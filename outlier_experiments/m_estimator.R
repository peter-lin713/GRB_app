#' m_estimator.R — robust-regression outlier cut.
#'
#' Removes catastrophic outliers from the imputed GRB feature set before model
#' training. A Huber M-estimator (MASS::rlm, method = "M") fits log10(z+1) on
#' the available linear predictors; its iteratively-reweighted least-squares
#' procedure assigns each GRB a weight in (0, 1] — points the robust fit cannot
#' accommodate get down-weighted toward 0. GRBs whose weight falls at or below
#' `weight_threshold` are dropped as outliers.
#'
#' Inputs (globals expected from the calling pipeline):
#'   GRBPred           imputed feature dataframe (rows = GRBs); must contain
#'                     log10z plus the predictor columns. Mutated in place: the
#'                     outlier-cut frame is written back to GRBPred.
#'   PLOTaddr          output directory prefix for diagnostic PNGs
#'   weight_threshold  numeric in (0, 1); GRBs with weight <= it are dropped
#'
#' Outputs:
#'   <PLOTaddr>M_estimator_plots/M_estimator_weights.png       weight histogram
#'   <PLOTaddr>M_estimator_plots/ScatterPlot_w_M_est_Weights.png  pairs plot
#'   OutputFiles/grb_xray_m_est.csv                            outlier-cut data
#'   GRBPred (global)  reassigned to the outlier-cut frame

require(MASS)
require(dplyr)

input_dir  <- if (exists("out_files_dir")) out_files_dir else "OutputFiles"
output_dir <- input_dir

# Create I/O directories if they don't exist.
if (!dir.exists(input_dir)) {
  dir.create(input_dir)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Directory for this stage's diagnostic plots.
MEstAddr <- paste(PLOTaddr, "M_estimator_plots", sep = "")
if (!dir.exists(MEstAddr)) {
  dir.create(MEstAddr)
}

# Stash the redshift columns, then strip everything the regression must not see
# as a predictor (the response and its transforms / the crosscheck flag).
rc <- GRBPred$Redshift_crosscheck
Response <- GRBPred$log10z
GRBPred <- subset(GRBPred, select = !colnames(GRBPred) %in% c("Redshift_crosscheck", "log10z", "invz", "X"))

# --- Build the robust-regression formula -----------------------------------
# We deliberately do NOT read Best_formula_GLM.txt here: it hardcodes squared
# terms (log10NHSqr, log10T90Sqr, ...) that don't exist in this frame yet
# (SqrTermGen runs later, inside the CV loop, on the top-7 LASSO vars only), so
# those formulas would need refitting. Instead we build a generic additive
# robust regression over the available linear predictors so outlier weights can
# still be computed.
numeric_cols <- names(GRBPred)[vapply(GRBPred, is.numeric, logical(1))]
predictor_cols <- setdiff(numeric_cols,
                          c("log10z", "Redshift_crosscheck", "invz", "X",
                            grep("Err$", colnames(GRBPred), value = TRUE)))
# Drop linearly dependent predictors so rlm gets a full-rank design (rlm does
# not implement singular fits). QR pivoting keeps a maximal independent subset.
.mm <- model.matrix(reformulate(predictor_cols), data = GRBPred)
.qr <- qr(.mm)
.indep <- setdiff(colnames(.mm)[.qr$pivot[seq_len(.qr$rank)]], "(Intercept)")
rlm_form <- list(reformulate(.indep, response = "log10z"))

grb_cols <- colnames(GRBPred)

# Fit the M-estimator. log10z was stripped from GRBPred above and stashed in
# Response; add it back into the fitting frame so the formula's response
# resolves locally.
M_est <- MASS::rlm(rlm_form[[1]], data = cbind(GRBPred, log10z = Response), method = "M")
summary(M_est)

# Per-GRB robustness weights: 1 = well fit, low = outlier. The higher the
# threshold, the more points get classified as outliers and removed.
weights <- M_est$w

pct <- if (exists("m_est_pct")) m_est_pct else 0.05
weight_threshold <- quantile(M_est$w, pct)

# Diagnostic: histogram of weights with the cutoff marked.
{
  png(filename = paste(MEstAddr, "/M_estimator_weights.png", sep = ""), width = 1000, height = 1000, res = 200)
  hist_weights <- hist(weights, breaks = 2 + length(unique(round(weights, 1))),
                       cex.main = 0.95,
                       main = paste0("M-estimator weights per GRB\n",
                                     "(red = cutoff ", weight_threshold, "; below = dropped outliers)"),
                       xlab = "M-estimator weight (1 = well fit, low = outlier)")
  abline(v = weight_threshold, col = "red")
  legend("topleft",
    legend = rownames(GRBPred)[weights < weight_threshold]
  )
  dev.off()
}

# Colour code: 1 = kept (weight above threshold), 2 = outlier (at/below).
color_weights <- vector(length = nrow(GRBPred))
names(color_weights) <- rownames(GRBPred)
color_weights <- ifelse(weights > weight_threshold, 1, 2)

# Diagnostic: pairwise feature scatter coloured by kept/outlier status.
{
  png(file = paste(MEstAddr, "/ScatterPlot_w_M_est_Weights.png", sep = ""), width = 3000, height = 3000, res = 200)
  pairs(as.matrix(GRBPred[, vapply(GRBPred, is.numeric, logical(1))]), # numeric cols only (GRB is character)
    horOdd = T,
    pch = 3,
    col = color_weights,
    cex = 0.5,
    cex.labels = 1.4,
    main = paste0("Pairwise feature scatter (", dim(GRBPred)[1],
                  " GRBs) colored by M-estimator weight\n",
                  "black = kept (weight > ", weight_threshold,
                  "), red = catastrophic outlier (weight <= ", weight_threshold, ")")
  )
  legend(
    title = "Weights",
    "bottomleft",
    inset = 0.1,
    fill = unique(color_weights),
    col = unique(color_weights),
    legend = c(
      paste("weight more than ", weight_threshold, " (", hist(color_weights, breaks = 2, plot = F)$counts[1], ")", sep = ""),
      paste("weight less than ", weight_threshold, " (", hist(color_weights, breaks = 2, plot = F)$counts[2], ")", sep = "")
    )
  )
  dev.off()
}

# Restore the redshift columns, then keep only the non-outlier rows.
GRBPred$Redshift_crosscheck <- rc
GRBPred$log10z <- Response

GRBCut <- GRBPred[M_est$w > weight_threshold, ]

nrow(GRBCut)

# Report which GRBs were dropped as outliers.
print(rownames(GRBPred[M_est$w <= weight_threshold, ]))

GRBPred <- GRBCut
write.csv(GRBCut, paste(output_dir, "/grb_xray_m_est.csv", sep = ""))
