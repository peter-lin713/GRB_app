#'
#' Generates the M Estimator for removing outliers in data
#' Should be run from main code

require(MASS)
require(dplyr)


input_dir <- "OutputFiles"
output_dir <- "OutputFiles"

# CREATE DIRECTORIES IF THEY DONT EXIST
if (!dir.exists(input_dir)) {
  dir.create(input_dir)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

#MEstAddr <- paste(output_dir, "/Plots", sep = "")
MEstAddr <- paste(PLOTaddr, "M_estimator_plots", sep = "")
# creating dir for plots
if (!dir.exists(MEstAddr)) {
  dir.create(MEstAddr)
}

#GRBPred <- read.csv(paste(input_dir, "/grb_xray_imputed.csv", sep = ""), header = TRUE, row.names = 1)

rc <- GRBPred$Redshift_crosscheck


Response = GRBPred$log10z

### save all the information into a stable dataframe
GRBPred_full <- GRBPred
GRBPred <- subset(GRBPred, select = !colnames(GRBPred) %in% c("Redshift_crosscheck","log10z","invz","X", "GRB"))

# # 80/20 HIGHEST CV CORRELATION

# rlm_formula = read.table("Best_formula_GLM.txt")
# rlm_form = apply(as.matrix(rlm_formula[,2]),1,as.formula)

# #rlm_form <- log10z ~ (log10NHSqr + log10T90Sqr + log10TaSqr + log10FaSqr + log10NH + PhotonIndex + log10T90 + log10Ta)^2 + log10Fa + log10PeakFlux + PhotonIndexSqr + log10PeakFluxSqr

# # RLM FORMULA FOR TASK 1
# # rlm_form <- log10z ~ (log10NHSqr + log10T90Sqr + log10NH + log10T90)^2 + log10Fa + log10Ta + PhotonIndex + log10PeakFlux + log10FaSqr + log10TaSqr + PhotonIndexSqr + log10PeakFluxSqr


# # 80/20 AFTER NEW LOG10NH CUT
# #rlm_form <- log10z ~ (PhotonIndex + log10T90 + log10Ta + log10Fa)^2 + log10NH + log10PeakFlux

# #rlm_form <- log10z ~ (log10T90Sqr + PhotonIndex + log10Ta + log10Fa)^2 + log10T90 + log10NH + log10PeakFlux + log10FaSqr + log10TaSqr + PhotonIndexSqr + log10NHSqr + log10PeakFluxSqr



# # 90/10 HIGHEST CV CORRELATION
# #rlm_form <- log10z ~ (log10NHSqr + log10T90Sqr + log10NH + log10T90)^2 + log10Fa + log10Ta + PhotonIndex + log10PeakFlux + log10FaSqr + log10TaSqr + PhotonIndexSqr + log10PeakFluxSqr

# grb_cols <- colnames(GRBPred)



# # creating a linear combo of all predictors
# # all fitting to linear redshift
# # for (i in 1:length(grb_cols)-1) {
# #   if (i == 1) {
# #     assign("rlm_form", paste("log10z ~", grb_cols[i]))
# #   } else {
# #     assign("rlm_form", paste(rlm_form, "+", grb_cols[i]))
# #   }
# # }

# # creating M Estimate object and showing summary
# M_est <- MASS::rlm(rlm_form[[1]], data = GRBPred, method = "M")
# summary(M_est)

# # saving weight values and threshold
# weights <- M_est$w

# Save predictor-only rownames for stable mapping
stopifnot(length(Response) == nrow(GRBPred))
full_rows <- rownames(GRBPred)

# Reconstruct "full-ish" frame for fitting (predictors + response)
GRBPred_fit <- GRBPred
GRBPred_fit$log10z <- Response

# Build GLM formula: response ~ all predictors currently in GRBPred
pred_cols <- colnames(GRBPred)  # predictors-only at this point
glm_form <- as.formula(paste("log10z ~", paste(pred_cols, collapse = " + ")))

# Fit GLM (drops NA rows)
glm_fit <- glm(glm_form, data = GRBPred_fit, family = gaussian(), na.action = na.omit)

# Get the exact rows used by glm (complete cases for the model)
mf <- model.frame(glm_fit)
used_rows <- rownames(mf)

# Build design matrix aligned to mf
X <- model.matrix(glm_fit)

# Reduce X to full rank (rlm requires non-singular X)
qrX <- qr(X)
keep <- qrX$pivot[seq_len(qrX$rank)]
X_reduced <- X[, keep, drop = FALSE]

# Drop intercept so we can use log10z ~ . cleanly
Xr <- X_reduced
if ("(Intercept)" %in% colnames(Xr)) {
  Xr <- Xr[, colnames(Xr) != "(Intercept)", drop = FALSE]
}

# Fit robust regression on the same rows, reduced design
GRBPred_fit_reduced <- data.frame(
  log10z = model.response(mf),
  Xr,
  check.names = FALSE
)

M_est <- MASS::rlm(log10z ~ ., data = GRBPred_fit_reduced, method = "M")
weights <- M_est$w  # weights for *used_rows* only

# Map weights back to full row space (length = nrow(GRBPred))
w_full <- rep(NA_real_, length(full_rows))
names(w_full) <- full_rows
w_full[used_rows] <- weights

keep_full <- !is.na(w_full) & (w_full > weight_threshold)

# Reconstruct full GRBPred with columns restored, then cut rows
cat("Total rows:", length(full_rows), "\n")
cat("Rows used by glm:", length(used_rows), "\n")
cat("Rows kept after weight threshold:", sum(keep_full), "\n")

print("printing the type of weights")
str(weights)
class(weights)
head(weights)

print("seeing which columns are non-numeric")
non_numeric <- names(GRBPred)[!sapply(GRBPred, is.numeric)]
print(non_numeric)

sapply(GRBPred[non_numeric], class)


cat("Smallest 10 weights:\n")
print(sort(weights)[1:10])

cat("\nLargest 10 weights:\n")
print(sort(weights, decreasing = TRUE)[1:10])



# the higher the threshold, the more likely a data point 
# will be classified as an outlier and is removed
#weight_threshold <- 0.5

# showing weight cutoff via histogram
{
  png(filename = paste(MEstAddr, "/M_estimator_weights.png", sep = ""), width = 1000, height = 1000, res = 200)
  hist_weights <- hist(weights, breaks = 2 + length(unique(round(weights, 1))))
  abline(v = weight_threshold, col = "red")
  legend("topleft",
    legend = rownames(GRBPred)[weights < weight_threshold]
  )
  dev.off()
}

color_weights <- vector(length = nrow(GRBPred))
names(color_weights) <- rownames(GRBPred)

color_weights <- ifelse(weights > weight_threshold, 1, 2)

# creating pairplot of data once weights have been assigned colors
{
  png(file = paste(MEstAddr, "/ScatterPlot_w_M_est_Weights.png", sep = ""), width = 3000, height = 3000, res = 200)
  pairs(as.matrix(GRBPred), # the excluded columns are the error columns
    horOdd = T,
    pch = 3,
    # col=round(weights*11),
    col = color_weights,
    cex = 0.5,
    cex.labels = 1.4,
    # cex.angle=45,
    main = paste("Scatter Plot of", dim(GRBPred)[1], " samples")
    # lower.panel = as.matrix(lowredshift[,3:6])
  )
  legend(
    title = "Weights",
    "bottomleft",
    inset = 0.1,
    # legend = unique(fullDatMat_WO$Label),
    # fill = sort(unique(round(weights*11))),#unique(fullDatMat_WO$LabelNo),
    # col = sort(unique(round(weights*11))),
    # legend = paste(sort(unique(round(weights,1))))

    fill = unique(color_weights), # unique(fullDatMat_WO$LabelNo),
    col = unique(color_weights),
    legend = c(
      paste("weight more than ", weight_threshold, " (", hist(color_weights, breaks = 2, plot = F)$counts[1], ")", sep = ""),
      paste("weight less than ", weight_threshold, " (", hist(color_weights, breaks = 2, plot = F)$counts[2], ")", sep = "")
    )

    # paste(plyr::count(Totdata,vars = c('LabelNo','Label'))$Label,
    #               ' (',plyr::count(Totdata,vars = c('LabelNo','Label'))$freq,
    #              ')',sep = ''),
    # fill = plyr::count(Totdata,vars = c('LabelNo','Label'))$LabelNo,
    # col = plyr::count(Totdata,vars = c('LabelNo','Label'))$LabelNo,
    # pch = 16
  )

  dev.off()
  # par(mar=c(10,10,10,10))
}

# # removing data with weights below the threshold (these are considered outliers)
# GRBPred$Redshift_crosscheck <- rc
# GRBPred$log10z <- Response

# # removing data with weights lower than weight threshold
# GRBCut <- GRBPred[M_est$w > weight_threshold, ]

# nrow(GRBCut)

# print(rownames(GRBPred[M_est$w < weight_threshold, ]))

# GRBPred <- GRBCut

# GRBPred_full <- GRBPred
GRBPred_full$Redshift_crosscheck <- rc
GRBPred_full$log10z <- Response

GRBPred <- GRBPred_full[keep_full, , drop = FALSE]

write.csv(GRBPred, paste(output_dir, "/grb_xray_m_est.csv", sep = ""))
