#!/usr/bin/env Rscript
# Apply the ALREADY-KNOWN emcee winning formula (idx 861, r=0.70 run) to a new
# combination dataset (RANSAC / Theil-Sen), skipping the expensive 100-split
# formula search -- justified because lasso_check.R confirmed both datasets'
# top-7 LASSO features are identical (same set, same rank order) to emcee's.
#
# Replicates GAM_analysis_8variables.R's M-estimator block (lines ~814-839)
# exactly, but with the winning formula supplied directly instead of derived
# from formula_win_frequency.csv.
#
# Usage (run from inside e.g. GAM_supercomputer/ransac_cleaned_formula_generation/):
#   Rscript ../apply_known_formula_mest.R superlearner_training_ransac.csv ransac

suppressMessages({
  require(mice)
  require(MASS)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
method_tag <- args[2]

raw_xray_data <- read.csv(input_file, header = TRUE, row.names = 1)

if ("T90" %in% colnames(raw_xray_data)) {
  raw_xray_data <- raw_xray_data[raw_xray_data$T90 > 2, ]
  raw_xray_data$log10T90 <- log10(raw_xray_data$T90)
} else {
  raw_xray_data <- raw_xray_data[raw_xray_data$log10T90 > log10(2), ]
}

features_for_mice_preds <- subset(raw_xray_data, select = c(log10T90,
                                                              log10Fa,
                                                              log10Ta,
                                                              Alpha,
                                                              Beta,
                                                              Gamma,
                                                              log10Fluence,
                                                              PhotonIndex,
                                                              log10NH,
                                                              log10PeakFlux))

source("../../generate_nonnegotiable_outliers.R")
.outliers <- update_nonnegotiable_outliers(raw_xray_data, out_dir = ".", method = method_tag)
to_drop <- .outliers$to_drop

grb_ids <- rownames(raw_xray_data)
dropped_mask <- grb_ids %in% to_drop
if (any(dropped_mask)) {
  cat("Removed", sum(dropped_mask), "GRB(s) via non-negotiable cuts\n")
  raw_xray_data <- raw_xray_data[!dropped_mask, ]
  features_for_mice_preds <- features_for_mice_preds[!dropped_mask, ]
}

features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20] <- NA
features_for_mice_preds$Beta[features_for_mice_preds$Beta > 2] <- NA
features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3] <- NA
features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3] <- NA
features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA

cat("N after non-negotiable cuts:", nrow(features_for_mice_preds), "\n")

set.seed(1)
mice_model_preds <- mice(data = features_for_mice_preds, m = 20, maxit = 20, method = 'midastouch', printFlag = FALSE)
features_for_mice_preds <- complete(mice_model_preds, 20)

# Confirmed via lasso_check.R: both RANSAC and Theil-Sen top-7 == emcee's top-7,
# same rank order -- so the emcee winning formula's referenced columns all exist.
lassovar <- c("log10NH", "log10PeakFlux", "PhotonIndex", "log10Ta",
              "log10Fa", "log10T90", "Alpha")

SqrTermGen <- function(inputData) {
  indVar <- colnames(inputData)
  for (i in 1:length(indVar)) {
    inputData[[paste(indVar[i], "Sqr", sep = "")]] <- inputData[, indVar[i]] * inputData[, indVar[i]]
  }
  return(inputData)
}

Predictor <- subset(features_for_mice_preds, select = lassovar)
SqrData   <- SqrTermGen(Predictor)
SqrData$Response <- log10(raw_xray_data$Redshift_crosscheck + 1)
rownames(SqrData) <- rownames(raw_xray_data)

Formula_for_outlier <- as.formula(
  "Response ~ (log10FaSqr + log10T90 + log10Fa + log10PeakFlux)^2 + log10NH + PhotonIndex + log10Ta + Alpha + log10NHSqr + log10PeakFluxSqr + PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr"
)
cat("\nUsing known emcee winning formula (idx 861, r=0.70 run):\n")
print(Formula_for_outlier)

M_est <- MASS::rlm(Formula_for_outlier, data = SqrData, method = "M", maxit = 50)
weights <- M_est$w
weight_threshold <- quantile(weights, 0.05)

kept_rows    <- rownames(SqrData)[weights > weight_threshold]
removed_rows <- rownames(SqrData)[weights <= weight_threshold]

cat("\nM-estimator outlier cut:", length(removed_rows), "of", nrow(SqrData),
    "GRBs removed (weight <=", round(weight_threshold, 4), ")\n")
writeLines(removed_rows, "removed_outliers.txt")

final_cut <- raw_xray_data[rownames(raw_xray_data) %in% kept_rows, ]
write.csv(final_cut, "final_outliers_removed.csv")
cat("Final outlier-removed data saved:", nrow(final_cut), "GRBs -> final_outliers_removed.csv\n")
