#!/usr/bin/env Rscript
# Builds a version of the emcee-recovered dataset where optical-only GRBs'
# plateau parameters (log10Fa, log10Ta, Alpha, Beta) are NOT filled by the
# emcee calibration/projection -- they're set back to NA and left for plain
# MICE to impute, pooling ALL rows (so the imputation model is fit on the
# X-ray-native rows' real relationships and applied to fill the optical
# rows' gaps -- no dedicated optical->X-ray calibration involved at all).
#
# Used to test: how much do Fa/Ta/Alpha/Beta contribute (LASSO ranking, and
# a downstream SuperLearner ablation) when optical GRBs get these features
# via plain cross-domain MICE instead of a calibrated projection?
suppressMessages(require(mice))

setwd("/Users/petelin/Desktop/test/GRB_app")
raw_xray_data <- read.csv("Data/superlearner_training_emcee_v2_recovered.csv", header = TRUE, row.names = 1)

# Scope filter: long GRBs only (matches every other script in this pipeline).
if ("T90" %in% colnames(raw_xray_data)) {
  raw_xray_data <- raw_xray_data[raw_xray_data$T90 > 2, ]
  raw_xray_data$log10T90 <- log10(raw_xray_data$T90)
} else {
  raw_xray_data <- raw_xray_data[raw_xray_data$log10T90 > log10(2), ]
}

cat("N after scope filter:", nrow(raw_xray_data), "\n")
cat("is_optical distribution:", paste(names(table(raw_xray_data$is_optical)), table(raw_xray_data$is_optical), sep = "=", collapse = ", "), "\n")

# ---- Non-negotiable outlier cuts (new datacut, own method tag) -------------
source("generate_nonnegotiable_outliers.R")
.outliers <- update_nonnegotiable_outliers(raw_xray_data, out_dir = "plateau_ablation_experiments", method = "no_calibration_ablation")
to_drop <- .outliers$to_drop
raw_xray_data <- raw_xray_data[!(rownames(raw_xray_data) %in% to_drop), ]
cat("N after non-negotiable cuts:", nrow(raw_xray_data), "\n")

# ---- Strip the emcee calibration for optical-only GRBs ---------------------
# This is the key manipulation: undo the calibrated projection entirely for
# is_optical==1 rows, so plain MICE (not emcee/theilsen/ransac) has to fill
# these in using only cross-feature structure learned from the pooled sample.
opt_mask <- raw_xray_data$is_optical == 1
cat("Optical-only GRBs whose Fa/Ta/Alpha/Beta are being reset to NA:", sum(opt_mask), "\n")
raw_xray_data$log10Fa[opt_mask] <- NA
raw_xray_data$log10Ta[opt_mask] <- NA
raw_xray_data$Alpha[opt_mask]   <- NA
raw_xray_data$Beta[opt_mask]    <- NA
# Their errors on these params come from the (now-unused) calibration too --
# blank them so MICE imputes plausible errors alongside the values rather
# than keeping calibration-derived error estimates for values MICE invented.
raw_xray_data$log10FaErr[opt_mask] <- NA
raw_xray_data$log10TaErr[opt_mask] <- NA
raw_xray_data$AlphaErr[opt_mask]   <- NA
raw_xray_data$BetaErr[opt_mask]    <- NA

features_for_mice_preds <- subset(raw_xray_data, select = c(log10T90, log10Fa, log10Ta, Alpha, Beta,
                                                              Gamma, log10Fluence, PhotonIndex,
                                                              log10NH, log10PeakFlux))
# Same safety-net rules as every other script in this pipeline.
features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20] <- NA
features_for_mice_preds$Beta[features_for_mice_preds$Beta > 2] <- NA
features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3] <- NA
features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3] <- NA
features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA

cat("Missingness per feature before MICE:\n")
print(round(colMeans(is.na(features_for_mice_preds)), 3))

set.seed(1)
mice_model <- mice(data = features_for_mice_preds, m = 20, maxit = 20, method = "midastouch", printFlag = FALSE)
features_imputed <- complete(mice_model, 20)

# Errors: T90Err/FluenceErr/PhotonIndexErr/PeakFluxErr are untouched; the four
# plateau errors are NA for optical rows (blanked above alongside their
# values) and need their own MICE pass, same pattern as the main pipeline.
features_for_mice_errs <- subset(raw_xray_data, select = c(T90Err, log10FaErr, log10TaErr,
                                                             AlphaErr, BetaErr, log10FluenceErr,
                                                             PhotonIndexErr, log10PeakFluxErr))
mice_model_errs <- mice(data = features_for_mice_errs, m = 20, maxit = 20, method = "midastouch", printFlag = FALSE)
errs_imputed <- complete(mice_model_errs, 20)

out <- cbind(features_imputed, errs_imputed)
out$Redshift_crosscheck <- raw_xray_data$Redshift_crosscheck
out$is_optical <- raw_xray_data$is_optical
out$GRB <- rownames(raw_xray_data)

write.csv(out, "plateau_ablation_experiments/no_calibration_mice_only.csv", row.names = FALSE)
cat("\nWrote plateau_ablation_experiments/no_calibration_mice_only.csv --", nrow(out), "rows,", ncol(out), "cols\n")
