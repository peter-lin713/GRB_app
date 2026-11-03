#!/usr/bin/env Rscript
# Data-prep stage of GAM_supercomputer/*/GAM_analysis_8variables.R, and nothing
# else: long-GRB filter -> non-negotiable cuts -> safety-net NA flags -> MICE
# (m=20, maxit=20, midastouch, seed 1). Output matches that script's
# grb_xray_imputed.csv (the pre-M-estimator frame).
#
# Deliberately skipped: the full-sample LASSO (selects features using every
# GRB's redshift), the formula search, and the M-estimator cut -- all of which
# leak redshift information into evaluation.
#
# Usage: Rscript prep_frame.R <input_csv> <out_dir> <method_label>

suppressPackageStartupMessages(library(mice))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]; out_dir <- args[2]; method <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
source(file.path(script_dir, "..", "generate_nonnegotiable_outliers.R"))

raw_xray_data <- read.csv(input_file, header = TRUE, row.names = 1)
raw_xray_data <- raw_xray_data[raw_xray_data$log10T90 > log10(2), ]

features_for_mice_preds <- subset(raw_xray_data, select = c(
  log10T90, log10Fa, log10Ta, Alpha, Beta, Gamma,
  log10Fluence, PhotonIndex, log10NH, log10PeakFlux))

.outliers <- update_nonnegotiable_outliers(raw_xray_data, out_dir = out_dir, method = method)
dropped_mask <- rownames(raw_xray_data) %in% .outliers$to_drop
raw_xray_data           <- raw_xray_data[!dropped_mask, ]
features_for_mice_preds <- features_for_mice_preds[!dropped_mask, ]

features_for_mice_preds$log10NH[features_for_mice_preds$log10NH < 20]        <- NA
features_for_mice_preds$Beta[features_for_mice_preds$Beta > 2]               <- NA
features_for_mice_preds$Gamma[features_for_mice_preds$Gamma > 3]             <- NA
features_for_mice_preds$Alpha[features_for_mice_preds$Alpha > 3]             <- NA
features_for_mice_preds$PhotonIndex[features_for_mice_preds$PhotonIndex < 0] <- NA
features_for_mice_preds$log10PeakFlux[is.infinite(features_for_mice_preds$log10PeakFlux)] <- NA

fluence_err_col  <- if ("FluenceErr"  %in% colnames(raw_xray_data)) "FluenceErr"  else "log10FluenceErr"
peakflux_err_col <- if ("PeakFluxErr" %in% colnames(raw_xray_data)) "PeakFluxErr" else "log10PeakFluxErr"
features_for_mice_errs <- raw_xray_data[, c("T90Err", "log10FaErr", "log10TaErr",
                                            "AlphaErr", "BetaErr", fluence_err_col,
                                            "PhotonIndexErr", peakflux_err_col)]

set.seed(1)
mice_model_preds <- mice(data = features_for_mice_preds, m = 20, maxit = 20,
                         method = 'midastouch', printFlag = F)
features_for_mice_preds <- complete(mice_model_preds, 20)
mice_model_errs <- mice(data = features_for_mice_errs, m = 20, maxit = 20,
                        method = 'midastouch', printFlag = F)
features_for_mice_errs <- complete(mice_model_errs, 20)

GRBPred <- cbind(features_for_mice_preds, features_for_mice_errs)
if (!"log10T90Err" %in% colnames(GRBPred)) {
  GRBPred$log10T90Err <- GRBPred$T90Err / (10^GRBPred$log10T90 * log(10))
}
GRBPred$Redshift_crosscheck <- raw_xray_data$Redshift_crosscheck
GRBPred$log10z <- log10(raw_xray_data$Redshift_crosscheck + 1)
GRBPred$is_optical <- raw_xray_data$is_optical

write.csv(GRBPred, file.path(out_dir, "grb_xray_imputed.csv"))
cat("Wrote", nrow(GRBPred), "GRBs ->", file.path(out_dir, "grb_xray_imputed.csv"), "\n")
