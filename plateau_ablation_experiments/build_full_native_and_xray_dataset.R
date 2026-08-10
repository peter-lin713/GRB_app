#!/usr/bin/env Rscript
# Builds a dataset where EVERY GRB has BOTH representations of the plateau
# params: the standard X-ray-scale (log10Fa/log10Ta/Alpha/Beta -- real for
# X-ray-native GRBs, emcee-projected for optical-only GRBs, already present
# in the recovered/native_opt data) AND the optical-native-scale versions
# (*_native -- real for optical-only GRBs, already present; and for X-ray-
# native GRBs, filled here via the INVERSE of the same emcee calibration
# line: optical = (xray - b) / m, using the posterior-median (m, b) per
# parameter from Data/emcee_chains_v2.npz -- the same calibration already
# used everywhere else in this pipeline, just applied in the other direction).
#
# Lets LASSO compare both representations directly for every GRB, rather
# than only for the ~73 optical-only GRBs that had a native measurement.
suppressMessages(require(mice))

setwd("/Users/petelin/Desktop/test/GRB_app")
d <- read.csv("Data/superlearner_training_emcee_v2_native_opt.csv", header = TRUE, row.names = 1)

# Posterior-median (m, b) per parameter, computed from Data/emcee_chains_v2.npz
# (same calibration this whole pipeline already uses for optical->X-ray).
calib <- list(
  log10Fa = c(m = 0.428340, b = -5.361624),
  log10Ta = c(m = 0.359438, b =  2.365699),
  Alpha   = c(m = 0.481452, b =  0.735096),
  Beta    = c(m = 0.093260, b =  0.831452)
)

xray_mask <- d$is_optical == 0
cat("X-ray-native GRBs getting an inverse-projected native-optical value:", sum(xray_mask), "\n")

for (core in names(calib)) {
  native_col     <- paste0(core, "_native")
  native_err_col <- paste0(core, "Err_native")  # log10FaErr_native / AlphaErr_native, matches naming already in the file

  m <- calib[[core]]["m"]; b <- calib[[core]]["b"]
  d[[native_col]][xray_mask] <- (d[[core]][xray_mask] - b) / m

  err_col <- paste0(core, "Err")
  if (err_col %in% colnames(d) && native_err_col %in% colnames(d)) {
    d[[native_err_col]][xray_mask] <- d[[err_col]][xray_mask] / abs(m)
  }
}

cat("Native-column missingness after inverse projection:\n")
print(colMeans(is.na(d[, c("log10Fa_native", "log10Ta_native", "Alpha_native", "Beta_native")])))

# ---- Scope filter + non-negotiable cuts ("new datacut") --------------------
if ("T90" %in% colnames(d)) {
  d <- d[d$T90 > 2, ]
  d$log10T90 <- log10(d$T90)
} else {
  d <- d[d$log10T90 > log10(2), ]
}
source("generate_nonnegotiable_outliers.R")
.outliers <- update_nonnegotiable_outliers(d, out_dir = "plateau_ablation_experiments", method = "full_native_and_xray")
d <- d[!(rownames(d) %in% .outliers$to_drop), ]
cat("N after cuts:", nrow(d), "\n")

# ---- MICE the full 14-feature set (10 core + 4 native) ---------------------
feat_cols <- c("log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
               "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux",
               "log10Fa_native", "log10Ta_native", "Alpha_native", "Beta_native")
feats <- d[, feat_cols]
feats$log10NH[feats$log10NH < 20] <- NA
feats$Beta[feats$Beta > 2] <- NA
feats$Gamma[feats$Gamma > 3] <- NA
feats$Alpha[feats$Alpha > 3] <- NA
feats$PhotonIndex[feats$PhotonIndex < 0] <- NA
feats$log10PeakFlux[is.infinite(feats$log10PeakFlux)] <- NA

cat("Missingness per feature before MICE:\n")
print(round(colMeans(is.na(feats)), 3))

set.seed(1)
mice_model <- mice(data = feats, m = 20, maxit = 20, method = "midastouch", printFlag = FALSE)
feats_imputed <- complete(mice_model, 20)
feats_imputed$Redshift_crosscheck <- d$Redshift_crosscheck
feats_imputed$is_optical <- d$is_optical
feats_imputed$GRB <- rownames(d)

write.csv(feats_imputed, "plateau_ablation_experiments/full_native_and_xray.csv", row.names = FALSE)
cat("\nWrote plateau_ablation_experiments/full_native_and_xray.csv --", nrow(feats_imputed), "rows,", ncol(feats_imputed), "cols\n")
