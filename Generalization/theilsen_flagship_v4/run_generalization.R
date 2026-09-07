#!/usr/bin/env Rscript
# Generalization analysis using the CONFIRMED-EXACT flagship Theil-Sen model.
suppressMessages({
  library(mice)
  library(SuperLearner)
  library(kSamples)
})
setwd("/Users/petelin/Desktop/test/GRB_app")
source("Generalization/Bias_Correction_function.R")
# Custom learners + their predict.* S3 methods -- required in this session
# before predict(sl_model, ...) can dispatch on SL.custom_glm etc.
for (f in list.files("Custom_SL", pattern = "\\.R$", full.names = TRUE)) source(f)

RUN_DIR <- "runs/theilsen_cleaned_linearz_5pct_mice20_RETRAIN_nesteddata"
sl_model <- readRDS(file.path(RUN_DIR, "superlearner_model"))
train_imputed <- read.csv(file.path(RUN_DIR, "OutputFiles", "grb_xray_imputed.csv"))
cv_results <- read.csv(file.path(RUN_DIR, "Results", "Results_wo_catout_with_catOutl_correlation_plot.csv"))

cat("varNames required by model:\n"); print(sl_model$varNames)

SqrTermGen <- function(inputData) {
  for (nm in colnames(inputData)) inputData[[paste0(nm, "Sqr")]] <- inputData[[nm]]^2
  inputData
}

# ---- 1. Load + remap TOTAL_GENERALIZATION_DATA_v4.csv -----------------------
gen_raw <- read.csv("TOTAL_GENERALIZATION_DATA_v4.csv", colClasses = "character")
cat("\nRaw generalization set:", nrow(gen_raw), "GRBs\n")

# Data-quality issue found this session: photon_index and logPeakFlux contain
# contaminated string values (trailing unit/model-tag suffixes like "1.81PL",
# stray trailing commas like "2.03,", and one literal unevaluated formula
# string "0.43429448190325176*Log[\"n/a\"]") that silently force the whole
# column to character on read.csv(). Strip to the leading numeric token;
# anything left non-numeric (the formula-string row, blank cells) becomes NA
# for MICE to impute, same as every other missing value in this pipeline.
clean_numeric <- function(x) as.numeric(sub("[^0-9.-].*$", "", x))
n_contaminated <- sum(is.na(suppressWarnings(as.numeric(gen_raw$photon_index)))) +
                    sum(is.na(suppressWarnings(as.numeric(gen_raw$logPeakFlux))))
for (col in c("Fbest","T_abest","F_min","F_max","T_amin","T_amax","logT90",
               "logPeakFlux","errorlogPeakFlux","photon_index","errorphotonindex",
               "logNH","T90_err","Gamma","Alpha","AlphaErr")) {
  gen_raw[[col]] <- clean_numeric(gen_raw[[col]])
}
cat("Cleaned", n_contaminated, "contaminated non-numeric string values (trailing unit/model-tag",
     "suffixes and stray commas in photon_index/logPeakFlux) -> numeric or NA for MICE\n")

gen_dat_preds <- data.frame(
  GRB           = gen_raw$X,
  log10Fa       = gen_raw$Fbest,
  log10Ta       = gen_raw$T_abest,
  log10T90      = gen_raw$logT90,
  log10PeakFlux = gen_raw$logPeakFlux,
  PhotonIndex   = gen_raw$photon_index,
  log10NH       = gen_raw$logNH,
  Gamma         = gen_raw$Gamma,
  Alpha         = gen_raw$Alpha,
  Beta          = NA_real_,
  log10Fluence  = NA_real_
)

# Null out physically-impossible values (same thresholds as training script)
gen_dat_preds$log10NH[gen_dat_preds$log10NH < 20]     <- NA
gen_dat_preds$PhotonIndex[gen_dat_preds$PhotonIndex < 0] <- NA
gen_dat_preds$Gamma[gen_dat_preds$Gamma > 3]          <- NA
gen_dat_preds$Alpha[gen_dat_preds$Alpha > 3]          <- NA

# ---- 2. Impute the partially-missing predictors via MICE ---------------------
# Beta and log10Fluence are 100% missing (0 of 299 rows) in this generalization
# set -- mice's midastouch/PMM has zero real donor values to draw from for a
# fully-empty column, so it cannot impute them (confirmed: left all-NA when
# included). This matches why Generalization_Aditya_v1.R's own MICE call for
# the earlier v2/v3 generalization sets never included Beta/Fluence as columns
# at all. MICE here only handles the genuinely-partially-missing columns
# (log10T90, log10PeakFlux, PhotonIndex, log10NH, Gamma, Alpha).
gen_dat_mice_cols <- setdiff(colnames(gen_dat_preds), c("GRB", "Beta", "log10Fluence"))
set.seed(1)
mice_gen <- mice(data = gen_dat_preds[, gen_dat_mice_cols], m = 20, maxit = 20,
                  method = "midastouch", printFlag = FALSE)
gen_imputed <- complete(mice_gen, 20)
gen_imputed$GRB <- gen_dat_preds$GRB

# Beta/log10Fluence: training-population mean substitution (a flagged
# limitation, not a silent one -- see report). Every generalization GRB gets
# the identical placeholder for these two features; this only damps the
# contribution of the two generic ensemble members that use them (SL.xgboost,
# SL.caret.rpart, ~38% combined weight) -- the formula-based GLM learner
# (~62% weight) never references Beta or Fluence at all.
gen_imputed$Beta         <- mean(train_imputed$Beta)
gen_imputed$log10Fluence <- mean(train_imputed$log10Fluence)
cat(sprintf("Beta/log10Fluence: 100%% missing in generalization set (0/%d observed) -- mice cannot\n", nrow(gen_imputed)),
     sprintf("  impute a fully-empty column; substituted training-population means (Beta=%.3f, log10Fluence=%.3f)\n",
              mean(train_imputed$Beta), mean(train_imputed$log10Fluence)))

# ---- 3. Remove out-of-parameter-space GRBs (avoid extrapolation) ------------
# Same feature set + convention as Generalization_Aditya_v1.R (lines ~481-500),
# extended with Alpha since v4 (unlike v2/v3) actually provides real Alpha.
range_cols <- c("log10NH", "PhotonIndex", "log10Fa", "log10Ta",
                 "log10PeakFlux", "log10T90", "Gamma", "Alpha")
in_range <- rep(TRUE, nrow(gen_imputed))
for (col in range_cols) {
  rng <- range(train_imputed[[col]])
  in_range <- in_range & gen_imputed[[col]] > rng[1] & gen_imputed[[col]] < rng[2]
}
n_before <- nrow(gen_imputed)
gen_filtered <- gen_imputed[in_range, ]
n_after <- nrow(gen_filtered)
cat(sprintf("Out-of-parameter-space filter: %d -> %d (%d removed, %.1f%%)\n",
             n_before, n_after, n_before - n_after, 100*(n_before-n_after)/n_before))

# ---- 4. Predict --------------------------------------------------------------
gen_sqr <- SqrTermGen(gen_filtered[, setdiff(colnames(gen_filtered), "GRB"), drop = FALSE])
newdata <- gen_sqr[, sl_model$varNames]
stopifnot(!anyNA(newdata))

pred <- predict(sl_model, newdata = newdata, onlySL = TRUE)$pred[, 1]
gen_filtered$InvZphot <- pred
gen_filtered$Zphot    <- 10^pred - 1
cat(sprintf("\nPoint predictions: N=%d, Zphot range [%.2f, %.2f], median %.2f\n",
             nrow(gen_filtered), min(gen_filtered$Zphot), max(gen_filtered$Zphot), median(gen_filtered$Zphot)))

# ---- 5. Bias correction (BC_3way, genuine mode: fit on Zspec, route by Zphot) -
png(filename = "/private/tmp/claude-501/-Users-petelin-Desktop-test-GRB-app/e05be0d8-9ae7-4510-8c50-2e6d5f9d1f7a/scratchpad/BC_3way_diagnostic_plots.png", width=1200, height=900)
par(mfrow=c(2,2))
gen_filtered$corrected_InvZphot <- BC_3way(crossvalidated_results = cv_results,
                                             generalization_set = gen_filtered,
                                             cut1 = 2, cut2 = 3.5)
dev.off()
gen_filtered$corrected_Zphot <- 10^gen_filtered$corrected_InvZphot - 1
cat(sprintf("Bias-corrected predictions: Zphot range [%.2f, %.2f], median %.2f\n",
             min(gen_filtered$corrected_Zphot), max(gen_filtered$corrected_Zphot), median(gen_filtered$corrected_Zphot)))

# ---- 6. Distribution validation ----------------------------------------------
cat("\n--- Distribution validation: generalization predicted-z vs training observed-z ---\n")
cat("Anderson-Darling (raw Zphot vs training Zspec):\n")
print(ad.test(gen_filtered$Zphot, cv_results$Zspec))
cat("\nAnderson-Darling (corrected Zphot vs training Zspec):\n")
print(ad.test(gen_filtered$corrected_Zphot, cv_results$Zspec))
cat("\nKS test (raw Zphot vs training Zspec): ")
print(ks.test(gen_filtered$Zphot, cv_results$Zspec))
cat("KS test (corrected Zphot vs training Zspec): ")
print(ks.test(gen_filtered$corrected_Zphot, cv_results$Zspec))

cat("\n--- Per-feature distribution checks (generalization vs training) ---\n")
for (col in c("log10T90", "log10NH", "log10Fa", "log10Ta", "PhotonIndex", "Alpha")) {
  kt <- ks.test(gen_filtered[[col]], train_imputed[[col]])
  cat(sprintf("%-14s KS p-value = %.4f %s\n", col, kt$p.value, if (kt$p.value < 0.05) "(DIFFERS)" else ""))
}

# ---- 7. Save outputs ----------------------------------------------------------
out_cols <- c("GRB", range_cols, "Beta", "log10Fluence", "InvZphot", "Zphot",
               "corrected_InvZphot", "corrected_Zphot")
write.csv(gen_filtered[, out_cols], "/private/tmp/claude-501/-Users-petelin-Desktop-test-GRB-app/e05be0d8-9ae7-4510-8c50-2e6d5f9d1f7a/scratchpad/generalization_predictions_v4.csv", row.names = FALSE)
cat("\nWrote generalization_predictions_v4.csv (", nrow(gen_filtered), "GRBs )\n")
