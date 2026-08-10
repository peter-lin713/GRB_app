#!/usr/bin/env Rscript
# Standalone LASSO feature-ranking check.
#
# Replicates the exact LASSO block from GAM_analysis_8variables.R (scope
# filter -> non-negotiable outlier cuts -> MICE(midastouch, m=20) on the 10
# core predictors -> cv.glmnet(alpha=1) x100 averaged at lambda.1se -> top-7
# by |mean coef|) without running the full 100-split formula search, so we
# can compare RANSAC's/Theil-Sen's top-7 against emcee's before deciding
# whether to reuse the existing winning formula.
#
# Usage (run from inside e.g. GAM_supercomputer/ransac_cleaned_formula_generation/):
#   Rscript ../lasso_check.R superlearner_training_ransac.csv ransac_lasso_check

suppressMessages({
  require(mice)
  library(glmnet)
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

cat("N after cuts:", nrow(features_for_mice_preds), "\n")

set.seed(1)
mice_model_preds <- mice(data = features_for_mice_preds, m = 20, maxit = 20, method = 'midastouch', printFlag = FALSE)
features_for_mice_preds <- complete(mice_model_preds, 20)

LASSO <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.vector(Y)
  cv.glmnet(X, Y, alpha = 1)
}

lasso_coef <- vector()
for (a in 1:100) {
  lasmod <- LASSO(Y = log10(raw_xray_data$Redshift_crosscheck + 1), X = features_for_mice_preds)
  lasso_coef <- cbind(lasso_coef, lasmod$glmnet.fit$beta[, lasmod$glmnet.fit$lambda == lasmod$lambda.1se])
}
lasso_coef_avg <- rowMeans(lasso_coef)

cat("\n--- Mean LASSO coefficients (", method_tag, ") ---\n", sep = "")
print(lasso_coef_avg)

lassovar <- names(lasso_coef_avg[order(abs(lasso_coef_avg), decreasing = TRUE)])
top7 <- head(lassovar, 7)

cat("\nTop-7 features (", method_tag, "):\n", sep = "")
print(top7)

saveRDS(list(coef = lasso_coef_avg, top7 = top7), paste0("lasso_check_", method_tag, ".rds"))
