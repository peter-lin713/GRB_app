#!/usr/bin/env Rscript
# LASSO ranking on the no-calibration (plain-MICE) dataset: same 100-rep
# cv.glmnet procedure used everywhere else in this pipeline, to see whether
# log10Fa/log10Ta/Alpha/Beta still rank as important predictors when optical
# GRBs get them via plain MICE instead of emcee/theilsen/ransac calibration.
suppressMessages(library(glmnet))

setwd("/Users/petelin/Desktop/test/GRB_app/plateau_ablation_experiments")
d <- read.csv("no_calibration_mice_only.csv", header = TRUE)

FEATURES <- c("log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
              "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux")
X <- as.matrix(d[, FEATURES])
Y <- log10(d$Redshift_crosscheck + 1)

LASSO <- function(X, Y) cv.glmnet(X, Y, alpha = 1)

lasso_coef <- vector()
for (a in 1:100) {
  lasmod <- LASSO(X, Y)
  lasso_coef <- cbind(lasso_coef, lasmod$glmnet.fit$beta[, lasmod$glmnet.fit$lambda == lasmod$lambda.1se])
}
lasso_coef_avg <- rowMeans(lasso_coef)

cat("\n--- Mean LASSO coefficients (no-calibration, plain-MICE dataset, N=", nrow(d), ") ---\n", sep = "")
print(lasso_coef_avg)

ranked <- names(sort(abs(lasso_coef_avg), decreasing = TRUE))
cat("\nFull ranking by |coefficient|:\n")
print(ranked)

cat("\nPlateau params (Fa/Ta/Alpha/Beta) rank position out of 10:\n")
for (v in c("log10Fa", "log10Ta", "Alpha", "Beta")) {
  cat(sprintf("  %-10s coef=%8.5f  rank=%d/10\n", v, lasso_coef_avg[v], which(ranked == v)))
}

saveRDS(list(coef = lasso_coef_avg, ranked = ranked), "lasso_no_calibration_result.rds")
