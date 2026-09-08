# cv_xray10f.R — X-ray-only subset, 10-fold CV (paper protocol), 3 repeats.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost)
})
source("Custom_SL/sl_xgboost_safe.R")

dat <- read.csv("outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv",
                row.names = 1, stringsAsFactors = FALSE)
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
VARS <- c(linear_vars, paste0(linear_vars, "Sqr"))
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)
keep <- nmiss_lookup[dat$GRB] < 6
dat <- dat[keep, ]
X <- dat[, VARS]; Y <- dat$log10z
cat("X-ray-only rows:", nrow(dat), "\n")
FAST_LIBS <- c("SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost_safe", "SL.mean")

res <- parallel::mclapply(1:3, function(rp) {
  set.seed(41000 + rp)
  folds <- caret::createFolds(Y, k = 10, returnTrain = FALSE)
  pred <- rep(NA_real_, length(Y))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_along(Y), te)
    fit <- suppressWarnings(SuperLearner(Y = Y[tr], X = X[tr, ], newX = X[te, ],
      family = gaussian(), SL.library = FAST_LIBS, cvControl = list(V = 5), verbose = FALSE))
    pred[te] <- as.numeric(fit$SL.predict)
  }
  cor(Y, pred)
}, mc.cores = 3L)
cat("xray-only 10-fold r:", sprintf("%.3f", unlist(res)), " mean:", round(mean(unlist(res)), 3), "\n")
