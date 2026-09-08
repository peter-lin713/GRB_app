# cv_permctrl.R — negative control (permuted target) + xray-only ceiling check.
# Run from GRB-Web-App repo root.
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
is_opt <- nmiss_lookup[dat$GRB] >= 6
FAST_LIBS <- c("SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost_safe", "SL.mean")

oof_r <- function(X, y, seed) {
  set.seed(seed)
  folds <- caret::createFolds(y, k = 5, returnTrain = FALSE)
  pred <- rep(NA_real_, length(y))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_along(y), te)
    fit <- suppressWarnings(SuperLearner(Y = y[tr], X = X[tr, , drop = FALSE],
      newX = X[te, , drop = FALSE], family = gaussian(), SL.library = FAST_LIBS,
      cvControl = list(V = 5), verbose = FALSE))
    pred[te] <- as.numeric(fit$SL.predict)
  }
  cor(y, pred)
}

X <- dat[, VARS]; Y <- dat$log10z
tasks <- lapply(1:6, function(p) list(p = p))
res <- parallel::mclapply(1:6, function(p) {
  set.seed(7000 + p)
  yp <- sample(Y)
  oof_r(X, yp, seed = 7100 + p)
}, mc.cores = 3L)
cat("permutation r values:", sprintf("%.3f", unlist(res)), "\n")
cat("perm mean:", round(mean(unlist(res)), 3), " sd:", round(sd(unlist(res)), 3), "\n")
