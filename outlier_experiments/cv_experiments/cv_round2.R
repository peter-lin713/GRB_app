# cv_round2.R — round 2: (a) combine round-1 winners (rich lib + CC_LS),
# (b) X-ray-only vs all-rows training scored on the same X-ray rows,
# (c) CV-safe piecewise-linear bias correction on the best config.
# Same fold protocol/seed as rounds 1/1b. Run from GRB-Web-App repo root.

options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(kernlab); library(MASS)
})
source("Custom_SL/sl_xgboost_safe.R")

out_prefix <- commandArgs(trailingOnly = TRUE)[1]

dat <- read.csv("outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv",
                row.names = 1, stringsAsFactors = FALSE)
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
sqr_vars <- paste0(linear_vars, "Sqr")
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)
is_opt <- nmiss_lookup[dat$GRB] >= 6

Y <- dat$log10z; Z <- dat$Redshift_crosscheck; N <- nrow(dat)
cat("rows:", N, " optical-projected:", sum(is_opt), "\n")

N_REP <- 5L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

rng_tuned <- create.Learner("SL.ranger", tune = list(
  mtry = c(3L, 8L), min.node.size = c(5L, 15L)), detailed_names = TRUE, name_prefix = "rng")
xgb_tuned <- create.Learner("SL.xgboost_safe", tune = list(
  max_depth = c(2L, 3L), shrinkage = c(0.03), ntrees = c(400L)),
  detailed_names = TRUE, name_prefix = "xgbT")
glmnet_ridge <- create.Learner("SL.glmnet", tune = list(alpha = c(0, 0.5)),
  detailed_names = TRUE, name_prefix = "gnet")
SL.ksvm.rbf <- function(Y, X, newX, family, obsWeights, id, ...) {
  fit <- kernlab::ksvm(as.matrix(X), Y, kernel = "rbfdot", C = 1, epsilon = 0.1)
  pred <- kernlab::predict(fit, as.matrix(newX))[, 1]
  out <- list(object = fit); class(out) <- "SL.ksvm.rbf"
  list(pred = pred, fit = out)
}
predict.SL.ksvm.rbf <- function(object, newdata, ...) kernlab::predict(object$object, as.matrix(newdata))[, 1]

FAST_LIBS <- c("SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost_safe", "SL.mean")
RICH_LIBS <- c("SL.glmnet", glmnet_ridge$names, "SL.earth", "SL.mean", "SL.ksvm.rbf",
               rng_tuned$names, xgb_tuned$names)
VARS <- c(linear_vars, sqr_vars)

CFG <- list(
  combo      = list(libs = RICH_LIBS, method = "method.CC_LS", train_sub = "all"),
  trainall   = list(libs = FAST_LIBS, method = "method.NNLS",  train_sub = "all"),
  trainxray  = list(libs = FAST_LIBS, method = "method.NNLS",  train_sub = "xray")
)

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  X <- dat[, VARS, drop = FALSE]
  pred <- rep(NA_real_, N)
  set.seed(3000 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    if (cfg$train_sub == "xray") tr <- tr[!is_opt[tr]]
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, , drop = FALSE], newX = X[te, , drop = FALSE],
      family = gaussian(), SL.library = cfg$libs, method = cfg$method,
      cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) {
      cat("ERROR", task$cfg, task$rep, i, substr(as.character(fit), 1, 120), "\n")
      next
    }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(cfg = task$cfg, rep = task$rep, row = seq_len(N),
             grb = dat$GRB, is_opt = is_opt, y = Y, z = Z, pred = pred,
             stringsAsFactors = FALSE)
}

tasks <- do.call(rbind, lapply(names(CFG), function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

# ---- metrics ----
metr <- function(d) c(r = cor(d$y, d$pred, use = "complete.obs"),
                      rmse = sqrt(mean((d$y - d$pred)^2, na.rm = TRUE)))
cat("\n== per-config, per-subset r(log10(z+1)), mean over repeats ==\n")
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp, ]
    c(all  = unname(metr(d)["r"]),
      xray = unname(metr(d[!d$is_opt, ])["r"]),
      opt  = unname(metr(d[d$is_opt, ])["r"]))
  })
  cat(sprintf("%-10s r_all=%.3f  r_xray=%.3f  r_opt=%.3f\n",
              cn, mean(rr["all", ]), mean(rr["xray", ]), mean(rr["opt", ], na.rm = TRUE)))
}

# ---- CV-safe piecewise-linear bias correction (2-fold cross-fitted) ----
# Paper (sec 4.6) fits Zc = B + A*Zp in 3 prediction ranges; we cross-fit so the
# correction is never fit on the rows it is applied to.
bc_crossfit <- function(d, n_bins = 3) {
  set.seed(99); grp <- sample(rep(1:2, length.out = nrow(d)))
  out <- rep(NA_real_, nrow(d))
  for (g in 1:2) {
    fit_d <- d[grp != g, ]; app <- d[grp == g, ]
    qs <- quantile(fit_d$pred, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
    qs[1] <- -Inf; qs[n_bins + 1] <- Inf
    for (b in seq_len(n_bins)) {
      inb_fit <- fit_d$pred >  qs[b] & fit_d$pred <= qs[b + 1]
      inb_app <- app$pred  >  qs[b] & app$pred  <= qs[b + 1]
      if (sum(inb_fit, na.rm = TRUE) >= 10) {
        cf <- coef(lm(y ~ pred, data = fit_d[inb_fit, ]))
        out[which(grp == g)[which(inb_app)]] <- cf[1] + cf[2] * app$pred[inb_app]
      } else {
        out[which(grp == g)[which(inb_app)]] <- app$pred[inb_app]
      }
    }
  }
  out
}
cat("\n== combo config with cross-fitted piecewise bias correction ==\n")
for (rp in seq_len(N_REP)) {
  d <- res[res$cfg == "combo" & res$rep == rp & is.finite(res$pred), ]
  d$bc <- bc_crossfit(d)
  cat(sprintf("rep %d: raw r=%.3f  BC r=%.3f | linear-z: raw r=%.3f BC r=%.3f\n",
              rp, cor(d$y, d$pred), cor(d$y, d$bc, use = "complete.obs"),
              cor(d$z, 10^d$pred - 1), cor(d$z, 10^d$bc - 1, use = "complete.obs")))
}
