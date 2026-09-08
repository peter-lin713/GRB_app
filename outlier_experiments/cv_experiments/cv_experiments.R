# cv_experiments.R — paired repeated-CV comparison of methods to improve CV r
# for GRB redshift prediction. Run from GRB-Web-App repo root:
#   Rscript cv_experiments.R <round>
# Data: post-MICE, post-5%-M-est frame from the combined_5pct run (n=291).
# Protocol: R repeats x K folds; the SAME fold assignments are used for every
# config, so per-repeat r values are paired across configs.

options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(kernlab); library(MASS)
})
source("Custom_SL/sl_xgboost_safe.R")

args   <- commandArgs(trailingOnly = TRUE)
round_id <- if (length(args) >= 1) args[1] else "1"
out_csv  <- if (length(args) >= 2) args[2] else paste0("cv_experiments_round", round_id, ".csv")

dat <- read.csv("outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv",
                row.names = 1, stringsAsFactors = FALSE)

linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
sqr_vars <- paste0(linear_vars, "Sqr")
err_vars <- c("log10T90Err","log10FaErr","log10TaErr","AlphaErr","BetaErr",
              "log10FluenceErr","PhotonIndexErr","log10PeakFluxErr")
Y  <- dat$log10z                      # log10(z+1), the standard target
Z  <- dat$Redshift_crosscheck
N  <- nrow(dat)

# ---- fixed folds shared by every config ----
N_REP <- 5L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# ---- learner variants ----
# tuned ranger: small-data settings
rng_tuned <- create.Learner("SL.ranger", tune = list(
  mtry = c(3L, 8L), min.node.size = c(5L, 15L)), detailed_names = TRUE, name_prefix = "rng")
# conservative xgboost variants for n~230 training rows
xgb_tuned <- create.Learner("SL.xgboost_safe", tune = list(
  max_depth = c(2L, 3L), shrinkage = c(0.03), ntrees = c(400L)),
  detailed_names = TRUE, name_prefix = "xgbT")
# glmnet ridge-ish variant
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

# ---- error-based observation weights: down-weight noisy measurements ----
err_mat <- as.matrix(dat[, err_vars])
err_mat[!is.finite(err_mat)] <- NA
tot_err <- rowMeans(scale(err_mat, center = FALSE,
                          scale = apply(err_mat, 2, function(c) stats::median(abs(c), na.rm = TRUE)))^2,
                    na.rm = TRUE)
w_err <- 1 / (1 + tot_err); w_err <- w_err / mean(w_err)

# ---- configs ----
# each: predictors (colnames), libs, obsWeights (NULL = uniform), method (SL metalearner),
#       y = target vector on model scale, to_log10z1 = fn mapping model-scale preds -> log10(z+1)
identity_map <- function(p) p
CFG <- list()
if (round_id == "1") {
  CFG$baseline      <- list(vars = c(linear_vars, sqr_vars), libs = FAST_LIBS)
  CFG$linear10      <- list(vars = linear_vars,              libs = FAST_LIBS)
  CFG$richlib       <- list(vars = c(linear_vars, sqr_vars), libs = RICH_LIBS)
  CFG$richlib_lin   <- list(vars = linear_vars,              libs = RICH_LIBS)
  CFG$errweight     <- list(vars = c(linear_vars, sqr_vars), libs = FAST_LIBS, w = w_err)
  CFG$errfeat       <- list(vars = c(linear_vars, sqr_vars, err_vars), libs = FAST_LIBS)
  CFG$target_logz   <- list(vars = c(linear_vars, sqr_vars), libs = FAST_LIBS,
                            y = log10(Z), map = function(p) log10(10^p + 1))
  CFG$meta_ccls     <- list(vars = c(linear_vars, sqr_vars), libs = FAST_LIBS,
                            method = "method.CC_LS")
} else {
  # round 2 configs are appended by editing this file
}

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  y_model <- if (!is.null(cfg$y)) cfg$y else Y
  map_fn  <- if (!is.null(cfg$map)) cfg$map else identity_map
  w       <- if (!is.null(cfg$w)) cfg$w else rep(1, N)
  X       <- dat[, cfg$vars, drop = FALSE]
  method  <- if (!is.null(cfg$method)) cfg$method else "method.NNLS"
  pred <- rep(NA_real_, N)
  set.seed(1000 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = y_model[tr], X = X[tr, , drop = FALSE], newX = X[te, , drop = FALSE],
      family = gaussian(), SL.library = cfg$libs, method = method,
      obsWeights = w[tr] / mean(w[tr]), cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) return(data.frame(
      cfg = task$cfg, rep = task$rep, log_r = NA, log_rmse = NA, z_r = NA,
      error = as.character(fit)))
    pred[te] <- map_fn(as.numeric(fit$SL.predict))
  }
  z_pred <- 10^pred - 1
  data.frame(cfg = task$cfg, rep = task$rep,
             log_r = cor(Y, pred), log_rmse = sqrt(mean((Y - pred)^2)),
             z_r = cor(Z, z_pred), error = "")
}

tasks <- do.call(rbind, lapply(names(CFG), function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks (", length(CFG), "configs x", N_REP, "repeats )\n")

res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, out_csv, row.names = FALSE)

cat("\n== per-config summary (mean over", N_REP, "paired repeats) ==\n")
agg <- aggregate(cbind(log_r, log_rmse, z_r) ~ cfg, data = res, FUN = mean)
agg <- agg[order(-agg$log_r), ]
print(agg, digits = 3)
cat("\n== paired delta vs baseline (log_r) ==\n")
base <- res[res$cfg == "baseline", c("rep", "log_r")]
for (cn in setdiff(unique(res$cfg), "baseline")) {
  m <- merge(base, res[res$cfg == cn, c("rep", "log_r")], by = "rep", suffixes = c("_b", "_x"))
  d <- m$log_r_x - m$log_r_b
  cat(sprintf("%-14s delta = %+0.4f (sd %0.4f)\n", cn, mean(d), sd(d)))
}
