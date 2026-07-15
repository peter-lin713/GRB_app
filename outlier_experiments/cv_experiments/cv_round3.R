# cv_round3.R — final round: combine the winners around the formula-learner
# library (paperform). Same paired folds as all prior rounds.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(MASS); library(mgcv); library(gam)
})
source("Custom_SL/sl_xgboost_safe.R")
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

out_prefix <- commandArgs(trailingOnly = TRUE)[1]

dat <- read.csv("outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv",
                row.names = 1, stringsAsFactors = FALSE)
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
sqr_vars <- paste0(linear_vars, "Sqr")
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)
is_opt <- nmiss_lookup[dat$GRB] >= 6
dat$is_optical <- as.numeric(is_opt); dat$nmiss <- as.numeric(nmiss_lookup[dat$GRB])

Y <- dat$log10z; Z <- dat$Redshift_crosscheck; N <- nrow(dat)
N_REP <- 5L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

fg <- read.table("Best_formula_GAM.txt"); gam_forms <- apply(as.matrix(fg[, 2]), 1, as.formula)
fl <- read.table("Best_formula_GLM.txt"); glm_forms <- apply(as.matrix(fl[, 2]), 1, as.formula)
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1:min(6, length(gam_forms))]),
                       detailed_names = FALSE, name_prefix = "gamF")
glmL <- create.Learner("SL.custom_glm", tune = list(glm.model = glm_forms[1:min(4, length(glm_forms))]),
                       detailed_names = FALSE, name_prefix = "glmF")
PAPER_LIBS <- c(gamL$names, glmL$names, "SL.randomForest", "SL.mean")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")

CFG <- list(
  pf          = list(vars = c(linear_vars, sqr_vars), libs = PAPER_LIBS, method = "method.NNLS", train_sub = "all"),
  pf_top7     = list(vars = c(linear_vars, paste0(top7, "Sqr")), libs = PAPER_LIBS, method = "method.NNLS", train_sub = "all"),
  pf_ccls     = list(vars = c(linear_vars, sqr_vars), libs = PAPER_LIBS, method = "method.CC_LS", train_sub = "all"),
  pf_opt      = list(vars = c(linear_vars, sqr_vars, "is_optical", "nmiss"), libs = PAPER_LIBS, method = "method.NNLS", train_sub = "all"),
  pf_xrayonly = list(vars = c(linear_vars, sqr_vars), libs = PAPER_LIBS, method = "method.NNLS", train_sub = "xray")
)

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  X <- dat[, cfg$vars, drop = FALSE]
  pred <- rep(NA_real_, N)
  set.seed(5000 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    if (cfg$train_sub == "xray") tr <- tr[!is_opt[tr]]
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, , drop = FALSE], newX = X[te, , drop = FALSE],
      family = gaussian(), SL.library = cfg$libs, method = cfg$method,
      cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", task$cfg, task$rep, i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(cfg = task$cfg, rep = task$rep, row = seq_len(N), is_opt = is_opt,
             y = Y, z = Z, pred = pred, stringsAsFactors = FALSE)
}

tasks <- do.call(rbind, lapply(names(CFG), function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 3: r(log10(z+1)) mean over 5 paired repeats ==\n")
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    c(all = cor(d$y, d$pred), xray = cor(d$y[!d$is_opt], d$pred[!d$is_opt]))
  })
  cat(sprintf("%-12s r_all=%.3f  r_xray=%.3f\n", cn, mean(rr["all", ]), mean(rr["xray", ])))
}
