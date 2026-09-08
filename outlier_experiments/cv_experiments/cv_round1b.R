# cv_round1b.R — paper-inspired configs (Narendra/Dainotti 2025): formula GAM/GLM
# learners, top-7 LASSO features, physical combo features, is_optical indicator,
# high-z upsampling. Same paired 5x5 CV protocol/folds as cv_experiments.R.
# Run from GRB-Web-App repo root.

options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(MASS); library(mgcv); library(gam)
})
source("Custom_SL/sl_xgboost_safe.R")
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

out_csv <- commandArgs(trailingOnly = TRUE)[1]

dat <- read.csv("outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv",
                row.names = 1, stringsAsFactors = FALSE)
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
sqr_vars <- paste0(linear_vars, "Sqr")
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)

dat$is_optical <- as.numeric(nmiss_lookup[dat$GRB] >= 6)
dat$nmiss      <- as.numeric(nmiss_lookup[dat$GRB])
dat$fluxdecay  <- dat$log10Fa + dat$log10Ta          # plateau energy proxy
dat$meanflux   <- dat$log10Fluence - dat$log10T90    # mean prompt flux

Y <- dat$log10z; Z <- dat$Redshift_crosscheck; N <- nrow(dat)

N_REP <- 5L; K <- 5L
set.seed(31415)  # SAME folds as round 1
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# formula learners (paper: 6 GAM + 4 GLM + RF)
fg <- read.table("Best_formula_GAM.txt"); gam_forms <- apply(as.matrix(fg[, 2]), 1, as.formula)
fl <- read.table("Best_formula_GLM.txt"); glm_forms <- apply(as.matrix(fl[, 2]), 1, as.formula)
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1:min(6, length(gam_forms))]),
                       detailed_names = FALSE, name_prefix = "gamF")
glmL <- create.Learner("SL.custom_glm", tune = list(glm.model = glm_forms[1:min(4, length(glm_forms))]),
                       detailed_names = FALSE, name_prefix = "glmF")

FAST_LIBS  <- c("SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost_safe", "SL.mean")
PAPER_LIBS <- c(gamL$names, glmL$names, "SL.randomForest", "SL.mean")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")

CFG <- list(
  paperform   = list(vars = c(linear_vars, sqr_vars), libs = PAPER_LIBS),
  paperform_p = list(vars = c(linear_vars, sqr_vars), libs = c(PAPER_LIBS, FAST_LIBS)),
  top7        = list(vars = c(top7, paste0(top7, "Sqr")), libs = FAST_LIBS),
  physfeat    = list(vars = c(linear_vars, sqr_vars, "fluxdecay", "meanflux"), libs = FAST_LIBS),
  optind      = list(vars = c(linear_vars, sqr_vars, "is_optical", "nmiss"), libs = FAST_LIBS),
  upsample_hz = list(vars = c(linear_vars, sqr_vars), libs = FAST_LIBS, upsample = TRUE)
)

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  X <- dat[, cfg$vars, drop = FALSE]
  pred <- rep(NA_real_, N)
  set.seed(2000 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    if (isTRUE(cfg$upsample)) {          # duplicate z>2.5 training rows once
      extra <- tr[Z[tr] > 2.5]
      tr2 <- c(tr, extra)
    } else tr2 <- tr
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr2], X = X[tr2, , drop = FALSE], newX = X[te, , drop = FALSE],
      family = gaussian(), SL.library = cfg$libs,
      cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) return(data.frame(
      cfg = task$cfg, rep = task$rep, log_r = NA, log_rmse = NA, z_r = NA,
      error = substr(as.character(fit), 1, 200)))
    pred[te] <- as.numeric(fit$SL.predict)
  }
  z_pred <- 10^pred - 1
  data.frame(cfg = task$cfg, rep = task$rep,
             log_r = cor(Y, pred), log_rmse = sqrt(mean((Y - pred)^2)),
             z_r = cor(Z, z_pred), error = "")
}

tasks <- do.call(rbind, lapply(names(CFG), function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 3L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, out_csv, row.names = FALSE)
agg <- aggregate(cbind(log_r, log_rmse, z_r) ~ cfg, data = res, FUN = mean)
print(agg[order(-agg$log_r), ], digits = 3)
