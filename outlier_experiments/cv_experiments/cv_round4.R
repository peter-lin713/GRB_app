# cv_round4.R — recover REAL prompt/X-ray features for the optical-projected rows
# from Data/optical_data.csv (units validated on 88 overlap GRBs), re-impute the
# small remainder with MICE, and re-run the winning formula-library configs on the
# SAME paired folds. Also: domain-weight control on the old frame.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(MASS); library(mgcv); library(gam); library(mice)
})
source("Custom_SL/sl_xgboost_safe.R")
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

out_prefix <- commandArgs(trailingOnly = TRUE)[1]

old <- read.csv("outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv",
                row.names = 1, stringsAsFactors = FALSE)
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
opt <- read.csv("Data/optical_data.csv", stringsAsFactors = FALSE); names(opt)[1] <- "GRB"

linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
sqr_vars <- paste0(linear_vars, "Sqr")
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)
is_opt <- nmiss_lookup[old$GRB] >= 6
cat("frame:", nrow(old), "rows;", sum(is_opt), "optical-projected\n")

# ---- rebuild feature block with real optical-catalog values ----
new <- old
mo <- match(new$GRB[is_opt], opt$GRB)
stopifnot(!anyNA(mo))
oc <- opt[mo, ]
put <- function(col, vals) { v <- new[[col]]; v[is_opt] <- vals; new[[col]] <<- v }
put("log10T90",      log10(oc$T90))
put("log10Fluence",  log10(oc$Fluence) - 7)
put("PhotonIndex",   oc$PhotonIndex)
put("log10NH",       log10(oc$NH) + 21)
put("log10PeakFlux", log10(ifelse(oc$PeakFlux > 0, oc$PeakFlux, NA)))
put("Gamma",         rep(NA_real_, sum(is_opt)))   # not measured for optical rows

# pipeline cleaning rules
new$log10NH[new$log10NH < 20] <- NA
new$PhotonIndex[new$PhotonIndex < 0] <- NA
new$log10PeakFlux[is.infinite(new$log10PeakFlux)] <- NA
cat("NAs to re-impute:", sum(is.na(new[, linear_vars])), "cells in",
    sum(!complete.cases(new[, linear_vars])), "rows\n")

set.seed(1)
imp <- mice(new[, linear_vars], m = 20, method = "midastouch", printFlag = FALSE)
new[, linear_vars] <- complete(imp, 20)
for (v in linear_vars) new[[paste0(v, "Sqr")]] <- new[[v]]^2
stopifnot(!anyNA(new[, c(linear_vars, sqr_vars)]))
write.csv(new, paste0(out_prefix, "_frame.csv"))

old$is_optical <- as.numeric(is_opt); new$is_optical <- as.numeric(is_opt)
Y <- old$log10z; Z <- old$Redshift_crosscheck; N <- nrow(old)
stopifnot(identical(old$GRB, new$GRB))

N_REP <- 5L; K <- 5L
set.seed(31415)  # SAME folds as all rounds
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

fg <- read.table("Best_formula_GAM.txt"); gam_forms <- apply(as.matrix(fg[, 2]), 1, as.formula)
fl <- read.table("Best_formula_GLM.txt"); glm_forms <- apply(as.matrix(fl[, 2]), 1, as.formula)
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1:6]), detailed_names = FALSE, name_prefix = "gamF")
glmL <- create.Learner("SL.custom_glm", tune = list(glm.model = glm_forms[1:4]), detailed_names = FALSE, name_prefix = "glmF")
PAPER_LIBS <- c(gamL$names, glmL$names, "SL.randomForest", "SL.mean")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
V7 <- c(linear_vars, paste0(top7, "Sqr"))

CFG <- list(
  new_top7     = list(frame = "new", vars = V7, w_opt = 1),
  new_top7_ind = list(frame = "new", vars = c(V7, "is_optical"), w_opt = 1),
  new_full     = list(frame = "new", vars = c(linear_vars, sqr_vars), w_opt = 1),
  old_w05      = list(frame = "old", vars = V7, w_opt = 0.5)
)

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  fr <- if (cfg$frame == "new") new else old
  X <- fr[, cfg$vars, drop = FALSE]
  w_all <- ifelse(is_opt, cfg$w_opt, 1)
  pred <- rep(NA_real_, N)
  set.seed(9000 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, , drop = FALSE], newX = X[te, , drop = FALSE],
      family = gaussian(), SL.library = PAPER_LIBS,
      obsWeights = w_all[tr] / mean(w_all[tr]),
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

cat("\n== round 4 (recovered real features): mean over 5 paired repeats ==\n")
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    c(all = cor(d$y, d$pred), xr = cor(d$y[!d$is_opt], d$pred[!d$is_opt]),
      op = cor(d$y[d$is_opt], d$pred[d$is_opt]))
  })
  cat(sprintf("%-14s r_all=%.3f  r_xray=%.3f  r_opt=%.3f\n",
              cn, mean(rr["all", ]), mean(rr["xr", ]), mean(rr["op", ])))
}
