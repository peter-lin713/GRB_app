# cv_round5.R — final: average predictions over 3 MICE completions on the
# recovered-features frame (new_top7_ind config). Anchor = single completion.
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
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)
is_opt <- nmiss_lookup[old$GRB] >= 6

new <- old
mo <- match(new$GRB[is_opt], opt$GRB); oc <- opt[mo, ]
put <- function(col, vals) { v <- new[[col]]; v[is_opt] <- vals; new[[col]] <<- v }
put("log10T90",      log10(oc$T90))
put("log10Fluence",  log10(oc$Fluence) - 7)
put("PhotonIndex",   oc$PhotonIndex)
put("log10NH",       log10(oc$NH) + 21)
put("log10PeakFlux", log10(ifelse(oc$PeakFlux > 0, oc$PeakFlux, NA)))
put("Gamma",         rep(NA_real_, sum(is_opt)))
new$log10NH[new$log10NH < 20] <- NA
new$PhotonIndex[new$PhotonIndex < 0] <- NA
new$log10PeakFlux[is.infinite(new$log10PeakFlux)] <- NA

set.seed(1)
imp <- mice(new[, linear_vars], m = 20, method = "midastouch", printFlag = FALSE)
COMPS <- c(20L, 10L, 5L)   # completion 20 == round-4 frame
frames <- lapply(COMPS, function(cc) {
  f <- new; f[, linear_vars] <- complete(imp, cc)
  for (v in linear_vars) f[[paste0(v, "Sqr")]] <- f[[v]]^2
  f$is_optical <- as.numeric(is_opt)
  f
})

Y <- old$log10z; N <- nrow(old)
N_REP <- 5L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

fg <- read.table("Best_formula_GAM.txt"); gam_forms <- apply(as.matrix(fg[, 2]), 1, as.formula)
fl <- read.table("Best_formula_GLM.txt"); glm_forms <- apply(as.matrix(fl[, 2]), 1, as.formula)
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1:6]), detailed_names = FALSE, name_prefix = "gamF")
glmL <- create.Learner("SL.custom_glm", tune = list(glm.model = glm_forms[1:4]), detailed_names = FALSE, name_prefix = "glmF")
PAPER_LIBS <- c(gamL$names, glmL$names, "SL.randomForest", "SL.mean")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")

run_rep <- function(rp) {
  folds <- fold_sets[[rp]]
  pred_by_comp <- matrix(NA_real_, nrow = N, ncol = length(frames))
  set.seed(70000 + rp)
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    for (fc in seq_along(frames)) {
      X <- frames[[fc]][, VARS]
      fit <- try(suppressWarnings(SuperLearner(
        Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
        SL.library = PAPER_LIBS, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
      if (!inherits(fit, "try-error")) pred_by_comp[te, fc] <- as.numeric(fit$SL.predict)
    }
  }
  data.frame(rep = rp, row = seq_len(N), is_opt = is_opt, y = Y,
             pred_single = pred_by_comp[, 1],
             pred_avg3   = rowMeans(pred_by_comp, na.rm = TRUE))
}

res <- parallel::mclapply(seq_len(N_REP), run_rep, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 5: single MICE completion vs 3-completion average ==\n")
for (col in c("pred_single", "pred_avg3")) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$rep == rp & is.finite(res[[col]]), ]
    cor(d$y, d[[col]])
  })
  cat(sprintf("%-12s r_all=%.4f (per-rep: %s)\n", col, mean(rr), paste(sprintf("%.3f", rr), collapse = " ")))
}
