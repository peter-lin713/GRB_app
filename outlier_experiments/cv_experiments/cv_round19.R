# cv_round19.R — final hypothesis screen: imputation engines + robust-loss learners.
#   anchor    round-10 semi-MICE frame, lean lib, daume
#   sbgcop    5 Gaussian-copula posterior completions (co-imputed w/ unlabeled), pred-averaged
#   missf     missForest co-imputed frame (deterministic, single)
#   huber     anchor + gbm Huber-loss + gbm median (quantile 0.5) learners
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(mice)
  library(sbgcop); library(missForest); library(gbm)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

SC <- "/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad"
out_prefix <- commandArgs(trailingOnly = TRUE)[1]

f0 <- read.csv(file.path(SC, "round10_frame.csv"), row.names = 1, stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")
Y <- f0$log10z; N <- nrow(f0); is_opt <- f0$is_optical == 1

# ---- rebuild pre-imputation corrected features + unlabeled pool (as round 10) ----
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
opt <- read.csv("Data/optical_data.csv", stringsAsFactors = FALSE); names(opt)[1] <- "GRB"
nmiss <- rowSums(is.na(raw[, linear_vars])); io <- nmiss >= 6
mo <- match(raw$GRB[io], opt$GRB); oc <- opt[mo, ]
raw$log10T90[io] <- log10(oc$T90); raw$log10Fluence[io] <- log10(oc$Fluence) - 7
raw$PhotonIndex[io] <- oc$PhotonIndex; raw$log10NH[io] <- log10(oc$NH) + 21
raw$log10PeakFlux[io] <- log10(ifelse(oc$PeakFlux > 0, oc$PeakFlux, NA))
raw$Gamma[io] <- NA_real_
raw <- raw[is.na(raw$log10T90) | raw$log10T90 > 0.301, ]
raw$log10PeakFlux[is.infinite(raw$log10PeakFlux)] <- NA
raw$log10NH[raw$log10NH < 20] <- NA
raw$Beta[raw$Beta > 3] <- NA; raw$Gamma[raw$Gamma > 3] <- NA; raw$Alpha[raw$Alpha > 3] <- NA
raw$PhotonIndex[raw$PhotonIndex < 0] <- NA
rownames(raw) <- raw$GRB
gen <- read.csv("Data/TOTAL_GENERALIZATION_DATA_v4.csv", row.names = 1, stringsAsFactors = FALSE)
gu <- data.frame(
  log10PeakFlux = suppressWarnings(as.numeric(gen$logPeakFlux)),
  log10NH = gen$logNH, log10Ta = gen$T_abest,
  PhotonIndex = suppressWarnings(as.numeric(gen$photon_index)),
  log10Fa = gen$Fbest, Alpha = gen$Alpha, log10T90 = gen$logT90,
  Beta = NA_real_, Gamma = gen$Gamma, log10Fluence = NA_real_)
gu$log10NH[gu$log10NH < 20] <- NA; gu$PhotonIndex[gu$PhotonIndex < 0] <- NA
gu$Gamma[gu$Gamma > 3] <- NA
tr_idx <- match(f0$GRB, raw$GRB)   # rows of raw matching our (post-cut) frame
comb <- rbind(raw[, linear_vars], gu)

mk_frame <- function(vals) {
  f <- f0
  f[, linear_vars] <- vals[tr_idx, ]
  for (v in linear_vars) f[[paste0(v, "Sqr")]] <- f[[v]]^2
  f
}

cat("Fitting sbgcop MCMC...\n")
set.seed(11)
sb <- sbgcop::sbgcop.mcmc(as.matrix(comb), nsamp = 2500, odens = 100, verb = FALSE)
sb_ids <- round(seq(10, dim(sb$Y.impute)[3], length.out = 5))
sb_frames <- lapply(sb_ids, function(i) mk_frame(sb$Y.impute[, , i]))
cat("sbgcop frames ready\n")

cat("Fitting missForest...\n")
set.seed(12)
mf <- missForest::missForest(comb)$ximp
mf_frame <- list(mk_frame(mf))
cat("missForest frame ready\n")

fg <- read.table("Best_formula_GAM.txt"); gam_forms <- unique(apply(as.matrix(fg[, 2]), 1, as.formula))
band_form <- as.formula(paste(
  "Response ~ (log10Fa + log10Ta + log10NH + PhotonIndex + log10PeakFlux)^2 +",
  "Alpha + Beta + log10T90 + Gamma + log10Fluence +",
  "is_optical:(log10Fa + log10Ta + Alpha + Beta) + is_optical"))
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1]), detailed_names = FALSE, name_prefix = "gamF")
glmB <- create.Learner("SL.custom_glm", tune = list(glm.model = list(band_form)), detailed_names = FALSE, name_prefix = "glmB")
SL.glm_all <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- glm(Y ~ ., data = df, family = family, weights = obsWeights)
  out <- list(object = fit); class(out) <- "SL.glm_all"
  list(pred = as.numeric(predict(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.glm_all <- function(object, newdata, ...) as.numeric(predict(object$object, newdata = newdata, type = "response"))
mk_gbm <- function(dist) {
  function(Y, X, newX, family, obsWeights, ...) {
    df <- data.frame(Y = Y, X, check.names = FALSE)
    fit <- gbm::gbm(Y ~ ., data = df, distribution = dist, n.trees = 500,
                    interaction.depth = 2, shrinkage = 0.05, bag.fraction = 0.8,
                    n.minobsinnode = 10, verbose = FALSE)
    out <- list(object = fit); class(out) <- paste0("SL.gbm_", if (is.character(dist)) dist else "q50")
    list(pred = as.numeric(predict(fit, newdata = data.frame(newX, check.names = FALSE), n.trees = 500)), fit = out)
  }
}
SL.gbm_huber <- mk_gbm(list(name = "huberized"))
SL.gbm_q50   <- mk_gbm(list(name = "quantile", alpha = 0.5))
predict.SL.gbm_huberized <- predict.SL.gbm_q50 <- function(object, newdata, ...)
  as.numeric(predict(object$object, newdata = data.frame(newdata, check.names = FALSE), n.trees = 500))

LEAN <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.randomForest", "SL.mean")

daume_aug <- function(X) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

N_REP <- 3L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# grid: (rep, cfg, sub) where sub indexes frames within a config
grid <- rbind(
  expand.grid(rep = 1:N_REP, cfg = "anchor", sub = 1),
  expand.grid(rep = 1:N_REP, cfg = "sbgcop", sub = 1:5),
  expand.grid(rep = 1:N_REP, cfg = "missf",  sub = 1),
  expand.grid(rep = 1:N_REP, cfg = "huber",  sub = 1)
)
tasks <- split(grid, seq_len(nrow(grid)))

run_task <- function(t) {
  f <- switch(as.character(t$cfg),
    anchor = f0, huber = f0,
    sbgcop = sb_frames[[t$sub]],
    missf  = mf_frame[[1]])
  libs <- if (t$cfg == "huber") c(LEAN, "SL.gbm_huber", "SL.gbm_q50") else LEAN
  X <- daume_aug(f[, VARS])
  folds <- fold_sets[[t$rep]]
  pred <- rep(NA_real_, N)
  set.seed(220000 + 1000 * t$rep + t$sub + 100 * match(as.character(t$cfg), c("anchor","sbgcop","missf","huber")))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = libs, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", as.character(t$cfg), t$rep, t$sub, i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(rep = t$rep, cfg = as.character(t$cfg), sub = t$sub,
             row = seq_len(N), is_opt = is_opt, y = Y, pred = pred)
}

cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 19 imputation engines + robust losses ==\n")
for (cn in c("anchor", "sbgcop", "missf", "huber")) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp, ]
    pa <- tapply(d$pred, d$row, mean, na.rm = TRUE)   # averages over subs for sbgcop
    yy <- tapply(d$y, d$row, mean)
    cor(yy, pa)
  })
  cat(sprintf("%-7s r_all=%.4f (per-rep: %s)\n", cn, mean(rr), paste(sprintf("%.3f", rr), collapse = " ")))
}
