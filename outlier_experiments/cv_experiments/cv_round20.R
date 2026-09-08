# cv_round20.R — FINAL consolidation run. All verified winners:
# corrected features + projection-draw bagging (4 param draws) x missForest
# co-imputation (3 seeds) = 12 frames; Daumé augmentation; nonneg-ridge stacking;
# 10-fold x 3 reps; predictions bagged over frames x reps.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(missForest)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

SC <- "/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad"
out_prefix <- commandArgs(trailingOnly = TRUE)[1]

f0 <- read.csv(file.path(SC, "round10_frame.csv"), row.names = 1, stringsAsFactors = FALSE)  # row set/cut
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")
Y <- f0$log10z; Z <- f0$Redshift_crosscheck; N <- nrow(f0); is_opt <- f0$is_optical == 1

# corrected raw features (pre-imputation, with NAs) + unlabeled pool
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
tr_idx <- match(f0$GRB, raw$GRB)
opt_in_raw <- which(io[match(rownames(raw), raw$GRB)])  # optical rows within raw (post-cut order)

# optical-native measurements for projection draws
mo2 <- match(f0$GRB[is_opt], opt$GRB)
onat <- data.frame(logFa = opt$logFa[mo2], logTa = opt$logT_a[mo2],
                   Alpha = opt$Alpha[mo2], Beta = opt$Beta[mo2])
tgt_map <- c(logFa = "log10Fa", logTa = "log10Ta", Alpha = "Alpha", Beta = "Beta")
draws <- read.csv(file.path(SC, "proj_draws12.csv"), stringsAsFactors = FALSE)
DRAWS <- 1:4; SEEDS <- c(12, 77, 301)

cat("Building", length(DRAWS) * length(SEEDS), "frames (draw x missForest seed)...\n")
frames <- list()
for (d in DRAWS) {
  rawd <- raw
  ri <- match(f0$GRB[is_opt], rawd$GRB)
  for (p in names(tgt_map)) {
    dr <- draws[draws$param == p & draws$draw == d - 1, ]
    val <- dr$m * onat[[p]] + dr$b
    tv <- tgt_map[[p]]
    ok <- is.finite(val)
    rawd[[tv]][ri[ok]] <- val[ok]
  }
  comb <- rbind(rawd[, linear_vars], gu)
  for (s in SEEDS) {
    set.seed(s)
    mfi <- missForest::missForest(comb)$ximp
    f <- f0
    f[, linear_vars] <- mfi[tr_idx, ]
    for (v in linear_vars) f[[paste0(v, "Sqr")]] <- f[[v]]^2
    frames[[length(frames) + 1]] <- f
    cat("frame", length(frames), "done\n")
  }
}

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
method.NNRidge <- function() {
  out <- SuperLearner::method.NNLS()
  out$computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights, ...) {
    fit <- glmnet::cv.glmnet(Z, Y, alpha = 0, lower.limits = 0, nfolds = 10)
    co <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    co[co < 0] <- 0
    if (sum(co) > 0) co <- co / sum(co) else co <- rep(1/ncol(Z), ncol(Z))
    list(cvRisk = apply(Z, 2, function(p) mean((p - Y)^2)), coef = co, optimizer = fit)
  }
  out
}
LEAN <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.randomForest", "SL.mean")
daume_aug <- function(X) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

N_REP <- 3L; K <- 10L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

grid <- expand.grid(rep = seq_len(N_REP), fr = seq_along(frames))
tasks <- split(grid, seq_len(nrow(grid)))
cat("Running", length(tasks), "tasks (", N_REP, "reps x", length(frames), "frames, K =", K, ")\n")

run_task <- function(t) {
  X <- daume_aug(frames[[t$fr]][, VARS])
  folds <- fold_sets[[t$rep]]
  pred <- rep(NA_real_, N)
  set.seed(240000 + 1000 * t$rep + t$fr)
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = LEAN, method = method.NNRidge(),
      cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR rep", t$rep, "frame", t$fr, "fold", i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(rep = t$rep, fr = t$fr, row = seq_len(N), is_opt = is_opt, y = Y, z = Z, pred = pred)
}

res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== ROUND 20 FINAL ==\n")
for (rp in seq_len(N_REP)) {
  d <- res[res$rep == rp, ]
  pa <- tapply(d$pred, d$row, mean, na.rm = TRUE)
  yy <- tapply(d$y, d$row, mean)
  cat(sprintf("rep %d (12-frame bag): r_all=%.4f\n", rp, cor(yy, pa)))
}
pa <- tapply(res$pred, res$row, mean, na.rm = TRUE)
yy <- tapply(res$y, res$row, mean); oo <- tapply(res$is_opt, res$row, mean) == 1
zz <- tapply(res$z, res$row, mean)
cat(sprintf("\nFULLY BAGGED (3 reps x 12 frames): r_all=%.4f  r_xray=%.4f  r_opt=%.4f | linear-z r=%.4f  RMSE(log)=%.4f\n",
            cor(yy, pa), cor(yy[!oo], pa[!oo]), cor(yy[oo], pa[oo]),
            cor(zz, 10^pa - 1), sqrt(mean((yy - pa)^2))))
