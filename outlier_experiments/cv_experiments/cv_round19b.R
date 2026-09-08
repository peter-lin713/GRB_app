# cv_round19b.R — lit-rev-4 try-verdicts on the missForest frame (new best):
#   anchor   missForest frame, lean lib, daume
#   quadEN   + SL.glmnet.quad (full pairwise-product basis, Daumé AFTER expansion,
#              alpha 0 / 0.25 as two candidates)
#   ridgeSt  anchor library but nonneg-ridge metalearner instead of NNLS
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(missForest); library(nnls)
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

# rebuild missForest frame (same recipe as round 19, seed 12)
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
comb <- rbind(raw[, linear_vars], gu)
set.seed(12)
mf <- missForest::missForest(comb)$ximp
fMF <- f0
fMF[, linear_vars] <- mf[tr_idx, ]
for (v in linear_vars) fMF[[paste0(v, "Sqr")]] <- fMF[[v]]^2
write.csv(fMF, file.path(SC, "missforest_frame.csv"))

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

# full-quadratic EN with Daumé applied AFTER expansion; operates on the linear vars only
quad_expand <- function(D) {
  mm <- model.matrix(~ .^2, data = as.data.frame(D))[, -1, drop = FALSE]
  cbind(mm, as.matrix(D)^2)
}
mk_quad <- function(alpha) {
  force(alpha)
  function(Y, X, newX, family, obsWeights, ...) {
    dom  <- X$is_optical
    Q    <- quad_expand(X[, linear_vars]);   Qn <- quad_expand(newX[, linear_vars])
    Xq   <- cbind(Q, Q * dom, is_optical = dom)
    Xqn  <- cbind(Qn, Qn * newX$is_optical, is_optical = newX$is_optical)
    fit  <- glmnet::cv.glmnet(as.matrix(Xq), Y, alpha = alpha, nfolds = 10)
    out  <- list(object = fit); class(out) <- if (alpha == 0) "SL.quad0" else "SL.quad25"
    list(pred = as.numeric(predict(fit, as.matrix(Xqn), s = "lambda.min")), fit = out)
  }
}
SL.quad0  <- mk_quad(0)
SL.quad25 <- mk_quad(0.25)
predict.SL.quad0 <- predict.SL.quad25 <- function(object, newdata, ...) {
  dom <- newdata$is_optical
  Qn  <- quad_expand(newdata[, linear_vars])
  as.numeric(predict(object$object, as.matrix(cbind(Qn, Qn * dom, is_optical = dom)), s = "lambda.min"))
}

# nonneg-ridge metalearner
method.NNRidge <- function() {
  out <- SuperLearner::method.NNLS()
  out$computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights, ...) {
    fit <- glmnet::cv.glmnet(Z, Y, alpha = 0, lower.limits = 0, nfolds = 10)
    co <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    co[co < 0] <- 0
    if (sum(co) > 0) co <- co / sum(co) else co <- rep(1/ncol(Z), ncol(Z))
    cvRisk <- apply(Z, 2, function(p) mean((p - Y)^2))
    list(cvRisk = cvRisk, coef = co, optimizer = fit)
  }
  out
}

LEAN <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.randomForest", "SL.mean")
daume_aug <- function(X) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

CFG <- list(
  anchor  = list(libs = LEAN, method = "method.NNLS"),
  quadEN  = list(libs = c(LEAN, "SL.quad0", "SL.quad25"), method = "method.NNLS"),
  ridgeSt = list(libs = LEAN, method = method.NNRidge())
)

N_REP <- 3L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  # quad learners need raw linear vars + is_optical (they self-expand); others get daume aug
  Xbase <- fMF[, VARS]
  X <- daume_aug(Xbase)
  X[, linear_vars] <- Xbase[, linear_vars]   # keep originals available for quad learners
  pred <- rep(NA_real_, N)
  set.seed(230000 + 100 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = cfg$libs, method = cfg$method,
      cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", task$cfg, task$rep, i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(cfg = task$cfg, rep = task$rep, row = seq_len(N), is_opt = is_opt, y = Y, pred = pred)
}

tasks <- do.call(rbind, lapply(names(CFG), function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 19b quad-EN + ridge stacking (missForest frame) ==\n")
base <- NULL
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    cor(d$y, d$pred)
  })
  if (cn == "anchor") base <- rr
  cat(sprintf("%-8s r_all=%.4f (d=%+.4f)\n", cn, mean(rr), mean(rr - base)))
}
