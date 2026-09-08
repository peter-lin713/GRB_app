# cv_round7.R — headline run: corrected pipeline end-to-end.
# Raw emcee frame -> recover real optical-catalog features -> clean -> MICE(m=20)
# -> fresh 5% M-estimator cut computed on CORRECTED features -> 10-fold CV x N reps,
# predictions averaged over 3 MICE completions, expanded learner library.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(MASS); library(mgcv); library(gam)
  library(kernlab); library(randomForest); library(mice)
})
source("Custom_SL/sl_xgboost_safe.R")
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

args <- commandArgs(trailingOnly = TRUE)
out_prefix <- args[1]
N_REP <- if (length(args) >= 2) as.integer(args[2]) else 3L
K     <- 10L

raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
opt <- read.csv("Data/optical_data.csv", stringsAsFactors = FALSE); names(opt)[1] <- "GRB"
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")

# recover real optical-catalog features (units validated on 88 overlap GRBs)
nmiss <- rowSums(is.na(raw[, linear_vars]))
io <- nmiss >= 6
mo <- match(raw$GRB[io], opt$GRB)
stopifnot(!anyNA(mo))
oc <- opt[mo, ]
raw$log10T90[io]      <- log10(oc$T90)
raw$log10Fluence[io]  <- log10(oc$Fluence) - 7
raw$PhotonIndex[io]   <- oc$PhotonIndex
raw$log10NH[io]       <- log10(oc$NH) + 21
raw$log10PeakFlux[io] <- log10(ifelse(oc$PeakFlux > 0, oc$PeakFlux, NA))
raw$Gamma[io]         <- NA_real_
raw$is_optical <- as.numeric(io)

# pipeline cleaning + T90 cut (now with real T90 for optical rows)
raw <- raw[is.na(raw$log10T90) | raw$log10T90 > 0.301, ]
raw$log10PeakFlux[is.infinite(raw$log10PeakFlux)] <- NA
raw$log10NH[raw$log10NH < 20] <- NA
raw$Beta[raw$Beta > 3] <- NA
raw$Gamma[raw$Gamma > 3] <- NA
raw$Alpha[raw$Alpha > 3] <- NA
raw$PhotonIndex[raw$PhotonIndex < 0] <- NA
rownames(raw) <- raw$GRB
cat("rows after T90 cut:", nrow(raw), "( optical:", sum(raw$is_optical), ")\n")

set.seed(1)
imp <- mice(raw[, linear_vars], m = 20, method = "midastouch", printFlag = FALSE)
COMPS <- c(20L, 10L, 5L)
frames <- lapply(COMPS, function(cc) {
  f <- raw
  f[, linear_vars] <- complete(imp, cc)
  for (v in linear_vars) f[[paste0(v, "Sqr")]] <- f[[v]]^2
  f$log10z <- log10(f$Redshift_crosscheck + 1)
  f
})

# fresh 5% M-estimator cut on corrected features (completion 20, pipeline logic)
f0 <- frames[[1]]
pred_cols <- c(linear_vars, paste0(linear_vars, "Sqr"))
mm <- model.matrix(reformulate(pred_cols), data = f0)
qrr <- qr(mm)
indep <- setdiff(colnames(mm)[qrr$pivot[seq_len(qrr$rank)]], "(Intercept)")
M_est <- MASS::rlm(reformulate(indep, response = "log10z"), data = f0, method = "M")
thr <- quantile(M_est$w, 0.05)
keep <- M_est$w > thr
cat("M-est cut drops", sum(!keep), "rows; kept:", sum(keep),
    "( optical dropped:", sum(!keep & f0$is_optical == 1), ")\n")
frames <- lapply(frames, function(f) f[keep, ])
f0 <- frames[[1]]

Y <- f0$log10z; Z <- f0$Redshift_crosscheck; N <- nrow(f0)
is_opt <- f0$is_optical == 1
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")

set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# ---- expanded library (round 6) ----
fg <- read.table("Best_formula_GAM.txt"); gam_forms <- unique(apply(as.matrix(fg[, 2]), 1, as.formula))
fl <- read.table("Best_formula_GLM.txt"); glm_forms <- unique(apply(as.matrix(fl[, 2]), 1, as.formula))
smooth_forms <- list(
  as.formula("Response ~ s(log10PeakFlux) + s(log10NH) + s(log10Ta) + s(log10Fa) + s(PhotonIndex) + Alpha + log10T90 + Beta + Gamma + log10Fluence + is_optical"),
  as.formula("Response ~ te(log10Fa, log10Ta) + s(log10NH) + s(log10PeakFlux) + s(PhotonIndex) + Alpha + Gamma + is_optical"),
  as.formula("Response ~ s(log10PeakFlux) + s(log10NH) + log10Ta + log10Fa + PhotonIndex + Alpha + Gamma + is_optical")
)
band_form <- as.formula(paste(
  "Response ~ (log10Fa + log10Ta + log10NH + PhotonIndex + log10PeakFlux)^2 +",
  "Alpha + Beta + log10T90 + Gamma + log10Fluence +",
  "is_optical:(log10Fa + log10Ta + Alpha + Beta) + is_optical"))
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms),    detailed_names = FALSE, name_prefix = "gamF")
gamS <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = smooth_forms), detailed_names = FALSE, name_prefix = "gamS")
glmL <- create.Learner("SL.custom_glm", tune = list(glm.model = glm_forms),    detailed_names = FALSE, name_prefix = "glmF")
glmB <- create.Learner("SL.custom_glm", tune = list(glm.model = list(band_form)), detailed_names = FALSE, name_prefix = "glmB")
SL.mgcv_sel <- function(Y, X, newX, family, obsWeights, ...) {
  sm <- intersect(c("log10PeakFlux","log10NH","log10Ta","log10Fa","PhotonIndex","Alpha","log10T90","Beta","Gamma","log10Fluence"), colnames(X))
  rest <- setdiff(colnames(X), sm)
  f <- as.formula(paste("Y ~", paste(c(sprintf("s(%s, k=5)", sm), rest), collapse = "+")))
  fit <- mgcv::gam(f, data = X, family = family, weights = obsWeights, select = TRUE, method = "REML")
  out <- list(object = fit); class(out) <- "SL.mgcv_sel"
  list(pred = as.numeric(mgcv::predict.gam(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.mgcv_sel <- function(object, newdata, ...)
  as.numeric(mgcv::predict.gam(object$object, newdata = newdata, type = "response"))
SL.gausspr <- function(Y, X, newX, family, obsWeights, ...) {
  fit <- kernlab::gausspr(as.matrix(X), Y, kernel = "rbfdot", var = 0.05)
  out <- list(object = fit); class(out) <- "SL.gausspr"
  list(pred = as.numeric(kernlab::predict(fit, as.matrix(newX))), fit = out)
}
predict.SL.gausspr <- function(object, newdata, ...) as.numeric(kernlab::predict(object$object, as.matrix(newdata)))
SL.glm_all <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- glm(Y ~ ., data = df, family = family, weights = obsWeights)
  out <- list(object = fit); class(out) <- "SL.glm_all"
  list(pred = as.numeric(predict(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.glm_all <- function(object, newdata, ...) as.numeric(predict(object$object, newdata = newdata, type = "response"))
LIBS <- c(gamL$names, gamS$names, "SL.mgcv_sel", glmL$names, glmB$names,
          "SL.glm_all", "SL.gausspr", "SL.glmnet", "SL.randomForest", "SL.mean")

run_task <- function(t) {   # t = (rep, comp)
  folds <- fold_sets[[t$rep]]
  X <- frames[[t$comp]][keeprows <- seq_len(N), VARS]
  pred <- rep(NA_real_, N)
  set.seed(91000 + 100 * t$rep + t$comp)
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = LIBS, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR rep", t$rep, "comp", t$comp, "fold", i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(rep = t$rep, comp = t$comp, row = seq_len(N), is_opt = is_opt, y = Y, z = Z, pred = pred)
}

tasks <- expand.grid(rep = seq_len(N_REP), comp = seq_along(COMPS))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks (", N_REP, "reps x", length(COMPS), "completions, K =", K, ")\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 7 headline (corrected pipeline, 10-fold, MICE-avg3) ==\n")
for (rp in seq_len(N_REP)) {
  d <- res[res$rep == rp, ]
  pa <- tapply(d$pred, d$row, mean, na.rm = TRUE)          # avg over completions
  yy <- tapply(d$y, d$row, mean); oo <- tapply(d$is_opt, d$row, mean) == 1
  cat(sprintf("rep %d: r_all=%.4f  r_xray=%.4f  r_opt=%.4f\n",
              rp, cor(yy, pa), cor(yy[!oo], pa[!oo]), cor(yy[oo], pa[oo])))
}
