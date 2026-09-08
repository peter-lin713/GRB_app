# cv_round6.R — expanded, deduplicated learner library on the recovered frame:
# distinct formula GAMs/GLMs, new smooth mgcv GAMs (incl. select=TRUE and tensor),
# explicit all-feature linear GLM, band-interaction GLM, Gaussian process, glmnet.
# Same paired folds. Frame = round4 recovered frame (MICE completion 20).
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(ranger)
  library(earth); library(xgboost); library(MASS); library(mgcv); library(gam)
  library(kernlab); library(randomForest)
})
source("Custom_SL/sl_xgboost_safe.R")
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

out_prefix <- commandArgs(trailingOnly = TRUE)[1]

dat <- read.csv("/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad/round4_frame.csv",
                row.names = 1, stringsAsFactors = FALSE)
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
nmiss_lookup <- setNames(rowSums(is.na(raw[, linear_vars])), raw$GRB)
is_opt <- nmiss_lookup[dat$GRB] >= 6
dat$is_optical <- as.numeric(is_opt)
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")

Y <- dat$log10z; N <- nrow(dat)
N_REP <- 5L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# ---- learners ----
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

gamL  <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms),  detailed_names = FALSE, name_prefix = "gamF")
gamS  <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = smooth_forms), detailed_names = FALSE, name_prefix = "gamS")
glmL  <- create.Learner("SL.custom_glm", tune = list(glm.model = glm_forms),  detailed_names = FALSE, name_prefix = "glmF")
glmB  <- create.Learner("SL.custom_glm", tune = list(glm.model = list(band_form)), detailed_names = FALSE, name_prefix = "glmB")

# mgcv with shrinkage smooths on everything (select=TRUE)
SL.mgcv_sel <- function(Y, X, newX, family, obsWeights, ...) {
  sm <- c("log10PeakFlux","log10NH","log10Ta","log10Fa","PhotonIndex","Alpha","log10T90","Beta","Gamma","log10Fluence")
  sm <- intersect(sm, colnames(X))
  rest <- setdiff(colnames(X), sm)
  f <- as.formula(paste("Y ~", paste(c(sprintf("s(%s, k=5)", sm), rest), collapse = "+")))
  fit <- mgcv::gam(f, data = X, family = family, weights = obsWeights, select = TRUE, method = "REML")
  pred <- as.numeric(mgcv::predict.gam(fit, newdata = newX, type = "response"))
  out <- list(object = fit); class(out) <- "SL.mgcv_sel"
  list(pred = pred, fit = out)
}
predict.SL.mgcv_sel <- function(object, newdata, ...)
  as.numeric(mgcv::predict.gam(object$object, newdata = newdata, type = "response"))

SL.gausspr <- function(Y, X, newX, family, obsWeights, ...) {
  fit <- kernlab::gausspr(as.matrix(X), Y, kernel = "rbfdot", var = 0.05)
  pred <- as.numeric(kernlab::predict(fit, as.matrix(newX)))
  out <- list(object = fit); class(out) <- "SL.gausspr"
  list(pred = pred, fit = out)
}
predict.SL.gausspr <- function(object, newdata, ...) as.numeric(kernlab::predict(object$object, as.matrix(newdata)))

# explicit linear GLM on all columns (the accidental round-3/4 winner, now on purpose)
SL.glm_all <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- glm(Y ~ ., data = df, family = family, weights = obsWeights)
  pred <- as.numeric(predict(fit, newdata = newX, type = "response"))
  out <- list(object = fit); class(out) <- "SL.glm_all"
  list(pred = pred, fit = out)
}
predict.SL.glm_all <- function(object, newdata, ...) as.numeric(predict(object$object, newdata = newdata, type = "response"))

LIBS <- c(gamL$names, gamS$names, "SL.mgcv_sel", glmL$names, glmB$names,
          "SL.glm_all", "SL.gausspr", "SL.glmnet", "SL.randomForest", "SL.mean")
cat("library size:", length(LIBS), "\n")

run_rep <- function(rp) {
  folds <- fold_sets[[rp]]
  X <- dat[, VARS]
  pred <- rep(NA_real_, N)
  coefs <- NULL
  set.seed(80000 + rp)
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = LIBS, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR rep", rp, "fold", i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
    coefs <- rbind(coefs, coef(fit))
  }
  list(df = data.frame(rep = rp, row = seq_len(N), is_opt = is_opt, y = Y, pred = pred),
       coefs = colMeans(coefs))
}

res <- parallel::mclapply(seq_len(N_REP), run_rep, mc.cores = 5L, mc.preschedule = FALSE)
df <- do.call(rbind, lapply(res, `[[`, "df"))
write.csv(df, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

rr <- sapply(seq_len(N_REP), function(rp) { d <- df[df$rep == rp & is.finite(df$pred), ]; cor(d$y, d$pred) })
cat(sprintf("\nround 6 expanded library: r_all=%.4f (per-rep: %s)\n",
            mean(rr), paste(sprintf("%.3f", rr), collapse = " ")))
cat("\nmean ensemble weights:\n")
w <- colMeans(do.call(rbind, lapply(res, `[[`, "coefs")))
print(round(sort(w[w > 0.01], decreasing = TRUE), 3))
