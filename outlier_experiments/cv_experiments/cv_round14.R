# cv_round14.R — formula regeneration with honest per-fold selection (lean screen,
# round-10 semi frame, daume applied everywhere):
#   anchor | +bestglm (150-formula pool, inner-CV select) | +bestglm+bestgam (20 GAM pool)
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

SC <- "/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad"
out_prefix <- commandArgs(trailingOnly = TRUE)[1]

f0 <- read.csv(file.path(SC, "round10_frame.csv"), row.names = 1, stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
top5 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")
Y <- f0$log10z; N <- nrow(f0); is_opt <- f0$is_optical == 1

# ---- candidate formula pools (structure fixed a priori, seed-deterministic) ----
pairs5 <- combn(top5, 2)
inter_terms <- apply(pairs5, 2, paste, collapse = ":")
term_pool <- c(linear_vars, paste0(top7, "Sqr"), "is_optical",
               paste0("is_optical:", c("log10Fa","log10Ta","Alpha","Beta")), inter_terms)
set.seed(777)
glm_pool <- lapply(seq_len(150), function(i) {
  k <- sample(8:18, 1)
  reformulate(sample(term_pool, k), response = "Y")
})
gam_pool <- lapply(seq_len(20), function(i) {
  ns <- sample(2:4, 1)
  sv <- sample(top5, ns)
  lin <- setdiff(sample(linear_vars, sample(4:8, 1)), sv)
  as.formula(paste("Y ~", paste(c(sprintf("s(%s, k=5)", sv), lin, "is_optical"), collapse = "+")))
})

inner_cv_mse <- function(form, df, fitter, predicter, V = 5) {
  n <- nrow(df); fold <- sample(rep(1:V, length.out = n))
  se <- 0
  for (v in 1:V) {
    tr <- df[fold != v, ]; te <- df[fold == v, ]
    fit <- try(suppressWarnings(fitter(form, tr)), silent = TRUE)
    if (inherits(fit, "try-error")) return(Inf)
    p <- try(suppressWarnings(predicter(fit, te)), silent = TRUE)
    if (inherits(p, "try-error") || anyNA(p)) return(Inf)
    se <- se + sum((te$Y - p)^2)
  }
  se / n
}

SL.bestglm <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  mses <- vapply(glm_pool, inner_cv_mse, numeric(1), df = df,
                 fitter = function(f, d) glm(f, data = d, family = gaussian()),
                 predicter = function(fit, d) as.numeric(predict(fit, newdata = d)))
  best <- glm_pool[[which.min(mses)]]
  fit <- glm(best, data = df, family = gaussian())
  out <- list(object = fit); class(out) <- "SL.bestglm"
  list(pred = as.numeric(predict(fit, newdata = data.frame(newX, check.names = FALSE))), fit = out)
}
predict.SL.bestglm <- function(object, newdata, ...)
  as.numeric(predict(object$object, newdata = data.frame(newdata, check.names = FALSE)))

SL.bestgam <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  mses <- vapply(gam_pool, inner_cv_mse, numeric(1), df = df,
                 fitter = function(f, d) mgcv::gam(f, data = d, method = "REML"),
                 predicter = function(fit, d) as.numeric(mgcv::predict.gam(fit, newdata = d)),
                 V = 3)
  best <- gam_pool[[which.min(mses)]]
  fit <- mgcv::gam(best, data = df, method = "REML")
  out <- list(object = fit); class(out) <- "SL.bestgam"
  list(pred = as.numeric(mgcv::predict.gam(fit, newdata = data.frame(newX, check.names = FALSE))), fit = out)
}
predict.SL.bestgam <- function(object, newdata, ...)
  as.numeric(mgcv::predict.gam(object$object, newdata = data.frame(newdata, check.names = FALSE)))

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
LEAN <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.randomForest", "SL.mean")

daume_aug <- function(X) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

CFG <- list(
  anchor  = LEAN,
  bglm    = c(LEAN, "SL.bestglm"),
  bglmgam = c(LEAN, "SL.bestglm", "SL.bestgam")
)

N_REP <- 3L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

run_task <- function(task) {
  folds <- fold_sets[[task$rep]]
  X <- daume_aug(f0[, VARS])
  pred <- rep(NA_real_, N)
  set.seed(170000 + 100 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = CFG[[task$cfg]], cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
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

cat("\n== round 14 formula regeneration (lean screen, daume everywhere) ==\n")
base <- NULL
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    cor(d$y, d$pred)
  })
  if (cn == "anchor") base <- rr
  cat(sprintf("%-8s r_all=%.4f (d=%+.4f)\n", cn, mean(rr), mean(rr - base)))
}
