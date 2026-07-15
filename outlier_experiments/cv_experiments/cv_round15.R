# cv_round15.R — protocol/cleaning screen on round-10 semi frame (daume everywhere):
#   anchor | infold (extra in-fold rlm training-row filter) | v10 (inner V=10)
#   | avg10 (10-completion MICE averaging — needs frames; approximated via
#     completions of the cached mice run: here we re-impute quickly with m=10)
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
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")
Y <- f0$log10z; N <- nrow(f0); is_opt <- f0$is_optical == 1

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

CFG <- c("anchor", "infold", "v10")
N_REP <- 3L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

run_task <- function(task) {
  folds <- fold_sets[[task$rep]]
  X <- daume_aug(f0[, VARS])
  pred <- rep(NA_real_, N)
  V <- if (task$cfg == "v10") 10 else 5
  set.seed(180000 + 100 * task$rep + match(task$cfg, CFG))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    if (task$cfg == "infold") {
      dfr <- data.frame(Y = Y[tr], f0[tr, c(linear_vars, paste0(top7, "Sqr"))])
      rf <- tryCatch(MASS::rlm(Y ~ ., data = dfr, method = "M", maxit = 50), error = function(e) NULL)
      if (!is.null(rf)) {
        r <- residuals(rf)
        keep <- abs(r) <= 3 * mad(r)
        if (sum(keep) >= 30) tr <- tr[keep]
      }
    }
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = LEAN, cvControl = list(V = V), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", task$cfg, task$rep, i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(cfg = task$cfg, rep = task$rep, row = seq_len(N), is_opt = is_opt, y = Y, pred = pred)
}

tasks <- do.call(rbind, lapply(CFG, function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 15 protocol screen ==\n")
base <- NULL
for (cn in CFG) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    cor(d$y, d$pred)
  })
  if (cn == "anchor") base <- rr
  cat(sprintf("%-7s r_all=%.4f (d=%+.4f)\n", cn, mean(rr), mean(rr - base)))
}
