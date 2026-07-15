# cv_round16.R — untied-variance Daumé (Finkel-Manning column scaling) + V=20 screen.
# Round-10 semi frame, lean library, 5f x 3 reps.
# Configs: c=1 (anchor, current daume) | c=0.5 | c=0.25 | c=0.1 | c2 (=2) | v20 (c=1, inner V=20)
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

daume_aug <- function(X, c_scale = 1) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical * c_scale
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

CFG <- list(anchor = list(c = 1), c05 = list(c = 0.5), c025 = list(c = 0.25),
            c01 = list(c = 0.1), c2 = list(c = 2), v20 = list(c = 1, V = 20))

N_REP <- 3L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  X <- daume_aug(f0[, VARS], cfg$c)
  V <- if (!is.null(cfg$V)) cfg$V else 5
  pred <- rep(NA_real_, N)
  set.seed(190000 + 100 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = LEAN, cvControl = list(V = V), verbose = FALSE)), silent = TRUE)
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

cat("\n== round 16 untied-Daumé + V20 screen ==\n")
base <- NULL
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    c(all = cor(d$y, d$pred), op = cor(d$y[d$is_opt], d$pred[d$is_opt]))
  })
  if (cn == "anchor") base <- rr
  cat(sprintf("%-7s r_all=%.4f (d=%+.4f)  r_opt=%.4f\n", cn,
              mean(rr["all", ]), mean(rr["all", ] - base["all", ]), mean(rr["op", ])))
}
