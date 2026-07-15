# cv_round13.R — interaction ablation on the round-10 semi-MICE frame (lean screen):
#   plain | daume | coral | daume+coral | daume+iso | daume+coral+iso
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(FNN)
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

N_REP <- 3L; K <- 5L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

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
shrink_cov <- function(X, lam = 0.2) { S <- cov(X); (1 - lam) * S + lam * diag(diag(S)) }
coral_align <- function(Xtr, Xte, dom_tr, dom_te) {
  xr <- Xtr[!dom_tr, linear_vars]; op <- Xtr[dom_tr, linear_vars]
  if (sum(dom_tr) < 15) return(list(tr = Xtr, te = Xte))
  Cs <- shrink_cov(op); Ct <- shrink_cov(xr)
  A <- solve(chol(Cs)) %*% chol(Ct)
  mu_s <- colMeans(op); mu_t <- colMeans(xr)
  fix <- function(X, dom) {
    if (!any(dom)) return(X)
    Xo <- sweep(as.matrix(X[dom, linear_vars]), 2, mu_s) %*% A
    Xo <- sweep(Xo, 2, mu_t, "+")
    X[dom, linear_vars] <- Xo
    for (v in top7) X[[paste0(v, "Sqr")]][dom] <- X[[v]][dom]^2
    X
  }
  list(tr = fix(Xtr, dom_tr), te = fix(Xte, dom_te))
}
add_iso <- function(Xtr, Xte) {
  sc <- scale(Xtr[, linear_vars])
  Xtr$outlier_score <- FNN::knn.dist(sc, k = 10)[, 10]
  Xte$outlier_score <- FNN::get.knnx(sc, scale(Xte[, linear_vars],
    center = attr(sc, "scaled:center"), scale = attr(sc, "scaled:scale")), k = 10)$nn.dist[, 10]
  list(tr = Xtr, te = Xte)
}

CFG <- list(
  plain      = c(),
  daume      = c("daume"),
  coral      = c("coral"),
  dm_coral   = c("coral", "daume"),
  dm_iso     = c("daume", "iso"),
  dm_cor_iso = c("coral", "daume", "iso")
)

run_task <- function(task) {
  mods <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  X <- f0[, VARS]
  pred <- rep(NA_real_, N)
  set.seed(160000 + 100 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- X[tr, ]; Xte <- X[te, ]
    if ("coral" %in% mods) { al <- coral_align(Xtr, Xte, is_opt[tr], is_opt[te]); Xtr <- al$tr; Xte <- al$te }
    if ("iso"   %in% mods) { ai <- add_iso(Xtr, Xte); Xtr <- ai$tr; Xte <- ai$te }
    if ("daume" %in% mods) { Xtr <- daume_aug(Xtr); Xte <- daume_aug(Xte) }
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = Xtr, newX = Xte, family = gaussian(),
      SL.library = LEAN, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
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

cat("\n== round 13 ablation (semi frame, 5f x 3 reps, lean) ==\n")
base <- NULL
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    c(all = cor(d$y, d$pred), op = cor(d$y[d$is_opt], d$pred[d$is_opt]))
  })
  if (cn == "plain") base <- rr
  cat(sprintf("%-11s r_all=%.4f (d=%+.4f)  r_opt=%.4f\n", cn,
              mean(rr["all", ]), mean(rr["all", ] - base["all", ]), mean(rr["op", ])))
}
