# cv_round11.R — SCREENING tier (5-fold x 3 reps, lean library, cached frame):
#   anchor    lean lib on corrected frame (round9_frame.csv)
#   coral     CORAL-align optical rows to X-ray covariance (fold-internal, shrunk cov)
#   wopt      optical obsWeights 0.5 (on corrected frame — retest post-correction)
#   coreg     COREG-style pseudo-labeling of unlabeled GRBs via 2 kNN views, downweighted
#   poolpca   +4 PCA components fit on pooled labeled+unlabeled (self-taught, y-blind)
#   isoscore  +isolation-proxy outlier score feature (kNN distance based)
#   knnenc    +OOF kNN-mean-target encoding (leakage-safe: encoder built per fold)
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(FNN)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

SC <- "/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad"
out_prefix <- commandArgs(trailingOnly = TRUE)[1]

f0 <- read.csv(file.path(SC, "round9_frame.csv"), row.names = 1, stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")
Y <- f0$log10z; N <- nrow(f0); is_opt <- f0$is_optical == 1

# unlabeled pool (for coreg/poolpca): impute simply (median) — embeddings/pseudo-labels only
gen <- read.csv("Data/TOTAL_GENERALIZATION_DATA_v4.csv", row.names = 1, stringsAsFactors = FALSE)
gu <- data.frame(
  log10PeakFlux = suppressWarnings(as.numeric(gen$logPeakFlux)),
  log10NH = gen$logNH, log10Ta = gen$T_abest,
  PhotonIndex = suppressWarnings(as.numeric(gen$photon_index)),
  log10Fa = gen$Fbest, Alpha = gen$Alpha, log10T90 = gen$logT90,
  Beta = NA_real_, Gamma = gen$Gamma, log10Fluence = NA_real_)
for (v in linear_vars) gu[[v]][!is.finite(gu[[v]])] <- median(f0[[v]])
gu <- gu[, linear_vars]

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

shrink_cov <- function(X, lam = 0.2) {
  S <- cov(X); (1 - lam) * S + lam * diag(diag(S))
}
coral_align <- function(Xtr, Xte, dom_tr, dom_te) {
  # align optical rows to X-ray covariance; transform fit on training rows only
  xr <- Xtr[!dom_tr, linear_vars]; op <- Xtr[dom_tr, linear_vars]
  if (sum(dom_tr) < 15) return(list(tr = Xtr, te = Xte))
  Cs <- shrink_cov(op); Ct <- shrink_cov(xr)
  A <- solve(chol(Cs)) %*% chol(Ct)
  mu_s <- colMeans(op); mu_t <- colMeans(xr)
  fix <- function(X, dom) {
    Xo <- sweep(as.matrix(X[dom, linear_vars]), 2, mu_s) %*% A
    Xo <- sweep(Xo, 2, mu_t, "+")
    X[dom, linear_vars] <- Xo
    for (v in top7) X[[paste0(v, "Sqr")]][dom] <- X[[v]][dom]^2
    X
  }
  list(tr = fix(Xtr, dom_tr), te = fix(Xte, dom_te))
}

knn_encode <- function(Xtr, ytr, Xq, k = 15) {
  sc <- scale(Xtr[, linear_vars])
  q  <- scale(Xq[, linear_vars], center = attr(sc, "scaled:center"), scale = attr(sc, "scaled:scale"))
  # OOF for training rows: exclude self via k+1 then drop first neighbor
  nnq <- FNN::get.knnx(sc, q, k = k)$nn.index
  rowMeans(matrix(ytr[nnq], nrow = nrow(q)))
}

CFG <- c("anchor", "coral", "wopt", "coreg", "poolpca", "isoscore", "knnenc")

run_task <- function(task) {
  folds <- fold_sets[[task$rep]]
  X <- f0[, VARS]
  pred <- rep(NA_real_, N)
  set.seed(140000 + 100 * task$rep + match(task$cfg, CFG))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- X[tr, ]; Xte <- X[te, ]; ytr <- Y[tr]
    w <- rep(1, length(tr))
    if (task$cfg == "coral") {
      al <- coral_align(Xtr, Xte, is_opt[tr], is_opt[te])
      Xtr <- al$tr; Xte <- al$te
    } else if (task$cfg == "wopt") {
      w <- ifelse(is_opt[tr], 0.5, 1)
    } else if (task$cfg == "coreg") {
      # two kNN views propose pseudo-labels for unlabeled pool; accept top-50 agreeing
      sc <- scale(Xtr[, linear_vars])
      qu <- scale(gu, center = attr(sc, "scaled:center"), scale = attr(sc, "scaled:scale"))
      v1 <- FNN::knn.reg(sc[, 1:5], qu[, 1:5], ytr, k = 5)$pred     # view 1: top-5 features
      v2 <- FNN::knn.reg(sc, qu, ytr, k = 10)$pred                  # view 2: all features
      agree <- abs(v1 - v2)
      pick <- order(agree)[1:50]
      Xps <- f0[rep(1, 50), VARS]; Xps[, linear_vars] <- gu[pick, ]
      for (v in top7) Xps[[paste0(v, "Sqr")]] <- Xps[[v]]^2
      Xps$is_optical <- 0
      Xtr <- rbind(Xtr, Xps); ytr <- c(ytr, (v1[pick] + v2[pick]) / 2)
      w <- c(w, rep(0.3, 50))
    } else if (task$cfg == "poolpca") {
      pool <- rbind(Xtr[, linear_vars], gu)
      sc <- scale(pool); pc <- prcomp(sc, rank. = 4)
      prj <- function(Xq) predict(pc, scale(Xq[, linear_vars], center = attr(sc, "scaled:center"),
                                            scale = attr(sc, "scaled:scale")))[, 1:4]
      Ptr <- prj(Xtr); Pte <- prj(Xte)
      colnames(Ptr) <- colnames(Pte) <- paste0("PC", 1:4)
      Xtr <- cbind(Xtr, Ptr); Xte <- cbind(Xte, Pte)
    } else if (task$cfg == "isoscore") {
      sc <- scale(Xtr[, linear_vars])
      d_tr <- FNN::knn.dist(sc, k = 10)[, 10]
      d_te <- FNN::get.knnx(sc, scale(Xte[, linear_vars], center = attr(sc, "scaled:center"),
                                      scale = attr(sc, "scaled:scale")), k = 10)$nn.dist[, 10]
      Xtr$outlier_score <- d_tr; Xte$outlier_score <- d_te
    } else if (task$cfg == "knnenc") {
      # leakage-safe: encode training rows via inner split, test rows via full train
      half <- sample(seq_along(tr), floor(length(tr) / 2))
      enc_tr <- rep(NA_real_, length(tr))
      enc_tr[half]  <- knn_encode(Xtr[-half, ], ytr[-half], Xtr[half, ])
      enc_tr[-half] <- knn_encode(Xtr[half, ],  ytr[half],  Xtr[-half, ])
      Xtr$knn_z <- enc_tr
      Xte$knn_z <- knn_encode(Xtr[, setdiff(colnames(Xtr), "knn_z")], ytr, Xte)
    }
    fit <- try(suppressWarnings(SuperLearner(
      Y = ytr, X = Xtr, newX = Xte, family = gaussian(),
      SL.library = LEAN, obsWeights = w / mean(w), cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", task$cfg, task$rep, i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(cfg = task$cfg, rep = task$rep, row = seq_len(N), is_opt = is_opt, y = Y, pred = pred)
}

tasks <- do.call(rbind, lapply(CFG, function(cn)
  data.frame(cfg = cn, rep = seq_len(N_REP), stringsAsFactors = FALSE)))
tasks <- split(tasks, seq_len(nrow(tasks)))
cat("Running", length(tasks), "screening tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 11 screening (5f x 3 reps, lean lib) ==\n")
base <- NULL
for (cn in CFG) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    c(all = cor(d$y, d$pred), op = cor(d$y[d$is_opt], d$pred[d$is_opt]))
  })
  if (cn == "anchor") base <- rr
  cat(sprintf("%-9s r_all=%.4f (d=%+.4f)  r_opt=%.4f\n", cn,
              mean(rr["all", ]), mean(rr["all", ] - base["all", ]), mean(rr["op", ])))
}
