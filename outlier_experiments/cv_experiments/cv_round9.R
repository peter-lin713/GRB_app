# cv_round9.R — literature round on the CORRECTED frame (round-7 pipeline),
# 10-fold x 3 reps, lean core library, single MICE completion (screening).
# Configs: anchor | daume (domain feature augmentation) | pls (+PLS learner)
#          | noise (measurement-error noise injection) | cmixup | rankgauss
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(mice); library(pls)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

SC <- "/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad"
out_prefix <- commandArgs(trailingOnly = TRUE)[1]

linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
err_map <- c(log10PeakFlux="log10PeakFluxErr", log10NH=NA, log10Ta="log10TaErr",
             PhotonIndex="PhotonIndexErr", log10Fa="log10FaErr", Alpha="AlphaErr",
             log10T90="log10T90Err", Beta="BetaErr", Gamma=NA, log10Fluence="log10FluenceErr")

# ---- rebuild corrected frame exactly as round 7 (deterministic) ----
raw <- read.csv("Data/superlearner_training_emcee_errcut_relative.csv", stringsAsFactors = FALSE)
opt <- read.csv("Data/optical_data.csv", stringsAsFactors = FALSE); names(opt)[1] <- "GRB"
nmiss <- rowSums(is.na(raw[, linear_vars])); io <- nmiss >= 6
mo <- match(raw$GRB[io], opt$GRB); oc <- opt[mo, ]
raw$log10T90[io] <- log10(oc$T90); raw$log10Fluence[io] <- log10(oc$Fluence) - 7
raw$PhotonIndex[io] <- oc$PhotonIndex; raw$log10NH[io] <- log10(oc$NH) + 21
raw$log10PeakFlux[io] <- log10(ifelse(oc$PeakFlux > 0, oc$PeakFlux, NA))
raw$Gamma[io] <- NA_real_; raw$is_optical <- as.numeric(io)
raw <- raw[is.na(raw$log10T90) | raw$log10T90 > 0.301, ]
raw$log10PeakFlux[is.infinite(raw$log10PeakFlux)] <- NA
raw$log10NH[raw$log10NH < 20] <- NA
raw$Beta[raw$Beta > 3] <- NA; raw$Gamma[raw$Gamma > 3] <- NA; raw$Alpha[raw$Alpha > 3] <- NA
raw$PhotonIndex[raw$PhotonIndex < 0] <- NA
rownames(raw) <- raw$GRB
set.seed(1)
imp <- mice(raw[, linear_vars], m = 20, method = "midastouch", printFlag = FALSE)
f0 <- raw; f0[, linear_vars] <- complete(imp, 20)
for (v in linear_vars) f0[[paste0(v, "Sqr")]] <- f0[[v]]^2
f0$log10z <- log10(f0$Redshift_crosscheck + 1)
pred_cols <- c(linear_vars, paste0(linear_vars, "Sqr"))
mm <- model.matrix(reformulate(pred_cols), data = f0); qrr <- qr(mm)
indep <- setdiff(colnames(mm)[qrr$pivot[seq_len(qrr$rank)]], "(Intercept)")
M_est <- MASS::rlm(reformulate(indep, response = "log10z"), data = f0, method = "M")
keep <- M_est$w > quantile(M_est$w, 0.05)
f0 <- f0[keep, ]
write.csv(f0, file.path(SC, "round9_frame.csv"))
cat("corrected frame:", nrow(f0), "rows\n")

# per-row measurement sigmas for noise injection (dex errors from raw where present)
sig <- matrix(0, nrow(f0), length(linear_vars), dimnames = list(NULL, linear_vars))
for (v in linear_vars) {
  ec <- err_map[[v]]
  if (!is.na(ec) && ec %in% names(raw)) {
    s <- raw[rownames(f0), ec]; s[!is.finite(s)] <- median(s[is.finite(s)], na.rm = TRUE)
    sig[, v] <- pmin(abs(s), 1)   # cap runaway errors at 1 dex
  }
}

Y <- f0$log10z; N <- nrow(f0)
is_opt <- f0$is_optical == 1
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")

N_REP <- 3L; K <- 10L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# ---- core library ----
fg <- read.table("Best_formula_GAM.txt"); gam_forms <- unique(apply(as.matrix(fg[, 2]), 1, as.formula))
smooth3 <- list(as.formula("Response ~ s(log10PeakFlux) + s(log10NH) + log10Ta + log10Fa + PhotonIndex + Alpha + Gamma + is_optical"))
band_form <- as.formula(paste(
  "Response ~ (log10Fa + log10Ta + log10NH + PhotonIndex + log10PeakFlux)^2 +",
  "Alpha + Beta + log10T90 + Gamma + log10Fluence +",
  "is_optical:(log10Fa + log10Ta + Alpha + Beta) + is_optical"))
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1]), detailed_names = FALSE, name_prefix = "gamF")
gamS <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = smooth3),      detailed_names = FALSE, name_prefix = "gamS")
glmB <- create.Learner("SL.custom_glm", tune = list(glm.model = list(band_form)), detailed_names = FALSE, name_prefix = "glmB")
SL.mgcv_sel <- function(Y, X, newX, family, obsWeights, ...) {
  sm <- intersect(linear_vars, colnames(X)); rest <- setdiff(colnames(X), sm)
  f <- as.formula(paste("Y ~", paste(c(sprintf("s(%s, k=5)", sm), rest), collapse = "+")))
  fit <- mgcv::gam(f, data = X, family = family, weights = obsWeights, select = TRUE, method = "REML")
  out <- list(object = fit); class(out) <- "SL.mgcv_sel"
  list(pred = as.numeric(mgcv::predict.gam(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.mgcv_sel <- function(object, newdata, ...)
  as.numeric(mgcv::predict.gam(object$object, newdata = newdata, type = "response"))
SL.glm_all <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- glm(Y ~ ., data = df, family = family, weights = obsWeights)
  out <- list(object = fit); class(out) <- "SL.glm_all"
  list(pred = as.numeric(predict(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.glm_all <- function(object, newdata, ...) as.numeric(predict(object$object, newdata = newdata, type = "response"))
SL.pls <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- pls::plsr(Y ~ ., data = df, ncomp = min(8, ncol(X)), validation = "CV")
  nc <- max(1, which.min(pls::RMSEP(fit)$val[1, 1, -1]))
  out <- list(object = fit, nc = nc); class(out) <- "SL.pls"
  list(pred = as.numeric(predict(fit, newdata = data.frame(newX, check.names = FALSE), ncomp = nc)), fit = out)
}
predict.SL.pls <- function(object, newdata, ...)
  as.numeric(predict(object$object, newdata = data.frame(newdata, check.names = FALSE), ncomp = object$nc))
CORE <- c(gamL$names, gamS$names, "SL.mgcv_sel", glmB$names, "SL.glm_all",
          "SL.glmnet", "SL.randomForest", "SL.mean")
FASTL <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.mean")

CFG <- list(
  anchor    = list(libs = CORE),
  daume     = list(libs = CORE, daume = TRUE),
  pls       = list(libs = c(CORE, "SL.pls")),
  noise     = list(libs = FASTL, noise_k = 5L),
  cmixup    = list(libs = FASTL, cmix_k = 2L),
  rankgauss = list(libs = CORE, rankg = TRUE)
)

daume_aug <- function(X) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

run_task <- function(task) {
  cfg <- CFG[[task$cfg]]
  folds <- fold_sets[[task$rep]]
  X <- f0[, VARS]
  if (isTRUE(cfg$daume)) X <- daume_aug(X)
  pred <- rep(NA_real_, N)
  set.seed(120000 + 100 * task$rep + match(task$cfg, names(CFG)))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- X[tr, , drop = FALSE]; ytr <- Y[tr]
    if (!is.null(cfg$noise_k)) {           # measurement-error noise injection
      reps <- lapply(seq_len(cfg$noise_k), function(k) {
        Xa <- Xtr
        Xa[, linear_vars] <- Xa[, linear_vars] + sig[tr, ] * matrix(rnorm(length(tr) * length(linear_vars)), length(tr))
        for (v in top7) Xa[[paste0(v, "Sqr")]] <- Xa[[v]]^2
        Xa
      })
      Xtr <- do.call(rbind, c(list(Xtr), reps))
      ytr <- rep(ytr, cfg$noise_k + 1)
    }
    if (!is.null(cfg$cmix_k)) {            # C-mixup augmentation
      bw <- sd(ytr) / 2
      Pm <- exp(-outer(ytr, ytr, "-")^2 / (2 * bw^2)); diag(Pm) <- 0
      idx <- seq_along(ytr)
      newX <- list(); newy <- c()
      for (k in seq_len(cfg$cmix_k)) {
        j <- vapply(idx, function(ii) sample(idx, 1, prob = Pm[ii, ]), integer(1))
        lam <- rbeta(length(idx), 2, 2)
        newX[[k]] <- Xtr * lam + Xtr[j, , drop = FALSE] * (1 - lam)
        newy <- c(newy, ytr * lam + ytr[j] * (1 - lam))
      }
      Xtr <- rbind(Xtr, do.call(rbind, newX)); ytr <- c(ytr, newy)
    }
    if (isTRUE(cfg$rankg)) {               # rank-gauss target (fold-internal map)
      yg <- qnorm((rank(ytr, ties.method = "average") - 0.5) / length(ytr))
      back <- function(p) approx(sort(yg), sort(ytr), xout = p, rule = 2)$y
      ymod <- yg
    } else { ymod <- ytr; back <- identity }
    fit <- try(suppressWarnings(SuperLearner(
      Y = ymod, X = Xtr, newX = X[te, , drop = FALSE], family = gaussian(),
      SL.library = cfg$libs, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", task$cfg, task$rep, i, "\n"); next }
    pred[te] <- back(as.numeric(fit$SL.predict))
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

cat("\n== round 9 (corrected frame, 10-fold, 3 reps) ==\n")
for (cn in names(CFG)) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$cfg == cn & res$rep == rp & is.finite(res$pred), ]
    c(all = cor(d$y, d$pred), op = cor(d$y[d$is_opt], d$pred[d$is_opt]))
  })
  cat(sprintf("%-10s r_all=%.4f  r_opt=%.4f (per-rep: %s)\n", cn,
              mean(rr["all", ]), mean(rr["op", ]), paste(sprintf("%.3f", rr["all", ]), collapse = " ")))
}
