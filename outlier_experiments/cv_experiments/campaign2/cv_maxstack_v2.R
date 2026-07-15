# cv_maxstack_v2.R — max-stack (round-25 config) ported to the new v2/v3 frames.
# Usage: Rscript cv_maxstack_v2.R <input_frame.csv> <out_prefix>
# Pipeline convention: clean -> impute leftovers -> fresh 5% M-est cut ->
# seed-42 20% holdout -> CV on train (3 reps x 10 folds, bagged over 8
# projection-draw frames) -> retrained 2-sigma catastrophic-outlier pass.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(missForest)
  library(mboost); library(FNN)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

a <- commandArgs(trailingOnly = TRUE)
input_csv <- a[1]; out_prefix <- a[2]

norm_grb <- function(s) {
  s <- trimws(gsub("GRB", "", as.character(s)))
  ifelse(grepl("[A-Za-z]$", s), s, paste0(s, "A"))
}

linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")

f <- read.csv(input_csv, stringsAsFactors = FALSE)
names(f)[1] <- "GRB"; f$GRB <- norm_grb(f$GRB)
xr <- read.csv("Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv",
               stringsAsFactors = FALSE)
xray_ids <- norm_grb(xr[[1]])
f$is_optical <- as.integer(!(f$GRB %in% xray_ids))

# clean (pipeline rules)
f <- f[!is.na(f$log10T90) & f$log10T90 > 0.301, ]
f$log10NH[f$log10NH < 20] <- NA
f$Beta[f$Beta > 3] <- NA; f$Gamma[f$Gamma > 3] <- NA; f$Alpha[f$Alpha > 3] <- NA
f$PhotonIndex[f$PhotonIndex < 0] <- NA
cat("after clean:", nrow(f), "rows,", sum(is.na(f[, linear_vars])), "NAs in features\n")
if (sum(is.na(f[, linear_vars])) > 0) {
  set.seed(12)
  f[, linear_vars] <- missForest::missForest(f[, linear_vars])$ximp
}
f$log10z <- log10(f$Redshift_crosscheck + 1)
rownames(f) <- f$GRB

# fresh 5% M-estimator cut (Huber weights, drop lowest 5%)
dfm <- data.frame(Y = f$log10z, f[, linear_vars])
rfit <- MASS::rlm(Y ~ ., data = dfm, method = "M", maxit = 100)
w <- rfit$w
cut_n <- ceiling(0.05 * nrow(f))
drop_ids <- rownames(f)[order(w)[seq_len(cut_n)]]
f <- f[!(rownames(f) %in% drop_ids), ]
cat("M-est 5% cut removed", cut_n, "->", nrow(f), "rows (", sum(f$is_optical), "optical )\n")

# LASSO top-7 for squared terms
set.seed(1)
cvl <- cv.glmnet(scale(as.matrix(f[, linear_vars])), f$log10z, alpha = 1)
co <- abs(as.numeric(coef(cvl, s = "lambda.min"))[-1])
top7 <- linear_vars[order(co, decreasing = TRUE)][1:7]
cat("top7:", paste(top7, collapse = ", "), "\n")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")

Y <- f$log10z; Z <- f$Redshift_crosscheck; N <- nrow(f); is_opt <- f$is_optical == 1

# projection-draw frames: redraw the 4 emcee-projected features for optical rows
opt_cat <- read.csv("Data/OnlyLGRBs_data_171_optical_corrected.csv", stringsAsFactors = FALSE)
names(opt_cat)[1] <- "GRB"; opt_cat$GRB <- norm_grb(opt_cat$GRB)
mo <- match(f$GRB[is_opt], opt_cat$GRB)
onat <- data.frame(logFa = opt_cat$log10Faopt[mo], logTa = opt_cat$log10Taopt[mo],
                   Alpha = opt_cat$Alpha_opt[mo], Beta = opt_cat$Beta_opt[mo])
tgt_map <- c(logFa = "log10Fa", logTa = "log10Ta", Alpha = "Alpha", Beta = "Beta")
draws <- read.csv("proj_draws_v2.csv", stringsAsFactors = FALSE)
N_DRAWS <- 8L
frames <- list()
for (d in seq_len(N_DRAWS)) {
  fd <- f
  ri <- which(is_opt)
  for (p in names(tgt_map)) {
    dr <- draws[draws$param == p & draws$draw == d - 1, ]
    val <- dr$m * onat[[p]] + dr$b
    tv <- tgt_map[[p]]
    ok <- is.finite(val)
    fd[[tv]][ri[ok]] <- val[ok]
  }
  for (v in linear_vars) fd[[paste0(v, "Sqr")]] <- fd[[v]]^2
  frames[[d]] <- fd
}
cat("built", length(frames), "projection-draw frames\n")

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
SL.gamboost <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  num <- intersect(linear_vars, colnames(X))
  fm <- as.formula(paste("Y ~", paste(c(sprintf("bbs(%s, df = 2)", num), "bols(is_optical)"), collapse = "+")))
  fit <- mboost::gamboost(fm, data = df, control = mboost::boost_control(mstop = 300, nu = 0.1))
  cvr <- try(mboost::cvrisk(fit, folds = mboost::cv(model.weights(fit), type = "kfold", B = 5),
                            papply = lapply), silent = TRUE)
  if (!inherits(cvr, "try-error")) fit <- fit[mboost::mstop(cvr)]
  out <- list(object = fit); class(out) <- "SL.gamboost"
  list(pred = as.numeric(predict(fit, newdata = data.frame(newX, check.names = FALSE))), fit = out)
}
predict.SL.gamboost <- function(object, newdata, ...)
  as.numeric(predict(object$object, newdata = data.frame(newdata, check.names = FALSE)))
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
LEAN <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.randomForest", "SL.gamboost", "SL.mean")

daume_aug <- function(X, c_scale = 0.25) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical * c_scale
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}
add_iso <- function(Xtr, Xte) {
  sc <- scale(Xtr[, linear_vars])
  Xtr$outlier_score <- FNN::knn.dist(sc, k = 10)[, 10]
  Xte$outlier_score <- FNN::get.knnx(sc, scale(Xte[, linear_vars],
    center = attr(sc, "scaled:center"), scale = attr(sc, "scaled:scale")), k = 10)$nn.dist[, 10]
  list(tr = Xtr, te = Xte)
}
infold_filter <- function(tr_idx, frame) {
  dfr <- data.frame(Y = Y[tr_idx], frame[tr_idx, c(linear_vars, paste0(top7, "Sqr"))])
  rf <- tryCatch(MASS::rlm(Y ~ ., data = dfr, method = "M", maxit = 50), error = function(e) NULL)
  if (is.null(rf)) return(tr_idx)
  r <- residuals(rf)
  keep <- abs(r) <= 3 * mad(r)
  if (sum(keep) >= 30) tr_idx[keep] else tr_idx
}

# pipeline convention: seed-42 20% holdout carved out BEFORE CV
set.seed(42)
holdout_idx <- sample(N, size = floor(0.20 * N))
train_all   <- setdiff(seq_len(N), holdout_idx)
cat("pipeline split: train", length(train_all), "/ holdout", length(holdout_idx), "\n")

N_REP <- 3L; K <- 10L

run_stage <- function(train_idx, tag) {
  set.seed(31415)
  fold_sets <- lapply(seq_len(N_REP), function(r)
    caret::createFolds(Y[train_idx], k = K, returnTrain = FALSE))
  grid <- expand.grid(rep = seq_len(N_REP), fr = seq_along(frames))
  tasks <- split(grid, seq_len(nrow(grid)))
  cat(tag, ": running", length(tasks), "tasks, train n =", length(train_idx), "\n")
  run_task <- function(t) {
    frame <- frames[[t$fr]]
    Xall <- daume_aug(frame[, VARS])
    folds <- fold_sets[[t$rep]]
    pred_cv <- rep(NA_real_, length(train_idx))
    set.seed(270000 + 1000 * t$rep + t$fr)
    for (i in seq_along(folds)) {
      te_loc <- folds[[i]]; tr_loc <- setdiff(seq_along(train_idx), te_loc)
      tr <- train_idx[tr_loc]; te <- train_idx[te_loc]
      tr <- infold_filter(tr, frame)
      ai <- add_iso(Xall[tr, ], Xall[te, ])
      fit <- try(suppressWarnings(SuperLearner(
        Y = Y[tr], X = ai$tr, newX = ai$te, family = gaussian(),
        SL.library = LEAN, method = method.NNRidge(),
        cvControl = list(V = 10), verbose = FALSE)), silent = TRUE)
      if (inherits(fit, "try-error")) { cat("ERR", tag, "rep", t$rep, "frame", t$fr, "fold", i, "\n"); next }
      pred_cv[te_loc] <- as.numeric(fit$SL.predict)
    }
    trf <- infold_filter(train_idx, frame)
    ai <- add_iso(Xall[trf, ], Xall[holdout_idx, ])
    fit_full <- try(suppressWarnings(SuperLearner(
      Y = Y[trf], X = ai$tr, newX = ai$te, family = gaussian(),
      SL.library = LEAN, method = method.NNRidge(),
      cvControl = list(V = 10), verbose = FALSE)), silent = TRUE)
    ph <- if (inherits(fit_full, "try-error")) rep(NA_real_, length(holdout_idx)) else as.numeric(fit_full$SL.predict)
    rbind(
      data.frame(rep = t$rep, fr = t$fr, set = "cv",      row = train_idx,   pred = pred_cv),
      data.frame(rep = t$rep, fr = t$fr, set = "holdout", row = holdout_idx, pred = ph)
    )
  }
  res <- parallel::mclapply(tasks, run_task, mc.cores = 14L, mc.preschedule = FALSE)
  do.call(rbind, res)
}

report <- function(res, tag) {
  out <- list()
  for (s in c("cv", "holdout")) {
    d <- res[res$set == s, ]
    pa <- tapply(d$pred, d$row, mean, na.rm = TRUE)
    rows <- as.integer(names(pa))
    yy <- Y[rows]; zz <- Z[rows]
    zp <- 10^pa - 1
    cat(sprintf("%s %-8s n=%d  log_r=%.4f  lin_r=%.4f\n",
                tag, s, length(pa), cor(yy, pa), cor(zz, zp)))
    out[[s]] <- data.frame(row = rows, GRB = f$GRB[rows], y = yy, z = zz, pred = as.numeric(pa), z_pred = zp)
  }
  out
}

res1 <- run_stage(train_all, "stage1")
write.csv(res1, paste0(out_prefix, "_stage1_preds.csv"), row.names = FALSE)
agg1 <- report(res1, "STAGE1 (pre cat-outlier removal)")

# retrained catastrophic-outlier pass: flag train rows with pooled CV resid > 2 sigma
cvp <- agg1$cv
resid <- cvp$y - cvp$pred
sdr <- sd(resid)
catoutl <- cvp$row[abs(resid) > 2 * sdr]
cat("catastrophic outliers (>2sd):", length(catoutl), ":", paste(f$GRB[catoutl], collapse = ", "), "\n")
train2 <- setdiff(train_all, catoutl)

res2 <- run_stage(train2, "stage2")
write.csv(res2, paste0(out_prefix, "_stage2_preds.csv"), row.names = FALSE)
agg2 <- report(res2, "FINAL (cat-outliers removed, retrained)")
write.csv(agg2$cv, paste0(out_prefix, "_final_cv.csv"), row.names = FALSE)
write.csv(agg2$holdout, paste0(out_prefix, "_final_holdout.csv"), row.names = FALSE)
writeLines(f$GRB[catoutl], paste0(out_prefix, "_catoutl.txt"))
cat("DONE_MAXSTACK\n")
