# harness_prep.R — shared data prep + configurable single-stage CV runner.
# Source this from the REPO ROOT (GRB-Web-App/). Leaves in the global env:
#   frames (list of projection-draw frames), Y, Z, N, train_all, holdout_idx,
#   linear_vars, top7, VARS, f, and functions run_config()/report_cfg().
# A "config" is a list of knobs (defaults = current max-stack best). run_config
# runs reduced-setting paired CV (single stage, no cat-outlier removal) on the
# 216-row train set and returns per-rep r + pooled r. Identical folds/frames
# across configs => deltas are paired.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest); library(missForest)
  library(mboost); library(FNN)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

if (!exists("INPUT_CSV")) INPUT_CSV <- "Data/superlearner_training_emcee_v2_errcut_relative.csv"
# screening defaults (override before sourcing for full verification)
if (!exists("N_REP"))   N_REP   <- 2L
if (!exists("K"))       K       <- 5L
if (!exists("V_INNER")) V_INNER <- 5L
if (!exists("N_FRAMES"))N_FRAMES<- 4L

norm_grb <- function(s) {
  s <- trimws(gsub("GRB", "", as.character(s)))
  ifelse(grepl("[A-Za-z]$", s), s, paste0(s, "A"))
}
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")

f <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
names(f)[1] <- "GRB"; f$GRB <- norm_grb(f$GRB)
xr <- read.csv("Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv",
               stringsAsFactors = FALSE)
xray_ids <- norm_grb(xr[[1]])
f$is_optical <- as.integer(!(f$GRB %in% xray_ids))

f <- f[!is.na(f$log10T90) & f$log10T90 > 0.301, ]
f$log10NH[f$log10NH < 20] <- NA
f$Beta[f$Beta > 3] <- NA; f$Gamma[f$Gamma > 3] <- NA; f$Alpha[f$Alpha > 3] <- NA
f$PhotonIndex[f$PhotonIndex < 0] <- NA
cat("after clean:", nrow(f), "rows,", sum(is.na(f[, linear_vars])), "NAs\n")
if (sum(is.na(f[, linear_vars])) > 0) {
  set.seed(12); f[, linear_vars] <- missForest::missForest(f[, linear_vars])$ximp
}
f$log10z <- log10(f$Redshift_crosscheck + 1)
rownames(f) <- f$GRB

dfm <- data.frame(Y = f$log10z, f[, linear_vars])
rfit <- MASS::rlm(Y ~ ., data = dfm, method = "M", maxit = 100)
cut_n <- ceiling(0.05 * nrow(f))
drop_ids <- rownames(f)[order(rfit$w)[seq_len(cut_n)]]
f <- f[!(rownames(f) %in% drop_ids), ]
cat("M-est 5% cut ->", nrow(f), "rows (", sum(f$is_optical), "optical )\n")

set.seed(1)
cvl <- cv.glmnet(scale(as.matrix(f[, linear_vars])), f$log10z, alpha = 1)
co <- abs(as.numeric(coef(cvl, s = "lambda.min"))[-1])
top7 <- linear_vars[order(co, decreasing = TRUE)][1:7]
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")

Y <- f$log10z; Z <- f$Redshift_crosscheck; N <- nrow(f); is_opt <- f$is_optical == 1

# projection-draw frames
opt_cat <- read.csv("Data/OnlyLGRBs_data_171_optical_corrected.csv", stringsAsFactors = FALSE)
names(opt_cat)[1] <- "GRB"; opt_cat$GRB <- norm_grb(opt_cat$GRB)
mo <- match(f$GRB[is_opt], opt_cat$GRB)
onat <- data.frame(logFa = opt_cat$log10Faopt[mo], logTa = opt_cat$log10Taopt[mo],
                   Alpha = opt_cat$Alpha_opt[mo], Beta = opt_cat$Beta_opt[mo])
tgt_map <- c(logFa = "log10Fa", logTa = "log10Ta", Alpha = "Alpha", Beta = "Beta")
draws <- read.csv("proj_draws_v2.csv", stringsAsFactors = FALSE)
build_frames <- function(n_draws) {
  fr <- list()
  for (d in seq_len(n_draws)) {
    fd <- f; ri <- which(is_opt)
    for (p in names(tgt_map)) {
      dr <- draws[draws$param == p & draws$draw == d - 1, ]
      val <- dr$m * onat[[p]] + dr$b; tv <- tgt_map[[p]]; ok <- is.finite(val)
      fd[[tv]][ri[ok]] <- val[ok]
    }
    for (v in linear_vars) fd[[paste0(v, "Sqr")]] <- fd[[v]]^2
    fr[[d]] <- fd
  }
  fr
}
frames <- build_frames(N_FRAMES)
cat("built", length(frames), "frames; N_REP", N_REP, "K", K, "V_INNER", V_INNER, "\n")

# ---- shared learner pieces ----
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
  cvr <- try(mboost::cvrisk(fit, folds = mboost::cv(model.weights(fit), type = "kfold", B = 5), papply = lapply), silent = TRUE)
  if (!inherits(cvr, "try-error")) fit <- fit[mboost::mstop(cvr)]
  out <- list(object = fit); class(out) <- "SL.gamboost"
  list(pred = as.numeric(predict(fit, newdata = data.frame(newX, check.names = FALSE))), fit = out)
}
predict.SL.gamboost <- function(object, newdata, ...) as.numeric(predict(object$object, newdata = data.frame(newdata, check.names = FALSE)))
method.NNRidge <- function() {
  out <- SuperLearner::method.NNLS()
  out$computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights, ...) {
    fit <- glmnet::cv.glmnet(Z, Y, alpha = 0, lower.limits = 0, nfolds = 10)
    co <- as.numeric(coef(fit, s = "lambda.min"))[-1]; co[co < 0] <- 0
    if (sum(co) > 0) co <- co / sum(co) else co <- rep(1/ncol(Z), ncol(Z))
    list(cvRisk = apply(Z, 2, function(p) mean((p - Y)^2)), coef = co, optimizer = fit)
  }
  out
}
band_form_default <- as.formula(paste(
  "Response ~ (log10Fa + log10Ta + log10NH + PhotonIndex + log10PeakFlux)^2 +",
  "Alpha + Beta + log10T90 + Gamma + log10Fluence +",
  "is_optical:(log10Fa + log10Ta + Alpha + Beta) + is_optical"))
fg <- read.table("Best_formula_GAM.txt"); gam_forms_default <- unique(apply(as.matrix(fg[, 2]), 1, as.formula))

daume_aug <- function(X, c_scale = 0.25) {
  if (c_scale <= 0) return(X)
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical * c_scale
  colnames(Xo) <- paste0(colnames(Xo), "_opt"); cbind(X, Xo)
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
  r <- residuals(rf); keep <- abs(r) <= 3 * mad(r)
  if (sum(keep) >= 30) tr_idx[keep] else tr_idx
}

# pipeline convention: seed-42 20% holdout carved out BEFORE CV
set.seed(42)
holdout_idx <- sample(N, size = floor(0.20 * N))
train_all   <- setdiff(seq_len(N), holdout_idx)
cat("split: train", length(train_all), "/ holdout", length(holdout_idx), "\n")

# ---- default config ----
default_cfg <- function() list(
  c_scale = 0.25, gam_forms = gam_forms_default, band_form = band_form_default,
  use_band_glm = TRUE, use_glm_all = TRUE, use_gamboost = TRUE, use_rf = TRUE,
  use_glmnet = TRUE, use_iso = TRUE, use_infold = TRUE, meta = "NNRidge",
  extra_vars = character(0), extra_learners = character(0), vars = VARS
)

build_library <- function(cfg) {
  lib <- character(0)
  if (is.null(cfg$use_gam) || cfg$use_gam) {
    gamL <- create.Learner("SL.mgcv_gam", tune = list(gam.model = cfg$gam_forms[1]),
                            detailed_names = FALSE, name_prefix = "gamF", env = globalenv())
    lib <- c(lib, gamL$names)
  }
  if (cfg$use_band_glm) {
    glmB <- create.Learner("SL.custom_glm", tune = list(glm.model = list(cfg$band_form)),
                           detailed_names = FALSE, name_prefix = "glmB", env = globalenv())
    lib <- c(lib, glmB$names)
  }
  if (cfg$use_glm_all) lib <- c(lib, "SL.glm_all")
  if (cfg$use_glmnet)  lib <- c(lib, "SL.glmnet")
  if (cfg$use_rf)      lib <- c(lib, "SL.randomForest")
  if (cfg$use_gamboost)lib <- c(lib, "SL.gamboost")
  lib <- c(lib, cfg$extra_learners, "SL.mean")
  lib
}

# single-stage paired CV for one config on train_idx; returns list(r_pool, r_by_rep, preds_by_row)
run_config <- function(cfg, train_idx, frame_list, label = "cfg") {
  meth <- if (cfg$meta == "NNLS") SuperLearner::method.NNLS() else method.NNRidge()
  set.seed(31415)
  fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y[train_idx], k = K, returnTrain = FALSE))
  grid <- expand.grid(rep = seq_len(N_REP), fr = seq_along(frame_list))
  tasks <- split(grid, seq_len(nrow(grid)))
  run_task <- function(t) {
    frame <- frame_list[[t$fr]]
    lib <- build_library(cfg)
    Xall <- daume_aug(frame[, cfg$vars], cfg$c_scale)
    folds <- fold_sets[[t$rep]]
    pred_cv <- rep(NA_real_, length(train_idx))
    set.seed(270000 + 1000 * t$rep + t$fr)
    for (i in seq_along(folds)) {
      te_loc <- folds[[i]]; tr_loc <- setdiff(seq_along(train_idx), te_loc)
      tr <- train_idx[tr_loc]; te <- train_idx[te_loc]
      if (cfg$use_infold) tr <- infold_filter(tr, frame)
      Xtr <- Xall[tr, ]; Xte <- Xall[te, ]
      if (cfg$use_iso) { ai <- add_iso(Xtr, Xte); Xtr <- ai$tr; Xte <- ai$te }
      fit <- try(suppressWarnings(SuperLearner(
        Y = Y[tr], X = Xtr, newX = Xte, family = gaussian(),
        SL.library = lib, method = meth, cvControl = list(V = V_INNER), verbose = FALSE)), silent = TRUE)
      if (inherits(fit, "try-error")) next
      pred_cv[te_loc] <- as.numeric(fit$SL.predict)
    }
    data.frame(rep = t$rep, fr = t$fr, row = train_idx, pred = pred_cv)
  }
  res <- parallel::mclapply(tasks, run_task, mc.cores = 14L, mc.preschedule = FALSE)
  res <- do.call(rbind, res)
  # pooled r: average over frames+reps per row
  pa <- tapply(res$pred, res$row, mean, na.rm = TRUE)
  rows <- as.integer(names(pa)); r_pool <- cor(Y[rows], as.numeric(pa))
  # per-rep r: average over frames within rep
  r_by_rep <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$rep == rp, ]; pr <- tapply(d$pred, d$row, mean, na.rm = TRUE)
    rr <- as.integer(names(pr)); cor(Y[rr], as.numeric(pr))
  })
  cat(sprintf("  [%s] r_pool=%.4f  r_by_rep=%s\n", label, r_pool, paste(sprintf("%.4f", r_by_rep), collapse=",")))
  list(label = label, r_pool = r_pool, r_by_rep = r_by_rep, preds = pa, rows = rows)
}

run_set <- function(cfg_list, train_idx = train_all, frame_list = frames) {
  out <- list()
  for (nm in names(cfg_list)) out[[nm]] <- run_config(cfg_list[[nm]], train_idx, frame_list, nm)
  # paired deltas vs first config (baseline)
  base <- out[[1]]
  cat("\n=== paired deltas vs", names(cfg_list)[1], "===\n")
  for (nm in names(cfg_list)) {
    d_pool <- out[[nm]]$r_pool - base$r_pool
    d_rep <- out[[nm]]$r_by_rep - base$r_by_rep
    cat(sprintf("%-20s r_pool=%.4f  d_pool=%+.4f  d_by_rep=%s  (mean %+.4f, %d/%d +)\n",
                nm, out[[nm]]$r_pool, d_pool, paste(sprintf("%+.3f", d_rep), collapse=","),
                mean(d_rep), sum(d_rep > 0), length(d_rep)))
  }
  out
}
cat("harness_prep loaded. train_all n=", length(train_all), "\n")
