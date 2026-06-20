#!/usr/bin/env Rscript
#' =============================================================================
#' translation_cv.R -- CV comparison of translation-encoder approaches
#' =============================================================================
#'
#' Four arms per fold, identical GRB-level folds (assigned over the union):
#'
#'   A  : x-ray baseline. SuperLearner (lean library: tuned GAM/GLM/bayesglm
#'        formulas + rf + glmnet + xgboost + ranger + mean) on x-ray GRBs,
#'        10 features + squares. The bar to beat on x-ray GRBs.
#'   O  : optical baseline. SuperLearner (generic library) on optical GRBs,
#'        5 optical features + squares. The bar to beat on optical-only GRBs.
#'   T1 : translation as augmentation. Ridge translation (optical -> x-ray
#'        feature space, LOO-tuned on training pairs) creates pseudo x-ray rows
#'        for optical-only training GRBs (weight TRANS_WEIGHT); arm-A library
#'        trains on real + pseudo rows; optical-only test GRBs are predicted
#'        through their translated features.
#'   T2 : end-to-end. Translation matrix + quadratic ridge head trained jointly
#'        (prediction loss backprops through the head into the translation;
#'        feature-match loss on pairs; warm-started from the ridge translation).
#'
#' Leakage control: MICE, feature scaling, translation fitting and SL training
#' all use only the fold's training GRBs (per-fold mice_complete_fold with
#' ignore = test, unlike concat_superlearner.R's full-data MICE).
#'
#' Usage:
#'   Rscript translation_cv.R [reps] [k_folds]            # defaults 3, 10
#'   SMOKE_TEST=true Rscript translation_cv.R 1           # structural check
#'   Env: TRANS_CORES (3), TRANS_WEIGHT (1.0), MATCH_WEIGHT (1.0)
#'
#' Outputs: OutputFiles/TranslationExperiment/{cv_predictions.csv, metrics.csv,
#'          fold_diagnostics.csv, progress.log}
#' =============================================================================

source("translation_model.R")
suppressPackageStartupMessages({
  library(SuperLearner)
  library(caret)
  library(parallel)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")
source("Custom_SL/sl_custom_bayesglm.R")
source("Custom_SL/sl_xgboost_safe.R")

args         <- commandArgs(trailingOnly = TRUE)
reps         <- if (length(args) >= 1) as.integer(args[1]) else 3L
k_folds      <- if (length(args) >= 2) as.integer(args[2]) else 10L
smoke        <- tolower(Sys.getenv("SMOKE_TEST", "false")) == "true"
n_cores      <- max(1L, as.integer(Sys.getenv("TRANS_CORES", "3")))
trans_weight <- as.numeric(Sys.getenv("TRANS_WEIGHT", "1.0"))
match_weight <- as.numeric(Sys.getenv("MATCH_WEIGHT", "1.0"))
out_dir      <- "OutputFiles/TranslationExperiment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data ---------------------------------------------------------------------
matched <- match_modalities(
  load_xray_modality("x-ray_data.csv"),
  load_optical_modality("optical_data.txt")
)
xray    <- matched$xray
optical <- matched$optical

sample_ids <- union(xray$sample_id, optical$sample_id)
xi <- match(sample_ids, xray$sample_id)
oi <- match(sample_ids, optical$sample_id)
units <- data.frame(
  sample_id    = sample_ids,
  GRB          = ifelse(!is.na(xi), xray$id[xi], optical$id[oi]),
  availability = ifelse(!is.na(xi) & !is.na(oi), "both",
                 ifelse(!is.na(xi), "xray", "optical")),
  y            = rowMeans(cbind(xray$y[xi], optical$y[oi]), na.rm = TRUE),
  stringsAsFactors = FALSE
)
stopifnot(!anyNA(units$y))
cat(sprintf("Units: %d (%d xray-only, %d both, %d optical-only)\n",
            nrow(units), sum(units$availability == "xray"),
            sum(units$availability == "both"), sum(units$availability == "optical")))

# ---- Learner libraries ----------------------------------------------------------
formula_table_GAM <- read.table("Best_formula_GAM.txt")
bestGAM1 <- apply(as.matrix(formula_table_GAM[, 2]), 1, as.formula)
learner1 <- create.Learner("SL.mgcv_gam", tune = list(gam.model = c(bestGAM1)),
                           detailed_names = FALSE, name_prefix = "gam")
formula_table_GLM <- read.table("Best_formula_GLM.txt")
best_lm3 <- apply(as.matrix(formula_table_GLM[, 2]), 1, as.formula)
sl_glm1  <- create.Learner("SL.custom_glm", tune = list(glm.model = c(best_lm3)),
                           detailed_names = FALSE, name_prefix = "cglm")
sl_bglm  <- create.Learner("SL.custom_bayesglm", tune = list(bglm.model = c(best_lm3)),
                           detailed_names = FALSE, name_prefix = "bglm")
caret_learner <- create.Learner("SL.caret",
                                tune = list(method = "rf", tuneLength = 1, verboseIter = FALSE),
                                detailed_names = TRUE)
libs_xray <- c(learner1$names, sl_glm1$names, sl_bglm$names, caret_learner$names,
               "SL.glmnet", "SL.xgboost_safe", "SL.ranger", "SL.mean")
# Optical features have different names, so the tuned formulas do not apply.
libs_optical <- c("SL.custom_glm", "SL.glmnet", "SL.xgboost_safe", "SL.ranger", "SL.mean")
if (smoke) {
  libs_xray <- libs_optical <- c("SL.mean", "SL.lm")
  k_folds <- 5L
  cat("SMOKE TEST: libraries reduced to SL.mean, SL.lm\n")
}

fit_sl <- function(Y, X, newX, lib, seed, weights = NULL) {
  set.seed(seed)
  suppressMessages(capture.output(
    fit <- SuperLearner(Y = Y, X = X, newX = newX, family = gaussian(),
                        SL.library = lib, cvControl = list(V = 5),
                        obsWeights = weights, verbose = FALSE),
    file = nullfile()
  ))
  drop(fit$SL.predict)
}

append_squares <- function(d) {
  d <- as.data.frame(d)
  for (nm in colnames(d)) d[[paste0(nm, "Sqr")]] <- d[[nm]]^2
  d
}

# ---- One fold -------------------------------------------------------------------
run_fold <- function(f, fold_of_unit, rep) {
  is_test  <- fold_of_unit == f
  test_ids <- units$sample_id[is_test]

  x_train <- !(xray$sample_id %in% test_ids)     # logical over x-ray modality rows
  o_train <- !(optical$sample_id %in% test_ids)  # logical over optical modality rows

  # Per-fold MICE (train rows only inform the imputation model) and scaling.
  x_complete <- mice_complete_fold(xray$imputation, x_train, m = 5L, maxit = 5L,
                                   seed = 100L * rep + f)
  o_complete <- mice_complete_fold(optical$imputation, o_train, m = 5L, maxit = 5L,
                                   seed = 100L * rep + f + 50L)
  x_scaled <- scale_fold_features(x_complete, x_train, xray$model_features)
  o_scaled <- scale_fold_features(o_complete, o_train, optical$model_features)
  stopifnot(identical(x_scaled$features, xray$model_features))

  # Optical inputs for translation: scaled features + squares.
  O_input <- cbind(o_scaled$matrix, o_scaled$matrix^2)
  colnames(O_input) <- c(colnames(o_scaled$matrix), paste0(colnames(o_scaled$matrix), "Sqr"))

  # Unscaled completed x-ray features (for the SuperLearner arms).
  X_raw <- as.data.frame(x_complete[, xray$model_features, drop = FALSE])
  unscale <- function(M) sweep(sweep(M, 2L, x_scaled$scale, "*"), 2L, x_scaled$center, "+")

  # Index helpers (unit -> modality row).
  u_xrow <- match(units$sample_id, xray$sample_id)
  u_orow <- match(units$sample_id, optical$sample_id)
  has_x  <- !is.na(u_xrow)
  has_o  <- !is.na(u_orow)

  train_units <- which(!is_test)
  test_units  <- which(is_test)

  # Training pairs (both modalities, training fold) for the translation fit.
  pair_units <- which(units$availability == "both" & !is_test)
  O_pair <- O_input[u_orow[pair_units], , drop = FALSE]
  X_pair <- x_scaled$matrix[u_xrow[pair_units], , drop = FALSE]
  translation <- fit_ridge_translation(O_pair, X_pair)

  # Translation quality on TEST pairs (out-of-sample feature-space fidelity).
  test_pairs <- which(units$availability == "both" & is_test)
  test_translation_mse <- if (length(test_pairs)) {
    mean((predict_translation(translation, O_input[u_orow[test_pairs], , drop = FALSE]) -
          x_scaled$matrix[u_xrow[test_pairs], , drop = FALSE])^2)
  } else NA_real_

  # ---- Arm A: x-ray baseline -------------------------------------------------
  a_train <- intersect(train_units, which(has_x))
  a_test  <- intersect(test_units,  which(has_x))
  XA_train <- append_squares(X_raw[u_xrow[a_train], , drop = FALSE])
  XA_test  <- append_squares(X_raw[u_xrow[a_test],  , drop = FALSE])
  pred_A <- setNames(fit_sl(units$y[a_train], XA_train, XA_test, libs_xray,
                            seed = rep * 1000L + f), units$sample_id[a_test])

  # ---- Arm O: optical baseline -----------------------------------------------
  o_train_units <- intersect(train_units, which(has_o))
  o_test_units  <- intersect(test_units,  which(has_o))
  O_raw <- as.data.frame(o_complete[, o_scaled$features, drop = FALSE])
  pred_O <- setNames(fit_sl(units$y[o_train_units],
                            append_squares(O_raw[u_orow[o_train_units], , drop = FALSE]),
                            append_squares(O_raw[u_orow[o_test_units],  , drop = FALSE]),
                            libs_optical, seed = rep * 1000L + f + 1L),
                     units$sample_id[o_test_units])

  # ---- Arm T1: translation as augmentation ------------------------------------
  oo_train <- intersect(train_units, which(units$availability == "optical"))
  oo_test  <- intersect(test_units,  which(units$availability == "optical"))
  pseudo_train <- as.data.frame(unscale(
    predict_translation(translation, O_input[u_orow[oo_train], , drop = FALSE])))
  colnames(pseudo_train) <- xray$model_features
  pseudo_test <- as.data.frame(unscale(
    predict_translation(translation, O_input[u_orow[oo_test], , drop = FALSE])))
  colnames(pseudo_test) <- xray$model_features

  T1_train <- rbind(XA_train, append_squares(pseudo_train))
  T1_y     <- c(units$y[a_train], units$y[oo_train])
  T1_w     <- c(rep(1, length(a_train)), rep(trans_weight, length(oo_train)))
  T1_test  <- rbind(XA_test, append_squares(pseudo_test))
  pred_T1 <- setNames(fit_sl(T1_y, T1_train, T1_test, libs_xray,
                             seed = rep * 1000L + f + 2L, weights = T1_w),
                      units$sample_id[c(a_test, oo_test)])

  # ---- Arm T2: end-to-end ------------------------------------------------------
  y_center <- mean(units$y[train_units]); y_scale <- sd(units$y[train_units])
  e2e <- fit_end_to_end_translation(
    X_real = x_scaled$matrix[u_xrow[a_train], , drop = FALSE],
    y_real = (units$y[a_train] - y_center) / y_scale,
    O_only = O_input[u_orow[oo_train], , drop = FALSE],
    y_only = (units$y[oo_train] - y_center) / y_scale,
    O_pair = O_pair, X_pair = X_pair,
    match_weight = match_weight, seed = rep * 1000L + f + 3L,
    init_translation = translation
  )
  t2_x <- predict(e2e, x_scaled$matrix[u_xrow[a_test], , drop = FALSE], "xray")$prediction
  t2_o <- predict(e2e, O_input[u_orow[oo_test], , drop = FALSE], "optical")$prediction
  pred_T2 <- setNames(c(t2_x, t2_o) * y_scale + y_center,
                      units$sample_id[c(a_test, oo_test)])

  msg <- sprintf("[%s] rep %d fold %d done: %d test (%d xray, %d optical-only); ridge lambda=%.3g LOO=%.3f testMSE=%.3f\n",
                 format(Sys.time(), "%H:%M:%S"), rep, f, length(test_units),
                 length(a_test), length(oo_test),
                 translation$lambda, translation$loo_mse, test_translation_mse)
  cat(msg)
  cat(msg, file = file.path(out_dir, "progress.log"), append = TRUE)

  ids <- units$sample_id[test_units]
  list(
    predictions = data.frame(
      rep = rep, fold = f, sample_id = ids,
      pred_A  = unname(pred_A[ids]),
      pred_O  = unname(pred_O[ids]),
      pred_T1 = unname(pred_T1[ids]),
      pred_T2 = unname(pred_T2[ids]),
      stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      rep = rep, fold = f, n_pairs_train = length(pair_units),
      ridge_lambda = translation$lambda, ridge_loo_mse = translation$loo_mse,
      test_translation_mse = test_translation_mse,
      e2e_loss = e2e$loss, e2e_pred_loss = e2e$pred_loss,
      e2e_match_loss = e2e$match_loss
    )
  )
}

# ---- CV loop -----------------------------------------------------------------------
all_preds <- list(); all_diag <- list()
for (rep in seq_len(reps)) {
  set.seed(2026L + rep)
  fold_list <- createFolds(units$y, k = k_folds)
  fold_of_unit <- integer(nrow(units))
  for (f in seq_along(fold_list)) fold_of_unit[fold_list[[f]]] <- f
  cat(sprintf("=== Repetition %d/%d (%d folds, %d cores) ===\n", rep, reps, k_folds, n_cores))
  RNGkind("L'Ecuyer-CMRG"); set.seed(rep)
  results <- mclapply(seq_len(k_folds), run_fold,
                      fold_of_unit = fold_of_unit, rep = rep, mc.cores = n_cores)
  bad <- !vapply(results, is.list, logical(1))
  if (any(bad)) stop("Fold failure in rep ", rep, ":\n",
                     paste(utils::capture.output(print(results[[which(bad)[1]]])), collapse = "\n"))
  all_preds[[rep]] <- do.call(rbind, lapply(results, `[[`, "predictions"))
  all_diag[[rep]]  <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  partial <- merge(do.call(rbind, all_preds), units, by = "sample_id", sort = FALSE)
  write.csv(partial, file.path(out_dir, "cv_predictions_partial.csv"), row.names = FALSE)
}

predictions <- merge(do.call(rbind, all_preds), units, by = "sample_id", sort = FALSE)
diagnostics <- do.call(rbind, all_diag)
write.csv(predictions, file.path(out_dir, "cv_predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "fold_diagnostics.csv"), row.names = FALSE)

# ---- Metrics --------------------------------------------------------------------------
scopes <- list(
  xray_grbs    = predictions$availability %in% c("xray", "both"),
  both_subset  = predictions$availability == "both",
  xronly_subset = predictions$availability == "xray",
  optical_only = predictions$availability == "optical",
  all_grbs     = rep(TRUE, nrow(predictions))
)
arms <- c(A = "pred_A", O = "pred_O", T1 = "pred_T1", T2 = "pred_T2")
metrics <- do.call(rbind, lapply(names(arms), function(arm) {
  do.call(rbind, lapply(names(scopes), function(scope) {
    d <- predictions[scopes[[scope]], ]
    if (all(!is.finite(d[[arms[[arm]]]]))) return(NULL)
    cbind(arm = arm, scope = scope, regression_metrics(d$y, d[[arms[[arm]]]]))
  }))
}))
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
cat("\n==== RESULTS (y = log10(z+1); _z = linear) ====\n")
print(metrics, row.names = FALSE, digits = 3)
cat("\nTranslation diagnostics (mean over folds): LOO MSE =",
    round(mean(diagnostics$ridge_loo_mse), 3),
    "| test MSE =", round(mean(diagnostics$test_translation_mse, na.rm = TRUE), 3),
    "(scaled feature units; 1.0 = predicting the mean)\n")
cat("Wrote", file.path(out_dir, "cv_predictions.csv"), ", metrics.csv, fold_diagnostics.csv\n")
