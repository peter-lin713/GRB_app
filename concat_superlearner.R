#!/usr/bin/env Rscript
#' =============================================================================
#' concat_superlearner.R -- Does appending optical features help the original
#' superlearner.R training scheme, vs x-ray features alone?
#' =============================================================================
#'
#' Two arms, identical training scheme (same learner library as superlearner.R,
#' same k-fold CV with inner V=5 ensemble weights, MICE midastouch m=20):
#'
#'   Arm A (baseline) : x-ray GRBs only, the 10 x-ray features (+ squares).
#'                      Replicates the original x-ray-only pipeline.
#'   Arm B (concat)   : union of x-ray and optical GRBs, 10 x-ray + 5 optical
#'                      features (+ squares) in ONE table; missing modality
#'                      columns are NA and are imputed by a single joint MICE
#'                      run, so the ~80 paired GRBs teach MICE the cross-modal
#'                      relationships used to fill the unpaired rows.
#'
#' Folds are assigned once per repetition at the GRB level over the UNION, and
#' both arms reuse them, so every x-ray GRB is out-of-fold in the same fold in
#' both arms. The primary comparison is r / RMSE of out-of-fold predictions on
#' the x-ray GRB subset (the only set both arms predict).
#'
#' Fidelity notes (deliberate, to mirror superlearner.R):
#'   - MICE runs once on the full table before CV (same train/test leakage the
#'     original scheme has; identical in both arms so the comparison is fair).
#'   - Cleaning rules (T90 cut, physical-value nulling) come from
#'     load_xray_modality / load_optical_modality, which replicate superlearner.R.
#'   - Learner library is the active library of superlearner.R (tuned GAM/GLM/
#'     bayesglm formulas + caret rf + generics).
#'
#' Usage:
#'   Rscript concat_superlearner.R [reps] [k_folds]      # defaults: 1, 10
#'   SMOKE_TEST=true Rscript concat_superlearner.R       # fast structural check
#'
#' Outputs: OutputFiles/ConcatExperiment/{cv_predictions.csv, metrics.csv}
#' =============================================================================

source("shared_latent_model.R")  # loaders, matching, regression_metrics
suppressPackageStartupMessages({
  library(SuperLearner)
  library(mice)
  library(caret)
  library(parallel)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")
source("Custom_SL/sl_custom_bayesglm.R")
source("Custom_SL/sl_xgboost_safe.R")

args    <- commandArgs(trailingOnly = TRUE)
reps    <- if (length(args) >= 1) as.integer(args[1]) else 1L
k_folds <- if (length(args) >= 2) as.integer(args[2]) else 10L
smoke   <- tolower(Sys.getenv("SMOKE_TEST", "false")) == "true"
n_cores <- max(1L, as.integer(Sys.getenv("CONCAT_CORES", "3")))
out_dir <- "OutputFiles/ConcatExperiment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load, clean, match -----------------------------------------------------
matched <- match_modalities(
  load_xray_modality("x-ray_data.csv"),
  load_optical_modality("optical_data.txt")
)
xray    <- matched$xray
optical <- matched$optical
cat(sprintf("Loaded %d x-ray, %d optical rows; %d pairs; %d conflicts dropped\n",
            length(xray$id), length(optical$id), nrow(matched$pairs), nrow(matched$conflicts)))

# ---- Union table (features raw, NA where modality missing) -------------------
x_feats <- xray$imputation[, xray$model_features, drop = FALSE]
o_feats <- optical$imputation[, optical$model_features, drop = FALSE]
names(o_feats) <- paste0("opt_", names(o_feats))

sample_ids <- union(xray$sample_id, optical$sample_id)
xi <- match(sample_ids, xray$sample_id)     # NA when GRB lacks x-ray
oi <- match(sample_ids, optical$sample_id)  # NA when GRB lacks optical

units <- data.frame(
  sample_id    = sample_ids,
  GRB          = ifelse(!is.na(xi), xray$id[xi], optical$id[oi]),
  availability = ifelse(!is.na(xi) & !is.na(oi), "both",
                 ifelse(!is.na(xi), "xray", "optical")),
  stringsAsFactors = FALSE
)
# Target: mean of the two (conflict-filtered, so they agree within 0.05).
yx <- xray$y[xi]; yo <- optical$y[oi]
units$y <- rowMeans(cbind(yx, yo), na.rm = TRUE)
stopifnot(!anyNA(units$y))

# data.frame[NA, ] yields an all-NA row, which is exactly what we want for the
# missing modality.
X_union <- cbind(x_feats[xi, , drop = FALSE], o_feats[oi, , drop = FALSE])
rownames(X_union) <- units$sample_id
cat(sprintf("Union table: %d GRBs (%d xray-only, %d both, %d optical-only), %d features\n",
            nrow(units), sum(units$availability == "xray"),
            sum(units$availability == "both"), sum(units$availability == "optical"),
            ncol(X_union)))

# ---- MICE (once, full data -- mirrors superlearner.R) ------------------------
mice_take_last <- function(data, m = 20L, seed = 1L) {
  set.seed(seed)
  fit <- mice(data, m = m, method = "midastouch", printFlag = FALSE)
  out <- complete(fit, m)
  stopifnot(!anyNA(out))
  out
}

is_xray_unit <- units$availability != "optical"

cat("MICE arm B (joint, union table)...\n")
Xb <- mice_take_last(X_union, seed = 1L)

cat("MICE arm A (x-ray rows/features only)...\n")
Xa_raw <- X_union[is_xray_unit, xray$model_features, drop = FALSE]
Xa <- mice_take_last(Xa_raw, seed = 1L)

append_squares <- function(d) {
  for (nm in colnames(d)) d[[paste0(nm, "Sqr")]] <- d[[nm]]^2
  d
}
Xa <- append_squares(Xa)   # 20 cols
Xb <- append_squares(Xb)   # 30 cols
rownames(Xa) <- units$sample_id[is_xray_unit]
rownames(Xb) <- units$sample_id

# ---- Learner library (identical to superlearner.R's active library) ----------
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
#' LIB_MODE=lean (default): tuned GAM/GLM/bayesglm formulas + rf + glmnet +
#'   xgboost + ranger + mean. Mirrors the small library of the historical ~0.63
#'   runs ("1GAM_1GLM_2algo") and avoids the stepwise learners, whose runtime
#'   explodes on arm B's 30 columns (~465-term interaction scope).
#' LIB_MODE=full: the complete generic library of current superlearner.R.
lib_mode <- tolower(Sys.getenv("LIB_MODE", "lean"))
generic_libs <- if (lib_mode == "full") c(
  "SL.glmnet", "SL.xgboost_safe",
  "SL.caret.rpart", "SL.earth", "SL.ipredbagg",
  "SL.mean", "SL.nnet", "SL.randomForest", "SL.ranger",
  "SL.rpart", "SL.step", "SL.step.forward",
  "SL.step.interaction", "SL.stepAIC"
) else c(
  "SL.glmnet", "SL.xgboost_safe", "SL.ranger", "SL.mean"
)
libs <- c(learner1$names, sl_glm1$names, sl_bglm$names, caret_learner$names, generic_libs)
cat("LIB_MODE:", lib_mode, "\n")
if (smoke) {
  libs    <- c("SL.mean", "SL.lm")
  k_folds <- 5L
  cat("SMOKE TEST: library reduced to", paste(libs, collapse = ", "), "\n")
}
cat("Learner library (", length(libs), "):", paste(libs, collapse = ", "), "\n")

fit_sl <- function(Y, X, newX, seed) {
  set.seed(seed)
  suppressMessages(capture.output(
    fit <- SuperLearner(Y = Y, X = X, newX = newX, family = gaussian(),
                        SL.library = libs, cvControl = list(V = 5), verbose = FALSE),
    file = nullfile()
  ))
  drop(fit$SL.predict)
}

# ---- CV ----------------------------------------------------------------------
run_fold <- function(f, folds, rep) {
  test_idx <- folds[[f]]
  train_idx <- setdiff(seq_len(nrow(units)), test_idx)

  # Arm B: union rows.
  pred_b <- fit_sl(units$y[train_idx],
                   Xb[train_idx, , drop = FALSE],
                   Xb[test_idx, , drop = FALSE],
                   seed = rep * 1000L + f)

  # Arm A: x-ray rows only, addressed by sample_id.
  tr_a_ids <- units$sample_id[intersect(train_idx, which(is_xray_unit))]
  te_a_ids <- units$sample_id[intersect(test_idx,  which(is_xray_unit))]
  pred_a <- if (length(te_a_ids)) {
    fit_sl(units$y[match(tr_a_ids, units$sample_id)],
           Xa[tr_a_ids, , drop = FALSE],
           Xa[te_a_ids, , drop = FALSE],
           seed = rep * 1000L + f)
  } else numeric(0)

  msg <- sprintf("[%s] rep %d fold %d done: %d test units (%d x-ray)\n",
                 format(Sys.time(), "%H:%M:%S"), rep, f, length(test_idx), length(te_a_ids))
  cat(msg)
  cat(msg, file = file.path(out_dir, "progress.log"), append = TRUE)
  data.frame(
    rep = rep, fold = f,
    sample_id = units$sample_id[test_idx],
    pred_b = pred_b,
    pred_a = unname(setNames(pred_a, te_a_ids)[units$sample_id[test_idx]]),
    stringsAsFactors = FALSE
  )
}

all_rows <- list()
for (rep in seq_len(reps)) {
  set.seed(2026L + rep)
  folds <- createFolds(units$y, k = k_folds)
  cat(sprintf("=== Repetition %d/%d: %d folds across %d cores ===\n",
              rep, reps, length(folds), n_cores))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(rep)
  fold_results <- mclapply(seq_along(folds), run_fold,
                           folds = folds, rep = rep, mc.cores = n_cores)
  bad <- !vapply(fold_results, is.data.frame, logical(1))
  if (any(bad)) {
    stop("Fold failure in rep ", rep, ":\n",
         paste(utils::capture.output(print(fold_results[[which(bad)[1]]])), collapse = "\n"))
  }
  all_rows[[rep]] <- do.call(rbind, fold_results)
  # Incremental save so partial results survive an interrupted run.
  partial <- merge(do.call(rbind, all_rows), units, by = "sample_id", sort = FALSE)
  write.csv(partial, file.path(out_dir, "cv_predictions_partial.csv"), row.names = FALSE)
}

predictions <- do.call(rbind, all_rows)
predictions <- merge(predictions, units, by = "sample_id", sort = FALSE)
write.csv(predictions, file.path(out_dir, "cv_predictions.csv"), row.names = FALSE)

# ---- Metrics -------------------------------------------------------------------
metric_row <- function(arm, scope, sel, pred_col) {
  d <- predictions[sel, ]
  cbind(arm = arm, scope = scope,
        regression_metrics(d$y, d[[pred_col]]))
}
xray_scope <- predictions$availability %in% c("xray", "both")
metrics <- rbind(
  metric_row("A_xray_only",  "xray_grbs",    xray_scope, "pred_a"),
  metric_row("B_concat",     "xray_grbs",    xray_scope, "pred_b"),
  metric_row("B_concat",     "all_grbs",     rep(TRUE, nrow(predictions)), "pred_b"),
  metric_row("A_xray_only",  "both_subset",  predictions$availability == "both", "pred_a"),
  metric_row("B_concat",     "both_subset",  predictions$availability == "both", "pred_b"),
  metric_row("A_xray_only",  "xronly_subset", predictions$availability == "xray", "pred_a"),
  metric_row("B_concat",     "xronly_subset", predictions$availability == "xray", "pred_b"),
  metric_row("B_concat",     "optical_only", predictions$availability == "optical", "pred_b")
)
# Per-repetition r on the primary scope, to see split-to-split variance.
for (rep in unique(predictions$rep)) {
  sel <- xray_scope & predictions$rep == rep
  metrics <- rbind(metrics,
                   metric_row(sprintf("A_xray_only_rep%d", rep), "xray_grbs", sel, "pred_a"),
                   metric_row(sprintf("B_concat_rep%d",   rep), "xray_grbs", sel, "pred_b"))
}
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
cat("\n==== RESULTS (log10(z+1) unless _z) ====\n")
print(metrics, row.names = FALSE, digits = 3)
cat("\nWrote", file.path(out_dir, "cv_predictions.csv"), "and metrics.csv\n")
