#!/usr/bin/env Rscript

# Nested experiment: learn modality encoders inside each outer fold, average
# paired GRBs into one aligned latent vector, then train SuperLearner on the
# deduplicated latent table. The outer holdout is never used by MICE, scaling,
# encoder fitting, latent alignment, or SuperLearner fitting.

source("shared_latent_model.R")
source("Custom_SL/sl_xgboost_safe.R")
suppressPackageStartupMessages(library(SuperLearner))

parse_latent_sl_arguments <- function(args) {
  defaults <- list(
    xray = "x-ray_data.csv", optical = "optical_data.txt",
    output = "OutputFiles/SharedLatentSuperLearner", folds = 5L,
    dimensions = c(1L, 2L, 3L, 5L), models = c("affine", "nonlinear"),
    mice_m = 5L, mice_maxit = 5L, epochs = 1500L, patience = 150L,
    hidden_dim = 8L, align_weight = 1, l2 = 1e-3,
    learning_rate = 0.01, inner_folds = 5L, sl_mode = "full", seed = 42L
  )
  for (arg in args) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L || !parts[[1L]] %in% names(defaults)) stop("Unknown argument: ", arg)
    name <- parts[[1L]]
    value <- parts[[2L]]
    defaults[[name]] <- switch(name,
      folds =, mice_m =, mice_maxit =, epochs =, patience =, hidden_dim =,
      inner_folds =, seed = as.integer(value),
      align_weight =, l2 =, learning_rate = as.numeric(value),
      dimensions = as.integer(strsplit(value, ",", fixed = TRUE)[[1L]]),
      models = strsplit(value, ",", fixed = TRUE)[[1L]],
      value
    )
  }
  defaults
}

latent_sl_library <- function(latent_dim, mode = c("full", "smoke")) {
  mode <- match.arg(mode)
  if (mode == "smoke") return(c("SL.mean", "SL.glm"))
  learners <- c(
    "SL.glm", "SL.glmnet", "SL.xgboost_safe",
    "SL.caret.rpart", "SL.earth", "SL.ipredbagg",
    "SL.mean", "SL.nnet", "SL.randomForest", "SL.ranger",
    "SL.rpart", "SL.step", "SL.step.forward",
    "SL.step.interaction", "SL.stepAIC"
  )
  # glmnet rejects a one-column design matrix. The remaining learners all
  # support a scalar latent representation.
  if (latent_dim == 1L) learners <- setdiff(learners, "SL.glmnet")
  learners
}

fit_latent_superlearner <- function(train_table, test_table, inner_folds,
                                    latent_dim, sl_mode, seed) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) stop("Package 'SuperLearner' is required")
  latent_names <- paste0("latent_", seq_len(latent_dim))
  train_x <- train_table[latent_names]
  test_x <- test_table[latent_names]
  stopifnot(!anyNA(train_x), !anyNA(test_x), !anyNA(train_table$y))
  set.seed(seed)
  warning_messages <- character()
  fit <- withCallingHandlers({
    capture.output(
      sl_fit <- SuperLearner::SuperLearner(
        Y = train_table$y,
        X = train_x,
        newX = test_x,
        family = gaussian(),
        SL.library = latent_sl_library(latent_dim, sl_mode),
        cvControl = list(V = min(as.integer(inner_folds), nrow(train_table))),
        verbose = FALSE
      ),
      file = nullfile()
    )
    sl_fit
  }, warning = function(warning) {
    warning_messages <<- c(warning_messages, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  list(
    prediction = drop(fit$SL.predict),
    coefficients = setNames(as.numeric(fit$coef), names(fit$coef)),
    risks = setNames(as.numeric(fit$cvRisk), names(fit$cvRisk)),
    warnings = unique(warning_messages)
  )
}

summarize_latent_sl_metrics <- function(predictions, prediction_column, predictor_name) {
  keys <- unique(predictions[c("encoder", "latent_dim", "availability")])
  by_availability <- lapply(seq_len(nrow(keys)), function(i) {
    selected <- predictions$encoder == keys$encoder[[i]] &
      predictions$latent_dim == keys$latent_dim[[i]] &
      predictions$availability == keys$availability[[i]]
    cbind(predictor = predictor_name, keys[i, , drop = FALSE], regression_metrics(
      predictions$y_true[selected], predictions[[prediction_column]][selected]
    ))
  })
  overall_keys <- unique(predictions[c("encoder", "latent_dim")])
  overall <- lapply(seq_len(nrow(overall_keys)), function(i) {
    selected <- predictions$encoder == overall_keys$encoder[[i]] &
      predictions$latent_dim == overall_keys$latent_dim[[i]]
    cbind(predictor = predictor_name, overall_keys[i, , drop = FALSE], availability = "all",
          regression_metrics(predictions$y_true[selected], predictions[[prediction_column]][selected]))
  })
  do.call(rbind, c(by_availability, overall))
}

run_latent_superlearner_cv <- function(config) {
  if (!all(config$models %in% c("affine", "nonlinear"))) stop("models must be affine and/or nonlinear")
  if (!config$sl_mode %in% c("full", "smoke")) stop("sl_mode must be full or smoke")
  dir.create(config$output, recursive = TRUE, showWarnings = FALSE)

  matched <- match_modalities(
    load_xray_modality(config$xray),
    load_optical_modality(config$optical)
  )
  xray <- matched$xray
  optical <- matched$optical
  all_ids <- c(xray$sample_id, optical$sample_id)
  all_y <- c(xray$y, optical$y)
  fold_map <- make_grouped_folds(all_ids, all_y, config$folds, config$seed)
  x_fold <- unname(fold_map[xray$sample_id])
  o_fold <- unname(fold_map[optical$sample_id])
  cat(sprintf("Loaded %d X-ray and %d optical rows -> %d unique GRBs (%d paired)\n",
              length(xray$id), length(optical$id), length(fold_map), nrow(matched$pairs)))

  predictions <- list()
  learner_rows <- list()
  diagnostics <- list()
  prediction_index <- learner_index <- diagnostic_index <- 0L

  for (fold in seq_len(max(fold_map))) {
    cat(sprintf("Outer fold %d/%d: modality-local MICE\n", fold, max(fold_map)))
    x_train <- x_fold != fold
    o_train <- o_fold != fold
    x_complete <- mice_complete_fold(
      xray$imputation, x_train, config$mice_m, config$mice_maxit,
      config$seed + fold * 100L + 1L
    )
    o_complete <- mice_complete_fold(
      optical$imputation, o_train, config$mice_m, config$mice_maxit,
      config$seed + fold * 100L + 2L
    )
    x_scaled <- scale_fold_features(x_complete, x_train, xray$model_features)
    o_scaled <- scale_fold_features(o_complete, o_train, optical$model_features)
    y_train_all <- c(xray$y[x_train], optical$y[o_train])
    y_center <- mean(y_train_all)
    y_scale <- sd(y_train_all)

    for (kind in config$models) for (dimension in config$dimensions) {
      cat(sprintf("  encoder=%s latent_dim=%d: fit encoder + SuperLearner\n", kind, dimension))
      encoder <- fit_shared_latent(
        x_scaled$matrix[x_train, , drop = FALSE],
        (xray$y[x_train] - y_center) / y_scale,
        o_scaled$matrix[o_train, , drop = FALSE],
        (optical$y[o_train] - y_center) / y_scale,
        xray$sample_id[x_train], optical$sample_id[o_train],
        latent_dim = dimension, kind = kind, hidden_dim = config$hidden_dim,
        align_weight = config$align_weight, l2 = config$l2,
        learning_rate = config$learning_rate, epochs = config$epochs,
        patience = config$patience,
        seed = config$seed + fold * 1000L + dimension * 10L + match(kind, config$models)
      )

      x_latent_train <- predict(encoder, x_scaled$matrix[x_train, , drop = FALSE], "xray")$latent
      o_latent_train <- predict(encoder, o_scaled$matrix[o_train, , drop = FALSE], "optical")$latent
      x_latent_test <- predict(encoder, x_scaled$matrix[!x_train, , drop = FALSE], "xray")$latent
      o_latent_test <- predict(encoder, o_scaled$matrix[!o_train, , drop = FALSE], "optical")$latent
      train_table <- deduplicate_latent_views(
        x_latent_train, xray$sample_id[x_train], xray$id[x_train], xray$y[x_train],
        o_latent_train, optical$sample_id[o_train], optical$id[o_train], optical$y[o_train]
      )
      test_table <- deduplicate_latent_views(
        x_latent_test, xray$sample_id[!x_train], xray$id[!x_train], xray$y[!x_train],
        o_latent_test, optical$sample_id[!o_train], optical$id[!o_train], optical$y[!o_train]
      )
      stopifnot(length(intersect(train_table$sample_id, test_table$sample_id)) == 0L)

      sl <- fit_latent_superlearner(
        train_table, test_table, config$inner_folds, dimension, config$sl_mode,
        config$seed + fold * 10000L + dimension * 100L + match(kind, config$models)
      )
      latent_names <- paste0("latent_", seq_len(dimension))
      head_y_pred <- (
        drop(as.matrix(test_table[latent_names]) %*% encoder$params$h_W) +
          drop(encoder$params$h_b)
      ) * y_scale + y_center
      prediction_index <- prediction_index + 1L
      predictions[[prediction_index]] <- data.frame(
        encoder = kind, latent_dim = dimension, fold = fold,
        sample_id = test_table$sample_id, GRB = test_table$GRB,
        availability = test_table$availability,
        y_true = test_table$y, y_pred = sl$prediction,
        head_y_pred = head_y_pred,
        z_true = 10^test_table$y - 1, z_pred = 10^sl$prediction - 1,
        head_z_pred = 10^head_y_pred - 1,
        view_distance = test_table$view_distance,
        stringsAsFactors = FALSE
      )
      learner_names <- union(names(sl$coefficients), names(sl$risks))
      learner_index <- learner_index + 1L
      learner_rows[[learner_index]] <- data.frame(
        encoder = kind, latent_dim = dimension, fold = fold,
        learner = learner_names,
        coefficient = unname(sl$coefficients[learner_names]),
        cv_risk = unname(sl$risks[learner_names]),
        stringsAsFactors = FALSE
      )
      diagnostic_index <- diagnostic_index + 1L
      diagnostics[[diagnostic_index]] <- data.frame(
        encoder = kind, latent_dim = dimension, fold = fold,
        train_n = nrow(train_table), test_n = nrow(test_table),
        train_both = sum(train_table$availability == "both"),
        test_both = sum(test_table$availability == "both"),
        encoder_train_loss = encoder$train_loss,
        encoder_alignment_mse = encoder$alignment_loss,
        test_view_distance = mean(test_table$view_distance, na.rm = TRUE),
        sl_warning_count = length(sl$warnings),
        sl_warnings = paste(sl$warnings, collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  }

  predictions <- do.call(rbind, predictions)
  learner_results <- do.call(rbind, learner_rows)
  diagnostics <- do.call(rbind, diagnostics)
  metrics <- rbind(
    summarize_latent_sl_metrics(predictions, "y_pred", "superlearner"),
    summarize_latent_sl_metrics(predictions, "head_y_pred", "shared_linear_head")
  )
  write.csv(predictions, file.path(config$output, "cv_predictions.csv"), row.names = FALSE)
  write.csv(metrics, file.path(config$output, "cv_metrics.csv"), row.names = FALSE)
  write.csv(learner_results, file.path(config$output, "learner_weights.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(config$output, "fold_diagnostics.csv"), row.names = FALSE)
  write.csv(matched$pairs, file.path(config$output, "matched_pairs.csv"), row.names = FALSE)
  if (nrow(matched$conflicts)) {
    write.csv(matched$conflicts, file.path(config$output, "target_conflicts.csv"), row.names = FALSE)
  }
  saveRDS(list(config = config, fold_map = fold_map), file.path(config$output, "cv_run.rds"))
  overall <- metrics[metrics$availability == "all", ]
  print(overall[order(overall$rmse_y), ], row.names = FALSE)
  invisible(list(predictions = predictions, metrics = metrics,
                 learner_results = learner_results, diagnostics = diagnostics))
}

if (sys.nframe() == 0L) {
  run_latent_superlearner_cv(parse_latent_sl_arguments(commandArgs(trailingOnly = TRUE)))
}
