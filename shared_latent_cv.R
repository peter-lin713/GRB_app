#!/usr/bin/env Rscript

# Grouped CV for affine and nonlinear modality encoders with a shared redshift
# head. Example:
#   Rscript shared_latent_cv.R --dimensions=1,2,3,5 \
#     --models=affine,nonlinear --align_weight=1
# Every option uses --name=value syntax; defaults are defined below.

source("shared_latent_model.R")

parse_arguments <- function(args) {
  defaults <- list(
    xray = "x-ray_data.csv", optical = "optical_data.txt",
    output = "OutputFiles/SharedLatent", folds = 5L,
    dimensions = c(1L, 2L, 3L, 5L), models = c("affine", "nonlinear"),
    mice_m = 5L, mice_maxit = 5L, epochs = 1500L, patience = 150L,
    hidden_dim = 8L, align_weight = 1, l2 = 1e-3,
    learning_rate = 0.01, seed = 42L
  )
  for (arg in args) {
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L || !parts[[1L]] %in% names(defaults)) stop("Unknown argument: ", arg)
    name <- parts[[1L]]
    value <- parts[[2L]]
    defaults[[name]] <- switch(name,
      folds =, mice_m =, mice_maxit =, epochs =, patience =, hidden_dim =, seed = as.integer(value),
      align_weight =, l2 =, learning_rate = as.numeric(value),
      dimensions = as.integer(strsplit(value, ",", fixed = TRUE)[[1L]]),
      models = strsplit(value, ",", fixed = TRUE)[[1L]],
      value
    )
  }
  defaults
}

metric_table <- function(predictions) {
  keys <- unique(predictions[c("model", "latent_dim", "modality")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    selected <- predictions$model == keys$model[[i]] &
      predictions$latent_dim == keys$latent_dim[[i]] &
      predictions$modality == keys$modality[[i]]
    cbind(keys[i, , drop = FALSE], regression_metrics(
      predictions$y_true[selected], predictions$y_pred[selected]
    ))
  })
  combined_keys <- unique(predictions[c("model", "latent_dim")])
  combined <- lapply(seq_len(nrow(combined_keys)), function(i) {
    selected <- predictions$model == combined_keys$model[[i]] &
      predictions$latent_dim == combined_keys$latent_dim[[i]]
    cbind(combined_keys[i, , drop = FALSE], modality = "combined",
          regression_metrics(predictions$y_true[selected], predictions$y_pred[selected]))
  })
  do.call(rbind, c(rows, combined))
}

run_shared_latent_cv <- function(config) {
  if (!all(config$models %in% c("affine", "nonlinear"))) stop("models must be affine and/or nonlinear")
  dir.create(config$output, recursive = TRUE, showWarnings = FALSE)

  xray <- load_xray_modality(config$xray)
  optical <- load_optical_modality(config$optical)
  matched <- match_modalities(xray, optical)
  xray <- matched$xray
  optical <- matched$optical
  cat(sprintf("Loaded %d X-ray and %d optical rows; %d valid pairs; %d target conflicts excluded\n",
              length(xray$id), length(optical$id), nrow(matched$pairs), nrow(matched$conflicts)))
  if (nrow(matched$conflicts)) {
    write.csv(matched$conflicts, file.path(config$output, "target_conflicts.csv"), row.names = FALSE)
  }

  all_ids <- c(xray$sample_id, optical$sample_id)
  all_y <- c(xray$y, optical$y)
  fold_map <- make_grouped_folds(all_ids, all_y, config$folds, config$seed)
  x_fold <- unname(fold_map[xray$sample_id])
  o_fold <- unname(fold_map[optical$sample_id])
  stopifnot(!anyNA(x_fold), !anyNA(o_fold))

  prediction_rows <- list()
  diagnostic_rows <- list()
  output_index <- 0L
  diagnostic_index <- 0L

  for (fold in seq_len(max(fold_map))) {
    cat(sprintf("Fold %d/%d: fitting modality-local MICE models\n", fold, max(fold_map)))
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
    yx_train <- (xray$y[x_train] - y_center) / y_scale
    yo_train <- (optical$y[o_train] - y_center) / y_scale
    Xx_train <- x_scaled$matrix[x_train, , drop = FALSE]
    Xo_train <- o_scaled$matrix[o_train, , drop = FALSE]
    Xx_test <- x_scaled$matrix[!x_train, , drop = FALSE]
    Xo_test <- o_scaled$matrix[!o_train, , drop = FALSE]

    for (kind in config$models) for (dimension in config$dimensions) {
      cat(sprintf("  %s latent_dim=%d\n", kind, dimension))
      model <- fit_shared_latent(
        Xx_train, yx_train, Xo_train, yo_train,
        xray$sample_id[x_train], optical$sample_id[o_train],
        latent_dim = dimension, kind = kind, hidden_dim = config$hidden_dim,
        align_weight = config$align_weight, l2 = config$l2,
        learning_rate = config$learning_rate, epochs = config$epochs,
        patience = config$patience,
        seed = config$seed + fold * 1000L + dimension * 10L + match(kind, config$models)
      )
      pred_x <- predict(model, Xx_test, "xray")
      pred_o <- predict(model, Xo_test, "optical")
      pred_x_y <- pred_x$prediction * y_scale + y_center
      pred_o_y <- pred_o$prediction * y_scale + y_center

      output_index <- output_index + 1L
      prediction_rows[[output_index]] <- data.frame(
        model = kind, latent_dim = dimension, fold = fold, modality = "xray",
        sample_id = xray$sample_id[!x_train], GRB = xray$id[!x_train],
        paired = grepl("^P:", xray$sample_id[!x_train]),
        y_true = xray$y[!x_train], y_pred = pred_x_y,
        z_true = xray$z[!x_train], z_pred = 10^pred_x_y - 1
      )
      output_index <- output_index + 1L
      prediction_rows[[output_index]] <- data.frame(
        model = kind, latent_dim = dimension, fold = fold, modality = "optical",
        sample_id = optical$sample_id[!o_train], GRB = optical$id[!o_train],
        paired = grepl("^P:", optical$sample_id[!o_train]),
        y_true = optical$y[!o_train], y_pred = pred_o_y,
        z_true = optical$z[!o_train], z_pred = 10^pred_o_y - 1
      )

      shared_test <- intersect(xray$sample_id[!x_train], optical$sample_id[!o_train])
      tx <- match(shared_test, xray$sample_id[!x_train])
      to <- match(shared_test, optical$sample_id[!o_train])
      test_alignment <- if (length(shared_test)) {
        mean((pred_x$latent[tx, , drop = FALSE] - pred_o$latent[to, , drop = FALSE])^2)
      } else NA_real_
      diagnostic_index <- diagnostic_index + 1L
      diagnostic_rows[[diagnostic_index]] <- data.frame(
        model = kind, latent_dim = dimension, fold = fold,
        train_loss = model$train_loss, supervised_loss = model$supervised_loss,
        train_alignment_mse = model$alignment_loss,
        train_pairs = model$paired_n, test_pairs = length(shared_test),
        test_alignment_mse = test_alignment
      )
    }
  }

  predictions <- do.call(rbind, prediction_rows)
  diagnostics <- do.call(rbind, diagnostic_rows)
  metrics <- metric_table(predictions)
  write.csv(predictions, file.path(config$output, "cv_predictions.csv"), row.names = FALSE)
  write.csv(metrics, file.path(config$output, "cv_metrics.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(config$output, "alignment_diagnostics.csv"), row.names = FALSE)
  write.csv(matched$pairs, file.path(config$output, "matched_pairs.csv"), row.names = FALSE)
  saveRDS(list(config = config, fold_map = fold_map), file.path(config$output, "cv_run.rds"))
  print(metrics[order(metrics$rmse_y), ], row.names = FALSE)
  invisible(list(predictions = predictions, metrics = metrics, diagnostics = diagnostics))
}

if (sys.nframe() == 0L) run_shared_latent_cv(parse_arguments(commandArgs(trailingOnly = TRUE)))
