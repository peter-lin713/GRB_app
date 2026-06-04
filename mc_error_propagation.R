#' mc_error_propagation.R — Monte Carlo error propagation for SL predictions.
#'
#' Propagates per-GRB measurement uncertainties through a trained SuperLearner
#' model. Each observed feature is treated as the mean of a Gaussian whose
#' standard deviation is the catalog error; we draw `n` perturbed feature sets,
#' predict log10(z+1) for each, and summarise the resulting predictive
#' distribution into a point estimate plus a confidence interval.
#'
#' Pipeline of helpers (top to bottom):
#'   linear_to_log10_error  scale-convert a linear error to log10 scale
#'   add_error              draw one noisy realisation of a measurement
#'   prepare_errors         assemble the per-feature error table from raw data
#'   prepare_features       assemble the model-ready feature table (+ squared terms)
#'   default_feature_error_map  feature -> error-column mapping
#'   mc_predict             run the MC loop, returns an (nGRB x n) prediction matrix
#'   summarize_mc           collapse that matrix into means / sds / quantile intervals
#'
#' Consumed downstream by the MC plotting code in mc_plots.R (mc_predict +
#' summarize_mc feed make_mc_plots()).

library(SuperLearner)

#' Convert linear-scale error to log10-scale error.
#' @param x numeric. Value(s) in linear scale.
#' @param x_err numeric. Error(s) in linear scale.
#' @return numeric. Error(s) expressed on the log10 scale.
linear_to_log10_error <- function(x, x_err) {
  x_err / (x * log(10))
}

#' Perturb a measurement by its error
#' Draws from N(measurement, |error|), treating the reported error column as a
#' 1-sigma standard deviation. NOTE: earlier versions used sd = |error/2| with no
#' documented justification -- that is only correct if the catalog errors are a
#' full-width / 2-sigma quantity. If the *Err columns are already 1-sigma (the
#' usual convention), |error/2| halves the noise and makes every MC interval too
#' narrow. Confirm the catalog's error convention; revert to error/2 only if the
#' reported errors are 2-sigma widths.
add_error <- function(measurement, error) {
  rnorm(length(measurement), mean = measurement, sd = abs(error))
}

#' Prepare error dataframe from raw data
#' Converts linear-scale errors to log10-scale where needed
#' and computes log10T90 error
#'
#' @param dat The raw dataframe (e.g. combined_data_with_redshift_V8.csv)
#' @return A dataframe with all errors in log10 or matching scale
prepare_errors <- function(dat) {
  errs <- data.frame(row.names = rownames(dat))

  # Already in matching scale — use directly
  errs$log10FaErr      <- dat$log10FaErr
  errs$log10TaErr      <- dat$log10TaErr
  errs$PhotonIndexErr  <- dat$PhotonIndexErr
  errs$AlphaErr        <- dat$AlphaErr
  errs$BetaErr         <- dat$BetaErr

  # Convert PeakFluxErr from linear to log10 scale
  errs$log10PeakFluxErr <- linear_to_log10_error(10^dat$log10PeakFlux, dat$PeakFluxErr)

  # Convert FluenceErr from linear to log10 scale
  errs$log10FluenceErr <- linear_to_log10_error(10^dat$log10Fluence, dat$FluenceErr)

  # Convert T90Err from linear to log10 scale
  errs$log10T90Err <- linear_to_log10_error(dat$T90, dat$T90Err)

  return(errs)
}

#' Prepare feature dataframe from raw data
#' Computes log10T90 and all squared terms
#'
#' @param dat The raw dataframe
#' @return A dataframe of model-ready features
prepare_features <- function(dat) {
  feats <- data.frame(row.names = rownames(dat))

  feats$log10Fa       <- dat$log10Fa
  feats$log10Ta       <- dat$log10Ta
  feats$log10NH       <- dat$log10NH
  feats$log10PeakFlux <- dat$log10PeakFlux
  feats$log10Fluence  <- dat$log10Fluence
  feats$PhotonIndex   <- dat$PhotonIndex
  feats$log10T90      <- log10(dat$T90)

  # Squared terms
  feats$log10FaSqr       <- feats$log10Fa^2
  feats$log10TaSqr       <- feats$log10Ta^2
  feats$log10NHSqr       <- feats$log10NH^2
  feats$log10PeakFluxSqr <- feats$log10PeakFlux^2
  feats$log10T90Sqr      <- feats$log10T90^2
  feats$PhotonIndexSqr   <- feats$PhotonIndex^2

  return(feats)
}

#' Default mapping from feature name to error column name.
#' Keys are feature column names in the prediction dataframe; values are the
#' corresponding error column names in the errors dataframe.
default_feature_error_map <- function() {
  list(
    log10Fa       = "log10FaErr",
    log10Ta       = "log10TaErr",
    PhotonIndex   = "PhotonIndexErr",
    log10PeakFlux = "log10PeakFluxErr",
    log10T90      = "log10T90Err",
    Alpha         = "AlphaErr",
    Beta          = "BetaErr",
    log10Fluence  = "log10FluenceErr"
  )
}

#' Monte Carlo error propagation through a trained SuperLearner model
#'
#' Perturbs each feature in `feature_error_map` by drawing from
#' N(feature, |error|) and regenerates the corresponding squared term
#' (convention: <feature>Sqr). Features in the map that are not present
#' in `features` are silently skipped so the function stays correct
#' regardless of which columns the upstream feature selection kept.
#'
#' @param model A trained SuperLearner object
#' @param features A dataframe of model-ready features (must match model's training columns)
#' @param errors A dataframe with error columns, indexed to match `features`
#' @param n Number of MC samples per GRB (default 100)
#' @param feature_error_map Named list: feature name -> error column name
#' @return A matrix of dimension nrow(features) x n, each column is one MC
#'         realization of predicted log10(z+1)
mc_predict <- function(model, features, errors, n = 100,
                       feature_error_map = default_feature_error_map()) {

  # Filter map to features actually present; warn about skipped ones
  present  <- names(feature_error_map) %in% colnames(features)
  skipped  <- names(feature_error_map)[!present]
  if (length(skipped) > 0) {
    message("mc_predict: skipping features not in data: ",
            paste(skipped, collapse = ", "))
  }
  active_map <- feature_error_map[present]

  # Validate error columns exist for the active features
  missing_err <- setdiff(unlist(active_map), colnames(errors))
  if (length(missing_err) > 0) {
    stop("mc_predict: error columns missing from errors df: ",
         paste(missing_err, collapse = ", "))
  }

  # Store base values for each active feature
  base_vals <- lapply(names(active_map), function(f) features[[f]])
  names(base_vals) <- names(active_map)

  prediction_matrix <- matrix(nrow = nrow(features), ncol = n)

  for (k in 1:n) {
    for (feat in names(active_map)) {
      err_col <- active_map[[feat]]
      features[[feat]] <- add_error(base_vals[[feat]], errors[[err_col]])

      # Regenerate squared term if it exists in the features frame
      sqr_name <- paste0(feat, "Sqr")
      if (sqr_name %in% colnames(features)) {
        features[[sqr_name]] <- features[[feat]]^2
      }
    }
    prediction_matrix[, k] <- predict(model, features)$pred
  }

  rownames(prediction_matrix) <- rownames(features)
  return(prediction_matrix)
}

#' Summarize MC prediction matrix into point estimates and intervals
#'
#' @param pred_matrix Output from mc_predict
#' @param level Confidence level for prediction intervals (default 0.95)
#' @return A dataframe with log10(z+1) and linear z summary statistics
summarize_mc <- function(pred_matrix, level = 0.95) {
  alpha <- (1 - level) / 2

  log_mean  <- rowMeans(pred_matrix)
  log_sd    <- apply(pred_matrix, 1, sd)
  log_lower <- apply(pred_matrix, 1, quantile, probs = alpha)
  log_upper <- apply(pred_matrix, 1, quantile, probs = 1 - alpha)

  linear_matrix <- 10^pred_matrix - 1
  lin_mean  <- rowMeans(linear_matrix)
  lin_sd    <- apply(linear_matrix, 1, sd)
  lin_lower <- apply(linear_matrix, 1, quantile, probs = alpha)
  lin_upper <- apply(linear_matrix, 1, quantile, probs = 1 - alpha)

  data.frame(
    log10z_mean  = log_mean,
    log10z_sd    = log_sd,
    log10z_lower = log_lower,
    log10z_upper = log_upper,
    z_mean       = lin_mean,
    z_sd         = lin_sd,
    z_lower      = lin_lower,
    z_upper      = lin_upper,
    row.names    = rownames(pred_matrix)
  )
}
