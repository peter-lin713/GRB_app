#' sl_mgcv_gam.R — generalized additive model (GAM) base learner for SuperLearner.
#'
#' A SuperLearner wrapper around mgcv::gam fitting degree-2 smoothing splines.
#' Continuous predictors (more than `cts.num` unique values) enter through
#' spline terms s(x, deg.gam); low-cardinality predictors enter linearly. A
#' pre-built formula may be supplied via `gam.model` instead of auto-generating
#' one. Easy to clone for other spline degrees, e.g.:
#'   SL.gam.3 <- function(..., deg.gam = 3) SL.mgcv_gam(..., deg.gam = deg.gam)
#'
#' Registered with SuperLearner under the name "SL.mgcv_gam"; the companion
#' predict.SL.mgcv_gam method scores new data.

#' Fit a degree-`deg.gam` GAM base learner.
#' @param Y numeric. Training response.
#' @param X data.frame. Training predictors.
#' @param newX data.frame. Predictors to score after fitting.
#' @param family family. GLM family object (e.g. gaussian()).
#' @param obsWeights numeric. Per-observation weights.
#' @param deg.gam integer. Spline degree for continuous terms (default 2).
#' @param cts.num integer. Min unique values for a predictor to be treated as
#'   continuous and given a spline term (default 4).
#' @param gam.model formula or NA. Optional pre-built formula; if NA (default)
#'   a formula is generated automatically.
#' @param verbose logical. If TRUE, print the formula and fitting progress.
#' @return list(pred = numeric, fit = <SL.mgcv_gam object>).
SL.mgcv_gam <- function(Y, X, newX, family, obsWeights, deg.gam = 2, cts.num = 4, gam.model = NA, verbose = F, ...) {
  # mgcv::s() must resolve as a special function when the formula is parsed, so
  # the package is loaded via require() (gam::s() is not recognized here).

  # A predictor is "continuous" (eligible for a spline) if it has more than
  # cts.num distinct values.
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if(is.atomic(gam.model)){
    if (sum(!cts.x) > 0) {
      # Mix of continuous (spline) and discrete (linear) predictors.
      gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam,")", sep=""), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop=FALSE]), collapse = "+")))
    } else {
      # All predictors continuous: spline terms only.
      gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""), collapse = "+")))
    }
    # All predictors discrete/binomial: purely linear model.
    if (sum(!cts.x) == length(cts.x)) {
      gam.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+"), sep = ""))
    }
  } else {
    # A formula was supplied: rewrite it so the response is the local Y.
    if(verbose){print(gam.model)}
    fn <- strsplit(as.character(gam.model),split = '~')
    gam.model <- as.formula(paste('Y ~',fn[[3]][1]))
  }

  if(verbose){print(gam.model)}

  fit.gam <- mgcv::gam(gam.model, data = X, family = family, control = mgcv::gam.control(maxit = 50), weights = obsWeights)

  if(verbose){print('fitted')}

  if(packageVersion('gam') >= 1.15) {
    pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response") # updated gam class in version 1.15
  } else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.mgcv_gam")
  return(out)
}

#' Prediction method for SL.mgcv_gam fits.
#' @param object SL.mgcv_gam fit (carries the fitted mgcv::gam in $object).
#' @param newdata data.frame. Predictors to score.
#' @return numeric. Predicted responses.
predict.SL.mgcv_gam <- function(object, newdata, ...){
  if(packageVersion('gam') >= 1.15) {
    pred <- mgcv::predict.gam(object = object$object, newdata = newdata, type = "response") # updated gam class in version 1.15
  } else {
    stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
  }

  return(pred)
}
