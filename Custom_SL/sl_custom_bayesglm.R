#' sl_custom_bayesglm.R — Bayesian GLM base learner for SuperLearner.
#'
#' A SuperLearner wrapper around arm::bayesglm (Bayesian generalized linear
#' regression with weakly-informative priors). By default it regresses the
#' response on all predictors additively; a pre-built formula may be supplied
#' via `bglm.model`, in which case its right-hand side is reused with the
#' response rebound to the local Y. Registered under the name
#' "SL.custom_bayesglm"; the companion predict method scores new data.
require('arm')

#' Fit a Bayesian GLM base learner.
#' @param Y numeric. Training response.
#' @param X data.frame. Training predictors.
#' @param newX data.frame. Predictors to score after fitting.
#' @param family family. GLM family object (e.g. gaussian()).
#' @param obsWeights numeric. Per-observation weights.
#' @param bglm.model formula or NA. Optional pre-built formula; if NA (default)
#'   an additive formula over all columns of X is generated.
#' @return list(pred = numeric, fit = <SL.custom_bayesglm object>).
SL.custom_bayesglm <- function(Y, X, newX, family, obsWeights, bglm.model = NA, ...){

  if(is.atomic(bglm.model)){
    bglm.model <- as.formula(paste("Y~", paste(colnames(X), collapse="+")))
  } else {
    # A formula was supplied: rewrite it so the response is the local Y.
    fn <- strsplit(as.character(bglm.model),split = '~')
    bglm.model <- as.formula(paste('Y ~',fn[[3]][1]))
  }

  fit.glm <- arm::bayesglm(bglm.model, data = X, family = family, weights = obsWeights)

  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.custom_bayesglm")
  return(out)
}

#' Prediction method for SL.custom_bayesglm fits.
#' @param object SL.custom_bayesglm fit (carries the fitted bayesglm in $object).
#' @param newdata data.frame. Predictors to score.
#' @return numeric. Predicted responses.
predict.SL.custom_bayesglm <- function(object, newdata, ...){
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  return(pred)
}
