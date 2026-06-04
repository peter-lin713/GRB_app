#' sl_custom_glm.R — plain GLM base learner for SuperLearner.
#'
#' A thin SuperLearner wrapper around stats::glm. By default it regresses the
#' response on all predictors additively; a pre-built formula may be supplied
#' via `glm.model`, in which case its right-hand side is reused with the
#' response rebound to the local Y. Registered under the name "SL.custom_glm";
#' the companion predict.SL.custom_glm method scores new data.
require('arm')

#' Fit a GLM base learner.
#' @param Y numeric. Training response.
#' @param X data.frame. Training predictors.
#' @param newX data.frame. Predictors to score after fitting.
#' @param family family. GLM family object (e.g. gaussian()).
#' @param obsWeights numeric. Per-observation weights.
#' @param glm.model formula or NA. Optional pre-built formula; if NA (default)
#'   an additive formula over all columns of X is generated.
#' @return list(pred = numeric, fit = <SL.custom_glm object>).
SL.custom_glm <- function(Y, X, newX, family, obsWeights, glm.model = NA, ...){

  if(is.atomic(glm.model)){
    glm.model <- as.formula(paste("Y~", paste(colnames(X), collapse="+")))
  } else {
    # A formula was supplied: rewrite it so the response is the local Y.
    fn <- strsplit(as.character(glm.model),split = '~')
    glm.model <- as.formula(paste('Y ~',fn[[3]][1]))
  }

  fit.glm <- glm(glm.model, data = X, family = family, weights = obsWeights)

  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.custom_glm")
  return(out)
}

#' Prediction method for SL.custom_glm fits.
#' @param object SL.custom_glm fit (carries the fitted glm in $object).
#' @param newdata data.frame. Predictors to score.
#' @return numeric. Predicted responses.
predict.SL.custom_glm <- function(object, newdata, ...){
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  return(pred)
}
