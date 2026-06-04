#' sl_xgboost_safe.R — version-proof xgboost base learner for SuperLearner.
#'
#' SuperLearner::SL.xgboost calls the high-level xgboost::xgboost(), whose
#' signature changed in xgboost >= 2.1.0 (the old data=/label= interface became
#' a required positional x, y). On those versions the stock wrapper dies with:
#'   argument "y" is missing, with no default
#' xgb.train() has a stable interface across xgboost 1.x and 2.x, so we call it
#' directly. The returned fit carries class "SL.xgboost" so the existing
#' predict.SL.xgboost method scores it unchanged.

#' Fit an xgboost base learner via the stable xgb.train() interface.
#' @param Y numeric. Training response.
#' @param X matrix or data.frame. Training predictors (auto-converted to a
#'   numeric design matrix if not already a matrix).
#' @param newX matrix or data.frame. Predictors to score after fitting.
#' @param family family. GLM family; supports gaussian() and binomial().
#' @param obsWeights numeric. Per-observation weights.
#' @param id ignored. Present for SuperLearner wrapper-signature compatibility.
#' @param ntrees integer. Number of boosting rounds (default 1000).
#' @param max_depth integer. Max tree depth (default 4).
#' @param shrinkage numeric. Learning rate / eta (default 0.1).
#' @param minobspernode integer. Min child weight (default 10).
#' @param params list. Extra xgboost params, merged with the defaults above.
#' @param nthread integer. Threads for xgboost (default 1).
#' @param verbose integer. xgboost verbosity (default 0, silent).
#' @return list(pred = numeric, fit = <SL.xgboost object>).
SL.xgboost_safe <- function(Y, X, newX, family, obsWeights, id,
                            ntrees = 1000, max_depth = 4, shrinkage = 0.1,
                            minobspernode = 10, params = list(), nthread = 1,
                            verbose = 0, ...) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("SL.xgboost_safe requires the xgboost package")
  }
  if (!is.matrix(X)) X <- model.matrix(~ . - 1, X)

  objective <- switch(family$family,
    gaussian = "reg:squarederror",
    binomial = "binary:logistic",
    stop("SL.xgboost_safe: unsupported family ", family$family))

  xgmat <- xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  params_full <- c(params, list(
    objective        = objective,
    max_depth        = max_depth,
    eta              = shrinkage,
    min_child_weight = minobspernode,
    nthread          = nthread
  ))
  model <- xgboost::xgb.train(params = params_full, data = xgmat,
                              nrounds = ntrees, verbose = verbose)

  if (!is.matrix(newX)) newX <- model.matrix(~ . - 1, newX)
  pred <- predict(model, newdata = newX)

  fit <- list(object = model)
  class(fit) <- "SL.xgboost"   # reuse SuperLearner's predict.SL.xgboost
  list(pred = pred, fit = fit)
}
