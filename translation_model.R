#' =============================================================================
#' translation_model.R -- optical -> x-ray feature-space translation models
#' =============================================================================
#'
#' Unlike the shared-latent models (symmetric encoders into an abstract
#' bottleneck), the translation approach is asymmetric and keeps the x-ray
#' branch untouched:
#'
#'   - X-ray GRBs keep their full feature set; no bottleneck, no tax.
#'   - Optical features are mapped INTO x-ray feature space by a translation
#'     model trained on the paired GRBs, for which the target (the GRB's real
#'     x-ray features) is observed. The mapping has a supervised, physical
#'     target rather than an abstract latent objective.
#'
#' Two variants:
#'
#'   1. fit_ridge_translation / predict_translation -- two-stage. Multi-target
#'      ridge regression from (scaled) optical inputs to the 10 (scaled) x-ray
#'      features, ridge penalty chosen by exact leave-one-out CV on the pairs.
#'      Output feeds the existing SuperLearner pipeline as pseudo x-ray rows.
#'
#'   2. fit_end_to_end_translation -- joint. The translation matrix and a
#'      quadratic ridge head ([x, x^2] %*% h + b) are trained together by Adam:
#'      prediction loss flows from y all the way through the head into the
#'      translation weights (the "backprop all the way" requirement), plus a
#'      feature-matching loss on the pairs and L2 regularization. Everything is
#'      differentiable; gradients are analytic and finite-difference checked in
#'      tests_translation.R.
#'
#' All functions operate in SCALED space (caller standardizes features and y).
#' =============================================================================

source("shared_latent_model.R")  # loaders, matching, mice_complete_fold, scale_fold_features, regression_metrics

# ---- Variant 1: two-stage ridge translation ---------------------------------

#' Fit a multi-target ridge translation O -> X on paired GRBs.
#'
#' @param O numeric matrix, n_pairs x p_o. Scaled optical inputs (caller may
#'   append squared columns for mild nonlinearity).
#' @param X numeric matrix, n_pairs x p_x. Scaled x-ray feature targets.
#' @param lambdas numeric vector of ridge penalties to score by exact LOO-CV.
#' @return list(W, lambda, loo_mse, loo_mse_by_lambda). W is (p_o + 1) x p_x
#'   including the intercept row; the LOO residual uses the shared hat matrix
#'   e_loo = e / (1 - h_ii), valid columnwise because all targets share O.
fit_ridge_translation <- function(O, X, lambdas = 10^seq(-3, 3, length.out = 31)) {
  O <- as.matrix(O); X <- as.matrix(X)
  stopifnot(nrow(O) == nrow(X), nrow(O) > ncol(O) + 1L)
  O1 <- cbind(intercept = 1, O)
  best <- NULL
  scores <- numeric(length(lambdas))
  for (k in seq_along(lambdas)) {
    penalty <- diag(ncol(O1)) * lambdas[[k]]
    penalty[1L, 1L] <- 0  # never penalize the intercept
    A_inv <- solve(crossprod(O1) + penalty)
    W <- A_inv %*% crossprod(O1, X)
    hat_diag <- rowSums((O1 %*% A_inv) * O1)
    hat_diag <- pmin(hat_diag, 1 - 1e-6)
    loo_resid <- (X - O1 %*% W) / (1 - hat_diag)
    scores[[k]] <- mean(loo_resid^2)
    if (is.null(best) || scores[[k]] < best$loo_mse) {
      best <- list(W = W, lambda = lambdas[[k]], loo_mse = scores[[k]])
    }
  }
  best$loo_mse_by_lambda <- setNames(scores, signif(lambdas, 3))
  best
}

#' Translate optical inputs into x-ray feature space.
#' @param fit result of fit_ridge_translation.
#' @param O numeric matrix, n x p_o (same input convention as the fit).
#' @return numeric matrix, n x p_x of translated (scaled) x-ray features.
predict_translation <- function(fit, O) {
  cbind(1, as.matrix(O)) %*% fit$W
}

# ---- Variant 2: end-to-end translation + quadratic ridge head ----------------

#' Quadratic basis: [x, x^2]. Differentiable; gives the head the squared terms
#' the tuned GAM/GLM formulas rely on.
quad_basis <- function(X) cbind(X, X^2)

#' Loss and analytic gradients for the end-to-end model.
#'
#' Parameters: W (p_o x p_x) + c (1 x p_x) translation; h (2*p_x) + b head.
#' Loss = mean squared prediction error over (x-ray rows on REAL features +
#'        optical-only rows on TRANSLATED features)
#'      + match_weight * mean((translate(O_pair) - X_pair)^2)
#'      + l2 * (sum(W^2) + sum(h^2)).
#' Prediction error backpropagates through the head into W via the optical-only
#' rows; the match loss anchors the translation to observed x-ray features.
end_to_end_loss_grad <- function(params, X_real, y_real, O_only, y_only,
                                 O_pair, X_pair, match_weight, l2) {
  W <- params$W; c_vec <- params$c; h <- params$h; b <- params$b
  p_x <- ncol(W)
  h1 <- h[seq_len(p_x)]; h2 <- h[p_x + seq_len(p_x)]

  X_hat <- sweep(O_only %*% W, 2L, c_vec, "+")          # translated features
  pred_real <- drop(quad_basis(X_real) %*% h) + b
  pred_only <- drop(quad_basis(X_hat) %*% h) + b
  n_pred <- length(y_real) + length(y_only)
  r_real <- pred_real - y_real
  r_only <- pred_only - y_only
  pred_loss <- (sum(r_real^2) + sum(r_only^2)) / n_pred

  M <- sweep(O_pair %*% W, 2L, c_vec, "+") - X_pair      # match residuals
  match_loss <- mean(M^2)

  loss <- pred_loss + match_weight * match_loss + l2 * (sum(W^2) + sum(h^2))

  # Head gradients.
  d_h <- 2 * (crossprod(quad_basis(X_real), r_real) +
              crossprod(quad_basis(X_hat),  r_only)) / n_pred + 2 * l2 * h
  d_b <- 2 * (sum(r_real) + sum(r_only)) / n_pred

  # Prediction loss into translated features: d pred / d x = h1 + 2 x h2.
  d_Xhat <- (2 / n_pred) * (tcrossprod(r_only, h1) +
                            2 * tcrossprod(r_only, h2) * X_hat)
  d_M <- (2 * match_weight / length(M)) * M

  d_W <- crossprod(O_only, d_Xhat) + crossprod(O_pair, d_M) + 2 * l2 * W
  d_c <- colSums(d_Xhat) + colSums(d_M)

  list(loss = loss, pred_loss = pred_loss, match_loss = match_loss,
       gradients = list(W = d_W, c = matrix(d_c, 1L), h = drop(as.matrix(d_h)), b = d_b))
}

#' Train the end-to-end model with Adam (full batch) and early stopping.
#'
#' @param X_real n_x x p_x scaled x-ray features (x-ray + paired GRBs), with
#'   y_real their scaled targets.
#' @param O_only n_o x p_o scaled optical inputs for optical-ONLY GRBs, with
#'   y_only their scaled targets.
#' @param O_pair,X_pair the paired GRBs' optical inputs and x-ray features
#'   (the feature-matching supervision).
#' @return object of class end_to_end_translation with $params and diagnostics.
fit_end_to_end_translation <- function(X_real, y_real, O_only, y_only,
                                       O_pair, X_pair,
                                       match_weight = 1, l2 = 1e-3,
                                       learning_rate = 0.01, epochs = 2000L,
                                       patience = 200L, seed = 1L,
                                       init_translation = NULL, verbose = FALSE) {
  X_real <- as.matrix(X_real); O_only <- as.matrix(O_only)
  O_pair <- as.matrix(O_pair); X_pair <- as.matrix(X_pair)
  p_x <- ncol(X_real); p_o <- ncol(O_only)
  stopifnot(ncol(O_pair) == p_o, ncol(X_pair) == p_x,
            nrow(O_pair) == nrow(X_pair), nrow(O_only) == length(y_only),
            nrow(X_real) == length(y_real))

  set.seed(seed)
  params <- list(
    W = if (is.null(init_translation)) {
      matrix(rnorm(p_o * p_x, sd = sqrt(2 / (p_o + p_x))), p_o, p_x)
    } else init_translation$W[-1L, , drop = FALSE],
    c = if (is.null(init_translation)) matrix(0, 1L, p_x)
        else matrix(init_translation$W[1L, ], 1L),
    h = rnorm(2L * p_x, sd = 0.1),
    b = 0
  )

  m1 <- lapply(params, function(x) x * 0)
  m2 <- lapply(params, function(x) x * 0)
  beta1 <- 0.9; beta2 <- 0.999; eps <- 1e-8
  best <- params; best_loss <- Inf; stale <- 0L

  for (epoch in seq_len(epochs)) {
    cur <- end_to_end_loss_grad(params, X_real, y_real, O_only, y_only,
                                O_pair, X_pair, match_weight, l2)
    if (cur$loss < best_loss - 1e-8) {
      best_loss <- cur$loss; best <- params; stale <- 0L
    } else {
      stale <- stale + 1L
      if (stale >= patience) break
    }
    g <- cur$gradients
    g_norm <- sqrt(sum(vapply(g, function(x) sum(x^2), numeric(1))))
    if (is.finite(g_norm) && g_norm > 10) g <- lapply(g, function(x) x * 10 / g_norm)
    for (nm in names(params)) {
      m1[[nm]] <- beta1 * m1[[nm]] + (1 - beta1) * g[[nm]]
      m2[[nm]] <- beta2 * m2[[nm]] + (1 - beta2) * g[[nm]]^2
      params[[nm]] <- params[[nm]] - learning_rate *
        (m1[[nm]] / (1 - beta1^epoch)) / (sqrt(m2[[nm]] / (1 - beta2^epoch)) + eps)
    }
    if (verbose && epoch %% 200L == 0L) {
      cat(sprintf("epoch=%d loss=%.5f pred=%.5f match=%.5f\n",
                  epoch, cur$loss, cur$pred_loss, cur$match_loss))
    }
  }
  final <- end_to_end_loss_grad(best, X_real, y_real, O_only, y_only,
                                O_pair, X_pair, match_weight, l2)
  structure(list(params = best, loss = final$loss, pred_loss = final$pred_loss,
                 match_loss = final$match_loss),
            class = "end_to_end_translation")
}

#' Predict from the end-to-end model.
#' @param newdata scaled features; x-ray feature matrix when modality="xray"
#'   (real features go straight to the head), optical input matrix when
#'   modality="optical" (translated first).
#' @return list(prediction, translated). `translated` is NULL for x-ray input.
predict.end_to_end_translation <- function(object, newdata, modality = c("xray", "optical"), ...) {
  modality <- match.arg(modality)
  newdata <- as.matrix(newdata)
  X <- if (modality == "xray") newdata else
    sweep(newdata %*% object$params$W, 2L, drop(object$params$c), "+")
  list(prediction = drop(quad_basis(X) %*% object$params$h) + object$params$b,
       translated = if (modality == "optical") X else NULL)
}
