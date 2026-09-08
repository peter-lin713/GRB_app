source("translation_model.R")

failed <- FALSE
check <- function(condition, message) {
  if (isTRUE(condition)) cat("PASS:", message, "\n") else {
    cat("FAIL:", message, "\n")
    failed <<- TRUE
  }
}

# ---- Ridge translation: recovers a known linear map --------------------------
set.seed(21)
n <- 90L; p_o <- 6L; p_x <- 10L
O <- matrix(rnorm(n * p_o), n, p_o)
T_true <- matrix(rnorm(p_o * p_x, sd = 0.6), p_o, p_x)
c_true <- rnorm(p_x, sd = 0.3)
X <- sweep(O %*% T_true, 2L, c_true, "+") + matrix(rnorm(n * p_x, sd = 0.1), n, p_x)
tr_fit <- fit_ridge_translation(O[1:60, ], X[1:60, ])
X_hat <- predict_translation(tr_fit, O[61:90, ])
oos_mse <- mean((X_hat - X[61:90, ])^2)
base_mse <- mean(sweep(X[61:90, ], 2L, colMeans(X[1:60, ]))^2)
check(all(dim(tr_fit$W) == c(p_o + 1L, p_x)), "ridge translation W has intercept + p_o rows")
check(oos_mse < 0.05 * base_mse, "ridge translation recovers a linear map out of sample")
check(tr_fit$loo_mse < base_mse, "LOO-selected lambda beats the mean predictor")

# ---- End-to-end gradients match finite differences ----------------------------
set.seed(33)
n_x <- 12L; n_o <- 7L; n_p <- 9L
X_real <- matrix(rnorm(n_x * p_x), n_x, p_x)
y_real <- rnorm(n_x)
O_only <- matrix(rnorm(n_o * p_o), n_o, p_o)
y_only <- rnorm(n_o)
O_pair <- matrix(rnorm(n_p * p_o), n_p, p_o)
X_pair <- matrix(rnorm(n_p * p_x), n_p, p_x)
params <- list(
  W = matrix(rnorm(p_o * p_x, sd = 0.4), p_o, p_x),
  c = matrix(rnorm(p_x, sd = 0.2), 1L),
  h = rnorm(2L * p_x, sd = 0.3),
  b = 0.15
)
analytic <- end_to_end_loss_grad(params, X_real, y_real, O_only, y_only,
                                 O_pair, X_pair, match_weight = 0.7, l2 = 1e-3)
eps <- 1e-6
max_err <- 0
probes <- list(c("W", 1L), c("W", 17L), c("c", 3L), c("h", 2L), c("h", 14L), c("b", 1L))
for (probe in probes) {
  name <- probe[[1L]]; idx <- as.integer(probe[[2L]])
  plus <- minus <- params
  plus[[name]][idx]  <- plus[[name]][idx] + eps
  minus[[name]][idx] <- minus[[name]][idx] - eps
  numeric_grad <- (
    end_to_end_loss_grad(plus,  X_real, y_real, O_only, y_only, O_pair, X_pair, 0.7, 1e-3)$loss -
    end_to_end_loss_grad(minus, X_real, y_real, O_only, y_only, O_pair, X_pair, 0.7, 1e-3)$loss
  ) / (2 * eps)
  max_err <- max(max_err, abs(numeric_grad - analytic$gradients[[name]][idx]))
}
check(max_err < 1e-5, sprintf("end-to-end gradients match finite differences (max err %.2e)", max_err))

# ---- End-to-end learns on synthetic data --------------------------------------
set.seed(44)
n2 <- 160L
O2 <- matrix(rnorm(n2 * p_o), n2, p_o)
X2 <- sweep(O2 %*% T_true, 2L, c_true, "+") + matrix(rnorm(n2 * p_x, sd = 0.1), n2, p_x)
h_true <- rnorm(p_x, sd = 0.4)
y2 <- drop(X2 %*% h_true) + 0.2 * drop(X2[, 1L]^2) + rnorm(n2, sd = 0.1)
xr <- 1:70; oo <- 71:110; pr <- 111:160                 # xray rows / optical-only / pairs
warm <- fit_ridge_translation(O2[pr, ], X2[pr, ])
e2e <- fit_end_to_end_translation(
  X_real = X2[xr, ], y_real = y2[xr],
  O_only = O2[oo, ], y_only = y2[oo],
  O_pair = O2[pr, ], X_pair = X2[pr, ],
  match_weight = 1, seed = 5L, init_translation = warm
)
test_o <- predict(e2e, O2[oo, ], "optical")$prediction
test_x <- predict(e2e, X2[xr, ], "xray")$prediction
base_rmse <- sqrt(mean((y2[oo] - mean(y2[xr]))^2))
check(sqrt(mean((y2[oo] - test_o)^2)) < 0.5 * base_rmse,
      "end-to-end model predicts optical-only rows well above baseline")
check(sqrt(mean((y2[xr] - test_x)^2)) < 0.5 * base_rmse,
      "end-to-end model fits x-ray rows through the quadratic head")
check(!is.null(predict(e2e, O2[oo, ], "optical")$translated) &&
      is.null(predict(e2e, X2[xr, ], "xray")$translated),
      "predict returns translated features only for optical input")
check(e2e$match_loss < 0.2, "warm-started translation keeps feature matching tight")

cat(if (failed) "\nSOME TESTS FAILED\n" else "\nALL TESTS PASSED\n")
if (failed) quit(status = 1L)
