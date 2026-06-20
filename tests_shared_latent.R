source("shared_latent_model.R")

failed <- FALSE
check <- function(condition, message) {
  if (isTRUE(condition)) cat("PASS:", message, "\n") else {
    cat("FAIL:", message, "\n")
    failed <<- TRUE
  }
}

optical <- read_optical_catalog("optical_data.txt")
check(nrow(optical) == 179L, "optical parser retains all 179 catalog rows")
check(abs(optical$z[[1L]] - 0.83) < 1e-10, "optical parser handles extra metadata fields")
check(abs(optical$logLumTa[[1L]] - 43.46001832) < 1e-8, "optical parser aligns trailing numeric fields")

xray_modality <- load_xray_modality("x-ray_data.csv")
optical_modality <- load_optical_modality("optical_data.txt")
check(!any(c("ktotal", "dk", "logTarescaled", "logLumTa") %in% optical_modality$model_features),
      "redshift-derived optical columns are excluded")
matched <- match_modalities(xray_modality, optical_modality)
check(nrow(matched$conflicts) == 2L, "two conflicting paired redshifts are detected")
check(all(c("081203A", "081029") %in% matched$conflicts$x_id),
      "GRB081203A and GRB081029 are the conflicting targets")
check(!anyDuplicated(matched$pairs$x_id) && !anyDuplicated(matched$pairs$o_id),
      "cross-modal pairs are one-to-one")

set.seed(7)
n <- 120L
latent <- matrix(rnorm(n * 2L), n, 2L)
Xx <- cbind(latent, matrix(rnorm(n * 3L, sd = 0.2), n, 3L))
Xo <- cbind(latent %*% matrix(c(1, 0.4, -0.2, 0.8), 2L),
             matrix(rnorm(n * 2L, sd = 0.2), n, 2L))
y <- drop(latent %*% c(0.7, -0.4) + rnorm(n, sd = 0.05))
ids <- paste0("P:", seq_len(n))
train <- seq_len(90L)
test <- 91:120
model <- fit_shared_latent(
  Xx[train, ], y[train], Xo[train, ], y[train], ids[train], ids[train],
  latent_dim = 2L, kind = "affine", align_weight = 1, l2 = 1e-3,
  learning_rate = 0.01, epochs = 1000L, patience = 150L, seed = 11L
)
px <- predict(model, Xx[test, ], "xray")$prediction
po <- predict(model, Xo[test, ], "optical")$prediction
baseline_rmse <- sqrt(mean((y[test] - mean(y[train]))^2))
check(sqrt(mean((y[test] - px)^2)) < baseline_rmse * 0.5, "affine X-ray branch learns shared target")
check(sqrt(mean((y[test] - po)^2)) < baseline_rmse * 0.5, "affine optical branch learns shared target")

gradient_params <- initialize_joint_model(5L, 4L, 2L, "nonlinear", hidden_dim = 3L, seed = 13L)
gradient_fit <- joint_loss_gradient(
  gradient_params, Xx[1:5, ], y[1:5], Xo[1:5, ], y[1:5],
  1:5, 1:5, "nonlinear", align_weight = 0.7, l2 = 1e-3
)
epsilon <- 1e-6
plus <- minus <- gradient_params
plus$x_W1[1, 1] <- plus$x_W1[1, 1] + epsilon
minus$x_W1[1, 1] <- minus$x_W1[1, 1] - epsilon
numeric_gradient <- (
  joint_loss_gradient(plus, Xx[1:5, ], y[1:5], Xo[1:5, ], y[1:5], 1:5, 1:5,
                      "nonlinear", 0.7, 1e-3)$loss -
  joint_loss_gradient(minus, Xx[1:5, ], y[1:5], Xo[1:5, ], y[1:5], 1:5, 1:5,
                      "nonlinear", 0.7, 1e-3)$loss
) / (2 * epsilon)
check(abs(numeric_gradient - gradient_fit$gradients$x_W1[1, 1]) < 1e-5,
      "nonlinear backpropagation matches a finite-difference gradient")

deduplicated <- deduplicate_latent_views(
  matrix(c(1, 3, 2, 4), 2L, 2L), c("P:A", "X:B"), c("A", "B"), c(0.2, 0.3),
  matrix(c(3, 5, 8, 10), 2L, 2L), c("P:A", "O:C"), c("A", "C"), c(0.2, 0.4)
)
check(nrow(deduplicated) == 3L && !anyDuplicated(deduplicated$sample_id),
      "latent views are deduplicated to one row per GRB")
check(all(as.numeric(deduplicated["P:A", c("latent_1", "latent_2")]) == c(2, 5)),
      "paired latent vectors are averaged")
check(identical(deduplicated["P:A", "availability"], "both"),
      "deduplicated paired rows retain availability metadata")

folds <- make_grouped_folds(c(ids, ids), c(y, y), k = 5L, seed = 3L)
check(length(folds) == n && all(folds >= 1L & folds <= 5L), "grouped fold map has one fold per GRB")

cat(if (failed) "\nSOME TESTS FAILED\n" else "\nALL TESTS PASSED\n")
if (failed) quit(status = 1L)
