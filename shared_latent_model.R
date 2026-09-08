# Shared-latent multimodal redshift models and preprocessing utilities.

normalize_grb_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("^GRB", "", x)
}

grb_date_root <- function(x) sub("^([0-9]{6})[A-Z]$", "\\1", normalize_grb_id(x))

read_optical_catalog <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2L) stop("Optical catalog has no data rows: ", path)

  fields <- strsplit(trimws(lines[-1L]), "[[:space:]]+")
  tail_names <- c(
    "z", "T90", "Class", "logFa", "logFaErr", "logTa", "logTaErr",
    "Alpha", "AlphaErr", "beta", "betaErr", "ktotal", "dk",
    "logTarescaled", "logTaErrRescaled", "logLumTa", "logLumTaErr"
  )
  if (any(lengths(fields) < length(tail_names) + 1L)) {
    stop("Optical catalog contains a row with too few fields")
  }

  rows <- lapply(fields, function(parts) {
    c(GRB = parts[[1L]], setNames(tail(parts, length(tail_names)), tail_names))
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  numeric_names <- setdiff(names(out), c("GRB", "Class"))
  out[numeric_names] <- lapply(out[numeric_names], function(x) suppressWarnings(as.numeric(x)))
  out$GRB <- normalize_grb_id(out$GRB)
  if (anyNA(out$GRB) || anyDuplicated(out$GRB)) stop("Optical GRB identifiers are missing or duplicated")
  out
}

to_dex_error <- function(linear_error, log10_value) {
  value <- linear_error / (10^log10_value * log(10))
  value[!is.finite(value)] <- NA_real_
  value
}

load_xray_modality <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  id <- normalize_grb_id(if ("GRB" %in% names(raw)) raw$GRB else raw[[1L]])
  if (anyNA(id) || anyDuplicated(id)) stop("X-ray GRB identifiers are missing or duplicated")
  if (!"log10T90" %in% names(raw)) raw$log10T90 <- log10(raw$T90)

  raw$log10PeakFlux[!is.finite(raw$log10PeakFlux)] <- NA_real_
  raw$log10NH[raw$log10NH < 20] <- NA_real_
  raw$Beta[raw$Beta > 3] <- NA_real_
  raw$Gamma[raw$Gamma > 3] <- NA_real_
  raw$Alpha[raw$Alpha > 3] <- NA_real_
  raw$PhotonIndex[raw$PhotonIndex < 0] <- NA_real_

  model_features <- c(
    "log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
    "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux"
  )
  required <- c("Redshift_crosscheck", model_features)
  missing <- setdiff(required, names(raw))
  if (length(missing)) stop("X-ray catalog is missing: ", paste(missing, collapse = ", "))

  imputation <- raw[model_features]
  error_columns <- list(
    log10T90Err = to_dex_error(raw$T90Err, raw$log10T90),
    log10FaErr = raw$log10FaErr,
    log10TaErr = raw$log10TaErr,
    AlphaErr = raw$AlphaErr,
    BetaErr = raw$BetaErr,
    log10FluenceErr = to_dex_error(raw$FluenceErr, raw$log10Fluence),
    PhotonIndexErr = raw$PhotonIndexErr,
    log10PeakFluxErr = to_dex_error(raw$PeakFluxErr, raw$log10PeakFlux)
  )
  for (name in names(error_columns)) imputation[[name]] <- error_columns[[name]]

  keep <- is.na(raw$log10T90) | raw$log10T90 > log10(2)
  list(
    modality = "xray",
    id = id[keep],
    z = as.numeric(raw$Redshift_crosscheck[keep]),
    y = log10(as.numeric(raw$Redshift_crosscheck[keep]) + 1),
    imputation = imputation[keep, , drop = FALSE],
    model_features = model_features
  )
}

load_optical_modality <- function(path) {
  raw <- read_optical_catalog(path)
  raw$log10T90 <- log10(raw$T90)
  raw$log10T90[!is.finite(raw$log10T90)] <- NA_real_

  # Rest-frame time, luminosity, and k-correction columns depend on redshift and
  # are deliberately excluded. Only observed-frame measurements are modeled.
  model_features <- c("log10T90", "logFa", "logTa", "Alpha", "beta")
  error_features <- c("logFaErr", "logTaErr", "AlphaErr", "betaErr")
  imputation <- raw[c(model_features, error_features)]
  keep <- is.na(raw$log10T90) | raw$log10T90 > log10(2)

  list(
    modality = "optical",
    id = raw$GRB[keep],
    z = as.numeric(raw$z[keep]),
    y = log10(as.numeric(raw$z[keep]) + 1),
    imputation = imputation[keep, , drop = FALSE],
    model_features = model_features
  )
}

match_modalities <- function(xray, optical, conflict_tolerance = 0.05) {
  exact_ids <- intersect(xray$id, optical$id)
  pairs <- data.frame(
    x_id = exact_ids, o_id = exact_ids, match_type = rep("exact", length(exact_ids)),
    stringsAsFactors = FALSE
  )

  x_left <- setdiff(xray$id, exact_ids)
  o_left <- setdiff(optical$id, exact_ids)
  x_root <- grb_date_root(x_left)
  o_root <- grb_date_root(o_left)
  shared_roots <- intersect(unique(x_root), unique(o_root))
  for (root in shared_roots) {
    xi <- x_left[x_root == root]
    oi <- o_left[o_root == root]
    if (length(xi) == 1L && length(oi) == 1L) {
      pairs <- rbind(pairs, data.frame(x_id = xi, o_id = oi, match_type = "date_root"))
    }
  }
  if (anyDuplicated(pairs$x_id) || anyDuplicated(pairs$o_id)) stop("Non-unique cross-modal match")

  pairs$x_row <- match(pairs$x_id, xray$id)
  pairs$o_row <- match(pairs$o_id, optical$id)
  pairs$x_z <- xray$z[pairs$x_row]
  pairs$o_z <- optical$z[pairs$o_row]
  pairs$target_difference <- abs(pairs$x_z - pairs$o_z)
  pairs$target_conflict <- !is.finite(pairs$target_difference) |
    pairs$target_difference > conflict_tolerance

  xray$sample_id <- paste0("X:", xray$id)
  optical$sample_id <- paste0("O:", optical$id)
  for (i in seq_len(nrow(pairs))) {
    key <- paste0("P:", pairs$x_id[[i]], "|", pairs$o_id[[i]])
    xray$sample_id[pairs$x_row[[i]]] <- key
    optical$sample_id[pairs$o_row[[i]]] <- key
  }

  conflicts <- pairs[pairs$target_conflict, , drop = FALSE]
  if (nrow(conflicts)) {
    bad_keys <- paste0("P:", conflicts$x_id, "|", conflicts$o_id)
    x_keep <- !xray$sample_id %in% bad_keys
    o_keep <- !optical$sample_id %in% bad_keys
    xray <- lapply_modality_rows(xray, x_keep)
    optical <- lapply_modality_rows(optical, o_keep)
    pairs <- pairs[!pairs$target_conflict, , drop = FALSE]
    pairs$x_row <- match(pairs$x_id, xray$id)
    pairs$o_row <- match(pairs$o_id, optical$id)
  }
  list(xray = xray, optical = optical, pairs = pairs, conflicts = conflicts)
}

lapply_modality_rows <- function(modality, keep) {
  modality$id <- modality$id[keep]
  modality$z <- modality$z[keep]
  modality$y <- modality$y[keep]
  modality$sample_id <- modality$sample_id[keep]
  modality$imputation <- modality$imputation[keep, , drop = FALSE]
  modality
}

mice_complete_fold <- function(data, train_rows, m = 5L, maxit = 5L, seed = 1L) {
  data <- as.data.frame(data)
  if (!all(vapply(data, is.numeric, logical(1)))) stop("MICE input must be numeric")
  if (!anyNA(data)) return(data)
  if (!requireNamespace("mice", quietly = TRUE)) stop("Package 'mice' is required")

  usable <- vapply(data, function(col) any(is.finite(col[train_rows])), logical(1))
  if (!all(usable)) data <- data[usable]
  if (!ncol(data)) stop("No features have observed training values")

  method <- mice::make.method(data)
  has_missing <- vapply(data, anyNA, logical(1))
  method[has_missing] <- "midastouch"
  method[!has_missing] <- ""
  predictor_matrix <- mice::make.predictorMatrix(data)
  diag(predictor_matrix) <- 0
  constants <- vapply(data, function(col) {
    values <- col[train_rows & is.finite(col)]
    length(unique(values)) < 2L
  }, logical(1))
  predictor_matrix[, constants] <- 0

  fit <- mice::mice(
    data, m = as.integer(m), maxit = as.integer(maxit), method = method,
    predictorMatrix = predictor_matrix, ignore = !train_rows,
    printFlag = FALSE, seed = as.integer(seed),
    remove.constant = FALSE, remove.collinear = FALSE
  )
  completed <- lapply(seq_len(m), function(i) mice::complete(fit, i))
  pooled <- Reduce(`+`, completed) / length(completed)
  pooled <- as.data.frame(pooled)
  if (anyNA(pooled)) stop("MICE left missing values in completed fold data")
  pooled
}

scale_fold_features <- function(completed, train_rows, feature_names) {
  feature_names <- intersect(feature_names, names(completed))
  means <- vapply(completed[train_rows, feature_names, drop = FALSE], mean, numeric(1))
  sds <- vapply(completed[train_rows, feature_names, drop = FALSE], sd, numeric(1))
  keep <- is.finite(means) & is.finite(sds) & sds > 1e-8
  if (!any(keep)) stop("No non-constant model features remain in this fold")
  feature_names <- feature_names[keep]
  matrix <- sweep(as.matrix(completed[feature_names]), 2L, means[keep], "-")
  matrix <- sweep(matrix, 2L, sds[keep], "/")
  list(matrix = matrix, features = feature_names, center = means[keep], scale = sds[keep])
}

make_grouped_folds <- function(sample_ids, outcomes, k = 5L, seed = 42L) {
  units <- unique(sample_ids)
  unit_y <- vapply(units, function(id) mean(outcomes[sample_ids == id], na.rm = TRUE), numeric(1))
  if (any(!is.finite(unit_y))) stop("Every CV unit must have a finite outcome")
  k <- min(as.integer(k), length(units))
  set.seed(seed)
  order_y <- order(unit_y, runif(length(unit_y)))
  unit_fold <- integer(length(units))
  chunks <- split(order_y, ceiling(seq_along(order_y) / k))
  for (chunk in chunks) unit_fold[chunk] <- sample(seq_along(chunk))
  setNames(unit_fold, units)
}

initialize_joint_model <- function(px, po, latent_dim, kind = c("affine", "nonlinear"),
                                   hidden_dim = 8L, seed = 1L) {
  kind <- match.arg(kind)
  set.seed(seed)
  init <- function(n_in, n_out) matrix(rnorm(n_in * n_out, sd = sqrt(2 / (n_in + n_out))), n_in, n_out)
  d <- as.integer(latent_dim)
  if (kind == "affine") {
    return(list(
      x_W = init(px, d), x_b = matrix(0, 1, d),
      o_W = init(po, d), o_b = matrix(0, 1, d),
      h_W = init(d, 1), h_b = matrix(0, 1, 1)
    ))
  }
  h <- as.integer(hidden_dim)
  list(
    x_W1 = init(px, h), x_b1 = matrix(0, 1, h), x_W2 = init(h, d), x_b2 = matrix(0, 1, d),
    o_W1 = init(po, h), o_b1 = matrix(0, 1, h), o_W2 = init(h, d), o_b2 = matrix(0, 1, d),
    h_W = init(d, 1), h_b = matrix(0, 1, 1)
  )
}

encoder_forward <- function(X, params, prefix, kind) {
  if (kind == "affine") {
    z <- sweep(X %*% params[[paste0(prefix, "_W")]], 2L,
               drop(params[[paste0(prefix, "_b")]]), "+")
    return(list(z = z, X = X))
  }
  h <- tanh(sweep(X %*% params[[paste0(prefix, "_W1")]], 2L,
                  drop(params[[paste0(prefix, "_b1")]]), "+"))
  z <- tanh(sweep(h %*% params[[paste0(prefix, "_W2")]], 2L,
                  drop(params[[paste0(prefix, "_b2")]]), "+"))
  list(z = z, h = h, X = X)
}

encoder_backward <- function(cache, dz, params, prefix, kind) {
  if (kind == "affine") {
    return(setNames(list(crossprod(cache$X, dz), matrix(colSums(dz), 1L)),
                    c(paste0(prefix, "_W"), paste0(prefix, "_b"))))
  }
  dz_pre <- dz * (1 - cache$z^2)
  dW2 <- crossprod(cache$h, dz_pre)
  db2 <- matrix(colSums(dz_pre), 1L)
  dh <- dz_pre %*% t(params[[paste0(prefix, "_W2")]])
  dh_pre <- dh * (1 - cache$h^2)
  setNames(
    list(crossprod(cache$X, dh_pre), matrix(colSums(dh_pre), 1L), dW2, db2),
    paste0(prefix, c("_W1", "_b1", "_W2", "_b2"))
  )
}

joint_loss_gradient <- function(params, Xx, yx, Xo, yo, pair_x, pair_o,
                                kind, align_weight, l2) {
  fx <- encoder_forward(Xx, params, "x", kind)
  fo <- encoder_forward(Xo, params, "o", kind)
  px <- drop(fx$z %*% params$h_W) + drop(params$h_b)
  po <- drop(fo$z %*% params$h_W) + drop(params$h_b)
  n_total <- length(yx) + length(yo)
  rx <- px - yx
  ro <- po - yo
  loss <- (sum(rx^2) + sum(ro^2)) / n_total
  dpx <- 2 * rx / n_total
  dpo <- 2 * ro / n_total
  dzx <- tcrossprod(dpx, drop(params$h_W))
  dzo <- tcrossprod(dpo, drop(params$h_W))

  align_loss <- 0
  if (length(pair_x)) {
    difference <- fx$z[pair_x, , drop = FALSE] - fo$z[pair_o, , drop = FALSE]
    align_loss <- mean(difference^2)
    align_gradient <- 2 * align_weight * difference / length(difference)
    dzx[pair_x, ] <- dzx[pair_x, , drop = FALSE] + align_gradient
    dzo[pair_o, ] <- dzo[pair_o, , drop = FALSE] - align_gradient
    loss <- loss + align_weight * align_loss
  }

  grads <- c(
    encoder_backward(fx, dzx, params, "x", kind),
    encoder_backward(fo, dzo, params, "o", kind)
  )
  grads$h_W <- crossprod(fx$z, dpx) + crossprod(fo$z, dpo)
  grads$h_b <- matrix(sum(dpx) + sum(dpo), 1L, 1L)

  weight_names <- grep("_W", names(params), value = TRUE)
  reg_loss <- sum(vapply(params[weight_names], function(x) sum(x^2), numeric(1)))
  loss <- loss + l2 * reg_loss
  for (name in weight_names) grads[[name]] <- grads[[name]] + 2 * l2 * params[[name]]
  list(loss = loss, supervised_loss = (sum(rx^2) + sum(ro^2)) / n_total,
       alignment_loss = align_loss, gradients = grads)
}

fit_shared_latent <- function(Xx, yx, Xo, yo, x_ids, o_ids, latent_dim,
                              kind = c("affine", "nonlinear"), hidden_dim = 8L,
                              align_weight = 1, l2 = 1e-3, learning_rate = 0.01,
                              epochs = 1500L, patience = 150L, seed = 1L,
                              verbose = FALSE) {
  kind <- match.arg(kind)
  paired <- intersect(x_ids, o_ids)
  pair_x <- match(paired, x_ids)
  pair_o <- match(paired, o_ids)
  params <- initialize_joint_model(ncol(Xx), ncol(Xo), latent_dim, kind, hidden_dim, seed)
  first_moment <- lapply(params, function(x) x * 0)
  second_moment <- lapply(params, function(x) x * 0)
  best <- params
  best_loss <- Inf
  stale <- 0L
  beta1 <- 0.9
  beta2 <- 0.999
  eps <- 1e-8

  for (epoch in seq_len(epochs)) {
    current <- joint_loss_gradient(params, Xx, yx, Xo, yo, pair_x, pair_o,
                                   kind, align_weight, l2)
    if (current$loss < best_loss - 1e-8) {
      best_loss <- current$loss
      best <- params
      stale <- 0L
    } else {
      stale <- stale + 1L
      if (stale >= patience) break
    }
    grad_norm <- sqrt(sum(vapply(current$gradients, function(x) sum(x^2), numeric(1))))
    if (is.finite(grad_norm) && grad_norm > 10) {
      current$gradients <- lapply(current$gradients, function(x) x * 10 / grad_norm)
    }
    for (name in names(params)) {
      gradient <- current$gradients[[name]]
      first_moment[[name]] <- beta1 * first_moment[[name]] + (1 - beta1) * gradient
      second_moment[[name]] <- beta2 * second_moment[[name]] + (1 - beta2) * gradient^2
      m_hat <- first_moment[[name]] / (1 - beta1^epoch)
      v_hat <- second_moment[[name]] / (1 - beta2^epoch)
      params[[name]] <- params[[name]] - learning_rate * m_hat / (sqrt(v_hat) + eps)
    }
    if (verbose && epoch %% 100L == 0L) {
      cat(sprintf("epoch=%d loss=%.6f supervised=%.6f alignment=%.6f\n",
                  epoch, current$loss, current$supervised_loss, current$alignment_loss))
    }
  }
  final <- joint_loss_gradient(best, Xx, yx, Xo, yo, pair_x, pair_o,
                               kind, align_weight, l2)
  structure(list(
    params = best, kind = kind, latent_dim = latent_dim, hidden_dim = hidden_dim,
    train_loss = final$loss, supervised_loss = final$supervised_loss,
    alignment_loss = final$alignment_loss, paired_n = length(paired)
  ), class = "shared_latent_model")
}

predict.shared_latent_model <- function(object, newdata, modality = c("xray", "optical"), ...) {
  modality <- match.arg(modality)
  prefix <- if (modality == "xray") "x" else "o"
  encoded <- encoder_forward(as.matrix(newdata), object$params, prefix, object$kind)$z
  prediction <- drop(encoded %*% object$params$h_W) + drop(object$params$h_b)
  list(prediction = prediction, latent = encoded)
}

deduplicate_latent_views <- function(x_latent, x_sample_id, x_grb, x_y,
                                     o_latent, o_sample_id, o_grb, o_y,
                                     target_tolerance = 0.05) {
  x_latent <- as.matrix(x_latent)
  o_latent <- as.matrix(o_latent)
  if (ncol(x_latent) != ncol(o_latent)) stop("Latent dimensions differ across modalities")
  if (anyDuplicated(x_sample_id) || anyDuplicated(o_sample_id)) {
    stop("Each modality must contain at most one row per sample_id")
  }

  sample_ids <- c(x_sample_id, setdiff(o_sample_id, x_sample_id))
  latent_names <- paste0("latent_", seq_len(ncol(x_latent)))
  rows <- lapply(sample_ids, function(sample_id) {
    xi <- match(sample_id, x_sample_id)
    oi <- match(sample_id, o_sample_id)
    has_x <- !is.na(xi)
    has_o <- !is.na(oi)
    views <- rbind(
      if (has_x) x_latent[xi, , drop = FALSE],
      if (has_o) o_latent[oi, , drop = FALSE]
    )
    targets <- c(if (has_x) x_y[[xi]], if (has_o) o_y[[oi]])
    if (length(targets) == 2L && abs(diff(targets)) > target_tolerance) {
      stop("Target disagreement remains after matching for ", sample_id)
    }
    grbs <- c(if (has_x) x_grb[[xi]], if (has_o) o_grb[[oi]])
    data.frame(
      sample_id = sample_id,
      GRB = paste(unique(grbs), collapse = "|"),
      availability = if (has_x && has_o) "both" else if (has_x) "xray" else "optical",
      y = mean(targets),
      view_distance = if (nrow(views) == 2L) mean((views[1L, ] - views[2L, ])^2) else NA_real_,
      as.data.frame(as.list(setNames(colMeans(views), latent_names))),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- result$sample_id
  if (anyDuplicated(result$sample_id) || anyNA(result[c("sample_id", "y", latent_names)])) {
    stop("Invalid deduplicated latent table")
  }
  result
}

regression_metrics <- function(observed_y, predicted_y) {
  keep <- is.finite(observed_y) & is.finite(predicted_y)
  observed_y <- observed_y[keep]
  predicted_y <- predicted_y[keep]
  observed_z <- 10^observed_y - 1
  predicted_z <- 10^predicted_y - 1
  safe_cor <- function(x, y) if (length(x) > 2L && sd(x) > 0 && sd(y) > 0) cor(x, y) else NA_real_
  data.frame(
    n = length(observed_y),
    rmse_y = sqrt(mean((observed_y - predicted_y)^2)),
    mae_y = mean(abs(observed_y - predicted_y)),
    correlation_y = safe_cor(observed_y, predicted_y),
    rmse_z = sqrt(mean((observed_z - predicted_z)^2)),
    mae_z = mean(abs(observed_z - predicted_z)),
    correlation_z = safe_cor(observed_z, predicted_z)
  )
}
