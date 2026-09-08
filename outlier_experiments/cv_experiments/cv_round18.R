# cv_round18.R — projection-uncertainty bagging: rebuild the optical rows' 4
# projected features (log10Fa, log10Ta, Alpha, Beta) under 6 emcee posterior
# calibration draws; average OOF predictions across draws.
#   anchor        median projection (current frame)
#   projbag_param 6 draws, parameter uncertainty only (m_k x + b_k)
#   projbag_full  6 draws incl. intrinsic scatter + measurement perturbation
# Round-10 semi frame base, daume, lean lib, 5f x 3 reps.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  library(SuperLearner); library(caret); library(glmnet); library(MASS)
  library(mgcv); library(gam); library(randomForest)
})
source("Custom_SL/sl_mgcv_gam.R")
source("Custom_SL/sl_custom_glm.R")

SC <- "/private/tmp/claude-501/-Users-spencergibson-Desktop-peter-code/2ffd2eec-28e2-4d34-93c6-737e2b069e0e/scratchpad"
out_prefix <- commandArgs(trailingOnly = TRUE)[1]

f0 <- read.csv(file.path(SC, "round10_frame.csv"), row.names = 1, stringsAsFactors = FALSE)
opt <- read.csv("Data/optical_data.csv", stringsAsFactors = FALSE); names(opt)[1] <- "GRB"
draws <- read.csv(file.path(SC, "proj_draws.csv"), stringsAsFactors = FALSE)
linear_vars <- c("log10PeakFlux","log10NH","log10Ta","PhotonIndex","log10Fa",
                 "Alpha","log10T90","Beta","Gamma","log10Fluence")
top7 <- c("log10NH","log10PeakFlux","PhotonIndex","log10Ta","log10Fa","Gamma","Alpha")
VARS <- c(linear_vars, paste0(top7, "Sqr"), "is_optical")
Y <- f0$log10z; N <- nrow(f0); is_opt <- f0$is_optical == 1

# optical-native measurements for the optical rows
mo <- match(f0$GRB[is_opt], opt$GRB)
onat <- data.frame(
  logFa = opt$logFa[mo], logFaErr = opt$logFaErr[mo],
  logTa = opt$logT_a[mo], logTaErr = opt$logTaErr[mo],
  Alpha = opt$Alpha[mo], AlphaErr = opt$AlphaErr[mo],
  Beta  = opt$Beta[mo],  BetaErr  = opt$betaErr[mo])
for (v in c("logFaErr","logTaErr","AlphaErr","BetaErr")) {
  onat[[v]][!is.finite(onat[[v]])] <- 1e-6
}
tgt_map <- c(logFa = "log10Fa", logTa = "log10Ta", Alpha = "Alpha", Beta = "Beta")

make_frame <- function(draw_id, full = FALSE, seed = 1) {
  f <- f0
  set.seed(seed)
  for (p in names(tgt_map)) {
    dr <- draws[draws$param == p & draws$draw == draw_id - 1, ]
    x <- onat[[p]]
    if (full) x <- x + rnorm(length(x)) * onat[[paste0(p, "Err")]]
    val <- dr$m * x + dr$b
    if (full) val <- val + rnorm(length(x)) * dr$s
    tv <- tgt_map[[p]]
    keepNA <- !is.finite(val)
    val[keepNA] <- f[[tv]][is_opt][keepNA]   # fall back to median-projection value
    f[[tv]][is_opt] <- val
  }
  for (v in top7) f[[paste0(v, "Sqr")]] <- f[[v]]^2
  f
}

fg <- read.table("Best_formula_GAM.txt"); gam_forms <- unique(apply(as.matrix(fg[, 2]), 1, as.formula))
band_form <- as.formula(paste(
  "Response ~ (log10Fa + log10Ta + log10NH + PhotonIndex + log10PeakFlux)^2 +",
  "Alpha + Beta + log10T90 + Gamma + log10Fluence +",
  "is_optical:(log10Fa + log10Ta + Alpha + Beta) + is_optical"))
gamL <- create.Learner("SL.mgcv_gam",   tune = list(gam.model = gam_forms[1]), detailed_names = FALSE, name_prefix = "gamF")
glmB <- create.Learner("SL.custom_glm", tune = list(glm.model = list(band_form)), detailed_names = FALSE, name_prefix = "glmB")
SL.glm_all <- function(Y, X, newX, family, obsWeights, ...) {
  df <- data.frame(Y = Y, X, check.names = FALSE)
  fit <- glm(Y ~ ., data = df, family = family, weights = obsWeights)
  out <- list(object = fit); class(out) <- "SL.glm_all"
  list(pred = as.numeric(predict(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.glm_all <- function(object, newdata, ...) as.numeric(predict(object$object, newdata = newdata, type = "response"))
LEAN <- c(gamL$names, glmB$names, "SL.glm_all", "SL.glmnet", "SL.randomForest", "SL.mean")

daume_aug <- function(X) {
  Xo <- X[, setdiff(colnames(X), "is_optical"), drop = FALSE] * X$is_optical
  colnames(Xo) <- paste0(colnames(Xo), "_opt")
  cbind(X, Xo)
}

N_REP <- 3L; K <- 5L; N_DRAW <- 6L
set.seed(31415)
fold_sets <- lapply(seq_len(N_REP), function(r) caret::createFolds(Y, k = K, returnTrain = FALSE))

# task grid: (rep, variant, draw); anchor has draw = 0
grid <- rbind(
  expand.grid(rep = seq_len(N_REP), variant = "anchor", draw = 0),
  expand.grid(rep = seq_len(N_REP), variant = "param",  draw = seq_len(N_DRAW)),
  expand.grid(rep = seq_len(N_REP), variant = "full",   draw = seq_len(N_DRAW))
)
tasks <- split(grid, seq_len(nrow(grid)))

run_task <- function(t) {
  f <- if (t$variant == "anchor") f0 else make_frame(t$draw, full = t$variant == "full",
                                                     seed = 500 + t$draw)
  X <- daume_aug(f[, VARS])
  folds <- fold_sets[[t$rep]]
  pred <- rep(NA_real_, N)
  set.seed(210000 + 1000 * t$rep + t$draw + 100 * (t$variant == "full"))
  for (i in seq_along(folds)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    fit <- try(suppressWarnings(SuperLearner(
      Y = Y[tr], X = X[tr, ], newX = X[te, ], family = gaussian(),
      SL.library = LEAN, cvControl = list(V = 5), verbose = FALSE)), silent = TRUE)
    if (inherits(fit, "try-error")) { cat("ERR", as.character(t$variant), t$rep, t$draw, i, "\n"); next }
    pred[te] <- as.numeric(fit$SL.predict)
  }
  data.frame(rep = t$rep, variant = as.character(t$variant), draw = t$draw,
             row = seq_len(N), is_opt = is_opt, y = Y, pred = pred)
}

cat("Running", length(tasks), "tasks\n")
res <- parallel::mclapply(tasks, run_task, mc.cores = 5L, mc.preschedule = FALSE)
res <- do.call(rbind, res)
write.csv(res, paste0(out_prefix, "_preds.csv"), row.names = FALSE)

cat("\n== round 18 projection-uncertainty bagging ==\n")
for (v in c("anchor", "param", "full")) {
  rr <- sapply(seq_len(N_REP), function(rp) {
    d <- res[res$variant == v & res$rep == rp, ]
    pa <- tapply(d$pred, d$row, mean, na.rm = TRUE)   # average over draws
    yy <- tapply(d$y, d$row, mean); oo <- tapply(d$is_opt, d$row, mean) == 1
    c(all = cor(yy, pa), op = cor(yy[oo], pa[oo]))
  })
  cat(sprintf("%-7s r_all=%.4f  r_opt=%.4f\n", v, mean(rr["all", ]), mean(rr["op", ])))
}
