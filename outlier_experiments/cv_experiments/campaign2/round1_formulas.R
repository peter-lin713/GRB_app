# Round 1: formula regeneration on the new frame. The supplied GAM formula has no
# s() terms, so the "GAM" learner is actually parametric. Test: (a) genuine smooth
# GAM replacing it, (b) smooth GAM ADDED as extra learner, (c) data-driven stepwise
# interaction formula, (d) larger-k smooths.
source("../rounds/harness_prep.R")

# genuine smooth GAM extra learner (penalized, select=TRUE)
SL.gam_smooth <<- function(Y, X, newX, family, obsWeights, ...) {
  num <- intersect(linear_vars, colnames(X))
  fm <- as.formula(paste("Y ~", paste(sprintf("s(%s, k=5)", num), collapse="+"), "+ is_optical"))
  fit <- mgcv::gam(fm, data = X, family = family, weights = obsWeights, select = TRUE,
                   control = mgcv::gam.control(maxit = 50))
  out <- list(object = fit); class(out) <- "SL.gam_smooth"
  list(pred = as.numeric(mgcv::predict.gam(fit, newdata = newX, type = "response")), fit = out)
}
predict.SL.gam_smooth <<- function(object, newdata, ...)
  as.numeric(mgcv::predict.gam(object$object, newdata = newdata, type = "response"))

# smooth-GAM formula string to REPLACE the interaction gam (fit without select in wrapper)
smooth_rhs <- paste(sprintf("s(%s, k=5)", linear_vars), collapse = " + ")
smooth_form <- as.formula(paste("Response ~", smooth_rhs))

# data-driven stepwise interaction formula on the train set (frame 1), BIC
dtr <- data.frame(Y = Y[train_all], frames[[1]][train_all, linear_vars])
null_m <- lm(Y ~ ., data = dtr)
scope_up <- as.formula(paste("~ (", paste(linear_vars, collapse = "+"), ")^2"))
step_m <- MASS::stepAIC(null_m, scope = list(lower = ~1, upper = scope_up),
                        direction = "both", k = log(nrow(dtr)), trace = FALSE)
dd_rhs <- paste(deparse(formula(step_m)[[3]]), collapse = " ")
cat("data-driven RHS:", dd_rhs, "\n")
dd_form <- as.formula(paste("Response ~", dd_rhs))

base <- default_cfg()
cfgs <- list(
  base            = base,
  gam_smooth_repl = modifyList(base, list(gam_forms = list(smooth_form))),
  gam_smooth_add  = modifyList(base, list(extra_learners = "SL.gam_smooth")),
  gam_dd_repl     = modifyList(base, list(gam_forms = list(dd_form))),
  glmB_dd_repl    = modifyList(base, list(band_form = dd_form)),
  gam_dd_both     = modifyList(base, list(gam_forms = list(dd_form), band_form = dd_form))
)
t0 <- Sys.time()
res <- run_set(cfgs)
cat("elapsed:", round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n")
cat("ROUND1_DONE\n")
