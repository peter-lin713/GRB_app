# Round 2: leave-one-out ablations. What actually carries weight?
source("../rounds/harness_prep.R")
base <- default_cfg()
cfgs <- list(
  base        = base,
  no_gam      = modifyList(base, list(use_gam = FALSE)),
  no_band_glm = modifyList(base, list(use_band_glm = FALSE)),
  no_glm_all  = modifyList(base, list(use_glm_all = FALSE)),
  no_glmnet   = modifyList(base, list(use_glmnet = FALSE)),
  no_rf       = modifyList(base, list(use_rf = FALSE)),
  no_gamboost = modifyList(base, list(use_gamboost = FALSE)),
  no_iso      = modifyList(base, list(use_iso = FALSE)),
  no_infold   = modifyList(base, list(use_infold = FALSE)),
  meta_NNLS   = modifyList(base, list(meta = "NNLS"))
)
t0 <- Sys.time()
res <- run_set(cfgs)
cat("elapsed:", round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n")
cat("ROUND2_DONE\n")
