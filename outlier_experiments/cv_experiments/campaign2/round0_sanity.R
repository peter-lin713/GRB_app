# Round 0: validate harness reproduces baseline under screening settings + permutation null.
source("../rounds/harness_prep.R")
t0 <- Sys.time()
cfgs <- list(base = default_cfg())
res <- run_set(cfgs)

# permutation null: shuffle Y, rerun default config
cat("\n=== permutation null (shuffled Y) ===\n")
Y_orig <- Y
set.seed(99); Y <<- Y_orig[sample(length(Y_orig))]
res_perm <- run_config(default_cfg(), train_all, frames, "perm")
Y <<- Y_orig
cat(sprintf("\nBASE r_pool=%.4f  PERM r_pool=%.4f\n", res$base$r_pool, res_perm$r_pool))
cat("elapsed:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
cat("ROUND0_DONE\n")
