# GAM 8-variable formula search — OT baseline dataset

Self-contained copy of `GAM_analysis_8variables.R`, adapted to run locally
against the OT baseline combined dataset instead of the original
`combined_data_with_redshift.csv`.

## What's here
- `GAM_analysis_8variables.R` — the script, with 3 changes from the original (see below)
- `Data/superlearner_training_ot_v3_errcut_relative.csv` — the OT baseline dataset (296 GRBs, error-cut), copied from `Data/` in the main repo

## Changes made to the original script
1. **Input file**: points to `Data/superlearner_training_ot_v3_errcut_relative.csv` instead of `combined_data_with_redshift.csv`.
2. **T90 filter**: the original did `raw_xray_data$T90 > 2` then derived `log10T90` from it. The OT dataset only ships `log10T90` (no raw linear T90 column), so the filter is now `raw_xray_data$log10T90 > 0.301` directly (same threshold, `T90 > 2s`, same convention `superlearner.R` uses elsewhere).
3. **`MyLaptop = TRUE`** (was `FALSE`): the original defaulted to spinning up a 71-worker **MPI** cluster (`makeCluster(71, type="MPI")`), clearly written for a supercomputer run. This machine has no MPI runtime installed at all (`mpirun`/`orted` not found) — with `MyLaptop=FALSE` the script would crash immediately at the cluster setup step. `TRUE` routes to a plain local cluster (`makeCluster(detectCores()-2)`), which works here.

Nothing else was changed — column names in the OT dataset (`log10Fa`, `log10Ta`, `Alpha`, `Beta`, `Gamma`, `log10Fluence`, `PhotonIndex`, `log10NH`, `log10PeakFlux`, `Redshift_crosscheck`, and the `*Err` columns) already match exactly what the script expects.

## Verified working (smoke-tested)
Ran the data-load → filter → MICE → LASSO portion standalone: 296 → 285 GRBs survive the T90 cut, MICE imputation completes, and LASSO ranks all 10 variables cleanly. The script keeps the top 8 by `|coefficient|`:
`log10PeakFlux, log10Ta, PhotonIndex, log10NH, log10Fa, Alpha, log10T90, Beta` (drops `Gamma` and `log10Fluence`, the two weakest).

## Before you run the full thing: this is a large compute job
The exhaustive formula search generates candidate GAM formulas over **every subset** of the 8 variables (quadratic/interaction, linear, and spline-smoothed variants combined), then for **each** candidate formula runs 50 repetitions of fresh 10-fold CV (`InnerLoop = 50`). Even parallelized across `detectCores()-2` local workers, this is realistically a multi-hour-to-multi-day job on a laptop — it was written for a cluster. Consider:
- reducing `InnerLoop` (currently 50) for a faster/rougher pass, or
- restricting `all_formula_list` to a subset (e.g. only formulas up to a certain size) before the big `foreach` loop, or
- just running it and letting it go in the background — nothing about it is incorrect, it's just large.

## To run
```bash
cd GAM_analysis_8var_ot_baseline
Rscript GAM_analysis_8variables.R
```
All outputs land in this same folder (relative paths, unchanged from the original): `grb_xray_imputed.csv`, `O1_formula_list`, `O2_formula_list`, `SmoothedO1_formula_list`, `Gam_train_set.csv`, `Gam_validation_set.csv`, `SuperLearner_complete_data.rds` (every candidate formula's CV + validation predictions), `Formulas_used.rds`.
