# Theil-Sen (opt→x-ray, the 0.707 main result), 1000-split formula search — supercomputer handoff

**What this is:** the same "more splits, not more variables or higher
interaction order" experiment as `../emcee_v2_1000splits_supercomputer/`
(see that package's README for the full reasoning), applied to the *other*
flagship result in this pipeline: the single-variable Theil-Sen
optical→X-ray calibration behind the reported **r=0.707** main result
(`GAM_supercomputer/theilsen_cleaned_formula_generation/`, whose CV output
lives in `runs/theilsen_cleaned_linearz_5pct_mice20/`). Same idea, same
mechanics, different source dataset — **1000 random 80/20 splits instead
of 100**, nothing else changed.

**Confirmed identical formula space:** reuses the exact `O1_formula_list` /
`O2_formula_list` / `SmoothedO1_formula_list` / `candidate_outliers_
theilsen.csv` / `confirmed_outliers_to_drop.txt` already committed from the
original 100-split Theil-Sen run, rather than letting the unseeded LASSO
ranking step possibly regenerate a different top-7 variable set.
`GAM_analysis_8variables.R` is patched (same as the emcee_v2 package) to
skip rewriting those files if they already exist. Original winning formula
for reference (14/100 votes):
```
Response ~ (log10FaSqr + log10Fa + log10PeakFlux)^2 + log10NH +
    PhotonIndex + log10Ta + log10T90 + Alpha + log10NHSqr +
    log10PeakFluxSqr + PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr
```

**Verified locally before packaging:** ran a 3-split smoke test
(`../theilsen_1000splits_SMOKETEST/`) with this exact setup (`MyLaptop=T`,
`N_SPLITS=3`). Reproduced the original pipeline's numbers exactly --
`Formulas: 16383 | GRBs: 281` -- confirming the reused candidate space and
non-negotiable-cut logic are wired correctly before spending supercomputer
time on it (same 4 non-negotiable-cut GRBs as every other combination in
this pipeline: 090516A, 100621A, 131117A, 230328B).

## Cost estimate

Same order of magnitude as the emcee_v2 package -- identical formula count
(16,383 O2 + 127 O1 = 16,510), identical GRB count (281 pre-M-estimator),
identical per-formula fitting cost (10-fold `mgcv::gam()` CV). Expect
**roughly the same ~6-7 hour optimistic / 12-20 hour budgeted wall-time**
with 71 MPI workers (see the emcee_v2 package's README for the scaling
derivation from local per-split timing). Run the two packages' smoke tests
side by side if you want a same-cluster confirmation before committing
both to a full submission.

## Requirements on the cluster

Identical to the emcee_v2 package: same R package list (`doParallel`,
`mice`, `ggplot2`, `lattice`, `stringr`, `dplyr`, `MASS`, `randomForest`,
`earth`, `glmnet`, `SuperLearner`, `mgcv`, `xgboost`, `gbm`, `caret`,
`biglasso`, `arm`), plus **Rmpi** + an MPI implementation (OpenMPI/MPICH)
for `parallel::makeCluster(..., type="MPI")` when `MyLaptop=F` (set here).

**No SLURM/PBS batch script is included** -- same situation as the emcee_v2
package; write the job script for your cluster's scheduler around:

```
Rscript GAM_analysis_8variables.R
```

requesting 72 MPI ranks/tasks (or adjust `slaves <- 72 - 1` on line 717 of
the script to match however many ranks you actually request).

## Resume support

Same checkpoint-and-skip logic as every formula-generation script in this
pipeline: `checkpoints/split_%03d.rds` per split, automatically skipped on
rerun if already present. Resubmitting an interrupted job picks up where
it left off.

## Output

`checkpoints/best_formulas_per_split.csv`, `formula_win_frequency.csv`
(full 1000-split vote tally), `final_outliers_removed.csv` (266-GRB
post-M-estimator data, ready for `superlearner_beta2.R`) -- same as every
other formula-generation run in this pipeline, directly comparable against
the original 100-split Theil-Sen result once complete.
