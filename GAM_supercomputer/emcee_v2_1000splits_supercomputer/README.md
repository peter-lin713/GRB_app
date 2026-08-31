# emcee v2, 1000-split formula search — supercomputer handoff

**What this is:** the exact same formula-generation search as the main
pipeline's flagship result (`GAM_supercomputer/emcee_v2_cleaned_formula_
generation/`, the one behind the reported r=0.720 result), but with **1000
random 80/20 splits instead of 100** — nothing else changed. This is
deliberately a "more splits" experiment, not a "more variables" or
"higher-order interactions" one; see the reasoning below.

**Why more splits, not more variables or higher interaction order:** with
this pipeline's discussed options for spending extra compute (more LASSO
variables at pairwise order, more splits, or higher polynomial degree
per formula), higher-degree interactions (3rd/4th/5th order) were ruled
out — with ~260-280 GRBs, formulas with many interaction terms overfit
badly and have no physical motivation (GRB plateau relations are empirical,
mostly pairwise). More splits, by contrast, directly strengthens the
*existing* selection procedure's statistical reliability (the 100-split
vote tally is already an implicit multiple-comparisons correction across
16,383 candidates; 1000 splits makes that correction far more robust)
without expanding the model space or increasing per-formula overfit risk.

**Confirmed identical formula space:** this package reuses the exact
`O1_formula_list` / `O2_formula_list` / `SmoothedO1_formula_list` /
`candidate_outliers_emcee_v2.csv` / `confirmed_outliers_to_drop.txt` files
already committed from the original 100-split run, rather than letting the
LASSO ranking step (unseeded, 100 reps of `cv.glmnet`) possibly regenerate
a different top-7 variable set. `GAM_analysis_8variables.R` is patched to
skip rewriting those files if they already exist (search for "1000-split
run" in the script) — so this run and the original 100-split run are
searching the identical 16,383-formula (O2) + 127-formula (O1) space of
16,510 total candidates, differing only in how many splits vote.

**Verified locally before packaging:** ran a 3-split smoke test
(`../emcee_v2_1000splits_SMOKETEST/`) with this exact setup (just
`MyLaptop=T`, `N_SPLITS=3` instead of `F`/`1000`). It reproduced the
original pipeline's numbers exactly — `Formulas: 16383 | GRBs: 281`,
`M-estimator outlier cut: 15 of 281 GRBs removed... 266 GRBs`, and 2 of
the 3 smoke-test splits picked the exact same formula as the real
100-split winner (`(log10FaSqr+log10T90+log10Fa+log10PeakFlux)^2 + ...`).
Confirms the reused candidate space and non-negotiable-cut logic are
wired correctly before spending supercomputer time on it.

## Cost estimate

Locally (12-core laptop, `MyLaptop=T`, `detectCores()-1` = 11 parallel
workers): **2.5 minutes/split**. At `N_SPLITS=1000` that would be
~42 hours on this laptop alone — the reason this needs a cluster.

With `MyLaptop=F` here configured for **71 MPI workers** (`slaves <- 72-1`,
matching this project's established convention for prior supercomputer
submissions), scaling by worker-count ratio (71/11 ≈ 6.5x) gives a rough
**optimistic estimate of ~6-7 hours** for all 1000 splits, assuming
near-linear parallel scaling. Real-world MPI communication overhead
(distributing 16,383 formula-evaluation tasks per split across ranks) and
this cluster's actual per-core speed (unknown, could differ from the
laptop's cores in either direction) will change this — treat 6-7 hours as
optimistic and budget for a submission wall-time closer to **12-20 hours**
to be safe. Resume support is built in (see below), so a wall-time cap that
turns out too short just means resubmitting to finish the remaining splits,
not losing progress.

## Requirements on the cluster

R packages (from the script's own `library()`/`require()` calls):
`doParallel`, `mice`, `ggplot2`, `lattice`, `stringr`, `dplyr`, `MASS`,
`randomForest`, `earth`, `glmnet`, `SuperLearner`, `mgcv`, `xgboost`,
`gbm`, `caret`, `biglasso`, `arm`.

**Not an explicit `library()` call but still required:** `parallel::
makeCluster(..., type = "MPI")` (used when `MyLaptop = F`, which this
package is set to) needs the **Rmpi** package installed and an MPI
implementation (OpenMPI or MPICH) available on the cluster — R loads Rmpi
internally for the MPI cluster type, so it won't show up as a `library()`
call anywhere in the script but the job will fail at `makeCluster()`
without it.

**No SLURM/PBS batch script is included.** Checked this project's full git
history — no supercomputer job-submission template was ever committed here
(prior supercomputer-bound experiments in this repo, e.g. the emcee_v2
native-opt/top-9-variable test, were prepared the same way — R script
configured `MyLaptop=F` and handed off — but the actual `sbatch`/`qsub`
script was written external to this repo, matching whatever the target
cluster's scheduler requires). You'll need to write the actual job script
for your specific cluster around:

```
Rscript GAM_analysis_8variables.R
```

as the payload, requesting 72 MPI ranks/tasks (or adjust `slaves <- 72 - 1`
on line 740 of the script to match however many ranks you actually
request) and enough wall-time per the estimate above.

## Resume support

The script checkpoints each split to `checkpoints/split_%03d.rds` and skips
any split already checkpointed on rerun (`done_ids` logic near the top of
the split loop). If the job gets killed partway through or the wall-time
cap is hit before all 1000 splits finish, just resubmit the same job — it
picks up where it left off rather than restarting from split 1.

## Output

Same as every other formula-generation run in this pipeline:
`checkpoints/best_formulas_per_split.csv` (per-split winners),
`formula_win_frequency.csv` (the full vote tally across all 1000 splits),
`final_outliers_removed.csv` (the 266-GRB post-M-estimator-cut data, ready
for the SuperLearner stage). Once complete, this can be plugged into
`superlearner_beta2.R` exactly like the 100-split winner was, for a direct
100-split-vs-1000-split comparison on the identical formula space and data.
