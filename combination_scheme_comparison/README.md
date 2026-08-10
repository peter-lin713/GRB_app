# Combination-scheme comparison

Everything from the multivariate/emcee/hard-zero-domain task batch, compiled in one place. Extended with a second batch: a multivariate Theil-Sen calibration, a no-plateau (6-var) ablation, and a corrected hard-zero-domain v2.

## Start here

`plots/master_scheme_comparison_table.png` — the headline comparison across all schemes tested, all using the "Without 2σ outliers" subset with corrected MICE (maxit=20). Every row states its **Formula** provenance (NEW own search vs. Reused standing) and whether it has **Opt cols** (the native-optical `_opt` columns added alongside the calibrated core columns) directly in the table, not just in prose.

## Folder guide

- `formulas/` — every winning formula found across both batches, with a README explaining provenance (which search found it, and whether it beat or lost to reusing the standing formula).
- `datasets/` — the new combination datasets built in the first batch:
  - `superlearner_training_multivariate_emcee.csv` — emcee-fit (MCMC) multivariate calibration, core columns only
  - `superlearner_training_hardzero_domain.csv` — hard-zero domain-separated v1 (Alpha_opt/Beta_opt/log10Fa_opt/log10Ta_opt columns, no regression fit at all)
  - (second-batch datasets live in `Data/` at the repo root: `superlearner_training_multivariate_theilsen.csv`, `superlearner_training_hardzero_domain_v2.csv` -- not copied in here to avoid duplicating large files)
- `scripts/` — modified copies of `superlearner.R` needed for datasets with the extra `_opt` columns (the stock script is blind to them and silently ignores them otherwise -- see "bugs found" below):
  - `superlearner_hardzero.R` — used for the hard-zero-domain dataset (v1 and v2)
  - `superlearner_optadd.R` — used for the "opt cols added" (Dataset B) run
  - `superlearner_optadd_theilsenv2.R` — a further-patched copy of `superlearner_optadd.R`: its internal errs MICE call uses `pmm` instead of `midastouch` for one column block (see bug #5 below); used for the Theil-Sen-multivariate + opt-cols-added run
- `hardzero_domain_run/` — full SuperLearner Plot_Output for the hard-zero-domain **v1** dataset (own fresh formula; shared GRBs' opt cols were hard-zeroed, discarding real optical data -- see "hard-zero-domain v2 fix" below)
- `hardzero_domain_v2_run/` — full SuperLearner Plot_Output for the hard-zero-domain **v2** dataset (own fresh formula, opt cols added; corrected shared-GRB handling)
- `multivariate_theilsen_formula_run/` — multivariate Theil-Sen calibration (all 4 optical params -> each X-ray-scale plateau param via `sklearn.TheilSenRegressor`), own fresh formula, no opt cols
- `theilsen_multivariate_optadded_reused_formula_run/` — same multivariate Theil-Sen calibration, with the v2-convention `_opt` columns added back in, reusing the standing formula instead of searching a fresh one (uses `scripts/superlearner_optadd_theilsenv2.R` -- see below)
- `no_plateau_formula_run/` — 6-variable ablation dataset (Alpha/Beta/log10Fa/log10Ta dropped entirely), own fresh formula, no opt cols
- `multivariate_emcee_formula_run/`
  - `dataset_A/` — multivariate-emcee core-only dataset, own fresh formula
  - `dataset_B_fixed/` — Dataset A + the 4 `_opt` columns added alongside (correct result)
  - `dataset_B_BUGGY_optcols_unused/` — kept for transparency: the first attempt at Dataset B, which silently didn't use the opt columns at all (identical numbers to Dataset A) because it used the stock `superlearner.R`. Superseded by `dataset_B_fixed/`.
- `results/` — plots from the earlier reused-formula runs (emcee_multivariate, multivariate_ridge, theilsen with/without is_optical, theilsen own formula)
- `theilsen_no_optid/` — full working folder for the Theil-Sen no-`is_optical` test
- `zips/` — the two fixed formula-generation zips (`no_plateau_formula_generation.zip`, `multivariate_formula_generation.zip`)

## Hard-zero-domain v2 fix

v1's hard-zero-domain dataset discarded real data: 78 GRBs have **both** a real X-ray-scale and a real optical-native measurement (the "shared" GRBs, also the calibration anchor sample for every other scheme). v1 followed the base dataset's X-ray-preferred `is_optical` convention and hard-zeroed these GRBs' `_opt` columns even though real optical values existed for them. v2 leaves those 78 GRBs' `_opt` columns as `NaN` instead, so MICE imputes them like any other missing value, while genuinely-X-ray-only GRBs (145 of them, no optical data at all) still get a structural zero. Core (`Alpha`/`Beta`/`log10Fa`/`log10Ta`) columns are unaffected -- still hard-zeroed only for the 73 optical-only GRBs, as in v1.

## Results summary (Without 2σ outliers, r on the z scale unless noted)

| Scheme | Formula | Opt cols? | N | r(z) | Note |
|---|---|---|---|---|---|
| Single-variable emcee (standing, reference) | Standing | No | 205 | 0.712 | current paper winner |
| Multivariate ridge | Reused | No | 207 | 0.678 | ridge shrinks Fa/Ta cross-terms, hurts a bit |
| Multivariate emcee (MCMC) | Reused | No | 202 | 0.714 | ties the reference almost exactly |
| Theil-Sen (single-var), with `is_optical` | Reused | No | 206 | 0.707 | |
| Theil-Sen (single-var), without `is_optical` | Reused | No | 206 | 0.679 | domain cols cost ~0.03 |
| Theil-Sen (single-var), own formula | NEW | No | 206 | 0.681 | own formula slightly worse than reused |
| Hard-zero domain **v1** | NEW | Yes | 208 | 0.683 | competitive; no calibration model at all |
| Multivariate emcee, Dataset A | NEW | No | 206 | 0.672 | own formula worse than reused (same pattern as Theil-Sen) |
| Multivariate emcee, Dataset B | Reused | Yes | 205 | 0.705 | adding native-optical cols alongside the calibrated ones helped |
| **Theil-Sen MULTIVARIATE** (4 opt params -> each xray param) | NEW | No | 203 | 0.692 | own formula slightly worse than reused standing (0.707) -- same recurring pattern |
| **No-plateau (6-var)**, drops Alpha/Beta/Fa/Ta entirely | NEW | No | 204 | 0.644 | ablation: dropping the plateau params entirely costs the most of any scheme tested |
| **Hard-zero domain v2** (shared GRBs' opt cols MICE-imputed, not hard-zeroed) | NEW | Yes | 207 | 0.691 | v1 fix -- essentially ties v1 (0.683), confirming the shared-GRB correction didn't change the picture much |
| **Theil-Sen MULTIVARIATE + opt cols added** (v2 convention) | Reused | Yes | 205 | **0.708** | best reused-formula result besides the two ~0.71 references; beats its own fresh formula (0.692) and edges out multivariate-emcee Dataset B (0.705) |

## Bugs found and fixed this batch (worth knowing before reusing any of this)

1. **Formula-installation gap**: the first version of the orchestration chain never installed each dataset's own freshly-found formula into `Best_formula_GAM.txt` before running its SuperLearner step -- it silently reused whatever formula was already installed. Fixed by adding explicit `install_formula()` calls before every SuperLearner invocation, and by re-running the affected steps (Theil-Sen's "own formula" result was corrected as a result).
2. **M-estimator singularity**: the hard-zero-domain dataset's top 7 ranked formulas all hit `'x' is singular` in `rlm()`, because any formula that interacts two `_opt` columns together is exactly collinear (both are zero for the same 223 X-ray-native rows). Fixed by falling back to the next-ranked formula that keeps all `_opt` terms as plain additive terms -- rank 8 (2 votes) fit cleanly in v1. This fallback was built proactively into the v2 formula-generation script itself (rather than patched on after a crash); it triggered again at v2's full scale and correctly fell back to rank 2, confirming the fix generalizes.
3. **Stock `superlearner.R` is blind to new columns**: its own `features_for_mice_preds` construction hardcodes the standard 10 features. Any dataset with extra columns (hard-zero's `_opt` columns, Dataset B's added `_opt` columns) needs a modified copy that explicitly includes them, or the extra columns are silently unused (as happened in the first, buggy Dataset B run -- caught because its numbers exactly matched Dataset A's).
4. **`$(date)` clobbers `$?`**: every orchestration script's `echo "... (exit $?)"` pattern was reporting `date`'s exit status, not the preceding command's -- harmless until the hard-zero-domain formula search actually crashed and still printed "(exit 0)". Worth fixing in any future chain scripts (capture `$?` into a variable immediately after the command, before any other command runs).
5. **`midastouch` singular on the `_opt` error columns under the v2 convention**: once a dataset's `_opt` columns follow the v2 convention (NaN for shared GRBs, not hard-zeroed), the four `*Err_opt` columns are identically 0 for every X-ray-only row (zero variance in that block) while genuinely NaN for the shared rows -- combined with `T90Err`'s separate NaN block (missing only for optical-only rows), some internal `midastouch` donor-regression subproblem is rank-deficient no matter the seed (tried 16 seeds, all failed identically). Confirmed this is v2-specific: Dataset B and hard-zero v1 never hit it because their `_opt` error columns were hard-zeroed for shared GRBs too (never actually missing), and hard-zero **v2**'s own formula-generation script didn't hit it either since its errs MICE call has a different column mix. Fixed by switching just that one `mice()` call to `method = "pmm"` instead of `"midastouch"` -- confirmed harmless since none of the `*Err` columns feed the M-estimator formula or the final outlier cut (`final_outliers_removed.csv` is written from the original un-imputed `raw_xray_data`, same as every other formula-generation script).
