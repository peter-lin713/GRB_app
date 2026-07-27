## OT v3, native-optical-columns + 9-variable experiment — supercomputer handoff

**This is a test, confined to this folder.** Nothing in the main pipeline (`emcee_v2_cleaned_formula_generation/`, `ot_v3_cleaned_formula_generation/`) was touched. Mirrors `../emcee_v2_native_opt_formula_generation/`, swapped onto the OT-fusion (optimal-transport) catalogue instead of the emcee-calibrated one.

**What this tests:** two changes from the standard corrected pipeline at once:
1. The 73 optically-sourced GRBs get 4 **extra** candidate columns --
   `log10Fa_native`, `log10Ta_native`, `Alpha_native`, `Beta_native` (+ their error columns) --
   holding their real native optical-band measurements (from the optical catalog's `*_opt`
   fields). These sit **alongside**, not in place of, the standard OT-projected
   Alpha/Beta/log10Fa/log10Ta -- the originals are untouched for all 296 GRBs. The 223
   X-ray-native GRBs get NA on the 4 new columns, imputed by MICE like any other missing
   predictor. This lets LASSO's own ranking decide whether the native version carries additional
   signal, rather than presupposing an answer by overwriting.
   Built by `build_native_opt_dataset_ot_v3.R` (the OT analogue of
   `combine_opt_xray_emcee/build_native_opt_dataset.R`), which reads
   `Data/superlearner_training_ot_v3_recovered.csv` + the optical catalog and writes
   `Data/superlearner_training_ot_v3_native_opt.csv` (296 GRBs x 28 cols, all 73
   `is_optical == 1` rows matched). A local copy of that output CSV is included in this
   folder so the two analysis scripts below are self-contained.
2. LASSO keeps the top **9** variables (from the pool of 14: the standard 10 + the 4 native
   columns) instead of the paper's top-7, same as the emcee experiment.

**Cost:** 9 vars -> 18 O2-columns -> 2^18-1 = 262,143 formulas (16x the standard 16,383) x 100
splits x 10-fold CV. Not attempted locally (~67 hours at laptop-scale parallelism) -- both
scripts are configured `MyLaptop = F` (72-1 MPI workers) for cluster submission.

**Two independent searches, matching the paper's Sec 4.2.3** ("The sets of best formulas for GAM
and GLM are determined independently"):
- `GAM_analysis_8variables.R` -- fits `mgcv::gam()`, writes to `checkpoints/`,
  `formula_win_frequency.csv`, `Formula_for_outlier.txt`, `final_outliers_removed.csv`.
- `GLM_analysis_8variables.R` -- identical data/LASSO/formula-generation, fits `stats::glm()`
  instead, writes to `checkpoints_GLM/`, `formula_win_frequency_GLM.csv`,
  `Formula_for_outlier_GLM.txt`, `final_outliers_removed_GLM.csv`.

**Non-negotiable outlier detection** is delegated to `../../generate_nonnegotiable_outliers.R`
(repo root) -- the single source of truth, not duplicated per folder. Both scripts source it and
call `update_nonnegotiable_outliers()` at startup, so `candidate_outliers_ot_v3_native_opt.csv`
is recomputed fresh from whatever data is currently loaded every run (never hand-edit it), and
`confirmed_outliers_to_drop.txt` gets the freshly-detected automatic candidates unioned in on top
of anything already there, so manual additions from 4D-plot visual review persist across reruns.
The err>value check is restricted to the standard Alpha/Beta/log10Fa/log10Ta only (the native
columns aren't subject to this cut; they're supplementary candidates, not core physical
measurements assumed always present) -- expect the same 4 GRBs (090516A, 100621A, 131117A,
230328B) that both the emcee and OT combinations flag independently on the standard columns,
though this run will recompute and confirm that rather than assume it.

**Caution:** both scripts regenerate `O2_formula_list`/`O1_formula_list`/`SmoothedO1_formula_list`/
`grb_xray_imputed.csv` independently (same content expected, since both start from the same data,
but LASSO's 100-rep ranking is unseeded). Don't launch both truly concurrently from this same
directory -- run sequentially, or copy the folder so each has its own working directory, to avoid
a write race on those shared intermediate files. (The outlier-detection files are deterministic,
so no race risk there even if run concurrently.)
