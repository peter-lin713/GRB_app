# Reproducing the OT v3 result

Snapshot updated 2026-07-15 after removing a target-leakage bug (see "Leakage
found and fixed" below) -- the r=0.73/r=0.70 numbers quoted in earlier
versions of this file are superseded and were inflated by that leakage. The
script and data in this folder now reflect the leakage-free version. See
`runs/ot_v3_noleak_5pct/run.log` for the corrected MC metrics.

## Files in this folder
- `combine_optical_xray_ot_v3.py` — the OT-fusion script, exact version that produced this result
- `superlearner.R` — main pipeline script, exact version (includes the dex-error
  reconstruction fix at the input-normalization step, needed to ingest datasets
  that ship log-space errors instead of linear ones)
- `OnlyLGRBs_data_171_optical_corrected.csv`, `Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv` — upstream catalog inputs to the OT script
- `emcee_chains_v2.npz` — cached emcee MCMC calibration chains the OT script projects optical -> X-ray scale with
- `superlearner_training_ot_v3.csv`, `superlearner_training_ot_v3_errcut_relative.csv` — the OT script's output (296 GRBs after error cut); this is what was fed to superlearner.R
- `ot_v3_loo_validation.csv` — leave-one-out validation table (per-anchor true vs. imputed values for the 5 OT-imputed features)
- `quick_ridge_check.py` — fast ridge-regression sanity gate

## To reproduce from scratch
```bash
cd GRB_app
cp runs/ot_v3_5pct/repro/combine_optical_xray_ot_v3.py ot_fusion/
cp "runs/ot_v3_5pct/repro/OnlyLGRBs_data_171_optical_corrected.csv" Data/
cp "runs/ot_v3_5pct/repro/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv" Data/
cp runs/ot_v3_5pct/repro/emcee_chains_v2.npz Data/
cp runs/ot_v3_5pct/repro/superlearner.R .

python3 ot_fusion/combine_optical_xray_ot_v3.py   # -> Data/superlearner_training_ot_v3_errcut_relative.csv
python3 quick_ridge_check.py Data/superlearner_training_ot_v3_errcut_relative.csv   # sanity gate, expect r~0.43

mkdir -p runs/ot_v3_5pct_repro_test
Rscript superlearner.R Data/superlearner_training_ot_v3_errcut_relative.csv TRUE FALSE TRUE FALSE 0.65 10 runs/ot_v3_5pct_repro_test 0.05 FALSE hard \
  > runs/ot_v3_5pct_repro_test/run.log 2>&1
```

Note: `superlearner.R`'s CV folds/MC draws are seeded (`RNGkind("L'Ecuyer-CMRG")`
+ `set.seed(20240603)` before the 10-rep `mclapply`), and the 80/20 holdout
split is seeded (`set.seed(42)`), so results should reproduce closely modulo
R package version drift (this machine's R was reinstalled from scratch this
session, so package versions here may not exactly match a future re-run).

## Bug found and fixed after the initial snapshot
The first version of `combine_optical_xray_ot_v3.py` had `xray.rename()` map
the X-ray catalog's log-space error columns (`logFluenceErr`,
`log10PeakFluxErr`) straight to `FluenceErr`/`PeakFluxErr` with no unit
conversion -- so those columns held dex-scale values mislabeled as linear.
`superlearner.R`'s `to_dex_err()` then divided by `10^log10Fluence` a second
time, a double conversion that exploded for GRBs with very negative
log10Fluence. GRB 090927A hit this exactly: a normal ~0.03 dex error came out
as `log10FluenceErr=141458` after the double conversion, which alone was
enough to send the MC-propagated RMSE into the billions even though the
underlying model (per the raw-CV number) was fine the whole time.

Fixed by converting `xray['FluenceErr']`/`['PeakFluxErr']` (understood as the
dex values they actually are) to real linear errors before they enter the
pipeline: `linear_err = dex_err * 10**log10_value * ln(10)`. Verified fix:
GRB 090927A's dex error re-derives to a sane 0.065, MC metrics went from
r≈0/RMSE=Inf to r=0.70 (log)/0.66 (z), and the sample size (216 rows) is
unchanged -- no data was sacrificed to get there.

## Leakage found and fixed (supersedes the fix above)
`MATCH_COLS` originally included `'z'` (= `Redshift_crosscheck`, the actual
training target) alongside the four calibrated Dainotti parameters. Since
`z` got by far the largest anchor-informed weight (~3.6x vs ~0.3-0.4x for
the others), the OT match for each optical-only GRB was driven mostly by
"find X-ray donors with the closest redshift" -- so the imputed Gamma/
PhotonIndex/NH/Fluence/PeakFlux values were contaminated with knowledge of
that GRB's own true redshift, baked in before the SuperLearner train/test
split ever happens. Two problems: (1) it inflates every CV/MC metric
reported for the pre-fix version, and (2) it's not reproducible at real
inference time, since a GRB with genuinely unknown redshift can't be
OT-matched on redshift in the first place -- which is the entire point of
the model.

Fixed by dropping `'z'` from `MATCH_COLS`, leaving only
`['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']`. Re-ran the bandwidth sweep,
LOO validation, and full SuperLearner pipeline after the fix -- see
`runs/ot_v3_noleak_5pct/run.log` for the corrected (lower, honest) numbers.
Quick ridge check dropped from r=0.470 to r=0.434 (log scale) after removing
the leak -- real signal lost, as expected, but not degenerate.
