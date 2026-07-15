# Improving the CV correlation for GRB redshift prediction — experiment report

*July 14, 2026. Scripts and per-repeat results in this folder; run scripts from the repo root.*

> **FINAL (rounds 4–25):** keeping both bands, the campaign finished at
> **CV r = 0.684 log / 0.712 linear** (n=227, pipeline seed-42 convention, retrained
> 2σ outlier removal; 0.661/0.688 before outlier removal on n=233) vs 0.51 at baseline
> and 0.646 published on the cleaner pure-X-ray sample. Untouched 58-GRB holdout:
> 0.415 log (old frame scored 0.324 on the same split; n=58 ⇒ SE≈0.11).
> Full history in `LEDGER.md`; final predictions and plots in `results/`
> (`final_cv_correlation.png`, `campaign_trajectory.png`). The winning configuration
> is integrated into the root `superlearner.R` on this branch.

## Setup

- **Data:** the emcee-combined training frame (`Data/superlearner_training_emcee_errcut_relative.csv`)
  after MICE and the 5% M-estimator hard cut — i.e. the
  `outlier_experiments/runs/combined_5pct/OutputFiles/grb_xray_m_est.csv` frame:
  **n = 291**, of which **64 are optical-projected GRBs** carrying 6/10 MICE-imputed
  features (they only have logFa, logTa, Alpha, Beta from the emcee projection).
- **Protocol:** 5 repeats × 5-fold CV with **identical fold assignments for every
  method** (seed 31415), so all comparisons are paired and split-luck cancels.
  Metric: Pearson r between out-of-fold predictions and observed log10(z+1).
- **Baseline:** generic library (glmnet, ranger, earth, xgboost, mean) on the
  10 linear + 10 squared predictors: **r = 0.510**.
- **Permutation null:** shuffled target gives r = −0.03 ± 0.07 (6 permutations), so
  the signal is real (~7σ) and method differences below ~±0.02 are noise.

## Results ladder (all on the same paired folds)

| Config | r (log10(z+1)) | paired Δ vs baseline | verdict |
|---|---|---|---|
| **Formula GAM/GLM + RF, top-7 features** | **0.547** | **+0.038 (t=4.9, 5/5 reps +)** | **best** |
| Formula GAM/GLM + RF, + is_optical flag | 0.547 | +0.037 (t=6.7, 5/5 reps +) | tied |
| Formula GAM/GLM + RF (all 20 vars) | 0.539 | +0.029 (t=3.6, 5/5 reps +) | real gain |
| top-7 LASSO features, generic library | 0.531 | +0.021 (t=2.9) | modest |
| CC_LS metalearner (generic lib) | 0.518 | +0.009 | noise-level |
| rich tuned library (SVM/ridge/tuned RF+XGB) | 0.516 | +0.006 | noise |
| physical combos (logFa+logTa, logFluence−logT90) | 0.516 | +0.006 | noise |
| log10(z) target instead of log10(z+1) | 0.515 | +0.005 | noise |
| error-based obsWeights | 0.508 | −0.002 | no |
| measurement errors as features | 0.498 | −0.012 | no |
| linear-only features (drop squares) | 0.502 | −0.008 | squares help |
| 2× upsampling of z>2.5 | 0.464 | −0.046 | hurts |
| CC_LS + formula lib | 0.524 | — | worse than NNLS |

**Headline: r 0.510 → 0.547** (log scale; linear-z r = 0.551, RMSE(log) = 0.149).
On the X-ray subset the best config reaches **r_xray = 0.582**.

## Key structural findings

1. **The formula learners are switched off in the outlier-experiment pipeline and
   that costs ~0.03–0.04 of r.** The sweep (`run_outlier_experiments.sh`) passes
   `use_formula_learners=FALSE`, so the tuned GAM/GLM formulas
   (`Best_formula_GAM.txt` / `Best_formula_GLM.txt`) never enter the library.
   Re-enabling them is the single biggest, replicable gain — consistent with the
   ApJS 271:22 finding that GAM+GLM carry essentially all SuperLearner weight
   (0.649 + 0.351) at this sample size, with trees contributing ~nothing.

2. **The 64 optical-projected rows slightly dilute the X-ray predictions.**
   Scored on the same X-ray rows with identical folds: train-on-all r_xray = 0.545 vs
   train-on-X-ray-only r_xray = 0.559 (generic lib); under the formula library the gap
   narrows (0.572 vs 0.580). The projected rows are themselves predicted much worse
   (r ≈ 0.41 vs 0.55–0.58). Since these rows have 6 of 10 features fabricated by MICE,
   any dataset-combination scheme (incl. the planned optimal transport) should aim to
   carry *real* information into those columns — e.g. project more optical features, or
   give projected rows their own noise model / obsWeights — rather than imputing them
   from the X-ray rows' joint distribution. A cheap robust fix that costs nothing:
   keep the rows but add the `is_optical` indicator (tied-best config above).

3. **The emcee frame is intrinsically harder than the paper's sample.** Even the pure
   X-ray subset (n=227) under the paper's 10-fold protocol gives r ≈ 0.557, well below
   the r = 0.646 reported on the 180-GRB pure X-ray sample (arXiv:2410.13985). The gap
   between "our pipeline" and "the paper" is mostly the dataset, not the method.

4. **The optimal-transport / piecewise bias correction does not survive leakage-free
   evaluation.** Cross-fitted (correction fit on half the OOF predictions, applied to
   the other half), it *lowers* r from ~0.51 to ~0.44. The published r = 0.89–0.93
   values fit the correction on the same predictions they evaluate; the honest
   apples-to-apples comparison to published work is their *pre-correction*
   r = 0.646–0.762. Recommend reporting raw CV r alongside any corrected number.

5. **There is an information ceiling near r ≈ 0.55 for this frame.** Everything tested
   lands in 0.46–0.55. The literature's uncorrected results on richer/cleaner samples
   have plateaued at 0.65–0.76 for 15 years (Ukwatta 2016 → Rácz → Dainotti/Narendra),
   with physics reasons: D_L(z) flattens beyond z ≈ 2, GRB luminosity-function scatter
   is several dex, and Malmquist bias thins the high-z training coverage.

## Recommendations (in order)

1. Re-enable `use_formula_learners=TRUE` (arg 10) in the sweep, with the top-7 design
   — or simply adopt `cv_round3.R`'s `pf_top7`/`pf_opt` config. Expected: r ≈ 0.55
   overall, ≈ 0.58 on X-ray rows.
2. Add the `is_optical` indicator column to the design matrix (free, robust).
3. For the dataset-combination effort: judge any new scheme (OT included) by the paired
   protocol here — same folds, scored on X-ray rows, against the train-X-ray-only
   control (r_xray = 0.580). A combination method only earns its place if it beats
   that. The current emcee projection does not.
4. Stop tuning libraries/metalearners/weights on this frame — all within noise.
   The leverage is in data quality: more real (non-imputed) features for the optical
   rows, or growing the X-ray sample.
5. Report raw CV r; treat bias-corrected r as a calibration diagnostic only.

## Part 2 (rounds 4–7): keeping both bands and fixing the optical rows

The project requires using X-ray **and** optical GRBs together, so instead of
dropping/downweighting the optical rows, rounds 4–7 attack their actual defect:
6 of 10 features were MICE-fabricated.

**Round 4 — recover real features (`cv_round4.R`).** `Data/optical_data.csv`
carries T90, Fluence, PhotonIndex, NH and PeakFlux for all 69 optical GRBs —
prompt/Swift quantities the emcee combination discarded. Unit conversions were
validated on the 88 overlap GRBs (r = 1.000 exact matches): `log10Fluence =
log10(Fluence) − 7` (10⁻⁷ erg/cm²), `log10NH = log10(NH) + 21` (10²¹ cm⁻²),
`log10PeakFlux = log10(PeakFlux)`, `log10T90 = log10(T90)`. Only Gamma (and ~30%
scattered gaps) still need imputation. Result, same paired folds:

| Frame / config | r_all | r_xray | r_opt |
|---|---|---|---|
| old frame (round-3 best) | 0.547 | 0.582 | 0.435 |
| recovered frame | 0.566 | 0.599 | 0.435 |
| recovered frame + is_optical | **0.572** (+0.025, t=5.0) | 0.599 | 0.416 |
| old frame, optical ½-weight (control) | 0.553 | 0.591 | 0.431 |

Recovery beats downweighting, and r_xray improves too — the corrected rows stop
contaminating the joint feature distribution.

**Round 5 — MICE-completion averaging (`cv_round5.R`).** Averaging predictions
over 3 MICE completions instead of using one: r = 0.575 → **0.579**.

**Round 6 — expanded library (`cv_round6.R`).** Found that `Best_formula_GAM.txt`
has only 3 lines (2 distinct) and `Best_formula_GLM.txt` only 3 — the "6 GAM/4 GLM"
indexing silently created broken learners, and the accidental 4th GLM (a plain
all-features linear GLM) carries the largest ensemble weight (0.34). A properly
deduplicated library with new smooth mgcv GAMs (incl. select=TRUE), a Gaussian
process, and a band-interaction GLM (`is_optical × plateau features`, weight 0.26)
gave r = 0.570 — no net gain; the frame is at its ceiling near 0.57–0.58 under
5-fold CV. Regenerating the formula files on the corrected frame is the open lever.

**Round 7 — corrected pipeline headline (`cv_round7.R`).** Full rebuild: recover
features on the raw frame → clean → MICE(m=20) → fresh 5% M-estimator cut computed
on corrected features → 10-fold CV (paper protocol) × 3 repeats × 3-completion
averaging, expanded library. (Results in `results/round7_preds.csv`.)

## Reproducing

```bash
cd GRB-Web-App
Rscript outlier_experiments/cv_experiments/cv_experiments.R 1 round1.csv  # method screen
Rscript outlier_experiments/cv_experiments/cv_round1b.R round1b.csv      # paper-inspired
Rscript outlier_experiments/cv_experiments/cv_round2.R round2            # dilution + BC
Rscript outlier_experiments/cv_experiments/cv_round3.R round3            # winners combined
Rscript outlier_experiments/cv_experiments/cv_permctrl.R                 # permutation null
Rscript outlier_experiments/cv_experiments/cv_xray10f.R                  # xray-only 10-fold
```

Per-repeat metrics and per-row OOF predictions are in `results/`.
