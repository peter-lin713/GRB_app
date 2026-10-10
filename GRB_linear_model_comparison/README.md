# Does the flagship protocol inflate a simple linear model too?

**Question.** The flagship GRB redshift estimator reports r ≈ 0.71 (linear z) via
a three-step evaluation: (1) an M-estimator "catastrophic-outlier" cut driven by
the winning GLM formula (drops the lowest-weight ~5%); (2) repeated 10-fold CV;
(3) a post-hoc "2σ cone" clip that removes CV points > 2σ from the 1:1 line
before recomputing r. Is that ~0.71 specific to the elaborate winning GLM, or
does a **plain 7-feature additive linear model** (no interactions, no squares)
reach nearly the same number under the *identical* protocol?

## Files
- `grb_imputed_281.csv` — the 281 long-GRB labeled feature matrix (pre-M-est-cut).
- `generalization_299_raw.csv` — the 299 GRBs with NO measured redshift.
- `simple_linear_flagship_check.py` — the numbers experiment. Run: `python3 simple_linear_flagship_check.py`
- `make_comparison_figure.py` — predicted-vs-observed side-by-side. Run: `python3 make_comparison_figure.py`
- `make_generalization_figure.py` — the no-redshift distribution comparison. Run: `python3 make_generalization_figure.py`
- `make_scatter_matrix.py` — scatter matrix with the M-est-removed bursts colored. Run: `python3 make_scatter_matrix.py`
- `make_luminosity_redshift.py` — plateau luminosity vs redshift with the flux
  limit curve, M-est-removed bursts colored. Run: `python3 make_luminosity_redshift.py`
- `make_luminosity_cut_glm.py` — winning GLM under a physical luminosity cut vs the
  statistical M-est cut (flagship protocol). Run: `python3 make_luminosity_cut_glm.py`
- `outlier_train_vs_eval.py` — do the 15 outliers hurt via training or evaluation?
  2x2 train/test decomposition + clean-trained prediction of the held-out 15.
  Run: `python3 outlier_train_vs_eval.py`
- `outlier_mahalanobis.py` — per-outlier Mahalanobis distributional test (each of
  the 15 vs the 266 used bursts, normalized for scale/correlation/dimension).
  Run: `python3 outlier_mahalanobis.py`
  (need numpy, pandas, patsy, statsmodels, matplotlib, scipy, seaborn.)

## Figures (for side-by-side comparison)
- `comparison_linear_vs_glm.png` — 2x2 predicted-vs-observed panels: columns =
  winning GLM vs plain linear model; top = flagship protocol (2σ-clipped points in
  red), bottom = honest CV. Same axes/style for direct comparison.
- `generalization_glm_vs_linear.png` — predicted-redshift **distribution** of the
  no-redshift GRBs (predict in log10(1+z), transform to z) vs the spectroscopic
  sample; winning GLM (left) vs simple linear (right). Both give the same
  compressed spike (KS D ≈ 0.25, p ≈ 1e-8) — the elaborate model adds nothing on
  generalization either. Both models use the same 7 features (no Beta/Fluence), so
  no missing-feature imputation is needed; an out-of-parameter-space filter
  (285 of 299 kept) prevents extrapolation blow-ups, applied identically to both.
- `luminosity_vs_redshift.png` — plateau-end luminosity L_a vs z, with the flux-limit
  ("limiting luminosity") curve. L_a = 4π D_L(z)² F_a (1+z)^(β−1) (spectral k-correction,
  β = the Beta column). Flat ΛCDM H0=70, Ωm=0.3. The 15 M-est-removed bursts (red) split
  into two groups: a low-z / low-L cluster (z<0.3, logL 42–44) hugging the flux limit —
  the classic Malmquist-bias underluminous nearby bursts — and a few of the highest-z
  bursts (z≈5). F_lim and β_ref are labeled constants at the top of the script; confirm
  these cosmology / flux-limit choices with Prof. Dainotti.
- `luminosity_cut_vs_mest_glm.png` — winning GLM, flagship CV+2σ clip, under the
  M-est cut (drop 15) vs a physical luminosity cut (drop 5). A matched-fraction 5%
  luminosity cut (14 bursts) still only reaches r≈0.62 vs the M-est's 0.66 — the
  M-est lift is not reproduced by a principled Malmquist cut.
- `outlier_train_vs_eval.png` — LEFT: 2x2 showing the outliers' evaluation effect
  (+0.088 r) dwarfs their training effect (+0.021). RIGHT: a model trained on the
  266 clean bursts still predicts the 15 held-out outliers with ~2.8x the RMSE and
  r≈−0.22 — they are intrinsically unpredictable, not fit-poisoning. So the flagship
  gain is mostly from NOT SCORING hard bursts, not from a cleaner fit.
- `outlier_mahalanobis.png` — per-outlier feature-space anomaly. In pure feature
  space only 6/15 are anomalous at p<0.05 (2/15 Bonferroni); 9/15 have ordinary
  features and are removed only because the formula mis-predicts their z (residual,
  not feature, outliers). Adding redshift to the space lifts this to 8/15 (5/15
  Bonferroni), i.e. several are anomalous mainly in z, not in their predictors.
- `flagship_pred_vs_obs_linear_WITHOUT_catoutl.png` — the flagship's OWN committed
  figure (the r≈0.71 "without catastrophic outliers" plot), for reference.
- `flagship_pred_vs_obs_linear_WITH_catoutl.png` — flagship's plot before the 2σ clip.
- `flagship_pred_vs_obs_log_WITHOUT_catoutl.png` — the log10(1+z) version.

## Reconciliation with the published r ≈ 0.71
The published 0.71 is a **single seed-42 split** on the flagship's own `mice20`
imputation (with the ensemble). Reproduced exactly (winning GLM: 0.713). The
`comparison_linear_vs_glm.png` figure reads ~0.66 for the winning GLM because it
(1) **averages over splits** (removing the seed-42 luck, ~−0.03) and (2) uses the
281-row imputation frame needed for the honest full-sample baseline (~−0.02).
Same protocol, and the conclusion is frame-independent: the simple linear tracks
the winning GLM throughout (within ~0.03–0.05), and both fall to r ≈ 0.5 honestly.

## Result (winning GLM shown for reference)

| Model | honest held-out r | flagship CV+2σ-clip r |
|---|---|---|
| Simple linear (7 features, additive) | **0.46** | **0.63 ± 0.03** |
| Winning GLM (Eq. winner) | 0.53 | 0.66 ± 0.03 |

- `honest held-out r`: mean over 100 random 80/20 splits, refit per split, scored
  on the untouched 20%. No cut, no clip.
- `flagship CV+2σ-clip r`: the flagship protocol, averaged over 20 random 80/20
  splits (see note on splits below).

## Conclusion

A plain additive linear model reaches **~0.63** under the flagship protocol,
statistically indistinguishable from the elaborate winning GLM's **~0.66**, while
its *honest* held-out performance is only **~0.46** (the winning GLM genuinely
edges it, 0.53 vs 0.46 — the interaction/square terms buy ~0.07 of real signal,
not ~0.20). So the ~0.15 lift is the **evaluation protocol** — the M-estimator
cut removes generically hard-to-predict bursts (which helps *any* model) and the
2σ clip drops each model's own worst residuals — **not** the specific model, and
not a circular reuse of the winning formula.

## Notes on reproducibility

- The flagship uses a **single** 80/20 split; that one split is arbitrary and
  swings the number by ~0.05–0.08 (a favorable draw gives ~0.67; the average is
  ~0.63). This script **averages over 20 splits** so the figure is stable and
  language-independent — an R port gives the same averages. That single-split
  instability is itself a weakness of the flagship's one-split protocol.
- Absolute levels also shift ~0.02 with the imputation frame; the
  simple-vs-winning comparison and the honest-vs-protocol gap are robust to both.
