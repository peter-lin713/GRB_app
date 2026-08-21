# Cross-validated results — predicted vs. observed redshift

Raw per-GRB CV output backing the three `paper_candidate_results/plots/*.png`,
copied directly from each result's `Results/10fCVResults_without_catOutl_
correlation_plot.csv` (the "without categorical/non-negotiable outliers"
sample — the full N used for every headline stat quoted in this repo, before
any 1σ/2σ cone filtering).

| File | N | Source run |
|---|---|---|
| `1_theilsen_singlevar_10fCVResults.csv` | 206 | `runs/theilsen_cleaned_linearz_5pct_mice20/` |
| `2_noplateau_10fCVResults.csv` | 204 | `combination_scheme_comparison/no_plateau_formula_run/` |
| `3_multivariate_emcee_10fCVResults.csv` | 202 | `plateau_ablation_experiments/runs/multivariate_emcee_reused_formula/` |

## Columns

- `InvZphot`, `InvZspec` — predicted / observed log10(z+1)
- `Zphot`, `Zspec` — predicted / observed redshift (linear scale)
- `Dlogz`, `Dz`, `normDz` — residuals: log-scale, linear-scale, and
  linear-scale normalized by `(1+Zspec)`
- `pred_min`/`pred_max`, `linpred_min`/`linpred_max` — per-GRB CV fold
  min/max prediction range (log scale and linear scale), used as the error
  bars in `plots/*.png`

## Cone filter (matches the plots)

The plots hide points outside each result's own 2σ cone
(`Result_plot_maker.R`'s convention): `L_Sigma = sd(Dlogz)`, a point is
"inside the 2σ cone" if `abs(Dlogz) < 2*L_Sigma`. All rows are included
here unfiltered — apply that formula directly if you need the same
outlier-free subset used for the plots.
