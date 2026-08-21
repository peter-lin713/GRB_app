# Bias correction — applied to the training/CV results directly

Uses the same method as `Generalization/Bias_Correction_function.R`'s `BC_1way`
(a quantile-quantile linear correction: sort predicted and observed log10(z+1)
independently, fit `observed ~ predicted` on the sorted pairs, apply that fit back
to the raw predictions) — but self-applied to each result's own cross-validated
predictions, not extended to a separate unknown-redshift generalization set.

## Source data caveat

The exact model/CV-results underlying the reported **0.707** Theil-Sen main result
were never saved to disk, so bias correction there uses the closest available
methodologically-equivalent run instead (`daume_fix_comparison/runs/
theilsen_single_WITHOUT_daume`, r=0.628 — same standing formula, same "no working
Daume" setup, built at a different point in the session so the exact M-estimator
cut/MICE draw differs). The other two match their `paper_candidate_results/`
plots exactly (no-plateau: r=0.644; multivariate emcee: r=0.714≈0.7136).

## Result: bias correction trades bias for variance, consistently

| Result | N | r before→after | RMSE before→after | Bias before→after | NMAD before→after |
|---|---|---|---|---|---|
| Theil-Sen single-var (proxy) | 215 | 0.628→0.626 | 0.919→1.034 | -0.135→**-0.002** | 0.225→0.231 |
| No-plateau ablation | 204 | 0.644→0.642 | 0.997→1.097 | -0.201→**-0.009** | 0.238→0.242 |
| Multivariate emcee | 202 | 0.714→0.702 | 0.925→1.008 | -0.207→**-0.006** | 0.195→0.196 |

Same pattern in all three: the correction all but **eliminates the systematic bias**
(mean signed error shrinks to near zero every time), but r and RMSE both get
slightly *worse*, and NMAD ticks up marginally. Visually (see the `*_before_after.png`
plots), the raw model under-predicts at high observed z — a compressed-spread
effect common to ML regressors. The linear correction stretches the whole
prediction distribution to match, which fixes the aggregate bias but also
amplifies scatter for the already-uncertain high-z tail (a few points overshoot
substantially post-correction).

**Takeaway for the paper**: bias correction is a real, standard bias-variance
tradeoff here, not a free improvement — worth mentioning if discussing bias
specifically, but it would not be honest to report the corrected numbers as a
better result than the raw ones without also disclosing the RMSE/r cost.

## Files

- `summary_table.png` — the before/after comparison table
- `1_theilsen_before_after.png`, `2_noplateau_before_after.png`,
  `3_multivariate_emcee_before_after.png` — predicted-vs-observed scatter,
  before and after, side by side
