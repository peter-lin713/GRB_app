# Bias correction — applied to the training/CV results directly

Uses the same method as `Generalization/Bias_Correction_function.R`'s `BC_1way`
(a quantile-quantile linear correction: sort predicted and observed log10(z+1)
independently, fit `observed ~ predicted` on the sorted pairs, apply that fit back
to the raw predictions) — but self-applied to each result's own cross-validated
predictions, not extended to a separate unknown-redshift generalization set.

## Source data

All three now use their exact, correct source runs. The Theil-Sen main result's
underlying run (`runs/theilsen_cleaned_linearz_5pct_mice20`) had been deleted
during an earlier cleanup pass (mistaken for a redundant duplicate of
`theilsen_cleaned_linearz_5pct`) — recovered from `backup-pre-cleanup-2026-08-10`
and restored to the working tree, confirmed as an exact match (r=0.7071, N=206,
Sigma=0.895, RMS=0.91, Bias=0.15, NMAD=1.22 — identical to the reported "0.707"
figure). The other two match their `paper_candidate_results/` plots directly
(no-plateau: r=0.644; multivariate emcee: r=0.714≈0.7136).

## Result: bias correction trades bias for variance, consistently

| Result | N | r before→after | RMSE before→after | Bias before→after | NMAD before→after |
|---|---|---|---|---|---|
| Theil-Sen single-var (MAIN RESULT) | 206 | 0.707→0.704 | 0.907→0.988 | -0.146→**-0.001** | 0.213→0.204 |
| No-plateau ablation | 204 | 0.644→0.642 | 0.997→1.097 | -0.201→**-0.009** | 0.238→0.242 |
| Multivariate emcee | 202 | 0.714→0.702 | 0.925→1.008 | -0.207→**-0.006** | 0.195→0.196 |

Same pattern in all three: the correction all but **eliminates the systematic bias**
(mean signed error shrinks to near zero every time), but r and RMSE both get
slightly *worse* — except NMAD, which actually improves for the main Theil-Sen
result (0.213→0.204) even though it ticks up marginally for the other two.
Visually (see the `*_before_after.png` plots), the raw model under-predicts at
high observed z — a compressed-spread effect common to ML regressors. The linear
correction stretches the whole prediction distribution to match, which fixes the
aggregate bias but also amplifies scatter for the already-uncertain high-z tail
(a few points overshoot substantially post-correction).

**Takeaway for the paper**: bias correction is a real, standard bias-variance
tradeoff here, not a free improvement — worth mentioning if discussing bias
specifically, but it would not be honest to report the corrected numbers as a
better result than the raw ones without also disclosing the RMSE/r cost.

## Files

- `summary_table.png` — the before/after comparison table
- `1_theilsen_before_after.png`, `2_noplateau_before_after.png`,
  `3_multivariate_emcee_before_after.png` — predicted-vs-observed scatter,
  before and after, side by side
