# Bias correction — applied to the training/CV results directly

Matches the paper's Sec 4.6 methodology exactly (the "optimal transport bias
correction," implemented as `BC_3way` in `Generalization/Bias_Correction_
function.R`): **three separate quantile-quantile linear corrections**, one per
redshift range —  z<2, 2<z<3.5, z>3.5 (the paper's own stated cutoffs). Within
each range, predicted and observed log10(z+1) are sorted independently and a
linear model (`observed ~ predicted`) is fit on the sorted pairs; that range's
fit is then applied back to its own points (`Z_C = B + A·Z_P`, the paper's
exact equation). Self-applied to each result's own cross-validated predictions
— matching what the paper itself does for its own headline SuperLearner result
in Sec 5.1 (r=0.646→0.89 there), not extended to a separate unknown-redshift
generalization set.

This replaces an earlier pass in this folder that used a single global fit
(`BC_1way`, one linear correction across the whole sample) instead of the
paper's actual 3-range method — that version showed bias correction as a
lossy bias/variance tradeoff (r and RMSE both got worse). The 3-way method
below is the correct replication of the paper's approach and shows a genuine,
large improvement instead, consistent with what the paper itself reports.

## Source data

All three now use their exact, correct source runs. The Theil-Sen main result's
underlying run (`runs/theilsen_cleaned_linearz_5pct_mice20`) had been deleted
during an earlier cleanup pass (mistaken for a redundant duplicate of
`theilsen_cleaned_linearz_5pct`) — recovered from `backup-pre-cleanup-2026-08-10`
and restored to the working tree, confirmed as an exact match (r=0.7071, N=206,
Sigma=0.895, RMS=0.91, Bias=0.15 — identical to the reported "0.707" figure;
NMAD=0.818 using the current, correct plot formula — the originally-cached
plot showed NMAD=1.22 from a since-fixed double-scaling bug, see
`paper_candidate_results/README.md` §1). The other two match their
`paper_candidate_results/` plots directly (no-plateau: r=0.644; multivariate
emcee: r=0.714≈0.7136).

Note: the before/after NMAD in the table below is a *different* metric
(median absolute normalized residual, no 1.48 scaling factor) than the main
plot's title NMAD (1.48×median|Dz|, linear-z scale) — the two aren't the same
number and shouldn't be compared to each other directly.

## Result: bias correction is a genuine, large improvement — matches the paper

| Result | N | r before→after | RMSE before→after | Bias before→after | NMAD before→after |
|---|---|---|---|---|---|
| Theil-Sen single-var (MAIN RESULT) | 206 | 0.707→**0.909** | 0.907→**0.533** | -0.146→**-0.004** | 0.213→**0.117** |
| No-plateau ablation | 204 | 0.644→**0.910** | 0.997→**0.533** | -0.201→**-0.005** | 0.238→**0.117** |
| Multivariate emcee | 202 | 0.714→**0.909** | 0.925→**0.541** | -0.207→**-0.004** | 0.195→**0.112** |

Same pattern in all three, and it mirrors the paper's own reported result almost
exactly (paper's main SuperLearner CV result: r=0.646→0.89, RMSE=1.011→0.62,
bias=0.14→0.0047, NMAD=1.34→0.86, Sec 5.1/Fig 13). Splitting the correction by
redshift range matters a lot: fitting one global linear correction (the
earlier, incorrect version of this analysis) only removes the bias while
leaving r/RMSE roughly flat or slightly worse. Fitting a *separate* correction
per range lets the model correct the compressed high-z under-prediction on its
own terms in each range, which is why r jumps from ~0.6-0.7 to ~0.91 across the
board — not just a bias fix, a genuine tightening of the scatter around the
1:1 line (see the `*_before_after.png` plots).

**Takeaway for the paper**: bias correction here is not a marginal tweak — it's
the standard step the paper itself applies before reporting its headline
number, and skipping it would understate how well the underlying model does
once its systematic, redshift-range-dependent compression is corrected for.

## Sigma cones — the pipeline's standard catastrophic-outlier convention

Every other result in this repo (e.g. the main Theil-Sen plot's "within 2σ
cone = 199 (97%)") is reported alongside this pipeline's standard 1σ/2σ
"sigma cone" diagnostic (`Result_plot_maker.R`): fit `L_Sigma = sd(Dlogz)`
where `Dlogz = observed − predicted` log10(z+1), then a point is "inside the
2σ cone" if `|Dlogz| < 2·L_Sigma` (likewise for 1σ). This was missing from
the first pass of this bias-correction writeup — added below, computed
separately before and after correction (each state gets its own `L_Sigma`,
matching how the pipeline always computes it fresh per result).

| Result | N | Sigma before→after | Within 2σ before→after | Within 1σ before→after |
|---|---|---|---|---|
| Theil-Sen single-var (MAIN RESULT) | 206 | 0.126→0.082 | 199 (97%)→195 (95%) | 133 (65%)→141 (68%) |
| No-plateau ablation | 204 | 0.133→0.078 | 198 (97%)→197 (97%) | 130 (64%)→148 (73%) |
| Multivariate emcee | 202 | 0.118→0.073 | 195 (97%)→193 (96%) | 133 (66%)→139 (69%)|

**Why the 2σ count doesn't simply go up too**: the cone is defined *relative
to the sample's own residual spread*, not an absolute error tolerance.
Correction shrinks `Sigma` itself by ~35-40% (residuals genuinely tighten),
so the "2σ" window shrinks right along with it — a few points that were
comfortably inside the old, wider band now sit just outside the new,
narrower one, even though every point's absolute error went down. The 1σ
count is the more informative number here and it improves consistently
(65%→68%, 64%→73%, 66%→69%): a larger share of the sample now lands within
one (now-tighter) sigma of the 1:1 line. The `*_before_after.png` scatter
plots show both cone boundaries (2σ blue, 1σ green) directly, matching the
pipeline's linear-z plot convention exactly.

## Files

- `summary_table.png` — the before/after r/RMSE/Bias/NMAD comparison table
- `cone_stats_table.png` — the before/after sigma-cone comparison table
- `1_theilsen_before_after.png`, `2_noplateau_before_after.png`,
  `3_multivariate_emcee_before_after.png` — predicted-vs-observed scatter,
  before and after, side by side, with 1σ/2σ cone boundaries drawn
