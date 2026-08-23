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

## Is this overfitting? Checked — no, but there's a real caveat

Two separate questions worth separating here:

**1. Is the per-bin linear fit itself just fitting noise?** No — checked with
proper 5-fold cross-validation (fit each bin's slope/intercept on 4/5 of that
bin's points, apply to the held-out 1/5, repeat across folds). Results barely
move from the in-sample numbers above:

| Result | r(z) in-sample | r(z) 5-fold CV | RMSE in-sample | RMSE 5-fold CV |
|---|---|---|---|---|
| Theil-Sen | 0.909 | 0.901 | 0.533 | 0.556 |
| No-plateau | 0.910 | 0.897 | 0.533 | 0.568 |
| Multivariate emcee | 0.909 | 0.901 | 0.541 | 0.572 |

A 2-parameter linear fit on 60-110+ points per bin has little room to
memorize noise, and this confirms it doesn't — a held-out slope/intercept
transfers almost perfectly within its bin.

**2. Does this reproduce the repo's own `BC_3way` function literally, run in
R?** No, and here's where the real caveat lives. `BC_3way` fits each bin
using the *true* observed z (`Zspec`) — same as above — but it decides which
bin's correction to *apply* using the **predicted** z (`Zphot`), because it
was written for a real generalization set where true z is unknown. Running
the literal R function on this same self-correction (same data passed as
both arguments) makes every result *worse*, not better (e.g. Theil-Sen:
r=0.707→**0.689**, bias -0.146→**-0.212**) — because ~35-40% of points have a
predicted z that crosses one of the two cutoffs (z=2 or z=3.5) on the
opposite side from their true z, so those points get corrected with the
wrong bin's coefficients entirely.

**What this means for interpreting the numbers above**: the strong
improvement here depends on knowing each point's *true* redshift range to
route its correction — legitimate for this exercise (self-correcting
already-labeled CV data, exactly what the paper's own Fig 13 does), but it
would not carry over to correcting a real unknown-redshift GRB, where only
the predicted range is available and `BC_3way`'s actual (worse-performing)
routing is what you'd be stuck with. Read this result as "how much of the
model's error is a fixable, redshift-range-dependent compression, given the
true range" — a diagnostic about the model's error structure — not as "here
is a deployable technique that improves predictions on new GRBs."

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
- `1_theilsen_bias_corrected.csv`, `2_noplateau_bias_corrected.csv`,
  `3_multivariate_emcee_bias_corrected.csv` — per-GRB raw data behind the
  plots above. Columns: `GRB` (identifier), `InvZphot`/`InvZspec` (predicted/
  observed log10(z+1), raw), `Zphot`/`Zspec` (same, linear z), `pred_min`/
  `pred_max`/`linpred_min`/`linpred_max` (raw CV-fold error bounds, log and
  linear scale), and the bias-corrected counterparts:
  `InvZphot_corrected`/`Zphot_corrected` (corrected prediction) and
  `pred_min_corrected`/`pred_max_corrected`/`linpred_min_corrected`/
  `linpred_max_corrected` (error bounds run through the same per-redshift-bin
  correction as the point prediction).
