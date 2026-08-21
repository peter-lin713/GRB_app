# Candidate results for the paper

Three results pulled out of the full `combination_scheme_comparison/` investigation,
specifically the ones flagged as worth publishing. All three use the "Without 2σ
outliers" subset and the pipeline described in the main methodology (MICE `midastouch`
m=20/maxit=20, non-negotiable + M-estimator outlier cuts, 17-learner SuperLearner
ensemble, 10×10 CV with 100-draw MC error propagation).

**Important context established since these were first proposed**: the `is_optical`
Daume domain-adaptation column present in the Theil-Sen dataset never actually reached
the model (a dead-code bug in `superlearner.R`, since fixed and separately investigated
in `daume_fix_comparison/` — genuinely enabling it doesn't help anyway). So "Theil-Sen
WITH is_optical" is not a domain-adaptation result; it's the standard 10-feature model
alone, fed by Theil-Sen-calibrated inputs. Frame it that way in the writeup.

## 1. Main result — Theil-Sen single-variable calibration

**r(z) = 0.707, N=206**, Sigma=0.895, RMS=0.91, Bias=0.15, NMAD=0.818.
`plots/1_MAIN_theilsen_singlevar_r0.707.png`

Rounds to the same 0.71 already in the paper, and has better bias/sigma/RMSE than
both the OT-fusion scheme and the standing single-variable emcee baseline (though
NMAD is not the best of the schemes tested — multivariate emcee's, 0.795, is
marginally lower). Reuses the standing formula (`formulas/standing_formula_
reused_by_theilsen_and_multivariate_emcee.txt`) rather than a freshly-searched
one — own-formula search on this same dataset scored lower (0.681).

**NMAD correction (this pass):** this result's plot previously showed NMAD=1.22.
That number came from a stale double-scaling bug in `Result_plot_maker.R`
(`1.48*mad(Dz)` — R's `mad()` already applies its own 1.4826 factor, so the
extra `*1.48` double-counted it), fixed in this repo's history (`git log -p --
Result_plot_maker.R` shows the exact diff) but never re-propagated to this
specific cached result. The other two candidates below were already regenerated
post-fix (their displayed NMAD already matched recomputation). Corrected value:
0.818, using the current, correct formula (`1.48*median(abs(Dz))`) — confirmed
by regenerating the plot directly from the CV results CSV.

Suggested framing: not "we improved the number" (0.707 is not higher than the
paper's existing 0.712/0.70), but "we replaced the optical→X-ray calibration
machinery (emcee MCMC) with a simpler, more standard robust estimator (Theil-Sen)
and lost nothing" — a robustness result, not an accuracy claim.

## 2. Ablation — predicting without the plateau parameters

**r(z) = 0.644, N=204**, Sigma=0.977, RMS=0.99, Bias=0.20, NMAD=0.947.
`plots/2_ablation_no_plateau_params_r0.644.png`

Alpha, Beta, log10Fa, log10Ta (the plateau-fit parameters) dropped entirely — only
the 6 shared/direct-measurement variables remain (log10T90, Gamma, log10Fluence,
PhotonIndex, log10NH, log10PeakFlux). Own fresh formula
(`formulas/no_plateau_own_fresh_formula.txt`), no opt cols.

**Be precise about what this shows**: the correlation *does* drop measurably
(~0.06-0.07 below baseline, not negligible) — the honest claim is "the model still
predicts reasonably without the plateau parameters, at a real but modest cost,"
not "the correlation doesn't go down." NMAD also gets modestly worse here
(0.947 vs the corrected baseline's 0.818) — an earlier version of this note
claimed NMAD improved, which was an artifact of comparing against the baseline's
since-corrected 1.22 figure (see note in §1); with the correct baseline, this
ablation costs a bit on both r and NMAD, consistently.

## 3. Best basic fit with no domain-adaptation columns at all

**r(z) = 0.714, N=202**, Sigma=0.901, RMS=0.92, Bias=0.19, NMAD=0.792.
`plots/3_best_basic_fit_multivariate_emcee_r0.714.png`

Multivariate emcee (MCMC) calibration — each X-ray-scale plateau parameter fit from
all four optical parameters jointly rather than just its own counterpart — reusing
the standing formula. The single best number of any scheme tested with no `is_optical`/
opt-columns involved at all; ties/edges the paper's own baseline (0.712).

## Files

- `plots/summary_table.png` — the 3-row comparison table
- `plots/1_*.png`, `plots/2_*.png`, `plots/3_*.png` — cross-validated predicted-vs-observed
  plots (linear z, without 2σ outliers), full stats embedded in each title
- `formulas/` — the two formulas involved (standing, reused by both #1 and #3; and
  the no-plateau ablation's own fresh formula for #2)

Full provenance, the complete 9-scheme comparison these were selected from, and the
bugs found/fixed along the way live in `combination_scheme_comparison/README.md`.
