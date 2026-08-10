# Daume domain-adaptation: bug, fix, and results

Everything from the "is the Daume `is_optical` domain-adaptation mechanism actually working"
investigation, compiled in one place.

## Start here

`plots/daume_comparison_table.png` — genuine Daume vs. no Daume, across all 4 calibration
schemes that carry an `is_optical` column, same standing formula throughout, same
"Without 2σ outliers" subset.

**Headline finding: genuine Daume domain-adaptation never meaningfully helps.** Across all
four datasets it's either a wash or a slight regression:

| Dataset | Daume | N | r(z) |
|---|---|---|---|
| Standing (single-var emcee) | Yes | 208 | **0.670** |
| Standing (single-var emcee) | No | 206 | 0.663 |
| Multivariate emcee | Yes | 204 | 0.711 |
| Multivariate emcee | No | 202 | **0.714** |
| Theil-Sen multivariate | Yes | 207 | 0.678 |
| Theil-Sen multivariate | No | 205 | **0.687** |
| Theil-Sen single-var | Yes | 216 | 0.614 |
| Theil-Sen single-var | No | 215 | **0.628** |

For context, none of these numbers are close to the *historically reported* ~0.70-0.71 for
these same datasets — that's because those historical numbers were themselves computed with
a **broken** Daume mechanism (see below), which turns out to be functionally identical to no
Daume at all, just via a different, more circuitous route to the same place. The small drop
from "historical" to "no Daume" here (e.g. standing: 0.712 historical → 0.663 here) reflects
the M-estimator cut / MICE draw differing between whenever each was originally run, not
anything about Daume.

## The bug: Daume was dead code the entire time

`superlearner.R`'s `is_optical` block computed `is_optical` + 10 columns of
`value × is_optical × 0.25` (the intended domain-adaptation mechanism), printed a confirming
message ("Added is_optical + N Daume-scaled domain columns"), and even wrote them into the
diagnostic `grb_xray_imputed.csv` snapshot — but never used them. `lassovar` (the predictor
name list) is computed once, *before* the Daume block runs, and never gets updated to include
the new columns; the actual design matrix fed to `SuperLearner()` (`O1Predictors <-
subset(GRBPred, select = lassovar)`, further down) rebuilds `GRBPred` from `lassovar` alone,
silently discarding everything the Daume block added.

**Confirmed empirically**, not just by code trace: re-ran the multivariate-emcee dataset with
`is_optical` present vs. absent through the *original* (buggy) script — the two runs produced
**byte-identical** predictions (r=0.713572 to 6 decimal places, same 202 GRBs). A real ensemble
of 17 stochastic learners cannot coincidentally converge to the same output unless the extra
columns genuinely never entered the fit.

This means **every "WITH is_optical" result reported anywhere earlier in this project** —
including the "Theil-Sen WITH is_optical (0.707) vs WITHOUT is_optical (0.679)" comparison —
was never actually testing Daume. Whatever real difference existed between those specific runs
came from something else (most likely which M-estimator cut / MICE draw each used), not from
`is_optical`'s presence.

**One exception, a different mechanism entirely**: hard-zero-domain (v1/v2) and
multivariate-emcee "Dataset B" used literal native `Alpha_opt`/`Beta_opt`/etc. *value* columns
via separately-modified scripts (`superlearner_hardzero.R`/`superlearner_optadd.R`) that
explicitly extend `features_for_mice_preds` themselves — not the `is_optical`-flag/Daume-
weighting path this investigation is about. Those genuinely did use their extra columns
(caught and fixed independently, earlier in the project, via the "buggy Dataset B" identical-
output bug).

## The fix

- **`superlearner.R`** (main script, used for all "no Daume" work going forward): the dead
  Daume block removed entirely. Confirmed via A/B test to be a pure no-op — re-ran the same
  `is_optical`-containing dataset through both the pre-fix checkpoint and the fixed script;
  byte-identical output (r=0.692030 to 6 decimals, N=203, same GRB set).
- **`superlearner_CHECKPOINT_pre_daume_fix.R`**: exact snapshot of `superlearner.R` from
  before the fix, preserved for provenance (every historical result before this batch traces
  to this version).
- **`superlearner_daume_working.R`**: a genuinely-working version. Extends `lassovar_final <-
  c(lassovar, "is_optical", colnames(opt_design))` right after building the Daume columns, and
  uses `lassovar_final` (not `lassovar`) at the design-matrix-rebuild step, so the columns
  actually reach `Predictors`. Also happens to fix a latent inefficiency: the formula file
  always contains the same formula 3x (an artifact of how it's installed), so the old code
  built 3 identical, redundant learners; the reader now dedupes to 1.

## A debugging detour worth keeping honest

Getting `superlearner_daume_working.R` to actually run cleanly took several failed attempts,
useful to record so the same time isn't lost twice:

1. First attempt crashed instantly on every dataset with `Error in attr(x, "formula") %||% {
   : invalid formula`, preceded by `read.table()` warnings ("incomplete final line",
   "EOF within quoted string") on `Best_formula_GAM.txt`.
2. **Misdiagnosed as an iCloud file-provider race** (this project directory lives under
   iCloud Desktop sync, confirmed via `com.apple.file-provider-domain-id` on `~/Desktop` and
   `brctl status` showing `com.apple.CloudDocs` actively syncing at the exact crash timestamp)
   — plausible on the surface, and a retry-with-sleep wrapper around `read.table()` was added.
   Retries didn't help even over 5+ seconds.
3. Switched to a `readLines()` + regex parser to bypass `read.table()`/`scan()` entirely —
   **also failed identically**, which ruled out "read.table's parser specifically" as the
   cause and should have been the signal to stop chasing environmental causes.
4. **True root cause**, found by directly inspecting the file mid-crash: it was
   double-wrapped — literal `[1] "[1] "Response ~ ...` — because the orchestration script
   (`run_daume_working_batch_v2.sh`) backed up the *already-wrapped* `Best_formula_GAM.txt`
   and then called `install_formula()` (which wraps its input in `[1] "..."`) on that backup,
   wrapping it a second time. The formula was already correctly installed and needed no
   reinstallation at all — the fix was to delete the erroneous `install_formula()` call, not
   to make the reader more robust to a genuinely corrupt file. The `readLines()`/regex reader
   is kept anyway (harmless, marginally simpler than `read.table()`), but it was never the
   actual fix.

## Folder guide

- `runs/<dataset>_WITH_daume/` — full SuperLearner `Plot_Output`/`Results` for each dataset
  with genuinely-working Daume
- `runs/<dataset>_WITHOUT_daume/` — same 4 datasets, `is_optical` dropped entirely before
  MICE/M-estimator (built via `scripts/prep_no_isoptical.R`)
- `scripts/` — `superlearner_no_daume_current.R` (= repo-root `superlearner.R`, fixed),
  `superlearner_daume_working.R` (genuine Daume), `superlearner_CHECKPOINT_pre_daume_fix.R`
  (pre-fix snapshot), `prep_no_isoptical.R` (the no-Daume prep pipeline)
- `plots/daume_comparison_table.png` — the headline table
- `plots/linear_predictor_plots/` — all 8 predicted-vs-observed plots, paired and numbered
  (1-2 = standing, 3-4 = multivariate emcee, 5-6 = Theil-Sen multivariate, 7-8 = Theil-Sen
  single-var; odd = with Daume, even = without)
- `build_comparison_table.py` — regenerates the table image

## Reproducing a "WITH Daume" run

```bash
# 1. Build the no-Daume cut (M-estimator using the standing formula):
Rscript scripts/prep_no_isoptical.R <source_csv> <method_tag> <out_dir>

# 2. Merge is_optical back in from the source CSV (see the merge logic used for
#    GAM_supercomputer/*_no_isoptical/final_outliers_removed_WITH_isoptical.csv)

# 3. Run with genuine Daume (standing formula already installed in Best_formula_GAM.txt):
Rscript scripts/superlearner_daume_working.R <final_outliers_removed_WITH_isoptical.csv> \
    TRUE FALSE FALSE FALSE 0.65 10 <out_dir> 0.05 TRUE hard
```

`Best_formula_GAM.txt`/`Best_formula_GLM.txt` at the repo root are left on the standing
formula by default — confirmed intact after every run in this batch.
