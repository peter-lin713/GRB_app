# GRB Redshift Prediction Pipeline — Methodology

Documents the current, fully-corrected pipeline for both combination methods
(emcee-projection and OT-fusion), the bugs found and fixed along the way, and
where to find the outputs. Written alongside the final production runs of
`emcee_v2_cleaned_formula_generation` and `ot_v3_cleaned_formula_generation`.

---

## 1. Raw catalogs

| File | Description | GRBs |
|---|---|---|
| `Data/xray_data.csv` | Unfiltered X-ray catalog | 251 |
| `Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv` | Filtered X-ray catalog (error-cut + MICE) | 223 |
| `Data/OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv` | Optical catalog, emcee-combination variant | 161 |
| `Data/OnlyLGRBs_data_171_optical_corrected.csv` | Optical catalog, OT-combination variant | 161 |
| `Data/optical_data.csv` | Raw/uncleaned optical catalog (not used in combination; historical) | 171 |

Overlap: filtered X-ray ∩ filtered optical = 78 GRBs. Unfiltered X-ray ∩
filtered optical = 86 GRBs. GRB ID suffixes matter (e.g. `A`/`B` distinguish
otherwise-identical designations) — both combination scripts normalize IDs by
stripping the `GRB` prefix and appending `A` if no letter suffix is present.

---

## 2. Combination methods

Two independent ways of building one combined (X-ray + optical) training set.
Both start from the same two catalogs (filtered X-ray + the relevant optical
variant) and produce a 296-GRB combined, error-cut dataset with an
`is_optical` flag marking which GRBs came from the optical catalog only.

### 2a. emcee-projection (`combine_opt_xray_emcee/combine_optical_xray_emcee_v2.py`)

For each of the four "prompt" light-curve parameters (log10Fa, log10Ta,
Alpha, Beta), fits a linear calibration (`y = m·x + b` plus intrinsic
scatter) between the optical-band and X-ray-band measurements, using only the
~74-90 GRBs with *both* measurements (the overlap set). The MCMC posterior
over (m, b, scatter) is then used to project every optical-only GRB's
measurement into X-ray-equivalent terms, with full posterior-predictive
uncertainty propagation (perturbing each optical measurement within its own
error bar, for every posterior sample, then taking the median and 16th-84th
percentile spread as the point estimate and 1σ error).

**Reproducibility fix (this session):** the MCMC fit previously ran as a
single unseeded chain — every invocation gave a different posterior. Fixed by
running `N_CHAINS = 4` independent chains per parameter, each with a fixed,
deterministic seed (`seed = 1000*(param_index+1) + chain_index`), and pooling
all chains' post-burn-in samples into one larger posterior sample. This is
standard MCMC practice (multiple chains are normally run anyway, to check
convergence) and, combined with the seeding, makes the whole calibration
exactly reproducible run to run. Verified: two consecutive runs of the
combination script now produce byte-identical output.

Output: `Data/superlearner_training_emcee_v2_errcut_relative.csv` (296 GRBs),
`Data/emcee_chains_v2.npz` (the pooled posterior, also used by OT — see
below), `Data/emcee_projection_draws.csv` (4 sampled (m, b) draws per
parameter, for the SuperLearner-stage recovery mechanism — see §4).

### 2b. OT-fusion (`ot_fusion/combine_optical_xray_ot_v3.py`)

Imputes five X-ray-only features (Gamma, PhotonIndex, log10NH, log10Fluence,
log10PeakFlux) for optical-only GRBs via donor-matching: each optical-only
GRB is matched against the X-ray-native GRBs in a standardized space of the
four calibrated Dainotti parameters (logFa_x, logTa_x, Alpha_x, Beta_x — the
*projected*, X-ray-scale versions, using the emcee calibration above), via a
Gibbs/softmax kernel (weight ∝ exp(-cost/bandwidth), no donor-marginal
constraint — a fix from earlier v1/v2 versions, which used balanced Sinkhorn
OT and had a structural degeneracy for single-query matching). The kernel
bandwidth is chosen by a leave-one-out sweep over candidate bandwidths,
picking whichever gives the best LOO correlation on the anchor (overlap)
GRBs. The weighted variance of the matched donor pool becomes each imputed
value's own uncertainty ("barycentric" error), capped at 3× the X-ray
population's own spread for that feature.

**Reproducibility:** OT's own matching and bandwidth-selection logic is
fully deterministic (no randomness anywhere in `ot_impute()` or the LOO
sweep) — but it depends on `Data/emcee_chains_v2.npz` for the coordinate
space it matches in, and that file used to change every time the emcee
script ran (see above). With the emcee-side fix, this file is now stable,
which makes OT's output stable too. Verified: two consecutive runs of the OT
combination script now produce byte-identical output.

Output: `Data/superlearner_training_ot_v3_errcut_relative.csv` (296 GRBs).

### `is_optical` flag

Added identically for both combination outputs: a GRB is flagged
`is_optical = 1` if its normalized ID is absent from the filtered X-ray
catalog's GRB list — a direct catalog-membership check (not a
missingness-count heuristic). 73 of 296 GRBs are flagged this way in both
combination datasets.

---

## 3. Cleaning: non-negotiable cuts (full removal)

Applied after the long-GRB scope filter (T90 > 2s: 296 → 285 for both
methods). Physically-impossible or unreliable GRBs are **removed from the
dataset entirely** — not nulled and imputed. Criteria:

- `Alpha > 3`
- `Beta > 2`
- `log10T90 > 6`
- `PhotonIndex < 0`
- Relative error exceeds half the value (`err/|value| > 0.5`), checked
  **only on Alpha, Beta, log10Fa, log10Ta**.

**The error-ratio check is deliberately restricted to Alpha/Beta/log10Fa/
log10Ta.** These are the four parameters emcee's MCMC calibration actually
fits **for the optically-sourced subset** (`is_optical == 1`, ~73 GRBs); T90,
PhotonIndex, Fluence, and NH are direct catalog measurements for every GRB,
so a large relative error there reflects genuine measurement uncertainty, not
an unreliable projection fit — applying the same err>value logic to them
conflates two different kinds of "uncertain." (`PeakFlux` was excluded
earlier still, for a separate reason: 94% of its flags traced to a numerical
artifact — `log10PeakFlux` crosses zero for roughly a third of the sample,
inflating the ratio near-zero-denominator regardless of actual data quality,
6 literal `Inf` values observed for emcee; OT's remaining flags had no clear
explanation. Excluded entirely rather than accepting a ~30% cut on weak
grounds.)

**Important scoping note:** the check runs on all 285 GRBs, not just the
`is_optical == 1` subset it was motivated by. For X-ray-native GRBs
(`is_optical == 0`), Alpha/Beta/Fa/Ta are also direct measurements, not
emcee-fit values — a flag there means "poorly-constrained measurement,"
not "unreliable emcee projection." Related but distinct kinds of
unreliability, on the same four parameters. This isn't just a hypothetical:
as of the current data, all 4 flagged GRBs are `is_optical == 0` — none of
the currently-flagged cases are actually catching an unreliable emcee fit.
The check is being kept as-is (broadly applied) rather than restricted to
`is_optical == 1` only, since a poorly-constrained direct measurement of
these parameters is still a reasonable thing to screen for — but the
rationale given above (emcee's fit) only strictly applies to part of what
this criterion actually catches.

**Self-regenerating, not a stale manual file.** `generate_nonnegotiable_
outliers.R` (repo root) is the single source of truth — a callable function
(`update_nonnegotiable_outliers()`), sourced fresh by each formula-generation
script at the start of every run, rather than trusting a pre-computed static
file. (This replaces an earlier design where `candidate_outliers_*.csv` was
generated once by hand and just read thereafter — it went stale without
anything noticing, at one point flagging only 1 GRB when a proper
recomputation on the same data found several more.) `confirmed_outliers_to_
drop.txt` unions the freshly-detected automatic candidates into whatever's
already there each run, so manual additions from visual review of the
Alpha/Beta/log10Fa/log10Ta 4D fundamental-plane plots persist across reruns,
but the automatic part can never silently drift out of sync again.

Current counts (emcee v2, corrected criteria): 4 GRBs flagged
(090516A, 100621A, 131117A, 230328B — one hard PhotonIndex<0 cut, three
err>value:Beta) → 285 → **281**. OT v3 has **not yet been rerun** with this
corrected (Alpha/Beta/Fa/Ta-only) criterion — its last run used the broader,
now-superseded err>value check across all seven columns; treat OT's numbers
below as stale pending a rerun.

---

## 4. Real-value recovery (SuperLearner-stage only)

A separate mechanism, used only in `superlearner_spencer.R` (not in the
formula-generation stage), confirmed to match spencer's actual
`cv_maxstack_v2.R` methodology on manual inspection:

1. **Point recovery**: for every `is_optical == 1` GRB, pull real measured
   values (log10T90, log10Fluence, PhotonIndex, log10NH, log10PeakFlux,
   Gamma) directly from the optical catalog
   (`OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv` — the
   single canonical source, matching what's used everywhere else in the
   pipeline) rather than letting MICE impute them from population
   correlations.
2. **Projection**: the four prompt parameters (Alpha, Beta, log10Fa,
   log10Ta) are projected from optical- to X-ray-scale using
   `Data/emcee_projection_draws.csv` (4 sampled posterior draws per
   parameter), building 5 candidate feature frames (1 unprojected + 4
   projected) that get pooled for MC uncertainty propagation downstream.
3. **Daume-style domain adaptation**: every predictor gets a scaled duplicate
   (`_opt` suffix, value × is_optical × 0.25) so penalized learners
   (glmnet/ridge) can learn a small, regularized optical-specific correction
   on top of the shared model, rather than needing separate per-domain
   models.
4. **Isolation feature**: a 10th-nearest-neighbor distance in standardized
   predictor space (`FNN::knn.dist`), fit on the training fold only.

This mechanism was found to be dead code until this session (the required
`emcee_projection_draws.csv` never existed) — now live and verified: all 73
optically-recovered GRBs get valid projected values across all 4 draws.

---

## 5. LASSO feature selection

100 repeated `cv.glmnet` fits (`log10(z+1) ~` the 10 base predictors:
log10T90, log10Fa, log10Ta, Alpha, Beta, Gamma, log10Fluence, PhotonIndex,
log10NH, log10PeakFlux), averaging `|coefficient|` across repeats to rank all
10 predictors. All 10 remain as linear terms in the formula search; only the
top 7 by this ranking get a squared term.

---

## 6. Formula generation (100-split search)

Matches the paper's stated methodology (Dainotti et al. 2025, Sec. 4.2.3)
closely, verified directly against the paper text:

1. Build the candidate formula space from the top-7 LASSO variables and
   their squares: O1 (127, linear-only) + O2 (16,383, includes all pairwise
   interactions within `(...)^2` blocks) = 16,510 total — matches the paper's
   stated formula count exactly. A separate 148-formula smoothing-term set
   (`SmoothedO1_formula_list`) is also built but not part of the O2 search
   itself.
2. 100 randomized 80:20 train/test splits. Per split: 10-fold CV on the
   training set for every formula, filter to the subset with r ≥ 99.9th
   percentile AND RMSE ≤ 2nd percentile (falling back to the top-10 by r if
   that intersection is empty), then evaluate the survivors on the held-out
   test set.
3. **Scoring is in LINEAR z-space throughout**, not log10(z+1) — predictions
   are converted via `z = 10^pred - 1` before computing r/RMSE, matching the
   paper's Fig 6 exactly: its own axes are labeled "Linear scale RMSE" /
   "Linear scale correlation," and its stated cutoffs (r=0.565, RMSE=1.167)
   only make sense there (log10(z+1) RMSE never exceeds ~1 for this sample).
   The GAM/GLM models still *fit* log10(z+1) internally (that's the point of
   the log link); only the selection metric changed. (An earlier version of
   this pipeline scored in log10(z+1) space throughout — a reasonable
   substitute at the time, but not what the paper actually does.)
4. Bias (`⟨z_pred − z_obs⟩`, paper Sec 4.4's definition, linear z-space) is a
   **gate, not a vote**: candidates worse than the split's median |bias| are
   dropped, then best-r and best-rmse among the survivors each cast one vote
   (200 total votes across 100 splits). Diagnosed empirically that voting on
   bias directly doesn't work: the best-|bias| criterion picked a
   near-different formula on ~90/100 splits (vs. ~60-70 for r/rmse) — a
   single split's mean signed residual is too noisy a statistic (small
   held-out test sets, ±errors cancel by luck) to vote on directly; it
   rewards lucky cancellation, not real calibration. Demoted to a pre-filter
   instead. (An even earlier version gave bias its own equal-weight vote,
   before this diagnosis.)
5. Tally the r/rmse votes across all 100 splits and take the most-frequent
   formula as the winner.

**Difference from the paper:** the paper runs GAM and GLM as two independent
tracks (6 winning GAM formulas + 4 winning GLM formulas); the main pipeline
here only fits `mgcv::gam()` for scoring, one track (an independent GLM
generator — `stats::glm()` instead of `mgcv::gam()`, otherwise identical —
exists as of this session, but only in the experimental
`emcee_v2_native_opt_formula_generation/` test folder, not yet ported to the
main pipeline). The paper's M-estimator step uses its GLM-track winner
specifically; here it uses whichever (single-track) formula wins the tally.

**Also confirmed:** the single-split "web app" version of this methodology
described in the paper (formulas generated without squared terms, one
90/10 split) is a real, documented fallback in the paper itself — not a bug
— though this implementation's single-split variant does include squared
terms, going slightly further than the paper's literal web-app version.

---

## 7. M-estimator outlier cut

Using the winning formula from step 6, fit `MASS::rlm(..., method="M")` on
the full (post-non-negotiable-cut) sample, then drop the bottom 5% by weight
(quantile-based, not the paper's fixed 0.65 absolute threshold — same
~5%-removal intent, different literal mechanism). Output:
`final_outliers_removed.csv`, ready for SuperLearner.

---

## 8. SuperLearner

Run with `do_m_estimator = FALSE` (the M-estimator cut already happened in
step 7 using the winning formula; re-cutting generically would be a
redundant second cut on a different criterion) and
`use_formula_learners = TRUE` (installs the winning formula as one of the
learners in the ensemble). 80/20 holdout split, 10-fold internal CV,
`loop = 10` external repetitions for Monte Carlo error bars, optional 2σ
catastrophic-outlier retrain pass.

---

## 9. Known issues found and fixed this session

1. **Real-value recovery was dead code** in `superlearner_spencer.R` — the
   projection mechanism (§4, step 2) silently never fired because its input
   file never existed. Fixed by wiring up the export from the combination
   script.
2. **Recovery pulled from the wrong optical source** — was reading the raw,
   uncleaned `optical_data.csv` instead of the single canonical
   `(1).csv` file used everywhere else. Fixed; confirmed the two sources
   agree on all but one GRB/column pair anyway (170405A, log10PeakFlux), and
   that GRB is excluded from the final sample regardless.
3. **`err>value` criterion unstable for PeakFlux** — see §3. Excluded from
   the non-negotiable cuts.
4. **emcee MCMC calibration was unseeded** — different posterior every run,
   which silently propagated into OT's donor-matching via the shared
   `emcee_chains_v2.npz` file (OT's own logic is deterministic; only the
   shared input varied). Fixed via seeded, pooled multi-chain fitting (§2a).
   Verified reproducible for both combination methods.
5. **`err>value` was too broad** — originally checked on T90/Fa/Ta/Fluence/
   Alpha/Beta/PhotonIndex; restricted to Alpha/Beta/Fa/Ta only (see §3) since
   the others are direct measurements, not emcee-fit parameters.
6. **`spencer_edits` `build_design()` bug**: two independently-computed LASSO
   rankings (the formula-search script's own, and `superlearner_spencer.R`'s
   own) aren't guaranteed to agree on which 7 variables rank top — when they
   disagreed (OT specifically), a needed `*Sqr` design-matrix column was
   never built, crashing every learner fit with "object not found." Fixed by
   unioning in whatever the fixed formula file actually references,
   regardless of the script's own ranking that run.
7. **Non-negotiable outlier list was a stale manual file** — see §3. Now
   self-regenerating (`generate_nonnegotiable_outliers.R`), single source of
   truth.
8. **Formula-selection scoring was in log10(z+1) space, not linear z** —
   see §6. The paper's own Fig 6 scores and thresholds are linear-z; fixed to
   match.
9. **Bias was an equal-weight vote, not a gate** — see §6. Diagnosed as too
   noisy a criterion to vote on directly (picks a near-different formula on
   ~90/100 splits); demoted to a pre-filter.

---

## 10. Final data-count summary

| Stage | emcee v2 (linear-z, corrected) | OT v3 (stale — not yet rerun) | Paper (Dainotti et al. 2025) |
|---|---|---|---|
| Combined dataset | 296 | 296 | — |
| Long-GRB scope filter (T90>2s) | 285 | 285 | — |
| Non-negotiable cuts (full removal) | **281** (−4) | 283 (stale) | — |
| Pre-M-estimator pool | 281 | 283 (stale) | 238 |
| M-estimator 5% cut | **266** (−15) | 268 (stale) | 226 (−12) |
| SuperLearner CV (with cat. outliers) | n=213, r(log)=0.672, r(z)=0.664 | not yet rerun | — |
| SuperLearner CV (without cat. outliers) | **n=204, r(log)=0.720, r(z)=0.700** | not yet rerun | — |
| MC, inside 2σ cone | n=55, r(log)=0.964, r(z)=0.967 | not yet rerun | — |

emcee rerun on the new 266-GRB final data (`runs/emcee_v2_cleaned_linearz_
5pct/`), plain `superlearner.R` (do_m_estimator=FALSE, use_formula_learners=
TRUE, loop=10) with the new linear-z winning formula (§11) installed. A real,
meaningful improvement over the previous (269-GRB, pre-linear-z) result of
r(log)=0.671/r(z)=0.660 — both the corrected non-negotiable cuts and the
linear-z-selected formula appear to genuinely help, not just shuffle the
sample. OT v3 has not been rerun with either the Alpha/Beta/Fa/Ta-only
err>value fix or linear-z scoring — its row reflects the last actual run and
should be treated as stale until rerun.

## 11. Winning formulas (final, on stable/reproducible data)

**emcee v2** (linear-z, corrected non-negotiable cuts; 8/200 votes, r+rmse
tally with bias as gate):
```
Response ~ (log10FaSqr + log10T90 + log10Fa + log10PeakFlux)^2 + log10NH +
    PhotonIndex + log10Ta + Alpha + log10NHSqr + log10PeakFluxSqr +
    PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr
```
Structurally new relative to every earlier version of this pipeline: log10T90
is now inside the core interaction block in place of log10Ta (which drops to
a linear-only term). Full formula list + plots:
`formula_results/best_formulas_emcee_LINEARZ_run.csv` and
`formula_results/plots/rmse_vs_correlation_emcee_LINEARZ_run.png`.

**OT v3** (21/100 splits — stale, from the pre-linear-z, broader-err>value
run; not yet rerun with either fix):
```
Response ~ (log10FaSqr + log10Fa + log10PeakFlux)^2 + log10Ta + PhotonIndex +
    log10NH + Alpha + log10T90 + log10TaSqr + log10PeakFluxSqr +
    PhotonIndexSqr + log10NHSqr + AlphaSqr + log10T90Sqr
```

---

## 12. Experimental (not part of the main pipeline)

`GAM_supercomputer/emcee_v2_native_opt_formula_generation/` — a supercomputer-
bound test, entirely separate from the main pipeline above (nothing in
§1-§11 was touched to build it). Tests two changes at once: (a) 4 extra
candidate columns (`*_native`) carrying the 73 optically-sourced GRBs' real
native optical-band Alpha/Beta/Fa/Ta, additive alongside the standard
emcee-projected columns, not overwriting them; (b) LASSO keeps top-9
variables instead of top-7. Cost ~67 hours at laptop-scale parallelism (16x
the standard formula count), hence supercomputer submission (`MyLaptop = F`,
71 MPI workers) rather than a local run. Includes both a GAM and an
independent GLM search (`stats::glm()`), matching the paper's Sec 4.2.3
"determined independently" framing — the first GLM-track implementation in
this codebase, though confined to this experimental folder for now. See that
folder's `README.md` for full detail and rationale.

Both share the same core interaction block (`log10FaSqr + log10Fa +
log10PeakFlux`, sometimes with `log10Ta` folded in) — consistent with the
paper's own published best formula (Eq. 1), which has an identical structure
except for one variable (the paper's 7th LASSO-selected variable was Gamma;
ours is log10T90, reflecting a slightly different top-7 ranking on this
data).

Fixing the PeakFlux artifact (§3) recovered most of the correlation that the
buggy criterion had cost: emcee went from r(log)=0.615 to **0.671**, OT from
0.561 to **0.643** — both now close to or matching their respective
pre-this-session baselines, while keeping the more defensible full-removal
cleaning, real recovery, and robustness-checked formula search intact.

Output locations: `GAM_supercomputer/{emcee_v2,ot_v3}_..._cleaned_formula_generation/`
for the formula search, `runs/{emcee_v2,ot_v3}_cleaned_5pct/` for the final
SuperLearner results.
