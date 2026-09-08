# Multimodal GRB Redshift Estimation: Fusion Experiments Summary

**Date:** June 2026 · **Question:** can combining the X-ray and optical afterglow
catalogs improve machine-learned redshift estimates over either modality alone?

**Answer:** no fusion architecture meaningfully beats the single-modality models
on their own samples — the modalities are largely redundant about redshift — but
the project gained something real anyway: optical-only GRBs went from effectively
unpredictable (r ≈ 0.09 under imputation) to r ≈ 0.50 with a dedicated optical
model, and a hybrid system now covers **all 294 GRBs at r = 0.527** versus the
previous state of r ≈ 0.52 on the ~70% of GRBs with X-ray data and nothing on the
rest. The win was coverage, not boost.

![Method comparison](Plot_Output/fusion_method_comparison.png)

---

## 1. Data

| Catalog | File | Rows | After T90 > 2s cut | Features |
|---|---|---|---|---|
| X-ray | `x-ray_data.csv` | 222 | 207 | log10T90, log10Fa, log10Ta, Alpha, Beta, Gamma, log10Fluence, PhotonIndex, log10NH, log10PeakFlux (+ 8 error cols) |
| Optical | `optical_data.txt` | 179 | 171 | log10T90, logFa, logTa, Alpha, beta (+ 4 error cols; redshift-derived cols excluded) |

Cross-matching (exact GRB id, then date-root, e.g. 050904 ↔ 050904A) yields the
working union: **123 X-ray-only, 84 paired, 87 optical-only = 294 unique GRBs.**
Two pairs (GRB081203A, GRB081029) have conflicting catalog redshifts (Δz > 0.05)
and are excluded. Loaders/matching: `shared_latent_model.R`
(`load_xray_modality`, `load_optical_modality`, `match_modalities`).

Target throughout: y = log10(z + 1). Metrics: Pearson r and RMSE of out-of-fold
CV predictions (`regression_metrics` in `shared_latent_model.R`).

**A note on the historical 0.63 baseline.** The often-quoted r ≈ 0.63–0.69 for
the X-ray SuperLearner (`Results10fCVResults_OG_*.csv`, n = 180) comes from the
older ~225-GRB `combined_data_with_redshift` sample. On the current 222-GRB
`x-ray_data.csv` the same training scheme achieves r ≈ 0.50–0.53. Much of the
apparent "combining made it worse" effect was (a) this dataset-vintage gap and
(b) pooled metrics being dragged down by unpredictable optical-only GRBs — not
model degradation (see §3).

---

## 2. Experiments and code paths

### 2.1 X-ray-only SuperLearner (baseline)
**Code:** `superlearner.R` (+ `Custom_SL/*`, `mc_error_propagation.R`)
The production pipeline: MICE imputation (midastouch, m = 20), squared-term
features, tuned GAM/GLM/bayesglm formula learners (`Best_formula_GAM.txt`,
`Best_formula_GLM.txt`) + generic learners, 10-fold CV, MC error propagation.
**Result:** r ≈ 0.63 (old sample) / **0.52–0.53** (current sample, lean library).

### 2.2 Shared latent space, end-to-end (negative result)
**Code:** `shared_latent_model.R`, `shared_latent_cv.R`,
`shared_latent_superlearner_cv.R`, `tests_shared_latent.R`
**Outputs:** `OutputFiles/SharedLatent*/`
Symmetric encoders (affine or tanh-MLP) map each modality into a shared latent
space; a linear head predicts y; an alignment loss pulls paired GRBs' embeddings
together; gradients flow end-to-end. A second script freezes the encoders and
trains a SuperLearner on the latent vectors.
**Result:** affine r ≈ 0.39–0.41, nonlinear 0.15–0.26. Diagnostic signature of
failure: metrics identical across latent dims 1/2/3/5, and the SuperLearner adds
nothing over the linear head (0.392 vs 0.397) — the latent collapses to rank-1
("the prediction itself"), leaving no structure to mine. The bottleneck also
taxes X-ray-only GRBs, which must squeeze 10 informative features through ≤5
dims that optical can also reach.

### 2.3 Feature concatenation + joint MICE
**Code:** `concat_superlearner.R` · **Outputs:** `OutputFiles/ConcatExperiment/`
One table over the union (10 X-ray + 5 optical columns, NA where a modality is
missing, single joint MICE), versus an X-ray-only arm on identical GRB-level
folds (3 reps × 10 folds). Paired bootstrap over GRBs:

| Scope | n | X-ray only (A) | Concat (B) | Δ(A−B), 95% CI |
|---|---|---|---|---|
| X-ray GRBs | 207 | 0.512 | 0.488 | +0.024 [−0.040, +0.088] |
| paired | 84 | 0.568 | 0.588 | −0.020 [−0.114, +0.071] |
| X-ray-only GRBs | 123 | 0.471 | 0.424 | +0.047 [−0.040, +0.136] |
| optical-only GRBs | 87 | — | **0.085** | |

Real optical features help slightly where they exist; MICE-imputed columns hurt
where they don't; optical-only GRBs (all 10 X-ray features imputed) are
essentially unpredictable and crush any pooled metric they enter.

### 2.4 Translation encoder (asymmetric)
**Code:** `translation_model.R`, `translation_cv.R`, `tests_translation.R`
**Outputs:** `OutputFiles/TranslationExperiment/`
Ridge regression maps optical features (+squares) into the 10-dimensional X-ray
feature space, trained on each fold's training pairs with exact-LOO penalty
selection — a supervised mapping with an *observable* target, unlike the
abstract latent. Four arms on identical folds, all preprocessing (MICE, scaling,
translation) fit per fold:

- **A** — X-ray baseline (as 2.1, lean library)
- **O** — optical-features SuperLearner (new baseline, previously missing)
- **T1** — translation as augmentation: pseudo X-ray rows for optical-only GRBs
  feed the full SuperLearner
- **T2** — end-to-end: translation + quadratic ridge head trained jointly by
  Adam; prediction loss backpropagates through the head into the translation
  (gradients finite-difference-verified in `tests_translation.R`)

| Scope | A | O | T1 | T2 |
|---|---|---|---|---|
| X-ray GRBs (207) | 0.528 | — | **0.539** | 0.510 |
| paired (84) | 0.605 | 0.334 | **0.638** | 0.570 |
| optical-only (87) | — | **0.502** | 0.405 | 0.494 |

Key paired-bootstrap contrasts: T1 vs A on X-ray GRBs **+0.011 [−0.011, +0.033]**
(paired subset +0.033 [−0.003, +0.071] — the only fusion effect that approaches
significance); T2 vs O on optical-only GRBs −0.008 [−0.094, +0.081].

**Translation fidelity diagnostic** (`fold_diagnostics.csv`): out-of-sample
feature-space MSE ≈ 0.62–0.92 (scaled units; 1.0 = predicting the mean), LOO-chosen
λ ≈ 40–63. Optical parameters recover only ~20% of X-ray feature variance —
the mapping itself is weak, which bounds what any translation-based method can do.

### 2.5 Late fusion (prediction averaging)
Computed from `OutputFiles/TranslationExperiment/cv_predictions.csv`: blends
w·A + (1−w)·O on the 84 paired GRBs. Best fixed blend (w ≈ 0.75–0.8) gives
r = 0.616 vs A's 0.605: **+0.011 [−0.031, +0.061]**. Blending O on top of T1:
+0.005. Nothing there either.

---

## 3. Synthesis

Every fusion strategy, from naive to engineered, lands within noise of the
single-modality baselines:

| Strategy | Effect vs best single-modality baseline |
|---|---|
| Shared latent (frozen or end-to-end) | −0.10 or worse |
| Feature concat + joint MICE | −0.02 (X-ray GRBs) |
| Frozen translation, optical-only GRBs | −0.10 vs direct optical |
| End-to-end translation (T2) | ≈ 0 |
| Pseudo-row augmentation (T1) | +0.01 to +0.03, n.s. |
| Late fusion | +0.01, n.s. |

The consistent verdict across five architectures is the signature of
**redundant information**, not of optimization or architecture failure: whatever
the optical afterglow parameters know about z, the X-ray parameters essentially
already know. Physically plausible — both feature sets are dominated by plateau
parameters (Fa, Ta) tracing the same emission mechanism, and their shared
Dainotti-relation structure is the z-sensitive part. The end-to-end experiments
(2.2, 2.4-T2) specifically rule out "the mapping was fit separately from the
predictor" as the explanation: full gradient flow recovered the frozen
translation's losses but never exceeded the direct single-modality model.

Two earlier failure modes are now understood and avoidable:
- **Composition artifacts:** pooled metrics over a union sample mix predictable
  and unpredictable GRBs; the historical "combined got much lower" numbers were
  dominated by this. Always report metrics scoped by data availability.
- **Imputation contamination:** rows whose features are mostly imputed (MICE
  across a missing modality) degrade training and are themselves unpredictable.
  Don't force optical-only GRBs through X-ray feature space.

### Recommended configuration (hybrid)

GRBs with X-ray data → X-ray SuperLearner (optionally T1's augmented variant);
optical-only GRBs → the optical-feature SuperLearner (arm O).
**CV performance: r = 0.527, RMSE(log10(z+1)) = 0.162 over all 294 GRBs.**

![Hybrid predictions](Plot_Output/hybrid_predictions_scatter.png)

### Caveats and open directions

- **Power:** with 84 pairs, a true fusion gain of +0.03 is at the edge of
  detectability. T1's paired-subset bump may be real; ~3× more pairs would
  settle it. Newer afterglow catalogs may already make this possible.
- **Feature parametrization:** redundancy is a statement about *these* features.
  Optical information orthogonal to the X-ray parametrization — light-curve
  morphology, colors, host extinction — is not in the catalog and could carry
  independent z signal.
- **Per-feature translation fidelity** (complete-case, n = 90 pairs, LOO ridge):
  log10T90 R² = 0.76 (a near-tautology — the same prompt-duration measurement in
  both catalogs), log10Fluence 0.21, all other X-ray features ≈ 0 (Alpha 0.09,
  log10Ta 0.05; Fa/PeakFlux/Beta/Gamma/NH/PhotonIndex at or below noise).
  Excluding T90, optical predicts essentially nothing about the X-ray afterglow
  parametrization. Note the in-pipeline fold diagnostic (~0.2 mean) is inflated
  by MICE-imputed values being intrinsically smooth; the honest figure is 5–10%.
  Together with the 0.83 residual correlation this pins down the mechanism: the
  two feature sets share exactly two common causes — the redshift imprint
  (signal) and intrinsic burst energetics (shared scatter) — and little else, so
  optical re-measures what X-ray already has and adds band-specific detail that
  is irrelevant to z.

---

## 4. Reproduction

```bash
# Tests (gradient checks, parsers, matching):
Rscript tests_shared_latent.R
Rscript tests_translation.R

# Experiments (~15–45 min each; SMOKE_TEST=true for a structural check):
CONCAT_CORES=4 Rscript concat_superlearner.R 3 10     # feature concat vs x-ray-only
TRANS_CORES=4  Rscript translation_cv.R 3 10          # A / O / T1 / T2

# Knobs: LIB_MODE=full (concat: full superlearner.R library),
#        TRANS_WEIGHT (pseudo-row weight), MATCH_WEIGHT (T2 feature-match loss).
```

Bootstrap CIs in this document were computed from the `cv_predictions.csv`
files (per-GRB predictions averaged over reps, 4000 resamples, paired by GRB).
