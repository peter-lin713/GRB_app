# Generalization analysis — confirmed-exact Theil-Sen flagship model, TOTAL_GENERALIZATION_DATA_v4

Replicates the predecessor paper's (Dainotti et al. 2024c, §3.4.3/§4.4)
generalization methodology, using the model whose exact provenance was
confirmed this session (see `PIPELINE_METHODOLOGY.md` §2c/§11): the
`superlearner_model` in
`runs/theilsen_cleaned_linearz_5pct_mice20_RETRAIN_nesteddata/`, r(z)=0.7071.

## Pipeline (`run_generalization.R`)

1. **Data-quality fix (new finding this session):** `TOTAL_GENERALIZATION_DATA_v4.csv`'s
   `photon_index` and `logPeakFlux` columns contain contaminated string values
   — trailing model-tag suffixes (`"1.81PL"`), stray trailing commas
   (`"2.03,"`, 82 of 299 rows), and one literal unevaluated formula string
   (`"0.43429448190325176*Log[\"n/a\"]"`) — that silently force the whole
   column to character on read, corrupting every downstream numeric
   comparison (confirmed: this alone had inflated the out-of-range filter
   from 21 removed to 85 removed before the fix). Cleaned by extracting the
   leading numeric token; genuinely unparseable values become NA for MICE.
2. **Column remap** to training convention (Fbest→log10Fa, T_abest→log10Ta,
   logT90→log10T90, logPeakFlux→log10PeakFlux, photon_index→PhotonIndex,
   logNH→log10NH; Gamma/Alpha/AlphaErr pass through), following the exact
   precedent set by `Generalization/Generalization_Aditya_v1.R` for the
   earlier v2/v3 generalization sets.
3. **MICE imputation** (m=20, maxit=20, midastouch) for the genuinely
   partially-missing columns (log10T90, log10PeakFlux, PhotonIndex, log10NH,
   Gamma, Alpha — each missing 1-4 of 299 rows).
4. **Beta / log10Fluence — flagged limitation:** the trained model needs
   these two features (used by 2 of the 3 non-zero-weight ensemble members,
   ~38% combined weight — the formula-based GLM learner, ~62% weight, never
   uses them), but they are **100% missing** in this generalization set (0/299
   rows). MICE cannot impute a fully-empty column (no real donor values to
   draw from — confirmed empirically: left all-NA when attempted). Substituted
   the **training-population mean** for both (Beta=0.922, log10Fluence=-5.578)
   as the least-assumption placeholder. This means Beta/Fluence carry zero
   discriminating information for every generalization GRB — a real caveat on
   the ~38%-weighted portion of the ensemble, not a cosmetic one.
5. **Out-of-parameter-space filter** (avoid extrapolation, matching
   `Generalization_Aditya_v1.R`'s established feature set, extended with Alpha
   since v4 — unlike v2/v3 — actually provides real Alpha values): drop any
   GRB outside the training set's [min, max] range on log10NH, PhotonIndex,
   log10Fa, log10Ta, log10PeakFlux, log10T90, Gamma, or Alpha.
   **299 → 278 GRBs (21 removed, 7.0%)**.
6. **Predict**: `predict(sl_model, newdata)` on the confirmed-exact ensemble
   (requires sourcing `Custom_SL/*.R` first, for the formula-based learners'
   S3 `predict` methods).
7. **Bias correction**: `BC_3way` (`Generalization/Bias_Correction_function.R`)
   in its genuine intended mode — bins fit on the training CV set's known
   Zspec (z<2, 2–3.5, z>3.5, matching this paper's own 3-way split), applied
   to the generalization set routed by its *own* predicted Zphot (the only
   thing available for a real unknown-z set).
8. **Distribution validation**: Anderson-Darling (`kSamples::ad.test`,
   matching the predecessor paper's method) and Kolmogorov-Smirnov, both
   generalization-predicted-z vs. training-observed-z, plus per-feature checks.

## Results

| | N | Zphot range | median |
|---|---|---|---|
| Point predictions (raw) | 278 | [0.71, 3.87] | 1.84 |
| Bias-corrected | 278 | [0.39, 5.03] | 1.42 |
| Training Zspec (reference) | 206 | [0.08, 8.26] | 1.76 |

**Distribution tests** (generalization predicted-z vs. training observed-z):

| Test | Raw | Bias-corrected |
|---|---|---|
| Anderson-Darling | AD=22.4, p=1.5e-12 | AD=8.4, p=6.6e-5 |
| KS | D=0.259, p=2.6e-7 | D=0.187, p=5.1e-4 |

**Both REJECTED at p<0.05** — same qualitative outcome as the predecessor
paper (Dainotti et al. 2024c also rejected this test). Bias correction
produces a large, genuine improvement (AD statistic drops ~2.7x) but doesn't
fully resolve the mismatch: the raw predictions are visibly compressed into a
narrow z~1-3 band relative to training's fuller spread (see
`generalization_z_distributions.png`), and correction re-spreads the
distribution but doesn't perfectly recover the training shape.

**Per-feature distribution checks** (KS test, generalization vs. training):
`log10NH` (p=0.0004) and `Alpha` (p=0.0003) differ significantly;
`log10T90` is borderline (p=0.056); `log10Fa`, `log10Ta`, `PhotonIndex` do
not differ significantly. The predecessor paper's own diagnostic similarly
found `log10T90` and `log10NH` failing — a partial, consistent replication of
their root-cause pattern (their sample, cuts, and exact feature set differ
from ours, so an exact match isn't expected).

## Not done (scope limits, flagged rather than silently skipped)

- **No per-point MC uncertainty propagation**: would require per-GRB
  measurement errors for every predictor, including `log10NH`, which this
  generalization set does not provide at all. Only point estimates are
  reported here.
- **No final "prediction error exceeds predicted redshift" cut**: the
  predecessor paper's last reporting step needs a per-point error estimate,
  which requires the MC propagation above. Not applied.
- **Beta/Fluence limitation** (see step 4): real information is missing for
  ~38% of the ensemble's weight; treat this generalization result as
  informative about the model's redshift-range behavior, not as a
  fully-independent validation on features the model was never truly tested
  against for this set.
