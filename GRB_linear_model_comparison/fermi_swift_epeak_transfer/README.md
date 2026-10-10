# Fermi → Swift transfer learning for GRB redshift (E_peak feasibility)

Feasibility study (with Prof. M. Dainotti): can Fermi-GBM **prompt-emission**
information — chiefly the spectral peak energy **E_peak** — improve the Swift
X-ray redshift estimator, by *transferring* it from Fermi into the Swift sample
using the GRBs observed by both instruments?

## Headline result

**No.** Adding E_peak does not improve redshift prediction, by any of the three
methods tried. The reason is physical: the **observed E_peak is uncorrelated with
redshift** (Pearson r = 0.02 over 163 Fermi GRBs). The useful E_peak relations
(Amati, Yonetoku) are *rest-frame* and require knowing z — circular for a
redshift predictor.

| test | result |
|---|---|
| observed E_peak vs redshift (163 GRBs) | r = 0.02 (p = 0.8) — flat |
| winning GLM + transfer-learned E_peak (281) | 0.518 → 0.515 (Δ = −0.004) |
| random forest + transfer-learned E_peak (281) | 0.413 → 0.412 (Δ = −0.002) |
| pooled Swift+Fermi on shared vars + dummy | 0.287 → 0.275 (Δ = −0.012) |

Full explanation and figures are in **`transfer_learning_feasibility.pdf`**.

## Contents

**Deliverables**
- `transfer_learning_feasibility.pdf` — the write-up (what we're transferring, how
  it's done, why it doesn't help). Read this first.
- `Epeak_vs_redshift_for_Maria.png` — observed E_peak vs z, flat (the "why").
- `redshift_pred_with_without_epeak.png` — true vs predicted z, winning GLM, with
  vs without E_peak — indistinguishable (the "so what").

**Data**
- `common67_transfer.csv` — the 67 bridge GRBs (Swift plateau features + Fermi
  FLNC_COMP spectral fits with uncertainties + redshift). 60 have well-constrained E_peak.
- `epeak_z_163.csv` — 163 Fermi GRBs with measured z and E_peak (the correlation test).
- `fermi_shared_vars.csv` — Fermi shared variables (index, T90, fluence, peak flux)
  for the pooled test; `is_common` flags the 65 that overlap the Swift-67.
- `source_GRB_redshift_table_SWIFTcompare176.xlsx` — Maria's original file (Sheet 4
  defines the common GRBs).
- `maria_brief.txt` — Maria's original instructions.

**Code (reproducible, offline)**
- `transfer_learning_analysis.py` — regenerates every headline number and both
  figures from the CSVs here + `../grb_imputed_281.csv`. No network needed.
- `build_note.py` — rebuilds the PDF.

## Reproduce

```bash
# from this folder, in the project's python env (needs: numpy pandas patsy
# scipy scikit-learn matplotlib reportlab pillow)
python3 transfer_learning_analysis.py     # prints the results, regenerates figures
python3 build_note.py                      # rebuilds the PDF
```

`transfer_learning_analysis.py` reads `../grb_imputed_281.csv` (the 281-burst Swift
redshift sample in the parent folder). The GBM spectra in the CSVs were originally
pulled from the HEASARC `fermigbrst` catalog via its TAP service; the saved CSVs
make the analysis reproducible without re-querying.

## How the transfer was done (short version)

- **Bridge:** 67 GRBs observed by both Swift and Fermi (Sheet 4) link the two instruments.
- **(4a) Direct:** add *measured* E_peak to the redshift model on the 60 good bridge GRBs → no gain.
- **(4b) Transfer map:** learn Swift→E_peak on the bridge, impute E_peak for all 281,
  add it → no gain (imputed E_peak is a function of the Swift features, so it's
  redundant — true for GLM, RF, or SuperLearner).
- **(4c) Pooled:** combine both samples on the shared variables + instrument dummy →
  slightly worse, and capped low because the shared variables exclude Swift's strong
  plateau features.
