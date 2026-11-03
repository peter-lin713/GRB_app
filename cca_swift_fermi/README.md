# Swift ↔ Fermi-GBM shared latent space (CCA family)

Item 3 of Dr. Dainotti's note asks whether a shared latent representation, such as CCA
built from the bursts both instruments observed, aligns Swift and Fermi more
principledly than simple pooling, and whether that helps predict redshift.

## Data
- **Swift view** (`Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv`,
  212 long GRBs): BAT prompt (log T90, log fluence, photon index, log peak flux) plus the
  XRT plateau (log Fa, log Ta, α, β) and log NH.
- **GBM view** (`data/`, copied from `spencer_most_recent` @ `d79a6f8`,
  `GRB_linear_model_comparison/fermi_swift_epeak_transfer/`). These are cutoff-power-law
  ("comp") fits over the fluence interval: log T90, log energy fluence, index,
  log photon flux, and log Epeak. These are the same columns `fermi_shared_vars.csv` was
  built from (checked value for value). Epeak uses Spencer's quality filter.
- **Bridge:** 64 long GRBs seen by both instruments (59 with a reliable Epeak).
  **Fermi-only:** 106 long GRBs with z (103 with Epeak). **Optical plateau:** 78 of the
  Swift GRBs; 26 GRBs have all three instruments.
- Caveat: the GBM "photon flux" is averaged over the whole burst, not the GBM 1-s peak
  flux, which isn't in any file we have.

## Methods (`cca_lib.py`, numpy only)
All methods are fit inside training splits only.
- Classical CCA
- Ridge CCA
- PLS-SVD
- Sparse CCA (Witten et al. 2009 PMD)
- Robust CCA (MCD covariance)
- Rank CCA (normal scores + ridge)
- RBF kernel CCA
- Probabilistic CCA / generalized multi-view PCCA, fitted by EM with whole views
  allowed to be missing, so Fermi-only and Swift-only bursts are both usable

Sanity checks:
- CCA matches sklearn exactly.
- Unpenalized sparse CCA reduces to PLS.
- 2-view PCCA reproduces the classical canonical correlations.
- PCCA's EM runs to convergence (tolerance 1e-7). k=1 needs ~1,500 iterations.

Deep CCA was not run. It needs thousands of samples, and torch isn't installed.

## Part A: is there real shared structure? (`part_a_shared_structure.py`)
Held-out canonical correlation from repeated 5-fold CV (10 repeats, hyperparameters tuned
by inner CV), checked against a permutation null (GBM rows shuffled, whole procedure
rerun, 100 shuffles). Full table: `results_part_a.csv`.

| View pair | N | held-out r, comp 1 (range over methods) | null 95th pct | p |
|---|---|---|---|---|
| Swift prompt ↔ GBM prompt | 64 | 0.67 – 0.79 | ≤ 0.27 | 0.01 (floor) |
| All Swift ↔ GBM + Epeak | 59 | 0.68 – 0.81 | ≤ 0.31 | 0.01 (floor) |
| Plateau + NH ↔ GBM + Epeak | 59 | 0.52 – 0.70 | ≤ 0.32 | 0.01 (floor) |

- All methods find strong, real shared structure. The variants agree within ~0.1. Rank
  CCA and robust CCA are usually best on component 1, and plain CCA is worst (it overfits
  the most).
- The shared directions are the obvious ones. In the Swift-all ↔ GBM pair,
  comp 1 = brightness (Swift fluence / peak flux / plateau Fa ↔ GBM fluence / photon flux)
  and comp 2 = duration (T90 ↔ T90).
- Plateau ↔ GBM: comp 1 is X-ray plateau flux Fa ↔ GBM fluence and flux (brighter prompt,
  brighter plateau). Comp 2 (Ta, β ↔ Epeak, index) looks suggestive in-sample (0.47) but
  is weak held out (0.03–0.42 depending on method).
- **Epeak cannot be predicted from Swift.** Held-out R² for log Epeak from the Swift
  features is ≤ 0.06 for every method (OLS, all CCA variants, PCCA). Epeak doesn't load on
  any shared direction, so it carries information Swift lacks. But Spencer showed it is
  uncorrelated with z (r = 0.02).

## Part B: does it help redshift? (`part_b_redshift.py`)
Spencer's honest protocol: 100 random 80/20 splits of the 212 Swift GRBs, refit on 80%,
r in linear z on the untouched 20%. Test GRBs are always predicted from the Swift view
alone. "Pooled" runs add the 106 Fermi-only GRBs to training only. Full table:
`results_part_b.csv` (per-split values in `results_part_b_per_split.csv`).

| Model | mean r | Δ vs SIMPLE7 (paired) | beats SIMPLE7 |
|---|---|---|---|
| **OLS, SIMPLE7 (baseline)** | **0.413** | — | — |
| OLS, all 9 Swift | 0.396 | −0.017 | 32% |
| OLS, Swift prompt only | 0.305 | −0.108 | 13% |
| Best CCA prompt latent (rank CCA), pooled | 0.334 | −0.079 | 23% |
| Other CCA variants, Swift-train or pooled | 0.25 – 0.32 | −0.09 to −0.17 | 7–25% |
| Unsupervised PCCA k=3, pooled | 0.366 | −0.047 | 29% |
| Supervised PCCA k=3, pooled | 0.404 | −0.009 | 45% |
| **Generalized PCCA (Swift + optical + GBM + z) k=2, pooled** | **0.417** | **+0.004 ± 0.052** | **56%** |

- Nothing beats plain OLS on Spencer's 7 Swift features. The best shared-latent model,
  4-view generalized PCCA, ties it (+0.004, wins 56% of splits).
- Pooling the 106 Fermi-only GRBs into training changes results by ≤ 0.01. A prompt-only
  shared space caps out around 0.33, because the plateau parameters (which Fermi lacks)
  carry most of the redshift signal. This matches Spencer's pooled result.
- Sparse CCA and kernel CCA are the worst for redshift. Sparse keeps only the
  single strongest variable per component.

## Bottom line
CCA does recover a real, significant shared Swift ↔ GBM space, and every variant agrees on
what it contains: brightness and duration. That shared content is already in Swift's own
features, though. The one GBM quantity Swift doesn't have, Epeak, sits outside the shared
space and doesn't track redshift. So no CCA variant improves honest redshift prediction
over a plain 7-feature linear model. The multi-view probabilistic models are the only
ones that don't lose ground.

## Reproduce
```bash
python3 part_a_shared_structure.py   # ~2 min on 12 cores  -> results_part_a.csv
python3 part_b_redshift.py           # ~1 min on 12 cores  -> results_part_b*.csv
```
