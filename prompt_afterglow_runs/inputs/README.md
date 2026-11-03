# Next-steps experiments (Swift BAT E_peak, Fermi GBM peak flux)

Hi Peter — this folder holds two of Dr. Dainotti's suggested experiments plus one
follow-up. Each adds a new quantity to the redshift estimator and checks whether it
improves the honest cross-validated r (Pearson r between predicted and true z, in
linear redshift, repeated 5-fold CV). The estimator is the winning GLM; for these
tests it is fit with ridge regression (stable on the smaller matched samples).

The PDFs `next_steps_figures.pdf` (figures only) and `next_steps_explained.pdf`
(figures + full method and results) summarize everything below.

## Task 1 — Swift BAT prompt E_peak
`bat_retrieve_match.py` -> `task1_bat_epeak.py`

Pull an E_peak directly from the Swift BAT catalogs (HEASARC `swiftgrbba` for the
1-second peak photon index / peak flux, `swiftgrb` for a direct `bat_epeak`), matched
to our 281-GRB redshift sample (274/281 matched; 255 with a 1-s peak photon index,
57 with a direct E_peak). Where no direct E_peak exists we estimate it from the BAT
photon index via a relation calibrated on BAT bursts that have both (never uses z).

- Directly-measured BAT E_peak vs z: r = +0.08 (not significant).
- Adding BAT E_peak / peak flux / photon index to the estimator (baseline r = 0.47,
  255 GRBs): Δr between −0.0001 and −0.014.
- **No improvement.** A clean, non-circular confirmation that observer-frame prompt
  E_peak carries essentially no redshift signal.

Figure: `task1_bat_epeak.png`. Data: `bat_match_281_v2.csv`, `task1_bat_sample.csv`.

## Task 2 — substitute the Fermi GBM peak flux for the Swift one
`gbm_crosswalk.py` -> `task2_gbm_substitution.py`

Train the estimator twice on the same GRBs — once with the Swift peak flux, once with
that one column swapped for the GBM peak flux — and compare CV r. To get past the 65
GRBs in the curated common list, a date + sky-position crosswalk between the 281-GRB
sample and the Fermi GBM catalog (`fermigbrst`) enlarged the common sample to 122
(116 with a usable GBM peak flux).

- Original 65 GRBs: Δr = +0.05 ± 0.11 (73% of resamples positive) — within the noise.
- Enlarged 116 GRBs (healthier baseline r = 0.41): Δr = +0.01 ± 0.09 (53% positive).
- **The +0.05 did not survive the larger sample.** GBM peak flux correlates with z
  (−0.32) about as well as Swift's (−0.30) — no real gain.

Figure: `task2_gbm_substitution.png`. Data: `gbm_crosswalk_281.csv`, `task2_enlarged_sample.csv`.

## Task 2 follow-up — impute GBM into the full 281-GRB sample
`task2b_fullsample_impute.py`

Keep all 281 GRBs; for the subset that has a GBM peak flux, substitute it (cross-
calibrated onto the Swift scale by OLS on the overlap — no z used) and leave everyone
else on Swift. Baseline (all Swift) r = 0.50.

- Impute the 65-set (63 of 281): 0.499 → 0.491, Δr = −0.008 ± 0.024 (36% positive).
- Impute the 116-set (116 of 281): 0.499 → 0.480, Δr = −0.019 ± 0.037 (28% positive).
- **No improvement** (marginally negative, within the noise) — same verdict from the
  other direction.

## Reproduce
```
conda activate xgbenv   # numpy, pandas, scipy, scikit-learn, patsy, astropy, matplotlib, reportlab, pillow
python3 bat_retrieve_match.py       # Task 1 data (HEASARC; network)
python3 task1_bat_epeak.py          # Task 1 result + figure
python3 gbm_crosswalk.py            # Task 2 crosswalk (HEASARC; network)
python3 task2_gbm_substitution.py   # Task 2 result + figure
python3 task2b_fullsample_impute.py # Task 2 follow-up
python3 build_task_pdfs.py          # both PDFs
```
