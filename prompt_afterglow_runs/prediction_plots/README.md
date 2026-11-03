# Prediction plots: prompt / afterglow / Fermi SuperLearner runs

These are copies of each run's predicted-vs-observed plots. The originals stay in
`runs/<run>/Plot_Output/`.

**`headline/`** has one plot per run: cross-validated predicted vs observed z in
linear z, after the 2σ catastrophic-outlier clip. That's the plot the flagship's
0.707 comes from, so the r in each title is the headline number.

| Plot | Run | r (linear z, after clip) |
|---|---|---|
| `00_flagship_0.707_reference` | Flagship 0.707 (earlier M-est snapshot, emcee formula) | 0.707 |
| `01_prompt_swiftpf` | Prompt only, Swift peak flux | 0.469 |
| `02_prompt_gbmpf` | Prompt only, GBM peak flux | 0.458 |
| `03_afterglow` | Afterglow only | 0.611 |
| `04_all_swiftpf` | All features, Swift peak flux (baseline) | 0.681 |
| `04b_all_swiftpf_emceeformula` | Run 4 data with the emcee formula | 0.679 |
| `05_all_gbmpf` | All features, GBM peak flux | 0.641 |
| `06_revmulti_all_swiftpf` | Reversed multivariate, Swift peak flux | 0.662 |
| `07_revmulti_all_gbmpf` | Reversed multivariate, GBM peak flux | 0.645 |
| `08_all_gbmpf_epeak0` | All + GBM peak flux + Epeak (0 if missing) | 0.612 |
| `09_all_gbmpf_epeakmice` | All + GBM peak flux + Epeak (MICE-imputed) | 0.671 |
| `10_prompt_nh_swiftpf` | Prompt + NH, Swift peak flux | 0.618 |
| `11_prompt_nh_gbmpf2` | Prompt + NH, GBM peak flux (corrected match) | 0.601 |
| `12_prompt_gbmpf2` | Prompt only, GBM peak flux (corrected match) | 0.486 |
| `13_all_gbmpf2` | All features, GBM peak flux (corrected match) | 0.684 |
| `14_all_gbmpf2cal` | All features, Theil-Sen-calibrated GBM peak flux (corrected match) | 0.662 |
| `15_revmulti_all_gbmpf2` | Reversed multivariate, GBM peak flux (corrected match) | 0.685 |
| `16_all_gbmpf2_epeak0` | Run 13 + Epeak (0 if missing) | 0.661 |
| `17_all_gbmpf2_epeakmice` | Run 13 + Epeak (MICE-imputed) | 0.688 |
| `13b_all_gbmpf2_emceeformula` | Run 13 data with the flagship emcee formula | 0.681 |
| `17b_all_gbmpf2_epeakmice_emceeformula` | Run 17 data with the flagship emcee formula | 0.676 |

Runs 2, 5, 7, 8 and 9 used the old date-based GBM crosswalk, which matched about
39 GRBs to the wrong GBM burst. Runs 11-17 use the corrected time+position match.

The per-run folders (`00_…` to `09_…`) each hold the full set of 11 plots for that run:

- `z_pred_v_obs_{linear,log}_{with,without}_catOutl_correlation_plot.png` is predicted vs
  observed z, in linear z or log10(z+1), before ("with") or after ("without") the 2σ clip.
- `For_proposal_*` are the same four plots in the proposal styling.
- `mc_pred_vs_obs_*` are the Monte-Carlo error-propagation versions.

All numbers include the full-sample M-estimator cut and the 2σ clip, so they share
that leakage caveat. The full numbers table is `../results_summary.csv`.
