# Multivariate Theil-Sen, X-ray → optical, under honest linear models

Rebuilds the multivariate Theil-Sen combination in the reversed direction and
scores it with plain linear models. There's no SuperLearner, M-estimator cut,
2σ clip, or full-sample LASSO, because those leak redshift into the evaluation.

## Pipeline
1. `../combine_opt_xray_emcee/combine_optical_xray_theilsen_multivariate_reversed.py`
   fits each optical plateau parameter (Fa, Ta, α, β) from **all four** X-ray plateau
   parameters with `sklearn.TheilSenRegressor` on the 69–73 overlap GRBs, then projects
   X-ray-only GRBs onto the optical scale with bootstrap (1000) uncertainty.
   Everything else (inputs, cuts, error propagation, output schema) is the same as the
   single-variable `combine_optical_xray_theilsen_reversed.py`.
   It writes `Data/superlearner_training_theilsen_multivariate_reversed_xray_to_opt.csv`
   (306 GRBs).
2. `prep_frame.R` runs only the data-prep stage of `GAM_analysis_8variables.R`: the
   long-GRB filter, non-negotiable cuts, and MICE. It writes `frame/grb_xray_imputed.csv`
   (278 GRBs). It reproduces the existing reversed frame exactly (max difference 0).
3. `eval_linear.py` runs the honest protocol from Spencer's `simple_linear_flagship_check.py`:
   100 random 80/20 splits, OLS refit on each, scored on the held-out 20%, and r reported
   in linear z. All four Theil-Sen schemes are scored on their own frames and on the
   275 GRBs common to all four. It writes `results_honest_linear.csv`.

## Result (SIMPLE7, r in linear z, mean over 100 splits)

| scheme | own frame | common 275 |
|---|---|---|
| forward single (0.707 family) | 0.466 | 0.475 |
| forward multivariate | 0.464 | 0.476 |
| reversed single | 0.455 | 0.455 |
| **reversed multivariate (new)** | **0.465** | **0.451** |

All four are the same within noise: split-to-split SD is about 0.11, and on the
same splits the new version differs from the reversed single version by −0.004 ± 0.028.
Neither the calibration direction nor going multivariate changes honest performance.
