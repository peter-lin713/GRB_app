#!/usr/bin/env python3
"""
Combine optical and x-ray GRB data using a MULTIVARIATE Theil-Sen calibration,
REVERSED direction (x-ray -> optical).

Same structure as combine_optical_xray_theilsen_reversed.py (same input
catalogs, same overlap-GRB anchor sample, same per-parameter outlier cuts,
same bootstrap uncertainty propagation, same relative-error cut on projected
values, same output schema). The only change is the calibration model:

  single-variable (reversed script):  Opt_p ~ X_p                  (one x-ray param)
  multivariate    (this script):      Opt_p ~ Fa_x + Ta_x + Alpha_x + Beta_x

i.e. each optical-scale plateau parameter is fit from ALL FOUR x-ray plateau
parameters jointly via sklearn.linear_model.TheilSenRegressor -- the reverse
of the forward multivariate run ("all 4 optical params -> each X-ray-scale
plateau param", combination_scheme_comparison/README.md), whose own build
script was never committed.

Anchor sample per target parameter: overlap GRBs where all four x-ray
predictors AND the optical target pass the same value / error-ratio cuts the
single-variable script applies to each parameter.

Uncertainty: N_BOOT bootstrap resamples of the anchors, each refit with
TheilSenRegressor, stand in for a posterior; each x-ray-only GRB's four x-ray
inputs are perturbed within their own errors, projected through each
bootstrap model plus a MAD-based intrinsic-scatter draw, and the median /
half 16-84 range become the value / 1-sigma error.
"""

import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import TheilSenRegressor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
DATA_DIR   = os.path.join(REPO_DIR, 'Data')
OPT_FILE   = os.path.join(DATA_DIR, 'OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv')
XRAY_FILE  = os.path.join(DATA_DIR, 'Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv')
OUT_FILE   = os.path.join(DATA_DIR, 'superlearner_training_theilsen_multivariate_reversed_xray_to_opt.csv')
PLOT_DIR   = os.path.join(REPO_DIR, 'Plot_Output', 'theilsen_calibration_multivariate_reversed')

N_BOOT  = 1000
MIN_ERR = 1e-6
SEED    = 2000

# Load optical data (identical to the single-variable reversed script)

opt = pd.read_csv(OPT_FILE)
opt_sel = opt[['GRB', 'z', 'log10T90',
               'log10Faopt', 'log10FaErr_opt', 'log10Taopt', 'logTaErr_opt',
               'Alpha_opt', 'AlphaErr_opt', 'Beta_opt', 'betaErr_opt',
               'Gamma', 'log10Fluence', 'log10FluenceErr',
               'PhotonIndex', 'PhotonIndexErr',
               'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']].copy()
opt_sel = opt_sel.rename(columns={
    'log10Faopt': 'logFa_opt', 'log10FaErr_opt': 'logFaErr_opt',
    'log10Taopt': 'logTa_opt', 'betaErr_opt':    'BetaErr_opt',
})
for col in opt_sel.columns:
    if col != 'GRB':
        opt_sel[col] = pd.to_numeric(opt_sel[col], errors='coerce')

print(f'Optical data loaded:  {opt_sel.shape[0]} GRBs')

# Load x-ray data

xray = pd.read_csv(XRAY_FILE)
xray = xray.rename(columns={
    'GRB_Name':            'GRB',
    'Redshift_crosscheck': 'z',
    'log10Fa_X':    'logFa_x',   'log10FaErr_X': 'logFaErr_x',
    'log10Ta_X':    'logTa_x',   'log10TaErr':   'logTaErr_x',
    'Alpha_X':      'Alpha_x',   'AlphaErr':     'AlphaErr_x',
    'Beta_X':       'Beta_x',    'BetaErr_X':    'BetaErr_x',
    'logFluenceErr': 'log10FluenceErr',
})
xray['GRB'] = xray['GRB'].str.replace('GRB', '', regex=False)
xray['GRB'] = xray['GRB'].apply(lambda g: g if (g and g[-1].isalpha()) else g + 'A')

print(f'X-ray data loaded:    {xray.shape[0]} GRBs')

xray_sel = xray[[
    'GRB', 'z', 'log10T90', 'T90Err', 'Gamma',
    'log10Fluence', 'log10FluenceErr', 'PhotonIndex', 'PhotonIndexErr',
    'log10NH', 'log10PeakFlux', 'log10PeakFluxErr',
    'logFa_x', 'logFaErr_x', 'logTa_x', 'logTaErr_x',
    'Alpha_x', 'AlphaErr_x', 'Beta_x', 'BetaErr_x',
]].copy()

PARAMS = ['logFa', 'logTa', 'Alpha', 'Beta']
X_COLS = [f'{p}_x' for p in PARAMS]
X_ERRS = [f'{p}Err_x' for p in PARAMS]

merged = pd.merge(xray_sel[['GRB'] + X_COLS + X_ERRS], opt_sel, on='GRB')
print(f'Overlapping GRBs:     {len(merged)}')

# Same per-parameter cuts as the single-variable script, applied to every
# variable that enters a fit (all four x-ray predictors + the optical target).
param_cuts = {
    'logTa': (6.0, False),
    'Alpha': (2.5, True),
    'Beta':  (1.5, True),
}
CALIB_ERR_RATIO = 0.5

def passes(df, param, vc, ec):
    ok = (df[ec].abs() / df[vc].abs().clip(lower=1e-10)) <= CALIB_ERR_RATIO
    if param in param_cuts:
        vmax, use_err = param_cuts[param]
        ok &= (df[vc] + df[ec] < vmax) if use_err else (df[vc] < vmax)
    return ok

x_ok = np.ones(len(merged), bool)
for p, xc, xe in zip(PARAMS, X_COLS, X_ERRS):
    x_ok &= passes(merged, p, xc, xe).values


def fit_ts(X, y, seed):
    m = TheilSenRegressor(random_state=seed).fit(X, y)
    return np.r_[m.intercept_, m.coef_]


results = {}
print()
for p_idx, p in enumerate(PARAMS):
    yc, ye = f'{p}_opt', f'{p}Err_opt'
    sub = merged[X_COLS + X_ERRS + [yc, ye]].apply(pd.to_numeric, errors='coerce')
    keep = x_ok & passes(sub, p, yc, ye).values & sub.notna().all(axis=1).values
    sub = sub[keep]
    X, y = sub[X_COLS].values, sub[yc].values
    n = len(y)

    coef = fit_ts(X, y, SEED + p_idx)
    resid = y - (coef[0] + X @ coef[1:])
    sigma_int = np.median(np.abs(resid)) * 1.4826  # robust scatter, as in the single-var script

    rng_p = np.random.default_rng(SEED + p_idx)
    idxs = [rng_p.integers(0, n, size=n) for _ in range(N_BOOT)]
    boot = np.array(Parallel(n_jobs=-1)(
        delayed(fit_ts)(X[ix], y[ix], SEED + 10 * p_idx + i) for i, ix in enumerate(idxs)))

    # Single-variable reference fit Opt_p ~ X_p on the same anchors, so the
    # multivariate gain over the plain reversed calibration is visible.
    j = PARAMS.index(p)
    coef_1 = fit_ts(X[:, [j]], y, SEED + p_idx)
    sig_1 = np.median(np.abs(y - (coef_1[0] + X[:, j] * coef_1[1]))) * 1.4826

    results[p] = dict(coef=coef, boot=boot, sigma_int=sigma_int, n=n,
                      fitted=coef[0] + X @ coef[1:], y=y)
    terms = '  '.join(f'{c}={v:+.3f}' for c, v in zip(['b'] + X_COLS, coef))
    print(f'{p}_opt ~ 4 x-ray params:  N={n}  {terms}  sigma_int={sigma_int:.4f}'
          f'  (single-var sigma_int={sig_1:.4f})')

# Calibration diagnostic: fitted (from all 4 x-ray params) vs actual optical

os.makedirs(PLOT_DIR, exist_ok=True)
fig, axes = plt.subplots(2, 2, figsize=(10, 9))
fig.suptitle('Multivariate Theil-Sen calibration (REVERSED, x-ray -> optical)\n'
             'Each optical parameter fit from all four x-ray plateau parameters (overlap GRBs)',
             fontsize=11)
for ax, p in zip(axes.flat, PARAMS):
    r = results[p]
    ax.scatter(r['fitted'], r['y'], s=12, color='steelblue', alpha=0.7,
               label=f"Overlap (N={r['n']})")
    lo, hi = min(r['fitted'].min(), r['y'].min()), max(r['fitted'].max(), r['y'].max())
    ax.plot([lo, hi], [lo, hi], 'k--', lw=1, alpha=0.5, label='1:1')
    rr = np.corrcoef(r['fitted'], r['y'])[0, 1]
    ax.set_title(f"{p}  (r={rr:.2f}, sigma_int={r['sigma_int']:.2f})", fontsize=10)
    ax.set_xlabel(f'{p} optical, Theil-Sen fitted from x-ray', fontsize=9)
    ax.set_ylabel(f'{p} optical, measured', fontsize=9)
    ax.legend(fontsize=7, loc='upper left')
plt.tight_layout()
out_path = os.path.join(PLOT_DIR, 'theilsen_multivariate_xray_vs_optical.png')
plt.savefig(out_path, dpi=150)
plt.close()
print(f'\nCalibration plot saved to {out_path}')

# Project every x-ray GRB onto the optical scale

xray_tr = xray_sel[['GRB', 'z', 'log10T90', 'Gamma',
                    'log10Fluence', 'log10FluenceErr',
                    'PhotonIndex', 'PhotonIndexErr',
                    'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']].copy()

Xv = xray_sel[X_COLS].apply(pd.to_numeric, errors='coerce').values
Xe = xray_sel[X_ERRS].apply(pd.to_numeric, errors='coerce').values
Xe = np.where(np.isfinite(Xe), np.maximum(Xe, MIN_ERR), MIN_ERR)
n_grb = len(Xv)

rng = np.random.default_rng(42)
for p in PARAMS:
    r = results[p]
    boot = r['boot']
    X_pert = Xv[None, :, :] + rng.standard_normal((N_BOOT, n_grb, 4)) * Xe[None, :, :]
    preds = (boot[:, None, 0] + np.einsum('bgk,bk->bg', X_pert, boot[:, 1:])
             + rng.standard_normal((N_BOOT, n_grb)) * r['sigma_int'])
    xray_tr[f'{p}_opt_x'] = np.nanmedian(preds, axis=0)
    xray_tr[f'{p}Err_opt_x'] = (np.nanpercentile(preds, 84, axis=0)
                                - np.nanpercentile(preds, 16, axis=0)) / 2.0
    print(f'{p}: median 1sigma err = {np.nanmedian(xray_tr[f"{p}Err_opt_x"]):.4f}'
          f'  (sigma_int={r["sigma_int"]:.4f})')

# Outer merge; real OPTICAL values take priority for overlap GRBs

combined = pd.merge(xray_tr, opt_sel, on='GRB', how='outer',
                    suffixes=('_from_xray', '_from_opt'))

def resolve(df, opt_col, xray_col):
    return (pd.to_numeric(df[opt_col], errors='coerce')
            .combine_first(pd.to_numeric(df[xray_col], errors='coerce')))

final = pd.DataFrame({'GRB': combined['GRB']})
final['Redshift_crosscheck'] = resolve(combined, 'z_from_opt', 'z_from_xray')
final['log10T90'] = resolve(combined, 'log10T90_from_opt', 'log10T90_from_xray')
for p, out_v in zip(PARAMS, ['log10Fa', 'log10Ta', 'Alpha', 'Beta']):
    final[out_v]         = resolve(combined, f'{p}_opt',    f'{p}_opt_x')
    final[out_v + 'Err'] = resolve(combined, f'{p}Err_opt', f'{p}Err_opt_x')
for col in ['Gamma', 'log10Fluence', 'log10FluenceErr', 'PhotonIndex', 'PhotonIndexErr',
            'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']:
    final[col] = resolve(combined, f'{col}_from_opt', f'{col}_from_xray')
final['T90Err'] = combined['GRB'].map(xray_sel.set_index('GRB')['T90Err'])

col_order = [
    'GRB', 'Redshift_crosscheck', 'log10T90',
    'log10Fa', 'log10Ta', 'Alpha', 'Beta',
    'Gamma', 'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux',
    'T90Err', 'log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr',
    'log10FluenceErr', 'PhotonIndexErr', 'log10PeakFluxErr',
]
out = final[col_order]
out = out[out['Redshift_crosscheck'].notna()].set_index('GRB')
out['is_optical'] = (~out.index.isin(set(xray_sel['GRB']))).astype(int)

# Relative-error cut on x-ray-only (projected) GRBs: original x-ray value OR
# projected optical value with err/|value| > 0.5 -- same as the reversed script.

RATIO_THRESH = 0.5
xray_only = out.index[~out.index.isin(set(opt_sel['GRB']))]
xo = xray_sel.set_index('GRB').loc[xray_only].apply(pd.to_numeric, errors='coerce')
po = out.loc[xray_only].apply(pd.to_numeric, errors='coerce')

def ratio_bad(v, e):
    return (v.abs() > 0) & (e / v.abs() > RATIO_THRESH)

any_cut = pd.Series(False, index=xray_only)
for xc, xe, oc in zip(X_COLS, X_ERRS, ['log10Fa', 'log10Ta', 'Alpha', 'Beta']):
    c_orig, c_proj = ratio_bad(xo[xc], xo[xe]), ratio_bad(po[oc], po[oc + 'Err'])
    any_cut |= c_orig | c_proj
    print(f'  {oc:8s} cut: xray {c_orig.sum():3d}  proj {c_proj.sum():3d}  combined {(c_orig | c_proj).sum():3d}')
out_cut = out.drop(index=any_cut[any_cut].index)

cut_file = OUT_FILE.replace('.csv', '_errcut_relative.csv')
out.to_csv(OUT_FILE)
out_cut.to_csv(cut_file)
print(f'\nRatio cut removed {any_cut.sum()} x-ray-only GRBs')
print(f'Saved {out.shape[0]} GRBs -> {OUT_FILE}')
print(f'Saved {out_cut.shape[0]} GRBs -> {cut_file}')
print(f'is_optical = 1 for {int(out["is_optical"].sum())} GRBs; '
      f'log10T90 missing for {int(out["log10T90"].isna().sum())}')
