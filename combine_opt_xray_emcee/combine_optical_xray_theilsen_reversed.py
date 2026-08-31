#!/usr/bin/env python3
"""
Combine optical and x-ray GRB data using Theil-Sen calibration -- REVERSED
direction (x-ray -> optical), mirroring combine_optical_xray_emcee_v2.py's
structure exactly but swapping MCMC for single-variable Theil-Sen and
swapping which side gets projected onto which scale.

No committed script for the original (optical -> x-ray, Theil-Sen) 0.707
combination was ever found in this repo's git history -- only its output
CSVs survived (Data/superlearner_training_theilsen*.csv). This script
replicates that method's structure (same overlap-GRB anchor sample, same
4 prompt parameters, same per-parameter outlier cuts, same relative-error
cut on projected values) but fits Optical ~ X-ray (instead of X-ray ~
Optical) via scipy.stats.theilslopes, then projects X-RAY-ONLY GRBs onto
OPTICAL scale (instead of projecting optical-only GRBs onto X-ray scale).
Final dataset carries all 4 prompt parameters on the OPTICAL scale, with
real optical measurements preferred over projected ones for overlap GRBs.

Uncertainty propagation: Theil-Sen has no natural posterior, so a bootstrap
of the calibration sample (1000 resamples, each refit with theilslopes)
stands in for the MCMC posterior -- same posterior-predictive-style
perturb-and-propagate logic as the emcee script, just with bootstrap draws
instead of MCMC samples.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import theilslopes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
DATA_DIR   = os.path.join(REPO_DIR, 'Data')
OUT_DIR    = REPO_DIR
OPT_FILE   = os.path.join(DATA_DIR, 'OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv')
XRAY_FILE  = os.path.join(DATA_DIR, 'Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv')
OUT_FILE   = os.path.join(DATA_DIR, 'superlearner_training_theilsen_reversed_xray_to_opt.csv')
PLOT_DIR   = os.path.join(OUT_DIR, 'Plot_Output', 'theilsen_calibration_reversed')

N_BOOT  = 1000
MIN_ERR = 1e-6

# Load optical data

opt = pd.read_csv(OPT_FILE)
opt_sel = opt[['GRB', 'z', 'log10T90',
               'log10Faopt', 'log10FaErr_opt', 'log10Taopt', 'logTaErr_opt',
               'Alpha_opt', 'AlphaErr_opt', 'Beta_opt', 'betaErr_opt',
               'Gamma', 'log10Fluence', 'log10FluenceErr',
               'PhotonIndex', 'PhotonIndexErr',
               'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']].copy()
opt_sel = opt_sel.rename(columns={
    'log10Faopt':   'logFa_opt',  'log10FaErr_opt': 'logFaErr_opt',
    'log10Taopt':   'logTa_opt',  'logTaErr_opt':   'logTaErr_opt',
    'AlphaErr_opt': 'AlphaErr_opt',
    'betaErr_opt':  'BetaErr_opt',
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
xray['GRB_norm'] = xray['GRB'].apply(
    lambda g: g if (g and g[-1].isalpha()) else g + 'A'
)

print(f'X-ray data loaded:    {xray.shape[0]} GRBs')

xray_sel = xray[[
    'GRB_norm', 'z', 'log10T90', 'T90Err', 'Gamma',
    'log10Fluence', 'log10FluenceErr', 'PhotonIndex', 'PhotonIndexErr',
    'log10NH', 'log10PeakFlux', 'log10PeakFluxErr',
    'logFa_x', 'logFaErr_x', 'logTa_x', 'logTaErr_x',
    'Alpha_x', 'AlphaErr_x', 'Beta_x', 'BetaErr_x',
]].rename(columns={'GRB_norm': 'GRB'}).copy()

# Overlapping GRBs

overlap = set(xray_sel['GRB']) & set(opt_sel['GRB'])
print(f'Overlapping GRBs:     {len(overlap)}  '
      f'(exact suffix matches + naming-convention fixes)')

merged = pd.merge(
    xray_sel[['GRB', 'logFa_x', 'logFaErr_x', 'logTa_x', 'logTaErr_x',
              'Alpha_x', 'AlphaErr_x', 'Beta_x', 'BetaErr_x']],
    opt_sel,
    on='GRB',
)
print(f'Merged overlap shape: {merged.shape}')

# REVERSED: fit Optical ~ X-ray (x=x-ray, y=optical), the mirror image of
# the emcee script's Optical ~ X-ray... no: emcee fits X-ray ~ Optical
# (x=optical, y=x-ray). Here we fit Optical ~ X-ray (x=x-ray, y=optical).
params_map = {
    'logFa': ('logFa_x',  'logFaErr_x',  'logFa_opt',  'logFaErr_opt'),
    'logTa': ('logTa_x',  'logTaErr_x',  'logTa_opt',  'logTaErr_opt'),
    'Alpha': ('Alpha_x',  'AlphaErr_x',  'Alpha_opt',  'AlphaErr_opt'),
    'Beta':  ('Beta_x',   'BetaErr_x',   'Beta_opt',   'BetaErr_opt'),
}

# Same per-parameter outlier cuts as the emcee script, applied on whichever
# axis now carries that parameter's values (x-ray axis is now the predictor,
# optical axis is now the response -- values themselves are unchanged).
param_cuts = {
    'logTa': (6.0, 6.0, False),
    'Alpha': (2.5, 2.5, True),
    'Beta':  (1.5, 1.5, True),
}
CALIB_ERR_RATIO = 0.5

def apply_calib_cuts(sub, xc, xe, yc, ye, param):
    if param in param_cuts:
        x_max, y_max, use_err = param_cuts[param]
        before = len(sub)
        if use_err:
            sub = sub[(sub[xc] + sub[xe] < x_max) & (sub[yc] + sub[ye] < y_max)]
            label = 'val+err'
        else:
            sub = sub[(sub[xc] < x_max) & (sub[yc] < y_max)]
            label = 'val'
        if len(sub) < before:
            print(f'  {param}: dropped {before - len(sub)} value-outlier(s) ({label} >= {x_max})')
    before = len(sub)
    sub = sub[
        (sub[xe].abs() / sub[xc].abs().clip(lower=1e-10) <= CALIB_ERR_RATIO) &
        (sub[ye].abs() / sub[yc].abs().clip(lower=1e-10) <= CALIB_ERR_RATIO)
    ]
    if len(sub) < before:
        print(f'  {param}: dropped {before - len(sub)} error-ratio outlier(s) (err/|val| > {CALIB_ERR_RATIO})')
    return sub

theilsen_results = {}
print()

rng_boot = np.random.default_rng(7)
BOOT_SEED_BASE = 2000

for p_idx, (param, (xc, xe, yc, ye)) in enumerate(params_map.items()):
    sub = merged[[xc, xe, yc, ye]].apply(pd.to_numeric, errors='coerce').dropna()
    sub = apply_calib_cuts(sub, xc, xe, yc, ye, param)
    x_vals = sub[xc].values
    y_vals = sub[yc].values
    n = len(x_vals)

    m_med, b_med, m_lo, m_hi = theilslopes(y_vals, x_vals)
    resid = y_vals - (m_med * x_vals + b_med)
    sigma_int = np.median(np.abs(resid)) * 1.4826  # robust residual scatter (MAD-based)

    # Bootstrap the calibration sample to approximate a posterior over (m, b),
    # standing in for the emcee script's MCMC chain.
    rng_p = np.random.default_rng(BOOT_SEED_BASE + p_idx)
    boot_m = np.empty(N_BOOT)
    boot_b = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = rng_p.integers(0, n, size=n)
        bm, bb, _, _ = theilslopes(y_vals[idx], x_vals[idx])
        boot_m[i] = bm
        boot_b[i] = bb

    theilsen_results[param] = dict(
        m=m_med, b=b_med, sigma_int=sigma_int,
        boot_m=boot_m, boot_b=boot_b, n=n,
    )
    print(f'{param}: m={m_med:.4f} (CI {m_lo:.4f}-{m_hi:.4f})  b={b_med:.4f}  '
          f'sigma_int={sigma_int:.4f}  N={n}  ({N_BOOT} bootstrap draws)')

# Calibration diagnostic plot

os.makedirs(PLOT_DIR, exist_ok=True)
axis_labels = {
    'logFa': (r'log$_{10}(F_a)$ X-ray', r'log$_{10}(F_a)$ optical'),
    'logTa': (r'log$_{10}(T_a)$ X-ray', r'log$_{10}(T_a)$ optical'),
    'Alpha': (r'$\alpha$ X-ray',         r'$\alpha$ optical'),
    'Beta':  (r'$\beta$ X-ray',          r'$\beta$ optical'),
}

fig, axes = plt.subplots(2, 2, figsize=(10, 9))
fig.suptitle('Theil-Sen calibration (REVERSED, x-ray -> optical): X-ray vs optical (overlapping GRBs)\n'
             'Line = Theil-Sen median fit; band = bootstrap 1$\\sigma$ envelope',
             fontsize=11)

for ax, (param, (xc, xe, yc, ye)) in zip(axes.flat, params_map.items()):
    sub = merged[[xc, xe, yc, ye]].apply(pd.to_numeric, errors='coerce').dropna()
    sub = apply_calib_cuts(sub, xc, xe, yc, ye, param)
    x_d, y_d, sx, sy = sub[xc].values, sub[yc].values, sub[xe].values, sub[ye].values

    r = theilsen_results[param]
    x_line = np.linspace(x_d.min(), x_d.max(), 200)
    y_line = r['m'] * x_line + r['b']
    y_band = r['boot_m'][:, None] * x_line[None, :] + r['boot_b'][:, None]
    y_lo = np.percentile(y_band, 16, axis=0)
    y_hi = np.percentile(y_band, 84, axis=0)

    ax.errorbar(x_d, y_d, xerr=sx, yerr=sy, fmt='o', ms=3, color='steelblue',
                alpha=0.6, elinewidth=0.7, capsize=2, label=f'Overlap (N={len(x_d)})')
    ax.plot(x_line, x_line, 'k--', lw=1, alpha=0.5, label='1:1')
    ax.plot(x_line, y_line, 'r-', lw=1.5, label=f"m={r['m']:.2f}, b={r['b']:.2f}")
    ax.fill_between(x_line, y_lo, y_hi, color='red', alpha=0.15, label='bootstrap 1$\\sigma$')

    xlab, ylab = axis_labels[param]
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel(ylab, fontsize=9)
    ax.set_title(param, fontsize=10)
    ax.legend(fontsize=7, loc='upper left')
    ax.tick_params(labelsize=8)

plt.tight_layout()
out_path = os.path.join(PLOT_DIR, 'theilsen_xray_vs_optical.png')
plt.savefig(out_path, dpi=150)
plt.close()
print(f'\nCalibration plot saved to {out_path}')

# Transform X-RAY-ONLY GRBs to OPTICAL scale (reversed direction)

xray_tr = xray_sel[['GRB', 'z', 'log10T90', 'Gamma',
                     'log10Fluence', 'log10FluenceErr',
                     'PhotonIndex', 'PhotonIndexErr',
                     'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']].copy()

rng = np.random.default_rng(42)

for param, (xc, xe, _, _) in params_map.items():
    r = theilsen_results[param]
    boot_m, boot_b, sigma_int = r['boot_m'], r['boot_b'], r['sigma_int']

    x_vals = pd.to_numeric(xray_sel[xc], errors='coerce').values
    x_errs = pd.to_numeric(xray_sel[xe], errors='coerce').values
    x_errs = np.where(np.isfinite(x_errs), np.maximum(x_errs, MIN_ERR), MIN_ERR)

    n_boot, n_grb = len(boot_m), len(x_vals)
    pick = rng.integers(0, n_boot, size=n_boot)  # use full bootstrap set as the "posterior"

    x_pert = x_vals[None, :] + rng.standard_normal((n_boot, n_grb)) * x_errs[None, :]
    preds = (boot_m[pick, None] * x_pert + boot_b[pick, None]
             + rng.standard_normal((n_boot, n_grb)) * sigma_int)

    xray_tr[f'{param}_opt']    = np.nanmedian(preds, axis=0)
    lo = np.nanpercentile(preds, 16, axis=0)
    hi = np.nanpercentile(preds, 84, axis=0)
    xray_tr[f'{param}Err_opt'] = (hi - lo) / 2.0

    print(f'{param}: median 1sigma err = {np.nanmedian(xray_tr[f"{param}Err_opt"]):.4f}'
          f'  (sigma_int={sigma_int:.4f})')

xray_tr = xray_tr.rename(columns={
    'logFa_opt': 'logFa_opt_x', 'logFaErr_opt': 'logFaErr_opt_x',
    'logTa_opt': 'logTa_opt_x', 'logTaErr_opt': 'logTaErr_opt_x',
    'Alpha_opt': 'Alpha_opt_x', 'AlphaErr_opt': 'AlphaErr_opt_x',
    'Beta_opt':  'Beta_opt_x',  'BetaErr_opt':  'BetaErr_opt_x',
})

# Outer merge; OPTICAL values now take priority for overlapping GRBs (reversed)

combined = pd.merge(
    xray_tr, opt_sel,
    on='GRB', how='outer', suffixes=('_from_xray', '_from_opt')
)

def resolve(df, opt_col, xray_col):
    oc = df.get(opt_col)
    xc = df.get(xray_col)
    if oc is None:
        return pd.to_numeric(xc, errors='coerce') if xc is not None else np.nan
    if xc is None:
        return pd.to_numeric(oc, errors='coerce')
    return pd.to_numeric(oc, errors='coerce').combine_first(pd.to_numeric(xc, errors='coerce'))

final = pd.DataFrame()
final['GRB'] = combined['GRB']

z_o = combined.get('z_from_opt', combined.get('z'))
z_x = combined.get('z_from_xray', combined.get('z'))
final['Redshift_crosscheck'] = (
    pd.to_numeric(z_o, errors='coerce').combine_first(pd.to_numeric(z_x, errors='coerce'))
)

final['log10T90'] = resolve(combined, 'log10T90_from_opt', 'log10T90_from_xray')

final['log10Fa']    = resolve(combined, 'logFa_opt',    'logFa_opt_x')
final['log10FaErr'] = resolve(combined, 'logFaErr_opt', 'logFaErr_opt_x')
final['log10Ta']    = resolve(combined, 'logTa_opt',    'logTa_opt_x')
final['log10TaErr'] = resolve(combined, 'logTaErr_opt', 'logTaErr_opt_x')
final['Alpha']      = resolve(combined, 'Alpha_opt',    'Alpha_opt_x')
final['AlphaErr']   = resolve(combined, 'AlphaErr_opt', 'AlphaErr_opt_x')
final['Beta']       = resolve(combined, 'Beta_opt',     'Beta_opt_x')
final['BetaErr']    = resolve(combined, 'BetaErr_opt',  'BetaErr_opt_x')

for col in ['Gamma', 'log10Fluence', 'log10FluenceErr', 'PhotonIndex', 'PhotonIndexErr',
            'log10NH', 'log10PeakFlux', 'log10PeakFluxErr']:
    final[col] = resolve(combined, f'{col}_from_opt', f'{col}_from_xray')

t90err_lookup = xray_sel.set_index('GRB')['T90Err']  # only ever tracked on the x-ray side
final['T90Err'] = combined['GRB'].map(t90err_lookup)

col_order = [
    'GRB', 'Redshift_crosscheck', 'log10T90',
    'log10Fa', 'log10Ta', 'Alpha', 'Beta',
    'Gamma', 'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux',
    'T90Err', 'log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr',
    'log10FluenceErr', 'PhotonIndexErr', 'log10PeakFluxErr',
]

out = final[col_order].copy()
out = out[out['Redshift_crosscheck'].notna()].copy()
out = out.set_index('GRB')

# is_optical flag: same data-provenance convention as combine_optical_xray_emcee_v2.py
# (1 if absent from the X-ray catalog) -- a fact about the GRB's source catalogs,
# independent of which calibration direction/scale this script chooses to combine on.
xray_ids_all = set(xray_sel['GRB'])
out['is_optical'] = (~out.index.isin(xray_ids_all)).astype(int)

# The projected-onto-optical-scale population is the X-RAY-only GRBs this time
# (the reverse of the emcee script, which projects optical-only GRBs).
opt_ids = set(opt_sel['GRB'])
xray_only_mask = ~out.index.isin(opt_ids)

RATIO_THRESH = 0.5
xray_only_grbs = out.index[xray_only_mask]
xray_orig = xray_sel[xray_sel['GRB'].isin(xray_only_grbs)].set_index('GRB')
xray_orig = xray_orig.apply(pd.to_numeric, errors='coerce')

fa_cut_orig    = ((xray_orig['logFa_x'].abs() > 0) & (xray_orig['logFaErr_x'] / xray_orig['logFa_x'].abs() > RATIO_THRESH)).reindex(xray_only_grbs, fill_value=False)
ta_cut_orig    = ((xray_orig['logTa_x'].abs() > 0) & (xray_orig['logTaErr_x'] / xray_orig['logTa_x'].abs() > RATIO_THRESH)).reindex(xray_only_grbs, fill_value=False)
alpha_cut_orig = ((xray_orig['Alpha_x'].abs() > 0) & (xray_orig['AlphaErr_x'] / xray_orig['Alpha_x'].abs() > RATIO_THRESH)).reindex(xray_only_grbs, fill_value=False)
beta_cut_orig  = ((xray_orig['Beta_x'].abs()  > 0) & (xray_orig['BetaErr_x']  / xray_orig['Beta_x'].abs()  > RATIO_THRESH)).reindex(xray_only_grbs, fill_value=False)

out_xray = out.loc[xray_only_grbs].apply(pd.to_numeric, errors='coerce')
fa_cut_proj    = (out_xray['log10Fa'].abs() > 0) & (out_xray['log10FaErr'] / out_xray['log10Fa'].abs() > RATIO_THRESH)
ta_cut_proj    = (out_xray['log10Ta'].abs() > 0) & (out_xray['log10TaErr'] / out_xray['log10Ta'].abs() > RATIO_THRESH)
alpha_cut_proj = (out_xray['Alpha'].abs()   > 0) & (out_xray['AlphaErr']   / out_xray['Alpha'].abs()   > RATIO_THRESH)
beta_cut_proj  = (out_xray['Beta'].abs()    > 0) & (out_xray['BetaErr']    / out_xray['Beta'].abs()    > RATIO_THRESH)

fa_cut    = fa_cut_orig    | fa_cut_proj
ta_cut    = ta_cut_orig    | ta_cut_proj
alpha_cut = alpha_cut_orig | alpha_cut_proj
beta_cut  = beta_cut_orig  | beta_cut_proj

any_cut  = fa_cut | ta_cut | alpha_cut | beta_cut
drop_ids = any_cut[any_cut].index
out_cut  = out.drop(index=drop_ids)

cut_file = OUT_FILE.replace('.csv', '_errcut_relative.csv')
out_cut.to_csv(cut_file, index=True)
print(f'\nRatio cut (err/|value| > {RATIO_THRESH}) on original x-ray + projected optical:')
print(f'  Removed {len(drop_ids)} x-ray-only GRBs ({out_cut.shape[0]} remain)  -> {cut_file}')
print(f'  logFa:  xray {fa_cut_orig.sum()}  proj {fa_cut_proj.sum()}  combined {fa_cut.sum()}')
print(f'  logTa:  xray {ta_cut_orig.sum()}  proj {ta_cut_proj.sum()}  combined {ta_cut.sum()}')
print(f'  Alpha:  xray {alpha_cut_orig.sum()}  proj {alpha_cut_proj.sum()}  combined {alpha_cut.sum()}')
print(f'  Beta:   xray {beta_cut_orig.sum()}  proj {beta_cut_proj.sum()}  combined {beta_cut.sum()}')

print(f'\nFinal dataset shape: {out.shape}')
out.to_csv(OUT_FILE, index=True)
print(f'Saved to {OUT_FILE}')

t90_missing_base = out['log10T90'].isna().sum()
t90_missing_cut  = out_cut['log10T90'].isna().sum()
print(f'\n=== T90 completeness check ===')
print(f'  Unfiltered ({out.shape[0]} GRBs):  log10T90 missing = {t90_missing_base}')
print(f'  Error-cut  ({out_cut.shape[0]} GRBs):  log10T90 missing = {t90_missing_cut}')
