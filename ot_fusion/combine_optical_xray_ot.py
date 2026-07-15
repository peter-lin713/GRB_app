#!/usr/bin/env python3
"""
combine_optical_xray_ot.py
OT-fusion of optical and X-ray GRB catalogs with anchor-informed metric learning.
Produces a training CSV in the same column format as
superlearner_training_emcee_errcut_relative.csv so it can be fed directly into
superlearner.R with identical settings.

Key differences from the emcee MCMC approach:
  - Uses POT (Sinkhorn OT) to impute the X-ray-only column Gamma for optical-only GRBs
  - Fills log10NH, log10Fluence, log10PeakFlux from optical measurements (unit-converted)
    instead of leaving them NaN for MICE
  - Matches catalogs via anchor-informed dimension weights learned from the ~88 GRBs
    present in both catalogs
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import ot

OPT_FILE  = 'Data/optical_data.csv'
XRAY_FILE = 'Data/xray_data.csv'
OUT_FILE  = 'Data/superlearner_training_ot.csv'

RATIO_THRESH = 0.5   # same as emcee script: drop if err/|value| > 0.5
OT_REG       = 0.3   # Sinkhorn entropic regularisation

# ── Load optical ────────────────────────────────────────────────────────────────
opt = pd.read_csv(OPT_FILE, index_col=0)
opt.index.name = 'GRB'
opt = opt.reset_index()
for col in opt.columns:
    if col not in ('GRB', 'class'):
        opt[col] = pd.to_numeric(opt[col], errors='coerce')

# Unit fixes: convert optical-native units to match xray column conventions
opt['log10NH']       = np.log10(opt['NH'] * 1e22)          # 10^22 cm^-2 -> log10(cm^-2)
opt['log10Fluence']  = np.log10(opt['Fluence'] * 1e-7)     # optical unit -> erg/cm^2 -> log10
opt['log10PeakFlux'] = np.log10(opt['PeakFlux'])            # direct log10 (same scale as xray)
opt['log10T90']      = np.log10(opt['T90'])

opt = opt.rename(columns={
    'logFa':   'log10Fa',   'logFaErr':  'log10FaErr',
    'logT_a':  'log10Ta',   'logTaErr':  'log10TaErr',
    'betaErr': 'BetaErr',
    'z':       'Redshift_crosscheck',
})

# Normalize GRB ID: append 'A' to suffix-less IDs to match xray convention
opt['GRB_norm'] = opt['GRB'].apply(lambda g: g if g[-1].isalpha() else g + 'A')
print(f'Optical data loaded: {len(opt)} GRBs')

# ── Load xray ───────────────────────────────────────────────────────────────────
xray = pd.read_csv(XRAY_FILE, index_col=0)
xray.index.name = 'GRB'
xray = xray.reset_index()
xray['GRB'] = xray['GRB'].str.replace('GRB', '', regex=False)
xray['GRB_norm'] = xray['GRB'].apply(lambda g: g if g[-1].isalpha() else g + 'A')
xray['log10T90'] = np.log10(pd.to_numeric(xray['T90'], errors='coerce'))
for col in xray.columns:
    if col not in ('GRB', 'GRB_norm'):
        xray[col] = pd.to_numeric(xray[col], errors='coerce')
print(f'X-ray data loaded:  {len(xray)} GRBs')

# ── Find anchors (GRBs in both catalogs) ────────────────────────────────────────
anchors = set(opt['GRB_norm']) & set(xray['GRB_norm'])
print(f'Anchor GRBs:        {len(anchors)}')

# ── Build 6-D shared feature space and standardize jointly ──────────────────────
SHARED = ['log10T90', 'PhotonIndex', 'Alpha', 'Beta', 'log10Fa', 'log10Ta']

opt_z_raw  = opt[SHARED].values.astype(float)
xray_z_raw = xray[SHARED].values.astype(float)

scaler = StandardScaler()
pooled = np.vstack([opt_z_raw, xray_z_raw])
ok_mask = ~np.isnan(pooled).any(axis=1)
scaler.fit(pooled[ok_mask])

def safe_transform(arr):
    out = np.full(arr.shape, np.nan)
    m = ~np.isnan(arr).any(axis=1)
    out[m] = scaler.transform(arr[m])
    return out

opt_z  = safe_transform(opt_z_raw)
xray_z = safe_transform(xray_z_raw)

# ── Anchor-informed dimension weights ───────────────────────────────────────────
common = sorted(anchors)
opt_idx  = {g: i for i, g in enumerate(opt['GRB_norm'])}
xray_idx = {g: i for i, g in enumerate(xray['GRB_norm'])}

opt_anch  = np.array([opt_z[opt_idx[g]]  for g in common])
xray_anch = np.array([xray_z[xray_idx[g]] for g in common])

usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
opt_anch  = opt_anch[usable]
xray_anch = xray_anch[usable]
print(f'Usable anchors for weighting: {usable.sum()}')

diffs   = opt_anch - xray_anch
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = 1.0 / dim_var
weights = weights / weights.mean()
print('Dimension weights:')
for dim, w in zip(SHARED, weights):
    print(f'  {dim}: {w:.3f}')

opt_zw  = opt_z  * np.sqrt(weights)
xray_zw = xray_z * np.sqrt(weights)

# ── Optical-only GRBs: impute Gamma via Sinkhorn OT ────────────────────────────
opt_only_mask = ~opt['GRB_norm'].isin(anchors)
opt_only      = opt[opt_only_mask].copy().reset_index(drop=True)
opt_only_zw   = opt_zw[opt_only_mask.values]
print(f'Optical-only GRBs:  {len(opt_only)}')

xray_gamma    = xray['Gamma'].values
xray_rows_ok  = ~np.isnan(xray_zw).any(axis=1) & ~np.isnan(xray_gamma)
opt_rows_ok   = ~np.isnan(opt_only_zw).any(axis=1)

n_opt  = opt_rows_ok.sum()
n_xray = xray_rows_ok.sum()

M = ot.dist(opt_only_zw[opt_rows_ok], xray_zw[xray_rows_ok])
M = M / M.max()
P = ot.sinkhorn(
    a=np.ones(n_opt)  / n_opt,
    b=np.ones(n_xray) / n_xray,
    M=M, reg=OT_REG
)
P_norm = P / P.sum(axis=1, keepdims=True)
gamma_vals = np.full(len(opt_only), np.nan)
gamma_vals[opt_rows_ok] = P_norm @ xray_gamma[xray_rows_ok]
opt_only['Gamma'] = gamma_vals
print(f'Gamma imputed for {opt_rows_ok.sum()} / {len(opt_only)} optical-only GRBs')

# ── Build final dataset ─────────────────────────────────────────────────────────
OUT_COLS = [
    'GRB', 'Redshift_crosscheck', 'log10T90',
    'log10Fa', 'log10Ta', 'Alpha', 'Beta',
    'Gamma', 'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux',
    'T90Err', 'log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr',
    'FluenceErr', 'PhotonIndexErr', 'PeakFluxErr',
]

xray_out = xray.copy()
xray_out['GRB'] = xray_out['GRB_norm']

opt_only_out = opt_only.copy()
opt_only_out['GRB'] = opt_only_out['GRB_norm']
# optical columns not in xray (T90Err, FluenceErr) will be NaN — MICE handles them

combined = pd.concat([xray_out, opt_only_out], ignore_index=True)

final = pd.DataFrame()
for col in OUT_COLS:
    if col in combined.columns:
        final[col] = combined[col].values
    else:
        final[col] = np.nan

print(f'\nBefore error cut: {len(final)} GRBs')

# ── Relative error cut (same as emcee script, RATIO_THRESH=0.5) ─────────────────
def rel_err_cut(df, val_col, err_col):
    v = pd.to_numeric(df[val_col], errors='coerce').abs()
    e = pd.to_numeric(df[err_col], errors='coerce').abs()
    return (v > 0) & (e / v.clip(lower=1e-10) > RATIO_THRESH)

any_cut = (
    rel_err_cut(final, 'log10Fa', 'log10FaErr') |
    rel_err_cut(final, 'log10Ta', 'log10TaErr') |
    rel_err_cut(final, 'Alpha',   'AlphaErr')   |
    rel_err_cut(final, 'Beta',    'BetaErr')
)
final = final[~any_cut].reset_index(drop=True)
print(f'After error cut:  {len(final)} GRBs (dropped {any_cut.sum()})')

# ── Feature nulling (same as superlearner.R) ────────────────────────────────────
final.loc[np.isinf(pd.to_numeric(final['log10PeakFlux'], errors='coerce')), 'log10PeakFlux'] = np.nan
final.loc[pd.to_numeric(final['log10NH'],      errors='coerce') < 20,  'log10NH']      = np.nan
final.loc[pd.to_numeric(final['Beta'],         errors='coerce') > 3,   'Beta']          = np.nan
final.loc[pd.to_numeric(final['Gamma'],        errors='coerce') > 3,   'Gamma']         = np.nan
final.loc[pd.to_numeric(final['Alpha'],        errors='coerce') > 3,   'Alpha']         = np.nan
final.loc[pd.to_numeric(final['PhotonIndex'],  errors='coerce') < 0,   'PhotonIndex']   = np.nan

# ── Save ────────────────────────────────────────────────────────────────────────
final.to_csv(OUT_FILE, index=True)
print(f'\nSaved to {OUT_FILE}')
print('\nMissing values per column:')
print(final[OUT_COLS[1:]].isnull().sum().to_string())
