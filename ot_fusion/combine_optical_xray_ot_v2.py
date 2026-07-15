#!/opt/anaconda3/envs/grb/bin/python3
"""
ot_fusion/combine_optical_xray_ot_v2.py

OT imputation of X-ray features for optical-only GRBs.

Inputs:
  - OnlyLGRBs_data_171_optical_corrected.csv:
      overlapping GRBs already have real X-ray values;
      optical-only GRBs have MICE-imputed X-ray values (nulled out before OT).
  - Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv

Matching space : emcee-projected [logFa_x, logTa_x, Alpha_x, Beta_x] + z
                 No T90 (same across catalogs). No PhotonIndex (X-ray only).
Target imputed : Gamma, PhotonIndex, log10NH, log10Fluence, log10PeakFlux
Anchor weights : learned from 78 overlapping GRBs (residual variance of
                 projected-optical vs real-xray in standardized space).
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import ot

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
REPO_DIR     = os.path.dirname(SCRIPT_DIR)
DATA_DIR     = os.path.join(REPO_DIR, 'Data')

OPT_FILE     = os.path.join(DATA_DIR, 'OnlyLGRBs_data_171_optical_corrected.csv')
XRAY_FILE    = os.path.join(DATA_DIR, 'Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv')
CHAINS_FILE  = os.path.join(DATA_DIR, 'emcee_chains_v2.npz')
OUT_FILE     = os.path.join(DATA_DIR, 'superlearner_training_ot_v2.csv')

OT_REG       = 0.1
MIN_ERR      = 1e-6
N_PROJ_SAMP  = 1000
RATIO_THRESH = 0.5

# ── Load optical (corrected) ──────────────────────────────────────────────────
opt = pd.read_csv(OPT_FILE, index_col=0)

def norm_grb(s):
    s = str(s).replace('GRB', '').strip()
    return s if s[-1].isalpha() else s + 'A'

opt.index = [norm_grb(g) for g in opt.index]
opt.index.name = 'GRB'

for col in opt.columns:
    if col != 'class':
        opt[col] = pd.to_numeric(opt[col], errors='coerce')

print(f'Optical: {len(opt)} GRBs')

# ── Load xray ─────────────────────────────────────────────────────────────────
xray = pd.read_csv(XRAY_FILE, index_col=0)
xray.index = [norm_grb(g) for g in xray.index]
xray.index.name = 'GRB'
xray = xray.rename(columns={
    'Redshift_crosscheck': 'z',
    'log10Fa_X':   'logFa_x',   'log10FaErr_X': 'logFaErr_x',
    'log10Ta_X':   'logTa_x',   'log10TaErr':   'logTaErr_x',
    'Alpha_X':     'Alpha_x',   'AlphaErr':     'AlphaErr_x',
    'Beta_X':      'Beta_x',    'BetaErr_X':    'BetaErr_x',
    'logFluenceErr': 'FluenceErr',
    'log10PeakFluxErr': 'PeakFluxErr',
})
for col in xray.columns:
    xray[col] = pd.to_numeric(xray[col], errors='coerce')

print(f'X-ray:   {len(xray)} GRBs')

# ── Identify overlap and optical-only ────────────────────────────────────────
overlap      = sorted(set(opt.index) & set(xray.index))
opt_only_idx = ~opt.index.isin(overlap)
print(f'Overlap: {len(overlap)}, Optical-only: {opt_only_idx.sum()}')

# ── Null out MICE-imputed X-ray features for optical-only GRBs ───────────────
# Overlapping GRBs already have real X-ray values; optical-only have MICE estimates.
# Replace with NaN so OT fills them cleanly.
XRAY_FEATURE_COLS = ['Gamma', 'PhotonIndex', 'PhotonIndexErr',
                      'log10NH', 'log10Fluence', 'log10FluenceErr',
                      'log10PeakFlux', 'log10PeakFluxErr']
opt_clean = opt.copy()
for col in XRAY_FEATURE_COLS:
    if col in opt_clean.columns:
        opt_clean.loc[opt_only_idx, col] = np.nan

print('Nulled MICE-imputed X-ray features for optical-only GRBs.')

# ── Emcee projection: optical → X-ray scale ──────────────────────────────────
# Column map: optical (1) file uses log10Faopt/log10Taopt/Alpha_opt/Beta_opt
cache = np.load(CHAINS_FILE, allow_pickle=False)

params_map = {
    'logFa': ('log10Faopt', 'log10FaErr_opt'),
    'logTa': ('log10Taopt', 'logTaErr_opt'),
    'Alpha': ('Alpha_opt',  'AlphaErr_opt'),
    'Beta':  ('Beta_opt',   'betaErr_opt'),
}
proj_out = {
    'logFa': ('logFa_x', 'logFaErr_x'),
    'logTa': ('logTa_x', 'logTaErr_x'),
    'Alpha': ('Alpha_x', 'AlphaErr_x'),
    'Beta':  ('Beta_x',  'BetaErr_x'),
}

rng      = np.random.default_rng(42)
opt_proj = pd.DataFrame({'z': opt_clean['z'], 'log10T90': opt_clean['log10T90']},
                         index=opt_clean.index)

print('\nEmcee projection (cached v2 chains):')
for param, (xc, xe) in params_map.items():
    flat   = cache[param]
    idx    = rng.integers(0, len(flat), N_PROJ_SAMP)
    flat_s = flat[idx]
    m_s, b_s = flat_s[:, 0], flat_s[:, 1]
    s_s = np.exp(flat_s[:, 2])

    x_v = opt_clean[xc].values.astype(float)
    x_e = opt_clean[xe].values.astype(float)
    x_e = np.where(np.isfinite(x_e), np.maximum(x_e, MIN_ERR), MIN_ERR)

    n      = len(x_v)
    x_pert = x_v[None, :] + rng.standard_normal((N_PROJ_SAMP, n)) * x_e[None, :]
    preds  = (m_s[:, None] * x_pert + b_s[:, None]
              + rng.standard_normal((N_PROJ_SAMP, n)) * s_s[:, None])

    yc, ye           = proj_out[param]
    opt_proj[yc]     = np.nanmedian(preds, axis=0)
    lo               = np.nanpercentile(preds, 16, axis=0)
    hi               = np.nanpercentile(preds, 84, axis=0)
    opt_proj[ye]     = (hi - lo) / 2.0
    print(f'  {param}: median 1σ projected err = {np.nanmedian(opt_proj[ye]):.4f}')

# Optical-only subset with projected features
opt_only_proj = opt_proj[opt_only_idx].copy()

# ── Matching space: [logFa_x, logTa_x, Alpha_x, Beta_x, z] ─────────────────
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x', 'z']

opt_anch_df  = opt_proj.loc[overlap]
xray_anch_df = xray.loc[overlap]

opt_anch  = opt_anch_df[MATCH_COLS].values.astype(float)
xray_anch = xray_anch_df[MATCH_COLS].values.astype(float)

usable    = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
print(f'\nUsable anchors: {usable.sum()} / {len(overlap)}')
opt_anch  = opt_anch[usable]
xray_anch = xray_anch[usable]

# ── Joint standardization ─────────────────────────────────────────────────────
opt_match  = opt_proj[MATCH_COLS].values.astype(float)
xray_match = xray[MATCH_COLS].values.astype(float)

scaler = StandardScaler()
pooled = np.vstack([opt_match, xray_match])
scaler.fit(pooled[~np.isnan(pooled).any(axis=1)])

def safe_scale(arr):
    out = np.full(arr.shape, np.nan)
    m   = ~np.isnan(arr).any(axis=1)
    out[m] = scaler.transform(arr[m])
    return out

opt_anch_z  = safe_scale(opt_anch)
xray_anch_z = safe_scale(xray_anch)
opt_z_all   = safe_scale(opt_match)
xray_z_all  = safe_scale(xray_match)

# ── Anchor-informed dimension weights ─────────────────────────────────────────
diffs   = opt_anch_z - xray_anch_z
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = 1.0 / dim_var
weights = weights / weights.mean()
print('\nDimension weights:')
for d, w in zip(MATCH_COLS, weights):
    print(f'  {d}: {w:.3f}')

# ── Sinkhorn OT ───────────────────────────────────────────────────────────────
opt_only_z = opt_z_all[np.array(opt_only_idx)] * np.sqrt(weights)
xray_z_w   = xray_z_all * np.sqrt(weights)

opt_ok  = ~np.isnan(opt_only_z).any(axis=1)
xray_ok = ~np.isnan(xray_z_w).any(axis=1)
n_opt, n_xray = opt_ok.sum(), xray_ok.sum()
print(f'\nOT: {n_opt} optical-only × {n_xray} xray GRBs')

M = ot.dist(opt_only_z[opt_ok], xray_z_w[xray_ok])
M = M / M.max()
P = ot.sinkhorn(
    a=np.ones(n_opt)  / n_opt,
    b=np.ones(n_xray) / n_xray,
    M=M, reg=OT_REG,
)
P_norm = P / P.sum(axis=1, keepdims=True)

# ── Impute target features ────────────────────────────────────────────────────
TARGET_COLS  = ['Gamma', 'PhotonIndex', 'log10NH', 'log10Fluence', 'log10PeakFlux']
xray_targets = xray[TARGET_COLS].values.astype(float)[xray_ok]

print('\nImputation:')
imputed = {}
for i, col in enumerate(TARGET_COLS):
    t    = xray_targets[:, i]
    t_ok = ~np.isnan(t)
    vals = np.full(opt_only_idx.sum(), np.nan)
    if t_ok.sum() > 0:
        P_sub        = P_norm[:, t_ok]
        P_sub        = P_sub / P_sub.sum(axis=1, keepdims=True)
        vals[opt_ok] = P_sub @ t[t_ok]
    imputed[col] = vals
    print(f'  {col}: imputed for {(~np.isnan(vals)).sum()} / {len(vals)} optical-only GRBs')

# ── Build final dataset ───────────────────────────────────────────────────────
# X-ray GRBs: real values directly from xray file
xray_out = pd.DataFrame({
    'Redshift_crosscheck': xray['z'],
    'log10T90':            xray['log10T90'],
    'log10Fa':             xray['logFa_x'],
    'log10FaErr':          xray['logFaErr_x'],
    'log10Ta':             xray['logTa_x'],
    'log10TaErr':          xray['logTaErr_x'],
    'Alpha':               xray['Alpha_x'],
    'AlphaErr':            xray['AlphaErr_x'],
    'Beta':                xray['Beta_x'],
    'BetaErr':             xray['BetaErr_x'],
    'Gamma':               xray['Gamma'],
    'log10Fluence':        xray['log10Fluence'],
    'PhotonIndex':         xray['PhotonIndex'],
    'log10NH':             xray['log10NH'],
    'log10PeakFlux':       xray['log10PeakFlux'],
    'T90Err':              xray.get('T90Err', np.nan),
    'FluenceErr':          xray.get('FluenceErr', np.nan),
    'PhotonIndexErr':      xray.get('PhotonIndexErr', np.nan),
    'PeakFluxErr':         xray.get('PeakFluxErr', np.nan),
}, index=xray.index)

# Median errors from X-ray GRBs — used as fill for OT-imputed GRBs that have no
# true error estimates. This prevents MICE from hallucinating extreme values and
# keeps MC error propagation in a physically reasonable range.
err_medians = {
    'T90Err':         xray_out['T90Err'].dropna().median(),
    'FluenceErr':     xray_out['FluenceErr'].dropna().median(),
    'PhotonIndexErr': xray_out['PhotonIndexErr'].dropna().median(),
    'PeakFluxErr':    xray_out['PeakFluxErr'].dropna().median(),
}
print('\nError fill medians (from X-ray population):')
for k, v in err_medians.items():
    print(f'  {k}: {v:.6f}')

# Optical-only GRBs: emcee-projected light curve params + OT-imputed X-ray features
opt_only_grbs = opt_clean.index[opt_only_idx]
opt_only_out = pd.DataFrame({
    'Redshift_crosscheck': opt_clean.loc[opt_only_idx, 'z'].values,
    'log10T90':            opt_clean.loc[opt_only_idx, 'log10T90'].values,
    'log10Fa':             opt_only_proj['logFa_x'].values,
    'log10FaErr':          opt_only_proj['logFaErr_x'].values,
    'log10Ta':             opt_only_proj['logTa_x'].values,
    'log10TaErr':          opt_only_proj['logTaErr_x'].values,
    'Alpha':               opt_only_proj['Alpha_x'].values,
    'AlphaErr':            opt_only_proj['AlphaErr_x'].values,
    'Beta':                opt_only_proj['Beta_x'].values,
    'BetaErr':             opt_only_proj['BetaErr_x'].values,
    'Gamma':               imputed['Gamma'],
    'log10Fluence':        imputed['log10Fluence'],
    'PhotonIndex':         imputed['PhotonIndex'],
    'log10NH':             imputed['log10NH'],
    'log10PeakFlux':       imputed['log10PeakFlux'],
    'T90Err':              err_medians['T90Err'],
    'FluenceErr':          err_medians['FluenceErr'],
    'PhotonIndexErr':      err_medians['PhotonIndexErr'],
    'PeakFluxErr':         err_medians['PeakFluxErr'],
}, index=opt_only_grbs)

col_order = [
    'Redshift_crosscheck', 'log10T90',
    'log10Fa', 'log10Ta', 'Alpha', 'Beta',
    'Gamma', 'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux',
    'T90Err', 'log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr',
    'FluenceErr', 'PhotonIndexErr', 'PeakFluxErr',
]

out = pd.concat([xray_out, opt_only_out])[col_order]
out = out[pd.to_numeric(out['Redshift_crosscheck'], errors='coerce').notna()].copy()
out.index.name = 'GRB'
print(f'\nBefore error cut: {len(out)} GRBs')

# ── Error ratio cut (same dual-sided logic as emcee v1) ──────────────────────
xray_ids          = set(xray.index)
opt_only_mask_out = ~out.index.isin(xray_ids)
opt_only_grbs_out = out.index[opt_only_mask_out]

opt_orig = opt_clean.loc[opt_clean.index.isin(opt_only_grbs_out)].apply(pd.to_numeric, errors='coerce')

def _cut_orig(col_v, col_e):
    v = opt_orig[col_v].abs()
    e = opt_orig[col_e].abs()
    return ((v > 0) & (e / v.clip(lower=1e-10) > RATIO_THRESH)).reindex(opt_only_grbs_out, fill_value=False)

def _cut_proj(col_v, col_e):
    sub = out.loc[opt_only_grbs_out].apply(pd.to_numeric, errors='coerce')
    return (sub[col_v].abs() > 0) & (sub[col_e].abs() / sub[col_v].abs().clip(lower=1e-10) > RATIO_THRESH)

fa_cut    = _cut_orig('log10Faopt', 'log10FaErr_opt') | _cut_proj('log10Fa', 'log10FaErr')
ta_cut    = _cut_orig('log10Taopt', 'logTaErr_opt')   | _cut_proj('log10Ta', 'log10TaErr')
alpha_cut = _cut_orig('Alpha_opt',  'AlphaErr_opt')   | _cut_proj('Alpha',   'AlphaErr')
beta_cut  = _cut_orig('Beta_opt',   'betaErr_opt')    | _cut_proj('Beta',    'BetaErr')

any_cut  = fa_cut | ta_cut | alpha_cut | beta_cut
drop_ids = any_cut[any_cut].index
out_cut  = out.drop(index=drop_ids)

print(f'Error cut: removed {len(drop_ids)} optical-only GRBs → {out_cut.shape[0]} remain')
print(f'  logFa: {fa_cut.sum()}  logTa: {ta_cut.sum()}  Alpha: {alpha_cut.sum()}  Beta: {beta_cut.sum()}')

# ── T90 check ─────────────────────────────────────────────────────────────────
t90_miss = out_cut['log10T90'].isna().sum()
print(f'\nlog10T90 missing: {t90_miss}')
if t90_miss == 0:
    print('  OK: all GRBs have log10T90')

# ── Save ──────────────────────────────────────────────────────────────────────
cut_file = OUT_FILE.replace('.csv', '_errcut_relative.csv')
out_cut.to_csv(cut_file)
out.to_csv(OUT_FILE)
print(f'\nSaved: {OUT_FILE}')
print(f'Saved: {cut_file}')
print(f'\nNon-null counts (error-cut dataset):')
print(out_cut[col_order].notna().sum().to_string())
