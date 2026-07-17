#!/usr/bin/env python3
"""
common.py -- shared data prep for the combination experiments.

Loads the optical/X-ray catalogs, applies the emcee calibration (cached
posterior, no refit) to project optical-only GRBs' shared Dainotti
parameters onto the X-ray scale, and hands back clean dataframes. Each
experiment script only has to supply ONE thing: how to fill in the 5
X-ray-only target features (Gamma, PhotonIndex, log10NH, log10Fluence,
log10PeakFlux) for the optical-only GRBs. Everything else (loading,
projection, error-cut, final CSV assembly) is identical across experiments
so results are only comparable, not doing genuinely different data prep by
diverging on details already fixed in the OT v3 script.

IMPORTANT: unlike the original OT v3 script, the matching/prediction
features made available here for optical-only GRBs are ONLY
[logFa_x, logTa_x, Alpha_x, Beta_x, log10T90] -- redshift ('z') is
deliberately excluded from anything an imputation method is allowed to see,
since it's the training target (see the leakage bug found and fixed in
ot_fusion/combine_optical_xray_ot_v3.py).
"""

import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
DATA_DIR   = os.path.join(REPO_DIR, 'Data')

OPT_FILE    = os.path.join(DATA_DIR, 'OnlyLGRBs_data_171_optical_corrected.csv')
XRAY_FILE   = os.path.join(DATA_DIR, 'Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv')
CHAINS_FILE = os.path.join(DATA_DIR, 'emcee_chains_v2.npz')

MIN_ERR      = 1e-6
N_PROJ_SAMP  = 1000
RATIO_THRESH = 0.5

TARGET_COLS = ['Gamma', 'PhotonIndex', 'log10NH', 'log10Fluence', 'log10PeakFlux']
# Features an imputation method is allowed to use -- no 'z', see docstring.
FEATURE_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x', 'log10T90']


def _norm_grb(s):
    s = str(s).replace('GRB', '').strip()
    return s if s[-1].isalpha() else s + 'A'


def load_data():
    """Returns (xray, opt_clean, opt_proj, opt_only_proj, overlap) where:
    - xray: 223 GRBs, all real measured values (feature + target cols + errors)
    - opt_clean: 161 GRBs, optical catalog with X-ray-only cols nulled for optical-only rows
    - opt_proj: 161 GRBs, emcee-projected shared params (FEATURE_COLS) + z + log10T90
    - opt_only_proj: the 83 optical-only rows of opt_proj
    - overlap: sorted list of the 78 anchor GRB ids
    """
    opt = pd.read_csv(OPT_FILE, index_col=0)
    opt.index = [_norm_grb(g) for g in opt.index]
    opt.index.name = 'GRB'
    for col in opt.columns:
        if col != 'class':
            opt[col] = pd.to_numeric(opt[col], errors='coerce')

    xray = pd.read_csv(XRAY_FILE, index_col=0)
    xray.index = [_norm_grb(g) for g in xray.index]
    xray.index.name = 'GRB'
    xray = xray.rename(columns={
        'Redshift_crosscheck': 'z',
        'log10Fa_X':   'logFa_x',   'log10FaErr_X': 'logFaErr_x',
        'log10Ta_X':   'logTa_x',   'log10TaErr':   'logTaErr_x',
        'Alpha_X':     'Alpha_x',   'AlphaErr':     'AlphaErr_x',
        'Beta_X':      'Beta_x',    'BetaErr_X':    'BetaErr_x',
        'logFluenceErr': 'FluenceErr',       # dex-scale despite the name, see below
        'log10PeakFluxErr': 'PeakFluxErr',   # dex-scale despite the name, see below
    })
    for col in xray.columns:
        xray[col] = pd.to_numeric(xray[col], errors='coerce')

    overlap      = sorted(set(opt.index) & set(xray.index))
    opt_only_idx = ~opt.index.isin(overlap)

    opt_clean = opt.copy()
    xray_feature_cols = ['Gamma', 'PhotonIndex', 'PhotonIndexErr',
                          'log10NH', 'log10Fluence', 'log10FluenceErr',
                          'log10PeakFlux', 'log10PeakFluxErr']
    for col in xray_feature_cols:
        if col in opt_clean.columns:
            opt_clean.loc[opt_only_idx, col] = np.nan

    # Emcee projection: optical -> X-ray scale (cached posterior, no refit)
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
    rng = np.random.default_rng(42)
    opt_proj = pd.DataFrame({'z': opt_clean['z'], 'log10T90': opt_clean['log10T90']},
                             index=opt_clean.index)
    for param, (xc, xe) in params_map.items():
        flat = cache[param]
        idx = rng.integers(0, len(flat), N_PROJ_SAMP)
        flat_s = flat[idx]
        m_s, b_s = flat_s[:, 0], flat_s[:, 1]
        s_s = np.exp(flat_s[:, 2])
        x_v = opt_clean[xc].values.astype(float)
        x_e = opt_clean[xe].values.astype(float)
        x_e = np.where(np.isfinite(x_e), np.maximum(x_e, MIN_ERR), MIN_ERR)
        n = len(x_v)
        x_pert = x_v[None, :] + rng.standard_normal((N_PROJ_SAMP, n)) * x_e[None, :]
        preds = (m_s[:, None] * x_pert + b_s[:, None]
                 + rng.standard_normal((N_PROJ_SAMP, n)) * s_s[:, None])
        yc, ye = proj_out[param]
        opt_proj[yc] = np.nanmedian(preds, axis=0)
        lo = np.nanpercentile(preds, 16, axis=0)
        hi = np.nanpercentile(preds, 84, axis=0)
        opt_proj[ye] = (hi - lo) / 2.0

    opt_only_proj = opt_proj[opt_only_idx].copy()
    return xray, opt_clean, opt_proj, opt_only_proj, overlap


def xray_linear_errors(xray):
    """FluenceErr/PeakFluxErr columns are dex-scale despite their linear-
    sounding names (see the double-conversion bug fixed in OT v3) -- convert
    to real linear errors here, once, for every experiment to reuse."""
    fluence_err_lin  = xray['FluenceErr']  * (10 ** xray['log10Fluence'])  * np.log(10)
    peakflux_err_lin = xray['PeakFluxErr'] * (10 ** xray['log10PeakFlux']) * np.log(10)
    return fluence_err_lin, peakflux_err_lin


def kernel_match_impute(query_df, donor_df, target_cols, match_cols, weights=None, bandwidth=0.2):
    """Shared OT-style kernel-weighted nearest-neighbor imputation (the
    fixed, non-leaky mechanism from combine_optical_xray_ot_v3.py): weight
    ∝ exp(-cost/bandwidth) per query row, cost is a weighted standardized
    Euclidean distance in match_cols space, no donor-side marginal
    constraint. Returns (means_df, stds_df) indexed like query_df.
    """
    from sklearn.preprocessing import StandardScaler
    import ot as pot

    scaler = StandardScaler()
    pooled = np.vstack([query_df[match_cols].values, donor_df[match_cols].values]).astype(float)
    ok = ~np.isnan(pooled).any(axis=1)
    scaler.fit(pooled[ok])

    def safe_scale(df):
        arr = df[match_cols].values.astype(float)
        out = np.full(arr.shape, np.nan)
        m = ~np.isnan(arr).any(axis=1)
        out[m] = scaler.transform(arr[m])
        return out

    q_z = safe_scale(query_df)
    d_z = safe_scale(donor_df)
    w = np.ones(len(match_cols)) if weights is None else np.asarray(weights)
    sqrt_w = np.sqrt(w)

    q_ok = ~np.isnan(q_z).any(axis=1)
    d_ok = ~np.isnan(d_z).any(axis=1)

    global_M = pot.dist(d_z[d_ok] * sqrt_w, d_z[d_ok] * sqrt_w)
    scale = np.median(global_M[global_M > 0])

    M = pot.dist(q_z[q_ok] * sqrt_w, d_z[d_ok] * sqrt_w)
    K = np.exp(-M / (bandwidth * scale))
    P = K / K.sum(axis=1, keepdims=True)

    donor_targets = donor_df[target_cols].values.astype(float)[d_ok]
    means = np.full((len(query_df), len(target_cols)), np.nan)
    stds  = np.full((len(query_df), len(target_cols)), np.nan)
    q_idx = np.where(q_ok)[0]
    for k in range(len(target_cols)):
        t = donor_targets[:, k]
        t_ok = ~np.isnan(t)
        if t_ok.sum() == 0:
            continue
        P_sub = P[:, t_ok]
        P_sub = P_sub / P_sub.sum(axis=1, keepdims=True)
        mu = P_sub @ t[t_ok]
        var = (P_sub * (t[t_ok][None, :] - mu[:, None]) ** 2).sum(axis=1)
        means[q_idx, k] = mu
        stds[q_idx, k] = np.sqrt(var)

    means_df = pd.DataFrame(means, index=query_df.index, columns=target_cols)
    stds_df  = pd.DataFrame(stds,  index=query_df.index, columns=target_cols)
    return means_df, stds_df


def assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path):
    """imputed / imputed_std: dict of {col: array-like over opt_only_proj.index}
    for the 5 TARGET_COLS. Builds the final training CSV (same schema as the
    OT v3 script) and applies the same relative-error quality cut."""
    fluence_err_lin, peakflux_err_lin = xray_linear_errors(xray)

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
        'FluenceErr':          fluence_err_lin,
        'PhotonIndexErr':      xray.get('PhotonIndexErr', np.nan),
        'PeakFluxErr':         peakflux_err_lin,
    }, index=xray.index)

    t90err_median = xray_out['T90Err'].dropna().median()
    opt_only_grbs = opt_only_proj.index

    def get_col(name):
        v = imputed.get(name)
        return pd.Series(v, index=opt_only_grbs) if v is not None else pd.Series(np.nan, index=opt_only_grbs)

    def get_err_lin(name, val_col):
        s = imputed_std.get(name)
        if s is None:
            return pd.Series(np.nan, index=opt_only_grbs)
        s = pd.Series(s, index=opt_only_grbs)
        return s * (10 ** get_col(val_col)) * np.log(10)

    opt_only_out = pd.DataFrame({
        'Redshift_crosscheck': opt_clean.loc[opt_only_grbs, 'z'].values,
        'log10T90':            opt_clean.loc[opt_only_grbs, 'log10T90'].values,
        'log10Fa':             opt_only_proj['logFa_x'].values,
        'log10FaErr':          opt_only_proj['logFaErr_x'].values,
        'log10Ta':             opt_only_proj['logTa_x'].values,
        'log10TaErr':          opt_only_proj['logTaErr_x'].values,
        'Alpha':               opt_only_proj['Alpha_x'].values,
        'AlphaErr':            opt_only_proj['AlphaErr_x'].values,
        'Beta':                opt_only_proj['Beta_x'].values,
        'BetaErr':             opt_only_proj['BetaErr_x'].values,
        'Gamma':               get_col('Gamma').values,
        'log10Fluence':        get_col('log10Fluence').values,
        'PhotonIndex':         get_col('PhotonIndex').values,
        'log10NH':             get_col('log10NH').values,
        'log10PeakFlux':       get_col('log10PeakFlux').values,
        'T90Err':              t90err_median,
        'FluenceErr':          get_err_lin('log10Fluence', 'log10Fluence').values,
        'PhotonIndexErr':      imputed_std.get('PhotonIndex', pd.Series(np.nan, index=opt_only_grbs)),
        'PeakFluxErr':         get_err_lin('log10PeakFlux', 'log10PeakFlux').values,
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

    # Relative error cut, same logic/threshold as OT v3.
    xray_ids = set(xray.index)
    opt_only_mask = ~out.index.isin(xray_ids)
    opt_only_grbs_out = out.index[opt_only_mask]
    opt_orig = opt_clean.loc[opt_clean.index.isin(opt_only_grbs_out)].apply(pd.to_numeric, errors='coerce')

    def _cut_orig(col_v, col_e):
        v = opt_orig[col_v].abs(); e = opt_orig[col_e].abs()
        return ((v > 0) & (e / v.clip(lower=1e-10) > RATIO_THRESH)).reindex(opt_only_grbs_out, fill_value=False)

    def _cut_proj(col_v, col_e):
        sub = out.loc[opt_only_grbs_out].apply(pd.to_numeric, errors='coerce')
        return (sub[col_v].abs() > 0) & (sub[col_e].abs() / sub[col_v].abs().clip(lower=1e-10) > RATIO_THRESH)

    any_cut = (
        (_cut_orig('log10Faopt', 'log10FaErr_opt') | _cut_proj('log10Fa', 'log10FaErr')) |
        (_cut_orig('log10Taopt', 'logTaErr_opt')   | _cut_proj('log10Ta', 'log10TaErr')) |
        (_cut_orig('Alpha_opt',  'AlphaErr_opt')   | _cut_proj('Alpha',   'AlphaErr')) |
        (_cut_orig('Beta_opt',   'betaErr_opt')    | _cut_proj('Beta',    'BetaErr'))
    )
    out_cut = out.drop(index=any_cut[any_cut].index)
    out_cut.to_csv(out_path)
    return out_cut
