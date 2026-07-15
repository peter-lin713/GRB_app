#!/usr/bin/env python3
"""
ot_fusion/combine_optical_xray_ot_v3.py

Same OT imputation as v2 (Sinkhorn-matched, anchor-weighted, in the
emcee-projected [logFa_x, logTa_x, Alpha_x, Beta_x, z] space), but fixes the
root cause of the v2 blowup (r=-0.045, RMSE in the billions): v2 imputed a
point estimate for the 5 X-ray-only features but no uncertainty, so the
downstream MICE/MC-propagation step in superlearner.R had to guess a sigma
for those rows and occasionally hallucinated an extreme one.

Fix: the Sinkhorn transport plan already gives a soft weighting over X-ray
donor GRBs for each optical-only GRB. Instead of only taking the weighted
mean (barycenter) of the donors' feature values, this version also computes
the weighted variance of the same donor pool around that mean -- a "how much
do the donors agree" uncertainty, native to the OT computation itself. A
tight match -> small uncertainty; a diffuse match -> large uncertainty. That
becomes the real per-GRB error estimate fed into the output error columns,
so MICE never sees a NaN error to hallucinate from.

A leave-one-out validation loop over the anchor GRBs (the ones with real
X-ray data) checks whether this is a good idea: for each anchor, pretend
it's optical-only, impute its 5 target features from the *other* anchors +
X-ray-only GRBs, and compare the imputed value (and its barycentric
uncertainty) against the true measured value.
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
OUT_FILE     = os.path.join(DATA_DIR, 'superlearner_training_ot_v3.csv')
LOO_REPORT   = os.path.join(DATA_DIR, 'ot_v3_loo_validation.csv')

#' Two compounding bugs in v1/v2's OT_REG=0.1 balanced-Sinkhorn matching:
#' (1) it over-smooths -- 0.1 spreads each optical-only GRB's match across
#'     ~190 of the 223 X-ray donors (effective donor count via inverse
#'     Simpson index), collapsing every imputed value toward the population
#'     mean regardless of which GRB is being matched.
#' (2) worse, *balanced* OT (ot.sinkhorn's hard donor-side marginal, which
#'     forces every donor to receive exactly 1/n_d total mass) is the wrong
#'     tool for "one query vs many donors" matching at all: with a single
#'     query row (as in the leave-one-out check below), that marginal
#'     constraint forces near-uniform weights *regardless of cost or reg*,
#'     since one source alone must satisfy the full uniform target marginal.
#'     This is a structural issue, not a tuning one -- confirmed with a toy
#'     example where a donor at cost~0 still only gets weight 1/n_d until reg
#'     is pushed low enough to be numerically unstable.
#' Fix: use a plain per-row Gibbs/softmax kernel (weight ∝ exp(-cost/bw),
#' normalized independently per query row) instead of balanced OT. This is
#' the standard tool for "impute from weighted similar neighbors" and has no
#' donor-marginal degeneracy. `OT_REG` here is the kernel bandwidth, applied
#' to cost distances normalized by a *fixed* global scale (the median
#' pairwise distance among the X-ray donors) so the effective bandwidth means
#' the same thing whether imputing one query or eighty-three at once.
OT_REG       = 0.15
MIN_ERR      = 1e-6
N_PROJ_SAMP  = 1000
RATIO_THRESH = 0.5

# Cap any barycentric error estimate at this multiple of the X-ray
# population's own spread for that feature -- a diffuse/low-confidence OT
# match shouldn't be allowed to inject an arbitrarily large sigma downstream.
ERR_CAP_MULT = 3.0

TARGET_COLS = ['Gamma', 'PhotonIndex', 'log10NH', 'log10Fluence', 'log10PeakFlux']
MATCH_COLS  = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x', 'z']

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

overlap      = sorted(set(opt.index) & set(xray.index))
opt_only_idx = ~opt.index.isin(overlap)
print(f'Overlap: {len(overlap)}, Optical-only: {opt_only_idx.sum()}')

XRAY_FEATURE_COLS = ['Gamma', 'PhotonIndex', 'PhotonIndexErr',
                      'log10NH', 'log10Fluence', 'log10FluenceErr',
                      'log10PeakFlux', 'log10PeakFluxErr']
opt_clean = opt.copy()
for col in XRAY_FEATURE_COLS:
    if col in opt_clean.columns:
        opt_clean.loc[opt_only_idx, col] = np.nan
print('Nulled MICE-imputed X-ray features for optical-only GRBs.')

# ── Emcee projection: optical -> X-ray scale ──────────────────────────────────
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

    yc, ye       = proj_out[param]
    opt_proj[yc] = np.nanmedian(preds, axis=0)
    lo           = np.nanpercentile(preds, 16, axis=0)
    hi           = np.nanpercentile(preds, 84, axis=0)
    opt_proj[ye] = (hi - lo) / 2.0
    print(f'  {param}: median 1σ projected err = {np.nanmedian(opt_proj[ye]):.4f}')

opt_only_proj = opt_proj[opt_only_idx].copy()

# ── Joint standardization over the matching space ─────────────────────────────
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

opt_z_all  = safe_scale(opt_match)
xray_z_all = safe_scale(xray_match)

opt_anch_df  = opt_proj.loc[overlap]
xray_anch_df = xray.loc[overlap]
opt_anch     = opt_anch_df[MATCH_COLS].values.astype(float)
xray_anch    = xray_anch_df[MATCH_COLS].values.astype(float)
usable       = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
print(f'\nUsable anchors: {usable.sum()} / {len(overlap)}')

opt_anch_z  = safe_scale(opt_anch[usable])
xray_anch_z = safe_scale(xray_anch[usable])

diffs   = opt_anch_z - xray_anch_z
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = 1.0 / dim_var
weights = weights / weights.mean()
print('\nDimension weights:')
for d, w in zip(MATCH_COLS, weights):
    print(f'  {d}: {w:.3f}')

sqrt_w = np.sqrt(weights)

#' Fixed global cost scale so `OT_REG` (the kernel bandwidth) means the same
#' thing for a single-row LOO query as for the full 83-row batch -- computed
#' once from the X-ray donor pool's own pairwise distances, not from
#' whatever query happens to be passed to a given ot_impute() call.
_xray_ok_global = ~np.isnan(xray_z_all).any(axis=1)
_M_global = ot.dist(xray_z_all[_xray_ok_global] * sqrt_w, xray_z_all[_xray_ok_global] * sqrt_w)
GLOBAL_COST_SCALE = np.median(_M_global[_M_global > 0])
print(f'\nGlobal cost scale (median pairwise donor distance): {GLOBAL_COST_SCALE:.4f}')


def ot_impute(query_z, donor_z, donor_targets, reg=OT_REG):
    """Gibbs-kernel-match `query_z` rows against `donor_z` rows (weight ∝
    exp(-cost/bandwidth), normalized per query row -- no donor-side marginal
    constraint, unlike balanced Sinkhorn OT). Returns per-row imputed value
    (weighted mean) and barycentric std (weighted spread of the donor pool
    around that mean) for each column of `donor_targets`.
    """
    q_ok = ~np.isnan(query_z).any(axis=1)
    d_ok = ~np.isnan(donor_z).any(axis=1)
    n_q, n_d = q_ok.sum(), d_ok.sum()

    n_feat = donor_targets.shape[1]
    means = np.full((len(query_z), n_feat), np.nan)
    stds  = np.full((len(query_z), n_feat), np.nan)
    if n_q == 0 or n_d == 0:
        return means, stds

    M = ot.dist(query_z[q_ok] * sqrt_w, donor_z[d_ok] * sqrt_w)
    K = np.exp(-M / (reg * GLOBAL_COST_SCALE))
    P_norm = K / K.sum(axis=1, keepdims=True)

    targets = donor_targets[d_ok]
    q_idx = np.where(q_ok)[0]
    for k in range(n_feat):
        t    = targets[:, k]
        t_ok = ~np.isnan(t)
        if t_ok.sum() == 0:
            continue
        P_sub = P_norm[:, t_ok]
        P_sub = P_sub / P_sub.sum(axis=1, keepdims=True)
        mu    = P_sub @ t[t_ok]
        var   = (P_sub * (t[t_ok][None, :] - mu[:, None]) ** 2).sum(axis=1)
        means[q_idx, k] = mu
        stds[q_idx, k]  = np.sqrt(var)
    return means, stds


anchor_ids   = np.array(overlap)[usable]
anchor_z_all = opt_z_all[np.array([opt.index.get_loc(g) for g in anchor_ids])]


def run_loo(reg):
    """Leave-one-out impute every anchor GRB's 5 target features from the
    other X-ray donors, using bandwidth `reg`. Returns the per-anchor result
    table and the mean correlation across features (the bandwidth-selection
    score).
    """
    rows = []
    for i, grb in enumerate(anchor_ids):
        donor_mask    = (xray.index != grb).values if hasattr(xray.index != grb, 'values') else (xray.index != grb)
        donor_z       = xray_z_all[donor_mask]
        donor_targets = xray[TARGET_COLS].values.astype(float)[donor_mask]

        query_z = anchor_z_all[i:i + 1]
        mu, sd  = ot_impute(query_z, donor_z, donor_targets, reg=reg)

        true_vals = xray.loc[grb, TARGET_COLS].values.astype(float)
        row = {'GRB': grb}
        for k, col in enumerate(TARGET_COLS):
            row[f'{col}_true']    = true_vals[k]
            row[f'{col}_imputed'] = mu[0, k]
            row[f'{col}_sigma']   = sd[0, k]
        rows.append(row)

    df = pd.DataFrame(rows).set_index('GRB')
    corrs = []
    for col in TARGET_COLS:
        t, p = df[f'{col}_true'].values, df[f'{col}_imputed'].values
        ok = ~(np.isnan(t) | np.isnan(p))
        if ok.sum() >= 3 and p[ok].std() > 1e-8:
            corrs.append(np.corrcoef(t[ok], p[ok])[0, 1])
    return df, (np.mean(corrs) if corrs else -1.0)


# ── Bandwidth selection: sweep reg, pick the one with best LOO correlation ────
print('\nBandwidth sweep (mean correlation across the 5 target features, LOO):')
sweep_results = {}
for reg in [1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03]:
    _, score = run_loo(reg)
    sweep_results[reg] = score
    print(f'  reg={reg:<6} mean LOO corr = {score:.4f}')

OT_REG = max(sweep_results, key=sweep_results.get)
print(f'\nSelected bandwidth: reg={OT_REG} (best mean LOO correlation = {sweep_results[OT_REG]:.4f})')

# ── Impute the 5 target features for optical-only GRBs, using the selected bandwidth ──
xray_targets_full = xray[TARGET_COLS].values.astype(float)
opt_only_z = opt_z_all[np.array(opt_only_idx)]

imp_mean, imp_std = ot_impute(opt_only_z, xray_z_all, xray_targets_full, reg=OT_REG)

print('\nImputation (with barycentric uncertainty):')
for i, col in enumerate(TARGET_COLS):
    n_ok = (~np.isnan(imp_mean[:, i])).sum()
    med_std = np.nanmedian(imp_std[:, i])
    print(f'  {col}: imputed for {n_ok} / {len(imp_mean)} optical-only GRBs, '
          f'median barycentric 1σ = {med_std:.4f}')

imputed     = {col: imp_mean[:, i] for i, col in enumerate(TARGET_COLS)}
imputed_std = {col: imp_std[:, i]  for i, col in enumerate(TARGET_COLS)}

# Safety cap: clip any barycentric sigma at ERR_CAP_MULT x the X-ray
# population's own spread for that feature, so a diffuse/low-confidence
# match can't inject an arbitrarily large error downstream.
pop_std = {col: np.nanstd(xray[col].values.astype(float)) for col in TARGET_COLS}
for col in TARGET_COLS:
    cap = ERR_CAP_MULT * pop_std[col]
    n_capped = (imputed_std[col] > cap).sum()
    if n_capped > 0:
        print(f'  capping {col}: {n_capped} GRB(s) had barycentric σ > {cap:.4f} '
              f'({ERR_CAP_MULT}x population std) -- clipped')
    imputed_std[col] = np.minimum(imputed_std[col], cap)

# ── Leave-one-out validation report at the selected bandwidth ─────────────────
print('\n' + '=' * 70)
print('Leave-one-out validation on anchor GRBs (selected bandwidth)')
print('=' * 70)

loo_df, _ = run_loo(OT_REG)
loo_df.to_csv(LOO_REPORT)

print(f'\n{"feature":<15}{"n":>5}{"corr":>8}{"RMSE":>10}{"|resid|/σ<1":>14}')
for col in TARGET_COLS:
    true_v = loo_df[f'{col}_true'].values
    imp_v  = loo_df[f'{col}_imputed'].values
    sig_v  = loo_df[f'{col}_sigma'].values
    ok = ~(np.isnan(true_v) | np.isnan(imp_v))
    if ok.sum() < 3:
        print(f'{col:<15}{ok.sum():>5}   (too few for stats)')
        continue
    corr = np.corrcoef(true_v[ok], imp_v[ok])[0, 1]
    rmse = np.sqrt(np.mean((true_v[ok] - imp_v[ok]) ** 2))
    within_1sig = np.abs(true_v[ok] - imp_v[ok]) < np.maximum(sig_v[ok], 1e-6)
    coverage = within_1sig.mean()
    print(f'{col:<15}{ok.sum():>5}{corr:>8.3f}{rmse:>10.4f}{coverage:>13.1%}')
print('\n(coverage is the fraction of held-out anchors whose true value fell')
print(' within ±1 barycentric σ of the imputed value -- well-calibrated ≈ 68%)')
print(f'\nFull per-anchor LOO table saved to {LOO_REPORT}')
print('=' * 70)

# ── Build final dataset ─────────────────────────────────────────────────────────
# FluenceErr/PeakFluxErr: the xray.rename() above maps the catalog's
# log-space columns (logFluenceErr, log10PeakFluxErr) straight to
# 'FluenceErr'/'PeakFluxErr' with NO unit conversion -- so xray['FluenceErr']
# here is actually still a dex-scale value, just mislabeled as linear. Left
# as-is, superlearner.R's to_dex_err() treats it as linear and divides by
# 10**log10Fluence again -- a double dex-conversion that blows up for any
# GRB with very negative log10Fluence (tiny denominator). GRB 090927A hit
# this exactly: a normal ~0.03 dex error came out as log10FluenceErr=141458
# after the double conversion, which then blew up MC error propagation.
# Fix: treat xray['FluenceErr']/['PeakFluxErr'] as the dex values they
# actually are and convert to real linear errors here, so to_dex_err's single
# conversion downstream lands back on the original, sane dex value.
xray_fluence_err_lin  = xray['FluenceErr']  * (10 ** xray['log10Fluence'])  * np.log(10)
xray_peakflux_err_lin = xray['PeakFluxErr'] * (10 ** xray['log10PeakFlux']) * np.log(10)

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
    'FluenceErr':          xray_fluence_err_lin,
    'PhotonIndexErr':      xray.get('PhotonIndexErr', np.nan),
    'PeakFluxErr':         xray_peakflux_err_lin,
}, index=xray.index)

# T90Err has no optical-side source at all (neither catalog carries a T90
# uncertainty for optical-only GRBs) -- fall back to the X-ray population
# median, same as v2. This is unrelated to the OT imputation fix.
t90err_median = xray_out['T90Err'].dropna().median()

# Convert the barycentric log-space std back to a linear-scale error for the
# two features whose output column is linear (FluenceErr, PeakFluxErr) --
# superlearner.R re-derives the dex error via to_dex_err(lin_err, log10_val),
# so this round-trips back to exactly imputed_std at ingestion time.
fluence_err_lin  = imputed_std['log10Fluence']  * (10 ** imputed['log10Fluence'])  * np.log(10)
peakflux_err_lin = imputed_std['log10PeakFlux'] * (10 ** imputed['log10PeakFlux']) * np.log(10)

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
    'T90Err':              t90err_median,
    'FluenceErr':          fluence_err_lin,
    'PhotonIndexErr':      imputed_std['PhotonIndex'],
    'PeakFluxErr':         peakflux_err_lin,
}, index=opt_only_grbs)

# Note: Gamma and log10NH have no error column anywhere in the pipeline's
# ingest contract (never did, even for X-ray-measured rows) -- their
# barycentric uncertainty is reported above for QC/LOO purposes only and
# isn't written out, since there's no slot for it downstream.

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

# ── Error ratio cut (same dual-sided logic as v1/v2) ──────────────────────────
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

print(f'Error cut: removed {len(drop_ids)} optical-only GRBs -> {out_cut.shape[0]} remain')
print(f'  logFa: {fa_cut.sum()}  logTa: {ta_cut.sum()}  Alpha: {alpha_cut.sum()}  Beta: {beta_cut.sum()}')

t90_miss = out_cut['log10T90'].isna().sum()
print(f'\nlog10T90 missing: {t90_miss}')
if t90_miss == 0:
    print('  OK: all GRBs have log10T90')

# Sanity check: no NaN error columns feeding into MICE/MC downstream.
err_cols = ['T90Err', 'log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr',
            'FluenceErr', 'PhotonIndexErr', 'PeakFluxErr']
nan_errs = out_cut[err_cols].isna().sum()
print('\nNaN counts in error columns (should be 0 to avoid the v2 failure mode):')
print(nan_errs.to_string())

# ── Save ──────────────────────────────────────────────────────────────────────
cut_file = OUT_FILE.replace('.csv', '_errcut_relative.csv')
out_cut.to_csv(cut_file)
out.to_csv(OUT_FILE)
print(f'\nSaved: {OUT_FILE}')
print(f'Saved: {cut_file}')
print(f'\nNon-null counts (error-cut dataset):')
print(out_cut[col_order].notna().sum().to_string())
