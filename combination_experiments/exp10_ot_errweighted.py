#!/usr/bin/env python3
"""exp10: OT matching (4 Dainotti params, as exp08) with an adaptive
per-query bandwidth: a GRB whose own calibrated parameters carry large
relative errors gets a wider/softer match (blends more donors), while a
precisely-measured GRB gets a sharper, more local match. Tests whether
accounting for each optical-only GRB's own measurement quality -- not just
the fixed anchor-informed dimension weights -- improves the imputation."""
import os
import numpy as np
from sklearn.preprocessing import StandardScaler
import ot as pot
from common import load_data, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']
ERR_COLS   = ['logFaErr_x', 'logTaErr_x', 'AlphaErr_x', 'BetaErr_x']

opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
scaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
diffs = scaler.transform(opt_anch[usable]) - scaler.transform(xray_anch[usable])
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()
sqrt_w = np.sqrt(weights)

xray_z = scaler.transform(xray[MATCH_COLS].values.astype(float))
d_ok = ~np.isnan(xray_z).any(axis=1)
global_M = pot.dist(xray_z[d_ok] * sqrt_w, xray_z[d_ok] * sqrt_w)
scale = np.median(global_M[global_M > 0])
donor_targets = xray[TARGET_COLS].values.astype(float)[d_ok]

# Per-query relative-error score: mean(|err| / (|val|+eps)) across match dims.
rel_err = (opt_only_proj[ERR_COLS].values / (np.abs(opt_only_proj[MATCH_COLS].values) + 1e-3)).mean(axis=1)
rel_err = np.nan_to_num(rel_err, nan=np.nanmedian(rel_err))
base_bw = 0.2
adaptive_bw = base_bw * (1.0 + np.clip(rel_err, 0, 3))  # cap so one huge error doesn't fully flatten a row
print('Adaptive bandwidth range:', adaptive_bw.min(), '-', adaptive_bw.max())

q_z = scaler.transform(opt_only_proj[MATCH_COLS].values.astype(float))
q_ok = ~np.isnan(q_z).any(axis=1)
M = pot.dist(q_z[q_ok] * sqrt_w, xray_z[d_ok] * sqrt_w)

means = np.full((len(opt_only_proj), len(TARGET_COLS)), np.nan)
stds  = np.full((len(opt_only_proj), len(TARGET_COLS)), np.nan)
q_idx = np.where(q_ok)[0]
for i, qi in enumerate(q_idx):
    bw = adaptive_bw[qi]
    K = np.exp(-M[i] / (bw * scale))
    P = K / K.sum()
    for k, col in enumerate(TARGET_COLS):
        t = donor_targets[:, k]
        t_ok = ~np.isnan(t)
        p_sub = P[t_ok] / P[t_ok].sum()
        mu = (p_sub * t[t_ok]).sum()
        var = (p_sub * (t[t_ok] - mu) ** 2).sum()
        means[qi, k] = mu
        stds[qi, k] = np.sqrt(var)

imputed = {col: means[:, k] for k, col in enumerate(TARGET_COLS)}
imputed_std = {col: stds[:, k] for k, col in enumerate(TARGET_COLS)}

from common import assemble_and_save
out_path = os.path.join(SCRIPT_DIR, 'results', 'exp10_ot_errweighted.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp10_ot_errweighted: {len(out)} GRBs -> {out_path}')
