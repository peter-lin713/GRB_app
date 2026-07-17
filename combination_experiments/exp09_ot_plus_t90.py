#!/usr/bin/env python3
"""exp09: same OT kernel matching as exp08, but adds log10T90 to the
matching space. T90 is measured directly and identically in both catalogs
(unlike the calibrated Dainotti params, no emcee projection needed) and is
not the training target, so it's safe to use -- tests whether burst
duration carries information the 4 calibrated params miss."""
import os
import numpy as np
from sklearn.preprocessing import StandardScaler
from common import load_data, kernel_match_impute, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x', 'log10T90']

opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
scaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
diffs = scaler.transform(opt_anch[usable]) - scaler.transform(xray_anch[usable])
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()
print('Dimension weights:', dict(zip(MATCH_COLS, np.round(weights, 3))))

anchor_ids = np.array(overlap)[usable]
anchor_df = opt_proj.loc[anchor_ids]
best_reg, best_score = 0.2, -np.inf
for reg in [0.5, 0.3, 0.2, 0.15, 0.1, 0.07]:
    mu_all, _ = kernel_match_impute(anchor_df, xray, TARGET_COLS, MATCH_COLS, weights=weights, bandwidth=reg)
    corrs = []
    for col in TARGET_COLS:
        t = xray.loc[anchor_ids, col].values
        p = mu_all[col].values
        ok = ~(np.isnan(t) | np.isnan(p))
        if ok.sum() >= 3 and p[ok].std() > 1e-8:
            corrs.append(np.corrcoef(t[ok], p[ok])[0, 1])
    score = np.mean(corrs) if corrs else -1
    if score > best_score:
        best_score, best_reg = score, reg

mu, sd = kernel_match_impute(opt_only_proj, xray, TARGET_COLS, MATCH_COLS, weights=weights, bandwidth=best_reg)
imputed = {col: mu[col].values for col in TARGET_COLS}
imputed_std = {col: sd[col].values for col in TARGET_COLS}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp09_ot_plus_t90.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp09_ot_plus_t90: bandwidth={best_reg} (LOO score={best_score:.3f}), {len(out)} GRBs -> {out_path}')
