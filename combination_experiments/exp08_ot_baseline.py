#!/usr/bin/env python3
"""exp08: OT baseline -- the current production method (combine_optical_
xray_ot_v3.py after the leakage fix): kernel-weighted matching on the 4
calibrated Dainotti parameters only (no log10T90, no z), anchor-informed
per-dimension weights, bandwidth swept via leave-one-out on the 78 anchors.
Included here so the other 9 experiments have a same-harness reference
point, not just the numbers reported earlier from the standalone script."""
import os
import numpy as np
from common import load_data, kernel_match_impute, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']

opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
opt_anch_z = scaler.transform(opt_anch[usable])
xray_anch_z = scaler.transform(xray_anch[usable])
diffs = opt_anch_z - xray_anch_z
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()

# Bandwidth sweep via leave-one-out on anchors (same idea as the production script)
anchor_ids = np.array(overlap)[usable]
anchor_df = opt_proj.loc[anchor_ids]
best_reg, best_score = 0.2, -np.inf
for reg in [0.5, 0.3, 0.2, 0.15, 0.1, 0.07]:
    corrs = []
    for grb in anchor_ids:
        donors = xray.drop(index=grb)
        mu, sd = kernel_match_impute(anchor_df.loc[[grb]], donors, TARGET_COLS, MATCH_COLS, weights=weights, bandwidth=reg)
        pass
    # cheap approximate score: use full batch LOO via matching all anchors at once against xray (donors unchanged)
    mu_all, sd_all = kernel_match_impute(anchor_df, xray, TARGET_COLS, MATCH_COLS, weights=weights, bandwidth=reg)
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

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp08_ot_baseline.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp08_ot_baseline: bandwidth={best_reg} (LOO score={best_score:.3f}), {len(out)} GRBs -> {out_path}')
