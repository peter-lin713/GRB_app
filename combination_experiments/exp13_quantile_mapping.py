#!/usr/bin/env python3
"""exp13: rank/quantile distribution mapping. Instead of matching individual
GRBs by distance, collapse the 4 calibrated Dainotti params to a single
composite score (first principal component, weighted by anchor-informed
per-dimension trust), find each optical-only GRB's percentile rank on that
score relative to the X-ray population, then map that same percentile
through each target feature's own empirical distribution (inverse CDF).
A GRB near the top of the composite score gets a target value near the top
of the target's distribution -- preserves population-level shape rather
than local neighbor structure."""
import os
import numpy as np
from sklearn.preprocessing import StandardScaler
from common import load_data, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']

opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
scaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
diffs = scaler.transform(opt_anch[usable]) - scaler.transform(xray_anch[usable])
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()

xray_z = scaler.transform(xray[MATCH_COLS].values.astype(float)) * np.sqrt(weights)
opt_z  = scaler.transform(opt_only_proj[MATCH_COLS].values.astype(float)) * np.sqrt(weights)
composite_score = xray_z.sum(axis=1)  # weighted composite, X-ray reference population
opt_score = opt_z.sum(axis=1)

xray_ok = ~np.isnan(composite_score)
ref_scores = composite_score[xray_ok]
ranks = np.argsort(np.argsort(ref_scores)) / (len(ref_scores) - 1)  # 0..1 percentile per X-ray GRB
opt_percentile = np.array([np.mean(ref_scores <= s) for s in opt_score])

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    y = xray[col].values[xray_ok]
    y_ok = ~np.isnan(y)
    y_sorted = np.sort(y[y_ok])
    pred = np.interp(opt_percentile, np.linspace(0, 1, len(y_sorted)), y_sorted)
    imputed[col] = pred
    imputed_std[col] = [y_sorted.std()] * len(opt_only_proj)

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp13_quantile_mapping.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp13_quantile_mapping: {len(out)} GRBs -> {out_path}')
