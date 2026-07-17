#!/usr/bin/env python3
"""exp20: stacked ensemble -- average the three best-performing methods from
the first round (OT baseline, random forest, and plain median), recomputed
here inline. Ensembling independent methods often reduces variance even
when no single method is clearly best; tests whether that holds here."""
import os
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from common import load_data, kernel_match_impute, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']

# -- Method 1: median --
median_pred = {col: np.full(len(opt_only_proj), xray[col].median()) for col in TARGET_COLS}

# -- Method 2: OT baseline --
opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
scaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
diffs = scaler.transform(opt_anch[usable]) - scaler.transform(xray_anch[usable])
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()
mu_ot, sd_ot = kernel_match_impute(opt_only_proj, xray, TARGET_COLS, MATCH_COLS, weights=weights, bandwidth=0.1)

# -- Method 3: random forest --
X_train_raw = xray[FEATURE_COLS].values
X_query_raw = opt_only_proj[FEATURE_COLS].values
imp = SimpleImputer(strategy='median').fit(X_train_raw)
X_train = imp.transform(X_train_raw)
X_query = imp.transform(X_query_raw)
rf_pred = {}
for col in TARGET_COLS:
    y = xray[col].values
    ok = ~np.isnan(y)
    model = RandomForestRegressor(n_estimators=300, max_depth=4, random_state=42, min_samples_leaf=5)
    model.fit(X_train[ok], y[ok])
    rf_pred[col] = model.predict(X_query)

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    stacked = np.vstack([median_pred[col], mu_ot[col].values, rf_pred[col]])
    imputed[col] = stacked.mean(axis=0)
    imputed_std[col] = stacked.std(axis=0) + xray[col].std() * 0.3  # ensemble spread + a floor

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp20_stacked_ensemble.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp20_stacked_ensemble: {len(out)} GRBs -> {out_path}')
