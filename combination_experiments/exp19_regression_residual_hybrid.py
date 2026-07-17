#!/usr/bin/env python3
"""exp19: regression + residual matching, a "kriging"-style hybrid. First
fit a global ridge regression per target (same as exp05) to capture the
broad trend. Then, instead of trusting that global prediction alone,
OT-match each optical-only GRB against the X-ray donors' RESIDUALS from
that same regression (donor_residual = donor_true - donor_regression_
prediction) and add the matched local residual correction back on top.
Combines a stable global trend with a local correction the trend misses --
distinct from pure OT (no global trend) and pure regression (no local
correction)."""
import os
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from common import load_data, kernel_match_impute, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']

X_train_raw = xray[FEATURE_COLS].values
X_query_raw = opt_only_proj[FEATURE_COLS].values
imp = SimpleImputer(strategy='median').fit(X_train_raw)
scaler = StandardScaler().fit(imp.transform(X_train_raw))
X_train = scaler.transform(imp.transform(X_train_raw))
X_query = scaler.transform(imp.transform(X_query_raw))

opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
mscaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
diffs = mscaler.transform(opt_anch[usable]) - mscaler.transform(xray_anch[usable])
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()

imputed, imputed_std = {}, {}
xray_resid = xray.copy()
for col in TARGET_COLS:
    y = xray[col].values
    ok = ~np.isnan(y)
    model = RidgeCV(alphas=np.logspace(-3, 3, 25)).fit(X_train[ok], y[ok])
    global_pred_train = np.full(len(y), np.nan)
    global_pred_train[ok] = model.predict(X_train[ok])
    xray_resid[col + '_resid'] = y - global_pred_train  # NaN where y is NaN
    global_pred_query = model.predict(X_query)
    imputed[col + '_global'] = global_pred_query

resid_cols = [c + '_resid' for c in TARGET_COLS]
mu_resid, sd_resid = kernel_match_impute(opt_only_proj, xray_resid, resid_cols, MATCH_COLS, weights=weights, bandwidth=0.15)

for col in TARGET_COLS:
    imputed[col] = imputed.pop(col + '_global') + mu_resid[col + '_resid'].values
    imputed_std[col] = sd_resid[col + '_resid'].values

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp19_regression_residual_hybrid.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp19_regression_residual_hybrid: {len(out)} GRBs -> {out_path}')
