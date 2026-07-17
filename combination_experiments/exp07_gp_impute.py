#!/usr/bin/env python3
"""exp07: Gaussian Process regression per target feature -- a smooth
nonlinear model class that natively returns a calibrated posterior std,
so uncertainty comes directly from the model rather than a workaround
(tree spread, donor-pool variance, etc)."""
import os
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

X_train_raw = xray[FEATURE_COLS].values
X_query_raw = opt_only_proj[FEATURE_COLS].values
imp = SimpleImputer(strategy='median').fit(X_train_raw)
scaler = StandardScaler().fit(imp.transform(X_train_raw))
X_train = scaler.transform(imp.transform(X_train_raw))
X_query = scaler.transform(imp.transform(X_query_raw))

kernel = ConstantKernel(1.0) * RBF(length_scale=np.ones(X_train.shape[1])) + WhiteKernel(noise_level=1.0)

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    y = xray[col].values
    ok = ~np.isnan(y)
    y_mean, y_std = y[ok].mean(), y[ok].std()
    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False, n_restarts_optimizer=2, random_state=42)
    gp.fit(X_train[ok], (y[ok] - y_mean) / y_std)
    pred_mean, pred_std = gp.predict(X_query, return_std=True)
    imputed[col] = pred_mean * y_std + y_mean
    imputed_std[col] = pred_std * y_std

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp07_gp_impute.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp07_gp_impute: {len(out)} GRBs -> {out_path}')
