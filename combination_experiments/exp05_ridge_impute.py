#!/usr/bin/env python3
"""exp05: same as exp04 but ridge regression (regularized) with alpha chosen
by cross-validation on the X-ray GRBs -- guards against overfitting the
global relationship on only 223 training rows."""
import os
import numpy as np
from sklearn.linear_model import RidgeCV
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

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    y = xray[col].values
    ok = ~np.isnan(y)
    model = RidgeCV(alphas=np.logspace(-3, 3, 25)).fit(X_train[ok], y[ok])
    pred = model.predict(X_query)
    resid_std = np.std(y[ok] - model.predict(X_train[ok]))
    imputed[col] = pred
    imputed_std[col] = [resid_std] * len(opt_only_proj)

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp05_ridge_impute.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp05_ridge_impute: {len(out)} GRBs -> {out_path}')
