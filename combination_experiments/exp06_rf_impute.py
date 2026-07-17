#!/usr/bin/env python3
"""exp06: random forest regressor per target feature, fit on the 223 X-ray
GRBs, predict for the 83 optical-only GRBs. Captures nonlinear relationships
the linear/ridge variants can't, while still using all X-ray GRBs jointly
(unlike the purely-local OT matching)."""
import os
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

X_train_raw = xray[FEATURE_COLS].values
X_query_raw = opt_only_proj[FEATURE_COLS].values
imp = SimpleImputer(strategy='median').fit(X_train_raw)
X_train = imp.transform(X_train_raw)
X_query = imp.transform(X_query_raw)

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    y = xray[col].values
    ok = ~np.isnan(y)
    model = RandomForestRegressor(n_estimators=300, max_depth=4, random_state=42, min_samples_leaf=5)
    model.fit(X_train[ok], y[ok])
    pred = model.predict(X_query)
    # per-tree spread at each query point as the uncertainty estimate
    tree_preds = np.stack([t.predict(X_query) for t in model.estimators_], axis=0)
    imputed[col] = pred
    imputed_std[col] = tree_preds.std(axis=0)

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp06_rf_impute.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp06_rf_impute: {len(out)} GRBs -> {out_path}')
