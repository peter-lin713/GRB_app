#!/usr/bin/env python3
"""exp04: fit one plain linear regression per target feature on the 223
X-ray GRBs (X = the 4 calibrated Dainotti params + log10T90), apply to the
83 optical-only GRBs. A global-relationship alternative to local matching:
assumes each X-ray-only feature varies smoothly/linearly with the shared
parameters across the whole population, rather than borrowing from
similar neighbors."""
import os
import numpy as np
from sklearn.linear_model import LinearRegression
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
    model = LinearRegression().fit(X_train[ok], y[ok])
    pred = model.predict(X_query)
    resid_std = np.std(y[ok] - model.predict(X_train[ok]))
    imputed[col] = pred
    imputed_std[col] = [resid_std] * len(opt_only_proj)

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp04_linreg_impute.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp04_linreg_impute: {len(out)} GRBs -> {out_path}')
