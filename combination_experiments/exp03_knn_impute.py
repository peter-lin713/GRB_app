#!/usr/bin/env python3
"""exp03: off-the-shelf scikit-learn KNNImputer (k=5) on the combined
feature+target matrix -- standard hard-nearest-neighbor-average imputation,
no custom weighting or bandwidth tuning. A simple, well-known alternative
to the custom OT kernel matching."""
import os
import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

xray_rows = xray[FEATURE_COLS + TARGET_COLS].copy()
opt_rows  = opt_only_proj[FEATURE_COLS].copy()
for c in TARGET_COLS:
    opt_rows[c] = np.nan
combined = pd.concat([xray_rows, opt_rows[FEATURE_COLS + TARGET_COLS]])

imputer = KNNImputer(n_neighbors=5, weights='distance')
scaled = (combined - combined.mean()) / combined.std()
filled_scaled = pd.DataFrame(imputer.fit_transform(scaled), index=combined.index, columns=combined.columns)
filled = filled_scaled * combined.std() + combined.mean()

opt_filled = filled.loc[opt_only_proj.index]
imputed = {col: opt_filled[col].values for col in TARGET_COLS}
imputed_std = {col: [xray[col].std()] * len(opt_only_proj) for col in TARGET_COLS}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp03_knn_impute.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp03_knn_impute: {len(out)} GRBs -> {out_path}')
