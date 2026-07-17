#!/usr/bin/env python3
"""exp11: proper multivariate chained-equations imputation (sklearn's
IterativeImputer, the real MICE algorithm -- BayesianRidge estimator,
cycles through each incomplete column using all others as predictors,
10 iterations) across the joint optical+X-ray table. Distinct from
exp02 (single-column median) and exp04/05 (single-pass regression):
MICE iterates and lets every column inform every other column, including
the 5 target columns informing each other."""
import os
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

xray_rows = xray[FEATURE_COLS + TARGET_COLS].copy()
opt_rows = opt_only_proj[FEATURE_COLS].copy()
for c in TARGET_COLS:
    opt_rows[c] = np.nan
combined = pd.concat([xray_rows, opt_rows[FEATURE_COLS + TARGET_COLS]])

imputer = IterativeImputer(estimator=BayesianRidge(), max_iter=15, sample_posterior=True, random_state=42)
filled = pd.DataFrame(imputer.fit_transform(combined), index=combined.index, columns=combined.columns)

opt_filled = filled.loc[opt_only_proj.index]
imputed = {col: opt_filled[col].values for col in TARGET_COLS}
imputed_std = {col: [xray[col].std()] * len(opt_only_proj) for col in TARGET_COLS}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp11_iterative_mice.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp11_iterative_mice: {len(out)} GRBs -> {out_path}')
