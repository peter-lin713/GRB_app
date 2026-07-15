#!/usr/bin/env python3
"""
quick_ridge_check.py -- fast sanity check for a new training CSV.

Fits a cross-validated ridge regression (10-fold, RidgeCV over a log-spaced
alpha grid) predicting log10(z+1) from the 10 core features, and reports the
out-of-fold Pearson correlation. Meant as a cheap (~seconds, vs. ~30 minutes
for the full SuperLearner run) gate before committing to a full run on a
newly built dataset: median-imputes remaining NaNs rather than running MICE,
so it is not a substitute for the real pipeline -- just a smoke test of
whether the dataset has any signal at all.

Usage: python3 quick_ridge_check.py <path/to/dataset.csv>
"""

import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

FEATURES = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma',
            'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux']
RESPONSE = 'Redshift_crosscheck'

def main(path):
    df = pd.read_csv(path)
    missing = [c for c in FEATURES + [RESPONSE] if c not in df.columns]
    if missing:
        print(f'ERROR: missing columns {missing}')
        sys.exit(1)

    df = df[df[RESPONSE].notna()].copy()
    X = df[FEATURES].apply(pd.to_numeric, errors='coerce').values
    y = np.log10(pd.to_numeric(df[RESPONSE], errors='coerce').values + 1)

    n = len(df)
    print(f'{path}')
    print(f'n = {n} GRBs, {len(FEATURES)} features')
    nan_frac = np.isnan(X).mean(axis=0)
    for f, frac in zip(FEATURES, nan_frac):
        if frac > 0:
            print(f'  {f}: {frac:.1%} missing (median-imputed for this quick check)')

    model = make_pipeline(
        SimpleImputer(strategy='median'),
        StandardScaler(),
        RidgeCV(alphas=np.logspace(-3, 3, 25)),
    )
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    pred = cross_val_predict(model, X, y, cv=kf)

    r_log = np.corrcoef(y, pred)[0, 1]
    rmse_log = np.sqrt(np.mean((y - pred) ** 2))
    z_true = 10 ** y - 1
    z_pred = 10 ** pred - 1
    r_z = np.corrcoef(z_true, z_pred)[0, 1]

    print(f'\n10-fold CV ridge regression (quick check, median-imputed, no MICE):')
    print(f'  log10(z+1): r = {r_log:.3f}   RMSE = {rmse_log:.3f}')
    print(f'  z scale:    r = {r_z:.3f}')
    return r_log, r_z

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python3 quick_ridge_check.py <path/to/dataset.csv>')
        sys.exit(1)
    main(sys.argv[1])
