#!/usr/bin/env python3
"""
quick_ridge_check_spencer.py -- fast sanity gate for the "real recovery +
is_optical + Daume domain columns" idea, before committing hours to the full
spencer-logic SuperLearner run.

Replicates (in Python, ~seconds) the data-side part of spencer's pipeline:
  - real T90/Fluence/PhotonIndex/NH/PeakFlux recovered from Data/optical_data.csv
    for optical-only GRBs where available (same unit conversions as the R version)
  - is_optical flag
  - Daume-style scaled per-domain copies (value * is_optical * 0.25)
Then runs the same 10-fold RidgeCV check as quick_ridge_check.py, with and
without the added columns, so the two can be compared directly.

Usage: python3 quick_ridge_check_spencer.py <base_dataset.csv>
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
DAUME_C  = 0.25


def norm(g):
    g = str(g).replace('GRB', '').strip()
    return g if g[-1].isalpha() else g + 'A'


def run_ridge(X, y, label):
    model = make_pipeline(
        SimpleImputer(strategy='median'),
        StandardScaler(),
        RidgeCV(alphas=np.logspace(-3, 3, 25)),
    )
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    pred = cross_val_predict(model, X, y, cv=kf)
    r_log = np.corrcoef(y, pred)[0, 1]
    z_true, z_pred = 10 ** y - 1, 10 ** pred - 1
    r_z = np.corrcoef(z_true, z_pred)[0, 1]
    print(f'  [{label}] n={len(y)} feats={X.shape[1]}  log10(z+1) r={r_log:.3f}   z r={r_z:.3f}')
    return r_log, r_z


def main(path):
    df = pd.read_csv(path, index_col=0)
    print(f'{path}  ({len(df)} GRBs)')

    xray = pd.read_csv('Data/Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv')
    xray_ids = set(norm(g) for g in xray['GRB_Name'])
    is_opt = ~df.index.to_series().apply(norm).isin(xray_ids)
    print(f'  optical-only GRBs: {is_opt.sum()}')

    # ---- Baseline: 10 core features, whatever is already in the dataset ----
    y = np.log10(pd.to_numeric(df[RESPONSE], errors='coerce').values + 1)
    X_base = df[FEATURES].apply(pd.to_numeric, errors='coerce').values
    r_log_base, r_z_base = run_ridge(X_base, y, 'baseline (as-is)')

    # ---- Recovery: real values from optical_data.csv where available ----
    opt_cat = pd.read_csv('Data/optical_data.csv', index_col=0)
    opt_cat.index = opt_cat.index.to_series().apply(norm)
    opt_cat = opt_cat[~opt_cat.index.duplicated()]

    recovered = df.copy()
    opt_rows = recovered.index[is_opt]
    matched = opt_cat.reindex([norm(i) for i in opt_rows])
    matched.index = opt_rows

    recovered.loc[opt_rows, 'log10T90']      = np.log10(matched['T90']).values
    recovered.loc[opt_rows, 'log10Fluence']  = (np.log10(matched['Fluence']) - 7).values
    recovered.loc[opt_rows, 'PhotonIndex']   = matched['PhotonIndex'].values
    recovered.loc[opt_rows, 'log10NH']       = (np.log10(matched['NH']) + 21).values
    recovered.loc[opt_rows, 'log10PeakFlux'] = np.log10(matched['PeakFlux'].where(matched['PeakFlux'] > 0)).values
    # Gamma has no optical analog -- stays whatever was already in the dataset/NaN

    n_recovered = matched['T90'].notna().sum()
    print(f'  recovered real values for {n_recovered}/{is_opt.sum()} optical-only GRBs (T90); '
          f'partial for Fluence/PhotonIndex/NH/PeakFlux (~19-21/73 missing at source)')

    X_recovered = recovered[FEATURES].apply(pd.to_numeric, errors='coerce').values
    r_log_rec, r_z_rec = run_ridge(X_recovered, y, 'recovered (no Daume cols)')

    # ---- + is_optical + Daume-scaled domain copies ----
    is_opt_vec = is_opt.astype(float).values
    opt_cols = recovered[FEATURES].apply(pd.to_numeric, errors='coerce').values * is_opt_vec[:, None] * DAUME_C
    X_daume = np.hstack([X_recovered, is_opt_vec[:, None], opt_cols])
    r_log_daume, r_z_daume = run_ridge(X_daume, y, 'recovered + is_optical + Daume')

    print(f'\n  Summary for {path}:')
    print(f'    baseline                      log r={r_log_base:.3f}  z r={r_z_base:.3f}')
    print(f'    + real recovery                log r={r_log_rec:.3f}  z r={r_z_rec:.3f}')
    print(f'    + is_optical/Daume columns      log r={r_log_daume:.3f}  z r={r_z_daume:.3f}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python3 quick_ridge_check_spencer.py <base_dataset.csv>')
        sys.exit(1)
    main(sys.argv[1])
