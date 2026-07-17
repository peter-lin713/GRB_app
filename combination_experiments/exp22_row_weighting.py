#!/usr/bin/env python3
"""exp22: heteroscedastic row weighting -- don't decide IF projected GRBs
enter training, decide HOW MUCH. X-ray rows get sample_weight 1; the 73
optical-only (projected+imputed) rows get either a swept global weight
w in [0,1] or a per-row inverse-variance weight from their projection
errors. Evaluated like round2 (RidgeCV, 10-fold CV x 10 seeds) but r is
reported separately on ALL rows and on the X-RAY SUBSET ONLY -- the
latter is the fair comparison against exp01 (does adding weighted
optical GRBs help or hurt prediction of the X-ray population?).

Uses exp08_ot_baseline.csv as the fused dataset (rerun exp08 first if
missing). w=0 reproduces exp01-like training; w=1 reproduces exp08."""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FEATURES = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma',
            'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux']
ERR_COLS = ['log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr']
RESPONSE = 'Redshift_crosscheck'

fused = pd.read_csv(os.path.join(SCRIPT_DIR, 'results', 'exp08_ot_baseline.csv'), index_col=0)
xray_only = pd.read_csv(os.path.join(SCRIPT_DIR, 'results', 'exp01_xray_only.csv'), index_col=0)

fused = fused[fused[RESPONSE].notna()].copy()
is_xray = fused.index.isin(xray_only.index)
print(f'fused: {len(fused)} rows ({is_xray.sum()} x-ray, {(~is_xray).sum()} projected)')

X = fused[FEATURES].apply(pd.to_numeric, errors='coerce').values
y = np.log10(pd.to_numeric(fused[RESPONSE], errors='coerce').values + 1)

# Per-row inverse-variance weights from the projected plateau errors,
# normalized so the x-ray rows' mean error gives weight 1 (capped at 1).
err = fused[ERR_COLS].apply(pd.to_numeric, errors='coerce')
row_var = (err ** 2).mean(axis=1)
ref_var = row_var[is_xray].median()
invvar_w = np.clip(ref_var / row_var.values, 0, 1)
invvar_w[is_xray] = 1.0
print(f'inv-var weights on projected rows: median={np.median(invvar_w[~is_xray]):.3f}, '
      f'range {invvar_w[~is_xray].min():.3f}-{invvar_w[~is_xray].max():.3f}')


def evaluate(weights, seeds=10):
    """10-fold CV x seeds with sample-weighted RidgeCV; returns
    (mean r all rows, mean r x-ray rows)."""
    r_all, r_xray = [], []
    for seed in range(seeds):
        pred = np.full(len(y), np.nan)
        kf = KFold(n_splits=10, shuffle=True, random_state=seed)
        for tr, te in kf.split(X):
            model = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(),
                                  RidgeCV(alphas=np.logspace(-3, 3, 25)))
            model.fit(X[tr], y[tr], ridgecv__sample_weight=weights[tr])
            pred[te] = model.predict(X[te])
        r_all.append(np.corrcoef(y, pred)[0, 1])
        r_xray.append(np.corrcoef(y[is_xray], pred[is_xray])[0, 1])
    return np.mean(r_all), np.std(r_all), np.mean(r_xray), np.std(r_xray)


rows = []
for w in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    weights = np.where(is_xray, 1.0, w)
    ra, sa, rx, sx = evaluate(weights)
    rows.append(('global_w=%.1f' % w, ra, sa, rx, sx))
    print(f'w={w:.1f}  r_all={ra:.3f}+-{sa:.3f}  r_xray={rx:.3f}+-{sx:.3f}')

ra, sa, rx, sx = evaluate(invvar_w)
rows.append(('invvar', ra, sa, rx, sx))
print(f'invvar  r_all={ra:.3f}+-{sa:.3f}  r_xray={rx:.3f}+-{sx:.3f}')

out = pd.DataFrame(rows, columns=['scheme', 'r_all', 'std_all', 'r_xray', 'std_xray'])
out_path = os.path.join(SCRIPT_DIR, 'results', 'exp22_row_weighting_summary.csv')
out.to_csv(out_path, index=False)
print(f'exp22_row_weighting: wrote {out_path}')
