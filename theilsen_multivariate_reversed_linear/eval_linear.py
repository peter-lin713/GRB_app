#!/usr/bin/env python3
"""
Honest linear-model comparison of the four Theil-Sen combination schemes:
forward (optical -> x-ray) vs reversed (x-ray -> optical), single-variable vs
multivariate calibration.

Protocol = the "honest held-out r" of
GRB_linear_model_comparison/simple_linear_flagship_check.py (spencer_most_recent):
100 random 80/20 splits, OLS refit on the 80%, scored on the untouched 20%.
No M-estimator cut, no 2-sigma clip, no LASSO / formula search. Target is
log10(z+1); r is Pearson in LINEAR z (z = 10**pred - 1), plus r in log10(z+1)
and linear-z RMSE.

Each scheme is scored on its own pre-M-estimator frame (grb_xray_imputed.csv,
built by the same prep code) and again on the GRBs common to all four frames,
using identical splits, so frame membership can't explain a difference.

Run from this folder: python3 eval_linear.py
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GAM = os.path.join(HERE, '..', 'GAM_supercomputer')
FRAMES = {
    'fwd single (0.707 family)': os.path.join(GAM, 'theilsen_cleaned_formula_generation', 'grb_xray_imputed.csv'),
    'fwd multivariate':          os.path.join(GAM, 'multivariate_calibration_theilsen_formula_generation', 'grb_xray_imputed.csv'),
    'rev single':                os.path.join(GAM, 'theilsen_reversed_xray_to_opt_formula_generation', 'grb_xray_imputed.csv'),
    'rev multivariate (NEW)':    os.path.join(HERE, 'frame', 'grb_xray_imputed.csv'),
}
MODELS = {
    'SIMPLE7': ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'log10NH', 'PhotonIndex', 'log10PeakFlux'],
    'ALL10':   ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma',
                'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux'],
}
N_SPLITS, SEED = 100, 1


def honest(d, feats, n_splits=N_SPLITS, seed=SEED):
    y = d['log10z'].values
    X = np.column_stack([np.ones(len(d)), d[feats].values])
    rng, n = np.random.default_rng(seed), len(y)
    r_lin, r_log, rmse = [], [], []
    for _ in range(n_splits):
        te = np.zeros(n, bool)
        te[rng.choice(n, size=int(0.20 * n), replace=False)] = True
        beta, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        p = X[te] @ beta
        zp, zo = 10 ** p - 1, 10 ** y[te] - 1
        r_lin.append(np.corrcoef(zp, zo)[0, 1])
        r_log.append(np.corrcoef(p, y[te])[0, 1])
        rmse.append(np.sqrt(np.mean((zp - zo) ** 2)))
    return dict(N=n, r_z=np.mean(r_lin), r_z_sd=np.std(r_lin),
                r_log=np.mean(r_log), rmse_z=np.mean(rmse))


frames = {k: pd.read_csv(v, index_col=0) for k, v in FRAMES.items()}
common = sorted(set.intersection(*(set(f.index) for f in frames.values())))

rows = []
for scope in ('own frame', 'common GRBs'):
    for name, d in frames.items():
        dd = d if scope == 'own frame' else d.loc[common]
        for m, feats in MODELS.items():
            rows.append(dict(scope=scope, scheme=name, model=m, **honest(dd, feats)))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(HERE, 'results_honest_linear.csv'), index=False, float_format='%.4f')

pd.set_option('display.width', 140)
for scope in ('own frame', 'common GRBs'):
    print(f'\n=== {scope}' + (f' (N={len(common)})' if scope == 'common GRBs' else '') + ' ===')
    t = res[res.scope == scope].drop(columns='scope')
    print(t.to_string(index=False, float_format=lambda v: f'{v:.3f}'))

# Split-to-split noise check: paired difference NEW - rev single on common
# GRBs, same splits, so the comparison isn't swamped by split variance.
print('\nPaired difference on common GRBs (same 100 splits), r in linear z:')
for m, feats in MODELS.items():
    a = frames['rev multivariate (NEW)'].loc[common]
    for ref in ['rev single', 'fwd single (0.707 family)', 'fwd multivariate']:
        b = frames[ref].loc[common]
        da, db = [], []
        rng = np.random.default_rng(SEED)
        n = len(common)
        ya = a['log10z'].values
        Xa = np.column_stack([np.ones(n), a[feats].values])
        Xb = np.column_stack([np.ones(n), b[feats].values])
        for _ in range(N_SPLITS):
            te = np.zeros(n, bool)
            te[rng.choice(n, size=int(0.20 * n), replace=False)] = True
            for X, acc in ((Xa, da), (Xb, db)):
                beta, *_ = np.linalg.lstsq(X[~te], ya[~te], rcond=None)
                acc.append(np.corrcoef(10 ** (X[te] @ beta) - 1, 10 ** ya[te] - 1)[0, 1])
        dd = np.array(da) - np.array(db)
        print(f'  {m:7s} NEW - {ref:26s}: {dd.mean():+.3f}  '
              f'(sd {dd.std():.3f}; NEW better on {np.mean(dd > 0):.0%} of splits)')
