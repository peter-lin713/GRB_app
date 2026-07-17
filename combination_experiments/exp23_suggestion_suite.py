#!/usr/bin/env python3
"""exp23: suite covering the deep-research suggestions + combinations,
all gated with the round2 ridge check, reporting r on ALL rows and on
the X-RAY subset only (the fair comparison vs exp01, see exp22).

Suggestion coverage:
  1. Hierarchical Bayesian calibration -- common.load_data()'s emcee
     projection IS this (per-param MCMC, x/y errors + intrinsic scatter);
     it underlies every dataset here. 'bayes_full' additionally gates the
     emcee_v2 production dataset, where the prompt features come from the
     optical catalog itself instead of being imputed.
  2. Deming errors-in-variables regression -- 'deming': per-param Deming
     fit on the 78 anchors replaces the emcee projection; prompt features
     via the shared kernel-match imputation.
  3. OT + bootstrap multiple imputation -- 'ot_boot_mi': B=20 reps of
     kernel-match imputation with donor resampling + error perturbation;
     pooled mean, total (within+between) variance.
  4. Heteroscedastic weighting -- applied as a combo to each dataset:
     w=1 (plain pooling), w=0.3, and inverse-variance per-row weights.
  5. Semi-supervised pseudo-labeling -- 'pseudolabel': train on X-ray,
     pseudo-label the projected rows, retrain with them at w=0.3.
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

from common import (load_data, assemble_and_save, kernel_match_impute,
                    TARGET_COLS, FEATURE_COLS, SCRIPT_DIR, DATA_DIR, MIN_ERR)

FEATURES = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'Beta', 'Gamma',
            'log10Fluence', 'PhotonIndex', 'log10NH', 'log10PeakFlux']
ERR_COLS = ['log10FaErr', 'log10TaErr', 'AlphaErr', 'BetaErr']
RESPONSE = 'Redshift_crosscheck'
RESULTS = os.path.join(SCRIPT_DIR, 'results')

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- deming
params_map = {   # optical col, optical err  ->  xray col, xray err
    'logFa': ('log10Faopt', 'log10FaErr_opt', 'logFa_x', 'logFaErr_x'),
    'logTa': ('log10Taopt', 'logTaErr_opt',   'logTa_x', 'logTaErr_x'),
    'Alpha': ('Alpha_opt',  'AlphaErr_opt',   'Alpha_x', 'AlphaErr_x'),
    'Beta':  ('Beta_opt',   'betaErr_opt',    'Beta_x',  'BetaErr_x'),
}
proj_out = {'logFa': ('logFa_x', 'logFaErr_x'), 'logTa': ('logTa_x', 'logTaErr_x'),
            'Alpha': ('Alpha_x', 'AlphaErr_x'), 'Beta':  ('Beta_x',  'BetaErr_x')}

deming_proj = pd.DataFrame({'z': opt_clean['z'], 'log10T90': opt_clean['log10T90']},
                           index=opt_clean.index)
for p, (oc, oe, xc, xe) in params_map.items():
    a = pd.DataFrame({
        'x': opt_clean.loc[overlap, oc], 'xe': opt_clean.loc[overlap, oe],
        'y': xray.loc[overlap, xc],      'ye': xray.loc[overlap, xe],
    }).apply(pd.to_numeric, errors='coerce').dropna()
    delta = np.maximum((a['ye'] ** 2).mean(), MIN_ERR) / np.maximum((a['xe'] ** 2).mean(), MIN_ERR)
    xb, yb = a['x'].mean(), a['y'].mean()
    sxx = ((a['x'] - xb) ** 2).mean()
    syy = ((a['y'] - yb) ** 2).mean()
    sxy = ((a['x'] - xb) * (a['y'] - yb)).mean()
    m = ((syy - delta * sxx) + np.sqrt((syy - delta * sxx) ** 2 + 4 * delta * sxy ** 2)) / (2 * sxy)
    b = yb - m * xb
    resid_sd = (a['y'] - (m * a['x'] + b)).std()
    yc, ye_out = proj_out[p]
    xv = opt_clean[oc].astype(float)
    xerr = opt_clean[oe].astype(float).fillna(MIN_ERR).clip(lower=MIN_ERR)
    deming_proj[yc] = m * xv + b
    deming_proj[ye_out] = np.sqrt((m * xerr) ** 2 + resid_sd ** 2)
    print(f'deming {p}: m={m:.3f} b={b:.3f} resid_sd={resid_sd:.3f} (N={len(a)})')

deming_only = deming_proj[~deming_proj.index.isin(overlap)].copy()
mu, sd = kernel_match_impute(deming_only, xray, TARGET_COLS, FEATURE_COLS)
assemble_and_save(xray, opt_clean, deming_proj, deming_only,
                  {c: mu[c].values for c in TARGET_COLS},
                  {c: sd[c].values for c in TARGET_COLS},
                  os.path.join(RESULTS, 'exp23_deming.csv'))

# ------------------------------------------------------------ ot_boot_mi
B = 20
boot_means, boot_vars = [], []
err_map = {'logFa_x': 'logFaErr_x', 'logTa_x': 'logTaErr_x',
           'Alpha_x': 'AlphaErr_x', 'Beta_x': 'BetaErr_x'}
for bi in range(B):
    donors = xray.sample(n=len(xray), replace=True, random_state=bi)
    donors = donors[~donors.index.duplicated()]
    q = opt_only_proj.copy()
    for fc, ec in err_map.items():
        e = q[ec].astype(float).fillna(MIN_ERR).clip(lower=MIN_ERR)
        q[fc] = q[fc].astype(float) + rng.standard_normal(len(q)) * e
    mu_b, sd_b = kernel_match_impute(q, donors, TARGET_COLS, FEATURE_COLS)
    boot_means.append(mu_b)
    boot_vars.append(sd_b ** 2)
mi_mean = sum(boot_means) / B
within = sum(boot_vars) / B
between = sum((m - mi_mean) ** 2 for m in boot_means) / (B - 1)
mi_sd = np.sqrt(within + (1 + 1 / B) * between)
assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj,
                  {c: mi_mean[c].values for c in TARGET_COLS},
                  {c: mi_sd[c].values for c in TARGET_COLS},
                  os.path.join(RESULTS, 'exp23_ot_boot_mi.csv'))
print(f'ot_boot_mi: B={B}, between-var share='
      f'{(between.values.mean() / (within.values.mean() + between.values.mean())):.2f}')

# --------------------------------------------------------------- gating
def load_ds(path):
    df = pd.read_csv(path, index_col=0)
    df = df[pd.to_numeric(df[RESPONSE], errors='coerce').notna()].copy()
    X = df[FEATURES].apply(pd.to_numeric, errors='coerce').values
    y = np.log10(pd.to_numeric(df[RESPONSE], errors='coerce').values + 1)
    is_x = df.index.isin(xray.index)
    err = df[ERR_COLS].apply(pd.to_numeric, errors='coerce')
    row_var = (err ** 2).mean(axis=1)
    ref = row_var[is_x].median()
    invvar = np.clip(ref / row_var.values, 0, 1)
    invvar[is_x] = 1.0
    return X, y, is_x, invvar

def ridge():
    return make_pipeline(SimpleImputer(strategy='median'), StandardScaler(),
                         RidgeCV(alphas=np.logspace(-3, 3, 25)))

def gate(X, y, is_x, weights, seeds=10):
    r_all, r_x = [], []
    for seed in range(seeds):
        pred = np.full(len(y), np.nan)
        for tr, te in KFold(10, shuffle=True, random_state=seed).split(X):
            m = ridge()
            m.fit(X[tr], y[tr], ridgecv__sample_weight=weights[tr])
            pred[te] = m.predict(X[te])
        r_all.append(np.corrcoef(y, pred)[0, 1])
        r_x.append(np.corrcoef(y[is_x], pred[is_x])[0, 1])
    return np.mean(r_all), np.mean(r_x), np.std(r_x)

def gate_pseudolabel(X, y, is_x, w_pseudo=0.3, seeds=10):
    """Train on x-ray fold, pseudo-label projected rows, retrain, score x-ray."""
    xi = np.where(is_x)[0]
    pi = np.where(~is_x)[0]
    r_x = []
    for seed in range(seeds):
        pred = np.full(len(y), np.nan)
        for tr_l, te_l in KFold(10, shuffle=True, random_state=seed).split(xi):
            tr, te = xi[tr_l], xi[te_l]
            m1 = ridge()
            m1.fit(X[tr], y[tr])
            y_pseudo = m1.predict(X[pi])
            X2 = np.vstack([X[tr], X[pi]])
            y2 = np.concatenate([y[tr], y_pseudo])
            w2 = np.concatenate([np.ones(len(tr)), np.full(len(pi), w_pseudo)])
            m2 = ridge()
            m2.fit(X2, y2, ridgecv__sample_weight=w2)
            pred[te] = m2.predict(X[te])
        r_x.append(np.corrcoef(y[xi], pred[xi])[0, 1])
    return np.nan, np.mean(r_x), np.std(r_x)

datasets = {
    'ot_baseline (exp08)': os.path.join(RESULTS, 'exp08_ot_baseline.csv'),
    'deming':              os.path.join(RESULTS, 'exp23_deming.csv'),
    'ot_boot_mi':          os.path.join(RESULTS, 'exp23_ot_boot_mi.csv'),
    'bayes_full (emceev2)': os.path.join(DATA_DIR, 'superlearner_training_emcee_v2_errcut_relative.csv'),
}

rows = []
for name, path in datasets.items():
    X, y, is_x, invvar = load_ds(path)
    for scheme in ['w=1', 'w=0.3', 'invvar']:
        w = (np.ones(len(y)) if scheme == 'w=1'
             else np.where(is_x, 1.0, 0.3) if scheme == 'w=0.3'
             else invvar)
        ra, rx, sx = gate(X, y, is_x, w)
        rows.append((name, scheme, len(y), ra, rx, sx))
        print(f'{name:22} {scheme:7} n={len(y)}  r_all={ra:.3f}  r_xray={rx:.3f}+-{sx:.3f}')

X, y, is_x, _ = load_ds(datasets['ot_baseline (exp08)'])
ra, rx, sx = gate_pseudolabel(X, y, is_x)
rows.append(('pseudolabel (on exp08)', 'w=0.3', len(y), ra, rx, sx))
print(f'{"pseudolabel (on exp08)":22} w=0.3   n={len(y)}  r_all=  -    r_xray={rx:.3f}+-{sx:.3f}')

out = pd.DataFrame(rows, columns=['dataset', 'scheme', 'n', 'r_all', 'r_xray', 'std_xray'])
out.to_csv(os.path.join(RESULTS, 'exp23_suggestion_suite_summary.csv'), index=False)
print('wrote results/exp23_suggestion_suite_summary.csv')
