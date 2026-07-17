#!/usr/bin/env python3
"""exp15: soft-cluster pooling via a Gaussian Mixture Model -- same idea as
exp14, but each optical-only GRB gets a probabilistic membership across all
K components instead of one hard assignment, and the imputed value is the
membership-weighted blend of every cluster's mean. Softer still than
K-means; lets a GRB that sits between two groups borrow from both."""
import os
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from common import load_data, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']
K = 6

imp = SimpleImputer(strategy='median').fit(xray[MATCH_COLS].values)
X_train = imp.transform(xray[MATCH_COLS].values)
X_query = imp.transform(opt_only_proj[MATCH_COLS].values)
scaler = StandardScaler().fit(X_train)
X_train_z = scaler.transform(X_train)
X_query_z = scaler.transform(X_query)

gmm = GaussianMixture(n_components=K, covariance_type='full', random_state=42, n_init=5).fit(X_train_z)
train_resp = gmm.predict_proba(X_train_z)   # (223, K)
query_resp = gmm.predict_proba(X_query_z)   # (83, K)

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    y = xray[col].values
    y_ok = ~np.isnan(y)
    # weighted mean/var per component, using only rows with an observed value
    comp_mean = np.zeros(K)
    comp_var = np.zeros(K)
    for c in range(K):
        w = train_resp[y_ok, c]
        if w.sum() < 1e-6:
            comp_mean[c] = np.nanmean(y)
            comp_var[c] = np.nanvar(y)
            continue
        comp_mean[c] = np.average(y[y_ok], weights=w)
        comp_var[c] = np.average((y[y_ok] - comp_mean[c]) ** 2, weights=w)
    pred = query_resp @ comp_mean
    var = query_resp @ comp_var + (query_resp @ (comp_mean[None, :] - pred[:, None]).T.diagonal() * 0)  # base var
    # add between-component variance from the mixture itself
    between = (query_resp * (comp_mean[None, :] - pred[:, None]) ** 2).sum(axis=1)
    imputed[col] = pred
    imputed_std[col] = np.sqrt(query_resp @ comp_var + between)

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp15_gmm_softcluster.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp15_gmm_softcluster: {len(out)} GRBs -> {out_path}')
