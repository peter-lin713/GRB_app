#!/usr/bin/env python3
"""exp21: CCA-reduced OT matching. The existing OT baseline (exp08) weights
each RAW matching dimension independently (1/variance of the optical-vs-
X-ray residual per dimension) -- it can't capture that some *combination*
of the 4 shared params might be more predictive of the target features than
any single axis. Canonical Correlation Analysis (CCA) instead learns the
linear combination(s) of the shared params that are maximally correlated
with the target features directly, using the 78 anchors (the only GRBs
where we know both sides). Matching then happens in this small (2-D)
canonical space instead of the raw 4-D standardized space.

LEAKAGE NOTE: CCA's fit uses the anchors' TRUE target values (Y side) to
learn the projection -- unlike the other experiments, simply excluding a
held-out anchor from the *donor pool* is not enough for a fair leave-one-out
test, since the anchor's own true values would still have shaped the CCA
projection used to match it. The LOO loop below refits CCA from scratch
for every held-out anchor, dropping it from the CCA fit AND the donor pool.
This is the correct (if slower) way to validate a supervised-projection
method without leaking the label into its own evaluation.
"""
import os
import numpy as np
import pandas as pd
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
from common import load_data, kernel_match_impute, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']
N_COMPONENTS = 2

overlap = np.array(overlap)
X_anch_raw = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
Y_anch_raw = xray.loc[overlap, TARGET_COLS].values.astype(float)
usable = ~(np.isnan(X_anch_raw).any(axis=1) | np.isnan(Y_anch_raw).any(axis=1))
overlap_u = overlap[usable]
print(f'Usable anchors for CCA fit: {usable.sum()} / {len(overlap)}')


def fit_cca_and_project(anchor_ids, donor_df, query_df):
    """Fit CCA on the given anchors (X=calibrated shared params, Y=true
    target features), then project donor_df and query_df's shared params
    into the learned canonical X-space. Returns (donor_scores, query_scores).
    """
    Xa = opt_proj.loc[anchor_ids, MATCH_COLS].values.astype(float)
    Ya = xray.loc[anchor_ids, TARGET_COLS].values.astype(float)
    x_scaler = StandardScaler().fit(Xa)
    y_scaler = StandardScaler().fit(Ya)
    Xa_z = x_scaler.transform(Xa)
    Ya_z = y_scaler.transform(Ya)

    n_comp = min(N_COMPONENTS, Xa_z.shape[0] - 1, Xa_z.shape[1], Ya_z.shape[1])
    cca = CCA(n_components=n_comp, max_iter=2000)
    cca.fit(Xa_z, Ya_z)

    donor_X = x_scaler.transform(donor_df[MATCH_COLS].values.astype(float))
    query_X = x_scaler.transform(query_df[MATCH_COLS].values.astype(float))
    donor_scores = cca.transform(donor_X)
    query_scores = cca.transform(query_X)
    return donor_scores, query_scores, cca


# ---- Fit CCA on all usable anchors, report canonical correlations ----
_, _, cca_full = fit_cca_and_project(overlap_u, xray, opt_only_proj)
Xa_z = StandardScaler().fit_transform(opt_proj.loc[overlap_u, MATCH_COLS].values.astype(float))
Ya_z = StandardScaler().fit_transform(xray.loc[overlap_u, TARGET_COLS].values.astype(float))
Xc, Yc = cca_full.transform(Xa_z, Ya_z)
for k in range(Xc.shape[1]):
    print(f'  canonical component {k+1}: corr(X_c, Y_c) = {np.corrcoef(Xc[:,k], Yc[:,k])[0,1]:.3f}')

# ---- Leave-one-out validation: refit CCA excluding each held-out anchor ----
print('\nLeave-one-out validation (CCA refit per anchor, no leakage):')
loo_true, loo_pred = {c: [] for c in TARGET_COLS}, {c: [] for c in TARGET_COLS}
for grb in overlap_u:
    train_anchors = overlap_u[overlap_u != grb]
    donors = xray.drop(index=grb)
    donor_scores, query_scores, _ = fit_cca_and_project(train_anchors, donors, opt_proj.loc[[grb]])
    # build a lightweight donor/query frame in canonical-score space for kernel_match_impute
    n_comp = donor_scores.shape[1]
    score_cols = [f'cc{k}' for k in range(n_comp)]
    donor_df = pd.DataFrame(donor_scores, index=donors.index, columns=score_cols)
    for c in TARGET_COLS:
        donor_df[c] = donors[c].values
    query_df = pd.DataFrame(query_scores, index=[grb], columns=score_cols)
    mu, sd = kernel_match_impute(query_df, donor_df, TARGET_COLS, score_cols, bandwidth=0.2)
    for c in TARGET_COLS:
        true_v = xray.loc[grb, c]
        if np.isfinite(true_v) and np.isfinite(mu.loc[grb, c]):
            loo_true[c].append(true_v)
            loo_pred[c].append(mu.loc[grb, c])

print(f'{"feature":<15}{"n":>5}{"corr":>8}{"RMSE":>10}')
for c in TARGET_COLS:
    t, p = np.array(loo_true[c]), np.array(loo_pred[c])
    if len(t) >= 3 and p.std() > 1e-8:
        corr = np.corrcoef(t, p)[0, 1]
        rmse = np.sqrt(np.mean((t - p) ** 2))
        print(f'{c:<15}{len(t):>5}{corr:>8.3f}{rmse:>10.4f}')

# ---- Final imputation for the real 83 optical-only GRBs ----
# (uses ALL usable anchors for the CCA fit -- correct here, since these
# query GRBs are genuinely unlabeled, not held-out anchors being scored)
donor_scores, query_scores, _ = fit_cca_and_project(overlap_u, xray, opt_only_proj)
n_comp = donor_scores.shape[1]
score_cols = [f'cc{k}' for k in range(n_comp)]
donor_df = pd.DataFrame(donor_scores, index=xray.index, columns=score_cols)
for c in TARGET_COLS:
    donor_df[c] = xray[c].values
query_df = pd.DataFrame(query_scores, index=opt_only_proj.index, columns=score_cols)

mu, sd = kernel_match_impute(query_df, donor_df, TARGET_COLS, score_cols, bandwidth=0.2)
imputed = {c: mu[c].values for c in TARGET_COLS}
imputed_std = {c: sd[c].values for c in TARGET_COLS}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp21_cca_ot.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'\nexp21_cca_ot: {len(out)} GRBs -> {out_path}')
