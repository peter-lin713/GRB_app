#!/usr/bin/env python3
"""exp14: hard cluster pooling. Cluster the 223 X-ray GRBs into k groups by
their calibrated shared parameters, assign each optical-only GRB to its
nearest cluster centroid, impute each target feature as a shrinkage blend
of that cluster's mean and the global mean (partial pooling) -- softer and
more stable than single-nearest-neighbor OT matching, at the cost of
losing fine-grained per-GRB resolution."""
import os
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']
K = 6
SHRINK = 0.5  # blend weight toward cluster mean vs global mean

imp = SimpleImputer(strategy='median').fit(xray[MATCH_COLS].values)
X_train = imp.transform(xray[MATCH_COLS].values)
X_query = imp.transform(opt_only_proj[MATCH_COLS].values)
scaler = StandardScaler().fit(X_train)
X_train_z = scaler.transform(X_train)
X_query_z = scaler.transform(X_query)

km = KMeans(n_clusters=K, n_init=10, random_state=42).fit(X_train_z)
xray_cluster = km.labels_
query_cluster = km.predict(X_query_z)

imputed, imputed_std = {}, {}
for col in TARGET_COLS:
    y = xray[col].values
    global_mean, global_std = np.nanmean(y), np.nanstd(y)
    cluster_means = {c: np.nanmean(y[xray_cluster == c]) if np.isfinite(np.nanmean(y[xray_cluster == c])) else global_mean
                      for c in range(K)}
    cluster_stds = {c: np.nanstd(y[xray_cluster == c]) if (xray_cluster == c).sum() > 3 else global_std
                     for c in range(K)}
    pred = np.array([SHRINK * cluster_means[c] + (1 - SHRINK) * global_mean for c in query_cluster])
    imputed[col] = pred
    imputed_std[col] = np.array([cluster_stds[c] for c in query_cluster])

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp14_kmeans_cluster.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp14_kmeans_cluster: {len(out)} GRBs -> {out_path}')
