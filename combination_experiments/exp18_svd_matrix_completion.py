#!/usr/bin/env python3
"""exp18: iterative low-rank matrix completion (same family as the
Netflix-prize-style recommender algorithms). Stack the 4 shared params + 5
targets for all 306 GRBs (223 X-ray complete, 83 optical-only missing the 5
targets) into one matrix, initialize missing entries at the column mean,
then repeatedly: standardize, take a rank-k SVD, reconstruct, re-fill ONLY
the originally-missing cells, and repeat to convergence. Treats the whole
feature table as a single structure to complete, rather than predicting
target columns from feature columns."""
import os
import numpy as np
from common import load_data, assemble_and_save, TARGET_COLS, FEATURE_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
ALL_COLS = FEATURE_COLS + TARGET_COLS

xray_mat = xray[ALL_COLS].values.astype(float)
opt_mat = opt_only_proj[FEATURE_COLS].values.astype(float)
opt_mat = np.hstack([opt_mat, np.full((len(opt_mat), len(TARGET_COLS)), np.nan)])
full = np.vstack([xray_mat, opt_mat])
missing = np.isnan(full)

col_mean = np.nanmean(full, axis=0)
col_std = np.nanstd(full, axis=0)
col_std[col_std == 0] = 1
filled = np.where(missing, np.tile(col_mean, (full.shape[0], 1)), full)

RANK = 4
for it in range(50):
    z = (filled - col_mean) / col_std
    U, S, Vt = np.linalg.svd(z, full_matrices=False)
    S_trunc = np.zeros_like(S); S_trunc[:RANK] = S[:RANK]
    recon = (U * S_trunc) @ Vt
    recon = recon * col_std + col_mean
    new_filled = np.where(missing, recon, full)
    delta = np.abs(new_filled[missing] - filled[missing]).max() if missing.any() else 0
    filled = new_filled
    if delta < 1e-5:
        break
print(f'SVD matrix completion converged after {it+1} iterations')

n_target = len(TARGET_COLS)
opt_filled = filled[len(xray_mat):, -n_target:]
imputed = {col: opt_filled[:, k] for k, col in enumerate(TARGET_COLS)}
imputed_std = {col: [xray[col].std()] * len(opt_only_proj) for col in TARGET_COLS}

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp18_svd_matrix_completion.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp18_svd_matrix_completion: {len(out)} GRBs -> {out_path}')
