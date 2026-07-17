#!/usr/bin/env python3
"""exp16: Gaussian copula conditional imputation. Rank-transform all 9
numeric columns (4 shared params + 5 targets) to a standard-normal scale via
each column's own empirical CDF on the 223 X-ray GRBs, fit the full 9x9
correlation matrix in that normal space (captures the whole joint
dependency structure at once, not just pairwise distance), then use the
standard multivariate-normal conditional-distribution formula to predict
the 5 (still-normal-scale) targets given the 4 known shared params, and
invert the transform back to the original scale. More principled than
matching on distance: models the actual joint distribution shape."""
import os
import numpy as np
from scipy.stats import norm
from common import load_data, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()
MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']
ALL_COLS = MATCH_COLS + TARGET_COLS

xray_mat = xray[ALL_COLS].values.astype(float)


def to_normal_scale(col_vals, ref_vals):
    """Rank-transform col_vals against ref_vals' empirical CDF, to N(0,1)."""
    ref_ok = ref_vals[~np.isnan(ref_vals)]
    ranks = np.searchsorted(np.sort(ref_ok), col_vals, side='right') / (len(ref_ok) + 1)
    ranks = np.clip(ranks, 1 / (len(ref_ok) + 2), 1 - 1 / (len(ref_ok) + 2))
    return norm.ppf(ranks)


def from_normal_scale(z_vals, ref_vals):
    ref_ok = np.sort(ref_vals[~np.isnan(ref_vals)])
    p = norm.cdf(z_vals)
    idx = (p * (len(ref_ok) - 1)).astype(int)
    idx = np.clip(idx, 0, len(ref_ok) - 1)
    return ref_ok[idx]


# Build normal-scale matrix on complete X-ray rows only (for the correlation fit)
Z = np.full_like(xray_mat, np.nan)
for j, col in enumerate(ALL_COLS):
    Z[:, j] = to_normal_scale(xray_mat[:, j], xray_mat[:, j])

complete = ~np.isnan(Z).any(axis=1)
corr = np.corrcoef(Z[complete].T)  # 9x9
n1 = len(MATCH_COLS)
S11 = corr[:n1, :n1]
S21 = corr[n1:, :n1]
S11_inv = np.linalg.pinv(S11)

# Query: optical-only GRBs' 4 known shared params, normal-scaled against the SAME X-ray reference
opt_match_raw = opt_only_proj[MATCH_COLS].values.astype(float)
Zq = np.column_stack([to_normal_scale(opt_match_raw[:, j], xray_mat[:, j]) for j in range(n1)])

cond_mean_z = Zq @ (S21 @ S11_inv).T   # (83, 5) in normal space
S22 = corr[n1:, n1:]
cond_cov_z = S22 - S21 @ S11_inv @ S21.T
cond_std_z = np.sqrt(np.diag(cond_cov_z))

imputed, imputed_std = {}, {}
for k, col in enumerate(TARGET_COLS):
    pred = from_normal_scale(cond_mean_z[:, k], xray_mat[:, n1 + k])
    # propagate the normal-space conditional std through the local slope of the inverse transform
    hi = from_normal_scale(cond_mean_z[:, k] + cond_std_z[k], xray_mat[:, n1 + k])
    lo = from_normal_scale(cond_mean_z[:, k] - cond_std_z[k], xray_mat[:, n1 + k])
    imputed[col] = pred
    imputed_std[col] = np.abs(hi - lo) / 2.0

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp16_copula.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp16_copula: {len(out)} GRBs -> {out_path}')
