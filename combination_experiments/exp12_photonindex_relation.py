#!/usr/bin/env python3
"""exp12: physically-motivated substitution for PhotonIndex. In X-ray
spectroscopy, photon index and energy spectral index are related by the
identity Gamma_photon = beta_energy + 1 (N(E) ~ E^-Gamma vs F_nu ~ nu^-beta).
Beta is measured (calibrated) for every GRB, optical-only included -- so
PhotonIndex needs no borrowing/matching at all for this one feature. The
other 4 targets (Gamma, NH, Fluence, PeakFlux) still use the OT baseline
(exp08) since there's no equivalent shortcut for them. Validated against
the 78 anchors before use."""
import os
import numpy as np
from common import load_data, kernel_match_impute, assemble_and_save, TARGET_COLS, SCRIPT_DIR

xray, opt_clean, opt_proj, opt_only_proj, overlap = load_data()

# Validate the relation on the 78 anchors first: true PhotonIndex vs Beta_x + 1
anchor_true = xray.loc[overlap, 'PhotonIndex']
anchor_beta = xray.loc[overlap, 'Beta_x']
ok = anchor_true.notna() & anchor_beta.notna()
resid = anchor_true[ok] - (anchor_beta[ok] + 1)
corr = np.corrcoef(anchor_true[ok], anchor_beta[ok] + 1)[0, 1]
print(f'PhotonIndex ~ Beta+1 relation check on {ok.sum()} anchors: corr={corr:.3f}, '
      f'residual std={resid.std():.3f}, residual mean={resid.mean():.3f}')

MATCH_COLS = ['logFa_x', 'logTa_x', 'Alpha_x', 'Beta_x']
opt_anch = opt_proj.loc[overlap, MATCH_COLS].values.astype(float)
xray_anch = xray.loc[overlap, MATCH_COLS].values.astype(float)
usable = ~(np.isnan(opt_anch).any(axis=1) | np.isnan(xray_anch).any(axis=1))
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(np.vstack([opt_proj[MATCH_COLS].values, xray[MATCH_COLS].values]).astype(float))
diffs = scaler.transform(opt_anch[usable]) - scaler.transform(xray_anch[usable])
dim_var = np.maximum(diffs.var(axis=0), 1e-3)
weights = (1.0 / dim_var); weights = weights / weights.mean()

other_targets = ['Gamma', 'log10NH', 'log10Fluence', 'log10PeakFlux']
mu, sd = kernel_match_impute(opt_only_proj, xray, other_targets, MATCH_COLS, weights=weights, bandwidth=0.1)

imputed = {col: mu[col].values for col in other_targets}
imputed_std = {col: sd[col].values for col in other_targets}
# PhotonIndex via the physical relation, using calibrated Beta_x (opt_only_proj)
imputed['PhotonIndex'] = (opt_only_proj['Beta_x'] + 1 - resid.mean()).values  # bias-correct by the anchor offset
imputed_std['PhotonIndex'] = [resid.std()] * len(opt_only_proj)

out_path = os.path.join(SCRIPT_DIR, 'results', 'exp12_photonindex_relation.csv')
out = assemble_and_save(xray, opt_clean, opt_proj, opt_only_proj, imputed, imputed_std, out_path)
print(f'exp12_photonindex_relation: {len(out)} GRBs -> {out_path}')
