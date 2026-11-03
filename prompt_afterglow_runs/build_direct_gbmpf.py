#!/usr/bin/env python3
"""Flagship dataset with the GBM 1-s peak photon flux swapped in directly
(no Theil-Sen calibration), for GRBs that have one. Everyone else keeps
Swift's log10PeakFlux. Uses the crosswalk saved by build_inputs.py.

  log10PeakFlux    <- log10(pflx_comp_phtflux)
  log10PeakFluxErr <- pflx_comp_phtflux_error / (pflx_comp_phtflux * ln 10)

Writes Data/superlearner_training_theilsen_errcut_relative_gbmpf_direct.csv."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), 'Data')
MIN_ERR = 1e-6

d = pd.read_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative.csv'), index_col=0)
cw = pd.read_csv(os.path.join(HERE, 'gbm_crosswalk_306.csv'), index_col=0)
pf = cw.pflx_comp_phtflux.reindex(d.index)
pfe = cw.pflx_comp_phtflux_error.reindex(d.index)
ok = (pf > 0) & np.isfinite(pf)

out = d.copy()
out.loc[ok, 'log10PeakFlux'] = np.log10(pf[ok])
out.loc[ok, 'log10PeakFluxErr'] = (pfe[ok] / (pf[ok] * np.log(10))).fillna(MIN_ERR).clip(lower=MIN_ERR)
out['gbm_peakflux_substituted'] = ok.astype(int).values
out.to_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative_gbmpf_direct.csv'))

diff = out.loc[ok, 'log10PeakFlux'] - d.loc[ok, 'log10PeakFlux']
print(f'{len(out)} GRBs, {int(ok.sum())} with GBM peak flux substituted directly; '
      f'GBM - Swift (log10): median {diff.median():+.3f}, MAD-sigma {1.4826 * (diff - diff.median()).abs().median():.3f}')
