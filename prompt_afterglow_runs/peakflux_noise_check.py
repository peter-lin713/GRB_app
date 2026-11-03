#!/usr/bin/env python3
"""
Does substituting the (Theil-Sen-calibrated) GBM peak flux add noise, or is the
problem mixing two instruments in one column?  Honest linear check on the
flagship prepared frame (281 GRBs; Spencer's protocol: 100 paired 80/20 splits,
OLS, r in linear z).
Variants of the peak-flux information, all else = Spencer's SIMPLE7 features:
  swift       Swift BAT peak flux for everyone (baseline)
  mixed_cal   GBM value, Theil-Sen-calibrated to the Swift scale, where GBM saw the
              burst; Swift elsewhere (what the full-pipeline GBM runs used)
  mixed_raw   same but uncalibrated raw log GBM peak flux (offset not removed)
  swift+gbm   Swift peak flux kept + calibrated GBM as an extra column
              (train-mean fill + missing indicator; fit on the training split)
Also on the GBM-overlap subset only: swift vs calibrated GBM, same GRBs.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
d = pd.read_csv(os.path.join(REPO, 'GAM_supercomputer', 'pa_04_all_swiftpf_formula_generation',
                             'grb_xray_imputed.csv'), index_col=0)
cal = pd.read_csv(os.path.join(REPO, 'Data', 'superlearner_training_theilsen_errcut_relative_gbmpf.csv'),
                  index_col=0).log10PeakFlux.reindex(d.index)
cw = pd.read_csv(os.path.join(HERE, 'gbm_crosswalk_306.csv'), index_col=0)
raw = np.log10(cw.pflx_comp_phtflux.where(cw.pflx_comp_phtflux > 0)).reindex(d.index)
has = raw.notna()
cal_gbm = cal.where(has)                       # calibrated GBM value, only where GBM measured
print(f'GBM peak flux on {int(has.sum())} of {len(d)} GRBs | Swift vs calibrated GBM on those: '
      f'r = {np.corrcoef(d.log10PeakFlux[has], cal_gbm[has])[0,1]:.3f}, '
      f'mean offset = {(cal_gbm[has] - d.log10PeakFlux[has]).mean():+.3f} dex, '
      f'scatter = {(cal_gbm[has] - d.log10PeakFlux[has]).std():.3f} dex')
for name, s in (('Swift', d.log10PeakFlux), ('GBM calibrated', cal_gbm)):
    print(f'  corr(log peak flux, log(1+z)) on the {int(has.sum())} overlap GRBs [{name:14s}] = '
          f'{np.corrcoef(s[has], d.log10z[has])[0,1]:+.3f}')

BASE = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'log10NH', 'PhotonIndex']


def honest(frame, pf, extra=None, n_splits=100, seed=1):
    y = frame.log10z.values
    X0 = np.column_stack([frame[BASE].values, pf.values])
    rng, n, rs = np.random.default_rng(seed), len(y), []
    for _ in range(n_splits):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.2 * n), replace=False)] = True
        X = X0
        if extra is not None:
            e = extra.values.copy(); miss = np.isnan(e); e[miss] = np.nanmean(e[~te & ~miss])
            X = np.column_stack([X0, e, miss.astype(float)])
        D = np.column_stack([np.ones(n), X])
        b = np.linalg.lstsq(D[~te], y[~te], rcond=None)[0]
        rs.append(np.corrcoef(10 ** (D[te] @ b) - 1, 10 ** y[te] - 1)[0, 1])
    return np.array(rs)


rows = []
base = honest(d, d.log10PeakFlux)
for label, pf, extra in (('mixed_cal (pipeline GBM runs)', d.log10PeakFlux.where(~has, cal_gbm), None),
                         ('mixed_raw (no calibration)',    d.log10PeakFlux.where(~has, raw), None),
                         ('swift + GBM as extra column',   d.log10PeakFlux, cal_gbm)):
    r = honest(d, pf, extra); dd = r - base
    rows.append(dict(sample=f'all {len(d)}', variant=label, swift_r=base.mean(), variant_r=r.mean(),
                     diff=dd.mean(), diff_sd=dd.std(), better_pct=100 * (dd > 0).mean()))
sub = d[has]
b2 = honest(sub, sub.log10PeakFlux); r2 = honest(sub, cal_gbm[has]); dd = r2 - b2
rows.append(dict(sample=f'GBM overlap {int(has.sum())}', variant='GBM calibrated instead of Swift',
                 swift_r=b2.mean(), variant_r=r2.mean(), diff=dd.mean(), diff_sd=dd.std(),
                 better_pct=100 * (dd > 0).mean()))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, 'peakflux_noise_check.csv'), index=False, float_format='%.4f')
pd.set_option('display.width', 200)
print('\n' + out.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
