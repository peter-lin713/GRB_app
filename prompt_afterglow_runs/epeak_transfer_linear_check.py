#!/usr/bin/env python3
"""
Quick honest-linear check: does Epeak help once BAT Epeak is transferred onto
the GBM scale with Theil-Sen?

Frame: flagship prepared frame (GAM_supercomputer/pa_04_all_swiftpf_formula_generation/
grb_xray_imputed.csv, 281 GRBs, Swift peak flux, MICE-imputed features).
Epeak sources (log10 keV):
  GBM  -- flnc_comp_epeak via our crosswalk (Spencer's quality filter), 109 GRBs
  BAT  -- direct bat_epeak from Spencer's swiftgrb pull (bat_match_281_v2.csv), 57 GRBs
  transfer: Theil-Sen  log GBM ~ log BAT  on the 25 GRBs with both (no z used);
            BAT-only GRBs get the calibrated value, GBM kept where it exists -> 141 GRBs
  MICE -- the Epeak column MICE imputed in run 9's prep (all 281; no z used)
Protocol: Spencer's honest held-out r -- 100 random 80/20 splits, OLS refit on
the 80%, Pearson r in linear z on the 20%; all variants on identical splits.
Missing Epeak in the full-sample variants: train-mean fill + missing indicator
(both computed on the training split only).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import theilslopes

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GAM = os.path.join(REPO, 'GAM_supercomputer')
SIMPLE7 = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'log10NH', 'PhotonIndex', 'log10PeakFlux']
ALL10 = SIMPLE7 + ['Beta', 'Gamma', 'log10Fluence']
N_SPLITS, SEED = 100, 1

d = pd.read_csv(os.path.join(GAM, 'pa_04_all_swiftpf_formula_generation', 'grb_xray_imputed.csv'), index_col=0)
gbm = pd.read_csv(os.path.join(REPO, 'Data', 'superlearner_training_theilsen_errcut_relative_gbmpf_epeakmice.csv'),
                  index_col=0).log10Epeak.reindex(d.index)
bat = pd.read_csv(os.path.join(HERE, 'inputs', 'bat_match_281_v2.csv')).set_index('key').bat_epeak
bat = np.log10(bat.where(bat > 0)).reindex(d.index)
mice = pd.read_csv(os.path.join(GAM, 'pa_09_all_gbmpf_epeakmice_formula_generation', 'grb_xray_imputed.csv'),
                   index_col=0).log10Epeak.reindex(d.index)

both = gbm.notna() & bat.notna()
m, b, lo, hi = theilslopes(gbm[both], bat[both])
r_cal = np.corrcoef(bat[both], gbm[both])[0, 1]
transfer = gbm.copy()
only_bat = bat.notna() & gbm.isna()
transfer[only_bat] = m * bat[only_bat] + b
print(f'Theil-Sen transfer  logEp_GBM = {m:.3f}*logEp_BAT + {b:.3f}  (slope CI {lo:.2f}-{hi:.2f}, '
      f'r={r_cal:.2f}, N={int(both.sum())})')
print(f'Epeak coverage of 281: GBM {int(gbm.notna().sum())} | BAT {int(bat.notna().sum())} | both {int(both.sum())} '
      f'| GBM+transferred BAT {int(transfer.notna().sum())}')
for name, s in (('GBM', gbm), ('GBM+BAT transfer', transfer), ('MICE', mice)):
    ok = s.notna()
    print(f'  corr(log Epeak, log(1+z)) [{name:16s}] r = {np.corrcoef(s[ok], d.log10z[ok])[0,1]:+.3f}  (N={int(ok.sum())})')


def honest(frame, feats, ep=None, seed=SEED):
    """Mean held-out r (linear z) over random 80/20 splits; optional Epeak
    column with train-mean fill + missing indicator."""
    y = frame.log10z.values
    X0 = frame[feats].values
    rng, n, rs = np.random.default_rng(seed), len(y), []
    for _ in range(N_SPLITS):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.2 * n), replace=False)] = True
        X = X0
        if ep is not None:
            e = ep.values.copy(); miss = np.isnan(e)
            e[miss] = np.nanmean(e[~te & ~miss])
            X = np.column_stack([X0, e, miss.astype(float)]) if miss.any() else np.column_stack([X0, e])
        D = np.column_stack([np.ones(n), X])
        beta = np.linalg.lstsq(D[~te], y[~te], rcond=None)[0]
        rs.append(np.corrcoef(10 ** (D[te] @ beta) - 1, 10 ** y[te] - 1)[0, 1])
    return np.array(rs)


rows = []
def report(label, base, new):
    dd = new - base
    rows.append(dict(test=label, base_r=base.mean(), with_epeak_r=new.mean(), diff=dd.mean(),
                     diff_sd=dd.std(), epeak_better_pct=100 * (dd > 0).mean()))

for feats, fname in ((SIMPLE7, 'SIMPLE7'), (ALL10, 'ALL10')):
    base = honest(d, feats)
    for name, s in (('GBM only (109)', gbm), ('GBM+BAT transfer (141)', transfer), ('MICE (281)', mice)):
        report(f'{fname} | all 281 | + {name}', base, honest(d, feats, s))
    # cleanest test: only GRBs with a measured / transferred Epeak, no fill at all
    for name, s in (('GBM only', gbm), ('GBM+BAT transfer', transfer)):
        sub = d[s.notna()]
        report(f'{fname} | {name} subset (N={len(sub)}) | + Epeak',
               honest(sub, feats), honest(sub, feats, s[s.notna()]))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, 'epeak_transfer_linear_check.csv'), index=False, float_format='%.4f')
pd.set_option('display.width', 200)
print('\n' + out.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
