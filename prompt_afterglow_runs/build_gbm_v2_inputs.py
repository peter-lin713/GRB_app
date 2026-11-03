#!/usr/bin/env python3
"""Inputs for the second series (runs 10-17), all on the corrected time+position
GBM crosswalk (build_gbm_swap.crosswalk_time). Peak flux = GBM flux_1024
(10-1000 keV, 1.024 s peak photon flux).

gbmpf2     Dainotti's scheme: GRBs in both Swift and GBM get GBM's measured peak
           flux directly; Swift-only keep Swift's; GBM-only (no Swift peak flux)
           get the GBM value transferred to the Swift scale with Theil-Sen
           (Swift ~ GBM fit on the overlap). In these datasets every GRB has a
           Swift peak flux, so the transfer branch fills 0 rows -- reported below.
gbmpf2cal  Old run-5 method on the corrected matches: Theil-Sen-calibrated GBM
           value (bootstrap median, build_inputs.substitute_peakflux) on the overlap.
gbmpf2_epeak0 / _epeakmice
           gbmpf2 + log10 GBM E_peak (flnc_comp_epeak, Spencer's quality filter);
           missing -> 0 or NA for MICE.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import theilslopes
from build_inputs import tap, substitute_peakflux, DATA, HERE, MIN_ERR
from build_gbm_swap import crosswalk_time

LN10 = np.log(10)
BASES = {'theilsen': 'superlearner_training_theilsen_errcut_relative.csv',
         'theilsen_multivariate_reversed_xray_to_opt':
             'superlearner_training_theilsen_multivariate_reversed_xray_to_opt_errcut_relative.csv'}


def dainotti_peakflux(d, cw, label):
    v = cw.flux_1024.reindex(d.index)
    e = cw.flux_1024_error.reindex(d.index)
    gbm = (v > 0) & np.isfinite(v)
    sw = d.log10PeakFlux.notna()
    both, gonly = gbm & sw, gbm & ~sw
    out = d.copy()
    out.loc[both, 'log10PeakFlux'] = np.log10(v[both])
    out.loc[both, 'log10PeakFluxErr'] = (e[both] / (v[both] * LN10)).fillna(MIN_ERR).clip(lower=MIN_ERR)
    m, b, _, _ = theilslopes(d.log10PeakFlux[both], np.log10(v[both]))
    if gonly.any():
        out.loc[gonly, 'log10PeakFlux'] = m * np.log10(v[gonly]) + b
        res = d.log10PeakFlux[both] - (m * np.log10(v[both]) + b)
        out.loc[gonly, 'log10PeakFluxErr'] = 1.4826 * np.median(np.abs(res))
    out['gbm_peakflux_substituted'] = gbm.astype(int).values
    print(f'[{label}] {len(d)} GRBs: {int(both.sum())} GBM measured value swapped in, '
          f'{int(gonly.sum())} GBM-only transferred (Swift ~ {m:.3f} GBM {b:+.3f}), '
          f'{int((sw & ~gbm).sum())} keep Swift')
    return out


def main():
    bases = {k: pd.read_csv(os.path.join(DATA, f), index_col=0) for k, f in BASES.items()}
    keys = sorted(set().union(*(set(d.index) for d in bases.values())))
    cw = crosswalk_time(keys)
    ex = tap('fermigbrst_extra', 'SELECT name,t50,flux_64,flux_256,flnc_comp_epeak,flnc_comp_epeak_pos_err,'
                                 'flnc_comp_index,pflx_comp_epeak,flnc_best_fitting_model FROM fermigbrst')
    ex['name'] = ex.name.str.strip()
    ex = ex.drop_duplicates('name').set_index('name')
    cw['flnc_comp_epeak'] = ex.flnc_comp_epeak.reindex(cw.name).values
    cw['flnc_comp_epeak_pos_err'] = ex.flnc_comp_epeak_pos_err.reindex(cw.name).values
    cw.to_csv(os.path.join(HERE, 'gbm_crosswalk_time_all.csv'))
    print(f'Corrected crosswalk: {len(cw)} of {len(keys)} GRBs matched')

    for k, d in bases.items():
        dainotti_peakflux(d, cw, k).to_csv(os.path.join(DATA, f'superlearner_training_{k}_errcut_relative_gbmpf2.csv'))

    cal_cw = pd.DataFrame({'pflx_comp_phtflux': cw.flux_1024, 'pflx_comp_phtflux_error': cw.flux_1024_error})
    substitute_peakflux(bases['theilsen'], cal_cw, 'theilsen calibrated').to_csv(
        os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative_gbmpf2cal.csv'))

    base = pd.read_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative_gbmpf2.csv'), index_col=0)
    ep = cw.flnc_comp_epeak.reindex(base.index)
    epe = cw.flnc_comp_epeak_pos_err.reindex(base.index)
    good = (ep > 8) & (ep < 1e4) & (epe < ep)
    log_ep = np.where(good, np.log10(ep.where(good)), np.nan)
    for tag, fill in (('epeak0', 0.0), ('epeakmice', np.nan)):
        e = base.copy()
        e['log10Epeak'] = np.where(np.isfinite(log_ep), log_ep, fill)
        e.to_csv(os.path.join(DATA, f'superlearner_training_theilsen_errcut_relative_gbmpf2_{tag}.csv'))
    print(f'E_peak: {int(good.sum())} of {len(base)} GRBs have a well-constrained GBM E_peak')


if __name__ == '__main__':
    main()
