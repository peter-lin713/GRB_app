#!/usr/bin/env python3
"""Flagship dataset with GBM's own measurements swapped in (no Theil-Sen) for the
prompt columns Swift BAT and Fermi GBM both measure, on GRBs GBM also detected.

Crosswalk (replaces the date-based one in build_inputs.py, which accepted any GBM
burst on the same day when there was only one -- 36 of its 130 matches were
other bursts): a Swift GRB matches the GBM trigger closest in time with
  |t_GBM - t_BAT| < 1000 s   and   sep < 3 * sqrt(err_GBM^2 + 3.7^2) deg,  err_GBM < 10 deg
(3.7 deg = GBM systematic localization error). Swift times/positions: swiftgrbba.

Swapped (GBM column -> dataset column; errors converted to the dataset's form):
  t90        (50-300 keV)            -> log10T90,      T90Err (s)
  fluence    (10-1000 keV, erg/cm^2) -> log10Fluence,  log10FluenceErr
  flux_1024  (10-1000 keV, 1.024 s peak photon flux) -> log10PeakFlux, log10PeakFluxErr
PhotonIndex is swapped only in the _pi variant: GBM's single power-law index for
the fluence spectrum (flnc_plaw_index, sign flipped to Swift's convention).

Writes Data/superlearner_training_theilsen_errcut_relative_gbmswap.csv and
_gbmswap_pi.csv, plus gbm_crosswalk_time.csv here."""
import os
import numpy as np
import pandas as pd
from build_inputs import tap, norm, stem, sep, HERE, DATA, MIN_ERR

LN10 = np.log(10)


def crosswalk_time(keys):
    sb = tap('swiftgrbba_time', 'SELECT name,time,ra,dec FROM swiftgrbba')
    sb['key'] = sb.name.map(norm); sb['stem'] = sb.key.map(stem)
    bf = {k: i for i, k in enumerate(sb.key)}
    sc = sb.stem.value_counts()
    bs = {s: i for i, s in enumerate(sb.stem) if sc[s] == 1}
    g = tap('fermigbrst_swap', 'SELECT name,trigger_time,ra,dec,error_radius,t90,t90_error,'
                               'fluence,fluence_error,flux_1024,flux_1024_error,'
                               'flnc_plaw_index,flnc_plaw_index_pos_err,flnc_plaw_index_neg_err FROM fermigbrst')
    g['name'] = g.name.str.strip()
    rows = []
    for k in keys:
        i = bf.get(k, bs.get(stem(k)))
        if i is None:
            continue
        dt = (g.trigger_time.values - sb.time.iloc[i]) * 86400
        c = np.flatnonzero(np.abs(dt) < 1000)
        if len(c) == 0:
            continue
        j = c[np.argmin(np.abs(dt[c]))]
        s = float(sep(sb.ra.iloc[i], sb.dec.iloc[i], g.ra.values[j], g.dec.values[j]))
        err = g.error_radius.values[j]
        if not (err < 10 and s < 3 * np.hypot(err, 3.7)):
            continue
        rows.append(dict(key=k, dt_s=round(dt[j], 2), sep_deg=round(s, 3), **g.iloc[j].to_dict()))
    return pd.DataFrame(rows).set_index('key')


def main():
    d = pd.read_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative.csv'), index_col=0)
    cw = crosswalk_time(d.index)
    cw.to_csv(os.path.join(HERE, 'gbm_crosswalk_time.csv'))
    old = pd.read_csv(os.path.join(HERE, 'gbm_crosswalk_306.csv'), index_col=0).reindex(d.index).dropna(subset=['fermi_name'])
    agree = (old.fermi_name.str.strip() == cw.fermi_name.reindex(old.index) if 'fermi_name' in cw else
             old.fermi_name.str.strip() == cw.name.reindex(old.index))
    print(f'Time+position crosswalk: {len(cw)} of {len(d)} GRBs matched to a GBM trigger '
          f'(old date-based crosswalk: {len(old)}, of which {int(agree.sum())} agree)')

    out = d.copy()
    swaps = [('log10T90', 'T90Err', cw.t90, cw.t90_error, lambda v, e: e),
             ('log10Fluence', 'log10FluenceErr', cw.fluence, cw.fluence_error, lambda v, e: e / (v * LN10)),
             ('log10PeakFlux', 'log10PeakFluxErr', cw.flux_1024, cw.flux_1024_error, lambda v, e: e / (v * LN10))]
    flag = pd.Series(0, index=d.index)
    for col, ecol, v, e, conv in swaps:
        ok = (v > 0) & np.isfinite(v)
        idx = ok[ok].index
        diff = np.log10(v[idx]) - d.loc[idx, col]
        print(f'  {col:14s} swapped {len(idx):3d} | GBM - Swift (log10): median {diff.median():+.3f}, '
              f'MAD-sigma {1.4826 * (diff - diff.median()).abs().median():.3f}, r {np.corrcoef(np.log10(v[idx]), d.loc[idx, col])[0, 1]:.3f}')
        out.loc[idx, col] = np.log10(v[idx])
        out.loc[idx, ecol] = conv(v[idx], e[idx]).fillna(MIN_ERR).clip(lower=MIN_ERR)
        flag[idx] = 1
    out['gbm_swapped'] = flag.values
    out.to_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative_gbmswap.csv'))

    pi = -cw.flnc_plaw_index
    ok = pi.notna() & np.isfinite(pi)
    idx = ok[ok].index
    diff = pi[idx] - d.loc[idx, 'PhotonIndex']
    print(f'  {"PhotonIndex":14s} swapped {len(idx):3d} | GBM - Swift: median {diff.median():+.3f}, '
          f'MAD-sigma {1.4826 * (diff - diff.median()).abs().median():.3f}, r {np.corrcoef(pi[idx], d.loc[idx, "PhotonIndex"])[0, 1]:.3f}  (_pi file only)')
    out_pi = out.copy()
    out_pi.loc[idx, 'PhotonIndex'] = pi[idx]
    out_pi.loc[idx, 'PhotonIndexErr'] = ((cw.flnc_plaw_index_pos_err.abs() + cw.flnc_plaw_index_neg_err.abs()) / 2)[idx].fillna(MIN_ERR).clip(lower=MIN_ERR)
    out_pi.to_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative_gbmswap_pi.csv'))


if __name__ == '__main__':
    main()
