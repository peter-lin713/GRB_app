#!/usr/bin/env python3
"""
Build the GBM-augmented input datasets for the prompt/afterglow SuperLearner runs.

1. Crosswalk every GRB in the two base datasets to the Fermi GBM catalog
   (fermigbrst), using Spencer's exact logic from
   GRB_linear_model_comparison/next_steps_bat_gbm/gbm_crosswalk.py
   (spencer_most_recent): date-primary match; if several GBM bursts share the
   day, nearest sky position (Swift position from swiftgrbba) within 8 deg.
   Re-queried here to also pull the GBM Fermi name, the 1-s peak photon flux
   error and E_peak, which Spencer's crosswalk did not save.
2. Peak flux: Theil-Sen fit of Swift log10PeakFlux ~ log10(GBM 1-s peak photon
   flux, pflx_comp_phtflux) on the GRBs with both; GRBs with a GBM value get the
   calibrated GBM value (bootstrap value/error propagation, as in
   combine_optical_xray_theilsen_reversed.py); everyone else keeps Swift's.
   Fit separately per base dataset against its own Swift column. No z used.
3. E_peak (flnc_comp_epeak, Spencer's quality filter 8 < Ep < 1e4 keV and
   pos_err < Ep), as log10Epeak: missing -> 0 (one file) or NA for MICE (other).

HEASARC is queried with curl (python here lacks SSL certs); raw VOTables are
saved under heasarc_raw/ and parsed with the standard-library XML parser.

Outputs (Data/):
  superlearner_training_theilsen_errcut_relative_gbmpf.csv
  superlearner_training_theilsen_multivariate_reversed_xray_to_opt_errcut_relative_gbmpf.csv
  superlearner_training_theilsen_errcut_relative_gbmpf_epeak0.csv
  superlearner_training_theilsen_errcut_relative_gbmpf_epeakmice.csv
(The flagship 0.707 formula search used the _errcut_relative file -- byte-identical
to GAM_supercomputer/theilsen_cleaned_formula_generation/superlearner_training_theilsen.csv.)
plus gbm_crosswalk_306.csv (this folder) for provenance.
"""
import base64
import os
import struct
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.stats import theilslopes

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(REPO, 'Data')
RAW = os.path.join(HERE, 'heasarc_raw')
BASES = {
    'theilsen': 'superlearner_training_theilsen_errcut_relative.csv',
    'theilsen_multivariate_reversed_xray_to_opt':
        'superlearner_training_theilsen_multivariate_reversed_xray_to_opt_errcut_relative.csv',
}
N_BOOT, MIN_ERR = 1000, 1e-6


def tap(name, query):
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, name + '.xml')
    if not os.path.exists(path):
        subprocess.run(['curl', '-s', '--max-time', '300', '-G',
                        'https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync',
                        '--data-urlencode', 'REQUEST=doQuery', '--data-urlencode', 'LANG=ADQL',
                        '--data-urlencode', 'FORMAT=votable', '--data-urlencode', 'QUERY=' + query,
                        '-o', path], check=True)
    ns = {'v': 'http://www.ivoa.net/xml/VOTable/v1.3'}
    root = ET.parse(path).getroot()
    status = root.find('.//v:INFO[@name="QUERY_STATUS"]', ns)
    if status is not None and status.get('value') != 'OK':
        raise RuntimeError(f'{name}: {status.text}')
    fields = [(f.get('name'), f.get('datatype'), f.get('arraysize'))
              for f in root.iter('{%s}FIELD' % ns['v'])]
    stream = root.find('.//v:BINARY/v:STREAM', ns)
    return _decode_votable_binary(base64.b64decode(stream.text), fields)


def _decode_votable_binary(buf, fields):
    """VOTable 1.x BINARY serialization: rows back to back; variable-length
    char ('*') = 4-byte big-endian length + bytes; numbers big-endian, NaN = null."""
    fmt = {'double': ('>d', 8), 'float': ('>f', 4), 'int': ('>i', 4),
           'long': ('>q', 8), 'short': ('>h', 2), 'unsignedByte': ('>B', 1)}
    rows, p = [], 0
    while p < len(buf):
        row = []
        for name, dt, size in fields:
            if dt == 'char':
                if size and size.endswith('*'):
                    n = struct.unpack_from('>i', buf, p)[0]; p += 4
                else:
                    n = int(size or 1)
                row.append(buf[p:p + n].decode('utf-8', 'replace').rstrip('\x00').strip()); p += n
            else:
                f, n = fmt[dt]
                row.append(struct.unpack_from(f, buf, p)[0]); p += n
        rows.append(row)
    return pd.DataFrame(rows, columns=[f[0] for f in fields])


norm = lambda s: str(s).replace('GRB', '').replace(' ', '').strip().upper()
stem = lambda k: str(k)[:6] if len(str(k)) >= 6 and str(k)[:6].isdigit() else str(k)


def sep(r1, d1, r2, d2):
    a, b, c, e = map(np.radians, [r1, d1, r2, d2])
    return np.degrees(np.arccos(np.clip(np.sin(b) * np.sin(e) + np.cos(b) * np.cos(e) * np.cos(a - c), -1, 1)))


def crosswalk(keys):
    sb = tap('swiftgrbba', 'SELECT name,ra,dec FROM swiftgrbba')
    sb['key'] = sb.name.map(norm); sb['stem'] = sb.key.map(stem)
    bf = {k: i for i, k in enumerate(sb.key)}
    sc = sb.stem.value_counts()
    bs = {s: i for i, s in enumerate(sb.stem) if sc.get(s, 0) == 1}
    g = tap('fermigbrst', 'SELECT name,ra,dec,flnc_comp_phtflux,pflx_comp_phtflux,'
                          'pflx_comp_phtflux_error,flnc_comp_epeak,flnc_comp_epeak_pos_err,'
                          'flnc_comp_epeak_neg_err FROM fermigbrst')
    g['key'] = g.name.map(norm); g['date'] = g.key.map(stem); g = g[g.ra.notna()]
    rows = []
    for k in keys:
        i = bf.get(k, bs.get(stem(k)))
        ra, dec = (sb.ra.iloc[i], sb.dec.iloc[i]) if i is not None else (np.nan, np.nan)
        cand = g[g.date == stem(k)]
        if len(cand) == 0:
            continue
        if len(cand) == 1:
            c, sp = cand.iloc[0], np.nan
        else:
            if not np.isfinite(ra):
                continue
            s = sep(ra, dec, cand.ra.values, cand.dec.values); j = int(np.argmin(s))
            if s[j] >= 8.0:
                continue
            c, sp = cand.iloc[j], round(float(s[j]), 2)
        rows.append(dict(key=k, fermi_name=c['name'], sep_deg=sp, n_sameday=len(cand),
                         flnc_comp_phtflux=c.flnc_comp_phtflux,
                         pflx_comp_phtflux=c.pflx_comp_phtflux,
                         pflx_comp_phtflux_error=c.pflx_comp_phtflux_error,
                         flnc_comp_epeak=c.flnc_comp_epeak,
                         flnc_comp_epeak_pos_err=c.flnc_comp_epeak_pos_err))
    return pd.DataFrame(rows).set_index('key')


def substitute_peakflux(d, cw, label):
    pf = cw.pflx_comp_phtflux.reindex(d.index)
    pfe = cw.pflx_comp_phtflux_error.reindex(d.index)
    ok = (pf > 0) & np.isfinite(pf)
    x = np.log10(pf[ok])
    xe = (pfe[ok] / (pf[ok] * np.log(10))).fillna(MIN_ERR).clip(lower=MIN_ERR)
    both = ok & d.log10PeakFlux.notna()
    xs, ys = np.log10(pf[both]).values, d.loc[both, 'log10PeakFlux'].values
    m, b, lo, hi = theilslopes(ys, xs)
    sigma_int = np.median(np.abs(ys - (m * xs + b))) * 1.4826
    rng = np.random.default_rng(2000)
    boot = np.array([theilslopes(ys[ix], xs[ix])[:2]
                     for ix in (rng.integers(0, len(xs), len(xs)) for _ in range(N_BOOT))])
    rng = np.random.default_rng(42)
    xp = x.values[None, :] + rng.standard_normal((N_BOOT, len(x))) * xe.values[None, :]
    pred = boot[:, :1] * xp + boot[:, 1:] + rng.standard_normal((N_BOOT, len(x))) * sigma_int
    out = d.copy()
    out.loc[x.index, 'log10PeakFlux'] = np.median(pred, 0)
    out.loc[x.index, 'log10PeakFluxErr'] = (np.percentile(pred, 84, 0) - np.percentile(pred, 16, 0)) / 2
    out['gbm_peakflux_substituted'] = out.index.isin(x.index).astype(int)
    r = np.corrcoef(xs, ys)[0, 1]
    print(f'[{label}] Theil-Sen Swift~GBM peak flux: m={m:.3f} ({lo:.3f}-{hi:.3f}) b={b:.3f} '
          f'sigma_int={sigma_int:.3f} r={r:.3f} N_overlap={len(xs)} | substituted {ok.sum()} GRBs')
    return out


def main():
    bases = {k: pd.read_csv(os.path.join(DATA, f), index_col=0) for k, f in BASES.items()}
    keys = sorted(set().union(*(set(d.index) for d in bases.values())))
    cw = crosswalk(keys)
    cw.to_csv(os.path.join(HERE, 'gbm_crosswalk_306.csv'))
    print(f'Crosswalk: {len(cw)} of {len(keys)} GRBs matched to fermigbrst; '
          f'{int((cw.pflx_comp_phtflux > 0).sum())} with a 1-s peak flux')

    # Sanity: reproduce Spencer's crosswalk on his keys (same GBM burst => same flnc flux)
    sp = pd.read_csv(os.path.join(HERE, 'inputs', 'gbm_crosswalk_281.csv')).set_index('key')
    shared = sp.index.intersection(cw.index)
    same = np.isclose(sp.loc[shared, 'flnc_comp_phtflux'], cw.loc[shared, 'flnc_comp_phtflux'],
                      rtol=1e-6, equal_nan=True)
    print(f'Spencer crosswalk check: {len(sp)} of his matches, {len(shared)} also matched here, '
          f'{int(same.sum())} identical GBM burst')

    for k, d in bases.items():
        out = substitute_peakflux(d, cw, k)
        out.to_csv(os.path.join(DATA, f'superlearner_training_{k}_errcut_relative_gbmpf.csv'))

    # E_peak variants (flagship + GBM peak flux)
    ep = cw.flnc_comp_epeak.reindex(bases['theilsen'].index)
    epe = cw.flnc_comp_epeak_pos_err.reindex(bases['theilsen'].index)
    good = (ep > 8) & (ep < 1e4) & (epe < ep)
    log_ep = np.where(good, np.log10(ep), np.nan)
    base = pd.read_csv(os.path.join(DATA, 'superlearner_training_theilsen_errcut_relative_gbmpf.csv'), index_col=0)
    for tag, fill in (('epeak0', 0.0), ('epeakmice', np.nan)):
        e = base.copy()
        e['log10Epeak'] = np.where(np.isfinite(log_ep), log_ep, fill)
        e.to_csv(os.path.join(DATA, f'superlearner_training_theilsen_errcut_relative_gbmpf_{tag}.csv'))
    print(f'E_peak: {int(good.sum())} of {len(base)} GRBs have a well-constrained GBM E_peak')


if __name__ == '__main__':
    main()
