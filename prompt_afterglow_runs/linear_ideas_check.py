#!/usr/bin/env python3
"""
Honest linear screen of ideas for using Fermi GBM data (corrected time+position
crosswalk, gbm_crosswalk_time.csv) plus a few feature ideas.

Frame: flagship prepared frame (pa_04 grb_xray_imputed.csv, 281 GRBs, MICE-imputed
Swift features). GBM values are put in after imputation, only where GBM measured.
Protocol: Spencer's -- 100 random 80/20 splits, OLS refit on the 80%, Pearson r in
linear z on the 20%. Every variant uses the same splits as its baseline, so the
per-split difference is paired. Extra columns with gaps: train-split-mean fill +
missing indicator. Offsets (GBM - Swift median on the overlap) use no redshift.
"""
import os
import numpy as np
import pandas as pd
from build_inputs import tap

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SIMPLE7 = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'log10NH', 'PhotonIndex', 'log10PeakFlux']
ALL10 = SIMPLE7 + ['Beta', 'Gamma', 'log10Fluence']
N_SPLITS, SEED = 100, 1

d = pd.read_csv(os.path.join(REPO, 'GAM_supercomputer', 'pa_04_all_swiftpf_formula_generation',
                             'grb_xray_imputed.csv'), index_col=0)
cw = pd.read_csv(os.path.join(HERE, 'gbm_crosswalk_time.csv'), index_col=0)
g = tap('fermigbrst_extra', 'SELECT name,t50,flux_64,flux_256,flnc_comp_epeak,flnc_comp_epeak_pos_err,'
                            'flnc_comp_index,pflx_comp_epeak,flnc_best_fitting_model FROM fermigbrst')
g['name'] = g.name.str.strip()
g = g.drop_duplicates('name').set_index('name').reindex(cw.name.values)
g.index = cw.index
G = pd.concat([cw, g], axis=1).reindex(d.index)
has = G.t90.notna()
print(f'{int(has.sum())} of {len(d)} GRBs in the frame have a GBM match')

lg = lambda s: np.log10(s.where(s > 0))
gbm = pd.DataFrame({'log10T90': lg(G.t90), 'log10Fluence': lg(G.fluence),
                    'log10PeakFlux': lg(G.flux_1024), 'PhotonIndex': -G.flnc_plaw_index}, index=d.index)
off = {c: (gbm[c] - d[c])[has].median() for c in gbm}
ep_ok = (G.flnc_comp_epeak > 8) & (G.flnc_comp_epeak < 1e4) & (G.flnc_comp_epeak_pos_err < G.flnc_comp_epeak)
log_ep = lg(G.flnc_comp_epeak.where(ep_ok))


def swap(cols, mode='raw'):
    f = d.copy()
    for c in cols:
        v = gbm[c] - (off[c] if mode == 'offset' else 0)
        if mode == 'avg':
            v = (d[c] + gbm[c] - off[c]) / 2
        f[c] = f[c].where(~v.notna(), v)
    return f


def honest(frame, feats, extra=None, rows=None):
    fr = frame if rows is None else frame[rows]
    y = fr.log10z.values
    X0 = fr[feats].values
    ex = [] if extra is None else [e.reindex(fr.index).values for e in extra]
    rng, n, rs = np.random.default_rng(SEED), len(y), []
    for _ in range(N_SPLITS):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.2 * n), replace=False)] = True
        cols = [X0]
        for e in ex:
            e = e.astype(float).copy(); miss = np.isnan(e)
            if miss.all():
                continue
            e[miss] = np.nanmean(e[~te & ~miss])
            cols += [e[:, None]] + ([miss[:, None].astype(float)] if miss.any() else [])
        D = np.column_stack([np.ones(n)] + cols)
        b = np.linalg.lstsq(D[~te], y[~te], rcond=None)[0]
        rs.append(np.corrcoef(10 ** (D[te] @ b) - 1, 10 ** y[te] - 1)[0, 1])
    return np.array(rs)


flag = has.astype(float).where(has, 0.0)
ideas = [  # (group, label, frame, extra columns, sample rows)
    ('swap', 'peak flux swapped (raw GBM)', swap(['log10PeakFlux']), None),
    ('swap', 'peak flux swapped, offset removed', swap(['log10PeakFlux'], 'offset'), None),
    ('swap', 'peak flux = mean(Swift, GBM-offset)', swap(['log10PeakFlux'], 'avg'), None),
    ('swap', 'T90 swapped (raw)', swap(['log10T90']), None),
    ('swap', 'T90 swapped, offset removed', swap(['log10T90'], 'offset'), None),
    ('swap', 'fluence swapped (raw)', swap(['log10Fluence']), None),
    ('swap', 'fluence swapped, offset removed', swap(['log10Fluence'], 'offset'), None),
    ('swap', 'T90+fluence+PF swapped (raw)', swap(['log10T90', 'log10Fluence', 'log10PeakFlux']), None),
    ('swap', 'T90+fluence+PF swapped, offset removed', swap(['log10T90', 'log10Fluence', 'log10PeakFlux'], 'offset'), None),
    ('swap', 'T90+fluence+PF = mean(Swift, GBM-offset)', swap(['log10T90', 'log10Fluence', 'log10PeakFlux'], 'avg'), None),
    ('swap', 'T90+fluence+PF swapped raw + GBM-flag column', swap(['log10T90', 'log10Fluence', 'log10PeakFlux']), [flag]),
    ('swap', 'T90+fluence+PF+PhotonIndex swapped, offset removed',
     swap(['log10T90', 'log10Fluence', 'log10PeakFlux', 'PhotonIndex'], 'offset'), None),
    ('extra', 'Swift kept + GBM PF as extra column', d, [gbm.log10PeakFlux]),
    ('extra', 'Swift kept + GBM T90, fluence, PF as extra columns', d, [gbm.log10T90, gbm.log10Fluence, gbm.log10PeakFlux]),
    ('extra', '+ GBM Epeak (corrected match)', d, [log_ep]),
    ('extra', '+ GBM Epeak + Comptonized index', d, [log_ep, G.flnc_comp_index.where(ep_ok)]),
    ('extra', '+ GBM 64-ms peak flux', d, [lg(G.flux_64)]),
    ('extra', '+ GBM T50', d, [lg(G.t50)]),
    ('extra', '+ GBM peak-flux / 64-ms ratio (variability)', d, [lg(G.flux_64) - lg(G.flux_1024)]),
    ('features', 'drop log10NH', d, None),
    ('features', '+ squares of all features', d, None),
    ('combo', 'peak flux swapped (offset removed) + squares', swap(['log10PeakFlux'], 'offset'), None),
    ('combo', 'squares + all pairwise interactions', d, None),
]

rows = []
for feats_name, base_feats in (('SIMPLE7', SIMPLE7), ('ALL10', ALL10)):
    base = honest(d, base_feats)
    rows.append(dict(features=feats_name, group='baseline', idea='Swift only (baseline)', r=base.mean(),
                     diff=0.0, diff_sd=0.0, better_pct=np.nan, N=len(d)))
    for group, label, frame, extra in ideas:
        feats = list(base_feats)
        if 'fluence' in label.lower() and 'log10Fluence' not in feats and group == 'swap':
            continue  # fluence isn't a SIMPLE7 feature
        ex = extra
        if label == 'drop log10NH':
            feats = [f for f in feats if f != 'log10NH']
        if 'squares' in label:
            ex = [frame[f] ** 2 for f in feats]
        if 'interactions' in label:
            ex += [frame[a] * frame[b] for i, a in enumerate(feats) for b in feats[i + 1:]]
        r = honest(frame, feats, ex); dd = r - base
        rows.append(dict(features=feats_name, group=group, idea=label, r=r.mean(), diff=dd.mean(),
                         diff_sd=dd.std(), better_pct=100 * (dd > 0).mean(), N=len(d)))
    # Same GRBs only: Swift values vs GBM values for the shared columns
    sub = has
    b2 = honest(d, base_feats, rows=sub)
    for label, frame in (('overlap only: GBM T90+fluence+PF instead of Swift', swap(['log10T90', 'log10Fluence', 'log10PeakFlux'])),
                         ('overlap only: GBM PF instead of Swift', swap(['log10PeakFlux']))):
        r = honest(frame, base_feats, rows=sub); dd = r - b2
        rows.append(dict(features=feats_name, group='overlap', idea=label, r=r.mean(), diff=dd.mean(),
                         diff_sd=dd.std(), better_pct=100 * (dd > 0).mean(), N=int(sub.sum()),
                         base_r=b2.mean()))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, 'linear_ideas_check.csv'), index=False, float_format='%.4f')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
print('offsets GBM - Swift:', {k: round(v, 3) for k, v in off.items()}, '| GBM Epeak on', int(log_ep.notna().sum()))
print(out.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
