#!/usr/bin/env python3
"""
Part B -- does a Swift <-> Fermi shared latent space improve redshift
prediction for Swift GRBs?

Protocol = Spencer's honest held-out r (GRB_linear_model_comparison/
simple_linear_flagship_check.py): 100 random 80/20 splits of the 212 long
Swift X-ray GRBs, everything refit on the 80%, scored on the untouched 20%,
r in linear z. Test GRBs are always Swift GRBs predicted from their Swift
view only. "pooled" variants also put the 106 Fermi-only GRBs (GBM view + z)
into training -- they are never in a test set. Unsupervised latent fits
(CCA, PCCA) never see z; supervised PCCA sees z for training rows only.

Models
  B0  baselines: OLS on SIMPLE7, on all 9 Swift features, on Swift prompt only
  B1  CCA-family latent of the PROMPT variables (Swift prompt <-> GBM prompt,
      fit on the bridge GRBs inside the training split, k=3), then
      z ~ latent; Swift-only training vs pooled with Fermi-only GRBs
  B2  unsupervised PCCA (views: Swift9 | GBM prompt | Epeak), z ~ latent
  B3  supervised PCCA (views: Swift9 | GBM prompt | Epeak | z), predict the
      z view from the Swift view
  B4  generalized PCCA (views: Swift9 | optical plateau | GBM prompt | Epeak
      | z), predict z from the Swift view

Run from this folder: python3 part_b_redshift.py
"""
import warnings
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import cca_lib as L

warnings.filterwarnings('ignore')
N_SPLITS, SEED = 100, 1

# Hyperparameters for B1: the typical inner-CV choice in Part A (A1 pair).
B1_METHODS = {
    'CCA':          lambda X, Y: L.fit_cca(X, Y, 3)[:2],
    'ridge CCA':    lambda X, Y: L.fit_cca(X, Y, 3, r=0.9)[:2],
    'PLS':          lambda X, Y: L.fit_pls(X, Y, 3)[:2],
    'sparse CCA':   lambda X, Y: L.fit_sparse(X, Y, 3, cfrac=0.3)[:2],
    'robust (MCD)': lambda X, Y: L.fit_robust_mcd(X, Y, 3)[:2],
    'rank CCA':     None,   # built below (needs the transform)
    'kernel CCA':   lambda X, Y: L.fit_kcca(X, Y, 3, 0.5, 0.03)[:2],
}


def _rank_cca(X, Y):
    tx, ty = L.normal_scores_fit(X), L.normal_scores_fit(Y)
    gx, gy = L.fit_cca(tx(X), ty(Y), 3, r=0.9)[:2]
    return (lambda A: gx(tx(A))), (lambda B: gy(ty(B)))


B1_METHODS['rank CCA'] = _rank_cca


def r_lin(pred_log, obs_log):
    return np.corrcoef(10 ** pred_log - 1, 10 ** obs_log - 1)[0, 1]


def ols(A, y, At):
    D = np.column_stack([np.ones(len(A)), A])
    b = np.linalg.lstsq(D, y, rcond=None)[0]
    return np.column_stack([np.ones(len(At)), At]) @ b


def one_split(d, s):
    sw = d[d.source == 'swift']
    fo = d[d.source == 'fermi_only']
    n = len(sw)
    rng = np.random.default_rng(SEED + s)
    te = np.zeros(n, bool)
    te[rng.choice(n, size=int(0.2 * n), replace=False)] = True
    tr_df, te_df = sw[~te], sw[te]
    y_te = te_df.log10z.values
    out = {}

    # ---- B0 baselines
    for name, cols in (('B0 OLS SIMPLE7', L.SIMPLE7), ('B0 OLS Swift9', L.S_ALL),
                       ('B0 OLS Swift prompt4', L.S_PROMPT)):
        out[name] = r_lin(ols(tr_df[cols].values, tr_df.log10z.values, te_df[cols].values), y_te)

    # ---- B1 CCA-family prompt latent
    br = tr_df[tr_df.bridge]
    sx, sy = L.standardizer(br[L.S_PROMPT].values), L.standardizer(br[L.G_PROMPT].values)
    for m, fit in B1_METHODS.items():
        fx, fy = fit(sx(br[L.S_PROMPT].values), sy(br[L.G_PROMPT].values))
        Ztr = fx(sx(tr_df[L.S_PROMPT].values))
        Zte = fx(sx(te_df[L.S_PROMPT].values))
        Zfo = fy(sy(fo[L.G_PROMPT].values))
        out[f'B1 {m} | Swift-train'] = r_lin(ols(Ztr, tr_df.log10z.values, Zte), y_te)
        out[f'B1 {m} | pooled'] = r_lin(ols(np.vstack([Ztr, Zfo]),
                                            np.r_[tr_df.log10z.values, fo.log10z.values], Zte), y_te)

    # ---- PCCA family (B2-B4)
    def pcca_block(cols_views, pooled, k, supervised):
        rows = pd.concat([tr_df, fo]) if pooled else tr_df
        allcols = sum(cols_views, [])
        if supervised:
            allcols = allcols + ['log10z']
        st = L.standardizer(rows[allcols].values)
        Xtr = st(rows[allcols].values)
        views, i = [], 0
        for cv in cols_views + ([['log10z']] if supervised else []):
            views.append(np.arange(i, i + len(cv))); i += len(cv)
        m = L.PCCA(k, seed=s).fit(Xtr, views)
        Xte = st(np.column_stack([te_df[allcols[:-1] if supervised else allcols].values,
                                  np.zeros(len(te_df))]) if supervised else te_df[allcols].values)
        if supervised:
            zs = m.predict_view(Xte, [0], len(views) - 1)[:, 0]
            mu, sd = np.nanmean(rows.log10z), np.nanstd(rows.log10z)
            return r_lin(zs * sd + mu, y_te)
        # unsupervised: z ~ latent. Swift rows -> latent from Swift view only
        # (as at test time); Fermi-only rows -> latent from their GBM views.
        sw_rows = rows.source.values == 'swift'
        Lsw = m.latent(Xtr[sw_rows], [0])
        Lte = m.latent(Xte, [0])
        A, y = Lsw, rows.log10z.values[sw_rows]
        if pooled:
            fo_rows = ~sw_rows
            Xfo = Xtr[fo_rows]
            has_ep = ~np.isnan(Xfo[:, views[2]]).any(1)
            Lfo = np.zeros((fo_rows.sum(), k))
            Lfo[has_ep] = m.latent(Xfo[has_ep], [1, 2])
            Lfo[~has_ep] = m.latent(np.nan_to_num(Xfo[~has_ep]), [1])
            A, y = np.vstack([A, Lfo]), np.r_[y, rows.log10z.values[fo_rows]]
        return r_lin(ols(A, y, Lte), y_te)

    two = [L.S_ALL, L.G_PROMPT, ['g_log10Epeak']]
    four = [L.S_ALL, L.G_PROMPT, ['g_log10Epeak'], L.O_PLATEAU]
    for k in (1, 2, 3):
        for pooled in (False, True):
            tag = 'pooled' if pooled else 'Swift-train'
            out[f'B2 PCCA k={k} | {tag}'] = pcca_block(two, pooled, k, False)
            out[f'B3 supervised PCCA k={k} | {tag}'] = pcca_block(two, pooled, k, True)
        out[f'B4 generalized PCCA k={k} | pooled'] = pcca_block(four, True, k, True)
    return out


def main():
    d = L.load_all()
    res = Parallel(n_jobs=-1)(delayed(one_split)(d, s) for s in range(N_SPLITS))
    R = pd.DataFrame(res)
    base = R['B0 OLS SIMPLE7']
    summ = pd.DataFrame({
        'mean_r': R.mean(), 'sd_r': R.std(),
        'diff_vs_SIMPLE7': R.sub(base, axis=0).mean(),
        'diff_sd': R.sub(base, axis=0).std(),
        'beats_SIMPLE7_pct': (R.sub(base, axis=0) > 0).mean() * 100,
    })
    summ.to_csv('results_part_b.csv', float_format='%.4f')
    R.to_csv('results_part_b_per_split.csv', index=False, float_format='%.4f')
    pd.set_option('display.width', 160)
    print(summ.to_string(float_format=lambda v: f'{v:.3f}'))


if __name__ == '__main__':
    main()
