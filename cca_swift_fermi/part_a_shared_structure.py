#!/usr/bin/env python3
"""
Part A -- is there real shared structure between the Swift and Fermi-GBM
views of the bridge GRBs, and which CCA variant finds it best?

For each view pair and each method:
  * held-out canonical correlation of components 1 and 2: repeated 5-fold CV
    (10 repeats); the method (and its hyperparameter, via inner 3-fold CV) is
    fit on the training folds only and the held-out GRBs are projected.
  * permutation p-value for component 1: the GBM rows are shuffled relative to
    the Swift rows and the whole CV (tuning included) is rerun; p = share of
    shuffles that match or beat the real held-out correlation.
  * held-out R^2 for predicting log Epeak from the Swift view (view pairs that
    include Epeak): OLS on the method's top-2 Swift canonical scores, vs plain
    OLS on all Swift features and vs PCCA's conditional mean.

Run from this folder: python3 part_a_shared_structure.py
"""
import warnings
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import cca_lib as L

warnings.filterwarnings('ignore')
K, N_REP, N_PERM = 2, 10, 100

VIEW_PAIRS = {
    'A1 prompt <-> GBM prompt':   (L.S_PROMPT, L.G_PROMPT),
    'A2 Swift all <-> GBM+Epeak': (L.S_ALL, L.G_ALL),
    'A3 plateau <-> GBM+Epeak':   (L.S_PLATEAU + ['log10NH'], L.G_ALL),
}

GRIDS = {
    'CCA':          [None],
    'ridge CCA':    [0.1, 0.3, 0.6, 0.9],
    'PLS':          [None],
    'sparse CCA':   [0.3, 0.5, 0.7, 1.0],
    'robust (MCD)': [None],
    'rank CCA':     [0.1, 0.3, 0.6, 0.9],
    'kernel CCA':   [(g, kap) for g in (0.5, 1.0) for kap in (0.03, 0.1, 0.3)],
}


def fitter(method, h):
    """-> function(Xtr, Ytr) returning (fx, fy) that take RAW standardized rows."""
    def f(X, Y):
        if method == 'CCA':          return L.fit_cca(X, Y, K)[:2]
        if method == 'ridge CCA':    return L.fit_cca(X, Y, K, r=h)[:2]
        if method == 'PLS':          return L.fit_pls(X, Y, K)[:2]
        if method == 'sparse CCA':   return L.fit_sparse(X, Y, K, cfrac=h)[:2]
        if method == 'robust (MCD)': return L.fit_robust_mcd(X, Y, K)[:2]
        if method == 'kernel CCA':   return L.fit_kcca(X, Y, K, *h)[:2]
        if method == 'rank CCA':
            tx, ty = L.normal_scores_fit(X), L.normal_scores_fit(Y)
            gx, gy = L.fit_cca(tx(X), ty(Y), K, r=h)[:2]
            return (lambda A: gx(tx(A))), (lambda B: gy(ty(B)))
    return f


def heldout_corr(X, Y, fit, folds):
    """Mean over folds of corr(u_k, v_k) on the held-out rows."""
    out = []
    for te in folds:
        tr = ~te
        sx, sy = L.standardizer(X[tr]), L.standardizer(Y[tr])
        fx, fy = fit(sx(X[tr]), sy(Y[tr]))
        u, v = fx(sx(X[te])), fy(sy(Y[te]))
        out.append([np.corrcoef(u[:, j], v[:, j])[0, 1] for j in range(K)])
    return np.nanmean(out, 0)


def kfold(n, k, rng):
    lab = rng.permutation(np.arange(n) % k)
    return [lab == i for i in range(k)]


def tuned_fit(method, Xtr, Ytr, rng):
    grid = GRIDS[method]
    if len(grid) == 1:
        return fitter(method, grid[0]), grid[0]
    inner = kfold(len(Xtr), 3, rng)
    scores = [heldout_corr(Xtr, Ytr, fitter(method, h), inner)[0] for h in grid]
    h = grid[int(np.nanargmax(scores))]
    return fitter(method, h), h


def cv_eval(X, Y, method, seed):
    """Repeated 5-fold CV with inner tuning; returns (comp1, comp2) means and
    the chosen hyperparameters."""
    rng = np.random.default_rng(seed)
    res, chosen = [], []
    for _ in range(N_REP):
        for te in kfold(len(X), 5, rng):
            tr = ~te
            sx, sy = L.standardizer(X[tr]), L.standardizer(Y[tr])
            Xs, Ys = sx(X[tr]), sy(Y[tr])
            f, h = tuned_fit(method, Xs, Ys, rng)
            fx, fy = f(Xs, Ys)
            u, v = fx(sx(X[te])), fy(sy(Y[te]))
            res.append([np.corrcoef(u[:, j], v[:, j])[0, 1] for j in range(K)])
            chosen.append(h)
    return np.nanmean(res, 0), chosen


def epeak_r2(X, Y, method, seed):
    """Held-out R^2 for log Epeak (last Y column) from the Swift view."""
    rng = np.random.default_rng(seed)
    r2 = []
    for _ in range(N_REP):
        pred = np.full(len(X), np.nan)
        for te in kfold(len(X), 5, rng):
            tr = ~te
            sx, sy = L.standardizer(X[tr]), L.standardizer(Y[tr])
            Xs, Ys = sx(X[tr]), sy(Y[tr])
            y = Y[tr, -1]
            if method == 'OLS (all Swift)':
                A, At = Xs, sx(X[te])
            elif method.startswith('PCCA'):
                k = int(method[-1])
                m = L.PCCA(k).fit(np.hstack([Xs, Ys]),
                                              [np.arange(X.shape[1]), X.shape[1] + np.arange(Y.shape[1])])
                Zt = np.hstack([sx(X[te]), np.full((te.sum(), Y.shape[1]), np.nan)])
                ys = m.predict_view(np.nan_to_num(Zt), [0], 1)[:, -1]
                pred[te] = ys * np.nanstd(Y[tr, -1]) + np.nanmean(Y[tr, -1])
                continue
            else:
                f, _ = tuned_fit(method, Xs, Ys, rng)
                fx, _ = f(Xs, Ys)
                A, At = fx(Xs), fx(sx(X[te]))
            D = np.column_stack([np.ones(len(A)), A])
            beta = np.linalg.lstsq(D, y, rcond=None)[0]
            pred[te] = np.column_stack([np.ones(len(At)), At]) @ beta
        r2.append(1 - np.sum((Y[:, -1] - pred) ** 2) / np.sum((Y[:, -1] - Y[:, -1].mean()) ** 2))
    return float(np.mean(r2))


def main():
    d = L.load_all()
    rows = []
    for pair, (sc, gc) in VIEW_PAIRS.items():
        b = d[d.bridge].dropna(subset=sc + gc)
        X, Y = b[sc].values, b[gc].values
        print(f'\n=== {pair}   (N={len(b)}, {X.shape[1]} vs {Y.shape[1]} variables)')
        for method in GRIDS:
            (c1, c2), chosen = cv_eval(X, Y, method, seed=0)
            null = Parallel(n_jobs=-1)(delayed(
                lambda s: cv_eval(X, Y[np.random.default_rng(1000 + s).permutation(len(Y))],
                                  method, seed=s)[0][0])(s) for s in range(N_PERM))
            p = (1 + np.sum(np.array(null) >= c1)) / (1 + N_PERM)
            hp = pd.Series([str(h) for h in chosen]).mode()[0]
            row = dict(pair=pair, method=method, N=len(b), heldout_r1=c1, heldout_r2=c2,
                       null_r1_mean=float(np.mean(null)), null_r1_95=float(np.quantile(null, 0.95)),
                       perm_p=p, typical_hparam=hp)
            if 'Epeak' in pair:
                row['epeak_R2'] = epeak_r2(X, Y, method, seed=0)
            rows.append(row)
            print(f'  {method:13s} held-out r1={c1:.3f} r2={c2:.3f} | null r1 mean={np.mean(null):.3f}'
                  f' 95%={np.quantile(null, .95):.3f}  p={p:.3f}  hparam~{hp}'
                  + (f' | Epeak R2={row["epeak_R2"]:.3f}' if 'epeak_R2' in row else ''))
        if 'Epeak' in pair:
            for m in ('OLS (all Swift)', 'PCCA k=1', 'PCCA k=2'):
                r2 = epeak_r2(X, Y, m, seed=0)
                rows.append(dict(pair=pair, method=m, N=len(b), epeak_R2=r2))
                print(f'  {m:13s} Epeak R2={r2:.3f}')
    pd.DataFrame(rows).to_csv('results_part_a.csv', index=False, float_format='%.4f')
    print('\nWrote results_part_a.csv')


if __name__ == '__main__':
    main()
