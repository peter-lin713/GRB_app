"""
Data assembly + CCA-family implementations for the Swift <-> Fermi-GBM
shared-latent experiments. Everything is plain numpy / scipy / sklearn (no
CCA package is installed), small enough to read in one sitting.

All methods take TRAIN-standardized matrices and return projection functions,
so every fit can be done strictly inside a training split.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.covariance import MinCovDet

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

# ---------------------------------------------------------------- data ------

# Swift view: BAT prompt (4) + XRT plateau (4) + XRT NH
S_PROMPT  = ['log10T90', 'log10Fluence', 'PhotonIndex', 'log10PeakFlux']
S_PLATEAU = ['log10Fa', 'log10Ta', 'Alpha', 'Beta']
S_ALL     = S_PROMPT + S_PLATEAU + ['log10NH']
# GBM view, cutoff-power-law ("comp") fits over the fluence interval -- the
# same columns fermi_shared_vars.csv was built from (verified value-for-value).
G_PROMPT  = ['g_log10T90', 'g_log10Fluence', 'g_Index', 'g_log10PhtFlux']
G_ALL     = G_PROMPT + ['g_log10Epeak']
O_PLATEAU = ['o_log10Fa', 'o_log10Ta', 'o_Alpha', 'o_Beta']
# Spencer's SIMPLE7 (GRB_linear_model_comparison/simple_linear_flagship_check.py)
SIMPLE7   = ['log10T90', 'log10Fa', 'log10Ta', 'Alpha', 'log10NH', 'PhotonIndex', 'log10PeakFlux']


def _norm_id(g):
    g = str(g).replace('GRB', '')
    return g if g[-1].isalpha() else g + 'A'


def _good_epeak(ep, ep_pos_err):
    """Spencer's well-constrained-Epeak filter (transfer_learning_analysis.py)."""
    ok = (ep > 8) & (ep < 1e4) & (ep_pos_err < ep)
    return np.where(ok, np.log10(ep), np.nan)


def load_all():
    """One row per GRB: long Swift X-ray GRBs (+ GBM view for the bridge
    GRBs, + optical plateau where the optical catalog has it) followed by the
    Fermi-only GRBs (GBM view only). Columns: S_ALL, G_ALL, O_PLATEAU, z,
    source in {'swift', 'fermi_only'}, bridge (bool)."""
    x = pd.read_csv(os.path.join(REPO, 'Data',
                    'Xray_data_with_redshift_V8-web-app_processed_filtered_MICE (1).csv'))
    sw = pd.DataFrame({
        'id': x.GRB_Name.map(_norm_id), 'z': x.Redshift_crosscheck,
        'log10T90': x.log10T90, 'log10Fluence': x.log10Fluence,
        'PhotonIndex': x.PhotonIndex, 'log10PeakFlux': x.log10PeakFlux,
        'log10Fa': x.log10Fa_X, 'log10Ta': x.log10Ta_X,
        'Alpha': x.Alpha_X, 'Beta': x.Beta_X, 'log10NH': x.log10NH,
    })
    sw = sw[sw.log10T90 > np.log10(2)].set_index('id')

    c = pd.read_csv(os.path.join(HERE, 'data', 'common67_transfer.csv'))
    num = lambda s: pd.to_numeric(c[s], errors='coerce').values
    gb = pd.DataFrame({
        'g_log10T90': np.log10(num('t90')),
        'g_log10Fluence': np.log10(num('flnc_comp_ergflnc')),
        'g_Index': num('flnc_comp_index'),
        'g_log10PhtFlux': np.log10(num('flnc_comp_phtflux')),
        'g_log10Epeak': _good_epeak(num('flnc_comp_epeak'), num('flnc_comp_epeak_pos_err')),
    }, index=c.GRB_Name.map(_norm_id))
    sw = sw.join(gb)

    o = pd.read_csv(os.path.join(REPO, 'Data',
                    'OnlyLGRBs_data_171_optical_processed_error-cut_MICE (1).csv'))
    op = pd.DataFrame({'o_log10Fa': o.log10Faopt.values, 'o_log10Ta': o.log10Taopt.values,
                       'o_Alpha': o.Alpha_opt.values, 'o_Beta': o.Beta_opt.values},
                      index=o.GRB.map(_norm_id))
    op = op[~op.index.duplicated()]
    sw = sw.join(op)
    sw['source'] = 'swift'
    sw['bridge'] = sw[G_PROMPT].notna().all(axis=1)

    f = pd.read_csv(os.path.join(HERE, 'data', 'fermi_shared_vars.csv'))
    e = pd.read_csv(os.path.join(HERE, 'data', 'epeak_z_163.csv'))
    f = f[~f.is_common & (f.gbm_log10T90 > np.log10(2))].merge(
        e[['Fermi_name', 'flnc_comp_epeak']], on='Fermi_name', how='left')
    fo = pd.DataFrame({
        'z': f.z.values,
        'g_log10T90': f.gbm_log10T90.values, 'g_log10Fluence': f.gbm_log10fluence.values,
        'g_Index': f.gbm_photindex.values, 'g_log10PhtFlux': f.gbm_log10peakflux.values,
        # epeak_z_163 carries no Epeak errors; apply the value-range part of the filter
        'g_log10Epeak': np.where((f.flnc_comp_epeak > 8) & (f.flnc_comp_epeak < 1e4),
                                 np.log10(f.flnc_comp_epeak), np.nan),
    }, index=f.Fermi_name.values)
    fo['source'] = 'fermi_only'
    fo['bridge'] = False

    d = pd.concat([sw, fo], axis=0)
    d['log10z'] = np.log10(d.z + 1)
    return d


# ------------------------------------------------------------ helpers -------

def standardizer(A):
    mu, sd = np.nanmean(A, 0), np.nanstd(A, 0)
    sd = np.where(sd > 0, sd, 1.0)
    return lambda B: (B - mu) / sd


def inv_sqrt(S, eps=1e-8):
    w, V = np.linalg.eigh((S + S.T) / 2)
    return V @ np.diag(1 / np.sqrt(np.maximum(w, eps))) @ V.T


def normal_scores_fit(A):
    """Rank -> normal-score transform fitted on train columns, applied to new
    rows by interpolating against the train order statistics."""
    n = len(A)
    srt = np.sort(A, 0)
    q = norm.ppf((np.arange(1, n + 1) - 0.5) / n)
    def tf(B):
        return np.column_stack([np.interp(B[:, j], srt[:, j], q) for j in range(A.shape[1])])
    return tf


# ---------------------------------------------- linear CCA family -----------
# Every fitter returns (fx, fy, rho): fx(X) / fy(Y) project new standardized
# rows onto the k canonical directions.

def cca_from_cov(Sxx, Syy, Sxy, k, rx=0.0, ry=0.0):
    px, py = len(Sxx), len(Syy)
    Kx = inv_sqrt((1 - rx) * Sxx + rx * np.eye(px))
    Ky = inv_sqrt((1 - ry) * Syy + ry * np.eye(py))
    U, s, Vt = np.linalg.svd(Kx @ Sxy @ Ky)
    Wx, Wy = Kx @ U[:, :k], Ky @ Vt.T[:, :k]
    return (lambda X: X @ Wx), (lambda Y: Y @ Wy), s[:k], Wx, Wy


def fit_cca(X, Y, k, r=0.0):
    n = len(X)
    C = np.cov(np.hstack([X, Y]).T)
    px = X.shape[1]
    return cca_from_cov(C[:px, :px], C[px:, px:], C[:px, px:], k, r, r)


def fit_pls(X, Y, k):
    """PLS-SVD: maximize covariance (ridge CCA with r = 1)."""
    return fit_cca(X, Y, k, r=1.0)


def fit_robust_mcd(X, Y, k, r=0.05):
    C = MinCovDet(random_state=0).fit(np.hstack([X, Y])).covariance_
    px = X.shape[1]
    return cca_from_cov(C[:px, :px], C[px:, px:], C[:px, px:], k, r, r)


def _soft_l1(a, c):
    """Soft-threshold a so that ||u||_2 = 1 and ||u||_1 <= c (Witten et al. 2009)."""
    u = a / (np.linalg.norm(a) + 1e-12)
    if np.abs(u).sum() <= c:
        return u
    lo, hi = 0.0, np.abs(a).max()
    for _ in range(60):
        lam = (lo + hi) / 2
        s = np.sign(a) * np.maximum(np.abs(a) - lam, 0)
        s = s / (np.linalg.norm(s) + 1e-12)
        lo, hi = (lam, hi) if np.abs(s).sum() > c else (lo, lam)
    s = np.sign(a) * np.maximum(np.abs(a) - hi, 0)
    return s / (np.linalg.norm(s) + 1e-12)


def fit_sparse(X, Y, k, cfrac=0.5):
    """Sparse CCA via penalized matrix decomposition of X'Y (Witten,
    Tibshirani & Hastie 2009), rank-k by deflation."""
    Z = X.T @ Y / len(X)
    cx = max(1.0, cfrac * np.sqrt(X.shape[1]))
    cy = max(1.0, cfrac * np.sqrt(Y.shape[1]))
    Wx, Wy, ds = [], [], []
    for _ in range(k):
        v = np.linalg.svd(Z)[2][0]
        for _ in range(200):
            u = _soft_l1(Z @ v, cx)
            v_new = _soft_l1(Z.T @ u, cy)
            if np.allclose(v_new, v, atol=1e-8):
                break
            v = v_new
        d = u @ Z @ v
        Z = Z - d * np.outer(u, v)
        Wx.append(u); Wy.append(v); ds.append(d)
    Wx, Wy = np.array(Wx).T, np.array(Wy).T
    return (lambda A: A @ Wx), (lambda B: B @ Wy), np.array(ds), Wx, Wy


# -------------------------------------------------------- kernel CCA --------

def _rbf(A, B, g):
    d2 = (A ** 2).sum(1)[:, None] + (B ** 2).sum(1)[None, :] - 2 * A @ B.T
    return np.exp(-g * np.maximum(d2, 0))


def fit_kcca(X, Y, k, gfrac=1.0, kappa=0.1):
    """Regularized RBF kernel CCA (Hardoon, Szedmak & Shawe-Taylor 2004)."""
    n = len(X)
    gx, gy = gfrac / X.shape[1], gfrac / Y.shape[1]
    H = np.eye(n) - 1.0 / n
    Kx0, Ky0 = _rbf(X, X, gx), _rbf(Y, Y, gy)
    Kx, Ky = H @ Kx0 @ H, H @ Ky0 @ H
    Lx = Kx @ Kx + kappa * n * Kx + 1e-6 * np.eye(n)
    Ly = Ky @ Ky + kappa * n * Ky + 1e-6 * np.eye(n)
    Ax, Ay = inv_sqrt(Lx), inv_sqrt(Ly)
    U, s, Vt = np.linalg.svd(Ax @ Kx @ Ky @ Ay)
    a, b = Ax @ U[:, :k], Ay @ Vt.T[:, :k]

    def centered(Kt, K0):
        m = Kt.shape[0]
        return Kt - Kt.mean(1, keepdims=True) - K0.mean(0)[None, :] + K0.mean()
    fx = lambda A: centered(_rbf(A, X, gx), Kx0) @ a
    fy = lambda B: centered(_rbf(B, Y, gy), Ky0) @ b
    return fx, fy, s[:k], None, None


# ------------------------------------- probabilistic CCA (multi-view) -------

class PCCA:
    """Probabilistic CCA / group factor analysis with full within-view noise
    covariance (Bach & Jordan 2005), fitted by EM with whole views allowed to
    be missing per row. x_i = W z_i + mu + eps,  z ~ N(0, I_k),
    eps ~ N(0, blockdiag(Psi_1..Psi_V)). With block-full noise the latent only
    captures covariance SHARED across views. Two views reduce to classical CCA;
    3+ views = probabilistic generalized CCA."""

    def __init__(self, k, n_iter=20000, ridge=1e-3, seed=0, tol=1e-7):
        # n_iter is a cap; EM stops once no loading moves by more than tol.
        # (k=1 fits can need several thousand iterations -- a fixed 200 left
        # predictions off by >1 sd.)
        self.k, self.n_iter, self.ridge, self.seed, self.tol = k, n_iter, ridge, seed, tol

    def fit(self, Xfull, views):
        """Xfull: n x p with NaN for missing; views: list of column-index arrays.
        A view counts as observed for a row only if all its columns are."""
        n, p = Xfull.shape
        k = self.k
        self.views = views
        obs = np.column_stack([~np.isnan(Xfull[:, v]).any(1) for v in views])
        rng = np.random.default_rng(self.seed)
        mu = np.nanmean(Xfull, 0)
        W = rng.normal(scale=0.1, size=(p, k))
        Psi = [np.diag(np.nanvar(Xfull[:, v], 0) + self.ridge) for v in views]
        X0 = np.where(np.isnan(Xfull), 0.0, Xfull)
        pats, inv = np.unique(obs, axis=0, return_inverse=True)
        inv = inv.ravel()
        for it in range(self.n_iter):
            W_prev = W.copy()
            Ez = np.zeros((n, k)); Ezz = np.zeros((n, k, k))
            for pi, pat in enumerate(pats):
                rows = np.where(inv == pi)[0]
                cols = np.concatenate([views[v] for v in np.where(pat)[0]]) if pat.any() else np.array([], int)
                if len(cols) == 0:
                    Ezz[rows] = np.eye(k); continue
                Pinv = np.linalg.inv(_blockdiag([Psi[v] for v in np.where(pat)[0]]))
                Wo = W[cols]
                M = np.linalg.inv(np.eye(k) + Wo.T @ Pinv @ Wo)
                Ez[rows] = (X0[np.ix_(rows, cols)] - mu[cols]) @ Pinv @ Wo @ M
                Ezz[rows] = M[None] + Ez[rows][:, :, None] * Ez[rows][:, None, :]
            for vi, v in enumerate(views):
                r = obs[:, vi]
                zt = np.hstack([Ez[r], np.ones((r.sum(), 1))])
                Zt = np.zeros((k + 1, k + 1))
                Zt[:k, :k] = Ezz[r].sum(0); Zt[:k, k] = Zt[k, :k] = Ez[r].sum(0); Zt[k, k] = r.sum()
                B = (X0[np.ix_(r, v)].T @ zt) @ np.linalg.inv(Zt)
                W[v], mu[v] = B[:, :k], B[:, k]
                Xr = X0[np.ix_(r, v)]
                P = (Xr.T @ Xr - B @ zt.T @ Xr) / r.sum()
                Psi[vi] = (P + P.T) / 2 + self.ridge * np.eye(len(v))
            if np.max(np.abs(W - W_prev)) < self.tol:
                break
        self.W, self.mu, self.Psi, self.n_iter_run = W, mu, Psi, it + 1
        return self

    def latent(self, Xrows, view_ids):
        """Posterior mean of z given only the listed views."""
        cols = np.concatenate([self.views[v] for v in view_ids])
        Pinv = np.linalg.inv(_blockdiag([self.Psi[v] for v in view_ids]))
        Wo = self.W[cols]
        M = np.linalg.inv(np.eye(self.k) + Wo.T @ Pinv @ Wo)
        return (Xrows[:, cols] - self.mu[cols]) @ Pinv @ Wo @ M

    def predict_view(self, Xrows, from_views, to_view):
        z = self.latent(Xrows, from_views)
        v = self.views[to_view]
        return z @ self.W[v].T + self.mu[v]


def _blockdiag(blocks):
    n = sum(len(b) for b in blocks)
    out = np.zeros((n, n)); i = 0
    for b in blocks:
        m = len(b); out[i:i + m, i:i + m] = b; i += m
    return out
