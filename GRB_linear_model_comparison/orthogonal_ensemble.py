#!/usr/bin/env python3
"""
An ORTHOGONAL ensemble: learners chosen for DECORRELATED errors, evaluated honestly.

The flagship "SuperLearner" is dominated by one member -- the winning GLM (weight
~0.62) -- which is also used for the M-est cut and was itself selected on the data.
That is the opposite of what makes an ensemble work: an ensemble only beats its
members when their errors are diverse (nearly orthogonal), so combining them
cancels error. Here we build genuinely different model families, MEASURE how
orthogonal their errors actually are, combine them, and ask the real question:

    does honest, non-circular diversity beat the winning GLM's honest r (~0.53),
    or does performance cap there -- i.e. is the ceiling the physics, not the model?

Learners (all predict log10(1+z); r is Pearson in linear z = 10**pred - 1):
  glm    : the winning GLM formula (OLS)          -- the flagship's dominant member
  lin7   : plain 7-feature additive OLS           -- simple parametric
  ridge  : L2 linear on 10 features (standardized)-- regularized parametric
  rf     : random forest                          -- bagged trees, nonlinear
  gbm    : hist gradient boosting                 -- boosted trees, nonlinear
  knn    : k-nearest neighbours (standardized)    -- local, non-parametric

Combination: non-negative least squares (NNLS) weights fit on inner out-of-fold
predictions of the TRAINING split only (SuperLearner-style), so the ensemble is
evaluated honestly on the untouched 20%. We also report a plain average.

Outputs the honest per-learner r, the residual-correlation matrix (the actual
"orthogonality"), the fitted ensemble weights, and the ensemble r.
Saves 'orthogonal_ensemble.png'.
"""
import warnings, numpy as np, pandas as pd, patsy
from scipy.optimize import nnls
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import KFold
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from simple_linear_flagship_check import d, WINNER

FEAT10 = ["log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
          "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux"]
FEAT7  = ["log10T90", "log10Fa", "log10Ta", "Alpha", "log10NH", "PhotonIndex", "log10PeakFlux"]
i7     = [FEAT10.index(c) for c in FEAT7]

y     = d["log10z"].to_numpy()
Xbase = d[FEAT10].to_numpy()
Xwin  = np.asarray(patsy.dmatrix(WINNER.split("~")[1], d, return_type="matrix"))   # winning design + intercept
Xlin7 = np.column_stack([np.ones(len(y)), Xbase[:, i7]])
r_lin = lambda p, o: np.corrcoef(10**p - 1, 10**o - 1)[0, 1]

# ---- learners as fit/predict closures over row indices ----------------------
def _lstsq(Xd):
    def f(tr, te):
        b = np.linalg.lstsq(Xd[tr], y[tr], rcond=None)[0]
        return Xd[te] @ b
    return f
def _sk(builder):
    def f(tr, te):
        m = builder().fit(Xbase[tr], y[tr])
        return m.predict(Xbase[te])
    return f

LEARNERS = {
    "glm":   _lstsq(Xwin),
    "lin7":  _lstsq(Xlin7),
    "ridge": _sk(lambda: make_pipeline(StandardScaler(), Ridge(alpha=10.0))),
    "rf":    _sk(lambda: RandomForestRegressor(n_estimators=300, min_samples_leaf=3, random_state=0, n_jobs=-1)),
    "gbm":   _sk(lambda: HistGradientBoostingRegressor(max_depth=3, learning_rate=0.05, max_iter=400, random_state=0)),
    "knn":   _sk(lambda: make_pipeline(StandardScaler(), KNeighborsRegressor(n_neighbors=15))),
}
names = list(LEARNERS)

# ---- one honest out-of-fold pass on the FULL sample (for residual corr) ------
def oof_full(K=10, seed=0):
    P = {k: np.full(len(y), np.nan) for k in names}
    for tr, te in KFold(K, shuffle=True, random_state=seed).split(y):
        for k in names:
            P[k][te] = LEARNERS[k](tr, te)
    return P
Pf = oof_full()
resid = pd.DataFrame({k: Pf[k] - y for k in names})           # out-of-fold residuals
Rmat  = resid.corr().to_numpy()                               # the "orthogonality" matrix
oof_r = {k: r_lin(Pf[k], y) for k in names}

# ---- honest 80/20 splits: per-learner r + NNLS ensemble + average -----------
def honest(n_splits=30, seed=1):
    rng = np.random.default_rng(seed); n = len(y)
    per = {k: [] for k in names}; ens_nnls = []; ens_avg = []; W = []
    for _ in range(n_splits):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.2 * n), replace=False)] = True
        tr = np.where(~te)[0]; teI = np.where(te)[0]
        # inner OOF on the training rows -> matrix to fit ensemble weights
        Poof = np.zeros((len(tr), len(names)))
        for j, (itr, ite) in enumerate(KFold(5, shuffle=True, random_state=1).split(tr)):
            for c, k in enumerate(names):
                Poof[ite, c] = LEARNERS[k](tr[itr], tr[ite])
        w, _ = nnls(Poof, y[tr]); W.append(w)
        # base learners refit on full train, predict test
        Pte = np.column_stack([LEARNERS[k](tr, teI) for k in names])
        for c, k in enumerate(names):
            per[k].append(r_lin(Pte[:, c], y[teI]))
        ens_nnls.append(r_lin(Pte @ w, y[teI]))
        ens_avg.append(r_lin(Pte.mean(1), y[teI]))
    return ({k: float(np.mean(v)) for k, v in per.items()},
            float(np.mean(ens_nnls)), float(np.mean(ens_avg)), np.mean(W, 0))

hon, r_nnls, r_avg, wbar = honest()

# ---- report -----------------------------------------------------------------
off = Rmat[~np.eye(len(names), dtype=bool)]
print("Residual-correlation matrix (out-of-fold, want off-diagonal near 0):")
print(pd.DataFrame(Rmat, index=names, columns=names).round(2))
print(f"\nmean |off-diagonal residual correlation| = {np.abs(off).mean():.2f}  "
      f"(1.0 = identical errors, 0 = orthogonal)\n")
print(f"{'learner':8s} {'honest r':>9s}   NNLS weight")
for k in names:
    print(f"{k:8s} {hon[k]:9.3f}   {wbar[names.index(k)]:.2f}")
print(f"\n{'ENSEMBLE (NNLS)':8s} honest r = {r_nnls:.3f}")
print(f"{'ENSEMBLE (avg) ':8s} honest r = {r_avg:.3f}")
print(f"\nbest single learner honest r = {max(hon.values()):.3f} ({max(hon, key=hon.get)})")
print(f"winning GLM honest r          = {hon['glm']:.3f}")

# ---- figure -----------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
im = ax[0].imshow(Rmat, vmin=0, vmax=1, cmap="RdYlBu_r")
ax[0].set_xticks(range(len(names))); ax[0].set_xticklabels(names)
ax[0].set_yticks(range(len(names))); ax[0].set_yticklabels(names)
for i in range(len(names)):
    for j in range(len(names)):
        ax[0].text(j, i, f"{Rmat[i,j]:.2f}", ha="center", va="center",
                   color="white" if abs(Rmat[i,j]) > .6 else "black", fontsize=9)
ax[0].set_title("Residual-correlation matrix\n(off-diagonal near 0 = orthogonal errors)")
fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04, label="corr of out-of-fold residuals")

order = names + ["ENS(avg)", "ENS(NNLS)"]
vals  = [hon[k] for k in names] + [r_avg, r_nnls]
cols  = ["#9ec3e6"]*len(names) + ["#f4a582", "#d7301f"]
cols[names.index("glm")] = "#4d4d4d"
ax[1].bar(range(len(order)), vals, color=cols)
for i, v in enumerate(vals): ax[1].text(i, v+.005, f"{v:.3f}", ha="center", fontsize=8.5)
ax[1].axhline(hon["glm"], color="#4d4d4d", ls="--", lw=1, alpha=.7)
ax[1].set_xticks(range(len(order))); ax[1].set_xticklabels(order, rotation=30, ha="right")
ax[1].set_ylabel("honest held-out r (avg 30x 80/20)"); ax[1].set_ylim(0.30, 0.62)
ax[1].set_title("Honest r: base learners vs orthogonal ensemble\n(dashed = winning GLM level)")
ax[1].grid(axis="y", alpha=0.25)
fig.tight_layout()
fig.savefig("orthogonal_ensemble.png", dpi=150)
print("\nwrote orthogonal_ensemble.png")
