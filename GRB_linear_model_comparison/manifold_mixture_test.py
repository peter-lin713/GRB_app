#!/usr/bin/env python3
"""
Test the mixture-of-manifolds idea HONESTLY.

Hypothesis (correct): distinct data manifolds exist in the joint (features, z)
space even though they are invisible in feature space alone -- so a single
function fit to all of them must do poorly. Fix: partition the TRAINING data in
the (flux, z) plane, fit a linear model per partition ("expert").

The catch this script measures: to score a NEW point you must route it to a
partition, and the partitions were defined using z -- the very thing you are
predicting. So we evaluate every mixture two ways, over repeated 80/20 splits:

  HONEST  : partitions are defined on the training z; a gate g(x) is trained to
            predict the partition from FEATURES ONLY; test points are routed by
            g(x_test). This is what you could actually deploy.
  ORACLE  : test points are routed to the expert whose prediction is closest to
            their TRUE z. This is an upper bound that secretly uses the answer --
            i.e. the leakage a z-based partition smuggles in.

K = number of manifolds is swept as a hyperparameter. If HONEST r rises with K,
the manifolds are real and usable. If only ORACLE rises while HONEST stays flat,
the manifold structure is real but NOT recoverable from x -> not deployable, and
any CV that reports a gain is routing on the target.

Partitions: KMeans on standardized (X, z) of the training rows (uses z, as the
idea requires). Experts: OLS of log10(1+z) on the 10 features within a partition,
falling back to the global fit for tiny partitions. r is Pearson in linear z.
Saves 'manifold_mixture_test.png'.
"""
import warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from simple_linear_flagship_check import d

FEAT = ["log10T90","log10Fa","log10Ta","Alpha","Beta","Gamma",
        "log10Fluence","PhotonIndex","log10NH","log10PeakFlux"]
X = d[FEAT].to_numpy(); y = d["log10z"].to_numpy(); n = len(y)
X1 = np.column_stack([np.ones(n), X])                       # design with intercept
r_lin = lambda p, o: np.corrcoef(10**p - 1, 10**o - 1)[0, 1]

def fit_expert(rows):
    if rows.sum() < 12:                                     # too small -> global fit
        b, *_ = np.linalg.lstsq(X1, y, rcond=None); return b
    b, *_ = np.linalg.lstsq(X1[rows], y[rows], rcond=None); return b

def run(K, n_splits=20, seed=0):
    rng = np.random.default_rng(seed); honest, oracle = [], []
    for _ in range(n_splits):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.2*n), replace=False)] = True
        tr = ~te; itr, ite = np.where(tr)[0], np.where(te)[0]
        if K == 1:
            b = fit_expert(tr)
            honest.append(r_lin(X1[ite] @ b, y[ite])); oracle.append(honest[-1]); continue
        # partition TRAIN in (X, z) space (uses z, as the idea requires)
        Z = StandardScaler().fit_transform(np.column_stack([X[itr], y[itr]]))
        lab = KMeans(K, n_init=5, random_state=0).fit_predict(Z)
        experts = [fit_expert_idx(itr[lab == k]) for k in range(K)]
        # gate: features -> partition label (deployable), trained on train only
        gate = RandomForestClassifier(n_estimators=200, random_state=0, n_jobs=-1).fit(X[itr], lab)
        g = gate.predict(X[ite])
        pred_h = np.array([X1[i] @ experts[g[j]] for j, i in enumerate(ite)])
        honest.append(r_lin(pred_h, y[ite]))
        # oracle: route each test point to the expert closest to its TRUE y
        P = np.column_stack([X1[ite] @ e for e in experts])
        pick = np.argmin(np.abs(P - y[ite][:, None]), axis=1)
        pred_o = P[np.arange(len(ite)), pick]
        oracle.append(r_lin(pred_o, y[ite]))
    return np.mean(honest), np.std(honest), np.mean(oracle), np.std(oracle)

def fit_expert_idx(idx):
    if len(idx) < 12:
        b, *_ = np.linalg.lstsq(X1, y, rcond=None); return b
    b, *_ = np.linalg.lstsq(X1[idx], y[idx], rcond=None); return b

Ks = [1, 2, 3, 4, 5]
res = {K: run(K) for K in Ks}
print(f"{'K':>2s}  {'HONEST r (gate on x)':>22s}  {'ORACLE r (routes on true z)':>28s}")
for K in Ks:
    h, hs, o, os_ = res[K]
    print(f"{K:>2d}  {h:>10.3f} ± {hs:.3f}        {o:>10.3f} ± {os_:.3f}")

fig, ax = plt.subplots(figsize=(8.5, 6))
h  = [res[K][0] for K in Ks]; hs = [res[K][1] for K in Ks]
o  = [res[K][2] for K in Ks]; os_= [res[K][3] for K in Ks]
ax.errorbar(Ks, h, yerr=hs, marker="o", lw=2, color="#2c7fb8", capsize=3,
            label="HONEST: route test points by a gate on x  (deployable)")
ax.errorbar(Ks, o, yerr=os_, marker="s", lw=2, color="#d7301f", capsize=3, ls="--",
            label="ORACLE: route test points by their true z  (leakage)")
ax.axhline(h[0], color="gray", ls=":", lw=1, label="single linear model (K=1)")
ax.set_xlabel("K  (number of manifolds / experts)"); ax.set_ylabel("held-out r (linear z), avg 20 splits")
ax.set_title("Mixture of z-partitioned linear experts: honest vs oracle routing\n"
             "manifolds are real, but membership is not recoverable from x")
ax.set_xticks(Ks); ax.legend(fontsize=9, loc="center right"); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig("manifold_mixture_test.png", dpi=150)
print("\nwrote manifold_mixture_test.png")
