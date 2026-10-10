#!/usr/bin/env python3
"""
Standard true-vs-predicted redshift, two panels side by side, using the SAME
honest out-of-fold CV predictions (winning GLM on all 281) in both -- only the
highlighted (removed) points differ:

  LEFT  : the 15 M-estimator-removed bursts in red   -> they lie OFF the 1:1 line
          (badly predicted), so removing them from the test set lifts r.
  RIGHT : the 14 luminosity-cut (faintest 5%) bursts in red -> they lie roughly ON
          the 1:1 line (well predicted, just faint), so removing them barely moves r.

Each panel annotates r on all points and r after dropping its red points.
Saves 'true_vs_pred_sidebyside.png'.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols, r_lin

# ---- honest out-of-fold CV predictions for all 281 (winning GLM) -------------
y, X = design(WINNER, d); n = len(y)
def oof(reps=25, K=10, seed=0):
    P = np.zeros(n)
    for r in range(reps):
        rng = np.random.default_rng(seed + r); fold = rng.permutation(np.arange(n) % K)
        Pr = np.full(n, np.nan)
        for k in range(K):
            te = fold == k; Pr[te] = _ols(X[~te], y[~te], X[te])
        P += Pr
    return P / reps
pred_log = oof()
z_true, z_pred = 10 ** y - 1, 10 ** pred_log - 1

# ---- the two removal sets ---------------------------------------------------
mest = d.index.isin(d.index.difference(d_cut.index))                 # 15 M-est
# luminosity ranking (same cosmology/k-correction as the L-z plot)
C_KMS, H0, OM, OL = 2.99792458e5, 70.0, 0.30, 0.70; MPC_CM = 3.0856775814913673e24
D_H = (C_KMS / H0) * MPC_CM
zz = d["Redshift_crosscheck"].to_numpy(); beta = d["Beta"].to_numpy(); bref = np.median(beta)
zg = np.linspace(0, zz.max() + 1e-6, 20001); Ei = 1/np.sqrt(OM*(1+zg)**3 + OL)
Dc = D_H * np.concatenate(([0], np.cumsum((Ei[1:]+Ei[:-1])*np.diff(zg)/2))); DL = (1+zz)*np.interp(zz, zg, Dc)
logL = np.log10(4*np.pi*DL**2 * 10**d["log10Fa"].to_numpy() * (1+zz)**(beta-1))
rank = logL - (2*np.log10(DL) + (bref-1)*np.log10(1+zz))
lum = np.zeros(n, bool); lum[np.argsort(rank)[:int(round(0.05*n))]] = True   # 14 faintest

# ---- plot -------------------------------------------------------------------
hi = max(z_true.max(), z_pred.max()) * 1.05
def panel(ax, red, title):
    ax.plot([0, hi], [0, hi], "k--", lw=1, alpha=0.6)
    ax.scatter(z_true[~red], z_pred[~red], s=16, c="#9ec3e6", alpha=0.8, edgecolor="none", label="kept")
    ax.scatter(z_true[red],  z_pred[red],  s=60, c="#d7301f", edgecolor="k", lw=0.4, zorder=5,
               label=f"removed ({red.sum()})")
    ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal")
    ax.set_xlabel("observed z"); ax.set_ylabel("predicted z"); ax.set_title(title, fontsize=11)
    ax.text(0.04, 0.96, f"r (all) = {r_lin(pred_log, y):.3f}\n"
                        f"r (after removing red) = {r_lin(pred_log[~red], y[~red]):.3f}",
            transform=ax.transAxes, va="top", fontsize=9.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.legend(loc="lower right", fontsize=8.5)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 6.2))
panel(ax[0], mest, "M-estimator cut (residual-based)\nremoves the worst-predicted points  ->  r: 0.53 -> 0.62")
panel(ax[1], lum,  "Luminosity cut (faintest 5%, physical)\nremoves faint, not worst-predicted  ->  r: 0.53 -> 0.54")
fig.suptitle("True vs predicted redshift (winning GLM, honest out-of-fold CV) — same predictions, different 5% removed",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("true_vs_pred_sidebyside.png", dpi=150)
print("wrote true_vs_pred_sidebyside.png")
print(f"r all = {r_lin(pred_log, y):.3f}")
print(f"r after M-est removal   = {r_lin(pred_log[~mest], y[~mest]):.3f}  (+{r_lin(pred_log[~mest], y[~mest])-r_lin(pred_log,y):.3f})")
print(f"r after luminosity cut  = {r_lin(pred_log[~lum],  y[~lum]):.3f}  (+{r_lin(pred_log[~lum], y[~lum])-r_lin(pred_log,y):.3f})")
print(f"overlap of the two removal sets: {(mest & lum).sum()} bursts")
