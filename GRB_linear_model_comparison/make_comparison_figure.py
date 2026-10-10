#!/usr/bin/env python3
"""
Side-by-side predicted-vs-observed redshift figures:
  columns : winning GLM  vs  plain 7-feature additive linear model
  top row : the FLAGSHIP protocol  (M-est cut + repeated 10-fold CV + 2-sigma clip)
  bottom  : an HONEST cross-validation (full sample, no cut, no clip)

Lets you see that under the flagship protocol both models look similar and good,
while honestly both are markedly weaker -- i.e. the lift is the protocol, not the
model. Saves 'comparison_linear_vs_glm.png' next to this script.

Reuses the data, the M-est cut and the model definitions from
simple_linear_flagship_check.py (imported, so that file's __main__ does not run).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, SIMPLE7, design, _ols, r_lin


def cv_oof(formula, frame, seed=0, reps=10, K=10):
    """Out-of-fold repeated-CV predictions (log10(1+z)) for every row of `frame`."""
    y, X = design(formula, frame)
    n, rng = len(y), np.random.default_rng(seed)
    P = np.full((n, reps), np.nan)
    for r in range(reps):
        fold = rng.permutation(np.arange(n) % K)
        for k in range(K):
            te = fold == k
            P[te, r] = _ols(X[~te], y[~te], X[te])
    return P.mean(axis=1), y


def panel(ax, pred_log, obs_log, title, clip=False):
    zp, zo = 10 ** pred_log - 1, 10 ** obs_log - 1
    if clip:                                              # flagship: show the 2-sigma clip
        L = np.std(obs_log - pred_log, ddof=1)
        inside = np.abs(pred_log - obs_log) < 2 * L
        ax.scatter(zo[inside],  zp[inside],  s=14, c="#2c7fb8", alpha=0.75, label="kept (inside 2σ)")
        ax.scatter(zo[~inside], zp[~inside], s=30, c="#d7301f", marker="x", label="clipped (>2σ)")
        r_all  = r_lin(pred_log, obs_log)
        r_clip = r_lin(pred_log[inside], obs_log[inside])
        txt = f"r = {r_all:.3f}  (all points)\nr = {r_clip:.3f}  (after 2σ clip)"
    else:                                                 # honest: no clip
        ax.scatter(zo, zp, s=14, c="#2c7fb8", alpha=0.75)
        txt = f"r = {r_lin(pred_log, obs_log):.3f}"
    hi = max(zo.max(), zp.max()) * 1.05
    ax.plot([0, hi], [0, hi], "k--", lw=1, alpha=0.6)     # 1:1 line
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    ax.set_title(title, fontsize=10)
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.set_xlabel("observed z"); ax.set_ylabel("predicted z")
    if clip:
        ax.legend(loc="lower right", fontsize=7, framealpha=0.9)


fig, ax = plt.subplots(2, 2, figsize=(11, 10))
# top: flagship protocol on the M-est-cut set
gp, go = cv_oof(WINNER,  d_cut); panel(ax[0, 0], gp, go, "Winning GLM  —  flagship protocol\n(M-est cut + CV + 2σ clip)", clip=True)
sp, so = cv_oof(SIMPLE7, d_cut); panel(ax[0, 1], sp, so, "Simple linear (7 feat)  —  flagship protocol\n(M-est cut + CV + 2σ clip)", clip=True)
# bottom: honest CV on the full sample
gp2, go2 = cv_oof(WINNER,  d);   panel(ax[1, 0], gp2, go2, "Winning GLM  —  honest CV\n(full sample, no cut, no clip)")
sp2, so2 = cv_oof(SIMPLE7, d);   panel(ax[1, 1], sp2, so2, "Simple linear (7 feat)  —  honest CV\n(full sample, no cut, no clip)")

fig.suptitle("Predicted vs observed redshift — winning GLM vs plain linear model", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig("comparison_linear_vs_glm.png", dpi=150)
print("wrote comparison_linear_vs_glm.png")
