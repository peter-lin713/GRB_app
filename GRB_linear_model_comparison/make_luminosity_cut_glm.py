#!/usr/bin/env python3
"""
Physically-motivated luminosity cut vs the M-estimator cut, under the flagship
protocol, for the WINNING GLM.

The flagship removes ~15 bursts with an M-estimator (a statistical down-weighting
of points the winning formula fits poorly) and that removal lifts the reported r.
Question: if we instead remove the same sort of "hard" bursts on PHYSICAL grounds
-- the low-luminosity bursts that fall below a flux/limiting-luminosity line (a
Malmquist cut) -- does the winning GLM's r rise the same way?

Cut: drop the N bursts that sit lowest relative to a limiting-luminosity line
L_lim(z) proportional to D_L(z)^2 (1+z)^(beta_ref-1). Ranking bursts by height
above such a line is independent of the line's absolute height, so this is just
"remove the N faintest bursts in the L-z plane". N=5 by default (~1.8%); set
N_DROP=14 for a 5% cut or 21 for 7.5% -- see the header table printed below.

Then evaluate the winning GLM under the IDENTICAL flagship protocol (repeated
10-fold CV + 2-sigma cone clip) on the surviving bursts, exactly as done for the
M-estimator cut, so the two cuts are compared apples-to-apples.

Saves 'luminosity_cut_vs_mest_glm.png' next to this script.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols, _cv_clip, r_lin, honest_r

N_DROP = 5                                        # bursts removed by the luminosity cut

# ---- luminosity + limiting-line ranking (same cosmology as the L-z plot) -----
C_KMS, H0, OM, OL = 2.99792458e5, 70.0, 0.30, 0.70
MPC_CM = 3.0856775814913673e24
D_H = (C_KMS / H0) * MPC_CM

def lum_dist_cm(z):
    zg = np.linspace(0.0, float(np.max(z)) + 1e-6, 20001)
    Ei = 1.0 / np.sqrt(OM * (1 + zg) ** 3 + OL)
    Dc = D_H * np.concatenate(([0.0], np.cumsum((Ei[1:] + Ei[:-1]) * np.diff(zg) / 2)))
    return (1 + np.asarray(z)) * np.interp(z, zg, Dc)

Fa   = 10.0 ** d["log10Fa"].to_numpy()
z    = d["Redshift_crosscheck"].to_numpy()
beta = d["Beta"].to_numpy()
bref = float(np.median(beta))
logL = np.log10(4 * np.pi * lum_dist_cm(z) ** 2 * Fa * (1 + z) ** (beta - 1))
rank = logL - (2 * np.log10(lum_dist_cm(z)) + (bref - 1) * np.log10(1 + z))   # height above line (+const)

drop_idx = d.index[np.argsort(rank)[:N_DROP]]
d_lum    = d.drop(drop_idx)
print(f"luminosity cut: dropped {N_DROP} lowest-L bursts ({100*N_DROP/len(d):.1f}%): {list(drop_idx)}")
print(f"  overlap with the 15 M-est-removed bursts: "
      f"{len(set(drop_idx) & set(d.index.difference(d_cut.index)))} of {N_DROP}\n")

# ---- winning GLM r under each training cut (flagship protocol) ---------------
def flagship(frame, n_splits=20, reps=5, K=10, seed=0):
    y, X = design(WINNER, frame)
    rng, n, cv, cl = np.random.default_rng(seed), len(y), [], []
    for _ in range(n_splits):
        ho = rng.choice(n, int(0.2 * n), replace=False)
        tr = np.setdiff1d(np.arange(n), ho)
        c, k = _cv_clip(y[tr], X[tr], rng, reps, K); cv.append(c); cl.append(k)
    return float(np.mean(cv)), float(np.mean(cl))

m_cv, m_cl = flagship(d_cut)
l_cv, l_cl = flagship(d_lum)
h_r = honest_r(WINNER)
print("Winning GLM, r in linear z:")
print(f"  M-estimator cut  (drop 15):   flagship CV r = {m_cv:.3f}   CV+2s-clip r = {m_cl:.3f}")
print(f"  luminosity cut   (drop {N_DROP:2d}):   flagship CV r = {l_cv:.3f}   CV+2s-clip r = {l_cl:.3f}")
print(f"  honest held-out  (no cut):    r = {h_r:.3f}\n")

# ---- predicted-vs-observed panels (same style as make_comparison_figure) -----
def cv_oof(frame, seed=0, reps=10, K=10):
    y, X = design(WINNER, frame)
    n, rng = len(y), np.random.default_rng(seed)
    P = np.full((n, reps), np.nan)
    for r in range(reps):
        fold = rng.permutation(np.arange(n) % K)
        for k in range(K):
            te = fold == k
            P[te, r] = _ols(X[~te], y[~te], X[te])
    return P.mean(axis=1), y

def panel(ax, pred_log, obs_log, title):
    zp, zo = 10 ** pred_log - 1, 10 ** obs_log - 1
    L = np.std(obs_log - pred_log, ddof=1)
    inside = np.abs(pred_log - obs_log) < 2 * L
    ax.scatter(zo[inside],  zp[inside],  s=14, c="#2c7fb8", alpha=0.75, label="kept (inside 2σ)")
    ax.scatter(zo[~inside], zp[~inside], s=30, c="#d7301f", marker="x", label="clipped (>2σ)")
    hi = max(zo.max(), zp.max()) * 1.05
    ax.plot([0, hi], [0, hi], "k--", lw=1, alpha=0.6)
    ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    ax.set_title(title, fontsize=10)
    ax.text(0.04, 0.96, f"r = {r_lin(pred_log, obs_log):.3f}  (all points)\n"
                        f"r = {r_lin(pred_log[inside], obs_log[inside]):.3f}  (after 2σ clip)",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.set_xlabel("observed z"); ax.set_ylabel("predicted z")
    ax.legend(loc="lower right", fontsize=7, framealpha=0.9)

fig, ax = plt.subplots(1, 2, figsize=(11, 5.4))
mp, mo = cv_oof(d_cut); panel(ax[0], mp, mo, "Winning GLM — M-estimator cut (drop 15)\nflagship CV + 2σ clip")
lp, lo = cv_oof(d_lum); panel(ax[1], lp, lo, f"Winning GLM — luminosity cut (drop {N_DROP})\nflagship CV + 2σ clip")
fig.suptitle("Winning GLM under the flagship protocol: statistical (M-est) cut vs physical (luminosity) cut",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("luminosity_cut_vs_mest_glm.png", dpi=150)
print("wrote luminosity_cut_vs_mest_glm.png")
