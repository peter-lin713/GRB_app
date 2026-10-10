#!/usr/bin/env python3
"""
2x2 training-vs-evaluation decomposition, but cutting the 5% by LUMINOSITY
(the faintest points relative to the limiting-luminosity line we discussed) rather
than by the M-estimator's residual. This cut is model-INDEPENDENT (it uses flux,
redshift and the spectral k-correction, not the GLM's fit), so it's the honest,
physically-motivated version of "drop the hardest 5%".

Question: does removing these faint bursts inflate the winning GLM's r through
EVALUATION (they're absent from the test set) the same way the M-est cut did, or
is a physical cut more benign?

Same design as outlier_train_vs_eval.py: winning GLM, many 80/20 splits,
train {with/without the cut points} x test {with/without}. r in linear z.
Saves 'luminosity_cut_train_vs_eval.png'.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols, r_lin

# ---- luminosity ranking (same cosmology/k-correction as the L-z plot) --------
C_KMS, H0, OM, OL = 2.99792458e5, 70.0, 0.30, 0.70
MPC_CM = 3.0856775814913673e24
D_H = (C_KMS / H0) * MPC_CM
def lum_dist_cm(z):
    zg = np.linspace(0.0, float(np.max(z)) + 1e-6, 20001)
    Ei = 1.0 / np.sqrt(OM * (1 + zg) ** 3 + OL)
    Dc = D_H * np.concatenate(([0.0], np.cumsum((Ei[1:] + Ei[:-1]) * np.diff(zg) / 2)))
    return (1 + np.asarray(z)) * np.interp(z, zg, Dc)

Fa   = 10.0 ** d["log10Fa"].to_numpy()
zz   = d["Redshift_crosscheck"].to_numpy()
beta = d["Beta"].to_numpy(); bref = float(np.median(beta))
DL   = lum_dist_cm(zz)
logL = np.log10(4 * np.pi * DL ** 2 * Fa * (1 + zz) ** (beta - 1))
rank = logL - (2 * np.log10(DL) + (bref - 1) * np.log10(1 + zz))   # height above limiting line

N = int(round(0.05 * len(d)))                                      # 5% -> 14 bursts
lum_out = d.index[np.argsort(rank)[:N]]
is_out  = d.index.isin(lum_out)
mest    = d.index.difference(d_cut.index)
print(f"luminosity cut: dropped {N} faintest bursts ({100*N/len(d):.1f}%)")
print(f"  overlap with the 15 M-est-removed: {len(set(lum_out) & set(mest))} of {N}")
print(f"  bursts: {sorted(lum_out)}\n")

# ---- 2x2 ---------------------------------------------------------------------
y, X = design(WINNER, d); n = len(y)
def two_by_two(n_splits=300, seed=0):
    rng = np.random.default_rng(seed); acc = {k: [] for k in ("ww","wo","cw","co")}
    for _ in range(n_splits):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.2*n), replace=False)] = True
        tr = ~te; tr_clean = tr & ~is_out
        te_clean = te & ~is_out
        b_with  = np.linalg.lstsq(X[tr],       y[tr],       rcond=None)[0]
        b_clean = np.linalg.lstsq(X[tr_clean], y[tr_clean], rcond=None)[0]
        if te_clean.sum() > 2:
            acc["ww"].append(r_lin(X[te]       @ b_with,  y[te]))
            acc["wo"].append(r_lin(X[te_clean] @ b_with,  y[te_clean]))
            acc["cw"].append(r_lin(X[te]       @ b_clean, y[te]))
            acc["co"].append(r_lin(X[te_clean] @ b_clean, y[te_clean]))
    return {k: float(np.mean(v)) for k, v in acc.items()}
r = two_by_two()

print("Winning GLM, 2x2 (r in linear z, avg 300 splits), LUMINOSITY 5% cut:\n")
print("                          test WITH cut pts    test WITHOUT cut pts")
print(f"  train WITH cut pts          {r['ww']:.3f}                {r['wo']:.3f}")
print(f"  train WITHOUT cut pts       {r['cw']:.3f}                {r['co']:.3f}\n")
ev = ((r['wo']-r['ww'])+(r['co']-r['cw']))/2
tr = ((r['cw']-r['ww'])+(r['co']-r['wo']))/2
print(f"  EVALUATION effect: {ev:+.3f}")
print(f"  TRAINING   effect: {tr:+.3f}")
print(f"  (M-est 5% cut, for comparison: evaluation +0.088, training +0.021)")

fig, ax = plt.subplots(figsize=(6.2, 5.3))
xp, w = np.arange(2), 0.36
ax.bar(xp - w/2, [r["ww"], r["cw"]], w, color="#d7301f", label="test WITH cut pts")
ax.bar(xp + w/2, [r["wo"], r["co"]], w, color="#2c7fb8", label="test WITHOUT cut pts")
for x0, v in zip(xp - w/2, [r["ww"], r["cw"]]): ax.text(x0, v+.006, f"{v:.3f}", ha="center", fontsize=9)
for x0, v in zip(xp + w/2, [r["wo"], r["co"]]): ax.text(x0, v+.006, f"{v:.3f}", ha="center", fontsize=9)
ax.set_xticks(xp); ax.set_xticklabels(["train WITH", "train WITHOUT"])
ax.set_ylabel("r (linear z), avg 300 splits"); ax.set_ylim(0.45, 0.67)
ax.set_title("Physical LUMINOSITY 5% cut: training vs evaluation effect\n(winning GLM)")
ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.25)
fig.tight_layout(); fig.savefig("luminosity_cut_train_vs_eval.png", dpi=150)
print("\nwrote luminosity_cut_train_vs_eval.png")
