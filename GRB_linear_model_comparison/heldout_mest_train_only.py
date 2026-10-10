#!/usr/bin/env python3
"""
Honest held-out protocol (the correct 'M-estimator on training only'):

  1. set aside 20 bursts as an untouched held-out set;
  2. on the remaining ~261 TRAINING bursts, fit the winning GLM as an M-estimator
     (Huber) and drop the lowest-5%-weight training bursts (the catastrophic
     outliers) -- TRAINING ONLY;
  3. retrain the winning GLM (OLS) on the cleaned training set;
  4. predict the 20 held-out bursts and report r -- the held-out set is scored in
     full, nothing removed from it.

Because r on 20 points is noisy, we do one representative split AND average over
many random 20-point holdouts. We also report the plain-OLS-training baseline (no
M-est) and the flagship number (remove outliers from the FULL set) for contrast.
r is Pearson in linear z = 10**pred - 1.
"""
import numpy as np
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols, r_lin

y, X = design(WINNER, d); n = len(y)

def one_split(rng):
    ho = rng.choice(n, 20, replace=False)
    te = np.zeros(n, bool); te[ho] = True; tr = np.where(~te)[0]
    # M-estimator on TRAINING only -> drop lowest 5% by weight
    w = sm.RLM(y[tr], X[tr], M=sm.robust.norms.HuberT()).fit().weights
    keep = tr[w > np.quantile(w, 0.05)]
    pred_mest = _ols(X[keep], y[keep], X[te])            # retrain GLM on cleaned training
    pred_ols  = _ols(X[tr],   y[tr],   X[te])            # baseline: no M-est
    return r_lin(pred_mest, y[te]), r_lin(pred_ols, y[te]), ho

# ---- one representative split -----------------------------------------------
rng = np.random.default_rng(7)
rm, ro, ho = one_split(rng)
print("ONE representative 20-point held-out split (seed 7):")
print(f"  M-est-on-training + retrained GLM:  held-out r = {rm:.3f}")
print(f"  plain OLS training (no M-est):      held-out r = {ro:.3f}")
print(f"  (held-out GRBs: {sorted(d.index[ho])})\n")

# ---- averaged over many random 20-point holdouts ----------------------------
rng = np.random.default_rng(0); M, O = [], []
NREP = 1000
for _ in range(NREP):
    a, b, _ = one_split(rng); M.append(a); O.append(b)
M, O = np.array(M), np.array(O)
print(f"Averaged over {NREP} random 20-point holdouts:")
print(f"  M-est-on-training + retrained GLM:  held-out r = {M.mean():.3f} ± {M.std():.3f}")
print(f"  plain OLS training (no M-est):      held-out r = {O.mean():.3f} ± {O.std():.3f}")
print(f"  M-est training-clean helps by:      {M.mean()-O.mean():+.3f}\n")

# ---- flagship contrast: remove outliers from the FULL set, then CV -----------
y2, X2 = design(WINNER, d[~d.index.isin(d.index.difference(d_cut.index))])
rs = []
for r in range(20):
    rng2 = np.random.default_rng(r); fold = rng2.permutation(np.arange(len(y2)) % 10); P = np.full(len(y2), np.nan)
    for k in range(10):
        te = fold == k; P[te] = _ols(X2[~te], y2[~te], X2[te])
    rs.append(r_lin(P, y2))
flagship = float(np.mean(rs))
print(f"CONTRAST -- flagship (remove 15 from FULL set, CV on 266): r = {flagship:.3f}")
print("  (the held-out numbers above are the honest ones; the flagship number")
print("   removes the hard points from the evaluation itself.)")

# ---- distribution figure ----------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.6))
bins = np.linspace(-0.1, 1.0, 45)
ax.hist(M, bins=bins, alpha=0.8, color="#2c7fb8",
        label=f"M-est on training + retrain  (mean {M.mean():.2f} ± {M.std():.2f})")
_hc, _edg = np.histogram(M, bins=bins)
mode_val = (_edg[np.argmax(_hc)] + _edg[np.argmax(_hc) + 1]) / 2          # tallest-bin centre
ax.axvline(M.mean(),      color="#08519c", lw=2, label=f"mean   = {M.mean():.2f}")
ax.axvline(np.median(M),  color="#08519c", lw=1.6, ls=":", label=f"median = {np.median(M):.2f}")
ax.axvline(mode_val,      color="#08519c", lw=1.6, ls="-.", label=f"mode   = {mode_val:.2f}")
ax.axvline(flagship, color="#d7301f", lw=2, ls="--",
           label=f"flagship (remove from full set) = {flagship:.2f}")
ax.set_xlabel("held-out r (20-point held-out set, linear z)")
ax.set_ylabel(f"count over {NREP} random holdouts")
ax.set_title("Distribution of honest held-out r  —  M-estimator on TRAINING only\n"
             f"single 20-point draws span ~0.1–0.9 (sd {M.std():.2f})")
ax.legend(fontsize=8.5, loc="upper left")
ax.grid(axis="y", alpha=0.25)
fig.tight_layout(); fig.savefig("heldout_r_distribution.png", dpi=150)
print("\nwrote heldout_r_distribution.png")
