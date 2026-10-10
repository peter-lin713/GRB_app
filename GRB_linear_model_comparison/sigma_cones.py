#!/usr/bin/env python3
"""
1-sigma and 2-sigma cones on the FULL honest sample -- no M-estimator removal.

Winning GLM, honest out-of-fold CV predictions for all 281 bursts. The residual
scatter in log10(1+z) is sigma = sd(pred - obs). The k-sigma "cone" is the band
|pred_log - obs_log| < k*sigma around the 1:1 line (this is exactly the flagship's
cone definition, L_Sigma = sd(Dlogz)).

There are TWO different things you can do with the cones -- keep them separate:

  (A) HONEST use = COVERAGE / uncertainty. Report what fraction of predictions land
      inside +/-1 sigma and +/-2 sigma. For a well-calibrated model these are ~68%
      and ~95%. This summarises the model WITHOUT deleting any point and WITHOUT
      touching r. This is the legitimate "1/2-sigma result".

  (B) LEAKY use = CLIP then recompute r. Drop the points OUTSIDE the cone and
      recompute r on what's left. This raises r -- but it is the SAME move as the
      M-estimator removal: it deletes the worst-predicted TEST points using their
      true redshift. It is the flagship's "second cut". Not less sketchy; identical
      mechanism, just thresholded on residual-sigma instead of M-est weight.

Saves 'sigma_cones.png'.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, WINNER, design, _ols, r_lin

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
P = oof()
res = P - y                                  # residual in log10(1+z)
sig = res.std(ddof=1)
in1, in2 = np.abs(res) < sig, np.abs(res) < 2 * sig

# (A) coverage
print(f"residual sigma (log10(1+z)) = {sig:.3f}\n")
print("(A) HONEST -- coverage (no points deleted, r untouched):")
print(f"    inside 1 sigma: {100*in1.mean():5.1f}%   (Gaussian ~68%)")
print(f"    inside 2 sigma: {100*in2.mean():5.1f}%   (Gaussian ~95%)")
print(f"    r on ALL 281 points = {r_lin(P, y):.3f}\n")
# (B) clip-and-recompute (leaky)
print("(B) LEAKY -- clip to the cone then recompute r (= flagship 'second cut'):")
print(f"    keep inside 2 sigma ({in2.sum()}/{n}, drop {n-in2.sum()}):  r = {r_lin(P[in2], y[in2]):.3f}")
print(f"    keep inside 1 sigma ({in1.sum()}/{n}, drop {n-in1.sum()}):  r = {r_lin(P[in1], y[in1]):.3f}")
print( "    (dropping the worst-predicted test points using their true z -- same")
print( "     leakage as the M-est cut, just thresholded on residual sigma.)")

# ---- figure: cones on the true-vs-predicted plane ---------------------------
z_true, z_pred = 10 ** y - 1, 10 ** P - 1
hi = max(z_true.max(), z_pred.max()) * 1.05
zz = np.linspace(0, hi, 300); base = np.log10(1 + zz)
def cone(kv): return 10 ** (base + kv * sig) - 1
fig, ax = plt.subplots(figsize=(8.4, 7.4))
ax.plot(zz, zz, "k--", lw=1, alpha=.6, label="1:1")
for kv, c, lab in [(2, "#d7301f", "2σ cone"), (1, "#f4a582", "1σ cone")]:
    ax.plot(zz, cone(kv), color=c, lw=1.4); ax.plot(zz, cone(-kv), color=c, lw=1.4, label=lab)
col = np.where(in1, "#2c7fb8", np.where(in2, "#f4a582", "#d7301f"))
ax.scatter(z_true, z_pred, s=20, c=col, alpha=.85, edgecolor="none")
ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal")
ax.set_xlabel("observed z"); ax.set_ylabel("predicted z")
ax.set_title(f"Winning GLM, honest CV (all 281, no M-est removal)\n"
             f"σ={sig:.3f} | inside 1σ {100*in1.mean():.0f}%, 2σ {100*in2.mean():.0f}% | "
             f"r(all)={r_lin(P,y):.2f}, r(≤2σ)={r_lin(P[in2],y[in2]):.2f}")
ax.legend(loc="upper left", fontsize=9)
fig.tight_layout(); fig.savefig("sigma_cones.png", dpi=150)
print("\nwrote sigma_cones.png")
