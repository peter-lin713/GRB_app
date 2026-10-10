#!/usr/bin/env python3
"""
Performance on the 10 residual-only outliers (the 15 M-est-removed minus the 5
verified different-class bursts), by name.

LEFT : true vs predicted redshift (honest out-of-fold CV, winning GLM). The other
       271 bursts in red, the 10 outliers in blue and labelled by GRB name.
RIGHT: density of absolute redshift error |z_pred - z_true|, others (red) vs the
       10 outliers (blue), with each group's mean absolute error (MAE) marked.

Saves 'outlier10_performance.png'.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols

# ---- honest out-of-fold CV predictions (winning GLM, all 281) ----------------
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
ae = np.abs(z_pred - z_true)                                   # absolute error in redshift

# ---- the 10 residual-only outliers ------------------------------------------
CONTAM = ["070809A", "070714B", "060614A", "171205A", "100316D"]
out10  = d.index.difference(d_cut.index).difference(CONTAM)
blue   = d.index.isin(out10)                                  # the 10
red    = ~blue                                                # the other 271
print(f"10 difficult-to-predict GRBs: {sorted(out10)}")
print(f"MAE  others = {ae[red].mean():.3f}   difficult = {ae[blue].mean():.3f}   ({ae[blue].mean()/ae[red].mean():.1f}x)")

# ---- figure -----------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(14, 6.3))

# LEFT: true vs predicted, names on the 10
hi = max(z_true.max(), z_pred.max()) * 1.05
ax[0].plot([0, hi], [0, hi], "k--", lw=1, alpha=.6, label="1:1")
ax[0].scatter(z_true[red],  z_pred[red],  s=16, c="#d7301f", alpha=.55, edgecolor="none",
              label=f"others ({red.sum()})  MAE={ae[red].mean():.2f}")
ax[0].scatter(z_true[blue], z_pred[blue], s=95, c="#2166ac", edgecolor="k", lw=.5, zorder=5,
              label=f"difficult to predict (10)  MAE={ae[blue].mean():.2f}")

# number the 10 by observed z, leader-line the numbers, and list the names in a key
diff_list = list(out10)
idx = [list(d.index).index(g) for g in diff_list]
gx = np.array([z_true[i] for i in idx]); gy = np.array([z_pred[i] for i in idx], float)
numlabel = {k: rank + 1 for rank, k in enumerate(np.argsort(gx))}          # 1..10 by observed z
lx, ly = gx.copy(), gy.copy()
for _ in range(500):                                                       # repel number positions
    moved = False
    for a in range(len(lx)):
        for b in range(a + 1, len(lx)):
            if abs(lx[a] - lx[b]) < 0.55 and abs(ly[a] - ly[b]) < 0.45:
                sh = (0.45 - abs(ly[a] - ly[b])) / 2 + 0.01
                if ly[a] >= ly[b]: ly[a] += sh; ly[b] -= sh
                else: ly[a] -= sh; ly[b] += sh
                moved = True
    if not moved: break
for k in range(len(diff_list)):
    ax[0].annotate(str(numlabel[k]), xy=(gx[k], gy[k]), xytext=(lx[k], ly[k]), zorder=7,
                   fontsize=8.5, fontweight="bold", color="#08306b", ha="center", va="center",
                   bbox=dict(boxstyle="circle,pad=0.12", fc="white", ec="#2166ac", lw=0.9),
                   arrowprops=dict(arrowstyle="-", color="#6baed6", lw=0.7))
inv = {numlabel[k]: diff_list[k] for k in range(len(diff_list))}
key = "difficult to predict:\n" + "\n".join(f"{nn:2d}  {inv[nn]}" for nn in range(1, len(diff_list) + 1))
ax[0].text(0.985, 0.97, key, transform=ax[0].transAxes, fontsize=8.2, family="monospace",
           ha="right", va="top", bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.95))
ax[0].set_xlim(0, hi); ax[0].set_ylim(0, hi); ax[0].set_aspect("equal")
ax[0].set_xlabel("observed z"); ax[0].set_ylabel("predicted z")
ax[0].set_title("True vs predicted redshift (winning GLM, honest CV)\n10 difficult-to-predict GRBs numbered")
ax[0].legend(loc="upper left", fontsize=9)

# RIGHT: density of absolute error
sns.kdeplot(ae[red],  ax=ax[1], color="#d7301f", fill=True, alpha=.35, lw=2, label=f"others ({red.sum()})")
sns.kdeplot(ae[blue], ax=ax[1], color="#2166ac", fill=True, alpha=.35, lw=2, label="difficult to predict (10)")
ax[1].plot(ae[blue], np.full(blue.sum(), -0.02), "|", color="#2166ac", markersize=12)   # rug for the 10
ax[1].axvline(ae[red].mean(),  color="#d7301f", ls="--", lw=1.5, label=f"MAE others = {ae[red].mean():.2f}")
ax[1].axvline(ae[blue].mean(), color="#2166ac", ls="--", lw=1.5, label=f"MAE difficult = {ae[blue].mean():.2f}")
ax[1].set_xlabel("absolute redshift error  |z$_{pred}$ - z$_{true}$|")
ax[1].set_ylabel("density")
ax[1].set_title("Density of absolute error: difficult points vs the rest")
ax[1].legend(fontsize=9)
fig.tight_layout()
fig.savefig("outlier10_performance.png", dpi=150)
print("wrote outlier10_performance.png")
