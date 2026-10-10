#!/usr/bin/env python3
"""
Re-run the training-vs-evaluation 2x2 AFTER removing the 5 bursts we verified are
a different physical class (they don't belong in a long-GRB plateau sample at all):

    070809A, 070714B  -> short / merger-class GRBs
    060614A           -> the no-supernova burst that forced a new GRB class
    171205A, 100316D  -> low-luminosity SN-associated GRBs

Those 5 are dropped from the sample entirely (physically justified). The remaining
10 M-estimator-removed bursts have NO such justification -- they were removed only
because the winning formula predicts their redshift poorly. Question: do those 10
still hurt mainly through EVALUATION (hard to score) rather than TRAINING (poison
the fit)? If yes, the circularity concern stands for them specifically.

Same design as outlier_train_vs_eval.py, now on the 276-burst sample with the 10
as the "outliers".  Saves 'outlier_train_vs_eval_drop5.png'.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols, r_lin

CONTAM = ["070809A", "070714B", "060614A", "171205A", "100316D"]     # verified different-class
all15  = d.index.difference(d_cut.index)
out10  = all15.difference(CONTAM)                                     # the 10 with no physical excuse

dd     = d.drop(CONTAM)                                               # 276 bursts: sample minus contaminants
is_out = dd.index.isin(out10)
y, X   = design(WINNER, dd)
n      = len(y)
print(f"dropped {len(CONTAM)} class contaminants -> {n} bursts; "
      f"outliers now = {is_out.sum()} (the 10 residual-only removals)\n")

def two_by_two(n_splits=300, seed=0):
    rng = np.random.default_rng(seed)
    acc = {k: [] for k in ("ww", "wo", "cw", "co")}
    out_res = []
    for _ in range(n_splits):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.20 * n), replace=False)] = True
        tr = ~te
        tr_clean = tr & ~is_out
        te_out, te_clean = te & is_out, te & ~is_out
        b_with  = np.linalg.lstsq(X[tr],       y[tr],       rcond=None)[0]
        b_clean = np.linalg.lstsq(X[tr_clean], y[tr_clean], rcond=None)[0]
        if te_clean.sum() > 2:
            acc["ww"].append(r_lin(X[te] @ b_with,  y[te]))
            acc["wo"].append(r_lin(X[te_clean] @ b_with,  y[te_clean]))
            acc["cw"].append(r_lin(X[te] @ b_clean, y[te]))
            acc["co"].append(r_lin(X[te_clean] @ b_clean, y[te_clean]))
        if te_out.any():
            out_res.append(X[te_out] @ b_clean - y[te_out])
    return {k: float(np.mean(v)) for k, v in acc.items()}, np.concatenate(out_res)

r, _ = two_by_two()
print("Winning GLM, 2x2 (r in linear z, avg 300 splits), 10 residual-only outliers:\n")
print("                         test WITH the 10     test WITHOUT the 10")
print(f"  train WITH the 10          {r['ww']:.3f}                {r['wo']:.3f}")
print(f"  train WITHOUT the 10       {r['cw']:.3f}                {r['co']:.3f}\n")
ev = ((r['wo']-r['ww'])+(r['co']-r['cw']))/2
tr = ((r['cw']-r['ww'])+(r['co']-r['wo']))/2
print(f"  EVALUATION effect: {ev:+.3f}")
print(f"  TRAINING   effect: {tr:+.3f}")
print( "  (for comparison, with all 15: evaluation +0.088, training +0.021)\n")

# sharp: train on 266 clean, predict the 10 held-out residual-only outliers
b_clean_all = np.linalg.lstsq(X[~is_out], y[~is_out], rcond=None)[0]
pred_out = X[is_out] @ b_clean_all
def cv_pred(mask, reps=10, K=10, seed=1):
    idx = np.where(mask)[0]; rng = np.random.default_rng(seed)
    P = np.full((len(idx), reps), np.nan)
    for r_ in range(reps):
        fold = rng.permutation(np.arange(len(idx)) % K)
        for k in range(K):
            te = fold == k
            P[te, r_] = _ols(X[idx[~te]], y[idx[~te]], X[idx[te]])
    return P.mean(1)
pred_clean = cv_pred(~is_out)
res_out, res_clean = pred_out - y[is_out], pred_clean - y[~is_out]
print("Sharp check -- model trained on the 266 clean bursts:")
print(f"  held-out 10:  RMSE(log)={np.sqrt(np.mean(res_out**2)):.3f}  r(linear z)={r_lin(pred_out, y[is_out]):.3f}")
print(f"  clean CV 266: RMSE(log)={np.sqrt(np.mean(res_clean**2)):.3f}  r(linear z)={r_lin(pred_clean, y[~is_out]):.3f}")
print(f"  RMSE ratio = {np.sqrt(np.mean(res_out**2))/np.sqrt(np.mean(res_clean**2)):.2f}x")

# ---- figure -----------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 5.3))
xpos, w = np.arange(2), 0.36
ax[0].bar(xpos - w/2, [r["ww"], r["cw"]], w, color="#d7301f", label="test WITH the 10")
ax[0].bar(xpos + w/2, [r["wo"], r["co"]], w, color="#2c7fb8", label="test WITHOUT the 10")
for xp, v in zip(xpos - w/2, [r["ww"], r["cw"]]): ax[0].text(xp, v+.006, f"{v:.3f}", ha="center", fontsize=9)
for xp, v in zip(xpos + w/2, [r["wo"], r["co"]]): ax[0].text(xp, v+.006, f"{v:.3f}", ha="center", fontsize=9)
ax[0].set_xticks(xpos); ax[0].set_xticklabels(["train WITH", "train WITHOUT"])
ax[0].set_ylabel("r (linear z), avg 300 splits"); ax[0].set_ylim(0.45, 0.67)
ax[0].set_title("After dropping the 5 class contaminants:\ntraining vs evaluation effect of the remaining 10")
ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(axis="y", alpha=0.25)

zo_out, zp_out = 10**y[is_out]-1,  10**pred_out-1
zo_cln, zp_cln = 10**y[~is_out]-1, 10**pred_clean-1
hi = max(zo_out.max(), zo_cln.max(), zp_out.max(), zp_cln.max())*1.05
ax[1].plot([0, hi], [0, hi], "k--", lw=1, alpha=.6)
ax[1].scatter(zo_cln, zp_cln, s=14, c="#9ec3e6", alpha=.8, label=f"clean, honest CV (266)  r={r_lin(pred_clean,y[~is_out]):.2f}")
ax[1].scatter(zo_out, zp_out, s=60, c="#d7301f", edgecolor="k", lw=.4, zorder=5,
              label=f"10 residual-only outliers  r={r_lin(pred_out,y[is_out]):.2f}")
ax[1].set_xlim(0, hi); ax[1].set_ylim(0, hi)
ax[1].set_xlabel("observed z"); ax[1].set_ylabel("predicted z")
ax[1].set_title("Trained on 266 clean bursts,\npredicting the 10 held-out residual-only outliers")
ax[1].legend(fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig("outlier_train_vs_eval_drop5.png", dpi=150)
print("\nwrote outlier_train_vs_eval_drop5.png")
