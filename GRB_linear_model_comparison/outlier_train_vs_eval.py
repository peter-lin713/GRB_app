#!/usr/bin/env python3
"""
Do the 15 M-estimator outliers hurt mainly through TRAINING or through EVALUATION?

The flagship removes them from BOTH the training fit and the scored test set, so
its r gain conflates two distinct effects:
  (T) training effect  -- outliers distort the fitted coefficients, so removing
      them from TRAINING gives a better model for everyone else;
  (E) evaluation effect -- outliers are intrinsically hard to predict, so removing
      them from the TEST set raises r regardless of the model.

We separate them with the winning GLM over many random 80/20 splits, filling a 2x2:

                         test = 20% WITH outliers   test = 20% WITHOUT outliers
    train 80% WITH        r_ww                        r_wo
    train 80% WITHOUT     r_cw                        r_co

  - Evaluation effect  = (fix training) r_*o - r_*w   > 0 if outliers are hard to score
  - Training effect    = (fix test)     r_c* - r_w*   > 0 if outliers poison the fit

Sharp version (bottom): train on ALL 266 non-outliers, then predict the 15
outliers -- which were NEVER in training -- and compare their errors to honest CV
errors on the 266. If the clean-trained model still predicts the outliers badly,
the outliers are an evaluation problem, not a training problem.

All r are Pearson r in linear z = 10**pred - 1.  Saves 'outlier_train_vs_eval.png'.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut, WINNER, design, _ols, r_lin

outliers = d.index.difference(d_cut.index)                 # the 15 M-est-removed
is_out   = d.index.isin(outliers)
y, X     = design(WINNER, d)                               # full design once
n        = len(y)
oidx     = np.where(is_out)[0]

# ---- 2x2: train {with,without} x test {with,without}, averaged over splits ----
def two_by_two(n_splits=300, seed=0):
    rng = np.random.default_rng(seed)
    acc = {k: [] for k in ("ww", "wo", "cw", "co")}
    out_res = []                                            # residuals of outliers that land in test (clean-trained)
    for _ in range(n_splits):
        te = np.zeros(n, bool); te[rng.choice(n, int(0.20 * n), replace=False)] = True
        tr = ~te
        tr_clean = tr & ~is_out                             # drop outliers from TRAIN only
        te_out   = te & is_out
        te_clean = te & ~is_out
        b_with  = np.linalg.lstsq(X[tr],       y[tr],       rcond=None)[0]
        b_clean = np.linalg.lstsq(X[tr_clean], y[tr_clean], rcond=None)[0]
        if te_clean.sum() > 2:
            acc["ww"].append(r_lin(X[te] @ b_with,  y[te]))
            acc["wo"].append(r_lin(X[te_clean] @ b_with,  y[te_clean]))
            acc["cw"].append(r_lin(X[te] @ b_clean, y[te]))
            acc["co"].append(r_lin(X[te_clean] @ b_clean, y[te_clean]))
        if te_out.any():                                   # clean-trained residuals on held-out outliers
            out_res.append(X[te_out] @ b_clean - y[te_out])
    r = {k: float(np.mean(v)) for k, v in acc.items()}
    return r, np.concatenate(out_res)

r, _ = two_by_two()
print("Winning GLM, 2x2 (r in linear z, avg over 300 random 80/20 splits):\n")
print("                         test WITH outliers   test WITHOUT outliers")
print(f"  train WITH outliers        {r['ww']:.3f}                {r['wo']:.3f}")
print(f"  train WITHOUT outliers     {r['cw']:.3f}                {r['co']:.3f}\n")
print(f"  EVALUATION effect (test without - with, avg over train rows): "
      f"{((r['wo']-r['ww'])+(r['co']-r['cw']))/2:+.3f}")
print(f"  TRAINING   effect (train without - with, avg over test cols): "
      f"{((r['cw']-r['ww'])+(r['co']-r['wo']))/2:+.3f}\n")

# ---- sharp version: train on all 266 clean, predict the 15 outliers ----------
b_clean_all = np.linalg.lstsq(X[~is_out], y[~is_out], rcond=None)[0]
pred_out    = X[is_out] @ b_clean_all                      # outliers were never in training
# honest CV residuals on the 266 clean (10x repeated 10-fold)
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
res_out   = pred_out   - y[is_out]                          # log10(1+z) residuals
res_clean = pred_clean - y[~is_out]
print("Sharp check -- model trained on the 266 clean bursts:")
print(f"  held-out outliers (15):   RMSE(log) = {np.sqrt(np.mean(res_out**2)):.3f}   "
      f"r(linear z) = {r_lin(pred_out, y[is_out]):.3f}")
print(f"  honest CV on clean (266): RMSE(log) = {np.sqrt(np.mean(res_clean**2)):.3f}   "
      f"r(linear z) = {r_lin(pred_clean, y[~is_out]):.3f}")
print(f"  outlier RMSE / clean RMSE = {np.sqrt(np.mean(res_out**2))/np.sqrt(np.mean(res_clean**2)):.2f}x")

# ---- figure -----------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 5.3))

# panel A: the 2x2 as grouped bars
labels = ["train WITH", "train WITHOUT"]
xpos = np.arange(2); w = 0.36
ax[0].bar(xpos - w/2, [r["ww"], r["cw"]], w, color="#d7301f", label="test WITH outliers")
ax[0].bar(xpos + w/2, [r["wo"], r["co"]], w, color="#2c7fb8", label="test WITHOUT outliers")
for xp, v in zip(xpos - w/2, [r["ww"], r["cw"]]): ax[0].text(xp, v + .006, f"{v:.3f}", ha="center", fontsize=9)
for xp, v in zip(xpos + w/2, [r["wo"], r["co"]]): ax[0].text(xp, v + .006, f"{v:.3f}", ha="center", fontsize=9)
ax[0].set_xticks(xpos); ax[0].set_xticklabels(labels)
ax[0].set_ylabel("r (linear z), avg 300 splits"); ax[0].set_ylim(0.45, 0.67)
ax[0].set_title("Training vs evaluation effect of the 15 outliers\n(winning GLM, 80/20)")
ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(axis="y", alpha=0.25)

# panel B: clean-trained predictions -- outliers (never trained on) vs clean CV
zo_out, zp_out     = 10**y[is_out]-1,  10**pred_out-1
zo_cln, zp_cln     = 10**y[~is_out]-1, 10**pred_clean-1
hi = max(zo_out.max(), zo_cln.max(), zp_out.max(), zp_cln.max())*1.05
ax[1].plot([0, hi], [0, hi], "k--", lw=1, alpha=.6)
ax[1].scatter(zo_cln, zp_cln, s=14, c="#9ec3e6", alpha=.8, label=f"clean, honest CV (266)  r={r_lin(pred_clean,y[~is_out]):.2f}")
ax[1].scatter(zo_out, zp_out, s=55, c="#d7301f", edgecolor="k", lw=.4, zorder=5,
              label=f"15 outliers, predicted from clean fit  r={r_lin(pred_out,y[is_out]):.2f}")
ax[1].set_xlim(0, hi); ax[1].set_ylim(0, hi)
ax[1].set_xlabel("observed z"); ax[1].set_ylabel("predicted z")
ax[1].set_title("Model trained on 266 clean bursts,\npredicting the 15 held-out outliers")
ax[1].legend(fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig("outlier_train_vs_eval.png", dpi=150)
print("\nwrote outlier_train_vs_eval.png")
