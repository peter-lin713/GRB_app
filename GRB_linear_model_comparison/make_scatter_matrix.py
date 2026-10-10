#!/usr/bin/env python3
"""
Scatter matrix (pairs plot) of the 281 imputed GRBs, with the 15 GRBs removed by
the M-estimator catastrophic-outlier cut highlighted in their own color.

The M-est cut here is the flagship one: fit the winning GLM formula with a Huber
M-estimator, give every burst an IRLS weight in (0,1], and drop the lowest 5% by
weight. Points are then plotted across all 10 predictors + log10(z+1).

Saves 'scatter_matrix_mest_removed.png' next to this script.
"""
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---- data + winning-formula design (for the M-est cut) ----------------------
d = pd.read_csv("grb_imputed_281.csv", index_col=0)
for b in ["log10Fa", "log10T90", "log10PeakFlux", "log10NH", "PhotonIndex", "log10Ta", "Alpha"]:
    d[b + "Sqr"] = d[b] ** 2
if "log10z" not in d.columns:
    d["log10z"] = np.log10(d["Redshift_crosscheck"] + 1)

WINNER = ("log10z ~ (log10FaSqr + log10T90 + log10Fa + log10PeakFlux)**2 "
          "+ log10NH + PhotonIndex + log10Ta + Alpha "
          "+ log10NHSqr + log10PeakFluxSqr + PhotonIndexSqr + log10TaSqr + log10T90Sqr + AlphaSqr")

# Huber M-estimator weights; drop the lowest 5% (the "catastrophic outliers").
y, X = patsy.dmatrices(WINNER, d, return_type="dataframe")
w = np.asarray(sm.RLM(y, X, M=sm.robust.norms.HuberT()).fit().weights)
thr = np.quantile(w, 0.05)
removed = d.index[w <= thr]
d["M-estimator"] = np.where(w <= thr, f"removed ({len(removed)})", f"kept ({len(d) - len(removed)})")
print(f"M-est removed {len(removed)} GRBs:\n{list(removed)}")

# ---- scatter matrix over the 10 predictors + log10(z+1) ---------------------
feats = ["log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
         "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux", "log10z"]
kept_lbl, rem_lbl = f"kept ({len(d) - len(removed)})", f"removed ({len(removed)})"

g = sns.pairplot(
    d, vars=feats, hue="M-estimator",
    hue_order=[kept_lbl, rem_lbl],                        # plot removed last -> on top
    palette={kept_lbl: "#9ec3e6", rem_lbl: "#d7301f"},
    plot_kws=dict(s=16, alpha=0.7, edgecolor="none"),
    diag_kind="hist", diag_kws=dict(alpha=0.6),
    corner=True, height=1.25,
)
g.fig.suptitle("Scatter matrix of the 281 imputed GRBs — the 15 M-estimator-removed bursts in red",
               y=1.005, fontsize=14)
g.savefig("scatter_matrix_mest_removed.png", dpi=130, bbox_inches="tight")
print("wrote scatter_matrix_mest_removed.png")
