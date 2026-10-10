#!/usr/bin/env python3
"""
Generalization distribution check, side by side: winning GLM vs plain linear.

The 299 GRBs in generalization_299_raw.csv have NO measured redshift. We fit each
model on the labeled sample, predict a redshift for all 299 (in log10(1+z), then
transform to z = 10**pred - 1), and compare the DISTRIBUTION of predicted z to the
distribution of the spectroscopic (labeled) sample. A well-behaved model should
reproduce a similar spread; regression-to-the-mean shows up as a narrow spike.

Both the winning GLM and the simple 7-feature linear model use the SAME seven
features (neither uses Beta or Fluence), so no missing-feature imputation is
needed -- the generalization catalog provides all seven directly.

Saves 'generalization_glm_vs_linear.png' next to this script.
"""
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import ks_2samp
from simple_linear_flagship_check import WINNER, SIMPLE7   # the two model formulas

# ---- labeled training data (add squared terms + log10z target) --------------
lab = pd.read_csv("grb_imputed_281.csv", index_col=0)
for b in ["log10Fa", "log10T90", "log10PeakFlux", "log10NH", "PhotonIndex", "log10Ta", "Alpha"]:
    lab[b + "Sqr"] = lab[b] ** 2
if "log10z" not in lab.columns:
    lab["log10z"] = np.log10(lab["Redshift_crosscheck"] + 1)
z_spec = lab["Redshift_crosscheck"].to_numpy()             # spectroscopic reference distribution

# ---- 299 generalization set: map raw columns -> the 7 model features --------
def clean_num(s):
    """Strip contaminated suffixes ('1.81PL', '2.03,') to the leading numeric token."""
    m = re.match(r"^\s*(-?\d+\.?\d*)", str(s))
    return float(m.group(1)) if m else np.nan

g = pd.read_csv("generalization_299_raw.csv", index_col=0)
gen = pd.DataFrame({
    "log10Fa":       g["Fbest"].map(clean_num),            # Fbest / T_abest are already log-scale
    "log10Ta":       g["T_abest"].map(clean_num),
    "log10T90":      g["logT90"].map(clean_num),
    "log10PeakFlux": g["logPeakFlux"].map(clean_num),
    "PhotonIndex":   g["photon_index"].map(clean_num),
    "log10NH":       g["logNH"].map(clean_num),
    "Alpha":         g["Alpha"].map(clean_num),
})
# physically-impossible values -> NaN, then fill the few remaining gaps with the
# training-set median (only a handful survive the string-cleaning as NaN).
gen.loc[gen["log10NH"] < 20, "log10NH"] = np.nan
gen.loc[gen["PhotonIndex"] < 0, "PhotonIndex"] = np.nan
gen.loc[gen["Alpha"] > 3, "Alpha"] = np.nan
for c in gen.columns:
    gen[c] = gen[c].fillna(lab[c].median())

# Out-of-parameter-space filter (Peter's step 5): drop generalization GRBs whose
# features fall outside the training [min, max] range, to avoid extrapolating the
# model (esp. the winning GLM's interaction/squared terms) far beyond its support.
# Applied identically to both models so they are compared on the same subset.
_n0 = len(gen)
_keep = np.ones(len(gen), bool)
for c in ["log10NH", "PhotonIndex", "log10Fa", "log10Ta", "log10PeakFlux", "log10T90", "Alpha"]:
    _keep &= (gen[c] > lab[c].min()) & (gen[c] < lab[c].max())
gen = gen[_keep].copy()
print(f"generalization set: {len(gen)} of {_n0} kept after out-of-range filter")

for b in ["log10Fa", "log10T90", "log10PeakFlux", "log10NH", "PhotonIndex", "log10Ta", "Alpha"]:
    gen[b + "Sqr"] = gen[b] ** 2

# ---- fit each model on labeled data, predict the 299 (log scale -> z) --------
def predict_z(formula):
    fit = smf.ols(formula, data=lab).fit()
    return 10 ** fit.predict(gen).to_numpy() - 1

z_glm = predict_z(WINNER)
z_lin = predict_z(SIMPLE7)

# ---- side-by-side distributions ---------------------------------------------
def panel(ax, z_pred, title):
    zp = z_pred[np.isfinite(z_pred)]
    lo = np.floor(min(zp.min(), z_spec.min()))
    hi = np.ceil(max(zp.max(), z_spec.max()))
    bins = np.arange(lo, hi + 0.5, 0.5)
    ax.hist(z_spec, bins=bins, density=True, color="#2c7fb8", alpha=0.55,
            label=f"spectroscopic (n={len(z_spec)})")
    ax.hist(zp, bins=bins, density=True, color="#d7301f", alpha=0.55,
            label=f"predicted (n={len(zp)})")
    ks = ks_2samp(zp, z_spec)
    ax.set_title(title, fontsize=11)
    ax.text(0.97, 0.97,
            f"KS D = {ks.statistic:.3f}, p = {ks.pvalue:.2g}\n"
            f"median: pred {np.median(zp):.2f}  vs  spec {np.median(z_spec):.2f}\n"
            f"IQR pred [{np.percentile(zp,25):.2f}, {np.percentile(zp,75):.2f}]  "
            f"spec [{np.percentile(z_spec,25):.2f}, {np.percentile(z_spec,75):.2f}]",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.set_xlabel("redshift z"); ax.set_ylabel("density"); ax.legend(fontsize=8, loc="upper left")

fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), sharex=True, sharey=True)
panel(ax[0], z_glm, "Winning GLM  —  299 no-redshift GRBs")
panel(ax[1], z_lin, "Simple linear (7 feat)  —  299 no-redshift GRBs")
fig.suptitle("Predicted-redshift distribution of the no-redshift GRBs vs the spectroscopic sample",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("generalization_glm_vs_linear.png", dpi=150)
print("wrote generalization_glm_vs_linear.png")
print(f"GLM:    median pred z = {np.median(z_glm):.3f}, range [{np.nanmin(z_glm):.2f}, {np.nanmax(z_glm):.2f}]")
print(f"Linear: median pred z = {np.median(z_lin):.3f}, range [{np.nanmin(z_lin):.2f}, {np.nanmax(z_lin):.2f}]")
