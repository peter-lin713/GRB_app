#!/usr/bin/env python3
"""
Is each of the 15 M-estimator outliers actually anomalous in FEATURE space, or
does it have ordinary features and get removed only because the winning FORMULA
fits its redshift badly (a residual outlier, not a feature outlier)?

For each outlier separately (never averaged -- they sit in different corners) we
measure the squared Mahalanobis distance to the distribution of the 266 "used"
(non-outlier) bursts:

    D^2 = (x - mu)^T  S^-1  (x - mu)

mu, S = mean and covariance of the 266 kept bursts over the 10 model features.
Mahalanobis is the right normalization for the request:
  - S^-1 divides out feature scales and CORRELATIONS, so a clump of correlated
    features is not counted several times -- we do not over-penalize the joint;
  - to not over-penalize the DIMENSION count, D^2 is compared to chi^2_df (whose
    mean is df) for a p-value, and also reported as D^2/df (≈1 for a typical
    burst). We add an assumption-free empirical p (fraction of the 266 kept whose
    own D^2 is >= the outlier's), which needs no normality.

Small p  => genuinely feature-space anomalous (leverage / extrapolation) -> removal
            is at least partly justified on physical/feature grounds.
Large p  => ordinary features; removed only because the formula mis-predicts its z
            -> a residual/circularity removal, not a feature outlier.

Saves 'outlier_mahalanobis.png'.
"""
import numpy as np
import pandas as pd
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simple_linear_flagship_check import d, d_cut

FEATURES = ["log10T90", "log10Fa", "log10Ta", "Alpha", "Beta", "Gamma",
            "log10Fluence", "PhotonIndex", "log10NH", "log10PeakFlux"]
outliers = d.index.difference(d_cut.index)
kept     = d.index.difference(outliers)

def mahalanobis_report(cols, tag):
    df  = len(cols)
    Xk  = d.loc[kept, cols].to_numpy()
    mu  = Xk.mean(0)
    S   = np.cov(Xk, rowvar=False)
    Si  = np.linalg.inv(S)
    def d2(x): return float((x - mu) @ Si @ (x - mu))
    kept_d2 = np.array([d2(r) for r in Xk])                     # in-sample reference
    rows = []
    for gid in outliers:
        D2 = d2(d.loc[gid, cols].to_numpy())
        rows.append((gid, D2, D2/df, chi2.sf(D2, df),
                     float(np.mean(kept_d2 >= D2))))            # empirical one-sided p
    out = (pd.DataFrame(rows, columns=["GRB", "D2", "D2_over_df", "p_chi2", "p_empirical"])
             .set_index("GRB").sort_values("D2", ascending=False))
    print(f"\n=== {tag}  (df={df}, reference = 266 kept bursts) ===")
    print(f"  kept D^2: median {np.median(kept_d2):.1f}, 95th pct {np.percentile(kept_d2,95):.1f}  "
          f"(chi^2_{df} mean={df}, 95th={chi2.ppf(.95,df):.1f})")
    with pd.option_context("display.float_format", lambda v: f"{v:.3g}"):
        print(out)
    n_sig = (out.p_chi2 < 0.05).sum()
    n_bonf = (out.p_chi2 < 0.05/len(outliers)).sum()
    print(f"  feature-space anomalous: {n_sig}/15 at p<0.05, {n_bonf}/15 at Bonferroni p<{0.05/15:.4f}")
    print(f"  ordinary features (p_chi2 > 0.05): {(out.p_chi2 > 0.05).sum()}/15")
    return out, kept_d2, df

feat_out,  kept_d2,  df  = mahalanobis_report(FEATURES, "Feature space (10 model predictors)")
featz_out, kept_d2z, dfz = mahalanobis_report(FEATURES + ["log10z"], "Feature space + log10z (11-D)")

# ---- figure: per-outlier -log10 p, features-only vs +redshift ----------------
order = feat_out.index                                          # sort by feature-space D^2
xf = -np.log10(feat_out.loc[order, "p_chi2"].clip(lower=1e-300))
xz = -np.log10(featz_out.loc[order, "p_chi2"].clip(lower=1e-300))
yy = np.arange(len(order))[::-1]

fig, ax = plt.subplots(figsize=(9.5, 7))
ax.barh(yy + 0.2, xf, height=0.4, color="#2c7fb8", label="features only (10-D)")
ax.barh(yy - 0.2, xz, height=0.4, color="#d7301f", label="features + redshift (11-D)")
ax.axvline(-np.log10(0.05), color="k", ls="--", lw=1, alpha=.7, label="p = 0.05")
ax.axvline(-np.log10(0.05/15), color="k", ls=":", lw=1, alpha=.7, label="Bonferroni p = 0.0033")
ax.set_yticks(yy); ax.set_yticklabels(order, fontsize=9)
ax.set_xlabel(r"$-\log_{10}\,p$  (Mahalanobis vs the 266 used bursts)")
ax.set_title("How anomalous is each removed burst in feature space?\n"
             "(bars left of the dashed line = ordinary features -> removed for residual, not feature, reasons)")
ax.legend(fontsize=8, loc="lower right"); ax.grid(axis="x", alpha=0.25)
fig.tight_layout()
fig.savefig("outlier_mahalanobis.png", dpi=150)
print("\nwrote outlier_mahalanobis.png")
