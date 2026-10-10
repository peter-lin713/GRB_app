#!/usr/bin/env python3
"""
Plateau luminosity vs redshift, with the 15 M-estimator-removed bursts colored,
plus the flux-limit ("limiting luminosity") curve.

For each GRB the isotropic plateau-end luminosity is

    L_a = 4 * pi * D_L(z)^2 * F_a * k(z),        k(z) = (1 + z)^(beta - 1)

where F_a = 10**log10Fa is the plateau flux [erg/cm^2/s], beta is the X-ray
spectral index (the 'Beta' column), and k(z) is Dainotti's spectral k-correction
(what Prof. Dainotti wrote as 1 / (1+z)^(1 - beta) == (1+z)^(beta-1)). D_L is the
luminosity distance for flat LCDM (H0=70, Om=0.3), integrated numerically.

The "limiting luminosity" curve is the same relation evaluated at a fixed flux
limit F_lim and a representative spectral index beta_ref -- i.e. the faintest
luminosity detectable at redshift z. Any flux-limited sample must sit above it,
so it shows the Malmquist truncation of the L-z plane.

    !! Research-judgment choices to confirm with Prof. Dainotti !!
      - cosmology (H0, Om);
      - F_lim = the faintest observed plateau flux in the sample (a common choice;
        she sometimes uses a fixed instrumental threshold instead);
      - beta_ref = median beta for the single limiting curve (k is per-object beta
        for the points, but the curve needs one representative value).
    All three are set as named constants at the top so they are easy to change.

Saves 'luminosity_vs_redshift.png' next to this script.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# the M-est cut and the 281-row frame come from ONE place (single source of truth)
from simple_linear_flagship_check import d, d_cut

# ---- cosmology + flux-limit choices (edit here) -----------------------------
H0     = 70.0                         # km/s/Mpc
OM, OL = 0.30, 0.70                   # flat LCDM
C_KMS  = 2.99792458e5                 # km/s
MPC_CM = 3.0856775814913673e24        # cm per Mpc
D_H_CM = (C_KMS / H0) * MPC_CM        # Hubble distance in cm

def lum_dist_cm(z):
    """Luminosity distance [cm] for flat LCDM, vectorized over z."""
    zg = np.linspace(0.0, float(np.max(z)) + 1e-6, 20001)
    Einv = 1.0 / np.sqrt(OM * (1 + zg) ** 3 + OL)
    Dc = D_H_CM * np.concatenate(([0.0], np.cumsum((Einv[1:] + Einv[:-1]) * np.diff(zg) / 2)))
    return (1 + np.asarray(z)) * np.interp(z, zg, Dc)          # D_L = (1+z) * D_C

def luminosity(Fa, z, beta):
    """Isotropic plateau luminosity [erg/s] with spectral k-correction."""
    DL = lum_dist_cm(z)
    k  = (1 + z) ** (beta - 1)                                # == 1/(1+z)^(1-beta)
    return 4 * np.pi * DL ** 2 * Fa * k

# ---- luminosity of every burst ----------------------------------------------
Fa   = 10.0 ** d["log10Fa"].to_numpy()
z    = d["Redshift_crosscheck"].to_numpy()
beta = d["Beta"].to_numpy()
logL = np.log10(luminosity(Fa, z, beta))

removed = d.index.difference(d_cut.index)                     # the 15 M-est-removed
is_rem  = d.index.isin(removed)
print(f"M-est removed {is_rem.sum()} of {len(d)} GRBs")

# ---- limiting-luminosity curve ----------------------------------------------
# A single flux-limit curve has a FIXED shape (L_lim ∝ D_L(z)^2, set by cosmology);
# only its HEIGHT is ours to choose via F_lim. Here the limit is DEFINED by one
# specific faint burst: we use its plateau flux as F_lim and its own spectral index
# as beta_ref, so the curve passes exactly through that GRB.
ANCHOR   = "050826A"
F_REP    = 10.0 ** d.loc[ANCHOR, "log10Fa"]                  # that burst's plateau flux
BETA_REF = float(d.loc[ANCHOR, "Beta"])                      # that burst's spectral index
z_anchor = float(d.loc[ANCHOR, "Redshift_crosscheck"])
zc       = np.linspace(0.02, z.max() * 1.02, 400)
logL_rep = np.log10(luminosity(np.full_like(zc, F_REP), zc, BETA_REF))

# second limit: the height that leaves ~5% of the sample below it -- i.e. a
# flux cut that removes the same ~5% fraction as the M-estimator does.
BETA_MED = float(np.median(beta))
F_5      = np.percentile(Fa, 5)
logL_5   = np.log10(luminosity(np.full_like(zc, F_5), zc, BETA_MED))

# ---- plot -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 6.5))
ax.scatter(z[~is_rem], logL[~is_rem], s=20, c="#9ec3e6", alpha=0.8,
           edgecolor="none", label=f"kept ({(~is_rem).sum()})")
ax.scatter(z[is_rem],  logL[is_rem],  s=55, c="#d7301f", edgecolor="k", lw=0.4,
           zorder=5, label=f"M-est removed ({is_rem.sum()})")
ax.plot(zc, logL_rep, "k--", lw=1.9,
        label=fr"limiting $L$  from {ANCHOR}  ($F_{{\rm lim}}=10^{{{np.log10(F_REP):.2f}}}$,  $\beta={BETA_REF:.2f}$)")
ax.plot(zc, logL_5, "--", color="#2166ac", lw=1.9,
        label=fr"limiting $L$  (5% cut: $F_{{\rm lim}}=10^{{{np.log10(F_5):.2f}}}$,  $\beta_{{\rm ref}}={BETA_MED:.2f}$)")
# mark the anchor burst
ia = d.index.get_loc(ANCHOR)
ax.scatter([z_anchor], [logL[ia]], s=140, facecolor="none", edgecolor="k",
           lw=1.8, zorder=6, label=f"{ANCHOR} (anchor)")

ax.set_xlabel("redshift  z")           # LINEAR z axis: L_lim(z) shows as a rising curve
ax.set_ylabel(r"$\log_{10}\, L_a$  [erg s$^{-1}$]  (plateau, k-corrected)")
ax.set_title("Plateau luminosity vs redshift — M-estimator-removed bursts in red")
ax.legend(loc="lower right", framealpha=0.95, fontsize=9)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("luminosity_vs_redshift.png", dpi=150)
print("wrote luminosity_vs_redshift.png")

# ---- quick numbers for the meeting ------------------------------------------
def summ(name, m):
    print(f"  {name:14s} n={m.sum():3d}  z median={np.median(z[m]):.2f}  "
          f"logL median={np.median(logL[m]):.2f}  logL range [{logL[m].min():.2f}, {logL[m].max():.2f}]")
print("\nsummary:")
summ("kept",    ~is_rem)
summ("removed", is_rem)
dist_dex = logL - np.interp(z, zc, logL_rep)                  # dex above the representative curve
print(f"\n  median distance above representative curve:  kept {np.median(dist_dex[~is_rem]):.2f} dex   "
      f"removed {np.median(dist_dex[is_rem]):.2f} dex")
print(f"  bursts on/below the {ANCHOR} curve: {(dist_dex < 0).sum()} of {len(z)} "
      f"({(dist_dex[is_rem] < 0).sum()} of them M-est-removed)")
dist5 = logL - np.interp(z, zc, logL_5)
print(f"  bursts on/below the 5%-cut curve:  {(dist5 < 0).sum()} of {len(z)} "
      f"({(dist5[is_rem] < 0).sum()} of them also M-est-removed)")
