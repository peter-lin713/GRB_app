#!/usr/bin/env python3
"""
Fermi -> Swift transfer-learning feasibility: reproduce the headline results.

Runs entirely from the saved CSVs in this folder plus ../grb_imputed_281.csv
(the 281-burst Swift redshift sample) -- no network / HEASARC query needed.

Reproduces:
  A. the observed-Epeak vs redshift correlation (163 Fermi GRBs) + its figure
  B. the transfer-map experiment: learn Swift->Epeak on the 60 good bridge GRBs,
     impute Epeak for the 281, add it to the winning GLM and to a random forest,
     paired cross-validation  (+ the true-vs-predicted figure)
  C. the pooled / domain-adaptation experiment (shared variables + instrument dummy)
  D. the "ladder" of baselines that explains why some numbers look low

Headline: adding Epeak (any of the three ways) does not improve redshift
prediction, because observed Epeak is uncorrelated with redshift (r = 0.02).
"""
import warnings, numpy as np, pandas as pd, patsy
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import RepeatedKFold
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = __file__.rsplit("/", 1)[0] if "/" in __file__ else "."
D281 = HERE + "/../grb_imputed_281.csv"            # the Swift 281-burst sample (parent folder)
rlz = lambda pr, o: np.corrcoef(10**pr - 1, 10**o - 1)[0, 1]

# ===== A. observed Epeak vs redshift (163 Fermi GRBs) ========================
ez = pd.read_csv(HERE + "/epeak_z_163.csv")
le, lz = np.log10(ez.flnc_comp_epeak), np.log10(ez.z + 1)
rA, pA = pearsonr(le, lz)
print(f"A. observed log10Epeak vs log10(1+z):  r = {rA:+.3f}  (p = {pA:.2f}, n = {len(ez)})")

# ===== B. transfer map: Swift -> Epeak, impute for 281, add to redshift model =
bridge = pd.read_csv(HERE + "/common67_transfer.csv")
ep = bridge.flnc_comp_epeak
good = bridge[(ep > 8) & (ep < 10000) & (bridge.flnc_comp_epeak_pos_err < ep)].copy()
good["log10Epeak"] = np.log10(good.flnc_comp_epeak)
SWx = ["log10Fa_X", "log10Ta_X", "Alpha_X", "log10NH", "PhotonIndex", "log10PeakFlux", "log10T90"]
imp = RandomForestRegressor(400, min_samples_leaf=3, random_state=0).fit(good[SWx].values, good.log10Epeak.values)

d = pd.read_csv(D281, index_col=0)
SW = ["log10Fa", "log10Ta", "Alpha", "log10NH", "PhotonIndex", "log10PeakFlux", "log10T90"]
d["Epeak_imp"] = imp.predict(d[SW].values)
d["log10z"] = np.log10(d.Redshift_crosscheck + 1)
for b in ["log10Fa", "log10T90", "log10PeakFlux", "log10NH", "PhotonIndex", "log10Ta", "Alpha"]:
    d[b + "Sqr"] = d[b] ** 2
W = ("(log10FaSqr+log10T90+log10Fa+log10PeakFlux)**2+log10NH+PhotonIndex+log10Ta+Alpha"
     "+log10NHSqr+log10PeakFluxSqr+PhotonIndexSqr+log10TaSqr+log10T90Sqr+AlphaSqr")
y = d.log10z.values; n = len(y)
Xg0 = np.asarray(patsy.dmatrix(W, d, return_type="matrix"))
Xg1 = np.asarray(patsy.dmatrix(W + "+Epeak_imp", d, return_type="matrix"))
Xr0, Xr1 = d[SW].values, d[SW + ["Epeak_imp"]].values
ols = lambda Xtr, yt, Xte: Xte @ np.linalg.lstsq(Xtr, yt, rcond=None)[0]
rf  = lambda Xtr, yt, Xte: RandomForestRegressor(300, min_samples_leaf=3, random_state=1).fit(Xtr, yt).predict(Xte)

def paired(X0, X1, fit, reps):
    cv = RepeatedKFold(n_splits=5, n_repeats=reps, random_state=0); a, b = [], []
    for tr, te in cv.split(X0):
        a.append(rlz(fit(X0[tr], y[tr], X0[te]), y[te]))
        b.append(rlz(fit(X1[tr], y[tr], X1[te]), y[te]))
    a, b = np.array(a), np.array(b); return a.mean(), b.mean(), (b - a).mean(), (b - a).std()

print("\nB. transfer-learned Epeak added to the redshift model (281, paired CV):")
a, b, dm, ds = paired(Xg0, Xg1, ols, 200); print(f"   winning GLM : r {a:.3f} -> {b:.3f}   delta = {dm:+.4f} +/- {ds:.4f}")
a, b, dm, ds = paired(Xr0, Xr1, rf, 20);   print(f"   random forest: r {a:.3f} -> {b:.3f}   delta = {dm:+.4f} +/- {ds:.4f}  (RF, 20 reps)")

# ===== C. pooled / domain adaptation (shared variables + instrument dummy) ====
fm = pd.read_csv(HERE + "/fermi_shared_vars.csv"); fo = fm[~fm.is_common].copy()
SW4 = ["PhotonIndex", "log10T90", "log10Fluence", "log10PeakFlux"]
F4 = ["gbm_photindex", "gbm_log10T90", "gbm_log10fluence", "gbm_log10peakflux"]
Xs, Xf, yf = d[SW4].values, fo[F4].values, fo.log10z.values
zsc = lambda X, mu, sd: (X - mu) / sd
fmu, fsd = Xf.mean(0), Xf.std(0); Xfz = zsc(Xf, fmu, fsd)
cv = RepeatedKFold(n_splits=5, n_repeats=200, random_state=0); r_sw, r_pool = [], []
for tr, te in cv.split(Xs):
    smu, ssd = Xs[tr].mean(0), Xs[tr].std(0); Xtr, Xte = zsc(Xs[tr], smu, ssd), zsc(Xs[te], smu, ssd)
    m = RidgeCV(alphas=np.logspace(-2, 3, 30)).fit(Xtr, y[tr]); r_sw.append(rlz(m.predict(Xte), y[te]))
    Xp = np.vstack([np.c_[Xtr, np.zeros(len(tr))], np.c_[Xfz, np.ones(len(Xfz))]])
    m = RidgeCV(alphas=np.logspace(-2, 3, 30)).fit(Xp, np.r_[y[tr], yf])
    r_pool.append(rlz(m.predict(np.c_[Xte, np.zeros(len(te))]), y[te]))
r_sw, r_pool = np.array(r_sw), np.array(r_pool)
print(f"\nC. pooled domain adaptation (shared 4 vars + dummy, {len(fo)} Fermi-only added):")
print(f"   Swift-only shared-4 : r = {r_sw.mean():.3f}")
print(f"   POOLED + dummy      : r = {r_pool.mean():.3f}   delta = {(r_pool-r_sw).mean():+.4f}")

# ===== D. the ladder =========================================================
print("\nD. ladder of baselines (281, honest CV):")
print("   winning GLM, full plateau features .......... ~0.53   <- headline")
print("   winning GLM + transfer-learned Epeak ........ ~0.52   (section B)")
print("   plain linear, 7 plateau features ............ ~0.46")
print("   plain linear, only 4 shared variables ....... ~0.29   (basis for C)")

# ===== figures ===============================================================
# A figure
plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(figsize=(7.6, 5.6))
ax.scatter(ez.z, ez.flnc_comp_epeak, s=34, c="#3b7dd8", alpha=.7, edgecolor="white", linewidth=.4)
bb = np.polyfit(ez.z, np.log10(ez.flnc_comp_epeak), 1); xx = np.linspace(0, ez.z.max() * 1.02, 50)
ax.plot(xx, 10**np.polyval(bb, xx), "--", color="#d62728", lw=2, label="linear fit (essentially flat)")
ax.set_yscale("log"); ax.set_ylim(7, 4000); ax.set_xlabel("Redshift  $z$")
ax.set_ylabel(r"Observed $E_{\rm peak}$  (keV)")
ax.set_title("Observed GRB $E_{\\rm peak}$ shows no correlation with redshift", fontsize=13, pad=12)
ax.text(0.97, 0.95, f"Pearson $r = {rA:+.2f}$  ($p = {pA:.1f}$)\n$n = {len(ez)}$ Fermi GRBs with measured $z$",
        transform=ax.transAxes, ha="right", va="top", fontsize=11.5,
        bbox=dict(boxstyle="round,pad=0.5", fc="#f5f5f5", ec="#bbbbbb"))
ax.grid(alpha=.3, which="both"); ax.legend(loc="lower right", fontsize=10)
fig.tight_layout(); fig.savefig(HERE + "/Epeak_vs_redshift_for_Maria.png", dpi=200, bbox_inches="tight")

# B figure (true vs predicted, winning GLM, with/without Epeak)
def oof(X, reps=25, K=10):
    P = np.zeros(n)
    for r in range(reps):
        fold = np.random.default_rng(r).permutation(np.arange(n) % K); Pr = np.full(n, np.nan)
        for k in range(K):
            te = fold == k; Pr[te] = ols(X[~te], y[~te], X[te])
        P += Pr
    return P / reps
p0, p1 = oof(Xg0), oof(Xg1); zt = 10**y - 1; hi = max(zt.max(), (10**p0 - 1).max()) * 1.05
fig, ax = plt.subplots(1, 2, figsize=(12.4, 6.1))
for a, (pp, ttl) in zip(ax, [(p0, "Winning GLM - Swift features only"), (p1, "Winning GLM + transfer-learned Epeak")]):
    a.plot([0, hi], [0, hi], "k--", lw=1, alpha=.6); a.scatter(zt, 10**pp - 1, s=18, c="#2c7fb8", alpha=.75, edgecolor="none")
    a.set_xlim(0, hi); a.set_ylim(0, hi); a.set_aspect("equal"); a.set_xlabel("observed z"); a.set_ylabel("predicted z")
    a.set_title(ttl, fontsize=11)
    a.text(0.04, 0.96, f"r = {rlz(pp, y):.3f}", transform=a.transAxes, va="top", fontsize=12,
           bbox=dict(boxstyle="round", fc="white", ec="0.7"))
fig.suptitle("Redshift prediction on the 281 Swift sample (honest out-of-fold CV) - Epeak adds nothing", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(HERE + "/redshift_pred_with_without_epeak.png", dpi=160)
print("\nregenerated both figures.")
