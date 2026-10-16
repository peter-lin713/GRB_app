#!/usr/bin/env python3
"""
Task 1 (Dr. Dainotti): get an E_peak from the Swift BAT catalog itself and add it
(plus the BAT 1-second peak photon index and peak flux) to the redshift estimator.

Run bat_retrieve_match.py first to produce bat_match_281_v2.csv (the 281-GRB
redshift sample matched to the HEASARC Swift BAT catalogs, with a date-stem
fallback: 274/281 matched, 255 with a 1-s peak photon index, 57 with a direct
bat_epeak).

E_peak: use the catalog bat_epeak where available (57 GRBs), otherwise estimate it
from the BAT photon index via log10(E_peak) = -0.402*Gamma + 2.882 (calibrated on
BAT bursts that have both; never uses redshift, so the CV stays honest). This is a
MEASURED Swift-side E_peak for the full sample -- no transfer/imputation, so it
avoids the circularity of the Fermi transfer approach.

We report the directly-measured E_peak separately from the index-derived proxy, and
avoid putting both the photon index and the derived E_peak in the model (nearly
collinear). Model: ridge on the winning-GLM feature set, honest repeated 5-fold CV,
r in linear z. Output: task1_bat_epeak.png
"""
import warnings, numpy as np, pandas as pd, patsy
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=__file__.rsplit("/",1)[0]
m=pd.read_csv(HERE+"/bat_match_281_v2.csv")
m=m[m.bat_photon_index.notna() & (m.bat_1s_peak_flux>0)].copy()
m["is_direct"]=m.bat_epeak.notna()&(m.bat_epeak>0)
m["bat_epeak_est"]=np.where(m.is_direct, m.bat_epeak, 10**(-0.402*m.bat_photon_index+2.882))
m["log10_bat_epeak"]=np.log10(m.bat_epeak_est); m["log10_bat_pf"]=np.log10(m.bat_1s_peak_flux)
m=m[np.isfinite(m.log10_bat_epeak)&np.isfinite(m.log10_bat_pf)].copy()
for b in ["log10Fa","log10T90","log10PeakFlux","log10NH","PhotonIndex","log10Ta","Alpha"]: m[b+"Sqr"]=m[b]**2
y=m.log10z.values; n=len(m); nd=int(m.is_direct.sum())
print(f"Task 1: {n} GRBs with BAT prompt data; {nd} with a DIRECT bat_epeak\n")
for c in ["bat_photon_index","log10_bat_pf","log10_bat_epeak"]:
    r,p=pearsonr(m[c],m.log10z); print(f"  corr {c:18s} r={r:+.3f} (p={p:.2g})")
dsub=m[m.is_direct]; rr,pp=pearsonr(np.log10(dsub.bat_epeak),dsub.log10z)
print(f"  corr DIRECT bat_epeak only (n={nd}): r={rr:+.3f} (p={pp:.2g})")
rlz=lambda pr,o:np.corrcoef(10**pr-1,10**o-1)[0,1]
W=("(log10FaSqr+log10T90+log10Fa+log10PeakFlux)**2+log10NH+PhotonIndex+log10Ta+Alpha"
   "+log10NHSqr+log10PeakFluxSqr+PhotonIndexSqr+log10TaSqr+log10T90Sqr+AlphaSqr")
Xw=np.asarray(patsy.dmatrix(W,m,return_type="matrix"))[:,1:]
def paired(add,reps=300):
    cv=RepeatedKFold(n_splits=5,n_repeats=reps,random_state=0); a,b=[],[]
    X1=np.column_stack([Xw]+add) if add else Xw
    for tr,te in cv.split(Xw):
        m0=make_pipeline(StandardScaler(),RidgeCV(alphas=np.logspace(-2,3,30))).fit(Xw[tr],y[tr])
        m1=make_pipeline(StandardScaler(),RidgeCV(alphas=np.logspace(-2,3,30))).fit(X1[tr],y[tr])
        a.append(rlz(m0.predict(Xw[te]),y[te])); b.append(rlz(m1.predict(X1[te]),y[te]))
    a,b=np.array(a),np.array(b); return a.mean(),b.mean(),(b-a).mean(),(b-a).std()
BP=m.bat_photon_index.values; PF=m.log10_bat_pf.values; EP=m.log10_bat_epeak.values
res={}; print("\nWinning-GLM feature set (ridge, honest CV):")
for name,cols in [("photon idx",[BP]),("peak flux",[PF]),("idx + peak flux",[BP,PF]),("+ derived Epeak",[EP,PF])]:
    a,bb,dm,ds=paired(cols); res[name]=(dm,ds); print(f"  +BAT {name:16s}: {a:.3f}->{bb:.3f}  d={dm:+.4f} ± {ds:.4f}")
base=paired([])[0]
fig,ax=plt.subplots(1,2,figsize=(12.6,5.3))
dd=m[m.is_direct]; de=m[~m.is_direct]
ax[0].scatter(de.Redshift_crosscheck,de.bat_epeak_est,s=24,c="#9ec3e6",alpha=.6,edgecolor="none",label=f"derived from index ({len(de)})")
ax[0].scatter(dd.Redshift_crosscheck,dd.bat_epeak_est,s=40,c="#d7301f",edgecolor="white",lw=.4,label=f"directly measured ({len(dd)})")
ax[0].set_yscale("log"); ax[0].set_xlabel("Redshift z"); ax[0].set_ylabel("Swift BAT E_peak (keV)")
ax[0].set_title(f"Swift BAT E_peak vs redshift  (n={n})")
ax[0].text(0.96,0.04,f"direct E_peak vs z: r={rr:+.2f} (n={nd})",transform=ax[0].transAxes,ha="right",va="bottom",fontsize=9,bbox=dict(boxstyle="round",fc="#f5f5f5",ec="#bbb"))
ax[0].legend(loc="upper left",fontsize=8.5); ax[0].grid(alpha=.3,which="both")
labels=list(res); dm=[res[k][0] for k in labels]; ds=[res[k][1] for k in labels]
ax[1].bar(range(len(labels)),dm,yerr=ds,capsize=5,color=["#9ec3e6","#9ec3e6","#2c7fb8","#2c7fb8"]); ax[1].axhline(0,color="k",lw=1)
ax[1].set_xticks(range(len(labels))); ax[1].set_xticklabels(labels,fontsize=8.5); ax[1].set_ylabel("change in honest CV r (ridge)")
ax[1].set_title(f"Adding BAT prompt variables  (baseline r={base:.2f})"); ax[1].grid(axis="y",alpha=.3)
fig.suptitle("Task 1: Swift BAT prompt E_peak / peak flux added to the redshift estimator",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(HERE+"/task1_bat_epeak.png",dpi=150); m.to_csv(HERE+"/task1_bat_sample.csv",index=False)
print(f"\nbaseline r={base:.3f}; wrote task1_bat_epeak.png")
