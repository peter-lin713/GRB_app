#!/usr/bin/env python3
"""
Task 2 (Dr. Dainotti): use the Fermi GBM peak flux in place of the Swift peak flux,
since the GBM measurement is more precise. Then follow it up with more data.

We compare the substitution (Swift peak flux -> GBM peak flux) on two samples:
  - the original 65 GRBs in Dr. Dainotti's curated common list (common67_transfer.csv)
  - an enlarged 116-GRB sample from a date+position crosswalk between the full
    281-GRB redshift sample and the Fermi GBM catalog (gbm_crosswalk.py ->
    gbm_crosswalk_281.csv), merged with the Swift features.
Model: ridge on the 7 Swift plateau features, honest repeated 5-fold CV (r in z).
Output: task2_gbm_substitution.png

Result: the small-sample +0.05 (n=65, 73% of resamples positive) regresses to
+0.01 (n=116, 53% positive) -- it does not survive the larger sample.
"""
import warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=__file__.rsplit("/",1)[0]; norm=lambda s:str(s).replace("GRB","").replace(" ","").strip().upper()
rlz=lambda pr,o:np.corrcoef(10**pr-1,10**o-1)[0,1]
def paired(data,SW,gbmcol,reps=400):
    y=data.log10z.values; X0=data[SW].values; X1=data[[c if c!=SW[-2] else gbmcol for c in SW]].values  # SW[-2] is the peak-flux col
    cv=RepeatedKFold(n_splits=5,n_repeats=reps,random_state=0); a,b=[],[]
    for tr,te in cv.split(X0):
        m0=make_pipeline(StandardScaler(),RidgeCV(alphas=np.logspace(-2,3,30))).fit(X0[tr],y[tr])
        m1=make_pipeline(StandardScaler(),RidgeCV(alphas=np.logspace(-2,3,30))).fit(X1[tr],y[tr])
        a.append(rlz(m0.predict(X0[te]),y[te])); b.append(rlz(m1.predict(X1[te]),y[te]))
    a,b=np.array(a),np.array(b); dl=b-a; return dl.mean(),dl.std(),100*(dl>0).mean()

# original 65-GRB sample
c=pd.read_csv("/Users/spencergibson/Desktop/peter_code/new_data/common67_transfer.csv")
c=c[(c.flnc_comp_phtflux>0)].copy(); c["log10z"]=np.log10(c.Redshift_crosscheck+1); c["gbm_pf"]=np.log10(c.flnc_comp_phtflux)
SW_c=["log10Fa_X","log10Ta_X","Alpha_X","log10NH","PhotonIndex","log10PeakFlux","log10T90"]
d1,s1,p1=paired(c,SW_c,"gbm_pf")

# enlarged crosswalk sample
cx=pd.read_csv(HERE+"/gbm_crosswalk_281.csv")
d=pd.read_csv(HERE+"/../grb_imputed_281.csv",index_col=0); d["key"]=[norm(i) for i in d.index]; d["log10z"]=np.log10(d.Redshift_crosscheck+1)
m=d.merge(cx,on="key",how="inner"); m=m[m.flnc_comp_phtflux>0].copy(); m["gbm_pf"]=np.log10(m.flnc_comp_phtflux)
SW_m=["log10Fa","log10Ta","Alpha","log10NH","PhotonIndex","log10PeakFlux","log10T90"]
d2,s2,p2=paired(m,SW_m,"gbm_pf")
print(f"original  (n={len(c)}): delta={d1:+.4f} ± {s1:.4f}  ({p1:.0f}% positive)")
print(f"enlarged  (n={len(m)}): delta={d2:+.4f} ± {s2:.4f}  ({p2:.0f}% positive)")

fig,ax=plt.subplots(figsize=(7.6,5.2))
dm=[d1,d2]; ds=[s1,s2]; pct=[p1,p2]; labels=[f"original\n(n = {len(c)})",f"enlarged crosswalk\n(n = {len(m)})"]
ax.bar(range(2),dm,yerr=ds,capsize=6,color=["#9ec3e6","#2c7fb8"]); ax.axhline(0,color="k",lw=1)
ax.set_xticks(range(2)); ax.set_xticklabels(labels,fontsize=10); ax.set_ylim(-0.12,0.18)
ax.set_ylabel("change in honest CV r  (Swift peak flux -> GBM peak flux)")
ax.set_title("Task 2: substituting the GBM peak flux\nThe small-sample +0.05 regresses to +0.01 on a larger common sample",fontsize=11)
for i,(v,e,p) in enumerate(zip(dm,ds,pct)): ax.text(i,v+e+.006,f"{v:+.3f} ± {e:.2f}\n{p:.0f}% positive",ha="center",fontsize=9)
ax.grid(axis="y",alpha=.3); fig.tight_layout(); fig.savefig(HERE+"/task2_gbm_substitution.png",dpi=150)
print("wrote task2_gbm_substitution.png")
