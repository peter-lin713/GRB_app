#!/usr/bin/env python3
"""
Task 2, follow-up: instead of restricting to the Swift-Fermi common sample, keep the
WHOLE 281-GRB redshift sample and, for the subset that has a Fermi GBM peak flux,
substitute that GBM value in place of the Swift peak flux (everyone else stays on Swift).
Then measure the honest CV r of the winning-GLM estimator.

The GBM peak flux (photon flux) is in different units from the Swift log10PeakFlux, so
we first cross-calibrate GBM onto the Swift scale with an OLS fit on the overlap
(Swift_pf ~ a + b*GBM_pf; never uses redshift, so the CV stays honest), then substitute.

Model: winning GLM, ridge, honest repeated 5-fold CV (r in linear z).
Result: no improvement -- marginally negative, within the noise
  65-set imputed  (63 of 281): 0.499 -> 0.491   d=-0.008 ± 0.024  (36% positive)
  116-set imputed (116 of 281): 0.499 -> 0.480   d=-0.019 ± 0.037  (28% positive)
"""
import warnings, numpy as np, pandas as pd, patsy
warnings.filterwarnings("ignore")
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
HERE=__file__.rsplit("/",1)[0]
norm=lambda s:str(s).replace("GRB","").replace(" ","").strip().upper()
rlz=lambda pr,o:np.corrcoef(10**pr-1,10**o-1)[0,1]

d=pd.read_csv(HERE+"/../grb_imputed_281.csv",index_col=0).copy()
d["key"]=[norm(i) for i in d.index]; d["log10z"]=np.log10(d.Redshift_crosscheck+1)
cx=pd.read_csv(HERE+"/gbm_crosswalk_281.csv"); cx=cx[cx.flnc_comp_phtflux>0]
d["gbm"]=d.key.map(dict(zip(cx.key,np.log10(cx.flnc_comp_phtflux))))
c65=pd.read_csv("/Users/spencergibson/Desktop/peter_code/new_data/common67_transfer.csv")
k65=set(norm(x) for x in c65.GRB_Name)|set(norm(x) for x in c65.GRBweb_name)
sub116=d.gbm.notna(); sub65=sub116 & d.key.isin(k65)

W=("(log10FaSqr+log10T90+log10Fa+log10PeakFlux)**2+log10NH+PhotonIndex+log10Ta+Alpha"
   "+log10NHSqr+log10PeakFluxSqr+PhotonIndexSqr+log10TaSqr+log10T90Sqr+AlphaSqr")
def design(frame):
    f=frame.copy()
    for b in ["log10Fa","log10T90","log10PeakFlux","log10NH","PhotonIndex","log10Ta","Alpha"]: f[b+"Sqr"]=f[b]**2
    return np.asarray(patsy.dmatrix(W,f,return_type="matrix"))[:,1:]
y=d.log10z.values
def cvr(X,reps=400):
    cv=RepeatedKFold(n_splits=5,n_repeats=reps,random_state=0); out=[]
    for tr,te in cv.split(X):
        m=make_pipeline(StandardScaler(),RidgeCV(alphas=np.logspace(-2,3,30))).fit(X[tr],y[tr])
        out.append(rlz(m.predict(X[te]),y[te]))
    return np.array(out)

base=cvr(design(d))
def imputed(mask):
    sw=d.loc[mask,"log10PeakFlux"].values; gb=d.loc[mask,"gbm"].values
    b1,b0=np.polyfit(gb,sw,1); f=d.copy(); f.loc[mask,"log10PeakFlux"]=b1*gb+b0
    return cvr(design(f)), np.corrcoef(b1*gb+b0,sw)[0,1]
print(f"full sample n={len(d)}; baseline (all Swift) r = {base.mean():.4f}\n")
for name,mask in [("65-set",sub65),("116-set",sub116)]:
    imp,rfit=imputed(mask); dl=imp-base
    print(f"impute GBM for {name} ({mask.sum()} of {len(d)}): {base.mean():.4f} -> {imp.mean():.4f}"
          f"  d={dl.mean():+.4f} ± {dl.std():.4f}  ({100*(dl>0).mean():.0f}% positive)  [GBM->Swift fit r={rfit:.2f}]")
