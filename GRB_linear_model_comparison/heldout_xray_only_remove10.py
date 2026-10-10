#!/usr/bin/env python3
"""
Same as heldout_xray_optical_remove10.py but on the X-RAY-ONLY subset (is_optical==0,
223 GRBs). Of the 10 difficult points, 8 are x-ray (050826A, 051109B, 080916A,
111008A, 130606A, 161129A, 220101A, 220521A); 100508A and 160314A are optical and
are not in this subset, so 8 are removed here.
"""
import warnings, numpy as np, pandas as pd, patsy
warnings.filterwarnings("ignore")
import statsmodels.api as sm
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SRC="/Users/spencergibson/Desktop/peter_code/grb_hogg_mixture/Data/superlearner_training_ot_v3_errcut_relative_isoptical.csv"
d=pd.read_csv(SRC); d.index=d[d.columns[0]].astype(str).str.replace("GRB","",regex=False).str.strip()
d=d[d.is_optical==0].copy()                                   # X-ray only
for b in ["log10Fa","log10T90","log10PeakFlux","log10NH","PhotonIndex","log10Ta","Alpha"]: d[b+"Sqr"]=d[b]**2
d["log10z"]=np.log10(d.Redshift_crosscheck+1)
W=("log10z ~ (log10FaSqr+log10T90+log10Fa+log10PeakFlux)**2+log10NH+PhotonIndex+log10Ta+Alpha"
   "+log10NHSqr+log10PeakFluxSqr+PhotonIndexSqr+log10TaSqr+log10T90Sqr+AlphaSqr")
yall,Xall=patsy.dmatrices(W,d,return_type="matrix"); yall=np.asarray(yall).ravel(); Xall=np.asarray(Xall)
OUT10=['050826A','051109B','080916A','100508A','111008A','130606A','160314A','161129A','220101A','220521A']
is_out=d.index.isin(OUT10)
print(f"X-ray only: {len(d)} GRBs; difficult points present: {is_out.sum()} "
      f"({sorted(d.index[is_out])})")
rlz=lambda pr,o:np.corrcoef(10**pr-1,10**o-1)[0,1]
def ols(Xtr,ytr,Xte): b,*_=np.linalg.lstsq(Xtr,ytr,rcond=None); return Xte@b
def heldout(rows, nsplits=1000, seed=0):
    idx=np.where(rows)[0]; rng=np.random.default_rng(seed); R=[]
    for _ in range(nsplits):
        perm=rng.permutation(idx); k=int(0.2*len(idx)); te=perm[:k]; tr=perm[k:]
        w=sm.RLM(yall[tr],Xall[tr],M=sm.robust.norms.HuberT()).fit().weights
        keep=tr[w>np.quantile(w,0.05)]
        R.append(rlz(ols(Xall[keep],yall[keep],Xall[te]),yall[te]))
    return np.array(R)
def flagship(rows, reps=20, K=10, seed=1):
    """M-estimator removes lowest 5% from the FULL sample, then repeated K-fold CV."""
    idx=np.where(rows)[0]
    w=sm.RLM(yall[idx],Xall[idx],M=sm.robust.norms.HuberT()).fit().weights
    keep=idx[w>np.quantile(w,0.05)]; rng=np.random.default_rng(seed); rs=[]
    for _ in range(reps):
        fold=rng.permutation(np.arange(len(keep))%K); P=np.full(len(keep),np.nan)
        for k in range(K):
            te=fold==k; P[te]=ols(Xall[keep[~te]],yall[keep[~te]],Xall[keep[te]])
        rs.append(rlz(P,yall[keep]))
    return float(np.mean(rs))
full=heldout(np.ones(len(d),bool)); cut=heldout(~is_out)
flag=flagship(np.ones(len(d),bool)); print(f"flagship (remove from full set, CV) = {flag:.3f}")
bins=np.linspace(-0.1,1.0,45)
def stats(R):
    h,e=np.histogram(R,bins=bins); return R.mean(),np.median(R),(e[np.argmax(h)]+e[np.argmax(h)+1])/2
mf=stats(full); mc=stats(cut)
print(f"\nWITH the 8 (n={len(d)}):   mean {mf[0]:.3f}  median {mf[1]:.3f}  mode {mf[2]:.2f}")
print(f"WITHOUT the 8 (n={int((~is_out).sum())}): mean {mc[0]:.3f}  median {mc[1]:.3f}  mode {mc[2]:.2f}")
print(f"shift in mean: {mc[0]-mf[0]:+.3f}")
fig,ax=plt.subplots(2,1,figsize=(9,8),sharex=True)
for a,(R,st,ttl,col) in zip(ax,[(full,mf,f"WITH the difficult points  (x-ray only, {len(d)})","#f0a202"),
                                 (cut,mc,f"WITHOUT the 8 x-ray difficult points  ({int((~is_out).sum())})","#2c7fb8")]):
    a.hist(R,bins=bins,color=col,alpha=0.85)
    a.axvline(st[0],color="#08519c",lw=2,label=f"mean   = {st[0]:.2f}")
    a.axvline(st[1],color="#08519c",lw=1.6,ls=":",label=f"median = {st[1]:.2f}")
    a.axvline(st[2],color="#08519c",lw=1.6,ls="-.",label=f"mode   = {st[2]:.2f}")
    a.axvline(flag,color="#d7301f",lw=2,ls="--",label=f"flagship (remove from full set) = {flag:.2f}")
    a.set_title(ttl,fontsize=11); a.set_ylabel("count over 1000 splits"); a.legend(fontsize=9,loc="upper left")
    a.grid(axis="y",alpha=0.25)
ax[1].set_xlabel("held-out r (20% held-out set, linear z)")
fig.suptitle("X-ray only: M-estimator on training only, 80/20 splits  —  effect of removing the difficult points",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig("heldout_xray_only_remove10.png",dpi=150)
print("\nwrote heldout_xray_only_remove10.png")
