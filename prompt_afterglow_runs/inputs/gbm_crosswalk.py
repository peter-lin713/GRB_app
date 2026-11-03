#!/usr/bin/env python3
"""Build a date+position crosswalk between the 281-GRB redshift sample and the
Fermi GBM catalog (fermigbrst), to enlarge the Swift-Fermi common sample beyond
the ~65 in Dr. Dainotti's curated list. Swift positions come from swiftgrbba;
matching is date-primary (unique GBM burst that day -> match; multiple that day ->
nearest sky position within 8 deg). Produces gbm_crosswalk_281.csv (122 matched,
~117 with a GBM peak flux). Requires network (HEASARC)."""
import urllib.request, urllib.parse, io, warnings, pandas as pd, numpy as np
warnings.filterwarnings("ignore")
from astropy.io.votable import parse_single_table
def tap(q):
    u="https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?"+urllib.parse.urlencode({"REQUEST":"doQuery","LANG":"ADQL","QUERY":q})
    return parse_single_table(io.BytesIO(urllib.request.urlopen(u,timeout=180).read())).to_table().to_pandas()
norm=lambda s:str(s).replace("GRB","").replace(" ","").strip().upper()
stem=lambda k:str(k)[:6] if len(str(k))>=6 and str(k)[:6].isdigit() else str(k)
sb=tap("SELECT name,ra,dec FROM swiftgrbba"); sb["key"]=sb.name.map(norm); sb["stem"]=sb.key.map(stem)
bf={k:i for i,k in enumerate(sb.key)}; sc=sb.stem.value_counts(); bs={s:i for i,s in enumerate(sb.stem) if sc.get(s,0)==1}
d=pd.read_csv("../grb_imputed_281.csv",index_col=0); d["key"]=[norm(i) for i in d.index]; d["stem"]=d.key.map(stem)
ra=[];dec=[]
for k in d.key:
    i=bf.get(k,bs.get(stem(k))); ra.append(sb.ra.iloc[i] if i is not None else np.nan); dec.append(sb.dec.iloc[i] if i is not None else np.nan)
d["ra"]=ra; d["dec"]=dec
g=tap("SELECT name,ra,dec,flnc_comp_phtflux,pflx_comp_phtflux,flnc_comp_index,flnc_comp_ergflnc FROM fermigbrst")
g["key"]=g.name.map(norm); g["date"]=g.key.map(stem); g=g[g.ra.notna()].copy()
sep=lambda r1,d1,r2,d2:(lambda a,b,c,e:np.degrees(np.arccos(np.clip(np.sin(b)*np.sin(e)+np.cos(b)*np.cos(e)*np.cos(a-c),-1,1))))(*map(np.radians,[r1,d1,r2,d2]))
rows=[]
for _,r in d.iterrows():
    cand=g[g.date==r.stem]
    if len(cand)==0: continue
    if len(cand)==1: c=cand.iloc[0]; sp=np.nan
    else:
        if not np.isfinite(r.ra): continue
        s=sep(r.ra,r.dec,cand.ra.values,cand.dec.values); j=int(np.argmin(s))
        if s[j]>=8.0: continue
        c=cand.iloc[j]; sp=round(float(s[j]),2)
    rows.append({"key":r.key,"sep_deg":sp,"n_sameday":len(cand),"flnc_comp_phtflux":c.flnc_comp_phtflux,
                 "pflx_comp_phtflux":c.pflx_comp_phtflux,"flnc_comp_index":c.flnc_comp_index,"flnc_comp_ergflnc":c.flnc_comp_ergflnc})
pd.DataFrame(rows).to_csv("gbm_crosswalk_281.csv",index=False); print(f"wrote gbm_crosswalk_281.csv ({len(rows)} matched)")
