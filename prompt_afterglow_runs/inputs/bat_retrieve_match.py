#!/usr/bin/env python3
"""Retrieve Swift BAT catalog data from HEASARC and match it to the 281-GRB
redshift sample, with a date-stem fallback so pre-2010 bursts that the BAT
catalog names without a trailing letter (e.g. '050826' vs '050826A') are not
dropped. Produces bat_match_281_v2.csv (274/281 matched; 260 with photon index,
57 with a direct bat_epeak)."""
import urllib.request, urllib.parse, io, warnings, pandas as pd, numpy as np
warnings.filterwarnings("ignore")
from astropy.io.votable import parse_single_table
def tap(adql):
    url="https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?"+urllib.parse.urlencode({"REQUEST":"doQuery","LANG":"ADQL","QUERY":adql})
    return parse_single_table(io.BytesIO(urllib.request.urlopen(url,timeout=120).read())).to_table().to_pandas()
norm=lambda s:str(s).replace("GRB","").replace(" ","").strip().upper()
stem=lambda k:k[:6] if len(k)>=6 and k[:6].isdigit() else k
bba=tap("SELECT name,bat_1s_peak_flux,bat_photon_index,bat_photon_index_error,bat_spectral_model FROM swiftgrbba")
sg =tap("SELECT name,bat_epeak,bat_peak_model FROM swiftgrb")
bba["key"]=bba.name.map(norm); sg["key"]=sg.name.map(norm)
cat=bba.merge(sg[["key","bat_epeak","bat_peak_model"]],on="key",how="outer"); cat["stem"]=cat.key.map(stem)
by_full={k:i for i,k in enumerate(cat.key)}
sc=cat.stem.value_counts(); by_stem={s:i for i,s in enumerate(cat.stem) if sc.get(s,0)==1}
d=pd.read_csv("../grb_imputed_281.csv",index_col=0)
d["key"]=[norm(i) for i in d.index]; d["log10z"]=np.log10(d.Redshift_crosscheck+1); rows=[]
for _,r in d.iterrows():
    k=r.key; idx=by_full.get(k, by_stem.get(stem(k))); rr=r.to_dict()
    if idx is not None:
        c=cat.iloc[idx]
        for col in ["bat_1s_peak_flux","bat_photon_index","bat_photon_index_error","bat_spectral_model","bat_epeak","bat_peak_model"]:
            rr[col]=c[col]
    rows.append(rr)
pd.DataFrame(rows).to_csv("bat_match_281_v2.csv",index=False); print("wrote bat_match_281_v2.csv")
