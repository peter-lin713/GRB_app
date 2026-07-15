#!/usr/bin/env python3
"""
make_4d_plot_ot_v3.py -- 3D scatter (log10Fa, log10Ta, Alpha) colored by
Beta, on the OT v3 dataset, with outliers flagged and highlighted.

Outlier flag: modified (MAD-based) z-score > 3.5 on any of
log10Fa/log10Ta/Alpha/Beta -- the standard Iglewicz-Hoaglin threshold, robust
to the heavy tails these afterglow parameters already show.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
DATA_FILE  = os.path.join(REPO_DIR, 'Data', 'superlearner_training_ot_v3_errcut_relative.csv')
OUT_FILE   = os.path.join(REPO_DIR, 'Plot_Output', 'fundamental_plane_4d_ot_v3.png')

df = pd.read_csv(DATA_FILE, index_col=0)
df = df[['Alpha', 'log10Fa', 'log10Ta', 'Beta']].apply(pd.to_numeric, errors='coerce').dropna()
print(f'Plotting {len(df)} GRBs with complete Alpha/log10Fa/log10Ta/Beta')

# Modified z-score per dimension: 0.6745*(x - median) / MAD
def modified_z(s):
    med = s.median()
    mad = (s - med).abs().median()
    if mad == 0:
        return pd.Series(0, index=s.index)
    return 0.6745 * (s - med) / mad

z_scores = df.apply(modified_z)
max_abs_z = z_scores.abs().max(axis=1)
is_outlier = max_abs_z > 3.5

print(f'Flagged {is_outlier.sum()} outliers (modified z-score > 3.5 on any of the 4 dims):')
for grb in df.index[is_outlier]:
    row = df.loc[grb]
    zrow = z_scores.loc[grb]
    worst_dim = zrow.abs().idxmax()
    print(f'  {grb}: log10Fa={row.log10Fa:.2f} log10Ta={row.log10Ta:.2f} '
          f'Alpha={row.Alpha:.2f} Beta={row.Beta:.2f}  (worst dim: {worst_dim}, z={zrow[worst_dim]:.2f})')

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

normal = df[~is_outlier]
outl   = df[is_outlier]

sc = ax.scatter(
    normal['log10Fa'], normal['log10Ta'], normal['Alpha'],
    c=normal['Beta'], cmap='viridis', s=35, alpha=0.8,
    edgecolors='k', linewidths=0.3, label=f'Normal (N={len(normal)})'
)
ax.scatter(
    outl['log10Fa'], outl['log10Ta'], outl['Alpha'],
    c='red', marker='X', s=140, edgecolors='black', linewidths=0.8,
    label=f'Outlier (N={len(outl)})', zorder=5
)

for grb in outl.index:
    row = outl.loc[grb]
    ax.text(row['log10Fa'], row['log10Ta'], row['Alpha'], f'  {grb}',
            fontsize=6.5, color='darkred')

ax.set_xlabel(r'log$_{10}(F_a)$', fontsize=10, labelpad=10)
ax.set_ylabel(r'log$_{10}(T_a)$', fontsize=10, labelpad=10)
ax.set_zlabel(r'$\alpha$ (post-plateau decay slope)', fontsize=10, labelpad=10)
ax.set_title('GRB afterglow fundamental plane (OT v3 dataset), colored by '
             r'$\beta$' '\nred X = outlier (modified z-score > 3.5)',
             fontsize=11, pad=20)

cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
cbar.set_label(r'$\beta$ (spectral index)', fontsize=10)
ax.legend(loc='upper left', fontsize=8)

plt.tight_layout()
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
plt.savefig(OUT_FILE, dpi=150)
plt.close()
print(f'\nSaved to {OUT_FILE}')
