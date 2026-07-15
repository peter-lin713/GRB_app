#!/usr/bin/env python3
"""
make_4d_plot.py -- 3D scatter (log10Fa, log10Ta, Alpha) colored by Alpha.

Intro-section figure: positions each GRB in (log10Fa, log10Ta, Alpha) space
(the Dainotti-style afterglow "fundamental plane" triplet) and color-codes by
Alpha as the 4th dimension.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
DATA_FILE  = os.path.join(REPO_DIR, 'Data', 'superlearner_training_emcee_errcut_relative.csv')
OUT_FILE   = os.path.join(REPO_DIR, 'Plot_Output', 'fundamental_plane_4d.png')

df = pd.read_csv(DATA_FILE)
df = df[['Alpha', 'log10Fa', 'log10Ta']].apply(pd.to_numeric, errors='coerce').dropna()
print(f'Plotting {len(df)} GRBs with complete Alpha/log10Fa/log10Ta')

fig = plt.figure(figsize=(9, 7.5))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(
    df['log10Fa'], df['log10Ta'], df['Alpha'],
    c=df['Alpha'], cmap='viridis', s=35, alpha=0.85,
    edgecolors='k', linewidths=0.3
)

ax.set_xlabel(r'log$_{10}(F_a)$', fontsize=10, labelpad=10)
ax.set_ylabel(r'log$_{10}(T_a)$', fontsize=10, labelpad=10)
ax.set_zlabel(r'$\alpha$ (post-plateau decay slope)', fontsize=10, labelpad=10)
ax.set_title('GRB afterglow fundamental plane, colored by ' r'$\alpha$',
             fontsize=12, pad=20)

cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
cbar.set_label(r'$\alpha$ (decay slope)', fontsize=10)

plt.tight_layout()
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
plt.savefig(OUT_FILE, dpi=150)
plt.close()
print(f'Saved to {OUT_FILE}')
