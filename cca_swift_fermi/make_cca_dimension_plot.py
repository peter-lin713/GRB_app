#!/usr/bin/env python3
"""
How do the CCA dimensions relate to the actual dataset columns?

Heatmap of structure correlations: Pearson r between each original column and
each CCA dimension, over the GRBs that have that column.
  Left:  classical CCA, Swift (9) <-> GBM + Epeak (5), bridge GRBs (Part A, A2).
         Redshift is shown for reference -- it is NOT used in the fit.
  Right: the best Part B model, multi-view probabilistic CCA, k = 2
         (views: Swift | optical plateau | GBM prompt | Epeak | z), fit on all
         318 GRBs; each GRB's latent = posterior mean from the views it has.
Latent signs are arbitrary (a flipped dimension is the same dimension).

Run from this folder: python3 make_cca_dimension_plot.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cca_lib as L

OUT = os.path.join(os.path.dirname(L.REPO), 'GRB_app', 'Plot_Output', 'cca_swift_fermi',
                   'cca_dimensions_vs_columns.png')

LABEL = {
    'log10T90': 'log T90', 'log10Fluence': 'log fluence', 'PhotonIndex': 'photon index',
    'log10PeakFlux': 'log peak flux', 'log10Fa': 'log Fa (X-ray)', 'log10Ta': 'log Ta (X-ray)',
    'Alpha': 'α (X-ray)', 'Beta': 'β (X-ray)', 'log10NH': 'log NH',
    'o_log10Fa': 'log Fa (opt)', 'o_log10Ta': 'log Ta (opt)', 'o_Alpha': 'α (opt)', 'o_Beta': 'β (opt)',
    'g_log10T90': 'log T90', 'g_log10Fluence': 'log fluence', 'g_Index': 'index',
    'g_log10PhtFlux': 'log photon flux', 'g_log10Epeak': 'log Epeak', 'log10z': 'log(1+z)',
}

# Diverging blue <-> gray <-> red (reference palette: poles + neutral midpoint)
CMAP = LinearSegmentedColormap.from_list('div', ['#184f95', '#6da7ec', '#f0efec', '#ec8a89', '#c4302f'])
INK, MUTED, SURFACE = '#1f1f1e', '#6b6a66', '#fcfcfb'


def corr_with(df, cols, Z):
    out = np.full((len(cols), Z.shape[1]), np.nan)
    for i, c in enumerate(cols):
        ok = df[c].notna().values
        for j in range(Z.shape[1]):
            out[i, j] = np.corrcoef(df[c].values[ok], Z[ok, j])[0, 1]
    return out


def classical_cca(d):
    b = d[d.bridge].dropna(subset=L.S_ALL + L.G_ALL)
    sx, sy = L.standardizer(b[L.S_ALL].values), L.standardizer(b[L.G_ALL].values)
    fx, fy, rho, _, _ = L.fit_cca(sx(b[L.S_ALL].values), sy(b[L.G_ALL].values), 2)
    U, V = fx(sx(b[L.S_ALL].values)), fy(sy(b[L.G_ALL].values))
    groups = [('Swift', L.S_ALL, corr_with(b, L.S_ALL, U)),
              ('Fermi GBM', L.G_ALL, corr_with(b, L.G_ALL, V)),
              ('not in fit', ['log10z'], corr_with(b, ['log10z'], U))]
    return groups, len(b), rho


def multiview_pcca(d):
    cols_views = [L.S_ALL, L.G_PROMPT, ['g_log10Epeak'], L.O_PLATEAU, ['log10z']]
    allcols = sum(cols_views, [])
    st = L.standardizer(d[allcols].values)
    X = st(d[allcols].values)
    views, i = [], 0
    for cv in cols_views:
        views.append(np.arange(i, i + len(cv))); i += len(cv)
    m = L.PCCA(2, seed=0).fit(X, views)
    obs = np.column_stack([~np.isnan(X[:, v]).any(1) for v in views])
    Z = np.zeros((len(X), 2))
    for pat in np.unique(obs, axis=0):
        rows = (obs == pat).all(1)
        vids = list(np.where(pat)[0])
        Z[rows] = m.latent(np.nan_to_num(X[rows]), vids)
    groups = [('Swift', L.S_ALL, corr_with(d, L.S_ALL, Z)),
              ('optical', L.O_PLATEAU, corr_with(d, L.O_PLATEAU, Z)),
              ('Fermi GBM', L.G_ALL, corr_with(d, L.G_ALL, Z)),
              ('redshift', ['log10z'], corr_with(d, ['log10z'], Z))]
    return groups, len(d)


def draw(ax, groups, title, sub):
    rows, labels, seps, glabels, y = [], [], [], [], 0
    for g, cols, M in groups:
        glabels.append((g, y, y + len(cols) - 1))
        rows.append(M); labels += [LABEL[c] for c in cols]
        y += len(cols); seps.append(y)
    M = np.vstack(rows)
    ax.imshow(M, cmap=CMAP, vmin=-1, vmax=1, aspect='auto')
    for (i, j), v in np.ndenumerate(M):
        ax.text(j, i, f'{v:+.2f}', ha='center', va='center', fontsize=8.5,
                color='white' if abs(v) > 0.6 else INK)
    for s in seps[:-1]:
        ax.axhline(s - 0.5, color=SURFACE, lw=4)
    ax.set_yticks(range(len(labels)), labels, fontsize=9, color=INK)
    ax.set_xticks([0, 1], ['dim 1', 'dim 2'], fontsize=9.5, color=INK)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for g, a, b in glabels:
        ax.annotate(g, xy=(1.04, (a + b) / 2), xycoords=('axes fraction', 'data'),
                    fontsize=9, color=MUTED, va='center', ha='left', rotation=0)
    ax.set_title(title, fontsize=11, color=INK, loc='left', pad=26)
    ax.text(0, 1.035, sub, transform=ax.transAxes, fontsize=8.5, color=MUTED, va='bottom')


def main():
    d = L.load_all()
    g1, n1, rho = classical_cca(d)
    g2, n2 = multiview_pcca(d)
    fig, axes = plt.subplots(1, 2, figsize=(11, 8.2), facecolor=SURFACE,
                             gridspec_kw=dict(wspace=0.85, height_ratios=[1]))
    for ax in axes:
        ax.set_facecolor(SURFACE)
    draw(axes[0], g1, 'Classical CCA (Swift ↔ GBM)',
         f'{n1} bridge GRBs · canonical r = {rho[0]:.2f}, {rho[1]:.2f} in-sample (0.68, 0.83 held out)')
    draw(axes[1], g2, 'Multi-view probabilistic CCA (best model)',
         f'{n2} GRBs · Swift + optical + GBM + z · 2 latent dims')
    fig.suptitle('What each CCA dimension is made of: correlation with the original columns',
                 x=0.06, ha='left', fontsize=12.5, color=INK, y=0.995)
    cax = fig.add_axes([0.30, 0.035, 0.40, 0.018])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(-1, 1)),
                      cax=cax, orientation='horizontal')
    cb.set_label('Pearson r between column and dimension (sign of a dimension is arbitrary)',
                 fontsize=8.5, color=MUTED)
    cb.ax.tick_params(labelsize=8, colors=MUTED, length=0)
    cb.outline.set_visible(False)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches='tight', facecolor=SURFACE)
    print('Wrote', OUT)


if __name__ == '__main__':
    main()
