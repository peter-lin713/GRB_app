#!/usr/bin/env python3
"""Epeak, GBM vs Swift BAT, on the corrected time+position GBM matches (flagship
281-GRB frame). Theil-Sen fit in the GBM -> Swift direction, no redshift used.
Error bars: GBM's 1-sigma flnc_comp_epeak errors. The Swift catalog has no error
column for bat_epeak, so the y values have none. Band: 68% bootstrap range of the fit.
Writes Plot_Output/epeak_transfer/epeak_gbm_to_swift_errorbars.png."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import theilslopes
from build_inputs import tap

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SURFACE, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
CPL, PL = '#2a78d6', '#eb6834'

d = pd.read_csv(os.path.join(REPO, 'GAM_supercomputer', 'pa_04_all_swiftpf_formula_generation',
                             'grb_xray_imputed.csv'), index_col=0)
gbm = pd.read_csv(os.path.join(REPO, 'Data', 'superlearner_training_theilsen_errcut_relative_gbmpf2_epeakmice.csv'),
                  index_col=0).log10Epeak.reindex(d.index)
B = pd.read_csv(os.path.join(HERE, 'inputs', 'bat_match_281_v2.csv')).set_index('key')
bat = np.log10(B.bat_epeak.where(B.bat_epeak > 0)).reindex(d.index)
model = B.bat_peak_model.reindex(d.index)
k = gbm.notna() & bat.notna()
x, y = gbm[k], bat[k]

# method='joint': intercept = median(y - m*x). scipy's default ('separate',
# median(y) - m*median(x)) leaves this line 0.10 dex below the median point.
m, b, lo, hi = theilslopes(y, x, method='joint')
r = np.corrcoef(x, y)[0, 1]
res = y - (m * x + b)
sig = 1.4826 * np.median(np.abs(res - np.median(res)))
m2, b2, lo2, hi2 = theilslopes(x, y, method='joint')
cw = pd.read_csv(os.path.join(HERE, 'gbm_crosswalk_time_all.csv'), index_col=0)
f = tap('fermigbrst', '')
f['name'] = f.name.str.strip()
f = f.drop_duplicates('name').set_index('name')
ep = f.flnc_comp_epeak.reindex(cw.name[x.index]).values
xerr = np.vstack([np.log10(ep) - np.log10(ep - f.flnc_comp_epeak_neg_err.reindex(cw.name[x.index]).values),
                  np.log10(ep + f.flnc_comp_epeak_pos_err.reindex(cw.name[x.index]).values) - np.log10(ep)])
assert np.allclose(np.log10(ep), x.values)
keep = x.index != '091024A'
r_no = np.corrcoef(x[keep], y[keep])[0, 1]
print(f'GBM -> Swift: log Ep(BAT) = {m:.3f} log Ep(GBM) {b:+.3f} | slope CI {lo:.2f}-{hi:.2f} | '
      f'r = {r:.3f} (without 091024A {r_no:.3f}) | scatter {sig:.3f} dex | N = {len(x)}')
print(f'Swift -> GBM: log Ep(GBM) = {m2:.3f} log Ep(BAT) {b2:+.3f} | slope CI {lo2:.2f}-{hi2:.2f}')
print('BAT model on the overlap:', model[k].value_counts().to_dict())

fig, ax = plt.subplots(figsize=(7, 6), dpi=200, facecolor=SURFACE)
ax.set_facecolor(SURFACE)
lim = np.array([1.4, 3.4])
ax.plot(lim, lim, color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=1)
grid = np.linspace(*lim, 100)
rng = np.random.default_rng(1)
curves = []
for _ in range(1000):
    i = rng.integers(0, len(x), len(x))
    mb, bb = theilslopes(y.values[i], x.values[i], method='joint')[:2]
    curves.append(mb * grid + bb)
ax.fill_between(grid, *np.percentile(curves, [16, 84], axis=0), color=GRID, lw=0, zorder=1, label='Fit, 68% bootstrap range')
ax.plot(lim, m * lim + b, color=INK, lw=2, zorder=2)
for lab, c, mk, sel in (('BAT model: cutoff power law', CPL, 'o', model[k] == 'CPL'),
                        ('BAT model: power law', PL, 's', model[k] == 'PL')):
    ax.errorbar(x[sel], y[sel], xerr=xerr[:, sel.values], fmt='none', ecolor=c, elinewidth=1.2, capsize=0, zorder=2.5)
    ax.scatter(x[sel], y[sel], s=62, c=c, marker=mk, edgecolors=SURFACE, linewidths=1.2, zorder=3,
               label=f'{lab} ({int(sel.sum())})')
ax.annotate('1:1', (3.28, 3.28), xytext=(-6, 6), textcoords='offset points', color=INK2, fontsize=10, ha='right')
ax.annotate('Theil–Sen fit', (2.95, m * 2.95 + b), xytext=(-4, -26), textcoords='offset points', color=INK, fontsize=10)
ax.annotate('091024A', (x['091024A'], y['091024A']), xytext=(0, 8), textcoords='offset points',
            color=INK2, fontsize=9.5, ha='center')
ax.text(0.03, 0.97, f'log Ep(Swift) = {m:.2f} · log Ep(GBM) {b:+.2f}\n'
                    f'slope 95% CI {lo:.2f}–{hi:.2f}   r = {r:.2f}\nscatter {sig:.2f} dex   N = {len(x)}',
        transform=ax.transAxes, va='top', fontsize=10.5, color=INK, linespacing=1.5)
ticks = [30, 100, 300, 1000, 2000]
ax.set_xticks(np.log10(ticks)); ax.set_xticklabels(ticks)
ax.set_yticks(np.log10(ticks)); ax.set_yticklabels(ticks)
ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect('equal')
ax.set_xlabel('GBM Epeak (keV), with 1σ errors', color=INK2, fontsize=11)
ax.set_ylabel('Swift BAT Epeak (keV), no catalog errors', color=INK2, fontsize=11)
ax.set_title('Epeak transfer, GBM to Swift', loc='left', color=INK, fontsize=13, fontweight='bold', pad=10)
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.tick_params(colors=INK2, length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(loc='lower right', frameon=False, fontsize=10, labelcolor=INK2)
out = os.path.join(REPO, 'Plot_Output', 'epeak_transfer')
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, 'epeak_gbm_to_swift_errorbars.png'), bbox_inches='tight', facecolor=SURFACE)
