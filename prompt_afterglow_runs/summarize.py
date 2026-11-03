#!/usr/bin/env python3
"""Headline numbers for each finished run, computed exactly as the pipeline's
correlation plots do (Result_plot_maker.R): Pearson r of predicted vs observed,
in log10(z+1) and in linear z, on the 10-fold CV results with and without the
2-sigma catastrophic-outlier clip. The flagship's reported 0.707 is the
'without catOutl' linear-z number (N = 206). Writes results_summary.csv."""
import os
import glob
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = {'flagship 0.707 (reference)':
       os.path.join(REPO, 'runs', 'theilsen_cleaned_linearz_5pct_mice20_RETRAIN_nesteddata')}
rows = []
dirs = dict(REF, **{os.path.basename(d): d for d in sorted(glob.glob(os.path.join(REPO, 'runs', 'pa_*')))})
for name, d in dirs.items():
    row = {'run': name}
    for tag in ('with', 'without'):
        f = os.path.join(d, 'Results', f'10fCVResults_{tag}_catOutl_correlation_plot.csv')
        if not os.path.exists(f):
            continue
        r = pd.read_csv(f, index_col=0)
        row[f'N_{tag}'] = len(r)
        row[f'r_z_{tag}'] = np.corrcoef(r.Zphot, r.Zspec)[0, 1]
        row[f'r_log_{tag}'] = np.corrcoef(r.InvZphot, r.InvZspec)[0, 1]
    fu = os.path.join(d, 'formula_used.txt')
    if os.path.exists(fu):
        row['formula'] = ' '.join(open(fu).read().split())
    if len(row) > 1:
        rows.append(row)
out = pd.DataFrame(rows)
out.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results_summary.csv'),
           index=False, float_format='%.4f')
pd.set_option('display.width', 200)
print(out.drop(columns=['formula'], errors='ignore').to_string(index=False, float_format=lambda v: f'{v:.3f}'))
