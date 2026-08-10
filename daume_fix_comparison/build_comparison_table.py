#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = '/Users/petelin/Desktop/test/GRB_app'
COLOR_HEAD = '#2a2a28'

# [Dataset, Daume?, N, r(z), Sigma, RMS, Bias, NMAD]
rows = [
    ['Standing (single-var emcee)', 'Yes (genuine)', '208', '0.670', '0.957', '0.97', '0.17', '0.862'],
    ['Standing (single-var emcee)', 'No',            '206', '0.663', '0.956', '0.97', '0.17', '0.812'],
    ['Multivariate emcee',          'Yes (genuine)', '204', '0.711', '0.906', '0.93', '0.19', '0.792'],
    ['Multivariate emcee',          'No',            '202', '0.714', '0.901', '0.92', '0.21', '0.795'],
    ['Theil-Sen multivariate',      'Yes (genuine)', '207', '0.678', '0.957', '0.97', '0.18', '0.812'],
    ['Theil-Sen multivariate',      'No',            '205', '0.687', '0.938', '0.96', '0.18', '0.817'],
    ['Theil-Sen single-var',        'Yes (genuine)', '216', '0.614', '0.920', '0.93', '0.15', '0.916'],
    ['Theil-Sen single-var',        'No',            '215', '0.628', '0.909', '0.92', '0.14', '0.924'],
]

fig, ax = plt.subplots(figsize=(13, 8.5), facecolor='white')
ax.axis('off')
col_labels = ['Dataset', 'Daume?', 'N', 'r (z)', 'Sigma', 'RMS', 'Bias', 'NMAD']
col_widths = [0.30, 0.15, 0.07, 0.09, 0.09, 0.08, 0.08, 0.09]

tbl = ax.table(cellText=rows, colLabels=col_labels, cellLoc='center', loc='center',
               colWidths=col_widths, bbox=[0.0, 0.0, 1.0, 1.0])
tbl.auto_set_font_size(False)
tbl.set_fontsize(11)

for (r, c), cell in tbl.get_celld().items():
    cell.set_edgecolor('#d8d7d1')
    if r == 0:
        cell.set_facecolor(COLOR_HEAD)
        cell.set_text_props(color='white', fontweight='bold')
        continue
    # Pair rows (Yes/No) share a light band so the with/without comparison reads as one unit
    pair_idx = (r - 1) // 2
    cell.set_facecolor('#f4f3ef' if pair_idx % 2 == 0 else 'white')
    if c == 0:
        cell.set_text_props(ha='left', fontweight='bold', color='#222')
        cell.PAD = 0.02
    elif c == 1:
        is_yes = rows[r - 1][1].startswith('Yes')
        cell.set_text_props(color='#1a6e3c' if is_yes else '#666', fontweight='bold' if is_yes else 'normal')
    else:
        cell.set_text_props(color='#222')

tbl.scale(1, 2.6)
ax.set_title('Genuine Daume domain-adaptation: with vs. without\n'
             '(same standing formula, same "Without 2σ outliers" subset; "Yes" rows use the fixed\n'
             'superlearner_daume_working.R where is_optical columns actually reach the model)',
             fontsize=12, color='#111', pad=16)

out = f'{REPO}/daume_fix_comparison/plots/daume_comparison_table.png'
plt.savefig(out, dpi=150, facecolor='white', bbox_inches='tight')
print(f'Saved {out}')
