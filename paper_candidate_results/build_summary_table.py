#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = '/Users/petelin/Desktop/test/GRB_app'
COLOR_HEAD = '#2a2a28'

# [Role, Scheme, N, r(z), Sigma, RMS, Bias, NMAD, Note]
rows = [
    ['MAIN RESULT', 'Theil-Sen single-var\n(reused standing formula)', '206', '0.707', '0.895', '0.91', '0.15', '0.818',
     'Rounds to the paper\'s existing 0.71; better\nbias/sigma/RMSE than OT and the standing\nemcee baseline. No Daume effect (dead code).'],
    ['ABLATION', 'No-plateau (6-var):\nAlpha/Beta/Fa/Ta dropped entirely', '204', '0.644', '0.977', '0.99', '0.20', '0.947',
     'Shows the model still predicts without the\nplateau-fit parameters, at a real but modest\ncost in both r and NMAD vs. baseline.'],
    ['BEST BASIC FIT\n(no opt cols)', 'Multivariate emcee\n(reused standing formula)', '202', '0.714', '0.901', '0.92', '0.19', '0.795',
     'Best result of any scheme with no domain-\nadaptation columns at all -- ties/edges the\npaper\'s own baseline.'],
]

fig, ax = plt.subplots(figsize=(14.5, 6.5), facecolor='white')
ax.axis('off')
col_labels = ['Role', 'Scheme', 'N', 'r(z)', 'Sigma', 'RMS', 'Bias', 'NMAD', 'Note']
col_widths = [0.12, 0.20, 0.05, 0.06, 0.06, 0.06, 0.06, 0.06, 0.33]

tbl = ax.table(cellText=rows, colLabels=col_labels, cellLoc='center', loc='center',
               colWidths=col_widths, bbox=[0.0, 0.0, 1.0, 1.0])
tbl.auto_set_font_size(False)
tbl.set_fontsize(10.5)

for (r, c), cell in tbl.get_celld().items():
    cell.set_edgecolor('#d8d7d1')
    if r == 0:
        cell.set_facecolor(COLOR_HEAD)
        cell.set_text_props(color='white', fontweight='bold')
        continue
    cell.set_facecolor('#f4f3ef' if r % 2 == 0 else 'white')
    if c == 0:
        cell.set_text_props(ha='center', fontweight='bold', color='#1a4e8c', fontsize=9.5)
    elif c == 1:
        cell.set_text_props(ha='left', fontweight='bold', color='#222')
        cell.PAD = 0.02
    elif c == 8:
        cell.set_text_props(ha='left', color='#555', fontsize=9)
    else:
        cell.set_text_props(color='#222')

tbl.scale(1, 4.2)
ax.set_title('Candidate results for the paper\n'
             '("Without 2σ outliers" subset; standing formula reused throughout, no Daume domain-adaptation effect)',
             fontsize=12.5, color='#111', pad=14)

out = f'{REPO}/paper_candidate_results/plots/summary_table.png'
plt.savefig(out, dpi=150, facecolor='white', bbox_inches='tight')
print(f'Saved {out}')
