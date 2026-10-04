from pathlib import Path
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = Path(__file__).resolve().parent
OUT = Path(os.environ.get('X24_OUTPUT_DIR', str(P)))
r1 = json.loads((OUT / 'results_R1.json').read_text())
r2 = json.loads((OUT / 'results_R2.json').read_text())
r3 = json.loads((OUT / 'results_R3.json').read_text())
b = json.loads((P / 'PHANTOM_BASELINE.json').read_text())
ref = {r['material']: r['gray_median'] for r in next((d for d in b if d['device'] == 'Varian'))['inserts']}
(fig, ax) = plt.subplots(2, 2, figsize=(12, 8.3), constrained_layout=True)
colors = ['#cc783b', '#2d8190', '#946bb0', '#425a73']
for (c, d) in zip(colors, b):
    x = [ref[r['material']] for r in d['inserts']]
    y = [r['gray_median'] for r in d['inserts']]
    ax[0, 0].plot(x, y, 'o-', label=d['device'].replace('Elekta_', ''), c=c)
ax[0, 0].plot([-1000, 1200], [-1000, 1200], '--', c='0.5', lw=1)
ax[0, 0].set(xlabel='Measured Varian reference [CT -tal]', ylabel='Measured grey value', title='A. Same Catphan material, different scanners')
ax[0, 0].legend(fontsize=8)
vv = r2['variants']
names = list(vv)
pos = np.arange(len(names))
ax[0, 1].bar(pos, [vv[v]['max_abs_HU'] for v in names], color=['#2d8190', '#946bb0', '#cc783b', '#6e7d89'])
ax[0, 1].axhline(40, color='#b64141', linestyle='--', label='Frozen gate 40 HU_ref')
ax[0, 1].set(xticks=pos, xticklabels=['3 kvadratiska', '3 styckvisa', '2 luft/Delrin', '2 luft/PS'], ylabel='Maximum holding errors [HU_ref ]', title='B. Material + 18 new scans')
ax[0, 1].legend(fontsize=8)
ax[0, 1].tick_params(axis='x', labelsize=8)
devs = list(r3['summary_by_device'])
sets = [[1 / a['E_naive_ratio'] for a in r3['answers'] if a['device'] == d] for d in devs]
ax[1, 0].boxplot(sets, labels=['LinacB', 'LinacD', 'LinacE'], showfliers=True)
ax[1, 0].set(yscale='log', ylabel='E-reference / E-raw-HU [×]', title='C. Conditional E-sensitivity, 34 sites')
ax[1, 0].text(0.02, 0.95, 'Density and E-teams are Closures.\nNo measured patient-E.', transform=ax[1, 0].transAxes, va='top', fontsize=9)
vals = [r3['summary_by_device'][d] for d in devs]
pp = np.arange(3)
ax[1, 1].bar(pp - 0.18, [v['u50_predicate_flips'] for v in vals], 0.36, label='50 µm - gate changed', color='#cc783b')
ax[1, 1].bar(pp + 0.18, [v['naive_density_invalid'] for v in vals], 0.36, label='Raw HU provides ρ_app ≤ 0', color='#2d8190')
ax[1, 1].set(xticks=pp, xticklabels=['LinacB', 'LinacD', 'LinacE'], ylabel='Sites by 34', title='D. Conditional mechanics; K42 abstains')
ax[1, 1].legend(fontsize=8)
for a in ax.ravel():
    a.grid(alpha=0.18)
    a.set_axisbelow(True)
fig.suptitle('X24 : Scanner information alters material base\nExternal Sensitometri; implant scenarios have no parade mechanical validation', fontsize=13)
fig.savefig(OUT / 'X24_CALIBRATION_FIGURE.png', dpi=170)
fig.savefig(OUT / 'X24_CALIBRATION_FIGURE.svg')
print(OUT / 'X24_CALIBRATION_FIGURE.png')
