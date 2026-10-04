"""System Python (-s) is used for Matplotlib ABI compatibility; no mesh dependency."""
from pathlib import Path
import json, os
ROOT = Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR'] = str(ROOT / 'raw/mpl-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r = json.loads((ROOT / 'rounds/R2.json').read_text())
s = json.loads((ROOT / 'raw/SUFFICIENCY_R1.json').read_text())
fig = plt.figure(figsize=(13, 4.4))
a = fig.add_subplot(131)
b = fig.add_subplot(132)
c = fig.add_subplot(133, projection='3d')
rows = [x for x in r['bond_lab_screens'] if x['with_cover']]
for (i, x) in enumerate(rows):
    mu = x['mean_MPa']
    (lo, hi) = x['pointwise_95pct_mean_CI_MPa']
    a.errorbar(i, mu, yerr=[[mu - lo], [hi - mu]], fmt='o', capsize=5, color='#26744d' if x['lab_mean_screen'] == 'PASS' else '#ad3636')
a.axhline(5, color='black', ls='--', label='Declared lab challenge, 5 MPa')
a.set_xticks(range(3))
a.set_xticklabels([x['system'] for x in rows])
a.set_ylabel('Nominal SBS [MPa], population mean')
a.set_title('Cement + adhesive + 1 mm cover')
a.legend(fontsize=7, loc='lower left')
a.set_ylim(0, 11)
for (j, substrate) in enumerate(['ZR', 'CR', 'MT']):
    rows = [x for x in r['optical_lab_screens'] if x['substrate'] == substrate and x['cement'] != 'CG']
    for (i, x) in enumerate(rows):
        mu = x['de00_of_mean_Lab']
        (lo, hi) = x['rigorous_formula_enclosure_de00']
        pos = i + (j - 1) * 0.18
        b.errorbar(pos, mu, yerr=[[mu - lo], [hi - mu]], fmt='o', capsize=3, color=['#6f68b0', '#4e9655', '#b96136'][j], label=substrate if i == 0 else None)
b.axhline(3.6, color='black', ls='--')
b.set_xticks(range(4))
b.set_xticklabels(['TR', 'BL', 'OP', 'A2'])
b.set_ylabel('CIEDE2000 to source A1')
b.set_title('1 mm measured stack\nBoxes: Lab mean ±2 SD, no probability')
b.legend(fontsize=8)
d = np.load(ROOT / 'raw/stress.npz')
xyz = d['cf'][::12]
c.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=1, color='#c6c6c6', alpha=0.55)
points = np.array(s['bond']['positions_mm'])
c.scatter(points[:, 0], points[:, 1], points[:, 2], s=50, c=['#af3535', '#337faa'])
c.set_xlabel('x [mm]')
c.set_ylabel('y [mm]')
c.set_zlabel('z [mm]')
c.set_title('Same STS D1 crown, simulated field\nSame mean strength: local screen ×2.48')
c.view_init(elev=22, azim=-48)
fig.suptitle('Two restorative information links; patient physical transfer remains UNKNOWN', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.91])
(ROOT / 'figures').mkdir(exist_ok=True)
fig.savefig(ROOT / 'figures/restorative_demo.png', dpi=180)
fig.savefig(ROOT / 'figures/restorative_demo.pdf')
print('Figure generated.')
