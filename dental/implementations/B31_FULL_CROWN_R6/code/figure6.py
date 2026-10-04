from dental_release.paths import expand as _release_expand
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
R = Path(__file__).resolve().parents[1]
D = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R6'))
P = D.parent / 'PROOF_LANE_FULL_CROWN_R5'
key = '079905ebf9504544_molar'
old = np.load(P / 'B/distance_local_thickening' / (key + '.npz'))
new = np.load(D / 'D_exterior.npz')
res = json.loads((R / 'RESULTS_D.json').read_text())
suff = json.loads((R / 'raw/SUFFICIENCY.json').read_text())

def section(tri):
    out = []
    for t in tri:
        pp = []
        for (a, b) in zip(t, np.roll(t, -1, axis=0)):
            if (a[1] < 0) != (b[1] < 0):
                u = -a[1] / (b[1] - a[1])
                pp.append((a + u * (b - a))[[0, 2]])
        if len(pp) == 2:
            out.append(pp)
    return out
(fig, ax) = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
for (tri, color, label) in [(old['vertices'][old['faces'][old['roles'] == 0]], '#a6a6a6', 'R5 outer'), (new['vertices'][new['faces']], '#1464a0', 'Repaired outer'), (old['vertices'][old['faces'][old['roles'] == 1]], '#b74738', 'R5 intaglio')]:
    ax[0].add_collection(LineCollection(section(tri), colors=color, lw=0.8, label=label))
ax[0].autoscale()
ax[0].set_aspect('equal')
ax[0].set(xlabel='x [mm]', ylabel='z [mm]', title='Priority molar: section y = 0')
ax[0].legend(fontsize=7)
x = np.arange(2)
width = 0.32
ax[1].bar(x - width / 2, [0.4367318076704878, 0.5948826396163763], width, label='R5', color='#a6a6a6')
ax[1].bar(x + width / 2, [res['shape']['p95_mm'], res['wall_against_R5_inner']['sampled_min_mm']], width, label='Repaired', color='#1464a0')
ax[1].hlines(0.437, -0.4, 0.4, color='black', linestyle='--')
ax[1].hlines(0.5, 0.6, 1.4, color='black', linestyle='--')
ax[1].set(xticks=x, xticklabels=['Shape p95\n(lower is better)', 'Sampled wall min\n(higher is better)'], ylabel='mm', title='Shape improves; wall fails')
ax[1].legend(fontsize=8)
vals = [suff['minimum_a_mm'] * 1000, suff['minimum_b_mm'] * 1000]
ax[2].bar(['State A', 'State B'], vals, color=['#1464a0', '#b74738'])
ax[2].axhline(25, color='black', ls='--', label='25 µm lab policy')
ax[2].set(ylabel='Minimum nominal gap [µm]', title='Same mean 31.25 µm exactly')
ax[2].legend(fontsize=8)
ax[2].text(0.5, 0.08, 'Mean identity error: 0\nDifferent downstream decisions', ha='center', transform=ax[2].transAxes, fontsize=9)
fig.suptitle('Research geometry: no qualified whole-crown export', fontsize=14)
(R / 'figures').mkdir(exist_ok=True)
fig.savefig(R / 'figures/demo.png', dpi=180)
fig.savefig(R / 'figures/demo.pdf')
plt.close(fig)
