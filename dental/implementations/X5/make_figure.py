from pathlib import Path
import csv, json, os
os.environ['MPLBACKEND'] = 'Agg'
import numpy as np
import matplotlib.pyplot as plt
P = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(P / 'DECISION_MARGINS.csv')))
r3 = json.load(open(P / 'R3_CONNECTED_COUNTEREXAMPLES.json'))
m = json.load(open(P / 'DATA_MANIFEST.json'))
(fig, ax) = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
for (reg, label, color) in [('canal_body', 'Canal, 2 mm rule', '#31688e'), ('lateral_bone', 'Bone wall, 1 mm rule', '#35b779'), ('adjacent_root', 'Adjacent root, 1.5 mm rule', '#b85b5b')]:
    v = np.sort([float(x['signed_margin_mm']) for x in rows if x['selection'] == 'ARCHIVED_RB' and x['region'] == reg])
    ax[0].plot(v, np.arange(1, len(v) + 1) / len(v), label=f'{label} (n={len(v)})', color=color)
ax[0].axvline(0.22, ls=':', color='k', label='0.22 mm scale only')
ax[0].set(xlim=(-0.2, 3), ylim=(0, 1), xlabel='Adverse edge-error margin (mm)', ylabel='Fraction of archived selected designs', title='A  Location changes tolerance')
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.2)
xx = np.arange(len(r3))
kk = [r['minimum_connected_added_voxels'] for r in r3]
ax[1].bar(xx, kk, color='#d95f02')
ax[1].set_xticks(xx)
ax[1].set_xticklabels([r['case'].replace('ToothFairy2', '') for r in r3], rotation=45, ha='right')
ax[1].set(ylabel='Minimum connected added voxels', title='B  Connected errors change all six decisions')
for (i, r) in enumerate(r3):
    ax[1].text(i, kk[i] + 0.12, f"{100 * r['Dice_adverse']:.4f}%", ha='center', fontsize=8)
ax[1].text(0.03, 0.98, 'Equal Dice and HD95=0 in each adverse/control pair', transform=ax[1].transAxes, va='top', fontsize=8)
ax[1].set_ylim(0, max(kk) + 3)
q = r3[-1]
a = np.load(m[-1]['path'])
X = a['body_voxel']
Y = a['hazard_surface_voxel']
s = a['spacing_zyx_mm']
path = np.array(q['adverse_added_voxel'])
seed = np.array(q['adverse_seed_voxel'])
allp = np.vstack([seed, path])
mid = allp.mean(0) * s
keepx = np.linalg.norm(X * s - mid, axis=1) < 5
keepy = np.linalg.norm(Y * s - mid, axis=1) < 5
delta = (path[-1] - seed) * s
i = int(np.argmax(abs(delta)))
j = max([z for z in range(3) if z != i], key=lambda z: float(np.ptp(X[:, z])))
ax[2].scatter((Y * s)[keepy, i], (Y * s)[keepy, j], s=5, alpha=0.35, color='#2686b2', label='Original canal')
ax[2].scatter((X * s)[keepx, i], (X * s)[keepx, j], s=5, alpha=0.2, color='gray', label='Fixed implant annotation')
ax[2].plot((allp * s)[:, i], (allp * s)[:, j], 'o-', color='#d95f02', ms=5, label='Minimum connected error')
ax[2].set_aspect('equal')
ax[2].set(xlabel=f'Grid axis {i} (mm)', ylabel=f'Grid axis {j} (mm)', title=f"C  {q['case'].replace('ToothFairy2', '')}: {q['minimum_connected_added_voxels']} added voxels")
ax[2].legend(fontsize=8, loc='best')
fig.suptitle('Segmentation error in decision units — annotated geometry, physical safety UNKNOWN', fontsize=13)
fig.savefig(P / 'DECISION_ERROR_FIGURE.png', dpi=180)
fig.savefig(P / 'DECISION_ERROR_FIGURE.pdf')
print('Figure saved')
