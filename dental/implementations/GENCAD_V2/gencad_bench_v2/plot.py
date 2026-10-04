"""System-only plotting process avoids the local NumPy/Matplotlib ABI mismatch."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
R = Path(__file__).resolve().parents[1]
d = json.loads((R / 'raw/FIGURE_DATA.json').read_text())
fig = plt.figure(figsize=(15, 9), constrained_layout=True)
ax = fig.add_subplot(221)
names = list(d['totals'])
left = np.zeros(len(names))
colors = {'PASS': '#19836b', 'FAIL': '#c45148', 'ABSTAIN': '#d59b32', 'UNKNOWN': '#8d94a5', 'UNKNOWN_SITE': '#b1b6c1', 'INVALID': '#4f465c'}
for (status, color) in colors.items():
    v = np.array([d['totals'][n].get(status, 0) for n in names])
    ax.barh(names, v, left=left, label=status, color=color)
    left += v
ax.set_xlabel('Requested tasks (576 per participant)')
ax.set_title('Scoped digital outcome; no composite score')
ax.legend(fontsize=8, ncol=3)
ax = fig.add_subplot(222)
families = list(dict.fromkeys((r['family'] for r in d['summaries'])))
matrix = []
for name in names:
    matrix.append([sum((r['counts'].get('PASS', 0) for r in d['summaries'] if r['family'] == f and r['participant'] == name)) / 64 for f in families])
im = ax.imshow(matrix, vmin=0, vmax=1, cmap='YlGnBu', aspect='auto')
ax.set_yticks(range(len(names)), names)
ax.set_xticks(range(len(families)), families, rotation=40, ha='right', fontsize=8)
ax.set_title('Digital PASS fraction by family')
fig.colorbar(im, ax=ax)
ax = fig.add_subplot(223)
regions = [r for r in ['marginal', 'axial', 'occlusal'] if r in d['measurements']]
x = np.arange(len(regions))
ax.bar(x - 0.18, [d['measurements'][r]['practice']['mae'] for r in regions], 0.36, label='CAD = measured gap', color='#98a5c3')
ax.bar(x + 0.18, [d['measurements'][r]['candidate']['mae'] for r in regions], 0.36, label='Held-study regional correction', color='#437c80')
ax.axhline(20, color='#bc4944', ls='--', label='Frozen 20 µm gate')
ax.set_xticks(x, regions)
ax.set_ylabel('Mean absolute error, µm')
ax.legend(fontsize=8)
ax.set_title('External cement measurements; transfer can fail')
ex = d['example']
ax = fig.add_subplot(224, projection='3d')
if ex:
    v = np.array(ex['outer'])
    f = np.array(ex['faces'])
    t = np.array(ex['antagonist'])
    ax.add_collection3d(Poly3DCollection(v[f], facecolor='#2f9a91', edgecolor='#376961', linewidth=0.2, alpha=0.9))
    if len(t):
        ax.add_collection3d(Poly3DCollection(t[::max(1, len(t) // 2000)], facecolor='#c8c8c8', edgecolor='none', alpha=0.15))
    points = np.concatenate([v, np.array(ex['inner'])])
    lo = points.min(0)
    hi = points.max(0)
    ax.set_xlim(lo[0] - 1, hi[0] + 1)
    ax.set_ylim(lo[1] - 1, hi[1] + 1)
    ax.set_zlim(lo[2] - 0.5, hi[2] + 2)
    ax.set_title('Generated roof in measured antagonist scene')
else:
    ax.text2D(0.05, 0.5, 'No fully passing molar example', transform=ax.transAxes)
ax.set_xlabel('mm')
ax.set_ylabel('mm')
ax.set_zlabel('mm')
fig.suptitle('DentalGenCAD-Bench v2 — research surface tasks, physical validity UNKNOWN', fontsize=15)
(R / 'figures').mkdir(exist_ok=True)
fig.savefig(R / 'figures/benchmark_v2.png', dpi=170)
fig.savefig(R / 'figures/benchmark_v2.pdf')
plt.close(fig)
if (R / 'raw/AXIAL_FIGURE_DATA.json').exists():
    d = json.loads((R / 'raw/AXIAL_FIGURE_DATA.json').read_text())
    fig = plt.figure(figsize=(11, 5), constrained_layout=True)
    for (i, (name, a)) in enumerate(d['arms'].items()):
        ax = fig.add_subplot(1, 2, i + 1, projection='3d')
        v = np.array(a['vertices'])
        f = np.array(a['faces'])
        tri = v[f]
        keep = tri[:, :, 1].mean(1) >= 0
        ax.add_collection3d(Poly3DCollection(tri[keep], facecolor='#c65c53' if name == 'unflared' else '#29927c', edgecolor='#555555', linewidth=0.25, alpha=0.8))
        inn = np.array(a['inner_vertices'])
        inf = np.array(a['inner_faces'])
        ax.add_collection3d(Poly3DCollection(inn[inf], facecolor='#d9cba1', edgecolor='#857955', linewidth=0.2, alpha=0.35))
        lo = v.min(0)
        hi = v.max(0)
        ax.set_xlim(lo[0] - 0.4, hi[0] + 0.4)
        ax.set_ylim(lo[1] - 0.4, hi[1] + 0.4)
        ax.set_zlim(lo[2] - 0.4, hi[2] + 0.4)
        ax.set_box_aspect(hi - lo)
        ax.view_init(20, -70)
        ax.set_xlabel('mm')
        ax.set_ylabel('mm')
        ax.set_zlabel('mm')
        ax.set_title(('Original axial wall: FAIL' if name == 'unflared' else 'Changed axial wall: PASS') + '\nComplete-cap exact distance, same roof/cavity')
    fig.suptitle('R3: 31/31 axial wall failures become digital passes; physical fit remains UNKNOWN')
    fig.savefig(R / 'figures/axial_crown_v2.png', dpi=170)
    fig.savefig(R / 'figures/axial_crown_v2.pdf')
    plt.close(fig)
