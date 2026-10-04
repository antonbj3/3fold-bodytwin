from dental_release.paths import expand as _release_expand
import pathlib, numpy as np, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_XBREAK_HUNT_5'))
fig = plt.figure(figsize=(14, 9), layout='constrained')
for (k, case) in enumerate(['6158-11', '6138-21']):
    z = np.load(DATA / (case + '_R1.npz'))
    v = z['v']
    f = z['f']
    lab = z['lab']
    ax = fig.add_subplot(2, 2, k + 1, projection='3d')
    for (label, color, alpha) in [(0, '#bcbcbc', 0.3), (1, '#376ca6', 0.85)]:
        facets = f[lab == label][::5]
        ax.add_collection3d(Poly3DCollection(v[facets], facecolor=color, edgecolor='none', alpha=alpha))
    lo = v.min(axis=0)
    hi = v.max(axis=0)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(hi - lo)
    ax.set_title(case + ' — original labelled die')
    ax.set_xlabel('source x')
    ax.set_ylabel('source y')
    ax.set_zlabel('source z')
r3 = json.loads((ROOT / 'raw/R3.json').read_text())
ax = fig.add_subplot(2, 2, 3)
for r in r3['rows']:
    h = [p['h_source_unit'] for p in r['grids']]
    lo = [p['gap']['lower_max'] for p in r['grids']]
    hi = [p['gap']['upper_max'] for p in r['grids']]
    ax.plot(h, lo, 'o-', label=r['case'] + ' lower witness')
    ax.plot(h, hi, 'x--', label=r['case'] + ' continuous upper bound')
ax.axhline(0.12, color='#ae3535', label='fixed geometry scenario .12')
ax.set_xlabel('Conservative cell width (source units)')
ax.set_ylabel('Euclidean blockout gap (source units)')
ax.legend(fontsize=7)
ax.set_ylim(0, 0.55)
ax.grid(alpha=0.2)
ax = fig.add_subplot(2, 2, 4)
ax.axis('off')
ax.text(0, 1, 'What this establishes\n\n11 measured die files imported; 2 reported files absent.\n0/11 pass an unsmoothed normal-cone test.\nA conservative cavity is exported for two pilots.\n6158-11 passes whole-facet .12 gap and opening scenarios.\n6138-21 fails gap; source-margin projection self-intersects.\n\nUnit calibration, true margin, cement, as-built and\na complete crown remain UNKNOWN.\n\nBlue = original Label 1. Rendering uses every fifth facet;\nall numerical tests use every original top facet.', va='top', fontsize=10)
fig.suptitle('Measured preparation → directional cavity → source margin', fontsize=15)
(ROOT / 'figures').mkdir(exist_ok=True)
fig.savefig(ROOT / 'figures/CONSTRUCTION.png', dpi=160)
fig.savefig(ROOT / 'figures/CONSTRUCTION.pdf')
