from dental_release.paths import expand as _release_expand
import pathlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/dental_sol_night/X69'))
r = json.loads((ROOT / 'raw/R2_RESULTS.json').read_text())
r3 = json.loads((ROOT / 'raw/R3_RESULTS.json').read_text())
w = np.load(DATA / 'figure_r2.npz')
fig = plt.figure(figsize=(14, 9))
ax = fig.add_subplot(221, projection='3d')
q = w['target']
ax.scatter(*q.T, s=0.5, c='#999999', alpha=0.35)
for jaw in ['lower', 'upper']:
    p = w[jaw + '_transformed_points'][::8]
    d = w[jaw + '_residual_mm'][::8]
    s = ax.scatter(*p.T, c=d, s=3, vmin=0, vmax=0.6, cmap='turbo')
ax.set(xlabel='CBCT native x (source mm)', ylabel='y', zlabel='z', title='Source-paired Demo_1: CBCT-only tooth mesh + IOS crowns\nSeparate jaw poses; metric scale assumed')
fig.colorbar(s, ax=ax, label='Nearest-triangle discrepancy (mm)', shrink=0.55)
ax = fig.add_subplot(222)
tooths = [x for jaw in r['jaws'] for x in jaw['tooth_results']]
names = [str(x['fdi']) for x in tooths]
v = [x['p95_mm'] for x in tooths]
ax.bar(names, v, color=['#d45232' if x['heldout_from_fit'] else '#207f9b' for x in tooths])
ax.axhline(0.5, color='black', ls='--', label='Frozen0.5 mm research gate')
ax.tick_params(axis='x', rotation=90)
ax.set(ylabel='p95 discrepancy (mm)', title='PER_TOOTH: orange teeth excluded from fit')
ax.legend(fontsize=8)
ax = fig.add_subplot(223)
for jaw in ['lower', 'upper']:
    d = np.sort(w[jaw + '_residual_mm'][np.isin(w[jaw + '_fdi'], [11, 21, 31, 41])])
    ax.plot(d, np.linspace(0, 1, len(d)), label=jaw)
ax.axvline(0.5, c='black', ls='--')
ax.set(xlim=(0, 1), xlabel='Heldout surface discrepancy (mm)', ylabel='Sample fraction', title='Four heldout incisors: consistency, not anatomical TRE')
ax.legend()
ax = fig.add_subplot(224)
ax.axis('off')
(A, B) = r3['designs']
ax.text(0, 1, f"What became possible\nReal IOS crowns in CBCT crown/root coordinate frame.\nHeldout p95: lower0.304 mm; upper0.311 mm.\nSource: 10.5281/zenodo.8027553 (CC BY4.0).\n\nWhat remains unknown\nIndependent pose/fiducial trueness; STL metric calibration.\nCBCT segmentation accuracy; force, material and bite state.\nNearest-surface witnesses are not independent landmarks.\n\nRoot-query conditional enclosure: error <= epsilon x L(x)\nCrown-centroid anchors: median L={A['amplification_median']:.2f}\nSpread crown anchors: median L={B['amplification_median']:.2f}\nepsilon UNKNOWN; exact digital bound, no clinical claim.\n\nX69 | PENDING INDEPENDENT REVIEW", va='top', fontsize=10)
fig.tight_layout()
fig.savefig(ROOT / 'PAIR_REGISTRATION.png', dpi=160)
