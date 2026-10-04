from dental_release.paths import expand as _release_expand
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(_release_expand('@DENTAL_WORK_ROOT@/X37_ALIGNER_MOVEMENT'))
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
pairs = json.loads((ROOT / 'raw/R2_PAIRS.json').read_text())
joined = json.loads((ROOT / 'raw/R3V2_JOINED.json').read_text())
upper = json.loads((ROOT / 'raw/R5_POSES.json').read_text())
(fig, ax) = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
for (a, p, label) in zip(ax[0], ['3485', '6457'], ['A', 'B']):
    worst = max((r for r in pairs if r['patient'] == p), key=lambda r: r['relative_tracking_error_mm'])
    (i, k) = (worst['fdi_i'], worst['fdi_j'])
    rr = sorted([r for r in pairs if r['patient'] == p and r['fdi_i'] == i and (r['fdi_j'] == k)], key=lambda r: r['week'])
    t = [0] + [r['week'] for r in rr]
    o = [0] + [r['observed_cumulative_separation_change_mm'] for r in rr]
    q = [0] + [r['planned_cumulative_separation_change_mm'] for r in rr]
    u = np.array([0.8] + [r['tracking_sensitivity_budget_mm'] for r in rr])
    a.fill_between(t, np.array(q) - u, np.array(q) + u, alpha=0.18, color='orange', label='Declared sensitivity budget (not CI)')
    a.plot(t, o, 'o-', label='Measured surfaces')
    a.plot(t, q, 's--', color='darkorange', label='Plan')
    a.set(title=f'Patient {label}, crown pair {i}–{k}', xlabel='Week', ylabel='Change in center separation [mm]')
    a.legend(fontsize=8)
for (p, color, label) in [('3485', '#285a96', 'A'), ('6457', '#b04337', 'B')]:
    rr = [r for r in joined if r['patient'] == p and r['scenario'] == 'formed334_zero_gap' and r['numerical_admissible']]
    ax[1, 0].scatter([r['force_norm_N'] for r in rr], [r['observed_weekly_mm'] for r in rr], s=12, alpha=0.5, label=f'Patient {label}', color=color)
ax[1, 0].set(xlabel='Conditional force scenario [N]', ylabel='Reference-dependent weekly motion [mm]', title='Force port: scenarios, no patient force validation')
ax[1, 0].legend()
for (p, color, label) in [('3485', '#285a96', 'A'), ('6457', '#b04337', 'B')]:
    rr = [r for r in upper if r['patient'] == p]
    ax[1, 1].scatter([r['week'] for r in rr], [r['reference_split_error_mm'] for r in rr], s=13, alpha=0.5, label=f'Patient {label}', color=color)
ax[1, 1].axhline(0.1, color='black', ls='--', label='Frozen split gate')
ax[1, 1].set(xlabel='Week', ylabel='Crown-center reference variation [mm]', title='Low surface residual does not fix the tooth frame')
ax[1, 1].legend(fontsize=8)
fig.suptitle('X37: n=2 pipeline test — pair tracking identifiable, absolute movement/force response UNKNOWN\nTop pairs selected after measurement for display; repeated teeth/weeks are dependent.', fontsize=12)
fig.savefig(ROOT / 'figures/chain_test.png', dpi=170)
fig.savefig(ROOT / 'figures/chain_test.pdf')
plt.close(fig)
(fig, axs) = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
for (a, p, label) in zip(axs, ['3485', '6457'], ['A', 'B']):
    d = np.load(CACHE / f'{p}_full_upper_T0_100k.npz')
    pts = d['points']
    clouds = {int(q.stem.split('_Z')[-1]): np.load(q) for q in CACHE.glob(f'{p}_Sirona_T0_Z*.npz') if int(q.stem.split('_Z')[-1]) // 10 in [1, 2]}
    cp = np.concatenate([d['points'][:4000] for d in clouds.values()])
    (dd, _) = cKDTree(cp).query(pts)
    zi = np.mean([d['centroid'][2] for (f, d) in clouds.items() if f % 10 == 1])
    zm = np.mean([d['centroid'][2] for (f, d) in clouds.items() if f % 10 == 6])
    lo = min(zm + 3, zi - 5)
    hi = max(zm + 3, zi - 5)
    r = pts[(np.abs(pts[:, 0]) <= 8) & (pts[:, 2] >= lo) & (pts[:, 2] <= hi) & (dd >= 3)]
    a.scatter(pts[::10, 0], pts[::10, 2], s=0.3, c='0.7')
    a.scatter(r[::2, 0], r[::2, 2], s=1, c='#2d7890')
    a.set(title=f'Patient {label}: untrimmed scan and candidate inner ROI', xlabel='Occlusal-frame x [mm]', ylabel='Occlusal-frame z [mm]', aspect='equal')
fig.suptitle('Candidate non-dental region is measured; third-ruga identity and stability remain unverified')
fig.savefig(ROOT / 'figures/reference_roi.png', dpi=170)
plt.close(fig)
