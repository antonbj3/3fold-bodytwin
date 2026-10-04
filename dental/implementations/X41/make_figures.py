"""Publication figures from frozen, replayed evidence. No simulated anatomy art."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parent
import platform
(ROOT / 'FIGURE_ENVIRONMENT.json').write_text(json.dumps({'python': platform.python_version(), 'numpy': np.__version__, 'matplotlib': matplotlib.__version__, 'PYTHONNOUSERSITE': '1', 'backend': 'Agg', 'threads': 1}, indent=2) + '\n')
plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
d = json.loads((ROOT / 'results.json').read_text())
rows = d['counterexamples']
y = np.arange(len(rows))
(fig, ax) = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={'width_ratios': [1.3, 1]})
for (i, r) in enumerate(rows):
    ax[0].plot([r['adverse_gap_mm'], r['reference_gap_mm']], [i, i], color='#666666', lw=1)
ax[0].scatter([r['reference_gap_mm'] for r in rows], y, label='Reference / location control', color='#2166ac', s=50, zorder=3)
ax[0].scatter([r['adverse_gap_mm'] for r in rows], y, label='Attached adverse addition', color='#b2182b', s=50, zorder=3)
ax[0].axvline(2, color='black', ls='--', lw=1, label='Declared threshold')
ax[0].set_yticks(y, [r['case'].replace('ToothFairy2', '') + f" / cc{r['component']}" for r in rows])
ax[0].invert_yaxis()
ax[0].set_xlabel('Complete-solid cube clearance (mm)')
ax[0].set_title('a  Decision changes under local additions', loc='left')
ax[0].legend(loc='lower right', fontsize=8)
ax[1].barh(y, [r['minimum_added_voxels'] for r in rows], color='#b2182b', height=0.65)
for (i, r) in enumerate(rows):
    ax[1].text(r['minimum_added_voxels'] + 0.15, i, f"Dice {100 * r['Dice']:.4f}%", va='center', fontsize=9)
ax[1].set_yticks([])
ax[1].invert_yaxis()
ax[1].set_xlim(0, 16)
ax[1].set_xlabel('Minimum added voxels (six-neighbour lattice)')
ax[1].set_title('b  Same Dice and HD95 in each pair', loc='left')
fig.text(0.5, 0.01, 'HD95 = 0 mm for both members of every pair. Six selected injections; no model-error prevalence or physical safety inference.', ha='center', fontsize=8)
fig.tight_layout(rect=(0, 0.045, 1, 1))
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(ROOT / 'figures' / ('FIGURE_1_DECISION_ERRORS.' + ext), dpi=220)
plt.close(fig)
raw = json.loads((ROOT / 'inputs/LANE_X26_DECIDABILITY/RAW_R2_CANAL.json').read_text())
floor = d['conditional_canal']['floor_mm']
margins = np.array([abs(r['margin']) for r in raw])
(fig, ax) = plt.subplots(figsize=(8.5, 4.3))
bins = np.linspace(0, 15, 61)
ax.hist([margins[margins <= floor], margins[margins > floor]], bins=bins, stacked=True, color=['#b2182b', '#2166ac'], label=['Class 2: no finite pitch in scenario', 'Class 1: finite conditional pitch'])
ax.axvline(floor, color='#b2182b', ls='--', label=f'Scenario floor = {floor:.3f} mm')
ax.axvline(floor + np.sqrt(3) * 0.3, color='black', ls=':', label='Floor + numerical bound at h = 0.3 mm')
ax.set_xlabel('Absolute signed minimum-clearance margin |d − 2 mm| (mm)')
ax.set_ylabel('Dependent decision queries (count)')
ax.legend(fontsize=8)
ax.set_title('Numerical reachability under an unvalidated physical-error scenario', loc='left')
fig.text(0.5, 0.015, '4,132 queries in 338 archive cases. No surface-level confidence or anatomical certificate; resampling does not improve acquisition.', ha='center', fontsize=8)
fig.tight_layout(rect=(0, 0.055, 1, 1))
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(ROOT / 'figures' / ('FIGURE_2_RESOLUTION_FLOOR.' + ext), dpi=220)
plt.close(fig)
