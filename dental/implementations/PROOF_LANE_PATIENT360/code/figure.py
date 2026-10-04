from dental_release.paths import expand as _release_expand
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
R = Path(__file__).resolve().parents[1]
run = Path(sys.argv[1])
D = R / 'data' / run.name
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
fig = plt.figure(figsize=(15, 9), constrained_layout=True)
gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 1])
m = np.load(D / _release_expand('@DENTAL_CASE_C@/contact_map.npz'))
ax = fig.add_subplot(gs[0, 0])
take = np.arange(len(m['gap']))[::max(1, len(m['gap']) // 30000)]
s = ax.scatter(m['xy'][take, 0], m['xy'][take, 1], c=np.clip(m['gap'][take], -1, 1), s=3, cmap='coolwarm', vmin=-1, vmax=1, rasterized=True)
ax.set_aspect('equal')
ax.set_title(_release_expand('@DENTAL_SURFACE_ID_C@ · same registered bite\n37 tooth pairs; FDI from X11'))
ax.set_xlabel('Shared local x [mm]')
ax.set_ylabel('Shared local y [mm]')
fig.colorbar(s, ax=ax, label='Projected gap [mm], clipped ±1\nPER_POINT; no contact pressure', shrink=0.68)
for (i, case) in enumerate([_release_expand('@DENTAL_CASE_A@'), _release_expand('@DENTAL_CASE_B@')], 1):
    ax = fig.add_subplot(gs[0, i], projection='3d')
    q = np.load(D / case / 'figure_geometry.npz')
    for (key, color, label, size, alpha) in [('tooth_xyz', '#bac6d4', 'Tooth36', 2, 0.14), ('pulp_xyz', '#d44730', 'Pulp', 8, 0.7), ('canal_xyz', '#c89810', 'Canal', 3, 0.35), ('sinus_xyz', '#43a0b5', 'Sinus', 2, 0.15)]:
        a = q[key]
        if len(a):
            ax.scatter(*a.T, s=size, c=color, alpha=alpha, label=label, rasterized=True)
    ax.set_title(case + '\nExactly the same CT as Pulpy3D')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')
    ax.view_init(21, -65)
    ax.legend(fontsize=7, loc='upper right')
    ax.text2D(0.01, -0.04, 'PER_POINT · labels, not histology', transform=ax.transAxes, fontsize=8)
r3 = json.loads((run / 'R3/RESULTS_R3.json').read_text())
rows = r3['theta_rows']
x = [r['shared_theta']['IOS_pose_delta_mm'] * 1000 for r in rows]
ax = fig.add_subplot(gs[1, 0])
ax.plot(x, [r['whole_bite']['absolute_near_area_mm2'] for r in rows], 'o-', label='Whole bite · PER_ARCH')
ax.plot(x, [r['crown']['contact_area_mm2'] for r in rows], 's-', label='Crown36 · PER_SURFACE_REGION')
ax.annotate('No load patch -> FE abstains', xy=(50, 0), xytext=(-40, 1.3), arrowprops={'arrowstyle': '->'}, fontsize=8)
ax.set_xlabel('ONE shared pose shift [micrometres]')
ax.set_ylabel('Geometric proximity area [mm^2]')
ax.set_title('Same theta changes several answers')
ax.legend(fontsize=8)
ax.grid(alpha=0.18)
r2 = json.loads((run / 'R2/RESULTS_R2.json').read_text())
ax = fig.add_subplot(gs[1, 1])
for t in r2['thermal']:
    if t['status'] == 'SIMULATED_CONDITIONAL':
        ax.plot([s['h_W_m2K'] for s in t['scenarios']], [s['peak_pulp_delta_C'] for s in t['scenarios']], 'o-', label=t['patient_id'].replace('ToothFairy2', ''))
ax.set_xlabel('Assumed cooling h [W/m²K]')
ax.set_ylabel('Maximum pulp temperature increase [degrees C]')
ax.set_title('Patient geometry, assumed heat source\nPER_POINT; physical accuracy UNKNOWN')
ax.legend()
ax.grid(alpha=0.18)
ax = fig.add_subplot(gs[1, 2])
labels = ['IOS / report', 'Pulp / canal / bone', 'Crown + design gate', 'Heat simulation', 'Export with sidecar', 'Individual clinical prediction']
matrix = np.array([[2, 0, 0], [0, 2, 2], [1, 0, 0], [0, 1, 1], [2, 0, 0], [0, 0, 0]])
ax.imshow(matrix, cmap=ListedColormap(['#d9dde2', '#e9bd65', '#60aa9a']), vmin=0, vmax=2, aspect='auto')
ax.set_xticks(range(3), [_release_expand('@DENTAL_SURFACE_ID_C@'), _release_expand('@DENTAL_VOLUME_ID_A@'), _release_expand('@DENTAL_VOLUME_ID_B@')])
ax.set_yticks(range(6), labels)
ax.set_title('Executed breadth with preserved stops')
ax.tick_params(axis='both', length=0)
for j in range(6):
    for i in range(3):
        ax.text(i, j, ['Abstain', 'Conditional/error', 'Executed'][matrix[j, i]], ha='center', va='center', fontsize=8)
fig.suptitle('PATIENT360 · three real cases, two separate patient tracks\nNo verified IOS-CBCT link between individuals; no clinical recommendation', fontsize=15)
fig.savefig(R / 'PATIENT360.png', dpi=170)
fig.savefig(R / 'PATIENT360.pdf')
plt.close(fig)
