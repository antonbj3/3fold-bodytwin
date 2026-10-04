import os
os.environ['MPLCONFIGDIR'] = str(__import__('pathlib').Path(__file__).resolve().parents[1] / '.mplconfig')
from pathlib import Path
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': 'white'})
a = json.load(open(R / 'raw/RESULTS_R1.json'))
b = json.load(open(R / 'raw/RESULTS_R2.json'))
c = json.load(open(R / 'raw/RESULTS_R3.json'))
g = json.load(open(R / 'inputs/ARCH_GEOMETRY.json'))
(fig, axs) = plt.subplots(2, 2, figsize=(13, 9.0))
ax = axs[0, 0]
row = a['forces_by_type_external']
pos = np.arange(3)
w = 0.24
ax.bar(pos - w, [r['nominal_median_N'] for r in row], w, label='Nominal 0.75 mm, E=2746 MPa', color='#aeb7c3')
ax.bar(pos, [r['formed_median_N'] for r in row], w, label='Formed literature h and E', color='#3d91b6')
ax.bar(pos + w, [r['force_N'] for r in row], w, yerr=[r['force_SD_N'] for r in row], label='External sensor mean ± SD', color='#c05c40', capsize=3)
ax.set_xticks(pos, ['Central incisor', 'Canine', 'Second premolar'])
ax.set_ylabel('Buccolingual force (N)')
ax.set_title('A  Real-arch ribbon simulation fails force transfer\n0.20 mm; different specimens and PETG brand unconfirmed', fontsize=11)
ax.legend(fontsize=8)
ax = axs[0, 1]
from series import stiffness
p = json.load(open(R / 'raw/PREDICTIONS_R2.json'))
cal = json.load(open(R / 'inputs/R2_CALIBRATION_ENDPOINTS.json'))
for (direct, col) in [('labial', '#487d9b'), ('palatal', '#c05c40')]:
    s = np.linspace(0.3, 0.75, 120)
    ax.plot(s, [stiffness(p['fits'][direct], t) for t in s], color=col, label=direct.capitalize() + ' series fit')
    ep = cal['calibration'][direct]
    ax.scatter([r['sheet_mm'] for r in ep], [r['K_N_per_mm'] for r in ep], s=42, color=col, marker='s')
    vv = [r for r in b['rows'] if r['direction'] == direct]
    ax.scatter([r['sheet_mm'] for r in vv], [r['measured_K_N_per_mm'] for r in vv], s=45, facecolors='white', edgecolors=col, linewidths=1.4)
ax.set_xlabel('Original sheet thickness (mm), not measured formed h')
ax.set_ylabel('Sensor stiffness (N/mm)')
ax.set_title('B  Two measured endpoints predict an unfit 0.625 mm group\nPrimary errors 0.8% / 6.4%; thinner-group max error 30.4%', fontsize=11)
ax.legend(fontsize=8)
ax.text(0.03, 0.04, 'Squares: fit endpoints; circles: unfit published groups\nElkholy 2016 Table 2; 37°C, saliva, 30 N securing load', transform=ax.transAxes, fontsize=8)
ax = axs[1, 0]
rr = c['groups']
p0 = np.arange(2)
ax.bar(p0 - w, [r['predicted_N'] for r in rr['between_study']['rows']], w, label='Frozen source-conditioned prediction', color='#3d91b6')
ax.bar(p0, [r['force_N'] for r in rr['within_study']['rows']], w, label='2016 median, Table 1', color='#88758f')
ax.bar(p0 + w, [r['force_N'] for r in rr['between_study']['rows']], w, yerr=[r['SD_N'] for r in rr['between_study']['rows']], label='2015 independent mean ± SD', color='#c05c40', capsize=3)
ax.set_xticks(p0, ['0.50 mm sheet', '0.75 mm sheet'])
ax.set_ylabel('Palatal force at 0.25 mm (N)')
ax.set_title('C  Calibration does not transfer between published fixtures\nAbsolute force errors 150% / 144%', fontsize=11)
ax.legend(fontsize=8)
ax = axs[1, 1]
arch = next((r for r in g['arches'] if r['dataset'] == 'Teeth3DS' and r['jaw'] == 'upper'))
teeth = arch['teeth']
fr = arch['frame']
U = np.array(fr['right_unit'])
V = np.array(fr['anterior_unit'])
cent = np.array(fr['center_mm'])
pr = json.load(open(R / 'raw/PREDICTIONS_R1.json'))
rs = {r['active_fdi']: r for r in pr['rows'] if r['case'] == arch['case'] and r['jaw'] == 'upper' and (r['activation_mm'] == 0.2) and (r['status'] == 'CONDITIONAL_RIBBON_MODEL')}
cs = [(np.array(t['centroid_mm']) - cent) @ [U, V] for t in []]
for t in teeth.values():
    xyz = np.array(t['centroid_mm']) - cent
    xx = xyz @ U
    yy = xyz @ V
    k = t['fdi']
    rrr = rs.get(k)
    F = rrr['formed']['active_force_norm_N'] if rrr else 0
    ax.scatter(xx, yy, c=[F], vmin=0, vmax=2, cmap='viridis', s=220, edgecolors='white')
    ax.text(xx, yy, str(k), ha='center', va='center', color='white', fontsize=8)
sm = plt.cm.ScalarMappable(norm=plt.Normalize(0, 2), cmap='viridis')
fig.colorbar(sm, ax=ax, label='Conditional simulated force (N)', shrink=0.8)
ax.set_xlabel('Arch lateral coordinate (mm)')
ax.set_ylabel('Arch anterior coordinate (mm)')
ax.set_aspect('equal')
ax.set_title('D  Source-labelled real upper arch; one tooth moved at a time\nNo specimen-matched thickness or force validation', fontsize=11)
fig.suptitle('X19 — Formed thickness matters, but contact-state calibration is still required', fontsize=15, y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(R / 'figures/aligner_force.png', dpi=180)
fig.savefig(R / 'figures/aligner_force.pdf')
plt.close(fig)
print('Saved figure PNG/PDF')
