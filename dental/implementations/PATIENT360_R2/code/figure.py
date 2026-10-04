"""Standalone publication/export plot; run with system NumPy/Matplotlib."""
from dental_release.paths import expand as _release_expand
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root = Path(sys.argv[1])
r = json.loads((root / 'results.json').read_text())
c = r['contact']
d = r['robust_design']
f = r['region_fields']
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
(fig, ax) = plt.subplots(2, 2, figsize=(12, 8.2), constrained_layout=True)
for row in c['partition']:
    if row['kind'] == 'OPEN_INTERVAL':
        ax[0, 0].plot(np.array([row['lower_mm'], row['upper_mm']]) * 1000, [row['contact_area_mm2']] * 2, lw=3, color='#147d92')
ax[0, 0].axvspan(c['no_contact_witness']['lower_mm'] * 1000, 50, color='#ba4944', alpha=0.2)
ax[0, 0].set(xlabel='Shared vertical pose (µm)', ylabel='Retained geometric area (mm²)', title=_release_expand('@DENTAL_CASE_ID@ tooth 36: contact support disappears'))
ax[0, 0].annotate('ABSTAIN: loaded pose missing', xy=(2, 0.1), xycoords='data', fontsize=10)
s = c['sufficient_summary_test']
ax[0, 1].bar(['Central point', 'Peripheral point'], s['downstream_tensile_MPa'], color=['#147d92', '#bd6257'])
ax[0, 1].set(ylabel='Conditional FE peak tensile stress (MPa)', title='Identical area and force, different response')
ax[0, 1].text(0.03, 0.93, 'Exact summary identity error = 0\n100 N scenario; no physical stress validation', transform=ax[0, 1].transAxes, va='top', fontsize=9)
cb = [v for v in f['descriptors'] if 'ToothFairy2P_048' in v][0]
meta = json.loads(Path(cb).read_text())
z = np.load(meta['arrays']['path'])
owner = z['owner']
sl = owner[:, :, owner.shape[2] // 2]
ax[1, 0].imshow(sl, origin='lower', cmap='viridis', vmin=0, vmax=2, interpolation='nearest', aspect='equal')
ax[1, 0].set(xlabel='Local crop y index', ylabel='Local crown-directed z index', title=_release_expand('@DENTAL_CASE_ID@: original tooth/pulp with signed crop affine'))
ax[1, 0].text(0.03, 0.98, 'PER_POINT\nSource verified\nLocal +z reverses\nCBCT array z', transform=ax[1, 0].transAxes, va='top', color='white', fontsize=8)
upper = d['nonpenetration_necessary_h_upper_mm'] * 1000
lower = d['any_point_contact_necessary_h_lower_mm'] * 1000
ax[1, 1].axvspan(-60, upper, color='#147d92', alpha=0.25, label='Necessary for no penetration')
ax[1, 1].axvspan(lower, 60, color='#bd6257', alpha=0.25, label='Necessary for any point contact')
ax[1, 1].axvline(upper, color='#147d92')
ax[1, 1].axvline(lower, color='#bd6257')
ax[1, 1].set(xlim=(-60, 60), ylim=(0, 1), yticks=[], xlabel='Rigid crown raise h (µm)', title='±50 µm pose: incompatible necessary bounds')
ax[1, 1].text(0, 0.52, f'No overlap\n{lower - upper:.3f} µm', ha='center')
ax[1, 1].legend(loc='lower center', fontsize=8)
fig.suptitle('Patient360 R2 — source geometry, shared pose, explicit abstention', fontsize=15)
fig.savefig(root / 'PATIENT360_R2.png', dpi=170)
fig.savefig(root / 'PATIENT360_R2.pdf')
plt.close(fig)
