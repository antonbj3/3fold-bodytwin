"""Use distro NumPy/Matplotlib: python3 -s code/figure.py."""
from dental_release.paths import expand as _release_expand
import pathlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/dental_sol_night/X69'))
g = np.load(DATA / 'figure_r1.npz')
j = json.loads((ROOT / 'raw/R1_RESULTS.json').read_text())
fig = plt.figure(figsize=(13, 8))
ax = fig.add_subplot(221, projection='3d')
q = g['target']
p = g['source']
d = g['residual']
ax.scatter(*q.T, s=0.2, c='gray', alpha=0.2)
s = ax.scatter(*p[::8].T, c=d[::8], s=2, vmin=0, vmax=3, cmap='turbo')
ax.set_title('R1: released CBCT isosurface + transformed IOS\nSame-subject source declaration; anatomy accuracy UNKNOWN')
ax.set(xlabel='CBCT RAS x (mm)', ylabel='y (mm)', zlabel='z (mm)')
fig.colorbar(s, ax=ax, label='distance to isosurface (mm)', shrink=0.55)
ax = fig.add_subplot(222)
rows = j['regions']
names = [r['jaw'][0].upper() + str(r['native_x_tercile']) + (' H' if r['heldout_from_fit'] else ' F') for r in rows]
vals = [r['p95_mm'] for r in rows]
ax.bar(names, vals, color=['#d45a36' if r['heldout_from_fit'] else '#296f91' for r in rows])
ax.axhline(0.5, c='black', ls='--', label='Frozen 0.5 mm research gate')
ax.set(ylabel='Regional p95 discrepancy (mm)', title='PER_SURFACE_REGION; H = held out, F = fitted')
ax.legend()
ax = fig.add_subplot(223)
for jaw in ['lower', 'upper']:
    a = np.sort(d[g['jaw'] == jaw])
    ax.plot(a, np.linspace(0, 1, len(a)), label=jaw)
ax.axvline(0.5, c='black', ls='--')
ax.set(xlim=(0, 10), xlabel='Distance to released intensity isosurface (mm)', ylabel='Sample fraction', title='Residual is not target registration error')
ax.legend()
ax = fig.add_subplot(224)
ax.axis('off')
ax.text(0, 1, 'Real pair acquired: Figshare case 001\nCBCT: 0.249 x 0.249 x 0.250 mm\nIOS: metric scale assumed (STL has no units)\n\nNo independent anatomical landmarks.\nNo per-tooth ownership from threshold.\nNo uniform registration error bound.\nNo physical/clinical release.\n\nSource: 10.6084/m9.figshare.26965903.v3\nCC BY 4.0 | X69 | PENDING INDEPENDENT REVIEW', va='top', fontsize=11)
fig.tight_layout()
fig.savefig(ROOT / 'PAIR_REGISTRATION.png', dpi=160)
