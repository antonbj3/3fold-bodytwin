"""Plot in the installed, compatible system NumPy/matplotlib pair (no user site)."""
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
x = json.loads(Path(sys.argv[1]).read_text())
out = Path(sys.argv[2])
out.mkdir(exist_ok=True)
ct = np.array(x['crown_triangles'])
at = np.array(x['antagonist_triangles'])
p = np.array(x['point'])
q = np.array(x['closest'])
fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(121, projection='3d')
ax.add_collection3d(Poly3DCollection(ct, facecolor='#47879a', edgecolor='none', alpha=0.45))
ax.add_collection3d(Poly3DCollection(at, facecolor='#d1a466', edgecolor='none', alpha=0.18))
ax.plot(*np.array([p, q]).T, color='#bd2444', linewidth=4)
ax.scatter(*p, c='#bd2444', s=35)
ax.scatter(*q, c='#bd2444', s=35)
verts = np.r_[ct.reshape(-1, 3), at.reshape(-1, 3)]
lo = verts.min(0)
hi = verts.max(0)
mid = (lo + hi) / 2
width = max(hi - lo)
for (f, z) in [(ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)]:
    f(mid[z] - width / 2, mid[z] + width / 2)
ax.view_init(23, 45)
ax.set_xlabel('x [mm]')
ax.set_ylabel('y [mm]')
ax.set_zlabel('z [mm]')
ax.set_title('X18 generated roof + original antagonist\nRed segment: independent nearest-facet witness')
ax2 = fig.add_subplot(122)
vals = [1.2, x['distance_mm']]
ax2.bar(['Declared axial shell\n(PER_POINT)', '3D witness separation\n(PER_POINT; global upper bound)'], vals, color=['#47879a', '#bd2444'])
ax2.axhline(0.5, color='#182e39', linestyle='--', label='KATANA HT posterior IFU: 0.5 mm')
for (i, v) in enumerate(vals):
    ax2.text(i, v + 0.035, f'{v:.3f} mm', ha='center')
ax2.set_ylabel('Digital geometric distance [mm]')
ax2.set_ylim(0, 1.4)
ax2.legend(loc='upper right', fontsize=9)
ax2.set_title('Axial thickness does not ensure 3D separation')
fig.text(0.51, 0.01, 'One derived crown; no measured preparation, manufactured part or adjustment time. Scanner/pose error UNKNOWN.', ha='center', fontsize=9)
fig.tight_layout(rect=[0, 0.06, 1, 1])
fig.savefig(out / 'designgate.png', dpi=170)
fig.savefig(out / 'designgate.pdf')
plt.close(fig)
