from fc_common import *
import matplotlib.pyplot as plt
from skimage.measure import marching_cubes
F = functional()
key = '8017101de9a27cf6_molar_crown'
p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
m = npz(DATA / 'diagnostic_oracles' / key / 'mesh.npz')
tri = a['source_triangles']
tri = tri[tri.mean(1)[:, 2] >= p['margin_z']]
ref = F.sample(tri, 4000)
ext = F.sample(m['vertices'][m['faces'][m['face_roles'] == 0]], 4000)
(v, f, _, _) = marching_cubes(p['cavity'], 0, spacing=(0.25,) * 3)
v += p['origin']
(fig, axs) = plt.subplots(1, 2, figsize=(12, 5))
for (i, ax) in enumerate(axs):
    ax.scatter(ref[:, i], ref[:, 2], s=2, label='native source')
    ax.scatter(ext[:, i], ext[:, 2], s=2, alpha=0.5, label='oracle height envelope')
    ax.scatter(v[:, i], v[:, 2], s=1, alpha=0.2, label='public cavity')
    ax.axhline(float(p['margin_z']), c='k')
    ax.set_aspect('equal')
    ax.legend()
    ax.set_xlabel(('x' if i == 0 else 'y') + ' mm')
    ax.set_ylabel('z mm')
fig.tight_layout()
fig.savefig(ROOT / 'figures/diagnostic_section.png', dpi=140)
print('source bounds', ref.min(0), ref.max(0), 'margin', p['margin_z'])
