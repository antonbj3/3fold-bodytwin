from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
fig = plt.figure(figsize=(15, 9))
frames = read(PREV / 'raw/CERVICAL_FRAME_DIAGNOSTIC.json')
for (k, rec) in enumerate(inputs()[:3]):
    (v, f, _) = native(rec)
    S = load(rec['private_path'])['target']
    M = np.array(next((r['M'] for r in frames if r['key'] == rec['key'])))
    T = S @ M.T
    for j in range(2):
        ax = fig.add_subplot(2, 3, j * 3 + k + 1, projection='3d')
        a = T if j == 0 else v[f] @ M.T
        ax.add_collection3d(Poly3DCollection(a[::max(1, len(a) // 7000)], facecolors='#b8d7e0', edgecolors='none', alpha=0.9))
        ax.auto_scale_xyz(*a.reshape(-1, 3).T)
        ax.set_box_aspect(np.ptp(a.reshape(-1, 3), axis=0))
        ax.set_title(rec['family'] + (' observed' if j == 0 else ' closed native'))
        ax.set_xlabel('x mm')
        ax.set_ylabel('y mm')
        ax.set_zlabel('z mm')
        ax.view_init(15, 30)
fig.tight_layout()
fig.savefig(R / 'figures/INITIAL_GEOMETRY.png', dpi=150)
