from geometry import *

def run():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    fig = plt.figure(figsize=(12, 4))
    norm = Normalize(0, 0.75)
    cmap = plt.get_cmap('viridis')
    inputs = {r['key']: r for r in read(ROOT / 'FROZEN_INPUTS_B.json')['records']}
    for (i, r) in enumerate(read(ROOT / 'FROZEN_EXPORTS.json')['exports']):
        m = npz(r['mesh_path'])
        t = m['vertices'][m['faces'][m['roles'] == 0]]
        target = npz(inputs[r['key']]['private_path'])['target']
        d = distance(target, t.mean(1))
        ax = fig.add_subplot(1, 3, i + 1, projection='3d')
        poly = Poly3DCollection(t, facecolors=cmap(norm(d)), linewidths=0)
        ax.add_collection3d(poly)
        v = t.reshape(-1, 3)
        lo = v.min(0)
        hi = v.max(0)
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(hi - lo)
        ax.view_init(20, 35)
        ax.set_xlabel('x mm')
        ax.set_ylabel('y mm')
        ax.set_zlabel('z mm')
        ax.set_title(r['family'] + '; p95 %.3f mm' % r['predictions']['in_situ']['p95_mm'])
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=fig.axes, shrink=0.55, pad=0.02, label='Facet-centroid distance to original (mm)')
    fig.suptitle('Frozen exterior candidates; inner wall fails — research geometry only')
    fig.savefig(ROOT / 'figures/crown_surfaces.png', dpi=160, bbox_inches='tight')
    plt.close(fig)
if __name__ == '__main__':
    run()
