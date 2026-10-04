from common_r3 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0, str(V6 / 'code'))
from contact import polygons
from memory_height import height

def run():
    r = next((r for r in read(ROOT / 'LAB_EXPORTS.json')['selected'] if r['family'] == 'premolar_crown'))
    key = r['key']
    p = npz(DATA / 'inputs' / (key + '.npz'))
    xy = p['contact_xy']
    f = p['contact_faces']
    native = p['contact_reference_gap_unclipped']
    ceiling = p['contact_ceiling_mm']
    (fig, axes) = plt.subplots(1, 4, figsize=(13, 4), sharex=True, sharey=True, constrained_layout=True)
    methods = [('R2', OLD_DATA / 'R2/outer_first_joint' / key / 'mesh.npz'), ('A', DATA / 'A/spatial_sheet' / key / 'mesh.npz'), ('B', DATA / 'B/explicit_sheet' / key / 'mesh.npz'), ('C', DATA / 'C/explicit_sheet' / key / 'mesh.npz')]
    (native_poly, _, _, _) = polygons(xy, f, native)
    for (ax, (tag, path)) in zip(axes, methods):
        source_shapes = [bary @ xy[f[idx]] for (idx, bary) in native_poly.items()]
        ax.add_collection(PolyCollection(source_shapes, facecolors='#359bce', alpha=0.55, edgecolors='none', label='native band'))
        if path.exists():
            m = npz(path)
            tri = m['vertices'][m['faces'][m['face_roles'] == 0]]
            gap = ceiling - height(tri, xy)
            (pp, _, _, _) = polygons(xy, f, gap)
            shapes = [bary @ xy[f[idx]] for (idx, bary) in pp.items()]
            ax.add_collection(PolyCollection(shapes, facecolors='#e07b27', alpha=0.5, edgecolors='none', label='constructed band'))
        else:
            ax.text(0.5, 0.5, 'construction rejected', transform=ax.transAxes, ha='center')
        ax.set_title(tag)
        ax.set_aspect('equal')
        ax.autoscale_view()
        ax.set_xlabel('site x (mm)')
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('site y (mm)')
    axes[0].legend(loc='upper left', fontsize=8)
    fig.suptitle('Selected premolar: blue = native, orange = constructed 0–0.1 mm band\nPER_SURFACE_REGION; supplied pose; source-derived contact prescription')
    fig.savefig(ROOT / 'figures/premolar_contact.png', dpi=170)
    plt.close(fig)
if __name__ == '__main__':
    run()
