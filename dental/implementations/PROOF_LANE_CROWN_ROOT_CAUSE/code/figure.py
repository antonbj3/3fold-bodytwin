from diagnose import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def lines(ax, tri, M, c, color, label, lw=1, ls='-'):
    mesh = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
    origin = np.r_[c, 0.0] @ M
    segments = trimesh.intersections.mesh_plane(mesh, plane_origin=origin, plane_normal=M[1])
    if not len(segments):
        return
    v = segments @ M.T
    ax.add_collection(LineCollection(v[:, :, [0, 2]], colors=color, linewidths=lw, linestyles=ls, label=label))

def run():
    p = R / 'RESULTS_D.json'
    if not p.exists():
        p = R / 'RESULTS_D_PARTIAL.json'
    new = {r['key']: r for r in read(p)['rows'] if r['status'] == 'GENERATED'}
    a = read(R / 'RESULTS_A.json')['rows']
    r4 = {r['key']: r for r in read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_PREDICTIONS_C.json')['rows'] if r['method'] == 'exact_margin'}
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, axs) = plt.subplots(3, 3, figsize=(13, 12.5), gridspec_kw={'width_ratios': [1, 1, 1.1]})
    for (row, (old, axes)) in enumerate(zip(a, axs)):
        rec = new[old['key']]
        m = load(rec['mesh_path'])
        (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
        base = load(r4[old['key']]['mesh_path'])
        M = np.array(rec['chosen']['world_to_local'])
        c = np.array(rec['chosen']['center_xy'])
        S = m['native_triangles']
        A = m['closure_triangles']
        (qv, qf) = (m['prep_vertices'], m['prep_faces'])
        for ax in axes[:2]:
            lines(ax, S, M, c, '#225ea8', 'Native measured surface', 1.4)
            lines(ax, A, M, c, '#9b9b9b', 'Artificial closure', 1, '--')
            ax.set_aspect('equal')
            sv = S.reshape(-1, 3) @ M.T
            ax.set_xlim(sv[:, 0].min() - 0.5, sv[:, 0].max() + 0.5)
            ax.set_ylim(sv[:, 2].min() - 0.5, sv[:, 2].max() + 0.5)
            ax.set_xlabel('Local transverse coordinate (mm)')
            ax.set_ylabel('Along selected axis (mm)')
        lines(axes[0], base['vertices'][base['faces'][base['roles'] == 1]], M, c, '#d7301f', 'Old one-apex preparation', 1)
        lines(axes[1], qv[qf], M, c, '#16876d', 'New virtual preparation', 1.3)
        lines(axes[1], v[f[roles == 1]], M, c, '#eea126', 'Cavity boundary', 0.8)
        axes[0].set_title(old['family'].capitalize() + ': old preparation')
        axes[1].set_title('New geometric witness')
        ax = axes[2]
        ax.axis('off')
        film = rec['film']
        red = rec['preparation_reduction_mm']
        txt = f"Exact old minimum wall: 0 mm\nOld gap error: {old['gap']['max_error_mm'] * 1000:.1f} µm\n\nNew nominal film:\n{film['min_mm'] * 1000:.3f}–{film['max_mm'] * 1000:.3f} µm\nWall cover bound: {rec['wall']['ray_cover_lower_mm']:.3f} mm\nNative facets kept: 100%\n\nReduction median / maximum:\n{red['median']:.2f} / {red['max']:.2f} mm\nAxis tilt: {rec['chosen']['tilt_degrees']}°\n\nGeometric conjunction: PASS\nClinical preparation: UNKNOWN"
        ax.text(0.02, 0.95, txt, va='top', linespacing=1.25, fontsize=10)
    axs[0, 0].legend(loc='lower left', fontsize=7)
    axs[0, 1].legend(loc='lower left', fontsize=7)
    fig.suptitle('The cone creates the conflict. A feasible digital pair requires a different preparation.\nThree prespecified Teeth3DS teeth; mm and µm are digital geometry, not measured cement.', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(R / 'figures/ROOT_CAUSE.png', dpi=170)
    fig.savefig(R / 'figures/ROOT_CAUSE.svg')
    plt.close(fig)
if __name__ == '__main__':
    run()
