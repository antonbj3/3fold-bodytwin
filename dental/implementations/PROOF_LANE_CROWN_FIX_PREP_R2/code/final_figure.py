from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def section(ax, v, f, M, y, color, label, lw=1.1):
    v = v @ M.T
    lines = trimesh.intersections.mesh_plane(trimesh.Trimesh(v, f, process=False), [0, 1, 0], [0, y, 0])
    ax.add_collection(LineCollection(lines[:, :, [0, 2]], colors=color, linewidths=lw, label=label))
    ax.autoscale_view()

def run():
    result = read(R / 'results.json')
    key = '06c45662b50389bd_premolar'
    c = next((x for x in read(R / 'RESULTS_C.json')['rows'] if x['key'] == key and x['material'] == 'KATANA_ML'))
    e = next((x for x in read(R / 'RESULTS_E.json')['rows'] if x['key'] == key and x['material'] == 'KATANA_ML'))
    (cm, em) = (load(c['mesh_path']), load(e['mesh_path']))
    M = cm['M']
    y = c['construction']['c'][1]
    (tv, tf, _) = native(next((x for x in inputs() if x['key'] == key)))
    (fig, axes) = plt.subplots(1, 3, figsize=(15, 5), gridspec_kw={'width_ratios': [1.2, 1, 1]})
    ax = axes[0]
    section(ax, tv, tf, M, y, '#69798a', 'Closed native T', 1.4)
    section(ax, em['prep_vertices'], em['prep_faces'], M, y, '#118a7e', 'Adaptive E core Q', 2)
    section(ax, cm['prep_vertices'], cm['prep_faces'], M, y, '#d56c26', 'Flat C core Q', 1.6)
    section(ax, em['vertices'], em['faces'], M, y, '#7567b2', 'E crown', 0.7)
    ax.axhline(cm['base'], color='#4c5963', ls=':', lw=1)
    ax.set_aspect('equal')
    ax.set_xlabel('Local x (mm)')
    ax.set_ylabel('Local z (mm)')
    ax.set_title('One premolar, same finish plane\nPER_POINT section')
    ax.legend(fontsize=8, loc='lower left')
    r = result['retained_substance']
    ax = axes[1]
    values = [r['flat_C_full_preparation_mm3'], r['adaptive_E_full_preparation_mm3']]
    bars = ax.bar(['Flat C', 'Adaptive E'], values, color=['#d56c26', '#118a7e'], width=0.55)
    ax.bar_label(bars, fmt='%.2f')
    ax.set_ylim(0, 220)
    ax.set_ylabel('Retained full preparation (mm³)')
    ax.set_title('Preserve more of this tooth\nPER_TOOTH digital volume')
    ax.text(0.5, 0.38, f"+{r['additional_retained_mm3']:.2f} mm³ retained\n{r['removed_volume_reduction_percent']:.1f}% less removed", transform=ax.transAxes, ha='center')
    ax.text(0.5, 0.02, 'Welded die approximation; biological location unknown', transform=ax.transAxes, ha='center', fontsize=7)
    ax = axes[2]
    x = np.arange(3)
    ax.bar(x - 0.18, [2, 2, 0], width=0.36, label='Assembly geometry', color='#118a7e')
    ax.bar(x + 0.18, [0, 0, 0], width=0.36, label='Complete crown chain', color='#7567b2')
    ax.set_xticks(x, ['Anterior', 'Premolar', 'Molar'])
    ax.set_ylim(0, 6.5)
    ax.set_yticks(range(7))
    ax.set_ylabel('Teeth out of six per type')
    ax.set_title('All 18 original teeth retained\nPER_TOOTH outcomes')
    ax.legend(fontsize=8)
    ax.text(1, 0.2, 'Complete chain: 0/18', ha='center', fontsize=9)
    fig.suptitle('Native-contained preparation: 4 digital assembly candidates, no qualified full crowns', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    (R / 'figures').mkdir(exist_ok=True)
    fig.savefig(R / 'figures/DEMO.png', dpi=180)
    fig.savefig(R / 'figures/DEMO.pdf')
    plt.close(fig)
if __name__ == '__main__':
    run()
