import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from occlusion_operator import H, DATA

def main():
    r1 = json.loads((H / 'round1/results.json').read_text())
    r4 = json.loads((H / 'round4/results.json').read_text())
    flow = json.loads((H / 'round4/force_feasibility_3d.json').read_text())
    models = json.loads((H / 'raw/segmentation_model.json').read_text())
    m2 = json.loads((H / 'raw/segmentation_r2_model.json').read_text())
    (fig, ax) = plt.subplots(2, 3, figsize=(14, 8))
    a = ax.ravel()
    x = np.arange(4)
    for (i, jaw) in enumerate(['upper', 'lower']):
        c = next((v for v in r4['comparisons'] if v['band_mm'] == 0.1 and v['jaw'] == jaw))
        a[i].bar(x - 0.24, c['reference_pp'], 0.24, label='Published independent reference', color='#555555')
        a[i].bar(x, c['predicted_area_pp'], 0.24, label='3D near-surface area', color='#2876a6')
        a[i].bar(x + 0.24, c['predicted_mechanics_pp'], 0.24, label='Conditional support model', color='#de8a42')
        a[i].set_xticks(x)
        a[i].set_xticklabels(['Incisor', 'Canine', 'Premolar', 'Molar'])
        a[i].set_ylabel('Relative share (pp)')
        a[i].set_title(('Upper: Ferrato T-Scan' if i == 0 else 'Lower: Hattori Prescale') + '\nUnmatched cohorts; provisional regions')
        a[i].set_ylim(0, 85)
    a[0].legend(fontsize=7)
    case = 1
    with np.load(DATA / f'{case:03d}_map.npz') as z:
        xy = z['xy']
        proj = z['gap'] <= 0.1
    with np.load(DATA / f'{case:03d}_clearance.npz') as z:
        near = z['distance_mm'] <= 0.1
    only = proj & ~near
    added = near & ~proj
    common = proj & near
    a[2].scatter(xy[common, 0], xy[common, 1], s=3, c='#2876a6', label='Both')
    a[2].scatter(xy[only, 0], xy[only, 1], s=8, c='#de8a42', label='Projection only')
    a[2].scatter(xy[added, 0], xy[added, 1], s=8, c='#7552a2', label='3D proximity only')
    a[2].set_aspect('equal')
    a[2].set_xlabel('RAS x (mm)')
    a[2].set_ylabel('RAS y (mm)')
    a[2].set_title('Case 1 fixed pose, 0.1mm band\nUnsigned sampled clearance')
    a[2].legend(fontsize=7)
    for (rid, (r, offset, color, label)) in enumerate([(r1, -0.07, '#999999', 'Projection'), (r4, 0.07, '#2876a6', '3D proximity')]):
        f = [c for c in r['contrasts'] if 'Open minus Normal' in c['quantity']]
        means = np.array([v['mean'] for v in f])
        intervals = np.array([v['interval'] for v in f])
        a[3].errorbar(means, np.arange(3) + offset, xerr=[means - intervals[:, 0], intervals[:, 1] - means], fmt='o', color=color, label=label)
    a[3].axvline(0, color='black', lw=0.7)
    a[3].set_yticks(range(3))
    a[3].set_yticklabels(['0.05mm', '0.10mm', '0.20mm'])
    a[3].set_xlabel('Open minus normal anterior area (pp)\nDescriptive bootstrap 95%; no force inference')
    a[3].legend(fontsize=8)
    a[3].set_title('Independent bite annotations')
    values1 = [models[j]['validation']['macro_recall'] for j in ['upper', 'lower']]
    values2 = [m2[j]['validation']['macro_recall'] for j in ['upper', 'lower']]
    a[4].bar(np.arange(2) - 0.15, values1, 0.3, label='2D prototype')
    a[4].bar(np.arange(2) + 0.15, values2, 0.3, label='Arch coordinate')
    a[4].axhline(0.85, color='red', ls='--', label='Frozen required 0.85')
    a[4].set_xticks(range(2))
    a[4].set_xticklabels(['Upper', 'Lower'])
    a[4].set_ylim(0.6, 1)
    a[4].set_ylabel('Heldout macro region recall')
    a[4].set_title('Four-region transfer failed\nSeparate heldout arches for the two models')
    a[4].legend(fontsize=7)
    dist = flow['rank_distribution']
    xx = sorted(map(int, dist))
    a[5].bar(xx, [dist[str(k)] for k in xx], color='#2876a6')
    a[5].set_xlabel('Minimum scalar shape channels\nNecessary reaction-flow outer set, provisional regions')
    a[5].set_ylabel('Scan pairs')
    a[5].set_title('Geometry does not identify force distribution')
    fig.suptitle('Bits2Bites: observed closed-bite geometry; force remains unmeasured', fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.text(0.02, 0.005, 'No virtual closure. Source triangles are preserved. Support, sensor transfer and anatomical region identity are not calibrated.', fontsize=9)
    (H / 'figures').mkdir(exist_ok=True)
    fig.savefig(H / 'figures/occlusion_demo.png', dpi=180)
    fig.savefig(H / 'figures/occlusion_demo.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
