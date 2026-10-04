"""Use /usr/bin/python3 -s to avoid the host NumPy2/system-matplotlib ABI clash."""
from dental_release.paths import expand as _release_expand
import json, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X30'))

def main():
    result = json.loads((ROOT / 'results.json').read_text())
    comp = result['external_population_comparison']
    (fig, axs) = plt.subplots(2, 2, figsize=(11.4, 8.6))
    ax = axs[0, 0]
    x = np.arange(len(comp))
    width = 0.35
    ax.bar(x - width / 2, [100 * r['annotation_second_path_fraction'] for r in comp], width, label='Pulpy annotation surrogate', color='#d4513c')
    ax.bar(x + width / 2, [100 * r['external_second_canal_fraction'] for r in comp], width, label='Independent CBCT Table 2', color='#245a86')
    ax.set_xticks(x, [str(r['fdi']) for r in comp])
    ax.set_ylabel('Second canal / section multiplicity (%)')
    ax.set_xlabel('FDI (mandibular incisors)')
    ax.legend(fontsize=8)
    ax.set_title('Population calibration fails the frozen 15-pp gate', fontsize=10)
    topo = list(csv.DictReader((ROOT / 'raw/R2_topology.csv').open()))
    ax = axs[0, 1]
    types = [('Incisors', [31, 32, 41, 42]), ('Canines', [33, 43]), ('Premolars', [34, 35, 44, 45]), ('Molars', [36, 37, 38, 46, 47, 48])]
    stable = []
    vanished = []
    for (_, fs) in types:
        rs = [r for r in topo if int(r['fdi']) in fs]
        stable.append(100 * np.mean([r['scenario_stable'] == 'True' for r in rs]))
        vanished.append(100 * np.mean([r['scenario_canal_vanished'] == 'True' for r in rs]))
    x = np.arange(4)
    ax.bar(x - 0.175, stable, 0.35, label='Count unchanged in all 3 scenarios', color='#407861')
    ax.bar(x + 0.175, vanished, 0.35, label='Canal disappears after erosion', color='#cb8a36')
    ax.set_xticks(x, [t[0] for t in types])
    ax.set_ylabel('Teeth (%)')
    ax.set_ylim(0, 100)
    ax.set_title('One-voxel edits expose topology uncertainty', fontsize=10)
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    s = result['R3']['summary']
    v = np.array([r['required_point_precision_P5_P50_P95_mm'] for r in s]) * 1000
    ax.errorbar([r['window_axial_mm'] for r in s], v[:, 1], yerr=[v[:, 1] - v[:, 0], v[:, 2] - v[:, 1]], fmt='o-', color='#245a86', capsize=4, label='Sufficient point budget: P5 / median / P95')
    ax.axhline(np.sqrt(3) * 0.3 / 2 * 1000, color='#d4513c', ls='--', label='0.3-mm voxel half diagonal')
    ax.set_xlabel('Axial window support (mm)')
    ax.set_ylabel('Sufficient point error for radius ±20% (µm)')
    ax.set_title('Conservative metrology budget, not a voxel-size limit', fontsize=10)
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    paths = [r for r in map(json.loads, (DATA / 'R1_paths.jsonl').open()) if r['case'] == 'P1' and r['fdi'] == 36]
    for (i, r) in enumerate(paths):
        p = np.array(r['smooth_points_mm'])
        p -= p[0]
        ax.plot(p[:, 2], p[:, 0], label='Saved branch ' + str(i + 1))
    ax.set_aspect('equal')
    ax.set_xlabel('X relative to lower graph endpoint (mm)')
    ax.set_ylabel('Z relative to lower graph endpoint (mm)')
    ax.set_title('P1 / FDI36: partial annotation branch paths', fontsize=10)
    ax.legend(fontsize=8)
    fig.suptitle('Root-canal annotation atlas: measurable geometry and unresolved anatomy', fontsize=14)
    fig.text(0.04, 0.018, 'External anatomy: Verma et al., DOI 10.12659/PJR.901840, Tables 1–4. Cohorts/operators differ.\nLengths are partial paths; radii are window circumcircles. Clinical working length, apex and fracture probability remain UNKNOWN.', fontsize=9)
    fig.tight_layout(rect=[0, 0.07, 1, 0.95])
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/atlas.png', dpi=160)
    fig.savefig(ROOT / 'figures/atlas.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
