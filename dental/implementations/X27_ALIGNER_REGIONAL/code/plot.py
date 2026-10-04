"""Standalone research figure. Run with PYTHONNOUSERSITE=1 system NumPy/MPL."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})

def main():
    regions = json.loads((R / 'inputs/PARK_REGIONS.json').read_text())
    r1 = json.loads((R / 'round1/results.json').read_text())
    r2 = json.loads((R / 'round2/results.json').read_text())
    r4 = json.loads((R / 'round4/results.json').read_text())
    r5 = json.loads((R / 'round5/results.json').read_text())
    (fig, axes) = plt.subplots(2, 2, figsize=(12.5, 8.5))
    ax = axes[0, 0]
    names = ['Buccogingival', 'Buccal', 'Incisal/occlusal', 'Palatal', 'Palatogingival']
    x = np.arange(5)
    h = [1000 * next((r['median_mm'] for r in regions if r['arm'] == 'TS' and r['region'] == n and (r['quantity'] == 'finished_thickness'))) for n in names]
    g = [1000 * next((r['median_mm'] for r in regions if r['arm'] == 'TS' and r['region'] == n and (r['quantity'] == 'passive_gap'))) for n in names]
    ax.bar(x - 0.18, h, 0.36, label='Formed thickness (measured)', color='#387ca3')
    ax.bar(x + 0.18, g, 0.36, label='Passive gap (measured)', color='#db9750')
    ax.axhline(750, color='gray', ls='--', label='Nominal sheet')
    ax.set_xticks(x)
    ax.set_xticklabels(['Buccal\ngingival', 'Buccal', 'Incisal /\nocclusal', 'Palatal', 'Palatal\ngingival'])
    ax.set_ylabel('Thickness / gap (µm)')
    ax.set_ylim(0, 900)
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title('A  Park TS: surface marginals, no joint tooth map')
    ax = axes[0, 1]
    x = np.arange(3)
    width = 0.19
    nom = [r['practice_N'] for r in r1['summary']]
    reg = [r['regional_N'] for r in r1['summary']]
    calibrated = [1.49] + [r['force_N'] for r in r5['rows']]
    ref = [r['external_N'] for r in r1['summary']]
    sd = [r['external_SD_N'] for r in r1['summary']]
    for (shift, values, label, c) in [(-1.5, nom, 'Nominal ribbon', '#9da9b0'), (-0.5, reg, 'Regional contact, no fit', '#db9750'), (0.5, calibrated, 'Regional + one force calibration', '#387ca3')]:
        ax.bar(x + shift * width, values, width, label=label, color=c)
    ax.errorbar(x + 1.5 * width, ref, yerr=sd, fmt='ko', capsize=3, label='Kaur mean ± SD')
    ax.set_xticks(x)
    ax.set_xticklabels(['Central incisor\ncalibration in R5', 'Canine\nunfit in R5', '2nd premolar\nunfit in R5'])
    ax.set_ylabel('Buccolingual force (N)')
    ax.set_ylim(0, 3.3)
    ax.legend(fontsize=8)
    ax.set_title('B  Type means: unfit errors 16.4% / 8.1%')
    ax = axes[1, 0]
    x = np.arange(6)
    force = [r['force_vector_error'] for r in r2['rows']]
    moment = [r['moment_vector_error'] for r in r2['rows']]
    ax.scatter(x - 0.1, force, label='Force-vector error', color='#387ca3')
    ax.scatter(x + 0.1, moment, label='Moment-vector error', color='#db9750')
    ax.axhline(0.25, color='red', ls='--', label='Frozen 25% gate')
    ax.set_yscale('log')
    ax.set_ylim(0.04, 50)
    ax.set_xticks(x)
    ax.set_xticklabels([r['material'][:3] + '\n#' + str(r['fdi']) for r in r2['rows']])
    ax.set_ylabel('Relative vector error (dimensionless)')
    ax.set_title('C  One 0.3 mm wrench does not transfer to 0.6 mm')
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    for (label, c) in [('R3 selected3', '#387ca3'), ('all16 frozen probes', '#db9750')]:
        points = []
        for error in [0.0, 0.001, 0.0025, 0.005, 0.01, 0.01745, 0.025, 0.0349, 0.05]:
            rs = [r for r in r4['rows'] if r['queries'] == label and r['compression_error_mm'] == error and (r['scaled_wrench_error_N'] == 0.001)]
            vals = [r['max_regional_relative_bound'] for r in rs]
            points.append(max(vals) if all((v is not None for v in vals)) else np.nan)
        ax.plot(np.array([0.0, 0.001, 0.0025, 0.005, 0.01, 0.01745, 0.025, 0.0349, 0.05]) * 1000, points, 'o-', label=label, color=c)
    ax.axhline(0.25, color='red', ls='--', label='±25% stiffness bound')
    ax.axvline(34.9, color='gray', ls=':', label='Park voxel size, not accuracy')
    ax.set_ylim(0, 5)
    ax.set_xlim(0, 51)
    ax.set_ylabel('Worst regional relative stiffness bound')
    ax.set_xlabel('Assumed compression error bound (µm)')
    ax.set_title('D  Metrology requirement, simulator fixture only')
    ax.legend(fontsize=8)
    fig.suptitle('Regional aligner calibration: a limited force-transfer result and an executable measurement design', fontsize=13)
    fig.text(0.02, 0.008, 'B: unmatched arches; known unfit type means, no held study. C: published independent fabrication groups. D: conditional five-patch model; no physical calibration.', fontsize=9)
    fig.tight_layout(rect=[0, 0.025, 1, 0.96])
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(R / f'figures/aligner_regional.{suffix}', dpi=180)
    r6 = json.loads((R / 'round6/results.json').read_text())
    r8 = json.loads((R / 'round8/results.json').read_text())
    r9 = json.loads((R / 'round9/results.json').read_text())
    plt.close(fig)
    (fig, axes) = plt.subplots(1, 2, figsize=(12.5, 5.1))
    x = np.arange(6)
    refs = [r['external_N'] for r in r6['neighbour_rows']]
    sd = [r['external_SD_N'] for r in r6['neighbour_rows']]
    labels = [f"{r['active_fdi']} → {r['response_fdi']}" for r in r6['neighbour_rows']]
    axes[0].errorbar(x, refs, yerr=sd, fmt='ko', capsize=4, label='Kaur Table 4.4 mean ± SD', zorder=5)
    for (delta, result, label, color) in [(-0.15, r6, 'R6 grounded: 1/6 pass', '#db9750'), (0, r8, 'R8 whole arch: 2/6 pass', '#387ca3'), (0.15, r9, 'R9 membrane: 1/6 pass', '#b95c74')]:
        pred = [r['predicted_N'] for r in result['neighbour_rows']]
        axes[0].scatter(x + delta, pred, label=label, color=color, marker='x', s=55, zorder=4)
        axes[1].scatter(x + delta, [100 * r['relative_error'] for r in result['neighbour_rows']], color=color, marker='x', s=55, label=label)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25)
        ax.set_xlabel('Active tooth → responding neighbour (type correspondence)')
    axes[0].set_ylabel('Signed buccolingual force (N)')
    axes[0].set_ylim(-1.48, 0.08)
    axes[0].legend(fontsize=8, loc='lower right')
    axes[0].set_title('A  Unfit neighbouring reactions refute force transfer')
    axes[1].axhline(25, color='red', ls='--', label='Frozen 25% gate')
    axes[1].set_yscale('log')
    axes[1].set_ylim(0.5, 900)
    axes[1].set_ylabel('Absolute relative error (%)')
    axes[1].legend(fontsize=8)
    axes[1].set_title('B  Fitted active-tooth magnitude hides wrong load path')
    fig.suptitle('Whole-arch aligner contact: numerically verified, physical neighbour-force claim fails', fontsize=13)
    fig.text(0.025, 0.025, 'PER_TOOTH / SIMULTANEOUS. One incisor calibration only. Unmatched Taglus sensor / Duran thickness / source arches. R9: 2/9 solver cycles excluded.', fontsize=9)
    fig.tight_layout(rect=[0, 0.07, 1, 0.94])
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(R / f'figures/neighbour_refutation.{suffix}', dpi=180)
    print('Research figures exported PNG/PDF/SVG')
if __name__ == '__main__':
    main()
