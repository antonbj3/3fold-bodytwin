"""Run with system Python -s to use its matching NumPy/Matplotlib ABI."""
from dental_release.paths import expand as _release_expand
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1])
    a = parser.parse_args()
    rows = list(csv.DictReader((a.input / 'raw/R8_teeth.csv').open()))
    rows = [r for r in rows if r['domain_truncated'] == 'False']
    order = ['central_incisor', 'lateral_incisor', 'canine', 'premolar1', 'premolar2', 'molar1', 'molar2', 'molar3']
    labels = ['I1', 'I2', 'C', 'P1', 'P2', 'M1', 'M2', 'M3']
    (fig, axes) = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle('Pulpy3D + ToothFairy2: calibrated annotated sample\nTotal hard tissue; physical boundary accuracy and clinical applicability unknown', fontsize=14)
    vals = [[float(r['horn_boundary_min_mm']) for r in rows if r['tooth_type'] == g] for g in order]
    axes[0, 0].boxplot(vals, labels=labels, showfliers=False)
    axes[0, 0].set_ylabel('Horn-region nearest whole-tooth boundary (mm)')
    axes[0, 0].set_title('New anatomical information, grouped by tooth type')
    axes[0, 0].text(0.02, 0.97, 'Digital two-boundary envelope: ±0.520 mm\nTop 2 mm of pulp; crown side is an area proxy', va='top', transform=axes[0, 0].transAxes, fontsize=9)
    for (i, g) in enumerate(order):
        axes[0, 0].text(i + 1, 0.04, str(len(vals[i])), ha='center', transform=axes[0, 0].get_xaxis_transform(), fontsize=8)
    cal = json.loads((a.input / 'raw/R8_calibration.json').read_text())
    r6 = json.loads((a.input / 'raw/R6_pairs.json').read_text())
    n = sum((r.get('paired', False) for r in r6))
    sw = sum((r.get('aligned', False) for r in r6))
    accepted = sum((r.get('accepted', False) for r in cal))
    axes[0, 1].bar(['Exact CT pairs', 'Fixed side swap', '2-anchor + held-out'], [n, sw, accepted], color=['#457b9d', '#e9c46a', '#2a9d8f'])
    axes[0, 1].set_ylabel('Cases')
    axes[0, 1].set_title('Independent annotation referent can reject pairing')
    for (i, v) in enumerate([n, sw, accepted]):
        axes[0, 1].text(i, v + 3, str(v), ha='center')
    axes[0, 1].set_ylim(0, max(n, sw, accepted) * 1.15)
    mats = list(csv.DictReader((a.input / 'material_summary.csv').open()))
    products = ['3Y_HT', '5Y_UTML', 'LS2_adhesive', 'LS2_conventional_lower_bound']
    x = np.arange(4)
    for (i, delta) in enumerate([0, 0.15, 0.3]):
        y = [100 * float(next((r['excluded_uniform_fraction'] for r in mats if r['product'] == p and float(r['buffer_mm']) == 0.5 and (float(r['annotation_delta_each_mm']) == delta)))) for p in products]
        axes[1, 0].bar(x + (i - 1) * 0.25, y, width=0.25, label=f'annotation sensitivity ±{delta:g} mm / boundary')
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(['HT/HTML', 'UTML', 'LS2 adhesive', 'LS2 conventional'], rotation=12)
    axes[1, 0].set_ylim(0, 100)
    axes[1, 0].set_ylabel('Excluded uniform-horn offset requests (%)')
    axes[1, 0].set_title('Illustrative buffer x=0.5 mm; product-specific wall cards')
    axes[1, 0].legend(fontsize=7, loc='upper left')
    axes[1, 0].text(0.02, 0.02, 'Necessary local proxy constraint; no full-crown or clinical veto', transform=axes[1, 0].transAxes, fontsize=8)
    q = np.load(_release_expand('@DENTAL_WORK_ROOT@/X12/R8_P1_36_paired.npz'))
    t = q['tooth']
    p = q['pulp']
    sp = q['spacing']
    points = np.argwhere(p)
    y = int(np.median(points[:, 1]))
    (z, xv) = t[:, y, :].shape
    axes[1, 1].imshow(t[:, y, :], origin='lower', cmap='Greys', extent=[0, xv * sp[2], 0, z * sp[0]])
    axes[1, 1].contour((np.arange(xv) + 0.5) * sp[2], (np.arange(z) + 0.5) * sp[0], p[:, y, :], levels=[0.5], colors=['#d62828'], linewidths=1)
    axes[1, 1].set_xlabel('Local x (mm)')
    axes[1, 1].set_ylabel('Local crown-oriented z (mm)')
    axes[1, 1].set_title('One paired lower molar: tooth exterior and pulp')
    axes[1, 1].text(0.02, 0.02, 'Voxel STL + SDF + 0.5 mm uniform-offset export\nRed = externally supplied pulp annotation', transform=axes[1, 1].transAxes, fontsize=8)
    for ax in axes.flat:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    (a.output / 'figures').mkdir(exist_ok=True, parents=True)
    fig.savefig(a.output / 'figures/atlas_demo.png', dpi=180)
    fig.savefig(a.output / 'figures/atlas_demo.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
