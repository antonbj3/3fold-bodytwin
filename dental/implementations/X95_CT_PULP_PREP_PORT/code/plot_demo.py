"""Standalone scientific plot; system interpreter avoids NumPy/Matplotlib ABI mismatch."""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]

def main():
    rows = json.loads((ROOT / 'raw/R1_ROWS.json').read_text())
    w = json.loads((ROOT / 'raw/R1_WITNESSES.json').read_text())
    r3 = json.loads((ROOT / 'raw/R3_ROWS.json').read_text())
    fdi = next((x['fdi'] for x in w if x['downstream_difference_mm'] > 0))
    (fig, axes) = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for (ax, name) in zip(axes[0], ['REGIONAL_A_LOW_X', 'REGIONAL_B_HIGH_X']):
        row = next((x for x in rows if x['fdi'] == fdi and x['operation'] == name))
        a = np.load(row['spatial_export'])
        depth = a['removed'].sum(axis=0) * 0.3
        depth = np.where(a['tooth'].any(axis=0), depth, np.nan)
        im = ax.imshow(depth, origin='lower', vmin=0, vmax=2.1, cmap='viridis', extent=[0, depth.shape[1] * 0.3, 0, depth.shape[0] * 0.3])
        projection = a['pulp'].any(axis=0)
        if projection.any():
            ax.contour((np.arange(depth.shape[1]) + 0.5) * 0.3, (np.arange(depth.shape[0]) + 0.5) * 0.3, projection, levels=[0.5], colors=['#e84545'], linewidths=1)
        ax.set_title(f"FDI {fdi}: {('A' if name.endswith('LOW_X') else 'B')} | V={row['removed_volume_mm3']:.3f} mm³\ncut gap={row['cut_exact_pulp_cube_min_mm']:.3f} mm")
        ax.set_xlabel('Local geometric X [mm]')
        ax.set_ylabel('Local geometric Y [mm]')
    fig.colorbar(im, ax=list(axes[0]), label='Removed axial column height [mm]; red=full pulp projection')
    ax = axes[1, 0]
    x = np.arange(len(w))
    width = 0.38
    ax.bar(x - width / 2, [v['cut_min_A_mm'] for v in w], width, label='Regional A')
    ax.bar(x + width / 2, [v['cut_min_B_mm'] for v in w], width, label='Regional B')
    ax.set_xticks(x)
    ax.set_xticklabels([str(v['fdi']) for v in w])
    ax.set_xlabel('Canonical FDI in one CT case')
    ax.set_ylabel('Exact voxel pulp-cube to cut gap [mm]')
    ax.legend()
    ax.set_title('Identical total volume, different local clearance')
    ax = axes[1, 1]
    x = np.arange(len(r3))
    ax.bar(x, [v['uniform_cap_deficit_voxels'] for v in r3], color='#c66a2d')
    ax.set_xticks(x)
    ax.set_xticklabels([str(v['fdi']) for v in r3])
    ax.set_xlabel('Canonical FDI in one CT case')
    ax.set_ylabel('Baseline cap cells that must remain [count]')
    ax.set_title('Exact minimum local depth-rule departure\n1mm research reserve, same total removal budget')
    fig.suptitle('Published P48 annotations → complete virtual preparation solids\nDigital voxel geometry; physical dentin, crown performance and clinical suitability UNKNOWN', fontsize=13)
    (ROOT / 'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT / 'figures/pulp_prep_port.png', dpi=180)
    fig.savefig(ROOT / 'figures/pulp_prep_port.pdf')
    plt.close(fig)
    (ROOT / 'raw/PLOT_ENVIRONMENT.json').write_text(json.dumps({'python': sys.version, 'numpy': np.__version__, 'matplotlib': matplotlib.__version__}, indent=2) + '\n')
if __name__ == '__main__':
    main()
