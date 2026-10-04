import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R = Path(__file__).resolve().parents[1]

def run():
    r1 = json.loads((R / 'R1_RESULTS.json').read_text())
    r2 = json.loads((R / 'R2_RESULTS.json').read_text())
    (fig, axes) = plt.subplots(1, 3, figsize=(13, 4.2), layout='constrained')
    keys = ['marginal', 'intaglio', 'occlusal', 'axial']
    x = np.arange(4)
    for mode in ['balanced', 'marginal', 'intaglio']:
        axes[0].plot(x, [r1['regions']['D1'][mode][k]['signed_bias_um'] for k in keys], marker='o', label=mode)
    axes[0].plot(x, [r2['datum_regions']['D1'][k]['signed_bias_um'] for k in keys], marker='o', label='external datum (simulated)')
    axes[0].set_xticks(x, keys, rotation=15)
    axes[0].set_ylabel('Signed normal bias [µm]')
    axes[0].set_title('Same part, different registration')
    axes[0].legend(fontsize=8)
    cells = json.loads((R / 'raw/LITERATURE_CELLS.json').read_text())['cells']
    crown = [c for c in cells if c['object'] in ['crown', 'crown_coping']]
    axes[1].scatter(np.arange(len(crown)), [c['mean_um'] for c in crown], c=['#f08b32' if c['process'] == 'printed' else '#2878b5' for c in crown])
    axes[1].errorbar(np.arange(len(crown)), [c['mean_um'] for c in crown], yerr=[c['sd_um'] for c in crown], fmt='none', color='gray', alpha=0.7)
    axes[1].set_title('Independent published crown RMS')
    axes[1].set_ylabel('Manufacturing RMS mean ± between-crown SD [µm]')
    axes[1].set_xlabel('2 CoCr groups + 9 zirconia groups; two studies')
    axes[1].text(0.03, 0.04, 'Aggregates do not validate raw scan fields', transform=axes[1].transAxes, fontsize=8)
    sf = r2['sufficiency']
    axes[2].bar([0, 1], [sf['minimum_a_um'], sf['minimum_b_um']], color=['#b13c3c', '#2878b5'])
    axes[2].axhline(0, color='black', lw=1)
    axes[2].set_xticks([0, 1], ['State A', 'State B'])
    axes[2].set_ylabel('Minimum local marginal clearance [µm]')
    axes[2].set_title('Identical regional bias and RMS')
    axes[2].text(0.03, 0.96, f"Identity error = {sf['identity_error_um']:g} µm\nConsumer difference = {sf['difference_um']:g} µm", transform=axes[2].transAxes, va='top', fontsize=9)
    fig.suptitle('K48: simulated measurement pipeline + independent literature context; physical validation UNKNOWN', fontsize=11)
    fig.savefig(R / 'figures/metrology.png', dpi=180)
    fig.savefig(R / 'figures/metrology.pdf')
    plt.close(fig)
if __name__ == '__main__':
    run()
