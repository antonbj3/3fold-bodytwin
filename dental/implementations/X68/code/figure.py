"""Publication-style standalone plot. Run with system Matplotlib: python3 -s."""
import json, pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]
read = lambda n: json.loads((ROOT / n).read_text())

def main():
    rows = {r['id']: r for r in read('raw/arms.json')}
    r1 = read('results_R1.json')
    r2 = read('results_R2.json')
    r3 = read('results_R3.json')
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
    (fig, axes) = plt.subplots(2, 2, figsize=(12, 8.6), layout='constrained')
    ax = axes[0, 0]
    kinds = ['practice', 'cortex', 'drilling', 'both']
    x = np.arange(4)
    for (shift, split, color, label) in [(-0.18, 'study', '#2563a5', 'Leave one study out'), (0.18, 'family', '#ea8c23', 'Leave one lab/implant family out')]:
        ax.bar(x + shift, [r1['scores'][split]['candidate_ridge'][k] for k in kinds], width=0.34, label=label, color=color)
    ax.axhline(np.log(1.5), color='gray', ls=':', label='Frozen absolute limit ln(1.5)')
    ax.set_xticks(x, ['Area + density', '+ cortex', '+ drilling', '+ both'])
    ax.set_ylim(0, 0.61)
    ax.set_ylabel('Mean of study medians |ln(prediction / mean torque)|')
    ax.set_title('A  Published PU means: added scalars do not transfer')
    ax.legend(fontsize=8, loc='upper left')
    for (ax, cal, p3, title) in zip([axes[0, 1], axes[1, 0]], r2['calibrations'], r3['predictions'], ['B  Same-regime calibration: no cortex', 'C  Same-regime calibration: 1.5 mm cortex']):
        rs = [rows[i] for i in cal['group']]
        d = np.array([r['Df_mm'] for r in rs])
        t = np.array([r['torque_Ncm'] for r in rs])
        sd = np.array([r['sd_Ncm'] for r in rs])
        D = rs[0]['D_mm']
        grid = np.linspace(d.min(), d.max(), 201)
        ax.axhspan(20, 30, color='#dbefe3', alpha=0.5, label='Research target 25 ±20% Ncm')
        ax.plot(grid, cal['a_Ncm'] + cal['b_Ncm_mm2'] * (D * D - grid * grid), color='#2563a5', label='2-point area proxy (closure)')
        ax.plot(grid, np.polyval(p3['poly_c_b_a'], D * D - grid * grid), ls='--', color='#b73b64', label='3-point quadratic (closure)')
        ax.errorbar(d, t, yerr=sd, fmt='o', color='black', capsize=3, label='Published means ± SD')
        ax.scatter(d[[0, -1]], t[[0, -1]], s=115, facecolors='none', edgecolors='#2563a5', zorder=4, label='2-point calibration endpoints')
        interval = p3['inverse']['Df_outer_interval_mm']
        ax.hlines(8, *interval, lw=5, color='#6b5194', label='25 Ncm root outer set*')
        ax.set_xlabel('Final bore diameter Df (mm)')
        ax.set_ylabel('Peak insertion torque (Ncm)')
        ax.set_title(title)
        ax.set_ylim(3, 44)
        ax.legend(fontsize=7.8, loc='upper right')
        ax.text(0.03, 0.2, '*Conditional exact mean + monotonicity\nSD is not a hard physical bound', transform=ax.transAxes, fontsize=8)
    ax = axes[1, 1]
    ax.bar([0, 1], [52.5, 56.8], yerr=[0.74, 1.33], capsize=4, color=['#2563a5', '#ea8c23'], width=0.6)
    ax.set_xticks([0, 1], ['No cortex\nDf = 3.6 mm', '1.5 mm cortex\nDf = 3.8 mm'])
    ax.set_ylabel('Published mean ISQ ± SD (dimensionless)')
    ax.set_ylim(0, 70)
    ax.set_title('D  Identical torque summary, different measured stability')
    ax.text(0.04, 0.88, '16.2 Ncm in both arms · bitwise identity error = 0\nISQ difference = 4.3 · radial pressure unidentified', transform=ax.transAxes, fontsize=9)
    fig.suptitle('X68: the missing information is a local response and regional contact state', fontsize=14)
    fig.savefig(ROOT / 'FIGURE_X68.png', dpi=180)
    fig.savefig(ROOT / 'FIGURE_X68.svg')
    print('FIGURE_X68.png and FIGURE_X68.svg written; all observations are POPULATION means.')
if __name__ == '__main__':
    main()
