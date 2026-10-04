from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = Path(__file__).resolve().parents[1]

def main():
    d = json.load(open(R / 'raw/EXTERNAL_VALIDATION.json'))['groups']['test']['metrics']['upper']
    t = json.load(open(R / 'raw/TRANSPORT_SUMMARY.json'))['aggregate']
    p = json.load(open(R / 'FROZEN_PREDICTIONS.json'))['profiles']
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 150})
    (fig, axs) = plt.subplots(2, 2, figsize=(11, 8), layout='constrained')
    for (ax, (tag, title)) in zip(axs[0], [('R', 'Right force share'), ('posterior', 'Posterior force share')]):
        r = d[tag]
        q = np.array(r['quantile_levels']) * 100
        ax.plot(q, np.array(r['external_quantiles_share']) * 100, 'o-', color='#145b7e', label='Measured T-Scan quantile knots (n=178)')
        ax.plot(q, np.array(r['predicted_quantiles_share']) * 100, 's-', color='#c44b39', label='Frozen area proxy (test n=188)')
        band = np.array(r['bootstrap95_quantiles_share']) * 100
        ax.fill_between(q, band[:, 0], band[:, 1], color='#c44b39', alpha=0.15, label='Proxy patient bootstrap 95%')
        default = 50 if tag == 'R' else 100 * 8 / 14
        ax.axhline(default, color='#666666', ls=':', label='No-geometry uniform14 default')
        ax.set(xlabel='Percentile across people (%)', ylabel='Force-share / proxy (%)', title=title, ylim=(-2, 102))
        ax.grid(alpha=0.15)
        ax.legend(fontsize=7, loc='lower right' if tag == 'R' else 'lower right')
    ax = axs[1, 0]
    x = np.arange(7)
    means = (np.array(p['lower']['area_mean'][:7]) + np.array(p['lower']['area_mean'][7:])) / 2 * 200
    ref = np.array([1.4, 1.2, 2.3, 2.5, 3.6, 15, 24]) * 2
    ax.bar(x - 0.17, means, 0.34, label='Frozen area proxy', color='#c44b39')
    ax.bar(x + 0.17, ref, 0.34, label='Hattori pooled-side digitization', color='#145b7e')
    ax.set(xticks=x, xticklabels=['I1', 'I2', 'C', 'P1', 'P2', 'M1', 'M2'], ylabel='Normalized pooled tooth-type shape (%)', title='Tooth shape differs; sides cluster in 42 people')
    ax.legend(fontsize=7)
    ax.grid(axis='y', alpha=0.15)
    ax = axs[1, 1]
    names = ['No force input', 'One upper region', 'Upper tooth vector', '+ crossing scalar']
    x = np.arange(4)
    for (tag, color, offset) in [('posterior', '#c44b39', -0.17), ('R', '#145b7e', 0.17)]:
        vals = [100, t[tag]['one_upper_region']['mean_width_share'] * 100, t[tag]['full_upper_tooth_marginal']['mean_width_share'] * 100, 0]
        ax.bar(x + offset, vals, 0.34, label=tag, color=color)
    ax.set(xticks=x, xticklabels=names, ylabel='Mean lower-region bound width (pp)', title='Conditional measurement consumer, no physical sensing')
    ax.tick_params(axis='x', labelrotation=18)
    ax.legend(fontsize=7)
    ax.grid(axis='y', alpha=0.15)
    fig.suptitle('X94 — similar central shares conceal wrong distribution tails\nGeometry and observation ports are conditional; cohorts / protocols are unmatched', fontsize=13)
    fig.savefig(R / 'figures/force_validation.png')
    fig.savefig(R / 'figures/force_validation.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
