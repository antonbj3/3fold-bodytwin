import thread_budget
import os
os.environ.setdefault('MPLCONFIGDIR', str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'raw/mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import ROOT, read

def main():
    a = read('raw/R1.json')
    d = read('raw/R2.json')
    g = read('raw/R3.json')
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, ax) = plt.subplots(2, 2, figsize=(12, 8))
    ns = np.arange(1, 101)
    u = 1 - 0.05 ** (1 / ns)
    ax[0, 0].plot(ns, 100 * u, color='#215b82', label='0 failures, fixed-n upper 95%')
    ax[0, 0].axhline(10, color='#b54b34', ls='--')
    ax[0, 0].scatter([3, 29, 59], 100 * (1 - 0.05 ** (1 / np.array([3, 29, 59]))), color='#b54b34')
    for n in [3, 29, 59]:
        ax[0, 0].annotate('n=' + str(n), (n, 100 * (1 - 0.05 ** (1 / n))), xytext=(5, 8), textcoords='offset points')
    ax[0, 0].set(xlabel='Specimens / fixed lab population (POPULATION)', ylabel='Failure probability upper bound (%)', title='A  |  Few runouts do not establish low risk', ylim=(0, 75))
    table = a['strength_groups']
    ys = np.arange(len(table))
    lab = ['Unpolished', 'Polished', 'Polished + chamfer', 'Urethane tape', 'Scotch tape']
    ax[0, 1].barh(ys + 0.15, [r['mean_profile']['lower'] for r in table], height=0.27, color='#215b82', label='Mean lower bound')
    ax[0, 1].barh(ys - 0.15, [r['q05_profile']['lower'] for r in table], height=0.27, color='#d08b31', label='5% quantile lower bound')
    ax[0, 1].axvline(300, color='#b54b34', ls='--')
    ax[0, 1].set_yticks(ys)
    ax[0, 1].set_yticklabels(lab)
    ax[0, 1].invert_yaxis()
    ax[0, 1].set(xlabel='MPa / published bar groups (POPULATION)', title='B  |  ISO mean and tail are different observables')
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].text(0.02, -0.24, 'Profile bounds: asymptotic, finite coverage UNKNOWN', transform=ax[0, 1].transAxes, fontsize=8)
    A = d['source_amplitudes']['sinter_size_um']
    x = np.arange(0, 101)
    colors = ['#b54b34', '#215b82', '#528d55']
    for (dur, col) in zip([14, 56, 168], colors):
        ax[1, 0].plot(x, 50 + A * np.minimum(x / dur, 1), color=col, label=str(dur) + '-day imposed ramp')
        row = next((r for r in d['scenarios'] if r['cause'] == 'sinter' and r['ramp_duration_days'] == dur))
        ad = row['median_cause_alarm_day']
        rd = row['first_rejection_day']
        ax[1, 0].scatter([ad], [50 + A * min(ad / dur, 1)], marker='o', color=col)
        ax[1, 0].scatter([rd], [50 + A * min(rd / dur, 1)], marker='x', s=50, color=col)
    ax[1, 0].axhline(120, color='#555', ls='--')
    ax[1, 0].set(xlabel='Day (PHENOMENOLOGICAL)', ylabel='Model signed marginal gap (um)', title='C  |  Weekly alarm o, first rejection x')
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].text(0.02, -0.24, 'Published spatial amplitude; time / gap transfer unmeasured', transform=ax[1, 0].transAxes, fontsize=8)
    z = d['sufficiency_test']
    i = np.arange(2)
    ax[1, 1].bar(i, z['true_gap_um'], color=['#215b82', '#b54b34'], width=0.5)
    ax[1, 1].axhline(120, color='#555', ls='--')
    ax[1, 1].set_xticks(i)
    ax[1, 1].set_xticklabels(['b=0, d=0', 'b=-80, d=80'])
    ax[1, 1].set(ylabel='Actual gap in additive fixture (um)', title='D  |  Identical finished-coupon scan: 0 um', ylim=(0, 165))
    ax[1, 1].text(0.02, 0.94, 'Identity error 0; downstream difference 80 um\nMinimum addition: independent signed reference', transform=ax[1, 1].transAxes, va='top', fontsize=8)
    fig.suptitle('X51: fixed lot evidence + observable drift + explicit physical UNKNOWN', fontsize=14)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.subplots_adjust(hspace=0.48, bottom=0.12)
    fig.savefig(ROOT / 'figures/LAB_ACCEPTANCE_DRIFT.png', dpi=150)
    fig.savefig(ROOT / 'figures/LAB_ACCEPTANCE_DRIFT.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
