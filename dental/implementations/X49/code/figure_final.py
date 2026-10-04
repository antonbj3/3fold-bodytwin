"""Standalone research figure using the compatible installed system matplotlib."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parents[1]

def main():
    r = json.loads((HERE / 'results.json').read_text())
    run = Path(r['run_directory'])
    (fig, ax) = plt.subplots(1, 3, figsize=(13, 5.2), gridspec_kw={'width_ratios': [1, 1.2, 1]})
    s = r['sufficiency']
    x = np.array([0.0, 1.0])
    ax[0].bar(x - 0.16, [s['summary_a_mm'], s['summary_b_mm']], 0.3, label='Mean', color='#688aab')
    ax[0].bar(x + 0.16, [s['a_min_mm'], s['b_min_mm']], 0.3, label='Regional minimum', color='#e79c53')
    ax[0].axhline(0.5, color='#b33', ls='--', label='Named IFU criterion')
    ax[0].set(xticks=x, xticklabels=['State A', 'State B'], ylabel='Thickness [mm]', title='Same mean, different decision')
    ax[0].legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.14))
    ax[0].text(0.5, 0.96, 'Identity error = 0 mm', transform=ax[0].transAxes, ha='center', va='top', fontsize=8)
    names = ['axial_valid', 'axial_penetration', 'axial_near']
    colors = {'PASS': '#27845c', 'FAIL': '#b33', 'UNKNOWN': '#b08017'}
    for (i, name) in enumerate(names):
        q = json.loads((run / (name + '.json')).read_text())['rules']['occlusal_contact']
        v = q['nominal_minimum_gap_mm']
        (lo, hi) = q['robust_interval_mm']
        ax[1].errorbar(v, i, xerr=[[v - lo], [hi - v]], fmt='o', color=colors[q['status']], capsize=4)
        ax[1].text(hi + 0.006, i, q['status'], va='center', fontsize=8, color=colors[q['status']])
    ax[1].axvspan(-0.03, 0.1, color='#27845c', alpha=0.1, label='Frozen lab band')
    ax[1].set(yticks=range(3), yticklabels=['Valid', 'Penetration', 'Larger uncertainty'], xlabel='Projected gap [mm]', title='Fixed-XY axial error envelope', xlim=(-0.25, 0.19))
    ax[1].legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.14))
    cases = [x['case'] for x in r['x1b']]
    health = [x['rules']['mesh_health'] for x in r['x1b']]
    for (i, (case, status)) in enumerate(zip(cases, health)):
        ax[2].text(0, i, case, va='center')
        ax[2].text(0.3, i, status, color=colors[status], va='center')
        ax[2].text(0.78, i, 'PASS', color=colors['PASS'], va='center')
    ax[2].set(xlim=(-0.05, 1.07), ylim=(-0.7, 4.8), xticks=[], yticks=[], title='Five actual X1b designs')
    ax[2].text(0.3, -0.45, 'Triple health', ha='left', fontsize=8)
    ax[2].text(0.78, -0.45, '3MF', ha='left', fontsize=8)
    ax[2].text(0, 4.6, '2 prep-mesh failures;\nphysical checks UNKNOWN', fontsize=8)
    ax[2].invert_yaxis()
    fig.suptitle('X49 — digital geometry capability; no fabricated parts or clinical validation', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(HERE / 'figures/design_gate.png', dpi=180)
    fig.savefig(HERE / 'figures/design_gate.pdf')
if __name__ == '__main__':
    main()
