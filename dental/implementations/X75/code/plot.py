"""Use system NumPy/Matplotlib: python3 -s code/plot.py."""
import json, os, pathlib
P = pathlib.Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR'] = str(P / 'raw/mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': 'white'})
    r1 = {(r['case'], r['site']): r for r in map(json.loads, (P / 'raw/PER_SITE_R1.jsonl').read_text().splitlines())}
    r2 = [json.loads(l) for l in (P / 'raw/PER_SITE_R2.jsonl').open()]
    s = json.loads((P / 'SUMMARY_SUFFICIENCY.json').read_text())
    h = s['heat']
    (fig, axes) = plt.subplots(2, 2, figsize=(11.3, 8.4))
    (a, b, c, d) = axes.flat
    x = [r1[r['case'], r['site']]['rays'][1]['envelope_chord_mm'] for r in r2]
    y = [r['rays'][1]['envelope_chord_mm'] for r in r2]
    a.scatter(x, y, s=18, c=['#4477aa' if r['kind'].startswith('X8') else '#cc6677' for r in r2], alpha=0.75)
    a.plot([0, 40], [0, 40], color='.7', lw=1)
    a.set(xlabel='Fixed image plane: envelope chord (mm)', ylabel='Local arch plane proxy: chord (mm)', title='A  Site frame changes the measured geometry')
    a.text(0.04, 0.96, '95 source sites; 94 local planes\nTrue anatomical plane: UNKNOWN', va='top', transform=a.transAxes, fontsize=9)
    retained = next(((r, q, side) for r in r2 for q in r['rays'] for (side, v) in q['profiles'].items() if v['rejection_reason'] is None))
    (r, q, side) = retained
    key = 'd' + str(int(q['depth_mm'])) + '_' + side
    with np.load(r['raw_file']['path']) as raw:
        xx = raw[key + '_depth_mm']
        gg = raw[key + '_gray_interp']
    b.plot(xx, gg, color='#4477aa')
    b.set(xlabel='Distance inward from annotated exterior (mm)', ylabel='Observed gray value (NOT HU)', title='B  A retained contrast profile')
    b.set_ylim(0, float(max(gg)) * 1.35)
    b.text(0.04, 0.96, r['case'] + ' / ' + r['site'] + '\n24 / 570 retained; cortex accuracy UNKNOWN', va='top', transform=b.transAxes, fontsize=9)
    with np.load(P / 'raw/EXACT_SUMMARY_WITNESS.npz') as raw:
        xx = (raw['bin_edges'][:-1] + raw['bin_edges'][1:]) / 2
        c.step(xx, raw['initial_A'], where='mid', label='A: ascending', color='#4477aa')
        c.step(xx, raw['initial_B'], where='mid', label='B: reversed', color='#cc6677')
    c.set(xlabel='Position (dimensionless)', ylabel='Constructed unit contrast', title='C  Exact histogram twins retain different positions', ylim=(0, 1.2))
    c.legend(frameon=False, fontsize=9)
    c.text(0.04, 0.96, 'Histogram and mean identity error = 0\nConstructed states; no measured temperature', va='top', transform=c.transAxes, fontsize=9)
    vals = [h['functional_A_unit_contrast'], h['functional_B_unit_contrast']]
    bounds = h['rigorous_model_intervals_A_B']
    d.bar([0, 1], vals, color=['#4477aa', '#cc6677'], width=0.55)
    d.errorbar([0, 1], vals, yerr=[[v - z[0] for (v, z) in zip(vals, bounds)], [z[1] - v for (v, z) in zip(vals, bounds)]], fmt='none', color='black', capsize=5)
    d.set(xticks=[0, 1], xticklabels=['A', 'B'], ylabel='Heat-kernel functional (unit contrast)', title='D  Spatial functional separates the same summary')
    d.text(0.04, 0.96, 'External closed-form heat kernel\nDifference lower bound = %.6f\nPhysical gray-to-heat map: UNKNOWN' % h['difference_lower_bound'], va='top', transform=d.transAxes, fontsize=9)
    d.set_ylim(0, 0.5)
    fig.suptitle('X75: observed site profiles; physical bone quality still uncalibrated', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(P / 'FIGURE_X75.png', dpi=170)
    fig.savefig(P / 'FIGURE_X75.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
