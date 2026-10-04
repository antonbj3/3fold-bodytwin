"""Standalone research figure, hypothetical fixtures only."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_r1 import ROOT
r5 = json.loads((ROOT / 'rounds/R5/results.json').read_text())
r6 = json.loads((ROOT / 'rounds/R6/results.json').read_text())
r4 = json.loads((ROOT / 'rounds/R4/results.json').read_text())
low = next((r for r in r5['records'] if r['world'] == 0 and r['support'] == 0 and (r['reaction_error_N'] == 0.005)))
query = next((r for r in r6['records'] if r['world'] == 0 and r['support'] == 0 and (r['law'] == 'nonlinear') and (r['reaction_error_N'] == 0.005)))
(fig, axes) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
ax = axes[0]
ax.axvline(query['facit_N'], color='#222222', ls='--', label='Nonlinear native reference')
for (y, label, lo, hi, c) in [(2, 'Elastic work', 205, 205, '#b23a48'), (1, 'Low-rate linear fit', low['lower_N'], low['upper_N'], '#cb8a28'), (0, 'Query-rate observation', query['lower_N'], query['upper_N'], '#197d68')]:
    ax.plot([lo, hi], [y, y], lw=5, color=c)
    ax.scatter([(lo + hi) / 2], [y], color=c, s=35)
ax.set_yticks([0, 1, 2], ['Query-rate observation', 'Low-rate linear fit', 'Elastic work'])
ax.set_xlabel('Local support force [N]')
ax.set_xlim(195, 300)
ax.set_ylim(-0.5, 2.5)
ax.legend(loc='upper left', fontsize=8)
ax.set_title('Same fixture: different observations')
ax.text(0.03, 0.08, 'Native velocity: 195 µm/s\nCalibration: 1–3 µm/s\nReaction error: hypothetical ±5 mN', transform=ax.transAxes, fontsize=8)
ax = axes[1]
ss = [s for s in r4['summaries'] if s['residual_error_mm'] == 1e-06]
ax.plot([s['intervals'] for s in ss], [s['fraction_top_two_le5pp'] for s in ss], 'o-', color='#436a9f')
ax.axhline(0.8, color='#b23a48', ls='--', label='Frozen 80% gate')
ax.set_xscale('log', base=2)
ax.set_xticks([512, 8192, 32768], ['512', '8192', '32768'])
ax.set_ylim(0, 1)
ax.set_xlabel('Force intervals per ramp')
ax.set_ylabel('Loaded queries with width ≤5 pp')
ax.set_title('Paired work still misses its gate')
ax.annotate('14/18 = 77.8%', (32768, 14 / 18), xytext=(-80, -35), textcoords='offset points', fontsize=9)
ax.text(0.03, 0.54, 'Hypothetical residual: 1 nm\nShared drift bound: ±1 µm', transform=ax.transAxes, fontsize=8)
ax.legend(loc='upper left', fontsize=8)
fig.suptitle('Conditional simulations; physical native force remains UNKNOWN', fontsize=11)
(ROOT / 'figures').mkdir(exist_ok=True)
fig.savefig(ROOT / 'figures/acquisition_results.png', dpi=180)
fig.savefig(ROOT / 'figures/acquisition_results.pdf')
