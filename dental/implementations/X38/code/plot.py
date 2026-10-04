"""Scientific figure. Run with existing system matplotlib (no user NumPy ABI)."""
from pathlib import Path
import json, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = Path(__file__).resolve().parents[1]
r = json.loads((R / 'rounds/R3.json').read_text())
s = json.loads((R / r['statistics_path']).read_text())
r1 = json.loads((R / 'rounds/R1.json').read_text())
(fig, ax) = plt.subplots(2, 2, figsize=(11, 7.5))
ax[0, 0].bar(['Facet-index map', 'Physical-triangle map'], [100 * r1['region_disagreement_fraction'], 100 * r['region_probe']['mismatches'] / r['region_probe']['facets']], color=['#c25a48', '#27847a'])
ax[0, 0].set(ylabel='Mislabeled facets (%)', title='Same crown geometry, reversed facet order\nSynthetic transport labels; no anatomical claim')
ax[0, 0].text(1, 2, '0 mismatches', ha='center')
n = list(range(1, 81))
bounds = [100 * (1 - math.exp(math.log(0.05) / i)) for i in n]
ax[0, 1].plot(n, bounds, color='#266a97')
ax[0, 1].axhline(10, color='#b55d30', ls=':')
ax[0, 1].scatter([3, 29], [100 * s['actual100N']['population_upper_one_sided'], 100 * s['plan29']['upper_if_acceptance_boundary']], color='#b55d30')
ax[0, 1].set(xlabel='Completed specimens, zero failures', ylabel='One-sided95% upper risk (%)', title='Fixed-plan risk certification (POPULATION)')
ax[0, 1].annotate('3:63.16%', (3, 63.16), xytext=(12, 65))
ax[0, 1].annotate('29:9.81%', (29, 9.81), xytext=(35, 22))
for (key, label) in [('binary', 'Failure risk30% vs10% (Fisher)'), ('continuous', 'Mean difference / SD =0.8 (t)')]:
    q = s[key]
    curve = q['search_curve']
    ax[1, 0].plot([v['n_per_group'] for v in curve], [100 * v['power'] for v in curve], label=label)
    ax[1, 0].scatter([q['n_per_group']], [100 * q['achieved_power']])
ax[1, 0].axhline(80, color='#888888', ls=':')
ax[1, 0].legend(fontsize=8)
ax[1, 0].set(xlabel='Specimens per design', ylabel='Power (%)', title='Two-sided alpha0.05; assumed effects (POPULATION)')
cells = s['actual_all_cells']
known = sum((x['n'] - x['unknown_early_censors'] for x in cells))
unknown = sum((x['unknown_early_censors'] for x in cells))
ax[1, 1].bar(['Known horizon endpoint', 'Early censored'], [known, unknown], color=['#27847a', '#d29047'])
ax[1, 1].set(ylabel='Published specimens', title='X25/Zenodo22306621 at5 million cycles\nKeep incomplete outcomes unknown (PER_TOOTH)')
for a in ax.flat:
    a.grid(axis='y', alpha=0.2)
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
fig.tight_layout()
(R / 'figures').mkdir(exist_ok=True)
fig.savefig(R / 'figures/export_lab.png', dpi=160)
fig.savefig(R / 'figures/export_lab.pdf')
plt.close(fig)
