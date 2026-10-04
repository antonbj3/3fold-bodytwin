import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
r6 = json.loads(Path('raw/R6_RESULTS.json').read_text())
world = json.loads(Path('raw/R7_TIMING_WORLDS.json').read_text())['normal_minus_lesion']['groups']
cats = ['normal', 'widening', 'lesion']
labels = ['Normal', 'Vidgning', 'Lesion']
(fig, ax) = plt.subplots(1, 2, figsize=(11, 4.5))
p = np.array([x['sample_risk'] for x in r6['observed_categories']])
ci = np.array([x['IID_simultaneous95_interval'] for x in r6['observed_categories']])
ax[0].errorbar(np.arange(3), p * 100, yerr=np.vstack([p - ci[:, 0], ci[:, 1] - p]) * 100, fmt='o', color='#145c75', capsize=5)
ax[0].set_xticks(range(3), labels)
ax[0].set_ylim(0, 82)
ax[0].set_ylabel("Switching to root treatment (%)")
ax[0].set_title('Publicerade parade grenantal')
for (i, row) in enumerate(r6['observed_categories']):
    ax[0].annotate(f"{row['fallback_RCT']}/{row['n']}", (i, p[i] * 100), xytext=(8, 3), textcoords='offset points')
means = []
for key in ['timing_world_lower', 'timing_world_upper']:
    values = []
    for c in cats:
        times = []
        for g in world.values():
            times.extend((t for (label, t) in zip(g['CBCT_labels'], g[key]) if label == c))
        values.append(np.mean(times))
    means.append(values)
for (vals, label, color) in zip(means, ["World A", "World B"], ['#a34239', '#28678f']):
    ax[1].plot(range(3), vals, 'o-', label=label, color=color)
ax[1].axhline(100, color='gray', ls='--', lw=1, label="100 min , scenario limit")
ax[1].set_xticks(range(3), labels)
ax[1].set_ylabel("Middle-term category, first visit ( min )")
ax[1].set_title("Two compatible, constructed time worlds")
ax[1].legend(fontsize=8)
for a in ax:
    a.spines[['top', 'right']].set_visible(False)
    a.grid(axis='y', alpha=0.2)
fig.suptitle("Diagnosis label → treatment course with intraoperative transition", fontsize=13)
fig.text(0.01, 0.02, "Left: PMID 39325552 ; 1 / 86 CBCT missing. Range: conditional IID , simultaneously 95% .\nRight: no measured patient time; both worlds preserve the same published branch agent/SD and image–branch count.", fontsize=8)
fig.tight_layout(rect=[0, 0.12, 1, 0.92])
fig.savefig('figures/PATHWAY_AND_TIME.png', dpi=170)
fig.savefig('figures/PATHWAY_AND_TIME.pdf')
plt.close(fig)
