import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
d = json.loads((ROOT / 'raw/X18.json').read_text())
fit = json.loads((ROOT / 'raw/X18_V2.json').read_text())
normalchanges = [dict(fdi=x['fdi'], planfit_outer_peak_MPa=y['coulomb_peak_interval_MPa'][1], facet_outer_peak_MPa=x['coulomb_peak_interval_MPa'][1]) for (x, y) in zip(d['rows'], fit['rows'])]
rankfails = sum((r['plane_fit_rank_failures'] for r in d['rows']))
patches = sum((r['patches'] for r in d['rows']))
(fig, ax) = plt.subplots(1, 2, figsize=(11, 4))
xx = np.arange(len(d['rows']))
before = [r['before_peak_MPa'] for r in d['rows']]
lo = [r['coulomb_peak_interval_MPa'][0] for r in d['rows']]
hi = [r['coulomb_peak_interval_MPa'][1] for r in d['rows']]
ax[0].bar(xx - 0.18, before, 0.36, label='100 N, vertikalt')
ax[0].bar(xx + 0.18, lo, 0.36, label='20 N local + gradient/friction')
ax[0].errorbar(xx + 0.18, lo, yerr=[np.zeros(len(xx)), np.array(hi) - lo], fmt='none', color='black', capsize=3)
ax[0].set_xticks(xx, [r['fdi'] for r in d['rows']])
ax[0].set_ylabel('Conditional peak traction voltage [MPa ]')
ax[0].set_xlabel('Predicted tooth ; separate local models')
ax[0].legend(fontsize=8)
ax[0].set_title('Scenario, no measured patient voltage')
xf = np.arange(len(normalchanges))
ax[1].bar(xf - 0.18, [x['planfit_outer_peak_MPa'] for x in normalchanges], 0.36, label='Planpassning i patch')
ax[1].bar(xf + 0.18, [x['facet_outer_peak_MPa'] for x in normalchanges], 0.36, label='Incidenta triangelfacetter')
ax[1].set_xticks(xf, [x['fdi'] for x in normalchanges])
ax[1].set_ylabel('Conditional outer boundary [ MPa ]')
ax[1].legend(fontsize=8)
ax[1].set_title(f'{rankfails}/{patches} patch plan under-defined')
for axis in ax:
    axis.grid(axis='y', alpha=0.2)
fig.tight_layout()
fp = ROOT / 'FIGURE.png'
fig.savefig(fp, dpi=170)
plt.close(fig)
