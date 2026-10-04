from pathlib import Path
import json
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = Path(__file__).resolve().parent
r = json.loads((P / 'results.json').read_text())
(fig, ax) = plt.subplots(1, 3, figsize=(12, 4.2))
n = r['harvest']['selected_reports']
yes = r['harvest']['reports_with_primary_facets']
no = n - yes
ax[0].bar(['Primary facets', 'Withheld'], [yes, no], color=['#277e8e', '#c65f56'])
ax[0].set_ylim(0, n + 3)
ax[0].set_ylabel('Adjudicated report candidates')
ax[0].set_title('Source verification')
for (i, v) in enumerate([yes, no]):
    ax[0].text(i, v + 0.5, str(v), ha='center')
a = [111.64, 23.96, 77.3, 36.22]
b = [23.96, 111.64, 77.3, 36.22]
ax[1].bar([0, 1], [62.28, 62.28], width=0.6, color='#999')
ax[1].scatter([0, 1], [a[0], b[0]], color='#c65f56', s=65, label='Vestibular mean')
ax[1].axhline(100, color='#777', ls=':', label='Illustrative lab query')
ax[1].set_xticks([0, 1], ['Region map A', 'Permuted B'])
ax[1].set_ylabel('Marginal gap (µm)')
ax[1].set_title('Identical full multiset, different region')
ax[1].legend(fontsize=8, loc='upper right')
ax[1].text(0.03, 0.04, 'Summary identity error: 0\nRegional difference: 87.68 µm', transform=ax[1].transAxes, fontsize=9)
ax[2].bar(['Scoped edges', 'Calibrated edges'], [r['mapping']['unique_existing_edges_with_scoped_facits'], 0], color=['#277e8e', '#999'])
ax[2].set_ylabel('Existing R4 edges')
ax[2].set_title('Reference availability ≠ calibration')
ax[2].text(1, 0.4, '0', ha='center')
fig.suptitle('Dental external reference harvest — source/protocol/region retined', fontsize=13)
fig.tight_layout()
fig.savefig(P / 'FACIT_HARVEST.png', dpi=160)
fig.savefig(P / 'FACIT_HARVEST.pdf', metadata={'CreationDate': datetime.datetime.fromisoformat(json.loads((P / 'PREREG_XFACIT_R1.json').read_text())['frozen_at'])})
plt.close(fig)
