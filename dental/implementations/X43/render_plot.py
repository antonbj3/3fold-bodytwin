"""Optional standard-tool export of the administrative figure; run_all emits SVG without dependencies."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p = Path(__file__).resolve().parent
d = json.loads((p / 'results.json').read_text())
labels = list(d['old_counts'])
x = list(range(len(labels)))
(fig, ax) = plt.subplots(figsize=(9, 4.2))
ax.bar([i - 0.18 for i in x], [d['old_counts'][s] for s in labels], 0.36, label='Frozen K table', color='#9da7af')
ax.bar([i + 0.18 for i in x], [d['proposed_counts'][s] for s in labels], 0.36, label='X43 proposal', color='#2776a8')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('Number of K chains (administrative)')
ax.set_title('52 Northstar chains: reviewed deliveries -> proposed statuses')
ax.set_ylim(0, 38)
ax.legend(frameon=False)
for (i, s) in enumerate(labels):
    ax.text(i - 0.18, d['old_counts'][s] + 0.5, str(d['old_counts'][s]), ha='center', fontsize=9)
    ax.text(i + 0.18, d['proposed_counts'][s] + 0.5, str(d['proposed_counts'][s]), ha='center', fontsize=9)
fig.text(0.12, 0.02, 'No new row becomes PRESENT. No added physical validation. K33 counts as PARTIAL.', fontsize=9)
fig.tight_layout(rect=(0, 0.045, 1, 1))
fig.savefig(p / 'status_counts.png', dpi=160)
fig.savefig(p / 'status_counts.pdf')
