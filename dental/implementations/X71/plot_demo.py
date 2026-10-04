import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
here = Path(__file__).resolve().parent
r = json.loads((here / 'raw/RESULTS_R2.json').read_text())
items = [p for p in r['individual_packages'] if p['measurement'] != 'M05']
labels = {'M01': 'scan/reference', 'M02': 'dry + cement film', 'M03': 'cement dynamics', 'M04': 'fracture pilot12', 'M05': 'full fracture72', 'M06': 'retention/aging24', 'M07': 'patch force + strain', 'M08': 'four-support force', 'M09': 'local implant slip', 'M10': 'HA + bone mechanics', 'M11': 'anatomy + pose', 'M12': 'orthodontic forces', 'M13': 'optics + resin', 'M14': 'spacer/batch40', 'M15': 'paired state/operator'}
(fig, (ax, bx)) = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={'width_ratios': [2, 1]})
for (i, p) in enumerate(items):
    (lo, hi) = p['value_per_operator_hour_interval']
    ax.plot([lo, hi], [i, i], color='#246a87', lw=3)
    ax.scatter([lo, hi], [i, i], s=18, color='#246a87')
ax.set_yticks(range(len(items)), [p['measurement'] + ' ' + labels[p['measurement']] for p in items])
ax.invert_yaxis()
ax.set_xlabel('Potential evaluable queries / operator hour\nAssumed boxes; full72 cost omitted:60 metrology records unpriced')
ax.set_title('Acquisition packages, including prerequisites')
ax.grid(axis='x', alpha=0.25)
bx.bar(['Original 72\nfracture plan', 'Retention 24\nalternative'], [72, 48], color=['#246a87', '#d07d37'])
bx.set_ylim(0, 80)
bx.set_ylabel('Available fracture endpoint slots')
bx.set_title('Identical manufactured summary\n72 specimens; 12 per cell')
bx.text(0.5, 75, 'Identity error 0; difference 24', ha='center', fontsize=10)
fig.suptitle('X71: a conditional first-week plan; no new physical measurements', fontsize=14)
fig.tight_layout()
fig.savefig(here / 'MEASUREMENT_PRIORITY.png', dpi=180)
fig.savefig(here / 'MEASUREMENT_PRIORITY.svg')
print('MEASUREMENT_PRIORITY.png')
