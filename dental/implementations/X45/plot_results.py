"""Standalone research figure; use system matplotlib with user packages disabled."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parent
raw = json.loads((ROOT / 'raw/PRIMARY_MEASUREMENTS.json').read_text())
ltd = json.loads((ROOT / 'raw/RESULT_K34_LTD_R1.json').read_text())
isq1 = json.loads((ROOT / 'raw/RESULT_K31_R1.json').read_text())
isq2 = json.loads((ROOT / 'raw/RESULT_K31_R2.json').read_text())
wear1 = json.loads((ROOT / 'raw/RESULT_K34_WEAR_R1.json').read_text())
wear2 = json.loads((ROOT / 'raw/RESULT_K34_WEAR_R2.json').read_text())
wear3 = json.loads((ROOT / 'raw/RESULT_K34_WEAR_R3.json').read_text())
(fig, axs) = plt.subplots(1, 3, figsize=(15, 4.4))
colors = {'OD': '#164c87', 'PISP': '#b85419'}
for (group, ps) in raw['ISQ'].items():
    c = colors[group]
    axs[0].plot([p['day'] for p in ps], [p['mean_ISQ'] for p in ps], '-o', color=c, label=group + ' published mean')
    ys = isq1['validation'][group]
    axs[0].plot([p['day'] for p in ys], [p['predicted'] for p in ys], ':', color=c, label=group + ' R1 decrease')
    ys = isq2['validation'][group]
    axs[0].scatter([p['day'] for p in ys], [p['predicted'] for p in ys], s=65, facecolors='none', edgecolors=c)
axs[0].set(title='K31: observed dip and recovery\nR2 held means:6/6 within1ISQ', xlabel='Healing time (days)', ylabel='Group-mean ISQ [POPULATION]')
axs[0].legend(fontsize=7)
axs[0].grid(alpha=0.2)
grades = [p['grade'] for p in ltd['query_results']]
x = list(range(4))
w = 0.35
axs[1].bar([a - w / 2 for a in x], [p['unaged_Q05_MPa'] for p in ltd['query_results']], width=w, color='#568db3', label='Unaged fit')
axs[1].bar([a + w / 2 for a in x], [p['aged_Q05_MPa'] for p in ltd['query_results']], width=w, color='#cf9443', label='After122C,2bar,8h')
axs[1].set_xticks(x, grades)
axs[1].set(title='K34: exposure changes scale and shape\nSampling confidence is UNKNOWN', ylabel='Fit-derived5% strength quantile (MPa) [POPULATION]')
axs[1].legend(fontsize=8)
axs[1].grid(axis='y', alpha=0.2)
g = 'Admira fusion flow'
ps = raw['WEAR'][g]
axs[2].errorbar([p['cycles'] / 1000 for p in ps], [p['mean_mg'] for p in ps], yerr=[p['SEM_mg'] for p in ps], fmt='ko-', label='Published mean +/-SEM')
for (record, style, label) in [(wear1, 'r:', 'R1 constant rate'), (wear2, 'm--', 'R2 run-in+rate')]:
    pred = record['predictions'][g]
    axs[2].plot([int(k) / 1000 for k in pred], [v['mean_mg'] for v in pred.values()], style, label=label)
for p in wear3['validation'][g]:
    (lo, hi) = p['interval_mg']
    axs[2].plot([p['cycles'] / 1000] * 2, [lo, hi], color='#168460', linewidth=6, alpha=0.35)
axs[2].plot([], [], color='#168460', linewidth=6, alpha=0.35, label='R3 monotone bracket')
axs[2].set(title='K34: later mass loss refutes rate transfer\nIntervals are conditional, not future predictions', xlabel='Lab wear cycles (thousands)', ylabel='Cumulative mass loss (mg) [POPULATION]')
axs[2].legend(fontsize=7)
axs[2].grid(alpha=0.2)
fig.suptitle('X45: measured history to typed research consumers — no patient or clinical-year promotion', fontsize=13)
fig.tight_layout()
out = ROOT / 'figures'
out.mkdir(exist_ok=True)
fig.savefig(out / 'chains.png', dpi=160)
fig.savefig(out / 'chains.svg')
plt.close(fig)
print('Figure: figures/chains.png and chains.svg')
