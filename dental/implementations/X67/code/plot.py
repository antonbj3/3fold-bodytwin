"""Plot-only system Python avoids changing the shared NumPy2 scientific environment."""
import json, pathlib, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]

def load(n):
    return json.loads((ROOT / n).read_text())
r1 = load('attempts/R1_results.json')
r2 = load('attempts/R2_results.json')
r3 = load('attempts/R3_results.json')
r4 = load('attempts/R4_results.json')
(fig, axs) = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
ax = axs[0, 0]
names = ['Legacy', 'eta=1', 'chi=0', 'chi=1', 'k=0.53 W/(m K)', 'h=25 W/(m² K)', 'unbounded fit']
vals = [d['metrics']['within_factor2_count'] for d in r1['diagnostic']]
ax.barh(names[::-1], vals[::-1], color=['#a47732'] + ['#668799'] * 5 + ['#a65252'])
ax.axvline(16 * 0.6, color='#222', ls='--', lw=1)
ax.set_xlim(0, 16)
ax.set_xlabel('conditions within factor 2 (out of 16)')
ax.set_title('A. Native K3 bisection; assumed process inputs', loc='left', fontsize=11)
ax.text(0.03, 0.02, 'Unbounded fit requires eta=1.738; old-table diagnosis only', transform=ax.transAxes, fontsize=8)
ax = axs[0, 1]
folds = r2['folds']
x = np.arange(3)
ax.bar(x - 0.16, [f['baseline']['within_factor2_fraction'] for f in folds], 0.3, label='uncalibrated K3', color='#829ba9')
ax.bar(x + 0.16, [f['bounded']['within_factor2_fraction'] for f in folds], 0.3, label='bounded spatial source', color='#b36759')
ax.axhline(0.6, color='#222', ls='--', lw=1)
ax.set_xticks(x)
ax.set_xticklabels(['Abuhajar\nprimary', 'Pourgiv\nprimary', 'Sindel\nlegacy extraction'])
ax.set_ylim(0, 1.1)
ax.set_ylabel('fraction within factor 2')
ax.legend(fontsize=8)
ax.set_title('B. Whole-study held out; nominal inputs, not physical validation', loc='left', fontsize=10)
ax.text(0.01, 0.94, f"MAE {r2['overall_baseline']['mae_C']:.3f} -> {r2['overall_bounded']['mae_C']:.3f} C", transform=ax.transAxes, fontsize=9)
ax = axs[1, 0]
sc = r4['scenarios']
colors = ['#a65252' if not d['acquisition_width_gate_pass'] else '#31836a' for d in sc]
for (i, (d, c)) in enumerate(zip(sc, colors)):
    ax.plot([d['eta_lower'], d['eta_upper']], [i, i], color=c, lw=5, solid_capstyle='butt')
ax.axvline(0.5, color='#777', ls=':', lw=1)
ax.set_yticks(range(6))
ax.set_yticklabels([f"{d['heater_power_W']} W / {d['scenario'].replace('target_', '').replace('published_spec', 'published ~0.85 C')}" for d in sc], fontsize=8)
ax.invert_yaxis()
ax.set_xlim(-0.02, 1.02)
ax.set_xlabel('effective heat fraction eta, interval')
ax.set_title('C. Certified LP bounds, declared two-volume interpolants', loc='left', fontsize=10)
ax = axs[1, 1]
s = r3['sufficiency']['population_pair']
ax.bar(['[48,48] C', '[37,59] C'], [s['mean_CEM43_A_min'], s['mean_CEM43_B_min']], color=['#668799', '#b36759'])
ax.set_yscale('log')
ax.set_ylabel('mean CEM43 min, constant 60 s')
ax.set_title('D. Identical mean peak = 48 C; dose differs 1024x', loc='left', fontsize=10)
ax.text(0.03, 0.93, 'Derived weighting-law example; no bone injury validation', transform=ax.transAxes, fontsize=8)
fig.suptitle('X67: source magnitude remains unvalidated; measurement requirements are executable', fontsize=13)
for ax in axs.flat:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
(ROOT / 'figures').mkdir(exist_ok=True)
fig.savefig(ROOT / 'figures/drill_heat.png', dpi=160, bbox_inches='tight', pad_inches=0.16)
fig.savefig(ROOT / 'figures/drill_heat.svg', bbox_inches='tight', pad_inches=0.16)
print('Saved PNG/SVG')
