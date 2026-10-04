"""System python -s gives ABI-compatible NumPy/matplotlib; no package edits."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 130})
targets = json.loads((R / 'raw/R2_PAPER_MOLAR_COMPARISON.json').read_text())
links = json.loads((R / 'raw/R3_TOOTH_LINKS.json').read_text())
rank = json.loads((R / 'raw/R3_RESULTS.json').read_text())
(fig, (a, b)) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.25, 1]}, layout='constrained')
y = np.arange(len(targets))
obs = np.array([r['observed_angular_deg'] for r in targets])
plan = np.array([r['planned_angular_deg'] for r in targets])
gap = np.array([r['angular_tracking_gap_deg'] for r in targets])
a.errorbar(obs, y - 0.11, xerr=4, fmt='o', color='#126e85', label='Observed relative rotation')
a.errorbar(plan, y + 0.11, xerr=4, fmt='s', color='#d38335', label='Planned relative rotation')
a.set_yticks(y, [f"{('A' if r['patient'] == '3485' else 'B')}  tooth {r['fdi']} / incisor {r['anchor_fdi']}" for r in targets])
a.invert_yaxis()
a.set_xlabel('Cumulative angular magnitude at week 8 (deg)')
a.set_title('Published patient scans versus plan')
a.legend(fontsize=8, loc='lower right')
a.grid(axis='x', alpha=0.2)
b.barh(y, gap, color='#704e78')
b.axvline(8, color='#a63e36', linestyle='--', label='Conditional error radius: 8 deg')
b.set_yticks(y, [])
b.invert_yaxis()
b.set_xlabel('SO(3) plan-tracking discrepancy (deg)')
b.set_title('Direction retained in relative pose')
b.legend(fontsize=8, loc='lower right')
b.grid(axis='x', alpha=0.2)
fig.suptitle('PER_TOOTH • reference = moving central incisor • two patient series')
fig.text(0.02, -0.025, 'Bars ±4 deg and dashed 8 deg assume 1 deg error per crown pose; these are sensitivity scenarios, not measured patient confidence intervals.', fontsize=8)
for ext in ['png', 'svg', 'pdf']:
    fig.savefig(R / f'figures/molar_plan_tracking.{ext}', bbox_inches='tight')
plt.close(fig)
typed = json.loads((R / 'raw/R4_TYPE_ENCLOSURES.json').read_text())
typed = [r for r in typed if r['motion'] == 'buccolingual_tip' and r['fdi'] in [17, 27] and (r['week'] == 8) and (r['interval'] == 'cumulative')]
(fig, ax) = plt.subplots(figsize=(9, 5), layout='constrained')
y = np.arange(len(typed))
sign = np.array([1 if r['planned'] >= 0 else -1 for r in typed])
ach = np.array([r['achieved'] for r in typed]) * sign
pl = np.array([abs(r['planned']) for r in typed])
ax.errorbar(ach, y - 0.12, xerr=[r['achieved_radius'] for r in typed], fmt='o', color='#126e85', label='Achieved, signed toward plan')
ax.errorbar(pl, y + 0.12, xerr=[r['planned_radius'] for r in typed], fmt='s', color='#d38335', label='Planned projected crown-axis tip')
ax.set_yticks(y, [f"{('A' if r['patient'] == '3485' else 'B')} tooth {r['fdi']} / incisor {r['anchor_fdi']}" for r in typed])
ax.invert_yaxis()
ax.axvline(0, color='gray', lw=0.8)
ax.set_xlabel('Relative buccolingual crown-coordinate-axis tip at week 8 (deg)')
ax.legend(loc='lower right', fontsize=8)
ax.grid(axis='x', alpha=0.2)
fig.suptitle('Finite signed type enclosure • root torque and absolute tip remain UNKNOWN')
fig.text(0.02, -0.03, 'Error bars follow the exact projection-disk bound assuming 1 deg per crown pose. This input is PHENOMENOLOGICAL; references are moving incisors.', fontsize=8)
for ext in ['png', 'svg']:
    fig.savefig(R / f'figures/typed_tip_tracking.{ext}', bbox_inches='tight')
plt.close(fig)
(fig, axs) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
for (p, ax) in zip(['3485', '6457'], axs):
    rows = [r for r in links if r['patient'] == p]
    x = [r['sum_force_proxy_N'] for r in rows]
    y = [r['lag_signed'] for r in rows]
    ax.scatter(x, y, color='#126e85')
    ax.axhline(0, color='gray', lw=0.7)
    for r in rows:
        if r['fdi'] in [17, 27]:
            ax.annotate(str(r['fdi']), (r['sum_force_proxy_N'], r['lag_signed']), xytext=(3, 3), textcoords='offset points', fontsize=8)
    ax.set_xlabel('Sum of 9 weekly X19 force proxies (N)')
    ax.set_ylabel('1 − signed angular projection fraction')
    rr = rank['by_patient'][p]
    ax.set_title(f"Patient {('A' if p == '3485' else 'B')}: n={rr['teeth']} teeth, rho={rr['rho_force_lag']:.3f}")
    ax.grid(alpha=0.2)
fig.suptitle('X19 geometric proxy fails the frozen descriptive gain gate in both patients')
fig.text(0.02, -0.03, 'PETG stiffness transferred to cellulose-coated Naturaligner; simulated forces are unvalidated. Sum is not force-time exposure. No population inference.', fontsize=8)
for ext in ['png', 'svg']:
    fig.savefig(R / f'figures/x19_force_lag.{ext}', bbox_inches='tight')
plt.close(fig)
print('Figures written')
