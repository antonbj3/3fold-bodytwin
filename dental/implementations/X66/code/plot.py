import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import *
r = read('R1_RESULTS.json')
data = read('raw/observations.json')
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 150})
(fig, axs) = plt.subplots(2, 2, figsize=(12, 8.5))
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.09, top=0.83, wspace=0.39, hspace=0.52)
ax = axs[0, 0]
x = np.arange(2)
meas = [a['measured_mean_ratio'] for a in r['L02']['paired_ratios']]
pred = [a['current_thickness_rule_ratio'] for a in r['L02']['paired_ratios']]
ax.bar(x - 0.18, meas, 0.36, label='Published surviving-mean ratio', color='#236b8e')
ax.bar(x + 0.18, pred, 0.36, label='Existing t² operation', color='#cc8355')
ax.set_xticks(x, ['3Y/5Y multilayer', '4Y/5Y multilayer'])
ax.set_ylim(0, 2.6)
ax.set_ylabel('1.0 mm / 0.7 mm fracture-load ratio')
ax.set_title('Cantilever FDP: thickness proxy fails transfer')
ax.legend(frameon=False, fontsize=8)
ax = axs[0, 1]
for (group, color) in zip('ABC', ['#236b8e', '#699d63', '#cc8355']):
    rows = [x for x in data['L03'] if x['group'] == group]
    ax.plot(range(4), [x['mean_um'] for x in rows], 'o-', label=group, color=color)
ax.axhline(120, color='#8a3042', ls='--', label='Descriptive screen: 120 µm')
ax.set_xticks(range(4), ['Distal', 'Mesial', 'Buccal', 'Palatal'])
ax.set_ylabel('Published post-cementation mean gap (µm)')
ax.set_title('Four repeated sites, 8 copings per group')
ax.legend(frameon=False, fontsize=8)
ax = axs[1, 0]
(lo, hi) = r['L04']['coefficient_rounding_interval']
ax.plot([100 * lo, 100 * hi], [1, 1], lw=9, color='#236b8e', solid_capstyle='butt')
ax.plot(r['L04']['reported_fraction'] * 100, 1, 'D', ms=9, color='#8a3042', label='Reported endpoint')
ax.plot(r['L04']['equation_value_from_printed_coefficients'] * 100, 1, 'o', color='white', mec='black', label='Printed Eq. 1 result')
ax.set_yticks([])
ax.set_xlim(9.5, 17.5)
ax.set_ylim(0.5, 1.5)
ax.set_xlabel('Nominal cyclic load / monotonic ultimate load (%)')
ax.set_title('Coefficient rounding box excludes 12.86%')
ax.legend(frameon=False, fontsize=8, loc='upper left')
ax.text(0.03, 0.04, 'Observed levels: 10% all 4 runout; 15% ≥1 failure/design\nLoads rounded to 10 N; intermediate/absolute force UNKNOWN', transform=ax.transAxes, fontsize=8)
ax = axs[1, 1]
upper = r['L05']['native_slab_eta0to1_deltaT_C_range'][1]
old = r['L05']['incident_dose_without_new_information_upper_C']
ax.barh([0, 1], [old, upper], color=['#cc8355', '#236b8e'], height=0.5)
ax.set_yticks([0, 1], ['Incident dose', 'Measured transmission'])
ax.set_xlabel('Radiation-only fixed-slab ΔT envelope width (°C)')
ax.set_title('Conditional range width falls by 75.7%')
ax.text(0.03, 0.98, 'Surface heat deposition + fixed thermal properties\nFull luting/real tooth temperature remains UNKNOWN', transform=ax.transAxes, va='top', fontsize=8)
ax.plot([4.1], [1], 'D', color='#8a3042', label='Published no-luting 4.1°C')
ax.legend(frameon=False, fontsize=8, loc='lower right')
ax.set_ylim(-0.6, 2.0)
fig.suptitle('X66 source-conditioned lab observations — PENDING INDEPENDENT REVIEW\nNo whole-chain or clinical predictive validation', fontsize=12)
for ext in ['png', 'svg', 'pdf']:
    fig.savefig(ROOT / ('figure.' + ext), bbox_inches='tight')
print('Figure written')
