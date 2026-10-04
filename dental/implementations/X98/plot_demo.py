"""Standalone publication/export figure; data input is the computed CSV."""
import csv
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
output = Path(sys.argv[1])
with (output / 'TABLE.csv').open() as f:
    plot = [r for r in csv.DictReader(f) if r['guide'] == 'fully_guided' and r['system'] == 'ProofLane_EV_Guided_2017']

def values(key):
    return [float(r[key]) for r in plot]
(fig, axs) = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
x = list(range(len(plot)))
axs[0].plot(x, values('nominal_gap_mm'), 'o', label='TF2 planned cylinder')
axs[0].plot(x, values('swept_gap_mm'), 'x', label='Revision union + max1mm full-radius tool envelope')
axs[0].axhline(2, color='black', ls='--', label='2mm reference')
axs[0].set_ylabel('Digital gap / mm')
axs[0].legend(fontsize=8)
axs[0].set_xticks(x)
axs[0].set_xticklabels([r['site_id'].replace('ToothFairy2P_', 'P').replace('/FDI', '/') for r in plot], rotation=70, fontsize=7)
r = plot[0]
axs[1].bar(['Entry bound', 'Apex bound', 'Rotation term', 'Combined guide'], [float(r[k]) for k in ['guide_entry_mm', 'guide_apex_mm', 'guide_rotation_mm', 'guide_combined_mm']], color=['#4c78a8', '#4c78a8', '#f58518', '#e45756'])
axs[1].set_ylabel('Conditional95% budget / mm')
axs[1].set_title('Guide = max(entry, apex) + rotation; physical joint coverage UNKNOWN', fontsize=10)
fig.suptitle('X98 digital research scenario: source-resolved clearance; physical nerve safety UNKNOWN', fontsize=11)
fig.savefig(output / 'FIGURE.png', dpi=160)
fig.savefig(output / 'FIGURE.pdf')
plt.close(fig)
