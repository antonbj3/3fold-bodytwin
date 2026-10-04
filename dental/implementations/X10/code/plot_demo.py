from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
r1 = load('RESULTS_R1.json')
r2 = load('RESULTS_R2.json')
r4 = load('RESULTS_R4.json')
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
(fig, ax) = plt.subplots(1, 3, figsize=(13.6, 4.6), constrained_layout=True)
rows = r1['clinical_rows']
x = np.arange(len(rows))
w = 0.35
ax[0].bar(x - w / 2, [r['predicted_group_gap_um'] for r in rows], w, label='40 + global scanner summary', color='#7392af')
ax[0].bar(x + w / 2, [r['value_um'] for r in rows], w, label='Independent published gap', color='#ba6b4f')
ax[0].axhline(120, color='#713b39', ls='--', lw=1)
ax[0].set_xticks(x, ['Omnicam*\nmedian', 'Cara TRIOS*\nmedian', 'i500†\nmean', 'CS3600†\nmean'])
ax[0].set_ylabel('Group marginal-gap summary (µm)')
ax[0].set_title('A  Global-to-gap shortcut fails')
ax[0].legend(fontsize=8, loc='upper left')
ax[0].set_ylim(0, 195)
ax[0].text(0.02, -0.23, '* Cross-year aliases; version match UNKNOWN\n† Local secondary table; primary numeric verification incomplete', transform=ax[0].transAxes, fontsize=7.7)
names = ['Omni', 'Prime', 'Trios 3', 'Trios 4']
commercial = [1, 2, 3, 7]
mat = np.array([[next((z for z in r2['rows'] if z['scanner'] == n and z['scanbody'] == j))['screens'][2]['composed_median_gap_screen_um'] for j in commercial] for n in names])
im = ax[1].imshow(mat, cmap='RdYlGn_r', vmin=70, vmax=310, aspect='auto')
ax[1].set_xticks(range(4), ['SB1', 'SB2', 'SB3', 'SB7'])
ax[1].set_yticks(range(4), names)
for i in range(4):
    for j in range(4):
        ax[1].text(j, i, f'{mat[i, j]:.0f}', ha='center', va='center', color='#172629', fontsize=9)
ax[1].set_title('B  Scanner + commercial scanbody')
fig.colorbar(im, ax=ax[1], label='Composed-median screen (µm)', shrink=0.7)
ax[1].text(0.0, -0.22, 'g0 = 40 µm, lever = 4 mm.\nSeparate medians; actual gap certificate UNKNOWN.', transform=ax[1].transAxes, fontsize=8)
eps = np.linspace(-65, 65, 300)
a = np.sin(np.deg2rad(3))
gap = 40 + np.maximum(0, -(40 + eps) / a)
ax[2].plot(eps, gap, color='#447c77', lw=2, label='Signed axial contact model')
ax[2].axhline(120, color='#713b39', ls='--', lw=1, label='Requested engineering limit')
ax[2].axvline(r4['axial_error_flip_um'], color='#888', ls=':', lw=1)
ax[2].scatter([-57.5, 57.5], [r4['case_bounds']['same_rms_contraction']['gap_upper_um'][0], 40], color=['#ba6b4f', '#447c77'], zorder=4)
ax[2].set_xlim(-65, 65)
ax[2].set_xlabel('Signed axial-wall clearance error (µm)')
ax[2].set_ylabel('Projected marginal gap (µm)')
ax[2].set_title('C  Same error size, opposite seating')
ax[2].legend(fontsize=8, loc='upper right')
ax[2].text(0.0, -0.22, 'Planar walls, 3° half taper, axial insertion.\nModel example; no physical crown-gap validation.', transform=ax[2].transAxes, fontsize=8)
fig.suptitle('Scanner information needed by the fit decision', fontsize=14)
fig.savefig(ROOT / 'FIGURE.png', dpi=190, bbox_inches='tight')
fig.savefig(ROOT / 'FIGURE.pdf', bbox_inches='tight')
plt.close(fig)
print('Figure saved: FIGURE.png and FIGURE.pdf')
