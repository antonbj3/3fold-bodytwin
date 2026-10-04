"""Plot summary with system NumPy/matplotlib; no data/model mutation."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = Path(__file__).resolve().parents[1]
j = json.load(open(P / 'raw/PLOT_DATA.json'))
fields = list(j['r1']['fields'])
names = ["Overbite", "Overjet", 'Korsbett', "M1 right", "M1 left"]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
(fig, ax) = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
a = np.load(j['case_map_path'])
g = a['gap']
xy = a['xy']
sel = g <= 0.1
p = ax[0, 0].scatter(xy[sel, 0], xy[sel, 1], c=np.clip(g[sel], -0.3, 0.1), s=5, cmap='coolwarm_r', vmin=-0.3, vmax=0.1, rasterized=True)
ax[0, 0].set(aspect='equal', xlabel="Assumed right axle ( mm )", ylabel='Antagen anterioraxel (mm)', title=f"{j['case']}: projected proximity at ≤0,10 mm")
fig.colorbar(p, ax=ax[0, 0], label='Projicerat signerat gap (mm)')
ax[0, 0].text(0.01, 0.98, "Negative gap = projected overlap\nPhysical contact and force UNKNOWN", va='top', transform=ax[0, 0].transAxes, fontsize=9)
up = [10 * q + t for q in [1, 2] for t in range(1, 9)]
lo = [10 * q + t for q in [4, 3] for t in range(1, 9)]
row = json.load(open(j['case_row_path']))
matrix = np.zeros((16, 16))
for r in row['pairs']:
    matrix[up.index(r['upper_fdi']), lo.index(r['lower_fdi'])] = r['projected_proximity_area_mm2']
im = ax[0, 1].imshow(matrix, aspect='auto', cmap='viridis')
ax[0, 1].set_xticks(range(16), lo, rotation=90)
ax[0, 1].set_yticks(range(16), up)
ax[0, 1].set(xlabel="Calculated FDI, lower jaw", ylabel="Calculated FDI , upper jaw", title="Dentised projected proximity area")
fig.colorbar(im, ax=ax[0, 1], label="mm² (including overlap)")
x = np.arange(5)
w = 0.25
for (offset, rd, label, color) in [(-w, j['r1'], 'Kontakt/gap R1', '#4985b4'), (0, j['r2'], 'Opposition R2', '#df8934')]:
    ax[1, 0].bar(x + offset, [rd['fields'][f]['first_report']['balanced_accuracy'] for f in fields], w, label=label, color=color)
ax[1, 0].bar(x + w, [j['r2']['fields'][f]['x7_first_report']['balanced_accuracy'] for f in fields], w, label='X7 geometri', color='#718270')
ax[1, 0].axhline(0.65, color='#b62929', ls='--', lw=1, label='Fryst 0,65-grind')
ax[1, 0].set_xticks(x, names)
ax[1, 0].set(ylim=(0, 1), ylabel="Balanced category of conformity", title="R1 / R2 : original labels before negative correction")
ax[1, 0].legend(fontsize=8, loc='upper left')
orders = {'negative': 0, 'reduced': 1, 'normal': 2, 'increased': 3}
points = []
for r in j['reportrows']:
    if r['field'] != 'overjet':
        continue
    measure = j['opposition'][r['case_id']]['incisor_measures_candidate_mm'].get('overjet_mm')
    if measure is None:
        continue
    vals = [orders[v] for v in r['reference_all'] if v in orders]
    if not vals:
        continue
    jitter = (int(r['case_id'][1:]) % 13 - 6) * 0.008
    ax[1, 1].plot([measure] * len(vals), np.array(vals) + jitter, 'o-', ms=3, lw=0.8, color='#b53333' if r['report_discordant'] else '#758a9a', alpha=0.7)
    points.append(measure)
ax[1, 1].set_yticks(range(4), ['Negativ', 'Reducerad', 'Normal', "Increased"])
ax[1, 1].axvline(0, ls=':', color='gray')
ax[1, 1].set(xlabel="Overjet of Unvalidated Incisive Candidates ( mm )", ylabel='Klinikerrapportens kategori', title="Red lines: different reports for the same case")
fig.suptitle(f"X21 · {j['geometry']['cases']} tandvisa kartor · {j['geometry']['cases_with_projected_penetration_beyond_0_03mm']} poser med gap <−0,03 mm\nProjected mesh geometry is conditional; clinical and physical certification is missing", fontsize=14)
(P / 'figures').mkdir(exist_ok=True)
fig.savefig(P / 'figures/contact_report_demo.png', dpi=160)
fig.savefig(P / 'figures/contact_report_demo.pdf')
plt.close(fig)
