import os
os.environ['MPLCONFIGDIR'] = '/tmp/proof_lane_combine_implant_matplotlib'
import pathlib, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = pathlib.Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
(fig, ax) = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
fig.suptitle('One CBCT site record: digital geometry is answerable; physical risk is not', fontsize=15)
pop = json.loads((P / 'raw/POPULATION_R1.json').read_text())
guides = ['freehand', 'pilot_guided', 'fully_guided', 'dynamic_navigation']
names = ['Freehand', 'Pilot', 'Fully guided', 'Dynamic']
for (off, seg, color, name) in [(-0.12, 'mental_endpoint_5mm_proxy', '#c95c38', 'Mental proxy'), (0.12, 'molar_nearest_FDI', '#28768f', 'Molar')]:
    rows = [next((r for r in pop if r['guide'] == g and r['segment'] == seg and (r['target'] == 0.95)))['fixed_accept_not_certified'] for g in guides]
    v = np.array([r['fraction'] for r in rows]) * 100
    ci = np.array([r['ci95'] for r in rows]) * 100
    ax[0, 0].errorbar(np.arange(4) + off, v, yerr=np.maximum(0, np.stack([v - ci[:, 0], ci[:, 1] - v])), fmt='o', color=color, capsize=4, label=name)
ax[0, 0].set(xticks=np.arange(4), xticklabels=names, ylim=(0, 105), ylabel='Sites (%), image-cluster bootstrap 95% CI', title='Fixed 2 mm accepts; conditional guide bound cannot certify')
ax[0, 0].legend(loc='center left')
ax[0, 0].text(0.02, 0.08, '95% model target; not an injury rate\nGuide moments treated as exact population moments', transform=ax[0, 0].transAxes, fontsize=8)
ch = json.loads((P / 'raw/RESEARCH_CHOICES_R2.json').read_text())
ch = [r for r in ch if r['guide'] == 'fully_guided' and r['target'] == 0.95]
old = [r['fixed_2mm_research_choice']['L_mm'] for r in ch]
new = [r['finite_research_choice']['L_mm'] if r['finite_research_choice'] else np.nan for r in ch]
x = np.arange(len(ch))
ax[0, 1].plot(x, old, 'o-', color='#777', label='TF2 fixed 2 mm')
ax[0, 1].plot(x, new, 's-', color='#28768f', label='Fully guided 95% scenario')
ix = np.where(np.isnan(new))[0]
ax[0, 1].scatter(ix, [5] * len(ix), marker='x', s=50, color='#c95c38', label='No candidate in menu')
ax[0, 1].set(xlabel='22 exact-pose sites in four CBCT cases', ylabel='Maximum menu length (mm)', yticks=[6, 8, 10, 12, 14], title='Dimensions are recomputed from the canal voxels')
ax[0, 1].legend(fontsize=8)
rr = list(csv.DictReader((P / 'raw/DRILL_ADVANCE_R3.csv').open()))
groups = {}
for r in rr:
    groups.setdefault((r['case'], r['fdi']), []).append(r)
for vals in groups.values():
    xx = [float(r['advance_mm']) for r in vals]
    yy = [float(r['gap_lower_mm']) for r in vals]
    ax[1, 0].plot(xx, yy, alpha=0.65, lw=1.4)
ax[1, 0].axhline(2, color='#333', ls='--')
ax[1, 0].axvline(1.5, color='#333', ls=':', alpha=0.4)
ax[1, 0].set(xlabel='Extra cylindrical drill advance (mm), prescribed scenario', ylabel='Distance to digital release union (mm)', title='15 of 21 initially passing sites lose the 2 mm margin')
ax[1, 0].text(0.02, 0.92, 'No real drilling or anatomical safety measurement', transform=ax[1, 0].transAxes, fontsize=8, bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'))
counts = [2581, 748, 80, 22, 0]
yy = np.arange(5)
ax[1, 1].barh(yy, [np.log10(c + 1) for c in counts], color=['#9aa3ac', '#637d8d', '#28768f', '#17646e', '#c95c38'])
ax[1, 1].set(yticks=yy, yticklabels=['Inherited virtual sites', 'Dense canal pairs', 'Same-pose bone profiles', 'Both inputs', 'Calibrated physical plans'], xticks=[], xlabel='Bar length: log10(count + 1)', title='Missing calibration is retained in the plan')
ax[1, 1].invert_yaxis()
for (y, n) in enumerate(counts):
    ax[1, 1].text(np.log10(n + 1) + 0.04, y, str(n), va='center')
ax[1, 1].set_xlim(0, 3.8)
fig.savefig(P / 'figures/site_plan.png', dpi=160)
fig.savefig(P / 'figures/site_plan.svg')
fig.savefig(P / 'figures/site_plan.pdf')
