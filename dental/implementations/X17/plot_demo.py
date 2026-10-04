from dental_release.paths import expand as _release_expand
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parent
r = [json.loads((ROOT / f'raw/round_R{i}.json').read_text()) for i in range(1, 5)]
cases = json.loads((ROOT / 'raw/cases_R1.json').read_text())
profiles = json.loads((ROOT / 'raw/profiles.json').read_text())
witness = json.loads((ROOT / 'raw/witnesses_R3.json').read_text())
case = cases[0]
name = case['case']
w = witness[name + '_nonresponders']
(fig, ax) = plt.subplots(2, 3, figsize=(15, 8), layout='constrained')
a0 = ax[0, 0]
for c in cases:
    pr = profiles[c['case']]
    e = np.array(pr['eligible'])
    z = np.array(pr['z_mm'])
    area = np.array(pr['area_mm2'])
    a0.plot(z[e] - z[e][0], area[e], alpha=0.6, lw=1, label=c['case'].split('_')[-1])
a0.set(xlabel='Distance from trimmed start (mm)', ylabel='Annotated section area (mm²)', title='12 measured ToothFairy2 pharynx profiles')
a1 = ax[0, 1]
rs = case['scenarios']
mm = sorted(set((s['advancement_mm'] for s in rs)))
lo = [min((s['minimum_mm2'] for s in rs if s['advancement_mm'] == m)) for m in mm]
hi = [max((s['minimum_mm2'] for s in rs if s['advancement_mm'] == m)) for m in mm]
a1.fill_between(mm, lo, hi, alpha=0.25, label='Unknown transfer + 0–6° rotation')
a1.plot(mm, [next((s['minimum_mm2'] for s in rs if s['advancement_mm'] == m and s['rotation_deg'] == 0 and (s['transfer'] == 0.5))) for m in mm], label='Illustration: transfer 0.5, rotation 0°')
a1.set(xlabel='Commanded incisor advancement (mm)', ylabel='Interior minimum area (mm²)', title=_release_expand('@DENTAL_CASE_ID@: conditional scenarios, not physical bounds'))
a1.legend(fontsize=8)
a2 = ax[0, 2]
external = r[1]['external_comparisons']
xx = np.arange(2)
b = 0.22
a2.bar(xx - b, [x['observed_area1_mm2'] for x in external], b, label='Published post area')
a2.bar(xx, [x['area1_predicted_mm2'] for x in external], b, label='Measured AP + lateral reconstruction')
a2.bar(xx + b, [x['proportional_volume_predicted_area1_mm2'] for x in r[0]['external_comparisons']], b, label='Measured volume proportional rule')
a2.set(xticks=xx, xticklabels=['Responders\nn=15', 'Nonresponders\nn=16'], ylabel='Minimum area (mm²)', title='External reference : Shi 2023 , Table 2 / 4')
a2.legend(fontsize=8)
a3 = ax[1, 0]
z = np.array(w['z_mm'])
a3.plot(z - z[0], w['area0_mm2'], label='Measured baseline')
a3.plot(z - z[0], w['area_low_mm2'], label='Conditional low-min witness')
a3.plot(z - z[0], w['area_high_mm2'], label='Conditional high-min witness')
a3.set(xlabel='Distance from trimmed start (mm)', ylabel='Section area (mm²)', title='Same exact profile volume, different minimum')
a3.legend(fontsize=8)
a4 = ax[1, 1]
sc = r[2]['scenarios']
a4.scatter([s['required_slice_fraction'] * 100 for s in sc], [s['width_mm2'] for s in sc])
a4.axhline(10, color='r', ls='--', label='Frozen 10 mm² width')
a4.axvline(25, color='gray', ls=':', label='Frozen 25% measurement budget')
a4.set(xlabel='Sufficient measured slices (%)', ylabel='Before measurement: min-area interval (mm²)', title='Volume alone leaves 48–213 mm² ambiguity')
a4.legend(fontsize=8)
a5 = ax[1, 2]
rr = r[3]['comparisons']
a5.bar(xx - b / 2, [s['predicted_area_mm2'] for s in rr], b, label='Slopes transferred from other group')
a5.bar(xx + b / 2, [s['observed_area_mm2'] for s in rr], b, label='Held-out published area')
a5.set(xticks=xx, xticklabels=['Responders', 'Nonresponders'], ylabel='Minimum area (mm²)', title='Shared tissue-response transfer fails')
a5.legend(fontsize=8)
fig.suptitle('X17 airway geometry: measurement gained; individual MAD prediction remains UNKNOWN', fontsize=14)
fig.savefig(ROOT / 'figure.png', dpi=160)
fig.savefig(ROOT / 'figure.pdf')
plt.close(fig)
