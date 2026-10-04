import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H = Path(__file__).resolve().parent
m = json.loads((H / 'results.json').read_text())['metrics']
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
(f, ax) = plt.subplots(2, 2, figsize=(11, 7.5), layout='constrained')
a = ax[0, 0]
n = m['L01']
a.bar(['Fracture sites', 'Test regions', 'Physical specimens'], [n['characterized_site_rows'], n['unique_test_regions'], n['unique_specimens']], color=['#4169a1', '#528a9d', '#b26945'])
a.set_title('Fatigue: preserve the sampling unit')
a.set_ylabel('Count [POPULATION]')
a.set_ylim(0, 160)
for bar in a.patches:
    a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2, f'{bar.get_height():.0f}', ha='center')
a = ax[0, 1]
v = [d['eigenvalues'][-1] for d in m['L05']['maps']]
a.bar(['Small Gaussian', 'Tailored beam'], v, color=['#4169a1', '#b26945'])
a.axhline(1 / 3, c='black', ls='--', label='Isotropic tensor')
a.set_title('EBSD: area-weighted texture')
a.set_ylabel('Largest tensor eigenvalue [POPULATION]')
a.set_ylim(0, 0.65)
a.legend(frameon=False)
a = ax[1, 0]
v = m['L08']
a.bar(['Largest observed pore', 'Median pore per sample'], [v['median_largest_sphericity'], v['median_sample_median_sphericity']], color=['#b26945', '#4169a1'])
a.axhline(1, c='black', ls='--', label='Spherical assumption')
a.set_ylabel('Median sphericity [POPULATION]')
a.set_ylim(0, 1.15)
a.set_title('XCT: extreme shape differs from average')
a.legend(frameon=False)
a = ax[1, 1]
v = m['L11']['by_topology']
a.bar([r['topology'] for r in v], [100 * r['relative_error'] for r in v], color='#a0a0a0')
a.axhline(10, c='black', ls='--', label='Historical 10% gate')
a.set_ylabel('Relative modulus error (%) [POPULATION]')
a.set_title('Rejected L11: native modulus unit missing')
a.set_xlabel('PLA topology; binary n=1/2, highest infill excluded from C fit')
a.legend(frameon=False)
f.suptitle('PROOF_LANE-infolink3 · local external observations and bounded comparisons', fontsize=14)
f.savefig(H / 'information_links_figure.svg')
f.savefig(H / 'information_links_figure.png', dpi=170)
f.savefig(H / 'information_links_figure.pdf')
