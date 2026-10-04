"""Standalone publication figure, system matplotlib/NumPy ABI via python3 -s."""
from dental_release.paths import expand as _release_expand
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H = Path(__file__).resolve().parents[1]
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X18b_prep_surface'))

def main():
    rounds = {t: json.loads((H / f'rounds/{t}.json').read_text()) for t in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']}
    (r1, r2, r3, r4, r5, r6) = [rounds[t] for t in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']]
    plt.rcParams.update({'font.size': 10})
    (fig, axs) = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    ex = next((r for r in r1['rows'] if r['case'] == 127 and r['fdi'] == 36))
    z = np.load(ex['array_file']['path'])
    xy = z['xy']
    good = np.isfinite(z['source_z'])
    mask = z['0.1_reference_contact']
    axs[0, 0].scatter(xy[good, 0], xy[good, 1], c=z['source_z'][good], s=4, cmap='Greys')
    axs[0, 0].scatter(xy[mask, 0], xy[mask, 1], s=6, c='crimson')
    axs[0, 0].set(title='Measured exterior + proximity regions\ncase127 / FDI36; .15mm samples', xlabel='x [mm]', ylabel='y [mm]')
    axs[0, 0].set_aspect('equal')
    rr = r6['rows']
    px = [r['practice_errors']['IoU'] for r in rr]
    cx = [r['candidate_errors']['IoU'] for r in rr]
    axs[0, 1].plot(px, 'o-', label='Generic donor proxy')
    axs[0, 1].plot(cx, 's-', label='Measured source retained')
    axs[0, 1].set(title='Same historical sample support\n12 tooth replacements', xlabel='Tooth case index', ylabel='Proximity-mask IoU', ylim=(-0.05, 1.05))
    axs[0, 1].legend()
    field = np.load(r4['prediction']['field_arrays']['path'])
    ix = field['phi'].shape[0] // 2
    sl = field['phi'][ix].T
    step = float(field['grid_step_mm'])
    origin = field['origin']
    yy = origin[1] + np.arange(sl.shape[1]) * step
    zz = origin[2] + np.arange(sl.shape[0]) * step
    axs[0, 2].contour(yy, zz, sl, levels=[-r5['depth_mm'], 0], colors=['steelblue', 'black'])
    axs[0, 2].set(title='Implicit virtual preparation section\nclosed generic erosion1.205mm', xlabel='y [mm]', ylabel='z [mm]')
    axs[0, 2].set_aspect('equal')
    hp = r2['source_tables']['Hattori']
    means = np.array(hp['mean_N'])
    sd = np.array(hp['sd_N'])
    lo = np.maximum(0, means - 2 * sd)
    hi = means + 2 * sd
    axs[1, 0].errorbar(np.arange(7), means, yerr=[means - lo, hi - means], fmt='o', capsize=4)
    axs[1, 0].axhline(100, color='gray', ls='--', label='Old100N scenario')
    axs[1, 0].set(xticks=np.arange(7), xticklabels=hp['types'], title='Published mandibular tooth-force scenarios\nmean±2SD; population, not patient bounds', ylabel='Force [N]')
    axs[1, 0].legend()
    old = np.array([r['old_shape_ratio'] for r in r2['stress_results']])
    new = np.array([r['verified_new_shape_ratio'] for r in r2['stress_results']])
    mom = np.array([r['conditional_shape_ratio'] for r in r3['rows']])
    axs[1, 1].plot(old, 'o-', label='Original unknown patch shares')
    axs[1, 1].plot(new, 'x', ms=10, label='With published tooth force')
    axs[1, 1].plot(mom, 's--', label='Hypothetical measured CoP')
    axs[1, 1].set(title='Fixed X18 roof: upper / area tensile peak\nTooth-force-only narrowing =0%', xlabel='Original X18 tooth case index', ylabel='Ratio [1]')
    axs[1, 1].legend(fontsize=8)
    ps = next((r for r in r2['pressure_predictions'] if r['case'] == 127 and r['fdi'] == 36))
    pp = ps['patches']
    up = np.array([p['average_pressure_hi'] for p in pp])
    axs[1, 2].vlines(np.arange(len(pp)), 0, up, color='purple')
    axs[1, 2].scatter(np.arange(len(pp)), up, c='purple')
    axs[1, 2].set(title='Per-region pressure scenario intervals\nUnmeasured fractions; not simultaneous maxima', xlabel='Proximity region, case127 / FDI36', ylabel='Average pressure scenario [MPa]')
    fig.suptitle('Measured surface retention succeeds; physical load/preparation validity remains unresolved', fontsize=14)
    (H / 'figures').mkdir(exist_ok=True)
    fig.savefig(H / 'figures/demo.png', dpi=180)
    fig.savefig(H / 'figures/demo.pdf')
    plt.close(fig)
if __name__ == '__main__':
    main()
