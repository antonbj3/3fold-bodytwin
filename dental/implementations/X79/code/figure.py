from dental_release.paths import expand as _release_expand
from common import *
import sys
sys.path.insert(0, str(X7 / 'code/vendor_plotting'))
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    st = time.perf_counter()
    cpu = time.process_time()
    r1 = read('raw/R1_RESULTS.json')
    r3 = read('raw/R3_RESULTS.json')
    r2 = read('raw/R2_RESULTS.json')
    r4 = read('raw/R4_RESULTS.json')
    obs = read('raw/REPORT_OBSERVABILITY.json')
    (fig, axes) = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    ax = axes[0]
    for (i, (label, r, methods)) in enumerate([('Restoration', r1, ['shape', 'tooth_position_prior']), ('Absent tooth site', r3, ['shape', 'any_vertex_rule'])]):
        for (j, name) in enumerate(methods):
            mm = r['metrics'][name]
            ci = mm['bootstrap']
            x = i * 3 + j
            vals = [mm['sensitivity'], mm['specificity']]
            ax.bar([x - 0.16, x + 0.16], vals, width=0.3, color=['#3379a8', '#ed9c42'], alpha=1 if j == 0 else 0.55)
            for (xx, y, ival) in zip([x - 0.16, x + 0.16], vals, [ci['sensitivity95'], ci['specificity95']]):
                ax.errorbar(xx, y, yerr=[[max(0, y - ival[0])], [max(0, ival[1] - y)]], color='black', capsize=2, fmt='none')
    ax.set_xticks([0, 1, 3, 4], ['Restoration\nshape', 'Position prior', 'Absence\nshape', 'Any FDI vertex'])
    ax.tick_params(axis='x', rotation=20)
    ax.set(ylim=(0, 1.06), ylabel='Reported-finding agreement, fraction', title='Patient-disjoint test, fixed calibration threshold')
    ax.text(0.02, 0.03, 'Blue: sensitivity; orange: specificity\nPatient-cluster descriptive 95% intervals', transform=ax.transAxes, fontsize=8)
    ax.grid(axis='y', alpha=0.2)
    ax = axes[1]
    counts = [obs['report_observability'][k]['positive_named_patients'] for k in ['ios', 'intraoral-photo']]
    bars = ax.bar(['IOS text', 'Photo text'], counts, color=['#aaaaaa', '#438574'])
    ax.bar_label(bars, padding=3)
    ax.set(ylabel='Patients with named treatment statements (any report)', title='Previously unused photo-report information', ylim=(0, 335))
    ax.text(0.03, 0.94, f"{r2['restored_teeth']} first-report teeth / {r2['restored_patients']} patients\nExact material/interface still UNKNOWN", transform=ax.transAxes, va='top', fontsize=9)
    ax = axes[2]
    ax.axis('off')
    lines = ['What the scan does not establish', '', f'Restoration: R1 numeric lift failed; only 13', 'explicitly negative test patients.', '', f"Spacing: {len(read('raw/R4_SPACING_GEOMETRY.json'))} positive arches; no negatives.", f"Positive-only point-set agreement {100 * r4['positive_only_agreement']:.1f}%.", 'Specificity UNKNOWN; no surface enclosure.', '', f"Wear: {obs['wear']['test_patients']} reported test patients.", 'Unknown unworn reference and rater identity.', _release_expand('One @DENTAL_CASE_ID@ native/site scope error corrected.')]
    ax.text(0, 1, '\n'.join(lines), va='top', fontsize=10)
    fig.savefig(ROOT / 'figures/report_geometry_findings.png', dpi=160)
    plt.close(fig)
    write('raw/FIGURE_COST.json', cost(st, cpu))
if __name__ == '__main__':
    main()
