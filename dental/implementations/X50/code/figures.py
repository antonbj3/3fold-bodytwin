from dental_release.paths import expand as _release_expand
from common import *
import sys, zipfile, time
sys.path.insert(0, str(X7 / 'code/vendor_plotting'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sample_stl(b):
    n = int.from_bytes(b[80:84], 'little')
    assert len(b) == 84 + 50 * n
    dt = np.dtype([('normal', '<f4', (3,)), ('vertices', '<f4', (3, 3)), ('attr', '<u2')])
    tri = np.frombuffer(b, dtype=dt, offset=84, count=n)['vertices']
    ids = np.linspace(0, n - 1, min(n, 8000), dtype=int)
    return np.asarray(tri[ids, 0, :], float)

def main():
    st = time.perf_counter()
    (P / 'figures').mkdir(exist_ok=True)
    rs = {r: read(P / f'raw/RESULTS_{r}.json') for r in ['R1', 'R2', 'R3']}
    colors = ['#3877aa', '#16826b', '#b26b26']
    (fig, axes) = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    for (j, sp) in enumerate(['0.9', '0.95']):
        for (k, rd) in enumerate(rs):
            vals = [rs[rd]['findings'][f]['methods']['candidate'][sp] for f in FINDINGS]
            x = np.arange(6) + (k - 1) * 0.18
            for (i, a) in enumerate(vals):
                e = a['empirical_test_roc']
                ci = e.get('bootstrap95')
                p = e.get('sensitivity')
                if ci is not None:
                    axes[0, j].errorbar(x[i], p, yerr=[[max(0, p - ci[0])], [max(0, ci[1] - p)]], fmt='o', ms=5, color=colors[k], capsize=3, label=rd if i == 0 else None)
                e = a['fixed_threshold_test']
                ci = e['sensitivity_wilson95']
                p = e['sensitivity']
                if ci is not None:
                    axes[1, j].errorbar(x[i], p, yerr=[[p - ci[0]], [ci[1] - p]], fmt='o', ms=5, color=colors[k], capsize=3, label=rd if i == 0 else None)
            axes[0, j].axhline(1 - float(sp), color='gray', ls=':', label='Oinformativ slumpordning' if k == 0 else '_nolegend_')
        for i in [0, 1]:
            axes[i, j].set_xticks(range(6), ['Korsbett', "Open\nanteriort", 'Djupt', 'Angle II', 'Angle III', 'Saxbett'])
            axes[i, j].set(ylim=(-0.02, 1.04), ylabel="Sensitivity, proportion of positive findings")
            axes[i, j].grid(axis='y', alpha=0.2)
        axes[0, j].set_title(f'Beskrivande test-ROC vid specificitet ≥ {float(sp):.2f}\n2000 patient-boot stairways; threshold is re-set in each move')
        axes[1, j].set_title(f'Threshold selected for specificity calibration {float(sp):.2f}\nThe specificity of the test may fall below the target; Wilson 95%')
        axes[0, j].legend(fontsize=8)
    fig.suptitle("X50 · fyndvis screening against Bite2Text-reports\nR1 : frozen arc model · R2 : dental dimensions · R3 : consistent training reports\nOpen bite 15 positive, Angle III 19 , Sax bet 9 : under frozen minimum requirements. Previously evaluated cohort.", fontsize=12)
    fig.savefig(P / 'figures/screening_sensitivity.png', dpi=150)
    fig.savefig(P / 'figures/screening_sensitivity.pdf')
    plt.close(fig)
    errors = read(P / 'raw/ERRORS_R2.json')
    selection = []
    for f in FINDINGS:
        candidates = [r for r in errors if r['finding'] == f and r['target_specificity'] == 0.95 and (r['error'] == 'FN')]
        reason = 'FN at cal-frozen Sp95 threshold'
        if not candidates:
            candidates = [r for r in errors if r['finding'] == f and r['target_specificity'] == 0.95]
            reason = 'No FN in evaluated positive cases; show FP at cal-frozen Sp95 threshold'
        if candidates:
            r = min(candidates, key=lambda x: hashlib.sha256(('X50-fixed-error:' + x['case_id']).encode()).hexdigest())
            selection.append(dict(**r, selection_rule=reason))
    write(P / 'raw/ERROR_PANEL.json', selection)
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    (fig, axes) = plt.subplots(3, 2, figsize=(15, 14), constrained_layout=True)
    with zipfile.ZipFile(m['zip']) as z:
        for (ax, r) in zip(axes.flat, selection):
            c = r['case_id']
            row = read(X21 / 'raw/cases' / (c + '.json'))
            R = np.array(row['frame']['R_columns'])
            center = np.array(row['frame']['center_xyz_mm'])
            source = []
            for (jaw, color) in [('upper', '#237bba'), ('lower', '#dd8433')]:
                src = row['source_zip_members'][jaw]
                blob = z.read(src['member'])
                assert hashlib.sha256(blob).hexdigest() == src['sha256']
                pts = (sample_stl(blob) - center) @ R
                ax.scatter(pts[:, 1], pts[:, 2], s=0.9, c=color, alpha=0.15, rasterized=True, label=jaw)
                source.append(src)
            ax.set(aspect='equal', xlabel='Antagen anterioraxel (mm)', ylabel='Antagen superioraxel (mm)', title=f"{NAMES[FINDINGS.index(r['finding'])]} · {c} · {r['error']}\nscore {r['score']:.3f}, threshold {r['threshold']:.3f}; rapporter {r['report_values']}")
            ax.legend(fontsize=8)
            r['mesh_sources'] = source
    fig.suptitle("Pre-determined selection: hash-first error rate per finding, R2 at frozen 0,95 - goal\nProjected mesh points, no patient photos. Physical pose, anatomical landmarks and cause of error UNKNOWN .", fontsize=12)
    fig.savefig(P / 'figures/error_cases.png', dpi=150)
    fig.savefig(P / 'figures/error_cases.pdf')
    plt.close(fig)
    write(P / 'raw/ERROR_PANEL.json', selection)
    rd = read(P / 'raw/RESULTS_R4.json')
    (fig, axes) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    labels = ['All-regionmodell R4', "Anterior-modell R2\non the same new reference", 'Logistisk, samma\nregionala information']
    for (j, sp) in enumerate(['0.9', '0.95']):
        for (i, meth) in enumerate(['candidate', 'anterior_scope_ablation', 'logistic_control']):
            q = rd['methods'][meth][sp]
            e = q['empirical_test_roc']
            v = e['sensitivity']
            ci = e['bootstrap95']
            axes[j].errorbar(i, v, yerr=[[v - ci[0]], [ci[1] - v]], fmt='o', ms=6, capsize=4, color=colors[i])
            axes[j].text(i, 0.05, f"Frozen cal threshold:\nSe {q['fixed_threshold_test']['sensitivity']:.3f}\nSp {q['fixed_threshold_test']['specificity']:.3f}", ha='center', fontsize=9)
        axes[j].set_xticks(range(3), labels, fontsize=9)
        axes[j].set(xlim=(-0.5, 2.5), ylim=(0, 1.05), ylabel="Sensitivity: test-ROC med bootstrap95", title='Specificitet ≥' + sp)
        axes[j].grid(axis='y', alpha=0.2)
    fig.suptitle("Open bite in any reported region · R4 : s explicit observation expansion\n194 test case, 21 positive; 6 changed the outcome case. The same new reference for all three scores.\nAnteriora R1 – R3 - result preserved. Previously used cohort, no screening gates pass.", fontsize=11)
    fig.savefig(P / 'figures/scope_extension.png', dpi=150)
    fig.savefig(P / 'figures/scope_extension.pdf')
    plt.close(fig)
    missed = [r for r in read(P / 'raw/ERRORS_R4.json') if r['error'] == 'FN' and r['target_specificity'] == 0.95]
    (fig, axes) = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    with zipfile.ZipFile(m['zip']) as z:
        for (ax, r) in zip(axes.flat, missed):
            c = r['case_id']
            row = read(X21 / 'raw/cases' / (c + '.json'))
            R = np.array(row['frame']['R_columns'])
            center = np.array(row['frame']['center_xyz_mm'])
            for (jaw, color) in [('upper', '#237bba'), ('lower', '#dd8433')]:
                src = row['source_zip_members'][jaw]
                b = z.read(src['member'])
                assert hashlib.sha256(b).hexdigest() == src['sha256']
                pts = (sample_stl(b) - center) @ R
                ax.scatter(pts[:, 1], pts[:, 2], s=0.9, c=color, alpha=0.18, rasterized=True)
            ax.set(aspect='equal', xlabel='Antagen anterioraxel (mm)', ylabel='Antagen superioraxel (mm)', title=f"{c}: laterally reported open bet missed\nscore {r['score']:.3f} < frozen threshold {r['threshold']:.3f};\nanterior etikett {r['source_reports'][0]['anterior_value']}")
    fig.suptitle(_release_expand("All four misses R4 at frozen 0,95 - goal are lateral source statements\n@DENTAL_CASE_ID@ combines anterior deep bite with lateral open bite. Physical /anatomical error cause UNKNOWN ."), fontsize=11)
    fig.savefig(P / 'figures/lateral_open_misses.png', dpi=150)
    fig.savefig(P / 'figures/lateral_open_misses.pdf')
    plt.close(fig)
    write(P / 'raw/FIGURE_COST.json', dict(wall_s=time.perf_counter() - st, mesh_decode_streamed=True, points_per_arch_at_most=8000, cases=len(selection), plotting_dependency=str(X7 / 'code/vendor_plotting')))
if __name__ == '__main__':
    main()
