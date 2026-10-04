import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .io import ROOT, dump

def report(r):

    def fmt(v):
        return 'UNKNOWN' if v is None else f'{v:.4f}'
    table = ['| Generator | PASS | FAIL | UNKNOWN | INVALID |', '|---|---:|---:|---:|---:|']
    for (n, s) in r['summary'].items():
        table.append(f"| {n} | {s['PASS']} | {s['FAIL']} | {s['UNKNOWN']} | {s['INVALID']} |")
    table += ['', 'Shape diagnostics (eight anatomical references; material replicas are not independent cases):', '', '| Generator | Material | Mean sampled RMS (mm) | Mean sampled p95 (mm) |', '|---|---|---:|---:|']
    for (n, s) in r['summary'].items():
        for (m, v) in s['by_material'].items():
            table.append(f"| {n} | {m} | {fmt(v['RMS_mm'])} | {fmt(v['p95_mm'])} |")
    (ROOT / 'BASELINE_RESULTS.md').write_text('\n'.join(table) + '\n')
    (fig, axs) = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
    names = list(r['summary'])
    x = np.arange(len(names))
    bottom = np.zeros(len(names))
    labels = [{'contralateral_mirror': 'Mirror', 'population_mean': 'Population', 'parametric_IFU': 'IFU'}.get(n, n.split(':')[0].split('.')[-1]) for n in names]
    for (status, color) in [('PASS', '#257447'), ('FAIL', '#b94432'), ('UNKNOWN', '#c49a34'), ('INVALID', '#6a4985')]:
        y = np.array([r['summary'][n][status] for n in names])
        axs[0].bar(x, y, bottom=bottom, label=status, color=color)
        bottom += y
    axs[0].set_xticks(x, labels)
    axs[0].set_ylabel('Tasks (24 per generator)')
    axs[0].set_title('Required geometric gates')
    axs[0].legend(fontsize=8)
    for (j, (m, label)) in enumerate([('KATANA_ML', 'ML'), ('KATANA_UTML', 'UTML'), ('IPS_e_max_CAD_adhesive_1mm', 'e.max')]):
        axs[1].bar(x + (j - 1) * 0.23, [r['summary'][n]['by_material'][m]['RMS_mm'] if r['summary'][n]['by_material'][m]['RMS_mm'] is not None else np.nan for n in names], width=0.23, label=label)
    axs[1].set_xticks(x, labels)
    axs[1].set_ylabel('Sampled surface RMS (mm)')
    axs[1].set_title('Withheld anatomy; not clinical quality')
    axs[1].legend(fontsize=8)
    cal = json.loads((ROOT / 'raw/calibration_results.json').read_text())
    for (name, color) in [('model', '#3578a3'), ('control', '#ad6944')]:
        rr = [v for v in cal['rows'] if v['method'] == name]
        xx = [v['measured_N'] for v in rr]
        yy = [v['prediction_N'] for v in rr]
        axs[2].scatter(xx, yy, label=name, color=color, s=22, alpha=0.7)
    maxv = max((v['measured_N'] for v in cal['rows'])) * 1.1
    axs[2].plot([0, maxv], [0, maxv], 'k--', linewidth=1)
    axs[2].set(xlabel='Published group mean (N)', ylabel='Held-study prediction (N)', title='Fracture transfer stress test')
    axs[2].legend(fontsize=8)
    fig.suptitle('DentalGenCAD-Bench foundation — digital proofs, measured references, explicit UNKNOWN')
    fig.savefig(ROOT / 'figures/benchmark.png', dpi=160)
    fig.savefig(ROOT / 'figures/benchmark.svg')
    plt.close(fig)
