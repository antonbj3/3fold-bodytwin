"""Run with python3 -s: standard system Matplotlib environment for standalone figures."""
import os, json, pathlib
if __name__ == '__main__':
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[k] = '1'
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    P = pathlib.Path(os.environ.get('X3_RUN_DIR', pathlib.Path(__file__).resolve().parent))
    rounds = {r: json.load(open(P / f'RESULTS_{r}.json')) for r in ['R1', 'R2'] if (P / f'RESULTS_{r}.json').exists()}
    (fig, axes) = plt.subplots(len(rounds), 2, figsize=(12, 4.4 * len(rounds)), squeeze=False, constrained_layout=True)
    types = ['incisor', 'canine', 'premolar', 'molar']
    colors = {'mirror': '#68737c', 'population_mean_sdf': '#d99524', 'neighbor_mesh': '#527dc3', 'neighbor_sdf': '#95b5e8', 'ring_mesh': '#278b71', 'ring_sdf': '#71b7a7'}
    for (ridx, (rn, res)) in enumerate(rounds.items()):
        methods = ['mirror', 'population_mean_sdf', 'neighbor_mesh', 'neighbor_sdf'] if rn == 'R1' else ['mirror', 'population_mean_sdf', 'neighbor_mesh', 'ring_mesh', 'ring_sdf']
        for (col, (metric, title)) in enumerate([('occlusal_p95_mm', 'Cusp / incisal surface p95'), ('proximal_patch_p95_mm', 'Proximal patch distance error p95')]):
            ax = axes[ridx, col]
            width = 0.8 / len(methods)
            x = np.arange(4)
            for (mindex, m) in enumerate(methods):
                vals = []
                for typ in types:
                    rr = next((t for t in res['tables'] if t['jaw'] == 'all' and t['tooth_type'] == typ and (t['method'] == m)), None)
                    vals.append(rr[metric] if rr else np.nan)
                ax.bar(x - 0.4 + width / 2 + mindex * width, vals, width, label=m, color=colors[m])
            ax.set_xticks(x)
            ax.set_xticklabels(types)
            ax.set_ylabel('Patient-balanced mean [provisional mm]')
            ax.set_title(f'{rn} | {title}')
            ax.grid(axis='y', alpha=0.2)
            ax.set_axisbelow(True)
            if col == 0:
                ax.axhline(0.5, color='#a83e3e', linestyle='--', linewidth=1, label='frozen0.50mm research target')
                ax.legend(fontsize=8, loc='upper left', ncol=2)
    fig.suptitle('Real Teeth3DS holdout geometry | Antagonist contact UNKNOWN', fontsize=14)
    fig.savefig(P / 'benchmark_figure.png', dpi=200)
    fig.savefig(P / 'benchmark_figure.pdf')
    print('FIGURE_WRITTEN')
