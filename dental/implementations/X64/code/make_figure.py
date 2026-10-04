import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
from pathlib import Path
import json, resource, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from freeze import ROOT, dump
OUT = Path(os.environ.get('X64_REPLAY_DIR', str(ROOT)))
(OUT / 'raw').mkdir(parents=True, exist_ok=True)

def main():
    tic = time.perf_counter()
    rays = np.load(ROOT / 'raw/OPTICAL_RAYS.npz')
    vis = rays['visible']
    pts = rays['points'][vis]
    h = rays['thickness_mm'][vis]
    r1 = json.loads((ROOT / 'raw/RESULTS_R1.json').read_text())
    r2 = json.loads((ROOT / 'raw/RESULTS_R2.json').read_text())
    source = json.loads((ROOT / 'raw/SOURCE_ROWS.json').read_text())
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8.5), layout='constrained')
    ax = axs[0, 0]
    sc = ax.scatter(pts[:, 0], pts[:, 2], c=h, cmap='viridis', s=24, vmin=0.5, vmax=3.0, edgecolor='none')
    fig.colorbar(sc, ax=ax, label='Normal material chord (mm), clipped at 3')
    ax.set(xlabel='D1 local x (mm)', ylabel='D1 local z (mm)', title='Same D1 anatomy: exposed surface samples')
    ax.text(0.02, 0.03, 'Ray length is a surrogate optical path.\nOnly 33/236 paths in source thickness range.', transform=ax.transAxes, fontsize=9, bbox=dict(facecolor='white', alpha=0.8))
    ax = axs[0, 1]
    colors = ['#397f68', '#cc5858', '#d9ab4c', '#a9b0bb']
    labels = ['Box pass', 'Box fail', 'Unresolved', 'Outside source support']
    bottom = np.zeros(2)
    for (label, color, key) in zip(labels, colors, ['MODEL_BOX_PASS', 'MODEL_BOX_FAIL', 'MODEL_BOX_UNRESOLVED', 'OUTSIDE_SOURCE_SUPPORT']):
        values = [r['counts'][key] for r in r2['TM02']]
        ax.bar(['Standard 420 min', 'Speed 105 min'], values, bottom=bottom, label=label, color=color)
        bottom += values
    ax.set(ylabel='Exposed sample sites (count)', title='TM02: optical source boxes + declared ±0.12 mm')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_ylim(0, 285)
    ax = axs[1, 0]
    times = [r['cure_min'] for r in source['cure']]
    mu = [r['mean_MPa'] for r in source['cure']]
    sd = [2 * r['SD_MPa'] for r in source['cure']]
    ax.errorbar(times, mu, yerr=sd, fmt='o-', color='#374d73', capsize=5, label='Source mean ±2 SD (scenario)')
    for (r, c) in zip(r2['TM05'], ['#397f68', '#cc5858']):
        ax.axhline(r['peak_demand_MPa'], color=c, ls='--', label=f"D1 {r['loadcase']}: {r['peak_demand_MPa']:.2f} MPa")
    ax.set(xlabel='Post-cure time (min)', ylabel='MPa', title='TM05: source bar strength vs frozen D1 demand, 150 N')
    ax.legend(fontsize=8, loc='lower right')
    ax = axs[1, 1]
    w = r1['sufficiency']
    x = np.arange(2)
    width = 0.35
    a = [w['optical']['downstream_min_visible_TP_A'] / 12, w['strength']['downstream_min_ratio_A']]
    b = [w['optical']['downstream_min_visible_TP_B'] / 12, w['strength']['downstream_min_ratio_B']]
    ax.bar(x - width / 2, a, width, label='State A', color='#d9ab4c')
    ax.bar(x + width / 2, b, width, label='State B', color='#374d73')
    ax.axhline(1, color='#cc5858', ls='--', label='Fixed query threshold')
    ax.set_xticks(x, ['Optical visible minimum /12', 'Strength /local demand'])
    ax.set(ylabel='Dimensionless local decision ratio', title='Identical mean + histogram, changed spatial association')
    ax.legend(fontsize=8)
    fig.suptitle('Two process chains on the same STS D1 crown\nConditional model queries; physical crown validation UNKNOWN', fontsize=15)
    fig.savefig(OUT / 'manufacturing_demo.png', dpi=160)
    fig.savefig(OUT / 'manufacturing_demo.pdf')
    dump(OUT / 'raw/FIGURE_COST.json', {'seconds': time.perf_counter() - tic, 'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
if __name__ == '__main__':
    main()
