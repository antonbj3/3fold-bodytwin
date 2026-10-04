import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT
from patient import load

def run():
    (points, vertical, rows, meshes) = load()
    r = json.loads((ROOT / 'raw/R1.json').read_text())
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, axs) = plt.subplots(2, 2, figsize=(11, 8), layout='constrained')
    ax = axs[0, 0]
    for (row, m) in zip(rows, meshes):
        v = m.vertices[::16]
        ax.scatter(v[:, 0], v[:, 1], s=0.3, c='#8799a7', alpha=0.35)
        p = np.array(row['lingual_surface_anchor_mm'])
        ax.text(p[0], p[1] - 0.4, str(row['FDI']), ha='center')
    ax.plot(points[:, 0], points[:, 1], '-o', c='#0d8068', lw=2, label='proposed six-bond retainer')
    ax.set_aspect('equal')
    ax.set(xlabel='source x (mm)', ylabel='source y (mm)', title='Same model: Open-Full-Jaw Patient 2')
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    tests = r['O18']['tests']
    x = np.arange(2)
    ax.bar(x - 0.18, [t['observed_N'] for t in tests], 0.36, label='published measured force', color='#29455c')
    ax.bar(x + 0.18, [t['predicted_N'] for t in tests], 0.36, label='nominal beam prediction', color='#c46640')
    ax.set_xticks(x, ['AP at 0.06 mm', 'Vertical at 0.02 mm'])
    ax.set(ylabel='force (N)', title='External fixture test: vertical closure fails')
    ax.legend(fontsize=8)
    for (i, t) in enumerate(tests):
        ax.text(i, max(t['observed_N'], t['predicted_N']) + 0.04, f"{100 * t['relative_error']:.1f}% error", ha='center')
    ax = axs[1, 0]
    pair = r['O18']['sufficiency']
    ax.bar([0, 1], pair['peak_reactions_N'], color=['#29455c', '#c46640'])
    ax.set_xticks([0, 1], ['common translation', 'alternating signed errors'])
    ax.set(ylabel='conditional peak model reaction (N)', title='Identical RMS, different deformation state')
    ax.text(0.04, 0.91, 'Both RMS = 0.015625 mm\nIdentity error = 0\nPatient forces are unvalidated', transform=ax.transAxes, va='top')
    ax = axs[1, 1]
    ax.bar([0, 1], [100, 93.75], color=['#29455c', '#c46640'])
    ax.axhline(95, c='#0d8068', ls='--', label='chosen lab release target')
    ax.set_xticks([0, 1], ['A: +/- 0.25 mm', 'B: rare +/- 1 mm'])
    ax.set(ylim=(80, 102), ylabel='within 0.5 mm (%)', title='Same mean and SD, different tolerance tails')
    ax.text(0.03, 0.16, 'Mean = 0; SD = 0.25 mm in both\nExact dyadic witness, not patient failures', transform=ax.transAxes)
    ax.legend(fontsize=8, loc='lower right')
    target = ROOT / 'figures'
    target.mkdir(exist_ok=True)
    fig.savefig(target / 'ortho_demo.png', dpi=180)
    fig.savefig(target / 'ortho_demo.pdf')
    plt.close(fig)
if __name__ == '__main__':
    run()
