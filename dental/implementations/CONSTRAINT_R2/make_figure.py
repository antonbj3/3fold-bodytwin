"""Standalone research figure using the installed, isolated system matplotlib."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

def main():
    n = read(HERE / 'CONSTRAINT_NET_DENTAL_R2.json')
    v = read(HERE / 'VALIDATION.json')
    s = read(HERE / 'CONSTRAINT_STRESS_MAP_DENTAL_R2.json')
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    fig.suptitle('Dental R2: sourced residuals and the limits of scalar stress', fontsize=16)
    ax = axs[0, 0]
    vals = [v['numeric_open_unknown_edges'], v['missing_numeric_open_unknown_edges']]
    ax.barh(['Scoped numeric residual', 'Unidentified — strict FAIL'], vals, color=['#207d70', '#b64941'])
    for (i, x) in enumerate(vals):
        ax.text(x + 0.5, i, str(x), va='center')
    ax.set_xlim(0, max(vals) * 1.18)
    ax.set_xlabel('OPEN / UNKNOWN edges (structural count)')
    ax.set_title('Completeness, not physical validation')
    ax = axs[0, 1]
    top = [r for r in s['all_variables_sorted'] if r['known_normalized_stress'] is not None][:5]
    ax.barh([r['variable'] for r in top][::-1], [r['known_normalized_stress'] for r in top][::-1], color='#51799d')
    ax.set_xlabel('Known sum: residual / source tolerance × decisions')
    ax.set_title('Partial display order — full rank UNKNOWN')
    ax = axs[1, 0]
    ax.bar([-0.18, 0.82], [0, 7], width=0.35, label='Completion A', color='#207d70')
    ax.bar([0.18, 1.18], [7, 0], width=0.35, label='Completion B', color='#ca934c')
    ax.set_xticks([0, 1], ['FDI18', 'FDI28'])
    ax.set_ylabel('Unobserved relative share (percentage points)')
    ax.set_ylim(0, 8.4)
    ax.legend(frameon=False)
    ax.set_title('Identical 7 pp remainder; different named tooth\nAlgebraic witnesses only; no measured allocation', fontsize=11)
    ax = axs[1, 1]
    ax.axis('off')
    ax.text(0, 1, 'Reviewed X71 pilot', fontsize=14, va='top', weight='bold')
    ax.text(0, 0.84, 'M01: registered as-built metrology\nM02: dry / cemented state, same crowns\nM04: fracture pilot + origin\n\n12 shared specimens; 5 evaluable questions\n18–30 assumed operator hours\n0 new physical decisions\n\nCost is phenomenological.\nGap sum alone does not select this plan.', va='top', fontsize=11, linespacing=1.5)
    for suffix in ['png', 'pdf']:
        fig.savefig(HERE / ('CONSTRAINT_NET_R2.' + suffix), dpi=160)
    plt.close(fig)
    print('Figure written: CONSTRAINT_NET_R2.png and .pdf')
if __name__ == '__main__':
    main()
