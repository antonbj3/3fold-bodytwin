from bench import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    r1 = json.loads((ROOT / 'raw/R1.json').read_text())
    r4 = json.loads((ROOT / 'raw/R4.json').read_text())
    r6 = json.loads((ROOT / 'raw/R6.json').read_text())
    (fig, ax) = plt.subplots(1, 3, figsize=(13.5, 4), constrained_layout=True)
    valid = [r for r in r1['rows'] if r.get('numerical_gate')]
    ax[0].bar(range(len(valid)), [r['max_local_response_ambiguity_N'] for r in valid], color='#7c4691')
    ax[0].axhline(1, color='black', ls='--', lw=1)
    ax[0].set(title='Same global K, different 10 µm edit response', xlabel='12 source layouts', ylabel='Largest ambiguity (N)')
    q = [t for r in r4['rows'] for t in r['queries']]
    ax[1].scatter([t['radius_N'] for t in q], [t['l2_force_error_N'] for t in q], s=8, alpha=0.5, color='#217a78')
    ax[1].plot([0, 1], [0, 1], color='black', ls='--', lw=1)
    ax[1].set(xlim=(0, 1), ylim=(0, 1), title='Coupled force response, declared sensor error', xlabel='R4 error radius (N)', ylabel='Actual simulated L2 error (N)')
    q = [t for r in r6['rows'] for t in r['queries']]
    a = np.sort([t['candidate_radius_N'] for t in q])
    b = np.sort([t['strong_control_radius_N'] for t in q])
    ax[2].plot(a, label='Triangle bound', color='#bb6f22')
    ax[2].plot(b, label='Exact corner control', color='#28679c')
    ax[2].axhline(1, color='black', ls='--', lw=1)
    ax[2].set(title='New 30/40 µm edits after contact activation', xlabel='Sorted new queries (1416)', ylabel='Error radius (N)')
    ax[2].legend(frameon=False, fontsize=8)
    fig.suptitle('Conditional mechanical bench • geometry positions from X21 • no physical force measurements', fontsize=11)
    fig.savefig(ROOT / 'figures/response_operator.png', dpi=180)
    fig.savefig(ROOT / 'figures/response_operator.pdf')
    print('Scientific PNG/PDF exported')
if __name__ == '__main__':
    main()
