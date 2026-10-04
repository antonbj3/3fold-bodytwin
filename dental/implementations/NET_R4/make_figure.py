from common import *

def make():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from collections import Counter
    n = read('CONSTRAINT_NET_DENTAL_R4.json')
    s = read('SUFFICIENCY.json')
    r = Counter((v['plausible_range']['status'] for v in n['variables'].values()))
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    ax = axs[0, 0]
    ax.bar(['1st', '10th'], [329.9, 253.7], yerr=[33.3, 36.8], color=['#17697d', '#76b7b2'], capsize=5)
    ax.set(ylabel='Axial preload, N (group mean ± SD)', title='Sagheb: source group, n = 25')
    ax.text(0.04, 0.91, '25 Ncm; dry; fixture separated 0.10 mm', transform=ax.transAxes, fontsize=9)
    ax.set_ylim(0, 430)
    ax = axs[0, 1]
    ax.bar(['Replica', 'Dry microCT'], [187.35, 126.53], yerr=[53.89, 18.12], capsize=5, color=['#d88436', '#e8bf8d'])
    ax.set(ylabel='Internal mid-occlusal gap, µm', title='Cunali: Amann group, n = 10')
    ax.text(0.03, 0.91, 'Group means ± SD; before cementation', transform=ax.transAxes, fontsize=9)
    ax.set_ylim(0, 300)
    ax = axs[1, 0]
    cats = ['SCOPED_SOURCE_ENVELOPE', 'MATHEMATICAL_DOMAIN', 'PARTIAL_DOMAIN', 'UNKNOWN', 'NOT_SCALAR']
    ax.barh(range(5), [r[c] for c in cats], color=['#17697d', '#76b7b2', '#d88436', '#b14548', '#808b96'])
    ax.set_yticks(range(5), ['Source context', 'Exact domain only', 'Upper bound missing', 'Range unknown', 'Structured value'])
    ax.invert_yaxis()
    ax.set(xlabel='Variables (100 = 99 inherited + 1 distinct gap port)', title='Range coverage: no universal bounds inferred')
    ax = axs[1, 1]
    ax.bar(['Pairing A', 'Pairing B'], [s['pairing']['minimum_individual_loss_A_N'], s['pairing']['minimum_individual_loss_B_N']], color=['#7560a5', '#baa6d6'])
    ax.set(ylabel='Minimum individual loss, N', title='Synthetic exact-summary counterexample')
    ax.text(0.02, 0.91, 'Identical marginal moments + ranges\nIdentity error = 0; downstream difference = 20 N', transform=ax.transAxes, fontsize=9)
    ax.set_ylim(0, 105)
    fig.suptitle('R4: audited reference facets and scope-dependent range admission\nNo new physical measurements or target-system validation', fontsize=13)
    fig.savefig(H / 'CONSTRAINT_NET_R4.png', dpi=160)
    fig.savefig(H / 'CONSTRAINT_NET_R4.pdf')
    plt.close(fig)
if __name__ == '__main__':
    make()
