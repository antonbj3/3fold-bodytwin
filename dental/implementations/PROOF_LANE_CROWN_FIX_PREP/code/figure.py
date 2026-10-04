from core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def run():
    rr = read(R / 'RESULTS_R3.json')['rows']
    fields = ['form', 'wall', 'nominal_gap', 'insertion', 'native_containment', 'margin', 'milling', 'STL']
    names = ['Form', 'Wall', 'Gap', 'Insertion¹', 'Native²', 'Marginal', 'Milling', 'Export']
    casekeys = sorted({r['case_key'] for r in rr})
    a = []
    labels = []
    for r in rr:
        a.append([1 if r.get('gates', {}).get(f) == 'PASS' else 0 if r.get('gates', {}).get(f) == 'FAIL' else 0.5 for f in fields])
        labels.append(f"{chr(65 + casekeys.index(r['case_key']))} / {r['fdi']} / {('Zr' if r['material'] == 'KATANA_ML' else 'LS2')}")
    fig = plt.figure(figsize=(14, 12), layout='constrained')
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1])
    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(a, cmap=ListedColormap(['#c66555', '#d5d9dc', '#438b77']), vmin=0, vmax=1, aspect='auto')
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_xticks(range(len(names)), names, rotation=30, ha='right', fontsize=9)
    ax.set_title('All 18 teeth × two material\nGreen: local PASS · Red: FAIL · Grey: UNKNOWN', fontsize=12)
    ax.set_ylabel('Jaw / FDI / material')
    ax.text(0.01, -0.08, '¹ Isolated model, not whole jaw .\n ²  Full Q⊕ B_w - enclosure  i sluten tandmodell.', transform=ax.transAxes, fontsize=9, va='top')
    ar = fig.add_subplot(gs[0, 1])
    row = next((r for r in rr if r['status'] == 'GENERATED' and r.get('native_containment_obstruction')))
    m = load(row['mesh_path'])
    M = m['M']
    t = sample(m['target'], 5000) @ M.T
    outer = m['vertices'][m['faces'][m['roles'] == 0]]
    s = sample(outer, 5000) @ M.T
    q = m['prep_vertices'] @ M.T
    ar.scatter(t[:, 0], t[:, 2], s=2, color='#aaa9a2', alpha=0.4, label='Original IOS')
    ar.scatter(s[:, 0], s[:, 2], s=1, color='#387c91', alpha=0.5, label='Mirrored crown shape')
    ar.scatter(q[:, 0], q[:, 2], s=2, color='#d49346', alpha=0.6, label='Taperad preparation')
    p = np.array(row['native_containment_obstruction']['point_mm']) @ M.T
    ar.plot(p[0], p[2], 'r*', ms=13, label='Precise point in Q, outside T')
    ar.set_aspect('equal')
    ar.set_xlabel('Local x ( mm )')
    ar.set_ylabel('Inserting axle z ( mm )')
    ar.set_title('Local wall not enough for full preparation')
    ar.legend(fontsize=8, loc='upper right')
    ar = fig.add_subplot(gs[1, 1])
    valid = [r for r in rr if r.get('shape')]
    limits = read(R / 'PREREG_A.json')['metrics']['shape_p95_mm']
    colors = {'anterior': '#ab5e69', 'premolar': '#438b77', 'molar': '#4d7da0'}
    for fam in colors:
        rows = [r for r in valid if r['family'] == fam]
        idx = [valid.index(r) for r in rows]
        val = [r['shape']['p95_mm'] / limits[fam] for r in rows]
        ar.scatter(idx, val, color=colors[fam], label=fam, s=30)
    ar.axhline(1, color='#333', ls='--', lw=1, label='Frozen form boundary')
    ar.set_xlabel('Generated tooth / material candidate')
    ar.set_ylabel('Shape error p95 / frozen tolerance')
    ar.set_title('Form against original external scan')
    ar.legend(fontsize=8)
    ar.grid(alpha=0.2)
    fig.suptitle('New Preparations: local proof, still 0/18 complete crowns', fontsize=16)
    fig.savefig(R / 'figures/DEMO.png', dpi=160)
    fig.savefig(R / 'figures/DEMO.svg')
    plt.close(fig)
if __name__ == '__main__':
    run()
