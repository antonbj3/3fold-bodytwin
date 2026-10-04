import json
from pathlib import Path
import numpy as np
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    P = Path(__file__).resolve().parents[1]
    rows = []
    for path in sorted((P / 'raw/refined_certificates').glob('*.json')):
        r = json.load(open(path))
        old = json.load(open(P / 'raw/certificates' / path.name))
        o = {(p['upper_fdi'], p['lower_fdi']): p for p in old['pairs']}
        for p in r['pairs']:
            rows.append((o[p['upper_fdi'], p['lower_fdi']]['interval_width_mm'], p['width_mm']))
    a = np.array(rows)
    (fig, ax) = plt.subplots(figsize=(8, 5.2))
    fig.subplots_adjust(bottom=0.2, top=0.9)
    ax.scatter(a[:, 0], a[:, 1], s=10, alpha=0.5, color='#477c99')
    ax.set(xscale='log', yscale='log', xlabel='Width before exact refinement , R3 ( mm )', ylabel='Width after exact refinement , R4 ( mm )', title='440 tooth pairs in 12 default cases: conditionally projected gap certificate')
    ax.axhline(0.005, ls='--', color='#b04433', label='Frozen 5 µm - boundary')
    ax.legend()
    fig.text(0.5, 0.03, 'Mathematical certificate for rounded mesh.\nPhysical contact, scanner pose and dental identity UNKNOWN.', ha='center', fontsize=9)
    fig.savefig(P / 'figures/certificate_width.png', dpi=160)
    plt.close(fig)
