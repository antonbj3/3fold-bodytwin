"""Standalone figure of actual replay costs and proposed laboratory observations."""
import json
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_replay import HERE, execution_manifest

def main():
    j = json.loads((HERE / 'raw/ALL_REPLAY.json').read_text())
    groups = defaultdict(lambda : [0.0, 0.0, 0, 0])
    for r in j['records']:
        i = 0 if r['status'] == 'PASS' else 1
        groups[r['batch']][i] += r['wall_seconds']
        groups[r['batch']][i + 2] += 1
    labels = sorted(groups, key=lambda x: int(x[5:]) if x.startswith('batch') else 0)
    (fig, (ax, bx)) = plt.subplots(1, 2, figsize=(13, 8), gridspec_kw={'width_ratios': [1.3, 1]})
    y = range(len(labels))
    a = [groups[k][0] for k in labels]
    b = [groups[k][1] for k in labels]
    ax.barh(y, a, color='#417d68', label="PASS : specified runcope")
    ax.barh(y, b, left=a, color='#b34b42', label="Run - or gate error")
    ax.set_yticks(list(y), ['Bas' if k == 'base' else 'Granskning ' + k[5:] for k in labels])
    for (i, k) in enumerate(labels):
        ax.text(a[i] + b[i] + max(a + b + [1]) * 0.01, i, f'{groups[k][2]} pass / {groups[k][3]} fail', va='center', fontsize=8)
    ax.set_xlim(0, max([a[i] + b[i] for i in range(len(a))] + [1]) * 1.45)
    ax.set_xlabel("Measured command time [s ], PHENOMENOLOGICAL")
    expected = len(json.loads(execution_manifest().read_text())['executions'])
    ax.set_title(f"{len(j['records'])}/{expected} active runs completed")
    ax.invert_yaxis()
    ax.legend(loc='lower right', fontsize=8)
    ax.grid(axis='x', alpha=0.2)
    ax.set_axisbelow(True)
    bx.axis('off')
    bx.set_title("Freezer before held out lab measurement")
    texts = [('1  Registrera samma prov', 'CAD, tillverkad form, die, antagonist, ram, batch'), ("2 Measure observation channels", "Height , regional Newton channels, total + moment"), ("3 Calibrate the actual change", "Two regional bases, 5 heights; repeat prehistory"), ('4  Frys prediktionen', 'Geometri-/kalibreringshashar, intervall, kriterier'), ("5 Measure an unused state", "Make/change and compare before the next fit"), ("Unknown before measurement", 'Fysisk kraft, nonlinearitetsrest, sparad justering')]
    for (i, (title, sub)) in enumerate(texts):
        h = 0.9 - i * 0.145
        bx.text(0.03, h, title, transform=bx.transAxes, fontsize=11, fontweight='bold', va='top')
        bx.text(0.03, h - 0.045, sub, transform=bx.transAxes, fontsize=8.5, va='top', wrap=True)
        if i < 4:
            bx.annotate('', xy=(0.5, h - 0.11), xytext=(0.5, h - 0.085), xycoords='axes fraction', arrowprops={'arrowstyle': '->', 'color': '#555'})
    fig.suptitle("Recycle digital operations and try a crown in the lab", fontsize=15)
    fig.text(0.03, 0.015, "PASS applies to specified digital/ source bound gate . The FAIL and UNKNOWN models are maintained. No physically specimen was performed.", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(HERE / 'PACKAGE_AND_LAB.png', dpi=180)
    fig.savefig(HERE / 'PACKAGE_AND_LAB.svg')
    plt.close(fig)
if __name__ == '__main__':
    main()
