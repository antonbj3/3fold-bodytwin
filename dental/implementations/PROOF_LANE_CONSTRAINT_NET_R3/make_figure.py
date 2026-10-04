from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    c = read(HERE / 'LAB_MEASUREMENT_COVERAGE.json')
    rows = sorted([r for r in c['groups'] if r['directly_touched_count']], key=lambda r: (-r['directly_touched_count'], r['id']))[:15]
    labels = [r['id'] + ' ' + r['name'][:49] for r in rows]
    y = list(range(len(rows)))
    (fig, (ax, bx)) = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={'width_ratios': [2.5, 1]})
    ax.barh(y, [r['directly_touched_count'] for r in rows], label='The relevant edges', color='#aac6dd')
    ax.barh(y, [r['sole_listed_acquisition_count'] for r in rows], height=0.42, label='Only type of measurement indicated', color='#29658d')
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Edge-ID (PHENOMENOLOGICAL Planning Dimensions)')
    ax.legend(loc='lower right')
    ax.set_title('Crossing of the measurement requirements with X71 and R3 suggestions')
    bx.bar(['Numeriskt\ngap', 'Physical\nmeasurement deficiency', 'Annat\nhinder'], [17, 51, 5], color=['#347c66', '#c48830', '#9e5361'])
    bx.set_ylabel('OPEN/UNKNOWN edges')
    bx.set_ylim(0, 60)
    for (i, n) in enumerate([17, 51, 5]):
        bx.text(i, n + 0.8, str(n), ha='center')
    bx.set_title('73 open or unknown edges')
    fig.suptitle('What measurements can make the gaps trialsome?', fontsize=17)
    fig.text(0.5, 0.025, '0 physically closed edges. The stacks show measurement type coverage, conditional on the same specimen/state and remaining model requirements.\nNumber not sufficient: identical number of measurements gives 12 and 0 matched triples (identity error 0 ).', ha='center', fontsize=10)
    fig.tight_layout(rect=[0, 0.085, 1, 0.95])
    fig.savefig(HERE / 'CONSTRAINT_NET_R3.png', dpi=150)
    plt.close(fig)
if __name__ == '__main__':
    main()
