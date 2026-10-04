import pathlib, json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from fixtures import fixture
from fractions import Fraction as Q
R = pathlib.Path(__file__).resolve().parents[1]

def run(root=R):
    root = pathlib.Path(root)
    r1 = json.load(open(root / 'rounds/R1/COHORT.json'))
    r3 = json.load(open(root / 'rounds/R3/FIRST_CONTACT.json'))
    r4 = json.load(open(root / 'rounds/R4/AXES.json'))
    fig = plt.figure(figsize=(12, 8), layout='constrained')
    gs = fig.add_gridspec(2, 2)
    a = fig.add_subplot(gs[0, 0])
    b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 0], projection='3d')
    d = fig.add_subplot(gs[1, 1])
    names = [x['key'][:8] for x in r3['rows']]
    heights = [x['first_surface_contact_height_mm'] for x in r3['rows']]
    a.barh(names, heights, color='#bd3b36')
    a.set_xlabel('First source-surface contact above seated pose (mm)')
    a.set_title('5 blocked crowns: exact digital first contact')
    a.invert_yaxis()
    vals = np.array([[0 if x['neighbours'][s]['endpoint']['status'] == 'SURFACE_CLEAR_CERTIFIED' else 2, 1 if x['neighbours'][s]['path']['status'] == 'COLLISION' else 0 if x['neighbours'][s]['path']['status'] == 'SURFACE_CLEAR_CERTIFIED' else 2] for x in r1['rows'] for s in ['mesial', 'distal']])
    from matplotlib.colors import ListedColormap
    b.imshow(vals, aspect='auto', cmap=ListedColormap(['#66a58a', '#bd3b36', '#b4b4b4']), vmin=0, vmax=2)
    b.set_xticks([0, 1], ['Seated surface', 'Whole path surface'])
    b.set_yticks([0, 12, 24, 35], ['pair1', 'pair13', 'pair25', 'pair36'])
    b.set_title('36 original pairs: green separated / red hit / gray missing')
    b.set_ylabel('All 18 sites × 2 neighbour regions')
    w = max((z['first_witness'] for z in r3['rows'][0]['all_pair_certificates'].values() if z['first_witness']), key=lambda x: Q(x['s_exact']))
    A = np.array(w['source_triangle_mm']) + np.array(w['translation_mm']) * float(Q(w['s_exact']))
    B = np.array(w['obstacle_triangle_mm'])
    pt = np.array(w['point_mm'])
    c.add_collection3d(Poly3DCollection([A], alpha=0.6, facecolor='#428cc4', edgecolor='k'))
    c.add_collection3d(Poly3DCollection([B], alpha=0.6, facecolor='#bd3b36', edgecolor='k'))
    c.scatter(*pt, color='black', s=40)
    allp = np.r_[A, B]
    lo = allp.min(0)
    hi = allp.max(0)
    cen = (lo + hi) / 2
    span = max(float(np.ptp(allp, axis=0).max()), 0.1) * 0.6
    c.set_xlim(cen[0] - span, cen[0] + span)
    c.set_ylim(cen[1] - span, cen[1] + span)
    c.set_zlim(cen[2] - span, cen[2] + span)
    c.set_xlabel('local x (mm)')
    c.set_ylabel('local y (mm)')
    c.set_zlabel('local z (mm)')
    c.set_title('Original facet witness at first contact; PER_POINT')
    d.axis('off')
    found = [x for x in r4['rows'] if x['selected_translation_mm']]
    text = 'Whole preparation / physical insertion: UNKNOWN\n\nEndpoint minimum witness: 0.1 mm in both states\nIdentity error: exactly 0\nIntermediate bulge: 0.2 mm -> COLLISION\n\nInside-only apex expansion: 0.2 mm\nProtected exterior/margin identity error: 0 mm\nCollision survives -> no feasible inside-only repair\n\nFinite straight-axis search: ' + str(len(found)) + '/5 surface paths found\nBoth fail virtual local normal condition\nComplete mesh and margin unchanged\nNo scanner or physical error band measured'
    d.text(0, 0.95, text, va='top', fontsize=10, linespacing=1.6)
    (root / 'figures').mkdir(exist_ok=True)
    fig.savefig(root / 'figures/demo.png', dpi=160)
    plt.close(fig)
if __name__ == '__main__':
    run()
