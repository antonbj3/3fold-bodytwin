from cadlib import *
from margin import *
import networkx as nx
from scipy.optimize import minimize_scalar
import resource

def detect_cut(m):
    fc = m.triangles_center[:, 2]
    (qlo, qhi) = np.quantile(fc, [0.2, 0.8])
    G = nx.DiGraph()
    s = len(fc)
    t = s + 1
    ea = m.face_adjacency_edges
    length = np.linalg.norm(m.vertices[ea[:, 0]] - m.vertices[ea[:, 1]], axis=1)
    cap = np.maximum(1, np.rint(1000000.0 * length / (0.05 + m.face_adjacency_angles) ** 2)).astype(np.int64)
    for ((i, j), w) in zip(m.face_adjacency, cap):
        G.add_edge(int(i), int(j), capacity=int(w))
        G.add_edge(int(j), int(i), capacity=int(w))
    inf = int(cap.sum()) + 1
    for i in range(len(fc)):
        if fc[i] >= qhi:
            G.add_edge(s, i, capacity=inf)
        if fc[i] <= qlo:
            G.add_edge(i, t, capacity=inf)
    (val, (A, B)) = nx.minimum_cut(G, s, t, flow_func=nx.algorithms.flow.preflow_push)
    cut = np.array([(int(i) in A) != (int(j) in A) for (i, j) in m.face_adjacency])
    actual = int(cap[cut].sum())
    if val != actual:
        raise ValueError('Mincut value differs from crossing capacities')
    (cs, rej) = cycles(ea[cut], m.vertices)
    info = dict(cut_capacity=actual, dual_minimum_cut_value=int(val), crossing_edges=int(cut.sum()), closed_components=len(cs), rejected_components=rej, automatic_terminal_z_mm=[qlo, qhi])
    if len(cs) != 1 or rej:
        return dict(status='UNKNOWN_MULTIPLE_OR_OPEN_CUT', diagnostics=info)
    return dict(status='PROPOSED', **cs[0], diagnostics=info)

def run():
    st = time.perf_counter()
    rows = []
    for r in read(ROOT / 'raw/R1_PREDICTIONS.json'):
        if r['status'] == 'UNAVAILABLE_SOURCE':
            rows.append(dict(key=r['key'], cohort=r['cohort'], status=r['status']))
            continue
        p = Path(r['source'])
        if sha(p) != r['source_sha256']:
            raise ValueError('source drift')
        if p.suffix == '.npz':
            a = np.load(p)
            m = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
        else:
            m = trimesh.load_mesh(p, process=True)
        t0 = time.perf_counter()
        d = detect_cut(m)
        rows.append(dict(key=r['key'], cohort=r['cohort'], fault=r.get('fault'), **d, seconds=time.perf_counter() - t0))
        dump(ROOT / 'raw/R2_PARTIAL.json', rows)
        print(r['key'], d['status'], round(time.perf_counter() - t0, 3), flush=True)
    dump(ROOT / 'raw/R2_PREDICTIONS.json', rows)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R2.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R2.json'), prediction_sha256=sha(ROOT / 'raw/R2_PREDICTIONS.json'), code_sha256=sha(Path(__file__)), seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    score()

def score():
    rows = read(ROOT / 'raw/R2_PREDICTIONS.json')
    refs = np.load(DATA / 'R1_references.npz')
    for r in rows:
        if r['status'] == 'PROPOSED' and r['key'] in refs:
            r['margin_error'] = errors(np.array(r['points']), refs[r['key']])
            r['margin_gate'] = 'PASS' if r['margin_error']['max_um'] <= 25 else 'FAIL'
        else:
            r['margin_gate'] = 'UNKNOWN'
    summary = {c: dict(requested=sum((r['cohort'] == c for r in rows)), proposed=sum((r['cohort'] == c and r['status'] == 'PROPOSED' for r in rows)), pass_count=sum((r['cohort'] == c and r['margin_gate'] == 'PASS' for r in rows)), fail_count=sum((r['cohort'] == c and r['margin_gate'] == 'FAIL' for r in rows)), unknown_count=sum((r['cohort'] == c and r['margin_gate'] == 'UNKNOWN' for r in rows))) for c in ['V2', 'X52']}
    dump(ROOT / 'raw/R2_SCORED.json', rows)
    dump(ROOT / 'RESULTS_R2.json', dict(claim_type='capability', summary=summary, external_referent=read(ROOT / 'PREREG_R2.json')['external_referent']))
    (ROOT / 'HANDOFF_R2.md').write_text('R2 global dual-graph separating cut: ' + json.dumps(summary) + '. Closed does not mean correct. Preserve 25um gate. Next: bind marginal ambiguity to real crown generation, refuse physical qualification and test complete geometry against independently measured IOS morphology.\n')
    state('R2_COMPLETE', summary, 'Whole crown generation and context constraints with explicit source semantics')
    print(summary)
if __name__ == '__main__':
    if not (ROOT / 'FROZEN_PREDICTIONS_R2.json').exists():
        run()
    else:
        score()
