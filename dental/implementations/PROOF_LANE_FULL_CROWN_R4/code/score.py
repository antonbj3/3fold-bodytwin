"""Held-surface scorer. This process runs only after predictions are frozen."""
from geometry import *
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
saved = sys.modules['common']
try:
    sys.modules['common'] = module('v6common_r4', BASE / 'PROOF_LANE_GENCAD_V6/code/common.py')
    contact = module('v6contact_r4', BASE / 'PROOF_LANE_GENCAD_V6/code/contact.py')
finally:
    sys.modules['common'] = saved
height = module('memory_height_r4', BASE / 'PROOF_LANE_FULL_CROWN_R3/code/memory_height.py').height

def proximal(t, ext):
    if not len(t):
        return dict(status='UNKNOWN_NO_NEIGHBOUR')
    lo = ext.min((0, 1)) - 0.2
    hi = ext.max((0, 1)) + 0.2
    tt = t[np.all(t.max(1) >= lo, axis=1) & np.all(t.min(1) <= hi, axis=1)]
    if not len(tt):
        return dict(status='NO_NEAR_CANDIDATE', area_mm2=0.0, count=0)
    (v, inv) = np.unique(tt.reshape(-1, 3), axis=0, return_inverse=True)
    f = inv.reshape(-1, 3)
    d = distance(ext, v)
    aa = []
    mm = []
    idx = []
    links = []
    edges = {}
    for (fi, ids) in enumerate(f):
        if d[ids].min() > 0.1:
            continue
        p = contact.Q.clip(np.eye(3), 0.1 - d[ids]) @ v[ids]
        if len(p) < 3:
            continue
        a = np.linalg.norm(np.cross(p[1:-1] - p[0], p[2:] - p[0]), axis=1) / 2
        ar = float(a.sum())
        if ar <= 1e-12:
            continue
        j = len(idx)
        idx.append(fi)
        aa.append(ar)
        mm.append(np.sum(a[:, None] * (p[0] + p[1:-1] + p[2:]) / 3, axis=0))
        for (x, y) in zip(ids, np.roll(ids, -1)):
            if min(d[x], d[y]) < 0.1 or d[x] == d[y] == 0.1:
                k = tuple(sorted((int(x), int(y))))
                if k in edges:
                    links.append((j, edges[k]))
                else:
                    edges[k] = j
    if idx:
        ij = np.array(links, int).reshape(-1, 2)
        (n, lab) = connected_components(coo_matrix((np.ones(len(ij)), (ij[:, 0], ij[:, 1])), shape=(len(idx), len(idx))).tocsr(), directed=False)
    else:
        n = 0
    return dict(status='SCORED_UNSIGNED_PL_DISTANCE', area_mm2=float(sum(aa)), count=int(n), centroid_mm=np.sum(mm, axis=0) / sum(aa) if aa else None, min_sampled_distance_mm=float(d.min()), signed_interference='UNKNOWN')

def wall(ext, inner):
    t = ext.copy()
    global_lower = np.inf
    min_sample = np.inf
    refined = 0
    unresolved = 0
    for level in range(5):
        q = t.mean(1)
        rad = np.linalg.norm(t - q[:, None], axis=2).max(1)
        d = distance(inner, q)
        lb = d - rad
        min_sample = min(min_sample, float(d.min()))
        refine = (lb < 0.5) & (d >= 0.5) & (rad > 0.04) & (level < 4)
        if np.any(~refine):
            global_lower = min(global_lower, float(lb[~refine].min()))
        if level == 4:
            unresolved = int(np.sum((lb < 0.5) & (d >= 0.5)))
        if not refine.any():
            break
        (a, b, c) = t[refine].transpose(1, 0, 2)
        ab = (a + b) / 2
        bc = (b + c) / 2
        ca = (c + a) / 2
        t = np.concatenate([np.stack([a, ab, ca], 1), np.stack([ab, b, bc], 1), np.stack([ca, bc, c], 1), np.stack([ab, bc, ca], 1)])
        refined += len(t)
    return dict(sampled_min_mm=min_sample, continuous_lower_mm=global_lower, adaptive_facets=refined, unresolved_facets=unresolved, scope='all exterior facets against explicit intaglio; annulus excluded', proof='1-Lipschitz triangle distance: d(centroid,inner)-max_vertex_radius. Exact-arithmetic geometric bound; rigorous floating/source enclosure MISSING')

def evaluate(tag='B', resume=False):
    st = time.perf_counter()
    fr = read(ROOT / f'FROZEN_PREDICTIONS_{tag}.json')
    lock = {r['key']: r for r in read(ROOT / 'FROZEN_INPUTS_B.json')['records']}
    pr = read(ROOT / 'PREREG_B.json')
    rows = []
    cache = {}
    old = read(ROOT / 'raw' / f'{tag}_SCORE.json') if resume and (ROOT / 'raw' / f'{tag}_SCORE.json').exists() else []
    for (index, rec) in enumerate(fr['rows']):
        if index < len(old):
            assert all((old[index][k] == rec[k] for k in ['key', 'method', 'mesh_sha256']))
            rows.append(old[index])
            continue
        r = dict(rec)
        t0 = time.perf_counter()
        key = r['key']
        rr = lock[key]
        if r['status'] == 'GENERATED':
            assert sha(r['mesh_path']) == r['mesh_sha256']
            try:
                if key not in cache:
                    assert sha(rr['private_path']) == rr['private_sha256']
                    p = npz(rr['public_path'])
                    target = npz(rr['private_path'])['target']
                    donor = p['donor'].copy()
                    donor[:, :, 0] *= -1
                    (R, u, errs) = fit(sample(donor, 1536), sample(target, 1536))
                    floor = metrics(donor @ R.T + u, target)
                    (xy, ff) = grid(target)
                    ceil = height(p['antagonist'], xy, True)
                    rg = ceil - height(target, xy)
                    native = {s: proximal(p[s], target) for s in ['mesial', 'distal']}
                    cache[key] = (p, target, floor, xy, ff, ceil, rg, native)
                (p, target, floor, xy, ff, ceil, rg, native) = cache[key]
                m = npz(r['mesh_path'])
                (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
                ext = v[f[roles == 0]]
                inner = v[f[roles == 1]]
                fixed = metrics(ext, target)
                (R, u, errs) = fit(sample(ext, 1536), sample(target, 1536))
                aligned = metrics(ext @ R.T + u, target)
                coarse = metrics(ext, target, 2048)
                mm = trimesh.Trimesh(v, f, process=False)
                outer = trimesh.Trimesh(v, f[roles == 0], process=False)
                outer.remove_unreferenced_vertices()
                ll = loops(outer)
                curve = outer.vertices[ll[0]]
                margin = curve_error(curve, p['margin_curve'])
                w = wall(ext, inner)
                pg = ceil - height(ext, xy)
                cm = contact.compare(xy, ff, pg, rg)
                cg = contact.gates(cm, pr['metrics']['v6'])
                encl = contact.enclosure(xy, ff, pg, rg, 0.05)
                prox = {s: proximal(p[s], ext) for s in ['mesial', 'distal']}
                threshold = pr['metrics']['shape_threshold_mm'][r['family']]
                gg = dict(shape=fixed['p95_mm'] <= threshold, registered_shape=aligned['p95_mm'] <= threshold, closed=bool(mm.is_watertight and mm.is_winding_consistent and (mm.volume > 0) and (len(mm.split(only_watertight=False)) == 1)), margin=margin <= 0.025, wall_sampled=w['sampled_min_mm'] >= 0.5, wall_geometric_bound=w['continuous_lower_mm'] >= 0.5, mesial=0.05 <= prox['mesial'].get('area_mm2', 0) <= 8, distal=0.05 <= prox['distal'].get('area_mm2', 0) <= 8, contact_nominal=all(cg.values()))
                r.update(status='SCORED', resolution='PER_TOOTH', in_situ=fixed, registered_shape=aligned, probe_sensitivity_p95_mm=abs(coarse['p95_mm'] - fixed['p95_mm']), natural_floor=floor, in_situ_excess_over_floor_mm=fixed['p95_mm'] - floor['p95_mm'], registered_excess_over_floor_mm=aligned['p95_mm'] - floor['p95_mm'], in_situ_to_floor_ratio=fixed['p95_mm'] / floor['p95_mm'], threshold_mm=threshold, margin_curve_sampled_max_mm=margin, wall=w, contact=cm, contact_gates=cg, contact_enclosure=encl, proximal=prox, native_proximal=native, gates=gg, technical_joint=all((gg[k] for k in ['shape', 'closed', 'margin', 'wall_geometric_bound'])), positive_native_band=cm.get('reference', {}).get('area_mm2', 0) > 0, full_functional_status='UNKNOWN_NO_VERIFIED_BITE_OR_ACTUAL_PREPARATION', self_intersection='UNKNOWN', source_uncertainty='UNKNOWN')
                dest = DATA / 'fields' / tag / r['method'] / (key + '.npz')
                dest.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(dest, xy=xy, faces=ff, predicted_gap=pg, reference_gap=rg)
                r.update(contact_field_path=dest, contact_field_sha256=sha(dest))
            except Exception as e:
                r.update(status='SCORE_FAILED', reason=repr(e))
        r['score_seconds'] = time.perf_counter() - t0
        rows.append(r)
        dump(ROOT / 'raw' / f'{tag}_SCORE.json', rows)
        print(tag, key, r['method'], r['status'], r.get('in_situ', {}).get('p95_mm'), r.get('gates'), flush=True)
        state(tag + '_SCORING', str(len(rows)) + '/' + str(len(fr['rows'])), 'Complete fixed test-set before deciding next construction')
    summary = {}
    for method in sorted({r['method'] for r in rows}):
        summary[method] = {}
        for family in ['anterior', 'premolar', 'molar']:
            ss = [r for r in rows if r['method'] == method and r['family'] == family and (r['status'] == 'SCORED')]
            summary[method][family] = dict(requested=6, scored=len(ss), shape_pass=sum((r['gates']['shape'] for r in ss)), registered_shape_pass=sum((r['gates']['registered_shape'] for r in ss)), technical_joint=sum((r['technical_joint'] for r in ss)), median_p95_mm=float(np.median([r['in_situ']['p95_mm'] for r in ss])) if ss else None, median_registered_p95_mm=float(np.median([r['registered_shape']['p95_mm'] for r in ss])) if ss else None, median_natural_floor_mm=float(np.median([r['natural_floor']['p95_mm'] for r in ss])) if ss else None, median_excess_mm=float(np.median([r['in_situ_excess_over_floor_mm'] for r in ss])) if ss else None, gate_counts={k: sum((r['gates'][k] for r in ss)) for k in ss[0]['gates']} if ss else {})
    result = dict(claim_type='capability', round=tag, summary=summary, rows=rows, external_referent=REFERENT, dropout=dict(requested=len(rows), rejected=sum((r['status'] != 'SCORED' for r in rows)), fraction=sum((r['status'] != 'SCORED' for r in rows)) / len(rows)), seconds=time.perf_counter() - st + sum((r['score_seconds'] for r in old)), resumed_prefix_rows=len(old), peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'rounds' / f'{tag}.json', result)
    state(tag + '_DECIDED', summary, 'Preserve outcome; change binding preparation/observation representation next')
    (ROOT / f'HANDOFF_{tag}.md').write_text(json.dumps(clean(summary), indent=2) + '\nNo verified bite or measured preparation. All failed gates retained.\n')
    return result
if __name__ == '__main__':
    evaluate()
