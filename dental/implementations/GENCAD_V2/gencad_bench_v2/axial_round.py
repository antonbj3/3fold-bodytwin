"""R2 full crown caps: new representation, immutable R1 roof predictions."""
import sys, time, math, collections
from pathlib import Path
import numpy as np
from .common import *
from .geometry import boundary_edges, crop, export_stl, parse_stl
from .tasks import load_task
sys.path.insert(0, str(ROOT / 'gencad_bench_v2/vendor/v1'))
sys.path.insert(0, str(ROOT / 'gencad_bench_v2/vendor'))
from gencad_bench.checks.exact import vec, sub, cross, dot, norm2, q, aabb_distance2, triangle_distance, strings
from gencad_bench.checks import cement, cement_max
from continuous import broad_phase, extremum

def cone(rb, rt, z0, z1, n=24):
    xy = np.round(np.c_[np.cos(np.arange(n) * 2 * np.pi / n), np.sin(np.arange(n) * 2 * np.pi / n)], 6)
    lo = np.c_[np.round(xy * round(rb, 6), 12), np.full(n, round(z0, 9))]
    hi = np.c_[np.round(xy * round(rt, 6), 12), np.full(n, round(z1, 9))]
    v = np.r_[lo, hi, [[0.0, 0.0, round(z1, 9)]]]
    f = []
    for i in range(n):
        j = (i + 1) % n
        f.extend([[i, j, n + j], [i, n + j, n + i], [n + i, n + j, 2 * n]])
    return (v, np.asarray(f, int))

def cap_planes(v, f):
    P = list(map(vec, v))
    planes = []
    seen = set()
    for tri in f:
        (a, b, c) = [P[k] for k in tri]
        n = cross(sub(b, a), sub(c, a))
        bb = dot(n, a)
        if not norm2(n):
            raise ValueError('degenerate conical cap')
        s = max((abs(x) for x in n))
        n = tuple((x / s for x in n))
        bb /= s
        if n[2] < 0 or any((dot(n, p) > bb for p in P)):
            raise ValueError('cap is not an exact convex extraction domain')
        row = tuple(strings(n + (bb,)))
        if row not in seen:
            seen.add(row)
            planes.append(row)
    return planes

def complete(t, design, flare):
    xy = np.asarray(t['xy'])
    faces = np.asarray(t['faces'])
    roof = np.asarray(design['outer_vertices'], float)
    roof[:, 2] = np.ceil(roof[:, 2] * 1000000000.0) / 1000000000.0
    inner_z = np.floor(float(np.asarray(design['inner_vertices'])[:, 2].min()) * 1000000000.0) / 1000000000.0
    z0 = inner_z - 3.0
    edges = boundary_edges(faces)
    ids = np.unique(edges)
    a = xy[edges[:, 0]]
    b = xy[edges[:, 1]]
    rin = float(np.min(np.abs(np.cross(a, b)) / np.linalg.norm(b - a, axis=1)))
    rb = rin - 0.35
    rt = rb - 3 * np.tan(np.radians(3))
    if rt < 0.25:
        raise ValueError('no positive conical cavity')
    (inn, inf) = cone(rb, rt, z0, inner_z)
    base = xy[ids] * (1 + flare / np.linalg.norm(xy[ids], axis=1)[:, None])
    ov = np.r_[roof, np.c_[np.round(base, 12), np.full(len(ids), round(z0, 9))]]
    lookup = {int(k): i + len(roof) for (i, k) in enumerate(ids)}
    of = faces.tolist()
    for (a, b) in edges:
        of.extend([[int(a), lookup[int(b)], int(b)], [int(a), lookup[int(a)], lookup[int(b)]]])
    of = np.asarray(of, int)
    ob = np.array([lookup[int(k)] for k in ids])
    oangle = np.mod(np.arctan2(ov[ob, 1], ov[ob, 0]), 2 * np.pi)
    ob = ob[np.argsort(oangle)]
    oa = np.sort(oangle)
    ib = np.arange(24)
    ia = np.mod(np.arctan2(inn[:24, 1], inn[:24, 0]), 2 * np.pi)
    ib = ib[np.argsort(ia)]
    ia = np.sort(ia)
    i = j = 0
    rim = []
    ni = len(ob)
    nj = len(ib)
    offset = len(ov)
    while i < ni or j < nj:
        no = oa[(i + 1) % ni] + (2 * np.pi if i + 1 >= ni else 0) if i < ni else np.inf
        nn = ia[(j + 1) % nj] + (2 * np.pi if j + 1 >= nj else 0) if j < nj else np.inf
        if no <= nn:
            rim.append([ob[i % ni], offset + ib[j % nj], ob[(i + 1) % ni]])
            i += 1
        else:
            rim.append([ob[i % ni], offset + ib[j % nj], offset + ib[(j + 1) % nj]])
            j += 1
    v = np.r_[ov, inn]
    f = np.r_[of, inf[:, ::-1] + offset, np.asarray(rim, int)]
    gap = 0.08
    prep0 = rb - gap / np.cos(np.radians(3))
    prep1 = rt + gap * np.tan(np.radians(3)) - gap / np.cos(np.radians(3))
    (pv, pf) = cone(prep0, prep1, z0, inner_z - gap)
    return dict(outer_vertices=ov, outer_faces=of, inner_vertices=inn, inner_faces=inf, vertices=v, faces=f, preparation_vertices=pv, preparation_faces=pf, flare_mm=flare, cavity_mouth_radius_mm=rb, cavity_top_radius_mm=rt, margin_z_mm=z0, serialization='homothetic decimal cone; outer roof ceil and inner roof floor at 1e-9 mm', scope='full digital caps with planar polygonal finish line; rounded chamfer not validated')

def wall_check(ov, of, iv, inf, w, max_pairs=2500):
    ov = np.asarray(ov, float)
    iv = np.asarray(iv, float)
    of = np.asarray(of)
    inf = np.asarray(inf)
    for (v, f) in [(ov, of), (iv, inf)]:
        if v.ndim != 2 or v.shape[1] != 3 or f.ndim != 2 or (f.shape[1] != 3) or (f.dtype.kind not in 'iu') or (not len(v)) or (not len(f)) or (not np.isfinite(v).all()) or (f.min() < 0) or (f.max() >= len(v)):
            return dict(status='INVALID')
    A = ov[of]
    B = iv[inf]
    lo = A.min(1)
    hi = A.max(1)
    blo = B.min(1)
    bhi = B.max(1)
    g = np.maximum(0, np.maximum(blo[None] - hi[:, None], lo[:, None] - bhi[None]))
    order = np.argsort(np.sum(g * g, axis=2).ravel(), kind='stable')
    Ar = [tuple(map(vec, t)) for t in A]
    Br = [tuple(map(vec, t)) for t in B]
    w = q(w)
    count = 0
    pruned = 0
    start = time.perf_counter()
    for index in order:
        (i, j) = divmod(int(index), len(B))
        (a, b) = (Ar[i], Br[j])
        if aabb_distance2(a, b) >= w * w:
            pruned += 1
            continue
        if count >= max_pairs:
            return dict(status='UNKNOWN', reason='exact pair budget reached', exact_pairs=count, pruned=pruned, seconds=time.perf_counter() - start)
        (dd, p1, p2) = triangle_distance(a, b)
        count += 1
        if dd < w * w:
            return dict(status='FAIL', witness=strings(dict(outer_triangle=i, inner_triangle=j, distance_squared_mm2=dd, outer_point=p1, inner_point=p2)), exact_pairs=count, seconds=time.perf_counter() - start)
    return dict(status='PASS', exact_pairs=count, pruned=pruned, seconds=time.perf_counter() - start, scope='all rational cap triangle pairs; cervical annular rim excluded')

def controls():
    v = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    f = np.array([[0, 1, 2]])
    good = wall_check(v + [0, 0, 2], f, v, f, 1)
    bad = wall_check(v + [0, 0, 0.5], f, v, f, 1)
    invalid = wall_check([[0, 0]], f, v, f, 1)
    if good['status'] != 'PASS' or bad['status'] != 'FAIL' or invalid['status'] != 'INVALID':
        raise ValueError('R2 wall controls failed')
    result = dict(parallel_triangle_distance_2mm=good, injected_half_mm=bad, injected_dimension=invalid, external_referent=dict(kind='closed_form', locator='Euclidean distance between congruent parallel triangles: distance equals perpendicular plane separation', compared_quantity='minimum separation mm', refutes_us=True))
    dump(ROOT / 'raw/CONTROL_R2.json', result)
    return result

def predict():
    controls()
    start = time.perf_counter()
    rows = []
    dest = DATA / 'axial_R2'
    dest.mkdir(exist_ok=True)
    selected = [p for p in sorted((DATA / 'public').glob('*.json')) if any((p.name.endswith('_' + f + '_normal.json') for f in ['molar_crown', 'premolar_crown']))]
    for p in selected:
        t = load_task(p)
        tid = t['task_id']
        d = read(DATA / 'predictions/constraint_optimizer' / (tid + '.json'))
        if t['status'] != 'READY' or d.get('status') != 'DESIGN':
            rows.append(dict(task_id=tid, status='UNKNOWN_NO_R1_FEASIBLE_ROOF'))
            continue
        for (arm, flare) in [('unflared', 0.0), ('axial_flare', max(0, t['requirements']['wall_mm'] - 0.35) + 0.15)]:
            try:
                design = complete(t, d, flare)
                path = dest / (tid + '_' + arm + '.npz')
                np.savez_compressed(path, **{k: v for (k, v) in design.items() if isinstance(v, np.ndarray)})
                export_stl(dest / (tid + '_' + arm + '.stl'), design['vertices'], design['faces'])
                rows.append(dict(task_id=tid, status='DESIGNED', arm=arm, file=str(path), sha256=sha(path), flare_mm=flare, scope=design['scope']))
            except ValueError as e:
                rows.append(dict(task_id=tid, arm=arm, status='UNKNOWN_CONSTRUCTION', reason=str(e)))
    dump(ROOT / 'raw/PREDICTIONS_R2.json', rows)
    files = {r['file']: r['sha256'] for r in rows if r['status'] == 'DESIGNED'}
    freeze(ROOT / 'FROZEN_PREDICTIONS_R2.json', dict(claim_type='capability', prereg_sha256=sha(ROOT / 'PREREG_R2.json'), rows_sha256=sha(ROOT / 'raw/PREDICTIONS_R2.json'), files=files, code_sha256=sha(Path(__file__)), whole_cap_wall_scored=False, seconds=time.perf_counter() - start))
    return rows

def score():
    import trimesh
    from .build import pair
    controls()
    freeze = read(ROOT / 'FROZEN_PREDICTIONS_R2.json')
    p = freeze['payload']
    if digest(p) != freeze['payload_sha256'] or sha(ROOT / 'raw/PREDICTIONS_R2.json') != p['rows_sha256'] or sha(Path(__file__)) != p['code_sha256']:
        raise ValueError('R2 freeze drift')
    cohort = {r['case_key']: r for r in read(ROOT / 'data/COHORT_LOCK.json')['payload']['cases']}
    references = {r['task_id']: r['reference'] for r in read(DATA / 'private/references.json')}
    rows = []
    pair_cache = None
    last_case = None
    start = time.perf_counter()
    for row in read(ROOT / 'raw/PREDICTIONS_R2.json'):
        if row['status'] != 'DESIGNED':
            rows.append(row)
            continue
        if sha(row['file']) != p['files'][row['file']]:
            raise ValueError('R2 mesh changed')
        t = load_task(DATA / 'public' / (row['task_id'] + '.json'))
        a = np.load(row['file'], allow_pickle=False)
        out = dict(**row)
        m = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
        out['mesh'] = dict(watertight=bool(m.is_watertight), consistent_winding=bool(m.is_winding_consistent), positive_volume=bool(m.volume > 0), volume_mm3=float(m.volume))
        if not m.is_watertight or not m.is_winding_consistent or m.volume <= 0:
            out['status'] = 'INVALID_MESH'
            rows.append(out)
            continue
        out['wall'] = wall_check(a['outer_vertices'], a['outer_faces'], a['inner_vertices'], a['inner_faces'], t['requirements']['wall_mm'])
        cp = cap_planes(a['inner_vertices'], a['inner_faces'])
        pp = cap_planes(a['preparation_vertices'], a['preparation_faces'])
        out['min_film'] = cement.check(a['preparation_vertices'].tolist(), cp, 0.04)
        maxpass = cement_max.verify_positive(pp, a['inner_vertices'].tolist(), a['preparation_vertices'].tolist(), 0.12)
        out['max_film'] = dict(status='PASS' if maxpass else 'UNKNOWN', scope='corresponding convex-preparation witness, nominal geometric distance')
        out['insertion'] = dict(status='PASS', scope='exact convex conical cap, all outward normals have nonnegative z; straight extraction; no remote lower-arch obstacles')
        om = trimesh.Trimesh(a['outer_vertices'], a['outer_faces'], process=False)
        samples = np.r_[a['inner_vertices'], a['inner_vertices'][a['inner_faces']].mean(1)]
        (_, dd, _) = trimesh.proximity.closest_point(om, samples)
        out['float_sample_control'] = dict(min_sample_distance_mm=float(dd.min()), scope='upper bound on true minimum; vertex/face-centroid samples cannot establish wall PASS')
        ref = np.load(DATA / 'private' / references[row['task_id']]['file'])
        (base, R) = (ref['base'], ref['R'])
        case = t['case_key']
        if case != last_case:
            (pair_cache, _) = pair(cohort[case])
            last_case = case
        upper = ((pair_cache['upper']['v'] - base) @ R)[pair_cache['upper']['f']]
        O = a['outer_vertices'][a['outer_faces']]
        U = crop(upper, a['outer_vertices'][:, :2])
        det = np.cross(O[:, 1, :2] - O[:, 0, :2], O[:, 2, :2] - O[:, 0, :2])
        O = O[abs(det) > 1e-09]
        if len(U) and len(O):
            pairs = broad_phase(U, O)
            occ = extremum(U, O, pairs, True) if len(pairs) else None
        else:
            occ = None
        out['occlusion'] = dict(status='PASS' if occ and occ['minimum_gap_mm'] >= -1e-06 else 'FAIL' if occ else 'UNKNOWN', measurement=occ, scope='projected full outer cap against native full upper arch; not force')
        stl = Path(row['file']).with_suffix('.stl')
        m2 = trimesh.load_mesh(stl, process=True)
        out['export'] = dict(sha256=sha(stl), watertight=bool(m2.is_watertight), volume_relative_error=float(abs(m2.volume - m.volume) / m.volume))
        rows.append(out)
        dump(ROOT / 'raw/R2_PARTIAL.json', rows)
        print('R2', row['task_id'], row['arm'], out['wall']['status'], flush=True)
    gains = []
    regressions = []
    for tid in sorted({r['task_id'] for r in rows}):
        arms = {r.get('arm'): r for r in rows if r['task_id'] == tid and r.get('wall')}
        if len(arms) != 2:
            continue
        (b, c) = (arms['unflared'], arms['axial_flare'])
        safe = all((c[k]['status'] == 'PASS' for k in ['min_film', 'max_film', 'occlusion']))
        if b['wall']['status'] != 'PASS' and c['wall']['status'] == 'PASS' and safe:
            gains.append(tid)
        if any((b[k]['status'] == 'PASS' and c[k]['status'] == 'FAIL' for k in ['min_film', 'max_film', 'occlusion'])):
            regressions.append(tid)
    out = dict(claim_type='capability', round='R2', requested_tasks=32, designed_meshes=sum((r.get('wall') is not None for r in rows)), new_wall_pass_without_film_or_occlusal_failure=gains, regressions=regressions, gate='PASS' if gains and (not regressions) else 'FAIL', wall_counts={arm: dict(collections.Counter((r['wall']['status'] for r in rows if r.get('arm') == arm and 'wall' in r))) for arm in ['unflared', 'axial_flare']}, seconds=time.perf_counter() - start, external_referent=read(ROOT / 'PREREG_R2.json')['external_referent'], limits=['Retrospective same-case construction; no new patient generalization', 'Conical rule-based cavity, unvalidated clinical margin and retention', 'Rim excluded from minimum cap separation', 'No physical milling, cement, fracture or tissue validation'])
    dump(ROOT / 'raw/SCORED_ROWS_R2.json', rows)
    dump(ROOT / 'rounds/R2.json', out)
    return out
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('phase', choices=['predict', 'score', 'all'])
    a = ap.parse_args()
    if a.phase in ['predict', 'all'] and (not (ROOT / 'FROZEN_PREDICTIONS_R2.json').exists()):
        predict()
    if a.phase in ['score', 'all']:
        print(score())
