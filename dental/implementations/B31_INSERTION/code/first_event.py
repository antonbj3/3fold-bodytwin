import json, pathlib, datetime, time, itertools, hashlib, warnings
import numpy as np
from sweep import *
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def freeze():
    rows = json.load(open(R / 'rounds/R1/COHORT.json'))['rows']
    keys = [r['key'] for r in rows if r['full_path_status'] == 'COLLISION']
    p = R / 'PREREG_R3.json'
    p.write_text(json.dumps(dict(round='R3', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), claim_type='capability', keys=keys, capability='Certify first contact height along fixed straight insertion for each of five blocked sites', obstacle='One intersection witness refutes path but does not locate first encounter', changed_operation='Exhaust all full-sweep triangle-pair candidates; maximize s with exact rational primal+dual certificate and rigorously upper-bound every remaining pair', metric='First-contact translation parameter exact rational, outward safe interval(s_contact,1]; exact reference boolean agreement; unresolved earlier candidates -> UNKNOWN', tolerance='no tolerance on certificates; margin0.01mm retained unchanged', strongest_equally_informed_control='Independent5variable barycentric HiGHS maximization, plus input-bound rational check', falsifier='Change dual coefficient by+1or claimed first height by1mm: exact upper-bound certificate must reject', full_cost={'fit': 0, 'prep': 'R1candidate AABBs', 'discovery': 'all pair optimization', 'validation': 'primal+dual reconstruction and independent LP', 'query': 'all5blocked sites', 'fallback': 'earliest UNKNOWN if uncovered', 'lab': 'UNKNOWN'}, resolution='PER_POINT', timescale='SIMULTANEOUS'), indent=2) + '\n')
    p.with_suffix('.sha256').write_text(sha(p) + '\n')
    print('R3FROZEN first-event5sites')

def dual_valid(A, B, D, w, y):
    (M, b) = primal_matrix(A, B, D)
    x = list(map(Q, w['weights']))
    y = list(map(Q, y))
    c = [Q(0)] * 8
    c[6] = Q(-1)
    return len(y) == 6 and all((sum((M[i][k] * y[i] for i in range(6))) <= c[k] for k in range(8))) and (sum((bi * yi for (bi, yi) in zip(b, y))) == sum((ci * xi for (ci, xi) in zip(c, x))))

def dual(A, B, D, w):
    (M, b) = primal_matrix(A, B, D)
    c = [Q(0)] * 8
    c[6] = Q(-1)
    x = list(map(Q, w['weights']))
    pos = [i for (i, v) in enumerate(x) if v > 0]
    extra = [i for i in range(8) if i not in pos]
    for count in range(len(pos), min(6, 8) + 1):
        for ex in itertools.combinations(extra, count - len(pos)):
            ids = pos + list(ex)
            y = rref_solve([[M[j][i] for j in range(6)] for i in ids], [c[i] for i in ids])
            if y is not None and dual_valid(A, B, D, w, y):
                return list(map(str, y))
    return None

def first(tri, roles, B, D):
    obs = Obstacle(B)
    (lo, hi) = swept_bounds(tri, D)
    pairs = []
    for (i, (a, b)) in enumerate(zip(lo, hi)):
        pairs.extend(((i, j) for j in obs.tree.intersection(tuple(np.r_[a, b])))) if obs.tree else None
    certs = []
    unresolved = []
    separated = 0
    lpqueries = 0
    for k in range(0, len(pairs), 4096):
        block = pairs[k:k + 4096]
        ii = np.array([q[0] for q in block])
        jj = np.array([q[1] for q in block])
        (sep, _) = separator_batch(tri[ii], B[jj], D)
        separated += int(sep.sum())
        for t in np.flatnonzero(~sep):
            (i, j) = block[t]
            lpqueries += 1
            (w, stat) = lp_witness(tri[i], B[j], D, True)
            upper = (rational(B[j, :, 2].max()) - rational(tri[i, :, 2].min())) / rational(D[2])
            if w:
                y = dual(tri[i], B[j], D, w)
                if y:
                    w.update(dual=y, source_face_id=int(i), obstacle_face_id=int(j), source_role=int(roles[i]), translation_mm=D.tolist(), source_triangle_mm=tri[i].tolist(), obstacle_triangle_mm=B[j].tolist())
                    certs.append(w)
                else:
                    unresolved.append(dict(source_face_id=int(i), obstacle_face_id=int(j), aabb_s_upper=str(upper), reason='exact dual reconstruction failed'))
            else:
                unresolved.append(dict(source_face_id=int(i), obstacle_face_id=int(j), aabb_s_upper=str(upper), reason='LP/intersction witness unavailable', LP_status=stat))
    best = max(certs, key=lambda w: Q(w['s_exact'])) if certs else None
    upper = Q(best['s_exact']) if best else Q(0)
    danger = [q for q in unresolved if Q(q['aabb_s_upper']) > upper]
    return dict(status='FIRST_SURFACE_CONTACT_CERTIFIED' if best and (not danger) else 'UNKNOWN', first_witness=best, candidate_pairs=len(pairs), separated_pairs=separated, LP_queries=lpqueries, intersection_certificates=certs, unresolved=unresolved, earlier_unresolved=len(danger), absent_surface=not len(B))

def run(outdir):
    pre = json.load(open(R / 'PREREG_R3.json'))
    co = json.load(open(R / 'FROZEN_COHORT.json'))
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    tick = time.perf_counter()
    for key in pre['keys']:
        rec = next((q for q in co['rows'] if q['key'] == key))
        with np.load(rec['mesh_path']) as m, np.load(rec['public_path']) as p:
            tri = m['vertices'][m['faces']]
            D = np.array([0.0, 0.0, 40.0])
            sides = {s: first(tri, m['roles'], p[s], D) for s in ['mesial', 'distal']}
            witnesses = [q['first_witness'] for q in sides.values() if q['first_witness']]
            best = max(witnesses, key=lambda w: Q(w['s_exact']))
            bad = [q for q in sides.values() if q['absent_surface'] or any((Q(x['aabb_s_upper']) > Q(best['s_exact']) for x in q['unresolved']))]
            A = np.array(best['source_triangle_mm'])
            B = np.array(best['obstacle_triangle_mm'])
            assert validate_witness(A, B, D, best) and dual_valid(A, B, D, best, best['dual'])
            yy = best['dual'].copy()
            yy[0] = str(Q(yy[0]) + 1)
            reject = not dual_valid(A, B, D, best, yy)
            assert reject
            independent = independent_lp(A, B, D)
            row = dict(key=key, status='FIRST_SURFACE_CONTACT_CERTIFIED' if not bad else 'UNKNOWN', first_surface_contact_height_mm=float(Q(best['s_exact']) * 40), height_exact_mm=str(Q(best['s_exact']) * 40), source_face_id=best['source_face_id'], point_mm=best['point_mm'], source_role=best['source_role'], safe_surface_interval_s=['(' + best['s_exact'] + ',1]'], all_pair_certificates=sides, dual_corruption_rejected=reject, independent_control=independent, resolution='PER_POINT', physical_first_contact='UNKNOWN: prep, registration, scene solids missing')
            rows.append(row)
            (outdir / 'FIRST_CONTACT.json').write_text(json.dumps(dict(rows=rows, seconds=time.perf_counter() - tick), indent=2) + '\n')
            print(key, row['status'], row['first_surface_contact_height_mm'], [(q['candidate_pairs'], len(q['intersection_certificates']), q['earlier_unresolved']) for q in sides.values()], flush=True)
    return rows
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--freeze', action='store_true')
    ap.add_argument('--out', default=str(R / 'rounds/R3'))
    a = ap.parse_args()
    freeze() if a.freeze else run(pathlib.Path(a.out))
