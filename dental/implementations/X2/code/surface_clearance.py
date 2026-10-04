"""Unsigned 3D surface clearance at fixed source-face sample points.
Distance derivation: Eberly, DistancePoint3Triangle3.pdf. Centroid-radius
screening is a triangle containment-ball lower bound, not a voxel SDF.
"""
import json, time, zipfile
import numpy as np
from scipy.spatial import cKDTree
from occlusion_operator import H, DATA, ZIP, read_stl, dump, sha, state, predict_regions

def closest(P, T):
    a = T[:, 0]
    e = T[:, 1] - a
    f = T[:, 2] - a
    w = P - a
    ee = np.einsum('ij,ij->i', e, e)
    ff = np.einsum('ij,ij->i', f, f)
    ef = np.einsum('ij,ij->i', e, f)
    we = np.einsum('ij,ij->i', w, e)
    wf = np.einsum('ij,ij->i', w, f)
    den = ee * ff - ef * ef
    ok = den > 1e-24
    u = np.divide(we * ff - wf * ef, den, out=np.zeros(len(P)), where=ok)
    v = np.divide(wf * ee - we * ef, den, out=np.zeros(len(P)), where=ok)
    Q = a + u[:, None] * e + v[:, None] * f
    d = np.linalg.norm(P - Q, axis=1)
    valid = ok & (u >= 0) & (v >= 0) & (u + v <= 1)
    best = np.where(valid, d, np.inf)
    bestQ = Q.copy()
    for (i, j) in [(0, 1), (1, 2), (2, 0)]:
        aa = T[:, i]
        bb = T[:, j]
        ed = bb - aa
        ed2 = np.einsum('ij,ij->i', ed, ed)
        t = np.clip(np.divide(np.einsum('ij,ij->i', P - aa, ed), ed2, out=np.zeros(len(P)), where=ed2 > 0), 0, 1)
        q = aa + t[:, None] * ed
        dd = np.linalg.norm(P - q, axis=1)
        take = dd < best
        best[take] = dd[take]
        bestQ[take] = q[take]
    return (best, bestQ)

def screened(P, T, hmax=0.2):
    centers = T.mean(1)
    radii = np.linalg.norm(T - centers[:, None, :], axis=2).max(1)
    caps = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, max(25.6, float(radii.max())) + 1e-10]
    best = np.full(len(P), np.inf)
    owner = np.full(len(P), -1, int)
    qp = np.full_like(P, np.nan)
    tested = 0
    lower = -1
    for cap in caps:
        ids = np.flatnonzero((radii > lower) & (radii <= cap))
        lower = cap
        if not len(ids):
            continue
        tree = cKDTree(centers[ids])
        for start in range(0, len(P), 128):
            pts = P[start:start + 128]
            neighbors = tree.query_ball_point(pts, hmax + cap + 1e-12, workers=1)
            count = np.array([len(x) for x in neighbors])
            if count.sum() == 0:
                continue
            pi = np.repeat(np.arange(len(pts)), count)
            ti = ids[np.concatenate([np.asarray(x, int) for x in neighbors if len(x)])]
            bound = np.linalg.norm(pts[pi] - centers[ti], axis=1) - radii[ti]
            keep = bound <= hmax + 1e-12
            pi = pi[keep]
            ti = ti[keep]
            if not len(pi):
                continue
            (dd, q) = closest(pts[pi], T[ti])
            tested += len(dd)
            order = np.lexsort((ti, dd, pi))
            opi = pi[order]
            first = np.r_[True, np.diff(opi) != 0]
            select = order[first]
            gi = start + pi[select]
            values = dd[select]
            take = values < best[gi]
            best[gi[take]] = values[take]
            owner[gi[take]] = ti[select[take]]
            qp[gi[take]] = q[select[take]]
    return (best, owner, qp, tested)

def brute(P, T):
    out = []
    for p in P:
        best = np.inf
        for start in range(0, len(T), 8192):
            tri = T[start:start + 8192]
            (d, _) = closest(np.broadcast_to(p, (len(tri), 3)), tri)
            best = min(best, float(d.min()))
        out.append(best)
    return np.asarray(out)

def main(max_cases=12):
    z = zipfile.ZipFile(ZIP)
    models = json.loads((H / 'raw/segmentation_model.json').read_text())
    rows = []
    for case in range(1, max_cases + 1):
        out = H / 'raw/clearance' / f'{case:03d}.json'
        if out.exists():
            rows.append(json.loads(out.read_text()))
            continue
        state('R4_3D_CLEARANCE_RUNNING', 'R2 fine-grid force gate failed; source gaps are only projected', f'3D clearance and full-triangle control case{case}')
        r = json.loads((H / 'raw/cases' / f'{case:03d}.json').read_text())
        st = time.perf_counter()
        (L, lh) = read_stl(z, r['upper_member'].replace('upper.stl', 'lower.stl'))
        (lc, _, _) = predict_regions(L.mean(1), models['lower'], True)
        with np.load(r['map_path']) as arr:
            gap = arr['gap']
            P = np.c_[arr['xy'], arr['upper_z']]
            ur = arr['upper_region']
        prep = time.perf_counter() - st
        st = time.perf_counter()
        (d, owner, q, tested) = screened(P, L)
        qs = time.perf_counter() - st
        ix = np.linspace(0, len(P) - 1, 8).astype(int)
        st = time.perf_counter()
        oracle = brute(P[ix], L)
        bs = time.perf_counter() - st
        parity = float(np.max(np.abs(np.minimum(d[ix], 0.2) - np.minimum(oracle, 0.2))))
        metrics = []
        for h in [0.05, 0.1, 0.2]:
            a = gap <= h
            b = d <= h
            upper = np.bincount(ur[b], minlength=8)
            lower = np.bincount(lc[owner[b]], minlength=8)
            shares = (np.r_[upper, lower] / b.sum() * 100).tolist() if b.any() else None
            metrics.append(dict(band_mm=h, projected_contact_pixels=int(a.sum()), unsigned_3D_near_pixels=int(b.sum()), projected_only_pixels=int((a & ~b).sum()), threeD_only_pixels=int((b & ~a).sum()), threeD_area_shares_pp=shares, negative_projected_gap_away_from_3D_band=int(((gap < 0) & ~b).sum())))
        path = DATA / f'{case:03d}_clearance.npz'
        np.savez_compressed(path, distance_mm=d, nearest_lower_face=owner, lower_point_mm=q, upper_point_mm=P, upper_region=ur)
        row = dict(case=case, samples=len(P), preparation_s=prep, screened_query_s=qs, brute_control_s=bs, tested_pairs=tested, brute_pairs=len(L) * 8, parity_mm=parity, parity_gate=parity <= 1e-09, oracle_indices=ix.tolist(), oracle_distance_mm=oracle.tolist(), candidate_clipped_distance_mm=np.minimum(d[ix], 0.2).tolist(), metrics=metrics, data_path=str(path), data_sha256=sha(path), stl_lower_sha256=lh, scope='Unsigned distance at sample points; all lower triangles screened, no signed penetration, no physical-force claim')
        dump(out, row)
        rows.append(row)
        print('3D clearance', case, 'parity', parity, 'query_s', round(qs, 2), flush=True)
    result = dict(round='R4', cases=list(range(1, max_cases + 1)), all_brute_parity_gate=all((x['parity_gate'] for x in rows)), maximum_parity_mm=max((x['parity_mm'] for x in rows)), central_band_projected_only_pixels=sum((x['metrics'][1]['projected_only_pixels'] for x in rows)), central_band_3D_only_pixels=sum((x['metrics'][1]['threeD_only_pixels'] for x in rows)), central_band_projected_pixels=sum((x['metrics'][1]['projected_contact_pixels'] for x in rows)), central_band_3D_pixels=sum((x['metrics'][1]['unsigned_3D_near_pixels'] for x in rows)), costs=dict(extra_preparation_s=sum((x['preparation_s'] for x in rows)), screened_queries_s=sum((x['screened_query_s'] for x in rows)), brute_checks_s=sum((x['brute_control_s'] for x in rows)), full_R1_R2_R3_discovery_and_preparation_charged=True), external_referent=dict(kind='closed_form', locator='https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf', compared_quantity='Unsigned point-to-triangle clearance mm; independently derived convex minimizer over triangle', refutes_us=True), outcome='FIXED_POSE_3D_CLEARANCE_OPERATOR; UNSIGNED_SAMPLED_ONLY; FORCE_UNKNOWN', rows=rows)
    dump(H / ('round4/results.json' if max_cases == 12 else 'round4/results_full200.json'), result)
    print(json.dumps({k: v for (k, v) in result.items() if k != 'rows'}, indent=2))
if __name__ == '__main__':
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 12)
