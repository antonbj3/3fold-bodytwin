"""Continuous regional affine minima. Patch areas remain an explicit closure."""
from dental_release.paths import expand as _release_expand
import os
os.environ['NUMBA_CACHE_DIR'] = _release_expand('@DENTAL_WORK_ROOT@/X54/numba')
import copy, json, time, zipfile, datetime
import numpy as np
from numba import njit
from scipy.optimize import linprog
from contact_model import P, D, BASE, ZIP, LABELS, sha, write, state, load_module, mesh_stl, solve
from run_r2 import force_control
cg = load_module('x54_x21_geometry', BASE / 'LANE_X21_CONTACT_MAP/code/contact.py')

def clip(poly, direction, center, sign):
    out = []
    for (i, q) in enumerate(poly):
        p = poly[i - 1]
        a = sign * ((p[:2] - center) @ direction)
        b = sign * ((q[:2] - center) @ direction)
        if (a >= 0) != (b >= 0):
            out.append(p + a / (a - b) * (q - p))
        if b >= 0:
            out.append(q)
    return out

def decoded(g):
    case = g['case']
    R = np.asarray(g['frame']['R'])
    center = np.asarray(g['frame']['center'])
    out = []
    with zipfile.ZipFile(ZIP) as z:
        for s in g['sources']:
            blob = z.read(s['zip_member'])
            (v, f) = mesh_stl(blob, fast=True)
            labs = np.load(s['label_path'])['labels']
            (a, b, c) = labs[f].T
            own = np.where((a == b) | (a == c), a, np.where(b == c, b, 0))
            allowed = [int(t) for t in g['centers'][s['jaw']]]
            keep = np.isin(own, allowed)
            T = ((v.astype(float) - center) @ R)[f[keep]]
            lab = own[keep].astype(np.int32)
            ids = np.flatnonzero(keep)
            norm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
            norm /= np.maximum(np.linalg.norm(norm, axis=1, keepdims=True), 1e-30)
            norm *= np.where(norm[:, 2] < 0, -1, 1)[:, None]
            ok = norm[:, 2] >= 0.35
            T = T[ok]
            lab = lab[ok]
            norm = norm[ok]
            ids = ids[ok]
            (co, valid) = cg.affine(T)
            out.append((T[valid], lab[valid], norm[valid], ids[valid], co[valid]))
    return out

@njit(cache=True)
def minima(U, L, ul, ll, un, ln, uc, lc, ub, lb, origin, h, shape, ptr, ids, prune):
    best = np.full((29, 49), np.inf)
    w = np.full((29, 49, 2), -1, np.int32)
    xy = np.zeros((29, 49, 2))
    seen = np.full(len(L), -1, np.int32)
    checked = 0
    for i in range(len(U)):
        a = ul[i]
        x0 = max(0, int(np.floor((ub[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((ub[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((ub[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((ub[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                cell = x * shape[1] + y
                for k in range(ptr[cell], ptr[cell + 1]):
                    j = ids[k]
                    if seen[j] == i:
                        continue
                    seen[j] = i
                    b = ll[j]
                    if ub[i, 0] > lb[j, 3] or ub[i, 3] < lb[j, 0] or ub[i, 1] > lb[j, 4] or (ub[i, 4] < lb[j, 1]):
                        continue
                    if np.dot(un[i], ln[j]) < 0.5:
                        continue
                    if prune and ub[i, 2] - lb[j, 5] >= best[a, b]:
                        continue
                    (value, px, py) = cg.intersection_min(U[i], L[j], uc[i] - lc[j])
                    checked += 1
                    if value < best[a, b]:
                        best[a, b] = value
                        w[a, b, 0] = i
                        w[a, b, 1] = j
                        xy[a, b, 0] = px
                        xy[a, b, 1] = py
    return (best, w, xy, checked)

def continuous(g):
    start = time.perf_counter()
    ((U, ul, un, uf, uc), (L, ll, ln, lf, lc)) = decoded(g)
    lb = cg.boxes(L)
    origin = np.minimum(U[:, :, :2].min((0, 1)), L[:, :, :2].min((0, 1))) - 0.8
    hi = np.maximum(U[:, :, :2].max((0, 1)), L[:, :, :2].max((0, 1)))
    shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
    (ptr, ids) = cg.buckets(lb, origin, 0.8, shape)
    out = copy.deepcopy(g)
    controls = []
    cost = []
    missing = []
    for region in range(4):
        ts = []
        labs = []
        normals = []
        src = []
        for (T, t, n, f) in zip(U, ul, un, uf):
            c = np.asarray(g['centers']['upper'][str(t)])[:2]
            rad = np.array([c[0], c[1] + 15.0])
            rad /= max(np.linalg.norm(rad), 1e-12)
            mes = np.array([-rad[1], rad[0]])
            pol = clip(list(T), rad, c, 1 if region // 2 else -1)
            if pol:
                pol = clip(pol, mes, c, 1 if region % 2 else -1)
            for k in range(1, len(pol) - 1):
                ts.append([pol[0], pol[k], pol[k + 1]])
                labs.append(t)
                normals.append(n)
                src.append(f)
        T = np.ascontiguousarray(ts, dtype=float)
        lab = np.array(labs, np.int32)
        nn = np.array(normals)
        ss = np.array(src)
        (co, valid) = cg.affine(T)
        T = T[valid]
        lab = lab[valid]
        nn = nn[valid]
        ss = ss[valid]
        co = co[valid]
        ub = cg.boxes(T)
        tick = time.perf_counter()
        (best, w, xy, ncheck) = minima(T, L, lab, ll, nn, ln, co, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
        cost.append(dict(region=region, upper_subtriangles=len(T), pairs_evaluated=ncheck, query_s=time.perf_counter() - tick))
        for p in out['patches']:
            if p['region'] != region:
                continue
            (a, b) = (p['upper_fdi'], p['lower_fdi'])
            (i, j) = w[a, b]
            if i < 0:
                missing.append((a, b, region))
                continue
            point = xy[a, b]
            n = nn[i] + ln[j]
            n /= np.linalg.norm(n)
            t = np.array([1.0, 0, 0])
            t -= t @ n * n
            t /= np.linalg.norm(t)
            p.update(gap_mm=float(best[a, b] * n[2]), vertical_gap_mm=float(best[a, b]), basis=np.vstack([n, t, np.cross(n, t)]).tolist(), upper_point_mm=np.r_[point, co[i, :2] @ point + co[i, 2]].tolist(), lower_point_mm=np.r_[point, lc[j, :2] @ point + lc[j, 2]].tolist(), upper_source_face=int(ss[i]), lower_source_face=int(lf[j]), continuous_upper_subtriangle=T[i].tolist(), continuous_lower_triangle=L[j].tolist(), minimum_scope='continuous affine projected model, regional halfplane clipping; float64, no formal interval certificate')
            cobj = np.r_[T[i, :, 2], -L[j, :, 2]]
            A = np.zeros((4, 6))
            A[0, :3] = 1
            A[1, 3:] = 1
            A[2, :3] = T[i, :, 0]
            A[2, 3:] = -L[j, :, 0]
            A[3, :3] = T[i, :, 1]
            A[3, 3:] = -L[j, :, 1]
            lp = linprog(cobj, A_eq=A, b_eq=[1, 1, 0, 0], bounds=(0, None), method='highs', options={'threads': 4})
            error = abs(lp.fun - best[a, b]) if lp.success else np.inf
            controls.append(dict(upper_fdi=a, lower_fdi=b, region=region, minimum_gap_mm=float(best[a, b]), independent_LP_gap_mm=float(lp.fun) if lp.success else None, error_mm=float(error), pass_gate=bool(error <= 1e-08), injected_001mm_error_rejected=bool(abs(best[a, b] + 0.001 - lp.fun) > 1e-08) if lp.success else False))
        regional_targets = [p for p in out['patches'] if p['region'] == region]
        if region == 0 and regional_targets:
            target = regional_targets[0]
            (a, b) = (target['upper_fdi'], target['lower_fdi'])
            iu = np.flatnonzero(lab == a)
            il = np.flatnonzero(ll == b)
            lsmall = L[il]
            bsmall = lb[il]
            (p0, i0) = cg.buckets(bsmall, origin, 0.8, shape)
            (true, _, _, count) = minima(T[iu], lsmall, lab[iu], ll[il], nn[iu], ln[il], co[iu], lc[il], ub[iu], bsmall, origin, 0.8, shape, p0, i0, False)
            controls.append(dict(kind='unpruned_regional_pair', upper_fdi=a, lower_fdi=b, region=region, error_mm=float(abs(true[a, b] - best[a, b])), pass_gate=bool(abs(true[a, b] - best[a, b]) <= 1e-08), pairs_evaluated=count))
    if missing:
        raise ValueError('MISSING_CONTINUOUS_WINNERS:' + str(missing))
    return (out, dict(controls=controls, cost=cost, wall_s=time.perf_counter() - start, area_status='original sampled patch area unchanged; finite elastic radius closure'))

def difference(a, b):
    return max((abs(a['force_shares_pp'][jaw].get(t, 0) - b['force_shares_pp'][jaw].get(t, 0)) for jaw in ('upper', 'lower') for t in set(a['force_shares_pp'][jaw]) | set(b['force_shares_pp'][jaw])))

def main():
    tick = time.perf_counter()
    rows = []
    manifest = []
    for case in [1, 2, 3]:
        state('R3_RUNNING', 'Raster force-share convergence failed in case1; PREREG_R3 frozen', 'Continuous source-triangle regional minimum case' + str(case))
        outfile = P / 'rounds/R3' / f'case{case:03d}.json'
        if outfile.exists():
            rows.append(json.loads(outfile.read_text()))
            continue
        g = json.loads((P / 'raw' / f'geometry_{case:03d}_h02.json').read_text())
        fine = json.loads((P / 'raw' / f'geometry_{case:03d}_h01.json').read_text())
        (con, info) = continuous(g)
        path = P / 'raw' / f'geometry_{case:03d}_continuous.json'
        write(path, con)
        (s, a) = solve(con, closure=0.05)
        ap = D / f'R3_case{case:03d}.npz'
        np.savez_compressed(ap, **a)
        (fs, fa) = force_control(con, 100.0)
        r1 = json.loads((P / 'rounds/R1' / f'case{case:03d}_d0.05.json').read_text())
        fmap = {(p['upper_fdi'], p['lower_fdi'], p['region']): p for p in fine['patches']}
        c2 = copy.deepcopy(con)
        omitted = []
        for p in c2['patches']:
            key = (p['upper_fdi'], p['lower_fdi'], p['region'])
            if key in fmap:
                p['area_mm2'] = fmap[key]['area_mm2']
            else:
                omitted.append(key)
        (sn, an) = solve(c2, closure=0.05)
        diff = difference(s, sn)
        row = dict(case=case, continuous_geometry_file=str(path), continuous_geometry_sha256=sha(path), checks=info, prediction=s, force_control_prediction=fs, raster_prediction_difference_pp=difference(s, r1), remaining_area_difference_pp=diff, remaining_area_gate=diff <= 2.0, omitted_regions=omitted, new_fine_regions=len(set(fmap) - {(p['upper_fdi'], p['lower_fdi'], p['region']) for p in con['patches']}), continuous_gap_shift_max_mm=max((abs(p['vertical_gap_mm'] - q['vertical_gap_mm']) for (p, q) in zip(con['patches'], g['patches']))), arrays=dict(path=str(ap), sha256=sha(ap)), scope='Named original regional patch set only; additional fine regions/area and physical model still incomplete')
        write(outfile, row)
        rows.append(row)
        manifest.append(dict(path=str(outfile), sha256=sha(outfile), arrays=row['arrays']))
        print('R3 case', case, 'LPmax', max((x['error_mm'] for x in info['controls'])), 'forcechange', row['raster_prediction_difference_pp'], 'remainingArea', diff, flush=True)
    fp = P / 'FROZEN_PREDICTIONS_R3.json'
    if not fp.exists():
        write(fp, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R3.json'), manifest=manifest, physical_measurement_status='NOT_RUN'))
        (P / 'FROZEN_PREDICTIONS_R3.json.sha256').write_text(sha(fp) + '\n')
    assert len(rows) > 0 and all((len(r['checks']['controls']) > 0 for r in rows)), 'EMPTY_CONTROL_SELECTION'
    result = dict(round='R3', claim_type='capability', answer='Continuous affine gap witnesses delivered at named surface regions; sampled minimum is removed, remaining area and topology debt exposed', rows=rows, all_LP_controls_pass=all((x['pass_gate'] for r in rows for x in r['checks']['controls'])), maximum_LP_error_mm=max((x['error_mm'] for r in rows for x in r['checks']['controls'])), all_remaining_area_gates_pass=all((r['remaining_area_gate'] for r in rows)), original_refinement_failed_gate_preserved=True, external_referent=dict(kind='closed_form', locator='barycentric affine minimum; scipy.optimize.linprog official implementation; same-triangle LP constraints stored in code/run_r3.py', compared_quantity='continuous projected regional gap minimum mm', refutes_us=True), physical_external_referent='Ferrato/Hattori R1/R2 force comparison remains failed/UNKNOWN; geometry LP is not patient-force validation', cost=dict(wall_s=time.perf_counter() - tick), rigorous_enclosure='MISSING: continuous affine evaluation float64; fixed patch area and topology, reduced elasticity and empirical uncertainty not enclosed')
    write(P / 'rounds/R3/results.json', result)
    state('R3_COMPLETE', 'Continuous geometry LP checks recorded; physical force validation still UNKNOWN', 'Package demo and frozen independent measurement port; local force/preload measurement missing')
if __name__ == '__main__':
    main()
