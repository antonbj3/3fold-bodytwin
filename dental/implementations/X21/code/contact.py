"""Tooth-resolved fixed-pose projected geometry. No physical or 3D certificate."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
os.environ.setdefault('NUMBA_CACHE_DIR', _release_expand('@DENTAL_WORK_ROOT@/X21/numba'))
import sys, json, time, zipfile, hashlib, datetime, resource, argparse
from pathlib import Path
import numpy as np
from numba import njit
P = Path(__file__).resolve().parents[1]
D = Path(_release_expand('@DENTAL_WORK_ROOT@/X21'))
X11 = P.parent / 'LANE_X11_TOOTH_SEG'
sys.path.insert(0, str(X11 / 'code'))
from morphology import mesh_stl
LABELS = Path(_release_expand('@DENTAL_WORK_ROOT@/X11/targets/Bite2Text'))
UP = [10 * q + t for q in [1, 2] for t in range(1, 9)]
LO = [10 * q + t for q in [4, 3] for t in range(1, 9)]

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x

def write(p, x):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(clean(x), ensure_ascii=False, allow_nan=False, indent=2) + '\n')
    t.replace(p)

def state(status, gate, next_op, **kw):
    write(P / 'CURRENT_WORK_STATE.json', dict(lane='X21-contact-map', status=status, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate=gate, next_operation=next_op, **kw))

@njit(cache=True)
def cross2(ax, ay, bx, by):
    return ax * by - ay * bx

@njit(cache=True)
def affine(T):
    n = len(T)
    out = np.empty((n, 3))
    valid = np.ones(n, np.bool_)
    for i in range(n):
        ax = T[i, 1, 0] - T[i, 0, 0]
        ay = T[i, 1, 1] - T[i, 0, 1]
        az = T[i, 1, 2] - T[i, 0, 2]
        bx = T[i, 2, 0] - T[i, 0, 0]
        by = T[i, 2, 1] - T[i, 0, 1]
        bz = T[i, 2, 2] - T[i, 0, 2]
        det = ax * by - ay * bx
        if abs(det) <= 1e-12:
            valid[i] = False
            out[i, :] = 0.0
            continue
        a = (az * by - bz * ay) / det
        b = (ax * bz - bx * az) / det
        out[i, 0] = a
        out[i, 1] = b
        out[i, 2] = T[i, 0, 2] - a * T[i, 0, 0] - b * T[i, 0, 1]
    return (out, valid)

@njit(cache=True)
def intersection_min(u, l, coeff):
    poly = np.empty((12, 2))
    tmp = np.empty((12, 2))
    for k in range(3):
        poly[k, 0] = u[k, 0]
        poly[k, 1] = u[k, 1]
    n = 3
    det = cross2(l[1, 0] - l[0, 0], l[1, 1] - l[0, 1], l[2, 0] - l[0, 0], l[2, 1] - l[0, 1])
    sgn = 1.0 if det > 0 else -1.0
    for e in range(3):
        ex = l[(e + 1) % 3, 0] - l[e, 0]
        ey = l[(e + 1) % 3, 1] - l[e, 1]
        m = 0
        if n == 0:
            return (np.inf, 0.0, 0.0)
        for k in range(n):
            prev = (k - 1) % n
            px = poly[prev, 0]
            py = poly[prev, 1]
            qx = poly[k, 0]
            qy = poly[k, 1]
            dp = sgn * cross2(ex, ey, px - l[e, 0], py - l[e, 1])
            dq = sgn * cross2(ex, ey, qx - l[e, 0], qy - l[e, 1])
            pin = dp >= -1e-12
            qin = dq >= -1e-12
            if pin != qin:
                den = dp - dq
                if abs(den) > 1e-30:
                    a = dp / den
                    tmp[m, 0] = px + a * (qx - px)
                    tmp[m, 1] = py + a * (qy - py)
                    m += 1
            if qin:
                tmp[m, 0] = qx
                tmp[m, 1] = qy
                m += 1
        n = m
        for k in range(n):
            poly[k, 0] = tmp[k, 0]
            poly[k, 1] = tmp[k, 1]
    best = np.inf
    x = 0.0
    y = 0.0
    for k in range(n):
        g = coeff[0] * poly[k, 0] + coeff[1] * poly[k, 1] + coeff[2]
        if g < best:
            best = g
            x = poly[k, 0]
            y = poly[k, 1]
    return (best, x, y)

@njit(cache=True)
def boxes(T):
    out = np.empty((len(T), 6))
    for i in range(len(T)):
        for a in range(3):
            out[i, a] = min(T[i, 0, a], T[i, 1, a], T[i, 2, a])
            out[i, a + 3] = max(T[i, 0, a], T[i, 1, a], T[i, 2, a])
    return out

@njit(cache=True)
def buckets(box, origin, h, shape):
    count = np.zeros(shape[0] * shape[1], np.int64)
    for i in range(len(box)):
        x0 = max(0, int(np.floor((box[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((box[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((box[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((box[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                count[x * shape[1] + y] += 1
    ptr = np.empty(len(count) + 1, np.int64)
    ptr[0] = 0
    for i in range(len(count)):
        ptr[i + 1] = ptr[i] + count[i]
    ids = np.empty(ptr[-1], np.int32)
    pos = ptr[:-1].copy()
    for i in range(len(box)):
        x0 = max(0, int(np.floor((box[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((box[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((box[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((box[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                k = x * shape[1] + y
                ids[pos[k]] = i
                pos[k] += 1
    return (ptr, ids)

@njit(cache=True)
def extrema(U, L, ul, ll, uc, lc, ub, lb, origin, h, shape, ptr, ids, prune=True):
    best = np.full((29, 49), np.inf)
    wi = np.full((29, 49, 2), -1, np.int32)
    xy = np.zeros((29, 49, 2))
    seen = np.full(len(L), -1, np.int32)
    checked = 0
    tested = 0
    for i in range(len(U)):
        a = ul[i]
        x0 = max(0, int(np.floor((ub[i, 0] - origin[0]) / h)))
        x1 = min(shape[0] - 1, int(np.floor((ub[i, 3] - origin[0]) / h)))
        y0 = max(0, int(np.floor((ub[i, 1] - origin[1]) / h)))
        y1 = min(shape[1] - 1, int(np.floor((ub[i, 4] - origin[1]) / h)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                bkt = x * shape[1] + y
                for k in range(ptr[bkt], ptr[bkt + 1]):
                    j = ids[k]
                    if seen[j] == i:
                        continue
                    seen[j] = i
                    b = ll[j]
                    if ub[i, 0] > lb[j, 3] or ub[i, 3] < lb[j, 0] or ub[i, 1] > lb[j, 4] or (ub[i, 4] < lb[j, 1]):
                        continue
                    checked += 1
                    if prune and ub[i, 2] - lb[j, 5] >= best[a, b]:
                        continue
                    (g, px, py) = intersection_min(U[i], L[j], uc[i] - lc[j])
                    tested += 1
                    if g < best[a, b]:
                        best[a, b] = g
                        wi[a, b, 0] = i
                        wi[a, b, 1] = j
                        xy[a, b, 0] = px
                        xy[a, b, 1] = py
    return (best, wi, xy, checked, tested)

@njit(cache=True)
def raster(T, lab, coeff, origin, h, shape, upper):
    z = np.full(shape[0] * shape[1], np.inf if upper else -np.inf)
    own = np.full(len(z), -1, np.int32)
    bb = boxes(T)
    for i in range(len(T)):
        det = cross2(T[i, 1, 0] - T[i, 0, 0], T[i, 1, 1] - T[i, 0, 1], T[i, 2, 0] - T[i, 0, 0], T[i, 2, 1] - T[i, 0, 1])
        sign = 1.0 if det > 0 else -1.0
        x0 = max(0, int(np.ceil((bb[i, 0] - origin[0]) / h - 0.5)))
        x1 = min(shape[0] - 1, int(np.floor((bb[i, 3] - origin[0]) / h - 0.5)))
        y0 = max(0, int(np.ceil((bb[i, 1] - origin[1]) / h - 0.5)))
        y1 = min(shape[1] - 1, int(np.floor((bb[i, 4] - origin[1]) / h - 0.5)))
        for x in range(x0, x1 + 1):
            for y in range(y0, y1 + 1):
                px = origin[0] + (x + 0.5) * h
                py = origin[1] + (y + 0.5) * h
                ok = True
                for e in range(3):
                    ex = T[i, (e + 1) % 3, 0] - T[i, e, 0]
                    ey = T[i, (e + 1) % 3, 1] - T[i, e, 1]
                    if sign * cross2(ex, ey, px - T[i, e, 0], py - T[i, e, 1]) < -1e-12:
                        ok = False
                        break
                if not ok:
                    continue
                zz = coeff[i, 0] * px + coeff[i, 1] * py + coeff[i, 2]
                k = x * shape[1] + y
                if upper and zz < z[k] or (not upper and zz > z[k]):
                    z[k] = zz
                    own[k] = i
    return (z, own)

def prepare(blob, meta, jaw, R, center):
    if hashlib.sha256(blob).hexdigest() != meta['input_sha256']:
        raise ValueError('SOURCE_HASH_MISMATCH')
    (v, f) = mesh_stl(blob, fast=True)
    label_path = Path(meta['labels']['path'])
    if sha(label_path) != meta['labels']['sha256']:
        raise ValueError('LABEL_HASH_MISMATCH')
    labs = np.load(label_path)['labels']
    if len(labs) != len(v):
        raise ValueError('LABEL_INDEX_LENGTH_MISMATCH')
    fl = labs[f]
    (a, b, c) = fl.T
    owner = np.where((a == b) | (a == c), a, np.where(b == c, b, 0)).astype(np.int32)
    allowed = UP if jaw == 'upper' else LO
    keep = np.isin(owner, allowed) & (owner > 0)
    faceids = np.flatnonzero(keep)
    T = np.ascontiguousarray(((v.astype(float) - center) @ R)[f[keep]])
    lab = owner[keep]
    (coeff, valid) = affine(T)
    stats = dict(vertices=len(v), faces=len(f), unassigned_faces=int((~keep).sum()), mixed_label_faces=int(np.any(fl != fl[:, 0, None], axis=1).sum()), degenerate_projections=int((~valid).sum()), labelled_face_fraction=float(keep.mean()))
    return (np.ascontiguousarray(T[valid]), lab[valid], coeff[valid], faceids[valid], stats)

def make_map(U, L, ul, ll, uc, lc, h=0.2, shift=(0.0, 0.0)):
    lo = np.minimum(U[:, :, :2].min((0, 1)), L[:, :, :2].min((0, 1)))
    hi = np.maximum(U[:, :, :2].max((0, 1)), L[:, :, :2].max((0, 1)))
    origin = np.floor(lo / h) * h - 2 * h + np.array(shift)
    shape = np.ceil((hi - origin) / h).astype(np.int64) + 3
    (uz, ui) = raster(U, ul, uc, origin, h, shape, True)
    (lz, li) = raster(L, ll, lc, origin, h, shape, False)
    ok = np.isfinite(uz) & np.isfinite(lz)
    ii = np.flatnonzero(ok)
    xy = origin + (np.c_[ii // shape[1], ii % shape[1]] + 0.5) * h
    return dict(gap=uz[ok] - lz[ok], xy=xy, upper_z=uz[ok], lower_z=lz[ok], upper_face=ui[ok], lower_face=li[ok], upper_fdi=ul[ui[ok]], lower_fdi=ll[li[ok]], origin=origin, shape=shape, grid_mm=h)

def summaries(m, best, wi, xy, U, L, uc, lc, uf, lf):
    rows = []
    features = {}
    h = m['grid_mm']
    gap = m['gap']
    u = m['upper_fdi']
    l = m['lower_fdi']
    features['total_signed_proximity_area_mm2'] = float(np.sum(gap <= 0.1) * h * h)
    features['total_absolute_near_area_mm2'] = float(np.sum(np.abs(gap) <= 0.1) * h * h)
    features['penetration_area_mm2'] = float(np.sum(gap < -0.03) * h * h)
    for a in UP:
        for b in LO:
            key = f'{a}_{b}'
            g = float(best[a, b])
            idx = (u == a) & (l == b)
            central = idx & (gap <= 0.1)
            ar = float(central.sum() * h * h)
            features['gap_' + key] = g if np.isfinite(g) else None
            features['area_' + key] = ar
            if not np.isfinite(g):
                continue
            (i, j) = wi[a, b]
            point = xy[a, b]
            w = dict(upper_source_face=int(uf[i]), lower_source_face=int(lf[j]), upper_triangle_frame_mm=U[i], lower_triangle_frame_mm=L[j], xy_frame_mm=point, upper_z_mm=float(uc[i, :2] @ point + uc[i, 2]), lower_z_mm=float(lc[j, :2] @ point + lc[j, 2]))
            if central.any():
                loc = m['xy'][central].mean(0)
                features['cx_' + key] = loc[0]
                features['cy_' + key] = loc[1]
            else:
                loc = None
                features['cx_' + key] = None
                features['cy_' + key] = None
            rows.append(dict(upper_fdi=a, lower_fdi=b, minimum_projected_gap_mm=g, numerical_tolerance_mm=1e-07, conditional_vertical_offset_interval_mm=[g - 0.05, g + 0.05], projected_proximity_area_mm2=ar, absolute_near_area_mm2=float(np.sum(idx & (np.abs(gap) <= 0.1)) * h * h), penetrating_projected_area_mm2=float(np.sum(idx & (gap < -0.03)) * h * h), projected_area_offset_scenario_mm2=[float(np.sum(idx & (gap <= 0.05)) * h * h), float(np.sum(idx & (gap <= 0.15)) * h * h)], area_uncertainty='offset scenario only; raster convergence separate; segmentation/pose error UNKNOWN', patch_centroid_xy_frame_mm=loc, near_contact_under_all_offset_scenarios=bool(g <= 0.05), near_contact_possible_under_offset_scenarios=bool(g <= 0.15), excess_projected_penetration=bool(g < -0.03), witness=w))
    for a in UP:
        ix = (u == a) & (gap <= 0.1)
        features['tooth_area_' + str(a)] = float(ix.sum() * h * h)
        features['tooth_degree_' + str(a)] = len(np.unique(l[ix]))
        for ax in range(2):
            features[f'tooth_centroid_{a}_{ax}'] = float(m['xy'][ix, ax].mean()) if ix.any() else None
    return (rows, features)

def load_case(z, case):
    meta = json.load(open(LABELS / case / 'labels+landmarks.json'))
    fr = meta['arches']['lower']['frame']
    R = np.column_stack([fr['right_unit'], fr['anterior_unit'], fr['superior_unit']])
    center = np.array(fr['center_mm'])
    if not np.allclose(R.T @ R, np.eye(3), atol=1e-10):
        raise ValueError('NON_ORTHONORMAL_FRAME')
    arches = []
    for jaw in ['upper', 'lower']:
        arches.append(prepare(z.read(case + '/ios/ios_' + jaw + '.stl'), meta['arches'][jaw], jaw, R, center))
    return (meta, R, center, arches)

def run_case(z, case, controls=False):
    st = time.perf_counter()
    cpu = time.process_time()
    (meta, R, center, arches) = load_case(z, case)
    (U, ul, uc, uf, us) = arches[0]
    (L, ll, lc, lf, ls) = arches[1]
    decode = time.perf_counter() - st
    ub = boxes(U)
    lb = boxes(L)
    origin = np.minimum(ub[:, :2].min(0), lb[:, :2].min(0)) - 0.8
    hi = np.maximum(ub[:, 3:5].max(0), lb[:, 3:5].max(0))
    shape = np.ceil((hi - origin) / 0.8).astype(np.int64) + 2
    (ptr, ids) = buckets(lb, origin, 0.8, shape)
    s = time.perf_counter()
    (best, wi, xy, checks, tests) = extrema(U, L, ul, ll, uc, lc, ub, lb, origin, 0.8, shape, ptr, ids, True)
    qs = time.perf_counter() - s
    s = time.perf_counter()
    m = make_map(U, L, ul, ll, uc, lc)
    patch_s = time.perf_counter() - s
    (rows, features) = summaries(m, best, wi, xy, U, L, uc, lc, uf, lf)
    output = D / 'maps' / (case + '.npz')
    np.savez_compressed(output, **m, minimum_gap=best, witness_faces=wi, witness_xy=xy)
    row = dict(case_id=case, claim_type='capability', schema='X21-projected-contact-map-v1', geometry_model='Conditional continuous projected affine minimum; floating geometry, not formal 3D or physical certificate', source_zip_members={jaw: dict(member=case + '/ios/ios_' + jaw + '.stl', sha256=meta['arches'][jaw]['input_sha256'], labels_sha256=meta['arches'][jaw]['labels']['sha256']) for jaw in ['upper', 'lower']}, frame=dict(R_columns=R, center_xyz_mm=center, anatomical_status='UNKNOWN_TARGET_FRAME'), segmentation_status='Predicted X11 FDI; Bite2Text correctness UNKNOWN', stats=dict(upper=us, lower=ls, upper_projected_triangles=len(U), lower_projected_triangles=len(L), broad_phase_aabb_pairs=checks, exact_pairs_evaluated=tests), pairs=rows, features=features, map_path=str(output), map_sha256=sha(output), landmark_source_sha256=sha(LABELS / case / 'labels+landmarks.json'), cost=dict(decode_s=decode, minimum_query_s=qs, patch_query_s=patch_s, wall_s=time.perf_counter() - st, cpu_s=time.process_time() - cpu, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024), error=None)
    if controls:
        s = time.perf_counter()
        finite = np.argwhere(np.isfinite(best))
        (a, b) = finite[0]
        ui = np.flatnonzero(ul == a)
        lj = np.flatnonzero(ll == b)
        uu = U[ui]
        lt = L[lj]
        uub = ub[ui]
        llb = lb[lj]
        (uo, lo) = (uc[ui], lc[lj])
        labu = np.full(len(uu), a, np.int32)
        labl = np.full(len(lt), b, np.int32)
        (p0, k0) = buckets(llb, origin, 0.8, shape)
        (oracle, _, _, c0, t0) = extrema(uu, lt, labu, labl, uo, lo, uub, llb, origin, 0.8, shape, p0, k0, False)
        parity = float(abs(oracle[a, b] - best[a, b]))
        pair_oracle_s = time.perf_counter() - s
        import importlib.util, types
        fake = types.ModuleType('occlusion')
        fake.HERE = P
        fake.DATA = D
        for k in ['TEETH', 'UP', 'LO', 'transformed', 'maps', 'sha', 'write_state']:
            setattr(fake, k, None)
        sys.modules['occlusion'] = fake
        spec = importlib.util.spec_from_file_location('inherited_continuous', P.parent / 'LANE_NEXT_P_OCCLUSION_VALIDATION/code/continuous_gap.py')
        cg = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cg)
        (i, j) = wi[a, b]
        (val, _, _, _) = cg.exact_batch(U, L, uc, lc, np.array([[i, j]], np.int32))
        werr = float(abs(val - best[a, b]))
        ixs = np.linspace(0, len(uu) - 1, min(300, len(uu)), dtype=int)
        jxs = np.linspace(0, len(lt) - 1, min(300, len(lt)), dtype=int)
        u0 = uu[ixs]
        l0 = lt[jxs]
        cu = uo[ixs]
        cl = lo[jxs]
        prod = np.stack(np.meshgrid(np.arange(len(u0)), np.arange(len(l0)), indexing='ij'), -1).reshape(-1, 2)
        truth = np.inf
        for f in range(0, len(prod), 4096):
            (v, _, _, _) = cg.exact_batch(u0, l0, cu, cl, prod[f:f + 4096])
            truth = min(truth, v)
        ub0 = boxes(u0)
        lb0 = boxes(l0)
        (p1, k1) = buckets(lb0, origin, 0.8, shape)
        (bm, _, _, _, _) = extrema(u0, l0, labu[:len(u0)], labl[:len(l0)], cu, cl, ub0, lb0, origin, 0.8, shape, p1, k1, True)
        subset_parity = float(abs(truth - bm[a, b])) if np.isfinite(truth) else 0.0 if not np.isfinite(bm[a, b]) else np.inf
        fine = make_map(U, L, ul, ll, uc, lc, h=0.1)
        ar = float((m['gap'] <= 0.1).sum() * 0.04)
        far = float((fine['gap'] <= 0.1).sum() * 0.01)
        delta = abs(ar - far)
        passed = delta <= max(1.0, 0.2 * far)
        row['numerical_control'] = dict(pair=[int(a), int(b)], triangle_counts=[len(ui), len(lj)], full_aabb_pairs=int(c0), all_pair_intersections_evaluated=int(t0), oracle_minimum_mm=float(oracle[a, b]), parity_mm=parity, witness_inherited_parity_mm=werr, full_subset_product=len(prod), full_subset_parity_mm=subset_parity, gate=bool(max(parity, werr, subset_parity) <= 1e-07), injected_plus_1mm_parity=abs(float(oracle[a, b]) - (best[a, b] + 1.0)), injected_rejected=bool(abs(float(oracle[a, b]) - (best[a, b] + 1.0)) > 1e-07), full_pair_query_s=pair_oracle_s, total_control_s=time.perf_counter() - s, area_refinement=dict(coarse_mm2=ar, fine_mm2=far, absolute_delta_mm2=delta, tolerance_mm2=max(1.0, 0.2 * far), gate=passed), area_injected_100mm2_rejected=bool(abs(ar + 100 - far) > max(1.0, 0.2 * far)))
    return row

def batch(limit=None, panel_only=False):
    manifest = json.load(open(P / 'raw/INPUT_MANIFEST.json'))
    cases = manifest['numerical_panel'] if panel_only else manifest['cases']
    cases = cases[:limit] if limit else cases
    done = {p.stem for p in (P / 'raw/cases').glob('*.json')} if (P / 'raw/cases').exists() else set()
    start = time.perf_counter()
    with zipfile.ZipFile(manifest['zip']) as z:
        for (n, c) in enumerate(cases):
            if c in done:
                continue
            state('R1_MAPPING', 'PREREG frozen; no report evaluation', f'Continuous projected minima and patches {c}', completed=len(done), expected=len(manifest['cases']))
            try:
                r = run_case(z, c, c in manifest['numerical_panel'])
            except Exception as e:
                import traceback
                r = dict(case_id=c, error=repr(e), traceback=traceback.format_exc())
            write(P / 'raw/cases' / (c + '.json'), r)
            done.add(c)
            print(c, 'elapsed', round(time.perf_counter() - start, 1), 's', 'query', r.get('cost', {}).get('minimum_query_s'), 'error', r['error'], flush=True)
            if n % 20 == 0:
                size = sum((p.stat().st_size for p in D.rglob('*') if p.is_file()))
                if size > 3000000000:
                    raise RuntimeError('X21 3GB limit exceeded')
    state('R1_MAP_BATCH_COMPLETED', 'Geometry output preserved; report gates pending', 'Freeze report models and predictions', completed=len(done), expected=len(manifest['cases']))
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int)
    ap.add_argument('--panel-only', action='store_true')
    a = ap.parse_args()
    batch(a.limit, a.panel_only)
