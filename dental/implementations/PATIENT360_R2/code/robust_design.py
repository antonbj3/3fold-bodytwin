"""Witness-based exact infeasibility for ONE declared vertical-offset family.

An exact local overlap gives a necessary nonpenetration limit. Every source
facet containing a contact query point is enumerated with rational predicates.
Thus rejection needs neither a floating global minimum nor a linearized bound.
"""
from fractions import Fraction as Q
import sys, time
from common import *
from region_field import *

def qpoint(p):
    return tuple((Q(float(v)) for v in p))

def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]

def sub(a, b):
    return tuple((x - y for (x, y) in zip(a, b)))

def det(a, b, c):
    return cross(sub(b, a), sub(c, a))

def qheight(tri, p):
    (a, b, c) = tri
    den = det(a, b, c)
    if not den:
        return None
    u = det(a, p, c) / den
    v = det(a, b, p) / den
    if u < 0 or v < 0 or u + v > 1:
        return None
    return a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2])

def overlap(a, b):
    poly = [p[:2] for p in a]
    ori = 1 if det(*b) > 0 else -1
    for j in range(3):
        (c, d) = (b[j], b[(j + 1) % 3])
        new = []
        for (k, q) in enumerate(poly):
            p = poly[k - 1]
            sp = ori * det(c, d, p)
            sq = ori * det(c, d, q)
            if (sp >= 0) != (sq >= 0):
                t = sp / (sp - sq)
                new.append(tuple((x + t * (y - x) for (x, y) in zip(p, q))))
            if sq >= 0:
                new.append(q)
        poly = new
        if not poly:
            break
    return poly

def exact_heights(tri, xy):
    lo = tri[:, :, :2].min(1)
    hi = tri[:, :, :2].max(1)
    out = []
    owners = []
    for p in xy:
        ids = np.flatnonzero(np.all(lo <= p, axis=1) & np.all(hi >= p, axis=1))
        q = qpoint(p)
        vals = []
        for i in ids:
            value = qheight(tuple((qpoint(v) for v in tri[i])), q)
            if value is not None:
                vals.append((value, int(i)))
        if not vals:
            out.append(None)
            owners.append(-1)
        else:
            (value, i) = min(vals)
            out.append(value)
            owners.append(i)
    return (out, owners)

def run(out, fields):
    t0 = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    field = next((f for f in fields if f.meta['field_id'] == 'Bite2Text_F4775_upper'))
    path = PARENT_DATA / 'Bite2Text_F4775/roof.npz'
    source(path)
    z = np.load(path)
    pr = load(ROOT / 'PREREG_C4_NATIVE_AND_ROBUST_DESIGN.json')
    (lo, hi) = pr['metrics']['pose_interval_mm']
    b = pr['metrics']['contact_band_mm']
    U = field.tri
    sel = np.all(U[:, :, :2].max(1) >= z['xy'].min(0), axis=1) & np.all(U[:, :, :2].min(1) <= z['xy'].max(0), axis=1)
    source_ids = np.flatnonzero(sel)
    U = U[sel]
    area = np.cross(U[:, 1, :2] - U[:, 0, :2], U[:, 2, :2] - U[:, 0, :2])
    valid = np.abs(area) > 1e-12
    source_ids = source_ids[valid]
    U = U[valid]
    sys.path.insert(0, str(RESULTS / 'LANE_X18_CROWN_ANTAGONIST/code'))
    import geometry as g
    coarse = g.signed_gap(U, z['xy'], z['z'], z['faces'])
    (upper_id, crown_id) = coarse['witness_triangle_indices']
    C = np.c_[z['xy'], z['z']][z['faces']]
    qa = tuple((qpoint(p) for p in U[upper_id]))
    qb = tuple((qpoint(p) for p in C[crown_id]))
    poly = overlap(qa, qb)
    candidates = [(qheight(qa, p) - qheight(qb, p), p) for p in poly]
    if not candidates:
        raise ValueError('EXACT_OVERLAP_WITNESS_MISSING')
    (witness_gap, witness_xy) = min(candidates)
    (heights, owners) = exact_heights(U, z['xy'])
    gaps = [v - Q(float(zz)) if v is not None else None for (v, zz) in zip(heights, z['z'])]
    finite = [(v, i) for (i, v) in enumerate(gaps) if v is not None]
    (minimum_point_gap, min_point) = min(finite)
    hmax = witness_gap + Q(lo)
    hmin = minimum_point_gap + Q(hi) - Q(b)
    deficit = hmin - hmax
    infeasible = hmin > hmax
    native_err = max((abs(float(v) - float(z['ceiling'][i])) for (i, v) in enumerate(heights) if v is not None))
    cert = dict(patient_id='Bite2Text_F4775', FDI=36, claim_type='capability', offset_definition='h>0 raises the entire fixed crown; source anatomy unchanged', pose_interval_mm=[lo, hi], nonpenetration_necessary_h_upper_mm=float(hmax), any_point_contact_necessary_h_lower_mm=float(hmin), exact_h_upper=str(hmax), exact_h_lower=str(hmin), incompatibility_mm=float(deficit), exact_deficit=str(deficit), status='NO_RIGID_VERTICAL_OFFSET_SATISFIES_BOTH' if infeasible else 'UNKNOWN_NEEDS_FULL_FEASIBILITY', witness=dict(source_facet=int(source_ids[upper_id]), crown_facet=int(crown_id), exact_gap=str(witness_gap), exact_xy=[str(v) for v in witness_xy], source_triangle=U[upper_id], crown_triangle=C[crown_id]), point_gap_witness=dict(point_id=min_point, exact_gap=str(minimum_point_gap), source_facet=int(source_ids[owners[min_point]])), exact_point_gap_field=[str(v) if v is not None else None for v in gaps], source=field.meta['source'], source_frame=field.meta['native_to_patient'], resolution='PER_POINT', time_scale='SIMULTANEOUS', enclosure='Exact rational inequalities on binary64 transformed source/crown geometry. Affine delta-h relation exact, zero remainder; no source-to-anatomy/FE enclosure claimed.', physical_measurement='NOT_RUN', physical_status='UNKNOWN', binding_quantity='Actual loaded relative pose/compliance; rigid offset cannot repair the declared geometric conflict', next_construction='Measured relative pose narrow enough for both constraints, or local surface deformation with independent preparation/load validation')
    freeze(out / 'FROZEN_PREDICTIONS.json', dict(prereg=source(ROOT / 'PREREG_C4_NATIVE_AND_ROBUST_DESIGN.json'), **cert))
    (h_ref, _) = g.query_height(U, z['xy'], True)
    referr = max((abs(float(v) - h_ref[i]) for (i, v) in enumerate(heights) if v is not None))
    direct_at_hmin = g.signed_gap(U + [0, 0, lo], z['xy'], z['z'] + float(hmin), z['faces'])
    direct_at_hmax = float(min((v + Q(hi) - hmax for v in gaps if v is not None)))

    def verify_offsets(candidate):
        if Q(candidate['exact_h_upper']) != witness_gap + Q(lo):
            raise ValueError('NONPENETRATION_LIMIT_CHANGED')
        if Q(candidate['exact_h_lower']) != minimum_point_gap + Q(hi) - Q(b):
            raise ValueError('CONTACT_LIMIT_CHANGED')
        if candidate['status'].startswith('NO_') and (not Q(candidate['exact_h_lower']) > Q(candidate['exact_h_upper'])):
            raise ValueError('NO_INCOMPATIBILITY')
        return True
    bad = dict(cert)
    bad['exact_h_upper'] = str(hmax + Q(1))
    checks = [check('native_source_height', native_err < 1e-06, native_err + 0.01 > 1e-06, dict(max_error_mm=native_err, injection='Add 0.01 mm to ceiling')), check('rational_height_vs_X18', referr < 1e-06, referr + 0.01 > 1e-06, dict(max_error_mm=referr, injection='Add 0.01 mm')), check('exact_rigid_offset_certificate', verify_offsets(cert), rejects(verify_offsets, bad), dict(injection='Raise h_upper by 1 mm')), check('direct_infeasibility_endpoints', infeasible and direct_at_hmin['minimum_gap_mm'] < 0 and (direct_at_hmax > b), not (0.0 < 0 and direct_at_hmax > b), dict(contact_offset_minimum_gap_mm=direct_at_hmin['minimum_gap_mm'], nonpenetrating_offset_minimum_point_gap_mm=direct_at_hmax, injection='Replace negative endpoint penetration by zero'))]
    result = dict(**cert, native_height_max_error_mm=native_err, checks=checks, external_referent=dict(kind='published_dataset', locator=field.meta['source']['archive'] + '::' + field.meta['source']['member'], compared_quantity='Original IOS antagonist heights and exact source facets in geometric offset infeasibility witness', refutes_us=not infeasible or native_err >= 1e-06), runtime_seconds=time.perf_counter() - t0)
    dump(out / 'ROBUST_DESIGN_RESULTS.json', result)
    return result
