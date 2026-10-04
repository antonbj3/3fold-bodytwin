"""Separate exact verifier: Gaussian elimination, no pose/design module imports.

Checks every point against all original source triangles selected by exact
coordinate comparisons, with no floating area cutoff. This can refute the
producer's near-degenerate-facet screening or overlap-witness calculation.
"""
from fractions import Fraction as F
import copy, time
import numpy as np
from common import dump, check, rejects, PARENT_DATA

def barycentric(triangle, point):
    t = [[F(float(x)) for x in vertex] for vertex in triangle]
    p = [F(x) if isinstance(x, (str, F)) else F(float(x)) for x in point]
    rows = [[t[j][i] for j in range(3)] + [p[i]] for i in range(2)] + [[F(1), F(1), F(1), F(1)]]
    for j in range(3):
        piv = next((i for i in range(j, 3) if rows[i][j]), None)
        if piv is None:
            return None
        (rows[j], rows[piv]) = (rows[piv], rows[j])
        den = rows[j][j]
        rows[j] = [x / den for x in rows[j]]
        for i in range(3):
            if i != j:
                mult = rows[i][j]
                rows[i] = [x - mult * y for (x, y) in zip(rows[i], rows[j])]
    weights = [row[3] for row in rows]
    if any((w < 0 for w in weights)):
        return None
    return sum((weights[i] * t[i][2] for i in range(3)), F(0))

def verify_certificate(cert, field, roof):
    w = cert['witness']
    source_tri = field.tri[w['source_facet']]
    crown_tri = np.c_[roof['xy'], roof['z']][roof['faces'][w['crown_facet']]]
    if not np.array_equal(source_tri, w['source_triangle']) or not np.array_equal(crown_tri, w['crown_triangle']):
        raise ValueError('CERTIFICATE_SOURCE_TRIANGLE_MISMATCH')
    a = barycentric(source_tri, w['exact_xy'])
    b = barycentric(crown_tri, w['exact_xy'])
    if a is None or b is None:
        raise ValueError('WITNESS_OUTSIDE_TRIANGLE')
    gap = a - b
    if gap != F(w['exact_gap']):
        raise ValueError('EXACT_WITNESS_GAP_MISMATCH')
    upper = gap + F(float(cert['pose_interval_mm'][0]))
    if upper != F(cert['exact_h_upper']):
        raise ValueError('CERTIFICATE_UPPER_MISMATCH')
    tri = field.tri
    lo = tri[:, :, :2].min(1)
    hi = tri[:, :, :2].max(1)
    gaps = []
    candidate_total = 0
    for (p, z) in zip(roof['xy'], roof['z']):
        candidates = np.flatnonzero(np.all(lo <= p, axis=1) & np.all(hi >= p, axis=1))
        candidate_total += len(candidates)
        vals = [v for i in candidates if (v := barycentric(tri[i], p)) is not None]
        gaps.append(min(vals) - F(float(z)) if vals else None)
    expected = [F(x) if x is not None else None for x in cert['exact_point_gap_field']]
    if gaps != expected:
        raise ValueError('EXACT_POINT_GAP_FIELD_MISMATCH')
    lower = min((g for g in gaps if g is not None)) + F(float(cert['pose_interval_mm'][1])) - F(0.1)
    if lower != F(cert['exact_h_lower']):
        raise ValueError('CERTIFICATE_LOWER_MISMATCH')
    if not lower > upper:
        raise ValueError('NO_INCOMPATIBILITY_CERTIFICATE')
    return dict(status='VERIFIED_EXACT_DIGITAL_INCOMPATIBILITY', points=len(gaps), source_candidate_triangle_queries=candidate_total, excluded_near_degenerate_facets=0, exact_gap_identity_error=0, exact_point_gap_identity_error=0, exact_upper=str(upper), exact_lower=str(lower), exact_deficit=str(lower - upper), scope='Exact binary digital source geometry; no physical/FE error enclosure')

def run(out, cert, fields):
    start = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    field = next((f for f in fields if f.meta['field_id'] == 'Bite2Text_F4775_upper'))
    roof = np.load(PARENT_DATA / 'Bite2Text_F4775/roof.npz')
    result = verify_certificate(cert, field, roof)
    checks = []
    mutations = []
    wrong = copy.deepcopy(cert)
    wrong['exact_h_upper'] = str(F(wrong['exact_h_upper']) + 1)
    mutations.append(('upper_plus_1mm', wrong))
    wrong = copy.deepcopy(cert)
    wrong['witness']['source_triangle'] = np.asarray(wrong['witness']['source_triangle']).copy()
    wrong['witness']['source_triangle'][0, 2] += 0.01
    mutations.append(('source_vertex_plus_0.01mm', wrong))
    wrong = copy.deepcopy(cert)
    i = next((i for (i, x) in enumerate(wrong['exact_point_gap_field']) if x is not None))
    wrong['exact_point_gap_field'][i] = str(F(wrong['exact_point_gap_field'][i]) + F(0.01))
    mutations.append(('point_gap_plus_0.01mm', wrong))
    for (name, wrong) in mutations:
        checks.append(check('independent_exact_certificate:' + name, True, rejects(verify_certificate, wrong, field, roof)))
    result.update(checks=checks, runtime_seconds=time.perf_counter() - start)
    dump(out / 'EXACT_CERTIFICATE_CHECK.json', result)
    return result
