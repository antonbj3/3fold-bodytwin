import numpy as np
from .exact import q, vec, norm2, sub, aabb_distance2, triangle_distance, strings
from .cement import check as offset_check
from ..geometry import cap, planes_from_cap
from ..io import result

def check(outer, inner, minimum_mm, max_exact_pairs=2500):
    scope = 'minimum Euclidean separation of complete outer/intaglio caps; cervical rim excluded'
    w = q(minimum_mm)
    if w <= 0:
        return result('INVALID', scope)
    (ov, of) = cap(outer)
    (iv, inf) = cap(inner)
    if np.any(np.asarray(inner['radii_mm'][0]) >= np.asarray(outer['radii_mm'][0])) or inner['apex_mm'] >= outer['apex_mm']:
        return result('FAIL', scope, reason='Intaglio protrudes outside outer radial cap at rim or apex')
    try:
        planes = planes_from_cap(outer)
        out = offset_check(iv.tolist(), planes, minimum_mm)
        if out['status'] == 'PASS':
            return result('PASS', scope, method='Exact convex outer supporting planes', certificate=out)
    except ValueError:
        pass
    A = ov[of]
    B = iv[inf]
    low = A.min(1)
    high = A.max(1)
    blo = B.min(1)
    bhi = B.max(1)
    gap = np.maximum(0, np.maximum(blo[None] - high[:, None], low[:, None] - bhi[None]))
    order = np.argsort(np.sum(gap * gap, axis=2).ravel(), kind='stable')
    Ar = [tuple(map(vec, t)) for t in A]
    Br = [tuple(map(vec, t)) for t in B]
    count = 0
    pruned = 0
    for index in order:
        (i, j) = divmod(int(index), len(B))
        (a, b) = (Ar[i], Br[j])
        if aabb_distance2(a, b) >= w * w:
            pruned += 1
            continue
        if count >= max_exact_pairs:
            return result('UNKNOWN', scope, reason='Exact pair budget reached', exact_pairs=count, exact_pruned=pruned)
        (dd, p, pp) = triangle_distance(a, b)
        count += 1
        if dd < w * w:
            return result('FAIL', scope, witness=strings(dict(outer_triangle=i, inner_triangle=j, squared_distance_mm2=dd, outer_point=p, inner_point=pp)), exact_pairs=count)
    return result('PASS', scope, method='Exhaustive exact triangle distances with exact AABB pruning', exact_pairs=count, exact_pruned=pruned)
