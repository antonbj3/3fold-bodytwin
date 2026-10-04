"""Optional maximum-film requirement, independent of the minimum-film gate."""
import numpy as np
import trimesh
from .exact import q, vec, dot, norm2, sub, strings
from ..geometry import cap, planes_from_cap
from ..io import result, digest

def verify_negative(prep_vertices, point, direction, maximum_mm):
    P = list(map(vec, prep_vertices))
    p = vec(point)
    v = vec(direction)
    g = q(maximum_mm)
    if not P or g < 0 or len(p) != 3 or (len(v) != 3) or any((len(x) != 3 for x in P)) or (not norm2(v)):
        return False
    slack = dot(v, p) - max((dot(v, x) for x in P))
    return slack > 0 and slack * slack > g * g * norm2(v)

def verify_positive(prep_planes, vertices, points_in_prep, maximum_mm):
    g = q(maximum_mm)
    if not prep_planes or not vertices or len(vertices) != len(points_in_prep) or (g < 0):
        return False
    if any((len(a) != 4 for a in prep_planes)) or any((len(p) != 3 for p in vertices)) or any((len(c) != 3 for c in points_in_prep)):
        return False
    if any((not norm2(vec(a[:3])) for a in prep_planes)):
        return False
    for (p, c) in zip(map(vec, vertices), map(vec, points_in_prep)):
        if any((dot(vec(a[:3]), c) > q(a[3]) for a in prep_planes)):
            return False
        if norm2(sub(p, c)) > g * g:
            return False
    return True

def check(preparation, intaglio, maximum_mm):
    (P, F) = cap(preparation)
    (V, _) = cap(intaglio)
    try:
        planes = planes_from_cap(preparation) + [[0, 0, -1, 0]]
    except ValueError:
        return result('UNKNOWN', 'maximum cement-space scenario', reason='Convex prep required')
    (closest, distance, _) = trimesh.proximity.closest_point(trimesh.Trimesh(P, F, process=False), V)
    order = np.argsort(-distance, kind='stable')
    for j in order:
        direction = [f'{x:.12f}' for x in V[j] - closest[j]]
        if verify_negative(P.tolist(), V[j].tolist(), direction, maximum_mm):
            return result('FAIL', 'maximum distance to convex preparation over required intaglio surface', witness=dict(vertex_index=int(j), point_mm=V[j].tolist(), direction=direction), discovery_gap_mm=float(distance[j]), maximum_mm=maximum_mm, requirement_status='MODELED SCENARIO, not a clinical threshold', input_sha256=digest([preparation, intaglio, maximum_mm]))
    centre = np.array([0, 0, preparation['apex_mm'] / 2])
    witnesses = [[f'{x:.12f}' for x in (1 - 1e-08) * c + 1e-08 * centre] for c in closest]
    ok = verify_positive(planes, V.tolist(), witnesses, maximum_mm)
    return result('PASS' if ok else 'UNKNOWN', 'maximum distance to convex preparation over required intaglio surface', certificate=dict(points_in_prep=witnesses) if ok else None, maximum_mm=maximum_mm, requirement_status='MODELED SCENARIO, not a clinical threshold', input_sha256=digest([preparation, intaglio, maximum_mm]))
