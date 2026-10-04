"""K/O set-monotonicity witness, independent of a generator's shape family."""
from .exact import vec, q, dot, norm2, sub, strings
from ..io import result

def point_conflict(required_point, outside_point, envelope_planes, wall_mm):
    (p, o) = (vec(required_point), vec(outside_point))
    w = q(wall_mm)
    if w < 0 or not envelope_planes:
        return result('INVALID', 'K/O necessary set inclusion')
    exterior = [i for (i, pl) in enumerate(envelope_planes) if dot(vec(pl[:3]), o) > q(pl[3])]
    if exterior and norm2(sub(o, p)) <= w * w:
        return result('FAIL', 'Q+B_w not subset T; no preparation K can meet both set requirements', witness=strings(dict(q=p, o=o, exterior_plane=exterior[0])), parents=['LANE_NEXT_K_MIN_PREP_PROOF/R11', 'LANE_NEXT_O_MATERIAL_FEASIBILITY'])
    return result('UNKNOWN', 'K/O necessary set inclusion', reason='This point pair does not prove impossibility')
