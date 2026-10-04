"""Full convex preparation + Euclidean ball containment in a cavity."""
from .exact import q, vec, dot, norm2, strings
from ..io import result, digest

def check(prep_vertices, planes, minimum_mm):
    g = q(minimum_mm)
    if g < 0 or not prep_vertices or (not planes) or any((len(p) != 3 for p in prep_vertices)) or any((len(a) != 4 for a in planes)):
        return result('INVALID', 'convex_offset_containment', reason='empty or negative requirement')
    P = list(map(vec, prep_vertices))
    limiting = None
    for (i, pl) in enumerate(planes):
        (a, b) = (vec(pl[:3]), q(pl[3]))
        nn = norm2(a)
        if not nn:
            return result('INVALID', 'convex_offset_containment', reason='zero normal')
        vals = [dot(a, p) for p in P]
        j = max(range(len(P)), key=vals.__getitem__)
        margin = b - vals[j]
        row = dict(plane=i, vertex=j, signed_numerator=margin, normal_squared=nn)
        if margin < 0 or margin * margin < g * g * nn:
            return result('FAIL', 'Q + Euclidean B_g subset declared convex intaglio', witness=strings(row), input_sha256=digest([prep_vertices, planes, minimum_mm]))
        if limiting is None or margin * margin / nn < limiting[0]:
            limiting = (margin * margin / nn, row)
    return result('PASS', 'Q + Euclidean B_g subset declared convex intaglio', checked_vertices=len(P), checked_planes=len(planes), minimum_gap_squared_mm2=str(limiting[0]), witness=strings(limiting[1]), arithmetic='rational, zero tolerance', input_sha256=digest([prep_vertices, planes, minimum_mm]))
