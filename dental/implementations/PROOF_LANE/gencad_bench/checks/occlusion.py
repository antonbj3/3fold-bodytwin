"""Exact affine gap on projected triangles. Arbitrary 3D scenes remain UNKNOWN."""
from .exact import q, vec, sub, cross, dot, solve, strings
from ..io import result

def _plane(t):
    return solve([[p[0], p[1], 1] for p in t], [p[2] for p in t])

def _inside(p, t):
    side = []
    for i in range(3):
        (a, b) = (t[i], t[(i + 1) % 3])
        side.append((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]))
    return all((x >= 0 for x in side)) or all((x <= 0 for x in side))

def pair(U, L, minimum_mm=0):
    (U, L) = (list(map(vec, U)), list(map(vec, L)))
    (pu, pl) = (_plane(U), _plane(L))
    scope = 'one projected affine upper/lower triangle pair'
    if pu is None or pl is None:
        return result('UNKNOWN', scope, reason='Degenerate projection')
    pts = [p[:2] for p in U if _inside(p, L)] + [p[:2] for p in L if _inside(p, U)]
    for i in range(3):
        (a, b) = (U[i], U[(i + 1) % 3])
        for j in range(3):
            (c, d) = (L[j], L[(j + 1) % 3])
            s = solve([[b[k] - a[k], c[k] - d[k]] for k in range(2)], [c[k] - a[k] for k in range(2)])
            if s is not None and all((0 <= v <= 1 for v in s)):
                pts.append(tuple((a[k] + s[0] * (b[k] - a[k]) for k in range(2))))
    if not pts:
        return result('PASS', scope, reason='Exact empty projected overlap', minimum_gap_mm=None)
    coeff = sub(pu, pl)
    values = [dot(coeff, tuple(p) + (q(1),)) for p in pts]
    i = min(range(len(values)), key=values.__getitem__)
    return result('PASS' if values[i] >= q(minimum_mm) else 'FAIL', scope, minimum_gap_mm=str(values[i]), witness_xy_mm=strings(pts[i]), arithmetic='exact rational')

def check_scene(upper, lower, *, single_valued_height_scope=False, minimum_mm=0):
    if not single_valued_height_scope:
        return result('UNKNOWN', 'full 3D nonpenetration', reason='Projected gap cannot certify arbitrary solids')
    if not upper or not lower:
        return result('UNKNOWN', 'projected gap', reason='Missing surface')
    rows = [pair(a, b, minimum_mm) for a in upper for b in lower]
    failures = [r for r in rows if r['status'] == 'FAIL']
    if failures:
        return result('FAIL', 'all projected pairs', witness=failures[0])
    if any((r['status'] == 'UNKNOWN' for r in rows)):
        return result('UNKNOWN', 'all projected pairs', reason='Degenerate projection')
    return result('PASS', 'all projected pairs', pairs=len(rows), physical_registration_error='UNKNOWN')
