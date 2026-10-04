"""Small exact rational kernel. Serialized decimals, not binary float residuals.

Rational square-root brackets are outward, with inequalities checked in integers.
These functions prove the encoded geometry, not a scanner's physical accuracy.
"""
from fractions import Fraction as F
from math import isqrt
from itertools import combinations

def q(x):
    if isinstance(x, F):
        return x
    return F(str(x))

def vec(x):
    return tuple((q(v) for v in x))

def dot(a, b):
    return sum((x * y for (x, y) in zip(a, b)), F(0))

def sub(a, b):
    return tuple((x - y for (x, y) in zip(a, b)))

def add(a, b):
    return tuple((x + y for (x, y) in zip(a, b)))

def mul(s, a):
    return tuple((s * x for x in a))

def norm2(a):
    return dot(a, a)

def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

def sqrt_bounds(x, digits=12):
    x = q(x)
    if x < 0:
        raise ValueError('Negative square')
    den = 10 ** digits
    low = F(isqrt(x.numerator * den * den // x.denominator), den)
    high = low if low * low == x else low + F(1, den)
    assert low * low <= x <= high * high
    return (low, high)

def solve(A, b):
    n = len(b)
    M = [list(map(q, row)) + [q(bb)] for (row, bb) in zip(A, b)]
    for j in range(n):
        piv = next((i for i in range(j, n) if M[i][j]), None)
        if piv is None:
            return None
        (M[j], M[piv]) = (M[piv], M[j])
        s = M[j][j]
        M[j] = [x / s for x in M[j]]
        for i in range(n):
            if i != j:
                s = M[i][j]
                M[i] = [x - s * y for (x, y) in zip(M[i], M[j])]
    return tuple((row[-1] for row in M))

def point_segment(p, a, b):
    e = sub(b, a)
    nn = norm2(e)
    t = max(F(0), min(F(1), dot(sub(p, a), e) / nn)) if nn else F(0)
    c = add(a, mul(t, e))
    return (norm2(sub(p, c)), c)

def point_triangle(p, t):
    (a, b, c) = t
    u = sub(b, a)
    v = sub(c, a)
    w = sub(p, a)
    uv = dot(u, v)
    sol = solve([[norm2(u), uv], [uv, norm2(v)]], [dot(u, w), dot(v, w)])
    if sol is not None and min(sol) >= 0 and (sum(sol) <= 1):
        pp = add(a, add(mul(sol[0], u), mul(sol[1], v)))
        return (norm2(sub(p, pp)), pp)
    return min((point_segment(p, t[i], t[(i + 1) % 3]) for i in range(3)), key=lambda x: x[0])

def segment_segment(a, b, c, d):
    (u, v, w) = (sub(b, a), sub(d, c), sub(a, c))
    uv = dot(u, v)
    sol = solve([[norm2(u), -uv], [-uv, norm2(v)]], [-dot(u, w), dot(v, w)])
    candidates = []
    if sol is not None and 0 <= sol[0] <= 1 and (0 <= sol[1] <= 1):
        (p, pp) = (add(a, mul(sol[0], u)), add(c, mul(sol[1], v)))
        candidates.append((norm2(sub(p, pp)), p, pp))
    for p in (a, b):
        (dd, pp) = point_segment(p, c, d)
        candidates.append((dd, p, pp))
    for pp in (c, d):
        (dd, p) = point_segment(pp, a, b)
        candidates.append((dd, p, pp))
    return min(candidates, key=lambda x: x[0])

def segment_triangle(a, b, tri):
    (t0, t1, t2) = tri
    (u, v, e) = (sub(t1, t0), sub(t2, t0), sub(b, a))
    sol = solve([[e[k], -u[k], -v[k]] for k in range(3)], sub(t0, a))
    if sol is not None and 0 <= sol[0] <= 1 and (min(sol[1:]) >= 0) and (sol[1] + sol[2] <= 1):
        return add(a, mul(sol[0], e))
    return None

def triangle_distance(t1, t2):
    """Exact minimum squared distance of two closed (possibly degenerate) triangles."""
    (t1, t2) = (tuple(map(vec, t1)), tuple(map(vec, t2)))
    for (a, b) in ((t1, t2), (t2, t1)):
        for i in range(3):
            p = segment_triangle(a[i], a[(i + 1) % 3], b)
            if p is not None:
                return (F(0), p, p)
    best = []
    for p in t1:
        (d, pp) = point_triangle(p, t2)
        best.append((d, p, pp))
    for pp in t2:
        (d, p) = point_triangle(pp, t1)
        best.append((d, p, pp))
    for i in range(3):
        for j in range(3):
            best.append(segment_segment(t1[i], t1[(i + 1) % 3], t2[j], t2[(j + 1) % 3]))
    return min(best, key=lambda x: x[0])

def aabb_distance2(a, b):
    out = F(0)
    for k in range(3):
        (lo1, hi1) = (min((p[k] for p in a)), max((p[k] for p in a)))
        (lo2, hi2) = (min((p[k] for p in b)), max((p[k] for p in b)))
        gap = max(F(0), lo2 - hi1, lo1 - hi2)
        out += gap * gap
    return out

def strings(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, (list, tuple)):
        return [strings(v) for v in x]
    if isinstance(x, dict):
        return {k: strings(v) for (k, v) in x.items()}
    return x
