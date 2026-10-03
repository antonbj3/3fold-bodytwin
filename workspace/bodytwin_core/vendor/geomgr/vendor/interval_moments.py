"""Outward binary64 intervals for the H5 tetrahedral moment formulas.

Inputs are finite binary64 coordinates/densities; their exact real values are
the values certified. No input discretisation or anatomical error is covered.
"""
from __future__ import annotations

import math


def down(x, n=1):
    for _ in range(n):
        x = math.nextafter(x, -math.inf)
    return x


def up(x, n=1):
    for _ in range(n):
        x = math.nextafter(x, math.inf)
    return x


class IV:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = float(lo)
        self.hi = self.lo if hi is None else float(hi)
        if not (math.isfinite(self.lo) and math.isfinite(self.hi) and self.lo <= self.hi):
            raise ValueError('invalid or nonfinite interval')

    @staticmethod
    def of(x):
        return x if isinstance(x, IV) else IV(x)

    def __add__(self, other):
        q = IV.of(other)
        return IV(down(self.lo + q.lo), up(self.hi + q.hi))

    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -IV.of(other)

    def __rsub__(self, other):
        return IV.of(other) - self

    def __mul__(self, other):
        q = IV.of(other)
        v = [a * b for a in (self.lo, self.hi) for b in (q.lo, q.hi)]
        return IV(down(min(v)), up(max(v)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        q = IV.of(other)
        if q.lo <= 0 <= q.hi:
            raise ZeroDivisionError('interval denominator includes zero')
        v = [a / b for a in (self.lo, self.hi) for b in (q.lo, q.hi)]
        return IV(down(min(v)), up(max(v)))

    def __rtruediv__(self, other):
        return IV.of(other) / self

    def pair(self):
        return [self.lo, self.hi]

    def contains(self, x):
        from fractions import Fraction
        x = Fraction(x)
        return Fraction(self.lo) <= x <= Fraction(self.hi)

    def width(self):
        return up(self.hi - self.lo)


def isum(xs):
    xs = list(xs)
    if not xs:
        return IV(0)
    # CPython fsum is near/correctly rounded on IEEE hardware. Two ulps
    # cover the documented possible one-ulp extended-precision double rounding.
    return IV(down(math.fsum(x.lo for x in xs), 2),
              up(math.fsum(x.hi for x in xs), 2))


def _vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def tetra_terms(a, b, c):
    """Return 13 local contributions: V, S[3], symmetric T[9]."""
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    six = ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
    u = _vadd(_vadd(a, b), c)
    v = six / 6
    s = [six * x / 24 for x in u]
    t = []
    for i in range(3):
        for j in range(3):
            q = a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + u[i] * u[j]
            t.append(six * q / 120)
    return [v] + s + t


def shell(vertices, faces, ref=None):
    """Certified V,S,T about an exact integer reference in mm."""
    if ref is None:
        ref = [int(round((min(row[i] for row in vertices) +
                          max(row[i] for row in vertices)) / 2)) for i in range(3)]
    if len(ref) != 3 or any(int(x) != x for x in ref):
        raise ValueError('reference must be integral')
    vv = [[IV(float(x)) - IV(int(r)) for x, r in zip(row, ref)] for row in vertices]
    terms = [[], [], [], [], [], [], [], [], [], [], [], [], []]
    for face in faces:
        if len(face) != 3:
            raise ValueError('nontriangle face')
        tri = tetra_terms(*(vv[int(k)] for k in face))
        for dst, x in zip(terms, tri):
            dst.append(x)
    return [isum(col) for col in terms], list(ref)


def mass_from_shells(shells, densities):
    """Nested shells, same reference; density[k] fills shell[k] less shell[k+1]."""
    if len(shells) != len(densities) or not shells:
        raise ValueError('shell/density mismatch')
    out = []
    for j in range(13):
        pieces = []
        for k, rho in enumerate(densities):
            if not math.isfinite(float(rho)) or rho <= 0:
                raise ValueError('invalid density')
            sh = shells[k][j] - shells[k+1][j] if k+1 < len(shells) else shells[k][j]
            pieces.append(sh * IV(float(rho)))
        out.append(isum(pieces))
    if out[0].lo <= 0:
        raise ValueError('nonpositive/uncertain mass')
    return out


def derive(moment, ref=(0, 0, 0)):
    """Mass, global COM, central inertia. Tensors are row major."""
    m, s, t = moment[0], moment[1:4], moment[4:]
    c_local = [x / m for x in s]
    com = [x + IV(int(r)) for x, r in zip(c_local, ref)]
    tc = [t[3*i+j] - s[i]*s[j]/m for i in range(3) for j in range(3)]
    trace = isum([tc[0], tc[4], tc[8]])
    inertia = [trace - tc[3*i+j] if i == j else -tc[3*i+j]
               for i in range(3) for j in range(3)]
    return {'m': m, 'com': com, 'I': inertia}


def affine(moment, F, trans):
    """Affine transport x'=Fx+t, including determinant and volume scaling."""
    f = [[IV.of(x) for x in row] for row in F]
    t = [IV.of(x) for x in trans]
    det = (f[0][0]*(f[1][1]*f[2][2]-f[1][2]*f[2][1])
           - f[0][1]*(f[1][0]*f[2][2]-f[1][2]*f[2][0])
           + f[0][2]*(f[1][0]*f[2][1]-f[1][1]*f[2][0]))
    if det.lo <= 0 <= det.hi:
        raise ValueError('singular or uncertain affine determinant')
    d = det if det.lo > 0 else -det
    m, s, T = moment[0], moment[1:4], moment[4:]
    fs = [isum(f[i][k]*s[k] for k in range(3)) for i in range(3)]
    mm = d*m
    ss = [d*(fs[i]+t[i]*m) for i in range(3)]
    tt = []
    for i in range(3):
        for j in range(3):
            ftf = isum(f[i][k]*T[3*k+l]*f[j][l] for k in range(3) for l in range(3))
            tt.append(d*(ftf+fs[i]*t[j]+t[i]*fs[j]+m*t[i]*t[j]))
    return [mm]+ss+tt
