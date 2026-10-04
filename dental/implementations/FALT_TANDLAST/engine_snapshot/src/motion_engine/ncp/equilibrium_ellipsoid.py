"""Exact reduced linear-equilibrium residual ellipsoid.

For symmetric SPD A, Ax*=b, and verified 0<M<=A, the error e=x*-x
obeys e'Me<=-r'e where r=Ax-b. Thus it belongs to the ellipsoid centered
at -M^-1 r/2 with squared radius r'M^-1 r/4 in metric M. Arbitrary linear
readouts are bounded by exact squared support and dyadic sqrt enclosures.
No eigentolerance, approximate inversion or nonlinear joint closure is trusted.
Rigid modes must have been eliminated explicitly; singular matrices are rejected.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt

def _q(v):
    if isinstance(v, bool) or not isinstance(v, (int, str, Q)):
        raise ValueError('exact rational input required; encode floats explicitly')
    return Q(v)

def _matrix(a):
    a = tuple((tuple((_q(v) for v in row)) for row in a))
    n = len(a)
    if not 1 <= n <= 32 or any((len(row) != n for row in a)):
        raise ValueError('square matrix, dimension 1..32 required')
    if any((a[i][j] != a[j][i] for i in range(n) for j in range(n))):
        raise ValueError('exact symmetry required')
    return a

def _psd(a, *, strict=False):
    """Exact Schur elimination; a zero PSD pivot must have a zero row."""
    a = [list(row) for row in a]
    for i in range(len(a)):
        d = a[i][i]
        if d < 0 or (strict and d == 0):
            return False
        if d == 0:
            if any((a[i][j] for j in range(i + 1, len(a)))):
                return False
            continue
        for j in range(i + 1, len(a)):
            for k in range(j, len(a)):
                a[j][k] -= a[j][i] * a[i][k] / d
                a[k][j] = a[j][k]
    return True

def _solve(a, b):
    a = [list(row) + [rhs] for (row, rhs) in zip(a, b)]
    n = len(a)
    for i in range(n):
        pivot = a[i][i]
        if not pivot:
            raise ValueError('singular system')
        a[i] = [v / pivot for v in a[i]]
        for j in range(i + 1, n):
            q = a[j][i]
            a[j] = [v - q * w for (v, w) in zip(a[j], a[i])]
    x = [Q(0)] * n
    for i in reversed(range(n)):
        x[i] = a[i][-1] - sum((a[i][j] * x[j] for j in range(i + 1, n)), Q(0))
    return tuple(x)

def _dot(a, b):
    return sum((x * y for (x, y) in zip(a, b)), Q(0))

def _sqrt_upper(q, bits):
    if q < 0 or type(bits) is not int or (not 1 <= bits <= 256):
        raise ValueError('invalid squared support or precision')
    den = 1 << bits
    k = isqrt(q.numerator * den * den // q.denominator)
    return Q(k if Q(k, den) ** 2 == q else k + 1, den)

@dataclass(frozen=True)
class EquilibriumBound:
    lower: Q
    upper: Q
    centre: Q
    support_squared: Q
    unit: str
    status: str = 'CERTIFIED_LINEAR_MODEL'

def equilibrium_readout(A, b, approximate, lower_matrix, readout, *, unit, bits=80):
    """Bound readout*x*, not readout error; model entries must be rational."""
    (A, M) = (_matrix(A), _matrix(lower_matrix))
    n = len(A)
    (b, x, ell) = (tuple((_q(v) for v in seq)) for seq in (b, approximate, readout))
    if len(M) != n or any((len(v) != n for v in (b, x, ell))) or (not isinstance(unit, str)) or (not unit):
        raise ValueError('dimension or unit mismatch')
    difference = [[A[i][j] - M[i][j] for j in range(n)] for i in range(n)]
    if not _psd(M, strict=True) or not _psd(difference):
        raise ValueError('exact 0<M<=A required')
    r = tuple((_dot(row, x) - rhs for (row, rhs) in zip(A, b)))
    (v, z) = (_solve(M, r), _solve(M, ell))
    centre = _dot(ell, x) - _dot(ell, v) / 2
    support_squared = _dot(r, v) * _dot(ell, z) / 4
    radius = _sqrt_upper(support_squared, bits)
    return EquilibriumBound(centre - radius, centre + radius, centre, support_squared, unit)
