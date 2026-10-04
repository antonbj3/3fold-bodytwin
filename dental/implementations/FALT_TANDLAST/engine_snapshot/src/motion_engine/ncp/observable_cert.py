"""Exact certificates for an affine *linearized* contact response A x = h.

Inputs are integers, Fraction or rational strings (floats are rejected). This
module does not establish branch existence, rank stability, or nonlinear NCP
uniqueness. C is the explicitly requested readout, e.g. J.T for body impulse.
An exact certificate applies to the supplied rational model, not rounded real data.
"""
from fractions import Fraction

def _q(x):
    if type(x) not in (int, str, Fraction):
        raise TypeError('use integers, rational strings or Fraction; no floats')
    return Fraction(x)

def _inputs(A, h, C):
    A = [[_q(x) for x in row] for row in A]
    h = [_q(x) for x in h]
    C = [[_q(x) for x in row] for row in C]
    if not A or not A[0] or (not C):
        raise ValueError('A and C must be nonempty matrices')
    n = len(A[0])
    if len(h) != len(A) or any((len(r) != n for r in A + C)):
        raise ValueError('incompatible dimensions')
    return (A, h, C)

def _mv(A, x):
    return [sum((a * b for (a, b) in zip(row, x))) for row in A]

def _mm(A, B):
    return [[sum((a * b for (a, b) in zip(row, col))) for col in zip(*B)] for row in A]

def _rref(A):
    (m, n) = (len(A), len(A[0]))
    R = [row[:] for row in A]
    E = [[Fraction(i == j) for j in range(m)] for i in range(m)]
    pivots = []
    for j in range(n):
        k = len(pivots)
        p = next((i for i in range(k, m) if R[i][j]), None)
        if p is None:
            continue
        (R[k], R[p]) = (R[p], R[k])
        (E[k], E[p]) = (E[p], E[k])
        v = R[k][j]
        R[k] = [x / v for x in R[k]]
        E[k] = [x / v for x in E[k]]
        for i in range(m):
            if i != k:
                v = R[i][j]
                R[i] = [a - v * b for (a, b) in zip(R[i], R[k])]
                E[i] = [a - v * b for (a, b) in zip(E[i], E[k])]
        pivots.append(j)
        if len(pivots) == m:
            break
    return (R, E, pivots)

def certify_observable(A, h, C):
    """Return INCOMPATIBLE, AMBIGUOUS, or UNIQUE with an exact witness.

    UNIQUE includes x0 with A*x0=h and L with L*A=C, hence C*x=L*h
    for every solution. AMBIGUOUS includes feasible x0 and z with A*z=0,
    C*z!=0. INCOMPATIBLE includes y with y*A=0 and y*h!=0.
    Discovery is standard rational row elimination; no speed advantage claimed.
    """
    (A, h, C) = _inputs(A, h, C)
    (m, n) = (len(A), len(A[0]))
    (R, E, piv) = _rref(A)
    rhs = _mv(E, h)
    for i in range(len(piv), m):
        if rhs[i]:
            return {'status': 'INCOMPATIBLE', 'left_null': E[i]}
    x = [Fraction(0) for _ in range(n)]
    for (i, j) in enumerate(piv):
        x[j] = rhs[i]
    for j in range(n):
        if j in piv:
            continue
        z = [Fraction(0) for _ in range(n)]
        z[j] = Fraction(1)
        for (i, p) in enumerate(piv):
            z[p] = -R[i][j]
        if any(_mv(C, z)):
            return {'status': 'AMBIGUOUS', 'particular': x, 'right_null': z}
    L = [[sum((row[p] * E[i][k] for (i, p) in enumerate(piv))) for k in range(m)] for row in C]
    return {'status': 'UNIQUE', 'particular': x, 'readout_map': L, 'value': _mv(C, x)}

def verify_observable(A, h, C, certificate):
    """Check certificate identities without elimination; invalid witness -> False."""
    (A, h, C) = _inputs(A, h, C)
    (m, n) = (len(A), len(A[0]))
    try:

        def vec(key, size):
            x = [_q(v) for v in certificate[key]]
            if len(x) != size:
                raise ValueError('certificate dimension')
            return x
        status = certificate['status']
        if status == 'INCOMPATIBLE':
            y = vec('left_null', m)
            return not any(_mv(list(zip(*A)), y)) and bool(sum((a * b for (a, b) in zip(y, h))))
        x = vec('particular', n)
        if _mv(A, x) != h:
            return False
        if status == 'AMBIGUOUS':
            z = vec('right_null', n)
            return not any(_mv(A, z)) and any(_mv(C, z))
        if status == 'UNIQUE':
            L = [[_q(v) for v in row] for row in certificate['readout_map']]
            value = vec('value', len(C))
            return len(L) == len(C) and all((len(row) == m for row in L)) and (_mm(L, A) == C) and (value == _mv(C, x))
        return False
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False
