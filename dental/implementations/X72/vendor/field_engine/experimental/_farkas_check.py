"""Exact producer-independent witness checking. No floating-point inputs admitted.

NEJ means infeasible under the supplied mathematical model, JA means an exact
primal witness, UNKNOWN means no verified witness. Geometry/model binding is
the caller's contract. All inequality variables are unrestricted. SOC blocks
are Lorentz cones (t,u,v); mixed products of LP and SOC are not yet supported.
"""
from fractions import Fraction as Q
from dataclasses import dataclass

def rational(v):
    if isinstance(v, bool) or isinstance(v, float):
        raise ValueError('exact rational input required')
    if isinstance(v, (int, str, Q)):
        return Q(v)
    raise ValueError('exact rational input required')

def vector(v):
    return tuple((rational(q) for q in v))

def dot(a, b):
    if len(a) != len(b):
        raise ValueError('dimension mismatch')
    return sum((x * y for (x, y) in zip(a, b)), Q(0))

def matrix(A, b):
    A = tuple((vector(r) for r in A))
    b = vector(b)
    if not A or not A[0] or len(A) != len(b) or any((len(r) != len(A[0]) for r in A)):
        raise ValueError('dimension mismatch')
    return (A, b)

def transpose_action(A, y):
    return tuple((sum((y[i] * A[i][j] for i in range(len(A))), Q(0)) for j in range(len(A[0]))))

def soc3(v):
    return len(v) == 3 and v[0] >= 0 and (v[0] * v[0] >= v[1] * v[1] + v[2] * v[2])

@dataclass(frozen=True)
class Decision:
    status: str
    reason: str
    separation: Q | None = None
    support: tuple = ()

def check_dual(kind, A, b, y, upper=None):
    try:
        (A, b) = matrix(A, b)
        y = vector(y)
        if len(y) != len(b):
            raise ValueError('dimension mismatch')
        r = transpose_action(A, y)
        c = dot(y, b)
        support = tuple((i for (i, v) in enumerate(y) if v))
        if kind == 'inequality':
            if any((v < 0 for v in y)):
                return Decision('UNKNOWN', 'negative inequality multiplier')
            if any(r):
                return Decision('UNKNOWN', 'unrestricted residual is nonzero')
            if upper is not None:
                u = vector(upper)
                if len(u) != len(b) or any((l > h for (l, h) in zip(b, u))):
                    raise ValueError('invalid box')
                c = dot(y, u)
        elif kind == 'orthant':
            if upper is not None:
                raise ValueError('box not supported for orthant')
            if any((v < 0 for v in r)):
                return Decision('UNKNOWN', 'transpose action not in dual orthant')
        elif kind == 'soc3':
            if upper is not None:
                raise ValueError('box not supported for SOC')
            if not soc3(r):
                return Decision('UNKNOWN', 'transpose action not in dual SOC')
        else:
            raise ValueError('unknown cone kind')
        if c >= 0:
            return Decision('UNKNOWN', 'separation is not strictly negative', c, support)
        return Decision('NEJ', 'exact strict dual certificate', c, support)
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as e:
        return Decision('UNKNOWN', str(e))

def check_primal(kind, A, b, x):
    try:
        (A, b) = matrix(A, b)
        x = vector(x)
        if len(x) != len(A[0]):
            raise ValueError('dimension mismatch')
        if kind == 'inequality':
            valid = all((dot(row, x) <= bi for (row, bi) in zip(A, b)))
        elif kind == 'orthant':
            valid = all((q >= 0 for q in x)) and all((dot(row, x) == bi for (row, bi) in zip(A, b)))
        elif kind == 'soc3':
            valid = soc3(x) and all((dot(row, x) == bi for (row, bi) in zip(A, b)))
        else:
            raise ValueError('unknown cone kind')
        return Decision('JA' if valid else 'UNKNOWN', 'exact primal certificate' if valid else 'invalid primal')
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as e:
        return Decision('UNKNOWN', str(e))

def decide_system(kind, A, b, *, dual=None, primal=None, upper=None):
    if dual is not None:
        d = check_dual(kind, A, b, dual, upper)
        if d.status == 'NEJ':
            return d
    if primal is not None and upper is None:
        d = check_primal(kind, A, b, primal)
        if d.status == 'JA':
            return d
    return Decision('UNKNOWN', 'no verified witness')

def support_reusable(old_A, old_b, new_A, new_b, y):
    """Full input binding is O(input size); only arithmetic is support-local."""
    try:
        (old_A, old_b) = matrix(old_A, old_b)
        (new_A, new_b) = matrix(new_A, new_b)
        y = vector(y)
        if len(y) != len(old_b) or len(old_A) != len(new_A) or len(old_A[0]) != len(new_A[0]):
            return False
        support = [i for (i, v) in enumerate(y) if v]
        if not support:
            return False
        return all((old_A[i] == new_A[i] and old_b[i] == new_b[i] for i in support))
    except (ValueError, TypeError, ZeroDivisionError):
        return False
