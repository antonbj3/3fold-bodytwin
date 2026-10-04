"""Input-bound certificates for a finite rational normal cone, in any dimension.

The 2d charts cover all nonzero directions by max-norm rescaling. Existing Farkas
arithmetic is reused. A normal-sign certificate is local geometry, not a swept
path certificate. `relaxation` bounds each row's action error at max-norm one;
establishing that bound from measurement/normalization belongs to the caller.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from ._farkas_check import rational, check_dual, check_primal

def _normals(normals):
    rows = tuple((tuple((rational(v) for v in r)) for r in normals))
    if not rows or not rows[0] or any((len(r) != len(rows[0]) for r in rows)):
        raise ValueError('nonempty rectangular rational normals required')
    if len(rows[0]) > 16:
        raise ValueError('dimension budget exceeded')
    return rows

def cone_chart(normals, axis, sign, *, relaxation=0):
    rows = _normals(normals)
    n = len(rows[0])
    relax = rational(relaxation)
    if relax < 0 or type(axis) is not int or (not 0 <= axis < n) or (type(sign) is not int) or (sign not in (-1, 1)):
        raise ValueError('invalid chart or relaxation')
    (A, b) = ([tuple((-v for v in r)) for r in rows], [relax] * len(rows))
    for j in range(n):
        for s in (-1, 1):
            A.append(tuple((Q(s if k == j else 0) for k in range(n))))
            b.append(Q(1))
    A.append(tuple((Q(-sign if k == axis else 0) for k in range(n))))
    b.append(Q(-1))
    return (tuple(A), tuple(b))

@dataclass(frozen=True)
class ConeDecision:
    status: str
    reason: str
    checked_charts: int = 0

def verify_no_direction(normals, certificates, *, relaxation=0):
    """Each entry carries axis, sign, sparse original row indices and y weights.

    Supplied matrices/labels cannot replace the chart reconstructed from input.
    Missing, duplicate or invalid charts return UNKNOWN, including empty lists.
    """
    try:
        rows = _normals(normals)
        n = len(rows[0])
        certs = tuple(certificates)
        if len(certs) != 2 * n:
            return ConeDecision('UNKNOWN', 'all 2d charts required')
        seen = set()
        for c in certs:
            (axis, sign) = (c['axis'], c['sign'])
            (A, b) = cone_chart(rows, axis, sign, relaxation=relaxation)
            key = (axis, sign)
            if key in seen:
                return ConeDecision('UNKNOWN', 'duplicate chart')
            (ids, y) = (tuple(c['rows']), tuple(c['y']))
            if len(ids) != len(y) or len(set(ids)) != len(ids) or any((type(i) is not int or not 0 <= i < len(A) for i in ids)):
                return ConeDecision('UNKNOWN', 'invalid sparse row binding')
            if check_dual('inequality', [A[i] for i in ids], [b[i] for i in ids], y).status != 'NEJ':
                return ConeDecision('UNKNOWN', 'unverified chart')
            seen.add(key)
        return ConeDecision('NO_NONZERO_DIRECTION', 'exact infeasibility of covering charts', len(seen))
    except (ValueError, TypeError, KeyError, ZeroDivisionError, OverflowError):
        return ConeDecision('UNKNOWN', 'invalid exact certificate')

def verify_direction(normals, direction, *, relaxation=0):
    try:
        rows = _normals(normals)
        d = tuple((rational(v) for v in direction))
        relax = rational(relaxation)
        if relax < 0 or len(d) != len(rows[0]) or (not any(d)):
            return ConeDecision('UNKNOWN', 'nonzero exact direction required')
        bound = relax * max((abs(v) for v in d))
        A = [tuple((-v for v in r)) for r in rows]
        if check_primal('inequality', A, [-bound] * len(rows), d).status != 'JA':
            return ConeDecision('UNKNOWN', 'invalid direction under row envelope')
        strict = all((sum((a * v for (a, v) in zip(r, d)), Q(0)) > bound for r in rows))
        return ConeDecision('YES_STRICT' if strict else 'YES_BOUNDARY', 'exact nonzero primal direction')
    except (ValueError, TypeError, ZeroDivisionError, OverflowError):
        return ConeDecision('UNKNOWN', 'invalid exact input')
