"""Exact checker for outward Krawczyk data from a complete true residual.

The caller must enclose R(c) and the TRUE dR/dx over the entire supplied box,
including every free neighbour and pair. Newton/PSD surrogate matrices do not
satisfy this contract. This checker proves only one root inside this box. It
never upgrades a box certificate into a global root count.

All input floats are treated as exact binary rational endpoints; a producer
must round interval endpoints outward before calling this function.
"""
from dataclasses import dataclass
from fractions import Fraction
import math

@dataclass(frozen=True)
class BoxRootCertificate:
    verdict: str
    contraction_upper: Fraction | None
    inclusion_upper: Fraction | None
    residual_id: str
    global_root_count: str = 'UNKNOWN'
    reason: str = ''

def _q(x):
    if isinstance(x, bool):
        raise ValueError('boolean is not an interval endpoint')
    if isinstance(x, float) and (not math.isfinite(x)):
        raise ValueError('nonfinite interval endpoint')
    q = Fraction(x)
    return Fraction(int(q.numerator), int(q.denominator))

def _interval(lo, hi):
    (lo, hi) = (_q(lo), _q(hi))
    if lo > hi:
        raise ValueError('reversed interval')
    return (lo, hi)

def _scale(v, a):
    (lo, hi) = v
    return (a * lo, a * hi) if a >= 0 else (a * hi, a * lo)

def _abs_upper(v):
    return max(abs(v[0]), abs(v[1]))

def certify_root_box(*, residual_id, radius, force_lower, force_upper, jacobian_lower, jacobian_upper, preconditioner, jacobian_role='TRUE_RESIDUAL'):
    """Check strict inclusion and contraction in a weighted infinity norm.

    The centre is implicit in the producer's force enclosure. `radius` must
    enclose the same box used for its Jacobian. The preconditioner is one fixed
    exact point matrix, not an interval-valued matrix or per-iterate surrogate.
    A successful result is conditional on the producer's enclosure contract.
    """

    def unknown(reason):
        return BoxRootCertificate('UNKNOWN', None, None, str(residual_id), reason=reason)
    if not residual_id or jacobian_role != 'TRUE_RESIDUAL':
        return unknown('complete identified true residual required')
    try:
        r = list(map(_q, radius))
        n = len(r)
        if n == 0 or any((v <= 0 for v in r)):
            raise ValueError('strictly positive box radii required')
        if len(force_lower) != n or len(force_upper) != n:
            raise ValueError('force dimensions')
        for a in [jacobian_lower, jacobian_upper, preconditioner]:
            if len(a) != n or any((len(row) != n for row in a)):
                raise ValueError('matrix dimensions')
        f = [_interval(a, b) for (a, b) in zip(force_lower, force_upper)]
        J = [[_interval(jacobian_lower[i][j], jacobian_upper[i][j]) for j in range(n)] for i in range(n)]
        C = [list(map(_q, row)) for row in preconditioner]
        (q, inclusion) = ([], [])
        for i in range(n):
            eta = [_scale(f[k], C[i][k]) for k in range(n)]
            eta = (sum((a for (a, _) in eta)), sum((b for (_, b) in eta)))
            qi = Fraction(0)
            for j in range(n):
                terms = [_scale(J[k][j], C[i][k]) for k in range(n)]
                (lo, hi) = (sum((a for (a, _) in terms)), sum((b for (_, b) in terms)))
                e = (int(i == j) - hi, int(i == j) - lo)
                qi += _abs_upper(e) * r[j] / r[i]
            q.append(qi)
            inclusion.append(_abs_upper(eta) / r[i] + qi)
        (qu, ku) = (max(q), max(inclusion))
        if qu < 1 and ku < 1:
            return BoxRootCertificate('UNIQUE_IN_BOX_IF_ENCLOSURES_VALID', qu, ku, str(residual_id))
        return BoxRootCertificate('UNKNOWN', qu, ku, str(residual_id), reason='strict Krawczyk test failed')
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as exc:
        return unknown(str(exc))
