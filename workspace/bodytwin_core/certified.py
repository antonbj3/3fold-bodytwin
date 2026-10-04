"""Numerical certificate for the K1 p=2 nonnegative equality QP.

The candidate is a one-sided frozen active-set solve. A strict KKT and
primal-dual gap check determines whether the exact fallback is needed.
"""
import numpy as np


def certificate(C, b, a, nu, force_scale=None, tolerance=0.005):
    C = np.asarray(C, float); b = np.asarray(b, float)
    a = np.asarray(a, float); nu = np.asarray(nu, float)
    if not (np.isfinite(C).all() and np.isfinite(b).all() and
            np.isfinite(a).all() and np.isfinite(nu).all()):
        return dict(accepted=False, reason='nonfinite', bound=np.inf)
    s = C.T @ nu
    projected = np.maximum(s, 0.0)
    primal = 0.5 * float(a @ a)
    lower = float(nu @ b - 0.5 * projected @ projected)
    # Subtraction and dot products introduce small signed roundoff in the gap.
    allowance = 64 * np.finfo(float).eps * max(1.0, abs(primal), abs(lower))
    gap = primal - lower
    gap_upper = max(0.0, gap) + allowance
    residual = float(np.linalg.norm(C @ a - b) / max(1.0, np.linalg.norm(b)))
    sign = float(max(np.max(np.maximum(-a, 0), initial=0),
                     np.max(np.abs(a - projected), initial=0)) /
                 max(1.0, np.max(np.abs(a), initial=0)))
    scale = np.ones_like(a) if force_scale is None else np.asarray(force_scale, float)
    f = a * scale
    bound = float(np.max(np.abs(scale)) * np.sqrt(2 * gap_upper) /
                  max(1e-12, np.linalg.norm(f)))
    accepted = bool(residual <= 1e-10 and sign <= 1e-10 and
                    gap >= -allowance and bound <= tolerance)
    return dict(accepted=accepted, reason='accepted' if accepted else 'uncertain',
                primal_residual=residual, dual_residual=sign,
                dual_gap=gap, relative_force_bound=bound, bound=bound)


def frozen_candidate(C, b, mask):
    C = np.asarray(C, float); b = np.asarray(b, float)
    mask = np.asarray(mask, bool)
    if mask.shape != (C.shape[1],) or mask.sum() < C.shape[0]:
        raise ValueError('invalid active-set mask')
    Cs = C[:, mask]
    nu = np.linalg.solve(Cs @ Cs.T, b)
    a = np.zeros(C.shape[1]); a[mask] = Cs.T @ nu
    return a, nu


def solve_step(C, b, exact, previous_mask=None, force_scale=None,
               tolerance=0.005, candidate_hook=None):
    """Return a, dual, receipt. exact(C,b) returns (a,dual)."""
    if previous_mask is not None:
        try:
            a, nu = frozen_candidate(C, b, previous_mask)
            if candidate_hook is not None:
                a, nu = candidate_hook(a, nu)
            cert = certificate(C, b, a, nu, force_scale, tolerance)
            if cert['accepted']:
                return a, nu, dict(cert, fallback=False, mask=a > 0)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            cert = dict(accepted=False, reason='candidate_failure')
    else:
        cert = dict(accepted=False, reason='no_seed')
    a, nu = exact(C, b)
    check = certificate(C, b, a, nu, force_scale, tolerance)
    if check['primal_residual'] > 1e-10 or check['dual_residual'] > 1e-10:
        raise RuntimeError('exact fallback failed KKT')
    return a, nu, dict(check, accepted=False, fallback=True,
                       fallback_reason=cert['reason'], mask=a > 0)


def warm_exact(C, b, seed=None, maxit=40):
    """B176-style dual warm start; cold source solver remains safety fallback."""
    from .emulator import xf4_model
    if seed is None:
        return xf4_model._dual_newton(C, b)
    nu = np.asarray(seed, float).copy()
    if not np.isfinite(nu).all():
        return xf4_model._dual_newton(C, b)
    eye = np.eye(C.shape[0])
    for _ in range(maxit):
        s = C.T @ nu
        a = np.maximum(s, 0)
        grad = C @ a - b
        if np.linalg.norm(grad) <= 1e-12 * max(1.0, np.linalg.norm(b)):
            if certificate(C, b, a, nu)['accepted']:
                return a, nu
            break
        H = C[:, s > 0] @ C[:, s > 0].T
        try:
            d = np.linalg.solve(H + 1e-10 * eye, grad)
        except np.linalg.LinAlgError:
            break
        obj = 0.5 * a @ a - nu @ b
        step = 1.0
        while step > 1e-6:
            nxt = nu - step * d
            pos = np.maximum(C.T @ nxt, 0)
            if 0.5 * pos @ pos - nxt @ b <= obj - 1e-4 * step * grad @ d:
                break
            step *= 0.5
        nu -= step * d
    return xf4_model._dual_newton(C, b)
