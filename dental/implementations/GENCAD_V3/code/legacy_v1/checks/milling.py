"""Exact certificates for a spherical tip in a declared convex access domain.

Positive: vertex centres + recession ray -> entire triangulated patch coverage.
Negative: nonnegative halfspace combination -> a point beyond r+allowance.
Float QP solves only discover candidates. See DESIGN.md for proof and limits.
"""
import numpy as np
from scipy.optimize import minimize, nnls
from .exact import q, vec, dot, norm2, sqrt_bounds, strings
from ..io import result, digest

def _inputs(planes, radius, allowance):
    (r, t) = (q(radius), q(allowance))
    if r <= 0 or t < 0 or (not planes) or any((len(p) != 4 for p in planes)):
        raise ValueError('invalid milling domain')
    A = [vec(p[:3]) for p in planes]
    b = [q(p[3]) for p in planes]
    if any((norm2(a) == 0 for a in A)):
        raise ValueError('zero plane normal')
    norms = [sqrt_bounds(norm2(a)) for a in A]
    return (r, t, A, b, norms)

def verify_positive(planes, vertices, radius, allowance, direction, centres):
    try:
        (r, t, A, b, ns) = _inputs(planes, radius, allowance)
        d = vec(direction)
        if not vertices or len(vertices) != len(centres) or len(d) != 3 or (norm2(d) == 0):
            return False
        if any((len(p) != 3 for p in vertices)) or any((len(c) != 3 for c in centres)):
            return False
        if any((dot(a, d) > 0 for a in A)):
            return False
        for (p, c) in zip(map(vec, vertices), map(vec, centres)):
            if norm2(tuple((x - y for (x, y) in zip(p, c)))) > (r + t) ** 2:
                return False
            if any((dot(a, c) > bb - r * n[1] for (a, bb, n) in zip(A, b, ns))):
                return False
        return True
    except (ValueError, TypeError, ZeroDivisionError):
        return False

def verify_negative(planes, point, radius, allowance, weights):
    try:
        (r, t, A, b, ns) = _inputs(planes, radius, allowance)
        p = vec(point)
        y = vec(weights)
        if len(p) != 3 or len(y) != len(A) or min(y) < 0:
            return False
        v = tuple((sum((yy * a[k] for (yy, a) in zip(y, A))) for k in range(3)))
        beta = sum((yy * (bb - r * n[0]) for (yy, bb, n) in zip(y, b, ns)))
        slack = dot(v, p) - beta
        if norm2(v) == 0:
            return beta < 0
        return slack > 0 and slack * slack > (r + t) ** 2 * norm2(v)
    except (ValueError, TypeError, ZeroDivisionError):
        return False

def check(planes, vertices, radius, allowance, direction=(0, 0, -1)):
    scope = 'ideal ball; convex cavity centre domain; continuous straight exit; residual <= allowance'
    try:
        (r, t, Ar, br, ns) = _inputs(planes, radius, allowance)
    except (ValueError, TypeError):
        return result('INVALID', scope)
    A = np.asarray([[float(v) for v in a] for a in Ar])
    b = np.array([float(v) for v in br])
    lengths = np.linalg.norm(A, axis=1)
    Au = A / lengths[:, None]
    bu = b / lengths - float(r)
    from .inherited import load
    parent = load('gencad_mill_farkas', 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
    bout = [bb - r * n[0] for (bb, n) in zip(br, ns)]
    lp = parent.highs(Ar, bout)
    if lp['ray'] is not None:
        cert = parent.reconstruct(Ar, bout, lp['ray'])
        if cert is not None:
            weights = ['0'] * len(Ar)
            for (index, value) in zip(cert['rows'], cert['y']):
                weights[index] = value
            if vertices and verify_negative(planes, vertices[0], radius, allowance, weights):
                return result('FAIL', scope, witness=dict(point=vertices[0], weights=weights), verified=True, reason='Exact Farkas certificate: even outer erosion is empty', erosion_infeasibility_certificate=cert, qp_calls=0, max_discovery_distance_mm=None, input_sha256=digest([planes, vertices, radius, allowance, list(direction)]))
    centres = []
    max_float_distance = 0.0
    calls = 0
    for p0 in vertices:
        p = np.array(p0, dtype=float)
        opt = minimize(lambda c: 0.5 * np.sum((c - p) ** 2), p - np.array([0, 0, float(r)]), jac=lambda c: c - p, constraints=[dict(type='ineq', fun=lambda c: bu - Au @ c, jac=lambda c: -Au)], method='SLSQP', options={'ftol': 1e-12, 'maxiter': 150})
        calls += 1
        c = opt.x
        max_float_distance = max(max_float_distance, float(np.linalg.norm(p - c)))
        active = np.flatnonzero(np.abs(Au @ c - bu) < 2e-06)
        if len(active):
            (yy, _) = nnls(Au[active].T, p - c, maxiter=max(100, 10 * len(active)))
            y = np.zeros(len(A))
            y[active] = yy / lengths[active]
            weights = [str(q(f'{v:.12f}')) for v in y]
            if verify_negative(planes, p0, radius, allowance, weights):
                return result('FAIL', scope, witness=dict(point=p0, weights=weights), verified=True, qp_calls=calls, max_discovery_distance_mm=max_float_distance, input_sha256=digest([planes, vertices, radius, allowance, list(direction)]))
        inner = np.array([0.0, 0.0, min(-10.0, float(c[2]) - 10.0)])
        if np.max(Au @ inner - bu) > 0:
            inneropt = minimize(lambda x: 0.0, inner, constraints=[dict(type='ineq', fun=lambda x: bu - 1e-05 - Au @ x)], method='SLSQP')
            inner = inneropt.x
        c = (1 - 1e-08) * c + 1e-08 * inner
        centres.append([f'{v:.12f}' for v in c])
    ok = verify_positive(planes, vertices, radius, allowance, direction, centres)
    return result('PASS' if ok else 'UNKNOWN', scope, certificate=dict(centres=centres, direction=list(direction)) if ok else None, reason=None if ok else 'No exact positive or negative certificate found; exit/domain may be unsupported', qp_calls=calls, max_discovery_distance_mm=max_float_distance, input_sha256=digest([planes, vertices, radius, allowance, list(direction)]))

def tangent_obstruction(point, inward_normal, material_point, radius):
    """Smooth-face exact zero-allowance impossibility, any approach direction.
    p and q must be independently bound to retained surface patches by caller.
    """
    (p, n, s) = (vec(point), vec(inward_normal), vec(material_point))
    r = q(radius)
    v = tuple((x - y for (x, y) in zip(p, s)))
    a = norm2(v)
    b = -2 * r * dot(v, n)
    return r > 0 and norm2(n) > 0 and (b > 0) and (a * a * norm2(n) < b * b)
