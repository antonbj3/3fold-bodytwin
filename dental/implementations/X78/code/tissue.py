"""Discrete tissue distances with exact half-voxel squared distances.

This is annotation geometry, not an anatomical error certificate. Inputs use
native array axes z,y,x; tick = 3/20 mm in the validation data.
"""
import math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import lsq_linear
from scipy.spatial import cKDTree

def sqrt_bounds(squared):
    squared = F(squared)
    if squared == 0:
        return [0.0, 0.0]
    x = math.sqrt(float(squared))
    lo = x
    hi = x
    while F.from_float(lo) ** 2 > squared:
        lo = math.nextafter(lo, -math.inf)
    while F.from_float(hi) ** 2 < squared:
        hi = math.nextafter(hi, math.inf)
    return [lo, hi]

def exact_box_query(query_ticks, centres_indices, tick=F(3, 20)):
    q = np.asarray(query_ticks, dtype=np.int64)
    c = np.asarray(centres_indices, dtype=np.int64) * 2
    if q.shape != (3,) or not len(c):
        raise ValueError('Invalid or empty point geometry')
    assert max(abs(q).max(), abs(c).max()) < 10000000
    delta = np.maximum(np.abs(c - q) - 1, 0)
    s = np.sum(delta * delta, axis=1)
    j = int(s.argmin())
    exact = F(int(s[j])) * tick * tick
    w = np.clip(q, c[j] - 1, c[j] + 1)
    return dict(squared_mm2=str(exact), distance_interval_mm=sqrt_bounds(exact), witness_voxel_index=j, witness_ticks=w.tolist())

def scipy_box_control(query_ticks, centres_indices, tick=F(3, 20)):
    q = np.asarray(query_ticks, float) * float(tick)
    c = np.asarray(centres_indices, float) * float(2 * tick)
    half = float(tick)
    tree = cKDTree(c)
    d0 = tree.query(q, workers=1)[0]
    candidates = tree.query_ball_point(q, math.nextafter(d0 + math.sqrt(3) * half + 1e-10, math.inf))
    distances = []
    for j in candidates:
        fit = lsq_linear(np.eye(3), q, bounds=(c[j] - half, c[j] + half), method='bvls', tol=1e-14, max_iter=200)
        if not fit.success:
            raise ValueError('Independent projection failed')
        distances.append(float(np.linalg.norm(fit.x - q)))
    return (min(distances), len(candidates))

def translated_interval(base, shift_norm, tissue_error=None, query_error=None):
    """Exact rational arithmetic plus outward sqrt, no linear sensitivity.

    Missing independently measured physical errors produce UNKNOWN. For
    discrete model sensitivity explicitly pass zero; do not call it calibrated.
    """
    if tissue_error is None or query_error is None:
        return None
    budget = F(shift_norm) + F(tissue_error) + F(query_error)
    lo = max(F(0), F.from_float(base[0]) - budget)
    hi = F.from_float(base[1]) + budget
    a = float(lo)
    b = float(hi)
    while F.from_float(a) > lo:
        a = math.nextafter(a, -math.inf)
    while F.from_float(b) < hi:
        b = math.nextafter(b, math.inf)
    return [a, b]

def rigid_displacement_enclosure(q_ticks, pivot_ticks, translation_bound, angle_bound, tick=F(3, 20)):
    delta = np.asarray(q_ticks, dtype=np.int64) - np.asarray(pivot_ticks, dtype=np.int64)
    r2 = F(int(delta @ delta)) * tick * tick
    ru = F.from_float(sqrt_bounds(r2)[1])
    budget = F(translation_bound) + ru * min(F(angle_bound), F(2))
    b = float(budget)
    while F.from_float(b) < budget:
        b = math.nextafter(b, math.inf)
    return b
