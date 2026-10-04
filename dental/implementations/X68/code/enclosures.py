"""Exact arithmetic bounds for a fitted quadratic, not for physical torque.

All supplied binary floats are interpreted as exact rational constants. Derivative
extrema are certified here only if the second derivative has one sign throughout
the positive bore interval. The fallback uses the convex hull of exact Bernstein
coefficients and is conservative. Outward float rounding encloses exact endpoints.
"""
from fractions import Fraction as F
import math

def outward(q, lower):
    v = float(q)
    if lower and F(v) > q or (not lower and F(v) < q):
        v = math.nextafter(v, -math.inf if lower else math.inf)
    return v

def certified_derivative(poly, diameter, bore_lo, bore_hi):
    (c, b, _) = [F(float(x)) for x in poly]
    (D, lo, hi) = map(lambda x: F(float(x)), [diameter, bore_lo, bore_hi])
    if not 0 < lo < hi:
        raise ValueError('Require 0 < bore_lo < bore_hi')
    A = b + 2 * c * D * D
    g = lambda d: -2 * A * d + 4 * c * d * d * d
    gp_ends = [-2 * A + 12 * c * d * d for d in [lo, hi]]
    if max(gp_ends) <= 0:
        (lower, upper) = (g(hi), g(lo))
        proof = 'T double-prime <= 0 on positive interval; exact endpoint extrema'
    elif min(gp_ends) >= 0:
        (lower, upper) = (g(lo), g(hi))
        proof = 'T double-prime >= 0 on positive interval; exact endpoint extrema'
    else:
        h = hi - lo
        a0 = g(lo)
        a1 = h * (-2 * A + 12 * c * lo * lo)
        a2 = 12 * c * lo * h * h
        a3 = 4 * c * h * h * h
        bernstein = [a0, a0 + a1 / 3, a0 + 2 * a1 / 3 + a2 / 3, a0 + a1 + a2 + a3]
        (lower, upper) = (min(bernstein), max(bernstein))
        proof = 'Exact cubic Bernstein convex hull on unit interval'
    return {'dT_dDf_Ncm_per_mm': [outward(lower, True), outward(upper, False)], 'exact_rational_bounds': [str(lower), str(upper)], 'proof': proof, 'certified_nonincreasing_fitted_torque': upper <= 0, 'scope': 'All bore values in interval; exact supplied binary-float polynomial only', 'physical_model_remainder_bound': 'MISSING', 'resolution': 'PHENOMENOLOGICAL'}
