"""Optimal myelin g-ratio: closed-form cable-geometry optimum, recomputed three ways.

Model: r_i ~ 1/g^2, r_m ~ ln(1/g), c_m ~ 1/ln(1/g) => tau = r_m*c_m is g-INVARIANT => conduction
velocity v ~ g*sqrt(-ln g). The optimum is the closed-form calculus result
g_opt = argmax_g[g*sqrt(-ln g)] = exp(-1/2), cross-checked with sympy, scipy and a 2M-point grid,
plus two variants: an adversary omitting the cable square root (g*(-ln g) -> 1/e = 0.36788) and a
g^2-exponent variant (g^2*sqrt(-ln g) -> 0.7788).

Reads: nothing. Writes: nothing (prints the comparison table).
Gate: each recomputed optimum must match the recorded number to < 1e-4; VERDICT is REPRODUCES
only if all of them do (exit 0), DIFFERS otherwise (exit 1).

Usage: python3 myelin_gratio_rebuild.py --selftest
"""

import sys
import sympy as sp
import numpy as np
from scipy.optimize import minimize_scalar


def sympy_optimum():
    g = sp.symbols('g', positive=True)
    v = g * sp.sqrt(-sp.log(g))
    dv = sp.diff(v, g)
    sol = sp.solve(sp.Eq(dv, 0), g)
    return float(sol[0])


def scipy_optimum(func):
    def neg(gg):
        if gg <= 0 or gg >= 1:
            return 1e9
        return -func(gg)
    res = minimize_scalar(neg, bounds=(1e-9, 1 - 1e-9), method='bounded',
                           options={'xatol': 1e-13})
    return res.x


def grid_optimum(func, n=2_000_000):
    gs = np.linspace(1e-9, 1 - 1e-9, n)
    vs = func(gs)
    return gs[np.argmax(vs)]


def main():
    results = {}

    g_sympy = sympy_optimum()
    g_scipy = scipy_optimum(lambda g: g * np.sqrt(-np.log(g)))
    g_grid = grid_optimum(lambda g: g * np.sqrt(-np.log(g)))
    results['baseline_sympy'] = g_sympy
    results['baseline_scipy'] = g_scipy
    results['baseline_grid_2M'] = g_grid

    g_adv_nosqrt = scipy_optimum(lambda g: g * (-np.log(g)))
    results['adversary_no_sqrt'] = g_adv_nosqrt

    g_g2 = scipy_optimum(lambda g: g**2 * np.sqrt(-np.log(g)))
    results['g_squared_variant'] = g_g2

    recorded = {
        'baseline_sympy': np.exp(-0.5),          # 0.6065306597...
        'baseline_scipy': 0.6065306590,
        'baseline_grid_2M': 0.6065306369,
        'adversary_no_sqrt': 1.0/np.e,            # 0.36787944...
        'g_squared_variant': 0.7788,
    }

    print("=== optimal myelin g-ratio: rebuild vs recorded numbers ===")
    all_ok = True
    for k, v in results.items():
        rec = recorded[k]
        diff = abs(v - rec)
        tol = 1e-4 if 'grid' not in k else 1e-4
        ok = diff < tol
        all_ok &= ok
        print(f"{k:22s} rebuild={v:.10f}  recorded={rec:.10f}  diff={diff:.2e}  {'PASS' if ok else 'FAIL'}")

    print()
    verdict = "REPRODUCES" if all_ok else "DIFFERS"
    print(f"VERDICT: {verdict} -- exact closed-form calculus optimum, no free parameters, no tuning possible.")
    return 0 if all_ok else 1


if __name__ == '__main__':
    if '--selftest' in sys.argv:
        sys.exit(main())
    sys.exit(main())
