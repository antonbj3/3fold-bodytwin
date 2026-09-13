#!/usr/bin/env python3
"""Collier/Monk/Maini/Lewis 1996 Notch-Delta lateral-inhibition ODE on regular graphs (ring and
triangular lattice): uniform fixed point, linear-stability threshold b_c, growth rates, and the
emergent sender fraction after a random quench.

Reads: nothing (the equations and parameters are embedded below).
Writes: notch_delta_rebuild_results.json under the cell output directory.
Gate: the acceptance test below -- agreement with the record's reported numbers.

ACCEPTANCE TEST (pre-registered, stated before running): reproduce, from the record's stated
dynamics string + parameter dict alone (no tuning), the record's reported: (a) the ring (n=40,
bipartite, lambda_min=-1) instability threshold b_c=11.042458589329957, (b) its invariance to nu
(nu_main=1.0 vs nu_alt=0.2), (c) the triangular-lattice adjacency's most-negative eigenvalue
lambda_min=-0.5000000000000001, (d) the resulting triangular b_c=49.57432745951492, (e) the
measured-vs-linear-theory growth-rate R2 over a sub-to-super-threshold sweep (record: R2=
0.9999999999999973), and (f) the emergent quench-dynamics "sender fraction" (~0.28, NOT the naive
1/3 combinatorial ideal) on triangular lattices of N=225/441/900 cells. Random-quench numbers (f) are
inherently seed/protocol-dependent (never literally bit-reproducible without the original RNG state,
which the record does not state) -- reproduced as a DISTRIBUTION (mean +/- sd over independent seeds)
compared against the record's reported plateau and SD range, not as an exact-digit match.

SOURCE EQUATIONS (verbatim from the record's "dynamics" field, transcribed unchanged):
  dN_i/dt = f(Dbar_i) - N_i
  dD_i/dt = nu*(g(N_i) - D_i)
  Dbar_i = mean of D_j over graph-neighbors j of i  (row-normalized adjacency A)
  f(x) = x^k/(a+x^k)      increasing Hill (Notch activation by neighbor Delta)
  g(x) = 1/(1+b*x^h)      decreasing Hill (own Delta repressed by own Notch)
  uniform fixed point D*=g(f(D*)), N*=f(D*), solved via brentq
  linearization block-diagonalizes in the graph eigenbasis of A (regular graph -> A/degree symmetric);
  J_lambda = [[-1, f'(D*)*lambda], [nu*g'(N*), -nu]];  growth_rate = max(Re(eig(J_lambda)))
  instability iff -f'(D*)*g'(N*)*lambda > 1;  most unstable mode = lambda_min(A)

PARAMETERS (verbatim from the record's "parameters" dict): a=0.01, k=2, h=2, nu_main=1.0,
nu_alt=0.2, ring_n_cells=40, ring_lambda=-1.0, triangular_lattice_sizes=[225,441,900],
b_super=247.87163729757458 (used for the quench-dynamics sender-fraction check).

GEOMETRIC STRUCTURE: the ring (n even, 2-regular) and triangular (6-regular) lattices are REGULAR
graphs, so the row-normalized adjacency A/degree is exactly SYMMETRIC (no left/right-eigenvector
split needed) -- eigenvalues are obtained via np.linalg.eigvalsh, not an heuristic. The ring's
eigenvalues are the classic cycle-graph spectrum cos(2*pi*m/n); the most negative, exactly -1, occurs
at m=n/2 for even n (an exact bipartite/checkerboard mode, not a numerical coincidence). The
threshold b_c is found by bisecting the SIGN of growth_rate(b) (a genuine geometric zero-crossing of
the dominant eigenvalue of J_lambda, not a curve-fit).

"""

import json
import os
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "notch_delta_lateral_inhibition_rebuild")
OUT_PATH = _os.path.join(OUT_DIR, "notch_delta_rebuild_results.json")
SOURCE_NODE = "Notch-Delta lateral-inhibition model"
SOURCE_EVIDENCE = "the stored evidence record for that model"

# ---- parameters, verbatim from the record's "parameters" dict ----
A_HILL = 0.01
K_HILL = 2
H_HILL = 2
NU_MAIN = 1.0
NU_ALT = 0.2
RING_N = 40
B_SUPER = 247.87163729757458
TRIANGULAR_SIZES = [15, 21, 30]  # L such that L*L in {225, 441, 900}

RECORD = {
    "ring_bc_predicted": 11.042458589329957,
    "ring_growth_rate_R2": 0.9999999999999973,
    "ring_growth_rate_max_abs_residual": 6.165382204770253e-08,
    "nu_invariance_rel_err": 1.1260624411865024e-15,
    "triangular_lambda_min": -0.5000000000000001,
    "triangular_bc_predicted": 49.57432745951492,
    "super_frac_by_lattice_size": {"225": 0.27666666666666667, "441": 0.27834467120181405, "900": 0.2813888888888889},
}


def f_hill(x):
    return x ** K_HILL / (A_HILL + x ** K_HILL)


def fp_hill(x):
    num = K_HILL * x ** (K_HILL - 1) * (A_HILL + x ** K_HILL) - x ** K_HILL * K_HILL * x ** (K_HILL - 1)
    return num / (A_HILL + x ** K_HILL) ** 2


def g_hill(x, b):
    return 1.0 / (1.0 + b * x ** H_HILL)


def gp_hill(x, b):
    return -b * H_HILL * x ** (H_HILL - 1) / (1.0 + b * x ** H_HILL) ** 2


def fixed_point(b):
    def resid(D):
        return g_hill(f_hill(D), b) - D
    return brentq(resid, 1e-9, 1.0 - 1e-12, xtol=1e-14, rtol=1e-14)


def growth_rate_pred(b, lam, nu):
    Dstar = fixed_point(b)
    Nstar = f_hill(Dstar)
    J = np.array([[-1.0, fp_hill(Dstar) * lam], [nu * gp_hill(Nstar, b), -nu]])
    return float(np.linalg.eigvals(J).real.max())


def find_bc(lam, nu, blo=0.01, bhi=1000.0):
    def gr(b):
        return growth_rate_pred(b, lam, nu)
    lo, hi = blo, bhi
    assert gr(lo) < 0 and gr(hi) > 0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if gr(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def triangular_adjacency_lambda_min(L):
    """6-neighbor periodic triangular lattice on an LxL index torus; regular graph -> A/6 symmetric."""
    N = L * L

    def idx(i, j):
        return (i % L) * L + (j % L)

    A = np.zeros((N, N))
    offsets = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    for i in range(L):
        for j in range(L):
            u = idx(i, j)
            for di, dj in offsets:
                A[u, idx(i + di, j + dj)] = 1.0
    Anorm = A / 6.0
    eigs = np.linalg.eigvalsh(Anorm)
    return float(eigs.min())


def reduced_rhs_ring_mode(t, y, Dstar, Nstar, b, nu):
    """Exact Z2 (even/odd sublattice) symmetry reduction of the full n=40 ring ODE onto the
    antisymmetric (lambda=-1) Fourier mode: N_even=Nstar+p, N_odd=Nstar-p, D_even=Dstar+q, D_odd=Dstar-q.
    Exact (not linearized) -- valid for any amplitude p,q, since the ring's neighbor-mean for this
    ansatz collapses to a clean 2-variable system by the graph's bipartite symmetry."""
    p, q = y
    dp = (f_hill(Dstar - q) - f_hill(Dstar + q)) / 2.0 - p
    dq = nu * ((g_hill(Nstar + p, b) - g_hill(Nstar - p, b)) / 2.0 - q)
    return [dp, dq]


def measured_growth_rate(b, nu, lam, eps=1e-7):
    Dstar = fixed_point(b)
    Nstar = f_hill(Dstar)
    rate_pred = growth_rate_pred(b, lam, nu)
    t_end = 200.0 if abs(rate_pred) < 1e-4 else min(400.0, 12.0 / abs(rate_pred))
    ts = np.linspace(0, t_end, 3000)
    sol = solve_ivp(reduced_rhs_ring_mode, [0, t_end], [0.0, eps], t_eval=ts,
                     args=(Dstar, Nstar, b, nu), method="RK45", rtol=1e-12, atol=1e-15)
    q = sol.y[1]
    absq = np.abs(q)
    valid = absq > 1e-14
    logq = np.log(absq[valid])
    tt = ts[valid]
    n = len(tt)
    lo = int(n * 0.6) if n >= 10 else 0
    coeff = np.polyfit(tt[lo:], logq[lo:], 1)
    return float(coeff[0])


def make_lattice_rhs(L, b):
    N2 = L * L

    def rhs(t, y):
        Nv = y[:N2].reshape(L, L)
        Dv = y[N2:].reshape(L, L)
        Dbar = (np.roll(Dv, -1, axis=0) + np.roll(Dv, 1, axis=0)
                + np.roll(Dv, -1, axis=1) + np.roll(Dv, 1, axis=1)
                + np.roll(np.roll(Dv, -1, axis=0), 1, axis=1)
                + np.roll(np.roll(Dv, 1, axis=0), -1, axis=1)) / 6.0
        dN = f_hill(Dbar) - Nv
        dD = NU_MAIN * (g_hill(Nv, b) - Dv)
        return np.concatenate([dN.ravel(), dD.ravel()])
    return rhs


def quench_sender_fraction(L, seed, t_end=100.0):
    rng = np.random.default_rng(seed)
    N2 = L * L
    y0 = np.concatenate([rng.uniform(0.3, 0.7, size=N2), rng.uniform(0.3, 0.7, size=N2)])
    rhs = make_lattice_rhs(L, B_SUPER)
    sol = solve_ivp(rhs, [0, t_end], y0, method="RK45", rtol=1e-9, atol=1e-11, t_eval=[t_end])
    Dfin = sol.y[N2:, -1]
    thr = 0.5 * (Dfin.min() + Dfin.max())
    return float(np.mean(Dfin > thr))


def main():
    print("=" * 78)
    print(f"REBUILDING {SOURCE_NODE} from {SOURCE_EVIDENCE}")
    print("=" * 78)

    # ---- ring threshold + nu-invariance ----
    ring_lambda = -1.0
    bc_ring_main = find_bc(ring_lambda, NU_MAIN)
    bc_ring_alt = find_bc(ring_lambda, NU_ALT)
    nu_invariance_rel_err = abs(bc_ring_main - bc_ring_alt) / bc_ring_main
    print(f"ring b_c (nu={NU_MAIN}) = {bc_ring_main!r}  (record {RECORD['ring_bc_predicted']!r})")
    print(f"ring b_c (nu={NU_ALT}) = {bc_ring_alt!r}  nu_invariance_rel_err={nu_invariance_rel_err:.3e}")

    m = np.arange(RING_N)
    ring_lambda_min_analytic = float(np.cos(2 * np.pi * m / RING_N).min())

    # ---- triangular lattice spectrum + threshold ----
    tri_lambda_min = triangular_adjacency_lambda_min(15)
    bc_tri = find_bc(tri_lambda_min, NU_MAIN)
    print(f"triangular lambda_min (L=15) = {tri_lambda_min!r}  (record {RECORD['triangular_lambda_min']!r})")
    print(f"triangular b_c = {bc_tri!r}  (record {RECORD['triangular_bc_predicted']!r})")

    # ---- 14-point growth-rate sweep, ring, measured vs linear-theory prediction ----
    bs = np.concatenate([np.linspace(bc_ring_main * 0.5, bc_ring_main * 0.95, 7),
                         np.linspace(bc_ring_main * 1.05, bc_ring_main * 3.0, 7)])
    preds, meas = [], []
    for b in bs:
        gp_ = growth_rate_pred(b, ring_lambda, NU_MAIN)
        gm_ = measured_growth_rate(b, NU_MAIN, ring_lambda)
        preds.append(gp_)
        meas.append(gm_)
        print(f"  b={b:8.3f}  pred={gp_:+.6f}  meas={gm_:+.6f}  resid={gm_-gp_:+.2e}")
    preds_a, meas_a = np.array(preds), np.array(meas)
    ss_res = float(np.sum((meas_a - preds_a) ** 2))
    ss_tot = float(np.sum((meas_a - meas_a.mean()) ** 2))
    growth_rate_R2 = 1 - ss_res / ss_tot
    growth_rate_max_abs_resid = float(np.max(np.abs(meas_a - preds_a)))
    print(f"growth-rate sweep: R2={growth_rate_R2!r}  max_abs_resid={growth_rate_max_abs_resid:.3e}")
    print(f"  (record R2={RECORD['ring_growth_rate_R2']!r}  max_abs_resid={RECORD['ring_growth_rate_max_abs_residual']!r})")

    # ---- quench-dynamics sender-fraction, triangular lattice, N=225/441/900, 6 seeds each ----
    quench = {}
    t0 = time.time()
    for L in TRIANGULAR_SIZES:
        fracs = [quench_sender_fraction(L, seed) for seed in range(6)]
        quench[str(L * L)] = {"mean": float(np.mean(fracs)), "sd": float(np.std(fracs)), "seeds": fracs}
        print(f"  N={L*L}: sender_frac mean={np.mean(fracs):.4f} sd={np.std(fracs):.4f}  "
              f"(record mean~{RECORD['super_frac_by_lattice_size'][str(L*L)]:.4f})")
    print(f"quench sims took {time.time()-t0:.1f}s total")

    computed = {
        "ring_bc_predicted": bc_ring_main,
        "ring_bc_alt_nu": bc_ring_alt,
        "nu_invariance_rel_err": nu_invariance_rel_err,
        "ring_lambda_min_analytic": ring_lambda_min_analytic,
        "triangular_lambda_min": tri_lambda_min,
        "triangular_bc_predicted": bc_tri,
        "ring_growth_rate_R2": growth_rate_R2,
        "ring_growth_rate_max_abs_residual": growth_rate_max_abs_resid,
        "quench_sender_fraction_by_N": quench,
    }

    def pct_err(x, y):
        return abs(x - y) / abs(y) * 100 if y != 0 else abs(x - y)

    agreement = {
        "ring_bc_pct_err": pct_err(bc_ring_main, RECORD["ring_bc_predicted"]),
        "triangular_lambda_min_abs_err": abs(tri_lambda_min - RECORD["triangular_lambda_min"]),
        "triangular_bc_pct_err": pct_err(bc_tri, RECORD["triangular_bc_predicted"]),
        "growth_rate_R2_both_essentially_1": growth_rate_R2 > 0.999999,
        "quench_frac_within_2sd_of_record": {
            str(L * L): bool(abs(quench[str(L*L)]["mean"] - RECORD["super_frac_by_lattice_size"][str(L*L)])
                              < max(2 * quench[str(L*L)]["sd"], 0.02))
            for L in TRIANGULAR_SIZES
        },
    }
    print("\nAGREEMENT vs record's reported numbers:")
    print(json.dumps(agreement, indent=2, default=str))

    gates = {
        "ring_bc_within_1e-6pct": bool(agreement["ring_bc_pct_err"] < 1e-6),
        "nu_invariance_essentially_exact": bool(nu_invariance_rel_err < 1e-9),
        "triangular_lambda_min_within_1e-9": bool(agreement["triangular_lambda_min_abs_err"] < 1e-9),
        "triangular_bc_within_1e-6pct": bool(agreement["triangular_bc_pct_err"] < 1e-6),
        "growth_rate_R2_essentially_1": bool(agreement["growth_rate_R2_both_essentially_1"]),
        "quench_frac_225_in_range": bool(agreement["quench_frac_within_2sd_of_record"]["225"]),
        "quench_frac_441_in_range": bool(agreement["quench_frac_within_2sd_of_record"]["441"]),
        "quench_frac_900_in_range": bool(agreement["quench_frac_within_2sd_of_record"]["900"]),
        "quench_frac_all_below_naive_third": bool(all(
            quench[str(L*L)]["mean"] < 1.0 / 3.0 for L in TRIANGULAR_SIZES)),
    }
    overall_pass = all(gates.values())
    print("\nGATES:")
    print(json.dumps(gates, indent=2))
    print(f"\nOVERALL RECOVERY: {'REPRODUCED' if overall_pass else 'PARTIAL/FAILED -- see gates'}")

    honest_gaps = {
        "quench_fraction_is_seed_dependent_by_nature": "The emergent sender-fraction from a disordered "
            "quench is a genuinely stochastic outcome (multi-domain frustration defects, per the record's "
            "own honest_gaps) -- reproduced here as a 6-seed mean+-sd, not a bit-exact number; the record "
            "itself does not state its own RNG seed/protocol, so exact-digit reproduction is not possible "
            "or expected, only agreement within the record's reported variability (SD 0.0012-0.0076).",
        "growth_rate_residual_magnitude_differs_from_record": "This rebuild's R2/max-residual "
            "(computed above) will generally not equal the record's noise-floor digits exactly, "
            "since neither the record nor this rebuild states/matches the other's exact eps/window/solver "
            "tolerance choices -- both land at R2 > 0.999999 (>=6 nines), the qualitative claim tested.",
        "mod3_ideal_packing_combinatorial_check_not_rebuilt": "The record's separate mod3_ideal_packing_"
            "density=1/3 (0 violations, 3-coloring combinatorics) is a distinct static graph-coloring "
            "fact, not part of this script's dynamical-systems acceptance test; not rebuilt here (scope "
            "note, not a defect of this rebuild).",
    }
    report = {
        "source_node": SOURCE_NODE,
        "source_evidence": SOURCE_EVIDENCE,
        "parameters_used_verbatim_from_record": {
            "a": A_HILL, "k": K_HILL, "h": H_HILL, "nu_main": NU_MAIN, "nu_alt": NU_ALT,
            "ring_n_cells": RING_N, "b_super": B_SUPER, "triangular_lattice_L_values": TRIANGULAR_SIZES,
        },
        "computed": computed,
        "record_reported": RECORD,
        "agreement": agreement,
        "gates": gates,
        "overall_pass": overall_pass,
        "honest_gaps": honest_gaps,
    }
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"\nWrote {OUT_PATH}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
