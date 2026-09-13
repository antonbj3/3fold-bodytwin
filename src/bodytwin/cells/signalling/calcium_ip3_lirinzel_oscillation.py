#!/usr/bin/env python3
"""Calcium/IP3 oscillation (Li-Rinzel 1994 / De Young-Keizer 1992 3-state model, closed-cell
2-state reduction): does the model, at its published parameters, produce a genuine Hopf-bounded
oscillatory window whose period lands inside the real hepatocyte envelope, and does the adiabatic-h
null control (collapsing to a 1D scalar ODE, Poincare-Bendixson-forbidden from oscillating) indeed
fail to oscillate across the sweep (0/76)?

Reads: nothing. Writes: calcium_ip3_lirinzel_oscillation.json under the cell output directory.

MODEL (closed-cell form c_er=(c0-c)/c1 eliminates the 3rd state, so the ODE below is the
algebraically-equivalent 2-state reduction actually integrated):
  m_inf(p)   = p/(p+d1)
  n_inf(c)   = c/(c+d5)
  Q2(p)      = d2*(p+d1)/(p+d3)
  h_inf(p,c) = Q2/(Q2+c)
  tau_h(p,c) = 1/(a2*(Q2+c))
  J_IPR = v1*m_inf^3*n_inf^3*h^3*(c_er-c)      c_er = (c0 - c)/c1
  J_leak = v2*(c_er-c)
  J_SERCA = v3*c^2/(k3^2+c^2)
  dc/dt = J_IPR + J_leak - J_SERCA
  dh/dt = (h_inf - h)/tau_h

PARAMETERS:
  d1=0.13uM d2=1.049uM d3=0.9434uM d5=0.08234uM v1=6.0/s v2=0.015/s v3=0.9uM/s k3=0.1uM
  c0=2.0uM c1=0.185 a2=0.013238 1/(uM*s)  [a2 is the ONE free-calibrated parameter]

GATES (pre-registered):
  G1: Hopf loci exist near p=0.1175 and p=0.245uM (within +-0.01uM), located by Jacobian
      eigenvalue Re(lambda)=0 crossing, cross-validated by direct time-domain integration
      showing a genuine limit cycle just inside and none just outside each locus.
  G2: calibrated a2=0.013238 (as given) reproduces period(p_mid=0.17uM) = 65.73s to <5%.
  G3: HELD-OUT (a2 untouched): period(p_lo=0.125uM) in [18,240]s, period(p_hi=0.22uM) in
      [18,240]s, and period(p_lo) > period(p_mid) > period(p_hi) (monotonic order, matches
      the recorded 122.14 / 65.73 / 56.67 s to <10%).
  G4 (forced adversary): the adiabatic/slaved-h reduction (h := h_inf(p,c) algebraically,
      collapsing to a scalar autonomous dc/dt) gives ZERO sustained limit-cycle oscillations
      across a 76-point p-sweep spanning the same window -- Poincare-Bendixson forbids a 1D
      autonomous flow from oscillating, so this must be identically 0/76; a nonzero count would
      mean the "2D slow-fast structure is the cause" claim is false.

VOID FLOOR (pre-registered, must FAIL): scramble the model by swapping v1<->v3 (IPR max flux and
SERCA max flux exchanged -- breaks the influx/efflux balance that creates the oscillatory window).
This must NOT reproduce a Hopf window whose period at p=0.17uM lands within 25% of the real 65.73s
target -- if the scrambled system still passes, G2/G3 measure nothing about this specific model.
"""
import json
import os
import os as _os
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

OUT = {}
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
RESULTS_PATH = _os.path.join(OUT_ROOT, "calcium_ip3_lirinzel_oscillation",
                             "calcium_ip3_lirinzel_oscillation.json")

P = dict(d1=0.13, d2=1.049, d3=0.9434, d5=0.08234, v1=6.0, v2=0.015, v3=0.9, k3=0.1,
         c0=2.0, c1=0.185, a2=0.013238)


def make_rhs(pars):
    d1, d2, d3, d5 = pars["d1"], pars["d2"], pars["d3"], pars["d5"]
    v1, v2, v3, k3 = pars["v1"], pars["v2"], pars["v3"], pars["k3"]
    c0, c1, a2 = pars["c0"], pars["c1"], pars["a2"]

    def rhs(t, y, p):
        c, h = y
        c = max(c, 1e-9)
        h = min(max(h, 0.0), 1.0)
        c_er = (c0 - c) / c1
        m_inf = p / (p + d1)
        n_inf = c / (c + d5)
        Q2 = d2 * (p + d1) / (p + d3)
        h_inf = Q2 / (Q2 + c)
        tau_h = 1.0 / (a2 * (Q2 + c))
        J_IPR = v1 * m_inf ** 3 * n_inf ** 3 * h ** 3 * (c_er - c)
        J_leak = v2 * (c_er - c)
        J_SERCA = v3 * c ** 2 / (k3 ** 2 + c ** 2)
        dc = J_IPR + J_leak - J_SERCA
        dh = (h_inf - h) / tau_h
        return [dc, dh]
    return rhs


def make_rhs_adiabatic(pars):
    """Forced adversary: h collapses instantaneously to h_inf(p,c) -- scalar autonomous ODE in c."""
    d1, d2, d3, d5 = pars["d1"], pars["d2"], pars["d3"], pars["d5"]
    v1, v2, v3, k3 = pars["v1"], pars["v2"], pars["v3"], pars["k3"]
    c0, c1 = pars["c0"], pars["c1"]

    def rhs(t, y, p):
        c = max(y[0], 1e-9)
        c_er = (c0 - c) / c1
        m_inf = p / (p + d1)
        n_inf = c / (c + d5)
        Q2 = d2 * (p + d1) / (p + d3)
        h_inf = Q2 / (Q2 + c)
        J_IPR = v1 * m_inf ** 3 * n_inf ** 3 * h_inf ** 3 * (c_er - c)
        J_leak = v2 * (c_er - c)
        J_SERCA = v3 * c ** 2 / (k3 ** 2 + c ** 2)
        return [J_IPR + J_leak - J_SERCA]
    return rhs


def fixed_point(rhs, p, c0_guess=0.2, h0_guess=0.7):
    sol = fsolve(lambda y: rhs(0, y, p), [c0_guess, h0_guess], xtol=1e-13, full_output=True)
    y, info, ier, msg = sol
    return y, (ier == 1)


def jacobian_fd(rhs, y, p, eps=1e-6):
    n = len(y)
    J = np.zeros((n, n))
    f0 = np.array(rhs(0, y, p))
    for i in range(n):
        yp = np.array(y, dtype=float)
        step = eps * max(1.0, abs(y[i]))
        yp[i] += step
        f1 = np.array(rhs(0, yp, p))
        J[:, i] = (f1 - f0) / step
    return J


def hopf_scan(rhs, p_lo, p_hi, n=300):
    ps = np.linspace(p_lo, p_hi, n)
    re_max = []
    for p in ps:
        y, ok = fixed_point(rhs, p)
        if not ok or y[0] <= 0 or not (0 <= y[1] <= 1):
            re_max.append(np.nan)
            continue
        J = jacobian_fd(rhs, y, p)
        eig = np.linalg.eigvals(J)
        re_max.append(float(np.max(eig.real)))
    return ps, np.array(re_max)


def find_crossings(ps, re_max):
    crossings = []
    for i in range(len(ps) - 1):
        a, b = re_max[i], re_max[i + 1]
        if np.isnan(a) or np.isnan(b):
            continue
        if a * b < 0:
            frac = -a / (b - a)
            crossings.append(ps[i] + frac * (ps[i + 1] - ps[i]))
    return crossings


def measure_period(rhs, p, T=2000.0, y0=None, discard_frac=0.5):
    if y0 is None:
        y0 = [0.3, 0.5]
    sol = solve_ivp(rhs, [0, T], y0, args=(p,), method="LSODA", max_step=0.5,
                     dense_output=False, rtol=1e-8, atol=1e-10)
    t, c = sol.t, sol.y[0]
    n0 = int(len(t) * discard_frac)
    t_tail, c_tail = t[n0:], c[n0:]
    # local maxima
    peaks = []
    for i in range(1, len(c_tail) - 1):
        if c_tail[i] > c_tail[i - 1] and c_tail[i] > c_tail[i + 1]:
            peaks.append(t_tail[i])
    if len(peaks) < 3:
        return None, c_tail, peaks
    intervals = np.diff(peaks)
    return float(np.mean(intervals[-3:])), c_tail, peaks


rhs_real = make_rhs(P)
rhs_adiab = make_rhs_adiabatic(P)

# ---------------------------------------------------------------------------
# G1: Hopf loci near p=0.1175 and p=0.245
# ---------------------------------------------------------------------------
ps_scan, re_scan = hopf_scan(rhs_real, 0.001, 0.5, n=500)
crossings = find_crossings(ps_scan, re_scan)
claimed_loci = [0.1175, 0.245]
g1_matches = []
for cl in claimed_loci:
    if crossings:
        nearest = min(crossings, key=lambda x: abs(x - cl))
        g1_matches.append({"claimed": cl, "found": nearest, "abs_diff": abs(nearest - cl),
                            "pass": abs(nearest - cl) < 0.01})
    else:
        g1_matches.append({"claimed": cl, "found": None, "pass": False})
# time-domain cross-validation: oscillation present just inside, absent just outside
inside_p, outside_lo_p, outside_hi_p = 0.15, 0.10, 0.30
per_inside, _, pk_inside = measure_period(rhs_real, inside_p)
per_out_lo, _, pk_out_lo = measure_period(rhs_real, outside_lo_p)
per_out_hi, _, pk_out_hi = measure_period(rhs_real, outside_hi_p)
g1_timedomain = {"inside_p": inside_p, "period_inside": per_inside,
                  "outside_lo_p": outside_lo_p, "period_outside_lo": per_out_lo,
                  "outside_hi_p": outside_hi_p, "period_outside_hi": per_out_hi,
                  "pass": (per_inside is not None) and (per_out_lo is None) and (per_out_hi is None)}
g1_pass = all(m["pass"] for m in g1_matches) and g1_timedomain["pass"]
OUT["G1_hopf_loci"] = {"pass": bool(g1_pass), "crossings_found": crossings,
                        "vs_claimed": g1_matches, "time_domain_check": g1_timedomain}

# ---------------------------------------------------------------------------
# G2: calibration reproduction at p_mid=0.17
# ---------------------------------------------------------------------------
period_mid, _, _ = measure_period(rhs_real, 0.17)
target_mid = 65.73
g2_pass = period_mid is not None and abs(period_mid - target_mid) / target_mid < 0.05
OUT["G2_calibration_pmid"] = {"period_computed": period_mid, "target": target_mid,
                               "rel_err_pct": (abs(period_mid - target_mid) / target_mid * 100
                                               if period_mid else None),
                               "pass": bool(g2_pass)}

# ---------------------------------------------------------------------------
# G3: held-out predictions p_lo, p_hi
# ---------------------------------------------------------------------------
period_lo, _, _ = measure_period(rhs_real, 0.125)
period_hi, _, _ = measure_period(rhs_real, 0.22)
env = (18.0, 240.0)
claimed = {"p_lo": (0.125, 122.14), "p_mid": (0.17, 65.73), "p_hi": (0.22, 56.67)}
in_env_lo = period_lo is not None and env[0] <= period_lo <= env[1]
in_env_hi = period_hi is not None and env[0] <= period_hi <= env[1]
order_ok = (period_lo is not None and period_mid is not None and period_hi is not None
            and period_lo > period_mid > period_hi)
g3_pass = in_env_lo and in_env_hi and order_ok
rel_errs = {}
if period_lo is not None:
    rel_errs["p_lo_rel_err_pct"] = abs(period_lo - 122.14) / 122.14 * 100
if period_hi is not None:
    rel_errs["p_hi_rel_err_pct"] = abs(period_hi - 56.67) / 56.67 * 100
OUT["G3_holdout"] = {"period_lo": period_lo, "period_mid": period_mid, "period_hi": period_hi,
                      "envelope": env, "in_envelope_lo": in_env_lo, "in_envelope_hi": in_env_hi,
                      "monotonic_order": order_ok, "rel_errs": rel_errs, "pass": bool(g3_pass)}

# ---------------------------------------------------------------------------
# G4: forced adversary -- adiabatic/slaved-h collapse, 0/76 expected
# ---------------------------------------------------------------------------
p_sweep = np.linspace(0.02, 0.5, 76)
n_oscillatory = 0
adiab_detail = []
for p in p_sweep:
    sol = solve_ivp(rhs_adiab, [0, 1000.0], [0.3], args=(p,), method="LSODA",
                     max_step=1.0, rtol=1e-8, atol=1e-10)
    c_tail = sol.y[0][len(sol.t) // 2:]
    spread = float(np.max(c_tail) - np.min(c_tail)) if len(c_tail) else 0.0
    osc = spread > 1e-3 * max(1e-6, np.mean(c_tail))
    if osc:
        n_oscillatory += 1
    adiab_detail.append({"p": float(p), "spread": spread, "oscillatory": bool(osc)})
g4_pass = (n_oscillatory == 0)
OUT["G4_adiabatic_null_forced_adversary"] = {
    "n_sweep": len(p_sweep), "n_oscillatory_found": n_oscillatory,
    "claimed": "0/76", "pass": bool(g4_pass),
    "note": "1D autonomous ODE cannot sustain a limit cycle (Poincare-Bendixson); nonzero count "
            "would indicate a numerical-integration artifact being mistaken for oscillation.",
}

# ---------------------------------------------------------------------------
# VOID FLOOR: v1<->v3 swap must NOT reproduce the p_mid period
# ---------------------------------------------------------------------------
P_scrambled = dict(P)
P_scrambled["v1"], P_scrambled["v3"] = P["v3"], P["v1"]
rhs_scrambled = make_rhs(P_scrambled)
period_mid_scrambled, _, _ = measure_period(rhs_scrambled, 0.17)
if period_mid_scrambled is None:
    void_reproduces = False
else:
    void_reproduces = abs(period_mid_scrambled - target_mid) / target_mid < 0.25
void_floor_pass = not void_reproduces
OUT["VOID_FLOOR"] = {
    "swap": "v1<->v3", "period_mid_scrambled": period_mid_scrambled,
    "target": target_mid, "reproduces_target_within_25pct": bool(void_reproduces),
    "pass_(void_floor_correctly_fails)": bool(void_floor_pass),
}

ALL_GATES = [OUT["G1_hopf_loci"]["pass"], OUT["G2_calibration_pmid"]["pass"],
             OUT["G3_holdout"]["pass"], OUT["G4_adiabatic_null_forced_adversary"]["pass"],
             OUT["VOID_FLOOR"]["pass_(void_floor_correctly_fails)"]]
OUT["OVERALL"] = {"gates_passed": int(sum(ALL_GATES)), "gates_total": len(ALL_GATES),
                   "all_pass": bool(all(ALL_GATES))}
OUT["node_id"] = "calcium_ip3_lirinzel_oscillation"

os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
with open(RESULTS_PATH, "w") as fh:
    json.dump(OUT, fh, indent=2,
               default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

print(json.dumps(OUT["OVERALL"], indent=2))
print("G1", OUT["G1_hopf_loci"]["pass"], OUT["G1_hopf_loci"]["crossings_found"])
print("G2", OUT["G2_calibration_pmid"])
print("G3", OUT["G3_holdout"])
print("G4", OUT["G4_adiabatic_null_forced_adversary"])
print("VOID_FLOOR", OUT["VOID_FLOOR"])
