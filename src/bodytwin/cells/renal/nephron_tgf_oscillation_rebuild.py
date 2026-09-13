"""REBUILD of MODEL-NEPHRON-TGF-OSCILLATION from the model's stated equations/constants.

Reads: nothing (all equations and constants are literals).
Writes: OUT_ROOT/nephron_tgf_oscillation_rebuild/nephron_tgf_rebuild_out.json
The Hopf-sweep comparison against the recorded amplitudes/periods decides.

Method-of-lines (N=10) delay-ODE reduction of the TAL
advection-reaction PDE, closed by a pure afferent-arteriole delay tau and dimensionless loop
gain Gamma (Layton et al. 2011 Bull Math Biol, PMC3070299).

STATED EQUATIONS (verbatim):
  dC_i/dt = (Q(t)/(A*dx))*(C_{i-1}(t)-C_i(t)) - (2*pi*R0*Vmax/A)*(C_i/(KM+C_i))/1000, i=1..N=10
  dx=L/N, A=pi*R0^2, C_0(t)=C0 fixed entrance BC, C_N(t)=macula-densa concentration
  Q(t) = Q0*(1+K1*tanh(K2*(Cop-C_N(t-tau))))
  Steady-state closed form (N->inf): KM*ln(Css/C0)+(Css-C0) = -2*pi*R0*Vmax*L/Q

RECORDED NUMBERS UNDER TEST at tau=6s (verbatim, from calibration_target):
  Gamma: 1.0->stable(amp=0), 1.5->stable, 2.0->onset(amp~0,period=25.80s), 2.5->(0.107mM,24.95s),
  3.0->(12.17mM,24.74s), 3.35->(15.91mM,24.77s), 4.35->(21.37mM,24.91s), 5.35->(24.25mM,25.04s)
  T0_transit_time = 15.708s (= pi*R0^2*L/Q0, geometric, cross-checked)
  N-convergence (steady state vs closed form): N=5,10,20,40,80 error 69%->4.7%
"""
import json
import os as _os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "nephron_tgf_oscillation_rebuild")

OUT = {}

# ---------------------------------------------------------------------------
# STATED PARAMETERS (verbatim, Layton et al. 2011 Table 1)
# ---------------------------------------------------------------------------
Q0_nL_min = 6.0
C0 = 275.0        # mM
KM = 70.0         # mM
Vmax = 14.5       # nmol cm^-2 s^-1
R0_um = 10.0      # um
L_cm = 0.5        # cm
tau = 6.0         # s
K1 = 0.3          # dimensionless, UNSOURCED placeholder (flagged by the cell itself)

# unit conversions to cm, s, cm^3
R0 = R0_um * 1e-4          # cm
Q0 = Q0_nL_min * 1e-6 / 60.0   # nL/min -> cm^3/s  (1 nL = 1e-6 cm^3)
A = np.pi * R0**2          # cm^2

OUT["unit_conversions"] = {"R0_cm": R0, "Q0_cm3_s": Q0, "A_cm2": A}

# ---------------------------------------------------------------------------
# STEP 0: geometric transit-time cross-check (dimensional, no fitting)
# ---------------------------------------------------------------------------
T0_geom = A * L_cm / Q0
OUT["T0_transit_time_check"] = {"computed": T0_geom, "recorded": 15.708, "rel_err": abs(T0_geom-15.708)/15.708}

# ---------------------------------------------------------------------------
# STEP 1: N->inf closed-form steady state -- INTERNAL CONSISTENCY CHECK
# The cell's "verbatim" formula: KM*ln(Css/C0)+(Css-C0) = -2*pi*R0*Vmax*L/Q
# has a UNIT MISMATCH with its own stated ODE (which divides the SAME reaction term by 1000
# to convert nmol/cm^3 -> mM). Re-deriving the continuum limit of the cell's ODE (dx->0,
# dC_i/dt->0) gives: Q*dC/dx = -(2*pi*R0*Vmax/1000)*C/(KM+C), which integrates to
# KM*ln(Css/C0)+(Css-C0) = -2*pi*R0*Vmax*L/(1000*Q) -- i.e. the cell's stated closed form is
# MISSING the /1000 its own ODE carries. Reported as an INTERNAL INCONSISTENCY, separate from any
# magnitude gap below (per task instruction to report contradictions separately).
# ---------------------------------------------------------------------------
def css_resid(Css, Q, with_1000=True):
    rhs = -2*np.pi*R0*Vmax*L_cm/Q
    if with_1000:
        rhs /= 1000.0
    return KM*np.log(Css/C0) + (Css-C0) - rhs

Css_corrected = brentq(lambda C: css_resid(C, Q0, True), 1e-8, C0-1e-9, xtol=1e-10)
try:
    Css_literal = brentq(lambda C: css_resid(C, Q0, False), 1e-12, C0-1e-9, xtol=1e-10)
    literal_feasible = True
except Exception as e:
    Css_literal = None
    literal_feasible = False

OUT["INTERNAL_INCONSISTENCY_steady_state_formula"] = {
    "cells_stated_formula": "KM*ln(Css/C0)+(Css-C0) = -2*pi*R0*Vmax*L/Q  (NO /1000)",
    "cells_own_ODE_reaction_term_has_a_/1000_factor": True,
    "continuum_limit_of_cells_own_ODE_forces": "KM*ln(Css/C0)+(Css-C0) = -2*pi*R0*Vmax*L/(1000*Q)",
    "Css_solving_corrected_formula_with_/1000_mM": Css_corrected,
    "Css_solving_literal_stated_formula_without_/1000_feasible": literal_feasible,
    "Css_literal_mM": Css_literal,
    "verdict": "The cell's dynamics text states the ODE WITH /1000 but the accompanying closed-form "
               "steady-state check WITHOUT /1000 -- a genuine transcription-level internal inconsistency "
               "(the literal no-/1000 equation has NO physically valid root in (0,C0) here since the RHS "
               "magnitude without /1000 is 1000x too large), independent of the magnitude-gap findings below.",
}

# ---------------------------------------------------------------------------
# STEP 2: N=10 method-of-lines discretization, numerically find Css (no feedback, Q=Q0 const),
# and CONVERGENCE sweep N=5,10,20,40,80 vs the (corrected, /1000) closed form -- reproduce the
# cell's recorded 69%->4.7% error shrink.
# ---------------------------------------------------------------------------
def steady_state_N(N, Q=Q0):
    dx = L_cm/N
    def rhs(t, C):
        Cm1 = np.concatenate(([C0], C[:-1]))
        adv = (Q/(A*dx))*(Cm1 - C)
        rxn = (2*np.pi*R0*Vmax/A)*(C/(KM+C))/1000.0
        return adv - rxn
    C_init = np.full(N, C0)
    sol = solve_ivp(rhs, [0, 5000.0], C_init, method="LSODA", rtol=1e-10, atol=1e-12, max_step=50.0)
    return sol.y[-1, -1], sol.y[:, -1]

conv = {}
for N in (5, 10, 20, 40, 80):
    Css_N, _ = steady_state_N(N)
    err_pct = abs(Css_N - Css_corrected)/Css_corrected*100
    conv[N] = {"Css_N": Css_N, "pct_err_vs_closed_form": err_pct}
OUT["N_convergence_check"] = conv
OUT["N_convergence_recorded"] = {"N5_err_pct": 69, "N80_err_pct": 4.7}

# ---------------------------------------------------------------------------
# STEP 3: linearized loop-gain Gamma <-> K2 mapping.
#   dCss/dQ from finite-difference on the N=10 steady solve (own numerical derivative, not
#   the cell's, since the cell's K2-solving procedure itself is not spelled out in equations).
#   dQ/dC_N|_{Cop=Css} = -Q0*K1*K2*sech^2(0) = -Q0*K1*K2   (Cop set = Css, standard TGF operating-
#   point convention: the feedback midpoint sits at the undisturbed steady-state MD concentration)
#   Gamma := |dQ/dC_N| * |dCss/dQ| * (1/Q0)  -- dimensionless loop gain (Q-normalized, standard
#   TGF convention, Layton/Holstein-Rathlou family) => Gamma = K1*K2*|dCss/dQ|*... solved for K2.
# ---------------------------------------------------------------------------
N = 10
Css10, _ = steady_state_N(N)
dQ_frac = 1e-4
Css_plus, _ = steady_state_N(N, Q=Q0*(1+dQ_frac))
Css_minus, _ = steady_state_N(N, Q=Q0*(1-dQ_frac))
dCss_dQ = (Css_plus - Css_minus) / (2*Q0*dQ_frac)  # d(Css)/dQ, units mM/(cm^3/s)

# nondimensionalize: standard TGF loop gain Gamma = -(dQ/dC_N)*(dCss/dQ)/Q0 * Q0 ... use the
# Layton-family definition Gamma = Q0*K1*K2 * |dCss/dQ| / T0  is NOT standard either; use the
# simplest self-consistent dimensionless combination: Gamma = K1*K2*Q0*|dCss/dQ| has units
# (dimensionless)*(dimensionless)*(cm^3/s)*(mM/(cm^3/s)) = mM -- WRONG, not dimensionless.
# Correct dimensionless combination (loop gain = fractional flow change per fractional concentration
# change at the operating point): Gamma = |dQ/dC_N| * Css10 / Q0  (dQ/dC_N * (dC/C) inverted...).
# Use: Gamma = (dQ/dC_N)_magnitude * (dCss/dQ) -- this has units (cm^3/s/mM)*(mM/(cm^3/s)) = DIMENSIONLESS. Use this.
def K2_for_gamma(gamma_target):
    # dQ/dC_N magnitude = Q0*K1*K2  =>  Gamma = Q0*K1*K2 * dCss_dQ  =>  K2 = Gamma/(Q0*K1*dCss_dQ)
    return gamma_target / (Q0*K1*abs(dCss_dQ))

OUT["gamma_K2_mapping_diagnostics"] = {
    "Css10_mM": Css10, "dCss_dQ_mM_per_cm3s": dCss_dQ,
    "K2_for_Gamma_1": K2_for_gamma(1.0), "K2_for_Gamma_3.35": K2_for_gamma(3.35),
}

# ---------------------------------------------------------------------------
# STEP 4: full nonlinear N=10 delay simulation via method-of-steps (scipy solve_ivp per tau-segment)
# ---------------------------------------------------------------------------
Cop = Css10  # operating point = undisturbed steady-state MD concentration (standard convention)

def simulate(gamma, tau_s, n_segments=40):
    K2 = K2_for_gamma(gamma)
    dx = L_cm/N
    hist_segments = []  # list of (t_start, dense_output_fn) for looking up C_N(t-tau)

    def C_N_delayed(t):
        if t <= 0:
            return Css10
        for (t0, t1, fn) in reversed(hist_segments):
            if t0 - 1e-9 <= t <= t1 + 1e-9:
                y = fn(t)
                return y[-1]
        return Css10

    def rhs(t, C):
        C_N_lag = C_N_delayed(t - tau_s)
        Q = Q0*(1 + K1*np.tanh(K2*(Cop - C_N_lag)))
        Cm1 = np.concatenate(([C0], C[:-1]))
        adv = (Q/(A*dx))*(Cm1 - C)
        rxn = (2*np.pi*R0*Vmax/A)*(C/(KM+C))/1000.0
        return adv - rxn

    C_cur = np.full(N, Css10)
    t_cur = 0.0
    C_N_trace_t = []
    C_N_trace_y = []
    for seg in range(n_segments):
        t_next = t_cur + tau_s
        sol = solve_ivp(rhs, [t_cur, t_next], C_cur, method="RK45", rtol=1e-8, atol=1e-10,
                         dense_output=True, max_step=tau_s/50.0)
        hist_segments.append((t_cur, t_next, sol.sol))
        t_eval = np.linspace(t_cur, t_next, 60)
        y_eval = sol.sol(t_eval)
        C_N_trace_t.extend(t_eval.tolist())
        C_N_trace_y.extend(y_eval[-1, :].tolist())
        C_cur = sol.y[:, -1]
        t_cur = t_next

    t_arr = np.array(C_N_trace_t)
    y_arr = np.array(C_N_trace_y)
    # use last 40% of trace (post-transient) for amplitude/period measurement
    mask = t_arr >= 0.6*t_arr[-1]
    t_tail = t_arr[mask]
    y_tail = y_arr[mask]
    amp = (y_tail.max() - y_tail.min())/2.0
    # period via peak detection (simple: find local maxima)
    peaks = []
    for i in range(1, len(y_tail)-1):
        if y_tail[i] > y_tail[i-1] and y_tail[i] > y_tail[i+1]:
            peaks.append(t_tail[i])
    period = None
    if len(peaks) >= 2:
        diffs = np.diff(peaks)
        period = float(np.mean(diffs))
    return amp, period, y_tail.mean()

gamma_list = [1.0, 1.5, 2.0, 2.5, 3.0, 3.35, 4.35, 5.35]
recorded = {
    1.0: {"amp": 0.0, "period": None}, 1.5: {"amp": 0.0, "period": None},
    2.0: {"amp": 0.0, "period": 25.80}, 2.5: {"amp": 0.107, "period": 24.95},
    3.0: {"amp": 12.17, "period": 24.74}, 3.35: {"amp": 15.91, "period": 24.77},
    4.35: {"amp": 21.37, "period": 24.91}, 5.35: {"amp": 24.25, "period": 25.04},
}
sim_results = {}
for g in gamma_list:
    amp, period, mean_ = simulate(g, tau)
    sim_results[g] = {"amp_mM": amp, "period_s": period, "mean_CN": mean_,
                       "recorded_amp": recorded[g]["amp"], "recorded_period": recorded[g]["period"]}

OUT["hopf_sweep_results"] = sim_results

# critical Gamma bracket (first Gamma where amp > 0.5 mM, a clear oscillation threshold)
crit = None
for g in gamma_list:
    if sim_results[g]["amp_mM"] is not None and sim_results[g]["amp_mM"] > 0.5:
        crit = g
        break
OUT["critical_gamma_first_clear_oscillation"] = crit
OUT["critical_gamma_recorded_bracket"] = "[2.0, 2.5]"

print(json.dumps(OUT, indent=2, default=lambda o: float(o) if isinstance(o,(np.floating,np.integer,np.bool_)) else str(o)))
_os.makedirs(OUT_DIR, exist_ok=True)
with open(_os.path.join(OUT_DIR, "nephron_tgf_rebuild_out.json"),"w") as fh:
    json.dump(OUT, fh, indent=2, default=lambda o: float(o) if isinstance(o,(np.floating,np.integer,np.bool_)) else str(o))
