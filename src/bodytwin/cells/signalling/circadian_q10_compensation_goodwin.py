#!/usr/bin/env python3
"""circadian_q10_compensation_goodwin.py.

NODE RESOLVED: MODEL-CIRCADIAN-TEMPERATURE-COMPENSATION
QUESTION: does a from-scratch Goodwin (1965) 3-variable negative-feedback oscillator (mRNA -> protein
-> repressor, Hill-function repression), with EVERY rate constant scaled by a UNIFORM Q10 factor
(the naive null hypothesis for a biochemical oscillator with no special compensation mechanism),
produce a period that is nearly temperature-INSENSITIVE (small effective Q10_period, in the claim's
0.85-1.2 band) purely from the NEGATIVE-FEEDBACK LOOP STRUCTURE itself -- i.e. is temperature
compensation partly a structural/geometric consequence of a delayed-negative-feedback oscillator
(where period is set by loop delay, and if EVERY step speeds up by the same Q10 factor, the loop
delay itself scales by 1/Q10, but higher-order nonlinear saturation and multi-step delay partially
cancel that scaling) rather than requiring a bespoke anti-Q10 mechanism per step?

DISTINCT from circadian_rhythm.py (tau/DLMO/light-PRC falsifiers, an anchor-graph-fit
model of the human SCN clock output) -- this cell is a from-scratch NUMERICAL Goodwin-oscillator
simulation with NO parameters borrowed from that file or from the biological literature cited in
the claim (per the claim's text: "no parameters borrowed from the biological literature").

METHOD: dx/dt=v0/(1+(z/Km)^n) - k1*x ; dy/dt=k2*x - k3*y ; dz/dt=k4*y - k5*z (classic Goodwin form,
all k_i, v0 uniformly Q10-scaled by exp(ln(Q10)*(T-T0)/10) at each of several T values spanning a
20C physiological range) -- measure period vs T, compute the empirical Q10_period =
(period(T0)/period(T0+10))^1 in the standard Q10 convention, and check it falls near 1 (compensated)
despite EVERY microscopic rate constant individually having Q10=Q10_rate>>1 (e.g. 2.5, a typical
enzymatic Q10).

VOID FLOOR: a version with the negative feedback loop BROKEN (Hill exponent n->0, i.e. z no longer
represses x -- pure linear cascade, no oscillation, so "period" becomes undefined / the settling
time to steady state) must NOT show the same compensation -- the settling-time-to-steady-state of
the broken linear cascade should scale directly with 1/Q10_rate (i.e. Q10-sensitive), confirming
the compensation is a property of the CLOSED LOOP, not just "any multi-step biochemical process
scaled uniformly happens to look compensated."
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from pathlib import Path as _Path
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = _os.path.join(OUT_ROOT, "circadian_q10_compensation_goodwin", "circadian_q10_compensation_goodwin.json")

T0 = 25.0       # reference temperature, C
Q10_RATE = 2.5  # typical enzymatic per-step Q10 (the "naive" per-reaction speedup)

# nominal Goodwin parameters at T0 (arbitrary units, chosen to sit in the oscillatory Hill-n regime)
V0, KM, N_HILL = 1.0, 1.0, 8
K1, K2, K3, K4, K5 = 0.1, 0.2, 0.1, 0.2, 0.1

def q10_scale(T, q10):
    return q10 ** ((T - T0) / 10.0)

def goodwin_rhs(t, s, scale, n_hill):
    x, y, z = np.maximum(s, 0)
    dx = scale * V0 / (1 + (z / KM) ** n_hill) - scale * K1 * x
    dy = scale * K2 * x - scale * K3 * y
    dz = scale * K4 * y - scale * K5 * z
    return [dx, dy, dz]

def measure_period(T, n_hill=N_HILL, t_max=2000):
    scale = q10_scale(T, Q10_RATE)
    y0 = [1.0, 1.0, 1.0]
    # rescale integration horizon so we always cover enough OSCILLATION CYCLES even though the
    # loop speeds up at higher T: use a generous fixed real-time window, adaptive step.
    t_eval = np.linspace(0, t_max, 20000)
    sol = solve_ivp(goodwin_rhs, [0, t_max], y0, args=(scale, n_hill), t_eval=t_eval,
                     method="RK45", rtol=1e-9, atol=1e-11)
    x = sol.x if hasattr(sol, "x") else sol.y[0]
    t = sol.t
    late = t > t_max * 0.5
    x_late, t_late = sol.y[0][late], t[late]
    amp = x_late.max() - x_late.min()
    peaks, _ = find_peaks(x_late, prominence=amp * 0.1 if amp > 1e-9 else 1e-12)
    if len(peaks) >= 2:
        period = float(np.mean(np.diff(t_late[peaks])))
    else:
        period = None
    return period, float(amp)

def measure_settling_time_broken_loop(T, thresh=0.01, t_max=2000):
    """n_hill=0 breaks the feedback (constant production v0/2, no repression) -> pure linear
    relaxation to steady state; settling time = time to get within `thresh` fraction of steady
    state, a Q10-sensitive-if-uncompensated quantity for comparison."""
    scale = q10_scale(T, Q10_RATE)
    y0 = [1.0, 1.0, 1.0]
    t_eval = np.linspace(0, t_max, 20000)
    sol = solve_ivp(goodwin_rhs, [0, t_max], y0, args=(scale, 1e-9), t_eval=t_eval,
                     method="RK45", rtol=1e-9, atol=1e-11)
    x = sol.y[0]
    x_ss = x[-1]
    # settling time: first t where |x-x_ss| < thresh*|x0-x_ss| and stays there
    dev = np.abs(x - x_ss)
    tol = thresh * abs(y0[0] - x_ss) if abs(y0[0] - x_ss) > 1e-9 else thresh
    below = dev < tol
    idx = np.argmax(below) if below.any() else -1
    return float(sol.t[idx]) if idx >= 0 else None

def main():
    temps = [T0 - 10, T0, T0 + 10]  # 15C, 25C, 35C -- a 20C physiological span
    periods = {}
    for T in temps:
        p, a = measure_period(T)
        periods[T] = {"period": p, "amplitude": a}

    p_lo = periods[T0 - 10]["period"]
    p_hi = periods[T0 + 10]["period"]
    # empirical Q10 of the PERIOD itself (standard convention: rate ratio over 10C; for a period,
    # "rate" = 1/period, so Q10_period = (1/p_hi)/(1/p_lo) = p_lo/p_hi)
    q10_period = (p_lo / p_hi) if (p_lo and p_hi) else None

    # void floor: broken-loop settling time at same two temperatures
    settle_lo = measure_settling_time_broken_loop(T0 - 10)
    settle_hi = measure_settling_time_broken_loop(T0 + 10)
    q10_settle = (settle_lo / settle_hi) if (settle_lo and settle_hi) else None

    claim = {"q10_period_band": [0.85, 1.2], "per_step_q10_rate": Q10_RATE}

    gates = {
        "G1_oscillation_found_all_T": bool(p_lo is not None and p_hi is not None and periods[T0]["period"] is not None),
        "G2_period_q10_in_claimed_band": bool(q10_period is not None and claim["q10_period_band"][0] <= q10_period <= claim["q10_period_band"][1]),
        "G3_perstep_rate_q10_is_naive_high": bool(Q10_RATE > 2.0),
        "G4_voidfloor_broken_loop_settling_q10_far_from_1": bool(q10_settle is not None and abs(q10_settle - 1.0) > 0.5),
    }
    all_pass = bool(all(gates.values()))

    result = {
        "node": "MODEL-CIRCADIAN-TEMPERATURE-COMPENSATION",
        "measured": {
            "periods_by_T": periods,
            "q10_period_measured": q10_period,
            "broken_loop_settling_time_by_T": {str(T0 - 10): settle_lo, str(T0 + 10): settle_hi},
            "q10_settling_broken_loop": q10_settle,
        },
        "claim": claim,
        "gates": gates,
        "all_gates_pass": all_pass,
        "void_floor_note": "G4: with the negative-feedback loop broken (Hill n->0, pure linear "
                            "relaxation), the settling-time Q10 must stay FAR from 1 (Q10-sensitive) "
                            "-- confirms compensation is a property of the closed feedback loop, not "
                            "an artifact of uniformly rescaling any multi-step process.",
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result

if __name__ == "__main__":
    main()
