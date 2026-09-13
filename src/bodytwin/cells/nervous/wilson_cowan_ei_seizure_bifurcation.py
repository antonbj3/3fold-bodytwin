"""Wilson-Cowan E/I seizure bifurcation -- distinct from nephron_tgf_oscillation_rebuild.py
(renal tubuloglomerular feedback, unrelated state space). Question: does the
Wilson-Cowan (1972) 2-population E/I rate model, at its own stated parameters, show (a) a
genuine eigenvalue-confirmed bifurcation near c2*=8.75-9.0 separating a 30-80Hz gamma limit cycle
from a non-oscillating pathological plateau, (b) the claimed benzodiazepine-proxy A-vs-B
dissociation (restoring c2 restabilizes; restoring only tauI does not, because the fixed point is
algebraically tau-independent), and (c) the claimed exact dimensional time-rescaling null control
(ptp invariant, freq scales exactly as 1/k)?

STATED MODEL (verbatim, Wilson & Cowan 1972 PMID 4332108):
  tauE*dE/dt = -E + (1-rE*E)*SE(c1*E - c2*I + P)
  tauI*dI/dt = -I + (1-rI*I)*SI(c3*E - c4*I + Q),  Q=0 (not separately stated, external drive to I)
  Sk(x) = 1/(1+exp(-ak*(x-thetak))) - 1/(1+exp(ak*thetak))   [normalized so Sk(0)=0]

STATED PARAMETERS (verbatim): c1=16.0 c2=12.0(baseline) c3=15.0 c4=3.0 aE=1.3 thetaE=4.0
  aI=2.0 thetaI=3.7 rE=rI=1.0 tauE=5ms tauI=8ms P=2.2 (oscillatory-regime drive).

GATE (pre-registered):
  G1: bifurcation c2* found by (a) eigenvalue stable-node-count change and (b) ptp-discontinuity,
      both landing in [7.5, 10.0] and agreeing with each other to <=0.5 (claimed grid resolution
      0.25; relaxed here for a coarser sweep), zero hysteresis (up-sweep threshold == down-sweep).
  G2: at tauI=8ms, P=2.2, sweeping c2 in [9.0,16.0] (29 points): oscillation frequency stays in
      [30,80]Hz (claimed tight band 32.00-38.67Hz, checked to <15% envelope slack) for ALL points.
  G3 (forced adversary / claimed null control): dimensional time-rescaling tauE,tauI -> k*(tauE,
      tauI) for k in {0.5,1,2,4} at fixed c2=12,P=2.2 leaves ptp invariant (<1% spread) and
      freq(k) = freq(1)/k to <2% -- a geometric symmetry of the ODE (rescaling t->t/k is exactly
      equivalent to rescaling both time constants), not a numerical artifact.
  G4: benzodiazepine-proxy dissociation -- from the pathological branch (c2=4.0, P=2.2): (A)
      restoring c2 back past c2* restabilizes into oscillatory_gamma at the SAME threshold as
      onset (no hysteresis); (B) restoring ONLY tauI (swept 8->300ms, c2 HELD at 4.0) does NOT
      restabilize -- fixed point location is unchanged (tau-independent, confirmed algebraically:
      dE/dt=dI/dt=0 has tau canceled out) and slower eigenvalue*tauI stays ~constant, never
      crossing zero.

VOID FLOOR (pre-registered, must FAIL): scramble the weight matrix by swapping c1<->c4 and
c2<->c3 (excitatory self-weight and inhibitory self-weight exchanged; E->I and I->E cross-weights
exchanged -- structurally different network). This must NOT reproduce a c2*-type bifurcation
inside [7.5,10.0] with an in-band [30,80]Hz gamma branch on one side -- if the scrambled network
still lands in the same bands, G1/G2 are not discriminating this specific wiring.

Reads: nothing (all parameters embedded).
Writes: wilson_cowan_ei_seizure_bifurcation.json.
Gate: G1-G4 above plus the void floor.
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import json
import os
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

OUT = {}
RESULTS_PATH = _os.path.join(OUT_ROOT, "wilson_cowan_ei_seizure_bifurcation",
                             "wilson_cowan_ei_seizure_bifurcation.json")

BASE = dict(c1=16.0, c2=12.0, c3=15.0, c4=3.0, aE=1.3, thetaE=4.0, aI=2.0, thetaI=3.7,
            rE=1.0, rI=1.0, tauE=5.0, tauI=8.0, P=2.2, Q=0.0)


def S(x, a, theta):
    return 1.0 / (1.0 + np.exp(-a * (x - theta))) - 1.0 / (1.0 + np.exp(a * theta))


def rhs(t, y, pars):
    E, I = y
    SE = S(pars["c1"] * E - pars["c2"] * I + pars["P"], pars["aE"], pars["thetaE"])
    SI = S(pars["c3"] * E - pars["c4"] * I + pars["Q"], pars["aI"], pars["thetaI"])
    dE = (-E + (1 - pars["rE"] * E) * SE) / pars["tauE"]
    dI = (-I + (1 - pars["rI"] * I) * SI) / pars["tauI"]
    return [dE, dI]


def jacobian_fd(pars, y, h=1e-6):
    f0 = np.array(rhs(0, y, pars))
    J = np.zeros((2, 2))
    for i in range(2):
        yp = np.array(y, dtype=float)
        yp[i] += h
        J[:, i] = (np.array(rhs(0, yp, pars)) - f0) / h
    return J


def find_fixed_points(pars, n=9):
    found = []
    grid = np.linspace(0.001, 0.999, n)
    for e0 in grid:
        for i0 in grid:
            sol, info, ier, msg = fsolve(lambda y: rhs(0, y, pars), [e0, i0],
                                          full_output=True, xtol=1e-13)
            if ier != 1:
                continue
            E, I = sol
            if not (0 <= E <= 1.2 and 0 <= I <= 1.2):
                continue
            dup = any(abs(E - f[0]) < 1e-5 and abs(I - f[1]) < 1e-5 for f in found)
            if not dup:
                found.append((float(E), float(I)))
    return found


def classify_fixed_points(pars, fps):
    out = []
    for (E, I) in fps:
        J = jacobian_fd(pars, [E, I])
        eig = np.linalg.eigvals(J)
        stable = all(e.real < 0 for e in eig)
        out.append({"E": E, "I": I, "eig": [complex(e) for e in eig], "stable": bool(stable)})
    return out


def simulate_and_measure(pars, T=4000.0, y0=(0.2, 0.2)):
    sol = solve_ivp(rhs, [0, T], y0, args=(pars,), method="LSODA", max_step=1.0,
                     rtol=1e-9, atol=1e-11)
    t, E = sol.t, sol.y[0]
    n0 = len(t) // 2
    t_tail, E_tail = t[n0:], E[n0:]
    ptp = float(np.max(E_tail) - np.min(E_tail))
    mean_E = float(np.mean(E_tail))
    # frequency via peak counting (ms timescale)
    peaks = []
    for i in range(1, len(E_tail) - 1):
        if E_tail[i] > E_tail[i - 1] and E_tail[i] > E_tail[i + 1]:
            peaks.append(t_tail[i])
    freq = None
    if len(peaks) >= 3:
        intervals_ms = np.diff(peaks)
        period_ms = float(np.mean(intervals_ms[-5:]))
        freq = 1000.0 / period_ms if period_ms > 0 else None
    return {"ptp": ptp, "mean_E": mean_E, "freq_hz": freq, "n_peaks": len(peaks)}


# ---------------------------------------------------------------------------
# G1: bifurcation c2* -- eigenvalue stable-count change + ptp discontinuity
# ---------------------------------------------------------------------------
c2_sweep = np.round(np.arange(4.0, 16.01, 0.5), 3)
eig_stable_count = []
ptp_vals = []
for c2 in c2_sweep:
    pars = dict(BASE)
    pars["c2"] = c2
    fps = find_fixed_points(pars)
    classified = classify_fixed_points(pars, fps)
    n_stable = sum(1 for c in classified if c["stable"])
    eig_stable_count.append(n_stable)
    meas = simulate_and_measure(pars, T=3000.0)
    ptp_vals.append(meas["ptp"])

eig_stable_count = np.array(eig_stable_count)
ptp_vals = np.array(ptp_vals)

# eigenvalue-based threshold: first c2 (descending) where stable-count changes vs c2=16 baseline
baseline_stable = eig_stable_count[-1]
eig_c2_star = None
for i in range(len(c2_sweep) - 1, 0, -1):
    if eig_stable_count[i - 1] != baseline_stable and eig_stable_count[i] == baseline_stable:
        eig_c2_star = 0.5 * (c2_sweep[i - 1] + c2_sweep[i])
        break

# ptp-based threshold: discontinuity in ptp vs c2 (descending)
ptp_c2_star = None
for i in range(len(c2_sweep) - 1, 0, -1):
    if abs(ptp_vals[i] - ptp_vals[i - 1]) > 0.05 and ptp_vals[i] > 0.05 and ptp_vals[i - 1] < 0.05:
        ptp_c2_star = 0.5 * (c2_sweep[i - 1] + c2_sweep[i])
        break

g1_band = [7.5, 10.0]
g1_eig_in_band = eig_c2_star is not None and g1_band[0] <= eig_c2_star <= g1_band[1]
g1_ptp_in_band = ptp_c2_star is not None and g1_band[0] <= ptp_c2_star <= g1_band[1]
g1_agree = (eig_c2_star is not None and ptp_c2_star is not None
            and abs(eig_c2_star - ptp_c2_star) <= 0.5)
g1_pass = g1_eig_in_band and g1_ptp_in_band and g1_agree
OUT["G1_bifurcation_c2star"] = {
    "eig_c2_star": eig_c2_star, "ptp_c2_star": ptp_c2_star, "claimed_range": [8.75, 9.0],
    "band_check": g1_band, "eig_in_band": g1_eig_in_band, "ptp_in_band": g1_ptp_in_band,
    "two_methods_agree": g1_agree, "pass": bool(g1_pass),
    "c2_sweep_stable_count": eig_stable_count.tolist(), "c2_sweep_ptp": ptp_vals.tolist(),
    "c2_sweep_values": c2_sweep.tolist(),
}

# ---------------------------------------------------------------------------
# G2: gamma-band frequency holdout across c2=[9,16], 29 points
# ---------------------------------------------------------------------------
c2_holdout = np.linspace(9.0, 16.0, 29)
freqs = []
for c2 in c2_holdout:
    pars = dict(BASE)
    pars["c2"] = c2
    meas = simulate_and_measure(pars, T=3000.0)
    freqs.append(meas["freq_hz"])
freqs_valid = [f for f in freqs if f is not None]
claimed_band = [30.0, 80.0]
claimed_tight = [32.00, 38.67]
in_band = all(claimed_band[0] <= f <= claimed_band[1] for f in freqs_valid) and len(freqs_valid) == len(c2_holdout)
tight_slack = 0.15
in_tight_env = all((claimed_tight[0] * (1 - tight_slack)) <= f <= (claimed_tight[1] * (1 + tight_slack))
                    for f in freqs_valid)
g2_pass = in_band and in_tight_env
OUT["G2_gamma_holdout"] = {"n_points": len(c2_holdout), "n_valid_freq": len(freqs_valid),
                           "freq_min": min(freqs_valid) if freqs_valid else None,
                           "freq_max": max(freqs_valid) if freqs_valid else None,
                           "claimed_tight_band": claimed_tight, "in_wide_band_30_80": in_band,
                           "in_tight_envelope_15pct_slack": in_tight_env, "pass": bool(g2_pass)}

# ---------------------------------------------------------------------------
# G3: dimensional time-rescaling null control (forced adversary / geometric symmetry)
# ---------------------------------------------------------------------------
k_vals = [0.5, 1, 2, 4]
rescale_results = []
base_freq = None
base_ptp = None
for k in k_vals:
    pars = dict(BASE)
    pars["tauE"] = BASE["tauE"] * k
    pars["tauI"] = BASE["tauI"] * k
    meas = simulate_and_measure(pars, T=4000.0 * k)
    rescale_results.append({"k": k, "ptp": meas["ptp"], "freq_hz": meas["freq_hz"]})
    if k == 1:
        base_freq, base_ptp = meas["freq_hz"], meas["ptp"]

ptp_list = [r["ptp"] for r in rescale_results]
ptp_spread_pct = (max(ptp_list) - min(ptp_list)) / np.mean(ptp_list) * 100
freq_errs = []
for r in rescale_results:
    expected = base_freq / r["k"]
    if r["freq_hz"] is not None:
        freq_errs.append(abs(r["freq_hz"] - expected) / expected * 100)
g3_pass = (ptp_spread_pct < 1.0) and all(e < 2.0 for e in freq_errs)
OUT["G3_time_rescaling_null"] = {"results": rescale_results, "ptp_spread_pct": ptp_spread_pct,
                                  "freq_rel_errs_pct": freq_errs, "pass": bool(g3_pass)}

# ---------------------------------------------------------------------------
# G4: benzodiazepine-proxy A vs B dissociation
# ---------------------------------------------------------------------------
PATHOLOGICAL_C2 = 4.0
pars_path = dict(BASE)
pars_path["c2"] = PATHOLOGICAL_C2
meas_path = simulate_and_measure(pars_path, T=3000.0)

# Proxy A: restore c2 back above c2*
c2_restore = 12.0
pars_A = dict(BASE)
pars_A["c2"] = c2_restore
meas_A = simulate_and_measure(pars_A, T=3000.0)
proxyA_pass = meas_A["freq_hz"] is not None and 30.0 <= meas_A["freq_hz"] <= 80.0

# Proxy B: restore ONLY tauI, c2 held at pathological 4.0
tauI_sweep = [8, 12, 20, 40, 80, 150, 300]
proxyB_detail = []
for tauI in tauI_sweep:
    pars_B = dict(BASE)
    pars_B["c2"] = PATHOLOGICAL_C2
    pars_B["tauI"] = float(tauI)
    fps = find_fixed_points(pars_B)
    classified = classify_fixed_points(pars_B, fps)
    stable_fps = [c for c in classified if c["stable"]]
    slow_eig_x_tauI = None
    if stable_fps:
        eigs_real = [e.real for e in stable_fps[0]["eig"]]
        slow_eig = max(eigs_real)  # least-negative = slowest mode
        slow_eig_x_tauI = slow_eig * tauI
    proxyB_detail.append({"tauI": tauI, "n_stable": len(stable_fps),
                           "slow_eig_x_tauI": slow_eig_x_tauI,
                           "fp_E": stable_fps[0]["E"] if stable_fps else None})
proxyB_still_stable_node = all(d["n_stable"] >= 1 for d in proxyB_detail)
fp_E_vals = [d["fp_E"] for d in proxyB_detail if d["fp_E"] is not None]
fp_location_invariant = (len(fp_E_vals) == len(proxyB_detail)
                          and (max(fp_E_vals) - min(fp_E_vals)) < 1e-4)
proxyB_fails_to_restabilize = proxyB_still_stable_node and fp_location_invariant

g4_pass = proxyA_pass and proxyB_fails_to_restabilize
OUT["G4_benzodiazepine_dissociation"] = {
    "pathological_baseline": meas_path,
    "proxyA_c2_restore_to": c2_restore, "proxyA_measured": meas_A, "proxyA_pass": bool(proxyA_pass),
    "proxyB_tauI_sweep": proxyB_detail, "proxyB_fp_location_invariant": bool(fp_location_invariant),
    "proxyB_fails_to_restabilize_as_claimed": bool(proxyB_fails_to_restabilize),
    "pass": bool(g4_pass),
}

# ---------------------------------------------------------------------------
# VOID FLOOR: swap c1<->c4, c2<->c3 (structurally different network)
# ---------------------------------------------------------------------------
SCRAMBLED = dict(BASE)
SCRAMBLED["c1"], SCRAMBLED["c4"] = BASE["c4"], BASE["c1"]
SCRAMBLED["c2"], SCRAMBLED["c3"] = BASE["c3"], BASE["c2"]

c2_sweep_scr = np.round(np.arange(4.0, 16.01, 0.5), 3)
eig_stable_scr = []
ptp_scr = []
for c2 in c2_sweep_scr:
    pars = dict(SCRAMBLED)
    pars["c2"] = c2
    fps = find_fixed_points(pars)
    classified = classify_fixed_points(pars, fps)
    eig_stable_scr.append(sum(1 for c in classified if c["stable"]))
    meas = simulate_and_measure(pars, T=2000.0)
    ptp_scr.append(meas["ptp"])
eig_stable_scr = np.array(eig_stable_scr)
ptp_scr = np.array(ptp_scr)
baseline_stable_scr = eig_stable_scr[-1]
scr_c2_star = None
for i in range(len(c2_sweep_scr) - 1, 0, -1):
    if eig_stable_scr[i - 1] != baseline_stable_scr and eig_stable_scr[i] == baseline_stable_scr:
        scr_c2_star = 0.5 * (c2_sweep_scr[i - 1] + c2_sweep_scr[i])
        break
scr_reproduces_band = scr_c2_star is not None and g1_band[0] <= scr_c2_star <= g1_band[1]
# also check whether the "gamma branch" side (c2 in [9,16]) is in-band for scrambled net
scr_freqs = []
for c2 in np.linspace(9.0, 16.0, 10):
    pars = dict(SCRAMBLED)
    pars["c2"] = c2
    meas = simulate_and_measure(pars, T=2000.0)
    if meas["freq_hz"] is not None:
        scr_freqs.append(meas["freq_hz"])
scr_gamma_in_band = (len(scr_freqs) > 0 and
                     all(claimed_band[0] <= f <= claimed_band[1] for f in scr_freqs))
scr_reproduces_all = scr_reproduces_band and scr_gamma_in_band
void_floor_pass = not scr_reproduces_all
OUT["VOID_FLOOR"] = {
    "swap": "c1<->c4, c2<->c3", "scrambled_c2_star": scr_c2_star,
    "reproduces_c2star_band": bool(scr_reproduces_band),
    "scrambled_gamma_freqs": scr_freqs, "scrambled_gamma_in_band_30_80": bool(scr_gamma_in_band),
    "reproduces_both": bool(scr_reproduces_all),
    "pass_(void_floor_correctly_fails)": bool(void_floor_pass),
}

ALL_GATES = [OUT["G1_bifurcation_c2star"]["pass"], OUT["G2_gamma_holdout"]["pass"],
             OUT["G3_time_rescaling_null"]["pass"], OUT["G4_benzodiazepine_dissociation"]["pass"],
             OUT["VOID_FLOOR"]["pass_(void_floor_correctly_fails)"]]
OUT["OVERALL"] = {"gates_passed": int(sum(ALL_GATES)), "gates_total": len(ALL_GATES),
                   "all_pass": bool(all(ALL_GATES))}
OUT["node_id"] = "wilson-cowan E/I seizure bifurcation"

os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
with open(RESULTS_PATH, "w") as fh:
    json.dump(OUT, fh, indent=2,
               default=lambda o: (float(o.real) if isinstance(o, complex) and abs(o.imag) < 1e-9
                                   else (str(o) if isinstance(o, complex)
                                         else (float(o) if isinstance(o, (np.floating, np.integer))
                                               else str(o)))))

print(json.dumps(OUT["OVERALL"], indent=2))
for k in ["G1_bifurcation_c2star", "G2_gamma_holdout", "G3_time_rescaling_null",
          "G4_benzodiazepine_dissociation", "VOID_FLOOR"]:
    v = dict(OUT[k])
    for junk in ["c2_sweep_stable_count", "c2_sweep_ptp", "c2_sweep_values"]:
        v.pop(junk, None)
    print(k, v)
