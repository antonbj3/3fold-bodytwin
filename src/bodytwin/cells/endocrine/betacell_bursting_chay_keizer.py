"""Beta-cell bursting (Chay-Keizer / Sherman-Rinzel-Keizer fast-slow ODE).

Answers: does the cell's stated
Chay-Keizer/Sherman-Rinzel-Keizer fast-slow ODE, at Bertram et al 2004 Table 1 conductances,
reproduce (a) the claimed 3-branch bifurcation structure (depolarized-block / oscillatory-spiking
/ silent-hyperpolarized) as Ca is scanned at fixed KATP-drive a, with the claimed spiking-window
bounds at a=0.001 and a=0.003, and (b) the cell's DISCLOSED NEGATIVE that the full 3-state
(V,n,Ca) system gives continuous modulated spiking rather than the clean discrete ~30s burst/
silence waveform real beta-cells show?

STATED MODEL (verbatim, node claim's DYNAMICS field): states V(mV), n(dimensionless, fast K+
gate), Ca(uM, slow). KATP-drive a held as an external bifurcation parameter (its own slow ODE,
tau_a=150000ms, is far slower than a burst period and is NOT integrated here -- consistent with
the claim's "Ca held fixed as parameter, fast V-n subsystem integrated" bifurcation-scan
method, extended one level to (V,n,Ca) for the burst-period test):
  Cm*dV/dt = -(ICa+IK+IKCa+IKATP)
  ICa = gCa*m_inf(V)*(V-VCa),  m_inf(V) = 1/(1+exp((vm-V)/sm))
  IK  = gK*n*(V-VK)
  IKCa = gKCa*(Ca/(Ca+KD_KCa))*(V-VK)
  IKATP = gKATP_bar*a*(V-VK)
  dn/dt = (n_inf(V)-n)/tau_n,  n_inf(V) = 1/(1+exp((vn-V)/sn))
  dCa/dt = f_cyt*(-alpha*ICa - k_PMCA*Ca)   [operational k_PMCA per honest_gaps]

STATED PARAMETERS (verbatim, Bertram2004 Table 1 lineage): Cm=5300fF gCa=1000pS gK=2700pS
  gKCa=600pS gKATP_bar=25000pS VCa=25mV VK=-75mV vm=-20mV sm=12mV vn=-16mV sn=5mV tau_n=20ms
  KD_KCa=0.5uM f_cyt=0.01 alpha=4.5e-6 uM/fA/ms k_PMCA_operational=0.05/ms.

GATE (pre-registered):
  G1: at BOTH a=0.001 and a=0.003, a Ca-scan (fast V-n subsystem, Ca held as parameter) shows the
      claimed 3 branches in order as Ca increases: depolarized-block (fixed pt, V high) ->
      oscillatory/spiking (ptp>=20mV limit cycle) -> silent-hyperpolarized (fixed pt, V low).
  G2: at a=0.001, spiking window bounds match claimed [0.10-0.11, 0.20]uM to <20% relative on
      each edge (a coarse-grid claim, checked at claim's resolution not over-precisely).
  G3: at a=0.003, spiking window bounds match claimed [0.05, 0.13]uM to <20% relative each edge.
  G4 (claim's DISCLOSED NEGATIVE, symmetric QC -- must ALSO reproduce as claimed, not tuned
      away): the full 3-state (V,n,Ca) system at a=0.001 gives CONTINUOUS modulated spiking, NOT
      a clean discrete burst/silent waveform -- measured by a duty-cycle discreteness statistic
      (fraction of inter-spike-intervals that fall in a long "silent" tail vs a short "active"
      cluster); the claim is REPRODUCED if this statistic shows NO clean bimodal separation
      (matching "gives continuous modulated spiking instead").

VOID FLOOR (pre-registered, must FAIL): scramble by swapping gCa<->gK (inward Ca conductance and
outward delayed-rectifier K+ conductance exchanged -- inverts which current dominates depolarizing
vs repolarizing). This must NOT reproduce a 3-branch structure with a spiking window inside
[0.02,0.30]uM at a=0.001 -- if the scrambled (inverted-current) cell still shows the claimed
window, G1/G2 are not discriminating this specific conductance assignment.
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

OUT = {}
RESULTS_PATH = os.path.join(OUT_ROOT, "betacell_bursting_chay_keizer",
                            "betacell_bursting_chay_keizer.json")

BASE = dict(Cm=5300.0, gCa=1000.0, gK=2700.0, gKCa=600.0, gKATP_bar=25000.0,
            VCa=25.0, VK=-75.0, vm=-20.0, sm=12.0, vn=-16.0, sn=5.0, tau_n=20.0,
            KD_KCa=0.5, f_cyt=0.01, alpha=4.5e-6, k_PMCA=0.05)


def m_inf(V, pars):
    return 1.0 / (1.0 + np.exp((pars["vm"] - V) / pars["sm"]))


def n_inf(V, pars):
    return 1.0 / (1.0 + np.exp((pars["vn"] - V) / pars["sn"]))


def currents(V, n, Ca, a, pars):
    ICa = pars["gCa"] * m_inf(V, pars) * (V - pars["VCa"])
    IK = pars["gK"] * n * (V - pars["VK"])
    IKCa = pars["gKCa"] * (Ca / (Ca + pars["KD_KCa"])) * (V - pars["VK"])
    IKATP = pars["gKATP_bar"] * a * (V - pars["VK"])
    return ICa, IK, IKCa, IKATP


def rhs_fast(t, y, Ca, a, pars):
    """Fast subsystem: Ca held as external bifurcation parameter."""
    V, n = y
    ICa, IK, IKCa, IKATP = currents(V, n, Ca, a, pars)
    dV = -(ICa + IK + IKCa + IKATP) / pars["Cm"]
    dn = (n_inf(V, pars) - n) / pars["tau_n"]
    return [dV, dn]


def rhs_full(t, y, a, pars):
    """Full 3-state (V,n,Ca) fast-slow system, a held fixed (slow-drive assumption)."""
    V, n, Ca = y
    Ca = max(Ca, 1e-6)
    ICa, IK, IKCa, IKATP = currents(V, n, Ca, a, pars)
    dV = -(ICa + IK + IKCa + IKATP) / pars["Cm"]
    dn = (n_inf(V, pars) - n) / pars["tau_n"]
    dCa = pars["f_cyt"] * (-pars["alpha"] * ICa - pars["k_PMCA"] * Ca)
    return [dV, dn, dCa]


def classify_fast(Ca, a, pars, T=12000.0, y0=(-50.0, 0.1)):
    """GEOMETRIC NOTE (OODA catch): a naive tail-ptp>threshold test is fooled by a slowly-decaying
    underdamped spiral near a stable focus -- its ptp rises smoothly with Ca just like a genuine
    limit cycle would, then cuts off, because it is decaying-oscillation amplitude, not sustained
    amplitude. Fixed by comparing ptp across TWO non-overlapping late windows: a real limit cycle
    has equal amplitude in both (ratio~1); a decaying spiral's second window is measurably smaller."""
    sol = solve_ivp(rhs_fast, [0, T], y0, args=(Ca, a, pars), method="LSODA",
                     max_step=1.0, rtol=1e-9, atol=1e-11)
    V = sol.y[0]
    n = len(V)
    win_a = V[int(n * 0.5):int(n * 0.75)]
    win_b = V[int(n * 0.75):]
    ptp_a = float(np.max(win_a) - np.min(win_a)) if len(win_a) else 0.0
    ptp_b = float(np.max(win_b) - np.min(win_b)) if len(win_b) else 0.0
    ptp = ptp_b
    mean_V = float(np.mean(win_b))
    sustained = (ptp_b > 1.0) and (ptp_a < 1e-6 or abs(ptp_b - ptp_a) / max(ptp_a, ptp_b) < 0.15)
    if ptp_b >= 20.0 and sustained:
        kind = "spiking"
    elif mean_V > -30.0:
        kind = "depolarized_block"
    else:
        kind = "silent_hyperpolarized"
    return kind, ptp, mean_V


def scan_window(a, pars, ca_lo=0.02, ca_hi=0.30, n=30):
    cas = np.linspace(ca_lo, ca_hi, n)
    results = []
    for ca in cas:
        kind, ptp, mean_V = classify_fast(ca, a, pars)
        results.append({"Ca": float(ca), "kind": kind, "ptp": ptp, "mean_V": mean_V})
    spiking = [r["Ca"] for r in results if r["kind"] == "spiking"]
    window = (min(spiking), max(spiking)) if spiking else None
    # 3-branch check: block region (low Ca) -> spiking -> silent region (high Ca), in that order
    kinds_seq = [r["kind"] for r in results]
    has_block_before = "depolarized_block" in kinds_seq[:kinds_seq.index("spiking")] if "spiking" in kinds_seq else False
    has_silent_after = ("silent_hyperpolarized" in kinds_seq[kinds_seq.index("spiking") + len([k for k in kinds_seq if k == "spiking"]):]
                         if "spiking" in kinds_seq else False)
    three_branch = ("spiking" in kinds_seq) and has_block_before and has_silent_after
    return {"scan": results, "window": window, "three_branch_confirmed": three_branch}


# ---------------------------------------------------------------------------
# G1 + G2 + G3: bifurcation scans at a=0.001 and a=0.003
# ---------------------------------------------------------------------------
scan_a001 = scan_window(0.001, BASE)
scan_a003 = scan_window(0.003, BASE)

g1_pass = scan_a001["three_branch_confirmed"] and scan_a003["three_branch_confirmed"]
OUT["G1_three_branch_structure"] = {
    "a_0.001_three_branch": scan_a001["three_branch_confirmed"],
    "a_0.003_three_branch": scan_a003["three_branch_confirmed"],
    "pass": bool(g1_pass),
}

claimed_window_a001 = (0.10, 0.20)  # main claim text says ~0.10-0.20; datapoints say 0.11-0.20
claimed_window_a003 = (0.05, 0.13)


def window_match(found, claimed, tol=0.20):
    if found is None:
        return False, None, None
    lo_err = abs(found[0] - claimed[0]) / claimed[0]
    hi_err = abs(found[1] - claimed[1]) / claimed[1]
    return (lo_err < tol and hi_err < tol), lo_err, hi_err


g2_ok, g2_lo_err, g2_hi_err = window_match(scan_a001["window"], claimed_window_a001)
g3_ok, g3_lo_err, g3_hi_err = window_match(scan_a003["window"], claimed_window_a003)
OUT["G2_window_a0.001"] = {"found_window": scan_a001["window"], "claimed": claimed_window_a001,
                           "lo_rel_err": g2_lo_err, "hi_rel_err": g2_hi_err, "pass": bool(g2_ok)}
OUT["G3_window_a0.003"] = {"found_window": scan_a003["window"], "claimed": claimed_window_a003,
                           "lo_rel_err": g3_lo_err, "hi_rel_err": g3_hi_err, "pass": bool(g3_ok)}

# ---------------------------------------------------------------------------
# G4: reproduce the claim's disclosed negative -- continuous modulated spiking,
# NOT a clean discrete burst/silent waveform, in the full 3-state system at a=0.001
# ---------------------------------------------------------------------------
sol_full = solve_ivp(rhs_full, [0, 60000.0], [-50.0, 0.1, 0.15], args=(0.001, BASE),
                      method="LSODA", max_step=2.0, rtol=1e-8, atol=1e-10)
t_full, V_full = sol_full.t, sol_full.y[0]
tail_idx = len(t_full) // 3
t_tail, V_tail = t_full[tail_idx:], V_full[tail_idx:]

# spike detection: local maxima above a threshold
spike_times = []
for i in range(1, len(V_tail) - 1):
    if V_tail[i] > V_tail[i - 1] and V_tail[i] > V_tail[i + 1] and V_tail[i] > -40.0:
        spike_times.append(t_tail[i])
isis = np.diff(spike_times) if len(spike_times) > 2 else np.array([])
if len(isis) > 3:
    isi_cv = float(np.std(isis) / np.mean(isis))
    # bimodality check: sort ISIs, look for a big gap separating a "short" cluster from a "long" cluster
    isis_sorted = np.sort(isis)
    gaps = np.diff(isis_sorted)
    max_gap = float(np.max(gaps)) if len(gaps) else 0.0
    max_gap_ratio = max_gap / (np.mean(isis) + 1e-9)
    # discrete bursting would show max_gap_ratio >> 1 (silent phase >> mean ISI); continuous
    # modulated spiking shows a smoothly varying ISI distribution, max_gap_ratio modest
    discrete_bimodal = max_gap_ratio > 3.0
else:
    isi_cv, max_gap_ratio, discrete_bimodal = None, None, False

claimed_continuous_not_discrete = True  # claim's disclosed finding
g4_reproduces_claim = (not discrete_bimodal)  # continuous modulated spiking, no clean burst gap
OUT["G4_disclosed_negative_check"] = {
    "n_spikes_detected_tail": len(spike_times), "isi_cv": isi_cv,
    "max_isi_gap_ratio": max_gap_ratio, "discrete_bimodal_bursting_found": bool(discrete_bimodal),
    "claim_says": "continuous modulated spiking, NOT clean discrete burst/silent waveform",
    "reproduces_claim_as_stated": bool(g4_reproduces_claim), "pass": bool(g4_reproduces_claim),
}

# ---------------------------------------------------------------------------
# VOID FLOOR: gCa <-> gK swap
# ---------------------------------------------------------------------------
SCRAMBLED = dict(BASE)
SCRAMBLED["gCa"], SCRAMBLED["gK"] = BASE["gK"], BASE["gCa"]
scan_scr = scan_window(0.001, SCRAMBLED)
scr_window = scan_scr["window"]
scr_reproduces = False
if scr_window is not None:
    ok, _, _ = window_match(scr_window, claimed_window_a001, tol=0.20)
    scr_reproduces = ok and scan_scr["three_branch_confirmed"]
void_floor_pass = not scr_reproduces
OUT["VOID_FLOOR"] = {
    "swap": "gCa<->gK", "scrambled_window": scr_window,
    "scrambled_three_branch": scan_scr["three_branch_confirmed"],
    "reproduces_claimed_window": bool(scr_reproduces),
    "pass_(void_floor_correctly_fails)": bool(void_floor_pass),
}

ALL_GATES = [OUT["G1_three_branch_structure"]["pass"], OUT["G2_window_a0.001"]["pass"],
             OUT["G3_window_a0.003"]["pass"], OUT["G4_disclosed_negative_check"]["pass"],
             OUT["VOID_FLOOR"]["pass_(void_floor_correctly_fails)"]]
OUT["OVERALL"] = {"gates_passed": int(sum(ALL_GATES)), "gates_total": len(ALL_GATES),
                   "all_pass": bool(all(ALL_GATES))}
OUT["node_id"] = "MODEL-BETACELL-BURSTING-SECRETION"
OUT["scope_note"] = ("Gates cover only the numerically-machine-verified core (3-branch fast-Ca "
                      "bifurcation structure + the claim's disclosed waveform negative). The "
                      "glucose-dose secretion curve, sulfonylurea/mutation predictions, and "
                      "cross-scale insulin-pulse anchors are qualitative/out-of-scope per the "
                      "claim's honest_gaps and are NOT gated here.")

os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
with open(RESULTS_PATH, "w") as fh:
    json.dump(OUT, fh, indent=2,
               default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

print(json.dumps(OUT["OVERALL"], indent=2))
for k in ["G1_three_branch_structure", "G2_window_a0.001", "G3_window_a0.003",
          "G4_disclosed_negative_check", "VOID_FLOOR"]:
    print(k, OUT[k])
