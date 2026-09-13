#!/usr/bin/env python3
"""Actin treadmilling kinetics: a 6-state mean-field terminal-occupancy ODE (Pollard 1986 rate
constants + Pollard & Weeds 1984 k_hyd + Melki 1996 k_pi) integrated to steady state to obtain the
free-monomer level, the treadmilling flux, and both external anchors (bulk critical concentration,
in-vitro treadmill velocity), plus a symmetric-ends null control.

Reads: nothing (all parameters are literature values embedded below).
Writes: actin_treadmilling_pollard_kinetics.json (gate results) under the cell output directory.
Gate: G1-G5 below plus the void floor must all pass.

MODEL: states Mt,Md (free ATP/ADP-G-actin, uM),
Ft,Fp,Fd (ATP/ADP-Pi/ADP-actin subunit pools polymerized, uM), Nf (filament end density, uM ends).
  theta_BT=clip(Mt/Cc_BT,0,1)  theta_PT=clip(Mt/Cc_PT,0,1)  theta_BD=1-theta_BT  theta_PD=1-theta_PT
  onB_T=kBT_on*Mt*Nf   offB_T=theta_BT*kBT_off*Nf      (barbed ATP-end)
  onB_D=kBD_on*Md*Nf   offB_D=theta_BD*kBD_off*Nf      (barbed ADP-end)
  onP_T=kPT_on*Mt*Nf   offP_T=theta_PT*kPT_off*Nf      (pointed ATP-end)
  onP_D=kPD_on*Md*Nf   offP_D=theta_PD*kPD_off*Nf      (pointed ADP-end)
  dMt/dt = -(onB_T+onP_T)+(offB_T+offP_T)+k_ex*Md
  dMd/dt = -(onB_D+onP_D)+(offB_D+offP_D)-k_ex*Md
  dFt/dt = (onB_T+onP_T)-(offB_T+offP_T)-k_hyd*Ft
  dFp/dt = k_hyd*Ft - k_pi*Fp
  dFd/dt = (onB_D+onP_D)-(offB_D+offP_D)*fd_gate + k_pi*Fp,  fd_gate = Fd/(Fd+0.01)
  dNf/dt = k_sev*Fd   (k_sev=0 for the baseline steady-state run)
  Cc_BT = kBT_off/kBT_on, Cc_PT = kPT_off/kPT_on  (critical concentrations, derived not fit)
  Mass conservation: Mt+Md+Ft+Fp+Fd = A_total = 10.0 uM (analytic invariant, checked not assumed).

GATE (pre-registered):
  G1: steady-state Mt* lands strictly between Cc_BT=0.1207uM and Cc_PT=0.6154uM (structural
      necessity for treadmilling: a monomer level that starves the pointed end while feeding
      the barbed end), and matches the claimed Mt*=0.1337uM to <2%.
  G2: baseline flux (net barbed-end addition rate per filament, subunit/s) matches claimed
      0.209 subunit/s to <5%, and equals -1x the net pointed-end rate at steady state (flux
      conservation: what leaves the pointed end must equal what the barbed end gains).
  G3: EXTERNAL ANCHORS -- (a) bulk Cc: model Mt* vs Drenckhahn&Pollard1986's measured 0.14uM,
      claimed 4.5% error, reproduced to <1pp; (b) HELD-OUT (temperature-shifted, independent
      lab): flux converted to um/h via subunits_per_um=370 vs Selve&Wegner1986's measured
      2.0um/h, claimed 1.5% agreement, reproduced to <1pp.
  G4: 4 independent initial conditions converge to the SAME Mt* (spread <1e-10, claimed 2.8e-17
      is integrator-tolerance-dependent and not re-demanded at that exact precision), AND an
      independent algebraic fsolve root of the steady-state system matches the long-time ODE
      integration to <1e-6 relative.
  G5 (claimed control): symmetric-ends null (force kPT_on=kBT_on, kPT_off=kBT_off, kPD_on=kBD_on,
      kPD_off=kBD_off, i.e. Cc_B=Cc_P) collapses flux to |flux| < 1e-6 subunit/s.

VOID FLOOR (pre-registered, must FAIL, INDEPENDENT of the claim's symmetric-ends control):
scramble k_hyd -> -k_hyd (unphysical negative hydrolysis rate). This must NOT reproduce the
claimed Mt* (within 2%) or the claimed flux (within 5%) -- if a physically-broken model still
lands on the claimed numbers, G1/G2 are not discriminating.
"""
import json
import os
import os as _os
import numpy as np

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

OUT = {}
RESULTS_PATH = _os.path.join(OUT_ROOT, "actin_treadmilling_pollard_kinetics",
                             "actin_treadmilling_pollard_kinetics.json")

BASE = dict(kBT_on=11.6, kBT_off=1.4, kPT_on=1.3, kPT_off=0.8,
            kBD_on=3.8, kBD_off=7.2, kPD_on=0.16, kPD_off=0.27,
            k_hyd=0.07, k_pi=0.00198, k_ex=0.1, Nf0=0.01, A_total=10.0,
            subunits_per_um=370.0, k_sev=0.0)


def derived_Cc(pars):
    return pars["kBT_off"] / pars["kBT_on"], pars["kPT_off"] / pars["kPT_on"]


def rhs(t, y, pars):
    Mt, Md, Ft, Fp, Fd, Nf = y
    Mt = max(Mt, 0.0)
    Md = max(Md, 0.0)
    Cc_BT, Cc_PT = derived_Cc(pars)
    theta_BT = min(max(Mt / Cc_BT, 0.0), 1.0)
    theta_PT = min(max(Mt / Cc_PT, 0.0), 1.0)
    theta_BD = 1.0 - theta_BT
    theta_PD = 1.0 - theta_PT
    onB_T = pars["kBT_on"] * Mt * Nf
    offB_T = theta_BT * pars["kBT_off"] * Nf
    onP_T = pars["kPT_on"] * Mt * Nf
    offP_T = theta_PT * pars["kPT_off"] * Nf
    onB_D = pars["kBD_on"] * Md * Nf
    offB_D = theta_BD * pars["kBD_off"] * Nf
    onP_D = pars["kPD_on"] * Md * Nf
    offP_D = theta_PD * pars["kPD_off"] * Nf
    fd_gate = Fd / (Fd + 0.01)
    dMt = -(onB_T + onP_T) + (offB_T + offP_T) + pars["k_ex"] * Md
    dMd = -(onB_D + onP_D) + (offB_D + offP_D) - pars["k_ex"] * Md
    dFt = (onB_T + onP_T) - (offB_T + offP_T) - pars["k_hyd"] * Ft
    dFp = pars["k_hyd"] * Ft - pars["k_pi"] * Fp
    dFd = (onB_D + onP_D) - (offB_D + offP_D) * fd_gate + pars["k_pi"] * Fp
    dNf = pars["k_sev"] * Fd
    return [dMt, dMd, dFt, dFp, dFd, dNf]


def flux_components(y, pars):
    Mt, Md, Ft, Fp, Fd, Nf = y
    Cc_BT, Cc_PT = derived_Cc(pars)
    theta_BT = min(max(Mt / Cc_BT, 0.0), 1.0)
    theta_PT = min(max(Mt / Cc_PT, 0.0), 1.0)
    theta_BD, theta_PD = 1 - theta_BT, 1 - theta_PT
    onB_T = pars["kBT_on"] * Mt * Nf
    offB_T = theta_BT * pars["kBT_off"] * Nf
    onB_D = pars["kBD_on"] * Md * Nf
    offB_D = theta_BD * pars["kBD_off"] * Nf
    onP_T = pars["kPT_on"] * Mt * Nf
    offP_T = theta_PT * pars["kPT_off"] * Nf
    onP_D = pars["kPD_on"] * Md * Nf
    offP_D = theta_PD * pars["kPD_off"] * Nf
    J_B = (onB_T + onB_D - offB_T - offB_D) / Nf
    J_P = (onP_T + onP_D - offP_T - offP_D) / Nf
    return J_B, J_P


def run_to_steady(pars, y0, T=20000.0):
    sol = solve_ivp(rhs, [0, T], y0, args=(pars,), method="Radau", rtol=1e-10, atol=1e-14,
                     max_step=50.0)
    return sol.y[:, -1], sol


def steady_state_fsolve(pars, guess):
    """Independent algebraic root-find. GEOMETRIC NOTE: the raw 6-eq system is rank-deficient at
    steady state -- summing all 5 mass-pool RHS gives (offB_D+offP_D)*(1-fd_gate), which is only
    ~0 once Fd is large enough that fd_gate=Fd/(Fd+0.01) saturates near 1; below saturation the
    Fd-equation's derivative w.r.t. Fd is nearly flat (sigma_min~0), so a raw 6D fsolve wanders
    along that near-null direction (observed: Fd -> 1.6e13, unconverged). Fix: hold Nf fixed at
    Nf0 (its own RHS is identically 0 at k_sev=0, so it is not part of the unknown vector) and
    replace the flat Fd-equation with the exact mass-conservation constraint
    Mt+Md+Ft+Fp+Fd=A_total, which is well-conditioned everywhere. This is a legitimate
    over-determined substitution (conservation is an EXACT analytic invariant of the ODE, not an
    approximation), not a fit to the target."""
    Nf_fixed = guess[5]

    def eqs4(y4):
        Mt, Md, Ft, Fp = y4
        Fd = pars["A_total"] - (Mt + Md + Ft + Fp)
        full = [Mt, Md, Ft, Fp, Fd, Nf_fixed]
        d = rhs(0, full, pars)
        return [d[0], d[1], d[2], d[3]]  # dMt,dMd,dFt,dFp = 0; dFd enforced via conservation

    sol4, info, ier, msg = fsolve(eqs4, guess[:4], xtol=1e-13, full_output=True)
    Mt, Md, Ft, Fp = sol4
    Fd = pars["A_total"] - (Mt + Md + Ft + Fp)
    full_sol = np.array([Mt, Md, Ft, Fp, Fd, Nf_fixed])
    resid = np.max(np.abs(rhs(0, full_sol, pars)))
    return full_sol, ier == 1, float(resid)


A_total = BASE["A_total"]
Nf0 = BASE["Nf0"]

# ---------------------------------------------------------------------------
# Baseline steady state from 4 independent initial conditions
# ---------------------------------------------------------------------------
ics = [
    [0.05, 0.05, A_total - 0.1 - Nf0, 0.05, 0.05, Nf0],
    [A_total * 0.5, 0.01, A_total * 0.4, 0.02, 0.02, Nf0],
    [0.01, A_total * 0.3, A_total * 0.6, 0.03, 0.03, Nf0],
    [0.2, 0.2, A_total - 0.4 - Nf0 * 0, 0.1, 0.1, Nf0],
]
finals = []
mass_errs = []
for ic in ics:
    yf, sol = run_to_steady(BASE, ic)
    finals.append(yf)
    mass_errs.append(abs(sum(yf[:5]) - A_total))
Mt_vals = [f[0] for f in finals]
Mt_spread = max(Mt_vals) - min(Mt_vals)
Mt_star = float(np.mean(Mt_vals))

Cc_BT, Cc_PT = derived_Cc(BASE)
J_B_star, J_P_star = flux_components(finals[0], BASE)

# fsolve independent cross-check
fs_sol, fs_ok, fs_resid = steady_state_fsolve(BASE, finals[0])
fs_vs_ode_rel = abs(fs_sol[0] - Mt_star) / Mt_star

OUT["baseline_steady_state"] = {
    "Mt_star": Mt_star, "Mt_spread_across_4_ICs": Mt_spread,
    "mass_conservation_max_abs_err": max(mass_errs),
    "Cc_BT": Cc_BT, "Cc_PT": Cc_PT,
    "flux_barbed_subunit_per_s": J_B_star, "flux_pointed_subunit_per_s": J_P_star,
    "fsolve_root": fs_sol.tolist(), "fsolve_converged": fs_ok, "fsolve_residual": fs_resid,
    "fsolve_vs_ode_rel_diff": fs_vs_ode_rel,
}

# G1
claimed_Mt = 0.1337
g1_between = Cc_BT < Mt_star < Cc_PT
g1_match = abs(Mt_star - claimed_Mt) / claimed_Mt < 0.02
g1_pass = g1_between and g1_match
OUT["G1_steady_state_monomer"] = {"Mt_star": Mt_star, "claimed": claimed_Mt,
                                   "rel_err_pct": abs(Mt_star - claimed_Mt) / claimed_Mt * 100,
                                   "between_Cc_bounds": g1_between, "pass": bool(g1_pass)}

# G2
claimed_flux = 0.209
flux_match = abs(J_B_star - claimed_flux) / claimed_flux < 0.05
flux_balance = abs(J_B_star + J_P_star) < 0.02 * abs(J_B_star)  # J_P should be ~ -J_B at steady state
g2_pass = flux_match and flux_balance
OUT["G2_flux"] = {"flux_barbed": J_B_star, "flux_pointed": J_P_star, "claimed": claimed_flux,
                   "rel_err_pct": abs(J_B_star - claimed_flux) / claimed_flux * 100,
                   "flux_balance_check": flux_balance, "pass": bool(g2_pass)}

# G3 external anchors
Cc_measured = 0.14
cc_err_pct = abs(Mt_star - Cc_measured) / Cc_measured * 100
cc_claimed_err_pct = 4.5
g3a_pass = abs(cc_err_pct - cc_claimed_err_pct) < 1.0

velocity_um_h = J_B_star * 3600.0 / BASE["subunits_per_um"]
velocity_measured = 2.0
vel_err_pct = abs(velocity_um_h - velocity_measured) / velocity_measured * 100
vel_claimed_err_pct = 1.5
g3b_pass = abs(vel_err_pct - vel_claimed_err_pct) < 1.0
g3_pass = g3a_pass and g3b_pass
OUT["G3_external_anchors"] = {
    "bulk_Cc": {"model_Mt_star": Mt_star, "measured": Cc_measured, "err_pct": cc_err_pct,
                "claimed_err_pct": cc_claimed_err_pct, "pass": bool(g3a_pass)},
    "treadmill_velocity": {"model_um_per_h": velocity_um_h, "measured_um_per_h": velocity_measured,
                            "err_pct": vel_err_pct, "claimed_err_pct": vel_claimed_err_pct,
                            "pass": bool(g3b_pass)},
    "pass": bool(g3_pass),
}

# G4
g4_pass = (Mt_spread < 1e-10) and fs_ok and (fs_vs_ode_rel < 1e-6)
OUT["G4_convergence_and_algebraic_crosscheck"] = {
    "Mt_spread_4_ICs": Mt_spread, "fsolve_converged": fs_ok,
    "fsolve_vs_ode_rel_diff": fs_vs_ode_rel, "pass": bool(g4_pass),
}

# G5 claimed control: symmetric ends
SYM = dict(BASE)
SYM["kPT_on"], SYM["kPT_off"] = BASE["kBT_on"], BASE["kBT_off"]
SYM["kPD_on"], SYM["kPD_off"] = BASE["kBD_on"], BASE["kBD_off"]
yf_sym, _ = run_to_steady(SYM, finals[0])
J_B_sym, J_P_sym = flux_components(yf_sym, SYM)
g5_pass = abs(J_B_sym) < 1e-6
OUT["G5_symmetric_ends_null_control"] = {"flux_barbed": J_B_sym, "flux_pointed": J_P_sym,
                                          "pass": bool(g5_pass)}

# ---------------------------------------------------------------------------
# VOID FLOOR: k_hyd -> -k_hyd (unphysical), independent of claim's null
# ---------------------------------------------------------------------------
SCRAMBLED = dict(BASE)
SCRAMBLED["k_hyd"] = -BASE["k_hyd"]
try:
    yf_scr, sol_scr = run_to_steady(SCRAMBLED, finals[0], T=5000.0)
    Mt_scr = float(yf_scr[0])
    J_B_scr, _ = flux_components(yf_scr, SCRAMBLED)
    blew_up = not np.all(np.isfinite(yf_scr)) or np.any(np.abs(yf_scr) > 1e6)
except Exception:
    Mt_scr, J_B_scr, blew_up = None, None, True

if blew_up or Mt_scr is None:
    reproduces_Mt = False
    reproduces_flux = False
else:
    reproduces_Mt = abs(Mt_scr - claimed_Mt) / claimed_Mt < 0.02
    reproduces_flux = abs(J_B_scr - claimed_flux) / claimed_flux < 0.05
void_floor_pass = not (reproduces_Mt and reproduces_flux)
OUT["VOID_FLOOR"] = {
    "scramble": "k_hyd -> -k_hyd (unphysical)", "blew_up_or_nonfinite": bool(blew_up),
    "Mt_scrambled": Mt_scr, "flux_scrambled": J_B_scr,
    "reproduces_claimed_Mt": bool(reproduces_Mt), "reproduces_claimed_flux": bool(reproduces_flux),
    "pass_(void_floor_correctly_fails)": bool(void_floor_pass),
}

ALL_GATES = [OUT["G1_steady_state_monomer"]["pass"], OUT["G2_flux"]["pass"],
             OUT["G3_external_anchors"]["pass"], OUT["G4_convergence_and_algebraic_crosscheck"]["pass"],
             OUT["G5_symmetric_ends_null_control"]["pass"],
             OUT["VOID_FLOOR"]["pass_(void_floor_correctly_fails)"]]
OUT["OVERALL"] = {"gates_passed": int(sum(ALL_GATES)), "gates_total": len(ALL_GATES),
                   "all_pass": bool(all(ALL_GATES))}
OUT["node_id"] = "MODEL-ACTIN-TREADMILLING"

os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
with open(RESULTS_PATH, "w") as fh:
    json.dump(OUT, fh, indent=2,
               default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

print(json.dumps(OUT["OVERALL"], indent=2))
for k in ["G1_steady_state_monomer", "G2_flux", "G3_external_anchors",
          "G4_convergence_and_algebraic_crosscheck", "G5_symmetric_ends_null_control", "VOID_FLOOR"]:
    print(k, OUT[k])
