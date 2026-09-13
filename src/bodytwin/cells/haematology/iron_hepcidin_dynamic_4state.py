"""Iron/hepcidin 4-state dynamical ODE (hours-to-decades timescale) -- untuned rebuild.

Reimplements the 4-state iron/hepcidin model from its claim-text specification (states, ODEs,
named rate constants) and checks whether the claimed headline numbers reproduce under an untuned
build. No parameter is tuned to hit a headline number except GAIN_IL6, which the specification
requires to be root-found once against the Nemeth 2004 calibration target (a calibration step, not
a fit to outcome) and is then frozen for every other check (Kemna, HFE, nuisance).

This is the fast dynamical (hours) companion to the iron_hepcidin cell, which is a slow
quasi-algebraic/mass-balance treatment of TSAT/recycling and tests different quantities.

Disclosed conventions (not given verbatim by the specification):
(A) Fe_p -> TSAT: TSAT(t) = 100 * Fe_p(t) / (TIBC_CONC_MG_PER_L * PLASMA_VOL_L), with
    TIBC_CONC_MG_PER_L = 3.0 mg/L (= 300 ug/dL, standard total-iron-binding-capacity reference)
    and PLASMA_VOL_L = BLOOD_VOLUME_L*(1-HCT) = 5.0*0.55 = 2.75 L. Hence Fe_p0 = 0.30*3.0*2.75 =
    2.475 mg at TSAT baseline 30% -- the right order of magnitude for circulating
    transferrin-bound iron (~2-4 mg), not independently literature-pinned.
(B) Fe_s0: the model has no separate Hb-bound-iron state, so Fe_s is "everything except the small
    circulating transferrin pool": Fe_s0 = 3500 - 2.475 = 3497.525 mg. Coarse lumping (Hb iron,
    which does not participate in these fast dynamics, is folded into storage for bookkeeping).
(C) IL6(t) input shapes (illustrative, not literature-digitized): Nemeth 2004 (continuous IV
    infusion) -> step-like rise IL6(t) = 1 - exp(-t/0.25 h); Kemna 2005 (single LPS bolus) ->
    Weibull-shaped pulse IL6(t) = (t/tp)*exp(1 - t/tp) with tp = 3.0 h, matching Kemna 2005's
    reported "IL-6 peaked within 3 h". Both dimensionless, amplitude-normalized to 1; GAIN_IL6
    carries the calibrated scale.
(D) HFE hemochromatosis is modeled as blunted hepcidin PRODUCTION capacity (k_synH divided by the
    Bridle 2003 measured 5.4-fold gain attenuation), not a change in FPN internalization kinetics.

Reads: nothing. Writes: iron_hepcidin_dynamic_4state_results.json under the cell output directory.
Gates (pre-registered, tolerances fixed before running, see the GATES dict): exact baseline fixed
point and four stable Jacobian eigenvalues; Nemeth calibration converged; emergent serum-iron
change same sign and within 2x of the measured -34%; Kemna held-out hepcidin lag positive and
within 3x of the measured 3 h; HFE TSAT crosses the EASL 45% threshold by 5 y; and the HFE
trajectory beats a +/-30% dietary/loss nuisance sweep at intact gain.
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "iron_hepcidin_dynamic_4state")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================================
# FIXED, GIVEN CONSTANTS (from the model's claim text -- not tuned)
# ============================================================================================
TOTAL_BODY_IRON_NORMAL_MG = 3500.0          # 3.5 g, band 3-4 g
TSAT_BASELINE = 0.30
BLOOD_VOLUME_L = 5.0
HCT = 0.45
TIBC_CONC_MG_PER_L = 3.0                    # disclosed convention (A): 300 ug/dL

RECYCLING_MG_DAY = 22.6
ABSORPTION_MG_DAY = 1.5
OBLIGATE_LOSS_MG_DAY = 1.5

K_DEGH = 0.0685                             # /h, Enculescu2017 best-fit
H_NORMAL_NM = 23.0                          # nM, Parmar&Mendes2019 WT steady state
FPN0 = 0.5
FPN_FAST_RELAX_RATE = 3.1                   # /h
BETA_IL6_FPN = 0.6
NTBI_FRACTION_OF_OBLIGATE_LOSS = 0.5        # disclosed convention, shared with the iron_hepcidin cell

HFE_GAIN_ATTENUATION = 5.4                  # Bridle2003

# ============================================================================================
# DERIVED SCALES (disclosed conventions A, B)
# ============================================================================================
PLASMA_VOL_L = BLOOD_VOLUME_L * (1.0 - HCT)                          # 2.75 L
Fe_p0 = TSAT_BASELINE * TIBC_CONC_MG_PER_L * PLASMA_VOL_L            # 2.475 mg
Fe_s0 = TOTAL_BODY_IRON_NORMAL_MG - Fe_p0                            # 3497.525 mg


def tsat_pct(fe_p):
    return 100.0 * fe_p / (TIBC_CONC_MG_PER_L * PLASMA_VOL_L)


# ============================================================================================
# RATE CONSTANTS SOLVED FROM THE BASELINE FIXED-POINT CONDITION (Fe_p0, H0, FPN0, Fe_s0),
# using ONLY the given fluxes/rates -- no free parameters left to tune.
# ============================================================================================
recyc_flux_mg_h = RECYCLING_MG_DAY / 24.0        # FPN0*v_recyc at baseline
abs_flux_mg_h = ABSORPTION_MG_DAY / 24.0          # FPN0*v_abs at baseline
obligate_loss_mg_h = OBLIGATE_LOSS_MG_DAY / 24.0  # k_trueloss*Fe_p0 at baseline

v_recyc = recyc_flux_mg_h / FPN0
v_abs = abs_flux_mg_h / FPN0
k_trueloss = obligate_loss_mg_h / Fe_p0

# remaining Fe_p sink (k_recycreturn + k_ntbi)*Fe_p0 must absorb the recycling+abs inflow
# exactly (dFe_p/dt=0 at baseline), with k_ntbi fixed as a disclosed fraction of k_trueloss
remaining_sink_mg_h = recyc_flux_mg_h + abs_flux_mg_h - obligate_loss_mg_h
k_ntbi = (NTBI_FRACTION_OF_OBLIGATE_LOSS * k_trueloss * Fe_p0) / Fe_p0
k_recycreturn = (remaining_sink_mg_h - k_ntbi * Fe_p0) / Fe_p0

# hepcidin synthesis constant from dH/dt=0 at baseline (IL6=0)
k_synH = K_DEGH * H_NORMAL_NM

# FPN kinetics from dFPN/dt=0 at baseline AND the given fast-relaxation rate as the
# linearized diagonal |d(dFPN/dt)/dFPN| at the fixed point (IL6=0):
#   k_synFPN + k_intFPN*H0 = FPN_FAST_RELAX_RATE   (relaxation)
#   k_synFPN*(1-FPN0) = k_intFPN*H0*FPN0            (steady state, FPN0=0.5 -> k_synFPN=H0*k_intFPN)
k_intFPN = FPN_FAST_RELAX_RATE / (2.0 * H_NORMAL_NM)
k_synFPN = H_NORMAL_NM * k_intFPN

# sanity: exact fixed point at (Fe_p0, H_NORMAL_NM, FPN0, Fe_s0)
Y0 = np.array([Fe_p0, H_NORMAL_NM, FPN0, Fe_s0])


def rhs(t, y, gain_il6, il6_func, k_synH_local, hdrug_func=None, antag_func=None):
    fe_p, h, fpn, fe_s = y
    il6 = il6_func(t)
    h_drug = hdrug_func(t) if hdrug_func is not None else 0.0
    antag = antag_func(t) if antag_func is not None else 0.0
    h_eff = (h + h_drug) * (1.0 - antag)

    dH = k_synH_local * (fe_p / Fe_p0) * np.exp(gain_il6 * il6) - K_DEGH * h
    dFPN = k_synFPN * (1.0 - fpn) * (1.0 - BETA_IL6_FPN * il6) - k_intFPN * h_eff * fpn
    dFe_p = (fpn * v_recyc * (fe_s / Fe_s0) + fpn * v_abs
             - k_trueloss * fe_p - k_recycreturn * fe_p - k_ntbi * fe_p)
    dFe_s = k_recycreturn * fe_p - fpn * v_recyc * (fe_s / Fe_s0) + k_ntbi * fe_p
    return [dFe_p, dH, dFPN, dFe_s]


def il6_zero(t):
    return 0.0


def il6_nemeth_step(t):
    return 1.0 - np.exp(-t / 0.25)


def il6_kemna_pulse(t, tp=3.0):
    x = t / tp
    return x * np.exp(1.0 - x)


# ============================================================================================
# CHECK 1: fixed point + Jacobian eigenvalues at baseline (IL6=0, gain irrelevant since IL6=0)
# ============================================================================================
f0 = rhs(0.0, Y0, gain_il6=1.0, il6_func=il6_zero, k_synH_local=k_synH)
max_residual = float(np.max(np.abs(f0)))

eps = 1e-6
J = np.zeros((4, 4))
for j in range(4):
    dy = np.zeros(4)
    dy[j] = eps * max(abs(Y0[j]), 1.0)
    f_plus = np.array(rhs(0.0, Y0 + dy, gain_il6=1.0, il6_func=il6_zero, k_synH_local=k_synH))
    f_minus = np.array(rhs(0.0, Y0 - dy, gain_il6=1.0, il6_func=il6_zero, k_synH_local=k_synH))
    J[:, j] = (f_plus - f_minus) / (2 * dy[j])

eigvals = np.linalg.eigvals(J)
eigvals_sorted = np.array(sorted(eigvals.real, key=lambda x: -abs(x)))
all_stable = bool(np.all(eigvals.real < 0))

CLAIMED_EIGS = [-3.105, -0.354, -0.118, -0.00007]

# ============================================================================================
# CHECK 2/3: Nemeth2004 -- root-find GAIN_IL6 so hepcidin is 7.5x higher at t=2h; report
# emergent (non-fit) serum-iron % change at the SAME timepoint, same run.
# ============================================================================================
NEMETH_T = 2.0
NEMETH_FOLD_TARGET = 7.5


def nemeth_run(gain):
    sol = solve_ivp(rhs, [0, NEMETH_T], Y0, method="Radau",
                     args=(gain, il6_nemeth_step, k_synH),
                     t_eval=[NEMETH_T], rtol=1e-9, atol=1e-12)
    fe_p_t, h_t, fpn_t, fe_s_t = sol.y[:, -1]
    return h_t, fe_p_t


def nemeth_fold_residual(gain):
    h_t, _ = nemeth_run(gain)
    return (h_t / H_NORMAL_NM) - NEMETH_FOLD_TARGET


GAIN_IL6 = brentq(nemeth_fold_residual, 0.5, 20.0, xtol=1e-8, rtol=1e-10)
h_2h, fe_p_2h = nemeth_run(GAIN_IL6)
nemeth_fold_achieved = h_2h / H_NORMAL_NM
nemeth_serum_iron_pctchange = (fe_p_2h / Fe_p0 - 1.0) * 100.0

CLAIMED_GAIN_IL6 = 4.157
CLAIMED_NEMETH_FOLD = 7.442
CLAIMED_SERUM_IRON_PCTCHANGE = -42.46
REAL_SERUM_IRON_PCTCHANGE = -34.0

# ============================================================================================
# CHECK 4: Kemna2005 HELD-OUT (same GAIN_IL6/K_DEGH, no re-fit; different IL6 pulse shape)
# ============================================================================================
KEMNA_TP = 3.0
t_grid = np.linspace(0, 24, 4801)
sol_kemna = solve_ivp(rhs, [0, 24], Y0, method="Radau",
                       args=(GAIN_IL6, lambda t: il6_kemna_pulse(t, KEMNA_TP), k_synH),
                       t_eval=t_grid, rtol=1e-9, atol=1e-12, dense_output=True)
h_traj = sol_kemna.y[1, :]
peak_idx = int(np.argmax(h_traj))
# refine peak with a finer local grid around the coarse max
lo = max(0.0, t_grid[peak_idx] - 0.5)
hi = min(24.0, t_grid[peak_idx] + 0.5)
t_fine = np.linspace(lo, hi, 2001)
h_fine = sol_kemna.sol(t_fine)[1, :]
kemna_peak_t = float(t_fine[int(np.argmax(h_fine))])
kemna_lag_h = kemna_peak_t - KEMNA_TP

CLAIMED_KEMNA_LAG = 2.8
REAL_KEMNA_LAG = 3.0  # Kemna2005: hepcidin peak ~6h, IL-6 peak within 3h

# ============================================================================================
# CHECK 5: HFE hemochromatosis -- chronic 5.4x hepcidin-gain attenuation, 25y trajectory
# ============================================================================================
YEAR_H = 365.0 * 24.0
k_synH_HFE = k_synH / HFE_GAIN_ATTENUATION
t_years = [1, 5, 10, 25]
t_eval_hours = [y * YEAR_H for y in t_years]

sol_hfe = solve_ivp(rhs, [0, 25 * YEAR_H], Y0, method="Radau",
                     args=(1.0, il6_zero, k_synH_HFE),
                     t_eval=t_eval_hours, rtol=1e-8, atol=1e-10)
hfe_tsat = {y: float(tsat_pct(sol_hfe.y[0, i])) for i, y in enumerate(t_years)}

CLAIMED_HFE_TSAT = {1: 58.5, 5: 87.0, 10: 107.4, 25: 136.3}

# ============================================================================================
# CHECK 6: nuisance-sweep adversary -- +/-30% dietary/loss wobble at INTACT (WT) gain,
# N draws, worst-case TSAT at 5y vs the HFE TSAT@5y computed above.
# ============================================================================================
rng = np.random.default_rng(20260728)
N_NUISANCE = 30
nuisance_tsat_5y = np.zeros(N_NUISANCE)
for i in range(N_NUISANCE):
    d_abs = rng.uniform(-0.30, 0.30)
    d_loss = rng.uniform(-0.30, 0.30)
    v_abs_pert = v_abs * (1.0 + d_abs)
    k_trueloss_pert = k_trueloss * (1.0 + d_loss)

    def rhs_pert(t, y, gain_il6, il6_func, k_synH_local,
                 v_abs_p=v_abs_pert, k_trueloss_p=k_trueloss_pert):
        fe_p, h, fpn, fe_s = y
        il6 = il6_func(t)
        h_eff = h
        dH = k_synH_local * (fe_p / Fe_p0) * np.exp(gain_il6 * il6) - K_DEGH * h
        dFPN = k_synFPN * (1.0 - fpn) * (1.0 - BETA_IL6_FPN * il6) - k_intFPN * h_eff * fpn
        dFe_p = (fpn * v_recyc * (fe_s / Fe_s0) + fpn * v_abs_p
                 - k_trueloss_p * fe_p - k_recycreturn * fe_p - k_ntbi * fe_p)
        dFe_s = k_recycreturn * fe_p - fpn * v_recyc * (fe_s / Fe_s0) + k_ntbi * fe_p
        return [dFe_p, dH, dFPN, dFe_s]

    sol_n = solve_ivp(rhs_pert, [0, 5 * YEAR_H], Y0, method="Radau",
                       args=(1.0, il6_zero, k_synH),
                       t_eval=[5 * YEAR_H], rtol=1e-8, atol=1e-10)
    nuisance_tsat_5y[i] = tsat_pct(sol_n.y[0, -1])

nuisance_worst_case = float(np.max(nuisance_tsat_5y))

CLAIMED_NUISANCE_CEILING = 40.0

# ============================================================================================
# GATES -- pre-registered BEFORE running (tolerances fixed above in comments), machine
# PASS/FAIL. Reproduces the pre-registered C/NOT-C thresholds where stated.
# ============================================================================================
GATES = {
    "gate1_fixed_point_exact": bool(max_residual < 1e-8),
    "gate1_all_4_eigs_stable": all_stable,
    "gate2_nemeth_calibration_converged": bool(abs(nemeth_fold_achieved - NEMETH_FOLD_TARGET) < 0.01),
    "gate3_serum_iron_emergent_same_sign_and_within_2x_of_real": bool(
        nemeth_serum_iron_pctchange < 0
        and abs(nemeth_serum_iron_pctchange) < 2.0 * abs(REAL_SERUM_IRON_PCTCHANGE)
    ),
    "gate4_kemna_lag_within_3x_real_per_cells_own_threshold": bool(
        0.0 < kemna_lag_h < 3.0 * REAL_KEMNA_LAG
    ),
    "gate5_hfe_crosses_easl_45pct_by_5y": bool(hfe_tsat[5] > 45.0),
    "gate6_hfe_beats_nuisance_specificity": bool(nuisance_worst_case < hfe_tsat[5]),
}
required_gates_overall_pass = bool(all(GATES.values()))

# ============================================================================================
# OUTPUT
# ============================================================================================
results = {
    "disclosed_conventions": {
        "tsat_conversion": "TSAT=100*Fe_p/(TIBC_CONC_MG_PER_L*PLASMA_VOL_L), "
                            "TIBC_CONC_MG_PER_L=3.0, PLASMA_VOL_L=2.75",
        "Fe_p0_mg": round(Fe_p0, 4),
        "Fe_s0_mg": round(Fe_s0, 3),
        "il6_shapes": "Nemeth: 1-exp(-t/0.25) step; Kemna: (t/3)*exp(1-t/3) Weibull pulse, peak@3h",
        "hfe_mechanism": "k_synH divided by 5.4 (blunted hepcidin PRODUCTION), FPN kinetics unchanged",
    },
    "solved_rate_constants_1_per_h_unless_noted": {
        "k_synH_nM_per_h": round(k_synH, 5), "k_degH": K_DEGH,
        "k_synFPN": round(k_synFPN, 5), "k_intFPN_per_nM_per_h": round(k_intFPN, 6),
        "v_recyc_mg_per_h": round(v_recyc, 5), "v_abs_mg_per_h": round(v_abs, 6),
        "k_trueloss": round(k_trueloss, 6), "k_recycreturn": round(k_recycreturn, 6),
        "k_ntbi": round(k_ntbi, 6),
    },
    "check1_fixed_point": {
        "max_residual": max_residual,
        "eigenvalues_1_per_h_computed": eigvals_sorted.tolist(),
        "eigenvalues_claimed": CLAIMED_EIGS,
        "all_4_stable_computed": all_stable,
    },
    "check2_nemeth_calibration": {
        "GAIN_IL6_root_found": round(GAIN_IL6, 5),
        "GAIN_IL6_claimed": CLAIMED_GAIN_IL6,
        "fold_change_2h_computed": round(nemeth_fold_achieved, 4),
        "fold_change_2h_claimed": CLAIMED_NEMETH_FOLD,
        "fold_change_2h_target": NEMETH_FOLD_TARGET,
    },
    "check3_nemeth_serum_iron_emergent": {
        "pctchange_2h_computed": round(nemeth_serum_iron_pctchange, 3),
        "pctchange_2h_claimed": CLAIMED_SERUM_IRON_PCTCHANGE,
        "pctchange_2h_real_nemeth2004": REAL_SERUM_IRON_PCTCHANGE,
    },
    "check4_kemna_holdout": {
        "il6_pulse_peak_h": KEMNA_TP,
        "hepcidin_peak_h_computed": round(kemna_peak_t, 3),
        "lag_h_computed": round(kemna_lag_h, 3),
        "lag_h_claimed": CLAIMED_KEMNA_LAG,
        "lag_h_real_kemna2005": REAL_KEMNA_LAG,
    },
    "check5_hfe_trajectory": {
        "tsat_pct_computed_by_year": {str(y): round(hfe_tsat[y], 2) for y in t_years},
        "tsat_pct_claimed_by_year": CLAIMED_HFE_TSAT,
    },
    "check6_nuisance_sweep": {
        "n_draws": N_NUISANCE,
        "worst_case_tsat_5y_computed": round(nuisance_worst_case, 3),
        "hfe_tsat_5y_computed": round(hfe_tsat[5], 3),
        "worst_case_ceiling_claimed": CLAIMED_NUISANCE_CEILING,
    },
    "gates": GATES,
    "required_gates_overall_pass": required_gates_overall_pass,
}

out_path = _os.path.join(OUT_DIR, "iron_hepcidin_dynamic_4state_results.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else o)

print(f"wrote {out_path}")
print(json.dumps(GATES, indent=2))
print(f"\nrequired_gates_overall_pass: {required_gates_overall_pass}")
print(json.dumps(results, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else o))
