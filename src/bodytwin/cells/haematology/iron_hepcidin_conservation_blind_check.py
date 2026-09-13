"""Iron-hepcidin homeostasis -- independent rebuild of the 4-state ODE plus a conservation falsifier.

Rebuilt only from the restated equations and constants of the model, without reading the first
rebuild (the iron_hepcidin_dynamic_4state cell).

Claim under test: summing the two iron-mass equations (dFe_p/dt, dFe_s/dt) shows total iron leaves
the system by exactly one route (k_trueloss*Fe_p), which forces a hard fixed-point ceiling on
plasma iron / TSAT, independent of integration or tuning -- contradicting the recorded long-horizon
HFE prediction of unbounded growth to 107.4% (10 y) / 136.3% (25 y).

Equations, y = [Fe_p, H, FPN, Fe_s]:
  dH/dt    = k_synH*(Fe_p/Fe_p0)*exp(GAIN_IL6*IL6(t)) - k_degH*H
  dFPN/dt  = k_synFPN*(1-FPN)*(1-beta_IL6*IL6(t)) - k_intFPN*H_eff*FPN
  dFe_p/dt = FPN*v_recyc*(Fe_s/Fe_s0) + FPN*v_abs - k_trueloss*Fe_p - k_recycreturn*Fe_p - k_ntbi*Fe_p
  dFe_s/dt = k_recycreturn*Fe_p - FPN*v_recyc*(Fe_s/Fe_s0) + k_ntbi*Fe_p

Stated baseline constants: k_deg_H = 0.0685 /h; H_normal = 23.0 nM; GAIN_IL6 = 4.15736; HFE
hepcidin-gain attenuation 5.4-fold; RBC-recycling iron flux 22.6 mg/day; dietary iron absorption
1.5 mg/day; total body iron 3.5 g (band 3-4); FPN0 = 0.5; FPN fast-relaxation rate 3.1 /h; NTBI-like
hepatocyte uptake 0.5x true obligate loss; blood volume / hematocrit 5.0 L / 0.45; TSAT0 = 30%.
Everything not literally handed a number (Fe_p0 in mg, TIBC concentration, plasma volume,
k_recycreturn, k_ntbi, k_synFPN, k_intFPN, k_synH) is solved algebraically from those baseline facts
(a fixed point at Fe_p0/H0/FPN0 with TSAT0 = 30%) -- disclosed algebra, not tuned.

Reads: nothing. Writes: blind_rebuild_results.json under the cell output directory.
Gate: the closed-form ceiling TSAT_ceiling = TSAT0/FPN0 must not be exceeded anywhere in the
numerically integrated HFE trajectory (key "ceiling_ever_exceeded"); the recorded 10 y / 25 y values
above that ceiling are thereby refuted.
"""
import json
import os as _os

import numpy as np
from scipy.integrate import solve_ivp

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "iron_hepcidin_conservation_blind_check")
OUT_PATH = _os.path.join(OUT_DIR, "blind_rebuild_results.json")

OUT = {}

# ---------------------------------------------------------------------------
# STEP A -- SYMBOLIC / CLOSED-FORM DERIVATION (done before any integration)
# ---------------------------------------------------------------------------
# Sum the two stated mass equations:
#   dFe_p/dt + dFe_s/dt
#     = [FPN*v_recyc*(Fe_s/Fe_s0) + FPN*v_abs - k_trueloss*Fe_p - k_recycreturn*Fe_p - k_ntbi*Fe_p]
#     + [k_recycreturn*Fe_p - FPN*v_recyc*(Fe_s/Fe_s0) + k_ntbi*Fe_p]
#   The FPN*v_recyc*(Fe_s/Fe_s0) terms cancel (+ / -).
#   The k_recycreturn*Fe_p terms cancel (- / +).
#   The k_ntbi*Fe_p terms cancel (- / +).
#   => d(Fe_p+Fe_s)/dt = FPN*v_abs - k_trueloss*Fe_p   (EXACT, algebraic, no approximation)
#
# Total body iron leaves the modeled system by exactly ONE route: k_trueloss*Fe_p.
# It enters by exactly one route: FPN*v_abs (dietary absorption gated by ferroportin).
# At any fixed point of the FULL 4-state system, d(Fe_p+Fe_s)/dt=0 individually forces:
#     Fe_p_ss = FPN_ss * v_abs / k_trueloss                       (*)
# FPN is STRUCTURALLY bounded above by 1: as FPN->1, dFPN/dt -> -k_intFPN*H_eff*FPN <= 0
# (the (1-FPN) synthesis-toward-max term vanishes, only the H-driven internalization/decay
# term remains, which is <=0), so FPN can never exceed 1 for any H_eff>=0.
# => Fe_p_ss <= v_abs / k_trueloss                                (hard ceiling, any hepcidin gain)
#
# At the WT (normal) fixed point itself, (*) gives:
#     Fe_p0 = FPN0 * v_abs / k_trueloss  =>  v_abs/k_trueloss = Fe_p0 / FPN0 = 2 * Fe_p0  (FPN0=0.5)
# So the ABSOLUTE ceiling (FPN_ss -> 1, i.e. total loss of hepcidin/FPN suppression, the HFE limit)
# is Fe_p_ss_max = v_abs/k_trueloss = 2*Fe_p0 -- exactly TWICE the baseline plasma iron.
# TSAT is proportional to Fe_p (TIBC/transferrin treated as a constant, non-dynamic quantity in
# this 4-state model -- it is not one of the 4 states), so:
#     TSAT_ceiling = 2 * TSAT0 = TSAT0 / FPN0
FPN0 = 0.5
TSAT0 = 0.30
tsat_ceiling_closed_form = TSAT0 / FPN0
OUT["closed_form_derivation"] = {
    "conserved_sum_ODE": "d(Fe_p+Fe_s)/dt = FPN*v_abs - k_trueloss*Fe_p  (all recycling/NTBI/recycreturn terms cancel exactly)",
    "fixed_point_identity": "Fe_p_ss = FPN_ss * v_abs / k_trueloss",
    "structural_FPN_bound": "FPN(t) < 1 for all t (from dFPN/dt = k_synFPN*(1-FPN)*(...) - k_intFPN*H_eff*FPN; as FPN->1 the synthesis term ->0 and the remaining term is <=0)",
    "baseline_calibration": "Fe_p0 = FPN0 * v_abs/k_trueloss  =>  v_abs/k_trueloss = Fe_p0/FPN0 = 2*Fe_p0 (using stated FPN0=0.5)",
    "TSAT_ceiling_formula": "TSAT_ceiling = TSAT0 / FPN0",
    "TSAT_ceiling_value_pct": tsat_ceiling_closed_form * 100.0,
    "note": "This closed-form ceiling needs ONLY TSAT0=30% and FPN0=0.5 (both stated baseline constants) -- it requires NO absolute mg/TIBC/plasma-volume unit convention at all, and is completely independent of v_abs, k_trueloss, v_recyc, k_recycreturn, k_ntbi individually (only their baseline RATIO matters, and that ratio is pinned by the baseline fixed point itself).",
}
print("CLOSED-FORM TSAT ceiling (%) =", tsat_ceiling_closed_form * 100.0)

# ---------------------------------------------------------------------------
# STEP B -- full numerical integration, to check the closed-form ceiling is respected
# and to get a comparable long-horizon trajectory number (unit conventions disclosed below;
# the closed-form result above does not depend on them).
# ---------------------------------------------------------------------------

# Disclosed unit convention:
# standard clinical TIBC = 300 ug/dL = 3.0 mg/L (textbook reference-range midpoint).
TIBC_conc_mg_per_L = 3.0
# plasma volume: FORCED (not a free choice) by the stated "blood volume/hematocrit = 5.0 L / 0.45"
# via plasma_vol = blood_volume * (1 - hematocrit).
blood_volume_L = 5.0
hematocrit = 0.45
plasma_vol_L = blood_volume_L * (1 - hematocrit)  # = 2.75 L

Fe_p0 = TSAT0 * TIBC_conc_mg_per_L * plasma_vol_L  # mg, baseline plasma iron mass
print("Fe_p0 (mg) =", Fe_p0, " plasma_vol_L =", plasma_vol_L)

# Baseline daily fluxes, stated directly:
absorption_loss_mg_per_day = 1.5   # stated "dietary iron absorption = 1.5 mg/day" (== obligate loss at steady state)
recycling_mg_per_day = 22.6        # stated "RBC-recycling iron flux = 22.6 mg/day"

# k_trueloss * Fe_p0 = absorption_loss (baseline balance of the SUM equation) => k_trueloss:
k_trueloss_per_day = absorption_loss_mg_per_day / Fe_p0
k_trueloss = k_trueloss_per_day / 24.0  # per hour

# FPN0*v_abs = k_trueloss*Fe_p0 = absorption_loss => v_abs:
v_abs_per_day = absorption_loss_mg_per_day / FPN0
v_abs = v_abs_per_day / 24.0  # mg/h

# NTBI fraction: stated "0.5x true obligate loss"
k_ntbi_per_day = 0.5 * k_trueloss_per_day
k_ntbi = k_ntbi_per_day / 24.0

# FPN0*v_recyc = (k_recycreturn+k_ntbi)*Fe_p0  =>  k_recycreturn:
v_recyc_per_day = recycling_mg_per_day
v_recyc = v_recyc_per_day / 24.0
k_recycreturn_per_day = (FPN0 * v_recyc_per_day) / Fe_p0 - k_ntbi_per_day
k_recycreturn = k_recycreturn_per_day / 24.0

# Fe_s0: not pinned by any Fe_p-ceiling-relevant equation (only the Fe_p/FPN/H subsystem
# matters for the TSAT ceiling); use standard reference storage-iron mass (~1 g) as a
# disclosed, non-load-bearing convention for running the full 4-state integration.
Fe_s0 = 1000.0  # mg

# Hepcidin / FPN kinetics, solved from the stated baseline fixed point + stated relaxation rate:
k_degH = 0.0685      # /h, stated
H0 = 23.0            # nM, stated
k_synH = k_degH * H0  # baseline balance of dH/dt at Fe_p=Fe_p0, IL6=0

# FPN baseline balance: k_synFPN*(1-FPN0) = k_intFPN*H0*FPN0  => k_synFPN = H0*k_intFPN (FPN0=0.5 cancels)
# stated "FPN fast-relaxation rate = 3.1 /h" = linearized decay rate at the fixed point
#   = k_synFPN + k_intFPN*H0  (standard linearization of a(1-x)-b*x type equation about x*)
# => k_synFPN + k_intFPN*H0 = 3.1, and k_synFPN = H0*k_intFPN => 2*H0*k_intFPN = 3.1
FPN_relax_rate = 3.1  # /h, stated
k_intFPN = FPN_relax_rate / (2.0 * H0)
k_synFPN = H0 * k_intFPN

GAIN_IL6 = 4.15736  # stated (irrelevant for chronic no-inflammation HFE run below)
HFE_ATTEN = 5.4     # stated

OUT["solved_rate_constants"] = {
    "k_trueloss_per_h": k_trueloss, "v_abs_mg_per_h": v_abs, "k_ntbi_per_h": k_ntbi,
    "v_recyc_mg_per_h": v_recyc, "k_recycreturn_per_h": k_recycreturn,
    "k_synH_nM_per_h": k_synH, "k_degH_per_h": k_degH,
    "k_synFPN_per_h": k_synFPN, "k_intFPN_per_nM_h": k_intFPN,
    "Fe_p0_mg": Fe_p0, "Fe_s0_mg": Fe_s0, "H0_nM": H0, "FPN0": FPN0,
}


def rhs(t, y, k_synH_local):
    Fe_p, H, FPN, Fe_s = y
    IL6 = 0.0  # chronic HFE scenario: no acute inflammation pulse
    H_eff = H  # no drug/antagonist in this run
    dH = k_synH_local * (Fe_p / Fe_p0) * np.exp(GAIN_IL6 * IL6) - k_degH * H
    dFPN = k_synFPN * (1 - FPN) * (1 - 0.6 * IL6) - k_intFPN * H_eff * FPN
    dFe_p = FPN * v_recyc * (Fe_s / Fe_s0) + FPN * v_abs - k_trueloss * Fe_p - k_recycreturn * Fe_p - k_ntbi * Fe_p
    dFe_s = k_recycreturn * Fe_p - FPN * v_recyc * (Fe_s / Fe_s0) + k_ntbi * Fe_p
    return [dFe_p, dH, dFPN, dFe_s]


y0 = [Fe_p0, H0, FPN0, Fe_s0]

# sanity: baseline (WT, k_synH unattenuated) must be an exact fixed point
r0 = rhs(0, y0, k_synH)
OUT["baseline_fixed_point_residual"] = r0
print("baseline residual (should be ~0):", r0)

# WT long run (should stay at baseline)
sol_wt = solve_ivp(rhs, [0, 25 * 365 * 24], y0, args=(k_synH,), method="Radau",
                    t_eval=[y * 365 * 24 for y in [1, 5, 10, 25]], rtol=1e-9, atol=1e-12)

# HFE run: hepcidin synthesis gain attenuated 5.4x (stated mechanism)
k_synH_hfe = k_synH / HFE_ATTEN
sol_hfe = solve_ivp(rhs, [0, 25 * 365 * 24], y0, args=(k_synH_hfe,), method="Radau",
                     t_eval=[y * 365 * 24 for y in [1, 5, 10, 25]], rtol=1e-9, atol=1e-12)

years = [1, 5, 10, 25]


def tsat_series(sol):
    Fe_p_vals = sol.y[0]
    return [100.0 * fp / Fe_p0 * TSAT0 for fp in Fe_p_vals]


tsat_wt = tsat_series(sol_wt)
tsat_hfe = tsat_series(sol_hfe)

OUT["integration"] = {
    "years": years,
    "TSAT_pct_WT": tsat_wt,
    "TSAT_pct_HFE": tsat_hfe,
    "TSAT_ceiling_closed_form_pct": tsat_ceiling_closed_form * 100.0,
    "recorded_cell_HFE_TSAT_pct": {"1y": 58.5, "5y": 87.0, "10y": 107.4, "25y": 136.3},
    "max_TSAT_HFE_over_horizon_pct": max(tsat_hfe),
    "ceiling_ever_exceeded": any(v > tsat_ceiling_closed_form * 100.0 + 1e-6 for v in tsat_hfe),
}

print("\nHFE TSAT%% trajectory (blind rebuild):")
for yv, t in zip(years, tsat_hfe):
    print(f"  {yv}y: {t:.2f}%")
print("Closed-form ceiling: %.2f%%" % (tsat_ceiling_closed_form * 100))
print("Recorded claim: 58.5% (1y) / 87.0% (5y) / 107.4% (10y) / 136.3% (25y)")
print("Ceiling ever exceeded in this integration?", OUT["integration"]["ceiling_ever_exceeded"])

_os.makedirs(OUT_DIR, exist_ok=True)
with open(OUT_PATH, "w") as f:
    json.dump(OUT, f, indent=2, default=float)

print("\nWROTE " + OUT_PATH)
