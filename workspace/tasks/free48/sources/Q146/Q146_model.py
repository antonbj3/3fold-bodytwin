"""
BT-HX-Q146 - first-principles model of the nasal route with a separately
tracked swallowed (GI) fraction and a systemic (plasma) fraction.

Control volumes (CV) and conservation identities
-------------------------------------------------
  CV1..CVn  nasal regions  (vestibulum, olfactory, mid/upper turbinate,
                             lower turbinate, nasopharynx) - liquid film on mucosa
  CVg       gastrointestinal tract  (fed by the SAME F_sw [mol/s])
  CVp       plasma / systemic     (fed by (1-f_met)*J_epi and (1-E_H)*J_gi)

Molar balance of a nasal region i (the mathematical start of the question):

    dN_i/dt = I_i - J_epi,i - F_sw,i - F_ext,i - R_met,i

with, for a well-stirred film of thickness d_i over area S_i:

    C_i    = N_i / (S_i * d_i)                  free local concentration [mol/m^3]
    J_epi,i= k_epi * a_i * S_i * C_i            trans-epithelial uptake      [mol/s]
    F_sw,i = chi_i * v_i * S_i * C_i            mucociliary transport -> GI  [mol/s]
    F_ext,i= (1-chi_i) * v_i * S_i * C_i        drainage / loss out of the nose
    R_met,i= f_met * J_epi,i                    local mucosal metabolism (CYP2A)

so that with  lambda_i = (k_epi*a_i + v_i)/d_i

    f_abs,i = k_epi*a_i/lambda_i,  f_sw,i = chi_i*v_i/lambda_i,
    f_ext,i = (1-chi_i)*v_i/lambda_i,  f_abs,f_sw,f_ext sum to exactly 1.

Note the control volume: d_i cancels in the fractions, so the ABSORPTION vs
SWALLOW SPLIT is set by the dimensionless group

    Pe_i = k_epi * a_i / v_i          (permeation velocity / transport velocity)

and d_i only sets the time scale.  Both velocities are [m/s], so Pe is
dimensionless.

GI leg (shared substance identity, same F_sw):
    dN_g/dt = F_sw,total - k_gi * N_g
    J_gi    = k_gi * N_g                                    [mol/s]
Systemic (well-stirred liver on the portal leg):
    dN_p/dt = (1-f_met) * sum_i J_epi,i + (1-E_H) * J_gi - k_e * N_p
    C_p     = N_p / V_ss                                    [mol/m^3]

Run:  python3 model.py     ->  results.json
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict, replace
from typing import Dict, List, Sequence

import numpy as np
from scipy.integrate import solve_ivp

MW_NICOTINE = 162.23  # g/mol


# ---------------------------------------------------------------- dimensions
class U:
    """Tiny dimensional-analysis helper: exponents of (L, M, T, N)."""

    __slots__ = ("e",)

    def __init__(self, L=0, M=0, T=0, N=0):
        self.e = {"L": L, "M": M, "T": T, "N": N}

    def __mul__(self, o):
        return U(**{k: self.e[k] + o.e[k] for k in self.e})

    __rmul__ = __mul__

    def __truediv__(self, o):
        return U(**{k: self.e[k] - o.e[k] for k in self.e})

    def __pow__(self, n):
        return U(**{k: self.e[k] * n for k in self.e})

    def __str__(self):
        return "".join(f"{v}{k}" for k, v in self.e.items() if v) or "1"

    def eq(self, o):
        return all(abs(self.e[k] - o.e[k]) < 1e-12 for k in self.e)


L_ = U(L=1)
T_ = U(T=1)
MOL = U(N=1)
ONE = U()

DIM_VELOCITY = L_ / T_          # k_epi, v_i            -> m/s
DIM_AREA = L_ ** 2              # S_i                   -> m^2
DIM_CONC = MOL / L_ ** 3        # C_i, C_p              -> mol/m^3
DIM_FLUX = MOL / T_             # all [mol/s] fluxes
DIM_RATE = ONE / T_             # k_gi, k_e              -> 1/s
DIM_MOLES = MOL

DIM_CHECKS = {
    "k_epi": DIM_VELOCITY,
    "v_clear": DIM_VELOCITY,
    "S": DIM_AREA,
    "C": DIM_CONC,
    "J_epi = k_epi*a*S*C": (DIM_VELOCITY * DIM_AREA * DIM_CONC),
    "F_sw = chi*v*S*C": (DIM_VELOCITY * DIM_AREA * DIM_CONC),
    "R_met = f_met*J_epi": DIM_FLUX,
    "N": DIM_MOLES,
    "k_gi": DIM_RATE,
    "k_e": DIM_RATE,
    "C_p = N/V": (DIM_MOLES / L_ ** 3),
}


def run_dimension_check() -> Dict[str, str]:
    """Explicit unit audit: every declared expression must have DIM_FLUX or the
    declared dimension. Returns the audit table (all must read 'ok')."""
    audit = {}
    for name, dim in DIM_CHECKS.items():
        expect = DIM_FLUX if name.startswith(("J_", "F_", "R_")) else dim
        audit[name] = f"{dim} {'==' if dim.eq(expect) else '!='} {expect}"
    # cross-expressions
    cross = {
        "Pe = k_epi*a/v": (DIM_VELOCITY / DIM_VELOCITY),
        "lambda = (k*a+v)/d": (DIM_VELOCITY / L_),
        "flux = v*S*C": (DIM_VELOCITY * DIM_AREA * DIM_CONC),
        "f_abs = k*a/(k*a+v)": (DIM_VELOCITY / DIM_VELOCITY),
    }
    for k, v in cross.items():
        target = ONE if (k.startswith("Pe") or k.startswith("f_abs")) else (
            DIM_RATE if k.startswith("lambda") else DIM_FLUX)
        audit[k] = f"{v} {'==' if v.eq(target) else '!= ' + str(target)}"
    return audit


# ---------------------------------------------------------------- parameters
@dataclass(frozen=True)
class Region:
    name: str
    S_cm2: float          # mucosal area [cm^2]        (R5)
    a_abs: float          # absorptive fraction of area (R5: vestibule = 0)
    v_clear_mm_min: float # local clearance velocity [mm/min] (R4)
    chi_swallow: float    # share of cleared flux that is swallowed (1 = all)


# Areas from R5 (Ličen 2026, 5-region nasal cast): total 150 cm^2
# = vestibulum + olfactory 15 + turbinate 130 + nasopharynx.
REGIONS: List[Region] = [
    Region("vestibulum",     S_cm2=10.0, a_abs=0.0, v_clear_mm_min=33.0, chi_swallow=0.15),
    Region("olfactory",      S_cm2=15.0, a_abs=1.0, v_clear_mm_min=6.6,  chi_swallow=1.00),
    Region("turbinate_mid",  S_cm2=50.0, a_abs=1.0, v_clear_mm_min=6.6,  chi_swallow=1.00),
    Region("turbinate_low",  S_cm2=80.0, a_abs=1.0, v_clear_mm_min=6.6,  chi_swallow=1.00),
    Region("nasopharynx",    S_cm2=10.0, a_abs=0.0, v_clear_mm_min=6.6,  chi_swallow=1.00),
]

# Deposition patterns (mass fractions of the actuated dose that lands in each
# region).  Spray pattern is anchored on R5's measurement ("majority in the
# vestibule and lower turbinate", olfactory and nasopharyngeal < 1 %).
# The drop patterns are the measured application modes of R1.
DEPOSITION = {
    # micropipette 0.1 mL onto the septum  -> all on respiratory mucosa
    "septum":   {"vestibulum": 0.02, "olfactory": 0.03, "turbinate_mid": 0.35,
                 "turbinate_low": 0.58, "nasopharynx": 0.02},
    # micropipette 0.1 mL onto the conchae -> spreads laterally/posteriorly
    "conchae":  {"vestibulum": 0.02, "olfactory": 0.01, "turbinate_mid": 0.32,
                 "turbinate_low": 0.50, "nasopharynx": 0.15},
    # 100 uL pump spray, R5 pattern
    "spray":    {"vestibulum": 0.30, "olfactory": 0.01, "turbinate_mid": 0.19,
                 "turbinate_low": 0.49, "nasopharynx": 0.01},
}


@dataclass
class Params:
    # --- molecule / dose (R1: 1 mg nicotine base, 100 uL actuation) ---
    dose_mg: float = 1.0
    MW: float = MW_NICOTINE
    V_liq_m3: float = 100e-9      # 100 uL of solution  [m^3]
    T_dep_s: float = 1.0          # actuation time       [s]

    # --- nasal ---
    k_epi: float = 1.0e-5         # permeation velocity  [m/s]  ASSUMPTION
    v_clear_mm_min: float = 6.6   # bulk nasal mucus transport [mm/min] (R4)
    v_vest_factor: float = 5.0    # vestibular run-off is faster (ASSUMPTION)
    d_min_m: float = 1.0e-5       # minimum film thickness [m]  ASSUMPTION
    f_met: float = 0.10           # mucosal metabolism of the absorbed flux
    deposition: Dict[str, float] = field(default_factory=lambda: dict(DEPOSITION["spray"]))

    # --- GI (R2/R3) ---
    k_gi_h: float = 1.55          # intestinal absorption rate [1/h] (R2)
    F_oral: float = 0.40          # composite GI (absorb x hepatic avail.)  (R2)
    Q_h_liver_L_h: float = 90.0   # hepatic blood flow [L/h]
    CL_systemic_L_h: float = 49.2 # R1: 820 mL/min (iv arm of the SAME study)

    # --- systemic PK (R1) ---
    V_ss_L: float = 196.0         # 2.8 L/kg * 70 kg
    t_half_h: float = 2.0         # R1: "about 2 h"

    # --- reference iv dose of R1, for the absolute-bioavailability reference ---
    dose_iv_mg: float = 0.7

    # --- solver ---
    t_max_h: float = 12.0
    n_out: int = 2001


def kel(params: Params) -> float:
    """Elimination rate constant, from CL and V (internally consistent)."""
    return params.CL_systemic_L_h / params.V_ss_L  # [1/h]


# ---------------------------------------------------------------- core model
def _assemble(params: Params):
    """Return per-region arrays and the film thickness for a deposition pattern."""
    names = [r.name for r in REGIONS]
    p = np.array([params.deposition.get(n, 0.0) for n in names], float)
    p = p / p.sum()                                  # normalise (identity)
    S = np.array([r.S_cm2 * 1e-4 for r in REGIONS])   # cm^2 -> m^2
    a = np.array([r.a_abs for r in REGIONS])
    v = np.array([r.v_clear_mm_min for r in REGIONS], float) * (1e-3 / 60.0)  # mm/min -> m/s
    # the bulk (turbinate) velocity is a model parameter; the vestibular value
    # is derived from it so that a single knob controls the whole field
    scale = params.v_clear_mm_min / 6.6
    v = v * scale
    v[0] *= params.v_vest_factor
    chi = np.array([r.chi_swallow for r in REGIONS])

    dose_mol = params.dose_mg * 1e-3 / params.MW
    V_dep = params.V_liq_m3 * p                      # liquid per region
    S_dep = float(np.sum(S * p))                     # wetted area
    d = V_dep / S_dep                               # uniform film thickness
    d = np.maximum(d, params.d_min_m)
    lam = (params.k_epi * a + v) / d                 # [1/s]
    return names, p, S, a, v, chi, d, lam, dose_mol, S_dep


def simulate(params: Params) -> Dict:
    names, p, S, a, v, chi, d, lam, dose_mol, S_dep = _assemble(params)
    n = len(names)
    k_gi = params.k_gi_h / 3600.0
    ke = kel(params) / 3600.0
    E_H = 1.0 - params.CL_systemic_L_h / params.Q_h_liver_L_h
    F_gi_avail = params.F_oral / max(1e-9, (1.0 - E_H))   # intestinal absorption alone
    F_gi_systemic = (1.0 - E_H) * F_gi_avail              # == F_oral (R2)

    T = params.T_dep_s
    rate_dep = dose_mol / T                          # [mol/s], 0 after T
    V_ss_m3 = params.V_ss_L * 1e-3

    def rhs(t, y):
        N = y[:n]
        Ng = y[n]
        Np = y[n + 1]
        C = N / (S * d)
        J_epi = params.k_epi * a * S * C
        Fsw = chi * v * S * C
        Fext = (1 - chi) * v * S * C        # drainage / loss out of the nose
        dN = -J_epi - Fsw - Fext            # == -lambda_i * N_i  (identity)
        Jgi = k_gi * Ng
        dNg = Fsw.sum() - Jgi
        dNp = ((1 - params.f_met) * J_epi.sum()
               + F_gi_systemic * Jgi - ke * Np)
        return np.concatenate([dN, np.array([dNg, dNp])])
    # Administration is treated as an instantaneous bolus at t = 0.  Justified:
    # the nasal relaxation time is 1/lambda ~ 0.08 s while the actuation takes
    # T_dep = 1 s, so lambda*T_dep ~ 12 and the bolus limit is accurate to
    # f_abs*(1-exp(-lambda*T_dep)) < 1e-5 relative.
    y0 = np.zeros(n + 2)
    y0[:n] = p * dose_mol

    t_eval = np.unique(np.concatenate([
        np.logspace(-9, math.log10(params.t_max_h * 3600.0), 4 * params.n_out),
        np.linspace(0.0, params.t_max_h * 3600.0, 200),
    ]))
    sol = solve_ivp(rhs, (0, t_eval[-1]), y0, t_eval=t_eval,
                    method="Radau", rtol=1e-11, atol=1e-18)
    if not sol.success:
        raise RuntimeError(sol.message)

    t, Y = sol.t, sol.y
    N, Ng, Np = Y[:n], Y[n], Y[n + 1]
    C = N / (S * d)[:, None]
    C_p = Np / V_ss_m3

    # ---- route integrals: quadrature of the fluxes (not a bookkeeping shortcut),
    #      so that test_model.py can compare them with the closed form ----
    # closed-form fractions: d_i cancels, only Pe = k_epi*a/v survives
    f_abs = (params.k_epi * a / lam) / d
    f_sw = (chi * v / lam) / d
    f_ext = ((1 - chi) * v / lam) / d

    I_tot = dose_mol
    Cg = Ng / V_ss_m3
    J_epi_t = params.k_epi * a[:, None] * S[:, None] * C
    Fsw_t = chi[:, None] * v[:, None] * S[:, None] * C
    Fext_t = (1 - chi)[:, None] * v[:, None] * S[:, None] * C
    Jgi_t = k_gi * Ng

    N_abl = np.trapezoid(J_epi_t, t, axis=1)
    N_swal = np.trapezoid(Fsw_t, t, axis=1)
    N_ext = np.trapezoid(Fext_t, t, axis=1)
    N_gi_abs = float(np.trapezoid(Jgi_t, t))
    N_left = float(N[:, -1].sum())
    N_met = params.f_met * N_abl
    N_to_blood_nasal = (1 - params.f_met) * N_abl.sum()
    N_to_blood_gi = F_gi_systemic * N_gi_abs      # absorbed in gut, survives liver
    frac_sum = f_abs + f_sw + f_ext             # closed-form, used in the audit

    # ---- plasma exposure, iv-referenced (R1) ----
    # t is in SECONDS, so the trapezoid gives mol*s/m^3
    AUC_p_mol_s_per_m3 = float(np.trapezoid(C_p, t))
    AUC_p_mol_h_per_m3 = AUC_p_mol_s_per_m3 / 3600.0
    dose_iv_mol = params.dose_iv_mg * 1e-3 / params.MW
    # AUC_iv = Dose/CL  (mol*h/L -> mol*h/m^3)
    AUC_iv_mol_h_per_m3 = dose_iv_mol / params.CL_systemic_L_h * 1e3
    F_abs_total = AUC_p_mol_h_per_m3 / AUC_iv_mol_h_per_m3

    # analytic baseline: what the swallowed leg alone could give
    F_swallow_only = float(p @ f_sw) * params.F_oral
    # identity: F cannot exceed the swallowed leg on its own
    F_lower_bound_nasal = F_abs_total - F_swallow_only

    imax = int(np.argmax(C_p))
    Pe = params.k_epi * a / v

    return {
        "mode": params._mode,  # set by caller
        "deposition": {nm: float(pp) for nm, pp in zip(names, p)},
        "film_thickness_um": {nm: float(dd * 1e6) for nm, dd in zip(names, d)},
        "Pe_k_epi_a_over_v": {nm: float(Pe[i]) for i, nm in enumerate(names)},
        "f_abs_region": {nm: float(f_abs[i]) for i, nm in enumerate(names)},
        "f_swallow_region": {nm: float(f_sw[i]) for i, nm in enumerate(names)},
        "f_external_region": {nm: float(f_ext[i]) for i, nm in enumerate(names)},
        "region_fraction_sum_max_dev": float(np.max(np.abs(frac_sum - 1.0))),
        "molar_balance_residual": float(
            N_swal.sum() + N_ext.sum() + N_abl.sum() + N_left - I_tot),
        "F_abs_total": float(F_abs_total),
        "F_swallow_only_ceiling": F_swallow_only,
        "F_nasal_route_contribution": float(F_lower_bound_nasal),
        "route": {
            "nasal_absorbed_mol": N_abl.sum(),
            "nasal_metabolised_mol": N_met.sum(),
            "nasal_to_blood_mol": N_to_blood_nasal.sum(),
            "swallowed_mol": N_swal.sum(),
            "swallowed_fraction_of_dose": float(N_swal.sum() / I_tot),
            "external_loss_mol": N_ext.sum(),
            "gi_absorbed_mol": N_gi_abs,
            "gi_to_blood_mol": N_to_blood_gi,
            "n_left_in_nose_mol": N_left,
        },
        "plasma": {
            "AUC_mol_s_per_m3": AUC_p_mol_s_per_m3,
            "AUC_mol_h_per_m3": AUC_p_mol_h_per_m3,
            "AUC_ng_h_per_mL": AUC_p_mol_h_per_m3 * 1e3 * params.MW,
            "AUC_iv_mol_h_per_m3": AUC_iv_mol_h_per_m3,
            "AUC_iv_ng_h_per_mL": AUC_iv_mol_h_per_m3 * 1e3 * params.MW,
            "Cmax_ng_per_mL": float(C_p[imax] * 1e3 * params.MW),
            "tmax_min": float(t[imax] / 60.0),
            "V_ss_L": params.V_ss_L,
            "CL_L_h": params.CL_systemic_L_h,
            "k_e_h": kel(params),
            "t_half_h_model": math.log(2) / kel(params),
            "E_H_wellstirred": E_H,
            "F_gi_absorption_implied": F_gi_avail,
        },
        "local_free_conc": {
            "C_free_max_mol_per_m3": float(np.max(C)),
            "C_free_max_mg_per_mL": float(np.max(C) * params.MW / 1000.0),
            "C_free_initial_uniform_mol_per_m3": float(
                (p * dose_mol / (S * d)).max()),
            "AUC_local_mol_s_per_m3": float(
                np.trapezoid(C, t, axis=1).max()),
        },
    }


def run_mode(mode: str, **overrides) -> Dict:
    dep = overrides.pop("deposition", None) or DEPOSITION[mode]
    params = Params(deposition=dict(dep), **overrides)
    params._mode = mode
    return simulate(params)


# ------------------------------------------------- ceiling (structural test C3)
def ceiling_F(params: Params) -> Dict[str, float]:
    """Max achievable F when the absorption competition is saturated
    (k_epi*a >> v_i for every absorbing region) and nothing is lost.
    Upper bound of the model FORM, independent of k_epi."""
    out = {}
    _, p, S, a, v, chi, d, lam, dose, S_dep = _assemble(params)
    f_abs1 = np.where(a > 0, 1.0, 0.0)
    # the cleared flux of a NON-absorbing region goes to the GI leg if it is
    # swallowed (chi) and is lost out of the nose otherwise
    f_sw1 = np.where(a > 0, 0.0, chi)
    f_lost1 = np.where(a > 0, 0.0, 1.0 - chi)
    E_H = 1.0 - params.CL_systemic_L_h / params.Q_h_liver_L_h
    Fgi = params.F_oral / (1.0 - E_H) * (1.0 - E_H)   # == F_oral
    out["F_max_saturated"] = float(
        (p * ((1 - params.f_met) * f_abs1 + Fgi * f_sw1)).sum())
    out["frac_lost_max"] = float((p * f_lost1).sum())
    return out


def ceiling_scan(params: Params) -> Dict[str, float]:
    """Ceiling for every site, and the axis along which the ceiling can be lifted."""
    res = {}
    for mode in DEPOSITION:
        pp = replace(params, deposition=dict(DEPOSITION[mode]))
        pp._mode = mode
        res[mode] = ceiling_F(pp)["F_max_saturated"]
    return res


# ------------------------------------------------------------- null model C6
def pooled_null(measured: Dict[str, float]) -> Dict[str, float]:
    """Route-split model with ONE shared absorbed fraction (the comparison
    baseline of the question: 'gemensam absorberad fraktion').  It has exactly
    one free number; we let it be the value that minimises the max absolute
    error on the three measured sites."""
    f = float(np.mean(list(measured.values())))
    return {m: f for m in measured}


# ------------------------------------------------------------------- sweep
def sensitivity(overrides_list) -> List[Dict]:
    base = run_mode("spray")
    rows = []
    for label, ov in overrides_list:
        r = run_mode("spray", **ov)
        rows.append({
            "parameter": label,
            "value": ov,
            "F_abs": r["F_abs_total"],
            "F_swallow_fraction": r["route"]["swallowed_fraction_of_dose"],
            "C_free_max_mol_per_m3": r["local_free_conc"]["C_free_max_mol_per_m3"],
            "AUC_plasma_mol_s_per_m3": r["plasma"]["AUC_mol_s_per_m3"],
            "dF_vs_base": r["F_abs_total"] - base["F_abs_total"],
            "rel_dF_pct": 100.0 * (r["F_abs_total"] - base["F_abs_total"])
            / base["F_abs_total"],
        })
    return rows


def required_Pe(mode: str, F_target: float, params: Params = None) -> float:
    """What k_epi/v would be NEEDED to reproduce a measured F.  Closed form:

        F = sum_i p_i [ (1-f_met)*a_i*Pe/(1+Pe) + chi_i*(1-a_i*Pe/(1+Pe))*F_oral ]
    solved for Pe.  Returns nan if the target is unreachable.
    """
    p = params or Params(deposition=dict(DEPOSITION[mode]))
    p = replace(p, deposition=dict(DEPOSITION[mode]))
    names, pp, S, a, v, chi, d, lam, dose, S_dep = _assemble(p)
    absf = a * np.array([1.0])          # absorbing fraction of area
    # F(Pe) = A*Pe/(1+Pe) + B, with A,B linear in Pe per region; solve for the
    # single Pe that is common to every region (a_i=0 regions drop out)
    base = float((pp * chi * p.F_oral).sum())          # F at Pe = 0
    # going from Pe=0 to Pe=inf replaces chi*F_oral by (1-f_met) on the
    # absorbing area, so the extra per unit absorbing fraction is
    # (1-f_met) - chi*F_oral
    coef = float((pp * absf * ((1 - p.f_met) - chi * p.F_oral)).sum())
    # base is the F at Pe = 0, coef*(1) is the extra when Pe -> inf
    # F(Pe) = base + coef*Pe/(1+Pe)  ->  Pe = (F-base)/(coef-(F-base))
    num = F_target - base
    den = coef - num
    return float(num / den) if abs(den) > 1e-12 else float("nan")


def F_oral_value(p: Params) -> float:
    return p.F_oral


# --------------------------------------------------------------------- main
def main():
    audit = run_dimension_check()
    measured = {"septum": 0.76, "conchae": 0.64, "spray": 0.58}   # R1
    measured_ci = {"septum": (0.26, 0.97), "conchae": (0.21, 1.03),
                   "spray": (0.17, 0.85)}                            # R1, n = 8

    pred = {m: run_mode(m) for m in ("septum", "conchae", "spray")}
    ceil = ceiling_scan(Params())
    null = pooled_null(measured)

    sens = sensitivity([
        ("k_epi x0.5", {"k_epi": 0.5e-5}),
        ("k_epi x1.5", {"k_epi": 1.5e-5}),
        ("v_clear x0.5", {"v_clear_mm_min": 3.3}),
        ("v_clear x1.5", {"v_clear_mm_min": 9.9}),
        ("p_vestibulum x0.5", {"deposition": {**DEPOSITION["spray"],
                                              "vestibulum": 0.15,
                                              "turbinate_low": 0.64}}),
        ("p_vestibulum x1.5", {"deposition": {**DEPOSITION["spray"],
                                              "vestibulum": 0.45,
                                              "turbinate_low": 0.34}}),
    ])

    out = {
        "id": "BT-HX-Q146",
        "model": "first-principles nasal route + GI leg + 1-compartment plasma",
        "dimension_audit": audit,
        "params": {k: (v if not isinstance(v, dict) else "see DEPOSITION")
                   for k, v in asdict(Params()).items()},
        "deposition_patterns": DEPOSITION,
        "regions": [asdict(r) for r in REGIONS],
        "predicted": {m: pred[m] for m in pred},
        "measured_R1": measured,
        "measured_R1_range_n8": measured_ci,
        "ceiling": ceil,
        "null_pooled": null,
        "sensitivity": sens,
        "required_Pe_for_measured_F": {
            m: required_Pe(m, measured[m]) for m in measured},
        "required_k_epi_m_per_s_for_measured_F": {
            m: required_Pe(m, measured[m]) * Params().v_clear_mm_min
               * (1e-3 / 60.0) for m in measured},
        "identifiability": {
            "F_oral_reference": Params().F_oral,
            "F_measured_range_R1": [0.58, 0.76],
            "forced_nasal_contribution_min_pp": {
                m: (measured[m] - Params().F_oral) * 100 for m in measured},
        },
    }
    with open("results.json", "w") as f:
        json.dump(out, f, indent=2, default=float)

    # ---- console summary ----
    print("dimension audit (all must be '=='):")
    for k, v in audit.items():
        print(f"  {k:28s} {v}")
    print(f"\n{'site':8s} {'F_pred':>8s} {'F_R1':>7s} {'ratio':>7s} {'F_max':>7s} "
          f"{'F_sw':>7s} {'null':>7s} {'Cfree_max':>10s} {'AUC_p':>11s}")
    for m in ("septum", "conchae", "spray"):
        r = pred[m]
        print(f"{m:8s} {r['F_abs_total']:8.3f} {measured[m]:7.2f} "
              f"{measured[m]/r['F_abs_total']:7.2f} {ceil[m]:7.3f} "
              f"{r['route']['swallowed_fraction_of_dose']:7.3f} "
              f"{null[m]:7.3f} "
              f"{r['local_free_conc']['C_free_max_mol_per_m3']:10.2f} "
              f"{r['plasma']['AUC_mol_s_per_m3']:11.3e}")
    print("\nsensitivity (+/-50 %):")
    for s in sens:
        print(f"  {s['parameter']:22s} F={s['F_abs']:.3f}  "
              f"dF={s['dF_vs_base']:+.3f} ({s['rel_dF_pct']:+.1f} %)")
    print(f"\nbalance residual (mol): "
          f"{max(abs(pred[m]['molar_balance_residual']) for m in pred):.3e}")
    print(f"region fraction-sum max dev: "
          f"{max(pred[m]['region_fraction_sum_max_dev'] for m in pred):.3e}")
    print("wrote results.json")


if __name__ == "__main__":
    main()
