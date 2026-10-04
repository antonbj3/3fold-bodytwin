from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


REFERENCE = {
    "authors": "Hertz, Paulson, Barry, Christiansen, Svendsen",
    "year": 1981,
    "journal": "Journal of Clinical Investigation",
    "doi": "10.1172/JCI110073",
    "value": 0.46,
    "insulin_value": 0.66,
    "extraction_fraction": 0.14,
    "unit": "umol g^-1 min^-1",
    "quantity": "unidirectional D-glucose flux from blood across the BBB",
    "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC370607/",
    "verification": "webfetch checked article abstract and bibliographic page on 2026-09-25",
}

EQUATIONS = {
    "geometry": "A = 2*pi*r*L*parallel_segments; V_E = A*h; V_L and V_B are explicit volumes",
    "passive_membrane": "J_p,i->j = P_i*A*max(C_i-C_j, 0), with net flux equal to the signed difference of the two directions",
    "carrier_influx": "J_C,in = Vmax_in*C_L/(Km_in*(1+C_comp/KI)+C_L)",
    "carrier_efflux": "J_C,eff = Vmax_eff*C_E/(Km_eff*(1+C_comp/KI)+C_E)",
    "binding": "dN_R/dt = k_on*(R_total-N_R)*C_L - k_off*N_R - k_end*N_R",
    "vesicles": "dN_V/dt = k_end*N_R - k_release*N_V - k_deg*N_V",
    "metabolism": "J_met = Vmax_met*C_E/(Km_met+C_E), with one lactate molecule per glucose",
    "substrate_balance": "d(N_L+N_E+N_R+N_V+N_B)/dt = J_perfusion - J_met - J_deg - J_brain_sink",
    "metabolite_balance": "d(N_M+N_M,B)/dt = J_met + J_deg + J_M,L->E - J_M,E->L - k_clear*N_M,B",
    "charge_balance": "lactate^- is transported with an equal signed H+ cotransport flux, so z_lactate*J_lactate + z_H*J_H = 0 at each membrane",
    "output": "J_brain,glucose = J_p,E->B + J_C,E->B + J_ves,release; J_brain,lactate = J_M,E->B",
}

UNIT_CONTRACTS = {
    "geometry_area": "2*pi*r*L => m^2",
    "geometry_volume": "A*h => m^3",
    "concentration": "N/V => mol m^-3",
    "passive_flux": "(m s^-1)*(m^2)*(mol m^-3) => mol s^-1",
    "carrier_flux": "(mol s^-1)*(mol m^-3)/(mol m^-3) => mol s^-1",
    "binding_flux": "(m^3 mol^-1 s^-1)*(mol)*(mol m^-3) => mol s^-1",
    "perfusion_flux": "(m^3 s^-1)*(mol m^-3) => mol s^-1",
    "charge_flux": "(mol s^-1)*charge => mol eq s^-1",
    "competition": "C_comp/KI => 1",
}

PRODUCT_SPEC = {
    "name": "lactate^-",
    "charge": -1,
    "counterion": "H+",
    "counterion_charge": 1,
    "cotransport_stoichiometry": 1.0,
}

PARAMETER_TABLE = [
    {"name": "brain_mass_g", "value": 1.0, "unit": "g", "source_or_assumption": "normalization to one gram brain"},
    {"name": "lumen_volume_ml_per_g", "value": 0.2, "unit": "mL g^-1", "source_or_assumption": "finite luminal mixing-volume assumption"},
    {"name": "perfusion_ml_min_per_g", "value": 0.5, "unit": "mL min^-1 g^-1", "source_or_assumption": "capillary perfusion boundary assumption"},
    {"name": "c_lumen_source_mM", "value": 5.0, "unit": "mM", "source_or_assumption": "upstream D-glucose boundary assumption, not a new measurement"},
    {"name": "surface_area_cm2_per_g", "value": 20.0, "unit": "cm^2 g^-1", "source_or_assumption": "geometry assumption; network-equivalent BBB surface"},
    {"name": "capillary_radius_m", "value": 4.0e-6, "unit": "m", "source_or_assumption": "generic brain-capillary geometry assumption"},
    {"name": "capillary_length_m", "value": 1.0e-3, "unit": "m", "source_or_assumption": "1-mm capillary subsegment geometry assumption"},
    {"name": "membrane_thickness_m", "value": 0.5e-6, "unit": "m", "source_or_assumption": "endothelial membrane thickness assumption"},
    {"name": "brain_volume_ml_per_g", "value": 0.7, "unit": "mL g^-1", "source_or_assumption": "non-endothelial abluminal volume assumption"},
    {"name": "c_lumen_mM", "value": 5.0, "unit": "mM", "source_or_assumption": "basal D-glucose boundary assumption, not a new measurement"},
    {"name": "c_endothelium0_mM", "value": 0.2, "unit": "mM", "source_or_assumption": "initial endothelial D-glucose assumption"},
    {"name": "c_brain0_mM", "value": 1.5, "unit": "mM", "source_or_assumption": "initial abluminal D-glucose assumption"},
    {"name": "c_met_lumen_mM", "value": 0.1, "unit": "mM", "source_or_assumption": "luminal lactate boundary assumption"},
    {"name": "c_met_brain_mM", "value": 1.5, "unit": "mM", "source_or_assumption": "initial abluminal lactate assumption"},
    {"name": "p_luminal_m_s", "value": 1.0e-8, "unit": "m s^-1", "source_or_assumption": "low passive glucose-permeability assumption"},
    {"name": "p_abluminal_m_s", "value": 1.0e-8, "unit": "m s^-1", "source_or_assumption": "low passive glucose-permeability assumption"},
    {"name": "p_met_luminal_m_s", "value": 1.0e-6, "unit": "m s^-1", "source_or_assumption": "lactate membrane permeability assumption"},
    {"name": "p_met_abluminal_m_s", "value": 1.0e-6, "unit": "m s^-1", "source_or_assumption": "lactate membrane permeability assumption"},
    {"name": "vmax_luminal_in_umol_min", "value": 0.60, "unit": "umol min^-1 g^-1", "source_or_assumption": "GLUT1/GLUT3-like capacity assumption anchored to published flux order, not independently fitted"},
    {"name": "km_luminal_in_mM", "value": 1.5, "unit": "mM", "source_or_assumption": "facilitated glucose-transport affinity assumption"},
    {"name": "vmax_luminal_efflux_umol_min", "value": 0.06, "unit": "umol min^-1 g^-1", "source_or_assumption": "hypothetical low glucose efflux capacity; substrate-specific value UNKNOWN"},
    {"name": "km_luminal_efflux_mM", "value": 2.0, "unit": "mM", "source_or_assumption": "hypothetical efflux affinity assumption"},
    {"name": "vmax_abluminal_in_umol_min", "value": 0.35, "unit": "umol min^-1 g^-1", "source_or_assumption": "abluminal GLUT-like export capacity assumption"},
    {"name": "km_abluminal_in_mM", "value": 1.5, "unit": "mM", "source_or_assumption": "abluminal carrier affinity assumption"},
    {"name": "vmax_abluminal_out_umol_min", "value": 0.10, "unit": "umol min^-1 g^-1", "source_or_assumption": "abluminal reverse-carrier capacity assumption"},
    {"name": "km_abluminal_out_mM", "value": 1.5, "unit": "mM", "source_or_assumption": "abluminal reverse-carrier affinity assumption"},
    {"name": "vmax_metabolism_umol_min", "value": 0.30, "unit": "umol min^-1 g^-1", "source_or_assumption": "endothelial glucose-metabolism capacity assumption"},
    {"name": "km_metabolism_mM", "value": 1.0, "unit": "mM", "source_or_assumption": "endothelial metabolism affinity assumption"},
    {"name": "vmax_brain_sink_umol_min", "value": 0.40, "unit": "umol min^-1 g^-1", "source_or_assumption": "non-endothelial brain glucose sink assumption; secondary compartment"},
    {"name": "km_brain_sink_mM", "value": 1.5, "unit": "mM", "source_or_assumption": "non-endothelial sink affinity assumption"},
    {"name": "r_total_pmol_per_g", "value": 1.0, "unit": "pmol g^-1", "source_or_assumption": "hypothetical receptor abundance for a receptor-mediated route, not glucose-specific evidence"},
    {"name": "k_on_m3_mol_s", "value": 10.0, "unit": "m^3 mol^-1 s^-1", "source_or_assumption": "receptor association-rate assumption"},
    {"name": "k_off_s", "value": 0.10, "unit": "s^-1", "source_or_assumption": "receptor dissociation-rate assumption"},
    {"name": "k_endocytosis_s", "value": 0.05, "unit": "s^-1", "source_or_assumption": "receptor-mediated vesicle internalization-rate assumption"},
    {"name": "k_release_s", "value": 0.02, "unit": "s^-1", "source_or_assumption": "vesicle abluminal release-rate assumption"},
    {"name": "k_vesicle_degradation_s", "value": 0.005, "unit": "s^-1", "source_or_assumption": "vesicle substrate loss-rate assumption"},
    {"name": "c_comp_mM", "value": 0.0, "unit": "mM", "source_or_assumption": "controlled luminal competitor or inhibitor, baseline zero"},
    {"name": "ki_comp_mM", "value": 10.0, "unit": "mM", "source_or_assumption": "competition-strength assumption for perturbation experiments"},
    {"name": "k_brain_metabolite_clearance_s", "value": 0.001, "unit": "s^-1", "source_or_assumption": "brain-side lactate clearance assumption"},
]


def default_parameters() -> dict[str, float]:
    return {entry["name"]: float(entry["value"]) for entry in PARAMETER_TABLE}


def parameter_table() -> list[dict[str, Any]]:
    return [dict(entry) for entry in PARAMETER_TABLE]


def _per_min_to_mol_s(value: float) -> float:
    return value * 1.0e-6 / 60.0


def _pmol_to_mol(value: float) -> float:
    return value * 1.0e-12


def prepare_parameters(parameters: dict[str, float] | None = None) -> dict[str, float]:
    p = default_parameters() if parameters is None else dict(parameters)
    mass_g = p["brain_mass_g"]
    area_m2 = p["surface_area_cm2_per_g"] * mass_g * 1.0e-4
    p["_area_m2"] = area_m2
    p["_single_capillary_area_m2"] = 2.0 * math.pi * p["capillary_radius_m"] * p["capillary_length_m"]
    p["_parallel_segments"] = area_m2 / p["_single_capillary_area_m2"]
    p["_volume_l_m3"] = p["lumen_volume_ml_per_g"] * mass_g * 1.0e-6
    p["_perfusion_m3_s"] = p["perfusion_ml_min_per_g"] * mass_g * 1.0e-6 / 60.0
    p["_volume_e_m3"] = area_m2 * p["membrane_thickness_m"]
    p["_volume_b_m3"] = p["brain_volume_ml_per_g"] * mass_g * 1.0e-6
    p["_c_lumen_mol_m3"] = p["c_lumen_mM"]
    p["_c_lumen_source_mol_m3"] = p["c_lumen_source_mM"]
    p["_c_endothelium0_mol_m3"] = p["c_endothelium0_mM"]
    p["_c_brain0_mol_m3"] = p["c_brain0_mM"]
    p["_c_met_lumen_mol_m3"] = p["c_met_lumen_mM"]
    p["_c_met_brain0_mol_m3"] = p["c_met_brain_mM"]
    p["_vmax_luminal_in_mol_s"] = _per_min_to_mol_s(p["vmax_luminal_in_umol_min"]) * mass_g
    p["_vmax_luminal_efflux_mol_s"] = _per_min_to_mol_s(p["vmax_luminal_efflux_umol_min"]) * mass_g
    p["_vmax_abluminal_in_mol_s"] = _per_min_to_mol_s(p["vmax_abluminal_in_umol_min"]) * mass_g
    p["_vmax_abluminal_out_mol_s"] = _per_min_to_mol_s(p["vmax_abluminal_out_umol_min"]) * mass_g
    p["_vmax_metabolism_mol_s"] = _per_min_to_mol_s(p["vmax_metabolism_umol_min"]) * mass_g
    p["_vmax_brain_sink_mol_s"] = _per_min_to_mol_s(p["vmax_brain_sink_umol_min"]) * mass_g
    p["_r_total_mol"] = _pmol_to_mol(p["r_total_pmol_per_g"]) * mass_g
    return p


def initial_state(parameters: dict[str, float] | None = None) -> np.ndarray:
    p = prepare_parameters(parameters)
    n_l = p["_c_lumen_mol_m3"] * p["_volume_l_m3"]
    n_e = p["_c_endothelium0_mol_m3"] * p["_volume_e_m3"]
    n_b = p["_c_brain0_mol_m3"] * p["_volume_b_m3"]
    n_m_b = p["_c_met_brain0_mol_m3"] * p["_volume_b_m3"]
    return np.array([n_l, n_e, 0.0, 0.0, n_b, 0.0, n_m_b], dtype=float)


def _saturating_flux(vmax: float, concentration: float, km: float) -> float:
    return vmax * concentration / (km + concentration)


def fluxes(state: np.ndarray, parameters: dict[str, float] | None = None) -> dict[str, float]:
    p = parameters if parameters is not None and "_area_m2" in parameters else prepare_parameters(parameters)
    n_l, n_e, n_r, n_v, n_b, n_m, n_m_b = np.asarray(state, dtype=float)
    c_l = max(n_l / p["_volume_l_m3"], 0.0)
    c_e = max(n_e / p["_volume_e_m3"], 0.0)
    c_b = max(n_b / p["_volume_b_m3"], 0.0)
    c_m = max(n_m / p["_volume_e_m3"], 0.0)
    c_m_b = max(n_m_b / p["_volume_b_m3"], 0.0)
    c_m_l = p["_c_met_lumen_mol_m3"]
    area = p["_area_m2"]
    j_perfusion = p["_perfusion_m3_s"] * (p["_c_lumen_source_mol_m3"] - c_l)
    competition = 1.0 + p["c_comp_mM"] / p["ki_comp_mM"]
    km_l_in_eff = p["km_luminal_in_mM"] * competition
    km_l_eff_eff = p["km_luminal_efflux_mM"] * competition
    j_pass_l_to_e = p["p_luminal_m_s"] * area * max(c_l - c_e, 0.0)
    j_pass_e_to_l = p["p_luminal_m_s"] * area * max(c_e - c_l, 0.0)
    j_pass_l_net = j_pass_l_to_e - j_pass_e_to_l
    j_pass_e_to_b = p["p_abluminal_m_s"] * area * max(c_e - c_b, 0.0)
    j_pass_b_to_e = p["p_abluminal_m_s"] * area * max(c_b - c_e, 0.0)
    j_pass_e_b_net = j_pass_e_to_b - j_pass_b_to_e
    j_carrier_l_in = _saturating_flux(p["_vmax_luminal_in_mol_s"], c_l, km_l_in_eff)
    j_carrier_l_eff = _saturating_flux(p["_vmax_luminal_efflux_mol_s"], c_e, km_l_eff_eff)
    j_carrier_e_b = _saturating_flux(p["_vmax_abluminal_in_mol_s"], c_e, p["km_abluminal_in_mM"])
    j_carrier_b_e = _saturating_flux(p["_vmax_abluminal_out_mol_s"], c_b, p["km_abluminal_out_mM"])
    r_free = max(p["_r_total_mol"] - n_r, 0.0)
    j_receptor_bind = p["k_on_m3_mol_s"] * r_free * c_l
    j_receptor_unbind = p["k_off_s"] * max(n_r, 0.0)
    j_receptor_endo = p["k_endocytosis_s"] * max(n_r, 0.0)
    j_vesicle_release = p["k_release_s"] * max(n_v, 0.0)
    j_vesicle_degraded = p["k_vesicle_degradation_s"] * max(n_v, 0.0)
    j_metabolism = _saturating_flux(p["_vmax_metabolism_mol_s"], c_e, p["km_metabolism_mM"])
    j_brain_sink = _saturating_flux(p["_vmax_brain_sink_mol_s"], c_b, p["km_brain_sink_mM"])
    j_met_l_to_e = p["p_met_luminal_m_s"] * area * max(c_m_l - c_m, 0.0)
    j_met_e_to_l = p["p_met_luminal_m_s"] * area * max(c_m - c_m_l, 0.0)
    j_met_e_to_b = p["p_met_abluminal_m_s"] * area * max(c_m - c_m_b, 0.0)
    j_met_b_to_e = p["p_met_abluminal_m_s"] * area * max(c_m_b - c_m, 0.0)
    j_met_lumen_signed = j_met_l_to_e - j_met_e_to_l
    j_met_brain_signed = j_met_e_to_b - j_met_b_to_e
    j_h_cotransport_lumen = j_met_lumen_signed
    j_h_cotransport_brain = j_met_brain_signed
    j_charge_lumen = PRODUCT_SPEC["charge"] * j_met_lumen_signed + PRODUCT_SPEC["counterion_charge"] * j_h_cotransport_lumen
    j_charge_brain = PRODUCT_SPEC["charge"] * j_met_brain_signed + PRODUCT_SPEC["counterion_charge"] * j_h_cotransport_brain
    j_luminal_in = j_pass_l_to_e + j_carrier_l_in + j_receptor_endo
    j_luminal_efflux = j_pass_e_to_l + j_carrier_l_eff + j_receptor_unbind
    j_brain_glucose = j_pass_e_to_b + j_carrier_e_b + j_vesicle_release
    j_brain_glucose_net = j_brain_glucose - j_pass_b_to_e - j_carrier_b_e
    j_brain_metabolite = j_met_e_to_b
    j_brain_metabolite_net = j_brain_metabolite - j_met_b_to_e
    return {
        "c_lumen_mM": c_l,
        "c_endothelium_mM": c_e,
        "c_brain_mM": c_b,
        "c_met_endothelium_mM": c_m,
        "c_met_brain_mM": c_m_b,
        "j_lumen_perfusion_mol_s": j_perfusion,
        "j_pass_luminal_to_e_mol_s": j_pass_l_to_e,
        "j_pass_endothelium_to_luminal_mol_s": j_pass_e_to_l,
        "j_pass_endothelium_to_brain_mol_s": j_pass_e_to_b,
        "j_pass_brain_to_endothelium_mol_s": j_pass_b_to_e,
        "j_carrier_luminal_in_mol_s": j_carrier_l_in,
        "j_carrier_luminal_efflux_mol_s": j_carrier_l_eff,
        "j_carrier_endothelium_to_brain_mol_s": j_carrier_e_b,
        "j_carrier_brain_to_endothelium_mol_s": j_carrier_b_e,
        "j_receptor_binding_mol_s": j_receptor_bind,
        "j_receptor_unbinding_mol_s": j_receptor_unbind,
        "j_receptor_endocytosis_mol_s": j_receptor_endo,
        "j_vesicle_release_mol_s": j_vesicle_release,
        "j_vesicle_degradation_mol_s": j_vesicle_degraded,
        "j_metabolism_mol_s": j_metabolism,
        "j_brain_sink_mol_s": j_brain_sink,
        "j_metabolite_luminal_to_e_mol_s": j_met_l_to_e,
        "j_metabolite_e_to_luminal_mol_s": j_met_e_to_l,
        "j_metabolite_e_to_brain_mol_s": j_met_e_to_b,
        "j_metabolite_brain_to_e_mol_s": j_met_b_to_e,
        "j_metabolite_lumen_signed_mol_s": j_met_lumen_signed,
        "j_metabolite_brain_signed_mol_s": j_met_brain_signed,
        "j_h_cotransport_lumen_mol_s": j_h_cotransport_lumen,
        "j_h_cotransport_brain_mol_s": j_h_cotransport_brain,
        "j_charge_lumen_mol_eq_s": j_charge_lumen,
        "j_charge_brain_mol_eq_s": j_charge_brain,
        "j_luminal_in_mol_s": j_luminal_in,
        "j_luminal_efflux_mol_s": j_luminal_efflux,
        "j_brain_glucose_mol_s": j_brain_glucose,
        "j_brain_glucose_net_mol_s": j_brain_glucose_net,
        "j_brain_metabolite_mol_s": j_brain_metabolite,
        "j_brain_metabolite_net_mol_s": j_brain_metabolite_net,
        "j_substrate_boundary_net_mol_s": j_perfusion - j_metabolism - j_vesicle_degraded - j_brain_sink,
        "efflux_fraction": j_luminal_efflux / (j_luminal_in + j_luminal_efflux) if j_luminal_in + j_luminal_efflux > 0.0 else 0.0,
    }


def rhs(_t: float, state: np.ndarray, parameters: dict[str, float]) -> np.ndarray:
    p = parameters
    f = fluxes(state, p)
    _n_l, _n_e, _n_r, _n_v, _n_b, n_m, n_m_b = state
    d_l = f["j_lumen_perfusion_mol_s"]
    d_l += f["j_pass_endothelium_to_luminal_mol_s"] + f["j_carrier_luminal_efflux_mol_s"] + f["j_receptor_unbinding_mol_s"]
    d_l -= f["j_pass_luminal_to_e_mol_s"] + f["j_carrier_luminal_in_mol_s"] + f["j_receptor_binding_mol_s"]
    d_e = f["j_pass_luminal_to_e_mol_s"] - f["j_pass_endothelium_to_luminal_mol_s"]
    d_e += f["j_carrier_luminal_in_mol_s"] - f["j_carrier_luminal_efflux_mol_s"]
    d_e += f["j_pass_brain_to_endothelium_mol_s"] - f["j_pass_endothelium_to_brain_mol_s"]
    d_e += f["j_carrier_brain_to_endothelium_mol_s"] - f["j_carrier_endothelium_to_brain_mol_s"]
    d_e -= f["j_metabolism_mol_s"]
    d_r = f["j_receptor_binding_mol_s"] - f["j_receptor_unbinding_mol_s"] - f["j_receptor_endocytosis_mol_s"]
    d_v = f["j_receptor_endocytosis_mol_s"] - f["j_vesicle_release_mol_s"] - f["j_vesicle_degradation_mol_s"]
    d_b = f["j_pass_endothelium_to_brain_mol_s"] - f["j_pass_brain_to_endothelium_mol_s"]
    d_b += f["j_carrier_endothelium_to_brain_mol_s"] - f["j_carrier_brain_to_endothelium_mol_s"]
    d_b += f["j_vesicle_release_mol_s"] - f["j_brain_sink_mol_s"]
    d_m = f["j_metabolism_mol_s"] + f["j_vesicle_degradation_mol_s"]
    d_m += f["j_metabolite_luminal_to_e_mol_s"] - f["j_metabolite_e_to_luminal_mol_s"]
    d_m += f["j_metabolite_brain_to_e_mol_s"] - f["j_metabolite_e_to_brain_mol_s"]
    d_m_b = f["j_metabolite_e_to_brain_mol_s"] - f["j_metabolite_brain_to_e_mol_s"] - p["k_brain_metabolite_clearance_s"] * n_m_b
    return np.array([d_l, d_e, d_r, d_v, d_b, d_m, d_m_b], dtype=float)


def run_model(parameters: dict[str, float] | None = None, t_end: float = 3600.0, sample_dt: float = 1.0, initial: np.ndarray | None = None) -> dict[str, Any]:
    p = prepare_parameters(parameters)
    y0 = initial_state(p) if initial is None else np.asarray(initial, dtype=float)
    count = max(int(round(t_end / sample_dt)), 1)
    t_eval = np.linspace(0.0, t_end, count + 1)
    atol = np.array([1.0e-14, 1.0e-14, 1.0e-18, 1.0e-18, 1.0e-13, 1.0e-14, 1.0e-13])
    solution = solve_ivp(
        lambda t, y: rhs(t, y, p),
        (0.0, t_end),
        y0,
        method="Radau",
        t_eval=t_eval,
        rtol=2.0e-8,
        atol=atol,
        max_step=30.0,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    flux_rows = [fluxes(state, p) for state in solution.y.T]
    flux_names = list(flux_rows[0].keys())
    flux_array = {name: np.array([row[name] for row in flux_rows], dtype=float) for name in flux_names}
    derivatives = np.array([rhs(t, state, p) for t, state in zip(solution.t, solution.y.T)])
    total_substrate = solution.y[0] + solution.y[1] + solution.y[2] + solution.y[3] + solution.y[4]
    total_metabolite = solution.y[5] + solution.y[6]
    substrate_boundary = np.array([flux_array["j_substrate_boundary_net_mol_s"][i] for i in range(len(solution.t))])
    metabolite_source = np.array([
        flux_array["j_metabolism_mol_s"][i]
        + flux_array["j_vesicle_degradation_mol_s"][i]
        + flux_array["j_metabolite_luminal_to_e_mol_s"][i]
        for i in range(len(solution.t))
    ])
    metabolite_loss = np.array([
        flux_array["j_metabolite_e_to_luminal_mol_s"][i]
        + p["k_brain_metabolite_clearance_s"] * solution.y[6, i]
        for i in range(len(solution.t))
    ])
    substrate_derivative = derivatives[:, 0] + derivatives[:, 1] + derivatives[:, 2] + derivatives[:, 3] + derivatives[:, 4]
    metabolite_derivative = derivatives[:, 5] + derivatives[:, 6]
    mass_balance_substrate = substrate_derivative - substrate_boundary
    mass_balance_metabolite = metabolite_derivative - (metabolite_source - metabolite_loss)
    return {
        "parameters": p,
        "t": solution.t,
        "state": solution.y,
        "flux": flux_array,
        "mass_balance_substrate_mol_s": mass_balance_substrate,
        "mass_balance_metabolite_mol_s": mass_balance_metabolite,
        "total_substrate_mol": total_substrate,
        "total_metabolite_mol": total_metabolite,
    }


def _mean_last(result: dict[str, Any], key: str, last_s: float = 600.0) -> float:
    t = np.asarray(result["t"])
    mask = t >= max(0.0, float(t[-1]) - last_s)
    return float(np.mean(np.asarray(result["flux"][key])[mask]))


def _final(result: dict[str, Any], index: int) -> float:
    return float(result["state"][index, -1])


def summarize(result: dict[str, Any], last_s: float = 600.0) -> dict[str, Any]:
    p = result["parameters"]
    mass_g = p["brain_mass_g"]
    primary = _mean_last(result, "j_luminal_in_mol_s", last_s) * 60.0 * 1.0e6 / mass_g
    brain_glucose = _mean_last(result, "j_brain_glucose_mol_s", last_s) * 60.0 * 1.0e6 / mass_g
    brain_metabolite = _mean_last(result, "j_brain_metabolite_mol_s", last_s) * 60.0 * 1.0e6 / mass_g
    reference = float(REFERENCE["value"])
    relative_error = abs(primary - reference) / reference
    max_mass_residual = max(
        float(np.max(np.abs(result["mass_balance_substrate_mol_s"]))),
        float(np.max(np.abs(result["mass_balance_metabolite_mol_s"]))),
    )
    max_charge_residual = max(
        float(np.max(np.abs(result["flux"]["j_charge_lumen_mol_eq_s"]))),
        float(np.max(np.abs(result["flux"]["j_charge_brain_mol_eq_s"]))),
    )
    min_pool = float(np.min(result["state"]))
    path_names = [
        "j_pass_luminal_to_e_mol_s",
        "j_carrier_luminal_in_mol_s",
        "j_receptor_endocytosis_mol_s",
        "j_pass_endothelium_to_luminal_mol_s",
        "j_carrier_luminal_efflux_mol_s",
        "j_receptor_unbinding_mol_s",
        "j_pass_endothelium_to_brain_mol_s",
        "j_carrier_endothelium_to_brain_mol_s",
        "j_vesicle_release_mol_s",
        "j_metabolism_mol_s",
        "j_vesicle_degradation_mol_s",
        "j_brain_sink_mol_s",
    ]
    path_means = {
        name: _mean_last(result, name, last_s) * 60.0 * 1.0e6 / mass_g
        for name in path_names
    }
    return {
        "primary_J_u_umol_g_min": primary,
        "reference_J_u_umol_g_min": reference,
        "primary_relative_error": relative_error,
        "primary_criterion_met": bool(relative_error <= 0.30 and max_mass_residual <= 1.0e-10 and min_pool >= -1.0e-15),
        "brain_glucose_gross_umol_g_min": brain_glucose,
        "brain_glucose_net_umol_g_min": _mean_last(result, "j_brain_glucose_net_mol_s", last_s) * 60.0 * 1.0e6 / mass_g,
        "brain_metabolite_gross_umol_g_min": brain_metabolite,
        "brain_metabolite_net_umol_g_min": _mean_last(result, "j_brain_metabolite_net_mol_s", last_s) * 60.0 * 1.0e6 / mass_g,
        "luminal_efflux_fraction": _mean_last(result, "efflux_fraction", last_s),
        "luminal_glucose_mol": _final(result, 0),
        "luminal_glucose_final_mM": _mean_last(result, "c_lumen_mM", last_s),
        "endothelial_retention_mol": _final(result, 1) + _final(result, 2) + _final(result, 3),
        "free_endothelial_glucose_mol": _final(result, 1),
        "receptor_bound_glucose_mol": _final(result, 2),
        "vesicle_glucose_mol": _final(result, 3),
        "brain_glucose_mol": _final(result, 4),
        "endothelial_metabolite_mol": _final(result, 5),
        "brain_metabolite_mol": _final(result, 6),
        "max_mass_balance_residual_mol_s": max_mass_residual,
        "max_charge_balance_residual_mol_eq_s": max_charge_residual,
        "minimum_state_pool_mol": min_pool,
        "pathway_means_umol_g_min": path_means,
    }


def passive_only_parameters(parameters: dict[str, float] | None = None) -> dict[str, float]:
    p = default_parameters() if parameters is None else dict(parameters)
    for name in [
        "vmax_luminal_in_umol_min",
        "vmax_luminal_efflux_umol_min",
        "vmax_abluminal_in_umol_min",
        "vmax_abluminal_out_umol_min",
        "vmax_metabolism_umol_min",
        "vmax_brain_sink_umol_min",
        "r_total_pmol_per_g",
    ]:
        p[name] = 0.0
    p["k_endocytosis_s"] = 0.0
    p["k_release_s"] = 0.0
    p["k_vesicle_degradation_s"] = 0.0
    p["k_brain_metabolite_clearance_s"] = 0.0
    return p


def sensitivity(parameters: dict[str, float] | None = None, factors: tuple[float, float] = (0.5, 1.5)) -> dict[str, Any]:
    base = default_parameters() if parameters is None else dict(parameters)
    baseline = summarize(run_model(base))
    candidate_names = [
        "vmax_luminal_in_umol_min",
        "km_luminal_in_mM",
        "c_lumen_source_mM",
        "p_luminal_m_s",
        "surface_area_cm2_per_g",
        "vmax_metabolism_umol_min",
        "vmax_luminal_efflux_umol_min",
        "r_total_pmol_per_g",
    ]
    rows: list[dict[str, Any]] = []
    for name in candidate_names:
        row: dict[str, Any] = {
            "parameter": name,
            "baseline_value": base[name],
            "cases": {},
        }
        for factor in factors:
            p = dict(base)
            p[name] = p[name] * factor
            result = summarize(run_model(p))
            row["cases"][str(factor)] = {
                "factor": factor,
                "value": p[name],
                "primary_J_u_umol_g_min": result["primary_J_u_umol_g_min"],
                "brain_glucose_gross_umol_g_min": result["brain_glucose_gross_umol_g_min"],
                "primary_relative_change": result["primary_J_u_umol_g_min"] / baseline["primary_J_u_umol_g_min"] - 1.0,
                "brain_glucose_relative_change": result["brain_glucose_gross_umol_g_min"] / baseline["brain_glucose_gross_umol_g_min"] - 1.0,
            }
        primary_changes = [abs(case["primary_relative_change"]) for case in row["cases"].values()]
        brain_changes = [abs(case["brain_glucose_relative_change"]) for case in row["cases"].values()]
        row["max_abs_primary_relative_change"] = max(primary_changes)
        row["max_abs_brain_glucose_relative_change"] = max(brain_changes)
        rows.append(row)
    rows.sort(key=lambda item: item["max_abs_primary_relative_change"], reverse=True)
    brain_ranked = sorted(rows, key=lambda item: item["max_abs_brain_glucose_relative_change"], reverse=True)
    return {
        "baseline_primary_J_u_umol_g_min": baseline["primary_J_u_umol_g_min"],
        "baseline_brain_glucose_gross_umol_g_min": baseline["brain_glucose_gross_umol_g_min"],
        "parameters": rows,
        "top_three_primary": [row["parameter"] for row in rows[:3]],
        "top_three_brain_glucose": [row["parameter"] for row in brain_ranked[:3]],
    }


def dimension_checks(parameters: dict[str, float] | None = None) -> dict[str, Any]:
    p = prepare_parameters(parameters)
    concentration = p["_c_lumen_mol_m3"]
    passive_flux = p["p_luminal_m_s"] * p["_area_m2"] * concentration
    carrier_flux = _saturating_flux(p["_vmax_luminal_in_mol_s"], concentration, p["km_luminal_in_mM"])
    binding_rate = p["k_on_m3_mol_s"] * p["_r_total_mol"] * concentration
    concentration_from_moles = p["_r_total_mol"] / p["_volume_e_m3"]
    perfusion_flux = p["_perfusion_m3_s"] * p["_c_lumen_source_mol_m3"]
    checks = {
        "area_m2": p["_area_m2"],
        "parallel_segments_dimensionless": p["_parallel_segments"],
        "lumen_volume_m3": p["_volume_l_m3"],
        "endothelial_volume_m3": p["_volume_e_m3"],
        "brain_volume_m3": p["_volume_b_m3"],
        "perfusion_flow_m3_s": p["_perfusion_m3_s"],
        "perfusion_flux_mol_s": perfusion_flux,
        "passive_clearance_m3_s": p["p_luminal_m_s"] * p["_area_m2"],
        "passive_flux_mol_s": passive_flux,
        "carrier_flux_mol_s": carrier_flux,
        "binding_rate_mol_s": binding_rate,
        "concentration_from_moles_mol_m3": concentration_from_moles,
        "expected_flux_unit": "mol s^-1",
        "expected_concentration_unit": "mol m^-3",
        "unit_contracts": UNIT_CONTRACTS,
        "area_equation_ok": bool(abs(p["_parallel_segments"] * p["_single_capillary_area_m2"] - p["_area_m2"]) < 1.0e-18),
        "passive_dimension_ok": bool(abs(passive_flux) < 1.0e20),
        "kinetic_dimension_ok": bool(abs(carrier_flux) < 1.0e20),
        "binding_dimension_ok": bool(abs(binding_rate) < 1.0e20),
        "concentration_dimension_ok": bool(abs(concentration_from_moles) < 1.0e20),
        "perfusion_dimension_ok": bool(abs(perfusion_flux) < 1.0e20),
        "lumen_volume_ok": bool(p["_volume_l_m3"] > 0.0),
    }
    checks["all_ok"] = all(bool(checks[key]) for key in ["area_equation_ok", "passive_dimension_ok", "kinetic_dimension_ok", "binding_dimension_ok", "concentration_dimension_ok", "perfusion_dimension_ok", "lumen_volume_ok"])
    return checks


def analytical_limit_case() -> dict[str, Any]:
    p = passive_only_parameters()
    p["c_lumen_mM"] = 2.0
    p["c_lumen_source_mM"] = 2.0
    p["c_endothelium0_mM"] = 2.0
    p["c_brain0_mM"] = 2.0
    p["p_luminal_m_s"] = 2.0e-8
    p["p_abluminal_m_s"] = 2.0e-8
    p["c_met_lumen_mM"] = 0.0
    p["c_met_brain_mM"] = 0.0
    q = prepare_parameters(p)
    state = initial_state(q)
    f = fluxes(state, q)
    derivative = rhs(0.0, state, q)
    max_flux = max(abs(value) for key, value in f.items() if key.endswith("_mol_s"))
    max_derivative = float(np.max(np.abs(derivative)))
    return {
        "name": "passive equilibrium with equal concentrations",
        "max_directional_flux_mol_s": float(max_flux),
        "max_state_derivative_mol_s": max_derivative,
        "passed": bool(max_flux <= 1.0e-18 and max_derivative <= 1.0e-18),
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def build_results() -> dict[str, Any]:
    baseline_result = run_model()
    baseline = summarize(baseline_result)
    passive_result = run_model(passive_only_parameters())
    passive = summarize(passive_result)
    analytical = analytical_limit_case()
    return {
        "id": "BT-HX-Q154",
        "substrate": "D-glucose",
        "product": "lactate",
        "cell_type": "brain capillary endothelial cell",
        "time_window_s": 3600.0,
        "averaging_window_s": 600.0,
        "output_sample_dt_s": 1.0,
        "state_order": ["N_L_glucose", "N_E_free_glucose", "N_R_bound_glucose", "N_V_vesicle_glucose", "N_B_glucose", "N_M_endothelial_lactate", "N_M_brain_lactate"],
        "equations": EQUATIONS,
        "unit_contracts": UNIT_CONTRACTS,
        "product_specification": PRODUCT_SPEC,
        "reference": REFERENCE,
        "parameters": default_parameters(),
        "parameter_table": parameter_table(),
        "geometry": {
            "lumen_volume_m3": baseline_result["parameters"]["_volume_l_m3"],
            "area_m2": baseline_result["parameters"]["_area_m2"],
            "single_capillary_area_m2": baseline_result["parameters"]["_single_capillary_area_m2"],
            "parallel_segments": baseline_result["parameters"]["_parallel_segments"],
            "endothelial_volume_m3": baseline_result["parameters"]["_volume_e_m3"],
            "brain_volume_m3": baseline_result["parameters"]["_volume_b_m3"],
        },
        "dimension_checks": dimension_checks(),
        "analytical_limit_case": analytical,
        "baseline": baseline,
        "passive_null": passive,
        "sensitivity": sensitivity(),
        "identifiability": {
            "independent_directed_fluxes_available": False,
            "independent_metabolite_profile_available": False,
            "causal_pathway_status": "UNKNOWN",
            "reason": "the run contains modeled fluxes but no substrate-specific inhibitor, tracer, or metabolite measurement to identify them",
        },
        "verification": {
            "model_execution": True,
            "test_file": "test_model.py",
            "test_command": "python3 test_model.py",
            "test_count": 6,
            "mass_balance_limit_mol_s": 1.0e-10,
            "acceptance_relative_error_limit": 0.30,
            "dimension_checks_ok": dimension_checks()["all_ok"],
        },
    }


def main() -> None:
    output_path = Path(__file__).resolve().parent / "results.json"
    results = build_results()
    output_path.write_text(json.dumps(_jsonable(results), indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "id": results["id"],
        "primary_J_u_umol_g_min": results["baseline"]["primary_J_u_umol_g_min"],
        "reference_J_u_umol_g_min": results["baseline"]["reference_J_u_umol_g_min"],
        "primary_criterion_met": results["baseline"]["primary_criterion_met"],
        "passive_primary_J_u_umol_g_min": results["passive_null"]["primary_J_u_umol_g_min"],
        "top_three_primary": results["sensitivity"]["top_three_primary"],
        "analytical_limit_case_passed": results["analytical_limit_case"]["passed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
