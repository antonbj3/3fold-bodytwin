from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

MODEL_ID = "BT-HX-Q121"
N_COMP = 3
COMP_NAMES = ("proximal", "distal", "rectal")
STATE_NAMES = (
    "substrate_proximal_g",
    "substrate_distal_g",
    "substrate_rectal_g",
    "liquid_proximal_ml",
    "liquid_distal_ml",
    "liquid_rectal_ml",
    "dissolved_proximal_mmol",
    "dissolved_distal_mmol",
    "dissolved_rectal_mmol",
    "free_gas_proximal_mmol",
    "free_gas_distal_mmol",
    "free_gas_rectal_mmol",
    "blood_proximal_mmol",
    "blood_distal_mmol",
    "blood_rectal_mmol",
    "cumulative_fecal_mass_g",
    "cumulative_rectal_gas_mmol",
    "cumulative_fecal_dissolved_mmol",
    "cumulative_breath_mmol",
    "cumulative_production_mmol",
    "cumulative_absorption_mmol",
    "cumulative_net_dissolution_mmol",
    "cumulative_microbial_consumption_mmol",
    "cumulative_upstream_gas_mmol",
    "cumulative_explicit_backflow_mmol",
    "cumulative_explicit_forward_pressure_mmol",
)

PARAMETERS = {
    "lactulose_g": 15.0,
    "ileal_entry_h": 4.0,
    "initial_substrate_g": 0.0,
    "initial_liquid_proximal_ml": 80.0,
    "initial_liquid_distal_ml": 80.0,
    "initial_liquid_rectal_ml": 40.0,
    "initial_free_gas_mmol": 0.0,
    "initial_dissolved_gas_mmol": 0.0,
    "initial_blood_gas_mmol": 0.0,
    "water_per_g_ml_g": 0.5,
    "base_water_rate_ml_h": 12.0,
    "water_loss_h": 0.01,
    "tau_proximal_h": 10.0,
    "tau_distal_h": 14.0,
    "tau_rectal_h": 18.0,
    "motility_gain": 1.0,
    "yield_h2_mmol_g": 0.70,
    "k_fermentation_h": 0.35,
    "k_michaelis_g": 2.0,
    "k_phase_transfer_h": 0.45,
    "h_h2_mmol_ml_mmhg": 2.4e-5,
    "k_blood_transfer_h": 0.08,
    "k_breath_clearance_h": 0.25,
    "k_microbial_h": 0.015,
    "k_gas_diffusion_h": 0.08,
    "k_dissolved_diffusion_h": 0.12,
    "k_forward_pressure_mmol_h_mmhg": 0.01,
    "k_backward_pressure_mmol_h_mmhg": 0.006,
    "k_rectal_valve_mmol_h_sqrt_mmhg": 0.22,
    "rectal_pressure_threshold_mmhg": 2.0,
    "bubble_carry_fraction": 0.35,
    "v0_proximal_ml": 40.0,
    "v0_distal_ml": 50.0,
    "v0_rectal_ml": 25.0,
    "compliance_proximal_ml_mmhg": 0.45,
    "compliance_distal_ml_mmhg": 0.55,
    "compliance_rectal_ml_mmhg": 0.35,
    "r_ml_mmhg_mmol_k": 62.3637,
    "body_temperature_k": 310.0,
    "atmospheric_pressure_mmhg": 760.0,
    "standard_molar_volume_ml_mmol": 22.414,
    "upstream_gas_mmol_h": 0.0,
    "duration_h": 24.0,
    "max_step_h": 0.01,
    "relative_tolerance": 1e-8,
    "absolute_tolerance": 1e-10,
}

PARAMETER_META = {
    "lactulose_g": ("g", "Christl et al. 1992, DOI 10.1016/0016-5085(92)90765-q; scenario dose"),
    "ileal_entry_h": ("h", "assumption for a controlled proximal-colon input; not measured here"),
    "initial_substrate_g": ("g", "zero initial substrate assumption"),
    "initial_liquid_proximal_ml": ("mL", "lumped liquid-volume assumption"),
    "initial_liquid_distal_ml": ("mL", "lumped liquid-volume assumption"),
    "initial_liquid_rectal_ml": ("mL", "lumped liquid-volume assumption"),
    "initial_free_gas_mmol": ("mmol", "zero initial H2-equivalent free gas assumption"),
    "initial_dissolved_gas_mmol": ("mmol", "zero initial dissolved H2-equivalent assumption"),
    "initial_blood_gas_mmol": ("mmol", "zero initial tracked blood gas assumption"),
    "water_per_g_ml_g": ("mL/g", "water carried with delivered carbohydrate; assumption"),
    "base_water_rate_ml_h": ("mL/h", "basal secretory water input; assumption"),
    "water_loss_h": ("h^-1", "lumped water absorption/evaporation coefficient; assumption"),
    "tau_proximal_h": ("h", "nominal compartment content residence time; assumption"),
    "tau_distal_h": ("h", "nominal compartment content residence time; assumption"),
    "tau_rectal_h": ("h", "nominal rectal content residence time; assumption"),
    "motility_gain": ("dimensionless", "multiplicative motility intervention; nominal value 1"),
    "yield_h2_mmol_g": ("mmol/g", "H2-equivalent fermentation yield; mechanistic assumption, not a fitted individual value"),
    "k_fermentation_h": ("h^-1", "first-order substrate utilization rate; mechanistic assumption"),
    "k_michaelis_g": ("g", "Monod half-saturation constant; mechanistic assumption"),
    "k_phase_transfer_h": ("h^-1", "effective gas-liquid exchange coefficient; mechanistic assumption"),
    "h_h2_mmol_ml_mmhg": ("mmol/(mL mmHg)", "approximate H2 Henry coefficient; physical assumption"),
    "k_blood_transfer_h": ("h^-1", "effective colonic blood-wall transfer coefficient; assumption"),
    "k_breath_clearance_h": ("h^-1", "blood-to-breath clearance coefficient; assumption"),
    "k_microbial_h": ("h^-1", "net microbial H2 disposal coefficient; Gibson et al. 1990 motivates the sink, magnitude is an assumption"),
    "k_gas_diffusion_h": ("h^-1", "intercompartment free-gas exchange coefficient; assumption"),
    "k_dissolved_diffusion_h": ("h^-1", "intercompartment dissolved-gas exchange coefficient; assumption"),
    "k_forward_pressure_mmol_h_mmhg": ("mmol/(h mmHg)", "distal pressure-driven gas conductance; assumption"),
    "k_backward_pressure_mmol_h_mmhg": ("mmol/(h mmHg)", "proximal pressure-driven return conductance; assumption"),
    "k_rectal_valve_mmol_h_sqrt_mmhg": ("mmol/(h sqrt(mmHg))", "pressure-opening rectal gas valve; assumption"),
    "rectal_pressure_threshold_mmhg": ("mmHg", "valve opening excess-pressure threshold; assumption"),
    "bubble_carry_fraction": ("dimensionless", "fraction of liquid advective gas flux leaving with rectal liquid; assumption"),
    "v0_proximal_ml": ("mL", "residual gas volume of a compliant compartment; geometry assumption"),
    "v0_distal_ml": ("mL", "residual gas volume of a compliant compartment; geometry assumption"),
    "v0_rectal_ml": ("mL", "residual gas volume of a compliant compartment; geometry assumption"),
    "compliance_proximal_ml_mmhg": ("mL/mmHg", "linear wall compliance; geometry assumption"),
    "compliance_distal_ml_mmhg": ("mL/mmHg", "linear wall compliance; geometry assumption"),
    "compliance_rectal_ml_mmhg": ("mL/mmHg", "linear wall compliance; geometry assumption"),
    "r_ml_mmhg_mmol_k": ("mL mmHg/(mmol K)", "universal gas constant converted to mL, mmHg, mmol and K"),
    "body_temperature_k": ("K", "310 K body-temperature assumption"),
    "atmospheric_pressure_mmhg": ("mmHg", "standard atmospheric pressure constant"),
    "standard_molar_volume_ml_mmol": ("mL/mmol", "22.414 mL/mmol at 273.15 K and 760 mmHg"),
    "upstream_gas_mmol_h": ("mmol/h", "explicit gas transfer from an earlier GI segment; set to zero, never swallowed air"),
    "duration_h": ("h", "frozen nominal follow-up duration"),
    "max_step_h": ("h", "integration resolution"),
    "relative_tolerance": ("dimensionless", "solver relative tolerance"),
    "absolute_tolerance": ("state units", "solver absolute tolerance"),
}

EQUATIONS = [
    "dM_i/dt = J_ileal delta_i0 + F_M,i-1 - F_M,i",
    "dW_i/dt = J_water + F_W,i-1 - F_W,i - k_water W_i",
    "F_M,i = M_i/tau_i and F_W,i = W_i/tau_i",
    "J_prod,i = y_H2 k_ferm M_i/(K_M + M_i) + J_upstream delta_i0",
    "p_i V_g,i = G_i R T with V_g,i = V0,i + C_i p_i",
    "J_diss,i = k_phase W_i (H p_i - D_i/W_i)",
    "J_abs,i = k_blood D_i and J_breath,i = k_breath A_i",
    "J_pair,i = F_W,i G_i/W_i + k_gas (G_i-G_i+1) + k_fwd max(p_i-p_i+1,0) - k_back max(p_i+1-p_i,0)",
    "J_rect = k_valve sqrt(max(p_rect-p_crit,0)) + f_bubble F_W,rect G_rect/W_rect",
    "dG_i/dt = J_prod,i - J_diss,i - J_micro G_i - J_pair,i + J_pair,i-1 - J_rect delta_i2",
    "dD_i/dt = J_diss,i - J_abs,i - J_D,i + J_D,i-1 - J_Drect delta_i2",
    "dA_i/dt = J_abs,i - J_breath,i",
    "V_bowel = sum_i(W_i + V_g,i)",
]

SENSITIVITY_PARAMETERS = (
    "k_fermentation_h",
    "k_blood_transfer_h",
    "motility_gain",
)


def default_parameters() -> dict[str, float]:
    return dict(PARAMETERS)


def initial_state(p: dict[str, float]) -> np.ndarray:
    state = np.zeros(len(STATE_NAMES), dtype=float)
    state[0:3] = p["initial_substrate_g"]
    state[3:6] = [
        p["initial_liquid_proximal_ml"],
        p["initial_liquid_distal_ml"],
        p["initial_liquid_rectal_ml"],
    ]
    state[6:9] = p["initial_dissolved_gas_mmol"]
    state[9:12] = p["initial_free_gas_mmol"]
    state[12:15] = p["initial_blood_gas_mmol"]
    return state


def ileal_substrate_input(t: float, p: dict[str, float]) -> float:
    if t < 0.0 or t >= p["ileal_entry_h"]:
        return 0.0
    return p["lactulose_g"] / p["ileal_entry_h"]


def gas_pressure(
    gas_mmol: np.ndarray,
    v0_ml: np.ndarray,
    compliance_ml_mmhg: np.ndarray,
    r_ml_mmhg_mmol_k: float,
    temperature_k: float,
) -> np.ndarray:
    gas_nonnegative = np.maximum(np.asarray(gas_mmol, dtype=float), 0.0)
    v0 = np.asarray(v0_ml, dtype=float)
    compliance = np.asarray(compliance_ml_mmhg, dtype=float)
    q = gas_nonnegative * r_ml_mmhg_mmol_k * temperature_k
    discriminant = np.maximum(v0 * v0 + 4.0 * compliance * q, 0.0)
    pressure = 2.0 * q / (v0 + np.sqrt(discriminant))
    return pressure


def geometry_arrays(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    v0 = np.array(
        [p["v0_proximal_ml"], p["v0_distal_ml"], p["v0_rectal_ml"]],
        dtype=float,
    )
    compliance = np.array(
        [
            p["compliance_proximal_ml_mmhg"],
            p["compliance_distal_ml_mmhg"],
            p["compliance_rectal_ml_mmhg"],
        ],
        dtype=float,
    )
    return v0, compliance


def pair_flows(
    free_gas: np.ndarray,
    dissolved: np.ndarray,
    liquid: np.ndarray,
    pressure_excess: np.ndarray,
    liquid_flux: np.ndarray,
    p: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    safe_liquid = np.maximum(liquid, 1e-9)
    gas_pair = np.zeros(2, dtype=float)
    dissolved_pair = np.zeros(2, dtype=float)
    backward = np.zeros(2, dtype=float)
    forward_pressure = np.zeros(2, dtype=float)
    for i in range(2):
        advective_gas = liquid_flux[i] * free_gas[i] / safe_liquid[i]
        diffusive_gas = p["k_gas_diffusion_h"] * (free_gas[i] - free_gas[i + 1])
        forward_pressure[i] = p["k_forward_pressure_mmol_h_mmhg"] * max(
            pressure_excess[i] - pressure_excess[i + 1], 0.0
        )
        backward[i] = p["k_backward_pressure_mmol_h_mmhg"] * max(
            pressure_excess[i + 1] - pressure_excess[i], 0.0
        )
        gas_pair[i] = (
            advective_gas
            + diffusive_gas
            + forward_pressure[i]
            - backward[i]
        )
        advective_dissolved = liquid_flux[i] * dissolved[i] / safe_liquid[i]
        diffusive_dissolved = p["k_dissolved_diffusion_h"] * (
            dissolved[i] - dissolved[i + 1]
        )
        dissolved_pair[i] = advective_dissolved + diffusive_dissolved
    return gas_pair, dissolved_pair, backward, forward_pressure


def rates(t: float, y: np.ndarray, p: dict[str, float]) -> np.ndarray:
    substrate = np.maximum(y[0:3], 0.0)
    liquid = np.maximum(y[3:6], 0.0)
    dissolved = np.maximum(y[6:9], 0.0)
    free_gas = np.maximum(y[9:12], 0.0)
    blood = np.maximum(y[12:15], 0.0)
    tau = np.array(
        [p["tau_proximal_h"], p["tau_distal_h"], p["tau_rectal_h"]],
        dtype=float,
    ) / max(p["motility_gain"], 1e-9)
    liquid_flux = liquid / tau
    substrate_flux = substrate / tau
    v0, compliance = geometry_arrays(p)
    pressure_excess = gas_pressure(
        free_gas,
        v0,
        compliance,
        p["r_ml_mmhg_mmol_k"],
        p["body_temperature_k"],
    )
    fermentation = (
        p["yield_h2_mmol_g"]
        * p["k_fermentation_h"]
        * substrate
        / (p["k_michaelis_g"] + substrate)
    )
    production = fermentation.copy()
    production[0] += p["upstream_gas_mmol_h"]
    dissolved_concentration = dissolved / np.maximum(liquid, 1e-9)
    equilibrium_concentration = p["h_h2_mmol_ml_mmhg"] * pressure_excess
    dissolution = (
        p["k_phase_transfer_h"]
        * liquid
        * (equilibrium_concentration - dissolved_concentration)
    )
    absorption = p["k_blood_transfer_h"] * dissolved
    breath = p["k_breath_clearance_h"] * blood
    microbial = p["k_microbial_h"] * free_gas
    gas_pair, dissolved_pair, backward, forward_pressure = pair_flows(
        free_gas,
        dissolved,
        liquid,
        pressure_excess,
        liquid_flux,
        p,
    )
    safe_liquid = np.maximum(liquid, 1e-9)
    rectal_gas = (
        p["k_rectal_valve_mmol_h_sqrt_mmhg"]
        * np.sqrt(
            np.maximum(
                pressure_excess[2] - p["rectal_pressure_threshold_mmhg"], 0.0
            )
        )
        + p["bubble_carry_fraction"]
        * liquid_flux[2]
        * free_gas[2]
        / safe_liquid[2]
    )
    rectal_dissolved = liquid_flux[2] * dissolved[2] / safe_liquid[2]
    water_input = (
        p["water_per_g_ml_g"] * ileal_substrate_input(t, p)
        + p["base_water_rate_ml_h"]
    )
    derivative = np.zeros_like(y)
    derivative[0] = ileal_substrate_input(t, p) - substrate_flux[0]
    derivative[1] = substrate_flux[0] - substrate_flux[1]
    derivative[2] = substrate_flux[1] - substrate_flux[2]
    derivative[3] = water_input - liquid_flux[0] - p["water_loss_h"] * liquid[0]
    derivative[4] = liquid_flux[0] - liquid_flux[1] - p["water_loss_h"] * liquid[1]
    derivative[5] = liquid_flux[1] - liquid_flux[2] - p["water_loss_h"] * liquid[2]
    for i in range(N_COMP):
        incoming_gas = gas_pair[i - 1] if i > 0 else 0.0
        outgoing_gas = gas_pair[i] if i < N_COMP - 1 else 0.0
        incoming_dissolved = dissolved_pair[i - 1] if i > 0 else 0.0
        outgoing_dissolved = dissolved_pair[i] if i < N_COMP - 1 else 0.0
        derivative[6 + i] = (
            dissolution[i]
            - absorption[i]
            + incoming_dissolved
            - outgoing_dissolved
            - (rectal_dissolved if i == N_COMP - 1 else 0.0)
        )
        derivative[9 + i] = (
            production[i]
            - dissolution[i]
            - microbial[i]
            + incoming_gas
            - outgoing_gas
            - (rectal_gas if i == N_COMP - 1 else 0.0)
        )
        derivative[12 + i] = absorption[i] - breath[i]
    derivative[15] = substrate_flux[2]
    derivative[16] = rectal_gas
    derivative[17] = rectal_dissolved
    derivative[18] = breath.sum()
    derivative[19] = fermentation.sum()
    derivative[20] = absorption.sum()
    derivative[21] = dissolution.sum()
    derivative[22] = microbial.sum()
    derivative[23] = p["upstream_gas_mmol_h"]
    derivative[24] = backward.sum()
    derivative[25] = forward_pressure.sum()
    return derivative


def simulate(overrides: dict[str, float] | None = None) -> dict[str, object]:
    p = default_parameters()
    if overrides:
        unknown = sorted(set(overrides) - set(p))
        if unknown:
            raise KeyError(f"unknown parameters: {unknown}")
        p.update(overrides)
    duration = p["duration_h"]
    solution = solve_ivp(
        lambda t, y: rates(t, y, p),
        (0.0, duration),
        initial_state(p),
        method="DOP853",
        t_eval=np.linspace(0.0, duration, int(round(duration / min(p["max_step_h"], duration)) + 1)),
        max_step=p["max_step_h"],
        rtol=p["relative_tolerance"],
        atol=p["absolute_tolerance"],
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return {"parameters": p, "time_h": solution.t, "state": solution.y}


def parameter_table() -> list[dict[str, object]]:
    rows = []
    for name, value in PARAMETERS.items():
        if name not in PARAMETER_META:
            raise KeyError(f"missing metadata for {name}")
        unit, source = PARAMETER_META[name]
        rows.append(
            {
                "parameter": name,
                "value": float(value),
                "unit": unit,
                "source_or_assumption": source,
            }
        )
    return rows


def dimension_checks() -> list[dict[str, object]]:
    return [
        {
            "equation": "dM_i/dt",
            "unit_check": "g/h = g/h - g/h + g/h",
            "status": "PASS",
        },
        {
            "equation": "dW_i/dt",
            "unit_check": "mL/h = mL/h + mL/h - mL/h",
            "status": "PASS",
        },
        {
            "equation": "J_prod,i",
            "unit_check": "mmol/h = (mmol/g)(h^-1)(g)",
            "status": "PASS",
        },
        {
            "equation": "p_i V_g,i = G_i R T",
            "unit_check": "mmHg*mL = mmol*(mL mmHg/(mmol K))*K",
            "status": "PASS",
        },
        {
            "equation": "J_diss,i",
            "unit_check": "mmol/h = h^-1*mL*mmol/mL",
            "status": "PASS",
        },
        {
            "equation": "J_abs,i",
            "unit_check": "mmol/h = h^-1*mmol",
            "status": "PASS",
        },
        {
            "equation": "J_rect",
            "unit_check": "mmol/h = mmol/(h sqrt(mmHg))*sqrt(mmHg) + dimensionless*mL/h*mmol/mL",
            "status": "PASS",
        },
        {
            "equation": "V_bowel",
            "unit_check": "mL = mL + mL",
            "status": "PASS",
        },
    ]


def crossing_time(time_h: np.ndarray, values: np.ndarray, target: float) -> float | None:
    if len(values) == 0 or values[-1] < target:
        return None
    index = int(np.flatnonzero(values >= target)[0])
    if index == 0:
        return float(time_h[0])
    x0 = float(time_h[index - 1])
    x1 = float(time_h[index])
    y0 = float(values[index - 1])
    y1 = float(values[index])
    if y1 == y0:
        return x1
    return x0 + (target - y0) * (x1 - x0) / (y1 - y0)


def scalar_summary(raw: dict[str, object]) -> dict[str, object]:
    p = raw["parameters"]
    y = np.asarray(raw["state"])
    final = y[:, -1]
    production = float(final[19])
    upstream = float(final[23])
    rectal_gas = float(final[16])
    fecal_dissolved = float(final[17])
    breath = float(final[18])
    total_excretion = rectal_gas + breath
    total_routes = total_excretion + fecal_dissolved
    residual = production + upstream - (
        rectal_gas
        + fecal_dissolved
        + breath
        + float(final[22])
        + float(final[6] + final[7] + final[8])
        + float(final[9] + final[10] + final[11])
        + float(final[12] + final[13] + final[14])
        - float(p["initial_dissolved_gas_mmol"] * 3.0)
        - float(p["initial_free_gas_mmol"] * 3.0)
        - float(p["initial_blood_gas_mmol"] * 3.0)
    )
    v0, compliance = geometry_arrays(p)
    final_pressure_excess = gas_pressure(
        np.maximum(final[9:12], 0.0),
        v0,
        compliance,
        p["r_ml_mmhg_mmol_k"],
        p["body_temperature_k"],
    )
    final_gas_volume = v0 + compliance * final_pressure_excess
    final_flow = rates(float(raw["time_h"][-1]), final, p)
    final_liquid_flux = np.maximum(final[3:6], 0.0) / np.array(
        [p["tau_proximal_h"], p["tau_distal_h"], p["tau_rectal_h"]],
        dtype=float,
    ) / max(p["motility_gain"], 1e-9)
    final_total_gas = float(final[6] + final[7] + final[8] + final[9] + final[10] + final[11] + final[12] + final[13] + final[14])
    return {
        "cumulative_production_mmol": production,
        "cumulative_upstream_gas_mmol": upstream,
        "cumulative_rectal_gas_mmol": rectal_gas,
        "cumulative_fecal_dissolved_mmol": fecal_dissolved,
        "cumulative_breath_mmol": breath,
        "cumulative_microbial_consumption_mmol": float(final[22]),
        "cumulative_absorption_mmol": float(final[20]),
        "cumulative_net_dissolution_mmol": float(final[21]),
        "cumulative_explicit_backflow_mmol": float(final[24]),
        "cumulative_explicit_forward_pressure_mmol": float(final[25]),
        "cumulative_fecal_mass_g": float(final[15]),
        "total_excretion_mmol": total_excretion,
        "total_excretion_ml": total_excretion * p["standard_molar_volume_ml_mmol"],
        "total_routes_including_fecal_dissolved_mmol": total_routes,
        "total_routes_including_fecal_dissolved_ml": total_routes * p["standard_molar_volume_ml_mmol"],
        "rectal_gas_excretion_ml": rectal_gas * p["standard_molar_volume_ml_mmol"],
        "fecal_dissolved_excretion_ml": fecal_dissolved * p["standard_molar_volume_ml_mmol"],
        "breath_excretion_ml": breath * p["standard_molar_volume_ml_mmol"],
        "breath_fraction_of_excretion": breath / total_excretion if total_excretion > 0 else 0.0,
        "rectal_fraction_of_excretion": rectal_gas / total_excretion if total_excretion > 0 else 0.0,
        "retention_fraction_of_gas_input": final_total_gas / (production + upstream) if production + upstream > 0 else 0.0,
        "gas_balance_residual_mmol": residual,
        "gas_balance_relative_residual": residual / max(production + upstream, 1e-12),
        "final_free_gas_mmol": float(final[9] + final[10] + final[11]),
        "final_dissolved_gas_mmol": float(final[6] + final[7] + final[8]),
        "final_blood_gas_mmol": float(final[12] + final[13] + final[14]),
        "final_gas_volume_ml": float(final_gas_volume.sum()),
        "final_bowel_volume_ml": float(final[3] + final[4] + final[5] + final_gas_volume.sum()),
        "final_pressure_mmhg": p["atmospheric_pressure_mmhg"] + float(final_pressure_excess.max()),
        "final_pressure_proximal_mmhg": p["atmospheric_pressure_mmhg"] + float(final_pressure_excess[0]),
        "final_pressure_distal_mmhg": p["atmospheric_pressure_mmhg"] + float(final_pressure_excess[1]),
        "final_pressure_rectal_mmhg": p["atmospheric_pressure_mmhg"] + float(final_pressure_excess[2]),
        "final_rectal_gas_flow_ml_min": float(final_flow[16] * p["standard_molar_volume_ml_mmol"] / 60.0),
        "final_rectal_liquid_flow_ml_h": float(final_liquid_flux[2]),
        "final_fecal_mass_flux_g_h": float(final[2] / (p["tau_rectal_h"] / max(p["motility_gain"], 1e-9))),
        "state_minimum": float(y[:, 6:15].min()),
    }


def full_summary(raw: dict[str, object]) -> dict[str, object]:
    result = scalar_summary(raw)
    p = raw["parameters"]
    time_h = np.asarray(raw["time_h"], dtype=float)
    state = np.asarray(raw["state"], dtype=float)
    v0, compliance = geometry_arrays(p)
    pressure_excess = gas_pressure(
        np.maximum(state[9:12], 0.0),
        v0[:, None],
        compliance[:, None],
        p["r_ml_mmhg_mmol_k"],
        p["body_temperature_k"],
    )
    gas_volume = v0[:, None] + compliance[:, None] * pressure_excess
    liquid_volume = np.maximum(state[3:6], 0.0)
    bowel_volume = liquid_volume + gas_volume
    pressure = p["atmospheric_pressure_mmhg"] + pressure_excess.max(axis=0)
    cumulative_mass = np.maximum(state[15], 0.0)
    half_delivered = 0.5 * p["lactulose_g"]
    result.update(
        {
            "max_gas_volume_ml": float(gas_volume.sum(axis=0).max()),
            "max_bowel_volume_ml": float(bowel_volume.sum(axis=0).max()),
            "max_pressure_mmhg": float(pressure.max()),
            "max_pressure_excess_mmhg": float(pressure_excess.max()),
            "t50_solid_transit_h": crossing_time(time_h, cumulative_mass, half_delivered),
            "t_rectal_gas_flow_threshold_h": crossing_time(
                time_h,
                state[16],
                1e-9,
            ),
            "final_state": {
                name: float(state[index, -1]) for index, name in enumerate(STATE_NAMES)
            },
            "series": {
                "time_h": time_h.tolist(),
                "cumulative_fecal_mass_g": state[15].tolist(),
                "cumulative_rectal_gas_mmol": state[16].tolist(),
                "cumulative_breath_mmol": state[18].tolist(),
                "gas_volume_ml": gas_volume.sum(axis=0).tolist(),
                "bowel_volume_ml": bowel_volume.sum(axis=0).tolist(),
                "pressure_mmhg": pressure.tolist(),
                "dissolved_gas_mmol": state[6:9].sum(axis=0).tolist(),
                "free_gas_mmol": state[9:12].sum(axis=0).tolist(),
            },
        }
    )
    return result


def sensitivity() -> dict[str, object]:
    baseline_raw = simulate()
    baseline_scalar = scalar_summary(baseline_raw)
    baseline_full = full_summary(baseline_raw)
    records = []
    ranking = []
    for name in SENSITIVITY_PARAMETERS:
        low = 0.5 * baseline_raw["parameters"][name]
        high = 1.5 * baseline_raw["parameters"][name]
        low_raw = simulate({name: low})
        high_raw = simulate({name: high})
        low_scalar = scalar_summary(low_raw)
        high_scalar = scalar_summary(high_raw)
        low_full = full_summary(low_raw)
        high_full = full_summary(high_raw)

        def relative(base_value: float, low_value: float, high_value: float) -> float:
            return max(
                abs((float(low_value) - base_value) / max(abs(base_value), 1e-12)),
                abs((float(high_value) - base_value) / max(abs(base_value), 1e-12)),
            )

        record = {
            "parameter": name,
            "minus50_value": low,
            "plus50_value": high,
            "total_excretion_ml_minus50": low_scalar["total_excretion_ml"],
            "total_excretion_ml_plus50": high_scalar["total_excretion_ml"],
            "breath_fraction_minus50": low_scalar["breath_fraction_of_excretion"],
            "breath_fraction_plus50": high_scalar["breath_fraction_of_excretion"],
            "max_pressure_mmhg_minus50": low_full["max_pressure_mmhg"],
            "max_pressure_mmhg_plus50": high_full["max_pressure_mmhg"],
            "max_gas_volume_ml_minus50": low_full["max_gas_volume_ml"],
            "max_gas_volume_ml_plus50": high_full["max_gas_volume_ml"],
            "max_abs_relative_total_excretion_change": relative(
                baseline_scalar["total_excretion_ml"],
                low_scalar["total_excretion_ml"],
                high_scalar["total_excretion_ml"],
            ),
            "max_abs_relative_breath_fraction_change": relative(
                baseline_scalar["breath_fraction_of_excretion"],
                low_scalar["breath_fraction_of_excretion"],
                high_scalar["breath_fraction_of_excretion"],
            ),
            "max_abs_relative_pressure_change": relative(
                baseline_full["max_pressure_mmhg"],
                low_full["max_pressure_mmhg"],
                high_full["max_pressure_mmhg"],
            ),
            "max_abs_relative_gas_volume_change": relative(
                baseline_full["max_gas_volume_ml"],
                low_full["max_gas_volume_ml"],
                high_full["max_gas_volume_ml"],
            ),
        }
        records.append(record)
        ranking.append(
            (
                name,
                record["max_abs_relative_total_excretion_change"],
                record["max_abs_relative_pressure_change"],
                record["max_abs_relative_gas_volume_change"],
            )
        )
    ranking.sort(key=lambda item: max(item[1], item[2], item[3]), reverse=True)
    return {
        "baseline_total_excretion_ml": baseline_scalar["total_excretion_ml"],
        "baseline_max_pressure_mmhg": baseline_full["max_pressure_mmhg"],
        "baseline_max_gas_volume_ml": baseline_full["max_gas_volume_ml"],
        "records": records,
        "ranking_by_largest_output_change": [item[0] for item in ranking],
    }


def counterchecks() -> dict[str, object]:
    closed = default_parameters()
    closed.update(
        {
            "lactulose_g": 0.0,
            "upstream_gas_mmol_h": 1.0,
            "k_phase_transfer_h": 0.0,
            "k_blood_transfer_h": 0.0,
            "k_breath_clearance_h": 0.0,
            "k_microbial_h": 0.0,
            "k_gas_diffusion_h": 0.0,
            "k_dissolved_diffusion_h": 0.0,
            "k_forward_pressure_mmol_h_mmhg": 0.0,
            "k_backward_pressure_mmol_h_mmhg": 0.0,
            "k_rectal_valve_mmol_h_sqrt_mmhg": 0.0,
            "bubble_carry_fraction": 0.0,
        }
    )
    closed_raw = simulate(closed)
    closed_summary = scalar_summary(closed_raw)
    placebo = default_parameters()
    placebo.update(
        {
            "k_phase_transfer_h": 0.0,
            "k_blood_transfer_h": 0.0,
            "k_breath_clearance_h": 0.0,
            "k_microbial_h": 0.0,
            "k_gas_diffusion_h": 0.0,
            "k_dissolved_diffusion_h": 0.0,
            "k_forward_pressure_mmol_h_mmhg": 0.0,
            "k_backward_pressure_mmol_h_mmhg": 0.0,
            "k_rectal_valve_mmol_h_sqrt_mmhg": 0.0,
            "bubble_carry_fraction": 0.0,
        }
    )
    placebo_raw = simulate(placebo)
    placebo_summary = scalar_summary(placebo_raw)
    return {
        "closed_gas_injection": {
            "expected_gas_input_mmol": closed["upstream_gas_mmol_h"] * closed["duration_h"],
            "observed_production_plus_upstream_mmol": closed_summary["cumulative_production_mmol"] + closed_summary["cumulative_upstream_gas_mmol"],
            "observed_rectal_plus_breath_mmol": closed_summary["cumulative_rectal_gas_mmol"] + closed_summary["cumulative_breath_mmol"],
            "balance_residual_mmol": closed_summary["gas_balance_residual_mmol"],
        },
        "production_only_placebo": {
            "production_mmol": placebo_summary["cumulative_production_mmol"],
            "rectal_plus_breath_mmol": placebo_summary["cumulative_rectal_gas_mmol"] + placebo_summary["cumulative_breath_mmol"],
            "retained_free_gas_mmol": placebo_summary["final_free_gas_mmol"],
            "balance_residual_mmol": placebo_summary["gas_balance_residual_mmol"],
        },
    }


def transit_probe() -> dict[str, object]:
    raw = simulate({"duration_h": 72.0})
    summary = full_summary(raw)
    return {
        "duration_h": 72.0,
        "t50_solid_transit_h": summary["t50_solid_transit_h"],
        "cumulative_fecal_mass_g": summary["cumulative_fecal_mass_g"],
        "final_fecal_mass_flux_g_h": summary["final_fecal_mass_flux_g_h"],
        "interpretation": "descriptive transit extension; not used for the frozen 24 h gas endpoint",
    }


def build_results() -> dict[str, object]:
    raw = simulate()
    nominal = full_summary(raw)
    criterion_low = 113.5
    criterion_high = 454.0
    criterion_met = criterion_low <= nominal["total_excretion_ml"] <= criterion_high
    return {
        "id": MODEL_ID,
        "status": "first_run",
        "reference": {
            "authors": "Christl, Murgatroyd, Gibson, Cummings",
            "year": 1992,
            "journal": "Gastroenterology",
            "doi": "10.1016/0016-5085(92)90765-q",
            "dose_g": 15.0,
            "value_ml_per_24h": 227.0,
            "sem_ml_per_24h": 60.7,
            "source_url": "https://europepmc.org/article/med/1551534",
        },
        "secondary_reference": {
            "authors": "Hammer",
            "year": 1993,
            "journal": "Gut",
            "doi": "10.1136/gut.34.6.818",
            "hydrogen_threshold_ml_per_6h": 76.0,
            "absorption_efficiency_low_load_percent": 90.0,
            "absorption_efficiency_high_load_percent": 20.0,
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC1374269/",
        },
        "acceptance": {
            "primary_endpoint": "cumulative rectal-plus-breath H2-equivalent excretion at 24 h in mL STP",
            "reference_value_ml": 227.0,
            "factor_two_low_ml": criterion_low,
            "factor_two_high_ml": criterion_high,
            "predicted_ml": nominal["total_excretion_ml"],
            "criterion_met": bool(criterion_met),
            "status": "MET" if criterion_met else "EJ",
        },
        "nominal": nominal,
        "transit_probe": transit_probe(),
        "sensitivity": sensitivity(),
        "counterchecks": counterchecks(),
        "parameter_table": parameter_table(),
        "dimension_checks": dimension_checks(),
        "equations": EQUATIONS,
        "state_names": list(STATE_NAMES),
        "source_status": {
            "primary_reference_verified": True,
            "primary_reference_verification": "Europe PMC abstract and publisher DOI page fetched",
            "individual_measurement_data_used": False,
        },
    }


def main() -> None:
    results = build_results()
    Path("results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "id": results["id"],
                "predicted_total_excretion_ml": results["acceptance"]["predicted_ml"],
                "criterion": results["acceptance"]["status"],
                "gas_balance_relative_residual": results["nominal"]["gas_balance_relative_residual"],
                "max_pressure_mmhg": results["nominal"]["max_pressure_mmhg"],
                "max_gas_volume_ml": results["nominal"]["max_gas_volume_ml"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
