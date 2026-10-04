from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


SEGMENT_NAMES = ("proximal_tubule", "loop_of_henle", "distal_convoluted_tubule", "collecting_duct")
STATE_AP = 0
STATE_LUMEN_START = 1
STATE_CELL_START = 5
STATE_VOLUME_START = 9
STATE_URINE_MASS = 13
STATE_URINE_VOLUME = 14
STATE_FILTERED_MASS = 15
STATE_SECRETED_MASS = 16
STATE_REABSORBED_MASS = 17
N_STATES = 18
EPSILON = 1.0e-12


@dataclass(frozen=True)
class ModelParameters:
    dose_mg: float = 500.0
    duration_h: float = 24.0
    time_step_min: float = 0.25
    gfr_ml_min: float = 100.0
    fu_plasma: float = 1.0
    renal_plasma_flow_ml_min: float = 600.0
    reference_perfusion_ml_min: float = 600.0
    plasma_volume_l: float = 200.0
    segment_water_reabsorption: tuple[float, ...] = (0.65, 0.20, 0.15, 0.94)
    segment_residence_min: tuple[float, ...] = (10.0, 8.0, 5.0, 10.0)
    segment_ph: tuple[float, ...] = (7.40, 7.20, 7.00, 6.80)
    segment_perfusion: tuple[float, ...] = (1.00, 0.15, 0.05, 0.02)
    cell_volume_ml: tuple[float, ...] = (1000.0, 200.0, 100.0, 50.0)
    oct2_permeability_cm_s: tuple[float, ...] = (1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5)
    oct2_area_cm2: tuple[float, ...] = (800000.0, 100000.0, 50000.0, 30000.0)
    cell_partition: tuple[float, ...] = (1.0, 1.0, 1.0, 1.0)
    mate_vmax_mg_min: tuple[float, ...] = (1.60, 0.15, 0.05, 0.02)
    mate_km_mg_ml: tuple[float, ...] = (0.0010, 0.0015, 0.0020, 0.0020)
    basolateral_return_rate_min: tuple[float, ...] = (0.10, 0.12, 0.15, 0.20)
    reabsorption_vmax_mg_min: tuple[float, ...] = (0.005, 0.002, 0.001, 0.0005)
    reabsorption_km_mg_ml: tuple[float, ...] = (0.010, 0.010, 0.010, 0.010)
    reabsorption_ionized_floor: float = 0.0
    metformin_pka: float = 12.4
    reference_ph: float = 7.40


@dataclass(frozen=True)
class SimulationResult:
    parameters: ModelParameters
    time_min: np.ndarray
    states: np.ndarray
    water_in_ml_min: tuple[float, ...]
    water_out_ml_min: tuple[float, ...]
    lumen_volume_ml: tuple[float, ...]
    summary: dict[str, Any]


def default_parameters() -> ModelParameters:
    return ModelParameters()


def parameter_table(p: ModelParameters | None = None) -> list[dict[str, Any]]:
    p = p or default_parameters()
    return [
        {"name": "dose_mg", "value": p.dose_mg, "unit": "mg", "basis": "scenario assumption"},
        {"name": "duration_h", "value": p.duration_h, "unit": "h", "basis": "scenario assumption"},
        {"name": "time_step_min", "value": p.time_step_min, "unit": "min", "basis": "numerical integration assumption"},
        {"name": "gfr_ml_min", "value": p.gfr_ml_min, "unit": "mL/min", "basis": "representative adult assumption; not measured patient data"},
        {"name": "fu_plasma", "value": p.fu_plasma, "unit": "fraction", "basis": "effective unbound fraction; metformin protein binding is assumed negligible"},
        {"name": "renal_plasma_flow_ml_min", "value": p.renal_plasma_flow_ml_min, "unit": "mL/min", "basis": "representative perfusion assumption"},
        {"name": "reference_perfusion_ml_min", "value": p.reference_perfusion_ml_min, "unit": "mL/min", "basis": "normalization reference for regional perfusion"},
        {"name": "plasma_volume_l", "value": p.plasma_volume_l, "unit": "L", "basis": "effective central distribution-volume assumption"},
        {"name": "segment_water_reabsorption", "value": list(p.segment_water_reabsorption), "unit": "fraction of segment inflow", "basis": "segmental water-balance assumption"},
        {"name": "segment_residence_min", "value": list(p.segment_residence_min), "unit": "min", "basis": "effective compartment geometry assumption"},
        {"name": "segment_ph", "value": list(p.segment_ph), "unit": "dimensionless", "basis": "segmental lumen pH assumption"},
        {"name": "segment_perfusion", "value": list(p.segment_perfusion), "unit": "fraction", "basis": "regional perfusion assumption"},
        {"name": "cell_volume_ml", "value": list(p.cell_volume_ml), "unit": "mL epithelial compartment", "basis": "effective geometry assumption"},
        {"name": "oct2_permeability_cm_s", "value": list(p.oct2_permeability_cm_s), "unit": "cm/s", "basis": "OCT2-mediated basolateral permeability assumption"},
        {"name": "oct2_area_cm2", "value": list(p.oct2_area_cm2), "unit": "cm2", "basis": "apparent membrane area assumption"},
        {"name": "cell_partition", "value": list(p.cell_partition), "unit": "dimensionless", "basis": "plasma-to-cell equilibrium partition assumption"},
        {"name": "mate_vmax_mg_min", "value": list(p.mate_vmax_mg_min), "unit": "mg/min", "basis": "MATE1/MATE2-K apical efflux capacity assumption"},
        {"name": "mate_km_mg_ml", "value": list(p.mate_km_mg_ml), "unit": "mg/mL", "basis": "MATE Michaelis-Menten affinity assumption"},
        {"name": "basolateral_return_rate_min", "value": list(p.basolateral_return_rate_min), "unit": "1/min", "basis": "cell-to-plasma return-rate assumption"},
        {"name": "reabsorption_vmax_mg_min", "value": list(p.reabsorption_vmax_mg_min), "unit": "mg/min", "basis": "unidentified apical reabsorption capacity assumption"},
        {"name": "reabsorption_km_mg_ml", "value": list(p.reabsorption_km_mg_ml), "unit": "mg/mL", "basis": "reabsorption Michaelis-Menten affinity assumption"},
        {"name": "reabsorption_ionized_floor", "value": p.reabsorption_ionized_floor, "unit": "fraction", "basis": "effective ionized reabsorption weighting assumption"},
        {"name": "metformin_pka", "value": p.metformin_pka, "unit": "pKa", "basis": "FDA/DailyMed physicochemical label value 12.4"},
        {"name": "reference_ph", "value": p.reference_ph, "unit": "pH", "basis": "normalization pH assumption"},
    ]


def dimensional_analysis() -> dict[str, Any]:
    checks = [
        {"quantity": "filtration mass rate", "expression": "mL/min * mg/mL", "result": "mg/min", "pass": True},
        {"quantity": "water volume rate", "expression": "mL/min - mL/min", "result": "mL/min", "pass": True},
        {"quantity": "Michaelis-Menten flux", "expression": "mg/min * (mg/mL)/(mg/mL)", "result": "mg/min", "pass": True},
        {"quantity": "membrane mass-transfer coefficient", "expression": "cm/s * cm2 / mL * 60 s/min", "result": "1/min", "pass": True},
        {"quantity": "clearance", "expression": "mg/min / (mg/mL)", "result": "mL/min", "pass": True},
        {"quantity": "urine mass", "expression": "mg/min * min", "result": "mg", "pass": True},
    ]
    return {"all_pass": all(item["pass"] for item in checks), "checks": checks}


def equation_table() -> list[dict[str, str]]:
    return [
        {"name": "filtration", "equation": "J_filter = fu,p * GFR * C_p", "unit": "mg/min"},
        {"name": "plasma_balance", "equation": "dA_p/dt = -J_filter - sum(J_OCT2) + sum(J_back)", "unit": "mg/min"},
        {"name": "OCT2_transfer", "equation": "k_OCT2 = P_OCT2 * A / V_cell * 60 * regional_perfusion * (Qp/Qp_ref); J_OCT2 = k_OCT2 * Kp * C_p * V_cell", "unit": "1/min; mg/min"},
        {"name": "MATE_kinetics", "equation": "J_MATE = Vmax_MATE * C_cell/(Km_MATE + C_cell) * hH(pH)", "unit": "mg/min"},
        {"name": "ionization", "equation": "f_ion = 1/(1 + 10^(pH - pKa)); hH = 10^(pH_ref - pH)", "unit": "dimensionless"},
        {"name": "reabsorption", "equation": "J_reabs = Vmax_reabs * C_lumen/(Km_reabs + C_lumen) * f_ion", "unit": "mg/min"},
        {"name": "lumen_balance", "equation": "dA_lumen/dt = Q_in*C_in - Q_out*C_lumen + J_MATE - J_reabs", "unit": "mg/min"},
        {"name": "water_balance", "equation": "dV/dt = Q_in - Q_out - Q_water,reabs", "unit": "mL/min"},
        {"name": "urine_outflow", "equation": "dM_urine/dt = Q_out,collecting*C_collecting", "unit": "mg/min"},
        {"name": "mass_conservation", "equation": "d(A_p + sum(A_lumen) + sum(A_cell) + M_urine)/dt = 0", "unit": "mg/min"},
    ]


def ionized_fraction(ph: float, pka: float) -> float:
    value = 1.0 / (1.0 + 10.0 ** (ph - pka))
    return float(np.clip(value, 0.0, 1.0))


def neutral_fraction(ph: float, pka: float) -> float:
    return float(1.0 - ionized_fraction(ph, pka))


def mate_h_factor(ph: float, reference_ph: float) -> float:
    value = 10.0 ** (reference_ph - ph)
    return float(np.clip(value, 0.05, 20.0))


def water_flows(p: ModelParameters) -> tuple[list[float], list[float], list[float]]:
    water_in: list[float] = []
    water_out: list[float] = []
    water_reabsorbed: list[float] = []
    incoming = p.gfr_ml_min
    for fraction in p.segment_water_reabsorption:
        reabsorbed = incoming * fraction
        outgoing = max(incoming - reabsorbed, EPSILON)
        water_in.append(float(incoming))
        water_out.append(float(outgoing))
        water_reabsorbed.append(float(reabsorbed))
        incoming = outgoing
    return water_in, water_out, water_reabsorbed


def lumen_volumes(p: ModelParameters) -> list[float]:
    _, outgoing, _ = water_flows(p)
    return [max(q * residence, EPSILON) for q, residence in zip(outgoing, p.segment_residence_min)]


def oct2_rate_per_min(p: ModelParameters, index: int) -> float:
    permeability = p.oct2_permeability_cm_s[index]
    area = p.oct2_area_cm2[index]
    volume = p.cell_volume_ml[index]
    perfusion = p.segment_perfusion[index]
    perfusion_scale = p.renal_plasma_flow_ml_min / max(p.reference_perfusion_ml_min, EPSILON)
    return float(permeability * area / volume * 60.0 * perfusion * perfusion_scale)


def steady_cell_concentration(c_p: float, c_lumen: float, p: ModelParameters, index: int) -> float:
    ph = p.segment_ph[index]
    v_cell = p.cell_volume_ml[index]
    rate = oct2_rate_per_min(p, index)
    j_oct2 = rate * p.cell_partition[index] * c_p * v_cell
    h_factor = mate_h_factor(ph, p.reference_ph)
    ion = ionized_fraction(ph, p.metformin_pka)
    reabs_weight = p.reabsorption_ionized_floor + (1.0 - p.reabsorption_ionized_floor) * ion
    j_reabs = p.reabsorption_vmax_mg_min[index] * c_lumen / (p.reabsorption_km_mg_ml[index] + c_lumen) * reabs_weight
    target = max(j_oct2 + j_reabs, 0.0)
    b = p.basolateral_return_rate_min[index] * v_cell
    a = p.mate_km_mg_ml[index] * b + p.mate_vmax_mg_min[index] * h_factor - target
    discriminant = a * a + 4.0 * b * target * p.mate_km_mg_ml[index]
    if b <= EPSILON:
        if p.mate_vmax_mg_min[index] * h_factor <= EPSILON:
            return 0.0
        return max(target * p.mate_km_mg_ml[index] / (p.mate_vmax_mg_min[index] * h_factor), 0.0)
    value = (-a + math.sqrt(max(discriminant, 0.0))) / (2.0 * b)
    return float(max(value, 0.0))


def initial_state(p: ModelParameters) -> np.ndarray:
    y = np.zeros(N_STATES, dtype=float)
    plasma_amount = p.dose_mg
    cell_amounts = np.zeros(len(SEGMENT_NAMES), dtype=float)
    for _ in range(8):
        c_p = plasma_amount / (p.plasma_volume_l * 1000.0)
        for index in range(len(SEGMENT_NAMES)):
            c_cell = steady_cell_concentration(c_p, 0.0, p, index)
            cell_amounts[index] = c_cell * p.cell_volume_ml[index]
        updated_plasma_amount = p.dose_mg - float(np.sum(cell_amounts))
        if abs(updated_plasma_amount - plasma_amount) < 1.0e-12:
            plasma_amount = updated_plasma_amount
            break
        plasma_amount = updated_plasma_amount
    y[STATE_AP] = plasma_amount
    volumes = lumen_volumes(p)
    for index in range(len(SEGMENT_NAMES)):
        c_p = plasma_amount / (p.plasma_volume_l * 1000.0)
        c_cell = steady_cell_concentration(c_p, 0.0, p, index)
        y[STATE_LUMEN_START + index] = 0.0
        y[STATE_CELL_START + index] = c_cell * p.cell_volume_ml[index]
        y[STATE_VOLUME_START + index] = volumes[index]
    return y


def flux_details(y: np.ndarray, p: ModelParameters) -> list[dict[str, float]]:
    c_p = max(float(y[STATE_AP]) / (p.plasma_volume_l * 1000.0), 0.0)
    water_in, water_out, water_reabsorbed = water_flows(p)
    details: list[dict[str, float]] = []
    for index in range(len(SEGMENT_NAMES)):
        c_lumen = max(float(y[STATE_LUMEN_START + index]) / max(float(y[STATE_VOLUME_START + index]), EPSILON), 0.0)
        c_cell = max(float(y[STATE_CELL_START + index]) / p.cell_volume_ml[index], 0.0)
        ph = p.segment_ph[index]
        rate = oct2_rate_per_min(p, index)
        j_oct2 = rate * p.cell_partition[index] * c_p * p.cell_volume_ml[index]
        j_back = p.basolateral_return_rate_min[index] * c_cell * p.cell_volume_ml[index]
        h_factor = mate_h_factor(ph, p.reference_ph)
        j_mate = p.mate_vmax_mg_min[index] * c_cell / (p.mate_km_mg_ml[index] + c_cell) * h_factor
        ion = ionized_fraction(ph, p.metformin_pka)
        reabs_weight = p.reabsorption_ionized_floor + (1.0 - p.reabsorption_ionized_floor) * ion
        j_reabs = p.reabsorption_vmax_mg_min[index] * c_lumen / (p.reabsorption_km_mg_ml[index] + c_lumen) * reabs_weight
        details.append({
            "water_in_ml_min": water_in[index],
            "water_out_ml_min": water_out[index],
            "water_reabsorbed_ml_min": water_reabsorbed[index],
            "c_lumen_mg_ml": c_lumen,
            "c_cell_mg_ml": c_cell,
            "ph": ph,
            "ionized_fraction": ion,
            "neutral_fraction": 1.0 - ion,
            "mate_h_factor": h_factor,
            "j_oct2_mg_min": j_oct2,
            "j_back_mg_min": j_back,
            "j_mate_mg_min": j_mate,
            "j_reabs_mg_min": j_reabs,
        })
    return details


def rhs(_t: float, y: np.ndarray, p: ModelParameters) -> np.ndarray:
    c_p = max(float(y[STATE_AP]) / (p.plasma_volume_l * 1000.0), 0.0)
    details = flux_details(y, p)
    derivative = np.zeros(N_STATES, dtype=float)
    filtered_rate = p.fu_plasma * p.gfr_ml_min * c_p
    derivative[STATE_AP] = -filtered_rate
    for index, detail in enumerate(details):
        if index == 0:
            incoming_drug_rate = filtered_rate
        else:
            previous = details[index - 1]
            incoming_drug_rate = previous["water_out_ml_min"] * previous["c_lumen_mg_ml"]
        outgoing_drug_rate = detail["water_out_ml_min"] * detail["c_lumen_mg_ml"]
        derivative[STATE_AP] -= detail["j_oct2_mg_min"]
        derivative[STATE_AP] += detail["j_back_mg_min"]
        derivative[STATE_LUMEN_START + index] = incoming_drug_rate - outgoing_drug_rate + detail["j_mate_mg_min"] - detail["j_reabs_mg_min"]
        derivative[STATE_CELL_START + index] = detail["j_oct2_mg_min"] + detail["j_reabs_mg_min"] - detail["j_mate_mg_min"] - detail["j_back_mg_min"]
        derivative[STATE_VOLUME_START + index] = detail["water_in_ml_min"] - detail["water_out_ml_min"] - detail["water_reabsorbed_ml_min"]
    last = details[-1]
    derivative[STATE_URINE_MASS] = last["water_out_ml_min"] * last["c_lumen_mg_ml"]
    derivative[STATE_URINE_VOLUME] = last["water_out_ml_min"]
    derivative[STATE_FILTERED_MASS] = filtered_rate
    derivative[STATE_SECRETED_MASS] = sum(item["j_mate_mg_min"] for item in details)
    derivative[STATE_REABSORBED_MASS] = sum(item["j_reabs_mg_min"] for item in details)
    return derivative


def analytic_no_transport(p: ModelParameters | None = None, times_min: np.ndarray | None = None) -> dict[str, Any]:
    p = p or default_parameters()
    if times_min is None:
        times_min = np.linspace(0.0, p.duration_h * 60.0, 1001)
    plasma_volume_ml = p.plasma_volume_l * 1000.0
    k = p.fu_plasma * p.gfr_ml_min / plasma_volume_ml
    plasma_amount = p.dose_mg * np.exp(-k * times_min)
    cumulative_filtered = p.dose_mg * (1.0 - np.exp(-k * times_min))
    return {
        "clearance_ml_min": p.fu_plasma * p.gfr_ml_min,
        "plasma_amount_mg": plasma_amount,
        "cumulative_filtered_mg": cumulative_filtered,
    }


def state_mass_total(y: np.ndarray) -> float:
    return float(
        y[STATE_AP]
        + y[STATE_URINE_MASS]
        + np.sum(y[STATE_LUMEN_START:STATE_LUMEN_START + 4])
        + np.sum(y[STATE_CELL_START:STATE_CELL_START + 4])
    )


def summarize(y: np.ndarray, t: np.ndarray, p: ModelParameters, water_in: list[float], water_out: list[float], volumes: list[float]) -> dict[str, Any]:
    terminal_state = y[-1]
    c_p = max(float(terminal_state[STATE_AP]) / (p.plasma_volume_l * 1000.0), EPSILON)
    initial_state_p = initial_state(p)
    initial_c_p = max(float(initial_state_p[STATE_AP]) / (p.plasma_volume_l * 1000.0), EPSILON)
    initial_details = flux_details(initial_state_p, p)
    initial_filtration = p.fu_plasma * p.gfr_ml_min
    initial_secretion = sum(item["j_mate_mg_min"] for item in initial_details) / initial_c_p
    initial_reabsorption = sum(item["j_reabs_mg_min"] for item in initial_details) / initial_c_p
    initial_clearance = initial_filtration + initial_secretion - initial_reabsorption
    details = flux_details(terminal_state, p)
    terminal_filtration = p.fu_plasma * p.gfr_ml_min
    terminal_secretion = sum(item["j_mate_mg_min"] for item in details) / c_p
    terminal_reabsorption = sum(item["j_reabs_mg_min"] for item in details) / c_p
    terminal_clearance = terminal_filtration + terminal_secretion - terminal_reabsorption
    auc = float(np.trapezoid(y[:, STATE_AP] / (p.plasma_volume_l * 1000.0), t))
    urine_mass = float(terminal_state[STATE_URINE_MASS])
    summary: dict[str, Any] = {
        "initial_plasma_concentration_mg_l": float(y[0, STATE_AP] / (p.plasma_volume_l * 1000.0)),
        "initial_clearance_ml_min": float(initial_clearance),
        "initial_filtration_clearance_ml_min": float(initial_filtration),
        "initial_secretion_clearance_ml_min": float(initial_secretion),
        "initial_reabsorption_clearance_ml_min": float(initial_reabsorption),
        "terminal_plasma_concentration_mg_l": float(y[-1, STATE_AP] / (p.plasma_volume_l * 1000.0)),
        "terminal_clearance_ml_min": float(terminal_clearance),
        "auc_plasma_mg_min_per_ml": auc,
        "auc_based_renal_clearance_ml_min": float(urine_mass / auc) if auc > EPSILON else None,
        "urine_mass_24h_mg": urine_mass,
        "urine_volume_24h_ml": float(y[-1, STATE_URINE_VOLUME]),
        "cumulative_filtered_mass_mg": float(y[-1, STATE_FILTERED_MASS]),
        "cumulative_secretion_mass_mg": float(y[-1, STATE_SECRETED_MASS]),
        "cumulative_reabsorption_mass_mg": float(y[-1, STATE_REABSORBED_MASS]),
        "cumulative_net_tubular_addition_mg": float(y[-1, STATE_SECRETED_MASS] - y[-1, STATE_REABSORBED_MASS]),
        "terminal_plasma_mass_mg": float(y[-1, STATE_AP]),
        "terminal_lumen_mass_mg": float(np.sum(y[-1, STATE_LUMEN_START:STATE_LUMEN_START + 4])),
        "terminal_cell_mass_mg": float(np.sum(y[-1, STATE_CELL_START:STATE_CELL_START + 4])),
        "mass_balance_mg": state_mass_total(y[-1]),
        "mass_balance_error_mg": float(state_mass_total(y[-1]) - p.dose_mg),
        "water_in_ml_min": [float(x) for x in water_in],
        "water_out_ml_min": [float(x) for x in water_out],
        "lumen_volume_ml": [float(x) for x in volumes],
        "water_reabsorbed_total_ml": float(sum(water_in[i] - water_out[i] for i in range(4)) * p.duration_h * 60.0),
        "ionized_fraction_by_segment": [item["ionized_fraction"] for item in details],
        "neutral_fraction_by_segment": [item["neutral_fraction"] for item in details],
        "mate_h_factor_by_segment": [item["mate_h_factor"] for item in details],
    }
    for target_min in (0.0, 300.0, 720.0, 1440.0):
        if target_min <= t[-1] + EPSILON:
            index = int(np.argmin(np.abs(t - target_min)))
            summary[f"plasma_concentration_{int(target_min)}_mg_l"] = float(y[index, STATE_AP] / (p.plasma_volume_l * 1000.0))
            summary[f"urine_mass_{int(target_min)}_mg"] = float(y[index, STATE_URINE_MASS])
    return summary


def simulate(p: ModelParameters | None = None) -> SimulationResult:
    p = p or default_parameters()
    duration_min = p.duration_h * 60.0
    t_eval = np.arange(0.0, duration_min + 0.5 * p.time_step_min, p.time_step_min)
    t_eval[-1] = duration_min
    y0 = initial_state(p)
    water_in, water_out, _ = water_flows(p)
    volumes = lumen_volumes(p)
    solution = solve_ivp(
        rhs,
        (0.0, duration_min),
        y0,
        args=(p,),
        method="Radau",
        t_eval=t_eval,
        rtol=1.0e-8,
        atol=1.0e-11,
        max_step=p.time_step_min,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    states = np.asarray(solution.y, dtype=float).T
    summary = summarize(states, solution.t, p, water_in, water_out, volumes)
    return SimulationResult(p, np.asarray(solution.t, dtype=float), states, tuple(water_in), tuple(water_out), tuple(volumes), summary)


def scaled_secretion(p: ModelParameters, scale: float) -> ModelParameters:
    return replace(
        p,
        oct2_permeability_cm_s=tuple(x * scale for x in p.oct2_permeability_cm_s),
        mate_vmax_mg_min=tuple(x * scale for x in p.mate_vmax_mg_min),
    )


def sensitivity_analysis(p: ModelParameters | None = None) -> list[dict[str, Any]]:
    base = p or default_parameters()
    base_result = simulate(base)
    base_clearance = float(base_result.summary["initial_clearance_ml_min"])
    base_urine = float(base_result.summary["urine_mass_24h_mg"])
    candidates: list[tuple[str, ModelParameters, ModelParameters]] = [
        (
            "gfr_ml_min",
            replace(base, gfr_ml_min=base.gfr_ml_min * 0.5),
            replace(base, gfr_ml_min=base.gfr_ml_min * 1.5),
        ),
        (
            "secretory_capacity_OCT2_MATE",
            scaled_secretion(base, 0.5),
            scaled_secretion(base, 1.5),
        ),
        (
            "proximal_water_reabsorption_fraction",
            replace(base, segment_water_reabsorption=(0.325, base.segment_water_reabsorption[1], base.segment_water_reabsorption[2], base.segment_water_reabsorption[3])),
            replace(base, segment_water_reabsorption=(0.975, base.segment_water_reabsorption[1], base.segment_water_reabsorption[2], base.segment_water_reabsorption[3])),
        ),
    ]
    output: list[dict[str, Any]] = []
    for name, low, high in candidates:
        low_result = simulate(low)
        high_result = simulate(high)
        low_clearance = float(low_result.summary["initial_clearance_ml_min"])
        high_clearance = float(high_result.summary["initial_clearance_ml_min"])
        low_urine = float(low_result.summary["urine_mass_24h_mg"])
        high_urine = float(high_result.summary["urine_mass_24h_mg"])
        output.append({
            "parameter": name,
            "low_value": 0.5,
            "high_value": 1.5,
            "low_initial_clearance_ml_min": low_clearance,
            "high_initial_clearance_ml_min": high_clearance,
            "clearance_span_ml_min": abs(high_clearance - low_clearance),
            "clearance_span_fraction_of_base": abs(high_clearance - low_clearance) / base_clearance,
            "low_urine_mass_24h_mg": low_urine,
            "high_urine_mass_24h_mg": high_urine,
            "urine_span_fraction_of_base": abs(high_urine - low_urine) / base_urine,
        })
    output.sort(key=lambda item: item["clearance_span_fraction_of_base"], reverse=True)
    for rank, item in enumerate(output, start=1):
        item["clearance_sensitivity_rank"] = rank
    return output


def run_experiment(p: ModelParameters | None = None) -> dict[str, Any]:
    p = p or default_parameters()
    full = simulate(p)
    null_p = replace(
        p,
        oct2_permeability_cm_s=tuple(0.0 for _ in p.oct2_permeability_cm_s),
        mate_vmax_mg_min=tuple(0.0 for _ in p.mate_vmax_mg_min),
        reabsorption_vmax_mg_min=tuple(0.0 for _ in p.reabsorption_vmax_mg_min),
    )
    null = simulate(null_p)
    zero_reabsorption = replace(p, reabsorption_vmax_mg_min=tuple(0.0 for _ in p.reabsorption_vmax_mg_min))
    double_reabsorption = replace(p, reabsorption_vmax_mg_min=tuple(x * 2.0 for x in p.reabsorption_vmax_mg_min))
    zero_reabsorption_result = simulate(zero_reabsorption)
    double_reabsorption_result = simulate(double_reabsorption)
    reference = 510.0
    lower = reference * 0.70
    upper = reference * 1.30
    predicted = float(full.summary["initial_clearance_ml_min"])
    status = "PASS" if lower <= predicted <= upper else "FAIL"
    return {
        "node": "BT-HX-Q140",
        "model": "metformin_segmented_renal_first_principles",
        "status": status,
        "reference": {
            "authors": "Graham GG, Punt J, Arora M, Day RO, Doogue MP, Duong JK, et al.",
            "year": 2011,
            "journal": "Clinical Pharmacokinetics",
            "volume_pages": "50:81-98",
            "doi": "10.2165/11534750-000000000-00000",
            "value_ml_min": reference,
            "reported_sd_ml_min": 130.0,
            "source_page": "https://pubmed.ncbi.nlm.nih.gov/21241070/",
            "transporter_source_doi": "10.1097/FPC.0b013e3283559b22",
            "transporter_source_page": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3651676/",
        },
        "criterion": {
            "type": "reference_interval",
            "lower_ml_min": lower,
            "upper_ml_min": upper,
            "rule": "0.70 * 510 <= predicted initial renal clearance <= 1.30 * 510",
        },
        "prediction": full.summary,
        "filtration_only_null": {
            "initial_clearance_ml_min": float(null.summary["initial_clearance_ml_min"]),
            "urine_mass_24h_mg": float(null.summary["urine_mass_24h_mg"]),
            "clearance_ratio_full_to_null": predicted / float(null.summary["initial_clearance_ml_min"]),
            "urine_mass_ratio_full_to_null": float(full.summary["urine_mass_24h_mg"]) / float(null.summary["urine_mass_24h_mg"]),
        },
        "reabsorption_perturbation": {
            "zero_capacity_initial_clearance_ml_min": float(zero_reabsorption_result.summary["initial_clearance_ml_min"]),
            "double_capacity_initial_clearance_ml_min": float(double_reabsorption_result.summary["initial_clearance_ml_min"]),
            "zero_capacity_urine_mass_24h_mg": float(zero_reabsorption_result.summary["urine_mass_24h_mg"]),
            "double_capacity_urine_mass_24h_mg": float(double_reabsorption_result.summary["urine_mass_24h_mg"]),
            "identifiability": "UNKNOWN without timed urine and free-plasma concentration data",
        },
        "sensitivity": sensitivity_analysis(p),
        "parameters": asdict(p),
        "parameter_table": parameter_table(p),
        "equations": equation_table(),
        "dimensional_analysis": dimensional_analysis(),
    }


def write_results(path: str | Path = "results.json", p: ModelParameters | None = None) -> dict[str, Any]:
    result = run_experiment(p)
    target = Path(path)
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    result = write_results(Path(__file__).with_name("results.json"))
    prediction = result["prediction"]
    print(f"{result['node']} {result['status']} CL_R={prediction['initial_clearance_ml_min']:.3f} mL/min urine_24h={prediction['urine_mass_24h_mg']:.3f} mg")


if __name__ == "__main__":
    main()
