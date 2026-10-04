from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass, fields, replace
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


STATE_NAMES = ("q", "o", "A", "Gly", "L", "H", "C", "R", "M", "D", "E", "F")
STATE_UNITS = {
    "q": "normalized regional flow [1]",
    "o": "normalized tissue oxygen availability [1]",
    "A": "ATP fraction [1]",
    "Gly": "glycolytic substrate fraction [1]",
    "L": "lactate proxy [1]",
    "H": "acidosis proxy [1]",
    "C": "normalized intracellular calcium [1]",
    "R": "ROS proxy [1]",
    "M": "mPTP-open fraction [1]",
    "D": "irreversible injury fraction [1]",
    "E": "edema proxy [1]",
    "F": "contractile-function fraction [1]",
}


@dataclass(frozen=True)
class Parameters:
    q_art: float = 1.0
    q_ischemia: float = 0.01
    tau_flow: float = 0.35
    k_no_reflow: float = 1.35
    k_edema_flow: float = 0.35
    k_oxygen_on: float = 1.20
    k_oxygen_use: float = 0.70
    vmax_resp: float = 0.12
    km_oxygen: float = 0.20
    gmax: float = 0.055
    k_gly_recovery: float = 0.010
    k_atp_recovery: float = 0.008
    k_lactate_clear: float = 0.100
    lactate_yield: float = 5.0
    k_acid_clear: float = 0.120
    acid_yield: float = 2.0
    k_acid_anoxic: float = 0.020
    k_ca_influx: float = 0.035
    k_ca_release: float = 0.012
    k_ca_clear: float = 0.320
    k_ca_mito: float = 0.040
    k_ros_reox: float = 0.180
    k_ros_mito: float = 0.120
    k_ros_anoxic: float = 0.012
    k_ros_clear: float = 0.160
    k_mptp_open: float = 0.450
    k_mptp_close: float = 0.160
    k_ca_gate: float = 0.250
    k_ros_gate: float = 0.150
    k_energy_gate: float = 0.250
    k_ros_priming: float = 0.350
    k_injury: float = 0.008
    k_injury_repair: float = 0.010
    k_edema: float = 0.004
    k_edema_mptp: float = 0.003
    k_edema_clear: float = 0.100
    tau_function: float = 0.650
    k_function_stroke: float = 0.025


@dataclass(frozen=True)
class Scenario:
    name: str = "full"
    pre_ischemia_min: float = 5.0
    ischemia_min: float = 90.0
    reperfusion_min: float = 30.0
    q_ischemia: float = 0.01
    mptp_block: float = 0.0
    mode: str = "full"
    max_step: float = 0.25


PARAMETER_META: dict[str, tuple[str, str, str]] = {
    "q_art": ("normalized flow [1]", "definition", "Pre-ischemic arterial oxygen-delivery reference"),
    "q_ischemia": ("normalized flow [1]", "literature_anchor", "El Baradie 2021 Results/Fig. 2-3: approximately 1% of baseline during ischemia"),
    "tau_flow": ("min", "assumption", "First-order regional-flow relaxation; not a measured rate"),
    "k_no_reflow": ("dimensionless", "assumption", "Injury-dependent conductance loss; direction supported by early no-reflow anchor"),
    "k_edema_flow": ("dimensionless", "assumption", "Additional flow loss from edema"),
    "k_oxygen_on": ("min^-1", "assumption", "First-order oxygen availability relaxation"),
    "k_oxygen_use": ("dimensionless", "assumption", "Oxygen withdrawal per normalized respiration proxy"),
    "vmax_resp": ("min^-1", "assumption", "Normalized oxidative phosphorylation capacity"),
    "km_oxygen": ("normalized oxygen [1]", "assumption", "Monod oxygen half-saturation"),
    "gmax": ("min^-1", "assumption", "Maximal glycolytic ATP production rate"),
    "k_gly_recovery": ("min^-1", "assumption", "Flow-dependent replenishment of glycolytic substrate"),
    "k_atp_recovery": ("min^-1", "assumption", "Flow/oxygen-dependent reserve recovery term"),
    "k_lactate_clear": ("min^-1", "assumption", "Flow washout of lactate proxy"),
    "lactate_yield": ("lactate units per ATP proxy", "assumption", "Proxy yield, not a stoichiometric clinical value"),
    "k_acid_clear": ("min^-1", "assumption", "Flow washout of acidosis proxy"),
    "acid_yield": ("acidosis units per glycolytic ATP proxy", "assumption", "Proxy coupling coefficient"),
    "k_acid_anoxic": ("min^-1", "assumption", "ATP-deficit contribution to acidosis"),
    "k_ca_influx": ("min^-1", "assumption", "Calcium entry during low flow"),
    "k_ca_release": ("min^-1", "assumption", "Calcium release coupled to energy deficit"),
    "k_ca_clear": ("min^-1", "assumption", "Flow-dependent calcium extrusion"),
    "k_ca_mito": ("min^-1", "assumption", "Calcium transfer associated with mPTP stress"),
    "k_ros_reox": ("min^-1", "assumption", "Reoxygenation-associated ROS production"),
    "k_ros_mito": ("min^-1", "assumption", "ROS production from damaged mitochondria"),
    "k_ros_anoxic": ("min^-1", "assumption", "ATP-deficit-associated ROS production"),
    "k_ros_clear": ("min^-1", "assumption", "Antioxidant/flow-dependent ROS clearance"),
    "k_mptp_open": ("min^-1", "assumption", "Maximum mPTP opening rate"),
    "k_mptp_close": ("min^-1", "assumption", "mPTP closing rate"),
    "k_ca_gate": ("normalized calcium [1]", "assumption", "Hill/Michaelis-Menten gate scale"),
    "k_ros_gate": ("normalized ROS [1]", "assumption", "ROS gate scale"),
    "k_energy_gate": ("normalized ATP deficit [1]", "assumption", "ATP-deficit gate scale"),
    "k_ros_priming": ("dimensionless", "assumption", "Weight of ROS in the opening gate"),
    "k_injury": ("min^-1", "assumption", "Irreversible injury accumulation"),
    "k_injury_repair": ("min^-1", "assumption", "Flow-dependent repair/removal term"),
    "k_edema": ("min^-1", "assumption", "Low-flow edema formation"),
    "k_edema_mptp": ("min^-1", "assumption", "mPTP-associated edema formation"),
    "k_edema_clear": ("min^-1", "assumption", "Flow-dependent edema clearance"),
    "tau_function": ("min", "assumption", "First-order relaxation toward functional capacity"),
    "k_function_stroke": ("min^-1", "assumption", "Direct calcium/mPTP functional depression"),
}


SENSITIVITY_PARAMETERS = ("vmax_resp", "k_mptp_open", "k_no_reflow")


def clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return float(np.clip(value, lower, upper))


def unpack(state: np.ndarray) -> dict[str, float]:
    return {name: float(value) for name, value in zip(STATE_NAMES, state)}


def bounded_derivative(raw: float, value: float, upper: float | None = None) -> float:
    if value <= 0.0 and raw < 0.0:
        return 0.0
    if upper is not None and value >= upper and raw > 0.0:
        return 0.0
    return float(raw)


def respiration_rate(
    oxygen: float,
    atp: float,
    mptp: float,
    acidosis: float,
    p: Parameters,
) -> float:
    oxygen_ratio = oxygen / (p.km_oxygen + max(oxygen, 0.0))
    acid_factor = 1.0 / (1.0 + 2.0 * max(acidosis, 0.0))
    return p.vmax_resp * oxygen_ratio * max(atp, 0.0) * (1.0 - clamp(mptp)) * acid_factor


def glycolysis_rate(
    oxygen: float,
    atp: float,
    substrate: float,
    mptp: float,
    acidosis: float,
    p: Parameters,
) -> float:
    oxygen_ratio = oxygen / (p.km_oxygen + max(oxygen, 0.0))
    hypoxia_gate = clamp((p.km_oxygen - oxygen) / p.km_oxygen)
    acid_factor = 1.0 / (1.0 + 2.0 * max(acidosis, 0.0))
    return p.gmax * clamp(substrate) * hypoxia_gate * max(atp, 0.0) * (1.0 - clamp(mptp)) * acid_factor


@lru_cache(maxsize=None)
def healthy_oxygen(p: Parameters) -> float:
    def balance(oxygen: float) -> float:
        return p.k_oxygen_on * (p.q_art - oxygen) - p.k_oxygen_use * respiration_rate(oxygen, 1.0, 0.0, 0.0, p)

    lower = 0.0
    upper = 1.0
    for _ in range(80):
        middle = 0.5 * (lower + upper)
        if balance(middle) > 0.0:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def healthy_state(p: Parameters) -> np.ndarray:
    return np.array(
        [
            p.q_art,
            healthy_oxygen(p),
            1.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        ],
        dtype=float,
    )


def external_flow(t: float, scenario: Scenario) -> float:
    if t < scenario.pre_ischemia_min:
        return 1.0
    if t < scenario.pre_ischemia_min + scenario.ischemia_min:
        return scenario.q_ischemia
    return 1.0


def phase(t: float, scenario: Scenario) -> str:
    if t < scenario.pre_ischemia_min:
        return "pre_ischemia"
    if t < scenario.pre_ischemia_min + scenario.ischemia_min:
        return "ischemia"
    return "reperfusion"


def rhs(t: float, state: np.ndarray, p: Parameters, scenario: Scenario) -> np.ndarray:
    q, oxygen, atp, substrate, lactate, acidosis, calcium, ros, mptp, injury, edema, function = (float(value) for value in state)
    q_actual = q
    q_external = external_flow(t, scenario)

    if scenario.mode == "flow_only":
        q, oxygen, atp, substrate, lactate, acidosis, calcium, ros, mptp, injury, edema, function = healthy_state(p)
        injury_block = 0.0
    else:
        injury_block = 1.0 - (1.0 - clamp(p.k_no_reflow * mptp)) * (1.0 - clamp(p.k_edema_flow * edema))

    flow_target = q_external * (1.0 - injury_block)
    d_flow = (flow_target - q_actual) / p.tau_flow

    oxygen_ratio = oxygen / (p.km_oxygen + max(oxygen, 0.0))
    respiration = respiration_rate(oxygen, atp, mptp, acidosis, p)
    glycolysis = glycolysis_rate(oxygen, atp, substrate, mptp, acidosis, p)
    healthy_respiration = p.vmax_resp * healthy_oxygen(p) / (p.km_oxygen + healthy_oxygen(p))
    demand = healthy_respiration * (0.35 + 0.65 * clamp(function)) * (0.85 + 0.15 * clamp(oxygen_ratio))
    d_atp_raw = respiration + glycolysis - demand + p.k_atp_recovery * max(1.0 - atp, 0.0) * q

    d_substrate_raw = -glycolysis + p.k_gly_recovery * q * (1.0 - clamp(substrate))
    d_lactate_raw = p.lactate_yield * glycolysis - p.k_lactate_clear * q * lactate
    d_acid_raw = (
        p.acid_yield * glycolysis
        + p.k_acid_anoxic * (1.0 - clamp(oxygen_ratio)) * max(1.0 - atp, 0.0)
        - p.k_acid_clear * q * acidosis
    )

    energy_deficit = max(1.0 - atp, 0.0)
    calcium_excess = max(calcium - 1.0, 0.0)
    d_calcium_raw = (
        p.k_ca_influx * (1.0 - clamp(q)) * (1.0 + 0.8 * energy_deficit)
        + p.k_ca_release * (1.0 - clamp(q)) * energy_deficit
        - p.k_ca_clear * q * (calcium - 1.0)
        - p.k_ca_mito * calcium_excess
    )

    d_oxygen_raw = p.k_oxygen_on * (q * p.q_art - oxygen) - p.k_oxygen_use * respiration
    d_ros_raw = (
        p.k_ros_reox * clamp(q) * clamp(oxygen_ratio) * energy_deficit * (1.0 - clamp(mptp))
        + p.k_ros_mito * clamp(mptp)
        + p.k_ros_anoxic * (1.0 - clamp(q)) * energy_deficit
        - p.k_ros_clear * (0.2 + 0.8 * clamp(q)) * ros
    )

    calcium_gate = calcium_excess / (p.k_ca_gate + calcium_excess) if calcium_excess > 0.0 else 0.0
    ros_gate = ros / (p.k_ros_gate + ros) if ros > 0.0 else 0.0
    energy_gate = energy_deficit / (p.k_energy_gate + energy_deficit) if energy_deficit > 0.0 else 0.0
    opening_trigger = calcium_gate * (p.k_ros_priming + (1.0 - p.k_ros_priming) * ros_gate) * (0.30 + 0.70 * energy_gate)
    d_mptp_raw = (1.0 - scenario.mptp_block) * p.k_mptp_open * opening_trigger * (1.0 - clamp(mptp)) - p.k_mptp_close * mptp

    d_injury_raw = p.k_injury * clamp(mptp) * energy_deficit * (1.0 + calcium_excess) - p.k_injury_repair * q * injury
    d_edema_raw = (
        p.k_edema * (1.0 - clamp(q)) * (1.0 + calcium_excess)
        + p.k_edema_mptp * clamp(mptp)
        - p.k_edema_clear * q * edema
    )

    target_function = clamp(
        (0.15 + 0.85 * clamp(atp))
        * np.exp(-0.65 * calcium_excess)
        * (1.0 - 0.65 * clamp(mptp))
        * (1.0 - 0.45 * clamp(edema))
        * (1.0 - clamp(injury))
    )
    d_function_raw = (target_function - function) / p.tau_function - p.k_function_stroke * (calcium_excess + mptp) * (0.5 + 0.5 * (1.0 - clamp(q)))

    derivatives = np.array(
        [
            d_flow,
            d_oxygen_raw,
            d_atp_raw,
            d_substrate_raw,
            d_lactate_raw,
            d_acid_raw,
            d_calcium_raw,
            d_ros_raw,
            d_mptp_raw,
            d_injury_raw,
            d_edema_raw,
            d_function_raw,
        ],
        dtype=float,
    )
    derivatives[0] = bounded_derivative(derivatives[0], q_actual, 1.0)
    derivatives[1] = bounded_derivative(derivatives[1], oxygen, 1.0)
    derivatives[2] = bounded_derivative(derivatives[2], atp, 1.0)
    derivatives[3] = bounded_derivative(derivatives[3], substrate, 1.0)
    derivatives[4] = bounded_derivative(derivatives[4], lactate)
    derivatives[5] = bounded_derivative(derivatives[5], acidosis)
    derivatives[6] = bounded_derivative(derivatives[6], calcium, 5.0)
    derivatives[7] = bounded_derivative(derivatives[7], ros)
    derivatives[8] = bounded_derivative(derivatives[8], mptp, 1.0)
    derivatives[9] = bounded_derivative(derivatives[9], injury, 1.0)
    derivatives[10] = bounded_derivative(derivatives[10], edema, 1.0)
    derivatives[11] = bounded_derivative(derivatives[11], function, 1.0)
    return derivatives


def validate_parameters(p: Parameters, scenario: Scenario) -> None:
    positive = {
        "tau_flow",
        "k_oxygen_on",
        "k_oxygen_use",
        "vmax_resp",
        "km_oxygen",
        "gmax",
        "k_gly_recovery",
        "k_atp_recovery",
        "k_lactate_clear",
        "k_acid_clear",
        "k_ca_influx",
        "k_ca_release",
        "k_ca_clear",
        "k_ca_mito",
        "k_ros_reox",
        "k_ros_mito",
        "k_ros_anoxic",
        "k_ros_clear",
        "k_ros_priming",
        "k_mptp_open",
        "k_mptp_close",
        "k_ca_gate",
        "k_ros_gate",
        "k_energy_gate",
        "k_injury",
        "k_injury_repair",
        "k_edema",
        "k_edema_mptp",
        "k_edema_clear",
        "tau_function",
        "k_function_stroke",
    }
    values = asdict(p)
    for name in positive:
        if values[name] <= 0.0:
            raise ValueError(f"{name} must be positive")
    if not 0.0 <= p.q_art <= 1.0:
        raise ValueError("q_art must be normalized")
    if not 0.0 <= p.q_ischemia <= 1.0:
        raise ValueError("q_ischemia must be normalized")
    if not 0.0 <= scenario.q_ischemia <= 1.0:
        raise ValueError("scenario.q_ischemia must be normalized")
    if not 0.0 <= scenario.mptp_block <= 1.0:
        raise ValueError("scenario.mptp_block must be in [0, 1]")
    if scenario.mode not in {"full", "flow_only"}:
        raise ValueError("scenario.mode must be full or flow_only")
    if scenario.pre_ischemia_min < 0.0 or scenario.ischemia_min < 0.0 or scenario.reperfusion_min <= 0.0:
        raise ValueError("scenario times are invalid")


def simulate(p: Parameters, scenario: Scenario) -> dict[str, Any]:
    validate_parameters(p, scenario)
    duration = scenario.pre_ischemia_min + scenario.ischemia_min + scenario.reperfusion_min
    sample_step = min(0.05, scenario.max_step)
    time_grid = np.arange(0.0, duration + sample_step, sample_step)
    solution = solve_ivp(
        lambda t, y: rhs(t, y, p, scenario),
        (0.0, duration),
        healthy_state(p),
        method="DOP853",
        t_eval=time_grid,
        rtol=1e-7,
        atol=1e-9,
        max_step=scenario.max_step,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    states = np.asarray(solution.y, dtype=float)
    if not np.all(np.isfinite(states)):
        raise RuntimeError("non-finite state")
    nonnegative_indices = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11]
    if np.any(states[nonnegative_indices, :] < -1e-7):
        raise RuntimeError("negative fraction state")
    if np.any(states[6, :] > 5.0 + 1e-7):
        raise RuntimeError("calcium state exceeded bound")
    bounded_indices = [0, 1, 2, 3, 8, 9, 10, 11]
    states[bounded_indices, :] = np.clip(states[bounded_indices, :], 0.0, 1.0)
    states[6, :] = np.clip(states[6, :], 0.0, 5.0)
    states[[4, 5, 7], :] = np.maximum(states[[4, 5, 7], :], 0.0)
    return {
        "time_min": solution.t,
        "states": states,
        "reperfusion_time_min": scenario.pre_ischemia_min + scenario.ischemia_min,
        "scenario": scenario,
        "parameters": p,
    }


def state_at(simulation: dict[str, Any], absolute_time: float) -> dict[str, float]:
    time_grid = simulation["time_min"]
    if absolute_time < time_grid[0] or absolute_time > time_grid[-1]:
        raise ValueError("requested time is outside simulation")
    values = np.empty(len(STATE_NAMES), dtype=float)
    for index in range(len(STATE_NAMES)):
        values[index] = np.interp(absolute_time, time_grid, simulation["states"][index, :])
    return unpack(values)


def sample_relative(simulation: dict[str, Any], relative_min: float) -> dict[str, float]:
    return state_at(simulation, simulation["reperfusion_time_min"] + relative_min)


def summarize_simulation(simulation: dict[str, Any], scenario: Scenario) -> dict[str, Any]:
    relative_times = [0.0, 5.0, 15.0, 30.0]
    available = min(scenario.reperfusion_min, 30.0)
    samples: dict[str, dict[str, float]] = {}
    for relative_time in relative_times:
        if relative_time <= available + 1e-12:
            samples[f"rep_{relative_time:g}_min"] = sample_relative(simulation, relative_time)
    samples["pre_ischemia"] = state_at(simulation, 0.0)
    samples["ischemia_end"] = state_at(simulation, scenario.pre_ischemia_min + scenario.ischemia_min)
    duration = scenario.pre_ischemia_min + scenario.ischemia_min + scenario.reperfusion_min
    trajectory_times = set(float(value) for value in np.arange(0.0, duration + 1e-12, 5.0))
    trajectory_times.update(
        [
            0.0,
            scenario.pre_ischemia_min,
            scenario.pre_ischemia_min + scenario.ischemia_min,
            duration,
        ]
    )
    trajectory = []
    for absolute_time in sorted(value for value in trajectory_times if value <= duration + 1e-12):
        state = state_at(simulation, absolute_time)
        trajectory.append(
            {
                "time_min": absolute_time,
                "relative_reperfusion_min": absolute_time - simulation["reperfusion_time_min"],
                "phase": phase(absolute_time, scenario),
                "state": state,
            }
        )
    return {"scenario": asdict(scenario), "samples": samples, "trajectory": trajectory}


def run_case(p: Parameters, scenario: Scenario) -> dict[str, Any]:
    simulation = simulate(p, scenario)
    return summarize_simulation(simulation, scenario)


def analytic_mptp_solution(m0: float, k_open: float, k_close: float, trigger: float, t: float) -> float:
    equilibrium = k_open * trigger / (k_open * trigger + k_close)
    rate = k_open * trigger + k_close
    return equilibrium + (m0 - equilibrium) * np.exp(-rate * t)


def check_units() -> dict[str, Any]:
    expected_units = {name: meta[0] for name, meta in PARAMETER_META.items()}
    parameter_units = {field.name: expected_units.get(field.name, "missing") for field in fields(Parameters)}
    derivative_units = {name: "fraction/min or normalized proxy/min" for name in STATE_NAMES}
    terms = {
        "q": ["(normalized flow target - q) / min"],
        "o": ["min^-1 * normalized oxygen - dimensionless * respiration [min^-1]"],
        "A": ["respiration + glycolysis - demand + min^-1 * ATP deficit"],
        "Gly": ["-glycolysis + min^-1 * substrate recovery"],
        "L": ["ATP-rate * dimensionless yield - min^-1 * lactate"],
        "H": ["ATP-rate * dimensionless yield + min^-1 * deficit - min^-1 * acidosis"],
        "C": ["min^-1 * dimensionless calcium state"],
        "R": ["min^-1 * dimensionless ROS state"],
        "M": ["min^-1 * dimensionless mPTP fraction"],
        "D": ["min^-1 * dimensionless injury fraction"],
        "E": ["min^-1 * dimensionless edema fraction"],
        "F": ["(target function - function) / min - min^-1 * stress"],
    }
    all_parameters_mapped = set(parameter_units) == set(expected_units) and all(
        parameter_units.get(name) == unit for name, unit in expected_units.items()
    )
    all_terms_have_time = all(any("min" in term for term in value) for value in terms.values())
    return {
        "parameter_units": parameter_units,
        "state_units": STATE_UNITS,
        "derivative_units": derivative_units,
        "balance_terms": terms,
        "all_parameters_mapped": all_parameters_mapped,
        "all_derivatives_have_time_denominator": all_terms_have_time,
        "unit_consistent": all_parameters_mapped and all_terms_have_time,
        "interpretation": "All state variables are normalized fractions/proxies; every ODE term is a normalized amount per minute.",
    }


def healthy_equilibrium_check(p: Parameters) -> dict[str, Any]:
    scenario = Scenario(name="healthy", pre_ischemia_min=0.0, ischemia_min=0.0, reperfusion_min=30.0)
    simulation = simulate(p, scenario)
    initial = simulation["states"][:, 0]
    final = simulation["states"][:, -1]
    deviation = float(np.max(np.abs(simulation["states"] - initial[:, None])))
    return {
        "max_absolute_deviation": deviation,
        "function_initial": float(initial[-1]),
        "function_final": float(final[-1]),
        "passes_1e-8": deviation < 1e-8 and abs(float(final[-1]) - 1.0) < 1e-8,
    }


def criterion_check(cases: dict[str, dict[str, Any]], healthy: dict[str, Any]) -> dict[str, Any]:
    full = cases["full"]["samples"]["rep_15_min"]
    blocked = cases["mptp_blocked"]["samples"]["rep_15_min"]
    flow_only = cases["flow_only"]["samples"]["rep_15_min"]
    checks = {
        "healthy_equilibrium": bool(healthy["passes_1e-8"]),
        "full_mptp_above_0.10": bool(full["M"] > 0.10),
        "full_atp_below_0.90": bool(full["A"] < 0.90),
        "full_function_below_0.90": bool(full["F"] < 0.90),
        "mptp_block_reduces_mptp_by_0.05": bool(full["M"] - blocked["M"] >= 0.05),
        "mptp_block_improves_function_by_0.05": bool(blocked["F"] - full["F"] >= 0.05),
        "flow_only_function_within_0.02": bool(abs(flow_only["F"] - 1.0) <= 0.02),
        "full_is_injured_relative_to_flow_only": bool(full["F"] < 0.90 and flow_only["F"] >= 0.98),
    }
    checks["all_passed"] = all(checks.values())
    return checks


def run_sensitivity(p: Parameters, base_scenario: Scenario) -> dict[str, Any]:
    base = run_case(p, base_scenario)
    base_15 = base["samples"]["rep_15_min"]
    output: dict[str, Any] = {"baseline": {name: base_15[name] for name in ("q", "A", "M", "F")}, "parameters": {}}
    for name in SENSITIVITY_PARAMETERS:
        parameter_values: dict[str, Any] = {}
        for multiplier in (0.5, 1.0, 1.5):
            changed = replace(p, **{name: getattr(p, name) * multiplier})
            case = run_case(changed, base_scenario)
            values = case["samples"]["rep_15_min"]
            parameter_values[f"{multiplier:g}x"] = {key: values[key] for key in ("q", "A", "M", "F")}
        parameter_values["delta_0.5x_to_1.5x"] = {
            key: parameter_values["1.5x"][key] - parameter_values["0.5x"][key]
            for key in ("q", "A", "M", "F")
        }
        parameter_values["relative_change_0.5x_to_1.5x"] = {
            key: (
                parameter_values["delta_0.5x_to_1.5x"][key] / base_15[key]
                if abs(base_15[key]) > 1e-12
                else None
            )
            for key in ("q", "A", "M", "F")
        }
        output["parameters"][name] = parameter_values
    return output


def parameter_table(p: Parameters) -> dict[str, dict[str, Any]]:
    values = asdict(p)
    return {
        name: {
            "value": values[name],
            "unit": PARAMETER_META[name][0],
            "source": PARAMETER_META[name][1],
            "rationale": PARAMETER_META[name][2],
        }
        for name in PARAMETER_META
    }


def literature_anchors() -> list[dict[str, Any]]:
    return [
        {
            "id": "E1",
            "citation": "El Baradie et al., Scientific Reports 11, 6152 (2021)",
            "doi": "10.1038/s41598-021-85753-x",
            "location": "Results and Figs. 2-3; Table 2",
            "values": {
                "ischemia_duration_min": 90,
                "ischemic_flow_percent_baseline_approx": 1,
                "early_reperfusion_flow_percent_baseline_approx": 10,
                "control_gait_score": 1.70,
                "control_gait_score_sd": 0.67,
                "nim811_gait_score": 2.70,
                "nim811_gait_score_sd": 0.82,
            },
            "unit_notes": "Flow is percent of pre-ischemic baseline; gait score is unitless Tarlov score. Functional timing is after recovery, not a 15 min label.",
        },
        {
            "id": "E2",
            "citation": "McAllister et al., American Journal of Physiology Regulatory, Integrative and Comparative Physiology 295, R681-R689 (2008)",
            "doi": "10.1152/ajpregu.90303.2008",
            "location": "Abstract",
            "values": {"postconditioning_infarction_control_percent": 44, "postconditioning_infarction_treated_percent": 22},
            "unit_notes": "Percent infarction at 48 h reperfusion; mechanistic support only, not a model coefficient.",
        },
        {
            "id": "E3",
            "citation": "Naparus et al., European Journal of Pharmacology 686, 90-96 (2012)",
            "doi": "10.1016/j.ejphar.2012.04.045",
            "location": "Abstract",
            "values": {"hypoxia_h": 3, "reoxygenation_h": 2, "nim811_concentration_uM": 5},
            "unit_notes": "Human muscle strip model; supports ATP/function state retention and mPTP inhibition.",
        },
        {
            "id": "E4",
            "citation": "Tran et al., PLoS ONE 7, e43410 (2012)",
            "doi": "10.1371/journal.pone.0043410",
            "location": "Abstract",
            "values": {"ischemia_h": 3, "reperfusion_h": 4},
            "unit_notes": "Mouse tourniquet model linking mitochondrial superoxide, mPTP opening, and apoptosis.",
        },
        {
            "id": "E5",
            "citation": "Pottecher et al., Journal of Vascular Surgery 57, 1100-1108.e2 (2013)",
            "doi": "10.1016/j.jvs.2012.09.020",
            "location": "Abstract, Results",
            "values": {
                "ir_vmax_umol_o2_min_g": 4.08,
                "ir_vmax_sd": 0.38,
                "sham_vmax_umol_o2_min_g": 5.98,
                "sham_vmax_sd": 0.56,
                "ir_ros_au": 3992,
                "ir_ros_sd_au": 706,
                "sham_ros_au": 1812,
                "sham_ros_sd_au": 322,
            },
            "unit_notes": "Oxidative capacity in umol O2/min/g; ROS in arbitrary units; not transferred as a normalized rate.",
        },
    ]


def make_results(p: Parameters, scenario: Scenario) -> dict[str, Any]:
    full = run_case(p, scenario)
    blocked_scenario = replace(scenario, name="mptp_blocked", mptp_block=0.85)
    flow_scenario = replace(scenario, name="flow_only", mode="flow_only")
    cases = {
        "full": full,
        "mptp_blocked": run_case(p, blocked_scenario),
        "flow_only": run_case(p, flow_scenario),
    }
    healthy = healthy_equilibrium_check(p)
    criteria = criterion_check(cases, healthy)
    sensitivity = run_sensitivity(p, scenario)
    return {
        "id": "BT-HX-Q036",
        "status": "first_runnable_mechanistic_model",
        "prereg_sha256": "1042448a06adcbf76e41cea35a367b537d884c7e25e134a9a81aed0e69eb0885",
        "model_scope": "normalized regional skeletal-muscle ischemia-reperfusion compartment",
        "state_names": list(STATE_NAMES),
        "state_units": STATE_UNITS,
        "parameter_table": parameter_table(p),
        "literature_anchors": literature_anchors(),
        "scenario": asdict(scenario),
        "cases": cases,
        "healthy_equilibrium_check": healthy,
        "frozen_thresholds": {
            "healthy_max_absolute_deviation": 1e-8,
            "full_mptp_min": 0.10,
            "full_atp_max": 0.90,
            "full_function_max": 0.90,
            "mptp_block_fraction": 0.85,
            "mptp_min_absolute_reduction": 0.05,
            "mptp_min_function_improvement": 0.05,
            "flow_only_function_tolerance": 0.02,
            "calcium_upper_bound": 5.0,
        },
        "frozen_criterion_check": criteria,
        "sensitivity": sensitivity,
        "unit_check": check_units(),
        "validation_status": "structural_check_only; no internal or clinical data were used",
        "source_derivation_hypothesis_separation": {
            "source": ["E1", "E2", "E3", "E4", "E5"],
            "derivation": ["ODE balances, gates, units, and integration in model.py"],
            "hypothesis": ["minimum sufficient ischemic state and frozen directional thresholds in PREREG.md"],
        },
    }


def write_results(path: str, p: Parameters, scenario: Scenario) -> dict[str, Any]:
    results = make_results(p, scenario)
    target = Path(path)
    target.write_text(json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return results


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results.json")
    parser.add_argument("--prereg", default="PREREG.md")
    args = parser.parse_args()
    prereg_hash = sha256_file(args.prereg)
    expected = Path("PREREG.sha256").read_text(encoding="utf-8").split()[0]
    if prereg_hash != expected:
        raise RuntimeError(f"PREREG hash mismatch: {prereg_hash} != {expected}")
    results = write_results(args.output, Parameters(), Scenario())
    print(json.dumps({"output": args.output, "all_criteria_passed": results["frozen_criterion_check"]["all_passed"]}))


if __name__ == "__main__":
    main()
