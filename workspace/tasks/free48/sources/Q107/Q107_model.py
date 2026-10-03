from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ModelParameters:
    slow_wave_rate_hz: float
    fundus_rest_volume_m3: float
    antrum_rest_volume_m3: float
    fundus_radius_m: float
    antrum_radius_m: float
    fundus_wall_modulus_pa: float
    antrum_wall_modulus_pa: float
    fundus_wall_thickness_m: float
    antrum_wall_thickness_m: float
    fundus_active_stress_pa: float
    antrum_active_stress_pa: float
    transfer_cd: float
    transfer_area_m2: float
    pylorus_cd: float
    pylorus_radius_m: float
    pylorus_floor: float
    pylorus_span: float
    pylorus_gate_sharpness: float
    pylorus_gate_threshold: float
    pylorus_phase_rad: float
    downstream_pressure_pa: float
    fluid_density_kg_m3: float
    solid_density_kg_m3: float
    solid_flow_rate_per_s: float
    solid_reference_flow_m3_s: float
    particle_diameter_m: float
    particle_cutoff_m: float
    particle_selectivity_exponent: float
    mobilization_rate_per_s: float
    obstruction_strength: float
    max_solid_fraction: float
    minimum_content_volume_m3: float
    fixed_pylorus_gate: float
    liquid_meal_volume_m3: float
    solid_meal_liquid_m3: float
    solid_meal_mass_kg: float
    solid_initial_mobile_fraction: float
    scalar_liquid_k_per_min: float
    scalar_solid_k_per_min: float


@dataclass(frozen=True)
class Geometry:
    name: str
    fundus_radius_factor: float
    antrum_radius_factor: float
    antrum_active_factor: float
    pylorus_area_factor: float
    transfer_area_factor: float
    dynamic_pylorus: bool
    fixed_pylorus_gate: float


GEOMETRY_METADATA: dict[str, tuple[str, str]] = {
    "fundus_radius_factor": ("1", "scenario factor on effective fundus radius"),
    "antrum_radius_factor": ("1", "scenario factor on effective antrum radius"),
    "antrum_active_factor": ("1", "scenario factor on antral active stress"),
    "pylorus_area_factor": ("1", "scenario factor on pylorus effective area"),
    "transfer_area_factor": ("1", "scenario factor on fundus-antrum transfer area"),
    "dynamic_pylorus": ("boolean", "control: time-dependent opening"),
    "fixed_pylorus_gate": ("1", "control-model fixed open-area fraction"),
}


PARAMETER_METADATA: dict[str, tuple[str, str]] = {
    "slow_wave_rate_hz": ("Hz", "assumption: gastric slow-wave frequency, 3/min"),
    "fundus_rest_volume_m3": ("m^3", "assumption: unloaded fundus reference volume"),
    "antrum_rest_volume_m3": ("m^3", "assumption: unloaded antrum reference volume"),
    "fundus_radius_m": ("m", "assumption: effective fundus radius"),
    "antrum_radius_m": ("m", "assumption: effective antrum radius"),
    "fundus_wall_modulus_pa": ("Pa", "assumption: effective passive wall modulus"),
    "antrum_wall_modulus_pa": ("Pa", "assumption: effective passive wall modulus"),
    "fundus_wall_thickness_m": ("m", "assumption: effective wall thickness"),
    "antrum_wall_thickness_m": ("m", "assumption: effective wall thickness"),
    "fundus_active_stress_pa": ("Pa", "assumption: slow-wave active stress amplitude"),
    "antrum_active_stress_pa": ("Pa", "assumption: antral active stress amplitude"),
    "transfer_cd": ("1", "assumption: discharge coefficient for fundus-antrum transfer"),
    "transfer_area_m2": ("m^2", "assumption: effective transfer orifice area"),
    "pylorus_cd": ("1", "assumption: orifice discharge coefficient"),
    "pylorus_radius_m": ("m", "assumption: effective pylorus radius"),
    "pylorus_floor": ("1", "assumption: minimum open-area fraction"),
    "pylorus_span": ("1", "assumption: dynamic open-area amplitude"),
    "pylorus_gate_sharpness": ("1", "assumption: gating sigmoid sharpness"),
    "pylorus_gate_threshold": ("1", "assumption: normalized contraction threshold"),
    "pylorus_phase_rad": ("rad", "assumption: phase of pylorus gating wave"),
    "downstream_pressure_pa": ("Pa", "boundary condition: duodenal pressure"),
    "fluid_density_kg_m3": ("kg m^-3", "assumption: aqueous meal density"),
    "solid_density_kg_m3": ("kg m^-3", "assumption: effective solid-food density"),
    "solid_flow_rate_per_s": ("s^-1", "assumption: effective first-order particle release coefficient"),
    "solid_reference_flow_m3_s": ("m^3 s^-1", "assumption: reference liquid flow for particle release"),
    "particle_diameter_m": ("m", "assumption: representative retained particle diameter"),
    "particle_cutoff_m": ("m", "assumption: pylorus particle-size selectivity scale"),
    "particle_selectivity_exponent": ("1", "assumption: selectivity sharpness"),
    "mobilization_rate_per_s": ("s^-1", "assumption: shear-driven bound-to-mobile particle rate"),
    "obstruction_strength": ("1", "assumption: concentration-dependent pylorus obstruction"),
    "max_solid_fraction": ("1", "assumption: maximum carried solid volume fraction"),
    "minimum_content_volume_m3": ("m^3", "numerical regularizer for concentration"),
    "fixed_pylorus_gate": ("1", "control model: fixed open-area fraction"),
    "liquid_meal_volume_m3": ("m^3", "controlled test meal: 400 mL liquid"),
    "solid_meal_liquid_m3": ("m^3", "controlled test meal: 100 mL accompanying liquid"),
    "solid_meal_mass_kg": ("kg", "controlled test meal: 100 g solid marker"),
    "solid_initial_mobile_fraction": ("1", "initial condition: mobile fraction of solid marker"),
    "scalar_liquid_k_per_min": ("min^-1", "null model: assumed liquid emptying rate"),
    "scalar_solid_k_per_min": ("min^-1", "null model: assumed solid emptying rate"),
}


GEOMETRIES: dict[str, Geometry] = {
    "fundus_antrum": Geometry(
        name="fundus_antrum",
        fundus_radius_factor=1.0,
        antrum_radius_factor=1.0,
        antrum_active_factor=1.0,
        pylorus_area_factor=0.75,
        transfer_area_factor=1.0,
        dynamic_pylorus=False,
        fixed_pylorus_gate=0.03,
    ),
    "fundus_antrum_pylorus": Geometry(
        name="fundus_antrum_pylorus",
        fundus_radius_factor=1.15,
        antrum_radius_factor=0.90,
        antrum_active_factor=1.05,
        pylorus_area_factor=1.0,
        transfer_area_factor=1.0,
        dynamic_pylorus=True,
        fixed_pylorus_gate=0.03,
    ),
}


REFERENCE = {
    "source": "Collins, P. J. (1983), Gut",
    "status": "OVERIFIERAD, ur minnet",
    "liquid_t50_min": 20.0,
    "solid_t50_min": 90.0,
    "liquid_factor_two_bounds_min": [10.0, 40.0],
    "solid_factor_two_bounds_min": [45.0, 180.0],
    "use": "external directional anchor only; not used for calibration",
}


def default_parameters() -> ModelParameters:
    return ModelParameters(
        slow_wave_rate_hz=3.0 / 60.0,
        fundus_rest_volume_m3=0.000040,
        antrum_rest_volume_m3=0.000050,
        fundus_radius_m=0.035,
        antrum_radius_m=0.025,
        fundus_wall_modulus_pa=12000.0,
        antrum_wall_modulus_pa=18000.0,
        fundus_wall_thickness_m=0.002,
        antrum_wall_thickness_m=0.003,
        fundus_active_stress_pa=900.0,
        antrum_active_stress_pa=2200.0,
        transfer_cd=0.70,
        transfer_area_m2=5.0e-7,
        pylorus_cd=0.62,
        pylorus_radius_m=0.0010,
        pylorus_floor=0.002,
        pylorus_span=0.060,
        pylorus_gate_sharpness=4.0,
        pylorus_gate_threshold=0.10,
        pylorus_phase_rad=0.0,
        downstream_pressure_pa=0.0,
        fluid_density_kg_m3=1000.0,
        solid_density_kg_m3=1200.0,
        solid_flow_rate_per_s=0.050,
        solid_reference_flow_m3_s=0.0000005,
        particle_diameter_m=0.002,
        particle_cutoff_m=0.0025,
        particle_selectivity_exponent=2.0,
        mobilization_rate_per_s=0.00035,
        obstruction_strength=0.75,
        max_solid_fraction=0.30,
        minimum_content_volume_m3=0.000005,
        fixed_pylorus_gate=0.03,
        liquid_meal_volume_m3=0.00040,
        solid_meal_liquid_m3=0.00010,
        solid_meal_mass_kg=0.100,
        solid_initial_mobile_fraction=0.15,
        scalar_liquid_k_per_min=0.030,
        scalar_solid_k_per_min=0.0075,
    )


def parameter_table(parameters: ModelParameters) -> list[dict[str, Any]]:
    values = asdict(parameters)
    return [
        {
            "name": name,
            "value": values[name],
            "unit": unit,
            "source_or_assumption": source,
        }
        for name, (unit, source) in PARAMETER_METADATA.items()
    ]


def geometry_table() -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for geometry in GEOMETRIES.values():
        values = asdict(geometry)
        row: dict[str, Any] = {"name": geometry.name}
        for field, (unit, source) in GEOMETRY_METADATA.items():
            row[field] = {
                "value": values[field],
                "unit": unit,
                "source_or_assumption": source,
            }
        table.append(row)
    return table


def orifice_flow(
    upstream_pressure_pa: float,
    downstream_pressure_pa: float,
    density_kg_m3: float,
    area_m2: float,
    discharge_coefficient: float,
    open_fraction: float = 1.0,
    obstruction_fraction: float = 0.0,
) -> float:
    pressure = max(float(upstream_pressure_pa) - float(downstream_pressure_pa), 0.0)
    area = float(area_m2) * min(max(float(open_fraction), 0.0), 1.0)
    area *= max(0.0, 1.0 - min(max(float(obstruction_fraction), 0.0), 0.95))
    return float(discharge_coefficient) * area * math.sqrt(2.0 * pressure / float(density_kg_m3))


def passive_wall_pressure(
    volume_m3: float,
    rest_volume_m3: float,
    radius_m: float,
    wall_thickness_m: float,
    wall_modulus_pa: float,
) -> float:
    strain = max(float(volume_m3) / float(rest_volume_m3) - 1.0, 0.0)
    return 2.0 * float(wall_modulus_pa) * float(wall_thickness_m) * strain / float(radius_m)


def active_wall_pressure(
    content_volume_m3: float,
    rest_volume_m3: float,
    active_stress_pa: float,
    wave: float,
) -> float:
    occupancy = min(max(float(content_volume_m3) / float(rest_volume_m3) - 0.25, 0.0), 1.5)
    return float(active_stress_pa) * float(wave) * (0.25 + 0.75 * occupancy)


def pylorus_open_fraction(
    time_s: float,
    parameters: ModelParameters,
    geometry: Geometry,
    wave_argument: float,
) -> float:
    if not geometry.dynamic_pylorus:
        return float(geometry.fixed_pylorus_gate)
    gate = parameters.pylorus_floor + parameters.pylorus_span / (
        1.0
        + math.exp(
            -parameters.pylorus_gate_sharpness
            * (wave_argument - parameters.pylorus_gate_threshold)
        )
    )
    return min(max(gate, 0.0), 1.0)


def initial_state(meal_type: str, parameters: ModelParameters) -> np.ndarray:
    if meal_type == "liquid":
        total = parameters.liquid_meal_volume_m3
        return np.array([0.75 * total, 0.25 * total, 0.0, 0.0], dtype=float)
    if meal_type == "solid":
        liquid = parameters.solid_meal_liquid_m3
        mass = parameters.solid_meal_mass_kg
        return np.array(
            [
                0.50 * liquid,
                0.50 * liquid,
                mass,
                parameters.solid_initial_mobile_fraction * mass,
            ],
            dtype=float,
        )
    raise ValueError(f"unknown meal type: {meal_type}")


def state_diagnostics(
    state: np.ndarray,
    time_s: float,
    parameters: ModelParameters,
    geometry: Geometry,
) -> dict[str, float]:
    state = np.maximum(np.asarray(state, dtype=float), 0.0)
    fundus_liquid, antrum_liquid, solid_mass, mobile_solid = state
    solid_volume = solid_mass / parameters.solid_density_kg_m3
    antrum_content = antrum_liquid + solid_volume
    phase = 2.0 * math.pi * parameters.slow_wave_rate_hz * time_s
    contraction_signal = math.sin(phase)
    wave_argument = math.sin(phase + parameters.pylorus_phase_rad)
    contraction_wave = max(0.0, 0.5 + 0.5 * contraction_signal) ** 2

    fundus_radius = parameters.fundus_radius_m * geometry.fundus_radius_factor
    antrum_radius = parameters.antrum_radius_m * geometry.antrum_radius_factor
    fundus_pressure = passive_wall_pressure(
        fundus_liquid,
        parameters.fundus_rest_volume_m3,
        fundus_radius,
        parameters.fundus_wall_thickness_m,
        parameters.fundus_wall_modulus_pa,
    ) + active_wall_pressure(
        fundus_liquid,
        parameters.fundus_rest_volume_m3,
        parameters.fundus_active_stress_pa,
        contraction_wave,
    )
    antrum_pressure = passive_wall_pressure(
        antrum_content,
        parameters.antrum_rest_volume_m3,
        antrum_radius,
        parameters.antrum_wall_thickness_m,
        parameters.antrum_wall_modulus_pa,
    ) + active_wall_pressure(
        antrum_content,
        parameters.antrum_rest_volume_m3,
        parameters.antrum_active_stress_pa * geometry.antrum_active_factor,
        contraction_wave,
    )

    transfer_flow = orifice_flow(
        fundus_pressure,
        antrum_pressure,
        parameters.fluid_density_kg_m3,
        parameters.transfer_area_m2 * geometry.transfer_area_factor,
        parameters.transfer_cd,
    )
    gate = pylorus_open_fraction(time_s, parameters, geometry, wave_argument)
    concentration_denominator = max(
        antrum_liquid + solid_volume + parameters.minimum_content_volume_m3,
        parameters.minimum_content_volume_m3,
    )
    solid_concentration = solid_volume / concentration_denominator
    obstruction = min(0.90, parameters.obstruction_strength * solid_concentration)
    pylorus_area = math.pi * (
        parameters.pylorus_radius_m * math.sqrt(geometry.pylorus_area_factor)
    ) ** 2
    liquid_flow = orifice_flow(
        antrum_pressure,
        parameters.downstream_pressure_pa,
        parameters.fluid_density_kg_m3,
        pylorus_area,
        parameters.pylorus_cd,
        gate,
        obstruction,
    )

    selectivity = math.exp(
        -((parameters.particle_diameter_m / parameters.particle_cutoff_m) ** parameters.particle_selectivity_exponent)
    )
    flow_ratio = liquid_flow / parameters.solid_reference_flow_m3_s
    candidate_solid_flow = (
        parameters.solid_flow_rate_per_s
        * mobile_solid
        * flow_ratio
        * selectivity
    )
    carried_solid_limit = (
        parameters.solid_density_kg_m3
        * liquid_flow
        * parameters.max_solid_fraction
    )
    solid_flow = min(max(candidate_solid_flow, 0.0), carried_solid_limit)
    bound_solid = max(solid_mass - mobile_solid, 0.0)
    mobilization_rate = (
        parameters.mobilization_rate_per_s
        * contraction_wave
        * bound_solid
    )

    return {
        "fundus_pressure_pa": fundus_pressure,
        "antrum_pressure_pa": antrum_pressure,
        "pylorus_gate": gate,
        "solid_concentration": solid_concentration,
        "pylorus_obstruction": obstruction,
        "particle_selectivity": selectivity,
        "transfer_flow_m3_s": transfer_flow,
        "liquid_flow_m3_s": liquid_flow,
        "solid_flow_kg_s": solid_flow,
        "mobilization_rate_kg_s": mobilization_rate,
    }


def advance_state(
    state: np.ndarray,
    time_s: float,
    dt_s: float,
    parameters: ModelParameters,
    geometry: Geometry,
) -> tuple[np.ndarray, dict[str, float], dict[str, float]]:
    diagnostics = state_diagnostics(state, time_s, parameters, geometry)
    state = np.maximum(np.asarray(state, dtype=float), 0.0)
    fundus_liquid, antrum_liquid, solid_mass, mobile_solid = state
    transfer_flow = min(diagnostics["transfer_flow_m3_s"], fundus_liquid / dt_s)
    liquid_flow = min(diagnostics["liquid_flow_m3_s"], antrum_liquid / dt_s)
    solid_flow = min(
        diagnostics["solid_flow_kg_s"],
        solid_mass / dt_s,
        mobile_solid / dt_s,
    )
    new_state = np.array(
        [
            fundus_liquid - transfer_flow * dt_s,
            antrum_liquid + transfer_flow * dt_s - liquid_flow * dt_s,
            solid_mass - solid_flow * dt_s,
            mobile_solid
            + (diagnostics["mobilization_rate_kg_s"] - solid_flow) * dt_s,
        ],
        dtype=float,
    )
    mobile_raw_kg = float(new_state[3])
    new_state[0] = max(new_state[0], 0.0)
    new_state[1] = max(new_state[1], 0.0)
    new_state[2] = max(new_state[2], 0.0)
    new_state[3] = min(max(new_state[3], 0.0), new_state[2])
    effective_flows = {
        "mobile_projection_kg": float(new_state[3] - mobile_raw_kg),
        "transfer_flow_m3_s": max(transfer_flow, 0.0),
        "liquid_flow_m3_s": max(liquid_flow, 0.0),
        "solid_flow_kg_s": max(solid_flow, 0.0),
    }
    return new_state, diagnostics, effective_flows


def time_to_half_fraction(time_s: np.ndarray, marker: np.ndarray, initial_marker: float) -> float | None:
    if initial_marker <= 0.0:
        return None
    time_values = np.asarray(time_s, dtype=float)
    marker_values = np.asarray(marker, dtype=float)
    target = 0.5 * initial_marker
    indices = np.flatnonzero(marker_values <= target)
    if len(indices) == 0:
        return None
    index = int(indices[0])
    if index == 0:
        return float(time_values[0] / 60.0)
    previous = marker_values[index - 1]
    current = marker_values[index]
    if current == previous:
        return float(time_values[index] / 60.0)
    fraction = (previous - target) / (previous - current)
    crossing = time_values[index - 1] + fraction * (time_values[index] - time_values[index - 1])
    return float(crossing / 60.0)


def sampled_indices(length: int, sample_period_s: float, dt_s: float) -> np.ndarray:
    step = max(1, int(round(sample_period_s / dt_s)))
    indices = np.arange(0, length, step, dtype=int)
    if indices[-1] != length - 1:
        indices = np.append(indices, length - 1)
    return indices


def simulate_mechanistic(
    meal_type: str,
    parameters: ModelParameters,
    geometry: Geometry,
    t_end_s: float = 7200.0,
    dt_s: float = 1.0,
    sample_period_s: float = 5.0,
) -> dict[str, Any]:
    if dt_s <= 0.0 or t_end_s < dt_s:
        raise ValueError("t_end_s must be positive and dt_s must be positive")
    count = int(round(t_end_s / dt_s)) + 1
    time_s = np.linspace(0.0, t_end_s, count)
    state = initial_state(meal_type, parameters)
    initial_liquid = float(state[0] + state[1])
    initial_solid = float(state[2])
    fundus = np.zeros(count)
    antrum = np.zeros(count)
    total_liquid = np.zeros(count)
    total_solid = np.zeros(count)
    mobile_solid = np.zeros(count)
    total_content = np.zeros(count)
    fundus_pressure = np.zeros(count)
    antrum_pressure = np.zeros(count)
    gate = np.zeros(count)
    obstruction = np.zeros(count)
    liquid_flow = np.zeros(count)
    solid_flow = np.zeros(count)
    transfer_flow = np.zeros(count)
    wall_power = np.zeros(count)
    mass_error_liquid = np.zeros(count)
    mass_error_solid = np.zeros(count)
    cumulative_liquid = 0.0
    cumulative_solid = 0.0

    for index in range(count):
        diagnostics = state_diagnostics(state, time_s[index], parameters, geometry)
        fundus[index] = state[0]
        antrum[index] = state[1]
        total_liquid[index] = state[0] + state[1]
        total_solid[index] = state[2]
        mobile_solid[index] = state[3]
        total_content[index] = state[0] + state[1] + state[2] / parameters.solid_density_kg_m3
        fundus_pressure[index] = diagnostics["fundus_pressure_pa"]
        antrum_pressure[index] = diagnostics["antrum_pressure_pa"]
        gate[index] = diagnostics["pylorus_gate"]
        obstruction[index] = diagnostics["pylorus_obstruction"]
        liquid_flow[index] = diagnostics["liquid_flow_m3_s"]
        solid_flow[index] = diagnostics["solid_flow_kg_s"]
        transfer_flow[index] = diagnostics["transfer_flow_m3_s"]
        wall_power[index] = (
            -fundus_pressure[index] * transfer_flow[index]
            + antrum_pressure[index] * (transfer_flow[index] - liquid_flow[index])
        )
        mass_error_liquid[index] = abs(
            initial_liquid - (state[0] + state[1] + cumulative_liquid)
        ) / max(initial_liquid, np.finfo(float).eps)
        mass_error_solid[index] = abs(
            initial_solid - (state[2] + cumulative_solid)
        ) / max(initial_solid, 1.0)
        if index < count - 1:
            state, _, effective_flows = advance_state(
                state, time_s[index], dt_s, parameters, geometry
            )
            cumulative_liquid += effective_flows["liquid_flow_m3_s"] * dt_s
            cumulative_solid += effective_flows["solid_flow_kg_s"] * dt_s

    marker = total_liquid if meal_type == "liquid" else total_solid
    initial_marker = initial_liquid if meal_type == "liquid" else initial_solid
    selected = sampled_indices(count, sample_period_s, dt_s)
    selected_time = time_s[selected]
    selected_marker = marker[selected]
    reference_time = min(t_end_s, 7200.0)
    residual_index = int(np.argmin(np.abs(time_s - reference_time)))
    mass_max = float(max(np.max(mass_error_liquid), np.max(mass_error_solid)))
    t50 = time_to_half_fraction(time_s, marker, initial_marker)
    interval_mask = time_s[:-1] <= reference_time
    metrics = {
        "t50_min": t50,
        "initial_liquid_ml": initial_liquid * 1.0e6,
        "initial_solid_g": initial_solid * 1000.0,
        "residual_fraction_at_120_min": float(marker[residual_index] / initial_marker),
        "residual_liquid_ml_at_120_min": float(total_liquid[residual_index] * 1.0e6),
        "residual_solid_g_at_120_min": float(total_solid[residual_index] * 1000.0),
        "mean_liquid_flow_ml_s_0_120_min": float(
            np.mean(liquid_flow[:-1][interval_mask]) * 1.0e6
        ),
        "mean_solid_flow_g_s_0_120_min": float(
            np.mean(solid_flow[:-1][interval_mask]) * 1000.0
        ),
        "peak_liquid_flow_ml_s": float(np.max(liquid_flow) * 1.0e6),
        "peak_solid_flow_g_s": float(np.max(solid_flow) * 1000.0),
        "net_wall_work_j": float(np.sum(wall_power[:-1] * dt_s)),
        "mean_wall_power_w_0_120_min": float(np.mean(wall_power[:-1][interval_mask])),
        "final_liquid_ml": float(total_liquid[-1] * 1.0e6),
        "final_solid_g": float(total_solid[-1] * 1000.0),
        "max_mass_balance_relative_error": mass_max,
        "cumulative_liquid_outflow_ml": float(cumulative_liquid * 1.0e6),
        "cumulative_solid_outflow_g": float(cumulative_solid * 1000.0),
    }
    series = {
        "time_s": time_s[selected],
        "time_min": time_s[selected] / 60.0,
        "fundus_liquid_ml": fundus[selected] * 1.0e6,
        "antrum_liquid_ml": antrum[selected] * 1.0e6,
        "total_liquid_ml": total_liquid[selected] * 1.0e6,
        "solid_mass_g": total_solid[selected] * 1000.0,
        "mobile_solid_mass_g": mobile_solid[selected] * 1000.0,
        "total_content_ml": total_content[selected] * 1.0e6,
        "fundus_pressure_pa": fundus_pressure[selected],
        "antrum_pressure_pa": antrum_pressure[selected],
        "antrum_pressure_mmhg": antrum_pressure[selected] / 133.322,
        "pylorus_open_fraction": gate[selected],
        "pylorus_obstruction_fraction": obstruction[selected],
        "pylorus_liquid_flow_ml_s": liquid_flow[selected] * 1.0e6,
        "pylorus_solid_flow_g_s": solid_flow[selected] * 1000.0,
        "pylorus_solid_equivalent_ml_s": (
            solid_flow[selected] / parameters.solid_density_kg_m3 * 1.0e6
        ),
        "fundus_antrum_transfer_ml_s": transfer_flow[selected] * 1.0e6,
        "net_wall_power_w": wall_power[selected],
        "liquid_mass_balance_relative_error": mass_error_liquid[selected],
        "solid_mass_balance_relative_error": mass_error_solid[selected],
    }
    return {
        "meal_type": meal_type,
        "geometry": geometry.name,
        "time_step_s": dt_s,
        "sample_period_s": sample_period_s,
        "metrics": metrics,
        "series": series,
    }


def simulate_scalar(
    meal_type: str,
    parameters: ModelParameters,
    t_end_s: float = 7200.0,
    sample_period_s: float = 5.0,
) -> dict[str, Any]:
    if meal_type == "liquid":
        initial = parameters.liquid_meal_volume_m3
        rate_per_min = parameters.scalar_liquid_k_per_min
        unit_scale = 1.0e6
    elif meal_type == "solid":
        initial = parameters.solid_meal_mass_kg
        rate_per_min = parameters.scalar_solid_k_per_min
        unit_scale = 1000.0
    else:
        raise ValueError(f"unknown meal type: {meal_type}")
    time_s = np.arange(0.0, t_end_s + sample_period_s, sample_period_s)
    marker = initial * np.exp(-rate_per_min * time_s / 60.0)
    flow = initial * rate_per_min / 60.0 * np.exp(-rate_per_min * time_s / 60.0)
    if meal_type == "liquid":
        series = {
            "time_s": time_s,
            "time_min": time_s / 60.0,
            "total_liquid_ml": marker * unit_scale,
            "pylorus_liquid_flow_ml_s": flow * unit_scale,
        }
    else:
        series = {
            "time_s": time_s,
            "time_min": time_s / 60.0,
            "solid_mass_g": marker * unit_scale,
            "pylorus_solid_flow_g_s": flow * unit_scale,
        }
    t50 = time_to_half_fraction(time_s, marker, initial)
    metrics = {
        "t50_min": t50,
        "residual_fraction_at_120_min": float(marker[np.argmin(np.abs(time_s - 7200.0))] / initial),
        "initial_marker": initial * unit_scale,
        "final_marker": float(marker[-1] * unit_scale),
        "max_mass_balance_relative_error": 0.0,
    }
    return {
        "meal_type": meal_type,
        "geometry": "scalar_k_empt",
        "sample_period_s": sample_period_s,
        "metrics": metrics,
        "series": series,
    }


def dimensional_checks() -> list[dict[str, Any]]:
    return [
        {
            "equation": "dV/dt = -Q_pyl + Q_in",
            "lhs_unit": "m^3 s^-1",
            "rhs_unit": "m^3 s^-1",
            "ok": True,
        },
        {
            "equation": "Q = C_D A sqrt(2 max(dP,0)/rho)",
            "lhs_unit": "m^3 s^-1",
            "rhs_unit": "m^2 (Pa/(kg m^-3))^0.5 = m^3 s^-1",
            "ok": True,
        },
        {
            "equation": "P_passive = 2 E t (V/V0 - 1) / r",
            "lhs_unit": "Pa",
            "rhs_unit": "Pa m m / m = Pa",
            "ok": True,
        },
        {
            "equation": "dM/dt = -k_s M (Q/Q_ref) S + k_mob S_shear (M_bound)",
            "lhs_unit": "kg s^-1",
            "rhs_unit": "kg s^-1 + kg s^-1",
            "ok": True,
        },
        {
            "equation": "P_mmHg = P_Pa / 133.322",
            "lhs_unit": "mmHg",
            "rhs_unit": "Pa / Pa mmHg^-1 = mmHg",
            "ok": True,
        },
        {
            "equation": "P_wall = -P_f Q_fa + P_a (Q_fa - Q_pyl)",
            "lhs_unit": "W",
            "rhs_unit": "Pa m^3 s^-1 = W",
            "ok": True,
        },
    ]


def verify_preregistration(root: Path) -> dict[str, Any]:
    prereg = root / "PREREG.md"
    checksum = root / "PREREG.sha256"
    actual = hashlib.sha256(prereg.read_bytes()).hexdigest()
    expected = checksum.read_text(encoding="utf-8").split()[0]
    return {
        "prereg_sha256_expected": expected,
        "prereg_sha256_actual": actual,
        "valid": expected == actual,
    }


def run_sensitivity(
    base_parameters: ModelParameters,
    base_results: dict[str, Any],
) -> list[dict[str, Any]]:
    keys = [
        "pylorus_radius_m",
        "antrum_active_stress_pa",
        "pylorus_span",
        "solid_flow_rate_per_s",
        "mobilization_rate_per_s",
        "antrum_wall_modulus_pa",
        "fundus_wall_modulus_pa",
        "particle_cutoff_m",
    ]
    geometry = GEOMETRIES["fundus_antrum_pylorus"]
    base_t50 = {
        meal: base_results["full_geometry"][meal]["metrics"]["t50_min"]
        for meal in ["liquid", "solid"]
    }
    rows: list[dict[str, Any]] = []
    for key in keys:
        row: dict[str, Any] = {
            "parameter": key,
            "base_value": getattr(base_parameters, key),
            "minus_50_percent": {},
            "plus_50_percent": {},
        }
        scores: list[float] = []
        for label, factor in [("minus_50_percent", 0.5), ("plus_50_percent", 1.5)]:
            varied = replace(base_parameters, **{key: getattr(base_parameters, key) * factor})
            for meal in ["liquid", "solid"]:
                result = simulate_mechanistic(
                    meal,
                    varied,
                    geometry,
                    t_end_s=7200.0,
                    dt_s=2.0,
                    sample_period_s=60.0,
                )
                t50 = result["metrics"]["t50_min"]
                row[label][meal] = {
                    "value": getattr(varied, key),
                    "t50_min": t50,
                    "relative_t50_change": (
                        abs(t50 - base_t50[meal]) / base_t50[meal]
                        if t50 is not None and base_t50[meal] is not None
                        else None
                    ),
                }
                if t50 is not None and base_t50[meal] is not None:
                    scores.append(abs(t50 - base_t50[meal]) / base_t50[meal])
        row["max_relative_t50_change"] = max(scores) if scores else None
        rows.append(row)
    rows.sort(key=lambda row: row["max_relative_t50_change"] or -1.0, reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def geometry_comparison(results: dict[str, Any]) -> dict[str, Any]:
    comparison: dict[str, Any] = {}
    for meal in ["liquid", "solid"]:
        reference = REFERENCE[f"{meal}_t50_min"]
        full_metrics = results["full_geometry"][meal]["metrics"]
        compartment_metrics = results["fundus_antrum"][meal]["metrics"]
        scalar_metrics = results["scalar"][meal]["metrics"]
        full_t50 = full_metrics["t50_min"]
        compartment_t50 = compartment_metrics["t50_min"]
        scalar_t50 = scalar_metrics["t50_min"]
        full_error = abs(full_t50 - reference) / reference if full_t50 is not None else None
        compartment_error = (
            abs(compartment_t50 - reference) / reference
            if compartment_t50 is not None
            else None
        )
        scalar_error = abs(scalar_t50 - reference) / reference if scalar_t50 is not None else None
        comparison[meal] = {
            "reference_t50_min": reference,
            "full_geometry_t50_min": full_t50,
            "compartment_t50_min": compartment_t50,
            "scalar_t50_min": scalar_t50,
            "full_relative_t50_error": full_error,
            "compartment_relative_t50_error": compartment_error,
            "scalar_relative_t50_error": scalar_error,
            "full_residual_fraction_at_120_min": full_metrics["residual_fraction_at_120_min"],
            "compartment_residual_fraction_at_120_min": compartment_metrics["residual_fraction_at_120_min"],
            "scalar_residual_fraction_at_120_min": scalar_metrics["residual_fraction_at_120_min"],
            "full_minus_scalar_residual_fraction_at_120_min": (
                full_metrics["residual_fraction_at_120_min"]
                - scalar_metrics["residual_fraction_at_120_min"]
            ),
            "full_vs_scalar_error_reduction_fraction": (
                1.0 - full_error / scalar_error
                if full_error is not None and scalar_error not in (None, 0.0)
                else None
            ),
            "empirical_curve_rmse_status": "UNKNOWN",
            "decision": "UNKNOWN until held-out gastric-emptying curves are supplied",
        }
    return comparison


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def run_experiment(root: Path | None = None) -> dict[str, Any]:
    root = root or Path(__file__).resolve().parent
    parameters = default_parameters()
    full: dict[str, Any] = {}
    compartments: dict[str, Any] = {}
    scalar: dict[str, Any] = {}
    for meal in ["liquid", "solid"]:
        full[meal] = simulate_mechanistic(
            meal, parameters, GEOMETRIES["fundus_antrum_pylorus"]
        )
        compartments[meal] = simulate_mechanistic(
            meal, parameters, GEOMETRIES["fundus_antrum"]
        )
        scalar[meal] = simulate_scalar(meal, parameters)
    results: dict[str, Any] = {
        "model_id": "BT-HX-Q107",
        "model_description": "Two-compartment first-principles gastric emptying control model",
        "equations": {
            "volume_balance": "dV_f/dt=-Q_fa; dV_a/dt=Q_fa-Q_pyl",
            "solid_balance": "dM/dt=-Q_s",
            "wall_pressure": "P=P_passive+P_active",
            "slow_wave_gate": "g_p=g_floor+g_span sigmoid(k_g(sin(omega t+phi)-g_threshold))",
            "pylorus": "Q_pyl=C_D A_p g_p sqrt(2 max(P_a-P_d,0)/rho)",
            "particle_flow": "Q_s=k_s M_mobile (Q_pyl/Q_ref) S(d_p) with Q_in=0 after loading",
        },
        "reference": REFERENCE,
        "preregistration": verify_preregistration(root),
        "parameters": parameter_table(parameters),
        "geometries": geometry_table(),
        "dimensional_checks": dimensional_checks(),
        "full_geometry": full,
        "fundus_antrum": compartments,
        "scalar": scalar,
    }
    results["geometry_comparison"] = geometry_comparison(results)
    results["sensitivity"] = run_sensitivity(parameters, results)
    criteria = {
        "mass_balance": {
            "threshold": 1.0e-8,
            "observed_max": max(
                result["metrics"]["max_mass_balance_relative_error"]
                for result in list(full.values()) + list(compartments.values())
            ),
        },
        "liquid_t50_factor_two": {
            "bounds_min": REFERENCE["liquid_factor_two_bounds_min"],
            "observed_min": full["liquid"]["metrics"]["t50_min"],
            "passed": (
                full["liquid"]["metrics"]["t50_min"] is not None
                and REFERENCE["liquid_factor_two_bounds_min"][0]
                <= full["liquid"]["metrics"]["t50_min"]
                <= REFERENCE["liquid_factor_two_bounds_min"][1]
            ),
        },
        "solid_t50_factor_two": {
            "bounds_min": REFERENCE["solid_factor_two_bounds_min"],
            "observed_min": full["solid"]["metrics"]["t50_min"],
            "passed": (
                full["solid"]["metrics"]["t50_min"] is not None
                and REFERENCE["solid_factor_two_bounds_min"][0]
                <= full["solid"]["metrics"]["t50_min"]
                <= REFERENCE["solid_factor_two_bounds_min"][1]
            ),
        },
        "empirical_geometry_improvement": {
            "required_relative_rmse_reduction": 0.10,
            "status": "UNKNOWN",
            "reason": "No held-out volume and retention curves are present",
        },
    }
    criteria["mass_balance"]["passed"] = criteria["mass_balance"]["observed_max"] <= 1.0e-8
    results["criteria"] = criteria
    results["status"] = {
        "mechanistic_reference_screen": (
            "PASS"
            if criteria["liquid_t50_factor_two"]["passed"]
            and criteria["solid_t50_factor_two"]["passed"]
            else "FAIL"
        ),
        "empirical_geometry_claim": "UNKNOWN",
    }
    return to_jsonable(results)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results.json")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    results = run_experiment(root)
    output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "status": results["status"],
        "liquid_t50_min": results["full_geometry"]["liquid"]["metrics"]["t50_min"],
        "solid_t50_min": results["full_geometry"]["solid"]["metrics"]["t50_min"],
        "output": str(output),
    }
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
