"""Reduced first-principles dynamic model for Q019."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

SEGMENT_NAMES = ("PCT", "DL", "TAL", "DCT", "CNT", "CCD", "IMCD")
SOLUTE_NAMES = ("na", "cl", "urea")
SOLUTE_LABELS = {"na": "NaCl-equivalent", "cl": "Cl-equivalent", "urea": "urea"}
OSMOLALITY_MMHG_PER_MM = 19.37
SECONDS_PER_MINUTE = 60.0
EPSILON = 1.0e-18


@dataclass(frozen=True)
class Parameters:
    sngfr_ml_min: float = 1.0e-4
    nephron_count: int = 1_000_000
    kf_nl_min_mmhg: float = 6.666666666666667
    p_gc_mmhg: float = 50.0
    p_bs_mmhg: float = 10.0
    pi_gc_mmhg: float = 25.0
    pi_bs_mmhg: float = 0.0
    temperature_c: float = 37.0
    fixed_luminal_osm_mM: float = 40.0
    avp_level: float = 1.0
    avp_floor: float = 0.10
    duration_min: float = 60.0
    output_window_min: float = 15.0
    rtol: float = 1.0e-7
    atol: float = 1.0e-11
    max_step_min: float = 0.5
    lengths_cm: tuple[float, ...] = (1.70, 0.32, 1.00, 0.20, 0.40, 0.40, 1.70)
    diameters_um: tuple[float, ...] = (37.0, 26.0, 26.0, 20.0, 24.0, 45.0, 50.0)
    residence_min: tuple[float, ...] = (0.85, 0.21, 1.25, 0.40, 0.80, 1.33, 5.10)
    # Legacy field name; values are synthetic hydraulic Lp [cm/(s mmHg)].
    lp_cm_s: tuple[float, ...] = (2.8e-6, 1.0e-6, 0.0, 0.0, 2.0e-9, 2.0e-9, 0.5e-9)
    vmax_na_umol_min: tuple[float, ...] = (0.01131, 0.0, 0.00151, 0.00030, 0.00020, 0.00010, 0.00005)
    vmax_cl_umol_min: tuple[float, ...] = (0.00761, 0.0, 0.00100, 0.00020, 0.00015, 0.00008, 0.00004)
    vmax_urea_umol_min: tuple[float, ...] = (0.00050, 0.0, 0.0, 0.0, 0.0, 0.00005, 0.00010)
    km_na_mM: tuple[float, ...] = (40.0, 30.0, 30.0, 20.0, 20.0, 20.0, 20.0)
    km_cl_mM: tuple[float, ...] = (30.0, 30.0, 30.0, 20.0, 20.0, 20.0, 20.0)
    km_urea_mM: tuple[float, ...] = (5.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0)
    cinter_osm_mM: tuple[float, ...] = (285.9, 480.0, 688.5, 688.5, 688.5, 688.5, 737.9)
    cinter_na_mM: tuple[float, ...] = (140.0, 299.0, 299.0, 299.0, 299.0, 299.0, 239.0)
    cinter_cl_mM: tuple[float, ...] = (100.0, 280.4, 280.4, 280.4, 280.4, 280.4, 236.5)
    cinter_urea_mM: tuple[float, ...] = (5.0, 60.0, 60.0, 60.0, 60.0, 60.0, 200.0)
    filtered_na_mM: float = 140.0
    filtered_cl_mM: float = 100.0
    filtered_urea_mM: float = 5.0

    def scaled_lp(self, segment: str, factor: float) -> "Parameters":
        if segment == "CD":
            indices = (SEGMENT_NAMES.index("CCD"), SEGMENT_NAMES.index("IMCD"))
        elif segment in SEGMENT_NAMES:
            indices = (SEGMENT_NAMES.index(segment),)
        else:
            raise ValueError(f"unknown segment: {segment}")
        values = list(self.lp_cm_s)
        for index in indices:
            values[index] = values[index] * factor
        return replace(self, lp_cm_s=tuple(values))

    def scaled_vmax(self, segment: str, factor: float) -> "Parameters":
        if segment not in SEGMENT_NAMES:
            raise ValueError(f"unknown segment: {segment}")
        values = list(self.vmax_na_umol_min)
        cl_values = list(self.vmax_cl_umol_min)
        urea_values = list(self.vmax_urea_umol_min)
        index = SEGMENT_NAMES.index(segment)
        values[index] = values[index] * factor
        cl_values[index] = cl_values[index] * factor
        urea_values[index] = urea_values[index] * factor
        return replace(
            self,
            vmax_na_umol_min=tuple(values),
            vmax_cl_umol_min=tuple(cl_values),
            vmax_urea_umol_min=tuple(urea_values),
        )


DEFAULT_PARAMETERS = Parameters()


def filtration_rate(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    net_pressure_mmhg = (
        parameters.p_gc_mmhg
        - parameters.p_bs_mmhg
        - (parameters.pi_gc_mmhg - parameters.pi_bs_mmhg)
    )
    return parameters.kf_nl_min_mmhg * net_pressure_mmhg * 1.0e-6


def geometry(parameters: Parameters = DEFAULT_PARAMETERS) -> tuple[np.ndarray, np.ndarray]:
    lengths_cm = np.asarray(parameters.lengths_cm, dtype=float)
    diameters_cm = np.asarray(parameters.diameters_um, dtype=float) * 1.0e-4
    surface_area_cm2 = np.pi * diameters_cm * lengths_cm
    return surface_area_cm2, diameters_cm


def water_conductance(parameters: Parameters = DEFAULT_PARAMETERS) -> np.ndarray:
    surface_area_cm2, _ = geometry(parameters)
    avp_multiplier = parameters.avp_floor + (1.0 - parameters.avp_floor) * parameters.avp_level
    multipliers = np.ones(len(SEGMENT_NAMES), dtype=float)
    multipliers[SEGMENT_NAMES.index("CCD")] = avp_multiplier
    multipliers[SEGMENT_NAMES.index("IMCD")] = avp_multiplier
    return (
        np.asarray(parameters.lp_cm_s, dtype=float)
        * multipliers
        * surface_area_cm2
        * SECONDS_PER_MINUTE
        * OSMOLALITY_MMHG_PER_MM
    )


def solute_arrays(parameters: Parameters) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vmax = np.column_stack(
        (
            np.asarray(parameters.vmax_na_umol_min, dtype=float),
            np.asarray(parameters.vmax_cl_umol_min, dtype=float),
            np.asarray(parameters.vmax_urea_umol_min, dtype=float),
        )
    )
    km = np.column_stack(
        (
            np.asarray(parameters.km_na_mM, dtype=float),
            np.asarray(parameters.km_cl_mM, dtype=float),
            np.asarray(parameters.km_urea_mM, dtype=float),
        )
    )
    cinter = np.column_stack(
        (
            np.asarray(parameters.cinter_na_mM, dtype=float),
            np.asarray(parameters.cinter_cl_mM, dtype=float),
            np.asarray(parameters.cinter_urea_mM, dtype=float),
        )
    )
    return vmax, km, cinter


def filtered_concentrations(parameters: Parameters) -> np.ndarray:
    return np.asarray(
        (parameters.filtered_na_mM, parameters.filtered_cl_mM, parameters.filtered_urea_mM),
        dtype=float,
    )


def initial_state(parameters: Parameters) -> np.ndarray:
    q_filtration = filtration_rate(parameters)
    concentrations = filtered_concentrations(parameters)
    volume = q_filtration * np.asarray(parameters.residence_min, dtype=float)
    amounts = volume[:, None] * concentrations[None, :]
    cumulative_transport = np.zeros(2, dtype=float)
    return np.concatenate((volume, amounts.reshape(-1), cumulative_transport))


def unpack_state(state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_segments = len(SEGMENT_NAMES)
    n_solutes = len(SOLUTE_NAMES)
    volumes = state[:n_segments]
    amounts = state[n_segments:n_segments + n_segments * n_solutes].reshape(n_segments, n_solutes)
    return volumes, amounts


def state_to_frame(state: np.ndarray, time_min: float, parameters: Parameters) -> dict[str, Any]:
    volumes, amounts = unpack_state(state)
    concentrations = amounts / np.maximum(volumes[:, None], EPSILON)
    terminal_volume = max(float(volumes[-1]), 0.0)
    terminal_amount = np.maximum(amounts[-1], 0.0)
    terminal_concentration = terminal_amount / max(terminal_volume, EPSILON)
    terminal_osm = float(np.sum(terminal_concentration) + parameters.fixed_luminal_osm_mM)
    return {
        "time_min": float(time_min),
        "segment_volume_ml_per_nephron": volumes.tolist(),
        "segment_concentration_mM": concentrations.tolist(),
        "terminal_flow_ml_min_per_kidney": terminal_volume
        / max(parameters.residence_min[-1], EPSILON)
        * parameters.nephron_count,
        "terminal_na_umol_min_per_kidney": terminal_amount[0]
        / max(parameters.residence_min[-1], EPSILON)
        * parameters.nephron_count,
        "terminal_cl_umol_min_per_kidney": terminal_amount[1]
        / max(parameters.residence_min[-1], EPSILON)
        * parameters.nephron_count,
        "terminal_urea_umol_min_per_kidney": terminal_amount[2]
        / max(parameters.residence_min[-1], EPSILON)
        * parameters.nephron_count,
        "terminal_osm_mM": terminal_osm,
    }


def _safe_flux(flux: float, amount: float, fraction_limit: float = 0.98) -> float:
    if amount <= EPSILON:
        return 0.0
    return float(np.clip(flux, -fraction_limit * amount, fraction_limit * amount))


def rhs(time_min: float, state: np.ndarray, parameters: Parameters) -> np.ndarray:
    del time_min
    volumes, amounts = unpack_state(state)
    volumes_safe = np.maximum(volumes, EPSILON)
    concentrations = amounts / volumes_safe[:, None]
    concentrations = np.maximum(concentrations, 0.0)
    q_filtration = filtration_rate(parameters)
    residence = np.asarray(parameters.residence_min, dtype=float)
    q_out = volumes_safe / residence
    q_in = np.empty(len(SEGMENT_NAMES), dtype=float)
    q_in[0] = q_filtration
    q_in[1:] = q_out[:-1]
    incoming_concentration = np.empty_like(concentrations)
    incoming_concentration[0] = filtered_concentrations(parameters)
    incoming_concentration[1:] = concentrations[:-1]
    vmax, km, _ = solute_arrays(parameters)
    active_flux = vmax * concentrations / (km + concentrations)
    active_flux = np.minimum(active_flux, amounts * 4.0)
    interstitium_osm = np.asarray(parameters.cinter_osm_mM, dtype=float)
    luminal_osm = np.sum(concentrations, axis=1) + parameters.fixed_luminal_osm_mM
    water_flux = water_conductance(parameters) * (interstitium_osm - luminal_osm)
    water_flux = np.minimum(np.maximum(water_flux, 0.0), q_in * 0.98)
    derivative_volume = q_in - q_out - water_flux
    derivative_amount = q_in[:, None] * incoming_concentration - amounts / residence[:, None] - active_flux
    # Donor-limited active efflux is above; do not limit net incoming mass.
    cumulative_water_reabsorbed = np.asarray([np.sum(water_flux)], dtype=float)
    cumulative_terminal_excretion = np.asarray([q_out[-1]], dtype=float)
    derivative_state = np.concatenate(
        (derivative_volume, derivative_amount.reshape(-1), cumulative_water_reabsorbed, cumulative_terminal_excretion)
    )
    return derivative_state


def no_transport_state(parameters: Parameters) -> Parameters:
    zero_lp = tuple(0.0 for _ in SEGMENT_NAMES)
    zero_vmax = tuple(0.0 for _ in SEGMENT_NAMES)
    return replace(
        parameters,
        lp_cm_s=zero_lp,
        vmax_na_umol_min=zero_vmax,
        vmax_cl_umol_min=zero_vmax,
        vmax_urea_umol_min=zero_vmax,
    )


def simulate(
    parameters: Parameters = DEFAULT_PARAMETERS,
    time_grid: np.ndarray | None = None,
) -> dict[str, Any]:
    if time_grid is None:
        time_grid = np.linspace(0.0, parameters.duration_min, 481)
    else:
        time_grid = np.asarray(time_grid, dtype=float)
    state0 = initial_state(parameters)
    solution = solve_ivp(
        rhs,
        (float(time_grid[0]), float(time_grid[-1])),
        state0,
        args=(parameters,),
        t_eval=time_grid,
        method="RK45",
        rtol=parameters.rtol,
        atol=parameters.atol,
        max_step=parameters.max_step_min,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    frames = [state_to_frame(state, time_min, parameters) for state, time_min in zip(solution.y.T, solution.t)]
    final_time = float(solution.t[-1])
    window_start = max(float(solution.t[0]), final_time - parameters.output_window_min)
    window_indices = [index for index, time_min in enumerate(solution.t) if time_min >= window_start]
    final_window = [frames[index] for index in window_indices]
    terminal_flow = float(np.mean([frame["terminal_flow_ml_min_per_kidney"] for frame in final_window]))
    terminal_na = float(np.mean([frame["terminal_na_umol_min_per_kidney"] for frame in final_window]))
    terminal_cl = float(np.mean([frame["terminal_cl_umol_min_per_kidney"] for frame in final_window]))
    terminal_urea = float(np.mean([frame["terminal_urea_umol_min_per_kidney"] for frame in final_window]))
    terminal_osm = float(np.mean([frame["terminal_osm_mM"] for frame in final_window]))
    mass_residual = mass_balance_residual(solution.y, solution.t, parameters)
    return {
        "time_min": solution.t.tolist(),
        "frames": frames,
        "terminal_output": {
            "urine_flow_ml_min_per_kidney": terminal_flow,
            "urinary_na_equiv_umol_min_per_kidney": terminal_na,
            "urinary_cl_equiv_umol_min_per_kidney": terminal_cl,
            "urinary_urea_umol_min_per_kidney": terminal_urea,
            "urine_osm_mM": terminal_osm,
        },
        "mass_balance": mass_residual,
        "solver": {
            "success": bool(solution.success),
            "message": solution.message,
            "n_steps": int(solution.t.size),
            "rtol": parameters.rtol,
            "atol": parameters.atol,
            "max_step_min": parameters.max_step_min,
        },
    }


def mass_balance_residual(
    states: np.ndarray,
    times: np.ndarray,
    parameters: Parameters,
    dense_solution: Any | None = None,
) -> dict[str, float]:
    del dense_solution
    n_segments = len(SEGMENT_NAMES)
    n_solutes = len(SOLUTE_NAMES)
    initial_volumes = initial_state(parameters)[:n_segments]
    final_volumes = states[:n_segments, -1]
    residence = np.asarray(parameters.residence_min, dtype=float)
    q_filtration = filtration_rate(parameters)
    filtered_water = q_filtration * (times[-1] - times[0])
    cumulative_start = n_segments + n_segments * n_solutes
    if states.shape[0] >= cumulative_start + 2:
        reabsorbed_water = float(states[cumulative_start, -1])
        excreted_water = float(states[cumulative_start + 1, -1])
    else:
        reabsorbed_water = 0.0
        excreted_water = float(np.trapezoid(states[n_segments - 1, :] / residence[-1], times))
    initial_water = float(np.sum(initial_volumes))
    final_water = float(np.sum(final_volumes))
    water_residual = initial_water + filtered_water - reabsorbed_water - excreted_water - final_water
    water_scale = max(abs(filtered_water), abs(reabsorbed_water), abs(excreted_water), abs(initial_water), 1.0e-12)
    return {
        "initial_stored_water_ml_per_nephron": initial_water,
        "final_stored_water_ml_per_nephron": final_water,
        "filtered_water_ml_per_nephron": filtered_water,
        "reabsorbed_water_ml_per_nephron": reabsorbed_water,
        "excreted_water_ml_per_nephron": excreted_water,
        "terminal_flow_ml_min_per_nephron_at_end": float(final_volumes[-1] / residence[-1]),
        "relative_total_water_balance_residual": water_residual / water_scale,
    }


def parameter_table(parameters: Parameters = DEFAULT_PARAMETERS) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = [
        {
            "name": "SNGFR",
            "value": parameters.sngfr_ml_min,
            "unit": "mL/min/nephron",
            "source_or_assumption": "Layton & Layton 2019, model-parameter text: 100 nL/min",
        },
        {
            "name": "nephron_count",
            "value": parameters.nephron_count,
            "unit": "nephrons/kidney",
            "source_or_assumption": "Layton & Layton 2019, approximately 10^6 nephrons",
        },
        {
            "name": "Kf",
            "value": parameters.kf_nl_min_mmhg,
            "unit": "nL/(min mmHg)",
            "source_or_assumption": "derived from 100 nL/min and the explicit net-pressure assumption",
        },
        {
            "name": "PCT length",
            "value": parameters.lengths_cm[0],
            "unit": "cm",
            "source_or_assumption": "Layton & Layton 2019, Table 1",
        },
        {
            "name": "PCT diameter",
            "value": parameters.diameters_um[0],
            "unit": "micrometre",
            "source_or_assumption": "Layton & Layton 2019, Table 1",
        },
        {
            "name": "PCT Lp",
            "value": parameters.lp_cm_s[0],
            "unit": "cm/(s mmHg)",
            "source_or_assumption": "assumption selected for the reduced first-principles baseline",
        },
        {
            "name": "TAL vmax",
            "value": parameters.vmax_na_umol_min[SEGMENT_NAMES.index("TAL")],
            "unit": "micromol/min/nephron",
            "source_or_assumption": "assumption for segment-specific active NaCl reabsorption",
        },
        {
            "name": "CD Lp",
            "value": parameters.lp_cm_s[SEGMENT_NAMES.index("IMCD")],
            "unit": "cm/(s mmHg)",
            "source_or_assumption": "assumption for AVP-regulated collecting-system water permeability",
        },
        {
            "name": "interstitial osmolality cortex",
            "value": parameters.cinter_osm_mM[0],
            "unit": "mosm/(kg H2O)",
            "source_or_assumption": "Layton & Layton 2019, Table 2",
        },
        {
            "name": "interstitial osmolality papillary tip",
            "value": parameters.cinter_osm_mM[-1],
            "unit": "mosm/(kg H2O)",
            "source_or_assumption": "Layton & Layton 2019, Table 2",
        },
        {
            "name": "reference urine flow",
            "value": 0.62,
            "unit": "ml/min/kidney",
            "source_or_assumption": "Layton & Layton 2019, Table 3, baseline; reporting anchor only",
        },
    ]
    for index, name in enumerate(SEGMENT_NAMES):
        rows.append(
            {
                "name": f"{name} residence time",
                "value": parameters.residence_min[index],
                "unit": "min",
                "source_or_assumption": "axial-transit assumption",
            }
        )
    return rows


def unit_check(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    area, diameter = geometry(parameters)
    conductance = water_conductance(parameters)
    expected_conversion = area * SECONDS_PER_MINUTE * OSMOLALITY_MMHG_PER_MM
    checks = {
        "filtration_rate_ml_min_per_nephron": filtration_rate(parameters),
        "filtration_rate_target_ml_min_per_nephron": parameters.sngfr_ml_min,
        "area_cm2": area.tolist(),
        "diameter_cm": diameter.tolist(),
        "water_conductance_ml_min_per_mM": conductance.tolist(),
        "conversion_factor_ml_min_per_mM_per_Lp": expected_conversion.tolist(),
        "solute_flux_unit": "micromol/min",
        "concentration_unit": "mM = micromol/mL",
        "water_flux_sign": "positive means lumen to interstitium",
    }
    checks["status"] = "PASS" if abs(filtration_rate(parameters) - parameters.sngfr_ml_min) < 1.0e-12 else "FAIL"
    return checks


def zero_transport_limit(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    result = simulate(no_transport_state(parameters), np.linspace(0.0, parameters.duration_min, 241))
    target = parameters.sngfr_ml_min * parameters.nephron_count
    measured = result["terminal_output"]["urine_flow_ml_min_per_kidney"]
    return {
        "target_flow_ml_min_per_kidney": target,
        "measured_flow_ml_min_per_kidney": measured,
        "relative_error": abs(measured - target) / max(abs(target), 1.0e-12),
        "water_mass_residual": result["mass_balance"]["relative_total_water_balance_residual"],
    }


def perturbation_label(segment: str, quantity: str) -> str:
    if quantity == "lp":
        return f"{segment}_Lp"
    if quantity == "vmax":
        return f"{segment}_vmax"
    raise ValueError(quantity)


def run_sensitivity(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    baseline = simulate(parameters)
    base_output = baseline["terminal_output"]
    perturbations = (
        ("PCT", "lp"),
        ("TAL", "vmax"),
        ("CD", "lp"),
    )
    results: list[dict[str, Any]] = []
    for segment, quantity in perturbations:
        if quantity == "lp":
            plus_parameters = parameters.scaled_lp(segment, 1.5)
            minus_parameters = parameters.scaled_lp(segment, 0.5)
        else:
            plus_parameters = parameters.scaled_vmax(segment, 1.5)
            minus_parameters = parameters.scaled_vmax(segment, 0.5)
        plus = simulate(plus_parameters)
        minus = simulate(minus_parameters)
        plus_output = plus["terminal_output"]
        minus_output = minus["terminal_output"]
        sensitivities: dict[str, float] = {}
        for output_name, baseline_value in base_output.items():
            if abs(baseline_value) > 1.0e-15:
                sensitivities[output_name] = (plus_output[output_name] - minus_output[output_name]) / (
                    2.0 * baseline_value
                )
            else:
                sensitivities[output_name] = 0.0
        results.append(
            {
                "parameter": perturbation_label(segment, quantity),
                "segment": segment,
                "quantity": quantity,
                "plus_50_percent": plus_output,
                "minus_50_percent": minus_output,
                "sensitivity": sensitivities,
            }
        )
    flow_ranking = sorted(results, key=lambda item: abs(item["sensitivity"]["urine_flow_ml_min_per_kidney"]), reverse=True)
    na_ranking = sorted(
        results,
        key=lambda item: abs(item["sensitivity"]["urinary_na_equiv_umol_min_per_kidney"]),
        reverse=True,
    )
    return {
        "baseline_terminal_output": base_output,
        "perturbations": results,
        "flow_ranking": [item["parameter"] for item in flow_ranking],
        "na_ranking": [item["parameter"] for item in na_ranking],
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def build_results(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    baseline = simulate(parameters)
    sensitivity = run_sensitivity(parameters)
    null_a = simulate(parameters)
    null_b = simulate(parameters)
    null_difference = max(
        abs(null_a["terminal_output"][key] - null_b["terminal_output"][key])
        for key in null_a["terminal_output"]
    )
    limit = zero_transport_limit(parameters)
    unit = unit_check(parameters)
    predicted_flow_winner = "CD_Lp"
    predicted_na_winner = "TAL_vmax"
    flow_winner = sensitivity["flow_ranking"][0]
    na_winner = sensitivity["na_ranking"][0]
    acceptance = {
        "unit_check": unit["status"] == "PASS",
        "zero_transport_limit": limit["relative_error"] < 1.0e-8,
        "null_perturbation": null_difference < 1.0e-12,
        "flow_ranking": flow_winner == predicted_flow_winner,
        "na_ranking": na_winner == predicted_na_winner,
        "finite_nonnegative_outputs": all(
            np.isfinite(value) and value >= 0.0
            for value in baseline["terminal_output"].values()
        ),
        "water_mass_residual": abs(baseline["mass_balance"]["relative_total_water_balance_residual"]) < 1.0e-6,
    }
    acceptance["all_pass"] = all(acceptance.values())
    results = {
        "id": "BT-HX-Q019",
        "status": "PASS" if acceptance["all_pass"] else "FAIL",
        "model_scope": "one superficial nephron scaled to one kidney; dynamic well-mixed segment compartments",
        "observable_definition": {
            "independent_measurements": [
                "terminal urine flow",
                "urinary NaCl-equivalent amount",
                "urinary urea amount",
                "terminal urine osmolality",
            ],
            "latent_states_not_used_as_measurements": list(SEGMENT_NAMES),
        },
        "source_anchor": {
            "citation": "Layton AT, Layton HE (2019), PLOS Computational Biology 15:e1006108",
            "doi": "10.1371/journal.pcbi.1006108",
            "table_or_figure": "Tables 1-3 and model-parameter text",
            "reference_flow_ml_min_per_kidney": 0.62,
            "reference_urea_mM": 280.0,
        },
        "parameters": parameter_table(parameters),
        "unit_check": unit,
        "baseline": {
            "terminal_output": baseline["terminal_output"],
            "solver": baseline["solver"],
            "mass_balance": baseline["mass_balance"],
            "time_course": [
                {
                    "time_min": frame["time_min"],
                    "urine_flow_ml_min_per_kidney": frame["terminal_flow_ml_min_per_kidney"],
                    "urinary_na_equiv_umol_min_per_kidney": frame["terminal_na_umol_min_per_kidney"],
                    "urinary_urea_umol_min_per_kidney": frame["terminal_urea_umol_min_per_kidney"],
                    "urine_osm_mM": frame["terminal_osm_mM"],
                }
                for frame in baseline["frames"][::40]
            ],
        },
        "sensitivity": sensitivity,
        "controls": {
            "null_max_absolute_output_difference": null_difference,
            "zero_transport_limit": limit,
        },
        "acceptance": acceptance,
        "next_resolution_step": "Calibrate segment-specific water and NaCl conductances against time-resolved urine volume, sodium, and osmolality measurements; then replace the well-mixed compartments with an axial countercurrent model and add nephron-population heterogeneity.",
    }
    return _jsonable(results)


def write_results(path: str | Path = "results.json", parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    output_path = Path(path)
    result = build_results(parameters)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    result = write_results()
    print(json.dumps({"status": result["status"], "acceptance": result["acceptance"]}, sort_keys=True))
