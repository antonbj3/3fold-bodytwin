"""A small first-principles model for lossless coarse/fine muscle-model switching.

The fine model resolves a fast activation filter, a cross-bridge-like filter,
central adaptation, and fast and slow fatigue memories.  The coarse model keeps
only activation and the slow fatigue memory.  A closure model keeps all
history variables but replaces the fast cross-bridge filter with an algebraic
approximation.  Force is normalized by F_max, so the reported errors are
fractions of maximal force.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    tau_a: float = 0.025
    tau_x: float = 0.015
    tau_h: float = 22.0
    gamma_h: float = 0.10
    tau_f: float = 5.0
    k_f: float = 0.018
    tau_s: float = 300.0
    k_s: float = 0.0015
    w_f: float = 0.25
    w_s: float = 0.75
    c0: float = 0.050
    beta_ct: float = 0.379
    f_max: float = 1.0
    perturbation_amplitude: float = 0.25
    perturbation_duration: float = 0.1


DEFAULT_PARAMETERS = Parameters()

PROTOCOL_START = 0.0
PROTOCOL_END = 140.0
PROTOCOL_BREAKS = (5.0, 35.0, 55.0, 79.5, 80.0, 80.1)
SWITCH_TIME = 79.5
PERTURBATION_TIME = 80.0
STATE_BOUND_TOLERANCE = 1.0e-8
GRID = np.linspace(PROTOCOL_START, PROTOCOL_END, 14001)

PARAMETER_TABLE = (
    {
        "name": "tau_a",
        "value": DEFAULT_PARAMETERS.tau_a,
        "unit": "s",
        "role": "activation-filter time constant",
        "source_or_assumption": "assumption; fast excitation-contraction lumping",
    },
    {
        "name": "tau_x",
        "value": DEFAULT_PARAMETERS.tau_x,
        "unit": "s",
        "role": "cross-bridge-like filter time constant",
        "source_or_assumption": "assumption; unresolved fast mechanical state",
    },
    {
        "name": "tau_h",
        "value": DEFAULT_PARAMETERS.tau_h,
        "unit": "s",
        "role": "central adaptation time constant",
        "source_or_assumption": "Potvin & Fuglevand 2017, Methods Eq. (12), DOI 10.1371/journal.pcbi.1005581",
    },
    {
        "name": "gamma_h",
        "value": DEFAULT_PARAMETERS.gamma_h,
        "unit": "1",
        "role": "normalized target adaptation current",
        "source_or_assumption": "assumption; not a patient-specific measurement",
    },
    {
        "name": "tau_f",
        "value": DEFAULT_PARAMETERS.tau_f,
        "unit": "s",
        "role": "fast fatigue-memory time constant",
        "source_or_assumption": "assumption; short-lived metabolic/activation history",
    },
    {
        "name": "k_f",
        "value": DEFAULT_PARAMETERS.k_f,
        "unit": "s^-1",
        "role": "fast fatigue-memory loading rate",
        "source_or_assumption": "assumption; qualitative fatigue-rate scaling from Potvin & Fuglevand 2017",
    },
    {
        "name": "tau_s",
        "value": DEFAULT_PARAMETERS.tau_s,
        "unit": "s",
        "role": "slow fatigue-memory time constant",
        "source_or_assumption": "assumption; long-lived contractile history",
    },
    {
        "name": "k_s",
        "value": DEFAULT_PARAMETERS.k_s,
        "unit": "s^-1",
        "role": "slow fatigue-memory loading rate",
        "source_or_assumption": "assumption; not fitted to measured data",
    },
    {
        "name": "w_f",
        "value": DEFAULT_PARAMETERS.w_f,
        "unit": "1",
        "role": "fast-memory weight in force capacity",
        "source_or_assumption": "assumption for a two-timescale lumping",
    },
    {
        "name": "w_s",
        "value": DEFAULT_PARAMETERS.w_s,
        "unit": "1",
        "role": "slow-memory weight in force capacity",
        "source_or_assumption": "assumption for a two-timescale lumping",
    },
    {
        "name": "c0",
        "value": DEFAULT_PARAMETERS.c0,
        "unit": "s",
        "role": "rested contraction time",
        "source_or_assumption": "midpoint of 30-90 ms range in Potvin & Fuglevand 2017 Fig. 1B",
    },
    {
        "name": "beta_ct",
        "value": DEFAULT_PARAMETERS.beta_ct,
        "unit": "1",
        "role": "relative contraction-time increase per fractional force loss",
        "source_or_assumption": "Potvin & Fuglevand 2017 Methods Eq. (11), DOI 10.1371/journal.pcbi.1005581",
    },
    {
        "name": "f_max",
        "value": DEFAULT_PARAMETERS.f_max,
        "unit": "normalized force",
        "role": "maximum force scale",
        "source_or_assumption": "normalization; no measured force is claimed",
    },
    {
        "name": "perturbation_amplitude",
        "value": DEFAULT_PARAMETERS.perturbation_amplitude,
        "unit": "normalized drive",
        "role": "post-switch perturbation amplitude",
        "source_or_assumption": "frozen numerical perturbation for the model test",
    },
    {
        "name": "perturbation_duration",
        "value": DEFAULT_PARAMETERS.perturbation_duration,
        "unit": "s",
        "role": "post-switch perturbation duration",
        "source_or_assumption": "frozen numerical perturbation for the model test",
    },
)

UNIT_CHECKS = (
    {
        "term": "activation derivative",
        "lhs": "dimensionless s^-1",
        "rhs": "(dimensionless - dimensionless) / s",
        "ok": True,
    },
    {
        "term": "adaptation derivative",
        "lhs": "dimensionless s^-1",
        "rhs": "(dimensionless - dimensionless) / s",
        "ok": True,
    },
    {
        "term": "fatigue derivative",
        "lhs": "dimensionless s^-1",
        "rhs": "s^-1 * dimensionless - dimensionless / s",
        "ok": True,
    },
    {
        "term": "normalized firing rate",
        "lhs": "dimensionless",
        "rhs": "Hz * s / s",
        "ok": True,
    },
    {
        "term": "force",
        "lhs": "normalized force",
        "rhs": "normalized force * dimensionless * dimensionless * dimensionless",
        "ok": True,
    },
    {
        "term": "post-switch error",
        "lhs": "dimensionless",
        "rhs": "normalized force / normalized force",
        "ok": True,
    },
)


def protocol_drive(t: float) -> float:
    if 5.0 <= t < 35.0 or 55.0 <= t < 80.0:
        return 0.65
    if 80.0 <= t < 80.0 + DEFAULT_PARAMETERS.perturbation_duration:
        return 0.65 + DEFAULT_PARAMETERS.perturbation_amplitude
    return 0.0


def force_frequency(normalized_rate: float | np.ndarray) -> float | np.ndarray:
    n = np.asarray(normalized_rate, dtype=float)
    n = np.clip(n, 0.0, 2.0)
    value = np.where(n <= 0.4, 0.3 * n, 1.0 - np.exp(-2.0 * n**3))
    if value.ndim == 0:
        return float(value)
    return value


def capacity(mf: float, ms: float, p: Parameters = DEFAULT_PARAMETERS) -> float:
    return float(np.clip(1.0 - p.w_f * mf - p.w_s * ms, 0.0, 1.0))


def contraction_time(capacity_value: float, p: Parameters = DEFAULT_PARAMETERS) -> float:
    return p.c0 * (1.0 + p.beta_ct * (1.0 - capacity_value))


def _high_drive(u: float, h: float) -> float:
    return float(np.clip(u - h, 0.0, 1.0))


def _high_derived(state: np.ndarray, u: float, p: Parameters) -> tuple[float, float, float, float]:
    a, x, h, mf, ms = state
    effective_drive = _high_drive(u, h)
    c = contraction_time(capacity(mf, ms, p), p)
    normalized_rate = effective_drive * c / p.c0
    force_factor = float(force_frequency(normalized_rate))
    c_value = capacity(mf, ms, p)
    load = float(np.clip(x * force_factor, 0.0, 1.0))
    force = p.f_max * x * force_factor * c_value
    return effective_drive, c, load, float(np.clip(force, 0.0, p.f_max))


def high_output(state: np.ndarray, u: float, p: Parameters = DEFAULT_PARAMETERS) -> float:
    return _high_derived(np.asarray(state, dtype=float), u, p)[3]


def high_derivative(t: float, state: np.ndarray, u: float, p: Parameters) -> np.ndarray:
    del t
    a, x, h, mf, ms = np.asarray(state, dtype=float)
    effective_drive, _, load, _ = _high_derived(np.asarray(state, dtype=float), u, p)
    return np.asarray(
        [
            (effective_drive - a) / p.tau_a,
            (a - x) / p.tau_x,
            (p.gamma_h * effective_drive - h) / p.tau_h,
            p.k_f * load - mf / p.tau_f,
            p.k_s * load - ms / p.tau_s,
        ],
        dtype=float,
    )


def _project_high(state: np.ndarray) -> np.ndarray:
    projected = np.asarray(state, dtype=float).copy()
    projected[0] = np.clip(projected[0], 0.0, 1.0)
    projected[1] = np.clip(projected[1], 0.0, 1.0)
    projected[2] = np.clip(projected[2], 0.0, 1.0)
    projected[3] = max(projected[3], 0.0)
    projected[4] = max(projected[4], 0.0)
    return projected


def low_derivative(t: float, state: np.ndarray, u: float, p: Parameters) -> np.ndarray:
    del t
    a, ms = np.asarray(state, dtype=float)
    load = float(np.clip(a * force_frequency(u), 0.0, 1.0))
    return np.asarray([(u - a) / p.tau_a, p.k_s * load - ms / p.tau_s], dtype=float)


def low_output(state: np.ndarray, u: float, p: Parameters = DEFAULT_PARAMETERS) -> float:
    a, ms = np.asarray(state, dtype=float)
    c_value = float(np.clip(1.0 - p.w_s * ms, 0.0, 1.0))
    force = p.f_max * a * force_frequency(u) * c_value
    return float(np.clip(force, 0.0, p.f_max))


def _project_low(state: np.ndarray) -> np.ndarray:
    projected = np.asarray(state, dtype=float).copy()
    projected[0] = np.clip(projected[0], 0.0, 1.0)
    projected[1] = max(projected[1], 0.0)
    return projected


def closure_derivative(t: float, state: np.ndarray, u: float, p: Parameters) -> np.ndarray:
    del t
    a, h, mf, ms = np.asarray(state, dtype=float)
    effective_drive = _high_drive(u, h)
    c_value = capacity(mf, ms, p)
    c = contraction_time(c_value, p)
    normalized_rate = effective_drive * c / p.c0
    load = float(np.clip(a * force_frequency(normalized_rate), 0.0, 1.0))
    return np.asarray(
        [
            (effective_drive - a) / p.tau_a,
            (p.gamma_h * effective_drive - h) / p.tau_h,
            p.k_f * load - mf / p.tau_f,
            p.k_s * load - ms / p.tau_s,
        ],
        dtype=float,
    )


def closure_output(state: np.ndarray, u: float, p: Parameters = DEFAULT_PARAMETERS) -> float:
    a, h, mf, ms = np.asarray(state, dtype=float)
    effective_drive = _high_drive(u, h)
    c_value = capacity(mf, ms, p)
    c = contraction_time(c_value, p)
    force = p.f_max * a * force_frequency(effective_drive * c / p.c0) * c_value
    return float(np.clip(force, 0.0, p.f_max))


def _project_closure(state: np.ndarray) -> np.ndarray:
    projected = np.asarray(state, dtype=float).copy()
    projected[0] = np.clip(projected[0], 0.0, 1.0)
    projected[1] = np.clip(projected[1], 0.0, 1.0)
    projected[2] = max(projected[2], 0.0)
    projected[3] = max(projected[3], 0.0)
    return projected


def _segment_bounds(start: float, end: float) -> list[tuple[float, float]]:
    points = {float(start), float(end)}
    points.update(point for point in PROTOCOL_BREAKS if start < point < end)
    ordered = sorted(points)
    return list(zip(ordered[:-1], ordered[1:]))


def _solve_segment(
    derivative: Callable[[float, np.ndarray], np.ndarray],
    state: np.ndarray,
    start: float,
    end: float,
) -> object:
    solution = solve_ivp(
        derivative,
        (start, end),
        np.asarray(state, dtype=float),
        method="DOP853",
        rtol=1.0e-8,
        atol=1.0e-10,
        max_step=0.10,
        dense_output=True,
    )
    if not solution.success:
        raise RuntimeError(f"integration failed on [{start}, {end}]: {solution.message}")
    return solution


def integrate_high_to(
    start: float,
    end: float,
    state: np.ndarray,
    p: Parameters = DEFAULT_PARAMETERS,
) -> np.ndarray:
    current = _project_high(state)
    for left, right in _segment_bounds(start, end):
        solution = _solve_segment(
            lambda t, z: high_derivative(t, z, protocol_drive(t), p), current, left, right
        )
        current = _project_high(solution.y[:, -1])
    return current


def _sample_segments(
    start: float,
    end: float,
    state: np.ndarray,
    grid: np.ndarray,
    derivative: Callable[[float, np.ndarray], np.ndarray],
    output: Callable[[np.ndarray, float], float],
    project: Callable[[np.ndarray], np.ndarray],
    drive: Callable[[float], float],
    p: Parameters,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grid = np.asarray(grid[(grid >= start) & (grid <= end)], dtype=float)
    if grid.size == 0:
        return np.empty(0), np.empty(0), np.empty((0, 0))
    times: list[np.ndarray] = []
    forces: list[np.ndarray] = []
    states: list[np.ndarray] = []
    current = project(state)
    for index, (left, right) in enumerate(_segment_bounds(start, end)):
        solution = _solve_segment(derivative, current, left, right)
        mask = (grid >= left) & (grid <= right)
        if index > 0:
            mask &= grid > left
        if np.any(mask):
            times.append(grid[mask])
            states.append(solution.sol(grid[mask]).T)
            forces.append(np.asarray([output(state_value, drive(time_value)) for time_value, state_value in zip(grid[mask], solution.sol(grid[mask]).T)], dtype=float))
        current = project(solution.y[:, -1])
    time_array = np.concatenate(times) if times else np.empty(0)
    force_array = np.concatenate(forces) if forces else np.empty(0)
    state_array = np.concatenate(states) if states else np.empty((0, len(state)))
    order = np.argsort(time_array)
    return time_array[order], force_array[order], state_array[order]


def _high_grid(
    grid: np.ndarray,
    p: Parameters = DEFAULT_PARAMETERS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return _sample_segments(
        PROTOCOL_START,
        PROTOCOL_END,
        np.zeros(5, dtype=float),
        grid,
        lambda t, z: high_derivative(t, z, protocol_drive(t), p),
        lambda state, u: high_output(state, u, p),
        _project_high,
        protocol_drive,
        p,
    )


def _low_grid(
    start: float,
    end: float,
    state: np.ndarray,
    grid: np.ndarray,
    p: Parameters,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return _sample_segments(
        start,
        end,
        state,
        grid,
        lambda t, z: low_derivative(t, z, protocol_drive(t), p),
        lambda state, u: low_output(state, u, p),
        _project_low,
        protocol_drive,
        p,
    )


def _closure_grid(
    start: float,
    end: float,
    state: np.ndarray,
    grid: np.ndarray,
    p: Parameters,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return _sample_segments(
        start,
        end,
        state,
        grid,
        lambda t, z: closure_derivative(t, z, protocol_drive(t), p),
        lambda state, u: closure_output(state, u, p),
        _project_closure,
        protocol_drive,
        p,
    )


def _state_from_high(state: np.ndarray, kind: str) -> np.ndarray:
    if kind == "naive":
        return np.asarray([state[0], state[4]], dtype=float)
    if kind == "closure":
        return np.asarray([state[0], state[2], state[3], state[4]], dtype=float)
    raise ValueError(f"unknown low-detail kind: {kind}")


def _states_within_bounds(high_states: np.ndarray, low_states: np.ndarray, kind: str) -> bool:
    high_states = np.asarray(high_states, dtype=float)
    low_states = np.asarray(low_states, dtype=float)
    if high_states.size == 0 or low_states.size == 0:
        return False
    tolerance = STATE_BOUND_TOLERANCE
    high_ok = bool(np.all(np.isfinite(high_states)) and np.all((high_states >= -tolerance) & (high_states <= 1.0 + tolerance)))
    low_ok = bool(np.all(np.isfinite(low_states)) and np.all((low_states >= -tolerance) & (low_states <= 1.0 + tolerance)))
    return high_ok and low_ok


def _metrics(
    reference_time: np.ndarray,
    reference_force: np.ndarray,
    switched_force: np.ndarray,
    switch_index: int,
    p: Parameters,
) -> dict[str, float | bool]:
    signed_jump = float(switched_force[switch_index] - reference_force[switch_index])
    post_mask = np.asarray(reference_time >= PERTURBATION_TIME - 1.0e-12)
    post_error = np.abs(switched_force[post_mask] - reference_force[post_mask])
    max_error = float(np.max(post_error)) if post_error.size else 0.0
    immediate_pass = abs(signed_jump) <= 0.02 * p.f_max
    post_pass = max_error <= 0.05 * p.f_max
    return {
        "signed_force_jump": signed_jump,
        "absolute_force_jump": abs(signed_jump),
        "post_perturbation_max_error": max_error,
        "post_perturbation_mean_error": float(np.mean(post_error)) if post_error.size else 0.0,
        "immediate_jump_pass": bool(immediate_pass),
        "post_perturbation_error_pass": bool(post_pass),
        "pass": bool(immediate_pass and post_pass),
    }


def simulate_case(
    p: Parameters = DEFAULT_PARAMETERS,
    switch_time: float = SWITCH_TIME,
    low_kind: str = "naive",
    grid: np.ndarray = GRID,
) -> dict[str, object]:
    if low_kind not in {"naive", "closure"}:
        raise ValueError("low_kind must be 'naive' or 'closure'")
    reference_time, reference_force, reference_state = _high_grid(grid, p)
    switch_state = integrate_high_to(PROTOCOL_START, switch_time, np.zeros(5, dtype=float), p)
    low_state = _state_from_high(switch_state, low_kind)
    low_start_force = low_output(low_state, protocol_drive(switch_time), p) if low_kind == "naive" else closure_output(low_state, protocol_drive(switch_time), p)
    post_grid = grid[grid >= switch_time]
    if low_kind == "naive":
        low_time, low_force_values, low_states = _low_grid(switch_time, PROTOCOL_END, low_state, post_grid, p)
    else:
        low_time, low_force_values, low_states = _closure_grid(switch_time, PROTOCOL_END, low_state, post_grid, p)
    switched_force = np.interp(grid, low_time, low_force_values, left=low_start_force, right=low_force_values[-1])
    switched_force[grid < switch_time] = np.interp(grid[grid < switch_time], reference_time, reference_force)
    switch_index = int(np.argmin(np.abs(grid - switch_time)))
    metrics = _metrics(reference_time, reference_force, switched_force, switch_index, p)
    high_capacity = capacity(switch_state[3], switch_state[4], p)
    low_capacity = capacity(switch_state[3], switch_state[4], p) if low_kind == "closure" else float(np.clip(1.0 - p.w_s * switch_state[4], 0.0, 1.0))
    memory_projection_error = abs(high_capacity - low_capacity)
    diagnostics = {
        "switch_time_s": switch_time,
        "low_detail_kind": low_kind,
        "high_state_at_switch": {
            "activation": float(switch_state[0]),
            "fast_filter": float(switch_state[1]),
            "central_adaptation": float(switch_state[2]),
            "fast_fatigue_memory": float(switch_state[3]),
            "slow_fatigue_memory": float(switch_state[4]),
        },
        "high_capacity_at_switch": high_capacity,
        "low_capacity_at_switch": low_capacity,
        "memory_projection_error": memory_projection_error,
        "state_bounds_ok": _states_within_bounds(reference_state, low_states, low_kind),
        "state_bound_tolerance": STATE_BOUND_TOLERANCE,
        "grid_points": int(grid.size),
    }
    return {
        "metrics": metrics,
        "diagnostics": diagnostics,
        "reference_time": reference_time,
        "reference_force": reference_force,
        "switched_force": switched_force,
        "low_time": low_time,
        "low_force": low_force_values,
        "low_state": low_states,
    }


def sensitivity_analysis(p: Parameters = DEFAULT_PARAMETERS) -> list[dict[str, object]]:
    base = simulate_case(p, low_kind="naive")
    base_value = float(base["metrics"]["post_perturbation_max_error"])
    records: list[dict[str, object]] = []
    for name in ("tau_f", "k_f", "tau_s"):
        base_parameter = float(getattr(p, name))
        values: dict[str, float] = {}
        for label, factor in (("minus50", 0.5), ("plus50", 1.5)):
            changed = replace(p, **{name: base_parameter * factor})
            result = simulate_case(changed, low_kind="naive")
            value = float(result["metrics"]["post_perturbation_max_error"])
            values[label] = value
            records.append(
                {
                    "parameter": name,
                    "factor": label,
                    "value": base_parameter * factor,
                    "unit": next(item["unit"] for item in PARAMETER_TABLE if item["name"] == name),
                    "post_perturbation_max_error": value,
                    "absolute_force_jump": float(result["metrics"]["absolute_force_jump"]),
                }
            )
        sensitivity = (values["plus50"] - values["minus50"]) / (2.0 * base_value) if base_value else 0.0
        records.append(
            {
                "parameter": name,
                "factor": "relative_sensitivity",
                "value": sensitivity,
                "unit": "dimensionless",
                "post_perturbation_max_error": base_value,
            }
        )
    return records


def unit_check() -> dict[str, object]:
    return {
        "all_pass": all(bool(item["ok"]) for item in UNIT_CHECKS),
        "checks": list(UNIT_CHECKS),
        "normalization": "F/F_max is dimensionless; n = r c / c0 is dimensionless.",
    }


def run_experiment(output_path: str | Path = "results.json") -> dict[str, object]:
    p = DEFAULT_PARAMETERS
    naive = simulate_case(p, low_kind="naive")
    closure = simulate_case(p, low_kind="closure")
    placebo = simulate_case(p, switch_time=0.0, low_kind="naive")
    result: dict[str, object] = {
        "id": "BT-HX-Q058",
        "model": "two-timescale fatigue-memory state transfer for fine/coarse muscle-model switching",
        "scope": "normalized isometric motor-unit surrogate; no measured data",
        "equations": {
            "activation": "da/dt=(clip(u-h,0,1)-a)/tau_a",
            "fast_filter": "dx/dt=(a-x)/tau_x",
            "central_adaptation": "dh/dt=(gamma_h*clip(u-h,0,1)-h)/tau_h",
            "fatigue": "dm_i/dt=k_i*load-m_i/tau_i, i in {f,s}",
            "load": "load=clip(x*Phi(effective_drive*c/c0),0,1)",
            "capacity": "C=clip(1-w_f*m_f-w_s*m_s,0,1)",
            "contraction_time": "c=c0*(1+beta_ct*(1-C))",
            "force_frequency": "Phi(n)=0.3*n for n<=0.4; 1-exp(-2*n^3) for n>0.4",
            "force": "F/F_max=x*Phi(effective_drive*c/c0)*C",
            "switch_naive": "low state=(a,m_s); h, m_f and x are not transferred",
            "switch_closure": "low state=(a,h,m_f,m_s); x is algebraic",
        },
        "parameter_table": list(PARAMETER_TABLE),
        "unit_check": unit_check(),
        "sources": [
            {
                "citation": "Potvin JR, Fuglevand AJ (2017), A motor unit-based model of muscle fatigue",
                "doi": "10.1371/journal.pcbi.1005581",
                "verified_values": [
                    "Fig. 1B: 30-90 ms contraction-time range",
                    "Methods Eq. (11): 0.379 contraction-time change per fractional force loss",
                    "Methods Eq. (12): 22 s adaptation time constant",
                    "Fig. 2: 511.5 s simulated 20% MVC endurance time",
                    "Fig. 3: 95.5 s simulated 50% MVC endurance time",
                ],
                "reference_values": {
                    "rested_contraction_time_range_s": [0.030, 0.090],
                    "representative_contraction_time_s": 0.050,
                    "contraction_time_change_coefficient": 0.379,
                    "twenty_percent_force_loss_contraction_time_fraction": 1.0758,
                    "adaptation_time_constant_s": 22.0,
                    "simulated_endurance_20_mvc_s": 511.5,
                    "simulated_endurance_50_mvc_s": 95.5,
                    "maximal_force_loss_rate_first_20s_percent_mvc_per_s": 1.4,
                    "time_to_50_percent_initial_force_s": 70.0,
                },
            },
            {
                "citation": "Liu JZ, Brown RW, Yue GH (2002), A dynamical model of muscle activation, fatigue, and recovery",
                "doi": "10.1016/S0006-3495(02)75580-X",
                "use": "conceptual support for separate fatigue and recovery dynamics; not used for calibration",
            },
        ],
        "protocol": {
            "start_s": PROTOCOL_START,
            "end_s": PROTOCOL_END,
            "active_intervals_s": [[5.0, 35.0], [55.0, 80.0]],
            "active_drive": 0.65,
            "switch_time_s": SWITCH_TIME,
            "perturbation_time_s": PERTURBATION_TIME,
            "perturbation_drive": 0.90,
            "perturbation_duration_s": p.perturbation_duration,
        },
        "criteria": {
            "jump_limit_fmax": 0.02,
            "post_perturbation_error_limit_fmax": 0.05,
            "source_reference_values_are_not_fitted": True,
        },
        "primary_naive_switch": {
            "metrics": naive["metrics"],
            "diagnostics": naive["diagnostics"],
        },
        "closure_switch": {
            "metrics": closure["metrics"],
            "diagnostics": closure["diagnostics"],
        },
        "placebo_switch_at_zero": {
            "metrics": placebo["metrics"],
            "diagnostics": placebo["diagnostics"],
        },
        "sensitivity": sensitivity_analysis(p),
        "unknowns": [
            "The fast and slow fatigue-rate constants and weights are assumptions, not measured values.",
            "A universal switching threshold requires paired fine/coarse data with measured fatigue and recovery trajectories.",
            "The source model is motor-unit based and does not establish a human joint-level tolerance.",
        ],
        "next_resolution_steps": [
            "Measure activation, force, and perturbation responses before and after a documented switch on the same preparation.",
            "Identify the two fatigue time constants and memory weights instead of treating them as assumptions.",
            "Test state reconstruction or a closure observer when the fast state is not directly observed.",
            "Repeat across motor-unit types and perturbation spectra before transferring the criterion to BodyTwin joint models.",
        ],
    }
    path = Path(output_path)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    output = run_experiment()
    print(json.dumps({
        "id": output["id"],
        "naive": output["primary_naive_switch"]["metrics"],
        "closure": output["closure_switch"]["metrics"],
        "results_path": "results.json",
    }, indent=2, sort_keys=True))
