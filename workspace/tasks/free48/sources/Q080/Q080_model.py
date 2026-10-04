from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Parameters:
    domain_length_um: float = 4.0
    grid_points: int = 21
    mean_receptor_density: float = 55.0
    patch_centers_um: tuple[tuple[float, float], ...] = (
        (1.2, 1.2),
        (2.8, 1.2),
        (1.2, 2.8),
        (2.8, 2.8),
    )
    patch_sigma_um: float = 0.30
    ligand_center_um: tuple[float, float] = (2.0, 2.0)
    ligand_sigma_um: float = 0.35
    pulse_duration_s: float = 0.10
    source_rate_nM_s: float = 15000.0
    ligand_diffusion_um2_s: float = 25.0
    clearance_rate_s: float = 0.5
    kon_s1_nM1: float = 0.02
    koff_s1: float = 0.66
    receptor_diffusion_um2_s: float = 0.04
    duration_s: float = 1.2
    dt_s: float = 0.00005


def parameter_table(p: Parameters | None = None) -> list[dict[str, Any]]:
    p = p or Parameters()
    dx = p.domain_length_um / p.grid_points
    return [
        {
            "name": "domain_length_um",
            "value": p.domain_length_um,
            "unit": "um",
            "source_or_assumption": "assumption: square effective membrane domain",
        },
        {
            "name": "grid_points",
            "value": p.grid_points,
            "unit": "cells per axis",
            "source_or_assumption": "assumption: finite-volume grid",
        },
        {
            "name": "dx",
            "value": dx,
            "unit": "um",
            "source_or_assumption": "derived: domain_length_um / grid_points",
        },
        {
            "name": "mean_receptor_density",
            "value": p.mean_receptor_density,
            "unit": "receptor um^-2",
            "source_or_assumption": "Briddon et al. 2004, Fig. 6b and text: approximately 55 receptor um^-2",
        },
        {
            "name": "patch_centers_um",
            "value": [list(center) for center in p.patch_centers_um],
            "unit": "um",
            "source_or_assumption": "assumption: four fixed patch centers",
        },
        {
            "name": "patch_sigma_um",
            "value": p.patch_sigma_um,
            "unit": "um",
            "source_or_assumption": "assumption; motivated by Suzuki et al. 2005 Fig. 4d 210 nm compartment scale",
        },
        {
            "name": "ligand_center_um",
            "value": list(p.ligand_center_um),
            "unit": "um",
            "source_or_assumption": "assumption: fixed pulse center",
        },
        {
            "name": "ligand_sigma_um",
            "value": p.ligand_sigma_um,
            "unit": "um",
            "source_or_assumption": "assumption: localized Gaussian ligand source",
        },
        {
            "name": "pulse_duration_s",
            "value": p.pulse_duration_s,
            "unit": "s",
            "source_or_assumption": "assumption: controlled temporal pulse",
        },
        {
            "name": "source_rate_nM_s",
            "value": p.source_rate_nM_s,
            "unit": "nM s^-1",
            "source_or_assumption": "assumption: pulse amplitude, not measured data",
        },
        {
            "name": "ligand_diffusion_um2_s",
            "value": p.ligand_diffusion_um2_s,
            "unit": "um^2 s^-1",
            "source_or_assumption": "effective 2D assumption; Briddon et al. free-ligand 251 um^2 s^-1 is a 3D/FCS anchor",
        },
        {
            "name": "clearance_rate_s",
            "value": p.clearance_rate_s,
            "unit": "s^-1",
            "source_or_assumption": "assumption: homogeneous post-pulse effective loss",
        },
        {
            "name": "kon_s1_nM1",
            "value": p.kon_s1_nM1,
            "unit": "s^-1 nM^-1",
            "source_or_assumption": "assumption; paired with koff to reproduce Briddon Kd anchor",
        },
        {
            "name": "koff_s1",
            "value": p.koff_s1,
            "unit": "s^-1",
            "source_or_assumption": "assumption; koff/kon equals 33 nM",
        },
        {
            "name": "receptor_diffusion_um2_s",
            "value": p.receptor_diffusion_um2_s,
            "unit": "um^2 s^-1",
            "source_or_assumption": "Calebiro et al. 2013 Fig. 2I beta2AR median 0.039 um^2 s^-1; rounded anchor",
        },
        {
            "name": "duration_s",
            "value": p.duration_s,
            "unit": "s",
            "source_or_assumption": "assumption: observation endpoint",
        },
        {
            "name": "dt_s",
            "value": p.dt_s,
            "unit": "s",
            "source_or_assumption": "assumption: explicit finite-difference step",
        },
    ]


def unit_check(p: Parameters | None = None) -> dict[str, Any]:
    p = p or Parameters()
    equations = {
        "ligand_diffusion": {
            "term": "D_L Laplacian(L)",
            "dimensions": "um^2 s^-1 * nM um^-2 = nM s^-1",
            "matches": "partial_t(L) = nM s^-1",
        },
        "receptor_diffusion": {
            "term": "D_R Laplacian(r)",
            "dimensions": "um^2 s^-1 * um^-2 = s^-1",
            "matches": "partial_t(r) = s^-1",
        },
        "association": {
            "term": "kon L r",
            "dimensions": "s^-1 nM^-1 * nM * 1 = s^-1",
            "matches": "partial_t(c) reaction term = s^-1",
        },
        "dissociation": {
            "term": "koff c",
            "dimensions": "s^-1 * 1 = s^-1",
            "matches": "partial_t(c) reaction term = s^-1",
        },
        "ligand_source": {
            "term": "source_rate Gaussian",
            "dimensions": "nM s^-1 * 1 = nM s^-1",
            "matches": "partial_t(L) source term = nM s^-1",
        },
        "ligand_clearance": {
            "term": "clearance_rate L",
            "dimensions": "s^-1 * nM = nM s^-1",
            "matches": "partial_t(L) loss term = nM s^-1",
        },
        "receptor_density": {
            "term": "Rbar (r+c)",
            "dimensions": "receptor um^-2 * 1 = receptor um^-2",
            "matches": "R = receptor um^-2",
        },
        "finite_volume_diffusion": {
            "term": "dt D field / dx^2",
            "dimensions": "s * um^2 s^-1 / um^2 = 1",
            "matches": "dimensionless finite-volume update",
        },
        "grid_spacing": {
            "term": "dx and dt",
            "dimensions": "um and s",
            "matches": "spatial and temporal coordinates",
        },
        "activation_density": {
            "term": "Rbar c",
            "dimensions": "receptor um^-2 * 1 = receptor um^-2",
            "matches": "A = receptor um^-2",
        },
        "receptor_count": {
            "term": "Rbar (r+c) dA",
            "dimensions": "receptor um^-2 * um^2 = receptor",
            "matches": "N_total = receptor",
        },
    }
    return {
        "consistent": True,
        "equations": equations,
        "grid_unit": "um",
        "time_unit": "s",
        "ligand_unit": "nM",
        "receptor_density_unit": "receptor um^-2",
    }


def stable_step_count(duration_s: float, dt_s: float) -> int:
    ratio = duration_s / dt_s
    nearest_integer = int(round(ratio))
    if np.isclose(ratio, nearest_integer, rtol=1e-12, atol=1e-12):
        return nearest_integer
    return int(np.ceil(ratio))


def validate_parameters(p: Parameters) -> tuple[float, int, int]:
    if p.domain_length_um <= 0 or p.grid_points < 5:
        raise ValueError("domain_length_um and grid_points are invalid")
    if p.grid_points % 2 == 0:
        raise ValueError("grid_points must be odd for a centered patch")
    if p.mean_receptor_density <= 0 or p.patch_sigma_um <= 0 or p.ligand_sigma_um <= 0:
        raise ValueError("density and spatial widths must be positive")
    if p.pulse_duration_s < 0 or p.duration_s <= 0 or p.dt_s <= 0:
        raise ValueError("time parameters are invalid")
    if p.kon_s1_nM1 <= 0 or p.koff_s1 <= 0 or p.clearance_rate_s < 0:
        raise ValueError("kinetic parameters are invalid")
    if p.ligand_diffusion_um2_s < 0 or p.receptor_diffusion_um2_s < 0:
        raise ValueError("diffusion parameters are invalid")
    dx = p.domain_length_um / p.grid_points
    cfl = 4.0 * max(p.ligand_diffusion_um2_s, p.receptor_diffusion_um2_s) * p.dt_s / dx**2
    if cfl > 0.95:
        raise ValueError(f"explicit diffusion CFL={cfl:.6g} exceeds 0.95")
    n_steps = stable_step_count(p.duration_s, p.dt_s)
    dt = p.duration_s / n_steps
    pulse_steps = (
        0
        if p.pulse_duration_s == 0.0
        else stable_step_count(p.pulse_duration_s, dt)
    )
    if pulse_steps > n_steps:
        raise ValueError("pulse duration exceeds simulation duration")
    return dt, n_steps, pulse_steps


def grid(p: Parameters) -> tuple[np.ndarray, np.ndarray, float]:
    dx = p.domain_length_um / p.grid_points
    coordinates = (np.arange(p.grid_points, dtype=float) + 0.5) * dx
    x, y = np.meshgrid(coordinates, coordinates, indexing="xy")
    return x, y, dx


def receptor_pattern(p: Parameters, pattern: str) -> np.ndarray:
    x, y, _ = grid(p)
    if pattern == "uniform":
        relative = np.ones_like(x)
    elif pattern == "patched":
        weights = np.zeros_like(x)
        for center_x, center_y in p.patch_centers_um:
            distance_sq = (x - center_x) ** 2 + (y - center_y) ** 2
            weights += np.exp(-distance_sq / (2.0 * p.patch_sigma_um**2))
        relative = weights / float(np.mean(weights))
    else:
        raise ValueError(f"unknown pattern: {pattern}")
    relative = np.maximum(relative, 0.0)
    relative /= float(np.mean(relative))
    return relative


def ligand_source_field(p: Parameters) -> np.ndarray:
    x, y, _ = grid(p)
    center_x, center_y = p.ligand_center_um
    distance_sq = (x - center_x) ** 2 + (y - center_y) ** 2
    return np.exp(-distance_sq / (2.0 * p.ligand_sigma_um**2))


def laplacian(field: np.ndarray, dx: float) -> np.ndarray:
    padded = np.pad(field, 1, mode="edge")
    return (
        padded[:-2, 1:-1]
        + padded[2:, 1:-1]
        + padded[1:-1, :-2]
        + padded[1:-1, 2:]
        - 4.0 * padded[1:-1, 1:-1]
    ) / dx**2


def receptor_total_count(p: Parameters, pattern: str) -> float:
    _, _, dx = grid(p)
    relative = receptor_pattern(p, pattern)
    return float(p.mean_receptor_density * np.sum(relative) * dx**2)


def equilibrium_bound(p: Parameters, pattern: str, ligand_nM: float) -> np.ndarray:
    if ligand_nM < 0:
        raise ValueError("ligand_nM must be non-negative")
    occupancy = p.kon_s1_nM1 * ligand_nM / (p.koff_s1 + p.kon_s1_nM1 * ligand_nM)
    return p.mean_receptor_density * receptor_pattern(p, pattern) * occupancy


def _spatial_summary(p: Parameters, bound: np.ndarray) -> dict[str, float]:
    _, _, dx = grid(p)
    activation = p.mean_receptor_density * bound
    mean_value = float(np.mean(activation))
    total = float(np.sum(activation) * dx**2)
    peak = float(np.max(activation))
    cv = float(np.std(activation) / mean_value) if mean_value > 0 else 0.0
    peak_to_mean = float(peak / mean_value) if mean_value > 0 else 0.0
    return {
        "peak_activation": peak,
        "mean_activation": mean_value,
        "total_bound_receptors": total,
        "spatial_cv": cv,
        "peak_to_mean": peak_to_mean,
    }


def run_simulation(
    p: Parameters,
    pattern: str,
    record_fields: bool = False,
) -> dict[str, Any]:
    dt, n_steps, pulse_steps = validate_parameters(p)
    _, _, dx = grid(p)
    relative = receptor_pattern(p, pattern)
    free_state = relative.copy()
    bound_state = np.zeros_like(relative)
    ligand = np.zeros_like(relative)
    source_field = ligand_source_field(p)
    expected_count = p.mean_receptor_density * p.domain_length_um**2
    initial_count = float(p.mean_receptor_density * np.sum(free_state + bound_state) * dx**2)
    max_mass_error = abs(initial_count - expected_count) / expected_count
    min_free = float(np.min(free_state))
    min_bound = float(np.min(bound_state))
    min_ligand = float(np.min(ligand))
    total_bound_integral = 0.0
    pulse_peak_bound = np.zeros_like(relative)
    pulse_peak_step = -1
    pulse_end_bound = np.zeros_like(relative)
    pulse_end_ligand = np.zeros_like(relative)
    final_ligand_total = 0.0

    for step in range(n_steps):
        association = p.kon_s1_nM1 * ligand * free_state
        dissociation = p.koff_s1 * bound_state
        free_state += dt * (-association + dissociation)
        bound_state += dt * (association - dissociation)
        source = p.source_rate_nM_s * source_field if step < pulse_steps else 0.0
        ligand += dt * (source - p.clearance_rate_s * ligand)
        free_state += dt * p.receptor_diffusion_um2_s * laplacian(free_state, dx)
        bound_state += dt * p.receptor_diffusion_um2_s * laplacian(bound_state, dx)
        ligand += dt * p.ligand_diffusion_um2_s * laplacian(ligand, dx)

        min_free = min(min_free, float(np.min(free_state)))
        min_bound = min(min_bound, float(np.min(bound_state)))
        min_ligand = min(min_ligand, float(np.min(ligand)))
        if min_free < -1e-10 or min_bound < -1e-10 or min_ligand < -1e-10:
            raise RuntimeError("negative state encountered in explicit integration")
        current_count = p.mean_receptor_density * float(np.sum(free_state + bound_state) * dx**2)
        current_error = abs(current_count - expected_count) / expected_count
        max_mass_error = max(max_mass_error, current_error)
        total_bound_integral += float(np.sum(bound_state) * p.mean_receptor_density * dx**2) * dt
        if step < pulse_steps and float(np.max(bound_state)) > float(np.max(pulse_peak_bound)):
            pulse_peak_bound = bound_state.copy()
            pulse_peak_step = step
        if step + 1 == pulse_steps:
            pulse_end_bound = bound_state.copy()
            pulse_end_ligand = ligand.copy()

    final_count = p.mean_receptor_density * float(np.sum(free_state + bound_state) * dx**2)
    final_ligand_total = float(np.sum(ligand) * dx**2)
    result: dict[str, Any] = {
        "pattern": pattern,
        "dt_s": dt,
        "n_steps": n_steps,
        "pulse_steps": pulse_steps,
        "pulse_end_time_s": pulse_steps * dt,
        "duration_s": p.duration_s,
        "expected_total_receptors": expected_count,
        "initial_total_receptors": initial_count,
        "final_total_receptors": final_count,
        "max_relative_receptor_mass_error": float(max_mass_error),
        "minimum_free_state": min_free,
        "minimum_bound_state": min_bound,
        "minimum_ligand_nM": min_ligand,
        "pulse_peak_time_s": (pulse_peak_step + 1) * dt if pulse_peak_step >= 0 else None,
        "activation_time_integral_receptor_s": total_bound_integral,
        "pulse_end": _spatial_summary(p, pulse_end_bound),
        "pulse_peak": _spatial_summary(p, pulse_peak_bound),
        "final": _spatial_summary(p, bound_state),
        "source_integral_nM_um2": float(
            np.sum(source_field) * dx**2 * p.source_rate_nM_s * pulse_steps * dt
        ),
    }
    if record_fields:
        result["free_state_field"] = free_state
        result["bound_state_field"] = bound_state
        result["pulse_end_bound_field_receptor_um2"] = p.mean_receptor_density * pulse_end_bound
        result["pulse_end_ligand_field_nM"] = pulse_end_ligand
        result["final_bound_field_receptor_um2"] = p.mean_receptor_density * bound_state
        result["final_ligand_field_nM"] = ligand
    return result


def compare_patterns(p: Parameters) -> dict[str, Any]:
    uniform = run_simulation(p, "uniform")
    patched = run_simulation(p, "patched")
    uniform_peak = uniform["pulse_end"]["peak_activation"]
    patched_peak = patched["pulse_end"]["peak_activation"]
    uniform_total = uniform["pulse_end"]["total_bound_receptors"]
    patched_total = patched["pulse_end"]["total_bound_receptors"]
    uniform_auc = uniform["activation_time_integral_receptor_s"]
    patched_auc = patched["activation_time_integral_receptor_s"]
    peak_ratio = patched_peak / uniform_peak if uniform_peak > 0 else 0.0
    total_ratio = patched_total / uniform_total if uniform_total > 0 else 0.0
    auc_ratio = patched_auc / uniform_auc if uniform_auc > 0 else 0.0
    cv_increment = patched["pulse_end"]["spatial_cv"] - uniform["pulse_end"]["spatial_cv"]
    return {
        "uniform": uniform,
        "patched": patched,
        "peak_ratio": float(peak_ratio),
        "total_response_ratio_at_pulse_end": float(total_ratio),
        "activation_integral_ratio": float(auc_ratio),
        "cv_increment": float(cv_increment),
        "primary_pass": bool(peak_ratio >= 1.20 and cv_increment >= 0.10),
        "frozen_criteria": {
            "peak_ratio_threshold": 1.20,
            "cv_increment_threshold": 0.10,
        },
    }


def run_constant_ligand_limit(
    p: Parameters,
    pattern: str,
    ligand_nM: float,
    duration_s: float = 20.0,
    integration_dt_s: float = 0.001,
) -> tuple[np.ndarray, dict[str, float]]:
    if ligand_nM < 0 or duration_s <= 0 or integration_dt_s <= 0:
        raise ValueError("constant-ligand control parameters are invalid")
    reaction_rate = p.kon_s1_nM1 * ligand_nM + p.koff_s1
    if reaction_rate * integration_dt_s >= 1:
        raise ValueError("constant-ligand integration step violates positivity")
    n_steps = stable_step_count(duration_s, integration_dt_s)
    dt = duration_s / n_steps
    free_state = receptor_pattern(p, pattern)
    bound_state = np.zeros_like(free_state)
    expected_relative_count = float(np.sum(free_state))
    min_free = float(np.min(free_state))
    min_bound = float(np.min(bound_state))
    for _ in range(n_steps):
        association = p.kon_s1_nM1 * ligand_nM * free_state
        dissociation = p.koff_s1 * bound_state
        free_state += dt * (-association + dissociation)
        bound_state += dt * (association - dissociation)
        min_free = min(min_free, float(np.min(free_state)))
        min_bound = min(min_bound, float(np.min(bound_state)))
    final_relative_count = float(np.sum(free_state + bound_state))
    diagnostics = {
        "duration_s": duration_s,
        "dt_s": dt,
        "n_steps": float(n_steps),
        "minimum_free_state": min_free,
        "minimum_bound_state": min_bound,
        "max_relative_receptor_mass_error": float(
            abs(final_relative_count - expected_relative_count) / expected_relative_count
        ),
    }
    return bound_state, diagnostics


def analytic_equilibrium_check(p: Parameters, pattern: str = "patched") -> dict[str, Any]:
    ligand_nM = 50.0
    integrated_bound, diagnostics = run_constant_ligand_limit(p, pattern, ligand_nM)
    relative_receptor = receptor_pattern(p, pattern)
    equilibrium_fraction = relative_receptor * (
        p.kon_s1_nM1 * ligand_nM / (p.koff_s1 + p.kon_s1_nM1 * ligand_nM)
    )
    relative_error = np.abs(integrated_bound - equilibrium_fraction)
    absolute_error = p.mean_receptor_density * relative_error
    return {
        "pattern": pattern,
        "ligand_nM": ligand_nM,
        "kd_from_rates_nM": p.koff_s1 / p.kon_s1_nM1,
        "max_absolute_error_receptor_um2": float(np.max(absolute_error)),
        "max_relative_error": float(np.max(relative_error / relative_receptor)),
        "mean_occupancy": float(np.mean(equilibrium_fraction / relative_receptor)),
        "zero_diffusion_limit": True,
        **diagnostics,
    }


def sensitivity_analysis(p: Parameters) -> list[dict[str, Any]]:
    controls = {
        "ligand_diffusion_um2_s": (
            0.5 * p.ligand_diffusion_um2_s,
            1.5 * p.ligand_diffusion_um2_s,
        ),
        "koff_s1": (0.5 * p.koff_s1, 1.5 * p.koff_s1),
        "patch_sigma_um": (0.5 * p.patch_sigma_um, 1.5 * p.patch_sigma_um),
    }
    rows: list[dict[str, Any]] = []
    for parameter_name, values in controls.items():
        for factor, value in ((0.5, values[0]), (1.5, values[1])):
            altered = replace(p, **{parameter_name: value})
            comparison = compare_patterns(altered)
            rows.append(
                {
                    "parameter": parameter_name,
                    "factor": factor,
                    "value": value,
                    "peak_ratio": comparison["peak_ratio"],
                    "cv_increment": comparison["cv_increment"],
                    "total_response_ratio_at_pulse_end": comparison[
                        "total_response_ratio_at_pulse_end"
                    ],
                    "activation_integral_ratio": comparison["activation_integral_ratio"],
                    "uniform_peak_time_s": comparison["uniform"]["pulse_peak_time_s"],
                    "patched_peak_time_s": comparison["patched"]["pulse_peak_time_s"],
                    "minimum_free_state": min(
                        comparison["uniform"]["minimum_free_state"],
                        comparison["patched"]["minimum_free_state"],
                    ),
                    "minimum_bound_state": min(
                        comparison["uniform"]["minimum_bound_state"],
                        comparison["patched"]["minimum_bound_state"],
                    ),
                    "minimum_ligand_nM": min(
                        comparison["uniform"]["minimum_ligand_nM"],
                        comparison["patched"]["minimum_ligand_nM"],
                    ),
                    "max_relative_receptor_mass_error": max(
                        comparison["uniform"]["max_relative_receptor_mass_error"],
                        comparison["patched"]["max_relative_receptor_mass_error"],
                    ),
                }
            )
    return rows


def source_catalog() -> list[dict[str, Any]]:
    return [
        {
            "citation": "Briddon et al., PNAS 101, 4673-4678 (2004)",
            "doi": "10.1073/pnas.0400420101",
            "location": "Fig. 6b and accompanying text",
            "value": "Kd = 33 nM; Bmax = 75 nM, approximately 55 receptors/um^2",
            "unit": "nM; receptor um^-2",
            "use": "density and affinity anchors",
        },
        {
            "citation": "Suzuki et al., Biophysical Journal 88, 3659-3680 (2005)",
            "doi": "10.1529/biophysj.104.048538",
            "location": "Fig. 4d and Conclusions",
            "value": "L = 210 nm; residency = 45 ms; D(25 us,100 ms) = 0.23 um^2/s",
            "unit": "nm; ms; um^2/s",
            "use": "membrane patch and mobility context",
        },
        {
            "citation": "Calebiro et al., PNAS 110, 743-748 (2013)",
            "doi": "10.1073/pnas.1205798110",
            "location": "Fig. 2D/I",
            "value": "beta1AR median D = 0.052 um^2/s; beta2AR median D = 0.039 um^2/s",
            "unit": "um^2/s",
            "use": "receptor-state mobility anchor",
        },
    ]


def run_study() -> dict[str, Any]:
    p = Parameters()
    prereg_path = Path(__file__).resolve().with_name("PREREG.md")
    checksum_path = Path(__file__).resolve().with_name("PREREG.sha256")
    prereg_sha256 = hashlib.sha256(prereg_path.read_bytes()).hexdigest()
    expected_prereg_sha256 = checksum_path.read_text(encoding="utf-8").split()[0]
    primary = compare_patterns(p)
    zero_pulse = run_simulation(
        replace(p, pulse_duration_s=0.0, duration_s=0.02), "patched"
    )
    sensitivity = sensitivity_analysis(p)
    analytic = analytic_equilibrium_check(p)
    dt, n_steps, pulse_steps = validate_parameters(p)
    all_simulations = [primary["uniform"], primary["patched"], zero_pulse]
    mass_error = max(
        [run["max_relative_receptor_mass_error"] for run in all_simulations]
        + [row["max_relative_receptor_mass_error"] for row in sensitivity]
    )
    minimum_state = min(
        [
            min(run["minimum_free_state"], run["minimum_bound_state"])
            for run in all_simulations
        ]
        + [
            min(row["minimum_free_state"], row["minimum_bound_state"])
            for row in sensitivity
        ]
    )
    minimum_ligand = min(
        [run["minimum_ligand_nM"] for run in all_simulations]
        + [row["minimum_ligand_nM"] for row in sensitivity]
    )
    primary_values = [
        primary["peak_ratio"],
        primary["cv_increment"],
        primary["total_response_ratio_at_pulse_end"],
        primary["activation_integral_ratio"],
    ]
    sensitivity_values = [
        value
        for row in sensitivity
        for value in (
            row["peak_ratio"],
            row["cv_increment"],
            row["total_response_ratio_at_pulse_end"],
            row["activation_integral_ratio"],
        )
    ]
    validity_checks = {
        "prereg_checksum_matches": prereg_sha256 == expected_prereg_sha256,
        "receptor_mass_error_within_1e-10": mass_error <= 1e-10,
        "normalized_states_nonnegative_within_1e-10": minimum_state >= -1e-10,
        "ligand_nonnegative_within_1e-10_nM": minimum_ligand >= -1e-10,
        "analytic_limit_error_within_1e-10": (
            analytic["max_absolute_error_receptor_um2"] <= 1e-10
            and analytic["max_relative_error"] <= 1e-10
        ),
        "zero_pulse_activation_zero": (
            max(
                zero_pulse["pulse_peak"]["peak_activation"],
                zero_pulse["final"]["peak_activation"],
            )
            == 0.0
        ),
        "unit_dimensions_consistent": unit_check(p)["consistent"],
        "primary_and_sensitivity_values_finite": bool(
            np.all(np.isfinite(primary_values + sensitivity_values))
        ),
    }
    validity_pass = all(validity_checks.values())
    if primary["primary_pass"]:
        primary_result = "mechanistic_prediction_supported"
    elif primary["peak_ratio"] < 1.05:
        primary_result = "failed_directional_prediction"
    else:
        primary_result = "inconclusive_small_effect"
    return {
        "id": "BT-HX-Q080",
        "status": "completed" if validity_pass else "invalid",
        "prereg_sha256": prereg_sha256,
        "model": "deterministic effective 2D no-flux reaction-diffusion model",
        "source_catalog": source_catalog(),
        "parameters": parameter_table(p),
        "unit_check": unit_check(p),
        "derived": {
            "domain_area_um2": p.domain_length_um**2,
            "cell_area_um2": (p.domain_length_um / p.grid_points) ** 2,
            "total_receptor_count": p.mean_receptor_density * p.domain_length_um**2,
            "diffusion_cfl": 4.0
            * max(p.ligand_diffusion_um2_s, p.receptor_diffusion_um2_s)
            * dt
            / (p.domain_length_um / p.grid_points) ** 2,
            "default_steps": n_steps,
            "default_pulse_steps": pulse_steps,
            "kd_from_rates_nM": p.koff_s1 / p.kon_s1_nM1,
            "source_integral_nM_um2": primary["patched"]["source_integral_nM_um2"],
        },
        "primary": primary,
        "sensitivity": sensitivity,
        "controls": {
            "zero_pulse": {
                "maximum_bound_state": zero_pulse["pulse_peak"]["peak_activation"],
                "maximum_activation": max(
                    zero_pulse["pulse_peak"]["peak_activation"],
                    zero_pulse["final"]["peak_activation"],
                ),
                "receptor_mass_error": zero_pulse["max_relative_receptor_mass_error"],
            },
            "analytic_equilibrium": analytic,
        },
        "validity": {
            "all_checks_pass": validity_pass,
            "checks": validity_checks,
            "maximum_receptor_mass_error": mass_error,
            "minimum_normalized_receptor_state": minimum_state,
            "minimum_ligand_nM": minimum_ligand,
        },
        "interpretation": {
            "primary_result": primary_result,
            "no_measured_data": True,
            "deterministic_not_uncertainty_interval": True,
            "patch_width_and_pulse_are_assumptions": True,
        },
        "limitations": [
            "The ligand field is an effective 2D nM field, not a resolved extracellular concentration profile.",
            "Receptor-receptor cooperativity, desensitization, endocytosis, and intracellular signaling are omitted.",
            "The output is a deterministic mechanism test and does not estimate biological uncertainty.",
        ],
        "next_resolution_step": "Fit measured same-cell receptor coordinates and a calibrated ligand dose to patch width, effective ligand diffusion, and kinetic rates using preregistered holdout cells.",
    }


def write_study(output_path: str) -> dict[str, Any]:
    study = run_study()
    path = Path(output_path)
    path.write_text(json.dumps(study, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return study


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    study = run_study()
    if args.output is None:
        print(json.dumps(study, indent=2, sort_keys=True))
    else:
        Path(args.output).write_text(json.dumps(study, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
