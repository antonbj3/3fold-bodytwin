from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable

import numpy as np


@dataclass(frozen=True)
class Layer:
    name: str
    thickness_m: float
    modulus_Pa: float
    permeability_m4_Ns: float
    provenance: str


@dataclass(frozen=True)
class Config:
    viscosity_Pa_s: float = 1.0e-3  # Compatibility metadata; supplied permeability is hydraulic.
    alpha: float = 1.0
    film_thickness_m: float = 1.0e-6
    film_modulus_Pa: float = 1.0e6
    film_conductance_m_per_Pa_s: float = 1.0e-8
    direct_conductance_m_per_Pa_s: float = 1.0e-6
    bone_conductance_m_per_Pa_s: float = 1.0e-6
    normalization_strain: float = 0.05
    reference_modulus_Pa: float = 0.73e6
    n_per_layer: int = 16
    dt_max_s: float = 0.05
    observation_times_s: tuple[float, ...] = (1.0, 10.0, 60.0)


@dataclass(frozen=True)
class Grid:
    faces_m: np.ndarray
    centers_m: np.ndarray
    widths_m: np.ndarray
    modulus_Pa: np.ndarray
    storage_Pa_inv: np.ndarray
    conductivity_m2_Pa_s: np.ndarray
    layer_index: np.ndarray
    layer_names: tuple[str, ...]
    thicknesses_m: tuple[float, ...]


REFERENCE_STRESS_PA = 0.05 * 0.73e6

BASE_LAYERS: tuple[Layer, ...] = (
    Layer(
        name="superficial",
        thickness_m=0.20e-3,
        modulus_Pa=0.32e6,
        permeability_m4_Ns=0.39e-15,
        provenance="Krishnan et al. 2003, Table 1, femur layer 4; thickness frozen assumption",
    ),
    Layer(
        name="bulk",
        thickness_m=1.80e-3,
        modulus_Pa=0.73e6,
        permeability_m4_Ns=0.58e-15,
        provenance="Krishnan et al. 2003, Table 1, femur layer 1; thickness frozen assumption",
    ),
    Layer(
        name="calcified_barrier",
        thickness_m=0.10e-3,
        modulus_Pa=500.0e6,
        permeability_m4_Ns=0.05e-15,
        provenance="frozen mechanistic assumption; no measured value claimed",
    ),
)

SOURCE_REFERENCES = {
    "park_2003": {
        "citation": "Park et al., Journal of Biomechanics 36, 1785-1796 (2003)",
        "doi": "10.1016/S0021-9290(03)00231-8",
        "location": "Table 1",
        "value": "bovine articular surface 94 +/- 4% and deep zone 71 +/- 8% fluid load support",
        "unit": "dimensionless",
    },
    "krishnan_2003": {
        "citation": "Krishnan et al., Journal of Biomechanical Engineering 125, 569-577 (2003)",
        "doi": "10.1115/1.1610018",
        "location": "Tables 1-2",
        "value": "femur H-A 0.32-0.73 MPa, kz 0.39-0.58e-15 m^4/(N s), contact Wp/W 0.93-0.98",
        "unit": "MPa, m^4/(N s), dimensionless",
    },
    "flannery_1999": {
        "citation": "Flannery et al., Biochemical and Biophysical Research Communications 254, 535-541 (1999)",
        "doi": "10.1006/bbrc.1998.0104",
        "location": "article abstract/title",
        "value": "surface-associated lubricating protein context",
        "unit": "not used as a numeric parameter",
    },
}


def parameter_table(config: Config) -> list[dict[str, Any]]:
    rows = [
        {
            "name": "superficial.thickness",
            "value": BASE_LAYERS[0].thickness_m,
            "unit": "m",
            "provenance": BASE_LAYERS[0].provenance,
        },
        {
            "name": "superficial.modulus",
            "value": BASE_LAYERS[0].modulus_Pa,
            "unit": "Pa",
            "provenance": BASE_LAYERS[0].provenance,
        },
        {
            "name": "superficial.permeability",
            "value": BASE_LAYERS[0].permeability_m4_Ns,
            "unit": "m^4/(N s)",
            "provenance": BASE_LAYERS[0].provenance,
        },
        {
            "name": "bulk.thickness",
            "value": BASE_LAYERS[1].thickness_m,
            "unit": "m",
            "provenance": BASE_LAYERS[1].provenance,
        },
        {
            "name": "bulk.modulus",
            "value": BASE_LAYERS[1].modulus_Pa,
            "unit": "Pa",
            "provenance": BASE_LAYERS[1].provenance,
        },
        {
            "name": "bulk.permeability",
            "value": BASE_LAYERS[1].permeability_m4_Ns,
            "unit": "m^4/(N s)",
            "provenance": BASE_LAYERS[1].provenance,
        },
        {
            "name": "calcified_barrier.thickness",
            "value": BASE_LAYERS[2].thickness_m,
            "unit": "m",
            "provenance": BASE_LAYERS[2].provenance,
        },
        {
            "name": "calcified_barrier.modulus",
            "value": BASE_LAYERS[2].modulus_Pa,
            "unit": "Pa",
            "provenance": BASE_LAYERS[2].provenance,
        },
        {
            "name": "calcified_barrier.permeability",
            "value": BASE_LAYERS[2].permeability_m4_Ns,
            "unit": "m^4/(N s)",
            "provenance": BASE_LAYERS[2].provenance,
        },
        {
            "name": "film.thickness",
            "value": config.film_thickness_m,
            "unit": "m",
            "provenance": "frozen mechanistic assumption",
        },
        {
            "name": "film.modulus",
            "value": config.film_modulus_Pa,
            "unit": "Pa",
            "provenance": "frozen mechanistic assumption",
        },
        {
            "name": "film.hydraulic_conductance",
            "value": config.film_conductance_m_per_Pa_s,
            "unit": "m/(Pa s)",
            "provenance": "frozen mechanistic assumption",
        },
        {
            "name": "direct_drain_conductance",
            "value": config.direct_conductance_m_per_Pa_s,
            "unit": "m/(Pa s)",
            "provenance": "frozen mechanistic assumption",
        },
        {
            "name": "bone_boundary_conductance",
            "value": config.bone_conductance_m_per_Pa_s,
            "unit": "m/(Pa s)",
            "provenance": "frozen mechanistic assumption",
        },
        {
            "name": "viscosity",
            "value": config.viscosity_Pa_s,
            "unit": "Pa s",
            "provenance": "frozen water-like first-order assumption",
        },
        {
            "name": "alpha",
            "value": config.alpha,
            "unit": "dimensionless",
            "provenance": "frozen reduced-storage assumption",
        },
        {
            "name": "normalization_strain",
            "value": config.normalization_strain,
            "unit": "dimensionless",
            "provenance": "frozen load normalization; not patient-specific",
        },
        {
            "name": "reference_modulus",
            "value": config.reference_modulus_Pa,
            "unit": "Pa",
            "provenance": "bulk anchor from Krishnan et al. 2003, Table 1",
        },
    ]
    return rows


def validate_parameters(layers: Iterable[Layer], config: Config) -> dict[str, Any]:
    layer_list = list(layers)
    checks: dict[str, bool] = {
        "positive_thickness": all(layer.thickness_m > 0.0 for layer in layer_list),
        "positive_modulus": all(layer.modulus_Pa > 0.0 for layer in layer_list),
        "nonnegative_permeability": all(layer.permeability_m4_Ns >= 0.0 for layer in layer_list),
        "positive_viscosity": config.viscosity_Pa_s > 0.0,
        "positive_alpha": config.alpha > 0.0,
        "positive_conductances": min(
            config.film_conductance_m_per_Pa_s,
            config.direct_conductance_m_per_Pa_s,
            config.bone_conductance_m_per_Pa_s,
        )
        >= 0.0,
        "positive_cell_count": config.n_per_layer >= 1,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError(f"parameter validation failed: {failed}")
    conductivity = [layer.permeability_m4_Ns for layer in layer_list]
    storage = [config.alpha * config.alpha / layer.modulus_Pa for layer in layer_list]
    return {
        "status": "pass",
        "checks": checks,
        "derived_conductivity_m2_Pa_s": conductivity,
        "derived_storage_Pa_inv": storage,
        "dimension_relations": {
            "Darcy_flux": "K [m^2/(Pa s)] * grad(p) [Pa/m] = m/s",
            "storage_rate": "S [1/Pa] * dp/dt [Pa/s] = 1/s",
            "boundary_flux": "G [m/(Pa s)] * p [Pa] = m/s",
        },
    }


def build_layers(scenario: str, config: Config) -> tuple[Layer, ...]:
    base = list(BASE_LAYERS)
    if scenario == "full" or scenario == "no_film":
        return tuple(base)
    if scenario == "no_barrier":
        return tuple(base[:2])
    if scenario == "merged_surface":
        bulk = base[1]
        return (
            Layer(
                name="bulk_merged",
                thickness_m=bulk.thickness_m + base[0].thickness_m,
                modulus_Pa=bulk.modulus_Pa,
                permeability_m4_Ns=bulk.permeability_m4_Ns,
                provenance="surface zone merged into bulk; bulk values from Krishnan et al. 2003 Table 1",
            ),
            base[2],
        )
    if scenario == "homogeneous":
        total_thickness = sum(layer.thickness_m for layer in base)
        mean_modulus = sum(layer.thickness_m * layer.modulus_Pa for layer in base) / total_thickness
        mean_permeability = sum(
            layer.thickness_m * layer.permeability_m4_Ns for layer in base
        ) / total_thickness
        return (
            Layer(
                name="homogeneous",
                thickness_m=total_thickness,
                modulus_Pa=mean_modulus,
                permeability_m4_Ns=mean_permeability,
                provenance="thickness-weighted arithmetic placebo from frozen layers",
            ),
        )
    raise ValueError(f"unknown scenario: {scenario}")


def apply_overrides(layers: Iterable[Layer], overrides: dict[str, float]) -> tuple[Layer, ...]:
    adjusted: list[Layer] = []
    for layer in layers:
        modulus = layer.modulus_Pa
        permeability = layer.permeability_m4_Ns
        if layer.name == "superficial":
            modulus = float(overrides.get("M_superficial", modulus))
        if layer.name in {"bulk", "bulk_merged"}:
            permeability = float(overrides.get("k_bulk", permeability))
        adjusted.append(
            replace(
                layer,
                modulus_Pa=modulus,
                permeability_m4_Ns=permeability,
            )
        )
    return tuple(adjusted)


def build_grid(layers: Iterable[Layer], config: Config, n_per_layer: int | None = None) -> Grid:
    layer_list = list(layers)
    count = config.n_per_layer if n_per_layer is None else n_per_layer
    if count < 1:
        raise ValueError("n_per_layer must be positive")
    faces = [0.0]
    current = 0.0
    for layer in layer_list:
        segment = np.linspace(current, current + layer.thickness_m, count + 1)
        faces.extend(segment[1:].tolist())
        current += layer.thickness_m
    faces_array = np.asarray(faces, dtype=float)
    centers = 0.5 * (faces_array[:-1] + faces_array[1:])
    widths = np.diff(faces_array)
    modulus = np.asarray([layer.modulus_Pa for layer in layer_list for _ in range(count)], dtype=float)
    permeability = np.asarray(
        [layer.permeability_m4_Ns for layer in layer_list for _ in range(count)],
        dtype=float,
    )
    storage = config.alpha * config.alpha / modulus
    conductivity = permeability  # m4/(N s) already equals m2/(Pa s)
    layer_index = np.asarray(
        [index for index, layer in enumerate(layer_list) for _ in range(count)],
        dtype=int,
    )
    return Grid(
        faces_m=faces_array,
        centers_m=centers,
        widths_m=widths,
        modulus_Pa=modulus,
        storage_Pa_inv=storage,
        conductivity_m2_Pa_s=conductivity,
        layer_index=layer_index,
        layer_names=tuple(layer.name for layer in layer_list),
        thicknesses_m=tuple(layer.thickness_m for layer in layer_list),
    )


def interface_conductivity(grid: Grid) -> np.ndarray:
    if len(grid.centers_m) == 1:
        return np.asarray([], dtype=float)
    left = grid.conductivity_m2_Pa_s[:-1]
    right = grid.conductivity_m2_Pa_s[1:]
    harmonic = 2.0 * left * right / np.maximum(left + right, np.finfo(float).tiny)
    return harmonic


def boundary_fluxes(
    pressure_Pa: np.ndarray,
    grid: Grid,
    top_conductance_m_per_Pa_s: float,
    bone_conductance_m_per_Pa_s: float,
) -> np.ndarray:
    flux = np.zeros(len(pressure_Pa) + 1, dtype=float)
    flux[0] = -top_conductance_m_per_Pa_s * pressure_Pa[0]
    flux[-1] = bone_conductance_m_per_Pa_s * pressure_Pa[-1]
    face_k = interface_conductivity(grid)
    for index in range(len(pressure_Pa) - 1):
        distance = grid.centers_m[index + 1] - grid.centers_m[index]
        flux[index + 1] = -face_k[index] * (
            pressure_Pa[index + 1] - pressure_Pa[index]
        ) / distance
    return flux


def pressure_derivative(
    pressure_Pa: np.ndarray,
    grid: Grid,
    top_conductance_m_per_Pa_s: float,
    bone_conductance_m_per_Pa_s: float,
) -> np.ndarray:
    flux = boundary_fluxes(
        pressure_Pa,
        grid,
        top_conductance_m_per_Pa_s,
        bone_conductance_m_per_Pa_s,
    )
    return -(flux[1:] - flux[:-1]) / grid.widths_m / grid.storage_Pa_inv


def implicit_pressure_step(
    pressure_Pa: np.ndarray,
    timestep_s: float,
    grid: Grid,
    top_conductance_m_per_Pa_s: float,
    bone_conductance_m_per_Pa_s: float,
) -> np.ndarray:
    count = len(pressure_Pa)
    face_k = interface_conductivity(grid)
    left_conductance = np.zeros(count, dtype=float)
    right_conductance = np.zeros(count, dtype=float)
    left_conductance[0] = top_conductance_m_per_Pa_s
    right_conductance[-1] = bone_conductance_m_per_Pa_s
    for index in range(count - 1):
        distance = grid.centers_m[index + 1] - grid.centers_m[index]
        conductance = face_k[index] / distance
        right_conductance[index] = conductance
        left_conductance[index + 1] = conductance
    storage_rate = grid.storage_Pa_inv / timestep_s
    diagonal = storage_rate + (left_conductance + right_conductance) / grid.widths_m
    lower = np.zeros(count - 1, dtype=float)
    upper = np.zeros(count - 1, dtype=float)
    for index in range(count - 1):
        lower[index] = -right_conductance[index] / grid.widths_m[index + 1]
        upper[index] = -right_conductance[index] / grid.widths_m[index]
    modified_diagonal = diagonal.copy()
    modified_upper = upper.copy()
    rhs = storage_rate * pressure_Pa
    for index in range(1, count):
        factor = lower[index - 1] / modified_diagonal[index - 1]
        modified_diagonal[index] -= factor * modified_upper[index - 1]
        rhs[index] -= factor * rhs[index - 1]
    solution = np.zeros(count, dtype=float)
    solution[-1] = rhs[-1] / modified_diagonal[-1]
    for index in range(count - 2, -1, -1):
        solution[index] = (
            rhs[index] - modified_upper[index] * solution[index + 1]
        ) / modified_diagonal[index]
    return solution


def stable_timestep(
    grid: Grid,
    top_conductance_m_per_Pa_s: float,
    bone_conductance_m_per_Pa_s: float,
    dt_max_s: float,
) -> float:
    return float(dt_max_s)


def fluid_content_m(pressure_Pa: np.ndarray, grid: Grid) -> float:
    return float(np.sum(grid.storage_Pa_inv * pressure_Pa * grid.widths_m))


def layer_means(pressure_Pa: np.ndarray, grid: Grid) -> dict[str, float]:
    result: dict[str, float] = {}
    for index, name in enumerate(grid.layer_names):
        mask = grid.layer_index == index
        result[name] = float(np.average(pressure_Pa[mask], weights=grid.widths_m[mask]))
    return result


def snapshot_metrics(
    pressure_Pa: np.ndarray,
    flux: np.ndarray,
    grid: Grid,
    total_stress_Pa: float,
    cumulative_outflow_m: float,
    initial_content_m: float,
    film_present: bool,
    film_thickness_m: float,
    film_modulus_Pa: float,
) -> dict[str, Any]:
    surface_pressure = float(pressure_Pa[0])
    bone_pressure = float(pressure_Pa[-1])
    surface_solid_stress = float(total_stress_Pa - surface_pressure)
    tissue_strain = float(
        np.sum(
            ((total_stress_Pa - pressure_Pa) / grid.modulus_Pa) * grid.widths_m
        )
        / sum(grid.thicknesses_m)
    )
    film_strain = total_stress_Pa / film_modulus_Pa if film_present else 0.0
    film_compliance = film_thickness_m / film_modulus_Pa if film_present else 0.0
    content = fluid_content_m(pressure_Pa, grid)
    return {
        "surface_pressure_Pa": surface_pressure,
        "bone_pressure_Pa": bone_pressure,
        "fluid_support_surface": surface_pressure / total_stress_Pa,
        "fluid_support_bone": bone_pressure / total_stress_Pa,
        "solid_stress_surface_Pa": surface_solid_stress,
        "mean_tissue_solid_strain": tissue_strain,
        "film_solid_strain": film_strain,
        "film_series_compliance_m_per_Pa": film_compliance,
        "top_exudation_flux_m_s": float(-flux[0]),
        "bone_exchange_flux_m_s": float(flux[-1]),
        "net_outward_flux_m_s": float(flux[-1] - flux[0]),
        "fluid_content_m": content,
        "cumulative_outflow_m": cumulative_outflow_m,
        "mass_balance_residual_m": initial_content_m - cumulative_outflow_m - content,
        "layer_pressure_mean_Pa": layer_means(pressure_Pa, grid),
    }


def simulate(
    layers: Iterable[Layer],
    config: Config,
    film_present: bool,
    top_conductance_m_per_Pa_s: float,
    bone_conductance_m_per_Pa_s: float | None = None,
    n_per_layer: int | None = None,
    observation_times_s: tuple[float, ...] | None = None,
) -> dict[str, Any]:
    layer_list = list(layers)
    if not layer_list:
        raise ValueError("at least one layer is required")
    grid = build_grid(layer_list, config, n_per_layer=n_per_layer)
    validate_parameters(layer_list, config)
    if bone_conductance_m_per_Pa_s is None:
        bone_conductance_m_per_Pa_s = config.bone_conductance_m_per_Pa_s
    times = tuple(config.observation_times_s if observation_times_s is None else observation_times_s)
    if any(time <= 0.0 for time in times):
        raise ValueError("observation times must be positive")
    if tuple(sorted(times)) != times:
        raise ValueError("observation times must be sorted")
    total_stress_Pa = config.normalization_strain * config.reference_modulus_Pa
    pressure = np.full(len(grid.centers_m), total_stress_Pa, dtype=float)
    initial_flux = boundary_fluxes(
        pressure,
        grid,
        top_conductance_m_per_Pa_s,
        bone_conductance_m_per_Pa_s,
    )
    initial_content = fluid_content_m(pressure, grid)
    cumulative_outflow = 0.0
    cumulative_top = 0.0
    cumulative_bone = 0.0
    snapshots: list[dict[str, Any]] = []
    target_index = 0
    time_s = 0.0
    initial_metric = snapshot_metrics(
        pressure,
        initial_flux,
        grid,
        total_stress_Pa,
        cumulative_outflow,
        initial_content,
        film_present,
        config.film_thickness_m,
        config.film_modulus_Pa,
    )
    initial_metric["cumulative_top_exudation_m"] = 0.0
    initial_metric["cumulative_bone_exchange_m"] = 0.0
    initial_metric["time_s"] = 0.0
    snapshots.append(initial_metric)
    while target_index < len(times):
        target_time = times[target_index]
        while time_s < target_time - 1.0e-12:
            timestep = min(
                stable_timestep(
                    grid,
                    top_conductance_m_per_Pa_s,
                    bone_conductance_m_per_Pa_s,
                    config.dt_max_s,
                ),
                target_time - time_s,
            )
            old_flux = boundary_fluxes(
                pressure,
                grid,
                top_conductance_m_per_Pa_s,
                bone_conductance_m_per_Pa_s,
            )
            pressure = implicit_pressure_step(
                pressure,
                timestep,
                grid,
                top_conductance_m_per_Pa_s,
                bone_conductance_m_per_Pa_s,
            )
            new_flux = boundary_fluxes(
                pressure,
                grid,
                top_conductance_m_per_Pa_s,
                bone_conductance_m_per_Pa_s,
            )
            new_net = new_flux[-1] - new_flux[0]
            cumulative_outflow += new_net * timestep
            cumulative_top += -new_flux[0] * timestep
            cumulative_bone += new_flux[-1] * timestep
            time_s += timestep
        metric = snapshot_metrics(
            pressure,
            boundary_fluxes(
                pressure,
                grid,
                top_conductance_m_per_Pa_s,
                bone_conductance_m_per_Pa_s,
            ),
            grid,
            total_stress_Pa,
            cumulative_outflow,
            initial_content,
            film_present,
            config.film_thickness_m,
            config.film_modulus_Pa,
        )
        metric["cumulative_top_exudation_m"] = cumulative_top
        metric["cumulative_bone_exchange_m"] = cumulative_bone
        metric["time_s"] = float(target_time)
        snapshots.append(metric)
        target_index += 1
    half_time = None
    half_pressure = 0.5 * total_stress_Pa
    for previous, current in zip(snapshots[:-1], snapshots[1:]):
        p0 = previous["surface_pressure_Pa"]
        p1 = current["surface_pressure_Pa"]
        if p1 <= half_pressure <= p0:
            if p0 == p1:
                half_time = float(current["time_s"])
            else:
                fraction = (half_pressure - p0) / (p1 - p0)
                half_time = float(previous["time_s"] + fraction * (current["time_s"] - previous["time_s"]))
            break
    final_metric = snapshots[-1]
    max_mass_residual = max(
        abs(metric["mass_balance_residual_m"]) for metric in snapshots
    )
    return {
        "snapshots": snapshots,
        "drainage_half_time_s": half_time,
        "max_mass_balance_residual_m": max_mass_residual,
        "normalized_mass_balance_error": max_mass_residual / max(initial_content, np.finfo(float).tiny),
        "total_stress_Pa": total_stress_Pa,
        "final_surface_pressure_Pa": final_metric["surface_pressure_Pa"],
        "final_top_exudation_flux_m_s": final_metric["top_exudation_flux_m_s"],
        "final_bone_exchange_flux_m_s": final_metric["bone_exchange_flux_m_s"],
        "n_cells": len(grid.centers_m),
        "n_per_layer": config.n_per_layer if n_per_layer is None else n_per_layer,
        "dt_max_s": config.dt_max_s,
    }


def scenario_configuration(
    scenario: str,
    config: Config,
    overrides: dict[str, float] | None = None,
) -> tuple[tuple[Layer, ...], bool, float, float]:
    overrides = {} if overrides is None else overrides
    layers = build_layers(scenario, config)
    if scenario == "no_film" or scenario == "homogeneous":
        film_present = False
        top_conductance = config.direct_conductance_m_per_Pa_s
    else:
        film_present = True
        top_conductance = config.film_conductance_m_per_Pa_s
    if "G_top" in overrides:
        if film_present:
            top_conductance = float(overrides["G_top"])
    layers = apply_overrides(layers, overrides)
    return layers, film_present, top_conductance, config.bone_conductance_m_per_Pa_s


def compact_scenario(
    scenario: str,
    config: Config,
    overrides: dict[str, float] | None = None,
    n_per_layer: int | None = None,
) -> dict[str, Any]:
    layers, film_present, top_conductance, bone_conductance = scenario_configuration(
        scenario,
        config,
        overrides,
    )
    result = simulate(
        layers,
        config,
        film_present,
        top_conductance,
        bone_conductance_m_per_Pa_s=bone_conductance,
        n_per_layer=n_per_layer,
    )
    return {
        "scenario": scenario,
        "film_present": film_present,
        "top_conductance_m_per_Pa_s": top_conductance,
        "bone_conductance_m_per_Pa_s": bone_conductance,
        "layers": [
            {
                "name": layer.name,
                "thickness_m": layer.thickness_m,
                "modulus_Pa": layer.modulus_Pa,
                "permeability_m4_Ns": layer.permeability_m4_Ns,
                "provenance": layer.provenance,
            }
            for layer in layers
        ],
        "total_stress_Pa": result["total_stress_Pa"],
        "drainage_half_time_s": result["drainage_half_time_s"],
        "max_mass_balance_error": result["normalized_mass_balance_error"],
        "n_cells": result["n_cells"],
        "n_per_layer": result["n_per_layer"],
        "snapshots": result["snapshots"],
    }


def metric_at(result: dict[str, Any], time_s: float) -> dict[str, Any]:
    for snapshot in result["snapshots"]:
        if abs(snapshot["time_s"] - time_s) < 1.0e-9:
            return snapshot
    raise KeyError(f"missing snapshot at {time_s}")


def layer_effects(results: dict[str, dict[str, Any]]) -> dict[str, float]:
    full_one = metric_at(results["full"], 1.0)
    merged_one = metric_at(results["merged_surface"], 1.0)
    nofilm_one = metric_at(results["no_film"], 1.0)
    merged_ten = metric_at(results["merged_surface"], 10.0)
    nobar_one = metric_at(results["no_barrier"], 10.0)
    full_ten = metric_at(results["full"], 10.0)
    full_bone = full_ten["bone_exchange_flux_m_s"]
    nofilm_top = nofilm_one["top_exudation_flux_m_s"]
    full_top = full_one["top_exudation_flux_m_s"]
    scale_top = max(abs(full_top), 1.0e-30)
    scale_bone = max(abs(full_bone), 1.0e-30)
    return {
        "superficial_fluid_support_delta_1s": full_one["fluid_support_surface"]
        - merged_one["fluid_support_surface"],
        "superficial_bone_flux_relative_delta_10s": (
            full_ten["bone_exchange_flux_m_s"] - merged_ten["bone_exchange_flux_m_s"]
        )
        / scale_bone,
        "film_fluid_support_delta_1s": full_one["fluid_support_surface"]
        - nofilm_one["fluid_support_surface"],
        "film_top_flux_relative_delta_1s": (full_top - nofilm_top) / scale_top,
        "barrier_bone_flux_relative_delta_10s": (
            full_ten["bone_exchange_flux_m_s"] - nobar_one["bone_exchange_flux_m_s"]
        )
        / scale_bone,
        "barrier_fluid_support_delta_10s": full_ten["fluid_support_surface"]
        - metric_at(results["no_barrier"], 10.0)["fluid_support_surface"],
    }


def sensitivity_effects(config: Config, overrides: dict[str, float]) -> dict[str, float]:
    results = {
        name: compact_scenario(name, config, overrides=overrides)
        for name in ("full", "merged_surface", "no_film", "no_barrier")
    }
    return layer_effects(results)


def sign_consistent(values: Iterable[float], reference: float) -> bool:
    reference_sign = np.sign(reference)
    return bool(reference_sign != 0.0 and all(np.sign(value) == reference_sign for value in values))


def analytical_limit(config: Config) -> dict[str, Any]:
    layer = Layer(
        name="sealed_homogeneous",
        thickness_m=1.0e-3,
        modulus_Pa=1.0e6,
        permeability_m4_Ns=1.0e-15,
        provenance="analytical test fixture",
    )
    result = simulate(
        (layer,),
        config,
        film_present=False,
        top_conductance_m_per_Pa_s=0.0,
        bone_conductance_m_per_Pa_s=0.0,
        n_per_layer=8,
        observation_times_s=(1.0,),
    )
    snapshot = metric_at(result, 1.0)
    expected_pressure = result["total_stress_Pa"]
    pressure_error = abs(snapshot["surface_pressure_Pa"] - expected_pressure) / expected_pressure
    return {
        "expected_surface_pressure_Pa": expected_pressure,
        "numerical_surface_pressure_Pa": snapshot["surface_pressure_Pa"],
        "expected_fluid_support": 1.0,
        "numerical_fluid_support": snapshot["fluid_support_surface"],
        "relative_pressure_error": pressure_error,
        "normalized_mass_balance_error": result["normalized_mass_balance_error"],
        "pass": pressure_error < 0.01 and result["normalized_mass_balance_error"] < 0.01,
    }


def resolution_check(config: Config) -> dict[str, Any]:
    coarse = compact_scenario("full", config, n_per_layer=16)
    fine = compact_scenario("full", config, n_per_layer=32)
    coarse_one = metric_at(coarse, 1.0)
    fine_one = metric_at(fine, 1.0)
    coarse_ten = metric_at(coarse, 10.0)
    fine_ten = metric_at(fine, 10.0)
    psi_error = abs(coarse_one["fluid_support_surface"] - fine_one["fluid_support_surface"])
    top_flux_error = abs(
        coarse_ten["cumulative_top_exudation_m"]
        - fine_ten["cumulative_top_exudation_m"]
    ) / max(abs(fine_ten["cumulative_top_exudation_m"]), 1.0e-30)
    bone_flux_error = abs(
        coarse_ten["cumulative_bone_exchange_m"]
        - fine_ten["cumulative_bone_exchange_m"]
    ) / max(abs(fine_ten["cumulative_bone_exchange_m"]), 1.0e-30)
    instantaneous_top_error = abs(
        coarse_one["top_exudation_flux_m_s"] - fine_one["top_exudation_flux_m_s"]
    ) / max(abs(fine_one["top_exudation_flux_m_s"]), 1.0e-30)
    instantaneous_bone_error = abs(
        coarse_ten["bone_exchange_flux_m_s"] - fine_ten["bone_exchange_flux_m_s"]
    ) / max(abs(fine_ten["bone_exchange_flux_m_s"]), 1.0e-30)
    return {
        "coarse_cells_per_layer": 16,
        "fine_cells_per_layer": 32,
        "fluid_support_absolute_difference_1s": psi_error,
        "cumulative_top_flux_relative_difference_10s": top_flux_error,
        "cumulative_bone_flux_relative_difference_10s": bone_flux_error,
        "instantaneous_top_flux_relative_difference_1s": instantaneous_top_error,
        "instantaneous_bone_flux_relative_difference_10s": instantaneous_bone_error,
        "pass": psi_error < 0.01 and top_flux_error < 0.05 and bone_flux_error < 0.05,
    }


def sensitivity_analysis(config: Config, baseline_effects: dict[str, float]) -> dict[str, Any]:
    parameters = {
        "k_bulk": {
            "base": BASE_LAYERS[1].permeability_m4_Ns,
            "minus50": BASE_LAYERS[1].permeability_m4_Ns * 0.5,
            "plus50": BASE_LAYERS[1].permeability_m4_Ns * 1.5,
        },
        "M_superficial": {
            "base": BASE_LAYERS[0].modulus_Pa,
            "minus50": BASE_LAYERS[0].modulus_Pa * 0.5,
            "plus50": BASE_LAYERS[0].modulus_Pa * 1.5,
        },
        "G_top": {
            "base": config.film_conductance_m_per_Pa_s,
            "minus50": config.film_conductance_m_per_Pa_s * 0.5,
            "plus50": config.film_conductance_m_per_Pa_s * 1.5,
        },
    }
    output: dict[str, Any] = {}
    for name, values in parameters.items():
        cases: dict[str, Any] = {}
        for case in ("minus50", "plus50"):
            effects = sensitivity_effects(config, {name: values[case]})
            cases[case] = effects
        effect_names = (
            "superficial_fluid_support_delta_1s",
            "film_fluid_support_delta_1s",
            "barrier_fluid_support_delta_10s",
        )
        robust = {
            effect_name: sign_consistent(
                [cases["minus50"][effect_name], cases["plus50"][effect_name]],
                baseline_effects[effect_name],
            )
            for effect_name in effect_names
        }
        output[name] = {
            "base_value": values["base"],
            "minus50_value": values["minus50"],
            "plus50_value": values["plus50"],
            "cases": cases,
            "sign_robustness": robust,
        }
    return output


def evaluate_decision(
    config: Config,
    scenarios: dict[str, dict[str, Any]],
    effects: dict[str, float],
    analytical: dict[str, Any],
    resolution: dict[str, Any],
) -> dict[str, Any]:
    full_one = metric_at(scenarios["full"], 1.0)
    full_ten = metric_at(scenarios["full"], 10.0)
    surface_pass = (
        abs(effects["superficial_fluid_support_delta_1s"]) >= 0.02
        or abs(effects["superficial_bone_flux_relative_delta_10s"]) >= 0.25
    )
    film_pass = (
        abs(effects["film_fluid_support_delta_1s"]) >= 0.02
        or abs(effects["film_top_flux_relative_delta_1s"]) >= 0.25
    )
    barrier_pass = (
        abs(effects["barrier_bone_flux_relative_delta_10s"]) >= 0.25
        or abs(effects["barrier_fluid_support_delta_10s"]) >= 0.02
    )
    plausibility_pass = 0.90 <= full_one["fluid_support_surface"] <= 1.0
    mass_pass = scenarios["full"]["max_mass_balance_error"] < 0.01
    core_pass = analytical["pass"] and resolution["pass"] and plausibility_pass and mass_pass
    return {
        "full_fluid_support_1s": full_one["fluid_support_surface"],
        "full_fluid_support_10s": full_ten["fluid_support_surface"],
        "plausibility_anchor_pass": plausibility_pass,
        "mass_balance_pass": mass_pass,
        "superficial_layer_resolved_as_necessary": surface_pass,
        "film_transport_relevant": film_pass,
        "barrier_transport_relevant": barrier_pass,
        "core_mechanistic_test_pass": core_pass,
        "interpretation_rule": "A false layer flag means not resolved as necessary by this model, not biological absence.",
    }


def run_experiment(config: Config | None = None) -> dict[str, Any]:
    active_config = Config() if config is None else config
    analytical = analytical_limit(active_config)
    scenarios = {
        name: compact_scenario(name, active_config)
        for name in ("full", "merged_surface", "no_film", "no_barrier", "homogeneous")
    }
    effects = layer_effects(scenarios)
    resolution = resolution_check(active_config)
    sensitivity = sensitivity_analysis(active_config, effects)
    decision = evaluate_decision(active_config, scenarios, effects, analytical, resolution)
    full_one = metric_at(scenarios["full"], 1.0)
    full_ten = metric_at(scenarios["full"], 10.0)
    return {
        "id": "BT-HX-Q049",
        "model": "one-dimensional layered poroelastic consolidation with interfacial and bone boundary resistances",
        "units": "SI",
        "status": "PASS" if decision["core_mechanistic_test_pass"] else "FAIL",
        "frozen_preregistration_sha256": "19f43a692263df5add62985912991197f7d0746c47dac182354cf988b60fc747",
        "parameter_table": parameter_table(active_config),
        "reference_anchors": SOURCE_REFERENCES,
        "stress_normalization": {
            "total_stress_Pa": scenarios["full"]["total_stress_Pa"],
            "definition": "normalization_strain * reference_modulus",
            "reference_modulus_Pa": active_config.reference_modulus_Pa,
        },
        "scenarios": scenarios,
        "layer_effects": effects,
        "sensitivity": sensitivity,
        "checks": {
            "unit_and_parameter_check": validate_parameters(BASE_LAYERS, active_config),
            "analytical_limit": analytical,
            "resolution": resolution,
        },
        "decision": decision,
        "next_step": "Replace the assumed film conductance and one-dimensional boundary laws with matched contact/transport experiments, then add lateral contact geometry and cyclic compression before assigning biological necessity.",
    }


def json_default(value: Any) -> Any:
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results.json")
    args = parser.parse_args()
    payload = run_experiment()
    output_path = Path(args.output)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": payload["status"], "output": str(output_path)}), flush=True)


if __name__ == "__main__":
    main()
