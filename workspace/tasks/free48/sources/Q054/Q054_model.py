from __future__ import annotations

import argparse
import csv
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np

BASE_DIR = Path(__file__).resolve().parent

REFERENCE = {
    "citation": "Anderson AE, Ellis BJ, Maas SA, Peters CL, Weiss JA. Validation of Finite Element Predictions of Cartilage Contact Pressure in the Human Hip Joint. Journal of Biomechanical Engineering. 2008;130(5):051008.",
    "doi": "10.1115/1.2953472",
    "pmcid": "PMC2840996",
    "location": "Abstract, Results and Figure 4",
    "experimental_peak_pressure_mpa": 10.0,
    "experimental_mean_pressure_mpa_range": [4.4, 5.0],
    "experimental_contact_area_mm2_range": [321.9, 425.1],
    "use": "scale and plausibility reference; not validation of synthetic geometry",
}

EQUATIONS = {
    "reduced_radius": "1/R_star = 1/R_head - 1/R_socket",
    "gap": "g(r) = clearance + r^2/(2 R_star) + form_error_amplitude sin(2 pi r/wavelength)",
    "contact_pressure": "p(r) = layer_stiffness max(travel - g(r), 0)",
    "layer_stiffness": "k = E_contact / layer_thickness",
    "load_balance": "P = integral_0^infinity p(r) 2 pi r dr",
    "analytic_parabolic_limit": "a = (4 R_star P/(pi k))^(1/4), p(0) = k a^2/(2 R_star)",
}


def effective_radius(head_radius_m: float, socket_radius_m: float) -> float:
    if head_radius_m <= 0 or socket_radius_m <= head_radius_m:
        raise ValueError("socket_radius_m must be greater than head_radius_m")
    return head_radius_m * socket_radius_m / (socket_radius_m - head_radius_m)


def layer_properties(shear_modulus_pa: float, poisson_ratio: float, thickness_m: float) -> tuple[float, float, float]:
    if shear_modulus_pa <= 0 or thickness_m <= 0 or not -1 < poisson_ratio < 0.5:
        raise ValueError("invalid layer properties")
    cartilage_modulus_pa = 2.0 * shear_modulus_pa * (1.0 + poisson_ratio)
    contact_modulus_pa = cartilage_modulus_pa / (2.0 * (1.0 - poisson_ratio**2))
    stiffness_pa_per_m = contact_modulus_pa / thickness_m
    return cartilage_modulus_pa, contact_modulus_pa, stiffness_pa_per_m


def default_cases() -> dict[str, dict[str, Any]]:
    common = {
        "load_n": 1800.0,
        "cartilage_shear_modulus_pa": 6.8e6,
        "poisson_ratio": 0.495,
        "layer_thickness_m": 1.55e-3,
        "form_error_wavelength_m": 8.0e-3,
        "n_cells": 1200,
        "radius_max_m": 30.0e-3,
    }
    nominal = {
        **common,
        "name": "nominal",
        "head_radius_m": 25.0e-3,
        "socket_radius_m": 50.0e-3,
        "clearance_m": 0.0,
        "form_error_amplitude_m": 0.0,
    }
    as_built = {
        **common,
        "name": "as_built",
        "head_radius_m": 24.75e-3,
        "socket_radius_m": 50.25e-3,
        "clearance_m": 0.05e-3,
        "form_error_amplitude_m": 0.04e-3,
    }
    return {"nominal": nominal, "as_built": as_built}


def reduced_radius(case: dict[str, Any]) -> float:
    if case.get("reduced_radius_m") is not None:
        value = float(case["reduced_radius_m"])
        if value <= 0:
            raise ValueError("reduced_radius_m must be positive")
        return value
    return effective_radius(float(case["head_radius_m"]), float(case["socket_radius_m"]))


def _grid(case: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    n_cells = int(case["n_cells"])
    radius_max_m = float(case["radius_max_m"])
    if n_cells < 8 or radius_max_m <= 0:
        raise ValueError("invalid radial grid")
    edges = radius_max_m * np.sqrt(np.linspace(0.0, 1.0, n_cells + 1))
    radii = 0.5 * (edges[:-1] + edges[1:])
    areas = math.pi * np.diff(edges**2)
    return radii, areas


def gap_profile(case: dict[str, Any], radii: np.ndarray) -> np.ndarray:
    measured_profile = case.get("measured_gap_profile_m")
    if measured_profile is not None:
        profile = np.asarray(measured_profile, dtype=float)
        if profile.ndim != 2 or profile.shape[1] != 2 or profile.shape[0] < 2:
            raise ValueError("measured_gap_profile_m must have shape (n, 2)")
        if np.any(np.diff(profile[:, 0]) <= 0):
            raise ValueError("measured profile radii must be strictly increasing")
        if np.any(radii < profile[0, 0]) or np.any(radii > profile[-1, 0]):
            raise ValueError("radial grid must lie inside measured profile")
        return np.interp(radii, profile[:, 0], profile[:, 1])
    radius_m = reduced_radius(case)
    clearance_m = float(case["clearance_m"])
    amplitude_m = float(case["form_error_amplitude_m"])
    wavelength_m = float(case["form_error_wavelength_m"])
    if wavelength_m <= 0:
        raise ValueError("form_error_wavelength_m must be positive")
    return clearance_m + radii**2 / (2.0 * radius_m) + amplitude_m * np.sin(2.0 * math.pi * radii / wavelength_m)


def _state(case: dict[str, Any], travel_m: float, target_load_n: float | None = None) -> dict[str, Any]:
    radii, areas = _grid(case)
    gap = gap_profile(case, radii)
    _, _, stiffness = layer_properties(
        float(case["cartilage_shear_modulus_pa"]),
        float(case["poisson_ratio"]),
        float(case["layer_thickness_m"]),
    )
    indentation = travel_m - gap
    pressure = np.where(indentation > 0.0, stiffness * indentation, 0.0)
    active = indentation > 0.0
    load_n = float(np.dot(pressure, areas))
    if load_n <= 0.0:
        raise ValueError("travel does not produce contact")
    contact_area_m2 = float(np.sum(areas[active]))
    equivalent_radius_m = math.sqrt(contact_area_m2 / math.pi)
    weighted_load = pressure * areas
    peak_index = int(np.argmax(pressure))
    inner_radius_m = 0.5 * equivalent_radius_m
    inner_mask = active & (radii <= inner_radius_m)
    outer_mask = active & (radii > inner_radius_m)
    bins = np.linspace(0.0, equivalent_radius_m, 6)
    radial_bins: list[float] = []
    for index in range(5):
        lower = bins[index]
        upper = bins[index + 1]
        if index == 4:
            mask = active & (radii >= lower) & (radii <= upper)
        else:
            mask = active & (radii >= lower) & (radii < upper)
        radial_bins.append(float(np.sum(weighted_load[mask]) / load_n))
    pressure_rms_pa = math.sqrt(float(np.sum(weighted_load * pressure) / np.sum(areas[active])))
    load_balance_error = abs(load_n - float(target_load_n)) / float(target_load_n) if target_load_n else 0.0
    return {
        "case_name": case.get("name", "unnamed"),
        "load_n": load_n,
        "target_load_n": target_load_n,
        "travel_m": float(travel_m),
        "clearance_m": float(case["clearance_m"]),
        "reduced_radius_m": reduced_radius(case),
        "layer_stiffness_pa_per_m": stiffness,
        "contact_area_m2": contact_area_m2,
        "contact_area_mm2": contact_area_m2 * 1.0e6,
        "equivalent_contact_radius_mm": equivalent_radius_m * 1.0e3,
        "mean_pressure_pa": load_n / contact_area_m2,
        "mean_pressure_mpa": load_n / contact_area_m2 / 1.0e6,
        "peak_pressure_pa": float(np.max(pressure)),
        "peak_pressure_mpa": float(np.max(pressure)) / 1.0e6,
        "pressure_rms_mpa": pressure_rms_pa / 1.0e6,
        "peak_radius_mm": float(radii[peak_index] * 1.0e3),
        "radial_load_centroid_mm": float(np.sum(weighted_load * radii) / load_n * 1.0e3),
        "inner_load_fraction": float(np.sum(weighted_load[inner_mask]) / load_n),
        "outer_load_fraction": float(np.sum(weighted_load[outer_mask]) / load_n),
        "radial_load_fraction_bins": radial_bins,
        "load_balance_relative_error": float(load_balance_error),
        "profile": {
            "radius_mm": (radii * 1.0e3).tolist(),
            "cell_area_mm2": (areas * 1.0e6).tolist(),
            "gap_mm": (gap * 1.0e3).tolist(),
            "pressure_mpa": (pressure / 1.0e6).tolist(),
        },
    }


def solve_travel_controlled(case: dict[str, Any], travel_m: float) -> dict[str, Any]:
    radii, _ = _grid(case)
    gap = gap_profile(case, radii)
    if travel_m <= float(np.min(gap)):
        raise ValueError("travel must exceed the minimum gap")
    return _state(case, float(travel_m))


def _load_at_travel(case: dict[str, Any], travel_m: float) -> float:
    radii, areas = _grid(case)
    gap = gap_profile(case, radii)
    _, _, stiffness = layer_properties(
        float(case["cartilage_shear_modulus_pa"]),
        float(case["poisson_ratio"]),
        float(case["layer_thickness_m"]),
    )
    return float(np.dot(np.where(travel_m - gap > 0.0, stiffness * (travel_m - gap), 0.0), areas))


def solve_load_controlled(case: dict[str, Any]) -> dict[str, Any]:
    target_load_n = float(case["load_n"])
    if target_load_n <= 0:
        raise ValueError("load_n must be positive")
    radii, areas = _grid(case)
    gap = gap_profile(case, radii)
    _, _, stiffness = layer_properties(
        float(case["cartilage_shear_modulus_pa"]),
        float(case["poisson_ratio"]),
        float(case["layer_thickness_m"]),
    )
    lower = float(np.min(gap))
    upper = float(np.max(gap) + target_load_n / (stiffness * float(np.sum(areas))))
    if _load_at_travel(case, lower) >= target_load_n:
        raise ValueError("invalid lower bracket")
    if _load_at_travel(case, upper) < target_load_n:
        raise ValueError("invalid upper bracket")
    for _ in range(120):
        midpoint = 0.5 * (lower + upper)
        if _load_at_travel(case, midpoint) > target_load_n:
            upper = midpoint
        else:
            lower = midpoint
    return _state(case, 0.5 * (lower + upper), target_load_n)


def analytic_parabolic_solution(case: dict[str, Any]) -> dict[str, float]:
    radius_m = reduced_radius(case)
    target_load_n = float(case["load_n"])
    _, _, stiffness = layer_properties(
        float(case["cartilage_shear_modulus_pa"]),
        float(case["poisson_ratio"]),
        float(case["layer_thickness_m"]),
    )
    contact_radius_m = (4.0 * radius_m * target_load_n / (math.pi * stiffness)) ** 0.25
    area_m2 = math.pi * contact_radius_m**2
    mean_pressure_pa = target_load_n / area_m2
    peak_pressure_pa = stiffness * contact_radius_m**2 / (2.0 * radius_m)
    return {
        "reduced_radius_m": radius_m,
        "contact_radius_m": contact_radius_m,
        "contact_area_m2": area_m2,
        "travel_m": float(case["clearance_m"]) + contact_radius_m**2 / (2.0 * radius_m),
        "mean_pressure_pa": mean_pressure_pa,
        "peak_pressure_pa": peak_pressure_pa,
        "load_n": target_load_n,
    }


def _dimension(unit: str) -> tuple[int, int, int, int]:
    mapping = {
        "1": (0, 0, 0, 0),
        "m": (1, 0, 0, 0),
        "mm": (1, 0, 0, 0),
        "m2": (2, 0, 0, 0),
        "N": (1, 1, -2, 0),
        "Pa": (-1, 1, -2, 0),
        "MPa": (-1, 1, -2, 0),
        "Pa/m": (-2, 1, -2, 0),
        "1/m": (-1, 0, 0, 0),
    }
    if unit not in mapping:
        raise ValueError(f"unsupported unit {unit}")
    return mapping[unit]


def _dim_mul(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(a + b for a, b in zip(left, right))


def _dim_div(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple(a - b for a, b in zip(left, right))


def _dim_add(left: tuple[int, int, int, int], right: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    if left != right:
        raise ValueError("dimensional addition requires identical dimensions")
    return left


def _dim_pow(value: tuple[int, int, int, int], exponent: float) -> tuple[int, int, int, int]:
    result = tuple(float(item) * exponent for item in value)
    if any(abs(item - round(item)) > 1e-12 for item in result):
        raise ValueError("non-integer dimensional exponent")
    return tuple(int(round(item)) for item in result)


def _format_dimension(value: tuple[int, int, int, int]) -> str:
    labels = ("m", "kg", "s", "K")
    powers = []
    for label, exponent in zip(labels, value):
        if exponent == 0:
            continue
        powers.append(label if exponent == 1 else f"{label}^{exponent}")
    return "1" if not powers else " ".join(powers)


def unit_check() -> dict[str, Any]:
    length = _dimension("m")
    force = _dimension("N")
    pressure = _dimension("Pa")
    stiffness = _dimension("Pa/m")
    area = _dimension("m2")
    curvature_term = _dim_div(_dim_pow(length, 2), length)
    inverse_length = _dim_div(_dimension("1"), length)
    inverse_curvature_difference = _dim_add(inverse_length, inverse_length)
    checks = [
        ("reduced radius", _dim_pow(inverse_curvature_difference, -1), length),
        ("layer stiffness", _dim_div(pressure, length), stiffness),
        ("curvature term", curvature_term, length),
        ("gap", length, length),
        ("pressure", _dim_mul(stiffness, length), pressure),
        ("load", _dim_mul(pressure, area), force),
        ("mean pressure", _dim_div(force, area), pressure),
        ("contact radius", _dim_pow(_dim_div(_dim_mul(length, force), stiffness), 0.25), length),
    ]
    rows = [
        {"name": name, "lhs": _format_dimension(lhs), "rhs": _format_dimension(rhs), "passed": lhs == rhs}
        for name, lhs, rhs in checks
    ]
    return {"passed": all(row["passed"] for row in rows), "checks": rows}


def parameter_table() -> list[dict[str, Any]]:
    cases = default_cases()
    nominal = cases["nominal"]
    as_built = cases["as_built"]
    nominal_radius = reduced_radius(nominal)
    as_built_radius = reduced_radius(as_built)
    cartilage_modulus, contact_modulus, stiffness = layer_properties(
        nominal["cartilage_shear_modulus_pa"], nominal["poisson_ratio"], nominal["layer_thickness_m"]
    )
    return [
        {"name": "load_n", "value": nominal["load_n"], "unit": "N", "kind": "assumption", "source_or_assumption": "synthetic fixed-load scenario; not patient measurement"},
        {"name": "head_radius_nominal_m", "value": nominal["head_radius_m"], "unit": "m", "kind": "assumption", "source_or_assumption": "synthetic local CAD-like geometry"},
        {"name": "socket_radius_nominal_m", "value": nominal["socket_radius_m"], "unit": "m", "kind": "assumption", "source_or_assumption": "synthetic local CAD-like geometry"},
        {"name": "head_radius_as_built_m", "value": as_built["head_radius_m"], "unit": "m", "kind": "synthetic perturbation", "source_or_assumption": "example of geometric deviation; not measurement data"},
        {"name": "socket_radius_as_built_m", "value": as_built["socket_radius_m"], "unit": "m", "kind": "synthetic perturbation", "source_or_assumption": "example of geometric deviation; not measurement data"},
        {"name": "reduced_radius_nominal_m", "value": nominal_radius, "unit": "m", "kind": "derived", "source_or_assumption": "(1/head_radius-1/socket_radius)^-1"},
        {"name": "reduced_radius_as_built_m", "value": as_built_radius, "unit": "m", "kind": "derived", "source_or_assumption": "(1/head_radius-1/socket_radius)^-1"},
        {"name": "clearance_as_built_m", "value": as_built["clearance_m"], "unit": "m", "kind": "synthetic perturbation", "source_or_assumption": "radial initial gap; not measurement data"},
        {"name": "form_error_amplitude_as_built_m", "value": as_built["form_error_amplitude_m"], "unit": "m", "kind": "synthetic perturbation", "source_or_assumption": "sinusoidal radial shape function; not measurement data"},
        {"name": "form_error_wavelength_m", "value": as_built["form_error_wavelength_m"], "unit": "m", "kind": "assumption", "source_or_assumption": "synthetic wavelength for shape deviation"},
        {"name": "cartilage_shear_modulus_pa", "value": nominal["cartilage_shear_modulus_pa"], "unit": "Pa", "kind": "source", "source_or_assumption": "Anderson et al. 2008, Methods; G=6.8 MPa"},
        {"name": "poisson_ratio", "value": nominal["poisson_ratio"], "unit": "1", "kind": "source/assumption", "source_or_assumption": "near-incompressible sensitivity value in Anderson et al. 2008"},
        {"name": "layer_thickness_m", "value": nominal["layer_thickness_m"], "unit": "m", "kind": "assumption", "source_or_assumption": "midpoint of reported femoral 1.5±0.5 mm and acetabular 1.6±0.4 mm"},
        {"name": "cartilage_modulus_pa", "value": cartilage_modulus, "unit": "Pa", "kind": "derived", "source_or_assumption": "E=2G(1+nu)"},
        {"name": "contact_modulus_pa", "value": contact_modulus, "unit": "Pa", "kind": "derived reduced model", "source_or_assumption": "E/[2(1-nu^2)]"},
        {"name": "layer_stiffness_pa_per_m", "value": stiffness, "unit": "Pa/m", "kind": "derived", "source_or_assumption": "E_contact/h"},
        {"name": "radial_cells", "value": nominal["n_cells"], "unit": "1", "kind": "numerical setting", "source_or_assumption": "equal-area annular discretization"},
        {"name": "radial_domain_radius_m", "value": nominal["radius_max_m"], "unit": "m", "kind": "numerical setting", "source_or_assumption": "integration domain"},
    ]


def _relative_change(value: float, baseline: float) -> float | None:
    if baseline == 0:
        return None
    return (value - baseline) / baseline


def compare_cases(nominal_case: dict[str, Any], as_built_case: dict[str, Any]) -> dict[str, Any]:
    nominal = solve_load_controlled(nominal_case)
    as_built = solve_load_controlled(as_built_case)
    as_built_same_travel = solve_travel_controlled(as_built_case, nominal["travel_m"])
    pressure_difference = np.asarray(as_built["profile"]["pressure_mpa"]) - np.asarray(nominal["profile"]["pressure_mpa"])
    pressure_rmse_mpa = float(np.sqrt(np.mean(pressure_difference**2)))
    return {
        "nominal": nominal,
        "as_built": as_built,
        "as_built_at_nominal_travel": as_built_same_travel,
        "deltas": {
            "contact_area_relative": _relative_change(as_built["contact_area_m2"], nominal["contact_area_m2"]),
            "mean_pressure_relative": _relative_change(as_built["mean_pressure_pa"], nominal["mean_pressure_pa"]),
            "peak_pressure_relative": _relative_change(as_built["peak_pressure_pa"], nominal["peak_pressure_pa"]),
            "travel_relative": _relative_change(as_built["travel_m"], nominal["travel_m"]),
            "radial_load_centroid_relative": _relative_change(as_built["radial_load_centroid_mm"], nominal["radial_load_centroid_mm"]),
            "outer_load_fraction_absolute": as_built["outer_load_fraction"] - nominal["outer_load_fraction"],
        },
        "pressure_field_rmse_mpa": pressure_rmse_mpa,
        "pressure_field_max_absolute_difference_mpa": float(np.max(np.abs(pressure_difference))),
        "same_travel_as_built_load_ratio": as_built_same_travel["load_n"] / nominal["load_n"],
        "same_travel_as_built_area_relative": _relative_change(as_built_same_travel["contact_area_m2"], nominal["contact_area_m2"]),
    }


def sensitivity(base_case: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    base = deepcopy(base_case or default_cases()["as_built"])
    baseline = solve_load_controlled(base)
    rows: list[dict[str, Any]] = []
    parameters = ("reduced_radius", "layer_modulus", "form_error_amplitude")
    for parameter in parameters:
        for factor in (0.5, 1.0, 1.5):
            case = deepcopy(base)
            if parameter == "reduced_radius":
                case["reduced_radius_m"] = reduced_radius(base) * factor
                input_value = case["reduced_radius_m"]
            elif parameter == "layer_modulus":
                case["cartilage_shear_modulus_pa"] = base["cartilage_shear_modulus_pa"] * factor
                input_value = case["cartilage_shear_modulus_pa"]
            else:
                case["form_error_amplitude_m"] = base["form_error_amplitude_m"] * factor
                input_value = case["form_error_amplitude_m"]
            state = solve_load_controlled(case)
            rows.append({
                "parameter": parameter,
                "factor": factor,
                "input_value": input_value,
                "input_unit": "m" if parameter != "layer_modulus" else "Pa",
                "contact_area_mm2": state["contact_area_mm2"],
                "mean_pressure_mpa": state["mean_pressure_mpa"],
                "peak_pressure_mpa": state["peak_pressure_mpa"],
                "outer_load_fraction": state["outer_load_fraction"],
                "travel_mm": state["travel_m"] * 1.0e3,
                "relative_area_change": _relative_change(state["contact_area_m2"], baseline["contact_area_m2"]),
                "relative_mean_pressure_change": _relative_change(state["mean_pressure_pa"], baseline["mean_pressure_pa"]),
                "relative_peak_pressure_change": _relative_change(state["peak_pressure_pa"], baseline["peak_pressure_pa"]),
                "load_balance_relative_error": state["load_balance_relative_error"],
            })
    return rows


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _write_json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_jsonable(value), handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def _write_sensitivity_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "parameter",
        "factor",
        "input_value",
        "input_unit",
        "contact_area_mm2",
        "mean_pressure_mpa",
        "peak_pressure_mpa",
        "outer_load_fraction",
        "travel_mm",
        "relative_area_change",
        "relative_mean_pressure_change",
        "relative_peak_pressure_change",
        "load_balance_relative_error",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _write_profile_csv(path: Path, state: dict[str, Any]) -> None:
    profile = state["profile"]
    fields = ["radius_mm", "cell_area_mm2", "gap_mm", "pressure_mpa"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in zip(*(profile[field] for field in fields)):
            writer.writerow(dict(zip(fields, row)))


def _validation(nominal_case: dict[str, Any], as_built_case: dict[str, Any], comparison: dict[str, Any]) -> dict[str, Any]:
    fine_case = deepcopy(nominal_case)
    fine_case["n_cells"] = 3000
    fine = solve_load_controlled(fine_case)
    analytic = analytic_parabolic_solution(nominal_case)
    analytic_errors = {
        "area_relative": abs(fine["contact_area_m2"] - analytic["contact_area_m2"]) / analytic["contact_area_m2"],
        "mean_pressure_relative": abs(fine["mean_pressure_pa"] - analytic["mean_pressure_pa"]) / analytic["mean_pressure_pa"],
        "peak_pressure_relative": abs(fine["peak_pressure_pa"] - analytic["peak_pressure_pa"]) / analytic["peak_pressure_pa"],
    }
    placebo_case = deepcopy(nominal_case)
    placebo_case["clearance_m"] = 0.1e-3
    placebo = solve_load_controlled(placebo_case)
    placebo_pressure_difference = np.max(np.abs(np.asarray(placebo["profile"]["pressure_mpa"]) - np.asarray(comparison["nominal"]["profile"]["pressure_mpa"])))
    placebo_area_difference = abs(placebo["contact_area_m2"] - comparison["nominal"]["contact_area_m2"]) / comparison["nominal"]["contact_area_m2"]
    placebo_travel_error = abs((placebo["travel_m"] - comparison["nominal"]["travel_m"]) - placebo_case["clearance_m"])
    unit_result = unit_check()
    source_area_low, source_area_high = REFERENCE["experimental_contact_area_mm2_range"]
    source_pressure_low, source_pressure_high = REFERENCE["experimental_mean_pressure_mpa_range"]
    source_scale_pass = source_area_low <= comparison["nominal"]["contact_area_mm2"] <= source_area_high and source_pressure_low <= comparison["nominal"]["mean_pressure_mpa"] <= source_pressure_high
    direction_pass = comparison["as_built"]["contact_area_m2"] < comparison["nominal"]["contact_area_m2"] and comparison["as_built"]["mean_pressure_pa"] > comparison["nominal"]["mean_pressure_pa"]
    form_distribution_pass = comparison["pressure_field_max_absolute_difference_mpa"] > 1.0e-9
    balance_pass = comparison["nominal"]["load_balance_relative_error"] < 1.0e-10 and comparison["as_built"]["load_balance_relative_error"] < 1.0e-10
    checks = {
        "unit_check": unit_result["passed"],
        "analytic_limit": max(analytic_errors.values()) < 0.01,
        "clearance_placebo": placebo_pressure_difference < 1.0e-9 and placebo_area_difference < 1.0e-9 and placebo_travel_error < 1.0e-12,
        "load_balance": balance_pass,
        "source_scale": source_scale_pass,
        "as_built_direction": direction_pass,
        "form_distribution": form_distribution_pass,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "analytic_reference": analytic,
        "analytic_relative_errors": analytic_errors,
        "placebo": {
            "added_clearance_mm": placebo_case["clearance_m"] * 1.0e3,
            "pressure_max_difference_mpa": float(placebo_pressure_difference),
            "area_relative_difference": float(placebo_area_difference),
            "travel_error_m": float(placebo_travel_error),
        },
        "source_scale_check": {
            "nominal_area_mm2": comparison["nominal"]["contact_area_mm2"],
            "nominal_mean_pressure_mpa": comparison["nominal"]["mean_pressure_mpa"],
            "reference_area_mm2": REFERENCE["experimental_contact_area_mm2_range"],
            "reference_mean_pressure_mpa": REFERENCE["experimental_mean_pressure_mpa_range"],
            "passed": source_scale_pass,
        },
    }


def run_experiment(output_dir: str | Path = BASE_DIR) -> dict[str, Any]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    cases = default_cases()
    comparison = compare_cases(cases["nominal"], cases["as_built"])
    sensitivity_rows = sensitivity(cases["as_built"])
    validation = _validation(cases["nominal"], cases["as_built"], comparison)
    hash_path = output_path / "PREREG.sha256"
    prereg_hash = hash_path.read_text(encoding="utf-8").split()[0] if hash_path.exists() else "MISSING"
    result = {
        "id": "BT-HX-Q054",
        "status": "first_runnable_mechanistic_model",
        "prereg_sha256": prereg_hash,
        "source_reference": REFERENCE,
        "equations": EQUATIONS,
        "parameter_table": parameter_table(),
        "unit_check": unit_check(),
        "cases": {
            "nominal": comparison["nominal"],
            "as_built": comparison["as_built"],
            "as_built_at_nominal_travel": comparison["as_built_at_nominal_travel"],
        },
        "comparison": {
            key: value for key, value in comparison.items() if key not in {"nominal", "as_built", "as_built_at_nominal_travel"}
        },
        "sensitivity": sensitivity_rows,
        "validation": validation,
        "limitations": [
            "The radial model is axisymmetric and does not predict a two-dimensional center of pressure.",
            "The Winkler layer omits elastic coupling, cartilage poroelasticity, fluid pressurisation, friction and sliding.",
            "Synthetic radii, clearance and form error are assumptions and are not measurements.",
            "The published range is a scale reference; no subject-specific validation was performed.",
        ],
        "next_resolution_step": [
            "Replace synthetic geometry with registered CAD and as-built radial or surface deviation maps.",
            "Measure load, travel and pressure or deformation on the same specimen in an independent contact test.",
            "Upgrade the layer foundation to a layered or finite-element contact model and test angular load transfer and shear.",
        ],
    }
    _write_json(output_path / "results.json", result)
    _write_sensitivity_csv(output_path / "sensitivity.csv", sensitivity_rows)
    _write_profile_csv(output_path / "profile_nominal.csv", comparison["nominal"])
    _write_profile_csv(output_path / "profile_as_built.csv", comparison["as_built"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=str(BASE_DIR))
    args = parser.parse_args()
    result = run_experiment(args.output_dir)
    print(json.dumps(_jsonable({
        "id": result["id"],
        "validation_passed": result["validation"]["passed"],
        "nominal_area_mm2": result["cases"]["nominal"]["contact_area_mm2"],
        "as_built_area_mm2": result["cases"]["as_built"]["contact_area_mm2"],
        "nominal_mean_pressure_mpa": result["cases"]["nominal"]["mean_pressure_mpa"],
        "as_built_mean_pressure_mpa": result["cases"]["as_built"]["mean_pressure_mpa"],
    }), ensure_ascii=False))


if __name__ == "__main__":
    main()
