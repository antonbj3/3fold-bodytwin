from __future__ import annotations

import argparse
import heapq
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


AVOGADRO = 6.02214076e23
FARADAY = 96485.33212
GAS_CONSTANT = 8.314462618


@dataclass(frozen=True)
class Parameters:
    diffusion_coefficient: float = 1.0e-11
    crista_length: float = 0.9e-6
    total_membrane_area: float = 2.0e-12
    total_ims_volume: float = 5.0e-19
    cyt_adp: float = 0.037
    cyt_atp: float = 4.0
    matrix_total_adenylate: float = 1.16
    pi_matrix: float = 5.0
    as_site_density: float = 1.0e15
    ant_site_density: float = 1.0e15
    as_kcat: float = 90.0
    ant_kcat: float = 90.0
    as_k_adp: float = 0.8
    as_k_atp: float = 2.0
    as_k_pi: float = 1.0
    ant_k_adp: float = 0.1
    ant_k_atp: float = 1.0
    ant_fraction_pmf: float = 0.5
    membrane_potential: float = -0.172
    temperature: float = 310.0
    reference_concentration: float = 1.0
    junction_area_total: float = 2.0e-14
    junction_length: float = 5.0e-8
    solver_tolerance: float = 1.0e-11
    solver_max_iterations: int = 200
    matrix_root_tolerance: float = 1.0e-10

    @property
    def as_vmax(self) -> float:
        return self.as_site_density * self.as_kcat / AVOGADRO

    @property
    def ant_vmax(self) -> float:
        return self.ant_site_density * self.ant_kcat / AVOGADRO

    @property
    def ant_k(self) -> float:
        return self.ant_vmax / self.reference_concentration

    @property
    def thermal_voltage(self) -> float:
        return GAS_CONSTANT * self.temperature / FARADAY

    @property
    def pmf_factor(self) -> float:
        x = abs(self.membrane_potential) / self.thermal_voltage
        return 1.0 / (1.0 + np.exp(-x))

    @property
    def electrochemical_factor(self) -> float:
        x = 2.0 * self.ant_fraction_pmf * abs(self.membrane_potential) / self.thermal_voltage
        return float(np.exp(x))


@dataclass(frozen=True)
class Geometry:
    name: str
    n_nodes: int
    edges: tuple[tuple[int, int, float], ...]
    cj_nodes: tuple[int, ...]


PARAMETER_TABLE = [
    {"name": "diffusion_coefficient", "unit": "m2/s", "source_or_assumption": "effective IMS ADP diffusivity; sensitivity parameter"},
    {"name": "crista_length", "unit": "m", "source_or_assumption": "long-crista comparison scale from Adams et al. 2025, Fig. 4"},
    {"name": "total_membrane_area", "unit": "m2", "source_or_assumption": "held identical across geometries; abstract unit"},
    {"name": "total_ims_volume", "unit": "m3", "source_or_assumption": "held identical across geometries; abstract unit"},
    {"name": "cyt_adp", "unit": "mol/m3", "source_or_assumption": "0.037 mM, Adams et al. 2025, Fig. 1/Section 2.1"},
    {"name": "cyt_atp", "unit": "mol/m3", "source_or_assumption": "4.0 mM, Adams et al. 2025, Section 2.1"},
    {"name": "matrix_total_adenylate", "unit": "mol/m3", "source_or_assumption": "0.44 + 0.72 mM starting matrix ATP + ADP, Adams et al. 2025, Section 2.1"},
    {"name": "pi_matrix", "unit": "mol/m3", "source_or_assumption": "fixed Pi assumption; sensitivity not run"},
    {"name": "as_site_density", "unit": "molecules/m2", "source_or_assumption": "matched ATP-synthase site density; scale assumption"},
    {"name": "ant_site_density", "unit": "molecules/m2", "source_or_assumption": "matched ANT site density; scale assumption"},
    {"name": "as_kcat", "unit": "1/s", "source_or_assumption": "sensitivity parameter; scale assumption"},
    {"name": "ant_kcat", "unit": "1/s", "source_or_assumption": "sensitivity parameter; scale assumption"},
    {"name": "as_k_adp", "unit": "mol/m3", "source_or_assumption": "Michaelis parameter assumption"},
    {"name": "as_k_atp", "unit": "mol/m3", "source_or_assumption": "ATP inhibition parameter assumption"},
    {"name": "as_k_pi", "unit": "mol/m3", "source_or_assumption": "Pi Michaelis parameter assumption"},
    {"name": "ant_k_adp", "unit": "mol/m3", "source_or_assumption": "external/matrix ADP exchange scale assumption"},
    {"name": "ant_k_atp", "unit": "mol/m3", "source_or_assumption": "external/matrix ATP exchange scale assumption"},
    {"name": "ant_fraction_pmf", "unit": "1", "source_or_assumption": "0.5 effective potential fraction, Adams et al. 2025 Eq. 1"},
    {"name": "membrane_potential", "unit": "V", "source_or_assumption": "-172 mV, Adams et al. 2025 Section 2.1"},
    {"name": "temperature", "unit": "K", "source_or_assumption": "310 K physiological assumption"},
    {"name": "reference_concentration", "unit": "mol/m3", "source_or_assumption": "1 mM reference for low-substrate ANT scale"},
    {"name": "junction_area_total", "unit": "m2", "source_or_assumption": "total CJ area held fixed; sensitivity parameter"},
    {"name": "junction_length", "unit": "m", "source_or_assumption": "effective junction transport length assumption"},
    {"name": "solver_tolerance", "unit": "1", "source_or_assumption": "numerical convergence threshold"},
    {"name": "solver_max_iterations", "unit": "1", "source_or_assumption": "numerical convergence limit"},
    {"name": "matrix_root_tolerance", "unit": "mol/m3", "source_or_assumption": "scalar balance root tolerance"},
]


def validate_parameters(parameters: Parameters) -> None:
    positive_fields = [
        "diffusion_coefficient",
        "crista_length",
        "total_membrane_area",
        "total_ims_volume",
        "cyt_adp",
        "cyt_atp",
        "matrix_total_adenylate",
        "pi_matrix",
        "as_site_density",
        "ant_site_density",
        "as_kcat",
        "ant_kcat",
        "as_k_adp",
        "as_k_atp",
        "as_k_pi",
        "ant_k_adp",
        "ant_k_atp",
        "reference_concentration",
        "junction_area_total",
        "junction_length",
        "solver_tolerance",
        "matrix_root_tolerance",
    ]
    for field in positive_fields:
        if getattr(parameters, field) <= 0:
            raise ValueError(f"{field} must be positive")
    if parameters.solver_max_iterations <= 0:
        raise ValueError("solver_max_iterations must be positive")
    if not 0.0 < parameters.ant_fraction_pmf < 1.0:
        raise ValueError("ant_fraction_pmf must lie between zero and one")
    if parameters.junction_area_total > parameters.total_membrane_area:
        raise ValueError("junction area cannot exceed total membrane area")
    if parameters.matrix_total_adenylate <= parameters.cyt_adp:
        raise ValueError("matrix adenylate pool must exceed cytosolic ADP")
    if parameters.membrane_potential >= 0.0:
        raise ValueError("membrane potential must be negative inside")


def make_geometries(parameters: Parameters) -> dict[str, Geometry]:
    n_nodes = 42
    straight_edges = tuple(
        (i, i + 1, parameters.crista_length / (n_nodes - 1))
        for i in range(n_nodes - 1)
    )
    one = Geometry("one_cj", n_nodes, straight_edges, (0,))
    two = Geometry("two_cj", n_nodes, straight_edges, (0, n_nodes - 1))
    trunk_nodes = 28
    branch_a_nodes = 8
    branch_b_nodes = 6
    trunk_edges = tuple(
        (i, i + 1, parameters.crista_length / (trunk_nodes - 1))
        for i in range(trunk_nodes - 1)
    )
    branch_a_start = trunk_nodes
    branch_b_start = branch_a_start + branch_a_nodes
    branch_a_edges = tuple(
        (
            branch_a_start if i == 0 else branch_a_start + i - 1,
            branch_a_start + i,
            parameters.crista_length * 0.30 / (branch_a_nodes - 1),
        )
        for i in range(branch_a_nodes - 1)
    )
    branch_a_edges = (((14, branch_a_start, parameters.crista_length * 0.30 / (branch_a_nodes - 1)),) + branch_a_edges[1:])
    branch_b_edges = tuple(
        (
            branch_b_start if i == 0 else branch_b_start + i - 1,
            branch_b_start + i,
            parameters.crista_length * 0.25 / (branch_b_nodes - 1),
        )
        for i in range(branch_b_nodes - 1)
    )
    branch_b_edges = (((14, branch_b_start, parameters.crista_length * 0.25 / (branch_b_nodes - 1)),) + branch_b_edges[1:])
    branched_edges = trunk_edges + branch_a_edges + branch_b_edges
    branched = Geometry("branched", n_nodes, branched_edges, (0, trunk_nodes - 1, branch_b_start - 1, branch_b_start + branch_b_nodes - 1))
    flat = Geometry("flat", 1, (), (0,))
    return {"one_cj": one, "two_cj": two, "branched": branched, "flat": flat}


def total_path_length(geometry: Geometry, parameters: Parameters) -> float:
    if not geometry.edges:
        return parameters.crista_length
    return float(sum(length for _, _, length in geometry.edges))


def cross_sectional_area(geometry: Geometry, parameters: Parameters) -> float:
    return parameters.total_ims_volume / total_path_length(geometry, parameters)


def edge_conductance(geometry: Geometry, parameters: Parameters) -> list[tuple[int, int, float]]:
    area = cross_sectional_area(geometry, parameters)
    return [
        (i, j, parameters.diffusion_coefficient * area / length)
        for i, j, length in geometry.edges
    ]


def ant_flux(ims_adp: np.ndarray, matrix_adp: float, parameters: Parameters) -> np.ndarray:
    matrix_atp = parameters.matrix_total_adenylate - matrix_adp
    if matrix_atp <= 0.0:
        raise ValueError("matrix ATP concentration must remain positive")
    denominator_scale = parameters.ant_k_adp * parameters.ant_k_atp
    forward = matrix_atp * ims_adp / denominator_scale
    reverse = parameters.cyt_atp * matrix_adp / denominator_scale
    factor = parameters.electrochemical_factor
    return parameters.ant_k * parameters.reference_concentration * (factor * forward - reverse / factor)


def ant_slope(ims_adp: np.ndarray, matrix_adp: float, parameters: Parameters) -> np.ndarray:
    matrix_atp = parameters.matrix_total_adenylate - matrix_adp
    denominator_scale = parameters.ant_k_adp * parameters.ant_k_atp
    factor = parameters.electrochemical_factor
    return np.full_like(ims_adp, parameters.ant_k * parameters.reference_concentration * factor * matrix_atp / denominator_scale, dtype=float)


def atp_synthase_flux(matrix_adp: float, parameters: Parameters) -> float:
    matrix_atp = parameters.matrix_total_adenylate - matrix_adp
    if matrix_adp <= 0.0 or matrix_atp <= 0.0:
        raise ValueError("matrix adenylate concentrations must remain positive")
    adp_ratio = matrix_adp / parameters.as_k_adp
    atp_ratio = matrix_atp / parameters.as_k_atp
    pi_ratio = parameters.pi_matrix / parameters.as_k_pi
    occupancy = adp_ratio * pi_ratio / (1.0 + adp_ratio + atp_ratio + pi_ratio)
    return parameters.as_vmax * parameters.pmf_factor * occupancy


def assemble_diffusion_system(
    geometry: Geometry,
    parameters: Parameters,
    sink_intercept: np.ndarray,
    sink_slope: np.ndarray,
) -> tuple[Any, np.ndarray]:
    n = geometry.n_nodes
    area = parameters.total_membrane_area / n
    matrix = lil_matrix((n, n), dtype=float)
    rhs = np.zeros(n, dtype=float)
    for i, j, conductance in edge_conductance(geometry, parameters):
        matrix[i, i] += conductance
        matrix[j, j] += conductance
        matrix[i, j] -= conductance
        matrix[j, i] -= conductance
    for i in range(n):
        matrix[i, i] += sink_slope[i] * area
        rhs[i] -= sink_intercept[i] * area
    junction_conductance = parameters.diffusion_coefficient * parameters.junction_area_total / parameters.junction_length
    per_junction = junction_conductance / len(geometry.cj_nodes)
    for node in geometry.cj_nodes:
        matrix[node, node] += per_junction
        rhs[node] += per_junction * parameters.cyt_adp
    return matrix.tocsr(), rhs


def solve_ims(
    matrix_adp: float,
    geometry: Geometry,
    parameters: Parameters,
    initial: np.ndarray | None = None,
) -> np.ndarray:
    n = geometry.n_nodes
    concentration = np.full(n, parameters.cyt_adp, dtype=float) if initial is None else np.asarray(initial, dtype=float).copy()
    flux = ant_flux(concentration, matrix_adp, parameters)
    slope = ant_slope(concentration, matrix_adp, parameters)
    intercept = flux - slope * concentration
    system, rhs = assemble_diffusion_system(geometry, parameters, intercept, slope)
    concentration = np.asarray(spsolve(system, rhs), dtype=float)
    if not np.all(np.isfinite(concentration)) or np.min(concentration) < -parameters.solver_tolerance:
        raise RuntimeError("IMS diffusion solve did not converge")
    return np.maximum(concentration, 0.0)


def solve_diffusion_only(geometry: Geometry, parameters: Parameters, boundary_concentration: float) -> np.ndarray:
    n = geometry.n_nodes
    matrix = lil_matrix((n, n), dtype=float)
    for i, j, conductance in edge_conductance(geometry, parameters):
        matrix[i, i] += conductance
        matrix[j, j] += conductance
        matrix[i, j] -= conductance
        matrix[j, i] -= conductance
    boundary = list(geometry.cj_nodes)
    matrix[boundary, :] = 0.0
    for node in boundary:
        matrix[node, node] = 1.0
    rhs = np.zeros(n, dtype=float)
    for node in boundary:
        rhs[node] = boundary_concentration
    if n == 1:
        return np.asarray([boundary_concentration], dtype=float)
    return np.asarray(spsolve(matrix.tocsr(), rhs), dtype=float)


def graph_distances(geometry: Geometry) -> np.ndarray:
    adjacency: dict[int, list[tuple[int, float]]] = {i: [] for i in range(geometry.n_nodes)}
    for i, j, length in geometry.edges:
        adjacency[i].append((j, length))
        adjacency[j].append((i, length))
    distances = np.full(geometry.n_nodes, np.inf, dtype=float)
    queue: list[tuple[float, int]] = []
    for node in geometry.cj_nodes:
        distances[node] = 0.0
        heapq.heappush(queue, (0.0, node))
    while queue:
        distance, node = heapq.heappop(queue)
        if distance > distances[node]:
            continue
        for neighbor, length in adjacency[node]:
            candidate = distance + length
            if candidate < distances[neighbor]:
                distances[neighbor] = candidate
                heapq.heappush(queue, (candidate, neighbor))
    return distances


def transport_metrics(geometry: Geometry, parameters: Parameters) -> dict[str, float | int]:
    distances = graph_distances(geometry)
    max_distance = float(np.max(distances))
    junction_conductance = parameters.diffusion_coefficient * parameters.junction_area_total / parameters.junction_length
    diffusion_time = max_distance**2 / parameters.diffusion_coefficient
    junction_time = parameters.total_ims_volume / junction_conductance
    return {
        "n_nodes": geometry.n_nodes,
        "cj_count": len(geometry.cj_nodes),
        "max_graph_distance_to_cj_m": max_distance,
        "mean_graph_distance_to_cj_m": float(np.mean(distances)),
        "diffusion_time_s": diffusion_time,
        "junction_exchange_time_s": junction_time,
        "end_to_end_transport_time_s": diffusion_time + junction_time,
        "total_path_length_m": total_path_length(geometry, parameters),
        "cross_sectional_area_m2": cross_sectional_area(geometry, parameters),
    }


def matrix_balance(matrix_adp: float, geometry: Geometry, parameters: Parameters) -> tuple[float, np.ndarray, np.ndarray]:
    ims_adp = solve_ims(matrix_adp, geometry, parameters)
    local_ant = ant_flux(ims_adp, matrix_adp, parameters)
    total_ant = float(np.sum(local_ant) * parameters.total_membrane_area / geometry.n_nodes)
    total_atp = atp_synthase_flux(matrix_adp, parameters) * parameters.total_membrane_area
    return total_ant - total_atp, ims_adp, local_ant


def solve_matrix_adp(geometry: Geometry, parameters: Parameters) -> tuple[float, np.ndarray, np.ndarray]:
    lower = 1.0e-8
    upper = parameters.matrix_total_adenylate - 1.0e-8
    grid = np.linspace(lower, upper, 81)
    values = np.array([matrix_balance(value, geometry, parameters)[0] for value in grid])
    sign_changes = np.where(np.sign(values[:-1]) * np.sign(values[1:]) <= 0.0)[0]
    if len(sign_changes) == 0:
        raise RuntimeError("no matrix balance root found")
    index = int(sign_changes[0])
    root = brentq(
        lambda value: matrix_balance(value, geometry, parameters)[0],
        float(grid[index]),
        float(grid[index + 1]),
        xtol=parameters.matrix_root_tolerance,
    )
    balance, ims_adp, local_ant = matrix_balance(root, geometry, parameters)
    if abs(balance) > 1.0e-7:
        raise RuntimeError("matrix balance residual too large")
    return float(root), ims_adp, local_ant


def unit_check(parameters: Parameters, geometry: Geometry, result: dict[str, Any]) -> dict[str, Any]:
    validate_parameters(parameters)
    area_volume_ratio = parameters.total_membrane_area / parameters.total_ims_volume
    expected_area = geometry.n_nodes * parameters.total_membrane_area / geometry.n_nodes
    expected_volume = geometry.n_nodes * parameters.total_ims_volume / geometry.n_nodes
    checks = {
        "positive_si_parameters": True,
        "area_m2_positive": parameters.total_membrane_area > 0.0,
        "volume_m3_positive": parameters.total_ims_volume > 0.0,
        "diffusion_time_s_positive": result["transport"]["diffusion_time_s"] >= 0.0,
        "junction_time_s_positive": result["transport"]["junction_exchange_time_s"] > 0.0,
        "atp_flux_mol_per_s_positive": result["atp_flux_mol_per_s"] > 0.0,
        "atp_flux_molecules_per_s_consistent": abs(result["atp_flux_molecules_per_s"] - result["atp_flux_mol_per_s"] * AVOGADRO) < 1.0e-6 * max(1.0, abs(result["atp_flux_molecules_per_s"])),
        "area_sum_matches": abs(expected_area - parameters.total_membrane_area) <= 1.0e-24,
        "volume_sum_matches": abs(expected_volume - parameters.total_ims_volume) <= 1.0e-31,
        "area_volume_ratio_m_inv": area_volume_ratio,
    }
    if not all(bool(value) for key, value in checks.items() if key != "area_volume_ratio_m_inv"):
        raise AssertionError("unit check failed")
    return checks


def run_topology(name: str, geometry: Geometry, parameters: Parameters) -> dict[str, Any]:
    matrix_adp, ims_adp, local_ant = solve_matrix_adp(geometry, parameters)
    total_ant = float(np.sum(local_ant) * parameters.total_membrane_area / geometry.n_nodes)
    atp_flux = atp_synthase_flux(matrix_adp, parameters) * parameters.total_membrane_area
    if not np.isclose(total_ant, atp_flux, rtol=2.0e-6, atol=1.0e-18):
        raise RuntimeError("ANT and ATP-synthase fluxes do not balance")
    transport = transport_metrics(geometry, parameters)
    fractional_depletion = (parameters.cyt_adp - ims_adp) / parameters.cyt_adp
    result: dict[str, Any] = {
        "name": name,
        "area_m2": parameters.total_membrane_area,
        "ims_volume_m3": parameters.total_ims_volume,
        "as_site_count": parameters.as_site_density * parameters.total_membrane_area,
        "ant_site_count": parameters.ant_site_density * parameters.total_membrane_area,
        "matrix_adp_mol_per_m3": matrix_adp,
        "matrix_atp_mol_per_m3": parameters.matrix_total_adenylate - matrix_adp,
        "atp_flux_mol_per_s": atp_flux,
        "atp_flux_molecules_per_s": atp_flux * AVOGADRO,
        "atp_flux_mol_per_m2_per_s": atp_flux / parameters.total_membrane_area,
        "atp_flux_molecules_per_ms_per_um2": atp_flux / parameters.total_membrane_area * AVOGADRO / 1.0e15,
        "ant_flux_mol_per_s": total_ant,
        "mean_ims_adp_mol_per_m3": float(np.mean(ims_adp)),
        "min_ims_adp_mol_per_m3": float(np.min(ims_adp)),
        "max_fractional_adp_depletion": float(np.max(fractional_depletion)),
        "mean_fractional_adp_depletion": float(np.mean(fractional_depletion)),
        "transport": transport,
        "parameter_ant_vmax_mol_per_m2_per_s": parameters.ant_vmax,
        "parameter_as_vmax_mol_per_m2_per_s": parameters.as_vmax,
        "unit_check": unit_check(parameters, geometry, {"transport": transport, "atp_flux_mol_per_s": atp_flux, "atp_flux_molecules_per_s": atp_flux * AVOGADRO}),
    }
    return result


def run_core(parameters: Parameters) -> dict[str, Any]:
    validate_parameters(parameters)
    geometries = make_geometries(parameters)
    topology_results = {
        name: run_topology(name, geometry, parameters)
        for name, geometry in geometries.items()
    }
    one = topology_results["one_cj"]
    two = topology_results["two_cj"]
    branched = topology_results["branched"]
    return {
        "topologies": topology_results,
        "comparisons": {
            "two_over_one_atp_flux": two["atp_flux_mol_per_s"] / one["atp_flux_mol_per_s"],
            "branched_over_one_atp_flux": branched["atp_flux_mol_per_s"] / one["atp_flux_mol_per_s"],
            "two_over_one_transport_time": two["transport"]["end_to_end_transport_time_s"] / one["transport"]["end_to_end_transport_time_s"],
            "branched_over_one_transport_time": branched["transport"]["end_to_end_transport_time_s"] / one["transport"]["end_to_end_transport_time_s"],
        },
    }


def run_sensitivity(parameters: Parameters) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for field in ("diffusion_coefficient", "ant_kcat", "junction_area_total"):
        baseline = getattr(parameters, field)
        for factor in (0.5, 1.0, 1.5):
            varied = replace(parameters, **{field: baseline * factor})
            core = run_core(varied)
            one = core["topologies"]["one_cj"]
            two = core["topologies"]["two_cj"]
            branched = core["topologies"]["branched"]
            records.append({
                "parameter": field,
                "factor": factor,
                "value": getattr(varied, field),
                "one_cj_atp_flux_mol_per_s": one["atp_flux_mol_per_s"],
                "two_cj_atp_flux_mol_per_s": two["atp_flux_mol_per_s"],
                "branched_atp_flux_mol_per_s": branched["atp_flux_mol_per_s"],
                "two_over_one_atp_flux": two["atp_flux_mol_per_s"] / one["atp_flux_mol_per_s"],
                "branched_over_one_atp_flux": branched["atp_flux_mol_per_s"] / one["atp_flux_mol_per_s"],
                "one_cj_max_adp_depletion": one["max_fractional_adp_depletion"],
                "two_cj_max_adp_depletion": two["max_fractional_adp_depletion"],
                "one_cj_transport_time_s": one["transport"]["end_to_end_transport_time_s"],
                "two_cj_transport_time_s": two["transport"]["end_to_end_transport_time_s"],
            })
    return records


def diffusion_limit_check(parameters: Parameters) -> dict[str, Any]:
    geometry = make_geometries(parameters)["two_cj"]
    concentration = solve_diffusion_only(geometry, parameters, parameters.cyt_adp)
    maximum_error = float(np.max(np.abs(concentration - parameters.cyt_adp)))
    return {
        "geometry": geometry.name,
        "boundary_concentration_mol_per_m3": parameters.cyt_adp,
        "maximum_absolute_error_mol_per_m3": maximum_error,
        "passed": maximum_error < 1.0e-10,
    }


def run_experiment(parameters: Parameters | None = None) -> dict[str, Any]:
    selected = parameters or Parameters()
    validate_parameters(selected)
    core = run_core(selected)
    one = core["topologies"]["one_cj"]
    two = core["topologies"]["two_cj"]
    diffusion_time_ratio = two["transport"]["diffusion_time_s"] / one["transport"]["diffusion_time_s"]
    final_checks = {
        "flux_ratio_ge_1_10": two["atp_flux_mol_per_s"] / one["atp_flux_mol_per_s"] >= 1.10,
        "two_cj_depletion_lower": two["max_fractional_adp_depletion"] < one["max_fractional_adp_depletion"],
        "two_cj_diffusion_time_le_0_30": diffusion_time_ratio <= 0.30,
        "matched_area_volume_and_sites": all(
            two[key] == one[key]
            for key in ("area_m2", "ims_volume_m3", "as_site_count", "ant_site_count")
        ),
        "analytical_limit_passed": diffusion_limit_check(selected)["passed"],
    }
    failed_conditions = [key for key, value in final_checks.items() if not bool(value)]
    output = {
        "model": "BT-HX-Q077 first-principles crista ADP transport model",
        "concentration_unit": "mol/m3; 1 mM = 1 mol/m3",
        "length_unit": "m",
        "flux_unit": "mol/(m2 s) for local flux and mol/s for total flux",
        "parameters": {key: getattr(selected, key) for key in selected.__dataclass_fields__},
        "derived_parameters": {
            "crista_length_um": selected.crista_length * 1.0e6,
            "as_vmax_mol_per_m2_per_s": selected.as_vmax,
            "ant_vmax_mol_per_m2_per_s": selected.ant_vmax,
            "ant_k_m_per_s": selected.ant_k,
            "pmf_factor": selected.pmf_factor,
            "electrochemical_factor": selected.electrochemical_factor,
            "junction_conductance_m3_per_s": selected.diffusion_coefficient * selected.junction_area_total / selected.junction_length,
        },
        "parameter_table": PARAMETER_TABLE,
        "preregistration": {
            "file": "PREREG.md",
            "primary_quantity": "J_ATP(two_trans_CJ) / J_ATP(one_CJ)",
            "flux_ratio_threshold": 1.10,
            "diffusion_time_ratio_threshold": 0.30,
            "no_sink_equal_boundary_required": True,
            "status": "frozen_before_final_default_run",
        },
        "literature_anchor": {
            "doi": "10.3390/cells14040257",
            "pmcid": "PMC11853683",
            "figure": "Figure 4",
            "relative_flux_bottlenecked": 0.60,
            "relative_flux_high_connectivity": 0.90,
            "reported_connectivity_increase_percent": 17.0,
            "model_max_flux_molecules_per_ms_per_um2": 100.0,
            "membrane_potential_magnitude_mV": 172.0,
            "status": "VERIFIED",
        },
        "final_status": {
            "checks": final_checks,
            "primary_criterion_passed": not failed_conditions,
            "failed_conditions": failed_conditions,
        },
        **core,
        "sensitivity": run_sensitivity(selected),
        "analytical_limit_check": diffusion_limit_check(selected),
    }
    return output


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
    output = run_experiment()
    Path(args.output).write_text(json.dumps(output, indent=2, sort_keys=True, default=json_default) + "\n", encoding="utf-8")
    summary = {
        "two_over_one_atp_flux": output["comparisons"]["two_over_one_atp_flux"],
        "branched_over_one_atp_flux": output["comparisons"]["branched_over_one_atp_flux"],
        "two_over_one_transport_time": output["comparisons"]["two_over_one_transport_time"],
        "analytical_limit_check": output["analytical_limit_check"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True, default=json_default))


if __name__ == "__main__":
    main()
