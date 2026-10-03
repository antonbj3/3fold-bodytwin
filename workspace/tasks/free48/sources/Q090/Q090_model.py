from __future__ import annotations

from dataclasses import dataclass, fields, replace
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root


@dataclass(frozen=True)
class FiberFamily:
    name: str
    area_fraction: float
    modulus_pa: float
    angle_rad: float


@dataclass(frozen=True)
class Parameters:
    pdl_fibers: tuple[FiberFamily, ...]
    molar_support_area_m2: float
    pdl_thickness_m: float
    pdl_strain_scale: float
    pdl_stiffening_coefficient: float
    alveolar_bone_modulus_pa: float
    alveolar_bone_thickness_m: float
    fluid_viscosity_pa_s: float
    pdl_permeability_m2: float
    tmj_joint_area_m2: float
    tmj_disc_thickness_m: float
    tmj_disc_modulus_pa: float
    tmj_condylar_cartilage_thickness_m: float
    tmj_condylar_cartilage_modulus_pa: float
    tmj_temporal_cartilage_thickness_m: float
    tmj_temporal_cartilage_modulus_pa: float
    tmj_displacement_scale_m: float
    tmj_stiffening_coefficient: float
    tmj_retrodiscal_tension_fraction: float
    molar_lever_m: float
    condyle_lever_m: float
    effective_condyle_radius_m: float
    close_load_n: float
    external_moment_nm: float
    mandible_mass_kg: float
    mandible_inertia_kgm2: float
    load_ramp_s: float
    simulation_time_s: float
    root_tolerance: float
    integration_max_step_s: float
    reference_mobility_m: float


DEFAULT_PARAMETERS = Parameters(
    pdl_fibers=(
        FiberFamily("vertical", 0.30, 2.0e6, 0.0),
        FiberFamily("oblique_45deg", 0.55, 1.0e6, math.radians(45.0)),
        FiberFamily("circumferential", 0.15, 3.0e6, math.radians(90.0)),
    ),
    molar_support_area_m2=2.0e-4,
    pdl_thickness_m=2.0e-4,
    pdl_strain_scale=0.20,
    pdl_stiffening_coefficient=2.0,
    alveolar_bone_modulus_pa=2.0e10,
    alveolar_bone_thickness_m=2.0e-3,
    fluid_viscosity_pa_s=1.0e-3,
    pdl_permeability_m2=1.0e-13,
    tmj_joint_area_m2=5.0e-4,
    tmj_disc_thickness_m=1.5e-3,
    tmj_disc_modulus_pa=1.0e6,
    tmj_condylar_cartilage_thickness_m=1.5e-4,
    tmj_condylar_cartilage_modulus_pa=2.0e6,
    tmj_temporal_cartilage_thickness_m=1.0e-4,
    tmj_temporal_cartilage_modulus_pa=2.0e6,
    tmj_displacement_scale_m=1.0e-4,
    tmj_stiffening_coefficient=1.0,
    tmj_retrodiscal_tension_fraction=0.20,
    molar_lever_m=3.0e-2,
    condyle_lever_m=4.5e-2,
    effective_condyle_radius_m=2.0e-2,
    close_load_n=100.0,
    external_moment_nm=0.0,
    mandible_mass_kg=4.0,
    mandible_inertia_kgm2=4.0e-3,
    load_ramp_s=2.0e-2,
    simulation_time_s=2.0e-1,
    root_tolerance=1.0e-12,
    integration_max_step_s=1.0e-3,
    reference_mobility_m=8.0e-5,
)


EXPECTED_PREREG_SHA256 = "a84ba74c54f80b296bccb00c26806236509770a3c24057c4e9eefe4f196c8ca4"


PARAMETER_META: dict[str, tuple[str, str]] = {
    "pdl_fibers": ("mixed", "Three assumed fibre families; homogenisation assumptions, not measurement data"),
    "molar_support_area_m2": ("m^2", "Assumed adult molar root support surface"),
    "pdl_thickness_m": ("m", "Antagen genomsnittlig axiell PDL-tjocklek"),
    "pdl_strain_scale": ("1", "Assumed reference strain for reversible stiffening"),
    "pdl_stiffening_coefficient": ("1", "Antagen fiber-crimp/engagemangsterm"),
    "alveolar_bone_modulus_pa": ("Pa", "Assumed short-term modulus in alveolar bone"),
    "alveolar_bone_thickness_m": ("m", "Assumed load-bearing bone length in series"),
    "fluid_viscosity_pa_s": ("Pa s", "Assumed pore fluid viscosity"),
    "pdl_permeability_m2": ("m^2", "Assumed PDL permeability; controls Darcy damping"),
    "tmj_joint_area_m2": ("m^2", "Antagen effektiv TMJ-kontaktarea"),
    "tmj_disc_thickness_m": ("m", "Antagen komprimerad skivtjocklek"),
    "tmj_disc_modulus_pa": ("Pa", "Antagen skivmodul"),
    "tmj_condylar_cartilage_thickness_m": ("m", "Antagen kondylkartiltjocklek"),
    "tmj_condylar_cartilage_modulus_pa": ("Pa", "Antagen kondylkartilmodul"),
    "tmj_temporal_cartilage_thickness_m": ("m", "Antagen temporal kartiltjocklek"),
    "tmj_temporal_cartilage_modulus_pa": ("Pa", "Antagen temporal kartilmodul"),
    "tmj_displacement_scale_m": ("m", "Assumed reference displacement for TMJ stiffening"),
    "tmj_stiffening_coefficient": ("1", "Assumed disc/tissue engagement term"),
    "tmj_retrodiscal_tension_fraction": ("1", "Assumed drag ratio for open-track control"),
    "molar_lever_m": ("m", "Assumed half molar support distance from the sagittal plane"),
    "condyle_lever_m": ("m", "Assumed half condylar distance from the sagittal plane"),
    "effective_condyle_radius_m": ("m", "Assumed geometric angle proxy for condylar translation"),
    "close_load_n": ("N", "Fryst normaliserande 100 N scenario, not measurement data"),
    "external_moment_nm": ("N m", "Frozen external moment; zero in the primary symmetric run"),
    "mandible_mass_kg": ("kg", "Assumed effective jaw mass for 1D control"),
    "mandible_inertia_kgm2": ("kg m^2", "Antagen effektiv rotationsinertsi"),
    "load_ramp_s": ("s", "Frozen time profile for checking quasi-statics"),
    "simulation_time_s": ("s", "Numeric control window"),
    "root_tolerance": ("1", "Numerisk tolerans"),
    "integration_max_step_s": ("s", "Numerisk tidsstegstak"),
    "reference_mobility_m": ("m", "UNVERIFIED, FROM MEMORY: comparison anchor 0.08 mm"),
}


EQUATION_REGISTER: dict[str, str] = {
    "fiber_homogenization": "E_ax = sum_i f_i E_i cos^2(alpha_i) [Pa]",
    "pdl_solid_force": "F_PDL = (E_ax A/t) delta [1 + a_PDL (delta/(t epsilon_ref))^2] [N]",
    "bone_series_force": "F_bone = (E_bone A/L_bone) delta [N]",
    "darcy_consolidation_drag": "F_drag = (eta t A/k_perm) delta_dot [N]",
    "tooth_series_force": "F_tooth = F_PDL + F_bone + F_drag [N]",
    "tmj_series_force": "F_TMJ = F_disc(F) + F_condylar_cartilage(F) + F_temporal_cartilage(F) [N]",
    "rigid_jaw_kinematics": "u_i = y + x_i theta; u_dot_i = y_dot + x_i theta_dot [m]",
    "static_momentum_balance": "sum_i F_i = F_close; sum_i x_i F_i = M_external [N, N m]",
    "energy": "U = sum_i integral_0^u_i F_i(s) ds; D_reduced_dot = sum_tooth((F(u,v)-F(u,0))*v); this is reduced resistive power, not internal Kelvin-Voigt heat; W_ext = integral Q(t)^T q(t) dt",
    "dynamic_momentum_balance": "diag(m,I) q_ddot + J(q)^T f(J q, J q_dot) = Q(t)",
    "dimensional_checks": "k t/A = E; C E A/L = 1; c k_perm/(eta t A) = 1; F x = N m",
    "damage_scope": "No irreversible damage law; reversible nonlinear elasticity and Darcy dissipation only",
}


def parameter_table(parameters: Parameters = DEFAULT_PARAMETERS) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for entry in fields(parameters):
        value = getattr(parameters, entry.name)
        unit, source = PARAMETER_META[entry.name]
        if entry.name == "pdl_fibers":
            for family in value:
                rows.extend(
                    [
                        {
                            "name": f"fiber_{family.name}_area_fraction",
                            "value_si": family.area_fraction,
                            "unit": "1",
                            "source_assumption": source,
                        },
                        {
                            "name": f"fiber_{family.name}_modulus",
                            "value_si": family.modulus_pa,
                            "unit": "Pa",
                            "source_assumption": "Antagen homogeniserad fibermodul",
                        },
                        {
                            "name": f"fiber_{family.name}_angle",
                            "value_si": family.angle_rad,
                            "unit": "rad",
                            "source_assumption": "Antagen familjeriktning relativt tandaxeln",
                        },
                    ]
                )
        else:
            rows.append(
                {
                    "name": entry.name,
                    "value_si": value,
                    "unit": unit,
                    "source_assumption": source,
                }
            )
    return rows


def validate_parameters(parameters: Parameters = DEFAULT_PARAMETERS) -> None:
    positive_fields = (
        "molar_support_area_m2",
        "pdl_thickness_m",
        "pdl_strain_scale",
        "alveolar_bone_modulus_pa",
        "alveolar_bone_thickness_m",
        "fluid_viscosity_pa_s",
        "pdl_permeability_m2",
        "tmj_joint_area_m2",
        "tmj_disc_thickness_m",
        "tmj_disc_modulus_pa",
        "tmj_condylar_cartilage_thickness_m",
        "tmj_condylar_cartilage_modulus_pa",
        "tmj_temporal_cartilage_thickness_m",
        "tmj_temporal_cartilage_modulus_pa",
        "tmj_displacement_scale_m",
        "molar_lever_m",
        "condyle_lever_m",
        "effective_condyle_radius_m",
        "close_load_n",
        "mandible_mass_kg",
        "mandible_inertia_kgm2",
        "load_ramp_s",
        "simulation_time_s",
        "root_tolerance",
        "integration_max_step_s",
        "reference_mobility_m",
    )
    for name in positive_fields:
        if getattr(parameters, name) <= 0.0:
            raise ValueError(f"{name} must be positive")
    if not math.isclose(sum(f.area_fraction for f in parameters.pdl_fibers), 1.0, abs_tol=1e-12):
        raise ValueError("PDL fiber area fractions must sum to one")
    if any(f.area_fraction < 0.0 or f.modulus_pa <= 0.0 for f in parameters.pdl_fibers):
        raise ValueError("Invalid PDL fiber family")
    if not 0.0 <= parameters.tmj_retrodiscal_tension_fraction <= 1.0:
        raise ValueError("Invalid retrodiscal tension fraction")


def effective_pdl_modulus_pa(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    return float(
        sum(
            family.area_fraction
            * family.modulus_pa
            * math.cos(family.angle_rad) ** 2
            for family in parameters.pdl_fibers
        )
    )


def pdl_tangent_stiffness_n_per_m(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    modulus = effective_pdl_modulus_pa(parameters)
    return modulus * parameters.molar_support_area_m2 / parameters.pdl_thickness_m


def bone_stiffness_n_per_m(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    return (
        parameters.alveolar_bone_modulus_pa
        * parameters.molar_support_area_m2
        / parameters.alveolar_bone_thickness_m
    )


def hydraulic_drag_n_s_per_m(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    return (
        parameters.fluid_viscosity_pa_s
        * parameters.pdl_thickness_m
        * parameters.molar_support_area_m2
        / parameters.pdl_permeability_m2
    )


def tmj_cartilage_stiffness_n_per_m(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    return (
        parameters.tmj_joint_area_m2
        * (
            parameters.tmj_condylar_cartilage_modulus_pa
            / parameters.tmj_condylar_cartilage_thickness_m
            + parameters.tmj_temporal_cartilage_modulus_pa
            / parameters.tmj_temporal_cartilage_thickness_m
        )
    )


def tmj_disc_stiffness_n_per_m(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    return parameters.tmj_disc_modulus_pa * parameters.tmj_joint_area_m2 / parameters.tmj_disc_thickness_m


def tmj_tangent_stiffness_n_per_m(parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    return tmj_cartilage_stiffness_n_per_m(parameters) + tmj_disc_stiffness_n_per_m(parameters)


def _signed_odd_nonlinear_force(
    displacement_m: float,
    stiffness_n_per_m: float,
    scale_m: float,
    coefficient: float,
) -> float:
    magnitude = abs(displacement_m)
    force = stiffness_n_per_m * magnitude * (
        1.0 + coefficient * (magnitude / scale_m) ** 2
    )
    return math.copysign(force, displacement_m)


def pdl_elastic_force_n(
    displacement_m: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> float:
    return _signed_odd_nonlinear_force(
        displacement_m,
        pdl_tangent_stiffness_n_per_m(parameters),
        parameters.pdl_thickness_m * parameters.pdl_strain_scale,
        parameters.pdl_stiffening_coefficient,
    )


def tooth_partition_m(
    total_displacement_m: float,
    total_displacement_rate_m_per_s: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> tuple[float, float, float]:
    lower = min(0.0, total_displacement_m)
    upper = max(0.0, total_displacement_m)
    expansion = max(
        abs(total_displacement_m),
        parameters.pdl_thickness_m,
        abs(hydraulic_drag_n_s_per_m(parameters) * total_displacement_rate_m_per_s)
        / max(bone_stiffness_n_per_m(parameters), 1.0),
        1.0e-12,
    )
    lower -= expansion
    upper += expansion

    def equation(pdl_displacement: float) -> float:
        bone_displacement = total_displacement_m - pdl_displacement
        return (
            pdl_elastic_force_n(pdl_displacement, parameters)
            + hydraulic_drag_n_s_per_m(parameters) * total_displacement_rate_m_per_s
            - bone_stiffness_n_per_m(parameters) * bone_displacement
        )

    pdl_displacement = float(
        brentq(
            equation,
            lower,
            upper,
            xtol=np.nextafter(0.0, 1.0),
            rtol=4.0 * np.finfo(float).eps,
        )
    )
    bone_displacement = total_displacement_m - pdl_displacement
    force = bone_stiffness_n_per_m(parameters) * bone_displacement
    return pdl_displacement, bone_displacement, force


def tooth_path_force_n(
    displacement_m: float,
    displacement_rate_m_per_s: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> float:
    _, _, force = tooth_partition_m(
        displacement_m, displacement_rate_m_per_s, parameters
    )
    return force


def tmj_path_force_n(
    displacement_m: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> float:
    compression = _signed_odd_nonlinear_force(
        displacement_m,
        tmj_disc_stiffness_n_per_m(parameters),
        parameters.tmj_displacement_scale_m,
        parameters.tmj_stiffening_coefficient,
    )
    tension_stiffness = (
        parameters.tmj_retrodiscal_tension_fraction * tmj_tangent_stiffness_n_per_m(parameters)
    )
    tension = tmj_cartilage_stiffness_n_per_m(parameters) * displacement_m
    return compression + tension + (tension_stiffness * displacement_m if displacement_m < 0.0 else 0.0)


def branch_layout(parameters: Parameters = DEFAULT_PARAMETERS) -> tuple[tuple[str, str, float], ...]:
    return (
        ("tooth_left", "tooth", -parameters.molar_lever_m),
        ("tooth_right", "tooth", parameters.molar_lever_m),
        ("tmj_left", "tmj", -parameters.condyle_lever_m),
        ("tmj_right", "tmj", parameters.condyle_lever_m),
    )


def generalized_load(parameters: Parameters = DEFAULT_PARAMETERS) -> np.ndarray:
    return np.array([parameters.close_load_n, parameters.external_moment_nm], dtype=float)


def branch_displacements(
    generalized_coordinates: np.ndarray,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> dict[str, float]:
    vertical = float(generalized_coordinates[0])
    angle = float(generalized_coordinates[1])
    return {
        name: vertical + lever * angle
        for name, _, lever in branch_layout(parameters)
    }


def branch_forces(
    generalized_coordinates: np.ndarray,
    rates: np.ndarray | None = None,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> dict[str, float]:
    displacements = branch_displacements(generalized_coordinates, parameters)
    if rates is None:
        rates = np.zeros(2, dtype=float)
    vertical_rate = float(rates[0])
    angle_rate = float(rates[1])
    forces: dict[str, float] = {}
    for name, kind, lever in branch_layout(parameters):
        local_rate = vertical_rate + lever * angle_rate
        if kind == "tooth":
            forces[name] = tooth_path_force_n(
                displacements[name], local_rate, parameters
            )
        else:
            forces[name] = tmj_path_force_n(displacements[name], parameters)
    return forces


def elastic_energy_j(
    generalized_coordinates: np.ndarray,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> float:
    displacements = branch_displacements(generalized_coordinates, parameters)
    total = 0.0
    pdl_k = pdl_tangent_stiffness_n_per_m(parameters)
    bone_k = bone_stiffness_n_per_m(parameters)
    disc_k = tmj_disc_stiffness_n_per_m(parameters)
    cartilage_k = tmj_cartilage_stiffness_n_per_m(parameters)
    tension_k = parameters.tmj_retrodiscal_tension_fraction * (
        disc_k + cartilage_k
    )
    for name, kind, _ in branch_layout(parameters):
        displacement = displacements[name]
        if kind == "tooth":
            pdl_displacement, bone_displacement, _ = tooth_partition_m(
                displacement, 0.0, parameters
            )
            pdl_magnitude = abs(pdl_displacement)
            total += 0.5 * pdl_k * pdl_magnitude**2
            total += (
                parameters.pdl_stiffening_coefficient
                * pdl_k
                * pdl_magnitude**4
                / (
                    4.0
                    * parameters.pdl_thickness_m**2
                    * parameters.pdl_strain_scale**2
                )
            )
            total += 0.5 * bone_k * bone_displacement**2
        elif displacement >= 0.0:
            total += 0.5 * disc_k * displacement**2
            total += 0.5 * cartilage_k * displacement**2
            total += (
                parameters.tmj_stiffening_coefficient
                * disc_k
                * displacement**4
                / (4.0 * parameters.tmj_displacement_scale_m**2)
            )
        else:
            total += 0.5 * disc_k * displacement**2
            total += 0.5 * (cartilage_k + tension_k) * displacement**2
            total += (
                parameters.tmj_stiffening_coefficient
                * disc_k
                * displacement**4
                / (4.0 * parameters.tmj_displacement_scale_m**2)
            )
    return float(total)


def dissipation_power_w(
    generalized_coordinates: np.ndarray,
    rates: np.ndarray,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> float:
    """Power conjugate to the condensed force and static U(q).

    This reduced resistance is not Kelvin-Voigt heat: internal PDL strain
    is not a dynamic state. The legacy fluid_viscous output is an alias.
    """
    vertical_rate = float(rates[0])
    angle_rate = float(rates[1])
    total = 0.0
    for name, kind, lever in branch_layout(parameters):
        if kind == "tooth":
            local_rate = vertical_rate + lever * angle_rate
            local_displacement = float(generalized_coordinates[0]) + lever * float(generalized_coordinates[1])
            total += (tooth_path_force_n(local_displacement, local_rate, parameters) - tooth_path_force_n(local_displacement, 0.0, parameters)) * local_rate
    return float(total)


def equilibrium_residual(
    generalized_coordinates: np.ndarray,
    rates: np.ndarray | None = None,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> np.ndarray:
    forces = branch_forces(generalized_coordinates, rates, parameters)
    force_sum = math.fsum(forces.values())
    moment_sum = math.fsum(
        lever * forces[name] for name, _, lever in branch_layout(parameters)
    )
    return np.array([force_sum, moment_sum], dtype=float) - generalized_load(parameters)


def pdl_local_deflection_m(
    total_force_n: float,
    total_displacement_m: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> float:
    def equation(displacement: float) -> float:
        return pdl_elastic_force_n(displacement, parameters) - total_force_n

    scale = max(
        parameters.pdl_thickness_m,
        abs(total_displacement_m),
        abs(total_force_n) / max(pdl_tangent_stiffness_n_per_m(parameters), 1.0) * 2.0,
    )
    return float(
        brentq(
            equation,
            -scale,
            scale,
            xtol=np.nextafter(0.0, 1.0),
            rtol=4.0 * np.finfo(float).eps,
        )
    )


def linear_stiffness_matrix(
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> np.ndarray:
    matrix = np.zeros((2, 2), dtype=float)
    for _, kind, lever in branch_layout(parameters):
        stiffness = (
            1.0
            / (
                1.0 / pdl_tangent_stiffness_n_per_m(parameters)
                + 1.0 / bone_stiffness_n_per_m(parameters)
            )
            if kind == "tooth"
            else tmj_tangent_stiffness_n_per_m(parameters)
        )
        jacobian = np.array([1.0, lever], dtype=float)
        matrix += stiffness * np.outer(jacobian, jacobian)
    return matrix


def linear_solution(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    stiffness = linear_stiffness_matrix(parameters)
    coordinates = np.linalg.solve(stiffness, generalized_load(parameters))
    forces = branch_forces(coordinates, parameters=parameters)
    return {
        "stiffness_matrix_n_per_m": stiffness.tolist(),
        "generalized_coordinates": coordinates.tolist(),
        "forces_n": forces,
        "residual": equilibrium_residual(coordinates, parameters=parameters).tolist(),
    }


def ramped_external_work_j(
    final_coordinates: np.ndarray,
    parameters: Parameters = DEFAULT_PARAMETERS,
    sample_count: int = 201,
) -> float:
    load = generalized_load(parameters)
    fractions = np.linspace(0.0, 1.0, sample_count)
    coordinates = np.empty((sample_count, 2), dtype=float)

    for index, fraction in enumerate(fractions):
        target = fraction * load
        if index == 0:
            coordinates[index] = 0.0
            continue
        if index == sample_count - 1:
            coordinates[index] = final_coordinates
            continue
        initial = fraction * final_coordinates

        def path_residual(candidate: np.ndarray) -> np.ndarray:
            candidate_forces = branch_forces(candidate, parameters=parameters)
            return np.array(
                [
                    math.fsum(candidate_forces.values()) - target[0],
                    math.fsum(
                        lever * candidate_forces[name]
                        for name, _, lever in branch_layout(parameters)
                    )
                    - target[1],
                ],
                dtype=float,
            )

        path_solution = root(
            path_residual,
            initial,
            method="hybr",
            tol=parameters.root_tolerance,
        )
        path_error = float(np.linalg.norm(path_residual(path_solution.x), ord=np.inf))
        if path_error > 1e-7:
            raise RuntimeError("Ramped static path failed")
        coordinates[index] = path_solution.x
    derivatives = np.gradient(coordinates, fractions, axis=0, edge_order=2)
    power = np.sum((fractions[:, None] * load[None, :]) * derivatives, axis=1)
    return float(np.sum(0.5 * (power[:-1] + power[1:]) * np.diff(fractions)))


def solve_static(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    solution = root(
        lambda coordinates: equilibrium_residual(coordinates, parameters=parameters),
        np.zeros(2, dtype=float),
        method="hybr",
        tol=parameters.root_tolerance,
    )
    coordinates = np.asarray(solution.x, dtype=float)
    residual = equilibrium_residual(coordinates, parameters=parameters)
    if not solution.success and float(np.linalg.norm(residual, ord=np.inf)) > 1e-7:
        raise RuntimeError(f"Static equilibrium failed: {solution.message}")
    forces = branch_forces(coordinates, parameters=parameters)
    displacements = branch_displacements(coordinates, parameters)
    pdl_deflections = {
        "tooth_left": pdl_local_deflection_m(
            forces["tooth_left"], displacements["tooth_left"], parameters
        ),
        "tooth_right": pdl_local_deflection_m(
            forces["tooth_right"], displacements["tooth_right"], parameters
        ),
    }
    tooth_load = forces["tooth_left"] + forces["tooth_right"]
    tmj_load = forces["tmj_left"] + forces["tmj_right"]
    total_load = tooth_load + tmj_load
    external_work = ramped_external_work_j(coordinates, parameters)
    stored_energy = elastic_energy_j(coordinates, parameters)
    condyle_angles = {
        "tmj_left": math.degrees(
            math.atan2(displacements["tmj_left"], parameters.effective_condyle_radius_m)
        ),
        "tmj_right": math.degrees(
            math.atan2(displacements["tmj_right"], parameters.effective_condyle_radius_m)
        ),
    }
    pressure_scale = 1.0 / (parameters.molar_support_area_m2 * 1000.0)
    result = {
        "converged": bool(solution.success),
        "message": str(solution.message),
        "iterations": int(solution.nfev),
        "jaw_translation_um": float(coordinates[0] * 1.0e6),
        "jaw_pitch_deg": float(math.degrees(coordinates[1])),
        "branch_displacements_um": {
            name: displacements[name] * 1.0e6 for name in displacements
        },
        "pdl_local_deflections_um": {
            name: pdl_deflections[name] * 1.0e6 for name in pdl_deflections
        },
        "contact_forces_n": {
            "tooth_left": forces["tooth_left"],
            "tooth_right": forces["tooth_right"],
        },
        "tmj_forces_n": {
            "tmj_left": forces["tmj_left"],
            "tmj_right": forces["tmj_right"],
        },
        "pdl_pressure_kpa": {
            "tooth_left": forces["tooth_left"] * pressure_scale,
            "tooth_right": forces["tooth_right"] * pressure_scale,
        },
        "tmj_angle_deg": condyle_angles,
        "tooth_load_fraction": tooth_load / total_load,
        "tmj_load_fraction": tmj_load / total_load,
        "equilibrium_residual": {
            "force_n": float(residual[0]),
            "moment_nm": float(residual[1]),
            "force_relative": abs(residual[0]) / parameters.close_load_n,
            "moment_relative": abs(residual[1])
            / max(abs(parameters.external_moment_nm), parameters.close_load_n * parameters.molar_lever_m),
        },
        "external_contact_work_mj": external_work * 1000.0,
        "stored_elastic_energy_mj": stored_energy * 1000.0,
        "energy_residual_mj": (external_work - stored_energy) * 1000.0,
    }
    return result


def solve_tooth_control(
    force_n: float = 100.0,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> dict[str, float]:
    if force_n <= 0.0:
        raise ValueError("Tooth control force must be positive")

    def equation(displacement: float) -> float:
        return tooth_path_force_n(displacement, 0.0, parameters) - force_n

    stiffness = 1.0 / (
        1.0 / pdl_tangent_stiffness_n_per_m(parameters)
        + 1.0 / bone_stiffness_n_per_m(parameters)
    )
    upper = max(
        parameters.pdl_thickness_m,
        4.0 * force_n / stiffness,
        2.0 * parameters.pdl_thickness_m,
    )
    displacement = float(brentq(equation, 0.0, upper, xtol=np.nextafter(0.0, 1.0), rtol=4.0 * np.finfo(float).eps))
    pdl_displacement = pdl_local_deflection_m(force_n, displacement, parameters)
    bone_displacement = displacement - pdl_displacement
    pressure = force_n / parameters.molar_support_area_m2
    return {
        "force_n": float(force_n),
        "tooth_displacement_um": displacement * 1.0e6,
        "pdl_deflection_um": pdl_displacement * 1.0e6,
        "bone_deflection_um": bone_displacement * 1.0e6,
        "pdl_stress_kpa": pressure / 1000.0,
    }


def linear_limit_parameters(parameters: Parameters = DEFAULT_PARAMETERS) -> Parameters:
    return replace(
        parameters,
        pdl_stiffening_coefficient=0.0,
        tmj_stiffening_coefficient=0.0,
        fluid_viscosity_pa_s=1.0e-3,
    )


def all_compliance_scale_parameters(
    parameters: Parameters, scale: float
) -> Parameters:
    return replace(
        parameters,
        pdl_thickness_m=parameters.pdl_thickness_m * scale,
        alveolar_bone_thickness_m=parameters.alveolar_bone_thickness_m * scale,
        tmj_disc_thickness_m=parameters.tmj_disc_thickness_m * scale,
        tmj_condylar_cartilage_thickness_m=(
            parameters.tmj_condylar_cartilage_thickness_m * scale
        ),
        tmj_temporal_cartilage_thickness_m=(
            parameters.tmj_temporal_cartilage_thickness_m * scale
        ),
        tmj_displacement_scale_m=parameters.tmj_displacement_scale_m * scale,
    )


def perturb_parameters(
    parameter: str, scale: float, parameters: Parameters = DEFAULT_PARAMETERS
) -> Parameters:
    if parameter == "E_PDL":
        fibers = tuple(
            FiberFamily(
                family.name,
                family.area_fraction,
                family.modulus_pa * scale,
                family.angle_rad,
            )
            for family in parameters.pdl_fibers
        )
        return replace(parameters, pdl_fibers=fibers)
    if parameter == "t_PDL":
        return replace(parameters, pdl_thickness_m=parameters.pdl_thickness_m * scale)
    if parameter == "E_TMJ_disc":
        return replace(
            parameters, tmj_disc_modulus_pa=parameters.tmj_disc_modulus_pa * scale
        )
    raise KeyError(parameter)


def _smoothstep(value: float) -> Any:
    clamped = np.clip(value, 0.0, 1.0)
    return clamped * clamped * (3.0 - 2.0 * clamped)


def transient_summary(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    def time_scaling(time_s: float) -> float:
        return _smoothstep(time_s / parameters.load_ramp_s)

    def rhs(time_s: float, state: np.ndarray) -> np.ndarray:
        coordinates = state[:2]
        rates = state[2:4]
        load = parameters.close_load_n * time_scaling(time_s)
        moment = parameters.external_moment_nm * time_scaling(time_s)
        forces = branch_forces(coordinates, rates, parameters)
        support_force = math.fsum(forces.values())
        support_moment = math.fsum(
            lever * forces[name] for name, _, lever in branch_layout(parameters)
        )
        acceleration = np.array(
            [
                (load - support_force) / parameters.mandible_mass_kg,
                (moment - support_moment) / parameters.mandible_inertia_kgm2,
            ],
            dtype=float,
        )
        return np.concatenate((rates, acceleration))

    sample_count = 201
    sample_times = np.linspace(0.0, parameters.simulation_time_s, sample_count)
    solution = solve_ivp(
        rhs,
        (0.0, parameters.simulation_time_s),
        np.zeros(4, dtype=float),
        method="BDF",
        t_eval=sample_times,
        max_step=parameters.integration_max_step_s,
        rtol=1e-9,
        atol=1e-11,
    )
    if not solution.success:
        raise RuntimeError(f"Transient failed: {solution.message}")
    coordinates = solution.y[:2].T
    rates = solution.y[2:4].T
    accelerations = np.empty_like(coordinates)
    dissipations = np.empty(solution.t.size, dtype=float)
    for index, time_s in enumerate(solution.t):
        accelerations[index] = rhs(time_s, solution.y[:, index])[2:4]
        dissipations[index] = dissipation_power_w(coordinates[index], rates[index], parameters)
    force = parameters.close_load_n * time_scaling(solution.t)
    displacement_power = force * rates[:, 0] + parameters.external_moment_nm * (
        time_scaling(solution.t) * rates[:, 1]
    )
    external_work = float(
        np.sum(0.5 * (displacement_power[:-1] + displacement_power[1:]) * np.diff(solution.t))
    )
    dissipated = float(
        np.sum(0.5 * (dissipations[:-1] + dissipations[1:]) * np.diff(solution.t))
    )
    post_ramp = solution.t >= parameters.load_ramp_s
    inertial_force = parameters.mandible_mass_kg * np.abs(accelerations[post_ramp, 0])
    rotational_inertial_moment = parameters.mandible_inertia_kgm2 * np.abs(
        accelerations[post_ramp, 1]
    )
    peak_scale = max(
        abs(parameters.close_load_n),
        abs(parameters.external_moment_nm) / max(parameters.molar_lever_m, 1e-12),
    )
    final_forces = branch_forces(coordinates[-1], rates[-1], parameters)
    final_displacements = branch_displacements(coordinates[-1], parameters)
    return {
        "converged": True,
        "message": solution.message,
        "sample_count": int(solution.t.size),
        "final_time_s": float(solution.t[-1]),
        "final_jaw_translation_um": float(coordinates[-1, 0] * 1.0e6),
        "final_jaw_pitch_deg": float(math.degrees(coordinates[-1, 1])),
        "final_tooth_forces_n": {
            "tooth_left": final_forces["tooth_left"],
            "tooth_right": final_forces["tooth_right"],
        },
        "final_tmj_forces_n": {
            "tmj_left": final_forces["tmj_left"],
            "tmj_right": final_forces["tmj_right"],
        },
        "final_tooth_displacements_um": {
            "tooth_left": final_displacements["tooth_left"] * 1.0e6,
            "tooth_right": final_displacements["tooth_right"] * 1.0e6,
        },
        "max_abs_inertial_force_n": float(np.max(inertial_force)),
        "max_abs_rotational_inertial_moment_nm": float(
            np.max(rotational_inertial_moment)
        ),
        "max_post_ramp_inertial_to_load": float(
            max(
                np.max(inertial_force) / parameters.close_load_n,
                np.max(rotational_inertial_moment) / peak_scale,
            )
        ),
        "external_contact_work_mj": external_work * 1000.0,
        "reduced_resistive_dissipation_mj": dissipated * 1000.0,
        "fluid_viscous_dissipation_mj": dissipated * 1000.0,  # legacy alias
        "dissipation_interpretation": "Reduced resistive work relative to static U(q); not internal Kelvin-Voigt heat. fluid_viscous_dissipation_mj is a legacy alias.",
    }


def sensitivity_analysis(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    baseline = solve_static(parameters)
    metrics = {
        "tooth_load_fraction": baseline["tooth_load_fraction"],
        "tmj_load_fraction": baseline["tmj_load_fraction"],
        "tooth_left_displacement_um": baseline["branch_displacements_um"]["tooth_left"],
        "pdl_pressure_left_kpa": baseline["pdl_pressure_kpa"]["tooth_left"],
        "tmj_left_load_n": baseline["tmj_forces_n"]["tmj_left"],
        "tmj_left_angle_deg": baseline["tmj_angle_deg"]["tmj_left"],
    }
    parameters_to_test = ("E_PDL", "t_PDL", "E_TMJ_disc")
    cases: list[dict[str, Any]] = []
    for name in parameters_to_test:
        relative_changes: list[float] = []
        perturbation_rows: list[dict[str, Any]] = []
        for factor in (0.5, 1.5):
            result = solve_static(perturb_parameters(name, factor, parameters))
            changed_metrics = {
                "tooth_load_fraction": result["tooth_load_fraction"],
                "tmj_load_fraction": result["tmj_load_fraction"],
                "tooth_left_displacement_um": result["branch_displacements_um"]["tooth_left"],
                "pdl_pressure_left_kpa": result["pdl_pressure_kpa"]["tooth_left"],
                "tmj_left_load_n": result["tmj_forces_n"]["tmj_left"],
                "tmj_left_angle_deg": result["tmj_angle_deg"]["tmj_left"],
            }
            changes = {
                metric: 100.0 * (changed_metrics[metric] / baseline_value - 1.0)
                for metric, baseline_value in metrics.items()
            }
            relative_changes.extend(abs(value) for value in changes.values())
            perturbation_rows.append(
                {
                    "factor": factor,
                    "relative_change_percent": changes,
                    "result": result,
                }
            )
        cases.append(
            {
                "parameter": name,
                "perturbations": perturbation_rows,
                "max_abs_relative_change_percent": max(relative_changes),
                "mean_abs_relative_change_percent": float(
                    np.mean(relative_changes)
                ),
            }
        )
    cases.sort(key=lambda item: item["max_abs_relative_change_percent"], reverse=True)
    return {
        "baseline": baseline,
        "ranked_parameters": [
            {
                "rank": index + 1,
                "parameter": item["parameter"],
                "max_abs_relative_change_percent": item["max_abs_relative_change_percent"],
                "mean_abs_relative_change_percent": item["mean_abs_relative_change_percent"],
            }
            for index, item in enumerate(cases)
        ],
        "cases": cases,
    }


def dimensional_checks(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, float]:
    pdl_k = pdl_tangent_stiffness_n_per_m(parameters)
    bone_c = parameters.alveolar_bone_thickness_m / (
        parameters.alveolar_bone_modulus_pa * parameters.molar_support_area_m2
    )
    drag_c = hydraulic_drag_n_s_per_m(parameters)
    disc_c = parameters.tmj_disc_thickness_m / (
        parameters.tmj_disc_modulus_pa * parameters.tmj_joint_area_m2
    )
    condylar_c = parameters.tmj_condylar_cartilage_thickness_m / (
        parameters.tmj_condylar_cartilage_modulus_pa * parameters.tmj_joint_area_m2
    )
    temporal_c = parameters.tmj_temporal_cartilage_thickness_m / (
        parameters.tmj_temporal_cartilage_modulus_pa * parameters.tmj_joint_area_m2
    )
    return {
        "pdl_k_times_t_over_A_recovered_E_ratio": pdl_k
        * parameters.pdl_thickness_m
        / parameters.molar_support_area_m2
        / effective_pdl_modulus_pa(parameters),
        "bone_C_times_E_times_A_over_L": bone_c
        * parameters.alveolar_bone_modulus_pa
        * parameters.molar_support_area_m2
        / parameters.alveolar_bone_thickness_m,
        "drag_C_times_k_over_eta_t_A": drag_c
        * parameters.pdl_permeability_m2
        / (
            parameters.fluid_viscosity_pa_s
            * parameters.pdl_thickness_m
            * parameters.molar_support_area_m2
        ),
        "disc_C_times_E_times_A_over_t": disc_c
        * parameters.tmj_disc_modulus_pa
        * parameters.tmj_joint_area_m2
        / parameters.tmj_disc_thickness_m,
        "condylar_C_times_E_times_A_over_t": condylar_c
        * parameters.tmj_condylar_cartilage_modulus_pa
        * parameters.tmj_joint_area_m2
        / parameters.tmj_condylar_cartilage_thickness_m,
        "temporal_C_times_E_times_A_over_t": temporal_c
        * parameters.tmj_temporal_cartilage_modulus_pa
        * parameters.tmj_joint_area_m2
        / parameters.tmj_temporal_cartilage_thickness_m,
        "force_times_lever_recovered_moment": 1.0
        * parameters.condyle_lever_m
        / parameters.condyle_lever_m,
    }


def verify_preregistration(path: Path | None = None) -> dict[str, Any]:
    prereg_path = path or Path(__file__).with_name("PREREG.md")
    checksum = sha256(prereg_path.read_bytes()).hexdigest()
    return {
        "path": str(prereg_path),
        "expected_sha256": EXPECTED_PREREG_SHA256,
        "actual_sha256": checksum,
        "matches": checksum == EXPECTED_PREREG_SHA256,
    }


def run_analysis(output_path: Path | None = None) -> dict[str, Any]:
    validate_parameters()
    preregistration = verify_preregistration()
    if not preregistration["matches"]:
        raise RuntimeError("PREREG.md changed after checksum freeze")
    baseline = solve_static(DEFAULT_PARAMETERS)
    control = solve_tooth_control(DEFAULT_PARAMETERS.close_load_n, DEFAULT_PARAMETERS)
    linear = linear_solution(linear_limit_parameters(DEFAULT_PARAMETERS))
    linear_static = solve_static(linear_limit_parameters(DEFAULT_PARAMETERS))
    sensitivity = sensitivity_analysis(DEFAULT_PARAMETERS)
    transient = transient_summary(DEFAULT_PARAMETERS)
    dimensions = dimensional_checks(DEFAULT_PARAMETERS)
    scaled = solve_static(all_compliance_scale_parameters(DEFAULT_PARAMETERS, 2.0))
    fraction_difference = max(
        abs(scaled["tooth_load_fraction"] - baseline["tooth_load_fraction"]),
        abs(scaled["tmj_load_fraction"] - baseline["tmj_load_fraction"]),
    )
    reference_ratio = (
        control["tooth_displacement_um"] * 1.0e-6
        / DEFAULT_PARAMETERS.reference_mobility_m
    )
    proxy_status = "PROXY_PASS" if 0.5 <= reference_ratio <= 2.0 else "PROXY_FAIL"
    equilibrium_ok = max(
        baseline["equilibrium_residual"]["force_relative"],
        baseline["equilibrium_residual"]["moment_relative"],
    ) < 1.0e-8
    linear_relative_error = max(
        abs(linear["generalized_coordinates"][0] - linear_static["jaw_translation_um"] * 1.0e-6)
        / max(abs(linear["generalized_coordinates"][0]), 1.0e-30),
        abs(linear["generalized_coordinates"][1] - math.radians(linear_static["jaw_pitch_deg"]))
        / max(abs(linear["generalized_coordinates"][1]), 1.0e-30),
    )
    analytic_ok = linear_relative_error < 1.0e-10
    mechanism_ok = (
        fraction_difference < 1.0e-10
        and perturb_parameters("t_PDL", 1.5).pdl_thickness_m
        > DEFAULT_PARAMETERS.pdl_thickness_m
        and solve_static(perturb_parameters("t_PDL", 1.5))["tooth_load_fraction"]
        < baseline["tooth_load_fraction"]
    )
    dimension_ok = all(abs(value - 1.0) < 1.0e-12 for value in dimensions.values())
    transient_ok = transient["max_post_ramp_inertial_to_load"] < 0.05
    result: dict[str, Any] = {
        "id": "BT-HX-Q090",
        "model": "rigid-jaw 2DOF nonlinear compliant support network",
        "equations": EQUATION_REGISTER,
        "preregistration": preregistration,
        "parameter_table": parameter_table(),
        "derived_baseline_parameters": {
            "effective_pdl_axial_modulus_pa": effective_pdl_modulus_pa(),
            "pdl_tangent_stiffness_n_per_m": pdl_tangent_stiffness_n_per_m(),
            "bone_stiffness_n_per_m": bone_stiffness_n_per_m(),
            "hydraulic_drag_n_s_per_m": hydraulic_drag_n_s_per_m(),
            "tmj_tangent_stiffness_n_per_m": tmj_tangent_stiffness_n_per_m(),
        },
        "static_baseline": baseline,
        "single_tooth_100n_control": control,
        "reference_comparison": {
            "reference_mobility_m": DEFAULT_PARAMETERS.reference_mobility_m,
            "reference_value": "0.08 mm",
            "reference_force": "nominal 100 N",
            "reference_source": "Bien & Topp (1970), Changes in the periodicity of the tooth mobility pattern during mastication, Journal of Periodology",
            "reference_source_status": "OVERIFIERAD, UR MINNET",
            "model_to_reference_ratio": reference_ratio,
            "frozen_numeric_interval": [0.5, 2.0],
            "numeric_proxy_status": proxy_status,
            "external_validation_status": "UNKNOWN",
            "reason": "No internet verification and no verified matching force protocol",
        },
        "linear_analytical_limit": {
            "matrix_solution": linear,
            "nonlinear_solver_in_linear_limit": linear_static,
            "relative_error": linear_relative_error,
            "status": "PASS" if analytic_ok else "FAIL",
        },
        "sensitivity": sensitivity,
        "transient_check": transient,
        "common_compliance_placebo": {
            "scale": 2.0,
            "scaled_solution": scaled,
            "max_fraction_difference": fraction_difference,
            "status": "PASS" if fraction_difference < 1.0e-10 else "FAIL",
        },
        "dimensional_checks": {
            "dimensionless_identities": dimensions,
            "status": "PASS" if dimension_ok else "FAIL",
        },
        "criteria": {
            "reference_proxy": proxy_status,
            "external_validation": "UNKNOWN",
            "equilibrium": "PASS" if equilibrium_ok else "FAIL",
            "mechanism": "PASS" if mechanism_ok else "FAIL",
            "analytical_limit": "PASS" if analytic_ok else "FAIL",
            "quasi_static_post_ramp": "PASS" if transient_ok else "FAIL",
            "overall_external_claim": "UNKNOWN",
        },
        "limitations": [
            "No individual geometry, microscopy, simultaneous force-motion data, or damage threshold was available",
            "The literature comparator is explicitly unverified and may not share the modeled load protocol",
            "The first model is sagittal and symmetric; fracture, contact loss, and muscle recruitment are not resolved",
            "The TMJ angle is a geometric translation-to-angle proxy, not a measured joint angle",
        ],
    }
    destination = output_path or Path(__file__).with_name("results.json")
    destination.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    analysis = run_analysis()
    print(
        json.dumps(
            {
                "id": analysis["id"],
                "tooth_load_fraction": analysis["static_baseline"]["tooth_load_fraction"],
                "single_tooth_100n_displacement_um": analysis["single_tooth_100n_control"]["tooth_displacement_um"],
                "criteria": analysis["criteria"],
            },
            ensure_ascii=False,
        )
    )
