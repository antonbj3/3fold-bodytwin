from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import expm_multiply, spsolve


EQUATIONS = {
    "solvent": "D_alpha(x) = D_alpha0 / (1 + k_visc (1 / V_r - 1))",
    "accessible": "chi_i(x) = 1 - phi_i(x)",
    "tortuosity": "tau_i(x) = 1 + 0.5 phi_i(x) / chi_i(x)",
    "permeability": "K_nano = a^2 chi_n^3 / (2 phi_n^2 + chi_n phi_n (1 + interaction_strength phi_n))",
    "hydrodynamic_hindrance": "H(x) = 1 + R / sqrt(K_nano) + (R / sqrt(K_nano))^2 / 9",
    "effective_diffusivity": "D(x) = D_alpha(x) (chi_m / tau_m) (chi_n / tau_n) / H(x)",
    "binding": "free -> bound: k_on(x) = k_on0 (1 + interaction_strength phi_n(x)) chi_n(x); bound -> free: k_off(x) = k_off0",
    "diffusion_equation": "partial_t p = partial_x (D(x) partial_x p) - k_on p_free + k_off p_bound",
    "first_passage": "-1 = sum_j Q_ij T_j; T_exit = 0",
    "adjusted_comparator": "D_adj = spatial mean of D(x) on transient cells",
}


@dataclass(frozen=True)
class Parameters:
    d_alpha0: float = 20.0
    d_water: float = 87.0
    length_um: float = 10.0
    n_cells: int = 161
    tracer_radius_m: float = 2.3e-9
    pore_radius_m: float = 4.0e-9
    phi_m: float = 0.08
    phi_n_mean: float = 0.18
    phi_n_amplitude: float = 0.08
    relative_volume: float = 1.0
    k_visc: float = 1.0
    tortuosity_coefficient: float = 0.5
    interaction_strength: float = 4.0
    k_on0: float = 0.12
    k_off0: float = 0.03
    observation_time_s: float = 10.0

    def validate(self) -> None:
        scalar_values = (
            self.d_alpha0,
            self.d_water,
            self.length_um,
            self.tracer_radius_m,
            self.pore_radius_m,
            self.phi_m,
            self.phi_n_mean,
            self.phi_n_amplitude,
            self.relative_volume,
            self.k_visc,
            self.tortuosity_coefficient,
            self.interaction_strength,
            self.k_on0,
            self.k_off0,
            self.observation_time_s,
        )
        if not all(np.isfinite(value) for value in scalar_values):
            raise ValueError("all parameters must be finite")
        if not isinstance(self.n_cells, (int, np.integer)) or self.n_cells < 3:
            raise ValueError("n_cells must be an integer at least 3")
        if self.d_alpha0 <= 0 or self.d_water <= 0 or self.length_um <= 0:
            raise ValueError("diffusivities and length must be positive")
        if self.tracer_radius_m < 0 or self.pore_radius_m <= 0:
            raise ValueError("radii must satisfy 0 <= tracer radius and pore radius > 0")
        if not 0 <= self.phi_m < 1:
            raise ValueError("phi_m must be in [0, 1)")
        if self.phi_n_amplitude < 0:
            raise ValueError("phi_n_amplitude must be nonnegative")
        phi_n_low = self.phi_n_mean - self.phi_n_amplitude
        phi_n_high = self.phi_n_mean + self.phi_n_amplitude
        if phi_n_low < 0 or phi_n_high >= 1:
            raise ValueError("nano excluded fraction must remain in [0, 1)")
        if self.relative_volume <= 0:
            raise ValueError("relative_volume must be positive")
        viscosity_denominator = 1 + self.k_visc * (1 / self.relative_volume - 1)
        if viscosity_denominator <= 0:
            raise ValueError("viscosity denominator must be positive")
        if self.tortuosity_coefficient < 0 or self.interaction_strength < 0:
            raise ValueError("dimensionless interaction coefficients must be nonnegative")
        if self.k_on0 < 0 or self.k_off0 < 0:
            raise ValueError("binding rates must be nonnegative")
        if self.observation_time_s <= 0:
            raise ValueError("observation_time_s must be positive")


@dataclass(frozen=True)
class LocalCoefficients:
    x_um: np.ndarray
    phi_m: np.ndarray
    phi_n: np.ndarray
    chi_m: np.ndarray
    chi_n: np.ndarray
    tau_m: np.ndarray
    tau_n: np.ndarray
    k_nano_m2: np.ndarray
    hindrance: np.ndarray
    d_alpha_um2_per_s: np.ndarray
    d_free_um2_per_s: np.ndarray
    k_on_per_s: np.ndarray
    k_off_per_s: np.ndarray


def coordinates(parameters: Parameters) -> np.ndarray:
    return np.linspace(0.0, parameters.length_um, parameters.n_cells)


def local_coefficients(parameters: Parameters) -> LocalCoefficients:
    parameters.validate()
    x_um = coordinates(parameters)
    phi_m = np.full(parameters.n_cells, parameters.phi_m, dtype=float)
    phi_n = parameters.phi_n_mean + parameters.phi_n_amplitude * np.cos(
        2.0 * np.pi * x_um / parameters.length_um
    )
    chi_m = 1.0 - phi_m
    chi_n = 1.0 - phi_n
    tau_m = 1.0 + parameters.tortuosity_coefficient * phi_m / chi_m
    tau_n = 1.0 + parameters.tortuosity_coefficient * phi_n / chi_n
    denominator = 2.0 * phi_n**2 + chi_n * phi_n * (
        1.0 + parameters.interaction_strength * phi_n
    )
    k_nano_m2 = np.full(parameters.n_cells, np.inf, dtype=float)
    positive_obstacle = phi_n > 0
    k_nano_m2[positive_obstacle] = (
        parameters.pore_radius_m**2
        * chi_n[positive_obstacle] ** 3
        / denominator[positive_obstacle]
    )
    particle_pore_ratio = np.zeros(parameters.n_cells, dtype=float)
    finite_k = np.isfinite(k_nano_m2)
    particle_pore_ratio[finite_k] = parameters.tracer_radius_m / np.sqrt(
        k_nano_m2[finite_k]
    )
    hindrance = 1.0 + particle_pore_ratio + particle_pore_ratio**2 / 9.0
    d_alpha = parameters.d_alpha0 / (
        1.0 + parameters.k_visc * (1.0 / parameters.relative_volume - 1.0)
    )
    d_alpha_um2_per_s = np.full(parameters.n_cells, d_alpha, dtype=float)
    d_free_um2_per_s = (
        d_alpha_um2_per_s
        * (chi_m / tau_m)
        * (chi_n / tau_n)
        / hindrance
    )
    k_on_per_s = (
        parameters.k_on0
        * (1.0 + parameters.interaction_strength * phi_n)
        * chi_n
    )
    k_off_per_s = np.full(parameters.n_cells, parameters.k_off0, dtype=float)
    if np.any(k_on_per_s > 0) != np.any(k_off_per_s > 0):
        raise ValueError("on and off binding rates must be enabled together")
    return LocalCoefficients(
        x_um=x_um,
        phi_m=phi_m,
        phi_n=phi_n,
        chi_m=chi_m,
        chi_n=chi_n,
        tau_m=tau_m,
        tau_n=tau_n,
        k_nano_m2=k_nano_m2,
        hindrance=hindrance,
        d_alpha_um2_per_s=d_alpha_um2_per_s,
        d_free_um2_per_s=d_free_um2_per_s,
        k_on_per_s=k_on_per_s,
        k_off_per_s=k_off_per_s,
    )


def transient_generator(
    d_free_um2_per_s: np.ndarray,
    k_on_per_s: np.ndarray,
    k_off_per_s: np.ndarray,
    dx_um: float,
) -> tuple[csr_matrix, bool, int]:
    n_total = len(d_free_um2_per_s)
    if n_total < 2:
        raise ValueError("at least two spatial nodes are required")
    if len(k_on_per_s) != n_total or len(k_off_per_s) != n_total:
        raise ValueError("coefficient arrays must have equal length")
    n = n_total - 1
    if dx_um <= 0 or np.any(d_free_um2_per_s <= 0):
        raise ValueError("diffusivity and grid spacing must be positive")
    if np.any(k_on_per_s < 0) or np.any(k_off_per_s < 0):
        raise ValueError("binding rates must be nonnegative")
    if np.any(k_on_per_s[:n] > 0) != np.any(k_off_per_s[:n] > 0):
        raise ValueError("on and off binding rates must be enabled together")
    has_binding = bool(np.any(k_on_per_s[:n] > 0))
    transient_size = 2 * n if has_binding else n
    # Fick concentration per fixed geometric volume; half-width reflecting cell.
    # Shared harmonic face diffusivity gives reciprocal conductance.
    face_diffusivity = 2.0 * d_free_um2_per_s[:-1] * d_free_um2_per_s[1:] / (
        d_free_um2_per_s[:-1] + d_free_um2_per_s[1:]
    )
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for i in range(n):
        width_weight = 0.5 if i == 0 else 1.0
        left_rate = float(face_diffusivity[i - 1]) / (width_weight * dx_um**2) if i > 0 else 0.0
        right_rate = float(face_diffusivity[i]) / (width_weight * dx_um**2)
        bind_rate = float(k_on_per_s[i]) if has_binding else 0.0
        rows.append(i)
        columns.append(i)
        values.append(-(left_rate + right_rate + bind_rate))
        if i > 0:
            rows.append(i)
            columns.append(i - 1)
            values.append(left_rate)
        if i < n - 1:
            rows.append(i)
            columns.append(i + 1)
            values.append(right_rate)
        if has_binding:
            rows.append(i)
            columns.append(n + i)
            values.append(bind_rate)
    if has_binding:
        for i in range(n):
            off_rate = float(k_off_per_s[i])
            rows.append(n + i)
            columns.append(n + i)
            values.append(-off_rate)
            rows.append(n + i)
            columns.append(i)
            values.append(off_rate)
    matrix = coo_matrix(
        (np.asarray(values), (np.asarray(rows), np.asarray(columns))),
        shape=(transient_size, transient_size),
        dtype=float,
    ).tocsr()
    matrix.sum_duplicates()
    return matrix, has_binding, n


def expected_passage_time(
    d_free_um2_per_s: np.ndarray,
    k_on_per_s: np.ndarray,
    k_off_per_s: np.ndarray,
    dx_um: float,
) -> float:
    generator, _, _ = transient_generator(
        d_free_um2_per_s, k_on_per_s, k_off_per_s, dx_um
    )
    solution = spsolve(-generator.tocsc(), np.ones(generator.shape[0], dtype=float))
    if np.any(~np.isfinite(solution)) or solution[0] <= 0:
        raise RuntimeError("first-passage linear solve did not converge")
    return float(solution[0])


def propagated_distribution(
    d_free_um2_per_s: np.ndarray,
    k_on_per_s: np.ndarray,
    k_off_per_s: np.ndarray,
    dx_um: float,
    time_s: float,
) -> tuple[np.ndarray, np.ndarray, bool, float]:
    if time_s <= 0:
        raise ValueError("time_s must be positive")
    generator, has_binding, n = transient_generator(
        d_free_um2_per_s, k_on_per_s, k_off_per_s, dx_um
    )
    initial = np.zeros(generator.shape[0], dtype=float)
    initial[0] = 1.0
    evolved = expm_multiply((generator.T * time_s).tocsr(), initial)
    evolved = np.asarray(evolved, dtype=float)
    if np.any(evolved < -1e-8):
        raise RuntimeError("propagated distribution has a material negative probability")
    evolved = np.maximum(evolved, 0.0)
    mass = float(np.sum(evolved))
    if not np.isfinite(mass) or mass <= 0:
        raise RuntimeError("propagated distribution has no surviving probability")
    free_profile = evolved[:n]
    bound_profile = evolved[n:] if has_binding else np.zeros(0, dtype=float)
    survival_probability = mass
    return free_profile, bound_profile, has_binding, survival_probability


def simulate_transport(
    d_free_um2_per_s: np.ndarray,
    k_on_per_s: np.ndarray,
    k_off_per_s: np.ndarray,
    parameters: Parameters,
) -> dict:
    parameters.validate()
    dx_um = parameters.length_um / (parameters.n_cells - 1)
    passage_time_s = expected_passage_time(
        d_free_um2_per_s, k_on_per_s, k_off_per_s, dx_um
    )
    free_profile, bound_profile, has_binding, survival_probability = propagated_distribution(
        d_free_um2_per_s,
        k_on_per_s,
        k_off_per_s,
        dx_um,
        parameters.observation_time_s,
    )
    if bound_profile.size:
        bound_mean = float(np.mean(bound_profile))
        bound_std = float(np.std(bound_profile))
        bound_cv = float(bound_std / bound_mean) if bound_mean > 0 else 0.0
        bound_range = float(np.max(bound_profile) - np.min(bound_profile))
    else:
        bound_mean = 0.0
        bound_cv = 0.0
        bound_range = 0.0
    return {
        "passage_time_s": passage_time_s,
        "free_profile": free_profile.tolist(),
        "bound_profile": bound_profile.tolist(),
        "has_binding": has_binding,
        "free_probability": float(np.sum(free_profile)),
        "bound_probability": float(np.sum(bound_profile)),
        "survival_probability": float(survival_probability),
        "absorption_probability": float(1.0 - survival_probability),
        "bound_occupancy_total": float(np.sum(bound_profile)),
        "bound_occupancy_mean": bound_mean,
        "bound_occupancy_cv": bound_cv,
        "bound_occupancy_range": bound_range,
        "d_free_mean_um2_per_s": float(np.mean(d_free_um2_per_s[:-1])),
        "d_free_min_um2_per_s": float(np.min(d_free_um2_per_s[:-1])),
        "d_free_max_um2_per_s": float(np.max(d_free_um2_per_s[:-1])),
    }


def compare_models(parameters: Parameters) -> dict:
    coefficients = local_coefficients(parameters)
    d_adjusted = float(np.mean(coefficients.d_free_um2_per_s[:-1]))
    adjusted_d = np.full(parameters.n_cells, d_adjusted, dtype=float)
    mechanistic = simulate_transport(
        coefficients.d_free_um2_per_s,
        coefficients.k_on_per_s,
        coefficients.k_off_per_s,
        parameters,
    )
    adjusted = simulate_transport(
        adjusted_d,
        coefficients.k_on_per_s,
        coefficients.k_off_per_s,
        parameters,
    )
    passage_ratio = mechanistic["passage_time_s"] / adjusted["passage_time_s"] - 1.0
    bound_l1 = float(
        np.sum(
            np.abs(
                np.asarray(mechanistic["bound_profile"], dtype=float)
                - np.asarray(adjusted["bound_profile"], dtype=float)
            )
        )
    )
    free_l1 = float(
        np.sum(
            np.abs(
                np.asarray(mechanistic["free_profile"], dtype=float)
                - np.asarray(adjusted["free_profile"], dtype=float)
            )
        )
    )
    return {
        "x_um": coefficients.x_um.tolist(),
        "x_profile_um": coefficients.x_um[:-1].tolist(),
        "phi_m": coefficients.phi_m.tolist(),
        "phi_n": coefficients.phi_n.tolist(),
        "d_mechanistic_um2_per_s": coefficients.d_free_um2_per_s.tolist(),
        "d_adjusted_um2_per_s": d_adjusted,
        "mechanistic_d_to_d_alpha0": float(
            np.mean(coefficients.d_free_um2_per_s[:-1]) / parameters.d_alpha0
        ),
        "adjusted_d_to_d_alpha0": float(d_adjusted / parameters.d_alpha0),
        "hindrance": coefficients.hindrance.tolist(),
        "k_on_per_s": coefficients.k_on_per_s.tolist(),
        "mechanistic": mechanistic,
        "adjusted": adjusted,
        "passage_time_ratio": float(passage_ratio),
        "passage_time_ratio_percent": float(100.0 * passage_ratio),
        "bound_profile_l1_difference": bound_l1,
        "free_profile_l1_difference": free_l1,
    }


def analytical_sanity(parameters: Parameters) -> dict:
    empty = replace(
        parameters,
        n_cells=max(3, min(parameters.n_cells, 81)),
        phi_m=0.0,
        phi_n_mean=0.0,
        phi_n_amplitude=0.0,
        interaction_strength=0.0,
        k_on0=0.0,
        k_off0=0.0,
    )
    coefficients = local_coefficients(empty)
    dx_um = empty.length_um / (empty.n_cells - 1)
    measured = expected_passage_time(
        coefficients.d_free_um2_per_s,
        coefficients.k_on_per_s,
        coefficients.k_off_per_s,
        dx_um,
    )
    expected = empty.length_um**2 / (2.0 * empty.d_alpha0)
    relative_error = abs(measured - expected) / expected
    mean_d = float(np.mean(coefficients.d_free_um2_per_s))
    return {
        "empty_space_discrete_passage_s": measured,
        "empty_space_analytic_passage_s": expected,
        "relative_error": float(relative_error),
        "relative_error_threshold": 0.02,
        "mean_empty_space_diffusivity_um2_per_s": mean_d,
        "input_diffusivity_um2_per_s": empty.d_alpha0,
        "diffusivity_absolute_error_um2_per_s": abs(mean_d - empty.d_alpha0),
        "passed": bool(
            relative_error <= 0.02
            and abs(mean_d - empty.d_alpha0) <= 1e-12
        ),
    }


def unit_check(parameters: Parameters) -> dict:
    coefficients = local_coefficients(parameters)
    finite = bool(
        np.all(np.isfinite(coefficients.d_free_um2_per_s))
        and np.all(coefficients.d_free_um2_per_s > 0)
        and np.all(np.isfinite(coefficients.k_on_per_s))
        and np.all(np.isfinite(coefficients.k_off_per_s))
        and np.all(coefficients.hindrance >= 1.0)
    )
    return {
        "coordinate": "um",
        "diffusivity": "um^2 s^-1",
        "length_for_diffusion_rates": "um",
        "binding_rate": "s^-1",
        "permeability": "m^2",
        "particle_radius": "m",
        "pore_radius": "m",
        "hindrance_ratio": "dimensionless",
        "passage_time": "s",
        "finite_positive_fields": finite,
    }


def parameter_table(parameters: Parameters) -> list[dict]:
    return [
        {
            "parameter": "D_alpha0",
            "value": parameters.d_alpha0,
            "unit": "um^2 s^-1",
            "role": "free iso-osmotic tracer diffusivity",
            "source_or_assumption": "assumption; context is source free-GFP range 16–20 um^2 s^-1",
        },
        {
            "parameter": "D_water",
            "value": parameters.d_water,
            "unit": "um^2 s^-1",
            "role": "water reference, not used in the effective law",
            "source_or_assumption": "Destrian et al. PNAS Fig. 3 discussion",
        },
        {
            "parameter": "L",
            "value": parameters.length_um,
            "unit": "um",
            "role": "transport compartment length",
            "source_or_assumption": "bounded model assumption",
        },
        {
            "parameter": "N",
            "value": parameters.n_cells,
            "unit": "cells",
            "role": "finite-volume spatial resolution",
            "source_or_assumption": "frozen numerical choice",
        },
        {
            "parameter": "t_obs",
            "value": parameters.observation_time_s,
            "unit": "s",
            "role": "late-time spatial-distribution readout",
            "source_or_assumption": "frozen numerical choice",
        },
        {
            "parameter": "R",
            "value": parameters.tracer_radius_m,
            "unit": "m",
            "role": "finite tracer radius in hydrodynamic hindrance",
            "source_or_assumption": "Destrian et al. PNAS text, GFP radius",
        },
        {
            "parameter": "a",
            "value": parameters.pore_radius_m,
            "unit": "m",
            "role": "characteristic nano-pore radius",
            "source_or_assumption": "frozen coarse-grained assumption",
        },
        {
            "parameter": "phi_m",
            "value": parameters.phi_m,
            "unit": "dimensionless",
            "role": "micro-obstacle excluded volume fraction",
            "source_or_assumption": "Destrian et al. PNAS model comparison choice",
        },
        {
            "parameter": "phi_n",
            "value": [parameters.phi_n_mean - parameters.phi_n_amplitude, parameters.phi_n_mean + parameters.phi_n_amplitude],
            "unit": "dimensionless",
            "role": "spatially varying nano-obstacle excluded fraction",
            "source_or_assumption": "frozen sinusoidal field; source uses local excluded fraction",
        },
        {
            "parameter": "k_visc",
            "value": parameters.k_visc,
            "unit": "dimensionless",
            "role": "relative-volume viscosity coefficient",
            "source_or_assumption": "ideal-solvent/Stokes–Einstein closure",
        },
        {
            "parameter": "V_r",
            "value": parameters.relative_volume,
            "unit": "dimensionless",
            "role": "cell volume after osmotic or crowding change",
            "source_or_assumption": "frozen nominal value; source equation uses relative cell volume",
        },
        {
            "parameter": "tortuosity_coefficient",
            "value": parameters.tortuosity_coefficient,
            "unit": "dimensionless",
            "role": "path-length detour closure coefficient",
            "source_or_assumption": "coarse-grained closure assumption",
        },
        {
            "parameter": "interaction_strength",
            "value": parameters.interaction_strength,
            "unit": "dimensionless",
            "role": "crowder-dependent permeability and encounter rate",
            "source_or_assumption": "frozen interaction closure",
        },
        {
            "parameter": "k_on0",
            "value": parameters.k_on0,
            "unit": "s^-1",
            "role": "reference reversible binding onset",
            "source_or_assumption": "frozen binding assumption; no internal data",
        },
        {
            "parameter": "k_off0",
            "value": parameters.k_off0,
            "unit": "s^-1",
            "role": "reference reversible binding offset",
            "source_or_assumption": "frozen binding assumption; no internal data",
        },
    ]


def sensitivity(parameters: Parameters) -> list[dict]:
    names = ("phi_n_mean", "interaction_strength", "k_on0")
    rows: list[dict] = []
    baseline = compare_models(parameters)
    for name in names:
        base_value = float(getattr(parameters, name))
        for multiplier in (0.5, 1.5):
            perturbed = replace(parameters, **{name: base_value * multiplier})
            result = compare_models(perturbed)
            rows.append(
                {
                    "parameter": name,
                    "multiplier": multiplier,
                    "base_value": base_value,
                    "value": float(getattr(perturbed, name)),
                    "passage_time_ratio_vs_adjusted": result["passage_time_ratio"],
                    "mechanistic_passage_time_s": result["mechanistic"]["passage_time_s"],
                    "adjusted_passage_time_s": result["adjusted"]["passage_time_s"],
                    "bound_occupancy_range": result["mechanistic"]["bound_occupancy_range"],
                    "bound_occupancy_cv": result["mechanistic"]["bound_occupancy_cv"],
                }
            )
    baseline_values = {
        "passage_time_ratio_vs_adjusted": baseline["passage_time_ratio"],
        "mechanistic_passage_time_s": baseline["mechanistic"]["passage_time_s"],
        "adjusted_passage_time_s": baseline["adjusted"]["passage_time_s"],
        "bound_occupancy_range": baseline["mechanistic"]["bound_occupancy_range"],
        "bound_occupancy_cv": baseline["mechanistic"]["bound_occupancy_cv"],
    }
    for row in rows:
        row["absolute_change_from_nominal"] = abs(
            row["passage_time_ratio_vs_adjusted"] - baseline_values["passage_time_ratio_vs_adjusted"]
        )
    return rows


def run(parameters: Parameters = Parameters()) -> dict:
    parameters.validate()
    comparison = compare_models(parameters)
    sanity = analytical_sanity(parameters)
    checks = {
        "mechanistic_separation": bool(abs(comparison["passage_time_ratio"]) >= 0.10),
        "spatial_signature": bool(
            comparison["mechanistic"]["bound_occupancy_range"] >= 0.10
            or comparison["mechanistic"]["bound_occupancy_cv"] >= 0.10
        ),
        "analytical_sanity": sanity["passed"],
    }
    return {
        "id": "BT-HX-Q084",
        "model_name": "hierarchical crowding with hydrodynamic hindrance and reversible binding",
        "provenance": {
            "source_reference": {
                "citation": "Destrian et al., PNAS 123(4), e2519599123 (2026)",
                "doi": "10.1073/pnas.2519599123",
                "pmcid": "PMC12846847",
                "figures": ["Fig. 3C", "Fig. 3D", "Fig. 4D-G"],
                "verified_values": {
                    "gfp_iso_osmotic_diffusivity_um2_per_s": [16.0, 20.0],
                    "low_porosity_relative_reduction_percent": 20.0,
                    "micro_obstacle_reduction_upper_percent": 10.0,
                    "nano_obstacle_reduction_upper_percent": 70.0,
                },
                "status": "VERIFIED",
            },
            "derivation": [
                "finite-volume reflecting-to-absorbing random walk",
                "accessible-volume and path-length tortuosity closure",
                "Kozeny-Carman-like permeability with particle hydrodynamic hindrance",
                "reversible two-state binding",
            ],
            "hypothesis": "local crowding and binding-state residence times differ from a mean-matched scalar diffusivity",
        },
        "parameters": asdict(parameters),
        "equations": EQUATIONS,
        "parameter_table": parameter_table(parameters),
        "unit_check": unit_check(parameters),
        "analytical_sanity": sanity,
        "comparison": comparison,
        "criteria": {
            "frozen_mechanistic_separation_threshold": 0.10,
            "frozen_spatial_range_or_cv_threshold": 0.10,
            "frozen_empty_space_relative_error_threshold": 0.02,
            "mechanistic_separation": checks["mechanistic_separation"],
            "spatial_signature": checks["spatial_signature"],
            "analytical_sanity": checks["analytical_sanity"],
            "overall": bool(all(checks.values())),
        },
        "sensitivity": sensitivity(parameters),
        "limitations": [
            "One-dimensional bounded domain and a deterministic sinusoidal obstacle field are not a measured cell geometry.",
            "The Kozeny-Carman-like closure and binding rates are mechanistic assumptions, not fitted internal data.",
            "A finite-volume absorbing outlet estimates first exit, not a biological target-search time.",
            "The scalar comparator is mean-matched and therefore tests spatial/binding mechanisms, not a claim that scalar diffusion is generally sufficient.",
        ],
    }


def write_results(path: str | Path, parameters: Parameters = Parameters()) -> dict:
    result = run(parameters)
    destination = Path(path)
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results.json")
    args = parser.parse_args()
    result = write_results(args.output)
    comparison = result["comparison"]
    summary = {
        "id": result["id"],
        "passage_time_ratio": comparison["passage_time_ratio"],
        "mechanistic_passage_time_s": comparison["mechanistic"]["passage_time_s"],
        "adjusted_passage_time_s": comparison["adjusted"]["passage_time_s"],
        "criteria": result["criteria"],
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
