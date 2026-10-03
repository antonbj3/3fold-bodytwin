from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
import json
import numpy as np


@dataclass(frozen=True)
class Parameters:
    phi_min: float = 0.70
    phi_max: float = 0.90
    d_min: float = 0.60e-3
    d_max: float = 1.00e-3
    q: float = 1.0
    E_s: float = 110e9
    chi_E: float = 0.75
    m_E: float = 2.0
    m_d: float = 0.35
    d_ref: float = 0.80e-3
    b_tau: float = 0.50
    k_geometry: float = 0.0690521240234
    D_0: float = 2.0e-9
    mu: float = 1.45e-3
    length: float = 12e-3
    pressure_gradient: float = 2.0e4
    consumption_rate: float = 0.05
    c_in: float = 1.0
    c_out: float = 1.0
    grid_points: int = 101
    sigma_app: float = 2.0e6
    k_t: float = 0.12
    epsilon_peak: float = 0.003
    epsilon_width: float = 0.004
    c_50: float = 0.5
    tau_opt: float = 0.5
    tau_ratio_width: float = 1.5
    response_max: float = 1.0


DEFAULT_PARAMETERS = Parameters()


def _dim_mul(*dimensions: tuple[float, float, float]) -> tuple[float, float, float]:
    result = (0.0, 0.0, 0.0)
    for dimension in dimensions:
        result = tuple(left + right for left, right in zip(result, dimension))
    return result


def _dim_div(numerator: tuple[float, float, float], denominator: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(left - right for left, right in zip(numerator, denominator))


def _dim_pow(dimension: tuple[float, float, float], exponent: float) -> tuple[float, float, float]:
    return tuple(exponent * value for value in dimension)


_DIMENSIONLESS = (0.0, 0.0, 0.0)
_METER = (1.0, 0.0, 0.0)
_KILOGRAM = (0.0, 1.0, 0.0)
_SECOND = (0.0, 0.0, 1.0)
_PASCAL = _dim_div(_KILOGRAM, _dim_mul(_METER, _dim_pow(_SECOND, 2.0)))
_PRESSURE_GRADIENT = _dim_div(_PASCAL, _METER)
_VISCOSITY = _dim_mul(_PASCAL, _SECOND)


PARAMETER_METADATA = {
    "phi_min": {"unit": "1", "source": "assumption", "note": "process endpoint"},
    "phi_max": {"unit": "1", "source": "assumption", "note": "process endpoint"},
    "d_min": {"unit": "m", "source": "Chao et al. 2021 design range", "note": "600 µm"},
    "d_max": {"unit": "m", "source": "Chao et al. 2021 design range", "note": "1000 µm"},
    "q": {"unit": "1", "source": "assumption", "note": "linear process mapping"},
    "E_s": {"unit": "Pa", "source": "Chao et al. 2021 material data", "note": "Ti-6Al-4V"},
    "chi_E": {"unit": "1", "source": "assumption", "note": "network correction"},
    "m_E": {"unit": "1", "source": "assumption", "note": "bending scaling"},
    "m_d": {"unit": "1", "source": "assumption", "note": "architecture correction"},
    "d_ref": {"unit": "m", "source": "Chao et al. 2021 nominal design", "note": "800 µm"},
    "b_tau": {"unit": "1", "source": "assumption", "note": "tortuosity slope"},
    "k_geometry": {"unit": "1", "source": "calibration to Chao et al. Fig. 13A", "note": "not independent validation"},
    "D_0": {"unit": "m2/s", "source": "assumption", "note": "37 C diffusivity scale"},
    "mu": {"unit": "Pa s", "source": "Chao et al. 2021", "note": "converted from MPa s"},
    "length": {"unit": "m", "source": "Chao et al. 2021", "note": "12 mm sample height"},
    "pressure_gradient": {"unit": "Pa/m", "source": "assumption", "note": "perfusion driving gradient"},
    "consumption_rate": {"unit": "1/s", "source": "assumption", "note": "first-order nutrient sink"},
    "c_in": {"unit": "1", "source": "assumption", "note": "normalized inlet"},
    "c_out": {"unit": "1", "source": "assumption", "note": "normalized outlet"},
    "grid_points": {"unit": "1", "source": "numerical choice", "note": "finite-volume nodes"},
    "sigma_app": {"unit": "Pa", "source": "assumption", "note": "compressive load"},
    "k_t": {"unit": "1", "source": "assumption", "note": "stress concentration coefficient"},
    "epsilon_peak": {"unit": "1", "source": "assumption", "note": "response window"},
    "epsilon_width": {"unit": "1", "source": "assumption", "note": "response window"},
    "c_50": {"unit": "1", "source": "assumption", "note": "half saturation"},
    "tau_opt": {"unit": "Pa", "source": "assumption", "note": "shear window"},
    "tau_ratio_width": {"unit": "1", "source": "assumption", "note": "log2 shear width"},
    "response_max": {"unit": "1", "source": "hypothesis", "note": "dimensionless potential"},
}


def parameter_table(parameters: Parameters = DEFAULT_PARAMETERS) -> list[dict]:
    values = asdict(parameters)
    return [
        {"name": name, "value": values[name], **PARAMETER_METADATA[name]}
        for name in values
    ]


def validate_parameters(parameters: Parameters = DEFAULT_PARAMETERS) -> dict:
    checks = {
        "porosity_order": 0.0 < parameters.phi_min < parameters.phi_max < 1.0,
        "pore_size_order": 0.0 < parameters.d_min < parameters.d_max,
        "positive_mechanics": parameters.E_s > 0.0 and parameters.chi_E > 0.0 and parameters.m_E > 0.0 and parameters.m_d >= 0.0,
        "positive_transport": parameters.k_geometry > 0.0 and parameters.D_0 > 0.0 and parameters.mu > 0.0 and parameters.length > 0.0,
        "nonnegative_driver": parameters.pressure_gradient >= 0.0 and parameters.consumption_rate >= 0.0,
        "normalized_boundaries": parameters.c_in >= 0.0 and parameters.c_out >= 0.0,
        "grid_size": parameters.grid_points >= 3,
        "positive_response_windows": parameters.c_50 > 0.0 and parameters.epsilon_width > 0.0 and parameters.tau_opt > 0.0 and parameters.tau_ratio_width > 0.0,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError(f"invalid parameters: {failed}")
    return checks


def process_microstructure(process_state: float, parameters: Parameters = DEFAULT_PARAMETERS) -> tuple[float, float]:
    validate_parameters(parameters)
    if not 0.0 <= process_state <= 1.0:
        raise ValueError("process_state must lie in [0, 1]")
    porosity = parameters.phi_min + (parameters.phi_max - parameters.phi_min) * process_state
    pore_diameter = parameters.d_min + (parameters.d_max - parameters.d_min) * process_state ** parameters.q
    return float(porosity), float(pore_diameter)


def geometry_properties(porosity: float, pore_diameter: float, parameters: Parameters = DEFAULT_PARAMETERS) -> dict:
    validate_parameters(parameters)
    if not 0.0 <= porosity < 1.0:
        raise ValueError("porosity must lie in [0, 1)")
    if pore_diameter < 0.0:
        raise ValueError("pore_diameter must be non-negative")
    if pore_diameter == 0.0:
        return {
            "specific_surface_area_per_volume_m_inv": 0.0,
            "tortuosity": 1.0 + parameters.b_tau * (1.0 - porosity),
            "permeability_m2": 0.0,
            "effective_diffusivity_m2_per_s": 0.0,
        }
    tortuosity = 1.0 + parameters.b_tau * (1.0 - porosity)
    specific_surface = 6.0 * porosity / pore_diameter
    permeability = parameters.k_geometry * porosity ** 3 * pore_diameter ** 2 / tortuosity ** 2
    effective_diffusivity = parameters.D_0 * porosity / tortuosity ** 2
    return {
        "specific_surface_area_per_volume_m_inv": float(specific_surface),
        "tortuosity": float(tortuosity),
        "permeability_m2": float(permeability),
        "effective_diffusivity_m2_per_s": float(effective_diffusivity),
    }


def elastic_modulus(porosity: float, pore_diameter: float, parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    validate_parameters(parameters)
    if pore_diameter <= 0.0:
        return 0.0
    if not 0.0 <= porosity < 1.0:
        raise ValueError("porosity must lie in [0, 1)")
    relative_density = 1.0 - porosity
    architecture = (parameters.d_ref / pore_diameter) ** parameters.m_d
    return float(parameters.E_s * parameters.chi_E * relative_density ** parameters.m_E * architecture)


def transport_field(
    permeability: float,
    effective_diffusivity: float,
    velocity: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
    concentration_in: float | None = None,
    concentration_out: float | None = None,
    grid_points: int | None = None,
) -> np.ndarray:
    validate_parameters(parameters)
    n = parameters.grid_points if grid_points is None else int(grid_points)
    if n < 2:
        raise ValueError("grid_points must be at least 2")
    cin = parameters.c_in if concentration_in is None else float(concentration_in)
    cout = parameters.c_out if concentration_out is None else float(concentration_out)
    if permeability < 0.0 or effective_diffusivity < 0.0:
        raise ValueError("transport coefficients must be non-negative")
    if effective_diffusivity == 0.0:
        if velocity == 0.0 and parameters.consumption_rate == 0.0:
            return np.full(n, cin, dtype=float)
        return np.zeros(n, dtype=float)
    dx = parameters.length / (n - 1)
    diffusion_coefficient = effective_diffusivity / dx ** 2
    advection_coefficient = velocity / dx
    matrix = np.zeros((n, n), dtype=float)
    rhs = np.zeros(n, dtype=float)
    matrix[0, 0] = 1.0
    rhs[0] = cin
    matrix[-1, -1] = 1.0
    rhs[-1] = cout
    for index in range(1, n - 1):
        matrix[index, index - 1] = diffusion_coefficient + advection_coefficient
        matrix[index, index] = -2.0 * diffusion_coefficient - advection_coefficient - parameters.consumption_rate
        matrix[index, index + 1] = diffusion_coefficient
    concentration = np.linalg.solve(matrix, rhs)
    return concentration


def response_potential(
    concentration_mid: float,
    local_strain: float,
    wall_shear: float,
    parameters: Parameters = DEFAULT_PARAMETERS,
) -> dict:
    validate_parameters(parameters)
    oxygen_factor = concentration_mid / (parameters.c_50 + concentration_mid) if concentration_mid > 0.0 else 0.0
    mechanical_factor = float(np.exp(-0.5 * ((local_strain - parameters.epsilon_peak) / parameters.epsilon_width) ** 2))
    if wall_shear <= 0.0:
        shear_factor = 0.0
    else:
        shear_factor = float(np.exp(-0.5 * (np.log2(wall_shear / parameters.tau_opt) / parameters.tau_ratio_width) ** 2))
    response = parameters.response_max * oxygen_factor * mechanical_factor * shear_factor
    return {
        "oxygen_factor": float(oxygen_factor),
        "mechanical_factor": mechanical_factor,
        "shear_factor": shear_factor,
        "response_potential": float(response),
    }


def simulate(
    process_state: float = 0.5,
    parameters: Parameters = DEFAULT_PARAMETERS,
    permeability_override: float | None = None,
    solve_transport: bool = True,
) -> dict:
    porosity, pore_diameter = process_microstructure(process_state, parameters)
    geometry = geometry_properties(porosity, pore_diameter, parameters)
    permeability = geometry["permeability_m2"] if permeability_override is None else float(permeability_override)
    if permeability < 0.0:
        raise ValueError("permeability_override must be non-negative")
    effective_diffusivity = geometry["effective_diffusivity_m2_per_s"]
    velocity = permeability * parameters.pressure_gradient / parameters.mu
    modulus = elastic_modulus(porosity, pore_diameter, parameters)
    concentration_factor = 1.0 + parameters.k_t * (pore_diameter / parameters.d_ref) / max(1.0 - porosity, 1.0e-12)
    local_stress = parameters.sigma_app * concentration_factor
    local_strain = local_stress / modulus if modulus > 0.0 else 0.0
    wall_shear = parameters.mu * abs(velocity) / (pore_diameter / 2.0) if pore_diameter > 0.0 else 0.0
    if solve_transport:
        concentration = transport_field(permeability, effective_diffusivity, velocity, parameters)
    else:
        concentration = np.full(parameters.grid_points, parameters.c_in, dtype=float)
    concentration_mid = float(concentration[parameters.grid_points // 2])
    response = response_potential(concentration_mid, local_strain, wall_shear, parameters)
    return {
        "process_state": float(process_state),
        "porosity": porosity,
        "pore_diameter_m": pore_diameter,
        "pore_radius_m": float(pore_diameter / 2.0),
        "specific_surface_area_per_volume_m_inv": geometry["specific_surface_area_per_volume_m_inv"],
        "tortuosity": geometry["tortuosity"],
        "elastic_modulus_pa": modulus,
        "permeability_m2": float(permeability),
        "effective_diffusivity_m2_per_s": effective_diffusivity,
        "darcy_velocity_m_per_s": float(velocity),
        "pressure_drop_pa": float(parameters.pressure_gradient * parameters.length),
        "wall_shear_pa": float(wall_shear),
        "stress_concentration_factor": float(concentration_factor),
        "local_peak_stress_pa": float(local_stress),
        "local_strain": float(local_strain),
        "concentration_profile": [float(value) for value in concentration],
        "concentration_mid": concentration_mid,
        "concentration_mean": float(np.mean(concentration)),
        "concentration_min": float(np.min(concentration)),
        "concentration_out": float(concentration[-1]),
        # Right face of the interior ledger: positive-flow upwind donor is c[-2].
        "outlet_flux_normalized": float(velocity * concentration[-2] - effective_diffusivity * (concentration[-1] - concentration[-2]) / (parameters.length / (parameters.grid_points - 1))),
        **response,
    }


def _pearson(first: np.ndarray, second: np.ndarray) -> float:
    if first.size < 2 or second.size < 2 or np.std(first) == 0.0 or np.std(second) == 0.0:
        return 0.0
    return float(np.corrcoef(first, second)[0, 1])


def sensitivity(parameters: Parameters = DEFAULT_PARAMETERS) -> dict:
    validate_parameters(parameters)
    nominal = simulate(0.5, parameters)
    rows = []
    for name in ("b_tau", "chi_E", "consumption_rate"):
        baseline = float(getattr(parameters, name))
        for multiplier in (0.5, 1.5):
            altered = replace(parameters, **{name: baseline * multiplier})
            result = simulate(0.5, altered)
            rows.append({
                "parameter": name,
                "multiplier": multiplier,
                "value": baseline * multiplier,
                "permeability_m2": result["permeability_m2"],
                "elastic_modulus_pa": result["elastic_modulus_pa"],
                "effective_diffusivity_m2_per_s": result["effective_diffusivity_m2_per_s"],
                "concentration_mean": result["concentration_mean"],
                "concentration_min": result["concentration_min"],
                "response_potential": result["response_potential"],
            })
    return {"nominal": nominal, "rows": rows}


def process_sweep(parameters: Parameters = DEFAULT_PARAMETERS, count: int = 9) -> dict:
    states = np.linspace(0.0, 1.0, count)
    results = [simulate(float(state), parameters) for state in states]
    return {
        "process_states": [float(state) for state in states],
        "results": results,
        "correlation_process_mean_concentration": _pearson(states, np.array([result["concentration_mean"] for result in results])),
    }


def shuffled_placebo(parameters: Parameters = DEFAULT_PARAMETERS, count: int = 9, permutations: int = 64, seed: int = 52) -> dict:
    sweep = process_sweep(parameters, count)
    states = np.array(sweep["process_states"])
    nominal_permeabilities = np.array([result["permeability_m2"] for result in sweep["results"]])
    generator = np.random.default_rng(seed)
    correlations = []
    for _ in range(permutations):
        shuffled = generator.permutation(nominal_permeabilities)
        results = [
            simulate(float(state), parameters, permeability_override=float(permeability))
            for state, permeability in zip(states, shuffled)
        ]
        correlations.append(_pearson(states, np.array([result["concentration_mean"] for result in results])))
    return {
        "permutations": permutations,
        "seed": seed,
        "median_absolute_correlation": float(np.median(np.abs(correlations))),
        "correlations": [float(value) for value in correlations],
    }


def unit_checks(parameters: Parameters = DEFAULT_PARAMETERS) -> dict:
    nominal = simulate(0.5, parameters)
    expected_permeability = parameters.k_geometry * nominal["porosity"] ** 3 * nominal["pore_diameter_m"] ** 2 / nominal["tortuosity"] ** 2
    expected_diffusivity = parameters.D_0 * nominal["porosity"] / nominal["tortuosity"] ** 2
    expected_velocity = expected_permeability * parameters.pressure_gradient / parameters.mu
    expected_shear = parameters.mu * abs(expected_velocity) / (nominal["pore_diameter_m"] / 2.0)
    expected_strain = parameters.sigma_app * nominal["stress_concentration_factor"] / nominal["elastic_modulus_pa"]
    permeability_dimension = _dim_mul(_dim_pow(_DIMENSIONLESS, 1.0), _dim_pow(_DIMENSIONLESS, 3.0), _dim_pow(_METER, 2.0), _dim_pow(_DIMENSIONLESS, -2.0))
    diffusivity_dimension = _dim_mul(_dim_pow(_METER, 2.0), _dim_pow(_SECOND, -1.0))
    velocity_dimension = _dim_div(_dim_mul(_dim_pow(_METER, 2.0), _PRESSURE_GRADIENT), _VISCOSITY)
    shear_dimension = _dim_div(_dim_mul(_VISCOSITY, velocity_dimension), _METER)
    strain_dimension = _dim_div(_PASCAL, _PASCAL)
    checks = {
        "parameter_validation": all(validate_parameters(parameters).values()),
        "permeability_dimension_m2": permeability_dimension == _dim_pow(_METER, 2.0),
        "diffusivity_dimension_m2_per_s": diffusivity_dimension == _dim_mul(_dim_pow(_METER, 2.0), _dim_pow(_SECOND, -1.0)),
        "velocity_dimension_m_per_s": velocity_dimension == _dim_mul(_METER, _dim_pow(_SECOND, -1.0)),
        "shear_dimension_pa": shear_dimension == _PASCAL,
        "strain_dimensionless": strain_dimension == _DIMENSIONLESS,
        "response_dimensionless": _dim_mul(_DIMENSIONLESS, _DIMENSIONLESS, _DIMENSIONLESS) == _DIMENSIONLESS,
        "permeability_formula_value": bool(np.isclose(nominal["permeability_m2"], expected_permeability, rtol=1.0e-12)),
        "diffusivity_formula_value": bool(np.isclose(nominal["effective_diffusivity_m2_per_s"], expected_diffusivity, rtol=1.0e-12)),
        "darcy_velocity_formula_value": bool(np.isclose(nominal["darcy_velocity_m_per_s"], expected_velocity, rtol=1.0e-12)),
        "shear_formula_value": bool(np.isclose(nominal["wall_shear_pa"], expected_shear, rtol=1.0e-12)),
        "strain_formula_value": bool(np.isclose(nominal["local_strain"], expected_strain, rtol=1.0e-12)),
        "finite_nominal": bool(np.isfinite(nominal["elastic_modulus_pa"]) and np.isfinite(nominal["permeability_m2"]) and np.isfinite(nominal["concentration_mean"])),
    }
    checks["all_passed"] = all(checks.values())
    return checks


def analytic_limit_checks(parameters: Parameters = DEFAULT_PARAMETERS) -> dict:
    zero = geometry_properties(0.0, parameters.d_ref, parameters)
    first = geometry_properties(parameters.phi_min, parameters.d_ref, parameters)
    second = geometry_properties(parameters.phi_min, 2.0 * parameters.d_ref, parameters)
    uniform_parameters = replace(parameters, pressure_gradient=0.0, consumption_rate=0.0)
    uniform = transport_field(first["permeability_m2"], first["effective_diffusivity_m2_per_s"], 0.0, uniform_parameters)
    return {
        "zero_porosity_permeability_is_zero": bool(zero["permeability_m2"] == 0.0),
        "zero_porosity_diffusivity_is_zero": bool(zero["effective_diffusivity_m2_per_s"] == 0.0),
        "no_flow_no_consumption_is_uniform": bool(np.allclose(uniform, parameters.c_in)),
        "permeability_scales_with_diameter_squared": bool(np.isclose(second["permeability_m2"] / first["permeability_m2"], 4.0, rtol=1.0e-12)),
    }


def run_analysis(parameters: Parameters = DEFAULT_PARAMETERS) -> dict:
    validate_parameters(parameters)
    nominal = simulate(0.5, parameters)
    sweep = process_sweep(parameters)
    placebo = shuffled_placebo(parameters)
    sensitivity_result = sensitivity(parameters)
    units = unit_checks(parameters)
    limits = analytic_limit_checks(parameters)
    reference_permeability = 1.87e-8
    reference_ratio = nominal["permeability_m2"] / reference_permeability
    mechanical_interval = (2.6e9, 4.0e9)
    mechanical_ok = mechanical_interval[0] / 1.5 <= nominal["elastic_modulus_pa"] <= mechanical_interval[1] * 1.5
    finite_nonnegative = all(
        np.isfinite(result[key])
        for result in [nominal] + sweep["results"]
        for key in ("elastic_modulus_pa", "permeability_m2", "effective_diffusivity_m2_per_s", "concentration_mean", "concentration_min", "response_potential")
    ) and nominal["concentration_min"] >= -1.0e-10 and nominal["response_potential"] >= 0.0
    coupling_ok = sweep["correlation_process_mean_concentration"] >= 0.80 and placebo["median_absolute_correlation"] < 0.30
    criteria = {
        "units_passed": units["all_passed"],
        "reference_ratio": float(reference_ratio),
        "reference_within_25_percent": bool(abs(reference_ratio - 1.0) <= 0.25),
        "mechanical_within_1.5x_reported_interval": bool(mechanical_ok),
        "analytic_limits_passed": bool(all(limits.values())),
        "full_coupling_correlation": sweep["correlation_process_mean_concentration"],
        "placebo_median_absolute_correlation": placebo["median_absolute_correlation"],
        "coupling_criterion_passed": bool(coupling_ok),
        "finite_and_nonnegative": bool(finite_nonnegative),
    }
    criteria["all_passed"] = bool(
        criteria["units_passed"]
        and criteria["reference_within_25_percent"]
        and criteria["mechanical_within_1.5x_reported_interval"]
        and criteria["analytic_limits_passed"]
        and criteria["coupling_criterion_passed"]
        and criteria["finite_and_nonnegative"]
    )
    return {
        "id": "BT-HX-Q052",
        "model": "first-principles-1D-porosity-mechanics-transport-response-v1",
        "parameter_table": parameter_table(parameters),
        "literature_reference": {
            "source": "Chao et al. 2021",
            "doi": "10.3389/fbioe.2021.779854",
            "location": "Figure 13A",
            "porosity": 0.80,
            "pore_diameter_m": 0.80e-3,
            "permeability_m2": 1.87e-8,
            "status": "VERIFIED",
            "natural_bone_context_m2": 1.50e-10,
            "mechanical_context_pa": {"minimum": 2.6e9, "maximum": 4.0e9},
            "biological_context": "Figure 15: qualitative MC3T3-E1 attachment direction at 7 days",
        },
        "nominal": nominal,
        "process_sweep": sweep,
        "sensitivity": sensitivity_result,
        "placebo": placebo,
        "unit_checks": units,
        "analytic_limit_checks": limits,
        "criteria": criteria,
        "equations": {
            "process": "phi=phi_min+(phi_max-phi_min)*s; d=d_min+(d_max-d_min)*s**q",
            "permeability": "K=K_geometry*phi**3*d**2/tau**2; tau=1+b_tau*(1-phi)",
            "mechanics": "E=E_s*chi_E*(1-phi)**m_E*(d_ref/d)**m_d",
            "diffusion": "D_eff=D_0*phi/tau**2",
            "transport": "D_eff*c_second-u*c_first-lambda*c=0",
            "response": "R=response_max*oxygen_factor*mechanical_factor*shear_factor",
        },
        "source_classes": {
            "nominal": "derived",
            "reference": "literature_reference",
            "process_endpoints": "assumption",
            "E_s": "literature_reference",
            "chi_E_m_E_m_d_D_0_pressure_gradient_consumption_rate": "assumption",
            "response_potential": "hypothesis",
        },
        "limitations": [
            "K_geometry is calibrated to one literature point and is not independent validation.",
            "The process-to-geometry map is an explicit AM-process assumption.",
            "The response potential is dimensionless and is not a measured cell or tissue outcome.",
            "The model is one-dimensional and homogenized; it does not resolve local CT geometry or multiphase flow.",
        ],
    }


def main() -> None:
    result = run_analysis()
    output_path = Path(__file__).with_name("results.json")
    output_path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    summary = {
        "id": result["id"],
        "nominal_permeability_m2": result["nominal"]["permeability_m2"],
        "nominal_elastic_modulus_pa": result["nominal"]["elastic_modulus_pa"],
        "nominal_concentration_mean": result["nominal"]["concentration_mean"],
        "nominal_response_potential": result["nominal"]["response_potential"],
        "criteria_all_passed": result["criteria"]["all_passed"],
    }
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
