from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from scipy.integrate import quad
from scipy.special import ellipe


@dataclass(frozen=True)
class Parameters:
    a_m: float = 0.023
    b_m: float = 0.0065
    thickness_m: float = 0.002
    prestress_pa: float = 20000.0
    modulus_pa: float = 20000.0
    nonlinearity_beta: float = 2.0
    reference_length_m: float = 0.020
    tip_radius_m: float = 0.0001
    pressure_delta_pa: float = 1000.0
    viscosity_pa_s: float = 0.001
    tissue_permeability_m2: float = 1.0e-15
    intact_path_length_m: float = 0.010
    wound_path_length_m: float = 0.002
    intact_tewl_g_m2_h: float = 1.0
    barrier_damage_multiplier: float = 5.9
    exudate_velocity_m_s: float = 1.0e-8

    def validate(self) -> None:
        positive = (
            "a_m",
            "thickness_m",
            "modulus_pa",
            "reference_length_m",
            "tip_radius_m",
            "viscosity_pa_s",
            "tissue_permeability_m2",
            "intact_path_length_m",
            "wound_path_length_m",
            "intact_tewl_g_m2_h",
            "barrier_damage_multiplier",
        )
        for name in positive:
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("b_m", "prestress_pa", "nonlinearity_beta", "exudate_velocity_m_s"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.b_m > self.a_m:
            raise ValueError("b_m must not exceed a_m")


@dataclass(frozen=True)
class ForceResult:
    prestress_force_n: float
    incremental_force_n: float
    total_force_n: float
    opening_width_m: float
    opening_strain_max: float
    opening_strain_mean: float
    nominal_stress_max_pa: float
    tip_stress_pa: float
    tip_stress_factor: float


@dataclass(frozen=True)
class TransportResult:
    intact_evap_flux_m_s: float
    cut_evap_flux_m_s: float
    intact_darcy_flux_m_s: float
    cut_darcy_flux_m_s: float
    intact_evap_flow_m3_s: float
    cut_evap_flow_m3_s: float
    intact_darcy_flow_m3_s: float
    cut_darcy_flow_m3_s: float
    exudate_flow_m3_s: float
    intact_total_flow_m3_s: float
    cut_total_flow_m3_s: float
    total_flow_ratio: float


DEFAULT_PARAMETERS = Parameters()


def parameter_table(parameters: Parameters = DEFAULT_PARAMETERS) -> list[dict[str, Any]]:
    parameters.validate()
    metadata = {
        "a_m": ("m", "standard geometry; 46 mm full length"),
        "b_m": ("m", "standard geometry; 13 mm full width"),
        "thickness_m": ("m", "assumption; load-bearing tissue thickness"),
        "prestress_pa": ("Pa", "assumption; natural membrane stress"),
        "modulus_pa": ("Pa", "assumption; small-strain effective modulus"),
        "nonlinearity_beta": ("1", "assumption; exponential stress exponent"),
        "reference_length_m": ("m", "assumption; opening strain gauge length"),
        "tip_radius_m": ("m", "assumption; rounded incision tip radius"),
        "pressure_delta_pa": ("Pa", "assumption; tissue-to-surface pressure drop"),
        "viscosity_pa_s": ("Pa s", "assumption; aqueous wound fluid viscosity"),
        "tissue_permeability_m2": ("m^2", "assumption; intrinsic tissue permeability"),
        "intact_path_length_m": ("m", "assumption; intact transport path"),
        "wound_path_length_m": ("m", "assumption; open wound transport path"),
        "intact_tewl_g_m2_h": ("g m^-2 h^-1", "assumption; normalized intact barrier flux"),
        "barrier_damage_multiplier": (
            "1",
            "Barthe et al. 2024, DOI 10.3389/fmed.2024.1481645, Fig. 2; 5.9-fold TEWL analog",
        ),
        "exudate_velocity_m_s": ("m s^-1", "assumption; wound-fluid source velocity"),
    }
    values = asdict(parameters)
    return [
        {
            "name": name,
            "value": float(values[name]),
            "unit": metadata[name][0],
            "source_or_assumption": metadata[name][1],
        }
        for name in values
    ]


def ellipse_geometry(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, float]:
    parameters.validate()
    a = parameters.a_m
    b = parameters.b_m
    area = math.pi * a * b
    if b == 0.0:
        perimeter = 4.0 * a
    else:
        eccentricity_squared = 1.0 - (b / a) ** 2
        perimeter = 4.0 * a * float(ellipe(eccentricity_squared))
    return {
        "a_m": float(a),
        "b_m": float(b),
        "full_length_m": float(2.0 * a),
        "full_width_m": float(2.0 * b),
        "exposed_area_m2": float(area),
        "intact_exposed_area_m2": 0.0,
        "cut_exposed_area_m2": float(area),
        "exposed_perimeter_m": float(perimeter),
        "opening_width_m": float(2.0 * b),
        "opening_strain_max": float(2.0 * b / parameters.reference_length_m),
        "opening_strain_mean": float(math.pi * b / (2.0 * parameters.reference_length_m)),
    }


def local_opening_width(y_m: float, parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    parameters.validate()
    if abs(y_m) > parameters.a_m:
        raise ValueError("y_m must lie on the cut interval")
    shape = math.sqrt(max(0.0, 1.0 - (y_m / parameters.a_m) ** 2))
    return float(2.0 * parameters.b_m * shape)


def opening_stress_pa(strain: float, parameters: Parameters = DEFAULT_PARAMETERS) -> float:
    parameters.validate()
    if strain < 0.0:
        raise ValueError("strain must be nonnegative")
    if parameters.nonlinearity_beta < 1.0e-12:
        return float(parameters.prestress_pa + parameters.modulus_pa * strain)
    value = parameters.modulus_pa * math.expm1(
        parameters.nonlinearity_beta * strain
    ) / parameters.nonlinearity_beta
    return float(parameters.prestress_pa + value)


def closure_force(parameters: Parameters = DEFAULT_PARAMETERS) -> ForceResult:
    parameters.validate()
    a = parameters.a_m
    b = parameters.b_m
    h = parameters.thickness_m
    beta = parameters.nonlinearity_beta
    reference_length = parameters.reference_length_m
    prestress_force = 2.0 * a * h * parameters.prestress_pa
    if b == 0.0:
        incremental_force = 0.0
    elif beta < 1.0e-12:
        incremental_force = parameters.modulus_pa * h * math.pi * a * b / reference_length
    else:
        k = beta * 2.0 * b / reference_length

        def integrand(u: float) -> float:
            shape = math.sqrt(max(0.0, 1.0 - u * u))
            return math.expm1(k * shape)

        integral, _ = quad(integrand, -1.0, 1.0, epsabs=1.0e-12, epsrel=1.0e-10)
        incremental_force = parameters.modulus_pa * h * a * integral / beta
    total_force = prestress_force + incremental_force
    opening_strain_max = 2.0 * b / reference_length
    nominal_stress_max = opening_stress_pa(opening_strain_max, parameters)
    tip_factor = 1.0 + 2.0 * math.sqrt(a / parameters.tip_radius_m)
    return ForceResult(
        prestress_force_n=float(prestress_force),
        incremental_force_n=float(incremental_force),
        total_force_n=float(total_force),
        opening_width_m=float(2.0 * b),
        opening_strain_max=float(opening_strain_max),
        opening_strain_mean=float(math.pi * b / (2.0 * reference_length)),
        nominal_stress_max_pa=float(nominal_stress_max),
        tip_stress_pa=float(tip_factor * nominal_stress_max),
        tip_stress_factor=float(tip_factor),
    )


def tewl_to_flux_m_s(tewl_g_m2_h: float) -> float:
    if not math.isfinite(tewl_g_m2_h) or tewl_g_m2_h < 0.0:
        raise ValueError("TEWL must be finite and nonnegative")
    return float(tewl_g_m2_h * 1.0e-3 / (3600.0 * 1000.0))


def darcy_flux_m_s(permeability_m2: float, pressure_delta_pa: float, viscosity_pa_s: float, path_length_m: float) -> float:
    if permeability_m2 < 0.0 or viscosity_pa_s <= 0.0 or path_length_m <= 0.0:
        raise ValueError("invalid Darcy parameters")
    return float(permeability_m2 * pressure_delta_pa / (viscosity_pa_s * path_length_m))


def transport_flow(parameters: Parameters = DEFAULT_PARAMETERS) -> TransportResult:
    parameters.validate()
    area = ellipse_geometry(parameters)["exposed_area_m2"]
    intact_evap_flux = tewl_to_flux_m_s(parameters.intact_tewl_g_m2_h)
    cut_evap_flux = intact_evap_flux * parameters.barrier_damage_multiplier
    intact_darcy_flux = darcy_flux_m_s(
        parameters.tissue_permeability_m2,
        parameters.pressure_delta_pa,
        parameters.viscosity_pa_s,
        parameters.intact_path_length_m,
    )
    cut_darcy_flux = darcy_flux_m_s(
        parameters.tissue_permeability_m2,
        parameters.pressure_delta_pa,
        parameters.viscosity_pa_s,
        parameters.wound_path_length_m,
    )
    intact_evap_flow = intact_evap_flux * area
    cut_evap_flow = cut_evap_flux * area
    intact_darcy_flow = intact_darcy_flux * area
    cut_darcy_flow = cut_darcy_flux * area
    exudate_flow = parameters.exudate_velocity_m_s * area
    intact_total_flow = intact_evap_flow + intact_darcy_flow
    cut_total_flow = cut_evap_flow + cut_darcy_flow + exudate_flow
    if intact_total_flow == 0.0:
        flow_ratio = 0.0
    else:
        flow_ratio = cut_total_flow / intact_total_flow
    return TransportResult(
        intact_evap_flux_m_s=float(intact_evap_flux),
        cut_evap_flux_m_s=float(cut_evap_flux),
        intact_darcy_flux_m_s=float(intact_darcy_flux),
        cut_darcy_flux_m_s=float(cut_darcy_flux),
        intact_evap_flow_m3_s=float(intact_evap_flow),
        cut_evap_flow_m3_s=float(cut_evap_flow),
        intact_darcy_flow_m3_s=float(intact_darcy_flow),
        cut_darcy_flow_m3_s=float(cut_darcy_flow),
        exudate_flow_m3_s=float(exudate_flow),
        intact_total_flow_m3_s=float(intact_total_flow),
        cut_total_flow_m3_s=float(cut_total_flow),
        total_flow_ratio=float(flow_ratio),
    )


def _asdict_force(force: ForceResult) -> dict[str, float]:
    return {key: float(value) for key, value in asdict(force).items()}


def _asdict_transport(transport: TransportResult) -> dict[str, float]:
    return {key: float(value) for key, value in asdict(transport).items()}


def simulate(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    parameters.validate()
    geometry = ellipse_geometry(parameters)
    force = closure_force(parameters)
    transport = transport_flow(parameters)
    return {
        "geometry": geometry,
        "before_after": {
            "before_exposed_area_m2": 0.0,
            "after_exposed_area_m2": geometry["exposed_area_m2"],
            "before_opening_width_m": 0.0,
            "after_opening_width_m": geometry["opening_width_m"],
        },
        "mechanics": _asdict_force(force),
        "transport": _asdict_transport(transport),
        "units": {
            "geometry": {
                "a_m": "m",
                "b_m": "m",
                "full_length_m": "m",
                "full_width_m": "m",
                "exposed_area_m2": "m^2",
                "intact_exposed_area_m2": "m^2",
                "cut_exposed_area_m2": "m^2",
                "exposed_perimeter_m": "m",
                "opening_width_m": "m",
                "opening_strain_max": "1",
                "opening_strain_mean": "1",
            },
            "before_after": {
                "before_exposed_area_m2": "m^2",
                "after_exposed_area_m2": "m^2",
                "before_opening_width_m": "m",
                "after_opening_width_m": "m",
            },
            "mechanics": {
                "prestress_force_n": "N",
                "incremental_force_n": "N",
                "total_force_n": "N",
                "opening_width_m": "m",
                "opening_strain_max": "1",
                "opening_strain_mean": "1",
                "nominal_stress_max_pa": "Pa",
                "tip_stress_pa": "Pa",
                "tip_stress_factor": "1",
            },
            "transport": {
                "intact_evap_flux_m_s": "m s^-1",
                "cut_evap_flux_m_s": "m s^-1",
                "intact_darcy_flux_m_s": "m s^-1",
                "cut_darcy_flux_m_s": "m s^-1",
                "intact_evap_flow_m3_s": "m^3 s^-1",
                "cut_evap_flow_m3_s": "m^3 s^-1",
                "intact_darcy_flow_m3_s": "m^3 s^-1",
                "cut_darcy_flow_m3_s": "m^3 s^-1",
                "exudate_flow_m3_s": "m^3 s^-1",
                "intact_total_flow_m3_s": "m^3 s^-1",
                "cut_total_flow_m3_s": "m^3 s^-1",
                "total_flow_ratio": "1",
            },
        },
    }


_DIMENSIONLESS = (0, 0, 0)
_M = (0, 1, 0)
_S = (0, 0, 1)
_PA = (1, -1, -2)
_M2 = (0, 2, 0)
_N = (1, 1, -2)


def _dim_mul(left: tuple[int, int, int], right: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(a + b for a, b in zip(left, right))


def _dim_div(left: tuple[int, int, int], right: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(a - b for a, b in zip(left, right))


def _dim_pow(dimension: tuple[int, int, int], exponent: int) -> tuple[int, int, int]:
    return tuple(value * exponent for value in dimension)


def unit_check() -> dict[str, Any]:
    checks = {
        "stress": _dim_pow(_PA, 1) == _PA,
        "force_from_stress_thickness_length": _dim_mul(_dim_mul(_PA, _M), _M) == _N,
        "strain": _dim_div(_M, _M) == _DIMENSIONLESS,
        "area": _dim_mul(_M, _M) == _M2,
        "darcy_flux": _dim_div(_dim_mul(_M2, _PA), _dim_mul(_dim_mul(_PA, _S), _M)) == _dim_div(_M, _S),
        "flow_from_flux_area": _dim_mul(_dim_div(_M, _S), _M2) == _dim_mul(_M2, _dim_div(_M, _S)),
        "tewl_conversion": tewl_to_flux_m_s(1.0) == 1.0e-3 / (3600.0 * 1000.0),
    }
    if not all(bool(value) for value in checks.values()):
        raise AssertionError(f"unit check failed: {checks}")
    return {"status": "pass", "checks": checks}


def prereg_checks(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    baseline = simulate(parameters)
    lower_b = replace(parameters, b_m=parameters.b_m * 0.5)
    upper_b = replace(parameters, b_m=parameters.b_m * 1.5)
    lower_a = replace(parameters, a_m=parameters.a_m * 0.5)
    upper_a = replace(parameters, a_m=parameters.a_m * 1.5)
    placebo = replace(parameters, b_m=0.0)
    placebo_result = simulate(placebo)
    lower_b_result = simulate(lower_b)
    upper_b_result = simulate(upper_b)
    lower_a_result = simulate(lower_a)
    upper_a_result = simulate(upper_a)
    force = baseline["mechanics"]["total_force_n"]
    transport = baseline["transport"]
    evap_ratio = transport["cut_evap_flux_m_s"] / transport["intact_evap_flux_m_s"]
    sensitivity_result = sensitivity(parameters)
    sensitivity_rows = sensitivity_result["rows"]
    sensitivity_map = {
        (row["parameter"], row["factor"]): row for row in sensitivity_rows
    }
    sensitivity_finite = all(
        math.isfinite(float(row["total_force_n"]))
        and math.isfinite(float(row["cut_total_flow_m3_s"]))
        and math.isfinite(float(row["evaporation_ratio"]))
        for row in sensitivity_rows
    )
    tolerance = 1.0e-12
    checks = {
        "unit_check": unit_check()["status"] == "pass",
        "force_sanity_0p5_to_6_n": bool(0.5 <= force <= 6.0),
        "force_increases_with_b": bool(
            lower_b_result["mechanics"]["total_force_n"]
            < upper_b_result["mechanics"]["total_force_n"]
        ),
        "area_increases_with_b": bool(
            lower_b_result["geometry"]["exposed_area_m2"]
            < upper_b_result["geometry"]["exposed_area_m2"]
        ),
        "area_increases_with_a": bool(
            lower_a_result["geometry"]["exposed_area_m2"]
            < upper_a_result["geometry"]["exposed_area_m2"]
        ),
        "placebo_incremental_force_zero": bool(
            abs(placebo_result["mechanics"]["incremental_force_n"]) <= tolerance
        ),
        "placebo_area_zero": bool(abs(placebo_result["geometry"]["exposed_area_m2"]) <= tolerance),
        "placebo_cut_flow_zero": bool(abs(placebo_result["transport"]["cut_total_flow_m3_s"]) <= tolerance),
        "evaporation_analog_ratio": bool(
            abs(evap_ratio - parameters.barrier_damage_multiplier) <= tolerance
        ),
        "sensitivity_finite": bool(sensitivity_finite),
        "modulus_sensitivity_direction": bool(
            sensitivity_map[("modulus_pa", 0.5)]["total_force_n"]
            < force
            < sensitivity_map[("modulus_pa", 1.5)]["total_force_n"]
        ),
        "prestress_sensitivity_direction": bool(
            sensitivity_map[("prestress_pa", 0.5)]["total_force_n"]
            < force
            < sensitivity_map[("prestress_pa", 1.5)]["total_force_n"]
        ),
        "barrier_sensitivity_direction": bool(
            sensitivity_map[("barrier_damage_multiplier", 0.5)]["evaporation_ratio"]
            < evap_ratio
            < sensitivity_map[("barrier_damage_multiplier", 1.5)]["evaporation_ratio"]
        ),
    }
    checks["all_pass"] = all(bool(value) for value in checks.values())
    return checks


def sensitivity(
    parameters: Parameters = DEFAULT_PARAMETERS,
    names: tuple[str, ...] = ("modulus_pa", "prestress_pa", "barrier_damage_multiplier"),
    factors: tuple[float, ...] = (0.5, 1.5),
) -> dict[str, Any]:
    parameters.validate()
    if not names or not factors:
        raise ValueError("sensitivity needs names and factors")
    baseline = simulate(parameters)
    baseline_force = baseline["mechanics"]["total_force_n"]
    baseline_flow = baseline["transport"]["cut_total_flow_m3_s"]
    baseline_evap_ratio = (
        baseline["transport"]["cut_evap_flux_m_s"]
        / baseline["transport"]["intact_evap_flux_m_s"]
    )
    rows: list[dict[str, Any]] = []
    for name in names:
        if not hasattr(parameters, name):
            raise ValueError(f"unknown sensitivity parameter: {name}")
        for factor in factors:
            if not math.isfinite(factor) or factor <= 0.0:
                raise ValueError("sensitivity factors must be positive")
            value = float(getattr(parameters, name)) * factor
            altered = replace(parameters, **{name: value})
            result = simulate(altered)
            force_value = result["mechanics"]["total_force_n"]
            flow_value = result["transport"]["cut_total_flow_m3_s"]
            evap_ratio = (
                result["transport"]["cut_evap_flux_m_s"]
                / result["transport"]["intact_evap_flux_m_s"]
            )
            rows.append(
                {
                    "parameter": name,
                    "factor": float(factor),
                    "value": value,
                    "total_force_n": force_value,
                    "incremental_force_n": result["mechanics"]["incremental_force_n"],
                    "exposed_area_m2": result["geometry"]["exposed_area_m2"],
                    "cut_total_flow_m3_s": flow_value,
                    "evaporation_ratio": evap_ratio,
                    "relative_force_change": (force_value - baseline_force) / baseline_force,
                    "relative_flow_change": (flow_value - baseline_flow) / baseline_flow,
                    "relative_evaporation_ratio_change": (evap_ratio - baseline_evap_ratio)
                    / baseline_evap_ratio,
                }
            )
    return {"baseline": baseline, "rows": rows}


def run_model(parameters: Parameters = DEFAULT_PARAMETERS) -> dict[str, Any]:
    result = simulate(parameters)
    return {
        "id": "BT-HX-Q033",
        "model_version": "0.1.0",
        "equations": {
            "opening_width": "g(y)=2*b*sqrt(1-(y/a)^2)",
            "opening_stress": "sigma=sigma0+E*expm1(beta*epsilon)/beta",
            "closure_force": "F=integral[-a,a] h*sigma(g(y)/L_ref) dy",
            "exposed_area": "A_exposed=pi*a*b",
            "darcy_flux": "J_d=k*delta_p/(mu*L)",
            "evaporation": "J_evap,cut=m_evap*J_evap,intact",
            "total_flow": "Q_total=Q_darcy+Q_evap+Q_exudate",
            "tip_stress_factor": "K_tip=1+2*sqrt(a/r_tip)",
        },
        "parameters": parameter_table(parameters),
        "reference_values": [
            {
                "doi": "10.1055/a-2150-0587",
                "location": "Table 1 and Fig. 7",
                "quantity": "experienced-surgeon direct-closure upper threshold",
                "value": 5.4,
                "unit": "N",
            },
            {
                "doi": "10.1055/a-2150-0587",
                "location": "Table 1 and Fig. 7",
                "quantity": "experienced-surgeon minimum tension without direct-closure recommendation",
                "value": 6.0,
                "unit": "N",
            },
            {
                "doi": "10.1055/a-2150-0587",
                "location": "Table 1",
                "quantity": "fresh-cadaver pooled direct-closure upper threshold",
                "value": 4.2,
                "unit": "N",
            },
            {
                "doi": "10.3389/fmed.2024.1481645",
                "location": "Fig. 2",
                "quantity": "TEWL increase after mechanical dermabrasion",
                "value": 5.9,
                "unit": "fold relative to control; TEWL measured in g m^-2 h^-1",
            },
            {
                "doi": "10.3389/fmed.2024.1481645",
                "location": "Fig. 4",
                "quantity": "Lucifer Yellow receptor concentration increase after dermabrasion",
                "value": 133.0,
                "unit": "fold relative to control",
            },
        ],
        "default_case": result,
        "prereg_checks": prereg_checks(parameters),
        "sensitivity": sensitivity(parameters),
        "sources": [
            {
                "doi": "10.1055/a-2150-0587",
                "use": "verified force reference; Table 1 and Fig. 7",
            },
            {
                "doi": "10.3389/fmed.2024.1481645",
                "use": "verified 5.9-fold TEWL barrier-damage analog; Fig. 2 and Fig. 4",
            },
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default=str(Path(__file__).with_name("results.json")),
    )
    args = parser.parse_args(argv)
    payload = run_model()
    output_path = Path(args.output)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "id": payload["id"],
                "total_force_n": payload["default_case"]["mechanics"]["total_force_n"],
                "exposed_area_m2": payload["default_case"]["geometry"]["exposed_area_m2"],
                "cut_total_flow_m3_s": payload["default_case"]["transport"]["cut_total_flow_m3_s"],
                "all_prereg_checks_pass": payload["prereg_checks"]["all_pass"],
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
