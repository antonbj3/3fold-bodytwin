from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_PARAMS: dict[str, float | int] = {
    "n_regions": 9,
    "pb_mmHg": 760.0,
    "ph2o_mmHg": 47.0,
    "fio2": 0.21,
    "rest_rq": 0.80,
    "heavy_rq": 0.90,
    "rest_va_Lmin": 4.2,
    "heavy_va_Lmin": 30.0,
    "rest_q_Lmin": 5.0,
    "heavy_q_Lmin": 25.0,
    "rest_vco2_Lmin": 0.21,
    "heavy_vco2_Lmin": 1.60,
    "rest_pv_o2_ml_dl": 12.0,
    "heavy_pv_o2_ml_dl": 7.5,
    "hb_g_dl": 15.0,
    "o2_binding_ml_g": 1.34,
    "o2_dissolved_ml_dl_mmHg": 0.003,
    "p50_mmHg": 26.6,
    "hill_n": 2.7,
    "oxygen_gas_coeff_ml_L_mmHg": 1000.0 / 760.0,
    "co2_pressure_scale_mmHg": 760.0,
    "lung_height_m": 0.30,
    "rho_blood_kg_m3": 1060.0,
    "gravity_m_s2": 9.81,
    "rest_art_mean_mmHg": 15.0,
    "heavy_art_mean_mmHg": 25.0,
    "rest_ven_mean_mmHg": 7.0,
    "heavy_ven_mean_mmHg": 15.0,
    "art_hydro_fraction": 1.0,
    "ven_hydro_fraction": 0.35,
    "rv_top_relative": 1.0,
    "rv_base_gradient": 0.35,
    "exercise_recruitment": 0.25,
    "hpv_gain": 0.80,
    "hpv_reference_pao2_mmHg": 100.0,
    "aw_resistance_top_relative": 1.25,
    "aw_resistance_gradient": 0.35,
    "compliance_base_gradient": 0.30,
    "exercise_bronchodilation": 0.55,
    "rest_driving_pressure_cm_h2o": 8.0,
    "heavy_driving_pressure_cm_h2o": 14.0,
    "max_fixed_point_iterations": 80,
    "fixed_point_tolerance": 1.0e-9,
}


PARAMETER_TABLE: list[dict[str, str]] = [
    {"name": "pb_mmHg", "unit": "mmHg", "source_or_assumption": "standard barometric pressure assumption"},
    {"name": "ph2o_mmHg", "unit": "mmHg", "source_or_assumption": "standard 37 C water-vapor pressure assumption"},
    {"name": "fio2", "unit": "fraction", "source_or_assumption": "room-air inspired fraction assumption"},
    {"name": "rest_rq, heavy_rq", "unit": "dimensionless", "source_or_assumption": "standard respiratory-exchange-ratio range; load interpolation assumption"},
    {"name": "rest_va_Lmin, heavy_va_Lmin", "unit": "L/min", "source_or_assumption": "alveolar ventilation state assumptions; heavy value chosen to illustrate compensatory hyperventilation"},
    {"name": "rest_q_Lmin, heavy_q_Lmin", "unit": "L/min", "source_or_assumption": "pulmonary blood-flow state assumptions; not a subject measurement"},
    {"name": "rest_vco2_Lmin, heavy_vco2_Lmin", "unit": "L/min gas", "source_or_assumption": "metabolic CO2-output state assumptions"},
    {"name": "rest_pv_o2_ml_dl, heavy_pv_o2_ml_dl", "unit": "mL O2/dL", "source_or_assumption": "mixed-venous content assumptions for the bounded demonstration"},
    {"name": "hb_g_dl", "unit": "g/dL", "source_or_assumption": "standard adult hemoglobin concentration assumption"},
    {"name": "o2_binding_ml_g", "unit": "mL O2/g Hb", "source_or_assumption": "standard oxygen-binding capacity assumption"},
    {"name": "o2_dissolved_ml_dl_mmHg", "unit": "mL O2/dL/mmHg", "source_or_assumption": "standard dissolved-oxygen solubility assumption"},
    {"name": "p50_mmHg", "unit": "mmHg", "source_or_assumption": "standard Hill-curve midpoint assumption"},
    {"name": "hill_n", "unit": "dimensionless", "source_or_assumption": "standard Hill-coefficient approximation"},
    {"name": "oxygen_gas_coeff_ml_L_mmHg", "unit": "mL gas/(L gas mmHg)", "source_or_assumption": "derived from 1000/760 at standard pressure"},
    {"name": "co2_pressure_scale_mmHg", "unit": "mmHg", "source_or_assumption": "ideal-gas conversion used for the reduced CO2 balance"},
    {"name": "lung_height_m", "unit": "m", "source_or_assumption": "upright lung-column geometric assumption"},
    {"name": "rho_blood_kg_m3", "unit": "kg/m3", "source_or_assumption": "blood-density assumption"},
    {"name": "gravity_m_s2", "unit": "m/s2", "source_or_assumption": "standard gravitational acceleration"},
    {"name": "rest_art_mean_mmHg, heavy_art_mean_mmHg", "unit": "mmHg", "source_or_assumption": "pulmonary arterial pressure assumptions"},
    {"name": "rest_ven_mean_mmHg, heavy_ven_mean_mmHg", "unit": "mmHg", "source_or_assumption": "pulmonary venous pressure assumptions"},
    {"name": "art_hydro_fraction, ven_hydro_fraction", "unit": "dimensionless", "source_or_assumption": "reduced hydrostatic pressure-distribution assumption"},
    {"name": "rv_top_relative, rv_base_gradient", "unit": "relative resistance", "source_or_assumption": "reduced extra-alveolar vascular-resistance assumption"},
    {"name": "exercise_recruitment", "unit": "dimensionless", "source_or_assumption": "exercise recruitment/dilation sensitivity assumption"},
    {"name": "hpv_gain, hpv_reference_pao2_mmHg", "unit": "dimensionless, mmHg", "source_or_assumption": "reduced HPV law; qualitative mechanism anchored to West 1969, DOI 10.1016/0034-5687(69)90071-1"},
    {"name": "aw_resistance_top_relative, aw_resistance_gradient", "unit": "relative resistance", "source_or_assumption": "reduced airway-resistance gradient assumption"},
    {"name": "compliance_base_gradient", "unit": "relative compliance", "source_or_assumption": "reduced regional compliance assumption"},
    {"name": "exercise_bronchodilation", "unit": "dimensionless", "source_or_assumption": "load-related airway-resistance flattening assumption"},
    {"name": "rest_driving_pressure_cm_h2o, heavy_driving_pressure_cm_h2o", "unit": "cmH2O", "source_or_assumption": "driving-pressure state assumptions; normalized regional fractions are reported"},
    {"name": "max_fixed_point_iterations, fixed_point_tolerance", "unit": "iteration, dimensionless", "source_or_assumption": "numerical solver settings frozen before the run"},
]


SOURCES: list[dict[str, str]] = [
    {
        "citation": "Hall ET et al., J Appl Physiol 2014;116:451-461",
        "doi": "10.1152/japplphysiol.00659.2013",
        "location": "Table 4",
        "anchor": "regional perfusion rest/exercise values in ml min^-1 ml^-1",
    },
    {
        "citation": "Tedjasaputra V et al., J Appl Physiol 2013;115:126-135",
        "doi": "10.1152/japplphysiol.00778.2012",
        "location": "Table 3",
        "anchor": "A-aDO2 rest/exercise values in Torr and arterial blood gases",
    },
    {
        "citation": "Harf A, Pratt T, Hughes JM, J Appl Physiol 1978;44:115-123",
        "doi": "10.1152/jappl.1978.44.1.115",
        "location": "indexed abstract",
        "anchor": "40-150% apical blood-flow increase at 50 W",
    },
    {
        "citation": "West JB, Respir Physiol 1969;7:88-110",
        "doi": "10.1016/0034-5687(69)90071-1",
        "location": "mechanistic context",
        "anchor": "ventilation-perfusion inequality and overall gas exchange",
    },
]


def _parameters(overrides: dict[str, float | int] | None = None) -> dict[str, float | int]:
    values = dict(DEFAULT_PARAMS)
    if overrides:
        unknown = sorted(set(overrides) - set(values))
        if unknown:
            raise ValueError(f"unknown parameters: {unknown}")
        values.update(overrides)
    if int(values["n_regions"]) < 3:
        raise ValueError("n_regions must be at least 3")
    if not 0.0 <= float(values["fio2"]) < 1.0:
        raise ValueError("fio2 must be in [0,1)")
    if float(values["rest_rq"]) <= 0.0 or float(values["heavy_rq"]) <= 0.0:
        raise ValueError("respiratory quotient must be positive")
    return values


def _lerp(a: float, b: float, fraction: float) -> float:
    return a + fraction * (b - a)


def _oxygen_content(po2_mmHg: float | np.ndarray, p: dict[str, float | int]) -> float | np.ndarray:
    po2 = np.asarray(po2_mmHg, dtype=float)
    p50 = float(p["p50_mmHg"])
    hill = float(p["hill_n"])
    saturation = np.power(po2, hill) / (np.power(po2, hill) + np.power(p50, hill))
    content = float(p["hb_g_dl"]) * float(p["o2_binding_ml_g"]) * saturation
    content = content + float(p["o2_dissolved_ml_dl_mmHg"]) * po2
    if np.ndim(po2_mmHg) == 0:
        return float(content)
    return content


def _oxygen_pressure_from_content(content_ml_dl: float, p: dict[str, float | int]) -> float:
    content = float(content_ml_dl)
    if content <= 0.0:
        return 0.0
    lo = 0.0
    hi = 300.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if float(_oxygen_content(mid, p)) < content:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _solve_alveolar_o2(
    alveolar_ventilation_Lmin: float,
    perfusion_Lmin: float,
    inspired_o2_mmHg: float,
    venous_o2_content_ml_dl: float,
    p: dict[str, float | int],
) -> tuple[float, float]:
    if perfusion_Lmin <= 0.0:
        return inspired_o2_mmHg, 0.0
    if alveolar_ventilation_Lmin <= 1.0e-12:
        alveolar_o2 = _oxygen_pressure_from_content(venous_o2_content_ml_dl, p)
        return alveolar_o2, 0.0
    lo = 0.0
    hi = inspired_o2_mmHg
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        blood_uptake = perfusion_Lmin * 10.0 * (float(_oxygen_content(mid, p)) - venous_o2_content_ml_dl)
        gas_delivery = float(p["oxygen_gas_coeff_ml_L_mmHg"]) * alveolar_ventilation_Lmin * (inspired_o2_mmHg - mid)
        if blood_uptake < gas_delivery:
            lo = mid
        else:
            hi = mid
    alveolar_o2 = 0.5 * (lo + hi)
    blood_uptake = perfusion_Lmin * 10.0 * (float(_oxygen_content(alveolar_o2, p)) - venous_o2_content_ml_dl)
    gas_delivery = float(p["oxygen_gas_coeff_ml_L_mmHg"]) * alveolar_ventilation_Lmin * (inspired_o2_mmHg - alveolar_o2)
    return alveolar_o2, blood_uptake - gas_delivery


def _state_values(p: dict[str, float | int], exercise_fraction: float) -> dict[str, float]:
    e = float(exercise_fraction)
    return {
        "va_total": _lerp(float(p["rest_va_Lmin"]), float(p["heavy_va_Lmin"]), e),
        "q_total": _lerp(float(p["rest_q_Lmin"]), float(p["heavy_q_Lmin"]), e),
        "vco2": _lerp(float(p["rest_vco2_Lmin"]), float(p["heavy_vco2_Lmin"]), e),
        "rq": _lerp(float(p["rest_rq"]), float(p["heavy_rq"]), e),
        "pv_o2": _lerp(float(p["rest_pv_o2_ml_dl"]), float(p["heavy_pv_o2_ml_dl"]), e),
        "art_mean": _lerp(float(p["rest_art_mean_mmHg"]), float(p["heavy_art_mean_mmHg"]), e),
        "ven_mean": _lerp(float(p["rest_ven_mean_mmHg"]), float(p["heavy_ven_mean_mmHg"]), e),
        "driving_pressure": _lerp(float(p["rest_driving_pressure_cm_h2o"]), float(p["heavy_driving_pressure_cm_h2o"]), e),
    }


def _regional_coordinates(p: dict[str, float | int]) -> np.ndarray:
    n = int(p["n_regions"])
    return (np.arange(n, dtype=float) + 0.5) / n


def _ventilation_distribution(
    p: dict[str, float | int], coordinates: np.ndarray, exercise_fraction: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    e = float(exercise_fraction)
    resistance_gradient = float(p["aw_resistance_gradient"]) * (1.0 - float(p["exercise_bronchodilation"]) * e)
    resistance = float(p["aw_resistance_top_relative"]) * np.exp(-resistance_gradient * coordinates)
    compliance = 1.0 + float(p["compliance_base_gradient"]) * coordinates
    raw = float(_state_values(p, e)["driving_pressure"]) * compliance / resistance
    return raw / np.sum(raw), resistance, compliance


def _perfusion_distribution(
    p: dict[str, float | int],
    coordinates: np.ndarray,
    exercise_fraction: float,
    alveolar_o2_mmHg: np.ndarray,
    hpv_enabled: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    s = _state_values(p, exercise_fraction)
    mmhg_per_m = float(p["rho_blood_kg_m3"]) * float(p["gravity_m_s2"]) / 133.322
    hydrostatic = mmhg_per_m * float(p["lung_height_m"]) * (coordinates - 0.5)
    arterial = s["art_mean"] + float(p["art_hydro_fraction"]) * hydrostatic
    venous = s["ven_mean"] + float(p["ven_hydro_fraction"]) * hydrostatic
    driving_pressure = np.maximum(arterial - venous, 1.0e-6)
    resistance = float(p["rv_top_relative"]) * (
        1.0 + float(p["rv_base_gradient"]) * coordinates
        - float(p["exercise_recruitment"]) * float(exercise_fraction) * (1.0 - coordinates)
    )
    resistance = np.maximum(resistance, 0.05)
    if hpv_enabled:
        hypoxia_fraction = np.clip(
            (float(p["hpv_reference_pao2_mmHg"]) - alveolar_o2_mmHg)
            / float(p["hpv_reference_pao2_mmHg"]),
            0.0,
            1.0,
        )
        resistance = resistance * (1.0 + float(p["hpv_gain"]) * hypoxia_fraction)
    raw = driving_pressure / resistance
    return raw / np.sum(raw), driving_pressure, resistance, hypoxia_fraction if hpv_enabled else np.zeros_like(coordinates)


def _local_paco2(
    ventilation_fractions: np.ndarray,
    perfusion_fractions: np.ndarray,
    va_total: float,
    q_total: float,
    global_paco2: float,
) -> np.ndarray:
    global_vq = 1.0  # local_vq below is a ratio of normalized fractions
    local_vq = np.divide(
        ventilation_fractions,
        perfusion_fractions,
        out=np.full_like(ventilation_fractions, np.inf),
        where=perfusion_fractions > 1.0e-14,
    )
    ratio = np.divide(
        global_vq,
        local_vq,
        out=np.ones_like(local_vq),
        where=np.isfinite(local_vq),
    )
    return np.minimum(global_paco2 * ratio, 200.0)


def simulate(
    params: dict[str, float | int] | None = None,
    exercise_fraction: float = 0.0,
    regional: bool = True,
    hpv_enabled: bool = True,
) -> dict[str, Any]:
    p = _parameters(params)
    e = float(exercise_fraction)
    if not 0.0 <= e <= 1.0:
        raise ValueError("exercise_fraction must be in [0,1]")
    s = _state_values(p, e)
    coordinates = _regional_coordinates(p)
    n = len(coordinates)
    ventilation_fractions, airway_resistance, compliance = _ventilation_distribution(p, coordinates, e)
    if regional:
        perfusion_fractions = np.full(n, 1.0 / n)
        alveolar_o2 = np.full(n, float(p["hpv_reference_pao2_mmHg"]))
        iterations = 0
        for iterations in range(1, int(p["max_fixed_point_iterations"]) + 1):
            perfusion_fractions, _, _, _ = _perfusion_distribution(p, coordinates, e, alveolar_o2, hpv_enabled)
            ventilation_Lmin = s["va_total"] * ventilation_fractions
            perfusion_Lmin = s["q_total"] * perfusion_fractions
            new_o2 = np.array(
                [
                    _solve_alveolar_o2(
                        float(ventilation_Lmin[i]),
                        float(perfusion_Lmin[i]),
                        float(p["fio2"]) * (float(p["pb_mmHg"]) - float(p["ph2o_mmHg"])),
                        s["pv_o2"],
                        p,
                    )[0]
                    for i in range(n)
                ]
            )
            if float(np.max(np.abs(new_o2 - alveolar_o2))) < float(p["fixed_point_tolerance"]):
                alveolar_o2 = new_o2
                break
            alveolar_o2 = 0.5 * alveolar_o2 + 0.5 * new_o2
        perfusion_fractions, _, _, hypoxia_fraction = _perfusion_distribution(p, coordinates, e, alveolar_o2, hpv_enabled)
    else:
        ventilation_fractions = np.full(n, 1.0 / n)
        perfusion_fractions = np.full(n, 1.0 / n)
        hypoxia_fraction = np.zeros(n)
        iterations = 0

    ventilation_Lmin = s["va_total"] * ventilation_fractions
    perfusion_Lmin = s["q_total"] * perfusion_fractions
    inspired_o2 = float(p["fio2"]) * (float(p["pb_mmHg"]) - float(p["ph2o_mmHg"]))
    roots: list[float] = []
    residuals: list[float] = []
    for i in range(n):
        root, residual = _solve_alveolar_o2(
            float(ventilation_Lmin[i]),
            float(perfusion_Lmin[i]),
            inspired_o2,
            s["pv_o2"],
            p,
        )
        roots.append(root)
        residuals.append(residual)
    alveolar_o2 = np.asarray(roots, dtype=float)
    alveolar_o2_content = np.asarray(_oxygen_content(alveolar_o2, p), dtype=float)
    mixed_arterial_content = float(np.dot(perfusion_fractions, alveolar_o2_content))
    mixed_arterial_o2 = _oxygen_pressure_from_content(mixed_arterial_content, p)
    global_paco2 = float(p["co2_pressure_scale_mmHg"]) * s["vco2"] / s["va_total"]
    alveolar_gas_equation_pao2 = inspired_o2 - global_paco2 * (
        float(p["fio2"]) + (1.0 - float(p["fio2"])) / s["rq"]
    )
    global_pao2 = float(np.dot(ventilation_fractions, alveolar_o2))
    paco2 = _local_paco2(
        ventilation_fractions,
        perfusion_fractions,
        s["va_total"],
        s["q_total"],
        global_paco2,
    )
    vq = np.divide(
        ventilation_fractions,
        perfusion_fractions,
        out=np.full(n, np.inf),
        where=perfusion_fractions > 1.0e-14,
    )
    top = slice(0, n // 3)
    middle = slice(n // 3, 2 * n // 3)
    bottom = slice(2 * n // 3, n)
    o2_uptake = float(
        np.dot(perfusion_Lmin, 10.0 * (alveolar_o2_content - s["pv_o2"]))
    )
    return {
        "regional": regional,
        "hpv_enabled": hpv_enabled,
        "exercise_fraction": e,
        "n_regions": n,
        "iterations": iterations,
        "coordinates": coordinates.tolist(),
        "ventilation_fraction": ventilation_fractions.tolist(),
        "perfusion_fraction": perfusion_fractions.tolist(),
        "ventilation_Lmin": ventilation_Lmin.tolist(),
        "perfusion_Lmin": perfusion_Lmin.tolist(),
        "vq": vq.tolist(),
        "airway_resistance_relative": airway_resistance.tolist(),
        "compliance_relative": compliance.tolist(),
        "hypoxia_fraction": hypoxia_fraction.tolist(),
        "alveolar_o2_mmHg": alveolar_o2.tolist(),
        "alveolar_o2_content_ml_dl": alveolar_o2_content.tolist(),
        "alveolar_co2_mmHg": paco2.tolist(),
        "global_paco2_mmHg": global_paco2,
        "global_pao2_mmHg": global_pao2,
        "alveolar_gas_equation_pao2_mmHg": alveolar_gas_equation_pao2,
        "gas_equation_consistency_gap_mmHg": alveolar_gas_equation_pao2 - global_pao2,
        "mixed_arterial_o2_mmHg": mixed_arterial_o2,
        "mixed_arterial_o2_content_ml_dl": mixed_arterial_content,
        "a_a_do2_mmHg": global_pao2 - mixed_arterial_o2,
        "o2_uptake_ml_min": o2_uptake,
        "gas_exchange_residual_ml_min": residuals,
        "regional_means": {
            "nondependent_vq_mean": float(np.mean(vq[top])),
            "middle_vq_mean": float(np.mean(vq[middle])),
            "dependent_vq_mean": float(np.mean(vq[bottom])),
            "nondependent_perfusion_mean": float(np.mean(perfusion_fractions[top])),
            "middle_perfusion_mean": float(np.mean(perfusion_fractions[middle])),
            "dependent_perfusion_mean": float(np.mean(perfusion_fractions[bottom])),
            "nondependent_alveolar_o2_mean_mmHg": float(np.mean(alveolar_o2[top])),
            "middle_alveolar_o2_mean_mmHg": float(np.mean(alveolar_o2[middle])),
            "dependent_alveolar_o2_mean_mmHg": float(np.mean(alveolar_o2[bottom])),
        },
    }


def unit_check(result: dict[str, Any]) -> dict[str, Any]:
    v = np.asarray(result["ventilation_fraction"], dtype=float)
    q = np.asarray(result["perfusion_fraction"], dtype=float)
    va = np.asarray(result["ventilation_Lmin"], dtype=float)
    qabs = np.asarray(result["perfusion_Lmin"], dtype=float)
    residual = np.asarray(result["gas_exchange_residual_ml_min"], dtype=float)
    checks = {
        "ventilation_fractions_sum_to_one": bool(abs(float(np.sum(v)) - 1.0) < 1.0e-12),
        "perfusion_fractions_sum_to_one": bool(abs(float(np.sum(q)) - 1.0) < 1.0e-12),
        "ventilation_nonnegative_Lmin": bool(np.all(va >= -1.0e-12)),
        "perfusion_nonnegative_Lmin": bool(np.all(qabs >= -1.0e-12)),
        "gas_residual_below_tolerance_ml_min": bool(np.max(np.abs(residual)) < 1.0e-7),
        "all_pressures_finite_mmHg": bool(np.all(np.isfinite(np.asarray(result["alveolar_o2_mmHg"], dtype=float)))),
        "o2_content_nonnegative_ml_dl": bool(np.all(np.asarray(result["alveolar_o2_content_ml_dl"], dtype=float) >= 0.0)),
    }
    checks["all_pass"] = bool(all(checks.values()))
    return checks


def _sensitivity(
    base_params: dict[str, float | int], exercise_fraction: float, homogeneous_result: dict[str, Any]
) -> dict[str, Any]:
    names = ["hpv_gain", "rv_base_gradient", "aw_resistance_gradient"]
    output: dict[str, Any] = {}
    for name in names:
        base = float(base_params[name])
        output[name] = {}
        for factor in (0.5, 1.0, 1.5):
            value = base * factor
            result = simulate({name: value}, exercise_fraction=exercise_fraction, regional=True, hpv_enabled=True)
            output[name][f"{factor:.1f}x"] = {
                "value": value,
                "a_a_do2_mmHg": result["a_a_do2_mmHg"],
                "o2_uptake_ml_min": result["o2_uptake_ml_min"],
                "nondependent_to_dependent_perfusion": (
                    result["regional_means"]["nondependent_perfusion_mean"]
                    / result["regional_means"]["dependent_perfusion_mean"]
                ),
                "regional_minus_homogeneous_a_a_do2_mmHg": (
                    result["a_a_do2_mmHg"] - homogeneous_result["a_a_do2_mmHg"]
                ),
            }
    return output


def run_all(overrides: dict[str, float | int] | None = None) -> dict[str, Any]:
    p = _parameters(overrides)
    regional_rest = simulate(p, exercise_fraction=0.0, regional=True, hpv_enabled=True)
    homogeneous_rest = simulate(p, exercise_fraction=0.0, regional=False, hpv_enabled=False)
    regional_heavy = simulate(p, exercise_fraction=1.0, regional=True, hpv_enabled=True)
    homogeneous_heavy = simulate(p, exercise_fraction=1.0, regional=False, hpv_enabled=False)
    placebo_parameters = _parameters(
        {
            "art_hydro_fraction": 0.0,
            "ven_hydro_fraction": 0.0,
            "rv_base_gradient": 0.0,
            "exercise_recruitment": 0.0,
            "aw_resistance_gradient": 0.0,
            "compliance_base_gradient": 0.0,
            "hpv_gain": 0.0,
        }
    )
    placebo_rest = simulate(placebo_parameters, exercise_fraction=0.0, regional=True, hpv_enabled=False)
    placebo_heavy = simulate(placebo_parameters, exercise_fraction=1.0, regional=True, hpv_enabled=False)
    rest_delta = regional_rest["a_a_do2_mmHg"] - homogeneous_rest["a_a_do2_mmHg"]
    heavy_delta = regional_heavy["a_a_do2_mmHg"] - homogeneous_heavy["a_a_do2_mmHg"]
    rest_perfusion_ratio = (
        regional_rest["regional_means"]["nondependent_perfusion_mean"]
        / regional_rest["regional_means"]["dependent_perfusion_mean"]
    )
    heavy_perfusion_ratio = (
        regional_heavy["regional_means"]["nondependent_perfusion_mean"]
        / regional_heavy["regional_means"]["dependent_perfusion_mean"]
    )
    placebo_rest_error = max(
        abs(float(a) - float(b))
        for a, b in zip(
            [homogeneous_rest["mixed_arterial_o2_mmHg"], homogeneous_rest["global_pao2_mmHg"]],
            [placebo_rest["mixed_arterial_o2_mmHg"], placebo_rest["global_pao2_mmHg"]],
        )
    )
    placebo_heavy_error = max(
        abs(float(a) - float(b))
        for a, b in zip(
            [homogeneous_heavy["mixed_arterial_o2_mmHg"], homogeneous_heavy["global_pao2_mmHg"]],
            [placebo_heavy["mixed_arterial_o2_mmHg"], placebo_heavy["global_pao2_mmHg"]],
        )
    )
    checks = {
        "regional_rest": unit_check(regional_rest),
        "homogeneous_rest": unit_check(homogeneous_rest),
        "regional_heavy": unit_check(regional_heavy),
        "homogeneous_heavy": unit_check(homogeneous_heavy),
    }
    criteria = {
        "heavy_a_a_increment_at_least_10_mmHg": bool(
            regional_heavy["a_a_do2_mmHg"] - regional_rest["a_a_do2_mmHg"] >= 10.0
        ),
        "heavy_a_a_in_10_to_40_mmHg": bool(10.0 <= regional_heavy["a_a_do2_mmHg"] <= 40.0),
        "regional_need_at_least_5_mmHg": bool(abs(heavy_delta) >= 5.0),
        "perfusion_ratio_tolerance_0_20": bool(abs(heavy_perfusion_ratio - 0.784) <= 0.20),
        "placebo_equivalent_within_1e_12": bool(max(placebo_rest_error, placebo_heavy_error) <= 1.0e-12),
    }
    return {
        "id": "BT-HX-Q021",
        "model": "nine-region pressure-resistance ventilation, perfusion, HPV, and O2 mass-balance model",
        "implementation_history": [
            "The first executable diagnostic exposed a non-zero homogeneous A-a offset when the independent alveolar-gas-equation PO2 was compared with the regional O2 mass-balance root.",
            "The final comparator uses the ventilation-weighted regional mass-balance PO2; the independent alveolar gas equation remains in results.json as a consistency check. No model parameter or frozen criterion was changed."
        ],
        "equations": {
            "oxygen_content": "C_O2 = Hb * 1.34 * P_O2^n/(P_O2^n + P50^n) + 0.003 P_O2",
            "perfusion": "Q_i proportional to max(P_art_i - P_ven_i, 0) / R_v,i",
            "ventilation": "V_i proportional to C_i * driving_pressure / R_aw,i",
            "hpv": "R_v,i multiplied by 1 + gain * clip((P_AO2_ref - P_AO2_i)/P_AO2_ref, 0, 1)",
            "regional_o2_balance": "Q_i(C_a,i - C_v) = 1000/760 * V_A,i(P_IO2 - P_AO2_i)",
            "global_alveolar_o2": "P_AO2_global = sum(V_i P_AO2_i)/sum(V_i)",
            "alveolar_gas_equation_check": "P_AO2_gas = P_IO2 - P_ACO2[F_IO2 + (1-F_IO2)/RQ]; reported as a consistency check, not used as the regional comparator",
        },
        "parameters": p,
        "parameter_table": PARAMETER_TABLE,
        "sources": SOURCES,
        "reference_anchors": [
            {
                "source": "Hall ET et al., DOI 10.1152/japplphysiol.00659.2013",
                "location": "Table 4",
                "quantity": "regional perfusion",
                "unit": "ml min^-1 ml^-1",
                "rest": {"nondependent": 2.9, "middle": 4.5, "dependent": 4.2, "sd": {"nondependent": 1.7, "middle": 1.7, "dependent": 1.7}},
                "exercise": {"nondependent": 4.0, "middle": 5.7, "dependent": 5.1, "sd": {"nondependent": 1.5, "middle": 2.1, "dependent": 1.6}},
                "derived_nondependent_to_dependent_ratio": {"rest": 0.690, "exercise": 0.784}
            },
            {
                "source": "Tedjasaputra V et al., DOI 10.1152/japplphysiol.00778.2012",
                "location": "Table 3",
                "quantity": "A-aDO2",
                "unit": "Torr",
                "rest": 6.3,
                "exercise": 23.3,
                "sd": {"rest": 3.7, "exercise": 5.3}
            },
            {
                "source": "Harf A et al., DOI 10.1152/jappl.1978.44.1.115",
                "location": "indexed abstract",
                "quantity": "apical blood-flow increase",
                "unit": "percent",
                "value": "40-150",
                "exercise": "50 W"
            }
        ],
        "cases": {
            "regional_rest": regional_rest,
            "homogeneous_rest": homogeneous_rest,
            "regional_heavy": regional_heavy,
            "homogeneous_heavy": homogeneous_heavy,
        },
        "comparisons": {
            "rest_regional_minus_homogeneous_a_a_do2_mmHg": rest_delta,
            "heavy_regional_minus_homogeneous_a_a_do2_mmHg": heavy_delta,
            "rest_nondependent_to_dependent_perfusion": rest_perfusion_ratio,
            "heavy_nondependent_to_dependent_perfusion": heavy_perfusion_ratio,
            "placebo_rest_max_gas_error": placebo_rest_error,
            "placebo_heavy_max_gas_error": placebo_heavy_error,
        },
        "criteria": criteria,
        "unit_checks": checks,
        "sensitivity": _sensitivity(p, 1.0, homogeneous_heavy),
    }


def write_results(path: str | Path = "results.json", overrides: dict[str, float | int] | None = None) -> dict[str, Any]:
    result = run_all(overrides)
    output = Path(path)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    result = write_results()
    summary = {
        "id": result["id"],
        "rest_a_a_mmHg": result["cases"]["regional_rest"]["a_a_do2_mmHg"],
        "heavy_a_a_mmHg": result["cases"]["regional_heavy"]["a_a_do2_mmHg"],
        "heavy_regional_minus_homogeneous_mmHg": result["comparisons"]["heavy_regional_minus_homogeneous_a_a_do2_mmHg"],
        "criteria": result["criteria"],
    }
    print(json.dumps(summary, indent=2))
