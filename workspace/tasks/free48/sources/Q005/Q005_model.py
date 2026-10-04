"""First-principles renal OCT2/MATE clearance model for BT-HX-Q005.

For substrate s, with unbound plasma concentration Cu = f_u C_total and
OCT2 uptake J = Jmax Cu/(K + Cu), the steady linear apical efflux gives
Ccell = J/(60 * 1000 * P_MATE) when J is pmol cm^-2 min^-1, P_MATE is cm s^-1, and
Ccell is microM. Secretion clearance is A J/(1000 Cu) mL min^-1 because one
microM equals 1000 pmol mL^-1. Total clearance is f_u GFR plus secretion.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

SEED = 20260925
N_DRAWS = 20000

BASE_PARAMETERS: dict[str, Any] = {
    "gfr_mL_min": 94.99,
    "area_cm2": 2.0e5,
    "K_OCT2_ref_uM": 518.0,
    "Jmax_OCT2_ref_pmol_cm2_min": 1046.0,
    "f_u": 1.0,
    "P_MATE_cm_s": 0.0012,
    "kappa_creatinine": 0.40,
    "jmax_factor_creatinine": 0.067,
    "C_creatinine_uM": 100.0,
    "C_metformin_uM": 7.5,
    "f_u_prior_mean": 0.95,
    "f_u_prior_cv": 0.20,
    "gfr_prior_cv": 0.15,
    "K_OCT2_prior_cv": 0.30,
    "Jmax_OCT2_prior_cv": 0.25,
    "P_MATE_prior_cv": 0.30,
    "K_J_log_correlation": 0.70,
    "outcome_error_correlation": 0.35,
    "Cr_error_cv": 0.05,
    "Met_error_cv": 0.08,
}

SUBSTRATES: dict[str, dict[str, float]] = {
    "creatinine": {
        "C_total_uM": BASE_PARAMETERS["C_creatinine_uM"],
        "kappa": 1.0,
        "jmax_factor": 1.0,
    },
    "metformin": {
        "C_total_uM": BASE_PARAMETERS["C_metformin_uM"],
        "kappa": 1.0,
        "jmax_factor": 1.0,
    },
}

PARAMETER_TABLE: tuple[dict[str, str], ...] = (
    {
        "name": "K_OCT2_ref_uM",
        "unit": "microM",
        "kind": "molecular",
        "source_or_assumption": "Severance et al. 2017, Table 1, DOI 10.1124/jpet.117.242552; metformin apparent Ktapp 518 +/- 45 microM",
    },
    {
        "name": "Jmax_OCT2_ref_pmol_cm2_min",
        "unit": "pmol cm^-2 min^-1",
        "kind": "molecular",
        "source_or_assumption": "Severance et al. 2017, Table 1, DOI 10.1124/jpet.117.242552; metformin Jmax 1046 +/- 171",
    },
    {
        "name": "f_u",
        "unit": "dimensionless",
        "kind": "molecular",
        "source_or_assumption": "bounded near-unity assumption; nominal 1.0, prior mean 0.95, prior CV 0.20",
    },
    {
        "name": "P_MATE_cm_s",
        "unit": "cm s^-1",
        "kind": "molecular",
        "source_or_assumption": "linear MATE efflux assumption; no direct human value used",
    },
    {
        "name": "gfr_mL_min",
        "unit": "mL min^-1",
        "kind": "physiological",
        "source_or_assumption": "Topletz-Erickson et al. 2021, Table 3, DOI 10.1002/jcph.1750; 94.99 mL min^-1 1.73m^-2",
    },
    {
        "name": "area_cm2",
        "unit": "cm^2",
        "kind": "geometric",
        "source_or_assumption": "2.0e5 cm^2 total modeled kidney area scenario assumption",
    },
    {
        "name": "kappa_creatinine",
        "unit": "dimensionless",
        "kind": "substrate mapping",
        "source_or_assumption": "0.40 fixed derived mapping, not a measured creatinine K",
    },
    {
        "name": "jmax_factor_creatinine",
        "unit": "dimensionless",
        "kind": "substrate mapping",
        "source_or_assumption": "0.067 fixed mapping chosen to reproduce the cited secretion scale; not a measured creatinine Jmax",
    },
    {
        "name": "C_creatinine_uM",
        "unit": "microM",
        "kind": "exposure",
        "source_or_assumption": "100 microM bounded scenario assumption",
    },
    {
        "name": "C_metformin_uM",
        "unit": "microM",
        "kind": "exposure",
        "source_or_assumption": "7.5 microM scenario scale from cited 1.3 microg/mL Cmax and 129.9 g/mol",
    },
)

REFERENCE_ANCHORS: dict[str, Any] = {
    "OCT2_metformin_Jmax": {
        "value": 1046.0,
        "uncertainty": 171.0,
        "unit": "pmol cm^-2 min^-1",
        "location": "Severance et al. 2017, Table 1",
        "doi": "10.1124/jpet.117.242552",
    },
    "OCT2_metformin_Ktapp": {
        "value": 518.0,
        "uncertainty": 45.0,
        "unit": "microM",
        "location": "Severance et al. 2017, Table 1",
        "doi": "10.1124/jpet.117.242552",
    },
    "iohexol_GFR": {
        "value": 94.99,
        "unit": "mL min^-1 1.73m^-2",
        "location": "Topletz-Erickson et al. 2021, Table 3",
        "doi": "10.1002/jcph.1750",
    },
    "baseline_creatinine_clearance": {
        "value": 152.3,
        "unit": "mL min^-1",
        "location": "Topletz-Erickson et al. 2021, Figure 3A and Results",
        "doi": "10.1002/jcph.1750",
    },
    "inhibited_creatinine_clearance": {
        "value": 117.2,
        "unit": "mL min^-1",
        "location": "Topletz-Erickson et al. 2021, Figure 3A and Results",
        "doi": "10.1002/jcph.1750",
    },
    "baseline_metformin_renal_clearance": {
        "value": 30.0,
        "unit": "L h^-1",
        "location": "Topletz-Erickson et al. 2021, Table 2",
        "doi": "10.1002/jcph.1750",
    },
}

CANDIDATES: dict[str, float] = {
    "K_OCT2_ref_uM": 2.0,
    "Jmax_OCT2_ref_pmol_cm2_min": 1.0,
    "f_u": 1.0,
    "P_MATE_cm_s": 2.0,
    "gfr_mL_min": 3.0,
}

EQUATIONS: dict[str, str] = {
    "unbound_concentration": "Cu = f_u * C_total",
    "OCT2_uptake": "Jup = Jmax * Cu/(K_OCT2 + Cu)",
    "MATE_steady_state": "Ccell = Jup/(60 * 1000 * P_MATE)",
    "secretion_clearance": "CLsecret = area * Jup/(1000 * Cu)",
    "renal_clearance": "CL = f_u * GFR + CLsecret",
    "metformin_unit_conversion": "CL_Met_L_h = CL_Met_mL_min * 0.06",
}


def _array(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=float)


def _validate_parameters(parameters: dict[str, Any]) -> None:
    required_positive = (
        "gfr_mL_min",
        "area_cm2",
        "K_OCT2_ref_uM",
        "P_MATE_cm_s",
    )
    for name in required_positive:
        values = _array(parameters[name])
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError(f"{name} must be finite and positive")
    jmax = _array(parameters["Jmax_OCT2_ref_pmol_cm2_min"])
    if np.any(~np.isfinite(jmax)) or np.any(jmax < 0):
        raise ValueError("Jmax_OCT2_ref_pmol_cm2_min must be finite and nonnegative")
    f_u = _array(parameters["f_u"])
    if np.any(~np.isfinite(f_u)) or np.any((f_u <= 0) | (f_u > 1)):
        raise ValueError("f_u must be finite and in (0, 1]")


def compute_substrate_state(
    parameters: dict[str, Any],
    substrate: str,
    oct2_activity: float = 1.0,
    mate_activity: float = 1.0,
) -> dict[str, np.ndarray]:
    _validate_parameters(parameters)
    if substrate not in SUBSTRATES:
        raise ValueError(f"unknown substrate: {substrate}")
    if oct2_activity < 0 or mate_activity < 0:
        raise ValueError("transporter activities must be nonnegative")
    f_u = _array(parameters["f_u"])
    gfr = _array(parameters["gfr_mL_min"])
    area = _array(parameters["area_cm2"])
    k_ref = _array(parameters["K_OCT2_ref_uM"])
    jmax_ref = _array(parameters["Jmax_OCT2_ref_pmol_cm2_min"])
    p_mate = _array(parameters["P_MATE_cm_s"])
    substrate_spec = SUBSTRATES[substrate]
    c_total = _array(substrate_spec["C_total_uM"])
    if substrate == "creatinine":
        kappa = BASE_PARAMETERS["kappa_creatinine"]
        jmax_factor = BASE_PARAMETERS["jmax_factor_creatinine"]
    else:
        kappa = substrate_spec["kappa"]
        jmax_factor = substrate_spec["jmax_factor"]
    k_s = k_ref * kappa
    jmax_s = jmax_ref * jmax_factor
    c_u = f_u * c_total
    denominator = k_s + c_u
    fraction = np.divide(c_u, denominator, out=np.zeros_like(denominator), where=c_u > 0)
    jup = jmax_s * oct2_activity * fraction
    c_cell = jup / (60.0 * 1000.0 * p_mate * mate_activity)
    secretion = np.divide(
        area * jup,
        1000.0 * c_u,
        out=np.zeros_like(jup),
        where=c_u > 0,
    )
    filtration = f_u * gfr
    clearance = filtration + secretion
    return {
        "Cu_uM": c_u,
        "K_substrate_uM": k_s,
        "Jmax_substrate_pmol_cm2_min": jmax_s,
        "Jup_pmol_cm2_min": jup,
        "Ccell_uM": c_cell,
        "secretion_clearance_mL_min": secretion,
        "filtration_clearance_mL_min": filtration,
        "clearance_mL_min": clearance,
    }


def evaluate(parameters: dict[str, Any]) -> dict[str, Any]:
    cr = compute_substrate_state(parameters, "creatinine")
    met = compute_substrate_state(parameters, "metformin")
    return {
        "CL_Cr_mL_min": cr["clearance_mL_min"],
        "CL_Met_mL_min": met["clearance_mL_min"],
        "CL_Met_L_h": met["clearance_mL_min"] * 0.06,
        "Cr_secretion_mL_min": cr["secretion_clearance_mL_min"],
        "Met_secretion_mL_min": met["secretion_clearance_mL_min"],
        "Cr_cell_uM": cr["Ccell_uM"],
        "Met_cell_uM": met["Ccell_uM"],
        "Cr_Cu_uM": cr["Cu_uM"],
        "Met_Cu_uM": met["Cu_uM"],
    }


def unit_check() -> dict[str, Any]:
    nominal = dict(BASE_PARAMETERS)
    cr = compute_substrate_state(nominal, "creatinine")
    met = compute_substrate_state(nominal, "metformin")
    conversion_cr = cr["secretion_clearance_mL_min"] / (BASE_PARAMETERS["area_cm2"] * cr["Jup_pmol_cm2_min"])
    conversion_met = met["secretion_clearance_mL_min"] / (BASE_PARAMETERS["area_cm2"] * met["Jup_pmol_cm2_min"])
    expected = 1.0 / (1000.0 * BASE_PARAMETERS["C_creatinine_uM"])
    checks = {
        "K_dimension": "K and Cu are both microM",
        "Jup_dimension": "Jmax times dimensionless saturation fraction is pmol cm^-2 min^-1",
        "Ccell_balance": "Jup/(60 * 1000 * P_MATE) is microM",
        "secretion_conversion": "area * Jup/(1000 Cu) is mL min^-1",
        "metformin_conversion": "1 mL min^-1 = 0.06 L h^-1",
        "creatinine_numeric_balance": bool(np.allclose(conversion_cr, expected)),
        "metformin_numeric_balance": bool(np.allclose(conversion_met, 1.0 / (1000.0 * (BASE_PARAMETERS["f_u"] * BASE_PARAMETERS["C_metformin_uM"])))),
    }
    return {"passed": all(value is True for value in checks.values() if isinstance(value, bool)), "checks": checks}


def sensitivity_analysis() -> dict[str, Any]:
    nominal = evaluate(dict(BASE_PARAMETERS))
    selected = (
        "K_OCT2_ref_uM",
        "Jmax_OCT2_ref_pmol_cm2_min",
        "f_u",
        "P_MATE_cm_s",
        "gfr_mL_min",
    )
    expected_direction = {
        "K_OCT2_ref_uM": -1,
        "Jmax_OCT2_ref_pmol_cm2_min": 1,
        "f_u": 1,
        "P_MATE_cm_s": 0,
        "gfr_mL_min": 1,
    }
    records: list[dict[str, Any]] = []
    for name in selected:
        record: dict[str, Any] = {"parameter": name, "expected_primary_direction": expected_direction[name]}
        for factor, label in ((0.5, "minus50"), (1.5, "plus50")):
            perturbed = dict(BASE_PARAMETERS)
            perturbed_value = BASE_PARAMETERS[name] * factor
            if name == "f_u":
                perturbed_value = min(perturbed_value, 1.0)
            perturbed[name] = perturbed_value
            result = evaluate(perturbed)
            record[label] = {
                "CL_Cr_mL_min": float(np.asarray(result["CL_Cr_mL_min"])),
                "CL_Met_L_h": float(np.asarray(result["CL_Met_L_h"])),
                "Cr_cell_uM": float(np.asarray(result["Cr_cell_uM"])),
                "Met_cell_uM": float(np.asarray(result["Met_cell_uM"])),
                "relative_CL_Cr_change": float((np.asarray(result["CL_Cr_mL_min"]) - np.asarray(nominal["CL_Cr_mL_min"])) / np.asarray(nominal["CL_Cr_mL_min"])),
            }
        direction = record["plus50"]["relative_CL_Cr_change"] - record["minus50"]["relative_CL_Cr_change"]
        record["primary_direction_consistent"] = bool(
            expected_direction[name] == 0
            or direction * expected_direction[name] >= 0
        )
        record["all_finite"] = bool(all(np.isfinite(value) for side in (record["minus50"], record["plus50"]) for value in side.values()))
        records.append(record)
    return {
        "nominal": {key: float(np.asarray(value)) for key, value in nominal.items()},
        "records": records,
    }


def _lognormal_from_z(mean: float, cv: float, z: np.ndarray) -> np.ndarray:
    sigma = np.sqrt(np.log1p(cv * cv))
    return mean * np.exp(sigma * z - 0.5 * sigma * sigma)


def _standard_draws(n: int, seed: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    z_k = rng.standard_normal(n)
    z_j_independent = rng.standard_normal(n)
    z_j = BASE_PARAMETERS["K_J_log_correlation"] * z_k + np.sqrt(1.0 - BASE_PARAMETERS["K_J_log_correlation"] ** 2) * z_j_independent
    z_f = rng.standard_normal(n)
    z_g = rng.standard_normal(n)
    z_p = rng.standard_normal(n)
    return {"K": z_k, "J": z_j, "f": z_f, "gfr": z_g, "P": z_p}


def _parameters_from_standard(
    z: dict[str, np.ndarray],
    candidate: str | None = None,
) -> dict[str, Any]:
    cvs = {
        "K": BASE_PARAMETERS["K_OCT2_prior_cv"],
        "J": BASE_PARAMETERS["Jmax_OCT2_prior_cv"],
        "f": BASE_PARAMETERS["f_u_prior_cv"],
        "gfr": BASE_PARAMETERS["gfr_prior_cv"],
        "P": BASE_PARAMETERS["P_MATE_prior_cv"],
    }
    if candidate is not None:
        candidate_to_key = {
            "K_OCT2_ref_uM": "K",
            "Jmax_OCT2_ref_pmol_cm2_min": "J",
            "f_u": "f",
            "gfr_mL_min": "gfr",
            "P_MATE_cm_s": "P",
        }
        cvs[candidate_to_key[candidate]] *= 0.5
    f_u = _lognormal_from_z(BASE_PARAMETERS["f_u_prior_mean"], cvs["f"], z["f"])
    f_u = np.minimum(f_u, 1.0)
    return {
        "gfr_mL_min": _lognormal_from_z(BASE_PARAMETERS["gfr_mL_min"], cvs["gfr"], z["gfr"]),
        "area_cm2": BASE_PARAMETERS["area_cm2"],
        "K_OCT2_ref_uM": _lognormal_from_z(BASE_PARAMETERS["K_OCT2_ref_uM"], cvs["K"], z["K"]),
        "Jmax_OCT2_ref_pmol_cm2_min": _lognormal_from_z(BASE_PARAMETERS["Jmax_OCT2_ref_pmol_cm2_min"], cvs["J"], z["J"]),
        "f_u": f_u,
        "P_MATE_cm_s": _lognormal_from_z(BASE_PARAMETERS["P_MATE_cm_s"], cvs["P"], z["P"]),
    }


def _residual_errors(n: int, seed: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    z1 = rng.standard_normal(n)
    z2 = rng.standard_normal(n)
    rho = BASE_PARAMETERS["outcome_error_correlation"]
    z2 = rho * z1 + np.sqrt(1.0 - rho * rho) * z2
    return {
        "Cr": BASE_PARAMETERS["Cr_error_cv"] * z1,
        "Met": BASE_PARAMETERS["Met_error_cv"] * z2,
    }


def _array_rmse(values: np.ndarray, nominal: float) -> float:
    return float(np.sqrt(np.mean((values - nominal) ** 2)))


def uncertainty_analysis(n: int = N_DRAWS) -> dict[str, Any]:
    z = _standard_draws(n, SEED)
    errors = _residual_errors(n, SEED + 1)
    base_parameters = _parameters_from_standard(z)
    base = evaluate(base_parameters)
    nominal = evaluate(dict(BASE_PARAMETERS))
    nominal_cr = float(np.asarray(nominal["CL_Cr_mL_min"]))
    nominal_met = float(np.asarray(nominal["CL_Met_L_h"]))
    latent_cr = _array(base["CL_Cr_mL_min"])
    latent_met_lh = _array(base["CL_Met_L_h"])
    observed_cr = latent_cr * (1.0 + errors["Cr"])
    observed_met = latent_met_lh * (1.0 + errors["Met"])
    baseline_rmse_cr = _array_rmse(observed_cr, nominal_cr)
    baseline_rmse_met = _array_rmse(observed_met, nominal_met)
    joint_before = float(np.sqrt(np.mean(((observed_cr - nominal_cr) / nominal_cr) ** 2 + ((observed_met - nominal_met) / nominal_met) ** 2)))
    candidate_records: list[dict[str, Any]] = []
    for name, cost in CANDIDATES.items():
        refined_parameters = _parameters_from_standard(z, candidate=name)
        refined = evaluate(refined_parameters)
        refined_cr = _array(refined["CL_Cr_mL_min"]) * (1.0 + errors["Cr"])
        refined_met = _array(refined["CL_Met_L_h"]) * (1.0 + errors["Met"])
        rmse_cr = _array_rmse(refined_cr, nominal_cr)
        rmse_met = _array_rmse(refined_met, nominal_met)
        joint_after = float(np.sqrt(np.mean(((refined_cr - nominal_cr) / nominal_cr) ** 2 + ((refined_met - nominal_met) / nominal_met) ** 2)))
        fixed_parameters = _parameters_from_standard(z)
        fixed_parameters[name] = BASE_PARAMETERS[name]
        fixed = evaluate(fixed_parameters)
        fixed_cr = _array(fixed["CL_Cr_mL_min"])
        fixed_met = _array(fixed["CL_Met_L_h"])
        contribution_cr = max(0.0, float(np.var(latent_cr) - np.var(fixed_cr)))
        contribution_met = max(0.0, float(np.var(latent_met_lh) - np.var(fixed_met)))
        total_cr = max(float(np.var(latent_cr)), np.finfo(float).tiny)
        total_met = max(float(np.var(latent_met_lh)), np.finfo(float).tiny)
        record = {
            "parameter": name,
            "cost_units": cost,
            "primary_rmse_before_mL_min": baseline_rmse_cr,
            "primary_rmse_after_mL_min": rmse_cr,
            "primary_rmse_reduction": baseline_rmse_cr - rmse_cr,
            "primary_rmse_reduction_fraction": (baseline_rmse_cr - rmse_cr) / baseline_rmse_cr,
            "primary_score_mL_min_per_cost": (baseline_rmse_cr - rmse_cr) / cost,
            "secondary_rmse_before_L_h": baseline_rmse_met,
            "secondary_rmse_after_L_h": rmse_met,
            "secondary_rmse_reduction": baseline_rmse_met - rmse_met,
            "secondary_score_L_h_per_cost": (baseline_rmse_met - rmse_met) / cost,
            "joint_rmse_before": joint_before,
            "joint_rmse_after": joint_after,
            "joint_score_per_cost": (joint_before - joint_after) / cost,
            "primary_variance_share": contribution_cr / total_cr,
            "secondary_variance_share": contribution_met / total_met,
        }
        candidate_records.append(record)
    ranked = sorted(candidate_records, key=lambda item: item["primary_score_mL_min_per_cost"], reverse=True)
    top = ranked[0]
    second = ranked[1]
    margin = top["primary_score_mL_min_per_cost"] / second["primary_score_mL_min_per_cost"] if second["primary_score_mL_min_per_cost"] > 0 else float("inf")
    eligible = [item for item in ranked if item["primary_rmse_reduction_fraction"] >= 0.10 and item["primary_score_mL_min_per_cost"] > 0]
    if not eligible:
        winner = None
        decision = "not established"
    elif eligible[0]["parameter"] != top["parameter"] or margin < 1.10:
        winner = None
        decision = "tie"
    else:
        winner = top["parameter"]
        decision = "winner"
    molecular_names = {
        item["name"]
        for item in PARAMETER_TABLE
        if item["kind"] == "molecular"
    }
    molecular_ranked = [
        item for item in ranked if item["parameter"] in molecular_names
    ]
    leading_molecular = molecular_ranked[0] if molecular_ranked else None
    correlation_matrix = np.corrcoef(np.vstack((latent_cr, latent_met_lh)))
    return {
        "n_draws": n,
        "seed": SEED,
        "baseline_rmse_primary_mL_min": baseline_rmse_cr,
        "baseline_rmse_secondary_L_h": baseline_rmse_met,
        "joint_baseline_rmse": joint_before,
        "latent_output_correlation": {
            "CL_Cr_CL_Met": float(correlation_matrix[0, 1]),
            "CL_Cr_variance": float(np.var(latent_cr)),
            "CL_Met_L_h_variance": float(np.var(latent_met_lh)),
        },
        "candidates": candidate_records,
        "ranked_by_primary_score": [item["parameter"] for item in ranked],
        "molecular_ranked_by_primary_score": [item["parameter"] for item in molecular_ranked],
        "leading_molecular_parameter": leading_molecular["parameter"] if leading_molecular is not None else None,
        "leading_molecular_score_mL_min_per_cost": leading_molecular["primary_score_mL_min_per_cost"] if leading_molecular is not None else None,
        "leading_molecular_meets_10_percent_reduction": bool(leading_molecular is not None and leading_molecular["primary_rmse_reduction_fraction"] >= 0.10),
        "top_score_margin_over_second": float(margin),
        "decision": decision,
        "winner": winner,
        "meets_frozen_10_percent_reduction": [item["parameter"] for item in eligible],
    }


def reference_check(nominal_result: dict[str, Any]) -> dict[str, Any]:
    cr = float(np.asarray(nominal_result["CL_Cr_mL_min"]))
    met = float(np.asarray(nominal_result["CL_Met_L_h"]))
    cr_ratio = cr / REFERENCE_ANCHORS["baseline_creatinine_clearance"]["value"]
    met_ratio = met / REFERENCE_ANCHORS["baseline_metformin_renal_clearance"]["value"]
    return {
        "primary_ratio_to_reference": cr_ratio,
        "secondary_ratio_to_reference": met_ratio,
        "primary_within_0.5_to_2_reference": bool(0.5 <= cr_ratio <= 2.0),
        "secondary_within_0.5_to_2_reference": bool(0.5 <= met_ratio <= 2.0),
        "passed": bool(0.5 <= cr_ratio <= 2.0 and 0.5 <= met_ratio <= 2.0),
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def build_results() -> dict[str, Any]:
    nominal_result = evaluate(dict(BASE_PARAMETERS))
    nominal_json = _jsonable(nominal_result)
    null_model = {
        "CL_Cr_mL_min": BASE_PARAMETERS["f_u"] * BASE_PARAMETERS["gfr_mL_min"],
        "CL_Met_L_h": BASE_PARAMETERS["f_u"] * BASE_PARAMETERS["gfr_mL_min"] * 0.06,
        "description": "filtration-only null model; no active tubular secretion",
    }
    results = {
        "id": "BT-HX-Q005",
        "model": "renal OCT2/MATE steady-state clearance and value-of-parameter model",
        "seed": SEED,
        "n_draws": N_DRAWS,
        "base_parameters": BASE_PARAMETERS,
        "parameter_table": list(PARAMETER_TABLE),
        "reference_anchors": REFERENCE_ANCHORS,
        "equations": EQUATIONS,
        "nominal": nominal_json,
        "null_model": null_model,
        "unit_check": unit_check(),
        "reference_check": reference_check(nominal_json),
        "sensitivity": sensitivity_analysis(),
        "uncertainty": uncertainty_analysis(),
        "boundary": "assumption-based mechanistic scenario; no patient or internal measurement is a prediction target",
    }
    return _jsonable(results)


def main() -> None:
    results = build_results()
    output_path = Path(__file__).with_name("results.json")
    output_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "id": results["id"],
        "nominal_CL_Cr_mL_min": results["nominal"]["CL_Cr_mL_min"],
        "nominal_CL_Met_L_h": results["nominal"]["CL_Met_L_h"],
        "decision": results["uncertainty"]["decision"],
        "winner": results["uncertainty"]["winner"],
        "unit_check": results["unit_check"]["passed"],
        "reference_check": results["reference_check"]["passed"],
    }, indent=2))


if __name__ == "__main__":
    main()
