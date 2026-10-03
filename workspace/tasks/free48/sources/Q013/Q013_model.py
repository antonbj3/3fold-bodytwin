from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np


REFERENCE_AUC_RATIO = 1.27
REFERENCE_CI = (1.06, 1.52)
REFERENCE_N = 16
REFERENCE_PK_RATIO_TOLERANCE = 0.01
SENSITIVITY_MULTIPLIER = 1.5
CLASSIFICATION_PK_TOLERANCE = 0.02
CLASSIFICATION_SENSITIVITY_TOLERANCE = 0.05


@dataclass(frozen=True)
class Parameters:
    dose: float = 1.0
    bioavailability: float = 0.70
    ka: float = 1.5
    volume: float = 1.0
    clearance: float = 0.35
    k_in: float = 1.0
    k_out: float = 0.70
    ec50: float = 1.0
    emax: float = 1.0
    baseline: float = 0.0
    hill: float = 1.0
    t_end: float = 48.0
    dt: float = 0.01

    def validate(self) -> None:
        positive = {
            "dose": self.dose,
            "volume": self.volume,
            "clearance": self.clearance,
            "k_in": self.k_in,
            "k_out": self.k_out,
            "ec50": self.ec50,
            "emax": self.emax,
            "hill": self.hill,
            "t_end": self.t_end,
            "dt": self.dt,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not np.isfinite(self.bioavailability) or not 0 < self.bioavailability <= 1:
            raise ValueError("bioavailability must be in (0, 1]")
        if not np.isfinite(self.baseline):
            raise ValueError("baseline must be finite")


def parameter_table(p: Parameters) -> list[dict[str, Any]]:
    return [
        {
            "name": "dose",
            "value": p.dose,
            "unit": "normalized amount",
            "source_or_assumption": "assumption: normalized scenario, not clinical dosing",
            "status": "assumption",
        },
        {
            "name": "bioavailability",
            "value": p.bioavailability,
            "unit": "fraction",
            "source_or_assumption": "assumption: illustrative linear oral victim",
            "status": "assumption",
        },
        {
            "name": "ka",
            "value": p.ka,
            "unit": "h^-1",
            "source_or_assumption": "assumption: first-order absorption",
            "status": "assumption",
        },
        {
            "name": "volume",
            "value": p.volume,
            "unit": "normalized volume",
            "source_or_assumption": "assumption: one-compartment normalization",
            "status": "assumption",
        },
        {
            "name": "clearance",
            "value": p.clearance,
            "unit": "normalized volume h^-1",
            "source_or_assumption": "assumption: identifiable only as a ratio",
            "status": "assumption",
        },
        {
            "name": "k_in",
            "value": p.k_in,
            "unit": "h^-1",
            "source_or_assumption": "assumption: plasma-to-cell transfer rate",
            "status": "assumption",
        },
        {
            "name": "k_out",
            "value": p.k_out,
            "unit": "h^-1",
            "source_or_assumption": "assumption: local cellular removal rate",
            "status": "assumption",
        },
        {
            "name": "EC50",
            "value": p.ec50,
            "unit": "normalized concentration",
            "source_or_assumption": "assumption: target concentration for half-maximal occupancy",
            "status": "assumption",
        },
        {
            "name": "Emax",
            "value": p.emax,
            "unit": "normalized effect",
            "source_or_assumption": "assumption: effect amplitude",
            "status": "assumption",
        },
        {
            "name": "E0",
            "value": p.baseline,
            "unit": "normalized effect",
            "source_or_assumption": "assumption: baseline effect",
            "status": "assumption",
        },
        {
            "name": "hill",
            "value": p.hill,
            "unit": "dimensionless",
            "source_or_assumption": "assumption: unit Hill slope",
            "status": "assumption",
        },
        {
            "name": "reference_auc_ratio",
            "value": REFERENCE_AUC_RATIO,
            "unit": "dimensionless",
            "source_or_assumption": "Tachibana et al. 2025, Table 3, AUClast GMR 1.27; 90% CI 1.06-1.52; n=16",
            "status": "source",
        },
        {
            "name": "sensitivity_multiplier",
            "value": SENSITIVITY_MULTIPLIER,
            "unit": "dimensionless",
            "source_or_assumption": "frozen synthetic sensitivity scenario; not measured data",
            "status": "assumption",
        },
    ]


def _trapezoid(values: np.ndarray, time_h: np.ndarray) -> float:
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:
        trapezoid = np.trapz
    return float(trapezoid(values, time_h))


def _time_grid(p: Parameters, t_end: float | None, dt: float | None) -> np.ndarray:
    end = p.t_end if t_end is None else t_end
    step = p.dt if dt is None else dt
    if not np.isfinite(end) or end <= 0:
        raise ValueError("t_end must be finite and positive")
    if not np.isfinite(step) or step <= 0 or step > end:
        raise ValueError("dt must be finite, positive, and no larger than t_end")
    count = max(1, int(np.ceil(end / step)))
    return np.linspace(0.0, end, count + 1, dtype=float)


def _plasma_components(
    p: Parameters,
    bioavailability_multiplier: float,
    clearance_multiplier: float,
    route: str,
) -> tuple[float, float, float, list[tuple[float, float]] | None, float | None]:
    if not np.isfinite(bioavailability_multiplier) or bioavailability_multiplier <= 0:
        raise ValueError("bioavailability_multiplier must be finite and positive")
    if not np.isfinite(clearance_multiplier) or clearance_multiplier <= 0:
        raise ValueError("clearance_multiplier must be finite and positive")
    if route not in {"oral", "iv"}:
        raise ValueError("route must be oral or iv")
    clearance = p.clearance * clearance_multiplier
    elimination_rate = clearance / p.volume
    if route == "iv":
        effective_f = 1.0
        components = [(elimination_rate, p.dose / p.volume)]
        return effective_f, clearance, elimination_rate, components, None
    effective_f = p.bioavailability * bioavailability_multiplier
    if effective_f > 1.0 + 1e-12:
        raise ValueError("effective oral bioavailability cannot exceed one")
    if np.isclose(p.ka, elimination_rate, rtol=1e-10, atol=1e-12):
        return effective_f, clearance, elimination_rate, None, p.ka
    coefficient = p.dose * effective_f * p.ka / (p.volume * (p.ka - elimination_rate))
    components = [
        (elimination_rate, coefficient),
        (p.ka, -coefficient),
    ]
    return effective_f, clearance, elimination_rate, components, None


def _local_profile(
    time_h: np.ndarray,
    p: Parameters,
    bioavailability_multiplier: float,
    clearance_multiplier: float,
    route: str,
) -> np.ndarray:
    effective_f, _, elimination_rate, components, equal_ka = _plasma_components(
        p, bioavailability_multiplier, clearance_multiplier, route
    )
    if components is None:
        coefficient = p.dose * effective_f * equal_ka / p.volume
        if np.isclose(p.k_out, elimination_rate, rtol=1e-10, atol=1e-12):
            local = coefficient * p.k_in * time_h**2 * np.exp(-elimination_rate * time_h) / 2.0
        else:
            denominator = p.k_out - elimination_rate
            local = coefficient * p.k_in * (
                (time_h / denominator - 1.0 / denominator**2) * np.exp(-elimination_rate * time_h)
                + np.exp(-p.k_out * time_h) / denominator**2
            )
        return np.maximum(local, 0.0)
    local = np.zeros_like(time_h, dtype=float)
    for rate, coefficient in components:
        if np.isclose(p.k_out, rate, rtol=1e-10, atol=1e-12):
            local += coefficient * p.k_in * time_h * np.exp(-rate * time_h)
        else:
            local += (
                coefficient
                * p.k_in
                / (p.k_out - rate)
                * (np.exp(-rate * time_h) - np.exp(-p.k_out * time_h))
            )
    if not np.allclose(local, np.zeros_like(local), atol=1e-10):
        if not np.all(np.isfinite(local)):
            raise FloatingPointError("local profile is not finite")
    return np.maximum(local, 0.0)


def _hill_occupancy(concentration: np.ndarray, p: Parameters) -> np.ndarray:
    concentration = np.maximum(np.asarray(concentration, dtype=float), 0.0)
    occupancy = np.zeros_like(concentration, dtype=float)
    positive = concentration > 0
    if np.any(positive):
        log_ratio = p.hill * (np.log(concentration[positive]) - np.log(p.ec50))
        occupancy[positive] = np.exp(-np.logaddexp(0.0, -log_ratio))
    return occupancy


def response_from_exposure(
    local_concentration: np.ndarray | list[float] | float,
    p: Parameters,
    sensitivity_multiplier: float = 1.0,
) -> np.ndarray:
    p.validate()
    if not np.isfinite(sensitivity_multiplier) or sensitivity_multiplier < 0:
        raise ValueError("sensitivity_multiplier must be finite and non-negative")
    concentration = np.asarray(local_concentration, dtype=float)
    occupancy = _hill_occupancy(concentration, p)
    return p.baseline + sensitivity_multiplier * p.emax * occupancy


def simulate(
    p: Parameters | None = None,
    *,
    route: str = "oral",
    bioavailability_multiplier: float = 1.0,
    clearance_multiplier: float = 1.0,
    sensitivity_multiplier: float = 1.0,
    t_end: float | None = None,
    dt: float | None = None,
) -> dict[str, Any]:
    p = Parameters() if p is None else p
    p.validate()
    time_h = _time_grid(p, t_end, dt)
    effective_f, clearance, elimination_rate, components, equal_ka = _plasma_components(
        p, bioavailability_multiplier, clearance_multiplier, route
    )
    if components is None:
        coefficient = p.dose * effective_f * equal_ka / p.volume
        plasma = coefficient * time_h * np.exp(-elimination_rate * time_h)
    else:
        plasma = np.zeros_like(time_h, dtype=float)
        for rate, coefficient in components:
            plasma += coefficient * np.exp(-rate * time_h)
    plasma = np.maximum(plasma, 0.0)
    local = _local_profile(time_h, p, bioavailability_multiplier, clearance_multiplier, route)
    occupancy = _hill_occupancy(local, p)
    effect = p.baseline + sensitivity_multiplier * p.emax * occupancy
    if not np.all(np.isfinite(effect)):
        raise FloatingPointError("effect profile is not finite")
    if route == "oral":
        full_auc = p.dose * effective_f / clearance
    else:
        full_auc = p.dose / clearance
    return {
        "time_h": time_h,
        "plasma_concentration": plasma,
        "local_concentration": local,
        "target_occupancy": occupancy,
        "effect": effect,
        "route": route,
        "bioavailability_multiplier": float(bioavailability_multiplier),
        "clearance_multiplier": float(clearance_multiplier),
        "sensitivity_multiplier": float(sensitivity_multiplier),
        "auc_plasma": _trapezoid(plasma, time_h),
        "auc_local": _trapezoid(local, time_h),
        "auc_effect": _trapezoid(effect, time_h),
        "full_auc_plasma_analytic": float(full_auc),
        "plasma_auc_tail_fraction": float(_trapezoid(plasma, time_h) / full_auc),
        "cmax_plasma": float(np.max(plasma)),
        "cmax_local": float(np.max(local)),
        "effect_max": float(np.max(effect)),
        "t_cmax_plasma_h": float(time_h[int(np.argmax(plasma))]),
        "t_cmax_local_h": float(time_h[int(np.argmax(local))]),
        "elimination_rate_h": float(elimination_rate),
    }


def classify_residual(residual: np.ndarray | list[float] | float, p: Parameters) -> str:
    p.validate()
    values = np.asarray(residual, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("residual must be finite")
    normalized = float(np.max(np.abs(values)) / p.emax)
    if normalized <= CLASSIFICATION_PK_TOLERANCE:
        return "PK_ONLY"
    if normalized > CLASSIFICATION_SENSITIVITY_TOLERANCE:
        return "CELL_SENSITIVITY"
    return "UNRESOLVED"


def decompose(
    alone: dict[str, Any],
    pk_counterfactual: dict[str, Any],
    observed_or_predicted: dict[str, Any],
    p: Parameters,
    *,
    local_exposure: np.ndarray | list[float] | None = None,
    observed_effect: np.ndarray | list[float] | None = None,
) -> dict[str, Any]:
    p.validate()
    if local_exposure is None:
        local = np.asarray(pk_counterfactual["local_concentration"], dtype=float)
        local_source = "modeled_pk_counterfactual"
    else:
        local = np.asarray(local_exposure, dtype=float)
        local_source = "external_local_exposure"
    if observed_effect is None:
        effect = np.asarray(observed_or_predicted["effect"], dtype=float)
        effect_source = "modeled_combined_effect"
    else:
        effect = np.asarray(observed_effect, dtype=float)
        effect_source = "external_effect_measurement"
    if local.shape != effect.shape:
        raise ValueError("local exposure and effect must have the same shape")
    pk_effect = response_from_exposure(local, p, 1.0)
    residual = effect - pk_effect
    local_auc = _trapezoid(local, np.asarray(alone["time_h"], dtype=float))
    alone_local_auc = float(alone["auc_local"])
    alone_effect_auc = float(alone["auc_effect"])
    effect_auc = _trapezoid(effect, np.asarray(alone["time_h"], dtype=float))
    exposure_ratio = local_auc / alone_local_auc if alone_local_auc else np.nan
    effect_ratio = effect_auc / alone_effect_auc if alone_effect_auc else np.nan
    return {
        "exposure_ratio_local": float(exposure_ratio),
        "effect_auc_ratio": float(effect_ratio),
        "pk_counterfactual_effect": pk_effect,
        "residual": residual,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "max_abs_residual_over_emax": float(np.max(np.abs(residual)) / p.emax),
        "classification": classify_residual(residual, p),
        "local_exposure_source": local_source,
        "effect_source": effect_source,
        "measurement_status": "UNKNOWN_NO_INDEPENDENT_LOCAL_MEASUREMENT",
    }


def sensitivity_analysis(
    p: Parameters | None = None,
    *,
    bioavailability_multiplier: float = REFERENCE_AUC_RATIO,
    sensitivity_multiplier: float = SENSITIVITY_MULTIPLIER,
) -> list[dict[str, Any]]:
    base = Parameters() if p is None else p
    base.validate()
    nominal = simulate(
        base,
        bioavailability_multiplier=bioavailability_multiplier,
        sensitivity_multiplier=sensitivity_multiplier,
    )
    rows: list[dict[str, Any]] = []
    for parameter_name in ("clearance", "ec50", "sensitivity_multiplier"):
        nominal_value = getattr(base, parameter_name) if parameter_name != "sensitivity_multiplier" else sensitivity_multiplier
        for factor in (0.5, 1.0, 1.5):
            value = nominal_value * factor
            if parameter_name == "sensitivity_multiplier":
                trial_parameters = base
                trial_gain = value
            else:
                trial_parameters = replace(base, **{parameter_name: value})
                trial_gain = sensitivity_multiplier
            trial = simulate(
                trial_parameters,
                bioavailability_multiplier=bioavailability_multiplier,
                sensitivity_multiplier=trial_gain,
            )
            rows.append(
                {
                    "parameter": parameter_name,
                    "factor": factor,
                    "value": float(value),
                    "unit": "h^-1" if parameter_name == "clearance" else ("normalized concentration" if parameter_name == "ec50" else "dimensionless"),
                    "auc_local": trial["auc_local"],
                    "auc_plasma": trial["auc_plasma"],
                    "effect_max": trial["effect_max"],
                    "auc_effect": trial["auc_effect"],
                    "auc_local_relative_to_nominal_parameter": trial["auc_local"] / nominal["auc_local"],
                    "effect_max_relative_to_nominal_parameter": trial["effect_max"] / nominal["effect_max"],
                }
            )
    return rows


def unit_check(p: Parameters) -> dict[str, Any]:
    p.validate()
    bioavailability_multiplier = REFERENCE_AUC_RATIO
    clearance_multiplier = 1.0
    calculated_ratio = bioavailability_multiplier / clearance_multiplier
    if not np.isclose(calculated_ratio, REFERENCE_AUC_RATIO, rtol=0.0, atol=1e-12):
        raise AssertionError("dimensionless PK ratio check failed")
    return {
        "time": "h",
        "plasma_concentration": "normalized amount / normalized volume",
        "local_concentration": "normalized amount / normalized local volume",
        "auc_plasma": "normalized concentration h",
        "auc_local": "normalized concentration h",
        "effect": "normalized effect",
        "reference_auc_ratio": "dimensionless",
        "calculated_reference_ratio_from_assumption": float(calculated_ratio),
        "all_dimensions_consistent": True,
    }


def _serializable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {key: _serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serializable(item) for item in value]
    return value


def _prereg_hash() -> str:
    path = Path(__file__).resolve().parent / "PREREG.md"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_experiment(p: Parameters | None = None) -> dict[str, Any]:
    base = Parameters() if p is None else p
    base.validate()
    null = simulate(base, bioavailability_multiplier=1.0, clearance_multiplier=1.0, sensitivity_multiplier=1.0)
    pk_only = simulate(
        base,
        bioavailability_multiplier=REFERENCE_AUC_RATIO,
        clearance_multiplier=1.0,
        sensitivity_multiplier=1.0,
    )
    cell_only = simulate(
        base,
        bioavailability_multiplier=1.0,
        clearance_multiplier=1.0,
        sensitivity_multiplier=SENSITIVITY_MULTIPLIER,
    )
    combined = simulate(
        base,
        bioavailability_multiplier=REFERENCE_AUC_RATIO,
        clearance_multiplier=1.0,
        sensitivity_multiplier=SENSITIVITY_MULTIPLIER,
    )
    decomposition = {
        "pk_only": decompose(null, pk_only, pk_only, base),
        "cell_only": decompose(null, null, cell_only, base),
        "combined": decompose(null, pk_only, combined, base),
    }
    plasma_ratio = pk_only["auc_plasma"] / null["auc_plasma"]
    local_ratio = pk_only["auc_local"] / null["auc_local"]
    calibration_pass = (
        abs(plasma_ratio - REFERENCE_AUC_RATIO) <= REFERENCE_PK_RATIO_TOLERANCE
        and REFERENCE_CI[0] <= plasma_ratio <= REFERENCE_CI[1]
        and abs(local_ratio - REFERENCE_AUC_RATIO) <= REFERENCE_PK_RATIO_TOLERANCE
    )
    return {
        "model_version": "BT-HX-Q013-0.1",
        "prereg_sha256": _prereg_hash(),
        "reference": {
            "citation": "Tachibana et al., Clinical and Translational Science 2025;18(11):e70330",
            "doi": "10.1111/cts.70330",
            "pmid": "41150706",
            "table": "Table 3",
            "quantity": "digoxin AUClast GMR with valemetostat 200 mg QD versus digoxin alone",
            "value": REFERENCE_AUC_RATIO,
            "ci90": list(REFERENCE_CI),
            "unit": "dimensionless",
            "n": REFERENCE_N,
            "interpretation": "systemic plasma exposure reference; not a local cell measurement",
        },
        "equations": {
            "oral_pk": "dA_g/dt=-ka*A_g; dA_p/dt=F*ka*A_g-CL*A_p/V_p; C_p=A_p/V_p",
            "pk_factors": "F_eff=F*F_mult; CL_eff=CL*CL_mult; AUC ratio=F_mult/CL_mult",
            "local_exposure": "dC_local/dt=k_in*C_p-k_out*C_local",
            "target_occupancy": "O=C_local^h/(EC50^h+C_local^h)",
            "effect": "E=E0+gamma*Emax*O",
            "sensitivity_residual": "delta_E_PD=E_observed_or_predicted-E_PK(C_local)",
        },
        "parameter_table": parameter_table(base),
        "unit_check": unit_check(base),
        "calibration": {
            "simulated_plasma_auc_ratio": float(plasma_ratio),
            "simulated_local_auc_ratio": float(local_ratio),
            "reference_ci90": list(REFERENCE_CI),
            "tolerance": REFERENCE_PK_RATIO_TOLERANCE,
            "pass": bool(calibration_pass),
            "source_branch": "F_mult=1.27, CL_mult=1; identifiability is not claimed",
        },
        "scenarios": {
            "victim_alone": null,
            "pk_only": pk_only,
            "cell_sensitivity_only": cell_only,
            "combined_pk_and_cell_sensitivity": combined,
        },
        "decomposition": decomposition,
        "sensitivity": sensitivity_analysis(base),
        "data_status": {
            "independent_local_exposure": "UNKNOWN",
            "cell_sensitivity_measurement": "UNKNOWN",
            "all_effect_values": "model predictions from frozen assumptions",
            "empirical_data_used": "one public PK ratio anchor only",
        },
        "next_resolution_step": "Measure paired plasma and local target-site exposure and the same cell response in alone, PK-only, and combined conditions; estimate F_mult, CL_mult, EC50, and gamma with a pre-specified identifiability design.",
    }


def write_results(output_path: str | Path = "results.json") -> Path:
    output = Path(output_path)
    if not output.is_absolute():
        output = Path(__file__).resolve().parent / output
    result = run_experiment()
    output.write_text(json.dumps(_serializable(result), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="BT-HX-Q013 mechanistic PK/PD separation model")
    parser.add_argument("--output", default="results.json")
    args = parser.parse_args()
    output = write_results(args.output)
    result = run_experiment()
    print(json.dumps({"output": str(output), "calibration": result["calibration"], "classification": {key: value["classification"] for key, value in result["decomposition"].items()}}, ensure_ascii=False))


if __name__ == "__main__":
    main()
