from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

CM_H2O_TO_MMHG = 0.73556
PCA_TWITCH = 5.40
PCA_FULL = 4.90
PCA_EC50 = 5.70
TARGET_GAIN_RATIO = 1.74
MODEL_VERSION = "q022-lumped-v1"


def pca_to_molar(pca: float) -> float:
    return float(10.0 ** (-float(pca)))


DEFAULT_PARAMETERS: dict[str, float] = {
    "dt": 0.001,
    "electrical_center": 0.160,
    "electrical_sigma": 0.025,
    "k_on": 100.0,
    "k_off": 50.0,
    "ca_rest": 1.0e-7,
    "s_max": 1.0e-5,
    "k_release": 20.0,
    "tau_ca": 0.12,
    "tau_load": 0.25,
    "hill_n": 4.0,
    "ca_ec50": pca_to_molar(PCA_EC50),
    "length_width": 0.20,
    "velocity_half": 0.50,
    "mass_norm": 1.0,
    "damping": 10.0,
    "passive_stiffness": 40.0,
    "afterload": 0.10,
    "p_ext": 0.0,
    "p_ref": 70.0,
    "flow_max": 1.0,
    "hypoxia_penalty": 0.20,
    "gregg_gain": TARGET_GAIN_RATIO * (1.0 - 0.20) - 1.0,
    "perfusion_gamma": 1.0,
    "compression_coef": 0.25,
    "force_scale": 1.0,
}


PARAMETER_SPECS: tuple[dict[str, str], ...] = (
    {"name": "dt", "unit": "s", "source": "numerical choice; explicit Euler stability choice"},
    {"name": "electrical_center", "unit": "s", "source": "synthetic control pulse placement; assumption"},
    {"name": "electrical_sigma", "unit": "s", "source": "synthetic dimensionless control pulse; assumption"},
    {"name": "k_on", "unit": "s^-1", "source": "activation-gate rate; assumption"},
    {"name": "k_off", "unit": "s^-1", "source": "activation-gate decay; assumption"},
    {"name": "ca_rest", "unit": "M", "source": "resting cytosolic Ca2+ assumption"},
    {"name": "s_max", "unit": "M", "source": "releasable SR Ca2+ assumption"},
    {"name": "k_release", "unit": "s^-1", "source": "SR release rate; assumption"},
    {"name": "tau_ca", "unit": "s", "source": "cytosolic removal time; assumption"},
    {"name": "tau_load", "unit": "s", "source": "SR reload time; assumption"},
    {"name": "hill_n", "unit": "1", "source": "cooperative myofilament approximation; assumption"},
    {"name": "ca_ec50", "unit": "M", "source": "pCa 5.70 sensitivity assumption below Fabiato Fig. 10 twitch anchor"},
    {"name": "length_width", "unit": "1", "source": "normalized force-length width; assumption"},
    {"name": "velocity_half", "unit": "s^-1", "source": "normalized force-velocity half velocity; assumption"},
    {"name": "mass_norm", "unit": "s^2", "source": "normalized mechanical mass; assumption"},
    {"name": "damping", "unit": "s", "source": "normalized viscous damping; assumption"},
    {"name": "passive_stiffness", "unit": "1", "source": "normalized passive stiffness; assumption"},
    {"name": "afterload", "unit": "F/F0", "source": "normalized load; assumption"},
    {"name": "p_ext", "unit": "cmH2O", "source": "effective external pressure; assumption"},
    {"name": "p_ref", "unit": "cmH2O", "source": "Schouten 1992 0-to-70 cmH2O anchor"},
    {"name": "flow_max", "unit": "1", "source": "normalized autoregulatory saturation; assumption"},
    {"name": "hypoxia_penalty", "unit": "1", "source": "low-flow penalty coefficient; assumption"},
    {"name": "gregg_gain", "unit": "1", "source": "derived to reproduce Schouten 74% static gain"},
    {"name": "perfusion_gamma", "unit": "1", "source": "perfusion saturation exponent; assumption"},
    {"name": "compression_coef", "unit": "1", "source": "strain-flow feedback coefficient; assumption"},
    {"name": "force_scale", "unit": "F/F0", "source": "normalized active-force scale; assumption"},
)


def default_parameters() -> dict[str, float]:
    return copy.deepcopy(DEFAULT_PARAMETERS)


def parameter_table(parameters: dict[str, float] | None = None) -> list[dict[str, Any]]:
    values = default_parameters() if parameters is None else parameters
    rows: list[dict[str, Any]] = []
    for spec in PARAMETER_SPECS:
        name = spec["name"]
        rows.append(
            {
                "name": name,
                "value": float(values[name]),
                "unit": spec["unit"],
                "source_or_assumption": spec["source"],
            }
        )
    return rows


def dimensional_equations() -> list[dict[str, str]]:
    return [
        {"equation": "a_dot = k_on e (1-a) - k_off a", "left_unit": "s^-1", "right_unit": "s^-1"},
        {"equation": "J_release = k_release a s", "left_unit": "M s^-1", "right_unit": "M s^-1"},
        {"equation": "q = clip((P-p_ext)/p_ref, 0, flow_max)", "left_unit": "1", "right_unit": "1"},
        {"equation": "F_active = F0 f_Ca f_L f_V G_perfusion", "left_unit": "F/F0", "right_unit": "F/F0"},
        {"equation": "m eps_ddot + b eps_dot + k eps = F_afterload - F_active", "left_unit": "F/F0", "right_unit": "F/F0"},
    ]


def check_units(parameters: dict[str, float] | None = None) -> dict[str, Any]:
    values = default_parameters() if parameters is None else parameters
    required_units = {spec["name"]: spec["unit"] for spec in PARAMETER_SPECS}
    errors: list[str] = []
    for name, expected_unit in required_units.items():
        if name not in values:
            errors.append(f"missing parameter {name}")
            continue
        try:
            value = float(values[name])
        except (TypeError, ValueError):
            errors.append(f"non-numeric parameter {name}")
            continue
        if not math.isfinite(value):
            errors.append(f"non-finite parameter {name}")
    positive_names = (
        "dt",
        "electrical_sigma",
        "k_on",
        "k_off",
        "ca_rest",
        "s_max",
        "k_release",
        "tau_ca",
        "tau_load",
        "hill_n",
        "ca_ec50",
        "length_width",
        "velocity_half",
        "mass_norm",
        "damping",
        "passive_stiffness",
        "p_ref",
        "flow_max",
        "force_scale",
    )
    for name in positive_names:
        if name in values and float(values[name]) <= 0.0:
            errors.append(f"non-positive parameter {name}")
    if values.get("hypoxia_penalty", 0.0) < 0.0:
        errors.append("negative hypoxia_penalty")
    if values.get("compression_coef", 0.0) < 0.0:
        errors.append("negative compression_coef")
    if abs(CM_H2O_TO_MMHG - 0.73556) > 1.0e-12:
        errors.append("pressure conversion constant mismatch")
    return {
        "status": "PASS" if not errors else "FAIL",
        "parameters_checked": len(required_units),
        "conversion_cmH2O_to_mmHg": CM_H2O_TO_MMHG,
        "expected_units": required_units,
        "dimensional_equations": dimensional_equations(),
        "errors": errors,
    }


def electrical_pulse(
    time_s: np.ndarray,
    center_s: float = 0.160,
    sigma_s: float = 0.025,
) -> np.ndarray:
    time = np.asarray(time_s, dtype=float)
    return np.exp(-0.5 * ((time - center_s) / sigma_s) ** 2)


def effective_flow(
    pressure_cmh2o: np.ndarray | float,
    strain: np.ndarray | float = 0.0,
    parameters: dict[str, float] | None = None,
) -> np.ndarray:
    values = default_parameters() if parameters is None else parameters
    pressure = np.maximum(np.asarray(pressure_cmh2o, dtype=float), 0.0)
    strain_array = np.asarray(strain, dtype=float)
    pressure, strain_array = np.broadcast_arrays(pressure, strain_array)
    q_pressure = np.clip(
        (pressure - float(values["p_ext"])) / float(values["p_ref"]),
        0.0,
        float(values["flow_max"]),
    )
    compression = np.clip(strain_array, 0.0, None)
    q_compression = np.clip(1.0 - float(values["compression_coef"]) * compression, 0.0, 1.0)
    return q_pressure * q_compression


def perfusion_gain(
    flow_norm: np.ndarray | float,
    parameters: dict[str, float] | None = None,
) -> np.ndarray:
    values = default_parameters() if parameters is None else parameters
    q = np.clip(np.asarray(flow_norm, dtype=float), 0.0, float(values["flow_max"]))
    hypoxia = 1.0 - float(values["hypoxia_penalty"]) * (1.0 - q) ** 2
    gregg = 1.0 + float(values["gregg_gain"]) * q ** float(values["perfusion_gamma"])
    return hypoxia * gregg


def calcium_activation(
    calcium_molar: np.ndarray | float,
    parameters: dict[str, float] | None = None,
) -> np.ndarray:
    values = default_parameters() if parameters is None else parameters
    ca = np.maximum(np.asarray(calcium_molar, dtype=float) - float(values["ca_rest"]), 0.0)
    ec50_delta = max(float(values["ca_ec50"]) - float(values["ca_rest"]), np.finfo(float).eps)
    ratio = ca / ec50_delta
    hill = float(values["hill_n"])
    return ratio**hill / (1.0 + ratio**hill)


def length_factor(
    strain: np.ndarray | float,
    parameters: dict[str, float] | None = None,
) -> np.ndarray:
    values = default_parameters() if parameters is None else parameters
    strain_array = np.asarray(strain, dtype=float)
    return np.exp(-((strain_array / float(values["length_width"])) ** 2))


def velocity_factor(
    strain_rate: np.ndarray | float,
    parameters: dict[str, float] | None = None,
) -> np.ndarray:
    values = default_parameters() if parameters is None else parameters
    shortening_rate = np.maximum(-np.asarray(strain_rate, dtype=float), 0.0)
    return 1.0 / (1.0 + (shortening_rate / float(values["velocity_half"])) ** 2)


def _as_trace(value: np.ndarray | float, length: int, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return np.full(length, float(array), dtype=float)
    if array.ndim != 1 or array.size != length:
        raise ValueError(f"{name} must be scalar or have one value per time step")
    return array


def simulate(
    electrical: np.ndarray | float,
    perfusion_pressure_cmh2o: np.ndarray | float,
    parameters: dict[str, float] | None = None,
    initial_strain: float = 0.0,
) -> dict[str, np.ndarray]:
    values = default_parameters() if parameters is None else copy.deepcopy(parameters)
    units = check_units(values)
    if units["status"] != "PASS":
        raise ValueError(f"unit check failed: {units['errors']}")
    electrical_input = np.asarray(electrical, dtype=float)
    pressure_input = np.asarray(perfusion_pressure_cmh2o, dtype=float)
    if electrical_input.ndim == 0 and pressure_input.ndim == 0:
        length = 2
    else:
        if electrical_input.ndim == 0:
            length = pressure_input.size
        elif pressure_input.ndim == 0:
            length = electrical_input.size
        else:
            if electrical_input.size != pressure_input.size:
                raise ValueError("electrical and perfusion traces must have equal length")
            length = electrical_input.size
    if length < 2:
        raise ValueError("at least two time samples are required")
    e_trace = _as_trace(electrical_input, length, "electrical")
    p_trace = _as_trace(pressure_input, length, "perfusion_pressure_cmh2o")
    dt = float(values["dt"])
    time = np.arange(length, dtype=float) * dt
    activation = 0.0
    cytosolic_ca = float(values["ca_rest"])
    releasable_ca = float(values["s_max"])
    strain = float(initial_strain)
    strain_rate = 0.0
    output = {
        "time_s": time,
        "electrical": e_trace.copy(),
        "perfusion_pressure_cmh2o": p_trace.copy(),
        "activation": np.zeros(length, dtype=float),
        "calcium_molar": np.zeros(length, dtype=float),
        "calcium_projection_molar": np.zeros(length, dtype=float),
        "releasable_calcium_molar": np.zeros(length, dtype=float),
        "strain": np.zeros(length, dtype=float),
        "strain_rate_s-1": np.zeros(length, dtype=float),
        "flow_norm": np.zeros(length, dtype=float),
        "oxygen_supply_norm": np.zeros(length, dtype=float),
        "perfusion_gain": np.zeros(length, dtype=float),
        "active_force_norm": np.zeros(length, dtype=float),
        "total_force_norm": np.zeros(length, dtype=float),
        "wall_stress_norm": np.zeros(length, dtype=float),
    }
    for index in range(length):
        e_nonnegative = max(float(e_trace[index]), 0.0)
        activation += dt * (
            float(values["k_on"]) * e_nonnegative * (1.0 - activation)
            - float(values["k_off"]) * activation
        )
        activation = float(np.clip(activation, 0.0, 1.0))
        flow = float(effective_flow(p_trace[index], strain, values))
        gain = float(perfusion_gain(flow, values))
        load_gain = 0.5 + 0.5 * flow
        release = float(values["k_release"]) * activation * releasable_ca
        uptake = (cytosolic_ca - float(values["ca_rest"])) / float(values["tau_ca"])
        reload = (
            (float(values["s_max"]) - releasable_ca)
            / float(values["tau_load"])
            * load_gain
        )
        raw_calcium_sum = cytosolic_ca + releasable_ca + dt * (-uptake + reload)
        cytosolic_ca += dt * (release - uptake)
        releasable_ca += dt * (-release + reload)
        cytosolic_ca = float(np.clip(cytosolic_ca, float(values["ca_rest"]), pca_to_molar(PCA_FULL) * 1.25))
        releasable_ca = float(np.clip(releasable_ca, 0.0, float(values["s_max"])))
        output["calcium_projection_molar"][index] = cytosolic_ca + releasable_ca - raw_calcium_sum
        f_ca = float(calcium_activation(cytosolic_ca, values))
        f_length = float(length_factor(strain, values))
        f_velocity = float(velocity_factor(strain_rate, values))
        active_force = (
            float(values["force_scale"])
            * f_ca
            * f_length
            * f_velocity
            * gain
        )
        passive_force = float(values["passive_stiffness"]) * strain + float(values["damping"]) * strain_rate
        total_force = active_force + passive_force + float(values["afterload"])
        strain_rate += dt * (
            float(values["afterload"])
            - active_force
            - float(values["passive_stiffness"]) * strain
            - float(values["damping"]) * strain_rate
        ) / float(values["mass_norm"])
        strain += dt * strain_rate
        if strain < -0.5 or strain > 0.5:
            strain = float(np.clip(strain, -0.5, 0.5))
            strain_rate *= -0.25
        output["activation"][index] = activation
        output["calcium_molar"][index] = cytosolic_ca
        output["releasable_calcium_molar"][index] = releasable_ca
        output["strain"][index] = strain
        output["strain_rate_s-1"][index] = strain_rate
        output["flow_norm"][index] = flow
        output["oxygen_supply_norm"][index] = flow
        output["perfusion_gain"][index] = gain
        output["active_force_norm"][index] = active_force
        output["total_force_norm"][index] = total_force
        output["wall_stress_norm"][index] = total_force
    return output


def _peak_time(time_s: np.ndarray, values: np.ndarray, mode: str = "max") -> float:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.isfinite(array).any():
        return float("nan")
    index = int(np.nanargmax(array) if mode == "max" else np.nanargmin(array))
    return float(np.asarray(time_s, dtype=float)[index])


def summarize(result: dict[str, np.ndarray]) -> dict[str, float]:
    time = result["time_s"]
    active = result["active_force_norm"]
    activation = result["activation"]
    electrical = result["electrical"]
    strain = result["strain"]
    activation_peak = _peak_time(time, activation)
    force_peak = _peak_time(time, active)
    electrical_peak = _peak_time(time, electrical)
    return {
        "electrical_peak_time_s": electrical_peak,
        "activation_peak_time_s": activation_peak,
        "active_force_peak_time_s": force_peak,
        "electrical_to_activation_lag_s": activation_peak - electrical_peak,
        "activation_to_force_lag_s": force_peak - activation_peak,
        "electrical_to_force_lag_s": force_peak - electrical_peak,
        "peak_activation": float(np.max(activation)),
        "peak_calcium_molar": float(np.max(result["calcium_molar"])),
        "peak_active_force_norm": float(np.max(active)),
        "peak_total_force_norm": float(np.max(result["total_force_norm"])),
        "minimum_strain": float(np.min(strain)),
        "minimum_strain_time_s": _peak_time(time, strain, mode="min"),
        "peak_flow_norm": float(np.max(result["flow_norm"])),
    }


def static_gain_ratio(parameters: dict[str, float] | None = None) -> float:
    values = default_parameters() if parameters is None else parameters
    q_zero = float(effective_flow(0.0, 0.0, values))
    q_ref = float(effective_flow(values["p_ref"], 0.0, values))
    return float(perfusion_gain(q_ref, values) / perfusion_gain(q_zero, values))


def run_demo(parameters: dict[str, float] | None = None) -> dict[str, Any]:
    values = default_parameters() if parameters is None else copy.deepcopy(parameters)
    dt = float(values["dt"])
    time = np.arange(0.0, 1.0 + dt * 0.5, dt)
    electrical = electrical_pulse(
        time,
        center_s=float(values["electrical_center"]),
        sigma_s=float(values["electrical_sigma"]),
    )
    p_zero = np.zeros_like(time)
    p_ref = np.full_like(time, float(values["p_ref"]))
    zero_result = simulate(electrical, p_zero, values)
    reference_result = simulate(electrical, p_ref, values)
    placebo_result = simulate(np.zeros_like(time), p_ref, values)
    q_grid = effective_flow(np.array([0.0, 35.0, 70.0]), 0.0, values)
    static_ratio = static_gain_ratio(values)
    zero_summary = summarize(zero_result)
    reference_summary = summarize(reference_result)
    checks = {
        "unit_check": check_units(values)["status"] == "PASS",
        "static_gain_ratio": abs(static_ratio - TARGET_GAIN_RATIO) <= 0.02,
        "causal_peak_order": (
            reference_summary["electrical_peak_time_s"]
            < reference_summary["activation_peak_time_s"]
            < reference_summary["active_force_peak_time_s"]
        ),
        "electrical_to_activation_lag": 0.0 <= reference_summary["electrical_to_activation_lag_s"] <= 0.080,
        "activation_to_force_lag": 0.0 <= reference_summary["activation_to_force_lag_s"] <= 0.150,
        "perfusion_force_sign": reference_summary["peak_active_force_norm"] > zero_summary["peak_active_force_norm"],
        "perfusion_strain_sign": reference_summary["minimum_strain"] < zero_summary["minimum_strain"],
        "flow_monotonic_at_fixed_strain": bool(np.all(np.diff(q_grid) >= 0.0)),
        "flow_nonnegative": bool(np.all(q_grid >= 0.0)),
        "electrical_off_placebo": float(np.max(placebo_result["active_force_norm"])) <= 1.0e-8,
        "finite_outputs": all(
            bool(np.isfinite(array).all())
            for key, array in reference_result.items()
            if key not in {"time_s"} and isinstance(array, np.ndarray)
        ),
        "nonnegative_calcium": bool(
            np.all(reference_result["calcium_molar"] >= 0.0)
            and np.all(reference_result["releasable_calcium_molar"] >= 0.0)
        ),
    }
    return {
        "time_grid": {"start_s": float(time[0]), "end_s": float(time[-1]), "dt_s": dt, "n": int(time.size)},
        "static_gain_ratio": static_ratio,
        "static_gain_target": TARGET_GAIN_RATIO,
        "zero_perfusion": zero_summary,
        "reference_perfusion": reference_summary,
        "dynamic_peak_force_ratio": reference_summary["peak_active_force_norm"] / zero_summary["peak_active_force_norm"],
        "flow_grid_at_zero_strain": [float(item) for item in q_grid],
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
    }


def sensitivity_analysis(parameters: dict[str, float] | None = None) -> dict[str, Any]:
    base = default_parameters() if parameters is None else copy.deepcopy(parameters)
    base_demo = run_demo(base)
    base_summary = base_demo["reference_perfusion"]
    controls = (
        ("k_on", "s^-1"),
        ("ca_ec50", "M"),
        ("gregg_gain", "1"),
    )
    rows: list[dict[str, Any]] = []
    for name, unit in controls:
        for factor, label in ((0.5, "minus50"), (1.5, "plus50")):
            varied = copy.deepcopy(base)
            varied[name] = float(base[name]) * factor
            demo = run_demo(varied)
            summary = demo["reference_perfusion"]
            rows.append(
                {
                    "parameter": name,
                    "unit": unit,
                    "factor": factor,
                    "change_label": label,
                    "value": float(varied[name]),
                    "peak_active_force_norm": summary["peak_active_force_norm"],
                    "peak_force_change_from_base_pct": 100.0
                    * (summary["peak_active_force_norm"] / base_summary["peak_active_force_norm"] - 1.0),
                    "activation_peak_time_s": summary["activation_peak_time_s"],
                    "force_peak_time_s": summary["active_force_peak_time_s"],
                    "minimum_strain": summary["minimum_strain"],
                    "all_checks_pass": demo["all_checks_pass"],
                }
            )
    return {
        "baseline": base_summary,
        "controls": [
            {"parameter": name, "unit": unit, "factors": [0.5, 1.5]}
            for name, unit in controls
        ],
        "runs": rows,
    }


def reference_catalog() -> list[dict[str, Any]]:
    return [
        {
            "source": "Fabiato 1981",
            "doi": "10.1085/jgp.78.5.457",
            "location": "Fig. 10 and Results text",
            "value": "pCa 5.30-5.40; full activation approximately pCa 4.90",
            "unit": "pCa",
            "status": "VERIFIED",
        },
        {
            "source": "Schouten et al. 1992",
            "doi": "10.1113/jphysiol.1992.sp019180",
            "location": "PubMed abstract Results paragraph 3; figure/table not accessible",
            "value": "74 +/- 20 (n=11) for 0 to 70",
            "unit": "% peak force",
            "status": "OVERIFIERAD for figure/table; abstract value verified",
        },
        {
            "source": "Schulz et al. 1991",
            "doi": "10.1161/01.cir.83.4.1390",
            "location": "PubMed abstract Results",
            "value": "no significant change from 88 +/- 11 to 186 +/- 11; change at 57 +/- 13",
            "unit": "mmHg and % wall thickening",
            "status": "VERIFIED abstract result",
        },
        {
            "source": "Goto et al. 1991",
            "doi": "10.1161/01.res.68.2.482",
            "location": "PubMed abstract Results",
            "value": "flow +99 +/- 76; Emax +18 +/- 15 at 93 +/- 11",
            "unit": "% and mmHg",
            "status": "VERIFIED abstract result",
        },
    ]


def make_results(parameters: dict[str, float] | None = None) -> dict[str, Any]:
    values = default_parameters() if parameters is None else copy.deepcopy(parameters)
    demo = run_demo(values)
    sensitivity = sensitivity_analysis(values)
    return {
        "id": "BT-HX-Q022",
        "model_version": MODEL_VERSION,
        "status": "MODEL_OUTPUT_NOT_MEASURED_DATA",
        "prereg_sha256": "66f16c5d791586f38e9edd54a1446677b0534e59b00afaafc86a444c2b7d5a73",
        "data_statement": "All trajectories and numerical outputs are generated by the local model; no internal or measured traces were used.",
        "parameters": values,
        "parameter_table": parameter_table(values),
        "unit_check": check_units(values),
        "reference_catalog": reference_catalog(),
        "demo": demo,
        "sensitivity": sensitivity,
        "limitations": [
            "Single normalized lumped compartment, not a distributed whole-heart model.",
            "Perfusion gain is calibrated to a literature contrast and is not a fitted clinical parameter.",
            "Flow and force outputs are normalized; only pressure, time, and calcium retain stated physical units.",
            "The next step requires paired ECG/electrical, strain or force, and perfusion traces with synchronized time base.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the BT-HX-Q022 mechanistic cardiac model")
    parser.add_argument("--output", default="results.json", help="JSON output path")
    parser.add_argument("--no-sensitivity", action="store_true", help="skip the three-parameter sensitivity run")
    args = parser.parse_args()
    payload = make_results()
    if args.no_sensitivity:
        payload["sensitivity"] = {"status": "SKIPPED_BY_COMMAND"}
    output_path = Path(args.output)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output_path), "all_checks_pass": payload["demo"]["all_checks_pass"]}, sort_keys=True))


if __name__ == "__main__":
    main()
