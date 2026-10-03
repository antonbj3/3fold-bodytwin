from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    tau_r: float = 5.0
    tau_a: float = 20.0
    tau_f: float = 2.2
    tau_total: float = 49.0
    delay_ra: float = 2.0
    delay_ac: float = 3.0
    delay_fa: float = 8.0
    feedback_gain: float = 0.9
    feedback_scale: float = 0.5
    hill_exponent: float = 2.0
    acth_ec50: float = 0.35
    max_secretion: float = 1.0
    pulse_amplitude: float = 1.0
    pulse_duration: float = 10.0
    t_end: float = 180.0
    dt: float = 0.1
    initial_state: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0)


@dataclass(frozen=True)
class DelaySpec:
    name: str
    length_min: float
    stages: int
    stage_tau_min: float


@dataclass
class Trajectory:
    time_min: np.ndarray
    values: np.ndarray
    success: bool
    message: str

    @property
    def drive(self) -> np.ndarray:
        return self.values[:, 0]

    @property
    def acth(self) -> np.ndarray:
        return self.values[:, 1]

    @property
    def free_cortisol(self) -> np.ndarray:
        return self.values[:, 2]

    @property
    def total_cortisol(self) -> np.ndarray:
        return self.values[:, 3]

    def column(self, name: str) -> np.ndarray:
        indices = {
            "drive": 0,
            "acth": 1,
            "free_cortisol": 2,
            "total_cortisol": 3,
        }
        return self.values[:, indices[name]]


def default_parameters() -> Parameters:
    return Parameters()


def input_drive(time_min: float, p: Parameters) -> float:
    if 0.0 <= time_min < p.pulse_duration:
        return p.pulse_amplitude
    return 0.0


def secretion_from_acth(acth: float, p: Parameters) -> float:
    positive = max(acth, 0.0)
    numerator = positive ** p.hill_exponent
    denominator = p.acth_ec50 ** p.hill_exponent + numerator
    return p.max_secretion * numerator / denominator


def feedback_from_cortisol(free_cortisol: float, p: Parameters) -> float:
    positive = max(free_cortisol, 0.0)
    return p.feedback_gain * positive / (p.feedback_scale + positive)


def _delay_spec(name: str, length_min: float) -> DelaySpec | None:
    if length_min == 0.0:
        return None
    stages = max(1, int(math.ceil(length_min / 0.5)))
    return DelaySpec(name, length_min, stages, length_min / stages)


def _delay_specs(p: Parameters) -> tuple[DelaySpec | None, DelaySpec | None, DelaySpec | None]:
    return (
        _delay_spec("drive_to_acth", p.delay_ra),
        _delay_spec("acth_to_cortisol", p.delay_ac),
        _delay_spec("cortisol_to_acth_feedback", p.delay_fa),
    )


def unit_check(p: Parameters) -> dict[str, Any]:
    numeric_fields = [
        "tau_r",
        "tau_a",
        "tau_f",
        "tau_total",
        "delay_ra",
        "delay_ac",
        "delay_fa",
        "feedback_gain",
        "feedback_scale",
        "hill_exponent",
        "acth_ec50",
        "max_secretion",
        "pulse_amplitude",
        "pulse_duration",
        "t_end",
        "dt",
    ]
    numeric_values = {name: float(getattr(p, name)) for name in numeric_fields}
    positive_fields = [
        "tau_r",
        "tau_a",
        "tau_f",
        "tau_total",
        "feedback_scale",
        "hill_exponent",
        "acth_ec50",
        "max_secretion",
        "pulse_amplitude",
        "t_end",
        "dt",
    ]
    nonnegative_fields = ["delay_ra", "delay_ac", "delay_fa", "feedback_gain", "pulse_duration"]
    state_values = np.asarray(p.initial_state, dtype=float)
    checks = {
        "finite_parameters": bool(np.all(np.isfinite(list(numeric_values.values())))),
        "positive_time_constants": bool(all(numeric_values[name] > 0.0 for name in positive_fields)),
        "nonnegative_delays_and_gain": bool(all(numeric_values[name] >= 0.0 for name in nonnegative_fields)),
        "finite_initial_state": bool(np.all(np.isfinite(state_values))),
        "initial_state_dimension": state_values.shape == (4,),
        "time_minutes": True,
        "state_unit": "relative concentration or relative drive",
        "derivative_unit": "relative concentration per minute",
        "delay_equation": "d q0/dt=(u-q0)/tau; d qi/dt=(q(i-1)-qi)/tau",
        "secretion_equation": "S(a)=Smax*a^n/(EC50^n+a^n), n>0, a>=0",
        "feedback_equation": "F(f)=g*f/(K+f), f>=0",
    }
    checks["status"] = "pass" if all(
        checks[name] for name in (
            "finite_parameters",
            "positive_time_constants",
            "nonnegative_delays_and_gain",
            "finite_initial_state",
            "initial_state_dimension",
        )
    ) else "fail"
    if checks["status"] != "pass":
        raise ValueError(checks)
    return checks


def parameter_table(p: Parameters | None = None) -> list[dict[str, Any]]:
    p = default_parameters() if p is None else p
    values = asdict(p)
    provenance = {
        "tau_r": ("min", "assumption; upstream hypothalamic drive relaxation"),
        "tau_a": ("min", "source-anchored ACTH half-life 20 min, Keenan et al. 2004 DOI 10.1152/ajpendo.00167.2004; tau=half-life/ln(2)"),
        "tau_f": ("min", "source-anchored free-cortisol half-life 2.2 min, Dorin et al. 2012 DOI 10.1210/jc.2011-2227 Fig. 5; tau=half-life/ln(2)"),
        "tau_total": ("min", "source-anchored total-cortisol half-life 49 min, Keenan et al. 2004 DOI 10.1152/ajpendo.00167.2004; tau=half-life/ln(2)"),
        "delay_ra": ("min", "assumption; distributed upstream drive-to-ACTH delay"),
        "delay_ac": ("min", "source-anchored ACTH-to-cortisol response lag 3.0 min, Dorin et al. 2012 DOI 10.1210/jc.2011-2227 Results/Fig. 2; represented as a distributed line"),
        "delay_fa": ("min", "assumption; distributed cortisol-to-ACTH negative-feedback delay"),
        "feedback_gain": ("relative drive", "assumption; dimensionless feedback strength"),
        "feedback_scale": ("relative cortisol", "assumption; half-saturation scale"),
        "hill_exponent": ("dimensionless", "assumption; positive secretion/feedback shape exponent"),
        "acth_ec50": ("relative ACTH", "assumption; normalized secretion half-saturation"),
        "max_secretion": ("relative cortisol", "assumption; normalized maximal secretion target"),
        "pulse_amplitude": ("relative drive", "frozen intervention amplitude"),
        "pulse_duration": ("min", "frozen intervention duration"),
        "t_end": ("min", "frozen observation horizon"),
        "dt": ("min", "numerical output and maximum solver step"),
        "initial_state": ("mixed relative state", "frozen zero-deviation initial condition"),
    }
    return [
        {"parameter": name, "value": values[name], "unit": unit, "provenance": text}
        for name, (unit, text) in provenance.items()
    ]


def simulate(p: Parameters | None = None, t_end: float | None = None, dt: float | None = None) -> Trajectory:
    p = default_parameters() if p is None else p
    unit_check(p)
    horizon = p.t_end if t_end is None else float(t_end)
    step = p.dt if dt is None else float(dt)
    if horizon <= 0.0 or step <= 0.0:
        raise ValueError("t_end and dt must be positive")
    t_eval = np.arange(0.0, horizon + 0.5 * step, step, dtype=float)
    t_eval = t_eval[t_eval <= horizon + 1.0e-12]
    specs = _delay_specs(p)
    base_size = 4
    offsets: list[int | None] = []
    cursor = base_size
    for spec in specs:
        if spec is None:
            offsets.append(None)
        else:
            offsets.append(cursor)
            cursor += spec.stages
    state_size = cursor
    y0 = np.zeros(state_size, dtype=float)
    y0[:4] = np.asarray(p.initial_state, dtype=float)

    def rhs(time_min: float, y: np.ndarray) -> np.ndarray:
        drive, acth, free_cortisol, _ = y[:4]
        if offsets[0] is None:
            drive_delayed = drive
        else:
            start = offsets[0]
            spec = specs[0]
            drive_delayed = y[start + spec.stages - 1]
        if offsets[1] is None:
            acth_delayed = acth
        else:
            start = offsets[1]
            spec = specs[1]
            acth_delayed = y[start + spec.stages - 1]
        if offsets[2] is None:
            feedback_delayed = free_cortisol
        else:
            start = offsets[2]
            spec = specs[2]
            feedback_delayed = y[start + spec.stages - 1]
        derivative = np.zeros_like(y)
        derivative[0] = (input_drive(time_min, p) - drive) / p.tau_r
        derivative[1] = (
            (drive_delayed - acth) / p.tau_a
            - feedback_from_cortisol(feedback_delayed, p) / p.tau_a
        )
        derivative[2] = (secretion_from_acth(acth_delayed, p) - free_cortisol) / p.tau_f
        derivative[3] = (free_cortisol - y[3]) / p.tau_total
        for index, spec in enumerate(specs):
            if spec is None or offsets[index] is None:
                continue
            start = offsets[index]
            source = (drive, acth, free_cortisol)[index]
            derivative[start] = (source - y[start]) / spec.stage_tau_min
            for stage in range(1, spec.stages):
                derivative[start + stage] = (y[start + stage - 1] - y[start + stage]) / spec.stage_tau_min
        return derivative

    solution = solve_ivp(
        rhs,
        (0.0, horizon),
        y0,
        method="DOP853",
        t_eval=t_eval,
        max_step=step,
        rtol=1.0e-8,
        atol=1.0e-10,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return Trajectory(solution.t.copy(), solution.y[:4, :].T.copy(), solution.success, solution.message)


def _onset(time_min: np.ndarray, signal: np.ndarray, fraction: float = 0.1) -> float | None:
    peak = float(np.max(signal)) if signal.size else 0.0
    if peak <= 0.0:
        return None
    indices = np.flatnonzero(signal >= fraction * peak)
    return float(time_min[indices[0]]) if indices.size else None


def _peak_time(time_min: np.ndarray, signal: np.ndarray) -> float:
    if signal.size == 0:
        return float("nan")
    return float(time_min[int(np.argmax(signal))])


def metrics(trajectory: Trajectory) -> dict[str, Any]:
    time_min = trajectory.time_min
    acth = trajectory.acth
    free_cortisol = trajectory.free_cortisol
    total_cortisol = trajectory.total_cortisol
    acth_onset = _onset(time_min, acth)
    free_onset = _onset(time_min, free_cortisol)
    free_peak_time = _peak_time(time_min, free_cortisol)
    after_free_peak = time_min >= free_peak_time
    acth_after = acth[after_free_peak]
    minimum_acth_after = float(np.min(acth_after)) if acth_after.size else float("nan")
    final_acth = float(acth[-1]) if acth.size else float("nan")
    phase_lag = None if acth_onset is None or free_onset is None else free_onset - acth_onset
    free_peak = float(np.max(free_cortisol)) if free_cortisol.size else 0.0
    total_peak = float(np.max(total_cortisol)) if total_cortisol.size else 0.0
    acth_peak_time = _peak_time(time_min, acth)
    total_peak_time = _peak_time(time_min, total_cortisol)
    finite = bool(np.all(np.isfinite(trajectory.values)))
    return {
        "acth_peak_time_min": acth_peak_time,
        "free_cortisol_peak_time_min": free_peak_time,
        "total_cortisol_peak_time_min": total_peak_time,
        "acth_to_free_peak_phase_lag_min": free_peak_time - acth_peak_time,
        "acth_onset_time_min": acth_onset,
        "free_cortisol_onset_time_min": free_onset,
        "acth_to_free_phase_lag_min": phase_lag,
        "acth_peak": float(np.max(acth)) if acth.size else 0.0,
        "free_cortisol_peak": free_peak,
        "total_cortisol_peak": total_peak,
        "minimum_acth_after_free_peak": minimum_acth_after,
        "final_acth": final_acth,
        "feedback_suppression_contrast": max(0.0, -minimum_acth_after),
        "rebound_rise_after_suppression": final_acth - minimum_acth_after,
        "free_to_total_peak_lag_min": total_peak_time - free_peak_time,
        "free_peak_to_total_peak_order_ok": total_peak_time > free_peak_time,
        "stable_finite": finite,
        "phase_criterion_1_to_5_min": bool(phase_lag is not None and 1.0 <= phase_lag <= 5.0),
        "transient_criterion_feedback_signal": bool(minimum_acth_after < -0.01 and final_acth > minimum_acth_after),
    }


def _metric_json(value: dict[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value_item in value.items():
        if isinstance(value_item, (np.floating, np.integer)):
            result[key] = value_item.item()
        elif value_item is None:
            result[key] = None
        else:
            result[key] = value_item
    return result


def sensitivity(p: Parameters, names: tuple[str, ...], factors: tuple[float, ...] = (0.5, 1.0, 1.5)) -> dict[str, dict[str, dict[str, Any]]]:
    output: dict[str, dict[str, dict[str, Any]]] = {}
    for name in names:
        output[name] = {}
        base_value = float(getattr(p, name))
        for factor in factors:
            changed = replace(p, **{name: base_value * factor})
            output[name][str(factor)] = _metric_json(metrics(simulate(changed)))
    return output


def parameter_scan(p: Parameters, name: str, values: tuple[float, ...]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for value in values:
        changed = replace(p, **{name: value})
        row = {"parameter": name, "value": value}
        row.update(_metric_json(metrics(simulate(changed))))
        rows.append(row)
    return rows


def delay_scan(p: Parameters, values: tuple[float, ...] = (0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0)) -> list[dict[str, Any]]:
    return parameter_scan(p, "delay_fa", values)


def prereg_hash() -> str:
    path = Path(__file__).with_name("PREREG.sha256")
    if not path.exists():
        raise FileNotFoundError("PREREG.sha256 must exist before analysis")
    token = path.read_text(encoding="utf-8").split()[0]
    return token


def run_analysis(p: Parameters | None = None) -> dict[str, Any]:
    p = default_parameters() if p is None else p
    default_trajectory = simulate(p)
    no_feedback_parameters = replace(p, feedback_gain=0.0)
    no_feedback_trajectory = simulate(no_feedback_parameters)
    default_metrics = metrics(default_trajectory)
    no_feedback_metrics = metrics(no_feedback_trajectory)
    controlling = sensitivity(p, ("delay_ac", "delay_fa", "feedback_gain"))
    feedback_contrast = {
        "default": default_metrics["feedback_suppression_contrast"],
        "no_feedback": no_feedback_metrics["feedback_suppression_contrast"],
        "increase_from_no_feedback": default_metrics["feedback_suppression_contrast"] - no_feedback_metrics["feedback_suppression_contrast"],
    }
    return {
        "id": "BT-HX-Q023",
        "axis": "human HPA: hypothalamic drive -> ACTH -> free cortisol -> total cortisol; cortisol negative feedback to ACTH",
        "equation_units": {
            "time": "min",
            "state": "relative concentration or relative drive",
            "rhs": "relative concentration/min",
            "input": "relative drive",
            "delay": "min",
        },
        "equations": {
            "drive": "dr/dt=(u(t)-r)/tau_r",
            "acth": "da/dt=(r_delay(t)-a)/tau_a-(g/tau_a)*f_delay(t)/(K_f+f_delay(t))",
            "free_cortisol": "df/dt=(S(a_delay(t))-f)/tau_f",
            "total_cortisol": "dc/dt=(f-c)/tau_total",
            "secretion": "S(a)=Smax*max(a,0)^n/(EC50^n+max(a,0)^n)",
            "delay_line": "dq0/dt=(u-q0)/delta, dqi/dt=(q(i-1)-qi)/delta",
        },
        "parameter_table": parameter_table(p),
        "unit_check": unit_check(p),
        "default": _metric_json(default_metrics),
        "no_feedback_counterfactual": _metric_json(no_feedback_metrics),
        "feedback_contrast": feedback_contrast,
        "sensitivity_plus_minus_50_percent": controlling,
        "delay_ac_scan": parameter_scan(p, "delay_ac", (0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)),
        "delay_fa_scan": delay_scan(p),
        "frozen_criteria": {
            "phase_order_positive": bool(default_metrics["acth_to_free_phase_lag_min"] is not None and default_metrics["acth_to_free_phase_lag_min"] > 0.0),
            "phase_lag_1_to_5_min": default_metrics["phase_criterion_1_to_5_min"],
            "free_before_total_peak": default_metrics["free_peak_to_total_peak_order_ok"],
            "feedback_transient": default_metrics["transient_criterion_feedback_signal"],
            "unit_check": unit_check(p)["status"],
        },
        "source_anchors": {
            "acth_to_cortisol_response_lag_min": {"mean": 3.0, "sd": 0.92, "source": "Dorin et al. 2012, Results with Fig. 2"},
            "free_cortisol_half_life_min": {"mean": 2.2, "sd": 1.3, "source": "Dorin et al. 2012, Fig. 5 and Results"},
            "maximal_cortisol_secretion_rate_nmol_l_s": {"mean": 0.44, "sd": 0.13, "source": "Dorin et al. 2012, abstract and Results"},
            "total_cortisol_half_life_min": {"mean": 49.0, "sd": 2.4, "source": "Keenan et al. 2004, primary 24 h cohort"},
            "acth_half_life_min": {"mean": 20.0, "sd": 1.3, "source": "Keenan et al. 2004, primary 24 h cohort"},
            "hcrf_elimination_half_life_10_min_iv_min": {"mean": 45.0, "sd": 7.0, "source": "Angst et al. 1998, abstract"},
        },
        "prereg_sha256": prereg_hash(),
        "provenance": {
            "data_status": "No internal hormone measurements were supplied; all trajectory values are predictions.",
            "unknowns": [
                "The internal observed feedback trace",
                "Individual CRH, ACTH, free-cortisol, and total-cortisol trajectories",
                "The true portal and feedback transport delays",
            ],
            "source_1": "Dorin et al. 2012 DOI 10.1210/jc.2011-2227, Fig. 2/5 and Results",
            "source_2": "Keenan et al. 2004 DOI 10.1152/ajpendo.00167.2004, primary 24 h sampling cohort",
            "source_3": "Angst et al. 1998 DOI 10.1016/S0009-9236(98)90133-3, hCRF PK intervention",
            "interpretation": "Mechanistic sensitivity analysis, not a fit to an observed internal trace.",
        },
        "next_resolution_step": "Fit the same equations to timestamped simultaneous CRH, ACTH, free/total cortisol measurements with intervention metadata and assay compartment labels.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results.json")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    result = run_analysis()
    encoded = json.dumps(result, indent=2, sort_keys=True)
    if not args.no_write:
        Path(args.output).write_text(encoded + "\n", encoding="utf-8")
    print(encoded)


if __name__ == "__main__":
    main()
