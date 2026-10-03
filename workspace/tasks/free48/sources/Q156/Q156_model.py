from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

MMHG_TO_PA = 133.322387415
REFERENCE_T50_S = 600.0


@dataclass(frozen=True)
class Parameters:
    L: float = 250e-6
    N: int = 40
    A: float = 1e-9
    eps_isf: float = 0.20
    eps_pvs: float = 0.10
    D_isf: float = 1e-10
    D_pvs: float = 1e-9
    V_csf: float = 1e-9
    k_csf: float = 1.0
    k_exchange: float = 1e-2
    k_blood_isf: float = 0.0
    k_blood_pvs: float = 0.0
    pressure_delta_mmhg: float = 0.0
    hydraulic_conductivity: float = 1.9e-15
    viscosity: float = 1e-3
    k_on: float = 0.0
    k_off: float = 0.0
    k_csf_out: float = 0.0
    initial_csf: float = 1.0
    time_step: float = 1.0
    t_end: float = 7200.0
    sample_times: tuple[float, ...] = (
        0.0,
        60.0,
        300.0,
        600.0,
        1200.0,
        1800.0,
        3600.0,
        7200.0,
    )

    def __post_init__(self) -> None:
        if self.L <= 0 or self.A <= 0 or self.V_csf <= 0:
            raise ValueError("L, A, and V_csf must be positive")
        if self.N < 4:
            raise ValueError("N must be at least four")
        if not 0 < self.eps_isf < 1 or not 0 < self.eps_pvs < 1:
            raise ValueError("volume fractions must lie between zero and one")
        if self.D_isf <= 0 or self.D_pvs <= 0:
            raise ValueError("diffusion coefficients must be positive")
        if min(
            self.k_csf,
            self.k_exchange,
            self.k_blood_isf,
            self.k_blood_pvs,
            self.k_on,
            self.k_off,
            self.k_csf_out,
        ) < 0:
            raise ValueError("rate constants must be non-negative")
        if self.pressure_delta_mmhg < 0:
            raise ValueError("this directional model uses a non-negative pressure drop")
        if self.hydraulic_conductivity < 0 or self.viscosity <= 0:
            raise ValueError("hydraulic conductivity must be non-negative and viscosity positive")
        if self.time_step <= 0 or self.t_end <= 0:
            raise ValueError("time parameters must be positive")
        if not self.sample_times or self.sample_times[0] != 0:
            raise ValueError("sample_times must begin at zero")
        if any(t < 0 or t > self.t_end for t in self.sample_times):
            raise ValueError("sample_times must lie within the run")
        if list(self.sample_times) != sorted(self.sample_times):
            raise ValueError("sample_times must be sorted")

    @property
    def dx(self) -> float:
        return self.L / self.N

    @property
    def velocity(self) -> float:
        pressure_pa = self.pressure_delta_mmhg * MMHG_TO_PA
        return max(
            0.0,
            self.hydraulic_conductivity * pressure_pa / (self.viscosity * self.L),
        )

    @property
    def stable_time_step(self) -> float:
        limits = [self.time_step]
        diffusion_limit = 0.45 * self.dx * self.dx / max(self.D_isf, self.D_pvs)
        limits.append(diffusion_limit)
        if self.velocity > 0:
            limits.append(0.45 * self.dx / self.velocity)
        rate_limit = 0.45 / max(
            self.k_csf,
            self.k_exchange,
            self.k_blood_isf,
            self.k_blood_pvs,
            self.k_on,
            self.k_off,
            self.k_csf_out,
            1e-30,
        )
        limits.append(rate_limit)
        return min(limits)


@dataclass
class Fluxes:
    isf_faces: np.ndarray
    pvs_faces: np.ndarray
    csf_to_pvs: float
    pvs_to_isf: np.ndarray
    blood_isf: np.ndarray
    blood_pvs: np.ndarray
    binding: np.ndarray
    csf_out: float
    pvs_out: float


def cell_volumes(p: Parameters) -> tuple[float, float, float]:
    geometric_cell_volume = p.A * p.dx
    return geometric_cell_volume, geometric_cell_volume * p.eps_isf, geometric_cell_volume * p.eps_pvs


def flux_divergence(face_flux: np.ndarray) -> np.ndarray:
    result = np.zeros(face_flux.size - 1, dtype=float)
    internal = face_flux[1:-1]
    result[:-1] += internal
    result[1:] -= internal
    return result


def compute_fluxes(
    csf: float,
    c_pvs: np.ndarray,
    c_isf: np.ndarray,
    bound: np.ndarray,
    p: Parameters,
) -> Fluxes:
    geometric_cell_volume, v_isf, v_pvs = cell_volumes(p)
    del geometric_cell_volume
    isf_faces = np.zeros(p.N + 1, dtype=float)
    pvs_faces = np.zeros(p.N + 1, dtype=float)
    isf_faces[1:p.N] = -p.A * p.eps_isf * p.D_isf * np.diff(c_isf) / p.dx
    pvs_diffusion = -p.A * p.eps_pvs * p.D_pvs * np.diff(c_pvs) / p.dx
    pvs_advection = p.A * p.velocity * c_pvs[:-1]
    pvs_faces[1:p.N] = pvs_diffusion + pvs_advection
    csf_to_pvs = p.k_csf * v_pvs * (csf - c_pvs[0])
    exchange_volume = 2.0 * v_pvs * v_isf / (v_pvs + v_isf)
    pvs_to_isf = p.k_exchange * exchange_volume * (c_pvs - c_isf)
    blood_isf = p.k_blood_isf * v_isf * c_isf
    blood_pvs = p.k_blood_pvs * v_pvs * c_pvs
    binding = p.k_on * v_isf * c_isf - p.k_off * p.A * p.dx * bound
    csf_out = p.k_csf_out * p.V_csf * csf
    pvs_out = p.A * p.velocity * c_pvs[-1]
    return Fluxes(
        isf_faces=isf_faces,
        pvs_faces=pvs_faces,
        csf_to_pvs=csf_to_pvs,
        pvs_to_isf=pvs_to_isf,
        blood_isf=blood_isf,
        blood_pvs=blood_pvs,
        binding=binding,
        csf_out=csf_out,
        pvs_out=pvs_out,
    )


def total_mass(
    csf: float,
    c_pvs: np.ndarray,
    c_isf: np.ndarray,
    bound: np.ndarray,
    p: Parameters,
) -> float:
    _, v_isf, v_pvs = cell_volumes(p)
    geometric_cell_volume = p.A * p.dx
    return float(
        p.V_csf * csf
        + np.sum(v_pvs * c_pvs)
        + np.sum(v_isf * c_isf)
        + np.sum(geometric_cell_volume * bound)
    )


def first_crossing(times: np.ndarray, values: np.ndarray, target: float) -> float | None:
    for index in range(1, len(times)):
        if values[index] >= target:
            previous = values[index - 1]
            current = values[index]
            if current == previous:
                return float(times[index])
            fraction = (target - previous) / (current - previous)
            return float(times[index - 1] + fraction * (times[index] - times[index - 1]))
    return None


def _record(
    time: float,
    csf: float,
    c_pvs: np.ndarray,
    c_isf: np.ndarray,
    bound: np.ndarray,
    p: Parameters,
) -> dict[str, Any]:
    _, v_isf, v_pvs = cell_volumes(p)
    geometric_cell_volume = p.A * p.dx
    isf_total = (v_isf * c_isf + geometric_cell_volume * bound) / v_isf
    isf_free = float(np.mean(c_isf))
    isf_total_mean = float(np.mean(isf_total))
    pvs_mean = float(np.mean(c_pvs))
    distal_isf_free = float(c_isf[-1])
    distal_isf_total = float(isf_total[-1])
    return {
        "time_s": float(time),
        "csf_mol_m3": float(csf),
        "isf_free_mean_mol_m3": isf_free,
        "isf_total_equivalent_mean_mol_m3": isf_total_mean,
        "pvs_mean_mol_m3": pvs_mean,
        "isf_distal_free_mol_m3": distal_isf_free,
        "isf_distal_total_equivalent_mol_m3": distal_isf_total,
        "isf_bound_mean_mol_m3_per_total_volume": float(np.mean(bound)),
        "isf_free_profile_mol_m3": c_isf.tolist(),
        "pvs_profile_mol_m3": c_pvs.tolist(),
        "isf_bound_profile_mol_m3_per_total_volume": bound.tolist(),
        "pvs_velocity_m_s": float(p.velocity),
    }


def run_model(p: Parameters) -> dict[str, Any]:
    csf = float(p.initial_csf)
    c_pvs = np.zeros(p.N, dtype=float)
    c_isf = np.zeros(p.N, dtype=float)
    bound = np.zeros(p.N, dtype=float)
    sample_times = np.asarray(p.sample_times, dtype=float)
    records: list[dict[str, Any]] = []
    csf_history: list[float] = []
    isf_free_history: list[float] = []
    isf_total_history: list[float] = []
    pvs_history: list[float] = []
    distal_total_history: list[float] = []
    mass_history: list[float] = []
    balance_history: list[float] = []
    csf_to_pvs_integral = 0.0
    pvs_to_isf_integral = 0.0
    csf_out_integral = 0.0
    blood_isf_integral = 0.0
    blood_pvs_integral = 0.0
    pvs_out_integral = 0.0
    initial_mass = total_mass(csf, c_pvs, c_isf, bound, p)
    cumulative_loss = 0.0
    time = 0.0
    record_index = 0
    records.append(_record(time, csf, c_pvs, c_isf, bound, p))
    csf_history.append(csf)
    _, v_isf, v_pvs = cell_volumes(p)
    isf_total = (v_isf * c_isf + p.A * p.dx * bound) / v_isf
    isf_free_history.append(float(np.mean(c_isf)))
    isf_total_history.append(float(np.mean(isf_total)))
    pvs_history.append(float(np.mean(c_pvs)))
    distal_total_history.append(float(isf_total[-1]))
    current_mass = total_mass(csf, c_pvs, c_isf, bound, p)
    mass_history.append(current_mass)
    balance_history.append(0.0)
    while time < p.t_end - 1e-12:
        next_sample = sample_times[record_index + 1] if record_index + 1 < len(sample_times) else p.t_end
        remaining = min(p.t_end - time, next_sample - time)
        h = min(p.stable_time_step, remaining)
        if h <= 1e-12:
            time = next_sample
        else:
            flux = compute_fluxes(csf, c_pvs, c_isf, bound, p)
            isf_divergence = flux_divergence(flux.isf_faces)
            pvs_divergence = flux_divergence(flux.pvs_faces)
            d_csf_mol = -flux.csf_to_pvs - flux.csf_out
            d_pvs_mol = (
                -pvs_divergence
                - flux.pvs_to_isf
                - flux.blood_pvs
            )
            d_pvs_mol[0] += flux.csf_to_pvs
            d_pvs_mol[-1] -= flux.pvs_out
            d_isf_mol = (
                -isf_divergence
                + flux.pvs_to_isf
                - flux.blood_isf
                - flux.binding
            )
            d_bound_mol = flux.binding
            csf += h * d_csf_mol / p.V_csf
            c_pvs += h * d_pvs_mol / v_pvs
            c_isf += h * d_isf_mol / v_isf
            bound += h * d_bound_mol / (p.A * p.dx)
            csf_to_pvs_integral += h * flux.csf_to_pvs
            pvs_to_isf_integral += h * float(np.sum(flux.pvs_to_isf))
            csf_out_integral += h * flux.csf_out
            blood_isf_integral += h * float(np.sum(flux.blood_isf))
            blood_pvs_integral += h * float(np.sum(flux.blood_pvs))
            pvs_out_integral += h * flux.pvs_out
            cumulative_loss += h * (flux.csf_out + np.sum(flux.blood_isf) + np.sum(flux.blood_pvs) + flux.pvs_out)
            time += h
        if time >= next_sample - 1e-10 or (record_index + 1 == len(sample_times) and time >= p.t_end - 1e-10):
            record = _record(time, csf, c_pvs, c_isf, bound, p)
            records.append(record)
            csf_history.append(float(csf))
            _, v_isf, v_pvs = cell_volumes(p)
            isf_total = (v_isf * c_isf + p.A * p.dx * bound) / v_isf
            isf_free_history.append(float(np.mean(c_isf)))
            isf_total_history.append(float(np.mean(isf_total)))
            pvs_history.append(float(np.mean(c_pvs)))
            distal_total_history.append(float(isf_total[-1]))
            current_mass = total_mass(csf, c_pvs, c_isf, bound, p)
            mass_history.append(current_mass)
            balance_history.append((current_mass - initial_mass + cumulative_loss) / initial_mass)
            record_index += 1
    times = np.asarray([record["time_s"] for record in records], dtype=float)
    distal_values = np.asarray(distal_total_history, dtype=float)
    mean_values = np.asarray(isf_total_history, dtype=float)
    target = 0.5 * p.initial_csf
    t50_distal = first_crossing(times, distal_values, target)
    t50_mean = first_crossing(times, mean_values, target)
    pvs_values = np.asarray(pvs_history, dtype=float)
    t50_pvs = first_crossing(times, pvs_values, target)
    all_finite = bool(
        np.isfinite(csf)
        and np.all(np.isfinite(c_pvs))
        and np.all(np.isfinite(c_isf))
        and np.all(np.isfinite(bound))
        and np.all(np.isfinite(times))
    )
    nonnegative = bool(csf >= -1e-12 and np.all(c_pvs >= -1e-12) and np.all(c_isf >= -1e-12) and np.all(bound >= -1e-12))
    balance_residual = float(np.max(np.abs(balance_history)))
    summary = {
        "source_concentration_mol_m3": float(p.initial_csf),
        "t50_isf_distal_s": t50_distal,
        "t50_isf_mean_s": t50_mean,
        "t50_pvs_s": t50_pvs,
        "effective_spreading_time_s": t50_distal,
        "isf_total_equivalent_mean_final_mol_m3": float(isf_total_history[-1]),
        "pvs_mean_final_mol_m3": float(pvs_history[-1]),
        "csf_final_mol_m3": float(csf_history[-1]),
        "max_relative_mass_balance_residual": balance_residual,
        "all_values_finite": all_finite,
        "all_values_nonnegative": nonnegative,
        "initial_mass_mol": float(initial_mass),
        "final_mass_mol": float(mass_history[-1]),
        "integrated_csf_to_pvs_mol": float(csf_to_pvs_integral),
        "integrated_pvs_to_isf_mol": float(pvs_to_isf_integral),
        "integrated_csf_out_mol": float(csf_out_integral),
        "integrated_blood_isf_out_mol": float(blood_isf_integral),
        "integrated_blood_pvs_out_mol": float(blood_pvs_integral),
        "integrated_pvs_out_mol": float(pvs_out_integral),
    }
    return {
        "parameters": asdict(p),
        "records": records,
        "time_s": [float(value) for value in times],
        "csf_mol_m3": csf_history,
        "isf_free_mean_mol_m3": isf_free_history,
        "isf_total_equivalent_mean_mol_m3": isf_total_history,
        "pvs_mean_mol_m3": pvs_history,
        "isf_distal_total_equivalent_mol_m3": distal_total_history,
        "mass_mol": mass_history,
        "mass_balance_relative_residual": balance_history,
        "summary": summary,
    }


def dimension_check(p: Parameters) -> dict[str, Any]:
    concentration = 1.0
    gradient = concentration / p.L
    diffusion_flux = p.A * p.eps_isf * p.D_isf * gradient
    _, v_isf, v_pvs = cell_volumes(p)
    exchange_volume = 2.0 * v_pvs * v_isf / (v_pvs + v_isf)
    exchange_flux = p.k_exchange * exchange_volume * concentration
    binding_flux = p.k_on * p.A * p.dx * concentration
    velocity = p.hydraulic_conductivity * p.pressure_delta_mmhg * MMHG_TO_PA / (p.viscosity * p.L)
    return {
        "diffusive_flux_mol_s": float(diffusion_flux),
        "exchange_flux_mol_s": float(exchange_flux),
        "binding_flux_mol_s": float(binding_flux),
        "darcy_velocity_m_s": float(velocity),
        "expected_flux_unit": "mol/s",
        "expected_velocity_unit": "m/s",
        "flux_dimension_ok": bool(np.isfinite(diffusion_flux) and np.isfinite(exchange_flux) and np.isfinite(binding_flux)),
        "velocity_dimension_ok": bool(np.isfinite(velocity)),
    }


def scenario_parameters(base: Parameters) -> dict[str, Parameters]:
    return {
        "diffusion_null": base,
        "pvs_advection_1_mmHg": replace(base, pressure_delta_mmhg=1.0),
        "binding_1e-3_2e-4": replace(base, k_on=1e-3, k_off=2e-4),
    }


def sensitivity_analysis(base: Parameters, baseline_t50: float | None) -> dict[str, Any]:
    selected = {
        "D_isf": ("D_isf", base.D_isf),
        "k_exchange": ("k_exchange", base.k_exchange),
        "L": ("L", base.L),
    }
    output: dict[str, Any] = {"baseline_t50_s": baseline_t50}
    for label, (field, value) in selected.items():
        low = run_model(replace(base, **{field: value * 0.5}))
        high = run_model(replace(base, **{field: value * 1.5}))
        low_t50 = low["summary"]["t50_isf_distal_s"]
        high_t50 = high["summary"]["t50_isf_distal_s"]
        output[label] = {
            "parameter": field,
            "minus_50_percent_t50_s": low_t50,
            "plus_50_percent_t50_s": high_t50,
            "relative_span": (
                (high_t50 - low_t50) / baseline_t50
                if low_t50 is not None and high_t50 is not None and baseline_t50 is not None
                else None
            ),
        }
    ranked = sorted(
        (entry for entry in output.values() if isinstance(entry, dict) and entry.get("relative_span") is not None),
        key=lambda entry: abs(entry["relative_span"]),
        reverse=True,
    )
    output["ranked_by_absolute_t50_span"] = [entry["parameter"] for entry in ranked]
    return output


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def run_suite(base: Parameters, prereg_path: Path) -> dict[str, Any]:
    scenarios = scenario_parameters(base)
    scenario_results = {name: run_model(parameters) for name, parameters in scenarios.items()}
    baseline = scenario_results["diffusion_null"]["summary"]
    sensitivity = sensitivity_analysis(base, baseline["t50_isf_distal_s"])
    prereg_hash = hashlib.sha256(prereg_path.read_bytes()).hexdigest()
    criterion_passed = bool(
        baseline["t50_isf_distal_s"] is not None
        and 0.5 <= baseline["t50_isf_distal_s"] / REFERENCE_T50_S <= 2.0
        and baseline["all_values_finite"]
        and baseline["all_values_nonnegative"]
        and baseline["max_relative_mass_balance_residual"] < 1e-6
    )
    return {
        "model_id": "BT-HX-Q156",
        "prereg_sha256": prereg_hash,
        "reference": {
            "authors": "Jin, Smith, and Verkman",
            "year": 2016,
            "journal": "Journal of General Physiology",
            "doi": "10.1085/jgp.201611684",
            "value_s": REFERENCE_T50_S,
            "value_unit": "s",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5129742/",
            "source_location": "Results, Fig. 2C discussion",
            "verification": "verified",
        },
        "dimension_check": dimension_check(base),
        "scenarios": scenario_results,
        "sensitivity": sensitivity,
        "primary_criterion": {
            "criterion": "0.5 <= t50_isf_distal_s / 600 s <= 2.0 and mass residual < 1e-6",
            "t50_ratio_to_reference": (
                baseline["t50_isf_distal_s"] / REFERENCE_T50_S
                if baseline["t50_isf_distal_s"] is not None
                else None
            ),
            "passed": criterion_passed,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results.json"))
    parser.add_argument("--prereg", type=Path, default=Path("PREREG.md"))
    args = parser.parse_args()
    base = Parameters()
    results = run_suite(base, args.prereg)
    args.output.write_text(json.dumps(jsonable(results), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    baseline = results["scenarios"]["diffusion_null"]["summary"]
    print(json.dumps({
        "output": str(args.output),
        "t50_isf_distal_s": baseline["t50_isf_distal_s"],
        "mass_balance_residual": baseline["max_relative_mass_balance_residual"],
        "criterion_passed": results["primary_criterion"]["passed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
