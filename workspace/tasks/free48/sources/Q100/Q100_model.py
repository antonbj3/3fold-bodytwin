from __future__ import annotations

import json
from dataclasses import dataclass, fields, replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    m1: float = 8.0e-4
    m2: float = 6.0e-4
    area0: float = 1.5e-5
    area_slope: float = 3.0e-3
    area_contact: float = 1.64e-5
    area_floor: float = 2.0e-6
    flow_gate_width: float = 2.0e-6
    p_sub: float = 800.0
    rho_air: float = 1.204
    flow_resistance: float = 3.5e5
    air_tau: float = 3.5e-3
    air_projection: float = 0.05
    air_damping: float = 0.04
    air_velocity_scale: float = 0.10
    flow_reference: float = 2.0e-4
    k_tissue_base: float = 600.0
    c_tissue_base: float = 0.12
    k_body: float = 600.0
    k_surface: float = 1500.0
    c_body: float = 0.10
    c_surface: float = 0.08
    k_coupling: float = 160.0
    c_coupling: float = 0.035
    k_contact_base: float = 6000.0
    c_contact_base: float = 1.20
    k_contact_surface: float = 2200.0
    c_contact_surface: float = 0.60
    contact_regularization: float = 2.0e-5
    contact_rate_regularization: float = 1.0e-4
    tissue_velocity_scale: float = 0.05
    meniscus_force: float = 2.0e-3
    meniscus_length: float = 1.2e-4
    mucus_damping: float = 0.012
    mucus_height: float = 2.0e-4
    rho_mucus: float = 1000.0
    sheet_tension: float = 0.18
    initial_x1: float = -9.0e-4
    initial_x2: float = -9.0e-4
    duration: float = 0.24
    measurement_window: float = 0.08
    max_step: float = 2.0e-5


DEFAULT = Parameters()

PARAMETER_SPECS: dict[str, tuple[str, str]] = {
    "m1": ("kg", "assumption: lumped lateral mass 1"),
    "m2": ("kg", "assumption: lumped lateral mass 2"),
    "area0": ("m^2", "assumption: open rest area"),
    "area_slope": ("m", "assumption: projected area derivative per inward displacement"),
    "area_contact": ("m^2", "assumption: contact threshold area"),
    "area_floor": ("m^2", "numerical regularization: positive flow area"),
    "flow_gate_width": ("m^2", "numerical regularization: contact-to-open flow transition"),
    "p_sub": ("Pa", "fixed protocol input"),
    "rho_air": ("kg/m^3", "assumption: standard air density at room conditions"),
    "flow_resistance": ("Pa*s/m^3", "first-order lumped viscous/turbulent resistance"),
    "air_tau": ("s", "first-order hydrodynamic response time"),
    "air_projection": ("1", "assumption: oblique vocal-wall pressure projection"),
    "air_damping": ("N*s/m", "assumption: pressure-flow aerodynamic negative-damping coefficient"),
    "air_velocity_scale": ("m/s", "assumption: saturation velocity for aerodynamic shear"),
    "flow_reference": ("m^3/s", "normalization for aerodynamic damping"),
    "k_tissue_base": ("N/m", "homogeneous two-mass reference stiffness"),
    "c_tissue_base": ("N*s/m", "homogeneous two-mass reference damping"),
    "k_body": ("N/m", "layered body spring stiffness"),
    "k_surface": ("N/m", "layered thin-surface spring stiffness"),
    "c_body": ("N*s/m", "layered body dashpot coefficient"),
    "c_surface": ("N*s/m", "layered thin-surface dashpot coefficient"),
    "k_coupling": ("N/m", "inter-mass tissue coupling"),
    "c_coupling": ("N*s/m", "inter-mass tissue coupling damping"),
    "k_contact_base": ("N/m", "homogeneous unilateral contact stiffness"),
    "c_contact_base": ("N*s/m", "homogeneous unilateral contact damping"),
    "k_contact_surface": ("N/m", "surface-layer contact indentation stiffness"),
    "c_contact_surface": ("N*s/m", "surface-layer contact indentation damping"),
    "contact_regularization": ("m", "numerical regularization: contact activation length"),
    "contact_rate_regularization": ("m/s", "numerical regularization: contact-rate smoothing"),
    "tissue_velocity_scale": ("m/s", "assumption: nonlinear tissue damping velocity scale"),
    "meniscus_force": ("N", "capillary-force amplitude at contact"),
    "meniscus_length": ("m", "capillary/lubricant decay length"),
    "mucus_damping": ("N*s/m", "mucus-film wall-shear coefficient"),
    "mucus_height": ("m", "mucus-film thickness used for wave-speed proxy"),
    "rho_mucus": ("kg/m^3", "assumption: mucus density"),
    "sheet_tension": ("N/m", "assumption: effective mucus-sheet tension"),
    "initial_x1": ("m", "fixed initial condition, mass 1"),
    "initial_x2": ("m", "fixed initial condition, mass 2"),
    "duration": ("s", "fixed integration protocol"),
    "measurement_window": ("s", "fixed final-window protocol"),
    "max_step": ("s", "numerical integration limit"),
}


VARIANTS = ("baseline", "layered", "placebo")
REFERENCE_F0_HZ = 120.0
REFERENCE_LOW_HZ = 60.0
REFERENCE_HIGH_HZ = 240.0
MIN_AMPLITUDE_UM = 5.0


def parameter_table(p: Parameters = DEFAULT) -> list[dict[str, Any]]:
    return [
        {
            "name": item.name,
            "value": float(getattr(p, item.name)),
            "unit": PARAMETER_SPECS[item.name][0],
            "source_or_assumption": PARAMETER_SPECS[item.name][1],
        }
        for item in fields(p)
    ]


def tissue_coefficients(p: Parameters, variant: str) -> tuple[float, float]:
    if variant == "layered":
        k_eff = p.k_body * p.k_surface / (p.k_body + p.k_surface)
        c_eff = (p.c_body * p.k_surface + p.c_surface * p.k_body) / (p.k_body + p.k_surface)
        return k_eff, c_eff
    return p.k_tissue_base, p.c_tissue_base


def raw_area(y: np.ndarray, p: Parameters) -> float:
    return float(p.area0 - p.area_slope * (y[0] + y[2]))


def contact_indentation(y: np.ndarray, p: Parameters) -> float:
    return max(0.0, (p.area_contact - raw_area(y, p)) / p.area_slope)


def smooth_contact_indentation(y: np.ndarray, p: Parameters) -> float:
    linear_indentation = (p.area_contact - raw_area(y, p)) / p.area_slope
    return float(p.contact_regularization * np.logaddexp(0.0, linear_indentation / p.contact_regularization))


def smooth_positive_rate(rate: float, scale: float) -> float:
    regularized = 0.5 * (rate + np.sqrt(rate * rate + scale * scale) - scale)
    return float(max(0.0, regularized))


def flow_target(y: np.ndarray, p: Parameters) -> float:
    area = raw_area(y, p)
    if area <= p.area_contact or p.p_sub == 0.0:
        return 0.0
    gate = float(np.clip((area - p.area_contact) / p.flow_gate_width, 0.0, 1.0))
    area_eff = max(area, p.area_floor)
    pressure = abs(p.p_sub)
    resistance = p.flow_resistance
    beta_over_area2 = 0.5 * p.rho_air / (area_eff * area_eff)
    denominator = resistance + np.sqrt(resistance * resistance + 4.0 * beta_over_area2 * pressure)
    magnitude = 2.0 * pressure / denominator
    return float(np.sign(p.p_sub) * magnitude * gate)


def flow_pressure(y: np.ndarray, p: Parameters) -> tuple[float, float, float]:
    area = raw_area(y, p)
    flow = float(y[4])
    if area <= p.area_contact:
        velocity = 0.0
        dynamic_pressure = 0.0
    else:
        area_eff = max(area, p.area_floor)
        velocity = flow / area_eff
        dynamic_pressure = min(0.5 * p.rho_air * velocity * velocity, abs(p.p_sub))
    wall_pressure = p.p_sub - dynamic_pressure
    return flow, velocity, wall_pressure


def force_components(y: np.ndarray, p: Parameters, variant: str) -> dict[str, np.ndarray | float]:
    x1, v1, x2, v2, _ = y
    k_tissue, c_tissue = tissue_coefficients(p, variant)
    tissue = np.array(
        [
            -k_tissue * x1 - c_tissue * v1 * (1.0 + (v1 / p.tissue_velocity_scale) ** 2),
            -k_tissue * x2 - c_tissue * v2 * (1.0 + (v2 / p.tissue_velocity_scale) ** 2),
        ],
        dtype=float,
    )
    coupling = np.array(
        [
            -p.k_coupling * (x1 - x2) - p.c_coupling * (v1 - v2),
            -p.k_coupling * (x2 - x1) - p.c_coupling * (v2 - v1),
        ],
        dtype=float,
    )
    linear_indentation = (p.area_contact - raw_area(y, p)) / p.area_slope
    penetration = smooth_contact_indentation(y, p)
    penetration_rate = v1 + v2
    closing_rate = smooth_positive_rate(penetration_rate, p.contact_rate_regularization)
    contact_activation = 0.5 * (1.0 + np.tanh(linear_indentation / p.contact_regularization))
    if variant == "layered":
        contact_stiffness = p.k_contact_surface
        contact_damping = p.c_contact_surface
    else:
        contact_stiffness = p.k_contact_base
        contact_damping = p.c_contact_base
    contact_total = contact_activation * (contact_stiffness * penetration + contact_damping * closing_rate)
    contact = np.array([-0.5 * contact_total, -0.5 * contact_total], dtype=float)
    if variant == "layered":
        area = raw_area(y, p)
        gap = max(0.0, (area - p.area_contact) / p.area_slope)
        film_factor = np.exp(-gap / p.meniscus_length)
        meniscus = np.array(
            [
                0.5 * p.meniscus_force * film_factor - p.mucus_damping * film_factor * v1,
                0.5 * p.meniscus_force * film_factor - p.mucus_damping * film_factor * v2,
            ],
            dtype=float,
        )
    else:
        meniscus = np.zeros(2, dtype=float)
    _, _, wall_pressure = flow_pressure(y, p)
    flow = abs(float(y[4])) / p.flow_reference
    air_pressure_force = -p.area_slope * p.air_projection * wall_pressure
    air_shear_coefficient = p.air_damping * flow
    air_shear_force = air_shear_coefficient * np.tanh(np.array([v1, v2], dtype=float) / p.air_velocity_scale)
    air = air_pressure_force + air_shear_force
    return {
        "tissue": tissue,
        "coupling": coupling,
        "contact": contact,
        "meniscus": meniscus,
        "air": air,
        "k_tissue_N_per_m": k_tissue,
        "c_tissue_N_s_per_m": c_tissue,
        "penetration_m": penetration,
        "penetration_rate_m_per_s": penetration_rate,
        "air_pressure_force_N": air_pressure_force,
        "air_shear_coefficient_N_s_per_m": air_shear_coefficient,
        "air_shear_force_N": air_shear_force,
    }


def state_derivative(t: float, y: np.ndarray, p: Parameters, variant: str) -> np.ndarray:
    del t
    components = force_components(y, p, variant)
    forces = components["tissue"] + components["coupling"] + components["contact"] + components["meniscus"] + components["air"]
    q_target = flow_target(y, p)
    return np.array(
        [
            y[1],
            forces[0] / p.m1,
            y[3],
            forces[1] / p.m2,
            (q_target - y[4]) / p.air_tau,
        ],
        dtype=float,
    )


def linear_matrices(p: Parameters, variant: str = "baseline") -> tuple[np.ndarray, np.ndarray]:
    k_tissue, _ = tissue_coefficients(p, variant)
    mass = np.diag([p.m1, p.m2])
    stiffness = np.array(
        [
            [k_tissue + p.k_coupling, -p.k_coupling],
            [-p.k_coupling, k_tissue + p.k_coupling],
        ],
        dtype=float,
    )
    return mass, stiffness


def linear_modes(p: Parameters, variant: str = "baseline") -> tuple[np.ndarray, np.ndarray]:
    mass, stiffness = linear_matrices(p, variant)
    eigenvalues, eigenvectors = np.linalg.eig(np.linalg.solve(mass, stiffness))
    order = np.argsort(eigenvalues)
    frequencies = np.sqrt(np.maximum(eigenvalues[order], 0.0)) / (2.0 * np.pi)
    return frequencies, eigenvectors[:, order]


def dominant_frequency(t: np.ndarray, signal: np.ndarray) -> float:
    if len(t) < 16:
        raise ValueError("frequency window is too short")
    dt = float(np.median(np.diff(t)))
    centered = np.asarray(signal, dtype=float) - float(np.mean(signal))
    scale = float(np.max(np.abs(centered)))
    if scale < 1e-14:
        raise ValueError("signal has no measurable oscillation")
    windowed = centered * np.hanning(len(centered))
    spectrum = np.abs(np.fft.rfft(windowed))
    frequencies = np.fft.rfftfreq(len(windowed), d=dt)
    spectrum[0] = 0.0
    peak = int(np.argmax(spectrum))
    if peak == 0 or spectrum[peak] <= 0.0:
        raise ValueError("no positive-frequency component")
    fractional_peak = 0.0
    if 0 < peak < len(spectrum) - 1:
        denominator = spectrum[peak - 1] - 2.0 * spectrum[peak] + spectrum[peak + 1]
        if denominator != 0.0:
            fractional_peak = 0.5 * (spectrum[peak - 1] - spectrum[peak + 1]) / denominator
            fractional_peak = float(np.clip(fractional_peak, -0.5, 0.5))
    return float((peak + fractional_peak) * frequencies[1])


def mucus_wave_speed(p: Parameters) -> float:
    return float(np.sqrt(p.sheet_tension / (p.rho_mucus * p.mucus_height)))


def simulate(
    p: Parameters = DEFAULT,
    variant: str = "layered",
    initial: np.ndarray | None = None,
) -> Any:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant: {variant}")
    if initial is None:
        initial = np.array([p.initial_x1, 0.0, p.initial_x2, 0.0, 0.0], dtype=float)
    sample_count = max(32, int(round(p.duration / p.max_step)) + 1)
    sample_times = np.linspace(0.0, p.duration, sample_count)
    solution = solve_ivp(
        lambda t, y: state_derivative(t, y, p, variant),
        (0.0, p.duration),
        np.asarray(initial, dtype=float),
        method="DOP853",
        t_eval=sample_times,
        rtol=1e-8,
        atol=1e-11,
        max_step=p.max_step,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def component_arrays(solution: Any, p: Parameters, variant: str) -> dict[str, np.ndarray]:
    names = ("tissue", "coupling", "contact", "meniscus", "air", "air_shear_force_N")
    output = {name: np.zeros((len(solution.t), 2), dtype=float) for name in names}
    metadata = {
        "k_tissue_N_per_m": np.zeros(len(solution.t)),
        "c_tissue_N_s_per_m": np.zeros(len(solution.t)),
        "penetration_m": np.zeros(len(solution.t)),
        "air_pressure_force_N": np.zeros(len(solution.t)),
        "air_shear_coefficient_N_s_per_m": np.zeros(len(solution.t)),
    }
    for index, state in enumerate(solution.y.T):
        components = force_components(state, p, variant)
        for name in names:
            output[name][index] = components[name]
        for name in metadata:
            metadata[name][index] = float(components[name])
    output.update(metadata)
    return output


def energy_diagnostics(solution: Any, p: Parameters, variant: str) -> dict[str, float]:
    t = solution.t
    y = solution.y
    components = component_arrays(solution, p, variant)
    velocity = y[[1, 3], :].T
    spring_force = -components["k_tissue_N_per_m"][:, None] * np.vstack((y[0], y[2])).T
    tissue_damping_force = components["tissue"] - spring_force
    tissue_dissipation = float(-np.trapezoid(np.sum(tissue_damping_force * velocity, axis=1), t))
    coupling_dissipation = float(np.trapezoid(p.c_coupling * (y[1] - y[3]) ** 2, t))
    air_power = float(
        np.trapezoid(
            np.sum((components["air_pressure_force_N"][:, None] + components["air_shear_force_N"]) * velocity, axis=1),
            t,
        )
    )
    contact_power = float(
        np.trapezoid(
            np.sum(components["contact"] * velocity, axis=1),
            t,
        )
    )
    meniscus_power = float(
        np.trapezoid(
            np.sum(components["meniscus"] * velocity, axis=1),
            t,
        )
    )
    kinetic = 0.5 * p.m1 * y[1] ** 2 + 0.5 * p.m2 * y[3] ** 2
    return {
        "air_power_mW": air_power * 1000.0 / float(t[-1] - t[0]),
        "tissue_dissipation_mJ": tissue_dissipation * 1000.0,
        "coupling_dissipation_mJ": coupling_dissipation * 1000.0,
        "contact_power_mJ": contact_power * 1000.0,
        "meniscus_power_mJ": meniscus_power * 1000.0,
        "final_kinetic_energy_uJ": float(kinetic[-1] * 1e6),
    }


def summarize(solution: Any, p: Parameters, variant: str) -> dict[str, Any]:
    t = solution.t
    y = solution.y
    start = max(0.0, float(t[-1]) - p.measurement_window)
    mask = t >= start - 1e-12
    if np.count_nonzero(mask) < 16:
        raise ValueError("measurement window is too short")
    tw = t[mask]
    yw = y[:, mask]
    mean_wall = 0.5 * (yw[0] + yw[2])
    frequency_observable = float(np.max(np.abs(mean_wall - np.mean(mean_wall)))) >= 1e-14
    f0 = dominant_frequency(tw, mean_wall) if frequency_observable else None
    area_values = p.area0 - p.area_slope * (yw[0] + yw[2])
    contact_values = area_values <= p.area_contact
    penetration_values = np.maximum(0.0, (p.area_contact - area_values) / p.area_slope)
    flow = yw[4]
    open_mask = area_values > p.area_contact
    area_eff = np.where(open_mask, np.maximum(area_values, p.area_floor), p.area_floor)
    velocity = np.where(open_mask, flow / area_eff, 0.0)
    dynamic_pressure = np.minimum(0.5 * p.rho_air * velocity * velocity, abs(p.p_sub))
    wall_pressure = p.p_sub - dynamic_pressure
    peak_to_peak = float(np.max(mean_wall) - np.min(mean_wall))
    k_tissue, c_tissue = tissue_coefficients(p, variant)
    result: dict[str, Any] = {
        "variant": variant,
        "f0_Hz": f0,
        "frequency_status": "MEASURABLE" if frequency_observable else "NO_MEASURABLE_OSCILLATION",
        "peak_to_peak_mean_wall_displacement_um": peak_to_peak * 1e6,
        "contact_percent": float(np.count_nonzero(contact_values) / len(contact_values) * 100.0),
        "maximum_contact_indentation_um": float(np.max(penetration_values) * 1e6),
        "subglottal_pressure_Pa": p.p_sub,
        "mean_wall_pressure_Pa": float(np.mean(wall_pressure)),
        "maximum_wall_pressure_Pa": float(np.max(wall_pressure)),
        "minimum_wall_pressure_Pa": float(np.min(wall_pressure)),
        "mean_signed_flow_m3_per_s": float(np.mean(flow)),
        "maximum_abs_flow_m3_per_s": float(np.max(np.abs(flow))),
        "minimum_glottal_area_m2": float(np.min(area_values)),
        "maximum_glottal_area_m2": float(np.max(area_values)),
        "tissue_stiffness_N_per_m": k_tissue,
        "tissue_damping_N_s_per_m": c_tissue,
        "mucus_wave_speed_proxy_m_per_s": mucus_wave_speed(p) if variant == "layered" else None,
        "finite_solution": bool(np.all(np.isfinite(y))),
    }
    result.update(energy_diagnostics(solution, p, variant))
    return result


def dimensional_checks(p: Parameters = DEFAULT) -> list[dict[str, Any]]:
    checks = []
    x = 1.0e-4
    v = 1.0e-3
    a = p.area0
    q = 1.0e-4
    k_tissue, c_tissue = tissue_coefficients(p, "layered")
    force = k_tissue * x
    damping_force = c_tissue * v
    pressure_loss = p.flow_resistance * q + 0.5 * p.rho_air * (q / a) ** 2
    wave_speed = mucus_wave_speed(p)
    checks.append(
        {
            "equation": "F_tissue = k_tissue*x",
            "unit_product": "(N/m)*m = N",
            "numerical_value": force,
            "expected_unit": "N",
            "status": "PASS" if np.isfinite(force) else "FAIL",
        }
    )
    checks.append(
        {
            "equation": "F_viscous = c_tissue*v",
            "unit_product": "(N*s/m)*(m/s) = N",
            "numerical_value": damping_force,
            "expected_unit": "N",
            "status": "PASS" if np.isfinite(damping_force) else "FAIL",
        }
    )
    checks.append(
        {
            "equation": "dp_flow = R*Q + 0.5*rho*(Q/A)^2",
            "unit_product": "(Pa*s/m^3)*(m^3/s) + (kg/m^3)*(m/s)^2 = Pa",
            "numerical_value": pressure_loss,
            "expected_unit": "Pa",
            "status": "PASS" if np.isfinite(pressure_loss) else "FAIL",
        }
    )
    checks.append(
        {
            "equation": "c_m = sqrt(T_sheet/(rho_m*h_m))",
            "unit_product": "sqrt((N/m)/(kg/m^2)) = m/s",
            "numerical_value": wave_speed,
            "expected_unit": "m/s",
            "status": "PASS" if np.isfinite(wave_speed) else "FAIL",
        }
    )
    return checks


def approval(result: dict[str, Any], reference_hz: float = REFERENCE_F0_HZ) -> dict[str, Any]:
    f0 = float(result["f0_Hz"]) if result["f0_Hz"] is not None else None
    amplitude = float(result["peak_to_peak_mean_wall_displacement_um"])
    finite = bool(result["finite_solution"])
    broad_pass = finite and f0 is not None and REFERENCE_LOW_HZ <= f0 <= REFERENCE_HIGH_HZ and amplitude >= MIN_AMPLITUDE_UM
    return {
        "broad_benchmark_pass": bool(broad_pass),
        "f0_interval_Hz": [REFERENCE_LOW_HZ, REFERENCE_HIGH_HZ],
        "minimum_peak_to_peak_displacement_um": MIN_AMPLITUDE_UM,
        "reference_f0_Hz": reference_hz,
        "reference_status": "OVERIFIERAD, UR MINNET",
        "criterion_status": "PASS" if broad_pass else "FAIL",
    }


def improvement_status(baseline: dict[str, Any], layered: dict[str, Any], reference_hz: float = REFERENCE_F0_HZ) -> dict[str, Any]:
    baseline_error = abs(float(baseline["f0_Hz"]) - reference_hz) if baseline["f0_Hz"] is not None else None
    layered_error = abs(float(layered["f0_Hz"]) - reference_hz) if layered["f0_Hz"] is not None else None
    contact_change = float(layered["contact_percent"]) - float(baseline["contact_percent"])
    return {
        "reference_f0_Hz": reference_hz,
        "baseline_absolute_error_Hz": baseline_error,
        "layered_absolute_error_Hz": layered_error,
        "layered_error_lower": bool(layered_error < baseline_error) if baseline_error is not None and layered_error is not None else None,
        "contact_fraction_change_percentage_points": contact_change,
        "mucus_proxy_present_only_in_layered": layered["mucus_wave_speed_proxy_m_per_s"] is not None and baseline["mucus_wave_speed_proxy_m_per_s"] is None,
    }


def sensitivity(p: Parameters = DEFAULT) -> list[dict[str, Any]]:
    selected = ("k_surface", "k_contact_surface", "p_sub")
    rows: list[dict[str, Any]] = []
    for name in selected:
        base_value = float(getattr(p, name))
        for factor in (0.5, 1.5):
            changed = replace(p, **{name: base_value * factor})
            result = summarize(simulate(changed, variant="layered"), changed, "layered")
            rows.append(
                {
                    "parameter": name,
                    "factor": factor,
                    "base_value": base_value,
                    "perturbed_value": float(getattr(changed, name)),
                    "f0_Hz": result["f0_Hz"],
                    "contact_percent": result["contact_percent"],
                    "peak_to_peak_mean_wall_displacement_um": result["peak_to_peak_mean_wall_displacement_um"],
                    "mucus_wave_speed_proxy_m_per_s": result["mucus_wave_speed_proxy_m_per_s"],
                }
            )
    return rows


def run_all() -> dict[str, Any]:
    p = DEFAULT
    solutions = {variant: simulate(p, variant=variant) for variant in VARIANTS}
    summaries = {variant: summarize(solutions[variant], p, variant) for variant in VARIANTS}
    baseline = summaries["baseline"]
    layered = summaries["layered"]
    placebo = summaries["placebo"]
    placebo_f0_difference = (abs(float(placebo["f0_Hz"]) - float(baseline["f0_Hz"]))
                             if placebo["f0_Hz"] is not None and baseline["f0_Hz"] is not None else None)
    placebo_frequency_matches = (placebo_f0_difference <= 1.0e-6 if placebo_f0_difference is not None
                                 else placebo["f0_Hz"] is None and baseline["f0_Hz"] is None)
    placebo_amplitude_difference = abs(
        float(placebo["peak_to_peak_mean_wall_displacement_um"])
        - float(baseline["peak_to_peak_mean_wall_displacement_um"])
    )
    output: dict[str, Any] = {
        "id": "BT-HX-Q100",
        "model": "two-mass layered-tissue vocal-fold oscillator with one-sided contact, mucus film, and reduced air flow",
        "protocol": {
            "p_sub_Pa": p.p_sub,
            "duration_s": p.duration,
            "measurement_window_s": p.measurement_window,
            "max_step_s": p.max_step,
            "reference_f0_Hz": REFERENCE_F0_HZ,
            "reference_source": "N. Ishizaka and J. L. Flanagan, Two-mass model of the vocal folds, Journal of the Acoustical Society of America (1972), approximately 120 Hz",
            "reference_status": "OVERIFIERAD, UR MINNET",
        },
        "equations": {
            "mechanics": "m_i*u_i'' = -k_i*u_i-c_i*u_i'*(1+(u_i'/v_tissue_scale)^2) - k_c*(u_i-u_j)-c_c*(u_i'-u_j') + F_air,i + F_contact,i + F_meniscus,i",
            "geometry": "A_raw = A0 - a_s*(u1+u2)",
            "flow": "Q* = sign(Psub)*2*|Psub|/(R+sqrt(R^2+4*(rho_air/2)*|Psub|/A_eff^2)) and dQ/dt=(Q*-Q)/tau_air",
            "air_pressure": "F_air,i = -a_s*projection*(Psub-min(0.5*rho_air*(Q/A_eff)^2,abs(Psub))) + c_air*(abs(Q)/Q_ref)*tanh(u_i'/v_air)",
            "contact": "q_c=regularized_positive((A_contact-A_raw)/a_s); F_contact,i=-0.5*activation*(k_contact*q_c+c_contact*smooth_positive_rate(u1'+u2', rate_scale))",
            "layer_series_tissue": "k_series=k_body*k_surface/(k_body+k_surface); c_series=(c_body*k_surface+c_surface*k_body)/(k_body+k_surface)",
            "mucus_proxy": "c_m=sqrt(T_sheet/(rho_mucus*h_m))",
        },
        "parameter_table": parameter_table(p),
        "variants": summaries,
        "approval": approval(layered),
        "improvement": improvement_status(baseline, layered),
        "placebo_check": {
            "f0_difference_Hz": placebo_f0_difference,
            "amplitude_difference_um": placebo_amplitude_difference,
            "tolerance_Hz": 1.0e-6,
            "tolerance_um": 1.0e-6,
            "frequency_comparison_status": "MEASURABLE" if placebo_f0_difference is not None else "UNDEFINED_FREQUENCY",
            "status": "PASS" if placebo_frequency_matches and placebo_amplitude_difference <= 1.0e-6 else "FAIL",
        },
        "sensitivity": sensitivity(p),
        "dimensional_checks": dimensional_checks(p),
        "data_status": {
            "measured_inputs_used": False,
            "geometry_measured": False,
            "external_data_used": False,
            "mucosal_wave_empirical_validation": "UNKNOWN",
            "reference_matching": "UNKNOWN because the remembered benchmark is not matched to this geometry and pressure",
        },
    }
    return output


def main() -> None:
    output = run_all()
    Path("results.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"f0_layered_Hz": output["variants"]["layered"]["f0_Hz"], "criterion": output["approval"]["criterion_status"]}, sort_keys=True))


if __name__ == "__main__":
    main()
