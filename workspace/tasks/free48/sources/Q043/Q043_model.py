from dataclasses import dataclass, replace
from pathlib import Path
import json
import math
import numpy as np

EQUATIONS = {
    "recruitment": "q_m(t)=sigmoid((d(t)-r_m)/s_r)",
    "rate": "nu_m(t)=q_m(t)[f_min+(f_peak,m-f_min)clip((d(t)-r_m)/(1-r_m),0,1)]",
    "spikes": "intervals_i~Gamma(1/CV^2, mean_interval/(1/CV^2)), mean_interval=1/nu_m",
    "force": "F_m(t)=G_m sum_j h_F(t-t_j), h_F(tau)=(tau/tau_p^2)exp(1-tau/tau_p) 1(tau>=0)",
    "emg": "E(t)=sum_m b_m g_m sum_j h_E(t-t_j), h_E=finite bipolar Gaussian derivative",
    "observation": "F_obs=F+noise_F, E_obs=lowpass(abs(E))+noise_E",
    "joint_error": "R_joint=sqrt(mean(((F_z-F_hat)^2+(E_z-E_hat)^2)/2))",
}

REFERENCES = [
    {
        "citation": "Kutch JJ, Kuo AD, Rymer WZ. Extraction of Individual Muscle Mechanical Action From Endpoint Force. Journal of Neurophysiology. 2010;103(6):3535-3546.",
        "doi": "10.1152/jn.00956.2009",
        "locations": [
            {
                "location": "Table 1",
                "values": {
                    "motor_units": 120,
                    "action_potential_duration_ms": [5, 20],
                    "longest_contraction_time_ms": 90,
                    "minimum_firing_rate_hz": 8,
                    "first_unit_peak_rate_hz": 45,
                    "peak_rate_difference_hz": 10,
                    "isi_cv": 0.2,
                    "recruitment_threshold_range_model_units": 30
                }
            },
            {
                "location": "Figure 5B",
                "values": {
                    "action_potential_duration_ms_to_correlation": {"5": 0.91, "10": 0.73, "15": 0.53, "20": 0.30}
                }
            },
            {
                "location": "Figure 5A",
                "value": "peak near 50% MVC"
            }
        ]
    }
]

@dataclass(frozen=True)
class Parameters:
    n_mus: int = 120
    dt: float = 0.001
    duration_s: float = 4.0
    seed: int = 43
    recruitment_center: float = 0.60
    recruitment_spread: float = 0.30
    recruitment_width: float = 0.035
    drive_mean: float = 0.60
    drive_amplitude: float = 0.18
    drive_frequency_hz: float = 0.45
    drive_noise: float = 0.012
    min_rate_hz: float = 8.0
    first_peak_rate_hz: float = 45.0
    peak_rate_step_hz: float = 10.0
    isi_cv: float = 0.20
    total_peak_force_n: float = 100.0
    reference_rate_hz: float = 20.0
    force_tau_min_s: float = 0.030
    force_tau_max_s: float = 0.090
    action_potential_duration_ms: float = 5.0
    emg_event_amplitude_mv: float = 0.035
    emg_lowpass_tau_s: float = 0.012
    muscle_half_width_mm: float = 30.0
    muscle_half_length_mm: float = 35.0
    electrode_half_separation_mm: float = 5.0
    geometry_length_mm: float = 12.0
    force_noise_fraction: float = 0.02
    emg_noise_fraction: float = 0.05
    train_fraction: float = 0.60
    pooled_error_threshold: float = 0.10
    gain_threshold: float = 0.20


def parameter_table(p: Parameters):
    return [
        {"name": "n_mus", "value": p.n_mus, "unit": "count", "source": "Kutch et al. 2010, Table 1"},
        {"name": "dt", "value": p.dt, "unit": "s", "source": "frozen synthetic protocol"},
        {"name": "duration_s", "value": p.duration_s, "unit": "s", "source": "frozen synthetic protocol"},
        {"name": "recruitment_spread", "value": p.recruitment_spread, "unit": "normalized drive", "source": "Kutch et al. 2010, Table 1, range 30 model units; scaled to 0.30"},
        {"name": "recruitment_width", "value": p.recruitment_width, "unit": "normalized drive", "source": "smooth recruitment assumption"},
        {"name": "min_rate_hz", "value": p.min_rate_hz, "unit": "Hz", "source": "Kutch et al. 2010, Table 1"},
        {"name": "first_peak_rate_hz", "value": p.first_peak_rate_hz, "unit": "Hz", "source": "Kutch et al. 2010, Table 1"},
        {"name": "peak_rate_step_hz", "value": p.peak_rate_step_hz, "unit": "Hz", "source": "Kutch et al. 2010, Table 1"},
        {"name": "isi_cv", "value": p.isi_cv, "unit": "dimensionless", "source": "Kutch et al. 2010, Table 1"},
        {"name": "total_peak_force_n", "value": p.total_peak_force_n, "unit": "N", "source": "normalization assumption; not a patient measurement"},
        {"name": "reference_rate_hz", "value": p.reference_rate_hz, "unit": "Hz", "source": "force normalization assumption"},
        {"name": "force_tau_min_s", "value": p.force_tau_min_s, "unit": "s", "source": "Kutch et al. 2010, Table 1, range 3 and longest 90 ms; mapped to 30–90 ms"},
        {"name": "force_tau_max_s", "value": p.force_tau_max_s, "unit": "s", "source": "Kutch et al. 2010, Table 1, longest contraction time 90 ms"},
        {"name": "action_potential_duration_ms", "value": p.action_potential_duration_ms, "unit": "ms", "source": "Kutch et al. 2010, Table 1 and Figure 5B"},
        {"name": "emg_event_amplitude_mv", "value": p.emg_event_amplitude_mv, "unit": "mV", "source": "amplitude assumption; geometry is relative"},
        {"name": "emg_lowpass_tau_s", "value": p.emg_lowpass_tau_s, "unit": "s", "source": "rectified EMG observation assumption"},
        {"name": "muscle_half_width_mm", "value": p.muscle_half_width_mm, "unit": "mm", "source": "geometry assumption"},
        {"name": "muscle_half_length_mm", "value": p.muscle_half_length_mm, "unit": "mm", "source": "geometry assumption"},
        {"name": "electrode_half_separation_mm", "value": p.electrode_half_separation_mm, "unit": "mm", "source": "differential electrode assumption"},
        {"name": "geometry_length_mm", "value": p.geometry_length_mm, "unit": "mm", "source": "volume-conduction length assumption"},
        {"name": "force_noise_fraction", "value": p.force_noise_fraction, "unit": "fraction of force SD", "source": "observation-noise assumption"},
        {"name": "emg_noise_fraction", "value": p.emg_noise_fraction, "unit": "fraction of EMG SD", "source": "observation-noise assumption"},
    ]


def sigmoid(x):
    x = np.clip(x, -40.0, 40.0)
    return 1.0 / (1.0 + np.exp(-x))


def make_drive(p: Parameters, rng: np.random.Generator):
    n = int(round(p.duration_s / p.dt))
    t = np.arange(n, dtype=float) * p.dt
    slow = p.drive_amplitude * np.sin(2.0 * np.pi * p.drive_frequency_hz * t)
    fast = 0.30 * p.drive_amplitude * np.sin(2.0 * np.pi * 2.3 * p.drive_frequency_hz * t + 0.7)
    noise = rng.normal(0.0, p.drive_noise, n)
    return np.clip(p.drive_mean + slow + fast + noise, 0.0, 1.0)


def make_thresholds(p: Parameters):
    if p.recruitment_spread <= 0.0:
        return np.full(p.n_mus, p.recruitment_center, dtype=float)
    rank = np.linspace(0.0, 1.0, p.n_mus)
    return p.recruitment_center + (rank - 0.5) * p.recruitment_spread


def make_rates(drive, thresholds, p: Parameters):
    excess = np.clip(drive[None, :] - thresholds[:, None], 0.0, 1.0)
    recruitment = sigmoid(excess / p.recruitment_width)
    rank = np.linspace(0.0, 1.0, p.n_mus)
    peak_rate = p.first_peak_rate_hz + rank * p.peak_rate_step_hz
    rate_fraction = excess / np.maximum(1.0 - thresholds[:, None], 1e-6)
    rate_fraction = np.clip(rate_fraction, 0.0, 1.0)
    return recruitment * (p.min_rate_hz + (peak_rate[:, None] - p.min_rate_hz) * rate_fraction)


def make_geometry(p: Parameters):
    rank = np.linspace(0.0, 1.0, p.n_mus)
    x = np.linspace(-p.muscle_half_length_mm, p.muscle_half_length_mm, p.n_mus)
    y = 0.65 * np.sin(2.0 * np.pi * rank + 0.4) * p.muscle_half_width_mm
    y += 0.25 * np.cos(5.0 * np.pi * rank) * p.muscle_half_width_mm
    x_pos = x + p.electrode_half_separation_mm
    x_neg = x - p.electrode_half_separation_mm
    d_pos = np.sqrt(x_pos * x_pos + y * y)
    d_neg = np.sqrt(x_neg * x_neg + y * y)
    raw = np.exp(-d_pos / p.geometry_length_mm) - np.exp(-d_neg / p.geometry_length_mm)
    return raw / max(float(np.mean(np.abs(raw))), 1e-12)


def force_weights(p: Parameters):
    rank = np.linspace(0.0, 1.0, p.n_mus)
    weights = 0.15 + 0.85 * rank ** 1.6
    return weights / np.sum(weights)


def force_taus(p: Parameters):
    rank = np.linspace(0.0, 1.0, p.n_mus)
    return p.force_tau_min_s + (p.force_tau_max_s - p.force_tau_min_s) * rank


def force_kernel(tau_s: float, dt: float):
    length = max(2, int(np.ceil(8.0 * tau_s / dt)))
    tau = np.arange(length + 1, dtype=float) * dt
    kernel = tau / (tau_s * tau_s) * np.exp(1.0 - tau / tau_s)
    area = np.trapz(kernel, dx=dt)
    return kernel / max(float(area), 1e-12)


def emg_kernel(duration_ms: float, dt: float):
    duration_s = max(duration_ms / 1000.0, 4.0 * dt)
    sigma = max(duration_s / 4.0, dt)
    half = max(2, int(np.ceil(2.0 * sigma / dt)))
    lag = np.arange(-half, half + 1, dtype=float) * dt
    kernel = -(lag / sigma) * np.exp(-0.5 * (lag / sigma) ** 2)
    peak = float(np.max(np.abs(kernel)))
    return lag, kernel / max(peak, 1e-12)


def generate_spikes(rates, p: Parameters, rng: np.random.Generator, homogeneous=False):
    n_mus, n_time = rates.shape
    spikes = np.zeros((n_mus, n_time), dtype=np.int8)
    if homogeneous:
        for t in range(1, n_time):
            if rates[0, t] > 0.0 and rng.random() < rates[0, t] * p.dt:
                spikes[:, t] = 1
        return spikes
    if p.isi_cv <= 0.0:
        for m in range(n_mus):
            next_time = rng.random() / max(float(rates[m, 0]), 1e-9)
            while next_time < p.duration_s:
                index = int(next_time / p.dt)
                if index >= n_time:
                    break
                spikes[m, index] = 1
                local_rate = max(float(rates[m, min(index, n_time - 1)]), 1e-9)
                next_time += 1.0 / local_rate
        return spikes
    shape = 1.0 / (p.isi_cv * p.isi_cv)
    for m in range(n_mus):
        local_rate = max(float(rates[m, 0]), 1e-9)
        next_time = rng.random() / local_rate
        while next_time < p.duration_s:
            index = int(next_time / p.dt)
            if index >= n_time:
                break
            spikes[m, index] = 1
            local_rate = max(float(rates[m, min(index, n_time - 1)]), 1e-9)
            mean_interval = 1.0 / local_rate
            interval = rng.gamma(shape, mean_interval / shape)
            next_time += max(interval, p.dt * 0.25)
    return spikes


def add_force_events(spikes, gains_n_s, taus_s, p: Parameters):
    n_time = spikes.shape[1]
    force = np.zeros(n_time, dtype=float)
    for m in range(spikes.shape[0]):
        events = np.flatnonzero(spikes[m])
        if events.size == 0:
            continue
        kernel = force_kernel(float(taus_s[m]), p.dt)
        for index in events:
            count = min(kernel.size, n_time - int(index))
            force[int(index):int(index) + count] += float(gains_n_s[m]) * kernel[:count]
    return force


def add_emg_events(spikes, gains_mv, p: Parameters):
    n_time = spikes.shape[1]
    lag, kernel = emg_kernel(p.action_potential_duration_ms, p.dt)
    half = kernel.size // 2
    emg = np.zeros(n_time, dtype=float)
    for m in range(spikes.shape[0]):
        events = np.flatnonzero(spikes[m])
        for index in events:
            start = max(0, int(index) - half)
            stop = min(n_time, int(index) + half + 1)
            offset = start - (int(index) - half)
            emg[start:stop] += float(gains_mv[m]) * kernel[offset:offset + (stop - start)]
    return emg


def lowpass_rectified(raw, p: Parameters):
    n_time = raw.size
    output = np.zeros(n_time, dtype=float)
    alpha = 1.0 - math.exp(-p.dt / p.emg_lowpass_tau_s)
    state = abs(float(raw[0])) if n_time else 0.0
    for i in range(n_time):
        state = (1.0 - alpha) * state + alpha * abs(float(raw[i]))
        output[i] = state
    return output


def simulate(p: Parameters, mode="nominal"):
    n_time = int(round(p.duration_s / p.dt))
    rng = np.random.default_rng(p.seed)
    drive = make_drive(p, rng)
    if mode == "zero":
        zeros = np.zeros(n_time, dtype=float)
        return {
            "drive": drive,
            "thresholds": make_thresholds(p),
            "rates": np.zeros((p.n_mus, n_time)),
            "spikes": np.zeros((p.n_mus, n_time), dtype=np.int8),
            "force_true": zeros.copy(),
            "emg_true": zeros.copy(),
            "force_observed": zeros.copy(),
            "emg_observed": zeros.copy(),
            "geometry": np.zeros(p.n_mus),
            "force_gains_n_s": np.zeros(p.n_mus),
            "emg_gains_mv": np.zeros(p.n_mus),
            "mode": mode,
        }
    homogeneous = mode == "homogeneous"
    thresholds = np.full(p.n_mus, p.recruitment_center, dtype=float) if homogeneous else make_thresholds(p)
    excess = np.clip(drive[None, :] - thresholds[:, None], 0.0, 1.0)
    recruitment = sigmoid(excess / p.recruitment_width)
    if homogeneous:
        rank = np.zeros(p.n_mus)
        peak_rate = np.full(p.n_mus, p.first_peak_rate_hz)
    else:
        rank = np.linspace(0.0, 1.0, p.n_mus)
        peak_rate = p.first_peak_rate_hz + rank * p.peak_rate_step_hz
    rate_fraction = np.clip(excess / np.maximum(1.0 - thresholds[:, None], 1e-6), 0.0, 1.0)
    rates = recruitment * (p.min_rate_hz + (peak_rate[:, None] - p.min_rate_hz) * rate_fraction)
    if homogeneous:
        rates = np.tile(rates[0:1, :], (p.n_mus, 1))
    spikes = generate_spikes(rates, p, rng, homogeneous=homogeneous)
    if homogeneous:
        force_gains = np.full(p.n_mus, p.total_peak_force_n / (p.n_mus * p.reference_rate_hz))
        geometry = np.ones(p.n_mus, dtype=float)
        emg_gains = np.full(p.n_mus, p.emg_event_amplitude_mv)
    else:
        weights = force_weights(p)
        force_gains = p.total_peak_force_n * weights / p.reference_rate_hz
        geometry = make_geometry(p)
        emg_gains = p.emg_event_amplitude_mv * (0.5 + rank) * geometry
    taus = np.full(p.n_mus, (p.force_tau_min_s + p.force_tau_max_s) * 0.5) if homogeneous else force_taus(p)
    force_true = add_force_events(spikes, force_gains, taus, p)
    emg_raw = add_emg_events(spikes, emg_gains, p)
    if homogeneous:
        emg_true = lowpass_rectified(force_true, p)
        emg_scale = max(float(np.std(force_true)), 1e-12)
        emg_true *= p.emg_event_amplitude_mv / emg_scale
    else:
        emg_true = lowpass_rectified(emg_raw, p)
    force_scale = float(np.std(force_true))
    emg_scale = float(np.std(emg_true))
    force_observed = force_true.copy()
    emg_observed = emg_true.copy()
    if mode != "zero" and force_scale > 0.0:
        force_observed += rng.normal(0.0, p.force_noise_fraction * force_scale, n_time)
    if mode != "zero" and emg_scale > 0.0:
        emg_observed += rng.normal(0.0, p.emg_noise_fraction * emg_scale, n_time)
    return {
        "drive": drive,
        "thresholds": thresholds,
        "rates": rates,
        "spikes": spikes,
        "force_true": force_true,
        "emg_true": emg_true,
        "force_observed": force_observed,
        "emg_observed": emg_observed,
        "geometry": geometry,
        "force_gains_n_s": force_gains,
        "emg_gains_mv": emg_gains,
        "mode": mode,
    }


def rank_one_joint_error(force, emg, train_fraction):
    force = np.asarray(force, dtype=float)
    emg = np.asarray(emg, dtype=float)
    n_time = force.size
    n_train = max(2, min(n_time - 1, int(round(train_fraction * n_time))))
    train = np.vstack((force[:n_train], emg[:n_train]))
    test = np.vstack((force[n_train:], emg[n_train:]))
    mean = np.mean(train, axis=1, keepdims=True)
    scale = np.std(train, axis=1, keepdims=True)
    scale[scale < 1e-12] = 1.0
    train_z = (train - mean) / scale
    test_z = (test - mean) / scale
    covariance = (train_z @ train_z.T) / max(train_z.shape[1] - 1, 1)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    observation_basis = eigenvectors[:, int(np.argmax(eigenvalues))]
    test_projection = observation_basis @ test_z
    prediction = np.outer(observation_basis, test_projection)
    residual = test_z - prediction
    return float(np.sqrt(np.mean(residual * residual))), prediction, mean, scale


def phase_shuffle(values, seed):
    values = np.asarray(values, dtype=float)
    if values.size < 4:
        return values.copy()
    centered = values - float(np.mean(values))
    spectrum = np.fft.rfft(centered)
    rng = np.random.default_rng(seed)
    phases = rng.uniform(0.0, 2.0 * np.pi, spectrum.size)
    phases[0] = 0.0
    if spectrum.size > 1:
        phases[-1] = 0.0
    shuffled = np.fft.irfft(spectrum * np.exp(1j * phases), n=values.size)
    return float(np.mean(values)) + shuffled


def cross_channel_lag_ms(force, emg, dt_s, max_lag_ms=200.0):
    force = np.asarray(force, dtype=float)
    emg = np.asarray(emg, dtype=float)
    f = force - np.mean(force)
    e = emg - np.mean(emg)
    if np.std(f) < 1e-12 or np.std(e) < 1e-12:
        return 0.0
    max_lag = int(round(max_lag_ms / 1000.0 / dt_s))
    best_lag = 0
    best_value = -np.inf
    denom = float(np.linalg.norm(f) * np.linalg.norm(e))
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            value = float(np.dot(f[-lag:], e[:lag])) / denom
        elif lag > 0:
            value = float(np.dot(f[:-lag], e[lag:])) / denom
        else:
            value = float(np.dot(f, e)) / denom
        if value > best_value:
            best_value = value
            best_lag = lag
    return float(best_lag * dt_s * 1000.0)


def finite_float(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("non-finite model output")
    return value


def evaluate_case(simulation, p: Parameters):
    force_true = simulation["force_true"]
    emg_true = simulation["emg_true"]
    force_observed = simulation["force_observed"]
    emg_observed = simulation["emg_observed"]
    force_scale = max(float(np.std(force_true)), 1e-12)
    emg_scale = max(float(np.std(emg_true)), 1e-12)
    pooled_error, _, _, _ = rank_one_joint_error(force_observed, emg_observed, p.train_fraction)
    mu_error = float(np.sqrt(np.mean(np.concatenate((((force_observed - force_true) / force_scale) ** 2, ((emg_observed - emg_true) / emg_scale) ** 2)) * 0.5)))
    mu_error = finite_float(mu_error)
    if pooled_error <= 1e-12 and mu_error <= 1e-12:
        gain = 0.0
    else:
        gain = finite_float(1.0 - mu_error / max(pooled_error, 1e-12))
    placebo_emg = np.maximum(phase_shuffle(emg_observed, p.seed + 701), 0.0)
    placebo_error, _, _, _ = rank_one_joint_error(force_observed, placebo_emg, p.train_fraction)
    correlation = 0.0
    if np.std(force_observed) > 1e-12 and np.std(emg_observed) > 1e-12:
        correlation = finite_float(np.corrcoef(force_observed, emg_observed)[0, 1])
    detail_needed = bool(pooled_error >= p.pooled_error_threshold and gain >= p.gain_threshold)
    geometry = np.asarray(simulation["geometry"], dtype=float)
    force_gains = np.asarray(simulation["force_gains_n_s"], dtype=float)
    visibility_cv = 0.0
    if np.mean(np.abs(geometry)) > 1e-12:
        visibility_cv = finite_float(np.std(np.abs(geometry), ddof=1) / np.mean(np.abs(geometry)))
    force_capacity_cv = 0.0
    if np.mean(force_gains) > 1e-12:
        force_capacity_cv = finite_float(np.std(force_gains, ddof=1) / np.mean(force_gains))
    return {
        "mode": simulation["mode"],
        "force_true_peak_n": finite_float(np.max(np.abs(force_true))),
        "emg_true_rms_mv": finite_float(np.sqrt(np.mean(np.square(emg_true)))),
        "force_emg_correlation": correlation,
        "force_lead_lag_ms": finite_float(cross_channel_lag_ms(force_observed, emg_observed, p.dt)),
        "pooled_joint_nrmse_test": finite_float(pooled_error),
        "mu_joint_nrmse_test": mu_error,
        "mu_gain": gain,
        "placebo_joint_nrmse_test": finite_float(placebo_error),
        "placebo_minus_aligned_pooled": finite_float(placebo_error - pooled_error),
        "geometry_visibility_cv": visibility_cv,
        "force_capacity_cv": force_capacity_cv,
        "detail_needed": detail_needed,
    }


def unit_check(p: Parameters):
    if p.n_mus < 2:
        raise ValueError("n_mus must be at least two")
    if p.dt <= 0.0 or p.duration_s <= 0.0:
        raise ValueError("time parameters must be positive")
    if p.recruitment_spread < 0.0 or p.recruitment_width <= 0.0:
        raise ValueError("recruitment parameters are invalid")
    if not 0.0 < p.train_fraction < 1.0:
        raise ValueError("train_fraction must lie between zero and one")
    if p.force_tau_min_s <= 0.0 or p.force_tau_max_s < p.force_tau_min_s:
        raise ValueError("force time constants are invalid")
    if p.electrode_half_separation_mm <= 0.0 or p.geometry_length_mm <= 0.0:
        raise ValueError("electrode geometry is invalid")
    kernel = force_kernel(p.force_tau_min_s, p.dt)
    area = float(np.trapz(kernel, dx=p.dt))
    if not math.isfinite(area) or abs(area - 1.0) > 0.02:
        raise ValueError("force kernel does not pass the unit-integral check")
    return {
        "force_kernel_integral_s": finite_float(area),
        "force_kernel_unit": "1/s",
        "force_gain_unit": "N s",
        "force_output_unit": "N",
        "emg_kernel_unit": "dimensionless per event",
        "emg_gain_unit": "mV per event",
        "emg_output_unit": "mV",
    }


def case_result(p: Parameters, mode="nominal"):
    unit_check(p)
    simulation = simulate(p, mode=mode)
    result = evaluate_case(simulation, p)
    result["parameters"] = {
        "n_mus": p.n_mus,
        "dt_s": p.dt,
        "duration_s": p.duration_s,
        "action_potential_duration_ms": p.action_potential_duration_ms,
        "geometry_length_mm": p.geometry_length_mm,
        "recruitment_spread": p.recruitment_spread,
        "seed": p.seed,
    }
    return result, simulation


def run_benchmark(p: Parameters | None = None):
    p = p or Parameters()
    unit_check(p)
    nominal, nominal_simulation = case_result(p, mode="nominal")
    source_control_p = replace(p, action_potential_duration_ms=20.0)
    source_control, _ = case_result(source_control_p, mode="nominal")
    homogeneous_control_p = replace(p, force_noise_fraction=0.0, emg_noise_fraction=0.0)
    homogeneous, _ = case_result(homogeneous_control_p, mode="homogeneous")
    zero, _ = case_result(p, mode="zero")
    sensitivity = []
    for name in ("geometry_length_mm", "action_potential_duration_ms", "recruitment_spread"):
        for factor in (0.5, 1.5):
            changed = replace(p, **{name: getattr(p, name) * factor})
            result, _ = case_result(changed, mode="nominal")
            result["changed_parameter"] = name
            result["factor"] = factor
            sensitivity.append(result)
    phase_shuffled_emg = np.maximum(phase_shuffle(nominal_simulation["emg_observed"], p.seed + 701), 0.0)
    placebo_error, _, _, _ = rank_one_joint_error(nominal_simulation["force_observed"], phase_shuffled_emg, p.train_fraction)
    return {
        "id": "BT-HX-Q043",
        "model_version": "first-principles-v1",
        "status": "synthetic mechanistic benchmark; not measured data",
        "observation_budget": {
            "force_channels": 1,
            "emg_channels": 1,
            "sample_rate_hz": 1.0 / p.dt,
            "duration_s": p.duration_s,
            "train_fraction": p.train_fraction,
        },
        "equations": EQUATIONS,
        "references": REFERENCES,
        "validity_notes": {
            "data_status": "All signals are synthetic predictions; no measured data are fitted.",
            "mu_error": "The full motor-unit result is a frozen forward-model oracle scored against the generated pre-noise channels.",
            "lag_sign": "Positive force_lead_lag_ms means EMG leads force; negative means force leads EMG.",
            "homogeneous_control": "Noise-free matched-temporal-response structural limit.",
            "phase_shuffle_scope": "Nominal aligned versus phase-shuffled EMG only."
        },
        "parameter_table": parameter_table(p),
        "unit_check": unit_check(p),
        "frozen_criterion": {
            "pooled_error_threshold": p.pooled_error_threshold,
            "mu_gain_threshold": p.gain_threshold,
            "detail_rule": "pooled_joint_nrmse_test >= 0.10 and mu_gain >= 0.20",
        },
        "nominal": nominal,
        "source_control_action_potential_20ms": source_control,
        "controls": {
            "homogeneous_visibility": homogeneous,
            "zero_drive": zero,
            "phase_shuffled_emg_pooled_joint_nrmse_test": finite_float(placebo_error),
            "phase_shuffled_emg_minus_aligned": finite_float(placebo_error - nominal["pooled_joint_nrmse_test"]),
        },
        "sensitivity": sensitivity,
    }


def write_results(result, path=None):
    output_path = Path(path) if path is not None else Path(__file__).with_name("results.json")
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return output_path


def main():
    result = run_benchmark()
    path = write_results(result)
    print(json.dumps({
        "output": str(path),
        "nominal_pooled_joint_nrmse": result["nominal"]["pooled_joint_nrmse_test"],
        "nominal_mu_joint_nrmse": result["nominal"]["mu_joint_nrmse_test"],
        "nominal_mu_gain": result["nominal"]["mu_gain"],
        "nominal_detail_needed": result["nominal"]["detail_needed"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
