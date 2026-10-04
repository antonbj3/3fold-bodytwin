from __future__ import annotations

import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import solve_banded


@dataclass(frozen=True)
class Parameters:
    D_um2_per_s: float = 55.0
    k_clear_per_s: float = 0.0015
    sigma_um: float = 25.0
    K_D_nM: float = 2.0
    s0_nM_per_s: float = 1.0
    half_width_um: float = 500.0
    candidate_dx_um: float = 10.0


EQUATIONS = {
    "ligand": "dL/dt = D d2L/dx2 - k_clear L + s0 exp(-x^2/(2 sigma^2))",
    "steady": "0 = D d2L/dx2 - k_clear L + s0 exp(-x^2/(2 sigma^2))",
    "occupancy": "theta = L/(K_D + L)",
    "boundary": "dL/dx = 0 at both finite-domain boundaries",
}
REFERENCE_C0_NM = 8.0


def parameter_table() -> list[dict[str, Any]]:
    return [
        {
            "parameter": "D_um2_per_s",
            "value": 55.0,
            "unit": "um^2 s^-1",
            "role": "fast extracellular diffusion coefficient",
            "source_or_assumption": "VERIFIED: Harish et al. 2023, Fig. 3C",
        },
        {
            "parameter": "k_clear_per_s",
            "value": 0.0015,
            "unit": "s^-1",
            "role": "linear effective clearance in the first free-ligand model",
            "source_or_assumption": "ASSUMPTION frozen in PREREG.md",
        },
        {
            "parameter": "sigma_um",
            "value": 25.0,
            "unit": "um",
            "role": "Gaussian source width",
            "source_or_assumption": "ASSUMPTION frozen in PREREG.md",
        },
        {
            "parameter": "K_D_nM",
            "value": 2.0,
            "unit": "nM",
            "role": "fast-equilibrium receptor half-occupancy concentration",
            "source_or_assumption": "ASSUMPTION frozen in PREREG.md",
        },
        {
            "parameter": "s0_nM_per_s",
            "value": 1.0,
            "unit": "nM s^-1",
            "role": "source amplitude; profile is normalized for threshold distances",
            "source_or_assumption": "ASSUMPTION frozen in PREREG.md",
        },
        {
            "parameter": "half_width_um",
            "value": 500.0,
            "unit": "um",
            "role": "distance from source to each no-flux outer boundary",
            "source_or_assumption": "ASSUMPTION frozen in PREREG.md",
        },
        {
            "parameter": "candidate_dx_um",
            "value": 10.0,
            "unit": "um",
            "role": "candidate spatial reporting grid",
            "source_or_assumption": "FROZEN criterion in PREREG.md",
        },
        {
            "parameter": "D_slow_um2_per_s",
            "value": 4.0,
            "unit": "um^2 s^-1",
            "role": "slow component reported by the anchor source, not used here",
            "source_or_assumption": "VERIFIED: Harish et al. 2023, Fig. 3C; excluded first model",
        },
        {
            "parameter": "slow_fraction",
            "value": 0.07,
            "unit": "dimensionless",
            "role": "slow fraction reported by the anchor source, not used here",
            "source_or_assumption": "VERIFIED: Harish et al. 2023, Fig. 3C; excluded first model",
        },
        {
            "parameter": "C0_reference_nM",
            "value": 8.0,
            "unit": "nM",
            "role": "source concentration reported by the anchor source, not imposed on model",
            "source_or_assumption": "VERIFIED: Harish et al. 2023, Fig. 4C",
        },
    ]


def unit_check(p: Parameters) -> dict[str, Any]:
    positive = {
        "D": p.D_um2_per_s,
        "k_clear": p.k_clear_per_s,
        "sigma": p.sigma_um,
        "K_D": p.K_D_nM,
        "s0": p.s0_nM_per_s,
        "half_width": p.half_width_um,
        "candidate_dx": p.candidate_dx_um,
    }
    all_positive = all(np.isfinite(value) and value > 0 for value in positive.values())
    diffusion_rate = p.D_um2_per_s / (p.candidate_dx_um**2)
    return {
        "status": "OK" if all_positive else "FAIL",
        "all_positive": bool(all_positive),
        "diffusion_rate_at_candidate_grid_per_s": float(diffusion_rate),
        "clearance_rate_per_s": float(p.k_clear_per_s),
        "rate_terms_have_same_unit": "s^-1" ,
        "source_term_unit": "nM s^-1",
        "concentration_unit": "nM",
        "coordinate_unit": "um",
    }


def _grid(p: Parameters, dx_um: float) -> tuple[np.ndarray, float]:
    if not np.isfinite(dx_um) or dx_um <= 0:
        raise ValueError("dx_um must be positive and finite")
    count = int(round(2.0 * p.half_width_um / dx_um))
    if count < 8 or count % 2 != 0:
        raise ValueError("2*half_width/dx must produce an even cell count of at least 8")
    dx = 2.0 * p.half_width_um / count
    x = -p.half_width_um + (np.arange(count) + 0.5) * dx
    return x, dx


def steady_profile(p: Parameters, dx_um: float | None = None) -> dict[str, np.ndarray | float]:
    selected_dx = p.candidate_dx_um if dx_um is None else dx_um
    x, dx = _grid(p, selected_dx)
    count = x.size
    source = p.s0_nM_per_s * np.exp(-0.5 * (x / p.sigma_um) ** 2)
    diffusion_rate = p.D_um2_per_s / (dx**2)
    diagonal = np.full(count, p.k_clear_per_s + 2.0 * diffusion_rate)
    diagonal[0] = p.k_clear_per_s + diffusion_rate
    diagonal[-1] = p.k_clear_per_s + diffusion_rate
    upper = np.full(count - 1, -diffusion_rate)
    lower = np.full(count - 1, -diffusion_rate)
    bands = np.zeros((3, count))
    bands[0, 1:] = upper
    bands[1] = diagonal
    bands[2, :-1] = lower
    concentration = solve_banded((1, 1), bands, source, check_finite=True)
    if np.any(concentration < -1e-10) or not np.all(np.isfinite(concentration)):
        raise FloatingPointError("steady solution is non-finite or negative")
    concentration = np.maximum(concentration, 0.0)
    source_index = int(np.argmin(np.abs(x)))
    source_value = float(concentration[source_index])
    scale = REFERENCE_C0_NM / source_value if source_value > 0.0 else 1.0
    reference_concentration = concentration * scale
    return {
        "x_um": x,
        "L_nM": concentration,
        "L_reference_nM": reference_concentration,
        "source_scale_to_C0": scale,
        "dx_um": float(dx),
    }


def _crossing(x: np.ndarray, values: np.ndarray, level: float) -> float:
    normalized = values / values[0]
    for index in range(1, x.size):
        if normalized[index] <= level:
            x0 = x[index - 1]
            x1 = x[index]
            y0 = normalized[index - 1]
            y1 = normalized[index]
            if y1 == y0:
                return float(x1)
            fraction = (level - y0) / (y1 - y0)
            return float(x0 + fraction * (x1 - x0))
    return float("inf")


def receptor_occupancy(concentration_nM: np.ndarray | float, K_D_nM: float) -> np.ndarray | float:
    values = np.asarray(concentration_nM, dtype=float)
    occupancy = values / (K_D_nM + values)
    if np.ndim(concentration_nM) == 0:
        return float(occupancy)
    return occupancy


def threshold_distances(profile: dict[str, np.ndarray | float], p: Parameters) -> dict[str, float]:
    x = np.asarray(profile["x_um"], dtype=float)
    values = np.asarray(profile.get("L_reference_nM", profile["L_nM"]), dtype=float)
    midpoint = int(np.argmin(np.abs(x)))
    x_positive = x[midpoint:]
    values_positive = values[midpoint:]
    x80 = _crossing(x_positive, values_positive, 0.8)
    x40 = _crossing(x_positive, values_positive, 0.4)
    occupancy = np.asarray(receptor_occupancy(values_positive, p.K_D_nM), dtype=float)
    x_half = _crossing(x_positive, values_positive / (p.K_D_nM + values_positive), 0.5)
    return {
        "x_0.8_um": x80,
        "x_0.4_um": x40,
        "threshold_span_um": x40 - x80,
        "x_KD_um": x_half,
        "occupancy_at_source": float(occupancy[0]),
        "occupancy_at_100um": float(
            values_positive[np.argmin(np.abs(x_positive - 100.0))]
            / (p.K_D_nM + values_positive[np.argmin(np.abs(x_positive - 100.0))])
        ),
    }


def solve_metrics(p: Parameters, dx_um: float | None = None) -> dict[str, Any]:
    profile = steady_profile(p, dx_um)
    metrics = threshold_distances(profile, p)
    x = np.asarray(profile["x_um"], dtype=float)
    raw_values = np.asarray(profile["L_nM"], dtype=float)
    reference_values = np.asarray(profile.get("L_reference_nM", raw_values), dtype=float)
    source_index = int(np.argmin(np.abs(x)))
    hundred_index = int(np.argmin(np.abs(x - 100.0)))
    metrics["source_L_nM"] = float(reference_values[source_index])
    metrics["raw_source_L_before_C0_scaling_nM"] = float(raw_values[source_index])
    metrics["L_at_0_um"] = float(reference_values[source_index])
    metrics["L_at_100um"] = float(reference_values[hundred_index])
    return {"metrics": metrics, "profile": profile}


def sensitivity_analysis(p: Parameters) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, field in (("D_um2_per_s", "D_um2_per_s"), ("k_clear_per_s", "k_clear_per_s"), ("K_D_nM", "K_D_nM")):
        for factor in (0.5, 1.0, 1.5):
            changed = replace(p, **{field: getattr(p, field) * factor})
            result = solve_metrics(changed, p.candidate_dx_um)["metrics"]
            rows.append(
                {
                    "parameter": name,
                    "factor": factor,
                    "D_um2_per_s": changed.D_um2_per_s,
                    "k_clear_per_s": changed.k_clear_per_s,
                    "K_D_nM": changed.K_D_nM,
                    "x_0.8_um": result["x_0.8_um"],
                    "x_0.4_um": result["x_0.4_um"],
                    "x_KD_um": result["x_KD_um"],
                    "threshold_span_um": result["threshold_span_um"],
                }
            )
    return rows


def _relative_change(coarse: float, fine: float) -> float:
    if not np.isfinite(coarse) or not np.isfinite(fine) or coarse == 0:
        return float("inf")
    return abs(fine - coarse) / abs(fine)


def run(p: Parameters | None = None) -> dict[str, Any]:
    selected = Parameters() if p is None else p
    coarse = solve_metrics(selected, selected.candidate_dx_um)["metrics"]
    fine = solve_metrics(selected, selected.candidate_dx_um / 2.0)["metrics"]
    dx80_change = _relative_change(coarse["x_0.8_um"], fine["x_0.8_um"])
    dx40_change = _relative_change(coarse["x_0.4_um"], fine["x_0.4_um"])
    span_grid_cells = coarse["threshold_span_um"] / selected.candidate_dx_um
    reference = {
        "doi": "10.1242/dev.201559",
        "source": "Harish et al., Development 150, dev201559 (2023)",
        "D_fast_um2_per_s": 55.0,
        "D_slow_um2_per_s": 4.0,
        "fast_fraction": 0.93,
        "slow_fraction": 0.07,
        "C0_nM": 8.0,
        "x_0.8_window_um": [40.0, 70.0],
        "x_0.4_window_um": [70.0, 140.0],
        "figure": "Fig. 3C and Fig. 4C",
    }
    criteria = {
        "x_0.8_in_40_70_um": bool(40.0 <= coarse["x_0.8_um"] <= 70.0),
        "x_0.4_in_70_140_um_diagnostic": bool(70.0 <= coarse["x_0.4_um"] <= 140.0),
        "grid_x_0.8_change_le_5_percent": bool(dx80_change <= 0.05),
        "grid_x_0.4_change_le_5_percent": bool(dx40_change <= 0.05),
        "threshold_span_at_least_4_grid_cells": bool(span_grid_cells >= 4.0),
        "zero_source_null_is_zero": bool(np.all(steady_profile(replace(selected, s0_nM_per_s=0.0))["L_nM"] == 0.0)),
    }
    return {
        "id": "BT-HX-Q026",
        "model": "one-dimensional free-ligand reaction-diffusion with equilibrium receptor readout",
        "equations": EQUATIONS,
        "parameters": asdict(selected),
        "parameter_table": parameter_table(),
        "unit_check": unit_check(selected),
        "reference": reference,
        "coarse_dx_um": selected.candidate_dx_um,
        "fine_dx_um": selected.candidate_dx_um / 2.0,
        "coarse_metrics": coarse,
        "fine_metrics": fine,
        "grid_convergence": {
            "relative_change_x_0.8": dx80_change,
            "relative_change_x_0.4": dx40_change,
            "threshold_span_grid_cells": span_grid_cells,
        },
        "criteria": criteria,
        "sensitivity": sensitivity_analysis(selected),
        "next_resolution_step": "Resolve L(x,t) and receptor occupancy on a grid no coarser than 10 um in the first pass, then add the measured slow/HSPG-bound pool and time-resolved source/receptor data before interpreting cell fates.",
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def main() -> None:
    output = run()
    Path("results.json").write_text(json.dumps(_json_safe(output), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
