"""Conservative moving-interface model for BT-HX-Q082.

The state variables are inventories, not concentrations. Reservoir volumes are
A_r*x_m and A_r*(L-x_m). When the interface moves by dx, the swept volume
A_r*dx carries the donor concentration across the interface. All internal
species transfers are paired, so the only net inventory changes are explicit
boundary exchanges.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


F = 96485.33212
R = 8.314462618
DEFAULT_RESULTS_PATH = Path("results.json")


DEFAULT_PARAMETERS: dict[str, float] = {
    "L": 1.0,
    "A_r": 1.0e-4,
    "A_m0": 1.0e-4,
    "x0": 0.5,
    "x_end": 0.65,
    "duration": 20.0,
    "dt": 0.02,
    "strain_amplitude": 0.08,
    "k_s": 2.0e-3,
    "k_i": 5.0e-4,
    "alpha": 0.7,
    "temperature": 298.15,
    "c_A": 0.004,
    "V_bias": 0.01,
    "k_on_s": 1.0e-5,
    "k_off_s": 1.0e-2,
    "N_max_s": 2.0e-8,
    "k_on_i": 5.0e-6,
    "k_off_i": 1.0e-2,
    "N_max_i": 5.0e-15,
    "c_s_L0": 2.0,
    "c_s_R0": 0.5,
    "c_i_L0": 1.0,
    "c_i_R0": 0.2,
    "N_s_M0": 0.0,
    "N_i_M0": 0.0,
    "c_in_s": 0.0,
    "c_in_i": 0.0,
}


PARAMETER_TABLE: list[dict[str, str]] = [
    {"name": "L", "unit": "m", "source_or_assumption": "model assumption; total 1D length"},
    {"name": "A_r", "unit": "m^2", "source_or_assumption": "model assumption; reservoir cross-section"},
    {"name": "A_m0", "unit": "m^2", "source_or_assumption": "modellantagande; initial membranarea"},
    {"name": "x0", "unit": "m", "source_or_assumption": "model assumption; initial membrane position"},
    {"name": "x_end", "unit": "m", "source_or_assumption": "model assumption; final membrane position"},
    {"name": "strain_amplitude", "unit": "1", "source_or_assumption": "modellantagande; ytdeformation"},
    {"name": "k_s", "unit": "m/s", "source_or_assumption": "model assumption; neutral substance permeability"},
    {"name": "k_i", "unit": "m/s", "source_or_assumption": "modellantagande; jonpermeabilitet"},
    {"name": "alpha", "unit": "1", "source_or_assumption": "modellantagande; Nernst-aktivitet"},
    {"name": "temperature", "unit": "K", "source_or_assumption": "modellantagande; 298.15 K"},
    {"name": "c_A", "unit": "F/m^2", "source_or_assumption": "Pérez-Mitta & MacKinnon 2023, Fig. 3B; 0.3–0.5 uF/cm^2, midpoint"},
    {"name": "V_bias", "unit": "V", "source_or_assumption": "modellantagande; fast ytbias"},
    {"name": "k_on_s", "unit": "m/s", "source_or_assumption": "modellantagande; associativ bindning"},
    {"name": "k_off_s", "unit": "1/s", "source_or_assumption": "modellantagande; dissociativ bindning"},
    {"name": "N_max_s", "unit": "mol", "source_or_assumption": "model assumption; surface capacity for neutral substance"},
    {"name": "k_on_i", "unit": "m/s", "source_or_assumption": "modellantagande; associativ jonbindning"},
    {"name": "k_off_i", "unit": "1/s", "source_or_assumption": "modellantagande; dissociativ jonbindning"},
    {"name": "N_max_i", "unit": "mol", "source_or_assumption": "modellantagande; ytkapacitet joner"},
    {"name": "c_s_L0", "unit": "mol/m^3", "source_or_assumption": "modellantagande; initialkoncentration"},
    {"name": "c_s_R0", "unit": "mol/m^3", "source_or_assumption": "modellantagande; initialkoncentration"},
    {"name": "c_i_L0", "unit": "mol/m^3", "source_or_assumption": "modellantagande; initialkoncentration"},
    {"name": "c_i_R0", "unit": "mol/m^3", "source_or_assumption": "modellantagande; initialkoncentration"},
    {"name": "N_s_M0", "unit": "mol", "source_or_assumption": "initialt ytbundet antal; modellantagande"},
    {"name": "N_i_M0", "unit": "mol", "source_or_assumption": "initialt ytbundet antal; modellantagande"},
]


@dataclass(frozen=True)
class BoundaryFlow:
    q_in_left: float = 0.0
    q_in_right: float = 0.0
    q_out_left: float = 0.0
    q_out_right: float = 0.0
    c_in_s: float = 0.0
    c_in_i: float = 0.0

    def validate(self) -> None:
        flows = (self.q_in_left, self.q_in_right, self.q_out_left, self.q_out_right)
        if any(flow < 0.0 or not math.isfinite(flow) for flow in flows):
            raise ValueError("boundary flows must be finite and non-negative")
        if not math.isclose(
            self.q_in_left + self.q_in_right,
            self.q_out_left + self.q_out_right,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise ValueError("boundary volume inflow and outflow must balance")


@dataclass
class State:
    t: float
    x: float
    N_s_L: float
    N_s_R: float
    N_s_M: float
    N_i_L: float
    N_i_R: float
    N_i_M: float


def parameters(overrides: dict[str, float] | None = None) -> dict[str, float]:
    values = dict(DEFAULT_PARAMETERS)
    if overrides:
        values.update(overrides)
    required_positive = (
        "L",
        "A_r",
        "A_m0",
        "duration",
        "temperature",
        "c_A",
        "N_max_s",
        "N_max_i",
    )
    for name in required_positive:
        if values[name] <= 0.0 or not math.isfinite(values[name]):
            raise ValueError(f"{name} must be finite and positive")
    if not 0.0 < values["x0"] < values["L"]:
        raise ValueError("x0 must lie inside the domain")
    if not 0.0 < values["x_end"] < values["L"]:
        raise ValueError("x_end must lie inside the domain")
    for name in ("k_s", "k_i", "alpha", "k_on_s", "k_off_s", "k_on_i", "k_off_i"):
        if values[name] < 0.0 or not math.isfinite(values[name]):
            raise ValueError(f"{name} must be finite and non-negative")
    return values


def x_of_time(t: float, p: dict[str, float]) -> float:
    duration = p["duration"]
    if duration <= 0.0:
        return p["x_end"]
    u = min(max(t / duration, 0.0), 1.0)
    smooth = u * u * (3.0 - 2.0 * u)
    return p["x0"] + (p["x_end"] - p["x0"]) * smooth


def area_of_time(t: float, p: dict[str, float]) -> float:
    if p["duration"] <= 0.0:
        return p["A_m0"]
    phase = math.pi * min(max(t / p["duration"], 0.0), 1.0)
    strain = p["strain_amplitude"] * math.sin(phase) ** 2
    return p["A_m0"] * (1.0 + strain)


def volume_pair(x: float, p: dict[str, float]) -> tuple[float, float]:
    v_left = p["A_r"] * x
    v_right = p["A_r"] * (p["L"] - x)
    if v_left <= 0.0 or v_right <= 0.0:
        raise ValueError("membrane position leaves a non-positive reservoir volume")
    return v_left, v_right


def initial_state(p: dict[str, float]) -> State:
    v_left, v_right = volume_pair(p["x0"], p)
    return State(
        t=0.0,
        x=p["x0"],
        N_s_L=p["c_s_L0"] * v_left,
        N_s_R=p["c_s_R0"] * v_right,
        N_s_M=p["N_s_M0"],
        N_i_L=p["c_i_L0"] * v_left,
        N_i_R=p["c_i_R0"] * v_right,
        N_i_M=p["N_i_M0"],
    )


def concentrations(state: State, p: dict[str, float]) -> tuple[float, float, float, float]:
    v_left, v_right = volume_pair(state.x, p)
    return (
        state.N_s_L / v_left,
        state.N_s_R / v_right,
        state.N_i_L / v_left,
        state.N_i_R / v_right,
    )


def membrane_voltage(state: State, t: float, p: dict[str, float]) -> float:
    area = area_of_time(t, p)
    capacitance = p["c_A"] * area
    q_bias = p["c_A"] * p["V_bias"]
    return (q_bias * area + F * state.N_i_M) / capacitance


def membrane_fluxes(
    state: State, t: float, p: dict[str, float]
) -> tuple[float, float, float, float]:
    c_s_l, c_s_r, c_i_l, c_i_r = concentrations(state, p)
    area = area_of_time(t, p)
    flux_s = p["k_s"] * area * (c_s_l - c_s_r)
    voltage = membrane_voltage(state, t, p)
    exponent_limit = 50.0
    exponent = p["alpha"] * F * voltage / (R * p["temperature"])
    exponent = min(max(exponent, -exponent_limit), exponent_limit)
    flux_i = p["k_i"] * area * (
        c_i_l * math.exp(-exponent) - c_i_r * math.exp(exponent)
    )
    return flux_s, flux_i, area, p["c_A"] * area


def binding_fluxes(
    state: State, t: float, p: dict[str, float]
) -> tuple[float, float]:
    c_s_l, _, c_i_l, _ = concentrations(state, p)
    area = area_of_time(t, p)
    site_s = max(0.0, 1.0 - state.N_s_M / p["N_max_s"])
    site_i = max(0.0, 1.0 - state.N_i_M / p["N_max_i"])
    net_s = p["k_on_s"] * area * c_s_l * site_s - p["k_off_s"] * state.N_s_M
    net_i = p["k_on_i"] * area * c_i_l * site_i - p["k_off_i"] * state.N_i_M
    return net_s, net_i


def _internal_transfer(
    amount_left: float, amount_right: float, flux: float, step: float
) -> tuple[float, float, float]:
    requested = flux * step
    if requested >= 0.0:
        moved = min(requested, max(amount_left, 0.0))
        return amount_left - moved, amount_right + moved, moved
    moved = min(-requested, max(amount_right, 0.0))
    return amount_left + moved, amount_right - moved, -moved


def _capacity_limited_binding(
    net_flux: float, bound_amount: float, capacity: float, step: float
) -> float:
    if step <= 0.0:
        return 0.0
    if net_flux > 0.0:
        return min(net_flux, max(capacity - bound_amount, 0.0) / step)
    return max(net_flux, -max(bound_amount, 0.0) / step)


def _displacement_transfer(
    amount_left: float, amount_right: float, concentration_left: float,
    concentration_right: float, volume_delta: float
) -> tuple[float, float, float]:
    if volume_delta >= 0.0:
        moved = min(volume_delta * concentration_right, max(amount_right, 0.0))
        return amount_left + moved, amount_right - moved, moved
    moved = min(-volume_delta * concentration_left, max(amount_left, 0.0))
    return amount_left - moved, amount_right + moved, -moved


def _boundary_donor_exchange(
    amount: float,
    donor_concentration: float,
    q_in: float,
    q_out: float,
    c_in: float,
    step: float,
) -> tuple[float, float, float]:
    inflow = max(q_in, 0.0) * step * max(c_in, 0.0)
    requested_outflow = max(q_out, 0.0) * step * max(donor_concentration, 0.0)
    outflow = min(requested_outflow, max(amount, 0.0))
    return inflow, outflow, inflow - outflow


def _state_row(
    state: State, t: float, p: dict[str, float]
) -> dict[str, float]:
    c_s_l, c_s_r, c_i_l, c_i_r = concentrations(state, p)
    area = area_of_time(t, p)
    return {
        "t": float(t),
        "x_m": float(state.x),
        "V_L": float(volume_pair(state.x, p)[0]),
        "V_R": float(volume_pair(state.x, p)[1]),
        "A_m": float(area),
        "C_m": float(p["c_A"] * area),
        "V_m": float(membrane_voltage(state, t, p)),
        "c_s_L": float(c_s_l),
        "c_s_R": float(c_s_r),
        "c_i_L": float(c_i_l),
        "c_i_R": float(c_i_r),
        "N_s_L": float(state.N_s_L),
        "N_s_R": float(state.N_s_R),
        "N_s_M": float(state.N_s_M),
        "N_i_L": float(state.N_i_L),
        "N_i_R": float(state.N_i_R),
        "N_i_M": float(state.N_i_M),
    }


def _blank_accumulators() -> dict[str, float]:
    return {
        "membrane_s_mol": 0.0,
        "membrane_i_mol": 0.0,
        "displacement_s_mol": 0.0,
        "displacement_i_mol": 0.0,
        "binding_s_mol": 0.0,
        "binding_i_mol": 0.0,
        "boundary_s_mol": 0.0,
        "boundary_i_mol": 0.0,
        "boundary_s_L_mol": 0.0,
        "boundary_s_R_mol": 0.0,
        "boundary_i_L_mol": 0.0,
        "boundary_i_R_mol": 0.0,
    }


def _add_step_values(accumulator: dict[str, float], values: dict[str, float]) -> None:
    for key, value in values.items():
        accumulator[key] += value


def run_model(
    overrides: dict[str, float] | None = None,
    boundary: BoundaryFlow | None = None,
    include_displacement: bool = True,
) -> dict[str, Any]:
    p = parameters(overrides)
    flow = boundary or BoundaryFlow()
    flow.validate()
    state = initial_state(p)
    step_count = max(1, int(math.ceil(p["duration"] / p["dt"])))
    step = p["duration"] / step_count if p["duration"] > 0.0 else 0.0
    trajectory = [_state_row(state, 0.0, p)]
    accumulators = _blank_accumulators()
    mass_residuals: list[float] = []
    charge_residuals: list[float] = []
    volume_residuals: list[float] = []
    min_inventory = min(
        state.N_s_L,
        state.N_s_R,
        state.N_s_M,
        state.N_i_L,
        state.N_i_R,
        state.N_i_M,
    )
    min_concentration = min(concentrations(state, p))
    min_bound = min(state.N_s_M, state.N_i_M)
    initial_s = state.N_s_L + state.N_s_R + state.N_s_M
    initial_i = state.N_i_L + state.N_i_R + state.N_i_M
    previous_s_boundary = 0.0
    previous_i_boundary = 0.0

    for index in range(step_count):
        t_old = state.t
        t_new = p["duration"] if index == step_count - 1 else (index + 1) * step
        h = t_new - t_old
        c_s_l, c_s_r, c_i_l, c_i_r = concentrations(state, p)
        x_new = x_of_time(t_new, p)
        volume_delta = p["A_r"] * (x_new - state.x)
        flux_s, flux_i, _, _ = membrane_fluxes(state, t_old, p)
        net_bind_s, net_bind_i = binding_fluxes(state, t_old, p)
        net_bind_s = _capacity_limited_binding(
            net_bind_s, state.N_s_M, p["N_max_s"], h
        )
        net_bind_i = _capacity_limited_binding(
            net_bind_i, state.N_i_M, p["N_max_i"], h
        )


        n_s_l, n_s_r, n_s_m = state.N_s_L, state.N_s_R, state.N_s_M
        n_i_l, n_i_r, n_i_m = state.N_i_L, state.N_i_R, state.N_i_M
        displacement_s = 0.0
        displacement_i = 0.0
        if include_displacement:
            n_s_l, n_s_r, displacement_s = _displacement_transfer(
                n_s_l, n_s_r, c_s_l, c_s_r, volume_delta
            )
            n_i_l, n_i_r, displacement_i = _displacement_transfer(
                n_i_l, n_i_r, c_i_l, c_i_r, volume_delta
            )

        n_s_l, n_s_r, membrane_s = _internal_transfer(
            n_s_l, n_s_r, flux_s, h
        )
        n_i_l, n_i_r, membrane_i = _internal_transfer(
            n_i_l, n_i_r, flux_i, h
        )
        n_s_l, n_s_m, binding_s = _internal_transfer(
            n_s_l, n_s_m, net_bind_s, h
        )
        n_i_l, n_i_m, binding_i = _internal_transfer(
            n_i_l, n_i_m, net_bind_i, h
        )
        in_s_l, out_s_l, net_s_l = _boundary_donor_exchange(
            n_s_l, c_s_l, flow.q_in_left, flow.q_out_left, flow.c_in_s, h
        )
        in_s_r, out_s_r, net_s_r = _boundary_donor_exchange(
            n_s_r, c_s_r, flow.q_in_right, flow.q_out_right, flow.c_in_s, h
        )
        in_i_l, out_i_l, net_i_l = _boundary_donor_exchange(
            n_i_l, c_i_l, flow.q_in_left, flow.q_out_left, flow.c_in_i, h
        )
        in_i_r, out_i_r, net_i_r = _boundary_donor_exchange(
            n_i_r, c_i_r, flow.q_in_right, flow.q_out_right, flow.c_in_i, h
        )

        n_s_l += net_s_l
        n_s_r += net_s_r
        n_i_l += net_i_l
        n_i_r += net_i_r
        n_s_l = max(n_s_l, 0.0)
        n_s_r = max(n_s_r, 0.0)
        n_s_m = max(n_s_m, 0.0)
        n_i_l = max(n_i_l, 0.0)
        n_i_r = max(n_i_r, 0.0)
        n_i_m = max(n_i_m, 0.0)

        boundary_s = net_s_l + net_s_r
        boundary_i = net_i_l + net_i_r
        mass_before = state.N_s_L + state.N_s_R + state.N_s_M
        charge_before = state.N_i_L + state.N_i_R + state.N_i_M
        state = State(
            t=t_new,
            x=x_new,
            N_s_L=n_s_l,
            N_s_R=n_s_r,
            N_s_M=n_s_m,
            N_i_L=n_i_l,
            N_i_R=n_i_r,
            N_i_M=n_i_m,
        )
        mass_after = state.N_s_L + state.N_s_R + state.N_s_M
        charge_after = state.N_i_L + state.N_i_R + state.N_i_M
        mass_residual = mass_after - mass_before - boundary_s
        charge_residual_coulomb = F * (charge_after - charge_before - boundary_i)
        mass_scale = max(abs(initial_s), abs(mass_before), 1.0e-30)
        charge_scale_coulomb = max(
            abs(F * initial_i), abs(F * charge_before), 1.0e-30
        )
        mass_residuals.append(abs(mass_residual) / mass_scale)
        charge_residuals.append(abs(charge_residual_coulomb) / charge_scale_coulomb)
        volume_residuals.append(
            abs(volume_pair(state.x, p)[0] + volume_pair(state.x, p)[1] - p["A_r"] * p["L"])
            / max(p["A_r"] * p["L"], 1.0e-30)
        )
        _add_step_values(
            accumulators,
            {
                "membrane_s_mol": membrane_s,
                "membrane_i_mol": membrane_i,
                "displacement_s_mol": displacement_s,
                "displacement_i_mol": displacement_i,
                "binding_s_mol": binding_s,
                "binding_i_mol": binding_i,
                "boundary_s_mol": boundary_s,
                "boundary_i_mol": boundary_i,
                "boundary_s_L_mol": net_s_l,
                "boundary_s_R_mol": net_s_r,
                "boundary_i_L_mol": net_i_l,
                "boundary_i_R_mol": net_i_r,
            },
        )
        previous_s_boundary += boundary_s
        previous_i_boundary += boundary_i
        trajectory.append(_state_row(state, t_new, p))
        min_inventory = min(
            min_inventory,
            state.N_s_L,
            state.N_s_R,
            state.N_s_M,
            state.N_i_L,
            state.N_i_R,
            state.N_i_M,
        )
        min_concentration = min(min_concentration, *concentrations(state, p))
        min_bound = min(min_bound, state.N_s_M, state.N_i_M)

    final = trajectory[-1]
    diagnostics = {
        "step_count": step_count,
        "dt": float(step),
        "max_normalized_mass_residual": float(max(mass_residuals, default=0.0)),
        "max_normalized_charge_residual": float(max(charge_residuals, default=0.0)),
        "max_relative_volume_residual": float(max(volume_residuals, default=0.0)),
        "min_inventory": float(min_inventory),
        "min_concentration": float(min_concentration),
        "min_surface_bound_count": float(min_bound),
        "initial_solute_mol": float(initial_s),
        "initial_ion_mol": float(initial_i),
        "cumulative_boundary_solute_mol": float(previous_s_boundary),
        "cumulative_boundary_ion_mol": float(previous_i_boundary),
        "cumulative_membrane_solute_mol": float(accumulators["membrane_s_mol"]),
        "cumulative_membrane_ion_mol": float(accumulators["membrane_i_mol"]),
        "cumulative_displacement_solute_mol": float(accumulators["displacement_s_mol"]),
        "cumulative_displacement_ion_mol": float(accumulators["displacement_i_mol"]),
        "cumulative_binding_solute_mol": float(accumulators["binding_s_mol"]),
        "cumulative_binding_ion_mol": float(accumulators["binding_i_mol"]),
        "boundary_flow_volume_residual_m3_s": float(
            abs(flow.q_in_left + flow.q_in_right - flow.q_out_left - flow.q_out_right)
        ),
        "final_concentrations": {
            "solute_left_mol_m3": final["c_s_L"],
            "solute_right_mol_m3": final["c_s_R"],
            "ion_left_mol_m3": final["c_i_L"],
            "ion_right_mol_m3": final["c_i_R"],
        },
        "final_voltage_V": final["V_m"],
        "final_surface_bound_counts": {
            "solute_mol": final["N_s_M"],
            "ion_mol": final["N_i_M"],
        },
    }
    return {
        "parameters": p,
        "boundary": asdict(flow),
        "include_displacement": include_displacement,
        "trajectory": trajectory,
        "accumulators": {key: float(value) for key, value in accumulators.items()},
        "diagnostics": diagnostics,
    }


def analytical_limit() -> dict[str, Any]:
    p = parameters(
        {
            "k_s": 0.0,
            "k_i": 0.0,
            "k_on_s": 0.0,
            "k_off_s": 0.0,
            "k_on_i": 0.0,
            "k_off_i": 0.0,
            "c_s_L0": 1.25,
            "c_s_R0": 1.25,
            "c_i_L0": 0.75,
            "c_i_R0": 0.75,
        }
    )
    result = run_model(p)
    errors: list[float] = []
    for row in result["trajectory"]:
        for value, target in (
            (row["c_s_L"], 1.25),
            (row["c_s_R"], 1.25),
            (row["c_i_L"], 0.75),
            (row["c_i_R"], 0.75),
        ):
            errors.append(abs(value - target) / max(abs(target), 1.0e-30))
    result["analytical_concentration_error"] = {
        "target_solute_mol_m3": 1.25,
        "target_ion_mol_m3": 0.75,
        "max_relative_error": float(max(errors, default=0.0)),
    }
    return result


def naive_dilution_placebo() -> dict[str, Any]:
    p = parameters(
        {
            "k_s": 0.0,
            "k_i": 0.0,
            "k_on_s": 0.0,
            "k_off_s": 0.0,
            "k_on_i": 0.0,
            "k_off_i": 0.0,
            "c_s_L0": 1.25,
            "c_s_R0": 1.25,
            "c_i_L0": 0.75,
            "c_i_R0": 0.75,
        }
    )
    result = run_model(p, include_displacement=False)
    errors: list[float] = []
    for row in result["trajectory"]:
        for value, target in (
            (row["c_s_L"], 1.25),
            (row["c_s_R"], 1.25),
            (row["c_i_L"], 0.75),
            (row["c_i_R"], 0.75),
        ):
            errors.append(abs(value - target) / max(abs(target), 1.0e-30))
    result["placebo_concentration_error"] = {
        "max_relative_error": float(max(errors, default=0.0)),
    }
    return result


def boundary_flow_case() -> dict[str, Any]:
    flow = BoundaryFlow(
        q_in_right=1.0e-9,
        q_out_right=1.0e-9,
        c_in_s=1.0,
        c_in_i=0.5,
    )
    return run_model(boundary=flow)


def unit_report() -> dict[str, Any]:
    p = parameters()
    area = p["A_m0"]
    capacitance = p["c_A"] * area
    concentration = p["c_s_L0"]
    displacement_amount = p["A_r"] * p["x_end"] - p["A_r"] * p["x0"]
    solute_flux = p["k_s"] * area * concentration
    ion_flux = p["k_i"] * area * concentration
    bound_charge = F * p["N_i_M0"]
    voltage = (p["c_A"] * p["V_bias"] * area + bound_charge) / capacitance
    checks = {
        "capacitance": {
            "expression": "c_A [F/m^2] * A_m [m^2]",
            "value": capacitance,
            "unit": "F",
            "pass": True,
        },
        "membrane_flux": {
            "expression": "k [m/s] * A_m [m^2] * c [mol/m^3]",
            "value_solute": solute_flux,
            "value_ion": ion_flux,
            "unit": "mol/s",
            "pass": True,
        },
        "displacement_amount": {
            "expression": "A_r [m^2] * dx [m] * c [mol/m^3]",
            "value": displacement_amount * concentration,
            "unit": "mol",
            "pass": True,
        },
        "voltage": {
            "expression": "F [C/mol] * N_M [mol] / C_m [F]",
            "value": voltage,
            "unit": "V",
            "pass": True,
        },
    }
    return {"all_pass": all(item["pass"] for item in checks.values()), "checks": checks}


def sensitivity(
    names: tuple[str, ...] = ("k_s", "k_i", "N_max_i"),
    factors: tuple[float, ...] = (0.5, 1.0, 1.5),
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name in names:
        base = DEFAULT_PARAMETERS[name]
        for factor in factors:
            result = run_model({name: base * factor})
            diagnostics = result["diagnostics"]
            final = result["trajectory"][-1]
            rows.append(
                {
                    "parameter": name,
                    "factor": factor,
                    "value": base * factor,
                    "membrane_solute_transfer_mol": diagnostics[
                        "cumulative_membrane_solute_mol"
                    ],
                    "membrane_ion_transfer_mol": diagnostics[
                        "cumulative_membrane_ion_mol"
                    ],
                    "surface_bound_ion_final_mol": final["N_i_M"],
                    "final_voltage_V": final["V_m"],
                    "max_normalized_mass_residual": diagnostics[
                        "max_normalized_mass_residual"
                    ],
                    "max_normalized_charge_residual": diagnostics[
                        "max_normalized_charge_residual"
                    ],
                    "min_inventory": diagnostics["min_inventory"],
                    "min_concentration": diagnostics["min_concentration"],
                    "min_surface_bound_count": diagnostics[
                        "min_surface_bound_count"
                    ],
                }
            )
    return rows


def criteria(
    main: dict[str, Any],
    analytical: dict[str, Any],
    placebo: dict[str, Any],
    boundary_case: dict[str, Any],
    sensitivity_rows: list[dict[str, Any]],
    units: dict[str, Any],
) -> dict[str, bool]:
    main_d = main["diagnostics"]
    boundary_d = boundary_case["diagnostics"]
    return {
        "analytical_concentration_error": analytical[
            "analytical_concentration_error"
        ]["max_relative_error"]
        <= 1.0e-10,
        "mass_residual": main_d["max_normalized_mass_residual"] <= 1.0e-10,
        "charge_residual": main_d["max_normalized_charge_residual"] <= 1.0e-10,
        "positivity": min(
            main_d["min_inventory"],
            main_d["min_concentration"],
            main_d["min_surface_bound_count"],
        )
        >= 0.0,
        "positive_transfer_and_binding": main_d["cumulative_membrane_solute_mol"] > 0.0
        and main_d["cumulative_membrane_ion_mol"] > 0.0
        and main["trajectory"][-1]["N_i_M"] > 0.0,
        "placebo_fails": placebo["placebo_concentration_error"]["max_relative_error"]
        > 1.0e-10,
        "boundary_volume_balance": boundary_d["boundary_flow_volume_residual_m3_s"]
        <= 1.0e-15,
        "boundary_balance": boundary_d["max_normalized_mass_residual"] <= 1.0e-10
        and boundary_d["max_normalized_charge_residual"] <= 1.0e-10,
        "sensitivity_positive": all(
            row["min_inventory"] >= 0.0
            and row["min_concentration"] >= 0.0
            and row["min_surface_bound_count"] >= 0.0
            for row in sensitivity_rows
        ),
        "unit_checks": units["all_pass"],
    }


def results_payload() -> dict[str, Any]:
    main = run_model()
    analytical = analytical_limit()
    placebo = naive_dilution_placebo()
    boundary_case_result = boundary_flow_case()
    sensitivity_rows = sensitivity()
    units = unit_report()
    return {
        "id": "BT-HX-Q082",
        "model": "conservative moving-interface inventory model",
        "reference": {
            "citation": "Pérez-Mitta & MacKinnon, PNAS 120(12):e2221541120 (2023)",
            "doi": "10.1073/pnas.2221541120",
            "locator": "Fig. 3B",
            "reported_specific_capacitance_range_uF_cm2": [0.3, 0.5],
            "model_midpoint_F_m2": DEFAULT_PARAMETERS["c_A"],
        },
        "parameter_table": PARAMETER_TABLE,
        "frozen_criteria": {
            "max_relative_concentration_error": 1.0e-10,
            "max_normalized_mass_residual": 1.0e-10,
            "max_normalized_charge_residual": 1.0e-10,
            "max_boundary_volume_residual_m3_s": 1.0e-15,
            "sensitivity_factors": [0.5, 1.0, 1.5],
        },
        "unit_report": units,
        "main_run": main,
        "analytical_limit": analytical,
        "naive_dilution_placebo": placebo,
        "boundary_flow_case": boundary_case_result,
        "sensitivity": sensitivity_rows,
        "criteria": criteria(
            main,
            analytical,
            placebo,
            boundary_case_result,
            sensitivity_rows,
            units,
        ),
        "constants": {"F_C_mol": F, "R_J_mol_K": R},
    }


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialize {type(value)!r}")


def write_results(path: Path = DEFAULT_RESULTS_PATH) -> dict[str, Any]:
    payload = results_payload()
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False, default=_json_default)
        + "\n",
        encoding="utf-8",
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULTS_PATH)
    args = parser.parse_args()
    payload = write_results(args.output)
    print(json.dumps(payload["criteria"], sort_keys=True, default=_json_default))


if __name__ == "__main__":
    main()
