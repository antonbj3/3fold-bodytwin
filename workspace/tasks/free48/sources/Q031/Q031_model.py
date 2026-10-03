from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class Parameters:
    F_C_mol: float = 96485.33212
    R_J_mol_K: float = 8.314462618
    temperature_K: float = 297.15
    water_molar_volume_m3_mol: float = 18.07e-6
    initial_volume_m3: float = 1.7e-15
    area_to_volume_initial_m_inv: float = 675000.0
    membrane_capacitance_F_m2: float = 1.0e-2
    hydraulic_permeability_m_s_Pa: float = 1.65e-11
    water_area_fraction: float = 0.025
    na_conductance_S_m2: float = 0.2
    k_conductance_S_m2: float = 1.0
    cl_conductance_S_m2: float = 0.5
    pump_density_mol_m2_s: float = 2.0e-7
    pump_na_half_saturation_mM: float = 10.0
    pump_k_half_saturation_mM: float = 4.0
    pump_baseline_activation: float = 0.25
    initial_na_mM: float = 10.0
    initial_k_mM: float = 140.0
    initial_cl_mM: float = 10.0
    initial_fixed_anion_mM: float = 140.0
    external_na_mM: float = 100.0
    external_k_mM: float = 4.0
    external_cl_mM: float = 196.0
    post_shock_scale: float = 0.5
    shock_time_s: float = 0.0
    t_end_s: float = 30.0


PARAMETER_TABLE = [
    {
        "name": "F_C_mol",
        "value": 96485.33212,
        "unit": "C mol^-1",
        "source_or_assumption": "CODATA 2018 Faraday constant",
    },
    {
        "name": "R_J_mol_K",
        "value": 8.314462618,
        "unit": "J mol^-1 K^-1",
        "source_or_assumption": "CODATA 2018 molar gas constant",
    },
    {
        "name": "water_molar_volume_m3_mol",
        "value": 18.07e-6,
        "unit": "m^3 mol^-1",
        "source_or_assumption": "Assumption: water molar volume at 24 degrees C",
    },
    {
        "name": "initial_volume_m3",
        "value": 1.7e-15,
        "unit": "m^3",
        "source_or_assumption": "Farinas & Verkman 1996, Results p. 3516: 1.7 pL",
    },
    {
        "name": "area_to_volume_initial_m_inv",
        "value": 675000.0,
        "unit": "m^-1",
        "source_or_assumption": "Farinas & Verkman 1996, Results p. 3516: 6750 cm^-1",
    },
    {
        "name": "hydraulic_permeability_m_s_Pa",
        "value": 1.65e-11,
        "unit": "m s^-1 Pa^-1",
        "source_or_assumption": "Derived from Farinas & Verkman 1996 Pf 6.1e-4 cm s^-1 over a 150 mM step; conversion assumes Pf is per osmolar concentration",
    },
    {
        "name": "water_area_fraction",
        "value": 0.025,
        "unit": "dimensionless",
        "source_or_assumption": "Assumption: effective hydraulically exposed fraction of the cell-layer geometry",
    },
    {
        "name": "temperature_K",
        "value": 297.15,
        "unit": "K",
        "source_or_assumption": "24 degrees C protocol temperature in Farinas & Verkman 1996",
    },
    {
        "name": "membrane_capacitance_F_m2",
        "value": 1.0e-2,
        "unit": "F m^-2",
        "source_or_assumption": "Assumption: 1 microF cm^-2 generic membrane value",
    },
    {
        "name": "na_conductance_S_m2",
        "value": 0.2,
        "unit": "S m^-2",
        "source_or_assumption": "Assumption: generic open Na conductance density",
    },
    {
        "name": "k_conductance_S_m2",
        "value": 1.0,
        "unit": "S m^-2",
        "source_or_assumption": "Assumption: generic open K conductance density",
    },
    {
        "name": "cl_conductance_S_m2",
        "value": 0.5,
        "unit": "S m^-2",
        "source_or_assumption": "Assumption: generic open Cl conductance density",
    },
    {
        "name": "pump_density_mol_m2_s",
        "value": 2.0e-7,
        "unit": "mol m^-2 s^-1",
        "source_or_assumption": "Assumption: pump density chosen to balance the initial generic currents",
    },
    {
        "name": "pump_na_half_saturation_mM",
        "value": 10.0,
        "unit": "mM",
        "source_or_assumption": "Assumption: simple pump activation scale",
    },
    {
        "name": "pump_k_half_saturation_mM",
        "value": 4.0,
        "unit": "mM",
        "source_or_assumption": "Assumption: simple pump activation scale",
    },
    {
        "name": "pump_baseline_activation",
        "value": 0.25,
        "unit": "dimensionless",
        "source_or_assumption": "Normalization so the initial generic pump activation equals one",
    },
    {
        "name": "initial_na_mM",
        "value": 10.0,
        "unit": "mM",
        "source_or_assumption": "Assumption chosen for a 300 mM electroneutral starting solution",
    },
    {
        "name": "initial_k_mM",
        "value": 140.0,
        "unit": "mM",
        "source_or_assumption": "Assumption chosen for a 300 mM electroneutral starting solution",
    },
    {
        "name": "initial_cl_mM",
        "value": 10.0,
        "unit": "mM",
        "source_or_assumption": "Assumption chosen for a 300 mM electroneutral starting solution",
    },
    {
        "name": "initial_fixed_anion_mM",
        "value": 140.0,
        "unit": "mM",
        "source_or_assumption": "Assumption chosen for electroneutrality and 300 mM osmolarity",
    },
    {
        "name": "external_na_mM",
        "value": 100.0,
        "unit": "mM",
        "source_or_assumption": "Assumption: baseline external ionic composition summing to 300 mM",
    },
    {
        "name": "external_k_mM",
        "value": 4.0,
        "unit": "mM",
        "source_or_assumption": "Assumption: baseline external ionic composition summing to 300 mM",
    },
    {
        "name": "external_cl_mM",
        "value": 196.0,
        "unit": "mM",
        "source_or_assumption": "Assumption: baseline external ionic composition summing to 300 mM",
    },
    {
        "name": "post_shock_scale",
        "value": 0.5,
        "unit": "dimensionless",
        "source_or_assumption": "Protocol control: 300 mM to 150 mM external osmolarity",
    },
    {
        "name": "shock_time_s",
        "value": 0.0,
        "unit": "s",
        "source_or_assumption": "Protocol definition: shock applied at integration start",
    },
    {
        "name": "t_end_s",
        "value": 30.0,
        "unit": "s",
        "source_or_assumption": "Frozen local simulation horizon",
    },
]

REFERENCE = {
    "doi": "10.1016/S0006-3495(96)79546-2",
    "citation": "Farinas J, Verkman AS, Biophysical Journal 71, 3511-3522 (1996)",
    "figure": "Figure 6 and Results pp. 3517-3518",
    "volume_ratio_dimensionless": 2.0,
    "half_time_s": 6.0,
    "water_permeability_cm_s": 6.1e-4,
    "water_permeability_uncertainty_cm_s": 2.0e-4,
    "area_to_volume_cm_inv": 6750.0,
    "initial_volume_pL": 1.7,
}

SPECIES = ("Na", "K", "Cl")
VALENCE = {"Na": 1, "K": 1, "Cl": -1}


def concentration_from_amount(amount_mol: float, volume_m3: float) -> float:
    return amount_mol / volume_m3


def amount_from_concentration(concentration_m3: float, volume_m3: float) -> float:
    return concentration_m3 * volume_m3


def initial_area(parameters: Parameters) -> float:
    return parameters.area_to_volume_initial_m_inv * parameters.initial_volume_m3


def scaled_area(volume_m3: float, parameters: Parameters) -> float:
    ratio = max(volume_m3 / parameters.initial_volume_m3, 1.0e-12)
    return initial_area(parameters) * ratio ** (2.0 / 3.0)


def nernst_potential(
    concentration_inside_m3: float,
    concentration_outside_m3: float,
    valence: int,
    parameters: Parameters,
) -> float:
    thermal_voltage = parameters.R_J_mol_K * parameters.temperature_K / parameters.F_C_mol
    return thermal_voltage / valence * np.log(concentration_outside_m3 / concentration_inside_m3)


def ionic_fluxes(
    concentration_inside_m3: dict[str, float],
    concentration_outside_m3: dict[str, float],
    membrane_potential_V: float,
    area_m2: float,
    parameters: Parameters,
) -> tuple[dict[str, float], dict[str, float]]:
    molar_flux_out_mol_s: dict[str, float] = {}
    electrical_current_A: dict[str, float] = {}
    conductances = {
        "Na": parameters.na_conductance_S_m2,
        "K": parameters.k_conductance_S_m2,
        "Cl": parameters.cl_conductance_S_m2,
    }
    for species in SPECIES:
        equilibrium_potential_V = nernst_potential(
            concentration_inside_m3[species],
            concentration_outside_m3[species],
            VALENCE[species],
            parameters,
        )
        flux_out_mol_s = (
            conductances[species]
            * area_m2
            / (VALENCE[species] * parameters.F_C_mol)
            * (membrane_potential_V - equilibrium_potential_V)
        )
        molar_flux_out_mol_s[species] = flux_out_mol_s
        electrical_current_A[species] = VALENCE[species] * parameters.F_C_mol * flux_out_mol_s
    return molar_flux_out_mol_s, electrical_current_A


def pump_flux(
    concentration_na_m3: float,
    concentration_k_out_m3: float,
    area_m2: float,
    parameters: Parameters,
) -> tuple[dict[str, float], float, float]:
    na_activation = concentration_na_m3 / (parameters.pump_na_half_saturation_mM + concentration_na_m3)
    k_activation = concentration_k_out_m3 / (parameters.pump_k_half_saturation_mM + concentration_k_out_m3)
    activation = na_activation * k_activation / parameters.pump_baseline_activation
    activation = float(np.clip(activation, 0.0, 1.5))
    cycle_flux_mol_s = parameters.pump_density_mol_m2_s * area_m2 * activation
    flux_out_mol_s = {"Na": 3.0 * cycle_flux_mol_s, "K": -2.0 * cycle_flux_mol_s, "Cl": 0.0}
    pump_current_A = parameters.F_C_mol * cycle_flux_mol_s
    return flux_out_mol_s, pump_current_A, activation


def initial_external_concentrations(parameters: Parameters, scale: float = 1.0) -> dict[str, float]:
    return {
        "Na": parameters.external_na_mM * scale,
        "K": parameters.external_k_mM * scale,
        "Cl": parameters.external_cl_mM * scale,
    }


def initial_internal_concentrations(parameters: Parameters) -> dict[str, float]:
    return {
        "Na": parameters.initial_na_mM,
        "K": parameters.initial_k_mM,
        "Cl": parameters.initial_cl_mM,
    }


def total_ionic_current(
    membrane_potential_V: float,
    concentration_inside_m3: dict[str, float],
    concentration_outside_m3: dict[str, float],
    area_m2: float,
    parameters: Parameters,
) -> tuple[float, dict[str, float], dict[str, float], float, float]:
    flux_ion, current_ion = ionic_fluxes(
        concentration_inside_m3,
        concentration_outside_m3,
        membrane_potential_V,
        area_m2,
        parameters,
    )
    flux_pump, pump_current_A, activation = pump_flux(
        concentration_inside_m3["Na"],
        concentration_outside_m3["K"],
        area_m2,
        parameters,
    )
    all_flux = {species: flux_ion[species] + flux_pump[species] for species in SPECIES}
    all_current = {species: current_ion[species] for species in SPECIES}
    total_A = sum(all_current.values()) + pump_current_A
    return total_A, all_flux, all_current, pump_current_A, activation


def resting_membrane_potential(parameters: Parameters) -> float:
    inside = initial_internal_concentrations(parameters)
    outside = initial_external_concentrations(parameters)
    area_m2 = initial_area(parameters)

    def residual(potential_V: float) -> float:
        return total_ionic_current(potential_V, inside, outside, area_m2, parameters)[0]

    return float(brentq(residual, -0.30, 0.30, xtol=1.0e-12))


def initial_state(parameters: Parameters) -> tuple[np.ndarray, float]:
    volume_m3 = parameters.initial_volume_m3
    inside = initial_internal_concentrations(parameters)
    amounts = {species: amount_from_concentration(inside[species], volume_m3) for species in SPECIES}
    fixed_anion_m3 = parameters.initial_fixed_anion_mM
    fixed_anion_amount_mol = amount_from_concentration(fixed_anion_m3, volume_m3)
    state = np.array(
        [
            resting_membrane_potential(parameters),
            volume_m3,
            amounts["Na"],
            amounts["K"],
            amounts["Cl"],
            volume_m3 / parameters.water_molar_volume_m3_mol,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )
    return state, fixed_anion_amount_mol


def external_scale_function(parameters: Parameters) -> Callable[[float], float]:
    def scale(time_s: float) -> float:
        if time_s < parameters.shock_time_s:
            return 1.0
        return parameters.post_shock_scale

    return scale


def state_concentrations(
    state: np.ndarray, parameters: Parameters, time_s: float = 0.0
) -> tuple[dict[str, float], dict[str, float]]:
    volume_m3 = state[1]
    inside = {
        "Na": concentration_from_amount(state[2], volume_m3),
        "K": concentration_from_amount(state[3], volume_m3),
        "Cl": concentration_from_amount(state[4], volume_m3),
    }
    return inside, outside_at_time(0.0, time_s, parameters)


def outside_at_time(_state_value: float, time_s: float, parameters: Parameters) -> dict[str, float]:
    return initial_external_concentrations(parameters, external_scale_function(parameters)(time_s))


def rhs(
    time_s: float,
    state: np.ndarray,
    parameters: Parameters,
    coupled_voltage: bool = True,
    fixed_voltage_V: float | None = None,
) -> np.ndarray:
    voltage_state_V = state[0]
    volume_m3 = max(state[1], 1.0e-18)
    inside = {
        "Na": concentration_from_amount(state[2], volume_m3),
        "K": concentration_from_amount(state[3], volume_m3),
        "Cl": concentration_from_amount(state[4], volume_m3),
    }
    outside = outside_at_time(0.0, time_s, parameters)
    area_m2 = scaled_area(volume_m3, parameters)
    voltage_for_flux_V = voltage_state_V if coupled_voltage else fixed_voltage_V
    if voltage_for_flux_V is None:
        raise ValueError("fixed_voltage_V is required when coupled_voltage is False")
    flux_ion, current_ion = ionic_fluxes(inside, outside, voltage_for_flux_V, area_m2, parameters)
    flux_pump, pump_current_A, _ = pump_flux(inside["Na"], outside["K"], area_m2, parameters)
    total_flux = {species: flux_ion[species] + flux_pump[species] for species in SPECIES}
    net_tracked_charge_flux_mol_s = total_flux["Na"] + total_flux["K"] - total_flux["Cl"]
    ionic_current_A = sum(current_ion.values())
    fixed_anion_concentration_m3 = (
        parameters.initial_fixed_anion_mM
        * parameters.initial_volume_m3
        / volume_m3
    )
    total_inside_osm_m3 = sum(inside.values()) + fixed_anion_concentration_m3 + abs(state[9]) / volume_m3
    total_outside_osm_m3 = sum(outside.values())
    osmotic_pressure_difference_Pa = (
        parameters.R_J_mol_K * parameters.temperature_K * (total_inside_osm_m3 - total_outside_osm_m3)
    )
    water_area_m2 = area_m2 * parameters.water_area_fraction
    water_flux_m3_s = parameters.hydraulic_permeability_m_s_Pa * water_area_m2 * osmotic_pressure_difference_Pa
    dvolume_dt = water_flux_m3_s
    capacitance_F = parameters.membrane_capacitance_F_m2 * area_m2
    dcapacitance_dt = (2.0 / 3.0) * capacitance_F * dvolume_dt / volume_m3
    dvoltage_dt = 0.0 if not coupled_voltage else -(pump_current_A + ionic_current_A + voltage_state_V * dcapacitance_dt) / capacitance_F
    dstate = np.zeros(10, dtype=float)
    dstate[0] = dvoltage_dt
    dstate[1] = dvolume_dt
    dstate[2] = -total_flux["Na"]
    dstate[3] = -total_flux["K"]
    dstate[4] = -total_flux["Cl"]
    dstate[5] = dvolume_dt / parameters.water_molar_volume_m3_mol
    dstate[6] = total_flux["Na"]
    dstate[7] = total_flux["K"]
    dstate[8] = total_flux["Cl"]
    dstate[9] = -net_tracked_charge_flux_mol_s
    return dstate


def solve_trajectory(
    parameters: Parameters = Parameters(),
    coupled_voltage: bool = True,
    time_end_s: float | None = None,
    points: int = 2001,
) -> dict[str, np.ndarray | float | bool]:
    state0, fixed_anion_m3 = initial_state(parameters)
    initial_voltage_V = state0[0]
    end_s = parameters.t_end_s if time_end_s is None else time_end_s
    time_grid_s = np.linspace(parameters.shock_time_s, end_s, points)
    absolute_atol = np.array(
        [1.0e-11, 1.0e-23, 1.0e-23, 1.0e-23, 1.0e-23, 1.0e-22, 1.0e-25, 1.0e-25, 1.0e-25, 1.0e-25],
        dtype=float,
    )

    def wrapped_rhs(time_s: float, state: np.ndarray) -> np.ndarray:
        return rhs(time_s, state, parameters, coupled_voltage, initial_voltage_V)

    solution = solve_ivp(
        wrapped_rhs,
        (time_grid_s[0], time_grid_s[-1]),
        state0,
        method="Radau",
        t_eval=time_grid_s,
        atol=absolute_atol,
        rtol=2.0e-9,
        max_step=0.005,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return {
        "time_s": solution.t,
        "state": solution.y,
        "fixed_anion_m3": fixed_anion_m3,
        "initial_voltage_V": initial_voltage_V,
        "success": solution.success,
    }


def _first_crossing(time_s: np.ndarray, values: np.ndarray, target: float) -> float:
    above = np.flatnonzero(values >= target)
    if len(above) == 0:
        return float("nan")
    index = int(above[0])
    if index == 0:
        return float(time_s[0])
    x0 = float(time_s[index - 1])
    x1 = float(time_s[index])
    y0 = float(values[index - 1])
    y1 = float(values[index])
    if y1 == y0:
        return x1
    return x0 + (target - y0) * (x1 - x0) / (y1 - y0)


def trajectory_metrics(
    full: dict[str, np.ndarray | float | bool],
    clamped: dict[str, np.ndarray | float | bool] | None = None,
    parameters: Parameters = Parameters(),
) -> dict[str, float | bool]:
    full_time = np.asarray(full["time_s"], dtype=float)
    full_state = np.asarray(full["state"], dtype=float)
    volume_ratio = full_state[1] / parameters.initial_volume_m3
    final_ratio = float(volume_ratio[-1])
    half_target = 0.5 * (1.0 + final_ratio)
    half_time_s = _first_crossing(full_time, volume_ratio, half_target)
    result: dict[str, float | bool] = {
        "initial_voltage_mV": float(full["initial_voltage_V"]) * 1.0e3,
        "final_volume_ratio": final_ratio,
        "peak_volume_ratio": float(np.max(volume_ratio)),
        "half_time_s": half_time_s,
    }
    if clamped is not None:
        clamped_state = np.asarray(clamped["state"], dtype=float)
        difference = np.abs(full_state[1] - clamped_state[1]) / parameters.initial_volume_m3
        result["chi_volume"] = float(np.max(difference))
        result["max_volume_difference_ratio"] = result["chi_volume"]
        result["clamped_final_volume_ratio"] = float(clamped_state[1, -1] / parameters.initial_volume_m3)
    else:
        result["chi_volume"] = float("nan")
        result["max_volume_difference_ratio"] = float("nan")
        result["clamped_final_volume_ratio"] = float("nan")
    return result


def trajectory_state_summary(
    trajectory: dict[str, np.ndarray | float | bool],
    parameters: Parameters = Parameters(),
) -> dict[str, object]:
    time_s = np.asarray(trajectory["time_s"], dtype=float)
    state = np.asarray(trajectory["state"], dtype=float)
    target_times_s = np.array([0.0, 0.01, 0.1, 1.0, 5.0, 10.0, float(time_s[-1])], dtype=float)
    target_times_s = target_times_s[(target_times_s >= time_s[0]) & (target_times_s <= time_s[-1])]
    indices = sorted(set(int(np.argmin(np.abs(time_s - target))) for target in target_times_s))
    samples: list[dict[str, float]] = []
    for index in indices:
        volume_m3 = float(state[1, index])
        fixed_anion_concentration_m3 = (
            parameters.initial_fixed_anion_mM * parameters.initial_volume_m3 / volume_m3
        )
        countercharge_concentration_m3 = float(state[9, index]) / volume_m3
        concentrations_m3 = {
            species: float(state[2 + species_index, index]) / volume_m3
            for species_index, species in enumerate(SPECIES)
        }
        samples.append(
            {
                "time_s": float(time_s[index]),
                "Vm_mV": float(state[0, index]) * 1.0e3,
                "volume_m3": volume_m3,
                "volume_ratio": volume_m3 / parameters.initial_volume_m3,
                "N_Na_mol": float(state[2, index]),
                "N_K_mol": float(state[3, index]),
                "N_Cl_mol": float(state[4, index]),
                "N_water_mol": float(state[5, index]),
                "N_countercharge_mol": float(state[9, index]),
                "C_Na_mM": concentrations_m3["Na"],
                "C_K_mM": concentrations_m3["K"],
                "C_Cl_mM": concentrations_m3["Cl"],
                "C_fixed_anion_mM": fixed_anion_concentration_m3,
                "C_countercharge_mM": countercharge_concentration_m3,
                "total_osmolarity_mM": sum(concentrations_m3.values())
                + fixed_anion_concentration_m3
                + abs(countercharge_concentration_m3),
                "net_bulk_charge_mM": concentrations_m3["Na"]
                + concentrations_m3["K"]
                - concentrations_m3["Cl"]
                - fixed_anion_concentration_m3
                - countercharge_concentration_m3,
            }
        )
    return {
        "samples": samples,
        "terminal": samples[-1],
        "max_abs_Vm_change_mV": float(np.max(np.abs(state[0] - state[0, 0])) * 1.0e3),
    }


def balance_metrics(
    trajectory: dict[str, np.ndarray | float | bool],
    parameters: Parameters = Parameters(),
) -> dict[str, float]:
    state = np.asarray(trajectory["state"], dtype=float)
    initial_state_array = state[:, 0]
    volume_m3 = np.maximum(state[1], 1.0e-18)
    fixed_anion_amount_mol = float(trajectory["fixed_anion_m3"])
    tracked_bulk_charge_C = parameters.F_C_mol * (
        state[2] + state[3] - state[4] - fixed_anion_amount_mol
    )
    integrated_charge_flux_C = parameters.F_C_mol * (state[6] + state[7] - state[8])
    charge_closure_C = tracked_bulk_charge_C + integrated_charge_flux_C
    charge_closure_C -= charge_closure_C[0]
    charge_scale_C = parameters.F_C_mol * (
        initial_state_array[2] + initial_state_array[3] + initial_state_array[4] + fixed_anion_amount_mol
    )
    species_balance_mol = np.array(
        [
            state[2] + state[6] - initial_state_array[2],
            state[3] + state[7] - initial_state_array[3],
            state[4] + state[8] - initial_state_array[4],
        ]
    )
    species_balance_rel = np.max(np.abs(species_balance_mol), axis=1) / np.array(
        [initial_state_array[2], initial_state_array[3], initial_state_array[4]]
    )
    water_reference_mol = initial_state_array[5] - initial_state_array[1] / parameters.water_molar_volume_m3_mol
    water_balance_mol = state[5] - state[1] / parameters.water_molar_volume_m3_mol - water_reference_mol
    water_balance_rel = np.max(np.abs(water_balance_mol)) / (
        parameters.initial_volume_m3 / parameters.water_molar_volume_m3_mol
    )
    total_bulk_charge_C = tracked_bulk_charge_C - parameters.F_C_mol * state[9]
    bulk_charge_density_m3 = total_bulk_charge_C / volume_m3
    cation_scale_m3 = (state[2] + state[3]) / volume_m3
    bulk_electroneutrality_rel = np.max(np.abs(bulk_charge_density_m3) / np.maximum(cation_scale_m3, 1.0e-30))
    countercharge_rel = np.max(np.abs(state[9])) / (
        initial_state_array[2] + initial_state_array[3] + initial_state_array[4] + fixed_anion_amount_mol
    )
    return {
        "max_species_balance_relative": float(np.max(species_balance_rel)),
        "max_water_balance_relative": float(water_balance_rel),
        "max_charge_closure_relative": float(np.max(np.abs(charge_closure_C)) / charge_scale_C),
        "max_bulk_electroneutrality_relative": float(bulk_electroneutrality_rel),
        "max_countercharge_relative": float(countercharge_rel),
        "max_charge_closure_C": float(np.max(np.abs(charge_closure_C))),
    }


def ideal_osmotic_final_volume(
    parameters: Parameters = Parameters(),
    post_shock_scale: float | None = None,
) -> float:
    scale = parameters.post_shock_scale if post_shock_scale is None else post_shock_scale
    inside_osm_m3 = (
        parameters.initial_na_mM + parameters.initial_k_mM + parameters.initial_cl_mM + parameters.initial_fixed_anion_mM
    )
    outside_osm_m3 = (
        parameters.external_na_mM + parameters.external_k_mM + parameters.external_cl_mM
    ) * scale
    return parameters.initial_volume_m3 * inside_osm_m3 / outside_osm_m3


def run_experiment(
    parameters: Parameters = Parameters(),
    time_end_s: float | None = None,
    points: int = 2001,
) -> dict[str, object]:
    full = solve_trajectory(parameters, coupled_voltage=True, time_end_s=time_end_s, points=points)
    clamped = solve_trajectory(parameters, coupled_voltage=False, time_end_s=time_end_s, points=points)
    metrics = trajectory_metrics(full, clamped, parameters)
    balances = balance_metrics(full, parameters)
    state_summary = trajectory_state_summary(full, parameters)
    reference_ratio = ideal_osmotic_final_volume(parameters)
    reference_pass = bool(
        1.8 <= float(metrics["final_volume_ratio"]) <= 2.2
        and np.isfinite(float(metrics["half_time_s"]))
        and 3.0 <= float(metrics["half_time_s"]) <= 10.0
    )
    coupling_required = bool(float(metrics["chi_volume"]) > 0.02)
    balance_pass = bool(
        float(balances["max_species_balance_relative"]) < 1.0e-6
        and float(balances["max_water_balance_relative"]) < 1.0e-6
        and float(balances["max_charge_closure_relative"]) < 1.0e-6
    )
    result: dict[str, object] = {
        "parameters": asdict(parameters),
        "reference": REFERENCE,
        "reference_ideal_final_volume_ratio": reference_ratio / parameters.initial_volume_m3,
        "metrics": metrics,
        "balances": balances,
        "trajectory_summary": state_summary,
        "reference_pass": reference_pass,
        "balance_pass": balance_pass,
        "coupling_required_chi_gt_0_02": coupling_required,
        "placebo_osmotic_final_volume_ratio": reference_ratio / parameters.initial_volume_m3,
    }
    return result


def sensitivity_analysis(
    parameters: Parameters = Parameters(),
    time_end_s: float | None = None,
    points: int = 1201,
) -> dict[str, dict[str, object]]:
    names = ("hydraulic_permeability_m_s_Pa", "k_conductance_S_m2", "pump_density_mol_m2_s")
    base = run_experiment(parameters, time_end_s=time_end_s, points=points)
    base_metrics = base["metrics"]
    output: dict[str, dict[str, object]] = {}
    for name in names:
        low = replace(parameters, **{name: getattr(parameters, name) * 0.5})
        high = replace(parameters, **{name: getattr(parameters, name) * 1.5})
        low_result = run_experiment(low, time_end_s=time_end_s, points=points)
        high_result = run_experiment(high, time_end_s=time_end_s, points=points)
        output[name] = {}
        for label, result in (("minus50", low_result), ("plus50", high_result)):
            output[name][label] = {
                "metrics": result["metrics"],
                "balances": result["balances"],
                "reference_pass": result["reference_pass"],
                "balance_pass": result["balance_pass"],
                "coupling_required_chi_gt_0_02": result["coupling_required_chi_gt_0_02"],
            }
        output[name]["sensitivity"] = {
            metric: float(
                (float(high_result["metrics"][metric]) - float(low_result["metrics"][metric]))
                / float(base_metrics[metric])
            )
            for metric in ("final_volume_ratio", "half_time_s", "chi_volume")
        }
    return output


def main() -> None:
    result = run_experiment()
    result["sensitivity"] = sensitivity_analysis()
    result["parameter_table"] = PARAMETER_TABLE
    result["species"] = list(SPECIES)
    result["valence"] = VALENCE
    result["balance_definition"] = {
        "species": "N_i(t) + integral(J_i dt) - N_i(0)",
        "water": "N_water(t) - V(t)/v_water - [N_water(0) - V(0)/v_water]",
        "charge": "F[(N_Na + N_K - N_Cl - N_X) + integral(J_Na + J_K - J_Cl) - (N_Na(0) + N_K(0) - N_Cl(0) - N_X(0))]",
        "countercharge": "N_delta keeps Q_bulk = F(N_Na + N_K - N_Cl - N_X - N_delta) at bulk electroneutrality",
    }
    result["prereg_sha256"] = Path("PREREG.sha256").read_text().split()[0]
    result["frozen_criteria"] = {
        "volume_ratio_interval": [1.8, 2.2],
        "half_time_interval_s": [3.0, 10.0],
        "chi_volume_threshold": 0.02,
        "species_water_charge_relative_threshold": 1.0e-6,
        "external_osmolarity_step_mOsm": [300.0, 150.0],
        "sensitivity_factors": [0.5, 1.5],
    }
    result["pump_stoichiometry"] = {"Na_out": 3, "K_in": 2}
    result["geometry_area_exponent"] = 2.0 / 3.0
    result["structural_source"] = {
        "doi": "10.1152/ajpcell.1980.238.5.c196",
        "citation": "Jakobsson E, American Journal of Physiology-Cell Physiology 238, C196-C206 (1980)",
    }
    Path("results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
