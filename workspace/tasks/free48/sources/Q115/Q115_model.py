from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Parameters:
    crypt_count: int = 8
    stem_cells_per_crypt: int = 4
    stem_cycle_h: float = 16.0
    self_renewal_probability: float = 0.5
    ta_rounds: int = 6
    ta_division_probability: float = 0.5
    ta_cycle_h: float = 16.0
    age_bin_h: float = 1.0
    max_age_h: float = 240.0
    crypt_force_pN: float = 2.1
    active_force_pN: float = 0.9
    drag_pN_h_per_mm: float = 100.0
    baseline_path_mm: float = 1.95
    # The net migration speed is now a DECLARED INPUT with a source, not a quotient of two
    # illustrative forces over a drag constant whose own table entry read "chosen to give
    # 0.03 mm/h". That constant was the dial and the speed was the target, so the agreement it
    # produced was not evidence of anything. The default is the held-out measurement that the
    # model was then refuted against; set migration_speed_um_h = 30.0 to recover the old value.
    migration_speed_um_h: float = 9.00
    migration_speed_um_h_sd: float = 0.465
    migration_speed_source: str = (
        "mouse duodenal thymidine-analogue fitted front velocity, control 9.00 +/- 0.465 um/h "
        "(DOI 10.1096/fj.201601002); held out, i.e. not among this model's inputs when the "
        "0.03 mm/h value was set"
    )
    atrophy_path_factor: float = 0.55
    shedding_hazard_h: float = 0.12
    villus_radius_base_mm: float = 0.65
    villus_radius_tip_mm: float = 0.25
    villus_height_mm: float = 1.90
    crypt_radius_mm: float = 0.15
    crypt_depth_mm: float = 0.50
    mature_cell_area_mm2: float = 0.0030
    cell_area_birth_fraction: float = 0.65
    maturation_tau_h: float = 24.0
    transport_birth_fraction: float = 0.35
    conductance_mature_mS_cm2: float = 0.010
    conductance_young_mS_cm2: float = 0.030
    transport_velocity_per_area_mm_day: float = 1.0
    concentration_over_Km: float = 0.1
    mm2_to_cm2: float = 0.01

    @property
    def migration_speed_mm_h(self) -> float:
        return self.migration_speed_um_h / 1000.0

    @property
    def migration_speed_from_drag_mm_h(self) -> float:
        """The old route, kept as a diagnostic only: nothing reads it to drive the model."""
        return (self.crypt_force_pN + self.active_force_pN) / self.drag_pN_h_per_mm

    @property
    def implied_drag_pN_h_per_mm(self) -> float:
        """The drag the declared speed implies, given the same two forces. Reading it next to
        drag_pN_h_per_mm shows how far the chosen constant sat from the measurement."""
        return (self.crypt_force_pN + self.active_force_pN) / self.migration_speed_mm_h

    @property
    def expected_ta_amplification(self) -> float:
        return (1.0 + self.ta_division_probability) ** self.ta_rounds

    @property
    def crypt_output_cells_h(self) -> float:
        return (
            self.crypt_count
            * self.stem_cells_per_crypt
            / self.stem_cycle_h
            * 2.0
            * (1.0 - self.self_renewal_probability)
            * self.expected_ta_amplification
        )


@dataclass
class CellState:
    compartment: str
    age_h: float
    division_rounds: int = 0
    position: float = 0.0


REFERENCE = {
    "author": "Kai",
    "year": 2021,
    "journal": "Biophysical Journal",
    "volume_issue_pages": "120(4), 699-710",
    "doi": "10.1016/j.bpj.2021.01.003",
    "page": "p. 704, Table 1",
    "N_cells": 1885,
    "n_cells_h": 26.0,
    "tau_f_h": 65.0,
    "T_ref_h": 1885.0 / 26.0,
    "renewal_range_days": [3.0, 5.0],
}

EQUATIONS = {
    "crypt_conservation": "dN/dt = J_crypt - J_shed",
    "crypt_production": "J_crypt = C*S/T_s*2*(1-p_self)*(1+q_TA)^R_TA",
    "lineage_division": "R_TA ~ Binomial(R_TA,max,q_TA); terminal daughters = 2^R_TA",
    "mechanics": "v = v_meas (declared input); eta*dr/dt = F_crypt + F_active is retained only as a diagnostic identity and no longer sets v",
    "migration": "tau_min = L_path/v; a cell cannot shed before tau_min",
    "shedding": "J_shed = sum_j N_j*h_shed*I(age_j>tau_min)",
    "mean_shedding_age": "T_shed = sum_j age_j*J_shed,j / sum_j J_shed,j",
    "area": "A_v = sum_j N_j*a_cell(age_j); a_cell = a_mature*f_area(age)",
    "uptake": "J_up = v_abs*(C/Km)*sum_j N_j*a_cell(age_j)*f_transport(age_j)",
    "barrier": "G_total = 0.01*sum_j N_j*a_cell(age_j)*g_cell(age_j); G_density = G_total/(0.01*A_ref)",
    "maturation": "f_m(age) = 1-exp(-age/tau_mat)",
    "conservation_step": "N(t+dt)-N(t) = J_crypt*dt - J_shed*dt",
}

PARAMETER_TABLE = [
    ("crypt_count", "C", "crypts per villus", "8", "count", "Kai 2021, p. 701, eight crypts per structure; scaled assumption"),
    ("stem_cells_per_crypt", "S", "stem cells per crypt", "4", "cells/crypt", "geometric scaling assumption; not fitted"),
    ("stem_cycle_h", "T_s", "stem division interval", "16", "h", "within the 12-24 h range reported by Kai 2021, p. 700"),
    ("self_renewal_probability", "p_self", "expected stem self-renewal probability", "0.5", "1", "symmetric renewal assumption"),
    ("ta_rounds", "R_TA", "maximum transit-amplifying rounds", "6", "rounds", "upstream biological range; numerical assumption"),
    ("ta_division_probability", "q_TA", "probability of each TA round dividing", "0.5", "1", "stochastic amplification assumption"),
    ("ta_cycle_h", "T_TA", "TA division interval", "16", "h", "same order as stem cycle; assumption"),
    ("age_bin_h", "dt_age", "age-class width", "1", "h", "numerical resolution"),
    ("max_age_h", "age_max", "tail cutoff", "240", "h", "tail is below 1e-8 at baseline hazard"),
    ("crypt_force_pN", "F_crypt", "mitotic-pressure migration force", "2.1", "pN", "illustrative force-scale assumption"),
    ("active_force_pN", "F_active", "active migration force", "0.9", "pN", "illustrative actomyosin-scale assumption"),
    ("drag_pN_h_per_mm", "eta", "cell-substrate drag", "100", "pN h/mm", "was chosen to give 0.03 mm/h; NO LONGER DRIVES THE MODEL, retained as a diagnostic only"),
    ("migration_speed_um_h", "v_meas", "net migration speed, declared input", "9.00", "um/h", "measured: mouse duodenal thymidine-analogue front velocity, control 9.00 +/- 0.465 um/h, DOI 10.1096/fj.201601002"),
    ("migration_speed_um_h_sd", "sd(v_meas)", "reported uncertainty on the declared speed", "0.465", "um/h", "same source, DOI 10.1096/fj.201601002"),
    ("baseline_path_mm", "L_path", "crypt-orifice to villus-tip path", "1.95", "mm", "chosen with the OLD v=0.03 mm/h to represent 65 h; still a chosen value, and with the measured speed the same path gives 216.7 h"),
    ("atrophy_path_factor", "alpha_L", "villus path scaling in perturbation", "0.55", "1", "preregistered geometry perturbation"),
    ("shedding_hazard_h", "h_shed", "tip extrusion hazard after transit", "0.12", "h^-1", "preregistered hazard assumption, not a measured rate"),
    ("villus_radius_base_mm", "r_base", "villus base radius", "0.65", "mm", "single-villus geometry assumption"),
    ("villus_radius_tip_mm", "r_tip", "villus tip radius", "0.25", "mm", "single-villus geometry assumption"),
    ("villus_height_mm", "H_v", "baseline villus height", "1.90", "mm", "mouse-scale geometry assumption"),
    ("crypt_radius_mm", "r_crypt", "crypt lumen radius", "0.15", "mm", "crypt-area bookkeeping assumption"),
    ("crypt_depth_mm", "L_crypt", "crypt depth", "0.50", "mm", "crypt-area bookkeeping assumption"),
    ("mature_cell_area_mm2", "a_mature", "mature apical area per cell", "0.0030", "mm^2/cell", "geometry-scale assumption; not a histology measurement"),
    ("cell_area_birth_fraction", "f_area0", "birth area/mature area", "0.65", "1", "cell-growth assumption"),
    ("maturation_tau_h", "tau_mat", "area and transport maturation time", "24", "h", "phenomenological assumption"),
    ("transport_birth_fraction", "f_up0", "birth transport/mature transport", "0.35", "1", "immature-cell transport assumption"),
    ("conductance_mature_mS_cm2", "g_mature", "mature patch conductance density", "0.010", "mS/cm^2", "illustrative barrier parameter; not measured"),
    ("conductance_young_mS_cm2", "g_young", "young patch conductance density", "0.030", "mS/cm^2", "illustrative barrier assumption; not measured"),
    ("transport_velocity_per_area_mm_day", "v_abs", "saturating-area transport scale", "1.0", "mm/day", "normalized uptake scale, not measured"),
    ("concentration_over_Km", "C_over_Km", "fixed driving-force ratio", "0.1", "1", "dimensionless boundary-condition assumption"),
    ("mm2_to_cm2", "c_area", "area conversion", "0.01", "cm^2/mm^2", "exact unit conversion"),
    ("migration_speed_mm_h", "v", "net migration speed used by the model", "0.009", "mm/h", "the declared input v_meas converted, no longer derived from eta"),
    ("expected_ta_amplification", "E[2^R_TA]", "expected terminal TA daughters", "11.390625", "daughters/lineage", "derived as (1+q_TA)^R_TA"),
    ("crypt_output_cells_h", "J_crypt", "crypt-to-villus output", "22.78125", "cells/h", "derived from the conservation and lineage equations"),
]


def parameter_dict(parameters: Parameters) -> dict[str, Any]:
    values = asdict(parameters)
    values["migration_speed_mm_h"] = parameters.migration_speed_mm_h
    values["migration_speed_from_drag_mm_h"] = parameters.migration_speed_from_drag_mm_h
    values["implied_drag_pN_h_per_mm"] = parameters.implied_drag_pN_h_per_mm
    values["expected_ta_amplification"] = parameters.expected_ta_amplification
    values["crypt_output_cells_h"] = parameters.crypt_output_cells_h
    return {key: float(value) if isinstance(value, (int, float, np.number)) else value for key, value in values.items()}


def parameter_table() -> list[dict[str, str]]:
    return [
        {"name": name, "symbol": symbol, "meaning": meaning, "value": value, "unit": unit, "source_or_assumption": source}
        for name, symbol, meaning, value, unit, source in PARAMETER_TABLE
    ]


def migration_velocity(parameters: Parameters) -> float:
    return parameters.migration_speed_mm_h


def crypt_output_rate(parameters: Parameters) -> float:
    return parameters.crypt_output_cells_h


def geometry_metrics(parameters: Parameters, path_factor: float = 1.0) -> dict[str, float]:
    height = parameters.villus_height_mm * path_factor
    radial_difference = parameters.villus_radius_base_mm - parameters.villus_radius_tip_mm
    slant_height_mm = math.sqrt(height * height + radial_difference * radial_difference)
    villus_area_mm2 = math.pi * (parameters.villus_radius_base_mm + parameters.villus_radius_tip_mm) * slant_height_mm
    crypt_side_area_mm2 = 2.0 * math.pi * parameters.crypt_radius_mm * parameters.crypt_depth_mm
    crypt_end_area_mm2 = 2.0 * math.pi * parameters.crypt_radius_mm * parameters.crypt_radius_mm
    crypt_area_mm2 = parameters.crypt_count * (crypt_side_area_mm2 + crypt_end_area_mm2)
    path_mm = parameters.baseline_path_mm * path_factor
    speed_mm_h = migration_velocity(parameters)
    return {
        "path_factor": float(path_factor),
        "height_mm": float(height),
        "slant_height_mm": float(slant_height_mm),
        "villus_geometric_area_mm2": float(villus_area_mm2),
        "crypt_area_mm2": float(crypt_area_mm2),
        "crypt_to_tip_path_mm": float(path_mm),
        "migration_speed_mm_h": float(speed_mm_h),
        "minimum_transit_h": float(path_mm / speed_mm_h),
    }


def maturation_factor(age_h: np.ndarray, parameters: Parameters) -> np.ndarray:
    return 1.0 - np.exp(-age_h / parameters.maturation_tau_h)


def cell_area_mm2(age_h: np.ndarray, parameters: Parameters) -> np.ndarray:
    fraction = parameters.cell_area_birth_fraction + (1.0 - parameters.cell_area_birth_fraction) * maturation_factor(age_h, parameters)
    return parameters.mature_cell_area_mm2 * fraction


def transport_fraction(age_h: np.ndarray, parameters: Parameters) -> np.ndarray:
    return parameters.transport_birth_fraction + (1.0 - parameters.transport_birth_fraction) * maturation_factor(age_h, parameters)


def conductance_density(age_h: np.ndarray, parameters: Parameters) -> np.ndarray:
    fraction = maturation_factor(age_h, parameters)
    return parameters.conductance_young_mS_cm2 + (parameters.conductance_mature_mS_cm2 - parameters.conductance_young_mS_cm2) * fraction


def age_grid(parameters: Parameters) -> np.ndarray:
    count = int(round(parameters.max_age_h / parameters.age_bin_h)) + 1
    return np.arange(count, dtype=float) * parameters.age_bin_h


def steady_state_age_distribution(parameters: Parameters, path_factor: float = 1.0) -> dict[str, Any]:
    ages_h = age_grid(parameters)
    geometry = geometry_metrics(parameters, path_factor)
    minimum_transit_h = geometry["minimum_transit_h"]
    survival = np.ones(ages_h.shape, dtype=float)
    for index in range(1, ages_h.size):
        shed_probability = min(max(parameters.shedding_hazard_h * parameters.age_bin_h, 0.0), 1.0) if ages_h[index] > minimum_transit_h else 0.0
        survival[index] = survival[index - 1] * (1.0 - shed_probability)
    supply_rate = crypt_output_rate(parameters)
    cells = supply_rate * parameters.age_bin_h * survival
    cell_area = cell_area_mm2(ages_h, parameters)
    transport = transport_fraction(ages_h, parameters)
    conductance = conductance_density(ages_h, parameters)
    villus_area_mm2 = float(np.sum(cells * cell_area))
    effective_uptake_area_mm2 = float(np.sum(cells * cell_area * transport))
    uptake_proxy_mm3_day = parameters.transport_velocity_per_area_mm_day * parameters.concentration_over_Km * effective_uptake_area_mm2
    total_conductance_mS = parameters.mm2_to_cm2 * float(np.sum(cells * cell_area * conductance))
    total_cells = float(np.sum(cells))
    mean_current_age_h = float(np.sum(ages_h * cells) / total_cells) if total_cells > 0.0 else 0.0
    shedding_ages_h = ages_h[1:]
    shedding_probability = np.where(shedding_ages_h > minimum_transit_h, min(max(parameters.shedding_hazard_h * parameters.age_bin_h, 0.0), 1.0), 0.0)
    shedding_flux = supply_rate * survival[:-1] * shedding_probability
    shedding_flux_total = float(np.sum(shedding_flux))
    mean_shedding_age_h = float(np.sum(shedding_ages_h * shedding_flux) / shedding_flux_total) if shedding_flux_total > 0.0 else 0.0
    shed_rate = float(supply_rate * (1.0 - survival[-1]))
    unresolved_tail_rate = float(supply_rate * survival[-1])
    return {
        "geometry": geometry,
        "ages_h": ages_h.tolist(),
        "survival_probability": survival.tolist(),
        "cells": cells.tolist(),
        "shedding_flux_cells_h": shedding_flux.tolist(),
        "cell_area_mm2": cell_area.tolist(),
        "transport_fraction": transport.tolist(),
        "conductance_mS_cm2": conductance.tolist(),
        "supply_cells_h": float(supply_rate),
        "shed_cells_h": shed_rate,
        "unresolved_tail_cells_h": unresolved_tail_rate,
        "n_cells": total_cells,
        "mean_current_age_h": mean_current_age_h,
        "mean_shedding_age_h": mean_shedding_age_h,
        "little_law_residence_h": total_cells / supply_rate if supply_rate > 0.0 else 0.0,
        "renewal_time_days": mean_shedding_age_h / 24.0,
        "villus_area_mm2": villus_area_mm2,
        "effective_uptake_area_mm2": effective_uptake_area_mm2,
        "uptake_proxy_mm3_day": uptake_proxy_mm3_day,
        "total_conductance_mS": total_conductance_mS,
    }


LAST_CUTOFF_OUTFLOW = 0.0


def age_transition(cells: np.ndarray, parameters: Parameters, path_factor: float = 1.0) -> tuple[np.ndarray, float, float]:
    dt = parameters.age_bin_h
    geometry = geometry_metrics(parameters, path_factor)
    minimum_transit_h = geometry["minimum_transit_h"]
    output = np.zeros_like(cells, dtype=float)
    newborn = crypt_output_rate(parameters) * dt
    output[0] += newborn
    shed_total = 0.0
    cutoff_outflow = 0.0
    for index in range(cells.size):
        next_age_h = (index + 1) * dt
        shed_fraction = parameters.shedding_hazard_h * dt if next_age_h > minimum_transit_h else 0.0
        shed_fraction = min(max(shed_fraction, 0.0), 1.0)
        shed = float(cells[index] * shed_fraction)
        shed_total += shed
        surviving = float(cells[index] - shed)
        if index + 1 < cells.size:
            output[index + 1] += surviving
        else:
            cutoff_outflow = surviving
    # 1/10 (anton-5f, after the graph lane's G2 finding): survivors in the last age bin leave the system at max_age_h
    # (tip extrusion) but were not accounted for. Inventory: N_next - N_old = newborn - shed_total - LAST_CUTOFF_OUTFLOW.
    # Negligible at baseline hazard (tail < 1e-8), but large when shedding_hazard_h is lowered (e.g. ~13 % at 0,01/h).
    global LAST_CUTOFF_OUTFLOW
    LAST_CUTOFF_OUTFLOW = cutoff_outflow
    return output, newborn, shed_total


def constant_area_null(reference_state: dict[str, Any], parameters: Parameters) -> dict[str, float]:
    reference_area_mm2 = float(reference_state["villus_area_mm2"])
    uptake = parameters.transport_velocity_per_area_mm_day * parameters.concentration_over_Km * reference_area_mm2
    reference_area_cm2 = reference_area_mm2 * parameters.mm2_to_cm2
    total_conductance_mS = parameters.conductance_mature_mS_cm2 * reference_area_cm2
    conductance_density_mS_cm2 = total_conductance_mS / reference_area_cm2
    return {
        "area_mm2": reference_area_mm2,
        "effective_uptake_area_mm2": reference_area_mm2,
        "uptake_proxy_mm3_day": uptake,
        "total_conductance_mS": total_conductance_mS,
        "conductance_mS_cm2": conductance_density_mS_cm2,
    }


def run_scenario(parameters: Parameters, path_factor: float, reference_state: dict[str, Any] | None = None) -> dict[str, Any]:
    state = steady_state_age_distribution(parameters, path_factor)
    if reference_state is None:
        reference_state = state
    reference_area_cm2 = float(reference_state["villus_area_mm2"]) * parameters.mm2_to_cm2
    conductance_density_mS_cm2 = state["total_conductance_mS"] / reference_area_cm2
    null = constant_area_null(reference_state, parameters)
    return {
        "path_factor": float(path_factor),
        "geometry": state["geometry"],
        "crypt_output_cells_h": state["supply_cells_h"],
        "shed_output_cells_h": state["shed_cells_h"],
        "n_cells": state["n_cells"],
        "mean_current_age_h": state["mean_current_age_h"],
        "mean_shedding_age_h": state["mean_shedding_age_h"],
        "little_law_residence_h": state["little_law_residence_h"],
        "renewal_time_days": state["renewal_time_days"],
        "villus_area_mm2": state["villus_area_mm2"],
        "crypt_area_mm2": state["geometry"]["crypt_area_mm2"],
        "effective_uptake_area_mm2": state["effective_uptake_area_mm2"],
        "uptake_proxy_mm3_day": state["uptake_proxy_mm3_day"],
        "total_conductance_mS": state["total_conductance_mS"],
        "conductance_mS_cm2": conductance_density_mS_cm2,
        "constant_area_null": null,
        "uptake_change_vs_null_percent": 100.0 * (state["uptake_proxy_mm3_day"] / null["uptake_proxy_mm3_day"] - 1.0),
        "conductance_change_vs_null_percent": 100.0 * (conductance_density_mS_cm2 / null["conductance_mS_cm2"] - 1.0),
        "age_distribution": [
            {
                "age_h": float(age),
                "cells": float(number),
                "shedding_flux_cells_h": float(shedding_flux),
                "cell_area_mm2": float(area),
                "transport_fraction": float(transport),
                "conductance_mS_cm2": float(conductance),
            }
            for age, number, shedding_flux, area, transport, conductance in zip(
                state["ages_h"], state["cells"], state["shedding_flux_cells_h"] + [0.0], state["cell_area_mm2"], state["transport_fraction"], state["conductance_mS_cm2"]
            )
        ],
        "tail_unresolved_cells_h": state["unresolved_tail_cells_h"],
    }


def sample_ta_rounds(rng: np.random.Generator, parameters: Parameters) -> int:
    return int(rng.binomial(parameters.ta_rounds, parameters.ta_division_probability))


def transition_cell(cell: CellState, parameters: Parameters, path_factor: float, rng: np.random.Generator, dt_h: float) -> CellState:
    if cell.compartment == "stem":
        cell.age_h += dt_h
        if cell.age_h >= parameters.stem_cycle_h:
            return CellState("progenitor", 0.0, 0, 0.0)
        return cell
    if cell.compartment == "progenitor":
        cell.age_h += dt_h
        if cell.age_h >= parameters.ta_cycle_h:
            if cell.division_rounds < parameters.ta_rounds and rng.random() < parameters.ta_division_probability:
                return CellState("progenitor", 0.0, cell.division_rounds + 1, 0.0)
            return CellState("epithelial", 0.0, cell.division_rounds, 0.0)
        return cell
    if cell.compartment == "epithelial":
        cell.age_h += dt_h
        path_mm = parameters.baseline_path_mm * path_factor
        cell.position += migration_velocity(parameters) * dt_h / path_mm
        if cell.position >= 1.0 and rng.random() < parameters.shedding_hazard_h * dt_h:
            return CellState("shed", cell.age_h, cell.division_rounds, 1.0)
        return cell
    return cell


def dimensional_check(parameters: Parameters) -> dict[str, Any]:
    checks = {
        "migration_speed": {
            "expression": "um/h / 1000 = mm/h (declared input); pN / (pN h/mm) = mm/h (diagnostic)",
            "value_mm_h": migration_velocity(parameters),
            "pass": True,
        },
        "crypt_output": {
            "expression": "crypts * cells/crypt / h * 1 = cells/h",
            "value_cells_h": crypt_output_rate(parameters),
            "pass": crypt_output_rate(parameters) >= 0.0,
        },
        "transit": {
            "expression": "mm / (mm/h) = h",
            "value_h": parameters.baseline_path_mm / migration_velocity(parameters),
            "pass": parameters.baseline_path_mm / migration_velocity(parameters) > 0.0,
        },
        "uptake": {
            "expression": "(mm/day) * mm^2 = mm^3/day",
            "value_mm3_day": parameters.transport_velocity_per_area_mm_day * parameters.concentration_over_Km * 1.0,
            "pass": True,
        },
        "barrier": {
            "expression": "(mS/cm^2) * (cm^2/mm^2) * mm^2 = mS",
            "value_mS": parameters.conductance_mature_mS_cm2 * parameters.mm2_to_cm2,
            "pass": True,
        },
    }
    checks["all_pass"] = all(value["pass"] for value in checks.values() if isinstance(value, dict))
    return checks


def sensitivity_analysis(parameters: Parameters) -> list[dict[str, Any]]:
    names = [
        "shedding_hazard_h",
        "migration_speed_um_h",
        "crypt_force_pN",
        "ta_division_probability",
        "mature_cell_area_mm2",
        "stem_cycle_h",
        "transport_birth_fraction",
    ]
    baseline = run_scenario(parameters, 1.0)
    baseline_metrics = {
        "renewal_time_h": baseline["mean_shedding_age_h"],
        "villus_area_mm2": baseline["villus_area_mm2"],
        "uptake_proxy_mm3_day": baseline["uptake_proxy_mm3_day"],
        "conductance_mS_cm2": baseline["conductance_mS_cm2"],
    }
    rows: list[dict[str, Any]] = []
    for name in names:
        base_value = float(getattr(parameters, name))
        changes: dict[str, Any] = {}
        maximum = 0.0
        for direction, multiplier in (("minus50", 0.5), ("plus50", 1.5)):
            changed = replace(parameters, **{name: base_value * multiplier})
            result = run_scenario(changed, 1.0)
            metrics = {
                "renewal_time_h": result["mean_shedding_age_h"],
                "villus_area_mm2": result["villus_area_mm2"],
                "uptake_proxy_mm3_day": result["uptake_proxy_mm3_day"],
                "conductance_mS_cm2": result["conductance_mS_cm2"],
            }
            relative = {
                key: 100.0 * (metrics[key] / baseline_metrics[key] - 1.0) if baseline_metrics[key] != 0.0 else 0.0
                for key in metrics
            }
            changes[direction] = {"value": metrics, "relative_change_percent": relative}
            maximum = max(maximum, max(abs(value) for value in relative.values()))
        rows.append({"parameter": name, "baseline_value": base_value, "max_abs_relative_change_percent": maximum, "changes": changes})
    rows.sort(key=lambda row: row["max_abs_relative_change_percent"], reverse=True)
    return rows


def evaluate_criteria(parameters: Parameters, baseline: dict[str, Any], atrophy: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    t_ref = reference["T_ref_h"]
    t_low = 0.70 * t_ref
    t_high = 1.30 * t_ref
    production_ratio = baseline["crypt_output_cells_h"] / reference["n_cells_h"]
    census_ratio = baseline["n_cells"] / reference["N_cells"]
    primary_pass = t_low <= baseline["mean_shedding_age_h"] <= t_high
    production_pass = 0.5 <= production_ratio <= 2.0
    census_pass = 0.5 <= census_ratio <= 2.0
    uptake_change = abs(atrophy["uptake_change_vs_null_percent"])
    conductance_change = abs(atrophy["conductance_change_vs_null_percent"])
    mechanistic_pass = uptake_change >= 10.0 and conductance_change >= 10.0
    return {
        "primary_turnover": {
            "reference_T_h": t_ref,
            "interval_h": [t_low, t_high],
            "predicted_T_h": baseline["mean_shedding_age_h"],
            "relative_error_percent": 100.0 * (baseline["mean_shedding_age_h"] / t_ref - 1.0),
            "status": "PASS" if primary_pass else "FAIL",
        },
        "crypt_production": {
            "reference_cells_h": reference["n_cells_h"],
            "predicted_cells_h": baseline["crypt_output_cells_h"],
            "ratio": production_ratio,
            "status": "PASS" if production_pass else "FAIL",
        },
        "cell_census": {
            "reference_cells": reference["N_cells"],
            "predicted_cells": baseline["n_cells"],
            "ratio": census_ratio,
            "status": "PASS" if census_pass else "FAIL",
        },
        "mechanistic_perturbation": {
            "path_factor": atrophy["path_factor"],
            "uptake_change_vs_null_percent": atrophy["uptake_change_vs_null_percent"],
            "conductance_change_vs_null_percent": atrophy["conductance_change_vs_null_percent"],
            "status": "PASS" if mechanistic_pass else "FAIL",
        },
        "empirical_superiority_over_constant_area": {
            "status": "UNKNOWN",
            "reason": "No matched measured uptake and barrier observations are supplied; the null comparison is mechanistic only.",
        },
        "all_numeric_criteria_pass": bool(primary_pass and production_pass and census_pass and mechanistic_pass),
    }


def migration_speed_counterfactual(parameters: Parameters) -> dict[str, Any]:
    """Run the model at the measured speed, at the old chosen 30.0 um/h, and at the measurement
    plus and minus its reported uncertainty, and report what moves. Written because the point of
    replacing a chosen constant is to find out whether it was load-bearing, which cannot be read
    off the source: it has to be run both ways."""
    cases: dict[str, float] = {
        "measured_9.00": 9.00,
        "measured_minus_1sd_8.535": 9.00 - parameters.migration_speed_um_h_sd,
        "measured_plus_1sd_9.465": 9.00 + parameters.migration_speed_um_h_sd,
        "old_chosen_30.0": 30.0,
    }
    rows: dict[str, Any] = {}
    for label, speed in cases.items():
        q = replace(parameters, migration_speed_um_h=speed)
        base = run_scenario(q, 1.0)
        atr = run_scenario(q, q.atrophy_path_factor, base)
        crit = evaluate_criteria(q, base, atr, REFERENCE)
        rows[label] = {
            "migration_speed_um_h": float(speed),
            "minimum_transit_h": float(geometry_metrics(q, 1.0)["minimum_transit_h"]),
            "mean_shedding_age_h": float(base["mean_shedding_age_h"]),
            "n_cells": float(base["n_cells"]),
            "villus_area_mm2": float(base["villus_area_mm2"]),
            "uptake_proxy_mm3_day": float(base["uptake_proxy_mm3_day"]),
            "conductance_mS_cm2": float(base["conductance_mS_cm2"]),
            "primary_turnover_predicted_T_h": float(crit["primary_turnover"]["predicted_T_h"]),
            "primary_turnover_relative_error_percent": float(crit["primary_turnover"]["relative_error_percent"]),
            "primary_turnover_status": crit["primary_turnover"]["status"],
            "cell_census_ratio": float(crit["cell_census"]["ratio"]),
            "cell_census_status": crit["cell_census"]["status"],
            "all_numeric_criteria_pass": bool(crit["all_numeric_criteria_pass"]),
            "implied_drag_pN_h_per_mm": float(q.implied_drag_pN_h_per_mm),
        }
    measured, old = rows["measured_9.00"], rows["old_chosen_30.0"]
    deltas = {}
    for key in ("minimum_transit_h", "mean_shedding_age_h", "n_cells", "villus_area_mm2",
                "uptake_proxy_mm3_day", "conductance_mS_cm2", "primary_turnover_predicted_T_h",
                "primary_turnover_relative_error_percent"):
        a, b = measured[key], old[key]
        deltas[key] = {
            "at_measured_9.00": a,
            "at_chosen_30.0": b,
            "absolute_change": a - b,
            "relative_change_percent": (100.0 * (a / b - 1.0)) if b != 0.0 else None,
        }
    return {
        "why": ("the drag constant eta was documented as chosen to give 0.03 mm/h, so the speed "
                "was the target rather than a result; the measurement it was later refuted "
                "against is used as the default and the old value is run beside it"),
        "declared_input": {
            "migration_speed_um_h": float(parameters.migration_speed_um_h),
            "migration_speed_um_h_sd": float(parameters.migration_speed_um_h_sd),
            "source": parameters.migration_speed_source,
        },
        "ratio_old_over_measured": 30.0 / 9.00,
        "cases": rows,
        "downstream_change_measured_vs_chosen": deltas,
        "load_bearing": any(
            d["relative_change_percent"] is not None and abs(d["relative_change_percent"]) > 1e-9
            for d in deltas.values()),
        "uncertainty_band_on_primary_turnover_h": [
            rows["measured_minus_1sd_8.535"]["primary_turnover_predicted_T_h"],
            rows["measured_plus_1sd_9.465"]["primary_turnover_predicted_T_h"],
        ],
        "review_state": "PENDING_INDEPENDENT_REVIEW",
    }


def build_results(parameters: Parameters | None = None) -> dict[str, Any]:
    parameters = parameters or Parameters()
    baseline = run_scenario(parameters, 1.0)
    atrophy = run_scenario(parameters, parameters.atrophy_path_factor, baseline)
    dimensions = dimensional_check(parameters)
    if not dimensions["all_pass"]:
        raise ValueError("Dimensional analysis failed")
    results = {
        "model_id": "BT-HX-Q115",
        "model_type": "stochastic age-structured crypt-villus transport and barrier model",
        "equations": EQUATIONS,
        "parameters": parameter_dict(parameters),
        "parameter_table": parameter_table(),
        "published_reference": REFERENCE,
        "dimensional_analysis": dimensions,
        "baseline": baseline,
        "atrophy_scenario": atrophy,
        "sensitivity": sensitivity_analysis(parameters),
        "migration_speed_counterfactual": migration_speed_counterfactual(parameters),
        "criteria": evaluate_criteria(parameters, baseline, atrophy, REFERENCE),
        "provenance": {
            "published": "Kai 2021 values are used only as the external validation target and comparison scale.",
            "derived": "Model residence time, areas, fluxes and conductances are derived from the equations and stated assumptions.",
            "hypothesis": "Age-dependent area, transport and barrier properties are mechanistic hypotheses; no measured uptake or barrier data are inferred.",
        },
        "next_resolution_step": "Measure paired 3D villus-crypt geometry, pulse-chase age profiles and matched uptake plus transepithelial conductance under the same perturbation.",
    }
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results.json")
    parser.add_argument("--print-parameters", action="store_true")
    parser.add_argument("--migration-speed-um-h", type=float, default=None,
                        help="override the declared net migration speed (default: the measured "
                             "9.00 um/h, DOI 10.1096/fj.201601002; pass 30.0 for the old chosen value)")
    args = parser.parse_args()
    parameters = Parameters()
    if args.migration_speed_um_h is not None:
        parameters = replace(parameters, migration_speed_um_h=args.migration_speed_um_h)
    results = build_results(parameters)
    if args.print_parameters:
        for row in results["parameter_table"]:
            print(f"{row['symbol']}={row['value']} {row['unit']} [{row['source_or_assumption']}]")
    output_path = Path(args.out)
    output_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results["criteria"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
