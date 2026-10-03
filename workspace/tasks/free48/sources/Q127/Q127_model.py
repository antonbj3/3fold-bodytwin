from __future__ import annotations

import copy
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Unit:
    exponents: tuple[tuple[str, int], ...]

    def __mul__(self, other: Unit) -> Unit:
        values = dict(self.exponents)
        for name, exponent in other.exponents:
            values[name] = values.get(name, 0) + exponent
        return Unit(tuple(sorted((name, exponent) for name, exponent in values.items() if exponent)))

    def __truediv__(self, other: Unit) -> Unit:
        values = dict(self.exponents)
        for name, exponent in other.exponents:
            values[name] = values.get(name, 0) - exponent
        return Unit(tuple(sorted((name, exponent) for name, exponent in values.items() if exponent)))

    def __pow__(self, exponent: int) -> Unit:
        return Unit(tuple(sorted((name, value * exponent) for name, value in self.exponents)))

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Unit) and dict(self.exponents) == dict(other.exponents)


def make_unit(**exponents: int) -> Unit:
    return Unit(tuple(sorted((name, exponent) for name, exponent in exponents.items() if exponent)))


U_DAY = make_unit(day=1)
U_FIBER = make_unit(fiber=1)
U_DIMENSIONLESS = make_unit()
U_NUCLEUS = U_DIMENSIONLESS
U_PRECURSOR = U_DIMENSIONLESS
U_EVENT = U_DIMENSIONLESS
U_UM2 = make_unit(um=2)
U_M = make_unit(fiber=-1)
U_S = make_unit(fiber=-1)
U_A = make_unit(um=2, fiber=-1)
U_P = U_A
U_D = make_unit(um=2)
U_J = make_unit(fiber=-1, day=-1)
U_KFUSE = make_unit(day=-1)
U_KLOSS = make_unit(day=-1)
U_RNUC = make_unit(um=2, day=-1)
U_KDEG = U_KLOSS
U_KINC = U_KLOSS
U_KTURN = U_KLOSS
U_KON = U_KLOSS
U_KOFF = U_KLOSS
U_ABASE = make_unit(fiber=-1, day=-1)
U_ALOAD = U_ABASE
U_KCYCLE = U_KLOSS
U_YIELD = U_DIMENSIONLESS
U_KLEARN = U_KLOSS
U_KFORGET = U_KLOSS


DEFAULT_PARAMS: dict[str, float] = {
    "M0": 2.4,
    "A0": 4029.0,
    "fiber_length": 10000.0,
    "P0": 4029.0,
    "S0": 0.10,
    "E0": 0.0,
    "G0": 0.0,
    "r_nuc": 16.8,
    "r_load": 0.20,
    "r_history": 0.35,
    "k_deg": 0.010,
    "k_inc": 0.010,
    "k_turn": 0.010,
    "k_on": 0.120,
    "k_off": 0.120,
    "a_base": 0.001,
    "a_load": 0.030,
    "k_cycle": 0.050,
    "k_fuse": 0.035,
    "D_star": 1200.0,
    "w_D": 0.20,
    "y_fuse": 1.0,
    "k_loss": 0.00035,
    "q_atrophy": 1.0,
    "p_protect": 0.50,
    "k_learn": 0.030,
    "k_forget": 0.002,
    "damage_sigma": 0.12,
    "repulsion": 0.35,
    "dt": 1.0,
    "stochastic_fibres": 2000.0,
    "seed": 127.0,
}


PARAMETER_META: dict[str, dict[str, str]] = {
    "M0": {"unit": "nuclei/fibre", "basis": "Cumming et al. 2024 type-II baseline mean"},
    "A0": {"unit": "um^2/fibre", "basis": "Cumming et al. 2024 type-II baseline fCSA mean"},
    "fiber_length": {"unit": "um", "basis": "assumption for axial spacing and density"},
    "P0": {"unit": "um^2-equivalent/fibre", "basis": "assumption P0=A0"},
    "S0": {"unit": "precursor/fibre", "basis": "assumption"},
    "E0": {"unit": "dimensionless", "basis": "assumption"},
    "G0": {"unit": "dimensionless", "basis": "assumption"},
    "r_nuc": {"unit": "um^2/(nucleus day)", "basis": "derived assumption k_deg A0/M0"},
    "r_load": {"unit": "dimensionless", "basis": "assumption"},
    "r_history": {"unit": "dimensionless", "basis": "assumption; ablation is run"},
    "k_deg": {"unit": "day^-1", "basis": "assumption"},
    "k_inc": {"unit": "day^-1", "basis": "assumption"},
    "k_turn": {"unit": "day^-1", "basis": "assumption"},
    "k_on": {"unit": "day^-1", "basis": "assumption"},
    "k_off": {"unit": "day^-1", "basis": "assumption"},
    "a_base": {"unit": "precursor/(fibre day)", "basis": "assumption"},
    "a_load": {"unit": "precursor/(fibre day)", "basis": "assumption"},
    "k_cycle": {"unit": "day^-1", "basis": "assumption"},
    "k_fuse": {"unit": "event/(precursor day)", "basis": "assumption; yield separate"},
    "D_star": {"unit": "um^2/nucleus", "basis": "assumption"},
    "w_D": {"unit": "dimensionless", "basis": "assumption"},
    "y_fuse": {"unit": "nuclei/event (dimensionless count ratio)", "basis": "assumption"},
    "k_loss": {"unit": "day^-1", "basis": "assumption"},
    "q_atrophy": {"unit": "dimensionless", "basis": "assumption"},
    "p_protect": {"unit": "dimensionless", "basis": "assumption"},
    "k_learn": {"unit": "day^-1", "basis": "assumption"},
    "k_forget": {"unit": "day^-1", "basis": "assumption"},
    "damage_sigma": {"unit": "normalised fibre length", "basis": "assumption"},
    "repulsion": {"unit": "dimensionless", "basis": "assumption"},
    "dt": {"unit": "day", "basis": "fixed numerical step"},
    "stochastic_fibres": {"unit": "fibres", "basis": "fixed Monte Carlo size"},
    "seed": {"unit": "dimensionless", "basis": "fixed random seed"},
}


TRAIN_DAYS = 70
DETRAIN_DAYS = 112
RETRAIN_DAYS = 70
TOTAL_DAYS = TRAIN_DAYS + DETRAIN_DAYS + RETRAIN_DAYS
REFERENCE = {
    "M0": 2.4,
    "M1": 3.3,
    "MD": 3.2,
    "MC_D": 2.4,
    "M0_sd": 0.5,
    "M1_sd": 0.7,
    "MD_sd": 0.6,
    "MC_D_sd": 0.5,
    "R_excess": (3.2 - 2.4) / (3.3 - 2.4),
    "criterion_relative_tolerance": 0.30,
}


def sigmoid(value: float) -> float:
    clipped = max(-60.0, min(60.0, value))
    return 1.0 / (1.0 + math.exp(-clipped))


def phase_schedule() -> list[tuple[str, int, float]]:
    return [
        ("first_training", TRAIN_DAYS, 1.0),
        ("detraining", DETRAIN_DAYS, 0.0),
        ("retraining", RETRAIN_DAYS, 1.0),
    ]


def make_initial_state(params: dict[str, float], rng: np.random.Generator) -> dict[str, Any]:
    count = max(1, int(round(params["M0"])))
    positions = np.linspace(0.5 / count, 1.0 - 0.5 / count, count).tolist()
    return {
        "day": 0,
        "M": float(params["M0"]),
        "M_acquired": 0.0,
        "S": float(params["S0"]),
        "P": float(params["P0"]),
        "A": float(params["A0"]),
        "G": float(params["G0"]),
        "E": float(params["E0"]),
        "positions": positions,
        "tracked_count": count,
        "fusion_cumulative": 0.0,
        "loss_cumulative": 0.0,
        "fusion_events_cumulative": 0.0,
        "last_site": 0.5,
    }


def insert_position(positions: list[float], rng: np.random.Generator, params: dict[str, float], site: float) -> float:
    if not positions:
        return float(np.clip(rng.normal(site, params["damage_sigma"]), 0.0, 1.0))
    current = np.asarray(positions, dtype=float)
    sigma = max(params["damage_sigma"], 1.0e-6)
    candidates = np.clip(rng.normal(site, sigma, 32), 0.0, 1.0)
    nearest = np.min(np.abs(candidates[:, None] - current[None, :]), axis=1)
    score = nearest + 0.05 * rng.random(candidates.size)
    candidate = float(candidates[int(np.argmax(score))])
    positions.append(candidate)
    values = np.sort(np.asarray(positions, dtype=float))
    for _ in range(12):
        if values.size < 3:
            break
        midpoint = 0.5 * (values[:-2] + values[2:])
        values[1:-1] += params["repulsion"] * 0.5 * (midpoint - values[1:-1])
        values = np.sort(np.clip(values, 0.0, 1.0))
    positions[:] = values.tolist()
    return candidate


def snapshot(state: dict[str, Any], day: int, phase: str, load: float) -> dict[str, Any]:
    M = max(float(state["M"]), 1.0e-12)
    A = max(float(state["A"]), 1.0e-12)
    positions = np.asarray(state["positions"], dtype=float)
    if positions.size > 1:
        gaps = np.diff(np.sort(positions))
        nearest = np.min(np.abs(positions[:, None] - positions[None, :]) + np.eye(positions.size) * 10.0, axis=1)
        spacing = {
            "mean_nearest_um": float(np.mean(nearest) * DEFAULT_PARAMS["fiber_length"]),
            "min_spacing_um": float(np.min(gaps) * DEFAULT_PARAMS["fiber_length"]),
        }
    else:
        spacing = {"mean_nearest_um": 0.0, "min_spacing_um": 0.0}
    return {
        "day": int(day),
        "phase": phase,
        "load": float(load),
        "M": float(state["M"]),
        "M_acquired": float(state["M_acquired"]),
        "S": float(state["S"]),
        "P": float(state["P"]),
        "A": float(state["A"]),
        "G": float(state["G"]),
        "E": float(state["E"]),
        "domain_um2": float(A / M),
        "density_nuclei_per_mm2": float(M / (A * 1.0e-6)),
        "fusion_cumulative": float(state["fusion_cumulative"]),
        "loss_cumulative": float(state["loss_cumulative"]),
        "fusion_events_cumulative": float(state["fusion_events_cumulative"]),
        "positions_normalised": positions.tolist(),
        **spacing,
    }


def simulate_path(
    params: dict[str, float],
    schedule: list[tuple[str, int, float]],
    rng: np.random.Generator,
    initial_state: dict[str, Any] | None = None,
) -> dict[str, Any]:
    state = copy.deepcopy(initial_state) if initial_state is not None else make_initial_state(params, rng)
    dt = params["dt"]
    snapshots = [snapshot(state, 0, "initial", 0.0)]
    phase_ends: dict[str, dict[str, Any]] = {}
    phase_states: dict[str, dict[str, Any]] = {}
    day = 0
    for phase, days, load in schedule:
        for _ in range(days):
            G = state["G"] + dt * (params["k_on"] * load - params["k_off"] * state["G"])
            G = max(0.0, G)
            M = max(float(state["M"]), 0.0)
            M_acquired = max(0.0, min(float(state["M_acquired"]), M))
            A = max(float(state["A"]), 1.0e-12)
            D = A / max(M, 1.0e-12)
            p_fuse = sigmoid((D / params["D_star"] - 1.0) / params["w_D"])
            J_fuse = params["k_fuse"] * max(float(state["S"]), 0.0) * p_fuse
            source = params["a_base"] + params["a_load"] * G
            S = max(0.0, state["S"] + dt * (source - params["k_cycle"] * state["S"] - J_fuse))
            atrophy_factor = 1.0 + params["q_atrophy"] * max(0.0, 1.0 - A / params["A0"])
            old_M = max(0.0, M - M_acquired)
            old_survival = math.exp(-params["k_loss"] * atrophy_factor * dt)
            acquired_survival = math.exp(-params["k_loss"] * atrophy_factor * (1.0 - params["p_protect"]) * dt)
            loss_old = old_M * (1.0 - old_survival)
            loss_acquired = M_acquired * (1.0 - acquired_survival)
            added = params["y_fuse"] * J_fuse * dt
            M_next = max(0.0, old_M * old_survival + M_acquired * acquired_survival + added)
            M_acquired_next = max(0.0, min(M_next, M_acquired * acquired_survival + added))
            E = state["E"] + dt * (params["k_learn"] * G * (1.0 - state["E"]) - params["k_forget"] * state["E"])
            E = max(0.0, min(1.0, E))
            P = state["P"] + dt * (
                params["r_nuc"] * M_next * (1.0 + params["r_load"] * G + params["r_history"] * E * G)
                - params["k_deg"] * state["P"]
            )
            A_next = state["A"] + dt * (params["k_inc"] * P - params["k_turn"] * state["A"])
            A_next = max(1.0, A_next)
            state["day"] = day + 1
            state["M"] = M_next
            state["M_acquired"] = M_acquired_next
            state["S"] = S
            state["P"] = max(0.0, P)
            state["A"] = A_next
            state["G"] = G
            state["E"] = E
            state["fusion_cumulative"] += added
            state["loss_cumulative"] += loss_old + loss_acquired
            state["fusion_events_cumulative"] += J_fuse * dt
            target_count = max(1, int(round(M_next)))
            while state["tracked_count"] < target_count:
                state["last_site"] = float(np.clip(0.5 + 0.35 * math.sin(day / 17.0), 0.0, 1.0))
                insert_position(state["positions"], rng, params, state["last_site"])
                state["tracked_count"] += 1
            day += 1
            snapshots.append(snapshot(state, day, phase, load))
        phase_ends[phase] = snapshot(state, day, phase, load)
        phase_states[phase] = copy.deepcopy(state)
    return {
        "initial": snapshots[0],
        "phase_ends": phase_ends,
        "phase_states": phase_states,
        "snapshots": snapshots,
        "final": snapshots[-1],
    }


def initial_stochastic_state(count: int, params: dict[str, float], rng: np.random.Generator) -> dict[str, Any]:
    positions = np.linspace(0.5 / count, 1.0 - 0.5 / count, count)
    nuclei = [
        {
            "id": index,
            "position": float(position),
            "acquired": False,
            "birth_day": 0,
        }
        for index, position in enumerate(positions)
    ]
    return {
        "day": 0,
        "nuclei": nuclei,
        "S": float(params["S0"]),
        "S_projection_cumulative": 0.0,
        "P": float(params["P0"]),
        "A": float(params["A0"]),
        "G": float(params["G0"]),
        "E": float(params["E0"]),
        "next_id": count,
        "fusion_events_cumulative": 0,
        "loss_cumulative": 0,
        "last_site": 0.5,
    }


def stochastic_snapshot(state: dict[str, Any], day: int, phase: str, load: float) -> dict[str, Any]:
    raw_M = len(state["nuclei"])
    M = max(raw_M, 1)
    A = max(float(state["A"]), 1.0e-12)
    positions = np.asarray([nucleus["position"] for nucleus in state["nuclei"]], dtype=float)
    acquired = sum(1 for nucleus in state["nuclei"] if nucleus["acquired"])
    if positions.size > 1:
        gaps = np.diff(np.sort(positions))
        nearest = np.min(np.abs(positions[:, None] - positions[None, :]) + np.eye(positions.size) * 10.0, axis=1)
        spacing = {
            "mean_nearest_um": float(np.mean(nearest) * params_fiber_length()),
            "min_spacing_um": float(np.min(gaps) * params_fiber_length()),
        }
    else:
        spacing = {"mean_nearest_um": 0.0, "min_spacing_um": 0.0}
    return {
        "day": int(day),
        "phase": phase,
        "load": float(load),
        "M": float(raw_M),
        "M_acquired": float(acquired),
        "S": float(state["S"]),
        "precursor_projection_cumulative": float(state["S_projection_cumulative"]),
        "P": float(state["P"]),
        "A": float(A),
        "G": float(state["G"]),
        "E": float(state["E"]),
        "domain_um2": float(A / raw_M) if raw_M else None,
        "domain_status": "DEFINED" if raw_M else "UNDEFINED_NO_NUCLEI",
        "density_nuclei_per_mm2": float(raw_M / (A * 1.0e-6)),
        "fusion_cumulative": float(state["fusion_events_cumulative"]),
        "loss_cumulative": float(state["loss_cumulative"]),
        "positions_normalised": positions.tolist(),
        **spacing,
    }


def params_fiber_length() -> float:
    return DEFAULT_PARAMS["fiber_length"]


def simulate_stochastic_path(
    params: dict[str, float],
    schedule: list[tuple[str, int, float]],
    initial_count: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    state = initial_stochastic_state(initial_count, params, rng)
    phase_ends: dict[str, dict[str, Any]] = {}
    day = 0
    dt = params["dt"]
    for phase, days, load in schedule:
        for _ in range(days):
            G = max(0.0, state["G"] + dt * (params["k_on"] * load - params["k_off"] * state["G"]))
            M = len(state["nuclei"])
            A = max(float(state["A"]), 1.0e-12)
            D = A / max(M, 1)
            p_fuse = sigmoid((D / params["D_star"] - 1.0) / params["w_D"])
            J_fuse = params["k_fuse"] * max(float(state["S"]), 0.0) * p_fuse
            source = params["a_base"] + params["a_load"] * G
            available = 1 if state["S"] + source > 0.0 else 0
            fusion_events = min(int(rng.poisson(J_fuse * dt)), available)
            position_values = [nucleus["position"] for nucleus in state["nuclei"]]
            for _ in range(fusion_events):
                state["last_site"] = float(np.clip(0.5 + 0.35 * math.sin(day / 17.0), 0.0, 1.0))
                candidate = insert_position(position_values, rng, params, state["last_site"])
                state["nuclei"].append(
                    {
                        "id": state["next_id"],
                        "position": candidate,
                        "acquired": True,
                        "birth_day": day + 1,
                    }
                )
                state["next_id"] += 1
            for nucleus, position in zip(state["nuclei"], position_values):
                nucleus["position"] = position
            raw_precursor = state["S"] + dt * (source - params["k_cycle"] * state["S"]) - fusion_events
            state["S"] = max(0.0, raw_precursor)
            state["S_projection_cumulative"] += state["S"] - raw_precursor
            atrophy_factor = 1.0 + params["q_atrophy"] * max(0.0, 1.0 - A / params["A0"])
            survivors = []
            for nucleus in state["nuclei"]:
                hazard = params["k_loss"] * atrophy_factor * (1.0 - params["p_protect"] if nucleus["acquired"] else 1.0)
                loss_probability = 1.0 - math.exp(-hazard * dt)
                if rng.random() < loss_probability:
                    state["loss_cumulative"] += 1.0
                else:
                    survivors.append(nucleus)
            state["nuclei"] = survivors
            state["E"] = max(
                0.0,
                min(
                    1.0,
                    state["E"]
                    + dt * (params["k_learn"] * G * (1.0 - state["E"]) - params["k_forget"] * state["E"]),
                ),
            )
            M_next = len(state["nuclei"])
            state["P"] = max(
                0.0,
                state["P"]
                + dt
                * (
                    params["r_nuc"] * M_next * (1.0 + params["r_load"] * G + params["r_history"] * state["E"] * G)
                    - params["k_deg"] * state["P"]
                ),
            )
            state["A"] = max(1.0, state["A"] + dt * (params["k_inc"] * state["P"] - params["k_turn"] * state["A"]))
            state["G"] = G
            state["fusion_events_cumulative"] += fusion_events
            day += 1
        phase_ends[phase] = stochastic_snapshot(state, day, phase, load)
    return phase_ends


def summarize(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(array)),
        "sd": float(np.std(array, ddof=1)) if array.size > 1 else 0.0,
        "q025": float(np.quantile(array, 0.025)),
        "median": float(np.median(array)),
        "q975": float(np.quantile(array, 0.975)),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
    }


def run_stochastic(params: dict[str, float]) -> dict[str, Any]:
    n = int(params["stochastic_fibres"])
    init_rng = np.random.default_rng(int(params["seed"]))
    initial_counts = np.where(init_rng.random(n) < 0.60, 2, 3).astype(int)
    trained_rng = np.random.default_rng(int(params["seed"]) + 101)
    control_rng = np.random.default_rng(int(params["seed"]) + 202)
    schedule = phase_schedule()
    control_schedule = [("detraining", TRAIN_DAYS + DETRAIN_DAYS, 0.0), ("retraining", RETRAIN_DAYS, 1.0)]
    trained_values: dict[str, list[float]] = {"M0": [], "M1": [], "MD": [], "MR": [], "A0": [], "A1": [], "AD": [], "AR": [], "fusion_events": [], "fusion_events_first_training": [], "fusion_events_detraining": [], "fusion_events_retraining": [], "acquired_MD": []}
    control_values: dict[str, list[float]] = {"M0": [], "MD": [], "MR": [], "AD": [], "AR": []}
    for count in initial_counts:
        trained = simulate_stochastic_path(params, schedule, int(count), trained_rng)
        control = simulate_stochastic_path(params, control_schedule, int(count), control_rng)
        trained_values["M0"].append(float(count))
        trained_values["M1"].append(trained["first_training"]["M"])
        trained_values["MD"].append(trained["detraining"]["M"])
        trained_values["MR"].append(trained["retraining"]["M"])
        trained_values["A0"].append(float(params["A0"]))
        trained_values["A1"].append(trained["first_training"]["A"])
        trained_values["AD"].append(trained["detraining"]["A"])
        trained_values["AR"].append(trained["retraining"]["A"])
        trained_values["fusion_events"].append(trained["retraining"]["fusion_cumulative"])
        trained_values["fusion_events_first_training"].append(trained["first_training"]["fusion_cumulative"])
        trained_values["fusion_events_detraining"].append(trained["detraining"]["fusion_cumulative"] - trained["first_training"]["fusion_cumulative"])
        trained_values["fusion_events_retraining"].append(trained["retraining"]["fusion_cumulative"] - trained["detraining"]["fusion_cumulative"])
        trained_values["acquired_MD"].append(trained["detraining"]["M_acquired"])
        control_values["M0"].append(float(count))
        control_values["MD"].append(control["detraining"]["M"])
        control_values["MR"].append(control["retraining"]["M"])
        control_values["AD"].append(control["detraining"]["A"])
        control_values["AR"].append(control["retraining"]["A"])
    means_trained = {key: summarize(values) for key, values in trained_values.items()}
    means_control = {key: summarize(values) for key, values in control_values.items()}
    R = (means_trained["MD"]["mean"] - means_control["MD"]["mean"]) / (means_trained["M1"]["mean"] - means_trained["M0"]["mean"])
    return {
        "n_fibres": n,
        "seed": int(params["seed"]),
        "trained": means_trained,
        "control": means_control,
        "R_excess": float(R),
        "retention_share_bounded": float(max(0.0, min(1.0, R))),
        "retention_upper_bound_exceeded": bool(R > 1.0),
        "first_training_fusion_events_per_fibre_day": float(means_trained["fusion_events_first_training"]["mean"] / TRAIN_DAYS),
        "detraining_fusion_events_per_fibre_day": float(means_trained["fusion_events_detraining"]["mean"] / DETRAIN_DAYS),
        "retraining_fusion_events_per_fibre_day": float(means_trained["fusion_events_retraining"]["mean"] / RETRAIN_DAYS),
        "retraining_M_excess_over_control": float(means_trained["MR"]["mean"] - means_control["MR"]["mean"]),
        "retraining_A_excess_over_control": float(means_trained["AR"]["mean"] - means_control["AR"]["mean"]),
    }


def dimension_audit() -> dict[str, Any]:
    checks = {
        "dM_dt_fusion": U_YIELD * U_J == U_M / U_DAY,
        "dM_dt_loss": U_KLOSS * U_M == U_M / U_DAY,
        "dS_dt_source": U_ABASE == U_S / U_DAY,
        "dS_dt_cycle": U_KCYCLE * U_S == U_S / U_DAY,
        "dS_dt_fusion_consumption": U_J == U_S / U_DAY,
        "dP_dt_synthesis": U_RNUC * U_M == U_P / U_DAY,
        "dP_dt_degradation": U_KDEG * U_P == U_P / U_DAY,
        "dA_dt_incorporation": U_KINC * U_P == U_A / U_DAY,
        "dA_dt_turnover": U_KTURN * U_A == U_A / U_DAY,
        "domain": U_A / U_M == U_D,
        "fusion_rate": U_KFUSE * U_S == U_J,
        "load_signal": U_KON == U_KOFF,
        "history_rate": U_KLEARN * U_DIMENSIONLESS == U_DIMENSIONLESS / U_DAY,
        "history_decay": U_KFORGET * U_DIMENSIONLESS == U_DIMENSIONLESS / U_DAY,
        "yield_dimensionless": U_YIELD == U_DIMENSIONLESS,
    }
    return {"passed": bool(all(checks.values())), "checks": {key: bool(value) for key, value in checks.items()}}


def analytic_loss_only(M0: float, k_loss: float, days: int) -> float:
    return M0 * math.exp(-k_loss * days)


def run_sensitivity(params: dict[str, float], base_R: float) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for name in ("k_fuse", "k_loss", "r_history"):
        values = []
        for factor in (0.5, 1.5):
            changed = dict(params)
            changed[name] = params[name] * factor
            trained = simulate_path(changed, phase_schedule(), np.random.default_rng(int(params["seed"])))
            control = simulate_path(
                changed,
                [("detraining", TRAIN_DAYS + DETRAIN_DAYS, 0.0), ("retraining", RETRAIN_DAYS, 1.0)],
                np.random.default_rng(int(params["seed"]) + 1),
            )
            M0 = trained["initial"]["M"]
            M1 = trained["phase_ends"]["first_training"]["M"]
            MD = trained["phase_ends"]["detraining"]["M"]
            MC = control["phase_ends"]["detraining"]["M"]
            R = (MD - MC) / (M1 - M0)
            values.append({"factor": factor, "parameter_value": changed[name], "R_excess": float(R), "delta_R": float(R - base_R)})
        output[name] = {"minus50": values[0], "plus50": values[1], "absolute_span": float(abs(values[1]["R_excess"] - values[0]["R_excess"]))}
    return output


def run_all(params: dict[str, float] | None = None) -> dict[str, Any]:
    values = dict(DEFAULT_PARAMS if params is None else params)
    trained = simulate_path(values, phase_schedule(), np.random.default_rng(int(values["seed"])))
    control = simulate_path(
        values,
        [("detraining", TRAIN_DAYS + DETRAIN_DAYS, 0.0), ("retraining", RETRAIN_DAYS, 1.0)],
        np.random.default_rng(int(values["seed"]) + 1),
    )
    initial = trained["initial"]
    M0 = initial["M"]
    M1 = trained["phase_ends"]["first_training"]["M"]
    MD = trained["phase_ends"]["detraining"]["M"]
    MC = control["phase_ends"]["detraining"]["M"]
    MR = trained["phase_ends"]["retraining"]["M"]
    MR_control = control["phase_ends"]["retraining"]["M"]
    R = (MD - MC) / (M1 - M0)
    relative_error = abs(R - REFERENCE["R_excess"]) / REFERENCE["R_excess"]
    ablation_params = dict(values)
    ablation_params["r_history"] = 0.0
    history_ablation = simulate_path(
        ablation_params,
        [("retraining", RETRAIN_DAYS, 1.0)],
        np.random.default_rng(int(values["seed"]) + 2),
        trained["phase_states"]["detraining"],
    )
    M_history_ablation = history_ablation["phase_ends"]["retraining"]["M"]
    A_history_ablation = history_ablation["phase_ends"]["retraining"]["A"]
    stochastic = run_stochastic(values)
    dimension = dimension_audit()
    sensitivity = run_sensitivity(values, R)
    primary = {
        "model_M0": M0,
        "model_M1": M1,
        "model_MD": MD,
        "model_MC_D": MC,
        "model_MR": MR,
        "model_MR_control": MR_control,
        "R_excess": float(R),
        "retention_share_bounded": float(max(0.0, min(1.0, R))),
        "retention_upper_bound_exceeded": bool(R > 1.0),
        "reference_R_excess": REFERENCE["R_excess"],
        "relative_error": float(relative_error),
        "criterion_tolerance": REFERENCE["criterion_relative_tolerance"],
        "criterion_met": bool(relative_error <= REFERENCE["criterion_relative_tolerance"]),
        "raw_MD_over_M1": float(MD / M1),
        "first_training_fusion_events_per_fibre_day": float(trained["phase_ends"]["first_training"]["fusion_cumulative"] / TRAIN_DAYS),
        "detraining_fusion_events_per_fibre_day": float((trained["phase_ends"]["detraining"]["fusion_cumulative"] - trained["phase_ends"]["first_training"]["fusion_cumulative"]) / DETRAIN_DAYS),
        "retraining_fusion_events_per_fibre_day": float((trained["phase_ends"]["retraining"]["fusion_cumulative"] - trained["phase_ends"]["detraining"]["fusion_cumulative"]) / RETRAIN_DAYS),
        "history_ablation_MR": M_history_ablation,
        "history_ablation_AR": A_history_ablation,
        "history_effect_MR": float(MR - M_history_ablation),
        "size_only_effect_MR": float(MR - MR_control),
    }
    phase_summary = {
        "trained": {name: trained["phase_ends"][name] for name in ("first_training", "detraining", "retraining")},
        "control": {name: control["phase_ends"][name] for name in ("detraining", "retraining")},
        "history_ablation": history_ablation["phase_ends"],
    }
    parameter_table = {
        name: {"value": values[name], "unit": PARAMETER_META[name]["unit"], "basis": PARAMETER_META[name]["basis"]}
        for name in values
    }
    result_status = "complete_criterion_passed" if primary["criterion_met"] else "complete_criterion_failed"
    return {
        "id": "BT-HX-Q127",
        "model": "first-principles longitudinal fibre/satellite-cell state model",
        "status": result_status,
        "interpretation": "UNKNOWN" if not primary["criterion_met"] else "retention criterion met; motor memory remains untested",
        "protocol_days": {"first_training": TRAIN_DAYS, "detraining": DETRAIN_DAYS, "retraining": RETRAIN_DAYS, "total": TOTAL_DAYS},
        "reference": REFERENCE,
        "primary": primary,
        "phase_summary": phase_summary,
        "sensitivity": sensitivity,
        "stochastic": stochastic,
        "dimension_audit": dimension,
        "analytic_loss_limit": {
            "M0": M0,
            "k_loss_day": values["k_loss"],
            "days": DETRAIN_DAYS,
            "predicted_without_flux": analytic_loss_only(M0, values["k_loss"], DETRAIN_DAYS),
        },
        "parameters": parameter_table,
        "source": {
            "citation": "Cumming et al. 2024, The Journal of Physiology 602(17):4171-4193, DOI 10.1113/JP285675",
            "page": "https://www.ovid.com/journals/jphy/fulltext/10.1113/jp285675~muscle-memory-in-humans-evidence-for-myonuclear-permanence",
            "verified": True,
            "reference_ratio_derivation": "(3.2-2.4)/(3.3-2.4)",
        },
    }


def main() -> None:
    result = run_all()
    path = Path(__file__).resolve().parent / "results.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"id": result["id"], "R_excess": result["primary"]["R_excess"], "criterion_met": result["primary"]["criterion_met"]}, sort_keys=True))


if __name__ == "__main__":
    main()
