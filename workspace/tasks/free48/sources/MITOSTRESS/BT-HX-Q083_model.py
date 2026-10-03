import copy
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


@dataclass
class Compartment:
    volume: float
    amounts: Dict[str, float]
    psi_mv: float
    label: str = ""

    def clone(self, label: Optional[str] = None) -> "Compartment":
        return Compartment(
            volume=float(self.volume),
            amounts={key: float(value) for key, value in self.amounts.items()},
            psi_mv=float(self.psi_mv),
            label=self.label if label is None else label,
        )

    def validate(self) -> None:
        if not np.isfinite(self.volume) or self.volume <= 0:
            raise ValueError("volume must be positive and finite")
        if not np.isfinite(self.psi_mv):
            raise ValueError("potential must be finite")
        for species, amount in self.amounts.items():
            if not np.isfinite(amount) or amount < 0:
                raise ValueError(f"invalid amount for {species}")

    def concentration(self, species: str) -> float:
        return float(self.amounts.get(species, 0.0) / self.volume)


@dataclass
class FusedPair:
    left: Compartment
    right: Compartment

    def validate(self) -> None:
        self.left.validate()
        self.right.validate()


PARAMETER_TABLE: Dict[str, Dict[str, Dict[str, object]]] = {
    "mitochondrion": {
        "V0": {"value": 1.0, "unit": "um^3", "source_or_assumption": "normaliserat initialfall"},
        "tau_substrate_s": {"value": 12.0, "unit": "s", "source_or_assumption": "Eisner et al. 2017, Fig. 3A; fast mixing upper boundary"},
        "tau_enzyme_s": {"value": 600.0, "unit": "s", "source_or_assumption": "assumption: protein memory is slower"},
        "tau_damage_s": {"value": 600.0, "unit": "s", "source_or_assumption": "antagande: skade-/reparationsminne"},
        "tau_opa1_s": {"value": 600.0, "unit": "s", "source_or_assumption": "antagande: OPA1-relaterat reglerande minne"},
        "tau_psi_s": {"value": 5.0, "unit": "s", "source_or_assumption": "antagande: elektrisk koppling under fusion"},
        "tau_post_s": {"value": 77.0, "unit": "s", "source_or_assumption": "Twig et al. 2008, Fig. 2B"},
        "delta_psi_mv": {"value": 6.0, "unit": "mV", "source_or_assumption": "Twig et al. 2008, Fig. 3B; mean magnitude"},
        "q": {"value": 0.5, "unit": "fraction", "source_or_assumption": "antagande: symmetriskt basfall"},
        "partition_copy_number": {"value": 10000.0, "unit": "copies", "source_or_assumption": "antagande for finite-copy partition"},
        "k_atp_nmol_um3_s": {"value": 0.08, "unit": "nmol um^-3 s^-1", "source_or_assumption": "antagande for normalized response"},
        "K_M_nmol_um3": {"value": 0.8, "unit": "nmol um^-3", "source_or_assumption": "antagande for Michaelis-Menten term"},
        "psi_ref_mv": {"value": -180.0, "unit": "mV", "source_or_assumption": "antagande for normalized response"},
        "psi_scale_mv": {"value": 25.0, "unit": "mV", "source_or_assumption": "antagande for potential sensitivity"},
        "k_leak_nmol_um3_s": {"value": 0.01, "unit": "nmol um^-3 s^-1", "source_or_assumption": "antagande for leak term"},
        "K_opa1_nmol_um3": {"value": 0.6, "unit": "nmol um^-3", "source_or_assumption": "antagande for OPA1-dependent hazard"},
        "n_opa1": {"value": 4.0, "unit": "dimensionless", "source_or_assumption": "antagande: cooperative selectivity"},
        "k_fusion_base_s": {"value": 0.01, "unit": "s^-1", "source_or_assumption": "antagande for next-event hazard"},
        "risk_window_s": {"value": 180.0, "unit": "s", "source_or_assumption": "Twig et al. 2008, Fig. 4; 3 min"},
        "substrate_perturbation": {"value": 0.5, "unit": "fraction", "source_or_assumption": "fryst kontrollfall"},
    },
    "vesicle": {
        "V0": {"value": 0.01, "unit": "um^3", "source_or_assumption": "normaliserad vesikelvolym"},
        "tau_cargo_s": {"value": 30.0, "unit": "s", "source_or_assumption": "antagande: loslig vesiclecargo"},
        "tau_activation_s": {"value": 600.0, "unit": "s", "source_or_assumption": "antagande: priming/SNARE state persists"},
        "tau_enzyme_s": {"value": 600.0, "unit": "s", "source_or_assumption": "antagande: protein state persists"},
        "tau_damage_s": {"value": 600.0, "unit": "s", "source_or_assumption": "antagande: quality state persists"},
        "tau_psi_s": {"value": 5.0, "unit": "s", "source_or_assumption": "antagande: shared electrical state"},
        "tau_post_s": {"value": 30.0, "unit": "s", "source_or_assumption": "antagande for vesicle partition window"},
        "delta_psi_mv": {"value": 0.0, "unit": "mV", "source_or_assumption": "vesicle potential is not used in release law"},
        "q": {"value": 0.5, "unit": "fraction", "source_or_assumption": "antagande: symmetriskt basfall"},
        "partition_copy_number": {"value": 10000.0, "unit": "copies", "source_or_assumption": "antagande for finite-copy partition"},
        "p_single_channel": {"value": 0.06, "unit": "fraction", "source_or_assumption": "Luo et al. 2015, Fig. 5C"},
        "k_fusion_vesicle_s": {"value": 0.0037, "unit": "s^-1", "source_or_assumption": "Kim et al. 2012, Fig. 3E"},
        "vesicle_dt_s": {"value": 1.0, "unit": "s", "source_or_assumption": "definerat prediktionsfonster"},
        "calcium_perturbation": {"value": 0.5, "unit": "fraction", "source_or_assumption": "fryst kontrollfall"},
    },
}


MITOCHONDRIAL_DEFAULTS: Dict[str, float] = {
    key: float(value["value"])
    for key, value in PARAMETER_TABLE["mitochondrion"].items()
}
VESICLE_DEFAULTS: Dict[str, float] = {
    key: float(value["value"])
    for key, value in PARAMETER_TABLE["vesicle"].items()
}

LOBE_SPECIES = {"enzyme", "damage", "opa1", "activation"}


SOURCES: Dict[str, Dict[str, object]] = {
    "twig_2008": {
        "citation": "Twig et al., EMBO Journal 27, 433-446 (2008)",
        "doi": "10.1038/sj.emboj.7601963",
        "verified": True,
        "values": {
            "post_fusion_connected_s_ins1": 77.0,
            "post_fusion_connected_sd_s_ins1": 71.0,
            "post_fusion_connected_s_cos7": 87.0,
            "post_fusion_connected_sd_s_cos7": 78.0,
            "depolarizing_daughter_mv": 5.9,
            "depolarizing_daughter_sd_mv": 3.7,
            "hyperpolarizing_daughter_mv": 6.5,
            "hyperpolarizing_daughter_sd_mv": 3.8,
            "hyperpolarized_fusion_ratio": 6.0,
        },
        "locations": ["Fig. 2B", "Fig. 3B", "Fig. 4"],
    },
    "eisner_2017": {
        "citation": "Eisner et al., PNAS 114, E859-E868 (2017)",
        "doi": "10.1073/pnas.1617288114",
        "verified": True,
        "values": {
            "fresh_avcm_fusion_events_per_min": 1.4,
            "fresh_avcm_fusion_events_sd_per_min": 0.1,
            "fast_mixing_upper_s": 12.0,
            "slow_mixing_example_s": 70.0,
            "pa_gfp_decay_at_300s_percent": 25.6,
        },
        "locations": ["Fig. 1C", "Fig. 2C", "Fig. 3A"],
    },
    "luo_2015": {
        "citation": "Luo et al., Journal of Neurophysiology 113, 2480-2489 (2015)",
        "doi": "10.1152/jn.00879.2014",
        "verified": True,
        "values": {"single_channel_release_probability": 0.06},
        "locations": ["Fig. 5C"],
    },
    "kim_2012": {
        "citation": "Kim et al., EMBO Journal 31, 2144-2155 (2012)",
        "doi": "10.1038/emboj.2012.57",
        "verified": True,
        "values": {
            "syt1_fusion_rate_s": 0.0037,
            "syt1_fusion_rate_sd_s": 0.0009,
            "snare_only_fusion_rate_s": 0.000116,
            "snare_only_fusion_rate_sd_s": 0.000025,
        },
        "locations": ["Fig. 3E", "Fig. 3F"],
    },
}


def _species_union(left: Compartment, right: Compartment) -> List[str]:
    return sorted(set(left.amounts) | set(right.amounts))


def _tau_for(params: Mapping[str, float], species: str) -> float:
    key = f"tau_{species}_s"
    if key in params:
        return float(params[key])
    return float(params.get("tau_default_s", 600.0))


def _exchange_coefficient(volume_left: float, volume_right: float, tau_s: float) -> float:
    if tau_s <= 0 or volume_left <= 0 or volume_right <= 0:
        raise ValueError("exchange volumes and time constants must be positive")
    return 1.0 / (tau_s * (1.0 / volume_left + 1.0 / volume_right))


def _weighted_potential(left: Compartment, right: Compartment) -> float:
    total_volume = left.volume + right.volume
    return float((left.volume * left.psi_mv + right.volume * right.psi_mv) / total_volume)


def _tau_map(params: Mapping[str, float], species: str) -> float:
    return _tau_for(params, species)


def fuse(
    left: Compartment,
    right: Compartment,
    params: Mapping[str, float],
    duration_s: float,
) -> FusedPair:
    left.validate()
    right.validate()
    if duration_s < 0:
        raise ValueError("duration must be non-negative")
    total_volume = left.volume + right.volume
    fused_amounts_left: Dict[str, float] = {}
    fused_amounts_right: Dict[str, float] = {}
    for species in _species_union(left, right):
        amount_left = float(left.amounts.get(species, 0.0))
        amount_right = float(right.amounts.get(species, 0.0))
        total_amount = amount_left + amount_right
        equilibrium_left = total_amount * left.volume / total_volume
        equilibrium_right = total_amount * right.volume / total_volume
        tau_s = _tau_map(params, species)
        decay = math.exp(-duration_s / tau_s)
        fused_amounts_left[species] = equilibrium_left + (amount_left - equilibrium_left) * decay
        fused_amounts_right[species] = equilibrium_right + (amount_right - equilibrium_right) * decay
    potential_equilibrium = _weighted_potential(left, right)
    potential_decay = math.exp(-duration_s / float(params["tau_psi_s"]))
    fused_left = Compartment(
        volume=left.volume,
        amounts=fused_amounts_left,
        psi_mv=potential_equilibrium + (left.psi_mv - potential_equilibrium) * potential_decay,
        label=f"{left.label or 'A'}+{right.label or 'B'}_left",
    )
    fused_right = Compartment(
        volume=right.volume,
        amounts=fused_amounts_right,
        psi_mv=potential_equilibrium + (right.psi_mv - potential_equilibrium) * potential_decay,
        label=f"{left.label or 'A'}+{right.label or 'B'}_right",
    )
    pair = FusedPair(fused_left, fused_right)
    pair.validate()
    return pair


def fission(
    pair: FusedPair,
    params: Mapping[str, float],
    rng: Optional[np.random.Generator] = None,
    stochastic: bool = False,
) -> Tuple[Compartment, Compartment]:
    pair.validate()
    total_volume = pair.left.volume + pair.right.volume
    q = float(np.clip(params["q"], 1.0e-9, 1.0 - 1.0e-9))
    total_amounts: Dict[str, float] = {}
    for species in _species_union(pair.left, pair.right):
        total_amounts[species] = float(
            pair.left.amounts.get(species, 0.0) + pair.right.amounts.get(species, 0.0)
        )
    amounts_left: Dict[str, float] = {}
    amounts_right: Dict[str, float] = {}
    copy_number = float(params.get("partition_copy_number", 10000.0))
    target_volume = q * total_volume
    left_allocation = min(1.0, target_volume / pair.left.volume)
    right_allocation = float(
        np.clip(
            (target_volume - left_allocation * pair.left.volume) / pair.right.volume,
            0.0,
            1.0,
        )
    )
    for species, total_amount in total_amounts.items():
        if species in LOBE_SPECIES:
            amount_left = (
                pair.left.amounts.get(species, 0.0) * left_allocation
                + pair.right.amounts.get(species, 0.0) * right_allocation
            )
            amount_right = total_amount - amount_left
        else:
            fraction = q
            if stochastic and rng is not None:
                standard_fraction = math.sqrt(q * (1.0 - q) / max(copy_number, 1.0))
                fraction = float(np.clip(rng.normal(q, standard_fraction), 0.0, 1.0))
            amount_left = total_amount * fraction
            amount_right = total_amount - amount_left
        amounts_left[species] = float(amount_left)
        amounts_right[species] = float(amount_right)
    potential_mean = _weighted_potential(pair.left, pair.right)
    delta_psi = float(params.get("delta_psi_mv", 0.0))
    potential_left = potential_mean - 2.0 * (1.0 - q) * delta_psi
    potential_right = potential_mean + 2.0 * q * delta_psi
    daughter_left = Compartment(
        volume=q * total_volume,
        amounts=amounts_left,
        psi_mv=potential_left,
        label="daughter_1",
    )
    daughter_right = Compartment(
        volume=(1.0 - q) * total_volume,
        amounts=amounts_right,
        psi_mv=potential_right,
        label="daughter_2",
    )
    daughter_left.validate()
    daughter_right.validate()
    return daughter_left, daughter_right


def initial_mitochondria() -> Tuple[Compartment, Compartment]:
    return (
        Compartment(
            volume=1.0,
            amounts={
                "substrate": 1.20,
                "enzyme": 0.80,
                "damage": 0.10,
                "opa1": 0.90,
            },
            psi_mv=-175.0,
            label="A",
        ),
        Compartment(
            volume=1.0,
            amounts={
                "substrate": 0.50,
                "enzyme": 0.40,
                "damage": 0.60,
                "opa1": 0.30,
            },
            psi_mv=-150.0,
            label="B",
        ),
    )


def initial_vesicles() -> Tuple[Compartment, Compartment]:
    return (
        Compartment(
            volume=0.01,
            amounts={
                "cargo": 0.0020,
                "activation": 0.0090,
                "enzyme": 0.0010,
                "damage": 0.0001,
                "substrate": 0.0020,
            },
            psi_mv=-40.0,
            label="primed",
        ),
        Compartment(
            volume=0.01,
            amounts={
                "cargo": 0.0010,
                "activation": 0.0020,
                "enzyme": 0.0005,
                "damage": 0.0005,
                "substrate": 0.0010,
            },
            psi_mv=-30.0,
            label="unprimed",
        ),
    )


def atp_response(
    compartment: Compartment,
    params: Mapping[str, float],
    substrate_scale: float = 1.0,
) -> float:
    substrate = compartment.concentration("substrate") * float(substrate_scale)
    enzyme = compartment.concentration("enzyme")
    damage = compartment.concentration("damage")
    potential_factor = math.exp(
        (compartment.psi_mv - float(params["psi_ref_mv"])) / float(params["psi_scale_mv"])
    )
    gross = (
        compartment.volume
        * float(params["k_atp_nmol_um3_s"])
        * enzyme
        * substrate
        / (float(params["K_M_nmol_um3"]) + substrate)
        * potential_factor
        * math.exp(-damage)
    )
    leak = compartment.volume * float(params["k_leak_nmol_um3_s"])
    return float(max(gross - leak, 0.0))


def fusion_hazard(compartment: Compartment, params: Mapping[str, float]) -> float:
    opa1_concentration = compartment.concentration("opa1")
    hill = float(params["n_opa1"])
    k_opa1 = float(params["K_opa1_nmol_um3"])
    opa1_factor = opa1_concentration**hill / (k_opa1**hill + opa1_concentration**hill)
    potential_factor = math.exp(
        (compartment.psi_mv - float(params["psi_ref_mv"])) / float(params["psi_scale_mv"])
    )
    damage_factor = math.exp(-compartment.concentration("damage"))
    return float(params["k_fusion_base_s"]) * opa1_factor * potential_factor * damage_factor


def fusion_probability(compartment: Compartment, params: Mapping[str, float]) -> float:
    window_s = float(params.get("risk_window_s", 180.0))
    return float(-math.expm1(-fusion_hazard(compartment, params) * window_s))


def vesicle_release_probability(
    compartment: Compartment,
    params: Mapping[str, float],
    calcium_scale: float = 1.0,
) -> float:
    activation = float(np.clip(compartment.concentration("activation"), 0.0, 1.0))
    single_channel = float(params["p_single_channel"]) * activation * calcium_scale**4
    rate_event = -math.expm1(
        -float(params["k_fusion_vesicle_s"])
        * float(params["vesicle_dt_s"])
        * activation
        * calcium_scale
    )
    return float(np.clip(1.0 - (1.0 - single_channel) * math.exp(-rate_event), 0.0, 1.0))


def vesicle_release_amount(
    compartment: Compartment,
    params: Mapping[str, float],
    calcium_scale: float = 1.0,
) -> float:
    return float(
        compartment.amounts.get("cargo", 0.0)
        * vesicle_release_probability(compartment, params, calcium_scale)
    )


def _weighted_mean_state(states: Sequence[Compartment]) -> Dict[str, float]:
    total_volume = sum(state.volume for state in states)
    pooled: Dict[str, float] = {}
    species = sorted({species for state in states for species in state.amounts})
    for name in species:
        pooled[name] = sum(state.amounts.get(name, 0.0) for state in states) / total_volume
    potential = sum(state.volume * state.psi_mv for state in states) / total_volume
    return {"volume": total_volume, "concentrations": pooled, "psi_mv": potential}


def total_only_projection(states: Sequence[Compartment]) -> List[Compartment]:
    pooled = _weighted_mean_state(states)
    projected: List[Compartment] = []
    for index, state in enumerate(states):
        amounts = {
            species: concentration * state.volume
            for species, concentration in pooled["concentrations"].items()
        }
        projected.append(
            Compartment(
                volume=state.volume,
                amounts=amounts,
                psi_mv=pooled["psi_mv"],
                label=f"total_only_{index + 1}",
            )
        )
    return projected


def no_potential_projection(states: Sequence[Compartment]) -> List[Compartment]:
    pooled_potential = sum(state.volume * state.psi_mv for state in states) / sum(
        state.volume for state in states
    )
    return [
        Compartment(
            volume=state.volume,
            amounts=dict(state.amounts),
            psi_mv=pooled_potential,
            label=f"no_potential_{index + 1}",
        )
        for index, state in enumerate(states)
    ]


def _response_vector(
    states: Sequence[Compartment],
    mode: str,
    params: Mapping[str, float],
    perturbation_scale: float,
) -> List[float]:
    if mode == "mitochondrion":
        return [atp_response(state, params, perturbation_scale) for state in states]
    if mode == "vesicle":
        return [vesicle_release_amount(state, params, perturbation_scale) for state in states]
    raise ValueError(f"unknown mode: {mode}")


def memory_audit(
    states: Sequence[Compartment],
    mode: str,
    params: Mapping[str, float],
    perturbation_scale: float,
) -> Dict[str, object]:
    full = _response_vector(states, mode, params, perturbation_scale)
    total_only = _response_vector(
        total_only_projection(states), mode, params, perturbation_scale
    )
    no_potential = _response_vector(
        no_potential_projection(states), mode, params, perturbation_scale
    )
    full_gap = abs(full[0] - full[1])
    total_gap = abs(total_only[0] - total_only[1])
    potential_gap = abs(no_potential[0] - no_potential[1])
    scale = max(float(np.mean(np.abs(full))), 1.0e-15)
    return {
        "full_response": [float(value) for value in full],
        "total_only_response": [float(value) for value in total_only],
        "no_potential_response": [float(value) for value in no_potential],
        "full_between_daughter_gap": float(full_gap),
        "total_only_between_daughter_gap": float(total_gap),
        "no_potential_between_daughter_gap": float(potential_gap),
        "total_only_gap_fraction_of_full": float(total_gap / full_gap) if full_gap > 0 else 0.0,
        "total_only_mean_absolute_error": float(
            np.mean(np.abs(np.asarray(total_only) - np.asarray(full))) / scale
        ),
        "no_potential_mean_absolute_error": float(
            np.mean(np.abs(np.asarray(no_potential) - np.asarray(full))) / scale
        ),
    }


def _state_record(state: Compartment, mode: str, params: Mapping[str, float]) -> Dict[str, object]:
    if mode == "mitochondrion":
        response = atp_response(state, params, 1.0)
        risk = fusion_probability(state, params)
    else:
        response = vesicle_release_amount(state, params, 1.0)
        risk = vesicle_release_probability(state, params, 1.0)
    return {
        "label": state.label,
        "volume_um3": float(state.volume),
        "psi_mv": float(state.psi_mv),
        "amounts_nmol": {key: float(value) for key, value in state.amounts.items()},
        "concentrations_nmol_um3": {
            key: float(state.concentration(key)) for key in sorted(state.amounts)
        },
        "baseline_response": float(response),
        "next_event_probability": float(risk),
    }


def _mass_error(
    before: Sequence[Compartment], after: Sequence[Compartment]
) -> Dict[str, float]:
    species = sorted(
        {species for state in before for species in state.amounts}
        | {species for state in after for species in state.amounts}
    )
    errors: Dict[str, float] = {}
    for name in species:
        initial = sum(state.amounts.get(name, 0.0) for state in before)
        final = sum(state.amounts.get(name, 0.0) for state in after)
        denominator = max(abs(initial), 1.0e-15)
        errors[name] = float(abs(final - initial) / denominator)
    volume_initial = sum(state.volume for state in before)
    volume_final = sum(state.volume for state in after)
    errors["volume"] = float(abs(volume_final - volume_initial) / max(volume_initial, 1.0e-15))
    return errors


def _potential_mean(states: Sequence[Compartment]) -> float:
    return float(
        sum(state.volume * state.psi_mv for state in states)
        / sum(state.volume for state in states)
    )


def canonical_mitochondrion(
    params: Optional[Mapping[str, float]] = None,
    stochastic: bool = False,
    seed: int = 830,
) -> Dict[str, object]:
    active_params = dict(MITOCHONDRIAL_DEFAULTS if params is None else params)
    initial = initial_mitochondria()
    fused = fuse(initial[0], initial[1], active_params, active_params["tau_post_s"])
    rng = np.random.default_rng(seed)
    daughters = fission(fused, active_params, rng=rng, stochastic=stochastic)
    before = [state.clone() for state in initial]
    errors = _mass_error(before, list(daughters))
    fused_errors = _mass_error([fused.left, fused.right], list(daughters))
    initial_potential_mean = _potential_mean(initial)
    final_potential_mean = _potential_mean(daughters)
    baseline = _response_vector(
        daughters, "mitochondrion", active_params, active_params["substrate_perturbation"]
    )
    risk = [fusion_probability(state, active_params) for state in daughters]
    high_index = int(np.argmin([state.psi_mv for state in daughters]))
    low_index = 1 - high_index
    return {
        "tau_post_s": float(active_params["tau_post_s"]),
        "potential_gap_mv": float(abs(daughters[0].psi_mv - daughters[1].psi_mv)),
        "initial_potential_mean_mv": float(initial_potential_mean),
        "final_potential_mean_mv": float(final_potential_mean),
        "potential_mean_error_mv": float(abs(final_potential_mean - initial_potential_mean)),
        "risk_ratio_hyper_to_depolarized": float(risk[high_index] / max(risk[low_index], 1.0e-15)),
        "response_after_perturbation": [float(value) for value in baseline],
        "response_gap_after_perturbation": float(abs(baseline[0] - baseline[1])),
        "daughters": [
            _state_record(state, "mitochondrion", active_params) for state in daughters
        ],
        "memory_audit": memory_audit(
            daughters, "mitochondrion", active_params, active_params["substrate_perturbation"]
        ),
        "mass_relative_error": errors,
        "fused_to_daughter_mass_relative_error": fused_errors,
        "all_finite": bool(
            all(np.isfinite(value) for state in daughters for value in state.amounts.values())
            and all(np.isfinite(state.psi_mv) for state in daughters)
        ),
    }


def canonical_vesicle(
    params: Optional[Mapping[str, float]] = None,
    stochastic: bool = False,
    seed: int = 831,
) -> Dict[str, object]:
    active_params = dict(VESICLE_DEFAULTS if params is None else params)
    initial = initial_vesicles()
    fused = fuse(initial[0], initial[1], active_params, active_params["tau_post_s"])
    rng = np.random.default_rng(seed)
    daughters = fission(fused, active_params, rng=rng, stochastic=stochastic)
    errors = _mass_error(initial, daughters)
    perturbed = _response_vector(
        daughters, "vesicle", active_params, active_params["calcium_perturbation"]
    )
    baseline = _response_vector(daughters, "vesicle", active_params, 1.0)
    return {
        "tau_post_s": float(active_params["tau_post_s"]),
        "single_channel_reference_probability": float(active_params["p_single_channel"]),
        "fusion_rate_reference_s": float(active_params["k_fusion_vesicle_s"]),
        "daughters": [
            _state_record(state, "vesicle", active_params) for state in daughters
        ],
        "baseline_expected_release_nmol": [float(value) for value in baseline],
        "perturbed_expected_release_nmol": [float(value) for value in perturbed],
        "perturbed_release_gap_nmol": float(abs(perturbed[0] - perturbed[1])),
        "memory_audit": memory_audit(
            daughters, "vesicle", active_params, active_params["calcium_perturbation"]
        ),
        "mass_relative_error": errors,
        "all_finite": bool(
            all(np.isfinite(value) for state in daughters for value in state.amounts.values())
            and all(np.isfinite(state.psi_mv) for state in daughters)
        ),
    }


def analytic_limit() -> Dict[str, float]:
    left = Compartment(
        volume=1.0,
        amounts={"substrate": 2.0, "enzyme": 1.0},
        psi_mv=-180.0,
        label="left",
    )
    right = Compartment(
        volume=1.0,
        amounts={"substrate": 0.0, "enzyme": 1.0},
        psi_mv=-180.0,
        label="right",
    )
    params = dict(MITOCHONDRIAL_DEFAULTS)
    params["tau_substrate_s"] = 1.0e-12
    params["tau_enzyme_s"] = 1.0e-12
    params["tau_psi_s"] = 1.0e-12
    params["q"] = 0.5
    params["delta_psi_mv"] = 0.0
    daughters = fission(fuse(left, right, params, 1.0), params, stochastic=False)
    pooled_substrate = (left.amounts["substrate"] + right.amounts["substrate"]) / (
        left.volume + right.volume
    )
    pooled_enzyme = (left.amounts["enzyme"] + right.amounts["enzyme"]) / (
        left.volume + right.volume
    )
    errors = [
        abs(daughter.concentration("substrate") - pooled_substrate) / pooled_substrate
        for daughter in daughters
    ] + [
        abs(daughter.concentration("enzyme") - pooled_enzyme) / pooled_enzyme
        for daughter in daughters
    ]
    return {
        "pooled_substrate_nmol_um3": float(pooled_substrate),
        "pooled_enzyme_nmol_um3": float(pooled_enzyme),
        "max_relative_concentration_error": float(max(errors)),
        "potential_gap_mv": float(abs(daughters[0].psi_mv - daughters[1].psi_mv)),
    }


def placebo() -> Dict[str, float]:
    state = initial_mitochondria()[0]
    params = dict(MITOCHONDRIAL_DEFAULTS)
    params["delta_psi_mv"] = 0.0
    daughters = fission(
        fuse(state, state.clone("B"), params, params["tau_post_s"]),
        params,
        stochastic=False,
    )
    responses = _response_vector(
        daughters, "mitochondrion", params, params["substrate_perturbation"]
    )
    return {
        "response_gap": float(abs(responses[0] - responses[1])),
        "potential_gap_mv": float(abs(daughters[0].psi_mv - daughters[1].psi_mv)),
    }


def unit_check() -> Dict[str, object]:
    left_volume = 1.0
    right_volume = 1.0
    tau_s = MITOCHONDRIAL_DEFAULTS["tau_substrate_s"]
    coefficient = _exchange_coefficient(left_volume, right_volume, tau_s)
    concentration_difference = 0.7
    flux = coefficient * concentration_difference
    dimensionless = coefficient * tau_s * (1.0 / left_volume + 1.0 / right_volume)
    return {
        "exchange_coefficient_value": float(coefficient),
        "exchange_coefficient_unit": "um^3 s^-1",
        "flux_value_nmol_s": float(flux),
        "flux_unit": "nmol s^-1",
        "dimensionless_exchange_check": float(dimensionless),
        "expected_flux_unit": "nmol s^-1",
        "valid": bool(
            np.isfinite(coefficient)
            and np.isfinite(flux)
            and abs(dimensionless - 1.0) < 1.0e-12
        ),
    }


def _sensitivity_case(params: Mapping[str, float]) -> Dict[str, float]:
    result = canonical_mitochondrion(params=params, stochastic=False)
    audit = result["memory_audit"]
    return {
        "response_gap_nmol_s": float(result["response_gap_after_perturbation"]),
        "potential_gap_mv": float(result["potential_gap_mv"]),
        "risk_ratio": float(result["risk_ratio_hyper_to_depolarized"]),
        "total_only_gap_fraction": float(audit["total_only_gap_fraction_of_full"]),
        "mass_error_max": float(max(result["mass_relative_error"].values())),
    }


def sensitivity() -> Dict[str, object]:
    cases: Dict[str, object] = {}
    baseline_params = dict(MITOCHONDRIAL_DEFAULTS)
    baseline = _sensitivity_case(baseline_params)
    definitions = {
        "tau_bound_group": ["tau_enzyme_s", "tau_damage_s", "tau_opa1_s"],
        "delta_psi": ["delta_psi_mv"],
        "q": ["q"],
    }
    for name, keys in definitions.items():
        rows: Dict[str, Dict[str, float]] = {}
        for factor in (0.5, 1.0, 1.5):
            params = dict(baseline_params)
            for key in keys:
                params[key] = baseline_params[key] * factor
            rows[f"{factor:.1f}x"] = _sensitivity_case(params)
        low = rows["0.5x"]
        high = rows["1.5x"]
        rows["elasticity_endpoint"] = {
            key: float((high[key] - low[key]) / max(abs(baseline[key]), 1.0e-15))
            for key in baseline
        }
        cases[name] = rows
    return {"baseline": baseline, "parameters": cases}


def _criteria(mito: Mapping[str, object], analytic: Mapping[str, float], placebo_result: Mapping[str, float], units: Mapping[str, object]) -> Dict[str, object]:
    mass_errors = list(mito["mass_relative_error"].values()) + list(
        mito["fused_to_daughter_mass_relative_error"].values()
    )
    audit = mito["memory_audit"]
    checks = {
        "mass_conservation": bool(max(mass_errors) <= 1.0e-10),
        "tau_post_50_to_100_s": bool(50.0 <= float(mito["tau_post_s"]) <= 100.0),
        "potential_gap_9_to_15_mv": bool(9.0 <= float(mito["potential_gap_mv"]) <= 15.0),
        "risk_ratio_4_to_8": bool(4.0 <= float(mito["risk_ratio_hyper_to_depolarized"]) <= 8.0),
        "placebo_response_gap": bool(float(placebo_result["response_gap"]) <= 1.0e-10),
        "total_null_reduction_at_least_half": bool(
            float(audit["total_only_gap_fraction_of_full"]) <= 0.5
        ),
        "analytic_limit": bool(float(analytic["max_relative_concentration_error"]) <= 1.0e-10),
        "unit_check": bool(units["valid"]),
        "finite": bool(mito["all_finite"]),
    }
    return {"checks": checks, "pass": bool(all(checks.values()))}


def run_all() -> Dict[str, object]:
    mito = canonical_mitochondrion(stochastic=False)
    mito_stochastic = canonical_mitochondrion(stochastic=True, seed=830)
    vesicle = canonical_vesicle(stochastic=False)
    analytic = analytic_limit()
    placebo_result = placebo()
    units = unit_check()
    sensitivity_result = sensitivity()
    criteria = _criteria(mito, analytic, placebo_result, units)
    prereg_hash = hashlib.sha256(Path(__file__).with_name("PREREG.md").read_bytes()).hexdigest()
    return {
        "id": "BT-HX-Q083",
        "model_status": "first_runnable_mechanistic_model",
        "prereg_sha256": prereg_hash,
        "source_registry": SOURCES,
        "parameter_table": PARAMETER_TABLE,
        "equations": {
            "exchange": "J_s = K_s (c_i,s - c_j,s); K_s = 1/(tau_s (1/V_i + 1/V_j))",
            "potential": "dpsi_i/dt = -K_psi (psi_i - psi_j); fission preserves volume-weighted mean",
            "fission": "soluble n_1=f n_total; lobe species use a*n_left+b*n_right with a*V_left+b*V_right=q*V_total; psi_1=psi_mean-2(1-q)delta_psi, psi_2=psi_mean+2q delta_psi",
            "mitochondrial_response": "J_ATP = V k_ATP c_enzyme c_substrate/(K_M+c_substrate) exp((psi-psi_ref)/psi_scale) exp(-c_damage)-V k_leak",
            "vesicle_response": "p_release combines single-channel trigger and k_fusion over dt; expected_release = cargo*p_release",
        },
        "unit_check": units,
        "analytic_limit": analytic,
        "placebo": placebo_result,
        "mitochondrion": mito,
        "mitochondrion_stochastic_partition": mito_stochastic,
        "vesicle": vesicle,
        "sensitivity": sensitivity_result,
        "frozen_criteria": criteria,
        "limitations": [
            "Nuclear mtDNA copy number and membrane protein topology are not inferred.",
            "The mitochondrial response law is a normalized first-pass model, not a fitted metabolic network.",
            "Vesicle branch is a reusable fusion/partition model; synaptic exocytosis does not imply a literal post-fusion daughter.",
        ],
        "next_resolution_step": [
            "Measure species-specific exchange and partition efficiencies in a controlled fusion/fission pair.",
            "Calibrate potential, OPA1 and damage trajectories with simultaneous time-lapse and perturbation.",
            "For vesicles, distinguish pre-fusion priming memory from post-fission daughter memory experimentally.",
        ],
    }


def main() -> None:
    result = run_all()
    output_path = Path(__file__).with_name("results.json")
    output_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
