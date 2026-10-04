from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

MODEL_EQUATIONS = """
P_v = P_v0 + f_stressed (V - V0) / C_sys
R_sys = R_sys0 (1 + k_tpr_a A + k_tpr_s S - k_tpr_n N)
Q = clip(Q0 + k_pre (P_v - P_v0) + k_con S + q_ext, 0, Q_max)
C_a dP_a/dt = Q - P_a / R_sys
P_r = P_a - P_v
Q_r = clip(P_r / (R_a + m_R R_e), 0, Q_r_max)
P_gc = P_a - Q_r R_a
G = K_f [P_gc - P_bow - pi_gc]_+
f_reabs = clip(f_reabs0 + k_reabs_a A + k_reabs_s S - k_reabs_n N, f_min, f_max)
U = G (1 - f_reabs)
R_a,target = R_a0 exp(k_myo (P_r - P_r0) + k_tgf (G/G0 - 1) + k_aff_s S - k_aff_n N)
tau_a dR_a/dt = R_a,target - R_a
tau_S dS/dt = S* - S
tau_A dA/dt = A* - A
tau_N dN/dt = N* - N
dV/dt = I_in(t) - U - L_other
dI_cum/dt = I_in(t)
dU_cum/dt = U
dL_cum/dt = L_other
""".strip()

STATE_NAMES = (
    "ecf_volume_ml",
    "arterial_pressure_mmHg",
    "afferent_resistance_mmHg_s_ml",
    "sympathetic_activity",
    "angiotensin_ii_activity",
    "natriuretic_activity",
    "cumulative_input_ml",
    "cumulative_urine_ml",
    "cumulative_other_loss_ml",
)

SENSITIVITY_PARAMETERS = (
    "preload_gain_ml_s_mmHg",
    "myogenic_gain_per_mmHg",
    "raas_time_constant_s",
)

PROTOCOLS = (
    "saline_load",
    "afterload_step",
    "renal_resistance_step",
    "combined",
    "baseline",
)


@dataclass(frozen=True)
class Parameters:
    arterial_pressure0_mmHg: float = 100.0
    cardiac_output0_l_min: float = 4.9
    arterial_compliance_ml_mmHg: float = 1.4
    venous_pressure0_mmHg: float = 5.0
    systemic_compliance_ml_mmHg: float = 68.0
    ecf_volume0_ml: float = 17690.0
    stressed_volume_fraction: float = 0.30
    gfr0_ml_min: float = 122.8
    renal_blood_flow0_ml_min: float = 1076.8
    glomerular_pressure0_mmHg: float = 60.0
    bowman_pressure_mmHg: float = 18.0
    glomerular_oncotic_pressure_mmHg: float = 32.0
    reabsorption_fraction0: float = 0.992
    other_loss_ml_s: float = 0.005
    preload_gain_ml_s_mmHg: float = 8.0
    contractility_gain_ml_s: float = 6.0
    external_contractility_ml_s: float = 0.0
    tpr_angiotensin_gain: float = 0.8
    tpr_sympathetic_gain: float = 0.3
    tpr_natriuretic_gain: float = 0.4
    myogenic_gain_per_mmHg: float = 0.015
    tgf_gain: float = 0.8
    afferent_sympathetic_gain: float = 0.08
    afferent_natriuretic_gain: float = 0.05
    renal_time_constant_s: float = 20.0
    sympathetic_time_constant_s: float = 5.0
    natriuretic_time_constant_s: float = 30.0
    raas_time_constant_s: float = 300.0
    sympathetic_pressure_gain: float = 0.015
    sympathetic_volume_gain: float = 0.8
    angiotensin_pressure_gain: float = 0.025
    angiotensin_volume_gain: float = 0.9
    angiotensin_sympathetic_gain: float = 0.5
    natriuretic_volume_gain: float = 0.8
    natriuretic_pressure_gain: float = 0.02
    reabsorption_angiotensin_gain: float = 0.015
    reabsorption_sympathetic_gain: float = 0.008
    reabsorption_natriuretic_gain: float = 0.010
    volume_scale_ml: float = 500.0
    reabsorption_min: float = 0.97
    reabsorption_max: float = 0.999
    afferent_factor_min: float = 0.30
    afferent_factor_max: float = 3.50
    pressure_floor_mmHg: float = 40.0
    pressure_ceiling_mmHg: float = 200.0
    maximum_cardiac_output_ml_s: float = 150.0
    maximum_renal_blood_flow_ml_s: float = 40.0
    maximum_gfr_ml_s: float = 5.0
    bolus_volume_ml: float = 3000.0
    bolus_duration_s: float = 10800.0
    perturbation_start_s: float = 300.0
    simulation_duration_s: float = 14400.0
    integration_rtol: float = 1.0e-8
    integration_atol: float = 1.0e-9
    maximum_step_s: float = 5.0
    output_points: int = 241
    renal_autoregulation_enabled: bool = True
    neural_hormonal_enabled: bool = True
    volume_pressure_enabled: bool = True
    filtration_enabled: bool = True

    def __post_init__(self) -> None:
        positive_fields = (
            "cardiac_output0_l_min",
            "arterial_compliance_ml_mmHg",
            "systemic_compliance_ml_mmHg",
            "ecf_volume0_ml",
            "gfr0_ml_min",
            "renal_blood_flow0_ml_min",
            "arterial_pressure0_mmHg",
            "bolus_volume_ml",
            "bolus_duration_s",
            "simulation_duration_s",
            "maximum_step_s",
            "volume_scale_ml",
        )
        for name in positive_fields:
            if getattr(self, name) <= 0.0:
                raise ValueError(f"{name} must be positive")
        if not 0.0 < self.reabsorption_fraction0 < 1.0:
            raise ValueError("reabsorption_fraction0 must be between zero and one")
        if not self.reabsorption_min <= self.reabsorption_fraction0 <= self.reabsorption_max:
            raise ValueError("reabsorption fraction is outside its bounds")
        if not self.reabsorption_min < self.reabsorption_max < 1.0:
            raise ValueError("reabsorption bounds are invalid")
        if not 0.0 < self.stressed_volume_fraction <= 1.0:
            raise ValueError("stressed_volume_fraction must be in (0, 1]")
        if not 0.0 < self.venous_pressure0_mmHg < self.arterial_pressure0_mmHg:
            raise ValueError("venous and arterial operating pressures are invalid")
        if self.glomerular_pressure0_mmHg <= self.bowman_pressure_mmHg + self.glomerular_oncotic_pressure_mmHg:
            raise ValueError("baseline glomerular net pressure must be positive")
        if self.perturbation_start_s < 0.0 or self.perturbation_start_s >= self.simulation_duration_s:
            raise ValueError("perturbation_start_s must lie inside the simulation")
        if self.bolus_duration_s <= 0.0 or self.bolus_duration_s > self.simulation_duration_s:
            raise ValueError("bolus_duration_s is invalid")
        if self.output_points < 3:
            raise ValueError("output_points must be at least three")
        if self.pressure_floor_mmHg >= self.pressure_ceiling_mmHg:
            raise ValueError("pressure bounds are invalid")
        if self.maximum_cardiac_output_ml_s <= 0.0 or self.maximum_renal_blood_flow_ml_s <= 0.0:
            raise ValueError("flow bounds must be positive")
        if self.maximum_gfr_ml_s <= 0.0:
            raise ValueError("GFR bound must be positive")

    @property
    def cardiac_output0_ml_s(self) -> float:
        return self.cardiac_output0_l_min * 1000.0 / 60.0

    @property
    def gfr0_ml_s(self) -> float:
        return self.gfr0_ml_min / 60.0

    @property
    def renal_blood_flow0_ml_s(self) -> float:
        return self.renal_blood_flow0_ml_min / 60.0

    @property
    def renal_pressure0_mmHg(self) -> float:
        return self.arterial_pressure0_mmHg - self.venous_pressure0_mmHg

    @property
    def efferent_resistance0_mmHg_s_ml(self) -> float:
        return (self.glomerular_pressure0_mmHg - self.venous_pressure0_mmHg) / self.renal_blood_flow0_ml_s

    @property
    def afferent_resistance0_mmHg_s_ml(self) -> float:
        return self.renal_pressure0_mmHg / self.renal_blood_flow0_ml_s - self.efferent_resistance0_mmHg_s_ml

    @property
    def filtration_coefficient_ml_s_mmHg(self) -> float:
        net_pressure = self.glomerular_pressure0_mmHg - self.bowman_pressure_mmHg - self.glomerular_oncotic_pressure_mmHg
        return self.gfr0_ml_s / net_pressure

    @property
    def urine0_ml_s(self) -> float:
        if not self.filtration_enabled:
            return 0.0
        return self.gfr0_ml_s * (1.0 - self.reabsorption_fraction0)

    @property
    def baseline_input_ml_s(self) -> float:
        return self.urine0_ml_s + self.other_loss_ml_s

    @property
    def systemic_resistance0_mmHg_s_ml(self) -> float:
        return self.arterial_pressure0_mmHg / self.cardiac_output0_ml_s

    def validate_derived_operating_point(self) -> dict[str, float]:
        q_r = self.renal_blood_flow0_ml_s
        p_gc = self.venous_pressure0_mmHg + q_r * self.efferent_resistance0_mmHg_s_ml
        q = self.cardiac_output0_ml_s
        return {
            "cardiac_output_ml_s": q,
            "renal_blood_flow_ml_s": q_r,
            "glomerular_pressure_mmHg": p_gc,
            "filtration_ml_s": self.filtration_coefficient_ml_s_mmHg * max(p_gc - self.bowman_pressure_mmHg - self.glomerular_oncotic_pressure_mmHg, 0.0) if self.filtration_enabled else 0.0,
            "urine_ml_s": self.urine0_ml_s,
        }


PARAMETER_UNITS = {
    "arterial_pressure0_mmHg": "mmHg",
    "cardiac_output0_l_min": "L/min",
    "arterial_compliance_ml_mmHg": "mL/mmHg",
    "venous_pressure0_mmHg": "mmHg",
    "systemic_compliance_ml_mmHg": "mL/mmHg",
    "ecf_volume0_ml": "mL",
    "stressed_volume_fraction": "1",
    "gfr0_ml_min": "mL/min",
    "renal_blood_flow0_ml_min": "mL/min",
    "glomerular_pressure0_mmHg": "mmHg",
    "bowman_pressure_mmHg": "mmHg",
    "glomerular_oncotic_pressure_mmHg": "mmHg",
    "reabsorption_fraction0": "1",
    "other_loss_ml_s": "mL/s",
    "preload_gain_ml_s_mmHg": "mL/(s mmHg)",
    "contractility_gain_ml_s": "mL/s per activity unit",
    "external_contractility_ml_s": "mL/s",
    "tpr_angiotensin_gain": "1/activity unit",
    "tpr_sympathetic_gain": "1/activity unit",
    "tpr_natriuretic_gain": "1/activity unit",
    "myogenic_gain_per_mmHg": "1/mmHg",
    "tgf_gain": "1/(GFR fraction)",
    "afferent_sympathetic_gain": "1/activity unit",
    "afferent_natriuretic_gain": "1/activity unit",
    "renal_time_constant_s": "s",
    "sympathetic_time_constant_s": "s",
    "natriuretic_time_constant_s": "s",
    "raas_time_constant_s": "s",
    "sympathetic_pressure_gain": "1/(mmHg activity unit)",
    "sympathetic_volume_gain": "1/activity unit",
    "angiotensin_pressure_gain": "1/(mmHg activity unit)",
    "angiotensin_volume_gain": "1/activity unit",
    "angiotensin_sympathetic_gain": "1/activity unit",
    "natriuretic_volume_gain": "1/activity unit",
    "natriuretic_pressure_gain": "1/(mmHg activity unit)",
    "reabsorption_angiotensin_gain": "1/activity unit",
    "reabsorption_sympathetic_gain": "1/activity unit",
    "reabsorption_natriuretic_gain": "1/activity unit",
    "volume_scale_ml": "mL",
    "reabsorption_min": "1",
    "reabsorption_max": "1",
    "afferent_factor_min": "1",
    "afferent_factor_max": "1",
    "pressure_floor_mmHg": "mmHg",
    "pressure_ceiling_mmHg": "mmHg",
    "maximum_cardiac_output_ml_s": "mL/s",
    "maximum_renal_blood_flow_ml_s": "mL/s",
    "maximum_gfr_ml_s": "mL/s",
    "bolus_volume_ml": "mL",
    "bolus_duration_s": "s",
    "perturbation_start_s": "s",
    "simulation_duration_s": "s",
    "integration_rtol": "1",
    "integration_atol": "state units",
    "maximum_step_s": "s",
    "output_points": "1",
    "renal_autoregulation_enabled": "boolean",
    "neural_hormonal_enabled": "boolean",
    "volume_pressure_enabled": "boolean",
    "filtration_enabled": "boolean",
}

PARAMETER_SOURCES = {
    "arterial_pressure0_mmHg": "Bikia et al. 2024, Table 2, in-vivo mean blood pressure 100 +/- 12 mmHg",
    "cardiac_output0_l_min": "Bikia et al. 2024, Table 2, in-vivo cardiac output 4.9 +/- 1.2 L/min",
    "arterial_compliance_ml_mmHg": "Bikia et al. 2024, Table 2, in-vivo total arterial compliance 1.4 +/- 0.5 mL/mmHg",
    "venous_pressure0_mmHg": "controlled operating-point assumption; not a Q017 measurement",
    "systemic_compliance_ml_mmHg": "Maas et al. 2012, systemic compliance approximately 0.97 mL/mmHg/kg for 70 kg, rounded",
    "ecf_volume0_ml": "Bhave and Neilson 2011, Table 3, 73-kg man ECF; review-tier initialization",
    "stressed_volume_fraction": "Magder value reused from existing venous_return.py; reduced-model assumption",
    "gfr0_ml_min": "Davies and Shock 1950, Table II, 20-29-year group, 122.8 mL/min/1.73m2",
    "renal_blood_flow0_ml_min": "Davies and Shock 1950, Table II, 20-29-year group, 1076.8 mL/min/1.73m2",
    "glomerular_pressure0_mmHg": "controlled local glomerular working point from existing nephron model",
    "bowman_pressure_mmHg": "controlled local pressure assumption from existing nephron model",
    "glomerular_oncotic_pressure_mmHg": "controlled local oncotic-pressure assumption from existing nephron model",
    "reabsorption_fraction0": "controlled model assumption; approximately 1 mL/min baseline urine",
    "other_loss_ml_s": "controlled nonrenal loss assumption",
    "preload_gain_ml_s_mmHg": "controlled reduced cardiac preload coefficient; not fitted to the primary reference",
    "contractility_gain_ml_s": "controlled reduced sympathetic pump coefficient",
    "external_contractility_ml_s": "zero in primary arm; reserved for external perturbation",
    "tpr_angiotensin_gain": "controlled reduced RAAS-to-TPR coefficient",
    "tpr_sympathetic_gain": "controlled reduced sympathetic-to-TPR coefficient",
    "tpr_natriuretic_gain": "controlled reduced natriuretic-to-TPR coefficient",
    "myogenic_gain_per_mmHg": "controlled myogenic feedback gain",
    "tgf_gain": "controlled TGF gain; direction from Ito and Abe 1996",
    "afferent_sympathetic_gain": "controlled sympathetic afferent-tone coefficient",
    "afferent_natriuretic_gain": "controlled natriuretic afferent-tone coefficient",
    "renal_time_constant_s": "controlled local feedback time constant",
    "sympathetic_time_constant_s": "controlled fast neural time constant",
    "natriuretic_time_constant_s": "controlled natriuretic time constant",
    "raas_time_constant_s": "controlled slow hormonal time constant",
    "sympathetic_pressure_gain": "controlled baroreflex pressure coefficient",
    "sympathetic_volume_gain": "controlled volume-to-sympathetic coefficient",
    "angiotensin_pressure_gain": "controlled low-pressure-to-AngII coefficient",
    "angiotensin_volume_gain": "controlled low-volume-to-AngII coefficient",
    "angiotensin_sympathetic_gain": "controlled sympathetic-to-AngII coefficient",
    "natriuretic_volume_gain": "controlled high-volume-to-natriuretic coefficient",
    "natriuretic_pressure_gain": "controlled high-pressure-to-natriuretic coefficient",
    "reabsorption_angiotensin_gain": "controlled AngII-to-reabsorption coefficient",
    "reabsorption_sympathetic_gain": "controlled sympathetic-to-reabsorption coefficient",
    "reabsorption_natriuretic_gain": "controlled natriuretic-to-reabsorption coefficient",
    "volume_scale_ml": "controlled normalization scale",
    "reabsorption_min": "numerical physiological bound",
    "reabsorption_max": "numerical physiological bound",
    "afferent_factor_min": "numerical resistance bound",
    "afferent_factor_max": "numerical resistance bound",
    "pressure_floor_mmHg": "numerical safety bound",
    "pressure_ceiling_mmHg": "numerical safety bound",
    "maximum_cardiac_output_ml_s": "numerical safety bound",
    "maximum_renal_blood_flow_ml_s": "numerical safety bound",
    "maximum_gfr_ml_s": "numerical safety bound",
    "bolus_volume_ml": "Kumar et al. 2004, Methods, 3 L saline; synthetic protocol dose",
    "bolus_duration_s": "Kumar et al. 2004, Methods, 3 h; synthetic protocol duration",
    "perturbation_start_s": "controlled protocol timing",
    "simulation_duration_s": "controlled protocol horizon",
    "integration_rtol": "numerical setting",
    "integration_atol": "numerical setting",
    "maximum_step_s": "numerical setting",
    "output_points": "numerical output resolution",
    "renal_autoregulation_enabled": "ablation switch",
    "neural_hormonal_enabled": "ablation switch",
    "volume_pressure_enabled": "ablation switch",
    "filtration_enabled": "zero-filtration control switch",
}


def parameter_table(parameters: Parameters | None = None) -> list[dict[str, Any]]:
    p = parameters or Parameters()
    rows = []
    for name, value in asdict(p).items():
        rows.append(
            {
                "name": name,
                "value": value,
                "unit": PARAMETER_UNITS.get(name, "declared model unit"),
                "source_or_assumption": PARAMETER_SOURCES.get(name, "controlled model assumption"),
            }
        )
    derived = p.validate_derived_operating_point()
    derived_rows = (
        ("cardiac_output0_ml_s", derived["cardiac_output_ml_s"], "mL/s", "derived from L/min reference"),
        ("renal_blood_flow0_ml_s", derived["renal_blood_flow_ml_s"], "mL/s", "derived from mL/min reference"),
        ("glomerular_pressure0_mmHg", derived["glomerular_pressure_mmHg"], "mmHg", "derived from resistor divider"),
        ("filtration_coefficient_ml_s_mmHg", p.filtration_coefficient_ml_s_mmHg, "mL/(s mmHg)", "derived from Starling balance"),
        ("urine0_ml_s", derived["urine_ml_s"], "mL/s", "derived from GFR and reabsorption fraction"),
        ("systemic_resistance0_mmHg_s_ml", p.systemic_resistance0_mmHg_s_ml, "mmHg s/mL", "derived from MAP/CO"),
    )
    for name, value, unit, source in derived_rows:
        rows.append({"name": name, "value": value, "unit": unit, "source_or_assumption": source})
    return rows


def initial_state(parameters: Parameters) -> np.ndarray:
    return np.array(
        [
            parameters.ecf_volume0_ml,
            parameters.arterial_pressure0_mmHg,
            parameters.afferent_resistance0_mmHg_s_ml,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )


def infusion_rate_ml_s(time_s: float, parameters: Parameters, protocol: str) -> float:
    rate = parameters.baseline_input_ml_s
    if protocol in {"saline_load", "combined"}:
        start = parameters.perturbation_start_s
        end = start + parameters.bolus_duration_s
        if start <= time_s < end:
            rate += parameters.bolus_volume_ml / parameters.bolus_duration_s
    return rate


def _external_multipliers(time_s: float, parameters: Parameters, protocol: str) -> tuple[float, float, float]:
    active = time_s >= parameters.perturbation_start_s
    afterload = 1.20 if active and protocol in {"afterload_step", "combined"} else 1.0
    renal_resistance = 1.25 if active and protocol in {"renal_resistance_step", "combined"} else 1.0
    contractility = parameters.external_contractility_ml_s
    return afterload, renal_resistance, contractility


def derived_quantities(time_s: float, state: np.ndarray, parameters: Parameters, protocol: str) -> dict[str, float]:
    volume, arterial_pressure, afferent_resistance, sympathetic, angiotensin, natriuretic = state[:6]
    arterial_pressure_eff = float(np.clip(arterial_pressure, parameters.pressure_floor_mmHg, parameters.pressure_ceiling_mmHg))
    stressed_fraction = parameters.stressed_volume_fraction if parameters.volume_pressure_enabled else 0.0
    venous_pressure = parameters.venous_pressure0_mmHg + stressed_fraction * (volume - parameters.ecf_volume0_ml) / parameters.systemic_compliance_ml_mmHg
    afterload_factor, renal_resistance_factor, external_contractility = _external_multipliers(time_s, parameters, protocol)
    if parameters.neural_hormonal_enabled:
        tpr_activity = (
            1.0
            + parameters.tpr_angiotensin_gain * angiotensin
            + parameters.tpr_sympathetic_gain * sympathetic
            - parameters.tpr_natriuretic_gain * natriuretic
        )
    else:
        tpr_activity = 1.0
    systemic_resistance = parameters.systemic_resistance0_mmHg_s_ml * max(0.25, tpr_activity) * afterload_factor
    cardiac_output = parameters.cardiac_output0_ml_s + parameters.preload_gain_ml_s_mmHg * (venous_pressure - parameters.venous_pressure0_mmHg) + parameters.contractility_gain_ml_s * sympathetic + external_contractility
    cardiac_output = float(np.clip(cardiac_output, 0.0, parameters.maximum_cardiac_output_ml_s))
    renal_pressure = arterial_pressure_eff - venous_pressure
    total_renal_resistance = max(afferent_resistance, 1.0e-9) + parameters.efferent_resistance0_mmHg_s_ml * renal_resistance_factor
    renal_blood_flow = float(np.clip(renal_pressure / total_renal_resistance, 0.0, parameters.maximum_renal_blood_flow_ml_s))
    glomerular_pressure = arterial_pressure_eff - renal_blood_flow * max(afferent_resistance, 0.0)
    if parameters.filtration_enabled:
        net_pressure = max(glomerular_pressure - parameters.bowman_pressure_mmHg - parameters.glomerular_oncotic_pressure_mmHg, 0.0)
        gfr = parameters.filtration_coefficient_ml_s_mmHg * net_pressure
        gfr = float(np.clip(gfr, 0.0, parameters.maximum_gfr_ml_s))
    else:
        gfr = 0.0
    if parameters.neural_hormonal_enabled:
        reabsorption = parameters.reabsorption_fraction0 + parameters.reabsorption_angiotensin_gain * angiotensin + parameters.reabsorption_sympathetic_gain * sympathetic - parameters.reabsorption_natriuretic_gain * natriuretic
        reabsorption = float(np.clip(reabsorption, parameters.reabsorption_min, parameters.reabsorption_max))
    else:
        reabsorption = parameters.reabsorption_fraction0
    urine = max(0.0, gfr * (1.0 - reabsorption))
    volume_deviation = (parameters.ecf_volume0_ml - volume) / parameters.volume_scale_ml
    if parameters.neural_hormonal_enabled:
        sympathetic_target = float(np.clip(-parameters.sympathetic_pressure_gain * (arterial_pressure_eff - parameters.arterial_pressure0_mmHg) + parameters.sympathetic_volume_gain * volume_deviation, -1.5, 1.5))
        angiotensin_target = float(np.clip(parameters.angiotensin_pressure_gain * (parameters.arterial_pressure0_mmHg - arterial_pressure_eff) + parameters.angiotensin_volume_gain * volume_deviation + parameters.angiotensin_sympathetic_gain * sympathetic, -0.5, 2.0))
        natriuretic_target = float(np.clip(parameters.natriuretic_volume_gain * (-volume_deviation) + parameters.natriuretic_pressure_gain * (arterial_pressure_eff - parameters.arterial_pressure0_mmHg), -0.5, 2.0))
    else:
        sympathetic_target = 0.0
        angiotensin_target = 0.0
        natriuretic_target = 0.0
    if parameters.renal_autoregulation_enabled:
        gfr_ratio = gfr / parameters.gfr0_ml_s if parameters.gfr0_ml_s > 0.0 else 0.0
        log_gain = parameters.myogenic_gain_per_mmHg * (renal_pressure - parameters.renal_pressure0_mmHg) + parameters.tgf_gain * (gfr_ratio - 1.0)
        if parameters.neural_hormonal_enabled:
            log_gain += parameters.afferent_sympathetic_gain * sympathetic - parameters.afferent_natriuretic_gain * natriuretic
    else:
        log_gain = 0.0
    factor_min = parameters.afferent_factor_min
    factor_max = parameters.afferent_factor_max
    afferent_factor = float(np.clip(math.exp(log_gain), factor_min, factor_max))
    afferent_target = parameters.afferent_resistance0_mmHg_s_ml * afferent_factor
    return {
        "ecf_volume_ml": float(volume),
        "arterial_pressure_mmHg": float(arterial_pressure),
        "venous_pressure_mmHg": float(venous_pressure),
        "systemic_resistance_mmHg_s_ml": float(systemic_resistance),
        "cardiac_output_ml_s": cardiac_output,
        "renal_pressure_mmHg": float(renal_pressure),
        "renal_blood_flow_ml_s": renal_blood_flow,
        "glomerular_pressure_mmHg": float(glomerular_pressure),
        "gfr_ml_s": gfr,
        "reabsorption_fraction": reabsorption,
        "urine_ml_s": urine,
        "input_ml_s": infusion_rate_ml_s(time_s, parameters, protocol),
        "other_loss_ml_s": parameters.other_loss_ml_s,
        "sympathetic_target": sympathetic_target,
        "angiotensin_target": angiotensin_target,
        "natriuretic_target": natriuretic_target,
        "afferent_target_mmHg_s_ml": float(afferent_target),
        "afferent_factor": afferent_factor,
    }


def rhs(time_s: float, state: np.ndarray, parameters: Parameters, protocol: str) -> np.ndarray:
    volume, arterial_pressure, afferent_resistance, sympathetic, angiotensin, natriuretic = state[:6]
    d = derived_quantities(time_s, state, parameters, protocol)
    d_volume = d["input_ml_s"] - d["urine_ml_s"] - d["other_loss_ml_s"]
    d_pressure = (d["cardiac_output_ml_s"] - arterial_pressure / d["systemic_resistance_mmHg_s_ml"]) / parameters.arterial_compliance_ml_mmHg
    d_afferent = (d["afferent_target_mmHg_s_ml"] - afferent_resistance) / parameters.renal_time_constant_s
    d_sympathetic = (d["sympathetic_target"] - sympathetic) / parameters.sympathetic_time_constant_s
    d_angiotensin = (d["angiotensin_target"] - angiotensin) / parameters.raas_time_constant_s
    d_natriuretic = (d["natriuretic_target"] - natriuretic) / parameters.natriuretic_time_constant_s
    return np.array(
        [
            d_volume,
            d_pressure,
            d_afferent,
            d_sympathetic,
            d_angiotensin,
            d_natriuretic,
            d["input_ml_s"],
            d["urine_ml_s"],
            d["other_loss_ml_s"],
        ],
        dtype=float,
    )


def _time_grid(parameters: Parameters, duration_s: float) -> np.ndarray:
    base = np.linspace(0.0, duration_s, parameters.output_points)
    keys = [0.0, duration_s, parameters.perturbation_start_s]
    if parameters.bolus_duration_s > 0.0:
        keys.append(min(duration_s, parameters.perturbation_start_s + parameters.bolus_duration_s))
    return np.unique(np.clip(np.concatenate((base, np.asarray(keys, dtype=float))), 0.0, duration_s))


def simulate(parameters: Parameters | None = None, protocol: str = "saline_load", duration_s: float | None = None) -> dict[str, Any]:
    p = parameters or Parameters()
    if protocol not in PROTOCOLS:
        raise ValueError(f"unknown protocol: {protocol}")
    duration = p.simulation_duration_s if duration_s is None else float(duration_s)
    if duration <= 0.0:
        raise ValueError("duration_s must be positive")
    times = _time_grid(p, duration)
    solution = solve_ivp(
        lambda time_s, state: rhs(time_s, state, p, protocol),
        (0.0, duration),
        initial_state(p),
        t_eval=times,
        method="DOP853",
        rtol=p.integration_rtol,
        atol=p.integration_atol,
        max_step=p.maximum_step_s,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    states = solution.y.T
    output_rows = [derived_quantities(time_s, state, p, protocol) for time_s, state in zip(solution.t, states)]
    output_names = tuple(output_rows[0].keys())
    outputs = {name: np.asarray([row[name] for row in output_rows], dtype=float) for name in output_names}
    mass_balance_error = states[:, 0] - p.ecf_volume0_ml - states[:, 6] + states[:, 7] + states[:, 8]
    return {
        "time_s": np.asarray(solution.t, dtype=float),
        "states": {name: states[:, index] for index, name in enumerate(STATE_NAMES)},
        "outputs": outputs,
        "mass_balance_error_ml": mass_balance_error,
        "success": bool(solution.success),
        "message": str(solution.message),
    }


def _index_at_time(times: np.ndarray, target_s: float) -> int:
    return int(np.argmin(np.abs(times - target_s)))


def _point(outputs: dict[str, np.ndarray], index: int) -> dict[str, float]:
    point = {name: float(values[index]) for name, values in outputs.items()}
    point["cardiac_output_l_min"] = point["cardiac_output_ml_s"] * 60.0 / 1000.0
    point["renal_blood_flow_l_min"] = point["renal_blood_flow_ml_s"] * 60.0 / 1000.0
    point["gfr_ml_min"] = point["gfr_ml_s"] * 60.0
    point["urine_ml_min"] = point["urine_ml_s"] * 60.0
    point["ecf_volume_l"] = point["ecf_volume_ml"] / 1000.0
    return point


def summarize(simulation: dict[str, Any], parameters: Parameters | None = None, protocol: str = "saline_load") -> dict[str, Any]:
    p = parameters or Parameters()
    times = np.asarray(simulation["time_s"], dtype=float)
    outputs = simulation["outputs"]
    states = simulation["states"]
    end_time = min(float(times[-1]), p.perturbation_start_s + p.bolus_duration_s) if protocol in {"saline_load", "combined"} else float(times[-1])
    baseline_index = 0
    end_index = _index_at_time(times, end_time)
    final_index = len(times) - 1
    baseline = _point(outputs, baseline_index)
    endpoint = _point(outputs, end_index)
    final = _point(outputs, final_index)
    peak_names = {
        "arterial_pressure_mmHg": "arterial_pressure_mmHg",
        "cardiac_output_l_min": "cardiac_output_l_min",
        "renal_blood_flow_l_min": "renal_blood_flow_l_min",
        "gfr_ml_min": "gfr_ml_min",
        "urine_ml_min": "urine_ml_min",
        "ecf_volume_l": "ecf_volume_l",
    }
    peaks = {}
    for label, name in peak_names.items():
        values = np.asarray([_point(outputs, i)[name] for i in range(len(times))])
        max_index = int(np.argmax(values))
        min_index = int(np.argmin(values))
        peaks[label] = {
            "maximum": float(values[max_index]),
            "maximum_time_s": float(times[max_index]),
            "minimum": float(values[min_index]),
            "minimum_time_s": float(times[min_index]),
        }
    mass_error = np.asarray(simulation["mass_balance_error_ml"], dtype=float)
    return {
        "protocol": protocol,
        "baseline": baseline,
        "end_of_perturbation": endpoint,
        "final": final,
        "relative_cardiac_output_change": (endpoint["cardiac_output_l_min"] - baseline["cardiac_output_l_min"]) / baseline["cardiac_output_l_min"],
        "peaks": peaks,
        "mass_balance": {
            "maximum_absolute_error_ml": float(np.max(np.abs(mass_error))),
            "final_error_ml": float(mass_error[-1]),
            "relative_final_error": float(abs(mass_error[-1]) / max(1.0, p.ecf_volume0_ml)),
        },
        "state_minima": {name: float(np.min(values)) for name, values in states.items()},
        "state_maxima": {name: float(np.max(values)) for name, values in states.items()},
    }


def analytic_pressure_flow_limit(parameters: Parameters | None = None) -> dict[str, float | bool]:
    p = parameters or Parameters()
    pressure = p.arterial_pressure0_mmHg
    resistance_sum = p.afferent_resistance0_mmHg_s_ml + p.efferent_resistance0_mmHg_s_ml
    epsilon = 1.0e-6
    flow_low = (pressure - epsilon) / resistance_sum
    flow_high = (pressure + epsilon) / resistance_sum
    numerical_derivative = (flow_high - flow_low) / (2.0 * epsilon)
    analytic_derivative = 1.0 / resistance_sum
    return {
        "pressure_mmHg": pressure,
        "resistance_sum_mmHg_s_ml": resistance_sum,
        "analytic_derivative_ml_s_mmHg": analytic_derivative,
        "numerical_derivative_ml_s_mmHg": numerical_derivative,
        "absolute_error_ml_s_mmHg": abs(numerical_derivative - analytic_derivative),
        "passed": abs(numerical_derivative - analytic_derivative) < 1.0e-8,
    }


def renal_pressure_sweep(parameters: Parameters | None = None, pressures_mmHg: np.ndarray | None = None) -> dict[str, Any]:
    p = parameters or Parameters()
    pressures = np.arange(80.0, 180.0 + 0.1, 5.0) if pressures_mmHg is None else np.asarray(pressures_mmHg, dtype=float)
    rows = []
    gfr0 = p.gfr0_ml_s if p.filtration_enabled else 0.0
    for pressure in pressures:
        afferent = p.afferent_resistance0_mmHg_s_ml
        state = initial_state(p)
        for _ in range(500):
            state[1] = pressure
            state[2] = afferent
            d = derived_quantities(0.0, state, p, "baseline")
            afferent = float(np.clip(0.5 * afferent + 0.5 * d["afferent_target_mmHg_s_ml"], p.afferent_factor_min * p.afferent_resistance0_mmHg_s_ml, p.afferent_factor_max * p.afferent_resistance0_mmHg_s_ml))
        state[1] = pressure
        state[2] = afferent
        d = derived_quantities(0.0, state, p, "baseline")
        rows.append(
            {
                "arterial_pressure_mmHg": float(pressure),
                "renal_pressure_mmHg": d["renal_pressure_mmHg"],
                "afferent_resistance_mmHg_s_ml": float(afferent),
                "renal_blood_flow_ml_s": d["renal_blood_flow_ml_s"],
                "gfr_ml_s": d["gfr_ml_s"],
                "gfr_ratio": d["gfr_ml_s"] / gfr0 if gfr0 > 0.0 else 0.0,
            }
        )
    max_deviation = max(abs(row["gfr_ratio"] - 1.0) for row in rows) if rows else float("inf")
    return {
        "pressures_mmHg": pressures.tolist(),
        "rows": rows,
        "maximum_absolute_gfr_fractional_deviation": max_deviation,
        "criterion_fractional_deviation": 0.20,
        "criterion_passed": max_deviation <= 0.20,
    }


def dimension_analysis(parameters: Parameters | None = None) -> dict[str, Any]:
    p = parameters or Parameters()
    checks = {
        "cardiac_output_conversion": "L/min * 1000 mL/L / 60 s/min = mL/s",
        "pressure_storage": "(mL/s - mmHg/(mmHg s/mL)) / (mL/mmHg) = mmHg/s",
        "renal_flow": "mmHg / (mmHg s/mL) = mL/s",
        "glomerular_filtration": "mL/(s mmHg) * mmHg = mL/s",
        "urine_flow": "mL/s * dimensionless reabsorption complement = mL/s",
        "volume_balance": "mL/s - mL/s - mL/s = mL/s",
        "resistance_state": "mmHg s/mL / s = mmHg/mL",
        "hormonal_activity": "dimensionless target and state",
    }
    positive = all(
        getattr(p, name) > 0.0
        for name in (
            "cardiac_output0_l_min",
            "arterial_compliance_ml_mmHg",
            "systemic_compliance_ml_mmHg",
            "renal_blood_flow0_ml_min",
        )
    )
    return {"passed": positive, "checks": checks, "operating_point": p.validate_derived_operating_point()}


def _metric_values(summary: dict[str, Any]) -> dict[str, float]:
    endpoint = summary["end_of_perturbation"]
    peaks = summary["peaks"]
    return {
        "relative_cardiac_output_change": summary["relative_cardiac_output_change"],
        "end_arterial_pressure_mmHg": endpoint["arterial_pressure_mmHg"],
        "end_cardiac_output_l_min": endpoint["cardiac_output_l_min"],
        "end_renal_blood_flow_l_min": endpoint["renal_blood_flow_l_min"],
        "end_gfr_ml_min": endpoint["gfr_ml_min"],
        "end_urine_ml_min": endpoint["urine_ml_min"],
        "end_ecf_volume_l": endpoint["ecf_volume_l"],
        "peak_arterial_pressure_mmHg": peaks["arterial_pressure_mmHg"]["maximum"],
        "peak_gfr_ml_min": peaks["gfr_ml_min"]["maximum"],
    }


def sensitivity_analysis(parameters: Parameters | None = None) -> dict[str, Any]:
    p = parameters or Parameters()
    base_summary = summarize(simulate(p, "saline_load"), p, "saline_load")
    base_metrics = _metric_values(base_summary)
    rows = []
    for name in SENSITIVITY_PARAMETERS:
        base_value = float(getattr(p, name))
        for label, factor in (("minus50", 0.5), ("plus50", 1.5)):
            value = base_value * factor
            variant = replace(p, **{name: value})
            summary = summarize(simulate(variant, "saline_load"), variant, "saline_load")
            metrics = _metric_values(summary)
            changes = {}
            for metric, base_metric in base_metrics.items():
                changes[metric] = 100.0 * (metrics[metric] - base_metric) / abs(base_metric) if base_metric != 0.0 else 0.0
            rows.append(
                {
                    "parameter": name,
                    "direction": label,
                    "value": value,
                    "metrics": metrics,
                    "relative_change_percent": changes,
                }
            )
    ranked = sorted(
        rows,
        key=lambda row: max(abs(value) for value in row["relative_change_percent"].values()),
        reverse=True,
    )
    return {"baseline_metrics": base_metrics, "runs": rows, "ranked_by_max_absolute_relative_change": [row["parameter"] + ":" + row["direction"] for row in ranked]}


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def run_experiments(parameters: Parameters | None = None) -> dict[str, Any]:
    p = parameters or Parameters()
    full_simulation = simulate(p, "saline_load")
    full_summary = summarize(full_simulation, p, "saline_load")
    control_parameters = {
        "full": p,
        "no_renal_autoregulation": replace(p, renal_autoregulation_enabled=False),
        "no_neural_hormonal": replace(p, neural_hormonal_enabled=False),
        "no_volume_pressure": replace(p, volume_pressure_enabled=False),
    }
    controls = {}
    for name, control_p in control_parameters.items():
        simulation = simulate(control_p, "saline_load")
        controls[name] = summarize(simulation, control_p, "saline_load")
    controls["afterload_step"] = summarize(simulate(p, "afterload_step"), p, "afterload_step")
    controls["renal_resistance_step"] = summarize(simulate(p, "renal_resistance_step"), p, "renal_resistance_step")
    controls["combined"] = summarize(simulate(p, "combined"), p, "combined")
    baseline_simulation = simulate(p, "baseline", duration_s=600.0)
    baseline_control = summarize(baseline_simulation, p, "baseline")
    zero_filtration_p = replace(p, filtration_enabled=False)
    zero_filtration_simulation = simulate(zero_filtration_p, "saline_load")
    zero_filtration_summary = summarize(zero_filtration_simulation, zero_filtration_p, "saline_load")
    sweep = renal_pressure_sweep(p)
    analytic = analytic_pressure_flow_limit(p)
    dimensions = dimension_analysis(p)
    sensitivity = sensitivity_analysis(p)
    primary_value = full_summary["relative_cardiac_output_change"]
    primary_reference = 0.147
    primary_error = abs(primary_value - primary_reference)
    criteria = {
        "primary_cardiac_output_change_within_0.05": primary_error <= 0.050 and primary_value > 0.0,
        "renal_pressure_sweep_within_20_percent": sweep["criterion_passed"],
        "mass_balance_relative_error_below_1e-9": full_summary["mass_balance"]["relative_final_error"] < 1.0e-9,
        "analytic_pressure_flow_limit": analytic["passed"],
        "dimension_analysis": dimensions["passed"],
        "zero_filtration_has_no_gfr_or_urine": zero_filtration_summary["final"]["gfr_ml_min"] == 0.0 and zero_filtration_summary["final"]["urine_ml_min"] == 0.0,
        "baseline_no_perturbation_stays_within_1_percent": all(
            abs(baseline_control["final"][name] / full_summary["baseline"][name] - 1.0) <= 0.01
            for name in ("arterial_pressure_mmHg", "cardiac_output_l_min", "gfr_ml_min", "urine_ml_min", "ecf_volume_l")
        ),
    }
    criteria["overall_pass"] = all(criteria.values())
    source_status = {
        "davies_shock_1950": {
            "doi": "10.1172/JCI102286",
            "pmid": "15415454",
            "location": "Table II, 20-29-year group",
            "values": {"gfr_ml_min": 122.8, "effective_renal_blood_flow_ml_min": 1076.8},
            "status": "VERIFIERAD primary source; table image/source cross-checked in PMC and existing repository transcription",
        },
        "kumar_2004": {
            "doi": "10.1186/cc2844",
            "pmid": "15153240",
            "location": "Table 1, group 1",
            "values": {"volume_l": 3.0, "duration_h": 3.0, "cardiac_index_change_percent": 14.7, "uncertainty_percentage_points": 2.4, "n": 24},
            "status": "VERIFIERAD primary open full text",
        },
        "bikia_2024": {
            "doi": "10.1038/s41598-024-56137-8",
            "pmid": "38467721",
            "location": "Table 2, in-vivo population",
            "values": {"mean_blood_pressure_mmHg": 100.0, "cardiac_output_l_min": 4.9, "arterial_compliance_ml_mmHg": 1.4, "peripheral_resistance_mmHg_s_ml": 1.3, "time_constant_s": 1.2, "n": 2263},
            "status": "VERIFIERAD primary open full text",
        },
        "maas_2012": {
            "doi": "10.1213/ANE.0b013e31825fb01d",
            "pmid": "22763909",
            "location": "abstract/systemic compliance report",
            "values": {"systemic_compliance_ml_mmHg_per_kg": 0.97},
            "status": "VERIFIERAD secondary physiological anchor used as assumption",
        },
        "bhave_neilson_2011": {
            "doi": "10.1681/ASN.2011080865",
            "pmid": "22034644",
            "location": "Table 3",
            "values": {"73kg_man_ecf_l": 17.69, "functional_ecf_l": 12.82},
            "status": "VERIFIERAD review-tier table; not a transient measurement",
        },
    }
    report = {
        "id": "BT-HX-Q017",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "equations": MODEL_EQUATIONS,
        "state_names": list(STATE_NAMES),
        "parameters": parameter_table(p),
        "derived_operating_point": p.validate_derived_operating_point(),
        "source_status": source_status,
        "primary_reference": {
            "quantity": "relative cardiac-output change after 3 L saline over 3 h",
            "value": primary_reference,
            "unit": "fraction",
            "uncertainty_percentage_points": 2.4,
            "source": "Kumar et al. 2004, DOI 10.1186/cc2844, Table 1",
        },
        "primary_prediction": {
            "relative_cardiac_output_change": primary_value,
            "absolute_error": primary_error,
            "criterion_interval": [0.097, 0.197],
        },
        "full_model": {
            "summary": full_summary,
            "trajectory": {
                "time_s": full_simulation["time_s"],
                "states": full_simulation["states"],
                "outputs": full_simulation["outputs"],
                "mass_balance_error_ml": full_simulation["mass_balance_error_ml"],
            },
        },
        "controls": controls,
        "baseline_no_perturbation": baseline_control,
        "zero_filtration_control": zero_filtration_summary,
        "renal_pressure_sweep": sweep,
        "analytic_limit": analytic,
        "dimension_analysis": dimensions,
        "sensitivity": sensitivity,
        "criteria": criteria,
        "interpretation": {
            "claim_status": "MODEL_RUN_COMPLETE",
            "clinical_claim": "UNKNOWN",
            "reason": "The perturbation is synthetic and no simultaneous patient pressure, flow, volume, GFR and organ-response dataset was available in the bounded task.",
        },
    }
    return _jsonable(report)


def main() -> None:
    report = run_experiments()
    output_path = Path(__file__).with_name("results.json")
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    criteria = report["criteria"]
    print(json.dumps({"output": str(output_path), "criteria": criteria, "primary_prediction": report["primary_prediction"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
