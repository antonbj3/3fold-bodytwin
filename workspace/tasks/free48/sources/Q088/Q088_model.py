from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np


EQUATIONS = r"""
C_t = D_O2 * (1/r) * d/dr(r * dC/dr) - V_O2 * C/(K_O2+C)
-D_O2 * dC/dr|_R = k_transfer * (C_boundary - C_R)
V_O2 = basal_O2 + oxygen_cost * J_ATP_oxidative
J_ATP_oxidative = Vmax_ATP * substrate_factor * mean(C)/(K_O2+mean(C))
dATP/dt = J_ATP_oxidative + J_CK - demand - J_PCr_resynthesis
dPCr/dt = J_PCr_resynthesis - J_CK
dI/dt = leak_rate(activity) - pump_rate * I * (ATP/ATP_rest) * oxygen_factor
F = 0.50*PCr/PCr_rest + 0.30*ATP/ATP_rest + 0.20*(1-I)
"""


@dataclass(frozen=True)
class Parameters:
    radius_m: float = 25e-6
    oxygen_boundary_kpa: float = 13.3
    # EDGE DRAWN 2026-10-03 to Q021, which CALCULATES this boundary instead of assuming it.
    # 13.3 kPa = 99.76 mmHg, which is Q021's HPV REFERENCE and close to its alveolar gas equation
    # (104.225 mmHg) -- thus alveolar oxygen at the cell boundary, with no gradient down to tissue.
    # Q021's own calculated global PAO2 with V/Q heterogeneity: 83.998 mmHg at rest (-15.8 %),
    # 82.407 homogeneous rest (-17.4 %), 76.162 during hard work (-23.7 %).
    # Measured consequence of switching to the resting value 11.1989 kPa: 5 001 of 7 407 fields move, no
    # sign reversals, BUT sensitivity[3] -- work demand 0.3 mM/s, thus HALF the nominal -- goes
    # from t90 functional recovery 46.568 s to CENSORED at the horizon 600 s. The case
    # thus does not recover within the window, and the true value is unknown.
    # The default is left unchanged; the alternatives are here for sweeps and consumers.
    oxygen_boundary_q021_kpa: tuple = (11.198909, 10.986690, 10.154142)
    work_demand_mM_s: float = 0.60
    recovery_demand_mM_s: float = 0.05
    exercise_duration_s: float = 250.0
    recovery_duration_s: float = 600.0
    dt_s: float = 0.25
    sample_stride_s: float = 1.0
    radial_nodes: int = 20
    oxygen_diffusivity_m2_s: float = 1.5e-9
    oxygen_boundary_transfer_m_s: float = 1.0e-5
    oxygen_solubility_mol_m3_kpa: float = 0.030
    oxygen_michaelis_mol_m3: float = 0.05
    atp_rest_mM: float = 8.2
    pcr_rest_mM: float = 32.0
    atp_capacity_mM_s: float = 0.70
    substrate_factor: float = 0.85
    oxygen_cost_mol_o2_mol_atp: float = 0.20
    # Source-backed interval for the above, from the older project's P/O grid. The default is kept
    # at 0.20 so this edge does not silently move any result; the band is reported and swept instead.
    oxygen_cost_band: tuple = (0.1833333333333333, 0.2037037037037037)
    basal_oxygen_mol_m3_s: float = 0.002
    ck_rate_s: float = 0.50
    ck_capacity_mM_s: float = 1.50
    pcr_resynthesis_rate_s: float = 0.05
    ion_leak_rate_s: float = 0.25
    ion_basal_leak_fraction: float = 0.0
    ion_pump_rate_s: float = 0.10
    functional_atp_weight: float = 0.30
    functional_pcr_weight: float = 0.50
    functional_ion_weight: float = 0.20

    def validate(self) -> None:
        positive = {
            "radius_m": self.radius_m,
            "oxygen_diffusivity_m2_s": self.oxygen_diffusivity_m2_s,
            "oxygen_boundary_transfer_m_s": self.oxygen_boundary_transfer_m_s,
            "oxygen_solubility_mol_m3_kpa": self.oxygen_solubility_mol_m3_kpa,
            "oxygen_michaelis_mol_m3": self.oxygen_michaelis_mol_m3,
            "atp_rest_mM": self.atp_rest_mM,
            "pcr_rest_mM": self.pcr_rest_mM,
            "atp_capacity_mM_s": self.atp_capacity_mM_s,
            "ck_capacity_mM_s": self.ck_capacity_mM_s,
            "dt_s": self.dt_s,
            "sample_stride_s": self.sample_stride_s,
            "ion_pump_rate_s": self.ion_pump_rate_s,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not np.isfinite(self.oxygen_boundary_kpa) or self.oxygen_boundary_kpa < 0:
            raise ValueError("oxygen_boundary_kpa must be finite and nonnegative")
        if not np.isfinite(self.work_demand_mM_s) or self.work_demand_mM_s < 0:
            raise ValueError("work_demand_mM_s must be finite and nonnegative")
        if not np.isfinite(self.recovery_demand_mM_s) or self.recovery_demand_mM_s < 0:
            raise ValueError("recovery_demand_mM_s must be finite and nonnegative")
        if self.exercise_duration_s < 0 or self.recovery_duration_s < 0:
            raise ValueError("durations must be nonnegative")
        if self.radial_nodes < 2:
            raise ValueError("radial_nodes must be at least two")
        if not 0 <= self.ion_basal_leak_fraction <= 1:
            raise ValueError("ion_basal_leak_fraction must lie in [0,1]")
        weights = np.array([
            self.functional_atp_weight,
            self.functional_pcr_weight,
            self.functional_ion_weight,
        ])
        if np.any(weights < 0) or not np.isclose(weights.sum(), 1.0):
            raise ValueError("functional weights must be nonnegative and sum to one")


PARAMETER_TABLE = [
    {
        "name": "radius_m",
        "value": 25e-6,
        "unit": "m",
        "source_or_assumption": "assumption; geometry sweep is a mechanistic hypothesis; context in Egginton et al. 2002, DOI 10.1242/jeb.205.6.769",
    },
    {
        "name": "oxygen_boundary_kpa",
        "value": 13.3,
        "unit": "kPa",
        "source_or_assumption": "controlled boundary condition; normoxic arterial-scale assumption, not a measured fiber value",
    },
    {
        "name": "oxygen_diffusivity_m2_s",
        "value": 1.5e-9,
        "unit": "m^2/s",
        "source_or_assumption": "assumption; representative effective tissue value, not fitted",
    },
    {
        "name": "oxygen_boundary_transfer_m_s",
        "value": 1.0e-5,
        "unit": "m/s",
        "source_or_assumption": "assumption; Robin capillary mass-transfer coefficient preventing an infinite-supply boundary",
    },
    {
        "name": "oxygen_solubility_mol_m3_kpa",
        "value": 0.030,
        "unit": "mol m^-3 kPa^-1",
        "source_or_assumption": "physical conversion assumption for dissolved oxygen",
    },
    {
        "name": "oxygen_michaelis_mol_m3",
        "value": 0.05,
        "unit": "mol m^-3",
        "source_or_assumption": "assumption for the oxygen dependence of respiration",
    },
    {
        "name": "atp_rest_mM",
        "value": 8.2,
        "unit": "mM",
        "source_or_assumption": "external anchor: Layec et al. 2013, Methods, DOI 10.1152/japplphysiol.00257.2013",
    },
    {
        "name": "pcr_rest_mM",
        "value": 32.0,
        "unit": "mM",
        "source_or_assumption": "external anchor: Kushmerick et al. 1992, abstract/first page, DOI 10.1073/pnas.89.16.7521",
    },
    {
        "name": "atp_capacity_mM_s",
        "value": 0.70,
        "unit": "mM/s",
        "source_or_assumption": "upper-scale anchor: Layec et al. 2013, Table 3, 41.2 mM/min converted to 0.687 mM/s and rounded",
    },
    {
        "name": "substrate_factor",
        "value": 0.85,
        "unit": "dimensionless",
        "source_or_assumption": "assumption for nonoxygen substrate limitation",
    },
    {
        "name": "oxygen_cost_mol_o2_mol_atp",
        "value": 0.20,
        "unit": "mol O2 mol^-1 ATP",
        "source_or_assumption": (
            "SOURCE-BACKED 2026-10-02, was 'assumption consistent with an order-of-magnitude P/O "
            "near five'. Edge drawn to the older project's own mitochondrial_oxphos cell "
            "(~/projects/bodytwin/data/mitochondrial_oxphos/, 17 of 17 gates true), which "
            "reproduces Hinkle's classic and revised P/O and the Watt c8 mammalian ring. "
            "Converting its P/O grid (ATP per oxygen ATOM) by 1/(2*P/O): mammalian c8 with classic "
            "complex I gives 0.18333, with revised complex I 0.20370, so the mammalian band is "
            "OXYGEN_COST_BAND = [0.18333, 0.20370], 10.5 percent wide rather than an order of "
            "magnitude. The value kept here, 0.20, is Hinkle's pre-1999 classic P/O = 2.5 to five "
            "digits -- the old note's 'P/O near five' counted ATP per O2 while the standard "
            "convention counts per atom, which is the same number under two conventions and "
            "exactly the kind of ambiguity that produces a factor of two elsewhere."
        ),
    },
    {
        "name": "basal_oxygen_mol_m3_s",
        "value": 0.002,
        "unit": "mol m^-3 s^-1",
        "source_or_assumption": "assumption for oxygen-independent maintenance consumption",
    },
    {
        "name": "ck_rate_s",
        "value": 0.50,
        "unit": "s^-1",
        "source_or_assumption": "assumption for rapid creatine-kinase buffering",
    },
    {
        "name": "ck_capacity_mM_s",
        "value": 1.50,
        "unit": "mM/s",
        "source_or_assumption": "assumption bounding rapid PCr-to-ATP transfer",
    },
    {
        "name": "pcr_resynthesis_rate_s",
        "value": 0.05,
        "unit": "s^-1",
        "source_or_assumption": "assumption converting oxidative ATP surplus into PCr recovery",
    },
    {
        "name": "ion_leak_rate_s",
        "value": 0.25,
        "unit": "s^-1",
        "source_or_assumption": "assumption for activity-dependent membrane leak; no direct calibration",
    },
    {
        "name": "ion_basal_leak_fraction",
        "value": 0.0,
        "unit": "dimensionless",
        "source_or_assumption": "assumption that membrane leak is activity-dependent and vanishes in the zero-demand sham",
    },
    {
        "name": "ion_pump_rate_s",
        "value": 0.10,
        "unit": "s^-1",
        "source_or_assumption": "assumption for ATP- and oxygen-dependent Na/K pump restoration",
    },
    {
        "name": "functional weights",
        "value": [0.30, 0.50, 0.20],
        "unit": "dimensionless for ATP, PCr, ion terms",
        "source_or_assumption": "frozen preregistered composite outcome; not a measured scale",
    },
    {
        "name": "exercise_duration_s",
        "value": 250.0,
        "unit": "s",
        "source_or_assumption": "frozen protocol choice",
    },
    {
        "name": "recovery_duration_s",
        "value": 600.0,
        "unit": "s",
        "source_or_assumption": "frozen protocol choice",
    },
    {
        "name": "dt_s",
        "value": 0.25,
        "unit": "s",
        "source_or_assumption": "frozen numerical protocol; implicit oxygen update",
    },
    {
        "name": "radial_nodes",
        "value": 20,
        "unit": "count",
        "source_or_assumption": "frozen numerical protocol",
    },
]


REFERENCE_ANCHORS = [
    {
        "source": "Layec et al. 2013",
        "doi": "10.1152/japplphysiol.00257.2013",
        "location": "Table 3 and Table 4",
        "values": {
            "pcr_tau_free_flow_s": "33 +/- 21",
            "pcr_tau_reactive_hyperemia_s": "27 +/- 10",
            "atp_vmax_free_flow_mM_min": "28.7 +/- 13.3",
            "atp_vmax_reactive_hyperemia_mM_min": "41.2 +/- 13.6",
            "reoxygenation_mean_response_free_flow_s": "70 +/- 15",
            "reoxygenation_mean_response_reactive_hyperemia_s": "24 +/- 15",
        },
    },
    {
        "source": "Heskamp et al. 2021",
        "doi": "10.1113/JP280771",
        "location": "Table 2",
        "values": {
            "distal_k_pcr_min_inv": "0.44 +/- 0.26",
            "proximal_k_pcr_min_inv": "1.50 +/- 0.57",
            "distal_v_pcr_mM_min": "5.2 +/- 3.1",
            "proximal_v_pcr_mM_min": "23.3 +/- 8.9",
            "distal_k_o2hb_min_inv": "5.4 +/- 3.8",
            "proximal_k_o2hb_min_inv": "7.8 +/- 4.4",
        },
    },
    {
        "source": "Kushmerick et al. 1992",
        "doi": "10.1073/pnas.89.16.7521",
        "location": "abstract/full-text first page",
        "values": {
            "fast_twitch_atp_mM": "8",
            "fast_twitch_pcr_mM": "32",
            "fast_twitch_total_creatine_mM": "39",
        },
    },
    {
        "source": "Piiper and Scheid 1986",
        "doi": "10.1016/0034-5687(86)90118-0",
        "location": "abstract",
        "values": {
            "capillary_to_fiber_number": "about 2",
            "capillary_to_fiber_radius": "about 0.1",
            "geometry_conclusion": "Krogh model more adequate for skeletal muscle",
        },
    },
    {
        "source": "Segal and Faulkner 1985",
        "doi": "10.1152/ajpcell.1985.248.3.C265",
        "location": "abstract",
        "values": {
            "critical_radius_20C_mm": "1.19",
            "critical_radius_40C_mm": "0.51",
        },
    },
]


@dataclass(frozen=True)
class Geometry:
    r_m: np.ndarray
    volumes_m3: np.ndarray
    weights: np.ndarray
    a_left_s: np.ndarray
    a_right_s: np.ndarray
    boundary_s: float


def build_geometry(p: Parameters) -> Geometry:
    p.validate()
    n = p.radial_nodes
    r = np.linspace(0.0, p.radius_m, n)
    edges = np.empty(n + 1)
    edges[0] = 0.0
    edges[-1] = p.radius_m
    edges[1:-1] = 0.5 * (r[:-1] + r[1:])
    volumes = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    weights = volumes / np.sum(volumes)
    faces = np.empty(n + 1)
    faces[0] = 0.0
    faces[1:-1] = 2.0 * np.pi * edges[1:-1]
    faces[-1] = 2.0 * np.pi * p.radius_m
    left_distance = np.zeros(n)
    right_distance = np.empty(n)
    left_distance[1:] = np.diff(r)
    right_distance[:-1] = np.diff(r)
    right_distance[-1] = 0.5 * (r[-1] - r[-2])
    a_left = np.divide(
        faces[:-1] * p.oxygen_diffusivity_m2_s,
        volumes * left_distance,
        out=np.zeros(n),
        where=left_distance > 0.0,
    )
    a_right_internal = np.divide(
        faces[1:-1] * p.oxygen_diffusivity_m2_s,
        volumes[:-1] * right_distance[:-1],
        out=np.zeros(n - 1),
        where=right_distance[:-1] > 0.0,
    )
    a_boundary = p.oxygen_boundary_transfer_m_s * faces[-1] / volumes[-1]
    a_right = np.concatenate((a_right_internal, np.array([a_boundary])))
    return Geometry(r, volumes, weights, a_left, a_right, float(a_boundary))


def solve_tridiagonal(
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    d = np.asarray(diagonal, dtype=float).copy()
    b = np.asarray(rhs, dtype=float).copy()
    c = np.asarray(upper, dtype=float).copy()
    a = np.asarray(lower, dtype=float)
    n = len(d)
    if len(b) != n or len(c) != n - 1 or len(a) != n - 1:
        raise ValueError("invalid tridiagonal dimensions")
    for i in range(1, n):
        pivot = d[i - 1]
        if abs(pivot) < 1e-30:
            raise FloatingPointError("singular tridiagonal system")
        factor = a[i - 1] / pivot
        d[i] -= factor * c[i - 1]
        b[i] -= factor * b[i - 1]
    if abs(d[-1]) < 1e-30:
        raise FloatingPointError("singular tridiagonal system")
    x = np.zeros(n)
    x[-1] = b[-1] / d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (b[i] - c[i] * x[i + 1]) / d[i]
    return x


OXYGEN_SINK_SCHEME = "nonlinear_implicit"


def advance_oxygen(
    concentration: np.ndarray,
    p: Parameters,
    geometry: Geometry,
    dt_s: float,
    vmax_oxygen_mol_m3_s: float,
    boundary_concentration_mol_m3: float,
) -> np.ndarray:
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")
    if vmax_oxygen_mol_m3_s < 0.0 or boundary_concentration_mol_m3 < 0.0:
        raise ValueError("oxygen rates and boundary must be nonnegative")
    c = np.asarray(concentration, dtype=float)
    if np.any(c < 0.0):
        raise ValueError("oxygen inventory must be nonnegative")
    k = p.oxygen_michaelis_mol_m3
    z = np.zeros_like(c)  # Newton subsolution for the concave Michaelis sink.
    boundary_source = np.zeros_like(c)
    boundary_source[-1] = geometry.boundary_s * boundary_concentration_mol_m3
    lower = -dt_s * geometry.a_left_s[1:]
    upper = -dt_s * geometry.a_right_s[:-1]
    for _ in range(60):
        diffusion = np.zeros_like(z)
        diffusion[1:] += geometry.a_left_s[1:] * (z[:-1] - z[1:])
        diffusion[:-1] += geometry.a_right_s[:-1] * (z[1:] - z[:-1])
        diffusion[-1] -= geometry.boundary_s * z[-1]
        sink = vmax_oxygen_mol_m3_s * z / (k + z)
        residual = z - c - dt_s * (diffusion + boundary_source - sink)
        # Include uncancelled face operands in the arithmetic residual scale.
        diffusion_scale = (geometry.a_left_s + geometry.a_right_s) * np.abs(z)
        diffusion_scale[1:] += geometry.a_left_s[1:] * np.abs(z[:-1])
        diffusion_scale[:-1] += geometry.a_right_s[:-1] * np.abs(z[1:])
        scale = np.max(np.abs(z) + np.abs(c) + dt_s * (diffusion_scale + boundary_source + sink))
        boundary_port = p.oxygen_boundary_transfer_m_s * 2.0 * np.pi * p.radius_m * (boundary_concentration_mol_m3 - z[-1])
        inventory_terms = np.concatenate((z * geometry.volumes_m3, -c * geometry.volumes_m3,
                                          np.asarray([-dt_s * boundary_port]), dt_s * sink * geometry.volumes_m3))
        inventory_budget = 64 * np.finfo(float).eps * max(float(np.sum(np.abs(inventory_terms))), np.finfo(float).tiny)
        if (np.max(np.abs(residual)) <= 32 * np.finfo(float).eps * max(scale, np.finfo(float).tiny)
                and abs(float(np.sum(inventory_terms))) <= inventory_budget):
            return z
        slope = vmax_oxygen_mol_m3_s * k / np.square(k + z)
        diagonal = 1.0 + dt_s * (geometry.a_left_s + geometry.a_right_s + slope)
        correction = solve_tridiagonal(lower, diagonal, upper, -residual)
        # With F(0)<=0 and a concave sink, Newton stays below the root.
        z = z + correction
        if np.any(z < 0.0):
            raise RuntimeError("implicit oxygen Newton positivity invariant failed")
    raise RuntimeError("positive implicit oxygen solve did not converge")


def oxygen_factor(concentration: np.ndarray, geometry: Geometry, p: Parameters) -> float:
    mean = float(np.dot(geometry.weights, np.maximum(concentration, 0.0)))
    return mean / (p.oxygen_michaelis_mol_m3 + mean) if mean > 0.0 else 0.0


def energy_rates(
    p: Parameters,
    atp_mM: float,
    pcr_mM: float,
    o2_factor: float,
    demand_mM_s: float,
) -> dict[str, float]:
    atp_deficit = max(p.atp_rest_mM - atp_mM, 0.0)
    pcr_fraction = max(pcr_mM, 0.0) / p.pcr_rest_mM
    j_ck = min(p.ck_capacity_mM_s, p.ck_rate_s * atp_deficit * pcr_fraction)
    j_oxidative = p.atp_capacity_mM_s * p.substrate_factor * max(o2_factor, 0.0)
    gross_atp = j_oxidative + j_ck
    excess = max(gross_atp - demand_mM_s, 0.0)
    pcr_capacity = max(p.pcr_rest_mM - pcr_mM, 0.0)
    can_resynthesize = atp_mM >= p.atp_rest_mM - 1e-10
    j_pcr = min(excess, p.pcr_resynthesis_rate_s * pcr_capacity) if can_resynthesize else 0.0
    return {
        "j_oxidative_mM_s": float(j_oxidative),
        "j_ck_mM_s": float(j_ck),
        "j_pcr_mM_s": float(j_pcr),
        "d_atp_mM_s": float(gross_atp - demand_mM_s - j_pcr),
        "d_pcr_mM_s": float(j_pcr - j_ck),
    }


def ion_rates(
    p: Parameters,
    ion_disturbance: float,
    atp_mM: float,
    o2_factor: float,
    demand_mM_s: float,
) -> dict[str, float]:
    activity = min(1.0, demand_mM_s / max(p.work_demand_mM_s, 1e-12))
    leak = p.ion_leak_rate_s * (
        p.ion_basal_leak_fraction + (1.0 - p.ion_basal_leak_fraction) * activity
    )
    pump = (
        p.ion_pump_rate_s
        * max(ion_disturbance, 0.0)
        * min(atp_mM / p.atp_rest_mM, 1.0)
        * max(o2_factor, 0.0)
    )
    return {"ion_leak_s": float(leak), "ion_pump_s": float(pump), "d_ion_s": float(leak - pump)}


def demand_at_time(p: Parameters, t_s: float) -> float:
    return p.work_demand_mM_s if t_s < p.exercise_duration_s else p.recovery_demand_mM_s


def functional_score(p: Parameters, atp_mM: float, pcr_mM: float, ion_disturbance: float) -> float:
    score = (
        p.functional_pcr_weight * pcr_mM / p.pcr_rest_mM
        + p.functional_atp_weight * atp_mM / p.atp_rest_mM
        + p.functional_ion_weight * (1.0 - ion_disturbance)
    )
    return float(np.clip(score, 0.0, 1.0))


def unit_check(p: Parameters) -> dict[str, Any]:
    p.validate()
    boundary_concentration = p.oxygen_solubility_mol_m3_kpa * p.oxygen_boundary_kpa
    oxygen_rate = (
        p.oxygen_cost_mol_o2_mol_atp
        * p.atp_capacity_mM_s
        * p.substrate_factor
    )
    boundary_transfer_rate = 2.0 * p.oxygen_boundary_transfer_m_s / p.radius_m
    score_weights = np.array([
        p.functional_atp_weight,
        p.functional_pcr_weight,
        p.functional_ion_weight,
    ])
    checks = {
        "boundary_concentration": {
            "value": boundary_concentration,
            "unit": "mol m^-3",
            "pass": bool(np.isfinite(boundary_concentration) and boundary_concentration >= 0.0),
        },
        "oxygen_consumption_rate": {
            "value": oxygen_rate,
            "unit": "mol O2 m^-3 s^-1",
            "pass": bool(np.isfinite(oxygen_rate) and oxygen_rate >= 0.0),
        },
        "boundary_transfer_rate": {
            "value": boundary_transfer_rate,
            "unit": "s^-1",
            "pass": bool(np.isfinite(boundary_transfer_rate) and boundary_transfer_rate >= 0.0),
        },
        "functional_score": {
            "value": float(score_weights.sum()),
            "unit": "dimensionless",
            "pass": bool(np.isclose(score_weights.sum(), 1.0)),
        },
    }
    checks["all_pass"] = bool(all(check["pass"] for check in checks.values()))
    return checks


def threshold_time(
    times_s: np.ndarray,
    values: np.ndarray,
    threshold: float,
    origin_s: float,
) -> float | None:
    for i in range(1, len(values)):
        if values[i - 1] < threshold <= values[i]:
            span = values[i] - values[i - 1]
            fraction = 0.0 if span == 0.0 else (threshold - values[i - 1]) / span
            return float(times_s[i] - origin_s + fraction * (times_s[i] - times_s[i - 1]))
    if len(values) > 0 and values[0] >= threshold:
        return 0.0
    return None


def estimate_pcr_tau_s(
    times_s: np.ndarray,
    pcr_mM: np.ndarray,
    pcr_rest_mM: float,
    origin_s: float,
) -> float | None:
    mask = times_s >= origin_s
    t = times_s[mask] - origin_s
    y = pcr_mM[mask]
    end_value = float(y[0]) if len(y) else pcr_rest_mM
    amplitude = pcr_rest_mM - end_value
    if amplitude <= 1e-9:
        return None
    fraction = np.clip((y - end_value) / amplitude, 0.0, 1.0)
    deficit_fraction = np.clip(1.0 - fraction, 1e-8, 1.0)
    valid = (t >= 2.0) & (deficit_fraction > 0.0) & (deficit_fraction < 0.999)
    if np.count_nonzero(valid) < 3:
        return None
    slope, _ = np.polyfit(t[valid], np.log(deficit_fraction[valid]), 1)
    if slope >= -1e-9:
        return None
    return float(-1.0 / slope)


def summarize_result(
    p: Parameters,
    times_s: np.ndarray,
    atp_mM: np.ndarray,
    pcr_mM: np.ndarray,
    ion: np.ndarray,
    oxygen: np.ndarray,
    pcr_rates: np.ndarray,
    geometry: Geometry,
    null_model: bool,
) -> dict[str, Any]:
    exercise_origin = p.exercise_duration_s
    scores = np.array([
        functional_score(p, float(atp_mM[i]), float(pcr_mM[i]), float(ion[i]))
        for i in range(len(times_s))
    ])
    end_index = int(np.searchsorted(times_s, exercise_origin, side="left"))
    end_index = min(end_index, len(times_s) - 1)
    recovery_mask = times_s >= exercise_origin
    recovery_times = times_s[recovery_mask] - exercise_origin
    recovery_scores = scores[recovery_mask]
    recovery_pcr_fraction = (pcr_mM / p.pcr_rest_mM)[recovery_mask]
    recovery_ion_fraction = (1.0 - ion)[recovery_mask]
    t90 = threshold_time(recovery_times, recovery_scores, 0.90, 0.0)
    t90_pcr = threshold_time(recovery_times, recovery_pcr_fraction, 0.90, 0.0)
    t90_ion = threshold_time(recovery_times, recovery_ion_fraction, 0.90, 0.0)
    recovery_rates = pcr_rates[recovery_mask]
    positive_rates = recovery_rates[recovery_rates > 1e-9]
    initial_pcr_rate = float(positive_rates[0]) if len(positive_rates) else 0.0
    finite = bool(
        np.all(np.isfinite(times_s))
        and np.all(np.isfinite(atp_mM))
        and np.all(np.isfinite(pcr_mM))
        and np.all(np.isfinite(ion))
        and np.all(np.isfinite(oxygen))
    )
    valid_bounds = bool(
        np.all(atp_mM >= -1e-9)
        and np.all(atp_mM <= p.atp_rest_mM + 1e-9)
        and np.all(pcr_mM >= -1e-9)
        and np.all(pcr_mM <= p.pcr_rest_mM + 1e-9)
        and np.all(ion >= -1e-9)
        and np.all(ion <= 1.0 + 1e-9)
        and np.all(oxygen >= -1e-9)
    )
    return {
        "null_model": null_model,
        "t90_functional_s": t90,
        "t90_functional_censored": t90 is None,
        "t90_functional_lower_bound_s": float(t90 if t90 is not None else p.recovery_duration_s),
        "t90_pcr_s": t90_pcr,
        "t90_pcr_censored": t90_pcr is None,
        "t90_pcr_lower_bound_s": float(t90_pcr if t90_pcr is not None else p.recovery_duration_s),
        "t90_ion_s": t90_ion,
        "t90_ion_censored": t90_ion is None,
        "t90_ion_lower_bound_s": float(t90_ion if t90_ion is not None else p.recovery_duration_s),
        "pcr_tau_s": estimate_pcr_tau_s(times_s, pcr_mM, p.pcr_rest_mM, exercise_origin),
        "initial_recovery_pcr_rate_mM_s": initial_pcr_rate,
        "initial_recovery_v_pcr_mM_min": initial_pcr_rate * 60.0,
        "score_end_exercise": float(scores[end_index]),
        "score_final": float(scores[-1]),
        "atp_end_exercise_mM": float(atp_mM[end_index]),
        "pcr_end_exercise_mM": float(pcr_mM[end_index]),
        "ion_end_exercise": float(ion[end_index]),
        "atp_final_mM": float(atp_mM[-1]),
        "pcr_final_mM": float(pcr_mM[-1]),
        "ion_final": float(ion[-1]),
        "oxygen_min_mol_m3": float(np.min(oxygen)),
        "oxygen_max_mol_m3": float(np.max(oxygen)),
        "oxygen_center_final_mol_m3": float(oxygen[-1, 0]) if oxygen.ndim == 2 else float(oxygen[-1]),
        "oxygen_boundary_mol_m3": float(
            p.oxygen_solubility_mol_m3_kpa * p.oxygen_boundary_kpa
        ),
        "oxygen_mean_final_mol_m3": float(np.dot(geometry.weights, oxygen[-1])),
        "finite": finite,
        "valid_bounds": valid_bounds,
        "valid": bool(finite and valid_bounds),
    }


def simulate(
    p: Parameters = Parameters(),
    include_trajectory: bool = True,
    null_model: bool = False,
) -> dict[str, Any]:
    p.validate()
    geometry = build_geometry(p)
    boundary = p.oxygen_solubility_mol_m3_kpa * p.oxygen_boundary_kpa
    concentration = np.full(p.radial_nodes, boundary, dtype=float)
    atp = p.atp_rest_mM
    pcr = p.pcr_rest_mM
    ion = 0.0
    total_duration = p.exercise_duration_s + p.recovery_duration_s
    n_steps = int(np.ceil(total_duration / p.dt_s))
    stride_steps = max(1, int(round(p.sample_stride_s / p.dt_s)))
    times: list[float] = [0.0]
    atps: list[float] = [atp]
    pcrs: list[float] = [pcr]
    ions: list[float] = [ion]
    oxygen_rows: list[np.ndarray] = [concentration.copy()]
    pcr_rates: list[float] = [0.0]
    high_energy_projection_mM = 0.0
    for step in range(1, n_steps + 1):
        t_old = (step - 1) * p.dt_s
        demand = demand_at_time(p, t_old)
        if null_model:
            o2_factor = 1.0
        else:
            o2_factor = oxygen_factor(concentration, geometry, p)
        rates = energy_rates(p, atp, pcr, o2_factor, demand)
        vmax_oxygen = 0.0 if null_model else p.basal_oxygen_mol_m3_s + p.oxygen_cost_mol_o2_mol_atp * rates["j_oxidative_mM_s"]
        concentration = advance_oxygen(
            concentration,
            p,
            geometry,
            p.dt_s,
            vmax_oxygen,
            boundary,
        )
        o2_factor_after = 1.0 if null_model else oxygen_factor(concentration, geometry, p)
        rates_after = energy_rates(p, atp, pcr, o2_factor_after, demand)
        ion_state = ion_rates(p, ion, atp, o2_factor_after, demand)
        atp_raw = atp + p.dt_s * rates_after["d_atp_mM_s"]
        pcr_raw = pcr + p.dt_s * rates_after["d_pcr_mM_s"]
        atp = float(np.clip(atp_raw, 0.0, p.atp_rest_mM))
        pcr = float(np.clip(pcr_raw, 0.0, p.pcr_rest_mM))
        high_energy_projection_mM += (atp - atp_raw) + (pcr - pcr_raw)
        ion = float(np.clip(ion + p.dt_s * ion_state["d_ion_s"], 0.0, 1.0))
        if step % stride_steps == 0 or step == n_steps:
            times.append(float(step * p.dt_s))
            atps.append(atp)
            pcrs.append(pcr)
            ions.append(ion)
            oxygen_rows.append(concentration.copy())
            pcr_rates.append(rates_after["j_pcr_mM_s"])
    time_array = np.asarray(times, dtype=float)
    atp_array = np.asarray(atps, dtype=float)
    pcr_array = np.asarray(pcrs, dtype=float)
    ion_array = np.asarray(ions, dtype=float)
    oxygen_array = np.asarray(oxygen_rows, dtype=float)
    pcr_rate_array = np.asarray(pcr_rates, dtype=float)
    summary = summarize_result(
        p,
        time_array,
        atp_array,
        pcr_array,
        ion_array,
        oxygen_array,
        pcr_rate_array,
        geometry,
        null_model,
    )
    summary["high_energy_projection_mM"] = float(high_energy_projection_mM)
    summary["high_energy_projection_interpretation"] = "signed homeostatic projection port; not measured ATP heat"
    summary["exercise_duration_s"] = p.exercise_duration_s
    summary["recovery_duration_s"] = p.recovery_duration_s
    summary["oxygen_boundary_kpa"] = p.oxygen_boundary_kpa
    summary["work_demand_mM_s"] = p.work_demand_mM_s
    summary["recovery_demand_mM_s"] = p.recovery_demand_mM_s
    summary["radius_m"] = p.radius_m
    if include_trajectory:
        summary["trajectory"] = {
            "time_s": time_array.tolist(),
            "atp_mM": atp_array.tolist(),
            "pcr_mM": pcr_array.tolist(),
            "ion_disturbance": ion_array.tolist(),
            "functional_score": [
                functional_score(p, float(atp_array[i]), float(pcr_array[i]), float(ion_array[i]))
                for i in range(len(time_array))
            ],
            "pcr_resynthesis_rate_mM_s": pcr_rate_array.tolist(),
            "oxygen_center_mol_m3": oxygen_array[:, 0].tolist(),
            "oxygen_mean_mol_m3": (oxygen_array @ geometry.weights).tolist(),
        }
    return summary


def compact_result(result: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "null_model",
        "high_energy_projection_mM",
        "high_energy_projection_interpretation",
        "t90_functional_s",
        "t90_functional_censored",
        "t90_functional_lower_bound_s",
        "t90_pcr_s",
        "t90_pcr_censored",
        "t90_pcr_lower_bound_s",
        "t90_ion_s",
        "t90_ion_censored",
        "t90_ion_lower_bound_s",
        "pcr_tau_s",
        "initial_recovery_v_pcr_mM_min",
        "score_end_exercise",
        "score_final",
        "pcr_end_exercise_mM",
        "ion_end_exercise",
        "oxygen_center_final_mol_m3",
        "oxygen_mean_final_mol_m3",
        "valid",
    ]
    return {key: result[key] for key in keys}


def sensitivity(p: Parameters) -> list[dict[str, Any]]:
    names = ["oxygen_boundary_kpa", "work_demand_mM_s", "radius_m"]
    output: list[dict[str, Any]] = []
    for name in names:
        center = float(getattr(p, name))
        for factor in (0.5, 1.0, 1.5):
            candidate = replace(p, **{name: center * factor})
            result = simulate(candidate, include_trajectory=False)
            output.append({
                "parameter": name,
                "factor": factor,
                "value": center * factor,
                "result": compact_result(result),
            })
    return output


def grid(p: Parameters) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for oxygen in (2.0, 6.0, 13.3):
        for work in (0.30, 0.60, 0.90):
            for radius in (20e-6, 40e-6):
                candidate = replace(
                    p,
                    oxygen_boundary_kpa=oxygen,
                    work_demand_mM_s=work,
                    radius_m=radius,
                )
                result = simulate(candidate, include_trajectory=False)
                output.append({
                    "oxygen_boundary_kpa": oxygen,
                    "work_demand_mM_s": work,
                    "radius_m": radius,
                    "result": compact_result(result),
                })
    return output


def evaluate_criteria(
    baseline: dict[str, Any],
    control: dict[str, Any],
    corner: dict[str, Any],
    null_corner: dict[str, Any],
    sham: dict[str, Any],
    anoxia: dict[str, Any],
    grid_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    def row_for(oxygen: float, work: float, radius: float) -> dict[str, Any]:
        for row in grid_rows:
            if (
                np.isclose(row["oxygen_boundary_kpa"], oxygen)
                and np.isclose(row["work_demand_mM_s"], work)
                and np.isclose(row["radius_m"], radius)
            ):
                return row["result"]
        raise ValueError("requested grid row is missing")

    monotonic: dict[str, bool] = {}
    oxygen_rows = [row_for(value, 0.90, 40e-6) for value in (2.0, 6.0, 13.3)]
    oxygen_rates = [row["initial_recovery_v_pcr_mM_min"] for row in oxygen_rows]
    oxygen_times = [row["t90_functional_lower_bound_s"] for row in oxygen_rows]
    monotonic["oxygen_rate_nondecreasing"] = bool(
        all(oxygen_rates[i] <= oxygen_rates[i + 1] + 1e-9 for i in range(2))
    )
    monotonic["oxygen_t90_nonincreasing"] = bool(
        all(oxygen_times[i] >= oxygen_times[i + 1] - 1e-9 for i in range(2))
    )
    work_rows = [row_for(13.3, value, 40e-6) for value in (0.60, 0.90)]
    work_rates = [row["initial_recovery_v_pcr_mM_min"] for row in work_rows]
    work_times = [row["t90_functional_lower_bound_s"] for row in work_rows]
    monotonic["work_rate_nondecreasing"] = bool(
        all(work_rates[i] <= work_rates[i + 1] + 1e-9 for i in range(len(work_rates) - 1))
    )
    monotonic["work_t90_nondecreasing"] = bool(
        all(work_times[i] <= work_times[i + 1] + 1e-9 for i in range(len(work_times) - 1))
    )
    radius_rows = [row_for(2.0, 0.90, radius) for radius in (20e-6, 40e-6)]
    radius_rates = [row["initial_recovery_v_pcr_mM_min"] for row in radius_rows]
    radius_times = [row["t90_pcr_lower_bound_s"] for row in radius_rows]
    geometry_check = bool(
        radius_rates[1] <= radius_rates[0] + 1e-9
        and radius_times[1] >= radius_times[0] - 1e-9
    )
    corner_lower_bound = corner["t90_functional_lower_bound_s"]
    control_t90 = control["t90_functional_s"]
    control_rate = control["initial_recovery_v_pcr_mM_min"]
    joint_prediction = bool(
        control_t90 is not None
        and corner_lower_bound - control_t90 >= 90.0
        and corner["initial_recovery_v_pcr_mM_min"] <= 0.5 * control_rate
    )
    source_scale = bool(5.2 <= baseline["initial_recovery_v_pcr_mM_min"] <= 41.2)
    sham_ok = bool(
        sham["score_end_exercise"] >= 0.999
        and sham["pcr_end_exercise_mM"] >= 31.999
        and sham["ion_end_exercise"] <= 1e-6
    )
    anoxia_ok = bool(anoxia["initial_recovery_v_pcr_mM_min"] < 1e-9)
    null_stress_ok = bool(
        corner_lower_bound >= null_corner["t90_functional_lower_bound_s"]
        and corner["initial_recovery_v_pcr_mM_min"] <= null_corner["initial_recovery_v_pcr_mM_min"]
    )
    return {
        "monotonic_checks": monotonic,
        "geometry_check": geometry_check,
        "joint_prediction": joint_prediction,
        "source_scale_descriptive": source_scale,
        "sham_control": sham_ok,
        "anoxia_control": anoxia_ok,
        "null_stress_check": null_stress_ok,
        "all_mechanistic_checks_pass": bool(
            all(monotonic.values())
            and geometry_check
            and joint_prediction
            and sham_ok
            and anoxia_ok
            and null_stress_ok
        ),
        "published_rate_range_mM_min": [5.2, 41.2],
        "baseline_rate_mM_min": baseline["initial_recovery_v_pcr_mM_min"],
        "control_rate_mM_min": control_rate,
        "control_t90_s": control_t90,
        "corner_rate_mM_min": corner["initial_recovery_v_pcr_mM_min"],
        "corner_t90_s": corner["t90_functional_s"],
        "corner_t90_lower_bound_s": corner_lower_bound,
        "corner_t90_censored": corner["t90_functional_censored"],
        "null_corner_t90_lower_bound_s": null_corner["t90_functional_lower_bound_s"],
    }


def run_all(p: Parameters = Parameters()) -> dict[str, Any]:
    baseline = simulate(p, include_trajectory=True)
    control_p = replace(p, radius_m=20e-6)
    control = simulate(control_p, include_trajectory=False)
    corner_p = replace(
        p,
        oxygen_boundary_kpa=2.0,
        work_demand_mM_s=0.90,
        radius_m=40e-6,
    )
    corner = simulate(corner_p, include_trajectory=False)
    null_corner = simulate(corner_p, include_trajectory=False, null_model=True)
    sham_p = replace(p, work_demand_mM_s=0.0, recovery_demand_mM_s=0.0)
    sham = simulate(sham_p, include_trajectory=False)
    anoxia_p = replace(p, oxygen_boundary_kpa=0.0)
    anoxia = simulate(anoxia_p, include_trajectory=False)
    grid_rows = grid(p)
    return {
        "id": "BT-HX-Q088",
        "model": "radial oxygen diffusion with ATP-PCr and ion-balance state dynamics",
        "equations": EQUATIONS,
        "parameters": asdict(p),
        "parameter_table": PARAMETER_TABLE,
        "reference_anchors": REFERENCE_ANCHORS,
        "unit_check": unit_check(p),
        "baseline": baseline,
        "control": compact_result(control),
        "corner": compact_result(corner),
        "null_corner": compact_result(null_corner),
        "sham_control": compact_result(sham),
        "anoxia_control": compact_result(anoxia),
        "sensitivity": sensitivity(p),
        "grid": grid_rows,
        "criteria": evaluate_criteria(
            compact_result(baseline),
            compact_result(control),
            compact_result(corner),
            compact_result(null_corner),
            compact_result(sham),
            compact_result(anoxia),
            grid_rows,
        ),
        "provenance": {
            "measured_data": False,
            "fit_to_individuals": False,
            "external_literature_used_as_anchor": True,
            "predictions_are_not_measurements": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results.json")
    args = parser.parse_args()
    payload = run_all()
    output = Path(args.output)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(payload["criteria"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
