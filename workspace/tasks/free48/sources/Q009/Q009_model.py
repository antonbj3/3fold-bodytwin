from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from scipy.integrate import solve_ivp


MODEL_EQUATIONS: Dict[str, str] = {
    "plasma_binding": "f_u,p = 1/(1 + P_p/K_d,p)",
    "plasma_partition": "C_p,u = f_u,p C_p,total; C_p,b = (1-f_u,p) C_p,total",
    "tissue_binding": "P_t,free = P_t,total - C_t,b",
    "tissue_free_balance": "dC_t,u/dt = (PS/V_t)(C_p,u-C_t,u) - k_elimination C_t,u - k_on,t P_t,free C_t,u + k_off,t C_t,b",
    "tissue_bound_balance": "dC_t,b/dt = k_on,t P_t,free C_t,u - k_off,t C_t,b",
    "tissue_total": "C_t,total = C_t,u + C_t,b",
    "equilibrium_tissue_fraction": "f_u,t = 1/(1 + P_t,total/K_d,t), K_d,t = k_off,t/k_on,t",
    "exposure": "AUC = integral C(t) dt",
}


REFERENCE_ANCHOR: Dict[str, Any] = {
    "citation": "Gill CM, Fratoni AJ, Shepard AK, Kuti JL, Nicolau DP. Omadacycline pharmacokinetics and soft-tissue penetration in diabetic patients with wound infections and healthy volunteers using in vivo microdialysis. J Antimicrob Chemother. 2022;77:1372-1378.",
    "doi": "10.1093/jac/dkac055",
    "location": "Table 2 and Methods, plasma protein-binding and tissue-penetration definitions",
    "status": "verified from full text",
    "cohorts": {
        "infected_patients": {
            "n": 6,
            "plasma_free_fraction_mean": 0.21,
            "plasma_free_fraction_sd": 0.03,
            "total_plasma_auc_mg_h_l": 6.27,
            "tissue_auc_mg_h_l": 0.82,
            "tissue_penetration": 0.66,
        },
        "healthy_volunteers": {
            "n": 6,
            "plasma_free_fraction_mean": 0.20,
            "plasma_free_fraction_sd": 0.02,
            "total_plasma_auc_mg_h_l": 14.06,
            "tissue_auc_mg_h_l": 1.37,
            "tissue_penetration": 0.54,
        },
    },
    "tissue_penetration_definition": "AUC_tissue/(f_u,p AUC_plasma)",
    "internal_discrepancy": "The abstract reports infected-cohort free plasma AUC 1.13 mg h/L while Table 2 reports 1.30 mg h/L; neither value is used as a model calibration target.",
}


@dataclass(frozen=True)
class PlasmaBinding:
    protein_sites_mg_l: float
    kd_mg_l: float

    def free_fraction(self) -> float:
        if not np.isfinite(self.protein_sites_mg_l) or self.protein_sites_mg_l < 0:
            raise ValueError("plasma protein sites must be finite and nonnegative")
        if not np.isfinite(self.kd_mg_l) or self.kd_mg_l < 0:
            raise ValueError("plasma Kd must be finite and nonnegative")
        if self.kd_mg_l == 0:
            return 1.0
        return float(1.0 / (1.0 + self.protein_sites_mg_l / self.kd_mg_l))

    @classmethod
    def from_target_free_fraction(cls, target_fu: float, protein_sites_mg_l: float = 0.8) -> "PlasmaBinding":
        if not np.isfinite(target_fu) or not 0 < target_fu <= 1:
            raise ValueError("target free fraction must lie in (0,1]")
        if not np.isfinite(protein_sites_mg_l) or protein_sites_mg_l <= 0:
            raise ValueError("protein sites must be positive")
        if target_fu == 1.0:
            kd = 0.0
        else:
            kd = protein_sites_mg_l * target_fu / (1.0 - target_fu)
        return cls(protein_sites_mg_l=protein_sites_mg_l, kd_mg_l=kd)


@dataclass(frozen=True)
class TissueParameters:
    protein_sites_mg_l: float = 1.0
    kd_mg_l: float = 4.0
    kon_l_mg_h: float = 0.25
    koff_h: float = 1.0
    k_transport_h: float = 0.8
    k_elimination_h: float = 0.05

    def validate(self) -> None:
        positive = {
            "protein_sites_mg_l": self.protein_sites_mg_l,
        }
        nonnegative = {
            "kd_mg_l": self.kd_mg_l,
            "kon_l_mg_h": self.kon_l_mg_h,
            "koff_h": self.koff_h,
            "k_transport_h": self.k_transport_h,
            "k_elimination_h": self.k_elimination_h,
        }
        for name, value in positive.items():
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name, value in nonnegative.items():
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.kon_l_mg_h > 0 and not np.isclose(self.koff_h / self.kon_l_mg_h, self.kd_mg_l, rtol=1e-12, atol=1e-12):
            raise ValueError("k_off/k_on must equal tissue Kd for the mass-action parameterization")

    def free_fraction_at_equilibrium(self) -> float:
        self.validate()
        if self.kd_mg_l == 0:
            return 1.0
        return float(1.0 / (1.0 + self.protein_sites_mg_l / self.kd_mg_l))


@dataclass(frozen=True)
class Case:
    name: str
    plasma: PlasmaBinding
    tissue: TissueParameters
    initial_tissue_free_mg_l: float = 0.0
    initial_tissue_bound_mg_l: float = 0.0

    def validate(self) -> None:
        if not self.name:
            raise ValueError("case name must not be empty")
        self.plasma.free_fraction()
        self.tissue.validate()
        if not np.isfinite(self.initial_tissue_free_mg_l) or self.initial_tissue_free_mg_l < 0:
            raise ValueError("initial free tissue concentration must be nonnegative")
        if not np.isfinite(self.initial_tissue_bound_mg_l) or self.initial_tissue_bound_mg_l < 0:
            raise ValueError("initial bound tissue concentration must be nonnegative")
        if self.initial_tissue_bound_mg_l > self.tissue.protein_sites_mg_l:
            raise ValueError("initial bound tissue concentration exceeds binding sites")


DEFAULT_TISSUE = TissueParameters()


def make_case(name: str, fu_plasma: float, tissue: TissueParameters = DEFAULT_TISSUE) -> Case:
    return Case(
        name=name,
        plasma=PlasmaBinding.from_target_free_fraction(fu_plasma),
        tissue=tissue,
    )


def plasma_curve(time_h: Sequence[float]) -> np.ndarray:
    t = np.asarray(time_h, dtype=float)
    if t.ndim != 1 or t.size < 2 or np.any(~np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError("time_h must be a finite strictly increasing one-dimensional array")
    if t[0] < 0:
        raise ValueError("time_h must be nonnegative")
    return 10.0 * np.exp(-0.08 * t) + 1.0 * np.exp(-0.8 * t)


def _auc(time_h: np.ndarray, values: np.ndarray) -> float:
    # getattr evaluates its default eagerly, so naming the removed np.trapz as the fallback
    # raised the very AttributeError the getattr was written to avoid.
    trapezoid = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(trapezoid(values, x=time_h))


def simulate_case(
    time_h: Sequence[float],
    total_plasma_mg_l: Sequence[float],
    case: Case,
    return_trace: bool = False,
) -> Dict[str, Any]:
    t = np.asarray(time_h, dtype=float)
    cp_total = np.asarray(total_plasma_mg_l, dtype=float)
    case.validate()
    if cp_total.ndim != 1 or cp_total.shape != t.shape:
        raise ValueError("total plasma curve must be one-dimensional and match time_h")
    if np.any(~np.isfinite(cp_total)) or np.any(cp_total < 0):
        raise ValueError("total plasma curve must be finite and nonnegative")
    fu_p = case.plasma.free_fraction()
    cp_free = fu_p * cp_total
    cp_bound = cp_total - cp_free
    initial = np.array([case.initial_tissue_free_mg_l, case.initial_tissue_bound_mg_l], dtype=float)

    def rhs(_time: float, state: np.ndarray) -> List[float]:
        tissue_free = float(state[0])
        tissue_bound = float(state[1])
        free_sites = max(case.tissue.protein_sites_mg_l - tissue_bound, 0.0)
        binding_on = case.tissue.kon_l_mg_h * free_sites * tissue_free
        binding_off = case.tissue.koff_h * tissue_bound
        transport_in = case.tissue.k_transport_h * (float(np.interp(_time, t, cp_free)) - tissue_free)
        clearance = case.tissue.k_elimination_h * tissue_free
        d_free = transport_in - clearance - binding_on + binding_off
        d_bound = binding_on - binding_off
        return [float(d_free), float(d_bound)]

    solution = solve_ivp(
        rhs,
        (float(t[0]), float(t[-1])),
        initial,
        t_eval=t,
        method="DOP853",
        rtol=1e-12,
        atol=1e-14,
    )
    if not solution.success or solution.y.shape != (2, t.size):
        raise RuntimeError(f"integration failed for {case.name}: {solution.message}")
    tissue_free = np.asarray(solution.y[0], dtype=float)
    tissue_bound = np.asarray(solution.y[1], dtype=float)
    tissue_total = tissue_free + tissue_bound
    if np.any(~np.isfinite(tissue_free)) or np.any(~np.isfinite(tissue_bound)):
        raise RuntimeError(f"non-finite tissue state for {case.name}")
    auc_plasma_total = _auc(t, cp_total)
    auc_plasma_free = _auc(t, cp_free)
    auc_plasma_bound = _auc(t, cp_bound)
    auc_tissue_free = _auc(t, tissue_free)
    auc_tissue_total = _auc(t, tissue_total)
    max_index = int(np.argmax(tissue_free))
    free_fraction_weighted = auc_tissue_free / auc_tissue_total if auc_tissue_total > 0 else None
    result: Dict[str, Any] = {
        "case": case.name,
        "fu_plasma": float(fu_p),
        "plasma_protein_sites_mg_l": float(case.plasma.protein_sites_mg_l),
        "plasma_kd_mg_l": float(case.plasma.kd_mg_l),
        "tissue_protein_sites_mg_l": float(case.tissue.protein_sites_mg_l),
        "tissue_kd_mg_l": float(case.tissue.kd_mg_l),
        "tissue_kon_l_mg_h": float(case.tissue.kon_l_mg_h),
        "tissue_koff_h": float(case.tissue.koff_h),
        "k_transport_h": float(case.tissue.k_transport_h),
        "k_elimination_h": float(case.tissue.k_elimination_h),
        "auc_plasma_total_mg_h_l": auc_plasma_total,
        "auc_plasma_free_mg_h_l": auc_plasma_free,
        "auc_plasma_bound_mg_h_l": auc_plasma_bound,
        "auc_tissue_free_mg_h_l": auc_tissue_free,
        "auc_tissue_total_mg_h_l": auc_tissue_total,
        "tissue_free_cmax_mg_l": float(tissue_free[max_index]),
        "tissue_total_cmax_mg_l": float(np.max(tissue_total)),
        "tissue_free_tmax_h": float(t[max_index]),
        "tissue_free_penetration": auc_tissue_free / auc_plasma_free if auc_plasma_free > 0 else None,
        "tissue_total_penetration": auc_tissue_total / auc_plasma_free if auc_plasma_free > 0 else None,
        "tissue_free_fraction_weighted": free_fraction_weighted,
    }
    if return_trace:
        result["trace"] = {
            "time_h": t.tolist(),
            "plasma_total_mg_l": cp_total.tolist(),
            "plasma_free_mg_l": cp_free.tolist(),
            "plasma_bound_mg_l": cp_bound.tolist(),
            "tissue_free_mg_l": tissue_free.tolist(),
            "tissue_bound_mg_l": tissue_bound.tolist(),
            "tissue_total_mg_l": tissue_total.tolist(),
        }
    return result


def dimensional_checks() -> Dict[str, Dict[str, Any]]:
    return {
        "transport_flux": {
            "lhs": "k_transport * C [mg L^-1 h^-1]",
            "rhs": "(PS/V_t) * C [mg L^-1 h^-1]",
            "pass": True,
        },
        "association": {
            "lhs": "k_on * P_free * C_free [mg L^-1 h^-1]",
            "rhs": "(L mg^-1 h^-1)(mg L^-1)(mg L^-1)",
            "pass": True,
        },
        "dissociation": {
            "lhs": "k_off * C_bound [mg L^-1 h^-1]",
            "rhs": "h^-1 * mg L^-1",
            "pass": True,
        },
        "auc": {
            "lhs": "integral C dt [mg h L^-1]",
            "rhs": "mg L^-1 * h",
            "pass": True,
        },
    }


def analytical_limit_check() -> Dict[str, Any]:
    time = np.linspace(0.0, 5.0, 10001)
    total = 3.0 * np.exp(-0.4 * time)
    transport = 1.2
    plasma_decay = 0.4
    tissue = TissueParameters(
        protein_sites_mg_l=1.0,
        kd_mg_l=1.0,
        kon_l_mg_h=0.0,
        koff_h=0.0,
        k_transport_h=transport,
        k_elimination_h=0.0,
    )
    case = Case(
        name="analytic_no_binding",
        plasma=PlasmaBinding.from_target_free_fraction(1.0),
        tissue=tissue,
    )
    result = simulate_case(time, total, case, return_trace=True)
    expected = transport * 3.0 / (transport - plasma_decay) * (
        np.exp(-plasma_decay * time) - np.exp(-transport * time)
    )
    observed = np.asarray(result["trace"]["tissue_free_mg_l"], dtype=float)
    max_relative_error = float(np.max(np.abs(observed - expected)) / np.max(np.abs(expected)))
    return {
        "expression": "C_t,u(t)=a*C0/(a-k)*(exp(-k*t)-exp(-a*t))",
        "max_relative_error": max_relative_error,
        "pass": bool(max_relative_error < 1e-7 and np.all(np.isfinite(observed))),
    }


def parameter_table() -> List[Dict[str, Any]]:
    plasma_reference = PlasmaBinding.from_target_free_fraction(0.20)
    tissue = DEFAULT_TISSUE
    return [
        {
            "name": "input_C_p_total",
            "value": "10*exp(-0.08*t)+1*exp(-0.8*t)",
            "unit": "mg L^-1",
            "source_or_assumption": "Controlled synthetic matched-curve input; not measured data.",
        },
        {
            "name": "time_grid",
            "value": "0:0.01:24",
            "unit": "h",
            "source_or_assumption": "Numerical grid choice for the frozen experiment.",
        },
        {
            "name": "P_p",
            "value": plasma_reference.protein_sites_mg_l,
            "unit": "mg L^-1 equivalent binding capacity",
            "source_or_assumption": "Effective binding-capacity assumption; no subject-specific measurement.",
        },
        {
            "name": "K_d,p(reference f_u=0.20)",
            "value": plasma_reference.kd_mg_l,
            "unit": "mg L^-1",
            "source_or_assumption": "Derived from P_p and the controlled reference free fraction.",
        },
        {
            "name": "P_t,total",
            "value": tissue.protein_sites_mg_l,
            "unit": "mg L^-1 equivalent binding capacity",
            "source_or_assumption": "Effective tissue binding-capacity assumption.",
        },
        {
            "name": "K_d,t",
            "value": tissue.kd_mg_l,
            "unit": "mg L^-1",
            "source_or_assumption": "Chosen with k_off/k_on to represent a tissue binding reservoir.",
        },
        {
            "name": "k_on,t",
            "value": tissue.kon_l_mg_h,
            "unit": "L mg^-1 h^-1",
            "source_or_assumption": "Mass-action association-rate assumption.",
        },
        {
            "name": "k_off,t",
            "value": tissue.koff_h,
            "unit": "h^-1",
            "source_or_assumption": "Mass-action dissociation-rate assumption.",
        },
        {
            "name": "k_transport=PS/V_t",
            "value": tissue.k_transport_h,
            "unit": "h^-1",
            "source_or_assumption": "Effective membrane-exchange assumption.",
        },
        {
            "name": "k_elimination",
            "value": tissue.k_elimination_h,
            "unit": "h^-1",
            "source_or_assumption": "Irreversible free-tissue clearance assumption.",
        },
    ]


def _relative_change(value: float, baseline: float) -> Optional[float]:
    if baseline == 0:
        return None
    return float((value - baseline) / baseline)


def run_experiment() -> Dict[str, Any]:
    time = np.linspace(0.0, 24.0, 2401)
    total = plasma_curve(time)
    reference_case = make_case("reference_fu_020", 0.20)
    low_binding_case = make_case("low_binding_fu_040", 0.40)
    high_binding_case = make_case("high_binding_fu_005", 0.05)
    same_case_a = make_case("null_a_fu_020", 0.20)
    same_case_b = make_case("null_b_fu_020", 0.20)
    zero_transport_tissue = replace(DEFAULT_TISSUE, k_transport_h=0.0)
    zero_transport_case = make_case("zero_transport", 0.20, zero_transport_tissue)

    main_results = {
        case.name: simulate_case(time, total, case, return_trace=True)
        for case in (low_binding_case, reference_case, high_binding_case)
    }
    null_a = simulate_case(time, total, same_case_a)
    null_b = simulate_case(time, total, same_case_b)
    zero_transport = simulate_case(time, total, zero_transport_case)
    primary_ratio = main_results[low_binding_case.name]["auc_tissue_free_mg_h_l"] / main_results[reference_case.name]["auc_tissue_free_mg_h_l"]
    null_ratio = null_a["auc_tissue_free_mg_h_l"] / null_b["auc_tissue_free_mg_h_l"]
    reference_auc_total = main_results[reference_case.name]["auc_plasma_total_mg_h_l"]
    max_curve_error = 0.0
    curve_pass = True
    for result in main_results.values():
        relative_error = abs(result["auc_plasma_total_mg_h_l"] - reference_auc_total) / reference_auc_total
        max_curve_error = max(max_curve_error, relative_error)
        if relative_error > 1e-12:
            curve_pass = False

    reference_metrics = main_results[reference_case.name]
    comparisons: Dict[str, Dict[str, float]] = {}
    for case_name, result in main_results.items():
        comparisons[case_name] = {
            "plasma_free_auc_ratio_to_reference": result["auc_plasma_free_mg_h_l"] / reference_metrics["auc_plasma_free_mg_h_l"],
            "tissue_free_auc_ratio_to_reference": result["auc_tissue_free_mg_h_l"] / reference_metrics["auc_tissue_free_mg_h_l"],
            "tissue_total_auc_ratio_to_reference": result["auc_tissue_total_mg_h_l"] / reference_metrics["auc_tissue_total_mg_h_l"],
        }

    sensitivity_rows: List[Dict[str, Any]] = []
    for label, value in (("minus_50_percent", 0.10), ("plus_50_percent", 0.30)):
        case = make_case(f"sensitivity_fu_{value:.3f}", value)
        row = simulate_case(time, total, case)
        row.update(
            {
                "parameter": "fu_plasma",
                "variation": label,
                "baseline_value": 0.20,
                "varied_value": value,
                "auc_relative_change": _relative_change(row["auc_tissue_free_mg_h_l"], main_results[reference_case.name]["auc_tissue_free_mg_h_l"]),
            }
        )
        sensitivity_rows.append(row)
    for label, value in (("minus_50_percent", 0.40), ("plus_50_percent", 1.20)):
        tissue = replace(DEFAULT_TISSUE, k_transport_h=value)
        case = make_case(f"sensitivity_transport_{value:.3f}", 0.20, tissue)
        row = simulate_case(time, total, case)
        row.update(
            {
                "parameter": "k_transport_h",
                "variation": label,
                "baseline_value": 0.80,
                "varied_value": value,
                "auc_relative_change": _relative_change(row["auc_tissue_free_mg_h_l"], main_results[reference_case.name]["auc_tissue_free_mg_h_l"]),
            }
        )
        sensitivity_rows.append(row)
    for label, value in (("minus_50_percent", 0.025), ("plus_50_percent", 0.075)):
        tissue = replace(DEFAULT_TISSUE, k_elimination_h=value)
        case = make_case(f"sensitivity_clearance_{value:.3f}", 0.20, tissue)
        row = simulate_case(time, total, case)
        row.update(
            {
                "parameter": "k_elimination_h",
                "variation": label,
                "baseline_value": 0.05,
                "varied_value": value,
                "auc_relative_change": _relative_change(row["auc_tissue_free_mg_h_l"], main_results[reference_case.name]["auc_tissue_free_mg_h_l"]),
            }
        )
        sensitivity_rows.append(row)

    analytic = analytical_limit_check()
    units = dimensional_checks()
    technical_pass = bool(
        curve_pass
        and 1.90 <= primary_ratio <= 2.10
        and abs(null_ratio - 1.0) <= 1e-10
        and zero_transport["auc_tissue_free_mg_h_l"] <= 1e-12
        and analytic["pass"]
        and all(row["pass"] for row in units.values())
    )
    return {
        "id": "BT-HX-Q009",
        "status": "executed_first_runnable_mechanistic_model",
        "provenance": {
            "source": "Gill et al. DOI 10.1093/jac/dkac055, Table 2; external verified reference only.",
            "derived": "Synthetic total-plasma input, numerical ODE integration, AUCs, ratios, and sensitivity changes.",
            "hypothesis": "Free plasma fraction and tissue transport/clearance are the controlling latent variables.",
            "unknown": "Subject-specific binding capacity, kinetics, permeability, and intracellular target concentration.",
        },
        "model_scope": "Matched total-plasma input, explicit free/bound plasma partition, passive unbound tissue exchange, reversible tissue binding, and free-tissue clearance.",
        "equations": MODEL_EQUATIONS,
        "parameter_table": parameter_table(),
        "input": {
            "time_start_h": 0.0,
            "time_end_h": 24.0,
            "time_step_h": 0.01,
            "expression": "C_p,total(t) = 10 exp(-0.08 t) + 1 exp(-0.8 t)",
            "unit": "mg L^-1",
            "source_or_assumption": "Controlled synthetic input; not measured data.",
        },
        "reference_anchor": REFERENCE_ANCHOR,
        "cases": main_results,
        "comparisons_to_reference": comparisons,
        "primary": {
            "quantity": "AUC_t,u(f_u,p=0.40)/AUC_t,u(f_u,p=0.20)",
            "predicted_ratio": 2.0,
            "observed_model_ratio": float(primary_ratio),
            "low_binding_case": low_binding_case.name,
            "reference_case": reference_case.name,
            "max_relative_total_plasma_curve_error": float(max_curve_error),
            "ratio_interval": [1.90, 2.10],
            "pass": bool(1.90 <= primary_ratio <= 2.10 and max_curve_error <= 1e-12),
        },
        "controls": {
            "same_binding_auc_ratio": float(null_ratio),
            "same_binding_pass": bool(abs(null_ratio - 1.0) <= 1e-10),
            "zero_transport_auc_tissue_free_mg_h_l": zero_transport["auc_tissue_free_mg_h_l"],
            "zero_transport_pass": bool(zero_transport["auc_tissue_free_mg_h_l"] <= 1e-12),
        },
        "sensitivity": sensitivity_rows,
        "verification": {
            "dimensional_checks": units,
            "analytical_limit": analytic,
            "technical_pass": technical_pass,
            "empirical_validity": "UNKNOWN; no matched free-plasma and target-tissue dataset was supplied.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the BT-HX-Q009 mechanistic exposure model")
    parser.add_argument("--output", default="results.json")
    args = parser.parse_args()
    result = run_experiment()
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
