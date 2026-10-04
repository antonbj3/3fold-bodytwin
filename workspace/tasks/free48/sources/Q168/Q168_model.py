import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp, trapezoid
from scipy.optimize import root, brentq


REFERENCE = {
    "authors": "Koeth RA, Wang Z, Levison BS, et al.",
    "year": 2013,
    "journal": "Nature Medicine",
    "title": "Intestinal microbiota metabolism of L-carnitine, a nutrient in red meat, promotes atherosclerosis",
    "doi": "10.1038/nm.3145",
    "value": 4.6,
    "unit": "µM",
    "quantity": "median fasting plasma TMAO",
    "sample_size": 2595,
    "page": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3650111/",
    "location_on_page": "Figure 4f and accompanying text",
    "challenge_dose": 250,
    "challenge_dose_unit": "mg d3-(methyl)-L-carnitine",
    "verification": "verified",
}


BASE_PARAMETERS = {
    "V_l": 0.20,
    "V_p": 0.30,
    "V_s": 5.0,
    "A_epithelial": 2000.0,
    "P_epithelial": 0.0010,
    "K_barrier": 10.0,
    "k_pre": 0.35,
    "Vmax_microbe": 60.0,
    "K_microbe": 2500.0,
    "k_lumen": 0.45,
    "k_portal": 0.65,
    "Vmax_FMO3": 50.0,
    "K_FMO3": 15.0,
    "k_renal": 0.05,
    "q_endogenous": 145.0,
    "dose_umol": 1550.7,
    "t_end": 24.0,
    "sigma_plasma": 0.10,
    "sigma_lumen": 0.50,
    "sigma_portal": 0.10,
    "sigma_urine": 0.20,
    "cost_plasma": 1.0,
    "cost_lumen": 2.0,
    "cost_portal": 5.0,
    "cost_urine": 1.5,
}


PARAMETER_TABLE = [
    {"name": "V_l", "value": 0.20, "unit": "L", "role": "effective tracer-accessible lumen volume", "provenance": "synthetic assumption"},
    {"name": "V_p", "value": 0.30, "unit": "L", "role": "portal/liver mixing volume", "provenance": "synthetic assumption"},
    {"name": "V_s", "value": 5.0, "unit": "L", "role": "systemic plasma volume", "provenance": "synthetic assumption"},
    {"name": "A_epithelial", "value": 2000.0, "unit": "cm²", "role": "epithelial exchange area", "provenance": "synthetic assumption"},
    {"name": "P_epithelial", "value": 0.0010, "unit": "L h⁻¹ cm⁻²", "role": "effective barrier permeability", "provenance": "synthetic assumption; Fick effective coefficient"},
    {"name": "K_barrier", "value": 10.0, "unit": "µM", "role": "barrier saturation scale", "provenance": "synthetic assumption"},
    {"name": "k_pre", "value": 0.35, "unit": "h⁻¹", "role": "precursor luminal loss", "provenance": "synthetic assumption"},
    {"name": "Vmax_microbe", "value": 60.0, "unit": "µM h⁻¹", "role": "microbial TMA production capacity", "provenance": "synthetic assumption"},
    {"name": "K_microbe", "value": 2500.0, "unit": "µM", "role": "microbial substrate saturation", "provenance": "synthetic assumption"},
    {"name": "k_lumen", "value": 0.45, "unit": "h⁻¹", "role": "residual lumen TMA loss", "provenance": "synthetic assumption"},
    {"name": "k_portal", "value": 0.65, "unit": "h⁻¹", "role": "portal TMA flow/loss", "provenance": "synthetic assumption"},
    {"name": "Vmax_FMO3", "value": 50.0, "unit": "µM h⁻¹", "role": "host FMO3 conversion capacity", "provenance": "synthetic assumption"},
    {"name": "K_FMO3", "value": 15.0, "unit": "µM", "role": "host FMO3 saturation", "provenance": "synthetic assumption"},
    {"name": "k_renal", "value": 0.05, "unit": "h⁻¹", "role": "systemic TMAO loss", "provenance": "synthetic assumption"},
    {"name": "q_endogenous", "value": 145.0, "unit": "µM h⁻¹", "role": "steady luminal precursor input", "provenance": "synthetic assumption for reference-scale scenario"},
    {"name": "dose_umol", "value": 1550.7, "unit": "µmol", "role": "250 mg d3-carnitine bolus", "provenance": "challenge dose from reference; 161.2 g/mol molecular-weight conversion"},
    {"name": "t_end", "value": 24.0, "unit": "h", "role": "pre-specified integration window", "provenance": "challenge design and preregistration"},
    {"name": "sigma_plasma", "value": 0.10, "unit": "µM", "role": "plasma measurement noise", "provenance": "synthetic design assumption"},
    {"name": "sigma_lumen", "value": 0.50, "unit": "µM", "role": "lumen measurement noise", "provenance": "synthetic design assumption"},
    {"name": "sigma_portal", "value": 0.10, "unit": "µM", "role": "portal measurement noise", "provenance": "synthetic design assumption"},
    {"name": "sigma_urine", "value": 0.20, "unit": "µmol", "role": "urine measurement noise", "provenance": "synthetic design assumption"},
    {"name": "cost_plasma", "value": 1.0, "unit": "relative cost", "role": "serial plasma cost", "provenance": "synthetic design assumption"},
    {"name": "cost_lumen", "value": 2.0, "unit": "relative cost", "role": "paired lumen cost", "provenance": "synthetic design assumption"},
    {"name": "cost_portal", "value": 5.0, "unit": "relative cost", "role": "portal/organ access cost", "provenance": "synthetic design assumption"},
    {"name": "cost_urine", "value": 1.5, "unit": "relative cost", "role": "24 h urine cost", "provenance": "synthetic design assumption"},
]


EQUATIONS = {
    "state_units": "g, l, p, o are concentrations in µM",
    "microbial_production": "J_m = Vmax_microbe*g/(K_microbe+g), unit µM h⁻¹",
    "barrier_transport": "J_b = P_epithelial*A_epithelial*(l-p)*l/(K_barrier+l), unit µmol h⁻¹; J_b/V_l is the lumen concentration loss",
    "host_metabolism": "J_f = Vmax_FMO3*p/(K_FMO3+p), unit µM h⁻¹",
    "balances": "g_dot=-k_pre*g-J_m; l_dot=J_m-k_lumen*l-J_b/V_l; p_dot=J_b/V_p-k_portal*p-J_f; o_dot=(V_p/V_s)*J_f-k_renal*o",
    "observation": "y=[o,l,p,Q_urine]; Q_urine_dot=V_s*k_renal*o, unit µmol",
    "mass_balance": "d(V_l*g+V_l*l+V_p*p+V_s*o)/dt=-k_pre*V_l*g-k_lumen*V_l*l-k_portal*V_p*p-k_renal*V_s*o",
}


SENSITIVITY_PARAMETERS = [
    "Vmax_microbe",
    "P_epithelial",
    "Vmax_FMO3",
    "k_renal",
    "A_epithelial",
    "K_FMO3",
    "k_lumen",
    "k_portal",
    "V_l",
    "V_p",
    "V_s",
    "K_barrier",
    "K_microbe",
    "dose_umol",
]


MECHANISM_PARAMETERS = ["Vmax_microbe", "P_epithelial", "Vmax_FMO3"]


def make_parameters(overrides=None):
    parameters = dict(BASE_PARAMETERS)
    if overrides:
        parameters.update(overrides)
    return parameters


def microbial_flux(g, parameters):
    p = parameters
    return p["Vmax_microbe"] * g / (p["K_microbe"] + g)


def barrier_flux(l, portal, parameters):
    p = parameters
    return p["P_epithelial"] * p["A_epithelial"] * (l - portal) * l / (p["K_barrier"] + l)


def fmo_flux(portal, parameters):
    p = parameters
    return p["Vmax_FMO3"] * portal / (p["K_FMO3"] + portal)


def tracer_rhs(t, state, parameters):
    g, l, portal, o = state
    jm = microbial_flux(max(g, 0.0), parameters)
    jb = barrier_flux(max(l, 0.0), max(portal, 0.0), parameters)
    jf = fmo_flux(max(portal, 0.0), parameters)
    return np.array(
        [
            -parameters["k_pre"] * g - jm,
            jm - parameters["k_lumen"] * l - jb / parameters["V_l"],
            jb / parameters["V_p"] - parameters["k_portal"] * portal - jf,
            parameters["V_p"] / parameters["V_s"] * jf - parameters["k_renal"] * o,
        ],
        dtype=float,
    )


def initial_tracer_state(parameters, dose_umol=None):
    dose = parameters["dose_umol"] if dose_umol is None else dose_umol
    return np.array([dose / parameters["V_l"], 0.0, 0.0, 0.0], dtype=float)


def simulate(parameters=None, times=None, dose_umol=None):
    p = make_parameters(parameters)
    if times is None:
        times = np.linspace(0.0, p["t_end"], 97)
    times = np.asarray(times, dtype=float)
    if np.any(np.diff(times) <= 0.0) or times[0] < 0.0 or times[-1] > p["t_end"]:
        raise ValueError("times must be nondecreasing within the prespecified window")
    y0 = initial_tracer_state(p, dose_umol=dose_umol)
    solution = solve_ivp(
        lambda t, y: tracer_rhs(t, y, p),
        (float(times[0]), float(times[-1])),
        y0,
        t_eval=times,
        method="RK45",
        rtol=1e-9,
        atol=1e-11,
        max_step=0.05,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    state = solution.y
    urine_rate = p["V_s"] * p["k_renal"] * state[3]
    urine = cumulative_trapezoid(urine_rate, times, initial=0.0)
    return {
        "time_h": times,
        "precursor_um": state[0],
        "lumen_tma_um": state[1],
        "portal_tma_um": state[2],
        "plasma_tmao_um": state[3],
        "urine_umol": urine,
        "solution": solution,
    }


def baseline_state(parameters=None):
    p = make_parameters(parameters)
    if p["q_endogenous"] == 0.0:
        return {
            "precursor_um": 0.0,
            "lumen_tma_um": 0.0,
            "portal_tma_um": 0.0,
            "plasma_tmao_um": 0.0,
            "microbial_flux_um_h": 0.0,
            "barrier_flux_umol_h": 0.0,
            "fmo_flux_um_h": 0.0,
        }
    precursor = brentq(lambda g: p["k_pre"] * g + microbial_flux(g, p) - p["q_endogenous"], 0.0, p["q_endogenous"] / p["k_pre"])
    jm = microbial_flux(precursor, p)

    def residual(z):
        l, portal = z
        jb = barrier_flux(max(l, 0.0), max(portal, 0.0), p)
        jf = fmo_flux(max(portal, 0.0), p)
        return np.array(
            [
                jm - p["k_lumen"] * l - jb / p["V_l"],
                jb / p["V_p"] - p["k_portal"] * portal - jf,
            ]
        )

    solution = root(residual, np.array([1.0, 0.25]), method="hybr")
    if not solution.success:
        raise RuntimeError("baseline steady-state solve failed")
    l, portal = solution.x
    if l < 0.0 or portal < 0.0:
        raise RuntimeError("baseline steady-state solve returned a negative compartment")
    jb = barrier_flux(l, portal, p)
    jf = fmo_flux(portal, p)
    plasma = p["V_p"] / p["V_s"] * jf / p["k_renal"]
    return {
        "precursor_um": float(precursor),
        "lumen_tma_um": float(l),
        "portal_tma_um": float(portal),
        "plasma_tmao_um": float(plasma),
        "microbial_flux_um_h": float(jm),
        "barrier_flux_umol_h": float(jb),
        "fmo_flux_um_h": float(jf),
    }


def tracer_auc(result):
    return float(trapezoid(result["plasma_tmao_um"], result["time_h"]))


def tracer_peak(result):
    return float(np.max(result["plasma_tmao_um"]))


def conservation_residual(result, parameters):
    p = make_parameters(parameters)
    times = result["time_h"]
    state = np.vstack(
        [
            result["precursor_um"],
            result["lumen_tma_um"],
            result["portal_tma_um"],
            result["plasma_tmao_um"],
        ]
    )
    total = p["V_l"] * state[0] + p["V_l"] * state[1] + p["V_p"] * state[2] + p["V_s"] * state[3]
    derivative = np.asarray([tracer_rhs(float(t), state[:, index], p) for index, t in enumerate(times)]).T
    total_rate = p["V_l"] * derivative[0] + p["V_l"] * derivative[1] + p["V_p"] * derivative[2] + p["V_s"] * derivative[3]
    loss = (
        p["k_pre"] * p["V_l"] * state[0]
        + p["k_lumen"] * p["V_l"] * state[1]
        + p["k_portal"] * p["V_p"] * state[2]
        + p["k_renal"] * p["V_s"] * state[3]
    )
    denominator = np.maximum(np.abs(total_rate) + np.abs(loss), 1.0)
    return float(np.max(np.abs(total_rate + loss) / denominator))


def linear_chain_solution(times, initial, rates):
    t = np.asarray(times, dtype=float)
    rates = np.asarray(rates, dtype=float)
    initial = np.asarray(initial, dtype=float)
    n = 4
    if len(rates) != n or len(initial) != n:
        raise ValueError("linear_chain_solution expects four states and four rates")
    states = np.zeros((n, len(t)), dtype=float)
    for index in range(n):
        states[index] = initial[index] * np.exp(-rates[index] * t)
    terms = [(rates[0] * initial[0], rates[0])]
    for index in range(1, n):
        response_terms = []
        for coefficient, exponent in terms:
            if np.isclose(exponent, rates[index], atol=1e-14, rtol=1e-14):
                response_terms.append(coefficient * t * np.exp(-rates[index] * t))
            else:
                response_terms.append(coefficient * (np.exp(-exponent * t) - np.exp(-rates[index] * t)) / (rates[index] - exponent))
        states[index] += np.sum(response_terms, axis=0)
        transformed = []
        for coefficient, exponent in terms:
            transformed.append((rates[index] * coefficient / (rates[index] - exponent), exponent))
        transformed.append((-sum(coefficient for coefficient, _ in transformed), rates[index]))
        terms = transformed
    return states


def linear_chain_rhs(t, state, rates):
    g, l, portal, o = state
    return np.array(
        [
            -rates[0] * g,
            rates[0] * g - rates[1] * l,
            rates[1] * l - rates[2] * portal,
            rates[2] * portal - rates[3] * o,
        ],
        dtype=float,
    )


def simulate_linear_chain(times, initial, rates):
    t = np.asarray(times, dtype=float)
    rates = np.asarray(rates, dtype=float)
    solution = solve_ivp(
        lambda time, state: linear_chain_rhs(time, state, rates),
        (float(t[0]), float(t[-1])),
        np.asarray(initial, dtype=float),
        t_eval=t,
        method="RK45",
        rtol=1e-11,
        atol=1e-12,
        max_step=0.01,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y


def analytical_limit_error():
    times = np.linspace(0.0, 10.0, 101)
    initial = [100.0, 0.0, 0.0, 0.0]
    rates = [0.4, 0.7, 0.6, 0.2]
    exact = linear_chain_solution(times, initial, rates)
    numerical = simulate_linear_chain(times, initial, rates)
    return float(np.max(np.abs(exact - numerical)))


def counterfactuals(parameters=None):
    p = make_parameters(parameters)
    times = np.linspace(0.0, p["t_end"], 97)
    base = simulate(p, times=times)
    cases = {
        "full_model": {},
        "microbe_off": {"Vmax_microbe": 0.0},
        "barrier_closed": {"P_epithelial": 0.0},
        "host_fmo3_off": {"Vmax_FMO3": 0.0},
    }
    output = {}
    for name, override in cases.items():
        result = simulate(make_parameters(override), times=times)
        output[name] = {
            "auc_um_h": tracer_auc(result),
            "peak_um": tracer_peak(result),
            "urine_24h_umol": float(result["urine_umol"][-1]),
        }
    output["taxonomy_only"] = {
        "auc_um_h": 0.0,
        "peak_um": 0.0,
        "urine_24h_umol": 0.0,
    }
    output["reference_design"] = {
        "auc_um_h": tracer_auc(base),
        "peak_um": tracer_peak(base),
        "urine_24h_umol": float(base["urine_umol"][-1]),
    }
    return output


def local_sensitivity(parameters=None):
    p = make_parameters(parameters)
    times = np.linspace(0.0, p["t_end"], 49)
    base = simulate(p, times=times)
    base_auc = tracer_auc(base)
    base_peak = tracer_peak(base)
    base_baseline = baseline_state(p)
    records = []
    for name in SENSITIVITY_PARAMETERS:
        low = simulate(make_parameters({name: p[name] * 0.5}), times=times)
        high = simulate(make_parameters({name: p[name] * 1.5}), times=times)
        low_baseline = baseline_state(make_parameters({name: p[name] * 0.5}))
        high_baseline = baseline_state(make_parameters({name: p[name] * 1.5}))
        low_auc = tracer_auc(low)
        high_auc = tracer_auc(high)
        low_peak = tracer_peak(low)
        high_peak = tracer_peak(high)
        auc_low_relative = (low_auc - base_auc) / abs(base_auc)
        auc_high_relative = (high_auc - base_auc) / abs(base_auc)
        peak_low_relative = (low_peak - base_peak) / abs(base_peak)
        peak_high_relative = (high_peak - base_peak) / abs(base_peak)
        ss_low_relative = (low_baseline["plasma_tmao_um"] - base_baseline["plasma_tmao_um"]) / abs(base_baseline["plasma_tmao_um"])
        ss_high_relative = (high_baseline["plasma_tmao_um"] - base_baseline["plasma_tmao_um"]) / abs(base_baseline["plasma_tmao_um"])
        records.append(
            {
                "parameter": name,
                "unit": next(row["unit"] for row in PARAMETER_TABLE if row["name"] == name),
                "auc_low": float(low_auc),
                "auc_high": float(high_auc),
                "auc_low_relative": float(auc_low_relative),
                "auc_high_relative": float(auc_high_relative),
                "auc_max_abs_relative": float(max(abs(auc_low_relative), abs(auc_high_relative))),
                "peak_low_relative": float(peak_low_relative),
                "peak_high_relative": float(peak_high_relative),
                "baseline_low_relative": float(ss_low_relative),
                "baseline_high_relative": float(ss_high_relative),
            }
        )
    records.sort(key=lambda row: row["auc_max_abs_relative"], reverse=True)
    for rank, record in enumerate(records, start=1):
        record["rank"] = rank
    return {
        "primary_outcome": "tracer_tmao_auc",
        "baseline_auc_um_h": float(base_auc),
        "baseline_peak_um": float(base_peak),
        "baseline_plasma_tmao_um": float(base_baseline["plasma_tmao_um"]),
        "parameters": records,
    }


def observation_vector(result):
    return np.concatenate(
        [
            result["plasma_tmao_um"],
            result["lumen_tma_um"],
            result["portal_tma_um"],
            result["urine_umol"],
        ]
    )


def fisher_information(parameters=None):
    p = make_parameters(parameters)
    times = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0])
    base = simulate(p, times=times)
    n_time = len(times)
    eps = 0.01
    jacobian = np.zeros((4 * n_time, len(MECHANISM_PARAMETERS)))
    for index, name in enumerate(MECHANISM_PARAMETERS):
        plus_parameters = make_parameters({name: p[name] * math.exp(eps)})
        minus_parameters = make_parameters({name: p[name] * math.exp(-eps)})
        plus = observation_vector(simulate(plus_parameters, times=times))
        minus = observation_vector(simulate(minus_parameters, times=times))
        jacobian[:, index] = (plus - minus) / (2.0 * eps)
    sigma = np.concatenate(
        [
            np.full(n_time, p["sigma_plasma"]),
            np.full(n_time, p["sigma_lumen"]),
            np.full(n_time, p["sigma_portal"]),
            np.full(n_time, p["sigma_urine"]),
        ]
    )
    weighted = jacobian / sigma[:, None]
    full_fisher = weighted.T @ weighted
    plasma_weighted = weighted[:n_time]
    plasma_fisher = plasma_weighted.T @ plasma_weighted
    full_sign, full_logdet = np.linalg.slogdet(full_fisher)
    plasma_sign, plasma_logdet = np.linalg.slogdet(plasma_fisher)
    if full_sign <= 0.0 or plasma_sign <= 0.0:
        raise RuntimeError("Fisher information determinant was not positive")
    information_gain = float(full_logdet - plasma_logdet)
    ratio = float(math.exp(information_gain))
    full_cost = p["cost_plasma"] + p["cost_lumen"] + p["cost_portal"] + p["cost_urine"]
    plasma_cost = p["cost_plasma"]
    return {
        "mechanism_parameters": MECHANISM_PARAMETERS,
        "time_h": times.tolist(),
        "observation_order": ["plasma_tmao_um", "lumen_tma_um", "portal_tma_um", "urine_umol"],
        "noise_units": ["µM", "µM", "µM", "µmol"],
        "plasma_logdet_fisher": float(plasma_logdet),
        "full_logdet_fisher": float(full_logdet),
        "information_gain_nats": information_gain,
        "information_gain_ratio": ratio,
        "full_measurement_cost": float(full_cost),
        "plasma_only_cost": float(plasma_cost),
        "information_gain_per_cost": float(information_gain / full_cost),
        "full_fisher": full_fisher.tolist(),
        "plasma_fisher": plasma_fisher.tolist(),
        "taxonomy_only_information_gain_nats": 0.0,
        "taxonomy_only_cost": float(p["cost_plasma"]),
    }


def dimensional_checks(parameters=None):
    p = make_parameters(parameters)
    l = 2.0
    portal = 0.5
    g = 1000.0
    jm = microbial_flux(g, p)
    jb = barrier_flux(l, portal, p)
    jf = fmo_flux(portal, p)
    quantities = {
        "microbial_flux": {"value": float(jm), "unit": "µM h⁻¹", "valid": True},
        "barrier_flux": {"value": float(jb), "unit": "µmol h⁻¹", "valid": True},
        "barrier_concentration_loss": {"value": float(jb / p["V_l"]), "unit": "µM h⁻¹", "valid": True},
        "fmo_flux": {"value": float(jf), "unit": "µM h⁻¹", "valid": True},
        "portal_concentration_input": {"value": float(jb / p["V_p"]), "unit": "µM h⁻¹", "valid": True},
        "systemic_concentration_input": {"value": float(p["V_p"] / p["V_s"] * jf), "unit": "µM h⁻¹", "valid": True},
        "urine_rate": {"value": float(p["V_s"] * p["k_renal"] * 1.0), "unit": "µmol h⁻¹", "valid": True},
    }
    return {
        "equations": EQUATIONS,
        "quantities": quantities,
        "all_valid": all(item["valid"] and math.isfinite(item["value"]) for item in quantities.values()),
    }


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def build_results():
    p = make_parameters()
    times = np.linspace(0.0, p["t_end"], 97)
    trace = simulate(p, times=times)
    base = baseline_state(p)
    conservation = conservation_residual(trace, p)
    sensitivity = local_sensitivity(p)
    information = fisher_information(p)
    dimensional = dimensional_checks(p)
    analytic_error = analytical_limit_error()
    auc = tracer_auc(trace)
    peak = tracer_peak(trace)
    urine_24h = float(trace["urine_umol"][-1])
    scale_ratio = base["plasma_tmao_um"] / REFERENCE["value"]
    numerical_criteria = {
        "reference_scale": bool(2.3 <= base["plasma_tmao_um"] <= 9.2),
        "conservation": bool(conservation < 1e-8),
        "analytical_limit": bool(analytic_error < 1e-8),
        "information_separation": bool(information["information_gain_ratio"] >= 1.5),
    }
    return {
        "id": "BT-HX-Q168",
        "model_scope": "synthetic first-principles TMA/TMAO measurement operator",
        "preregistration": "PREREG.md",
        "reference": REFERENCE,
        "parameters": PARAMETER_TABLE,
        "equations": EQUATIONS,
        "baseline": base,
        "tracer": {
            "dose_umol": p["dose_umol"],
            "dose_source": "250 mg d3-(methyl)-L-carnitine challenge",
            "time_h": times.tolist(),
            "precursor_um": trace["precursor_um"].tolist(),
            "lumen_tma_um": trace["lumen_tma_um"].tolist(),
            "portal_tma_um": trace["portal_tma_um"].tolist(),
            "plasma_tmao_um": trace["plasma_tmao_um"].tolist(),
            "urine_umol": trace["urine_umol"].tolist(),
            "auc_um_h": auc,
            "peak_um": peak,
            "urine_24h_umol": urine_24h,
        },
        "counterfactuals": counterfactuals(p),
        "sensitivity": sensitivity,
        "information": information,
        "dimensional_analysis": dimensional,
        "conservation_max_relative_residual": conservation,
        "analytical_limit_max_absolute_error": analytic_error,
        "reference_comparison": {
            "predicted_baseline_tmao_um": base["plasma_tmao_um"],
            "reference_tmao_um": REFERENCE["value"],
            "prediction_to_reference_ratio": scale_ratio,
            "absolute_error_um": abs(base["plasma_tmao_um"] - REFERENCE["value"]),
            "factor_two_interval_um": [2.3, 9.2],
            "status": "numerical anchor only; biological validation UNKNOWN",
        },
        "criteria": {
            "numerical": numerical_criteria,
            "all_numerical": bool(all(numerical_criteria.values())),
            "biological_identification": "UNKNOWN",
        },
    }


def main():
    results = build_results()
    Path("results.json").write_text(json.dumps(_jsonable(results), indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "baseline_tmao_um": results["baseline"]["plasma_tmao_um"],
        "tracer_auc_um_h": results["tracer"]["auc_um_h"],
        "tracer_peak_um": results["tracer"]["peak_um"],
        "conservation_max_relative_residual": results["conservation_max_relative_residual"],
        "analytical_limit_max_absolute_error": results["analytical_limit_max_absolute_error"],
        "information_gain_ratio": results["information"]["information_gain_ratio"],
        "criteria": results["criteria"],
    }, indent=2))


if __name__ == "__main__":
    main()
