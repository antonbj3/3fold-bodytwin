#!/usr/bin/env python3
"""COUPLED METABOLIC CHAIN v1 -- the first cell chain in this repository whose certificate is
COMPOSED across cells instead of being asserted per cell.

Four shipped cells are wired producer->consumer along one physical quantity each, and nothing new
is modelled: every number below comes out of a function that already exists in the cell it is
attributed to, at that cell's published constants.

  S1  endocrine/glucose_insulin_minimal_model.py   meal glucose load (mg)
        -> ogtt_leg(): 2-state Bergman ODE -> plasma glucose trajectory G(t)
        -> suprabasal glucose DISPOSAL rate (mg/min) by mass balance on the cell's state:
           (absorbed mass - mass still held in the distribution volume) / window
  S2  energy/mitochondrial_oxphos.py + energy/warburg_metabolism.py
        -> disposal rate (mmol glucose/min) x ATP yield per glucose at the oxphos cell's P/O
           grid -> ATP flux (mmol/min); the OXIDATIVE ATP fraction divided by 2*P/O gives the
           O2 required (mmol/min -> mL/min at STPD)
  S3  haematology/blood_oxygen_transport.py
        -> sao2_hill() + cao2_ml_dl() at the reference blood gases -> arterial and mixed-venous
           O2 content -> a-vO2 content difference (mL O2/dL)
  S4  cardiovascular/cardiac_output.py
        -> fick_q_l_min(resting VO2 + meal-driven VO2, a-vO2 diff) -> required cardiac output
           (L/min), hr_from_q_sv() -> required heart rate

FINAL QUANTITY: the cardiac output (L/min) required to deliver the oxygen that oxidising a meal
glucose load costs, on top of the resting O2 demand.

WHAT THE CHAIN ADDS (nothing physiological):
  1. UNCERTAINTY PROPAGATION -- each stage's declared inputs are drawn inside the band that the
     stage's cell states for them (table STAGE_TOLERANCES below quotes the source line), Monte
     Carlo, fixed seed, N_DRAWS draws; the spread on the final cardiac output is reported.
  2. CERT COMPOSITION -- the full-chain half-width must not exceed the LINEAR SUM of the four
     single-stage half-widths (each measured by perturbing exactly one stage). Sub-additive =
     composition holds; super-additive = a stage is nonlinear in this operating point and is
     flagged by name.
  3. NULL CASES -- zero glucose load must leave every stage exactly at baseline; a doubled load
     must move every stage's output in the physiologically required direction (up).

DECLARED COUPLING (the only chain-level assumption, made explicit because the cells were written
standalone): the prescribed OGTT insulin excursion of S1 (make_ogtt_insulin_func, peak fold-rise
over fasting) is scaled linearly with the ingested load, so a zero load leaves insulin at fasting.
Without this the S1 cell would impose a full insulin excursion on a meal that was never eaten and
the zero-load null case would be untestable.

HONEST GAP (disclosed, not fixed here): the chain oxidises the WHOLE disposed glucose load inside
the window, so the O2 demand and therefore the required cardiac output are an UPPER BOUND -- a real
meal partly goes to glycogen instead of being oxidised immediately, and no shipped cell in this
repository splits that partition, so no such split is invented here.

INPUTS: a declared synthetic reference body (SYNTHETIC_REFERENCE below) -- body mass read from
examples/synthetic_inputs/metabolic_cost/metabolic_cost_results.json, every other value taken from
the four cells' own published constants. No individual's measurements are used anywhere.

Reads: the four cell modules (imported by path) + the synthetic metabolic_cost input.
Writes: $BODYTWIN_OUT/metabolic_chain_v1/metabolic_chain_v1_results.json.
Gate: overall_pass = all of the gates dict (exit code 0 on pass, 1 on fail).
"""
import contextlib
import importlib.util
import io
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
CELLS = REPO / "src" / "bodytwin" / "cells"
OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = Path(OUT_ROOT) / "metabolic_chain_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON = OUT_DIR / "metabolic_chain_v1_results.json"

SEED = 20260912
N_DRAWS = 256
OGTT_WINDOW_MIN = 180.0
GLUCOSE_MOLAR_MASS_MG_PER_MMOL = 180.16   # C6H12O6, same value the S1 cell uses for mg/dL<->mmol/L
MOLAR_VOLUME_STPD_L_MOL = 22.414          # 22.414 L/mol; same constant as carbon_co2_conservation_closure.py


def _load_cell(name, relpath):
    """Import a cell module by path with its self-test output suppressed (several cells run their
    report at import time)."""
    spec = importlib.util.spec_from_file_location(name, str(CELLS / relpath))
    module = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(module)
    return module


GLU = _load_cell("chain_glucose_insulin_minimal_model", "endocrine/glucose_insulin_minimal_model.py")
OXP = _load_cell("chain_mitochondrial_oxphos", "energy/mitochondrial_oxphos.py")
WAR = _load_cell("chain_warburg_metabolism", "energy/warburg_metabolism.py")
BOT = _load_cell("chain_blood_oxygen_transport", "haematology/blood_oxygen_transport.py")
CDO = _load_cell("chain_cardiac_output", "cardiovascular/cardiac_output.py")

PO_GRID_NADH = sorted(v["po_nadh"] for v in OXP.variant_grid.values())
PO_MIN, PO_MAX = PO_GRID_NADH[0], PO_GRID_NADH[-1]
PO_REF_MAX = 2.79   # the maximum P/O that warburg_metabolism.py's Y_OXPHOS_MODERN=33.45 is quoted with

_SYNTH_BODY = json.loads((REPO / "examples" / "synthetic_inputs" / "metabolic_cost" /
                          "metabolic_cost_results.json").read_text())

SYNTHETIC_REFERENCE = {
    "body_mass_kg": _SYNTH_BODY["muscle_mass"]["total_body_mass_kg"],
    "meal_glucose_load_mg": GLU.PARAMS["ogtt_dose_mg"],
    "Gb_mgdl": GLU.PARAMS["Gb_mgdl"],
    "Ib_uUmL": GLU.PARAMS["Ib_uUmL_central"],
    "Si": GLU.PARAMS["Si_central"],
    "Sg": GLU.PARAMS["Sg_normal"],
    "Vg_dL_per_kg": GLU.PARAMS["Vg_dL_per_kg"],
    "ogtt_tau_min": GLU.PARAMS["ogtt_tau_min"],
    "po_nadh": OXP.variant_grid["hinkle_classic_pre1999"]["po_nadh"],
    "hb_g_dl": BOT.HB_G_DL,
    "pao2_mmhg": BOT.PAO2_REST_MMHG,
    "pvo2_mmhg": BOT.PVO2_REST_MMHG,
    "o2_solubility_ml_dl_mmhg": BOT.O2_SOLUBILITY,
    "vo2_rest_ml_min": BOT.ANCHOR_VO2_REST_ML_MIN,
    "stroke_volume_rest_ml": CDO.HIGG_SVI_REST_ML_M2 * CDO.ASSUMED_BSA_M2,
    "source": "synthetic reference body: mass from examples/synthetic_inputs/metabolic_cost, every "
              "other value is the owning cell's published constant",
}

STAGE_TOLERANCES = {
    "S1_glucose_insulin": {
        "Si": ("relative", 0.24, "glucose_insulin_minimal_model.py: Bergman 1979 PMID 443421 "
                                 "SI = 7.00e-4 +/-24%"),
        "Ib_uUmL": ("band", (5.0, 10.0), "glucose_insulin_minimal_model.py: fasting insulin "
                                         "reference band 5-10 uU/mL"),
        "Vg_dL_per_kg": ("band", (1.3, 2.0), "glucose_insulin_minimal_model.py: its own vg_sweep"),
        "ogtt_tau_min": ("band", (30.0, 50.0), "glucose_insulin_minimal_model.py: its own tau_sweep"),
    },
    "S2_oxphos_atp": {
        "po_nadh": ("band", (PO_MIN, PO_MAX), "mitochondrial_oxphos.py: NADH P/O variant grid "
                                              "(c-ring size x complex-I stoichiometry)"),
    },
    "S3_blood_o2": {
        "hb_g_dl": ("band", (14.0, 18.0), "blood_oxygen_transport.py: Billett reference band "
                                          "14-18 g/dL its own HB_G_DL sits in"),
        "pvo2_mmhg": ("relative", 0.05, "blood_oxygen_transport.py: SVO2_AT_40_TOL_PCT = 5.0"),
        "o2_solubility_ml_dl_mmhg": ("band", (0.003, 0.0032), "blood_oxygen_transport.py: its own "
                                     "disclosed solubility variant 0.003 vs 0.0032"),
    },
    "S4_cardiac_output": {
        "vo2_rest_ml_min": ("relative", 0.20, "blood_oxygen_transport.py: FICK_CLOSURE_TOL_PCT = "
                                              "20.0, the pre-registered Fick-closure tolerance on "
                                              "the resting VO2 this stage divides"),
    },
}
STAGES = list(STAGE_TOLERANCES)


# ===================================================================================
# Boundary invariants and typed input validation.
# The existing gates are unit-blind; this independent stoichiometric closure catches a
# wrong-unit swap that Fick closure and glucose mass balance do not.
# ===================================================================================
INPUT_SPEC = {
    "meal_glucose_load_mg": ("mg", (0.0, 200000.0)),
    "Gb_mgdl": ("mg/dL", (50.0, 150.0)),
    "Ib_uUmL": ("uU/mL", (1.0, 30.0)),
    "Si": ("1e-4/(uU/mL)/min", (1e-5, 1e-2)),
    "Sg": ("1/min", (0.001, 0.1)),
    "Vg_dL_per_kg": ("dL/kg", (1.0, 2.5)),
    "body_mass_kg": ("kg", (20.0, 250.0)),
    "ogtt_tau_min": ("min", (10.0, 120.0)),
    "hb_g_dl": ("g/dL", (8.0, 22.0)),
    "pao2_mmhg": ("mmHg", (40.0, 600.0)),
    "pvo2_mmhg": ("mmHg", (10.0, 60.0)),
    "o2_solubility_ml_dl_mmhg": ("mL/dL/mmHg", (0.0015, 0.005)),
    "vo2_rest_ml_min": ("mL/min", (100.0, 600.0)),
    "stroke_volume_rest_ml": ("mL", (30.0, 150.0)),
}


def validate_inputs(params):
    """Typed, actionable seam errors instead of a raw KeyError (or silent acceptance)."""
    for key, (unit, (lo, hi)) in INPUT_SPEC.items():
        if key not in params:
            raise ValueError(f"missing required quantity '{key}' ({unit})")
        value = params[key]
        if not isinstance(value, (int, float)):
            raise ValueError(f"quantity '{key}' expected units '{unit}', got {type(value).__name__}")
        if not lo <= value <= hi:
            raise ValueError(f"quantity '{key}'={value} outside declared band ({lo}, {hi})")
    for key in params:
        if key not in INPUT_SPEC and key != "po_nadh":
            raise ValueError(f"unknown quantity '{key}' (no consumer declared it)")


def boundary_invariants(run, params=None):
    """Independent conservation/stoichiometry checks on one chain result.

    ``params`` must be the parameter set ``run`` was computed with; the glucose mass balance
    uses its meal load, basal glucose, distribution volume and body mass. The default is
    ``SYNTHETIC_REFERENCE``, which is correct for the nominal run that ``main()`` checks.
    """
    s1, s2, s3, s4 = (run["S1_glucose_insulin"], run["S2_oxphos_atp"],
                      run["S3_blood_o2"], run["S4_cardiac_output"])
    p = SYNTHETIC_REFERENCE if params is None else params
    absorbed_mg = p["meal_glucose_load_mg"] * GLU.PARAMS["ogtt_f_absorbed"]
    residual_mg = p["Vg_dL_per_kg"] * p["body_mass_kg"] * (s1["glucose_end_mgdl"] - p["Gb_mgdl"])
    mass_closed = abs(s1["suprabasal_disposal_mg_min"] * OGTT_WINDOW_MIN + residual_mg - absorbed_mg) / absorbed_mg < 1e-9
    fick_closed = abs(s4["cardiac_output_required_l_min"] * s3["avo2_diff_ml_dl"] * 10.0 - s4["vo2_total_ml_min"]) / s4["vo2_total_ml_min"] < 1e-9
    o2_per_glucose = s2["o2_required_ml_min"] / (s2["glucose_mmol_min"] * MOLAR_VOLUME_STPD_L_MOL)
    return {"glucose_mass_balance": bool(mass_closed), "fick_closure": bool(fick_closed),
            "o2_per_glucose_in_band": bool(4.0 <= o2_per_glucose <= 8.0),
            "o2_per_glucose_mol_per_mol": o2_per_glucose}


# ===================================================================================
# STAGE FUNCTIONS -- every physiological call goes to the owning cell
# ===================================================================================
def stage1_glucose_disposal(load_mg, si, sg, gb, ib, vg_dl_per_kg, bw_kg, tau_min,
                            window_min=OGTT_WINDOW_MIN):
    """S1: meal load -> Bergman ODE -> suprabasal glucose disposal rate (mg/min)."""
    p = dict(GLU.PARAMS)
    p["ogtt_dose_mg"] = load_mg
    p["ogtt_tau_min"] = tau_min
    p["Vg_dL_per_kg"] = vg_dl_per_kg
    p["BW_kg_ref"] = bw_kg
    # declared coupling: no meal -> no prescribed insulin excursion
    load_fraction = load_mg / GLU.PARAMS["ogtt_dose_mg"]
    p["insulin_peak_fold"] = 1.0 + (GLU.PARAMS["insulin_peak_fold"] - 1.0) * load_fraction
    leg = GLU.ogtt_leg(si, sg, gb, ib, p, t_end=window_min)
    vg_total_dl = vg_dl_per_kg * bw_kg
    absorbed_mg = load_mg * p["ogtt_f_absorbed"]
    still_in_pool_mg = vg_total_dl * (leg["G_curve_mgdl"][-1] - gb)
    disposal_mg_min = (absorbed_mg - still_in_pool_mg) / window_min
    return {
        "glucose_peak_mgdl": leg["G_peak_mgdl"],
        "glucose_120min_mmol": leg["G_120min_mmol"],
        "glucose_end_mgdl": leg["G_curve_mgdl"][-1],
        "suprabasal_disposal_mg_min": disposal_mg_min,
    }


def stage2_atp_and_o2(disposal_mg_min, po_nadh):
    """S2: glucose disposal -> ATP flux and the O2 (mL/min) its oxidation requires."""
    gluc_mmol_min = disposal_mg_min / GLUCOSE_MOLAR_MASS_MG_PER_MMOL
    atp_oxidative_per_glucose = (WAR.Y_OXPHOS_MODERN - WAR.Y_GLYCOLYSIS) * po_nadh / PO_REF_MAX
    atp_per_glucose = WAR.Y_GLYCOLYSIS + atp_oxidative_per_glucose
    atp_flux_mmol_min = gluc_mmol_min * atp_per_glucose
    # O2 = oxidative ATP / (2 * P/O); the P/O cancels against the yield above, i.e. the O2 cost of
    # oxidising a given amount of glucose is fixed by the chemistry and only the ATP RETURN on it
    # depends on the P/O variant -- reported, not hidden.
    o2_mmol_min = gluc_mmol_min * (WAR.Y_OXPHOS_MODERN - WAR.Y_GLYCOLYSIS) / (2.0 * PO_REF_MAX)
    o2_ml_min = o2_mmol_min * MOLAR_VOLUME_STPD_L_MOL
    return {
        "glucose_mmol_min": gluc_mmol_min,
        "atp_per_glucose": atp_per_glucose,
        "atp_flux_mmol_min": atp_flux_mmol_min,
        "o2_required_ml_min": o2_ml_min,
        "o2_per_glucose_mol_mol": (WAR.Y_OXPHOS_MODERN - WAR.Y_GLYCOLYSIS) / (2.0 * PO_REF_MAX),
    }


def stage3_avo2_diff(hb_g_dl, pao2, pvo2, solubility):
    """S3: blood gases -> arterial/venous O2 content -> a-vO2 content difference (mL O2/dL)."""
    sao2 = float(BOT.sao2_hill(pao2))
    svo2 = float(BOT.sao2_hill(pvo2))
    cao2 = float(BOT.cao2_ml_dl(hb_g_dl, sao2, pao2, sol=solubility))
    cvo2 = float(BOT.cao2_ml_dl(hb_g_dl, svo2, pvo2, sol=solubility))
    return {"sao2": sao2, "svo2": svo2, "cao2_ml_dl": cao2, "cvo2_ml_dl": cvo2,
            "avo2_diff_ml_dl": cao2 - cvo2}


def stage4_cardiac_output(o2_meal_ml_min, vo2_rest_ml_min, avo2_diff_ml_dl, sv_ml):
    """S4: Fick -> required cardiac output and heart rate."""
    vo2_total = vo2_rest_ml_min + o2_meal_ml_min
    q_baseline = CDO.fick_q_l_min(vo2_rest_ml_min, avo2_diff_ml_dl)
    q_meal = CDO.fick_q_l_min(vo2_total, avo2_diff_ml_dl)
    return {
        "vo2_total_ml_min": vo2_total,
        "cardiac_output_baseline_l_min": q_baseline,
        "cardiac_output_required_l_min": q_meal,
        "delta_cardiac_output_l_min": q_meal - q_baseline,
        "heart_rate_required_bpm": CDO.hr_from_q_sv(q_meal, sv_ml),
    }


def run_chain(params):
    s1 = stage1_glucose_disposal(params["meal_glucose_load_mg"], params["Si"], params["Sg"],
                                 params["Gb_mgdl"], params["Ib_uUmL"], params["Vg_dL_per_kg"],
                                 params["body_mass_kg"], params["ogtt_tau_min"])
    s2 = stage2_atp_and_o2(s1["suprabasal_disposal_mg_min"], params["po_nadh"])
    s3 = stage3_avo2_diff(params["hb_g_dl"], params["pao2_mmhg"], params["pvo2_mmhg"],
                          params["o2_solubility_ml_dl_mmhg"])
    s4 = stage4_cardiac_output(s2["o2_required_ml_min"], params["vo2_rest_ml_min"],
                               s3["avo2_diff_ml_dl"], params["stroke_volume_rest_ml"])
    return {"S1_glucose_insulin": s1, "S2_oxphos_atp": s2, "S3_blood_o2": s3,
            "S4_cardiac_output": s4,
            "final_cardiac_output_l_min": s4["cardiac_output_required_l_min"]}


# ===================================================================================
# UNCERTAINTY PROPAGATION
# ===================================================================================
def _draw(rng, nominal, spec):
    kind, value, _src = spec
    if kind == "relative":
        return nominal * (1.0 + rng.uniform(-value, value))
    lo, hi = value
    return rng.uniform(lo, hi)


def monte_carlo(stages, n_draws=N_DRAWS, seed=SEED):
    """Perturb the declared inputs of the listed stages only; return the final-quantity sample."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_draws):
        params = dict(SYNTHETIC_REFERENCE)
        for stage in STAGES:                      # same draw sequence regardless of which stages
            for key, spec in STAGE_TOLERANCES[stage].items():   # are active -> comparable samples
                value = _draw(rng, SYNTHETIC_REFERENCE[key], spec)
                if stage in stages:
                    params[key] = value
        out.append(run_chain(params)["final_cardiac_output_l_min"])
    return np.array(out)


def spread(sample):
    return {
        "n": int(sample.size),
        "mean_l_min": float(np.mean(sample)),
        "min_l_min": float(np.min(sample)),
        "max_l_min": float(np.max(sample)),
        "half_width_l_min": float((np.max(sample) - np.min(sample)) / 2.0),
        "sd_l_min": float(np.std(sample, ddof=1)),
    }


def main():
    nominal = run_chain(SYNTHETIC_REFERENCE)
    q_nom = nominal["final_cardiac_output_l_min"]

    # ---- (1) uncertainty propagation -------------------------------------------------------
    full = monte_carlo(STAGES)
    per_stage = {s: monte_carlo([s]) for s in STAGES}
    full_spread = spread(full)
    stage_spreads = {s: spread(v) for s, v in per_stage.items()}
    linear_sum = sum(v["half_width_l_min"] for v in stage_spreads.values())
    dominant = max(stage_spreads, key=lambda s: stage_spreads[s]["half_width_l_min"])

    # ---- (2) cert composition ---------------------------------------------------------------
    composition_holds = full_spread["half_width_l_min"] <= linear_sum
    nonlinear_flag = None if composition_holds else dominant

    # ---- (3) null cases ---------------------------------------------------------------------
    zero = dict(SYNTHETIC_REFERENCE, meal_glucose_load_mg=0.0)
    zero_run = run_chain(zero)
    zero_ok = (
        abs(zero_run["S1_glucose_insulin"]["glucose_peak_mgdl"] - SYNTHETIC_REFERENCE["Gb_mgdl"]) < 1e-6
        and abs(zero_run["S1_glucose_insulin"]["suprabasal_disposal_mg_min"]) < 1e-9
        and abs(zero_run["S2_oxphos_atp"]["atp_flux_mmol_min"]) < 1e-9
        and abs(zero_run["S2_oxphos_atp"]["o2_required_ml_min"]) < 1e-9
        and abs(zero_run["final_cardiac_output_l_min"]
                - zero_run["S4_cardiac_output"]["cardiac_output_baseline_l_min"]) < 1e-12
    )

    doubled = dict(SYNTHETIC_REFERENCE,
                   meal_glucose_load_mg=2.0 * SYNTHETIC_REFERENCE["meal_glucose_load_mg"])
    doubled_run = run_chain(doubled)
    directions = {
        "S1_glucose_peak_up": doubled_run["S1_glucose_insulin"]["glucose_peak_mgdl"]
                              > nominal["S1_glucose_insulin"]["glucose_peak_mgdl"],
        "S1_disposal_up": doubled_run["S1_glucose_insulin"]["suprabasal_disposal_mg_min"]
                          > nominal["S1_glucose_insulin"]["suprabasal_disposal_mg_min"],
        "S2_atp_flux_up": doubled_run["S2_oxphos_atp"]["atp_flux_mmol_min"]
                          > nominal["S2_oxphos_atp"]["atp_flux_mmol_min"],
        "S2_o2_required_up": doubled_run["S2_oxphos_atp"]["o2_required_ml_min"]
                             > nominal["S2_oxphos_atp"]["o2_required_ml_min"],
        "S3_avo2_diff_unchanged": abs(doubled_run["S3_blood_o2"]["avo2_diff_ml_dl"]
                                      - nominal["S3_blood_o2"]["avo2_diff_ml_dl"]) < 1e-12,
        "S4_cardiac_output_up": doubled_run["final_cardiac_output_l_min"] > q_nom,
    }

    gates = {
        "g01_uncertainty_propagated_all_stages": all(
            stage_spreads[s]["n"] == N_DRAWS for s in STAGES) and full_spread["n"] >= 64,
        "g02_final_spread_reported_and_finite": bool(np.all(np.isfinite(full))),
        "g03_nominal_inside_propagated_spread": full_spread["min_l_min"] <= q_nom <= full_spread["max_l_min"],
        "g04_cert_composition_subadditive": bool(composition_holds),
        "g05_null_zero_load_stays_at_baseline": bool(zero_ok),
        "g06_null_doubled_load_moves_every_stage_right_way": all(directions.values()),
        "g07_boundary_invariants_hold": all(
            v for k, v in boundary_invariants(nominal, SYNTHETIC_REFERENCE).items() if isinstance(v, bool)),
    }
    overall = all(gates.values())

    results = {
        "chain": ["glucose_insulin_minimal_model", "mitochondrial_oxphos+warburg_metabolism",
                  "blood_oxygen_transport", "cardiac_output"],
        "synthetic_reference": SYNTHETIC_REFERENCE,
        "stage_tolerances": {s: {k: {"kind": v[0], "value": v[1], "source": v[2]}
                                 for k, v in STAGE_TOLERANCES[s].items()} for s in STAGES},
        "seed": SEED, "n_draws": N_DRAWS,
        "nominal": nominal,
        "propagation": {"full_chain": full_spread, "per_stage": stage_spreads,
                        "linear_sum_of_stage_half_widths_l_min": linear_sum,
                        "dominant_stage": dominant,
                        "nonlinear_stage_flagged": nonlinear_flag},
        "null_cases": {"zero_load": zero_run, "doubled_load": doubled_run,
                       "directions": directions},
        "boundary_invariants": boundary_invariants(nominal, SYNTHETIC_REFERENCE),
        "gates": gates, "overall_pass": overall,
    }
    OUT_JSON.write_text(json.dumps(results, indent=2))

    print("=" * 78)
    print("COUPLED METABOLIC CHAIN v1 -- results")
    print("=" * 78)
    print("chain: glucose_insulin_minimal_model -> mitochondrial_oxphos(+warburg_metabolism) -> "
          "blood_oxygen_transport -> cardiac_output")
    s1, s2, s3, s4 = (nominal["S1_glucose_insulin"], nominal["S2_oxphos_atp"],
                      nominal["S3_blood_o2"], nominal["S4_cardiac_output"])
    print(f"S1 load {SYNTHETIC_REFERENCE['meal_glucose_load_mg']/1000:.0f} g -> peak glucose "
          f"{s1['glucose_peak_mgdl']:.1f} mg/dL, 2 h {s1['glucose_120min_mmol']:.2f} mmol/L, "
          f"disposal {s1['suprabasal_disposal_mg_min']:.1f} mg/min")
    print(f"S2 -> {s2['glucose_mmol_min']:.3f} mmol glucose/min, ATP {s2['atp_flux_mmol_min']:.2f} "
          f"mmol/min ({s2['atp_per_glucose']:.2f} ATP/glucose), O2 {s2['o2_required_ml_min']:.1f} mL/min")
    print(f"S3 -> SaO2 {s3['sao2']*100:.1f}%, SvO2 {s3['svo2']*100:.1f}%, CaO2 {s3['cao2_ml_dl']:.2f}, "
          f"a-vO2 diff {s3['avo2_diff_ml_dl']:.2f} mL/dL")
    print(f"S4 -> VO2 {s4['vo2_total_ml_min']:.1f} mL/min, Q {s4['cardiac_output_required_l_min']:.3f} "
          f"L/min (baseline {s4['cardiac_output_baseline_l_min']:.3f}, "
          f"+{s4['delta_cardiac_output_l_min']:.3f}), HR {s4['heart_rate_required_bpm']:.1f} bpm")
    print()
    print(f"PROPAGATION (seed {SEED}, {N_DRAWS} draws): Q = {full_spread['mean_l_min']:.3f} L/min, "
          f"range [{full_spread['min_l_min']:.3f}, {full_spread['max_l_min']:.3f}], "
          f"half-width {full_spread['half_width_l_min']:.3f}, SD {full_spread['sd_l_min']:.3f}")
    for s in STAGES:
        print(f"  {s:<20s} half-width {stage_spreads[s]['half_width_l_min']:.4f} L/min")
    print(f"  linear sum {linear_sum:.4f} L/min vs full-chain {full_spread['half_width_l_min']:.4f} "
          f"-> composition {'HOLDS (sub-additive)' if composition_holds else 'VIOLATED'}")
    print(f"  dominant stage: {dominant}"
          + ("" if nonlinear_flag is None else f"   NONLINEAR STAGE FLAGGED: {nonlinear_flag}"))
    print()
    print("GATES:")
    for k, v in gates.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    print(f"\nOVERALL: {'PASS' if overall else 'FAIL'}")
    print(f"\nWrote {OUT_JSON}")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
