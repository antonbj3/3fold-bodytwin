"""
Aqueous humor dynamics / Goldmann equation cert.

    IOP = (F - U) / C + Pev            [modern "expanded" Goldmann equation]
    F   = C * (IOP - Pev) + U          [rearranged: inflow = pressure-dependent + pressure-independent outflow]

Every numeric input below is a LIVE-VERIFIED (NCBI eutils) measured quantity from a named
primary source (the PMID for each number is given inline). This script does NOT
fit any parameter to the population IOP anchor — F, C, U, Pev are each sourced independently of the
tonometric population-IOP studies used to check the result (over-determination, not a tautology gate).

Reads: nothing (all inputs embedded). Writes: aqueous_humor_iop_evidence.json.
Gate: the eleven pre-registered gates in gates_summary; "overall" counts how many pass.
Pure Python/numpy, no scipy, deterministic.
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import json
import math
import numpy as np

_OUT_DIR = _os.path.join(OUT_ROOT, "aqueous_humor_iop")
OUT_EVIDENCE = _os.path.join(_OUT_DIR, "aqueous_humor_iop_evidence.json")

# ----------------------------------------------------------------------------------------------
# 1. LIVE-VERIFIED MEASURED INPUTS (independently sourced — see doc for PMID/DOI per number)
# ----------------------------------------------------------------------------------------------

# F: aqueous flow, healthy adults, fluorophotometry. Brubaker RF 1991 (Friedenwald Lecture), PMID 1748546.
F_MEAN, F_SD = 2.75, 0.63
F_95_LO, F_95_HI = 1.8, 4.3          # Brubaker's quoted 95% normal range

# F by age strata (fluorophotometry), Toris et al 1999, PMID 10218693 (n=51 young 20-30y, n=53 old 60y+)
F_YOUNG_MEAN, F_YOUNG_SD, F_YOUNG_N = 2.8, 0.8, 51
F_OLD_MEAN, F_OLD_SD, F_OLD_N = 2.4, 0.6, 53

# C: tonographic/trabecular outflow facility. Derived from Goel et al 2010 review (PMID 21293732):
# "resistance of the conventional aqueous drainage tissues is approximately 3-4 mmHg/(uL/min)" -> C=1/R.
R_LO, R_HI = 3.0, 4.0                 # mmHg per (uL/min)
C_LO, C_HI = 1.0 / R_HI, 1.0 / R_LO   # uL/min/mmHg  ->  0.25 .. 0.333
C_MID = 1.0 / 3.5                     # 0.2857

# U: uveoscleral (pressure-insensitive) outflow, Toris et al 1999, PMID 10218693.
U_YOUNG_MEAN, U_YOUNG_SD, U_YOUNG_N = 1.52, 0.81, 51
U_OLD_MEAN, U_OLD_SD, U_OLD_N = 1.10, 0.81, 53

# Anterior chamber volume, same cohorts (Toris 1999) -- used only for turnover-time cross-check.
ACV_YOUNG_MEAN, ACV_YOUNG_SD = 247.0, 39.0
ACV_OLD_MEAN, ACV_OLD_SD = 160.0, 39.0

# Pev: episcleral venous pressure, venomanometry. Phelps CD, Armaly MF 1978, PMID 619684.
PEV_MEAN, PEV_SD = 9.0, 1.6

# ----------------------------------------------------------------------------------------------
# 2. POPULATION IOP ANCHORS (measured independently of F/C/U/Pev above -- the falsifier target)
# ----------------------------------------------------------------------------------------------

# Beaver Dam Eye Study, applanation tonometry, Klein BE et al 1992, PMID 1607232.
IOP_BEAVERDAM_WOMEN, N_BEAVERDAM_WOMEN = 15.5, 2721
IOP_BEAVERDAM_MEN, N_BEAVERDAM_MEN = 15.3, 2135

# Secondary/textbook-tier population mean+-SD, via Goel 2010 review (PMID 21293732), citing
# Schottenstein/Ritch & Shields "The Glaucomas" 1989 textbook chapter (not independently live-fetched
# when this cell was written -- textbook-grade citation, disclosed as such).
IOP_POP_MEAN, IOP_POP_SD = 15.5, 2.6

# Distribution SHAPE anchor: approx-normal <21 mmHg, right-skew above, skew grows with age.
# Colton T, Ederer F 1980, PMID 7466591.
IOP_SHAPE_BREAKPOINT = 21.0

# POAG / ocular-hypertensive baseline (pre-treatment) diurnal IOP, real RCT, n=267.
# Alm A, Stjernschantz J 1995 (Scandinavian Latanoprost Study Group), PMID 9098273.
IOP_POAG_LO, IOP_POAG_HI = 24.6, 25.5

GLAUCOMA_THRESHOLD = 21.0   # task's pre-registered ">21 mmHg" glaucomatous-range threshold

# ----------------------------------------------------------------------------------------------
# 3. DRUG TRIAL DATA (real matched pre/post pairs, independent of the population-IOP anchor)
# ----------------------------------------------------------------------------------------------

# Carbonic anhydrase inhibitors, Maus TL et al 1997, PMID 9006424, n=40.
CAI = {
    "acetazolamide": {
        "F_pre": 3.18, "F_post": 2.23, "F_sd_pre": 0.70, "F_sd_post": 0.48,
        "IOP_pre": 12.5, "IOP_post": 10.1, "IOP_sd_pre": 2.2, "IOP_sd_post": 2.2,
    },
    "dorzolamide": {
        "F_pre": 3.18, "F_post": 2.65, "F_sd_pre": 0.70, "F_sd_post": 0.64,
        "IOP_pre": 12.5, "IOP_post": 10.8, "IOP_sd_pre": 2.2, "IOP_sd_post": 2.1,
    },
}

# Prostaglandin F2-alpha analog (PhXA41/latanoprost precursor), Toris CB et al 1993, PMID 8371915, n=22.
PG = {
    "IOP_drop_mean": 5.5, "IOP_drop_sem": 0.6,           # measured directly, day 8 vs baseline
    "U_baseline": 0.39, "U_baseline_sd": 0.20,
    "U_treated": 0.87, "U_treated_sd": 0.22,
    "U_contralateral_vehicle": 0.14, "U_contralateral_sd": 0.30,
    "F_and_C_significantly_changed": False,               # paper's explicit conclusion
}

# Independent full-RCT fractional confirmation, Alm & Stjernschantz 1995, PMID 9098273, n=267.
LATANOPROST_FRACTIONAL_DROP = {"morning": 0.31, "evening": 0.35}
TIMOLOL_FRACTIONAL_DROP = 0.27

# ----------------------------------------------------------------------------------------------
# 4. THE MODEL (pure algebra on the geometry: 1-node steady-state flow balance)
# ----------------------------------------------------------------------------------------------

def goldmann_iop(F, C, U, Pev):
    """IOP = (F-U)/C + Pev. The 1-node steady-state balance of a current source F, a conductance-C
    path to reference Pev, and a pressure-independent sink U."""
    return (F - U) / C + Pev


def goldmann_F(IOP, C, U, Pev):
    """Inverse direction: F = C*(IOP-Pev) + U."""
    return C * (IOP - Pev) + U


def implied_C(F, U, IOP, Pev):
    """Backward-solve C from independently measured F, U, IOP, Pev (over-determination direction)."""
    return (F - U) / (IOP - Pev)


def d_iop_d_C(F, C, U):
    """Analytic sensitivity: d/dC[(F-U)/C + Pev] = -(F-U)/C^2 (Pev drops out -- pure algebra)."""
    return -(F - U) / (C ** 2)


def d_iop_d_F(C):
    return 1.0 / C


def d_iop_d_U(C):
    return -1.0 / C


# ================================================================================================
# GATE 1 -- nested-model identity: modern 3-param equation at U=0 EXACTLY recovers the classical
# 1950 Goldmann 2-param form F=C(IOP-Pev), IOP=F/C+Pev. Machine-checked, not asserted.
# ================================================================================================
rng = np.random.default_rng(20260722)
_F = rng.uniform(1.0, 5.0, 5000)
_C = rng.uniform(0.05, 0.6, 5000)
_Pev = rng.uniform(5.0, 13.0, 5000)
modern_at_U0 = goldmann_iop(_F, _C, 0.0, _Pev)
classical = _F / _C + _Pev
gate1_nested_identity_exact = bool(np.allclose(modern_at_U0, classical, atol=1e-12))

# ================================================================================================
# GATE 2/3 -- FORWARD test: independently-sourced F, C, U, Pev sweep -> reproduce measured population IOP
# ================================================================================================
F_grid = np.array([F_MEAN - F_SD, F_MEAN, F_MEAN + F_SD, F_YOUNG_MEAN, F_OLD_MEAN])
C_grid = np.array([C_LO, C_MID, C_HI])
U_grid = np.array([U_OLD_MEAN, (U_OLD_MEAN + U_YOUNG_MEAN) / 2.0, U_YOUNG_MEAN])
Pev_grid = np.array([PEV_MEAN - PEV_SD, PEV_MEAN, PEV_MEAN + PEV_SD])

FF, CC, UU, PP = np.meshgrid(F_grid, C_grid, U_grid, Pev_grid, indexing="ij")
iop_sweep = goldmann_iop(FF, CC, UU, PP)
n_total = iop_sweep.size

# Falsifiable band: population mean +/- 2 SD (Goel/textbook-tier 15.5+/-2.6 -> [10.3, 20.7]), the
# pre-registered, EXTERNALLY-anchored (not tautological) target band.
band_lo, band_hi = IOP_POP_MEAN - 2 * IOP_POP_SD, IOP_POP_MEAN + 2 * IOP_POP_SD
frac_in_band = float(np.mean((iop_sweep >= band_lo) & (iop_sweep <= band_hi)))

# Non-degeneracy / void-floor disclosure: how wide IS the raw sweep, and how does the pass fraction
# change under a much STRICTER band (mean +/- 1 SD)? Reported, not gated -- honesty check on gate 3.
sweep_min, sweep_max, sweep_mean, sweep_std = (
    float(iop_sweep.min()), float(iop_sweep.max()), float(iop_sweep.mean()), float(iop_sweep.std())
)
band1_lo, band1_hi = IOP_POP_MEAN - 1 * IOP_POP_SD, IOP_POP_MEAN + 1 * IOP_POP_SD
frac_in_band_1sd = float(np.mean((iop_sweep >= band1_lo) & (iop_sweep <= band1_hi)))

# Central-estimate point (single best-guess plug-in, not cherry-picked from the sweep)
iop_central = goldmann_iop(F_MEAN, C_MID, (U_OLD_MEAN + U_YOUNG_MEAN) / 2.0, PEV_MEAN)

gate2_central_in_band = bool(band_lo <= iop_central <= band_hi)
gate3_majority_in_band = bool(frac_in_band >= 0.5)

# ================================================================================================
# GATE 4 -- forced adversary: P_ev forced to 0 must give a physiologically WRONG IOP (outside band)
# ================================================================================================
iop_pev0 = goldmann_iop(F_MEAN, C_MID, (U_OLD_MEAN + U_YOUNG_MEAN) / 2.0, 0.0)
gate4_pev0_outside_band = bool(not (band_lo <= iop_pev0 <= band_hi))
pev0_relative_error_pct = 100.0 * (iop_central - iop_pev0) / iop_central

# ================================================================================================
# GATE 5 -- BACKWARD: does a real measured POAG IOP + normal F (Larsson 1995: F not elevated) +
# measured Pev imply a REDUCED facility C, landing below the independently-sourced normal C_LO?
# Two age-blend variants shown (disclosed, not cherry-picked to the better-fitting one).
# ================================================================================================
poag_iop_mid = (IOP_POAG_LO + IOP_POAG_HI) / 2.0

variant_overall = dict(
    F=F_MEAN, U=(U_OLD_MEAN + U_YOUNG_MEAN) / 2.0, label="overall-population F,U average"
)
variant_age_matched = dict(F=F_OLD_MEAN, U=U_OLD_MEAN, label="age-matched-older F,U (Toris 60y+ arm)")

backward_variants = {}
for key, v in (("overall", variant_overall), ("age_matched_older", variant_age_matched)):
    C_poag = implied_C(v["F"], v["U"], poag_iop_mid, PEV_MEAN)
    pct_of_normal_lo = 100.0 * C_poag / C_LO
    backward_variants[key] = dict(
        label=v["label"], F_used=v["F"], U_used=v["U"],
        C_poag_implied=C_poag, pct_of_C_LO=pct_of_normal_lo,
        below_normal_range=bool(C_poag < C_LO),
    )

gate5_poag_facility_below_normal = bool(all(v["below_normal_range"] for v in backward_variants.values()))

# ================================================================================================
# GATE 6 -- forced adversary: PURE OVERPRODUCTION hypothesis (F up, C normal) explaining POAG IOP.
# Steelmanned: use NORMAL-central C (not cherry-picked low) and the LOWER (older-normal) U, which
# MINIMIZES the required F -- the most generous possible treatment of the adversary.
# ================================================================================================
F_needed_overproduction = goldmann_F(poag_iop_mid, C_MID, U_OLD_MEAN, PEV_MEAN)
overproduction_exceeds_95pct_ceiling = bool(F_needed_overproduction > F_95_HI)
overproduction_pct_above_ceiling = 100.0 * (F_needed_overproduction - F_95_HI) / F_95_HI
overproduction_fold_vs_measured_mean = F_needed_overproduction / F_MEAN

# Direct empirical kill: Larsson 1995 (PMID 7887840) measured F in real POAG (n=20) vs control (n=20)
# and found NO significant daytime difference -- directly contradicts the F_needed=5.7ish requirement.
larsson_1995_poag_f_elevated = False   # paper's explicit finding (quoted in doc)

gate6_overproduction_rejected = bool(overproduction_exceeds_95pct_ceiling and not larsson_1995_poag_f_elevated)

# ================================================================================================
# GATE 7 -- sensitivity/governor: |dIOP/dC| blows up as C falls (hyperbolic, ~1/C^2) -- evaluate at
# normal C_MID vs the backward-implied POAG C (age-matched-older variant), machine-computed ratio.
# ================================================================================================
C_poag_for_sensitivity = backward_variants["age_matched_older"]["C_poag_implied"]
sens_normal = d_iop_d_C(F_OLD_MEAN, C_MID, U_OLD_MEAN)
sens_poag = d_iop_d_C(F_OLD_MEAN, C_poag_for_sensitivity, U_OLD_MEAN)
sensitivity_ratio = abs(sens_poag) / abs(sens_normal)
predicted_ratio_from_1_over_C2 = (C_MID / C_poag_for_sensitivity) ** 2
gate7_governor_blowup = bool(sensitivity_ratio > 1.0)
gate7_matches_analytic_1_over_c2 = bool(math.isclose(sensitivity_ratio, predicted_ratio_from_1_over_C2, rel_tol=1e-9))

# ================================================================================================
# GATE 8 -- CAI (aqueous suppressant) dIOP/dF check: apply independently-sourced C to Maus 1997's
# OWN measured delta-F, compare predicted delta-IOP to Maus 1997's measured delta-IOP.
# ================================================================================================
cai_results = {}
for name, d in CAI.items():
    dF = d["F_pre"] - d["F_post"]
    dIOP_measured = d["IOP_pre"] - d["IOP_post"]
    dIOP_predicted = dF * d_iop_d_F(C_MID)   # = dF / C_MID
    ratio = dIOP_predicted / dIOP_measured
    cai_results[name] = dict(
        dF=dF, dIOP_measured=dIOP_measured, dIOP_predicted=dIOP_predicted,
        predicted_over_measured_ratio=ratio,
        same_order_of_magnitude=bool(0.4 <= ratio <= 2.5),
        frac_drop_measured=dIOP_measured / d["IOP_pre"],
    )
gate8_cai_same_order = bool(all(v["same_order_of_magnitude"] for v in cai_results.values()))

# ================================================================================================
# GATE 9 -- Prostaglandin dIOP/dU check: direction + qualitative mechanism (F,C unchanged, only U
# moves) is the PRE-REGISTERED falsifier; magnitude is reported honestly as an open tension.
# ================================================================================================
dU_vs_baseline = PG["U_treated"] - PG["U_baseline"]
dU_vs_contralateral = PG["U_treated"] - PG["U_contralateral_vehicle"]
dIOP_pred_vs_baseline = -dU_vs_baseline * (1.0 / C_MID)
dIOP_pred_vs_contralateral = -dU_vs_contralateral * (1.0 / C_MID)
dIOP_measured_pg = -PG["IOP_drop_mean"]

gate9_direction_correct = bool(dIOP_pred_vs_baseline < 0 and dIOP_pred_vs_contralateral < 0)
gate9_mechanism_isolated = bool(not PG["F_and_C_significantly_changed"])
pg_magnitude_ratio_baseline = dIOP_pred_vs_baseline / dIOP_measured_pg
pg_magnitude_ratio_contralateral = dIOP_pred_vs_contralateral / dIOP_measured_pg

# ================================================================================================
# GATE 10 -- aqueous turnover time: ACV / F, real matched Toris-1999 pairs, vs task's ~100min prior.
# ================================================================================================
turnover_young_min = ACV_YOUNG_MEAN / F_YOUNG_MEAN
turnover_old_min = ACV_OLD_MEAN / F_OLD_MEAN
turnover_heuristic_100pct = 250.0 / (0.01 * 250.0)   # Goel's "1.0%/min of ACV" heuristic -> 100 min exactly
turnover_heuristic_150pct = 250.0 / (0.015 * 250.0)  # Goel's "1.5%/min" upper heuristic -> 66.7 min

gate10_turnover_order_of_magnitude = bool(50.0 <= turnover_young_min <= 150.0 and 50.0 <= turnover_old_min <= 150.0)

# ================================================================================================
# GATE 11 -- void-floor necessity checks: is EACH term doing real, non-decorative work? (same
# discipline as the arterial_pressure notes sec 9's R/C void floors)
# ================================================================================================
U_mid_for_voidfloor = (U_OLD_MEAN + U_YOUNG_MEAN) / 2.0
iop_U_ignored = goldmann_iop(F_MEAN, C_MID, 0.0, PEV_MEAN)               # a model that forgets U entirely
iop_C_huge = goldmann_iop(F_MEAN, 1.0e6, U_mid_for_voidfloor, PEV_MEAN)  # C->inf: perfect drainage, no resistance
iop_C_tiny = goldmann_iop(F_MEAN, 1.0e-3, U_mid_for_voidfloor, PEV_MEAN) # C->0: angle-closure-like extreme

U_ignored_error_mmHg = iop_U_ignored - iop_central
C_necessary_work_mmHg = iop_central - iop_C_huge   # how much of central IOP is C "holding up" above Pev
gate11_U_necessary = bool(abs(U_ignored_error_mmHg) > 1.0)   # >1 mmHg is clinically non-trivial
gate11_C_necessary = bool(C_necessary_work_mmHg > 1.0)
gate11_C_tiny_blows_up = bool(iop_C_tiny > 1000)  # sanity: near-zero facility -> unbounded modeled IOP

# ================================================================================================
# ASSEMBLE EVIDENCE JSON
# ================================================================================================
evidence = {
    "model": "IOP = (F-U)/C + Pev (modern expanded Goldmann equation); F=C(IOP-Pev)+U",
    "inputs_independently_measured": {
        "F_uL_per_min": {"mean": F_MEAN, "sd": F_SD, "range95": [F_95_LO, F_95_HI],
                          "young": [F_YOUNG_MEAN, F_YOUNG_SD, F_YOUNG_N],
                          "old": [F_OLD_MEAN, F_OLD_SD, F_OLD_N]},
        "C_uL_per_min_per_mmHg": {"lo": C_LO, "mid": C_MID, "hi": C_HI},
        "U_uL_per_min": {"young": [U_YOUNG_MEAN, U_YOUNG_SD, U_YOUNG_N],
                         "old": [U_OLD_MEAN, U_OLD_SD, U_OLD_N]},
        "Pev_mmHg": {"mean": PEV_MEAN, "sd": PEV_SD},
    },
    "population_iop_anchor": {
        "beaverdam_women": [IOP_BEAVERDAM_WOMEN, N_BEAVERDAM_WOMEN],
        "beaverdam_men": [IOP_BEAVERDAM_MEN, N_BEAVERDAM_MEN],
        "textbook_mean_sd": [IOP_POP_MEAN, IOP_POP_SD],
        "band_used_mean_2sd": [band_lo, band_hi],
    },
    "gate1_nested_identity_exact_at_U0": gate1_nested_identity_exact,
    "gate2_central_estimate": {
        "iop_central_mmHg": iop_central, "in_band": gate2_central_in_band, "band": [band_lo, band_hi],
    },
    "gate3_sweep": {
        "n_total_combinations": int(n_total), "fraction_in_band": frac_in_band,
        "pass": gate3_majority_in_band,
        "sweep_raw_range": [sweep_min, sweep_max], "sweep_mean_sd": [sweep_mean, sweep_std],
        "fraction_in_tighter_1sd_band": frac_in_band_1sd, "band_1sd": [band1_lo, band1_hi],
    },
    "gate4_pev_forced_zero_adversary": {
        "iop_with_pev0": iop_pev0, "iop_central": iop_central,
        "relative_error_pct": pev0_relative_error_pct, "falls_outside_band": gate4_pev0_outside_band,
    },
    "gate5_poag_backward_solve": backward_variants | {"pass_both_variants_below_normal": gate5_poag_facility_below_normal},
    "gate6_overproduction_adversary": {
        "poag_iop_used": poag_iop_mid, "F_needed_uL_per_min": F_needed_overproduction,
        "F_95_hi_ceiling": F_95_HI, "pct_above_ceiling": overproduction_pct_above_ceiling,
        "fold_vs_measured_mean_F": overproduction_fold_vs_measured_mean,
        "larsson_1995_directly_measured_poag_F_elevated": larsson_1995_poag_f_elevated,
        "rejected": gate6_overproduction_rejected,
    },
    "gate7_sensitivity_governor": {
        "C_normal_mid": C_MID, "C_poag_implied": C_poag_for_sensitivity,
        "dIOP_dC_normal": sens_normal, "dIOP_dC_poag": sens_poag,
        "sensitivity_ratio": sensitivity_ratio, "matches_analytic_1_over_C2": gate7_matches_analytic_1_over_c2,
        "blowup_confirmed": gate7_governor_blowup,
    },
    "gate8_cai_dIOP_dF": cai_results,
    "gate8_pass": gate8_cai_same_order,
    "gate9_prostaglandin_dIOP_dU": {
        "dU_vs_baseline": dU_vs_baseline, "dU_vs_contralateral": dU_vs_contralateral,
        "dIOP_predicted_vs_baseline": dIOP_pred_vs_baseline,
        "dIOP_predicted_vs_contralateral": dIOP_pred_vs_contralateral,
        "dIOP_measured": dIOP_measured_pg,
        "ratio_vs_baseline": pg_magnitude_ratio_baseline, "ratio_vs_contralateral": pg_magnitude_ratio_contralateral,
        "direction_correct": gate9_direction_correct, "mechanism_isolated_to_U": gate9_mechanism_isolated,
        "magnitude_note": "predicted undershoots measured -- open tension, disclosed, not laundered",
    },
    "gate10_turnover_time": {
        "young_min": turnover_young_min, "old_min": turnover_old_min,
        "heuristic_1pct_min": turnover_heuristic_100pct, "heuristic_1p5pct_min": turnover_heuristic_150pct,
        "pass": gate10_turnover_order_of_magnitude,
        "note": "heuristic_* is algebraically 1/rate_pct (ACV cancels) -- NOT independent of the "
                "young/old direct measurement; both trace to the same Mayo/Nebraska fluorophotometry "
                "lineage, disclosed, not double-counted as two independent anchors.",
    },
    "gate11_void_floor_necessity": {
        "iop_central": iop_central,
        "iop_if_U_ignored": iop_U_ignored, "U_ignored_error_mmHg": U_ignored_error_mmHg,
        "iop_if_C_to_infinity": iop_C_huge, "C_necessary_work_mmHg": C_necessary_work_mmHg,
        "iop_if_C_near_zero": iop_C_tiny,
        "U_necessary": gate11_U_necessary, "C_necessary": gate11_C_necessary,
        "C_tiny_blows_up_unbounded": gate11_C_tiny_blows_up,
    },
    "gates_summary": {
        "gate1_nested_identity_exact": gate1_nested_identity_exact,
        "gate2_central_estimate_in_band": gate2_central_in_band,
        "gate3_sweep_majority_in_band": gate3_majority_in_band,
        "gate4_pev0_falls_outside_band": gate4_pev0_outside_band,
        "gate5_poag_facility_below_normal_both_variants": gate5_poag_facility_below_normal,
        "gate6_overproduction_rejected": gate6_overproduction_rejected,
        "gate7_governor_blowup_confirmed": gate7_governor_blowup,
        "gate8_cai_same_order_of_magnitude": gate8_cai_same_order,
        "gate9_pg_direction_and_mechanism_correct": bool(gate9_direction_correct and gate9_mechanism_isolated),
        "gate10_turnover_order_of_magnitude": gate10_turnover_order_of_magnitude,
        "gate11_U_and_C_both_necessary_nondecorative": bool(gate11_U_necessary and gate11_C_necessary),
    },
}

n_pass = sum(1 for v in evidence["gates_summary"].values() if v)
n_gate = len(evidence["gates_summary"])
evidence["overall"] = f"{n_pass}/{n_gate} PASS"

if __name__ == "__main__":
    _os.makedirs(_OUT_DIR, exist_ok=True)
    with open(OUT_EVIDENCE, "w") as fh:
        json.dump(evidence, fh, indent=2)
    print(json.dumps(evidence["gates_summary"], indent=2))
    print("OVERALL:", evidence["overall"])
    print()
    print("central IOP:", round(iop_central, 2), "mmHg  (band", [round(band_lo,2), round(band_hi,2)], ")")
    print("Pev=0 IOP:", round(iop_pev0, 2), "mmHg  (rel error", round(pev0_relative_error_pct,1), "%)")
    print("sweep fraction in band:", round(frac_in_band, 3), "  n=", n_total)
    print("sweep raw range:", [round(sweep_min,2), round(sweep_max,2)], " mean/sd:",
          round(sweep_mean,2), round(sweep_std,2))
    print("sweep fraction in TIGHTER 1sd band", [round(band1_lo,2), round(band1_hi,2)], ":",
          round(frac_in_band_1sd,3))
    print("backward C_poag (overall):", round(backward_variants["overall"]["C_poag_implied"], 4),
          " (age-matched):", round(backward_variants["age_matched_older"]["C_poag_implied"], 4),
          " vs C_LO=", C_LO)
    print("overproduction F_needed:", round(F_needed_overproduction, 3), " vs F_95_hi=", F_95_HI,
          " (+", round(overproduction_pct_above_ceiling,1), "% over ceiling,",
          round(overproduction_fold_vs_measured_mean,2), "x measured mean)")
    print("sensitivity ratio (governor blowup):", round(sensitivity_ratio, 2), "x")
    for k, v in cai_results.items():
        print(f"  CAI[{k}]: dF={v['dF']:.2f} -> predicted dIOP={v['dIOP_predicted']:.2f}"
              f" vs measured dIOP={v['dIOP_measured']:.2f}  ratio={v['predicted_over_measured_ratio']:.2f}")
    print("PG: predicted dIOP(vs baseline)=", round(dIOP_pred_vs_baseline,2), " (vs contralateral)=",
          round(dIOP_pred_vs_contralateral,2), " measured=", dIOP_measured_pg)
    print("turnover young/old (min):", round(turnover_young_min,1), round(turnover_old_min,1),
          " heuristic 1%/1.5%:", round(turnover_heuristic_100pct,1), round(turnover_heuristic_150pct,1))
    print("void floor: U ignored -> IOP=", round(iop_U_ignored,2), " (error", round(U_ignored_error_mmHg,2), "mmHg)")
    print("void floor: C->inf -> IOP=", round(iop_C_huge,2), " (C doing", round(C_necessary_work_mmHg,2), "mmHg of work)")
    print("void floor: C->~0 -> IOP=", round(iop_C_tiny,1), "mmHg (unbounded blowup, sanity)")
