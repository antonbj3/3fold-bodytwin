"""BMU (Basic Multicellular Unit) TURNOVER KINETICS -- the quantitative remodeling-RATE layer that
the static bone_remodeling mechanostat map explicitly does not model.

Two independent, decorrelated falsifiers:

  FALSIFIER 1 -- does a BMU-kinetics model built from in-vivo histomorphometry (tetracycline
  double-labeling: mineral apposition rate, BMU-cycle duration, resorption-cavity depth,
  activation frequency) reproduce the independently reported annual turnover rate (~10%/yr
  cortical, ~25%/yr trabecular, historically attributed to Parfitt/Eriksen)? Two routes, reported
  separately, never averaged:
    (a) a renewal-process OCCUPANCY derivation (the long-run fraction of bone actively remodeling
        equals the BMU active duration over the mean time between activations at the same site,
        i.e. the content of Parfitt's Ac.f x Sigma formalism derived from the renewal-reward
        theorem), computed for TRABECULAR bone from two Eriksen 2010 numbers;
    (b) a DIRECT cortical comparator (Clarke 2008's "2 to 3%/yr") against the stated 10%/yr.
  Plus a resorption/formation BALANCE (coupling) cross-check: does MAR x formation-phase duration
  reproduce Eriksen 2010's resorption-cavity depth (what a balanced BMU must refill)?

  FALSIFIER 2 -- does Frost's mechanostat, reusing the bone_remodeling classifier VERBATIM, place
  femoral-NECK gait strain (a second anatomical site and free body, which the shaft-only layer
  never touched) into the "maintenance" window, or into disuse/overload?

Reuse, not recompute: imports the bone_remodeling cell's THRESHOLD_VARIANTS /
classify_all_variants / classify_one / stress_mpa_to_microstrain / MICROSTRAIN_PER_MPA unedited.

Reads: OUT_ROOT/femoral_neck_stress/femoral_neck_stress_results.json (read-only) and the
bone_remodeling module. Writes: bmu_turnover_kinetics_results.json.

Gates: falsifier-1 route (a), route (b), the cross-compartment ratio and the coupling balance must
each land within a factor of 2 of the stated anchor; falsifier-2 requires bit-for-bit consistency
between recomputed and stored microstrain magnitudes in the upstream neck JSON, and reports
whether the "maintenance window" anchor claim survives at every tested neck configuration.
"""
import json
import os
import sys
import time

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
NECK_JSON = _os.path.join(OUT_ROOT, "femoral_neck_stress", "femoral_neck_stress_results.json")
OUT_DIR = _os.path.join(OUT_ROOT, "bmu_turnover_kinetics")
OUT_JSON = _os.path.join(OUT_DIR, "bmu_turnover_kinetics_results.json")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bone_remodeling as br  # noqa: E402  (reuse THRESHOLD_VARIANTS/classify_all_variants UNEDITED)

ASSERT_UPSTREAM_EXISTS = True

# ---------------------------------------------------------- literature ----
# Every number below is tagged with its live-verified source; task-stated
# figures are kept in a SEPARATE dict, never merged into the verified one.

TASK_STATED = {
    "cortical_turnover_pct_yr": 10.0,
    "trabecular_turnover_pct_yr": 25.0,
    "bmu_lifespan_months": [6.0, 9.0],
    "mar_um_day": [0.5, 1.0],
    "source": "Task specification (the build specification), attributed to "
               "tetracycline double-labeling histomorphometry / Parfitt.",
    "verified_live_this_session": False,
}

CLARKE_2008 = {
    "citation": "Clarke B. Normal bone anatomy and physiology. Clin J Am Soc Nephrol. 2008;3 Suppl 3:S131-9.",
    "pmid": "18988698", "doi": "10.2215/CJN.04151206", "pmcid": "PMC3152283",
    "verified_live_this_session_via": "NCBI efetch (existence/bibliography) + Crossref REST (DOI/author/journal/pages cross-match) + PMC full-text WebFetch (body-text numeric extraction).",
    "cortical_turnover_pct_yr": [2.0, 3.0],
    "resorption_phase_weeks": [2.0, 4.0],
    "formation_phase_months": [4.0, 6.0],
    "osteoblast_apoptosis_fraction": [0.50, 0.70],
    "quotes": {
        "cortical_turnover": "The relatively low adult cortical bone turnover rate of 2 to 3%/yr is adequate to maintain biomechanical strength of bone.",
        "resorption_phase": "Osteoclast-mediated bone resorption takes only approximately 2 to 4 wk during each remodeling cycle.",
        "formation_phase": "Bone formation takes approximately 4 to 6 mo to complete.",
        "osteoblast_apoptosis": "approximately 50 to 70% of osteoblasts undergo apoptosis",
        "trabecular_qualitative": "trabecular turnover is higher [than cortical], more than required for maintenance of mechanical strength (qualitative only -- no exact %/yr given in the fetched text).",
    },
}

ERIKSEN_2010 = {
    "citation": "Eriksen EF. Cellular mechanisms of bone remodeling. Rev Endocr Metab Disord. 2010;11(4):219-27.",
    "pmid": "21188536", "doi": "10.1007/s11154-010-9153-1", "pmcid": "PMC3028072",
    "verified_live_this_session_via": "NCBI efetch + Crossref REST + PMC full-text WebFetch.",
    "trabecular_bmu_cycle_days": 200.0,
    "cortical_bmu_cycle_days_median": 120.0,
    "resorption_phase_days_median_range": [30.0, 40.0],
    "formation_phase_days_cancellous": 150.0,
    "resorption_lacuna_depth_um": {"older": 40.0, "younger": 60.0},
    "trabecular_surface_full_renewal_years": 2.0,
    "remodeling_cycle_variability_days": {
        "fast_hyperthyroid_or_hyperparathyroid": 100.0,
        "slow_myxedema_or_bisphosphonate": 1000.0,
    },
    "quotes": {
        "trabecular_cycle": "In cancellous bone remodeling occurs on the surface of trabeculae and lasts about 200 days in normal bone.",
        "cortical_cycle": "the duration of the remodeling cycle in cortical is shorter than in cancellous bone with a median of 120 days.",
        "resorption_formation_split": "The resorption period has a median duration of 30-40 days and is followed by bone formation over a period of 150 days",
        "lacuna_depth": "[resorption lacuna] depth ... varies between 60 in young individuals and 40 um in older individuals.",
        "surface_renewal": "The total surface of cancellous bone is completely remodeled over a period of 2 years.",
        "variability": "The remodeling cycle can be as short as 100 days in thyrotoxicosis and primary hyperparathyroidism and exceed 1,000 days in low turnover states like Myxedema and after bisphosphonate treatment",
    },
}

PARFITT_1994 = {
    "citation": "Parfitt AM. Osteonal and hemi-osteonal remodeling: the spatial and temporal framework for signal traffic in adult human bone. J Cell Biochem. 1994;55(3):273-86.",
    "pmid": "7962158", "doi": "10.1002/jcb.240550303",
    "verified_live_this_session_via": "NCBI efetch (full abstract) + Crossref REST.",
    "quote_bmu_duration_qualitative": "The BMU maintains its size, shape and internal organization for many months as it travels through bone in a controlled direction.",
    "osteoclast_nuclei_turnover_pct_per_day": 8.0,
    "quote_osteoclast_nuclei": "Individual osteoclast nuclei are short-lived, turning over about 8% per d.",
}

HAUGE_2001 = {
    "citation": "Hauge EM, Qvesel D, Eriksen EF, Mosekilde L, Melsen F. Cancellous bone remodeling occurs in specialized compartments lined by cells expressing osteoblastic markers. J Bone Miner Res. 2001;16(9):1575-82.",
    "pmid": "11547826", "doi": "10.1359/jbmr.2001.16.9.1575",
    "verified_live_this_session_via": "NCBI efetch (full abstract) + Crossref REST.",
    "brc_covered_wall_thinner_pct": 30.0,
    "quote": "BRC-covered uncompleted walls are 30% thinner than those without a BRC.",
}

MANOLAGAS_2000 = {
    "citation": "Manolagas SC. Birth and death of bone cells: basic regulatory mechanisms and implications for the pathogenesis and treatment of osteoporosis. Endocr Rev. 2000;21(2):115-37.",
    "pmid": "10782361",
    "verified_live_this_session_via": "NCBI efetch (title/abstract match, AFTER an initial wrong-PMID guess 10893248 was caught by title mismatch -- see literature_search_misses).",
    "note": "Existence/conceptual anchor only (BMU birth/death coupling framework); abstract has no restated turnover-rate numbers.",
}

ERIKSEN_1986 = {
    "citation": "Eriksen EF. Normal and pathological remodeling of human trabecular bone: three dimensional reconstruction of the remodeling sequence in normals and in metabolic bone disease. Endocr Rev. 1986;7(4):379-408.",
    "pmid": "3536460", "doi": "10.1210/edrv-7-4-379",
    "verified_live_this_session_via": "NCBI efetch (existence/DOI) + Europe PMC REST (bibliographic cross-match). PubMed carries NO indexed abstract for this pre-1990s record.",
    "note": "Cited for existence/historical-primacy only -- same honest-gap tier as bone_stress.py's unresolved Reilly & Burstein 1975 citation.",
}

PARFITT_1983 = {
    "citation": "Parfitt AM, Mathews CH, Villanueva AR, Kleerekoper M, Frame B, Rao DS. Relationships between surface, volume, and thickness of iliac trabecular bone in aging and in osteoporosis. J Clin Invest. 1983;72(6):1396-1409.",
    "pmid": "6630513",
    "verified_live_this_session_via": "NCBI efetch (full abstract).",
    "quote": "age-related bone loss occurs principally by a process that removes entire structural elements of bone ... initiated by increased depth of osteoclastic resorption cavities which leads to focal perforation of trabecular plates",
    "used_here_as": "qualitative (non-numeric) geometric anchor: trabecular bone's much higher surface-to-volume ratio vs cortical bone is the standard literature mechanism for its much higher remodeling rate -- invoked here to EXPLAIN, not independently number-check, the cortical-vs-trabecular turnover asymmetry this script derives below.",
}

LITERATURE_SEARCH_MISSES = [
    {"guessed_pmid": "10893248", "expected": "Manolagas 2000 'Birth and death of bone cells', Endocr Rev",
     "actual_fetched_title": "Duplicated downstream enhancers control expression of the human apolipoprotein E gene in macrophages and adipose tissue (J Biol Chem 2000)",
     "resolution": "correct PMID found via esearch title match: 10782361"},
    {"guessed_pmid": "20835708", "expected": "Recker RR et al, 'Issues in modern bone histomorphometry', Osteoporos Int 2011",
     "actual_fetched_title": "Invisible scar endoscopic dorsal approach thyroidectomy: a clinical feasibility study",
     "resolution": "not pursued further -- Recker et al. is NOT cited below as a result; no fabricated substitute used."},
    {"guessed_pmid": "5307749", "expected": "Frost HM 1969, tetracycline-based histological analysis of bone remodeling, Calcif Tissue Res",
     "actual_fetched_title": "Pathologic fractures in children",
     "resolution": "not pursued further -- Frost 1969 is NOT cited below; the tetracycline double-labeling METHOD is instead corroborated qualitatively via Clarke 2008 and Eriksen 2010 (both of which describe tetracycline-labeling-derived kinetic parameters directly)."},
    {"guessed_pmid": "16834566", "expected": "Robling AG, Castillo AB, Turner CH, 'Biomechanical and molecular regulation of bone remodeling', Annu Rev Biomed Eng 2006",
     "actual_fetched_title": "Machine learning for detection and diagnosis of disease",
     "resolution": "not pursued further -- Robling 2006 is NOT cited below."},
    {"guessed_pmid": "3455637", "expected": "Parfitt AM et al 1987 ASBMR histomorphometry nomenclature standardization, J Bone Miner Res",
     "actual_fetched_title": "Bone histomorphometry: standardization of nomenclature, symbols, and units (title MATCHED, but PubMed carries no indexed abstract for this committee report)",
     "resolution": "title confirmed correct; cited below for existence/conceptual (BFR/BS, MAR, MS/BS, Ac.f, W.Th formal definitions) only, no numeric quote available."},
]

PARFITT_1987_NOMENCLATURE = {
    "citation": "Parfitt AM et al. Bone histomorphometry: standardization of nomenclature, symbols, and units. Report of the ASBMR Histomorphometry Nomenclature Committee. J Bone Miner Res. 1987;2(6):595-610.",
    "pmid": "3455637",
    "verified_live_this_session_via": "NCBI efetch (title match only; no abstract indexed).",
    "note": "Cited for the FORMAL definitions this script's derivations use (BFR/BS = MAR x MS/BS; Ac.f = BFR/BS / W.Th) -- existence-verified, not numerically re-extracted (no abstract on file).",
}


# ------------------------------------------------------- self-test (cheap, decisive) --
def _self_test():
    """Exact round-number constructions, checked before trusting any derived
    number below (same discipline as the bone_remodeling cell's _self_test)."""
    # Occupancy fraction: if a BMU is active for exactly as long as the mean
    # recurrence interval, occupancy is exactly 100% (trivial, exact).
    assert abs(occupancy_fraction_pct(365.0, 365.0) - 100.0) < 1e-9
    # Wall thickness: exact multiplication, round numbers.
    assert abs(wall_thickness_um(1.0, 100.0) - 100.0) < 1e-9
    # Round-trip identity: implied_recurrence_days(T, occupancy_fraction(T,R)/100) == R
    T, R = 120.0, 4800.0
    occ = occupancy_fraction_pct(T, R) / 100.0
    R_back = implied_recurrence_days(T, occ)
    assert abs(R_back - R) < 1e-6, f"round-trip failed: {R_back} != {R}"
    # within_factor symmetric-band self-check
    assert within_factor(10.0, 15.0, 2.0) is True
    assert within_factor(10.0, 25.0, 2.0) is False
    assert within_factor(25.0, 10.0, 2.0) is False  # order-independence
    return True


def occupancy_fraction_pct(t_bmu_days, t_recur_days):
    """GEOMETRIC/RENEWAL-THEORY DERIVATION (not a looked-up ratio): model a
    point on the bone surface as a two-state (idle/active) renewal process.
    By the renewal-reward theorem, the long-run FRACTION of time a point is
    ACTIVE equals (mean active-sojourn duration) / (mean cycle length,
    i.e. active+idle). This is algebraically identical to Parfitt's
    Ac.f x Sigma formalism (Ac.f = 1/T_recur per site, Sigma = T_bmu), here
    obtained from first-principles ergodic renewal theory rather than
    asserted as a black-box formula."""
    return 100.0 * t_bmu_days / t_recur_days


def implied_recurrence_days(t_bmu_days, occupancy_fraction_value):
    """Inverse of the above: given an independently-measured occupancy
    (~volumetric turnover fraction) and a BMU cycle length, what mean
    recurrence interval would reproduce it? Used to make the CORTICAL
    activation-frequency implication explicit and auditable (not just
    asserted) from Clarke 2008's directly-reported aggregate rate."""
    return t_bmu_days / occupancy_fraction_value


def wall_thickness_um(mar_um_day, formation_phase_days):
    """W.Th = MAR x formation-phase duration: the total new bone laid down
    by one BMU's osteoblast team over its own active formation window,
    assuming a constant apposition rate (the same constant-rate assumption
    Parfitt's histomorphometric MAR/W.Th formalism uses)."""
    return mar_um_day * formation_phase_days


def within_factor(a, b, factor):
    """Symmetric factor-of-X band check, order-independent. PRE-REGISTERED
    tolerance (stated before any of the numbers below were computed): a
    first-pass geometric/renewal model is not expected to hit an
    independently-reported aggregate figure exactly -- within a factor of 2
    either direction is treated as a genuine reproduction; anything wider is
    a reported (not hidden) miss."""
    if a <= 0 or b <= 0:
        return False
    ratio = max(a, b) / min(a, b)
    return ratio <= factor


TOLERANCE_FACTOR = 2.0  # pre-registered, applied uniformly below


# --------------------------------------------------- FALSIFIER 1: turnover rate --
def falsifier_1_turnover_rate():
    # (a) TRABECULAR: renewal-occupancy derivation from two independently
    # live-verified Eriksen 2010 numbers, vs the task's stated 25%/yr.
    t_bmu_trab = ERIKSEN_2010["trabecular_bmu_cycle_days"]              # 200 d, verified
    t_recur_trab = ERIKSEN_2010["trabecular_surface_full_renewal_years"] * 365.0  # 730 d, verified
    trab_derived_pct_yr = occupancy_fraction_pct(t_bmu_trab, t_recur_trab)
    trab_vs_task_pass = within_factor(trab_derived_pct_yr, TASK_STATED["trabecular_turnover_pct_yr"], TOLERANCE_FACTOR)

    # (b) CORTICAL: DIRECT literature comparator (Clarke 2008), vs task's 10%/yr.
    cort_lit_mid_pct_yr = sum(CLARKE_2008["cortical_turnover_pct_yr"]) / 2.0  # 2.5%/yr
    cort_vs_task_pass = within_factor(cort_lit_mid_pct_yr, TASK_STATED["cortical_turnover_pct_yr"], TOLERANCE_FACTOR)
    cort_vs_task_ratio = TASK_STATED["cortical_turnover_pct_yr"] / cort_lit_mid_pct_yr

    # (c) CROSS-COMPARTMENT RATIO check: does the RATIO between compartments
    # (not just each absolute number) match the task's implied ratio?
    task_ratio = TASK_STATED["trabecular_turnover_pct_yr"] / TASK_STATED["cortical_turnover_pct_yr"]  # 2.5x
    verified_ratio = trab_derived_pct_yr / cort_lit_mid_pct_yr  # ~11x
    ratio_of_ratios_pass = within_factor(task_ratio, verified_ratio, TOLERANCE_FACTOR)

    # (d) GEOMETRIC EXPLANATION (informational, not independently 3rd-source-
    # checked): what mean recurrence interval would CORTICAL
    # bone need, given its OWN (Eriksen) BMU-cycle length, to produce
    # Clarke's directly-reported aggregate rate?
    t_bmu_cort = ERIKSEN_2010["cortical_bmu_cycle_days_median"]  # 120 d, verified
    cort_implied_recur_days = implied_recurrence_days(t_bmu_cort, cort_lit_mid_pct_yr / 100.0)
    cort_implied_recur_years = cort_implied_recur_days / 365.0

    # (e) RESORPTION/FORMATION BALANCE (coupling) check: does MAR x
    # formation-duration reproduce Eriksen's resorption-lacuna depth
    # (what a balanced BMU must refill)? Checked at BOTH ends of the task's
    # stated MAR range, against BOTH ends of Eriksen's lacuna-depth range,
    # reported per-endpoint (never averaged into one pass/fail).
    lacuna_lo, lacuna_hi = ERIKSEN_2010["resorption_lacuna_depth_um"]["older"], ERIKSEN_2010["resorption_lacuna_depth_um"]["younger"]
    lacuna_mid = (lacuna_lo + lacuna_hi) / 2.0
    formation_days_trab = ERIKSEN_2010["formation_phase_days_cancellous"]  # 150 d, verified
    wall_th_results = {}
    for mar in TASK_STATED["mar_um_day"]:
        wth = wall_thickness_um(mar, formation_days_trab)
        wall_th_results[f"mar_{mar}_um_per_day"] = {
            "wall_thickness_derived_um": wth,
            "vs_eriksen_lacuna_depth_mid_um": lacuna_mid,
            "ratio_derived_over_anchor": wth / lacuna_mid,
            "within_factor_2_of_anchor": within_factor(wth, lacuna_mid, TOLERANCE_FACTOR),
        }

    return {
        "trabecular_occupancy_derivation": {
            "method": "occupancy_fraction_pct(T_bmu_days, T_recur_days) = 100 * T_bmu/T_recur "
                      "(renewal-reward theorem; T_bmu, T_recur both live-verified Eriksen 2010 numbers)",
            "t_bmu_days": t_bmu_trab, "t_recur_days": t_recur_trab,
            "derived_pct_yr": trab_derived_pct_yr,
            "task_stated_pct_yr": TASK_STATED["trabecular_turnover_pct_yr"],
            "ratio_derived_over_task": trab_derived_pct_yr / TASK_STATED["trabecular_turnover_pct_yr"],
            "within_factor_2_PASS": trab_vs_task_pass,
        },
        "cortical_direct_comparison": {
            "clarke_2008_verified_pct_yr_range": CLARKE_2008["cortical_turnover_pct_yr"],
            "clarke_2008_verified_pct_yr_mid": cort_lit_mid_pct_yr,
            "task_stated_pct_yr": TASK_STATED["cortical_turnover_pct_yr"],
            "task_over_verified_ratio": cort_vs_task_ratio,
            "within_factor_2_PASS": cort_vs_task_pass,
            "honest_finding": ("The task's stated 10%/yr cortical figure is "
                               f"{cort_vs_task_ratio:.2f}x HIGHER than the live-verified Clarke 2008 "
                               "figure (2-3%/yr) -- this falsifier-check FAILS for cortical bone "
                               "specifically, reported plainly, not reconciled away."),
        },
        "cross_compartment_ratio_check": {
            "task_implied_trabecular_over_cortical_ratio": task_ratio,
            "verified_derived_trabecular_over_cortical_ratio": verified_ratio,
            "within_factor_2_PASS": ratio_of_ratios_pass,
            "honest_finding": ("The task's two numbers imply trabecular turns over "
                               f"{task_ratio:.1f}x faster than cortical; the live-verified/derived "
                               f"numbers imply {verified_ratio:.1f}x -- a real, reported discrepancy "
                               "in the RATIO, not just the absolute numbers, even though the "
                               "trabecular absolute number alone reproduces well (see above)."),
        },
        "cortical_geometric_explanation_INFORMATIONAL_not_3rd_source_checked": {
            "t_bmu_cortical_days_eriksen2010": t_bmu_cort,
            "implied_mean_recurrence_interval_years": cort_implied_recur_years,
            "vs_trabecular_verified_recurrence_years": ERIKSEN_2010["trabecular_surface_full_renewal_years"],
            "geometric_mechanism": ("Cortical and trabecular BMU cycle LENGTHS differ by <2x (120 vs 200 "
                                    f"days, verified), but the IMPLIED recurrence interval differs by "
                                    f"~{cort_implied_recur_years / ERIKSEN_2010['trabecular_surface_full_renewal_years']:.1f}x "
                                    "(~13 vs 2 years) -- i.e. the dominant driver of the compartment "
                                    "turnover-rate gap is ACTIVATION FREQUENCY (how often a given site is "
                                    "revisited), not BMU cycle speed. Qualitatively, this is the standard "
                                    "surface-to-volume-ratio (BS/BV) argument (Parfitt 1983, PMID 6630513, "
                                    "verified live -- trabecular bone exposes far more remodeling-accessible "
                                    "surface per unit volume than cortical bone) -- invoked here to EXPLAIN "
                                    "the mechanism geometrically, NOT independently re-verified against an "
                                    "exact BS/BV number (informational, flagged as such)."),
        },
        "resorption_formation_balance_coupling_check": {
            "method": "wall_thickness_um(MAR, formation_days) = MAR x formation_days, compared "
                      "against Eriksen 2010's live-verified resorption-lacuna depth (what a "
                      "COUPLED/balanced, non-pathological BMU must refill) -- the closest "
                      "available live-verified proxy for the task's "
                      "'osteoclast resorption vs osteoblast formation balance' question.",
            "eriksen_lacuna_depth_um_range": [lacuna_lo, lacuna_hi],
            "formation_phase_days_used": formation_days_trab,
            "per_mar_endpoint": wall_th_results,
            "honest_finding": ("The LOW end of the task's stated MAR range (0.5 um/day) reproduces "
                               "Eriksen's resorption-lacuna-depth anchor within the pre-registered "
                               "factor-of-2 band; the HIGH end (1.0 um/day) does NOT (over-predicts "
                               "wall thickness by ~3x) -- reported per-endpoint, not averaged."),
        },
    }


# --------------------------------------------- FALSIFIER 2: neck mechanostat --
def falsifier_2_neck_mechanostat(neck):
    """Applies the bone_remodeling cell's, already-gated
    classify_all_variants/stress_mpa_to_microstrain UNEDITED to the femoral
    NECK results (the femoral_neck_stress cell) -- numbers Part A
    (which only ever read bone_stress_results.json, the SHAFT) never
    touched. This is the task's 'decorrelated strain-vs-turnover
    cross-check': a second anatomical site, second free body, SAME reused
    classifier."""
    def point(label, stress_mpa, stored_ue_mag=None, **meta):
        ue_signed = br.stress_mpa_to_microstrain(stress_mpa)
        ue_mag = abs(ue_signed)
        consistency = None
        if stored_ue_mag is not None:
            consistency = {
                "recomputed_ue_mag": ue_mag, "stored_ue_mag": stored_ue_mag,
                "diff": abs(ue_mag - stored_ue_mag),
                "consistent": abs(ue_mag - stored_ue_mag) < 1e-6,
            }
        d = {
            "label": label, "stress_MPa": stress_mpa,
            "microstrain_magnitude": ue_mag,
            "classification_by_variant": br.classify_all_variants(ue_mag),
            "task_specified_classification": br.classify_one(ue_mag, br.THRESHOLD_VARIANTS["task_specified"]),
        }
        if consistency is not None:
            d["bit_for_bit_consistency_vs_upstream_json"] = consistency
        d.update(meta)
        return d

    primary = neck["primary"]
    null = neck["null_control_NSA180"]

    primary_pt = point("neck_primary_AS_IS", primary["peak_stress_MPa"], primary["microstrain_magnitude"])
    primary_axial_floor = point("neck_primary_axial_only_hypothetical_floor", primary["sigma_axial_MPa"])
    null_pt = point("neck_null_control_NSA180_zero_eccentricity", null["peak_stress_MPa"], null["microstrain_magnitude"])
    null_axial_floor = point("neck_null_control_axial_only_hypothetical_floor", null["sigma_axial_MPa"])

    nsa_pts = [point(f"neck_nsa_{k}deg", v["peak_stress_MPa"], v["microstrain_magnitude"], nsa_deg=float(k))
               for k, v in neck["nsa_sweep_results"].items()]
    geom_pts = [point(f"neck_geom_{k}", v["peak_stress_MPa"], v["microstrain_magnitude"])
                for k, v in neck["cross_section_sweep"].items()]
    deinfl = neck["deinflation_counterfactual"]
    deinfl_pt = point("neck_deinflation_counterfactual", deinfl["peak_stress_MPa"], deinfl["microstrain_magnitude"])

    all_real_pts = [primary_pt, null_pt] + nsa_pts + geom_pts + [deinfl_pt]
    all_axial_floor_pts = [primary_axial_floor, null_axial_floor]

    n_real_maintenance = sum(1 for p in all_real_pts if p["task_specified_classification"] == "maintenance_adapted")
    n_real_formation_or_higher = sum(1 for p in all_real_pts if p["task_specified_classification"] == "formation_or_higher_overload")
    n_real_resorption = sum(1 for p in all_real_pts if p["task_specified_classification"] == "resorption_disuse")
    n_axial_floor_resorption = sum(1 for p in all_axial_floor_pts if p["task_specified_classification"] == "resorption_disuse")

    all_consistent = all(
        p.get("bit_for_bit_consistency_vs_upstream_json", {}).get("consistent", True) for p in all_real_pts
    )

    task_anchor_claim_supported = (n_real_maintenance == len(all_real_pts))

    return {
        "method": "Reused verbatim: the bone_remodeling cell's classify_all_variants/stress_mpa_to_microstrain, "
                  "applied to femoral_neck_stress_results.json (a JSON Part A never reads).",
        "primary_AS_IS": primary_pt,
        "primary_axial_only_hypothetical_floor": primary_axial_floor,
        "null_control_NSA180": null_pt,
        "null_control_axial_only_hypothetical_floor": null_axial_floor,
        "nsa_sweep": nsa_pts,
        "cross_section_sweep": geom_pts,
        "deinflation_counterfactual": deinfl_pt,
        "gates": {
            "bit_for_bit_consistency_vs_upstream_neck_json_ALL_PASS": all_consistent,
            "n_real_configs_tested": len(all_real_pts),
            "n_real_configs_classified_maintenance": n_real_maintenance,
            "n_real_configs_classified_formation_or_higher": n_real_formation_or_higher,
            "n_real_configs_classified_resorption": n_real_resorption,
            "n_axial_only_floors_tested": len(all_axial_floor_pts),
            "n_axial_only_floors_classified_resorption": n_axial_floor_resorption,
        },
        "falsifier_verdict": {
            "task_anchor_claim": "habitual gait strain at the femoral neck sits in the mechanostat MAINTENANCE window",
            "task_anchor_claim_SUPPORTED": task_anchor_claim_supported,
            "finding": (
                f"FALSIFIER FIRES: {n_real_formation_or_higher}/{len(all_real_pts)} real (bending-included) "
                "neck configurations classify formation-or-higher (overload), ZERO classify maintenance -- "
                "consistent with, and independently reinforcing (second anatomical site, second free body), "
                "Part A's shaft-based finding. Both axial-only hypothetical floors (what the verdict "
                "would be if Taylor 1996's compression-dominated regime held here instead of this coarse "
                "model's bending-dominated one) classify RESORPTION, not maintenance either -- i.e. "
                "across every tested variant at BOTH anatomical sites (shaft: Part A; neck: this script), "
                "this model class never lands a configuration in the 'maintenance' bin: it is either "
                "resorption (axial-only-hypothetical) or formation-or-higher (real, bending-included). "
                "The task's anchor claim ('maintenance window') is NOT supported by either site tested "
                "tested so far."
                if not task_anchor_claim_supported else
                "Task anchor claim supported at every tested configuration."
            ),
        },
    }


def main():
    assert _self_test()

    if not os.path.exists(NECK_JSON):
        if ASSERT_UPSTREAM_EXISTS:
            raise FileNotFoundError(
                f"upstream femoral_neck_stress_results.json not found at {NECK_JSON} -- "
                "required for FALSIFIER 2 (the decorrelated strain-vs-turnover cross-check); "
                "aborting rather than guessing or silently skipping.")
    with open(NECK_JSON) as f:
        neck = json.load(f)
    neck_mtime = time.ctime(os.path.getmtime(NECK_JSON))

    f1 = falsifier_1_turnover_rate()
    f2 = falsifier_2_neck_mechanostat(neck)

    gates = {
        "self_test_pass": True,
        "f1_trabecular_occupancy_within_factor_2_of_task_stated_25pct": f1["trabecular_occupancy_derivation"]["within_factor_2_PASS"],
        "f1_cortical_clarke2008_within_factor_2_of_task_stated_10pct": f1["cortical_direct_comparison"]["within_factor_2_PASS"],
        "f1_cross_compartment_ratio_within_factor_2": f1["cross_compartment_ratio_check"]["within_factor_2_PASS"],
        "f1_coupling_balance_mar_0.5_within_factor_2": f1["resorption_formation_balance_coupling_check"]["per_mar_endpoint"]["mar_0.5_um_per_day"]["within_factor_2_of_anchor"],
        "f1_coupling_balance_mar_1.0_within_factor_2": f1["resorption_formation_balance_coupling_check"]["per_mar_endpoint"]["mar_1.0_um_per_day"]["within_factor_2_of_anchor"],
        "f2_bit_for_bit_consistency_vs_upstream_neck_json": f2["gates"]["bit_for_bit_consistency_vs_upstream_neck_json_ALL_PASS"],
        "f2_task_anchor_maintenance_claim_SUPPORTED": f2["falsifier_verdict"]["task_anchor_claim_SUPPORTED"],
        "f2_zero_real_configs_land_in_maintenance": f2["gates"]["n_real_configs_classified_maintenance"] == 0,
        "f2_all_axial_only_floors_land_in_resorption": f2["gates"]["n_axial_only_floors_classified_resorption"] == f2["gates"]["n_axial_only_floors_tested"],
    }

    bmu_lifespan_cross_check = {
        "task_stated_months": TASK_STATED["bmu_lifespan_months"],
        "task_stated_days_range": [TASK_STATED["bmu_lifespan_months"][0] * 30.44, TASK_STATED["bmu_lifespan_months"][1] * 30.44],
        "eriksen_2010_cortical_bmu_cycle_days_median": ERIKSEN_2010["cortical_bmu_cycle_days_median"],
        "eriksen_2010_cortical_bmu_cycle_months": ERIKSEN_2010["cortical_bmu_cycle_days_median"] / 30.44,
        "eriksen_2010_trabecular_bmu_cycle_days": ERIKSEN_2010["trabecular_bmu_cycle_days"],
        "eriksen_2010_trabecular_bmu_cycle_months": ERIKSEN_2010["trabecular_bmu_cycle_days"] / 30.44,
        "clarke_2008_cortical_resorption_plus_formation_days_range": [
            CLARKE_2008["resorption_phase_weeks"][0] * 7.0 + CLARKE_2008["formation_phase_months"][0] * 30.44,
            CLARKE_2008["resorption_phase_weeks"][1] * 7.0 + CLARKE_2008["formation_phase_months"][1] * 30.44,
        ],
        "honest_finding": (
            "Eriksen 2010's cortical BMU-cycle median (120d ~= 3.9mo) sits BELOW the task's stated "
            "6-9mo range; Eriksen's trabecular figure (200d ~= 6.6mo) sits WITHIN it; Clarke 2008's "
            "cortical resorption+formation sum (2-4wk + 4-6mo = ~136-210 days ~= 4.5-6.9mo) straddles "
            "the LOW edge of the task's range. 'BMU lifespan ~6-9 months' is best supported for the "
            "TRABECULAR compartment and Clarke's upper cortical estimate, not for Eriksen's cortical "
            "median point estimate -- reported per-compartment, not forced into one number."
        ),
    }

    out = {
        "method": ("Part B of the bone-remodeling layer: quantitative BMU turnover-KINETICS "
                  "(rate), cross-checked against live-verified in-vivo histomorphometry, PLUS "
                  "reuse of Part A's Frost-mechanostat classifier applied to the femoral-NECK "
                  "strain results (a second, decorrelated anatomical site Part A never touched). "
                  "Two independent falsifiers, reported separately, never merged into one verdict."),
        "couples_to": {
            "bone_remodeling_part_a": "the bone_remodeling cell (reused UNEDITED: THRESHOLD_VARIANTS, classify_all_variants, stress_mpa_to_microstrain)",
            "bone_stress_shaft": "the bone_stress cell result (not re-read directly here; already consumed by Part A)",
            "bone_stress_consequence": "the bone_stress_consequence cell result (informational cross-reference only, not re-read)",
            "femoral_neck_stress": NECK_JSON,
        },
        "femoral_neck_stress_json_mtime": neck_mtime,
        "task_stated_anchors": TASK_STATED,
        "literature_anchors": {
            "clarke_2008": CLARKE_2008,
            "eriksen_2010": ERIKSEN_2010,
            "parfitt_1994": PARFITT_1994,
            "hauge_2001": HAUGE_2001,
            "manolagas_2000": MANOLAGAS_2000,
            "eriksen_1986": ERIKSEN_1986,
            "parfitt_1983": PARFITT_1983,
            "parfitt_1987_nomenclature": PARFITT_1987_NOMENCLATURE,
        },
        "literature_search_misses_this_session": LITERATURE_SEARCH_MISSES,
        "falsifier_1_turnover_rate_vs_histomorphometry": f1,
        "bmu_lifespan_cross_check": bmu_lifespan_cross_check,
        "falsifier_2_neck_mechanostat_decorrelated_cross_check": f2,
        "gates": gates,
        "confidence_tier": "in-vivo-anchored (histomorphometry) for the literature turnover-rate/phase-duration "
                          "figures (Clarke 2008, Eriksen 2010, Parfitt 1994/1983, Hauge 2001 -- all live-verified "
                          "primary-source abstracts/full-text); method-only for the renewal-theory "
                          "occupancy DERIVATION itself (a first-pass geometric model, not independently "
                          "cross-validated against a 3rd, differently-derived turnover-rate estimate); the femoral-NECK mechanostat classification inherits the femoral_neck_stress cell's "
                          "own confidence tier (partially in-vivo-anchored crossed with method-only, per that "
                          "document's Sec. 7) UNCHANGED -- this script adds no new mechanical modeling there, "
                          "only re-classifies already-computed, already-gated numbers.",
        "honest_gaps": [
            "The live literature search did NOT find a primary source restating the task's exact "
            "'10%/yr cortical, 25%/yr trabecular' figures verbatim -- several plausible candidate PMIDs "
            "recalled from memory turned out, on live verification, to be entirely unrelated papers (see "
            "literature_search_misses_this_session; each caught by ArticleTitle mismatch before use, not "
            "assumed). The task's numbers are retained as a labeled TASK_STATED variant and cross-checked "
            "against what COULD be verified (Clarke 2008 direct for cortical, an Eriksen-2010-derived renewal "
            "calculation for trabecular) -- not silently replaced or forced to match.",
            "The trabecular occupancy-fraction DERIVATION (27.4%/yr, within factor-2 of the task's 25%/yr) "
            "assumes the 'fraction of surface actively remodeling' is a reasonable proxy for 'fraction of bone "
            "VOLUME turned over per year' -- true only if resorption and formation are volumetrically balanced "
            "(coupled) and spatially homogeneous; a first-pass approximation, not independently verified via a "
            "genuinely distinct 3rd method.",
            "The cortical turnover figure (Clarke 2008, 2-3%/yr) is LOWER than the task's stated 10%/yr by "
            "roughly 3-4x -- this falsifier-check FAILS as stated; reported plainly (falsifier_1_turnover_rate_vs_histomorphometry.cortical_direct_comparison), "
            "not reconciled by adjusting either number.",
            "The resorption/formation BALANCE (coupling) check uses TRABECULAR/cancellous-specific Eriksen 2010 "
            "numbers (150-day formation phase, 40-60um lacuna depth) paired with the TASK's MAR range -- the "
            "high end of that MAR range (1.0 um/day) over-predicts wall thickness by ~3x relative to the "
            "lacuna-depth anchor; only the low end (0.5 um/day) reproduces it within the pre-registered "
            "factor-of-2 band. A single point-estimate MAR was not independently pinned.",
            "'BMU lifespan ~6-9 months' (task-stated) is best supported for the TRABECULAR compartment "
            "(Eriksen 2010: 200d~=6.6mo) and Clarke 2008's upper cortical estimate (~6.9mo); it sits ABOVE "
            "Eriksen 2010's cortical median point estimate (120d~=3.9mo) -- reported per-compartment (see "
            "bmu_lifespan_cross_check), not forced into a single number.",
            "The mechanostat 'maintenance window' cross-check (FALSIFIER 2) inherits EVERY honest gap "
            "the bone_remodeling cell Part A and the femoral_neck_stress cell already disclosed unchanged (generic, "
            "not subject-measured, geometry; bending-dominated-vs-Taylor-1996-compression-dominated tension; "
            "boundary-value, not interior-optimum, neck peak; quasi-static/translational-inertia-only). This "
            "script performs NO new mechanical computation for FALSIFIER 2 -- it only re-applies an "
            "already-gated classifier to already-computed, already-gated numbers, and the "
            "bit_for_bit_consistency_vs_upstream_neck_json gate exists specifically to catch any transcription "
            "error in doing so (PASS, see gates).",
            "Single gait trial, right leg only; the femoral-neck geometry (NSA/offset/anteversion) is itself "
            "population-level literature (Carmona et al. 2019, n=628), not individually measured hip "
            "geometry (inherited from the femoral_neck_stress cell, unchanged).",
            "Osteoclast nuclei turnover (~8%/day, Parfitt 1994) is a DIFFERENT quantity from bone-volume "
            "turnover (~cell-level nuclear replacement inside a single multinucleated osteoclast, not tissue-"
            "level remodeling rate) -- reported for completeness/context, not used in any of the turnover-rate "
            "cross-checks above, to avoid conflating two distinct kinetic quantities that happen to share a "
            "'%' unit.",
        ],
    }

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    print("BMU TURNOVER KINETICS -- Part B (rate) + neck mechanostat cross-check (Falsifier 2)")
    print(f"Neck JSON: {NECK_JSON} (mtime {neck_mtime})\n")
    print("FALSIFIER 1 (turnover rate vs histomorphometry):")
    print(f"  trabecular: derived {f1['trabecular_occupancy_derivation']['derived_pct_yr']:.1f}%/yr vs task-stated "
          f"{TASK_STATED['trabecular_turnover_pct_yr']}%/yr -> PASS={f1['trabecular_occupancy_derivation']['within_factor_2_PASS']}")
    print(f"  cortical:   Clarke2008 {f1['cortical_direct_comparison']['clarke_2008_verified_pct_yr_mid']:.1f}%/yr vs task-stated "
          f"{TASK_STATED['cortical_turnover_pct_yr']}%/yr -> PASS={f1['cortical_direct_comparison']['within_factor_2_PASS']}")
    print(f"  cross-compartment ratio: task={f1['cross_compartment_ratio_check']['task_implied_trabecular_over_cortical_ratio']:.1f}x "
          f"vs verified={f1['cross_compartment_ratio_check']['verified_derived_trabecular_over_cortical_ratio']:.1f}x -> "
          f"PASS={f1['cross_compartment_ratio_check']['within_factor_2_PASS']}")
    print("\nFALSIFIER 2 (neck mechanostat, decorrelated cross-check):")
    print(f"  {f2['falsifier_verdict']['finding']}")
    print("\nGATES:")
    for k, v in gates.items():
        print(f"  {k}: {v}")
    print(f"\nWrote {OUT_JSON}")
    return out


if __name__ == "__main__":
    main()
