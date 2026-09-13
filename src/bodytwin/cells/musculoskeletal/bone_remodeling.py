"""Static strain-to-remodeling-regime map for the femur (Frost mechanostat), the first step above
a bone-stress layer: it reads a bone_stress result JSON (Euler-Bernoulli beam femur stress),
converts sigma -> epsilon via epsilon = sigma/E with E = 17 GPa (Reilly & Burstein 1975 lineage,
the same constant the stress layer uses), and classifies every upstream-swept geometry/offset
configuration against three independently sourced Frost mechanostat threshold variants
(a commonly repeated simplified set, Frost 1987 primary-verified, and a secondary summary).

It is NOT a time-integrated bone-mass simulation: no bone-mass/BMD state variable is advanced
through time; only an instantaneous strain -> regime classification per configuration. Femur
geometry is generic (literature-range), not individually measured.

Reads: OUT_ROOT/bone_stress/bone_stress_results.json.
Writes: bone_remodeling_results.json.

Gates: every real swept configuration must sit at or above the 1500 microstrain formation
threshold, no configuration may classify as maintenance or resorption, and the verdict must hold
across all three threshold variants (with the disclosed caveat that the three variants share the
identical 1500 microstrain formation threshold and differ only on the disuse threshold). An
axial-only hypothetical floor (bending zeroed) is reported alongside, because a bending-dominated
primary configuration contradicts the Taylor et al. 1996 anchor (PMID 8673318: the intact femur is
loaded primarily in compression); that tension is reported, not adjudicated.
"""
import json
import os
import time

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
BONE_STRESS_JSON = _os.path.join(OUT_ROOT, "bone_stress", "bone_stress_results.json")
OUT_DIR = _os.path.join(OUT_ROOT, "bone_remodeling")
OUT_JSON = _os.path.join(OUT_DIR, "bone_remodeling_results.json")

# Reused verbatim from scripts/msk/bone_stress.py's CORTICAL_E_GPA -- NOT
# re-derived here, so the sigma->epsilon conversion is identical between the
# two layers. Same honest gap carried forward (see module docstring).
CORTICAL_E_GPA = 17.0
MICROSTRAIN_PER_MPA = 1.0e6 / (CORTICAL_E_GPA * 1000.0)  # = 1e6 / (E in MPa)

# --------------------------------------------------------------- thresholds --
# Three independently-sourced threshold variants (microstrain). Each is run
# against every computed configuration below (the "force the adversary"
# step: is the qualitative verdict robust to WHICH published threshold set is
# used, not just to the one the task happened to specify?).
#
# Binning convention (half-open, pre-registered before computing any strain
# number below): resorption < r_max <= maintenance < m_max <= formation.
# Values >= m_max are reported as "formation_or_higher": Frost 1987's
# verified abstract phrasing ("in or above the 1500-3000 microstrain range")
# does not itself split that band into two further sub-zones, so this script
# does NOT invent an unverified 4th "pathologic overload" category as a
# primary bin. A hedged, explicitly-labeled secondary split at 3000 ue is
# reported ONLY as a common-secondary-literature annotation, not a verified
# 4th threshold.
THRESHOLD_VARIANTS = {
    "task_specified": {
        "resorption_max_ue": 1000.0,
        "maintenance_max_ue": 1500.0,
        "common_secondary_split_ue": 3000.0,  # hedge only, see docstring
        "source": "Task specification (the build specification).",
        "verified_live_this_session": False,
        "note": ("Round-number simplified version widely repeated in "
                 "secondary/clinical literature; does NOT match Frost 1987's "
                 "own primary-text lower threshold (100-300 ue) -- see "
                 "frost_1987 variant below and honest_gaps."),
    },
    "frost_1987_verified": {
        "resorption_max_ue": 300.0,   # conservative (upper) end of the stated 100-300 range
        "resorption_max_ue_lower_bound": 100.0,
        "maintenance_max_ue": 1500.0,
        "source": ("Frost HM 1987, Anat Rec 219(1):1-9, PMID 3688455, DOI "
                    "10.1002/ar.1092190104 -- abstract fetched verbatim via "
                    "NCBI efetch."),
        "verified_live_this_session": True,
        "quote": ("bone strains in or above the 1500-3000 microstrain range "
                  "cause bone modelling to increase cortical bone mass, "
                  "while strains below the 100-300 microstrain range "
                  "release BMU-based remodeling which then removes existing "
                  "cortical-endosteal and trabecular bone."),
        "note": ("Primary source does not subdivide the >=1500 band further "
                 "-- classified here as a single formation_or_higher zone."),
    },
    "wikipedia_secondary": {
        "resorption_max_ue": 800.0,
        "maintenance_max_ue": 1500.0,
        "source": ("Wikipedia 'Mechanostat' article, raw wikitext fetched "
                    "live (secondary/tertiary, NOT peer-"
                    "reviewed -- used only to triangulate threshold values, "
                    "not as the anchor)."),
        "verified_live_this_session": True,
        "note": ("Same source states habitual/'typical load' for a long "
                 "bone such as the tibia is 2000-3000 microstrain -- i.e. "
                 "BY ITS OWN NUMBERS this sits inside 'overload', not "
                 "inside its own 800-1500 'adapted' band. Quoted verbatim "
                 "in honest_gaps."),
    },
}

ASSERT_UPSTREAM_EXISTS = True  # hard-fail, never silently substitute synthetic data


def stress_mpa_to_microstrain(stress_mpa):
    """epsilon (microstrain) = sigma / E, E in the SAME units as bone_stress.py
    uses (17 GPa cortical bone). Sign preserved: negative = compressive."""
    return stress_mpa * MICROSTRAIN_PER_MPA


def classify_one(microstrain_magnitude, variant):
    r_max = variant["resorption_max_ue"]
    m_max = variant["maintenance_max_ue"]
    if microstrain_magnitude < r_max:
        return "resorption_disuse"
    if microstrain_magnitude < m_max:
        return "maintenance_adapted"
    return "formation_or_higher_overload"


def classify_all_variants(microstrain_magnitude):
    return {name: classify_one(microstrain_magnitude, v)
            for name, v in THRESHOLD_VARIANTS.items()}


def secondary_severity_hedge(microstrain_magnitude):
    """Non-authoritative annotation only (see module docstring): common
    secondary/clinical convention further splits formation_or_higher at
    ~3000 ue into 'mild overload (classic formation)' vs 'severe/pathologic
    overload (microdamage-risk zone)'. NOT independently verified against a
    primary Frost source -- reported as a hedge, not a bin."""
    if microstrain_magnitude < 1500.0:
        return "n/a (below formation threshold)"
    if microstrain_magnitude <= 3000.0:
        return "mild overload (classic formation band, hedge-consistent)"
    return "beyond classic 1500-3000 formation band (hedge: commonly-cited " \
           "'pathologic overload' territory, NOT independently verified " \
           "live as a distinct 4th Frost threshold)"


def make_point(label, stress_mpa, **meta):
    ue_signed = stress_mpa_to_microstrain(stress_mpa)
    ue_mag = abs(ue_signed)
    point = {
        "label": label,
        "stress_MPa": stress_mpa,
        "microstrain_signed": ue_signed,
        "microstrain_magnitude": ue_mag,
        "sign": "compressive" if stress_mpa < 0 else "tensile",
        "classification_by_variant": classify_all_variants(ue_mag),
        "secondary_severity_hedge": secondary_severity_hedge(ue_mag),
    }
    point.update(meta)
    return point


# ------------------------------------------------------------- self-test --
# NOTE on the self-test below (behaviour-preserving rewrite): an earlier version hardcoded the
# literal 17.0 (MPa) and asserted its conversion equals exactly 1000.0, which is only a tautology
# of the sigma/E formula WHEN 17.0 equals CORTICAL_E_GPA; for any other CORTICAL_E_GPA it is a
# false assertion about a fixed literal, so it silently pinned the elastic-modulus parameter and
# crashed every perturbed draw. Likewise the boundary-classification asserts hardcoded
# 999.999/1000.0/1499.999/1500.0/3000.0, pinning the task_specified resorption/maintenance bounds.
# Fix: reference the module's live constants (CORTICAL_E_GPA, v["resorption_max_ue"],
# v["maintenance_max_ue"]). At the baseline (17.0 GPa, 1000.0, 1500.0) this produces identical
# IDENTICAL returned value as the original.
def _self_test():
    """Cheap, decisive, machine-checked correctness gate on the ARITHMETIC
    and BINNING LOGIC only (not on the science) -- run every time, before
    trusting any number below. Uses an exact round-number construction:
    (stress numerically equal to CORTICAL_E_GPA) is ALWAYS exactly 1000
    microstrain, for any positive E -- an algebraic tautology of sigma/E,
    disclosed as such, not a physiological check."""
    ue = stress_mpa_to_microstrain(CORTICAL_E_GPA)
    assert abs(ue - 1000.0) < 1e-9, f"unit conversion self-test FAILED: {ue} != 1000.0"
    v = THRESHOLD_VARIANTS["task_specified"]
    r_max = v["resorption_max_ue"]
    m_max = v["maintenance_max_ue"]
    assert classify_one(r_max - 0.001, v) == "resorption_disuse"
    assert classify_one(r_max, v) == "maintenance_adapted"
    assert classify_one(m_max - 0.001, v) == "maintenance_adapted"
    assert classify_one(m_max, v) == "formation_or_higher_overload"
    assert classify_one(m_max * 2.0, v) == "formation_or_higher_overload"
    return True


def main():
    assert _self_test()

    if not os.path.exists(BONE_STRESS_JSON):
        if ASSERT_UPSTREAM_EXISTS:
            raise FileNotFoundError(
                f"upstream bone_stress_results.json not found at {BONE_STRESS_JSON} -- "
                "per task spec the fallback is the corrected hip/knee contact-force "
                "JSONs, NOT a silent synthetic substitute; this script does not "
                "implement that fallback path because the primary path exists in "
                "this repo. Aborting rather than guessing.")

    with open(BONE_STRESS_JSON) as f:
        bs = json.load(f)

    upstream_mtime = time.ctime(os.path.getmtime(BONE_STRESS_JSON))

    primary = bs["primary"]
    # Medial (compression) side is bs['primary']['peak_stress_MPa'] itself
    # (confirmed: peak_side == 'medial' for this trial/geometry).
    assert primary["peak_side"] == "medial", "unexpected peak side -- re-check sign convention"
    medial_mpa = primary["peak_stress_MPa"]
    axial_mpa = primary["sigma_axial_MPa"]
    bend_amp_mpa = primary["sigma_bend_amp_MPa"]
    lateral_mpa = axial_mpa + bend_amp_mpa  # opposite cortex: axial + bend (not - bend)
    # cross-check the medial reconstruction reproduces the stored peak exactly
    medial_reconstructed = axial_mpa - bend_amp_mpa
    assert abs(medial_reconstructed - medial_mpa) < 1e-6, (
        f"medial reconstruction mismatch: {medial_reconstructed} vs {medial_mpa}")

    bending_dominated_primary = abs(axial_mpa) < bend_amp_mpa

    primary_meta = dict(region=primary["region"], s=primary["s"], t_s=primary["t_s"],
                         pct_gait_cycle=primary["pct_gait_cycle"],
                         offset_mm=primary["offset_mm"], outer_diam_mm=primary["outer_diam_mm"],
                         thickness_mm=primary["thickness_mm"])

    primary_medial = make_point("primary_medial_compression_cortex", medial_mpa, **primary_meta)
    primary_lateral = make_point("primary_lateral_tension_cortex", lateral_mpa, **primary_meta)
    primary_axial_only = make_point(
        "primary_axial_only_hypothetical_floor_if_bending_were_zero", axial_mpa, **primary_meta)

    # ---- offset sweep (each entry carries its own axial/bend split) ----
    offset_map = []
    for off_mm, entry in bs["offset_sweep_results"].items():
        meta = dict(offset_mm=float(off_mm), region=entry["region"], s=entry["s"],
                    t_s=entry["t_s"], pct_gait_cycle=entry["pct_gait_cycle"])
        pk = make_point(f"offset_{off_mm}mm_peak_{entry['peak_side']}", entry["peak_stress_MPa"], **meta)
        pk["bending_dominated"] = abs(entry["sigma_axial_MPa"]) < entry["sigma_bend_amp_MPa"]
        pk["consistent_with_taylor_1996_axial_dominated"] = not pk["bending_dominated"]
        offset_map.append(pk)
    offset0_axial_only = make_point(
        "offset_0mm_axial_only_hypothetical_floor",
        bs["offset_sweep_results"]["0.0"]["sigma_axial_MPa"],
        offset_mm=0.0, region=bs["offset_sweep_results"]["0.0"]["region"])

    # ---- cross-section (geometry) sweep -- NOTE: no s/t_s stored per-entry;
    # inherited from the primary offset per bone_stress.py's source
    # comment ("Cross-section sensitivity (at primary offset ...)"), not
    # independently re-confirmed per-point in this JSON -- stated, not hidden.
    geom_map = []
    for key, entry in bs["geom_sweep_results"].items():
        meta = dict(outer_diam_mm=entry["outer_diam_mm"], thickness_mm=entry["thickness_mm"],
                    offset_mm_inherited_from_primary=primary["offset_mm"],
                    s_and_t_s_inherited_from_primary_NOT_independently_stored=True)
        geom_map.append(make_point(f"geom_{key}_peak_{entry['peak_side']}",
                                    entry["peak_stress_MPa"], **meta))

    # ---- out-of-valid-domain diagnostic (explicitly NOT headline upstream) --
    diag = bs["unrestricted_whole_domain_diagnostic_NOT_headline"]
    out_of_domain_point = make_point(
        "out_of_valid_domain_femoral_head_neck_diagnostic_NOT_HEADLINE",
        diag["peak_stress_MPa"], region=diag["region"], s=diag["s"], t_s=diag["t_s"])

    # ------------------------------------------------------------- gates --
    all_real_configs = [primary_medial, primary_lateral] + offset_map + geom_map
    all_ue_mag = [p["microstrain_magnitude"] for p in all_real_configs]
    min_ue, max_ue = min(all_ue_mag), max(all_ue_mag)
    # offset=40mm (== primary) and geom od27.0_th6.5 (== primary) are the SAME
    # underlying computation reappearing across sweeps -- dedupe by stress
    # value (6dp) for an honest "distinct configurations" count, separate
    # from the raw (duplicate-including) sweep-entry count.
    n_distinct_configs = len({round(p["stress_MPa"], 6) for p in all_real_configs})

    gate_all_at_or_above_task_formation_threshold = (
        min_ue >= THRESHOLD_VARIANTS["task_specified"]["maintenance_max_ue"])
    gate_zero_configs_maintenance_or_resorption_task = all(
        p["classification_by_variant"]["task_specified"] == "formation_or_higher_overload"
        for p in all_real_configs)
    gate_robust_across_all_3_variants = all(
        all(p["classification_by_variant"][v] == "formation_or_higher_overload" for v in THRESHOLD_VARIANTS)
        for p in all_real_configs)
    # HONEST CAVEAT ON THE ABOVE GATE: all 3 variants share the identical
    # maintenance_max_ue=1500.0 -- they were NOT found to disagree about the
    # upper (formation-onset) threshold, only about the lower (disuse)
    # threshold. So "robust across 3 variants" here specifically means
    # "robust to the disuse-threshold disagreement", which is irrelevant
    # since no configuration is anywhere near that lower boundary -- it does
    # NOT mean 3 independently-conflicting formation-threshold hypotheses
    # were all tried and all agreed. Stated plainly, not oversold.
    variants_share_identical_formation_threshold = len(
        {v["maintenance_max_ue"] for v in THRESHOLD_VARIANTS.values()}) == 1
    gate_primary_bending_dominated_vs_taylor1996 = bending_dominated_primary
    gate_axial_only_floor_would_be_resorption_task = (
        primary_axial_only["classification_by_variant"]["task_specified"] == "resorption_disuse")

    gates = {
        "self_test_pass": True,
        "medial_reconstruction_exact_match": True,
        "ALL_real_swept_configs_at_or_above_task_formation_threshold_1500ue": gate_all_at_or_above_task_formation_threshold,
        "ZERO_real_configs_classify_maintenance_or_resorption_under_task_thresholds": gate_zero_configs_maintenance_or_resorption_task,
        "verdict_ROBUST_across_all_3_independent_threshold_variants": gate_robust_across_all_3_variants,
        "CAVEAT_all_3_variants_share_the_identical_1500ue_formation_threshold": variants_share_identical_formation_threshold,
        "CAVEAT_note": ("the 3 variants were found to DISAGREE only on the lower (disuse) threshold "
                         "(100-300 vs 800 vs 1000 ue) -- none disputed the ~1500 ue formation-onset "
                         "threshold, so 'robust across variants' means robust to that disuse-threshold "
                         "disagreement (which turned out to be moot here, since no configuration is "
                         "anywhere near it), NOT evidence from 3 independently-conflicting "
                         "formation-threshold hypotheses that all happened to agree."),
        "min_microstrain_in_sweep": min_ue,
        "max_microstrain_in_sweep": max_ue,
        "n_raw_sweep_entries_with_duplicates": len(all_real_configs),
        "n_distinct_configurations_deduped_by_stress_value": n_distinct_configs,
        "note_on_duplicates": ("offset=40mm (offset sweep) and OD=27.0mm/th=6.5mm (geometry sweep) are "
                               "both exactly the primary configuration re-appearing across sweeps -- "
                               "counted once in n_distinct_configurations_deduped_by_stress_value."),
        "primary_config_bending_dominated_axial_lt_bend_amp": gate_primary_bending_dominated_vs_taylor1996,
        "primary_config_CONSISTENT_with_Taylor1996_compression_dominated": not gate_primary_bending_dominated_vs_taylor1996,
        "axial_only_hypothetical_floor_would_classify_resorption_under_task_thresholds": gate_axial_only_floor_would_be_resorption_task,
    }

    headline = {
        "which_window_does_gait_strain_sit_in": (
            f"FORMATION_OR_HIGHER (mild-overload / modeling), NOT maintenance and NOT "
            f"resorption -- robust across all {n_distinct_configs} distinct upstream-swept "
            f"geometry/offset configurations ({len(all_real_configs)} raw sweep entries "
            f"including 2 exact duplicates, see gates.note_on_duplicates) and across all "
            f"3 independently-sourced threshold variants (task-specified, "
            f"Frost-1987-verified-primary, Wikipedia-secondary)."
        ),
        "headline_number": (
            f"primary configuration (subtrochanteric shaft, s={primary['s']:.3f}, "
            f"{primary['pct_gait_cycle']:.0f}% gait cycle): medial(compression) = "
            f"{primary_medial['microstrain_magnitude']:.0f} ue, "
            f"lateral(tension) = {primary_lateral['microstrain_magnitude']:.0f} ue -- "
            f"both >= the 1500 ue formation threshold used by every variant checked."
        ),
        "sweep_range": f"{min_ue:.0f} to {max_ue:.0f} ue across the full upstream sweep; "
                       f"the single most conservative real configuration (offset=0mm, "
                       f"i.e. zero hip eccentricity) still reads "
                       f"{offset_map[[o['offset_mm'] for o in offset_map].index(0.0)]['microstrain_magnitude']:.0f} ue"
                       f" -- ~5% ABOVE the 1500 ue formation threshold, not below it.",
        "CRITICAL_CAVEAT": (
            "This 'formation, not maintenance' verdict is likely biased HIGH, not a "
            "settled biological finding: the SAME upstream model is bending-dominated "
            f"(|sigma_axial|={abs(axial_mpa):.2f} MPa vs sigma_bend_amp={bend_amp_mpa:.2f} MPa "
            "at the primary configuration), which directly contradicts the verified "
            "Taylor et al. 1996 anchor (PMID 8673318: physiological muscle+joint loading "
            "loads the intact femur PRIMARILY IN COMPRESSION, not bending). If the "
            "axial-dominated regime Taylor 1996 found held here instead, the pure-axial "
            f"component alone ({abs(axial_mpa)*MICROSTRAIN_PER_MPA:.0f} ue) would classify "
            "as RESORPTION under every threshold variant checked -- i.e. the true, "
            "unresolved answer plausibly sits somewhere between this model's 'formation' "
            "estimate and a lower 'maintenance' estimate closer to Duda et al. 1998's "
            "subject-specific finite-element result (<2000 ue, all-muscle case), and this "
            "first-step layer cannot adjudicate between them -- it inherits, rather than "
            "resolves, the upstream layer's already-flagged reliance on a "
            "literature-range (not subject-measured) neck-shaft offset."
        ),
        "on_the_tasks_own_anchor_claim": (
            "'Habitual gait strain sits in the maintenance window' is NOT unambiguously "
            "supported by what could be verified live. Frost 1987's "
            "abstract does not make a habitual-gait-specific claim (it states the two "
            "threshold RANGES, not where gait sits). The Wikipedia secondary summary "
            "explicitly states typical/habitual peak long-bone strain is 2000-3000 ue -- "
            "which, by that SAME source's 4-window table, is inside its 'overload' "
            "band, not its 'adapted'(maintenance) band. This is consistent with the "
            "mechanostat's logic (a cost-minimizing feedback controller should tune "
            "habitual peak strain to sit NEAR the modeling threshold, not deep inside a "
            "quiescent zone, to avoid wastefully over-building bone) but it means this "
            "task's stated anchor should be treated as a common simplification, not a "
            "verified empirical fact -- flagged rather than silently assumed."
        ),
    }

    out = {
        "method": ("Static strain-to-remodeling-signal map (Frost mechanostat). Reads "
                   "bone_stress_results.json (Euler-Bernoulli beam femur stress) "
                   "verbatim, converts sigma->epsilon via epsilon=sigma/E "
                   "(E=17 GPa, same constant as bone_stress.py), classifies every "
                   "upstream-swept configuration against 3 independently-sourced Frost "
                   "mechanostat threshold variants. NOT a time-integrated bone-mass "
                   "simulation; generic (not subject-measured) femur geometry; first "
                   "step."),
        "upstream_bone_stress_json": BONE_STRESS_JSON,
        "upstream_bone_stress_json_mtime": upstream_mtime,
        "upstream_reused_not_recomputed": True,
        "elastic_modulus_GPa": CORTICAL_E_GPA,
        "elastic_modulus_source_note": (
            "Reused verbatim from bone_stress.py's CORTICAL_E_GPA constant -- a "
            "widely-cited figure attributed to Reilly & Burstein 1975 (PMID 1206042, "
            "existence-verified live, DOI 10.1016/0021-9290(75)90075-5) "
            "but not independently re-extracted from that primary source (PubMed has "
            "no abstract on file for it) -- same honest gap as upstream, not resolved "
            "here."),
        "threshold_variants": THRESHOLD_VARIANTS,
        "primary": {
            "medial_compression_cortex": primary_medial,
            "lateral_tension_cortex": primary_lateral,
            "axial_only_hypothetical_floor": primary_axial_only,
        },
        "offset_sweep_map": offset_map,
        "offset0_axial_only_hypothetical_floor": offset0_axial_only,
        "geom_sweep_map": geom_map,
        "out_of_valid_domain_diagnostic_NOT_headline": out_of_domain_point,
        "gates": gates,
        "headline": headline,
        "literature_anchors": {
            "frost_1987": {
                "citation": "Frost HM. Bone \"mass\" and the \"mechanostat\": a proposal. Anat Rec. 1987 Sep;219(1):1-9.",
                "pmid": "3688455", "doi": "10.1002/ar.1092190104",
                "verified_live_this_session_via": "NCBI efetch (abstract text) + Crossref REST (bibliographic fields) -- both queried directly, cross-matching title/journal/vol/issue/pages/year/author.",
                "quote": THRESHOLD_VARIANTS["frost_1987_verified"]["quote"],
            },
            "frost_2003_update": {
                "citation": "Frost HM. Bone's mechanostat: a 2003 update. Anat Rec A Discov Mol Cell Evol Biol. 2003 Dec;275(2):1081-101.",
                "pmid": "14613308", "doi": "10.1002/ar.a.10119",
                "verified_live_this_session_via": "NCBI efetch + Crossref REST.",
                "note": "Conceptual corroboration only; abstract has no restated microstrain numbers.",
            },
            "frost_2004_angle_orthod_clinicians": {
                "citation": "Frost HM. A 2003 update of bone physiology and Wolff's Law for clinicians. Angle Orthod. 2004 Feb;74(1):3-15.",
                "pmid": "15038485", "doi": "10.1043/0003-3219(2004)074<0003:AUOBPA>2.0.CO;2",
                "verified_live_this_session_via": "NCBI efetch (abstract text).",
                "note": "Existence/context only; abstract has no restated microstrain numbers; full text not accessed.",
            },
            "wikipedia_mechanostat_secondary": THRESHOLD_VARIANTS["wikipedia_secondary"],
            "duda_1998": {
                "citation": "Duda GN, Heller M, Albinger J, Schulz O, Schneider E, Claes L. Influence of muscle forces on femoral strain distribution. J Biomech. 1998 Sep;31(9):841-6.",
                "pmid": "9802785", "doi": "10.1016/s0021-9290(98)00080-3",
                "verified_live_this_session_via": "NCBI efetch (full abstract text re-fetched and re-confirmed).",
                "quote": ("the model with all thigh muscles showed peak surface strains below 2000 mu epsilon (45% gait cycle). "
                          "Under simplified load regimes surface strains reached values close to 3000 mu epsilon."),
                "carried_from": "bone_stress_results.json literature_anchors.duda_1998",
            },
            "taylor_1996": {
                "citation": "Taylor ME, Tanner KE, Freeman MA, Yettram AL. Stress and strain distribution within the intact femur: compression or bending? Med Eng Phys. 1996 Mar;18(2):122-31.",
                "pmid": "8673318", "doi": "10.1016/1350-4533(95)00031-3",
                "verified_live_this_session_via": "NCBI efetch (full abstract text re-fetched and re-confirmed).",
                "quote": ("The results of this investigation strongly support the hypothesis that the femur is loaded "
                          "primarily in compression, and not bending as previously thought."),
                "used_here_as": "the CONSISTENCY CHECK anchor against this model's bending-dominated primary result.",
            },
            "reilly_burstein_1975": {
                "citation": "Reilly DT, Burstein AH. The elastic and ultimate properties of compact bone tissue. J Biomech. 1975;8(6):393-405.",
                "pmid": "1206042", "doi": "10.1016/0021-9290(75)90075-5",
                "verified_live_this_session_via": "NCBI efetch (existence/title/journal/year only -- no abstract on file).",
                "note": "Source of the E=17 GPa figure per bone_stress.py's comment; GPa value itself not independently re-extracted (no abstract available).",
            },
            "turner_1998_three_rules": {
                "citation": "Turner CH. Three rules for bone adaptation to mechanical stimuli. Bone. 1998 Nov;23(5):399-407.",
                "pmid": "9823445", "doi": "10.1016/s8756-3282(98)00118-5",
                "verified_live_this_session_via": "NCBI efetch (full abstract text fetched).",
                "quote": ("(1) It is driven by dynamic, rather than static, loading. (2) Only a short duration of "
                          "mechanical loading is necessary to initiate an adaptive response. (3) Bone cells accommodate "
                          "to a customary mechanical loading environment, making them less responsive to routine loading signals."),
                "used_here_as": "honest-gap anchor: this script's static peak-magnitude-only map omits rate- and novelty-dependence.",
            },
        },
        "honest_gaps": [
            "STATIC map, not a time-integrated remodeling simulation: no bone-mass/BMD/cortical-area state variable is "
            "advanced through simulated time; only an instantaneous strain->regime classification at each upstream-swept "
            "configuration's peak. No remodeling RATE is modeled.",
            "Generic (literature-range, not this subject's CT/DXA-measured) femur geometry and neck-shaft offset -- "
            "inherited unchanged from bone_stress_results.json's honest_gaps; this IS the parameter the "
            "remodeling verdict is most sensitive to (see gates.min/max_microstrain_in_sweep spanning "
            f"{min_ue:.0f}-{max_ue:.0f} ue, a {max_ue/min_ue:.1f}x range, from geometry/offset assumptions alone).",
            "The task's specified thresholds (1000/1500/3000 ue) are a commonly-repeated simplified/secondary-"
            "literature version, NOT verbatim Frost 1987 (whose own verified abstract states a 100-300 ue disuse "
            "threshold, 3-10x lower) -- both are reported side by side rather than silently reconciled.",
            "The single-source '4th pathologic-overload bin above 3000 ue' some clinical literature uses was NOT "
            "independently verified against a primary Frost source (both fetched Frost abstracts, "
            "1987 and 2003, describe only the two-threshold/three-zone picture) -- reported only as a hedge "
            "annotation (secondary_severity_hedge), never as a verified bin boundary.",
            "'Habitual gait sits in the maintenance window' (the task's anchor claim) could not be confirmed "
            "live -- if anything the Wikipedia secondary source's numbers place typical/habitual "
            "long-bone strain in the overload band, not the adapted/maintenance band (see headline."
            "on_the_tasks_own_anchor_claim).",
            "CRITICAL: this model's primary configuration is bending-dominated, which directly contradicts the "
            "verified Taylor et al. 1996 anchor (real femur under physiological loading is compression-dominated). "
            "This is the single largest reliability caveat on the 'formation, not maintenance' headline -- flagged "
            "as unresolved, not adjudicated, by this first-step layer (see headline.CRITICAL_CAVEAT).",
            "Mechanostat strain-MAGNITUDE thresholds omit loading RATE, cycle count, and rest-interval effects that "
            "the literature (Turner 1998, PMID 9823445, verified live) shows also govern the real adaptive response; "
            "this map only classifies peak magnitude.",
            "Single gait trial, right leg only -- no inter-trial or between-person variability is "
            "characterized.",
            "The E=17 GPa cortical modulus is reused, not re-verified, from bone_stress.py (Reilly & Burstein 1975 "
            "citation exists; the specific GPa figure was not independently re-extracted from primary text, no "
            "abstract on file at PubMed for that 1975 record).",
            "geom_sweep_map entries have no independently-stored s/t_s (gait-cycle location) of their own in "
            "bone_stress_results.json; they are inherited from the primary offset configuration per bone_stress.py's "
            "own source comment, not re-verified per-point in this script.",
        ],
    }

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    # ------------------------------------------------------------- stdout --
    print("BONE REMODELING (Frost mechanostat) -- first-step static map")
    print(f"Upstream: {BONE_STRESS_JSON} (mtime {upstream_mtime})")
    print(f"E = {CORTICAL_E_GPA} GPa -> {MICROSTRAIN_PER_MPA:.4f} microstrain per MPa\n")
    print(f"PRIMARY config: {primary['region']}, s={primary['s']:.3f}, "
          f"{primary['pct_gait_cycle']:.1f}% gait cycle, offset={primary['offset_mm']}mm, "
          f"OD={primary['outer_diam_mm']}mm, th={primary['thickness_mm']}mm")
    print(f"  medial (compression): {medial_mpa:.2f} MPa -> {primary_medial['microstrain_magnitude']:.0f} ue "
          f"-> {primary_medial['classification_by_variant']}")
    print(f"  lateral (tension):    {lateral_mpa:.2f} MPa -> {primary_lateral['microstrain_magnitude']:.0f} ue "
          f"-> {primary_lateral['classification_by_variant']}")
    print(f"  axial-only floor:     {axial_mpa:.2f} MPa -> {primary_axial_only['microstrain_magnitude']:.0f} ue "
          f"-> {primary_axial_only['classification_by_variant']}")
    print(f"  bending-dominated (|axial|<bend_amp)? {bending_dominated_primary} "
          f"({'CONTRADICTS' if bending_dominated_primary else 'consistent with'} Taylor 1996 PMID 8673318)\n")

    print("Offset sweep (0=null/no hip eccentricity -> 60mm literature upper bound):")
    for p in offset_map:
        print(f"  offset={p['offset_mm']:>5.1f}mm  {p['stress_MPa']:>8.2f} MPa  "
              f"{p['microstrain_magnitude']:>7.0f} ue  {p['classification_by_variant']['task_specified']:<28s} "
              f"taylor1996-consistent={p['consistent_with_taylor_1996_axial_dominated']}")

    print("\nCross-section (geometry) sweep, at primary offset:")
    for p in geom_map:
        print(f"  OD={p['outer_diam_mm']:>5.1f}mm th={p['thickness_mm']:>4.1f}mm  "
              f"{p['stress_MPa']:>8.2f} MPa  {p['microstrain_magnitude']:>7.0f} ue  "
              f"{p['classification_by_variant']['task_specified']}")

    print(f"\nSweep range: {min_ue:.0f}-{max_ue:.0f} ue "
          f"(all {len(all_real_configs)} real configs)")
    print("\nGATES:")
    for k, v in gates.items():
        print(f"  {k}: {v}")

    print("\nHEADLINE:")
    for k, v in headline.items():
        print(f"  [{k}]\n    {v}\n")

    print(f"Wrote {OUT_JSON}")
    return out


if __name__ == "__main__":
    main()
