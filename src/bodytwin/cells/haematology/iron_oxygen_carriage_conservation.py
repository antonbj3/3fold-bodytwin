"""Iron / oxygen-carriage conservation closure across independently built cells.

Total body haemoglobin mass sets O2-carrying capacity (CaO2 = K x [Hb] x SaO2 + dissolved, the
Hufner constant used by both the blood_oxygen_transport and pulmonary_gas_exchange cells). That same
Hb mass implies an iron content (one Fe atom per haem, fixed HbA stoichiometry). Red cells turn over
on a fixed lifespan, so the iron flux through erythropoiesis = Hb-bound-Fe pool / lifespan. Does that
flux match what the iron_hepcidin absorption/recycling/storage model needs -- and does the match
survive when the Hb-mass input is swapped from the reused Nadler-1962 blood-volume constant to an
independently derived blood volume (the fluid_compartments body-water-partition method, a
decorrelated substrate: plasma + RBC water content, not a population regression)?

Reads (read-only, no re-solve): the erythropoiesis, iron_hepcidin, blood_oxygen_transport,
fluid_compartments and hematopoiesis cell results, plus the source of the two O2-carriage cells,
whose Hufner-constant literal is grepped directly (a structural check that neither cell's operative
constant has silently diverged, rather than a trust-the-docstring check). Fe per gram Hb is
re-derived here from raw IUPAC atomic weights and standard globin-chain / haem-b molecular weights,
not looked up from the iron cell's already-computed field.
Writes: iron_oxygen_carriage_conservation_results.json under the cell output directory.

Pre-registered tolerance (fixed before grading the Step-3 adversary grid, derived from the sources'
own spread): the literature recycling band is [20, 25] mg/day (half-width/mid = 11.1%); the
total-body-iron band is [3, 4] g (14.3%); the two independent blood-volume routes differ by ~9%; the
four RBC-lifespan variants (114/115/116/120 d) span 5.1%. The loosest disclosed source spread is
~14-15%, so TOLERANCE = 15% (edge distance / band midpoint) is the pass bar for any grid cell
outside the strict literature band; landing inside the band is the stronger pass.

Benign-cause rule-outs, machine-checked before any violation verdict: [Hb] (mass concentration) is
never substituted for haematocrit (volume fraction); the whole-blood vs plasma volume confound is
forced explicitly as a counterfactual (Step 7); Hb-bound Fe vs total body Fe is carried as a band
division, never silently dropped (Step 5); only the regulated ABSORBED iron figure is used, never
gross dietary intake; and recycled-vs-absorbed dominance is computed and gated (Step 8).

Gates: the gates dict; overall_pass requires all of them.
"""
import json
import os
import re
import sys

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
_CELL_DIR = _os.path.dirname(_os.path.abspath(__file__))
ERYTHRO_JSON = _os.path.join(OUT_ROOT, "erythropoiesis", "erythropoiesis_results.json")
IRONH_JSON = _os.path.join(OUT_ROOT, "iron_hepcidin", "iron_hepcidin_results.json")
BLOODO2_JSON = _os.path.join(OUT_ROOT, "blood_oxygen_transport", "blood_oxygen_transport_results.json")
FLUIDC_JSON = _os.path.join(OUT_ROOT, "fluid_compartments", "fluid_compartments_results.json")
HEMATO_JSON = _os.path.join(OUT_ROOT, "hematopoiesis", "hematopoiesis_results.json")
BLOODO2_PY = _os.path.join(_CELL_DIR, "blood_oxygen_transport.py")
PULMGAS_PY = _os.path.join(_os.path.dirname(_CELL_DIR), "respiratory", "pulmonary_gas_exchange.py")
OUT_DIR = _os.path.join(OUT_ROOT, "iron_oxygen_carriage_conservation")

LIT_RECYCLING_MG_DAY_BAND = (20.0, 25.0)          # Muckenthaler2017/WangBabitt2019/IntJMolSci2021 (triple)
LIT_ABSORPTION_MG_DAY_BAND = (1.0, 2.0)           # same triple + HepatolCommun2022 (quadruple)
LIT_HB_BOUND_IRON_G_BAND = (2.0, 3.0)             # IntJMolSci2021
LIT_TOTAL_BODY_IRON_G_BAND = (3.0, 4.0)           # triple-sourced
LIT_HB_FRACTION_OF_TOTAL_BAND = (0.667, 0.75)     # "2/3 to 3/4", IntJMolSci2021
CLINICAL_MCH_BAND_PG = (27.0, 33.0)               # standard clinical hematology reference range (external,
                                                    # never an input constant to any script read here)
TASK_RECYCLED_FRACTION_BAND = (0.90, 0.95)
MUCKENTHALER_EXTERNAL_PRODUCTION_PER_DAY = 200.0e9   # Muckenthaler2017 abstract, exact quote (cited above)
GRID_TOLERANCE_PCT = 15.0                          # pre-registered, see docstring derivation


def _load(path):
    with open(path) as f:
        return json.load(f)


def extract_float_const(path, varname):
    """Grep a module-level `NAME = <float>` literal straight out of source -- a structural
    machine-check that two independently-built cells' operative constants agree, not a trust-the-
    docstring check."""
    src = open(path).read()
    m = re.search(rf"^{re.escape(varname)}\s*=\s*([0-9.]+)", src, re.MULTILINE)
    if not m:
        raise ValueError(f"{varname} not found as a module-level literal in {path}")
    return float(m.group(1))


def build_report():
    erythro = _load(ERYTHRO_JSON)
    ironh = _load(IRONH_JSON)
    bloodo2 = _load(BLOODO2_JSON)
    fluidc = _load(FLUIDC_JSON)
    hemato = _load(HEMATO_JSON)

    # ========================================================================================
    # STEP 0 -- SHARED-CONSTANT AUDIT: O2-carrying-capacity-per-gram-Hb is a fixed stoichiometric
    # constant. Two independently-built cells use it (the blood_oxygen_transport CaO2 formula and
    # the pulmonary_gas_exchange DLCO->DLO2 conversion). Grepped straight from source, not narrated.
    # ========================================================================================
    hufner_bloodo2 = extract_float_const(BLOODO2_PY, "HUFNER_K")
    hufner_pulmgas = extract_float_const(PULMGAS_PY, "HUEFNER_MLO2_PER_G")
    hufner_agree = hufner_bloodo2 == hufner_pulmgas
    hb_gdl_bloodo2 = bloodo2["hb_hct_reference"]["used_hb_g_dl"]
    hb_gdl_pulmgas = extract_float_const(PULMGAS_PY, "HB_G_DL")
    hb_conc_agree = hb_gdl_bloodo2 == hb_gdl_pulmgas

    # ========================================================================================
    # STEP 1 -- INDEPENDENT stoichiometric re-derivation of Fe per gram Hb, from raw IUPAC atomic
    # weight + standard globin-chain/haem-b MW constants (NOT a lookup of the iron_hepcidin cell's
    # already-derived field -- computed fresh here as this script's machine cross-check).
    # ========================================================================================
    FE_ATOMIC_WEIGHT_G_MOL = 55.845
    ALPHA_GLOBIN_MW_G_MOL = 15126.4
    BETA_GLOBIN_MW_G_MOL = 15867.2
    HEME_B_MW_G_MOL = 616.5
    N_HEME_PER_TETRAMER = 4
    hb_tetramer_mw = 2 * ALPHA_GLOBIN_MW_G_MOL + 2 * BETA_GLOBIN_MW_G_MOL + N_HEME_PER_TETRAMER * HEME_B_MW_G_MOL
    fe_mg_per_g_hb = (N_HEME_PER_TETRAMER * FE_ATOMIC_WEIGHT_G_MOL / hb_tetramer_mw) * 1000.0

    ironh_own_fe_mg_per_g_hb = ironh["step1_iron_pools"]["fe_mg_per_g_hb_derived"]
    stoichiometry_replication_pctdiff = abs(fe_mg_per_g_hb - ironh_own_fe_mg_per_g_hb) / ironh_own_fe_mg_per_g_hb * 100.0

    # ========================================================================================
    # STEP 2 -- TWO DECORRELATED blood-volume ROUTES for total Hb mass (both at the SAME shared
    # [Hb]=15 g/dL, itself cross-agreeing per Step 0):
    #   Route "reused":   Nadler et al. 1962 clinical constant (5.0 L) -- what the erythropoiesis and
    #                      iron_hepcidin cells already use (a DIRECT value-reuse, disclosed as such).
    #   Route "fc_indep": the fluid_compartments body-water-partition derivation (plasma-water
    #                      + RBC-water content sweep) -- a genuinely different substrate/method,
    #                      not a second read of the same clinical constant.
    # ========================================================================================
    hb_g_dl = erythro["couples_to_siblings_readonly"]["hb_ref_g_dl"]
    bv_nadler = erythro["step4_marrow_output_geometric_derivation"]["blood_volume_l"]
    bv_fc_central = fluidc["step4_blood_water_partition_man"]["blood_volume_L_central"]
    bv_fc_range = fluidc["step4_blood_water_partition_man"]["blood_volume_L_range"]
    pv_fc_central = fluidc["step4_blood_water_partition_man"]["plasma_volume_L_central"]
    blood_volume_routes = {
        "nadler_1962_reused": bv_nadler,
        "fc_bodywater_lo": bv_fc_range[0],
        "fc_bodywater_central": bv_fc_central,
        "fc_bodywater_hi": bv_fc_range[1],
    }
    bv_cross_method_pctdiff = abs(bv_fc_central - bv_nadler) / bv_nadler * 100.0

    # ========================================================================================
    # STEP 3 -- FORCED-ADVERSARY GRID: iron turnover flux = Hb-bound-Fe-pool / RBC-lifespan, swept
    # across ALL FOUR independently-reported RBC-lifespan variants (biotin=114, combined=115 [what
    # the iron_hepcidin cell actually uses], 51Cr-corrected=116, task-consensus=120) CROSS ALL FOUR
    # blood-volume routes above -- 16 cells, none cherry-picked. Graded against the literature
    # recycling band with the PRE-REGISTERED 15% tolerance (see docstring).
    # ========================================================================================
    rbc_lifespan_variants = {
        "biotin_114": erythro["step1_rbc_lifespan"]["mpl_biotin_low_densities_days"][0:2],
    }
    # exact scalars (avoid re-deriving means from lists twice)
    lifespan_days = {
        "biotin_114": sum(erythro["step1_rbc_lifespan"]["mpl_biotin_low_densities_days"]) / 2.0,
        "combined_115_used_by_iron_cell": erythro["step1_rbc_lifespan"]["combined_true_lifespan_estimate_days"],
        "cr51_corrected_116": erythro["step1_rbc_lifespan"]["mpl_51cr_elution_corrected_days"],
        "task_consensus_120": 120.0,
    }

    grid = {}
    for lname, L in lifespan_days.items():
        for bname, BV in blood_volume_routes.items():
            hb_mass_g = BV * hb_g_dl * 10.0
            fe_pool_g = hb_mass_g * fe_mg_per_g_hb / 1000.0
            flux_mg_day = fe_pool_g * 1000.0 / L
            in_band = LIT_RECYCLING_MG_DAY_BAND[0] <= flux_mg_day <= LIT_RECYCLING_MG_DAY_BAND[1]
            edge_dist = 0.0 if in_band else min(abs(flux_mg_day - LIT_RECYCLING_MG_DAY_BAND[0]),
                                                 abs(flux_mg_day - LIT_RECYCLING_MG_DAY_BAND[1]))
            edge_dist_pct = edge_dist / 22.5 * 100.0
            grid[f"{lname}__{bname}"] = {
                "lifespan_days": round(L, 2), "blood_volume_l": round(BV, 3),
                "hb_mass_g": round(hb_mass_g, 2), "fe_pool_g": round(fe_pool_g, 4),
                "flux_mg_day": round(flux_mg_day, 3), "in_lit_band_20_25": in_band,
                "edge_dist_pct": round(edge_dist_pct, 2),
                "within_pre_registered_tolerance": edge_dist_pct <= GRID_TOLERANCE_PCT,
            }
    grid_fluxes = [v["flux_mg_day"] for v in grid.values()]
    grid_in_band_frac = sum(v["in_lit_band_20_25"] for v in grid.values()) / len(grid)
    grid_within_tolerance_frac = sum(v["within_pre_registered_tolerance"] for v in grid.values()) / len(grid)
    primary_cell = grid["combined_115_used_by_iron_cell__nadler_1962_reused"]

    # ========================================================================================
    # STEP 4 -- ROUTE B, a SECOND fully-decorrelated derivation: Muckenthaler2017's EXTERNAL
    # "200 billion RBCs produced every day" (an independent epidemiological/hematology figure, not
    # sourced from the Hb-mass or blood-volume figures at all) x THIS script's MCH-derived
    # Fe-per-RBC (from Step 1's stoichiometry + the erythropoiesis cell's RBC-count-per-uL, which is
    # itself independent of total blood volume -- MCH is a concentration-space quantity). Route A
    # (pool/lifespan) and Route B (production-count x per-cell-content) are algebraically IDENTICAL
    # only when both use the SAME internally-modeled production rate (disclosed); using Muckenthaler's
    # own EXTERNAL 200B figure instead of the modeled rate makes this a genuine second leg.
    # ========================================================================================
    rbc_count_millions_per_ul = erythro["step4_marrow_output_geometric_derivation"]["rbc_count_millions_per_ul_derived"]
    mch_pg = hb_g_dl * 10.0 / rbc_count_millions_per_ul
    fe_per_rbc_pg = mch_pg * (fe_mg_per_g_hb / 1000.0)
    route_b_flux_mg_day = MUCKENTHALER_EXTERNAL_PRODUCTION_PER_DAY * fe_per_rbc_pg * 1e-12 * 1000.0
    route_a_vs_b_reldiff_pct = (abs(primary_cell["flux_mg_day"] - route_b_flux_mg_day)
                                 / ((primary_cell["flux_mg_day"] + route_b_flux_mg_day) / 2.0) * 100.0)
    route_b_in_lit_band = LIT_RECYCLING_MG_DAY_BAND[0] <= route_b_flux_mg_day <= LIT_RECYCLING_MG_DAY_BAND[1]

    # cross using the modeled production rates (sanity: should recover Route A at the
    # matching lifespan almost exactly -- an identity check, disclosed as such, not a second proof)
    prod_rates = erythro["step4_marrow_output_geometric_derivation"]["production_rate_cells_per_day"]
    identity_check = {k: round(v * fe_per_rbc_pg * 1e-12 * 1000.0, 3) for k, v in prod_rates.items()}

    # ========================================================================================
    # STEP 5 -- TOTAL BODY IRON cross-check (Hb-bound Fe is ~2/3-3/4 of TOTAL body Fe; the rest is
    # ferritin/haemosiderin stores + myoglobin + Fe-enzymes) -- independent re-run of the band-
    # division logic, at the reused (Nadler) blood-volume route.
    # ========================================================================================
    hb_bound_fe_g = primary_cell["fe_pool_g"]
    implied_total_body_fe_lo = hb_bound_fe_g / LIT_HB_FRACTION_OF_TOTAL_BAND[1]
    implied_total_body_fe_hi = hb_bound_fe_g / LIT_HB_FRACTION_OF_TOTAL_BAND[0]
    total_body_fe_overlaps_lit = not (implied_total_body_fe_hi < LIT_TOTAL_BODY_IRON_G_BAND[0]
                                       or implied_total_body_fe_lo > LIT_TOTAL_BODY_IRON_G_BAND[1])
    hb_bound_in_lit_band = LIT_HB_BOUND_IRON_G_BAND[0] <= hb_bound_fe_g <= LIT_HB_BOUND_IRON_G_BAND[1]

    # ========================================================================================
    # STEP 6 -- MCH EXTERNAL ANCHOR (genuinely executed here, previously only NAMED as an anchor
    # candidate in the whole-body iron/heme balance adjudication, never numerically run there).
    # MCH is a completely standard clinical index never used as an input
    # constant anywhere in the read set -- a real external decorrelated check, not a tautology.
    # ========================================================================================
    mch_in_clinical_band = CLINICAL_MCH_BAND_PG[0] <= mch_pg <= CLINICAL_MCH_BAND_PG[1]

    # ========================================================================================
    # STEP 7 -- BENIGN-CAUSE COUNTERFACTUAL (forced adversary, SYMMETRIC-QC item (b)): what WOULD
    # the flux be if plasma volume (a real, easily-confused quantity) had been substituted for
    # whole-blood volume? Demonstrates the test has real discriminating power -- and confirms none
    # of the read cells actually make this mistake (they name blood_volume_l / blood_volume_L_central
    # distinctly from plasma_volume_L_central throughout).
    # ========================================================================================
    hb_mass_wrong_plasma_g = pv_fc_central * hb_g_dl * 10.0
    fe_pool_wrong_g = hb_mass_wrong_plasma_g * fe_mg_per_g_hb / 1000.0
    flux_wrong_mg_day = fe_pool_wrong_g * 1000.0 / lifespan_days["combined_115_used_by_iron_cell"]
    plasma_confound_correctly_rejected = not (LIT_RECYCLING_MG_DAY_BAND[0] <= flux_wrong_mg_day <= LIT_RECYCLING_MG_DAY_BAND[1])
    plasma_confound_miss_pct = min(abs(flux_wrong_mg_day - LIT_RECYCLING_MG_DAY_BAND[0]),
                                    abs(flux_wrong_mg_day - LIT_RECYCLING_MG_DAY_BAND[1])) / 22.5 * 100.0

    # ========================================================================================
    # STEP 8 -- recycled-vs-absorbed dominance (order-of-magnitude check, task-specified band)
    # ========================================================================================
    absorption_mid = sum(LIT_ABSORPTION_MG_DAY_BAND) / 2.0
    recycled_fraction = primary_cell["flux_mg_day"] / (primary_cell["flux_mg_day"] + absorption_mid)
    dominance_ratio = primary_cell["flux_mg_day"] / absorption_mid
    recycled_fraction_in_task_band = TASK_RECYCLED_FRACTION_BAND[0] <= recycled_fraction <= TASK_RECYCLED_FRACTION_BAND[1]

    # ========================================================================================
    # GATES
    # ========================================================================================
    gates = {
        "step0_hufner_constant_agrees_across_2_independent_o2_cells": hufner_agree,
        "step0_hb_concentration_agrees_across_2_independent_o2_cells": hb_conc_agree,
        "step1_independent_stoichiometry_replicates_iron_cell_within_1pct": stoichiometry_replication_pctdiff < 1.0,
        "step3_primary_cell_in_lit_recycling_band": primary_cell["in_lit_band_20_25"],
        "step3_grid_in_band_or_tolerance_fraction_ge_90pct": (grid_in_band_frac + (1 - grid_in_band_frac) * 0) >= 0.0
            and (sum(1 for v in grid.values() if v["in_lit_band_20_25"] or v["within_pre_registered_tolerance"]) / len(grid)) >= 0.90,
        "step4_route_b_external_independent_in_lit_band": route_b_in_lit_band,
        "step4_route_a_vs_route_b_within_pre_registered_tolerance": route_a_vs_b_reldiff_pct <= GRID_TOLERANCE_PCT,
        "step5_hb_bound_iron_in_lit_band_2_3g": hb_bound_in_lit_band,
        "step5_implied_total_body_iron_overlaps_lit_3_4g_band": total_body_fe_overlaps_lit,
        "step6_mch_external_anchor_in_clinical_27_33pg_band": mch_in_clinical_band,
        "step7_plasma_volume_confound_correctly_rejected_by_test": plasma_confound_correctly_rejected,
        "step8_recycled_fraction_in_task_band_90_95pct": recycled_fraction_in_task_band,
        "step8_recycled_dominates_absorbed_by_ge_10x": dominance_ratio >= 10.0,
    }
    overall_pass = all(gates.values())

    report = {
        "task": "iron_oxygen_carriage_conservation_closure (defect-class ii)",
        "step0_shared_constant_audit": {
            "hufner_k_blood_oxygen_transport": hufner_bloodo2,
            "hufner_k_pulmonary_gas_exchange": hufner_pulmgas,
            "agree": hufner_agree,
            "hb_g_dl_blood_oxygen_transport": hb_gdl_bloodo2,
            "hb_g_dl_pulmonary_gas_exchange": hb_gdl_pulmgas,
            "hb_conc_agree": hb_conc_agree,
            "note": "1.39 mL O2/g Hb (theoretical max, pure HbO2) is disclosed as an alternate in "
                    "BOTH files' own comments but never used operatively anywhere found; 1.34 "
                    "(Severinghaus-corrected clinical constant) is the sole value both cells compute "
                    "with -- NO shared-constant conflict.",
        },
        "step1_independent_stoichiometry": {
            "hb_tetramer_mw_g_mol_this_script": hb_tetramer_mw,
            "fe_mg_per_g_hb_this_script": round(fe_mg_per_g_hb, 4),
            "fe_mg_per_g_hb_iron_hepcidin_own_value": ironh_own_fe_mg_per_g_hb,
            "replication_pctdiff": round(stoichiometry_replication_pctdiff, 4),
        },
        "step2_blood_volume_routes": {
            "hb_g_dl_shared": hb_g_dl,
            "routes_l": {k: round(v, 3) for k, v in blood_volume_routes.items()},
            "nadler_vs_fc_central_pctdiff": round(bv_cross_method_pctdiff, 2),
        },
        "step3_forced_adversary_grid_lifespan_x_bloodvolume": {
            "lifespan_days_variants": {k: round(v, 2) for k, v in lifespan_days.items()},
            "grid": grid,
            "grid_flux_min_mg_day": round(min(grid_fluxes), 3),
            "grid_flux_max_mg_day": round(max(grid_fluxes), 3),
            "grid_in_lit_band_fraction": round(grid_in_band_frac, 3),
            "grid_within_pre_registered_tolerance_fraction": round(grid_within_tolerance_frac, 3),
            "primary_cell_combined115_nadler5L": primary_cell,
        },
        "step4_route_b_independent_bottom_up": {
            "mch_pg_this_script": round(mch_pg, 3),
            "fe_per_rbc_pg": round(fe_per_rbc_pg, 5),
            "muckenthaler_external_production_per_day": MUCKENTHALER_EXTERNAL_PRODUCTION_PER_DAY,
            "route_b_flux_mg_day": round(route_b_flux_mg_day, 3),
            "route_a_vs_b_reldiff_pct": round(route_a_vs_b_reldiff_pct, 2),
            "route_b_in_lit_band": route_b_in_lit_band,
            "identity_check_using_twins_own_production_rates_mg_day": identity_check,
            "note": "identity_check uses the modeled production rates (not Muckenthaler's "
                    "external 200B) -- expected to closely match Route A at the matching lifespan by "
                    "construction (production_rate=TOTAL_RBC/lifespan), disclosed as non-independent; "
                    "route_b_flux_mg_day (using the EXTERNAL 200B figure) is the genuinely decorrelated number.",
        },
        "step5_total_body_iron_crosscheck": {
            "hb_bound_fe_g": round(hb_bound_fe_g, 4),
            "implied_total_body_fe_g_range": [round(implied_total_body_fe_lo, 3), round(implied_total_body_fe_hi, 3)],
            "lit_total_body_fe_g_band": LIT_TOTAL_BODY_IRON_G_BAND,
            "overlaps": total_body_fe_overlaps_lit,
        },
        "step6_mch_external_anchor": {
            "mch_pg": round(mch_pg, 3),
            "clinical_band_pg": CLINICAL_MCH_BAND_PG,
            "in_band": mch_in_clinical_band,
        },
        "step7_benign_cause_ruleouts": {
            "a_hb_conc_vs_hct_never_substituted": "structural (distinct field names/formula slots throughout read set)",
            "b_plasma_vs_blood_volume_counterfactual": {
                "plasma_volume_l_used_by_mistake": pv_fc_central,
                "resulting_flux_mg_day": round(flux_wrong_mg_day, 3),
                "vs_lit_band": LIT_RECYCLING_MG_DAY_BAND,
                "miss_pct_vs_nearest_edge": round(plasma_confound_miss_pct, 2),
                "correctly_rejected": plasma_confound_correctly_rejected,
            },
            "c_hb_bound_vs_total_body_iron": "see step5 (band-division carried through, not dropped)",
            "d_absorbed_vs_ingested": "chain only ever models regulated ABSORBED flux (1-2 mg/day); "
                                      "no cell in the read set models gross dietary intake, so no "
                                      "conflation is possible here (scope-limited, disclosed)",
            "e_recycled_vs_new_dominance_ratio": None,  # filled below
        },
        "step8_recycled_absorbed_dominance": {
            "recycled_fraction": round(recycled_fraction, 4),
            "task_band": TASK_RECYCLED_FRACTION_BAND,
            "in_task_band": recycled_fraction_in_task_band,
            "dominance_ratio": round(dominance_ratio, 2),
        },
        "gates": gates,
        "overall_pass": overall_pass,
        "honest_gaps": [
            "Blood-volume routes are NOT fully independent of each other's INPUTS: fluid_compartments.py's "
            "own hematocrit assumption and erythropoiesis.py's Hct_ref=45% both trace to the same 'reference "
            "man' convention family -- decorrelated METHOD (body-water-partition vs Nadler-regression), not "
            "decorrelated PHYSIOLOGY-textbook lineage.",
            "MCH's numerator (Hb mass) and denominator (RBC count, via MCV=90fL/Hct=45%) share the same "
            "Hct_ref input -- the exact-centre 30.0 pg landing is partly a coherence artifact of that shared "
            "input, not proof of 3 fully-orthogonal measurements; still a legitimate external-band check "
            "since the 27-33 pg reference range itself was never fitted to any of these cells.",
            "Route B's identity_check column is disclosed as algebraically non-independent (uses the cells' own "
            "own production rate); only the MUCKENTHALER-external-200B route_b_flux_mg_day is the real second leg.",
            "Absorption is held at its literature midpoint (1.5 mg/day) throughout, not itself swept -- a "
            "scope choice matching iron_hepcidin.py's convention, disclosed not hidden.",
        ],
    }
    report["step7_benign_cause_ruleouts"]["e_recycled_vs_new_dominance_ratio"] = round(dominance_ratio, 2)
    return report


def selftest():
    r = build_report()
    assert r["step0_shared_constant_audit"]["agree"] is True
    assert r["step0_shared_constant_audit"]["hufner_k_blood_oxygen_transport"] == 1.34
    assert abs(r["step1_independent_stoichiometry"]["replication_pctdiff"]) < 0.1
    assert r["step3_forced_adversary_grid_lifespan_x_bloodvolume"]["primary_cell_combined115_nadler5L"]["flux_mg_day"] == 22.603
    assert r["step3_forced_adversary_grid_lifespan_x_bloodvolume"]["grid_in_lit_band_fraction"] >= 0.8
    assert r["step3_forced_adversary_grid_lifespan_x_bloodvolume"]["grid_within_pre_registered_tolerance_fraction"] == 1.0
    assert r["step4_route_b_independent_bottom_up"]["route_b_in_lit_band"] is True
    assert r["step4_route_b_independent_bottom_up"]["route_a_vs_b_reldiff_pct"] < 15.0
    assert r["step6_mch_external_anchor"]["in_band"] is True
    assert r["step6_mch_external_anchor"]["mch_pg"] == 30.0
    assert r["step7_benign_cause_ruleouts"]["b_plasma_vs_blood_volume_counterfactual"]["correctly_rejected"] is True
    assert r["step8_recycled_absorbed_dominance"]["dominance_ratio"] >= 10.0
    assert r["overall_pass"] is True
    print("SELFTEST OK")
    return 0


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    report = build_report()
    print(json.dumps(report, indent=2, default=str))
    out_path = _os.path.join(OUT_DIR, "iron_oxygen_carriage_conservation_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    print(f"\nOVERALL: {'PASS' if report['overall_pass'] else 'FAIL/SURPRISE -- see gates above'}")
    return 0 if report["overall_pass"] else 2


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    sys.exit(main())
