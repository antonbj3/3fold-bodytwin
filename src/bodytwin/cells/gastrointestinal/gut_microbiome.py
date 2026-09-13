"""Gut microbiome: quantitative CENSUS + FERMENTATION-OUTPUT layer.

Bacterial density / bacteria:host-cell ratio, SCFA production/ratio/absorption, and butyrate as
colonocyte fuel -- explicitly NOT the unfalsifiable "dysbiosis" framing that dominates the popular
microbiome literature. Questions answered by computation: does the directly-measured colonic
bacterial census reproduce the corrected ~1:1 bacteria:human-cell ratio (refuting the entrenched
10:1 myth), does directly-measured SCFA production/ratio/absorption arithmetic self-consistently
close, and is butyrate really the dominant colonocyte fuel? Population-parametrized arithmetic
model; no individual data.

FALSIFIER (pre-registered, stated before any number below is computed):
  F1 (PRIMARY, census) the directly-reported bacteria:human-cell ratio (Sender/Fuchs/Milo 2016,
     PMID 27541692) lands in [0.3, 3.3] (same order of magnitude as 1:1) AND explicitly EXCLUDES
     [7, 13] (the entrenched 10:1 band) -- a conjunctive refutation, not just a point estimate.
  F2 (arithmetic cross-check) colonic bacterial density (1e11/mL, Sender Table 1) x colon content
     volume (0.4 L, Sender Box 1) reproduces Sender's independently-stated headline total
     (3.8e13) within 20% -- an internal geometric consistency check, not asserted.
  F3 (PRIMARY, SCFA) the directly-reported molar ratio acetate:propionate:butyrate (den Besten 2013,
     PMID 23821742, citing Cummings 1987 PMID 3678950 among its sources) is 60:20:20, and this
     falls inside the independently-reported cross-species/cross-diet band of Bergman 1990 (PMID
     2181501, 75:15:10 to 40:40:20) component-wise.
  F4 (PRIMARY, SCFA) directly-measured luminal SCFA concentration (Cummings 1987 autopsy intestinal-
     content sampling: 131 mmol/kg caecum, 80 mmol/kg descending colon; den Besten's independent
     70-140 mM proximal / 20-70 mM distal band) brackets the commonly-quoted ~100 mM point estimate.
  F5 (arithmetic cross-check) total SCFA production (400-600 mmol/day, den Besten) minus fecal
     excretion (5-30 mmol/day, diet-dependent, den Besten) reproduces the independently-STATED ~95%
     colonic absorption fraction (den Besten) to within 5 percentage points -- 3 numbers that were
     each separately reported must be mutually arithmetically consistent, not just individually true.
  F6 (decorrelated corroboration) butyrate accounts for ~70% of colonocyte oxidative energy
     (Roediger 1980, human colonocytes, PMID 7429343: 73%/75% ascending/descending, pure-substrate)
     and this DOMINANCE SURVIVES a forced adversary (competing 10 mM glucose present: 59%/72%, still
     >=50%) -- not an artifact of a substrate-starved assay.
  F7 (decorrelated corroboration, germ-free-vs-colonized) Backhed 2004 (PMID 15505215): direct
     germ-free-vs-conventionalized manipulation in mice shows a 60% body-fat increase in 14 days
     DESPITE REDUCED food intake -- a non-trivial, causally-forced caloric-harvest effect,
     independent methodologically from Bergman 1990's ~10%-of-human-calories mass-balance estimate.

SYMMETRIC QC, HELD OPEN PER TASK INSTRUCTION (not resolved away):
  - Microbiome composition (not just density) varies ENORMOUSLY between individuals and is diet-
    driven (Sender's Table 2: 24-78% CV across 12 independent stool-density studies; HMP 2012,
    PMID 22699609: community STRUCTURE varies widely across healthy subjects even though metabolic-
    PATHWAY carriage is comparatively stable; den Besten: diet shifts fecal SCFA output 2-3x).
  - Most human disease associations are CORRELATIONAL, not causal, and often drug/diet-confounded:
    Forslund et al 2015 (PMID 26633628, n=784 metagenomes) shows a T2D-microbiome "signature" was
    substantially a METFORMIN effect, not a disease effect, until stratified -- the field's
    clean confound-collapse case, freshly verified (not re-used from a sibling cell).
  - "Psychobiotic" / gut-brain causal-mood claims and most 16S disease-association literature are
    NOT re-litigated here -- out of scope for this census/fermentation layer.
  - SCFA numbers are frequently fecal-PROXY, not luminal/production, in the wider field (already
    flagged for Schwiertz 2010). This cell's
    headline numbers are NOT fecal-proxy (Cummings 1987 is direct autopsy intestinal-content
    sampling; den Besten/Bergman are production-flux syntheses) -- disclosed explicitly, Section D.

GEOMETRIC STRUCTURE (derive, don't assert):
  - The bacteria:human-cell ratio is a RATIO OF TWO INDEPENDENT VOLUME-INTEGRALS (bacteria: density
    x colon content volume; human: dominated by hematopoietic-lineage cell density x blood/marrow
    volume) -- the historical 10:1 myth is a GEOMETRIC (volumetric) error, not a density error: this
    script reconstructs Luckey 1972's back-of-envelope arithmetic and shows it used the SAME
    bacterial density this paper uses (1e11/g) but the WRONG integration volume (1 L whole
    alimentary tract vs the correct 0.4 L colon-restricted content) -- upstream-GI bacterial density
    is 3-8 orders of magnitude lower than colonic density (Sender Table 1: stomach/duodenum/jejunum
    1e3-1e4/mL, ileum 1e8/mL, colon 1e11/mL) because gastric acid + bile + fast transit suppress
    it -- a real anatomical fact, not an arbitrary refit.
  - SCFA absorption is a MASS-BALANCE CLOSURE: production - fecal_excretion = absorbed, and the
    resulting fraction must equal the independently-reported ~95% -- an over-determined system (3
    separately-measured quantities constrained by 1 conservation law), not a free parameter.

TWO CORRECTIONS MADE DURING THE BUILD (disclosed, not hidden):
  Fix #1: the common framing bundles "~10% of host energy from SCFA" together with "Turnbaugh
  2006 germ-free-transplant adiposity" as if one finding. Turnbaugh et al 2006 (PMID 17183312)
  reports NO percentage-of-host-energy figure at all -- it is a donor-COMPOSITION-dependent
  (ob/ob vs lean microbiota) adiposity-transfer result. The ~10% figure is Bergman 1990's
  (PMID 2181501) separate, independently-methodologized (species-comparative mass-balance)
  estimate. FIX: disentangled below (Section E) -- reported as two separate citations, not merged
  into a false single-source number.
  Fix #2 (a caught near-double-count): den Besten 2013's "~10% of human caloric requirements"
  sentence was initially going to be counted as a THIRD independent corroborating source alongside
  Bergman 1990 and Backhed 2004. Tracing den Besten's reference list showed its citation for
  that number IS Bergman 1990 (den Besten's ref 31) -- i.e. den Besten REPEATS Bergman's estimate,
  it does not independently re-derive it. FIX: this cell counts only 2 genuinely decorrelated
  angles for the caloric-harvest triangulation (Backhed's direct causal germ-free-vs-colonized
  manipulation; Bergman's independent mass-balance/stoichiometric estimate), not 3 -- avoiding a
  shared-source evidence-inflation false positive.

Citations: every PMID/DOI below was verified via NCBI eutils (esearch+efetch), several via full
PMC HTML fetch (not abstract-only) to pull exact in-text numbers. Backhed 2004 is reported
elsewhere as "+57% total body fat, +61% epididymal fat"; the abstract headline is a single "60%
increase in body fat content" figure -- the 57/61 sub-breakdown is not verified here and is NOT
repeated as if independently confirmed.

Reads: nothing.
Writes: gut_microbiome.json
Gate: F1-F7 (primary falsifiers F1/F3/F4/F6 plus the arithmetic closures) -> verdict field.
"""
import json
from pathlib import Path

import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_JSON = Path(OUT_ROOT) / "gut_microbiome" / "gut_microbiome.json"


def _to_native(obj):
    """Recursively convert numpy scalars/bools/arrays to native Python types so gates persist
    as REAL JSON booleans (machine-checkable), not stringified via a default=str fallback."""
    if isinstance(obj, dict):
        return {k: _to_native(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_native(v) for v in obj]
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return _to_native(obj.tolist())
    return obj

CITATIONS = {
    "sender_2016_bacteria_human_census": {
        "pmid": "27541692",
        "doi": "10.1371/journal.pbio.1002533",
        "cite": "Sender R, Fuchs S, Milo R (2016). Revised Estimates for the Number of Human and "
                "Bacteria Cells in the Body. PLoS Biol 14(8):e1002533.",
        "role": "PRIMARY census anchor. Verbatim (abstract): 'total number of bacteria in the 70 "
                "kg reference man to be 3.8x10^13 ... revise past estimates to 3.0x10^13 human "
                "cells ... updates the widely-cited 10:1 ratio, showing that the number of "
                "bacteria in the body is actually of the same order as the number of human cells, "
                "and their total mass is about 0.2 kg.' Full text (PMC4991899): B/H = 1.3, "
                "uncertainty 25%, population CV 53%. Colon density Table 1: 1e11 bacteria/mL "
                "content; colon content volume 0.4 L (Box 1). Table 2 compiles 12 independent "
                "empirical stool bacterial-density studies (1966-2008), range ~0.35-3.2x10^11/g "
                "wet stool, CV 24-78% across studies. Upper-GI density (Table 1): stomach/"
                "duodenum/jejunum 1e3-1e4/mL, ileum 1e8/mL -- 3-8 orders of magnitude below colon."
    },
    "luckey_1972_origin_of_myth": {
        "pmid": "4639749",
        "doi": "10.1093/ajcn/25.12.1292",
        "cite": "Luckey TD (1972). Introduction to intestinal microecology. Am J Clin Nutr "
                "25(12):1292-4.",
        "role": "Sender 2016's ref [3], explicitly named as the traceable 'single back-of-"
                "the-envelope estimate' origin of the ~1e14-bacteria (hence 10:1) figure: 'assumes "
                "the volume of the alimentary tract to be 1 liter ... multiplies this volume by "
                "the number density of bacteria, known to be about 1e11 bacteria per gram of wet "
                "content' -- SAME density this paper uses, WRONG (whole-tract, not colon-only) "
                "volume."
    },
    "savage_1977_propagation": {
        "pmid": "334036",
        "doi": "10.1146/annurev.mi.31.100177.000543",
        "cite": "Savage DC (1977). Microbial ecology of the gastrointestinal tract. Annu Rev "
                "Microbiol 31:107-33.",
        "role": "Sender 2016's ref [2], cited as a/the source of the 'previous estimate of "
                "about 1e14 bacteria' that propagated the pre-2016 consensus figure."
    },
    "cummings_1987_luminal_scfa": {
        "pmid": "3678950",
        "doi": "10.1136/gut.28.10.1221",
        "cite": "Cummings JH, Pomare EW, Branch WJ, Naylor CP, Macfarlane GT (1987). Short chain "
                "fatty acids in human large intestine, portal, hepatic and venous blood. Gut "
                "28(10):1221-7.",
        "role": "PRIMARY DIRECT luminal (NOT fecal) SCFA anchor -- autopsy sampling of sudden-"
                "death victims within 4h, all large-intestine regions. Verbatim: total SCFA "
                "(mmol/kg) 13+/-6 terminal ileum, 131+/-9 caecum, 80+/-11 descending colon; blood "
                "SCFA (umol/L) portal 375+/-70, hepatic 148+/-42, peripheral 79+/-22; 'acetate was "
                "the principal anion'; molar ratios shift colon->portal->hepatic (greater butyrate "
                "uptake by colonic epithelium, propionate by liver). One of den Besten 2013's "
                "3 cited sources (refs 9-11) for the 60:20:20 molar-ratio figure (confirmed by "
                "live full-text reference-list fetch)."
    },
    "den_besten_2013_scfa_review": {
        "pmid": "23821742",
        "doi": "10.1194/jlr.R036012",
        "cite": "den Besten G, van Eunen K, Groen AK, Venema K, Reijngoud DJ, Bakker BM (2013). The "
                "role of short-chain fatty acids in the interplay between diet, gut microbiota, "
                "and host energy metabolism. J Lipid Res 54(9):2325-40.",
        "role": "PRIMARY SCFA-quantity review, full text (PMC3735932) fetched live, not abstract-"
                "only. Verbatim: 'Acetate, propionate, and butyrate are present in an approximate "
                "molar ratio of 60:20:20 in the colon and stool (9-11)'; 'the total concentration "
                "of SCFAs decreases from 70 to 140 mM in the proximal colon to 20 to 70 mM in the "
                "distal colon (12)'; 'In the cecum and large intestine, 95% of the produced SCFAs "
                "are rapidly absorbed by the colonocytes while the remaining 5% are secreted in "
                "the feces (12-15)'; 'Fermentation of the carbohydrates reaching the cecum yield "
                "400-600 mmol SCFAs/day ... equivalent to ~10% of the human caloric requirements "
                "(31)'; 'Fecal secretion rates for SCFA are in the range of 10-30 mmol/day for "
                "diets with high fiber content compared with 5-15 mmol/day for control diets.' "
                "Live reference-list fetch confirms ref [31] = Bergman 1990 (below) -- this "
                "review's ~10% figure is NOT an independent re-derivation, it repeats Bergman's."
    },
    "bergman_1990_vfa_energy_review": {
        "pmid": "2181501",
        "doi": "10.1152/physrev.1990.70.2.567",
        "cite": "Bergman EN (1990). Energy contributions of volatile fatty acids from the "
                "gastrointestinal tract in various species. Physiol Rev 70(2):567-90.",
        "role": "PRIMARY independent caloric-contribution + cross-species SCFA-ratio anchor. "
                "Verbatim: 'principal VFA in either the rumen or large intestine are acetate, "
                "propionate, and butyrate and are produced in a ratio varying from approximately "
                "75:15:10 to 40:40:20'; 'VFA contribute approximately 70% to the caloric "
                "requirements of ruminants ... approximately 10% for humans, and approximately "
                "20-30% for several other omnivorous or herbivorous animals.'"
    },
    "roediger_1980_human_colonocytes": {
        "pmid": "7429343",
        "doi": "10.1136/gut.21.9.793",
        "cite": "Roediger WE (1980). Role of anaerobic bacteria in the metabolic welfare of the "
                "colonic mucosa in man. Gut 21(9):793-8.",
        "role": "PRIMARY HUMAN colonocyte-fuel anchor (isolated human colonocytes, ascending n=7 + "
                "descending n=7 surgical specimens). Verbatim: butyrate alone accounted for 73% "
                "(ascending) / 75% (descending) of O2 consumption; WITH 10 mM glucose present: 59% "
                "/ 72% -- FORCED ADVERSARY (does dominance survive substrate competition?): yes, "
                "attenuated not abolished. Glucose alone: 85% (ascending) / 30% (descending) O2 "
                "consumption; with 10mM butyrate present, glucose's share falls to 41%/16%."
    },
    "roediger_1982_rat_crossspecies": {
        "pmid": "7084619",
        "doi": None,
        "cite": "Roediger WE (1982). Utilization of nutrients by isolated epithelial cells of the "
                "rat colon. Gastroenterology 83(2):424-9.",
        "role": "Cross-species (RAT) replication, different assay batch/species than the 1980 "
                "human paper. Verbatim: butyrate (10 mM) alone accounted for 86% of total O2 "
                "consumption, suppressed endogenous-fuel oxidation by 82%; substrate preference "
                "order butyrate > acetate > propionate; glucose alone = 30% of O2 consumption. No "
                "DOI on record for this pre-DOI-era article (PMID sufficient, confirmed via "
                "esummary)."
    },
    "clausen_1994_rat_kinetics": {
        "pmid": "8299908",
        "doi": "10.1016/0016-5085(94)90601-7",
        "cite": "Clausen MR, Mortensen PB (1994). Kinetic studies on the metabolism of short-chain "
                "fatty acids and glucose by isolated rat colonocytes. Gastroenterology "
                "106(2):423-32.",
        "role": "DECORRELATED METHOD (Michaelis-Menten kinetics, not %-of-O2-consumption): Km "
                "order butyrate (0.184 mM) < propionate (0.339) < acetate (0.487) < glucose "
                "(0.777) -- butyrate has the HIGHEST affinity (lowest Km), supporting preferential "
                "use at physiological sub-saturating concentrations, independent of the O2-"
                "consumption-fraction metric above."
    },
    "clausen_1995_human_kinetics": {
        "pmid": "8549946",
        "doi": "10.1136/gut.37.5.684",
        "cite": "Clausen MR, Mortensen PB (1995). Kinetic studies on colonocyte metabolism of "
                "short chain fatty acids and glucose in ulcerative colitis. Gut 37(5):684-9.",
        "role": "HUMAN replication of the kinetic (Km) preference-ordering result: 'considerably "
                "lower Km value of butyrate ... supports a specific role of butyrate as an energy "
                "source for the colonic mucosa in both health and ulcerative colitis' (n=14 UC + "
                "n=8 control, no significant Vmax/Km difference by disease status)."
    },
    "donohoe_2011_mouse_mechanism": {
        "pmid": "21531334",
        "doi": "10.1016/j.cmet.2011.02.018",
        "cite": "Donohoe DR, Garge N, Zhang X, Sun W, O'Connell TM, Bunger MK, Bultman SJ (2011). "
                "The microbiome and butyrate regulate energy metabolism and autophagy in the "
                "mammalian colon. Cell Metab 13(5):517-26.",
        "role": "DECORRELATED METHOD + SPECIES (mouse; transcriptomic/metabolic-state readout, NOT "
                "O2-consumption tracer or kinetics). Verbatim: germ-free-mouse colonocytes are in "
                "an energy-DEPRIVED state (decreased NADH/NAD+, oxidative phosphorylation, ATP -> "
                "AMPK activation -> autophagy); adding butyrate rescues mitochondrial respiration "
                "and prevents autophagy, acting 'as an energy source rather than as an HDAC "
                "inhibitor' -- confirms NECESSITY of butyrate as fuel via a completely different "
                "experimental readout than Roediger's %-O2-consumption or Clausen's Km."
    },
    "backhed_2004_germfree_manipulation": {
        "pmid": "15505215",
        "doi": "10.1073/pnas.0407076101",
        "cite": "Backhed F, Ding H, Wang T, Hooper LV, Koh GY, Nagy A, Semenkovich CF, Gordon JI "
                "(2004). The gut microbiota as an environmental factor that regulates fat storage. "
                "Proc Natl Acad Sci U S A 101(44):15718-23.",
        "role": "PRIMARY causal germ-free-vs-colonized DECORRELATED anchor (direct binary "
                "presence/absence manipulation, not donor-composition-dependent). Verbatim: "
                "'conventionalization of adult germ-free (GF) C57BL/6 mice with a normal "
                "microbiota harvested from the distal intestine (cecum) of conventionally raised "
                "animals produces a 60% increase in body fat content and insulin resistance within "
                "14 days DESPITE REDUCED FOOD INTAKE' -- rules out an 'ate more' confound; "
                "mechanism = increased monosaccharide absorption + hepatic de novo lipogenesis, "
                "Fiaf-suppression-dependent."
    },
    "turnbaugh_2006_donor_composition": {
        "pmid": "17183312",
        "doi": "10.1038/nature05414",
        "cite": "Turnbaugh PJ, Ley RE, Mahowald MA, Magrini V, Mardis ER, Gordon JI (2006). An "
                "obesity-associated gut microbiome with increased capacity for energy harvest. "
                "Nature 444(7122):1027-31.",
        "role": "The standard citation for 'germ-free-transplant adiposity' -- re-read "
                "shows this paper's abstract reports NO percentage-of-host-energy figure; it "
                "is a DONOR-COMPOSITION-dependent (not germ-free-vs-colonized presence/absence) "
                "result: 'colonization of germ-free mice with an obese microbiota results in a "
                "significantly greater increase in total body fat than colonization with a lean "
                "microbiota' -- qualitative directional corroboration of Backhed's quantitative "
                "result, NOT the source of the '~10% of host energy' figure (that is "
                "Bergman 1990, disentangled explicitly in Section E, an OODA catch)."
    },
    "forslund_2015_metformin_confound": {
        "pmid": "26633628",
        "doi": "10.1038/nature15766",
        "cite": "Forslund K et al (2015). Disentangling type 2 diabetes and metformin treatment "
                "signatures in the human gut microbiota. Nature 528(7581):262-6.",
        "role": "Independently-verified confound anchor for "
                "symmetric QC: n=784 human gut metagenomes. Verbatim: 'antidiabetic medication "
                "confounds' the T2D-microbiome association; 'controlling for metformin treatment, "
                "we report a unified signature of gut microbiome shifts in T2D with a depletion of "
                "butyrate-producing taxa' -- disease association was substantially a DRUG effect "
                "until stratified. Directly ties back to this doc's butyrate theme."
    },
    "hmp_2012_interindividual_variation": {
        "pmid": "22699609",
        "doi": "10.1038/nature11234",
        "cite": "Human Microbiome Project Consortium (2012). Structure, function and diversity of "
                "the healthy human microbiome. Nature 486(7402):207-14.",
        "role": "Symmetric-QC anchor for inter-individual composition variation. Verbatim: 'healthy "
                "individuals differ remarkably in the microbes that occupy habitats ... diversity "
                "and abundance of each habitat's signature microbes [vary] widely even among "
                "healthy subjects'; BUT 'metagenomic carriage of metabolic pathways was stable "
                "among individuals despite variation in community structure' -- composition is "
                "individual, function is conserved (consistent with this repo's already-built "
                "MOL-MICROBIOME-GUT-BRAIN-HUB framing: 'taxonomy is the confound, function is the "
                "survivor')."
    },
}

TASK_BANDS = {
    "ratio_same_order_band": (0.3, 3.3),
    "ratio_old_myth_band": (7.0, 13.0),
    "arithmetic_tol_frac": 0.20,
    "scfa_ratio_target": (60.0, 20.0, 20.0),
    "bergman_ratio_endpoints": ((75.0, 15.0, 10.0), (40.0, 40.0, 20.0)),
    "task_production_floor_mmol_day": 50.0,
    "task_luminal_point_mM": 100.0,
    "absorption_fraction_stated": 0.95,
    "absorption_tol_pct_points": 5.0,
    "butyrate_task_target_pct": 70.0,
    "butyrate_tol_pct_points": 10.0,
    "butyrate_competition_floor_pct": 50.0,
    "caloric_harvest_nontrivial_floor_pct": 10.0,
    "bergman_human_task_adjacent_band": (5.0, 20.0),
}


def census_section():
    density_per_mL = 1e11
    colon_volume_mL = 400.0
    total_bacteria_headline = 3.8e13
    total_human_cells = 3.0e13
    bh_ratio_stated = 1.3
    bh_uncertainty_frac = 0.25
    bh_population_cv = 0.53
    table2_density_range_e11 = (0.35, 3.2)
    table2_cv_range_pct = (24, 78)
    upper_gi_density_range = (1e3, 1e8)  # stomach/duodenum/jejunum..ileum, per mL

    # F2: arithmetic cross-check (density x volume reproduces the independently-stated headline)
    arithmetic_total = density_per_mL * colon_volume_mL
    arithmetic_rel_err = abs(arithmetic_total - total_bacteria_headline) / total_bacteria_headline

    # F1: ratio checks
    ratio_lo, ratio_hi = TASK_BANDS["ratio_same_order_band"]
    myth_lo, myth_hi = TASK_BANDS["ratio_old_myth_band"]
    gate_ratio_same_order = ratio_lo <= bh_ratio_stated <= ratio_hi
    gate_ratio_excludes_myth = not (myth_lo <= bh_ratio_stated <= myth_hi)
    gate_arithmetic = arithmetic_rel_err <= TASK_BANDS["arithmetic_tol_frac"]

    # Myth-mechanism reconstruction (Luckey 1972's arithmetic, same density, wrong volume)
    luckey_volume_mL = 1000.0  # 1 L whole alimentary tract (Luckey's stated assumption)
    luckey_reconstructed_total = density_per_mL * luckey_volume_mL
    orders_of_magnitude_suppression_upper_gi = (
        np.log10(density_per_mL / upper_gi_density_range[1]),  # smallest gap (ileum, 1e8)
        np.log10(density_per_mL / upper_gi_density_range[0]),  # largest gap (stomach, 1e3)
    )
    gate_density_cross_validated = (
        table2_density_range_e11[0] <= (density_per_mL / 1e11) <= table2_density_range_e11[1]
    )
    gate_upper_gi_suppressed = orders_of_magnitude_suppression_upper_gi[0] >= 3.0

    return {
        "citations": ["sender_2016_bacteria_human_census", "luckey_1972_origin_of_myth",
                      "savage_1977_propagation"],
        "inputs": {
            "colon_bacterial_density_per_mL": density_per_mL,
            "colon_content_volume_mL": colon_volume_mL,
            "total_bacteria_headline": total_bacteria_headline,
            "total_human_cells": total_human_cells,
            "bh_ratio_stated": bh_ratio_stated,
            "bh_uncertainty_frac": bh_uncertainty_frac,
            "bh_population_cv": bh_population_cv,
            "total_bacteria_mass_kg": 0.2,
            "bacteria_mass_pct_bodyweight": 0.3,
            "table2_independent_stool_density_studies_n": 12,
            "table2_density_range_x1e11_per_g_wet_stool": table2_density_range_e11,
            "table2_cv_range_pct_across_studies": table2_cv_range_pct,
        },
        "derived": {
            "arithmetic_total_from_density_x_volume": arithmetic_total,
            "arithmetic_rel_err_vs_headline": round(arithmetic_rel_err, 4),
            "luckey_1972_reconstructed_total_bacteria": luckey_reconstructed_total,
            "luckey_used_same_density_wrong_volume_mL": luckey_volume_mL,
            "orders_of_magnitude_upper_gi_density_below_colon": [
                round(x, 2) for x in orders_of_magnitude_suppression_upper_gi
            ],
        },
        "gates": {
            "ratio_same_order_as_1to1": gate_ratio_same_order,
            "ratio_excludes_10to1_myth_band": gate_ratio_excludes_myth,
            "density_x_volume_reproduces_headline_lt20pct": gate_arithmetic,
            "density_cross_validated_by_12_independent_studies": gate_density_cross_validated,
            "upper_gi_density_suppressed_ge3_orders": gate_upper_gi_suppressed,
        },
    }


def scfa_section():
    ratio_target = np.array(TASK_BANDS["scfa_ratio_target"])
    bergman_lo_end, bergman_hi_end = TASK_BANDS["bergman_ratio_endpoints"]
    production_band = (400.0, 600.0)
    production_central = float(np.mean(production_band))
    task_floor = TASK_BANDS["task_production_floor_mmol_day"]
    luminal_proximal_mM = (70.0, 140.0)
    luminal_distal_mM = (20.0, 70.0)
    cummings_direct_mmol_per_kg = {"terminal_ileum": 13.0, "caecum": 131.0,
                                    "descending_colon": 80.0}
    fecal_band = {"control_diet": (5.0, 15.0), "high_fiber_diet": (10.0, 30.0)}
    stated_absorption_fraction = TASK_BANDS["absorption_fraction_stated"]

    # F3: ratio sums to 100, and component-wise brackets Bergman's cross-species/cross-diet band
    ratio_sum_ok = abs(float(ratio_target.sum()) - 100.0) < 1e-9
    bracket_lo = np.minimum(bergman_lo_end, bergman_hi_end)
    bracket_hi = np.maximum(bergman_lo_end, bergman_hi_end)
    component_in_band = [(bracket_lo[i] <= ratio_target[i] <= bracket_hi[i]) for i in range(3)]
    butyrate_at_edge = abs(ratio_target[2] - bracket_hi[2]) < 1e-9

    # F5: production vs task band (disclose undershoot, don't hide behind the "+")
    task_undershoot_factor = production_central / TASK_BANDS.get("task_lenient_ceiling", 100.0) \
        if False else production_central / 100.0
    gate_production_floor = production_central >= task_floor

    # F4: luminal concentration bracketing
    task_point = TASK_BANDS["task_luminal_point_mM"]
    in_proximal = luminal_proximal_mM[0] <= task_point <= luminal_proximal_mM[1]
    in_distal = luminal_distal_mM[0] <= task_point <= luminal_distal_mM[1]
    in_union = min(luminal_distal_mM[0], luminal_proximal_mM[0]) <= task_point <= \
        max(luminal_distal_mM[1], luminal_proximal_mM[1])

    # F5 (absorption arithmetic closure): production - fecal = absorbed; check vs stated 95%
    implied_fraction = {}
    for diet, (flo, fhi) in fecal_band.items():
        fecal_mid = float(np.mean([flo, fhi]))
        implied_fraction[diet] = round(1.0 - fecal_mid / production_central, 4)
    max_abs_dev_pct_points = max(
        abs(v * 100.0 - stated_absorption_fraction * 100.0) for v in implied_fraction.values()
    )
    gate_absorption_consistent = max_abs_dev_pct_points <= TASK_BANDS["absorption_tol_pct_points"]

    return {
        "citations": ["cummings_1987_luminal_scfa", "den_besten_2013_scfa_review",
                      "bergman_1990_vfa_energy_review"],
        "inputs": {
            "molar_ratio_acetate_propionate_butyrate": ratio_target.tolist(),
            "bergman_cross_species_diet_band_endpoints": [bergman_lo_end, bergman_hi_end],
            "total_production_band_mmol_day": production_band,
            "task_band_floor_mmol_day": task_floor,
            "luminal_proximal_colon_mM": luminal_proximal_mM,
            "luminal_distal_colon_mM": luminal_distal_mM,
            "cummings_1987_direct_autopsy_mmol_per_kg": cummings_direct_mmol_per_kg,
            "fecal_excretion_band_mmol_day": fecal_band,
            "stated_absorption_fraction": stated_absorption_fraction,
        },
        "derived": {
            "production_central_mmol_day": production_central,
            "task_band_undershoot_factor_vs_measured": round(task_undershoot_factor, 2),
            "component_in_bergman_band": component_in_band,
            "butyrate_component_at_bergman_band_edge": bool(butyrate_at_edge),
            "implied_absorption_fraction_by_diet": implied_fraction,
            "max_abs_dev_pct_points_vs_stated_95pct": round(max_abs_dev_pct_points, 2),
        },
        "gates": {
            "ratio_sums_to_100": ratio_sum_ok,
            "ratio_60_20_20_within_bergman_crossspecies_band": all(component_in_band),
            "production_central_ge_task_floor": gate_production_floor,
            "task_100mM_within_luminal_union_range": in_union,
            "task_100mM_within_proximal_range": in_proximal,
            "task_100mM_within_distal_range_ONLY": in_distal,
            "absorption_arithmetic_closes_within_5pts": gate_absorption_consistent,
        },
        "disclosed_mismatch": (
            f"The commonly-quoted '~50-100+ mmol/day' band undershoots the directly-verified literature "
            f"central estimate (den Besten 2013, citing Bergman 1990) of "
            f"{production_central:.0f} mmol/day by a factor of "
            f"~{task_undershoot_factor:.1f}x at its upper edge (100). The lenient "
            f"'+' qualifier is technically not violated (500 >= 50), but the headline number "
            f"reported here is the measured {production_band[0]:.0f}-{production_band[1]:.0f} "
            f"mmol/day, not the looser band."
        ),
    }


def colonocyte_fuel_section():
    human_butyrate_only = {"ascending": 73.0, "descending": 75.0}
    human_with_glucose = {"ascending": 59.0, "descending": 72.0}
    rat_butyrate_only = 86.0
    rat_glucose_only = 30.0
    task_target = TASK_BANDS["butyrate_task_target_pct"]
    tol = TASK_BANDS["butyrate_tol_pct_points"]
    competition_floor = TASK_BANDS["butyrate_competition_floor_pct"]

    human_mean_butyrate_only = float(np.mean(list(human_butyrate_only.values())))
    human_mean_with_glucose = float(np.mean(list(human_with_glucose.values())))
    abs_diff_from_task = abs(human_mean_butyrate_only - task_target)
    relative_attenuation_pct = round(
        100.0 * (1.0 - human_mean_with_glucose / human_mean_butyrate_only), 1
    )

    gate_matches_task = abs_diff_from_task <= tol
    gate_adversary_survives = human_with_glucose["ascending"] >= competition_floor and \
        human_with_glucose["descending"] >= competition_floor

    return {
        "citations": ["roediger_1980_human_colonocytes", "roediger_1982_rat_crossspecies",
                      "clausen_1994_rat_kinetics", "clausen_1995_human_kinetics",
                      "donohoe_2011_mouse_mechanism"],
        "inputs": {
            "human_pct_O2_consumption_butyrate_only": human_butyrate_only,
            "human_pct_O2_consumption_with_10mM_glucose": human_with_glucose,
            "rat_pct_O2_consumption_butyrate_only": rat_butyrate_only,
            "rat_pct_O2_consumption_glucose_only": rat_glucose_only,
            "task_target_pct": task_target,
        },
        "derived": {
            "human_mean_butyrate_only_pct": human_mean_butyrate_only,
            "human_mean_with_glucose_competition_pct": human_mean_with_glucose,
            "abs_diff_from_task_target_pct_points": abs_diff_from_task,
            "relative_attenuation_under_glucose_competition_pct": relative_attenuation_pct,
            "rat_vs_human_ratio": round(rat_butyrate_only / human_mean_butyrate_only, 2),
            "kinetic_Km_ordering_rat_mM": {"butyrate": 0.184, "propionate": 0.339,
                                             "acetate": 0.487, "glucose": 0.777},
            "donohoe_mouse_orthogonal_confirmation": (
                "germ-free colonocytes energy-deprived (autophagy/AMPK state); butyrate "
                "specifically rescues -- decorrelated method (transcriptomic/metabolic-state, "
                "not O2-consumption-fraction or Km), qualitative not numeric corroboration."
            ),
        },
        "gates": {
            "human_mean_within_10pts_of_task_70pct": gate_matches_task,
            "dominance_survives_glucose_competition_forced_adversary": gate_adversary_survives,
        },
    }


def caloric_harvest_decorrelated_section():
    backhed_body_fat_increase_pct = 60.0
    bergman_human_pct = 10.0
    nontrivial_floor = TASK_BANDS["caloric_harvest_nontrivial_floor_pct"]
    band_lo, band_hi = TASK_BANDS["bergman_human_task_adjacent_band"]

    gate_backhed_nontrivial = backhed_body_fat_increase_pct >= nontrivial_floor
    gate_bergman_in_band = band_lo <= bergman_human_pct <= band_hi

    return {
        "citations": ["backhed_2004_germfree_manipulation", "turnbaugh_2006_donor_composition",
                      "bergman_1990_vfa_energy_review", "den_besten_2013_scfa_review"],
        "inputs": {
            "backhed_2004_body_fat_increase_pct_14days": backhed_body_fat_increase_pct,
            "backhed_2004_food_intake_change": "REDUCED (rules out an 'ate more' confound)",
            "bergman_1990_human_caloric_contribution_pct": bergman_human_pct,
            "turnbaugh_2006_finding": "obese-donor > lean-donor adiposity transfer, SAME diet "
                                       "(donor-composition-dependent, qualitative, no % figure)",
        },
        "gates": {
            "backhed_germfree_manipulation_nontrivial_ge10pct": gate_backhed_nontrivial,
            "bergman_estimate_within_task_adjacent_5_20pct_band": gate_bergman_in_band,
        },
        "ooda_catches_disclosed": [
            "Fix #1: the common framing bundled '~10% of host energy' with 'Turnbaugh 2006 "
            "germ-free-transplant adiposity' as if one finding. Turnbaugh 2006's abstract "
            "(re-read live) reports NO percentage-of-host-energy figure -- it is a donor-"
            "composition result. The ~10% figure is Bergman 1990's separate estimate. "
            "Disentangled, not merged.",
            "Fix #2: den Besten 2013's '~10% of human caloric requirements' sentence was "
            "nearly counted as a 3rd independent corroborating source. Live full-text reference-"
            "list fetch shows den Besten's citation for that number IS Bergman 1990 (ref 31) "
            "-- a repeated citation, not an independent re-derivation. This section counts only "
            "2 genuinely decorrelated angles (Backhed's causal manipulation; Bergman's "
            "independent mass-balance estimate), not 3 -- avoiding a shared-source evidence-"
            "inflation false positive.",
        ],
    }


def symmetric_qc_section():
    return {
        "citations": ["hmp_2012_interindividual_variation", "forslund_2015_metformin_confound",
                      "den_besten_2013_scfa_review", "sender_2016_bacteria_human_census"],
        "held_open_not_resolved": {
            "individual_variation": (
                "Sender 2016's B/H ratio has 53% population CV; its own Table 2 compiles 12 "
                "independent stool bacterial-density studies with 24-78% CV EACH. HMP 2012 "
                "(PMID 22699609): microbial community STRUCTURE 'vary[ies] widely even among "
                "healthy subjects' -- but metabolic-PATHWAY carriage is comparatively stable. "
                "Composition varies enormously; this document's census/production NUMBERS are "
                "population medians/typical-orders-of-magnitude, not individual predictions."
            ),
            "diet_driven": (
                "den Besten 2013: fecal SCFA secretion 10-30 mmol/day (high-fiber) vs 5-15 "
                "mmol/day (control diet) -- ~2-3x diet-dependent swing; luminal concentration "
                "band itself is diet-dependent (70-140 vs 20-70 mM by colon region AND diet)."
            ),
            "disease_association_mostly_correlational": (
                "Forslund et al 2015 (PMID 26633628, n=784, freshly verified): a "
                "T2D gut-microbiome 'signature' reported by two prior unstratified studies was "
                "substantially a METFORMIN drug effect, not disease -- the field's clean "
                "confound-collapse case. This document does NOT extend its own census/SCFA/fuel "
                "numbers into any disease-association or causal claim."
            ),
            "psychobiotic_and_causal_overreach_NOT_relitigated": (
                "Already covered, more skeptically and at greater depth, by this repo's "
                "(read-only, NOT re-verified) graph cells: SYS-GUT-MICROBIOME-HOST-"
                "AXIS (flags PMID 29197739's null psychobiotic-mood RCT meta-analysis as "
                "open/contested) and MOL-MICROBIOME-GUT-BRAIN-HUB (its own Edge C explicitly "
                "'LOW confidence ... explicitly marked hypothesis'). MOL-MICROBIOME-METAGENOMICS-"
                "HOST-COUPLING's IBD/CGM anchor is only MARGINAL (AUROC 0.5825, z=1.45). This "
                "document does not re-derive or contradict those; it is a narrower, more basic "
                "quantitative layer underneath all of them."
            ),
            "scfa_fecal_proxy_caveat": (
                "This repo's MICROBIOME-HOST-METABOLIC cell already flags 'fecal SCFA "
                "concentration != production flux' (re: Schwiertz 2010) as a field-wide hazard. "
                "THIS document's headline numbers are NOT fecal-proxy: Cummings 1987 is "
                "direct autopsy intestinal-content sampling (not voided stool), and den Besten/"
                "Bergman's production-flux + 95%-absorption figures are explicit flux estimates, "
                "not concentration stand-ins. Fecal SCFA (the ~5% non-absorbed remainder) is "
                "reported here SEPARATELY (Section B) and never substituted for production. "
                "The caveat remains valid as a general field-level warning even though this "
                "document's core numbers clear that bar."
            ),
        },
    }


def main():
    census = census_section()
    scfa = scfa_section()
    fuel = colonocyte_fuel_section()
    harvest = caloric_harvest_decorrelated_section()
    qc = symmetric_qc_section()

    all_gates = {}
    all_gates.update({"census." + k: v for k, v in census["gates"].items()})
    all_gates.update({"scfa." + k: v for k, v in scfa["gates"].items()})
    all_gates.update({"fuel." + k: v for k, v in fuel["gates"].items()})
    all_gates.update({"harvest." + k: v for k, v in harvest["gates"].items()})

    n_pass = sum(1 for v in all_gates.values() if bool(v))
    n_total = len(all_gates)

    # Primary falsifier (pre-registered conjunction): census ratio (F1,F2) AND
    # SCFA ratio+concentration (F3,F4). Production-floor (F5-lite) and absorption-closure (F5)
    # are reported but the production-vs-task-band MISMATCH is disclosed, not gated as a hard
    # fail, per the lenient '+' qualifier -- explicitly flagged, not hidden.
    primary_falsifier_gates = [
        census["gates"]["ratio_same_order_as_1to1"],
        census["gates"]["ratio_excludes_10to1_myth_band"],
        scfa["gates"]["ratio_60_20_20_within_bergman_crossspecies_band"],
        scfa["gates"]["task_100mM_within_luminal_union_range"],
    ]
    primary_falsifier_pass = all(primary_falsifier_gates)

    secondary_corroboration_gates = [
        census["gates"]["density_x_volume_reproduces_headline_lt20pct"],
        scfa["gates"]["absorption_arithmetic_closes_within_5pts"],
        fuel["gates"]["human_mean_within_10pts_of_task_70pct"],
        fuel["gates"]["dominance_survives_glucose_competition_forced_adversary"],
        harvest["gates"]["backhed_germfree_manipulation_nontrivial_ge10pct"],
        harvest["gates"]["bergman_estimate_within_task_adjacent_5_20pct_band"],
    ]
    secondary_pass = all(secondary_corroboration_gates)

    report = {
        "citations": CITATIONS,
        "census": census,
        "scfa": scfa,
        "colonocyte_fuel": fuel,
        "caloric_harvest_decorrelated": harvest,
        "symmetric_qc_held_open": qc,
        "gates_all": all_gates,
        "gates_summary": {"n_pass": n_pass, "n_total": n_total},
        "verdict": {
            "primary_falsifier_ratio_and_scfa_conjunction": "PASS" if primary_falsifier_pass
                                                              else "FAIL",
            "secondary_corroboration_fuel_and_harvest": "PASS" if secondary_pass else "PARTIAL",
            "overall": "PASS" if (primary_falsifier_pass and secondary_pass) else
                       ("PARTIAL" if primary_falsifier_pass else "FAIL"),
        },
        "couples_to_PROSE_ONLY_not_graph_edited": [
            "GI (fermentation substrate): the gi_absorption_transit cell's small-"
            "bowel-transit + (unmodeled here) colonic-transit layer sets the rate at which "
            "fermentable substrate (resistant starch, NSP) reaches the colon -- an input this "
            "script takes as a disclosed fixed given (Western-diet fiber intake 20-25 g/day, "
            "den Besten), not itself modeled.",
            "Immune (SCFA -> Treg): mechanistically causal in gnotobiotic mice (Furusawa 2013, "
            "Smith 2013 GPR43/FFAR2 -- already live-verified by the existing, read-only "
            "SYS-GUT-MICROBIOME-HOST-AXIS graph cell, NOT re-verified in; cross-"
            "referenced, not re-derived).",
            "Metabolic (caloric harvest): Section D (this script) -- Backhed 2004 causal "
            "germ-free manipulation (+60% body fat/14d) and Bergman 1990's independent ~10%-of-"
            "human-calories mass-balance estimate, 2 genuinely decorrelated angles (an initial "
            "3rd angle, den Besten's repeat-citation of Bergman, was caught and NOT double-"
            "counted -- see ooda_catches_disclosed).",
        ],
        "honest_gaps": [
            "No subject-specific data anywhere -- population-parametrized census/production-"
            "arithmetic, matching every other systemic layer's disclosed scope.",
            "The '10:1 myth' correction's NUMERATOR mechanism (bacteria: ~1e14 -> 3.8e13, a "
            "volume-misattribution error) is directly reconstructed and verified; "
            "the myth's DENOMINATOR side (what 'old' human-cell figure was implicitly used to "
            "produce exactly '10:1') is NOT independently re-derived from a single live-"
            "verified primary source -- Sender 2016's paper is cited directly "
            "for the corrected ratio and its own explicit statement that it 'updates the widely-"
            "cited 10:1 ratio,' which is treated as sufficient without reconstructing every "
            "historical step of the denominator's provenance.",
            "The frequently-quoted SCFA total-production band ('~50-100+ mmol/day') is a substantial "
            "(~5x at its upper edge) undershoot of the directly-verified literature figure "
            "(400-600 mmol/day) -- disclosed prominently in scfa_section, not silently absorbed "
            "by the lenient '+' qualifier.",
            "Colonocyte-fuel % figures (Roediger 1980/1982) come from ISOLATED-CELL suspensions "
            "ex vivo, not in-vivo flux measurement -- the forced glucose-competition adversary "
            "(Section C) partially addresses the 'substrate-starved-assay-inflates-butyrate' "
            "concern but does not fully substitute for an in-vivo tracer study (none was found "
            "live with a directly quotable colonocyte-specific in-vivo % figure).",
            "Bacteria:human-cell ratio and SCFA numbers are BOTH population-level syntheses "
            "(Sender's 'reference man,' den Besten/Bergman's review-level synthesis) -- "
            "individual variation is large and diet-dependent (Section E / symmetric_qc), "
            "reported not resolved.",
        ],
        "confidence_tier": (
            "in-vivo-anchored (Sender/Fuchs/Milo 2016 census, PMID 27541692, cross-validated "
            "against 12 independent empirical stool-density studies in their own Table 2; "
            "Cummings 1987 direct autopsy luminal SCFA measurement, PMID 3678950; Roediger 1980 "
            "human isolated-colonocyte O2-consumption, PMID 7429343; Backhed 2004 causal germ-"
            "free-vs-colonized mouse manipulation, PMID 15505215) -- population-parametrized "
            "arithmetic, NOT subject-specific, matching gi_absorption_transit.py / "
            "renal_filtration.py's stated confidence tier. Weaker for the production-total "
            "and colonocyte-fuel-% numbers specifically (review-synthesized / ex-vivo isolated-"
            "cell respectively), disclosed in honest_gaps."
        ),
    }

    report = _to_native(report)

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(report, f, indent=2)

    print(json.dumps(_to_native({
        "verdict": report["verdict"],
        "gates_summary": report["gates_summary"],
        "gates_all": all_gates,
    }), indent=2))


if __name__ == "__main__":
    main()
