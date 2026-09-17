#!/usr/bin/env python3
"""MSK build: CLOSE THE CROSS-SCALE HANDOFF between the signalling model
(bone_rankl_opg_lemaire2004.py, osteoclast "concentration" C in
arbitrary pM units), the tissue-turnover model (bmu_turnover_kinetics.py,
annual %/yr activation-frequency), and the classic whole-skeleton calcium-flux
anchor (~500 mg/day, previously asserted from memory, never re-sourced).

TASK (per handoff): (1) build osteoclast-activity -> resorbed-mineral-mass
conversion from an INDEPENDENTLY MEASURED per-osteoclast/per-BMU rate (sourced,
not back-solved from the Lemaire ODE -- that would be a fit, not a conversion);
(2) test whether the classic ~500 mg/day anchor's literature already
includes RAPID bone-fluid/osteocytic calcium exchange, a method-family
distinct from histomorphometric remodeling, before treating a mg/day gap as a
model defect; (3) RE-FETCH the ~500 anchor itself -- it was flagged as
asserted from memory; (4) name the method family behind every number.

TWO DECORRELATED LEGS for whole-skeleton resorption flux (never averaged):

  LEG A -- TOP-DOWN, annual-%/yr activation-frequency route. Reuses this
  repo's LIVE-VERIFIED annual turnover anchors: Clarke 2008 (PMID
  18988698) cortical ~2-3%/yr (the repo's live-verified figure, NOT the
  folk "10%/yr" bmu_turnover_kinetics.py already flagged as an over-estimate
  vs Clarke), Eriksen 2010 trabecular ~25%/yr, weighted by a standard
  cortical:trabecular mass split (80:20, Clarke 2008 / standard skeletal
  anatomy), applied to ICRP Publication 70 Reference Man total skeletal
  calcium (~1000 g).

  LEG B -- BOTTOM-UP, single-BMU erosion-geometry route. INDEPENDENT
  measured parameters: Jaworski & Lok 1972 (Calcif Tissue Int, canine rib
  Haversian remodeling, cited as applicable to human lamellar dynamics)
  longitudinal osteoclastic erosion rate 39+-14 um/day; Qiu et al. 2003
  (PMID 12740946, human rib) osteon/cutting-cone cross-sectional area
  0.02-0.07 mm^2; ~1,000,000 concurrently active BMUs in the adult skeleton
  (Eriksen/Parfitt consensus, e.g. Jilka 2003 review); a derived cortical-bone
  calcium density from standard mineral-fraction (~65% ash by wet weight,
  Ca ~39.8% of hydroxyapatite by mass, apparent density ~1.9 g/cm^3) --
  this last figure is a STANDARD textbook combination, flagged, not
  independently re-measured when this cell was written.

EXTERNAL ANCHOR (decorrelated from both legs, isotope tracer-kinetics, NOT
histomorphometry): Bauer et al. (long-term whole-body 47Ca/85Sr tracer study,
PMC1935462) report accretion to the "fixed" bone pool of 100-210 mg Ca/day in
normal adults, and state this is ~1/3 of the accretion to the TOTAL bone
pool -> implied total tracer-kinetic accretion ~300-630 mg/day. This is the
RE-FETCHED anchor (the ~500 mg/day figure could NOT be re-confirmed as a
specific primary-sourced number when this cell was written -- Guyton & Hall's chapter
text, fetched directly, does not state it; it is retired as an unverified
memory-asserted number, per handoff instruction 3).

METHOD-FAMILY CANDIDATE FOR ANY RESIDUAL GAP: rapid bone-fluid / osteocytic
("perilacunar") calcium exchange is qualitatively reported (Parfitt 1993,
Handbook of Exp. Pharmacology, cited via PMC4377030) as "many fold higher"
than remodeling-based flux. This is a PHYSIOLOGICALLY DISTINCT, faster
process (labile surface pool <-> ECF, minutes-hours) that histomorphometry
structurally CANNOT see (it only resolves net BMU-based bone volume change).
Isotope tracer-kinetics (Leg-anchor above) CAN partially pick up this pool
(rapid isotope redistribution into exchangeable bone Ca), which is the
candidate mechanism for tracer-kinetic estimates running higher than
histomorphometric ones. HONEST GAP: no numeric fold-factor for the rapid
exchange flux could be machine-extracted from primary text when this cell was written
(only the qualitative "many fold" claim, itself a secondary citation) --
this piece is NOT independently quantified, flagged, not fabricated.

THE LEMAIRE ODE's C VARIABLE: deliberately NOT bridged to mg/day here.
Inverting the ODE's abstract pM-unit C to a mass flux would require fitting
an extra calibration constant against these same target numbers -- circular
by the handoff's instruction ("source it rather than back-solving, or
the conversion is a fit"). The correct, non-circular position is: C remains
an uncalibrated RELATIVE signalling index; the mass conversion is built
bottom-up from Legs A/B instead, fully decorrelated from the ODE.
"""
import json
from pathlib import Path

from pathlib import Path as _Path
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

REPO = _Path(OUT_ROOT)
OUT = REPO / "bone_calcium_flux_reconciliation" / "bone_calcium_flux_reconciliation_results.json"

# ---------------------------------------------------------------- Leg A ----
ICRP_SKELETAL_CA_G = 1000.0          # ICRP Publication 70 Reference Man, standard
CORTICAL_FRACTION = 0.80             # standard adult cortical:trabecular mass split
TRABECULAR_FRACTION = 0.20
CORTICAL_PCT_YR_RANGE = (0.02, 0.03)   # Clarke 2008 PMID 18988698, live-verified in repo
TRABECULAR_PCT_YR = 0.25               # Eriksen 2010, live-verified in repo

def leg_a_mg_day():
    lo_cort, hi_cort = CORTICAL_PCT_YR_RANGE
    out = []
    for cort_pct in (lo_cort, hi_cort):
        annual_g = (CORTICAL_FRACTION * ICRP_SKELETAL_CA_G * cort_pct
                    + TRABECULAR_FRACTION * ICRP_SKELETAL_CA_G * TRABECULAR_PCT_YR)
        mg_day = annual_g * 1000.0 / 365.0
        out.append(mg_day)
    return min(out), max(out)

# ---------------------------------------------------------------- Leg B ----
EROSION_LONGITUDINAL_UM_DAY = (25.0, 53.0)   # Jaworski & Lok 1972, 39+-14 um/day, +-1SD band
OSTEON_AREA_MM2 = (0.02, 0.07)                # Qiu et al. 2003, PMID 12740946
ACTIVE_BMUS = 1.0e6                            # Eriksen/Parfitt consensus (Jilka 2003)
# derived cortical-bone Ca density: apparent density ~1.9 g/cm3 * mineral(ash)
# fraction ~0.65 * Ca fraction of hydroxyapatite (Ca10(PO4)6(OH)2, M=1004.6,
# Ca mass fraction = 40.08*10/1004.6 = 0.3989) -- standard textbook combination
APPARENT_DENSITY_G_CM3 = 1.9  # g/cm3
MINERAL_ASH_FRACTION = 0.65
CA_FRACTION_OF_HA = 40.08 * 10 / 1004.6
CA_DENSITY_MG_MM3 = (APPARENT_DENSITY_G_CM3 * MINERAL_ASH_FRACTION * CA_FRACTION_OF_HA  # mg/mm3
                      * 1000.0 / 1000.0)  # g/cm3 -> mg/mm3 (1 g/cm3 = 1 mg/mm3)

def leg_b_mg_day():
    out = []
    for area in OSTEON_AREA_MM2:
        for eros_um in EROSION_LONGITUDINAL_UM_DAY:
            vol_mm3_day_per_bmu = area * (eros_um / 1000.0)
            mg_day_per_bmu = vol_mm3_day_per_bmu * CA_DENSITY_MG_MM3
            out.append(mg_day_per_bmu * ACTIVE_BMUS)
    return min(out), max(out)

# ------------------------------------------------------------ external ----
# Bauer et al. (PMC1935462): accretion to "fixed" bone pool 100-210 mg/day,
# stated as ~1/3 of accretion to TOTAL bone pool -> implied total range.
BAUER_FIXED_POOL_MG_DAY = (100.0, 210.0)
FIXED_POOL_FRACTION_OF_TOTAL = 1.0 / 3.0

def tracer_kinetic_anchor_mg_day():
    lo, hi = BAUER_FIXED_POOL_MG_DAY
    return lo / FIXED_POOL_FRACTION_OF_TOTAL, hi / FIXED_POOL_FRACTION_OF_TOTAL

def overlap(a, b):
    lo = max(a[0], b[0]); hi = min(a[1], b[1])
    return (lo, hi) if lo <= hi else None

def gap_ratio(a, b):
    """max/min across the two ranges' nearest edges -- the residual multiplicative gap."""
    if overlap(a, b) is not None:
        return 1.0
    if a[1] < b[0]:
        return b[0] / a[1]
    return a[0] / b[1]

def main():
    a_lo, a_hi = leg_a_mg_day()
    b_lo, b_hi = leg_b_mg_day()
    ext_lo, ext_hi = tracer_kinetic_anchor_mg_day()

    classic_500_reverified = False  # explicit: could NOT re-confirm a primary
                                     # source stating exactly "~500 mg/day"
                                     # (Guyton & Hall chapter text fetched
                                     # directly does not contain it)

    results = {
        "leg_a_topdown_annual_pct_route_mg_day": [round(a_lo, 1), round(a_hi, 1)],
        "leg_a_sources": ["Clarke 2008 PMID 18988698 (cortical 2-3%/yr, live-verified in repo)",
                          "Eriksen 2010 (trabecular 25%/yr, live-verified in repo)",
                          "ICRP Publication 70 Reference Man skeletal Ca ~1000 g"],
        "leg_b_bottomup_bmu_geometry_mg_day": [round(b_lo, 1), round(b_hi, 1)],
        "leg_b_sources": ["Jaworski & Lok 1972 (Calcif Tissue Int) erosion rate 39+-14 um/day",
                          "Qiu et al. 2003 PMID 12740946 osteon cross-section 0.02-0.07 mm^2",
                          "~1e6 active BMUs, Eriksen/Parfitt consensus (Jilka 2003 review)",
                          "derived Ca density from standard apparent-density x ash-fraction x "
                          "hydroxyapatite-Ca-mass-fraction (textbook combination, flagged)"],
        "external_anchor_tracer_kinetic_mg_day": [round(ext_lo, 1), round(ext_hi, 1)],
        "external_anchor_source": ("Bauer et al. long-term whole-body 47Ca/85Sr tracer study "
                                    "(PMC1935462): fixed-bone-pool accretion 100-210 mg/day, "
                                    "stated ~1/3 of total-bone-pool accretion"),
        "classic_500mg_day_reverified_as_primary_sourced": classic_500_reverified,
        "classic_500_reverify_note": ("Guyton & Hall chapter text (doctorlib.org, fetched "
                                       "directly when this cell was written) does NOT contain a '500 mg/day' "
                                       "statement; no other primary source located when this cell was written "
                                       "either. RETIRED as an unverified memory-asserted anchor "
                                       "per handoff instruction 3 -- replaced by the tracer-"
                                       "kinetic external anchor above."),
        "leg_a_vs_external_gap_ratio": round(gap_ratio((a_lo, a_hi), (ext_lo, ext_hi)), 2),
        "leg_a_vs_external_overlap_mg_day": overlap((a_lo, a_hi), (ext_lo, ext_hi)),
        "leg_b_vs_external_overlap_mg_day": overlap((b_lo, b_hi), (ext_lo, ext_hi)),
        "rapid_bone_fluid_exchange": {
            "claim": ("rapid osteocytic/bone-fluid ('perilacunar') calcium exchange is a "
                      "physiologically distinct, faster (minutes-hours) process vs BMU-based "
                      "remodeling (weeks-months); qualitatively reported 'many fold higher' "
                      "flux than remodeling-based turnover"),
            "source": "Parfitt 1993 (Handbook of Experimental Pharmacology), cited via PMC4377030",
            "numeric_fold_factor_quantified_this_session": False,
            "honest_gap": ("no absolute mg/day or a specific numeric fold-factor for the rapid "
                            "exchange flux could be machine-extracted from primary text this "
                            "run -- only the qualitative 'many fold' claim, itself a "
                            "secondary citation. NOT fabricated as a number."),
            "mechanism_role": ("candidate explanation for why isotope tracer-kinetic estimates "
                               "(external anchor) run higher than histomorphometry-only "
                               "estimates (Leg A): tracer methods can pick up rapid isotope "
                               "exchange into the labile surface pool in addition to true net "
                               "remodeling accretion; histomorphometry structurally cannot see "
                               "this pool at all (it resolves BMU bone-volume change only). This "
                               "is a METHOD-FAMILY difference, not a numeric defect in either "
                               "model, IF the residual gap is of the right sign and rough "
                               "magnitude -- checked below."),
        },
        "lemaire_ode_osteoclast_C_bridge": {
            "bridged_to_mg_day": False,
            "reason": ("Lemaire2004's C is an abstract pM-unit signalling proxy for the active "
                       "osteoclast pool, never calibrated against histomorphometric osteoclast "
                       "density or resorbed volume in the original paper or the project. Inverting "
                       "it to mg/day would require fitting a new calibration constant against "
                       "these same target numbers -- circular per the handoff's instruction. "
                       "Legs A/B instead build the mass conversion bottom-up from independently "
                       "measured parameters, fully decorrelated from the ODE's C."),
        },
    }

    # verdict: does the gap SHRINK once the anchor is properly re-sourced, and
    # is the residual gap's SIGN consistent with the rapid-exchange-flux
    # mechanism (external/tracer-kinetic should run >= Leg A/histomorphometric,
    # since tracer methods can additionally capture some of the faster pool)?
    sign_consistent = ext_lo <= a_hi * 3.0 and ext_lo >= a_lo  # anchor not below leg A floor,
    # and not so far above as to need a >3x mechanism to close (would falsify
    # "many fold" as a full explanation, only a partial one)
    results["verdict"] = {
        "gap_vs_original_handoff_2x_claim": ("SHRUNK: prior unwired estimate was 205-266 vs an "
                                              "unverified ~500 (ratio ~1.9-2.4x). Leg A recomputed "
                                              f"fresh here ({round(a_lo,1)}-{round(a_hi,1)}) vs the "
                                              "re-sourced tracer-kinetic anchor "
                                              f"({round(ext_lo,1)}-{round(ext_hi,1)}) gives ratio "
                                              f"{round(gap_ratio((a_lo,a_hi),(ext_lo,ext_hi)),2)}x, "
                                              "and Leg B (wide, geometry-uncertain) OVERLAPS the "
                                              "anchor entirely."),
            "sign_consistent_with_rapid_exchange_mechanism": sign_consistent,
            "classification": ("PARTIAL METHOD-FAMILY EXPLANATION, not a full quantitative "
                                "closure: the direction (tracer-kinetic > histomorphometric) "
                                "matches the rapid-exchange-adds-to-isotope-uptake mechanism, "
                                "and Leg B's independent bottom-up geometry route is fully "
                                "consistent with the anchor. But the rapid-exchange flux itself "
                                "is NOT independently quantified when this cell was written, so the residual "
                                "~1.1-1.6x gap between Leg A and the anchor is NOT closed by a "
                                "number, only rendered PLAUSIBLE and same-sign. Honest partial "
                                "negative, not a full confirmation."),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
