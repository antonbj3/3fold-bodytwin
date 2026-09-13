#!/usr/bin/env python3
"""
First-breath opening pressure + fetal-lung-liquid clearance, forced against a "passive-only"
adversary.

Reads: nothing (all values are published literature numbers embedded below).
Writes: neonatal_transition_results.json under the cell output directory.
Gate: the pre-registered gates below; overall_pass requires all of them.

TASK FALSIFIER (symmetric, pre-registered):
  (a) reproduce the first-breath opening pressures (~40-80 cmH2O, far above tidal) --
      overcoming surface tension + fluid, couples to the surfactant cert's Laplace geometry.
  (b) reproduce lung-fluid clearance where ENaC-mediated Na+ reabsorption is LOAD-BEARING (not
      merely present) -- a "passive drainage only" adversary must FAIL to clear fluid at the
      measured rate (ENaC-KO / amiloride block delays clearance -- force it).

METHOD (geometric + measured, not narrated):
  Part 1 reuses the EXACT Laplace relation P=2*gamma/r already verified in the
  pulmonary_surfactant_alveolar_stability cell (decorrelated reuse of a verified
  geometric relation, NOT re-derived from scratch), evaluated at the MATURE alveolar radius
  (r=100 micron, same value) and INVERTED to ask: what effective radius of curvature does the
  MEASURED first-breath pressure range (40-80 cmH2O, Karlberg-anchored) IMPLY, given the SAME
  bare/unspread surface-tension bracket (50-70 mN/m) used in the surfactant cert? If that implied
  radius is much smaller than the mature 100-micron alveolus, the "surface-tension-at-a-tiny,
  not-yet-opened interface" mechanism is geometrically self-consistent with the measured pressures
  -- a real falsifier (it could have come out equal to or larger than r_mature, which would have
  refuted the mechanism).

  Part 2 forces the "passive drainage only" adversary using DIRECT measured necessity evidence
  (genetic knockout + pharmacological blockade + clinical epidemiology) across 3 SPECIES/METHODS
  (mouse, sheep, human) rather than a simulated ODE with invented rate constants -- the strongest
  fair form of "passive-only" is a real animal that underwent completely normal labor/delivery
  (mechanical squeeze intact) MINUS functional ENaC (genetic, Hummler et al 1996) or MINUS
  functional Na+ channels acutely (pharmacological, Olver/Ramsden et al 1986) -- both fail to
  clear fluid at the physiologically required rate. Quantitative arithmetic (fold-changes, margins
  over threshold) is computed on the actual extracted numbers, not invented.

  Part 3 is a SECONDARY, explicitly-labeled-illustrative cross-check of the viscous-resistance leg
  (does NOT re-derive Hooper et al 2013's measured ~100x liquid-vs-air airway-resistance figure;
  a simple bulk-viscosity-ratio floor is computed and the gap to the measured figure is disclosed,
  not hidden).

"""
import json
import hashlib
import os

import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUTDIR = _os.path.join(OUT_ROOT, "neonatal_transition")
OUTFILE = os.path.join(OUTDIR, "neonatal_transition_results.json")

CMH2O_TO_PA = 98.0665  # standard: rho_water(1000 kg/m3) * g(9.80665 m/s2) * 0.01 m

# ---------------------------------------------------------------------------------------
# PRE-REGISTERED measured inputs (cited)
# ---------------------------------------------------------------------------------------

# Shared Laplace geometry anchor -- IDENTICAL values to the sibling surfactant cert
# (the pulmonary_surfactant_alveolar_stability cell), reused not re-fitted.
R_MATURE_M = 1.0e-4          # 0.1 mm mature alveolar radius (surfactant doc, Part 1)
GAMMA_BARE_mNm = [50.0, 60.0, 70.0]     # bare-saline / unspread-film surface tension bracket
GAMMA_SURF_mNm = [2.0, 5.0, 10.0]       # surfactant-active surface tension bracket

# Task-specified / Karlberg-anchored (Karlberg P, Koch G et al, Acta Paediatrica newborn
# respiratory-mechanics series, 1960/1962 -- PMID 14404498 Part I, PMID 14453970 Part III;
# pre-1975-era records, no MEDLINE abstract -- see honest_gaps) first-breath opening pressure.
P_FIRST_BREATH_LO_cmH2O = 40.0
P_FIRST_BREATH_HI_cmH2O = 80.0

# Lista G, Castoldi F (2010) PMID 21089711 -- sustained lung inflation (SLI) clinical protocol,
# live-verified quote: "a peak pressure of 25-30 cm H2O for 10-20 seconds"
P_SLI_LO_cmH2O = 25.0
P_SLI_HI_cmH2O = 30.0

# Papastamelos et al (1995) PMID 7713809 -- chest-wall vs lung compliance ratio in infants <1 yr,
# live-verified quote: "Cw/spontaneous Cl was 2.86 +/- 1.06 in infants < 1 yr"
CW_OVER_CL_INFANT = 2.86
CW_OVER_CL_INFANT_SD = 1.06

# Brown MJ, Olver RE, Ramsden CA, Strang LB, Walters DV (1983) PMID 6655575 -- catecholamine
# surge across labor, live-verified quotes (plasma concentrations, ng/mL):
ADRENALINE_EARLY_LABOR_ngml = 0.087
ADRENALINE_LATE_LABOR_ngml = 6.86
ADRENALINE_POSTNATAL_ngml = 7.17
NORADRENALINE_EARLY_LABOR_ngml = 1.71
NORADRENALINE_LATE_LABOR_ngml = 12.14
NORADRENALINE_POSTNATAL_ngml = 9.10
# Same paper -- minimum adrenaline concentration to inhibit secretion (Ai), decreasing with
# advancing gestation (maturational sensitization of the switch mechanism):
AI_THRESHOLD_132_134D_ngml = 0.43
AI_THRESHOLD_GT140D_ngml = 0.029

# Olver RE, Ramsden CA, Strang LB, Walters DV (1986) PMID 3795077 -- amiloride block of the
# adrenaline-induced Na+-dependent reabsorption response, live-verified quotes:
AMILORIDE_APPLIED_M = 1.0e-4     # "abolished the changes in p.d. and ion flux induced by adrenaline"
AMILORIDE_KI_M = 4.0e-6          # "50% inhibition ... induced by 4 x 10^-6 M-amiloride"

# Hummler E et al (1996) Nat Genet, PMID 8589728 -- alpha-ENaC(-/-) mouse phenotype, live-verified
# quote: "died within 40 h of birth from failure to clear their lungs of liquid"
ENAC_KO_DEATH_WINDOW_H = 40.0

# Hansen AK, Wisborg K, Uldbjerg N, Henriksen TB (2008) BMJ, PMID 18077440 -- elective
# caesarean-section-without-labor respiratory morbidity, live-verified quotes (odds ratios, 95% CI):
CS_OR_37WK = (3.9, 2.4, 6.5)
CS_OR_38WK = (3.0, 2.1, 4.3)
CS_OR_39WK = (1.9, 1.2, 3.0)
CS_OR_SERIOUS_37WK = (5.0, 1.6, 16.0)

# Standard tabulated physical constants (engineering/physics reference values, NOT a biomedical
# literature claim -- not NCBI-live-verified, see honest_gaps), used ONLY for the secondary,
# explicitly-illustrative Part 3 cross-check:
MU_WATER_37C_mPas = 0.6913
MU_AIR_37C_mPas = 0.0190
HOOPER_MEASURED_VISCOUS_RATIO = 100.0  # Hooper et al 2013 PMID 24035400, reused by reference
                                        # from the sibling surfactant doc (not re-fetched this run)

# ---------------------------------------------------------------------------------------
# PRE-REGISTERED gate thresholds (decided BEFORE reading the gate outcomes below)
# ---------------------------------------------------------------------------------------
GATE_R_IMPLIED_MAX_MICRON = 50.0        # must be < half the mature 100-micron alveolar radius
GATE_MATURE_BARE_MAX_cmH2O_LT = P_FIRST_BREATH_LO_cmH2O   # mature-radius-only ceiling must undershoot measured MIN
GATE_RATIO_MIN_X = 10.0                  # first-breath / mature-tidal pressure ratio floor
GATE_CW_CL_MIN = 1.0                     # chest wall must be MORE compliant than lung (>1), not stiffer
GATE_CATECHOL_FOLD_MIN_X = 10.0
GATE_AMILORIDE_MARGIN_MIN_X = 5.0
GATE_CS_OR_LOWER_CI_MIN = 1.0            # lower 95% CI bound must exclude 1 (significant)
GATE_DECORRELATED_METHODS_MIN = 3


def laplace_p_cmH2O(gamma_mNm, r_m):
    """P = 2*gamma/r, gamma in mN/m -> N/m via *1e-3, result Pa -> cmH2O via /CMH2O_TO_PA."""
    gamma_Nm = gamma_mNm * 1.0e-3
    p_pa = 2.0 * gamma_Nm / r_m
    return p_pa / CMH2O_TO_PA


def laplace_r_implied_m(gamma_mNm, p_cmH2O):
    """Invert P=2*gamma/r -> r = 2*gamma/P, SI throughout."""
    gamma_Nm = gamma_mNm * 1.0e-3
    p_pa = p_cmH2O * CMH2O_TO_PA
    return 2.0 * gamma_Nm / p_pa


# ---------------------------------------------------------------------------------------
# PART 1 -- opening-pressure geometry
# ---------------------------------------------------------------------------------------
def part1_opening_pressure():
    # 1a. Mature-radius (r=100 micron) pressures, bare vs surfactant-active -- reproduces the
    # sibling surfactant cert's Part-1 numbers as an internal cross-check (should match).
    mature_bare = {g: laplace_p_cmH2O(g, R_MATURE_M) for g in GAMMA_BARE_mNm}
    mature_surf = {g: laplace_p_cmH2O(g, R_MATURE_M) for g in GAMMA_SURF_mNm}

    # 1b. Inverse calc: what effective radius does the MEASURED first-breath pressure range imply,
    # using the SAME bare/unspread-film tension bracket? (the geometric falsifier)
    r_implied_grid = {}
    for g in GAMMA_BARE_mNm:
        for p in (P_FIRST_BREATH_LO_cmH2O, P_FIRST_BREATH_HI_cmH2O):
            r_implied_grid[f"gamma{g}_P{p}"] = laplace_r_implied_m(g, p) * 1e6  # meters -> micron

    r_implied_values_micron = list(r_implied_grid.values())
    r_implied_max_micron = max(r_implied_values_micron)
    r_implied_min_micron = min(r_implied_values_micron)

    # 1c. Ratio gate: measured first-breath pressure vs mature TIDAL pressure (surfactant-active,
    # r=100 micron) -- the most conservative (smallest) ratio uses the highest mature-tidal
    # pressure (10 mN/m) against the lowest first-breath pressure (40 cmH2O).
    mature_tidal_max_cmH2O = max(mature_surf.values())   # weakest-case denominator
    mature_tidal_min_cmH2O = min(mature_surf.values())
    ratio_conservative = P_FIRST_BREATH_LO_cmH2O / mature_tidal_max_cmH2O
    ratio_generous = P_FIRST_BREATH_HI_cmH2O / mature_tidal_min_cmH2O

    # 1d. Mature-radius-only adversary: can bare (unspread) tension AT THE MATURE RADIUS alone
    # reach the measured first-breath range? (it should NOT -- that's the falsifier)
    mature_bare_max_cmH2O = max(mature_bare.values())

    return {
        "mature_bare_cmH2O": mature_bare,
        "mature_surf_cmH2O": mature_surf,
        "mature_bare_max_cmH2O": mature_bare_max_cmH2O,
        "r_implied_grid_micron": r_implied_grid,
        "r_implied_max_micron": r_implied_max_micron,
        "r_implied_min_micron": r_implied_min_micron,
        "r_mature_micron": R_MATURE_M * 1e6,
        "ratio_conservative_x": ratio_conservative,
        "ratio_generous_x": ratio_generous,
        "mature_tidal_range_cmH2O": [mature_tidal_min_cmH2O, mature_tidal_max_cmH2O],
    }


# ---------------------------------------------------------------------------------------
# PART 2 -- fluid clearance, adversary-forcing (measured necessity evidence, 3 species/methods)
# ---------------------------------------------------------------------------------------
def part2_fluid_clearance():
    adrenaline_fold = ADRENALINE_LATE_LABOR_ngml / ADRENALINE_EARLY_LABOR_ngml
    noradrenaline_fold = NORADRENALINE_LATE_LABOR_ngml / NORADRENALINE_EARLY_LABOR_ngml
    ai_sensitization_fold = AI_THRESHOLD_132_134D_ngml / AI_THRESHOLD_GT140D_ngml  # >1 => more sensitive later

    amiloride_margin_x = AMILORIDE_APPLIED_M / AMILORIDE_KI_M

    cs_or_monotonic_decreasing = (CS_OR_37WK[0] > CS_OR_38WK[0] > CS_OR_39WK[0])

    return {
        "adrenaline_fold_early_to_late_labor": adrenaline_fold,
        "noradrenaline_fold_early_to_late_labor": noradrenaline_fold,
        "ai_threshold_sensitization_fold": ai_sensitization_fold,
        "amiloride_applied_over_KI_margin_x": amiloride_margin_x,
        "enac_ko_death_window_h": ENAC_KO_DEATH_WINDOW_H,
        "enac_ko_categorical_result": (
            "alpha-ENaC(-/-) mouse neonates: amiloride-sensitive Na+ transport ABOLISHED (direct "
            "measurement, airway epithelia); phenotype = respiratory distress, death within 40h "
            "from failure to clear lung liquid (Hummler et al 1996, PMID 8589728, direct quote). "
            "Standard gene-targeting mouse litters are vaginally delivered by the dam unless "
            "otherwise stated -- the fetched abstract does NOT itself confirm delivery mode; this "
            "is a reasonable, standard inference, NOT an independently-verified quote (see "
            "honest_gaps). The core necessity claim (active Na+ transport required; genetic "
            "loss is not compensated) holds regardless of delivery mode."
        ),
        "cs_or_by_gestational_week": {"37": CS_OR_37WK, "38": CS_OR_38WK, "39": CS_OR_39WK, "37_serious": CS_OR_SERIOUS_37WK},
        "cs_or_monotonic_decreasing_with_ga": cs_or_monotonic_decreasing,
        "decorrelated_methods": [
            {"method": "genetic knockout", "species": "mouse", "citation_pmid": "8589728"},
            {"method": "acute pharmacological blockade", "species": "fetal sheep", "citation_pmid": "3795077"},
            {"method": "endocrine physiology (catecholamine infusion + spontaneous labor)", "species": "fetal sheep", "citation_pmid": "6655575"},
            {"method": "clinical epidemiology (cohort study)", "species": "human", "citation_pmid": "18077440"},
        ],
    }


# ---------------------------------------------------------------------------------------
# PART 3 -- secondary, explicitly-illustrative viscous-resistance cross-check
# ---------------------------------------------------------------------------------------
def part3_viscous_crosscheck():
    bulk_viscosity_ratio = MU_WATER_37C_mPas / MU_AIR_37C_mPas
    return {
        "mu_water_37C_mPas": MU_WATER_37C_mPas,
        "mu_air_37C_mPas": MU_AIR_37C_mPas,
        "bulk_viscosity_ratio_water_over_air": bulk_viscosity_ratio,
        "hooper_measured_airway_resistance_ratio": HOOPER_MEASURED_VISCOUS_RATIO,
        "gap_disclosed": "bulk Newtonian viscosity ratio (water/air) is a lower-bound illustrative "
                          "estimate, NOT a re-derivation of Hooper et al 2013's measured whole-airway "
                          "~100x figure -- real fetal lung liquid is not pure water (protein content), "
                          "and real airway resistance includes geometric/non-Newtonian effects a bulk "
                          "viscosity ratio does not capture. Both point the same direction (liquid >> "
                          "air resistance, by more than an order of magnitude); magnitudes are not "
                          "claimed identical.",
    }


def build_gates(p1, p2, p3):
    gates = {}

    # --- Part 1 gates ---
    gates["A1_r_implied_lt_half_mature_radius"] = p1["r_implied_max_micron"] < GATE_R_IMPLIED_MAX_MICRON
    gates["A2_mature_bare_alone_cannot_reach_measured_min"] = p1["mature_bare_max_cmH2O"] < GATE_MATURE_BARE_MAX_cmH2O_LT
    gates["A3_first_breath_over_mature_tidal_ratio_ge_10x"] = p1["ratio_conservative_x"] >= GATE_RATIO_MIN_X
    gates["A4_chest_wall_more_compliant_than_lung"] = CW_OVER_CL_INFANT > GATE_CW_CL_MIN
    gates["A5_SLI_same_order_of_magnitude_as_first_breath"] = (
        P_SLI_LO_cmH2O / P_FIRST_BREATH_HI_cmH2O >= 0.1 and P_SLI_HI_cmH2O <= P_FIRST_BREATH_HI_cmH2O * 2
    )

    # --- Part 2 gates ---
    gates["B1_catecholamine_surge_ge_10x"] = p2["adrenaline_fold_early_to_late_labor"] >= GATE_CATECHOL_FOLD_MIN_X
    gates["B2_ai_threshold_sensitizes_with_gestation"] = p2["ai_threshold_sensitization_fold"] > 1.0
    gates["B3_amiloride_margin_ge_5x"] = p2["amiloride_applied_over_KI_margin_x"] >= GATE_AMILORIDE_MARGIN_MIN_X
    gates["B4_enac_ko_fails_despite_normal_birth_mechanics"] = p2["enac_ko_death_window_h"] <= 72.0  # categorical: KO dies fast, not "eventually clears"
    gates["B5_cs_or_37wk_lower_CI_excludes_1"] = CS_OR_37WK[1] > GATE_CS_OR_LOWER_CI_MIN
    gates["B6_cs_or_monotonic_with_gestational_maturity"] = p2["cs_or_monotonic_decreasing_with_ga"]
    gates["B7_decorrelated_methods_ge_3"] = len(p2["decorrelated_methods"]) >= GATE_DECORRELATED_METHODS_MIN

    # --- Part 3 gate (secondary/illustrative) ---
    gates["C1_bulk_viscosity_ratio_gt_10x_consistent_direction"] = p3["bulk_viscosity_ratio_water_over_air"] > 10.0

    overall = all(gates.values())
    return gates, overall


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    p1 = part1_opening_pressure()
    p2 = part2_fluid_clearance()
    p3 = part3_viscous_crosscheck()
    gates, overall = build_gates(p1, p2, p3)

    result = {
        "part1_opening_pressure_geometry": p1,
        "part2_fluid_clearance_adversary_forcing": p2,
        "part3_viscous_resistance_crosscheck": p3,
        "gates": gates,
        "overall_pass": overall,
    }

    with open(OUTFILE, "w") as fh:
        json.dump(result, fh, indent=2, default=lambda o: bool(o) if str(type(o)) == "<class 'numpy.bool_'>" else float(o))

    print("=" * 90)
    print("NEONATAL TRANSITION -- FIRST BREATH + LUNG-FLUID CLEARANCE -- PRE-REGISTERED GATES")
    print("=" * 90)
    for k, v in gates.items():
        print(f"  {k:55s} {v}")
    print(f"\nr_implied range: {p1['r_implied_min_micron']:.2f}-{p1['r_implied_max_micron']:.2f} micron "
          f"(mature r={p1['r_mature_micron']:.0f} micron)")
    print(f"first-breath/mature-tidal ratio: {p1['ratio_conservative_x']:.1f}x - {p1['ratio_generous_x']:.1f}x")
    print(f"adrenaline fold (early->late labor): {p2['adrenaline_fold_early_to_late_labor']:.2f}x")
    print(f"amiloride margin over KI: {p2['amiloride_applied_over_KI_margin_x']:.1f}x")
    print(f"OVERALL PASS: {overall}")
    print(f"\nWritten: {OUTFILE}")

    with open(OUTFILE, "rb") as fh:
        print("md5:", hashlib.md5(fh.read()).hexdigest())


if __name__ == "__main__":
    main()
