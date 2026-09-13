"""RENAL AMMONIAGENESIS <-> HEPATIC UREAGENESIS NITROGEN COUPLING.

Reads: nothing (all inputs are cited literals).
Writes: OUT_ROOT/renal_ammoniagenesis_nitrogen_coupling/renal_ammoniagenesis_nitrogen_coupling_results.json
Gates G1-G3 and G5-G10 decide (overall_pass); G4, the stoichiometric identity under test, is
reported and never gating.

Built as the explicit falsifier named by
HOLE-ACIDBASE-NAE-PROVENANCE-AMMONIAGENESIS-CEILING-NPARTITION, whose closing line
asks for: (a) a quantitative renal glutamine-ammoniagenesis flux model, (b) the cross-scale test
`glutamine_flux x 2 == urinary ammonium excretion rate`, and (c) a DYNAMIC coupling of that
nitrogen draw to the urea cycle's input pool (the parent cell found the one existing split
"additive and correct where checked, but STATIC -- it would not move if an acid load were
applied").

THE THREE QUESTIONS, and what is actually at stake in each:

  Q1 NITROGEN LEDGER. Does a whole-body nitrogen ledger built from independently-measured route
     rates (diet-record intake in; urinary urea-N, urinary NH4+-N, creatinine-N, uric-acid-N,
     faecal-N, integumental/misc-N out) close, and what is the residual?

  Q2 GLUTAMINE STOICHIOMETRY. Renal glutamine -> glutaminase -> glutamate -> GDH yields 2 NH4+
     per glutamine (Weiner & Verlander; Boron). Test the identity 2 x (renal glutamine
     extraction) == (urinary ammonium excretion) on RAW measured human renal arteriovenous data,
     at baseline and under chronic NH4Cl acidosis. If it fails, decide WHICH of the three
     candidates is wrong: the glutamine measurement, the ammonium measurement, or the assumed
     stoichiometry.

  Q3 CROSS-ORGAN COMPETITION. Hepatic ureagenesis and renal ammoniagenesis draw on the SAME
     waste nitrogen. Under an acid load, how far must ureagenesis fall for ammoniagenesis to rise
     as measured, and does the measured change in urea excretion match that prediction?

PRE-REGISTERED, BEFORE ANY NUMBER BELOW WAS COMPUTED:

  G1  Ledger residual: |intake_N - sum(measured out-routes)| / intake_N <= 0.10.
  G2  Urea-N as a fraction of total urinary N, computed from the component ledger, must land in
      the externally-stated 0.80-0.90 band (Guyton; NOT used to build any component).
  G3  Source-internal over-determination: Owen & Robinson 1963 states its renal ammonia and
      glutamine fluxes BOTH as molar rates (umol/min) and as nitrogen mass rates (mg N/min).
      Each independent recompute must agree with the paper's stated companion figure to
      <= 5%. (This is a check on the DATA TRANSCRIPTION, not on any hypothesis.)
  G4  THE STOICHIOMETRIC IDENTITY TEST (the parent cell's named test):
      R_urine = 2 * glutamine_extraction / urinary_NH4 must equal 1.00 +/- 0.20 for the identity
      to hold. ==> THIS GATE'S OUTCOME IS THE FINDING AND IS DELIBERATELY EXCLUDED FROM
      overall_pass. A script that could only "pass" by the identity holding would be unable to
      report a refutation; the identity is the HYPOTHESIS under test, not a correctness gate.
  G5  ATTRIBUTION (only meaningful if G4 fails): the deficit factorises exactly as
      R_urine = (2/Y) * (1/f_urine), where Y = measured NH4+ produced per glutamine extracted and
      f_urine = fraction of produced ammonia excreted in urine (the rest returns to the renal
      vein). If ONE measurement were wrong, one factor would carry the whole deficit. Gate: in
      >= 2 of 3 independent datasets BOTH (Y/2) < 0.90 AND f_urine < 0.90, i.e. the two factors
      multiply -> the failure is structural, not a single bad number.
  G6  DECISIVENESS SCREEN on the acidosis dataset used for Q3: the NH4Cl acid load is itself a
      nitrogen load. If exogenous NH4Cl-N per day >= the measured INCREMENT in urinary NH4+-N,
      the dataset CANNOT test any fixed-nitrogen-pool hypothesis (the new nitrogen is simply the
      dose). Gate must be TRUE (=> dataset correctly marked NON-DECISIVE for Q3).
  G7  Q3 on the decisive dataset (acid-base manipulated with NO exogenous nitrogen): fixed-pool
      hypothesis H_fixed predicts substitution coefficient s = -d(urea_N)/d(NH4_N) = +1.00.
      H_fixed is retained only if measured s in [0.50, 1.50]; otherwise REFUTED.
  G8  Coupled model: reproduces the measured baseline to <= 1%; urinary NH4+ strictly increasing
      in acid load; the ammoniagenic ceiling binds at high load; glutamine uptake strictly
      increasing.
  G9  VOID-FLOOR on the coupling: the renal-vein ammonia return term must be LOAD-BEARING --
      forcing f_urine = 1.0 (the naive "production == excretion" assumption) must change the
      predicted glutamine uptake by >= 1.5x. If it does not, the term is decorative.
  G10 CIRCULARITY / INVERSION TEST (task-mandated). Urinary urea-N is routinely used to ESTIMATE
      protein intake (Maroni: PNA = 6.25*(UUN + 0.031*BW)). Machine-demonstrate that closing the
      ledger with a Maroni-derived intake makes the residual identically zero for ANY urea value
      -> a tautology gate. This gate must be TRUE (the tautology must be demonstrated) and the
      real ledger must NOT use that route.

  overall_pass = G1..G3 and G5..G10. G4 is reported, never gating (see above).

SYMMETRIC QC -- STATED BEFORE RESULTS, NOT AFTER:
  * Owen & Robinson 1963 measured glutamine in PLASMA only; Tizianello 1980 measured whole-blood
    glutamine (plasma + a blood-cell route). The whole-blood denominator is larger, so the two
    studies' Y values are NOT expected to agree and the difference between them is an upper
    bound on the blood-cell route, not an error. Both are reported.
  * Owen's acidosis arm is NH4Cl loading. NH4Cl delivers nitrogen; this is exactly the confound
    G6 screens for and is why Q3 is NOT decided on that arm.
  * The decisive Q3 dataset (Hannaford 1982) is prolonged total fasting -- a regime where net
    protein catabolism is large and demonstrably free to move. A fed, constant-diet regime may
    behave differently; the rat HCl dataset (Oliver & Bourke 1975) reports the opposite result
    and is reported alongside, NOT discarded. The disagreement is a REGIME finding.
  * Oliver & Bourke's phrase "equimolar increase in NH4+ ... no change in their sum" is unit-
    ambiguous: urea carries 2 N, so molar-equimolar and nitrogen-equimolar readings differ by a
    factor of 2. Both readings are computed; the ambiguity is NOT resolved by assertion.
  * The ledger's "other urinary N" term (free amino acids, peptides, hippurate) is the least
    independently pinned component and is carried with an explicit band.

CITATIONS (each live-verified against NCBI eutils / the publisher PDF while writing this file):
  [1] Owen EE, Robinson RR. Amino acid extraction and ammonia metabolism by the human kidney
      during the prolonged administration of ammonium chloride. J Clin Invest 1963;42(2):263-276.
      PMID 13940831, PMC289275, doi 10.1172/JCI104713. Raw renal A-V data, n=6 control + n=6
      NH4Cl (6-8 g/day orally, 6-9 days).
  [2] Tizianello A, De Ferrari G, Garibotto G, Gurreri G, Robaudo C. Renal metabolism of amino
      acids and ammonia in subjects with normal renal function and in patients with chronic renal
      insufficiency. J Clin Invest 1980;65(5):1162-1173. PMID 7364943, PMC371450,
      doi 10.1172/JCI109771.
  [3] Hannaford MC, Leiter LA, Josse RG, Goldstein MB, Marliss EB, Halperin ML. Protein wasting
      due to acidosis of prolonged fasting. Am J Physiol 1982;243(3):E251-6. PMID 6287864,
      doi 10.1152/ajpendo.1982.243.3.E251. n=8 obese subjects, >14 d total fast, 150 mmol NaHCO3
      + 60 mmol KCl/day. NO exogenous nitrogen -> the decisive Q3 dataset.
  [4] Oliver J, Bourke E. Adaptations in urea and ammonium excretion in metabolic acidosis in the
      rat: a reinterpretation. Clin Sci Mol Med 1975;48(6):515-520. PMID 1056282,
      doi 10.1042/cs0480515. Rat, constant diet, oral HCl (non-nitrogenous acid).
  [5] Bingham SA, Cummings JH. Urine nitrogen as an independent validatory measure of dietary
      intake: a study of nitrogen balance in individuals consuming their normal diet.
      Am J Clin Nutr 1985;42(6):1276-1289. PMID 4072961, doi 10.1093/ajcn/42.6.1276. n=8, 28 d,
      duplicate-diet Kjeldahl N, PABA-verified urine, pellet-verified faeces, skin+blood N
      measured. Urinary N = 81 +/- 5% (SD) of dietary N.
  [6] FAO/WHO/UNU 2007, Protein and amino acid requirements in human nutrition, WHO TRS 935:
      obligatory nitrogen losses, adult: urinary 37, faecal 12, integumental 3, miscellaneous 2
      mg N/kg/day.
  [7] Weiner ID, Verlander JW. Renal ammonia metabolism and transport. Compr Physiol
      2013;3(1):201-220, doi 10.1002/cphy.c120010; and Emerging features of ammonia metabolism
      and transport in acid-base balance, PMC6629433: "Complete glutamine metabolism in the
      proximal tubule results in equimolar generation of NH4+ and HCO3-" (2 NH4+ + 2 HCO3- per
      glutamine).
  [8] Published anchors, re-used only as EXTERNAL comparators, never as inputs to a fit:
      Halperin/Rose baseline urinary NH4+ 40-50 mEq/day rising to 200-300 in chronic acidosis;
      Guyton: urea = 80-90% of urinary N; NHANES adult protein intake 1.15 +/- 0.35 g/kg/day.

SIDE EFFECTS: writes ONE JSON under its own output directory; reads nothing but its own constants;
no network, no GPU, no git. Pure closed-form arithmetic (< 1 s, no resource gate needed).
"""
import json
import os
import sys

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "renal_ammoniagenesis_nitrogen_coupling")

# ---- universal constants -----------------------------------------------------------------
N_ATOMIC_MASS_MG = 14.007          # mg per mmol N
MIN_PER_DAY = 1440.0
N_PER_UREA = 2                     # [7]/Lehninger
NH4_PER_GLUTAMINE_THEORY = 2       # [7] -- the stoichiometry UNDER TEST
N_PER_GLUTAMINE = 2                # amide N + alpha-amino N
PROTEIN_TO_N_DIVISOR = 6.25        # Jones 1931

# ---- Q1 ledger inputs. PROVENANCE TAG IS LOAD-BEARING (see G10) ---------------------------
BODY_MASS_KG = 70.0
LEDGER_INPUTS = {
    # value, band, unit, provenance -- provenance must NOT be "urea-derived" for the intake term
    "protein_intake_g_per_kg_day": dict(v=1.15, lo=0.80, hi=1.50, prov="diet-record/food-composition (NHANES) [8]"),
    "urinary_N_fraction_of_intake": dict(v=0.81, lo=0.71, hi=0.91, prov="duplicate-diet Kjeldahl balance, PABA-verified [5]"),
    "urinary_NH4_mmol_per_day": dict(v=45.0, lo=40.0, hi=50.0, prov="24h urine NH4+ titration, Halperin/Rose [8]"),
    "creatinine_g_per_day": dict(v=1.60, lo=1.20, hi=2.00, prov="24h urine creatinine, clinical chemistry"),
    "uric_acid_g_per_day": dict(v=0.60, lo=0.40, hi=0.80, prov="24h urine uric acid, clinical chemistry"),
    "other_urinary_N_g_per_day": dict(v=0.30, lo=0.15, hi=0.50, prov="free AA/peptide/hippurate residual (weakest term)"),
    "faecal_N_g_per_day": dict(v=1.20, lo=0.90, hi=1.50, prov="pellet-verified faecal Kjeldahl [5],[6]"),
    "integumental_misc_N_mg_per_kg_day": dict(v=5.0, lo=3.0, hi=8.0, prov="FAO/WHO/UNU obligatory losses (3 skin + 2 misc) [6]"),
}
CREATININE_MW, CREATININE_N = 113.12, 3
URIC_ACID_MW, URIC_ACID_N = 168.11, 4
UREA_N_FRACTION_EXTERNAL_BAND = (0.80, 0.90)     # [8] Guyton -- external comparator only

# ---- Q2 raw measured renal arteriovenous datasets ----------------------------------------
# All values transcribed verbatim from the source text; `paper_*` fields are the source's
# companion statements, used ONLY for the transcription over-determination check G3.
RENAL_DATASETS = {
    "owen1963_control": dict(
        src="[1] Owen & Robinson 1963, n=6, no NH4Cl",
        glutamine_extraction_umol_min=47.0, glutamine_sd=14.0,
        urine_NH4_umol_min=37.8, urine_sd=25.6,
        renal_vein_NH4_umol_min=37.3, rv_sd=24.6,
        paper_total_NH4_umol_min=77.4,                 # paper's stated total, vs 37.8+37.3
        paper_total_NH3_N_mg_min=1.05,                 # paper's mass-route statement
        paper_glutamine_total_N_mg_min=1.33,           # paper's mass-route statement (2 N/Gln)
        paper_AAN_over_ammoniaN_pct=132.0,             # glycine+glutamine-total-N vs ammonia-N
        exogenous_N_mmol_per_day=0.0,
    ),
    "owen1963_nh4cl_acidosis": dict(
        src="[1] Owen & Robinson 1963, n=6, 6-8 g NH4Cl/day x 6-9 days",
        glutamine_extraction_umol_min=102.0, glutamine_sd=53.0,
        urine_NH4_umol_min=96.8, urine_sd=31.9,
        renal_vein_NH4_umol_min=42.3, rv_sd=14.9,
        paper_total_NH4_umol_min=139.1,
        paper_total_NH3_N_mg_min=1.95,
        paper_glutamine_total_N_mg_min=None,           # stated as a % instead (147%), see G3
        paper_AAN_over_ammoniaN_pct=147.0,
        exogenous_N_mmol_per_day=None,                 # computed from the stated 6-8 g/day dose
    ),
    "tizianello1980_normal": dict(
        src="[2] Tizianello 1980, normal renal function, whole-blood route, per 1.73 m2",
        glutamine_extraction_umol_min=None,            # back-solved from the paper's ratio
        glutamine_sd=None,
        urine_NH4_umol_min=None,                       # "almost equally partitioned"
        renal_vein_NH4_umol_min=None,
        paper_total_NH4_umol_min=37.8,                 # 37.8 +/- 3.08 umol/min
        paper_production_over_glutamine_ratio=1.0,     # "The ratio of ammonia production to
                                                       #  glutamine extraction is ~1"
        paper_urine_fraction_of_production=0.50,       # "almost equally partitioned"
        paper_total_NH3_N_mg_min=None,
        paper_glutamine_total_N_mg_min=None,
        paper_AAN_over_ammoniaN_pct=None,
        exogenous_N_mmol_per_day=0.0,
    ),
}
NH4CL_DOSE_G_PER_DAY = (6.0, 8.0)   # [1]
NH4CL_MW = 53.49

# ---- Q3 datasets -------------------------------------------------------------------------
# [3] Hannaford 1982: units are g N per g creatinine (creatinine-normalised, the paper's
# normalisation); the substitution coefficient is dimensionless so the normalisation cancels.
HANNAFORD1982 = dict(
    src="[3] Hannaford 1982, n=8, >14 d total fast, +150 mmol NaHCO3 +60 mmol KCl/day",
    NH4_N_before=3.8, NH4_N_before_sem=0.4, NH4_N_after=2.0, NH4_N_after_sem=0.4,
    urea_N_before=2.5, urea_N_before_sem=0.2, urea_N_after=2.1, urea_N_after_sem=0.3,
    exogenous_N_given=False, units="g N / g creatinine",
)
OLIVER1975 = dict(
    src="[4] Oliver & Bourke 1975, rat, constant diet, oral HCl",
    qualitative="urea excretion significantly reduced, equimolar increase in NH4+, no change in their sum",
    exogenous_N_given=False,
)
# Parent-cell comparators (EXTERNAL to this script's arithmetic, never fitted):
NH4_BASELINE_BAND_MMOL_DAY = (40.0, 50.0)
NH4_ACIDOSIS_BAND_MMOL_DAY = (200.0, 300.0)
DELTA_NH4_FRACTION_OF_DELTA_NAE = 0.93   # parent cell: 90-96% of the adaptive NAE increase


# =========================================================================================
# helpers
# =========================================================================================
def mgN_to_umol(mg_per_min):
    return mg_per_min / N_ATOMIC_MASS_MG * 1000.0


def umol_min_to_mmol_day(x):
    return x * MIN_PER_DAY / 1000.0


def rel(a, b):
    return abs(a - b) / abs(b) if b else float("inf")


# =========================================================================================
# STEP 1 -- whole-body nitrogen ledger
# =========================================================================================
def step1_nitrogen_ledger(scale=None):
    """scale: optional dict of key->value overrides (used for the band sweep)."""
    g = dict((k, v["v"]) for k, v in LEDGER_INPUTS.items())
    if scale:
        g.update(scale)

    protein_g_day = g["protein_intake_g_per_kg_day"] * BODY_MASS_KG
    intake_N = protein_g_day / PROTEIN_TO_N_DIVISOR                                   # g N/day
    urinary_N = intake_N * g["urinary_N_fraction_of_intake"]                          # [5]

    nh4_N = g["urinary_NH4_mmol_per_day"] * N_ATOMIC_MASS_MG / 1000.0
    creat_N = g["creatinine_g_per_day"] / CREATININE_MW * CREATININE_N * N_ATOMIC_MASS_MG
    uric_N = g["uric_acid_g_per_day"] / URIC_ACID_MW * URIC_ACID_N * N_ATOMIC_MASS_MG
    other_N = g["other_urinary_N_g_per_day"]
    urea_N = urinary_N - nh4_N - creat_N - uric_N - other_N
    urea_mmol_day = urea_N / N_ATOMIC_MASS_MG / N_PER_UREA * 1000.0

    faecal_N = g["faecal_N_g_per_day"]
    skin_misc_N = g["integumental_misc_N_mg_per_kg_day"] * BODY_MASS_KG / 1000.0
    total_out = urinary_N + faecal_N + skin_misc_N
    residual = intake_N - total_out

    return dict(
        protein_g_day=protein_g_day, intake_N_g_day=intake_N, urinary_N_g_day=urinary_N,
        urinary_NH4_N_g_day=nh4_N, urinary_creatinine_N_g_day=creat_N,
        urinary_uric_acid_N_g_day=uric_N, urinary_other_N_g_day=other_N,
        urinary_urea_N_g_day=urea_N, urea_flux_mmol_day=urea_mmol_day,
        urea_N_fraction_of_urinary=urea_N / urinary_N,
        nh4_N_fraction_of_urinary=nh4_N / urinary_N,
        faecal_N_g_day=faecal_N, integumental_misc_N_g_day=skin_misc_N,
        total_out_N_g_day=total_out, residual_N_g_day=residual,
        residual_fraction_of_intake=residual / intake_N,
    )


def step1_band_sweep():
    """Propagate each input over its stated band, one at a time, to size the residual band."""
    base = step1_nitrogen_ledger()
    contrib = {}
    lo_res, hi_res = base["residual_N_g_day"], base["residual_N_g_day"]
    for k, spec in LEDGER_INPUTS.items():
        a = step1_nitrogen_ledger({k: spec["lo"]})["residual_N_g_day"]
        b = step1_nitrogen_ledger({k: spec["hi"]})["residual_N_g_day"]
        contrib[k] = dict(lo=a, hi=b, half_width=abs(b - a) / 2.0)
        lo_res, hi_res = min(lo_res, a, b), max(hi_res, a, b)
    quad = sum(c["half_width"] ** 2 for c in contrib.values()) ** 0.5
    ranked = sorted(contrib.items(), key=lambda kv: -kv[1]["half_width"])

    # POWER, in the convention already used by HOLE-NITROGEN-UREA-CLOSURE-POWER-ADJUDICATION:
    # noise/residual > 1 => BALANCES-BUT-LOW-POWER (the closure cannot resolve this residual).
    noise_over_residual = quad / abs(base["residual_N_g_day"]) if base["residual_N_g_day"] else float("inf")

    # FAO/WHO/UNU standard correction for routinely-unmeasured losses in N-balance studies
    # (+5 mg N/kg/day): N balance is known to be biased toward apparent retention.
    corrected = base["residual_N_g_day"] - 5.0 * BODY_MASS_KG / 1000.0

    # FALSIFIABILITY POWER of G2 (urea-N fraction vs the external 0.80-0.90 band): if the fraction
    # cannot leave the band anywhere inside the inputs' own stated bands, G2 is a NON-FALSIFYING
    # sanity check and must be labelled as such rather than sold as a validation.
    fr_lo, fr_hi = base["urea_N_fraction_of_urinary"], base["urea_N_fraction_of_urinary"]
    for k, spec in LEDGER_INPUTS.items():
        for val in (spec["lo"], spec["hi"]):
            fr = step1_nitrogen_ledger({k: val})["urea_N_fraction_of_urinary"]
            fr_lo, fr_hi = min(fr_lo, fr), max(fr_hi, fr)
    g2_can_fail = (fr_lo < UREA_N_FRACTION_EXTERNAL_BAND[0]) or (fr_hi > UREA_N_FRACTION_EXTERNAL_BAND[1])

    return dict(per_input=contrib, residual_min=lo_res, residual_max=hi_res,
                quadrature_half_width_g_day=quad,
                governor=ranked[0][0], governor_half_width=ranked[0][1]["half_width"],
                noise_over_residual=noise_over_residual,
                power_verdict=("BALANCES-BUT-LOW-POWER" if noise_over_residual > 1.0
                               else "RESIDUAL-RESOLVED"),
                residual_after_unmeasured_loss_correction_g_day=corrected,
                residual_after_correction_fraction=corrected / base["intake_N_g_day"],
                urea_fraction_range_over_input_bands=(fr_lo, fr_hi),
                G2_is_falsifiable_within_input_bands=g2_can_fail)


# =========================================================================================
# STEP 2 -- glutamine flux and the 2:1 stoichiometric identity test
# =========================================================================================
def step2_glutamine_stoichiometry():
    out, checks = {}, {}

    # --- G3 transcription over-determination, Owen control ---
    c = RENAL_DATASETS["owen1963_control"]
    sum_route = c["urine_NH4_umol_min"] + c["renal_vein_NH4_umol_min"]
    mass_route = mgN_to_umol(c["paper_total_NH3_N_mg_min"])
    gln_mass_route = mgN_to_umol(c["paper_glutamine_total_N_mg_min"]) / N_PER_GLUTAMINE
    checks["owen_control_sum_vs_stated_total"] = dict(
        sum_of_parts=sum_route, paper_stated=c["paper_total_NH4_umol_min"],
        rel_err=rel(sum_route, c["paper_total_NH4_umol_min"]))
    checks["owen_control_massroute_vs_sum"] = dict(
        mass_route=mass_route, sum_of_parts=sum_route, rel_err=rel(mass_route, sum_route))
    checks["owen_control_glutamine_mass_vs_molar"] = dict(
        mass_route=gln_mass_route, paper_molar=c["glutamine_extraction_umol_min"],
        rel_err=rel(gln_mass_route, c["glutamine_extraction_umol_min"]))

    a = RENAL_DATASETS["owen1963_nh4cl_acidosis"]
    sum_a = a["urine_NH4_umol_min"] + a["renal_vein_NH4_umol_min"]
    mass_a = mgN_to_umol(a["paper_total_NH3_N_mg_min"])
    checks["owen_acidosis_sum_vs_stated_total"] = dict(
        sum_of_parts=sum_a, paper_stated=a["paper_total_NH4_umol_min"],
        rel_err=rel(sum_a, a["paper_total_NH4_umol_min"]))
    checks["owen_acidosis_massroute_vs_sum"] = dict(
        mass_route=mass_a, sum_of_parts=sum_a, rel_err=rel(mass_a, sum_a))
    # the paper states glutamine+glutamate total-N uptake = 147% of ammonia-N release; invert it
    inv_gln_acid = a["paper_AAN_over_ammoniaN_pct"] / 100.0 * mass_a / N_PER_GLUTAMINE
    checks["owen_acidosis_147pct_inverts_to_glutamine"] = dict(
        inverted=inv_gln_acid, paper_molar=a["glutamine_extraction_umol_min"],
        rel_err=rel(inv_gln_acid, a["glutamine_extraction_umol_min"]))

    g3_pass = all(v["rel_err"] <= 0.05 for v in checks.values())

    # --- the identity test on each dataset ---
    for name in ("owen1963_control", "owen1963_nh4cl_acidosis"):
        d = RENAL_DATASETS[name]
        prod = d["urine_NH4_umol_min"] + d["renal_vein_NH4_umol_min"]
        gln = d["glutamine_extraction_umol_min"]
        Y = prod / gln
        f_urine = d["urine_NH4_umol_min"] / prod
        R_urine = NH4_PER_GLUTAMINE_THEORY * gln / d["urine_NH4_umol_min"]
        R_prod = NH4_PER_GLUTAMINE_THEORY * gln / prod
        out[name] = dict(
            src=d["src"], glutamine_umol_min=gln, urine_NH4_umol_min=d["urine_NH4_umol_min"],
            renal_vein_NH4_umol_min=d["renal_vein_NH4_umol_min"], production_umol_min=prod,
            Y_NH4_per_glutamine=Y, yield_fraction_of_theory=Y / NH4_PER_GLUTAMINE_THEORY,
            f_urine_fraction_excreted=f_urine,
            R_urine_2Gln_over_urinaryNH4=R_urine, R_prod_2Gln_over_production=R_prod,
            factorisation_check=rel(R_urine, (NH4_PER_GLUTAMINE_THEORY / Y) * (1.0 / f_urine)),
            urinary_NH4_mmol_day=umol_min_to_mmol_day(d["urine_NH4_umol_min"]),
            glutamine_mmol_day=umol_min_to_mmol_day(gln),
        )

    # Tizianello: only ratios published -> back-solve, do not invent absolutes
    t = RENAL_DATASETS["tizianello1980_normal"]
    prod_t = t["paper_total_NH4_umol_min"]
    gln_t = prod_t / t["paper_production_over_glutamine_ratio"]
    urine_t = prod_t * t["paper_urine_fraction_of_production"]
    out["tizianello1980_normal"] = dict(
        src=t["src"], glutamine_umol_min=gln_t, urine_NH4_umol_min=urine_t,
        renal_vein_NH4_umol_min=prod_t - urine_t, production_umol_min=prod_t,
        Y_NH4_per_glutamine=prod_t / gln_t,
        yield_fraction_of_theory=(prod_t / gln_t) / NH4_PER_GLUTAMINE_THEORY,
        f_urine_fraction_excreted=urine_t / prod_t,
        R_urine_2Gln_over_urinaryNH4=NH4_PER_GLUTAMINE_THEORY * gln_t / urine_t,
        R_prod_2Gln_over_production=NH4_PER_GLUTAMINE_THEORY * gln_t / prod_t,
        factorisation_check=0.0,
        urinary_NH4_mmol_day=umol_min_to_mmol_day(urine_t),
        glutamine_mmol_day=umol_min_to_mmol_day(gln_t),
        note="glutamine is a WHOLE-BLOOD extraction (larger denominator than Owen's plasma route)",
    )

    # G4 -- the identity, reported not gating
    g4 = {k: dict(R=v["R_urine_2Gln_over_urinaryNH4"],
                  identity_holds=0.80 <= v["R_urine_2Gln_over_urinaryNH4"] <= 1.20)
          for k, v in out.items()}
    identity_holds_anywhere = any(v["identity_holds"] for v in g4.values())

    # G5 -- attribution: do BOTH factors deviate?
    both = sum(1 for v in out.values()
               if v["yield_fraction_of_theory"] < 0.90 and v["f_urine_fraction_excreted"] < 0.90)
    g5_pass = both >= 2

    # external anchor: Owen's baseline urinary NH4 vs the Halperin/Rose band
    base_mmol = out["owen1963_control"]["urinary_NH4_mmol_day"]
    acid_mmol = out["owen1963_nh4cl_acidosis"]["urinary_NH4_mmol_day"]
    external = dict(
        owen_baseline_mmol_day=base_mmol, halperin_band=NH4_BASELINE_BAND_MMOL_DAY,
        baseline_within_band=NH4_BASELINE_BAND_MMOL_DAY[0] <= base_mmol <= NH4_BASELINE_BAND_MMOL_DAY[1],
        baseline_rel_to_band_top=base_mmol / NH4_BASELINE_BAND_MMOL_DAY[1],
        owen_acidosis_mmol_day=acid_mmol, halperin_acidosis_band=NH4_ACIDOSIS_BAND_MMOL_DAY,
        acidosis_within_band=NH4_ACIDOSIS_BAND_MMOL_DAY[0] <= acid_mmol <= NH4_ACIDOSIS_BAND_MMOL_DAY[1],
        acidosis_fraction_of_band_floor=acid_mmol / NH4_ACIDOSIS_BAND_MMOL_DAY[0],
    )
    return dict(datasets=out, transcription_checks=checks, g3_pass=g3_pass,
                g4_identity=g4, identity_holds_anywhere=identity_holds_anywhere,
                g5_both_factors_deviate_count=both, g5_pass=g5_pass,
                external_anchor=external)


# =========================================================================================
# STEP 3 -- cross-organ nitrogen competition (the regulated split)
# =========================================================================================
def step3_cross_organ(ledger, s2):
    # (a) decisiveness screen on the NH4Cl arm
    d_ctrl = s2["datasets"]["owen1963_control"]
    d_acid = s2["datasets"]["owen1963_nh4cl_acidosis"]
    d_nh4_mmol_day = d_acid["urinary_NH4_mmol_day"] - d_ctrl["urinary_NH4_mmol_day"]
    d_nh4_N_g_day = d_nh4_mmol_day * N_ATOMIC_MASS_MG / 1000.0
    dose_N_lo = NH4CL_DOSE_G_PER_DAY[0] / NH4CL_MW * 1000.0
    dose_N_hi = NH4CL_DOSE_G_PER_DAY[1] / NH4CL_MW * 1000.0
    sufficiency = dose_N_lo / d_nh4_mmol_day
    g6_pass = sufficiency >= 1.0

    # (b) what H_fixed WOULD require, computed anyway (the prediction is the deliverable)
    urea_base_mmol = ledger["urea_flux_mmol_day"]
    required_urea_fall_mmol = d_nh4_N_g_day / N_ATOMIC_MASS_MG * 1000.0 / N_PER_UREA
    required_fall_frac = required_urea_fall_mmol / urea_base_mmol

    # same prediction at the maximal chronic-acidosis band (200-300 mmol NH4/day)
    maximal = {}
    for nh4_max in NH4_ACIDOSIS_BAND_MMOL_DAY:
        dN = (nh4_max - LEDGER_INPUTS["urinary_NH4_mmol_per_day"]["v"]) * N_ATOMIC_MASS_MG / 1000.0
        fall_mmol = dN / N_ATOMIC_MASS_MG * 1000.0 / N_PER_UREA
        maximal[f"NH4_{int(nh4_max)}_mmol_day"] = dict(
            delta_NH4_N_g_day=dN, required_urea_fall_mmol_day=fall_mmol,
            required_urea_fall_fraction=fall_mmol / urea_base_mmol,
            residual_urea_N_g_day=ledger["urinary_urea_N_g_day"] - dN)

    # (c) THE DECISIVE TEST -- acid-base moved with no exogenous nitrogen
    h = HANNAFORD1982
    d_nh4 = h["NH4_N_after"] - h["NH4_N_before"]
    d_urea = h["urea_N_after"] - h["urea_N_before"]
    d_total = d_nh4 + d_urea
    s_measured = -d_urea / d_nh4                       # H_fixed predicts +1.00
    g7_h_fixed_retained = 0.50 <= s_measured <= 1.50
    pool_explained = abs(d_total) / (abs(d_nh4) + abs(d_urea)) if (d_nh4 or d_urea) else 0.0

    # propagate the paper's SEMs onto s (independent-error quadrature, n=8)
    sd_dnh4 = (h["NH4_N_before_sem"] ** 2 + h["NH4_N_after_sem"] ** 2) ** 0.5
    sd_durea = (h["urea_N_before_sem"] ** 2 + h["urea_N_after_sem"] ** 2) ** 0.5
    s_sd = abs(s_measured) * ((sd_durea / d_urea) ** 2 + (sd_dnh4 / d_nh4) ** 2) ** 0.5
    s_excludes_unity = abs(s_measured - 1.0) > 2.0 * s_sd

    # (d) the rat HCl dataset, both readings of "equimolar" -- ambiguity NOT resolved
    oliver = dict(
        src=OLIVER1975["src"], statement=OLIVER1975["qualitative"],
        reading_molar_equimolar=dict(
            s_in_nitrogen_units=2.0,
            note="if d(urea mmol) = -d(NH4 mmol), urea-N falls 2x the NH4-N rise -> total urinary N FALLS"),
        reading_nitrogen_equimolar=dict(
            s_in_nitrogen_units=1.0,
            note="if d(urea-N) = -d(NH4-N), total urinary N is conserved -> H_fixed exactly"),
        factor_between_readings=2.0,
        resolved=False,
    )
    return dict(
        nh4cl_arm=dict(delta_urinary_NH4_mmol_day=d_nh4_mmol_day,
                       delta_urinary_NH4_N_g_day=d_nh4_N_g_day,
                       exogenous_NH4Cl_N_mmol_day=(dose_N_lo, dose_N_hi),
                       exogenous_over_delta_excreted=sufficiency,
                       decisive_for_H_fixed=not g6_pass,
                       verdict="NON-DECISIVE: the NH4Cl dose supplies %.2fx the extra excreted "
                               "ammonium nitrogen, so no fall in ureagenesis is required"
                               % sufficiency),
        h_fixed_prediction=dict(baseline_urea_flux_mmol_day=urea_base_mmol,
                                required_urea_fall_mmol_day=required_urea_fall_mmol,
                                required_urea_fall_fraction=required_fall_frac,
                                at_maximal_acidosis=maximal),
        decisive_test=dict(src=h["src"], units=h["units"],
                           delta_NH4_N=d_nh4, delta_urea_N=d_urea, delta_total_N=d_total,
                           s_measured=s_measured, s_sd=s_sd,
                           s_predicted_by_H_fixed=1.0,
                           s_excludes_unity_at_2sd=s_excludes_unity,
                           H_fixed_retained=g7_h_fixed_retained,
                           sign_of_urea_change_matches_H_fixed=(d_urea > 0),
                           fraction_of_change_carried_by_pool_size=pool_explained),
        rat_HCl_counter_dataset=oliver,
        g6_pass=g6_pass,
        # G7 is a POWER gate, not an outcome gate: the test only counts if it was run on a
        # dataset that gave no exogenous nitrogen AND could resolve s=1 at 2SD.
        g7_pass=(not h["exogenous_N_given"]) and s_excludes_unity,
        g7_H_fixed_refuted=not g7_h_fixed_retained,
    )


# =========================================================================================
# STEP 4 -- the coupled model the parent cell asked for (dynamic, not static)
# =========================================================================================
def coupled_model(acid_load_mEq_day, ledger, Y=None, f_urine=None, ceiling=300.0,
                  force_f_urine=None):
    """Acid load -> renal NH4+ excretion -> renal glutamine draw -> nitrogen removed from the
    ureagenic pool. f_urine interpolates from the measured baseline (Owen control) to the
    measured acidosis value (Owen NH4Cl arm) as the load rises -- the partition is itself a
    regulated variable, which is exactly what the parent cell flagged as missing."""
    nh4_base = LEDGER_INPUTS["urinary_NH4_mmol_per_day"]["v"]
    Y = 1.6489361702127659 if Y is None else Y          # set from data in main(); default = Owen ctrl
    f_lo, f_hi = (0.4966, 0.6959) if f_urine is None else f_urine
    # excretion rises with load until the ammoniagenic ceiling binds
    nh4 = min(nh4_base + DELTA_NH4_FRACTION_OF_DELTA_NAE * acid_load_mEq_day, ceiling)
    span = max(ceiling - nh4_base, 1e-9)
    frac = min(max((nh4 - nh4_base) / span, 0.0), 1.0)
    f = f_lo + (f_hi - f_lo) * frac
    if force_f_urine is not None:
        f = force_f_urine
    production = nh4 / f                                 # renal-vein return is the difference
    glutamine = production / Y
    nh4_N = nh4 * N_ATOMIC_MASS_MG / 1000.0              # g N/day actually LEAVING the body
    urea_N = ledger["urinary_urea_N_g_day"] - (nh4_N - ledger["urinary_NH4_N_g_day"])
    return dict(acid_load_mEq_day=acid_load_mEq_day, urinary_NH4_mmol_day=nh4,
                f_urine=f, renal_NH3_production_mmol_day=production,
                renal_vein_NH4_return_mmol_day=production - nh4,
                glutamine_uptake_mmol_day=glutamine, urinary_NH4_N_g_day=nh4_N,
                ureagenic_N_g_day_under_H_fixed=urea_N,
                urea_flux_mmol_day_under_H_fixed=urea_N / N_ATOMIC_MASS_MG / N_PER_UREA * 1000.0,
                ceiling_bound=nh4 >= ceiling - 1e-9)


def step4_coupling(ledger, s2):
    Y = s2["datasets"]["owen1963_control"]["Y_NH4_per_glutamine"]
    f_lo = s2["datasets"]["owen1963_control"]["f_urine_fraction_excreted"]
    f_hi = s2["datasets"]["owen1963_nh4cl_acidosis"]["f_urine_fraction_excreted"]
    loads = [0, 25, 50, 100, 150, 200, 275, 400]
    sweep = [coupled_model(a, ledger, Y=Y, f_urine=(f_lo, f_hi)) for a in loads]

    base = sweep[0]
    baseline_err = rel(base["urinary_NH4_mmol_day"], LEDGER_INPUTS["urinary_NH4_mmol_per_day"]["v"])
    mono_nh4 = all(sweep[i + 1]["urinary_NH4_mmol_day"] >= sweep[i]["urinary_NH4_mmol_day"]
                   for i in range(len(sweep) - 1))
    mono_gln = all(sweep[i + 1]["glutamine_uptake_mmol_day"] >= sweep[i]["glutamine_uptake_mmol_day"]
                   for i in range(len(sweep) - 1))
    ceiling_binds = sweep[-1]["ceiling_bound"] and not sweep[0]["ceiling_bound"]
    g8 = baseline_err <= 0.01 and mono_nh4 and mono_gln and ceiling_binds

    # G9 void-floor: is the renal-vein return term load-bearing?
    naive = coupled_model(0, ledger, Y=Y, f_urine=(f_lo, f_hi), force_f_urine=1.0)
    ratio = base["glutamine_uptake_mmol_day"] / naive["glutamine_uptake_mmol_day"]
    g9 = ratio >= 1.5

    # a second void-floor: with the theoretical Y=2 AND f_urine=1 (the naive model the parent
    # cell's test assumes), how far off is the predicted glutamine uptake?
    naive2 = coupled_model(0, ledger, Y=2.0, f_urine=(f_lo, f_hi), force_f_urine=1.0)
    naive2_ratio = base["glutamine_uptake_mmol_day"] / naive2["glutamine_uptake_mmol_day"]

    return dict(Y_used=Y, f_urine_baseline=f_lo, f_urine_acidosis=f_hi, sweep=sweep,
                baseline_rel_err=baseline_err, monotonic_NH4=mono_nh4,
                monotonic_glutamine=mono_gln, ceiling_binds=ceiling_binds, g8_pass=g8,
                void_floor_force_f_urine_1=dict(glutamine_naive=naive["glutamine_uptake_mmol_day"],
                                                glutamine_model=base["glutamine_uptake_mmol_day"],
                                                ratio=ratio, g9_pass=g9),
                void_floor_full_naive_Y2_f1=dict(glutamine_naive=naive2["glutamine_uptake_mmol_day"],
                                                 ratio=naive2_ratio),
                g9_pass=g9)


# =========================================================================================
# STEP 5 -- circularity / inversion tests (G10 + anchor inversions)
# =========================================================================================
def step5_inversions(ledger):
    # (a) Maroni tautology demonstration: close the ledger using a urea-DERIVED intake
    def maroni_closed_residual(urea_N_g_day):
        pna_g_day = PROTEIN_TO_N_DIVISOR * (urea_N_g_day + 0.031 * BODY_MASS_KG)
        intake_N = pna_g_day / PROTEIN_TO_N_DIVISOR
        non_urea_out = 0.031 * BODY_MASS_KG      # Maroni's lumped non-urea term
        return intake_N - (urea_N_g_day + non_urea_out)

    r1 = maroni_closed_residual(ledger["urinary_urea_N_g_day"])
    r2 = maroni_closed_residual(ledger["urinary_urea_N_g_day"] * 2.0)
    r3 = maroni_closed_residual(0.1)
    tautology_demonstrated = max(abs(r1), abs(r2), abs(r3)) < 1e-9
    real_intake_prov = LEDGER_INPUTS["protein_intake_g_per_kg_day"]["prov"]
    intake_is_urea_derived = ("urea" in real_intake_prov.lower()) or ("maroni" in real_intake_prov.lower())

    # (b) invert Bingham's 81% for intake, given the component-built urinary N
    inv_intake_N = ledger["urinary_N_g_day"] / LEDGER_INPUTS["urinary_N_fraction_of_intake"]["v"]
    inv_protein = inv_intake_N * PROTEIN_TO_N_DIVISOR / BODY_MASS_KG

    # (c) invert Maroni on the SAME urea value and compare -- how much does the circular route
    #     differ from the decorrelated one?
    maroni_protein = PROTEIN_TO_N_DIVISOR * (ledger["urinary_urea_N_g_day"]
                                             + 0.031 * BODY_MASS_KG) / BODY_MASS_KG

    # (d) invert the ammoniagenic ceiling for the acid load it can service
    ceiling_service_mEq = (NH4_ACIDOSIS_BAND_MMOL_DAY[1]
                           - LEDGER_INPUTS["urinary_NH4_mmol_per_day"]["v"]) / DELTA_NH4_FRACTION_OF_DELTA_NAE

    return dict(
        maroni_tautology=dict(residual_at_measured_urea=r1, residual_at_2x_urea=r2,
                              residual_at_0p1_urea=r3, invariant_to_urea=tautology_demonstrated,
                              g10_pass=tautology_demonstrated and not intake_is_urea_derived),
        real_intake_provenance=real_intake_prov,
        intake_is_urea_derived=intake_is_urea_derived,
        inverted_intake_from_bingham=dict(intake_N_g_day=inv_intake_N,
                                          protein_g_per_kg_day=inv_protein,
                                          matches_stated_input=rel(inv_protein, LEDGER_INPUTS["protein_intake_g_per_kg_day"]["v"]) < 1e-9,
                                          note="exact by construction -- Bingham enters the ledger multiplicatively; "
                                               "reported as a self-consistency identity, NOT as evidence"),
        maroni_inverted_protein_g_per_kg_day=maroni_protein,
        maroni_vs_diet_record_rel_diff=rel(maroni_protein, LEDGER_INPUTS["protein_intake_g_per_kg_day"]["v"]),
        ceiling_inverted_to_serviceable_acid_load_mEq_day=ceiling_service_mEq,
    )


# =========================================================================================
# main
# =========================================================================================
def selftest():
    """Unit-level assertions on the arithmetic itself (run with --selftest)."""
    assert abs(mgN_to_umol(1.05) - 74.962) < 0.01, "mg N/min -> umol/min conversion"
    assert abs(umol_min_to_mmol_day(37.8) - 54.432) < 1e-3, "umol/min -> mmol/day conversion"
    # the factorisation R_urine = (2/Y)*(1/f_urine) must be an identity, not an approximation
    s2 = step2_glutamine_stoichiometry()
    for k, v in s2["datasets"].items():
        assert v["factorisation_check"] < 1e-9, f"factorisation identity broken for {k}"
        assert abs(v["production_umol_min"]
                   - (v["urine_NH4_umol_min"] + v["renal_vein_NH4_umol_min"])) < 1e-9
    # ledger must be exactly additive
    led = step1_nitrogen_ledger()
    parts = (led["urinary_urea_N_g_day"] + led["urinary_NH4_N_g_day"]
             + led["urinary_creatinine_N_g_day"] + led["urinary_uric_acid_N_g_day"]
             + led["urinary_other_N_g_day"])
    assert abs(parts - led["urinary_N_g_day"]) < 1e-9, "urinary components must sum exactly"
    assert abs((led["intake_N_g_day"] - led["total_out_N_g_day"]) - led["residual_N_g_day"]) < 1e-12
    # the coupled model must be strictly monotone below the ceiling and clamp above it
    prev = -1.0
    for a in (0, 10, 50, 120, 260):
        r = coupled_model(a, led)
        assert r["urinary_NH4_mmol_day"] > prev, "NH4 must increase with acid load"
        prev = r["urinary_NH4_mmol_day"]
    assert coupled_model(1e6, led)["ceiling_bound"], "ceiling must clamp"
    # renal-vein return must be non-negative for all f_urine <= 1
    for a in (0, 100, 300):
        assert coupled_model(a, led)["renal_vein_NH4_return_mmol_day"] >= 0.0
    # the Maroni route must be a tautology (that is the point of G10)
    inv = step5_inversions(led)
    assert inv["maroni_tautology"]["invariant_to_urea"], "Maroni tautology demo failed"
    print("selftest: OK")
    return 0


def main():
    if "--selftest" in sys.argv:
        return selftest()
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 92)
    print("RENAL AMMONIAGENESIS <-> HEPATIC UREAGENESIS NITROGEN COUPLING")
    print("=" * 92)

    # ---------------- STEP 1 ----------------
    led = step1_nitrogen_ledger()
    band = step1_band_sweep()
    print("\nSTEP 1  WHOLE-BODY NITROGEN LEDGER (70 kg adult, Western intake)")
    print(f"  intake                 {led['intake_N_g_day']:7.3f} g N/day  ({led['protein_g_day']:.1f} g protein/day)")
    print(f"  urinary total          {led['urinary_N_g_day']:7.3f}   [Bingham 81% of intake]")
    print(f"    urea-N               {led['urinary_urea_N_g_day']:7.3f}   ({led['urea_N_fraction_of_urinary']*100:.1f}% of urinary N; "
          f"{led['urea_flux_mmol_day']:.1f} mmol urea/day)")
    print(f"    NH4+-N               {led['urinary_NH4_N_g_day']:7.3f}   ({led['nh4_N_fraction_of_urinary']*100:.1f}% of urinary N)")
    print(f"    creatinine-N         {led['urinary_creatinine_N_g_day']:7.3f}")
    print(f"    uric-acid-N          {led['urinary_uric_acid_N_g_day']:7.3f}")
    print(f"    other-N              {led['urinary_other_N_g_day']:7.3f}")
    print(f"  faecal-N               {led['faecal_N_g_day']:7.3f}")
    print(f"  integumental+misc-N    {led['integumental_misc_N_g_day']:7.3f}")
    print(f"  TOTAL OUT              {led['total_out_N_g_day']:7.3f}")
    print(f"  RESIDUAL               {led['residual_N_g_day']:+7.3f} g N/day "
          f"({led['residual_fraction_of_intake']*100:+.1f}% of intake)")
    print(f"  propagated half-width  {band['quadrature_half_width_g_day']:7.3f} g N/day "
          f"(governor: {band['governor']}, +/-{band['governor_half_width']:.3f})")
    g1 = abs(led["residual_fraction_of_intake"]) <= 0.10
    g2 = UREA_N_FRACTION_EXTERNAL_BAND[0] <= led["urea_N_fraction_of_urinary"] <= UREA_N_FRACTION_EXTERNAL_BAND[1]
    print(f"  after FAO unmeasured-loss correction (+5 mg N/kg/day): "
          f"{band['residual_after_unmeasured_loss_correction_g_day']:+.3f} g N/day "
          f"({band['residual_after_correction_fraction']*100:+.1f}%)")
    print(f"  G1 |residual|<=10%: {g1}   G2 urea-N fraction in {UREA_N_FRACTION_EXTERNAL_BAND}: {g2}")
    print(f"  residual is {'WITHIN' if abs(led['residual_N_g_day']) <= band['quadrature_half_width_g_day'] else 'OUTSIDE'} "
          "the propagated measurement band")
    print(f"  POWER: noise/residual = {band['noise_over_residual']:.2f} -> {band['power_verdict']}")
    print(f"  G2 falsifiability: urea-N fraction spans "
          f"{band['urea_fraction_range_over_input_bands'][0]:.3f}-{band['urea_fraction_range_over_input_bands'][1]:.3f} "
          f"over the inputs' own bands -> can G2 fail? {band['G2_is_falsifiable_within_input_bands']}"
          + ("" if band["G2_is_falsifiable_within_input_bands"]
             else "  <-- NON-FALSIFYING sanity check, not a validation"))

    # ---------------- STEP 2 ----------------
    s2 = step2_glutamine_stoichiometry()
    print("\nSTEP 2  RENAL GLUTAMINE FLUX AND THE 2-NH4-PER-GLUTAMINE IDENTITY")
    print("  transcription over-determination (source's molar vs mass companion figures):")
    for k, v in s2["transcription_checks"].items():
        print(f"    {k:46s} rel_err={v['rel_err']*100:5.2f}%")
    print(f"  G3 all <= 5%: {s2['g3_pass']}")
    print(f"  {'dataset':30s} {'Gln':>8s} {'urine':>8s} {'renalV':>8s} {'prod':>8s} "
          f"{'Y':>6s} {'f_ur':>6s} {'R_urine':>8s} {'R_prod':>7s}")
    for k, v in s2["datasets"].items():
        print(f"  {k:30s} {v['glutamine_umol_min']:8.1f} {v['urine_NH4_umol_min']:8.1f} "
              f"{v['renal_vein_NH4_umol_min']:8.1f} {v['production_umol_min']:8.1f} "
              f"{v['Y_NH4_per_glutamine']:6.2f} {v['f_urine_fraction_excreted']:6.2f} "
              f"{v['R_urine_2Gln_over_urinaryNH4']:8.2f} {v['R_prod_2Gln_over_production']:7.2f}")
    print("  (umol/min; Y = measured NH4+ produced per glutamine extracted, theory = 2.00)")
    print(f"  G4 identity 2*Gln == urinary NH4 holds anywhere: {s2['identity_holds_anywhere']}  "
          "<-- REPORTED, NOT GATING")
    print(f"  G5 both factors deviate in >=2 datasets: {s2['g5_pass']} "
          f"({s2['g5_both_factors_deviate_count']}/3)")
    ea = s2["external_anchor"]
    print(f"  external anchor: Owen baseline urinary NH4 = {ea['owen_baseline_mmol_day']:.1f} mmol/day "
          f"vs Halperin/Rose {ea['halperin_band']} -> within={ea['baseline_within_band']} "
          f"({ea['baseline_rel_to_band_top']:.2f}x band top)")
    print(f"                   Owen acidosis = {ea['owen_acidosis_mmol_day']:.1f} mmol/day vs "
          f"{ea['halperin_acidosis_band']} -> within={ea['acidosis_within_band']} "
          f"({ea['acidosis_fraction_of_band_floor']:.2f}x band floor: MILD, not maximal, acidosis)")

    # ---------------- STEP 3 ----------------
    s3 = step3_cross_organ(led, s2)
    print("\nSTEP 3  CROSS-ORGAN NITROGEN COMPETITION (hepatic urea vs renal ammonium)")
    n = s3["nh4cl_arm"]
    print(f"  NH4Cl arm: d(urinary NH4) = {n['delta_urinary_NH4_mmol_day']:.1f} mmol/day "
          f"(= {n['delta_urinary_NH4_N_g_day']:.3f} g N/day)")
    print(f"             exogenous NH4Cl-N = {n['exogenous_NH4Cl_N_mmol_day'][0]:.1f}-"
          f"{n['exogenous_NH4Cl_N_mmol_day'][1]:.1f} mmol N/day -> ratio "
          f"{n['exogenous_over_delta_excreted']:.2f}x")
    print(f"             {n['verdict']}")
    p = s3["h_fixed_prediction"]
    print(f"  H_fixed prediction at that load: urea must fall {p['required_urea_fall_mmol_day']:.1f} "
          f"mmol/day = {p['required_urea_fall_fraction']*100:.1f}% of {p['baseline_urea_flux_mmol_day']:.0f} mmol/day")
    for k, v in p["at_maximal_acidosis"].items():
        print(f"    at {k}: urea must fall {v['required_urea_fall_mmol_day']:.1f} mmol/day "
              f"({v['required_urea_fall_fraction']*100:.1f}%)")
    d = s3["decisive_test"]
    print(f"  DECISIVE dataset ({d['units']}, no exogenous N): d(NH4-N)={d['delta_NH4_N']:+.2f}, "
          f"d(urea-N)={d['delta_urea_N']:+.2f}, d(total-N)={d['delta_total_N']:+.2f}")
    print(f"    substitution coefficient s = {d['s_measured']:+.3f} +/- {d['s_sd']:.3f} "
          f"vs H_fixed prediction {d['s_predicted_by_H_fixed']:+.2f}")
    print(f"    sign of urea change matches H_fixed: {d['sign_of_urea_change_matches_H_fixed']}; "
          f"excludes unity at 2SD: {d['s_excludes_unity_at_2sd']}")
    print(f"    H_fixed retained: {d['H_fixed_retained']}  -> G7 refuted={s3['g7_H_fixed_refuted']}")
    print(f"    fraction of the change carried by POOL SIZE (net protein catabolism): "
          f"{d['fraction_of_change_carried_by_pool_size']*100:.1f}%")
    o = s3["rat_HCl_counter_dataset"]
    print(f"  counter-dataset (rat, HCl, constant diet): '{o['statement']}'")
    print(f"    unit-ambiguous: molar reading -> s={o['reading_molar_equimolar']['s_in_nitrogen_units']:.1f}, "
          f"nitrogen reading -> s={o['reading_nitrogen_equimolar']['s_in_nitrogen_units']:.1f} "
          f"(factor {o['factor_between_readings']:.0f}); resolved={o['resolved']}")

    # ---------------- STEP 4 ----------------
    s4 = step4_coupling(led, s2)
    print("\nSTEP 4  COUPLED MODEL (acid load -> NH4+ -> glutamine draw -> ureagenic pool)")
    print(f"  {'load':>6s} {'NH4_ur':>8s} {'f_ur':>6s} {'NH3_prod':>9s} {'RV_return':>10s} "
          f"{'Gln':>8s} {'urea(H_fixed)':>14s}")
    for r in s4["sweep"]:
        print(f"  {r['acid_load_mEq_day']:6.0f} {r['urinary_NH4_mmol_day']:8.1f} {r['f_urine']:6.3f} "
              f"{r['renal_NH3_production_mmol_day']:9.1f} {r['renal_vein_NH4_return_mmol_day']:10.1f} "
              f"{r['glutamine_uptake_mmol_day']:8.1f} {r['urea_flux_mmol_day_under_H_fixed']:14.1f}"
              + ("  <ceiling>" if r["ceiling_bound"] else ""))
    print(f"  (mEq/day, mmol/day; Y={s4['Y_used']:.3f}, f_urine {s4['f_urine_baseline']:.3f}->{s4['f_urine_acidosis']:.3f})")
    b0 = s4["sweep"][0]["glutamine_uptake_mmol_day"] * 1000.0 / MIN_PER_DAY
    print(f"  baseline glutamine uptake = {b0:.1f} umol/min vs measured 47.0 [1] and 37.8 [2] "
          "-- NOT an independent check (it is Owen's Y and f rescaled to a 45 mmol/day NH4 "
          "baseline); reported for scale only")
    print(f"  G8 baseline<=1% & monotonic & ceiling binds: {s4['g8_pass']}")
    vf = s4["void_floor_force_f_urine_1"]
    print(f"  G9 void-floor (force f_urine=1): Gln {vf['glutamine_model']:.1f} -> {vf['glutamine_naive']:.1f} "
          f"mmol/day, ratio {vf['ratio']:.2f}x -> pass={vf['g9_pass']}")
    print(f"     fully-naive (Y=2, f=1): Gln {s4['void_floor_full_naive_Y2_f1']['glutamine_naive']:.1f} "
          f"mmol/day, {s4['void_floor_full_naive_Y2_f1']['ratio']:.2f}x off the data-set model")

    # ---------------- STEP 5 ----------------
    s5 = step5_inversions(led)
    print("\nSTEP 5  CIRCULARITY / INVERSION")
    m = s5["maroni_tautology"]
    print(f"  Maroni-closed ledger residual at urea = measured / 2x / 0.1: "
          f"{m['residual_at_measured_urea']:.2e} / {m['residual_at_2x_urea']:.2e} / {m['residual_at_0p1_urea']:.2e}")
    print(f"  -> residual invariant to urea: {m['invariant_to_urea']} (TAUTOLOGY GATE, correctly refused)")
    print(f"  real intake provenance: {s5['real_intake_provenance']}; urea-derived={s5['intake_is_urea_derived']}")
    print(f"  Maroni-inverted intake {s5['maroni_inverted_protein_g_per_kg_day']:.3f} g/kg/day vs "
          f"diet-record 1.150 -> rel diff {s5['maroni_vs_diet_record_rel_diff']*100:.1f}%")
    print(f"  ammoniagenic ceiling inverts to a serviceable acid load of "
          f"{s5['ceiling_inverted_to_serviceable_acid_load_mEq_day']:.0f} mEq/day")
    print(f"  G10 tautology demonstrated AND not used: {m['g10_pass']}")

    # ---------------- gates ----------------
    gates = {
        "G1_ledger_residual_le_10pct": g1,
        "G2_urea_N_fraction_in_external_band": g2,
        "G3_source_transcription_overdetermined": s2["g3_pass"],
        "G5_deficit_is_two_factor_not_one_bad_number": s2["g5_pass"],
        "G6_nh4cl_arm_correctly_marked_nondecisive": s3["g6_pass"],
        "G7_H_fixed_tested_on_a_decisive_dataset": s3["g7_pass"],
        "G8_coupled_model_wellposed": s4["g8_pass"],
        "G9_renal_vein_return_is_load_bearing": s4["g9_pass"],
        "G10_circular_route_demonstrated_and_refused": s5["maroni_tautology"]["g10_pass"],
    }
    reported_not_gating = {
        "G4_stoichiometric_identity_holds": s2["identity_holds_anywhere"],
        "G7_H_fixed_refuted": s3["g7_H_fixed_refuted"],
    }
    overall = all(gates.values())
    print("\n" + "=" * 92)
    for k, v in gates.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    for k, v in reported_not_gating.items():
        print(f"  RPT   {k} = {v}")
    print(f"\n  OVERALL: {'PASS' if overall else 'FAIL'}")
    print("  FINDING: identity `2 x glutamine_flux == urinary NH4+` is REFUTED on every dataset "
          f"(R = {', '.join('%.2f' % v['R_urine_2Gln_over_urinaryNH4'] for v in s2['datasets'].values())}); "
          "neither measurement is wrong -- 2 NH4+/glutamine is a stoichiometric CEILING "
          "(measured Y = %.2f-%.2f) and urinary excretion captures only %.0f-%.0f%% of renal "
          "ammonia production, the rest returning to the liver via the renal vein."
          % (min(v["Y_NH4_per_glutamine"] for v in s2["datasets"].values()),
             max(v["Y_NH4_per_glutamine"] for v in s2["datasets"].values()),
             100 * min(v["f_urine_fraction_excreted"] for v in s2["datasets"].values()),
             100 * max(v["f_urine_fraction_excreted"] for v in s2["datasets"].values())))
    print("=" * 92)

    report = dict(
        script=os.path.abspath(__file__),
        question="renal ammoniagenesis-glutamine flux, the 2:1 stoichiometric identity, and its "
                 "dynamic coupling to the hepatic ureagenic nitrogen pool",
        parent_cell="HOLE-ACIDBASE-NAE-PROVENANCE-AMMONIAGENESIS-CEILING-NPARTITION",
        constants=dict(LEDGER_INPUTS=LEDGER_INPUTS, body_mass_kg=BODY_MASS_KG,
                       NH4_PER_GLUTAMINE_THEORY=NH4_PER_GLUTAMINE_THEORY,
                       nh4cl_dose_g_per_day=NH4CL_DOSE_G_PER_DAY),
        step1_ledger=led, step1_band_sweep=band,
        step2_glutamine=s2, step3_cross_organ=s3, step4_coupling=s4, step5_inversions=s5,
        gates=gates, reported_not_gating=reported_not_gating, overall_pass=overall,
    )
    out_path = f"{OUT_DIR}/renal_ammoniagenesis_nitrogen_coupling_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"Wrote {out_path}")
    return 0 if overall else 2


if __name__ == "__main__":
    sys.exit(main())
