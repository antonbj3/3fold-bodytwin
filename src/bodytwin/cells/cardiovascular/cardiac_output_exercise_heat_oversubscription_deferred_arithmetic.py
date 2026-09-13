"""
Combined maximal exercise + heat stress: cardiac-output budget oversubscription and how much
of it published flow-reallocation actually closes.

At maximal exercise the muscle flow fraction of CO_max plus the passive-heat-stress skin-flow
ceiling sum to MORE than CO_max (a face-value oversubscription of 17.0-25.5 percentage points).
This cell computes the deferred quantity: how much skin/muscle-flow attenuation would be
REQUIRED to close the gap, and whether external, decorrelated exercise-in-heat literature
supports that much attenuation.

PRE-REGISTERED GATE (stated BEFORE computing the literature cross-check below):
  C  = the literature-reported flow reallocation under combined maximal
       exercise+heat (Gonzalez-Alonso 2008 J Physiol review: total CO decline
       ~1.2 L/min under heat+dehydration, of which "two-thirds" is attributed to
       combined leg+skin flow reduction) closes >=50% of the required
       flow-reduction gap computed from the cited numbers.
  not-C = it closes <50%: the face-value oversubscription remains a real,
       largely-unresolved structural gap even after crediting the one concrete
       literature-supported compensation mechanism found.
  Falsifier for C: required_reduction_L_min * 0.5 <= literature_credit_L_min

ADVERSARY forced: the naive resolution is "skin flow must just be lower during
exercise than during passive heat" -- but Rowell 1986 (restated in multiple
secondary sources) reports SkBF "as high as 7 L/min" IS SUSTAINED during
exercise-heat stress, not merely a passive-heat artifact -- so the adversary "the
6-8 L/min ceiling doesn't really apply once you're exercising" does NOT trivially
resolve this; forced to its strongest form it must specify how MUCH of a CO decline
is real and independently measured (~1.2 L/min, 2/3 muscle+skin) -- used here as the
anchor, not asserted from recall.

This is a single-scenario budget-arithmetic problem, not a permutation/shuffle
setting -- there is no population to scramble. The substitute for a void-floor here
is an explicit uncertainty SWEEP over every cited band (muscle fraction 88.1-94.0%,
skin ceiling 6-8 L/min) to confirm the <50%-closure conclusion is not an artifact of
picking convenient point estimates.

Reads:  nothing (all constants embedded).
Writes: results.json under the cell output directory.
Gate:   closure fraction >= 0.50 in the most favorable cited band combination.
"""
import json
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "cardiac_output_exercise_heat_oversubscription_deferred_arithmetic")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Inputs, all cited (not invented):
# ---------------------------------------------------------------------------
CO_MAX_L_MIN = 23.6208  # Higginbotham 1986-anchored, from the cardiac_output_geometric cell
MUSCLE_FRAC_CENTRAL = 0.91647
MUSCLE_FRAC_BAND = (0.88067, 0.94033)  # sensitivity band, cited
NON_MUSCLE_BUDGET_L_MIN = 1.973  # = CO_MAX*(1-MUSCLE_FRAC_CENTRAL), cited & cross-checked below
SKIN_CEILING_BAND_L_MIN = (6.0, 8.0)  # Rowell1974/Charkoudian2003 PMID12744548, passive heat-stress

# cross-check the cited non-muscle budget against direct arithmetic (must reproduce)
non_muscle_recomputed = CO_MAX_L_MIN * (1.0 - MUSCLE_FRAC_CENTRAL)
non_muscle_selfcheck_relerr_pct = abs(non_muscle_recomputed - NON_MUSCLE_BUDGET_L_MIN) / NON_MUSCLE_BUDGET_L_MIN * 100.0

# ---------------------------------------------------------------------------
# STEP 1: reproduce the cited face-value oversubscription (17.0-25.5pp)
# ---------------------------------------------------------------------------
def joint_total_L_min(muscle_frac, skin_L_min):
    return CO_MAX_L_MIN * muscle_frac + skin_L_min

# cited 17.0-25.5pp band holds MUSCLE_FRAC at its central point estimate and
# varies ONLY the skin ceiling (verified: pairing both bands independently
# widens it to 13.5-27.9pp, which does NOT match the cited figure -- this
# mismatch was caught by running it both ways, not assumed)
recon_lo = joint_total_L_min(MUSCLE_FRAC_CENTRAL, SKIN_CEILING_BAND_L_MIN[0])
recon_hi = joint_total_L_min(MUSCLE_FRAC_CENTRAL, SKIN_CEILING_BAND_L_MIN[1])
over_pp_lo = (recon_lo - CO_MAX_L_MIN) / CO_MAX_L_MIN * 100.0
over_pp_hi = (recon_hi - CO_MAX_L_MIN) / CO_MAX_L_MIN * 100.0
cited_reproduced = (16.5 <= over_pp_lo <= 17.5) and (25.0 <= over_pp_hi <= 26.0)

# ---------------------------------------------------------------------------
# STEP 2: THE DEFERRED QUANTITY -- required skin-flow value/reduction to fit CO_max,
# swept over the full cited muscle-frac and skin-ceiling bands (not just central pt).
# ---------------------------------------------------------------------------
grid_muscle = [MUSCLE_FRAC_BAND[0], MUSCLE_FRAC_CENTRAL, MUSCLE_FRAC_BAND[1]]
grid_skin_ceiling = [SKIN_CEILING_BAND_L_MIN[0], 7.0, SKIN_CEILING_BAND_L_MIN[1]]

sweep_rows = []
for mf in grid_muscle:
    muscle_flow = CO_MAX_L_MIN * mf
    skin_allowed_L_min = CO_MAX_L_MIN - muscle_flow  # what skin flow the budget can afford
    for sc in grid_skin_ceiling:
        required_reduction_L_min = sc - skin_allowed_L_min
        required_reduction_frac = required_reduction_L_min / sc
        sweep_rows.append({
            "muscle_frac": mf, "muscle_flow_L_min": round(muscle_flow, 4),
            "skin_ceiling_L_min": sc, "skin_allowed_L_min": round(skin_allowed_L_min, 4),
            "required_reduction_L_min": round(required_reduction_L_min, 4),
            "required_reduction_frac_of_ceiling": round(required_reduction_frac, 4),
        })

required_reduction_min = min(r["required_reduction_L_min"] for r in sweep_rows)
required_reduction_max = max(r["required_reduction_L_min"] for r in sweep_rows)

# ---------------------------------------------------------------------------
# STEP 3: external, decorrelated literature anchor -- not the reference body's
# numbers, not fetched to make this pass.
#   Gonzalez-Alonso 2008 J Physiol review ("The cardiovascular challenge of
#   exercising in the heat"): CO decline under heat+dehydration ~1.2 L/min;
#   "lowering in exercising leg blood flow and skin blood flow account for
#   two-thirds of the decline in cardiac output."
# ---------------------------------------------------------------------------
LIT_CO_DECLINE_L_MIN = 1.2
LIT_FRACTION_FROM_MUSCLE_SKIN = 2.0 / 3.0
literature_credit_L_min = LIT_CO_DECLINE_L_MIN * LIT_FRACTION_FROM_MUSCLE_SKIN  # = 0.8 L/min

closure_fraction_worst = literature_credit_L_min / required_reduction_max
closure_fraction_best = literature_credit_L_min / required_reduction_min

GATE_THRESHOLD = 0.50
gate_C_holds_at_best_case = closure_fraction_best >= GATE_THRESHOLD
gate_C_holds_at_worst_case = closure_fraction_worst >= GATE_THRESHOLD
# Robust-negative test (symmetric-QC "adversary must fall across the instance
# space"): C is only accepted if it holds even in the LEAST-favorable cited
# band combination; here we check whether it holds even in the MOST favorable,
# which is the fair steelman for accepting C.
verdict = "C-HOLDS-RESOLVED" if gate_C_holds_at_best_case else "NOT-C-STRUCTURAL-GAP-REMAINS"

residual_gap_after_credit_min_L_min = required_reduction_min - literature_credit_L_min
residual_gap_after_credit_max_L_min = required_reduction_max - literature_credit_L_min
residual_pp_of_CO_max_range = (
    round(residual_gap_after_credit_min_L_min / CO_MAX_L_MIN * 100.0, 2),
    round(residual_gap_after_credit_max_L_min / CO_MAX_L_MIN * 100.0, 2),
)

OUT = {
    "cell": "cardiac-output combined exercise-heat oversubscription",
    "pre_registered_gate": {
        "C": "literature-supported flow-reallocation credit closes >=50% of the required reduction",
        "threshold": GATE_THRESHOLD,
    },
    "step1_reproduce_cited_facevalue_oversubscription": {
        "non_muscle_budget_selfcheck_relerr_pct": round(non_muscle_selfcheck_relerr_pct, 4),
        "over_pp_band_recomputed_muscle_fixed_skin_varies": (round(over_pp_lo, 2), round(over_pp_hi, 2)),
        "cited_band_17.0_to_25.5pp_reproduced": cited_reproduced,
        "wider_band_if_muscle_AND_skin_both_vary_independently": (
            round((joint_total_L_min(MUSCLE_FRAC_BAND[0], SKIN_CEILING_BAND_L_MIN[0]) - CO_MAX_L_MIN) / CO_MAX_L_MIN * 100.0, 2),
            round((joint_total_L_min(MUSCLE_FRAC_BAND[1], SKIN_CEILING_BAND_L_MIN[1]) - CO_MAX_L_MIN) / CO_MAX_L_MIN * 100.0, 2),
        ),
    },
    "step2_deferred_quantity_required_skin_reduction": {
        "sweep_rows": sweep_rows,
        "required_reduction_L_min_range": (round(required_reduction_min, 3), round(required_reduction_max, 3)),
        "required_reduction_frac_of_ceiling_range": (
            round(min(r["required_reduction_frac_of_ceiling"] for r in sweep_rows), 3),
            round(max(r["required_reduction_frac_of_ceiling"] for r in sweep_rows), 3),
        ),
    },
    "step3_external_literature_anchor": {
        "source": "Gonzalez-Alonso 2008 J Physiol 586(1):45-53 'The cardiovascular challenge of exercising in the heat' (decorrelated from the CO/flow-redistribution primaries used here)",
        "lit_CO_decline_L_min": LIT_CO_DECLINE_L_MIN,
        "lit_fraction_attributed_to_muscle_plus_skin": LIT_FRACTION_FROM_MUSCLE_SKIN,
        "literature_credit_L_min": round(literature_credit_L_min, 3),
    },
    "step4_gate_verdict": {
        "closure_fraction_best_case": round(closure_fraction_best, 4),
        "closure_fraction_worst_case": round(closure_fraction_worst, 4),
        "gate_C_holds_at_best_case": gate_C_holds_at_best_case,
        "gate_C_holds_at_worst_case": gate_C_holds_at_worst_case,
        "verdict": verdict,
        "residual_unresolved_gap_pp_of_CO_max_range": residual_pp_of_CO_max_range,
    },
    "honest_note": (
        "This is a single-scenario budget-arithmetic problem (no population to "
        "permute), so the void-floor requirement is satisfied here by an "
        "explicit uncertainty sweep over every cited band (muscle-frac x "
        "skin-ceiling), not a scramble. The literature anchor (Gonzalez-Alonso "
        "2008) was used specifically to test the claim that 'consistent with "
        "SOME resolution existing' "
        "-- the result QUANTIFIES that the one concrete literature-supported "
        "mechanism found covers only 12-25% of the required gap across the "
        "full cited band, i.e. the honest_gap's hedge ('some resolution exists') "
        "is directionally correct but the gap is NOT closed by it -- most of "
        "the oversubscription is a genuine unresolved structural property of "
        "this model (fixed CO_max shared across regimes that real exercise-heat "
        "physiology does not actually hold fixed)."
    ),
}

print(json.dumps(OUT, indent=2))
with open(os.path.join(OUT_DIR, "results.json"), "w") as fh:
    json.dump(OUT, fh, indent=2)
