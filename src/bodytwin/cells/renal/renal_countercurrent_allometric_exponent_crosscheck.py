"""renal_countercurrent_allometric_exponent_crosscheck.py -- cross-species allometric exponents.

Resolves the allometric-exponent sub-computation of MODEL-KIDNEY-COUNTERCURRENT-MULTIPLIER
(electrolyte/fluid-balance cluster): does the node's cited Beuchat 1990/1996 cross-species
power-law exponents reproduce its claimed "RMT and Uosm mass-exponents match each other at ratio
0.90, while the unrelated Kleiber whole-organism metabolic exponent (-0.25) diverges from them by
~2.6x" -- a geometry-over-metabolism triangulation from two decorrelated data families (renal
concentrating-capacity scaling vs whole-body metabolic scaling)?

DISTINCT from MODEL-RENAL-COUNTERCURRENT-MULTIPLIER (a sibling narrative-only node with an
overlapping topic stem, id NOT built here): that sibling's claim is about the Greenwald
metabolic-CAPACITY-density hypothesis (mitochondrial infoldings per unit volume) crosschecked
against a from-scratch Na-K-ATPase electrochemical calculation -- a completely different
computation (thermodynamic driving force vs kinetic capacity) from this script's allometric-
exponent-ratio arithmetic. Both nodes are left with couples_to unresolved between them; this script
only claims MODEL-KIDNEY-COUNTERCURRENT-MULTIPLIER.

SCOPE DECISION (avoiding a too-minimal-rebuild false refutation): the node's honest_gaps
explicitly flags that (a) the absolute cross-species RMT->Uosm point predictions (e.g. "human
0.72x, kangaroo rat 0.31x of measured") depend on a reference/normalization convention the claim
text does not pin down numerically (which species/mass anchors the power-law's leading
coefficient), and (b) R^2=0.59 means 41% of variance is unexplained -- reconstructing those exact
per-species point-prediction ratios without the paper's fitted leading coefficients would risk
inventing a proxy. This script therefore gates ONLY the two fully-specified, coefficient-free
ratios the node states directly: exponent-ratio (0.90) and Kleiber-divergence (2.6x), both of which
depend on the exponents ALONE, not on any leading-coefficient convention. The absolute cross-
species point predictions are computed and reported as informational, not gated.

STATED NUMBERS (verbatim, Beuchat 1990 AJP258:R298 / Beuchat 1996 AJP271:R157):
  RMT = 5.408 * M^-0.108   (relative medullary thickness vs body mass)
  Uosm = 2564 * M^-0.097   (max urine osmolality vs body mass)
  Kleiber whole-organism metabolic exponent: -0.25

PRE-REGISTERED GATES:
  G1 exponent_ratio = |{-0.097}| / |{-0.108}| must reproduce the claimed 0.90 within 2%.
  G2 kleiber_divergence = 0.25 / |{-0.097}| must reproduce the claimed 2.6x within 2%.
  G3 VOID FLOOR: swap which exponent belongs to which quantity (RMT<->Uosm exponents exchanged)
     -> exponent_ratio is a symmetric quotient so this ALONE is a weak floor; the true floor is
     recomputing kleiber_divergence against the WRONG (swapped-in) Uosm exponent (-0.108 instead
     of -0.097) -> must NOT reproduce 2.6x within 2% (must differ by >5%), proving G2 depends on
     using the correct, non-interchangeable exponent, not "any exponent near 0.1 gives ~2.6x".

Run: python3 renal_countercurrent_allometric_exponent_crosscheck.py
Reads: nothing.
Writes: OUT_ROOT/renal_countercurrent_allometric_exponent_crosscheck/renal_countercurrent_allometric_exponent_crosscheck.json
Gates: G1-G3 below decide.
"""
import json
import math
import os

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT = os.path.join(OUT_ROOT, "renal_countercurrent_allometric_exponent_crosscheck",
                   "renal_countercurrent_allometric_exponent_crosscheck.json")

RMT_EXPONENT = -0.108
UOSM_EXPONENT = -0.097
KLEIBER_EXPONENT = -0.25

CLAIMED_EXPONENT_RATIO = 0.90
CLAIMED_KLEIBER_DIVERGENCE = 2.6

# informational cross-species check (coefficient-dependent, NOT gated -- see scope_note)
RMT_COEF, UOSM_COEF = 5.408, 2564.0
HUMAN_MASS_KG = 70.0
KANGAROO_RAT_MASS_KG = 0.042
HUMAN_RMT_MEASURED, KANGAROO_RMT_MEASURED = 4.5, 8.5
HUMAN_UOSM_MEASURED, KANGAROO_UOSM_MEASURED_MID = 1200.0, 5750.0


def power_law(coef, mass_kg, exponent):
    return coef * (mass_kg ** exponent)


def main():
    gates = {}

    exponent_ratio = abs(UOSM_EXPONENT) / abs(RMT_EXPONENT)
    g1 = abs(exponent_ratio - CLAIMED_EXPONENT_RATIO) / CLAIMED_EXPONENT_RATIO < 0.02
    gates["G1_exponent_ratio_reproduces_0.90"] = bool(g1)

    kleiber_divergence = abs(KLEIBER_EXPONENT) / abs(UOSM_EXPONENT)
    g2 = abs(kleiber_divergence - CLAIMED_KLEIBER_DIVERGENCE) / CLAIMED_KLEIBER_DIVERGENCE < 0.02
    gates["G2_kleiber_divergence_reproduces_2.6x"] = bool(g2)

    # G3 void floor: use the WRONG (RMT) exponent in the Kleiber-divergence formula
    kleiber_divergence_void = abs(KLEIBER_EXPONENT) / abs(RMT_EXPONENT)
    void_fails = abs(kleiber_divergence_void - CLAIMED_KLEIBER_DIVERGENCE) / CLAIMED_KLEIBER_DIVERGENCE >= 0.05
    gates["G3_void_floor_wrong_exponent_fails_2.6x"] = bool(void_fails)

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "CONFIRMED" if all(gates.values()) else "DISAGREE"

    # informational-only cross-species prediction (NOT gated, coefficient/normalization dependent)
    rmt_human_pred = power_law(RMT_COEF, HUMAN_MASS_KG, RMT_EXPONENT)
    rmt_kangaroo_pred = power_law(RMT_COEF, KANGAROO_RAT_MASS_KG, RMT_EXPONENT)
    uosm_human_pred = power_law(UOSM_COEF, HUMAN_MASS_KG, UOSM_EXPONENT)
    uosm_kangaroo_pred = power_law(UOSM_COEF, KANGAROO_RAT_MASS_KG, UOSM_EXPONENT)

    result = {
        "node_id": "MODEL-KIDNEY-COUNTERCURRENT-MULTIPLIER",
        "exponent_ratio_computed": exponent_ratio,
        "kleiber_divergence_computed": kleiber_divergence,
        "kleiber_divergence_void_floor_wrong_exponent": kleiber_divergence_void,
        "gates": gates,
        "verdict": verdict,
        "informational_not_gated": {
            "rmt_predicted_human": rmt_human_pred, "rmt_measured_human": HUMAN_RMT_MEASURED,
            "rmt_predicted_kangaroo_rat": rmt_kangaroo_pred, "rmt_measured_kangaroo_rat": KANGAROO_RMT_MEASURED,
            "uosm_predicted_human_mOsm": uosm_human_pred, "uosm_measured_human_mOsm": HUMAN_UOSM_MEASURED,
            "uosm_predicted_kangaroo_rat_mOsm": uosm_kangaroo_pred,
            "uosm_measured_kangaroo_rat_mOsm_mid": KANGAROO_UOSM_MEASURED_MID,
            "note": ("absolute point predictions from the bare power laws using the paper's "
                     "leading coefficients (fitted across the FULL cross-species dataset, not just "
                     "these 2 points) -- expected to be order-of-magnitude-consistent but not "
                     "point-exact per the node's honest_gaps (R^2=0.59, n=2 species used here "
                     "is illustrative only); NOT pre-registered as a gate."),
        },
        "scope_note": ("The absolute per-species Uosm/RMT prediction ratios cited in the node's "
                        "own honest_gaps (e.g. 'human 0.72x, kangaroo rat 0.31x of measured') are "
                        "NOT gated here -- reconstructing them exactly requires a reference-point "
                        "convention the claim text does not specify; only the coefficient-free "
                        "exponent-ratio and Kleiber-divergence numbers, which the claim states "
                        "directly as this node's decorrelated-anchor triangulation, are gated."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict, "exponent_ratio": exponent_ratio,
                       "kleiber_divergence": kleiber_divergence}, indent=2))
    return result


if __name__ == "__main__":
    main()
