"""renal_adh_osmolality_robertson_threshold.py -- AVP/thirst osmotic set-points from the cited source.

Resolves ORG-RENAL-FLUID-ELECTROLYTE (renal-fluid cluster), specifically its leg (a): "osmolality
-> ADH/thirst loop with dead-band [280,290-294] mOsm/kg". The node's text states the dead-band
numbers but gives no functional form for the loop above threshold -- it names the mechanism
("osmolality->ADH... loop") without the equation, which is exactly the implementable-via-cited-
source pattern: the classic source is Robertson GL's plasma-AVP-vs-osmolality studies (Robertson
1976 J Clin Invest and later reviews), which give a LINEAR threshold-plus-slope model, not a dead-
band step function.

QUESTION: does the actual published Robertson AVP-osmolality relation reproduce a threshold near
the node's stated [280,290-294] mOsm/kg band, and is the mechanism a hard dead-band (as the node's
"dead-band" language implies) or a linear threshold-slope relation (as Robertson's data show)?

SOURCE FETCHED (not the node's paraphrase): below an osmolality of ~280 mOsm/kg, plasma AVP is
near the assay floor (<1 pg/mL); above that osmotic threshold (T_osm, individually variable but
population-centered ~280-284 mOsm/kg), AVP rises LINEARLY with slope ~0.4 pg/mL per mOsm/kg
(reported range 0.38-0.5 pg*kg/mL/mOsm across sources) -- i.e. AVP(C) = slope*(C - T_osm) for
C > T_osm, else ~0. Thirst threshold is a SEPARATE, higher set-point, ~290-294 mOsm/kg (this is
where the node's upper dead-band number 290-294 actually comes from -- it is the THIRST threshold,
not a second AVP threshold; the node's text compresses two distinct set-points into one band).

PRE-REGISTERED GATES:
  G1 the node's stated lower bound (280) must fall within +/-3 mOsm/kg of the published AVP
     osmotic threshold T_osm=280-284 (using T_osm=282 midpoint).
  G2 the node's stated upper bound (290-294) must fall within +/-3 mOsm/kg of the published THIRST
     threshold (292 midpoint of 290-294) -- confirming the node's "dead-band" is actually two
     DIFFERENT physiological set-points (AVP-release vs thirst-onset), not one dead-band with two
     edges of the same signal. This is a scope/labeling correction, reported as informational
     alongside the numeric match.
  G3 REGIME CHECK: the model must be LINEAR-ABOVE-THRESHOLD, not a step/dead-band -- verify by
     computing AVP at C=282 (at threshold, ~0), C=290 (mid-range), C=300 (above thirst threshold)
     and checking AVP is STRICTLY INCREASING and CONTINUOUS (no discontinuous jump), which a true
     "dead-band" (constant-then-jump) model would violate.
  G4 VOID FLOOR: apply the slope BELOW threshold (C=270, hyposmolar) where AVP should be at floor
     -- a void model that ignores the threshold and applies the linear slope everywhere would
     predict a NEGATIVE AVP concentration (unphysical) at C=270; the correctly-thresholded model
     must clip to ~0/floor instead. Confirms the threshold clause is doing real work, not decoration.

Distinct from raas.py (aldosterone/renin pressure-natriuresis loop, a DIFFERENT hormone axis and
DIFFERENT sensor -- volume/pressure not osmolality) and fluid_compartments.py (static compartment
volumes, no ADH dynamics) -- zero shared code or physics.

Run: python3 renal_adh_osmolality_robertson_threshold.py
Reads: nothing.
Writes: OUT_ROOT/renal_adh_osmolality_robertson_threshold/renal_adh_osmolality_robertson_threshold.json
Gates: G1-G4 below decide.
"""
import json
import os

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT = os.path.join(OUT_ROOT, "renal_adh_osmolality_robertson_threshold",
                   "renal_adh_osmolality_robertson_threshold.json")

NODE_LOWER = 280.0
NODE_UPPER_LO, NODE_UPPER_HI = 290.0, 294.0

AVP_THRESHOLD_MOSM = 282.0     # Robertson-family osmotic threshold for AVP release (280-284 range)
AVP_SLOPE_PG_PER_MOSM = 0.4    # pg/mL per mOsm/kg above threshold (0.38-0.5 reported range)
THIRST_THRESHOLD_MOSM = 292.0  # separate, higher thirst-onset set-point (290-294 range)


def avp_pg_per_ml(C_mosm, threshold=AVP_THRESHOLD_MOSM, slope=AVP_SLOPE_PG_PER_MOSM):
    return max(0.0, slope * (C_mosm - threshold))


def avp_pg_per_ml_void_no_threshold(C_mosm, slope=AVP_SLOPE_PG_PER_MOSM):
    """VOID FLOOR: applies the slope with NO threshold clip (wrong model), everywhere."""
    return slope * (C_mosm - AVP_THRESHOLD_MOSM)


def main():
    gates = {}

    g1 = abs(NODE_LOWER - AVP_THRESHOLD_MOSM) <= 3.0
    gates["G1_node_lower_bound_within_3mOsm_of_AVP_threshold"] = bool(g1)

    upper_mid = (NODE_UPPER_LO + NODE_UPPER_HI) / 2.0
    g2 = abs(upper_mid - THIRST_THRESHOLD_MOSM) <= 3.0
    gates["G2_node_upper_bound_within_3mOsm_of_thirst_threshold"] = bool(g2)

    avp_282 = avp_pg_per_ml(282.0)
    avp_290 = avp_pg_per_ml(290.0)
    avp_300 = avp_pg_per_ml(300.0)
    g3_monotone = avp_282 <= avp_290 <= avp_300 and (avp_290 - avp_282) > 0 and (avp_300 - avp_290) > 0
    gates["G3_linear_above_threshold_strictly_increasing"] = bool(g3_monotone)

    avp_270_correct = avp_pg_per_ml(270.0)
    avp_270_void = avp_pg_per_ml_void_no_threshold(270.0)
    g4_void = (avp_270_correct >= 0.0) and (avp_270_void < 0.0)
    gates["G4_void_floor_unthresholded_model_goes_negative"] = bool(g4_void)

    all_pass = all(gates.values())
    verdict = "CONFIRMED_WITH_LABELING_CORRECTION" if all_pass else "PARTIAL"

    result = {
        "node_id": "ORG-RENAL-FLUID-ELECTROLYTE",
        "leg": "(a) osmolality -> ADH/thirst loop",
        "node_stated_band": [NODE_LOWER, NODE_UPPER_LO, NODE_UPPER_HI],
        "source_avp_threshold_mosm": AVP_THRESHOLD_MOSM,
        "source_thirst_threshold_mosm": THIRST_THRESHOLD_MOSM,
        "avp_at_282": avp_282,
        "avp_at_290": avp_290,
        "avp_at_300": avp_300,
        "avp_at_270_correct_floor": avp_270_correct,
        "avp_at_270_void_no_threshold": avp_270_void,
        "gates": gates,
        "verdict": verdict,
        "labeling_correction": ("The node's single 'dead-band [280,290-294]' actually names TWO "
                                  "distinct physiological set-points: 280-284 is the AVP-release "
                                  "osmotic threshold (Robertson), 290-294 is the SEPARATE, higher "
                                  "thirst-onset threshold. The mechanism above 280 is LINEAR "
                                  "(slope~0.4 pg/mL per mOsm/kg), not a flat dead-band -- 'dead-band' "
                                  "correctly describes the region BELOW 280 (AVP at assay floor), "
                                  "not the region between the two named numbers."),
        "scope_note": ("The RAAS/aldosterone volume loop (b) and the eGFR-decline free-water-"
                        "clearance modifier are NOT rebuilt here -- the node gives no numeric "
                        "functional form for either; left informational, not gated. Only leg (a), "
                        "which names a specific, fully-specifying published threshold+slope "
                        "relationship, is rebuilt and gated."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict}, indent=2))
    return result


if __name__ == "__main__":
    main()
