"""Deferred arithmetic for HOLE-WATER-OSMOTIC-LEDGER-CROSSCELL-CONTRADICTION.

Closes the daily water-output ledger (urine + skin + respiratory + sweat + stool) against the
2.5 L/day anchor and extends the node's single void-floor case (drop urine) to a full
LEAVE-ONE-OUT sweep over all 5 components.

Reads: OUT_ROOT/renal_filtration/renal_filtration_results.json (GFR).
Writes: OUT_ROOT/water_osmotic_ledger_leaveoneout_deferred_arithmetic/results.json
Gate: the 15% closure band, applied to the baseline sum and to every leave-one-out case.

Inputs, cited verbatim from the node's leg text:
  urine = 1.415 L/day   (independently re-derived below from
                          renal_filtration.py's GFR=122.8 mL/min output
                          x (1-0.992) reabsorption x 1440 min/day -- NOT
                          just copied, machine-checked to reproduce)
  skin  = 0.353 L/day   (skin_barrier_tewl.py)
  resp  = 0.325 L/day   (external constant, cited)
  sweat = 0.10  L/day   (eccrine_sweat_gland.py family, resting)
  stool = 0.15  L/day   (external constant, cited)
  anchor = 2.5 L/day    (StatPearls/Guyton total daily water output)
  tolerance = 15%       (the node's pre-registered PASS band)
"""
import json
import os

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "water_osmotic_ledger_leaveoneout_deferred_arithmetic")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Independent re-derivation of the urine figure (the only component with a
# traceable mechanistic source) from renal_filtration.py's output,
# rather than trusting the node's "1.415L" number by name alone.
# ---------------------------------------------------------------------------
gfr = json.load(open(os.path.join(OUT_ROOT, "renal_filtration", "renal_filtration_results.json")))
gfr_ml_min = gfr["step1_human_primary_anchors"]["davies_shock_1950_20_29yo"]["gfr"]  # 122.8
REABSORPTION_FRAC = 0.992  # StatPearls-cited, hybrid (not itself re-derived here, disclosed)
urine_ml_day = gfr_ml_min * (1.0 - REABSORPTION_FRAC) * 1440.0
urine_L_day = urine_ml_day / 1000.0
urine_selfcheck_relerr_pct = abs(urine_L_day - 1.415) / 1.415 * 100.0

COMPONENTS = {
    "urine": round(urine_L_day, 4),  # independently re-derived above, not copied
    "skin": 0.353,
    "resp": 0.325,
    "sweat": 0.10,
    "stool": 0.15,
}
ANCHOR_L_DAY = 2.5
TOLERANCE_FRAC = 0.15

total = sum(COMPONENTS.values())
resid_pct = (total - ANCHOR_L_DAY) / ANCHOR_L_DAY * 100.0
baseline_pass = abs(resid_pct) <= TOLERANCE_FRAC * 100.0

# ---------------------------------------------------------------------------
# THE DEFERRED EXTENSION: full leave-one-out void floor (node only ran
# drop-urine; this runs all 5 -- does the closure test have teeth for EVERY
# component, or only the one the node happened to test?)
# ---------------------------------------------------------------------------
loo_results = {}
for name in COMPONENTS:
    remaining = {k: v for k, v in COMPONENTS.items() if k != name}
    remaining_total = sum(remaining.values())
    remaining_resid_pct = (remaining_total - ANCHOR_L_DAY) / ANCHOR_L_DAY * 100.0
    remaining_pass = abs(remaining_resid_pct) <= TOLERANCE_FRAC * 100.0
    loo_results[name] = {
        "dropped_component_L_day": round(COMPONENTS[name], 4),
        "remaining_total_L_day": round(remaining_total, 4),
        "remaining_resid_pct": round(remaining_resid_pct, 2),
        "still_passes_15pct_band": remaining_pass,
    }

n_loo_pass = sum(1 for r in loo_results.values() if r["still_passes_15pct_band"])
n_loo_total = len(loo_results)
closure_has_teeth_for_all_components = (n_loo_pass == 0)  # node's criterion: dropping SHOULD fail

# ---------------------------------------------------------------------------
# reproduce the node's single tested case exactly (drop urine -> 0.928,
# -62.9%, FAILS) as a machine cross-check that this rebuild matches its prose
# ---------------------------------------------------------------------------
drop_urine = loo_results["urine"]
drop_urine_matches_cited = (
    abs(drop_urine["remaining_total_L_day"] - 0.928) < 0.005
    and abs(drop_urine["remaining_resid_pct"] - (-62.9)) < 0.5
    and drop_urine["still_passes_15pct_band"] is False
)

OUT = {
    "cell": "HOLE-WATER-OSMOTIC-LEDGER-CROSSCELL-CONTRADICTION",
    "urine_independent_rederivation": {
        "gfr_ml_min_source": "renal_filtration.py:davies_shock_1950_20_29yo",
        "gfr_ml_min": gfr_ml_min, "reabsorption_frac_cited_StatPearls": REABSORPTION_FRAC,
        "urine_L_day_rederived": round(urine_L_day, 4),
        "urine_selfcheck_relerr_pct_vs_cited_1.415": round(urine_selfcheck_relerr_pct, 3),
    },
    "baseline_closure": {
        "components_L_day": {k: round(v, 4) for k, v in COMPONENTS.items()},
        "total_L_day": round(total, 4), "anchor_L_day": ANCHOR_L_DAY,
        "resid_pct": round(resid_pct, 2), "pass_15pct_band": baseline_pass,
    },
    "full_leave_one_out_voidfloor": loo_results,
    "reproduces_nodes_own_drop_urine_case": drop_urine_matches_cited,
    "n_of_5_dropouts_still_pass_15pct": n_loo_pass,
    "closure_has_teeth_for_ALL_5_components": closure_has_teeth_for_all_components,
    "honest_note": (
        "The node ran only ONE leave-one-out case (drop urine). This extension "
        "runs all 5 and finds " + (
            "the closure test has real teeth regardless of WHICH single "
            "component is dropped (0/5 pass), i.e. the -6.3% baseline match is "
            "not carried disproportionately by one term."
            if closure_has_teeth_for_all_components else
            f"{n_loo_pass}/5 single-component dropouts STILL pass the 15% band "
            "-- meaning the closure test's 'teeth' claim from the node's single "
            "urine-drop case does NOT generalize to every component; some "
            "components are individually small enough that the ledger is "
            "insensitive to them, a genuine weakening of the node's "
            "generality claim that it did not test."
        )
    ),
}

print(json.dumps(OUT, indent=2))
with open(os.path.join(OUT_DIR, "results.json"), "w") as fh:
    json.dump(OUT, fh, indent=2)
