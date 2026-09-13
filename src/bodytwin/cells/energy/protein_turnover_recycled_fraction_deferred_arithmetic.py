"""
Deferred-arithmetic execution of the protein-turnover RECYCLED-FRACTION estimate.

A prior estimate gave only the 4 corner combinations of "1 - intake/turnover" (3-4 g/kg/d turnover
x 0.8-1.0 g/kg/d intake = 66.7-80.0 pct, mid 74.3 pct) against the cited 75-90 pct anchor. This cell
rebuilds that 4-point grid as real code, then extends it to a FULL fine grid (the deferred part: the
interior was never checked, and the FRACTION of the plausible parameter rectangle that actually
lands inside the cited 75-90% Waterlow/Young anchor band was never reported).

PRE-REGISTERED GATE:
  C = the recycled-fraction estimate (1 - intake/turnover) lands inside the
      cited 75-90% anchor band for a MAJORITY (>50%) of the plausible
      (turnover in [3,4] g/kg/d, intake in [0.8,1.0] g/kg/d) parameter
      rectangle, i.e. the "mid 74.3 pct vs 75-90 pct" near-miss is not a
      fragile artifact of averaging over a mostly-failing grid.
  not-C = <50% of the grid lands inside the band -- the "near-miss" framing
      obscures that MOST of the physiologically plausible parameter space
      actually falls short of the cited anchor.

Reads: nothing.
Writes: <OUT_ROOT>/protein_turnover_recycled_fraction_deferred_arithmetic/results.json
Gate: pre_registered_gate.C_holds (frac_in_band > 0.50).
"""
import subprocess
import json
import os
import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "protein_turnover_recycled_fraction_deferred_arithmetic")
os.makedirs(OUT_DIR, exist_ok=True)

CELLS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

grep = subprocess.run(
    ["grep", "-rl", "-i", "recycled.fraction\\|recycled_fraction", "--include=*.py",
     CELLS_DIR],
    capture_output=True, text=True
)
pre_existing_hits = [f for f in grep.stdout.strip().split("\n")
                      if f and "protein_turnover_recycled_fraction" not in f]

ANCHOR_BAND = (75.0, 90.0)  # Waterlow/Young, cited

# ---------------------------------------------------------------------------
# STEP 1: reproduce the cited 4-corner grid exactly
# ---------------------------------------------------------------------------
corners_turnover = [3.0, 4.0]
corners_intake = [0.8, 1.0]
corner_vals = []
for t in corners_turnover:
    for i in corners_intake:
        recycled_pct = (1.0 - i / t) * 100.0
        corner_vals.append({"turnover_g_kg_d": t, "intake_g_kg_d": i, "recycled_pct": round(recycled_pct, 2)})

corner_min = min(c["recycled_pct"] for c in corner_vals)
corner_max = max(c["recycled_pct"] for c in corner_vals)
corner_mean = float(np.mean([c["recycled_pct"] for c in corner_vals]))
cited_range_reproduced = (abs(corner_min - 66.7) < 0.1) and (abs(corner_max - 80.0) < 0.1)
cited_mid_reproduced = abs(corner_mean - 74.3) < 0.2

# ---------------------------------------------------------------------------
# STEP 2: THE DEFERRED EXTENSION -- fine interior grid + fraction-in-band
# ---------------------------------------------------------------------------
turnover_grid = np.round(np.arange(3.0, 4.001, 0.1), 2)
intake_grid = np.round(np.arange(0.8, 1.001, 0.02), 2)

fine_vals = []
for t in turnover_grid:
    for i in intake_grid:
        recycled_pct = (1.0 - i / t) * 100.0
        in_band = ANCHOR_BAND[0] <= recycled_pct <= ANCHOR_BAND[1]
        fine_vals.append({"turnover": float(t), "intake": float(i),
                           "recycled_pct": round(float(recycled_pct), 3), "in_band": bool(in_band)})

n_total = len(fine_vals)
n_in_band = sum(1 for v in fine_vals if v["in_band"])
frac_in_band = n_in_band / n_total

# where exactly does the boundary fall? (turnover value at which the band
# just starts being reachable, at each intake level)
boundary_by_intake = {}
for i in intake_grid:
    row = [v for v in fine_vals if abs(v["intake"] - i) < 1e-9]
    lo_t = min(v["turnover"] for v in row if v["in_band"]) if any(v["in_band"] for v in row) else None
    boundary_by_intake[str(i)] = {
        "min_turnover_entering_band": lo_t,
        "recycled_pct_at_turnover_3.0": next(v["recycled_pct"] for v in row if abs(v["turnover"] - 3.0) < 1e-9),
        "recycled_pct_at_turnover_4.0": next(v["recycled_pct"] for v in row if abs(v["turnover"] - 4.0) < 1e-9),
    }

C_holds = frac_in_band > 0.50

OUT = {
    "cell": "protein_turnover_recycled_fraction_deferred_arithmetic",
    "preexisting_backing_script_check": {
        "grep_hits": pre_existing_hits, "confirmed_no_preexisting_script": len(pre_existing_hits) == 0,
        "manual_disambiguation": (
            "protein_turnover_atp_energy_ledger.py:step2_recycled_fraction_crux computes a "
            "DIFFERENT ratio (urea-N/intake-N vs turnover-N/intake-N, a G3 2-4x scale-crux "
            "gate) -- inspected line-by-line, it does NOT compute the recycled_pct=1-intake/"
            "turnover fine-grid-fraction-in-band check this cell performs. Not a duplicate; "
            "a genuinely distinct sub-question."
        ),
    },
    "step1_reproduce_cited_4corner_grid": {
        "corner_vals": corner_vals,
        "range_pct": (round(corner_min, 2), round(corner_max, 2)), "mean_pct": round(corner_mean, 2),
        "cited_range_66.7_80.0_reproduced": cited_range_reproduced,
        "cited_mid_74.3_reproduced": cited_mid_reproduced,
    },
    "step2_deferred_finegrid_fraction_in_band": {
        "anchor_band_pct": ANCHOR_BAND,
        "grid_n": n_total, "n_in_band": n_in_band, "frac_in_band": round(frac_in_band, 4),
        "boundary_by_intake_level": boundary_by_intake,
    },
    "pre_registered_gate": {"threshold": "frac_in_band > 0.50", "C_holds": C_holds},
    "verdict": (
        f"{round(frac_in_band*100,1)}% of the fine (turnover x intake) grid lands inside "
        f"the cited 75-90% Waterlow/Young anchor band -- "
        + ("CONFIRMS the 'near-miss, mid 74.3 vs 75-90' framing is representative "
           "of most of the plausible parameter space, not just a lucky corner average."
           if C_holds else
           "REFUTES the implicit 'near-miss' framing: MOST of the physiologically plausible "
           "parameter rectangle actually falls SHORT of the cited anchor band -- recycled_pct "
           "= 1-intake/turnover is monotonically INCREASING in turnover and DECREASING in "
           "intake, so only the HIGH-turnover (4g/kg/d, the upper/less-commonly-cited edge of "
           "the band) / LOW-intake corner clears 75%; the more commonly-cited lower turnover "
           "estimates (~3g/kg/d) fall short across the ENTIRE intake range tested (66.7-73.3%, "
           "machine-verified above), meaning the 'mid 74.3 pct vs 75-90 pct near-miss' "
           "framing is driven by averaging in the less-representative high-turnover corner.")
    ),
}

print(json.dumps(OUT, indent=2))
with open(os.path.join(OUT_DIR, "results.json"), "w") as fh:
    json.dump(OUT, fh, indent=2)
