#!/usr/bin/env python3
"""RESTING CARDIAC OUTPUT SPREAD GATE: a pre-registered bound on a deliberate three-anchor fork.

Three resting cardiac-output values circulate in this cell graph:
  (a) the population-level GEOMETRIC route (CMR-measured sex-averaged stroke volume x HR=65 bpm),
      read by fluid_compartments, arterial_pressure, renal_filtration, venous_return and
      exercise_bloodflow_redistribution;
  (b) Higginbotham et al. 1986 (PMID 3948345), a real direct right-heart catheterization
      measurement, n=24, used as a literal by coronary_absolute_flow;
  (c) the reference-trial FICK-CHAIN estimate (VO2 / a-vO2 difference = 5.0) of the cardiac_output
      cell, read by pulmonary_gas_exchange and blood_oxygen_transport.
The fork is deliberate: (b) is a decorrelated real-measurement cross-check against the two modelled
routes, and (a) vs (c) is a scope split (population-generic vs trial-specific), each consumer
reading the value that matches its own scope. This cell does not collapse the anchors; it adds the
machine-checked bound on their drift that the design previously lacked.

PRE-REGISTERED THRESHOLD (fixed before the first run): MAX_PAIRWISE_SPREAD_PCT = 5.0. Independent
cardiac-output techniques (thermodilution vs direct Fick) are clinically considered to agree within
10-20% (Narang et al. 2022, PMID 35613956); requiring these three decorrelated anchors to sit
inside a tighter 5% band is a stronger, failable test.

Reads: <BODYTWIN_OUT>/cardiac_output_geometric/cardiac_output_geometric_results.json and
<BODYTWIN_OUT>/cardiac_output/cardiac_output_results.json (both required, never written to).
Writes: <BODYTWIN_OUT>/cardiac_output_resting_spread_gate/cardiac_output_resting_spread_gate_results.json
Gate: resting_co_spread_within_preregistered_bound; exit 0 on pass, 2 on fail.
"""
import json
import os
import sys

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
CARDIAC_GEO_JSON = os.path.join(OUT_ROOT, "cardiac_output_geometric", "cardiac_output_geometric_results.json")
CARDIAC_SUBJ_JSON = os.path.join(OUT_ROOT, "cardiac_output", "cardiac_output_results.json")
OUT_DIR = os.path.join(OUT_ROOT, "cardiac_output_resting_spread_gate")

# ---- pre-registered threshold (declared before measuring) ---------------------------------------
MAX_PAIRWISE_SPREAD_PCT = 5.0


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def pct_spread(a, b):
    """Symmetric percent spread, denominator = pairwise mean (same convention as
    cardiac_output.py's route_diff_pct)."""
    return abs(a - b) / ((a + b) / 2.0) * 100.0


def main():
    for p in (CARDIAC_GEO_JSON, CARDIAC_SUBJ_JSON):
        if not os.path.exists(p):
            print(f"FAIL: required input missing: {p} -- run its producing script first.")
            return 1
    os.makedirs(OUT_DIR, exist_ok=True)

    geo = load_json(CARDIAC_GEO_JSON)
    subj = load_json(CARDIAC_SUBJ_JSON)

    co_geometric = geo["rest_co_flow_leg"]["central_estimate_sex_avg_hr65_l_min"]           # (a) 5.5575
    co_higginbotham_real = geo["rest_co_flow_leg"]["higginbotham_real_rest_co_l_min"]        # (b) 5.70
    co_subject_fick = subj["fick_cardiac_output"]["q_rest_l_min"]                            # (c) 5.653

    anchors = {
        "geometric_population_central_estimate": co_geometric,
        "higginbotham_1986_real_catheterization": co_higginbotham_real,
        "reference_walking1_fick_chain": co_subject_fick,
    }
    print("=" * 78)
    print("GATE 2 -- three-way resting cardiac output spread")
    print("=" * 78)
    for name, v in anchors.items():
        print(f"  {name:40s} = {v:.4f} L/min")

    names = list(anchors.keys())
    pairwise = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            key = f"{names[i]}__vs__{names[j]}"
            spread = pct_spread(anchors[names[i]], anchors[names[j]])
            pairwise[key] = spread
            print(f"  spread[{names[i]} vs {names[j]}] = {spread:.3f}%")

    max_spread_pct = max(pairwise.values())
    max_pair = max(pairwise, key=pairwise.get)
    gate_pass = max_spread_pct < MAX_PAIRWISE_SPREAD_PCT
    print(f"\nMax pairwise spread = {max_spread_pct:.3f}% ({max_pair}) vs pre-registered threshold "
          f"<{MAX_PAIRWISE_SPREAD_PCT}%: {'PASS' if gate_pass else 'FAIL'}")
    print("(This gate PASSING confirms the three-anchor DESIGN is healthy -- decorrelated cross-check, "
          "not an unlabelled drift. It does NOT collapse the anchors; it bounds their disagreement.)")

    report = {
        "anchors_l_min": anchors,
        "pairwise_spread_pct": pairwise,
        "max_pairwise_spread_pct": max_spread_pct,
        "max_pairwise_pair": max_pair,
        "max_pairwise_spread_threshold_pct_PREREGISTERED": MAX_PAIRWISE_SPREAD_PCT,
        "gate_spread_within_threshold": bool(gate_pass),
        "adjudication": "deliberate 3-anchor decorrelated cross-check (1 real measurement + 2 models "
                         "at 2 different scopes: population-generic vs trial-specific); NOT "
                         "collapsed, NOT staleness -- see docstring.",
        "consumers": {
            "geometric_population_central_estimate": [
                "fluid_compartments.py", "arterial_pressure.py", "renal_filtration.py",
                "venous_return.py", "exercise_bloodflow_redistribution.py",
            ],
            "higginbotham_1986_real_catheterization": ["coronary_absolute_flow.py"],
            "reference_walking1_fick_chain": ["pulmonary_gas_exchange.py", "blood_oxygen_transport.py"],
        },
    }
    gates = {"resting_co_spread_within_preregistered_bound": bool(gate_pass)}
    report["gates"] = gates
    report["overall_pass"] = gate_pass

    out_path = f"{OUT_DIR}/cardiac_output_resting_spread_gate_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    return 0 if gate_pass else 2


if __name__ == "__main__":
    sys.exit(main())
