#!/usr/bin/env python3
"""Gut crypt renewal time by Little's Law: tau = Ts / LI (steady-state population balance
applied to a labeling-index / S-phase-duration pair) from Potten et al 1992's two human colonic
BrdU-labeling cohorts, compared against a fitted compartmental ODE model's reported renewal times.

The compartmental ODE model (PMID:17360468, cross-checked against BioModels BIOMD0000000520) is
NOT reimplemented here: its published summary gives the model OUTPUTS (tau_all=1.179d,
tau_shed=1.453d) but not its transition-rate parameters, so those are used as given.

Reads: nothing (all numbers are published cohort values embedded below).
Writes: gut_crypt_littles_law_renewal.json under the cell output directory.
Gate: G1-G5 below.

STATED NUMBERS (verbatim, Potten CS et al. 1992 Gut):
  Cohort 1 (PMID 1740282): LI=10.0%, Ts=8.6h (same-cohort flow cytometry)
  Cohort 2 (PMID 1316306, companion cohort): LI=10.3% (colon)
  Model (reported, from the fitted compartmental ODE): tau_all=1.179d, tau_shed=1.453d

PRE-REGISTERED GATES:
  G1 tau1 = Ts / LI1 (in days) must reproduce the claimed 3.58d within 1%.
  G2 tau2 = Ts / LI2 must reproduce the claimed 3.48d within 1%.
  G3 mean(tau1, tau2) must reproduce the claimed 3.53d within 1%.
  G4 ratio_all = tau_all / mean(tau1,tau2) must reproduce the claimed 0.334x within 5%; ratio_shed
     must reproduce 0.412x within 5% (verifies the DISAGREEMENT arithmetic itself -- symmetric QC
     requires the same evidence burden for a confirmed disagreement as for a match).
  G5 VOID FLOOR: compute the inverted ratio LI/Ts (wrong Little's-Law direction: labeling index per
     hour instead of hours per labeling-fraction) -> must NOT land anywhere near a 3-4 day renewal
     time (off by >100x), proving G1-G3's near-day-scale answer is not a coincidence of "any
     ratio of these two small numbers looks like a few days."

"""
import json
import os

import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "gut_crypt_littles_law_renewal",
                    "gut_crypt_littles_law_renewal.json")

TS_HOURS = 8.6
LI1 = 0.100   # cohort 1, PMID 1740282
LI2 = 0.103   # cohort 2, PMID 1316306

CLAIMED_TAU1_D = 3.58
CLAIMED_TAU2_D = 3.48
CLAIMED_MEAN_D = 3.53
CLAIMED_TAU_ALL_D = 1.179
CLAIMED_TAU_SHED_D = 1.453
CLAIMED_RATIO_ALL = 0.334
CLAIMED_RATIO_SHED = 0.412


def tau_days(Ts_hours, LI):
    """Little's Law population-balance renewal time: mean sojourn time = (measured transit-phase
    duration) / (fraction of the population observed in that phase at steady state)."""
    return (Ts_hours / LI) / 24.0


def main():
    gates = {}

    tau1 = tau_days(TS_HOURS, LI1)
    tau2 = tau_days(TS_HOURS, LI2)
    mean_tau = (tau1 + tau2) / 2.0

    g1 = abs(tau1 - CLAIMED_TAU1_D) / CLAIMED_TAU1_D < 0.01
    g2 = abs(tau2 - CLAIMED_TAU2_D) / CLAIMED_TAU2_D < 0.01
    g3 = abs(mean_tau - CLAIMED_MEAN_D) / CLAIMED_MEAN_D < 0.01
    gates["G1_tau1_reproduces_3.58d"] = bool(g1)
    gates["G2_tau2_reproduces_3.48d"] = bool(g2)
    gates["G3_mean_reproduces_3.53d"] = bool(g3)

    ratio_all = CLAIMED_TAU_ALL_D / mean_tau
    ratio_shed = CLAIMED_TAU_SHED_D / mean_tau
    g4a = abs(ratio_all - CLAIMED_RATIO_ALL) / CLAIMED_RATIO_ALL < 0.05
    g4b = abs(ratio_shed - CLAIMED_RATIO_SHED) / CLAIMED_RATIO_SHED < 0.05
    gates["G4_disagreement_ratio_all_reproduces_0.334x"] = bool(g4a)
    gates["G4_disagreement_ratio_shed_reproduces_0.412x"] = bool(g4b)

    # G5 void floor: inverted (wrong-direction) Little's Law
    tau1_void = (LI1 / TS_HOURS) / 24.0
    void_far_off = (mean_tau / tau1_void) > 100.0
    gates["G5_void_floor_inverted_law_off_by_100x"] = bool(void_far_off)

    gates = {k: bool(v) for k, v in gates.items()}
    # The published claim is an ALREADY-DISCLOSED disagreement (model under-predicts renewal
    # time vs the human-anchor Little's-Law estimate) -- so "matching" means reproducing that
    # disagreement's arithmetic exactly, not finding agreement. Verdict reflects arithmetic fidelity.
    verdict = "CONFIRMED_DISAGREEMENT_REPRODUCED" if all(gates.values()) else "DISAGREE_WITH_NODES_OWN_ARITHMETIC"

    result = {
        "node_id": "MODEL-GUT-CRYPT-DYNAMICS",
        "tau1_days_computed": tau1,
        "tau2_days_computed": tau2,
        "mean_tau_days_computed": mean_tau,
        "ratio_all_computed": ratio_all,
        "ratio_shed_computed": ratio_shed,
        "tau1_void_floor_days": tau1_void,
        "gates": gates,
        "verdict": verdict,
        "scope_note": ("The compartmental ODE model's tau_all/tau_shed (1.179d/1.453d) are "
                        "taken from the published summary as given, NOT independently re-derived from "
                        "rate constants -- rebuilding the full ODE from a text summary would risk "
                        "a too-minimal reimplementation. This script certifies only that the "
                        "Little's-Law anchor arithmetic and the resulting disagreement ratios are "
                        "correctly computed from Potten 1992's directly-reported (LI,Ts)."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict, "mean_tau_days": mean_tau,
                       "ratio_all": ratio_all, "ratio_shed": ratio_shed}, indent=2))
    return result


if __name__ == "__main__":
    main()
