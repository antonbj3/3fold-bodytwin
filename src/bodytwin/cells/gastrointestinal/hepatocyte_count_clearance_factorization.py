"""Hepatocyte-count / clearance factorization check.

Tests whether the whole-liver factorization CL ~ hepatocellularity (HPGL, count/g) x liver mass (g)
x per-cell intrinsic clearance (CLint) arithmetically reproduces the stated total-hepatocyte-count
range and the Wilson 2003 individual-variation fold, and whether the two independent measurement
METHODS (direct counting-chamber count vs meta-analytic pooled estimate) converge within a modest
tolerance.

STATED NUMBERS (verbatim):
  Wilson2003 (PMID12968989, n=7, direct collagenase-perfusion + counting-chamber count):
    geometric mean 107e6 cells/g, individual range 65e6-185e6 cells/g (claimed 2.85x)
  Barter2007 (meta-analysis): weighted mean 99e6 cells/g
  Adult liver mass: 1400-1800g
  Claimed total count: 1.39e11-1.78e11 (~140-180 billion) hepatocytes
  Retrospective IVIVE clearance underprediction (PMID21058916 + corroborating 52/65-drug
  compilations): ~4-fold average, high-clearance CYP probes

PRE-REGISTERED GATES:
  G1 total_count = HPGL_Barter_mean(99e6) * liver_mass_range[1400,1800]g must reproduce the claimed
     [1.39e11, 1.78e11] range within 1%.
  G2 Wilson2003 individual fold = 185e6/65e6 must reproduce the claimed 2.85x within 1%.
  G3 CONVERGENCE (not a tautology gate): the two independent methods (Wilson2003 direct count
     107e6/g vs Barter2007 meta-analytic 99e6/g) must agree within a pre-registered 15% band --
     tighter than the within-Wilson2003 individual spread (2.85x=185%) would allow by chance,
     demonstrating the two decorrelated methods are not merely both landing in a wide plausible
     range.
  G4 VOID FLOOR: recompute total_count using a liver-mass typo one decimal order too small
     (140-180g instead of 1400-1800g, a real unit-confusion error class: g vs hg) -> must NOT
     overlap the claimed [1.39e11,1.78e11] range (proves G1's pass isn't "any product of two
     O(1e8)/O(1e3)-scale numbers lands near 1e11").

Reads: nothing.
Writes: hepatocyte_count_clearance_factorization.json
Gate: verdict CONFIRMED iff G1-G4 all pass.
"""
import json
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "hepatocyte_count_clearance_factorization")
OUT = _os.path.join(OUT_DIR, "hepatocyte_count_clearance_factorization.json")

HPGL_WILSON_GEOMEAN = 107e6
HPGL_WILSON_LO, HPGL_WILSON_HI = 65e6, 185e6
HPGL_BARTER_MEAN = 99e6

LIVER_MASS_LO_G, LIVER_MASS_HI_G = 1400.0, 1800.0

CLAIMED_TOTAL_LO, CLAIMED_TOTAL_HI = 1.39e11, 1.78e11
CLAIMED_WILSON_FOLD = 2.85
CONVERGENCE_BAND = 0.15  # 15%


def main():
    gates = {}

    total_lo = HPGL_BARTER_MEAN * LIVER_MASS_LO_G
    total_hi = HPGL_BARTER_MEAN * LIVER_MASS_HI_G
    err_lo = abs(total_lo - CLAIMED_TOTAL_LO) / CLAIMED_TOTAL_LO
    err_hi = abs(total_hi - CLAIMED_TOTAL_HI) / CLAIMED_TOTAL_HI
    g1 = (err_lo < 0.01) and (err_hi < 0.01)
    gates["G1_total_count_reproduces_claimed_range_1pct"] = bool(g1)

    wilson_fold = HPGL_WILSON_HI / HPGL_WILSON_LO
    g2 = abs(wilson_fold - CLAIMED_WILSON_FOLD) / CLAIMED_WILSON_FOLD < 0.01
    gates["G2_wilson_individual_fold_reproduces_2.85x"] = bool(g2)

    rel_diff_methods = abs(HPGL_WILSON_GEOMEAN - HPGL_BARTER_MEAN) / HPGL_BARTER_MEAN
    g3 = rel_diff_methods < CONVERGENCE_BAND
    gates["G3_two_independent_methods_converge_within_15pct"] = bool(g3)

    # G4 void floor: unit-confusion decade error on liver mass
    total_void_lo = HPGL_BARTER_MEAN * (LIVER_MASS_LO_G / 10.0)
    total_void_hi = HPGL_BARTER_MEAN * (LIVER_MASS_HI_G / 10.0)
    overlaps = not (total_void_hi < CLAIMED_TOTAL_LO or total_void_lo > CLAIMED_TOTAL_HI)
    g4 = not overlaps
    gates["G4_void_floor_decade_error_does_not_overlap"] = bool(g4)

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "CONFIRMED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "COUNT-HEPATOCYTE-NUMBER-CLEARANCE-FACTORIZATION",
        "total_count_computed_range": [total_lo, total_hi],
        "wilson_fold_computed": wilson_fold,
        "rel_diff_between_methods": rel_diff_methods,
        "void_floor_total_count_range": [total_void_lo, total_void_hi],
        "gates": gates,
        "verdict": verdict,
        "scope_note": ("The ~4-fold retrospective IVIVE clearance-underprediction figure "
                        "(PMID21058916) is taken as a cited external validation-literature number, "
                        "not independently re-derived here (no per-drug microsome/hepatocyte "
                        "concentration-time data is used here) -- this cell certifies "
                        "the count x mass factorization arithmetic and the cross-method count "
                        "convergence only."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict,
                       "total_count_computed_range": [total_lo, total_hi],
                       "wilson_fold": wilson_fold}, indent=2))
    return result


if __name__ == "__main__":
    main()
