"""Postprandial bile-acid response: sigma-normalized residual / z-score cross-check.

Checks whether the sigma-normalized residual arithmetic of a fitted enterohepatic bile-acid ODE --
comparing its postprandial fold-rise and time-to-peak against two held-out human serum studies --
holds on raw numbers, and whether the pre-registered PASS/FAIL split (fold-rise FAILS, time-to-peak
PASSES) follows correctly from that arithmetic.

SCOPE DECISION (avoiding a too-minimal-rebuild false refutation): the full hepatobiliary/upper-SI/
ileum/portal/serum ODE (calibrated to Voronova/Woodhead/Guiastrennec model papers, with a fitted
meal-pulse, saturable hepatic extraction, and a gallbladder-concentration step) is NOT
reimplemented here -- the model OUTPUT numbers (fold-rise=1.06x, time-to-peak=72min) are given but
not the rate constants needed to regenerate them from scratch; inventing those constants would
manufacture a disagreement that is not independently testable at that resolution. This cell takes
the model's reported outputs as given and rebuilds ONLY the fully-specified statistical comparison:
combining two independent anchor studies (Linnet1983 N=9, LaRusso1974 N=8) into a band, and
z-scoring the model output against it.

STATED NUMBERS (verbatim):
  Model fold-rise = 1.06x;               model time-to-peak = 72 min
  Linnet1983 (PMID6826110, N=9): fold-rise 2.4-4.7x (mean 3.55, sigma 1.15); peak 60-90min
  LaRusso1974 (PMID4851463, N=8): peak 90-120min
  Pre-registered bands: fold-rise [1.5,6.0]x; time-to-peak [45,150]min
  Combined time-to-peak anchor: mean=90min, sigma=30min (spanning both studies' peak windows,
  60-90 and 90-120 -> center ~90, half-range ~30 covering both)

PRE-REGISTERED GATES:
  G1 z_fold = (model_fold - mean_fold)/sigma_fold must reproduce the claimed -2.17 within 5%.
  G2 fold-rise PASS/FAIL: model_fold=1.06x must fall OUTSIDE [1.5,6.0]x (reproducing the claimed
     FAIL verdict).
  G3 z_peak = (model_peak - mean_peak)/sigma_peak must reproduce the claimed -0.60 within 5%.
  G4 time-to-peak PASS/FAIL: model_peak=72min must fall INSIDE [45,150]min (reproducing the
     claimed PASS verdict).
  G5 VOID FLOOR: swap which anchor-band feeds which z-score (score time-to-peak against the
     fold-rise band's mean/sigma, and vice versa, both unit-mismatched) -> must produce z-scores
     that do NOT reproduce either claimed value within 5%, proving G1/G3 are not an artifact of
     "any two mean/sigma numbers give some z-score."

Reads: nothing.
Writes: bile_acid_postprandial_zscore_crosscheck.json
Gate: verdict CONFIRMED iff G1-G5 all pass.
"""
import json
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "bile_acid_postprandial_zscore_crosscheck",
                    "bile_acid_postprandial_zscore_crosscheck.json")

MODEL_FOLD = 1.06
MODEL_PEAK_MIN = 72.0

FOLD_MEAN, FOLD_SIGMA = 3.55, 1.15
FOLD_BAND = (1.5, 6.0)
CLAIMED_Z_FOLD = -2.17

PEAK_MEAN, PEAK_SIGMA = 90.0, 30.0
PEAK_BAND = (45.0, 150.0)
CLAIMED_Z_PEAK = -0.60


def zscore(x, mean, sigma):
    return (x - mean) / sigma


def main():
    gates = {}

    z_fold = zscore(MODEL_FOLD, FOLD_MEAN, FOLD_SIGMA)
    g1 = abs(z_fold - CLAIMED_Z_FOLD) / abs(CLAIMED_Z_FOLD) < 0.05
    gates["G1_z_fold_reproduces_claimed_-2.17"] = bool(g1)

    fold_in_band = FOLD_BAND[0] <= MODEL_FOLD <= FOLD_BAND[1]
    g2 = not fold_in_band  # reproduces the claimed FAIL
    gates["G2_fold_rise_correctly_fails_band"] = bool(g2)

    z_peak = zscore(MODEL_PEAK_MIN, PEAK_MEAN, PEAK_SIGMA)
    g3 = abs(z_peak - CLAIMED_Z_PEAK) / abs(CLAIMED_Z_PEAK) < 0.05
    gates["G3_z_peak_reproduces_claimed_-0.60"] = bool(g3)

    peak_in_band = PEAK_BAND[0] <= MODEL_PEAK_MIN <= PEAK_BAND[1]
    g4 = peak_in_band  # reproduces the claimed PASS
    gates["G4_time_to_peak_correctly_passes_band"] = bool(g4)

    # G5 void floor: swap the anchor stats between the two quantities (unit-mismatched sigma)
    z_fold_swapped = zscore(MODEL_FOLD, PEAK_MEAN, PEAK_SIGMA)
    z_peak_swapped = zscore(MODEL_PEAK_MIN, FOLD_MEAN, FOLD_SIGMA)
    swap_matches_fold = abs(z_fold_swapped - CLAIMED_Z_FOLD) / abs(CLAIMED_Z_FOLD) < 0.05
    swap_matches_peak = abs(z_peak_swapped - CLAIMED_Z_PEAK) / abs(CLAIMED_Z_PEAK) < 0.05
    g5 = (not swap_matches_fold) and (not swap_matches_peak)
    gates["G5_void_floor_swapped_anchors_fail_to_reproduce"] = bool(g5)

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "CONFIRMED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "MODEL-BILE-ACID-ENTEROHEPATIC",
        "z_fold_computed": z_fold,
        "z_fold_claimed": CLAIMED_Z_FOLD,
        "z_peak_computed": z_peak,
        "z_peak_claimed": CLAIMED_Z_PEAK,
        "fold_in_band": fold_in_band,
        "peak_in_band": peak_in_band,
        "void_floor_z_fold_swapped": z_fold_swapped,
        "void_floor_z_peak_swapped": z_peak_swapped,
        "gates": gates,
        "verdict": verdict,
        "scope_note": ("The full hepatobiliary/ileum/portal ODE that PRODUCES the 1.06x/72min "
                        "model outputs is taken as given (not re-derived from rate constants, "
                        "which the claim text does not fully specify) -- this cell certifies "
                        "only that the sigma-normalized residual arithmetic and its "
                        "resulting FAIL/PASS split are computed correctly from raw anchor numbers."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict, "z_fold": z_fold, "z_peak": z_peak}, indent=2))
    return result


if __name__ == "__main__":
    main()
