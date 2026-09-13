"""
Node: RENAL-GLOMERULAR-FILTRATION

Claim under test: "Gansevoort 2011 (PMID21289597, n=1,019,017 pooled general+high-risk cohorts):
eGFR and ACR independently/multiplicatively predict incident ESRD up to HR454/HR28; replicated
in CKD-diagnosed cohorts (Astor 2011, PMID21289598, n=21,688, HR6.24)". No script anywhere
This cell computes whether Astor's CONTINUOUS per-unit hazard (HR=6.24 per 15 mL/min/1.73m2 lower
eGFR, below 45) actually predicts Gansevoort's directly-observed CATEGORICAL step (HR29 at eGFR45 ->
HR454 at eGFR15, both vs an eGFR~95 reference).

Reads: nothing. Writes: OUT_ROOT/renal_egfr_hr_compounding_crosscheck/renal_egfr_hr_compounding_crosscheck.json
plus its void-floor draws in the same directory. The [1/3,3] convergence gate decides.

PRE-REGISTERED (before computing): eGFR45 -> eGFR15 is a drop of 30 mL/min/1.73m2 = exactly 2
units of Astor's 15-unit step. If the per-unit hazard is roughly constant across this range
(the node's claim that eGFR and ACR act multiplicatively/independently, no eGFR x ACR
interaction), Astor's continuous rate COMPOUNDED over 2 units should predict a HR ratio in the
same order of magnitude as Gansevoort's directly-reported categorical ratio
(454/29 = 15.66x), even though Astor (n=21,688, CKD-diagnosed only) and Gansevoort (n=1,019,017,
pooled general+high-risk) are disjoint cohorts/populations/study designs.

GATE (fixed before computing): predicted/observed ratio (Astor-compounded HR-ratio divided by
Gansevoort's HR29->HR454 ratio) must fall in [1/3, 3] -- the same tolerance band used
for the tendon CAGR magnitude-convergence check, chosen for the same reason: these are
different cohorts/eras/populations, so exact identity is not expected, but order-of-magnitude
agreement is the real test of "not an artifact of one equation/cohort."

FORCED ADVERSARY: a per-unit hazard estimated in a narrow, sicker (CKD-diagnosed) cohort could
extrapolate very differently in a broader general-risk population purely from nonlinearity /
selection, with no real cross-cohort agreement. The void-floor below asks how often two
INDEPENDENT log-uniform-drawn per-unit/categorical hazard pairs would land in the same [1/3,3]
band purely by chance, given the wide range of hazard magnitudes plausible in nephrology
epidemiology.

SECONDARY leg (equation-independence, reported as-is, not gated pass/fail since the node itself
flags it as borderline): Shlipak 2013 (PMID24004120, n=93,710) NRI swapping creatinine-eGFR for
cystatin-C-eGFR = 0.23 for death (95%CI 0.18-0.28, robust) but only 0.10 for ESRD (95%CI
0.00-0.21, touches null) -- reported honestly, not folded into the primary gate.
"""
import json
import os
import random
import math

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "renal_egfr_hr_compounding_crosscheck")
SCRATCH_DIR = OUT_DIR
OUT_PATH = os.path.join(OUT_DIR, "renal_egfr_hr_compounding_crosscheck.json")

ASTOR_HR_PER_15UNIT = 6.24  # PMID21289598, n=21688, CKD-diagnosed
ASTOR_UNITS_45_TO_15 = 2  # (45-15)/15 = 2

GANSEVOORT_HR_EGFR45 = 29.0  # PMID21289597, n=1,019,017
GANSEVOORT_HR_EGFR15 = 454.0

SHLIPAK_NRI_DEATH = 0.23
SHLIPAK_NRI_DEATH_CI = (0.18, 0.28)
SHLIPAK_NRI_ESRD = 0.10
SHLIPAK_NRI_ESRD_CI = (0.00, 0.21)

GATE_LO, GATE_HI = 1.0 / 3.0, 3.0


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    astor_predicted_ratio = ASTOR_HR_PER_15UNIT ** ASTOR_UNITS_45_TO_15
    gansevoort_observed_ratio = GANSEVOORT_HR_EGFR15 / GANSEVOORT_HR_EGFR45
    convergence_ratio = astor_predicted_ratio / gansevoort_observed_ratio
    gate = GATE_LO <= convergence_ratio <= GATE_HI

    result = {
        "node": "RENAL-GLOMERULAR-FILTRATION",
        "astor_2011_pmid21289598": {
            "hr_per_15unit_egfr_decline": ASTOR_HR_PER_15UNIT,
            "units_compounded_egfr45_to_egfr15": ASTOR_UNITS_45_TO_15,
            "predicted_ratio": round(astor_predicted_ratio, 3),
        },
        "gansevoort_2011_pmid21289597": {
            "hr_egfr45": GANSEVOORT_HR_EGFR45,
            "hr_egfr15": GANSEVOORT_HR_EGFR15,
            "observed_ratio": round(gansevoort_observed_ratio, 3),
        },
        "convergence_ratio_predicted_over_observed": round(convergence_ratio, 4),
        "GATE_convergence_in_[1/3,3]": gate,
        "verdict": "MAGNITUDE-CONVERGENT (cross-cohort compounding survives)" if gate else "DIVERGED",
        "secondary_leg_shlipak_2013_pmid24004120_NRI_equation_independence": {
            "death_NRI": SHLIPAK_NRI_DEATH, "death_95CI": SHLIPAK_NRI_DEATH_CI,
            "esrd_NRI": SHLIPAK_NRI_ESRD, "esrd_95CI": SHLIPAK_NRI_ESRD_CI,
            "note": "robust for death (CI excludes 0), borderline for ESRD alone (CI touches 0) -- reported as-is, not gated",
        },
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result, convergence_ratio


def void_floor(n_draws=4000, seed=20260728, real_ratio=None):
    """Null: two INDEPENDENT log-uniform-drawn hazard magnitudes for the same nominal
    2-unit-eGFR-decline comparison, spanning the wide plausible range seen in nephrology
    epi (per-unit HR in [1.5,10], categorical ratio in [3,50] i.e. compounded-2-unit vs a
    directly observed step) -- what fraction of unrelated draws land in [1/3,3]?"""
    rng = random.Random(seed)
    passes = 0
    per_unit_draws = []
    for _ in range(n_draws):
        per_unit_hr = math.exp(rng.uniform(math.log(1.5), math.log(10.0)))
        predicted = per_unit_hr ** ASTOR_UNITS_45_TO_15
        observed = math.exp(rng.uniform(math.log(3.0), math.log(50.0)))
        per_unit_draws.append(per_unit_hr)
        r = predicted / observed
        if GATE_LO <= r <= GATE_HI:
            passes += 1
    moved = (max(per_unit_draws) - min(per_unit_draws)) > 1.0
    out = {
        "n_draws": n_draws,
        "per_unit_hr_sweep": [1.5, 10.0], "observed_ratio_sweep": [3.0, 50.0],
        "substitution_landed": moved,
        "void_pass_rate": round(passes / n_draws, 4),
        "real_ratio": real_ratio,
    }
    with open(os.path.join(SCRATCH_DIR, "void_floor.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    res, ratio = main()
    void_floor(real_ratio=res["convergence_ratio_predicted_over_observed"])
