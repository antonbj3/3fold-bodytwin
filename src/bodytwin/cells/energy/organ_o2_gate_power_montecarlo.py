"""Statistical power check for the combined organ-O2 gate (sum-ratio AND heart-specific) against
two null models, mirroring the void-floor sweep methodology (wired flows held fixed, ERO2 drawn
from the null distribution, N=200000 per null).

Reads: <OUT_ROOT>/cerebral_autoregulation/, /renal_filtration/,
/exercise_bloodflow_redistribution/ and /blood_oxygen_transport/ result JSONs.
Writes: nothing (stdout only).
Gate: the reported pass rates -- a useful gate must pass rarely under null 1 (fully uninformed
ERO2) and null 2 (heart-blind resting-tissue ERO2) while passing on the measured bands.
"""
import json, os, random

ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))

def load(rel):
    with open(os.path.join(ROOT, rel)) as f:
        return json.load(f)

cereb = load("cerebral_autoregulation/cerebral_autoregulation_results.json")
renal = load("renal_filtration/renal_filtration_results.json")
exber = load("exercise_bloodflow_redistribution/exercise_bloodflow_redistribution_results.json")
bot = load("blood_oxygen_transport/blood_oxygen_transport_results.json")

cao2_frac = bot["chemistry_rest"]["cao2_ml_dl"] / 100.0
brain_flow = cereb["resting_cbf_baseline"]["task_targets"]["whole"] / 100.0 * 1350.0
kidney_flow = renal["step5_cardiovascular_coupling"]["rbf_l_min"] * 1000.0
liver_flow = exber["rest_state_distribution"]["splanchnic"]["rest_l_min_central"] * 1000.0
muscle_flow = exber["rest_state_distribution"]["muscle"]["central_l_min"] * 1000.0
heart_flow_lo, heart_flow_hi = 198.6, 331.0
fick_ff = bot["fick_closure_falsifier"]
whole_body_fick = fick_ff["vo2_reconstructed_ml_min"]

# The three literals below carry the provenance of the organ_o2_consumption_partition and coronary
# flow cells, which cite the same numbers to two named primary sources. Same physical quantity,
# same units, same regime (heart-mass-normalized resting MVO2 anchors) at both sites.
HEART_MASS_G = 331.0  # Molina & DiMaio 2012, PMID 22182983, Am J Forensic Med Pathol 33(4):362-7 (n=232, mean 331+-56.7g)
anchor_direct_fick = (6.0 * HEART_MASS_G / 100.0, 10.0 * HEART_MASS_G / 100.0)  # Gobel1978 direct N2O-Fick PMID624164 + wholebodyfrac alloc, mL/100g/min
anchor_mechanical = (10.0 * HEART_MASS_G / 100.0, 20.0 * HEART_MASS_G / 100.0)  # Suga pressure-volume-area (PVA) mechanical MVO2 family, mL/100g/min

def overlap(lo, hi, alo, ahi):
    return max(0.0, min(hi, ahi) - max(lo, alo))

N = 200000
rng = random.Random(12345)

def run_null(brain_rng, kidney_rng, liver_rng, muscle_rng, heart_rng, label):
    pass_sum_only = 0
    pass_combined = 0
    for _ in range(N):
        e_brain = rng.uniform(*brain_rng)
        e_kidney = rng.uniform(*kidney_rng)
        e_liver = rng.uniform(*liver_rng)
        e_muscle = rng.uniform(*muscle_rng)
        e_heart = rng.uniform(*heart_rng)
        heart_flow = rng.uniform(heart_flow_lo, heart_flow_hi)

        vo2_brain = brain_flow * cao2_frac * e_brain
        vo2_kidney = kidney_flow * cao2_frac * e_kidney
        vo2_liver = liver_flow * cao2_frac * e_liver
        vo2_muscle = muscle_flow * cao2_frac * e_muscle
        vo2_heart = heart_flow * cao2_frac * e_heart

        total = vo2_brain + vo2_kidney + vo2_liver + vo2_muscle + vo2_heart
        ratio = total / whole_body_fick
        sum_pass = 0.75 <= ratio <= 1.15

        ov1 = overlap(vo2_heart, vo2_heart, *anchor_direct_fick) if False else max(0.0, min(vo2_heart, anchor_direct_fick[1]) - max(vo2_heart, anchor_direct_fick[0]))
        # point value vo2_heart vs anchor RANGE: pass iff point falls inside either anchor band
        heart_pass = (anchor_direct_fick[0] <= vo2_heart <= anchor_direct_fick[1]) or (anchor_mechanical[0] <= vo2_heart <= anchor_mechanical[1])

        if sum_pass:
            pass_sum_only += 1
        if sum_pass and heart_pass:
            pass_combined += 1

    print(f"{label}: sum_only_pass_rate={100*pass_sum_only/N:.2f}%  combined_gate_pass_rate={100*pass_combined/N:.2f}%")
    return 100*pass_sum_only/N, 100*pass_combined/N

print("=== NULL 1: fully-uninformed wide ERO2 (each organ incl heart) Uniform[0.05,0.85] ===")
run_null((0.05,0.85),(0.05,0.85),(0.05,0.85),(0.05,0.85),(0.05,0.85), "null1_wide_uninformed")

print("\n=== NULL 2: heart-blind resting-tissue guess, ALL organs incl heart Uniform[0.15,0.45] ===")
run_null((0.15,0.45),(0.15,0.45),(0.15,0.45),(0.15,0.45),(0.15,0.45), "null2_heart_blind_resting")

print("\n=== SANITY: measured/correct bands (brain .30-.34, kidney .10-.16, liver .28-.32, muscle .24-.28, heart .70-.80) ===")
run_null((0.30,0.34),(0.10,0.16),(0.28,0.32),(0.24,0.28),(0.70,0.80), "measured_correct_bands")
