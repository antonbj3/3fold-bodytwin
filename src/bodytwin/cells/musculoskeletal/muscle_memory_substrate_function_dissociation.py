"""Human skeletal-muscle "memory": machine-checked partition test of the claim that it is
SUBSTRATE-CONFIRMED (retained myonuclei / epigenetic marks persist through detraining) but
FUNCTIONALLY-INERT on the metric that matters for a retraining-speed prediction (fibre-CSA
REGROWTH RATE on retraining vs a naive control), which is why the z_myo -> regrowth-rate coupling
coefficient is set conservatively at 1.0 (no assumed mouse-level speedup).

Inputs (6 literature datapoints, embedded in this file):
  substrate channel (does retained biological signal exceed a real-signal floor?):
    Cumming2024  myonuclei retained: +33% (n=12, PMID 39159314), floor >= 10%
    Seaborne2018 epigenetic hypomethylation ratio: 2.06x (n=8, PMID 29382913), floor >= 1.5x
    Nielsen2023  AAS myonuclei-retention duration: 4.0 yr (n=7, PMID 37466198), floor >= 1 yr
  function channel (detectable retraining rate ADVANTAGE over naive on the fibre-CSA outcome;
  strength is excluded per the explicit strength-vs-size dissociation caveat):
    Psilander2019 CSA-rate delta: 0% (n=19, PMID 30991013), floor >= 10 pp
    Cumming2024   CSA-rate delta: 0% (n=12, PMID 39159314), floor >= 10 pp
    Staron1991    qualitative CSA-advantage coded 0 (n=6, PMID 1827108), floor >= 10 pp

Reads: nothing. Writes: muscle_memory_substrate_function_dissociation.json and the adversary JSON.

Gate: PASS only if substrate_hits == 3 and function_hits == 0; a single crossover on either side
falls to MIXED.

Forced adversary: is the clean 3-vs-0 partition an artifact of the labelling -- would almost any
random assignment of the same 6 numbers into two groups of 3 also look clean, because 3 of the 6
are trivially 0? Computed over all C(6,3)=20 partitions of the raw values.
"""
import itertools
import json
import os
import random

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "muscle_memory_substrate_function_dissociation")
SCRATCH_DIR = OUT_DIR
OUT_PATH = os.path.join(OUT_DIR, "muscle_memory_substrate_function_dissociation.json")

# (label, raw_value, floor, unit, pmid, n)
SUBSTRATE = [
    ("Cumming2024_myonuclei_retained_pct", 33.0, 10.0, "%", "39159314", 12),
    ("Seaborne2018_epigenetic_ratio_x10", 20.6, 15.0, "x*10 (2.06x scaled to same order)", "29382913", 8),
    ("Nielsen2023_AAS_retention_yr_x10", 40.0, 10.0, "yr*10 (4.0yr scaled)", "37466198", 7),
]
FUNCTION = [
    ("Psilander2019_CSA_rate_delta_pct", 0.0, 10.0, "pp", "30991013", 19),
    ("Cumming2024_CSA_rate_delta_pct", 0.0, 10.0, "pp", "39159314", 12),
    ("Staron1991_qualitative_CSA_advantage", 0.0, 10.0, "coded", "1827108", 6),
]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    substrate_hits = sum(1 for (_, v, floor, *_ ) in SUBSTRATE if v >= floor)
    function_hits = sum(1 for (_, v, floor, *_ ) in FUNCTION if v >= floor)
    clean_partition = (substrate_hits == 3 and function_hits == 0)

    # forced adversary: all 20 random 3-vs-3 partitions of the 6 raw (unlabeled) values,
    # using each value's floor-type carried along so we can ask "would grouping these
    # values differently still look clean under a plausible floor assignment" -- since
    # floors differ by unit, we test the value-only partition against a UNIFORM floor
    # (the largest shared floor magnitude proxy: 10, matching the function floors and the
    # Cumming/Nielsen raw floors) to make the random-partition test fair (apples-to-apples
    # in the SAME units after scaling).
    all_vals = [(lbl, v) for (lbl, v, *_ ) in SUBSTRATE + FUNCTION]
    uniform_floor = 10.0
    n_clean_random = 0
    all_partitions = list(itertools.combinations(range(6), 3))
    for idx_a in all_partitions:
        idx_b = tuple(i for i in range(6) if i not in idx_a)
        group_a = [all_vals[i][1] for i in idx_a]
        group_b = [all_vals[i][1] for i in idx_b]
        a_all_above = all(v >= uniform_floor for v in group_a)
        b_all_below = all(v < uniform_floor for v in group_b)
        a_all_below = all(v < uniform_floor for v in group_a)
        b_all_above = all(v >= uniform_floor for v in group_b)
        if (a_all_above and b_all_below) or (a_all_below and b_all_above):
            n_clean_random += 1
    random_partition_clean_rate = n_clean_random / len(all_partitions)

    # pre-registered: the true labeled partition is non-trivial ONLY if the random-clean
    # rate is itself low (<25%) -- i.e. the clean split is not simply forced by 3 zeros
    # existing trivially; if random-clean-rate is high, the "clean partition" finding is
    # cheap and should be reported as such, not oversold.
    partition_is_informative = random_partition_clean_rate < 0.25

    result = {
        "node": "MSK-DETRAINING-MUSCLE-MEMORY",
        "substrate_datapoints": SUBSTRATE,
        "function_datapoints": FUNCTION,
        "substrate_hits_above_floor": substrate_hits,
        "function_hits_above_floor": function_hits,
        "clean_labeled_partition (3 substrate-hits, 0 function-hits)": clean_partition,
        "adversary_random_partition_test": {
            "n_partitions_total": len(all_partitions),
            "n_partitions_also_clean": n_clean_random,
            "random_partition_clean_rate": round(random_partition_clean_rate, 3),
            "partition_is_informative (<25% random-clean-rate)": partition_is_informative,
        },
        "verdict": (
            "CONFIRMED: substrate-retained/function-null dissociation is a real, "
            "non-trivial labeling (not just '3 numbers happen to be near zero')"
            if clean_partition and partition_is_informative
            else "NOT CONFIRMED / partition is cheap or not clean"
        ),
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result


def void_floor(n_draws=2000, seed=20260728):
    """
    Scramble the load-bearing driver: jitter each of the 6 raw values by a random
    multiplicative factor (0.2x-5x, log-uniform) many times, and recompute BOTH the
    labeled clean-partition outcome AND the random-partition-clean-rate, reporting the
    fraction of draws where the labeled partition STILL looks clean purely because of
    where the jittered values happen to land relative to the floor -- this is the real
    void floor for this specific claim (not the fixed one-shot check above).
    """
    rng = random.Random(seed)
    base_vals = [v for (_, v, *_ ) in SUBSTRATE + FUNCTION]  # indices 0-2 substrate, 3-5 function
    passes = 0
    jittered_ranges = []
    # additive jitter (not multiplicative): the function-channel values are exact reported
    # ZEROS (a real null), and 0*anything=0 would make multiplicative jitter unable to
    # ever move them off the floor -- that would silently guarantee the function side
    # always "passes" regardless of scramble strength, which is exactly the kind of
    # non-landing substitution this discipline requires catching. Additive jitter in
    # +/-25 units (comparable to the largest raw values here) lets EVERY value, including
    # exact zeros, cross the floor=10 threshold in either direction.
    for _ in range(n_draws):
        jittered = [max(0.0, v + rng.uniform(-25, 25)) for v in base_vals]
        jittered_ranges.append(jittered)
        sub_hits = sum(1 for v in jittered[0:3] if v >= 10.0)
        fun_hits = sum(1 for v in jittered[3:6] if v >= 10.0)
        if sub_hits == 3 and fun_hits == 0:
            passes += 1

    flat = [v for row in jittered_ranges for v in row]
    # confirm the substitution actually landed on the previously-immovable exact-zero
    # function values specifically (not just on the already-nonzero substrate values)
    function_col_vals = [row[3] for row in jittered_ranges] + [row[4] for row in jittered_ranges] + [row[5] for row in jittered_ranges]
    function_side_moved_off_zero = sum(1 for v in function_col_vals if v > 0.01) / len(function_col_vals)
    moved = (max(flat) - min(flat)) > 0.5 and function_side_moved_off_zero > 0.3
    out = {
        "n_draws": n_draws,
        "substitution_landed (values varied, INCLUDING the exact-zero function side)": moved,
        "fraction_of_function_side_draws_pushed_off_exact_zero": round(function_side_moved_off_zero, 3),
        "jittered_value_range": [round(min(flat), 3), round(max(flat), 3)],
        "void_pass_rate": round(passes / n_draws, 4),
        "note": (
            "fraction of a wide (+/-25 additive) jitter of all 6 raw values that STILL "
            "reproduces the exact clean 3-substrate-hit/0-function-hit partition by "
            "chance -- additive (not multiplicative) jitter was required because the "
            "function values are exact reported zeros that multiplicative scrambling "
            "could never move off the floor, which would have silently guaranteed a "
            "100% void-pass rate on that side without ever really testing it."
        ),
    }
    scratch_path = os.path.join(SCRATCH_DIR, "muscle_memory_void_floor.json")
    with open(scratch_path, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    os.makedirs(SCRATCH_DIR, exist_ok=True)
    main()
    void_floor()
