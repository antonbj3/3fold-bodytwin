"""Deferred full-grid version of
HOLE-RENAL-NEPHRON-SEGMENT-FRACTION-CROSSMODEL-DISAGREEMENT-CONTRADICTION's "segment-sum
forced-adversary" leg, which only ran a 1-D sweep: "hold TAL25/DCT5/CD2.5 fixed, sweep
PCT{55,67,70} vs 99.4pct anchor". That leaves the TAL-must-compensate finding untested against the
alternative that DCT/CD (not TAL) absorb the slack -- a genuinely different, unexecuted question.

PRE-REGISTERED (before any grid cell computed):
  Anchor: total fractional Na reabsorption must land in [99.0, 99.8]% (FE_Na 0.2-1.0%, standard
          euvolemic clinical range; anchor independent of this script).
  C  : across the FULL joint grid of segment fractions consistent with independent per-segment
       literature ranges, the anchor-satisfying subset requires TAL > 28% (i.e. materially above
       the classical ~25% textbook figure) in a CLEAR MAJORITY (>=80%) of anchor-passing grid
       cells -- i.e. the 1-D finding ("TAL must compensate, ~36.9%") is a robust, structural
       feature of the joint space, not an artifact of holding DCT/CD fixed at their 1-D values.
  NOT-C: if a large minority/majority (>20%) of anchor-passing cells achieve the anchor via
       classical TAL (<=28%) compensated instead by DCT/CD/PCT variation, the "TAL-must-compensate"
       framing overstates TAL's uniqueness -- DCT/CD are equally-valid escape routes the 1-D sweep
       never tested because it held them fixed.
  Per-segment literature ranges (independent, not fit to the anchor):
    PCT  : 50-70%  (Boron/Guyton textbook spread; Kiil&Sejersted 2003 vol-expanded low end ~50-55%)
    TAL  : 15-40%  (25% classical; Kiil&Sejersted 2003 PMID12713517 dog vol-expanded up to 40%)
    DCT  : 3-8%    (Boron/Guyton ~5% central)
    CD   : 1-4%    (Boron/Guyton ~2.5% central, final fine-tuning segment)
  Void-floor: for each anchor-passing cell, additionally draw N_DRAWS scrambled re-assignments of
    which segment carries "above-classical" share (permute the four fraction values across the
    four segment LABELS while keeping the total fixed) -- confirms the substitution actually lands
    (label reassigned, sum unchanged) and moves a downstream number (whether TAL>28pct becomes
    PCT>72pct, etc).

Deterministic exhaustive grid (no RNG needed for the main gate; RNG only for the disclosed
void-floor label-permutation, seeded for reproducibility). Pure stdlib+numpy, <2s.
Reads: nothing (all inputs are literals).
Writes: OUT_ROOT/renal_nephron_segment_fraction_grid_deferred_arithmetic/
  renal_nephron_segment_fraction_grid_results.json plus the void-floor draws (.jsonl) beside it.
Gate C (>=80% of anchor-passing cells require TAL>28%) decides.

Run:            python3 renal_nephron_segment_fraction_grid_deferred_arithmetic.py
Self-test:      python3 renal_nephron_segment_fraction_grid_deferred_arithmetic.py --selftest
"""
import itertools
import json
import os
import random
import sys

OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "renal_nephron_segment_fraction_grid_deferred_arithmetic")
OUT_PATH = os.path.join(OUT_DIR, "renal_nephron_segment_fraction_grid_results.json")
NONCE_SCRATCH = OUT_DIR
VOIDFLOOR_OUT = os.path.join(NONCE_SCRATCH, "renal_segment_voidfloor_draws.jsonl")

ANCHOR_LO, ANCHOR_HI = 99.0, 99.8   # pct, FE_Na 0.2-1.0pct euvolemic clinical range

PCT_RANGE = (50.0, 70.0)
TAL_RANGE = (15.0, 40.0)   # Kiil & Sejersted 2003 PMID 12713517 up to 40pct (dog, vol-expanded)
DCT_RANGE = (3.0, 8.0)
CD_RANGE = (1.0, 4.0)

TAL_CLASSICAL_CEILING = 28.0   # "materially above classical ~25pct" threshold, pre-registered
STEP = 1.0                     # pct, grid resolution
N_DRAWS = 500
SEED = 20260728


def build_grid():
    pct_vals = [PCT_RANGE[0] + i * STEP for i in range(int((PCT_RANGE[1] - PCT_RANGE[0]) / STEP) + 1)]
    tal_vals = [TAL_RANGE[0] + i * STEP for i in range(int((TAL_RANGE[1] - TAL_RANGE[0]) / STEP) + 1)]
    dct_vals = [DCT_RANGE[0] + i * STEP for i in range(int((DCT_RANGE[1] - DCT_RANGE[0]) / STEP) + 1)]
    cd_vals = [CD_RANGE[0] + i * STEP for i in range(int((CD_RANGE[1] - CD_RANGE[0]) / STEP) + 1)]
    return pct_vals, tal_vals, dct_vals, cd_vals


def run(selftest=False):
    pct_vals, tal_vals, dct_vals, cd_vals = build_grid()
    total_cells = 0
    passing = []
    for pct, tal, dct, cd in itertools.product(pct_vals, tal_vals, dct_vals, cd_vals):
        total_cells += 1
        s = pct + tal + dct + cd
        if ANCHOR_LO <= s <= ANCHOR_HI:
            passing.append((pct, tal, dct, cd, s))

    n_pass = len(passing)
    n_tal_above = sum(1 for p in passing if p[1] > TAL_CLASSICAL_CEILING)
    frac_tal_above = n_tal_above / n_pass if n_pass else float("nan")

    # symmetric check: also test whether DCT or CD alone (not TAL) can equally carry the "excess"
    n_dct_above_central = sum(1 for p in passing if p[2] > 6.0)   # >6pct = above its own central ~5pct
    n_cd_above_central = sum(1 for p in passing if p[3] > 2.5)     # >2.5pct = above its own central
    n_pct_below_central = sum(1 for p in passing if p[0] < 65.0)  # below typical 65-67pct central estimate

    frac_dct_above_central = n_dct_above_central / n_pass if n_pass else float("nan")
    frac_cd_above_central = n_cd_above_central / n_pass if n_pass else float("nan")
    frac_pct_below_central = n_pct_below_central / n_pass if n_pass else float("nan")

    gate_C_pass = bool(n_pass > 0 and frac_tal_above >= 0.80)

    # 1-D reproduction check (node's original claim): hold TAL=25,DCT=5,CD=2.5 fixed, PCT swept
    orig_sweep = {}
    for pct_test in (55.0, 67.0, 70.0):
        s = pct_test + 25.0 + 5.0 + 2.5
        orig_sweep[pct_test] = {"total": s, "in_anchor": ANCHOR_LO <= s <= ANCHOR_HI}
    # TAL-need to rescue a PCT=55 scenario (node's stated "50-60pct-PCT scenario", midpoint 55)
    # to 99.4 with DCT/CD fixed at 5/2.5 -- the node's cited 36.9pct arithmetic
    tal_need_for_pct55 = 99.4 - 55.0 - 5.0 - 2.5

    # ---------------- void-floor: permute which LABEL gets which passing-cell's fraction --------
    os.makedirs(NONCE_SCRATCH, exist_ok=True)
    rng = random.Random(SEED)
    landed = 0
    moved = 0
    with open(VOIDFLOOR_OUT, "w", encoding="utf-8") as vf:
        sample_source = passing if passing else [(60.0, 25.0, 5.0, 2.5, 92.5)]
        for i in range(N_DRAWS):
            cell = sample_source[rng.randrange(len(sample_source))]
            vals = list(cell[:4])
            perm = vals[:]
            rng.shuffle(perm)
            substitution_landed = (perm != vals) and (abs(sum(perm) - sum(vals)) < 1e-9)
            # downstream number: does TAL-label now carry a value that flips the >28 classification?
            orig_tal_above = vals[1] > TAL_CLASSICAL_CEILING
            new_tal_above = perm[1] > TAL_CLASSICAL_CEILING
            downstream_moved = orig_tal_above != new_tal_above
            if perm != vals:
                landed += 1
            if downstream_moved:
                moved += 1
            vf.write(json.dumps({
                "draw": i, "orig": vals, "perm": perm,
                "substitution_landed": substitution_landed,
                "downstream_tal_classification_moved": downstream_moved,
            }) + "\n")
    substitution_landed_rate = landed / N_DRAWS
    downstream_moved_rate = moved / N_DRAWS

    report = {
        "anchor_band_pct": [ANCHOR_LO, ANCHOR_HI],
        "grid": {
            "pct_range": PCT_RANGE, "tal_range": TAL_RANGE, "dct_range": DCT_RANGE, "cd_range": CD_RANGE,
            "step_pct": STEP, "total_cells": total_cells, "n_anchor_passing": n_pass,
        },
        "gate_C": {
            "description": "frac of anchor-passing cells requiring TAL>28pct",
            "frac_tal_above_28pct": frac_tal_above,
            "threshold": 0.80,
            "PASS": gate_C_pass,
        },
        "symmetric_escape_routes": {
            "frac_dct_above_central_6pct": frac_dct_above_central,
            "frac_cd_above_central_2_5pct": frac_cd_above_central,
            "frac_pct_below_central_65pct": frac_pct_below_central,
            "interpretation": (
                "if these fractions are ALSO high, DCT/CD/PCT are equally-valid alternative escape "
                "routes and TAL is not uniquely privileged by the anchor alone -- the original 1-D "
                "sweep's TAL-specific framing came from holding DCT/CD fixed by construction, not "
                "from the anchor forcing TAL specifically."
            ),
        },
        "original_1d_sweep_reproduction": {
            "fixed_TAL25_DCT5_CD2_5_sweep_PCT": orig_sweep,
            "tal_need_for_pct55_scenario": tal_need_for_pct55,
            "node_claimed_tal_need": 36.9,
            "reproduction_diff_pp": abs(tal_need_for_pct55 - 36.9),
        },
        "void_floor": {
            "n_draws": N_DRAWS, "seed": SEED,
            "substitution_landed_rate": substitution_landed_rate,
            "downstream_tal_classification_moved_rate": downstream_moved_rate,
            "output_path_nonce_scratch": VOIDFLOOR_OUT,
        },
        "verdict": "C" if gate_C_pass else "NOT-C",
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(json.dumps(report, indent=2))
    print(f"\nWrote {OUT_PATH}")
    print(f"Void-floor draws written to: {VOIDFLOOR_OUT}")

    if selftest:
        assert total_cells > 1000, "grid too coarse"
        assert n_pass >= 0
        assert substitution_landed_rate > 0.9, "void-floor permutation barely ever changes labels -- broken shuffle"
        print("SELFTEST OK")
    return report


if __name__ == "__main__":
    run(selftest="--selftest" in sys.argv)
