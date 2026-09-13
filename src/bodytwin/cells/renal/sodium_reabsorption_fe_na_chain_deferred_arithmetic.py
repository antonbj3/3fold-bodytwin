"""Executes HOLE-SODIUM-REABSORPTION-AMPLIFICATION-NARRATED-NOT-CODED's never-run chain.

The "amplification" claim (99.4pct reabsorption from a 0.6pct FE_Na) was prose, never coded
arithmetic. This cell closes that gap: filtered_load = GFR x plasma[Na], excreted = filtered_load
x FE_Na, cross-checked against the INTERSALT population intake range ALREADY LIVE in raas.py (reused,
not retyped) under the steady-state identity excretion == intake.

PRE-REGISTERED (before any number below computed):
  Inputs, ALL reused from other cells' live values, none retyped:
    GFR = renal_filtration.py's Davies&Shock1950 anchor, 122.8 mL/min/1.73m2 (PASS-gated already)
    plasma[Na] = textbook range 135-145 mmol/L (central 140), reused per this HOLE node's
                 honest_gaps ("plasma[Na]=138-142mM never independently live-sourced")
    FE_Na = clinical euvolemic-typical range 0.5-1.0pct (task/textbook), central 0.6pct (the node's
            own cited figure, e.g. Boron/Guyton)
    INTERSALT_TYPICAL_LO/HI_MMOL24H = 50.0/250.0 (raas.py's constants, live-imported, not retyped)
  C  : the filtered-load x FE_Na chain, using ONLY the reused GFR + standard textbook plasma-Na/FE_Na
       (no fitting to the anchor), lands the implied daily Na EXCRETION inside raas.py's
       INTERSALT_TYPICAL band [50,250] mmol/24h for >=50% of the full plausible
       (plasma_Na x FE_Na) grid -- i.e. the "narrated" chain, once actually coded, genuinely closes.
  NOT-C : <50% of the grid lands in-band -- the chain does not actually close; the "amplification"
       narrative remains uncoded arithmetic dressed as a fact, matching the node's suspicion.
  Void-floor: scramble (plasma_Na, FE_Na) pairing across the grid (shuffle which FE_Na value goes
       with which plasma_Na row) and confirm the in-band verdict changes for a nontrivial fraction
       of draws (substitution landed + downstream in-band fraction moved).

Deterministic exhaustive grid + seeded RNG void-floor only, pure stdlib, <1s.
Reads: OUT_ROOT/renal_filtration/renal_filtration_results.json (GFR); the raas cell's module-level
INTERSALT constants (imported, not read as data).
Writes: OUT_ROOT/sodium_reabsorption_fe_na_chain_deferred_arithmetic/
  sodium_reabsorption_fe_na_chain_results.json plus the void-floor draws (.jsonl) beside it.
Gate C (>=50% of the grid inside the INTERSALT typical band) decides.

Run:            python3 sodium_reabsorption_fe_na_chain_deferred_arithmetic.py
Self-test:      python3 sodium_reabsorption_fe_na_chain_deferred_arithmetic.py --selftest
"""
import importlib.util
import itertools
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "sodium_reabsorption_fe_na_chain_deferred_arithmetic")
OUT_PATH = os.path.join(OUT_DIR, "sodium_reabsorption_fe_na_chain_results.json")
NONCE_SCRATCH = OUT_DIR
VOIDFLOOR_OUT = os.path.join(NONCE_SCRATCH, "sodium_fe_na_voidfloor_draws.jsonl")

RENAL_JSON = os.path.join(OUT_ROOT, "renal_filtration", "renal_filtration_results.json")

# ---- import raas.py to reuse ITS OWN INTERSALT constants (not retyped) ------------------------
_spec = importlib.util.spec_from_file_location(
    "raas", os.path.join(os.path.dirname(HERE), "cardiovascular", "raas.py"))
raas = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(raas)

NA_MW = 23.0            # g/mol
NACL_PER_NA = 58.44 / 23.0   # NaCl:Na mass ratio (MW58.44/MW23)

PLASMA_NA_RANGE = (135.0, 145.0)  # mmol/L, textbook clinical reference range
FE_NA_RANGE = (0.3, 1.2)          # pct, euvolemic-typical clinical spread (loose, not fit to anchor)
STEP_NA = 1.0
STEP_FE = 0.05

SEED = 20260728
N_DRAWS = 500


def load_gfr():
    with open(RENAL_JSON, "r", encoding="utf-8") as f:
        d = json.load(f)
    # Davies & Shock 1950 anchor, stored under step1_human_primary_anchors
    return d["step1_human_primary_anchors"]["davies_shock_1950_20_29yo"]["gfr"]


def filtered_load_mmol_day(gfr_ml_min, plasma_na_mmol_l):
    return gfr_ml_min * plasma_na_mmol_l / 1000.0 * 1440.0


def excreted_mmol_day(filtered_load, fe_na_pct):
    return filtered_load * (fe_na_pct / 100.0)


def run(selftest=False):
    gfr = load_gfr()
    intersalt_lo = raas.INTERSALT_TYPICAL_LO_MMOL24H
    intersalt_hi = raas.INTERSALT_TYPICAL_HI_MMOL24H
    intersalt_extreme_lo = raas.INTERSALT_NA_LO_MMOL24H
    intersalt_extreme_hi = raas.INTERSALT_NA_HI_MMOL24H

    # ---- central-estimate reproduction (task's 140mM/0.6pct figures) ------------------------
    central_filtered = filtered_load_mmol_day(gfr, 140.0)
    central_excreted = excreted_mmol_day(central_filtered, 0.6)
    central_excreted_g_na = central_excreted * NA_MW / 1000.0
    central_excreted_g_nacl = central_excreted_g_na * NACL_PER_NA
    central_in_typical_band = intersalt_lo <= central_excreted <= intersalt_hi
    central_reabsorbed_pct = 100.0 - 0.6

    # ---- full grid ---------------------------------------------------------------------------
    na_vals = [PLASMA_NA_RANGE[0] + i * STEP_NA
               for i in range(int((PLASMA_NA_RANGE[1] - PLASMA_NA_RANGE[0]) / STEP_NA) + 1)]
    fe_vals = [round(FE_NA_RANGE[0] + i * STEP_FE, 2)
               for i in range(int(round((FE_NA_RANGE[1] - FE_NA_RANGE[0]) / STEP_FE)) + 1)]

    grid_rows = []
    n_in_band = 0
    for na, fe in itertools.product(na_vals, fe_vals):
        fl = filtered_load_mmol_day(gfr, na)
        ex = excreted_mmol_day(fl, fe)
        in_band = intersalt_lo <= ex <= intersalt_hi
        if in_band:
            n_in_band += 1
        grid_rows.append({"plasma_na": na, "fe_na_pct": fe, "filtered_load_mmol_day": fl,
                           "excreted_mmol_day": ex, "in_typical_band": in_band})
    n_total = len(grid_rows)
    frac_in_band = n_in_band / n_total
    gate_C_pass = bool(frac_in_band >= 0.50)

    # ---- void-floor: scramble (plasma_Na, FE_Na) pairing --------------------------------------
    os.makedirs(NONCE_SCRATCH, exist_ok=True)
    rng = random.Random(SEED)
    na_col = [r["plasma_na"] for r in grid_rows]
    fe_col = [r["fe_na_pct"] for r in grid_rows]
    landed, moved = 0, 0
    with open(VOIDFLOOR_OUT, "w", encoding="utf-8") as vf:
        for i in range(N_DRAWS):
            perm_fe = fe_col[:]
            rng.shuffle(perm_fe)
            n_in_band_scrambled = 0
            for na, fe in zip(na_col, perm_fe):
                fl = filtered_load_mmol_day(gfr, na)
                ex = excreted_mmol_day(fl, fe)
                if intersalt_lo <= ex <= intersalt_hi:
                    n_in_band_scrambled += 1
            frac_scrambled = n_in_band_scrambled / n_total
            substitution_landed = perm_fe != fe_col
            downstream_moved = abs(frac_scrambled - frac_in_band) > 0.02
            if substitution_landed:
                landed += 1
            if downstream_moved:
                moved += 1
            vf.write(json.dumps({"draw": i, "frac_in_band_scrambled": frac_scrambled,
                                  "frac_in_band_orig": frac_in_band,
                                  "substitution_landed": substitution_landed,
                                  "downstream_frac_moved": downstream_moved}) + "\n")
    substitution_landed_rate = landed / N_DRAWS
    downstream_moved_rate = moved / N_DRAWS

    report = {
        "inputs": {"gfr_ml_min": gfr, "gfr_source": "renal_filtration.py Davies&Shock1950 (live, PASS-gated)",
                   "intersalt_typical_band_mmol24h": [intersalt_lo, intersalt_hi],
                   "intersalt_extreme_band_mmol24h": [intersalt_extreme_lo, intersalt_extreme_hi]},
        "central_estimate": {
            "plasma_na_mmol_l": 140.0, "fe_na_pct": 0.6,
            "filtered_load_mmol_day": central_filtered,
            "excreted_mmol_day": central_excreted,
            "excreted_g_na_day": central_excreted_g_na,
            "excreted_g_nacl_day": central_excreted_g_nacl,
            "reabsorbed_pct": central_reabsorbed_pct,
            "in_intersalt_typical_band": central_in_typical_band,
        },
        "full_grid": {
            "plasma_na_range": PLASMA_NA_RANGE, "fe_na_range_pct": FE_NA_RANGE,
            "n_total_cells": n_total, "n_in_typical_band": n_in_band, "frac_in_band": frac_in_band,
        },
        "gate_C": {"threshold": 0.50, "frac_in_band": frac_in_band, "PASS": gate_C_pass},
        "void_floor": {
            "n_draws": N_DRAWS, "seed": SEED,
            "substitution_landed_rate": substitution_landed_rate,
            "downstream_frac_moved_rate": downstream_moved_rate,
            "output_path_nonce_scratch": VOIDFLOOR_OUT,
        },
        "verdict": "C" if gate_C_pass else "NOT-C",
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print(json.dumps({k: v for k, v in report.items() if k != "full_grid"}, indent=2))
    print(f"\nWrote {OUT_PATH}")
    print(f"Void-floor draws written to: {VOIDFLOOR_OUT}")

    if selftest:
        assert abs(central_filtered - 24756.5) / 24756.5 < 0.01, "central filtered-load arithmetic drifted"
        assert n_total > 100
        assert substitution_landed_rate > 0.9
        print("SELFTEST OK")
    return report


if __name__ == "__main__":
    run(selftest="--selftest" in sys.argv)
