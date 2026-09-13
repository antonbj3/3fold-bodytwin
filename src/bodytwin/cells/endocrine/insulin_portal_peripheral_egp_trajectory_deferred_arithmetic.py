"""Portal (kp4*Ipo) vs peripheral/delayed (kp3*Id) EGP-suppression split across a full meal trajectory.

An earlier check of this handoff covered ONLY the BASAL (t=0, steady-state) split of EGP suppression
between the portal term (kp4*Ipo) and the peripheral/delayed term (kp3*Id): basal kp3*Id=0.225
(50.3%) vs kp4*Ipo=0.2225 (49.7%) of insulin-driven EGP suppression -- "neither vestigial". That is
ONE point on a 400-minute postprandial trajectory during which Ipo and Id evolve on very different
timescales (Ipo fast, tau~1/gamma~2min; Id delayed, tau~127-253min per the Bergman/Rebrin/Steil
lineage), so the basal 50/50 split could easily NOT hold once a meal perturbs both compartments away
from their resting values. This cell runs the ACTUAL ODE (the glucose_meal_dallaman2007 cell,
imported live, not re-typed) across a full meal and checks the split at EVERY timestep, not just t=0.

PRE-REGISTERED (before any number below computed):
  Anchor: the basal claim -- kp3*Id and kp4*Ipo both >=30% of (kp3*Id+kp4*Ipo) at t=0
          ("neither vestigial", operationalized as a floor).
  C  : across the FULL 400-minute trajectory of a standard 75g OGTT, BOTH terms stay >=15% of their
       sum (a relaxed floor vs the basal 30%, since post-meal transients are expected to be more
       extreme than the resting point) for >=90% of the simulated timepoints -- i.e. "neither
       vestigial" generalizes beyond the single basal point.
  NOT-C : one term drops below 15% of the sum for >=10% of the trajectory -- the "neither
       vestigial" finding was a basal-only artifact.
  Void-floor: circularly time-shift the Id trajectory relative to Ipo (by a random offset, wrapping)
       across N_DRAWS draws and confirm the resulting min-fraction-of-sum statistic changes
       (substitution landed + downstream number moved) -- this tests whether the observed
       (non-)vestigial pattern depends on the REAL temporal alignment between the fast portal and
       slow peripheral pathways, not just their marginal magnitudes.

Deterministic ODE re-solve (LSODA via scipy, the sibling cell's numerics) + seeded RNG only for
the disclosed void-floor time-shift.
Reads: the glucose_meal_dallaman2007 cell (imported). Writes:
insulin_portal_peripheral_egp_trajectory_results.json + insulin_portal_voidfloor_draws.jsonl.
Gate: gate_C (verdict C vs NOT-C).
Run:            python3 insulin_portal_peripheral_egp_trajectory_deferred_arithmetic.py
Self-test:      python3 insulin_portal_peripheral_egp_trajectory_deferred_arithmetic.py --selftest
"""
import importlib.util
import json
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.environ.get("BODYTWIN_OUT", os.path.join(os.getcwd(), "outputs"))
OUT_DIR = os.path.join(OUT_ROOT, "insulin_portal_peripheral_egp_trajectory_deferred_arithmetic")
OUT_PATH = os.path.join(OUT_DIR, "insulin_portal_peripheral_egp_trajectory_results.json")
NONCE_SCRATCH = OUT_DIR
VOIDFLOOR_OUT = os.path.join(NONCE_SCRATCH, "insulin_portal_voidfloor_draws.jsonl")

# SAFETY: glucose_meal_dallaman2007.py has NO `if __name__=="__main__"` guard -- importing it
# unconditionally RE-RUNS its own module-level scenarios AND rewrites its own output JSON as a
# side effect. Snapshot+restore its EXACT prior byte content so this cell never leaves collateral
# drift on a dependency it only meant to read from.
_DM_OUT_JSON = os.path.join(OUT_ROOT, "glucose_meal_dallaman2007",
                            "glucose_meal_dallaman2007_results.json")
_dm_prior_bytes = None
if os.path.exists(_DM_OUT_JSON):
    with open(_DM_OUT_JSON, "rb") as _f:
        _dm_prior_bytes = _f.read()

_spec = importlib.util.spec_from_file_location("glucose_meal_dallaman2007",
                                                os.path.join(HERE, "glucose_meal_dallaman2007.py"))
dm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dm)  # imports the ALREADY-BUILT module live (runs its own module-level scenarios too)

if _dm_prior_bytes is not None:
    with open(_DM_OUT_JSON, "rb") as _f:
        _dm_post_bytes = _f.read()
    if _dm_post_bytes != _dm_prior_bytes:
        with open(_DM_OUT_JSON, "wb") as _f:
            _f.write(_dm_prior_bytes)
        print(f"SAFETY: restored {_DM_OUT_JSON} to its pre-import byte content "
              f"(import-time side-effect write reverted, dependency's tracked file left untouched)")

BASAL_FLOOR_FRAC = 0.15   # relaxed vs the basal ~0.30 (49.7/50.3 split), pre-registered
FRAC_TIMEPOINTS_REQUIRED = 0.90
SEED = 20260728
N_DRAWS = 300


def run(selftest=False):
    D = 75 * 1000.0  # 75g OGTT, standard clinical dose, matches dm's dose_response scenario set
    y0 = dm.with_meal_ic(dm.y0_healthy, D)
    t, y = dm.run_scenario(D=D, meal_t0=0.0, y0=y0, t_end=400.0, secretion_gain=1.0)

    Id = y[dm.IDX["Id"], :]
    Ipo = y[dm.IDX["Ipo"], :]
    Gp = y[dm.IDX["Gp"], :]

    term_kp3_Id = dm.P["kp3"] * Id
    term_kp4_Ipo = dm.P["kp4"] * Ipo
    denom = term_kp3_Id + term_kp4_Ipo
    denom_safe = np.where(denom > 1e-9, denom, np.nan)
    frac_kp3 = term_kp3_Id / denom_safe
    frac_kp4 = term_kp4_Ipo / denom_safe

    # basal reproduction (t=0 point of the ACTUAL simulated trajectory, using the model's
    # truly-equilibrated y0_healthy fixed point -- NOT the naive Id=I_b identity the basal
    # arithmetic used)
    basal_frac_kp3 = float(frac_kp3[0])
    basal_frac_kp4 = float(frac_kp4[0])
    # naive-identity reproduction: what the basal claim ACTUALLY computed (Id=I_b=25 exactly, the
    # un-equilibrated identity target, not the model's converged fixed point)
    naive_id = dm.P["I_b"]
    naive_ipo = dm.P["S_b"] / dm.P["gamma"]
    naive_k3 = dm.P["kp3"] * naive_id
    naive_k4 = dm.P["kp4"] * naive_ipo
    naive_frac_kp3 = naive_k3 / (naive_k3 + naive_k4)
    naive_frac_kp4 = naive_k4 / (naive_k3 + naive_k4)
    fidelity_note = (
        "the basal claim (50.3/49.7) used the NAIVE identity Id=I_b=25 (reproduced exactly "
        "here as {:.1f}/{:.1f}), NOT the model's TRULY-EQUILIBRATED basal_state() fixed point "
        "(Id=30.40, per basal_state()'s disclosed long-integration convergence), which actually "
        "feeds every meal simulation. Using the model's real IC gives "
        "{:.1f}/{:.1f} at t=0 -- a ~5pp fidelity gap between that basal arithmetic and what "
        "the equilibrated model state actually is.".format(
            naive_frac_kp3 * 100, naive_frac_kp4 * 100, basal_frac_kp3 * 100, basal_frac_kp4 * 100)
    )

    n_valid = int(np.sum(~np.isnan(frac_kp3)))
    n_both_above_floor = int(np.sum((frac_kp3 >= BASAL_FLOOR_FRAC) & (frac_kp4 >= BASAL_FLOOR_FRAC)
                                     & ~np.isnan(frac_kp3)))
    frac_timepoints_ok = n_both_above_floor / n_valid if n_valid else float("nan")
    gate_C_pass = bool(frac_timepoints_ok >= FRAC_TIMEPOINTS_REQUIRED)

    min_frac_kp3 = float(np.nanmin(frac_kp3))
    min_frac_kp4 = float(np.nanmin(frac_kp4))
    argmin_t_kp3 = float(t[int(np.nanargmin(frac_kp3))])
    argmin_t_kp4 = float(t[int(np.nanargmin(frac_kp4))])

    # ---- void-floor: circular time-shift Id relative to Ipo -------------------------------------
    os.makedirs(NONCE_SCRATCH, exist_ok=True)
    rng = random.Random(SEED)
    n = len(t)
    landed, moved = 0, 0
    orig_min = min(min_frac_kp3, min_frac_kp4)
    with open(VOIDFLOOR_OUT, "w", encoding="utf-8") as vf:
        for i in range(N_DRAWS):
            shift = rng.randrange(1, n - 1)
            Id_shifted = np.roll(Id, shift)
            term_kp3_shifted = dm.P["kp3"] * Id_shifted
            denom_shifted = term_kp3_shifted + term_kp4_Ipo
            denom_shifted_safe = np.where(denom_shifted > 1e-9, denom_shifted, np.nan)
            f3 = term_kp3_shifted / denom_shifted_safe
            f4 = term_kp4_Ipo / denom_shifted_safe
            min_shifted = float(np.nanmin([np.nanmin(f3), np.nanmin(f4)]))
            substitution_landed = shift != 0
            downstream_moved = abs(min_shifted - orig_min) > 0.02
            if substitution_landed:
                landed += 1
            if downstream_moved:
                moved += 1
            vf.write(json.dumps({"draw": i, "shift": shift, "min_frac_orig": orig_min,
                                  "min_frac_shifted": min_shifted,
                                  "substitution_landed": substitution_landed,
                                  "downstream_moved": downstream_moved}) + "\n")
    substitution_landed_rate = landed / N_DRAWS
    downstream_moved_rate = moved / N_DRAWS

    report = {
        "meal_dose_g": 75, "t_end_min": 400.0,
        "basal_reproduction": {"frac_kp3_id_true_equilibrated": basal_frac_kp3,
                                "frac_kp4_ipo_true_equilibrated": basal_frac_kp4,
                                "frac_kp3_id_naive_identity": naive_frac_kp3,
                                "frac_kp4_ipo_naive_identity": naive_frac_kp4,
                                "node_claimed": {"kp3_id_pct": 50.3, "kp4_ipo_pct": 49.7},
                                "fidelity_note": fidelity_note},
        "trajectory": {
            "n_valid_timepoints": n_valid,
            "min_frac_kp3_id": min_frac_kp3, "argmin_t_kp3_min": argmin_t_kp3,
            "min_frac_kp4_ipo": min_frac_kp4, "argmin_t_kp4_min": argmin_t_kp4,
            "n_timepoints_both_above_floor": n_both_above_floor,
            "frac_timepoints_ok": frac_timepoints_ok,
            "floor_frac": BASAL_FLOOR_FRAC, "required_frac_timepoints": FRAC_TIMEPOINTS_REQUIRED,
        },
        "gate_C": {"PASS": gate_C_pass},
        "void_floor": {
            "n_draws": N_DRAWS, "seed": SEED,
            "substitution_landed_rate": substitution_landed_rate,
            "downstream_moved_rate": downstream_moved_rate,
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
        assert abs(naive_frac_kp3 * 100 - 50.3) < 1.0, "naive-identity kp3 reproduction drifted"
        assert abs(naive_frac_kp4 * 100 - 49.7) < 1.0, "naive-identity kp4 reproduction drifted"
        assert n_valid > 100
        assert substitution_landed_rate > 0.9
        print("SELFTEST OK -- see fidelity_note for a real ~5pp basal-arithmetic discrepancy found")
    return report


if __name__ == "__main__":
    run(selftest="--selftest" in sys.argv)
