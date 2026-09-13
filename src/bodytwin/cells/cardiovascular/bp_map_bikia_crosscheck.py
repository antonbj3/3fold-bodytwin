"""
Cross-model mean-arterial-pressure check against a held-out cohort anchor.

Claim under test: 3 MAP models (values in [84,93.33] mmHg) sit within |z|<2 of an
independently-measured, held-out anchor -- Bikia et al 2024 (Sci Rep 14:5913, PMID 38467721,
PMCID PMC10928153, Table 2), real in-vivo Asklepios cohort n=2263, "Mean blood pressure" =
100 +/- 12 mmHg -- never used to calibrate any of the 3 models.

PRE-REGISTERED (before running):
  Reuse the arterial_pressure cell's already-computed MAP values (do not recompute the
  physiology fresh -- the gap here is the missing cross-check step, not the MAP model itself).
  Take RAP=0 (the low end of that cell's 0-8 mmHg RAP sweep):
    model_1 = classic_map (DBP+PP/3)                = 93.3333 mmHg
    model_2 = CO*TPR Ohm's law, CO_geometric route   = 83.3696 mmHg  (rap_0/map_at_tpr1200)
    model_3 = CO*TPR Ohm's law, CO_subject route     = 84.8024 mmHg  (rap_0/map_at_tpr1200)
  GATE: |100 - model_i| / 12 < 2.0 for all 3 (Bikia's reported mean+/-SD, n=2263).
  FALSIFIER: any |z| >= 2 -> the "all 3 models land near a real independent cohort" claim
  does not hold for that model.

FORCED ADVERSARY: a wide +/-12 mmHg SD anchor could pass almost any physiologically-sane MAP
value (tautology-risk / a wide gate is not discriminating). Force it: void-floor a uniform
sweep of MAP over the full plausible resting-adult physiological envelope [55,145] mmHg (a
much wider envelope than any model or textbook resting range) and measure what FRACTION of
that envelope ALSO clears |z|<2 by chance alone. If a large majority of the envelope passes,
the gate is weak/uninformative (report, don't hide); if a minority passes, landing there is
informative. A second, independent structural adversary: the 3 real models are NOT a random
draw from that envelope -- they come from a physiology-constrained band; report both.

Reads:  the arterial_pressure cell result JSON.
Writes: bp_map_bikia_crosscheck.json and bp_map_void_floor.json under the cell output directory.
Gate:   |z| < 2 against the Bikia 2024 anchor for all 3 models.
"""
import json
import math
import os
import random

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "bp_map_bikia_crosscheck")
SCRATCH_DIR = OUT_DIR
OUT_PATH = _os.path.join(OUT_DIR, "bp_map_bikia_crosscheck.json")
ARTERIAL_RESULTS = _os.path.join(OUT_ROOT, "arterial_pressure", "arterial_pressure_results.json")

BIKIA_MEAN, BIKIA_SD, BIKIA_N = 100.0, 12.0, 2263  # PMID38467721 Table 2, Asklepios cohort
ENVELOPE_LO, ENVELOPE_HI = 55.0, 145.0  # mmHg, wide resting-adult physiological envelope


def ensure_arterial_results():
    if not os.path.exists(ARTERIAL_RESULTS):
        raise SystemExit(f"missing {ARTERIAL_RESULTS} -- run the arterial_pressure cell first")
    with open(ARTERIAL_RESULTS) as f:
        return json.load(f)


def extract_models(d):
    approx = d["map_approximation"]
    sweep = d["tpr_forward_backward"]["forward_sweep"]
    return {
        "classic_map (DBP+PP/3)": approx["classic_map_mmhg"],
        "CO*TPR (CO_geometric, RAP=0, TPR=1200)": sweep["co_geometric"]["rap_0"]["map_at_tpr1200"],
        "CO*TPR (CO_subject, RAP=0, TPR=1200)": sweep["co_subject"]["rap_0"]["map_at_tpr1200"],
    }


def zscore(x):
    return (BIKIA_MEAN - x) / BIKIA_SD


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SCRATCH_DIR, exist_ok=True)
    art = ensure_arterial_results()
    models = extract_models(art)

    per_model = {}
    all_pass = True
    for name, val in models.items():
        z = zscore(val)
        p = abs(z) < 2.0
        all_pass = all_pass and p
        per_model[name] = {"map_mmhg": round(val, 4), "z_vs_bikia": round(z, 4), "PASS_|z|<2": p}

    result = {
        "node": "cross-model MAP check vs a held-out cohort anchor",
        "anchor": {"source": "Bikia2024 PMID38467721 Table2, n=2263 Asklepios",
                   "MAP_mean_sd_mmhg": [BIKIA_MEAN, BIKIA_SD]},
        "models_from_arterial_pressure_py": per_model,
        "z_range_measured": [round(min(v["z_vs_bikia"] for v in per_model.values()), 3),
                             round(max(v["z_vs_bikia"] for v in per_model.values()), 3)],
        "node_claimed_z_range": [0.56, 1.33],
        "ALL_3_MODELS_PASS": all_pass,
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result, models


def void_floor(models, n_draws=2000, seed=20260728):
    rng = random.Random(seed)
    passes = 0
    draws = []
    for _ in range(n_draws):
        x = rng.uniform(ENVELOPE_LO, ENVELOPE_HI)
        draws.append(x)
        if abs(zscore(x)) < 2.0:
            passes += 1
    envelope_pass_rate = passes / n_draws

    moved = (max(draws) - min(draws)) > 1.0
    out = {
        "n_draws": n_draws,
        "envelope_mmhg": [ENVELOPE_LO, ENVELOPE_HI],
        "substitution_landed (uniform draw spans envelope)": moved,
        "draw_range": [round(min(draws), 2), round(max(draws), 2)],
        "envelope_void_pass_rate": round(envelope_pass_rate, 4),
        "real_models_z": {k: round(v, 4) for k, v in
                           {kk: zscore(vv) for kk, vv in models.items()}.items()},
        "note": ("fraction of a UNIFORM draw over the full wide resting-adult MAP envelope "
                 "[55,145] mmHg that would ALSO clear |z|<2 by chance -- measures whether the "
                 "Bikia SD=12 gate is discriminating or a near-tautology; real models are NOT "
                 "drawn from this envelope (physiology-constrained), reported separately."),
    }
    scratch_path = os.path.join(SCRATCH_DIR, "bp_map_void_floor.json")
    with open(scratch_path, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    result, models = main()
    void_floor(models)
