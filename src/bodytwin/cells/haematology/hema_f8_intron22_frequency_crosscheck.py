"""F8 intron-22 inversion frequency -- cross-replication check between two independent cohorts.

Computes the replication ratio between Antonarakis 1995 (PMID 7662970, n=2093 severe hemophilia A,
22 labs / 14 countries, 740 type-1 + 140 type-2 inversions = 880/2093 = 42.04%) and the earlier
Lakich 1993 discovery-era estimate (~45%, PMID 8275087, point value only -- no exact n or CI is
recoverable from the public abstract, so it is not treated as more precise than it is). A Wilson
95% CI is reported for the Antonarakis proportion so the point-vs-point comparison is not
overstated.

Forced adversary: inversion frequency could differ by ascertainment era (small early cohorts
overestimate recurrent lesions via referral/founder bias) or by population (the discovery cohort
was single-country). The gate is therefore a magnitude-convergence band, not identity, and a
void-floor redraws two UNRELATED mutation-frequency estimates uniformly over [10%, 70%] to measure
how often random pairs would clear the same band by chance.

Reads: nothing. Writes: hema_f8_intron22_frequency_crosscheck.json and void_floor.json under the
cell output directory.
Gate: ratio = antonarakis_pct / lakich_pct must fall in [0.70, 1.43] (+/-30%); both numbers
describe the same lesion in the same disease, hence the tight band.
"""
import json
import os
import random

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "hema_f8_intron22_frequency_crosscheck")
SCRATCH_DIR = OUT_DIR
OUT_PATH = os.path.join(OUT_DIR, "hema_f8_intron22_frequency_crosscheck.json")

ANTONARAKIS_N = 2093
ANTONARAKIS_POS = 880  # 740 type1 + 140 type2 per node text (740+140=880)
LAKICH_PCT = 45.0  # point estimate only, no recoverable CI, PMID8275087

GATE_LO, GATE_HI = 0.70, 1.43


def wilson_ci(k, n, z=1.96):
    phat = k / n
    denom = 1 + z * z / n
    center = phat + z * z / (2 * n)
    half = z * ((phat * (1 - phat) / n + z * z / (4 * n * n)) ** 0.5)
    return ((center - half) / denom, (center + half) / denom)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SCRATCH_DIR, exist_ok=True)
    antonarakis_pct = 100.0 * ANTONARAKIS_POS / ANTONARAKIS_N
    lo, hi = wilson_ci(ANTONARAKIS_POS, ANTONARAKIS_N)
    ratio = antonarakis_pct / LAKICH_PCT
    gate = GATE_LO <= ratio <= GATE_HI

    result = {
        "node": "HEMA-F8-ACTIVITY-DOSE-INTRON22-LESION",
        "antonarakis_1995_pmid7662970": {
            "n": ANTONARAKIS_N, "positive": ANTONARAKIS_POS,
            "pct": round(antonarakis_pct, 3),
            "wilson_95ci_pct": [round(100 * lo, 3), round(100 * hi, 3)],
        },
        "lakich_1993_pmid8275087_point_pct": LAKICH_PCT,
        "ratio_antonarakis_over_lakich": round(ratio, 4),
        "GATE_ratio_in_[0.70,1.43]": gate,
        "verdict": "REPLICATED (magnitude convergence across era/cohort/assay)" if gate else "DIVERGED",
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result, ratio


def void_floor(n_draws=4000, seed=20260728, real_ratio=None):
    """Null: two independent mutation-frequency estimates for a rare recurrent monogenic
    lesion, log-uniform over a wide plausible range [10%,70%] (covers everything from
    uncommon-but-recurrent to majority-cause), redrawn. What fraction of random UNRELATED
    pairs would also land in [0.70,1.43] purely by chance?"""
    rng = random.Random(seed)
    lo_b, hi_b = 10.0, 70.0
    passes = 0
    draws = []
    for _ in range(n_draws):
        a = rng.uniform(lo_b, hi_b)
        b = rng.uniform(lo_b, hi_b)
        draws.append((a, b))
        r = a / b
        if GATE_LO <= r <= GATE_HI:
            passes += 1
    a_vals = [d[0] for d in draws]
    moved = (max(a_vals) - min(a_vals)) > 5.0
    out = {
        "n_draws": n_draws, "sweep_band_pct": [lo_b, hi_b],
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
    void_floor(real_ratio=res["ratio_antonarakis_over_lakich"])
