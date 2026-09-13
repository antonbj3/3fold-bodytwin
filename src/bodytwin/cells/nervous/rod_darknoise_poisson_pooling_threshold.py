"""Rod thermal dark noise, Poisson pooling and the absolute visual threshold.

CLAIM UNDER TEST: "Rod thermal dark-noise rate from single-cell suction-electrode recordings ...
predicts via Poisson pooling (500 rods x 0.1-0.3s) a false-alarm-limited k*=3-9 overlapping the
psychophysically-measured k=4-10." This cell is an independent rebuild from the claim's stated
numbers (dark rate 0.0063-0.01 events/s/rod, N=500 rods, integration window 0.1-0.3s) plus a
pre-registered, disclosed modeling choice for the false-alarm criterion (the one free parameter the
claim text does not pin a single value for -- swept over a standard psychophysical range 1%-10%
rather than hand-picked to fit).

GEOMETRY: pooled dark events over N independent rods integrating for time t is Poisson with
mean mu = N*lambda_dark*t (sum of iid Poisson is Poisson; independence is the modeling
assumption the claim itself makes by "pooling"). An ideal detector sets an integer count
threshold k such that the false-alarm rate P(Poisson(mu) >= k) does not exceed an accepted
psychophysical criterion alpha. k*(mu, alpha) = smallest k with P(X>=k) <= alpha. This is the
SAME machinery as any matched-filter/shot-noise detection threshold (Rose-de Vries / Barlow
1956 "dark light" argument) -- not a curve fit, a direct tail-probability inversion.

PRE-REGISTERED GATES (stated before the run below is executed):
  G1: sweeping mu over the claim's parameter box (lambda_dark in [0.0063,0.01]/s,
      t in [0.1,0.3]s -> mu in [N*0.0063*0.1, N*0.01*0.3] = [0.315, 1.5]) and alpha over a
      standard psychophysical false-alarm band [0.01, 0.10], k* must fall inside [2, 12]
      (a generous envelope around the claim's stated 3-9, wide enough to not be a
      tautology gate but tight enough to be falsifiable).
  G2: the resulting k* range must OVERLAP the independently-measured psychophysical
      Hecht-Shlaer-Pirenne-family k=4-10 (i.e. max(k*_lo,4) <= min(k*_hi,10)).
  G3 (VOID FLOOR, pre-registered to FAIL): lambda_dark=0 (no thermal dark noise) must give
      k*=1 for EVERY alpha (any single absorbed-photon event is already improbable under a
      noiseless floor) -- i.e. the void floor's k* must NOT overlap [4,10], proving the
      dark-noise term (not the pooling/statistics alone) is what produces the match.

Reads: nothing (all parameters embedded).
Writes: rod_darknoise_poisson_pooling_threshold.json (raw numeric results).
"""

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

import os
import json
import numpy as np
from scipy.stats import poisson

OUT_PATH = _os.path.join(OUT_ROOT, "rod_darknoise_poisson_pooling_threshold",
                         "rod_darknoise_poisson_pooling_threshold.json")

N_RODS = 500  # claim's stated pooling geometry
LAMBDA_RANGE = (0.0063, 0.01)  # Baylor-Nunn-Schnapf 1984 primate 0.0063/s; recent human ~0.01/s
T_RANGE = (0.1, 0.3)  # claim's stated integration window, seconds
ALPHA_RANGE = (0.01, 0.10)  # standard psychophysical false-alarm criterion band (not fit)

PSYCHOPHYSICAL_K_BAND = (4, 10)  # HSP1942-family independently measured range, claim's anchor
CLAIM_K_BAND = (3, 9)  # claim's stated result band


def k_star(mu, alpha):
    """Smallest integer k with P(Poisson(mu) >= k) <= alpha. Direct tail inversion, no fitting."""
    k = 0
    while poisson.sf(k - 1, mu) > alpha:  # sf(k-1,mu) = P(X>=k)
        k += 1
        if k > 10000:
            return None
    return k


def sweep(lambda_lo, lambda_hi, t_lo, t_hi, alpha_lo, alpha_hi, n_grid=25):
    lambdas = np.linspace(lambda_lo, lambda_hi, n_grid)
    ts = np.linspace(t_lo, t_hi, n_grid)
    alphas = np.linspace(alpha_lo, alpha_hi, n_grid)
    ks = []
    mus = []
    for lam in lambdas:
        for t in ts:
            mu = N_RODS * lam * t
            mus.append(mu)
            for a in alphas:
                k = k_star(mu, a)
                if k is not None:
                    ks.append(k)
    return ks, mus


def main():
    # --- living model: dark-noise-driven pooled Poisson detector ---
    ks_live, mus_live = sweep(*LAMBDA_RANGE, *T_RANGE, *ALPHA_RANGE)
    k_lo, k_hi = int(min(ks_live)), int(max(ks_live))
    mu_lo, mu_hi = float(min(mus_live)), float(max(mus_live))

    g1_pass = (2 <= k_lo) and (k_hi <= 12)
    g2_overlap_lo = max(k_lo, PSYCHOPHYSICAL_K_BAND[0])
    g2_overlap_hi = min(k_hi, PSYCHOPHYSICAL_K_BAND[1])
    g2_pass = g2_overlap_lo <= g2_overlap_hi

    # claim-band cross-check (informational, not a gate by itself -- G1/G2 are the real gates)
    claim_overlap_lo = max(k_lo, CLAIM_K_BAND[0])
    claim_overlap_hi = min(k_hi, CLAIM_K_BAND[1])
    claim_band_overlaps = claim_overlap_lo <= claim_overlap_hi

    # --- VOID FLOOR: dark rate = 0, same pooling/statistics machinery, same alpha sweep ---
    ks_void = []
    for t in np.linspace(*T_RANGE, 10):
        mu0 = N_RODS * 0.0 * t  # = 0 identically
        for a in np.linspace(*ALPHA_RANGE, 10):
            ks_void.append(k_star(mu0, a))
    k_void_lo, k_void_hi = int(min(ks_void)), int(max(ks_void))
    void_overlap_lo = max(k_void_lo, PSYCHOPHYSICAL_K_BAND[0])
    void_overlap_hi = min(k_void_hi, PSYCHOPHYSICAL_K_BAND[1])
    void_should_fail = not (void_overlap_lo <= void_overlap_hi)  # PASS if it correctly does NOT overlap

    result = {
        "node_id": "vision quantum limit",
        "distinct_from": "none (first computation of this cell)",
        "mu_dark_range_pooled_500rods": [round(mu_lo, 4), round(mu_hi, 4)],
        "k_star_range_live": [k_lo, k_hi],
        "psychophysical_k_band_anchor": list(PSYCHOPHYSICAL_K_BAND),
        "claim_stated_k_band": list(CLAIM_K_BAND),
        "gates": {
            "G1_k_in_generous_envelope_2_12": g1_pass,
            "G2_overlaps_psychophysical_4_10": g2_pass,
            "informational_overlaps_claim_stated_3_9": claim_band_overlaps,
        },
        "void_floor": {
            "description": "lambda_dark forced to 0 (no thermal dark noise), same pooling/alpha machinery",
            "k_star_range_void": [k_void_lo, k_void_hi],
            "void_should_NOT_overlap_psychophysical_band": True,
            "void_correctly_failed_to_overlap": void_should_fail,
        },
        "verdict": "PASS" if (g1_pass and g2_pass and void_should_fail) else "DISAGREE",
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
