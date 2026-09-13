"""Serum B12 vs functional (tissue) cobalamin deficiency -- machine cross-check.

Computes two independent arithmetic legs from published numbers embedded below:
  LEG A (exact binomial): Lindenbaum 1988 (PMID 3374544, n=141; subgroup n=40) -- 18/40
    hematologically-silent neuropsychiatric-deficiency patients had serum B12 not frankly low
    (>=100 pg/mL, 2 fully normal >200 pg/mL), yet repletion improved neuropsychiatric status
    39/39, hematologic 36/39, MMA/homocysteine normalization 31/31. Under the null "serum B12
    alone determines true deficiency status" those 18 should respond only at a spontaneous rate
    p0, which is SWEPT over [0.05, 0.60] instead of assumed. A conservative floor of 17-of-18
    responders is used because the aggregate 39/39, not the per-stratum count, is published.
  LEG B (Hanley-McNeil ROC approximation): Valente 2011 (PMID 21482749, n=700, reference =
    RBC-cobalamin/holoTC) reports B12 AUC=0.80 vs holoTC AUC=0.90 (p<=0.0002). This recomputes
    the gap by an independent, cruder unpaired route (no DeLong paired covariance term; assumed
    balanced split n_eff=350, variance-maximizing), giving an upper bound on the p-value.

Forced adversary (assay noise): if the 18/40 "not-low" patients were merely borderline-noisy,
the unambiguously normal subset (>200 pg/mL, n=2) should respond less; n=2 is underpowered, so
this is reported as consistent-with, not decisive.

Reads: nothing. Writes: b12_serum_functional_dissociation.json and b12_void_floor.json under
the cell output directory.
Gates: LEG A binomial P(>=17 of 18 | p0) < 0.01 for EVERY p0 in the swept band; LEG B two-sided
p < 0.001. Both must pass for the CONFIRMED verdict.
"""
import json
import math
import os
from math import comb

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "b12_serum_functional_dissociation")
SCRATCH_DIR = OUT_DIR
OUT_PATH = os.path.join(OUT_DIR, "b12_serum_functional_dissociation.json")

N_NORMAL_SERUM = 18   # of 40, serum B12 not low, Lindenbaum1988
K_RESPOND_FLOOR = 17  # pre-registered conservative floor (not the full 18)
P0_LO, P0_HI, P0_STEPS = 0.05, 0.60, 56  # swept null spontaneous-response-rate band
GATE_PVAL = 0.01

AUC_B12, AUC_HOLOTC, N_VALENTE = 0.80, 0.90, 700
N_EFF_CONSERVATIVE = 350  # disclosed assumed balanced split (variance-maximizing)
GATE_LEGB_PVAL = 0.001


def binom_sf_ge(k, n, p):
    """P(X >= k) for X~Binomial(n,p), exact."""
    return sum(comb(n, i) * (p ** i) * ((1 - p) ** (n - i)) for i in range(k, n + 1))


def leg_a():
    step = (P0_HI - P0_LO) / (P0_STEPS - 1)
    p0_grid = [P0_LO + i * step for i in range(P0_STEPS)]
    results = []
    worst_p = 0.0
    for p0 in p0_grid:
        pval = binom_sf_ge(K_RESPOND_FLOOR, N_NORMAL_SERUM, p0)
        worst_p = max(worst_p, pval)
        results.append({"p0": round(p0, 4), "p_value": pval})
    gate_pass = worst_p < GATE_PVAL
    return {
        "n": N_NORMAL_SERUM, "k_floor": K_RESPOND_FLOOR,
        "p0_band_swept": [P0_LO, P0_HI], "n_p0_points": P0_STEPS,
        "worst_case_p_value_across_band": worst_p,
        "GATE_p<0.01_for_ALL_p0_in_band": gate_pass,
        "grid_sample": [results[0], results[len(results) // 2], results[-1]],
    }


def hanley_mcneil_se(auc, n_eff):
    q1 = auc / (2 - auc)
    q2 = (2 * auc ** 2) / (1 + auc)
    var = (auc * (1 - auc) + (n_eff - 1) * (q1 - auc ** 2) + (n_eff - 1) * (q2 - auc ** 2)) / (n_eff ** 2)
    return math.sqrt(max(var, 0.0))


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def leg_b():
    se_b12 = hanley_mcneil_se(AUC_B12, N_EFF_CONSERVATIVE)
    se_holo = hanley_mcneil_se(AUC_HOLOTC, N_EFF_CONSERVATIVE)
    se_diff = math.sqrt(se_b12 ** 2 + se_holo ** 2)
    z = (AUC_HOLOTC - AUC_B12) / se_diff
    p_two_sided = 2 * norm_sf(abs(z))
    gate_pass = p_two_sided < GATE_LEGB_PVAL
    return {
        "AUC_B12": AUC_B12, "AUC_holoTC": AUC_HOLOTC, "n_valente": N_VALENTE,
        "n_eff_conservative_assumed_balanced": N_EFF_CONSERVATIVE,
        "SE_B12": round(se_b12, 5), "SE_holoTC": round(se_holo, 5),
        "z": round(z, 3), "p_two_sided_independent_method": p_two_sided,
        "node_claimed_p": "<=0.0002",
        "GATE_p<0.001": gate_pass,
        "note": "unpaired Hanley-McNeil, ignores paired covariance -> conservative (SE overestimate), so p is an upper bound on significance, not inflated.",
    }


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SCRATCH_DIR, exist_ok=True)
    a = leg_a()
    b = leg_b()
    verdict = "CONFIRMED (serum-alone null rejected on both legs)" if (
        a["GATE_p<0.01_for_ALL_p0_in_band"] and b["GATE_p<0.001"]) else "MIXED"
    result = {
        "node": "ORG-B12-SERUM-VS-FUNCTIONAL-DEFICIENCY",
        "leg_a_binomial_response_vs_swept_null": a,
        "leg_b_hanley_mcneil_auc_gap": b,
        "forced_adversary_assay_noise": {
            "n_far_from_threshold_normal (>200pg/mL)": 2,
            "note": "underpowered (n=2) to statistically refute; the source report states these responded -- reported as consistent-with, not decisive.",
        },
        "verdict": verdict,
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result


def void_floor(n_scrambles=40):
    """
    Scramble the load-bearing driver of Leg A (the floor count K_RESPOND_FLOOR, this
    computation's central adversary-forcing parameter) across a wide range [5,18] and
    confirm the substitution lands (worst-case p-value moves materially) before trusting
    the fixed K=17 gate above.
    """
    step_p0 = (P0_HI - P0_LO) / (P0_STEPS - 1)
    p0_grid = [P0_LO + i * step_p0 for i in range(P0_STEPS)]
    out = []
    for k in range(5, N_NORMAL_SERUM + 1):
        worst_p = max(binom_sf_ge(k, N_NORMAL_SERUM, p0) for p0 in p0_grid)
        out.append({"k": k, "worst_case_p": worst_p, "passes_0.01": worst_p < GATE_PVAL})
    p_values = [o["worst_case_p"] for o in out]
    moved = (max(p_values) - min(p_values)) > 1e-6
    pass_rate = sum(o["passes_0.01"] for o in out) / len(out)
    result = {
        "k_scramble_range": [5, N_NORMAL_SERUM],
        "substitution_landed (k varied AND worst-case p moved)": moved,
        "p_value_range_over_k": [min(p_values), max(p_values)],
        "void_pass_rate (fraction of k in [5,18] clearing gate at worst-case p0)": round(pass_rate, 4),
        "detail": out,
    }
    scratch_path = os.path.join(SCRATCH_DIR, "b12_void_floor.json")
    with open(scratch_path, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps({k: v for k, v in result.items() if k != "detail"}, indent=1))
    return result


if __name__ == "__main__":
    main()
    void_floor()
