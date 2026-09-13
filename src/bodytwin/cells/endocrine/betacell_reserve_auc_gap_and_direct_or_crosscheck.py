"""Beta-cell functional reserve: AUC-gap significance test + independent remission-OR cross-check.

Claim under test: Fan 2026 (PMID41588568, n=2421, Hong Kong/Chinese T2D registry, 10.2y median
follow-up) reports that raw fasting C-peptide alone discriminates future insulin-requirement at
C-index=0.532 (near chance=0.5) while an IR-adjusted clinical risk score reaches C-index=0.729
(n=2421, 1143/2421 events). This cell tests whether that gap is STATISTICALLY significant rather
than just numerically larger, and cross-checks it against a disjoint population.

LEG A (Hanley-McNeil AUC-gap significance): compute the Hanley & McNeil (1982) standard error for
each reported C-index given the paper's n_events=1143 / n_nonevents=1278, then z-test the gap
between the composite score (0.729) and raw C-peptide (0.532). PRE-REGISTERED gate: p<0.001.
Treating the two AUCs as INDEPENDENT (no reported within-subject correlation) is deliberately
CONSERVATIVE: the true DeLong-correlated test on the same patients would give an even LARGER z, so
this design cannot manufacture significance out of noise -- if anything it understates it.

LEG B (independent, disjoint-population magnitude gate): DiRECT/Taylor's remission OR=19.7 (95%CI
7.8-49.8, n=298, UK cluster-RCT, PMID29221645/30078554) is a DIFFERENT outcome (remission, not
insulin-dependence), DIFFERENT population (UK primary-care T2D, not Hong Kong registry), DIFFERENT
design (RCT vs cohort). PRE-REGISTERED gate: OR's 95% CI lower bound must exceed 10 (a
"large-effect" floor fixed before computing, not tuned to the reported 7.8-49.8 CI). This is a
magnitude, not identity, over-determination check.

FORCED ADVERSARY: with n=2421 and this event count, even a modest true AUC gap could look
"significant" by chance under an underpowered or biased SE formula -- Hanley-McNeil is a
closed-form, non-tuned formula (not fit to this claim), and the void floor below scrambles the
event count and the two AUCs to show the gate is not vacuous.

Reads: nothing. Writes: betacell_reserve_auc_gap_and_direct_or_crosscheck.json and void_floor.json.
Gates: leg A p<0.001, leg B CI lower bound > 10.
"""
import json
import os
import random
import math

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

OUT_DIR = os.path.join(OUT_ROOT, "betacell_reserve_auc_gap_and_direct_or_crosscheck")
SCRATCH_DIR = OUT_DIR
OUT_PATH = os.path.join(OUT_DIR, "betacell_reserve_auc_gap_and_direct_or_crosscheck.json")

N_TOTAL = 2421
N_EVENTS = 1143  # progressed to insulin requirement
N_NONEVENTS = N_TOTAL - N_EVENTS

AUC_COMPOSITE = 0.729  # clinical risk score, PMID41588568
AUC_CPEPTIDE = 0.532   # raw fasting C-peptide alone, PMID41588568 (forced-adversary leg in-node)

DIRECT_OR = 19.7
DIRECT_OR_CI = (7.8, 49.8)  # PMID29221645/30078554, n=298

GATE_P = 0.001
GATE_OR_CI_LO_FLOOR = 10.0


def hanley_mcneil_se(auc, n_pos, n_neg):
    q1 = auc / (2.0 - auc)
    q2 = (2.0 * auc * auc) / (1.0 + auc)
    var = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc * auc) + (n_neg - 1) * (q2 - auc * auc)) / (n_pos * n_neg)
    return math.sqrt(max(var, 0.0))


def norm_sf(z):
    # survival function of standard normal via erfc (no scipy dependency)
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    se_composite = hanley_mcneil_se(AUC_COMPOSITE, N_EVENTS, N_NONEVENTS)
    se_cpeptide = hanley_mcneil_se(AUC_CPEPTIDE, N_EVENTS, N_NONEVENTS)
    se_diff = math.sqrt(se_composite ** 2 + se_cpeptide ** 2)  # conservative: assumes independence
    z = (AUC_COMPOSITE - AUC_CPEPTIDE) / se_diff
    p_two_sided = 2 * norm_sf(abs(z))
    gate_a = p_two_sided < GATE_P

    gate_b = DIRECT_OR_CI[0] > GATE_OR_CI_LO_FLOOR

    result = {
        "node": "ENDO-BETA-CELL-FUNCTION",
        "leg_A_hanley_mcneil_auc_gap": {
            "n_events": N_EVENTS, "n_nonevents": N_NONEVENTS,
            "auc_composite": AUC_COMPOSITE, "se_composite": round(se_composite, 5),
            "auc_cpeptide": AUC_CPEPTIDE, "se_cpeptide": round(se_cpeptide, 5),
            "se_diff_conservative_independent": round(se_diff, 5),
            "z": round(z, 4), "p_two_sided": p_two_sided,
            "GATE_p<0.001": gate_a,
        },
        "leg_B_direct_taylor_or_magnitude": {
            "or": DIRECT_OR, "ci95": list(DIRECT_OR_CI),
            "GATE_ci_lower>10": gate_b,
        },
        "verdict": "BOTH LEGS PASS (composite-reserve-score beats raw marker significantly, AND independent disjoint-population large-effect anchor clears its floor)" if (gate_a and gate_b) else "AT LEAST ONE LEG FAILED",
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result, p_two_sided, z


def void_floor(n_draws=4000, seed=20260728, real_p=None):
    """Null A: redraw n_events (hence n_nonevents) uniformly over a wide plausible cohort-size
    range [200, 5000] and redraw the two AUCs independently and log-uniformly around
    [0.5, 0.9] (chance to strong) -- what fraction of random (auc_hi>auc_lo) pairs would ALSO
    clear p<0.001 purely from sample-size/AUC-gap combinations, without any real reserve
    mechanism behind them? Demonstrates the gate is driven by the actual n and actual gap size,
    not automatically satisfied by any two distinct AUCs."""
    rng = random.Random(seed)
    passes = 0
    gaps = []
    for _ in range(n_draws):
        n_tot = rng.randint(200, 5000)
        n_ev = rng.randint(int(0.1 * n_tot), int(0.9 * n_tot))
        n_non = n_tot - n_ev
        auc_lo = rng.uniform(0.50, 0.65)
        auc_hi = rng.uniform(auc_lo, 0.90)
        gaps.append(auc_hi - auc_lo)
        se1 = hanley_mcneil_se(auc_hi, n_ev, n_non)
        se2 = hanley_mcneil_se(auc_lo, n_ev, n_non)
        sed = math.sqrt(se1 ** 2 + se2 ** 2)
        if sed <= 0:
            continue
        zz = (auc_hi - auc_lo) / sed
        pp = 2 * norm_sf(abs(zz))
        if pp < GATE_P:
            passes += 1
    moved = (max(gaps) - min(gaps)) > 0.05
    out = {
        "n_draws": n_draws,
        "n_total_sweep": [200, 5000], "auc_sweep": [0.50, 0.90],
        "substitution_landed": moved,
        "void_pass_rate (random AUC-gap+n combos clearing p<0.001)": round(passes / n_draws, 4),
        "real_p": real_p,
    }
    with open(os.path.join(SCRATCH_DIR, "void_floor.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    res, p, z = main()
    void_floor(real_p=p)
