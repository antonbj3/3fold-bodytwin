"""Endocrine receptor-number reserve: Hill-saturation elasticity and its departure threshold.

Claim under test: a closed-form Hill-saturation elasticity S(g)=n/(1+g^n) (g = normalized
receptor-reserve parameter, proportional to receptor number N) has a reserve-departure
threshold g*=9 (n=1) / g*=3 (n=2) fixed by generic math ALONE (10%-of-max elasticity-departure
criterion), and 4 real, decorrelated physiological/pathological beta-adrenoceptor states
(PMID1376792: rat heart normal, human right atria normal, human left papillary normal, human
end-stage heart-failure) land on this ONE curve with the reported elasticity S_el ranging
0.99-2.91% / 7.41-9.09% / 16.67-20.0% / 47.37-50.0%, monotonically increasing with disease
severity (shrinking receptor reserve), with human right atria described as landing "AT the
edge" (g*=9).

METHOD (re-derived independently from the disclosed formula description): the Black & Leff (1983)
operational model of a saturating agonist gives, at full agonist saturation (A -> infinity),
Emax_frac(g) = g^n / (1 + g^n) where g = tau = k*N/K_E is proportional to receptor number N (K_E =
fixed tissue transducer constant) and is STRUCTURALLY INDEPENDENT of agonist affinity Kd at that
limit (Kd only shifts EC50, not Emax) -- this is the real mechanism behind "Emax blind to
affinity", not an assumption.
  Elasticity: S(g) = d(ln Emax_frac)/d(ln g) = n / (1+g^n)  [derived below, verified vs
  finite-difference]
  Reserve-departure threshold g*: S(g*) = 0.10 * n  =>  g* = 9^(1/n)  (n=1 -> 9, n=2 -> 3) --
  FIXED BEFORE touching the PMID1376792 data (pure algebra, no fit).

CROSS-CHECK (not tautological: the threshold above was fixed from math alone; here the reported
S_el percentages are inverted, assuming n=1 (conservative single-site default, the
Goldbeter-Koshland convention of the insulin_pi3k_akt_signaling cell), to recover the IMPLIED g
per physiological state: g = 1/S_el - 1. PRE-REGISTERED gates (fixed before running): (1) implied
g must be MONOTONICALLY DECREASING from rat-heart-normal -> human-HF (receptor reserve shrinks
with disease severity); (2) human-right-atria-normal's implied g range must BRACKET g*=9 (the "AT
the edge" claim); (3) analytic elasticity must match a numerical finite-difference derivative of
ln(Emax_frac) w.r.t. ln(g) to better than 1e-6 relative error.

Reads: nothing. Writes: endo_receptor_reserve_elasticity_rebuild.json and void_floor.json.
"""
import json
import os
import random

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

OUT_DIR = os.path.join(OUT_ROOT, "endo_receptor_reserve_elasticity_rebuild")
SCRATCH_DIR = OUT_DIR
OUT_PATH = os.path.join(OUT_DIR, "endo_receptor_reserve_elasticity_rebuild.json")

N_HILL = 1  # conservative default (single-site), per insulin_pi3k_akt_signaling.py convention
DEPARTURE_FRAC = 0.10  # 10%-of-max-elasticity departure, fixed before touching data

# PMID1376792, 4 physiological/pathological states, S_el (%) ranges as disclosed in-node
STATES = {
    "rat_heart_normal": (0.99, 2.91),
    "human_right_atria_normal": (7.41, 9.09),
    "human_left_papillary_normal": (16.67, 20.0),
    "human_heart_failure": (47.37, 50.0),
}
STATE_ORDER = ["rat_heart_normal", "human_right_atria_normal",
               "human_left_papillary_normal", "human_heart_failure"]


def emax_frac(g, n=N_HILL):
    return g ** n / (1.0 + g ** n)


def elasticity_analytic(g, n=N_HILL):
    return n / (1.0 + g ** n)


def elasticity_finite_diff(g, n=N_HILL, h=1e-6):
    g_hi = g * (1 + h)
    g_lo = g * (1 - h)
    import math
    dlnE = math.log(emax_frac(g_hi, n)) - math.log(emax_frac(g_lo, n))
    dlng = math.log(g_hi) - math.log(g_lo)
    return dlnE / dlng


def gstar(departure_frac=DEPARTURE_FRAC, n=N_HILL):
    # S(g*) = departure_frac * n  =>  n/(1+g*^n) = departure_frac*n  =>  g*^n = 1/departure_frac - 1
    return (1.0 / departure_frac - 1.0) ** (1.0 / n)


def implied_g_from_Sel_pct(Sel_pct, n=N_HILL):
    S = Sel_pct / 100.0
    if n == 1:
        return 1.0 / S - 1.0
    return (n / S - 1.0) ** (1.0 / n)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    g_star = gstar(DEPARTURE_FRAC, N_HILL)
    g_star_n2 = gstar(DEPARTURE_FRAC, 2)
    gate_gstar = abs(g_star - 9.0) < 1e-9 and abs(g_star_n2 - 3.0) < 1e-9

    # analytic vs finite-difference elasticity cross-check, several g values
    fd_check = []
    max_rel_err = 0.0
    for g in [1.5, 3.0, 9.0, 16.4, 50.0, 200.0]:
        a = elasticity_analytic(g)
        f = elasticity_finite_diff(g)
        rel_err = abs(a - f) / a
        max_rel_err = max(max_rel_err, rel_err)
        fd_check.append({"g": g, "analytic": round(a, 8), "finite_diff": round(f, 8), "rel_err": rel_err})
    gate_fd = max_rel_err < 1e-6

    # invert each state's reported S_el range to implied g range
    implied = {}
    for name, (lo, hi) in STATES.items():
        g_lo = implied_g_from_Sel_pct(hi, N_HILL)  # higher S_el -> lower g
        g_hi = implied_g_from_Sel_pct(lo, N_HILL)
        implied[name] = (round(g_lo, 4), round(g_hi, 4))

    # gate 1: monotonically decreasing implied-g midpoint across STATE_ORDER
    midpoints = [(implied[s][0] + implied[s][1]) / 2.0 for s in STATE_ORDER]
    gate_monotonic = all(midpoints[i] > midpoints[i + 1] for i in range(len(midpoints) - 1))

    # gate 2: human_right_atria_normal implied-g range brackets g*=9
    ra_lo, ra_hi = implied["human_right_atria_normal"]
    gate_bracket = ra_lo <= g_star <= ra_hi
    # honest disclosure: even when the strict bracket gate fails, report HOW close (ratio of
    # nearest implied-g bound to g*) since the node's language ("AT the edge") is a
    # qualitative closeness claim, not a literal straddle -- the gate above is a stricter,
    # pre-registered operationalization of that language, and a near-miss should be shown, not hidden.
    nearest_bound = ra_lo if abs(ra_lo - g_star) < abs(ra_hi - g_star) else ra_hi
    closeness_ratio = nearest_bound / g_star

    all_pass = gate_gstar and gate_fd and gate_monotonic and gate_bracket

    result = {
        "node": "ENDO-RECEPTOR-NUMBER-RESERVE-ELASTICITY-REGISTRY",
        "gstar_n1": round(g_star, 6), "gstar_n2": round(g_star_n2, 6),
        "GATE_gstar_matches_9_and_3": gate_gstar,
        "analytic_vs_finite_diff": fd_check, "max_rel_err": max_rel_err,
        "GATE_fd_match_<1e-6": gate_fd,
        "implied_g_per_state": implied,
        "implied_g_midpoints_in_state_order": [round(m, 4) for m in midpoints],
        "GATE_monotonic_decreasing_reserve": gate_monotonic,
        "human_right_atria_implied_g_range": [ra_lo, ra_hi],
        "GATE_right_atria_brackets_gstar=9": gate_bracket,
        "nearest_bound_to_gstar_closeness_ratio (1.0=exact edge, disclosed even on gate-fail)": round(closeness_ratio, 4),
        "verdict": "REBUILD CONFIRMS node's disclosed numbers (independent re-derivation converges)" if all_pass else "REBUILD DIVERGES from node's disclosed numbers",
    }
    with open(OUT_PATH, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))
    return result, midpoints


def void_floor(n_draws=3000, seed=20260728, real_both_gates=None):
    """Null: randomly PERMUTE which S_el range is assigned to which of the 4 physiological
    states (breaking the real disease-severity ordering) -- what fraction of random
    permutations would ALSO satisfy both (a) monotonic implied-g ordering in STATE_ORDER and
    (b) the right-atria-analog slot bracketing g*=9? Demonstrates the two gates are not
    automatically satisfied by any assignment of these 4 ranges."""
    rng = random.Random(seed)
    ranges = list(STATES.values())
    g_star = gstar(DEPARTURE_FRAC, N_HILL)
    passes = 0
    perms_tried = []
    for _ in range(n_draws):
        perm = ranges[:]
        rng.shuffle(perm)
        perms_tried.append(tuple(perm))
        mids = []
        for (lo, hi) in perm:
            g_lo = implied_g_from_Sel_pct(hi, N_HILL)
            g_hi = implied_g_from_Sel_pct(lo, N_HILL)
            mids.append((g_lo + g_hi) / 2.0)
        mono = all(mids[i] > mids[i + 1] for i in range(len(mids) - 1))
        slot2_lo, slot2_hi = perm[1]
        g2_lo = implied_g_from_Sel_pct(slot2_hi, N_HILL)
        g2_hi = implied_g_from_Sel_pct(slot2_lo, N_HILL)
        bracket = g2_lo <= g_star <= g2_hi
        if mono and bracket:
            passes += 1
    distinct_perms = len(set(perms_tried))
    out = {
        "n_draws": n_draws,
        "distinct_permutations_sampled": distinct_perms,
        "substitution_landed": distinct_perms > 1,
        "void_pass_rate (random permutation also satisfies both gates)": round(passes / n_draws, 4),
        "real_both_gates_passed": bool(real_both_gates),
    }
    with open(os.path.join(SCRATCH_DIR, "void_floor.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    res, mids = main()
    void_floor(real_both_gates=(res["GATE_monotonic_decreasing_reserve"] and res["GATE_right_atria_brackets_gstar=9"]))
