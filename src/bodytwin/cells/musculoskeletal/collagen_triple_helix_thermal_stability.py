"""Collagen type-I triple-helix melting temperature (Tm): is it over-determined by two non-circular
routes -- (A) a same-sequence causal hydroxylation intervention (Berg & Prockop 1973) and (B) a
composition regression fit ONLY on 3 fish species and then extrapolated out-of-sample to mammalian
imino-acid content -- both converging on the measured mammalian anchor (~37-40 C)?

Inputs (literature, embedded in this file):
  Route A: Tm(0% Y-Hyp)=24 C -> Tm(full Y-Hyp)=39 C, dTm=+15 C  [Berg & Prockop 1973]
  Route B: 3-species same-lab DSC fit (ImAc residues/1000, Td degC):
           cod(157.3,15.2) salmon(161.5,20.6) catfish(193.5,29.3) -> Tm=-37.22+0.345*ImAc, R2=0.92;
           extrapolated to mammalian ImAc 215-225/1000 -> Tm 37.0-40.4 C
  Adversary (2 species, different lab/units): chum salmon(14.8%,12.01 C) tilapia(17.4%,31.31 C)
           -> naive extrapolation overshoots to 61.7-69.2 C
  Anchor: Leikina et al 2002 PNAS 99:1314-1318 (PMID 11805290), mammalian collagen Tm ~37-40 C,
          never used to fit either route.

Reads: nothing. Writes: collagen_triple_helix_thermal_stability.json.

Gates (fixed before computing):
  G1 Route-A causal dTm == 15.0 +/- 0.5 C
  G2 Route-B 3-species fit R2 >= 0.85
  G3 Route-B extrapolated mammalian Tm range overlaps the anchor band [37, 40] C
  G4 monotonic rank-order: ImAc order == Td order in the 3-species dataset
  G5 void floor: the same 3-point fit on a scrambled (ImAc, Td) pairing must fail to overlap
     [37, 40] C for at least one of the 5 non-identity permutations
  G6 adversary (2-species) extrapolation overshoots the anchor by > 20 C, reproducing the
     disclosed non-universality of the composition law
Verdict CONFIRMED iff G1-G4 and G6 pass and G5's void floor shows at least one failure.
"""
import itertools
import json
import os

import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "collagen_triple_helix_thermal_stability", "collagen_triple_helix_thermal_stability.json")

# ---- Route A: causal same-sequence hydroxylation intervention (Berg & Prockop 1973) ----
TM_0_HYP = 24.0
TM_FULL_HYP = 39.0
ROUTE_A_DTM = TM_FULL_HYP - TM_0_HYP

# ---- Route B: 3-species same-lab DSC dataset (well-powered) ----
IMAC_3 = np.array([157.3, 161.5, 193.5])   # cod, salmon, catfish (residues/1000)
TD_3 = np.array([15.2, 20.6, 29.3])        # degC

MAMMAL_IMAC_LOW, MAMMAL_IMAC_HIGH = 215.0, 225.0
ANCHOR_LOW, ANCHOR_HIGH = 37.0, 40.0


def fit_line(x, y):
    """Least-squares slope/intercept + R^2 for a linear fit y = a + b*x."""
    b, a = np.polyfit(x, y, 1)
    yhat = a + b * x
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return a, b, r2


def extrapolate_range(a, b, x_low, x_high):
    vals = [a + b * x_low, a + b * x_high]
    return min(vals), max(vals)


def overlaps(lo1, hi1, lo2, hi2):
    return lo1 <= hi2 and lo2 <= hi1


def main():
    gates = {}

    # G1
    gates["G1_routeA_dTm_15C"] = abs(ROUTE_A_DTM - 15.0) <= 0.5

    # Route B real fit
    a3, b3, r2_3 = fit_line(IMAC_3, TD_3)
    tm_mammal_lo, tm_mammal_hi = extrapolate_range(a3, b3, MAMMAL_IMAC_LOW, MAMMAL_IMAC_HIGH)
    gates["G2_routeB_R2_ge_0.85"] = r2_3 >= 0.85
    gates["G3_routeB_extrapolation_overlaps_anchor"] = overlaps(
        tm_mammal_lo, tm_mammal_hi, ANCHOR_LOW, ANCHOR_HIGH
    )

    # G4 monotonicity
    order_imac = np.argsort(IMAC_3)
    order_td = np.argsort(TD_3)
    gates["G4_monotonic_rank_order"] = bool(np.array_equal(order_imac, order_td))

    # G5 void floor: scramble the (ImAc, Td) pairing (non-identity permutations of Td against ImAc)
    perms = [p for p in itertools.permutations(range(3)) if list(p) != [0, 1, 2]]
    scrambled_results = []
    any_scrambled_fails_to_overlap = False
    for p in perms:
        td_scrambled = TD_3[list(p)]
        a_s, b_s, r2_s = fit_line(IMAC_3, td_scrambled)
        lo_s, hi_s = extrapolate_range(a_s, b_s, MAMMAL_IMAC_LOW, MAMMAL_IMAC_HIGH)
        ov = overlaps(lo_s, hi_s, ANCHOR_LOW, ANCHOR_HIGH)
        scrambled_results.append(
            {"perm": p, "r2": r2_s, "mammal_Tm_range": [lo_s, hi_s], "overlaps_anchor": ov}
        )
        if not ov:
            any_scrambled_fails_to_overlap = True
    gates["G5_void_floor_scrambled_pairing_fails_at_least_once"] = any_scrambled_fails_to_overlap

    # G6 adversary: 2-species cross-lab dataset (different units: % converted to residues/1000 by *10)
    imac_2 = np.array([148.0, 174.0])  # chum salmon 14.8%, tilapia 17.4%, *10 to residues/1000
    td_2 = np.array([12.01, 31.31])
    a2, b2, r2_2 = fit_line(imac_2, td_2)
    tm_adv_lo, tm_adv_hi = extrapolate_range(a2, b2, MAMMAL_IMAC_LOW, MAMMAL_IMAC_HIGH)
    overshoot_low = tm_adv_lo - ANCHOR_HIGH
    gates["G6_adversary_overshoots_by_gt_20C"] = overshoot_low > 20.0

    gates = {k: bool(v) for k, v in gates.items()}
    overall_confirmed = all(gates.values())

    result = {
        "node_id": "MODEL-COLLAGEN-TRIPLE-HELIX-THERMAL-STABILITY",
        "route_A": {"Tm_0Hyp": TM_0_HYP, "Tm_fullHyp": TM_FULL_HYP, "dTm": ROUTE_A_DTM},
        "route_B": {
            "fit_intercept": a3,
            "fit_slope": b3,
            "R2": r2_3,
            "mammal_ImAc_range": [MAMMAL_IMAC_LOW, MAMMAL_IMAC_HIGH],
            "predicted_mammal_Tm_range": [tm_mammal_lo, tm_mammal_hi],
        },
        "anchor_band": [ANCHOR_LOW, ANCHOR_HIGH],
        "void_floor_scrambled_pairings": scrambled_results,
        "adversary_2species": {
            "fit_intercept": a2,
            "fit_slope": b2,
            "R2": r2_2,
            "predicted_mammal_Tm_range": [tm_adv_lo, tm_adv_hi],
            "overshoot_above_anchor_C": overshoot_low,
        },
        "gates": gates,
        "claimed_by_narrative": {
            "route_A_dTm": 15.0,
            "route_B_R2": 0.92,
            "route_B_Tm_range": [37.0, 40.4],
            "adversary_Tm_range": [61.7, 69.2],
        },
        "verdict": "CONFIRMED" if overall_confirmed else "DISAGREE",
    }

    def _default(o):
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.ndarray,)):
            return o.tolist()
        raise TypeError(f"Object of type {o.__class__.__name__} is not JSON serializable")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2, default=_default)

    print(json.dumps({"gates": gates, "verdict": result["verdict"]}, indent=2))
    return result


if __name__ == "__main__":
    main()
