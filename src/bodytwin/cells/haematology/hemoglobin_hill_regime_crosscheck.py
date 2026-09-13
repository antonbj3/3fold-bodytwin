"""Hemoglobin O2 cooperativity across a held-out regime -- Severinghaus Hill slope vs Adair anchor.

Evaluates, with NO tuning, the Severinghaus 1979 whole-blood empirical dissociation equation
(PMID 35496, 37C):  S(p) = [ (p^3 + 150p)^-1 * 23400 + 1 ]^-1,  p = pO2 in mmHg,
and computes its Hill coefficient at half-saturation, then compares it with the Yuan et al 2015
Adair-model calibration target n = 3.0 +/- 0.1 (PMID 26244770, purified stripped Hb-A). The two
differ in equation family (empirical rational fit vs mechanistic Adair/MWC), source (clinical
whole blood vs purified solution) and conditions (37C physiological vs 25C stripped), so agreement
tests whether the cooperativity CLASS is conserved across regimes.

Scope: a full 4-site Adair K1-K4 refit is deliberately NOT attempted -- K2, K3 are not individually
identifiable from (P50, n) alone (Jacobian sigma_min ~ 2e-26, numerically degenerate), so only the
independently-derivable Severinghaus side is recomputed here. The MWC two-state rebuild lives in
the hemoglobin_mwc_allostery_rebuild cell and is untouched.

Reads: nothing. Writes: hemoglobin_hill_regime_crosscheck.json under the cell output directory.
Gates (pre-registered):
  G1 P50 solved from the closed form must lie in the physiological band [20, 35] mmHg.
  G2 Hill n at P50 (centered finite difference in ln p, Richardson-extrapolated, step-halved to
     confirm convergence) must reproduce the claimed n_Severinghaus = 2.656 within 1%.
  G3 decisive: |n_Severinghaus - n_Adair_target| < 0.5 (band reflecting the measured literature
     spread of Hill coefficients, 2.3-3.2 across methods).
  G4 void floor: the SAME finite-difference machinery evaluated far from P50 at S = 0.05 (low
     shoulder, pO2 ~ 6.4 mmHg), where a sigmoidal binding curve's local log-log slope collapses
     toward 1, MUST fail the |delta n| < 0.5 test -- so G3 is not a tautology. (S = 0.95 does not
     discriminate; n stays > 2.85 well into the high shoulder.)
"""
import json
import os
from scipy.optimize import brentq

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "hemoglobin_hill_regime_crosscheck", "hemoglobin_hill_regime_crosscheck.json")

N_SEVERINGHAUS_CLAIMED = 2.656
N_ADAIR_TARGET = 3.000
N_ADAIR_SIGMA = 0.100
DELTA_N_THRESHOLD = 0.5
P50_PHYS_BAND = (20.0, 35.0)


def S_of_p(p):
    return 1.0 / (((p ** 3 + 150.0 * p) ** -1) * 23400.0 + 1.0)


def hill_n_at(p, rel_step=1e-4):
    """Hill coefficient n = d ln(Y/(1-Y)) / d ln(p), via centered finite difference in ln(p),
    Richardson-extrapolated over a halved step to confirm convergence (not a single-step guess)."""
    def logit_Y(pp):
        Y = S_of_p(pp)
        return (pp, __import__("math").log(Y / (1.0 - Y)))

    import math

    def n_at_step(h):
        lp = math.log(p)
        p_hi = math.exp(lp + h)
        p_lo = math.exp(lp - h)
        _, ly_hi = logit_Y(p_hi)
        _, ly_lo = logit_Y(p_lo)
        return (ly_hi - ly_lo) / (2.0 * h)

    n1 = n_at_step(rel_step)
    n2 = n_at_step(rel_step / 2.0)
    # Richardson extrapolation (finite-difference error ~O(h^2))
    n_extrap = n2 + (n2 - n1) / 3.0
    converged = abs(n2 - n1) < 1e-4
    return n_extrap, converged


def main():
    gates = {}

    # Solve P50: S(p) = 0.5
    p50 = brentq(lambda p: S_of_p(p) - 0.5, 1.0, 200.0, xtol=1e-8)
    g1 = P50_PHYS_BAND[0] <= p50 <= P50_PHYS_BAND[1]
    gates["G1_p50_in_physiological_band"] = bool(g1)

    n_sev, converged = hill_n_at(p50)
    g2 = abs(n_sev - N_SEVERINGHAUS_CLAIMED) / N_SEVERINGHAUS_CLAIMED < 0.01
    gates["G2_reproduces_claimed_n_severinghaus_1pct"] = bool(g2 and converged)

    delta_n = abs(n_sev - N_ADAIR_TARGET)
    g3 = delta_n < DELTA_N_THRESHOLD
    gates["G3_regime_invariant_cooperativity_class_conserved"] = bool(g3)

    # G4 void floor: evaluate the SAME machinery at S=0.05 (low shoulder, far from P50)
    p_shoulder = brentq(lambda p: S_of_p(p) - 0.05, 0.1, 500.0, xtol=1e-8)
    n_shoulder, _ = hill_n_at(p_shoulder)
    delta_n_void = abs(n_shoulder - N_ADAIR_TARGET)
    g4 = delta_n_void >= DELTA_N_THRESHOLD
    gates["G4_void_floor_offregime_point_fails_band"] = bool(g4)

    z = delta_n / N_ADAIR_SIGMA

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "CONFIRMED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "MODEL-HEMOGLOBIN-O2-COOPERATIVITY",
        "p50_computed_mmHg": p50,
        "n_severinghaus_computed": n_sev,
        "n_severinghaus_claimed": N_SEVERINGHAUS_CLAIMED,
        "n_adair_target": N_ADAIR_TARGET,
        "delta_n": delta_n,
        "z_score": z,
        "n_at_shoulder_S0.05_void_floor": n_shoulder,
        "gates": gates,
        "verdict": verdict,
        "scope_note": ("Adair K1-K4 refit intentionally NOT attempted (the claim's stated gap: "
                        "numerically degenerate, sigma_min~2e-26 from (P50,n) alone); this script "
                        "verifies only the independently-computable Severinghaus side."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict, "n_severinghaus_computed": n_sev,
                       "delta_n": delta_n, "z_score": z}, indent=2))
    return result


if __name__ == "__main__":
    main()
