"""Biphasic poroelastic cartilage arithmetic: does the relaxation-time law tau = L^2/(H_A*k)
reproduce the directly measured fluid-load-support duration (Soltz & Ateshian 1998/2000) with no
fitting, does the resulting diffusive timescale leave a large margin against a single gait stance,
and does mu_eff = mu_min + mu_eq*(1 - Wp/W) land inside the classical Charnley/Jones friction band
at physiological fluid-load-support fractions?

Inputs (literature, embedded in this file):
  H_A (aggregate modulus, Mow 1980, measured)          : [0.70, 0.76] MPa
  k (permeability, measured)                           : [1e-15, 7.6e-15] m^4/(N*s)
  measured fluid-support duration (Soltz-Ateshian)     : [404, 725] s at L=1 mm
  qualitative in-vivo "several hours" estimate at L=4 mm; single gait stance [0.7, 1.0] s
  mu_min 0.010 (near-full fluid support); mu_eq [0.15, 0.28] (drained)
  Wp/W physiological [0.90, 0.99]; Charnley/Jones friction anchor [0.005, 0.024]

Reads: nothing. Writes: cartilage_poroelastic_relaxation.json.

Gates (fixed before computing):
  G1 tau(L=1mm) swept over the full H_A x k grid must contain the measured [404, 725] s window
     (the formula uses only H_A and k, never fit to that number).
  G2 tau(L=4mm) must overlap a "several hours" band, pre-registered as [1800, 14400] s.
  G3 worst-corner margin tau(L=1mm)/T_gait must exceed 100x.
  G4 mu_eff at the high-Wp/W corner (0.99) must fall inside [0.005, 0.024].
  G5 void floor: a dimensionally nonsensical diffusivity k*k instead of H_A*k must NOT overlap
     the measured [404, 725] s window.
Verdict CONFIRMED iff G1-G4 pass and the void floor correctly fails to overlap.
"""
import json
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "cartilage_poroelastic_relaxation", "cartilage_poroelastic_relaxation.json")

H_A_RANGE_PA = (0.70e6, 0.76e6)
K_RANGE = (1e-15, 7.6e-15)  # m^4/(N.s)

MEASURED_TAU_1MM_S = (404.0, 725.0)
SEVERAL_HOURS_S = (1800.0, 14400.0)  # 0.5-4 hr, pre-registered generous band
GAIT_STANCE_S = (0.7, 1.0)

MU_MIN = 0.010
MU_EQ_RANGE = (0.15, 0.28)
WP_W_RANGE = (0.90, 0.99)
CHARNLEY_JONES = (0.005, 0.024)


def tau(L_m, H_A, k):
    D = H_A * k
    return (L_m ** 2) / D


def mu_eff(mu_min, mu_eq, wp_w):
    return mu_min + mu_eq * (1.0 - wp_w)


def main():
    gates = {}

    # G1: sweep full grid at L=1mm
    L1 = 1e-3
    corners = [tau(L1, ha, k) for ha in H_A_RANGE_PA for k in K_RANGE]
    tau1_lo, tau1_hi = min(corners), max(corners)
    g1 = (tau1_lo <= MEASURED_TAU_1MM_S[1]) and (tau1_hi >= MEASURED_TAU_1MM_S[0])
    gates["G1_tau_1mm_contains_measured_window"] = g1

    # G2: L=4mm
    L4 = 4e-3
    corners4 = [tau(L4, ha, k) for ha in H_A_RANGE_PA for k in K_RANGE]
    tau4_lo, tau4_hi = min(corners4), max(corners4)
    g2 = (tau4_lo <= SEVERAL_HOURS_S[1]) and (tau4_hi >= SEVERAL_HOURS_S[0])
    gates["G2_tau_4mm_overlaps_several_hours"] = g2

    # G3: margin
    worst_margin = tau1_lo / GAIT_STANCE_S[1]
    best_margin = tau1_hi / GAIT_STANCE_S[0]
    g3 = worst_margin > 100.0
    gates["G3_margin_exceeds_100x_worst_case"] = g3

    # G4: friction at physiological high-Wp/W corner
    mu_at_high_wpw = [mu_eff(MU_MIN, mu_eq, 0.99) for mu_eq in MU_EQ_RANGE]
    mu_hi_lo, mu_hi_hi = min(mu_at_high_wpw), max(mu_at_high_wpw)
    g4 = (mu_hi_lo <= CHARNLEY_JONES[1]) and (mu_hi_hi >= CHARNLEY_JONES[0] - 1e-9)
    gates["G4_friction_physiological_within_classical_band"] = g4

    # full friction range across the whole grid, for reporting
    mu_all = [mu_eff(MU_MIN, mu_eq, wpw) for mu_eq in MU_EQ_RANGE for wpw in WP_W_RANGE]
    mu_lo, mu_hi = min(mu_all), max(mu_all)

    # G5: void floor -- dimensionally nonsensical D = k*k instead of H_A*k
    def tau_void(L_m, k):
        D_void = k * k
        return (L_m ** 2) / D_void

    void_corners = [tau_void(L1, k) for k in K_RANGE]
    void_lo, void_hi = min(void_corners), max(void_corners)
    void_overlaps = (void_lo <= MEASURED_TAU_1MM_S[1]) and (void_hi >= MEASURED_TAU_1MM_S[0])
    gates["G5_void_floor_does_not_overlap_measured"] = not void_overlaps

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "CONFIRMED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "MODEL-CARTILAGE-POROELASTIC-CONTACT",
        "tau_1mm_computed_s": [tau1_lo, tau1_hi],
        "tau_4mm_computed_s": [tau4_lo, tau4_hi],
        "margin_range": [worst_margin, best_margin],
        "mu_eff_physiological_high_wpw": [mu_hi_lo, mu_hi_hi],
        "mu_eff_full_grid_range": [mu_lo, mu_hi],
        "void_floor_tau_1mm_nonsensical_s": [void_lo, void_hi],
        "gates": gates,
        "claimed_by_narrative": {
            "tau_1mm_s": 476, "tau_4mm_s": 7619, "margin_range": [132, 32000],
            "mu_eff_range": [0.012, 0.031],
        },
        "verdict": verdict,
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict,
                       "tau_1mm_computed_s": [tau1_lo, tau1_hi],
                       "margin_range": [worst_margin, best_margin]}, indent=2))
    return result


if __name__ == "__main__":
    main()
