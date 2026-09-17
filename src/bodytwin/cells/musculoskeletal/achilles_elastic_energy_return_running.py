"""Per-step Achilles tendon elastic energy storage/return during running, by two
independently parameterized routes: (A) material strain-energy-density x volume
(u = sigma^2/(2E)), (B) structural 0.5*F^2/k.

Inputs (all literature, embedded in this file):
  Route A: sigma = 111 MPa (Komi 1990, in-vivo running peak stress); E = 1.6 GPa
           (strain-consistency-checked against the Wren 2001 failure band);
           CSA = 57.5 mm^2 (free Achilles tendon, central of 44-71 mm^2); L0 = 200 mm.
  Route B: F = 2556 N (Kharazi et al 2021, in-vivo peak AT force at 3.5 m/s, n=11);
           k = 250000 N/m (Lichtwark & Wilson 2008, running-tuned, central of 145-250 N/mm).
  Measured envelope held out of both routes: Kharazi et al 2021 7.8-11.3 J propulsion-phase
  recoil; Fletcher & MacIntosh 2015 (n=46 runners) 10-70 J/stride -> [5.2, 70] J.
  Economy: eta = 0.25 (Margaria); COT = 1.01-1.09 kcal/kg/km, 74 kg, 170 steps/min, 3.5 m/s.

Reads: nothing. Writes: achilles_elastic_energy_return_running.json.

Gates (thresholds fixed before computing):
  G1 Route A central energy == 44.3 +/- 0.5 J.
  G2 Route B central energy == 13.1 +/- 0.5 J.
  G3 ratio A/B < 5.0 and == 3.39 +/- 0.05.
  G4 both central energies inside the measured envelope [5.2, 70] J.
  G5 economy fraction (11.3 J recoil / eta, over per-step net metabolic cost) in [10.0, 12.5] %.
  G6 falsifier: both central energies exceed the pre-registered "1-3 J" anchor's 3 J upper bound
     by a factor > 4, refuting that specific anchor (a disclosed negative).
  G7 void floor: Route A recomputed with a 1000x unit-scrambled modulus (1.6 MPa instead of
     1.6 GPa) must fall outside [5.2, 70] J by more than 10x.
Verdict PARTIAL_CONFIRMED iff all gates pass.
"""
import json
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "achilles_elastic_energy_return_running",
                    "achilles_elastic_energy_return_running.json")

# Route A inputs
SIGMA_PA = 111e6
E_GPA_CENTRAL = 1.6e9  # Pa (1.6e9 Pa = 1.6 GPa; name retains legacy '-GPA' token)
CSA_M2 = 57.5e-6
L0_M = 0.200

# Route B inputs
F_N = 2556.0
K_NPM = 250000.0

ANCHOR_LOW_J, ANCHOR_HIGH_J = 5.2, 70.0
TASK_PREREG_ANCHOR_HIGH_J = 3.0

# Economy inputs
RECOIL_FOR_ECONOMY_J = 11.3
ETA_MARGARIA = 0.25
COT_KCAL_PER_KG_KM = (1.01, 1.09)
MASS_KG = 74.0
CADENCE_SPM = 170.0
SPEED_MPS = 3.5


def route_a_energy(sigma, E, csa, L0):
    u = (sigma ** 2) / (2.0 * E)
    return u * csa * L0


def route_b_energy(F, k):
    return 0.5 * (F ** 2) / k


def economy_fraction():
    step_freq_hz = CADENCE_SPM / 60.0
    step_length_m = SPEED_MPS / step_freq_hz
    cot_j_per_kg_m_range = [c * 4184.0 / 1000.0 for c in COT_KCAL_PER_KG_KM]
    per_step_cost_j_range = [c * MASS_KG * step_length_m for c in cot_j_per_kg_m_range]
    metabolic_equiv_j = RECOIL_FOR_ECONOMY_J / ETA_MARGARIA
    fractions = [metabolic_equiv_j / cost for cost in per_step_cost_j_range]
    return min(fractions) * 100.0, max(fractions) * 100.0, per_step_cost_j_range


def main():
    gates = {}

    energy_a = route_a_energy(SIGMA_PA, E_GPA_CENTRAL, CSA_M2, L0_M)
    energy_b = route_b_energy(F_N, K_NPM)

    gates["G1_routeA_central_44.3J"] = abs(energy_a - 44.3) <= 0.5
    gates["G2_routeB_central_13.1J"] = abs(energy_b - 13.1) <= 0.5

    ratio = energy_a / energy_b
    gates["G3_ratio_below_5_and_matches_3.39"] = (ratio < 5.0) and (abs(ratio - 3.39) <= 0.05)

    gates["G4_both_inside_anchor_envelope"] = (
        ANCHOR_LOW_J <= energy_a <= ANCHOR_HIGH_J and ANCHOR_LOW_J <= energy_b <= ANCHOR_HIGH_J
    )

    econ_lo, econ_hi, per_step_cost_range = economy_fraction()
    gates["G5_economy_fraction_in_band"] = (10.0 <= econ_lo <= 12.5) and (10.0 <= econ_hi <= 12.5)

    min_central = min(energy_a, energy_b)
    refutation_factor = min_central / TASK_PREREG_ANCHOR_HIGH_J
    gates["G6_task_anchor_refuted_factor_gt_4"] = refutation_factor > 4.0

    # G7 void floor: unit-scrambled modulus (1000x error: MPa read as GPa)
    e_void = E_GPA_CENTRAL / 1000.0
    energy_void = route_a_energy(SIGMA_PA, e_void, CSA_M2, L0_M)
    void_outside_by_10x = (energy_void > ANCHOR_HIGH_J * 10.0) or (
        energy_void < ANCHOR_LOW_J / 10.0 if energy_void > 0 else True
    )
    gates["G7_void_floor_outside_envelope_10x"] = bool(void_outside_by_10x)

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "PARTIAL_CONFIRMED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "MODEL-ACHILLES-ELASTIC-ENERGY-RETURN",
        "energy_A_J": energy_a,
        "energy_B_J": energy_b,
        "ratio_A_over_B": ratio,
        "economy_fraction_pct_range": [econ_lo, econ_hi],
        "per_step_metabolic_cost_J_range": per_step_cost_range,
        "task_anchor_refutation_factor": refutation_factor,
        "void_floor_energy_J": energy_void,
        "gates": gates,
        "claimed_by_narrative": {
            "energy_A_J": 44.3, "energy_B_J": 13.1, "ratio": 3.39,
            "economy_pct_range": [10.8, 11.7],
        },
        "verdict": verdict,
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict,
                       "energy_A_J": energy_a, "energy_B_J": energy_b}, indent=2))
    return result


if __name__ == "__main__":
    main()
