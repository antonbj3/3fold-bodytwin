"""Red blood cell osmotic volume -- three closed-form cross-checks from raw geometry and constants.

Re-derives, with no fitting, three stated closed-form relations and compares each with an
independently stated external anchor:
  (1) geometry: maximum isovolumetric-sphere volume from the measured RBC membrane area
      (136 um^2), vs the stated ~150 fL maximum pre-lysis swelling volume.
  (2) van't Hoff osmotic pressure Pi = C*R*T at plasma osmolality 280 mOsm/kg, vs the
      StatPearls/Guyton-family reference 19.32-19.77 mmHg per mOsm.
  (3) Nernst decade slope R*T*ln(10)/F at 37C, vs the textbook 61.5 mV/decade.

Scope: the Boyle-van't-Hoff non-osmotic fraction (b = 0.41-0.48, historically fitted against an
osmotic-fragility %NaCl curve) and the Tosteson-Hoffman pump-leak/Donnan-instability argument are
NOT rebuilt -- mapping %NaCl fragility thresholds onto a normalized volume ratio needs a
reference-tonicity convention the source does not pin numerically, and inventing one would be the
"invented proxy" failure mode. Those legs are informational, not gated.

Reads: nothing. Writes: rbc_osmotic_volume_crosscheck.json under the cell output directory.
Gates (pre-registered, closed-form physics only):
  G1 V_max = A^1.5 / (6*sqrt(pi)) from A = 136 um^2 within 2% of the ~150 fL anchor.
  G2 Pi at C = 280 mOsm/kg, T = 310.15 K, 1 mmHg = 133.322 Pa, inside
     [19.32, 19.77] * 280 = [5409.6, 5535.6] mmHg.
  G3 Nernst decade slope at 310.15 K within 0.5% of 61.5 mV.
  G4 void floor: recomputing G2 with the Boltzmann constant k_B in place of the molar gas constant
     R (no Avogadro correction) must land many orders of magnitude outside the band, proving G2 is
     not an artifact of "any O(10) constant times 280 lands somewhere plausible".
"""
import json
import math
import os

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT = _os.path.join(OUT_ROOT, "rbc_osmotic_volume_crosscheck", "rbc_osmotic_volume_crosscheck.json")

MEMBRANE_AREA_UM2 = 136.0          # RBC membrane surface area
MAX_SPHERE_ANCHOR_FL = 150.0       # independently-stated max pre-lysis volume

R_GAS = 8.314462618                # J/(mol K)
K_B = 1.380649e-23                 # J/K
N_A = 6.02214076e23                # 1/mol
T_37C_K = 310.15
F_FARADAY = 96485.33212            # C/mol
MMHG_PER_PA = 1.0 / 133.322

PLASMA_OSM_MOSM = 280.0            # mOsm/kg ~ mol/m^3 (dilute-solution approx)
STATPEARLS_MMHG_PER_MOSM = (19.32, 19.77)

NERNST_ANCHOR_MV = 61.5


def max_sphere_volume_fl(area_um2):
    """A sphere of given surface area has the maximum volume-to-area ratio of any closed surface
    (isoperimetric inequality): r = sqrt(A/4pi), V=(4/3)pi r^3 = A^1.5/(6*sqrt(pi)).
    Units: um^2 in -> um^3 (=fL) out, since 1 fL = 1 um^3."""
    return (area_um2 ** 1.5) / (6.0 * math.sqrt(math.pi))


def vant_hoff_pressure_mmHg(C_mosm, R, T):
    C_mol_m3 = C_mosm  # 1 mOsm/kg water ~ 1 mmol/L ~ 1 mol/m^3 (dilute approx, as the source uses)
    Pi_pa = C_mol_m3 * R * T
    return Pi_pa * MMHG_PER_PA


def nernst_slope_mV(R, T, F):
    return 1000.0 * R * T * math.log(10.0) / F


def main():
    gates = {}

    v_max = max_sphere_volume_fl(MEMBRANE_AREA_UM2)
    rel_err_sphere = abs(v_max - MAX_SPHERE_ANCHOR_FL) / MAX_SPHERE_ANCHOR_FL
    g1 = rel_err_sphere < 0.02
    gates["G1_max_sphere_volume_within_2pct_of_150fL_anchor"] = bool(g1)

    pi_mmHg = vant_hoff_pressure_mmHg(PLASMA_OSM_MOSM, R_GAS, T_37C_K)
    lo, hi = (STATPEARLS_MMHG_PER_MOSM[0] * PLASMA_OSM_MOSM, STATPEARLS_MMHG_PER_MOSM[1] * PLASMA_OSM_MOSM)
    g2 = lo <= pi_mmHg <= hi
    gates["G2_vant_hoff_pressure_inside_statpearls_band"] = bool(g2)

    nernst = nernst_slope_mV(R_GAS, T_37C_K, F_FARADAY)
    rel_err_nernst = abs(nernst - NERNST_ANCHOR_MV) / NERNST_ANCHOR_MV
    g3 = rel_err_nernst < 0.005
    gates["G3_nernst_slope_within_0.5pct_of_61.5mV"] = bool(g3)

    # G4 void floor: use k_B (per-molecule) in place of R (per-mole) with NO Avogadro correction --
    # a real, tempting unit-confusion error (k_B*N_A = R is the correct relation; using k_B alone
    # is off by ~6e23).
    pi_void_mmHg = vant_hoff_pressure_mmHg(PLASMA_OSM_MOSM, K_B, T_37C_K)
    g4 = not (lo <= pi_void_mmHg <= hi)
    gates["G4_void_floor_wrong_constant_lands_outside_band"] = bool(g4)

    gates = {k: bool(v) for k, v in gates.items()}
    verdict = "CONFIRMED" if all(gates.values()) else "DISAGREE"

    result = {
        "node_id": "MODEL-OSMOTIC-CELL-VOLUME",
        "max_sphere_volume_fL_computed": v_max,
        "max_sphere_volume_fL_anchor": MAX_SPHERE_ANCHOR_FL,
        "rel_err_sphere": rel_err_sphere,
        "vant_hoff_pressure_mmHg_computed": pi_mmHg,
        "statpearls_band_mmHg": [lo, hi],
        "nernst_slope_mV_computed": nernst,
        "nernst_anchor_mV": NERNST_ANCHOR_MV,
        "void_floor_pressure_mmHg_kB_no_avogadro": pi_void_mmHg,
        "gates": gates,
        "verdict": verdict,
        "scope_note": ("b (non-osmotic volume fraction, Boyle-van't-Hoff plot) and the "
                        "Tosteson-Hoffman pump-leak/Donnan leg are NOT rebuilt -- the claim text "
                        "does not pin a numeric reference-tonicity convention needed to map "
                        "%NaCl fragility thresholds onto a normalized volume ratio without "
                        "inventing one; left informational per the false-refutation "
                        "caution, not gated."),
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)

    print(json.dumps({"gates": gates, "verdict": verdict, "max_sphere_volume_fL": v_max,
                       "vant_hoff_pressure_mmHg": pi_mmHg, "nernst_slope_mV": nernst}, indent=2))
    return result


if __name__ == "__main__":
    main()
