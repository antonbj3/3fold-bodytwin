"""Two structurally distinct composite-modulus routes for cortical bone -- (A) pure volume-fraction
bounds algebra (Voigt/Reuss/Hashin-Shtrikman, no geometry) and (B) Halpin-Tsai shear-lag mechanics
(adds mineral-platelet aspect ratio and matrix shear transfer) -- tested against the measured
cortical axial modulus (Reilly-Burstein 1975, 17 GPa) at the measured mineral volume fraction
(~43%), using two decorrelated measurement channels (mechanical E and compositional ash fraction).

Reads: nothing. Writes: bone_composite_modulus_bounds.json. Gates G1-G4 plus the void floor decide.

STATED MODEL (verbatim structure from the claim's legs):
  Two-phase composite: mineral (hydroxyapatite, E_m~100GPa) + organic/collagen-water matrix (E_c~2GPa),
  volume fraction Phi (mineral) ~ 0.40-0.45 (central 0.43).
  Route A -- bounds (no geometry):
    Voigt (upper, parallel/iso-strain):  E_V = Phi*E_m + (1-Phi)*E_c
    Reuss (lower, series/iso-stress):    E_R = 1 / (Phi/E_m + (1-Phi)/E_c)
    Hashin-Shtrikman (tighter two-phase bound, using shear moduli G_m,G_c via standard HS formula)
  Route B -- Halpin-Tsai shear-lag (adds mineral-platelet aspect ratio rho=L/t, a distinct mechanism):
    E_HT = E_c * [1 + xi*eta*Phi] / [1 - eta*Phi],   eta = (E_m/E_c - 1) / (E_m/E_c + xi),  xi = 2*rho
    (rho=30, the aspect-ratio convention the claim's decorrelated sibling cell MODEL-BONE-
    MECHANICAL-SAFETY-FACTOR independently landed on with a different modulus formula)

GATE (pre-registered):
  G1: Voigt >= Hashin-Shtrikman-upper >= Reuss (standard bound ORDERING, an algebraic identity of the
      HS construction -- checked, not assumed) at Phi=0.43.
  G2: the measured E=17 GPa (Reilly-Burstein axial, held OUT of both derivations) falls INSIDE the
      [Reuss, Voigt] bracket at Phi=0.43 (necessary condition for the bounds to be non-vacuous).
  G3 (PRE-REGISTERED, REPORTED EVEN IF IT FAILS): Halpin-Tsai at rho=30, Phi=0.43 reproduces E=17GPa
      to within 15% (a genuinely different, geometry-dependent mechanism landing near the SAME
      independently-measured anchor Route A only brackets). ORIENT (forced before accepting a
      negative): tried the standard Halpin-Tsai closed form (xi=2*rho) AND the classic Cox/Fratzl
      staggered shear-lag formula (a structurally different equation, same physical picture) -- BOTH
      over-predict E at phi=0.43,rho=30 (30.9 and 35.9 GPa respectively, vs measured 17 GPa) and
      BOTH require rho~7-11 (not 30) to hit 17GPa exactly at this Em/Ef choice. This is reported as a
      genuine Route-B disagreement, not resolved by silently substituting a smaller rho.
  G4 (over-determination, decorrelated anchor, ALSO REPORTED EVEN IF IT DISAGREES): using the
      corroborating ultrasonic modulus (Rho/Ashman/Turner 1993, 20.7 GPa) INSTEAD of the
      Reilly-Burstein value, Halpin-Tsai's PREDICTED mineral fraction (solving E_HT(Phi)=20.7GPa for
      Phi at rho=30, inverting the SAME formula) is checked against the independently-measured
      ash-fraction range [0.40, 0.45] -- this inherits the SAME rho-mismatch as G3 and is expected,
      pre-registered, to also disagree (predicted Phi ~0.30, below the measured band) as a direct
      consequence, not an independent second failure.

VOID-FLOOR (pre-registered, must FAIL): scramble the phase assignment (swap E_m and E_c, i.e. treat
the SOFT organic matrix as if it were the STIFF mineral phase). This must NOT bracket or reproduce
the measured 17 GPa at Phi=0.43 (predicts a modulus far too low).
"""
import json
import os
import numpy as np
from scipy.optimize import brentq

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
RESULTS_PATH = _os.path.join(OUT_ROOT, "bone_composite_modulus_bounds", "bone_composite_modulus_bounds.json")
OUT = {}

E_MINERAL = 100.0   # GPa, hydroxyapatite (order-of-magnitude literature value)
E_ORGANIC = 2.0      # GPa, collagen-water organic matrix (order-of-magnitude literature value)
NU_MINERAL, NU_ORGANIC = 0.28, 0.35  # Poisson ratios, for shear-modulus conversion (HS bound)
PHI = 0.43
RHO = 30.0  # mineral platelet aspect ratio (length/thickness), dimensionless -- NOT a density/litre; cross-validated convention


def shear_modulus(E, nu):
    return E / (2 * (1 + nu))


def voigt(phi, Em=E_MINERAL, Ec=E_ORGANIC):
    return phi * Em + (1 - phi) * Ec


def reuss(phi, Em=E_MINERAL, Ec=E_ORGANIC):
    return 1.0 / (phi / Em + (1 - phi) / Ec)


def hashin_shtrikman_upper(phi, Em=E_MINERAL, Ec=E_ORGANIC, num=NU_MINERAL, nuc=NU_ORGANIC):
    """HS upper bound for effective bulk+shear -> effective E, standard two-phase formula,
    using the stiffer phase (mineral) as the 'matrix' reference for the upper bound."""
    Gm, Gc = shear_modulus(Em, num), shear_modulus(Ec, nuc)
    Km = Em / (3 * (1 - 2 * num))
    Kc = Ec / (3 * (1 - 2 * nuc))
    G_ref = max(Gm, Gc)  # stiffer phase reference for upper bound
    K_ref = max(Km, Kc)
    f_soft = 1 - phi if Km > Kc else phi
    K_eff = K_ref + f_soft / (1 / (min(Km, Kc) - K_ref) + (1 - f_soft) / (K_ref + 4 / 3 * G_ref))
    G_eff = G_ref + f_soft / (1 / (min(Gm, Gc) - G_ref) + 2 * f_soft * (K_ref + 2 * G_ref) /
                               (5 * G_ref * (K_ref + 4 / 3 * G_ref)))
    E_eff = 9 * K_eff * G_eff / (3 * K_eff + G_eff)
    return E_eff


def halpin_tsai(phi, Em=E_ORGANIC, Ef=E_MINERAL, rho=RHO):
    """Halpin-Tsai shear-lag: Em/Ef roles are matrix/filler -- here matrix=organic, filler=mineral."""
    xi = 2 * rho
    eta = (Ef / Em - 1) / (Ef / Em + xi)
    return Em * (1 + xi * eta * phi) / (1 - eta * phi)


def halpin_tsai_scrambled(phi, Em=E_MINERAL, Ef=E_ORGANIC, rho=RHO):
    """VOID FLOOR: phases swapped -- treats mineral as the compliant 'matrix' and organic as filler."""
    xi = 2 * rho
    eta = (Ef / Em - 1) / (Ef / Em + xi)
    return Em * (1 + xi * eta * phi) / (1 - eta * phi)


# ---------------------------------------------------------------------------
# G1: bound ordering
# ---------------------------------------------------------------------------
E_voigt = voigt(PHI)
E_reuss = reuss(PHI)
E_hs = hashin_shtrikman_upper(PHI)
g1_pass = bool(E_voigt >= E_hs >= E_reuss)
OUT["G1_bound_ordering"] = {"E_voigt": E_voigt, "E_hs_upper": E_hs, "E_reuss": E_reuss, "pass": g1_pass}

# ---------------------------------------------------------------------------
# G2: measured E=17GPa falls inside [Reuss, Voigt]
# ---------------------------------------------------------------------------
E_measured = 17.0
g2_pass = bool(E_reuss <= E_measured <= E_voigt)
OUT["G2_measured_inside_bracket"] = {"E_measured": E_measured, "bracket": [E_reuss, E_voigt], "pass": g2_pass}

# ---------------------------------------------------------------------------
# G3: Halpin-Tsai at rho=30 reproduces 17GPa within 15%
# ---------------------------------------------------------------------------
E_ht = halpin_tsai(PHI)
err_pct = abs(E_ht - E_measured) / E_measured * 100
g3_pass = bool(err_pct < 15.0)
OUT["G3_halpin_tsai_match"] = {"E_halpin_tsai": E_ht, "E_measured": E_measured,
                                "error_pct": err_pct, "pass": g3_pass}

# ---------------------------------------------------------------------------
# G4: over-determination -- invert Halpin-Tsai at E=20.7GPa (ultrasonic anchor), predicted Phi
# should fall in the independently-measured ash-fraction range [0.40, 0.45]
# ---------------------------------------------------------------------------
E_ultrasonic = 20.7


def ht_minus_target(phi, target):
    return halpin_tsai(phi) - target


phi_pred = brentq(ht_minus_target, 0.01, 0.99, args=(E_ultrasonic,))
g4_pass = bool(0.40 <= phi_pred <= 0.45)
OUT["G4_over_determination_ultrasonic"] = {
    "E_ultrasonic_anchor": E_ultrasonic, "phi_predicted_from_inversion": phi_pred,
    "measured_ash_fraction_range": [0.40, 0.45], "pass": g4_pass,
}

# ---------------------------------------------------------------------------
# VOID FLOOR: phase-swapped Halpin-Tsai
# ---------------------------------------------------------------------------
E_ht_void = halpin_tsai_scrambled(PHI)
err_void_pct = abs(E_ht_void - E_measured) / E_measured * 100
void_reproduces = err_void_pct < 15.0
void_floor_pass = not void_reproduces
OUT["VOID_FLOOR"] = {
    "E_halpin_tsai_scrambled": E_ht_void, "error_pct_scrambled": err_void_pct,
    "reproduces_measured_within_15pct": bool(void_reproduces),
    "pass_(void_floor_correctly_fails)": bool(void_floor_pass),
    "note": "mineral and organic phase roles swapped in the Halpin-Tsai formula (compliant matrix "
            "treated as the stiff filler and vice versa)",
}

CORE_GATES = [g1_pass, g2_pass, void_floor_pass]  # Route A (bounds) + void floor
ALL_GATES = [g1_pass, g2_pass, g3_pass, g4_pass, void_floor_pass]
OUT["OVERALL"] = {
    "gates_passed": int(sum(ALL_GATES)), "gates_total": len(ALL_GATES), "all_pass": bool(all(ALL_GATES)),
    "route_A_bounds_core_all_pass": bool(all(CORE_GATES)),
    "route_B_shear_lag_disagrees": bool(not (g3_pass and g4_pass)),
    "verdict": ("DISAGREE (Route B specifically) -- Route A composite bounds (Voigt/Reuss/HS) MATCH "
                "and correctly bracket the measured 17GPa, and the void floor correctly fails; but "
                "Route B (Halpin-Tsai shear-lag) at the claim's cited rho=30 over-predicts the "
                "measured modulus by 82% (30.9 vs 17GPa) and the derived over-determination check "
                "(G4) inherits the same mismatch -- per symmetric-QC this node should NOT be linked "
                "as a full match."),
}
OUT["node_id"] = "MODEL-BONE-COMPOSITE-MODULUS"

os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
with open(RESULTS_PATH, "w") as fh:
    json.dump(OUT, fh, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

print(json.dumps(OUT["OVERALL"], indent=2))
for k in ["G1_bound_ordering", "G2_measured_inside_bracket", "G3_halpin_tsai_match",
          "G4_over_determination_ultrasonic", "VOID_FLOOR"]:
    print(k, OUT[k])
