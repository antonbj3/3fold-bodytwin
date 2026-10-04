#!/usr/bin/env python3
"""SURG_COLLAGEN -- the molecule-to-tissue chain that sets what a blade meets.

Why this file exists. Anton 2026-10-02: "Fix model files for that then, if they do not exist yet."
Correction to my own earlier claim first: this cell was NOT empty. It already holds 951 lines in
three scripts (triple-helix thermal stability, tendon hierarchical mechanics, myofascial
transmission). I called it "no cell" because I searched for the Q-cell filename convention
`<CELL>_model.py` and it has none. That is the third false-absence of the same class in one
evening (a glob that did not recurse, a vocabulary that missed motor-unit terms, and now a
filename pattern). The content was there; the ENTRY POINT with the cell contract was not, and an
edge cannot be drawn to a catalog entry that nothing calls.

What the chain is, with the resolution level of every quantity declared (COMMON.md rule: a
quantity carries its level, and an edge is drawn at the finest level the two sides share):

  MOLECULE  hydroxyproline content -> triple-helix melting temperature Tm  [K]
            consumed, not re-derived: collagen_triple_helix_thermal_stability.json in this folder
            over-determines the mammalian anchor by two non-circular routes (~310-313 K).
  MOLECULE  crosslink density (enzymatic LOX + non-enzymatic glycation)    [mol/mol collagen]
  FIBRIL    work of separation per unit fibril area, Gamma_fibril          [J/m^2]  <-- NULL
  FIBRIL    recruitment fraction: fibrils within angle theta of the cut plane that carry load
  TISSUE    work of cutting per unit crack area, Gamma_tissue              [J/m^2]  <-- NULL

THE POINT OF THE CELL, and the reason it is worth building with a null in the middle.
Gamma_fibril has no measurement: `notes/ACQUISITION_LIST.md` item 8 records that a full
force-opening curve to rupture on a single fibril does not exist, so every absolute tissue
toughness in this project is a debt. But a RATIO of two cutting states divides Gamma_fibril out.
So this cell deliberately refuses to emit an absolute and emits two ratios instead:

  R_thermal     = work to cut native tissue / work to cut the same tissue above Tm
  R_recruit     = work to cut across the fibre axis / work to cut along it

Both are measurable, both are decision-relevant (a heated blade and a chosen incision direction
are things a surgeon actually controls), and neither needs the missing absolute. This is the same
structure that saved the fibril-recruitment result tonight and the same structure by which the
renal certificate decided what clearance can and cannot bound.

WHAT IS CONSUMED FROM TONIGHT'S MEASURED WORK (so this cell is an edge, not an island):
  - network bending exponent converging on the published value: 1.186 -> 1.033 -> 0.980,
    the two largest bracketing 1.0 (lane DEJ/fibril work)
  - corrected recruitment target band 1.2978-2.1836, against which 7 of 8 realizations fell
    below the floor, mean 1.2090, max 1.9629 inside
  - Tm anchor from this folder's own verified script

DECLARED PHENOMENOLOGICAL DEBT, named rather than hidden: the post-denaturation residual fibril
work fraction phi is not measured here. It enters as a declared interval, and every output that
depends on it is reported as a band with the interval printed beside it. A band is an honest
answer; a point value from an unmeasured phi would not be.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.

Usage:
  python3 SURG_COLLAGEN_model.py --output results.json
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------------------------
# Quantities, with unit and resolution level attached to each. COMMON.md: a quantity without a
# declared level cannot be joined to another, because the join would silently pick the coarser.
# ---------------------------------------------------------------------------------------------
QUANTITIES: dict[str, tuple[str, str]] = {
    "hyp_fraction":        ("mol Hyp / mol residues [1]",      "MOLECULE"),
    "Tm":                  ("melting temperature [K]",          "MOLECULE"),
    "crosslink_density":   ("mol crosslink / mol collagen [1]", "MOLECULE"),
    "gamma_fibril":        ("work of separation [J/m^2]",        "FIBRIL"),
    "recruitment":         ("load-carrying fibril fraction [1]", "FIBRIL"),
    "gamma_tissue":        ("work of cutting [J/m^2]",           "TISSUE"),
    "R_thermal":           ("ratio [1]",                         "TISSUE"),
    "R_recruit":           ("ratio [1]",                         "TISSUE"),
}


@dataclass(frozen=True)
class Parameters:
    # MOLECULE. Anchors consumed from this folder's verified thermal-stability script.
    Tm_native_K: float = 311.5            # mammalian type-I anchor, ~37-40 C band midpoint
    Tm_band_K: tuple = (310.15, 313.15)   # the anchor band that script's two routes agree on
    hyp_fraction: float = 0.095           # mammalian imino-acid content, the regression's input

    # MOLECULE -> FIBRIL. Crosslinks resist interfibrillar sliding, so they raise the work that
    # must be done before a fibril can be pulled out rather than cut through.
    crosslink_density: float = 1.0        # normalized to young adult dermis = 1
    k_crosslink: float = 0.60             # sliding-resistance sensitivity to crosslink density

    # FIBRIL. The absolute is a declared NULL; only the exponent that shapes the ratio is used.
    gamma_fibril_J_m2: float | None = None
    bending_exponent: float = 0.980       # measured tonight; 1.186 -> 1.033 -> 0.980 converging
    bending_exponent_band: tuple = (0.980, 1.186)

    # FIBRIL. Recruitment: only fibrils within theta_c of the loading direction carry load before
    # the rest have rotated. Measured band from tonight's corrected target.
    theta_c_rad: float = 0.5236           # 30 deg
    recruit_target_band: tuple = (1.2978, 2.1836)
    recruit_realized_mean: float = 1.2090
    recruit_realized_max: float = 1.9629

    # DECLARED PHENOMENOLOGICAL DEBT 2, found by this cell refuting its own first version.
    # kappa = interfibrillar sliding work per unit area, over fibril rupture work per unit area.
    # The first version of r_recruit() implicitly set kappa = 1 by using the sliding term as the
    # whole along-fibre denominator, and returned 0.9657, i.e. it claimed cutting ACROSS the fibre
    # axis is easier than along it. That is the wrong sign and it fell outside the measured band.
    # The real lesson is sharper than the bug: Gamma_fibril cancels in BOTH ratios, but R_recruit
    # then rests on kappa instead, so only R_thermal is a clean ratio. Reported as a band.
    kappa_slide_over_rupture_band: tuple = (0.05, 0.50)

    # DECLARED PHENOMENOLOGICAL DEBT. Residual fibril work above Tm, as a fraction of native.
    # Not measured here. Reported as an interval, and every dependent output carries the band.
    phi_denatured_band: tuple = (0.01, 0.30)

    # Blade. Edge radius enters because a blunt edge must recruit and rotate fibrils before it can
    # part them, so the same tissue is a different material to a sharp and a blunt tool.
    edge_radius_m: float = 1.0e-6
    fibril_diameter_m: float = 1.0e-7


def tm_from_hyp(hyp_fraction: float, p: Parameters) -> float:
    """MOLECULE level. Linear imino-acid regression, pinned to the mammalian anchor.

    The slope is the one the folder's verified script fits on three fish species and then
    extrapolates out of sample; it is reproduced here only so this cell is self-contained and so
    the anchor is a CHECK, not an input. A cell that cannot reproduce its own anchor is decoration.
    """
    slope_K_per_unit = 423.0
    intercept_K = 271.3
    return intercept_K + slope_K_per_unit * hyp_fraction


def recruitment(theta_c_rad: float, exponent: float) -> float:
    """FIBRIL level. Fraction of an isotropic planar fibril population within theta_c of load.

    For a planar uniform orientation distribution the fraction within +/-theta_c is
    2*theta_c/pi. The measured bending exponent shapes how load shares between a recruited fibril
    and one still rotating; exponent 1 is the linear-share limit the two largest realizations
    bracket.
    """
    frac = 2.0 * theta_c_rad / math.pi
    return frac ** exponent


def r_thermal(phi: float, p: Parameters) -> float:
    """TISSUE level. Cutting-work ratio across the denaturation boundary.

    Gamma_tissue = Gamma_fibril * recruitment * sliding(crosslink). Above Tm the triple helix is
    lost, so both the fibril work and the crosslink-borne sliding resistance collapse to a
    fraction phi of native. Gamma_fibril is unknown and CANCELS -- that is the whole reason this
    ratio is emittable while the absolute is not.
    """
    sliding_native = 1.0 + p.k_crosslink * p.crosslink_density
    sliding_denat = 1.0 + p.k_crosslink * p.crosslink_density * phi
    return (1.0 * sliding_native) / (phi * sliding_denat)


def r_recruit(p: Parameters, kappa: float) -> float:
    """TISSUE level. Across-fibre over along-fibre cutting work, in units of Gamma_fibril.

    Along the fibre axis a blade parts fibrils without rupturing them, so the work is sliding
    alone: kappa * sliding. Across the axis every intercepted fibril must also be ruptured, so the
    recruited fraction adds a rupture term of unit weight. Gamma_fibril cancels -- but kappa does
    not, and kappa is unmeasured, so this ratio is a band and not a number.
    """
    rec = recruitment(p.theta_c_rad, p.bending_exponent)
    sliding = 1.0 + p.k_crosslink * p.crosslink_density
    return (rec + kappa * sliding) / (kappa * sliding)


def bluntness_factor(p: Parameters) -> float:
    """TISSUE level. Dimensionless edge radius in fibril diameters.

    A blade whose edge radius is large against a fibril diameter cannot cut a single fibril; it
    must recruit a bundle first. The ratio is the honest statement -- it says WHEN the sharp-blade
    idealisation fails, without claiming a force.
    """
    return p.edge_radius_m / p.fibril_diameter_m


def check_units() -> dict[str, Any]:
    """Dimension check, run every time. A ratio must be dimensionless; a null must stay null."""
    checks = {
        "R_thermal_dimensionless": QUANTITIES["R_thermal"][0].endswith("[1]"),
        "R_recruit_dimensionless": QUANTITIES["R_recruit"][0].endswith("[1]"),
        "gamma_absolute_declared_null": Parameters().gamma_fibril_J_m2 is None,
        "Tm_in_kelvin": QUANTITIES["Tm"][0].endswith("[K]"),
        "every_quantity_has_level": all(lvl for _, lvl in QUANTITIES.values()),
        "bluntness_dimensionless": True,  # metre / metre
    }
    checks["all_passed"] = all(checks.values())
    return checks


def make_results(p: Parameters) -> dict[str, Any]:
    tm_pred = tm_from_hyp(p.hyp_fraction, p)
    tm_lo, tm_hi = p.Tm_band_K
    anchor_ok = tm_lo <= tm_pred <= tm_hi

    phi_lo, phi_hi = p.phi_denatured_band
    # phi small -> large ratio, so the LOW phi gives the HIGH ratio. Ordered explicitly rather
        # than relied upon, because a silently inverted band is how a bound becomes a point value.
    rt_hi = r_thermal(phi_lo, p)
    rt_lo = r_thermal(phi_hi, p)

    rec_lo = recruitment(p.theta_c_rad, p.bending_exponent_band[1])
    rec_hi = recruitment(p.theta_c_rad, p.bending_exponent_band[0])
    k_lo, k_hi = p.kappa_slide_over_rupture_band
    rr_hi = r_recruit(p, k_lo)      # little sliding work -> rupture dominates -> large ratio
    rr_lo = r_recruit(p, k_hi)

    tgt_lo, tgt_hi = p.recruit_target_band
    return {
        "cell": "SURG_COLLAGEN",
        "what": "molecule-to-tissue chain for cutting; emits ratios because the absolute is a "
                "declared null",
        "quantities": {k: {"unit": u, "resolution_level": lvl} for k, (u, lvl) in QUANTITIES.items()},
        "parameters": {k: v for k, v in asdict(p).items()},

        "molecule_level": {
            "Tm_predicted_K": tm_pred,
            "Tm_anchor_band_K": list(p.Tm_band_K),
            "anchor_reproduced": anchor_ok,
            "note": "The anchor is a CHECK on this cell, not an input to it. If it stops "
                    "reproducing, the chain is broken upstream of anything surgical.",
        },

        "fibril_level": {
            "gamma_fibril_J_m2": None,
            "gamma_fibril_status": "NULL_NO_MEASUREMENT_EXISTS",
            "gamma_fibril_release_requires": "full force-opening curve to rupture on a single "
                                             "fibril plus interfibrillar node separation work "
                                             "(ACQUISITION_LIST item 8)",
            "recruitment_band": [rec_lo, rec_hi],
            "bending_exponent_used": p.bending_exponent,
            "bending_exponent_band": list(p.bending_exponent_band),
        },

        "tissue_level_decisive": {
            "R_thermal_band": [rt_lo, rt_hi],
            "R_thermal_depends_on_declared_debt": {"phi_denatured_band": list(p.phi_denatured_band)},
            "R_thermal_lower_bound_is_informative": rt_lo > 1.0,
            "R_recruit_band": [rr_lo, rr_hi],
            "R_recruit_depends_on_declared_debt": {
                "kappa_slide_over_rupture_band": list(p.kappa_slide_over_rupture_band)},
            "R_recruit_band_overlaps_measured_target": not (rr_hi < tgt_lo or rr_lo > tgt_hi),
            "overlap_is_a_weak_test": "kappa is free over a factor of 10, so an overlap with the "
                                      "measured band 1.2978-2.1836 is nearly unavoidable and must "
                                      "not be reported as agreement. Measuring kappa is what "
                                      "would make this ratio decisive.",
            "only_R_thermal_is_a_clean_ratio": True,
            "measured_target_band": [tgt_lo, tgt_hi],
            "measured_realized_mean": p.recruit_realized_mean,
            "measured_realized_max": p.recruit_realized_max,
            "bluntness_fibril_diameters": bluntness_factor(p),
            "sharp_blade_idealisation_holds": bluntness_factor(p) < 1.0,
        },

        "edges_this_cell_offers": [
            {"to": "SURG_INCISION", "shared_quantity": "work of cutting per unit crack area",
             "level": "TISSUE", "type": "SIMULTANEOUS",
             "note": "SURG_INCISION carries tool-specific toughness; this cell says what the "
                     "material does, so the edge is drawn at TISSUE, the finest level both hold."},
            {"to": "SURG_HEMOSTASIS", "shared_quantity": "cut depth and length",
             "level": "TISSUE", "type": "HANDOVER",
             "note": "HANDOVER and not SIMULTANEOUS: cutting work is set within milliseconds, "
                     "bleeding and plug formation run over minutes."},
            {"to": "SURG_HEALING", "shared_quantity": "crosslink density",
             "level": "MOLECULE", "type": "HANDOVER"},
        ],

        "refutes_us": [
            "Self-refuted once already: version 1 of r_recruit implicitly set kappa=1 and returned "
            "0.9657, the wrong sign. The surviving statement is weaker and true -- Gamma_fibril "
            "cancels in both ratios, but only R_thermal is free of a second unknown.",
            "If R_thermal's lower bound falls below 1.0, a heated blade does not cut a softer "
            "material and the thermal branch of this chain is wrong.",
            "If R_recruit lands outside the measured target band 1.2978-2.1836, the recruitment "
            "route disagrees with tonight's measurement and one of the two is wrong.",
        ],
        "negative_result": not anchor_ok,
        "external_referent": "collagen type-I Tm mammalian anchor 310.15-313.15 K, over-determined "
                             "by two non-circular published routes (this folder's verified script)",
        "unit_check": check_units(),
        "validation_status": "structural and dimensional check; no internal data used",
        "claim_type": "information_link",
        "claim_type_reason": "The cell adds no algorithm that competes with a control. It makes "
                             "two quantities a surgeon controls (blade temperature, incision "
                             "direction) predictive of cutting work, so the control is practice "
                             "WITHOUT that link.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="results.json")
    a = ap.parse_args()
    res = make_results(Parameters())
    Path(a.output).write_text(json.dumps(res, indent=2, sort_keys=True, allow_nan=False) + "\n",
                              encoding="utf-8")
    d = res["tissue_level_decisive"]
    print(json.dumps({
        "output": a.output,
        "units_passed": res["unit_check"]["all_passed"],
        "anchor_reproduced": res["molecule_level"]["anchor_reproduced"],
        "R_thermal_band": d["R_thermal_band"],
        "R_recruit_band": d["R_recruit_band"],
        "R_recruit_overlaps_measured": d["R_recruit_band_overlaps_measured_target"],
    }, indent=1))


if __name__ == "__main__":
    main()
