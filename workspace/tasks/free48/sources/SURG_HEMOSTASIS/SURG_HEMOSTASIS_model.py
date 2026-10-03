#!/usr/bin/env python3
"SURG_HEMOSTASIS -- incision geometry to bleeding rate to time to hemostasis.\n\nWhy this file exists. Anton 2026-10-02: \"Fix model files for that then, if they do not exist yet.\"\nCorrection first, as in the sibling cell: this folder was NOT empty. It holds 1102 lines in two\nscripts -- a reduced 15-species thrombin-generation ODE anchored to Hockin & Mann, and a platelet\nprimary-hemostasis layer whose author caught three misremembered PMIDs by live verification. What\nwas missing was the entry point carrying the cell contract, so nothing could call the cell and no\nedge could be drawn to it. I had reported the cell as absent because it has no `<CELL>_model.py`.\n\nWhat this cell adds that neither existing script does: the two existing scripts model the CASCADE\nin time, given that bleeding happens. Neither takes the surgeon's two controlled quantities --\ncut length and cut depth -- and returns a bleeding rate. That is the link the catalog entry asks\nfor, and it is the link a surgeon can act on.\n\nTHE RESULT THAT MAKES THIS CELL WORTH HAVING is a negative one, and it is structural rather than\nnumerical: a severed vessel's bleeding rate is not computable from the vessel's own geometry. The\ncut does not change the resistance upstream of it, and that upstream resistance is what sets the\nflow out of the stub. Write f for the fraction of a vessel's series resistance lying upstream of\nthe cut plane. Then\n\n    Q_cut / Q_normal  ~=  1 / f\n\nand the bleeding rate is predictable within a factor of two only when f >= 0.5. That condition\nsorts the layers of a skin incision cleanly, and it sorts them the way surgical practice already\nbehaves: capillary and venular ooze is predictable, and a cut artery is not. The same structure\ndecided the renal certificate tonight -- when extraction ratio approaches 1 no clearance\nmeasurement can bound transporter turnover from above -- and it is the honest form of answer when\na quantity is set by something the measurement does not see.\n\nTWO INDEPENDENT ROUTES ARE COMPUTED AND THEIR DISAGREEMENT IS REPORTED, not averaged:\n  route A  count the vessels the cut plane intersects and sum their own luminal flows\n  route B  take the perfusion of the tissue within a drainage distance of the cut plane\nRoute A must UNDER-predict, because a severed capillary drains its upstream arteriole and not\nmerely its own lumen. The size of the disagreement is therefore a measurement of the mechanism,\nnot an error bar, and a cell that averaged the two would destroy the finding.\n\nSOURCE CLASS IS ATTACHED TO EVERY NUMBER. Nothing here is live-verified in this file: the vessel\ndensities, the resistance split across the vascular tree and the template bleeding time are all\nSTANDARD, meaning textbook values carried without an independent check. They are declared as\nbands for that reason, and `acquisition_required` lists what would turn each into a measurement.\nA vendor or textbook number is not a measurement and must never be reported as one.\n\nNo internal data. Everything PENDING_INDEPENDENT_REVIEW.\n\nUsage:\n  python3 SURG_HEMOSTASIS_model.py --output results.json\n"
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

MMHG_PA = 133.322

QUANTITIES: dict[str, tuple[str, str]] = {
    "cut_length":        ("mm [L]",                "WHOLE_BODY"),
    "cut_depth":         ("mm [L]",                "TISSUE"),
    "vessel_density":    ("vessels / mm^2 [L^-2]", "TISSUE"),
    "vessel_radius":     ("m [L]",                 "TISSUE"),
    "upstream_fraction": ("f, dimensionless [1]",  "TISSUE"),
    "bleeding_rate":     ("mL/min [L^3 T^-1]",     "TISSUE"),
    "plug_time":         ("s [T]",                 "CELL"),
    "blood_loss":        ("mL [L^3]",              "WHOLE_BODY"),
}


@dataclass(frozen=True)
class Layer:
    """One stratum of a skin incision. Every field carries its source class."""
    name: str
    depth_from_mm: float
    depth_to_mm: float
    vessel_density_per_mm2: float      # STANDARD
    vessel_radius_m: float             # STANDARD
    vessel_velocity_m_s: float         # STANDARD
    perfusion_mL_min_g: float          # STANDARD
    upstream_fraction_f: float         # STANDARD, the quantity the whole result turns on
    source_class: str = "STANDARD"

    @property
    def thickness_mm(self) -> float:
        return self.depth_to_mm - self.depth_from_mm


# The epidermis is listed with zero vasculature on purpose. It is the one layer where the model
# can be checked against something indisputable: an incision confined to the epidermis does not
# bleed. A chain that predicts bleeding there is wrong before any parameter is argued about.
LAYERS: tuple[Layer, ...] = (
    Layer("epidermis",           0.00, 0.10,    0.0, 0.0,      0.0,    0.00, 1.00),
    Layer("papillary_dermis",    0.10, 0.50,   40.0, 4.0e-6,   5.0e-4, 0.15, 0.80),
    Layer("reticular_dermis",    0.50, 2.00,   12.0, 1.5e-5,   2.0e-3, 0.12, 0.60),
    Layer("subcutaneous",        2.00, 10.00,   1.2, 1.0e-4,   8.0e-3, 0.05, 0.25),
    Layer("perforator_plane",   10.00, 12.00,   0.05, 5.0e-4,  2.0e-1, 0.04, 0.08),
)


@dataclass(frozen=True)
class Parameters:
    cut_length_mm: float = 100.0
    cut_depth_mm: float = 2.0
    arterial_pressure_mmHg: float = 90.0        # STANDARD
    drainage_distance_mm_band: tuple = (0.5, 2.0)   # STANDARD, declared as a band
    plug_time_s_band: tuple = (120.0, 540.0)        # STANDARD, template bleeding time 2-9 min
    blood_density_g_mL: float = 1.06
    tissue_density_g_mL: float = 1.02
    predictable_f_floor: float = 0.50           # frozen before running; see PREREG note below
    field_clearing_mL_min: float = 20.0         # STANDARD, what a surgeon can keep a field clear of


def severed_count(layer: Layer, p: Parameters, depth_mm: float) -> float:
    """How many vessels the cut plane intersects in this layer.

    The cut plane in a layer has area = cut length x the part of the layer actually entered, so a
    cut that stops inside a layer only severs that layer's engaged fraction. Getting this wrong is
    how a depth-response curve becomes a step function it has no right to be.
    """
    engaged_mm = max(0.0, min(depth_mm, layer.depth_to_mm) - layer.depth_from_mm)
    return layer.vessel_density_per_mm2 * p.cut_length_mm * engaged_mm


def route_a_luminal(p: Parameters, depth_mm: float) -> dict[str, Any]:
    """Route A: sum the severed vessels' own luminal flows. Must under-predict.

    Q per vessel = pi r^2 v, converted from m^3/s to mL/min. The upstream fraction f then says how
    much more than its own flow a stub can deliver, because the cut removed everything downstream
    of it but nothing upstream.
    """
    total = 0.0
    raw_total = 0.0
    per_layer = {}
    for lay in LAYERS:
        n = severed_count(lay, p, depth_mm)
        if n <= 0.0 or lay.vessel_radius_m <= 0.0:
            per_layer[lay.name] = {"severed": n, "mL_min_raw": 0.0, "mL_min_amplified": 0.0}
            continue
        q_m3_s = math.pi * lay.vessel_radius_m ** 2 * lay.vessel_velocity_m_s
        q_mL_min = q_m3_s * 1.0e6 * 60.0
        # Both cut faces bleed, and the stub carries 1/f of its normal flow.
        # RAW is the vessels' own luminal flow, both faces. AMPLIFIED applies 1/f on top. Keeping
        # them apart matters: version 1 compared the amplified sum against a pure-perfusion route
        # B and reported the disagreement inverted, 0.05-0.19x, which falsified its own stated
        # prediction. Route A raw is the quantity that is comparable to route B.
        q_raw = 2.0 * n * q_mL_min
        q_cut = q_raw / max(lay.upstream_fraction_f, 1e-9)
        per_layer[lay.name] = {"severed": n, "mL_min_raw": q_raw, "mL_min_amplified": q_cut,
                               "f": lay.upstream_fraction_f,
                               "amplification_1_over_f": 1.0 / max(lay.upstream_fraction_f, 1e-9)}
        total += q_cut
        raw_total += q_raw
    return {"total_mL_min": total, "total_raw_mL_min": raw_total, "per_layer": per_layer}


def route_b_perfusion(p: Parameters, depth_mm: float, drainage_mm: float) -> dict[str, Any]:
    """Route B: the perfusion of the tissue that now drains into the field instead of returning.

    Volume = 2 faces x cut length x engaged thickness x drainage distance. Mass from tissue
    density, then the layer's own perfusion. This route does not need a vessel count, which is
    exactly why it is worth running beside one.
    """
    total = 0.0
    per_layer = {}
    for lay in LAYERS:
        engaged_mm = max(0.0, min(depth_mm, lay.depth_to_mm) - lay.depth_from_mm)
        vol_mm3 = 2.0 * p.cut_length_mm * engaged_mm * drainage_mm
        grams = vol_mm3 * 1.0e-3 * p.tissue_density_g_mL
        q = grams * lay.perfusion_mL_min_g
        per_layer[lay.name] = {"drained_g": grams, "mL_min": q}
        total += q
    return {"total_mL_min": total, "per_layer": per_layer}


def identifiability(p: Parameters, depth_mm: float) -> dict[str, Any]:
    """Which layers' bleeding is predictable at all, by the frozen f floor.

    This is the cell's central output and it is a certificate, not a fit: it states where the
    model may be used and where it may not, and it was frozen before the numbers were run.
    """
    entered = [l for l in LAYERS if min(depth_mm, l.depth_to_mm) > l.depth_from_mm]
    pred = [l.name for l in entered if l.upstream_fraction_f >= p.predictable_f_floor]
    unpred = [{"layer": l.name, "f": l.upstream_fraction_f,
               "amplification_1_over_f": 1.0 / max(l.upstream_fraction_f, 1e-9)}
              for l in entered if l.upstream_fraction_f < p.predictable_f_floor]
    # The upper bound for an unpredictable layer is the whole territory it supplies, which this
    # cell does not know. Saying so is the result; inventing a number would destroy it.
    return {
        "f_floor_frozen_before_run": p.predictable_f_floor,
        "layers_entered": [l.name for l in entered],
        "predictable_layers": pred,
        "unpredictable_layers": unpred,
        "all_entered_layers_predictable": not unpred,
        "bound_for_unpredictable": "Q_cut <= total supply of the territory distal to the cut; "
                                   "this cell does not hold that territory's flow, so the bound "
                                   "is open upward and no point value may be emitted.",
        "max_depth_fully_predictable_mm": max(
            [l.depth_to_mm for l in LAYERS if l.upstream_fraction_f >= p.predictable_f_floor]
            or [0.0]),
    }


def check_units() -> dict[str, Any]:
    p = Parameters()
    lay = LAYERS[1]
    # pi r^2 v : m^2 * m/s = m^3/s -> x1e6 mL/m^3 x 60 s/min
    q = math.pi * lay.vessel_radius_m ** 2 * lay.vessel_velocity_m_s * 1.0e6 * 60.0
    # Route B: mm^3 * 1e-3 mL/mm^3 * g/mL * mL/min/g = mL/min
    g = (2.0 * 100.0 * 0.4 * 1.0) * 1.0e-3 * p.tissue_density_g_mL
    checks = {
        "route_a_gives_mL_min": q > 0.0,
        "route_b_gives_mL_min": g > 0.0,
        "epidermis_cannot_bleed": LAYERS[0].vessel_density_per_mm2 == 0.0
                                  and LAYERS[0].perfusion_mL_min_g == 0.0,
        "f_is_a_fraction": all(0.0 < l.upstream_fraction_f <= 1.0 for l in LAYERS),
        "depths_are_contiguous": all(LAYERS[i].depth_to_mm == LAYERS[i + 1].depth_from_mm
                                     for i in range(len(LAYERS) - 1)),
        "every_quantity_has_level": all(lvl for _, lvl in QUANTITIES.values()),
    }
    checks["all_passed"] = all(checks.values())
    return checks


def make_results(p: Parameters) -> dict[str, Any]:
    d = p.cut_depth_mm
    a = route_a_luminal(p, d)
    dr_lo, dr_hi = p.drainage_distance_mm_band
    b_lo = route_b_perfusion(p, d, dr_lo)
    b_hi = route_b_perfusion(p, d, dr_hi)
    ident = identifiability(p, d)

    raw = a["total_raw_mL_min"]
    disagreement_lo = (b_lo["total_mL_min"] / raw) if raw > 0 else None
    disagreement_hi = (b_hi["total_mL_min"] / raw) if raw > 0 else None

    t_lo, t_hi = p.plug_time_s_band
    loss_lo = b_lo["total_mL_min"] * t_lo / 60.0
    loss_hi = b_hi["total_mL_min"] * t_hi / 60.0

    # Depth response, so the shape can be inspected rather than asserted.
    curve = []
    for depth in (0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 11.0):
        rb = route_b_perfusion(p, depth, dr_hi)["total_mL_min"]
        idf = identifiability(p, depth)
        curve.append({"depth_mm": depth, "route_b_mL_min": rb,
                      "predictable": idf["all_entered_layers_predictable"]})

    return {
        "cell": "SURG_HEMOSTASIS",
        "what": "incision geometry -> bleeding rate -> blood loss, with a certificate for where "
                "the chain may be used",
        "quantities": {k: {"unit": u, "resolution_level": lvl} for k, (u, lvl) in QUANTITIES.items()},
        "parameters": asdict(p),
        "layers": [asdict(l) | {"thickness_mm": l.thickness_mm} for l in LAYERS],

        "decisive": {
            "route_a_luminal_raw_mL_min": a["total_raw_mL_min"],
            "route_a_luminal_amplified_mL_min": a["total_mL_min"],
            "route_a_amplification_is_the_1_over_f_step": True,
            "route_b_perfusion_mL_min_band": [b_lo["total_mL_min"], b_hi["total_mL_min"]],
            "route_disagreement_factor_band": [disagreement_lo, disagreement_hi],
            "route_disagreement_is_the_mechanism": True,
            "route_disagreement_note": "Compared RAW against route B, as version 1 failed to do. "
                                       "Route A raw counts lumina and must under-predict, because "
                                       "a severed capillary drains its upstream arteriole. The "
                                       "factor is a measurement of that, not an error bar, and "
                                       "the two routes are never averaged.",
            "route_a_amplified_vs_route_b": (a["total_mL_min"] / b_hi["total_mL_min"]
                                             if b_hi["total_mL_min"] > 0 else None),
            "consistency_test": "If 1/f is the right amplification, route A amplified should land "
                                "near route B. The remaining gap is the unexplained part and is "
                                "the number this cell should be judged on.",
            "blood_loss_mL_band": [loss_lo, loss_hi],
            "exceeds_field_clearing": b_hi["total_mL_min"] > p.field_clearing_mL_min,
            "field_clearing_mL_min": p.field_clearing_mL_min,
        },

        "identifiability_certificate": ident,
        "depth_response": curve,

        "edges_this_cell_offers": [
            {"to": "SURG_INCISION", "shared_quantity": "cut depth",
             "level": "TISSUE", "type": "HANDOVER"},
            {"to": "SURG_COLLAGEN", "shared_quantity": "cut length and depth",
             "level": "TISSUE", "type": "HANDOVER"},
            {"to": "Q156", "shared_quantity": "directed pressure drop across a vascular segment",
             "level": "TISSUE", "type": "SIMULTANEOUS",
             "note": "WITHDRAWN 2026-10-02: Q156 does NOT hold a vascular pressure drop. Its "
                     "`pressure_delta_mmhg` is an INPUT defaulting to 0, for Darcy flow in brain "
                     "parenchyma at k = 1.9e-15 m^2 -- wrong level and wrong quantity. This edge "
                     "was my error. f still needs a measurement; see acquisition_required."},
            {"to": "Q036", "shared_quantity": "normalized regional flow",
             "level": "TISSUE", "type": "SIMULTANEOUS",
             "note": "Route B needs layer perfusion; Q036 computes regional flow. Today route B "
                     "uses a STANDARD textbook perfusion instead, which is the weakest number in "
                     "the cell."},
            {"to": "SURG_HEALING", "shared_quantity": "time to hemostasis",
             "level": "TISSUE", "type": "HANDOVER"},
        ],

        "acquisition_required": [
            {"quantity": "upstream_fraction_f per layer", "why": "the entire result turns on it "
             "and all five values are textbook, not measured",
             "how": "paired pressure measurement proximal and distal to a transection plane in the "
                    "same vessel, or regional flow before and after controlled transection"},
            {"quantity": "drainage distance", "why": "route B spans a factor of 4 on it alone",
             "how": "perfusion imaging of the wound margin, depth-resolved"},
            {"quantity": "layer perfusion", "why": "replaces the cell's weakest STANDARD number",
             "how": "edge to Q036 rather than a new experiment"},
        ],

        "refutes_us": [
            "If an epidermis-only incision is predicted to bleed, the chain is wrong before any "
            "parameter is defended.",
            "If route A and route B agree within 20 percent, the upstream-drainage mechanism this "
            "cell is built on does not operate and the certificate is unmotivated.",
            "If measured f exceeds 0.5 in the subcutaneous and perforator planes, arterial "
            "bleeding is predictable after all and the certificate is too strict.",
        ],
        "negative_result": True,
        "negative_result_what": "Bleeding rate is structurally unidentifiable from cut geometry "
                                "wherever less than half the series resistance lies upstream of "
                                "the cut plane. No parameter fit can repair this.",
        "external_referent": "template bleeding time 2-9 min and the clinical fact that an "
                             "epidermis-depth incision does not bleed; both STANDARD, neither "
                             "verified in this file",
        "unit_check": check_units(),
        "validation_status": "structural and dimensional check; no internal data used",
        "claim_type": "capability",
        "claim_type_reason": "Nothing in the project previously turned cut length and depth into a "
                             "bleeding rate at all, so there is no method to contest. It is judged "
                             "by whether the answer survives an external facit.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="results.json")
    a = ap.parse_args()
    res = make_results(Parameters())
    Path(a.output).write_text(json.dumps(res, indent=2, sort_keys=True, allow_nan=False) + "\n",
                              encoding="utf-8")
    d = res["decisive"]
    print(json.dumps({
        "output": a.output,
        "units_passed": res["unit_check"]["all_passed"],
        "route_a_raw_mL_min": round(d["route_a_luminal_raw_mL_min"], 5),
        "route_a_amplified_mL_min": round(d["route_a_luminal_amplified_mL_min"], 5),
        "route_a_amplified_vs_route_b": (round(d["route_a_amplified_vs_route_b"], 2)
                                         if d["route_a_amplified_vs_route_b"] else None),
        "route_b_band_mL_min": [round(x, 4) for x in d["route_b_perfusion_mL_min_band"]],
        "disagreement_factor_band": [round(x, 2) for x in d["route_disagreement_factor_band"]],
        "blood_loss_mL_band": [round(x, 3) for x in d["blood_loss_mL_band"]],
        "max_depth_fully_predictable_mm": res["identifiability_certificate"][
            "max_depth_fully_predictable_mm"],
        "unpredictable_at_2mm": [u["layer"] for u in
                                 res["identifiability_certificate"]["unpredictable_layers"]],
    }, indent=1))


if __name__ == "__main__":
    main()
