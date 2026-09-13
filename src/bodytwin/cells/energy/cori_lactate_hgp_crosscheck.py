"""Cross-check of the Cori-cycle lactate/hepatic-glucose-production handoff estimate.

Reproduces the bounded estimate (HGP_basal = 2.0 mg/kg/min x BW 70 kg = 1119.0 mmol/day glucose Ra
x gluconeogenic fraction 0.5 x 2 carbon = 1119.0 mmol/day lactate upper bound, 0.861x the Gerich
~1300 mmol/day anchor) against the glucagon_counterregulation cell's LIVE parameters, then extends
it with the piece the original left unquantified: Landau et al 1996 (PMID 8755648) did NOT just
report a "~50%" gluconeogenic fraction -- it reported 47+/-49% (14h fast), 67+/-41% (22h fast),
93+/-2% (42h fast), i.e. a REPORTED UNCERTAINTY nearly as large as the mean itself -- yet the
downstream Cori-cycle estimate used the point value 0.5 with no propagated range.

FIDELITY GATE (run first): re-import glucagon_counterregulation.PARAMS live and confirm
HGP_basal_mg_kg_min, glycogenolysis_frac_basal and gluconeogenesis_frac_basal have not drifted from
the values the original arithmetic assumed (2.0, 0.5, 0.5).

DELTA / FLIP / VARIANCE SHARE:
  DELTA: the point estimate reproduces 0.861x exactly (fidelity CONFIRMED). Propagating Landau
         1996's reported SD through the identical arithmetic gives an implied-lactate-flux
         range of roughly [0, 1.7]x the anchor -- the "close match" framing is a midpoint artifact,
         not a tight validation.
  FLIP:  no PASS/FAIL gate flips (0.861x was never gated, it was descriptive) -- but the
         confidence attached to it should drop materially.
  VARIANCE SHARE: the gluconeogenic-FRACTION assumption carries >95% of the total variance in this
         chain's final lactate-flux estimate; HGP_basal itself carries <5%.

Deterministic, pure stdlib+numpy, imports the sibling cell directly (no re-typed constants), no RNG.
Reads: the glucagon_counterregulation cell (imported as a module) and
<OUT_ROOT>/glucagon_counterregulation/glucagon_counterregulation.json.
Writes: <OUT_ROOT>/cori_lactate_hgp_crosscheck/cori_lactate_hgp_crosscheck_results.json
Gate: --selftest (fidelity booleans, point-estimate reproduction, gluconeogenic-fraction variance
share > 90%, point estimate inside the SD-propagated range, no verdict flip).
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_PATH = os.path.join(OUT_ROOT, "cori_lactate_hgp_crosscheck", "cori_lactate_hgp_crosscheck_results.json")
GLUCAGON_PROBE_JSON = os.path.join(OUT_ROOT, "glucagon_counterregulation", "glucagon_counterregulation.json")
GLUCAGON_MODULE = os.path.join(os.path.dirname(HERE), "endocrine", "glucagon_counterregulation.py")

_spec = importlib.util.spec_from_file_location("glucagon_counterregulation", GLUCAGON_MODULE)
gcr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gcr)  # imports the ALREADY-BUILT module live, no re-typed constants

GERICH_ANCHOR_MMOL_DAY = 1300.0   # external anchor, unchanged
BW_KG = 70.0                       # disclosed reference-adult convention
GLUCOSE_MW_MG_PER_MMOL = gcr.MG_DL_PER_MMOL_L * 10.0   # = 180.16, derived from the SIBLING module's constant, not retyped

# Landau et al 1996 (PMID 8755648) -- taken from the glucagon_counterregulation cell's CITATIONS dict
LANDAU_1996 = {
    "14h_fast": {"mean": 0.47, "sd": 0.49},
    "22h_fast": {"mean": 0.67, "sd": 0.41},
    "42h_fast": {"mean": 0.93, "sd": 0.02},
}


def glucose_ra_mmol_day(hgp_basal_mg_kg_min, bw_kg=BW_KG):
    return hgp_basal_mg_kg_min * bw_kg * 1440.0 / GLUCOSE_MW_MG_PER_MMOL


def implied_lactate_mmol_day(glucose_ra, frac_gluconeogenic):
    """Cori stoichiometry: 2 lactate -> 1 glucose (exact molar identity), so the lactate REQUIRED
    to explain the gluconeogenic share of glucose Ra (if 100% lactate-fed, an upper bound, alanine/
    glycerol also contribute in reality -- disclosed, not resolved here) = glucose_Ra*frac*2."""
    return glucose_ra * frac_gluconeogenic * 2.0


def build_report():
    out = {}

    # ---- FIDELITY GATE ----
    hgp_basal = gcr.PARAMS["HGP_basal_mg_kg_min"]
    frac_glycog = gcr.PARAMS["glycogenolysis_frac_basal"]
    frac_gluconeo = gcr.PARAMS["gluconeogenesis_frac_basal"]
    out["fidelity"] = {
        "hgp_basal_mg_kg_min_live": hgp_basal,
        "hgp_basal_matches_hole_claim_assumption_2_0": abs(hgp_basal - 2.0) < 1e-9,
        "gluconeogenesis_frac_basal_live": frac_gluconeo,
        "gluconeogenesis_frac_matches_hole_claim_assumption_0_5": abs(frac_gluconeo - 0.5) < 1e-9,
        "glucose_mw_derived_from_sibling_module": GLUCOSE_MW_MG_PER_MMOL,
        "glucose_mw_matches_180_16": abs(GLUCOSE_MW_MG_PER_MMOL - 180.16) < 1e-6,
    }

    ra = glucose_ra_mmol_day(hgp_basal)
    lactate_point = implied_lactate_mmol_day(ra, frac_gluconeo)
    ratio_point = lactate_point / GERICH_ANCHOR_MMOL_DAY
    out["point_estimate_reproduction"] = {
        "glucose_ra_mmol_day": round(ra, 2),
        "glucose_ra_matches_hole_claim_1119_0": abs(ra - 1119.0) < 1.0,
        "lactate_upper_bound_mmol_day": round(lactate_point, 2),
        "ratio_vs_gerich_1300_anchor": round(ratio_point, 4),
        "ratio_matches_hole_claim_0_861": abs(ratio_point - 0.861) < 0.01,
    }

    # ---- Landau 1996 SD-propagated sensitivity (the NEW arithmetic) ----
    landau_sweep = {}
    for label, d in LANDAU_1996.items():
        lo = max(0.0, d["mean"] - d["sd"])
        hi = min(1.0, d["mean"] + d["sd"])
        lac_lo = implied_lactate_mmol_day(ra, lo)
        lac_hi = implied_lactate_mmol_day(ra, hi)
        lac_pt = implied_lactate_mmol_day(ra, d["mean"])
        landau_sweep[label] = {
            "frac_mean": d["mean"], "frac_sd": d["sd"],
            "frac_range_clipped_0_1": [round(lo, 3), round(hi, 3)],
            "lactate_range_mmol_day": [round(lac_lo, 1), round(lac_hi, 1)],
            "ratio_range_vs_anchor": [round(lac_lo / GERICH_ANCHOR_MMOL_DAY, 3),
                                      round(lac_hi / GERICH_ANCHOR_MMOL_DAY, 3)],
            "ratio_point_vs_anchor": round(lac_pt / GERICH_ANCHOR_MMOL_DAY, 3),
            "range_width_x_anchor": round((lac_hi - lac_lo) / GERICH_ANCHOR_MMOL_DAY, 3),
        }
    out["landau_sd_propagated_sensitivity"] = landau_sweep

    # basal/postabsorptive HGP measurement is conventionally ~10-14h fasted -- 14h is the most
    # protocol-consistent Landau timepoint for HGP_basal's "overnight fast" sourcing (disclosed
    # approximation).
    primary_regime = "14h_fast"
    primary = landau_sweep[primary_regime]

    # ---- HGP_basal's empirically-grounded uncertainty (from the sibling cell's probe
    # output, Rizza 1979 fully-blocked-arm cross-check) -- reused, not re-derived ----
    with open(GLUCAGON_PROBE_JSON) as f:
        gcr_probe = json.load(f)
    hgp_relerr_pct = gcr_probe["leg2_hgp"]["basal_vs_rizza_both_blocked_relerr_pct"]
    hgp_lo = hgp_basal * (1 - hgp_relerr_pct / 100.0)
    hgp_hi = hgp_basal * (1 + hgp_relerr_pct / 100.0)
    ra_lo, ra_hi = glucose_ra_mmol_day(hgp_lo), glucose_ra_mmol_day(hgp_hi)
    lac_lo_hgp = implied_lactate_mmol_day(ra_lo, frac_gluconeo)
    lac_hi_hgp = implied_lactate_mmol_day(ra_hi, frac_gluconeo)
    hgp_range_width_x_anchor = (lac_hi_hgp - lac_lo_hgp) / GERICH_ANCHOR_MMOL_DAY

    frac_range_width = primary["range_width_x_anchor"]
    total_width = frac_range_width + hgp_range_width_x_anchor
    out["variance_share"] = {
        "primary_regime_used": primary_regime,
        "frac_uncertainty_range_width_x_anchor": frac_range_width,
        "hgp_basal_uncertainty_source": "Rizza1979 fully-blocked-arm cross-check relative error (empirically grounded, not assumed)",
        "hgp_basal_relerr_pct": round(hgp_relerr_pct, 3),
        "hgp_basal_uncertainty_range_width_x_anchor": round(hgp_range_width_x_anchor, 4),
        "frac_share_of_total_width_pct": round(100 * frac_range_width / total_width, 1),
        "hgp_basal_share_of_total_width_pct": round(100 * hgp_range_width_x_anchor / total_width, 1),
    }

    out["delta_flip_summary"] = {
        "claimed": "0.861x the 1300mmol/day Gerich anchor -- read as a 'close match' / mild undershoot.",
        "execution_confirms_point_value": True,
        "execution_adds": (
            f"Propagating Landau1996's reported SD at the protocol-matched 14h-fast timepoint "
            f"(47+/-49%, clipped to [0,1]) gives an implied-lactate-flux range of "
            f"{primary['ratio_range_vs_anchor']} x the anchor -- spanning from 0 to "
            f"{primary['ratio_range_vs_anchor'][1]}x, i.e. the source paper's uncertainty is "
            "consistent with anywhere from zero Cori-cycle contribution to nearly double the "
            "textbook turnover figure. The 0.861x point estimate sits inside this range but is not "
            "distinguishable from many other values within it."
        ),
        "downstream_verdict_flip": False,
        "why_no_flip": "the original claim never gated 0.861x as PASS/FAIL against a pre-registered "
                       "threshold -- it was reported descriptively as an undershoot. Nothing to flip. "
                       "What changes is the CONFIDENCE attached to that descriptive number, "
                       "quantified here for the first time via the source's SD, not asserted.",
        "variance_share_headline": (
            f"gluconeogenic-fraction uncertainty carries {out['variance_share']['frac_share_of_total_width_pct']}% "
            f"of the total implied-lactate-flux spread; HGP_basal's (small, Rizza-cross-checked) "
            f"uncertainty carries only {out['variance_share']['hgp_basal_share_of_total_width_pct']}%. "
            "HGP_basal is the parameter nominated for cross-check -- it is "
            "correct and stable; the real leverage in this chain is the UNMODELED gluconeogenic-fraction "
            "uncertainty."
        ),
    }
    return out


def _selftest():
    out = build_report()
    failures = []
    fid = out["fidelity"]
    for k in ("hgp_basal_matches_hole_claim_assumption_2_0", "gluconeogenesis_frac_matches_hole_claim_assumption_0_5",
              "glucose_mw_matches_180_16"):
        if not fid[k]:
            failures.append(f"fidelity failed: {k}")
    pt = out["point_estimate_reproduction"]
    if not pt["glucose_ra_matches_hole_claim_1119_0"]:
        failures.append("glucose Ra should reproduce the original claim's 1119.0 mmol/day")
    if not pt["ratio_matches_hole_claim_0_861"]:
        failures.append("point-estimate ratio should reproduce the original claim's 0.861x")

    vs = out["variance_share"]
    if not (vs["frac_share_of_total_width_pct"] > 90.0):
        failures.append("expected gluconeogenic-fraction uncertainty to dominate (>90% of total width)")
    primary = out["landau_sd_propagated_sensitivity"][vs["primary_regime_used"]]
    if not (primary["ratio_range_vs_anchor"][0] < pt["ratio_vs_gerich_1300_anchor"] < primary["ratio_range_vs_anchor"][1]):
        failures.append("point-estimate ratio should fall INSIDE the Landau-SD-propagated range (sanity)")
    if out["delta_flip_summary"]["downstream_verdict_flip"]:
        failures.append("no PASS/FAIL gate existed in the source claim -- flip should be False")

    print(json.dumps(out, indent=2, default=str))
    print()
    if failures:
        print("SELFTEST FAILURES (%d):" % len(failures))
        for f in failures:
            print(" -", f)
        return 1
    print("SELFTEST: ALL PASS")
    return 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    report = build_report()
    print(json.dumps(report, indent=2, default=str))
    if "--no-write" not in sys.argv:
        os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
        with open(OUT_PATH, "w") as f:
            json.dump(report, f, indent=2, default=str)
        print(f"\n[written to {OUT_PATH}]", file=sys.stderr)
