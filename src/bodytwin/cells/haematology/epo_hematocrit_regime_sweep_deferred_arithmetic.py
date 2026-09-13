"""EPO-hematocrit loop: does the RBC-production identity's top-down match generalize, and is it
hematocrit-specific?

The underlying claim was checked at a SINGLE calibration point (Hct = 45%, blood volume = 5.0 L,
MCV = 90 fL) taken from the erythropoiesis and hematopoiesis cells. This cell:
  (1) reproduces that single point bit-exactly from those two cells' result JSONs (reading the
      written results rather than importing the scripts, which execute at import time);
  (2) runs the full Hct x MCV x blood-volume grid over the physiological range, to see whether
      "RBC + neutrophil bottom-up lands within 9.6% of the Sender & Milo top-down anchor" is a
      generic feature of the identity (i.e. the tolerance is so wide that almost any Hct passes)
      or a real, Hct-specific signal;
  (3) a discriminating void floor: the pass rate inside the physiological Hct band [20, 60]% must
      exceed the pass rate over the surrounding non-physiological range, otherwise the match
      carries no information about which Hct is correct.

Reads: the erythropoiesis, hematopoiesis and renal_filtration cell results.
Writes: epo_hematocrit_regime_sweep_deferred_arithmetic.json under the cell output directory.

Gates (pre-registered before any grid was computed):
  G1 bit-exact reproduction of the cited single point (RBC P_ss, neutrophil primary rate,
     subtotal and the percentage difference against the Sender & Milo anchor), relative
     tolerance 1e-6 (1e-2 for the percentage difference, which is cited rounded to 1 decimal).
  G2 generalization: at least 50% of the physiological grid (Hct in [20, 60]%, MCV in
     [70, 105] fL, blood volume in [4.0, 6.5] L) must land within +/-15% of the anchor.
  G3 void floor: the pass rate inside the physiological Hct band must exceed the
     non-physiological-band rate by at least 1.5x, and that non-physiological rate must itself be
     <= 50%. If G3 fails, the 9.6% match is evidentially vacuous.
"""
import json
import numpy as np

import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "epo_hematocrit_regime_sweep_deferred_arithmetic")
OUT_PATH = _os.path.join(OUT_DIR, "epo_hematocrit_regime_sweep_deferred_arithmetic.json")

# ---- STEP 0: load the producer cells' results (read-only, no import/side-effect) ----
erythro = json.load(open(_os.path.join(OUT_ROOT, "erythropoiesis", "erythropoiesis_results.json")))
hemato = json.load(open(_os.path.join(OUT_ROOT, "hematopoiesis", "hematopoiesis_results.json")))
renal = json.load(open(_os.path.join(OUT_ROOT, "renal_filtration", "renal_filtration_results.json")))

HCT_REF_PCT = erythro["couples_to_siblings_readonly"]["hct_ref_pct"]          # 45.0
HB_REF_G_DL = erythro["couples_to_siblings_readonly"]["hb_ref_g_dl"]          # 15.0
MCV_FL_CITED = erythro["step4_marrow_output_geometric_derivation"]["mcv_fl"]  # 90.0
BV_L_CITED = erythro["step4_marrow_output_geometric_derivation"]["blood_volume_l"]  # 5.0
LIFESPAN_DAYS = 120.0  # LIFESPAN_DAYS_CONSENSUS, cited in the erythropoiesis cell step4/5

P_SS_CITED = erythro["step5_feedback_loop_parameters"]["p_ss_cells_per_day"]  # 2.08333e11
NEUTROPHIL_RATE_CITED = hemato["part3_amplification_decorrelated_check"]["neutrophil_leg"]["primary_rate_per_day"]
SUBTOTAL_CITED_PCTDIFF = hemato["part3_amplification_decorrelated_check"]["sendermilo2021_external_anchor"]["pctdiff_vs_rbc_neutrophil_only"]

# Sender & Milo top-down blood-cell anchor: total_cells_per_day * blood_fraction (from source code
# constants, reconfirmed here rather than re-imported: 0.33e12 * 0.90)
SENDERMILO_TOTAL_CELLS_PER_DAY = 0.33e12
SENDERMILO_BLOOD_FRACTION = 0.90
SENDERMILO_BLOOD_CELLS_PER_DAY = SENDERMILO_TOTAL_CELLS_PER_DAY * SENDERMILO_BLOOD_FRACTION  # 2.97e11

RENAL_IMPLIED_HCT = renal["step1_human_primary_anchors"]["implied_hematocrit"]
RENAL_HCT_BAND_PCT = (RENAL_IMPLIED_HCT["decade_20_29"] * 100.0, RENAL_IMPLIED_HCT["decade_30_39"] * 100.0)  # (43.03,45.04)


def total_rbc(hct_pct, mcv_fl, bv_l):
    rbc_count_millions_per_ul = hct_pct * 10.0 / mcv_fl
    return rbc_count_millions_per_ul * 1e12 * bv_l


def p_ss(hct_pct, mcv_fl, bv_l):
    return total_rbc(hct_pct, mcv_fl, bv_l) / LIFESPAN_DAYS


# ============================================================================
# GATE G1 -- bit-exact reproduction of the cited single point
# ============================================================================
p_ss_repro = p_ss(HCT_REF_PCT, MCV_FL_CITED, BV_L_CITED)
subtotal_repro = p_ss_repro + NEUTROPHIL_RATE_CITED
pctdiff_repro = (subtotal_repro - SENDERMILO_BLOOD_CELLS_PER_DAY) / SENDERMILO_BLOOD_CELLS_PER_DAY * 100.0

g1_p_ss_relerr = abs(p_ss_repro - P_SS_CITED) / P_SS_CITED
g1_pctdiff_relerr = abs(pctdiff_repro - SUBTOTAL_CITED_PCTDIFF) / abs(SUBTOTAL_CITED_PCTDIFF)
G1_PASS = bool(g1_p_ss_relerr < 1e-6 and g1_pctdiff_relerr < 1e-2)  # pctdiff cited rounded to 1dp

# ============================================================================
# GATE G2 -- full physiological grid: Hct x MCV x BV, does the +/-15% Sender&Milo match
# generalize, or was it a single lucky point?
# ============================================================================
HCT_PHYS = np.arange(20.0, 61.0, 1.0)         # anemia floor to polycythemia ceiling, 41 pts
MCV_PHYS = np.arange(70.0, 106.0, 2.0)        # microcytic to macrocytic, 18 pts
BV_PHYS = np.arange(4.0, 6.6, 0.25)           # 11 pts (Nadler-formula adult range)

HH, MM, BB = np.meshgrid(HCT_PHYS, MCV_PHYS, BV_PHYS, indexing="ij")
P_SS_GRID = total_rbc(HH, MM, BB) / LIFESPAN_DAYS
SUBTOTAL_GRID = P_SS_GRID + NEUTROPHIL_RATE_CITED  # neutrophil rate itself Hct-independent (disclosed)
PCTDIFF_GRID = (SUBTOTAL_GRID - SENDERMILO_BLOOD_CELLS_PER_DAY) / SENDERMILO_BLOOD_CELLS_PER_DAY * 100.0
WITHIN_15PCT = np.abs(PCTDIFF_GRID) <= 15.0

n_grid = WITHIN_15PCT.size
frac_grid_pass = float(np.mean(WITHIN_15PCT))
G2_PASS = bool(frac_grid_pass >= 0.50)

# joint constraint: of the grid cells that pass Sender&Milo, what fraction ALSO have Hct inside
# the renal-independent implied band?
hct_in_renal_band = (HH >= RENAL_HCT_BAND_PCT[0]) & (HH <= RENAL_HCT_BAND_PCT[1])
n_pass = int(np.sum(WITHIN_15PCT))
n_pass_and_renal = int(np.sum(WITHIN_15PCT & hct_in_renal_band))
frac_passing_cells_also_renal_consistent = (n_pass_and_renal / n_pass) if n_pass else float("nan")
# what fraction of the grid is renal-band at all (base rate, for comparison)
frac_grid_renal_band = float(np.mean(hct_in_renal_band))

# ============================================================================
# GATE G3 -- void floor: does the match have any Hct-specific discriminating power, or would
# almost any Hct pass (Sender&Milo band too wide)? Hold MCV/BV at their cited central values,
# sweep Hct over the FULL plausible-looking numeric range [1,99]%.
# ============================================================================
HCT_FULL = np.arange(1.0, 100.0, 1.0)  # 99 points, most NOT real human hematocrits
SUBTOTAL_FULL = p_ss(HCT_FULL, MCV_FL_CITED, BV_L_CITED) + NEUTROPHIL_RATE_CITED
PCTDIFF_FULL = (SUBTOTAL_FULL - SENDERMILO_BLOOD_CELLS_PER_DAY) / SENDERMILO_BLOOD_CELLS_PER_DAY * 100.0
WITHIN_15PCT_FULL = np.abs(PCTDIFF_FULL) <= 15.0

physiological_mask = (HCT_FULL >= 20.0) & (HCT_FULL <= 60.0)
nonphys_mask = ~physiological_mask

rate_physiological = float(np.mean(WITHIN_15PCT_FULL[physiological_mask]))
rate_nonphysiological = float(np.mean(WITHIN_15PCT_FULL[nonphys_mask]))

G3_ratio = (rate_physiological / rate_nonphysiological) if rate_nonphysiological > 0 else float("inf")
G3_PASS = bool(rate_physiological > 0 and (rate_nonphysiological <= 0.50) and (G3_ratio >= 1.5 or rate_nonphysiological == 0.0))

# ============================================================================
# Verdict assembly
# ============================================================================
overall_generalizes = bool(G1_PASS and G2_PASS and G3_PASS)

results = {
    "node": "EPO-hematocrit loop consistency check",
    "script": "epo_hematocrit_regime_sweep_deferred_arithmetic",
    "g1_bitexact_reproduction": {
        "p_ss_repro": p_ss_repro, "p_ss_cited": P_SS_CITED, "relerr": g1_p_ss_relerr,
        "subtotal_repro": subtotal_repro,
        "pctdiff_repro": round(pctdiff_repro, 1), "pctdiff_cited": SUBTOTAL_CITED_PCTDIFF,
        "pass": G1_PASS,
    },
    "g2_full_grid_generalization": {
        "grid_shape": list(P_SS_GRID.shape), "n_grid_cells": n_grid,
        "hct_range_pct": [float(HCT_PHYS.min()), float(HCT_PHYS.max())],
        "mcv_range_fl": [float(MCV_PHYS.min()), float(MCV_PHYS.max())],
        "bv_range_l": [float(BV_PHYS.min()), float(BV_PHYS.max())],
        "frac_grid_within_15pct_sendermilo": frac_grid_pass,
        "gate_threshold": 0.50,
        "pass": G2_PASS,
        "n_pass_cells": n_pass,
        "frac_passing_cells_also_renal_hct_consistent": frac_passing_cells_also_renal_consistent,
        "frac_grid_that_is_renal_hct_band_at_all_baserate": frac_grid_renal_band,
    },
    "g3_void_floor_hct_specificity": {
        "rate_within_15pct_physiological_hct_20_60": rate_physiological,
        "rate_within_15pct_nonphysiological_hct_1_19_and_61_99": rate_nonphysiological,
        "ratio_physiological_over_nonphysiological": G3_ratio,
        "pass_discriminating": G3_PASS,
        "interpretation": (
            "Sender&Milo's 15pct band DOES discriminate physiological from non-physiological Hct"
            if G3_PASS else
            "Sender&Milo's 15pct band does NOT discriminate -- most numeric Hct values in [1,99] "
            "would ALSO pass, so the cited 9.6pct match carries little Hct-specific information "
            "(an evidentially-weak leg, analogous to a scrambled-pairing null landing on the same number)"
        ),
    },
    "renal_independent_hct_anchor_band_pct": list(RENAL_HCT_BAND_PCT),
    "cited_hct_ref_pct": HCT_REF_PCT,
    "cited_hct_inside_renal_band": bool(RENAL_HCT_BAND_PCT[0] <= HCT_REF_PCT <= RENAL_HCT_BAND_PCT[1]),
    "overall_verdict": {
        "claim": "C: the identity's Sender&Milo match at the cited single point generalizes across "
                 "physiological (Hct,MCV,BV) space AND is Hct-specific (not a wide-tolerance artifact)",
        "G1_bitexact": G1_PASS, "G2_generalizes_grid": G2_PASS, "G3_hct_specific_not_vacuous": G3_PASS,
        "overall_pass_C": overall_generalizes,
    },
}

import os
os.makedirs(OUT_DIR, exist_ok=True)
with open(OUT_PATH, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results["overall_verdict"], indent=2))
print(json.dumps(results["g2_full_grid_generalization"], indent=2))
print(json.dumps(results["g3_void_floor_hct_specificity"], indent=2))
print(f"cited Hct={HCT_REF_PCT}% inside renal-implied band {RENAL_HCT_BAND_PCT}: "
      f"{results['cited_hct_inside_renal_band']}")
