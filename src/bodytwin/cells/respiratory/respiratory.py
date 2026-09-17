"""RESPIRATORY / VENTILATION-O2 layer: 0-D steady-state conversion of an already-computed whole-body
metabolic rate into oxygen uptake (VO2) and minute ventilation (VE), for REST and WALKING, compared
against published human exercise-physiology anchors.

METHOD (0-D steady-state stoichiometric conversion -- explicitly NOT a breath-dynamics model):
  1. VO2 from metabolic rate via the caloric equivalent of oxygen, E_O2 = 20.9 kJ/L O2, cross-checked
     (Step 2) against the RQ-boundary values (carbohydrate 21.13 kJ/L at RQ=1.0, fat 19.62 kJ/L at
     RQ=0.7): solving the linear-interpolation formula for RQ gives RQ=~0.95, physiologically
     reasonable for moderate, carbohydrate-leaning exercise metabolism.
       VO2 [mL/(kg*min)] = P [W/kg] * 60 / E_O2 [kJ/L]   (derivation in vo2_ml_per_kg_min())
  2. VE from VO2 via the ventilatory equivalent for oxygen VE/VO2, swept low/mid/high = 25/27.5/30:
       VE [L/min] = VO2 [L/min] * (VE/VO2)
  3. Applied to FOUR metabolic-rate inputs read from the metabolic_cost cell: REST (basal), and the
     three walking configurations (Umberger-primary, Bhargava-primary, combined-corrected), carried
     forward end-to-end rather than collapsed to one.
  4. FORCED, non-degeneracy adversary: is this just a linear rescaling that makes any plausible input
     look reasonable? Tested via (a) a broad sweep of the input metabolic rate with an
     analytical-vs-numerical derivative match (Step 8), and (b) whether the published-VO2-range gate
     (Step 5) actually discriminates between the three configurations.

SCOPE (stated before any number is computed): steady state only -- no breath-by-breath dynamics, no
VO2 kinetics time constant, no respiratory-rate/tidal-volume decomposition (only their product), no
dead-space or alveolar-ventilation separation, no diffusion limitation, no altitude/inspired-O2
coupling; generic constants (E_O2, VE/VO2), not subject-measured RER or lung function. Every number
is a linear re-expression of the SAME upstream metabolic rate the thermoregulation cell also consumes:
the agreement of three different external literature bases on which configuration is more consistent
is convergent external evidence, not three independent re-derivations of the upstream number.

CITATIONS (PMIDs/DOIs verified against NCBI eutils, not recalled):
  [1] Ainsworth BE et al (2011), "2011 Compendium of Physical Activities: a second update of codes
      and MET values", Med Sci Sports Exerc 43(8):1575-81, PMID 21681120,
      DOI 10.1249/MSS.0b013e31821ece12. The Compendium's public companion database gives, for
      walking: code 17255 "Walking, self-selected speed, indoor track or outdoors, firm surface" =
      4.0 METs; code 17190 "Walking, 2.8 to 3.4 mph, level, moderate pace, firm surface" = 3.8 METs;
      code 17170 "Walking, 2.5 mph, firm, level surface" = 3.0 METs. Via 1 MET = 3.5 mL O2/kg/min
      these are VO2 = 14.0, 13.3, 10.5 mL/kg/min.
  [2] Waters RL, Mulroy S (1999), "The energy expenditure of normal and pathologic gait", Gait
      Posture 9(3):207-31, PMID 10575082, DOI 10.1016/s0966-6362(99)00009-0 -- cited for topical
      grounding (paywalled; no number re-extracted here).
  [3] Browning RC, Kram R (2005), Obes Res 13(5):891-9, PMID 15919843, DOI 10.1038/oby.2005.103;
      Browning RC, Baker EA, Herron JA, Kram R (2006), J Appl Physiol 100(2):390-8, PMID 16210434,
      DOI 10.1152/japplphysiol.00767.2005 -- source of the preferred/self-selected walking speed
      ~1.40-1.47 m/s context figure.
  [4] Wasserman K, Whipp BJ, Koyl SN, Beaver WL (1973), J Appl Physiol 35(2):236-43, PMID 4723033,
      DOI 10.1152/jappl.1973.35.2.236 -- the classical reference establishing the
      ventilatory-equivalent (VE/VO2, VE/VCO2) method; the 25-30 moderate-exercise figure used here
      is the textbook-consensus value, not re-extracted from its primary text.
  [5] ATS/ACCP Statement on Cardiopulmonary Exercise Testing (2003), Am J Respir Crit Care Med
      167(2):211-77, PMID 12524257, DOI 10.1164/rccm.167.2.211 -- standard clinical reference for
      CPET normal values including ventilatory equivalents.
  [6] Peronnet F, Massicotte D (1991), "Table of nonprotein respiratory quotient: an update", Can J
      Sport Sci 16(1):23-9, PMID 1645211 (no DOI on record) -- provenance of the
      caloric-equivalent-of-oxygen concept.
  [7] "Indirect calorimetry" reference article (textbook-grade, flagged): 21.13 kJ/L for
      carbohydrate oxidation (RQ=1.0), 19.62 kJ/L for fat (RQ=0.7), plus the linear-interpolation
      formula between them, used to machine-check the 20.9 kJ/L constant.
  [8] "Minute ventilation" reference article (textbook-grade, flagged): resting minute volume about
      5-8 L/min; light activity around 12 L/min; moderate exercise between 40 and 60 L/min.

READS: <BODYTWIN_OUT>/metabolic_cost/metabolic_cost_results.json
WRITES: <BODYTWIN_OUT>/respiratory/respiratory_results.json
GATE: overall_pass = all ten gates in the final GATES block; exit 0 on pass, 2 otherwise.

RESULT ENVELOPE: this cell is a consumer of the shared metabolic-cost result
envelope. Before any number is read, `consumer_preflight_v1.preflight_metabolic_cost` checks
provenance / region / frame / time / regime / staleness / producer pin on
`metabolic_cost_results.json`, and labels the result scientific evidence, an explicit synthetic
demo (never promoted), or a refusal (exit 2, no numbers used). The shared quantity contract
(`METABOLIC_COST_QUANTITY_SPECS`, nine keys) is a SUPERSET of what this cell reads (six keys:
total body mass, the four headline w/kg values and the combined tendon+mass correction); the
producer always emits all nine, so the full contract is reused unchanged. Its own ten gates and
thresholds are NOT touched by the migration.
"""
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

# --------------------------------------------------------------------------- I1 preflight --
# Result-envelope preflight wiring: make the shared framework modules importable when this cell is
# run standalone (python src/bodytwin/cells/respiratory/respiratory.py). The import MUST NOT be
# silently skipped: a missing framework module is a hard error, never a downgrade to the legacy
# unchecked cached-result read.
_SRC_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if _SRC_ROOT not in sys.path:
    sys.path.insert(0, _SRC_ROOT)
try:
    from bodytwin.framework.consumer_preflight_v1 import preflight_metabolic_cost
    from bodytwin.framework.result_envelope_v1 import Verdict, get_dotted
except ImportError as _exc:  # pragma: no cover - exercised only on a broken checkout
    raise ImportError(
        "respiratory.py requires bodytwin.framework."
        "result_envelope_v1 and bodytwin.framework.consumer_preflight_v1; ensure the "
        f"src root is on sys.path (expected {_SRC_ROOT!r}). Original error: {_exc}"
    ) from _exc

# --------------------------------------------------------------------------- paths / consts --
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
METCOST_JSON = _os.path.join(OUT_ROOT, "metabolic_cost", "metabolic_cost_results.json")
OUT_DIR = _os.path.join(OUT_ROOT, "respiratory")

# ---- verified-live physiological constants (see docstring CITATIONS for the fetch trail) -----
E_O2_KJ_PER_L = 20.9         # this task's given caloric equivalent of O2 (moderate exercise,
                             # mixed/carb-leaning fuel) -- cross-checked below, not just asserted.
CARB_E_O2_KJ_PER_L = 21.13   # Wikipedia "Indirect calorimetry" [7]: pure carbohydrate, RQ=1.0
FAT_E_O2_KJ_PER_L = 19.62    # Wikipedia "Indirect calorimetry" [7]: pure fat, RQ=0.7
ONE_MET_ML_O2_PER_KG_MIN = 3.5   # universal DEFINITION of 1 MET (not a separately-cited number)

VE_VO2_LOW, VE_VO2_MID, VE_VO2_HIGH = 25.0, 27.5, 30.0   # this task's given ventilatory-equivalent
                                                          # range at moderate exercise [4][5]

# MET/Compendium-derived walking-VO2 anchor [1] (pacompendium.com, live-fetched codes):
COMPENDIUM_SELF_SELECTED_METS = 4.0   # code 17255 "self-selected speed... firm surface"
COMPENDIUM_MODERATE_2P8_3P4MPH_METS = 3.8   # code 17190
COMPENDIUM_2P5MPH_METS = 3.0   # code 17170
PUBLISHED_WALKING_VO2_LOW_ML_KG_MIN, PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN = 12.0, 15.0  # this
                                              # task's stated preferred-speed-walking anchor

# Wikipedia "Minute ventilation" [8] textbook-grade bands:
WIKI_VE_REST_LOW_L_MIN, WIKI_VE_REST_HIGH_L_MIN = 5.0, 8.0
WIKI_VE_LIGHT_L_MIN = 12.0
WIKI_VE_MODERATE_LOW_L_MIN, WIKI_VE_MODERATE_HIGH_L_MIN = 40.0, 60.0

# ---- pre-registered thresholds (set from the external literature, BEFORE reading whether the
# reference body's specific number passes) -----------------------------------------------------------
REST_VO2_TOLERANCE_FRAC = 0.10       # rest VO2 should be within +-10% of the 1-MET DEFINITION
                                     # (3.5 mL/kg/min) -- a self-consistency/reasonableness check
                                     # on E_O2's choice, not an independent empirical discovery.
SENSITIVITY_SWEEP_LOW_W, SENSITIVITY_SWEEP_HIGH_W = 90.0, 1200.0   # near-rest to near-maximal
                                                                    # sustained human exercise,
                                                                    # matching thermoregulation.py's
                                                                    # own upper bound for
                                                                    # comparability.


def load_metabolic_cost():
    if not os.path.exists(METCOST_JSON):
        return None
    with open(METCOST_JSON, "r", encoding="utf-8") as f:
        return json.load(f)


def _sha256_file(path):
    """sha256 of an input file's bytes -- the ``current_input_sha256`` the consumer uses."""
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _load_registry_entry(path, result_name="metabolic_cost_results.json"):
    """Load the I1 registry file and return the entry for the metabolic-cost result."""
    from bodytwin.framework.consumer_preflight_v1 import load_registry, registry_entry
    return registry_entry(load_registry(path), result_name)


def _write_refused_report(mode, verdict, input_path):
    """Write the I1 refusal record and print each failure prefixed ``REFUSED:``; return 2."""
    record = {
        "mode": mode,
        "failures": list(verdict.failures),
        "flags": list(verdict.flags),
        "input": input_path,
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    # A REFUSE MUST NOT leave a stale ACCEPT behind: remove any prior accepted result so a
    # downstream reader keyed on respiratory_results.json cannot pick up an older run.
    stale_results_path = os.path.join(OUT_DIR, "respiratory_results.json")
    if os.path.exists(stale_results_path):
        os.remove(stale_results_path)
    refused_path = os.path.join(OUT_DIR, "respiratory_refused.json")
    with open(refused_path, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2, default=str)
    for failure in verdict.failures:
        print(f"REFUSED: {failure}")
    print(f"REFUSED: wrote {refused_path}")
    return 2


def vo2_ml_per_kg_min(power_w_per_kg, e_o2_kj_per_l=E_O2_KJ_PER_L):
    """VO2 [mL/(kg*min)] from metabolic power [W/kg]:
    P [W/kg] = P [J/(s*kg)]  ->  *60 -> [J/(min*kg)]  -> /1000 -> [kJ/(min*kg)]
    VO2 [L/(min*kg)] = that / e_o2_kj_per_l   ->  *1000 -> [mL/(kg*min)]
    Net: VO2[mL/kg/min] = P[W/kg] * 60 / e_o2_kj_per_l."""
    return power_w_per_kg * 60.0 / e_o2_kj_per_l


def implied_rq_from_e_o2(e_o2_kj_per_l, carb=CARB_E_O2_KJ_PER_L, fat=FAT_E_O2_KJ_PER_L):
    """Invert Wikipedia's linear RQ-interpolation formula (e_o2 = fat + (carb-fat)*x, x=(RQ-0.7)/0.3)
    to recover the RQ this task's E_O2 choice implicitly assumes -- a machine cross-check, not an assertion."""
    x = (e_o2_kj_per_l - fat) / (carb - fat)
    return 0.7 + 0.3 * x


def ve_l_per_min(vo2_l_per_min, ve_vo2_ratio):
    return vo2_l_per_min * ve_vo2_ratio


def main():
    if not os.path.exists(METCOST_JSON):
        print(f"FAIL: required input missing: {METCOST_JSON} -- run the metabolic_cost cell first.")
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    report = {}

    print("=" * 78)
    print("STEP 1/8 -- load the already-computed metabolic rate (no re-solve)")
    print("=" * 78)
    # --- I1 preflight: decide evidence class / REFUSE before any number is used -----------------
    mode = os.environ.get("BODYTWIN_EVIDENCE_MODE", "auto").strip().lower()
    registry_entry = None
    _registry_path = os.environ.get("BODYTWIN_I1_REGISTRY")
    if _registry_path:
        registry_entry = _load_registry_entry(_registry_path)
    raw = load_metabolic_cost()
    try:
        payload, verdict = preflight_metabolic_cost(
            raw, mode=mode, registry_entry=registry_entry,
            current_input_sha256=_sha256_file(METCOST_JSON), producer_path=None,
        )
    except ValueError as exc:
        # An unknown BODYTWIN_EVIDENCE_MODE is a refusal, not a traceback.
        return _write_refused_report(
            mode,
            Verdict("REFUSE", (f"invalid evidence mode: {exc}",), (), "refused"),
            METCOST_JSON,
        )
    if verdict.state == "REFUSE":
        return _write_refused_report(mode, verdict, METCOST_JSON)
    mc = payload
    if verdict.state == "ACCEPT_SYNTHETIC_DEMO":
        print("EVIDENCE: SYNTHETIC DEMO -- not a scientific result")
    elif verdict.state == "ACCEPT_SCIENTIFIC":
        print(f"EVIDENCE: SCIENTIFIC provenance={verdict.evidence_class}")
    mass_kg = get_dotted(mc, "muscle_mass.total_body_mass_kg")
    umb_gross_w_per_kg = get_dotted(mc, "headline.umberger2010_gross_w_per_kg")
    umb_net_w_per_kg = get_dotted(mc, "headline.umberger2010_net_w_per_kg")
    bhg_gross_w_per_kg = get_dotted(mc, "headline.bhargava2004_gross_w_per_kg")
    bhg_net_w_per_kg = get_dotted(mc, "headline.bhargava2004_net_w_per_kg")
    basal_w_per_kg = umb_gross_w_per_kg - umb_net_w_per_kg
    basal_w_per_kg_check = bhg_gross_w_per_kg - bhg_net_w_per_kg
    basal_consistent = abs(basal_w_per_kg - basal_w_per_kg_check) < 1e-9
    combined_corrected_net_w_per_kg = get_dotted(
        mc, "sensitivity.combined_tendon_and_mass_correction_w_per_kg")
    combined_corrected_gross_w_per_kg = combined_corrected_net_w_per_kg + basal_w_per_kg

    configs_w_per_kg = {
        "rest_basal": basal_w_per_kg,
        "umberger_primary": umb_gross_w_per_kg,
        "bhargava_primary": bhg_gross_w_per_kg,
        "combined_corrected": combined_corrected_gross_w_per_kg,
    }
    print(f"Subject mass={mass_kg:.2f} kg  basal={basal_w_per_kg:.4f} W/kg  "
          f"(cross-model self-consistency umb vs bhg: {'PASS' if basal_consistent else 'FAIL'})")
    for name, w_kg in configs_w_per_kg.items():
        print(f"  {name:20s} {w_kg:7.4f} W/kg -> {w_kg*mass_kg:7.1f} W total")
    print("(combined_corrected = this task brief's named '~400W gross rate': "
          f"{combined_corrected_gross_w_per_kg*mass_kg:.1f} W)")

    print("\n" + "=" * 78)
    print("STEP 2/8 -- E_O2 caloric-equivalent choice, machine cross-checked against RQ bounds")
    print("=" * 78)
    implied_rq = implied_rq_from_e_o2(E_O2_KJ_PER_L)
    rq_plausible = 0.70 <= implied_rq <= 1.05
    print(f"E_O2 = {E_O2_KJ_PER_L} kJ/L O2 (this task's given constant). Wikipedia 'Indirect "
          f"calorimetry' RQ bounds: fat(RQ=0.7)={FAT_E_O2_KJ_PER_L} kJ/L, carb(RQ=1.0)="
          f"{CARB_E_O2_KJ_PER_L} kJ/L. Inverting the linear-interpolation formula for RQ: "
          f"implied RQ = {implied_rq:.3f}  ({'PLAUSIBLE for moderate carb-leaning exercise metabolism, in [0.70,1.05]' if rq_plausible else 'IMPLAUSIBLE -- outside the physiological RQ range'})")
    print(f"Gate: {'PASS' if rq_plausible else 'FAIL'}")

    print("\n" + "=" * 78)
    print("STEP 3/8 -- VO2 (mL/kg/min, L/min) for REST and all 3 WALKING configs")
    print("=" * 78)
    vo2_ml_kg_min = {name: vo2_ml_per_kg_min(w_kg) for name, w_kg in configs_w_per_kg.items()}
    vo2_l_min = {name: v * mass_kg / 1000.0 for name, v in vo2_ml_kg_min.items()}
    for name in configs_w_per_kg:
        print(f"  {name:20s} VO2 = {vo2_ml_kg_min[name]:7.3f} mL/kg/min  = {vo2_l_min[name]*1000:7.1f} mL/min "
              f"= {vo2_l_min[name]:.4f} L/min")

    print("\n" + "=" * 78)
    print("STEP 4/8 -- external anchor #1 (REST): 1-MET DEFINITION cross-check")
    print("=" * 78)
    rest_vo2_ratio_to_1met = vo2_ml_kg_min["rest_basal"] / ONE_MET_ML_O2_PER_KG_MIN
    rest_vo2_1met_ok = abs(rest_vo2_ratio_to_1met - 1.0) <= REST_VO2_TOLERANCE_FRAC
    print(f"Rest VO2 (from the basal rate={basal_w_per_kg:.3f} W/kg, via E_O2={E_O2_KJ_PER_L} kJ/L) = "
          f"{vo2_ml_kg_min['rest_basal']:.3f} mL/kg/min  vs  1 MET DEFINITION = {ONE_MET_ML_O2_PER_KG_MIN} mL/kg/min  "
          f"(ratio={rest_vo2_ratio_to_1met:.3f}, {'within' if rest_vo2_1met_ok else 'OUTSIDE'} "
          f"+-{REST_VO2_TOLERANCE_FRAC*100:.0f}%)")
    print("(NOTE: this is primarily a self-consistency check on this pipeline's arithmetic and "
          "the E_O2 constant's reasonableness -- 1 MET is a DEFINITION, not an independent empirical "
          "discovery -- reported plainly as such, not oversold.)")
    print(f"Gate: {'PASS' if rest_vo2_1met_ok else 'FAIL'}")

    print("\n" + "=" * 78)
    print("STEP 5/8 -- external anchor #2 (WALKING): MET/Compendium of Physical Activities [1]")
    print("=" * 78)
    compendium_vo2 = {
        "self_selected_speed_4p0_MET": COMPENDIUM_SELF_SELECTED_METS * ONE_MET_ML_O2_PER_KG_MIN,
        "moderate_2p8_3p4mph_3p8_MET": COMPENDIUM_MODERATE_2P8_3P4MPH_METS * ONE_MET_ML_O2_PER_KG_MIN,
        "2p5mph_3p0_MET": COMPENDIUM_2P5MPH_METS * ONE_MET_ML_O2_PER_KG_MIN,
    }
    for name, v in compendium_vo2.items():
        print(f"  Compendium {name:28s} -> VO2 = {v:.2f} mL/kg/min")
    print(f"Task-brief published anchor range: [{PUBLISHED_WALKING_VO2_LOW_ML_KG_MIN}, "
          f"{PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN}] mL/kg/min -- Compendium's 'self-selected "
          f"speed' code lands at {compendium_vo2['self_selected_speed_4p0_MET']:.1f} mL/kg/min, "
          f"inside/near-center of that range: independently corroborates the task's stated anchor.")
    walking_gate_results = {}
    for name in ("umberger_primary", "bhargava_primary", "combined_corrected"):
        v = vo2_ml_kg_min[name]
        in_range = PUBLISHED_WALKING_VO2_LOW_ML_KG_MIN <= v <= PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN
        ratio_to_range_mid = v / ((PUBLISHED_WALKING_VO2_LOW_ML_KG_MIN + PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN) / 2)
        walking_gate_results[name] = in_range
        print(f"  {name:20s} VO2={v:7.3f} mL/kg/min  in-published-range={in_range}  "
              f"ratio-to-range-midpoint={ratio_to_range_mid:.2f}")

    print("\n" + "=" * 78)
    print("STEP 6/8 -- VE (minute ventilation) via VE/VO2 sweep (25/27.5/30), REST + WALKING")
    print("=" * 78)
    ve_l_min = {}
    for name, vo2 in vo2_l_min.items():
        ve_l_min[name] = {
            "low_25": ve_l_per_min(vo2, VE_VO2_LOW),
            "mid_27p5": ve_l_per_min(vo2, VE_VO2_MID),
            "high_30": ve_l_per_min(vo2, VE_VO2_HIGH),
        }
        print(f"  {name:20s} VE = {ve_l_min[name]['low_25']:6.2f} (VE/VO2=25)  "
              f"{ve_l_min[name]['mid_27p5']:6.2f} (=27.5)  {ve_l_min[name]['high_30']:6.2f} (=30) L/min")
    walking_to_rest_ve_ratio_mid = ve_l_min["combined_corrected"]["mid_27p5"] / ve_l_min["rest_basal"]["mid_27p5"]
    print(f"combined_corrected walking VE / rest VE (mid ratio) = {walking_to_rest_ve_ratio_mid:.2f}x "
          f"(a derived internal-consistency number, not itself an independently-verified external "
          f"citation -- plausible order of magnitude for light-moderate walking exercise hyperpnea)")

    print("\n" + "=" * 78)
    print("STEP 7/8 -- external anchor #3 (VE): Wikipedia 'Minute ventilation' [8] textbook bands "
          "+ discriminating-power / triangulation check (the forced adversary)")
    print("=" * 78)
    rest_ve_mid = ve_l_min["rest_basal"]["mid_27p5"]
    rest_ve_in_wiki_band = WIKI_VE_REST_LOW_L_MIN <= rest_ve_mid <= WIKI_VE_REST_HIGH_L_MIN
    print(f"Rest VE (mid) = {rest_ve_mid:.2f} L/min  vs  Wikipedia resting band ["
          f"{WIKI_VE_REST_LOW_L_MIN},{WIKI_VE_REST_HIGH_L_MIN}] L/min: "
          f"{'PASS' if rest_ve_in_wiki_band else 'FAIL'}")
    for name in ("umberger_primary", "bhargava_primary", "combined_corrected"):
        v = ve_l_min[name]["mid_27p5"]
        band = ("MODERATE-EXERCISE band [40,60]" if WIKI_VE_MODERATE_LOW_L_MIN <= v <= WIKI_VE_MODERATE_HIGH_L_MIN
                else "LIGHT-ACTIVITY-to-moderate transition (~12 to 40)" if WIKI_VE_LIGHT_L_MIN <= v < WIKI_VE_MODERATE_LOW_L_MIN
                else "outside named textbook bands")
        print(f"  {name:20s} VE(mid)={v:6.2f} L/min  -> falls in Wikipedia's {band}")
    print("\nDISCRIMINATING-POWER / FORCED-ADVERSARY check: does this pipeline actually DISCRIMINATE "
          "the reference body's specific number, or would ANY plausible metabolic-rate input pass "
          "(a degenerate, non-informative mapping)? Evidence: the published-VO2-range gate "
          f"(Step 5) PASSES ONLY for combined_corrected ({walking_gate_results['combined_corrected']}) "
          f"and FAILS for both primary configs (umberger={walking_gate_results['umberger_primary']}, "
          f"bhargava={walking_gate_results['bhargava_primary']}) -- both primary configs read "
          f"{vo2_ml_kg_min['umberger_primary']/PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN:.1f}x-"
          f"{vo2_ml_kg_min['bhargava_primary']/PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN:.1f}x the anchor's "
          "own upper bound, and their VE lands in the 'moderate-to-vigorous exercise' Wikipedia "
          "band rather than a walking-appropriate one. This is the SAME discrimination pattern "
          "metabolic_cost.py's Koelewijn-COT anchor and thermoregulation.py's Saltin-"
          "Hermansen-validity-domain anchor already found -- now via a THIRD, differently-sourced "
          "external literature base (MET/Compendium classification + textbook ventilation ranges), "
          "not a re-run of the same comparison. SHARED-INPUT CAVEAT (stated once, applies "
          "throughout, see docstring): all three layers still judge the SAME single upstream "
          "metabolic-rate number -- this is convergent EXTERNAL evidence, not three independent "
          "re-derivations of the reference body's muscle model.")

    print("\n" + "=" * 78)
    print("STEP 8/8 -- void-floor / non-degeneracy sweep + pre-registered gates")
    print("=" * 78)
    m_sweep_w_per_kg = np.linspace(SENSITIVITY_SWEEP_LOW_W, SENSITIVITY_SWEEP_HIGH_W, 25) / mass_kg
    vo2_sweep = vo2_ml_per_kg_min(m_sweep_w_per_kg)
    ve_sweep = ve_l_per_min(vo2_sweep * mass_kg / 1000.0, VE_VO2_MID)
    d_vo2_d_p = np.gradient(vo2_sweep, m_sweep_w_per_kg)
    d_ve_d_p = np.gradient(ve_sweep, m_sweep_w_per_kg)
    analytical_d_vo2_d_p = 60.0 / E_O2_KJ_PER_L
    analytical_d_ve_d_p = analytical_d_vo2_d_p * mass_kg / 1000.0 * VE_VO2_MID
    vo2_deriv_match = bool(np.allclose(d_vo2_d_p, analytical_d_vo2_d_p, rtol=1e-6))
    ve_deriv_match = bool(np.allclose(d_ve_d_p, analytical_d_ve_d_p, rtol=1e-6))
    monotonic = bool(np.all(np.diff(vo2_sweep) > 0) and np.all(np.diff(ve_sweep) > 0))
    print(f"Sweep P in [{SENSITIVITY_SWEEP_LOW_W:.0f},{SENSITIVITY_SWEEP_HIGH_W:.0f}] W total "
          f"({SENSITIVITY_SWEEP_LOW_W/mass_kg:.2f}-{SENSITIVITY_SWEEP_HIGH_W/mass_kg:.2f} W/kg), n={len(m_sweep_w_per_kg)}: "
          f"VO2 range [{vo2_sweep.min():.2f},{vo2_sweep.max():.2f}] mL/kg/min, "
          f"VE(mid) range [{ve_sweep.min():.2f},{ve_sweep.max():.2f}] L/min")
    print(f"Strictly monotonic (real, non-pinned function of the reference body's number): {'PASS' if monotonic else 'FAIL'}")
    print(f"Analytical-vs-numerical derivative match: VO2 d/dP={analytical_d_vo2_d_p:.4f} mL/kg/min per W/kg "
          f"({'PASS' if vo2_deriv_match else 'FAIL'}); VE d/dP={analytical_d_ve_d_p:.4f} L/min per W/kg "
          f"({'PASS' if ve_deriv_match else 'FAIL'})")

    all_configs_positive_and_finite = bool(np.all(np.isfinite(list(vo2_ml_kg_min.values()))) and
                                            np.all(np.array(list(vo2_ml_kg_min.values())) > 0))
    gates = {
        "basal_self_consistent_across_models": basal_consistent,
        "e_o2_implied_rq_plausible": rq_plausible,
        "rest_vo2_within_10pct_of_1MET_definition": rest_vo2_1met_ok,
        "combined_corrected_walking_vo2_within_published_range": walking_gate_results["combined_corrected"],
        "primary_configs_exceed_published_range_expected_direction": (
            not walking_gate_results["umberger_primary"] and not walking_gate_results["bhargava_primary"]
            and vo2_ml_kg_min["umberger_primary"] > PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN
            and vo2_ml_kg_min["bhargava_primary"] > PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN
        ),
        "rest_ve_within_wikipedia_band": rest_ve_in_wiki_band,
        "walking_ve_exceeds_rest_ve": bool(ve_l_min["combined_corrected"]["mid_27p5"] > 2 * ve_l_min["rest_basal"]["mid_27p5"]),
        "sensitivity_sweep_monotonic": monotonic,
        "sensitivity_sweep_derivative_match": bool(vo2_deriv_match and ve_deriv_match),
        "all_vo2_finite_and_positive": all_configs_positive_and_finite,
    }
    gates = {k: bool(v) for k, v in gates.items()}
    overall_pass = all(gates.values())
    print(f"\nGATES: {json.dumps(gates, indent=2)}")
    print(f"OVERALL: {'PASS' if overall_pass else 'FAIL -- see gates above'}")

    report.update({
        "citations_verified_live": {
            "ainsworth_compendium_2011_pmid": "21681120", "ainsworth_compendium_2011_doi": "10.1249/MSS.0b013e31821ece12",
            "waters_mulroy_1999_pmid": "10575082", "waters_mulroy_1999_doi": "10.1016/s0966-6362(99)00009-0",
            "browning_kram_2005_pmid": "15919843", "browning_kram_2005_doi": "10.1038/oby.2005.103",
            "browning_baker_herron_kram_2006_pmid": "16210434", "browning_baker_herron_kram_2006_doi": "10.1152/japplphysiol.00767.2005",
            "wasserman_whipp_koyl_beaver_1973_pmid": "4723033", "wasserman_whipp_koyl_beaver_1973_doi": "10.1152/jappl.1973.35.2.236",
            "ats_accp_cpet_statement_2003_pmid": "12524257", "ats_accp_cpet_statement_2003_doi": "10.1164/rccm.167.2.211",
            "peronnet_massicotte_1991_pmid": "1645211", "peronnet_massicotte_1991_doi": None,
        },
        "inputs": {"mass_kg": mass_kg, "configs_w_per_kg": configs_w_per_kg,
                   "basal_self_consistent": basal_consistent},
        "e_o2_kj_per_l": E_O2_KJ_PER_L, "e_o2_implied_rq": implied_rq,
        "vo2_ml_per_kg_min": vo2_ml_kg_min, "vo2_l_per_min": vo2_l_min,
        "compendium_walking_vo2_ml_per_kg_min": compendium_vo2,
        "published_walking_vo2_range_ml_per_kg_min": [PUBLISHED_WALKING_VO2_LOW_ML_KG_MIN, PUBLISHED_WALKING_VO2_HIGH_ML_KG_MIN],
        "walking_vo2_in_published_range": walking_gate_results,
        "ve_l_per_min": ve_l_min,
        "walking_to_rest_ve_ratio_mid": walking_to_rest_ve_ratio_mid,
        "wikipedia_ve_bands_l_per_min": {
            "rest": [WIKI_VE_REST_LOW_L_MIN, WIKI_VE_REST_HIGH_L_MIN], "light": WIKI_VE_LIGHT_L_MIN,
            "moderate": [WIKI_VE_MODERATE_LOW_L_MIN, WIKI_VE_MODERATE_HIGH_L_MIN],
        },
        "sensitivity_sweep": {
            "power_w_per_kg": m_sweep_w_per_kg.tolist(), "vo2_ml_per_kg_min": vo2_sweep.tolist(),
            "ve_l_per_min_mid": ve_sweep.tolist(),
            "analytical_d_vo2_d_p": analytical_d_vo2_d_p, "analytical_d_ve_d_p": analytical_d_ve_d_p,
        },
        "gates": gates,
        "overall_pass": overall_pass,
        "reference_body": {
            "inherited_from": "metabolic_cost_results.json (mass_kg = mc['muscle_mass']"
                              "['total_body_mass_kg'], the reference body's scaled musculoskeletal-model mass, "
                              "read-only, not re-derived here)",
            "mass_kg": mass_kg, "name": None, "body_fat_fraction": None,
            "class": "reference_body",
        },
    })
    if verdict.state == "ACCEPT_SYNTHETIC_DEMO":
        report["evidence_class"] = "synthetic_demo"
        report["scientific_claim"] = False
    else:
        report["evidence_class"] = verdict.evidence_class
        report["scientific_claim"] = True
    report["preflight"] = {
        "mode": mode,
        "state": verdict.state,
        "failures": list(verdict.failures),
        "flags": list(verdict.flags),
        "evidence_class": verdict.evidence_class,
    }
    out_path = f"{OUT_DIR}/respiratory_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
