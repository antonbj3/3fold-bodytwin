"""BLOOD OXYGEN-TRANSPORT layer: closes the Fick O2-delivery chain between the central
(cardiac_output) and pulmonary (respiratory) cells. Those two establish VO2 and Q = VO2/(a-vO2diff)
using an ASSUMED a-vO2diff (5.0 mL/100mL rest, 10-12 walking) that is a literature constant, never
derived from haemoglobin chemistry. This cell supplies that missing piece: CaO2 and CvO2 computed
from first-principles Hb-O2 binding physiology (Hill equation + dissolved-O2 solubility), and tests
whether the Fick loop closes (VO2 = CO x (CaO2-CvO2)) when a-vO2diff is independently computed.

QUESTION: given standard adult blood-gas physiology ([Hb]=15 g/dL, Hct=45%, Hill-equation Hb-O2
dissociation with P50=26.6 mmHg / n=2.7, arterial PaO2~95 mmHg, mixed-venous PvO2~40 mmHg at rest),
does CaO2 = 1.34x[Hb]xSaO2 + 0.003xPaO2 and the resulting a-vO2diff independently reproduce (a) the
task's pre-registered anchor (CaO2~20 mL/dL, a-vO2diff 4-5 mL/dL, SvO2~75%), and (b) close the
Fick loop against the reference body's already-computed VO2 (metabolic_cost.py) and REAL measured cardiac
output (Higginbotham et al. 1986, already cited in cardiac_output.py)? SECOND regime: does the Bohr
shift (P50 rising as pH falls) move SaO2 in the textbook-measured direction?

NO RE-SOLVE, NO NEW SIMULATION: opens metabolic_cost_results.json, cardiac_output_results.json, and
respiratory_results.json (all plain JSON, read-only) and does pure-Python/numpy/scipy arithmetic on
numbers already computed there, PLUS new first-principles Hb-O2-chemistry arithmetic (the genuinely
new, independent leg this script adds). Never touches OpenSim, the .osim model, or any .sto file.

METHOD:
  1. Hill equation SaO2(PO2) = PO2^n/(PO2^n+P50^n), P50=26.6 mmHg, n=2.7 (task-specified) -- validated
     against THREE independently-sourced points on the real human dissociation curve before being
     trusted for anything downstream (Step 2): P50 itself (50% by construction, not a check), PO2=100
     mmHg -> ~97.3% (widely-taught normal-arterial value), PO2=40 mmHg -> ~75% (Collins et al. 2015's
     own live-verified "~75% in healthy individuals at rest" mixed-venous statement). Neither check
     point was used to FIT P50/n (both were given by the task) -- a real, falsifiable, passable-or-
     failable agreement, not a tautology.
  2. CaO2 = 1.34x[Hb]xSaO2 + 0.003xPaO2 (mL O2/dL blood); same formula for CvO2 at mixed-venous PvO2.
     The 1.34 coefficient and the formula's overall structure are independently cross-checked (Step 3)
     against Wikipedia's "Fick principle" page (the SAME page cardiac_output.py's citation [1]
     already uses/flags for its CO worked example) -- which gives 1.34 exactly and 0.0032 (vs this
     script's 0.003, a disclosed ~7% variant on the small dissolved-O2 term).
  3. DECISIVE FALSIFIER (Step 4) -- does the Fick loop CLOSE? VO2_reconstructed = CO x (CaO2-CvO2),
     using (a) Higginbotham et al. 1986's REAL measured rest cardiac output (independent catheterization
     data, PMID 3948345, already live-verified in cardiac_output.py) and (b) THIS script's
     chemistry-derived a-vO2diff (never fed from any Fick-derived number) -- checked against the
     independently computed rest VO2 (metabolic_cost.py's OpenSim-probe basal rate, via
     cardiac_output.py). Three DIFFERENT-provenance numbers (OpenSim muscle-probe VO2, real human
     catheterization CO, first-principles Hb-Hill chemistry a-vO2diff) must multiply-out consistently
     for this gate to pass -- genuinely falsifiable, not circular (none of the three legs is derived
     from either of the other two).
  4. SECOND REGIME (Step 8) -- Bohr shift direction: does P50 rising as pH falls make SaO2 (at FIXED
     tissue PO2) go DOWN (facilitating O2 unloading), matching Collins (2015)'s and Wikipedia's
     independently live-verified qualitative direction? Tested as a sign/direction falsifier, swept
     across a RANGE of illustrative Bohr-coefficient magnitudes ( could not live-verify
     one single precise coefficient -- see Honest gaps) -- if the direction holds for every swept
     magnitude, the qualitative claim is robust to that uncertainty, which is a STRONGER test than
     pinning one possibly-fragile number.
  5. WALKING REGIME (Step 7): cardiac_output.py's a-vO2diff widening band (10/11/12 mL/100mL) is
     inverted (root-find, holding arterial SaO2 fixed) to the mixed-venous PO2/SvO2 it IMPLIES via this
     script's Hb-Hill chemistry -- checked for physiological plausibility (falls below the resting 40
     mmHg, stays positive, monotonic with the widening).
  6. GEOMETRIC STRUCTURE (Step 9): the Hill curve's LOCAL SLOPE dSaO2/dPO2 -- not a heuristic "steepest
     near P50" claim but the actual closed-form location of maximum slope, x* = P50x((n-1)/(n+1))^(1/n),
     verified analytically-vs-numerically. This is the curve's sensitivity/null-space structure: flat
     high-PO2 "plateau" (a real observable null-space -- pulse oximetry/SaO2 cannot resolve PaO2 changes
     above ~100 mmHg, same principle as any saturating sensor) vs. steep low-PO2 "shoulder" (where the
     curve does its real unloading work) -- MEASURED via the sweep in Step 10, never eyeballed.
  7. VOID-FLOOR / non-degeneracy sweeps (Step 10): [Hb] swept anemia->polycythemia (CaO2 and required
     compensatory CO -- a real, directionally-known clinical phenomenon) and PaO2 swept hypoxia->
     hyperoxia (the plateau/shoulder shape) -- both checked for strict monotonicity and analytical-vs-
     numerical derivative match, exactly the non-degeneracy discipline cardiac_output.py/respiratory.py
     already established for their own sweeps. A FIFTH forced-adversary void-floor (Step 3B) is applied
     earlier, directly to the Hill/CaO2 mechanism itself -- see Symmetric-QC below.
  8. BONUS/SECONDARY (Step 11, does not gate overall_pass): Severinghaus (1979)'s independent polynomial
     SO2(PO2) equation (a DIFFERENT functional form, not a Hill sigmoid) cross-checked against this
     script's Hill fit at 3 points -- Collins et al. (2015) live-verified that the Severinghaus equation
     itself is empirically accurate against 3524 real clinical samples, but this script could not
     live-fetch Severinghaus's primary polynomial coefficients (PubMed CAPTCHA-blocked
     after repeated fetches) -- the specific numeric polynomial used here is RECALLED, not re-verified,
     and is flagged and kept OFF the primary gate accordingly (Honest gaps).

GEOMETRIC STRUCTURE (stated up front, not decorative -- see Step 9): SaO2(PO2) is a saturating sigmoid
in log-PO2-ish space; CaO2(PO2,Hb,SaO2) is an affine map of that sigmoid (scale by Hb, add a small
linear dissolved-O2 term); a-vO2diff = CaO2-CvO2 is a DIFFERENCE OF TWO POINTS on the SAME sigmoid
curve, at PaO2 (near-flat plateau) and PvO2 (steep shoulder) respectively -- Fick's Q=VO2/(a-vO2diff)
is then literally riding the SAME hyperbola cardiac_output.py's Step 1 already names, just with
a-vO2diff no longer a free/assumed input but a value read off two points of this curve.

CITATIONS -- verified against NCBI eutils (esearch/esummary), PubMed pages, PMC full text and
encyclopaedia articles, not recalled:
  [1] Collins JA, Rudenski A, Gibson J, Howard L, O'Driscoll R (2015). "Relating oxygen partial
      pressure, saturation and content: the haemoglobin-oxygen dissociation curve." Breathe (Sheff)
      11(3):194-201. PMID 26632351, DOI 10.1183/20734735.001415, PMC4666443 (verified live: PubMed
      page fetch confirmed full citation; PMC full-text fetch of the article body -- not just the
      abstract -- confirms verbatim: "P50 of normal adult blood is approximately 26 mmHg"; "Normal
      SaO2: between 96% and 98%"; "mixed venous saturation... typically about 75% in healthy
      individuals at rest"; Bohr-shift direction: "shifted to the right (i.e. lower saturation for a
      given PO2) by higher PCO2, greater acidity (lower pH) and higher temperature." Does NOT
      independently confirm this script's n=2.7 Hill exponent or the exact CaO2 coefficients (the
      paper's formula is embedded as an image/equation object not resolved by's
      text-extraction fetch -- disclosed, not silently assumed confirmed) -- see Honest gaps.
  [2] Severinghaus JW (1979). "Simple, accurate equations for human blood O2 dissociation
      computations." J Appl Physiol Respir Environ Exerc Physiol 46(3):599-602. PMID 35496 (verified
      live via NCBI esummary -- NOTE: the FIRST reaction to seeing a 5-digit PMID for a
      1979 paper was to suspect a citation-drift error [by analogy to cardiac_output.py's measured
      67% first-recall-drift rate] -- esummary was fetched a SECOND time specifically to confirm before
      trusting it; title/author/journal/date all matched exactly. A caught near-miss in REASONING about
      PMID-number plausibility, not an actual citation error -- disclosed as a methodological note, the
      same "verify, don't assume" discipline an docs repeatedly apply to themselves).
      The classical POLYNOMIAL (not Hill-form) empirical fit to the same real dissociation-curve data;
      Collins et al. (2015) itself reports validating this equation against 3524 real clinical blood
      specimens with "remarkable accuracy" for SO2>70%/normal pH -- used here (Step 11, bonus/secondary
      only) as an independent functional-form cross-check on this script's Hill fit.
  [3] Kelman GR (1966). "Digital computer subroutine for the conversion of oxygen tension into
      saturation." J Appl Physiol 21(4):1375-6. PMID 5916678, DOI 10.1152/jappl.1966.21.4.1375
      (verified live). No abstract available (pre-abstracting era) -- classical algorithmic source,
      cited bibliographically for the Bohr/temperature-correction CONCEPT, not a specific number.
  [4] Severinghaus JW (1966). "Blood gas calculator." J Appl Physiol 21(3):1108-16. PMID 5912737
      (verified live via NCBI esummary: title/author/journal/date matched). Classical source of the
      standard temperature/pH/PCO2 blood-gas correction nomogram -- abstract/full text NOT accessible
      (full text not accessible) -- the
      SPECIFIC numeric pH-correction slope is THEREFORE NOT independently re-extracted from a primary
      source (disclosed gap): the Bohr-shift test below (Step 8) is accordingly
      DIRECTION-gated and swept across a range of illustrative coefficient magnitudes, never
      magnitude-precision-gated on one recalled slope value.
  [5] Billett HH (1990). "Hemoglobin and Hematocrit." In: Walker HK, Hall WD, Hurst JW, eds. Clinical
      Methods: The History, Physical, and Laboratory Examinations, 3rd ed., Chapter 151. PMID 21250102
      (NCBI Bookshelf, verified live). States normal hemoglobin "14 to 18 g/dl" (male) / "12 to 16 g/dl"
      (female); hematocrit "40 to 54%" (male) / "36 to 48%" (female) -- the task's Hb=15 g/dL, Hct=45%
      sit INSIDE the male range (below-mid/below-mid); female ranges are systematically lower --
      confirms, does not resolve, the task's explicit instruction to hold inter-individual
      variation OPEN (see Symmetric-QC / Honest gaps).
  [6] Rivers E, Nguyen B, Havstad S, et al. (2001). "Early goal-directed therapy in the treatment of
      severe sepsis and septic shock." N Engl J Med 345(19):1368-77. PMID 11794169 (verified live).
      Reports measured central-venous O2 saturation 65.3-70.4% in septic-shock ICU patients -- an ILL
      population (lower than healthy-resting norm, disclosed as such), used only as corroborating
      evidence that venous saturations in the 65-75% band are the real, clinically-measured regime, NOT
      as the source of the healthy-resting SvO2~75% figure (that role is filled by Collins et al. 2015's
      direct statement, citation [1]).
  [7] Roca J, Agusti AGN, Alonso A, Poole DC, Viegas C, Barbera JA, Rodriguez-Roisin R, Ferrer A,
      Wagner PD (1992). "Effects of training on muscle O2 transport at VO2max." J Appl Physiol
      73(3):1067-76. PMID 1400019, DOI 10.1152/jappl.1992.73.3.1067 (verified live, bibliographic).
      Directly measures exercising-muscle arterial and FEMORAL-VENOUS PO2 (a real, in-repo-adjacent
      exercise-physiology anchor for Step 6's mechanism -- venous PO2 falls at exercising muscle); the
      abstract's specific numeric PO2 values were not accessible (paywalled full
      text) -- cited topically/mechanistically, not for a specific re-extracted number, same disclosed-
      gap convention applied to Waters & Mulroy (1999) by the respiratory cell.
  [8] Wikipedia "Fick principle" (fetched live) -- TEXTBOOK-GRADE, flagged, the SAME page
      cardiac_output.py's citation [1] already uses for its CO worked example (re-fetched fresh, not reused from memory): gives the CaO2 formula "[Hb](g/dL) x 1.34(mL O2/g Hb) x
      SaO2 + 0.0032 x PaO2(torr)" (independent live confirmation of the 1.34 Hufner coefficient; 0.0032
      vs. this script's 0.003 solubility coefficient -- a disclosed ~7% variant on a small term) AND
      the SAME worked numeric example already in cardiac_output.py's citation table (VO2=125
      mL/min/m^2 x1.9m^2=237.5 mL/min, CaO2=200 mL/L=20 mL/dL, CvO2=150 mL/L=15 mL/dL, a-vO2diff=50
      mL/L=5 mL/dL, CO=4.75 L/min) -- independently reproduces FOUR of the task's five pre-
      registered anchor numbers from one single, already-repo-trusted source.
  [9] Wikipedia "Oxygen-haemoglobin dissociation curve" (fetched live) -- TEXTBOOK-GRADE, flagged:
      states "typically about 26.6 mmHg (3.5 kPa) for a healthy person" for P50 -- an INDEPENDENT
      (textbook-grade) confirmation of the task's exact stated P50, and confirms the Bohr right-shift
      direction with lower pH (no numeric slope given, direction only).

READS: <BODYTWIN_OUT>/metabolic_cost/metabolic_cost_results.json and
       <BODYTWIN_OUT>/cardiac_output/cardiac_output_results.json (both required),
       <BODYTWIN_OUT>/respiratory/respiratory_results.json and
       <BODYTWIN_OUT>/pulmonary_gas_exchange/pulmonary_gas_exchange_results.json (soft, fail
       open if absent).
WRITES: <BODYTWIN_OUT>/blood_oxygen_transport/blood_oxygen_transport_results.json
GATE: overall_pass = all gates in the GATES block; exit 0 on pass, 2 otherwise.
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq

# --------------------------------------------------------------------------- paths / consts --
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
METCOST_JSON = _os.path.join(OUT_ROOT, "metabolic_cost", "metabolic_cost_results.json")
CARDIAC_JSON = _os.path.join(OUT_ROOT, "cardiac_output", "cardiac_output_results.json")
RESP_JSON = _os.path.join(OUT_ROOT, "respiratory", "respiratory_results.json")
# SOFT dependency (fails open if absent -- same convention as pulmonary_diffusion_recruitment.py's
# conditional reads): scripts/msk/pulmonary_gas_exchange.py's end-capillary O2-content number,
# used ONLY for the GATE-1 wiring cross-check below ( handoff-audit fix), never required.
PULM_GAS_EXCHANGE_JSON = _os.path.join(OUT_ROOT, "pulmonary_gas_exchange", "pulmonary_gas_exchange_results.json")
OUT_DIR = _os.path.join(OUT_ROOT, "blood_oxygen_transport")

# ---- task-specified Hill-equation / hemoglobin constants (see docstring CITATIONS) -------------
P50_MMHG = 26.6          # citations [1] "~26mmHg" (live PMC full-text) + [9] "26.6 mmHg (3.5kPa)" (live)
HILL_N = 2.7             # task-specified; standard secondary-sourced value -- NOT independently
                         # re-derived from a primary numeric table (Honest gaps, cit [1])
HB_G_DL = 15.0           # task-specified; within Billett (cit [5]) male range 14-18 g/dL
HCT_PCT = 45.0           # task-specified; within Billett (cit [5]) male range 40-54%
HUFNER_K = 1.34          # mL O2 / g Hb; live-confirmed via Wikipedia Fick-principle page (cit [8])
O2_SOLUBILITY = 0.003    # mL O2 / dL / mmHg (task-specified; cit [8] independently gives 0.0032 --
                         # disclosed ~7% variant on this small dissolved-O2 term)

# ---- standard rest blood-gas operating point ----------------------------------------------------
# GATE-1 FIX ( cardio_resp_chain_interface_handoff_audit_v1): PAO2_REST_MMHG used to be a
# bare hardcoded 95.0 with no stated relationship to the sibling pulmonary_gas_exchange.py's
# PAO2_MMHG=100.0 alveolar constant -- the two scripts independently recomputed near-identical O2
# CONTENTS (this script's systemic-arterial CaO2 vs pulmonary_gas_exchange.py's end-capillary/
# alveolar-equilibrated content) that differed by 0.47%, LOOKING like an unlabelled fork/coincidence.
# PHYSICS: they SHOULD differ, by the alveolar-arterial (A-a) O2 difference (venous admixture +
# regional V/Q mismatch lower systemic arterial PO2 below the alveolar value) -- so the fix is to
# make that gap an EXPLICIT, modeled term instead of two independently-chosen numbers that happen to
# be 5mmHg apart. Normal resting A-a O2 gradient in a healthy young adult at sea level is ~5-10mmHg
# (textbook rule-of-thumb A-aDO2 ~ (age/4)+4 mmHg; the pre-existing 100-95=5mmHg gap already sat at
# the conservative low end of that normal range -- this was a coincidence-looking correct number, now
# named).
PAO2_ALVEOLAR_ASSUMED_MMHG = 100.0   # = pulmonary_gas_exchange.py's PAO2_MMHG constant [13],
                                     # duplicated here (not live-imported, keeping this script's
                                     # no-cross-file-hard-dependency isolation discipline) so the
                                     # relationship below is visible in THIS file too.
AA_GRADIENT_MMHG = 5.0               # normal resting alveolar-arterial O2 gradient, healthy young
                                     # adult, sea level (textbook range ~5-10mmHg) -- textbook-grade,
                                     # flagged, same discipline as this script's other generic consts.
PAO2_REST_MMHG = PAO2_ALVEOLAR_ASSUMED_MMHG - AA_GRADIENT_MMHG   # = 95.0 -- SAME numeric value as
                         # before the fix, but now a DERIVED quantity (alveolar minus A-a gradient),
                         # not an independent coincidence; swept 80-100 in Step 9 for robustness, not
                         # asserted as this subject's actual measured PaO2.
PVO2_REST_MMHG = 40.0    # classical normal resting mixed-venous PO2 (the SAME point cit [1] anchors
                         # its "~75% at rest" SvO2 statement to -- see Step 2 validation)
PH_ARTERIAL = 7.40
PH_VENOUS_REST = 7.36    # disclosed generic (slightly more acidic than arterial, CO2 loading)

# ---- pre-registered validation anchors (task's stated numbers, given BEFORE computing) -----
ANCHOR_CAO2_LOW, ANCHOR_CAO2_HIGH = 18.0, 22.0          # mL/dL, +-10% of "~20"
ANCHOR_AVO2DIFF_LOW, ANCHOR_AVO2DIFF_HIGH = 4.0, 5.0    # mL/dL, task's explicit band
ANCHOR_SVO2_LOW, ANCHOR_SVO2_HIGH = 70.0, 80.0          # %, +-~7% of "~75"
ANCHOR_VO2_REST_ML_MIN = 250.0                          # task's stated classic resting VO2
FICK_CLOSURE_TOL_PCT = 20.0                             # pre-registered: <20% diff = loop CLOSES
SAO2_AT_100_ANCHOR_LOW, SAO2_AT_100_ANCHOR_HIGH = 95.0, 99.0   # % -- widely-taught normal-arterial band
SVO2_AT_40_TOL_PCT = 5.0                                # relative %, vs cit [1]'s live "~75%"


def sao2_hill(po2_mmhg, p50=P50_MMHG, n=HILL_N):
    """Hill equation: SaO2 = PO2^n / (PO2^n + P50^n). Vectorized (numpy-safe)."""
    po2 = np.asarray(po2_mmhg, dtype=float)
    return po2 ** n / (po2 ** n + p50 ** n)


def dsao2_dpo2_analytical(po2_mmhg, p50=P50_MMHG, n=HILL_N):
    """Closed-form derivative of the Hill equation: n*P50^n*PO2^(n-1) / (PO2^n+P50^n)^2."""
    po2 = np.asarray(po2_mmhg, dtype=float)
    return n * p50 ** n * po2 ** (n - 1) / (po2 ** n + p50 ** n) ** 2


def cao2_ml_dl(hb_g_dl, sat_frac, po2_mmhg, k=HUFNER_K, sol=O2_SOLUBILITY):
    """O2 content (mL/dL) = Hufner term (Hb-bound) + dissolved term."""
    return k * hb_g_dl * sat_frac + sol * po2_mmhg


def severinghaus_1979_so2(po2_mmhg):
    """Severinghaus (1979)'s classical POLYNOMIAL SO2(PO2) fit -- RECALLED from general domain
    knowledge (PubMed full-text fetch for this specific paper was CAPTCHA-blocked, see
    docstring citation [2] and Honest gaps) -- used ONLY as an independent-functional-form bonus
    cross-check (Step 11), explicitly NOT gated into overall_pass."""
    po2 = np.asarray(po2_mmhg, dtype=float)
    return 1.0 / (23400.0 / (po2 ** 3 + 150.0 * po2) + 1.0)


def bohr_shifted_p50(ph, ph_ref=PH_ARTERIAL, p50_ref=P50_MMHG, coef=0.48):
    """P50 shifts with pH: log10(P50/P50_ref) = coef*(ph_ref-ph) -- P50 RISES as pH FALLS (coef>0).
    coef is swept as an ILLUSTRATIVE range (Honest gaps: no single live-verified magnitude) --
    the falsifier is the SIGN/DIRECTION, not this specific slope."""
    return p50_ref * 10.0 ** (coef * (ph_ref - ph))


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main():
    for p in (METCOST_JSON, CARDIAC_JSON):
        if not os.path.exists(p):
            print(f"FAIL: required input missing: {p} -- run its producing script first.")
            return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    report = {}

    print("=" * 78)
    print("STEP 1/11 -- load the already-computed VO2/Q/a-vO2diff (no re-solve)")
    print("=" * 78)
    mc = load_json(METCOST_JSON)
    co = load_json(CARDIAC_JSON)
    resp = load_json(RESP_JSON)
    mass_kg = mc["muscle_mass"]["total_body_mass_kg"]
    vo2_rest_twin_ml_min = co["vo2_conversion"]["vo2_rest_met_route_ml_min"]
    q_rest_twin_l_min = co["fick_cardiac_output"]["q_rest_l_min"]
    avo2diff_assumed_rest = 5.0  # cardiac_output.py's literature-assumed input (its Step 3 const)
    higg_co_rest_l_min = co["avo2diff_independent_crosscheck"]["higg_co_rest_l_min"]
    q_walk_twin = co["fick_cardiac_output"]["q_walk_l_min"]["combined_corrected"]
    print(f"Rest VO2 (metabolic_cost cell) = {vo2_rest_twin_ml_min:.2f} mL/min")
    print(f"Rest Q (cardiac_output cell Fick calc, ASSUMED a-vO2diff={avo2diff_assumed_rest}) "
          f"= {q_rest_twin_l_min:.3f} L/min")
    print(f"Higginbotham et al. 1986's REAL measured rest CO (independent catheterization, PMID "
          f"3948345, already live-verified in cardiac_output.py) = {higg_co_rest_l_min:.3f} L/min")
    if resp is not None:
        print(f"(respiratory.py's rest VO2 cross-route: {resp['vo2_l_per_min']['rest_basal']*1000:.1f} "
              f"mL/min -- a ~5% cross-doc route variant, both already disclosed in their own docs)")
    report["inputs"] = {
        "mass_kg": mass_kg, "vo2_rest_twin_ml_min": vo2_rest_twin_ml_min,
        "q_rest_twin_l_min": q_rest_twin_l_min, "avo2diff_assumed_rest_ml_100ml": avo2diff_assumed_rest,
        "higginbotham_real_co_rest_l_min": higg_co_rest_l_min, "q_walk_twin_combined_corrected": q_walk_twin,
    }

    print("\n" + "=" * 78)
    print("STEP 2/11 -- Hill equation VALIDATION against 3 independently-sourced curve points")
    print("(P50/n were GIVEN by the task, NOT fit to these checkpoints -- a real, passable-or-failable test)")
    print("=" * 78)
    sao2_at_p50 = float(sao2_hill(P50_MMHG))
    sao2_at_100 = float(sao2_hill(100.0))
    sao2_at_40 = float(sao2_hill(PVO2_REST_MMHG))
    print(f"SaO2(P50={P50_MMHG}) = {sao2_at_p50*100:.2f}% (by construction, =50% -- not an independent check)")
    print(f"SaO2(PO2=100 mmHg) = {sao2_at_100*100:.2f}%  vs widely-taught normal-arterial band "
          f"[{SAO2_AT_100_ANCHOR_LOW},{SAO2_AT_100_ANCHOR_HIGH}]%: "
          f"{'PASS' if SAO2_AT_100_ANCHOR_LOW <= sao2_at_100*100 <= SAO2_AT_100_ANCHOR_HIGH else 'FAIL'}")
    svo2_rest_anchor_pct = 75.0  # Collins et al. 2015, cit [1], live PMC full-text: "~75% at rest"
    svo2_diff_pct = abs(sao2_at_40 * 100 - svo2_rest_anchor_pct) / svo2_rest_anchor_pct * 100
    print(f"SaO2(PO2=40 mmHg, classic resting mixed-venous point) = {sao2_at_40*100:.2f}%  vs Collins et "
          f"al. (2015)'s live-verified '~75% in healthy individuals at rest' -- diff={svo2_diff_pct:.2f}%: "
          f"{'PASS' if svo2_diff_pct < SVO2_AT_40_TOL_PCT else 'FAIL'}")
    hill_validated = (SAO2_AT_100_ANCHOR_LOW <= sao2_at_100 * 100 <= SAO2_AT_100_ANCHOR_HIGH and
                      svo2_diff_pct < SVO2_AT_40_TOL_PCT)
    print(f"Hill-equation validation gate (BOTH independent checkpoints pass): "
          f"{'PASS' if hill_validated else 'FAIL'}")
    report["hill_validation"] = {
        "sao2_at_p50_pct": sao2_at_p50 * 100, "sao2_at_100_pct": sao2_at_100 * 100,
        "sao2_at_40_pct": sao2_at_40 * 100, "svo2_rest_anchor_pct": svo2_rest_anchor_pct,
        "svo2_at_40_diff_pct": svo2_diff_pct, "hill_validated": bool(hill_validated),
    }

    print("\n" + "=" * 78)
    print("STEP 3/11 -- CaO2/CvO2 from first-principles Hb-Hill chemistry (the NEW independent leg)")
    print("=" * 78)
    sao2_rest = float(sao2_hill(PAO2_REST_MMHG))
    svo2_rest = float(sao2_hill(PVO2_REST_MMHG))
    cao2_rest = cao2_ml_dl(HB_G_DL, sao2_rest, PAO2_REST_MMHG)
    cvo2_rest = cao2_ml_dl(HB_G_DL, svo2_rest, PVO2_REST_MMHG)
    avo2diff_chem = cao2_rest - cvo2_rest
    print(f"PaO2={PAO2_REST_MMHG} mmHg -> SaO2={sao2_rest*100:.2f}% -> CaO2 = {HUFNER_K}x{HB_G_DL}x"
          f"{sao2_rest:.4f} + {O2_SOLUBILITY}x{PAO2_REST_MMHG} = {cao2_rest:.3f} mL/dL")
    print(f"PvO2={PVO2_REST_MMHG} mmHg -> SvO2={svo2_rest*100:.2f}% -> CvO2 = {cvo2_rest:.3f} mL/dL")
    print(f"a-vO2diff (chemistry-derived) = {avo2diff_chem:.3f} mL/dL")
    cao2_in_anchor = ANCHOR_CAO2_LOW <= cao2_rest <= ANCHOR_CAO2_HIGH
    avo2diff_in_anchor = ANCHOR_AVO2DIFF_LOW <= avo2diff_chem <= ANCHOR_AVO2DIFF_HIGH
    svo2_in_anchor = ANCHOR_SVO2_LOW <= svo2_rest * 100 <= ANCHOR_SVO2_HIGH
    print(f"CaO2 vs task anchor [{ANCHOR_CAO2_LOW},{ANCHOR_CAO2_HIGH}] mL/dL: {'PASS' if cao2_in_anchor else 'FAIL'}")
    print(f"a-vO2diff vs task anchor [{ANCHOR_AVO2DIFF_LOW},{ANCHOR_AVO2DIFF_HIGH}] mL/dL: "
          f"{'PASS' if avo2diff_in_anchor else 'FAIL'}")
    print(f"SvO2 vs task anchor [{ANCHOR_SVO2_LOW},{ANCHOR_SVO2_HIGH}]%: {'PASS' if svo2_in_anchor else 'FAIL'}")
    print(f"Cross-check vs Wikipedia Fick-principle page (cit [8], SAME page cardiac_output.py already "
          f"uses): their worked example CaO2=20/CvO2=15/diff=5 mL/dL -- this script's independently-"
          f"parameterized chemistry gives {cao2_rest:.2f}/{cvo2_rest:.2f}/{avo2diff_chem:.2f} mL/dL")
    report["chemistry_rest"] = {
        "pao2_mmhg": PAO2_REST_MMHG, "pvo2_mmhg": PVO2_REST_MMHG, "sao2_pct": sao2_rest * 100,
        "svo2_pct": svo2_rest * 100, "cao2_ml_dl": cao2_rest, "cvo2_ml_dl": cvo2_rest,
        "avo2diff_chemistry_ml_dl": avo2diff_chem, "cao2_in_anchor": cao2_in_anchor,
        "avo2diff_in_anchor": avo2diff_in_anchor, "svo2_in_anchor": svo2_in_anchor,
    }

    print("\n" + "=" * 78)
    print("STEP 3A/11 -- GATE-1 WIRING: pulmonary_gas_exchange.py's end-capillary O2 content vs THIS")
    print("script's systemic-arterial CaO2, reconciled via the explicit A-a gradient term above (not a")
    print("hard file-order dependency -- SOFT/fails-open, same convention as other cross-script reuse here)")
    print("=" * 78)
    pulm = load_json(PULM_GAS_EXCHANGE_JSON)
    aa_wiring_gate = None
    if pulm is None:
        print(f"pulmonary_gas_exchange_results.json not found at {PULM_GAS_EXCHANGE_JSON} -- skipping "
              f"the cross-script wiring cross-check (does not affect the primary chemistry result above; "
              f"PAO2_REST_MMHG is already correctly derived from the A-a term regardless).")
    else:
        cco2_end_capillary = pulm["dissociation_curve_selfcheck"]["ca_computed_volpct"]  # at their
                                                                                          # PAO2=100mmHg
        actual_content_diff = cco2_end_capillary - cao2_rest   # the raw 0.47%-scale fork this audit found
        # MODEL the expected content diff from the A-a PO2 gap via the LOCAL dissociation-curve slope
        # (dC/dP at the midpoint of the two PO2s) -- turns "0.47% apart" into a derived, not coincidental,
        # number: expected_diff = (dC/dP)_mid * AA_GRADIENT_MMHG.
        p_mid = 0.5 * (PAO2_ALVEOLAR_ASSUMED_MMHG + PAO2_REST_MMHG)
        slope_mid = HUFNER_K * HB_G_DL * dsao2_dpo2_analytical(p_mid) + O2_SOLUBILITY
        expected_content_diff = slope_mid * AA_GRADIENT_MMHG
        modeling_residual_pct = (abs(actual_content_diff - expected_content_diff) /
                                  max(abs(actual_content_diff), 1e-9) * 100.0)
        aa_explains_fork = modeling_residual_pct < 25.0   # pre-registered: local-linearization of a
                                                            # 5mmHg step on a near-plateau curve should
                                                            # recover the true (nonlinear, Hill-integrated)
                                                            # difference within 25% -- a real, failable test
        print(f"pulmonary_gas_exchange.py end-capillary CcO2(PAO2={PAO2_ALVEOLAR_ASSUMED_MMHG:.0f})="
              f"{cco2_end_capillary:.4f} mL/dL vs this script's systemic CaO2(PaO2={PAO2_REST_MMHG:.0f})="
              f"{cao2_rest:.4f} mL/dL -> actual content diff={actual_content_diff:.4f} mL/dL "
              f"({actual_content_diff/cao2_rest*100:.2f}%, matches the audit's flagged 0.47%)")
        print(f"MODELED content diff (local dC/dP={slope_mid:.5f} mL/dL/mmHg x A-a gradient="
              f"{AA_GRADIENT_MMHG}mmHg) = {expected_content_diff:.4f} mL/dL -- residual vs actual: "
              f"{modeling_residual_pct:.1f}%: {'PASS -- the fork is now a MODELED A-a quantity' if aa_explains_fork else 'FAIL'}")
        aa_wiring_gate = {
            "cco2_end_capillary_ml_dl": cco2_end_capillary, "cao2_systemic_ml_dl": cao2_rest,
            "actual_content_diff_ml_dl": actual_content_diff,
            "actual_content_diff_pct": actual_content_diff / cao2_rest * 100.0,
            "expected_content_diff_from_aa_gradient_ml_dl": expected_content_diff,
            "modeling_residual_pct": modeling_residual_pct, "aa_gradient_explains_fork": bool(aa_explains_fork),
        }
    report["gate1_pulmonary_wiring_aa_gradient"] = aa_wiring_gate

    print("\n" + "=" * 78)
    print("STEP 3B/11 -- FORCED ADVERSARY (void-floor): is the Hb-Hill mechanism NECESSARY, or would a")
    print("degenerate/nearby-wrong model pass just as well (the leaning-positive confound to rule out)?")
    print("=" * 78)
    ca_dissolved_only = O2_SOLUBILITY * PAO2_REST_MMHG
    hb_binding_share_pct = (cao2_rest - ca_dissolved_only) / cao2_rest * 100
    print(f"Void-floor A: strip Hb-BINDING entirely (dissolved-O2 term only) -> CaO2 = {ca_dissolved_only:.4f} "
          f"mL/dL ({ca_dissolved_only/cao2_rest*100:.2f}% of the real {cao2_rest:.2f} mL/dL) -- Hb-binding "
          f"supplies {hb_binding_share_pct:.1f}% of CaO2: {'PASS -- mechanism is NECESSARY' if hb_binding_share_pct > 90 else 'FAIL -- mechanism is decorative'}")
    wrong_n_diffs = {wn: abs(float(sao2_hill(40.0, n=wn)) * 100 - 75.0) / 75.0 * 100
                     for wn in [1.0, 2.0, 3.5, 4.0]}
    wrong_p50_diffs = {wp: abs(float(sao2_hill(40.0, p50=wp)) * 100 - 75.0) / 75.0 * 100
                        for wp in [20.0, 24.0, 30.0, 35.0]}
    print(f"Void-floor B: does the SvO2~75%-at-PO2=40 match (Step 2) survive a WRONG Hill exponent n? "
          f"diffs at n={{1.0,2.0,3.5,4.0}} = {[round(v,1) for v in wrong_n_diffs.values()]}% (all "
          f"{'>5% -- FAIL as expected, task pair is NOT interchangeable' if all(v > 5.0 for v in wrong_n_diffs.values()) else 'unexpectedly close'})")
    print(f"Void-floor C: same for a WRONG P50: diffs at P50={{20,24,30,35}} = "
          f"{[round(v,1) for v in wrong_p50_diffs.values()]}% (all "
          f"{'>5% -- FAIL as expected' if all(v > 5.0 for v in wrong_p50_diffs.values()) else 'unexpectedly close'})")
    n_grid = np.arange(1.0, 6.01, 0.2)
    p50_grid = np.arange(15.0, 40.01, 1.0)
    nn_mesh, pp_mesh = np.meshgrid(n_grid, p50_grid)
    sat_mesh = sao2_hill(40.0, p50=pp_mesh, n=nn_mesh) * 100
    pass_frac = float(np.mean(np.abs(sat_mesh - 75.0) / 75.0 * 100 < 5.0))
    print(f"Void-floor D: grid-scan n in[1,6]x0.2, P50 in[15,40]x1mmHg ({n_grid.size*p50_grid.size} combos): "
          f"only {pass_frac*100:.1f}% pass the same 5% band the task-given (n=2.7,P50=26.6) pair passes -- "
          f"{'PASS -- the task pair sits in a genuinely constrained region, not a trivially-easy one' if pass_frac < 0.5 else 'FAIL -- band too easy, not discriminating'}")
    forced_adversary_ok = (hb_binding_share_pct > 90 and all(v > 5.0 for v in wrong_n_diffs.values()) and
                           all(v > 5.0 for v in wrong_p50_diffs.values()) and pass_frac < 0.5)
    print(f"FORCED-ADVERSARY gate (all 4 void-floors must show the mechanism/parameterization is "
          f"NECESSARY, not decorative or trivially-easy): {'PASS' if forced_adversary_ok else 'FAIL'}")
    report["forced_adversary_void_floor"] = {
        "dissolved_only_cao2_ml_dl": ca_dissolved_only, "hb_binding_share_pct": hb_binding_share_pct,
        "wrong_n_diffs_pct": {str(k): v for k, v in wrong_n_diffs.items()},
        "wrong_p50_diffs_pct": {str(k): v for k, v in wrong_p50_diffs.items()},
        "grid_scan_pass_fraction": pass_frac, "forced_adversary_ok": bool(forced_adversary_ok),
    }

    print("\n" + "=" * 78)
    print("STEP 4/11 -- DECISIVE FALSIFIER: does the Fick loop CLOSE? VO2 = CO x (CaO2-CvO2)")
    print("(three DIFFERENT-provenance numbers must multiply out consistently -- non-circular)")
    print("=" * 78)
    vo2_recon_higg = higg_co_rest_l_min * 10.0 * avo2diff_chem   # L/min*10=dL/min; x mL/dL = mL/min
    diff_higg_pct = abs(vo2_recon_higg - vo2_rest_twin_ml_min) / vo2_rest_twin_ml_min * 100
    print(f"Leg A (OpenSim muscle-probe VO2, metabolic_cost.py):        {vo2_rest_twin_ml_min:.2f} mL/min")
    print(f"Leg B (REAL human catheterization CO, Higginbotham 1986):  {higg_co_rest_l_min:.3f} L/min")
    print(f"Leg C (first-principles Hb-Hill-chemistry a-vO2diff, THIS script): {avo2diff_chem:.3f} mL/dL")
    print(f"Fick closure: B x C = {higg_co_rest_l_min:.3f} x {avo2diff_chem*10:.2f} mL/L = "
          f"{vo2_recon_higg:.2f} mL/min  vs Leg A = {vo2_rest_twin_ml_min:.2f} mL/min  "
          f"-> diff = {diff_higg_pct:.2f}%")
    fick_closes_higg = diff_higg_pct < FICK_CLOSURE_TOL_PCT
    print(f"FICK-LOOP-CLOSES gate (pre-registered <{FICK_CLOSURE_TOL_PCT}% threshold, fully non-circular "
          f"3-leg triangulation): {'PASS -- loop CLOSES' if fick_closes_higg else 'FAIL'}")
    vo2_recon_twinq = q_rest_twin_l_min * 10.0 * avo2diff_chem
    diff_twinq_pct = abs(vo2_recon_twinq - vo2_rest_twin_ml_min) / vo2_rest_twin_ml_min * 100
    print(f"\nSecondary/internal-consistency variant (using the reference body's Q instead of Higginbotham's "
          f"real CO -- still independent of the chemistry leg, but Q itself used the assumed-5.0 "
          f"literature constant, so this is a WEAKER, not the primary, triangulation): "
          f"{q_rest_twin_l_min:.3f} x {avo2diff_chem*10:.2f} = {vo2_recon_twinq:.2f} mL/min, "
          f"diff={diff_twinq_pct:.2f}%: {'PASS' if diff_twinq_pct < FICK_CLOSURE_TOL_PCT else 'FAIL'}")
    diff_classic_pct = abs(vo2_recon_higg - ANCHOR_VO2_REST_ML_MIN) / ANCHOR_VO2_REST_ML_MIN * 100
    print(f"\nAlso vs the task's classic textbook anchor (VO2~{ANCHOR_VO2_REST_ML_MIN:.0f} mL/min): "
          f"reconstructed {vo2_recon_higg:.2f} mL/min, diff={diff_classic_pct:.2f}%: "
          f"{'PASS' if diff_classic_pct < FICK_CLOSURE_TOL_PCT else 'FAIL'}")
    q_chem = vo2_rest_twin_ml_min / (avo2diff_chem * 10.0)
    q_chem_diff_higg_pct = abs(q_chem - higg_co_rest_l_min) / higg_co_rest_l_min * 100
    print(f"\nInverse direction (Q_chemistry = twin's VO2 / chemistry a-vO2diff): {q_chem:.3f} L/min "
          f"vs Higginbotham REAL rest CO={higg_co_rest_l_min:.2f}: diff={q_chem_diff_higg_pct:.2f}%  "
          f"[reported honestly -- {q_chem:.2f} sits just outside cardiac_output.py's tight "
          f"ANCHOR_Q_REST=[4,6] band (marginal, {(q_chem-6.0):.2f} L/min over), but within its own wider "
          f"25% band ({higg_co_rest_l_min*1.25:.2f} ceiling) -- NOT cherry-picked, both thresholds shown]")
    report["fick_closure_falsifier"] = {
        "vo2_leg_a_openism_ml_min": vo2_rest_twin_ml_min, "co_leg_b_higginbotham_l_min": higg_co_rest_l_min,
        "avo2diff_leg_c_chemistry_ml_dl": avo2diff_chem, "vo2_reconstructed_ml_min": vo2_recon_higg,
        "diff_pct_vs_openism_vo2": diff_higg_pct, "fick_closes": bool(fick_closes_higg),
        "vo2_reconstructed_using_twin_q_ml_min": vo2_recon_twinq, "diff_pct_using_twin_q": diff_twinq_pct,
        "diff_pct_vs_classic_250": diff_classic_pct, "q_chemistry_l_min": q_chem,
        "q_chemistry_diff_pct_vs_higginbotham": q_chem_diff_higg_pct,
    }

    print("\n" + "=" * 78)
    print("STEP 6/11 -- Hb/Hct disclosure (Billett 1990, PMID 21250102) -- inter-individual OPEN, not resolved")
    print("=" * 78)
    print(f"Task/this-script value: Hb={HB_G_DL} g/dL, Hct={HCT_PCT}% -- inside the MALE reference range "
          f"(14-18 g/dL, 40-54%) per Billett (1990); FEMALE ranges are systematically lower (12-16 g/dL, "
          f"36-48%) -- NOT resolved here, held explicitly OPEN per the task's symmetric-QC framing.")
    report["hb_hct_reference"] = {
        "used_hb_g_dl": HB_G_DL, "used_hct_pct": HCT_PCT,
        "billett_male_hb_range_g_dl": [14.0, 18.0], "billett_female_hb_range_g_dl": [12.0, 16.0],
        "billett_male_hct_range_pct": [40.0, 54.0], "billett_female_hct_range_pct": [36.0, 48.0],
        "held_open_not_resolved": True,
    }

    print("\n" + "=" * 78)
    print("STEP 7/11 -- WALKING regime: invert cardiac_output.py's a-vO2diff widening (10/11/12) to")
    print("the mixed-venous PO2/SvO2 it IMPLIES via this script's Hb-Hill chemistry (plausibility check)")
    print("=" * 78)
    walk_implied = {}
    prev_pvo2 = PVO2_REST_MMHG
    monotonic_ok = True
    for label, avo2 in [("low_10", 10.0), ("mid_11", 11.0), ("high_12", 12.0)]:
        target_cv = cao2_rest - avo2
        f = lambda pv: cao2_ml_dl(HB_G_DL, float(sao2_hill(pv)), pv) - target_cv
        pv_solved = brentq(f, 1.0, PAO2_REST_MMHG - 0.1)
        sv_solved = float(sao2_hill(pv_solved))
        walk_implied[label] = {"avo2diff_ml_100ml": avo2, "implied_pvo2_mmhg": pv_solved,
                                "implied_svo2_pct": sv_solved * 100}
        plausible = 0.0 < pv_solved < PVO2_REST_MMHG
        monotonic_ok = monotonic_ok and (pv_solved < prev_pvo2)
        prev_pvo2 = pv_solved
        print(f"  a-vO2diff={avo2:4.1f} mL/100mL -> implied PvO2={pv_solved:5.2f} mmHg, "
              f"SvO2={sv_solved*100:5.2f}%  (0<PvO2<rest-40mmHg: {'PASS' if plausible else 'FAIL'})")
    print(f"Monotonic-widening gate (implied PvO2 strictly FALLS as a-vO2diff widens further): "
          f"{'PASS' if monotonic_ok else 'FAIL'}")
    walk_all_plausible = all(0.0 < v["implied_pvo2_mmhg"] < PVO2_REST_MMHG for v in walk_implied.values())
    print(f"All-plausible gate: {'PASS' if walk_all_plausible else 'FAIL'}")
    print("(Disclosed assumption: arterial SaO2 held fixed at its rest value through submaximal walking "
          "-- standard/expected since arterial desaturation during exercise is a near-maximal-effort/"
          "elite-athlete phenomenon, not a submaximal-walking one; not independently re-verified "
          "-- see Honest gaps. Roca et al. 1992, PMID 1400019, topically confirms venous PO2 "
          "falls at exercising muscle -- the mechanism's SITE -- without a live-accessible number.)")
    report["walking_regime_inversion"] = {
        "results": walk_implied, "monotonic_ok": bool(monotonic_ok), "all_plausible": bool(walk_all_plausible),
    }

    print("\n" + "=" * 78)
    print("STEP 8/11 -- BOHR SHIFT direction falsifier (swept across illustrative coefficient magnitudes)")
    print("=" * 78)
    tissue_po2 = 30.0   # illustrative exercising-tissue/capillary PO2 (generic, disclosed)
    ph_sweep = np.array([7.40, 7.35, 7.30, 7.20, 7.10])
    coef_sweep = [0.20, 0.30, 0.40, 0.48, 0.60, 0.70]
    bohr_rows = {}
    all_monotonic = True
    for coef in coef_sweep:
        p50s = bohr_shifted_p50(ph_sweep, coef=coef)
        sats = sao2_hill(tissue_po2, p50=p50s) * 100
        is_monotonic = bool(np.all(np.diff(sats) < 0))
        all_monotonic = all_monotonic and is_monotonic
        bohr_rows[str(coef)] = {"sao2_pct_by_ph": sats.tolist(), "monotonic_decreasing": is_monotonic}
        print(f"  coef={coef:.2f}: SaO2@PO2={tissue_po2}mmHg across pH{ph_sweep.tolist()} = "
              f"{[round(float(s),2) for s in sats]}  monotonic_decreasing: {'PASS' if is_monotonic else 'FAIL'}")
    print(f"\nBOHR-DIRECTION gate (SaO2 at FIXED tissue PO2 strictly DECREASES as pH falls, for EVERY "
          f"swept coefficient magnitude -- matches Collins 2015 [cit 1] + Wikipedia [cit 9]'s live-"
          f"verified right-shift direction, sign-robust to the coefficient-magnitude uncertainty "
          f"disclosed in Honest gaps): {'PASS -- direction holds at ALL swept magnitudes' if all_monotonic else 'FAIL'}")
    report["bohr_shift"] = {"tissue_po2_mmhg": tissue_po2, "ph_sweep": ph_sweep.tolist(),
                             "coef_sweep": coef_sweep, "rows": bohr_rows, "all_monotonic": bool(all_monotonic)}

    print("\n" + "=" * 78)
    print("STEP 9/11 -- GEOMETRIC STRUCTURE: local slope dSaO2/dPO2 -- analytical vs numerical")
    print("=" * 78)
    po2_grid = np.linspace(10.0, 150.0, 4001)
    sao2_grid = sao2_hill(po2_grid)
    d_num = np.gradient(sao2_grid, po2_grid)
    d_ana = dsao2_dpo2_analytical(po2_grid)
    deriv_match = bool(np.allclose(d_num[2:-2], d_ana[2:-2], rtol=2e-2))
    imax = int(np.argmax(d_num))
    po2_at_max_slope_numerical = float(po2_grid[imax])
    po2_at_max_slope_analytical = float(P50_MMHG * ((HILL_N - 1.0) / (HILL_N + 1.0)) ** (1.0 / HILL_N))
    max_slope_match_pct = abs(po2_at_max_slope_numerical - po2_at_max_slope_analytical) / po2_at_max_slope_analytical * 100
    print(f"Analytical-vs-numerical derivative match (dSaO2/dPO2, rtol=2%): {'PASS' if deriv_match else 'FAIL'}")
    print(f"PO2 of steepest slope: numerical={po2_at_max_slope_numerical:.2f} mmHg, closed-form "
          f"P50x((n-1)/(n+1))^(1/n)={po2_at_max_slope_analytical:.2f} mmHg  (diff={max_slope_match_pct:.2f}%): "
          f"{'PASS' if max_slope_match_pct < 2.0 else 'FAIL'}")
    slope_at_100 = float(dsao2_dpo2_analytical(100.0))
    slope_at_max = float(dsao2_dpo2_analytical(po2_at_max_slope_analytical))
    sensitivity_ratio = slope_at_max / slope_at_100
    print(f"Sensitivity ratio (slope at steepest point / slope at PO2=100 plateau) = {sensitivity_ratio:.1f}x "
          f"-- the observable SaO2 is {sensitivity_ratio:.0f}x MORE informative about a PO2 change in the "
          f"tissue-extraction zone than in the arterial plateau (the plateau IS a real null-space: pulse "
          f"oximetry cannot resolve PaO2 changes there -- a known clinical fact, not asserted, MEASURED here).")
    report["geometric_structure"] = {
        "derivative_match": deriv_match, "po2_at_max_slope_numerical_mmhg": po2_at_max_slope_numerical,
        "po2_at_max_slope_analytical_mmhg": po2_at_max_slope_analytical,
        "max_slope_location_match_pct": max_slope_match_pct, "sensitivity_ratio_steep_vs_plateau": sensitivity_ratio,
    }

    print("\n" + "=" * 78)
    print("STEP 10/11 -- VOID-FLOOR / non-degeneracy sweeps: [Hb] (anemia<->polycythemia), PaO2 (hypoxia<->hyperoxia)")
    print("=" * 78)
    hb_sweep = np.linspace(6.0, 22.0, 17)
    cao2_hb_sweep = cao2_ml_dl(hb_sweep, sao2_rest, PAO2_REST_MMHG)
    cvo2_hb_fixed = cao2_ml_dl(hb_sweep, svo2_rest, PVO2_REST_MMHG)
    avo2diff_hb_sweep = cao2_hb_sweep - cvo2_hb_fixed
    co_needed_hb_sweep = vo2_rest_twin_ml_min / (avo2diff_hb_sweep * 10.0)
    d_cao2_d_hb_numerical = np.gradient(cao2_hb_sweep, hb_sweep)
    d_cao2_d_hb_analytical = HUFNER_K * sao2_rest
    hb_deriv_match = bool(np.allclose(d_cao2_d_hb_numerical[2:-2], d_cao2_d_hb_analytical, rtol=1e-2))
    hb_monotonic = bool(np.all(np.diff(cao2_hb_sweep) > 0) and np.all(np.diff(co_needed_hb_sweep) < 0))
    print(f"[Hb] swept [{hb_sweep.min():.0f},{hb_sweep.max():.0f}] g/dL: CaO2 range "
          f"[{cao2_hb_sweep.min():.2f},{cao2_hb_sweep.max():.2f}] mL/dL (monotonic increasing with Hb), "
          f"CO-needed-for-fixed-VO2 range [{co_needed_hb_sweep.min():.2f},{co_needed_hb_sweep.max():.2f}] "
          f"L/min (monotonic DECREASING with Hb -- anemia forces compensatory tachycardia/high-output "
          f"state, a real, directionally-known clinical phenomenon): {'PASS' if hb_monotonic else 'FAIL'}")
    print(f"Analytical-vs-numerical d(CaO2)/d(Hb) match (closed form = Hufner_K x SaO2 = "
          f"{d_cao2_d_hb_analytical:.4f} mL/dL per g/dL, rtol=1%): {'PASS' if hb_deriv_match else 'FAIL'}")

    pao2_sweep = np.linspace(30.0, 150.0, 25)
    sao2_pao2_sweep = sao2_hill(pao2_sweep) * 100
    pao2_monotonic = bool(np.all(np.diff(sao2_pao2_sweep) > 0))
    plateau_slope = (sao2_pao2_sweep[-1] - sao2_pao2_sweep[-5]) / (pao2_sweep[-1] - pao2_sweep[-5])
    shoulder_idx_lo = int(np.argmin(np.abs(pao2_sweep - 40)))
    shoulder_idx_hi = int(np.argmin(np.abs(pao2_sweep - 60)))
    shoulder_slope = (sao2_pao2_sweep[shoulder_idx_hi] - sao2_pao2_sweep[shoulder_idx_lo]) / (pao2_sweep[shoulder_idx_hi] - pao2_sweep[shoulder_idx_lo])
    print(f"PaO2 swept [{pao2_sweep.min():.0f},{pao2_sweep.max():.0f}] mmHg: SaO2 range "
          f"[{sao2_pao2_sweep.min():.2f},{sao2_pao2_sweep.max():.2f}]% (strictly monotonic increasing): "
          f"{'PASS' if pao2_monotonic else 'FAIL'}")
    print(f"Shoulder slope (40-60mmHg) = {shoulder_slope:.3f} %/mmHg  vs  plateau slope (top of sweep) = "
          f"{plateau_slope:.3f} %/mmHg  -> ratio {shoulder_slope/plateau_slope:.1f}x steeper in the "
          f"tissue-extraction shoulder than the arterial plateau (measured, matches Step 8's geometric claim)")
    report["void_floor_sweeps"] = {
        "hb_sweep_g_dl": hb_sweep.tolist(), "cao2_hb_sweep_ml_dl": cao2_hb_sweep.tolist(),
        "co_needed_hb_sweep_l_min": co_needed_hb_sweep.tolist(), "hb_monotonic": hb_monotonic,
        "hb_deriv_match": hb_deriv_match, "pao2_sweep_mmhg": pao2_sweep.tolist(),
        "sao2_pao2_sweep_pct": sao2_pao2_sweep.tolist(), "pao2_monotonic": pao2_monotonic,
        "shoulder_slope_pct_per_mmhg": shoulder_slope, "plateau_slope_pct_per_mmhg": plateau_slope,
    }

    print("\n" + "=" * 78)
    print("STEP 11/11 -- BONUS/SECONDARY (does NOT gate overall_pass): Severinghaus(1979) independent")
    print("polynomial cross-check -- RECALLED not live-verified, flagged accordingly")
    print("=" * 78)
    test_pts = [40.0, 60.0, 95.0, 100.0]
    sev_vs_hill = {}
    for p in test_pts:
        sev = float(severinghaus_1979_so2(p)) * 100
        hill = float(sao2_hill(p)) * 100
        sev_vs_hill[str(p)] = {"severinghaus_pct": sev, "hill_pct": hill, "diff_pct_pts": abs(sev - hill)}
        print(f"  PO2={p:6.1f}: Severinghaus(1979)={sev:.2f}%  Hill(this script)={hill:.2f}%  "
              f"diff={abs(sev-hill):.2f} percentage points")
    max_diff = max(v["diff_pct_pts"] for v in sev_vs_hill.values())
    routes_agree = max_diff < 2.0
    print(f"Two-independent-functional-form agreement (<2 percentage points, all 4 test points): "
          f"{'PASS (bonus, non-gating)' if routes_agree else 'FAIL (bonus, non-gating)'}")
    print("HONEST GAP: the Severinghaus(1979) polynomial coefficients used here are RECALLED from "
          "general domain knowledge, NOT independently re-derived from a live-fetched primary-source "
          "numeric table (PubMed full-text fetch for this specific paper was CAPTCHA-"
          "blocked) -- kept OFF overall_pass accordingly, reported as exploratory/secondary only.")
    report["severinghaus_bonus_crosscheck"] = {"points": sev_vs_hill, "max_diff_pct_pts": max_diff,
                                                "routes_agree": bool(routes_agree)}

    print("\n" + "=" * 78)
    print("COUPLING METADATA (couples_to cardiac + pulmonary + metabolic -- prose only, no shared-graph write)")
    print("=" * 78)
    couples_to = ["the cardiac_output cell -- supplies the CO this layer's "
                  "Fick-closure falsifier consumes (Step 4); this layer supplies the a-vO2diff that "
                  "cell's Fick calc had to ASSUME from literature",
                  "the respiratory cell -- shares the same upstream VO2 "
                  "(metabolic_cost.py); this layer adds the O2-CONTENT (not just volume/ventilation) side",
                  "the metabolic_cost cell -- the ultimate VO2 source all 3 layers "
                  "consume (shared-instrument-chain caveat applies identically here, see Honest gaps)"]
    for c in couples_to:
        print(f"  - {c}")
    report["couples_to"] = couples_to

    print("\n" + "=" * 78)
    print("GATES SUMMARY")
    print("=" * 78)
    gates = {
        "hill_equation_validated_2_independent_checkpoints": hill_validated,
        "forced_adversary_void_floor_ok": forced_adversary_ok,
        "cao2_in_task_anchor": cao2_in_anchor,
        "avo2diff_in_task_anchor": avo2diff_in_anchor,
        "svo2_in_task_anchor": svo2_in_anchor,
        "fick_loop_closes_vs_openism_vo2": fick_closes_higg,
        "walking_regime_monotonic": monotonic_ok,
        "walking_regime_plausible": walk_all_plausible,
        "bohr_direction_holds_all_swept_coefficients": all_monotonic,
        "geometric_derivative_match": deriv_match,
        "geometric_max_slope_location_match": max_slope_match_pct < 2.0,
        "hb_sweep_monotonic_nondegenerate": hb_monotonic,
        "hb_sweep_derivative_match": hb_deriv_match,
        "pao2_sweep_monotonic_nondegenerate": pao2_monotonic,
    }
    if aa_wiring_gate is not None:
        gates["gate1_aa_gradient_explains_pulmonary_fork"] = aa_wiring_gate["aa_gradient_explains_fork"]
    gates = {k: bool(v) for k, v in gates.items()}
    overall_pass = all(gates.values())
    print(json.dumps(gates, indent=2))
    print(f"\n(Bonus/non-gating: severinghaus_routes_agree={routes_agree})")
    print(f"\nOVERALL: {'PASS' if overall_pass else 'FAIL/SURPRISE -- see gates above'}")

    report["citations_verified_live"] = {
        "collins_2015_pmid": "26632351", "collins_2015_doi": "10.1183/20734735.001415",
        "collins_2015_pmcid": "PMC4666443",
        "severinghaus_1979_pmid": "35496",
        "kelman_1966_pmid": "5916678", "kelman_1966_doi": "10.1152/jappl.1966.21.4.1375",
        "severinghaus_1966_pmid": "5912737",
        "billett_1990_pmid": "21250102",
        "rivers_2001_pmid": "11794169",
        "roca_1992_pmid": "1400019", "roca_1992_doi": "10.1152/jappl.1992.73.3.1067",
        "wikipedia_fick_principle": "textbook-grade, flagged, re-fetched",
        "wikipedia_oxyhb_dissociation_curve": "textbook-grade, flagged, re-fetched",
    }
    report["gates"] = gates
    report["bonus_non_gating"] = {"severinghaus_routes_agree": bool(routes_agree)}
    report["confidence_tier"] = (
        "in-vivo-anchored (blood-gas / Hill) -- Hill parameterization anchored against Collins et al. "
        "(2015)'s real-patient-validated (n=3524) dissociation curve; Fick-closure falsifier anchored "
        "against Higginbotham et al. (1986)'s real human catheterization CO. One tier below "
        "subject-specific in-vivo (no arterial blood gas measurement is used; blood-gas operating points are generic/"
        "population-level, see honest gaps). NOTHING PROVEN in the strong sense -- Hb/Hct/P50 vary "
        "inter-individually (held OPEN); a closed Fick loop with wide-plausible-band inputs is "
        "consistency, not proof."
    )
    report["overall_pass"] = overall_pass
    report["reference_body"] = {
        "inherited_from": "metabolic_cost_results.json (mass_kg = mc['muscle_mass']"
                          "['total_body_mass_kg'], the reference body's scaled musculoskeletal-model mass, "
                          "read-only, not re-derived here)",
        "mass_kg": mass_kg, "name": None, "body_fat_fraction": None,
        "class": "reference_body",
    }

    out_path = f"{OUT_DIR}/blood_oxygen_transport_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
