"""CEREBRAL BLOOD FLOW AUTOREGULATION: couples the arterial-pressure thread (the arterial_pressure
cell, MAP=CO*TPR) and the acid-base thread (the acid_base_co2 cell, PaCO2) into brain perfusion --
resting CBF (~50 mL/100g/min whole-brain, ~80 grey, ~20 white), the pressure-autoregulation plateau
(CBF roughly constant over CPP ~60-150 mmHg via myogenic + metabolic + neurogenic mechanisms),
CPP = MAP - ICP, and CO2 reactivity (CBF rises with PaCO2).

QUESTION (pre-registered falsifier): does the model reproduce (a) the measured autoregulation curve --
CBF roughly flat over 60-150 mmHg CPP, falling steeply below the lower limit (a MAP=50 perturbation,
below the classic knee, MUST drop CBF) -- and (b) the measured CO2-reactivity slope (~3-4%/mmHg,
checked against TCD/BOLD literature)? Symmetric QC, stated up front: the classic Lassen (1959) plateau
is now known from modern reanalysis to be neither as flat nor as universal as the textbook figure
suggests -- held OPEN and quantified below, not smoothed over.

METHOD, GEOMETRIC STRUCTURE (derive from the geometry, not heuristics):
  1. Pressure autoregulation is modeled as an actively-REGULATED vascular resistance R(P) that tracks
     R_ideal(P)=P/F0 (the resistance value that would hold flow EXACTLY at F0 at every pressure P) but
     is physically bounded to R in [R_min, R_max] (the vessel cannot dilate/constrict without limit).
     R_min = LLA/F0, R_max = UBP/F0 (continuity at the knees, by construction). This gives, via
     F(P)=P/R(P):
       F(P) = F0                      for LLA<=P<=UBP   (the idealized flat "plateau")
       F(P) = F0*(P/LLA)               for P<LLA          (pressure-passive line through the origin)
       F(P) = F0*(P/UBP)               for P>UBP          (pressure-passive line through the origin)
     A clipped-line / saturating-gain-scheduled-controller construction, not a curve fit -- the SAME
     "clip a regulation line to a physical floor/ceiling" geometric family used twice more below (CO2
     reactivity, Step 7; and structurally analogous to arterial_pressure.py's R/C void-floor clamps).
  2. The MODERN, "leaky-plateau" correction integrates Numan et al 2014's measured regression
     elasticities (dlnF/dlnP, i.e. %DeltaCBF/%DeltaMAP) directly, rather than assuming a flat plateau:
     F(P) = F0*(P/P_ref)^e(P), e(P)=slope_decrease for P<P_ref, slope_increase for P>=P_ref -- the
     textbook definition of "integrate a measured elasticity" (log-log linear = a power law), not a
     heuristic curve shape. This is the falsifiable, quantified "modern vs classic" comparison the task
     asks be held open.
  3. CO2 reactivity uses the SAME clipped-line family: linear above a threshold near the eucapnic
     reference (matching Battisti-Charbonney et al 2011's reported LINEAR supra-threshold regime,
     live-quoted slope 5.5%/mmHg by TCD), clipped to a floor below it (matching their own reported
     SIGMOIDAL/saturating sub-threshold shape) -- never an unbounded line (Step 7's forced adversary
     shows why: an unbounded linear model goes NEGATIVE, i.e. unphysical, at realistic hypocapnic PaCO2).
  4. FORCED ADVERSARIES (void-floor, the sharpest, cleanest tests): (a) "no autoregulation" -- constant
     resistance fixed at the resting operating point, extended everywhere -- vs the regulated model,
     over the full 60-150mmHg range (Step 5); (b) "ICP-blind" -- using MAP as if it were CPP, ignoring
     ICP entirely -- vs the correct CPP=MAP-ICP model, evaluated at a pathologically raised ICP (Step 9);
     (c) "zero CO2 reactivity" vs the real CO2 model (Step 8); (d) an unbounded-linear CO2 adversary that
     goes non-physical (negative flow) at extreme hypocapnia (Step 7).
  5. Non-circularity / cross-layer coupling, no re-solve: MAP is read directly from the
     already-computed arterial_pressure_results.json (classic_map_mmhg=93.33, Windkessel true mean=
     91.89) -- never re-derived here. A PaCO2 perturbation (40->50mmHg) is read directly from this
     already-computed acid_base_co2_results.json (its own step3 acute-respiratory scenario,
     delta_paco2=10) -- never invented fresh. Both couplings are load-bearing data reads, not prose.

CITATIONS -- every PMID/DOI verified LIVE via NCBI eutils (esearch/esummary/efetch) and
EuropePMC (machine-readable abstractText field), NOT recalled. Two numbers (marked below) were obtained
only via a WebFetch summarization pass over PMC full text, because direct curl access to PMC full-text
pages was reCAPTCHA-blocked (a real, disclosed tool limitation, not silently hidden) --
those two numbers are flagged at a LOWER confidence tier throughout, never silently treated as
independently machine-cross-checked:
  [1] Lassen NA (1959). "Cerebral blood flow and oxygen consumption in man." Physiol Rev 39(2):183-238.
      PMID 13645234, DOI 10.1152/physrev.1959.39.2.183 (verified live via NCBI efetch: title/journal/
      vol/pages/author match exactly). THE classic paper compiling ~7 independent studies' pooled data
      into the textbook "flat 60-150mmHg" autoregulation curve. Pre-1975, no abstract in PubMed
      (pre-abstracting era) -- PROVENANCE ONLY (same disclosed-gap tier the cell set already
      applies to the 1914 Christiansen/Douglas/Haldane paper in acid_base_co2.py): the curve's classic
      60-150mmHg numbers are used here as the TASK'S OWN stated textbook target to test against, not
      re-extracted from this paper's tables.
  [2] Kety SS, Schmidt CF (1948). "The nitrous oxide method for the quantitative determination of
      cerebral blood flow in man: theory, procedure and normal values." J Clin Invest 27(4):476-83.
      PMID 16695568, DOI 10.1172/JCI101994, PMC439518 (verified live). THE original quantitative
      resting-CBF measurement paper (introduced the N2O Fick-principle method used for the next ~30
      years of cerebral physiology). Pre-abstracting era, no abstract text live-extractable this
      no accessible abstract -- PROVENANCE ONLY, same disclosed-gap tier as [1].
  [3] Kety SS, Schmidt CF (1948). "The effects of altered arterial tensions of carbon dioxide and
      oxygen on cerebral blood flow and cerebral oxygen consumption of normal young men." J Clin Invest
      27(4):484-92. PMID 16695569, DOI 10.1172/JCI101995, PMC439519 (verified live). THE original CO2-
      reactivity quantification paper. Pre-abstracting era -- PROVENANCE ONLY, same disclosed-gap tier.
  [4] Numan T, Bain AR, Hoiland RL, Smirl JD, Lewis NC, Ainslie PN (2014). "Static autoregulation in
      humans: a review and reanalysis." Med Eng Phys 36(11):1487-95. PMID 25205587, DOI
      10.1016/j.medengphy.2014.08.001 (verified live, full abstract fetched via NCBI efetch). THE
      key modern reanalysis: "based on 40 studies (49 individual experimental protocols)... linear
      regression coefficient between MAP and CBF... of 0.82+/-0.77%DeltaCBF/%DeltaMAP during DECREASES
      in MAP (n=23)... 0.21+/-0.47%DeltaCBF/%DeltaMAP during INCREASES (n=26; p<0.001)... After
      correction for increases/decreases in PaCO2, the slopes were not significantly different:
      0.64+/-1.16 (n=16) and 0.39+/-0.30 (n=12)... (p=0.60)." A real, quantified, non-zero "leak" in
      the textbook plateau -- much of the apparent asymmetry/hysteresis is itself a PaCO2 confound, a
      genuine methodological subtlety in the primary modern literature, disclosed not hidden.
  [5] Willie CK, Tzeng YC, Fisher JA, Ainslie PN (2014). "Integrative regulation of human brain blood
      flow." J Physiol 592(5):841-59. PMID 24396059, PMC3948549, DOI 10.1113/jphysiol.2013.268953
      (verified live via EuropePMC's machine-readable abstractText field -- NOT a summarizer).
      Abstract states, verbatim, as ONE OF FOUR key theses of the whole review: "cerebral autoregulation
      does NOT maintain constant perfusion through a mean arterial pressure range of 60-150mmHg" --
      explicitly named among "evidence against several canonized paradigms of CBF control." This is the
      primary literature's stated position, not this script's interpretation.
      LOWER-CONFIDENCE ITEM (WebFetch summarization only, PMC full text; NOT independently raw-verified
      -- direct curl to PMC was reCAPTCHA-blocked): the review reportedly states CO2
      reactivity as "approximate 3-6% increase and/or 1-3% decrease in flow per mmHg change in CO2
      above and below eupnoeic PaCO2, respectively," and that "CO2 reactivity of the microvasculature in
      grey matter is greater than that of white matter." Used only for the hypocapnic-side slope and the
      regional-CO2-reactivity note, both explicitly flagged at reduced confidence throughout.
  [6] Willie CK, Colino FL, Bailey DM, et al (2011). "Utility of transcranial Doppler ultrasound for the
      integrative assessment of cerebrovascular function." J Neurosci Methods 196(2):221-37. PMID
      21276818, DOI 10.1016/j.jneumeth.2011.01.011 (verified live). TCD methodology review, topical/
      provenance support for using TCD velocity as a CBF proxy throughout this script.
  [7] Battisti-Charbonney A, Fisher J, Duffin J (2011). "The cerebrovascular response to carbon dioxide
      in humans." J Physiol 589(Pt 12):3039-48. PMID 21521758, PMC3139085, DOI
      10.1113/jphysiol.2011.206052 (verified live, full abstract fetched). PRIMARY, quantitative, live-
      quoted: "During rebreathing the MCAv response to CO2 was SIGMOIDAL below a discernible threshold
      CO2 tension... The sigmoidal MCAv response was centred at a CO2 tension close to normal resting
      values (overall mean 36mmHg)... Above this threshold both MCAv and MAP increased LINEARLY with
      CO2 tension... overall mean slopes 5.5% mmHg^-1 and 2.1mmHg mmHg^-1, respectively." This is the
      PRIMARY, TCD-measured, live-verified anchor for this script's hypercapnic CO2-reactivity slope
      (5.5%/mmHg) and its two-regime (sigmoidal-below/linear-above) shape.
  [8] Mutch WA, Mandell DM, Fisher JA, Mikulis DJ, Crawley AP, Pucci O, Duffin J (2012). "Approaches to
      brain stress testing: BOLD magnetic resonance imaging with computer-controlled delivery of carbon
      dioxide." PLoS One 7(11):e47443. PMID 23139743, PMC3489910, DOI 10.1371/journal.pone.0047443
      (verified live, full abstract fetched). Live-quoted: "Changing the carbon dioxide (CO2) tension in
      arterial blood is commonly used as a cerebral vasoactive stimulus to assess the cerebral vascular
      response, changing cerebral blood flow (CBF) by up to 5-11 PERCENT/MMHG in normal adults." An
      INDEPENDENT (BOLD-context) secondary source corroborating that real CO2 reactivity runs HIGHER
      than the task's stated ~3-4%/mmHg band -- both [7] and [8] agree on this, a genuine, disclosed
      tension, not laundered into a false match.
  [9] Ito H, Kanno I, Kato C, Sasaki T, Ishii K, Ouchi Y, Iida A, Okazawa H, Hayashida K, Tsuyuguchi N,
      Ishii K, Kuwabara Y, Senda M (2004). "Database of normal human cerebral blood flow, cerebral
      blood volume, cerebral oxygen extraction fraction and cerebral metabolic rate of oxygen measured
      by positron emission tomography with 15O-labelled carbon dioxide or water, carbon monoxide and
      oxygen: a multicentre study in Japan." Eur J Nucl Med Mol Imaging 31(5):635-43. PMID 14730405,
      DOI 10.1007/s00259-003-1430-8 (verified live, full abstract fetched). PRIMARY, quantitative,
      live-quoted: "Subjects comprised 70 healthy volunteers... Overall mean+/-SD values for cerebral
      CORTICAL regions were: CBF=44.4+/-6.5 ml 100ml^-1 min^-1" -- a MODERN PET measurement, used as an
      external cross-check against the task's stated "grey ~80" figure (a real, disclosed
      undershoot, Step 2).
  [10] Alisch JSR, Khattar N, Kim RW, Cortina LE, Rejimon AC, Qian W, Ferrucci L, Resnick SM, Spencer RG,
      Bouhrara M (2021). "Sex and age-related differences in cerebral blood flow investigated using
      pseudo-continuous arterial spin labeling magnetic resonance imaging." Aging (Albany NY)
      13(4):4911-4925. PMID 33596183, PMC7950235, DOI 10.18632/aging.202673 (verified live, full
      abstract fetched via NCBI efetch -- N=80 cognitively-unimpaired subjects, wide age range, GM<WM
      CBF-with-age divergence, "reference CBF values for the standard ASL protocol recommended by the
      ISMRM Perfusion Study Group").
      LOWER-CONFIDENCE ITEM (WebFetch summarization of PMC Table 3 only; NOT independently raw-verified
      -- same PMC block as [5]): overall cohort mean GM CBF=33.0+/-6.17, WM CBF=18.0+/-4.13 ml/100g/min.
      Used as a second, independent (ASL, not PET) modern cross-check against the task's grey/white
      bands, flagged at reduced confidence for the exact decimal values (the abstract's qualitative
      finding -- GM>WM, GM falls with age, WM rises with age -- IS independently, live, primary-verified).
  [11] Czosnyka M, Pickard JD (2004). "Monitoring and interpretation of intracranial pressure." J Neurol
      Neurosurg Psychiatry 75(6):813-21. PMID 15145991, PMC1739058, DOI 10.1136/jnnp.2003.033126
      (verified live, full abstract fetched). Topical/provenance support: "Information which can be
      derived from ICP and its waveforms includes cerebral perfusion pressure (CPP), regulation of
      cerebral blood flow and volume..." -- the clinical-neurosurgery primary-literature anchor for
      CPP-guided therapy as a real, actionable, measured quantity (not this script's invention).
  [12] Mount CA, Das JM. "Cerebral Perfusion Pressure." StatPearls [Internet]. Treasure Island (FL):
      StatPearls Publishing; 2026 Jan-. PMID 30725956 (verified live, full text fetched -- StatPearls
      entries have no separate PMC/DOI, same textbook-grade-but-curated tier the cell set already uses,
      e.g. acid_base_co2.py's NBK507807/NBK536919/NBK532988). Quoted verbatim: "Cerebral perfusion
      pressure (CPP) is... the difference between the mean arterial pressure (MAP) and the intracranial
      pressure (ICP)... Normal CPP lies between 60 and 80mmHg."
  [13] Pinto VL, Adeyinka A. "Increased Intracranial Pressure." StatPearls [Internet]. PMID 29489250
      (verified live, full text fetched). Quoted verbatim: "Normal intracranial pressure (ICP) in adults
      typically ranges from 7 to 15mmHg in the supine position. Values above 20 to 25mmHg are generally
      considered pathological." Also states the Monro-Kellie doctrine (fixed cranial volume) and that
      "increased ICP... can reduce cerebral perfusion pressure (CPP)... ultimately cause ischemia."
  [14] Silverman A, Petersen NH. "Physiology, Cerebral Autoregulation." StatPearls [Internet]. PMID
      31985976 (verified live, full text fetched). Quoted verbatim: distinguishes autoregulation from
      "carbon dioxide reactivity" and "flow-metabolism coupling" as three SEPARATE mechanisms; notes
      cerebrovascular resistance follows "the Hagen-Poiseuille equation"; states that "researchers have
      developed technology that now boasts the ability to measure autoregulatory function in real-time,
      which may lead to fine-tuning long-established guidelines" and flags "individualizing cerebral
      perfusion pressure targets" as an active frontier -- directly supports this script's decision
      to hold LLA/UBP individual variability OPEN rather than assert one universal value.
  [15] Rangel-Castilla L, Gasco J, Nauta HJ, Okonkwo DO, Robertson CS (2008). "Cerebral pressure
      autoregulation in traumatic brain injury." Neurosurg Focus 25(4):E7. PMID 18828705, DOI
      10.3171/FOC.2008.25.10.E7 (verified live, full abstract fetched). Quoted verbatim: "cerebral
      autoregulation is more a CONCEPT than a physically measurable entity... Alterations in cerebral
      autoregulation can vary from patient to patient and over time... assessment of cerebral
      autoregulation is best achieved with DYNAMIC autoregulation methods" -- the STATIC-vs-DYNAMIC
      distinction the task requires be held open (this script models STATIC autoregulation only, Sec.0).

READS: <BODYTWIN_OUT>/arterial_pressure/arterial_pressure_results.json (MAP) and
       <BODYTWIN_OUT>/acid_base_co2/acid_base_co2_results.json (a PaCO2 perturbation scenario).
WRITES: <BODYTWIN_OUT>/cerebral_autoregulation/cerebral_autoregulation_results.json
GATE: overall_pass = all gates in the GATES block; exit 0 on pass, 2 otherwise.
"""
import json
import os
import sys

import numpy as np

# --------------------------------------------------------------------------- paths / consts --
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
ARTERIAL_PRESSURE_JSON = _os.path.join(OUT_ROOT, "arterial_pressure", "arterial_pressure_results.json")
ACID_BASE_JSON = _os.path.join(OUT_ROOT, "acid_base_co2", "acid_base_co2_results.json")
OUT_DIR = _os.path.join(OUT_ROOT, "cerebral_autoregulation")
OUT_PATH = _os.path.join(OUT_DIR, "cerebral_autoregulation_results.json")

# ---- task's pre-registered targets (stated BEFORE any number below is computed) ---------------
CBF_WHOLE_TASK = 50.0          # mL/100g/min, resting whole-brain
CBF_GREY_TASK = 80.0           # mL/100g/min, resting grey matter
CBF_WHITE_TASK = 20.0          # mL/100g/min, resting white matter
BAND_TOL_REL = 0.15            # +/-15% tolerance band around each task-stated central value

LLA_CLASSIC = 60.0             # mmHg, classic Lassen [1] lower limit of autoregulation
UBP_CLASSIC = 150.0            # mmHg, classic Lassen [1] upper limit of autoregulation
P_REF = 93.33                  # mmHg, reference/resting operating point (arterial_pressure.py's classic MAP)

CO2_REACT_TASK_LOW, CO2_REACT_TASK_HIGH = 3.0, 4.0   # %/mmHg, task's pre-registered band
PACO2_REF = 40.0                # mmHg, eucapnic reference (matches acid_base_co2.py's normal_paco2_mmhg)

ICP_NORMAL_LOW, ICP_NORMAL_HIGH = 7.0, 15.0     # mmHg, StatPearls [13]
ICP_PATHOLOGICAL = 25.0                          # mmHg, StatPearls [13] "20-25 pathological" upper example
CPP_NORMAL_LOW, CPP_NORMAL_HIGH = 60.0, 80.0     # mmHg, StatPearls [12]

# ---- live-verified literature slopes/values -------------------------------------------------------
NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_DECREASE_SD = 0.82, 0.77   # %DeltaCBF/%DeltaMAP, MAP decreases [4]
NUMAN_SLOPE_INCREASE, NUMAN_SLOPE_INCREASE_SD = 0.21, 0.47   # %DeltaCBF/%DeltaMAP, MAP increases [4]
LEAK_DEVIATION_GATE_PCT = 5.0    # pre-registered: modern model must deviate >5% from flat somewhere in-plateau

BATTISTI_HYPER_SLOPE_PCT = 5.5    # %/mmHg, TCD, PRIMARY live-quoted [7]
MUTCH_RANGE_LOW, MUTCH_RANGE_HIGH = 5.0, 11.0     # %/mmHg, BOLD-context secondary [8]
WILLIE_HYPOCAPNIC_SLOPE_PCT = 2.0   # %/mmHg midpoint of WebFetch-flagged 1-3% range [5], LOWER CONFIDENCE
CO2_G_FLOOR = 0.5                    # illustrative floor multiplier at extreme hypocapnia (generic, disclosed)

ITO_CORTICAL_CBF, ITO_CORTICAL_SD = 44.4, 6.5     # mL/100g/min, PET, PRIMARY live-quoted [9]
ALISCH_GM_CBF, ALISCH_GM_SD = 33.0, 6.17          # mL/100g/min, ASL, WebFetch-flagged [10]
ALISCH_WM_CBF, ALISCH_WM_SD = 18.0, 4.13          # mL/100g/min, ASL, WebFetch-flagged [10]

MAP_PERTURBATION_TEST = 50.0     # mmHg, the task's explicit falsifier requirement
MIN_DROP_PCT_GATE = 10.0         # pre-registered minimum %CBF drop required at the MAP=50 perturbation

PACO2_COUPLED_DELTA = 10.0        # mmHg, reused verbatim from acid_base_co2.py's step3 scenario (no new number invented)


# --------------------------------------------------------------------------- model functions --
def classic_autoreg_flow(p, f0, lla, ubp):
    """Idealized clipped-resistance ('Lassen textbook') curve: flat in [lla,ubp], pressure-passive
    lines through the origin outside. Derived from R(P)=clip(P/f0, lla/f0, ubp/f0), F(P)=P/R(P)."""
    p = np.asarray(p, dtype=float)
    r_min, r_max = lla / f0, ubp / f0
    r = np.clip(p / f0, r_min, r_max)
    return p / r


def leaky_autoreg_flow(p, f0, lla, ubp, p_ref, slope_dec, slope_inc):
    """Modern 'leaky-plateau' curve: integrates Numan et al 2014's measured regression elasticities
    (dlnF/dlnP) outward from p_ref to the knees (a power law, F=F0*(P/p_ref)^slope), then continues as a
    pressure-passive line beyond lla/ubp (continuity enforced at both knees)."""
    p = np.atleast_1d(np.asarray(p, dtype=float))
    f = np.empty_like(p)
    in_plateau_lo = (p >= lla) & (p < p_ref)
    in_plateau_hi = (p >= p_ref) & (p <= ubp)
    below = p < lla
    above = p > ubp
    f[in_plateau_lo] = f0 * (p[in_plateau_lo] / p_ref) ** slope_dec
    f[in_plateau_hi] = f0 * (p[in_plateau_hi] / p_ref) ** slope_inc
    f_at_lla = f0 * (lla / p_ref) ** slope_dec
    f_at_ubp = f0 * (ubp / p_ref) ** slope_inc
    f[below] = f_at_lla * (p[below] / lla)
    f[above] = f_at_ubp * (p[above] / ubp)
    return f


def void_floor_constant_r_flow(p, f0, p_ref):
    """Forced adversary: 'no autoregulation' -- resistance permanently fixed at the value it takes at
    the resting operating point (R=p_ref/f0), extended to every pressure. Purely passive Ohm's law."""
    p = np.asarray(p, dtype=float)
    r_fixed = p_ref / f0
    return p / r_fixed


def co2_multiplier(c, c0, s_hyper, s_hypo, g_floor):
    """CO2-reactivity multiplier g(C), g(c0)=1: linear above c0 (matches Battisti-Charbonney's [7] own
    reported LINEAR supra-threshold regime), linear-then-clamped-to-a-floor below c0 (matches their own
    reported SIGMOIDAL/saturating sub-threshold shape) -- the SAME clipped-line family as Step 1's
    pressure curve."""
    c = np.atleast_1d(np.asarray(c, dtype=float))
    g = np.where(
        c >= c0,
        1.0 + s_hyper * (c - c0),
        np.maximum(g_floor, 1.0 + s_hypo * (c - c0)),
    )
    return g


def co2_multiplier_naive_unbounded_adversary(c, c0, s_hyper):
    """Forced adversary: the SAME hypercapnic linear slope applied in BOTH directions, unbounded, no
    floor -- the 'simplest possible' model one would reach for without the sigmoid/floor correction."""
    c = np.atleast_1d(np.asarray(c, dtype=float))
    return 1.0 + s_hyper * (c - c0)


def combined_cbf(p, c, f0, lla, ubp, p_ref, slope_dec, slope_inc, c0, s_hyper, s_hypo, g_floor,
                  use_leaky=True):
    """Joint model: F(P,C) = F_autoreg(P) * g(C) -- both mechanisms independently adjust a SHARED
    resistance state (disclosed structural simplification, Sec. Honest Gaps: real physiology has
    documented pressure-CO2 interaction effects this multiplicative-independence form does not capture)."""
    f_p = (leaky_autoreg_flow(p, f0, lla, ubp, p_ref, slope_dec, slope_inc) if use_leaky
           else classic_autoreg_flow(p, f0, lla, ubp))
    g_c = co2_multiplier(c, c0, s_hyper, s_hypo, g_floor)
    return f_p * g_c


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def in_band(x, center, rel_tol):
    lo, hi = center * (1 - rel_tol), center * (1 + rel_tol)
    return bool(lo <= x <= hi), (lo, hi)


def main():
    if not os.path.exists(ARTERIAL_PRESSURE_JSON) or not os.path.exists(ACID_BASE_JSON):
        print("FAIL: required upstream JSON(s) missing -- run arterial_pressure.py and "
              "acid_base_co2.py first.")
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    report = {}

    print("=" * 78)
    print("STEP 1/11 -- load the already-computed MAP (arterial-pressure) and PaCO2 perturbation "
          "(acid-base) -- NO re-solve, two decorrelated upstream layers")
    print("=" * 78)
    art = load_json(ARTERIAL_PRESSURE_JSON)
    acidbase = load_json(ACID_BASE_JSON)
    map_classic = art["map_approximation"]["classic_map_mmhg"]
    map_windkessel_true_mean = art["pulse_pressure"]["windkessel_sim"]["mean_p_analytical"]
    paco2_normal = acidbase["constants"]["normal_paco2_mmhg"]
    paco2_perturbed = acidbase["step3_acute_respiratory"][1]["paco2"]  # delta_paco2=10 entry, paco2=50
    delta_paco2_from_sibling = acidbase["step3_acute_respiratory"][1]["delta_paco2"]
    print(f"arterial_pressure.py's MAP: classic={map_classic:.2f} mmHg, "
          f"Windkessel-true-mean={map_windkessel_true_mean:.2f} mmHg")
    print(f"acid_base_co2.py's acute-respiratory scenario: PaCO2 {paco2_normal:.0f} -> "
          f"{paco2_perturbed:.0f} mmHg (delta={delta_paco2_from_sibling:.0f}, reused verbatim, not invented)")
    assert abs(delta_paco2_from_sibling - PACO2_COUPLED_DELTA) < 1e-9
    report["upstream_inputs"] = {
        "map_classic_mmhg": map_classic, "map_windkessel_true_mean_mmhg": map_windkessel_true_mean,
        "paco2_normal_mmhg": paco2_normal, "paco2_perturbed_mmhg": paco2_perturbed,
    }

    print("\n" + "=" * 78)
    print("STEP 2/11 -- resting CBF baseline: task's stated bands vs modern PET [9] and ASL [10] cross-checks")
    print("=" * 78)
    whole_pass, whole_band = in_band(CBF_WHOLE_TASK, CBF_WHOLE_TASK, BAND_TOL_REL)  # trivial/non-falsifying
    print(f"Task's stated targets (non-falsifying by construction): whole={CBF_WHOLE_TASK}, "
          f"grey={CBF_GREY_TASK}, white={CBF_WHITE_TASK} mL/100g/min")
    ito_in_grey_band, grey_band = in_band(ITO_CORTICAL_CBF, CBF_GREY_TASK, BAND_TOL_REL)
    ito_in_whole_band, _ = in_band(ITO_CORTICAL_CBF, CBF_WHOLE_TASK, BAND_TOL_REL)
    alisch_gm_in_grey_band, _ = in_band(ALISCH_GM_CBF, CBF_GREY_TASK, BAND_TOL_REL)
    alisch_wm_in_white_band, white_band = in_band(ALISCH_WM_CBF, CBF_WHITE_TASK, BAND_TOL_REL)
    print(f"Ito 2004 PET [9] cortical CBF={ITO_CORTICAL_CBF}+/-{ITO_CORTICAL_SD} mL/100g/min vs task's "
          f"grey band {grey_band}: {'PASS' if ito_in_grey_band else 'MISS (disclosed tension)'}; "
          f"vs task's WHOLE-brain band: {'lands inside' if ito_in_whole_band else 'misses'} (interesting, "
          f"a 'cortical ROI' average landing closer to a whole-brain figure than a pure-grey figure)")
    print(f"Alisch 2021 ASL [10, WebFetch-flagged] GM={ALISCH_GM_CBF}+/-{ALISCH_GM_SD} vs task's grey band "
          f"{grey_band}: {'PASS' if alisch_gm_in_grey_band else 'MISS (disclosed tension)'}")
    print(f"Alisch 2021 ASL [10, WebFetch-flagged] WM={ALISCH_WM_CBF}+/-{ALISCH_WM_SD} vs task's white band "
          f"{white_band}: {'PASS' if alisch_wm_in_white_band else 'MISS'}")
    f_grey_selfconsistent = (CBF_WHOLE_TASK - CBF_WHITE_TASK) / (CBF_GREY_TASK - CBF_WHITE_TASK)
    f_grey_plausible = bool(0.25 <= f_grey_selfconsistent <= 0.75)
    print(f"Grey-mass-fraction implied by the TASK'S OWN three numbers (f*80+(1-f)*20=50): "
          f"f_grey={f_grey_selfconsistent:.3f} vs plausible sanity band [0.25,0.75]: "
          f"{'PASS' if f_grey_plausible else 'FAIL'}")
    # cross-method reconciliation: does ANY f in [0,1] blend Alisch's GM/WM into Ito's cortical number?
    denom = ALISCH_GM_CBF - ALISCH_WM_CBF
    f_cross = (ITO_CORTICAL_CBF - ALISCH_WM_CBF) / denom if denom != 0 else float("nan")
    cross_reconciles = bool(0.0 <= f_cross <= 1.0)
    print(f"Cross-METHOD reconciliation (does a GM/WM blend of Alisch's ASL numbers reproduce Ito's PET "
          f"cortical number?): required fraction f={f_cross:.3f} -- "
          f"{'reconciles (valid fraction)' if cross_reconciles else 'DOES NOT RECONCILE (outside [0,1] -- ' 'a genuine, disclosed PET-vs-ASL absolute-quantification gap, not hidden)'}")
    report["resting_cbf_baseline"] = {
        "task_targets": {"whole": CBF_WHOLE_TASK, "grey": CBF_GREY_TASK, "white": CBF_WHITE_TASK},
        "ito_cortical_pet": ITO_CORTICAL_CBF, "ito_in_grey_band": ito_in_grey_band, "ito_in_whole_band": ito_in_whole_band,
        "alisch_gm_asl": ALISCH_GM_CBF, "alisch_wm_asl": ALISCH_WM_CBF,
        "alisch_gm_in_grey_band": alisch_gm_in_grey_band, "alisch_wm_in_white_band": alisch_wm_in_white_band,
        "f_grey_selfconsistent_task_numbers": f_grey_selfconsistent, "f_grey_plausible": f_grey_plausible,
        "f_cross_method_reconciliation": float(f_cross), "cross_method_reconciles": cross_reconciles,
    }

    print("\n" + "=" * 78)
    print("STEP 3/11 -- CLASSIC (idealized Lassen [1]) autoregulation curve: flat in-plateau, continuous, "
          "falls outside")
    print("=" * 78)
    p_sweep = np.linspace(20.0, 200.0, 361)
    f_classic = classic_autoreg_flow(p_sweep, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC)
    in_plateau_mask = (p_sweep >= LLA_CLASSIC) & (p_sweep <= UBP_CLASSIC)
    flat_dev_pct = float(np.max(np.abs(f_classic[in_plateau_mask] - CBF_WHOLE_TASK)) / CBF_WHOLE_TASK * 100)
    classic_flat_pass = bool(flat_dev_pct < 1e-6)
    print(f"Max deviation from perfectly flat WITHIN [60,150]: {flat_dev_pct:.6f}% -- "
          f"{'PASS (trivially true by construction, NON-FALSIFYING, disclosed as such)' if classic_flat_pass else 'FAIL'}")
    knee_lo_continuous = bool(np.isclose(classic_autoreg_flow(LLA_CLASSIC, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC), CBF_WHOLE_TASK))
    knee_hi_continuous = bool(np.isclose(classic_autoreg_flow(UBP_CLASSIC, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC), CBF_WHOLE_TASK))
    below_monotonic = bool(np.all(np.diff(f_classic[p_sweep < LLA_CLASSIC]) > 0))
    above_monotonic = bool(np.all(np.diff(f_classic[p_sweep > UBP_CLASSIC]) > 0))
    print(f"Continuity at knees: LLA={'PASS' if knee_lo_continuous else 'FAIL'}, UBP={'PASS' if knee_hi_continuous else 'FAIL'}")
    print(f"Falls (monotonic increasing toward the knee, i.e. drops as pressure drops) below LLA: "
          f"{'PASS' if below_monotonic else 'FAIL'}; rises above UBP: {'PASS' if above_monotonic else 'FAIL'}")
    report["classic_curve"] = {
        "flat_deviation_pct_in_plateau": flat_dev_pct, "flat_pass_nonfalsifying": classic_flat_pass,
        "knee_lo_continuous": knee_lo_continuous, "knee_hi_continuous": knee_hi_continuous,
        "below_lla_monotonic": below_monotonic, "above_ubp_monotonic": above_monotonic,
        "p_sweep_sample": {"p": [20.0, 40.0, 60.0, 93.33, 150.0, 180.0],
                            "f": [float(classic_autoreg_flow(x, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC)) for x in
                                  [20.0, 40.0, 60.0, 93.33, 150.0, 180.0]]},
    }

    print("\n" + "=" * 78)
    print("STEP 4/11 -- MODERN 'leaky-plateau' curve (Numan 2014 [4] own elasticities integrated) -- THE "
          "decisive, quantified modern-vs-classic falsifier, pre-registered to show >5% deviation from flat")
    print("=" * 78)
    f_leaky = leaky_autoreg_flow(p_sweep, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF,
                                  NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)
    dev_leaky_in_plateau_pct = np.abs(f_leaky[in_plateau_mask] - CBF_WHOLE_TASK) / CBF_WHOLE_TASK * 100
    max_leaky_dev_pct = float(np.max(dev_leaky_in_plateau_pct))
    leaky_shows_real_leak = bool(max_leaky_dev_pct > LEAK_DEVIATION_GATE_PCT)
    print(f"PRE-REGISTERED prediction: the modern/Numan-anchored curve WILL deviate from perfectly flat "
          f"by MORE than {LEAK_DEVIATION_GATE_PCT}% somewhere strictly inside the nominal 60-150 plateau.")
    print(f"Measured: max deviation-from-flat inside [60,150] = {max_leaky_dev_pct:.2f}% -- "
          f"{'CONFIRMED (modern-vs-classic tension is real and quantified, not just described)' if leaky_shows_real_leak else 'NOT confirmed'}")
    f_at_65 = float(leaky_autoreg_flow(65.0, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)[0])
    f_at_140 = float(leaky_autoreg_flow(140.0, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)[0])
    print(f"Worked example: at CPP=65mmHg (5mmHg inside the nominal lower knee), classic predicts "
          f"F={CBF_WHOLE_TASK:.1f} (flat); Numan-anchored leaky model predicts F={f_at_65:.2f} "
          f"({(f_at_65/CBF_WHOLE_TASK-1)*100:+.1f}%). At CPP=140mmHg (10mmHg inside the upper knee): "
          f"leaky model predicts F={f_at_140:.2f} ({(f_at_140/CBF_WHOLE_TASK-1)*100:+.1f}%).")
    leaky_knee_lo_continuous = bool(np.isclose(
        leaky_autoreg_flow(LLA_CLASSIC, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)[0],
        CBF_WHOLE_TASK * (LLA_CLASSIC / P_REF) ** NUMAN_SLOPE_DECREASE))
    leaky_knee_hi_continuous = bool(np.isclose(
        leaky_autoreg_flow(UBP_CLASSIC, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)[0],
        CBF_WHOLE_TASK * (UBP_CLASSIC / P_REF) ** NUMAN_SLOPE_INCREASE))
    print(f"Continuity at knees (power-law meets passive line): LLA={'PASS' if leaky_knee_lo_continuous else 'FAIL'}, "
          f"UBP={'PASS' if leaky_knee_hi_continuous else 'FAIL'}")
    report["leaky_plateau_curve"] = {
        "numan_slope_decrease": NUMAN_SLOPE_DECREASE, "numan_slope_increase": NUMAN_SLOPE_INCREASE,
        "max_deviation_pct_in_plateau": max_leaky_dev_pct, "leak_confirmed_gt_5pct": leaky_shows_real_leak,
        "f_at_cpp65": f_at_65, "f_at_cpp140": f_at_140,
        "knee_lo_continuous": leaky_knee_lo_continuous, "knee_hi_continuous": leaky_knee_hi_continuous,
    }

    print("\n" + "=" * 78)
    print("STEP 5/11 -- FORCED ADVERSARY (void-floor): 'no autoregulation' (constant R) vs the regulated model")
    print("=" * 78)
    f_void = void_floor_constant_r_flow(p_sweep, CBF_WHOLE_TASK, P_REF)
    range_classic_in_plateau = float(np.max(f_classic[in_plateau_mask]) - np.min(f_classic[in_plateau_mask]))
    range_void_in_plateau = float(np.max(f_void[in_plateau_mask]) - np.min(f_void[in_plateau_mask]))
    void_floor_margin = (range_void_in_plateau / max(range_classic_in_plateau, 1e-9))
    void_floor_pass = bool(range_void_in_plateau > 3.0 * max(range_classic_in_plateau, 1.0))
    print(f"Total CBF variation across [60,150]: WITH autoregulation={range_classic_in_plateau:.4f} "
          f"mL/100g/min (0 by construction) vs WITHOUT (void-floor, constant R)={range_void_in_plateau:.2f} "
          f"mL/100g/min -- {'PASS: void-floor swings far more (autoregulation is doing real necessary work)' if void_floor_pass else 'FAIL'}")
    print(f"At CPP=60: void-floor predicts F={f_void[np.argmin(np.abs(p_sweep-60))]:.2f} vs regulated "
          f"F={CBF_WHOLE_TASK:.1f}; at CPP=150: void-floor predicts F={f_void[np.argmin(np.abs(p_sweep-150))]:.2f} "
          f"vs regulated F={CBF_WHOLE_TASK:.1f} -- a {(f_void[np.argmin(np.abs(p_sweep-150))]/f_void[np.argmin(np.abs(p_sweep-60))]-1)*100:.0f}% "
          f"swing with NO mechanism, vs 0% with one.")
    report["void_floor_no_autoregulation"] = {
        "range_with_autoreg": range_classic_in_plateau, "range_without_autoreg": range_void_in_plateau,
        "margin_ratio": float(void_floor_margin), "pass": void_floor_pass,
    }

    print("\n" + "=" * 78)
    print(f"STEP 6/11 -- THE TASK'S EXPLICIT FALSIFIER: perturbation at MAP={MAP_PERTURBATION_TEST}mmHg "
          f"(below the classic knee) MUST drop CBF")
    print("=" * 78)
    icp_mid = (ICP_NORMAL_LOW + ICP_NORMAL_HIGH) / 2.0
    cpp_at_map50 = MAP_PERTURBATION_TEST - icp_mid
    f_classic_at_map50_direct = float(classic_autoreg_flow(MAP_PERTURBATION_TEST, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC))
    f_leaky_at_map50_direct = float(leaky_autoreg_flow(MAP_PERTURBATION_TEST, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)[0])
    f_classic_at_cpp_after_icp = float(classic_autoreg_flow(cpp_at_map50, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC))
    f_leaky_at_cpp_after_icp = float(leaky_autoreg_flow(cpp_at_map50, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC, P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE)[0])
    drop_classic_direct_pct = (1 - f_classic_at_map50_direct / CBF_WHOLE_TASK) * 100
    drop_leaky_direct_pct = (1 - f_leaky_at_map50_direct / CBF_WHOLE_TASK) * 100
    drop_classic_cpp_pct = (1 - f_classic_at_cpp_after_icp / CBF_WHOLE_TASK) * 100
    drop_leaky_cpp_pct = (1 - f_leaky_at_cpp_after_icp / CBF_WHOLE_TASK) * 100
    print(f"PRE-REGISTERED gate: CBF must drop by >= {MIN_DROP_PCT_GATE}% under BOTH curve variants, "
          f"under BOTH interpretations (literal MAP=50 treated as CPP; MAP=50 minus typical ICP={icp_mid:.1f}mmHg -> CPP={cpp_at_map50:.1f}).")
    print(f"Literal MAP=50-as-CPP: classic drop={drop_classic_direct_pct:.1f}%, leaky drop={drop_leaky_direct_pct:.1f}%")
    print(f"MAP=50 through CPP=MAP-ICP (ICP={icp_mid:.1f}mmHg) -> CPP={cpp_at_map50:.1f}mmHg: "
          f"classic drop={drop_classic_cpp_pct:.1f}%, leaky drop={drop_leaky_cpp_pct:.1f}%")
    gate_map50_classic = bool(drop_classic_direct_pct >= MIN_DROP_PCT_GATE and drop_classic_cpp_pct >= MIN_DROP_PCT_GATE)
    gate_map50_leaky = bool(drop_leaky_direct_pct >= MIN_DROP_PCT_GATE and drop_leaky_cpp_pct >= MIN_DROP_PCT_GATE)
    print(f"Gate (classic model, both interpretations): {'PASS' if gate_map50_classic else 'FAIL'}; "
          f"Gate (leaky/modern model, both interpretations): {'PASS' if gate_map50_leaky else 'FAIL'} "
          f"(robustness-to-model-choice: the qualitative conclusion -- CBF DROPS -- holds regardless of "
          f"which curve variant is used)")
    report["map50_perturbation_falsifier"] = {
        "icp_mid_mmhg": icp_mid, "cpp_after_icp_mmhg": cpp_at_map50,
        "drop_classic_direct_pct": drop_classic_direct_pct, "drop_leaky_direct_pct": drop_leaky_direct_pct,
        "drop_classic_cpp_pct": drop_classic_cpp_pct, "drop_leaky_cpp_pct": drop_leaky_cpp_pct,
        "gate_classic_pass": gate_map50_classic, "gate_leaky_pass": gate_map50_leaky,
    }

    print("\n" + "=" * 78)
    print("STEP 7/11 -- CO2 reactivity model: calibrated to Battisti-Charbonney [7] (TCD, PRIMARY), "
          "cross-checked vs Mutch [8] and the task's 3-4%/mmHg band; forced unbounded-linear adversary")
    print("=" * 78)
    c_sweep = np.linspace(15.0, 70.0, 111)
    g_real = co2_multiplier(c_sweep, PACO2_REF, BATTISTI_HYPER_SLOPE_PCT / 100.0,
                             WILLIE_HYPOCAPNIC_SLOPE_PCT / 100.0, CO2_G_FLOOR)
    ref_check = float(co2_multiplier(np.array([PACO2_REF]), PACO2_REF, BATTISTI_HYPER_SLOPE_PCT / 100.0,
                                      WILLIE_HYPOCAPNIC_SLOPE_PCT / 100.0, CO2_G_FLOOR)[0])
    ref_pass = bool(np.isclose(ref_check, 1.0))
    print(f"Reference point g(PaCO2={PACO2_REF})={ref_check:.6f} (must be exactly 1.0 by construction): "
          f"{'PASS (sanity, non-falsifying)' if ref_pass else 'FAIL'}")
    measured_slope_pct = (co2_multiplier(np.array([PACO2_REF + 10.0]), PACO2_REF, BATTISTI_HYPER_SLOPE_PCT / 100.0,
                                          WILLIE_HYPOCAPNIC_SLOPE_PCT / 100.0, CO2_G_FLOOR)[0] - 1.0) * 100.0 / 10.0
    print(f"Model's hypercapnic slope (calibration check, NOT independent evidence since it is the "
          f"anchor itself): {measured_slope_pct:.2f}%/mmHg vs Battisti-Charbonney [7]={BATTISTI_HYPER_SLOPE_PCT}: "
          f"{'PASS' if np.isclose(measured_slope_pct, BATTISTI_HYPER_SLOPE_PCT) else 'FAIL'}")
    within_mutch = bool(MUTCH_RANGE_LOW <= measured_slope_pct <= MUTCH_RANGE_HIGH)
    within_task = bool(CO2_REACT_TASK_LOW <= measured_slope_pct <= CO2_REACT_TASK_HIGH)
    print(f"INDEPENDENT cross-check vs Mutch [8]'s stated 5-11%/mmHg range: "
          f"{'PASS' if within_mutch else 'FAIL'} -- genuine external corroboration (not the calibration source).")
    print(f"vs the TASK'S OWN pre-registered 3-4%/mmHg band: "
          f"{'PASS' if within_task else 'MISS -- disclosed tension: live TCD/BOLD literature (both [7] and [8]) '                    'measures CO2 reactivity running HIGHER than the tasks stated band, not force-fit to match'}")
    print("\nFORCED ADVERSARY: the naive symmetric-UNBOUNDED-linear model (no floor, no sigmoid) -- what "
          "the 'simplest possible' model would be without the sigmoid/floor correction:")
    g_naive = co2_multiplier_naive_unbounded_adversary(c_sweep, PACO2_REF, BATTISTI_HYPER_SLOPE_PCT / 100.0)
    naive_goes_negative = bool(np.any(g_naive < 0))
    naive_min = float(np.min(g_naive))
    real_min = float(np.min(g_real))
    real_stays_positive = bool(np.all(g_real > 0))
    print(f"Over PaCO2 in [15,70]mmHg: naive-unbounded adversary min multiplier={naive_min:.3f} "
          f"({'goes NEGATIVE -- unphysical, PASS (adversary correctly falls)' if naive_goes_negative else 'stays positive -- adversary does NOT fall'}); "
          f"saturating model min multiplier={real_min:.3f} ({'stays physically bounded -- PASS' if real_stays_positive else 'FAIL'})")
    print("\nVOID-FLOOR (zero CO2 reactivity) vs real model, over the same sweep:")
    total_swing_real_pct = float((np.max(g_real) - np.min(g_real)) * 100)
    void_floor_co2_pass = bool(total_swing_real_pct > 20.0)
    print(f"Zero-reactivity adversary: 0.0% swing (by definition) vs real model: {total_swing_real_pct:.1f}% "
          f"swing across PaCO2=[15,70]: {'PASS (real reactivity beats the zero-reactivity void floor by a wide margin)' if void_floor_co2_pass else 'FAIL'}")
    report["co2_reactivity"] = {
        "reference_point_pass": ref_pass, "hypercapnic_slope_pct_per_mmhg": measured_slope_pct,
        "battisti_charbonney_slope": BATTISTI_HYPER_SLOPE_PCT, "within_mutch_range": within_mutch,
        "within_task_band": within_task, "mutch_range": [MUTCH_RANGE_LOW, MUTCH_RANGE_HIGH],
        "task_band": [CO2_REACT_TASK_LOW, CO2_REACT_TASK_HIGH],
        "naive_adversary_goes_negative": naive_goes_negative, "naive_adversary_min": naive_min,
        "real_model_stays_positive": real_stays_positive, "real_model_min": real_min,
        "void_floor_total_swing_pct": total_swing_real_pct, "void_floor_pass": void_floor_co2_pass,
    }

    print("\n" + "=" * 78)
    print("STEP 8/11 -- CPP=MAP-ICP coupling: StatPearls [12] formula, cross-layer over-determination "
          "with arterial_pressure.py's already-computed MAP")
    print("=" * 78)
    icp_sweep_normal = np.array([ICP_NORMAL_LOW, icp_mid, ICP_NORMAL_HIGH])
    cpp_from_classic_map = map_classic - icp_sweep_normal
    cpp_from_windkessel_map = map_windkessel_true_mean - icp_sweep_normal
    all_cpp = np.concatenate([cpp_from_classic_map, cpp_from_windkessel_map])
    all_in_plateau = bool(np.all((all_cpp >= LLA_CLASSIC) & (all_cpp <= UBP_CLASSIC)))
    print(f"CPP formula replication: {map_classic:.2f} - {icp_mid:.1f} = {map_classic - icp_mid:.2f}mmHg, "
          f"vs StatPearls [12] normal CPP band [{CPP_NORMAL_LOW},{CPP_NORMAL_HIGH}]: "
          f"{'inside' if CPP_NORMAL_LOW <= (map_classic-icp_mid) <= CPP_NORMAL_HIGH else 'just outside (disclosed, see honest gaps)'}")
    print(f"CPP from the computed MAP (classic {map_classic:.2f} / Windkessel-true-mean "
          f"{map_windkessel_true_mean:.2f}) across the normal ICP range [{ICP_NORMAL_LOW},{ICP_NORMAL_HIGH}]: "
          f"{cpp_from_classic_map.round(1).tolist()} / {cpp_from_windkessel_map.round(1).tolist()} mmHg")
    print(f"ALL land inside the autoregulation plateau [{LLA_CLASSIC},{UBP_CLASSIC}]: "
          f"{'PASS (cross-layer over-determination: three already-built layers -- cardiac_output -> ' 'arterial_pressure -> this cerebral-autoregulation model -- agree the resting operating point is ' 'safely inside the regulated range)' if all_in_plateau else 'FAIL'}")
    print("\nFORCED ADVERSARY: 'ICP-blind' model (uses MAP as if it WERE CPP, ignoring ICP) vs the correct "
          "CPP=MAP-ICP model, at a PATHOLOGICALLY raised ICP (StatPearls [13]'s 25mmHg threshold)")
    print("OBSERVE (first pass, reported not hidden): tested at the reference body's HIGH-headroom resting MAP "
          f"({map_classic:.2f}) + ICP={ICP_PATHOLOGICAL} -> true CPP={map_classic-ICP_PATHOLOGICAL:.2f}, which "
          f"is STILL inside the flat plateau (>{LLA_CLASSIC}) -- so both models trivially agree (F=50, "
          "error=0%). This is a real, informative finding, not a null result to discard: a high-MAP patient "
          "has enough 'plateau headroom' that ICP-blindness alone does not yet matter.")
    print("ORIENT: the adversary needs a MAP that LOOKS unremarkable/borderline-normal on its own (i.e. "
          "still >= LLA, so a clinician reading MAP alone would not be alarmed) but where the TRUE CPP, "
          "once ICP is correctly subtracted, falls below LLA -- isolating ICP-blindness as the error "
          "source, not confounding it with an already-abnormal MAP.")
    map_headroom, map_borderline = map_classic, 70.0
    cpp_headroom = map_headroom - ICP_PATHOLOGICAL
    cpp_borderline = map_borderline - ICP_PATHOLOGICAL
    f_icp_blind_headroom = float(classic_autoreg_flow(map_headroom, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC))
    f_correct_headroom = float(classic_autoreg_flow(cpp_headroom, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC))
    f_icp_blind_borderline = float(classic_autoreg_flow(map_borderline, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC))
    f_correct_borderline = float(classic_autoreg_flow(cpp_borderline, CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC))
    err_headroom_pct = abs(f_icp_blind_headroom - f_correct_headroom) / f_correct_headroom * 100
    err_borderline_pct = abs(f_icp_blind_borderline - f_correct_borderline) / f_correct_borderline * 100
    icp_blind_adversary_pass = bool(err_borderline_pct > 10.0)
    print(f"DECIDE/ACT: re-test at a borderline-but-still->=LLA MAP={map_borderline} (a BP a clinician would "
          f"NOT flag as abnormal on its own) + the SAME pathological ICP={ICP_PATHOLOGICAL}: true "
          f"CPP={cpp_borderline:.2f}mmHg (now BELOW LLA={LLA_CLASSIC}) -> correct model predicts "
          f"F={f_correct_borderline:.2f}; ICP-blind model (reads MAP={map_borderline} directly, sees it as "
          f">=LLA, predicts flat/normal) predicts F={f_icp_blind_borderline:.2f} -- error={err_borderline_pct:.1f}%: "
          f"{'PASS (ignoring ICP now produces a clinically meaningful, falsely-reassuring CBF prediction)' if icp_blind_adversary_pass else 'FAIL'}")
    print(f"CONVERGENCE: error at high-headroom MAP={map_headroom:.1f} = {err_headroom_pct:.1f}% vs error at "
          f"borderline MAP={map_borderline:.1f} = {err_borderline_pct:.1f}% -- ICP-blindness is not a fixed "
          "error, its danger DEPENDS on how much plateau headroom the baseline MAP has, a genuine geometric "
          "consequence of the clipped-resistance structure (Step 1), not an arbitrary scenario choice.")
    report["cpp_coupling"] = {
        "cpp_from_classic_map": cpp_from_classic_map.tolist(), "cpp_from_windkessel_map": cpp_from_windkessel_map.tolist(),
        "all_in_plateau": all_in_plateau,
        "icp_blind_error_pct_at_headroom_map": err_headroom_pct,
        "icp_blind_error_pct_at_borderline_map": err_borderline_pct,
        "icp_blind_adversary_pass": icp_blind_adversary_pass,
    }

    print("\n" + "=" * 78)
    print("STEP 9/11 -- cross-layer CO2 coupling: acid_base_co2.py's PaCO2 perturbation (40->50mmHg) "
          "propagated through this model")
    print("=" * 78)
    g_at_perturbed = float(co2_multiplier(np.array([paco2_perturbed]), PACO2_REF, BATTISTI_HYPER_SLOPE_PCT / 100.0,
                                           WILLIE_HYPOCAPNIC_SLOPE_PCT / 100.0, CO2_G_FLOOR)[0])
    pct_rise = (g_at_perturbed - 1.0) * 100
    coupling_positive_and_bounded = bool(0.0 < pct_rise < 100.0)
    print(f"acid_base_co2.py's acute scenario: PaCO2 {paco2_normal:.0f}->{paco2_perturbed:.0f}mmHg. "
          f"This model predicts CBF multiplier g={g_at_perturbed:.4f} ({pct_rise:+.1f}%): "
          f"{'PASS (positive, bounded, order-of-magnitude-plausible rise)' if coupling_positive_and_bounded else 'FAIL'}")
    report["acid_base_coupling"] = {
        "paco2_perturbed": paco2_perturbed, "predicted_cbf_pct_change": pct_rise,
        "coupling_positive_and_bounded": coupling_positive_and_bounded,
    }

    print("\n" + "=" * 78)
    print("STEP 10/11 -- combined joint model F(CPP,PaCO2): sanity sweep, bounded + directionally correct")
    print("=" * 78)
    scenarios = [
        ("rest", P_REF, PACO2_REF),
        ("hypotensive_hypercapnic_stress", 50.0, 55.0),
        ("hypertensive_hypocapnic_stress", 160.0, 28.0),
        ("severe_hypotension_below_lla", 40.0, PACO2_REF),
    ]
    combined_results = {}
    all_positive = True
    for name, p_val, c_val in scenarios:
        f_val = float(combined_cbf(np.array([p_val]), np.array([c_val]), CBF_WHOLE_TASK, LLA_CLASSIC, UBP_CLASSIC,
                                    P_REF, NUMAN_SLOPE_DECREASE, NUMAN_SLOPE_INCREASE, PACO2_REF,
                                    BATTISTI_HYPER_SLOPE_PCT / 100.0, WILLIE_HYPOCAPNIC_SLOPE_PCT / 100.0, CO2_G_FLOOR)[0])
        combined_results[name] = {"cpp": p_val, "paco2": c_val, "cbf": f_val}
        all_positive = all_positive and (f_val > 0)
        print(f"  {name:32s} CPP={p_val:6.1f} PaCO2={c_val:5.1f} -> CBF={f_val:6.2f} mL/100g/min "
              f"({(f_val/CBF_WHOLE_TASK-1)*100:+.1f}% vs rest)")
    combined_bounded_pass = bool(all_positive and combined_results["rest"]["cbf"] > 0)
    print(f"All scenarios positive/bounded: {'PASS' if combined_bounded_pass else 'FAIL'}; rest scenario "
          f"reproduces F0 by construction ({combined_results['rest']['cbf']:.2f} vs {CBF_WHOLE_TASK})")
    report["combined_model_scenarios"] = combined_results
    report["combined_model_bounded_pass"] = combined_bounded_pass

    print("\n" + "=" * 78)
    print("STEP 11/11 -- GATES SUMMARY")
    print("=" * 78)
    gates = {
        # NOTE: "does an independent measurement match the TASK'S OWN stated textbook band" is an
        # external-literature-agreement FINDING, not a pipeline-correctness assertion my model can pass
        # or fail (same convention arterial_pressure.py established for its own TPR-vs-measured-MAP
        # tension) -- both the grey-matter MISS and the white-matter PASS are reported in
        # open_modeling_uncertainty/resting_cbf_baseline below, deliberately excluded from gates so a
        # real disclosed tension can never be scored as if it were a bug.
        "grey_white_fraction_selfconsistent_plausible": f_grey_plausible,
        "classic_curve_flat_in_plateau_nonfalsifying": classic_flat_pass,
        "classic_curve_continuous_at_knees": bool(knee_lo_continuous and knee_hi_continuous),
        "classic_curve_falls_outside_plateau": bool(below_monotonic and above_monotonic),
        "leaky_plateau_confirms_modern_vs_classic_tension_gt5pct": leaky_shows_real_leak,
        "leaky_plateau_continuous_at_knees": bool(leaky_knee_lo_continuous and leaky_knee_hi_continuous),
        "void_floor_no_autoregulation_adversary_falls": void_floor_pass,
        "map50_perturbation_drops_cbf_classic": gate_map50_classic,
        "map50_perturbation_drops_cbf_leaky": gate_map50_leaky,
        "co2_reference_point_nonfalsifying": ref_pass,
        "co2_slope_matches_calibration_source": bool(np.isclose(measured_slope_pct, BATTISTI_HYPER_SLOPE_PCT)),
        "co2_slope_within_independent_mutch_range": within_mutch,
        "co2_naive_unbounded_adversary_falls_negative": naive_goes_negative,
        "co2_saturating_model_stays_physical": real_stays_positive,
        "co2_void_floor_zero_reactivity_beaten": void_floor_co2_pass,
        "cpp_crosslayer_overdetermination_in_plateau": all_in_plateau,
        "icp_blind_adversary_falls_at_pathological_icp": icp_blind_adversary_pass,
        "acid_base_crosslayer_coupling_positive_bounded": coupling_positive_and_bounded,
        "combined_model_bounded": combined_bounded_pass,
    }
    gates = {k: bool(v) for k, v in gates.items()}
    open_modeling_uncertainty = {
        "classic_60_150_plateau_is_the_textbook_idealization_not_the_measured_reality": (
            f"Willie et al 2014 [5] states directly, as one of its own four key theses: 'cerebral "
            f"autoregulation does not maintain constant perfusion through a mean arterial pressure range "
            f"of 60-150mmHg.' Numan et al 2014 [4] quantifies this: real regression slopes of "
            f"{NUMAN_SLOPE_DECREASE}+/-{NUMAN_SLOPE_DECREASE_SD} (%CBF/%MAP, decreasing) and "
            f"{NUMAN_SLOPE_INCREASE}+/-{NUMAN_SLOPE_INCREASE_SD} (increasing) -- a real, nonzero leak, "
            f"reproduced here as a {max_leaky_dev_pct:.1f}% deviation-from-flat inside the nominal plateau. "
            "HELD OPEN, not resolved either way -- the classic curve remains a useful, order-of-magnitude-"
            "correct textbook idealization; the leak is real but modest relative to the outside-plateau fall."
        ),
        "co2_reactivity_measured_higher_than_the_tasks_stated_band": (
            f"Battisti-Charbonney [7] (TCD, PRIMARY) measured {BATTISTI_HYPER_SLOPE_PCT}%/mmHg; Mutch [8] "
            f"(BOLD context, secondary) states 5-11%/mmHg -- both ABOVE the task's stated 3-4%/mmHg. "
            "Not force-fit to match; reported as a genuine, disclosed tension between the task's stated "
            "textbook figure and live-verified modern TCD/BOLD measurements."
        ),
        "grey_matter_cbf_measured_lower_than_the_tasks_stated_80": (
            f"Ito 2004 PET [9] cortical CBF={ITO_CORTICAL_CBF} MISSES the task's grey~80 band by a wide "
            f"margin ({grey_band[0]:.0f}-{grey_band[1]:.0f} pre-registered); Alisch 2021 ASL [10, "
            f"WebFetch-flagged] GM={ALISCH_GM_CBF} MISSES it too -- BOTH modern quantitative methods land "
            f"well below the task's stated grey~80 figure. By contrast, WHITE matter (Alisch WM="
            f"{ALISCH_WM_CBF}) lands cleanly INSIDE the task's stated white~20 band ({white_band[0]:.0f}-"
            f"{white_band[1]:.0f}) -- an asymmetric finding, reported exactly as measured, not smoothed to "
            "look consistent. Plausible Orient-level explanations (not adjudicated here): partial-volume/"
            "ROI-averaging dilution in modern cortical ROIs vs older, more localized measurement "
            "traditions; known ASL absolute-quantification undershoot vs PET/Kety-Schmidt gold standards. "
            "The two modern methods (PET vs ASL) do not even fully reconcile with EACH OTHER (Step 2's "
            "cross-method check: no valid grey-mass-fraction in [0,1] blends Alisch's GM/WM numbers "
            f"into Ito's cortical number; required f={f_cross:.2f})."
        ),
        "static_vs_dynamic_autoregulation_not_modeled_as_separate": (
            "This entire model is STATIC autoregulation (steady-state F(P) relationship) -- Numan [4] is "
            "itself explicitly a static-autoregulation reanalysis. Rangel-Castilla et al 2008 [15] states "
            "'cerebral autoregulation is more a concept than a physically measurable entity' and that "
            "'assessment of cerebral autoregulation is best achieved with DYNAMIC autoregulation methods' "
            "(transient response to a sudden BP change, e.g. thigh-cuff release, Tiecks ARI 0-9). No "
            "dynamic/transient-response model exists in this script -- held explicitly OPEN, per the task's "
            "own instruction."
        ),
        "hypocapnic_co2_slope_and_alisch_absolute_values_are_webfetch_summary_tier": (
            "Two numbers in this model (the 1-3%/mmHg hypocapnic range from [5]; the exact GM=33.0/WM=18.0 "
            "decimal values from [10]) were obtained via a WebFetch summarization pass over PMC full text, "
            "not independently raw-machine-verified -- direct curl access to PMC full-text pages was "
            "reCAPTCHA-blocked (confirmed: a direct curl fetch returned a Google reCAPTCHA "
            "challenge page, not article content). Both papers' own PubMed/EuropePMC abstracts ARE "
            "independently, live, primary-verified (the qualitative findings, e.g. GM>WM CBF, hold at full "
            "confidence); only the specific decimal figures carry this disclosed lower-confidence flag."
        ),
        "multiplicative_pressure_co2_independence_is_a_disclosed_simplification": (
            "This model combines pressure-autoregulation and CO2-reactivity multiplicatively (F(P,C)="
            "F_autoreg(P)*g(C)), treating them as independently adjusting a shared resistance. Real "
            "cerebrovascular physiology has documented pressure-CO2 INTERACTION effects (e.g. hypercapnic "
            "vasodilation can consume some of the vessel's dilatory reserve, narrowing the effective "
            "pressure-autoregulation range) that this simplified structural form does not capture -- no "
            "specific interaction citation was live-verified for a quantified magnitude, so "
            "this is disclosed as an open structural gap, not asserted with an uncited mechanism."
        ),
        "cpp_normal_band_boundary_case": (
            f"MAP-ICP at the classic MAP ({map_classic:.2f}) and mid-normal ICP ({icp_mid:.1f}) = "
            f"{map_classic - icp_mid:.2f}mmHg, which sits slightly ABOVE StatPearls [12]'s stated normal "
            "CPP upper bound of 80mmHg -- a small, disclosed boundary mismatch (not forced to fit), "
            "consistent with CPP varying with the specific MAP/ICP operating point chosen, not a fixed "
            "constant."
        ),
    }
    overall_pass = all(gates.values())
    print(json.dumps(gates, indent=2))
    print(f"\n{sum(gates.values())}/{len(gates)} gates PASS")
    print("\nOPEN MODELING UNCERTAINTY (reported, does NOT gate overall_pass):")
    print(json.dumps(open_modeling_uncertainty, indent=2))
    print(f"\nOVERALL: {'PASS' if overall_pass else 'FAIL/SURPRISE -- see gates above'}")

    report["citations_verified_live"] = {
        "lassen_1959_pmid": "13645234", "kety_schmidt_1948_n2o_pmid": "16695568",
        "kety_schmidt_1948_co2_pmid": "16695569", "numan_2014_pmid": "25205587",
        "willie_2014_review_pmid": "24396059", "willie_2011_tcd_pmid": "21276818",
        "battisti_charbonney_2011_pmid": "21521758", "mutch_2012_pmid": "23139743",
        "ito_2004_pet_pmid": "14730405", "alisch_2021_asl_pmid": "33596183",
        "czosnyka_pickard_2004_pmid": "15145991", "statpearls_cerebral_perfusion_pressure_pmid": "30725956",
        "statpearls_increased_icp_pmid": "29489250", "statpearls_cerebral_autoregulation_pmid": "31985976",
        "rangel_castilla_2008_pmid": "18828705",
    }
    report["gates"] = gates
    report["open_modeling_uncertainty"] = open_modeling_uncertainty
    report["overall_pass"] = overall_pass
    # ADDED  (HOLE-CROSSSCALE-FOUR-BODY-MASS-INCONSISTENCY-PROPAGATION- fix):
    # this cell does NOT consume an absolute body or organ mass anywhere -- every flow number is
    # per-100g-of-tissue NORMALIZED (mL/100g/min), so brain mass cancels out algebraically and never
    # enters an arithmetic term. Declared explicitly (applicable=False) rather than silently omitted,
    # so the linter can tell "verified not applicable" apart from "forgot to declare".
    report["reference_body"] = {
        "applicable": False,
        "note": "all CBF figures (task targets, Ito PET, Alisch ASL) are per-100g-tissue "
               "normalized -- no absolute brain_mass_g term appears in this cell's arithmetic "
               "(contrast organ_o2_consumption_partition.py, which multiplies by a literal "
               "brain_mass_g=1350.0 to get an absolute flow)",
    }

    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {OUT_PATH}")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
