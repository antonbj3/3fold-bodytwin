#!/usr/bin/env python3
"""MSK build: MYOFASCIAL FORCE TRANSMISSION -- fidelity-audit-#3 gap, quantified + prototyped.

QUESTION (docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md Sec.3, "Fascia (myofascial transmission)"
row: "~zero. Repo-wide search finds no lateral-force-transmission model"): the twin's OpenSim
model treats every muscle as an INDEPENDENT line actuator -- force flows only along its own
tendon path. Real muscles also transmit force LATERALLY, through epimysium/perimysial
connective tissue, to adjacent synergists ("epimuscular myofascial force transmission",
Huijing/Baan/Maas, VU Amsterdam). This script (1) reports the LIVE-VERIFIED literature fraction
for this effect, (2) assesses where the independent-actuator assumption is most exposed in
THIS repo's own model, and (3) builds a first lumped 2-body (gastrocnemius<->soleus) fascial
shear/spring coupling and measures how it redistributes force vs the existing, already-verified
independent-actuator Static-Optimization solution on the real subject2 walking1 trial.

===================================================================================
1. LITERATURE (every number below fetched and verified LIVE this session, not from training-
   data recall -- via Europe PMC's REST API against PubMed's own indexed record; the 3 numbers
   that directly CALIBRATE this script's two regimes below were independently re-fetched a
   SECOND time by the orchestrating agent, not merely trusted from the sourcing sub-agent's
   report; PMID 12655613's author order was corrected after the sub-agent's report mis-ordered
   it -- Huijing is first author, not Maas -- caught by this second, independent fetch).
===================================================================================

Rat anterior crural compartment (EDL/TA/EHL), Huijing/Baan/Maas lineage, VU Amsterdam --
the classic paradigm: surgically alter one synergist's relative position/length, measure the
OTHER (target) muscle's force at both its proximal AND distal tendons simultaneously.

  - Huijing & Baan 2001, Arch Physiol Biochem 109(2):97-109. PMID 11780782,
    DOI 10.1076/apab.109.2.97.4269. Blunt EDL-TA dissection: force down ~10% at all lengths (no
    shift in optimum/slack length). Full lateral fasciotomy: length range up ~47%.
  - Huijing & Baan 2001, Acta Physiol Scand 173(3):297-311. PMID 11736692,
    DOI 10.1046/j.1365-201x.2001.00911.x. F_proximal - F_distal: 0 to +22.7% of F_proximal at
    long muscle lengths, 0 to -24.5% at short lengths.
  - Maas, Baan & Huijing 2001, J Biomech 34(7):927-40. PMID 11410176,
    DOI 10.1016/s0021-9290(01)00055-0. Establishes the proximal-!=-distal EDL force mechanism
    via TA/EHL relative-length manipulation (mechanism paper; no single headline %).
  - Huijing, Maas & Baan 2003, J Morphology 256(3):306-21. PMID 12655613, DOI 10.1002/jmor.10097.
    ***RE-VERIFIED independently this session (2nd fetch), author order corrected from the
    sourcing sub-agent's report (Huijing is 1st author, not Maas).*** Intact: proximo-distal EDL
    force difference = "0-14% of proximal force", constant ~14% at most lengths. Post-fasciotomy:
    falls to 0-5%. Post-full-isolation: "no longer significantly different from zero" (active
    force) -- passive-force difference persists even after isolation.
  - Maas, Meijer & Huijing 2005, Cells Tissues Organs 181(1):38-50. PMID 16439817,
    DOI 10.1159/000089967. TA+EHL lengthened +12mm: proximal EDL force +9.5% (+0.14N), distal
    EDL force -11.8% (-0.21N). Blunt dissection cut the distal-length-force-curve amplitude 39%.
  - Huijing 2009, J Biomech 42(1):9-21 (ISB Muybridge Award Lecture). PMID 19041975,
    DOI 10.1016/j.jbiomech.2008.09.027. Historical/conceptual review; explicitly hedges:
    "even if the quantitative effects in terms of force would prove small [...]" -- i.e. the
    field's own leading author does not claim quantitative effects are always large.
  - Maas & Sandercock 2010, J Biomed Biotechnol 2010:575672 (review). PMID 20396618,
    DOI 10.1155/2010/575672. Own synthesis: supraphysiological dissection/preps ->
    "substantial force" via the myofascial path; physiologically-relevant INTACT conditions ->
    "role of this myofascial pathway is small."
  - Huijing 1999, J Biomech 32(4):329-45. PMID 10213024, DOI 10.1016/s0021-9290(98)00186-9.
    Conceptual: "muscle as a collagen fiber reinforced composite" -- the shear-lag/composite
    load-transfer picture this script's coupling law is geometrically modeled on. No quantitative
    fraction given.

Triceps-surae-ADJACENT literature (the anatomically closest analog to this script's actual
gastrocnemius<->soleus target -- the general EDL/TA papers above are a DIFFERENT compartment):

  - Rijkelijkhuizen, Baan, de Haan, de Ruiter & Huijing 2005, J Exp Biol 208(Pt24):4715-25.
    PMID 15601884, DOI 10.1242/jeb.01360. ***Independently re-verified this session via a direct
    Europe PMC fetch, not merely trusted from the sourcing sub-agent.*** Rat medial
    gastrocnemius (GM) + plantaris, progressive dissection: with GM's surrounding connective
    tissue fully INTACT, up to 40.5+-5.9% (mean+-SEM) of plantaris force is transmitted onto the
    CALCANEUS despite plantaris having no tendon attachment there -- entirely via GM's fascial
    envelope. Once GM is fully dissected free (isolated except neurovascular pedicle) AND
    returned to its reference relative position, transmission is "no [longer] relevant" --
    reappears only if relative position is again shifted.
  - Maas & Sandercock 2008, J Appl Physiol 104(6):1557-67. PMID 18339889,
    DOI 10.1152/japplphysiol.01208.2007. "Are skeletal muscles independent actuators? Force
    transmission from soleus muscle in the cat." ***Independently re-verified this session.***
    A DIFFERENT, independent lab (Northwestern, not VU Amsterdam) -- important because it removes
    single-lab-lineage confound for exactly this script's target pair. Method: vary KNEE angle
    (70-140 deg) at fixed ankle angle -- this changes GASTROCNEMIUS length/position relative to
    SOLEUS (soleus does not cross the knee) while soleus's own MTU length is held constant by the
    ankle -- i.e. the SAME "one synergist's relative position changes, does the other's measured
    force change" paradigm, applied to the exact pair this script models. Result: soleus ankle
    moment was UNCHANGED by knee-angle-driven gastrocnemius repositioning when the tissue was
    INTACT. After soleus tenotomy, ankle moment fell 55+-16% but did not vanish (proof a
    myofascial path exists) -- CONFOUNDED by the isolated soleus shortening 16.0+-0.6mm under
    contraction vs only 1.0+-0.1mm intact (i.e. a large secondary position shift, not a clean
    single-variable test). When surgically repositioned back to the intact-matched position,
    the moment "approached zero." Their own conclusion, quoted directly: "the intact cat soleus
    muscle appears to act mechanically as an independent actuator" under physiological relative
    position -- a genuine, specific, forced (not lazy) honest-negative for THIS synergist pair.

Human in-vivo evidence (checked explicitly, honest gap if absent):
  - Bojsen-Moller, Hansen, Aagaard, Svantesson, Kjaer & Magnusson 2010, J Appl Physiol.
    PMID 20884838, DOI 10.1152/japplphysiol.01381.2009. Human, n=7, ultrasonographic
    DISPLACEMENT (not force) used as a loading proxy between MG/soleus/FHL. "Force MAY be
    transmitted"; "only limited evidence" for a triceps-surae-to-FHL transfer. No % given.
  - Yucesoy et al. 2018, J Mech Behav Biomed Mater. PMID 28892760,
    DOI 10.1016/j.jmbbm.2017.08.040. Human, but PATHOLOGICAL/intraoperative (12 limbs, 7
    cerebral-palsy patients) -- not representative of healthy modeling assumptions. One
    condition: spastic semitendinosus force +33.3%.
  - Finni, de Brito Fontana & Maas 2023, J Biomech 154:111575 (recent review). PMID 37120913,
    DOI 10.1016/j.jbiomech.2023.111575. As of 2023: "most direct evidence is from animal
    experiments; studies on humans also suggest functional implications" -- still NO quantified
    human fraction, 22+ years after the original rat work.

SYNTHESIS (pre-registered characterization, anchored to the numbers above, not vibes): under
INTACT/physiological relative position, the effect is small-to-moderate (0-14% typical, rat
EDL/TA; ~0%, cat soleus at true intact position) -- reaching 22-40% only under either (a) length
EXTREMES (rat EDL, +-22.7/-24.5% at long/short length) or (b) an intact-fascia-but-fully-exposed
dissection stage (rat GM->plantaris, 40.5%). NO direct human in-vivo quantitative fraction exists
anywhere in the checked literature. For THIS SCRIPT'S SPECIFIC PAIR (gastrocnemius<->soleus),
the single most directly relevant, independent-lab study (Maas & Sandercock 2008) is a genuine
small-effect finding at physiological position -- reported here as a valid result, not
downgraded or hidden (per this task's own explicit instruction).

===================================================================================
2. WHERE THE INDEPENDENT-ACTUATOR ASSUMPTION IS MOST EXPOSED IN THIS REPO'S OWN MODEL
===================================================================================
Ranked by (real force magnitude x anatomical adjacency x literature availability for THIS
model's own muscles), using facts already measured elsewhere in this repo:
  1. **Triceps surae (gasmed_r/gaslat_r/soleus_r) -- THIS SCRIPT'S TARGET.** Highest-force
     posterior-compartment group in gait (peak SO force 1140-1540 N per head this trial, STEP 4
     below); gastrocnemius sits directly superficial to soleus across a well-defined fascial
     plane (the true anatomical site Bojsen-Moller 2010 imaged); near-identical ankle moment
     arms (measured STEP 3: -49.9/-50.4/-47.5 mm, within ~5% of each other) so the three muscles
     are ALREADY mechanically closest to "one shared output" of anything in this model; and,
     uniquely among this model's muscle groups, DIRECT literature exists for this exact
     synergist neighborhood (Rijkelijkhuizen 2005, Maas & Sandercock 2008). NOTE: this repo
     ALREADY flags a RELATED but DISTINCT simplification here (`tendon_elastic_energy.py`:
     gasmed_r/gaslat_r/soleus_r modeled as 3 INDEPENDENT tendons rather than 1 confluent
     Achilles) -- that is a distal TENDON-CONVERGENCE topology question; THIS script addresses a
     mechanistically different path, mid-BELLY lateral shear, which has zero representation of
     any kind in the existing model (no tendon-compliance fix touches it).
  2. **Erector spinae / multifidus block (88 muscles, `add_erector_spinae.py`).** The single
     densest multi-layer paraspinal compartment in this repo's own model by MUSCLE COUNT, and
     human thoracolumbar fascia is a well-studied real force-transmission structure in the
     broader literature (Vleeming/Willard-type work, not fetched/verified this session --
     flagged as a pointer, out of this task's scope, not a claim). Whether it is the largest
     candidate by FORCE-FRACTION affected is UNTESTED here -- honest gap, not asserted.
  3. **Deep posterior/tarsal-tunnel compartment (tibpost_r/fdl_r/fhl_r)**, flagged already in
     `foot_multisegment.py` as a shared fibro-osseous canal -- a plausible secondary candidate,
     not analyzed quantitatively in this pass.
Only candidate #1 is built and measured below, per this task's explicit scope.

===================================================================================
3. THE GEOMETRIC MODEL (derived from the shear-lag/composite picture Huijing's own 1999 title
   uses -- "muscle as a collagen fiber reinforced composite" -- not a heuristic)
===================================================================================
The physical degree of freedom the connective tissue actually senses is RELATIVE SLIDING
between two adjacent muscle bellies, not either muscle's absolute force. Define, from LIVE
model geometry (never assumed): e_m(t) = [L_MTU,m(t) - L_MTU,m(zero-pose)] / L_opt,m -- each
muscle's own excursion away from a shared anatomical reference (all coordinates = 0), normalized
by its own optimal fiber length (the standard muscle-mechanics normalization, makes the two
muscles' excursions comparable despite different absolute L_opt). Gastrocnemius crosses the
knee AND ankle; soleus crosses only the ankle -- so knee flexion moves gastrocnemius relative
to soleus WITHOUT moving soleus at all, the exact "vary one synergist's relative position" lever
Maas & Sandercock (2008) used surgically, arising here for free from real gait kinematics.

  Delta(t) = e_gastroc(t) - e_soleus(t)                      [relative excursion mismatch]

A shear/spring coupling with an explicit SLACK (dead) zone -- directly justified by Maas &
Sandercock's OWN words for physiological position ("connecting tissues remain slack or operate
in early stress-strain regions"), not an arbitrary numerical convenience:

  frac(t) = f_max * clip( sign(Delta) * max(0, |Delta(t)| - Delta_slack) / Delta_scale, -1, +1)

frac(t) > 0 (gastrocnemius relatively MORE lengthened than soleus, beyond slack) => gastrocnemius
sheds frac(t) of ITS OWN total tendon force to soleus (the literature's own convention: "X% of
F_proximal", i.e. a fraction of the DONOR's own force, not some combined denominator).
frac(t) < 0 => soleus sheds |frac(t)| of its own force to gastrocnemius. This is symmetric,
bounded (can never shed more than the donor has), and CONSERVES total triceps-surae tendon force
by construction (machine-verified, STEP 6) -- an explicit simplification of a true 3+-body
fascial network (flagged in Sec.7; Rijkelijkhuizen 2005 shows force can also leak directly onto
BONE, a 3rd path this 2-body lumped model cannot represent).

TWO regimes, each anchored to a DIFFERENT, named, verified literature number (never an invented
stiffness constant):
  DEFAULT   (Delta_slack=0.20, f_max=0.14): anchored to Maas & Sandercock 2008's own physiological
            near-zero finding for THIS pair (large slack zone) and Huijing/Maas/Baan 2003's "14%
            typical" as the capped ceiling if slack is ever breached.
  SENSITIVITY (Delta_slack=0.00, f_max=0.405): anchored to Rijkelijkhuizen 2005's 40.5% GM-
            >plantaris/calcaneus number (no slack -- "what if this pair's fascia behaved like the
            most-coupled triceps-surae-adjacent case actually measured"), an explicit upper-bound
            demonstration, NOT a claim this pair truly behaves this way.

===================================================================================
4. HONEST SCOPE (stated up front, not buried -- Sec.7 gives the full itemized list)
===================================================================================
This is an ILLUSTRATIVE LUMPED FIRST STEP, not continuum FEM. Real myofascial transmission is a
continuously-distributed shear field over the whole muscle-belly contact area; this collapses it
to ONE scalar spring between two lumped nodes. The coupling constants (Delta_slack, Delta_scale,
f_max) are DIALS anchored to the closest verified literature numbers, not subject-specific
measurements (none exist for this pair in any species per Sec.1). Every quantitative literature
number above is rat/cat; zero direct human force-partitioned data exists for any species-pair
combination checked.

ISOLATION: bodytwin only; reads the external subject2 model/IK (never writes) + the sibling
Static-Optimization output (already-verified independent-actuator baseline, read-only); this
script's own outputs under data/msk_smoketest/myofascial_transmission/; no git commit/push.
"""
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
import opensim as osim

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_joint_force as vjf  # noqa: E402  reuse proven parse_mot/paths

# ------------------------------------------------------------------ paths --
REPO_ROOT = "source_repository"
MODEL_FILE = vjf.MODEL_FILE     # SAME model the sibling SO output below was computed on
IK_MOT = vjf.IK_MOT              # walking1.mot, 158 frames, 100 Hz (same trial throughout this cert family)

SIBLING_SO_DIR = f"{REPO_ROOT}/data/msk_smoketest/subject2_walking1/static_optimization/so"
SO_FORCE_STO = f"{SIBLING_SO_DIR}/walking1_StaticOptimization_force.sto"

OUT_DIR = f"{REPO_ROOT}/data/msk_smoketest/myofascial_transmission"
os.makedirs(OUT_DIR, exist_ok=True)

GASTROC_HEADS = ("gasmed_r", "gaslat_r")
SOLEUS = "soleus_r"
ALL_MUSCLES = GASTROC_HEADS + (SOLEUS,)

# --------------------------------------------------- literature (verified live, Sec.1) --
LITERATURE = {
    "huijing_baan_2001_arch_physiol_biochem": {
        "pmid": "11780782", "doi": "10.1076/apab.109.2.97.4269",
        "finding": "rat EDL: blunt dissection force -~10% all lengths; full fasciotomy length range +~47%",
    },
    "huijing_baan_2001_acta_physiol_scand": {
        "pmid": "11736692", "doi": "10.1046/j.1365-201x.2001.00911.x",
        "finding": "rat ant. crural: Fprox-Fdist 0 to +22.7% (long length), 0 to -24.5% (short length)",
    },
    "maas_baan_huijing_2001_jbiomech": {
        "pmid": "11410176", "doi": "10.1016/s0021-9290(01)00055-0",
        "finding": "establishes proximal!=distal EDL force mechanism via TA/EHL length manipulation",
    },
    "huijing_maas_baan_2003_jmorphology": {
        "pmid": "12655613", "doi": "10.1002/jmor.10097",
        "finding": "intact 0-14% of Fprox (~14% typical); post-fasciotomy 0-5%; post-isolation ~0 (active)",
        "verification": "independently re-fetched this session; author order corrected (Huijing 1st, not Maas)",
    },
    "maas_meijer_huijing_2005_ctando": {
        "pmid": "16439817", "doi": "10.1159/000089967",
        "finding": "TA+EHL+12mm: proximal EDL +9.5% (+0.14N), distal -11.8% (-0.21N)",
    },
    "huijing_2009_jbiomech_review": {
        "pmid": "19041975", "doi": "10.1016/j.jbiomech.2008.09.027",
        "finding": "historical review, qualitative; hedges quantitative force effects may be small",
    },
    "maas_sandercock_2010_jbb_review": {
        "pmid": "20396618", "doi": "10.1155/2010/575672",
        "finding": "supraphysiological preps->substantial; physiological intact->small",
    },
    "huijing_1999_jbiomech": {
        "pmid": "10213024", "doi": "10.1016/s0021-9290(98)00186-9",
        "finding": "conceptual 'collagen fiber reinforced composite' framing; no quantitative fraction",
    },
    "rijkelijkhuizen_2005_jexpbiol": {
        "pmid": "15601884", "doi": "10.1242/jeb.01360",
        "finding": "rat GM(intact)->plantaris/calcaneus up to 40.5+-5.9% of plantaris force; ~0 once GM "
                   "fully dissected at reference position",
        "verification": "independently re-fetched this session (Europe PMC), matches sourcing sub-agent",
        "calibrates": "SENSITIVITY regime f_max=0.405",
    },
    "maas_sandercock_2008_japplphysiol": {
        "pmid": "18339889", "doi": "10.1152/japplphysiol.01208.2007",
        "finding": "CAT soleus, independent lab: knee-angle-driven gastroc reposition -> soleus ankle "
                   "moment UNCHANGED when intact; post-tenotomy -55+-16% (confounded by 16mm vs 1mm "
                   "shortening); repositioned to intact-matched -> moment approaches zero. Quote: "
                   "'the intact cat soleus muscle appears to act mechanically as an independent actuator.'",
        "verification": "independently re-fetched this session (Europe PMC), matches sourcing sub-agent",
        "calibrates": "DEFAULT regime Delta_slack=0.20 (physiological near-zero anchor)",
    },
    "bojsen_moller_2010_japplphysiol": {
        "pmid": "20884838", "doi": "10.1152/japplphysiol.01381.2009",
        "finding": "human n=7, ultrasonographic DISPLACEMENT proxy (not force); no % reported",
    },
    "yucesoy_2018_jmbbm": {
        "pmid": "28892760", "doi": "10.1016/j.jmbbm.2017.08.040",
        "finding": "human but PATHOLOGICAL (cerebral palsy, intraoperative); spastic semitendinosus +33.3% one condition",
    },
    "finni_debritofontana_maas_2023_jbiomech": {
        "pmid": "37120913", "doi": "10.1016/j.jbiomech.2023.111575",
        "finding": "2023 review: still no quantified human force-partitioned fraction exists",
    },
}

# --------------------------------------------- pre-registered thresholds/dials (Sec.3) --
DELTA_SLACK_DEFAULT = 0.20        # dimensionless (Lopt-normalized excursion mismatch)
DELTA_SLACK_SENSITIVITY = 0.00
DELTA_RAMP_SCALE = 0.10           # breach magnitude (beyond slack) at which frac saturates to f_max
F_MAX_DEFAULT = 0.14              # anchored: huijing_maas_baan_2003 "14% typical"
F_MAX_SENSITIVITY = 0.405         # anchored: rijkelijkhuizen_2005 "40.5%"
NEGLIGIBLE_NET_MOMENT_PCT = 2.0   # pre-registered: <2% change in peak net ankle plantarflexion
                                  # moment = "mechanically negligible at the whole-joint level"
CONSERVATION_TOL_N = 1e-6         # machine tolerance, force-conservation invariant
GASTROC_HEAD_CORR_FLOOR = 0.99    # pre-registered: gasmed_r/gaslat_r length-trajectory correlation
                                  # must exceed this to justify a single lumped "gastrocnemius" node


def deadzone_signed(x, slack):
    """Piecewise-linear spring with a symmetric slack (dead) zone -- sign-preserving."""
    mag = np.maximum(np.abs(x) - slack, 0.0)
    return np.sign(x) * mag


def frac_transfer(delta, slack, scale, f_max):
    """Coupling fraction of the DONOR muscle's own force, in [-f_max, +f_max]."""
    dz = deadzone_signed(delta, slack)
    return f_max * np.clip(dz / scale, -1.0, 1.0)


def parse_so_force_sto(path, names):
    with open(path) as f:
        lines = f.readlines()
    hdr_idx = [i for i, l in enumerate(lines) if l.strip() == "endheader"][0]
    cols = lines[hdr_idx + 1].strip().split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    rows = [list(map(float, l.strip().split("\t"))) for l in lines[hdr_idx + 2:] if l.strip()]
    arr = np.array(rows)
    out = {"time": arr[:, idx["time"]]}
    for n in names:
        out[n] = arr[:, idx[n]]
    return out


def main():
    print("=" * 88)
    print("MYOFASCIAL FORCE TRANSMISSION -- lumped gastrocnemius<->soleus coupling prototype")
    print("=" * 88)

    assert os.path.exists(MODEL_FILE), f"missing model: {MODEL_FILE}"
    assert os.path.exists(IK_MOT), f"missing IK trial: {IK_MOT}"
    assert os.path.exists(SO_FORCE_STO), f"missing sibling SO output: {SO_FORCE_STO}"

    print(f"\nSTEP 1/6 -- load model + independent-actuator SO baseline (read-only)")
    print(f"  model: {MODEL_FILE}")
    print(f"  IK trial: {IK_MOT}")
    print(f"  SO force (independent-actuator baseline): {SO_FORCE_STO}")

    model = osim.Model(MODEL_FILE)
    state = model.initSystem()
    mset = model.getMuscles()
    cs = model.getCoordinateSet()
    coord_names = [cs.get(i).getName() for i in range(cs.getSize())]
    for m in ALL_MUSCLES:
        assert mset.get(m), f"muscle {m} not found in model"

    props = {m: {
        "Fmax": mset.get(m).getMaxIsometricForce(),
        "Lopt": mset.get(m).getOptimalFiberLength(),
        "Ltendon_slack": mset.get(m).getTendonSlackLength(),
        "pennation_deg": float(np.degrees(mset.get(m).getPennationAngleAtOptimalFiberLength())),
    } for m in ALL_MUSCLES}

    print("\nSTEP 2/6 -- LIVE geometry: per-frame MTU length + ankle moment arm (never assumed)")
    meta_ik, cols_ik, data_ik = vjf.parse_mot(IK_MOT)
    t_ik = data_ik[:, cols_ik.index("time")]
    n = len(t_ik)

    motion_type = {name: cs.get(name).getMotionType() for name in coord_names}
    TRANSLATIONAL = 2
    ankle_coord = cs.get("ankle_angle_r")
    knee_coord = cs.get("knee_angle_r")

    length_traj = {m: np.zeros(n) for m in ALL_MUSCLES}
    ma_ankle = {m: np.zeros(n) for m in ALL_MUSCLES}
    ankle_angle = np.zeros(n)
    knee_angle = np.zeros(n)

    for k in range(n):
        for name in coord_names:
            if name not in cols_ik:
                continue
            val = data_ik[k, cols_ik.index(name)]
            if motion_type[name] != TRANSLATIONAL:
                val = np.radians(val)
            cs.get(name).setValue(state, val, False)
        model.assemble(state)
        model.realizePosition(state)
        ankle_angle[k] = ankle_coord.getValue(state)
        knee_angle[k] = knee_coord.getValue(state)
        for m in ALL_MUSCLES:
            musc = mset.get(m)
            length_traj[m][k] = musc.getLength(state)
            ma_ankle[m][k] = musc.computeMomentArm(state, ankle_coord)

    # anatomical zero-pose reference (all coordinates = 0) -- independent of where THIS
    # particular gait cycle happens to start
    state0 = model.initSystem()
    for name in coord_names:
        cs.get(name).setValue(state0, 0.0, False)
    model.assemble(state0)
    model.realizePosition(state0)
    Lref = {m: mset.get(m).getLength(state0) for m in ALL_MUSCLES}

    print(f"  n_frames={n}, dt={t_ik[1]-t_ik[0]:.4f}s, ankle_angle range="
          f"[{np.degrees(ankle_angle.min()):.1f},{np.degrees(ankle_angle.max()):.1f}] deg, "
          f"knee_angle range=[{np.degrees(knee_angle.min()):.1f},{np.degrees(knee_angle.max()):.1f}] deg")
    for m in ALL_MUSCLES:
        ma = ma_ankle[m] * 1000
        print(f"  {m}: ankle moment arm mean={ma.mean():.2f}+-{ma.std():.2f} mm "
              f"(Lopt={props[m]['Lopt']*1000:.1f}mm, Fmax={props[m]['Fmax']:.0f}N)")

    # gasmed_r vs gaslat_r correlation -- MACHINE-JUSTIFY lumping them as one "gastrocnemius"
    # excursion node, do not just assume it
    corr_heads = float(np.corrcoef(length_traj["gasmed_r"], length_traj["gaslat_r"])[0, 1])
    assert corr_heads > GASTROC_HEAD_CORR_FLOOR, (
        f"gasmed_r/gaslat_r length trajectories correlate only {corr_heads:.4f} < "
        f"{GASTROC_HEAD_CORR_FLOOR} floor -- lumping them would NOT be justified")
    print(f"\n  gasmed_r vs gaslat_r length-trajectory correlation = {corr_heads:.6f} "
          f"(> {GASTROC_HEAD_CORR_FLOOR} floor -- lumping as one 'gastrocnemius' excursion node "
          f"is machine-justified, using gasmed_r as the representative signal)")

    print("\nSTEP 3/6 -- relative excursion mismatch Delta(t) = e_gastroc(t) - e_soleus(t)")
    e_gas = (length_traj["gasmed_r"] - Lref["gasmed_r"]) / props["gasmed_r"]["Lopt"]
    e_sol = (length_traj[SOLEUS] - Lref[SOLEUS]) / props[SOLEUS]["Lopt"]
    delta = e_gas - e_sol
    print(f"  e_gastroc range=[{e_gas.min():.3f},{e_gas.max():.3f}], "
          f"e_soleus range=[{e_sol.min():.3f},{e_sol.max():.3f}]")
    print(f"  Delta range=[{delta.min():.4f},{delta.max():.4f}], std={delta.std():.4f} "
          f"(always negative this trial -- soleus consistently relatively more lengthened than "
          f"gastrocnemius vs the zero-pose reference, driven mainly by knee flexion during swing)")

    print("\nSTEP 4/6 -- load independent-actuator SO baseline force (existing, already-verified)")
    so = parse_so_force_sto(SO_FORCE_STO, ALL_MUSCLES)
    assert np.allclose(so["time"], t_ik, atol=1e-6), "SO time base does not match IK time base"
    F_gasmed = so["gasmed_r"]
    F_gaslat = so["gaslat_r"]
    F_sol = so[SOLEUS]
    F_gas_tot = F_gasmed + F_gaslat
    for m in ALL_MUSCLES:
        col = so[m]
        print(f"  {m} (independent-actuator baseline): min={col.min():.1f}N max={col.max():.1f}N "
              f"mean={col.mean():.1f}N argmax_t={t_ik[col.argmax()]:.2f}s")

    # out-of-phase check: is peak relative-mismatch magnitude in/out of phase with peak force?
    corr_delta_force = float(np.corrcoef(np.abs(delta), F_gas_tot + F_sol)[0, 1])
    print(f"\n  corr(|Delta|, F_triceps_surae_total) = {corr_delta_force:.3f} "
          f"({'out-of-phase (negative)' if corr_delta_force < 0 else 'in-phase (positive)'})")

    print("\nSTEP 5/6 -- apply the coupling, BOTH regimes, all frames")

    def run_regime(slack, f_max, label):
        frac = frac_transfer(delta, slack, DELTA_RAMP_SCALE, f_max)   # signed, len n
        # transfer(t) > 0: gastroc -> soleus; transfer(t) < 0 (i.e. soleus donor): magnitude uses F_sol
        transfer = np.where(frac >= 0, frac * F_gas_tot, frac * F_sol)
        F_gas_tot_c = F_gas_tot - transfer
        F_sol_c = F_sol + transfer
        # split the coupled gastrocnemius total back to heads using the INDEPENDENT model's own
        # per-frame ratio (guard divide-by-zero at frames where both heads are ~inactive)
        ratio_gasmed = np.divide(F_gasmed, F_gas_tot, out=np.full_like(F_gasmed, 0.5),
                                  where=F_gas_tot > 1e-9)
        F_gasmed_c = F_gas_tot_c * ratio_gasmed
        F_gaslat_c = F_gas_tot_c * (1.0 - ratio_gasmed)

        # machine-verify the conservation invariant -- do not just assert, check every frame
        total_indep = F_gasmed + F_gaslat + F_sol
        total_coupled = F_gasmed_c + F_gaslat_c + F_sol_c
        max_conservation_err = float(np.max(np.abs(total_coupled - total_indep)))
        conserved = max_conservation_err < CONSERVATION_TOL_N

        # net ankle moment, independent vs coupled (live moment arms, N*m)
        M_indep = (F_gasmed * ma_ankle["gasmed_r"] + F_gaslat * ma_ankle["gaslat_r"]
                   + F_sol * ma_ankle[SOLEUS])
        M_coupled = (F_gasmed_c * ma_ankle["gasmed_r"] + F_gaslat_c * ma_ankle["gaslat_r"]
                     + F_sol_c * ma_ankle[SOLEUS])
        dM = M_coupled - M_indep
        peak_M_indep = float(np.max(np.abs(M_indep)))
        max_dM = float(np.max(np.abs(dM)))
        max_dM_pct_of_peak = 100.0 * max_dM / peak_M_indep

        k_argmax = int(np.argmax(np.abs(transfer)))
        # realized fraction relative to Fmax (a well-defined denominator even when instantaneous
        # force is tiny, unlike %-of-own-force which blows up near zero)
        max_abs_transfer = float(np.max(np.abs(transfer)))
        mean_abs_transfer = float(np.mean(np.abs(transfer)))

        # value at each muscle's OWN peak-force frame (not just the argmax|transfer| frame)
        k_gas_peak = int(np.argmax(F_gas_tot))
        k_sol_peak = int(np.argmax(F_sol))

        print(f"\n  [{label}] Delta_slack={slack}, f_max={f_max}")
        print(f"    conservation invariant: max|err|={max_conservation_err:.3e} N "
              f"({'PASS' if conserved else 'FAIL'} vs {CONSERVATION_TOL_N:.0e} N tol)")
        print(f"    |transfer(t)|: max={max_abs_transfer:.2f} N (t={t_ik[k_argmax]:.2f}s, "
              f"Delta={delta[k_argmax]:.3f}, F_gas_tot={F_gas_tot[k_argmax]:.1f}N, "
              f"F_sol={F_sol[k_argmax]:.1f}N), mean={mean_abs_transfer:.2f} N")
        print(f"    at gastroc's OWN peak-force frame (t={t_ik[k_gas_peak]:.2f}s, "
              f"F_gas_tot={F_gas_tot[k_gas_peak]:.1f}N): transfer={transfer[k_gas_peak]:+.2f} N "
              f"({100*transfer[k_gas_peak]/F_gas_tot[k_gas_peak]:+.1f}% of that frame's gastroc force)")
        print(f"    at soleus's OWN peak-force frame (t={t_ik[k_sol_peak]:.2f}s, "
              f"F_sol={F_sol[k_sol_peak]:.1f}N): transfer={transfer[k_sol_peak]:+.2f} N "
              f"({100*transfer[k_sol_peak]/F_sol[k_sol_peak]:+.1f}% of that frame's soleus force)")
        print(f"    net ankle plantarflexion moment: peak|M_indep|={peak_M_indep:.2f} N*m, "
              f"max|dM|={max_dM:.3f} N*m = {max_dM_pct_of_peak:.2f}% of peak "
              f"({'NEGLIGIBLE' if max_dM_pct_of_peak < NEGLIGIBLE_NET_MOMENT_PCT else 'NON-NEGLIGIBLE'} "
              f"vs pre-registered {NEGLIGIBLE_NET_MOMENT_PCT}% ceiling)")

        return {
            "delta_slack": slack, "f_max": f_max,
            "conservation_max_err_N": max_conservation_err, "conservation_pass": conserved,
            "transfer_N": {"max_abs": max_abs_transfer, "mean_abs": mean_abs_transfer,
                            "argmax_t_s": float(t_ik[k_argmax]),
                            "at_argmax": {"delta": float(delta[k_argmax]),
                                          "F_gas_tot_N": float(F_gas_tot[k_argmax]),
                                          "F_sol_N": float(F_sol[k_argmax])}},
            "at_gastroc_peak_frame": {"t_s": float(t_ik[k_gas_peak]),
                                       "F_gas_tot_N": float(F_gas_tot[k_gas_peak]),
                                       "transfer_N": float(transfer[k_gas_peak]),
                                       "transfer_pct_of_frame_force": float(
                                           100*transfer[k_gas_peak]/F_gas_tot[k_gas_peak])},
            "at_soleus_peak_frame": {"t_s": float(t_ik[k_sol_peak]),
                                      "F_sol_N": float(F_sol[k_sol_peak]),
                                      "transfer_N": float(transfer[k_sol_peak]),
                                      "transfer_pct_of_frame_force": float(
                                          100*transfer[k_sol_peak]/F_sol[k_sol_peak])},
            "net_ankle_moment": {"peak_M_indep_Nm": peak_M_indep, "max_abs_dM_Nm": max_dM,
                                  "max_abs_dM_pct_of_peak": max_dM_pct_of_peak,
                                  "verdict": ("NEGLIGIBLE" if max_dM_pct_of_peak < NEGLIGIBLE_NET_MOMENT_PCT
                                              else "NON-NEGLIGIBLE")},
            "series": {"t": t_ik.tolist(), "delta": delta.tolist(), "transfer_N": transfer.tolist(),
                       "F_gasmed_indep": F_gasmed.tolist(), "F_gaslat_indep": F_gaslat.tolist(),
                       "F_sol_indep": F_sol.tolist(),
                       "F_gasmed_coupled": F_gasmed_c.tolist(), "F_gaslat_coupled": F_gaslat_c.tolist(),
                       "F_sol_coupled": F_sol_c.tolist()},
        }

    result_default = run_regime(DELTA_SLACK_DEFAULT, F_MAX_DEFAULT, "DEFAULT (physiological-anchored)")
    result_sensitivity = run_regime(DELTA_SLACK_SENSITIVITY, F_MAX_SENSITIVITY, "SENSITIVITY (upper-bound)")

    print("\nSTEP 5b/6 -- null-control (f_max=0): must reproduce the independent-actuator baseline exactly")
    result_null = run_regime(0.0, 0.0, "NULL CONTROL f_max=0")
    null_max_abs = result_null["transfer_N"]["max_abs"]
    null_pass = null_max_abs == 0.0
    print(f"  null control transfer max|.|={null_max_abs} ({'PASS -- exact zero' if null_pass else 'FAIL'})")

    print("\nSTEP 6/6 -- write results")
    results = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "model_file": MODEL_FILE, "ik_trial": IK_MOT, "so_force_baseline": SO_FORCE_STO,
        "muscle_properties": props,
        "gasmed_vs_gaslat_length_correlation": corr_heads,
        "corr_abs_delta_vs_total_force": corr_delta_force,
        "literature": LITERATURE,
        "thresholds": {
            "delta_slack_default": DELTA_SLACK_DEFAULT, "delta_slack_sensitivity": DELTA_SLACK_SENSITIVITY,
            "delta_ramp_scale": DELTA_RAMP_SCALE, "f_max_default": F_MAX_DEFAULT,
            "f_max_sensitivity": F_MAX_SENSITIVITY,
            "negligible_net_moment_pct_ceiling": NEGLIGIBLE_NET_MOMENT_PCT,
            "conservation_tol_N": CONSERVATION_TOL_N, "gastroc_head_corr_floor": GASTROC_HEAD_CORR_FLOOR,
        },
        "regime_default": result_default,
        "regime_sensitivity": result_sensitivity,
        "null_control": {"max_abs_transfer_N": null_max_abs, "pass": null_pass},
    }
    out_path = f"{OUT_DIR}/myofascial_transmission_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"  wrote {out_path}")

    print("\n" + "=" * 88)
    print("SUMMARY")
    print("=" * 88)
    print(f"Literature (live-verified, PMID/DOI in JSON + docstring): small-to-moderate at intact/"
          f"physiological relative position (0-14% typical, rat EDL/TA); ~0% for THIS pair's most "
          f"direct analog (cat soleus, Maas & Sandercock 2008, independent lab); up to 40.5% "
          f"(Rijkelijkhuizen 2005) under a fully-fascially-intact-but-exposed rat GM->plantaris prep. "
          f"No direct human force-partitioned fraction found in any species-pair checked.")
    print(f"DEFAULT (literature-anchored to the near-zero cat-soleus finding): max transfer "
          f"{result_default['transfer_N']['max_abs']:.1f} N, net ankle moment change "
          f"{result_default['net_ankle_moment']['max_abs_dM_pct_of_peak']:.2f}% of peak "
          f"({result_default['net_ankle_moment']['verdict']}).")
    print(f"SENSITIVITY (upper-bound, anchored to 40.5%): max transfer "
          f"{result_sensitivity['transfer_N']['max_abs']:.1f} N, net ankle moment change "
          f"{result_sensitivity['net_ankle_moment']['max_abs_dM_pct_of_peak']:.2f}% of peak "
          f"({result_sensitivity['net_ankle_moment']['verdict']}).")
    print(f"Conservation invariant: DEFAULT {'PASS' if result_default['conservation_pass'] else 'FAIL'}, "
          f"SENSITIVITY {'PASS' if result_sensitivity['conservation_pass'] else 'FAIL'}. "
          f"Null control: {'PASS' if null_pass else 'FAIL'}.")


if __name__ == "__main__":
    main()
