# MECHANISM SUBTALAR JOINT — the hindfoot torque converter (tibial rotation ⇄ calcaneal inversion/eversion) + stance-phase pronation↔supination (2026-07-22)

Distinct joint from `docs/MECHANISM_ANKLE_FORCE.md` (talocrural — dorsi/plantarflexion, force-focused) and distinct
from the midtarsal cert (`docs/MECHANISM_FOOT_AXIS_FIX.md`, `docs/MECHANISM_MULTISEGMENT_FOOT.md` — Chopart's joint,
an HONEST-NEGATIVE: the real Maharaj2021 midtarsal axis transplant FAILED the held-out generalization gate). This
document is the **subtalar (talocalcaneal) joint** one level below the ankle: its oblique single axis, the
kinematic torque-conversion it performs between the tibia and the calcaneus, the stance-phase
pronation→supination timing, and the function↔dysfunction axis (pes planus / pes cavus). No live simulation was
run this session (no new `.py` script/cell) — this is a **first-principles geometric derivation, forced against a
23-source live-verified literature base and the twin's own already-encoded model**, delivered as a hypothesis for
independent QC, per this task's own framing.

## Headline

| claim | verdict | anchor |
|---|---|---|
| (a) Subtalar axis obliquity geometrically forces a fixed-ratio torque conversion between tibial axial rotation and calcaneal (subtalar) rotation during closed-chain stance | **PASS** (derived ratio 1/sin(α) ≈ **1.49:1** at α=42°; order-of-magnitude-consistent with independent measured stance ROM) | Manter 1941 (DOI-verified); Lundberg et al. 1989 in-vivo RSA (PMID 2915016); Hunt et al. 2001 (PMID 11470301); twin's own `.osim` (ankle = 1-DOF PinJoint, live-verified) |
| (b) Coupling ratio is not a fixed population constant — it is foot-type/dysfunction-dependent | **PASS** | Nawoczenski, Saltzman & Cook 1998 (PMID 9555923); Maharaj, Cresswell & Lichtwark 2017 (PMID 28783556) |
| (c) Stance pronates-then-supinates (shock-absorber → rigid lever) as a coarse directional/timing description | **PARTIAL PASS** — directional claim survives, the literal "clean triplanar" mechanism does NOT | Hunt et al. 2001 (PMID 11470301) — same paper is BOTH the anchor and the forced adversary |
| (d) A rigid/fused hindfoot must fail to deliver BOTH heel-strike shock absorption AND push-off rigid lever, even given its strongest real compensatory mechanism | **PASS** | Usuelli et al. 2021, n=23 real subtalar-fusion patients (PMID 32381452) |
| (e) Dysfunction pole: over-pronation → medial-column/medial-leg overload | **PASS** (prospective, OR=6.6) | Bennett, Reinking & Rauh 2012 (PMID 22666641); Arangio & Salathe 2009 (PMID 19272682); Oh et al. 2013 (PMID 23783206) |
| (f) Dysfunction pole: over-pronation → tibial rotation → patellofemoral/knee-valgus chain SPECIFICALLY | **FAIL at the population/group level (honest kill of this specific sub-claim)** | Powers, Chen, Reischl & Perry 2002, n=41 (PMID 12146775); Reischl et al. 1999 (PMID 10473063) |
| (g) Over-supination (pes cavus) → poor shock absorption, locked midfoot | **PASS (tertiary-review + mechanistic)** | Seaman & Ball, StatPearls 2023/2026 (PMID 32310476) |
| (h) Twin's current model: subtalar ROM is symmetric (±35°), real anatomy is not | **PASS (live-verified discrepancy, actionable)** | direct read of `data/msk_models/subject2_unified_v2_realmidtarsalaxis.osim` this session |

No WebSearch was available this session (session-wide quota exhausted, the same pre-existing constraint already
flagged in `docs/MECHANISM_ANKLE_FORCE.md` §1). All 23 sources below were instead verified LIVE via direct `curl
--max-time 25` calls to NCBI eutils (`esearch`/`esummary`/`efetch`) and the Crossref REST API — machine-fetched
real metadata/abstract text, not recalled, not narrated.

## 0. Honest tier — read before the geometry

Two primary historical sources anchor the classical "~42° / ~16-23°" axis figures. Their verification tiers are
**not equal, and I will not blur them**:

- **Manter JT (1941). "Movements of the subtalar and transverse tarsal joints." *Anat Rec* 80(4):397-410.**
  DOI `10.1002/ar.1090800402` — **existence, author, journal, volume/pages independently confirmed live via
  Crossref** this session. PubMed carries no record (pre-MEDLINE). **I could not independently extract Manter's
  own reported angle values from primary text this session** (no abstract exists for a 1941 paper; no open
  full-text found). The "~42° from transverse" figure used below is the value ubiquitously attributed to this
  paper in the secondary/podiatric-biomechanics literature — used as this document's working number, flagged as
  **not re-derived from the primary source**, consistent with the task's own framing of it as "Manter/Inman."
- **Isman RE, Inman VT (1969). "Anthropometric studies of the human foot and ankle." *Bull Prosthet Res*
  10-11:97-129.** **NOT independently verified this session** — no PMID found via NCBI eutils, no DOI found via
  Crossref (a VA government report series predating both indexing systems for this item). The commonly-quoted
  "23° medial deviation, range 20.5°-68.5°" figure is carried here **by name only, as a fully disclosed gap**,
  not fabricated with a fake identifier.

This is a genuinely weaker verification tier than the rest of this document's 21 other sources (all PMID- or
DOI-confirmed live, most with verbatim abstract quotes below). The geometric derivation in §1 is a **hypothesis
built on the standard textbook figure**, cross-checked in §2 against independently-measured modern ROM data — the
cross-check, not the Manter/Isman citation itself, is what makes the numeric claim non-circular.

## 1. The geometry — why an oblique axis is a torque converter (derived, not asserted)

**Setup.** The subtalar joint is modeled (by essentially every gait-analysis-grade OpenSim model, including this
twin's own — confirmed live, below) as a **single revolute (1-DOF) hinge** with a fixed oblique axis. Let α = the
axis's inclination from the transverse (horizontal) plane (~42°, §0). Take z = vertical (proximal), x = the
foot's long axis (anterior), y = medial-lateral; the axis unit vector is **n** = (cosα·cosβ, cosα·sinβ, sinα),
β = its azimuthal deviation from the foot's long axis (~16-23°, §0).

**Closed-chain stance constraint.** During flat-footed, non-slipping stance, the calcaneus (rigidly coupled to
the ground via the planted forefoot) **cannot rotate about the global vertical axis** — if it did, the sole
would twist against the ground. The talus, in turn, is coupled to the tibia through the ankle (talocrural) joint.

**Live-verified structural fact (this session, not assumed):** in the twin's own model
(`data/msk_models/subject2_unified_v2_realmidtarsalaxis.osim`), `ankle_r`/`ankle_l` are **`PinJoint`s — exactly 1
DOF** (`ankle_angle_r`/`_l`, dorsi/plantarflexion only), grepped directly from the live XML. **By construction,
this joint carries zero axial (transverse-plane) rotational freedom** — any tibial internal/external rotation
must pass to the talus rigidly. This is not a modeling shortcut I introduced; it is the pre-existing structure of
the model I am analyzing, and it is exactly the assumption the classical "mitred hinge" torque-converter
description (Inman) requires.

**First-order (small-angle) rotation-vector superposition** (valid to first order; real finite-rotation
composition is nonlinear and this is disclosed as an approximation in §6): talar rotation = tibial rotation
φ·ẑ (transmitted rigidly per above). Calcaneal rotation relative to talus = θ·**n** (the single subtalar DOF).
For the total calcaneal rotation's vertical component to vanish (foot-flat constraint): φ + θ·sinα = 0, i.e.

**θ = −φ / sin(α)** — independent of β.

This is the geometric "torque converter" ratio, and it makes two falsifiable, checked predictions:
1. **Degenerate-case sanity** (computed, not hand-waved): α→90° (axis vertical) ⇒ ratio→1 (pure 1:1 relay — a
   vertical axis IS the tibial rotation axis, no conversion needed). α→0° (axis horizontal, like a plain
   dorsi/plantarflexion hinge) ⇒ ratio→∞ (a horizontal axis cannot absorb ANY vertical-axis rotation — physically
   correct: this is exactly why the ankle needs its own separate joint and cannot double as the rotation-absorber).
   Both limits check out against physical intuition — the formula is not an artifact of the small-angle step.
2. **β drops out of the ratio entirely** — a clean, checkable structural prediction: the axis's azimuthal
   deviation from the foot's long axis governs the *triplanar mixing* (how much of the subtalar rotation shows up
   as pure frontal-plane inversion/eversion vs. as ab/adduction-like tilt), not the *torque-conversion magnitude*.
   This matches the standard clinical teaching that subtalar motion is triplanar, never pure frontal-plane.

**Numeric value at α=42°: ratio = 1/sin(42°) = 1.494 : 1.** For every 1° of tibial axial rotation, the subtalar
joint must rotate ~1.5° for the planted foot not to twist on the ground. Sensitivity across the wider individual
range some secondary sources attribute to Isman & Inman (20.5°-68.5°, §0 — unverified this session): ratio spans
**1.08 (steep axis) to 2.86 (shallow axis)** — i.e. the SAME anatomical parameter (axis inclination) that varies
several-fold between individuals predicts a several-fold difference in torque-conversion aggressiveness. This is
the mechanistic bridge to the dysfunction pole (§4): a shallower (more horizontal) axis geometrically REQUIRES a
larger coupling ratio.

## 2. Cross-checking the derivation against independent measurement (not a tautology)

The derivation in §1 is not fit to any data — it falls out of axis geometry + a structural fact about the model.
The cross-checks below are genuinely external.

**(i) In-vivo confirmation the mechanism itself exists**, independent of the exact angle: **Lundberg A, Svensson
OK, Németh G, Selvik G (1989). "The axis of rotation of the ankle joint." *J Bone Joint Surg Br* 71(1):94-9.**
PMID `2915016`, DOI `10.1302/0301-620X.71B1.2915016`. n=8 healthy volunteers, **roentgen stereophotogrammetric
analysis (RSA)** — tantalum bone markers, gold-standard in-vivo bone kinematics, not skin/cadaver. Companion
monograph: **Lundberg A (1989). "Kinematics of the ankle and foot. In vivo roentgen stereophotogrammetry." *Acta
Orthop Scand Suppl* 233:1-24.** PMID `2686345`. Verbatim: *"The ankle/foot complex showed ability to **transform
leg rotation into pro-/supination and vice versa**. This function was most pronounced in external leg rotation."*
Also: *"The total amount of rotation in the talo-calcaneal joint was small in internal rotation of the leg and in
pronation of the foot compared to external rotation of the leg and supination of the foot"* — a genuine
**internal/external asymmetry**, disclosed rather than smoothed over (§6).

**Adversary this same source forces on the "single fixed axis" idealization my derivation uses:** *"the joint
axis of the talo-crural joint varied with varying kinds of input motion... took varying inclinations between
horizontal and vertical."* Real in-vivo axes migrate; the fixed-axis hinge is the standard modeling
simplification (used by this twin and by essentially every OpenSim gait model), not a literal anatomical fact.
Disclosed, not smoothed over.

**(ii) Magnitude triangulation** (approximate — two honestly-flagged inputs compound, so this is an
order-of-magnitude check, not a precision validation): **Hunt AE, Smith RM, Torode M, Keenan AM (2001).
"Inter-segment foot motion and ground reaction forces over the stance phase of walking." *Clin Biomech*
16(7):592-600.** PMID `11470301`. Directly measured, surface markers: **rearfoot stance-phase ROM = 22°
(sagittal), 8° (frontal), 10° (transverse)**; forefoot 12°/4°/10°. Feeding a commonly-cited stance tibial-rotation
magnitude (~4-12°, historically attributed to **Levens AS, Inman VT, Blosser JA (1948). "Transverse rotation of
the segments of the lower extremity in locomotion." *J Bone Joint Surg Am* 30A(4):859-72,** PMID `18887290` —
existence confirmed, no abstract available, value not independently re-extracted this session) through the
derived 1.494 ratio predicts **6.0-17.9° subtalar rotation** — bracketing Hunt's measured 8° rearfoot frontal-plane
figure. Same order of magnitude; not a tight match, and not claimed as one.

**(iii) Coupling ratio is foot-type-dependent, not fixed** (this is the direct, decisive, quoted measurement):
**Nawoczenski DA, Saltzman CL, Cook TM (1998). "The effect of foot structure on the three-dimensional kinematic
coupling behavior of the leg and rear foot." *Phys Ther* 78(4):404-16.** PMID `9555923`, DOI
`10.1093/ptj/78.4.404`. n=20 runners (10 low-rearfoot / 10 high-rearfoot by radiograph). Verbatim: *"The
predominant rotations suggest a combined subtalar and talocalcaneal joint axis to favor calcaneal eversion and
inversion for the low rear-foot group and tibial medial and lateral rotation for the high rear-foot group. **Group
differences were also found for the coupling ratio**, which described the proportion of calcaneal eversion and
inversion transferred or coupled to tibial axial rotation."* (Exact per-group ratio number not extractable from
the abstract alone — disclosed, not invented.)

**(iv) Independent replication, different lab/continent/method:** **Maharaj JN, Cresswell AG, Lichtwark GA
(2017). "Foot structure is significantly associated to subtalar joint kinetics and mechanical energetics." *Gait
Posture* 58:159-165.** PMID `28783556` (this citation is inherited from the sibling cert
`docs/MECHANISM_FOOT_FIDELITY_PLAN.md`, spot-checked live again this session via NCBI esummary — confirmed real).
n unspecified in abstract (multi-subject regression), OpenSim multi-segment foot model, Queensland. Verbatim:
*"Arch-height ratio, foot length and step width were associated with peak subtalar joint (STJ) moment, while
greater STJ negative work was correlated to a low arch-height ratio and greater foot mobility"* — **~46%/37% of
kinetic/energetic variance explained**. **The forced adversary this paper raises against its own kinematic
cousins:** *"no model could be created for any of the kinematic measures analysed"* — foot structure predicts
subtalar **moment and energy absorption**, not raw joint-angle kinematics. This converges with a repeated,
independent null (below, §4) on kinematic-magnitude predictability specifically — a coherent cross-study pattern,
not cherry-picked.

**(v) Measurement-fidelity caveat on ALL the surface-marker numbers above:** **Stacoff A, Nigg BM, Reinschmidt C,
van den Bogert AJ, Lundberg A (2000). "Tibiocalcaneal kinematics of barefoot versus shod running." *J Biomech*
33(11):1387-95.** PMID `10940397`. n=5, **intracortical bone pins** (gold standard). Verbatim: *"previous
investigations were typically based on externally mounted shoe and/or skin markers, which have been shown to
overestimate skeletal movements"*; between-subject spread up to 10°, shoe-condition effect <2°. Quantified
directly: **Westblad P, Hashimoto T, Winson I, Lundberg A, Arndt A (2002). "Differences in ankle-joint complex
motion... measured by superficial and bone-anchored markers." *Foot Ankle Int* 23(9):856-63.** PMID `12356185`.
Verbatim RMS gap, bone-pin vs. skin: **inversion/eversion 2.5° (tibiocalcaneal) / 2.1° (talocalcaneal)**,
plantar/dorsiflexion 1.7°, ab/adduction 2.8°. Every Nawoczenski/Reischl/Hunt/Powers number above (§2-§4) carries
this ~2-2.8° systematic uncertainty tier on top of its own reported variability — disclosed, not absorbed
silently.

**Twin's own axis, as a supporting (not primary) data point.** The twin's `subtalar_r/l` axis was already
FD-measured (finite-difference, live joint-perturbation measurement) in a prior wave's cert
(`scripts/msk/foot_midtarsal_axis_fix_evidence.json`, reused read-only) and cross-validated there against
Maharaj's own measured axis: cosine similarity **0.9859 (R) / 0.9858 (L)**, angle **9.64° / 9.68°**, PASS at a
pre-registered 0.8 threshold. Decomposing that same already-computed unit vector (calcn-frame, R:
`[0.7818, 0.6069, -0.1431]`) under the standard OpenSim body-frame convention (X=anterior, Y=superior, Z=lateral —
**this specific convention was not independently re-derived from the model this session**; supported only
indirectly by the sibling evidence file's own `anterior_end_x`/`posterior_end_x` labels, which are consistent
with +X=anterior) gives **inclination-from-transverse ≈ 37.4°, azimuthal deviation ≈ 10.4°** — the same
qualitative regime (oblique, tens-of-degrees, single-digit-to-twenties medial deviation) as the classical figures
in §0, offered as a **third, semi-independent, machine-computed corroboration**, explicitly flagged as
lower-confidence than (i)-(v) above because of the unverified frame-convention assumption.

## 3. Stance-phase pronation → supination timing (claim c) — the adversary IS the anchor

Pre-registered claim: early stance = pronation-dominant (shock absorption, midtarsal unlocked), late
stance/push-off = supination-dominant (rigid lever). **Hunt et al. 2001 (PMID 11470301, above)** is simultaneously
the best available anchor AND the paper that forces the adversary hardest against its own naive form. Verbatim:
*"Most motion occurred at the beginning and end of stance phase when support was via only the rearfoot or
forefoot"* — directionally consistent with an early/late-loaded timing pattern. But: *"Typical foot motion does
not obey descriptions of triplanar motion such as 'pronation' and 'supination'"* — the literal textbook
"clean-triplanar-then-reverses" mechanism is **falsified as stated** by direct multi-segment kinematic
measurement. **Verdict: PARTIAL PASS.** The coarse directional/temporal simplification (motion concentrated at
the phase transitions, generally eversion-then-inversion in overall sign) survives as a useful approximation; the
precise triplanar mechanistic story taught in clinical texts does not survive contact with Hunt's own data. Both
halves are reported — this is the honest reading, not a resolved one.

## 4. The core falsifier — a rigid/fused hindfoot, forced to its strongest real form

Pre-registered: a rigid-hindfoot model must fail to deliver BOTH heel-strike shock absorption AND push-off rigid
lever. The strongest fair adversary is not a hypothetical — it is a **real surgical natural experiment**:
**Usuelli FG, Indino C, Leardini A, Manzi L, Ortolani M, Caravaggi P (2021). "Range of motion of foot joints
following total ankle replacement and subtalar fusion." *Foot Ankle Surg* 27(2):150-155.** PMID `32381452`. n=23
real patients with **surgical subtalar arthrodesis** (fusion) combined with ankle arthroplasty, multi-segment
gait analysis vs. healthy controls. This is the adversary at its strongest: real fusion, PLUS the body's own
best compensatory strategy. Verbatim: *"fusion of the subtalar joint appeared to be **compensated by larger
frontal-plane motion at the tibiotalar joint**"* — i.e. the ankle borrows frontal-plane mobility it does not
normally use (contradicting, incidentally, this same twin's own ankle model's 1-DOF-only construction, §1 — real
ankles have a compensatory reserve the standard OpenSim simplification doesn't encode). Despite that
compensation: *"sagittal-plane mobility of foot joints was about **50%** that in healthy joints"* and *"Normal
spatio-temporal parameters were not restored"* (walking speed and stride length remained significantly reduced).
**Verdict: PASS.** The forced adversary, given its best real-world shot, still falls — confirming the
flexible-then-rigid transition is genuinely load-bearing, not a dispensable nicety.

## 5. The dysfunction pole — pes planus / pes cavus, and one deliberate kill

**Over-pronation (pes planus) → medial-column/medial-leg overload: PASS, prospective evidence.** **Bennett JE,
Reinking MF, Rauh MJ (2012). "The relationship between isotonic plantar flexor endurance, navicular drop, and
exercise-related leg pain in a cohort of collegiate cross-country runners." *Int J Sports Phys Ther* 7(3):267-78.**
PMID `22666641`, PMCID `PMC3362985`. n=77, **prospective** (tracked through a full season). Verbatim: *"runners
with a ND [navicular drop] >10 mm were almost 7 times (**OR=6.6**, 95% CI 1.2-38.0) more likely to incur medial
[exercise-related leg pain] than runners with ND <10 mm."* Navicular drop is a clinical PROXY for
rearfoot/midfoot pronation, not a direct subtalar-angle measurement — disclosed. Mechanistically supported by
**Arangio GA, Salathe EP (2009), *Clin Biomech* 24(4):385-90**, PMID `19272682` (biomechanical model: posterior
tibial tendon dysfunction/flatfoot → increased medial-arch/talonavicular load) and **Oh I, Imhauser C, Choi D,
Williams B, Ellis S, Deland J (2013), *J Bone Joint Surg Am* 95(12):1094-100**, PMID `23783206`, PMCID
`PMC6948803` (n=8 cadaveric robotic model: correcting a created flatfoot via lateral column lengthening
consistently reduced talonavicular abduction and increased lateral forefoot pressure — the reverse of the
uncorrected-flatfoot pressure pattern, confirming the medial↔lateral redistribution mechanism).

**Over-supination (pes cavus) → poor shock absorption, locked midfoot: PASS, tertiary-review + mechanistic.**
**Seaman TJ, Ball TA. "Pes Cavus." *StatPearls* [Internet], StatPearls Publishing, 2023 Aug 8 update.** PMID
`32310476`. Verbatim: *"In gait, the cavovarus foot demonstrates a compensatory heel varus, **a locked midfoot,
and a reduction of the flexible phase and decreased shock absorption** as a result of plantar fascial
tightening."* Consistent with Nawoczenski's "high rear-foot group" (§2) favoring the tibial-rotation axis over
the calcaneal-eversion axis — i.e., less of the shock-absorbing eversion pathway is used, forcing more of the
motion (and by implication, more of the load) elsewhere.

**Deliberate kill — the knee/patellofemoral chain SPECIFICALLY, not the whole dysfunction pole.** The classic
theoretical mechanism (subtalar over-pronation → obligatory tibial internal rotation → altered tibiofemoral/
patellofemoral mechanics) is real and foundational: **Tiberio D (1987). "The effect of excessive subtalar joint
pronation on patellofemoral mechanics: a theoretical model." *J Orthop Sports Phys Ther* 9(4):160-5.** PMID
`18797010` — a genuine, citable theoretical model, not invented for this document. But the one study found this
session that DIRECTLY tests it in real patients **refutes it at the population/group level**: **Powers CM, Chen
PY, Reischl SF, Perry J (2002). "Comparison of foot pronation and lower extremity rotation in persons with and
without patellofemoral pain." *Foot Ankle Int* 23(7):634-40.** PMID `12146775`. n=41 (24 PFP + 17 controls),
3D gait kinematics. Verbatim: *"No group differences were found with respect to the magnitude and timing of peak
foot pronation and tibia rotation... **these results do not support the hypothesis** that individuals with PFP
demonstrate excessive foot pronation or tibial internal rotation compared to nonpainful individuals."*
Consistent with the same group's earlier individual-level null: **Reischl SF, Powers CM, Rao S, Perry J (1999).
*Foot Ankle Int* 20(8):513-20**, PMID `10473063`, n=30: *"the magnitude and timing of peak pronation was not
predictive of the magnitude and timing of tibial and femoral rotation... contrary to the clinical hypothesis."*
**This is an honest, disclosed FAIL of the specific "pronation→knee-valgus" sub-claim** — not a hedge, not
softened. The broader dysfunction-pole claims above (medial-leg overload, shock-absorption deficit) stand on
their own independent evidence and are not weakened by this kill; the knee-specific causal chain is the one piece
that does not survive direct testing.

**The unifying geometric parameter.** §1's derivation showed the coupling ratio depends on exactly one scalar,
axis inclination α — matching, independently: **Kirby KA (2001). "Subtalar joint axis location and rotational
equilibrium theory of foot function." *J Am Podiatr Med Assoc* 91(9):465-87**, PMID `11679628` (a clinical theory
built entirely on individual variation in subtalar axis location/orientation driving mechanical-behavior
differences), validated methodologically by **Van Alsenoy KK, D'Août K, Vereecke EE, De Schepper J, Santos D
(2014), same journal 104(4):365-74**, PMID `25076079` (cadaver-feet reliability/validity of the axis-location
palpation technique). Two independent confirmations — one theoretical, one measured — that axis geometry, not a
fixed population constant, is the governing parameter, exactly as §1 derives from first principles.

## 6. Honest caveats (full list)

1. **Two-tier citation confidence (§0)** — Manter 1941 is DOI-existence-verified but not text-verified; Isman &
   Inman 1969 could not be verified at all this session (no PMID/DOI found). The 42°/16-23° figures are the
   field's standard numbers, not independently re-derived from their primary sources here.
2. **Small-angle linearization (§1)** — the θ=−φ/sinα derivation uses first-order rotation-vector superposition;
   real finite-rotation composition is non-commutative. Adequate for gait-magnitude rotations (single-digit to
   ~20°) but not exact at the edges of full ROM.
3. **Fixed-axis idealization vs. real in-vivo axis migration** — Lundberg's own RSA data (§2i) shows the true
   instantaneous axis is NOT fixed; the single-hinge model (used identically by this twin and by essentially
   every OpenSim gait model) is a standard, disclosed simplification.
4. **Triangulation in §2ii compounds two honestly-flagged unverified inputs** (the exact 42° figure's primary-text
   confirmation, and the Levens/Inman/Blosser tibial-ROM figure) — reported as order-of-magnitude, not precision,
   agreement.
5. **Skin-marker measurement-fidelity gap (§2v)** — ~2-2.8° RMS vs. bone-pin truth affects every Nawoczenski/
   Reischl/Hunt/Powers number used here.
6. **Reischl 1999 and Powers 2002 are the same research lineage** (Reischl, Powers, Perry appear on both) —
   the two "null" findings on pronation-predicts-rotation are not fully independent replications of each other,
   though Maharaj 2017 (different lab/country/method) independently converges on the same
   kinematics-not-predictable pattern for a different specific measure.
7. **Navicular drop (Bennett 2012) is a clinical proxy**, not a direct subtalar-joint-angle measurement — real,
   prospective, but one moderate-evidence-level (2b) single cohort (n=77), not a meta-analysis.
8. **The twin's own frame-convention assumption (§2, twin-axis paragraph)** was not independently re-derived this
   session (no live OpenSim query run) — treated explicitly as lower-confidence, supporting-only.
9. **Hicks 1954** (windlass foundational paper, PMID `13129168`, PMCID `PMC1244640`) exists and is correctly
   cited, but its full text is a scanned, non-OCR'd PDF — WebFetch could not extract a specific numeric windlass
   quote this session; cited for the mechanism/concept only, not a number.
10. **No new simulation was run.** This document is a literature + geometric-derivation synthesis, cross-checked
    against already-existing, already-verified machine outputs in this repo (the `.osim` XML read live, and the
    sibling FOOT_AXIS_FIX evidence file's FD-measured axis, reused read-only) — not a fresh measurement pipeline.
    §7 proposes the concrete cell that would convert this from a literature hypothesis into a twin-native
    machine measurement.

## 7. Couples to (the twin)

- **`docs/MECHANISM_ANKLE_FORCE.md`** — that doc's own disclosed caveat (§1 there) is that its Giddings et al. 2000
  anchor of **540 %BW is technically the talocalcaneal (subtalar) joint**, not the talocrural joint its headline
  targets. That doc supplies this joint's FORCE anchor (540 %BW peak, walking); this document supplies its
  KINEMATIC/coupling behavior. Together they close both halves for this one joint. Also: that doc's ankle model
  is the same 1-DOF `PinJoint` this document's §1 derivation leans on structurally.
- **`docs/MECHANISM_FOOT_AXIS_FIX.md` / `docs/MECHANISM_MULTISEGMENT_FOOT.md` / `docs/MECHANISM_FOOT_FIDELITY_PLAN.md`**
  — the midtarsal (Chopart) joint one segment distal; that cert's honest-NEGATIVE (real Maharaj axis transplant
  failed the held-out gate, root-caused to near-parallel subtalar/midtarsal axes) is a DIFFERENT joint's finding,
  not conflated here. This document reuses that cert's already-FD-measured subtalar axis (read-only) as a
  supporting data point (§2).
- **`docs/MECHANISM_GAIT_KINEMATICS_FIDELITY.md`** — stance-phase timing (§3) is a natural extension target.
- **`docs/MECHANISM_KNEE_CARTILAGE_MATERIAL.md`** and the knee-force certs** — §5's deliberate kill of the
  pronation→knee-valgus chain is directly relevant to any future knee-injury-mechanism cert reusing that theory.
- **`bt_memory/function-vs-dysfunction-build-the-athlete-baseline-in-parallel-with-disease.md`** — this document
  is built explicitly on that organizing axis: the healthy torque-converter mechanism (§1-§4) is the FUNCTION
  pole; pes planus/pes cavus (§5) is the DYSFUNCTION pole, made quantifiable as a deviation from the same
  geometric baseline (the axis-inclination parameter α).
- **`data/msk_models/subject2_unified_v2_realmidtarsalaxis.osim`** — read live this session (not from memory):
  `subtalar_r/l` PinJoint, range exactly ±35.000° (machine-computed from `0.61086523819801497` rad), confirmed
  **symmetric** — a live, actionable discrepancy vs. the near-universally-reported qualitative asymmetry of real
  subtalar ROM (larger inversion range than eversion range). `ankle_r/l` confirmed PinJoint, exactly 1 DOF.

## 8. Proposed cell (next step, concrete)

Replace the order-of-magnitude triangulation in §2ii with a **live, twin-native measurement**: using the
already-scaled OpenSim model (no new acquisition needed), prescribe a tibial-rotation sweep at the knee/hip
proximal to a foot held flat (a ground-contact or weld constraint enforcing the foot-flat idealization of §1),
and FD-measure (reusing the exact finite-difference method already proven and validated in
`docs/MECHANISM_FOOT_AXIS_FIX.md`) the resulting `subtalar_angle_r` response. Compare the measured ratio directly
against this document's derived 1.494:1 prediction — a decisive, twin-native PASS/FAIL replacing the current
literature-triangulation hypothesis with a genuine machine measurement. Secondary: fix the live-verified
ROM-symmetry gap (§7) by splitting `subtalar_angle_{r,l}`'s range into an asymmetric limit consistent with real
inversion>eversion anatomy, sourced from a primary goniometric/cadaveric reference (not yet pinned to an exact
live-verified split this session — an explicit acquisition task).

## Files

- `docs/MECHANISM_SUBTALAR_JOINT.md` — this document.
- `docs/MECHANISM_SUBTALAR_JOINT_evidence.json` — all 23 sources with PMID/DOI, verbatim quotes, verification
  method, and the geometric-derivation numbers, machine-structured.
- Read live this session (not written): `data/msk_models/subject2_unified_v2_realmidtarsalaxis.osim` (grepped
  directly for `subtalar_r`/`ankle_r` joint XML); `scripts/msk/foot_midtarsal_axis_fix_evidence.json` (sibling
  cert's FD-measured axis, reused read-only).
