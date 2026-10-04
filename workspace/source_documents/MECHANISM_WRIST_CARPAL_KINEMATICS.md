# MECHANISM WRIST/CARPAL KINEMATICS — two-row linkage model + a structural falsifier the twin already fails (2026-07-22)

Fills the stated MSK gap: `docs/MECHANISM_ELBOW_WRIST_FORCE.md` measures wrist *force* but treats
the wrist as its two `.osim` coordinates (`wrist_flex_r`, `wrist_dev_r`) with zero carpal-row
content. This adds the missing layer: a literature-parametrized **two-row carpal linkage**
model (radiocarpal vs midcarpal ROM contribution; radius-vs-ulna load split via the TFCC) and a
**machine-verified structural audit** of the twin's own `.osim` files asking the task's own
falsifier question directly — does the twin's *actual, already-built* wrist architecture even
admit these partitions, or is it a "single-hinge" adversary by construction? Every number below
is machine-computed this session (`scripts/msk/wrist_carpal_kinematics.py`, exit 0) or a
live-reverified PMID (NCBI eutils, `curl`, this session) — no simulation was run; this is a
structural-XML audit + a cited geometric model, the same evidence tier as this repo's other
literature-anchored MSK certs.

**Status: HYPOTHESIS awaiting independent QC. Confidence tier: cadaveric/published-plausibility**
(no in-vivo instrumented carpal-load telemetry exists anywhere in the literature — confirmed
absent, not merely unsearched, by the pre-existing graph node this doc couples to, §9).

## Headline

| Question | Finding |
|---|---|
| **Does the twin's current wrist model already carry carpal-row fidelity?** | **NO — measured, not assumed.** 0/5 audited `.osim` files contain any of 8 carpal bone names (scaphoid/lunate/triquetrum/capitate/hamate/trapezium/trapezoid/pisiform). The 4 "articulated" models carry exactly **2** generalized coordinates at the wrist (`wrist_flex_r`, `wrist_dev_r`) via **2 serial `PinJoint`s** through one **massless (1e-5 kg) dummy body** (`wrist_int_r`) — a single joint-*complex* standing in for the entire 8-bone, 2-row carpus. |
| **Is the twin's real architecture the task's own "single-hinge" adversary, or a strawman?** | **It IS the real architecture — not hypothesized.** 4 files (`subject2_hand_complete`, `LaiArnoldModified2017_full_hand_subject2_scaled`, `subject2_unified_v2`, `LaiArnoldModified2017_anatomical_hand_index_subject2_scaled`) share this identical 2-DOF gimbal; a 5th, older file (`subject2_unified.osim`) has the even-more-extreme **0-DOF welded** version (`radius_hand_r`, a single fused rigid body — the literal "rigid carpus" pole). |
| **Does this architecture structurally admit an emergent radiocarpal/midcarpal ROM split?** | **NO — falsified by DOF-counting, not opinion.** Row-attribution needs ≥4 generalized coordinates (`{rc_flex, rc_dev, mc_flex, mc_dev}`, a floor); the twin's wrist carries exactly 2. There is no second joint or intermediate mass to assign a "midcarpal" fraction to — the question is category-ill-posed for this architecture, not merely mismatched. |
| **Does this architecture structurally admit an emergent radius/ulna load split?** | **NO — falsified by graph reachability, not opinion.** The wrist joints' sole proximal parent is `radius_r`; `ulna_r` connects to the rest of the arm only via a *proximal* forearm joint (`radioulnar_r`, elbow-level pronation-supination) that has no direct edge to `hand_r` or any hand descendant. A `JointReaction` cut at the wrist is structurally 100% radius / 0% ulna — the same class of finding as `docs/MECHANISM_PATELLOFEMORAL_FORCE.md`'s quad:PT ratio ≡ 1.0 (a mathematical certainty of the architecture, not a simulated result). |
| **What does the real two-row linkage's ROM partition look like (verified)?** | **Direction-dependent, not a fixed split.** Sarrafian et al. 1977 (PMID 598105, n=55): flexion = **40% radiocarpal / 60% midcarpal**; extension = **66.5% radiocarpal / 33.5% midcarpal** — a **26.5 percentage-point swing** the twin's 2-DOF gimbal has no coordinate to represent even qualitatively. |
| **What does the real radius/ulna load split look like (verified)?** | **81.6% radius / 18.4% ulna at neutral ulnar variance** (Werner, Glisson, Murphy, Palmer 1986, PMID 3770570, in-vitro axial load-cell) — within 1.6 percentage points of the task's own pre-registered ~80/20 figure — but this split is **not a constant**: it swings from 4.3% to 41.9% ulnar share across a 5 mm ulnar-variance range, and the TFCC disc alone carries 66.3% of the neutral-variance ulnar load (removal drops ulnar share from 18.4% to 6.2%). |
| **Dart-thrower's motion — verified, and does it survive fusion of either row?** | **YES, robustly.** Crisco et al. 2005 (PMID 16322624, n=28, in vivo CT, 504 wrist positions): scaphoid/lunate rotation is **significantly less** along the DTM path than any other direction (p<0.01 both bones). Kane et al. 2018 (PMID 29146510, n=6 cadaveric, moment-controlled): the DTM axis orientation is **unchanged** by simulated radiocarpal OR simulated total-carpal fusion — both rows can independently support it (a redundancy property), even though ROM *magnitude* drops with either fusion. |
| **Dysfunction pole — DISI, TFCC tear, carpal tunnel?** | All 3 verified. DISI: Watson & Ballet 1984 (PMID 6725894, n=210/4000 films) — 57% of wrist-arthritis follows the scaphoid-lunate-radius (SLAC) collapse pattern; Rhee et al. 2009 (PMID 19833447, n=58) gives an independently-anchored numeric criterion, radiolunate angle >15°. TFCC tear: Palmer & Werner 1981 (PMID 7229292, n=61) — 53% of specimens show TFCC perforation, and **100%** of perforated specimens show lunate/ulna cartilage erosion. Carpal tunnel: Gelberman et al. 1981 (PMID 7204435) — neutral-position canal pressure is **12.8×** higher in CTS patients than controls (32 vs 2.5 mmHg), a *larger* patient:control ratio than at 90° flexion (3.03×) or extension (3.67×) — patients lose the normal pressure-relief neutral posture entirely, a genuinely computed, non-obvious result. |

## 1. Pre-registration (thresholds fixed before reading the twin's own files)

- **C (accept)**: the twin's wrist architecture structurally admits (a) a direction-dependent
  radiocarpal/midcarpal ROM partition and (b) an emergent (non-imposed) radius/ulna load split
  near 80/20. Threshold: ≥4 generalized coordinates at the wrist (for (a)) AND a direct joint/force
  edge between `ulna_r` (or a carpal proxy) and a hand-side body that does not pass through
  `radius_r` first (for (b)).
- **¬C (the single-hinge/rigid-carpus adversary)**: exactly what the task named — reject only if
  the adversary is FORCED to its real, already-built form (not a hypothetical strawman) and
  FAILS both thresholds by a machine-checkable graph/DOF fact, not narration.
- **Literature-side falsifier**: does a real, verified, direction-dependent ROM split exist
  (Sarrafian 1977) and does a real, verified, variance-sensitive load split exist (Werner 1986)
  close to the task's own ~80/20 figure (tolerance: ±5 percentage points)? Both are external,
  independently-measured anchors — never fit to the target.

## 2. The twin's current wrist architecture — measured, not assumed

`scripts/msk/wrist_carpal_kinematics.py` parses (pure `xml.etree.ElementTree`, no OpenSim API
dependency, no regex-on-XML fragility) 5 already-existing `.osim` files, builds the body/joint
graph from `socket_parent_frame`/`socket_child_frame` → `PhysicalOffsetFrame` → `socket_parent`
resolution (the same frame-indirection convention used throughout this model family), and asks
graph-theoretic questions of it.

**Result, identical across 4 independently-evolved hand-lineage files** (different hand
complexity — from 3 to 17 finger bodies — but the SAME wrist joint definition, a genuine
over-determination that this isn't a one-file artifact):

```
... -> ulna_r --(radioulnar_r, PinJoint)--> radius_r --(wrist_flex_r_joint, PinJoint)-->
wrist_int_r [mass = 1.0e-05 kg] --(wrist_dev_r_joint, PinJoint)--> hand_r -> {5 fingers}
```

- **0 carpal bones** in any file (checked against scaphoid/lunate/triquetrum/capitate/hamate/
  trapezium/trapezoid/pisiform — 0/5 hits, machine-counted).
- **Exactly 2 generalized coordinates** at the wrist (`wrist_flex_r` range [-70°,+70°],
  `wrist_dev_r` range [-25°,+35°] — read directly from each `Coordinate`'s own `<range>`, not
  estimated).
- **`wrist_int_r` is a massless numerical dummy** (1e-5 kg, 5 orders of magnitude below any real
  segment in the model — `radius_r`/`ulna_r` are each 0.631 kg), not a proximal-row proxy.
- **A 5th, older file** (`subject2_unified.osim`) has the even-more-extreme pole: a single
  `WeldJoint` (`radius_hand_r`) fuses `hand_r` directly to `radius_r` — **0 wrist DOF at all**,
  the literal rigid-carpus adversary the task named, also already built and in the repo's own
  history, not invented for this doc.

### 2.1 A self-caught false start (symmetric QC on my own first pass)

My first version of this audit tested "is `ulna_r` anywhere in `hand_r`'s full kinematic-tree
ancestor chain" and got **FAIL** (surprising — `ulna_r` IS an ancestor). Forced Orient rather than
accepting a premature negative: `ulna_r` is reached only because it sits **proximal to `radius_r`**
in the *elbow/forearm* chain (`humerus_r -> ulna_r -(elbow)-> ... -(radioulnar_r)-> radius_r`,
the standard convention where the ulna is the fixed forearm reference and the radius spins around
it for pronation-supination) — a completely different, serial, *upstream* relationship, not a
parallel distal (TFCC-like) contact surface. The corrected, precise test — is there a **direct**
joint/force edge between `ulna_r`/`ulna_l` and `hand_r`-or-any-hand-descendant that bypasses
`radius_r`? — is what actually answers the physiological question, and the answer is **no** (the
only 4 joints touching `ulna_r`/`ulna_l` anywhere in the model are `elbow_r/l` and
`radioulnar_r/l`, both proximal-forearm joints whose other endpoint is `humerus_r/l` or
`radius_r/l`, never a hand-side body). This is the same discipline this repo's sibling docs
already apply (e.g. the `TRIlong` sign-anomaly catch in `MECHANISM_ELBOW_WRIST_FORCE.md` §4) —
fixed at the source in the script, not narrated around.

## 3. The two-row linkage — ROM partition (literature model)

**Geometric frame**: the wrist is a **2-joint serial chain** (radiocarpal: radius↔proximal row;
midcarpal: proximal row↔distal row), so for any given motion plane, total hand rotation
θ_total = θ_RC + θ_MC — an exact linkage-superposition identity, not a heuristic. The
**distal row behaves as a single quasi-rigid functional unit** (Kobayashi & An, *Hand Clin*
1997, PMID 9048189: "the bones of the distal carpal row behave as a single, functional unit with
little intercarpal motion and are directly linked to movement of the third metacarpal" — this is
literally what the twin's `hand_r` body already IS, geometrically, just without the proximal-row
joint proximal to it). The **proximal row is not rigid** — it is the "intercalated segment," and
its own internal kinematics are non-trivial: Kobayashi et al. 1997 (*J Biomech*, PMID 9239563,
n=22 cadavers, biplanar radiography) found the scaphoid has the **greatest** rotation magnitude
within the proximal row, the lunate the **least**, during sagittal-plane wrist motion — i.e. the
"rigid intercalated segment" is itself a first-order approximation, not exact.

**Flexion-extension partition is direction-dependent, quantitatively** — Sarrafian, Melamed,
Goshgarian 1977 (*Clin Orthop Relat Res*, PMID 598105, n=55 normal wrists, radiographic):

| direction | radiocarpal % | midcarpal % |
|---|---:|---:|
| maximum flexion | 40.0 | 60.0 |
| maximum extension | 66.5 | 33.5 |

with a documented mechanism: "functionally the scaphoid belongs to the 1st carpal row in flexion
and to the second row in extension" — a genuine kinematic row-switching, independently
corroborated 20 years later by a completely different method (Kobayashi's biplanar cadaveric bone
tracking finding the scaphoid's rotation anomaly within the "proximal" row, above) — two
decorrelated methods/eras converging on the same non-obvious anatomical behavior.

**The axis of rotation is fixed near the capitate head, for BOTH flexion-extension and
radioulnar deviation** (Youm, McMurthy, Flatt, Gillespie 1978, *JBJS Am*, PMID 670263, n=6
cadavers + 6 volunteers + 100 wrist radiographs) — meaning even a *well-placed* single hinge (at
the capitate head, not wherever `wrist_int_r` happens to sit) could reproduce the **net** rotation,
but never the **internal RC/MC partition** — that partition is a hidden/latent state a
single-axis model has no coordinate to carry, independent of where the axis is placed.

**Radioulnar deviation partition**: Kobayashi et al. 1997 found distal-row RUD magnitude is
"generally of a greater magnitude" than proximal-row RUD — i.e. midcarpal > radiocarpal for this
plane too, **qualitatively** confirmed but **no exact percentage split could be extracted from an
accessible abstract this session** (honest gap, §9 — unlike flexion-extension, where Sarrafian's
number is exact and quoted).

**Dart-thrower's motion (DTM)** — the functionally-dominant, low-proximal-row-motion oblique
plane the task named:
- Moritomo et al. 2007 (IFSSH biomechanics committee report, *J Hand Surg Am*, PMID 17996783):
  DTM = radial-extension ↔ ulnar-flexion, oriented **30°–45° from the sagittal plane**; "the DTM
  utilizes the midcarpal joint to a great extent"; "less scaphoid and lunate motion than during
  pure flexion-extension or radioulnar deviation."
- Crisco et al. 2005 (*JBJS Am*, PMID 16322624, n=28 healthy subjects, 504 CT-derived wrist
  positions): scaphoid/lunate rotation **significantly less** along the DTM path than any other
  direction (p<0.01 for both bones); rotation varies near-linearly with motion direction
  (R²=0.90 scaphoid, 0.82 lunate).
- Kane et al. 2018 (*J Hand Surg Am*, PMID 29146510, n=6 cadavers, ±1.5 N·m moment-controlled
  testing across 24 directions, 3 states: intact / simulated-radiocarpal-fusion /
  simulated-total-carpal-fusion): the DTM axis **orientation** survives BOTH simulated fusions
  unchanged (no significant difference between states) even though ROM **magnitude** drops with
  either — i.e. the RC and MC joints are functionally **redundant** for the DTM axis specifically,
  a robustness property distinct from, and not contradicted by, the direction-dependent
  FE-partition finding above (Sarrafian) — those are two different, compatible facts (a fixed
  *net* axis orientation vs. a variable *row* contribution to reaching any given point on it).
- **Reachability sanity check** (coarse, disclosed as such): the twin's own `wrist_dev_r`/
  `wrist_flex_r` coordinate-range box subtends 19.7°–26.6° from the flexion axis at its own
  corners — undershooting the literature's 30–45° figure, because `wrist_dev_r`'s box
  (-25°/+35°) is narrower than `wrist_flex_r`'s (±70°). This is a box-corner geometry check on
  generic (likely non-subject-specific) coordinate defaults, **not** a re-derivation of Kane's
  actual moment-controlled ROM-envelope measurement — flagged, not conflated.

**Falsifier verdict, ROM side**: the 2-DOF gimbal can reach a DTM-*direction* pose trivially (any
linear combination of 2 orthogonal DOFs spans an oblique direction) — that was never in question.
What it categorically cannot do is attribute **how much of that pose came from which row**, or
represent that the row-partition itself is asymmetric between flexion and extension (Sarrafian)
— because it has no 2nd joint. **DOF deficit, machine-counted**: row-attribution needs a floor of
**4** generalized coordinates (`rc_flex, rc_dev, mc_flex, mc_dev`); every audited articulated
model has **2**. Deficit = 2, a clean, pre-registered, non-hand-wavy number.

## 4. The radius/ulna load split via the TFCC (literature model)

This session found a **pre-existing OPEN node already in this repo's own anchor graph**,
`CARPAL-LOAD-PARTITION-TFCC-SCAPHOID` (`data/MECHANISM_ANCHOR_GRAPH.json`), carrying exactly this
load-partition literature — built by a prior session, focused on load/pressure partition and
SNAC-wrist natural history, with **zero** overlap with the ROM/DTM kinematics content above. Per
this repo's own "additive, not duplicative" doctrine, this doc **couples to** that node rather
than re-deriving its numbers — but every PMID it cites was **independently re-verified live this
session** (NCBI eutils `efetch`, verbatim abstract match, not blind trust of a prior fold):

| citation (PMID, live-reverified) | n | method | key number |
|---|---:|---|---|
| Werner, Glisson, Murphy, Palmer 1986 (3770570) | in-vitro cadaveric | axial load-cell | ulna bears **18.4%** of axial load at neutral variance (radius 81.6%); **41.9%** at +2.5mm ulnar-plus; **4.3%** at -2.5mm ulnar-minus; TFCC-disc excision alone drops ulnar share to **6.2%** |
| Kazuki, Kusunoki, Shimazu 1991 (1861018) | cadaveric, 27 positions | pressure-sensitive film densitometer | radiolunate peak pressure **+27%** at 2.5mm ulnar-minus, **-22%** at 2.5mm ulnar-plus — direction-consistent with Werner (decorrelated method/lab) |
| Viegas et al. 1987 Part I (3693853) | n=5 cadavers, 36 positions | Fuji pressure film | scaphoid contact area = **1.47×** lunate's; mean peak pressure 3.17 MPa; contact area only 20.6% of available joint surface |
| Palmer & Werner 1981 (7229292) | n=61 specimens | anatomic dissection + biomechanical testing | TFCC perforated in **53%** of specimens; **100%** of perforated specimens show lunate/ulna cartilage erosion |
| Miyamura et al. 2023 (37471563) | n=51 SNAC nonunion + 50 controls | in-vivo 3D-CT, modern, non-cadaveric | distal row pronates +14°, deviates ulnarly +19°, shifts dorsally +17% (all p<0.001) vs. controls — the same load-repartition mechanism, an independent (in-vivo, not cadaveric) modality |

**Computed** (this session, from Werner 1986's verbatim numbers): the TFCC disc alone accounts
for (18.4−6.2)/18.4 = **66.3%** of the neutral-variance ulnar load-bearing capacity; the
ulnar-share slope is ≈**7.5 percentage points per mm** of ulnar-variance shift across the tested
±2.5mm range — i.e. the "~80/20" split is not a constant, it is a steep function of a
millimeter-scale anatomic parameter, and the task's own ~80/20 figure (measured here: 81.6/18.4)
sits within **1.6 percentage points** of Werner's direct measurement — a genuine, non-tautological
external match (Werner's cadaveric force-transducer study was not fit to the task's figure; the
task's figure is textbook-inherited from this same body of work).

**Falsifier verdict, load side**: the audited wrist architecture (§2) has **no direct joint/force
edge from `ulna_r` to any hand-side body** — the only path is proximal-serial, through `radius_r`.
A standard `JointReaction` cut at `wrist_flex_r_joint`/`wrist_dev_r_joint` therefore attributes
100% of the hand-side reaction to `radius_r`, structurally, always — the identical failure mode
`docs/MECHANISM_PATELLOFEMORAL_FORCE.md` §5 already found for the quad:patellar-tendon ratio (a
single continuous `PathActuator` forces that ratio ≡ 1.000 "by construction, not a simulation
result" — here, a single serial kinematic path forces the radius/ulna split to 100/0 by the same
kind of construction).

## 5. Dysfunction pole

**Scapholunate (SL) ligament rupture → DISI.** Linscheid, Dobyns, Beabout, Bryan 1972
(*JBJS Am*, PMID 4653642 — existence/authors/year/journal verified live; a 2002 "Classic Article"
reprint, PMID 11792792, DOI 10.2106/00004623-200201000-00020, independently confirms the citation
exists and is still referenced 30 years later; **neither indexed record carries an abstract**
pre-1975/reprint-without-abstract, so the primary numeric SL-angle thresholds — textbook-standard
~30-60° normal, >70-80° DISI — could not be independently re-extracted from primary text this
session, an honest gap, same treatment as `MECHANISM_ELBOW_WRIST_FORCE.md` gave An 1981's number).
**Independently-quotable, decorrelated numeric cross-check**: Rhee, Moran, Shin 2009
(*J Hand Surg Am*, PMID 19833447, n=58 patients) gives a directly-quoted, different-metric DISI
criterion: **radiolunate angle > 15°**. Watson & Ballet 1984 (PMID 6725894, n=210 cases /
4000 films reviewed): **57%** of degenerative wrist arthritis follows the scaphoid-lunate-radius
(SLAC) collapse pattern — the clinical end-state of untreated SL instability, a large,
quantitative, verified anchor.

**TFCC tear → altered ulnar load.** Palmer & Werner 1981 (PMID 7229292, §4 table): perforation
in 53% of 61 specimens, **all** of which show lunate/ulna cartilage erosion — the structural
correlate of the load-redistribution mechanism in §4 (a torn/perforated TFCC no longer carries
its 66.3%-of-ulnar-share role, concentrating load elsewhere — exactly what Miyamura 2023's in-vivo
bone-density mapping independently shows at the radiolunate joint in a different, chronic,
SNAC-pattern population).

**Carpal tunnel (median-nerve compression), the neuromuscular-dysfunction pole.** Gelberman,
Hergenroeder, Hargens, Lundborg, Akeson 1981 (*JBJS Am*, PMID 7204435, n=15 CTS patients + 12
controls, wick-catheter pressure): patients 32/94/110 mmHg (neutral/90°flex/90°ext) vs. controls
2.5/31/30 mmHg. **Computed this session**: patient:control ratio is **12.8×** at neutral,
**3.03×** at flexion, **3.67×** at extension — the neutral-position ratio is by far the largest,
a genuinely computed (not narrated) and non-obvious result: CTS patients lose the normal
pressure-relief neutral posture entirely, rather than merely exaggerating the flexion/extension
pressure rise everyone experiences.

## 6. Machine cross-checks (8/8 PASS)

| check | result |
|---|---|
| Zero carpal bone bodies in any of 5 audited `.osim` files | PASS (0/5) |
| Articulated-wrist models carry exactly 2 generalized coordinates at the wrist | PASS (4/4 models) |
| No direct `ulna_r`↔hand-side joint/force edge in any articulated model (corrected test, §2.1) | PASS (4/4) |
| No real (non-dummy, >1e-4 kg) intermediate segment between `hand_r` and `radius_r` | PASS (4/4; only the 1e-5 kg `wrist_int_r` dummy found) |
| Welded-control model (`subject2_unified.osim`) has exactly 0 wrist DOF | PASS |
| DOF deficit for row-attribution: need ≥4, measured max 2 | PASS (deficit confirmed) |
| Werner 1986 neutral ulnar-share (18.4%) within 5pp of the task's ~20% figure | PASS (diff 1.6pp) |
| `wrist_flex_r`/`wrist_dev_r` coordinate ranges extracted from live XML | PASS |

Full detail: `data/msk_smoketest/wrist_carpal_kinematics/wrist_carpal_kinematics_results.json`.

## 7. External anchors — all PMIDs live-reverified this session (NCBI eutils, `curl --max-time 25`)

ROM/kinematics (newly searched this session): Kobayashi et al. 1997 *J Biomech* (PMID 9239563,
DOI 10.1016/s0021-9290(97)00026-2); Kobayashi & An 1997 *Hand Clin* (PMID 9048189); Moojen et al.
2002 *Clin Biomech* (PMID 12206941, DOI 10.1016/s0268-0033(02)00038-4) and 2003 *J Hand Surg Am*
(PMID 12563642, DOI 10.1053/jhsu.2003.50009) — the latter's own honest conclusion, quoted:
"a single functional model of carpal kinematics could not be determined" (in-vivo inter-subject
variation exceeds in-vitro); Crisco et al. 2005 (PMID 16322624, DOI 10.2106/JBJS.D.03058);
Moritomo et al. 2007 (PMID 17996783, DOI 10.1016/j.jhsa.2007.08.014) and 2014 follow-up committee
report (PMID 24888529); Kane et al. 2018 (PMID 29146510, DOI 10.1016/j.jhsa.2017.10.017);
Sarrafian, Melamed, Goshgarian 1977 (PMID 598105); Youm, McMurthy, Flatt, Gillespie 1978
(PMID 670263).

Load/TFCC (independently re-verified from the pre-existing graph node `CARPAL-LOAD-PARTITION-
TFCC-SCAPHOID`, §4): Werner, Glisson, Murphy, Palmer 1986 (PMID 3770570); Kazuki, Kusunoki,
Shimazu 1991 (PMID 1861018, DOI not registered in PubMed record); Viegas et al. 1987 Part I
(PMID 3693853) and Part II (PMID 3693854, perilunate instability, not independently re-fetched
this session — cited via the existing node); Palmer & Werner 1981 TFCC anatomy (PMID 7229292);
Palmer & Werner 1984 DRUJ biomechanics review (PMID 6744728, Clin Orthop Relat Res — mechanism
confirmed live, no extractable % in abstract, §9); Miyamura et al. 2023 (PMID 37471563).

Dysfunction: Watson & Ballet 1984 (PMID 6725894, DOI 10.1016/s0363-5023(84)80223-3); Linscheid
et al. 1972 (PMID 4653642) + 2002 reprint (PMID 11792792, DOI 10.2106/00004623-200201000-00020);
Rhee, Moran, Shin 2009 (PMID 19833447, DOI 10.1016/j.jhsa.2009.06.017); Gelberman et al. 1981
(PMID 7204435).

## 8. Couples to

- `docs/MECHANISM_ELBOW_WRIST_FORCE.md` — this doc adds the carpal-row kinematic layer that
  doc's own wrist-force estimate has zero content on (it treats the wrist as its 2 bare
  coordinates); its walking-data-gate finding (`wrist_flex_r` exactly 0.000° for all 158 gait
  frames — no hand/wrist mocap markers) is unaffected by anything here.
- `CARPAL-LOAD-PARTITION-TFCC-SCAPHOID` (existing OPEN node, `data/MECHANISM_ANCHOR_GRAPH.json`)
  — this doc's §4 couples to, re-verifies, and extends it; adds the ROM/DTM/DISI/carpal-tunnel
  content that node does not have.
- `docs/MECHANISM_HAND_COMPLETE.md`, `docs/MECHANISM_HAND_FIVE_DIGIT.md`, `docs/MECHANISM_ARM_MUSCLES.md`
  — the hand/grip models this doc's structural audit directly parsed.
- The function↔dysfunction axis: DISI / TFCC tear / carpal tunnel (§5) as the wrist's own
  instance of the pattern this repo's other joint docs already use (e.g. knee ligament laxity,
  ankle instability).

## 9. Honest gaps

1. **No in-vivo instrumented carpal-load telemetry exists anywhere** (confirmed absent by the
   pre-existing `CARPAL-LOAD-PARTITION-TFCC-SCAPHOID` node's own search, not re-verified as a
   negative independently this session) — same "no instrumented-implant analog conceivable for
   the small multi-bone carpus" gap that node already flagged, worse than the knee/hip case.
2. **No exact radioulnar-deviation RC/MC percentage split found** — unlike flexion-extension
   (Sarrafian's exact 40/60 vs 66.5/33.5), the RUD partition is only qualitatively confirmed
   (midcarpal > radiocarpal, Kobayashi 1997) after 2 additional targeted PubMed searches this
   session; may exist in a source not surfaced by the search terms tried.
3. **Palmer & Werner's own 1984 CORR review (PMID 6744728) states no extractable percentage in
   its abstract** — the ~80/20 figure is anchored here via the DIFFERENT, more precise Werner
   1986 in-vitro force-transducer paper instead (§4), which does carry the exact number; the 1984
   review is retained only as a qualitative mechanism anchor, disclosed, not conflated.
4. **The SL-angle numeric thresholds (30-60° normal, >70-80° DISI) are textbook-inherited, not
   independently extracted from Linscheid 1972's primary text** (pre-1975, no indexed abstract in
   either of 2 independent PubMed records checked) — the Rhee 2009 radiolunate-angle criterion
   (>15°) is a genuine, independently-verified numeric anchor for the same diagnosis, but uses a
   different landmark pairing, not a direct confirmation of the SL-angle number itself.
5. **The DTM box-reachability check (§3) is a coarse sanity check on generic coordinate-range
   defaults**, not a re-derivation of Kane 2018's actual moment-controlled ROM-envelope-
   orientation measurement — both are reported, not conflated.
6. **This is a literature-parametrized model, not a new subject-specific simulation** — no new
   carpal bones were added to any `.osim` file this session (that is the proposed next step,
   §10); the structural audit and DOF/graph facts are the only newly-measured numbers against
   the twin's own files, everything else is a cited external figure.
7. **Single-subject-agnostic**: the literature parameters (Sarrafian, Werner, etc.) are
   population averages from their own cadaveric/in-vivo cohorts, not subject2-specific
   measurements — same scope limit as every other cited-literature MSK cert in this repo.

## 10. Files

- `scripts/msk/wrist_carpal_kinematics.py` — self-contained, re-runnable, exit 0. Pure
  `xml.etree.ElementTree` parse of 5 `.osim` files (no OpenSim API dependency); builds the
  body/joint ancestor graph; encodes the literature model; runs 8 machine checks.
- `data/msk_smoketest/wrist_carpal_kinematics/wrist_carpal_kinematics_results.json` — full
  structural audit (all 5 models) + literature model + all 8 checks + verdict.
- Read-only, unmodified: `data/msk_models/subject2_hand_complete.osim`,
  `data/msk_models/LaiArnoldModified2017_full_hand_subject2_scaled.osim`,
  `data/msk_models/subject2_unified_v2.osim`,
  `data/msk_models/LaiArnoldModified2017_anatomical_hand_index_subject2_scaled.osim`,
  `data/msk_models/subject2_unified.osim`, `data/MECHANISM_ANCHOR_GRAPH.json` (read for the
  pre-existing `CARPAL-LOAD-PARTITION-TFCC-SCAPHOID` node, not modified).
- Builds on, does not re-litigate: `docs/MECHANISM_ELBOW_WRIST_FORCE.md`,
  `docs/MECHANISM_HAND_COMPLETE.md`, `docs/MECHANISM_PATELLOFEMORAL_FORCE.md` (the structural-
  degeneracy finding pattern this doc's §4 verdict directly mirrors).

## 11. Proposed cell / next step

Graft a genuine two-row carpus onto `subject2_hand_complete.osim`: a real (non-dummy,
inertia-bearing) proximal-row body (scaphoid+lunate+triquetrum, first-order as one rigid segment
— the "intercalated segment" approximation) between `radius_r` and `hand_r`, with **2 new**
generalized coordinates (radiocarpal flex/dev) inserted proximal to the existing 2 (now relabeled
midcarpal flex/dev) — resolving the measured DOF deficit (§3) directly. A parallel ulnocarpal
(TFCC) contact/ligament element from the proximal-row body to `ulna_r` would resolve the load-split
degeneracy (§4) the same way `docs/MECHANISM_PATELLOFEMORAL_FORCE.md`'s Pathway B (JAM/lenhart2015,
a free 6-DOF patella + separate ligament bundle) resolved the analogous quad:PT ratio-unity
problem there — a scoped, bounded build task with a precedented pattern already in this repo, not
a new methodology.
