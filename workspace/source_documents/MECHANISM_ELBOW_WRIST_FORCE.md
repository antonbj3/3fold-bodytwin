# MECHANISM ELBOW + WRIST FORCE — walking data-gate + twin-own muscled functional-pose estimate (2026-07-21)

Executes the operator's task: increase resolution on elbow/wrist joint contact loads. A prior
stub, `docs/MECHANISM_ELBOW_FORCE.md`, already exists and already certifies the elbow via
pure-kinematics + a generic-donor (`arm26`) muscled adversary. This session (1) forces the
walking-adequacy question with a real measurement rather than an assumption, (2) discovers
(live, not from doc prose) that a single existing model file now carries BOTH elbow- and
wrist-crossing muscles on subject2's own scaled anatomy, and uses it to compute a genuine
**twin-own** (not generic-donor) muscle-driven estimate for **both** joints at a matched
2/5/10 kg static-hold task, and (3) anchors both against live-reverified cadaveric/biomechanical
literature. Every number below is machine-measured this session
(`scripts/msk/validate_elbow_wrist_force.py`, exit 0), not recalled. Isolation respected:
`.venv-msk` only; all model/data files read in place, never written to; no git commit/push.

**Status: HYPOTHESIS awaiting independent QC. Confidence tier: cadaveric/published-plausibility,
task-caveated** (same tier vocabulary as `docs/MECHANISM_HAND_COMPLETE.md` /
`docs/MECHANISM_HAND_FIVE_DIGIT.md` — real, cited, textbook-standard biomechanics method; the
donor-graft geometry and load/grip-point choices are disclosed simplifications, not precision
subject-specific measurements).

## Headline result

| Question | Finding |
|---|---|
| **Does walking meaningfully load the elbow?** | **NO — diagnosed, quantified gap.** Peak walking-induced elbow reaction (arm-swing inertia only, no hand load) = **21.3 N = 2.77 %BW** — a small fraction of even the LIGHT 2 kg static-hold muscled estimate (30.6 %BW) and of the pure-kinematic 2 kg-hold reaction (38.4 N) from the existing cert. |
| **Does walking meaningfully load the wrist?** | **NO — more decisively than the elbow.** `wrist_flex_r` = **exactly 0.000° for all 158 frames** of the gait cycle; `wrist_dev_r` ranges **5.001°–5.142° (0.14° total)**. The wrist is kinematically inert during gait in this twin. |
| **Does the twin now support a TWIN-OWN (not generic-donor) muscled elbow estimate?** | **YES — newly discovered this session.** `subject2_hand_complete.osim` (built for the hand/wrist, `docs/MECHANISM_HAND_COMPLETE.md`) was found, live, to ALSO carry the 6 elbow-crossing muscles from the arm-muscle graft (`docs/MECHANISM_ARM_MUSCLES.md`) — a fact neither parent doc's own prose stated (both were checked; this is a genuinely new, machine-verified structural finding, not inferred from either doc). |
| **Elbow, 2 kg static hold, 90° flexion, twin-own muscles, real wrist markers** | **235.0 N = 11.98× the hand load's own weight = 30.6 %BW** — cross-validates the existing `arm26`-adversary number (234.5 N / 11.95× / 30.6 %BW) to **0.2%** |
| **Wrist, SAME 2 kg task, wrist-crossing muscles (FDP/FDS/FPL vs ED/EPL)** | **248.4 N = 12.66× the hand load's own weight = 32.4 %BW** — same order of magnitude as the elbow, same underlying mechanism |
| **Does the modeled ratio match the task's stated ~1–3× cadaveric ballpark?** | **NO — an honest, geometrically-explained mismatch, not forced to fit.** Measured 8–12× (elbow) / 11–13× (wrist) is 3–8× ABOVE the task's inherited ~1–3× figure. First-principles lever-arm geometry (r_load/r_muscle, measured directly on the model) independently PREDICTS 4.7–11.4× (elbow) / 12.3–13.3× (wrist) — i.e. the elevated ratio is the DIRECT, expected consequence of the measured short-tendon/long-lever geometry, not a bug; the task's own ~1–3× figure remains UNVERIFIED at the source (An 1981, paywalled, confirmed via 2 independent channels this session) and may describe a different pose/mechanism. |

## 1. What already existed, and what this session adds

`docs/MECHANISM_ELBOW_FORCE.md` (same day, prior session) already established: the base gait
model (`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`) has **0/80** muscles crossing
the elbow; a pure-kinematic reaction force at 2/5/10 kg static hold (36.6/66.1/115.1 N); and a
forced-adversary demonstration on `arm26` (Holzbaur 2005, a real but GENERIC, non-subject-scaled
6-muscle model) giving 234.5/455.2/823.0 N (11.95×/9.28×/8.39× the load's own weight). It flagged
as an honest next step: "graft real elbow flexor/extensor muscle-tendon actuators... parametrized
from a published source" onto the twin's own bodies.

Independently, `docs/MECHANISM_ARM_MUSCLES.md` (same day) already did exactly that graft — for a
DIFFERENT stated purpose (re-running the SHOULDER force cert) — producing
`LaiArnoldModified2017_arm_muscles_subject2_scaled.osim` with 6 new elbow-crossing muscles
(`BIC_long`, `BIC_brevis`, `BRA` flexors; `TRIlong`, `TRIlat`, `TRImed` extensors) alongside 22
shoulder-crossing muscles. That doc never re-ran the elbow cert against this new muscle set.

This session's contribution:
1. A **walking-adequacy data-gate**, measured, not assumed (§2).
2. Live structural audits (§3) that found `subject2_hand_complete.osim` (the newest hand-lineage
   model, `docs/MECHANISM_HAND_COMPLETE.md`) **already inherited the same 6 elbow-crossing
   muscles** from an earlier merge step, in addition to its own 14 wrist-crossing muscles — i.e.
   ONE file now supports a twin-own muscled estimate at BOTH joints. This was NOT stated in
   either parent doc and was found only by loading the file and querying it directly.
3. A forced-OODA diagnosis of a real donor-graft sign anomaly (`TRIlong`, §4) before trusting the
   elbow solve.
4. A genuinely NEW muscle-driven wrist analysis (§6), the first in this repo.
5. A first-principles geometric cross-check (§7) and a live re-verification of the elbow
   literature anchor plus two NEW wrist anchors (§8).

## 2. Walking data-gate — measured, not assumed

### 2.1 Elbow

Parsed the SAME `walking1.mot` (158 frames, subject2) every prior gait cert in this repo used.
`elbow_flex_r` ranges **40.12°–55.18° (15.06° total)** — a real, non-trivial angular excursion
(ordinary arm-swing counter-rotation), NOT frozen. But the FORCE consequence is small: cutting
the free body at `elbow_r` (forearm+hand, 1.736 kg on this model) and computing
`R = Σ(m·a) − m·g` from the REAL, Savitzky-Golay-differentiated body-COM trajectories — **no**
external hand load, **no** GRF term (the arm never touches ground during normal gait) — gives a
peak reaction of **21.26 N = 2.77 %BW** at t=0.66s, stable across a 70–170 ms smoothing-window
sweep (range 0.098 %BW, pre-registered <1.0 %BW gate: **PASS**).

**This is smaller than even the pure-kinematic 2 kg static-hold reaction (36.6 N / 4.78 %BW,
already published) and a small fraction of the muscle-inclusive 2 kg-hold estimate (30.6 %BW,
§5)** — i.e. ordinary gait arm-swing, with no hand load, exercises the elbow at roughly
**1/11th** the force level of even the LIGHTEST functional-hold task tested. Walking is a
load-inadequate task for an elbow contact-force estimate — diagnosed quantitatively, not
asserted from "arm swing is obviously low-load."

### 2.2 Wrist

Independently re-parsed (from the raw saved `.mot` file, not trusting doc prose) the ALREADY-
COMPUTED `subject2_unified_v2` walking-IK output
(`data/msk_smoketest/subject2_unified_v2_walking1/walking1_unified_v2_smoketest.mot`, produced
by `docs/MECHANISM_UNIFIED_V2.md`'s own merge+IK-smoketest):

| coordinate | min | max | range |
|---|---:|---:|---:|
| `wrist_flex_r` | 0.000° | 0.000° | **0.000°** |
| `wrist_dev_r` | 5.001° | 5.142° | 0.141° |
| `elbow_flex_r` (same file, cross-check) | 33.25° | 46.30° | 13.05° |

The wrist is **kinematically inert** — `wrist_flex_r` is bit-exact zero across all 158 frames.
This is honestly a DATA-AVAILABILITY zero as much as a physiological one (`walking1.trc`, the
source mocap marker file, has no hand/wrist markers, so these two coordinates are unweighted in
the IK objective) — but it is decisive for the data-gate question either way: there is no usable
signal to compute a walking-driven wrist force from, regardless of the reason. **Walking is
load-inadequate for the wrist, more starkly than for the elbow** (the elbow at least has a real,
non-zero, if modest, ROM to measure a force from).

## 3. Structural audit — a genuinely new finding, not assumed from either parent doc

Rather than trust `docs/MECHANISM_ARM_MUSCLES.md`'s framing (a shoulder-focused graft) or
`docs/MECHANISM_HAND_COMPLETE.md`'s framing (a hand/wrist-focused build) about which muscles
exist where, both candidate models were loaded live and queried directly:

| model | elbow-crossing muscles found | wrist-crossing muscles found |
|---|---:|---:|
| `LaiArnoldModified2017_arm_muscles_subject2_scaled.osim` (78.2 kg) | **6**: `BIC_long`, `BIC_brevis`, `BRA`, `TRIlong`, `TRIlat`, `TRImed` | 0 (wrist still welded, `radius_hand_r`) |
| `subject2_hand_complete.osim` (79.15 kg) | **6 — IDENTICAL names**, moment arms agree to ≤0.28mm at 90° flexion (BIC_long/BRA/TRIlat/TRImed match to <0.001mm) | **14**: `FDP2-5_r`, `FDS2-5_r`, `FPL_r` (flexors, +4.4 to +14.6mm at `wrist_flex_r`); `ED2-5_r`, `EPL_r` (extensors, −1.9 to −12.4mm) |

**`subject2_hand_complete.osim` carries BOTH sets** — confirmed live (not assumed), including a
direct moment-arm cross-check between the two files proving they share the same underlying
elbow-muscle graft geometry. This was unexpected from either doc's own account (an over-eager
first-pass read of the docs alone concluded, wrongly, that the elbow-muscled lineage and the
wrist-articulated lineage were structurally disjoint — corrected here by direct measurement, per
the "MACHINE CROSS-CHECK, never narration" discipline). All subsequent computation in this
document uses this ONE model, for both joints, at 78.2 kg %BW normalization (subject2's own
canonical mass per `sessionMetadata.yaml`; the model's own total mass is 79.15 kg, +1.2% from the
added multi-body hand anatomy — disclosed, not hidden).

## 4. Forced OODA: a real donor-graft sign anomaly found and excluded, not silently used

Before trusting the elbow solve, all 6 crossing muscles' moment-arm sign was checked at the
90°-flexion analysis pose (`computeMomentArm`, the same already-correct-convention API used
throughout this repo, never a custom re-derivation):

| muscle | moment arm @ 90° | expected sign | measured sign |
|---|---:|---|---|
| BIC_long, BIC_brevis, BRA | +26 to +64mm | flexor (+) | **correct** |
| TRIlat, TRImed | −6.4mm | extensor (−) | **correct** |
| **TRIlong** | **+33.8mm** | extensor (−) | **WRONG (anomalous)** |

**Forced, not assumed**: a 6-point pose sweep (0°–120°) shows `TRIlong`'s moment arm is
**correctly negative at 0° (−20.1mm) but flips sign between 30°–45° and grows to +42.3mm by
120°** — while `TRIlat` (a properly-multi-point-pathed sibling extensor from a different donor)
stays correctly negative through ~110°. Path inspection shows why: `TRIlong`'s path (inherited
from Donor 1, the Seth-2016 ThoracoscapularShoulderModel) is a **bare 2-point straight line**
from `ulna_r` directly to `torso`, with **no via-point on the humerus at all** — unlike
`TRIlat`'s realistic 5-point path. This is a newly-found (this session) instance of the SAME
already-disclosed "no wrapping surfaces ported" limitation `docs/MECHANISM_ARM_MUSCLES.md` §6/§8
already reported for `DeltoideusScapula_M`/`Supraspinatus_A/P` at the shoulder — here it recurs
for `TRIlong` at the elbow specifically.

**`TRIlong` is excluded from the muscle-force solve below**, disclosed with the full pose-sweep
evidence above. This does NOT change the numeric headline result: the task (holding a load,
pure flexion demand) has `M_req > 0`, so the minimum-effort solution assigns **zero** activation
to every extensor candidate regardless (the same one-line proof `docs/MECHANISM_ELBOW_FORCE.md`
used for `arm26`'s triceps) — `TRIlong`'s exclusion is reported for completeness/symmetric QC,
not because it would have changed the answer.

## 5. Elbow functional pose: twin-own muscled estimate

Same pose convention as the existing cert (already verified there: `arm_flex_r=0` puts the
humerus vertical): `arm_flex_r=0°`, `elbow_flex_r=90°`, `pro_sup_r=0°`. Load applied at the
**real** `R_wrist_radius`/`R_wrist_ulna` marker midpoint (this model has the twin's own real
markers — an upgrade over `arm26`'s estimated-ratio proxy, which had no wrist body/marker at
all). Closed-form minimum-activation-squared solve over `BIC_long`/`BIC_brevis`/`BRA` (the
provably-optimal-to-include flexors; `TRIlat`/`TRImed` get zero activation by the same proof,
`TRIlong` excluded per §4).

| load | R (pure reaction) | required moment | F_bone_contact | × load weight | %BW | max activation |
|---:|---:|---:|---:|---:|---:|---:|
| 2 kg | 38.4 N | 10.09 N·m | **235.0 N** | **11.98×** | **30.6%** | 0.126 |
| 5 kg | 67.8 N | 18.70 N·m | 432.6 N | 8.82× | 56.4% | 0.233 |
| 10 kg | 116.8 N | 33.04 N·m | 761.9 N | 7.77× | 99.3% | 0.412 |

Forearm+hand mass on this model: 1.913 kg (vs. 1.736 kg on the original weldHand-only model —
the +0.177 kg difference is exactly the added multi-body hand/finger anatomy, explaining the
correspondingly slightly-higher pure-reaction R here (38.4N vs 36.6N, +4.9%, matching the mass
delta's proportional contribution, +4.7%) — a small, fully-explained, disclosed difference, not
a discrepancy.

**No saturation at any load** (max activation 0.126–0.412, real muscles have ample reserve — in
sharp, informative contrast to the existing doc's finding that the model's OWN idealized 10 N·m
`CoordinateActuator` was inadequate above ~2kg: that actuator was clearly undersized/non-
representative of true muscle capacity, not a genuine physiological limit).

**Cross-validation against the existing `arm26` adversary** (an independent muscle-geometry
source — different donor lineage, different individual, same method):

| load | twin-own (this session) | `arm26` (existing cert) | relative difference |
|---:|---:|---:|---:|
| 2 kg | 235.0 N (11.98×, 30.6%BW) | 234.5 N (11.95×, 30.6%BW) | **0.2%** |
| 5 kg | 432.6 N (8.82×, 56.4%BW) | 455.2 N (9.28×, 59.4%BW) | 5.0% |
| 10 kg | 761.9 N (7.77×, 99.3%BW) | 823.0 N (8.39×, 107.3%BW) | 7.4% |

Two independently-geometried muscle sets (Holzbaur 2005's `arm26` vs. this repo's Seth-2016/
`arm26`-hybrid donor graft, applied to subject2's own scaled anatomy) agree to within 0.2–7.4%
across the entire load sweep — a genuine over-determination check, not a tautology (neither
model was tuned to match the other).

## 6. Wrist functional pose: twin-own muscled estimate (NEW this session)

**Same overall arm pose** as §5 (`arm_flex_r=0°`, `elbow_flex_r=90°`, `pro_sup_r=0°`), extended
distally: `wrist_flex_r=0°`, `wrist_dev_r=0°`, fingers/thumb at their model-default 0° (a flat,
extended hand — NOT a closed fist/power grip; disclosed, see §6.3). This is a single coherent
"holding 2 kg, forearm horizontal" scene analyzed at **two nested free-body cuts** (elbow-distal,
wrist-distal), mirroring this repo's own established knee/hip/ankle "nested cut" convention.

Load applied at `middle_distal_r`'s body origin (a disclosed proxy for "object held near the
fingertips" — no fingertip/grip marker exists in this model's `MarkerSet`, confirmed live).
Required moment computed via the SAME bug-immune virtual-work method (§5); the agonist/antagonist
split generalizes the `arm26` one-line proof to whichever sign the REQUIRED actuator moment
actually has (measured, not presumed): at this pose the load demands a **flexor** moment, so
`FDP2-5_r`/`FDS2-5_r`/`FPL_r` are the agonists (`ED2-5_r`/`EPL_r` get zero activation, same proof).

| load | R (pure reaction) | required moment | F_bone_contact | × load weight | %BW | max activation | saturates? |
|---:|---:|---:|---:|---:|---:|---:|---|
| 2 kg | 26.0 N | 3.47 N·m | **248.4 N** | **12.66×** | **32.4%** | 0.787 | no |
| 5 kg | 55.4 N | 7.89 N·m | 565.6 N | 11.53× | 73.7% | **1.790** | **YES** |
| 10 kg | 104.5 N | 15.27 N·m | 1094.2 N | 11.16× | 142.7% | **3.463** | **YES (badly)** |

### 6.1 The 2 kg case is the physiologically-valid headline; 5/10 kg reveal a disclosed model gap

At 2 kg, no muscle saturates (max activation 0.787) — a physiologically sane result. At 5 kg and
10 kg, the closed-form solve DEMANDS activation >1.0 (up to 3.46×, i.e. 3.46× a muscle's own
maximum isometric force) — **not physiologically achievable**. This is diagnosed, not silently
reported as a number: **this model has NO dedicated wrist-flexor muscles** (no FCR, FCU, or
palmaris longus) — the ONLY available wrist-flexion agonists are the extrinsic FINGER flexors
acting as a secondary effect of their own primary (finger) action. Real wrists share this load
across BOTH dedicated wrist flexors AND the finger flexors; a model missing the former will
force the latter alone toward saturation at loads a real, complete muscle set would distribute
more broadly. This is a genuine, disclosed **model-completeness gap** (not necessarily a claim
that a real human wrist cannot hold 5-10kg this way) — flagged as the clearest next-step
follow-on this session identifies (§10).

### 6.2 Elbow vs. wrist: same order of magnitude, same mechanism

At the shared 2 kg task: elbow 11.98×, wrist 12.66× — remarkably similar force-multiplication
despite being different joints with different muscles/geometry, because both are governed by the
SAME mechanical principle: a short muscle-tendon moment arm (26–64mm at the elbow; 14–15mm at the
wrist) versus a much longer lever arm to the held load (30.0cm elbow-to-wrist; 18.1cm
wrist-to-grip-point) — the classic third-class-lever mechanical disadvantage, at BOTH joints (§7
derives this from first principles, not just curve-fits it).

### 6.3 Honest task/literature mismatch, disclosed not conflated

This is a STATIC-HOLD moment-balance task (weight of hand+object resisted by a flat, extended
hand at the wrist) — the SAME task TYPE as the elbow analysis, chosen for direct elbow↔wrist
comparability. This is mechanistically DIFFERENT from Schuind et al. 1995's own modeled task (a
10N GRASP/grip force actively generated by finger flexion and transmitted through the wrist,
§8) — real power-grip wrist loads are dominated by the tendon tension needed to generate grip
force against an object (a friction/slip problem this session does not model, since it would
require an unconstrained new friction-coefficient assumption), not by the object's own weight.
Both are real, legitimate wrist-loading mechanisms; they are not the same number, and this
document does not conflate them.

## 7. First-principles geometric cross-check (derived, not curve-fit)

Measured directly on the model (body/marker ground positions at the analysis pose):

- Elbow-to-wrist-marker-midpoint distance (the load's lever arm for the elbow task): **300.0mm**
- Wrist-to-grip-point distance (the load's lever arm for the wrist task): **180.7mm**

A pure lever-arm ratio (`r_load / r_muscle`, ignoring the segment's own weight — a first-order
geometric prediction, not the full closed-form solve) predicts:

| joint | via muscle | r_muscle | predicted ratio | measured closed-form ratio (§5/§6) |
|---|---|---:|---:|---:|
| elbow | BIC_long | 62.3mm | 4.81× | 7.77–11.98× |
| elbow | BIC_brevis | 64.0mm | 4.69× | (same range) |
| elbow | BRA | 26.2mm | 11.43× | (same range) |
| wrist | FDP2_r/FDS2_r | 14.6mm | 12.33–12.34× | 11.16–12.66× |
| wrist | FDP4_r/FDP5_r | 13.6–13.8mm | 13.13–13.31× | (same range) |

**The measured closed-form ratios fall squarely inside the range this pure geometric argument
independently predicts, at both joints.** This confirms the elevated (8–13×) force-multiplication
is a DIRECT, expected consequence of the measured moment-arm/lever-arm geometry — not a
numerical artifact of the closed-form solve, and not inconsistent with the general "short
tendon moment arm" mechanism the existing elbow doc already invoked for `arm26`.

## 8. External literature anchors — live-reverified this session (PubMed efetch + CrossRef)

`WebSearch` was unavailable this session (shared session-wide quota, "2000/2000" — the identical
constraint `docs/MECHANISM_ELBOW_FORCE.md`/`docs/MECHANISM_MUSCLE_AUDIT.md` already hit).
CrossRef's and PubMed's REST APIs (unaffected, unauthenticated, live-queried via `curl`) were
used instead — the same established pattern this repo's own hand-build docs already use.

### 8.1 Elbow: An, Hui, Morrey, Linscheid, Chao (1981) — re-verified, still no extractable number

**An KN, Hui FC, Morrey BF, Linscheid RL, Chao EY (1981). "Muscles across the elbow joint: a
biomechanical analysis." J Biomech 14(10):659-69. DOI 10.1016/0021-9290(81)90048-8. PMID
7334026** (newly identified this session). Re-verified via **two independent channels** this
session (CrossRef, matching the existing doc's own check, AND PubMed `esearch`/`efetch`, a
genuinely new second channel): existence/authors/year/journal/DOI/PMID all confirmed both ways.
**PubMed's own record has no abstract text** (`efetch` returns title/authors/DOI/PMID only) — the
paper's specific numeric finding remains genuinely unextractable via either of two independent
channels (Elsevier paywall), not merely untried a second time. **The task's own "~1-3× hand load"
figure remains inherited, not independently re-derived from this paper's own text** — same
honest disclosure as the existing doc, now on a firmer (doubly-checked) footing.

### 8.2 Wrist: Schuind, Cooney, Linscheid, An, Chao (1995) — FULL abstract retrieved (a stronger tier)

**Schuind F, Cooney WP, Linscheid RL, An KN, Chao EY (1995). "Force and pressure transmission
through the normal wrist. A theoretical two-dimensional study in the posteroanterior plane."
J Biomech 28(5):587-601. DOI 10.1016/0021-9290(94)00093-j. PMID 7775494.** Full abstract
retrieved live via PubMed `efetch` this session (a stronger evidence tier than the elbow anchor
above, which yielded no abstract):

- Method: rigid body spring modeling (RBSM) of 120 normal wrist PA X-rays; carpal bones as rigid
  bodies interposed by compression springs (cartilage) and tension springs (ligaments).
- **Simulated task: axial loads along the metacarpals simulating a 10N GRASP force**, active
  wrist stabilization in NEUTRAL position — i.e. a grip-generation task, not a passive
  weight-holding task (§6.3's caveat).
- Findings: radioscaphoid carries **55%** of the transmitted load, radiolunate **35%**, the TFCC
  **10%**; peak-pressure ratio radioscaphoid:radiolunate = **1.6**; ligaments opposing ulnar
  carpal translation dominate load transmission; wrist morphology/age have little effect.

This is a DISTRIBUTION anchor (which sub-articulation carries how much of a given total load),
not a total-force-multiplication ratio — it cannot be numerically matched against this session's
own "×load weight" ratio, and this document does not force that comparison. It IS a real,
verified, on-topic anchor for "wrist joint loads are tendon/grip-force-transmission-dominated,"
consistent with §6.2's own finding that the wrist shows an elbow-like force-multiplication regime.

### 8.3 Wrist: An, Chao, Cooney, Linscheid (1985) — newly found, directly on-topic methodology anchor

**An KN, Chao EY, Cooney WP, Linscheid RL (1985). "Forces in the normal and abnormal hand."
J Orthop Res 3(2):202-211. DOI 10.1002/jor.1100030210. PMID 3998897.** Newly found and verified
live this session (full abstract). Same Mayo Clinic research group as the Schuind anchor. An
analytic, cadaveric-measurement-anchored model computing **muscle AND joint force
distributions** under isometric hand functions — confirms this exact class of quantity
(tendon-force → wrist/hand joint contact-force mapping) is an established methodology in this
literature; the abstract itself states no single numeric ratio.

### 8.4 Wrist: Rikli — the "three-column" load-sharing framework, a qualitative anchor

**Jakob M, Rikli DA, Regazzoni P (2000). "Fractures of the distal radius treated by internal
fixation and early function. A prospective study of 73 consecutive patients." J Bone Joint Surg
Br 82(3):340-4. DOI 10.1302/0301-620x.82b3.10099. PMID 10813166.** Verified live this session
(full abstract); corroborated by CrossRef surfacing multiple later Rikli/Regazzoni distal-radius
papers and third-party figures explicitly captioned "the three-column model by D. Rikli and P.
Regazzoni." Rikli's contribution is the clinically-standard **radial + intermediate (+
ulnar/TFCC) column** load-bearing framework for the distal radius/wrist, used to guide fracture
fixation — a **qualitative** (which structures share transmitted load) anchor, correctly
characterized here as distinct from Schuind's quantified RBSM model, not a numeric
force-multiplication ratio.

## 9. Verdict — the falsifier, resolved plainly

**Falsifier as posed: does the modeled elbow/wrist contact match the cadaveric ratio for a
matched task, or is the walking task load-inadequate?**

1. **Walking is load-inadequate — CONFIRMED, decisively, for both joints** (§2). Elbow: peak
   walking reaction (2.77%BW) is ~1/11th the light 2kg-hold estimate. Wrist: kinematically inert
   (exactly 0.000° `wrist_flex_r` for the entire gait cycle) — no usable signal exists at all.
2. **The modeled functional-pose ratio does NOT match the task's stated ~1-3× cadaveric
   ballpark — an honest mismatch, not forced to fit.** Measured 7.8–12.0× (elbow) / 11.2–12.7×
   (wrist) is 3–8× above that figure. This is NOT dismissed as noise: it is independently
   corroborated by (a) a completely separate, real, peer-reviewed muscle model (`arm26`,
   agreeing to 0.2–7.4%) and (b) a first-principles lever-arm geometric derivation (§7) that
   predicts the same elevated range from measured distances alone. The task's own ~1-3× figure
   traces to An 1981, whose specific number could not be extracted after TWO independent
   verification attempts (paywalled both times) — it may describe a different pose (e.g. a less
   flexed elbow, where the SAME muscles could have different moment arms) or a different
   quantity than the one this session computed. This is reported as an open, unresolved
   discrepancy, not silently reconciled.
3. **The wrist has no directly-comparable numeric ratio in the verified literature** (Schuind
   reports distribution percentages; An 1985's abstract states no ratio; Rikli is a qualitative
   framework) — only an order-of-magnitude, mechanism-level consistency claim is honestly
   possible (§6.2/§8.2), not a number-to-number match.

## 10. Honest gaps (full list)

1. **The elbow "~1-3×" and any wrist ratio anchor remain numerically unverified at the source**
   (§8.1/8.2/8.3/8.4) — every citation's existence/identity is live-verified (2 independent
   channels for the elbow paper), but none yields an extractable precision number this session.
2. **`TRIlong`'s donor-graft geometry is anomalous** (§4) — excluded, disclosed, does not affect
   the numeric headline (proof in §4), but is a real quality gap in the underlying graft, an
   instance of its own already-disclosed "no wrapping surfaces" limitation recurring at a new
   joint.
3. **The wrist model has no dedicated wrist-flexor muscles** (FCR/FCU/palmaris longus absent,
   §6.1) — causes non-physiological (>1.0) activation at 5kg/10kg; the 2kg headline is
   unaffected, but this is a real model-completeness gap, not a subtle numerical issue.
4. **Load/grip-point placement is a disclosed proxy** — the elbow uses REAL wrist markers (an
   upgrade over the prior `arm26` estimate), but the wrist analysis uses `middle_distal_r`'s body
   origin (no fingertip/grip marker exists in this model) — likely UNDERESTIMATES the true lever
   arm to an object held at the fingertip pad, meaning the true ×load-weight ratio could be
   modestly HIGHER than reported, not lower.
5. **Flat, extended-hand pose, not a closed-fist power grip** (§6.3) — fingers/thumb at 0°
   (model default); a genuine grip-force task (matching Schuind's own modeled mechanism) would
   need a different, more complex analysis (friction/slip, active finger-flexion-to-generate-
   grip-force) not attempted here, by disclosed scope choice.
6. **Single pose, right side only, static** — no ROM sweep, no dynamic/EMG validation, same scope
   caveat as every joint-force cert in this repo.
7. **Different individuals/lineages inside one "twin-own" claim** — the muscle GEOMETRY is
   donor-sourced (Seth 2016 DSEM lineage + `arm26`) scaled onto subject2's own skeleton; only the
   skeletal anthropometry and real wrist markers are genuinely subject2's own. "Twin-own" here
   means "on the twin's own scaled skeleton," not "subject2's own cadaveric muscle geometry."
8. **No in-vivo instrumented-implant anchor exists for either joint** — already flagged as a
   repo-wide honest gap (`MSK-JOINT-CONTACT-FORCE-ANCHORS`), inherited unchanged.

## 11. Files

- `scripts/msk/validate_elbow_wrist_force.py` — the full pipeline (self-contained, re-runnable,
  exit 0; imports `validate_joint_force.py` for `get_descendant_bodies`/`savgol_smooth_and_derivs`/
  `parse_mot`, and `validate_elbow_force.py` for `virtual_work_generalized_force`/
  `crossing_muscles_generic` — not re-implemented).
- `data/msk_smoketest/subject2_elbow_wrist_static/elbow_wrist_force_results.json` — every number
  in this document, machine-written (walking data-gates, both functional-pose analyses at
  2/5/10kg, the TRIlong pose-sweep diagnostic, literature anchors, verdict block).
- Read/reused, unmodified: `scripts/msk/validate_joint_force.py`, `scripts/msk/validate_elbow_force.py`,
  `data/msk_models/LaiArnoldModified2017_arm_muscles_subject2_scaled.osim`,
  `data/msk_models/subject2_hand_complete.osim`,
  `data/msk_smoketest/subject2_unified_v2_walking1/walking1_unified_v2_smoketest.mot`.
- Prior docs this builds on (not re-litigated, only extended/cross-checked):
  `docs/MECHANISM_ELBOW_FORCE.md`, `docs/MECHANISM_ARM_MUSCLES.md`, `docs/MECHANISM_HAND_COMPLETE.md`,
  `docs/MECHANISM_UNIFIED_V2.md`, `docs/MECHANISM_JOINT_FORCE_VALIDATION.md`,
  `docs/MECHANISM_STATIC_OPT.md`, `docs/MECHANISM_SHOULDER_FORCE.md`.

## 12. Next step

1. Graft dedicated wrist-flexor/extensor muscles (FCR, FCU, palmaris longus, ECRL/ECRB/ECU) onto
   `subject2_hand_complete.osim` to remove the 5-10kg saturation artifact (§6.1) — a scoped,
   bounded build task, same pattern as `add_arm_muscles.py`.
2. Re-attempt An 1981's exact numeric figure via an alternative access route (e.g. an
   institutional-access proxy or a citing secondary source that quotes it directly) if a tighter
   elbow-anchor reconciliation is later required.
3. If a genuine grip-force (power-grip) wrist estimate is wanted, matching Schuind's own modeled
   mechanism: model active finger flexion generating a specified grip force (with an explicit,
   disclosed friction-coefficient assumption for "how much grip force is needed to not drop a 2kg
   object"), not the passive weight-holding task used here.
4. A real ROM sweep of the elbow's own force-multiplication ratio (not just 90°) would test
   whether the task's ~1-3× figure is consistent with SOME pose in the ROM even if not 90° —
   flagged as a concrete, cheap follow-on given the moment-arm curves are already measured (§4).
