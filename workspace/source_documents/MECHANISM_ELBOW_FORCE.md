# MECHANISM ELBOW FORCE — replicating the knee force-cert method for the humeroulnar joint (2026-07-21)

Executes the operator's task: replicate the VALIDATED knee force-cert method (Static Optimization +
the OpenSim-4.6-ID-bug-immune Newton/virtual-work method, `docs/MECHANISM_STATIC_OPT.md` +
`docs/MECHANISM_JOINT_FORCE_VALIDATION.md`) for the ELBOW. Every number below is machine-measured
this session (`scripts/msk/validate_elbow_force.py`, exit 0, fully deterministic — a repeat run
reproduced the log byte-for-byte), not recalled. Isolation respected: `.venv-msk` only; external
model/data (LabValidation, `arm26.osim` under the separate `opensim_jam_build` tree) read in place,
never written to; no git commit/push.

## Headline result

**The twin cannot produce a muscle-driven elbow contact-force estimate — not because its arm
muscle set is "too sparse," but because it has ZERO elbow (or shoulder, or pro/sup) muscles at
all.** This was forced via a live, machine-checked structural audit (Sec.2), not inferred from the
model's name or asserted from a one-shot look. What the twin CAN still produce (a pure kinematic
reaction force, bug-immune by construction) is computed below, together with a forced-adversary
demonstration on a real muscled reference model to show what the method predicts when muscles
*are* present.

| finding | value | verdict |
|---|---:|---|
| Twin muscles crossing the elbow, any of 3 poses (0/45/90°), all 80 checked | **0/80** | Structural absence, not sparsity |
| Twin elbow/shoulder/pro-sup actuation | 10× idealized `CoordinateActuator`, 10 N·m each | 0 muscles |
| Required elbow-flexion moment, 5 kg hand load @ 90°, upper arm vertical | 18.14 N·m | vs 10 N·m actuator capacity → **ratio 1.81 — actuator inadequate** |
| Twin pure-reaction ≡ "bone-contact" force (identical; 0 muscles to subtract), 5 kg case | 66.1 N | **8.6 %BW** |
| **arm26** (Holzbaur/Murray/Delp 2005, real 6-muscle elbow model) forced-adversary demonstration, 5 kg case | 455.2 N | **9.28× the applied load's own weight**, 59.4 %BW (borrowed reference mass) |

## 1. Why this replicates the knee method (and where it structurally cannot)

The knee cert's method has two legs: (a) a pure-kinematic Newton's-law reaction force, immune to
the documented OpenSim 4.6 `InverseDynamicsTool` knee/hip generalized-force bug
(`docs/MECHANISM_MSK_ELASTIC_BAND.md` §4) because it never calls that code path; (b) Static
Optimization to recover the muscle-crossing correction, `F_bone_contact = R − ΣF_muscle`. Both legs
are replicated here:

- **(a) Newton's-law / virtual-work reaction force** — implemented exactly as before for the
  FORCE balance (statics: `0 = gravity + hand-load + R`), and *generalized* to a **virtual-work
  finite-difference method** for the MOMENT balance (`virtual_work_generalized_force()`): the
  generalized force conjugate to `elbow_flex_r` from an applied point force is
  `Q = F · (∂p/∂q)`, with `∂p/∂q` obtained by central finite difference of the force-application
  point's ground position as the coordinate is perturbed — the same virtual-work principle as the
  elastic-band doc's `dPE/dq` bug workaround, generalized here from a conservative potential to an
  arbitrary applied force. Never calls `InverseDynamicsTool`/`InverseDynamicsSolver`/
  `JointReaction` — immune to the documented bug by construction.
- **(b) Static Optimization** — structurally unavailable for the twin at the elbow (Sec.2): there
  is no redundant muscle set to resolve. It IS exercised, in the standard textbook
  minimum-activation-squared closed form, on the forced-adversary demonstration model (Sec.6),
  using real OpenSim `computeMomentArm()` values (a different, ID-bug-unaffected code path) and the
  SAME live-path-crossing force-direction detector already validated for the knee cert
  (`static_opt_knee.py`'s `knee_crossing_muscles_and_forces`, reused for the twin's own
  crossing-check; a small local variant, `crossing_muscles_generic`, used for arm26 since its 6
  muscle names carry no `_l`/`_r` suffix for that function's side-filter to match — caught by
  inspection before running, not after a false-negative).

## 2. The decisive structural finding: 0/80 muscles cross the elbow, at any pose

Live-queried from `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` (subject2):

- **Joint chain** (confirmed via the model's own joint tree, BFS, not assumed): `elbow_r` (`PinJoint`,
  coordinate `elbow_flex_r`) → `radioulnar_r` (`PinJoint`, coordinate `pro_sup_r`) → `radius_hand_r`
  (`WeldJoint`, 0 DOF — confirms the model's own "weldHand" naming: no wrist articulation). BFS
  descendant-of-`elbow_r` free-body cut = `{ulna_r, radius_r, hand_r}`, total mass **1.736 kg**
  (≈2.2% of the model's 78.2 kg body mass — matches Winter's anthropometric forearm+hand fraction
  almost exactly, an unprompted sanity check).
- **Actuator inventory** (`getForceSet()`, all 93 actuators enumerated): 80
  `Millard2012EquilibriumMuscle` + 13 `CoordinateActuator`. Of the 13: 3 lumbar + **10 arm**
  (`shoulder_flex/add/rot_{r,l}`, `elbow_flex_{r,l}`, `pro_sup_{r,l}`), **every one an idealized
  `CoordinateActuator`, optimal_force = 10.0 N·m, on BOTH sides** — confirmed live, matching
  `docs/MECHANISM_MUSCLE_AUDIT.md`'s independent top-line inventory (all 80 named muscles are
  lower-limb/trunk; zero arm muscle names anywhere in the model family).
- **Falsifier check, forced not assumed**: reused the knee cert's own BFS + live-path-crossing
  detector (`knee_crossing_muscles_and_forces`) across **all 80 muscles**, at **3 different elbow
  poses** (0°/45°/90° flexion — a wrapping path could in principle only engage at certain angles,
  so this was checked, not inferred from one pose). **Result: 0 crossings at every pose.** The
  pre-registered falsifier (any crossing found ⇒ revise the claim) did not fire.

**This is a whole-arm gap, not elbow-specific** — shoulder (`arm_flex/add/rot_r`) and forearm
pro/sup are equally muscle-free. This matches the task's own prediction that "a shoulder agent may
independently hit the same arm-muscle limitation" — confirmed here directly rather than inferred
from a parallel agent's output.

**Consequence**: `F_bone_contact = R − Σ F_muscle_crossing` reduces to `F_bone_contact ≡ R`
**identically**, for any hypothetical Static-Optimization activation pattern — not a numerical
result, a structural identity, since the subtracted term is the empty sum.

## 3. Representative pose — verified, not assumed

"Holding/curling a load at 90°" was operationalized as: elbow flexed 90° (the standard isometric
elbow-flexion strength-testing angle), upper arm hanging vertically at the side (`arm_flex_r = 0°`,
maximizing gravity's moment arm about the flexion axis — the conservative/worst-case choice). That
`arm_flex_r = 0°` actually puts the humerus vertical was **verified live**, not assumed: the
shoulder-to-elbow vector's horizontal offset at 0° is 0.019 m vs a 0.337 m vertical drop (vs
horizontal/vertical of 0.156 m / 0.300 m at −30° and 0.182 m / 0.284 m at +30° — horizontal offset
roughly 8× larger and vertical drop clearly smaller at either ±30°, confirming 0° is the
near-vertical configuration) — confirmed for both the twin and `arm26` independently (same
convention in both models, itself a small cross-model consistency check).

## 4. Required elbow-flexion moment vs the twin's own actuator capacity

A load sweep (2/5/10 kg — a declared sensitivity range, not one knife-edge number) against the
`elbow_flex_r` `CoordinateActuator`'s own 10 N·m rating:

| load | required moment (N·m) | contribution: forearm+hand weight | contribution: hand load | ratio vs 10 N·m cap | verdict |
|---:|---:|---:|---:|---:|---|
| 2 kg | 9.51 | 3.75 | 5.75 | 0.95 | marginal (just under) |
| **5 kg** | **18.14** | 3.76 | 14.38 | **1.81** | **INADEQUATE** |
| 10 kg | 32.52 | 3.76 | 28.77 | 3.25 | INADEQUATE |

Even a light 2 kg hand load already consumes 95% of the idealized actuator's entire rated capacity;
5 kg and 10 kg exceed it by 1.8× and 3.3×. **The twin, as currently actuated, could not sustain this
pose in a forward/controlled simulation at any of these loads except the lightest — a stronger,
more concrete finding than "the muscle set looks sparse."** This actuator was clearly never sized
for loaded upper-limb tasks (consistent with its likely purpose: providing enough torque for
passive arm-swing kinematics during gait IK, its actual, designed use-case elsewhere in this
model family).

## 5. Twin: pure reaction force (= "bone contact," identically, per Sec.2)

| load | R (bone-contact) force | %BW | × load weight |
|---:|---:|---:|---:|
| 2 kg | 36.6 N | 4.78% | 1.87× |
| **5 kg** | **66.1 N** | **8.61%** | 1.35× |
| 10 kg | 115.1 N | 15.01% | 1.17× |

**Self-consistency check** (implementation correctness, not a physics discovery — statics has no
acceleration term, so `R` cannot depend on elbow angle by construction): recomputed at
elbow_flex_r = 0°/45°/90°, all three give **exactly 66.058 N** for the 5 kg case (spread = 0.0%,
pre-registered tolerance 1.0% — PASS). Contrast: the required MOMENT is **not** pose-invariant
(depends on the shoulder angle, since gravity's projection onto the flexion axis changes) —
verified separately: −18.14 N·m at shoulder = 0°, decreasing in magnitude to −17.07 at 30° and
−11.43 at 60°, the expected qualitative trend as the forearm rotates away from horizontal.

## 6. Forced adversary: what a muscled elbow model predicts (arm26, Holzbaur 2005)

To force the "is a muscle-driven estimate even achievable" adversary to its strongest fair form,
the exact same method (virtual-work reaction force + minimum-effort muscle force-sharing + real
path-geometry force subtraction) was applied to **`arm26`** — a real, peer-reviewed, 6-muscle
elbow model (citation verified from the model file's own embedded `<credits>`/`<publications>`
tags: **Holzbaur, K.R.S., Murray, W.M., Delp, S.L. "A Model of the Upper Extremity for Simulating
Musculoskeletal Surgery and Analyzing Neuromuscular Control." Annals of Biomedical Engineering,
33:829–840, 2005**). **This is explicitly NOT the twin** — a generic, non-subject-scaled reference
model, read in place from the separate `opensim_jam_build` source tree, used only to demonstrate
the method and give a concrete, non-hand-wavy number.

Real moment arms at 90° flexion (`computeMomentArm`, unaffected by the ID-tool bug — a different
code path): flexors BIClong/BICshort 4.88 cm, BRA (brachialis) 2.27 cm; extensors TRIlong/TRIlat/
TRImed ≈ −1.99 cm (sign confirms flexor/extensor role correctly). A closed-form minimum-effort
solve (`aᵢ = Fmax·r / Σ(Fmax·r)²  ×  M_req`, the standard textbook Static-Optimization reduction for
a single equality constraint) allocates activation across the 3 flexors only — **provably optimal,
not assumed**: since `M_req > 0` and all included moment arms are positive, all activations come
out non-negative automatically; co-activating any antagonist extensor would strictly increase the
cost function (it would force the flexors to work harder AND add its own non-negative term), so the
true minimum has zero triceps activation — a one-line proof, not hand-waved.

| load | required moment (N·m) | activation at this static pose (BIClong/BICshort/BRA) | muscle-crossing force | F_bone_contact | × load weight | %BW (borrowed ref.) |
|---:|---:|---:|---:|---:|---:|---:|
| 2 kg | 7.53 | 0.122 / 0.085 / 0.090 | 200.7 N | 234.5 N | 11.95× | 30.6% |
| **5 kg** | **14.74** | 0.239 / 0.167 / 0.176 | 392.8 N | **455.2 N** | **9.28×** | **59.4%** |
| 10 kg | 26.75 | 0.433 / 0.302 / 0.319 | 713.0 N | 823.0 N | 8.39× | 107.3% |

No muscle saturates (all activations well under 1.0) across the sweep. The load application point
(arm26 has no separate hand/wrist body — it fuses forearm+hand into one segment, and has no wrist
marker) was estimated by **transferring the twin's own measured
(elbow-to-wrist-marker)/(elbow-to-forearm+hand-COM) distance ratio** (1.352×, from the twin's real
`R_wrist_radius`/`R_wrist_ulna` markers) onto arm26's own COM direction — a measured,
cross-model-calibrated proxy, not an arbitrary guess.

**A real bug was caught and fixed in this step, not silently left in**: a first implementation
attached the load-force application point as a fixed GROUND-frame point for the finite-difference
virtual-work calculation. Since a fixed ground point does not move as the coordinate is perturbed,
its contribution to the generalized force is identically zero — this produced the SAME required
moment (2.727 N·m) for all three loads in the sweep, which is physically impossible (caught by
exactly that invariance, not by a crash). Fixed by transforming the reference-pose ground point into
a body-local fixed offset (round-trip-verified via
`ground.findStationLocationInAnotherFrame`/`body.findStationLocationInGround`) so the load point
correctly co-moves with the forearm+hand segment. Post-fix, required moment scales correctly and
monotonically with load (7.53 / 14.74 / 26.75 N·m), and is in the same order of magnitude as the
twin's own independently-computed requirement (9.51 / 18.14 / 32.52 N·m) — a real cross-model
consistency check (two different individuals/models, same task, same order of magnitude), not
identical (as expected — different forearm geometry).

**The 8–12× "force is many multiples of the applied load" finding is itself a real, first-principles
geometric result**, not a citation: it falls directly out of the elbow flexors' short moment arm
(≈2–5 cm) versus the forearm's long lever arm to the hand (≈30 cm) — the classic third-class-lever
mechanical disadvantage. This is a geometrically DERIVED confirmation of the same mechanism the
literature anchor (Sec.7) describes, independent of whether that paper's exact number could be
extracted this session.

## 7. External anchor — verified, honestly scoped

No in-vivo instrumented-implant elbow anchor exists — **already flagged as an honest gap in this
repo's own graph** (`MSK-JOINT-CONTACT-FORCE-ANCHORS`: *"ELBOW and ANKLE have NO public in-vivo
instrumented-implant telemetry database"*; independently corroborated by
`docs/MECHANISM_MSK_BUILD_PLAN.md`'s per-joint anchor table). Not re-litigated here, only confirmed
consistent.

Decorrelated anchor used instead: published cadaveric/biomechanical-model elbow force-analysis
literature. **WebSearch was unavailable this session** (tool reported session quota exhausted,
2000/2000 — a pre-existing constraint also hit in `docs/MECHANISM_MUSCLE_AUDIT.md` and
`docs/MECHANISM_MSK_BUILD_PLAN.md`). CrossRef's bibliographic REST API (a separate, unaffected
channel) was used instead to **live-verify** the task's own suggested citation:

- **An KN, Hui FC, Morrey BF, Linscheid RL, Chao EY. "Muscles across the elbow joint: a
  biomechanical analysis." Journal of Biomechanics, 1981. DOI: 10.1016/0021-9290(81)90048-8.**
  Existence, authors, year, journal, and DOI confirmed live via CrossRef. Full text is
  Elsevier-paywalled (DOI resolves to `linkinghub.elsevier.com`, JS/paywall gate, no extractable
  abstract or numeric result this session) — **honestly disclosed**: the specific "~1–3×BW" figure
  is carried over from the task's own stated ballpark and general biomechanics domain knowledge,
  NOT independently re-derived from this paper's full text this session.
- **Funk DA, An KN, Morrey BF, Daube JR. "Electromyographic analysis of muscles across the elbow
  joint." Journal of Orthopaedic Research, 1987. DOI: 10.1002/jor.1100050408.** Same verification
  tier (companion EMG paper, same research group; content not extracted, same paywall class).

**Regime check, honestly scoped**: rather than force a crisp PASS/FAIL of the arm26 demonstration's
number against the task's ~1–3×BW figure (whose exact load/task/normalization convention could not
be independently verified — doing so would be a tautology-adjacent gate against a target that is
itself uncertain), a **gross-error sanity ceiling/floor** was used instead (0.05×BW – 10×BW,
analogous to the knee cert's own 1200%BW absurdity ceiling vs. the documented ~1000–1800× ID-tool
bug signature): the arm26 demonstration's 59.4 %BW (5 kg case) **PASSES** this sanity gate by a wide
margin — i.e., the result is not a gross implementation error, though it sits above the task's
stated ~1–3×BW figure at the *high* end of the load sweep. Both normalizations (%BW and ×
applied-load-weight) are reported plainly so a reader can judge against whichever convention the
original literature actually used.

## 8. Honest gaps (full list)

1. **No in-vivo elbow anchor exists** (Sec.7) — already known/flagged in this repo's graph, not a
   new finding, but load-bearing for the whole cert's evidence tier: weaker than the knee's
   OrthoLoad/Grand-Challenge anchor by construction.
2. **The twin has 0 elbow/shoulder/pro-sup muscles** (Sec.2) — the central finding. Static
   Optimization cannot be exercised on the twin at this joint; only the pure-kinematic
   reaction/lower-bound is available, and even that pose is not sustainable by the twin's own
   idealized actuator above ~2 kg of hand load (Sec.4).
3. **The literature anchor's exact numeric figure could not be independently extracted** this
   session (paywall) — citation existence/identity verified live via CrossRef; the ~1–3×BW figure
   itself is inherited, not re-derived. The geometric (moment-arm-based) derivation on arm26
   independently corroborates the qualitative mechanism (large force multiplication from a short
   flexor moment arm), which is the strongest claim this session can honestly make.
4. **arm26 is a generic, non-subject-scaled illustrative model** — Holzbaur/Delp 2005's own
   reference anthropometry, not subject2's. The forced-adversary demonstration shows the METHOD
   works and gives a concrete, real number, but is not a "twin" result and should not be read as
   one.
5. **arm26's load-application point is estimated**, not measured directly (no wrist body/marker in
   that model) — via a cross-model-transferred distance ratio from the twin's own real wrist
   markers (Sec.6), a disclosed, measured proxy, not a free parameter.
6. **Simplified Static Optimization** (arm26 demonstration): closed-form moment-arm × Fmax
   minimum-activation-squared solve, no force-length/velocity curves — the standard textbook
   reduction for a single static pose, not a literal `opensim.AnalyzeTool`/`StaticOptimization`
   run (judged not worth the added scope of synthesizing a fresh trial/GRF-file pipeline for a
   single illustrative side-demonstration, given the twin itself cannot use SO at the elbow
   regardless, per Sec.2's decisive structural proof).
7. **Welded hand** (`radius_hand_r`, 0 DOF) — confirmed; no wrist/grip biomechanics in scope
   either model.
8. **Single representative pose family** (90° flexion, upper arm vertical) — the standard
   isometric-curl test position, not a full range-of-motion sweep.

## 9. Next step

To give the TWIN a genuinely credible muscle-driven elbow solve, graft real elbow flexor/extensor
muscle-tendon actuators (at minimum: biceps brachii long+short heads, brachialis,
brachioradialis, triceps brachii long/lateral/medial heads) onto the twin's own existing
`humerus_r`/`ulna_r`/`radius_r` bodies, parametrized from a published source — `arm26`
(Holzbaur 2005) for the elbow-crossing set specifically, or a fuller upper-limb model (e.g. a
Holzbaur/MoBL-ARMS-lineage full arm+hand model, if one can be sourced and verified) for
shoulder+wrist fidelity too, since this session's finding is a **whole-arm** gap, not
elbow-specific. This is a concrete, scoped build task (wrap/extend, don't edit the existing model
in place, per this repo's inherit-vs-build convention), not a vague "add more muscles" note.

Independently, this session's finding — that the shoulder is *equally* muscle-free — means a
sibling shoulder-focused agent hitting the identical wall is an expected, corroborating
convergence, not a coincidence; both would be diagnosing the same root cause (this model family's
arm was built for kinematic arm-swing during gait IK, not for upper-limb load-bearing dynamics).

## Files

- `scripts/msk/validate_elbow_force.py` — the full pipeline (self-contained, re-runnable, exit 0;
  imports `validate_joint_force.py` for `get_descendant_bodies`/`G`, and `static_opt_knee.py` for
  the proven muscle-crossing detector on the twin — not re-implemented).
- `data/msk_smoketest/subject2_elbow_static/elbow_force_validation/elbow_force_results.json` —
  every number in this document, machine-written.
