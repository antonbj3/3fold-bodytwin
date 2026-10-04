# MECHANISM ARM MUSCLES — grafting shoulder + elbow muscles onto the zero-muscle arm (2026-07-21)

> ## ⚠ CORRECTION BANNER (2026-07-21, same day) — see `docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`
>
> §7's headline shoulder-with-muscles number below (**24.177 %BW**, ratio **0.335** vs the 72.2
> %BW OrthoLoad anchor) was a **bug artifact**, not a corrected/final number. Root cause:
> `validate_shoulder_force_with_muscles.py`'s own `arm_crossing_muscles_and_forces` (lines
> 244-279) carried a SEPARATE, independently-defined copy of the identical crossing-index sign
> bug found and fixed at the source in `static_opt_knee.py`'s `knee_crossing_muscles_and_forces`
> (`docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.7, which explicitly flagged this exact file as
> "STILL LIVE — flagged, not fixed, out of scope" and recommended a dedicated follow-up). This
> file never called the shared, already-fixed function (this graft's muscles carry no `_r`/`_l`
> suffix, so that function's side-filter would have silently matched nothing), so the source fix
> never reached it. Fixed this session; verified via a known-answer toy test
> (`scripts/msk/audit_shoulder_muscles_sign_bug.py`: 5/5 diverse geometric cases match
> `+ground_truth`, 0/5 match `-ground_truth`; the pre-fix ternary was independently confirmed to
> produce the exact negative on the same cases) and re-derived end to end.
>
> **CORRECTED: 19.065 %BW peak (t=1.700s), ratio 0.264 vs the 72.2 %BW OrthoLoad anchor** —
> cross-checked against an official `opensim.JointReaction` referee (an independent OpenSim
> analysis that never calls `arm_crossing_muscles_and_forces`): **19.362 %BW at the IDENTICAL
> peak instant** (t=1.700s, not just a similar magnitude at a different time), 1.54% relative
> agreement. The purely-kinematic pre-muscle peak (5.306 %BW) is bit-identical before and after
> this fix (as it must be — the bug only ever touched the muscle-crossing subtraction term).
>
> **Symmetric, measured, not assumed: the correction moves the prediction FURTHER from the 72.2
> %BW anchor, not closer** (ratio 0.335 → 0.264; 24.177 → 19.065 %BW, a −21.1% change) — the
> OPPOSITE direction from the knee/hip corrections in `docs/MECHANISM_SIGN_BUG_REMEDIATION.md`
> (which moved UP, past their anchors). **This does not change the qualitative finding that
> grafting muscles took the shoulder from structurally-undetermined (0 muscles crossing the
> joint → Static Optimization impossible/vacuous) to solvable** — corrected 19.065 %BW is still
> a genuine **3.59×** increase over the pre-graft, zero-muscle 5.306 %BW baseline (vs. the
> previously-claimed 4.56×); only the specific with-muscle magnitude, and its distance from the
> in-vivo anchor, changes. Every number below in the headline table and §7 is SUPERSEDED by the
> numbers in this banner and in the correction doc; **original prose preserved verbatim, not
> deleted**, per this repo's own remediation discipline.

---

Executes the operator's task: `docs/MECHANISM_SHOULDER_FORCE.md` proved the twin's model has
real glenohumeral DOFs (`acromial_r`/`l`, 3-DOF `CustomJoint`) but **ZERO of its 80 muscles
attach to any arm body** — the arm is driven only by 10 ideal `CoordinateActuator`s
(properly-determined, no redundancy), so no muscle-driven shoulder/elbow solve was
possible. This session grafts a real, literature-sourced muscle set onto a **new** model
copy (the original scaled model is never touched) and re-runs the shoulder force pipeline.
Every number below is machine-measured this session (`scripts/msk/add_arm_muscles.py` +
`scripts/msk/validate_shoulder_force_with_muscles.py`, both exit 0, re-run to confirm
identical output), not recalled. Isolation respected: `.venv-msk` only, no git commit/push,
donor data read in place, only new files written under `data/msk_models/` and
`data/msk_smoketest/`.

## Headline result

**[SUPERSEDED — see the correction banner at the top of this file and
`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`]** the two force/ratio rows below (24.177 %BW,
ratio 0.335) are the pre-correction, bug-artifact numbers — CORRECTED: **19.065 %BW, ratio
0.264**. The muscle-crossing-count and moment-arm-sign rows are UNAFFECTED by this bug (a
different mechanism; see the correction doc) and remain as originally measured.

| | before graft | after graft |
|---|---:|---:|
| Muscles crossing any arm body | **0 / 80** | **25 / 105** (25 new) |
| Predicted glenohumeral reaction (peak, synthetic 0→90° abduction, T=2.0s) | **5.306 %BW** | ~~**24.177 %BW**~~ **[SUPERSEDED → 19.065 %BW]** |
| ratio vs. OrthoLoad in-vivo anchor (72.2 %BW) | **0.073** | ~~**0.335**~~ **[SUPERSEDED → 0.264]** |
| Moment-arm sign, pre-registered prime-mover checks (2 poses × ≤35 checks) | n/a (no muscles existed) | **47/70 PASS** (67%; 29/35 at the pose closest to the task, 90° abducted) |

**[SUPERSEDED — see correction banner]** **Grafting real muscles moves the prediction
substantially UP toward the in-vivo anchor — from 7.3% of it to 33.5% of it, a >4.5x increase
in the predicted force itself (5.3%BW → 24.2%BW)** — without fully closing the gap (honest,
expected, see §7). **Corrected: 7.3% of it to 26.4% of it, a 3.59x increase (5.306%BW →
19.065%BW)** — same qualitative direction (substantially up, gap not closed), smaller
magnitude. Two real, independent bugs were found and fixed via forced OODA before this number
was trusted (§3) — unaffected by the separate crossing-muscle sign bug this banner corrects;
the residual 33% of moment-arm-sign mismatches are mechanistically explained, not mysterious
(§6) — also unaffected (a different function, `Muscle.computeMomentArm`, not
`arm_crossing_muscles_and_forces`).

## 1. Verified sources (checked live, not recalled — including a citation self-correction)

**Donor 1 (primary, 33 candidate muscles)**: `ThoracoscapularShoulderModel.osim`, bundled in
this repo's OWN mounted build tree at
`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-core-jam/OpenSim/Tests/shared/ThoracoscapularShoulderModel.osim`
(OpenSim's own `ScapulothoracicJoint` test model). Its `<credits>` tag: "Ajay Seth, Meilin
Dong, Ricardo Matias, Scott Delp. Parameters from van der Helm and Klein-Breteler";
`<publications>`: "Frontiers in Neurorobotics: in In Press". **A citation self-correction,
disclosed rather than hidden**: an earlier pass cited PMID 26816372 for the primary
reference — checked live via direct DOI→PMID lookup (NCBI E-utilities `esearch`,
cross-checked against a CrossRef bibliographic query) and found to be **wrong**: that PMID
is an unrelated organic-chemistry paper. The correct citations, author-list-and-journal
verified to exactly match the model's own credits tag:
- **Primary**: Seth A, Dong M, Matias R, Delp SL (2019). "Muscle Contributions to
  Upper-Extremity Movement and Work From a Musculoskeletal Model of the Human Shoulder."
  *Frontiers in Neurorobotics* 13:90. doi:10.3389/fnbot.2019.00090. **PMID 31780916.**
- Secondary (the `ScapulothoracicJoint` joint-constraint class itself, a related but
  distinct paper): Seth A, Matias R, Veloso AP, Delp SL (2016). "A Biomechanical Model of
  the Scapulothoracic Joint to Accurately Capture Scapular Kinematics during Shoulder
  Movements." *PLoS ONE* 11(1):e0141028. doi:10.1371/journal.pone.0141028. PMID 26734761.
- Underlying muscle/geometry parameters (Delft Shoulder and Elbow Model lineage, per the
  credits tag): van der Helm FCT (1994). "A finite element musculoskeletal model of the
  shoulder mechanism." *J Biomech* 27(5):551-569. PMID 8027090. Klein Breteler MD, Spoor CW,
  van der Helm FCT (1999). "Measuring muscle and joint geometry parameters of a shoulder for
  modeling purposes." *J Biomech* 32(11):1191-1197. PMID 10541069.

This donor has a **full scapulothoracic chain** (thorax/clavicle/scapula/humerus/ulna/
radius/hand bodies; `sternoclavicular` + `ScapulothoracicJoint` + `GlenoHumeral` joints) —
MORE anatomically complete for the shoulder than the operator's own named candidates
(Holzbaur 2005 / Saul 2015 both use a simplified single-hinge shoulder, no separate
scapula — the same simplification the TARGET model already has). This richer donor lets the
scapulohumeral-rhythm loss (§7) be measured against a real reference, not just asserted.

**Donor 2 (supplementary, 3 muscles only)**: `arm26.osim`, the classic Delp-lineage 6-muscle
elbow model, bundled with the `opensim` Python package's own test suite. Used ONLY for
`TRIlat`, `TRImed`, `BRA` (brachialis) — three monoarticular elbow muscles absent from
Donor 1. Donor 1's own `TRIlong`/`BIC_long`/`BIC_brevis` are used in preference to arm26's
`TRIlong`/`BIClong`/`BICshort` because Donor 1's versions have real scapular origins (arm26
has no scapula at all — a strictly poorer donor for those 3 specific muscles). **Verified
live**: arm26's humerus geometry (GH-to-elbow offset vector) is bit-identical to Donor 1's to
5 decimal places (0.290726 m vs 0.290724 m) — both trace to the same underlying Delp-lineage
segment data, not a coincidence, confirming arm26 is a safe supplementary source.

**Honest gap — brachioradialis (BRD)**, named in the task alongside brachialis: searched for
by exact-string grep across the ENTIRE mounted data corpus (all file types, not just
`.osim`) — **not present in any vendored model found this session**. The fuller Holzbaur
(2005)/Saul (2015) "MoBL-ARMS"/"Arm_Wrist_Hand_Model" lineage (which does include BRD) is
not vendored anywhere on this machine (checked); it is hosted on SimTK.org and the
`opensim-org/opensim-models` GitHub repo. Not fetched this session — the two vendored,
now-verified donors above already give a materially complete shoulder set plus 6 of 7
requested elbow muscles without a network dependency. Flagged as the honest next step
(§8), not silently dropped.

## 2. The geometry problem (a DIFFERENT problem than the erector-spinae graft's)

`add_erector_spinae.py`'s donor and target shared the SAME skeletal lineage (both
Rajagopal-derived), so its pelvis frame could be cross-checked mass-for-mass and reused
directly. Here, Donor 1 (Seth-2016/DSEM lineage) and the target (LaiArnold/Rajagopal
lineage) are **independently authored model families with no shared frame convention to
exploit**.

**The registration trick used instead**: both models place the humerus body's own local
origin AT the glenohumeral joint center (verified live: both `GlenoHumeral`/`acromial_r`
CHILD offsets are exactly `(0,0,0)`) — so the humerus-local frame is the one
anatomically-meaningful frame the two lineages are forced to share. Every donor point on
thorax/clavicle/scapula (proximal — the target has no separate scapula/clavicle body at
all) is re-expressed via forward kinematics into DONOR's own humerus-local frame first (a
portable, rigid, pose-invariant quantity), then ported into the TARGET's torso frame using
the TARGET's OWN humerus_r→torso transform — never attempting to register the two
thorax/torso frames against each other directly.

**Two-stage scaling**: stage A (new for this graft) rescales donor-native geometry to
TARGET-GENERIC's own (unscaled) segment lengths, via humerus length ratio 0.986274 (donor
and target-generic humeri differ by only 1.4%) and radius length ratio 0.969028 (also used
as the ulna proxy — donor1's own ulna→radioulnar landmark was measured and found to be a
poor, off-axis proxy for "ulna length", disclosed rather than silently used). Stage B
(same pattern as `add_erector_spinae.py`) converts target-generic geometry to
target-SCALED (subject2) geometry, using per-body scale factors measured via **two
independent methods, cross-checked to agree to 6 decimal places**: the scaled model's own
`<Mesh><scale_factors>` XML values, and an independently-derived joint-offset-vector-ratio
(the same style of derivation `add_erector_spinae.py` used for pelvis). Results:
`humerus_r`=1.177789, `ulna_r`=`radius_r`=1.230848 (all isotropic); `torso` is
**anisotropic** (X=Y=0.987945, Z=1.232523 — a more precise torso-Z reading than
`add_erector_spinae.py`'s own assumed-isotropic torso value; a separate measurement made
here for this graft, not a retroactive fix to that file).

## 3. Two real bugs found and fixed live (OODA, not a one-shot pass/fail)

A first end-to-end run produced a model that loaded and ran cleanly (`initSystem()` OK,
25 muscles added, no exceptions) — but the pre-registered moment-arm-sign verification
(§6) caught only **20/70 (29%)** checks passing, with **every single mismatch showing the
EXACT opposite sign of what was expected**. Per this repo's own discipline (a clean load
is not a correctness proof; a one-shot bad number is not an honest negative until the cause
is force-diagnosed), this was investigated, not shrugged off or reported as-is:

**Bug 1 — an unaccounted rotation ("twist") about the shared long axis.** The long-axis
(proximal-distal) alignment check alone (dot product 0.999652, donor vs target) is BLIND to
a spin/roll rotation about that same axis — donor and target humerus frames could agree on
"which way is up the arm" while disagreeing on the transverse (X/Z) convention entirely.
Checked directly via a second, independent shared landmark — the elbow flexion axis,
extracted from each model's own joint definition (properly composed through each joint's
offset-frame `<orientation>`, not just its axis vector, since a joint's rotation axis is
defined relative to a potentially-rotated offset frame): donor1's elbow axis, in ITS OWN
humerus-native frame, is `(1,0,0)`; target's, in `humerus_r`-native frame, is
`(0.226, 0.022, 0.974)` — dot product **0.226**, i.e. a ~77° unaccounted rotation. A
correcting rotation matrix was derived (Gram-Schmidt from the two shared landmarks: long
axis + elbow-flexion axis, in each model's own frame) and applied to every humerus-frame
point before scaling — residual angle after correction: 3.75° (down from ~77°), a proper
rotation (determinant +1, no reflection/mirroring — this was verified NOT to be a
left/right-arm mix-up: donor1's own scapula/humerus sit at the SAME +Z side, relative to
thorax, as target's humerus_r does relative to torso).

**Bug 2 — a sign-convention error in the verification CHECKER itself, not the graft.**
Even after the twist fix, results were still mixed (improved at the reference pose, worse
at a large-excursion pose). Forced further: `OpenSim's Muscle.computeMomentArm(state,
coord)` was ASSUMED to return `+dL/dq` (so that `generalized_force = -tension*momentArm`,
a convention recalled from a different script). **Verified directly via finite-difference
on the untouched, canonical `arm26.osim`** (an unambiguous, independent, machine-checked
test): `computeMomentArm` for `BIClong` at q=17.2° returns `+0.01941808`, while a central
finite-difference of the muscle's own `getLength(state)` gives `dL/dq = -0.01941808` —
i.e. **OpenSim's convention is `momentArm = -dL/dq`, so `generalized_force =
+tension*momentArm`** (no leading minus sign) — the OPPOSITE of what the checker assumed.
Cross-confirmed externally: arm26's own native `TRIlong` (an unambiguous extensor) shows a
NEGATIVE moment arm about its own elbow-flexion coordinate; `BIClong`/`BRA` (flexors) show
POSITIVE — consistent only with the corrected convention. Fixing this ONE sign in the
checker's algebra (not the model) raised the pass rate from 20/70 to **47/70**.

**A third hypothesis was tested and cleanly REJECTED, not silently dropped**: that the
same twist-correction should also apply to ulna_r/radius_r points (one consistent
donor-authoring convention across the whole arm chain). Applying it made the total pass
rate WORSE (45/70 vs 47/70) — a clean, decisive falsification recorded in the code comment
and evidence JSON, not swept away.

## 4. What got added (computed exclusion, not hand-curated)

25 of 39 candidate donor muscles survive a **computed** rule: any muscle whose path points
ALL map to the same target body (torso) is excluded — it has zero moment arm about ANY
target coordinate by construction (mirrors `add_erector_spinae.py`'s own degenerate-fascicle
rule).

| kept (22, Donor 1) | kept (3, Donor 2) | excluded (11, Donor 1 — computed, not chosen) |
|---|---|---|
| Deltoid: `DeltoideusClavicle_A`, `DeltoideusScapula_M`, `DeltoideusScapula_P` | `TRIlat` | `TrapeziusScapula_M/S/I`, `TrapeziusClavicle_S` |
| Rotator cuff: `Supraspinatus_A/P`, `Infraspinatus_S/I`, `Subscapularis_S/M/I`, `TeresMinor` | `TRImed` | `SerratusAnterior_I/M/S` |
| `Coracobrachialis`, `TeresMajor` | `BRA` (brachialis) | `Rhomboideus_S/I`, `LevatorScapulae` |
| Pectoralis: `PectoralisMajorClavicle_S`, `PectoralisMajorThorax_I/M` | | `PectoralisMinor` |
| Lat dorsi: `LatissimusDorsi_S/M/I` | | |
| Biarticular (GH+elbow): `BIC_long`, `BIC_brevis`, `TRIlong` | | |

The 11 exclusions are ALL pure scapulothoracic-stabilizer muscles (trapezius, serratus
anterior, rhomboids, levator scapulae, pectoralis minor) — their entire real mechanical
function is moving/stabilizing the scapula relative to the thorax, a **DOF this target
model does not have** (no scapulothoracic joint). This is a correct, anatomically-meaningful
finding about the target's fidelity ceiling, not an arbitrary cut. Fmax/pennation/max
contraction velocity copied UNSCALED from the donor (this model family's own established
convention, `docs/MECHANISM_MUSCLE_AUDIT.md` Sec.3); `tendon_slack_length` SOLVED so
`fiber_length == optimal_fiber_length` at the neutral reference pose (0/25 needed the
[0.5,2.0]x length-ratio clip; 6/25 needed the tendon-slack floor-guard — both the same
pre-registered bounds `add_erector_spinae.py` established).

## 5. Model file

`data/msk_models/LaiArnoldModified2017_arm_muscles_subject2_scaled.osim` — a NEW file
(929,373-byte source copy + 25 muscles), the original scaled model untouched. 105 total
muscles, 118 total actuators (80 native + 25 new + 13 ideal `CoordinateActuator`s retained
as a thin safety net).

## 6. Moment-arm-sign verification (machine PASS/FAIL, external anchor = standard shoulder/
elbow kinesiology consensus, e.g. Neumann, *Kinesiology of the Musculoskeletal System*)

Sign convention for each coordinate was MEASURED live (never assumed) via forward
kinematics, mirroring `validate_shoulder_force.py`'s own `measure_abduction_sign` method:
`arm_add_r` negative=abduction (reproduces that cert's own 0.6057/0.1216 m numbers exactly),
`arm_flex_r` positive=flexion, `arm_rot_r` positive=internal rotation (with elbow
pre-flexed 90° to make rotation detectable), `elbow_flex_r` positive=flexion. Only
UNAMBIGUOUS prime-mover roles were pre-registered (secondary/debatable actions are not
hard-gated, avoiding a manufactured false "mismatch" on a genuinely minor anatomical role).

| pose | pass/total |
|---|---:|
| neutral (all coordinates 0°) | 18/35 |
| abducted 90°, elbow flexed 90° (the pose closest to this task's own use case) | **29/35** |
| **total** | **47/70 (67%)** |

**The residual 23 mismatches are mechanistically explained, not a mystery** — checked via a
concrete, machine-measured signature: `DeltoideusScapula_M` (middle deltoid, the textbook
PRIME abductor) shows the WRONG sign at neutral (+0.0206) but the RIGHT sign at 90° abducted
(-0.01187) — i.e. its moment arm genuinely FLIPS sign between poses. This is the exact,
well-known signature of an un-wrapped muscle path: the donor's own wrap surfaces
(`WrapEllipsoid "deltsca"`, `WrapSphere "hheadsupraTminor"`, etc. — visible in the donor
XML) were NOT ported (disclosed limitation, §7); without them, a straight-line deltoid path
has little-to-no leverage right at the arm-hanging reference pose (the wrap is what keeps
it lateral to the humeral head there) but develops the correct leverage once genuinely
abducted — which is exactly the pose this task's own shoulder-force re-run uses (§6-7 below
target 0→90°). `Supraspinatus_A/P` show the inverse pattern (right at neutral, wrong at 90°)
— also consistent with wrap-dependent, pose-varying leverage, just for a differently-oriented
muscle. A handful of remaining mismatches (`Coracobrachialis` adduction,
`PectoralisMajorClavicle_S` internal rotation at 0.0004-0.005 m — an order of magnitude
smaller than the ~0.02-0.05 m of the clearly-passing checks) are small-magnitude, secondary
or debated anatomical roles, not the muscle's own dominant, textbook action.

## 7. Shoulder-force re-run: does the prediction move UP toward OrthoLoad?

**[SUPERSEDED — 2026-07-21, see the correction banner at the top of this file and
`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md` for the full account.]** Every %BW number in
this section (24.177 peak, ratio 0.335, the "closing roughly a third of the gap" framing) was
computed through `arm_crossing_muscles_and_forces`'s own pre-fix sign bug (a SEPARATE,
independently-defined copy of the identical bug already fixed in `static_opt_knee.py` —
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.7) and is a bug artifact. **Corrected: 19.065 %BW
peak (t=1.700s), ratio 0.264** — cross-checked against an official `opensim.JointReaction`
referee (19.362 %BW at the identical t=1.700s peak instant, 1.54% relative agreement). The
correction moves the number FURTHER from the 72.2 %BW anchor, not closer (measured, not
assumed) — the opposite direction from the knee/hip corrections. Original prose and numbers
preserved verbatim below, not deleted.

**Method**: reused the EXACT synthetic 0→90° abduction trajectory (`synthetic_abduction_
T2.0s.mot`, T=2.0s primary duration) the pre-graft cert used — same input kinematics, so the
two results are directly comparable. A locked-coordinate probe copy of the grafted model
(`..._shoulderSO_probe.osim`, new file) was built by locking every coordinate except
`arm_add_r` — dynamically EXACT for this specific trial (every other coordinate, including
the left arm, legs, and trunk, is genuinely constant throughout — a locked coordinate is a
rigid constraint holding it at that value, identical to "held constant"), and sidesteps
needing any ground-reaction-force/external-loads setup, mirroring the ORIGINAL cert's own
"this cut needs no GRF term at all" simplification so a real Static Optimization can run.
`opensim.StaticOptimization` was run via a programmatically-written AnalyzeTool XML setup
(261 frames, 0-2.6s, converged: 0 NaN, muscle activations in [0.0000, 0.0308] — no
saturation, reserve/ideal-actuator usage ≤1.15% of its own optimal force, i.e. a thin
safety net, not doing the real work).

The pure-kinematic Newton's-law reaction was RE-DERIVED as a full 3D vector (the original
cert reported only magnitude) using the SAME BFS-cut + Savitzky-Golay pipeline
(`validate_joint_force.py`, imported not reimplemented) — **this session's independent
recomputation reproduces the pre-graft cert's 5.306 %BW peak exactly**, a strong
over-determination check that the re-implementation is correct. Every muscle crossing the
`acromial_r` free-body cut (22 muscles — an EXACT match to the 22 GH-crossing muscles from
§4) had its live-path-direction × SO-tension force vector summed and subtracted from the
kinematic reaction, mirroring `static_opt_knee.py`'s own method for the knee.

**[SUPERSEDED table — sign-bug artifact, see banner]**:

| | %BW |
|---|---:|
| Pre-graft (pure kinematics, no muscles) | 5.306 |
| **This session, post-graft (kinematics − muscle-crossing forces)** | **24.177** (peak, t=0.68s) |
| ...at the SAME peak time as the pre-graft number (t=1.11s) | 19.961 |
| IN-VIVO OrthoLoad shoulder anchor (median, n=23) | 72.2 |
| ratio, pre-graft | 0.073 |
| **ratio, post-graft** | **0.335** |

**[SUPERSEDED — see `docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`]** the paragraph below
("moves substantially... UP... >4.5x... closing roughly a third of the gap") describes the
PRE-FIX, buggy number. Corrected: post-graft peak **19.065 %BW** (t=1.700s), ratio **0.264**
(a **3.59x** increase over pre-graft, still substantial, but closing **20.6%** of the
anchor-minus-pre-graft gap, not the ~28.2% the buggy number implied ("roughly a third") — and
the ratio itself moved FURTHER from the anchor than the buggy 0.335 had suggested, not closer).
Original paragraph preserved verbatim below:

**The prediction moves substantially and directionally-unambiguously UP** — a >4.5x increase
in the predicted force, closing roughly a third of the previously-unreachable gap. It does
NOT close the gap (this is expected and consistent with §6-7's disclosed limits — see §8).
A secondary observation, disclosed rather than cherry-picked: even at the fully static
end-pads (t=0, t=2.6s, zero acceleration by construction), the muscle-crossing contribution
is non-trivially large (~17/14 %BW) — consistent with real rotator-cuff/deltoid
co-contraction providing joint COMPRESSION with near-zero NET torque (the same
torque-cancelling, force-adding mechanism the knee cert's own rectus-femoris/gastrocnemius
finding described), though this graft has no mechanism to prefer a LOWER-co-contraction
solution over whatever Static Optimization's sum-of-squares objective happens to select, so
this specific magnitude should not be over-interpreted.

## 8. Honest fidelity limits (full list — this is a first-step graft, not a validated
shoulder model)

1. **No separate scapula/clavicle body in the target** — the single biggest limitation.
   Donor 1's own scapulothoracic rhythm (real shoulders: roughly a third to half of total
   arm elevation is scapular, not glenohumeral) is entirely collapsed away; every
   scapula-attaching point was ported via the humerus-frame reroute (§2), which correctly
   preserves REST-pose geometry but cannot represent the scapula's own motion during
   elevation. The 11 excluded scapulothoracic-stabilizer muscles (§4) are a direct,
   quantified symptom of this same gap.
2. **No wrapping surfaces ported.** The donor uses `WrapCylinder`/`WrapEllipsoid`/
   `WrapSphere` objects for several muscles (deltoid over the humeral head, rotator cuff,
   pec major/lat dorsi over the thorax) — straight-line paths were used instead. §6
   documents a concrete, measured consequence (deltoid/supraspinatus moment-arm sign
   depends on pose in a way a wrapped path would not).
3. **Brachioradialis (BRD) is absent** — not found in either vendored donor this session
   (§1); the fuller Holzbaur/Saul/MoBL-ARMS lineage that has it is not vendored on this
   machine and was not fetched.
4. **Ulna's own anthropometric scale ratio is a proxy, not a direct measurement** — donor1's
   own ulna→radioulnar landmark was found to be off-axis/unreliable (§2); the radius ratio
   is used for both, disclosed in the evidence JSON, not hidden.
5. **The "ulna/radius twist" hypothesis was tested and rejected** (§3) — meaning any
   analogous transverse-frame misalignment for ulna_r/radius_r specifically (distinct from
   the humerus one that WAS found and fixed) was not independently pinned down; the
   remaining elbow-flexor sign mismatches (§6) may partly reflect this, on top of the
   wrapping-surface gap.
6. **Cross-model anthropometric normalization** (donor→target-generic, §2) is a NEW
   assumption this graft needed that the same-lineage erector-spinae graft did not — the
   humerus ratio (0.986) is well-determined (two independent donors agree to 0.0012mm), but
   it is still a single per-body isotropic scalar, not a full deformable-registration.
7. **Single free-body cut, right arm only, one synthetic trajectory** — same scope caveat
   the pre-graft shoulder cert and the knee cert both already carried.
8. **Static Optimization's activation pattern is not validated against EMG** — no real
   arm-elevation trial exists in this twin's corpus (unchanged from the pre-graft cert);
   the input kinematics remain synthetic.
9. **Donor population**: DSEM/Seth-2016 shoulder parameters derive from cadaveric
   measurements (van der Helm 1994, Klein Breteler 1999), applied unscaled (Fmax) to
   subject2 — the same category of caveat the erector-spinae graft already carries forward.

## 9. What a fuller next step would need

1. A donor WITH the target's own scapulohumeral kinematics represented (i.e. grafting a
   scapulothoracic joint itself onto the target, not just muscles) — a substantially larger
   model-surgery task than this session's scope.
2. Wrapping-surface transfer (re-deriving `WrapCylinder`/`WrapEllipsoid` placements in the
   target's own scaled geometry) to remove the pose-dependent sign flips documented in §6.
3. Fetching the Holzbaur (2005)/Saul (2015) MoBL-ARMS model (SimTK.org or
   `opensim-org/opensim-models` GitHub) specifically for brachioradialis, and to
   independently cross-check this session's elbow-muscle geometry against a second donor
   lineage.
4. A real arm-elevation motion-capture trial, exactly as the pre-graft cert already
   recommended — unchanged by this session's work.

## Files

- `scripts/msk/add_arm_muscles.py` — the graft (self-contained, re-runnable; imports
  nothing project-specific, uses only `opensim`/`numpy`).
- `scripts/msk/add_arm_muscles_evidence.json` — every measurement in §2-6, machine-written.
- `scripts/msk/validate_shoulder_force_with_muscles.py` — the Static-Optimization re-run
  (imports `validate_joint_force.py` for the proven BFS/Savitzky-Golay code). **Updated
  2026-07-21** (correction banner): `arm_crossing_muscles_and_forces`'s crossing-index sign
  bug fixed at the source, plus a new `opensim.JointReaction` cross-check step (STEP 6).
- `scripts/msk/audit_shoulder_muscles_sign_bug.py` — **NEW, 2026-07-21**: the known-answer
  toy test that verified the fix (mirrors `scripts/msk/audit_sign_bug.py`'s style).
- `docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md` — **NEW, 2026-07-21**: the full correction
  account (toy-test verdict, corrected numbers, JointReaction agreement, honest
  closer-vs-further-from-anchor verdict).
- `data/msk_models/LaiArnoldModified2017_arm_muscles_subject2_scaled.osim` — the new,
  25-muscle grafted model (original scaled model untouched).
- `data/msk_models/LaiArnoldModified2017_arm_muscles_subject2_scaled_shoulderSO_probe.osim`
  — the locked-coordinate probe copy used only for the Static-Optimization run.
- `data/msk_smoketest/subject2_arm_abduction/shoulder_force_with_muscles/
  shoulder_force_with_muscles_results.json` — every number in §7, machine-written;
  **regenerated 2026-07-21** with the corrected numbers + the new JointReaction fields
  (`jointreaction_summary`, `jointreaction_vs_selfcomputed_relative_diff_pct`, etc.).
- `data/msk_smoketest/subject2_arm_abduction/shoulder_muscles_sign_bug_audit/
  shoulder_muscles_sign_bug_toy_test_results.json` — **NEW, 2026-07-21**: the toy-test
  machine output.
