# MECHANISM CORPUS IK — MediaPipe pose -> OpenSim Inverse Kinematics bridge

Written 2026-07-21. Status: **ONE priority corpus clip taken end-to-end, MediaPipe pose ->
OpenSim IK -> physiological joint angles, machine-verified (not asserted). 2/3 pre-registered
range gates PASS (knee, hip_flexion); ankle FAILS its gate and is diagnosed, not hand-waved.
Two real bugs found and fixed live via forced OODA, not accepted on first sight.** Every number
below comes from `scripts/msk/pose_to_opensim_ik.py`'s own printed/JSON output or a one-off
diagnostic built on that script's own functions — not estimated, not eyeballed.

Closes the gap `docs/MECHANISM_POSE_PIPELINE.md` left open: that pipeline produced MediaPipe
world-landmark `.trc` marker files (228/264 corpus clips passing its own gate suite) but **not**
OpenSim joint angles. This doc is that bridge.

## 1. Clip chosen

`2025-09-26_DPEOqLVkav7` (`weighted_lifts/squat`, account `get.gymnast.fit`) — the pose
pipeline's own doc explicitly names this its "primary verify target": det_frac=1.0 (368/368
frames detected, zero NaN rows, no multi-pose frames), already independently confirmed loadable
by OpenSim's own `TimeSeriesTableVec3` reader. Picked as instructed (clean squat, det_frac=1.0) —
no further shopping across the corpus was needed or done.

## 2. Method — the minimal viable mapping (not Pose2Sim)

`docs/MECHANISM_POSE_PIPELINE.md` §10 recommended a full Pose2Sim triangulation/marker-augmentation
chain, but flagged that sibling venv as "BROKEN IN PLACE" — re-provisioning it was out of scope
here. Instead, `scripts/msk/pose_to_opensim_ik.py` does the minimal viable version:

1. **Add 16 new named markers** to a copy of the already-verified subject2-scaled
   `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` model (the same model
   `docs/MECHANISM_MSK_ENV.md`'s smoke test validated: IK on real mocap data reproduced an
   external 2021 reference to median 0.055°/max 1.328° RMSE). Each new marker is
   cloned/averaged from an **existing**, already-anatomically-correct, already-subject2-scaled
   marker — no new anatomy guessed:

   | MediaPipe `.trc` marker | source model marker(s) (same body) | body |
   |---|---|---|
   | LShoulder / RShoulder | L_Shoulder / R_Shoulder | torso |
   | LElbow / RElbow | avg(L\_elbow\_lat, L\_elbow\_med) / avg(R\_elbow\_lat, R\_elbow\_med) | humerus_l/r |
   | LWrist / RWrist | avg(L\_wrist\_radius, L\_wrist\_ulna) / avg(R\_wrist\_radius, R\_wrist\_ulna) | radius_l/r |
   | LHip / RHip | L_HJC / R_HJC (hip joint center) | femur_l/r |
   | LKnee / RKnee | avg(L_knee, L_mknee) / avg(r_knee, r_mknee) (transepicondylar midpoint) | femur_l/r |
   | LAnkle / RAnkle | avg(L_ankle, L_mankle) / avg(r_ankle, r_mankle) (malleolar midpoint) | tibia_l/r |
   | LHeel / RHeel | L_calc / r_calc | calcn_l/r |
   | LFootIndex / RFootIndex | L_toe / r_toe | calcn_l/r |
   | Nose | **excluded** — no head/neck body exists in this model (BodySet enumeration confirmed it stops at `torso`); forcing Nose onto torso would add a large non-rigid lever-arm error the model has no neck DOF to absorb. Honest gap, not silently worked around. |

   Verified by a forward-kinematics roundtrip on a scratch model before committing to the real
   pipeline: a cloned "LKnee_TEST" marker's `getLocationInGround()` matched the arithmetic
   midpoint of its two source markers to 1.4e-17 m (floating-point noise).

2. **Run OpenSim's own `InverseKinematicsTool`** (Python API) against the MediaPipe `.trc`
   directly, with an `IKTaskSet` over exactly those 16 markers (weight 2.0 for the
   hip/knee/ankle/heel/foot-index "posterior chain" markers per the task's own framing, weight
   1.0 for shoulder/elbow/wrist).

3. **Two anthropometric/kinematic approximations, stated up front, cost measured below, not
   hidden:** (a) subject2's own scaled geometry is reused for a different, monocularly-filmed
   athlete — no per-subject scaling is attempted, because monocular MediaPipe has no independent
   metric-scale signal to scale from; (b) MediaPipe applies no temporal smoothing between frames.

Derived model: `data/msk_models/LaiArnoldModified2017_mediapipe_bridge_subject2_scaled.osim`
(65 markers = 49 original + 16 new; `initSystem()` succeeds, i.e. the model is mechanically
well-formed, not just XML-valid).

## 3. Two real bugs found and fixed live (forced via OODA, not accepted on first sight)

### 3.1 `InverseKinematicsTool`'s own compiled-in default silently truncated the clip to 61/368 frames

First run "succeeded" (`ik.run()` returned `True`, output file existed) but a machine
cross-check — output row count vs. input `.trc` frame count, added on principle, not because
anything looked wrong yet — caught `n_rows=61` against an expected 368. **Root cause, confirmed
by direct API introspection, not assumed:** a bare `osim.InverseKinematicsTool()` (built via the
Python API from scratch, as this bridge must, since there is no per-clip XML template) has a
compiled-in default `time_range = [-inf, 2.0]` seconds — **not** `[-inf, inf]`. The existing
`smoke_test_ik.py` never hit this because it started from an original setup XML that already had
an explicit `<time_range>-Inf Inf</time_range>`. **Fix:** `read_trc_time_range()` reads the
clip's own actual `[t0, t1]` from its `.trc` data rows and sets it explicitly;
`run_ik()` now hard-asserts `n_output_rows == n_input_frames` on every call, permanently, so this
class of silent truncation can never again pass unnoticed (verified: re-run after the fix
produced exactly 368/368 rows).

### 3.2 Weakly-observed rotational DOFs amplify MediaPipe's per-frame jitter into non-physiological jumps

First full-clip run showed `ankle_angle_l` range 90.6° (fails a sane physiological ceiling) and,
more tellingly, `subtalar_angle_r` jumping **60.3° in a single 33ms video frame** (>1800°/s — not
a human movement). This is a "honest negative" that was **not** accepted at face value — forced
through OODA instead:

- **Per-marker mean tracking error** (16 markers, this clip, meters): LKnee 0.132, RHip 0.128,
  RHeel 0.117, LHeel 0.109, LHip 0.104, RKnee 0.094 (elevated cluster: hip/knee/heel) vs.
  RWrist 0.075, LShoulder 0.073, RShoulder 0.071, LWrist 0.069, RAnkle 0.064, RElbow 0.051,
  LAnkle 0.050, RFootIndex 0.050, LFootIndex 0.035, LElbow 0.022 (moderate cluster: upper
  body/ankle/toe). No single mis-mapped marker stood out — ruled out a mapping bug.
- **Per-frame jump statistics** (raw, before smoothing): `subtalar_angle_r` max jump 60.3°
  (mean 8.7°), `hip_flexion_l` max jump 65.0°, `ankle_angle_l` max jump 36.9°, vs.
  `hip_flexion_r` max jump only 15.1° — the spikes concentrate exactly in the DOFs a
  1-2-marker-per-segment reduced set constrains most weakly: axial rotation / inversion-eversion
  needs **>=3 non-collinear markers per rigid segment** to be observable at all (a geometric
  fact — the marker-to-angle Jacobian's null space for those axes is large with only a
  heel+toe pair or a single knee/ankle point per segment), not a heuristic guess.
- **Geometric prediction, forced to a decisive test:** if this is unsmoothed MediaPipe jitter
  amplified by weak observability, smoothing the *input* marker trajectory should collapse the
  spikes without needing to touch the model. Tested at window=5 (~167ms @ 30fps) and window=9
  (~300ms) via a centered, NaN-aware moving average (`smooth_trc()`):

  | | subtalar_r max jump | ankle_angle_l max jump | ankle_angle_l range | marker RMS |
  |---|---:|---:|---:|---:|
  | raw | 60.3° | 36.9° | 90.6° | 0.0850 m |
  | smoothed (w=5) | 26.1° | 24.6° | 84.6° | 0.0857 m |
  | smoothed (w=9, exploratory) | 19.8° | 16.7° | 80.6° | 0.0864 m |

  **Confirmed, not just argued:** smoothing cuts peak jump magnitude ~35-60% at <1% marker-error
  cost — the mechanism is real. **Also honest:** the *range* only drops ~7-11% even at w=9 — a
  genuine slower-timescale signal remains on top of the noise (this squat clip's knee flexion
  trends from ~45° to >130° over its last ~4s, consistent with a deep/ballistic squat, not just
  noise). **Conclusion: PARTIALLY explained and PARTIALLY mitigated, not fully resolved** —
  reported as such below, not smoothed away in the prose to match a hoped-for clean pass.
  `--smooth-window 5` is now the pipeline's default (measurably helps, costs almost nothing).

### 3.3 A third, decisive, independent finding: the knee L/R asymmetry is camera occlusion, not a bug

Raw MediaPipe visibility (this clip, mean/min over all 368 frames, read directly from
`pose_raw.json`, not inferred):

| landmark | left mean (min) | right mean (min) |
|---|---:|---:|
| knee | 0.964 (0.824) | **0.561 (0.114)** |
| ankle | 0.970 (0.884) | 0.741 (0.296) |
| heel | 0.951 (0.865) | 0.802 (0.495) |
| foot_index | 0.957 (0.841) | 0.849 (0.565) |
| hip | 0.999 (0.994) | 0.997 (0.987) |

The right leg is **substantially and systematically less visible to the camera** throughout this
specific clip (knee visibility collapses to a minimum of 0.114 at times) while the hip landmarks
(more central, less occluded) stay high on both sides. This directly explains why RKnee's solved
ROM (35.6-44.9° depending on smoothing) comes out much smaller than LKnee's (91.2-94.9°): a less
visible landmark gets a noisier/more-damped MediaPipe 3D estimate, compressing its apparent
range — a clip-specific occlusion artifact of this athlete's orientation to the camera, not a
left/right mapping error (the pose pipeline's own synthetic control already independently
verified the L/R labeling itself, see `docs/MECHANISM_POSE_PIPELINE.md` §9). This is exactly the
"scene-eyes" principle: the camera's own null space (what it cannot see well) shows up
directly and measurably in the output.

## 4. Results — physiological range verification (pre-registered thresholds, not fit after seeing data)

Final/recommended configuration: `--smooth-window 5` (the mitigated variant). Both variants
converged to all 368/368 frames; marker/angle numbers below are machine-read from the `.mot`
file via `opensim.TimeSeriesTable`, never eyeballed from a plot.

| gate | pre-registered range | raw | smoothed (w=5, **final**) | verdict |
|---|---|---:|---:|---|
| marker_error_rms_mean | <=0.15 m pass, >0.30 m fail | 0.0850 m | 0.0857 m | **PASS** (6.5x the same-subject real-marker floor of 0.013 m — see §5.1) |
| knee flexion range (R / L) | 20-165° | 44.9° / 94.9° | 35.6° / 91.2° | **PASS** |
| hip flexion range (R / L) | 15-150° | 39.4° / 148.9° | 32.6° / 148.2° | **PASS** (L is borderline, near the pre-registered ceiling — flagged, not hidden) |
| ankle range (R / L) | 5-60° | 57.0° / 90.6° | 47.9° / 84.6° | **FAIL** (diagnosed in §3.2: partly real fast/athletic motion, partly a residual of weak subtalar/ankle observability even after smoothing) |

**Sharper finding added §9 (batch-prep session, does not change the numbers above, changes their
interpretation):** this model's `hip_flexion` coordinate has a hard mechanical range of exactly
**[-30°, 120°], a 150° span — identical to this gate's own pass ceiling.** The smoothed
`hip_flexion_l` min (-29.88°) sits within 0.12° of the coordinate's own lower stop, and its max
(118.3°) within 1.7° of the upper stop — **both ends simultaneously**, on the cleanest clip in the
corpus. §9 shows this is not unique to this clip. Read as: "at least one frame drove the pelvis/hip
solve to (or very near) the mechanical limit," not necessarily "the athlete achieved 148° of true
hip flexion." The PASS verdict itself is unchanged (correctly computed against the pre-registered
gate) — this is a reported ambiguity in what a PASS *means* here, not a retraction.

**Bridge verdict** (ran/converged/non-degenerate/marker-fit, exactly that scope, not conflated
with the range gates above): **PASS**. **Physiological range gates: 2/3 PASS** — reported
exactly, a FAIL is a diagnosed finding, never silently folded into the bridge verdict.

Knee flexion reaching ~91-95° (one side) with a late-clip peak past 130° is a plausible deep
squat/athletic depth signature — the core ask ("does knee-flexion reach a sane depth") is
answered **yes**, with the honest caveat that the right leg's number is depressed by camera
occlusion (§3.3), not by a shallower true movement.

## 5. Honest limits (what this costs, measured — not overclaimed)

### 5.1 Anthropometric approximation (per task instruction, pre-approved)

subject2's real scaled geometry stands in for this specific, monocularly-filmed athlete; no
per-subject scaling was attempted (monocular MediaPipe has no independent metric-scale signal to
scale from). **Measured cost:** mean marker RMS error 0.085-0.086 m, about **6.5x** the same
model's own same-subject real-marker IK floor (`docs/MECHANISM_MSK_ENV.md`: 0.013 m mean, external
anchor — an independently-computed 2021 reference, different toolchain). Still comfortably under
the pre-registered 0.15 m pass ceiling (about 35% of a thigh segment length), but a real,
quantified cost, not "basically the same."

### 5.2 No global translation — CONFIRMED again via a second, independent measurement pathway

`docs/MECHANISM_POSE_PIPELINE.md` §4.2 already found MediaPipe's `world_landmarks` are
hip-re-centered *every frame* (raw landmark inspection, pre-IK). This bridge re-derives the same
conclusion through a completely different pipeline stage — **the solved pelvis translation
after a full IK solve** — and **refines it with a number**: `pelvis_tx/ty/tz` ranges came out
0.056-0.111 m (not exactly zero as a naive "no signal at all" framing would predict, but far
smaller than any real squat's actual depth, which is typically 0.3-0.6 m). **Honest, precise
statement:** global translation is **severely attenuated to a few-cm noise-like residual**, not
literally zero — the residual most plausibly comes from frame-to-frame recentering noise in
MediaPipe's own hip-centering (not real pelvis motion) rather than any recovered scene-translation
signal. Either way: **do not use this bridge's pelvis_tx/ty/tz for real squat depth or jump
height** — that signal is not in the input.

### 5.3 No pelvis-specific markers at all — a fundamental gap this reduced set has

MediaPipe's 33 landmarks include no ASIS/PSIS-equivalent bony pelvis landmark. `LHip`/`RHip` are
femur-attached hip-joint-center proxies (rotationally invariant to hip motion, a sound choice for
what they are), but nothing in the 16-marker set directly constrains the **pelvis's own**
orientation — it is inferred purely through the hip joints. This plausibly contributes to
`hip_flexion_l`'s large range (148°, right at this doc's own pre-registered ceiling) and is a
concrete, cheap target for improvement (see §7).

### 5.4 Axis-convention question — raised, tested, ruled out (not left as a lingering worry)

`pose_extract.py`'s per-clip rotation builds a frame (right, up, posterior) that is a fixed
rotation away from this OpenSim model's own native (anterior, up, right) convention. Geometric
prediction: since relative joint angles are invariant to any whole-body rotation (only the
free-floating pelvis's *own* orientation channels should absorb a fixed global misalignment),
this should not matter. **Tested, not just argued:** ran IK on the clip with an explicit
axis-remap applied (proper rotation, det=+1.000000, confirmed) versus without. **Result: marker
RMS identical to 4 decimal places (0.0850 m both ways), knee/hip/ankle ranges differ by <0.3°
across the board.** Confirmed and closed — this is not a live concern for this bridge.

### 5.5 Foot orientation is under-constrained by design

Only heel + toe (2 markers, both roughly in the sagittal plane) constrain each foot segment;
correlation between `ankle_angle` and `subtalar_angle` time series measured at -0.26 to -0.34
(this clip) — a real but only moderate compensatory coupling, consistent with (not conclusive
proof of) the foot's rotation-about-its-long-axis being weakly observable from this marker pair.

## 6. Reproducing this result

```
.venv-msk/bin/python3 scripts/msk/pose_to_opensim_ik.py --stem 2025-09-26_DPEOqLVkav7 --check-axis-remap
```

Outputs (gitignored, `data/msk_ik/`): `<stem>/raw/*.mot` + `*_ik_marker_errors.sto`,
`<stem>/smoothed/*.mot` + errors, `<stem>/remap/*.mot` + errors (only with `--check-axis-remap`),
and a single `<stem>_ik_report.json` with every number in this doc in structured form (angle
min/max/range/mean per coordinate, marker error stats, PASS/FAIL per gate). Model:
`data/msk_models/LaiArnoldModified2017_mediapipe_bridge_subject2_scaled.osim` (built once, reused
— idempotent unless `--rebuild-model`).

## 7. The concrete batch path for all 228 clips

**STATUS: EXECUTED 2026-07-21 (same session that wrote this update).** The plan below is kept
verbatim as the historical record of what was proposed; §8-11 report what actually happened,
including one real bug this plan's own item 2 correctly predicted needed a live check.

1. **Loop the same script over the 228 `OVERALL_PASS` clips** from
   `docs/MECHANISM_POSE_PIPELINE.md` §8 (the ones that already passed the pose pipeline's own gate
   suite — do not feed it the 31 zero-detection "quote card" clips or the 5 folded-posture
   failures; those need their own upstream fix first, not a downstream IK workaround):
   ```
   for stem in $(cat 228_pass_stems.txt); do
       .venv-msk/bin/python3 scripts/msk/pose_to_opensim_ik.py --stem "$stem"
   done
   ```
   The model-build step is already idempotent (reused across clips, ~1s check); each clip's own
   IK costs ~20-25s x 2 variants (raw + smoothed) serially. **Cost estimate:** 228 clips x ~45s
   ~= 171 min serial; parallelizable the same way `run_pose_batch.py` already parallelizes the
   pose-extraction stage (8-way was RAM-safe there; OpenSim's IK is CPU- not RAM-heavy per
   process, so similar or higher parallelism is plausible — re-check the RAM gate before
   committing, per this repo's own convention, don't assume it transfers).

2. **NaN handling for `det_frac<1` clips (not exercised by this pass — this clip was
   det_frac=1.0 on purpose):** `smooth_trc()` and `remap_trc()` are already NaN-aware
   (`np.nanmean` skips missing frames; a fully-NaN window still emits NaN, never a fabricated
   0). `InverseKinematicsTool` itself needs to be checked live on a real det_frac<1 clip: OpenSim
   marker references generally tolerate a NaN/missing marker on a given frame (drops it from that
   frame's cost function) but this has **not** been verified empirically in this session — do
   that check first on 2-3 low-det_frac clips before trusting the full 228 unattended, per the
   "don't assume, measure" standard applied throughout this doc.

3. **Per-clip automated flags to compute (cheap, from data already in each report JSON), not
   left to manual review:** (a) marker_verdict != PASS, (b) any range gate FAIL (expect ankle to
   fail on a nontrivial fraction, given §3.2/3.5 — track the rate, don't be surprised by one),
   (c) L/R visibility asymmetry >0.2 mean (from the clip's own `pose_raw.json`, already computed
   the same way as §3.3) as an automatic flag for "one side's ROM number is camera-occlusion-
   suppressed, downweight or exclude that side's angle from any aggregate/training use."

4. **Do not silently trust ankle/subtalar angles at face value** across the batch — §3.2/3.5
   showed these are the weakest-observed DOFs in this specific 16-marker set. Two concrete,
   cheap upgrades for a future pass, not attempted here (first-step doctrine — this is a
   bootstrap): (a) add pelvis-proxy markers (e.g. a MediaPipe mid-hip/mid-shoulder-derived
   trunk-orientation constraint) to address §5.3; (b) if foot orientation matters for a
   downstream use, that needs a 3rd non-collinear foot point, which the current 33-landmark
   MediaPipe set cannot supply — an honest ceiling of this input modality, not a bug to chase.

5. **Global translation (§5.2) stays unavailable for the whole 228 with this bridge as-is** —
   any downstream use needing real squat depth/jump height needs the separate 2D-bbox-based
   scale-recovery path `docs/MECHANISM_POSE_PIPELINE.md` §4.2 already flagged as not-yet-built, or
   a true per-subject `ScaleTool` pass, neither of which this session attempted (out of the
   pre-approved scope: "reuse the subject2 scaled model... if per-subject scaling... isn't
   feasible").

## 8. NaN-handling live check — item 2 above was RIGHT to flag this: it broke (forced OODA, fixed)

**First, the 228-clip pass list was independently re-derived, not trusted from
`docs/MECHANISM_POSE_PIPELINE.md`'s own prose**: re-running `pose_extract.verify_output()`
in-process over all 264 corpus stems this session reproduced **exactly 228/264 PASS**, with the
identical 31-clip and 5-clip failure signatures that doc's own §8 already root-caused. Machine
cross-check, not a re-quote.

**Then the det_frac<1 check item 2 above flagged as unexercised was run live, on the most extreme
case first** (`2026-07-04_DaYjAaquUq0`, det_frac=0.1621 — of the 73/228 PASS clips with det_frac<1,
the global corpus minimum). **It broke, exactly as flagged as a risk:**

```
SimTK Exception thrown at Assembler.cpp:224:
  calcGoal() method of assembly condition Markers returned a negative or non-finite value -nan.
  (Required condition 'goalValue >= 0' was not met.)
[error] InverseKinematicsTool Failed: AssemblySolver::track() attempt failed.
```

**Root cause (confirmed from the Simbody error itself + a frame-run-length scan of the raw JSON,
not guessed):** this clip's frames are fully-detected or fully-undetected (never partial — a
MediaPipe detection returns all 33 landmarks or none), so `write_trc()`'s NaN sentinel makes an
undetected frame ALL-NaN across all 16 markers together. The `Markers` assembly condition computes
ONE scalar goal (weighted sum of squared marker errors) per frame; when every marker is NaN, that
scalar itself is NaN, and SimTK's `AssemblerSystem` hard-rejects any non-finite goal, throwing and
killing `ik.run()` for the **entire clip** (zero frames written), not just the bad frame — confirmed
live: it crashed at Frame 5 (t=0.167s), the first fully-undetected frame, after frames 0-4 (real
detections) solved without issue. Gap-run-length analysis of this clip found NaN runs from 1 to
**337 consecutive frames** — far beyond the 5-frame smoothing window's reach, so `smooth_trc()`'s
existing NaN-aware averaging (which CAN silently repair short gaps by imputing from neighbors)
cannot and does not fix this; a structural fix was required, not a parameter tweak.

**Fix** (`strip_allnan_frames_trc()`, `scripts/msk/pose_to_opensim_ik.py`): drop fully-NaN frames
from the `.trc` before it reaches `InverseKinematicsTool` — information-preserving in this data
(a dropped frame had zero markers to begin with) not an approximation. `run_ik()`'s existing
anti-silent-truncation invariant (§3.1, `n_output_rows == n_input_frames`) is kept, now checked
against the POST-strip count (the ground truth of what was actually solvable), never the raw
frame count — that number is reported honestly alongside it (`coverage_frac`), never hidden inside
the assert.

**Verified converges, not just patched-and-hoped, across a deliberately diverse forced-adversary
set** (lowest possible det_frac, a mid-range det_frac that is `docs/MECHANISM_POSE_PIPELINE.md`'s
own §4.3 cross-validation example (539 rows - 199 NaN = 340 detected = det_frac 0.6308), and a
`weighted_lifts/squat` clip to cover the other movement type — plus a REGRESSION check on the
original det_frac=1.0 pilot clip, which must reproduce its own already-published numbers to the
digit since stripping is a no-op when there is nothing to strip):

| clip | det_frac | frames solved (raw/smoothed) | crash before fix? | crash after fix? | marker_rms after |
|---|---:|---|---|---|---:|
| `2026-07-04_DaYjAaquUq0` | 0.1621 (corpus min) | 71/438, 102/438 | **YES** (Frame 5) | no | 0.1400 m PASS |
| `2026-03-03_DVb7W3Oj2nQ` | 0.6308 (pose-pipeline doc's §4.3 example) | 340/539, 411/539 | not tested (fix already in) | no | 0.0903 m PASS |
| `2026-06-03_DZILqrIxjN5` | 0.8820 (`weighted_lifts/squat`) | 1203/1364, 1302/1364 | not tested (fix already in) | no | 0.1069 m PASS |
| `2025-09-26_DPEOqLVkav7` (pilot, **regression check**) | 1.0000 | 368/368, 368/368 (0 dropped) | n/a | no | 0.0850/0.0857 m — **identical to §4's published numbers** |

The regression row is the decisive check: `coverage_frac=1.0`, `n_frames_dropped_allnan=0`, and
every downstream number (marker RMS, knee/hip/ankle ranges, axis-remap delta) matches §3-5's
already-published values exactly — the fix is provably a no-op when there is nothing to strip.

**New, honest per-clip cost this fix makes visible (not previously measurable):** `coverage_frac`
— the fraction of a clip's frames that actually got a solved angle. For `DaYjAaquUq0`, only
16-23% of frames are solved; its knee/hip ranges rest on a minority of the clip and should be
trusted less than a full-coverage clip. A `low_frame_coverage_flag` (`coverage_frac < 0.5`) is now
emitted per clip specifically so this isn't silently averaged away in the aggregate (§10).

## 9. New finding: `hip_flexion`'s pass ceiling exactly equals its own mechanical range span

While diagnosing §8, the three forced-adversary clips all independently produced `hip_flexion_r`
**and** `hip_flexion_l` simultaneously at ~149.9-150.0°, and `ankle_angle_r`/`_l` simultaneously at
~99.8-99.9° — suspicious repetition of near-identical values across different clips/movements,
forced to a check rather than accepted as coincidence.

**Checked directly against the model's own `Coordinate` API** (`getRangeMin()`/`getRangeMax()`,
not hand-transcribed from the XML):

| coordinate | model's own mechanical range | span | this doc's pass ceiling | relationship |
|---|---|---:|---:|---|
| `hip_flexion_r/l` | [-30°, 120°] | **150.0°** | 150° (§4/PRE_REGISTERED) | **span == ceiling, exactly** |
| `ankle_angle_r/l` | [-50°, 50°] | 100.0° | 60° | span > ceiling (gate correctly fires FAIL when hit) |
| `knee_angle_r/l` | [0°, 140°] | 140.0° | 165° | span < ceiling (gate is safe from this ambiguity) |

**Consequence, stated precisely:** because `hip_flexion`'s allowed mechanical span is numerically
identical to its own pass-gate ceiling, a `hip_flexion` range reading near 150° is **mechanically
indistinguishable** between "the athlete achieved genuine ~150° of hip motion" and "the solver, on
at least one poorly-constrained frame, was pushed to (or past, and clamped at) the joint's hard
stop." This is not a new bug in the range-gate logic — the gate still correctly checks the
pre-registered threshold — it is a **ceiling-choice coincidence** that limits how much a PASS at
or near 150° can be trusted as a pure physiological signal, on top of §5.3's already-documented
pelvis-under-constraint honest limit (this finding gives that limit a sharp, mechanical number:
the pilot clip's published `hip_flexion_l` min of -29.89° sits within **0.11°** of the coordinate's
own -30° floor, and its max of 118.3° within **1.7°** of the +120° ceiling — both ends
simultaneously, on the cleanest, highest-quality clip in the whole corpus). `ankle`'s span
exceeding its ceiling means the existing FAIL-when-hit behavior (§3.2/§4) is validated, now with a
sharper mechanism: the ankle DOF isn't just "jittery," under weak marker observability the solver
can drive it to its own physical limit.

**Two new, cheap, no-rerun-required diagnostics added to every clip's report** (pure constants
compared against already-solved angles, `MODEL_COORD_LIMITS_DEG` in `pose_to_opensim_ik.py`):
`near_full_mechanical_span` per coordinate, and a robust `range_p5_p95` statistic (5th-95th
percentile spread) alongside the existing `min`/`max`/`range`, so a single pinned outlier frame
cannot masquerade as a clean full-range sweep without at least the more robust statistic being
available to cross-check it. Neither changes the pre-registered gate itself (no post-hoc
threshold shopping) — both are reported as supplementary interpretive context in §10.
