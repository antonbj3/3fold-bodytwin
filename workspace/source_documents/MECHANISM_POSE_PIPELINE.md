# MECHANISM POSE PIPELINE — IG corpus → MediaPipe pose → OpenSim .trc

Written 2026-07-21, updated same day after scale-up. Status: **pilot batch of 5 clips verified end-to-end
(§1-5) → left/right handedness open item resolved via synthetic control (§9) → full 264-clip priority corpus
extracted, 264/264 processed, 228/264 pass the full gate suite (§8).** Every number below was measured this
session (printed by the script, read directly from its JSON output, or computed by a machine-checked gate) —
not estimated.

## 1. What worked (verified live, exact commands)

- **Venv:** `~/projects/cad-to-simulation-I/.venv-humancap/bin/python3` — confirmed `mediapipe==0.10.32`,
  `cv2==4.13.0`, `numpy==2.4.1`, and `mediapipe.tasks.python.vision.PoseLandmarker`/`PoseLandmarkerOptions`
  all import cleanly. bodytwin itself has no such venv (by design, per the handoff).
- **Model:** `/media/anton/183E48713E484A48/cad-to-simulation_video/geoprobe_scratch/weights/pose_landmarker_full.task`
  (9.4 MB) — present exactly where `BODYTWIN_HANDOFF.md` said.
- **Reused, not reinvented:** the exact `PoseLandmarkerOptions`/`BaseOptions`/`detect_for_video` call shape
  from `cad-to-simulation-I/scripts/human_capture_chain_v1.py` (read-only reference — not modified). That script's
  own docstring had already flagged (2026-07-16) "world-landmarks are hip-centered relative → no metric scene
  scale"; this session's own measurement (§4.2) sharpened that finding considerably.
- **Cross-check tool:** this repo's own `.venv-msk` (`opensim==4.6-2026-06-22-85aaf64`, built by the sibling
  mechanism lane working the OpenSim/band side concurrently) — used **read-only**, only to load-test the `.trc`
  output via `opensim.TimeSeriesTableVec3`. Not modified.

## 2. Pipeline

`scripts/msk/pose_extract.py` (bodytwin-only, gitignored output in `data/msk_pose/`):

```
cad-to-simulation-I/.venv-humancap/bin/python3 scripts/msk/pose_extract.py \
    --video <clip.mp4> --group <corpus group> --quality <corpus quality>
```

Per clip, decodes frames with `cv2.VideoCapture`, runs `PoseLandmarker.detect_for_video` (VIDEO mode,
`num_poses=2`; if >1 person detected — e.g. a coach in a technique_tutorial clip — keeps the **larger 2D-bbox**
pose as a documented prominence heuristic, not a verified re-identification), and emits two artifacts per clip
in `data/msk_pose/<stem>/`:

1. **`<stem>.pose_raw.json`** — the primary, most-defensible artifact: per-frame timestamp, detection count,
   and ALL 33 MediaPipe landmarks (both 2D image-normalized and 3D "world") with visibility, completely
   untransformed. Nothing downstream can be more trustworthy than this file.
2. **`<stem>.trc`** — an OpenSim marker file, 17 named markers (Nose/Shoulders/Elbows/Wrists/Hips/Knees/
   Ankles/Heels/FootIndex — the same joint chain the sibling `MECHANISM_MSK_BUILD_PLAN.md`'s LaiArnold model
   pick uses), built from the world landmarks via a **per-clip geometric rotation** (§4.1) into OpenSim's
   Y-up convention.

`--verify <stem>` re-reads an already-produced output and runs machine-checkable (never eyeballed) gates;
`--rebuild-trc <stem>` regenerates just the `.trc` from the saved raw JSON without re-running MediaPipe
(used heavily below to iterate the axis fix cheaply).

## 3. Batch processed (5 clips, quality≥4, group ∈ {plyometrics_jumps/*, weighted_lifts/squat})

Picked for group+account diversity (3 IG accounts, both target groups, durations 7–19s):

| clip | group | account | frames | det_frac | multi-pose frames | wall (s, CPU) |
|---|---|---|---:|---:|---:|---:|
| 2025-09-26_DPEOqLVkav7 | weighted_lifts/squat | get.gymnast.fit | 368 | 1.0000 | 0 | 11.79 |
| 2025-09-22_DO6XLmikQvI | plyometrics_jumps/technique_tutorial | get.gymnast.fit | 581 | 1.0000 | 0 | 19.00 |
| 2026-03-03_DVb7W3Oj2nQ | plyometrics_jumps/plyometric_drill | theathleticpt | 539 | 0.6308 | 9 | 16.12 |
| 2026-05-30_DY-yVTeuXyR | plyometrics_jumps/plyometric_drill | samuelwatson__ | 560 | 0.8036 | 12 | 15.97 |
| 2026-03-22_DWNGMkfj7g7_5 | plyometrics_jumps/technique_tutorial | theathleticpt | 213 | 0.9343 | 0 | 6.49 |

CPU-only (XNNPACK), ~30 fps effective throughput. Total output 5.8 MB (`data/msk_pose/`, gitignored).
`det_frac<1` clips are genuinely harder shots (fast motion / partial framing), not a pipeline defect —
handled correctly via the NaN-sentinel gap convention (§4.3).

**Primary verify target** (explicitly named in the task, `weighted_lifts/squat`): `2025-09-26_DPEOqLVkav7`.
Sample `.trc` rows (`Frame#, Time, Nose(X,Y,Z), LShoulder(X,Y,Z), ...`, units m, Y-up):
```
1	0.000000	0.030377	0.768350	-0.056520	-0.151383	0.510524	-0.067183	...
2	0.033000	0.039392	0.769498	-0.047106	-0.144524	0.512999	-0.055612	...
3	0.067000	0.034202	0.765899	-0.061835	-0.148615	0.510942	-0.064670	...
```
Raw JSON sample (frame 0, nose world landmark `[x,y,z,visibility]`): `[-0.751377, 0.019966, -0.171792, 0.9858]`.

## 4. Findings this session (machine-measured, not assumed — the actual work, not just plumbing)

### 4.1 A fixed MediaPipe→OpenSim axis remap does NOT generalize — replaced with a per-clip geometric rotation

First attempt: `OpenSim_Y = -mediapipe_y` (assumed from the 2D image-landmark convention). Checked against a
physically-grounded external anchor (head above feet, gap must be >0.3 m) on the squat clip: **passed the sign
check but the magnitude was only 0.24 m** (should be ~1.4 m for a standing adult) — investigated rather than
accepted (per the standing "don't declare success without forcing it" discipline). Root cause: measured the
hip→shoulder unit vector (the spine direction, which must be close to "up") directly — it was **99% along the
MediaPipe *X* axis**, not Y, for this clip. Second attempt: assume X is vertical instead. **Tested across all
5 batch clips** (the correct diverse-instance-space test) — **falsified**: 2/5 clips X-dominant (~0.99), 3/5
Y-dominant (~0.96–0.98), one clip unstable even within itself (dominant axis flipped frame-to-frame 60% of the
time). MediaPipe's `pose_world_landmarks` frame orientation relative to gravity is **not fixed across clips**
(varies with the detected body/ROI orientation per clip) — confirmed from data, not from (recalled, and in
this case wrong) documentation; the installed package's own docstrings (`mediapipe/tasks/python/components/
containers/landmark.py`) don't specify the coordinate system at all.

**Fix** (`compute_up_rotation()` in `pose_extract.py`): derive the "up" direction from each clip's own measured
anatomy — the same principle OpenSim's own static-calibration trial already uses. `e_y` = mean(mid_shoulder −
mid_hip) direction; `e_x` = mean(right_hip − left_hip) direction with `e_y` projected out (Gram-Schmidt);
`e_z = e_x × e_y`. Applied as a per-clip rotation matrix (not a coordinate swap) to every landmark.
**Re-tested on the same 5-clip falsifier**: `head_above_feet_gap_m` now **0.888–1.457 m across all 5 clips**
(physically plausible for nose-to-ankle distance in varied athletic postures) — the fix generalizes where the
guess did not. Honest residual limit: assumes clip-mean spine ≈ gravity-vertical (degrades for sustained
non-upright postures, e.g. floor exercises); the in-plane X/Z (left-right/forward-back) assignment is an
orthonormal, right-handed construction but **not independently checked against a left/right ground truth**.

### 4.2 MediaPipe world landmarks are hip-origin-centered PER FRAME — no global translation signal at all

Investigating why the squat clip's hip vertical range was only 0.033 m (implausibly small for a squat rep)
led to a bigger and more important fact than the axis issue: **mid-hip position measured at |position| < 0.002 m
from the coordinate origin in literally every frame**, confirmed on both the squat clip and a jump clip (raw,
pre-rotation data). MediaPipe's `pose_world_landmarks` are **re-centered on the hip in every single frame
independently** — this removes ALL whole-body translation (squat depth, jump height in real space) from the
representation *by construction*. What survives is body pose *relative to the hip* — which is why knee/ankle
ranges are large and real (e.g. ankle-relative-to-hip range up to 1.49 m on jump clips, correctly reflecting
leg tuck/extension) while hip-relative-to-itself is necessarily ~0.

**Consequence for downstream OpenSim use:** this `.trc` supports relative joint-ANGLE kinematics correctly, but
**cannot** recover global pelvis trajectory (how high the person actually jumped, how deep the squat went in
real space) — that signal isn't present in `pose_world_landmarks` for any clip, not just these 5. A separate
global-scale/translation recovery (e.g. from the 2D image-space bounding box trajectory, in `image_landmarks`
in the raw JSON) would be needed, and per the sibling `human_capture_chain_v1.py`'s own already-verified
finding, is gauge-limited without a metric anchor in view for an arbitrary uncalibrated clip.

### 4.3 A real bug caught by an EXTERNAL cross-check (not caught by self-validation)

Built a `--verify` gate suite (file existence, row/marker-count consistency, non-degenerate-motion void-floor
check, head-above-feet anchor) — all passed on my own hand-rolled `.trc` parsing. Then loaded the same files
through **this repo's own `.venv-msk` OpenSim install** (`opensim.TimeSeriesTableVec3`) as a genuinely
independent parser: **3 of 5 clips FAILED** — `"Unexpected number of columns... Expected = 53. Received = 52"`
— exactly the 3 clips with `det_frac < 1.0`. Root cause, confirmed with a minimal synthetic 2-marker test file:
a fully-undetected frame's row ends in a long run of blank tab-separated fields; OpenSim's C++ TRC reader
trims trailing whitespace before splitting, silently dropping the row's final column. Fix: write an explicit
`"NaN"` sentinel instead of a blank cell for undetected markers — confirmed via the same synthetic test that
`NaN` round-trips as an actual NaN in the resulting OpenSim table (correct, machine-detectable missing-data
semantics), while an all-blank row throws. **5/5 clips now load cleanly** through OpenSim's own reader, and the
per-clip NaN-row counts cross-validate exactly against the independently-computed `det_frac` (e.g.
539 rows − 199 NaN = 340 detected ≡ det_frac 0.6308, matching to 4 decimal places). This external-tool check
is now a permanent, non-skippable-when-available gate in `verify_output()` (`opensim_cross_check`), not a
one-off manual step — run `--verify` under `.venv-msk/bin/python3` to exercise it live (it degrades to a
reported "skipped" status, not a silent pass, under an interpreter without opensim).

## 5. Final verification (all 5 clips, full gate set incl. live OpenSim cross-check)

Run via `.venv-msk/bin/python3 scripts/msk/pose_extract.py --verify <stem>` (opensim import succeeds there):

| clip | trc rows/cols via OpenSim | head_above_feet_gap_m | ankle_Y_range_m | OVERALL_PASS |
|---|---|---:|---:|---|
| 2025-09-26_DPEOqLVkav7 (squat) | 368/17 | 1.4565 | 0.1287 | **PASS** |
| 2025-09-22_DO6XLmikQvI | 581/17 | 1.2057 | 0.8053 | **PASS** |
| 2026-03-03_DVb7W3Oj2nQ | 539/17 | 0.9688 | 1.2987 | **PASS** |
| 2026-05-30_DY-yVTeuXyR | 560/17 | 0.9481 | 1.4933 | **PASS** |
| 2026-03-22_DWNGMkfj7g7_5 | 213/17 | 0.8884 | 0.3143 | **PASS** |

5/5 PASS on every gate: file existence, `.trc` header/body/raw-JSON row-count consistency, marker-count
consistency, non-degenerate ankle-relative motion (void-floor 0.02 m), head-above-feet physical anchor
(>0.3 m), and the live external OpenSim-reader cross-check.

## 6. Honest open items (not solved by this pass — do not imply otherwise)

1. **No global/scene-metric scale or translation** (§4.2) — joint angles yes, pelvis trajectory / real jump
   height / real squat depth no, for any clip processed by this pipeline as-is.
2. ~~Left/right and forward/back (in-plane X/Z) handedness of the per-clip rotation is unverified~~ **RESOLVED
   2026-07-21 — see §9.** Synthetic known-pose control: `compute_up_rotation`/`apply_rotation` correctly and
   invariantly recover left/right and up/down across 8 diverse input orientations (max residual 3e-6 m,
   wide-margin labeled-sign checks 0.18-1.02 m, det(R)=+1 always), extended to a det(R)=+1.0-exactly check on
   all 264 real corpus clips. Residual, DECORRELATED, NOT closed by this control: MediaPipe's own raw L/R
   semantic convention (anatomical vs. mirrored) is a property of the upstream model, not of this rotation —
   see §9 for the precise scope and the derived (and confirmed) reflection-case prediction.
3. **Multi-person "prominence = larger bbox" heuristic remains unverified against ground truth — now measured
   at full-corpus scale (§8)**: 60/264 clips (22.7%) had ≥1 multi-pose frame, 2115 multi-pose frames total (vs.
   the pilot's 2/5 clips, 21 frames) — a much larger exposure surface than the 5-clip pilot suggested; still
   not independently checked that the selected pose is the athlete rather than e.g. a closer-framed coach.
4. **Corpus group/quality labels are auto-labeled and unverified — CONFIRMED MATERIAL at scale (§8)**: 31/264
   quality≥4-labeled clips (11.7%, all `get.gymnast.fit`, all exactly 4.0s/120 frames) contain NO visible
   athlete at all — confirmed by direct pixel measurement (mean brightness 3.5-5.1/255, frame-to-frame diff
   <0.005: static, near-black "quote card" content), not a detector failure (null even at
   `min_pose_detection_confidence=0.1`). Upstream corpus-labeling defect, not a `pose_extract.py` bug — the
   pipeline correctly reports `det_frac=0.0` rather than fabricating a detection — but it means "264 priority
   clips" overstates usable athletic content by ~11.7%; true usable count is 233.
5. **MediaPipe's own coordinate-system documentation was not found in the installed package**; the axis
   findings above are from this session's direct measurement, cross-clip, not from a citable spec. A live fetch
   of Google's current MediaPipe Pose Landmarker docs this session DID confirm the 33-landmark index→name
   ordering exactly matches `LANDMARK_NAMES` (0=nose … 11=left_shoulder … 23=left_hip … 32=right_foot_index) —
   ruling out a transcription bug there — but confirmed the docs are silent on the anatomical-vs-mirrored L/R
   question (see §9); still not independently citable.
6. **NEW, only visible at full-corpus scale (§8): the "clip-mean spine ≈ vertical" assumption in
   `compute_up_rotation` (already documented as an honest limit in the code's own docstring) breaks for
   SUSTAINED folded/inverted postures** — 5/264 clips (1.9%, all `get.gymnast.fit` core-compression / mobility
   drills: hip compression, shoulder compression, spine/hip mobility) have real detections (det_frac 0.91-1.0)
   but FAIL the head-above-feet external anchor (gap as low as -0.17 m, i.e. head measured below feet on
   average). Not a bug — the pre-existing documented degradation case, now quantified for the first time.

## 7. Scale-up to the full priority set — COMPLETED 2026-07-21 (see §8 for results, §9 for handedness)

264 clips in the corpus match `quality≥4 AND group ∈ {plyometrics_jumps/*, weighted_lifts/squat}` (out of 893
total rows; the `weighted_lifts/squat`-only definition, not all of `weighted_lifts/*`, was confirmed empirically
against the live corpus file — it is the one that reproduces exactly 264 = 236 plyometrics_jumps/* + 28
weighted_lifts/squat) across 5 accounts (`realgame.athletics` 127, `get.gymnast.fit` 114, `theathleticpt` 20,
`samuelwatson__` 2, `iambburns` 1 — confirmed to match this session's original prediction exactly). Per the
recommendation below, honest-item #2 (L/R handedness) was resolved FIRST via a synthetic control (§9, PASS) —
no code fix was needed, so no bulk re-run risk materialized. The full 264-clip set was then extracted (§8).
Honest item #1 (no global translation) remains open and out of scope for this pass; see §10 for the current
recommended next step.

## 8. Full corpus batch results (264/264 processed, 2026-07-21)

**Orchestration** (`scripts/msk/run_pose_batch.py`, new this session): a thin subprocess-per-clip harness around
the UNMODIFIED, already-verified `pose_extract.py` command line — no pipeline logic reimplemented. 8 parallel
workers (RAM is the binding constraint per this repo's `COORDINATOR.md` §4 shared-machine convention; measured peak
RSS 294 MB/process on a real clip → even 20-way parallelism would use <6 GB, nowhere near the 12 GB floor;
`resource_gate.py check --ram-gb 4` → `OK` before launch). Sanity-checked on 3 then 20 clips before committing
to the full remaining set (`--verify`-style discipline applied to the orchestrator itself, not just the pipeline).

**Completion:** 264/264 manifest clips have a valid `.trc` + `.pose_raw.json` (264/264, zero missing). Every
individual subprocess invocation exited 0 — **zero crashes, zero timeouts, zero exceptions** across the whole
corpus. 337 MB total output. Wall-clock: ~13 min this session for the 259 newly-run clips at 8-way parallelism
(sum of individual clip wall_s = 5656s ≈ 94.3 min single-clip-equivalent CPU time, matching this doc's original
~70-90 min single-threaded estimate; 8x parallel brought elapsed wall time down to ~660s for the bulk of it).
Account breakdown confirmed exactly as predicted: `realgame.athletics` 127, `get.gymnast.fit` 114,
`theathleticpt` 20, `samuelwatson__` 2, `iambburns` 1.

**Detection quality:**
| metric | value |
|---|---|
| clips with det_frac > 0 (genuine detections) | 233/264 (88.3%) |
| — median det_frac (non-zero) | 1.0000 |
| — mean det_frac (non-zero) | 0.9717 |
| — ≥0.95 det_frac | 200/233 (85.8%) |
| — min det_frac (non-zero; still passes full gate) | 0.1621 (`2026-07-04_DaYjAaquUq0`, hard/fast clip) |
| clips with det_frac == 0 | 31/264 (11.7%) — see below, content issue not pipeline defect |
| clips with ≥1 multi-pose frame | 60/264 (22.7%), 2115 multi-pose frames total |

**Full gate-suite pass rate** (`pose_extract.py --verify`, run under `.venv-msk` for the live external OpenSim
cross-check, machine-checked not eyeballed): **228/264 (86.4%) OVERALL_PASS.** The independent OpenSim reader
(`opensim.TimeSeriesTableVec3`) loaded **264/264 (100%)** `.trc` files without error — every artifact is
structurally well-formed, including the 31 degenerate all-NaN ones (confirms the NaN-sentinel fix from §4.3
holds at full scale). The 36 non-passing clips resolve into exactly two known, root-caused, non-overlapping
failure signatures — no unexplained failures:
- **31 clips** fail `det_frac_gt_0` + `ankle_Y_range_above_noise_floor` + `head_above_feet_pass` together — the
  zero-detection "quote card" content clips (root-caused below), correctly fail-closed (void-floor sweep: the
  gate does NOT silently pass a degenerate all-NaN output).
- **5 clips** fail `head_above_feet_pass` ONLY (real detections, det_frac 0.91-1.0) — the folded/inverted-posture
  case, see open-item #6 above.

**Root cause of the 31 zero-detection clips** (forced via OODA, not assumed from captions): all 31 are
`get.gymnast.fit`, all exactly 120 frames / 4.0s. Direct pixel measurement (mean grayscale brightness,
frame-to-frame absolute difference) on multiple of these clips: mean brightness 3.5-5.1 / 255 (≈1.4-2.0%),
frame-to-frame diff max 0.0045 — i.e. static, near-black video for the full clip, not a detector failure (null
result persists even at `min_pose_detection_confidence=0.1`, 5x below production). This is an **upstream
corpus auto-labeling defect** (these are not athletic-movement footage despite their `quality=4` /
`plyometrics_jumps/technique_tutorial` label) — confirms the pre-existing honest item #4 concretely, at
exactly the 11.7% level, for the first time. **Not a `pose_extract.py` bug**: the pipeline's behavior here
(report `det_frac=0.0`, still emit a structurally valid all-NaN `.trc`, never fabricate a detection) is exactly
correct.

## 9. Handedness synthetic control — RESOLVES open item #2 (verdict: PASS, no code fix needed)

**Method** (`scripts/msk/handedness_synthetic_control.py`, new this session; imports and exercises the REAL
`pose_extract.py` functions — `write_trc`, `compute_up_rotation` — and parses the actual `.trc` artifact they
write, not an in-memory shortcut): a synthetic 33-landmark skeleton was built in a DECLARED ground-truth frame
(GT_X = subject's right, GT_Y = up, GT_Z = subject's front), deliberately ASYMMETRIC on BOTH X and Z (right arm
raised forward+up, left arm down at side; right leg forward+lifted, left leg planted) so an axis swap or sign
flip on either axis would be caught, not just Y.

**Geometric prediction, derived before measuring** (§4.1's own finding — MediaPipe's raw per-clip world-frame
orientation is arbitrary — means the "raw" input must be modeled as `Q @ GT_point` for some per-clip orthogonal
`Q`). Since `(Qa)×(Qb) = det(Q)·Q(a×b)`, `compute_up_rotation`'s cross-product construction implies its output
should be **invariant to any proper rotation Q** (det=+1), and a **reflection** (det=-1) should flip only the
recovered Z (front/back), leaving X (left/right) untouched.

**Measured, against pre-registered thresholds:**
| gate | result |
|---|---|
| proper-rotation invariance (8 orientations: identity, 2 axis-cycles, 5 random SO(3)) | max residual 3×10⁻⁶ m (< 10⁻⁴ m threshold; residual is rounding noise from two independent 6-decimal quantization steps in the harness itself, not a geometric error — confirmed by threshold-recalibration, not by loosening past the real-signal scale, which is 4-5 orders of magnitude larger) |
| det(R) == +1 for every proper-Q synthetic clip | exactly 1.0 in all 8 cases |
| labeled L/R sign checks (RHip>LHip on X, RShoulder>LShoulder on X, RWrist>LWrist on Y+Z, RAnkle>LAnkle on Z) | correct with wide margins (0.18-1.02 m), IDENTICAL across all 8 orientations |
| deliberate reflection case (det(Q)=-1): predicted X/Y preserved, Z flipped | confirmed to 0.0 m error (exact) |
| independent OpenSim reader cross-check on the synthetic `.trc` | `loaded_ok`, 5 rows / 17 cols |
| **extended to all 264 REAL corpus clips**: det(R) | exactly 1.0 (within 1e-6) on every single one, including the 31 identity-fallback (zero-detected-frames) cases |

**Verdict: PASS.** No flip found; no code change was needed or made to `pose_extract.py`. The per-clip rotation
correctly and invariantly recovers left/right and up/down regardless of MediaPipe's arbitrary per-clip raw
orientation, confirmed both synthetically (with known ground truth) and at full real-corpus scale (with the
det(R)=+1 structural invariant).

**Precise residual scope (do not overclaim beyond this):** this control verifies `pose_extract.py`'s OWN
geometric transform, which is the item docs §6 flagged. It does NOT independently verify MediaPipe's upstream
raw semantic convention (whether its `left_hip`/`right_hip` labels are anatomical or image-mirrored) — a
live fetch of Google's current MediaPipe Pose Landmarker documentation this session found the landmark
index→name ordering exactly matches this file's `LANDMARK_NAMES` (external cross-check, ruling out a
transcription bug) but found the docs explicitly silent on the anatomical-vs-mirrored question. That is a
property of the upstream model shared identically by the already-in-production `human_capture_chain_v1.py`
reference, decorrelated from anything this rotation could introduce or fix. The reflection-case test above
additionally shows that IF MediaPipe's raw per-clip frame were ever an improper (mirrored) transform — which
would require a monocular depth-sign ambiguity internal to MediaPipe, not any real physical camera motion —
the resulting corruption would be to front/back (Z) specifically, not left/right (X). Closing this fully would
require an end-to-end test with real images through the real MediaPipe detector (a known asymmetric marker in
a synthetic photo) — judged disproportionate effort this session relative to the doc's own flagged scope, and
not attempted.

## 10. Next step (revised after §8-9)

1. **Address the 31 mislabeled zero-content clips (§8)** before any downstream consumer treats "264" as "264
   usable athletic clips" — either fix the corpus auto-labeler upstream or filter `det_frac==0` clips out of
   the priority set (cheap: already flagged in each clip's own `summary.json`).
2. **Marker augmentation + OpenSim Scale+IK** (per `docs/MECHANISM_MSK_BUILD_PLAN.md` §5's already-proven chain:
   pose → Pose2Sim triangulation/augmentation → Scale+IK → `.mot`) is the natural consumer of this `.trc`
   corpus — note the sibling `video2kin/venv_pose2sim` Pose2Sim environment is separately documented as
   "BROKEN IN PLACE" (that doc's own finding, not this session's) and would need re-provisioning first (that
   same doc's §5 step 1 already has the exact `uv venv` recipe).
3. **Before wiring into IK**, decide scope on honest item #1 (no global translation/scale) — IK on the current
   `.trc` will silently produce a non-moving pelvis regardless of real jump/squat depth, for every one of the
   264 clips, not just the original 5.
4. **Multi-pose re-ID heuristic (open item #3)** is now known to engage on 22.7% of the corpus (60/264 clips,
   2115 frames) — worth a targeted spot-check (e.g. cross-referencing 2D bbox track continuity) before trusting
   those clips' `.trc` in any per-athlete analysis, given the exposure is ~3x larger than the pilot suggested.
