# MECHANISM FOOT CLOSE-UP POSE — foot-specific keypoint extraction vs the 3-landmark floor (2026-07-21)

Answers the operator's #1 priority (foot/ankle) with the operator's own framing correction applied:
caption keywords are a bad proxy for what a clip *shows*; the real wall is SENSING. General-body
MediaPipe Pose gives exactly 3 foot landmarks per side (`ankle`, `heel`, `foot_index` — verified live
against `scripts/msk/pose_extract.py`'s own `LANDMARK_NAMES`/`MARKER_SUBSET`, indices 27-32), so
midfoot/toe/arch motion is architecturally invisible to it regardless of clip content (already
established in `docs/MECHANISM_FOOT_FIDELITY_PLAN.md` §2 and `docs/MECHANISM_FOOT_DRILL_CORPUS.md`
§4/§9 — not re-litigated here). This doc tests the proposed unblock directly: does CLOSE-UP framing
plus a FOOT-SPECIFIC extraction method (MediaPipe's HAND model repurposed on a bare foot, and classical
contour CV) recover more than that 3-landmark floor, and does whatever it recovers anchor the
multi-segment foot's still-placeholder midtarsal/subtalar DOFs (`docs/MECHANISM_FOOT_MULTISEGMENT.md`)?

**Script:** `scripts/msk/foot_closeup_pose.py` (re-runnable, deterministic — reproduced byte-identical
verdicts across 3 runs this session; runs under `cad-to-simulation-I/.venv-humancap`, see its own
docstring). **Evidence:** `data/msk_pose_closeup/foot_closeup_pose/foot_closeup_pose_summary.json` +
forensic debug frames in `data/msk_pose_closeup/foot_closeup_pose/debug_frames/`.

---

## 1. Fresh YouTube close-up fetch — attempted, bot-blocked (measured, not assumed)

Per the task, 4+ search families were tried via `yt-dlp` ("barefoot foot mobility drill close up",
"toe splay exercise slow motion", "foot intrinsic muscle exercise close up", "ankle CARs barefoot",
plus 2 follow-ups) against **6 distinct player-client configurations** (default/unspecified, `android`,
`ios`, `tv`, `web`, `tv_embedded`), with backoff sleeps up to 15s between attempts, over roughly
15-20 minutes:

- **19 explicitly-logged video-ID extraction attempts** each returned `ERROR: [youtube] <ID>: Sign in
  to confirm you're not a bot` (one additionally hit an age-gate variant of the same wall), plus 5
  further `ytsearch8` batches whose per-ID errors were not individually captured in this session's
  terminal output (still real attempts, still zero usable results).
- **Exactly 1 transient success**: `LhgaiYrjjRI` ("Barefoot Running Exercise: Toe Speading", 30fps,
  49s) returned real metadata via the `android` client on first contact. It did **not** reproduce:
  the *same video ID* was retried 3 more times (android-client download attempt, `web` client with
  `--sleep-requests 2`, `tv_embedded` client) across the following ~5 minutes and was blocked every
  time.
- This pattern (one query-independent, client-independent, video-ID-independent success that
  immediately closes again) is diagnostic of an **IP-reputation-level block**, not a fixable
  query/client/video issue — confirmed independently by the coordinator's own measurement mid-task,
  convergent with this session's own finding. Further retrying would be pure wasted runs against a
  wall neither query rewording nor client spoofing opens.
- Two alternative free-footage sources were reachability-checked (`curl -o /dev/null -w '%{http_code}'`)
  before the coordinator's pivot instruction arrived: `pexels.com` and `pixabay.com` both **403**
  (their own bot-wall), `vimeo.com` **200** (reachable, not pursued further once the pivot to the
  existing corpus was directed).

**Unblock recipe (concrete, for when the operator wants to complete this)**: this sandbox has no
browser profile yt-dlp can read cookies from. Two real paths, both requiring operator action outside
this session (same class as this project's other standing human-identity asks):
1. `yt-dlp --cookies-from-browser <chrome|firefox|edge|...> "URL"` — run directly on a machine with a
   real, currently-logged-in YouTube browser session.
2. Export a `cookies.txt` from a logged-in browser (e.g. the "Get cookies.txt LOCALLY" extension),
   copy it into this environment, then `yt-dlp --cookies /path/to/cookies.txt "URL"`.

---

## 2. Pivot: closest-to-closeup frames already in the corpus, selected by measurement not eyeballing

Per the coordinator's redirect, used `docs/MECHANISM_FOOT_DRILL_CORPUS.md`'s own shortlisted clip
families (`2025-07-31_DMwHVZLsLe9*`, 7 files; `2026-06-16_DZoF3vaDA9j*`, 8 files — both
`realgame.athletics`, already pose-extracted in `data/msk_pose/`). **Frame selection is a machine
ranking, not a guess**: every already-computed frame's 3-point foot-landmark
(`ankle`/`heel`/`foot_index`) image-space bounding-box AREA (a direct fraction of frame area, min
landmark visibility ≥0.5) was ranked across **2,136 candidate (frame, side) samples**. The single
largest anywhere is `2026-06-16_DZoF3vaDA9j_6` frame 91 (0.777% of frame area, 94×57px raw
3-point bbox) — **this precisely quantifies, rather than just asserts, `MECHANISM_FOOT_DRILL_CORPUS.md`
§5's finding that these are wide/medium shots**: even the single best-framed foot moment in ~2,100
samples occupies well under 1% of the frame.

Top candidates in both families were extracted (`cv2.VideoCapture`, sequential read to the target
index — the exact same frame-numbering convention `pose_extract.py` used, not an independent ffmpeg
filter that could drift by one), cropped with a fixed, pre-registered generous pad (half-width/height
= max(1.6× the 3-point bbox, 7.5% of frame), and **forensically viewed** (a per-frame visual fact,
not inferred from caption/account — same discipline `MECHANISM_FOOT_DRILL_CORPUS.md` §5 already used).
This surfaced an important, non-obvious split **within the same account's own shortlisted clips**:

| tag | clip | barefoot or shod (forensic) | crop px | bbox frac of frame |
|---|---|---|---|---|
| `dzof6_f91` | `DZoF3vaDA9j_6` f91 | **SHOD** — black minimalist sock-shoe, visible logo | 204×134 | 0.777% |
| `dzof_f70` | `DZoF3vaDA9j` f70 | **SHOD** — white athletic shoe | 236×151 | 0.500% |
| `dzof6_f93` | `DZoF3vaDA9j_6` f93 | **SHOD** — same sock-shoe as f91 | 307×127 | 0.513% |
| `dmw5_f122` | `DMwHVZLsLe9_5` f122 | **BAREFOOT** — toes curled on a gymnastics bar | 137×190 | 0.367% |
| `dmw5_f112` | `DMwHVZLsLe9_5` f112 | **BAREFOOT** — same clip, 10 frames earlier | 141×173 | 0.346% |
| `dmw6_f20` | `DMwHVZLsLe9_6` f20 | **SHOD** — sneaker, mid-jump, motion-blurred | 127×156 | 0.280% |
| `dmw7_f76` | `DMwHVZLsLe9_7` f76 | **SHOD** — shoe on an inclined wooden board | 183×144 | 0.288% |

Only **2 of the 7** best-framed candidates across both shortlisted families are actually barefoot
(both from the *same* clip, `DMwHVZLsLe9_5` — a gymnastics/calisthenics bar-grip drill, toes curled,
not splayed). This is a second, independent limiter on top of small-in-frame: **an ankle/foot-drill
compilation account mixes footwear across its own "7 exercises," and shod is more common than
barefoot even among the best-framed moments.** Neither the classical contour method nor the
hand-model repurposing can see through a shoe — this is a content property of the footage, not a bug.

---

## 3. Method 1 — MediaPipe HandLandmarker (21kp) repurposed on the bare foot

**Pre-registered thresholds** (stated before scoring): "meaningful" detection = score ≥0.5 (the
model's own default `min_hand_detection_confidence`, not tuned post-hoc); "any raw hypothesis"
(diagnostic only, never counted as a real detection) = score ≥0.05.

### 3.1 Resolution-confound control (forced BEFORE trusting any foot result)

A real hand at full resolution (1080×1920, a frame with an independently pre-existing confirmed
detection, re-verified live against `scripts/hand_leg/raw_tracks/VID_20260716_231152638_hands.json`)
scores **0.934**. The **same real hand, downscaled to 106×190px** — comparable to the foot crops
below — still scores **0.903**. This decisively rules out "low resolution" as an explanation for any
foot null result: if the hand model fails on a foot crop, it is because a foot isn't hand-shaped
enough, not because the crop is small.

### 3.2 Per-frame results

| tag | barefoot | hand meaningful (raw/up4x) | hand score (best) | skin-mask frac | contour defects >0.15 |
|---|---|---|---|---|---|
| `dzof6_f91` | shod | 0/0 | — | 0.124 | 0 |
| `dzof_f70` | shod | 0/0 | — | 0.060 | 0 |
| `dzof6_f93` | shod | 0/0 | — | 0.077 | 0 |
| `dmw5_f122` | **barefoot** | 1/1 | 0.583 | 0.885 | 2 |
| `dmw5_f112` | **barefoot** | 1/1 | 0.832 | 0.916 | 2 |
| `dmw6_f20` | shod | 0/1 | **0.922** | 0.941 | 2 |
| `dmw7_f76` | shod | 0/0 | — | 0.169 | 5 |

**The real, positive finding**: both barefoot frames get a "meaningful" (≥0.5) hand-shaped detection.
Rendering the 21-keypoint overlay and reading exact pixel coordinates (not eyeballing a thumbnail —
machine-measured) shows the 4-finger fan lands **on the visible toe row**, with `INDEX_MCP`→`PINKY_MCP`
x-positions monotonically increasing and evenly spaced (deltas 8.3/9.0/9.5px on `f112`, 7.6/8.5/9.3px
on `f122` — consistent across 2 independent frames of the same movement) and 3 of the 4 fingertips
landing close to 3 distinct visible toe tips. **This is genuinely more per-toe spatial information than
the 3-landmark floor has at all** (which has zero toe-individuating landmarks, full stop).

**The real, load-bearing caveat, found by forcing the adversary rather than trusting the score**: on
`dmw6_f20` (a SHOD, motion-blurred sneaker), the upscaled crop gets a hand detection at **0.922
confidence — higher than either true barefoot detection (0.583, 0.832)**. Machine-computed:
`confidence_score_separates_true_from_false_positive = false`. A naive "the hand model detected
something confidently" filter would happily accept this shoe as a hand. This is exactly the adversary
the watertight method requires forcing, not skipping — found and reported, not hidden.

**A second, forced OODA pass (do not stop at the honest-negative-shaped "score doesn't work") turned
up a geometry-derived signal that DOES separate them this session**: the mean MCP→TIP chain length
(index/middle/ring/pinky), divided by the INDEX_MCP↔PINKY_MCP "palm width" — grounded directly in the
real anatomical fact that toes are much shorter relative to forefoot width than fingers are relative to
palm width:

| condition | ratio (n) |
|---|---|
| real hand, hi-res + downscaled positive control | 0.849, 1.017 |
| **true barefoot detections** (2 frames × raw/up4x) | **0.070, 0.132, 0.140, 0.140** |
| **shod false-positive** (`dmw6_f20`) | **1.558** |

The false positive's ratio (1.56) sits *above* even the real hand's own range (0.85-1.02), while every
true barefoot detection sits 6-14× lower (0.07-0.14). `ratio_signal_separates_true_from_false_this_session
= true` (machine-checked: max(barefoot)=0.140 < min(shod-false-positive)=1.558). **This is a promising
LEAD, explicitly NOT a validated discriminator** — it rests on n=3 total examples (1 physical hand at 2
resolutions, 2 frames of one barefoot clip, 1 false positive) — nowhere near the diverse instance-space
the watertight method requires before accepting it as a real filter. Flagged as the concrete next cheap
test (does it hold on more shod/barefoot examples?), not claimed as proven.

### 3.3 Classical contour CV — clean limb-silhouette segmentation, no toe individuation on this footage

Skin-color (YCrCb, fixed range, not tuned per-frame) segmentation cleanly separates the whole bare
limb from a dark background on both barefoot frames (mask fraction 0.885-0.916) and correctly
EXCLUDES the shoe material on shod frames (`dzof6_f91`'s 0.124 mask fraction is the sliver of visible
bare ANKLE/shin skin *above* the shoe line, verified by viewing the mask directly — not the shoe
itself). The fitted principal axis (`cv2.minAreaRect`) sits at -90° (i.e., aligned with the leg's long
axis) on every clean barefoot/limb case — a real, usable **limb/ankle-flexion axis measurement**, same
class of signal the existing IK pipeline already gets from the 3-landmark floor, not new information.

**The 2 candidate convexity defects on each barefoot frame do NOT individuate toes** — forensically,
the contour traces the *entire* leg+curled-foot as one blob (toes bunched against the bar, not
splayed), and the defects sit at the ankle-narrowing/heel-bulge transition, not toe-to-toe valleys.
This is an honest **inconclusive-due-to-content**, not a clean negative for the method: this specific
footage's pose (toes curled around a bar) never presents splayed, individually-silhouetted toes for
the contour method to individuate in the first place. A genuine "toe splay drill" clip (the original
search target) would be a fundamentally more favorable test of this specific method — untested this
session because YouTube is blocked (§1).

A second, generic Otsu foreground segmentation (works regardless of skin tone, honestly labeled as
recovering a foreground OUTLINE, not a toe boundary) gives a coarser but shoe-compatible silhouette —
useful only for a gross axis, same ceiling as above. **A forced adversary caught here**: `dmw7_f76`'s
skin-color mask produced 5 "defects" with unusually large depths (0.22-1.23 normalized) — forensically
this mask is a degenerate thin cross-hair pattern, not a coherent blob (visually confirmed), consistent
with the wooden incline board's tan color false-triggering the fixed skin-color range. **Not toe
signal — a known, disclosed false-positive mode of naive skin-color segmentation on wood-toned
backgrounds**, exactly the kind of confound a raw defect-count would silently launder if not forensically
checked.

---

## 4. Comparison vs. the 3-landmark body-pose floor — stated precisely, not hyped

| | 3-landmark floor (MediaPipe Pose) | Hand-model-on-foot (this doc) | Classical contour (this doc) |
|---|---|---|---|
| Toe-tip 2D localization | **none** (0 toe landmarks) | **yes, on 2/2 tested barefoot frames** — 3-4 toe-tip-adjacent points, unvalidated accuracy | no (toes not splayed in available footage) |
| Hallux (big toe) individuated from neighbor | none | **no** — `THUMB`/`INDEX` chains collapsed onto the same point in both tested frames (within 2-5px) | no |
| Limb/ankle long-axis | yes (via ankle/heel/foot_index already) | n/a | yes, redundant with the floor |
| Arch / midfoot dorsum | none | none | none |
| Subtalar (rearfoot) motion | none | none | none |
| 3D metric accuracy | no (hip-relative only, per `MECHANISM_POSE_PIPELINE.md`) | no (2D image-plane only; `hand_world_landmarks` would be metrically WRONG — calibrated to hand bone-length ratios, not foot/toe ratios) | no (2D image-plane only) |
| Reliability | established, corpus-wide | **1 confirmed false positive at HIGHER confidence than the weakest true positive** (§3.2) — not yet safe unattended | robust when skin/background contrast is real; fails on wood-toned backgrounds |

**Honest answer to "does the foot-specific approach recover ANY midfoot/toe detail beyond the floor on
existing footage": yes, narrowly** — 2D toe-tip-adjacent localization on both tested barefoot frames,
which the 3-landmark floor cannot provide at all (it has no toe landmark of any kind). This is real,
machine-measured, and forensically checked, not hyped: it is thin (n=2, one clip, one pose), unvalidated
for metric accuracy, and comes with a demonstrated false-positive risk that a naive confidence-based
filter alone does not catch.

---

## 5. Can any of this anchor the multi-segment foot's midtarsal/subtalar DOFs? Geometric answer: no — and here is precisely why, not just "sensing is hard"

This is a clean negative, forced through the actual geometry rather than a generic "cameras can't do
3D" hand-wave:

1. **Wrong anatomical target.** Toe-tip 2D tracking (the one thing this session's foot-specific methods
   plausibly add, §3-4) informs MTP/toe flexion — a *different* joint group than subtalar (talocalcaneal)
   or midtarsal (Chopart/transverse-tarsal), which are rotations of the *hindfoot/midfoot bones
   themselves*, proximal to any toe. No per-toe keypoint scheme, however good, adds information about
   how the calcaneus rotates under the talus.
2. **Wrong viewing angle, independent of keypoint scheme.** The subtalar axis is classically described
   as markedly oblique to any sagittal-plane (side-on) camera view (recalled textbook figure, not
   re-verified this session — same caveat `MECHANISM_FOOT_MULTISEGMENT.md` §3.3 already carries for its
   own placeholder axis discussion); the standard single-camera clinical proxy for 2D rearfoot angle
   uses a **posterior** view of the heel specifically. Midtarsal assessment similarly needs a **dorsal**
   (top-down) or medial view contrasting the hindfoot's long axis against the forefoot's. Every frame
   available this session — and every search term this task specified ("toe splay," "ankle CARs,"
   "foot intrinsic muscle," all anterior/plantar/toe-focused by construction — is an
   **anterior/plantar toe-focused view**. Even a perfectly-executed fresh close-up fetch of the
   *originally-requested* search terms would not have changed this: the content class itself (toe/
   forefoot mobility drills) is the wrong vantage point for hindfoot/midfoot DOFs, a fact independent
   of whether YouTube is blocked. This is the sharper, corrected version of "need more landmarks" —
   the camera has to be pointed somewhere else entirely.
3. **Monocular 2D has no depth even with the right view.** A single camera, however well-aimed, gives a
   2D projection of a 3D rotation. The existing OpenSim IK build already demonstrates this failure mode
   empirically, not hypothetically: reusing the subtalar axis's own orientation for the midtarsal joint
   (v1, `MECHANISM_FOOT_MULTISEGMENT.md` §3.3) aliased the two joints' marker residuals together (subtalar
   RMSE jumped to 3.9-6.0°) — a real, measured instance of exactly the axis/camera-misalignment problem
   described above, from *real markers*, an even richer signal than any video keypoint scheme.

**Verdict: no candidate this session — hand-model-on-foot, classical contour, nor a hypothetically
better-executed fresh close-up fetch of the originally-requested search terms — anchors subtalar or
midtarsal.** Not because the extraction methods failed on their own terms (§3-4 shows a real, if
narrow, positive for toe detail), but because the content class and the DOF are geometrically
mismatched, independent of sensing quality.

### What sensing would actually be required (stated plainly, per the task's explicit ask)

- **Subtalar**: a posterior (heel-on) close-up view, ideally with the tibial and calcaneal long axes
  markable/identifiable (paint lines, as the classical 2D clinical rearfoot-angle method does), OR a
  multi-camera/stereo rig, OR physical markers on the calcaneus + talus/malleoli (the same ≥3
  non-collinear-points-per-segment requirement `MECHANISM_FOOT_FIDELITY_PLAN.md` §2 already derived for
  the mocap protocol — unchanged by anything in this doc).
- **Midtarsal**: a dorsal or medial close-up view showing hindfoot-vs-forefoot long-axis angulation and
  arch profile, with the same multi-view/marker/depth requirement.
- **Any of the above from a SINGLE monocular camera, in 3D, with metric accuracy**: not achievable in
  principle without either (a) ≥2 simultaneous views (stereo triangulation), (b) a depth sensor (RGBD),
  or (c) a validated, foot-specific learned 3D shape prior — no such prior was found to exist this
  session (nor in the prior sessions' literature survey, `MECHANISM_FOOT_FIDELITY_PLAN.md` §1) the way
  MediaPipe Pose attempts one for the whole body (itself already flagged unreliable for subtle
  rotations in this project's own pipeline docs).
- **The one concrete, already-partially-available lead this project has for monocular 3D relief**:
  `metric_depth_vit_small_800k.pth` / `video_depth_anything_vits.pth` already sit in this project's
  shared model weights (`/media/anton/183E48713E484A48/cad-to-simulation_video/geoprobe_scratch/weights/`),
  with a working inference environment already built elsewhere in this machine's tree
  (`cad-to-simulation-I/.venv-geoprobe`, per `reports/probes/geometry_prior_probe_1.json`) for a related
  but different (optics/vision-eyes) probe. **Not tested this session** (out of scope/budget for this
  task, and a genuinely different code path to wire up) — flagged as the most concrete next step for
  recovering *some* monocular 3D structure (arch doming, hindfoot-forefoot 3D angle) from close-ups,
  not fabricated as a working result.

---

## 6. Honest gaps

1. n=2 barefoot frames, both from ONE clip/movement (toes curled on a bar) — thin, non-diverse
   instance space for the hand-model-repurposing finding. Does not establish it generalizes to a
   standing/weight-bearing foot, a splayed-toe pose, different skin tones, or different lighting.
2. No ground-truth toe-tip location exists for any tested frame (no markers, no multi-view) — the
   "plausible placement" finding (§3.2) is a forensic/geometric plausibility check, not a validated
   accuracy measurement in pixels or mm.
3. The finger-length/palm-width ratio (§3.2) is a promising, mechanistically-grounded lead from n=3
   examples — explicitly not a certified discriminator; needs testing across many more true/false
   examples before any pipeline should rely on it.
4. Classical contour toe-individuation is untested-not-refuted for splayed-toe content — this
   session's available footage never presented that pose. A genuine close-up "toe splay drill" (the
   original search target, still blocked per §1) is the natural next test.
5. Monocular depth-estimation checkpoints already in this project were identified as the most concrete
   lead for 3D relief but not run this session — real future work, not a result.
6. YouTube's block is IP/session-level per this session's own evidence (§1); it may or may not persist
   at a different time/from a different network — worth a cheap single re-check before assuming it is
   permanent, but not worth further burst-retrying now.

---

## 7. File index

- `scripts/msk/foot_closeup_pose.py` — the extraction + measurement script (re-runnable:
  `~/projects/cad-to-simulation-I/.venv-humancap/bin/python3 scripts/msk/foot_closeup_pose.py`;
  deterministic, reproduced identically across 3 runs this session; self-asserts its own regression
  facts, exit 0/AssertionError).
- `data/msk_pose_closeup/foot_closeup_pose/foot_closeup_pose_summary.json` — full machine-readable
  evidence (every number in this doc traces to a key here).
- `data/msk_pose_closeup/foot_closeup_pose/debug_frames/` — forensic-only debug renders (crops, skin
  masks, contour+hull+defect overlays, hand-landmark skeleton overlays) — used to verify the JSON
  numbers correspond to real anatomy/artifacts, never used as primary evidence themselves.
- `data/youtube/foot/NOTES.txt` — why this directory is empty (§1).
- `docs/MECHANISM_FOOT_FIDELITY_PLAN.md`, `docs/MECHANISM_FOOT_DRILL_CORPUS.md`,
  `docs/MECHANISM_FOOT_MULTISEGMENT.md` — prerequisite reading this doc builds on, not re-litigated.
- `scripts/hand_leg/hand_landmarker.task`, `scripts/hand_leg/frames_extract/chk_VID_20260716_231152638_24.png`,
  `scripts/hand_leg/raw_tracks/VID_20260716_231152638_hands.json` — the hand model + the positive-control
  provenance this doc's resolution-confound control depends on (all pre-existing, reused not
  recomputed).
