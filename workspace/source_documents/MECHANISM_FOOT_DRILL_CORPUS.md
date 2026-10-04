# MECHANISM FOOT DRILL CORPUS — mining the IG movement corpus for foot/ankle content (2026-07-21)

Answers the operator's #1 flagged priority ("all joints incl. feet") and closes an admitted gap:
`docs/MECHANISM_FOOT_MULTISEGMENT.md` (dial-turn 2, the midtarsal joint) was built **geometry-first**
off the LaiArnold LabValidation mocap trials, never driven by this project's own 893-clip IG
athletic corpus — even though the corpus handoff flags foot/ankle/mobility content as high-value.
This doc mines that corpus, cross-references it against `data/msk_pose/`/`data/msk_ik/`, and gives
an honest verdict on whether any clip can actually anchor the placeholder axis that failed to
generalize (held-out DJ2 max-abs RMSE 5.76° vs the 5.0° ceiling, driven by `subtalar_angle_l`).

**Script:** `scripts/msk/foot_drill_corpus_mine.py` (re-runnable, self-checked, exit 0/1).
**Evidence:** `scripts/msk/foot_drill_corpus_mine_evidence.json`,
`scripts/msk/foot_drill_corpus_mine_shortlist.csv`.

Isolation: read-only on `~/ig_downloads/*`, `data/msk_pose/`, `data/msk_ik/`, and this project's own
`scripts/msk/*.py` (imports `pose_extract.py` for its live landmark constants, never re-typed).
This doc + the mining script are the only writes. bodytwin only. No push/commit.

---

## 0. Verifying the premise, not assuming it

The task handed off two claims: "a whole `closeup-low-framemotion-likely-isolation-KEEP` category
= foot/ankle/mobility drills" and "~132 clips match foot/ankle/toe/heel/calf/mobility keywords."
Both were checked against the raw `ig_athletics__semantic.jsonl` (893 rows) rather than taken at
face value:

| Claim | Checked value | Verdict |
|---|---|---|
| "a whole isolation-KEEP category" | the literal `idea_flags` entry `closeup-low-framemotion-likely-isolation-KEEP` hits **3/893** clips, and their captions (generic "joint stability"/"end-range strength"/"hip mobility" copy) are not obviously foot-specific on inspection | **overstated** — it's a sparse, noisy flag, not a category |
| "~132 clips match keywords" | naive substring match of exactly `foot\|feet\|ankle(s)\|toe(s)\|heel(s)\|calf\|calves\|mobility` on the semantic jsonl's own (truncated) `what` field: **134** | **reconciles** — the number itself is real |

The reconciled 134 is still the wrong number to build on: it is dominated by the single generic
word "mobility" (which matches shoulder/hip/spine mobility content identically) and by idiomatic
"foot"/"feet" mentions ("getting back on your feet" = a full-body ground-recovery drill, not a
foot/ankle isolation exercise). A tiered, context-gated reclassification (Sec.1) is needed before
the count is anatomically trustworthy.

A second, independent structured signal exists in `~/ig_downloads/per_clip_metadata.jsonl`:
`movement_tags_secondary` contains `"foot"` on **360/893** clips — even broader than the naive
keyword count. Checked directly (not assumed reliable): **65/360 (18.1%)** of those have **zero**
foot-related keyword anywhere in their own caption text (e.g. "Spinal control isn't just about
posture", "Core strength is not about crunches" — genuinely about the trunk, not the foot). This
auto-tag is *also* not a trustworthy foot-content filter on its own.

## 1. Tiered, context-gated classification (this doc's own)

`scripts/msk/foot_drill_corpus_mine.py` classifies every clip's `what` + full un-truncated caption
into STRICT (anatomically specific: ankle, achilles, subtalar, dorsiflex/plantarflex, calf/heel
raise, plantar fascia, big toe, tibialis, peroneal/us, …) / MODERATE (generic body-part word:
foot, feet, toe, heel, calf, arch, barefoot, …) / WEAK (mobility/balance/proprioception/
footwork/stability *alone*, the words that inflate the naive 134) / NONE.

**A forced adversary, caught, not skipped:** the first pass counted bare "inversion" as STRICT.
It hit `2025-10-17_DP6Op-MkWoK.mp4` — caption about "handstands to L-sits and pikes," i.e. a
gymnastics **inversion** (an upside-down skill), not subtalar inversion. This is exactly the kind
of confound the watertight method requires forcing: inversion/eversion now only counts as STRICT
if a foot/ankle/subtalar/rearfoot word co-occurs in the same caption; otherwise it is logged under
`context_gate_rejected_clips` (auditable), not silently dropped. Result:

| Tier | Count |
|---|---|
| STRICT | 28 |
| MODERATE | 48 |
| WEAK (mobility/balance/etc. alone) | 351 |
| NONE | 466 |
| **STRICT + MODERATE (the trustworthy foot/ankle set)** | **76** |

76, not 132 or 360, is the number this doc builds the shortlist on.

### Movement breakdown (caption-text-derived — explicitly unverified visually beyond the Sec.5 spot-check)

| Bucket | Count | What it means |
|---|---|---|
| `ankle_dorsi_plantarflexion` | 23 | ankle-specific drills, "elite ankles," dorsiflexion/plantarflexion callouts |
| `general_foot_mention` | 23 | mentions foot/feet/heel but names no specific drill/DOF |
| `single_leg_balance` | 18 | single-leg stance/balance/proprioception work |
| `toe_hallux_articulation` | 10 | toe/hallux/plantar-fascia/big-toe content (3 of these also independently mention "arch" — e.g. "the arch lifts, and the toes press into the floor" — co-classified here because toe-pattern matching takes priority; both buckets carry the same CANNOT-anchor verdict, Sec.4) |
| `calf_heel_raise` | 2 | explicit calf/heel raise |
| `subtalar_inversion_eversion` / `midtarsal_forefoot_articulation` / `footwork_agility` (primary) | **0 each** | see Sec.2 — decisive |

## 2. The decisive content-level finding

Direct grep of `subtalar\|midtarsal\|chopart\|tarsometatarsal\|forefoot\|midfoot` across **all 893
full captions**: **0 hits.** Not one clip in this corpus — at the caption/content-description level,
independent of any landmark argument — even claims to target the joint the failed placeholder axis
needs. This settles the task's item (3) before the sensing argument (Sec.4) is even invoked: there
is currently no candidate anchor clip in this corpus for the midtarsal joint specifically, full stop.

## 3. Cross-reference vs `data/msk_pose/` and `data/msk_ik/`

Of the 76 STRICT+MODERATE clips: **43 have pose** (`data/msk_pose/<stem>/`), **42 have an IK
directory with a readable evidence record** (`data/msk_ik/`, joined against the already-built,
already-self-checked `scripts/msk/corpus_ik_aggregate_evidence.json` — reused, not recomputed).
33/76 have neither yet.

## 4. The hard sensing limit — machine-verified, not assumed

Imported live from `scripts/msk/pose_extract.py` (its own `LANDMARK_NAMES`/`MARKER_SUBSET`
constants, never re-typed by hand, so this can't silently drift from the actual pipeline):
MediaPipe BlazePose has **33 total landmarks**; exactly **3 per foot** survive into the `.trc`
driving IK — `Ankle`, `Heel`, `FootIndex`. **No midfoot landmark. No per-toe landmark.** This
matches and independently re-confirms `docs/MECHANISM_FOOT_FIDELITY_PLAN.md` §2's own prior finding
(which made the same point about the real marker-mocap protocol's 3-marker foot).

Geometric consequence (not a heuristic): 3 points can only ever trivially decompose into rigid
segments with **zero residual**, regardless of true anatomy — the same argument the Fidelity Plan
already made for the mocap protocol applies identically here. A real midtarsal break and an
apparent one (from foot orientation change alone) are **architecturally indistinguishable** with
this landmark set. This is the model's null space, not a data-quantity problem more clips can fix.

### CAN-anchor / CANNOT-anchor, by movement bucket (derived from the landmark topology above)

| Bucket | Verdict |
|---|---|
| `ankle_dorsi_plantarflexion` | **YES (qualitative/regression)** — the one DOF with a real 3-marker-per-side analog; already resolved in the existing corpus IK, ROM comparable to the real mocap anchor (Sec.6) |
| `calf_heel_raise` | **PARTIAL** — same ankle DOF shows plantarflexion excursion + rep timing; no tendon/muscle-specific detail |
| `single_leg_balance` | **PARTIAL** — whole-body COM/hip-ankle sway is within monocular-pose reach in principle, not yet computed by this pipeline |
| `footwork_agility` | **PARTIAL** — gross stance/step timing only |
| `subtalar_inversion_eversion` | **NO** — frontal-plane rearfoot motion has no independent landmark; none of the 3 foot points captures medial-lateral (width) deformation |
| `midtarsal_forefoot_articulation` | **NO** — the exact null space above |
| `toe_hallux_articulation` | **NO** — `foot_index` is one point approximating the whole forefoot tip; no per-digit or MTP landmark |
| `arch_intrinsic` | **NO** — arch height/doming needs a dorsal mid-foot point; none exists |
| `general_foot_mention` | **NO SPECIFIC ANCHOR** — caption names no specific drill/DOF to check against |

**Do not overclaim**: of the 76 candidate clips, only the 23 in `ankle_dorsi_plantarflexion` (plus
partial credit for `calf_heel_raise`/`single_leg_balance`/`footwork_agility`, 45 clips combined)
touch a DOF this pipeline can actually inform. The other 33 clips (`toe_hallux_articulation` +
`arch_intrinsic`-co-tagged + `general_foot_mention`) are real, on-topic foot content by caption, but
below the current sensing floor — valuable for *future* capture prioritization (Sec.7), not for
driving today's model.

## 5. Forensic visual spot-check (pixels, not caption trust)

Per the discipline that a caption match is a hypothesis, not a measurement, 1 mid-clip frame was
extracted (`ffmpeg`, saved to scratch, not committed) from the top 4 ranked candidates and viewed
directly:

- `2025-07-31_DMwHVZLsLe9_5` / `_7`, `2026-06-16_DZoF3vaDA9j`: real single-leg-balance / split-stance
  ankle-loading drills, barefoot or shod, on a step/plate — **on-screen text on `DZoF3vaDA9j`
  literally reads "HIGH PERFORMANCE FOOT & ANKLE"**, confirming the caption is not a mismatch.
- `2026-01-11_DTYBc1HEZHG`: on-screen text reads **"✅DORSIFLEXION ✅EXTENSION"** — direct visual
  confirmation of the specific DOF being demonstrated, athlete standing on a round balance disc.

**Earned, not assumed, additional finding**: all 4 are filmed **whole-body wide/medium** shots (to
show stance + balance context), not literal foot close-ups — the foot occupies a small fraction of
the frame even in the best candidates. (A shot-scale check across all STRICT-tier clips found the
mix is 2 wide / 12 medium / 7 close-up, so close framing does exist elsewhere in the corpus, just
not in these top 4.) This means the effective pixel resolution on the foot region is *coarser* than
the landmark-count argument alone suggests — a second, independent reason fine foot detail is out
of reach here, on top of (not instead of) the landmark topology limit in Sec.4.

## 6. Artifact check: is the video-IK's subtalar/mtp output real signal or noise?

Reused directly from `scripts/msk/corpus_ik_aggregate_evidence.json`/`_per_clip.csv` (no IK re-run),
anchored against the REAL marker-mocap ROM table already published in
`docs/MECHANISM_FOOT_FIDELITY_PLAN.md` §2 (a genuinely independent data source — 6 trials, 4 subjects):

| Coordinate | Video-IK mean range (all 228 clips) | Real mocap anchor mean | Ratio / flag |
|---|---|---|---|
| `ankle_angle` | 58.0° | 52.0° | comparable — plausible physiological signal |
| `subtalar_angle` | 48.2° | 23.3° | **2.1x the anchor**; 36.8% of sides pinned near the ±35° mechanical limit |
| `mtp_angle` | 1.0° | 2.2° | near-null both — but **L/R EXACT-TIE rate = 100.0% (228/228 clips)** |

The mtp finding is decisive on its own: **every single one of 228 clips** has `mtp_angle_r` and
`mtp_angle_l` bit-identical (min and max match exactly). No real bilateral human motion produces
that across 228 different people/clips — this proves `mtp_angle` carries **zero** independent
per-side signal in this video pipeline, consistent with (and a stronger, corpus-wide version of)
the pipeline's own existing code-comment finding (`pose_to_opensim_ik.py`: "subtalar_angle jumps up
to 60 deg in a single 33ms frame… axial rotation/inversion needs ≥3 non-collinear markers per
segment to be observable, a geometric fact, not a heuristic").

**Falsifier check, run not skipped**: does the foot-keyword-matched subset (42 clips with IK) look
any less artifact-dominated than the corpus-wide average — i.e. does having foot-relevant *content*
buy better *observability*? Answer: **no** — mtp tie rate 100.0% (42/42), subtalar pinned-side rate
34.5% (vs 36.8% corpus-wide), ankle range 58.2° (vs 58.0°) — statistically indistinguishable from
the whole corpus. This is the expected, honest result: observability is a landmark/geometry
property of the pipeline, not a property of what the clip's caption says it's about.

## 7. Geometry-first confirmation (machine-verified, not assumed)

Grepped `scripts/msk/foot_multisegment.py` + `scripts/msk/foot_detail_proto.py` for any
`ig_downloads`/`msk_pose`/`msk_ik`/`.mp4`/`mediapipe` reference: **0 hits, both files.** Confirms
the task's premise directly: the multi-segment foot model was built entirely off the LaiArnold
LabValidation mocap trials (`walking1.trc`, `DJ2.mot`), never touched by this IG corpus.

## 8. Prioritized shortlist

Score = tier (STRICT +3/MODERATE +1) + pose/IK availability (+2 IK, +1 pose-only) + duration
sanity (+1 if ≥3s, −1 if <1.5s — many IG "clips" are shot-split sub-2s fragments of one post) +
same-post multi-shot compilation bonus (+1) + rare-bucket flag (+1, informational — Sec.4 still
says NO for these) − **3.0 if this specific clip's own `ankle_gate` fails** (reused directly from
`corpus_ik_aggregate_per_clip.csv`'s already-computed per-joint gate, not recomputed — a clip can be
`classification_strict==DEGENERATE` for an unrelated knee/hip reason and still have a perfectly
good ankle number, so the ankle-specific gate, not the whole-clip verdict, is what's scored). 21/76
of the full candidate set fail their own ankle gate (mostly `ankle_angle` pinned at ±50°, i.e. 99+°
of "ROM" that is solver artifact, not motion) — these are excluded from the top ranks below, not
presented as validated.

| Rank | Clip | Tier | Bucket | Pose | IK | Ankle ROM (deg) | Ankle gate | Dur (s) | Score |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `2025-07-31_DMwHVZLsLe9_5` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 25.5 | OK | 5.8 | 7.0 |
| 2 | `2025-07-31_DMwHVZLsLe9_7` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 18.6 | OK | 3.5 | 7.0 |
| 3 | `2026-01-11_DTYBc1HEZHG` | STRICT | toe_hallux_articulation | Y | Y | 18.3 | OK | 9.4 | 7.0 |
| 4 | `2026-06-16_DZoF3vaDA9j` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 47.7 | OK | 3.5 | 7.0 |
| 5 | `2026-06-16_DZoF3vaDA9j_2` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 19.6 | OK | 5.8 | 7.0 |
| 6 | `2026-06-16_DZoF3vaDA9j_5` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 49.4 | OK | 3.3 | 7.0 |
| 7 | `2026-06-16_DZoF3vaDA9j_6` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 22.5 | OK | 3.2 | 7.0 |
| 8 | `2026-06-16_DZoF3vaDA9j_7` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 18.1 | OK | 6.3 | 7.0 |
| 9-11 | `2025-07-31_DMwHVZLsLe9_{2,4,6}` | STRICT | ankle_dorsi_plantarflexion | Y | Y | 24-27 | OK | 2.2-3.0 | 6.0 |
| 12 | `2025-07-31_DMwHVZLsLe9_3` | STRICT | ankle_dorsi_plantarflexion | N | N | — | — | 4.6 | 5.0 |
| 13 | `2025-09-13_DOjHs_PEXnt` | STRICT | toe_hallux_articulation | N | N | — | — | 13.8 | 5.0 |
| 14 | `2026-01-29_DUGxp4NgeIM` | STRICT | toe_hallux_articulation | N | N | — | — | 31.6 | 5.0 |
| 19-20 | `2025-07-24_DMeFzibMRh__{2,4}` | MODERATE | single_leg_balance | Y | Y | 45-51 | OK | 1.6-2.4 | 4.0 |

**What each cluster is actually good for**, being specific rather than hyping a "shortlist":
- **`DMwHVZLsLe9` (7 shots, "7 High Performance Ankle Exercises") and `DZoF3vaDA9j` (7 shots,
  "Exceptional foot and ankle function," on-screen-text-confirmed)**: the best available material
  for **validating/extending the ankle_angle DOF regression** already in the corpus IK — a real
  compilation of distinct ankle drills from one specialty account (`realgame.athletics`), several
  already IK'd and ankle-gate-clean. This is the closest thing to a "held-out ankle validation set"
  this corpus has — genuinely useful, but note it validates the *existing* 3-segment model's ankle
  DOF, not the failed midtarsal placeholder.
- **`DTYBc1HEZHG` ("DORSIFLEXION/EXTENSION" on-screen)**: same ankle DOF, clean IK, good duration.
- **`DOjHs_PEXnt`/`DUGxp4NgeIM`/`DWJd9T5EZVd` (toe/hallux/plantar-fascia captions, no pose/IK yet)**:
  the most anatomically specific toe/arch content in the corpus — worth a pose-extraction pass for
  completeness, but per Sec.4, running them through the existing pipeline will **not** produce a
  midtarsal/toe anchor; flagged instead as the best candidates for a **future finer-grained
  pose model or marker capture**, should one become available (Sec.9).
- **`DMeFzibMRh` (7 shots, "Get-Up Faster")**: `single_leg_balance`/general ground-recovery content;
  PARTIAL value only (whole-body COM proxy), not a foot-joint-specific anchor.

Full 76-row ranked list with all fields: `scripts/msk/foot_drill_corpus_mine_shortlist.csv`.

## 9. Honest verdict: observable now vs. needs better sensing

**Observable with the current MediaPipe-33-landmark + existing IK pipeline, right now:**
- Ankle dorsiflexion/plantarflexion (sagittal-plane shank-vs-foot rotation) — real signal, ROM
  comparable to the independent real-mocap anchor (Sec.6).
- Calf/heel-raise excursion and rep timing (same ankle DOF).
- Gross single-leg balance / weight-bearing timing (whole-body COM proxy; not yet computed by this
  pipeline but architecturally within reach).
- Gross footwork/stance timing.

**NOT observable — below the sensing floor, needs multi-camera, physical markers, or a
purpose-shot future capture (not "more of the same clips"):**
- Subtalar inversion/eversion (frontal-plane rearfoot motion — no landmark captures foot width/
  medial-lateral deformation, at any camera angle this pipeline currently uses).
- Midtarsal/transverse-tarsal (Chopart) articulation — **exactly the placeholder axis's own gap**;
  zero landmarks exist between heel/ankle and foot_index, so real vs. apparent midtarsal motion is
  architecturally indistinguishable, not just under-sampled.
- Tarsometatarsal motion.
- Individual toe/hallux articulation (MTP flexion/extension, hallux abduction) — `foot_index` is one
  point for the whole forefoot, no per-digit landmark.
- Arch height / doming / windlass mechanism — needs a dorsal mid-foot point; none exists.

**What would actually close the gap** (per `docs/MECHANISM_FOOT_MULTISEGMENT.md` §5's own honest-gaps
list, unchanged by this mining pass): either (a) a real cited midtarsal axis from the literature
(Malaquias/Postolka full text, still paywalled), (b) a from-scratch CT-derived axis, or (c) new
capture — physical markers on the midfoot specifically (≥3 non-collinear points, the same geometric
requirement §2 of the Fidelity Plan already derived for the mocap protocol), a multi-camera rig, or
a foot-specific pose model with more than 3 landmarks/foot. **No amount of additional IG clips run
through the existing pipeline substitutes for any of these three** — this is the corpus's null
space, not a data-quantity problem. If the operator does plan a future filming session, the highest-
leverage single addition (per this section's geometry) would be 2-3 markers/keypoints specifically
on the midfoot dorsum, independent of heel/ankle/toe.

## 10. Honest gaps

1. Movement-bucket labels (Sec.1) are caption-text-derived, same `auto-labeled-ig-import-unverified`
   caveat the corpus itself already carries — Sec.5's 4-frame spot-check corroborates the *top*
   candidates but is not a full visual audit of all 76.
2. `single_leg_balance`/`footwork_agility` PARTIAL verdicts (Sec.4) assume a whole-body COM proxy
   pipeline that is architecturally plausible but **not yet built** — flagged as a capability gap,
   not a working feature.
3. The 33/76 clips with neither pose nor IK yet were not run through the pipeline in this session
   (mining/characterization scope only, per the task) — Sec.8 flags which of those are worth a
   future pose-extraction pass.
4. This mining pass covers the 893-clip snapshot as of 2026-07-21; `foot_drill_corpus_mine.py`
   deliberately does not hard-fail if a future corpus refresh adds a clip naming a midtarsal-family
   term (Sec.2) — that would be new information to re-surface, not a regression.

## 11. File index

- `scripts/msk/foot_drill_corpus_mine.py` — the mining script (re-runnable: `python3
  scripts/msk/foot_drill_corpus_mine.py [--top N]`; plain `python3`+numpy, no mediapipe/opensim
  needed since it only reads already-produced JSON/JSONL).
- `scripts/msk/foot_drill_corpus_mine_evidence.json` — full machine-readable evidence (all sections
  above, every number in this doc traces to a key here).
- `scripts/msk/foot_drill_corpus_mine_shortlist.csv` — all 76 STRICT/MODERATE candidates, ranked.
- `docs/MECHANISM_FOOT_MULTISEGMENT.md` / `docs/MECHANISM_FOOT_FIDELITY_PLAN.md` — the model + the
  marker-observability finding this doc extends (not re-litigated here).
- `docs/MECHANISM_CORPUS_IK_AGGREGATE.md` / `scripts/msk/corpus_ik_aggregate.py` — the already-built,
  reused (not recomputed) corpus-wide IK aggregate this doc's Sec.3/6/8 depend on.
