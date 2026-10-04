# MECHANISM CORPUS SLOWMO CASCADE — does the corpus's genuine slow-motion footage resolve the ankle→knee→hip force cascade at its TRUE effective fps? (2026-07-21)

Executes the operator's flagged next step: `docs/MECHANISM_PLYO_JUMP_CASCADE.md` picked its one
clip by flight-run-length alone and never checked `ig_capture.db`'s own `genuine_highspeed`/
`slowmo` shot flags — i.e. it (and every sibling MSK build) implicitly treated every corpus clip
as real-time 30fps. Script: `scripts/msk/corpus_slowmo_cascade.py`
(`.venv-msk/bin/python3 -u scripts/msk/corpus_slowmo_cascade.py`). Every number below is this
script's own printed/JSON output.

## Bottom line

**A minority of this corpus's `genuine_highspeed=1`-flagged jump clips are physics-confirmed as
genuinely slowed; on those, the ankle→knee→hip force cascade MOSTLY still does not resolve, with
one real, non-replicated, method-fragile exception.** Of 12 db-flagged ballistic candidates with
available pose data, an independent gravity-based physics check **confirms** clear slow-down for
only **2/12 (17%)**, is merely **suggestive** for 3/12, and **does not corroborate** the flag at
all for 7/12 — indistinguishable from 4/4 known-real-time negative controls run through the
identical pipeline (s = 0.77–1.27, correctly clustered near 1). Of the 5 clips where cascade
timing could be re-run at the measured effective fps (up to 251fps effective, vs. the original
30fps-assumed floor): **4/5 still show the exact-same-frame tie** the original analysis found —
proving directly that finer effective resolution does not automatically break a tie (a tie at
identical frame INDEX is invariant to any time-rescaling). **Exactly 1/5** — our best-quality
confirmed-slowed clip — shows a real, moderately robust (SavGol- and takeoff-instant-stable)
**~20ms-real** stagger (ankle/knee≈tied, then hip, then back), but it does not survive an
independent cross-check (joint angular-velocity peak timing gives a different order on the same
clip/window) and is a single, non-replicated instance. **Net: mostly still tied, one intriguing
exception — reported as inconclusive/leaning-tied, not forced into either a clean resolve or a
blanket negative.**

**A major, load-bearing discovery made getting here**: the existing OpenSim/`.trc` pipeline's own
"pelvis height"/"freefall" signal — used by `force_scenes_batch.py` and
`plyo_jump_cascade.py` for flight-phase detection — **cannot be reused to verify slow-motion**,
because it is built entirely from MediaPipe `world_landmarks`, which are hip-centered PER FRAME
by construction (the pipeline's own pre-registered gate predicts `pelvis_ty` range ≈0, verified
live below at 0.056m on the original doc's own clip). This build instead measures directly on the
RAW, previously-unused `image_landmarks` channel already sitting in `pose_raw.json`.

## 1. Why the existing freefall/pelvis-height signal is the wrong tool for this job (verified, not assumed)

`pose_extract.py`'s own docstring and `pose_to_opensim_ik.py`'s own pre-registered gate
(`pelvis_translation_void_floor_m: 0.05`) already establish, from an earlier session's own
measurement, that MediaPipe `world_landmarks` — 100% of what feeds the `.trc`/IK/
`force_scenes_batch.py` pipeline — have mid-hip pinned at `|pos| < 0.002m` from the coordinate
origin **in every frame**. Whole-body translation (jump height) is removed from that channel by
construction. Re-verified live here, directly on the ORIGINAL cascade doc's own chosen clip
(`2026-06-30_DaMKH8aDBxu_6`), not merely quoted from the docstring:

| quantity | value | pre-registered expectation | verdict |
|---|---:|---|---|
| `pelvis_ty` range | **0.0557 m** | ≤0.05m void floor ("SMALL, confirms hip-centered/no-global-translation") | borderline-small; **far too small to be a real jump's flight-height excursion** |

**Consequence, derived (not guessed):** `com_raw["calcn_r"]` (the anchor-foot trajectory
`force_scenes_batch.py`'s freefall diagnostic thresholds at `<-0.6g`) = a near-constant
`pelvis_ty` **+** a kinematic-chain term driven entirely by hip/knee/ankle joint angles. During
real free flight this channel cannot see the true whole-body `-g` acceleration; during a fast
voluntary leg articulation (tuck/knee-drive) it readily produces `-0.6g`-scale apparent
acceleration from articulation alone. The existing "flight phase" detector is therefore likely
measuring **fast leg articulation**, not verified airborne time — and cannot be used, without
begging the question, as this build's own slow-motion-verification flight detector. (This also
offers a candidate explanation, not re-litigated further here, for why the original doc found
this corpus's "jump" content 10/11–15/21 tuck-dominated: a detector keyed to hip-relative
acceleration is structurally biased toward finding fast-articulation tuck events over
legs-extended flight.)

This build instead measures on `pose_raw.json`'s **`image_landmarks`** (normalized image-plane
x/y/z/visibility) — the other channel MediaPipe already emits, unused by the existing `.trc`/IK
pipeline, which *does* preserve whole-body position within the camera frame.

## 2. Method — gravity-parabola fit, anchored externally, forced through an adversary

**Geometry**: during true free flight, gravity is the only force, so any body-fixed point's
vertical position is an exact parabola in REAL time: `y(t_real) = y0 + v0·t_real − 0.5·g·t_real²`.
If clip-clock time is a stretched version of real time by slow-mo factor `s` (`t_real = t_clip/s`,
`s>1` = slowed), fitting the SAME parabola in clip-clock time gives curvature
`C_fit = −0.5·g/s²`, so **`s = sqrt(g / (2·|C_fit|))`** falls out directly — anchored by
`g = 9.80665 m/s²` (a physical constant) and a standard external anthropometric ratio
(Drillis & Contini / Winter tables: leg length ≈0.491×stature, trunk ≈0.288×stature), **never**
by the db's own `genuine_highspeed` flag. `effective_fps = container_fps × s`.

**Forced adversary** (an unconstrained single-landmark scan was tried FIRST and happily found
high-R² "parabolas" that were squats/dips/articulation — caught, not glossed over): the final
method REQUIRES hip, shoulder, AND nose to independently fit **positive** curvature of similar
magnitude in the same window (a genuinely falling body moves all its points together; a squat
does not move hip and ankle the same way; pure limb articulation does not move the trunk at all).
Anthropometric scale uses a **foreshortening-resistant estimator** (90th-percentile projected
segment length over a wide span — any camera angle/joint bend can only shrink a segment's
apparent length relative to its true fronto-parallel value, never lengthen it); this fixed an
early version of the calibration that was unreliable on 12/14 test clips, improved to 10/14 (and
similar in the final corpus run).

**Validation (forced before trusting any positive result)**: 4 clips from the same account/
category, **not** flagged `slowmo` or `genuine_highspeed` in the db, run through the identical
pipeline:

| stem | s (stature-swept range) | tier |
|---|---:|---|
| `2025-07-24_DMeFzibMRh__2` | 0.91–0.99 | NOT_CONFIRMED |
| `2025-07-24_DMeFzibMRh__4` | 1.22–1.33 | NOT_CONFIRMED |
| `2025-07-24_DMeFzibMRh__6` | 1.13–1.23 | NOT_CONFIRMED |
| `2025-07-24_DMeFzibMRh__7` | 0.74–0.81 | NOT_CONFIRMED |

**All 4 cluster near s≈1 (range 0.77–1.27) and all correctly clear as NOT_CONFIRMED** — the
method has the specificity to not spuriously report large slowdown on ordinary real-time footage.

**Pre-registered thresholds** (set before the final corpus-wide run, applied mechanically,
including to the negative controls): `s_lo` (low end of the stature sweep) `≥1.5` to count as
elevated at all; **CONFIRMED** additionally requires reliable leg/trunk calibration (ratio
∈[0.7,1.4]) and detection in the primary (non-relaxed) shot-window scope; anything elevated but
missing one of those is **SUGGESTIVE**; `s_lo<1.5` is **NOT_CONFIRMED**.

## 3. Candidate discovery (live `ig_capture.db` query)

`ig/`-prefixed (IG-sourced) shots with `genuine_highspeed=1 OR slowmo=1`, joined to
`ig_athletics__semantic.jsonl` for `movement_group`, restricted to ballistic groups
(jump/plyo/sprint/run):

| | count |
|---|---:|
| Ballistic shots flagged `genuine_highspeed=1` OR `slowmo=1` | 41 |
| …of those, strictly `genuine_highspeed=1` | 32 |
| …of those 32, already have `pose_raw.json` on disk (reused, not re-extracted) | 12 |
| Fresh MediaPipe extraction run this build (`.venv-humancap`, read-only on -I) | 1 (`2026-06-04_DZJMf91s7uA`) |

## 4. Slow-mo verification results (12 db candidates + 4 negative controls + 1 context clip)

| stem | shot window (s) | R²(hip/sh/nose) | calib | s range | eff. fps | **tier** |
|---|---|---|---|---:|---:|---|
| `2026-05-30_DY-yVTeuXyR` | [0.73,1.90] | 1.00/1.00/1.00 | OK (1.02) | 5.42–5.92 | **170** | **CONFIRMED** |
| `2026-04-10_DW89psJkf5d` | [7.53,7.77] | 0.94/1.00/0.98 | OK (0.86) | 1.61–1.76 | **50** | **CONFIRMED** |
| `2026-07-04_DaYjAaquUq0` | [1.13,2.10] | 0.93/0.96/0.92 | flag (1.44) | 7.04–7.69 | 220 | SUGGESTIVE |
| `2026-01-21_DTxsrw9EfAN` | [38.7,39.5] | 0.82/0.92/0.84 | flag (1.51) | 8.03–8.77 | 251 | SUGGESTIVE |
| `2026-06-04_DZJMf91s7uA` | [28.3,28.5] | 0.98/0.99/1.00 | flag (1.68) | 1.86–2.03 | 58 | SUGGESTIVE |
| `2025-12-30_DS5ULQvEUXu` | [8.10,8.33] | 0.95/0.98/0.87 | flag (0.55) | 1.09–1.19 | 34 | NOT_CONFIRMED |
| `2026-05-09_DYH8d74RhGi` | [12.5,12.7] | 0.99/1.00/0.99 | OK (1.40) | 1.07–1.17 | 34 | NOT_CONFIRMED |
| `2026-07-06_DaeDohKuJK2` | [14.5,14.6] | 1.00/1.00/1.00 | flag (1.61) | 1.34–1.47 | 42 | NOT_CONFIRMED |
| `2026-07-13_DavKz_lkXU5_6` | [2.20,2.57] | 0.93/0.99/1.00 | flag (0.38) | 0.89–0.97 | 28 | NOT_CONFIRMED |
| `2026-05-01_DXxohp9Mr65` | [21.1,21.3] | 1.00/1.00/1.00 | OK (1.29) | 0.85–0.93 | 27 | NOT_CONFIRMED |
| `2026-06-26_DaDba5vRhzO` | [37.9,38.1] | 0.85/0.96/0.96 | OK (1.35) | 0.60–0.65 | 19 | NOT_CONFIRMED |
| `2026-06-18_DZusbAbFwyV` | [49.5,49.8] | 0.84/0.86/0.84 | OK (1.29) | 0.52–0.57 | 16 | NOT_CONFIRMED |
| `2026-06-30_DaMKH8aDBxu_6`\* | [2.73,3.70] | 0.80/0.81/0.74 | flag (0.59) | 3.10–3.39 | 97 | SUGGESTIVE\* |
| 4× negative controls | — | 0.99–1.00 (all) | OK (all) | 0.74–1.33 | 23–38 | **all NOT_CONFIRMED (correct)** |

\* `2026-06-30_DaMKH8aDBxu_6` (the ORIGINAL cascade doc's own clip) is **not** flagged
`slowmo`/`genuine_highspeed` in the db at all (verified live: `slowmo=0`, `genuine_highspeed`
NULL for all 7 shots of this video) — carried here purely for continuity, not as a "found"
candidate. Its own SUGGESTIVE reading is low-confidence (weak fit, found only via the relaxed
scan tier, unreliable calibration) — reported as an honest, low-confidence secondary observation
that the original doc's own 33ms/"resolution floor" framing might deserve a second look, **not**
a claim that it was mislabeled.

**Tier tally**: of the 12 actual db-flagged candidates tested — **CONFIRMED: 2 (17%), SUGGESTIVE:
3 (25%), NOT_CONFIRMED: 7 (58%)**. The `genuine_highspeed` flag — "the optical pipeline's own
guess" per the task — **holds up under an independent physics check for a minority, not the
majority, of tested candidates on this corpus.**

## 5. Cascade-timing re-run at measured effective fps

Reused, unmodified: `plyo_jump_cascade.compute_core`/`propulsion_window`/
`cascade_order_in_window`, `force_scenes_batch`'s Newton-reaction-force machinery,
`validate_joint_force`'s `parse_mot`/`savgol_smooth_and_derivs`. New: the propulsion window's
"takeoff" instant is THIS build's own physics-verified flight-onset time (nearest-time match into
the `.mot`'s own frame index), not the confound-contaminated freefall detector (§1); frame-gaps
convert to REAL ms via measured `dt_real = dt_container/s`. Run on all 5 CONFIRMED/SUGGESTIVE
clips that have existing IK data (`2026-06-04_DZJMf91s7uA` has none yet — skipped, not
extrapolated).

| clip | tier | eff. fps | FORCE-peak order | gaps (real ms) | magnitude vs. 400%BW ceiling | ANGVEL-peak order (cross-check) |
|---|---|---:|---|---|---|---|
| `DW89psJkf5d` | CONFIRMED | 50 | ankle<knee<hip<back | **[19.6, 19.6, 20.2]** | **all under** (112–185%BW) | knee<ankle<hip (**disagrees**) |
| `DY-yVTeuXyR` | CONFIRMED | 170 | ankle=knee=hip=back | **[0, 0, 0] (TIED)** | over (2096–3576%BW) | window too short (n<23) |
| `DaYjAaquUq0` | SUGGESTIVE | 220 | ankle=knee=hip=back | **[0, 0, 0] (TIED)** | over (7479–10482%BW) | knee<hip<ankle (n/a, tied case) |
| `DTxsrw9EfAN` | SUGGESTIVE | 251 | hip<ankle=knee=back | [111.4, 0, 0] | over (2680–4013%BW) | hip<ankle<knee (partial agreement: hip first) |
| `DaMKH8aDBxu_6`\* | SUGGESTIVE(weak) | 97 | ankle=knee=hip=back | **[0, 0, 0] (TIED)** | over (603–840%BW) | knee=hip<ankle (n/a, tied case) |

**Geometric point proven directly, not just asserted**: rescaling `dt` cannot by itself break a
tie already found at identical frame INDEX (the `argmax` of a curve is invariant to a positive
rescaling of its time axis). `DY-yVTeuXyR`, `DaYjAaquUq0`, and `DaMKH8aDBxu_6` remain **exactly**
tied (0.0/0.0/0.0 ms) at effective sampling rates of 97–220fps — direct, measured proof that
"the original 30fps analysis just didn't have enough resolution" is **not** a universally correct
diagnosis: on 3 of these 5 instances, correctly interpreting the SAME already-computed peak
frames as far-finer-grained samples changes nothing, because the underlying peaks genuinely
coincide at the same discrete sample.

### 5.1 The one non-tied result — forced through robustness adversaries before trusting it

`DW89psJkf5d` (our best-quality CONFIRMED clip: reliable calibration, primary-tier detection,
and the only clip whose s²-corrected force magnitudes stay under the plausibility ceiling) shows
a genuinely non-tied order. Forced adversary checks, same discipline the original doc used:

- **SavGol-window sweep** (5 windows): **4/5 agree** on `ankle<knee<hip<back` with the same
  ~19.6–20.2ms real gaps between knee→hip and hip→back (ankle/knee sometimes tie within this).
  Window=31 alone reorders to `back<ankle<knee<hip` — a single non-adjacent-reproducing outlier,
  same "does not survive the sweep" signature the original doc itself used to dismiss its own
  one outlier window.
- **Takeoff-instant sensitivity** (±3 frames around the detected flight onset): **identical**
  order and gaps for shifts −3 through +1 (6/7 shifts); at +2/+3 the ankle<knee<hip leg is
  unchanged but the hip→back gap grows (`back`'s own parallel-branch ambiguity, inherited
  unresolved from the original doc's §5.4).

**But**: the independent joint-angular-velocity-peak method (same window, same clip) gives
`knee<ankle<hip` — a **different order** than the force-peak method's `ankle<knee<hip<back`.
No clip in this sample shows an order that is both non-tied AND cross-method-consistent.

## 6. Does the cascade resolve? — symmetric verdict

**No, not robustly, across this corpus's verified-slowed clips.** In 4 of 5 testable instances —
including our *second*-best slow-mo confirmation (`DY-yVTeuXyR`, the highest-R²/most-reliable-
calibration clip of all 12) — the exact-frame tie the original 30fps-assumed analysis found
**persists** even at 97–220fps effective resolution, directly demonstrating that (for those
instances) more resolution did not break the tie: this supports, rather than overturns, the
original doc's own geometric argument that a shared whole-body driver can genuinely dominate,
not merely appear to because of a frame-rate floor. Exactly 1 of 5 (`DW89psJkf5d`) shows a real,
~20ms, moderately-robust-to-SavGol-and-takeoff-choice stagger — but it is a **single, non-
replicated instance**, it does **not** survive an independent (joint angular-velocity)
operational definition of "peak timing," and it (like every other candidate in this sample)
carries a disclosed `ankle_range_gate`-FAIL / joint-pinning pose-quality caveat (§7). This does
**not** clear the diverse-instance-space bar needed to accept a general "the cascade resolves"
claim. Reported as **mostly-still-tied, with one intriguing, non-replicated, method-fragile
exception** — not forced into either a clean positive or a lazy blanket negative.

## 7. Honest limits

1. **Pose/IK quality caveats are pervasive and correlate inversely with slow-mo confirmation
   strength in this sample.** Every one of the 8 `CLEAN_BORDERLINE` gh=1 candidates checked
   against the corpus catalog shows `classification_strict=DEGENERATE` with an
   `ankle_range_gate FAIL`, and most also show hip/knee pinned at the model's mechanical limit.
   `DW89psJkf5d` (our one non-tied result) is not exempt from this — its ankle_range_gate also
   fails and hip is pinned. `DaYjAaquUq0` is outright `DEGENERATE` (`coverage_frac=0.23`). Only
   `DaMKH8aDBxu_6` (the original doc's clip) is `CLEAN_NO_CAVEAT`, and it has the *weakest*
   slow-mo confirmation of the five. No clip in this sample is simultaneously cleanly-classified
   AND clearly-and-reliably confirmed slowed — an intersection this corpus's available data does
   not (yet) populate.
2. **Anthropometric scale is an assumed external constant** (stature 1.55–1.85m sweep, standard
   Drillis&Contini/Winter leg/trunk fractions), not this athlete's real height — the slow-mo
   factor's dependence on this is a favorable, damped `sqrt(1/scale)`, confirmed numerically (the
   full 1.55–1.85m sweep moves `s` by ≤10% in every case above), but a true anthropometry error
   beyond that range would shift the exact `s` value, not just its precision.
3. **Single-camera monocular foreshortening**: the leg/trunk calibration-consistency check
   (§4 "calib" column) flags roughly half the tested clips (7/13) as unreliable — those `s`
   values are reported but should be weighted less; this is disclosed per-clip, not averaged
   away.
4. **`back`'s parallel-branch ambiguity is inherited, unresolved** (no GRF term; its exact
   peak-time coincidence with the GRF-bearing cuts remains open, same as both prior builds).
5. **The "takeoff instant" is this build's own detected common-mode-ballistic-window onset**,
   not an independently-verified biomechanical toe-off frame — a reasonable, geometrically-
   motivated proxy (§2), cross-checked for stability (§5.1), but not ground-truthed against a
   force plate or a second camera.
6. **`s²`-corrected force magnitudes are frequently wildly implausible** (up to ~10,000%BW) on
   the higher-`s`, lower-pose-quality clips — read as a joint finding: either those specific
   `s` estimates are inflated (plausible given their flagged-unreliable calibration) or the
   underlying pose/IK data is not fit to support kinematics at their implied sub-5ms effective
   resolution, or both. Not resolved further here; reported as a disclosed complication, not
   hidden by only reporting the timing (order) numbers.
7. **GRF is a kinematically-estimated, Tier-1 symmetric 50/50 L/R split** and subject2's
   anthropometry stands in for every athlete — both inherited, unchanged limits from every prior
   MSK build in this repo.
8. **n=12 db-flagged candidates tested** (of 32 in the corpus flagged `genuine_highspeed=1`
   overall) — a live, re-runnable sample restricted to clips that already had pose data on disk
   plus one fresh extraction, not an exhaustive sweep of all 32.

## 8. Reproducing this result

```
.venv-msk/bin/python3 -u scripts/msk/corpus_slowmo_cascade.py
```

Runs the live db candidate discovery (~1s), the 17-clip slow-mo physics measurement (~15s), and
the cascade-timing re-run including the SavGol-window and takeoff-instant sensitivity sweeps on
the two CONFIRMED clips (~90s, OpenSim model load dominates). Exit code reflects the
method-trustworthiness gate (negative controls correctly NOT_CONFIRMED + at least one cascade
re-run completed) — **not** whether a cascade was "found" (this repo's own
honest-negative-is-PASS convention).

## Files

- `scripts/msk/corpus_slowmo_cascade.py` — the full pipeline (self-contained, re-runnable;
  imports `plyo_jump_cascade.py`'s `compute_core`/`propulsion_window`/`cascade_order_in_window`,
  `force_scenes_batch.py`'s Newton-reaction primitives + constants, and
  `validate_joint_force.py`'s `parse_mot`/`savgol_smooth_and_derivs` — none re-implemented).
- `data/msk_smoketest/corpus_slowmo_cascade/corpus_slowmo_cascade_results.json` — every number
  above, machine-written (all 17 slow-mo measurements incl. tiers, all 5 cascade-timing results
  incl. both pre-margin conventions, the SavGol sweep, and the takeoff-instant sweep).
- Fresh MediaPipe extraction this build added (`.venv-humancap`, read-only on
  `cad-to-simulation-I`): `data/msk_pose/2026-06-04_DZJMf91s7uA/` (pose only; not yet IK'd, so not
  part of the cascade-timing table).
- Reused, unmodified: `scripts/msk/plyo_jump_cascade.py`, `scripts/msk/force_scenes_batch.py`,
  `scripts/msk/validate_joint_force.py`, `data/msk_ik/*/smoothed/*.mot`,
  `data/msk_pose/*/*.pose_raw.json` (the `image_landmarks` channel specifically — previously
  unused by any script in this repo).
