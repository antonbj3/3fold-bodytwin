# MECHANISM CORPUS IK AGGREGATE — per-clip scene catalog over the corpus IK batch

Written 2026-07-21. Status: **228/228 `data/msk_ik/` clips aggregated, machine-classified
CLEAN+PHYSIOLOGICAL vs DEGENERATE against thresholds inherited verbatim from the pipeline's
own pre-registered gates (no new numbers invented, none threshold-shopped after seeing the
corpus).** Every number below comes from running `scripts/msk/corpus_ik_aggregate.py` (prints
the full report + writes a per-clip CSV + a structured evidence JSON); nothing here is
estimated or eyeballed from a plot. Builds on `docs/MECHANISM_CORPUS_IK.md` (the pipeline this
aggregates) without re-litigating it.

Isolation: read-only on `data/msk_ik/`, `data/msk_pose/`, `~/ig_downloads/`; this script and
doc are its only writes; no push/commit performed.

## 0. The one mistake this aggregate is built not to repeat

A prior quick aggregate pass over this same corpus mistakenly reported each coordinate's
**fixed model mechanical limit** (`knee=[0,140]`, `hip_flexion=[-30,120]`, `ankle=[-50,50]` —
literally the same 2 numbers for every clip in the corpus, from
`pose_to_opensim_ik.MODEL_COORD_LIMITS_DEG`) as if it were that clip's **observed** range of
motion. This aggregate reads only `angle_stats[coord]["min"/"max"/"range"]` — the per-clip,
per-frame statistic `analyze_results()` computes via `np.min`/`np.max` over that clip's own
solved `.mot` column — and carries the model limit back out under an unambiguously separate,
labeled key (`model_mech_limit_deg`) purely as a diagnostic, never conflated with the observed
value. Three machine self-checks prove this, run and printed **before** any summary claim:

| self-check | result |
|---|---|
| **[1] Pilot-clip regression**: this aggregate's own extraction of `2025-09-26_DPEOqLVkav7` must reproduce `docs/MECHANISM_CORPUS_IK.md` §4's already-published, hand-verified numbers (knee 35.6°/91.2°, hip_flexion 32.6°/148.2°, ankle 47.9°/84.6°, marker RMS 0.0857 m) — an external anchor written before this aggregate existed. | **PASS**, exact match to the doc's stated precision (0/7 mismatches) |
| **[2] Not-a-constant-column**: if the bug above recurred, every clip's reported range would collapse to one identical number per coordinate (stdev **exactly** 0.0). Measured corpus-wide stdev instead: | **PASS** — knee σ=42.03°, hip_flexion σ=47.90°, ankle σ=30.89° (all ≫ the 5° pass floor; a constant column gives exactly 0.0) |
| **[3] det_frac == coverage_frac(raw)**: the upstream pose-extraction detection fraction and this IK stage's own independently-recomputed raw-frame coverage must agree (both count non-NaN `.trc` rows off the same file, computed by two different scripts). | **PASS** — 228/228 checked, 0 mismatches |

All three PASS. Script exits 1 if any self-check fails, so a future regression cannot pass
silently.

## 1. What was extracted, per clip (the 5 requested fields)

For every directory under `data/msk_ik/` with a `<stem>_ik_report.json` (228 directories present,
228 with a report — 0 missing; the 4 non-directory files at that path's top level,
`_aggregate_summary.json`/`_ik_batch_log.jsonl`/`_overall_pass_manifest.json`/
`_overall_pass_stems.txt`, are the batch-runner's own bookkeeping — the first is a **stale**
mid-batch checkpoint from `n_processed=19/228`, superseded by the now-complete 228; the other
three are the upstream POSE-stage's 264-clip gate manifest — none are per-clip IK output and
none were iterated as if they were):

1. **Observed per-coordinate ROM** — `min`/`max`/`range`/`range_p5_p95` for all 22
   `COORDS_OF_INTEREST` (imported directly from `pose_to_opensim_ik`, not retyped), focused for
   classification on the 3 gated pairs: `knee_angle_{r,l}`, `hip_flexion_{r,l}`,
   `ankle_angle_{r,l}`. Read from the `chosen` variant (`"smoothed"` in 227/228 reports, the
   pipeline's own recommended final configuration; 1 report additionally carries a
   `axis_remap_check` block from the `--check-axis-remap` pilot run — harmless, `chosen` is
   still `"smoothed"` there too).
2. **IK RMS + max marker error** — `marker_error_rms_mean_m`, `marker_error_rms_max_m`,
   `marker_error_max_mean_m`, `marker_error_max_max_m`, plus the pipeline's own
   `marker_verdict` (PASS/MARGINAL/FAIL against the pre-registered 0.15 m / 0.30 m ceilings).
3. **det_frac / pose quality** — `det_frac` from `data/msk_pose/<stem>/<stem>.summary.json`
   (the upstream, TRUE detection fraction) as the primary figure, **kept distinct** from
   `coverage_frac` (this IK report's own frame-solve fraction). Measured, not assumed: these
   two normally agree (verified, self-check [3] above), but for clips with any dropped
   all-NaN frames, the **smoothed** variant's `coverage_frac` can run **higher** than the true
   `det_frac` by up to +0.25 (39/228 clips differ by >0.01) — because centered-window
   smoothing can partially fill a previously-all-NaN frame from a detected neighbor before
   that frame reaches `strip_allnan_frames_trc()`. Both numbers are reported per clip, never
   merged into one.
4. **Frame count + duration** — `n_frames_original` (pose-extraction total), `n_frames_solved`
   (chosen variant, post-NaN-strip), `fps_used`, and `duration_s` (preferring the semantic
   jsonl's own independently-recorded `duration_s`, cross-checked against
   `n_frames_original/fps_used`: **0 mismatches >0.5 s across all 228 clips** with both
   available — total corpus duration agrees to within 12 s either way: 3977.0 s via the
   semantic field vs 3965.2 s via frames/fps).
5. **Movement group** — exact video-basename join against
   `~/ig_downloads/ig_athletics__semantic.jsonl` (893 unique basenames, **0 duplicates**,
   checked not assumed). **228/228 clips matched, 0 unmatched.**

## 2. Classification: CLEAN+PHYSIOLOGICAL vs DEGENERATE

**Criteria are inherited verbatim from `pose_to_opensim_ik.PRE_REGISTERED`/`verify()`** — no
new thresholds were invented and none were tuned after looking at the corpus-wide pass rate.

**DEGENERATE** if **any** of:
- `marker_verdict == FAIL` (mean marker RMS > 0.30 m)
- `non_degenerate_all == False` (≥1 of the 6 sagittal knee/hip/ankle DOFs has range ≤1° — a
  frozen/near-frozen solve)
- `coverage_frac (chosen variant) < 0.5` (majority of frames never solved)
- `knee` range gate FAIL (observed range outside [20°,165°]) **or** `hip_flexion` range gate
  FAIL (outside [15°,150°])

**`ankle`'s range gate is deliberately excluded** from this hard rule — folded instead into a
soft caveat on an otherwise-CLEAN clip. This is not a lenient carve-out invented to inflate the
CLEAN count: `docs/MECHANISM_CORPUS_IK.md` §3.2/3.5/9 already forced an OODA loop on this
*before* this aggregate was written and found ankle/subtalar is the single weakest-observed DOF
in this 16-marker set for a geometric reason (a 1–2-marker-per-segment reduced set cannot
observe axial/inversion rotation at all), predicting in advance that this specific gate would
fail often for reasons other than the clip's true physiology. **The STRICT variant (ankle
folded into the hard gate) is computed and reported in full below regardless**, so this choice
is never hidden.

Else **CLEAN+PHYSIOLOGICAL**, sub-tagged **no-caveat** (marker PASS, ankle gate PASS, not
pinned-at-mechanical-limit, no L/R visibility asymmetry) vs **borderline** (≥1 soft caveat,
specific caveat always reported, never silently absorbed).

### Results (228 clips total)

| | LENIENT (knee+hip_flexion hard; ankle=soft) | STRICT (knee+hip_flexion+ankle all hard) |
|---|---:|---:|
| **CLEAN+PHYSIOLOGICAL** | **181** (48 no-caveat + 133 borderline) | **67** |
| **DEGENERATE** | **47** | **161** |

Headline number to actually use downstream is the **lenient 181/47** split — but report both,
always, since 133/181 (73%) of "CLEAN" clips carry ≥1 caveat and the strict/lenient gap is
large (181 vs 67) entirely because of the one already-diagnosed weak ankle DOF. **Only 48/228
(21%) are clean with zero caveats of any kind.**

Degenerate reason breakdown (a clip can trip more than one, so these sum to more than 47):

| reason | count |
|---|---:|
| `knee_range_gate=FAIL` | 34 |
| `hip_flexion_range_gate=FAIL` | 26 |
| `non_degenerate_all=False` (a frozen DOF, range≤1°) | 11 |
| `coverage_frac<0.5` (majority of frames unsolved) | 1 |
| `marker_verdict=FAIL` | 0 |

Caveat patterns within the 181 lenient-CLEAN clips (top rows):

| pattern | count |
|---|---:|
| none (fully clean) | 48 |
| ankle_fail + lr_visibility_asymmetry | 21 |
| ankle_fail + knee_near_full_span + pinned | 21 |
| ankle_fail only | 20 |
| ankle_fail + knee_near_full_span + lr_asymmetry + pinned | 19 |
| lr_visibility_asymmetry only | 18 |
| ankle_fail + pinned | 17 |
| ankle_fail + lr_asymmetry + pinned | 12 |
| (remaining 4 rarer combinations, incl. `marker_verdict=MARGINAL`) | 5 |

### Movement-group × classification (lenient)

| movement group | CLEAN | DEGEN | total | degen% |
|---|---:|---:|---:|---:|
| `plyometrics_jumps/technique_tutorial` | 113 | 33 | 146 | 22.6% |
| `plyometrics_jumps/plyometric_drill` | 32 | 10 | 42 | 23.8% |
| `weighted_lifts/squat` | 25 | 3 | 28 | 10.7% |
| `plyometrics_jumps/plyometrics_jumps_misc` | 11 | 1 | 12 | 8.3% |

Descriptive, not causally over-claimed: squat/misc clips run a much lower degenerate rate than
technique_tutorial/drill clips — consistent with squats forcing a large, reliable dynamic
knee/hip excursion every rep, while technique-tutorial clips more often include static holds or
close-up partial-body demonstration segments with genuinely small real ROM (see the
`DMwHVZLsLe9_4` example below — this is a correct catch of real low-motion content, not a
pipeline artifact).

## 3. Observed-ROM sanity check — concrete clips, proving these are per-frame values

**(a) The pilot regression itself** (`2025-09-26_DPEOqLVkav7`, `weighted_lifts/squat`) — this
aggregate's own numbers, reproducing the doc's already-published table exactly:

| joint | observed min | observed max | observed range | model's fixed mech. limit (NOT the same field) |
|---|---:|---:|---:|---|
| knee_angle_r | 22.37° | 57.94° | **35.56°** | [0°, 140°] |
| knee_angle_l | 41.08° | 132.30° | **91.22°** | [0°, 140°] |
| hip_flexion_r | -20.93° | 11.67° | **32.60°** | [-30°, 120°] |
| hip_flexion_l | -29.89° | 118.35° | **148.23°** | [-30°, 120°] (near_full_mechanical_span=True) |
| ankle_angle_r | 1.31° | 49.21° | **47.90°** | [-50°, 50°] |
| ankle_angle_l | -34.85° | 49.78° | **84.63°** | [-50°, 50°] |

marker RMS mean 0.0857 m; range gates knee=PASS, hip_flexion=PASS, ankle=FAIL — all identical
to `docs/MECHANISM_CORPUS_IK.md` §4's table.

**(b) A clean, no-caveat clip** (`2025-07-24_DMeFzibMRh__2`, `plyometrics_jumps/technique_tutorial`,
1.55 s / 46 frames, det_frac=1.0): knee 45.2°/85.7°, hip_flexion 68.0°/136.6°, ankle
41.0°/49.6° — all 3 range gates PASS, marker RMS 0.0821 m, no pinned/asymmetry/marginal flags.

**(c) A genuinely degenerate clip — small true ROM, not a pipeline bug**
(`2025-07-31_DMwHVZLsLe9_4`, `plyometrics_jumps/plyometric_drill`): knee 18.7°/13.6° (both
below the 20° floor), hip_flexion 14.3°/18.4° (right side below the 15° floor), ankle
24.7°/23.5° (within bounds). Marker fit is fine (RMS 0.0780 m, PASS) — the clip simply doesn't
contain enough real joint excursion to be a usable force-transmission scene, correctly caught
by the range gates rather than the marker-error gate.

**(d) A genuinely frozen-DOF clip** (`2025-12-10_DSFMpkNkYru`): `knee_angle_r` range **0.11°**,
`ankle_angle_l` range **0.09°** — solved but essentially motionless on those two DOFs, despite
marker RMS still reading 0.1395 m PASS. This is the concrete proof that marker-fit-PASS alone
is **not** sufficient — `non_degenerate_all` is a necessary, independent gate, not redundant
with the marker-error check.

**(e) A pinned-at-mechanical-limit clip, all three joint pairs simultaneously**
(`2025-10-08_DPh0fvwDLex`): knee range 139.91°/139.85° (this coordinate's own **mechanical**
ceiling is 140°), hip_flexion range 149.98°/149.94° (mechanical ceiling 150°), ankle range
99.86°/99.85° (mechanical ceiling 100°) — **every one of the 6 sagittal DOFs sits within 0.2°
of its own model's mechanical hard stop, on both sides, at once.** Two different ceilings are
in play here, and it matters which: knee's *pass-gate* ceiling is a separate, looser 165°, so
139.9° clears it with a comfortable 25° of margin — not a near-miss on the gate itself, only on
the model's mechanical stop. hip_flexion's pass-gate ceiling, by contrast, is exactly 150°
(§9 of `docs/MECHANISM_CORPUS_IK.md` already found this is a coincidence, not a deliberate
design choice — the pass-gate ceiling happens to equal the coordinate's own full mechanical
span), so 149.98°/149.94° really is passing *at* the edge of both numbers simultaneously.
ankle FAILs its (60°) pass-gate outright. Classified
`CLEAN_BORDERLINE` here, correctly caveated (`pinned_at_mechanical_limit`,
`knee_near_full_mechanical_span`, `ankle_range_gate=FAIL`) rather than silently trusted — a
strong illustration of why "passes the literal machine gate" and "physiologically trustworthy"
are not the same claim, and why every caveat is surfaced per-clip in the CSV rather than
collapsed into a single boolean.

## 4. A new finding this aggregate surfaced (not in `docs/MECHANISM_CORPUS_IK.md`)

The upstream pipeline's own `pinned_at_mechanical_limit_flag`/`pinned_joints`
(`pose_to_opensim_ik.py` `main()`) is hardcoded to check only
`("hip_flexion_r","hip_flexion_l","ankle_angle_r","ankle_angle_l")` — **it never checks knee**,
even though `knee_angle_{r,l}`'s own `near_full_mechanical_span` diagnostic is computed and
stored in every report right alongside the others (example (e) above shows knee at 139.9°/140°,
clearly pinned, yet absent from that clip's own `pinned_joints` list). Measured directly from
already-available `angle_stats` data (no upstream pipeline change, no re-run needed): **44/228
clips (19.3%) have knee within 5° of its 140° mechanical ceiling on ≥1 side, silently uncaught
by the existing flag.** Folded into this aggregate's own caveat set
(`knee_near_full_mechanical_span`) so it isn't lost; does not change any CLEAN/DEGENERATE
verdict (checked: none of the 48 no-caveat clips are affected — a genuine coincidence, not
guaranteed to hold on a future, larger corpus). A one-line fix for the source pipeline (add
`"knee_angle_r","knee_angle_l"` to that tuple) is flagged, not applied here — out of scope for
a read-only aggregate.

## 5. Honest gaps / limitations

- **The lenient/strict gap is doing a lot of work.** 181 "CLEAN" vs 67 "strict-CLEAN" — a
  2.7x difference hinging entirely on whether ankle is a hard or soft gate. Anyone consuming
  the 181 number for downstream force work on the ankle/subtalar joint specifically should use
  the strict column instead; anyone working the knee/hip chain (as
  `scripts/msk/force_transmission_scene.py` currently does, on the left leg of the pilot clip)
  is on firmer ground with the lenient definition.
- **73% of lenient-CLEAN clips carry ≥1 caveat** — "CLEAN" here means "passable for
  knee+hip_flexion physiological range, not frozen, well-tracked, majority of frames solved,"
  not "flawless." Only 48/228 (21%) clear every soft flag too.
- **Duration**: corpus spans 1.21 s–172.64 s (median 5.98 s, total 66.3 min). 13 clips run
  under 2.0 s (36-71 frames), 9 of which classify CLEAN. A sub-2-second clip is a thin basis
  for any downstream finite-difference velocity/acceleration estimate even if its ROM and
  marker-fit numbers pass — this aggregate reports duration as a column but does **not** gate
  on it (that wasn't a requested classification criterion and adding an undocumented extra gate
  post-hoc would be exactly the kind of threshold-shopping this method forbids); flagging it
  here for whoever consumes the catalog next.
- **quality is a constant in this batch**: both `corpus_quality` (msk_pose summary.json) and
  the semantic jsonl's own `quality` field read **4** for all 228 clips, and agree with each
  other 228/228 times. The full semantic jsonl (893 records, all 6 IG accounts) has real spread
  ({3: 517, 4: 376}) — this 228-clip batch happens to be entirely drawn from the quality=4
  subset upstream (an input-selection fact, not a computation of this script, and not further
  investigated here).
- **det_frac vs coverage_frac is a real, if currently non-decisive, distinction.** Documented
  in §1 item 3 — no clip's classification currently flips on this (checked: `det_frac<0.5`
  agrees 1-for-1 with the pipeline's own `low_frame_coverage_flag` on the current 228), but a
  future clip with e.g. true `det_frac`≈0.45 and smoothed `coverage_frac`≈0.55 would pass the
  coverage gate on a number that overstates real detection density. Worth a tightened check if
  the corpus grows into that range.
- **No independent ground truth for "physiological."** The range gates are the pipeline's own
  pre-registered anatomical plausibility bounds (backed by the external subject2 real-marker
  RMS floor for marker error), not a per-clip external validation against, e.g., real mocap for
  these specific athletes — that anchor doesn't exist for monocular IG footage and isn't
  claimed here.
- **This is a snapshot of 228 dirs, not a hardcoded count.** `clip_stems()` iterates whatever
  directories currently exist under `data/msk_ik/` — if the batch grows (e.g. toward the
  remaining 36 upstream pose-stage failures, should those get fixed) or shrinks, re-running
  `corpus_ik_aggregate.py` picks up the new set without any code change.

## 6. Reproducing this result

```
python3 scripts/msk/corpus_ik_aggregate.py               # prints full report, writes both artifacts
python3 scripts/msk/corpus_ik_aggregate.py --full-table  # also prints all 228 rows, not just non-clean
python3 scripts/msk/corpus_ik_aggregate.py --no-write    # print only
```

No OpenSim/`.venv-msk` needed — every input number already lives in each clip's own
`_ik_report.json`; this script only does arithmetic over already-computed JSON (plain system
`python3` + numpy, the latter only pulled in transitively via importing
`pose_to_opensim_ik`'s module-level constants for DRY threshold reuse).

Outputs (gitignored, this repo's own `<name>_evidence.json` convention):
`scripts/msk/corpus_ik_aggregate_per_clip.csv` (228 rows, 1 per clip — the scene catalog
itself) and `scripts/msk/corpus_ik_aggregate_evidence.json` (full structured record + summary,
including every self-check's raw output).
