# OrthoLoad AKF force/moment index — per-joint validation anchor (2026-07-21)

Status: BUILT + machine-sanity-checked. Parser: `scripts/msk/index_orthoload_forces.py`
(stdlib-only Python, no new deps). Run over the full corpus: **3,942/3,942 AKF files
parsed, 0 parse errors.** Index lives under `data/external/orthoload/_index/` (gitignored,
same as the source data — regenerate any time by re-running the script; it is fully
deterministic).

This turns the raw OrthoLoad in-vivo force/moment trial corpus (fetched per
`docs/MECHANISM_ORTHOLOAD_FETCH.md`) into a queryable per-trial index plus a joint x
activity summary, and sanity-checks the result against two independent, pre-registered
literature bands (knee level-walking, hip walking) rather than trusting the parse blind.

## 1. What was built

| Output | Path | Rows |
|---|---|---|
| Per-trial index (JSONL, one record/trial) | `data/external/orthoload/_index/orthoload_akf_trials.jsonl` | 3,942 |
| Flagged-trial subset (non-empty QC flags) | `data/external/orthoload/_index/orthoload_akf_flagged.jsonl` | 19 |
| Parse-error log (enumerated, not silently dropped) | `data/external/orthoload/_index/orthoload_akf_parse_errors.jsonl` | 0 |
| Joint x sub_group x activity summary | `data/external/orthoload/_index/orthoload_joint_activity_summary.{csv,json}` | 90 |
| Per-joint distinct-subject index | `data/external/orthoload/_index/orthoload_subject_index.json` | 6 dirs |
| Run manifest (counts, cross-checks, sanity verdicts) | `data/external/orthoload/_index/orthoload_index_run_manifest.json` | — |

Each per-trial record carries: `subject_code`, `months_post_op`, `side`, `implant_type`,
`activity_raw` (verbatim source text) + `activity_bucket` (coarse classification),
`bodyweight_N`, `peak_force_N` (+ 3 independent cross-check values, see §3),
`peak_force_pctBW`, `mean_force_N`, `impulse_N_s`, `duration_s`, `n_data_rows_parsed`,
and a `flags` list — nothing is silently dropped; every one of the 3,942 input files
produced exactly one output record.

## 2. Format facts this parser is built on (measured, not assumed — see script docstring for full detail)

- **Encoding is cp1252, not utf-8.** 240/3,942 files contain a byte>127 (e.g. the degree
  sign). A shell `grep -rl "^BodyWeight"` check **silently failed on exactly those 240
  files** (grep's binary-file heuristic), falsely suggesting 240 files were missing a
  BodyWeight line. Re-checked in Python with an explicit cp1252 decode (which can never
  raise): **BodyWeight is present and parseable in all 3,942/3,942 files.** Recorded here
  because it's a real trap for anyone tempted to grep this corpus directly.
- **Two data-column layouts**, verified by scanning every file's column-header line (not
  a sample): 3,232 files (knee/hip_gen2/shoulder/spine_vbr/spine_fixator) carry
  `Time Fx Fy Fz F Mx My Mz Marker`; the other 710 (**exactly** hip_gen1's file count)
  carry `Time -Fx -Fy -Fz Fres - - - Marker` (moments not recorded for that older
  telemetry generation). The resultant-force column is index 4 in both.
- **Decimal separator is inconsistent within one file**: the "Max. Force (N,sec)" summary
  line uses a comma in newer-format files while the data table below it, in the *same*
  file, uses a period. Handled by normalizing per numeric token, never globally across a
  structural line (a global comma→period replace would break the "CODE, N Months PO"
  parse, which needs that comma to stay a comma).
- Every file has `Number of analog channels = 0` and `Number of Data Sets > 0` (verified,
  full corpus) — no stray extra columns, no empty trials anywhere in the corpus.

## 3. Peak force — three independent cross-checks per trial, not one trusted number

For every trial: (1) the file's own `Max. Force (N,sec)` summary line, (2) the max of the
declared resultant-force data column, (3) an independent recompute of
`sqrt(Fx²+Fy²+Fz²)` from the raw components. A trial is flagged `peak_mismatch` if any
pair disagrees by >2%. **Result: 3,941/3,942 trials agree to well under 2%; exactly one
does not** (see §5). The per-trial index carries all three values plus the disagreement
%, so this is auditable per-trial, not just asserted in aggregate.

## 4. Subject + activity enumeration

Subject code and months-post-op come from the AKF's own `Comment #2` field (two dialects:
`"<CODE>, N Months PO"` for knee/hip_gen2/shoulder/spine, and hip_gen1's own
`"Pat.: <CODE>, Side: <L|R>, N Months PO"`). Activity is a **34-bucket, priority-ordered
keyword classifier** over `Comment #1` (the trial's own free-text description), designed
against a >150-file-per-joint vocabulary sweep rather than guessed upfront, and revised
through 3 rounds of "observe the actual unclassified residual → fix the specific gap →
re-run" (recovered a 16-trial "stand up, no explicit *sit*" gap; separated "Aqua Gym"
into its own bucket after finding it — correctly, buoyancy — was dragging the Walking
bucket's low tail down to ~105% BW). Residual `Other/Unclassified`: **23/3,942 (0.6%)**,
a genuine long tail of one-off activities (agriculture sub-tasks not already bucketed,
bowling, table tennis, dance, calibration movements) — the raw text is always preserved
regardless of bucket, so nothing is lost even where the coarse label is "Other".

Per-joint subject counts (from `orthoload_subject_index.json`, cross-checked against the
live-crawl counts in `docs/MECHANISM_ORTHOLOAD_FETCH.md` §4):

| joint dir | AKF files | subjects found | subjects expected (fetch doc) | match? |
|---|---:|---:|---:|---|
| knee | 609 | 9 (K1L,K2L,K3R,K4R,K5R,K6L,K7L,K8L,K9L) | 9 | ✅ |
| hip_gen1 | 710 | **9** (EBL,EBR,HSR,IBL,JB,KWL,KWR,PFL,RHR) | 2 | ❌ — **see correction below** |
| hip_gen2 | 1,240 | 10 (H1L,H2R,H3L,H4L,H5L,H6R,H7R,H8L,H9L,H10R) | 10 | ✅ |
| shoulder | 140 | 7 (S1R,S2R,S3L,S4R,S5R,S7R,S8R) | 7 | ✅ |
| spine_vbr | 869 | 5 (WP1-WP5) | 5 | ✅ |
| spine_fixator | 374 | 10 (AGL,BBL,FJL,HBL,HSL,JTL,JWL,LGL,MSL,NFL) | 10 | ✅ |

**Correction to `docs/MECHANISM_ORTHOLOAD_FETCH.md`'s hip_gen1 subject count (2 → 9):**
the fetch doc's table says 2 subjects for hip_gen1; this parser independently finds 9,
converging from two directions: (a) direct regex parse of `Comment #2` across all 710
files, and (b) per-code trial counts sum to exactly 710 (185+30+109+68+33+22+132+105+26)
with only 2 files needing a fallback (one filename-prefix recovery for a file whose
`Comment #2` code field was genuinely blank in the source; one initials-with-periods
variant, "E.B." for the same physical patient as the other 184 "EBL" trials, reconciled
via the file's own stated side). Both fallbacks are flagged transparently in the index
(`subject_code_from_filename_fallback`, `subject_code_side_suffix_merge`), not silently
absorbed. This is reported as a correction, not swept under either number.

## 5. Anomalies found and how they were handled (not silently averaged in)

1. **Subject WP4 (spine_vbr), BodyWeight = 100 N in 16/192 of its own trials** — vs
   580 N (126 trials) / 630 N (50 trials) in its other 176. 100 N (~10.2 kgf) is
   implausible for an adult and is almost certainly an OrthoLoad-side placeholder, not a
   remeasurement — recurs identically across 3 different session dates for the same
   subject, never appears for anyone else. **Flagged `bodyweight_implausible` and
   excluded from every %BW aggregate** in the summary table (still present in the
   full per-trial index with the flag set and the raw N values intact).
2. **One peak-mismatch trial**: `hip_gen1/jb4541a.akf`, a "HIP JOINT Stumbling" trial.
   Header says peak 4085.25 N; the full-resolution data column (and an independent
   sqrt(Fx²+Fy²+Fz²) recompute, which agrees with the data column to 0.01 N) says
   4436.58 N, at the *same* timestamp (t=3.095s) — i.e. OrthoLoad's own one-line summary
   statistic disagrees with their own full trace for this one file. Diagnosed, not just
   flagged: this is the corpus's only "Stumbling" (fall-recovery perturbation) trial with
   this issue, plausibly from their export tool applying different filtering to a sharp
   transient; the parser's two independently-computed values agree with each other, so
   the *data column* max (869.9% BW) is used downstream, and the discrepancy is fully
   visible in the record rather than hidden.
3. Two subject-code edge cases in hip_gen1's older (`Measurement Programm 5.0.9`) format
   — see §4 correction above.

**19/3,942 trials (0.5%) carry any flag at all**; the rest parsed clean.

## 6. Sanity check against known literature values (pre-registered bands, from the task brief)

| Check | n trials | median peak %BW | range | literature band | verdict |
|---|---:|---:|---|---|---|
| Knee level walking | 169 | **254.8%** BW | 136.2–353.6% | 250–300% BW | **PASS** (median in band) |
| Hip walking (gen1+gen2) | 382 | **302.4%** BW | 191.8–550.8% | 210–330% BW | **PASS** |
| Hip walking, gen1 only | 201 | 296.7% BW | 191.8–550.8% | 210–330% BW | PASS |
| Hip walking, gen2 only | 181 | 305.5% BW | 218.0–469.2% BW | 210–330% BW | PASS |

Both headline medians land inside the task-specified bands. Being symmetric about it:
only 46% of individual knee-"Walking" trials and 66% of hip-"Walking" trials fall inside
the tight literature band — the median-PASS is real but the band is narrower than the
full observed spread, and that spread is not noise. Splitting the knee "Walking" bucket
by sub-condition shows exactly why:

| knee walking sub-condition | n | median %BW |
|---|---:|---:|
| level walking (free) | 73 | 250.9% |
| brace/wedge/orthosis | 48 | 247.8% |
| treadmill (no device) | 35 | 268.7% |
| barefoot | 7 | 288.0% |
| other (gym-floor/MBT/etc.) | 9 | 229.3% |

All five sub-conditions cluster in a tight 230–290% band — it's the pooled "Walking"
label that looks wide, not the underlying free-gait data. The trials actually dragging
the pooled range down/up are physically explained, not anomalous: **"Aqua Gym" walking
trials** (buoyancy offloads body weight — correctly given their own bucket after this was
found, see §4) at the low end, and, on the hip side, **"walking on treadmill + carrying
5 kg in the hands"** at the high end (added mass → added joint load, exactly the expected
direction). Both are real, physiologically-sensible signals, not parser noise.

## 7. Full joint x activity table (top-level joint, pooling sub_group — see the CSV for the gen1/gen2 and vbr/fixator split)

### knee
| Activity bucket | #trials | #subjects | peak %BW median | peak %BW range |
|---|---:|---:|---:|---:|
| Walking | 169 | 9 | 254.8% | 136.2-353.6% |
| Whole-Body Vibration | 62 | 6 | 201.4% | 144.5-374.2% |
| Stairs Up | 54 | 9 | 329.6% | 221.9-484.8% |
| Stairs Down | 52 | 9 | 348.1% | 284.1-458.9% |
| Aqua Gym (buoyancy-supported) | 50 | 5 | 137.7% | 40.6-414.5% |
| Knee Bend/Squat | 47 | 9 | 249.8% | 158.5-359.6% |
| Sit-to-Stand/Stand-to-Sit | 47 | 9 | 268.0% | 133.1-408.7% |
| Cycling | 45 | 9 | 132.0% | 58.2-211.4% |
| One-Leg Stance | 32 | 9 | 280.3% | 202.6-376.6% |
| Swimming | 21 | 4 | 105.6% | 69.4-189.5% |
| Standing/Stance | 16 | 8 | 91.6% | 56.6-154.6% |
| Stairs (Up+Down) | 6 | 5 | 334.0% | 312.8-388.3% |
| Jogging/Running | 4 | 2 | 621.5% | 577.1-633.8% |
| Trampoline/Jumping | 4 | 2 | 547.5% | 472.4-621.8% |

### hip
| Activity bucket | #trials | #subjects | peak %BW median | peak %BW range |
|---|---:|---:|---:|---:|
| Walking | 382 | 19 | 302.4% | 191.8-550.8% |
| Lying | 294 | 14 | 137.8% | 22.5-438.9% |
| Exercise/Gym | 215 | 9 | 251.0% | 66.9-492.1% |
| Sit-to-Stand/Stand-to-Sit | 119 | 17 | 184.2% | 74.2-403.1% |
| One-Leg Stance | 111 | 15 | 281.6% | 83.5-523.9% |
| Cycling | 111 | 15 | 99.7% | 55.9-211.4% |
| Standing/Stance | 106 | 17 | 223.1% | 69.7-512.1% |
| Sitting | 85 | 7 | 197.2% | 23.5-420.3% |
| Walking (crutch-assisted) | 57 | 12 | 220.2% | 75.0-324.9% |
| Aqua Gym (buoyancy-supported) | 57 | 4 | 158.8% | 91.9-439.4% |
| Lifting | 56 | 14 | 180.8% | 57.0-601.7% |
| Whole-Body Vibration | 47 | 4 | 131.7% | 92.6-375.3% |
| Swimming | 44 | 6 | 187.9% | 101.7-306.9% |
| Stairs (Up+Down) | 41 | 12 | 305.5% | 171.3-565.3% |
| Jogging/Running | 38 | 9 | 467.8% | 373.4-634.1% |
| Knee Bend/Squat | 37 | 12 | 202.0% | 100.0-587.6% |
| Stairs Up | 36 | 16 | 300.2% | 162.4-555.4% |
| Stairs Down | 34 | 15 | 309.4% | 164.6-529.6% |
| Other/Unclassified | 23 | 10 | 402.1% | 249.4-673.7% |
| Occupational/Agricultural Task | 20 | 1 (single-subject bucket) | 474.8% | 295.5-520.5% |
| Trampoline/Jumping | 14 | 8 | 428.4% | 248.8-554.1% |
| Car Transfer | 12 | 3 | 223.0% | 198.6-297.1% |
| Stumbling (perturbation) | 5 | 3 | 338.1% | 300.3-869.9% (max = the flagged peak-mismatch trial, §5) |
| Rowing | 4 | 4 | 308.6% | 223.3-380.8% |
| Bed Transfer | 2 | 1 | 163.2% | 152.1-174.3% |

### shoulder
| Activity bucket | #trials | #subjects | peak %BW median | peak %BW range |
|---|---:|---:|---:|---:|
| Walking (crutch-assisted) | 52 | 7 | 79.3% | 39.9-197.1% |
| Standing/Stance | 36 | 6 | 88.0% | 36.3-238.2% |
| Wheelchair | 25 | 6 | 88.7% | 20.7-115.3% |
| ADL (Daily Living Tasks) | 20 | 4 | 94.6% | 58.2-183.1% |
| Lifting | 6 | 3 | 106.0% | 96.1-113.2% |
| Lying | 1 | 1 | 29.2% | 29.2-29.2% |

Shoulder %BW is naturally much lower than knee/hip — it is not a primary weight-bearing
joint; load is transferred through the arm only during crutch-support/ADL tasks, exactly
the pattern above (crutch-walking and standing/ADL dominate the trial count).

### spine (vbr + fixator pooled)
| Activity bucket | #trials | #subjects | peak %BW median | peak %BW range |
|---|---:|---:|---:|---:|
| Trunk Flexion/Extension/Rotation | 271 | 15 | 29.8% | 5.2-170.6% |
| Lifting | 210 | 13 | 37.3% | 2.5-252.3% |
| Standing/Stance | 149 | 12 | 72.2% | 8.0-238.9% |
| Walking | 120 | 11 | 58.0% | 8.3-169.5% |
| Sit-to-Stand/Stand-to-Sit | 76 | 6 | 89.0% | 10.8-159.4% |
| Whole-Body Vibration | 75 | 4 | 49.7% | 11.1-94.8% |
| Lying | 71 | 13 | 22.8% | 1.4-130.0% |
| Arm Elevation | 68 | 10 | 54.9% | 7.2-146.4% |
| Sitting | 61 | 10 | 42.3% | 6.3-139.6% |
| Kneeling | 42 | 7 | 36.1% | 5.5-107.9% |
| Cycling | 29 | 6 | 39.9% | 6.8-62.2% |
| ADL (Daily Living Tasks) | 23 | 4 | 110.8% | 29.5-173.8% |
| Changing Position (transfer) | 17 | 4 | 60.7% | 19.0-68.8% |
| Exercise/Gym | 14 | 1 (single-subject bucket) | 15.5% | 5.8-52.9% |
| Jogging/Running | 5 | 2 | 79.3% | 63.5-81.8% |
| Walking (crutch-assisted) | 5 | 3 | 59.8% | 10.0-124.4% |
| Knee Bend/Squat | 4 | 2 | 89.8% | 78.8-103.5% |
| Wheelchair | 3 | 1 | 97.3% | 72.4-110.6% |

Spine %BW values are an order of magnitude lower than hip/knee across every activity —
mechanically expected (an instrumented vertebral segment does not carry anything close
to full body weight the way a lower-limb joint does in single-leg stance; axial spinal
load is shared across the whole column and paraspinal musculature). **No literature band
was pre-registered for spine/shoulder in this task** (only knee/hip were specified) — the
numbers above are reported descriptively; the directional pattern (lifting/trunk-flexion
> standing > lying, and standing > lying) is biomechanically sensible but was not
checked against an external spine-load reference. Flagged as an honest gap, §8.

## 8. Honest gaps / not done

- **No literature anchor for shoulder or spine** — only knee level-walking and hip
  walking were pre-registered bands in this task. The spine/shoulder numbers above are
  real, machine-computed, and internally cross-checked, but not externally verified
  against a published reference the way knee/hip were.
- **`Other/Unclassified` residual: 23/3,942 (0.6%)** — a long tail of low-count one-off
  activities (occupational sub-tasks, bowling, table tennis, dance, calibration
  movements). Raw text is fully preserved; further bucket-splitting would have single or
  low-double-digit n per new bucket, judged not worth the added taxonomy complexity here.
- **IOF/COF/EOF files (996 more force-like trial files, same header convention per the
  fetch doc) were explicitly out of scope** for this task (which asked for the AKF corpus
  specifically) and were not touched — a natural, cheap follow-on if the twin needs them,
  since the same parser's header/data-block logic should apply directly.
- **Tier-A `standard_loads/*.xlsx` files** (the separate, already-tidied per-subject-x-
  activity knee/hip spreadsheets) were not folded in — different format, different
  (already-averaged) granularity, out of this task's AKF-specific scope.
- **The `Occupational/Agricultural Task` (hip, 20 trials) and `Exercise/Gym` (spine, 14
  trials) buckets are single-subject** — real data, but not cross-subject-generalizable;
  flagged here so nobody mistakes n_trials for subject diversity without checking
  n_subjects (which the summary table always reports alongside it for exactly this
  reason).
- **hip_gen1 vs hip_gen2 are structurally different** (gen1 has no moment channels, older
  export format) — kept as separate `sub_group` values throughout the index rather than
  silently merged, so anyone using this for a physics fit can choose to pool or not.

## 9. Next step

The per-trial JSONL is ready to serve as the held-out validation anchor for the reduced-
representation generative-design generality claim (peak %BW by activity, per subject,
per joint) — e.g. compare a generated/reduced-order knee model's predicted peak contact
force under "Level Walking" against the 169-trial, 9-subject empirical distribution above
(254.8% BW median) rather than a single textbook number. Natural follow-on: fold in
IOF/COF/EOF once/if the twin needs the additional force-rig variants, and pull a hip/knee
literature-band equivalent for spine and shoulder if those joints become load-bearing
validation targets too.
