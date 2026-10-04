# MECHANISM WAVEFORM CROSS-SUBJECT — does the push-off-concentration reframe generalize? (2026-07-21)

`docs/MECHANISM_CONTACT_WAVEFORM.md` found, on subject2/`walking1` alone, that the knee/hip
over-prediction is **concentrated at push-off**: SO-knee's weight-acceptance hump (245.4 %BW) is
statistically indistinguishable from real OrthoLoad in-vivo data (224.6±21.6 %BW, Z=0.97) while its
push-off hump (390.3 %BW vs 266.9±46.2 %BW, Z=2.67) is a clear over-prediction. n=1 subject was the
open gap that doc itself flagged. This doc runs the identical waveform pipeline
(`scripts/msk/contact_waveform_analysis.py`) on **subject3** and **subject4** (whose SO solves
already existed, `docs/MECHANISM_CROSS_SUBJECT.md`) to test whether that split is a cohort-wide
property or subject2-specific.

**Headline (pre-registered falsifier from the task, applied literally): the cohort-wide "generalizes"
criterion is NOT met.** Of the 4 new series tested (subject3 knee/hip, subject4 knee/hip), only
**1 of 4** (subject3 SO-knee) reproduces subject2's clean pattern (weight-acceptance Z<1.5, push-off
Z>2); the other **3 of 4** (subject3 hip, subject4 knee, subject4 hip) over-predict at
weight-acceptance TOO (Z=2.77–4.49, all clearing the 1.5 bar) — Outcome B, "subject2-specific,"
per the task's own pre-registration. **But the weaker, DIRECTIONAL claim generalizes perfectly**:
push-off's ratio/Z exceeds weight-acceptance's own ratio/Z in **8 of 8** series now measured across
all 3 subjects (4 from subject2 + 4 here) — a mechanism that always shows up, just not always big
enough to leave weight-acceptance statistically invisible. The refined, honest claim: *the
"push-off is relatively worse than weight-acceptance" tendency is cohort-wide and robust; the
stronger "weight-acceptance is a genuine, statistically-invisible validation point" is fragile —
so far seen only for knee (2 of 3 subjects: subject2, subject3; NOT subject4), never for hip (0 of
3 subjects, though subject2's own hip Z=1.72 was already the mildest/most-borderline case even
there).**

A forced adversary surfaced and resolved along the way: subject3's official knee "hump1"
(250.8 %BW @ 29%GC) turned out to be a boundary-clipped ramp value, not a genuine peak (confirmed:
no real dip between it and hump2) — but the genuine local-peak alternate reading (225.1 %BW @
22%GC, Z=0.74) **still lands inside the "matches" verdict**, just as cleanly as the flawed official
number (Z=1.22) — the artifact did not manufacture subject3-knee's positive result, it was
independently reconfirmed.

## Headline numbers

| Series | hump1 (weight-accept.) model/ortho±SD | ratio₁ | Z₁ | hump2 (push-off) model/ortho±SD | ratio₂ | Z₂ | rel.spread | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| subject2 SO knee (reference, established) | 245.4 / 224.6±21.6 | 1.09× | 0.97 | 390.3 / 266.9±46.2 | 1.46× | 2.67 | 28.9% | CONCENTRATED |
| subject2 SO hip (reference, established) | 326.3 / 249.9±44.4 | 1.31× | 1.72 | 386.9 / 232.2±34.0 | 1.67× | 4.56 | 24.3% | UNIFORM (borderline) |
| **subject3 SO knee** (official, boundary-artifact*) | 250.8 / 224.6±21.6 | 1.12× | **1.22** | 512.1 / 266.9±46.2 | 1.92× | **5.30** | 52.8% | CONCENTRATED |
| *subject3 SO knee (genuine-peak ALT reading, 22%GC)* | *225.1 / 205.2±26.9* | *1.10×* | *0.74* | *(same)* | *(same)* | *(same)* | *64%* | *CONCENTRATED* |
| subject3 SO hip | 426.2 / 249.9±44.4 | 1.71× | **3.97** | 455.2 / 232.2±34.0 | 1.96× | 6.57 | 13.9% | UNIFORM |
| subject4 SO knee | 284.3 / 224.6±21.6 | 1.27× | **2.77** | 433.1 / 266.9±46.2 | 1.62× | 3.59 | 24.7% (borderline) | UNIFORM |
| subject4 SO hip | 449.1 / 249.9±44.4 | 1.80× | **4.49** | 418.4 / 232.2±34.0 | 1.80× | 5.48 | 0.25% | UNIFORM |

\* see §3.2 — no real dip exists between subject3-knee's reported hump1 and hump2 (the curve rises
monotonically through the fixed 30%GC split); flagged, not hidden, and independently reconfirmed
via the genuine-local-peak alternate reading on the row below.

**Bold Z-values ≥1.5** are the ones that fail the task's own pre-registered weight-acceptance-match
bar; only subject2-knee (0.97) and subject3-knee (1.22 / 0.74 alt) clear it.

**Confidence tier: in-vivo-anchored for knee** (OrthoLoad's official `Dat_BW` per-subject curves, 8
subjects, calibration r=0.9278483314012478 — reproduces subject2's own stored 0.9278483314012478
**exactly**, confirming the re-import introduced no drift) — **in-vivo-anchored-but-self-built for
hip** (same disclosed limitation as `docs/MECHANISM_CONTACT_WAVEFORM.md` §6: the official OrthoLoad
hip archive is a corrupted download; the hip ensemble is this repo's own raw-AKF build, calibrated
against official knee only, shared unchanged across all 3 subjects so it cannot itself explain any
between-subject difference).

## 1. Pre-registration (from the task, stated before any subject3/4 hump number was computed)

**Outcome A ("generalizes")**: weight-acceptance Z<~1.5 AND push-off Z>2 at each new subject/joint
→ the twin is validated at weight-acceptance cohort-wide, a genuine positive cross-subject result.
**Outcome B ("subject2-specific")**: subject3/4 over-predict at weight-acceptance too → subject2's
clean split was subject-specific, the reframe is weaker. The task named subject3's 1.99× knee
ratio as "a strong candidate to break it — it may over-predict everywhere."
**Symmetric-QC clause (also pre-registered)**: subject3/4's SO solves failed the raw
`convergence_pass` gate on shorter trials (`docs/MECHANISM_CROSS_SUBJECT.md` §4) — verify each has
enough gait-cycle coverage for a REAL weight-acceptance AND push-off hump before trusting the split;
if too short, say so as a diagnosed gap, not a forced number.

**What actually happened, measured**: subject3-knee → Outcome A. subject3-hip, subject4-knee,
subject4-hip → Outcome B. Interestingly, the task's own named "strong candidate to break it"
(subject3's knee) is the ONE series that does NOT break it — the break instead comes from hip (both
subjects) and subject4's knee, which the task did not specifically flag. A real, not
pattern-matched, result.

## 2. Method (reuse discipline)

`scripts/msk/contact_waveform_cross_subject.py` (new) imports
`scripts/msk/contact_waveform_analysis.py` (cwa) as a module and calls its functions **unedited** —
md5 `9d54780f5146b8fea3821b3f92bac18d` before AND after this session's runs (verified in-script,
asserted, not just narrated). Every subject-independent piece (OrthoLoad ensembles, calibration,
`find_two_humps`, `compare`, `smooth5`) is called completely unchanged and computed **once**, shared
across subjects (lean — the external cohort does not depend on which LabValidation subject is being
tested). No CMC comparison exists for subject3/4 (confirmed: no `cmc_second_solve/` output under
either subject's `data/msk_smoketest/` tree, unlike subject2) — SO only, disclosed, not forced.

**Two forced adversaries caught by reading cwa's full call chain before trusting it on new
subjects, not assumed safe because it worked on subject2:**

**Adversary #1 (a silent-corruption landmine).** `cwa.model_curve_on_grid()` hardcodes a call to
`joint_force_mag_pct_bw(reaction_sto, colprefix)` with no `bw_n` passed through. Python binds a
default-parameter value (`bw_n=BW_N`) **once at function-definition time**, not at call time —
unlike a bare-name lookup inside a function body (the mechanism `cross_subject_validation.py`
correctly relies on for `sok.SUBJECT_DIR` etc.), monkey-patching `cwa.BW_N` after import is
**inert** for this call chain. Verified live with a standalone toy closure before trusting this
conclusion (confirmed: patching the module constant after def-time left an already-bound default
unchanged; only an explicit keyword override picked up the new value). Reusing
`model_curve_on_grid()` verbatim would have silently normalized subject3/subject4 by subject2's
BW_N=766.88N instead of their own 622.72N/613.90N — an 18.8%/19.6% systematic %BW
UNDER-statement, corrupting every downstream ratio and Z-score with no error raised. **Fix (not an
edit to cwa.py):** call `cwa.joint_force_mag_pct_bw(reaction_sto, colprefix, bw_n=subject_bw_n)`
directly with an explicit override; `build_model_curve()` in the new script reimplements the
remaining ~6 lines (grid resample + peak metadata) inline. Verified via self-consistency: all 4
new series reproduce the already-published `cross_subject_validation.py` JR peak numbers to
**0.0000% relative difference** (513.75/477.93 %BW subject3 knee/hip; 434.15/497.44 %BW subject4
knee/hip) — confirms the correct per-subject bw_n is actually being used.

**Adversary #2 (a hard assert that blocks subject3 entirely).**
`cwa.model_pct_gc_mapping()` asserts exactly one each of {R heel-strike, R toe-off, L heel-strike, L
toe-off} inside the GRF window. Measured: subject4 satisfies this (all 4 present) — called **truly
unchanged**, zero adaptation, confirmed by running it directly. subject3 does not: L toe-off never
crosses back below the 20N threshold inside the 1.319s window
(`L_ground_force_vy(t=1.319)=21.95N`, still above threshold — toe-off is imminent but uncaptured).
Diagnosed via OODA, not silently patched: the assert demands more than the transform itself needs —
`T_stride` only requires R_HS+L_HS; L_TO is consumed only by the later, separate out-of-sample
verification check. `robust_model_pct_gc_mapping()` (new) reimplements the identical
contralateral-symmetry recipe (same 20N threshold, same edge detector, same `T_stride` formula)
requiring only R_HS+L_HS — **cross-validated to reproduce `cwa.model_pct_gc_mapping()`'s own diag
dict to a max absolute difference of exactly 0.0** on subject2 (published ground truth) and
subject4 (where cwa's strict function also runs), *before* being trusted standalone on subject3
where no ground truth exists to check it against.

## 3. Results in detail

### 3.1 Is each trial long enough to contain a REAL instance of both humps? (the task's own gate)

Measured directly (SCENE-EYES — a number, not an assumption from trial duration): for every model
frame, map it to %GC and check how many real frames land inside each hump's search band
(`find_two_humps`' own split: hump1∈[0,30)%GC, hump2∈[30,62]%GC).

| Subject | n frames | max %GC gap (real-data coverage hole) | hump1-band real frames | hump2-band real frames |
|---|---:|---:|---:|---:|
| subject3 | 132 | 0.91%GC (starts@23.7%GC) | 55 | 35 |
| subject4 | 131 | 0.94%GC (starts@32.1%GC) | 56 | 34 |

Both trials have dense real-frame support in BOTH hump bands and no coverage hole bigger than 1%GC
— **not** a "too short" case. This works out because `T_stride` (1.095s/1.063s) is almost as long
as the trial itself (1.31s/1.30s) with the captured heel-strike sitting well inside the window, so a
single continuous recording sweeps almost the entire 0–100%GC range even though (like subject2) only
one true heel-strike-to-heel-strike cycle is ever directly observed.

### 3.2 Phase-alignment verification: the narrow, subject2-tuned band misses, but two independent anchors hold

The out-of-sample textbook-band check (`cwa.verify_phase_alignment`, reused UNCHANGED, same 55-68%/
8-20% thresholds as subject2 — no re-tuning):

| Subject | R toe-off %GC (band 55-68%) | L toe-off %GC (band 8-20%) | raw `all_pass` |
|---|---:|---:|---|
| subject2 (established) | 66.9% — PASS | 17.0% — PASS | **True** |
| subject3 | 70.1% — **FAIL** | unavailable (event not captured in-window) | **False** (1 of 2 checks available, fails) |
| subject4 | 69.6% — **FAIL** | 19.7% — PASS | **False** |

Reported raw, not hidden. Forced via OODA rather than accepted as an honest negative — two
independent, decorrelated anchors were checked before trusting (or discarding) the phase mapping:

1. **T_stride plausibility vs the OrthoLoad raw-AKF corpus** (n=371 knee-derived cycles,
   mean=1.0947s SD=0.1252s; n=757 hip-derived cycles, mean=1.0997s SD=0.1251s — external data,
   nothing to do with this trial's own phase mapping): subject3's `T_stride`=1.095s sits **0.002 SD**
   from the corpus mean; subject4's 1.063s sits **0.25-0.29 SD** away. Subject2's own already-trusted
   `T_stride`=1.327s sits **1.82-1.86 SD** away — i.e., by this external, independent metric,
   **subject2's own mapping is the least-typical of the three**, not subject3/4's.
2. **Independent peak-time cross-check**: `cross_subject_validation.py` (a separate script,
   computed before this task, using none of this doc's phase-mapping code) already established
   subject3's knee peak at t=0.49s and subject4's at t=0.46s. Mapped through THIS doc's own
   `T_stride`, those land at **47.4%GC** and **46.2%GC** respectively — squarely inside subject2's
   own established push-off window (44.8-47.8%GC), exactly where the known gastroc/soleus-driven
   push-off mechanism predicts the global maximum should sit.

Both anchors support the mapping being fundamentally sound; the narrow textbook band (fit to n=1
subject2) simply doesn't transfer — a modest (1.6-2.1 percentage point) miss in the SAME direction
for both new subjects, most plausibly ordinary inter-subject stance-fraction variation (a
well-documented source of normal gait variability) rather than a broken transform. Diagnosed gap,
not swept under the rug: subject3's L-toe-off cross-check is **structurally unavailable**, not
failed — one of the two independent checks subject2 had is simply not observable in this trial.

### 3.3 Hump-shape validity: a forced adversary, caught and resolved

`cwa.find_two_humps`' fixed-30%GC-split argmax hump-finder was verified, per its own docstring,
only against subject2's 6 curves. That verification does not automatically transfer. Machine check
added (`has_real_dip_between`, new): is there a real dip (value below hump1's own) strictly between
the reported hump1 and hump2? Answer for all 4 new series:

| Series | Real dip between humps? |
|---|---|
| subject3 SO knee | **NO** — monotonic ramp, hump1 is a boundary-clipped value |
| subject3 SO hip | YES — genuine double-hump |
| subject4 SO knee | YES — genuine double-hump |
| subject4 SO hip | YES — genuine double-hump |

Only subject3-knee fails this check. Diagnostic alternate reading (`find_genuine_local_peak`, new):
the highest genuine local maximum (real decline on both sides) in subject3-knee's curve is
225.1 %BW @ 22%GC — vs OrthoLoad 205.2±26.9 %BW there, ratio 1.097×, **Z=0.74**. This is EVEN
CLEANER evidence of a weight-acceptance match than the flawed official reading (Z=1.22) — the
artifact did not manufacture the positive result; forcing it through an independent, decorrelated
detection method (genuine local-max-with-decline, vs. plain fixed-window argmax) reconfirms it.

### 3.4 The core question: does weight-acceptance stay invisible while push-off doesn't?

Applying the task's own pre-registered bar (Z₁<~1.5, Z₂>2) to all 4 new series plus the 4
already-established subject2 series (8 total, the full cross-subject-x-joint-x-solver space tested
to date):

- **Series clearing Z₁<1.5 (weight-acceptance "matches")**: subject2-SO-knee (0.97), subject3-SO-knee
  (1.22 official / 0.74 alt). **2 of 8** — both knee, both SO solver, 2 of the 3 subjects.
- **Series NOT clearing it**: subject2-CMC-knee (3.30), subject2-SO-hip (1.72 — already over the
  bar even for the reference subject), subject2-CMC-hip (3.23), subject3-SO-hip (3.97),
  subject4-SO-knee (2.77), subject4-SO-hip (4.49). **6 of 8.**
- **Series with Z₂>Z₁ (push-off relatively worse than weight-acceptance)**: **8 of 8** — perfect,
  including all 4 new series (subject3 knee 5.30>1.22; subject3 hip 6.57>3.97; subject4 knee
  3.59>2.77; subject4 hip 5.48>4.49, the last one nearly tied at 5.48 vs 4.49 but still directionally
  correct).

**Reading this honestly**: the pre-registered cohort-level "generalizes" criterion (Z₁<1.5 at
EVERY new subject/joint) fails — 3 of 4 new series over-predict at weight-acceptance too, Outcome
B. But it is not a clean negative either: hip NEVER passed the Z₁<1.5 bar in ANY subject including
subject2 itself (1.72, 3.23), so hip was never a fair test of the ORIGINAL "weight-acceptance
matches" headline (that headline was always specifically a **knee, SO-solver** finding) — and for
knee-SO specifically, 2 of 3 subjects DO replicate it. The directional mechanism (push-off is
always relatively worse) is the part that is robustly cohort-wide (8/8); the specific "weight-acceptance
is a genuine validation point" claim from `docs/MECHANISM_CONTACT_WAVEFORM.md` should be narrowed
to "sometimes true for knee via Static Optimization, not reliable across joints or subjects" rather
than treated as an established, generalized property of the twin.

## 4. Honest gaps

1. **n=3 subjects, n=4 new series** — still a small, descriptive instance-space (same scope caveat
   `docs/MECHANISM_CROSS_SUBJECT.md` already carries). Sufficient to show the "generalizes" criterion
   fails at the cohort level and to identify the knee/hip split; insufficient to fit a population
   model of WHY knee-SO sometimes clears the bar and hip never does (candidate explanations —
   per-joint muscle-strength-to-bodyweight scaling, hip's self-built/uncalibrated-on-itself OrthoLoad
   ensemble, genuine anatomical/mechanical differences — are not disentangled here).
2. **subject3's L-toe-off cross-check is structurally unavailable**, not failed — the trial ends
   while L-foot loading is still above the 20N stance threshold. This is a genuine data-length
   limit on ONE of the two independent phase-verification events, not on the hump comparison itself
   (§3.1 already shows dense real-frame support at both humps regardless).
3. **The R-toe-off textbook band (55-68%GC) misses for BOTH new subjects** (70.1%, 69.6%) — same
   direction, similar magnitude. Two independent anchors (§3.2) support the mapping being sound
   regardless, but the root cause of the consistent 1.6-2.1-point miss (ordinary stance-fraction
   variation vs. a small systematic T_stride bias in the contralateral-symmetry convention itself)
   is not distinguished here — would need a subject with a directly-observed second heel-strike to
   settle definitively.
4. **Hip's OrthoLoad ensemble is self-built and shared, unchanged, across all 3 subjects** — the
   same limitation `docs/MECHANISM_CONTACT_WAVEFORM.md` §4 already discloses for subject2 (calibrated
   against official KNEE data only, r=0.9278, never independently checked against an official hip
   curve because none exists). Since it is IDENTICAL across all 3 subjects' hip comparisons, it
   cannot manufacture a between-subject difference, but any absolute bias in it (over- or
   under-smoothing the true hip waveform) applies to all three alike.
5. **subject3/4's own raw `convergence_pass` gate reads False** (the symmetric-QC clause the task
   pre-registered) — already diagnosed in `docs/MECHANISM_CROSS_SUBJECT.md` §4 as a subject2-specific
   hardcoded `t≥1.50s` literal that fails BY CONSTRUCTION for every other subject in this cohort
   (subject2's trial is the longest of ten), not a real non-convergence; the subject-agnostic
   solve-quality read (frames>100, no NaN, in-bounds activations) is True for both. Carried forward
   unchanged here, not re-litigated.
6. **No CMC comparison for subject3/4** (no second-solve output exists) — SO only, half the series
   subject2 had. Disclosed, not forced.
7. **Only right side, only `walking1`** — same scope caveat inherited from every cert in this family.

## 5. Verification (machine-checked, not narrated)

- `scripts/msk/contact_waveform_analysis.py` integrity: md5 `9d54780f5146b8fea3821b3f92bac18d`
  before AND after this session's runs (asserted in-script; a mismatch would raise).
- Robust-mapping-vs-strict-cwa cross-validation: max absolute difference over 7 diagnostic keys =
  **0.0 exactly** on both subject2 and subject4 (asserted in-script; `gate_pass: true` in the JSON).
- OrthoLoad knee calibration reproduces subject2's own already-published value **exactly**
  (0.9278483314012478 both places, not just close).
- Self-consistency: all 4 new series' re-extracted overall-peak values match the already-published
  `cross_subject_validation.py` JR numbers to **0.0000%** relative difference (asserted in-script).
- Hump-shape validity: `has_real_dip_between` is a hard, pre-specifiable, machine-evaluated
  falsifier (does a value below hump1 exist strictly between the two reported hump locations?), not
  a visual/narrated judgment.
- All numbers in this document are pulled directly from
  `data/msk_smoketest/contact_waveform_cross_subject/contact_waveform_cross_subject_results.json`,
  cross-checked against `data/msk_smoketest/subject2_walking1/contact_waveform/contact_waveform_results.json`
  and `data/msk_smoketest/cross_subject_validation/subject{3,4}_walking1_cross_subject_results.json`
  for the established reference numbers — never transcribed from console prose.

## 6. Files

- `scripts/msk/contact_waveform_cross_subject.py` (new) — the orchestration script described above.
  Run with `/usr/bin/python3 scripts/msk/contact_waveform_cross_subject.py` (needs numpy+openpyxl,
  same as `contact_waveform_analysis.py`; deliberately not `.venv-msk`; never imports OpenSim).
- `docs/MECHANISM_WAVEFORM_CROSS_SUBJECT.md` — this doc.
- `data/msk_smoketest/contact_waveform_cross_subject/contact_waveform_cross_subject_results.json`
  (new) — every number in this document: robust-mapping cross-validation, OrthoLoad ensembles
  (computed once), per-subject mapping diagnostics, phase-alignment checks, coverage diagnostics,
  hump comparisons + hump-shape-validity diagnostics, for both subject3 and subject4.
- **Not touched** (md5-verified, §5): `scripts/msk/contact_waveform_analysis.py`. Also not touched
  or relied upon: `scripts/msk/contact_waveform_fmax_corrected.py` and
  `data/msk_smoketest/subject3_walking1/fmax_pcsa_correction_test/` — confirmed (via docstring +
  mtime, 2026-07-21, same day) to be concurrent, independent work on a different question
  (Fmax/PCSA correction), not read in detail or depended on here.
- Inputs read in place, never modified: `data/msk_smoketest/subject{3,4}_walking1/static_optimization/jr/`
  (model force time series), `data/external/orthoload/...` (OrthoLoad, subject-independent, shared),
  `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject{3,4}/ForceData/walking1_forces.mot`
  (GRF, heel-strike detection), `data/msk_smoketest/cross_subject_validation/subject{3,4}_walking1_cross_subject_results.json`
  (established reference numbers for self-consistency).

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
