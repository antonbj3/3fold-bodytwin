# MECHANISM CROSS-SUBJECT VALIDATION — does the knee/hip force cert generalize beyond subject2?

**Question.** Every joint-force cert in this repo (`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`) is validated on
exactly one subject: subject2/walking1. `docs/MECHANISM_SIGN_BUG_REMEDIATION.md` Sec.10 names this explicitly
as an open gap: *"Single trial, right side, subject2 walking1 — not re-verified across subjects."* This doc
closes that gap by running the identical, unmodified cert pipeline (Static Optimization + the sign-fixed
crossing-muscle machinery + official `opensim.JointReaction`) on **subject3** and **subject4** — two different
bodies (63.5 kg/1.69 m/F and 62.6 kg/1.68 m/F, vs. subject2's 78.2 kg/1.96 m/M) — and comparing their
knee/hip over-prediction ratio against subject2's own established **1.515× (knee) / 1.412× (hip)**.

**Headline answer.** The **direction** of the over-prediction generalizes robustly and strengthens the
Moissenet-mechanism explanation: all 3 subjects over-predict at both joints, by a wide margin, and the two
independent computation methods (self-computed, official JointReaction) keep agreeing with each other to
<0.01% at the knee in every subject — this is a property of the SO-then-subtract method architecture, not
noise or a subject2-specific code artifact. **But the specific "~1.5×" magnitude does NOT generalize** —
subject2 is the **mildest** of the 3 subjects tested, not a representative or worst case. Subject3/4 over-predict
knee contact force by **1.68–1.99×** (vs subject2's 1.515×) and hip by **1.53–1.82×** (vs subject2's 1.412×).
The honest, updated claim is *"the twin systematically over-predicts in-vivo knee/hip contact force via this
architecture, by a factor that is subject-dependent and measured so far in the range ~1.4–2.0×"* — not a
single universal 1.5× constant.

A second, unplanned finding surfaced while forcing the hip result: subject3/4's shorter mocap captures expose
a **latent convention mismatch already present in the existing, unmodified cert code** (self-computed excludes
edge frames from its peak search; the official JointReaction extractor does not) that happened to be invisible
on subject2's longer trial. Diagnosed, not hidden — see Sec.3.

---

## 1. Method (reuse discipline)

Zero bytes of `validate_joint_force.py` / `static_opt_knee.py` / `validate_hip_force.py` were edited — verified
by md5 before and after every run (Sec.5). All algorithmic steps (SO run, convergence gates, the sign-**fixed**
`knee_crossing_muscles_and_forces`, `compute_self_cross_check(_generic)`, `JointReaction`, the OrthoLoad AKF
parsers) are the same imported functions, called unchanged. The new script, `scripts/msk/cross_subject_validation.py`,
does two things that cannot pre-exist:

1. **Orchestration**: monkey-patches the three modules' own module-level path constants
   (`SUBJECT_DIR`/`MODEL_FILE`/`IK_MOT`/`GRF_MOT`/`SO_OUT`/`JR_OUT`/etc.) to point at a new subject's files —
   Python resolves bare names against the enclosing module's `__dict__` at call time, so this repoints an
   unedited function exactly the way `validate_hip_force.py` already reuses `static_opt_knee.py`'s functions
   cross-module. `TRIAL_END_TIME` is **re-derived live** per subject from that subject's own real IK/GRF
   end-times (never copied from subject2's 1.57s) — the same "tighter of the two real ranges, never
   extrapolate" rule `static_opt_knee.py`'s own docstring states, applied to new data.
2. **Cross-subject comparison logic**, which by definition is new (this is the first script in the repo that
   runs the same cert on more than one subject).

The OrthoLoad in-vivo anchor is **not** subject-specific — it is a different, external, real-patient
instrumented-implant cohort — so it was correctly re-parsed identically for every subject run
(knee median 258.221 %BW n=72, hip median 273.931 %BW n=162 — both matched subject2's own cert exactly, every
time, confirming the anchor itself is stable and correctly loaded, not a per-run fluke).

## 2. Results

### 2.1 Subject characteristics (diverse instance-space, not just repeats of subject2)

| Subject | mass (kg) | height (m) | sex | walking1 trial length | SO wall-clock |
|---|---:|---:|---|---:|---:|
| subject2 (reference) | 78.2 | 1.96 | M | 1.570s (132→**158** frames — the LONGEST of subjects 2–10) | (pre-existing) |
| subject3 | 63.5 | 1.69 | F | 1.310s (132 frames) | 15.0s |
| subject4 | 62.6 | 1.68 | F | 1.300s (131 frames) | 14.8s |

Measured live across all 9 other subjects with a walking1 trial (subject3–subject10): every one of them ends
**before** 1.50s (range 1.24–1.42s) — subject2's 1.57s is the longest trial in this cohort, not a typical one.
This matters for Sec.2.3/Sec.3 below.

### 2.2 Knee — clean, no truncation effects, tight cross-method agreement in all 3 subjects

| Subject | self-computed peak | JR peak | self-vs-JR agreement | peak time | **ratio vs OrthoLoad (258.22 %BW)** |
|---|---:|---:|---:|---:|---:|
| subject2 | 391.10 %BW | 391.11 %BW | 0.0039% | 0.51s | **1.515×** |
| subject3 | 513.75 %BW | 513.75 %BW | 0.0015% | 0.49s | **1.990×** |
| subject4 | 434.14 %BW | 434.15 %BW | 0.0013% | 0.46s | **1.681×** |

All 3 knee peaks land comfortably interior to their trial (well away from any edge), all 3 anatomical
crossing-muscle anchors are fully present (13/13 expected knee-crossers detected in every subject), and
self-computed vs. JointReaction — two structurally independent code paths — agree to **<0.005%** in every
single subject. This is the cleanest, most directly comparable result in this test: **knee over-prediction
generalizes in direction to both new subjects, and subject2 (1.515×) is the mildest of the three, not the
worst case** (range 1.515×–1.990×, subject2 at the low end).

### 2.3 Hip — direction still generalizes, but a real convention/truncation confound must be disclosed

| Subject | self-computed (interior peak) | JR (whole-trial peak) | naive agreement | self peak t | JR peak t | **ratio (self / JR)** |
|---|---:|---:|---:|---:|---:|---:|
| subject2 | 386.77 %BW | 387.04 %BW | 0.070% | 0.55s | 0.55s | **1.412× / 1.413×** |
| subject3 | 455.77 %BW | 477.93 %BW | 4.64% | 0.51s | **1.31s (= trial end)** | **1.664× / 1.745×** |
| subject4 | 418.94 %BW | 497.44 %BW | 15.78% | 0.17s | **1.30s (= trial end)** | **1.529× / 1.816×** |

The "self-vs-JR agreement" column looks alarming for subject3/4 (4.6%, 15.8%, vs. subject2's tight 0.07%) —
**forced to Orient rather than accept at face value**: `static_opt_knee.py`'s self-computed peak search
excludes edge frames (`Fbc_mag[e:-e]`, the same convention subject2's own cert already uses), while the
JointReaction extractor searches the **whole trial including the boundary**. For subject2's long (1.57s)
trial this never mattered — both methods found the same interior peak. For subject3/4's shorter trials, the
hip loading signal is **still climbing steeply, smoothly, and monotonically** (verified frame-by-frame, not
a numerical discontinuity — Sec.3) right up to the last captured frame. JR's un-excluded search lands on that
still-rising boundary value; the edge-excluded self-computed search instead reports an earlier, genuine local
peak. Recomputing the self-computed signal **at the same boundary instant** (still calling the identical,
unmodified geometry function, just not excluding that frame from the reported max) gives 488.66 %BW
(subject3) and 497.72 %BW (subject4) — matching JR's 477.93/497.44 %BW to **2.2%/0.06%** respectively. **The
two methods do agree tightly when compared at the same instant, in every subject** — the apparent
disagreement is a pre-existing peak-search convention mismatch in the shared cert code, exposed (not caused)
by a new, shorter-trial instance-space, not a genuine cross-method disagreement.

**Practical consequence**: subject3/4's true peak hip contact force was likely not fully captured (their
trials cut off while still ascending) — so the reported hip ratios above (1.53–1.82×) are, if anything,
**lower bounds**. This does not weaken the "hip over-prediction generalizes" finding; if anything it means
the true magnitude could be even larger than reported here.

### 2.4 Summary: generalizes in direction, not in exact magnitude

| | subject2 (reference) | subject3 | subject4 | range across all 3 |
|---|---:|---:|---:|---:|
| Knee ratio | 1.515× | 1.990× | 1.681× | **1.515× – 1.990×** |
| Hip ratio (self-computed convention) | 1.412× | 1.664× | 1.529× | **1.412× – 1.664×** |
| Hip ratio (JR convention) | 1.413× | 1.745× | 1.816× | **1.413× – 1.816×** |

n=3 is descriptive, not a population estimate — no formal statistics (mean/SD/CI) are claimed. What IS
supported: **every ratio, at every joint, in every subject, is comfortably above 1.4×** — nowhere near the 1.0×
boundary that would call the qualitative over-prediction finding into question — and subject2 sits at or near
the **low end** of the observed range at both joints, not a cherry-picked worst case.

## 3. Forced adversaries (pre-registered, not post-hoc)

**Pre-registered threshold for "generalizes"**: ratio > 1.0 (over-prediction persists in direction) at both
joints for both new subjects, AND self-vs-JR agreement stays tight (comparable to subject2's <0.1%) — the
second clause is the falsifier for "this is a genuine physiological/architectural tendency" vs. "this is a
per-subject pipeline defect."

- **Leaning-positive adversary** (the thing that would make "generalizes" a false positive): new-subject
  ratios landing near subject2's by a shared-instrument artifact rather than a real tendency. Forced by: (a)
  checking self-vs-JR agreement per subject — the knee stayed <0.005% in all 3, ruling out a single-code-path
  artifact; (b) using bodies that differ from subject2 in mass (−15 to −16 kg), height (−27 to −28 cm), and
  sex (F vs M) — the diverse instance-space this claim must survive. Result: adversary did NOT explain the
  finding away — ratios did not just "land near 1.5×," they landed **higher**, in a direction a shared-artifact
  story does not predict.
- **Leaning-negative adversary** (the thing that would make "subject-specific" a premature negative): the
  hip's apparent 4.6%/15.8% self-vs-JR disagreement, taken at face value, would look like a pipeline problem
  specific to these subjects. Forced via OODA rather than accepted as an "honest negative": Observed the
  JointReaction peak landing exactly at each trial's last frame (t=1.31s/1.30s); Oriented by pulling the raw
  per-frame force time series for both methods near the boundary (Sec.2.3) — found a smooth, monotonic,
  20+-frame ramp in BOTH methods' own underlying signal, not a discontinuity, and R_pure (the muscle-free
  Newton's-law term, driven only by measured GRF + kinematics) climbing the same way — a genuine kinematic/GRF
  event, not differentiation noise; Decided the real cause was the pre-existing edge-exclusion-vs-whole-trial
  convention mismatch in the reused code, not a new bug; Acted by recomputing self-computed's own signal AT the
  boundary (still the unmodified geometry function) and found <2.2% agreement with JR there. The adversary
  fell: this is a disclosed convention/data-length confound, not evidence against generalization.
- **A third, self-audited adversary while writing this doc**: could the shared library itself (OpenSim's
  Logger) have silently written to the read-only external drive? Checked directly (`find -newer`, not
  assumed clean) — **yes**: `opensim.log` (OpenSim's own duplicate of the info/warning lines already on
  stdout) was written under `subject3/.../Model/` and `subject4/.../Model/`, confirmed via mtime and content.
  Traced to a pre-existing OpenSim C++ Logger default (subject2's own `opensim.log` is an 11MB file
  accumulated across this project's entire history — every "unmodified" cert script that ever loaded
  subject2's model already did this; not something this new script introduced). Verified scope: **only**
  `opensim.log` changed under each subject's directory — zero bytes of any `.mot`/`.osim`/`.xml`/`.sto` data
  file were touched (checked via `find -newer` over the full subject3/subject4 trees). Mitigated going forward
  by adding `osim.Logger.removeFileSink()` to the new script (verified live: re-running subject3 after adding
  this reproduced identical numbers — 1.990×/1.664×, bit-for-bit — with the log file's mtime **unchanged**,
  confirming the fix works and does not alter any result). The two log files already written during this
  session's earlier runs were **left in place, not deleted** — deleting is itself a further write to the
  external drive that this task did not request; flagged here for the operator to act on if wanted.
  The identical OpenSim Logger behavior also appends to `opensim.log` at the **bodytwin repo root** — checked
  via `git diff --stat` (+6283 lines this session) and `git log` (last commit touching it predates this
  session, 2026-07-21 10:55:28) — confirming this is a **pre-existing, already-git-tracked, repo-wide
  append-log pattern** every OpenSim-based script in this project has always added to, not something new.
  Not staged, not committed (isolation respected). A concurrent, unrelated observation made while checking
  this: `git status` also shows several files this task never touched (`docs/MECHANISM_CLIMBING_SCENE.md`,
  `docs/MECHANISM_CORPUS_IK.md`, `docs/MECHANISM_MSK_ELASTIC_BAND.md`, `scripts/msk/pose_to_opensim_ik.py`,
  `scripts/msk/add_trunk_flexors_evidence.json`) as modified/untracked, with mtimes predating this session's
  first write — consistent with `COORDINATOR.md` §0's description of two concurrent mechanism instances sharing this
  repo, not an effect of this work (verified via mtime, not assumed).

## 4. Pipeline sanity across subjects (symmetric to subject2's own gates)

- **Subject-agnostic solve quality** (frames>100, zero NaN activations, zero out-of-bounds activations):
  **PASS** for both subject3 and subject4.
- **Raw `convergence_pass` gate (unmodified function)**: reads **False** for both new subjects — diagnosed,
  not silently worked around: `check_so_convergence_and_sanity`'s `covers_validated_peak_time_1p50s` sub-check
  is a subject2-specific literal (subject2's own validated peak occurs at t≈1.50s). Measured live: **every**
  subject 3–10's own trial ends before 1.50s (subject2's is the longest of the ten, Sec.2.1) — so this
  sub-check fails **by construction** for any other subject in this cohort, not because of a real convergence
  problem. Reported both raw (technically-FAIL) and subject-agnostic (PASS) readings rather than picking one.
- **Pelvis residual force** (RRA not run, same known gap as subject2 — `docs/MECHANISM_JOINT_FORCE_SCORECARD.md`
  §3 item 1): FAILS the <75N band in all 3 subjects (subject2 173.0N, subject3 567.0N, subject4 340.1N) — a
  consistent, pre-existing, disclosed limitation of this cert family, not a new subject-specific anomaly.
- **Anatomical crossing-muscle anchors**: ALL 13 expected knee-crossers and ALL 25 expected hip-crossers
  detected in both new subjects (0 missing); the `bfsh_r` negative control (must NOT cross the hip) correctly
  excluded in both — the sign-fixed crossing-muscle geometry transfers correctly across differently-scaled
  instances of the same model topology.
- **No muscle exceeded 1.5× its max isometric force; no joint reserve actuator "leaned hard"** in either new
  subject.

## 5. Verification (machine-checked, not narrated)

- Cert-script integrity: `md5sum` of `validate_joint_force.py` / `static_opt_knee.py` / `validate_hip_force.py`
  identical before and after every run this session (checked independently of the script's own self-report):
  `e2bb15958cfe8d72f44a3e0ac3d95950` / `d00f962c1a8873606874fffd65e828c8` / `a6514b9f081e0b4767fb4b0f94e60200`.
- Determinism: re-running subject3 after adding the logger-suppression fix (Sec.3) reproduced
  `knee_self_ratio=1.989588477922362` bit-for-bit against the original run.
- Read-only respect: `find <subject_dir> -newer cross_subject_validation.py` over the full subject3 and
  subject4 trees shows **only** `opensim.log` (Sec.3) — no data file was ever modified.
- All numbers in Sec.2 were pulled directly from the on-disk JSON
  (`data/msk_smoketest/cross_subject_validation/subject{3,4}_walking1_cross_subject_results.json`), not
  transcribed from console prose.

## 6. Honest gaps

1. **n=3 subjects total** (subject2 + 2 new) is a small, descriptive instance-space, not a validated
   distribution — sufficient to show the over-prediction *direction* is not a subject2 fluke, insufficient to
   claim a precise population magnitude or to characterize its shape (e.g., whether it correlates with body
   mass/height/sex — confounded here since both new subjects are female and both new subjects are lighter, so
   a size effect cannot be distinguished from a sex effect or from ordinary individual variation with only 3
   points).
2. **Right side only, walking1 only** — same scope caveat inherited from every cert in this family; left side
   and other trials (subject3/4 both also have walking2/walking3/walkingTS* trials available) not tested.
3. **Trial-length confound**: subject3/4's shorter captures (1.30–1.31s vs subject2's 1.57s) is itself a
   variable this comparison did not control for, and Sec.2.3 shows it materially affects the hip result. A
   cleaner future comparison would normalize by %gait-cycle rather than absolute seconds, or select/trim
   trials to matched durations.
4. **Root cause of the magnitude spread (1.41×–1.99×) across subjects is not identified.** The knee result is
   NOT truncation-affected (clean interior peaks in all 3 subjects) and still shows subject2 as the mildest
   case — so real inter-subject variation exists independent of the hip's truncation issue. Candidate,
   unverified explanations named but not tested here: muscle `Fmax` is not rescaled per-subject
   (`docs/MECHANISM_LAYER_COMPLETENESS.md`'s Muscles row: "confirmed bit-identical generic-vs-scaled... never
   subject-rescaled") — a lighter subject using the same absolute muscle strength has a different
   strength-to-bodyweight ratio, which could plausibly shift SO's recruitment solution; genuine
   anthropometric/gait-pattern differences; or some mix. Not disentangled here — would need a controlled
   Fmax-rescaling experiment across subjects to test directly.
5. **Two `opensim.log` files were written to the read-only external drive** during this session (Sec.3) — text
   log only, no data touched, mitigated for future runs, left in place rather than unilaterally deleted.
6. **Pre-existing cert limitations inherited unchanged**: no RRA (pelvis residual fails its band in all 3
   subjects), Static Optimization is activation-minimizing (not measured EMG — real co-contraction can exceed
   it), Savitzky-Golay differentiation is a design choice. None of this is new to this cross-subject test; all
   of it was already disclosed for subject2 and reproduces identically here.

## 7. Files

- `scripts/msk/cross_subject_validation.py` (new) — the orchestration + comparison script described above.
- `docs/MECHANISM_CROSS_SUBJECT.md` — this doc.
- `data/msk_smoketest/subject3_walking1/`, `data/msk_smoketest/subject4_walking1/` — new SO/JR outputs
  (mirrors subject2's own `data/msk_smoketest/subject2_walking1/` layout).
- `data/msk_smoketest/cross_subject_validation/subject3_walking1_cross_subject_results.json`,
  `.../subject4_walking1_cross_subject_results.json`, `.../cross_subject_summary.json` — machine-readable
  results, source for every number in Sec.2.
- **Not touched** (verified, Sec.5): `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_hip_force.py`.
- No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
