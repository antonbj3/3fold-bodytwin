# MECHANISM GLENOHUMERAL FORCE — increasing resolution on the shoulder now that arm muscles + scapula exist (2026-07-21)

Executes the operator's task: `docs/MECHANISM_SHOULDER_FORCE.md` was an earlier honest-negative
(0/80 arm muscles, no muscle-driven solve possible). The unified model now HAS arm muscles
(`docs/MECHANISM_ARM_MUSCLES.md`) and a real scapula/clavicle chain
(`docs/MECHANISM_SCAPULA_CLAVICLE.md`), merged into `data/msk_models/subject2_unified_v2.osim`
(`docs/MECHANISM_UNIFIED_V2.md`) — but neither the joint-force method nor a real gait trial had
ever been run on this exact substrate. This session: (1) FIRST enumerates whether subject2's
real `walking1` trial exercises the glenohumeral (GH) joint meaningfully, measured on real
mocap, not assumed; (2) re-runs the muscle-driven capacity pipeline at a functional
arm-elevation pose on the CURRENT (scapula-equipped) model; (3) verifies the OrthoLoad
Bergmann/Westerhoff PMIDs live; (4) forces a real adversary against a ~7x magnitude jump the
capacity re-run produced, before reporting anything. Every number below is machine-measured
this session (`scripts/msk/validate_glenohumeral_force_unified.py`, exit 0), isolation
respected (`.venv-msk` only, no git operations, only new files written).

## Headline result

| | | |
|---|---:|---|
| **Walking1 peak GH reaction (real mocap, pure kinematics)** | **6.72 %BW** | at t=0.71s |
| Walking1 ROM: peak\|arm_add_r\| (abduction-relevant DOF) | **19.60°** | vs. the anchor's 90° pose |
| Data-gate verdict (pre-registered: peak<20%BW AND ROM<45°) | **PASS — CONFIRMED DATA-GATE** | walking under-loads the GH joint |
| Task-matched functional pose found (FK-verified, this session) | arm_add_r = **−60.45°** | gives 90.0° TOTAL humerothoracic elevation |
| Capacity at that pose, self-computed (muscle-corrected) | 132.9 %BW | **NOT TRUSTED — see §4** |
| Capacity, official `JointReaction` referee | **NaN/unusable** (69–76% NaN rows, some ±inf) | referee broken on this topology |
| OrthoLoad shoulder in-vivo anchor (median, n=23, live re-derived this session) | **72.2 %BW** | Bergmann/Westerhoff, PMIDs verified §5 |
| Most recent TRUSTWORTHY capacity estimate (pre-scapula model) | **19.065/19.362 %BW** (ratio 0.264/0.268) | `docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`, unchanged by this session |

**CONFIDENCE TIER: walking claim = high-confidence DIAGNOSED DATA-GATE (in-vivo comparison not
applicable — the task itself under-loads the joint). Capacity claim = DIAGNOSED-GAP, not
in-vivo-anchored** — the pose is correctly task-matched (a genuine improvement, §3), but the
resulting force MAGNITUDE is not trustworthy at this model resolution for a separate,
mechanistically-diagnosed reason (§4), and the independent referee that would normally certify
it is broken. The pre-scapula corrected number remains the best trustworthy capacity estimate
this twin currently has.

## 0. Model substrate — structural facts confirmed live, not assumed

`data/msk_models/subject2_unified_v2.osim` (34 bodies, 34 joints, 240 muscles, 79.0341 kg).
Diagnostic queried directly via the OpenSim API this session:

- The joint named **`acromial_r`** now sits **`scapula_r` → `humerus_r`** (a `CustomJoint`) — its
  parent frame was re-anchored from torso onto scapula_r by `add_scapula_clavicle.py`. This is a
  genuine resolution IMPROVEMENT over the pre-scapula model's own `acromial_r` (which was a
  combined torso→humerus proxy): the free-body cut at this joint now isolates the TRUE
  glenohumeral joint, a strictly better match to what OrthoLoad's prosthesis-mounted transducer
  actually measures.
- An enforced `CoordinateCouplerConstraint` named **`scapulohumeral_rhythm_r`** ties
  `scapula_upward_rot_r` (dependent) to `arm_add_r` (independent), confirmed `isEnforced=True`
  live. Consequence: `arm_add_r` is no longer the same physical angle it was pre-scapula — it is
  now GH-relative-to-SCAPULA, not GH-relative-to-torso (§3).
- BFS descendants of `acromial_r` (the free body distal to the GH cut): `humerus_r, ulna_r,
  radius_r, hand_r, wrist_int_r, index_{proximal,medial,distal}_r, middle_{proximal,medial,distal}_r`
  — chain mass **3.904 kg = 4.940 %BW-mass-equivalent** (of 79.0341 kg total).
  BW = 775.06 N (live-measured, this model's own mass — NOT the older 78.2 kg used pre-merge).
- **22 muscles cross the `acromial_r` cut** at rest pose — an EXACT match to the 22 GH-crossing
  muscles identified in the pre-scapula graft (`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`):
  `BIC_brevis, BIC_long, Coracobrachialis, DeltoideusClavicle_A, DeltoideusScapula_{M,P},
  Infraspinatus_{I,S}, LatissimusDorsi_{I,M,S}, PectoralisMajorClavicle_S,
  PectoralisMajorThorax_{I,M}, Subscapularis_{I,M,S}, Supraspinatus_{A,P}, TRIlong, TeresMajor,
  TeresMinor`. None of the 11 scapular-stabilizer muscles re-routed onto `scapula_r`/`clavicle_r`
  cross this specific cut (they attach thorax↔scapula, both proximal to it) — expected, not a gap.

## 1. Corpus enumeration — does ANY real trial in subject2's corpus load the shoulder?

Every real IK trial for subject2 was enumerated live (not just `walking1`):
`DJ1-3, DJAsym1/4/5, squats1, squatsAsym1, STS1, STSweakLegs1, walking1-3, walkingTS1/2/4`.
**Zero** trials target arm elevation/abduction/overhead reach — this is a drop-jump/squat/
sit-to-stand/gait lab corpus, confirming and generalizing the prior honest gap
(`docs/MECHANISM_SHOULDER_FORCE.md` §2) beyond `walking1` alone: no real motion-capture trial
anywhere in this subject's corpus was designed to load the GH joint.

## 2. PART A — walking1 data-gate (measured on REAL mocap, not assumed)

Real IK (`data/msk_smoketest/subject2_unified_v2_walking1/walking1_unified_v2_smoketest.mot`,
158 frames, 100 Hz, already computed by a prior merge's own smoke test — reused, not
regenerated) → full 51-coordinate FK per frame → pure-kinematic Newton's-law GH reaction (no
GRF term needed — the arm chain never touches the ground, same simplification
`docs/MECHANISM_SHOULDER_FORCE.md` used).

| coordinate | min (deg) | max (deg) | **ROM (deg)** | peak\|.\| (deg) |
|---|---:|---:|---:|---:|
| arm_flex_r | −15.33 | −4.06 | 11.26 | 15.33 |
| **arm_add_r** (abduction-relevant) | −19.60 | −18.00 | **1.60** | 19.60 |
| arm_rot_r | 34.15 | 43.54 | 9.40 | 43.54 |
| scapula_upward_rot_r | 9.00 | 9.80 | 0.80 | 9.80 |

(ROM — max minus min — is the primary, offset-invariant readout: `docs/MECHANISM_UNIFIED_V2.md`
§4c already disclosed a ~55° absolute registration shift for `arm_rot_r` vs. the pre-merge
reference model, from the scapula's own 19.77° rest-tilt; ROM is robust to that constant
offset, absolute peak values are reported for context only.)

**Reaction-force sensitivity sweep** (Savitzky-Golay window, right arm):

| window (ms) | 70 | 90 | **110 (primary)** | 130 | 170 |
|---|---:|---:|---:|---:|---:|
| peak (%BW) | 6.785 | 6.675 | **6.717** | 6.640 | 6.564 |

Range across the first 4 windows = 0.145 %BW → **stable**. Static (gravity-only) floor =
4.940 %BW (the chain's own weight); the measured peak is 1.36× this floor — a small, expected
inertial arm-swing contribution, not a load event.

**Gates (pre-registered before measuring):**

| gate | threshold | measured | verdict |
|---|---|---:|---|
| ROM: peak\|arm_add_r\| during walking | < 45° (half the anchor's 90° pose) | 19.60° | **PASS** |
| Load: peak reaction during walking | < 20 %BW (well under 72–110 %BW anchor regime) | 6.717 %BW | **PASS** |
| Window-stability | range < 2.0 pct-BW | 0.145 | **PASS** |

**Verdict: DATA-GATE CONFIRMED, not assumed.** Real gait loads the GH joint to ~6.7 %BW, an
order of magnitude below the ~72–110 %BW in-vivo anchor regime for a 90° abduction/elevation
task. Walking is the wrong task to validate glenohumeral capacity against this anchor — a
measured, falsifiable finding (the gate could have failed; it did not, by a wide margin).

## 3. PART B0 — task-matching the functional pose on the CURRENT (scapula-equipped) model

Because `arm_add_r` is now GH-relative-to-scapula rather than GH-relative-to-torso (§0), simply
reusing the pre-scapula 0°→90° `arm_add_r` sweep would no longer represent a 90°
TOTAL-elevation pose (which is what OrthoLoad's "90 deg abduction" task label — an externally
observable, clinically-measured arm-vs-trunk angle — actually refers to). A fresh FK sweep was
run live on `subject2_unified_v2.osim` to re-verify and re-derive this:

| arm_add_r | scapula_upward_rot_r | **total elevation (independent FK)** |
|---:|---:|---:|
| 0° | 0.000° | 0.000° |
| −10° | 5.000° | 14.879° |
| −30° | 15.000° | 44.653° |
| −50° | 25.000° | 74.437° |
| −70° | 35.000° | 104.213° |
| −90° | 45.000° | 133.943° |

The independent total-elevation FK numbers reproduce the standalone `add_scapula_clavicle.py`
build's own values to 3 decimal places (e.g. 133.943° at −90°) — confirming the coupler and its
geometry survived the two additional forks (hip/ankle ligaments, full hand) merged on top of it
since. (The scapular:GH RATIO computed as `scapula_upward_rot_r / -arm_add_r` reads exactly
0.5000 at every point — this is the SAME already-adjudicated circularity
`docs/MECHANISM_SCAPULA_CLAVICLE.md`'s own correction banner flagged: it is reading the
constraint's own hard-coded slope back out, not independent evidence. Only the total-elevation
FK column above is the non-circular, independent check, and it is what is used below.)

**Interpolating to find the arm_add_r that gives 90.0° TOTAL elevation: −60.45°** — a NEW
trajectory was generated (the original synthetic 0°→90° `arm_add_r` `.mot` file rescaled by
0.6717×, preserving its exact smootherstep shape/timing) so the capacity test below is properly
task-matched to the anchor's actual stated pose, not the ~134°-total pose a naive reuse would
have produced.

**A real bug was caught and fixed in this session's own script before trusting this number**: a
first version of the interpolation reversed both the x/y arrays under a wrong assumption about
which one needed increasing order, silently returning 0.0° (a degenerate no-motion target)
instead of −60.45°. Caught via a standalone known-answer check
(`np.interp(90,[0,14.879,...,133.943],[0,-10,...,-90])`) before the pipeline was allowed to
build the rescaled trajectory from it — forced per this repo's own symmetric-QC discipline, not
discovered after the fact.

## 4. PART B1/B2 — capacity re-run, and why the headline magnitude is NOT trusted

Locked-coordinate probe built from the live model's OWN `ConstraintSet` (not hardcoded): 47
coordinates locked, 4 left free (`arm_add_r` plus the 3 constraint-DEPENDENT coordinates —
`scapula_upward_rot_r`, `knee_angle_{r,l}_beta` — found by querying every
`CoordinateCouplerConstraint` live). Two legs run: **B1** reuses the OLD 0°→90° `arm_add_r`
trajectory verbatim (now representing ~134° total elevation — a reproducibility check, NOT
task-matched); **B2** uses the NEW task-matched trajectory (~90° total elevation, §3).

| | B1 (old traj., ~134° total) | B2 (task-matched, ~90° total) |
|---|---:|---:|
| R_old (pure kinematic) peak | 5.766 %BW | 5.346 %BW |
| Real-muscle activation range (post gate-fix, §4a) | [0.0000, 0.5641] | [0.0100, 0.3522] |
| Reserve actuator max usage (`shoulder_add_r`) | 0.361 (36%) | 0.160 (16%) |
| **F_new (muscle-corrected) peak** | **138.73 %BW** (t=1.50s) | **132.88 %BW** (t=2.31s) |
| Official `JointReaction` referee | **NaN** (181/261 rows, 69.3%) | **NaN** (197/261 rows, 75.5%) |
| ratio vs. OrthoLoad anchor (72.2%) | 1.92 | 1.84 |

### 4a. A checker bug found and fixed first (before the magnitude was even inspected)

The reused `validate_shoulder_force_with_muscles.check_so_convergence`'s
`IDEAL_ACTUATOR_COORD_NAMES` tuple lists only 5 RIGHT-side reserve-actuator names — written for
the older, arm-only pre-scapula graft model. This full-body unified model has **13**
`CoordinateActuator`s (adds 5 left-side + 3 lumbar reserves). The gate FAILED on first run
because it misclassified reserve-actuator CONTROLS (legitimately signed, e.g.
`shoulder_add_r=-0.361`) as if they were `[0,1]`-bounded muscle activations — direct inspection
showed **zero** real muscles among the offending values; all 8 were `CoordinateActuator` names.
**Fixed** by querying the model's own `ForceSet` live for actual `CoordinateActuator` names
instead of extending the hardcoded tuple blindly. After the fix, both legs show clean,
well-bounded real-muscle activations (no saturation, no NaN) — this fix is necessary and
correct, but does **not** explain the magnitude jump below.

### 4b. The real diagnosis: a single muscle's moment arm is near-degenerate

F_new (132–139 %BW) is **~7× the pre-scapula corrected capacity number** (19.065/19.362 %BW).
Rather than report this at face value, the per-muscle force decomposition at the B2 peak frame
was inspected directly:

| muscle | tension (N) | % of \|sum\| |
|---|---:|---:|
| **DeltoideusScapula_M** | **914.73** | **74.0%** |
| DeltoideusClavicle_A | 146.65 | 11.9% |
| Infraspinatus_S | 55.70 | 4.5% |
| (19 others, each ≤17.6N) | ≤86.2 combined | 9.6% |

`DeltoideusScapula_M` alone contributes 914.73 N = **118 %BW in one muscle**. Its own parameters
and geometry, measured live on this exact model:

- `max_isometric_force` = **2597.8 N** (large for a single deltoid head vs. typical Delp-lineage
  values of ~500–1200 N; copied unscaled from the donor per this graft's own disclosed
  convention, `docs/MECHANISM_MUSCLE_AUDIT.md` §3 — a pre-existing property, not new this session).
- Moment arm about `arm_add_r`, swept live: **+5.98mm (0°) → −2.77mm (−10°) → −13.65mm (−30°) →
  −13.36mm (−50°) → −8.50mm (−60.45°, the B2 target) → −1.71mm (−70°) → +15.99mm (−90°)** — TWO
  sign flips (crossing near-zero between 0/−10° and again between −70/−90°).

**Mechanism**: near a sign-flip, this muscle's leverage collapses toward zero, forcing Static
Optimization to recruit large tension through a near-useless lever to satisfy the net-torque
requirement; with an oversized `Fmax` available, it can do this at only moderate activation
(~0.35), producing a non-physiological single-muscle compressive contribution. **This is not a
new mystery** — it is the SAME muscle and the SAME qualitative signature already disclosed in
`docs/MECHANISM_ARM_MUSCLES.md` §6 for the pre-scapula model ("DeltoideusScapula_M ... its moment
arm genuinely FLIPS sign between poses ... the exact, well-known signature of an un-wrapped
muscle path"), now shown — by live measurement on the scapula-equipped model — to have **gained
a second zero-crossing**, consistent with two already-disclosed, not-yet-fixed gaps: no
wrapping surfaces ported (`docs/MECHANISM_ARM_MUSCLES.md` §8), and this muscle's proximal
attachment still collapsed onto "torso" rather than re-registered onto the new `scapula_r`
(`docs/MECHANISM_SCAPULA_CLAVICLE.md` §8 item 4, explicitly flagged there as "an explicit,
disclosed opportunity NOT taken").

B1's F_new additionally shows genuine numerical jaggedness (frame-to-frame jumps up to 20 %BW
per 0.01s clustered exactly at t=1.59–1.63s and t=2.27–2.29s) — consistent with the sweep
passing directly through this muscle's zero-crossing region. B2's F_new is comparatively smooth
(largest jump 2.1 %BW, stable plateau at the trajectory's static hold) — smoothness alone does
not certify the magnitude as physiological; §4b shows it remains dominated by the same
oversized-Fmax/small-moment-arm muscle, just evaluated away from its exact zero-crossing.

### 4c. The independent referee is itself broken on this topology

The official `opensim.JointReaction` analysis — the bug-immune referee that gave 1.54%
agreement with the self-computed number on the pre-scapula model — returns **NaN in 69.3%
(B1) / 75.5% (B2) of frames**, and at least one **infinite** value among the nominally "valid"
remaining rows in both legs. This is a clean, decisive, machine-checked failure (not ambiguous):
the independent over-determination check that would normally validate or refute the
self-computed number is unusable here, most plausibly because the locked probe leaves 3
constraint-dependent coordinates unlocked-but-constrained, and `JointReaction`'s own per-step
state reconstruction from a 2-column coordinates file appears fragile under that configuration.
Not chased further this session (reported as a diagnosed limitation, not silently worked
around).

**Verdict: the B1/B2 capacity magnitudes are measured but NOT TRUSTED.** Two independent,
compounding, diagnosed reasons: (a) the dominant muscle's Fmax+moment-arm combination is a
mechanistically-confirmed continuation of an already-disclosed model gap, and (b) the
cross-check referee is numerically broken on this topology. The pre-scapula corrected number
(19.065 %BW self-computed / 19.362 %BW JointReaction-verified, 1.54% agreement,
`docs/MECHANISM_SHOULDER_MUSCLES_CORRECTION.md`) remains this twin's most recent TRUSTWORTHY
glenohumeral capacity estimate — unchanged and not superseded by this session.

## 5. External anchor — OrthoLoad PMIDs verified live (NCBI E-utilities, not recalled)

An initial recalled guess for the Bergmann citation's PMID (17188287) was **wrong** (an
unrelated organic-chemistry paper) — caught by live verification via `curl` against
`eutils.ncbi.nlm.nih.gov`, not trusted from memory:

- **Bergmann G, Graichen F, Bender A, Kääb M, Rohlmann A, Westerhoff P (2007).** "In vivo
  glenohumeral contact forces—measurements in the first patient 7 months postoperatively." *J
  Biomech* 40(10):2139-49. doi:10.1016/j.jbiomech.2006.10.037. **PMID 17169364.** (n=1
  hemiarthroplasty patient; abstract: contact force stayed <100%BW for most ADL, up to 130%BW
  near ROM limits/resistance, ~150%BW at maximum resisted effort.)
- **Westerhoff P, Graichen F, Bender A, Halder A, Beier A, Rohlmann A, Bergmann G (2009).** "In
  vivo measurement of shoulder joint loads during activities of daily living." *J Biomech*
  42(12):1840-9. doi:10.1016/j.jbiomech.2009.05.035. **PMID 19643418.** (n=4 patients, 8
  weighted ADL tasks; e.g. lifting a 1.5kg coffee pot = 105.0%BW mean (range 90-124.6); lifting
  2kg to head height = 98.3%BW (93-103.6); highest forces at ≥90° abduction/elevation WITH a
  hand load.)

**Cross-check vs. this repo's own AKF-derived anchor**: this session independently
re-*derived* (not just re-cited) the repo's shoulder AKF anchor live —
`select_orthoload_shoulder_trials()` re-run fresh this session: **n=23 unloaded 90°
abduction/elevation trials, 6 distinct subjects (s1r, s2r, s3l, s4r, s5r, s8r), min=36.28
median=72.21 mean=78.16 max=128.06 %BW** — matching the previously-published 72.2%BW to
rounding. This unloaded-task median sits below the two papers' own WEIGHTED-task figures
(91–131%BW), the expected direction (unloaded < weighted), and within Bergmann 2007's own
stated "up to 130%BW near ROM limits" ceiling — consistent, no red flag; both draw from the same
underlying OrthoLoad telemetry program. **Population caveat carried forward, unchanged**: both
papers' patients are instrumented hemiarthroplasty/reconstruction patients; the twin (subject2)
is a healthy OpenCap participant.

## 6. Falsifier verdict

*"Does the modeled glenohumeral contact match the OrthoLoad shoulder in-vivo range for a
matched task, or is the walking task load-inadequate (diagnosed-gap)?"*

- **Walking task: DIAGNOSED-GAP, confirmed not assumed.** Peak real-gait GH reaction = 6.72%BW,
  ~11x below the 72.2%BW anchor and ~1 order of magnitude below the general 72-130%BW in-vivo
  regime; ROM of the abduction-relevant coordinate is 1.60° (peak 19.60°) vs. the anchor's 90°
  pose. Walking does not exercise the shoulder meaningfully in this model/trial — the
  pre-registered gates (ROM<45°, load<20%BW) both PASS by a wide margin.
- **Functional-pose capacity comparison: a SEPARATE diagnosed-gap, for a different reason.** The
  pose itself is now correctly task-matched (90.0° total humerothoracic elevation, FK-verified
  live, §3) — a genuine resolution improvement over naively reusing the pre-scapula trajectory.
  But the resulting force magnitude (132-139%BW) is not a trustworthy match-or-mismatch signal
  against the 72.2%BW anchor: it is dominated by a single muscle's diagnosed
  Fmax/moment-arm-degeneracy artifact (§4b) and lacks a working independent referee (§4c). This
  is a MODEL-fidelity gap, not a task-mismatch gap — orthogonal to the walking finding above.

## 7. Honest caveats carried forward (unchanged from prior certs) + new this session

1. No wrapping surfaces ported for any of the 25 arm muscles (unchanged, `docs/MECHANISM_ARM_MUSCLES.md` §8) — now shown to matter MORE at the scapula-equipped topology (§4b), not less.
2. The 25 original arm muscles' proximal attachments remain collapsed onto "torso" rather than re-registered onto the new `scapula_r`/`clavicle_r` bodies (disclosed, not fixed, `docs/MECHANISM_SCAPULA_CLAVICLE.md` §8 item 4) — directly implicated in §4b's diagnosis.
3. Acromioclavicular joint is not a real moving joint (clavicle locked at rest pose) — unchanged, does not affect this session's GH-only free-body cut.
4. Only one scapular DOF (`scapula_upward_rot_r`) is kinematically driven; the 2:1 rhythm ratio is a literature constant, not independently re-derived from this twin's own data — unchanged.
5. Single free-body cut, right arm only, one synthetic (task-matched) trajectory for the capacity leg — same scope caveat every joint-force cert in this family carries.
6. `JointReaction`'s fragility on the multi-constraint locked-probe topology (§4c) is NEW this session and not yet root-caused beyond the plausible mechanism named — a genuine open item for any future re-attempt at this exact leg.
7. Different populations: OrthoLoad's shoulder subjects are hemiarthroplasty/reconstruction patients; subject2 is a healthy participant — unchanged.
8. Total-elevation FK measurement (§3) uses a 2-point body-origin vector proxy for "humerus long axis," not a formally-defined anatomical axis — adequate for the 6-point sweep and interpolation used here, not validated as a general-purpose measurement.

## Files

- `scripts/msk/validate_glenohumeral_force_unified.py` — the full pipeline (self-contained,
  re-runnable; imports `validate_joint_force.py` and `validate_shoulder_force_with_muscles.py`
  for proven BFS/Savitzky-Golay/crossing-muscle code, not reimplemented; includes the
  live-queried ideal-actuator-name fix, §4a).
- `data/msk_models/subject2_unified_v2_shoulderSO_probe.osim` — the new locked-coordinate probe
  (47 locked / 4 free: `arm_add_r` + 3 live-queried constraint-dependent coordinates).
- `data/msk_smoketest/subject2_unified_v2_glenohumeral/glenohumeral_force_unified_results.json`
  — every number in §2-4's headline tables, machine-written by the pipeline.
- `data/msk_smoketest/subject2_unified_v2_glenohumeral/glenohumeral_capacity_diagnosis_addendum.json`
  — the forced-adversary diagnosis (§4's per-muscle decomposition, moment-arm sweep, JointReaction
  NaN/inf statistics, PMID verification detail, corpus enumeration), machine-written this session.
- `data/msk_smoketest/subject2_unified_v2_glenohumeral/synthetic_taskmatched_abduction.mot` — the
  new, rescaled (0.6717×) task-matched trajectory (90.0° total elevation).
- `data/msk_smoketest/subject2_unified_v2_glenohumeral/so_output_b{1,2}_*/`,
  `jr_output_b{1,2}_*/` — raw Static Optimization and JointReaction AnalyzeTool outputs for both
  legs.
- Reused, unchanged: `scripts/msk/validate_joint_force.py`,
  `scripts/msk/validate_shoulder_force_with_muscles.py`,
  `data/msk_smoketest/subject2_unified_v2_walking1/walking1_unified_v2_smoketest.mot` (real
  walking1 IK on the current unified model, from `docs/MECHANISM_UNIFIED_V2.md`'s own smoke test),
  `data/external/orthoload/shoulder/database_api/akf/*.akf` (OrthoLoad raw data, read in place).
