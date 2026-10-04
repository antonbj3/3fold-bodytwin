# MECHANISM MZ PELVIS-MOMENT FLOOR — isolating the source (2026-07-21)

Executes the operator's task: `docs/MECHANISM_RRA_TASK_GAINS.md` found a ~20-30 N·m pelvis-moment
residual, labeled "MZ" in `run_rra.py`'s own residual-actuator convention, resists ALL task-gain
tuning (up to 50× weight ratio, `docs/MECHANISM_RRA_TASK_GAINS.md`) and actuator boosting (150 vs 300
N·m identical, `docs/MECHANISM_ANKLE_RESERVE_FIX.md`), and was left unisolated among 3 candidates:
(1) GRF-model consistency, (2) transverse-plane observability, (3) intrinsic kinematic-GRF
inconsistency. Every number below is machine-measured this session
(`scripts/msk/mz_pelvis_moment.py`, exit 0), cross-checked against the pre-existing OpenCap-computed
Inverse Dynamics output already on disk (exact, bit-for-bit reproduction, §2) and an independent,
from-scratch Newton's-law computation (§4, a different code path than any OpenSim tool). Isolation
respected: `.venv-msk` only, original model/IK/GRF read in place and never written (verified: this
session opened them read-only via `osim.Model`/`parse_mot`, and separately reused, never modified,
existing RRA outputs from `ankle_reserve_fix/` and `rra_task_gain_tune/`), all outputs under a NEW
tree, `data/msk_smoketest/subject2_walking1/mz_pelvis_moment/`, no git operations.

## Headline result

**The operator's own framing needed a correction before anything else could be isolated (§1): "MZ"
(`run_rra.py`'s residual-actuator label) is the SAGITTAL plane (pelvis_tilt-like), not the transverse
plane. The genuinely transverse-plane residual actuator is "MY" — and MY was already the SMALLEST,
best-behaved axis throughout every configuration tested in the upstream docs.** With that correction
made, the floor traces cleanly to **candidate (3): an intrinsic kinematic-vs-GRF inconsistency in the
SAGITTAL plane**, present with ZERO optimizer in the loop and robust to filter-cutoff choice.
Candidates (1) and (2), tested directly on their own terms, are both **REJECTED**.

| candidate | pre-registered test | result |
|---|---|---|
| (1) GRF-model consistency | CoP stays inside the model's own tracked foot span during stance; free moment inside a friction-plausibility bound | **REJECTED** — CoP violations negligible (≤1.7 cm, ≤0.9% of stance frames, both feet); free moment (≤5.08 N·m) is ~17% of a generous 29 N·m friction bound |
| (2) transverse-plane observability | the stubborn axis IS transverse; transverse-plane DOFs show anomalous noise vs an external reliability anchor | **REJECTED** — the stubborn axis is SAGITTAL, not transverse (§1); transverse (MY/pelvis_rotation) is the smallest, best-behaved component in BOTH the raw data and every RRA configuration; the external anchor (McGinley et al. 2009) says sagittal is the MOST reliable plane, the opposite of what candidate 2 needs |
| (3) intrinsic kinematic-GRF inconsistency | same axis dominates a ZERO-optimizer plain-ID computation, robustly across filter cutoffs, and RRA's own residual converges to match it as tracking tightens | **SUPPORTED** — sagittal dominates plain ID (zero RRA/CMC/task-gain) at every tested lowpass cutoff (4-20 Hz, and unfiltered); RRA's own MX/MY/MZ correlate with the ID-based sagittal/frontal/transverse decomposition along the predicted diagonal, cleanly (r=0.73-0.97) once tracking is tight, messily (r=-0.51 to 0.37) when it is loose — explained, not just observed (§5) |

**Fixable or intrinsic?** Symmetric answer, not forced either way: the floor is **NOT fixable by any
RRA/CMC parameter (task gains, actuator magnitude)** — that is now established across three sessions.
But "intrinsic to this kinematics-GRF pairing" is NOT the same claim as "intrinsic to the physical
world" — three distinct root-mechanism sub-hypotheses (mass-property/scaling error, sagittal marker/
soft-tissue error, rigid-body model-fidelity limit) remain open and were NOT distinguished this
session (§7, honest gap). What IS answered: it is not a marker-observability problem in the classic
sense (§1, §6), and it is not a GRF/CoP data problem (§5).

## 1. Axis-convention correction — forced BEFORE anything else, not assumed

The operator's task brief calls the stubborn residual "MZ-axis (transverse-plane)". This is checked,
not accepted, against **two independent parses that agree exactly**:

- **Live OpenSim API**: `model.getGravity() = (0, -9.80665, 0)` → global **Y is vertical**. The
  `ground_pelvis` `CustomJoint`'s own `SpatialTransform` (`CustomJoint.safeDownCast(...).getSpatialTransform().getTransformAxis(i)`)
  reports: `pelvis_tilt` → axis `(0,0,1)` = global **Z**; `pelvis_list` → axis `(1,0,0)` = global **X**;
  `pelvis_rotation` → axis `(0,1,0)` = global **Y**.
- **Independent raw-XML regex parse** of the same `.osim` file's `ground_pelvis` block (a different
  code path than the SWIG API call above): identical result, `two_independent_parses_agree: true`.
- **Cross-check via the GRF file's own torque columns** (`ForceData/walking1_forces.mot`):
  `R_ground_torque_x` and `R_ground_torque_z` are **exactly 0.0** for the entire trial (both feet,
  machine precision, not "small") — the standard force-plate convention (CoP is *defined* as the
  point where the horizontal-plane moment vanishes). Only `torque_y` (the vertical/transverse free
  moment) is nonzero (≤5.08 N·m R, ≤4.71 N·m L) — independently confirming Y is vertical/transverse.

Since a rotation about a global axis produces motion in the plane PERPENDICULAR to that axis:

| `run_rra.py` residual-actuator label | global axis | pelvis coordinate | **anatomical plane** |
|---|---|---|---|
| MX | (1,0,0) | pelvis_list | **frontal** (obliquity/list) |
| MY | (0,1,0) | pelvis_rotation | **transverse** (axial, internal/external rotation) |
| MZ | (0,0,1) | pelvis_tilt | **sagittal** (anterior/posterior tilt) |

**So "MZ" is SAGITTAL, not transverse.** The genuinely transverse-plane residual actuator is MY —
and the existing evidence on disk (`rra_task_gain_tune.py`'s own docstring) already recorded MY as
the smallest RMS component throughout (MX=46.5 > MZ=29.4 > **MY=7.6** N·m at baseline). This
correction changes which mechanism/literature is even relevant to candidate (2) — see §6.

## 2. The zero-optimizer floor (candidate 3, central test)

A completely fresh `opensim.InverseDynamicsTool` run — muscles excluded, **zero RRA, zero CMC, zero
task-gain optimizer anywhere in the computation** — on the ORIGINAL, untouched `walking1.mot` (IK) +
`walking1_forces.mot` (GRF). This is pure Newton-Euler: `pelvis_tilt/list/rotation_moment` are the
generalized forces literally required to reconcile the given (real, marker-based) kinematics with
the given (real, measured) GRF at the pelvis's own free, unactuated 6 DOFs.

**Reproducibility, checked not assumed**: this session's path-patched run (reusing
`static_opt_knee.sub_tag`, never `run_rra.py`'s or `validate_joint_force.py`'s own unsafe non-
anchored variant) reproduces the pre-existing OpenCap-computed `OpenSimData/Mocap/ID/walking1.sto`
**exactly** — `max_abs_diff = 0.0` on all 6 pelvis columns (tilt/list/rotation moment, tx/ty/tz
force). The known OpenSim 4.6 `InverseDynamicsTool` bug (`docs/MECHANISM_MSK_ELASTIC_BAND.md` §4,
~1000-1800× corruption confirmed CONFINED to `knee_angle_r/l`/`hip_flexion_r/l`) does not touch any
column read here — independently re-confirmed there via a from-scratch virtual-work method
(<0.03% agreement) for all 6 `ground_pelvis` coordinates, and this session's own sanity-bound guard
(<500 N·m) never fires.

**An edge-filter artifact was caught and corrected, not glossed over**: the raw full-trial peak for
`pelvis_tilt_moment` (57.86 N·m at t=0.000s) is a boundary transient — excluding the first/last 5%
of the trial (matching this repo's existing edge-artifact-exclusion convention,
`rra_task_gain_tune.py`), the interior peak drops to 17.71 N·m (3.3× smaller) and interior RMS drops
from 10.85 to 8.12 N·m. `pelvis_list`/`pelvis_rotation` show NO such edge effect (their peaks sit at
t=0.19s/0.48s, unchanged by edge exclusion) — this artifact is specific to `pelvis_tilt`, plausibly
because it alone carries a large, dynamically-sensitive baseline level near the trial start (a
representative static pose already showed pelvis_tilt has an unusually large moment sensitivity in
this exact model, `docs/MECHANISM_MSK_ELASTIC_BAND.md` §6, -605 N·m). All numbers below use the
edge-excluded (interior) values.

**Result — even with ZERO optimizer, the ranking already matches the direction RRA's own tuning was
pushing toward:**

| plane | coordinate | interior peak (N·m) | interior RMS (N·m) |
|---|---:|---:|---:|
| **sagittal** | pelvis_tilt_moment | 17.71 | **8.12** |
| frontal | pelvis_list_moment | 13.63 | 6.09 |
| transverse | pelvis_rotation_moment | 11.48 | 5.09 |

Sagittal dominates by RMS; transverse is smallest — with **no RRA, no CMC, no task-gain lever
touching this computation at all.**

## 3. Robustness to filter-cutoff choice (differentiation-noise vs genuine-inconsistency)

If the sagittal dominance were an artifact of one arbitrary lowpass-filter choice (double-
differentiation noise amplified a specific way), the ranking should be unstable across cutoffs. It
is not — swept over a 5× range plus a fully unfiltered case:

| lowpass cutoff | sagittal (tilt) interior RMS | frontal (list) | transverse (rotation) |
|---:|---:|---:|---:|
| 4 Hz | 6.40 | 6.06 | 5.14 |
| 6 Hz (primary, matches this pipeline's own RRA/SO/JR convention) | 8.12 | 6.09 | 5.09 |
| 10 Hz | 10.20 | 6.25 | 5.35 |
| 20 Hz | 12.29 | 7.51 | 5.88 |
| unfiltered (-1, no lowpass) | 60.60 | 36.62 | 18.71 |

Absolute magnitudes grow with cutoff (expected: less smoothing passes more high-frequency
differentiation noise) but **transverse stays smallest at every single cutoff, and sagittal is
largest or tied-largest at every cutoff** (nearly tied with frontal only at 4 Hz: 6.40 vs 6.06). Even
fully unfiltered — a totally different noise regime — the same ordering holds. This is the falsifier
this test was pre-registered against (a cutoff-dependent ranking flip would have rejected "intrinsic"
in favor of "differentiation-noise artifact"); it did not flip.

## 4. Independent cross-check: from-scratch Newton's-law computation (decorrelated code path)

A SEPARATE, from-scratch computation — reusing only primitives already validated elsewhere in this
repo (`vjf.savgol_smooth_and_derivs`, the proven Savitzky-Golay differentiator; the same COM-from-
kinematics loop pattern as `run_rra.compute_peak_net_grf_and_com_height`) — computes whole-body COM
acceleration directly from the RAW Cartesian COM trajectory, then applies Newton's 2nd law:
`F_residual(t) = M_total·a_COM(t) − F_gravity − F_GRF_R(t) − F_GRF_L(t)`. This never calls
`opensim.InverseDynamicsTool` at all — a fully decorrelated angle on the SAME translational-residual
question:

| axis | ID tool (Stage 2) | independent Newton (this check) | agreement |
|---|---:|---:|---|
| pelvis_tx (AP) | peak 143.1 N, RMS 24.8 N | peak 155.4 N, RMS 26.4 N | close (~8% peak, ~6% RMS) |
| pelvis_ty (vertical) | peak 172.0 N, RMS 37.5 N | peak 164.2 N, RMS 40.8 N | close (~5% peak, ~9% RMS) |
| pelvis_tz (ML) | peak 19.3 N, RMS 8.2 N | peak 69.6 N, RMS 11.0 N | diverges (~3.6× peak) |

The two dominant axes (AP, vertical) cross-validate the ID tool's translational-residual computation
via a completely independent method. The ML (tz) divergence is a KNOWN, already-diagnosed mechanism
in this repo (`reconcile_diff_scheme.py`: raw-Cartesian-SG-differentiation vs filtered-joint-angle-
spline-differentiation give genuinely different acceleration estimates, most visible on the smallest-
amplitude signal — ML sway is the smallest-amplitude COM motion during walking, so noise is
proportionally largest there) — disclosed, not a new mystery, and not load-bearing for the sagittal-
MOMENT finding this document is centrally about (this check is translational, a different quantity).
A full independent ROTATIONAL momentum cross-check (differentiating whole-body angular momentum
directly) was scope-cut this session — see §7 gap 1.

## 5. Reconciling with RRA's OWN numbers — why baseline looks different, and why that's expected

RRA's own residual-actuator MX/MY/MZ (a DIFFERENT computation: task-gain-optimized tracking control,
not plain ID) were correlated against this session's zero-optimizer ID decomposition, at two existing
configurations already on disk:

| RRA config | predicted-diagonal r (MX↔list, MY↔rotation, MZ↔tilt) | off-diagonal \|r\| |
|---|---|---|
| baseline (uniform gains, boost=150) | MX↔list=0.22, MY↔rotation=0.37, **MZ↔tilt=−0.51** | up to 0.51 |
| v4_extreme (heaviest proximal tuning, weight=50) | **MX↔list=0.82, MY↔rotation=0.97, MZ↔tilt=0.73** | all ≤0.28 |

At v4_extreme the correlation is clean and strong on exactly the predicted diagonal, with small
off-diagonal terms — a genuine, numerical (not just geometric) confirmation of the axis mapping in
§1. At baseline it is weak/mixed (even sign-flipped for MZ↔tilt) — **explained, not just observed**:
a direct re-measurement of pelvis-coordinate tracking error (this session) shows baseline lets
`pelvis_list` drift 11.56° (max\|diff\| vs the real IK) while v4_extreme holds it to 0.34°
(pelvis_tilt: 2.25°→0.37°; pelvis_rotation: 2.65°→0.43°). Under baseline's uniform gains, RRA's own
optimizer is tracking a MEANINGFULLY DIFFERENT trajectory than my zero-optimizer ID used (the real
IK) — so a weak correlation is expected, not a contradiction. As task-gain tuning forces RRA's
tracked kinematics to converge toward the real IK, RRA's own residual converges toward reproducing
the SAME axis-resolved pattern the raw, zero-optimizer computation already shows. **This is also why
the frontal (MX) component was the tunable one** (RRA_TASK_GAINS.md: MX shrank from 92.84→~0 N·m
under heavy proximal weighting) **while sagittal (MZ) was not**: a large fraction of baseline's
frontal residual was plausibly self-inflicted by RRA's own loose pelvis_list tracking under uniform
gains (a fixable allocation choice), whereas the sagittal residual is present even in the zero-
optimizer computation using the REAL trajectory throughout — a floor, not a tracking-slack artifact.

## 6. Candidate (1): GRF-model (free-moment + CoP) consistency — tested and rejected

- **The GRF file's torque_x/z ≡ 0 convention (§1) is standard, not a red flag** — CoP is defined as
  the point nulling the horizontal-plane moment; this is expected for any 6-DOF force-plate wrench.
- **Direct geometric test**: at every stance frame (`Fy > 20 N`), the model's OWN tracked
  ankle-joint-to-mtp-joint span (in the global X=AP/Z=ML axes confirmed in §1) was compared against
  the recorded CoP (px, pz). Violations are negligible:

  | foot | stance frames | max AP violation | frac frames >1cm AP | max ML violation | frac frames >1cm ML |
  |---|---:|---:|---:|---:|---:|
  | R | 114 | 1.69 cm | 0.9% | 0.0 cm | 0.0% |
  | L | 89 | 0.0 cm | 0.0% | 0.33 cm | 0.0% |

  The CoP stays essentially inside the physically-tracked foot throughout stance — no evidence of a
  GRF-vs-kinematics spatial registration mismatch.
- **Friction-plausibility bound** (generous, μ=0.7, r_eff=5cm): recorded free moment peaks at 5.08 N·m
  (R) / 4.71 N·m (L) vs a bound of ~29 N·m at peak vertical force — comfortably inside, i.e.
  physically unremarkable (neither anomalously large nor suspiciously exactly zero).
- **Most direct rejection**: candidate (1) is specifically about the free-moment/CoP data explaining
  transverse-plane (vertical-axis) kinematics — but the transverse-plane residual (MY/pelvis_rotation)
  is the SMALLEST of the three throughout (§2, §5, and the pre-existing RRA evidence) — there is no
  large transverse residual for a GRF/CoP inconsistency to explain in the first place.

## 7. Candidate (2): transverse-plane observability — tested and rejected for THIS residual

**External anchor, fetched live this session**: McGinley JL, Baker R, Wolfe R, Morris ME (2009),
"The reliability of three-dimensional kinematic gait measurements: a systematic review," *Gait &
Posture*, PMID 19013070, DOI 10.1016/j.gaitpost.2008.09.003. Verbatim abstract (fetched via
WebFetch): *"The highest reliability indices occurred in the hip and knee in the sagittal plane, with
lowest errors in pelvic rotation and obliquity and hip abduction. Lowest reliability and highest
error frequently occurred in the hip and knee transverse plane."* This CONFIRMS the general textbook
claim (transverse-plane hip/knee rotation is marker-based gait analysis's classic worst-observed
plane) — **but the stubborn residual here is SAGITTAL (§1), and this same citation states sagittal
hip/knee has the HIGHEST reliability** — the opposite of what candidate (2), so framed, would need.

**Direct, own-data measurement** (not just deferring to literature): ROM, raw 2nd-derivative RMS
normalized by ROM, and spectral high-frequency (≥4 Hz) power fraction, computed on the RAW 100 Hz IK
trajectory for every pelvis+hip rotational DOF plus `knee_angle_r` as a known-good reference:

| coordinate | plane | ROM (deg) | HF power frac (≥4Hz) |
|---|---|---:|---:|
| knee_angle_r (reference) | sagittal | 60.07 | **0.0013** (cleanest by far) |
| hip_flexion_r | sagittal | 39.51 | 0.015 |
| hip_rotation_r | transverse | 13.64 | 0.021 |
| hip_adduction_r | frontal | 7.45 | 0.021 |
| pelvis_rotation | transverse | 6.39 | 0.036 |
| hip_flexion_l | sagittal | 40.87 | 0.042 |
| hip_adduction_l | frontal | 10.16 | 0.043 |
| pelvis_list | frontal | 3.91 | 0.044 |
| **pelvis_tilt** | **sagittal** | 2.83 | **0.045** |
| **hip_rotation_l** | **transverse** | 8.55 | **0.134** (clear outlier, ~3× pelvis_tilt, ~102× knee_r) |

By the cleanest metric (spectral HF fraction, which does not conflate "small genuine ROM" with
"noisy"), the ONE genuine outlier in this trial is `hip_rotation_l` (transverse) — consistent with
the general literature. **`pelvis_tilt`'s HF fraction (0.045) is unremarkable**, similar to
`pelvis_list`/`hip_adduction_l`/`hip_flexion_l`, nowhere near the outlier. By a second metric
(accel-RMS/ROM), `pelvis_tilt` DOES rank #1 (159.4) — but this is explained, not just reported: its
absolute angular-acceleration RMS (450.5 deg/s²) is middling (lower than `hip_flexion_r`,
`hip_rotation_l`, `knee_angle_r`); the high RATIO is driven by its unusually SMALL genuine ROM in
this trial (2.83°, the smallest of all 10 coordinates tested) — a real, small-amplitude sagittal
oscillation during gait, not a large-absolute-noise signature. **Candidate (2), as literally framed
(transverse-plane observability), does not explain a SAGITTAL floor** — pelvis_tilt is not
anomalously noisy by the cleaner measure, and the axis that IS the literature's classic weak point
(hip_rotation_l) is not the one that resisted RRA's tuning.

## 8. Verdict — which candidate, and is it fixable

**The floor traces to candidate (3): an intrinsic inconsistency between the tracked kinematics and
the measured GRF, specifically in the SAGITTAL plane (pelvis_tilt), not the transverse plane the
operator's brief assumed.** It is present with zero optimizer in the loop, robust to a 5× range of
filter cutoffs plus the unfiltered case, and explains (via the tracking-tightness correlation, §5)
why task-gain tuning can shrink the frontal component (partly a tracking-slack artifact) but
structurally cannot touch this one.

**Fixable in software, or an intrinsic data/GRF limit?** Symmetric, not forced: **not fixable by
RRA/CMC parameter tuning** (established across this doc + the two upstream docs it builds on) — but
"intrinsic to this kinematics-GRF pairing" does not by itself mean "intrinsic to the physical
measurement" in the strongest sense. Three distinct root-mechanism sub-hypotheses remain open,
**not distinguished this session**:

1. **Mass-property/scaling error** — `docs/MECHANISM_MSK_ELASTIC_BAND.md` §6 independently found
   `pelvis_tilt` has an unusually LARGE static moment-sensitivity in this exact model (-605 N·m at a
   representative pose, order of the whole system's characteristic body-weight×length scale) —
   consistent with a torso/HAT mass-center or scaling error being dynamically amplified specifically
   through this coordinate. RRA's own auto-COM-position adjustment (`docs/MECHANISM_RRA.md` §3) is a
   partial version of this fix and does NOT resolve the floor (this session's zero-optimizer test
   used the ORIGINAL, pre-RRA-adjustment model and found the same floor) — but the separately-
   reported, never-applied +1.52% uniform mass-magnitude recommendation was not tested. **This is the
   most concrete, ready-to-test next step** and would be software-fixable (a scaling/personalization
   correction) if confirmed.
2. **Sagittal-plane marker/soft-tissue-artifact or scaling error specific to this subject/trial's
   pelvis marker cluster** — a genuine data-quality limit, not fixable downstream of data collection.
3. **Rigid-body model-fidelity limit** (no soft-tissue/footwear compliance, the "rolling knee"
   polynomial coupling, etc.) — an intrinsic modeling-assumption limit.

What this session DOES answer with direct evidence: it is **not** primarily a marker-observability
problem in the classic transverse-plane sense (§7), and it is **not** a GRF/CoP registration or
free-moment-data problem (§6).

## 9. Honest gaps (full list)

1. **No independent ROTATIONAL momentum cross-check was built** (§4 only did the translational
   Newton check) — a from-scratch whole-body angular-momentum-rate computation was scope-cut given
   time budget; the rotational floor's magnitude currently rests on the ID tool alone (itself
   cross-validated via exact reproduction of the pre-existing OpenCap computation and via
   `docs/MECHANISM_MSK_ELASTIC_BAND.md`'s independent virtual-work check on pelvis columns, but that
   check was static-pose, not this dynamic trial).
2. **The 3 root-mechanism sub-hypotheses in §8 were not distinguished** — this session isolates
   WHICH plane (sagittal) and WHICH broad candidate (intrinsic, not GRF-registration, not classic
   marker-observability), not the ultimate physical cause among mass-property error / sagittal
   marker error / rigid-body model-fidelity limit. Testing sub-hypothesis 1 (re-run RRA/ID with the
   already-computed but never-applied +1.52% mass-magnitude recommendation) is the natural, concrete
   next step and was not attempted here (time-boxed).
3. **The RRA-vs-ID correlation (§5) is diagnostic, not a rigorous causal estimate** — both are
   single, non-independent time series from one coupled-dynamics trial (the same caveat
   `docs/MECHANISM_RRA_TASK_GAINS.md` §7 already flagged for its own correlation analysis).
4. **The friction-plausibility bound (§6) uses assumed, not measured, μ and effective contact
   radius** — a generous plausibility check, not a tight physical derivation.
5. **CoP-vs-foot-geometry margins (§6) are heuristic** (heel/toe/width margins chosen for physical
   plausibility, not derived from this model's exact foot mesh, which is unavailable — the `.vtp`
   geometry files are missing from this model's directory, confirmed live via the harmless-but-
   numerous `Couldn't find file *.vtp` warnings every run in this family already produces).
6. **Single trial, right-leg-primary** (subject2 `walking1`) — same scope caveat as every upstream
   cert in this family; no claim of generality across subjects/trials/gait speeds.
7. **The lowpass sweep (§3) did not extend below 4 Hz** — a lower cutoff (e.g. 2-3 Hz) might show a
   different regime; not tested (reasoned scope cut: 4-20 Hz plus unfiltered already spans the
   plausible-to-degenerate range this pipeline would ever actually use).
8. **The observability probe (§7) covers pelvis+hip rotational DOFs plus one sagittal reference
   (knee_angle_r) only** — ankle/subtalar/mtp (independently known to mistrack badly,
   `docs/MECHANISM_ANKLE_RESERVE_FIX.md`) were not included since they are not pelvis/hip and not the
   axis in question here.

## Files

- `scripts/msk/mz_pelvis_moment.py` — the full, re-runnable pipeline (axis-convention verification,
  fresh zero-optimizer Inverse Dynamics + lowpass sweep, independent Newton translational cross-
  check, RRA-axis correlation, CoP-foot-geometry check, friction-plausibility bound, observability
  probe). Reuses `validate_joint_force.py` (`parse_mot`, paths, `savgol_smooth_and_derivs`, `G`) and
  `static_opt_knee.py` (`sub_tag`, `TRIAL_END_TIME`) throughout — never reimplements proven
  parsing/patching logic. CLI: `python3 scripts/msk/mz_pelvis_moment.py` (single `main()`, ~2-3 min).
- `data/msk_smoketest/subject2_walking1/mz_pelvis_moment/` — all outputs (NEW tree): `id_runs/`
  (one subdirectory per lowpass-cutoff ID run: `lp6`, `lp4p0`, `lp6p0`, `lp10p0`, `lp20p0`,
  `lpneg1p0`), `mz_pelvis_moment_results.json` (every number in this document traces back to this
  file).
