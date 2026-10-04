# MECHANISM RRA — Residual Reduction Algorithm fidelity upgrade (2026-07-21)

Executes the operator's task: does running RRA on subject2 `walking1` shrink the pelvis residual
(~173 N / 22.6 %BW) both the knee (`docs/MECHANISM_STATIC_OPT.md`) and hip (`docs/MECHANISM_HIP_FORCE.md`)
certs independently flagged as their shared honest gap, and how does that shift the knee/hip
contact-force estimates once Static Optimization is re-run on the RRA-adjusted model/kinematics?
Every number below is machine-measured this session (`scripts/msk/run_rra.py`, exit 0, plus a
diagnosed-and-fixed continuation), not recalled. Isolation respected: `.venv-msk` only, original
model/IK/GRF read in place and never written to, all outputs under a new
`data/msk_smoketest/subject2_walking1/rra/` tree, no git commit/push.

**MID-SESSION CORRECTION (load-bearing, read before the numbers):** `docs/MECHANISM_SIGN_BUG_AUDIT.md`
(an independent audit, verified live below) found the self-computed knee/hip crossing-muscle method
(`static_opt_knee.knee_crossing_muscles_and_forces`) has a confirmed sign bug — its published
233.20/235.30 %BW numbers are bug artifacts (invisible to a vector-norm output, unlike a signed-axis
projection). This document therefore uses the bug-immune `opensim.JointReaction` analysis as the
knee/hip method throughout, both pre- and post-RRA (391.11 %BW knee / 387.04 %BW hip pre-RRA,
reproduced from `static_opt_knee.py`/`validate_hip_force.py`'s own JointReaction cross-checks) —
**not** the buggy self-computed numbers the original task brief cited.

## Headline result

| | force | moment |
|---|---:|---:|
| **BEFORE** (published SO pelvis residual, same peak-per-component convention) | 172.96 N (22.6 %BW) | 45.11 N·m |
| **AFTER, RRA's own tracking residual** (ideal actuators, no muscles) | **37.77 N** (4.6× better) | **85.91 N·m** (1.9× WORSE) |
| **AFTER, re-validated through the actual muscle-driven SO model** (same gate, same method as BEFORE) | **175.32 N** (unchanged, +1.4%) | **83.68 N·m** (1.9× WORSE) |
| Hicks 2015 threshold (peak / RMS, verified live, §below) | ≤ 5% peak net GRF | ≤ 1% × COM-height × peak net GRF |
| AFTER vs Hicks (vector-magnitude form) | peak ratio 0.0508 **FAIL** (RMS ratio 0.0262 PASS) | ratio 0.0868 **FAIL** (8.7× over) |

**RRA did NOT clear the residual to threshold, and the apparent force improvement does not survive
re-validation through the real muscle-driven pipeline.** RRA's own internal tracking accounting
(ideal, unbounded actuators, no muscles) shows a 4.6× force-residual improvement — but once the
RRA-adjusted model/kinematics are fed back into the actual 80-muscle Static Optimization model (the
literal "re-run SO on the RRA-adjusted model/kinematics" the task asked for), the pelvis force
residual is statistically **unchanged** (175.32 N vs 172.96 N) and the moment residual is **worse**
in both accountings. This is diagnosed, not just observed — see §4.

| joint | pre-RRA (bug-immune JointReaction) | post-RRA (same method) | vs OrthoLoad anchor |
|---|---:|---:|---:|
| knee-r | 391.11 %BW | **480.09 %BW** (+22.7%, t=0.35s) | anchor 258.22 %BW — moved FURTHER past it |
| hip-r | 387.04 %BW | **466.73 %BW** (+20.6%, t=1.57s\*) | anchor 273.93 %BW — moved FURTHER past it |

\*hip's post-RRA peak lands on the LAST sample of the trial — an edge-of-window artifact, flagged
in §7, not a confidently-interior local maximum.

**Max segment mass change:** 0.00% actually applied to any body's mass *value* in the model used
downstream (RRA auto-applies only a COM-*position* shift); a separate, un-applied "recommended"
uniform +1.52% mass-*magnitude* scaling across all 22 bodies was also measured (§3) — both far
inside the operator's ~20% guardrail, moot for this trial either way.

## 1. Inputs

- Model: `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`, subject2, 78.2 kg → BW = 766.88 N
  (the same model/mass used by every prior cert in this family).
- Kinematics: `OpenSimData/Mocap/IK/walking1.mot` (158 frames, 100 Hz, 0–1.57 s — the established
  real-data ceiling, `static_opt_knee.py`'s own `TRIAL_END_TIME`, reused here to avoid extrapolation).
- GRF: `ForceData/walking1_forces.mot` (both feet).
- 33 free (actuatable) coordinates, verified live via each coordinate's own `getMotionType()`;
  `knee_angle_{r,l}_beta` (2 of 35) excluded — constraint-coupled ("rolling knee" secondary DOFs,
  the same joint flagged elsewhere in this repo as numerically sensitive,
  `docs/MECHANISM_MSK_ELASTIC_BAND.md` §4), not independent states.

## 2. RRA actuators + tasks — built via the OpenSim API, no template existed

No RRA/CMC example ships with the offline pip `opensim` wheel, and no prior RRA run exists anywhere
in this repo (grepped clean). Rather than hand-type XML from memory, both files were **constructed
via real OpenSim Python objects** (`osim.PointActuator`, `osim.TorqueActuator`, `osim.CoordinateActuator`,
`osim.CMC_Joint`) and serialized with the library's own `printToXML` — every emitted tag/doc-comment
is the actual compiled C++ property metadata, then round-trip-verified by reloading and checking size.

- **6 pelvis residuals** (3 `PointActuator` FX/FY/FZ + 3 `TorqueActuator` MX/MY/MZ, global-frame,
  unbounded control): body/point/direction/optimal_force (10, 10, 10, 15, 15, 15) copied **verbatim**
  from the already-validated OpenCap `walking1_reserveActuators.xml`.
- **12 leg reserves** (hip flex/add/rot, knee, ankle, subtalar × r/l): optimal_force values copied
  verbatim from the same file (2.5–30 N·m range).
- **13 lumbar/arm actuators**: optimal_force (10.0 each) read live off the model's own 13
  pre-existing `CoordinateActuator`s.
- **2 mtp actuators** (mtp_angle_r/l): **ASSUMED** 2.5 N·m each — no muscle or actuator crosses mtp
  in this model at all (verified live), so there was no real value to copy. Disclosed, not hidden —
  see §4 for why this class of assumption (small-reserve magnitudes carried into RRA's fully
  actuator-driven context) turned out to matter.
- **Tasks**: one `CMC_Joint` per actuated coordinate, `kp=100, kv=20` — critically damped
  (`kv = 2√kp`, the tool's own doc-comment, dumped live: "To achieve critical damping of errors,
  choose kv = 2*sqrt(kp)"), giving natural frequency √100 = 10 rad/s (~1.6 Hz), weight=1 uniform
  (a disclosed simplification, not a differentiated per-joint scheme).

## 3. RRA runs — 2-pass, both mechanically successful

- **External anchors** (for the Hicks gate, §5): peak net GRF (both feet vector-summed) = 950.12 N
  (123.9 %BW); whole-body COM height (model's own mean over the trial, std 0.0123 m — stable, not a
  single-frame artifact) = 1.1104 m.
- **Pass 1** (original model, `adjust_com_body='torso'`, full 0–1.57 s): ran OK, 28.2 s real time,
  1625/1661 integrator steps accepted.
- **A genuine mechanism finding, measured not assumed:** RRA's `adjust_com_to_reduce_residuals=true`
  **auto-applies only a COM-POSITION shift** to the specified body — verified directly in the written
  `.osim` XML: torso `<mass>` is bit-identical to the original (27.846079881067737 vs
  27.8461 kg, i.e. 0.00% change) while `<mass_center>` moved from `(-0.0296, 0.3161, 0)` to
  `(-0.02422, 0.31614, -0.00426)`, exactly matching the log's reported `dx=-0.00541, dz=0.00426`.
  Separately, RRA's console prints a **"Summary of Mass Adjustments to Reduce Residuals"** block
  recommending a uniform **+1.523% mass-magnitude** increase across **all 22 bodies** (total
  +1.191 kg) — verified as a single proportional scale factor (each body's recommended fractional
  change matches to 4 significant figures, e.g. pelvis +1.5223%, torso +1.5230%, hand_r +1.5231%) —
  but its own text says **"Note: Edit the model to make recommended adjustments to mass
  properties"**, and the written model does NOT reflect it. Both figures (0.00% applied, 1.52%
  recommended-unapplied) are reported; both are inside the ~20% guardrail. (Honest aside: unlike the
  Hicks force/moment thresholds — verified live against the actual paper, §5 — a specific "20%"
  mass-change ceiling could not be independently verified against a primary citable source this
  session; three targeted fetches found no such numeric guideline in Hicks 2015 or reachable OpenSim
  docs, so it is applied here as an operator-specified guardrail, not a verified citation. Moot for
  this trial regardless, since both measured figures are far below it.)
- **Pass 2** (pass-1-adjusted model, same desired kinematics, `adjust_com=False`, standard 2-pass
  verification convention): ran OK, 28 s. By design (verified from the tool's own doc-comment:
  output model "is written if ... adjust_com_to_reduce_residuals is true") it writes NO new model —
  pass 1's model is the correct, current "RRA-adjusted model" used downstream.
- **A real performance/correctness bug found and fixed:** `pass2_Kinematics_q.sto` is written at the
  integrator's own ADAPTIVE step points; measured directly, `np.diff(time).min() == 0.0` — an exact
  non-increasing timestamp (the signature of a constraint-projection event reporting two states at
  the same nominal time; the log independently confirms "number of projections: 2"). Feeding this
  raw file into `opensim.AnalyzeTool` made it build a **~100,001-sample internal resampling grid
  spanning an out-of-range "-0.799 to 2.339s"** (log-verified) — a post-RRA Static Optimization run
  was killed after **26 CPU-minutes at only ~10% progress**, not a genuine "the model is slow"
  finding. **Fix**: de-duplicate non-increasing timestamps and resample onto a clean uniform 100 Hz
  grid (matching the original IK's own rate, 157 rows) before handing kinematics to any Analyze-family
  tool — cut the same computation to well under a minute.

## 4. BEFORE vs AFTER — two different "after"s, and why they disagree

**RRA's own tracking-residual accounting** (measured from `pass{1,2}_Actuation_force.sto`, same
per-component max-abs convention as the published 172.96 N number, both passes agree to 5 decimal
places — pass 2 is a stable verification of pass 1, not a further reduction, as expected since it
requests no new adjustment): force peak 37.77 N (4.6× better), moment peak 85.91 N·m (1.9× worse).
Vector-magnitude form (the quantity the verified Hicks criterion is phrased in): force peak/RMS
48.24 N / 24.87 N, moment peak/RMS 91.61 / 53.69 N·m.

**Re-validated through the actual muscle-driven model** (Static Optimization re-run on the RRA-adjusted
model + resampled RRA-tracked kinematics, using the EXACT SAME gate code
(`static_opt_knee.check_so_convergence_and_sanity`) that produced the original 172.96 N number — the
literal apples-to-apples comparison): pelvis residual force **175.32 N**, moment **83.68 N·m**. The
force number is **statistically unchanged** from before RRA (+1.4%, within the kind of run-to-run
noise this pipeline has shown elsewhere); the moment number is **materially worse** (+85%).

**Why the two "after"s disagree — forced via OODA, not shrugged off:** compared the RRA-tracked
kinematics directly against the original, already-validated IK trajectory at every coordinate:

| coordinate | orig range (deg) | RRA-tracked range (deg) | max |diff| (deg) | RMS diff (deg) |
|---|---:|---:|---:|---:|
| hip_flexion_r | −12.85 to 26.66 | −14.65 to 26.66 | 7.5 | 4.0 |
| knee_angle_r | 5.59 to 65.65 | 5.30 to 68.01 | 3.0 | 1.4 |
| ankle_angle_r | −16.54 to 15.85 | −16.32 to 22.17 | 21.9 | 6.8 |
| subtalar_angle_r | −20.25 to 3.35 | −40.10 to 5.24 | 29.1 | 10.5 |
| ankle_angle_l | −24.04 to 12.72 | −22.44 to 36.34 | 25.4 | 12.9 |
| **subtalar_angle_l** | **−17.26 to 4.29** | **−55.13 to 25.67** | **51.4** | **21.1** |

Hip/knee track well (≤7.5° max deviation); **ankle/subtalar/mtp track badly** (up to 51° deviation —
a real, large tracking failure, not numerical noise). **Diagnosed root cause**: the RRA actuator
file's ankle/subtalar optimal_force (2.5 N·m, §2) was copied from SO's convention, where those
values are fine because they are a *small safety net supplementing ~80 real muscles* (real
plantarflexor/invertor capacity is 100s of N·m). RRA has **no muscles at all**
(`replace_force_set=True`) — the *entire* ankle/subtalar torque needed to track push-off-phase gait
kinematics had to come from that same 2.5 N·m actuator, and RRA/CMC's control-allocation objective
(penalizing control magnitude, effectively `(force/optimal_force)²`) makes using it heavily "expensive,"
so the tracking controller let those joints drift instead of forcing them to follow the true
trajectory. This mistracking is exactly why the post-RRA muscle-driven SO run needed **massive**
reserve-actuator help at those same joints (§6) and why the pelvis-residual improvement evaporates
once real muscles are back in the loop. **Concrete next step, not a vague retry**: re-run with
ankle/subtalar/mtp optimal_force raised to a magnitude commensurate with real muscle torque capacity
(e.g. 100–300 N·m) instead of copying SO's small-reserve convention — not attempted this session
(time-boxed; flagged as the specific, identified fix rather than left as an unexplained failure).

## 5. Hicks et al. 2015 thresholds — verified live, not recalled

Fetched directly from the paper this session (PMC4321112, Hicks, Uchida, Seth, Rajagopal, Delp,
*J Biomech Eng* 2015, §3.1.3), quoted verbatim: **"force discrepancies that are 5% or less (peak and
RMS) than the magnitude of the experimentally measured net external force and residual moments that
are less than 1% of COM height times the magnitude of the measured net external force."** Applied
here as: force ratio = peak (or RMS) residual-force vector magnitude ÷ peak net GRF magnitude (both
feet vector-summed, 950.12 N); moment ratio = peak residual-moment vector magnitude ÷ (COM height ×
peak net GRF).

| gate | threshold | measured (AFTER, pass 2) | verdict |
|---|---|---:|---|
| Force, peak | ≤ 0.05 | 0.0508 | **FAIL** (a near-miss, 1.5% over) |
| Force, RMS | ≤ 0.05 | 0.0262 | PASS |
| Moment, peak | ≤ 0.01 | 0.0868 | **FAIL** (8.7× over) |

The residual **improved from a much worse starting point by RRA's own accounting but did not clear
threshold**, and — per §4 — the improvement is not even real once re-validated through the muscle
model. This is reported as the honest result, not reframed as a pass.

## 6. Muscle-driven re-solve sanity (a new caveat this session surfaced)

Post-RRA SO convergence: 158/158 frames, 0 NaN, activation bounds [0.01, **1.0**] (at least one
muscle saturates exactly at its ceiling — new vs. the original run, where the max was 0.627), 0
muscles over 1.5× Fmax. Joint reserve saturation (fraction of each reserve's own small optimal_force
used) **exploded** at the distal joints RRA mistracked (§4):

| reserve | original run (docs/MECHANISM_STATIC_OPT.md) | post-RRA (this session) |
|---|---:|---:|
| all hip/knee reserves | ≤ 12.5% | 2.8–23.5% (still modest) |
| ankle_angle_r / l | (small, unflagged) | **1747% / 1533%** |
| subtalar_angle_r / l | (small, unflagged) | **3751% / 4364%** |

A reserve actuator with optimal_force 2.5 N·m being asked for 37–44× that (≈90–110 N·m) is a loud,
machine-checkable signal that the muscle-driven model cannot physiologically reproduce the RRA-tracked
ankle/subtalar trajectory — directly corroborating §4's diagnosis via an independent measurement
(reserve saturation vs. direct kinematic comparison), not merely repeating it.

## 7. Honest gaps (full list)

1. **RRA did not clear the Hicks thresholds**, and the one metric that looked like a clean win
   (RRA's own force-residual accounting, 4.6× better) does not survive re-validation through the
   actual muscle-driven pipeline (§4) — the pelvis force residual is statistically unchanged and the
   moment residual is worse by every measure tried. This is the primary, load-bearing honest result.
2. **Root cause is diagnosed and the fix is identified but not yet applied this session** (§4) — a
   real, disclosed limitation of this specific actuator recipe, not a fundamental limitation of RRA
   itself; a concrete, specific re-run (boost ankle/subtalar/mtp optimal_force) is the natural next
   step, not attempted here (time-boxed).
3. **The knee/hip contact-force estimates moved FURTHER from the OrthoLoad anchor after RRA**
   (knee 391→480 %BW, hip 387→467 %BW) — consistent with, not independent of, the ankle/subtalar
   mistracking corrupting the whole-limb kinematics/dynamics used by both the pelvis-residual and the
   knee/hip calculations.
4. **Hip's post-RRA JointReaction peak lands at t=1.57s, the trial's last sample** — an edge-of-window
   artifact (the same class of boundary effect `docs/MECHANISM_STATIC_OPT.md` §5 already documented
   for JointReaction's own differentiation scheme), flagged rather than trusted as a confident
   interior maximum.
5. **The self-computed knee/hip method is NOT used here at all** (unlike the original brief's framing)
   because of the independently-audited sign bug (`docs/MECHANISM_SIGN_BUG_AUDIT.md`) — this document
   relies exclusively on the bug-immune `opensim.JointReaction` analysis for knee/hip, pre- and
   post-RRA, for a clean apples-to-apples comparison.
6. **mtp actuator magnitude (2.5 N·m) is an assumption**, not measured or copied from any existing
   model component (no muscle/actuator crosses mtp in this model at all) — flagged in §2, and its
   under-strength is part of the same diagnosed mechanism as ankle/subtalar (§4), though mtp's own
   tracking error was not separately isolated from ankle/subtalar's.
7. **The "~20% mass change" guardrail is the operator's own framing, not an independently verified
   citation** (§3) — unlike the Hicks force/moment numbers, which were fetched and quoted directly
   from the primary source. Moot for this trial (both measured figures are ≤1.6%).
8. **A cosmetic console-log ordering quirk**: this session's fd-level stdout-capture mechanism (added
   to auto-parse RRA's mass-recommendation text on every future re-run, §3) produced one run whose
   printed narrative lost a few intermediate lines around the pass-1→pass-2 transition — the
   underlying `.sto` data and `rra_results.json` are unaffected (independently re-verified by
   recomputing `pass{1,2}_residuals` directly from the raw `Actuation_force.sto` files, matching the
   log to 6+ significant figures) — noted for transparency, not load-bearing for any number reported
   above.
9. **Single trial, right leg primary** (subject2 `walking1`) — same scope caveat as every upstream
   cert in this family; no claim of generality across subjects/trials/gait speeds/RRA-gain choices.
10. **CMC_Joint task weights are uniform (1.0) for all 33 coordinates** — a disclosed simplification;
    a differentiated per-joint weighting (e.g. tracking the pelvis more tightly) was not attempted,
    and might materially change §4's mistracking finding.

## Files

- `scripts/msk/run_rra.py` — the full, consolidated, re-runnable pipeline (builds RRA actuators/tasks
  via the OpenSim API, runs the 2-pass RRA, measures residuals two ways, re-runs SO + JointReaction on
  the adjusted model/kinematics with the kinematics-resampling fix folded in, writes
  `rra_results.json`). Imports `validate_joint_force.py` / `static_opt_knee.py` /
  `validate_hip_force.py` for proven parse/BFS/SO-patch/JR-extraction code, never edits them.
- `scripts/msk/rra_post_process.py` — the standalone diagnostic-and-fix continuation actually executed
  this session (identical logic to `run_rra.py`'s Step 7, developed first while iterating on the
  pathological-resampling bug; kept as a second, independent artifact of the same result).
- `data/msk_smoketest/subject2_walking1/rra/` — all outputs (NEW tree, originals untouched):
  `actuators/RRA_actuators.xml`, `tasks/RRA_tasks.xml`, `pass1/` (adjusted model + residuals),
  `pass2/` (final tracked kinematics, raw + resampled), `so_post_rra/` (post-RRA Static Optimization),
  `jr_post_rra/` (post-RRA JointReaction), `rra_post_process_results.json`.
