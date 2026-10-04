# MECHANISM RRA TASK-GAIN TUNING — does differentiated CMC_Joint weighting clear the Hicks moment gate? (2026-07-21)

Executes the operator's task: `docs/MECHANISM_ANKLE_RESERVE_FIX.md` ruled OUT actuator-capacity as the
cause of RRA's still-failing Hicks 2015 pelvis-**moment** gate (boost=150 vs 300 N·m identical to 4
sig figs) and flagged two untested candidates — (a) the uniform CMC_Joint task gains (`kp=100,
kv=20, weight=1` for all 33 coordinates) or (b) genuine ankle/subtalar/mtp IK-trajectory
inconsistency. This session tests (a): raise tracking weight/gains on the well-observed coords
(pelvis/hip/knee) relative to the poorly-observed (ankle/subtalar/mtp), or down-weight the latter,
holding actuator magnitude FIXED at the already-validated boost=150 N·m set (isolating the
task-gain lever from the actuator-magnitude lever the prior session already exhausted), and
re-runs the 2-pass RRA. Every number below is machine-measured this session
(`scripts/msk/rra_task_gain_tune.py`, reusing `run_rra.py`'s/`ankle_reserve_fix.py`'s own proven
functions throughout — `run_rra_pass`, `measure_residuals`, `resample_kinematics_to_uniform_grid`,
`get_or_compute_anchors` — never reimplemented), cross-checked against an independently-fetched
primary source (Hicks 2015) and independently re-verified via raw-XML regex (not just the same
code path that built the files). Isolation respected: `.venv-msk` only, original model/IK/GRF and
the existing boost=150 actuator file read in place and NEVER written (confirmed via unchanged
md5/mtime, §7), all outputs under a NEW tree,
`data/msk_smoketest/subject2_walking1/rra_task_gain_tune/`, no git operations.

## Headline result

| | force-peak ratio | Hicks force gate | moment-peak ratio | Hicks moment gate |
|---|---:|---|---:|---|
| **BEFORE** (boost=150 actuators, uniform tasks, `MECHANISM_ANKLE_RESERVE_FIX.md`) | 0.0491 | PASS | **0.0957** | **FAIL** |
| **BEST tested** (proximal-only task weight ×20) | 0.0838 | **FAIL** (newly broken) | **0.0220** | **FAIL** (2.2× over) |
| Hicks 2015 threshold (re-verified live, §1) | ≤ 0.05 | | ≤ 0.01 | |

**The moment gate does NOT clear at any of the 9 configurations tested (4 pre-registered variants +
a 5-point dose-response sweep), but the reason is not "task gains don't matter" — they matter a
great deal (up to −77% on the moment ratio) — the reason is a genuine, forced-and-confirmed
force/moment TRADEOFF plus a SATURATION, both directly measured, not assumed.** The exercise also
overturned the operator's/prior-doc's leading candidate direction: down-weighting the "noisy"
ankle/subtalar/mtp coordinates (the (b) option) is nearly INERT (−2.6%) despite tripling their own
tracking error; up-weighting pelvis/hip/knee (the (a) option) is the real, powerful lever, but
cannot be pushed far enough to clear Hicks without first breaking the (previously-passing) force
gate, and it plateaus at 2.2× over threshold regardless. See §5 for the forced, axis-resolved
diagnosis of why.

## 1. External anchor re-verified live (not trusted from the prior session's docstring)

`run_rra.py`'s existing constants (`HICKS_FORCE_FRAC_OF_PEAK_NET_GRF=0.05`,
`HICKS_MOMENT_FRAC_OF_HEIGHT_TIMES_PEAK_GRF=0.01`) were re-confirmed this session via an independent
WebFetch of the primary source, PMC4321112 (Hicks, Uchida, Seth, Rajagopal, Delp, 2015, J Biomech
Eng), Sec 3.1.3, which returns verbatim: *"force discrepancies that are 5% or less (peak and RMS)
than the magnitude of the experimentally measured net external force and residual moments that are
less than 1% of COM height times the magnitude of the measured net external force."* Moment
threshold = 1% (0.01) exactly, force = 5% (0.05) exactly — no discrepancy found; both constants
carried forward unchanged.

## 2. ORIENT, forced before picking a gain recipe (OODA, not skipped)

Before running anything new, the EXISTING boost=150/uniform-tasks data (already on disk) was
machine-cross-checked to test the operator's implicit mechanism ("ankle/subtalar/mtp mistracking
forces the pelvis moment residual") *before* betting a gain recipe on it:

- The moment-residual peak (100.94 N·m at t=1.424s) is a genuine, **sustained** ~200ms plateau, not
  a spline/boundary edge artifact: excluding the last 5% of the trial from the argmax search returns
  the IDENTICAL peak value/timestamp, and 13.7% of the whole trial sits ≥80% of the peak magnitude.
- The peak is 92% dominated by the **MX** component (92.84 of 100.94 N·m).
- Correlating |M(t)| against each coordinate's own |RRA-tracked − real-IK| tracking-error time
  series over the whole trial gives a **mixed, side-asymmetric** picture, not a clean confirmation:
  mean r = 0.483 for the well-observed (hip/knee) group vs only 0.204 for the poorly-observed
  (ankle/subtalar/mtp) group; strong **positive** r for LEFT-side distal coords (ankle_l 0.714,
  subtalar_l 0.539, mtp_l 0.689) but **negative** r for RIGHT-side distal coords (subtalar_r −0.510,
  ankle_r −0.108, mtp_r −0.099).

This tempered (not ruled out) the task-gain hypothesis *before* any variant was run — disclosed
here, not cherry-picked after seeing the fix's result.

## 3. Method — variants (kp/kv/weight fully specified before running)

33 free coordinates partitioned into 3 groups (14 well-observed + 6 poorly-observed + 13
unchanged-trunk/arm = 33, verified to exactly cover the model's free-coordinate set):

| group | coordinates | n |
|---|---|---|
| well-observed | pelvis ×6, hip ×6, knee ×2 | 14 |
| poorly-observed | ankle/subtalar/mtp ×2 sides ×3 | 6 |
| unchanged | lumbar ×3, arm/elbow/pro_sup ×2 sides ×5 | 13 |

Baseline (unchanged group, and every group in the "BEFORE" run): `kp=100, kv=20, weight=1`.

Pre-registered variants (actuators held fixed at the existing, round-trip-verified boost=150 N·m
set — confirmed byte-identical/untouched throughout, §7):

| variant | well-observed | poorly-observed |
|---|---|---|
| v1_downweight_distal | unchanged | `kp=25, kv=10` (half bandwidth, still critically damped), `weight=0.1` |
| v2_upweight_proximal | `weight=10` (kp/kv unchanged) | unchanged |
| v3_combined | v2's proximal change | v1's distal change |
| v4_extreme (escalation) | `weight=50` | `kp=10, kv=6.32, weight=0.02` |

Post-hoc (added after v2/v3 showed a large, dose-dependent, insufficient effect — disclosed as
exploratory, not relabeled as pre-registered): a 5-point proximal-weight-only sweep, `weight` ∈
{2, 3, 5, 20, 30}, kp/kv unchanged, poorly-observed group unchanged — isolating the ONE lever v3
showed actually matters (v2 vs v3 differ by <0.5% — the distal down-weight adds essentially nothing
on top of the proximal up-weight).

Every task file was round-trip verified TWICE: once internally (reload via the OpenSim API, compare
every coordinate's kp/kv/weight against the intended value) and once independently this session via
raw-XML regex (a different code path than the one that built/verified the file) — both confirm
exact, correct differentiation (e.g., `v2_upweight_proximal`: `pelvis_tilt weight=10`,
`ankle_angle_r weight=1`, unchanged).

## 4. Results — all 9 tested configurations

| variant | force ratio | F-gate | moment ratio | M-gate | Δmoment vs before |
|---|---:|---|---:|---|---:|
| BEFORE (uniform) | 0.0491 | PASS | 0.0957 | FAIL | — |
| V1 distal-down (w=0.1) only | 0.0487 | PASS | 0.0932 | FAIL | −2.6% |
| proximal-up w=2 only | 0.0700 | FAIL | 0.0588 | FAIL | −38.5% |
| proximal-up w=3 only | 0.0759 | FAIL | 0.0450 | FAIL | −53.0% |
| proximal-up w=5 only | 0.0689 | FAIL | 0.0331 | FAIL | −65.4% |
| V2 proximal-up (w=10) only | 0.0647 | FAIL | 0.0249 | FAIL | −73.9% |
| **proximal-up w=20 only (BEST)** | 0.0838 | FAIL | **0.0220** | FAIL | **−77.0%** |
| proximal-up w=30 only | 0.0973 | FAIL | 0.0227 | FAIL | −76.3% (worse than w=20) |
| V4 extreme (prox=50 + distal=0.02) | 0.1168 | FAIL | 0.0282 | FAIL | −70.5% (worse than w=20) |
| V3 combined (prox=10 + distal=0.1) | 0.0647 | FAIL | 0.0250 | FAIL | −73.9% (≈V2) |

**Zero of 9 configurations clear the moment gate; zero clear BOTH gates simultaneously.** Force
stays PASS only at baseline and V1 (the near-inert distal-only change) — every configuration that
meaningfully moves the moment ratio breaks the force gate, because baseline force (0.0491) already
sits at 98% of its 0.05 budget with almost no headroom. The moment-ratio improvement itself
**plateaus and reverses** past weight≈20 (w=20: 0.0220 → w=30: 0.0227 → v4 (w=50 combined): 0.0282,
all WORSE than w=20) — a saturation pattern structurally analogous to the actuator-magnitude
saturation `docs/MECHANISM_ANKLE_RESERVE_FIX.md` found (150 vs 300 N·m identical).

## 5. Mechanism, forced via the variants themselves (geometric, machine-measured, not narrated)

**V1 falsifies the "distal mistracking directly forces the pelvis moment" mechanism.** Down-weighting
ankle/subtalar/mtp lets their OWN tracking error roughly TRIPLE (max\|diff\|, resampled-vs-real-IK
convention matching `MECHANISM_RRA.md`'s Sec.4 table):

| coordinate | BEFORE max\|diff\|° | V1 (distal down) | V2 (proximal up) |
|---|---:|---:|---:|
| ankle_angle_r | 18.79 | **59.56** | 2.64 |
| subtalar_angle_r | 25.04 | **102.16** | 4.63 |
| subtalar_angle_l | 52.47 | **138.05** | 4.91 |
| mtp_angle_l | 28.35 | **69.48** | 2.09 |
| hip_flexion_r | 7.95 | 7.85 | **1.97** |
| knee_angle_r | 3.08 | 2.93 | **0.90** |

Yet V1's moment ratio barely moves (−2.6%). If poor distal tracking were directly forcing the
pelvis moment residual, RELIEVING that tracking demand should have shrunk the residual far more than
2.6% — it does not. **This cleanly refutes the "down-weight the noisy distal-foot coords" half of
the operator's hypothesis as the operative mechanism**, via a direct manipulation check, not
assumption.

**V2 reveals the real lever, and a striking side effect.** Up-weighting pelvis/hip/knee (untouching
ankle/subtalar/mtp weights at all) cuts the moment ratio 74% AND, as an unrequested side effect,
distal tracking *also* improves dramatically (subtalar_angle_l: 52.47°→4.91°; mtp_angle_l:
28.35°→2.09°) — better than V1 achieves by directly touching those coordinates' own gains. Read
geometrically: under uniform weighting, the optimizer's weighted least-squares tracking objective
treats a unit of hip/knee tracking error as worth the same as a unit of ankle/mtp tracking error,
but the DYNAMIC cost of correcting them differs (the doc's own low-inertia/no-passive-restraint
argument for mtp, `MECHANISM_ANKLE_RESERVE_FIX.md` §3) — so the optimizer's cheapest solution
under uniform weights lets the whole distal chain drift while nominally satisfying the objective.
Forcing the PROXIMAL chain to anchor tightly to the true kinematics removes that degree of freedom
for the optimizer to exploit, and the distal chain's own tracking improves as a geometric
consequence of a more kinematically-consistent whole-limb solution — not because its own weight
changed. The cost: more pelvis residual FORCE is spent to hit the tighter proximal targets
(35.8N→53.9N), which is what breaks the force gate.

**Axis-resolved diagnosis (geometric, not scalar) — WHY the plateau happens:**

| config | \|M\| peak (N·m) | MX | MY | MZ | \|MX\|/\|M\| |
|---|---:|---:|---:|---:|---:|
| BEFORE | 100.94 | 92.84 | 3.20 | −39.49 | 0.92 |
| V2 (w=10) | 26.31 | 24.22 | 6.82 | 7.66 | 0.92 |
| sweep w=20 | 23.20 | 16.87 | −0.61 | 15.92 | 0.73 |
| sweep w=30 | 23.90 | **−1.69** | −0.87 | **23.83** | **0.07** |
| V4 (w=50 combined) | 29.77 | **2.75** | −0.82 | **29.63** | **0.09** |

At the highest tested weights, the MX component (92% of the ORIGINAL problem) collapses to nearly
zero, but a previously-secondary MZ-axis residual (~24-30 N·m) does NOT shrink correspondingly and
becomes the new, resistant bottleneck — this is precisely why the vector-magnitude ratio plateaus
around 0.022-0.028 instead of continuing toward 0.01. **This task-gain lever fully addresses one
axis of the original problem and is structurally blind to the other** — a specific, falsifiable,
machine-measured finding, not a hand-wave. Which single coordinate's dynamics govern the residual
MZ-axis moment was NOT isolated this session (§7, honest gap for future work).

## 6. Honest answer to the pre-registered claim

**Pre-registered:** C = "≥1 variant drives moment ratio ≤0.01 while force ratio stays ≤0.05." **C is
FALSE** — confirmed across 9 configurations spanning both directions of the operator's stated OR,
a combined variant, an extreme escalation, and a 5-point dose-response sweep on the one lever that
matters (not a single arbitrary attempt). The adversary ("genuine kinematic inconsistency,
independent of anything gain-tunable") does NOT fully hold either, though — a substantial fraction
of the original residual (the MX-axis, 92% of the original 100.94 N·m) IS gain-tunable and
responds with a real, large, monotonic-then-plateauing effect (−77% at the optimum). The honest,
falsifiable, and re-scoped next cause: **a residual MZ-axis pelvis moment (~20-30 N·m, plateaued,
resistant to every tested task-gain configuration up to a 50× weight ratio) is the genuine
next-cause candidate** — narrower and more specific than the operator's original "ankle/subtalar/mtp
IK inconsistency" framing (which V1 directly weakens: relieving distal-tracking demand barely
changes the residual), more likely now a whole-body proximal-chain force/moment cost-allocation
limit compounded with a not-yet-isolated MZ-axis dynamic inconsistency (candidates: GRF-vs-model
consistency, scaling, or a genuine data limit specific to that axis) that this session's tests
cannot further localize without a new, targeted investigation (§7).

## 7. Honest gaps (full list)

1. **No configuration clears the Hicks moment gate** — the headline is a forced, dose-characterized,
   axis-resolved honest-negative on the pre-registered claim, not a partial win dressed up.
2. **The force/moment tradeoff was characterized only along the "proximal weight" axis** (kp/kv held
   at baseline for up-weighted coordinates); raising kp/kv (bandwidth) instead of/alongside weight
   was not sweept — disclosed as a different, untested lever, not assumed equivalent.
3. **Fine-grained weight ∈ (1, 2) was not tested** (reasoned, not measured, scope-cut): at weight=2
   the moment ratio has only reached −38.5% (needs −89.5% to clear 0.01) while force has ALREADY
   failed — even a hypothetical force-preserving weight in (1,2) could not plausibly reach the
   still-far-off 0.01 threshold, so this boundary region was not separately measured.
4. **The MZ-axis residual floor's anatomical/mechanical source was not isolated** — which single
   coordinate(s) or which specific dynamic inconsistency governs it is an open, well-motivated
   follow-on question this session surfaces but does not answer (§5).
5. **Only the RRA leg was re-tested** — the post-RRA muscle-driven SO/JointReaction re-solve and the
   CMC full-trial (both already resolved/characterized in `MECHANISM_ANKLE_RESERVE_FIX.md`) were not
   re-run against these differentiated task files; the operator's task scope was specifically "tune
   task gains, re-run RRA (2-pass)," and this is what was executed.
6. **Single trial, right-leg-primary** (subject2 `walking1`) — same scope caveat as every upstream
   cert in this family; no claim of generality across subjects/trials/gait speeds/models.
7. **The correlation analysis in §2 is diagnostic/orienting, not a rigorous causal estimate** — time
   series from a single coupled-dynamics trial are non-independent and can show spurious
   correlation/anti-correlation from shared gait-cycle phase; it correctly flagged the eventual
   mixed picture but should not be over-read as a precise attribution.
8. **Isolation verified for the reused boost=150 actuator file** (byte-identical md5/mtime,
   confirmed live this session, not assumed) and for all pre-existing tracked files (`git status`
   clean aside from this session's one new script) — but the new `data/` output tree is large
   (9 full 2-pass RRA runs) and was not further pruned.

## Files

- `scripts/msk/rra_task_gain_tune.py` — the full, re-runnable pipeline (`v1_downweight_distal`,
  `v2_upweight_proximal`, `v3_combined`, `v4_extreme` (`escalate`), `sweep` (post-hoc dose-response),
  `all`, `summary` CLI actions). Builds differentiated `CMC_TaskSet` XML fresh via the OpenSim API
  (never `.clone()`), round-trip verified twice (internal + independent raw-XML regex spot-check,
  §3). Reuses `run_rra.py`'s `run_rra_pass`/`measure_residuals`/`resample_kinematics_to_uniform_grid`
  and `ankle_reserve_fix.py`'s `get_or_compute_anchors`/boost=150 actuators/`DIAG_COORDS` throughout
  — no reimplementation.
- `data/msk_smoketest/subject2_walking1/rra_task_gain_tune/` — all outputs (NEW tree): one
  subdirectory per variant (`v1_downweight_distal/`, `v2_upweight_proximal/`, `v3_combined/`,
  `v4_extreme/`, `sweep_w2/`, `sweep_w3/`, `sweep_w5/`, `sweep_w20/`, `sweep_w30/`), each with
  `tasks/`, `rra/pass{1,2}/`, and its own `results.json` (gains, gains_verify, pass1/pass2 residuals,
  tracking_diagnostic, moment_diagnostic) — every number in this document traces back to one of
  these files.
