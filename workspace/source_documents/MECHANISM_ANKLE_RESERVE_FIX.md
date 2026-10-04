# MECHANISM ANKLE RESERVE FIX — boosting ankle/subtalar/mtp reserve optimal_force (2026-07-21)

Executes the operator's task: two independent threads (`docs/MECHANISM_RRA.md` §4, `docs/
MECHANISM_CONTRACTION_DYNAMICS.md` §8) converged on the SAME diagnosed cause — RRA's/CMC's
ankle/subtalar/mtp reserve actuators inherited Static Optimization's tiny 2.5 N·m
muscle-*supplementing* convention into contexts with no (RRA) or insufficient (CMC) muscle
support at those joints — and both flagged the identical, not-yet-applied fix: raise
ankle/subtalar/mtp reserve `optimal_force` to ~100-300 N·m. This session applies that fix (a new,
boosted actuator/reserve-file pair, `scripts/msk/ankle_reserve_fix.py`) and re-runs both RRA (2-pass)
and the full 1.57 s CMC trial. **The fix is a clean, decisive, disk-verified WIN for CMC (question 2)
but only a genuine, partial, mechanistically-diagnosed win for RRA's own pelvis residual (question 1)
— an asymmetry that is itself explained below, not left as a loose end.** CMC now runs the ENTIRE
1.57 s trial to completion (`ran_ok=True`, 518 s wall time; verified directly on the raw output files,
not just the boolean: data span 0.03-1.56 s, zero NaN/Inf, `Ipopt: Restoration failed` no longer
appears anywhere in the console capture) — a direct, disk-confirmed fix of the t=0.67s/56.1%GC
infeasibility. RRA's own force-residual Hicks gate flips FAIL to PASS and the post-RRA muscle-driven
reserve-saturation pathology resolves (1533-4364% → 63-93% of capacity), but the Hicks **moment** gate
still FAILS (and gets *worse*), and the ankle/subtalar/mtp kinematic-tracking error is only partially
resolved (mixed by side). The single most decisive finding this session, forced via an OODA follow-up
rather than accepted at the first data point: **the RRA-side effect completely SATURATES between 150
and 300 N·m** (every measured quantity is numerically identical to 3-4 significant figures at both
boost levels), proving the residual moment-gate failure and residual mistracking are **not**
actuator-capacity-limited within this fix family — the reason CMC fully resolves while RRA only
partially does is diagnosed in §5-§6. Every number below is
machine-measured this session (`scripts/msk/ankle_reserve_fix.py`, reusing `run_rra.py`/
`static_opt_knee.py`/`validate_hip_force.py`/`contraction_dynamics_cmc_attempt.py`'s own proven
functions throughout, never reimplemented), cross-checked live against the existing docs' cited
numbers before being trusted (§1). Isolation respected: `.venv-msk` only, original model/IK/GRF/
`RRA_tasks.xml` read in place and never written to, all outputs under a NEW tree,
`data/msk_smoketest/subject2_walking1/ankle_reserve_fix/`, no git operations.

## Headline result

| question | answer |
|---|---|
| **(1) Does RRA's pelvis residual now clear Hicks 2015** (vs prior force 0.0508/moment 0.0868 FAIL)? | **PARTIAL.** Force-peak ratio **0.0491 — PASS** (narrowly, was 0.0508 FAIL). Moment-peak ratio **0.0957 — still FAIL**, and *worse* than before (was 0.0868). |
| **(2) Does CMC now converge past the 56%GC push-off window to completion?** | **YES.** Full 1.57 s trial `ran_ok=True`, 518 s wall time, verified on the raw output files (span 0.03-1.56 s, 2786 rows, 0 NaN, 0 Inf). No Ipopt restoration failure anywhere in the console capture (vs the prior weak-reserve run's hard failure at t=0.67s/56.1%GC after only 240 s / 43% of the trial). §5. |
| **(3) If either still fails (the RRA moment gate does), the honest next cause?** | **Not actuator capacity.** Boost=150 and boost=300 N·m (the full width of the operator-identified range) give numerically IDENTICAL results on every RRA metric (§2-§3) — the mechanism has already saturated by 150 N·m, and CMC's own reserve usage during the successful full trial is modest (9-24% of the 150 N·m budget at ankle/subtalar, 0.2% at mtp, §5) — nowhere near saturated either. The honest next-cause candidates for the remaining RRA moment-gate failure (§6): the uniform CMC_Joint task gains (`kp=100/kv=20`, `weight=1` for all 33 coordinates, already disclosed as a simplification in `docs/MECHANISM_RRA.md` §7.10), or genuine kinematic inconsistency in the ankle/subtalar/mtp IK trajectory itself (skin-marker IK is known to poorly resolve these specific rotations) — RRA has no muscles at all so it bears 100% of that inconsistency, while CMC's real muscles absorb most of it and the small reserves only mop up the residual gap (§5-§6). |

| RRA's own tracking residual | force (N) | moment (N·m) | Hicks force-peak | Hicks force-RMS | Hicks moment-peak |
|---|---:|---:|---|---|---|
| **BEFORE** (weak, 2.5 N·m ankle/subtalar/mtp — `docs/MECHANISM_RRA.md`, cross-checked live §1) | 37.77 | 85.91 | 0.0508 **FAIL** | 0.0262 PASS | 0.0868 **FAIL** |
| **AFTER, boost=150 N·m** | 35.80 | 95.52 | 0.0491 **PASS** | 0.0265 PASS | 0.0957 **FAIL** (worse) |
| **AFTER, boost=300 N·m** | 35.80 | 95.53 | 0.0491 **PASS** | 0.0265 PASS | 0.0957 **FAIL** (worse) |

150 vs 300 N·m agree to 4 significant figures on every column — see §2.

## 1. Baseline cross-check (measured live before trusting any comparison)

Before running anything new, the existing weak-reserve (2.5 N·m) RRA result already on disk
(`data/msk_smoketest/subject2_walking1/rra/pass2/`) was **re-measured this session** with the exact
same reused functions (`run_rra.measure_residuals`, same `peak_net_grf_n`=950.12 N /
`com_height_m`=1.1104 m anchors) that will be applied to the boosted run, to confirm the comparison
is apples-to-apples before trusting it:

| quantity | doc-cited (`MECHANISM_RRA.md`) | live-recomputed this session | rel. error |
|---|---:|---:|---:|
| force peak (N) | 37.77 | 37.7729 | 0.008% |
| moment peak (N·m) | 85.91 | 85.9093 | 0.001% |
| Hicks force-peak ratio | 0.0508 | 0.0508 | 0.054% |
| Hicks moment-peak ratio | 0.0868 | 0.0868 | 0.031% |

A first attempt at cross-checking the ankle/subtalar tracking-error table (§4 of `MECHANISM_RRA.md`)
against the **raw** (adaptive-step) `pass2_Kinematics_q.sto` reproduced 2/4 coordinates to <0.1% but
was 2-8% off for `ankle_angle_r`/`subtalar_angle_r` (22.37/31.39° live vs 21.9/29.1° cited) — forced
via a cheap OODA check (tried the alternative interpolation direction, then tried the RESAMPLED-100Hz
kinematics file instead of raw): the doc's own table used the **resampled** file (what actually flows
downstream into the SO re-solve), not the raw adaptive-step one. Switching closed all 4 coordinates
to <0.2% agreement (21.87/29.15/25.39/51.37° vs cited 21.9/29.1/25.4/51.4°). `ankle_reserve_fix.py`
uses the resampled convention throughout for this reason.

## 2. Mechanism verified live before any expensive run

A first attempt to build the boosted SO-style reserve file via `osim.ForceSet.clone()` **segfaulted**
(OpenSim 4.6 Python SWIG bindings) — caught by a cheap, throwaway pre-check before wiring it into the
real pipeline, not discovered mid-run. Fix: build every actuator **fresh** via the OpenSim API
(`osim.PointActuator`/`TorqueActuator`/`CoordinateActuator` constructors), the same pattern
`run_rra.py`'s own `build_rra_actuators` already uses successfully — never clone an existing `Force`
object. Both new builders (`build_boosted_rra_actuators_xml`, `build_boosted_so_reserve_xml`) were
verified round-trip-correct (reload from disk, check every targeted coordinate's `optimal_force`
against the intended value, check every untouched coordinate is unchanged) before use.

Two boosted files were built, matching the two consumers already established in this repo:

- **Boosted RRA actuators** (33 elements, `RRA_actuators_boosted.xml`): identical structure to the
  existing, validated `rra/actuators/RRA_actuators.xml`, `ankle_angle_{r,l}`/`subtalar_angle_{r,l}`/
  `mtp_angle_{r,l}` optimal_force raised 2.5 → BOOST N·m, every other actuator (6 pelvis residuals,
  hip/knee reserves, 13 lumbar/arm) **unchanged**. Used only for the RRA passes.
- **Boosted SO-style reserve** (20 elements, `SO_reserve_boosted.xml`): mirrors the existing
  `walking1_reserveActuators.xml`'s 18 elements verbatim except the same 4 ankle/subtalar coords
  boosted, **plus 2 NEW mtp `CoordinateActuator`s** (the original file has none at all — verified live,
  no muscle or actuator crosses mtp in this model). Used for BOTH the post-RRA muscle-driven SO
  re-solve and the CMCTool full-trial attempt (both are "small reserve alongside real muscles,
  `replace_force_set=False`" contexts, exactly like the existing pre-fix convention already shared
  between `run_rra.py`'s post-RRA SO step and `contraction_dynamics_cmc_attempt.py`'s CMC step).
  Existing `RRA_tasks.xml` (task gains, orthogonal to actuator magnitude) reused **unchanged** for
  both RRA and CMC, per the operator's instruction.

**The decisive OODA follow-up — forced, not skipped, because the first (boost=150) result was a
mixed win, not a clean pass:** rather than stop at one data point, a second run at boost=300 N·m (the
top of the operator-identified range) was run for the cheap RRA leg. The two are **numerically
indistinguishable**:

| | boost=150 | boost=300 | difference |
|---|---:|---:|---:|
| pass2 force peak (N) | 35.796284 | 35.795796 | 0.0014% |
| pass2 moment peak (N·m) | 95.523962 | 95.526121 | 0.0023% |
| Hicks force-peak ratio | 0.049053 | 0.049053 | <0.001% |
| Hicks moment-peak ratio | 0.095672 | 0.095674 | <0.001% |
| `ankle_angle_l` max\|diff\| | 25.115° | 25.115° | <0.01° |
| `subtalar_angle_l` max\|diff\| | 52.472° | 52.473° | <0.01° |
| `mtp_angle_l` max\|diff\| | 28.355° | 28.354° | <0.01° |

Verified this is not a caching/no-op bug: the two actuator XML files were diffed directly and confirmed
to genuinely differ (`ankle_angle_r`/`subtalar_angle_r`/`mtp_angle_r` optimal_force = 150 vs 300 in
the respective files; every other actuator, e.g. `knee_angle_r` = 2.5, byte-identical between them),
and both RRA runs show fresh console/integration output (not a reused cache).

**Geometric read of the saturation:** every actuator in both files has `min_control=-Inf,
max_control=Inf` — `optimal_force` is not a hard torque ceiling here, it only rescales the *quadratic
control-effort cost* `(force/optimal_force)²` RRA's/CMC's allocation objective penalizes. Complete
insensitivity between a 2× change in that cost-scaling (150→300) means the allocator is no longer
being limited by that cost term at all in this regime — whatever residual mistracking and
moment-residual failure remain is being set by something else entirely (task gains or genuine
kinematic/dynamic inconsistency, §6), not by how expensive it is, in the optimizer's soft penalty, to
use these actuators.

## 3. Did boosting fix the actual tracking MECHANISM, or just move an aggregate number?

Direct, per-coordinate comparison of RRA-tracked kinematics vs the real IK trajectory (same method
`docs/MECHANISM_RRA.md` §4 used, resampled-kinematics convention, §1):

| coordinate | weak (2.5 N·m) max\|diff\| | boost=150 max\|diff\| | boost=300 max\|diff\| | change vs weak |
|---|---:|---:|---:|---|
| `ankle_angle_r` | 21.87° | 18.79° | 18.79° | **-14%** (improved) |
| `subtalar_angle_r` | 29.15° | 25.04° | 25.04° | **-14%** (improved) |
| `ankle_angle_l` | 25.39° | 25.12° | 25.12° | -1% (~flat) |
| `subtalar_angle_l` | 51.37° | 52.47° | 52.47° | **+2% (worse)** |
| `mtp_angle_r` | 27.24°\* | 24.68° | 24.68° | -9% (modest) |
| `mtp_angle_l` | 28.51°\* | 28.35° | 28.35° | -0.6% (~flat) |

\*mtp was not in the original doc's Sec.4 table; live-measured this session from the same on-disk
weak-reserve kinematics for a fair before/after (mtp already had a 2.5 N·m actuator in RRA's own
actuator set pre-fix — only the CMC-side reserve file lacked one entirely, added in §2).

**Honest read: a real but small, side-asymmetric improvement, not a fix.** The right leg improves
modestly (~14% at ankle_r/subtalar_r); the left leg is flat-to-slightly-worse. None approach the
hip/knee tracking fidelity (3.1-7.9° max deviation, unaffected by this fix as expected since those
actuators were not touched). `mtp` is the most surprising: its **desired** trajectory has almost no
range at all (-7.2° to -4.9°, a 2.3° span — consistent with a near-rigid toe during normal shod gait)
yet its RRA-tracked trajectory swings ~54° regardless of a 60× actuator boost (2.5→150 N·m).

**A specific alternative hypothesis (mtp coordinate-limit violation) was tested and ruled out**: mtp's
hard range is -45° to +30° (`getRangeMin/Max`, verified live); the tracked excursions (-35.1° to
+33.8°) stay inside that range — not a clamping artifact. The toe body itself is small (mass 0.225 kg,
`toes_r`/`toes_l`, verified live) with **no passive `CoordinateLimitForce`** restoring it toward
neutral — a plausible (disclosed, not independently proven) reason mtp shows large excursions
regardless of its own actuator's optimal_force: a low-inertia, passively-unrestrained DOF is
dynamically "cheap" to perturb by any leftover whole-body force/moment inconsistency, a property of
its *mass*, not of its *actuator cost*, which boosting optimal_force cannot address.

## 4. Muscle-driven re-solve (post-RRA Static Optimization + JointReaction)

Boosted SO-reserve file (§2) used for the SO re-solve on the RRA-pass1-adjusted model + resampled
pass-2 kinematics — the literal "re-run SO on the RRA-adjusted model" comparison, same gate code
(`static_opt_knee.check_so_convergence_and_sanity`) as every prior cert in this family:

| | pelvis force peak (N) | pelvis moment peak (N·m) |
|---|---:|---:|
| weak-reserve (2.5 N·m), `docs/MECHANISM_RRA.md` | 175.32 | 83.68 |
| **boost=150 N·m** | 174.39 (~unchanged, -0.5%) | 94.37 (**worse, +12.8%**) |

Same qualitative pattern as RRA's own residual (§Headline): force ~flat, moment worse. **But the
mechanism-level reserve-SATURATION pathology `docs/MECHANISM_RRA.md` §6 flagged is genuinely
resolved** — the muscle-driven model's own demand on these reserves, which was asking for 15-44× the
weak reserves' capacity (a physically-impossible ask), now sits at a plausible 63-93% of the new,
boosted capacity:

| reserve | weak-reserve saturation (docs/MECHANISM_RRA.md §6) | boost=150 saturation |
|---|---:|---:|
| `ankle_angle_r_reserve` | 1747% | **67.3%** |
| `ankle_angle_l_reserve` | 1533% | **63.3%** |
| `subtalar_angle_r_reserve` | 3751% | **89.0%** |
| `subtalar_angle_l_reserve` | 4364% | **92.7%** |

This is the fix's cleanest, most unambiguous win: it directly confirms the originally-diagnosed
mechanism (reserves being asked for far more torque than physically available) was real and is now
resolved at the muscle-driven-model level, even though the aggregate pelvis-residual number does not
show a matching clean improvement (the "extra" torque these joints now supply changes the whole-body
force balance in a way that helps the FORCE residual marginally but costs the MOMENT residual more,
§Headline).

**Knee/hip JointReaction (bug-immune method, vs OrthoLoad in-vivo anchor, `docs/MECHANISM_STATIC_OPT.md`/
`MECHANISM_HIP_FORCE.md`):**

| joint | no-RRA baseline | weak-RRA | boost=150-RRA | OrthoLoad anchor |
|---|---:|---:|---:|---:|
| knee-r | 391.11 %BW | 480.09 %BW | 475.85 %BW | 258.22 %BW |
| hip-r | 387.04 %BW | 466.73 %BW | 470.88 %BW | 273.93 %BW |

Essentially unchanged by this fix (knee -0.9%, hip +0.9% vs weak-RRA) — still substantially further
from the OrthoLoad anchor than the no-RRA baseline, consistent with `docs/MECHANISM_RRA.md`'s own §7.3
finding that this movement is driven by the same whole-limb kinematic/dynamic disturbance, not
independently fixable by this specific actuator change.

## 5. CMC full 1.57 s trial — YES, now converges to completion

Same `CMCTool` configuration as `contraction_dynamics_cmc_attempt.py`'s already-proven recipe
(`replace_force_set=False`, real 80 muscles kept, existing `RRA_tasks.xml` reused unchanged) but with
the boosted SO-style reserve file (§2) in place of the weak 2.5 N·m one, plus fd-level stdout+stderr
capture added this session (the prior weak-reserve session's own `opensim.log` for this exact run was
not actually present in its results directory when checked this session — a stray, unrelated
root-level `opensim.log` from interactive snippets was found instead and is unrelated to any CMC run;
the fd-capture wrapper added here guarantees a reliable, grep-able record regardless).

| | weak-reserve (2.5 N·m), `docs/MECHANISM_CONTRACTION_DYNAMICS.md` §8 | boost=150 N·m (this session) |
|---|---|---|
| `ran_ok` | **False** | **True** |
| wall time | 239.5 s | 517.9 s |
| outcome | `Ipopt: Restoration failed (status -2)`, "Unable to find a feasible solution at time = 0.67" | completes the full trial, no restoration-failure text anywhere in console capture |
| trial coverage reached | t=0.67s / 1.57s = 42.7% (56.1%GC) | t=1.56s / 1.57s = **99.4%** |

**Verified directly on the raw output files, not trusted from the boolean alone** (machine
cross-check, not narration): `subject2_walking1_cmc_full_1p57s_Actuation_force.sto` parsed with
numpy — shape (2786, 114), time column spans **[0.03, 1.56] s** (the requested window is 0-1.57s;
CMC's own internal integration grid starts a few steps in and the file legitimately ends at the last
completed step), **0 NaN, 0 Inf**, max\|value\| across all non-time columns = 1437.6 (a plausible
muscle-force magnitude in Newtons, not a divergence artifact). `grep -c "CMC::computeControls"` on the
console capture = 155 (vs an implied ~68 for the old failure point at t=0.67s) — the run genuinely
integrated almost the full window, not a truncated/silently-early-exit success.

**Reserve-actuator usage during the successful run** (peak\|force\| as % of each reserve's own
`optimal_force`, confirming the fix worked with real margin, not a knife-edge pass):

| reserve | peak\|force\| (N) | % of 150 N·m budget |
|---|---:|---:|
| `ankle_angle_r_reserve` | 31.8 | 21.2% |
| `ankle_angle_l_reserve` | 35.7 | 23.8% |
| `subtalar_angle_r_reserve` | 14.2 | 9.4% |
| `subtalar_angle_l_reserve` | 10.5 | 7.0% |
| `mtp_angle_r_reserve` | 0.28 | 0.2% |
| `mtp_angle_l_reserve` | 0.29 | 0.2% |

Comfortable headroom at every one of the previously-flagged joints (well under 25% of the boosted
capacity) — CMC's real muscles are doing almost all the work, exactly as `replace_force_set=False`
intends; the reserves only needed to supply a modest top-up. This is the mechanistic reason CMC
resolves so much more cleanly than RRA (§2-§4): RRA has **no muscles at all**, so its actuators (even
boosted) must supply the ENTIRE ankle/subtalar/mtp torque and visibly saturate the tracking objective
(§2-§3); CMC's muscles absorb the bulk of the demand, leaving the reserves a genuinely small,
comfortably-resourced safety net. (By contrast, the 6 pelvis-residual actuators — never touched by
this fix — are heavily leaned-on in both contexts: 213-1271% of their small 10/15 N optimal_force
during this same CMC run, consistent with this whole cert family's established finding that the
pelvis residual, not the leg reserves, is where this trial's dynamic inconsistency concentrates.)

**Smoke test (0-0.3s):** `ran_ok=True`, 110.2 s (vs the weak-reserve run's 118.1 s — no regression),
run first as the cheap decisive gate before committing to the ~9-minute full trial, per this
pipeline's own convention.

## 6. Honest answers to the operator's 3 questions

1. **Does RRA's pelvis residual now clear Hicks 2015 (vs prior 0.0508/0.0868 FAIL)?** PARTIAL, not a
   clean pass. Force-peak ratio flips to PASS (0.0491, narrowly under the 0.05 threshold). Moment-peak
   ratio is still FAIL and gets WORSE (0.0868→0.0957, now 9.6× over the 0.01 threshold vs 8.7× before).
2. **Does CMC now converge past 56%GC to completion?** **YES** — verified on the raw output files, not
   just `ran_ok`: the full 1.57s trial completes (`ran_ok=True`, 518s wall time), data span 0.03-1.56s
   (99.4% of the trial), 0 NaN/0 Inf, and no `Ipopt: Restoration failed` text anywhere in the console
   capture (§5). This is a clean, decisive fix of the specific failure `docs/
   MECHANISM_CONTRACTION_DYNAMICS.md` §8 diagnosed.
3. **Honest next cause for the ONE thing that still fails (RRA's moment gate), forced via OODA rather
   than asserted:** the 150-vs-300 N·m saturation (§2) is the decisive datum — this rules out "needs an
   even bigger boost" as the next step. §5's reserve-usage numbers reinforce the same conclusion from
   the other side: CMC's boosted reserves peak at only 7-24% of their 150 N·m budget (0.2% at mtp)
   during a SUCCESSFUL run, so there is headroom to spare there too — capacity was never the binding
   constraint once real muscles share the load. Two concrete, disclosed-but-untested candidates remain
   for RRA specifically (which has no muscles and so cannot spread the load the way CMC does), both
   already flagged as simplifications in `docs/MECHANISM_RRA.md` §7.10 and not touched by this session's
   fix: (a) the uniform CMC_Joint task gains (`kp=100, kv=20`, `weight=1` for all 33 coordinates) may be
   trading off ankle/subtalar/mtp tracking fidelity against other coordinates in a way a differentiated,
   per-joint weighting (e.g. heavier weight on ankle/subtalar) could change; (b) the ankle/subtalar/mtp
   desired kinematics themselves may carry genuine skin-marker IK noise/inconsistency at a magnitude no
   torque reallocation can track — a well-documented measurement-fidelity limitation for these specific
   rotations in marker-based gait analysis, independent of anything about this model's actuators.
   Neither is attempted here (time-boxed; identified, not vague).

## 7. Honest gaps (full list)

1. **The Hicks moment gate still fails and gets measurably worse**, not just "does not improve" — a
   real, disclosed regression on that specific metric, not hidden behind the force-gate's narrow pass.
2. **150-vs-300 saturation was measured only for the RRA leg** (cheap, ~1 min/run); it was not
   separately re-confirmed for the muscle-driven SO re-solve or CMC at boost=300 (skipped deliberately,
   lean: RRA's own numbers already show complete saturation, so re-running the ~10-17 CPU-minute CMC
   leg a second time at 300 N·m was judged very unlikely to show a different qualitative outcome than
   150 N·m — a reasoned scope cut, not an oversight, but not literally measured either).
3. **Intermediate boost values (e.g. 100, 200, 250 N·m) were not tested** — the 150-vs-300 bracket
   already saturates, so an interior value is very unlikely to behave differently, but this is
   reasoned, not measured.
4. **The mtp low-inertia/no-passive-restoring-force explanation (§3) is a plausible, disclosed
   hypothesis, not independently proven** — the coordinate-limit-violation alternative was explicitly
   tested and ruled out, but the "cheap to perturb" account was not itself isolated (e.g. by adding a
   test passive stiffness and re-measuring).
5. **The task-gain (kp/kv/weight) hypothesis (§6) is identified, not tested** — this session's fix
   only touched actuator `optimal_force`, per the operator's specific request; re-running with
   differentiated task weights is the natural next experiment, not attempted here (time-boxed).
6. **Knee/hip contact-force estimates remain far from the OrthoLoad anchor** and are not meaningfully
   changed by this fix (§4) — consistent with, not independent of, `docs/MECHANISM_RRA.md`'s own
   finding that this is a whole-limb kinematic/dynamic issue, not fixable by this actuator change alone.
7. **Single trial, right-leg-primary** (subject2 `walking1`) — same scope caveat as every upstream
   cert in this family; no claim of generality across subjects/trials/gait speeds.
8. **CMC_Joint task weights remain uniform (1.0)** for all 33 coordinates — unchanged by this session,
   flagged in §6 as the leading honest-next-cause candidate.
9. **The CMC full-trial fix was only run/confirmed at boost=150 N·m**, not separately re-confirmed at
   300 (deliberately skipped, lean, per gap 2 above) — the successful run already shows comfortable
   reserve headroom (7-24%, §5), so a 300 N·m repeat is expected to also pass, but this is reasoned from
   the RRA-leg saturation evidence, not independently measured for CMC itself.
10. **CMC's own forward-dynamics accuracy/tracking-error against the real IK was not separately
    re-measured** the way RRA's was (§3) — this session's CMC question was specifically "does it
    converge," which is now answered cleanly; whether the resulting motion tracks the real gait as
    tightly as RRA's (mixed) result is a distinct, unasked, unmeasured question.
11. **The 6 pelvis-residual actuators remain heavily saturated during CMC** (213-1271% of their small
    10/15 N optimal_force, §5) — untouched by this fix (only ankle/subtalar/mtp were boosted, per the
    operator's specific request) and consistent with, not resolving, this whole cert family's
    established pelvis-residual/dynamic-consistency finding.

## Files

- `scripts/msk/ankle_reserve_fix.py` — the full, re-runnable, staged pipeline (`baseline_check`,
  `actuators`, `rra`, `tracking_diagnostic`, `so_jr`, `cmc_smoke`, `cmc_full`, `fast_all` stages;
  `--boost <N.m>` selects the magnitude). Builds both boosted actuator files fresh via the OpenSim API
  (never `.clone()`, which segfaulted in a pre-check, §2), reuses `run_rra.py`'s/`static_opt_knee.py`'s/
  `validate_hip_force.py`'s/`validate_emg_timing.py`'s own proven functions throughout (via direct
  reuse for path-parametrized functions, via a save/restore monkeypatch for the few module-level
  constants those functions read at call time — verified correct and side-effect-free before use).
- `data/msk_smoketest/subject2_walking1/ankle_reserve_fix/` — all outputs (NEW tree, the existing
  `rra/` and `contraction_dynamics/cmc/` trees from the weak-reserve run are untouched, preserving
  their evidence): `anchors.json`, `baseline_check_results.json`; `boost150/` (the primary, fully-run
  boost — `actuators/`, `rra/pass{1,2}/`, `tracking_diagnostic` (in `results.json`),
  `so_jr/{so_post,jr_post}/`, `cmc/{smoke_0p3s,full_1p57s}/`); `boost300/` (the OODA follow-up,
  RRA-leg only, deliberately not extended to `so_jr`/`cmc` — gap 2/9 above — `actuators/`,
  `rra/pass{1,2}/`, `tracking_diagnostic`); each with its own `results.json` — every number in this
  document is machine-written into one of these files, traceable back to the exact run that produced it.
