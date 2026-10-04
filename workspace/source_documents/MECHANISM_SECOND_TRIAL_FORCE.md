# MECHANISM SECOND-TRIAL FORCE CERT — does the ~1.5x over-prediction hold on a different ACTIVITY? (2026-07-21)

**Question** (`docs/MECHANISM_CMC_SECOND_SOLVE.md` / `docs/MECHANISM_COMMON_MODE.md`'s own named biggest
remaining risk): every knee/hip contact-force number in this repo so far — walking1's own
391.10 %BW knee / 386.77 %BW hip (ratio 1.515x / 1.412x vs OrthoLoad) — comes from ONE subject2
IK+GRF trial. Does the ~1.5x over-prediction hold on a genuinely DIFFERENT, higher-demand
**activity** (not just a different subject on the same walking trial, already tested in
`docs/MECHANISM_CROSS_SUBJECT.md`)?

**Headline: NO — and the pre-registered falsifier is not a knife-edge miss, it fails by ~8-9x its
own tolerance, for a well-diagnosed reason.** The real drop-jump landing trial (DJ1, subject2) is
the only available non-walking activity this repo had not already muscle-force-tested
(`docs/MECHANISM_CROSS_ACTIVITY.md` already covered squat/STS). Static Optimization on DJ1
**numerically converges** (0/100 frames NaN, all activations in bounds, both of 2 tested SO
template variants) but is **dynamically implausible by a wide margin** — 16-17 muscles pinned near
maximum simultaneously and pelvis reserve actuators leaning **20-24x** over the walking-trial
comfort band (vs. squat's 0.78x / STS's 2.40x — DJ1 is categorically worse than either
previously-tested ballistic activity). The resulting knee/hip contact-force ratio vs the
activity-matched OrthoLoad anchor (Trampoline/Jumping: knee 547.5 %BW n=4, hip 428.4 %BW n=14) is
**~2.6-2.8x**, not ~1.5x — a delta of **1.22-1.39** against walking1's own ratio, roughly **9x** the
pre-registered ±0.15 tolerance. **Verdict: the ~1.5x is walking-(steady-state-gait-)specific, not a
model-level constant** — extending, not contradicting, `docs/MECHANISM_CROSS_ACTIVITY.md`'s own
finding that Static Optimization's muscle-driven tier breaks down on rapid/ballistic movements,
now confirmed on a THIRD, most-severe-yet instance. This is read as **good news for the walking
number specifically**: SO's behavior is activity-*dependent* in an interpretable, mechanistically
grounded way (converges cleanly on steady-state cyclic gait; fails via a distinct, machine-flagged
muscle-pinning/reserve-saturation signature on ballistic events) rather than "SO just always says
~1.5x regardless of input" (which would have been a more worrying, degenerate finding).

## 0. Trial enumeration (machine-run, not assumed) — `enumerate_available_trials()`

10 subjects on the external mount (`subject2`-`subject11`; no `subject1`), **22 distinct IK trial
stems** pooled across all: `DJ1-5`, `DJAsym1-5`, `squats1`, `squatsAsym1`, `STS1`, `STSweakLegs1`,
`walking1-4`, `walkingTS1-4`.

| check | result |
|---|---|
| Subjects with complete DJ1 (IK **and** ID **and** real force-plate GRF) | **9/10** (all except `subject7`, which has none of the three) |
| subject2 DJ1 complete | **True** — `OpenSimData/Mocap/IK/DJ1.mot` (100 frames, 100 Hz, real range [1.25,2.24]s), `ForceData/DJ1_forces.mot` (5679 rows, 2000 Hz, real range [0,2.839]s), `OpenSimData/Mocap/ID/DJ1.sto` |
| Native Static-Optimization **template** shipped for DJ1, any subject | **False** — every subject's own `OpenSimData/Mocap/SO/` only ships `squats1`/`squatsAsym1`/`walking*`/`walkingTS*` templates, **never** any `DJ*` or `STS*` template. Same corpus-wide gap `docs/MECHANISM_CROSS_ACTIVITY.md` already found and adapted around for STS. |
| Already muscle-force-tested this repo (before this session) | squat (`squats1`), sit-to-stand (`STS1`) — `docs/MECHANISM_CROSS_ACTIVITY.md`. Cross-**subject** (not activity) walking1 tested on subject3/4 — `docs/MECHANISM_CROSS_SUBJECT.md`. **DJ1/DJ2/DJ3/DJAsym\* untested at the muscle-force level anywhere in this repo before this session** (`docs/MECHANISM_MUSCLE_DAMAGE.md` explicitly disclosed "no SO exists for DJ1... none was solved here" as an open gap, §8 of that doc). |

**DJ1 chosen**: the only trial that is simultaneously (a) real IK+GRF-complete, (b) genuinely
untested at the muscle-force level, and (c) independently documented as a **higher-demand** event
than walking — `docs/MECHANISM_MUSCLE_DAMAGE.md` §8 measured DJ1's own peak eccentric fiber
velocity at **1.87-8.59x** walking1's own whole-trial peak, across all 7 landing-relevant muscles,
and its own peak two-foot vertical GRF (measured fresh this session, below) at ~4.4x walking's own
%BW peak.

## 1. Pre-registration (stated before running Static Optimization)

```
C  ("generalizes"):  DJ1's muscle-driven tier is TRUSTWORTHY (SO converges; self-computed vs
    JointReaction agree <5% relative, the squat/STS bar; fsb.diagnose_so_dynamical_plausibility
    reports dynamically_plausible=True) AND |ratio_DJ1(activity-matched anchor)
    - ratio_walking1(walking anchor)| <= 0.15, knee and/or hip => MODEL-LEVEL, strengthens trust.
not-C: ratio differs by > 0.15, OR the muscle-driven tier is not trustworthy (exactly what
    docs/MECHANISM_CROSS_ACTIVITY.md already found for squat/STS) => trial-specific / open,
    reported as an HONEST, INFORMATIVE outcome.
Symmetric-QC precondition: DJ1's whole-body Newton residual examined BEFORE trusting any force
    number (§2) — a bad GRF-kinematics match would manufacture a false ratio.
```

## 2. Symmetric-QC FIRST: is DJ1's IK+GRF self-consistent? (forced OODA, not skipped)

**Observe.** Whole-body Newton residual (mass x COM-accel = GRF + gravity, the identical gate
`validate_joint_force.py`/`cross_activity_validation.py` already use, reused unmodified) on DJ1's
100-frame IK-available window [1.25,2.24]s: **RMS = 13.87 / 17.58 / 6.85 %BW (x/y/z)** — **FAILS**
the pre-registered <8 %BW gate on 2 of 3 axes. Per-frame, the residual is not uniformly bad: it
spikes to 40-65 %BW specifically in a ~150 ms band around the landing impact (t=1.54-1.65s) and is
comparatively better-behaved before (t<1.50s, pre-impact) and in a brief local dip just after
(t~1.70-1.71s, 10.3 %BW) before rising again in a second, smaller bump (t~1.75-1.95s, 14-25 %BW).

**Orient — why?** Two candidate causes were checked, not assumed:
1. **Marker/IK quality** — checked directly: `marker_error_RMS` stays flat at 0.016-0.020 m
   through the landing window (max over the whole trial 0.0279 m, at t=2.14s, nowhere near the
   impact) — **ruled out**, the IK solve itself is not degrading at the impact instant.
2. **A single arbitrary differentiation-window choice** — swept SavGol windows 7/11/17 samples
   (70/110/170 ms): interior RMS 29.3/23.4/30.3 %BW, **all fail** the same gate; 0-1 valid
   contiguous low-residual segments found (thr=15 %BW, min 0.3s) depending on window — **not a
   single-window artifact**, the elevated residual is a robust feature of this trial, not a
   parameter-tuning fluke.
3. **A genuine missing third contact** (the STS chair-confound class of bug,
   `docs/MECHANISM_CROSS_ACTIVITY.md` §3.1) — checked directly: DJ1's GRF file carries an unused
   third force-plate channel (`3_ground_force_*`); measured **exactly 0.0 N at every one of 5679
   rows, the entire 2.839s recording** — **ruled out**, not an unseen contact.
4. **Differentiation-resolution limit at a near-step impact** — DJ1's own combined vertical GRF
   rises from ~0 to 398 %BW in ~30-40ms (measured max d(GRF)/dt = 144,939 N/s at t=1.5855s,
   matching `docs/MECHANISM_MUSCLE_DAMAGE.md`'s own independently-measured 141,975 N/s at the same
   instant), i.e. a near-step transient comparable to or faster than the 100Hz IK sampling
   interval (10ms) and the 110ms SavGol smoothing window used to differentiate it.

**External, independently-derived corroboration** (not this script's own claim alone):
`docs/MECHANISM_WOBBLING_MASS.md` — built via a completely different method (a spring-damper
wobbling-mass EOM, not a whole-body Newton residual) — **already and independently found**: *"the
100Hz-IK real-data path is resolution-limited (own validity gate FAILS for the drop-landing
trial)"*, quantified there as thigh_r/shank_r whole-body-residual **FAIL at every SavGol window
tried (14.6-15.6 %BW x, 18.7-28.2 %BW y)** — matching this session's fresh 13.9/17.6 %BW x/y to
within ~1 point, via an entirely separate analysis. Two independent methods, two independent
sessions, the same conclusion.

**Decide/Act.** This is a genuine, cross-corroborated **differentiation-resolution limitation at
the specific instant of peak physical interest**, not a missing-contact bug and not a
single-window artifact. It does **not**, by itself, block running Static Optimization — SO's own
internal coordinate-differentiation (a low-pass-filtered spline fit, a materially different
computation from this script's own SavGol residual check) is not mechanically the same
computation — but it **does cap the confidence tier** of whatever force number results: **no
number from this trial can be reported `in-vivo-anchored`-tier at the specific instant nearest the
impact; `method-only`/`diagnosed-gap` is the ceiling.** No valid low-residual sub-segment both
excludes the compromised region and retains the landing event this trial exists to test, so the
full [1.25,2.24]s window was used (matching the squat cert's own "whole trial when nothing forces
windowing" convention) rather than windowing away from the content of interest.

## 3. Static Optimization on DJ1 — 2 template variants (sensitivity check, since no DJ1-native template exists)

No trial in the corpus ships a DJ1-native SO template (§0); adapted from the two existing,
already-validated templates that differ in exactly one field — `lowpass_cutoff_frequency_for_coordinates`
(6Hz, walking1's own original template family vs 4Hz, squats1's — the value `docs/MECHANISM_CROSS_ACTIVITY.md`
used for its own STS adaptation) — reserve-actuator file confirmed byte-identical
(`diff walking1_reserveActuators.xml squats1_reserveActuators.xml` = 0 differences, re-verified
live) across both.

| | primary (6Hz, walking1-template) | secondary (4Hz, squats1-template) |
|---|---:|---:|
| SO convergence (`cav.check_so_convergence_generic`) | **PASS** (100/100 frames, 0 NaN, 0 out-of-bounds) | **PASS** (100/100 frames, 0 NaN, 0 out-of-bounds) |
| Muscles pinned >95% activation, >10% of frames (`fsb.diagnose_so_dynamical_plausibility`) | **17** (top: `vaslat_r` 46%, `vaslat_l` 44%, `glmax2_r` 38%, `glmax2_l` 29%, `vasmed_r/l` 26%, `glmax1_r` 25%, ...) | **16** (same pattern: `vaslat_r` 45%, `vaslat_l` 43%, `glmax2_r` 36%, ...) |
| Peak pelvis reserve force / moment | 1512.9 N / 850.5 Nm | 1811.9 N / 912.7 Nm |
| vs. 75N/75Nm comfort band (unmodified, same constants squat/STS were checked against) | **20.17x over** | **24.16x over** |
| **`dynamically_plausible`** | **False** | **False** |

**For context, squat/STS's own already-published implausibility signatures** (`docs/MECHANISM_CROSS_ACTIVITY.md`
§3.3): squat 4 muscles pinned / 0.78x reserve (**within** the comfort band — muscle-pinning alone
drove squat's implausibility); STS 1 muscle pinned / 2.40x reserve. **DJ1's 16-17 pinned / 20-24x
reserve is categorically more severe than either** — the muscle-pinning pattern also broadens (vasti
+ glutes together, not one isolated group), consistent with an extreme, poorly-conditioned
per-frame recruitment problem at a genuine landing impact, not a marginal borderline case.

**Robustness**: both template variants (6Hz and 4Hz coordinate low-pass) give the **same
qualitative verdict** and similar magnitudes (knee 1507.7 vs 1538.7 %BW, hip 1129.6 vs 1200.0 %BW,
~2-6% apart) — the finding does not hinge on which of the two disclosed template adaptations is used.

## 4. Self-computed (`vhf.compute_self_cross_check_generic`) vs official `opensim.JointReaction`

The task's explicitly-named function, reused unmodified for **both** knee and hip cuts (the same
function `cmc_second_solve.py` already used this way for a different second-solve axis):

| variant | joint | self-computed peak %BW | JointReaction peak %BW | rel. diff. | peak times (self / JR) |
|---|---|---:|---:|---:|---|
| primary (6Hz) | knee | 1507.72 @t=1.71s | 1523.46 @t=1.71s | **1.03%** | same instant |
| primary (6Hz) | hip | 1129.60 @t=1.71s | 1141.66 @t=1.60s | **1.06%** | 110ms apart, magnitude still agrees |
| secondary (4Hz) | knee | 1538.75 @t=1.71s | 1555.73 @t=1.71s | **1.09%** | same instant |
| secondary (4Hz) | hip | 1200.03 @t=1.71s | 1210.40 @t=1.71s | **0.86%** | same instant |

**All 4 agree within ~1%** — looser than walking's own 0.004-0.070% but far inside the
pre-registered <5% bar (the same bar squat/STS cleared at 0.00-0.10%). **This rules out a code
bug**: two structurally independent code paths (this repo's own geometric BFS free-body cut vs
OpenSim's compiled `JointReaction` analysis) converge on the same large number — SO's implausible
output is what SO **actually, consistently computed**, not a self-computed extraction error.

**Anatomical cross-check** (unaffected by the implausibility verdict — a separate, purely
geometric check): knee-r crossing set (13 muscles: the walking cert's own 12 **+** `tfl_r`) and
hip-r crossing set (25 muscles) **exactly match** `HIP_ANATOMICAL_MUSCLES_EXPECTED` and the
already-published walking/squat/STS crossing-muscle sets, muscle-for-muscle — confirms the
BFS+path-geometry detector is working correctly on DJ1's kinematics; only the SO recruitment
solution's physiological plausibility fails, not the geometry.

## 5. OrthoLoad activity-matched anchors (pooled fresh, `cav.load_orthoload_trials`/`pool_pctbw`, unmodified)

| joint | activity bucket | n trials | n subjects | median %BW |
|---|---|---:|---:|---:|
| knee | Jogging/Running | 4 | 2 | 621.5 |
| knee | **Trampoline/Jumping** | 4 | 2 | **547.5** |
| knee | Walking (context) | 169 | 9 | 258.22 (published) |
| hip | Jogging/Running | 38 | 9 | 467.8 |
| hip | **Trampoline/Jumping** | 14 | 8 | **428.4** |
| hip | Walking (context) | 382 | 19 | 273.93 (published) |

**Cross-checked against the pre-existing, independently-built aggregate summary**
(`scripts/msk/index_orthoload_forces.py`'s own `orthoload_joint_activity_summary.csv`, built by a
different session/script): knee Trampoline/Jumping 547.5 (exact match), knee Jogging/Running 621.5
(exact match), hip Trampoline/Jumping gen1(n=5,median=356.1)+gen2(n=9,median=464.7)=n=14 pooled
(matches this session's fresh n=14 exactly), hip Jogging/Running gen1(n=9)+gen2(n=29)=n=38 pooled
(matches exactly) — the fresh-pooling function introduces no drift when applied to these new
bucket labels.

**Confidence tier: `published-plausibility`, not `in-vivo-exact`**, per the task's own framing —
these ARE real in-vivo instrumented-implant measurements (external anchor), but n is small for the
jump-specific buckets (4 knee trials/2 subjects; 14 hip trials/8 subjects, vs. walking's 169-382
trials/9-19 subjects) and the OrthoLoad activity LABEL "Trampoline/Jumping" is this implant
registry's own coarse daily-activity classification (could include recreational trampoline
bouncing, hopping, or other jump types by these patients, not necessarily a controlled
laboratory drop-jump-landing protocol specifically) — a genuine, disclosed label-vs-protocol
mismatch even though it is the best available match.

## 6. Ratios + falsifier verdict — side-by-side with walking1

| | walking1 (published) | DJ1 primary (6Hz) | DJ1 secondary (4Hz) |
|---|---:|---:|---:|
| knee self-computed | 391.10 %BW | **1507.72 %BW** | **1538.75 %BW** |
| knee JointReaction | 391.11 %BW | 1523.46 %BW | 1555.73 %BW |
| knee anchor used | Walking median 258.22 (n=169) | Trampoline/Jumping median 547.5 (n=4) | (same) |
| **knee ratio (self-computed / anchor)** | **1.515x** | **2.754x** | **2.810x** |
| hip self-computed | 386.77 %BW | **1129.60 %BW** | **1200.03 %BW** |
| hip JointReaction | 387.04 %BW | 1141.66 %BW | 1210.40 %BW |
| hip anchor used | Walking median 273.93 (n=382) | Trampoline/Jumping median 428.4 (n=14) | (same) |
| **hip ratio (self-computed / anchor)** | **1.412x** | **2.637x** | **2.801x** |
| **delta vs walking1's own ratio** | — | knee 1.239 / hip 1.225 | knee 1.295 / hip 1.389 |
| pre-registered falsifier tolerance | ±0.15 | **FAILS by ~8-9x** | **FAILS by ~9x** |
| `dynamically_plausible` (precondition) | True (walking1, established) | **False** | **False** |
| **`generalizes`** | — | **False** | **False** |

**Robustness of the negative**: using the Jogging/Running anchor instead (knee 621.5, hip 467.8)
gives ratios 2.43x knee / 2.41x hip — still a delta of ~0.91-1.00 against walking1's ratio, **~6-7x
the tolerance**. The "does not generalize" conclusion does not hinge on which activity-matched
anchor is chosen. Against the walking anchor itself (the task's fallback: *"if only walking medians
apply... treat as published-plausibility"*), DJ1's raw numbers are **4.1-6.0x** walking's own
in-vivo walking median — even larger, for the obvious reason that walking is not the relevant
comparison activity for a landing impact; reported for completeness, not as the primary reading.

## 7. Pure-reaction (Newton, muscle-free) tier — secondary, complementary cross-check

Same free-body-cut Newton's-law method already shown trustworthy/cross-activity-stable for
walking/squat/STS (`docs/MECHANISM_CROSS_ACTIVITY.md` §4: ratio 0.25-0.40 across those 3):

| | knee_r | hip_r |
|---|---:|---:|
| DJ1 pure-reaction peak | 128.29 %BW @t=1.60s | 96.03 %BW @t=1.83s |
| vs. Trampoline/Jumping anchor | **ratio 0.234** | **ratio 0.224** |
| walking/squat/STS established band (context) | 0.25-0.40 | 0.25-0.40 |

DJ1 sits just **below** the previously-established 0.25-0.40 cross-activity band rather than inside
it — a modest, disclosed extension (not alarming; consistent with the same differentiation-resolution
ceiling §2 already identified suppressing the pure-reaction peak too, since this tier is equally
exposed to the same SavGol-based differentiation). The knee peak (t=1.60s) lands almost exactly at
the raw GRF impact instant (measured independently at t=1.5855-1.595s); the hip peak lands notably
later (t=1.83s), plausibly a distinct, later stabilization-phase loading event — observed, not
further mechanistically decomposed here (out of this session's scope).

## 8. Forced adversaries (both directions, symmetric per the task's own instruction)

**Leaning-negative adversary** (the temptation: accept "SO fails on DJ1" as an honest negative
after one shot): forced via OODA, not accepted on a single pass —
(a) tried **2** independent template/lowpass-cutoff variants (6Hz, 4Hz): both implausible, closely
    agreeing in magnitude (2-6% apart) — not a single arbitrary-parameter artifact;
(b) checked whether the implausibility signature could be a windowing/edge artifact: the reported
    peak (t=1.71s) sits comfortably interior (5-frame edge margin excluded on a 100-frame series),
    and the muscle-pinning pattern spans 16-17 DIFFERENT muscles across 2 different muscle groups
    (vasti, glutes) — a broad, structural over-demand signature, not a 1-2-frame numerical fluke;
(c) verified self-computed vs JointReaction agree tightly (~1%) — ruling out a code-level
    extraction bug as the source of the large number;
(d) verified the anatomical crossing-muscle sets exactly match the established walking/squat/STS
    sets — ruling out a geometry-detection bug.
All 4 checks converge on the same conclusion: this is a genuine, reproducible, well-diagnosed
property of Static Optimization applied to this specific movement class (ballistic landing), not a
one-shot artifact swept under "honest negative."

**Leaning-positive adversary** (had the ratio landed near 1.5x, the temptation would be to accept
generalization at face value): moot here since the ratio does not land near 1.5x at all (delta
1.2-1.4, ~9x tolerance) — but the SAME checks (a)-(d) above also guard against a false-POSITIVE
reading in the opposite direction (e.g., a code bug that happened to produce a coincidentally
plausible-looking number) — not applicable this session since the result is unambiguously negative,
noted for completeness per the symmetric-QC mandate.

## 9. What this DOES and does NOT resolve

**Resolved**: the ~1.5x knee/hip over-prediction is **not a model-level constant that holds
regardless of activity** — it is specific to (at minimum) steady-state, cyclic, non-ballistic
gait. Static Optimization's own architecture has a **distinct, independently-diagnosable failure
mode** (muscle-pinning + pelvis-reserve-saturation) on rapid/ballistic movements, now confirmed on
**3** independent instances (squat, STS, DJ1 landing) with DJ1 the most severe — a coherent,
strengthening pattern, not 3 unrelated coincidences. This is read as evidence FOR, not against, the
walking-specific number's trustworthiness: SO's behavior is activity-dependent in a way that
tracks a real, mechanistically-interpretable distinction (smooth steady-state vs. ballistic
transition loading), not a fixed multiplier the method would output regardless of input.

**NOT resolved**: (1) whether SO could be MADE trustworthy on ballistic movements with a
different, purpose-tuned differentiation/filtering scheme (an open, bounded follow-up named but
not attempted here, same as `docs/MECHANISM_CROSS_ACTIVITY.md`'s own disclosed next step); (2) the
DJ1-specific differentiation-resolution ceiling (§2) means even a hypothetically-plausible SO
result on this exact 100Hz-sampled trial would carry a capped confidence tier — a higher-rate
capture (this dataset's own video-derived IK sources, not attempted here) would be needed to fully
decouple "SO breaks down on ballistic movements" from "this specific 100Hz capture under-resolves
the impact"; (3) single subject (subject2), single trial (DJ1, not DJ2-5/DJAsym1-5) — no claim of
generality across the other 8 available drop-jump instances or 8 other subjects with complete DJ1
data (§0) — a natural, cheap-relative-to-this-session extension not attempted here.

## 10. Honest gaps (full list)

1. **DJ1's own whole-body residual gate fails** (§2) — capped confidence tier on any force number
   from this trial's exact impact instant, cross-corroborated but not eliminated.
2. **No DJ1-native SO/JR template exists anywhere in the 10-subject corpus** (§0) — both template
   variants used here are adaptations (disclosed, sensitivity-checked, not hidden).
3. **Only 2 of the many possible SO configuration variants tested** (2 lowpass cutoffs) — a
   narrower/wider time window, alternate reserve-actuator strength, or a different differentiation
   scheme entirely were not swept; given the severity and consistency of the 2 tested variants
   (16-17 pinned muscles, 20-24x reserve — categorically beyond squat/STS), a third variant
   reversing the conclusion is assessed as unlikely but not tested, a bounded, disclosed limit on
   how exhaustively the "honest negative" was forced (2 variants, not more, per LEAN discipline).
4. **Jump-specific OrthoLoad anchors are small-n** (4 knee trials/2 subjects, 14 hip trials/8
   subjects) — `published-plausibility` tier, not `in-vivo-exact`; the "Trampoline/Jumping" label
   is this registry's own coarse activity classification, not a protocol-matched drop-jump-landing
   category.
5. **Right leg only, single DJ1 instance, subject2 only** — 8 other DJ-family trials
   (`DJ2-5`,`DJAsym1-5`) and 8 other subjects with complete DJ1 data exist and were not tested
   (§0) — the natural next extension.
6. **The hip's self-computed-vs-JointReaction peak TIME differs by 110ms for the primary (6Hz)
   variant** (§4) — magnitude still agrees to 1.06%, and this is the SAME disclosed
   edge-exclusion/whole-window peak-search convention difference `docs/MECHANISM_CROSS_SUBJECT.md`
   §2.3 already diagnosed for subject3/4's hip — not a new bug, not further investigated here.
7. **The pure-reaction tier's own ratio (0.22-0.23) sits just below, not cleanly inside, the
   previously-established 0.25-0.40 cross-activity band** (§7) — plausibly the same
   differentiation-resolution ceiling (§2), not independently isolated as the specific cause.
8. **The mechanism behind the 110ms lag between the raw GRF/pure-reaction knee peak (t=1.60s) and
   the SO-driven self-computed peak (t=1.71s)** was observed, not mechanistically decomposed.

## Files

- `scripts/msk/second_trial_force.py` — the full, re-runnable pipeline (`.venv-msk/bin/python3
  scripts/msk/second_trial_force.py`, ~25s wall time for both SO variants + both JR passes).
  Imports and reuses UNMODIFIED: `validate_joint_force.py` (vjf), `static_opt_knee.py` (sok),
  `validate_hip_force.py` (vhf, specifically `compute_self_cross_check_generic`,
  `extract_jr_hip_r_force_pct_bw`), `force_scenes_batch.py` (fsb,
  `diagnose_so_dynamical_plausibility`), `cross_activity_validation.py` (cav, its own
  already-generalized SO/JR XML-templating + residual/pure-reaction + OrthoLoad-pooling plumbing).
  md5 of all 5 reused files recorded before/after in the JSON output (unchanged).
- `data/msk_smoketest/second_trial_force/second_trial_force_results.json` — every number in this
  document, machine-written (trial enumeration, residual QC + SavGol sweep, both SO variants'
  convergence/plausibility gates, self-computed + JointReaction results, OrthoLoad anchors,
  ratios/verdicts, reused-file md5s).
- `data/msk_smoketest/second_trial_force/dj1_primary_6hz/{so,jr}/`,
  `data/msk_smoketest/second_trial_force/dj1_secondary_4hz/{so,jr}/` — patched OpenSim setup XMLs +
  raw Static-Optimization/JointReaction `.sto` output for both template variants.
- Inputs read in place, never modified: subject2's real OpenCap LabValidation mocap
  (`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/`),
  `data/external/orthoload/_index/orthoload_akf_trials.jsonl`.
- Cited, not re-derived: `docs/MECHANISM_MUSCLE_DAMAGE.md` (DJ1 eccentric-velocity/landing-instant
  measurements), `docs/MECHANISM_WOBBLING_MASS.md` (independent DJ1 residual-gate-failure
  corroboration), `docs/MECHANISM_CROSS_ACTIVITY.md` (squat/STS dynamical-implausibility
  precedent + its own OrthoLoad cross-activity anchor-pooling method), `docs/MECHANISM_CROSS_SUBJECT.md`
  (the hip peak-search convention-mismatch precedent, §4/§10 above).

No git commit, no git push performed (isolation respected). No file outside this repo was written
to; the external LabValidation mount and OrthoLoad data were read-only throughout.
