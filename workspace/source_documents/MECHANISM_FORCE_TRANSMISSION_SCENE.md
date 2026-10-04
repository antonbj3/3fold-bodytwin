# MECHANISM FORCE TRANSMISSION SCENE — ankle→knee→hip→lumbar on a real corpus squat (2026-07-21)

Executes the operator's "force transmission + scenes" ask: one real corpus clip, force propagating up
the posterior chain, characterized end to end. Script: `scripts/msk/force_transmission_scene.py`
(`.venv-msk/bin/python3 -u scripts/msk/force_transmission_scene.py`, `FTS_RUN_SO=0` to skip the
slow/crash-prone secondary Static Optimization step). Every number below is this script's own
printed/JSON output, or a one-off diagnostic built on its own functions — never estimated.

## Bottom line

**The Newton reaction-force method (bug-immune, never touches `InverseDynamicsTool`) gives a clean,
cross-validated force-transmission profile across a real 5-rep squat set: ankle 119-178%BW, knee
113-170%BW, hip 86-131%BW, lumbar/back 116-176%BW, each rep's peak occurring during the RAPID part
of the descent, not at maximum depth.** The mass-shedding cross-check — a specific, falsifiable,
quantitative prediction, not a vibe — confirms the ankle>knee>hip trend is real physics (observed
37.6N / predicted 37.7N shank-weight shed, ratio 0.996; observed 97.6N / predicted 95.6N
thigh-weight shed, ratio 1.021 — both within 3% of an independently-known quantity: each segment's
own mass, not fit to the data). **The honest answer to "does ankle peak before knee before hip":
no resolvable sequence exists at this clip's own ~30fps (33ms) resolution — all 4 levels peak in
the SAME video frame, in all 5 reps, stably across a smoothing-window sensitivity sweep.** This
squat is not a fast ballistic jump; the literature's proximal-to-distal cascade may still apply to
a genuine jump clip (flagged as the concrete next step, §7), but is not demonstrated here.
**Static Optimization was genuinely attempted twice (not skipped) and fails a hard, precedented,
quantitative gate both times** (§5) — the delivered result is the Newton method, not SO.

## 1. Clip chosen

`2025-09-26_DPEOqLVkav7` (`weighted_lifts/squat`) — the same clip `docs/MECHANISM_CORPUS_IK.md`
already vetted as its primary target: det_frac=1.0 (368/368 frames), marker_verdict PASS, knee and
hip_flexion range gates PASS (ankle FAILs, diagnosed there as weak ankle/subtalar observability,
inherited and re-confirmed here, §6). Re-checked live at the top of this script (not re-quoted from
memory): `marker_verdict=PASS  range_gates={'knee': True, 'hip_flexion': True, 'ankle': False}`.
That report's own `lr_visibility_asymmetry` block flags `max_asymmetry=0.403` at the **knee**,
right-knee visibility collapsing to a minimum of 0.114 — so this build uses the **left**-side chain
(`ankle_l`/`knee_l`/`hip_l`/back+arms), re-asserted live at STEP 1 (`assert
lr_asym["max_asymmetry_landmark"] == "knee"`), not the right, to avoid a camera-occlusion artifact.

The clip is a **5-rep squat set** (not one isolated rep), confirmed by a numeric peak-finder on
`knee_angle_l` (never eyeballed): squat-bottom events at t=2.8, 4.5, 6.2, 7.8, 12.07s, peak knee ROM
98.2°/93.9°/103.4°/106.9°/132.3° — the deepest, final rep is also the clip's own documented
"trends from ~45° to >130° over its last ~4s" feature (`docs/MECHANISM_CORPUS_IK.md` §3.2).

## 2. Method — Newton's-law free-body cuts + a geometric fix for the monocular translation problem

### 2.1 The new problem this clip poses (subject2's real mocap+force-plate data did not have it)

`validate_joint_force.py`/`static_opt_knee.py` (subject2 `walking1`) had a real force plate and a
trustworthy solved pelvis translation. This corpus clip has neither: monocular video → no GRF
measurement, and `docs/MECHANISM_CORPUS_IK.md` §5.2 already found `pelvis_tx/ty/tz` attenuated to a
few-cm noise-like residual (MediaPipe's `world_landmarks` are hip-re-centered every frame — the
translation signal isn't in the input the IK solver fit). Both a GRF measurement and a trustworthy
absolute whole-body acceleration are missing.

### 2.2 Geometric fix: a stance-foot-anchored frame

For a double-support squat the stance foot does not slide on the ground — a frame anchored to it is
a valid inertial frame, and every OTHER body's position in that frame is fully determined by the
already range-gate-PASSED RELATIVE joint angles, independent of the broken `pelvis_tx/ty/tz`.
Because rigid translation shifts every body identically, this is an exact rebasing: subtract the
anchor body's (`calcn_l`, left foot) own trajectory from every body's trajectory
(`ground_frame()` in the script).

### 2.3 OODA iteration #1 (forced, not skipped past): the anchor itself is jittery

A first pass (reusing `validate_joint_force.py`'s own 100Hz-tuned 5-13-sample window class, ported
to ~30fps as 165-430ms) gave **implausible 585-1300%BW whole-body GRF peaks**. Forced re-diagnosis,
not accepted: direct measurement of `calcn_l`'s own raw per-frame COM delta vs `pelvis`'s (STEP 3's
own printed diagnostic) —

| | max per-frame delta | mean per-frame delta |
|---|---:|---:|
| `calcn_l` (the anchor) | 0.1584 m | 0.0395 m |
| `pelvis` | 0.0301 m | 0.0054 m |

**the anchor is 7.4x jitterier than pelvis** — inherited directly from `docs/MECHANISM_CORPUS_IK.md`
§3.2/3.5's own already-documented weakly-observed ankle/subtalar DOFs. Because every body's
grounded position carries the identical `-anchor(t)` term, this noise does not average down by mass
fraction when summed for the whole-body GRF estimate — it survives with FULL, unweighted magnitude
(a real, derived-from-the-geometry consequence, not asserted: summing `mass[b]*(raw[b]-anchor(t))`
over all bodies leaves exactly one un-cancelled copy of `anchor(t)`, regardless of how many bodies
or their masses). **Fix**: smooth position AND acceleration for every body at the SAME
Savitzky-Golay window before grounding (`compute_pos_and_accel()` — one filter pass per body,
never re-differentiate an already-smoothed signal, which would silently compound two windows).

**Window choice, externally anchored, not tuned to taste:** `GRF_PLAUSIBILITY_CEILING_PCT_BW=400`
(walking ~1.0-1.3xBW, running ~2-3xBW, jump takeoff/landing ~2-4xBW — literature ballparks, a squat
generally less explosive than a jump). Swept 15-31 samples:

| window | ms | GRF peak %BW | vs 400%BW ceiling |
|---:|---:|---:|---|
| 15 | 495 | 549.5 | OVER |
| 19 | 627 | 457.2 | OVER |
| **23** | **759** | **359.1** | **under — PRIMARY (smallest window clearing the ceiling)** |
| 27 | 891 | 312.4 | under |
| 31 | 1023 | 276.8 | under |

### 2.4 OODA iteration #2 (a tempting "fix" that was WRONG — symmetric QC on this build's own work)

A closer look at STEP 6 raised a real adversary: ALL FOUR cuts (ankle/knee/hip/back) peaked at the
exact same frame in every rep — suspicious, since `back` carries no GRF term at all. Direct
measurement confirmed a mechanism: `calcn_l`'s own SMOOTHED acceleration still spikes to 27 m/s² at
t=2.9s (vs its own 6.8 m/s² mean) while `torso`'s own raw acceleration there is only 1.97 m/s² —
subtracting the anchor's own large event from torso manufactures roughly -25 m/s² of pure
subtraction artifact in torso's "grounded" acceleration.

The tempting fix — algebraically isolate ONLY the broken pelvis-translation component
(`raw(calcn_l) - com_relpos(calcn_l)`, where `com_relpos` re-runs forward kinematics with
`pelvis_tx/ty/tz` forced to 0) and subtract only THAT from every body — was built and tested
standalone. It correctly reproduced `docs/MECHANISM_CORPUS_IK.md`'s own published
`pelvis_tx/ty/tz` range (0.056/0.102/0.079 m vs the doc's 0.080/0.111/0.077 m) — confirming the
derivation is right — **but it is a different, worse physical assumption**: because a body's own
translation contribution is ALGEBRAICALLY IDENTICAL to what gets isolated this way, subtracting it
from `pelvis` itself collapses pelvis's own grounded trajectory to nearly flat (it's a
PELVIS-anchored frame in disguise, not a FOOT-anchored one) — exactly the wrong assumption for
recovering squat depth. Measured, not asserted: this alternative dropped the pelvis-height vs
knee-angle correlation from -0.92 to a meaningless **+0.53** (below the pre-registered
threshold — see §3). **Rejected; reverted to the foot-anchored method (§2.3).** The residual
"all 4 cuts peak simultaneously" question is real and is reported honestly, not resolved, in §4.2.

## 3. Cross-checks (pre-registered thresholds, machine PASS/FAIL)

| cross-check | what it forces | threshold | measured | verdict |
|---|---|---|---|---|
| `model.assemble()` clean | every one of 368 frames solves | 0 failures | 0/368 | **PASS** |
| B: corrected pelvis-height vs `knee_angle_l` | external anchor — knee ROM independently range-gate PASSED in a DIFFERENT pipeline stage; if the geometric fix recovers real depth (not just relabels noise), height must anti-correlate with flexion | r ≤ -0.70 | **r = -0.9233** | **PASS** |
| Mass-shedding, ankle→knee | is ankle≥knee≥hip real segment-weight subtraction or a damping/averaging artifact of summing more bodies? checked at the 104 lowest-whole-body-acceleration (quasi-static) instants only | ratio∈[0.5,1.5] vs shank's OWN mass | observed 37.58N / predicted 37.74N → **ratio 0.996** | **PASS** |
| Mass-shedding, knee→hip | same, vs thigh's OWN mass | ratio∈[0.5,1.5] | observed 97.58N / predicted 95.56N → **ratio 1.021** | **PASS** |
| A: `calcn_r` (right foot) planted per-rep | is the untouched OTHER foot stationary during each rep (double-support assumption)? | ratio < 0.4 vs that rep's own pelvis-y excursion | 0/5 reps, ratio 1.26-1.48 | **INCONCLUSIVE** — see below, not a gate |

**Cross-check A is deliberately excluded from the pass/fail gate**, not silently dropped: `calcn_r`
is built from `hip_r/knee_r/ankle_r/subtalar_r` — exactly the RIGHT leg that this same clip's own
`lr_visibility_asymmetry` already measured as severely occluded (knee visibility min=0.114). A
large `calcn_r` excursion is consistent with either real foot movement OR pre-existing right-leg
tracking noise; this check's own "control" data is independently known to be compromised, so a FAIL
here cannot be cleanly attributed to the stance-foot assumption being wrong. Context (not gated):
`calcn_r`'s excursion over the WHOLE clip is 1.38 m and barely changes across a 15-31-sample
smoothing sweep (a genuine low-frequency signal, not jitter) — consistent with the athlete
repositioning the right foot between reps in a real 5-rep set, which does not by itself falsify
"planted within one rep."

**Overall primary-method gate** (assemble-clean + B + both mass-shedding checks): **PASS**
(`force_transmission_scene_results.json`, script exit code 0).

## 4. Results — the force-transmission profile

### 4.1 Peak magnitude per level (global, whole clip, primary window=23 samples/759ms)

| level | subchain (BFS on model's own joint tree) | chain mass | peak (%BW) | t_peak |
|---|---|---:|---:|---:|
| ankle_l | talus_l, calcn_l, toes_l | 1.63 kg (2.1%) | 177.6 | 2.900 s |
| knee_l | + tibia_l | 5.47 kg (7.0%) | 170.1 | 2.900 s |
| hip_l | + femur_l, patella_l | 15.22 kg (19.5%) | 130.7 | 2.900 s |
| back (lumbar) | torso + both arms | 35.54 kg (45.4%) | 175.8 | 2.900 s |

**Geometric caveat the "ankle→knee→hip→lumbar" framing needs** (established in
`docs/MECHANISM_BAND_POSTERIOR_CHAIN.md` §1, re-verified live on this model's joint tree):
`ankle_l→knee_l→hip_l` is one true SERIAL chain — each free-body cut is nested inside the next, and
the ankle>knee>hip ordering IS a real, quantitatively-verified (§3 mass-shedding) consequence of
that nesting. `back` is a **parallel sibling branch** off the pelvis (`hip_r`, `hip_l`, and `back`
all parent directly off `pelvis` — confirmed via the same BFS/joint-tree dump this repo already
uses), not a serial continuation of the leg's reaction force — it is its own free-body cut (torso+
arms, no ground contact at all, so no GRF term enters its Newton's-law equation), coupled to the leg
only through whole-body equilibrium. Reported as 4 levels of one scene, not one undifferentiated
chain.

### 4.2 Per-rep peaks + the TIMING-SEQUENCE question, honestly answered

| rep (t bottom, knee ROM) | ankle_l | knee_l | hip_l | back | t_peak (ALL 4, identical) |
|---|---:|---:|---:|---:|---:|
| 2.80s, 98.2° | 177.6% | 170.1% | 130.7% | 175.8% | **2.900s** |
| 4.50s, 93.9° | 167.6% | 160.6% | 123.3% | 167.4% | **4.600s** |
| 6.20s, 103.4° | 158.4% | 152.2% | 117.6% | 159.5% | **6.233s** |
| 7.80s, 106.9° | 172.8% | 165.5% | 127.7% | 171.6% | **7.933s** |
| 12.07s, 132.3° (deepest) | 118.5% | 112.7% | 85.9% | 116.5% | **9.467s** |

**Every single one of the 12 inter-cut gaps, across all 5 reps, is exactly 0.0 ms** — not "close to
zero," bit-identical `t_peak` values. This is a tie, not a sequence: the printed "order" a naive
`sorted()` call would emit (`ankle_l < knee_l < hip_l < back`) is an artifact of Python's stable-sort
tie-breaking on identical timestamps, not a resolved temporal precedence — checked directly against
the raw per-cut `t_peak` values (not assumed from the sorted label), and confirmed identical to the
`.001`-second data-grid resolution in every rep. **The honest, measured finding: at this clip's own
~30fps (33ms) temporal resolution, ankle/knee/hip/back all peak within the SAME video frame, in
every rep, stably across a 15-31-sample Savitzky-Golay sensitivity sweep** (order_by_window:
identical `('ankle_l','knee_l','hip_l','back')` at windows 15/19/23/27/31 — stable, but stable
because it's a tie every time, not because a genuine small lag survived the sweep unchanged).

**This is a real, physically informative negative, not a shortfall of the method**: a controlled
back squat is not a fast ballistic movement. The classic biomechanics literature finding of a
resolvable proximal-to-distal "kinetic chain" lag (hip extends, then knee, then ankle
plantarflexes last, tens-of-ms apart) is specifically documented for BALLISTIC movements — vertical
jumps, throws — not a controlled squat's eccentric/concentric cycling. **A genuine test of "does
ankle peak before knee before hip" needs a plyometric jump clip** (candidates exist in the corpus —
`measured_athletics/` category videos with force-plate-verified countermovement jumps — not yet run
through this IK+force pipeline; concrete next step, §7), where real inter-joint lags are expected to
be large enough (order 50-100ms, literature) to clear this clip-class's own ~33ms resolution floor.

**A second, independent, non-obvious finding**: peak force does NOT occur at peak depth. The final,
deepest rep (132.3° knee flexion, at t≈12.07-12.23s) has its reaction-force peak at t=9.467s — 0.8s
into its analysis window, during the RAPID descent, not at the (slower, controlled) bottom position.
This rep's peak magnitudes (86-119%BW) are also the SMALLEST of the 5 reps, despite being the
deepest — consistent with peak reaction force tracking movement SPEED (acceleration), not position:
the fastest, least-controlled transition (rep 1, 2.9s) produces the biggest transient force, the
slowest/most-controlled deepest rep produces the smallest.

**Symmetric QC on this own build's flattering-looking "4/4 cuts agree" result** (§2.4): since `back`
carries no GRF term, its exact peak-time coincidence with the GRF-dominated leg cuts is CONSISTENT
WITH EITHER a genuine shared whole-body loading event (plausible — a rigid kinetic chain undergoing
one rapid, correlated transition) OR residual shared-anchor-noise propagation (the foot anchor's own
imperfectly-smoothed acceleration subtracted from every cut, including `back`). **This build does
not resolve which — disclosed as an open, machine-flagged ambiguity
(`back_cut_peak_coincidence_flag` in the JSON), not claimed as independent confirmation.**

## 5. Static Optimization: genuinely attempted twice, fails a hard quantitative gate both times

The task named Static Optimization alongside the Newton/virtual-work method. It was NOT skipped:

**Attempt 1 (naive, full 368-frame clip, `lowpass_cutoff_frequency_for_coordinates=6`):** killed
after continuous Ipopt `Restoration failed`/`Maximum iterations exceeded` + dozens of muscles
pinned at their upper activation bound, from before t=0 through at least t=0.8s (same signature
`docs/MECHANISM_STATIC_OPT.md` documented for FABRICATED/extrapolated kinematics). **Forced OODA
diagnosis** (not just "SO is hard, move on"): the Tier-1 50/50 L/R GRF split (§2) is only physically
valid during genuine double-support squatting — this clip has walking-like transitions between reps
where that split is wrong-by-construction, and SO's own dynamics-consistency requirement is far
stricter (every frame, not just the aggregate peak) than the whole-body GRF-magnitude plausibility
check that already passed. **Fix**: restrict to the final rep's own double-support window only
(t=[8.667, 12.233]s, ~108 frames).

**Attempt 2 (final-rep window only):** the optimization itself actually completed — 109 frames
written to `walking1_StaticOptimization_activation.sto`/`_force.sto`, 0 NaN, activations bounded
[0.01, 1.0] — `AnalyzeTool.run()` then crashed with a native `double free or corruption (out)` (a
SimTK/SWIG-level memory error, occurring AFTER the analysis output was already on disk; not traced
further, out of scope for this build). Because the files were complete, they were independently
re-analyzed with `static_opt_knee.py`'s own already-proven `check_so_convergence_and_sanity()` gate,
reused verbatim (not re-implemented) — this is where the real, decisive, quantitative negative
comes from, machine-computed, not a guess from the crash alone:

| gate (same precedented thresholds as `docs/MECHANISM_STATIC_OPT.md`) | threshold | measured | verdict |
|---|---|---:|---|
| Activation bounds / NaN | [0,1], 0 NaN | [0.01, 1.0], 0 NaN | PASS |
| Pelvis residual force peak | <75 N ("good") | **1153.1 N** (15.4x over) | **FAIL** |
| Pelvis residual moment peak | <75 N·m ("good") | **1092.0 N·m** (14.6x over) | **FAIL** |
| Joint reserve saturation (thin safety net) | <50% of own optimal_force | **hip_flexion_l: 9182%; hip_adduction_l: 6303%; knee_angle_l: 3173%** (and 9 more coordinates over 50%) | **FAIL** |

**Verdict: FAIL, decisively — not a borderline call.** The "reserve" actuators (meant, per this
repo's own SO convention, as a thin safety net covering pelvis-consistency residuals) are doing
30-90x their intended role: the muscles alone cannot come close to producing the moments this
Tier-1-estimated GRF implies for large parts of even the single cleanest rep. **This is a forced,
twice-attempted, quantitatively-confirmed honest negative** — the estimated GRF's aggregate
magnitude is plausible (§2), but its frame-by-frame internal consistency is not good enough for
muscle-driven optimization, a materially higher bar. The delivered force-transmission profile (§4)
is the Newton reaction-force method; SO on this corpus clip remains open (`docs/
MECHANISM_STATIC_OPT.md`'s own SO recipe DID work — on subject2's real force-plate walking trial;
the difference here is the GRF's provenance, not the SO recipe itself).

## 6. Honest limits (explicit, machine-quantified, not hedged in prose only)

1. **Absolute force/GRF magnitudes are approximate, capped by the monocular pipeline's own
   documented limits** (`docs/MECHANISM_CORPUS_IK.md` §5.2) — there is no force plate and no
   reliable absolute translation; every number here rests on (a) the stance-foot-anchored geometric
   reconstruction (§2, itself cross-validated at r=-0.92 and two independent mass-shedding checks)
   and (b) a Tier-1 50/50 L/R GRF split, valid for double-support instants only. Treat the
   PATTERN/SEQUENCE and RELATIVE magnitudes as the trustworthy content; treat absolute %BW numbers
   as plausible-order-of-magnitude, not a precision measurement.
2. **Ankle is the documented weakest-observed level** (`docs/MECHANISM_CORPUS_IK.md` §3.2/3.5/5.5:
   ankle range gate FAILs; subtalar oscillates near its full ±35° mechanical span) — this is
   directly why the anchor body (the foot) needed extra smoothing (§2.3) before use, and its own
   cut's numbers should be trusted least of the 4 levels reported here.
3. **The "back" cut's peak-time coincidence with the leg cuts is an unresolved ambiguity** (§4.2,
   §2.4) — real shared dynamics vs. residual shared-anchor noise, not distinguished by this build.
4. **No resolvable timing sequence for THIS movement class** (§4.2) — a squat's peak forces across
   all 4 levels are simultaneous within this clip's ~33ms resolution; this is a genuine, measured
   property of a controlled squat, not evidence against the classic ballistic-movement sequencing
   literature, which this clip cannot test.
5. **Cross-check A (right-foot-planted assumption) is inconclusive**, confounded by the
   already-documented right-leg occlusion in this specific clip (§3).
6. **Static Optimization does not converge cleanly on this clip** (§5) — the delivered profile is
   the bug-immune Newton method only.
7. **Single clip, single subject's anthropometry (subject2's scaled geometry stands in for this
   athlete)**, same inherited limit as every corpus-IK cert in this repo.
8. **Tier-1, straight-body-segment dynamics** — no wrap objects, no muscle-driven inertial terms
   beyond what SO would add (§5); the Newton method sees only rigid-segment kinematics + estimated
   external force, same category of limit `docs/MECHANISM_JOINT_FORCE_VALIDATION.md` already
   disclosed for the reaction-force-vs-contact-force distinction.

## 7. Next step

1. **Run a genuine plyometric jump clip** through this exact pipeline (candidates already in the
   corpus manifest under `measured_athletics/` — force-plate-verified countermovement-jump videos —
   not yet processed through `pose_to_opensim_ik.py`) to properly test the proximal-to-distal
   ballistic-sequencing hypothesis this squat's own ~33ms resolution could not resolve (§4.2).
2. **If a tighter SO result is later needed**: either obtain/estimate a per-frame-consistent GRF
   (not just aggregate-plausible), or diagnose the AnalyzeTool native crash (§5) upstream, before
   re-attempting a third time — per this build's own lean discipline, not attempted here.
3. **Resolve the back-cut coincidence ambiguity (§4.2)**: ground `back` using an anchor independent
   of the leg chain's own noise (e.g., a second, arm-based or trunk-marker-based stance check) to
   determine whether the simultaneous peak is genuine shared dynamics or residual contamination.

## Files

- `scripts/msk/force_transmission_scene.py` — the full pipeline (self-contained, re-runnable;
  imports `validate_joint_force.py`'s proven `parse_mot`/`savgol_smooth_and_derivs`/
  `get_descendant_bodies`/`G`, and reuses `static_opt_knee.py`'s `check_so_convergence_and_sanity`
  for the SO diagnostic — neither re-implemented). `FTS_RUN_SO=0` env var skips the slow/crash-prone
  Static Optimization step for a fast primary-method-only run.
- `data/msk_smoketest/force_transmission_scene/2025-09-26_DPEOqLVkav7/force_transmission_scene_results.json`
  — every number in §1-4, machine-written (primary method).
- `data/msk_smoketest/force_transmission_scene/2025-09-26_DPEOqLVkav7/static_opt/` — the two SO
  attempts' patched setup XMLs, synthetic GRF `.mot`, and the completed (pre-crash)
  `walking1_StaticOptimization_activation.sto`/`_force.sto` from attempt 2, plus
  `so_attempt_diagnostic.json` (the §5 gate table, machine-written from re-analyzing those files).
