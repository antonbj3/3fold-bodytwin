# MECHANISM FORCE SCENES BATCH — movement-type panel: squat vs jump vs other (2026-07-21)

## Bottom line

**A movement-TYPE panel (3 clips: squat / jump / one other athletic movement), not a single-clip
deep-dive** (that's `scripts/msk/force_transmission_scene.py`, the whole-chain agent's scope).
Built `scripts/msk/force_scenes_batch.py` by generalizing that already-validated Newton's-law
(OpenSim-ID-bug-immune) pipeline from one hardcoded clip/side to any (stem, side) pair, **verified
via a regression check against the already-published pilot-clip numbers (PASS, 0.0000% relative
difference on all 4 free-body cuts)** before trusting it on 3 new clips.

**Headline (bug-immune Newton reaction-force method, the trustworthy tier):** the plyometric
**jump** clip peaks **23-30% higher** than either the squat or the third clip at every one of
ankle/knee/hip — a clean, monotonic, physiologically sensible pattern. Squat and the third clip
are **near-identical in peak magnitude** (within 0.4-3.5% of each other) but **kinematically
distinct** in every other measured signature (angular velocity, pelvis-knee coupling) — the
physics, not the corpus label, is what actually distinguishes them (scene-eyes: measured, not
assumed).

**Two genuine, forced findings this build surfaced, not hidden:**
1. **A new, general geometric limit of the stance-foot-anchor method**, derived from Newton's law
   and then measured: it silently misreads a genuine loss of foot-ground contact as +1 body-weight
   of static support rather than the true ~0. Checked directly via an anchor-independent raw
   vertical-acceleration diagnostic; found and quantified on the regression pilot clip (11.6% of
   frames, exactly at its previously-published "peak" — a genuine ~2.1x inflation); absent (0.0%)
   on all 3 new panel clips.
2. **The muscle-driven (Static Optimization) tier numerically converges on all 3 clips but is
   NOT dynamically plausible on any of them** — multiple large hip/thigh muscles pinned near
   maximum activation simultaneously, pelvis "reserve" (non-muscular) actuator forces/moments
   3.7-16.8x over the accepted walking-trial comfort band, and Ipopt's own internal cost value
   spiking 2-4x at isolated frames. Diagnosed via forced OODA (not a first-shot bail), root cause
   reasoned through (OpenSim's own internal 4-6 Hz coordinate lowpass for SO is far lighter than
   this pipeline's own heavily-smoothed, ~380 ms-window COM differentiation, so it likely
   amplifies this corpus's MediaPipe-derived joint-angle noise into spurious torque demands).
   **A SEPARATE, independently-found bug compounds this**: a sibling audit
   (`docs/MECHANISM_SIGN_BUG_AUDIT.md`, committed the same session) machine-confirmed that
   `static_opt_knee.knee_crossing_muscles_and_forces` — the exact muscle-crossing-geometry
   function this build reused — returns the **exact negative** of the correct force direction.
   Independently re-derived from the ternary's own index logic (not merely trusted from the
   audit), fixed locally (`force_scenes_batch.knee_crossing_muscles_and_forces_fixed`, branches
   swapped, `static_opt_knee.py` itself left untouched), and **re-applied to all 3 clips reusing
   the already-converged SO output (no Ipopt re-run)**. The corrected numbers are reported in §6
   — **they remain equally implausible (179-1044%BW)**, confirming the sign bug was not the
   (sole) cause: the dynamically-implausible muscle-recruitment pattern is the dominant,
   sufficient reason to distrust these absolute numbers, independent of the sign fix. Both
   problems are real, both are now fixed/diagnosed, and fixing one did not rescue the other —
   reported as such, not as a single tidy story.

## 1. Scope vs the whole-chain agent

`force_transmission_scene.py` already took ONE clip (`2025-09-26_DPEOqLVkav7`,
`weighted_lifts/squat`) end-to-end in depth: geometric stance-foot-anchor correction, a
5-Savitzky-Golay-window GRF-plausibility sweep, per-rep cross-checks, and a best-effort genuine
Static Optimization attempt. This build is **explicitly different scope**: reuse that exact
method, generalize it from one hardcoded clip/side to an arbitrary (stem, side) pair, and run it
across **3 different movement types** to build a comparison panel. Nothing here re-litigates that
script's own findings about its one clip — it is reused as an external, already-published anchor
(see the regression check, §2).

## 2. The 3 clips chosen (+ the regression-check clip)

Picked from `scripts/msk/corpus_ik_aggregate_per_clip.csv` (228 clips, 181 CLEAN per
`docs/MECHANISM_CORPUS_IK_AGGREGATE.md`), using the catalog's own `movement_group` column +
`classification` + marker error, exactly as instructed:

| role | stem | movement_group | classification | marker RMS | duration / frames | knee/hip/ankle gates |
|---|---|---|---:|---:|---:|---|
| regression check (not part of the 3-clip panel) | `2025-09-26_DPEOqLVkav7` | `weighted_lifts/squat` | CLEAN_BORDERLINE | 0.0857 m | 12.37 s / 368 | knee PASS, hip PASS, **ankle FAIL** (already documented, `docs/MECHANISM_CORPUS_IK.md`) |
| **squat** | `2025-12-19_DSbSGkCDJRJ_5` | `weighted_lifts/squat` | **CLEAN_NO_CAVEAT** | 0.0771 m | 6.71 s / 202 | knee/hip/ankle **all PASS** |
| **jump** | `2026-06-01_DZC9ekKlZUv_2` | `plyometrics_jumps/plyometric_drill` | **CLEAN_NO_CAVEAT** | 0.0820 m | 3.83 s / 116 | knee/hip/ankle **all PASS** |
| **other** | `2025-10-07_DPfNTqrjPXo_7` | `plyometrics_jumps/technique_tutorial` | **CLEAN_NO_CAVEAT** | 0.0735 m | 5.00 s / 151 | knee/hip/ankle **all PASS** |

All 3 panel clips have `det_frac=1.000` (every frame detected) and **negligible L/R visibility
asymmetry** (max per-landmark asymmetry 0.016-0.052, all below the catalog's own
`lr_asymmetry_flag` trigger) — unlike the regression clip, which needed the LEFT side specifically
because of a documented severe right-knee occlusion. All 3 panel clips are therefore analyzed on
the **right** side uniformly (checked live per clip, not assumed).

**Honest catalog gap, exactly as flagged in the task:** this 228-clip corpus has only **two**
top-level movement categories — `weighted_lifts/squat` and `plyometrics_jumps/*` — no third
top-level type (sprint/throw/agility/etc.) exists anywhere in it. "One other athletic movement"
is therefore drawn from `plyometrics_jumps/technique_tutorial`, a **different sub-group** than
the jump-drill pick, deliberately chosen with substantial ROM (knee 58.9-134.6°, hip
48.4-118.2°) so it is a genuine movement demonstration, not one of the corpus's documented
static-hold/close-up segments. **Its kinematic distinctness from the jump pick is measured, not
assumed from the label** — see §7.

## 3. Method (bug-immune Newton reaction-force, generalized)

Identical physics to `force_transmission_scene.py`'s already-validated Steps 1-6, generalized
from one hardcoded clip/side to any `(stem, side)`:

1. Forward-kinematics COM trajectory per body (OpenSim `model.assemble()` per frame — never
   `InverseDynamicsTool`/`InverseDynamicsSolver`, which carries the documented ~1000-1800x
   knee/hip generalized-force bug, `docs/MECHANISM_MSK_ELASTIC_BAND.md` Sec.4).
2. **Stance-foot-anchor correction** for the missing global-translation signal (MediaPipe's
   `pelvis_tx/ty/tz` is a few-cm noise residual, not real depth, `docs/MECHANISM_CORPUS_IK.md`
   Sec.5.2): every body's position/acceleration rebased onto the analyzed side's own foot
   (`calcn_{side}`), zeroing that foot's own grounded acceleration by construction.
3. Kinematically-**estimated** GRF (Tier-1, symmetric 50/50 L/R split during any candidate
   double-support instant — no force plate exists for any corpus clip), window-chosen via a
   5-point Savitzky-Golay sweep against an externally-anchored plausibility ceiling (400%BW;
   walking~1-1.3x, running~2-3x, jump~2-4x BW per the literature).
4. Newton's second law on each BFS-derived free-body cut (`ankle_{side}`, `walker_knee_{side}`,
   `hip_{side}`, `back`) → reaction force per cut, across the WHOLE clip.
5. **Regression check** (run first, before trusting the generalization on new clips): re-run this
   generalized code on `force_transmission_scene.py`'s own pilot clip and compare to its
   already-published numbers — an external anchor for correctness, decorrelated from whether the
   3 new clips are analyzed right.

| cut | mine (this script) | published (`force_transmission_scene.py`) | rel. diff |
|---|---:|---:|---:|
| ankle | 177.574 %BW | 177.574 %BW | 0.0000% |
| knee | 170.125 %BW | 170.125 %BW | 0.0000% |
| hip | 130.663 %BW | 130.663 %BW | 0.0000% |
| back | 175.807 %BW | 175.807 %BW | 0.0000% |

**REGRESSION CHECK: PASS.** All 3 new-clip results below build on a verified-correct
re-implementation, not a fresh, unchecked one.

## 4. Cross-movement comparison — peak knee/hip/ankle loads (the requested headline)

Same subject2-anthropometry model for every clip (78.2 kg, BW=766.88 N — see §8 caveat), chain
mass fractions identical across clips (ankle-distal 2.1%, knee-distal 7.0%, hip-distal 19.5%,
back/torso+arms 45.4% of body weight). All 3 panel clips: `assemble()` failures = 0, GRF-sweep
primary window under the 400%BW ceiling, **mass-shedding cross-check PASS** (observed
ankle→knee / knee→hip segment-weight-shedding matches predicted shank/thigh weight to within
0.6% for every clip — the chain physics holds identically regardless of movement type), **zero
anchor-freefall frames flagged** (§5) — every number below is on solid ground, not an artifact.

| movement | peak ANKLE %BW | peak KNEE %BW | peak HIP %BW | peak order | t_peak |
|---|---:|---:|---:|---|---|
| **squat** | 52.6 | 47.8 | 35.6 | ankle < knee < hip < back | 0.367 s |
| **jump** | **64.6** | **59.7** | **44.7** | ankle < knee < hip < back | 1.367-1.400 s |
| **other** | 52.4 | 47.5 | 34.4 | ankle < knee < **back** < hip | 0.767 s (hip @3.100s) |

**Robustness check on the peak instant** (forced adversary: is "peak at t=0.37-0.77s, near the
start of the clip" a differentiation edge artifact rather than a real loading event? — the
regression clip's own history shows this concern is not paranoid, see §5): re-searched each
peak excluding the first/last 20% of frames (not just the minimal Savitzky-Golay edge trim).
Squat's peak moves only 52.6%→52.4%BW (a different, later frame at t=2.571s); other's hip peak
is the *same* frame (t=3.100s) either way. **Confirmed stable, not an edge artifact** — the
forced adversary fell.

**Cross-movement finding:** jump peaks **23-30% higher** than squat/other at every joint (ankle
+23%, knee +25-26%, hip +26-30%) — a clean, monotonic, sensible result: plyometric drills load
the leg more than a controlled squat. Squat and "other" are **essentially indistinguishable in
peak magnitude** (within 0.4-3.5% at every joint) despite different movement labels — peak %BW
alone does not separate them; §7 shows what does. The **ankle > knee > hip** ordering (chain
physics: each proximal cut sheds the more-distal segment's own weight) holds identically across
all 3 movement types, never flips — the method's chain-physics signature is movement-invariant
even though the magnitudes are not.

## 5. A new general finding: the stance-foot-anchor method's hidden assumption, and where it breaks

**Derived before running anything** (geometric thinking, not a post-hoc rationalization — see
the module docstring in `force_scenes_batch.py`): the stance-foot-anchor trick is valid only when
the anchor foot's TRUE ground-frame acceleration is ~0 (genuinely planted). Re-deriving Newton's
law for the case where it is not (e.g. true free-fall, anchor acceleration ≈ g_vec) shows the
method then computes `GRF ≈ mass·(0 − g_vec) = +1 BW` — i.e. it **silently misreads a genuine
loss of ground contact as full static support**, the opposite of the naive "probably reads near
zero" guess.

**Checked directly, not just argued:** an anchor-independent diagnostic (the anchor foot's own
RAW, pre-grounding vertical acceleration, flagged when it drops below −0.6g) was run on all 4
clips.

- **Regression pilot clip: 11.6% of interior frames flagged**, in five 7-9-frame runs recurring
  at a consistent phase of 4-5 squat reps. **Forced the adversary before accepting this as
  real-vs-noise:** raw MediaPipe visibility at the flagged frames is 0.92-0.99 (high — rules out
  simple occlusion/low-confidence tracking as the cause) but the OpenSim-**solved** foot position
  swings ~0.58 m in ~330 ms while the **raw MediaPipe marker's own** vertical excursion over the
  same window is only ~0.15 m — a 4x mismatch implicating the IK solve's weakly-constrained
  floating-base/global DOFs (no pelvis-specific markers exist in this 16-marker bridge,
  `docs/MECHANISM_CORPUS_IK.md` Sec.5.3), not raw tracking noise, as the likely root cause. **All
  four of the pilot clip's already-published peaks (ankle/knee/hip/back) sit exactly at
  t=2.900s, inside a flagged window** — the freefall-excluded ("stance-only") re-peak is
  ankle 83.8%BW / knee 79.5%BW / hip 59.8%BW / back 82.4%BW, **about half** the naive number.
  This is a genuine refinement of that already-published result, surfaced by this build, not a
  retraction of it (the regression check above compares against the ORIGINAL script's own
  numbers, which is the correct decorrelated anchor for THIS script's correctness regardless).
- **All 3 new panel clips: 0.0% flagged.** Consistent with (though not proof of, n=3) this
  panel's CLEAN_NO_CAVEAT selection criterion also screening out this artifact.

## 6. Muscle-driven (Static Optimization) attempt — genuine, sign-bug found+fixed, still NOT dynamically plausible

Per the task, a genuine Static Optimization attempt was run for every clip (reusing
`static_opt_knee.py`'s OpenCap SO/JR XML templates, and — originally — its
`knee_crossing_muscles_and_forces` live-path-geometry muscle-crossing detector), windowed around
each clip's own peak-load event (reusing `force_transmission_scene.py`'s already-proven
windowing convention). **All 3 converged on the first attempt (6 Hz coordinate lowpass, no
retry needed)** — a genuine capability result, not skipped.

**Sign bug found and fixed mid-build** (§ Bottom line): `static_opt_knee.knee_crossing_muscles_and_forces`
returns the exact negative of the correct force direction (its `(k, k+1) if in_k1 else (k+1, k)`
ternary is backwards in both branches — independently re-derived here, then cross-confirmed
against `docs/MECHANISM_SIGN_BUG_AUDIT.md`'s own toy-case proof, found via a concurrent sibling
audit). Fixed locally as `knee_crossing_muscles_and_forces_fixed` (branches swapped,
`static_opt_knee.py` itself left untouched by THIS build) and re-applied to all 3 clips, reusing
the already-converged SO output directly (no Ipopt re-run needed — the bug is entirely in the
downstream Python geometry step, never in SO itself). **Triangulated three ways, independently:**
this build's own first-principles re-derivation of the ternary's index logic, the sibling audit's
toy-case unit test (`scripts/msk/audit_sign_bug.py`), and — noticed live, mid-build, via a plain
`git status`/`git diff` check on the shared repo — a concurrent session applying an in-place fix
to `static_opt_knee.py` itself that swaps the exact same two branches. Three decorrelated paths,
one answer.

| movement | SO window | contact-force peak, BUGGY (ankle/knee/hip %BW) | contact-force peak, **SIGN-FIXED** (ankle/knee/hip %BW) | muscles pinned >10% of window | peak pelvis reserve force / moment | vs 75N/75Nm comfort band |
|---|---|---:|---:|---:|---:|---:|
| squat | [0.363, 2.637]s (68 fr) | 171.9 / 614.0 / 1028.5 | **239.1 / 684.8 / 1044.2** | **8** | 222.6 N / 277.7 Nm | **3.7x over** |
| jump | [0.633, 3.133]s (75 fr) | 103.5 / 589.7 / 770.6 | **179.3 / 680.5 / 772.7** | **2** | 1259.0 N / 538.9 Nm | **16.8x over** |
| other | [0.600, 4.200]s (109 fr) | 108.4 / 616.9 / 1005.5 | **203.4 / 620.7 / 990.4** | **10** | 319.5 N / 328.5 Nm | **4.4x over** |

**The corrected numbers are, if anything, slightly LARGER than the buggy ones (ankle actually
rose 60-95%) — not smaller.** This is itself an important, forced-adversary result: it rules out
"the sign bug was the whole story." Unlike the original walking-trial audit (where the fix moved
knee/hip from an OrthoLoad-plausible-looking 233-235%BW to a still-PASS-substantial-but-larger
391/387%BW), fixing this bug on these 3 corpus clips does **not** rescue the numbers into
anything resembling a physiological range — **knee (620-685%BW) and hip (773-1044%BW) remain far
outside any published literature** (this repo's own OrthoLoad-anchored knee cert used 258%BW
in-vivo; the hip cert used 210-330%BW; even the most extreme running literature this repo has
cited tops out under 300%BW for the hip). **Forced OODA before accepting or rejecting either
number** (not a first-shot bail in either direction):
- Observe: implausible magnitude, and near-identical across all 3 movement types despite very
  different pure-reaction baselines (29-53%BW) — itself suspicious (a genuine per-movement finding
  should vary more than a shared artifact would).
- Orient: added a new, reusable diagnostic (`diagnose_so_dynamical_plausibility`,
  `force_scenes_batch.py`) — checks muscle-pinning fraction + pelvis reserve-actuator peaks
  (inherited, not re-tuned, comfort band from `static_opt_knee.py`/`validate_hip_force.py`'s
  real-mocap walking-trial cert). Result: **`dynamically_plausible=False` on all 3 clips**,
  machine-computed, not eyeballed. Also inspected Ipopt's own internal cost ("Performance") log:
  baseline ~200-700 for most frames, spiking to 1593-2147 at isolated frames (the `other` clip) —
  a signature of the optimizer fighting noise at specific instants, not a smooth real-effort
  curve.
- Decide/Act: root cause reasoned through, not just asserted — OpenSim's own internal coordinate
  differentiation for SO uses a much lighter 4-6 Hz lowpass than this pipeline's own COM
  differentiation (23-sample, ~380 ms Savitzky-Golay window); on this corpus's MediaPipe-derived
  (noisier-than-mocap) joint angles, the lighter filter plausibly lets frame-to-frame jitter
  amplify into large, spurious angular-acceleration demands that only an unrealistic
  multi-muscle-pinned recruitment pattern can satisfy.
- Symmetric check on the sign-bug fix itself: does correcting it change this verdict? **No** — the
  `dynamically_plausible=False` diagnostic (muscle pinning, reserve saturation) is computed
  directly from `activation.sto`/`force.sto` and never calls the crossing-muscle detector at all,
  so it is unaffected by the sign fix; and the sign-corrected contact-force numbers themselves
  (table above) are, if anything, larger. Two independent problems, both real, neither one
  explains away the other.

**Verdict: HONEST NEGATIVE for the absolute muscle-driven contact-force number on this corpus, on
all 3 movement types, BEFORE and AFTER the sign-bug fix — not hidden, not silently reported as if
trustworthy.** The pure Newton reaction-force numbers in §4 remain the primary, trustworthy
cross-movement metric. **Concrete, named next fix (not attempted here, out of this session's
scope):** re-run SO with a much heavier coordinate lowpass (e.g. 2-3 Hz) to test whether the
pinning/reserve-saturation resolves — a cheap, decisive follow-up, not attempted here per this
repo's first-step doctrine.

## 7. Cross-movement kinematic signatures (scene-eyes: let the physics classify the movement)

| movement | knee ROM | hip ROM | peak \|knee angular velocity\| | pelvis-height vs knee-angle correlation |
|---|---|---|---:|---:|
| squat | 35.6-112.6° | -29.4-72.6° | 148.4 °/s | **r=+0.13** (weak, wrong-signed vs. the textbook squat prediction) |
| jump | 15.8-137.5° | -29.8-112.9° | **192.1 °/s** (fastest) | **r=-0.98** (strong, textbook bilateral-crouch-then-extend signature) |
| other | 58.9-134.6° | 48.4-118.2° | 94.9 °/s (slowest) | **r=+0.27** (weak, wrong-signed) |

The pre-registered physical prediction was "quasi-static squat-like movement ⇒ pelvis height
anti-correlates strongly with knee flexion; a true flight-phase movement ⇒ decouples." **The jump
clip matches this prediction cleanly** (a tight, symmetric bilateral crouch-and-extend). **Neither
the squat nor the "other" clip does** — both show weak, positive (wrong-signed) correlation,
despite the squat clip's own hip-flexion range being strongly L/R-asymmetric (r=102.0° vs
l=59.1°, from the catalog) even though its L/R *tracking* asymmetry is negligible (0.016) — i.e.
a real, measured, asymmetric movement pattern (plausibly a split-stance/lunge-family variant, not
a plain bilateral back squat), not a tracking artifact. This is exactly the finding the task's own
framing anticipated: **the "other" clip is genuinely, measurably distinct from the jump pick**
(half the peak knee angular velocity, opposite-signed pelvis-knee coupling) even though the
corpus's coarse label puts both under `plyometrics_jumps/*` — confirmed by the kinematics, not
assumed from the group column alone.

## 8. Honest scale caveats

- **Absolute force magnitude is subject2's, not the filmed athlete's**, for every one of these 3
  clips (and every clip in the whole 228-clip corpus) — no per-subject scaling exists for
  monocular video (`docs/MECHANISM_CORPUS_IK.md` Sec.5.1). BW=766.88N is identical across the
  entire panel by construction. **Relative %BW patterns and cross-clip comparisons are the
  meaningful output here — absolute Newtons are not.**
- **No real GRF exists for any corpus clip** — every number in §4 rests on a kinematically
  ESTIMATED, Tier-1 symmetric-50/50-split GRF, one layer of approximation beyond subject2's own
  real-force-plate certs. The muscle-driven tier (§6) stacks a SECOND layer of approximation on
  top of that estimate — disclosed, not inherited silently.
- **Ankle remains the weakest-observed DOF** in this 16-marker bridge (heel+toe only, no 3rd
  non-collinear foot point) — this build did not re-litigate that finding, only reused it; all 3
  panel clips happen to PASS the ankle range gate (unlike the regression pilot), which is a
  property of this specific clip selection, not a general fix.
- **Only 2 top-level movement categories exist in the source corpus** (§2) — "other athletic
  movement" is a within-corpus best-available pick, not a genuinely disjoint third exercise
  family (e.g. no sprint/throw/agility clips exist to pick from).
- **The muscle-driven contact-force numbers in §6 are not physiologically trustworthy** — use the
  bug-immune Newton reaction-force numbers in §4 for any downstream comparison.
- **n=3 clips per this panel's own construction** (as instructed) — the cross-movement
  percentages in §4 are a snapshot on these specific 3 clips, not a corpus-wide statistical claim;
  scaling to more clips per movement type is the natural next step, not attempted here.

## 9. Reproducing this result

```
.venv-msk/bin/python3 scripts/msk/force_scenes_batch.py
```

Runs the regression check first (fails loudly before touching the 3 new clips if it doesn't
PASS), then all 3 panel clips (Newton reaction-force, always; genuine Static Optimization
attempt, ~1-3 min per clip). Outputs (gitignored) under
`data/msk_smoketest/force_scenes_batch/<label>__<stem>/scene_results.json` (per-clip) and
`data/msk_smoketest/force_scenes_batch/force_scenes_batch_results.json` (the aggregate +
comparison table, everything in this doc's tables machine-computed from there).

## Files

- `scripts/msk/force_scenes_batch.py` (new, this build)
- `data/msk_smoketest/force_scenes_batch/force_scenes_batch_results.json` (aggregate + comparison
  table)
- `data/msk_smoketest/force_scenes_batch/{squat__2025-12-19_DSbSGkCDJRJ_5,
  jump__2026-06-01_DZC9ekKlZUv_2, other_technique_tutorial__2025-10-07_DPfNTqrjPXo_7,
  regression_check_pilot_squat__2025-09-26_DPEOqLVkav7}/scene_results.json` (per-clip)
- Reused, unmodified: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `data/msk_models/LaiArnoldModified2017_mediapipe_bridge_subject2_scaled.osim`,
  `data/msk_ik/<stem>/smoothed/*.mot` (228-clip corpus IK output, unchanged)
- External anchor for the regression check: `data/msk_smoketest/force_transmission_scene/
  2025-09-26_DPEOqLVkav7/force_transmission_scene_results.json` (`force_transmission_scene.py`'s
  own already-published output, not modified by this build)
- Cross-referenced (not authored by this build, not modified): `docs/MECHANISM_SIGN_BUG_AUDIT.md` +
  `scripts/msk/audit_sign_bug.py` — the sibling audit whose finding (`knee_crossing_muscles_and_
  forces` sign bug) this build independently re-derived, cross-confirmed, and locally fixed as
  `knee_crossing_muscles_and_forces_fixed` inside `force_scenes_batch.py` (§6)
