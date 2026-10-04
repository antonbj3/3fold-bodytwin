# MECHANISM PLYOMETRIC JUMP CASCADE — does a ballistic jump resolve ankle→knee→hip→lumbar force timing? (2026-07-21)

Executes the whole-chain agent's flagged next step (`docs/MECHANISM_FORCE_TRANSMISSION_SCENE.md`
§7): that build's controlled squat found all 4 chain levels peaking in the SAME 30fps frame, every
rep — a tie, not a cascade, attributed to "a squat isn't ballistic." Script:
`scripts/msk/plyo_jump_cascade.py` (`.venv-msk/bin/python3 -u scripts/msk/plyo_jump_cascade.py`).
Every number below is this script's own printed/JSON output — never estimated, never eyeballed
from a figure.

## Bottom line

**A genuine plyometric jump, with a real, physics-cross-checked, multi-frame flight phase, STILL
does not resolve a proximal-to-distal force cascade at this corpus's ~30fps (33ms) resolution.**
Ankle, knee, hip, and back reaction-force peaks occur at the IDENTICAL frame during the propulsion
(on-ground, pre-takeoff) phase in **both** of the chosen clip's independent jump instances — stable
across a 3-value pre-margin sensitivity sweep and 4 of 5 Savitzky-Golay smoothing windows (the one
exception did not survive its own cross-event/cross-window re-check, see §5.3). This **extends**,
rather than overturns, the squat build's finding: even a ballistic movement's own concentric
propulsion sub-phase can be too fast for this pipeline's temporal resolution to stagger.

**Ground-contact-loss handled explicitly, as required**: this clip's own airborne phase (10.3% of
interior frames) triggers the documented free-fall misread (`force_scenes_batch.py`'s own
derivation: true GRF≈0 read as ≈+1 bodyweight) — the naive global peak sits exactly at a flagged
frame and is **1.29–1.37× inflated** relative to the stance-only (freefall-excluded) peak for this
clip (smaller than the ~2.1× found on `force_scenes_batch.py`'s own regression pilot clip — a
real, clip-specific number, not forced to match).

**Clip selection was machine-screened, not eyeballed, and surfaced an honest, disclosed corpus
finding**: of the corpus's CLEAN plyometric-jump clips with a resolvable flight phase, the large
majority show a knee-angle **rising** (tucking) through takeoff, not the classic crouch-then-EXTEND
countermovement-jump shape — this Instagram-sourced corpus's "jump" content skews toward
tuck/knee-drive drills. Disclosed in full in §1, not hidden behind the task's "ideally a
countermovement or drop jump" framing.

## 1. Clip chosen — machine-screened across the whole corpus, forced OODA mid-screen

Candidates: every `plyometrics_jumps/*` clip in the 228-clip corpus classified `CLEAN_NO_CAVEAT`
in `scripts/msk/corpus_ik_aggregate_per_clip.csv` (marker PASS, knee/hip/ankle range gates ALL
pass, no L/R visibility asymmetry, not pinned at a mechanical limit) — **37 clips**, all
independently re-verified live (`lr_asymmetry_flag=False` for every one) to justify using the
**right** side uniformly, same convention `force_scenes_batch.py` established.

**Screen (Step 1 of the script, re-run live, not quoted from memory)**: each candidate's anchor
foot (`calcn_r`) raw vertical acceleration is checked against the SAME `FREEFALL_FRAC_OF_G=0.6`
threshold `force_scenes_batch.py` already uses for its ground-contact-loss diagnostic, at a
15-sample (noisier, "first look") screening window, requiring a **contiguous run ≥3 frames
(~100ms)** to count as a real flight candidate (not a 1-frame blip):

| | count |
|---|---:|
| Candidates screened | 37 |
| Show ≥1 resolvable (≥3-frame) flight episode | 21 |
| …of those, show ≥1 genuine on-ground knee-**extension** frame immediately into takeoff (classic CMJ) | 6 |
| …of those, show **only** a knee-**flexion** ("tuck") signature through takeoff | 15 |

**Forced OODA, not skipped past**: the single longest-flight-run candidate
(`2025-08-05_DM8_eSNsE7p_3`, 363ms) was inspected frame-by-frame (raw knee angle / foot height /
foot acceleration, printed and read as numbers, never a figure) **before** being accepted, because
its own printed timeline showed knee angle **rising** continuously from 44.9° to 126.6° across the
exact frames leading into and past its own takeoff — a tuck-jump signature, not a crouch-then-
extend jump. A second classifier (walking backward from each candidate's own detected takeoff,
checking the SIGN of the knee-angle velocity in the immediately-preceding on-ground frames) was
built to check this systematically across all 37 clips, not just the one — table above is its
output. (A first draft of that classifier had the extension/tuck sign backwards — caught by
re-deriving the walk-back condition from the physical definition rather than trusting the first
result; both directions are now reported for every clip, §1's table and the script's own printed
per-clip detail.)

**Ranked by longest contiguous flight run** (screening window):

| stem | longest flight run (ms) | n flight events | freefall frac |
|---|---:|---:|---:|
| `2025-08-05_DM8_eSNsE7p_3` | 363.0 | 2 | 29.7% |
| `2025-07-24_DMeFzibMRh__2` | 306.0 | 1 | 28.1% |
| **`2026-06-30_DaMKH8aDBxu_6`** (chosen) | **264.0** | **2** | **7.9%** |
| `2025-07-24_DMeFzibMRh__4` | 234.5 | 3 | 33.3% |
| `2025-08-05_DM8_eSNsE7p_5` | 201.0 | 1 | 15.1% |

**Why the chosen clip, not the top-ranked ones**: the top 2 by raw flight-run-length were run
through the full physics battery during selection (same pipeline, `force_scenes_batch.run_scene`)
and both are — like nearly the whole corpus — tuck-style. The **chosen** clip,
`2026-06-30_DaMKH8aDBxu_6` (`plyometrics_jumps/plyometrics_jumps_misc`), was picked instead because
it is the strongest **all-around** candidate on the criteria the task actually specifies: low
marker error (0.074m RMS), a GRF-plausibility sweep that clears the ceiling at **every** one of 5
smoothing windows with the **tightest** spread of any candidate checked (199–210%BW — see §2), an
almost-exact mass-shedding cross-check (ratios 1.001/1.014, essentially matching predicted
segment weight to the percent), the strongest pelvis-height/knee-angle correlation (r=-0.989), and
**2 independent jump instances** in one clip for a repeatability check (mirroring the squat
build's own 5-rep-aggregation convention) — all while still clearing the ≥3-frame flight-resolution
bar (264ms, rank 3 of 21). It is TUCK-style like most of the corpus (disclosed, not hidden) — the
propulsion-window definition (§5) is deliberately built generic (on-ground frames before takeoff)
so it does not presuppose either kinematic shape.

## 2. Method — 100% reused Newton reaction-force pipeline, one new addition

Everything physics-side is `scripts/msk/force_scenes_batch.py` (`fsb`)'s already-validated,
bug-immune Newton's-second-law free-body-cut method (itself generalized from
`scripts/msk/force_transmission_scene.py`) — **never** `InverseDynamicsTool`/`InverseDynamicsSolver`
(the documented ~1000–1800× knee/hip generalized-force bug, `docs/MECHANISM_MSK_ELASTIC_BAND.md`
§4). **Static Optimization (muscle-driven tier) is deliberately NOT run** — `force_scenes_batch.py`
already found it dynamically implausible on this exact corpus (muscles pinned near-max
simultaneously, pelvis reserve actuators 3.7–16.8× over the walking-trial comfort band, on all 3 of
its own clips) — per the task's explicit instruction, only the trustworthy Newton-reaction tier is
used here.

**The one genuinely new piece** (neither sibling script built it): an explicit **propulsion-phase
window** — the on-ground (non-freefall-flagged) frames immediately preceding each detected takeoff
— in which each of the 4 free-body cuts' own local peak reaction-force TIME is measured and
ordered (§5). The sibling scripts' own knee-angle-threshold event windows are reproduced here too
(§5.1, "coarse" cross-reference) for comparability, but they were built for a squat's crouch-bottom
event and do not isolate the specific on-ground-to-takeoff sub-window the task's "during the
propulsion phase" wording asks for.

**Regression check** (same convention `force_scenes_batch.py` itself used against
`force_transmission_scene.py`'s published numbers): this script's own from-scratch recomputation of
the 4 global %BW peaks is checked against `force_scenes_batch.run_scene()`'s own independently
produced output for this exact clip, from an earlier, separate interactive run made during clip
selection — **before** this script's own pipeline existed.

| cut | this script | independent anchor | rel. diff | verdict |
|---|---:|---:|---:|---|
| ankle | 103.3132 %BW | 103.3132 %BW | 0.0000% | **PASS** |
| knee | 99.2832 %BW | 99.2832 %BW | 0.0000% | **PASS** |
| hip | 79.0357 %BW | 79.0357 %BW | 0.0000% | **PASS** |
| back | 106.4937 %BW | 106.4937 %BW | 0.0000% | **PASS** |

**REGRESSION CHECK: PASS** (exact match — this is a checked reuse of the method, not a silent
fork).

## 3. Base-battery cross-checks (pre-registered thresholds, machine PASS/FAIL)

| cross-check | threshold | measured | verdict |
|---|---|---:|---|
| `model.assemble()` clean (178 frames) | 0 failures | 0/178 | **PASS** |
| GRF-plausibility sweep, primary window (23) | <400%BW | 210.3%BW (full sweep: 204.8/209.3/210.3/209.6/199.3 across windows 15–31) | **PASS** (all 5 windows) |
| Pelvis-height vs `knee_angle_r` correlation | \|r\|≥0.70 | **r = -0.9889** | **PASS** |
| Mass-shedding, ankle→knee | ratio∈[0.5,1.5] vs shank's own mass | **1.001** | **PASS** |
| Mass-shedding, knee→hip | ratio∈[0.5,1.5] vs thigh's own mass | **1.014** | **PASS** |

**Overall method-trustworthiness gate** (assemble-clean + GRF-plausibility + mass-shedding ×2 +
regression-check): **PASS** — this gate is on whether the MEASUREMENT is trustworthy, not on
whether a cascade was "found" (honest-negative-is-PASS, this repo's own convention).

## 4. Ground-contact-loss (freefall) artifact — checked and handled explicitly, as required

`force_scenes_batch.py`'s module docstring derives (Newton's law, not assumed) that the
stance-foot-anchor correction silently misreads genuine free-fall as ≈+1 bodyweight of static
support. This clip's own airborne phase (10.3% of interior frames flagged) triggers it: the naive
global peak for **every** cut sits exactly at the flagged frame (t=4.467s):

| cut | naive (contaminated) | stance-only (freefall-excluded) | inflation |
|---|---:|---:|---:|
| ankle | 103.3 %BW | 80.1 %BW | **1.29×** |
| knee | 99.3 %BW | 75.5 %BW | **1.31×** |
| hip | 79.0 %BW | 57.5 %BW | **1.37×** |
| back | 106.5 %BW | 78.6 %BW | **1.36×** |

Smaller than the ~2.1× `force_scenes_batch.py` found on its own regression pilot clip — a genuine,
clip-specific measurement, not forced to reproduce that figure. **Every downstream cascade-timing
number in §5 explicitly excludes freefall-flagged frames from its search window** — this is not
just a magnitude caveat, it is applied structurally to the timing claim itself.

## 5. Propulsion-phase cascade timing — the task's core question

Two resolvable flight episodes detected (freefall run ≥3 frames): takeoffs at **t=1.400s** and
**t=4.367s**, each with 264ms of flight — two independent jump instances in one clip (the
diverse-instance-space check, mirroring the squat build's own 5-rep convention, on a smaller n).

### 5.1 Coarse cross-reference (whole knee-angle-threshold event window, established convention)

| event (knee-bottom t) | order (earliest→latest) | gaps (ms) | simultaneous? |
|---|---|---|---|
| t=1.53s | ankle < knee < hip < back | [0.0, 0.0, 0.0] | **ALL 4 TIED** |
| t=4.50s | ankle < knee < hip < back | [0.0, 0.0, 0.0] | **ALL 4 TIED** |

Same tie the squat found, at the coarser window — consistent with, not yet distinguishing from,
§5.2's finer test.

### 5.2 Primary result — propulsion-specific window (on-ground frames immediately before takeoff)

Pre-margin = 0.5s (this repo's own precedented range), on-ground frames only:

| takeoff | window | ankle peak | knee peak | hip peak | back peak | order | gaps (ms) |
|---|---|---:|---:|---:|---:|---|---|
| t=1.400s | 15 frames (495ms) | 75.7%BW @t=1.367s | 70.5%BW @t=1.367s | 52.0%BW @t=1.367s | 71.8%BW @t=1.367s | ankle<knee<hip<back | **[0.0, 0.0, 0.0]** |
| t=4.367s | 15 frames (495ms) | 74.7%BW @t=4.333s | 69.7%BW @t=4.333s | 51.6%BW @t=4.333s | 71.6%BW @t=4.333s | ankle<knee<hip<back | **[0.0, 0.0, 0.0]** |

**All four cuts peak at the exact same frame — one frame (33ms) before takeoff — in BOTH
independent jump instances.** This holds restricted to just the true serial ankle→knee→hip
subchain too (all three share the identical peak index, `back`'s parallel-branch ambiguity,
§5.4, is not driving this). Smallest inter-cut gap = 0.0ms vs the 33.0ms one-frame floor →
**AT/BELOW resolution floor — NOT a resolvable cascade at this frame rate**, in either instance.

**Pre-margin sensitivity** (0.3s / 0.5s / 0.6s): identical `ankle<knee<hip<back`-tied order at all
3 margins, both events — not an artifact of the specific 0.5s choice.

### 5.3 Forced adversary: Savitzky-Golay smoothing-window stability sweep

The literature's own reported inter-joint lags (tens of ms) are comparable to this pipeline's own
smoothing-window scale — a real adversary, not paranoia (per the squat build's own precedent of
checking this). Swept all 5 established windows:

| window | event 1 (t=1.40s) order | event 2 (t=4.37s) order |
|---:|---|---|
| 15 | ankle<knee<hip<back (tied) | ankle<knee<hip<back (tied) |
| 19 | **hip<back<ankle<knee (NOT tied)** | ankle<knee<hip<back (tied) |
| **23 (primary)** | ankle<knee<hip<back (tied) | ankle<knee<hip<back (tied) |
| 27 | ankle<knee<hip<back (tied) | ankle<knee<hip<back (tied) |
| 31 | ankle<knee<hip<back (tied) | ankle<knee<hip<back (tied) |

**4 of 5 windows tie for both events; only window=19/event-1 shows a non-tied order.** Forced
re-check (symmetric QC on this one "positive-looking" result, exactly as the discipline requires
before accepting ANY interesting-looking number): it does not reproduce at the immediately
adjacent windows (15 or 23), and event 2 stays tied even at window=19 itself — a single
non-robust outlier at one window for one of two events is the signature of a smoothing-window
sensitivity fluke, not a genuine signal. **Verdict: does not survive the sweep — the tie at the
established primary window (23) is the trustworthy reading.**

### 5.4 Symmetric QC carried forward from the squat build

`back` carries no GRF term (parallel sibling branch off the pelvis, per
`docs/MECHANISM_BAND_POSTERIOR_CHAIN.md` §1, re-affirmed here) — its exact peak-time coincidence
with the GRF-dominated leg cuts is, as in the squat build, consistent with either genuine shared
whole-body loading or shared-anchor-noise propagation, not independently resolved by this build
either. **What this build adds**: restricting the tie check to just `ankle`/`knee`/`hip` (the true
nested serial chain, no ambiguity about GRF terms) still shows all 3 peaking at the identical
frame — so the core finding does not depend on `back`'s own ambiguous status.

## 6. Why a genuine flight phase still didn't produce a resolvable cascade (geometric reasoning, not just a shrug)

All 4 cuts' reaction-force series share the SAME dominant driver right before takeoff: the
whole-body acceleration signal entering `GRF_side` (common to the ankle/knee/hip cuts) and the
shared anchor-correction term (common to all 4, including `back`). A genuine inter-segment TIMING
stagger would need each cut's own distinct, chain-specific acceleration term to diverge in time by
more than one frame period from this shared driver — for this clip's specific propulsion event,
it does not. This is the same geometric point the squat build made (§2.4 there): a rigid,
tightly-coupled kinetic chain undergoing one rapid, correlated transition does not, by itself,
guarantee a resolvable lag — the lag has to be large enough relative to the sampling period, and
here (as in the squat) it is not.

## 7. Honest limits

1. **Absolute force magnitudes remain approximate** (subject2's anthropometry stands in for this
   athlete; GRF is a kinematically ESTIMATED, Tier-1 symmetric 50/50 L/R split, valid for
   double-support instants only) — relative patterns and TIMING are the trustworthy content,
   exactly as the task asked to treat them.
2. **This corpus's plyometric-jump content is predominantly tuck/knee-drive-style** (§1) — a
   classic lab countermovement/drop jump with a clean on-ground crouch-then-extend signature is
   NOT well represented (only 1/11 resolvable-flight clips at the established primary window
   showed it in an earlier manual check during clip selection; 6/21 at the noisier 15-sample
   screening window used by the script's own Step 1 — the two counts differ because a smaller,
   noisier smoothing window flags more marginal candidate flight events, itself an honest,
   disclosed methodological detail, not a discrepancy to paper over). The chosen clip is
   tuck-style; the propulsion-window definition does not presuppose either shape, but a true
   lab-style CMJ was not available in this corpus to test separately.
3. **n=2 independent jump instances** in the chosen clip (fewer than the squat build's 5 reps) —
   both agree, but this is a thinner repeatability base than the squat's own.
4. **The window=19/event-1 non-tied result (§5.3) is disclosed, not swept under the rug**, even
   though the balance of evidence (4/5 windows, both events, the established primary window)
   says it does not survive as a real finding.
5. **`back`'s peak-time coincidence remains an open, disclosed ambiguity** (§5.4), inherited
   unchanged from the squat build — not re-resolved here.
6. **Single clip, single subject's anthropometry** — same inherited limit as every corpus-IK cert
   in this repo.
7. **No Static Optimization / muscle-driven tier** — deliberate, per the task's own instruction and
   `force_scenes_batch.py`'s own prior finding that it is not dynamically plausible on this corpus.

## 8. Reproducing this result

```
.venv-msk/bin/python3 -u scripts/msk/plyo_jump_cascade.py
```

Re-runs the full 37-clip corpus screen (Step 1, ~15s), then the chosen clip's base battery,
regression check, ground-contact-loss handling, and propulsion-phase cascade analysis (Steps 2–6,
~25s, includes a 5-window SavGol robustness sweep). Exit code reflects the method-trustworthiness
gate (§3), not whether a cascade was found.

## Files

- `scripts/msk/plyo_jump_cascade.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py`'s proven `parse_mot`/`savgol_smooth_and_derivs`/`get_descendant_bodies`/
  `G` and `force_scenes_batch.py`'s proven Newton-reaction primitives + established constants —
  neither re-implemented).
- `data/msk_smoketest/plyo_jump_cascade/plyo_jump_cascade_results.json` — every number in this doc,
  machine-written (corpus screen rows, base battery, regression check, ground-contact-loss
  inflation, flight events, both cascade-timing conventions, the full SavGol + pre-margin
  sensitivity sweeps).
- External anchor (not modified by this build): `data/msk_smoketest/force_scenes_batch/
  jump_candidate__2026-06-30_DaMKH8aDBxu_6/scene_results.json` — produced by
  `force_scenes_batch.run_scene()` in an earlier, separate interactive run during clip selection.
- Reused, unmodified: `scripts/msk/force_scenes_batch.py`, `scripts/msk/validate_joint_force.py`,
  `data/msk_models/LaiArnoldModified2017_mediapipe_bridge_subject2_scaled.osim`,
  `data/msk_ik/2026-06-30_DaMKH8aDBxu_6/smoothed/*.mot`.
