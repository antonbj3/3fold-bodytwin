# MECHANISM EMG TIMING VALIDATION — decorrelated muscle-activation anchor (2026-07-21)

Closes the honest gap every prior MSK force cert flagged identically: `docs/MECHANISM_STATIC_OPT.md`,
`MECHANISM_HIP_FORCE.md`, and `MECHANISM_ANKLE_FORCE.md` all disclose *"the mechanism is SO's
activation-minimizing solution, not measured EMG"* as an unclosed limitation. This session adds the
independent check: does the twin's **already-computed** Static-Optimization (SO) muscle-activation
time series (subject2 `walking1`, reused from `data/msk_smoketest/subject2_walking1/static_optimization/
so/` — **zero re-solving**) turn the right leg muscles ON/OFF at the right point in the gait cycle,
compared to published surface EMG for walking? Every number below is machine-computed this session
(`scripts/msk/validate_emg_timing.py`, exit 0). Isolation respected: `.venv-msk` only, reads existing
data in place, no git commit/push.

## Headline result

**Pre-registered gate (symmetric Jaccard vs a circular-shift null, §5): 0 PASS / 3 PARTIAL / 5 FAIL**
out of 8 muscle groups. Read alongside the required decomposition (§6), the picture is more specific
than that count alone suggests:

| finding | muscles | what it means |
|---|---|---|
| **Correct phase, shorter duration than EMG** (precision = 1.00 — twin is *never* active at an implausible phase, just for less of the cycle) | gluteus maximus, hamstrings, vasti | This is the pre-registered, EXPECTED signature of SO's activation-minimizing objective (§0) — not a phase error. |
| **Independently corroborated against the literature's OWN simulation-vs-EMG mismatch** | gluteus maximus, gluteus medius | Rajagopal et al. 2016's own CMC-based simulation showed the SAME qualitative disagreement vs their own measured EMG for these two exact muscles (quoted verbatim, §3) — this twin's SO reproduces it independently. |
| **Best match, boundary-line miss** | gastrocnemius | rank = 94.9% (tied for the single BEST possible phase-alignment among all 100 possible shifts; pre-registered bar was 95.0%). |
| **Specific literature-predicted failure mode, explicitly checked, found ABSENT** | hamstrings (mid-swing) | Rajagopal 2016 flagged "compensatory hamstring activity during swing in simulations not observed in EMG" as a known defect of activation-minimizing methods — measured directly in this twin's mid-swing window: raw activation never exceeds 0.017 (floor = 0.010), i.e., essentially OFF. Does not replicate. |
| **Genuine partial mismatch, moderate confidence in the literature anchor itself** | soleus, tibialis anterior, rectus femoris | Real, disclosed disagreements; rectus femoris carries the WEAKEST literature anchor of the eight (§3). |

Mean precision across all 8 muscles = **0.809**; mean recall = **0.631**; mean Jaccard = **0.514**; mean
null-rank = **77.0%** (all machine-computed, `emg_timing_results.json`). Consistent with the
pre-registered expectation: when the twin's SO fires a muscle, it is in a literature-plausible phase
81% of the time on average — it just doesn't sustain across as much of the cycle as real EMG does.

## 0. Why TIMING, not amplitude, is the fair test (pre-registered, stated before running)

SO minimizes `sum(activation^2)` subject to reproducing the net joint moments — it structurally
**under-predicts co-contraction and antagonist tone** (this is the exact gap being closed). An
amplitude mismatch is therefore expected, not a failure. ON/OFF **timing** is the part SO's optimality
argument does not get a free pass on: if the net joint moment requires a knee-extensor moment in early
stance, *some* knee extensor must turn on then, regardless of the objective function's amplitude bias.
This script scores **timing overlap only** — amplitude is never part of the verdict.

## 1. Method (pointer — full derivation in the script's own docstring)

1. **Reuse** `walking1_StaticOptimization_activation.sto` (158 frames, t∈[0,1.57s]) — no new SO solve.
2. **Detect the gait cycle from the RAW GRF**, not assumed (§2) — a genuine sub-problem this session
   had to solve, since the 1.57s trial does not contain two full right-heel-strike-to-heel-strike
   events (it starts mid-stance).
3. **Bin** each muscle group's activation onto a 100-bin (1%-GC) circular grid; normalize per muscle to
   its own `[floor, peak]` range (floor ≈ 0.01, SO's numerical lower bound, not physiological zero);
   "ON" = normalized activation ≥ a pre-registered fraction of that muscle's own dynamic range.
4. **Score** via circular-bin Jaccard overlap against the literature ON-mask, tested against an
   **exhaustive 99-shift circular-permutation null** of the twin's own mask (preserves its own duty
   cycle exactly — directly answers "is this overlap more than the two windows just both being big?").
5. **Sensitivity**: repeat at 3 activation-ON thresholds (15/25/35% of dynamic range) × 3 GRF
   stance-thresholds (2/5/10% BW) = 9 combinations; report which verdicts are threshold-stable.

## 2. Gait-cycle detection — measured, then cross-checked (not assumed)

The trial starts already in right-foot stance, so the right heel-strike that began it is **before**
t=0 (unobserved). Recovered via the standard contralateral-symmetry convention (left heel-strike = 50%
GC of the right-referenced cycle), then **cross-checked against an independent relation** — this is
the forced adversary for the %GC axis itself, not assumed correct:

| quantity | value (primary, 5%BW threshold) |
|---|---:|
| Right heel-strike (measured) | t = 1.2535 s |
| Right toe-off (measured) | t = 0.7925 s |
| Left heel-strike (measured) | t = 0.5885 s |
| Left toe-off (measured) | t = 1.4570 s |
| Cycle duration T (from left-HS-at-50%GC) | **1.3300 s** |
| Inferred prior right heel-strike (right_hs − T) | t = −0.0765 s |
| Right stance fraction (derived from T) | 65.34 %GC |
| Left stance fraction (**independently** derived, different event pair, same T) | 65.30 %GC |
| **Consistency check: right vs left stance-fraction agreement** | **diff = 0.04 percentage points** (gate ≤ 5.0) → **PASS** |

This is a genuine, non-tautological over-determination: two DIFFERENT pairs of measured GRF events
(right HS + right TO; left HS + left TO) independently triangulate the same stance-fraction to within
0.04 points. Stable across all 3 GRF thresholds tested (T = 1.327–1.333 s, <0.5% spread; §4 sensitivity
table in the JSON).

**Independent physiological sanity check on the %GC axis itself** (a hardwired law of human gait, not
a soft convention): peak plantarflexor activation (soleus 55%GC, gastrocnemius 46%GC, §6) occurs
immediately **before** the measured right toe-off (65.3%GC) — exactly where ankle push-off power must
peak in every human gait cycle. This could not happen by coincidence if the %GC axis were mis-scaled,
shifted, or phase-flipped; it is independent supporting evidence the axis itself is correctly built,
separate from the muscle-by-muscle literature-agreement scoring below.

## 3. Literature anchor — verified live this session, not recalled

**WebSearch quota was already exhausted this session** (0/2000 remaining, pre-existing per
`docs/MECHANISM_MSK_BUILD_PLAN.md`). Citations below were verified via `curl` to NCBI E-utilities
(`esearch`/`efetch`/`esummary` — a database API, not a "web search" call) and `WebFetch` on PubMed;
the Rajagopal 2016 full text was fetched and its actual prose read (never eyeballed as a figure).

- **PRIMARY**: Rajagopal A, Dembia CL, DeMers MS, Delp DD, Hicks JL, Delp SL. "Full-Body
  Musculoskeletal Model for Muscle-Driven Simulation of Human Gait." *IEEE Trans Biomed Eng.*
  2016;63(10):2068-79. **PMID 27392337, DOI 10.1109/TBME.2016.2586891, PMC5507211** (full text
  fetched: `efetch.fcgi?db=pmc&id=5507211`). This is the **direct ancestor model lineage** of the
  `LaiArnoldModified2017` model this repo's SO run uses, and its own Results section reports a
  walking EMG-vs-simulated-activation comparison (their Fig. 5, their own recorded 10-muscle surface
  EMG) **in prose**, quoted verbatim per muscle in the script (`EMG_LITERATURE` dict). Direct quote:
  *"the timing of the simulated muscle activations represented some of the major features of the
  measured EMG signals... hip extensor (gluteus medius and biceps femoris long head) and knee
  extensor (vastus lateralis) activity in early stance, plantarflexor (gastrocnemius and soleus)
  activity in late stance prior to toe-off, and tibialis anterior activity in early stance and swing,
  and hamstring (biceps femoris long head) activity in terminal swing prior to foot strike."*
- **SECONDARY** (cited by Rajagopal 2016 itself as ref [50]): Perry J. *Gait Analysis: Normal and
  Pathological Function.* SLACK Incorporated; 1992 — used only for the standard
  stance(~60%)/swing(~40%) phase-percentage convention (Rancho Los Amigos 8-phase model) needed to
  turn Rajagopal's qualitative phase words into numeric windows, since their prose does not give exact
  percentages. Disclosed per-muscle "confidence" field distinguishes direct quotes from this
  convention-based operationalization.
- **CORROBORATING** (existence + abstract verified, PMID/DOI-checked; NOT used for numeric windows —
  no per-muscle percentages were accessible from the fetched text): Winter DA, Yack HJ. "EMG profiles
  during normal human walking: stride-to-stride and inter-subject variability." *Electroencephalogr
  Clin Neurophysiol.* 1987;67(5):402-11. **PMID 2444408, DOI 10.1016/0013-4694(87)90003-4.** This
  confirms "Winter's EMG atlas" (named in the task) is real; its abstract's claim that soleus/tibialis
  anterior/gastrocnemii are the most active, least-variable distal muscles is consistent with (not
  contradicted by) the windows used here.
- Kadaba et al. 1989 (PMID 2795325, *J Orthop Res* 7(6):849-60) was also checked — real, verified, but
  its abstract (the only text this session could access) reports EMG *repeatability* statistics, not
  per-muscle on/off percentages, so it was **not usable** as a numeric anchor and is not relied upon.

**Honesty note on rectus femoris** (weakest-anchored muscle in this report): Rajagopal 2016's own
prose only states RF timing explicitly for the RUNNING condition ("quadriceps (rectus femoris and
vastus lateralis) activity in early stance"); RF is shown in their Fig. 5a for walking but not
separately described in text for that condition. The window used here (early stance + pre-swing/
initial-swing) is the STANDARD Perry-1992 two-burst convention, not a walking-specific verified quote
— flagged "low-moderate confidence" throughout, not silently treated as equal-strength evidence.

## 4. Per-muscle results (primary settings: 25% activation-ON threshold, 5%BW stance threshold)

| muscle | twin ON (%GC) | literature ON (%GC) | Jaccard | null rank | precision | recall | twin peak | verdict | lit. confidence |
|---|---|---|---:|---:|---:|---:|---:|---|---|
| Gluteus maximus | 5–23 | 85–100 + 0–50 | 0.277 | 52.5% | 1.00 | 0.28 | 12%GC | FAIL | moderate |
| Gluteus medius | 5–57 | 0–30 | 0.439 | 67.7% | 0.48 | 0.83 | 16%GC | FAIL | high |
| Hamstrings | 89–100 + 0–13 | 85–100 + 0–30 | 0.533 | 78.8% | 1.00 | 0.53 | 0%GC | FAIL | high (bflh) |
| Vasti/quadriceps | 6–22 | 0–30 | 0.533 | 85.9% | 1.00 | 0.53 | 12%GC | **PARTIAL** | high (vaslat) |
| Rectus femoris | 17–66, 67–77 | 0–30, 50–73 | 0.455 | 77.8% | 0.59 | 0.66 | 23%GC | FAIL | **low-moderate** |
| Gastrocnemius | 30–65 | 30–60 | 0.857 | **94.9%** | 0.86 | 1.00 | 46%GC | **PARTIAL** | high |
| Soleus | 22–24, 46–62 | 30–60 | 0.412 | 81.8% | 0.78 | 0.47 | 55%GC | **PARTIAL** | high |
| Tibialis anterior | 5–11, 17–46, 67–100 | 60–100 + 0–30 | 0.605 | 76.8% | 0.76 | 0.74 | 6%GC | FAIL | high |

`precision` = of the twin's own ON-time, the fraction that falls in a literature-plausible phase (a
narrowing-only mismatch scores 1.00 here). `recall` = of the literature's claimed window, the fraction
the twin actually covers (this is the one sensitive to SO's expected narrowing bias). Both are reported
**in addition to**, not instead of, the pre-registered Jaccard gate — added after observing that
Jaccard alone conflates "narrower duration than EMG" (expected, §0) with "active in the wrong phase"
(a real disagreement); the gate itself was **not** changed post-hoc (full derivation + the proof that
Jaccard/precision/recall share an identical null-rank ordering for fixed mask sizes: script
`circular_null` docstring).

## 5. Pre-registered gates (stated before Step 3 of the first run)

| gate | threshold | verdict rule |
|---|---|---|
| Gait-cycle consistency (right vs left stance-fraction) | diff ≤ 5.0 percentage points | PASS (0.04 pts, §2) |
| Muscle has real dynamic range | peak − floor ≥ 0.02 | 8/8 muscles pass (no NO-SIGNAL case) |
| Per-muscle verdict | Jaccard ≥ 0.15 **and** null-rank ≥ 95% → PASS; Jaccard ≥ 0.075 and rank ≥ 80% → PARTIAL; else FAIL | 0 PASS / 3 PARTIAL / 5 FAIL |

## 6. Two specific, falsifiable checks — both machine-measured, not narrated

1. **Does the twin replicate the literature's OWN documented simulation-vs-EMG disagreement for
   gluteus maximus/medius, or invent a different one?** Rajagopal 2016's CMC simulation showed
   "gluteus maximus only in early stance" (vs EMG's "relatively constant") and "gluteus medius...
   throughout stance" (vs EMG's "primarily early stance"). This twin's SO independently shows: glmax
   ON only 5–23%GC (a narrow early-stance burst, precision=1.00 vs the broad literature window) and
   glmed ON 5–57%GC, i.e. **87.2% of the measured stance phase** (0.57/0.6534) — both **replicate the
   same qualitative direction of mismatch**, from a completely independent generative method (SO here
   vs CMC there) and an independent dataset. This is evidence the mismatch is a structural property of
   activation-minimizing muscle-driven simulation in general, not a bug specific to this twin or this
   trial.
2. **Does the twin show the SPECIFIC "compensatory mid-swing hamstring activity not seen in EMG"
   defect Rajagopal 2016 flagged in their own CMC results?** Measured directly: hamstrings' raw SO
   activation across the entire mid-stance-through-mid-swing window (bins 15–79%GC, clear of both
   literature bursts) never exceeds **0.017** (floor = 0.010, i.e., essentially pinned at the
   numerical floor). **Does NOT replicate** — a specific, literature-predicted failure mode was
   checked for and found absent. Plausible (not further tested this session) mechanistic reason: SO
   optimizes each time frame independently with no forward-dynamics state to track, unlike CMC's
   tracking-controller machinery, which may need anticipatory muscle activity that SO's per-instant
   formulation has no analogous need for.

## 7. Sensitivity / robustness (9 threshold combinations: 3 activation-ON × 3 GRF-stance)

| muscle | primary verdict | across all 9 combos |
|---|---|---|
| Gluteus maximus | FAIL | **STABLE** (FAIL in all 9) |
| Gluteus medius | FAIL | **STABLE** (FAIL in all 9) |
| Vasti | PARTIAL | **STABLE** (PARTIAL in all 9) |
| Hamstrings | FAIL | varies: FAIL/PARTIAL |
| Rectus femoris | FAIL | varies: FAIL/PARTIAL |
| Gastrocnemius | PARTIAL | varies: PARTIAL/**PASS** (some combos clear the 95% bar) |
| Soleus | PARTIAL | varies: FAIL/PARTIAL |
| Tibialis anterior | FAIL | varies: FAIL/PARTIAL |

The two clearest-cut, robust findings (stable across every threshold tested) are: gluteus maximus's
narrow-early-stance-only pattern (a real, literature-corroborated mismatch), and vasti's clean
early-stance PARTIAL match. Gastrocnemius is a genuine boundary case — its primary rank (94.9%) misses
the pre-registered 95.0% bar by five hundredths of a percentage point (tied for the single best
possible alignment among all 100 shifts), and several nearby sensitivity settings do cross it.

## 8. Honest gaps (full list)

1. **Rectus femoris's literature anchor is the weakest of the eight** (§3) — the window used is a
   standard convention, not a walking-specific verified quote from the fetched source. Its FAIL
   verdict should be read with lower confidence than the others.
2. **The activation-ON threshold (25% of each muscle's own [floor,peak] range) is a convention, not an
   MVC-calibrated EMG onset criterion** — no real EMG amplitude exists to calibrate a model activation
   against. Swept 15/25/35% (§7); most verdicts are stable, several are not (disclosed per-muscle).
3. **The stance/swing phase-percentage convention (60%/40%) is a population average, not fit to this
   subject** — only the gait-cycle DURATION (T) and its origin (right_hs_0) are subject-specific/
   measured; the mapping from "early stance"/"late stance"/"swing" words to numeric windows uses the
   standard convention. A sensitivity check using this subject's own measured stance fraction (65.3%,
   vs the 60% convention) instead was not run this session (small, ~5-point shift; unlikely to flip
   the STABLE verdicts in §7, not verified).
4. **Jaccard as the primary gate structurally penalizes SO's known narrowing/co-contraction-
   suppression bias more than a phase-only test would** (§4) — mitigated by reporting precision/recall
   alongside, but the pre-registered PASS/PARTIAL/FAIL gate itself was deliberately NOT swapped to a
   precision-based rule after seeing the data (that would be post-hoc metric shopping); readers should
   weight the §4 decomposition, not the raw verdict label alone, for muscles with precision=1.00.
5. **Multi-compartment muscles (glmax1/2/3, glmed1/2/3) are averaged into one signal per group**, to
   match what a single surface-EMG electrode over one anatomical muscle belly would see — sub-
   compartment-specific (e.g., fine-wire) EMG could disagree with this averaged comparison; not tested.
6. **Single trial, right leg only** (subject2 `walking1`) — same scope caveat as every prior MSK force
   cert in this repo; no claim of generality across subjects, trials, speeds, or the left leg.
7. **This is a TIMING-only validation, by design** (§0) — amplitude/co-contraction magnitude is
   explicitly not scored; a full closure of the original "not measured EMG" gap would also want a
   co-contraction-INDEX-style amplitude comparison, which is a different (harder, MVC-calibration-
   dependent) question not attempted here.
8. **Kadaba et al. 1989 (PMID 2795325) could not be used quantitatively** — only its abstract was
   accessible this session (WebSearch quota exhausted, full text is Wiley-paywalled and pre-dates the
   PMC mandate); it is cited as a real, verified, existing reference only, not as a numeric source.

## 9. Next step

If sub-compartment-level EMG or a real MVC-calibrated amplitude comparison is later wanted, this
script's `MUSCLE_GROUPS`/`EMG_LITERATURE` dicts are the natural extension point. The clearest
remaining ambiguity is rectus femoris (weakest literature anchor, §3/§8.1) — a targeted fetch of a
walking-specific (not running-only) verified RF timing source would sharpen that one muscle's verdict
without re-running anything else.

## Files

- `scripts/msk/validate_emg_timing.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py` for the proven `.mot`/`.sto` parser and paths only, not re-implemented).
- `data/msk_smoketest/subject2_walking1/emg_timing_validation/emg_timing_results.json` — every number
  in this document, machine-written (gait-cycle detection at all 3 GRF thresholds, per-muscle scores,
  full 9-combination sensitivity sweep, all pre-registered thresholds).
- Reused, not re-solved: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_activation.sto`.
