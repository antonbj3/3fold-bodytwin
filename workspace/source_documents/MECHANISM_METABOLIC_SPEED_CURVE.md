# MECHANISM METABOLIC COST-OF-TRANSPORT vs SPEED CURVE — does the twin reproduce the U-shape, not just the level? (2026-07-21)

**Question** (`docs/MECHANISM_METABOLIC_CALORIMETRY.md` found the twin over-predicts net cost of
transport (COT) by **1.945×(SO)/2.106×(CMC)** (Umberger2010) / **1.569×/1.550×** (Bhargava2004) vs
Koelewijn 2019's in-vivo indirect-calorimetry anchor, at ONE speed, ~1.06 m/s). A single-speed
comparison cannot distinguish "the model is a wrong constant" from "the model has the wrong
physics" — a twin that reproduces the U-shaped net-COT-vs-speed curve (and its speed-of-minimum)
with an offset LEVEL is a fundamentally different, more diagnostic result than one that gets the
shape wrong too. **This doc builds the actual curve and tests its shape**, not just its level.

**Falsifier (stated by the task before this doc was assembled)**: does the twin reproduce the
U-shape + the ~1.3 m/s optimum, even if the absolute level over-predicts ~2× (already established)?

**Headline: NO on shape — forced through four decorrelated checks, all converging — but YES on a
different, genuinely positive finding: the ~1.6-1.9× over-prediction factor is remarkably STABLE
(CV 6-7%) across a 65%-wide speed range and 10 different simulated bodies.** Neither half of this
is hidden inside the other.

## Pre-registration

**C** ("twin reproduces the literature U-shape"): a pooled quadratic fit `COT(v) = a·v² + b·v + c`
across all gate-passing points is **convex** (`a>0`) for **both** Umberger2010 and Bhargava2004,
its vertex speed falls in **[0.8, 1.8] m/s** (a deliberately generous band covering both Ralston
1958's own net-COT-minimizing speed ~1.05 m/s and the task's cited ~1.3 m/s — not a knife-edge
single-point target), **and** at least half of subjects with a well-defined per-subject vertex
(`a>0`, vertex inside [0.3, 2.5] m/s) individually land in **[0.7, 2.0] m/s**.

**¬C**: fails any leg above.

**Falsifier restated as a machine gate**: `C_umberger_supported` / `C_bhargava_supported`, computed
directly in `scripts/msk/metabolic_speed_curve.py`, fixed before the fit was run.

## 0. Enumeration first (measured, not assumed)

Per the task's own instruction, walking-speed availability was **measured directly from IK
kinematics** (pelvis/whole-body-COM net displacement ÷ duration) before committing to any analysis
plan:

| check | result |
|---|---|
| subject2's own repeated walking-family trials (`walking1/2/3`, `walkingTS1/2/4`) | **6 trials, speeds 0.974–1.114 m/s** (±14% relative range) — repeats of the SAME self-selected comfortable pace, **not** a genuine speed sweep |
| Extended to all 10 subjects (`subject2`–`subject11`) × 8 candidate trial stems (`walking1-4`, `walkingTS1-4`), keeping only stems with BOTH an IK `.mot` AND a Static-Optimization activation file | **60 usable points**, speed range **0.875–1.445 m/s** (65% relative range) — a genuine, wide spread, because `walkingTS*` trials are consistently SLOWER than each subject's own `walking*` trials (cadence 1.09-1.17 Hz vs 1.26-1.32 Hz for subject2 specifically) — a real within-subject manipulation, not repeat noise, layered on top of between-subject differences in self-selected pace |
| Vendor-native Static Optimization | Every subject's `OpenSimData/Mocap/SO/<trial>_StaticOptimization_activation.sto` is dated **Jan 2022** (`stat`-verified) — predates this entire project; genuinely vendor-shipped with the external LabValidation dataset, not written by any script here |
| Vendor-native vs this repo's own independently-resolved subject2/walking1 activation | Same (158,112) shape/columns; flattened correlation **0.9999997**, mean abs diff **0.000102**, max abs diff **0.0056** across 80 muscles × 158 frames — functionally identical, safe to use directly (no re-solve) |

This is a materially different, richer picture than subject2 alone: the task's own conditional
("if multiple speeds exist... else use the model at one speed and compare shape") resolves to the
**first branch**, but via a cross-subject, not purely within-subject, route — disclosed as a
confound throughout, not hidden.

## 1. Pipeline (reuse, not re-derive)

New script `scripts/msk/metabolic_speed_curve.py` imports `metabolic_cost.py`'s own
`build_model`/`compute_kinematics`/`run_frames`/`trapz_mean`/`parse_mot` **unchanged** and
monkey-patches its module-level path constants per point — the same technique
`cross_subject_validation.py` already uses on the force-layer scripts. Zero bytes of
`metabolic_cost.py` are edited. Every point runs the SAME primary configuration as the walking1
anchor (rigid tendon, `ratio_slow_twitch=0.5`, both probes) — the force-reconstruction cross-check
and sensitivity sweeps (walking1-specific, already exhaustively run once) are not repeated per
point; this is a controlled, single-configuration speed sweep.

**Machine cross-check before trusting the new (vendor-activation) pipeline on 59 new points**: run
it on subject2/walking1 and diff against the ALREADY-PUBLISHED `metabolic_cost_results.json`
(which used our own re-solved activation):

| | new (vendor activation) | published (our re-solve) | rel. diff |
|---|---:|---:|---:|
| Umberger COT_net (J/kg/m) | 7.68346 | 7.68336 | 0.0013% |
| Bhargava COT_net (J/kg/m) | 6.19741 | 6.19768 | 0.0043% |
| speed (m/s) | 1.06473 | 1.06473 | 0.0000% |

**PASS** (tolerance 5%). All 60/60 points pass every pipeline-correctness gate (time-grid match,
no missing coordinates, peak-rotation-speed sanity <1200°/s, speed sanity 0.5-2.0 m/s, no NaN,
zero `equilibrateMuscles` exceptions). `muscle_mass_frac_ok` fails for all 60 (summed Fmax-derived
muscle mass implausibly exceeds whole-body skeletal-muscle literature bounds) — this is the SAME
already-disclosed, systemic issue `metabolic_cost.py`'s own docstring names for subject2/walking1
(confirmed here to be systemic, not newly discovered); per that script's own established
convention it is reported as an **open modeling uncertainty**, not a pipeline-correctness gate, and
does not exclude points from this analysis.

## 2. Curve-shape result — forced through four decorrelated checks

| check | Umberger2010 | Bhargava2004 |
|---|---|---|
| **(a) Raw pooled quadratic** (n=60, all subjects pooled) | **concave** (a=−3.19), vertex 0.76 m/s | **concave** (a=−4.26), vertex 1.06 m/s |
| **(b) Forced adversary — subject-demeaned** (removes each subject's own mean COT before pooling, isolating the within-subject speed effect from the between-subject economy-offset confound) | nominally convex (a=+0.093) but **quadratic term explains only 0.008% more variance than a straight line** — noise, not signal | **robustly concave** (a=−1.71), hump peaks near 1.09 m/s |
| **(c) Leave-one-subject-out stability of (b)'s vertex** | **UNSTABLE**: 7/10 refits stay convex but the vertex swings **1.82–4.58 m/s** (physiologically implausible, faster than any human walks) depending on which single subject is dropped | **STABLE in the wrong direction**: **0/10** refits stay convex — every leave-one-out refit is concave |
| **(d) Per-subject fits** (each subject its own control) | 4/10 subjects show any interior vertex at all (1.25, 1.47, 1.60, 2.37 m/s); only **3 of those 4** (1.25, 1.47, 1.60) land in the tighter [0.7,2.0] m/s corroboration band, 2.37 is an edge-case outside it | 2/10 show any interior vertex (1.23, 1.57 m/s); **both** land in-band |
| **(e) Non-parametric valley check** (no functional-form assumption: is each subject's own middle-speed COT below the mean of their slowest+fastest COT?) | **0/10 valleys**, 2/10 humps, 8/10 ties | **0/10 valleys**, 2/10 humps, 8/10 ties |
| **(f) Simple monotonicity** (per-subject speed-vs-COT correlation sign) | **8/10 negative** (COT falls as speed rises, throughout the WHOLE 0.875-1.445 m/s window) | 6/10 negative |
| **(g) Empirical ranking** (assumption-free: which points are actually cheapest?) | 8 cheapest points average **1.374 m/s** (near the fast end); dominated by WHICH SUBJECT (subject9/subject10 own 6/8 of the cheapest slots across all their own speeds) more than by speed itself | 8 cheapest average 1.282 m/s; same subject-dominance pattern |

**Pre-registered verdict**: `C_umberger_supported = False`, `C_bhargava_supported = False` — **¬C,
NOT SUPPORTED**, on the strict, pre-registered test. This is not a one-shot/lazy negative: check
(b)+(c) is the forced, adversary-strengthened variant (the OODA "Orient" step — the obvious
adversary to a messy pooled fit is the between-subject confound, and subject-demeaning is the
standard fixed-effects fix for exactly that) with an explicit robustness re-check (leave-one-out),
and checks (d)-(g) are three further, mutually-decorrelated angles (per-subject parametric,
non-parametric rank-based, and raw empirical). All four independent angles converge on the same
answer: **no robust U-shaped valley is detectable in the twin's net-COT-vs-speed relationship
within the 0.875-1.445 m/s range tested here.** Where a consistent signal exists at all (Umberger,
majority of subjects), it is a still-falling COT across the ENTIRE tested window — i.e., if this
twin has a true economical-speed minimum, the evidence here points to it lying at or above 1.445
m/s (faster than the fastest trial tested), not at the ~1.05-1.3 m/s literature range.

## 3. The level-stability result — a different, genuinely positive finding

| | Umberger2010 | Bhargava2004 |
|---|---:|---:|
| ratio to Koelewijn net-COT anchor (3.95 J/kg/m @ 1.3 m/s, reused not re-derived) — mean | **1.915** | **1.630** |
| range across all 60 points | [1.672, 2.228] | [1.462, 1.881] |
| coefficient of variation (std/mean) | **0.069** | **0.062** |

The ~1.9×/1.6× over-prediction the calorimetry doc found at ONE speed/subject **holds up, almost
unchanged, across a 65%-relative speed range and 10 different simulated bodies** — under 7%
coefficient of variation. This is itself a real, machine-checked, cross-subject/cross-speed-
replicated finding, distinct from the shape question: the twin's absolute error behaves as a
**stable, multiplicative bias** (consistent with a fixed structural cause — e.g. the already-
disclosed missing-co-contraction and Fmax-scaling-inflated-muscle-mass effects, both of which would
be expected to scale roughly proportionally, not idiosyncratically, with the overall activation
level) rather than a wildly speed-dependent or subject-dependent one. **Reported honestly as its
own outcome, not smoothed into either "shape confirmed" or "everything is noise."**

## 4. Forced adversary on the falsifier's OWN premise: is "~1.3 m/s" really the net-COT optimum?

Checked, not assumed — literature is **not unanimous** here:

| source | verified | metric | speed |
|---|---|---|---|
| Ralston HJ (1958), PMID 13610523 (secondary-sourced numbers, reused from the calorimetry doc) | title/journal/year via NCBI esummary; no MEDLINE abstract | **net**-COT-minimizing | **~1.05 m/s** |
| Ralston HJ (1958), same | same | **gross**-COT-minimizing | ~1.23 m/s |
| **Abe D, Fukuoka Y, Horiuchi M (2015)**, PLoS ONE 10(9):e0138154, PMID **26383249**, PMCID **PMC4575035** — **NEW this session**, full-text PMC-verified | live Europe PMC full-text XML fetch | **net** Economical Speed (level grade, n=11 male trained athletes) | **1.136 m/s** (4.09±0.31 km/h) |
| Abe et al. 2015, same | same | **gross** Economical Speed (level) | 1.431 m/s (5.15±0.18 km/h) |
| Browning RC, Kram R (2005), PMID 15919843 — re-verified this session (protocol detail is NEW) | live Europe PMC abstract re-fetch | preferred speed (**gross**-cost framing, their own wording) | 1.40 (obese) / 1.47 (normal-weight) m/s |
| Koelewijn et al. 2019 (primary magnitude anchor, reused) | PMC full text, prior session | fixed protocol speed (not a derived optimum) | 1.3 m/s |

**Two independent net-COT-specific estimates (Ralston secondary + Abe primary, PMC-verified this
session) cluster at ~1.05-1.14 m/s — materially below 1.3 m/s.** The ~1.3-1.47 m/s range is better
supported as the GROSS-cost-minimizing / preferred-speed regime (Ralston gross 1.23, Browning-Kram
preferred 1.40-1.47, Abe gross-ES 1.431 — three independent sources converging near 1.4-1.47).
Koelewijn's 1.3 m/s was a chosen protocol speed, not a literature-derived optimum. **This does not
change the shape-verdict** (no valley was found regardless of which target speed is used — see
§2), but it is a disclosed, symmetric-QC correction to the falsifier's own framing, exactly the same
discipline this repo's other docs already apply to task-brief citation imprecisions. Bonus
methodological cross-check: Abe et al. 2015's own Eq. 3-5 define `CoT(v)=a·v²+b·v+c` and
`ES=−b/(2a)` — **the identical quadratic-fit-plus-vertex method** this doc independently chose,
confirming the approach against the field's own established convention.

## 5. Confidence tier

Per this repo's own established ledger convention (`docs/MECHANISM_TRUST_LEDGER.md`), matching the
existing `METABOLIC_COST.md` row for this identical comparison type: **cadaveric-or-published-
plausibility** (a genuinely external, decorrelated, published human dataset; not a direct in-vivo
force/telemetry measurement) — **muscle-energetics-model-limited** (per the already-disclosed
Miller 2014 ~3× model-formulation spread) **and cross-subject-confound-limited** (54/60 points are
cross-subject, not within-subject; the subject-demeaning fix is standard but not a full substitute
for a true within-subject sweep at each of 10 bodies).

## 6. Honest gaps

1. **Range may not extend far enough past the optimum on the high-speed side.** The fastest tested
   points (1.4-1.445 m/s) are still the CHEAPEST in the raw data — if this twin's true minimum lies
   above 1.445 m/s, this dataset cannot see the rising (right) arm of its own U at all. A genuinely
   new, falsifiable prediction, not tested here: running faster trials (if any exist in this or a
   companion dataset) should show twin COT rising again.
2. **54/60 points are cross-subject**, confounded by anthropometry/individual model quirks
   (subject3/subject6 read expensive across ALL their own tested speeds; subject9/subject10 read
   cheap across all theirs) — subject-demeaning is the standard fixed-effects fix but cannot rule
   out subject × speed interaction effects a pure additive-offset model would miss.
3. **Vendor-native Static Optimization was validated for equivalence on exactly ONE point**
   (subject2/walking1, the only trial where both a vendor-native and an independently-resolved
   version exist) — assumed, not independently re-verified per subject, that the same near-identity
   holds for the other 9 subjects' vendor files.
4. **No force.sto cross-check or sensitivity sweeps repeated per point** (Steps 3+5 of
   `metabolic_cost.py`, run exhaustively once for the walking1 anchor only) — this speed curve's
   shape-verdict is conditional on the SAME single (rigid-tendon, ratio=0.5) configuration; whether
   curvature would emerge under the walking1 anchor's own corrected configuration (elastic tendon +
   literature-anchored muscle mass, shown there to shift the ratio from 1.95 to 0.93) is untested
   here — a natural next step, not run this session (would double the compute for a question not
   yet asked).
5. **`muscle_mass_frac_ok` fails for all 60 points** — the same already-disclosed, systemic
   Fmax-scaling issue named in `metabolic_cost.py`'s own docstring, confirmed here to recur across
   all 10 subjects, not newly discovered; affects the absolute LEVEL, not directly the
   curvature/shape question this doc is centrally about.
6. **Abe et al. 2015's population (n=11 male trained athletes)** differs from Koelewijn's (n=12,
   mixed-sex, general) and from this twin's own cross-subject cohort — a population-mismatch
   caveat specific to the Abe-2015-sourced economical-speed numbers.
7. **"TS" trial-naming's exact meaning was never confirmed via metadata** (`sessionMetadata.yaml`/
   `Setup.yaml` carry no explicit speed field) — inferred only from the measured, consistent
   within-subject pattern (TS trials systematically slower + lower cadence than each subject's own
   "walking" trials). A reasonable, evidence-based inference, not a confirmed protocol fact.
8. **Ralston 1958's specific numbers remain secondary-sourced** (Wikipedia, not the primary
   MEDLINE-unindexed abstract) — unchanged from the calorimetry doc's own disclosed gap, not
   re-resolved this session.

## Files

- `scripts/msk/metabolic_speed_curve.py` — new, re-runnable script (reuses `metabolic_cost.py`'s
  `build_model`/`compute_kinematics`/`run_frames`/`trapz_mean`/`parse_mot` unchanged, zero bytes of
  that file edited); enumerates 60 subject×trial points, cross-checks the new vendor-activation
  pipeline against the already-published subject2/walking1 result, runs the primary configuration
  on all 60, and performs the pooled/per-subject/subject-demeaned+leave-one-out/non-parametric
  curve-shape analysis plus the ratio-to-anchor stability check.
- `data/msk_smoketest/metabolic_speed_curve/metabolic_speed_curve_results.json` — every number in
  this doc, machine-written: all 60 raw points (subject, trial, speed, both COTs, all gates), the
  crosscheck, pooled + per-subject + subject-demeaned + leave-one-out-vertex + non-parametric-valley
  results, the ratio-to-anchor statistics, the pre-registered test verdict, and the full literature
  citation table with PMIDs/PMCIDs.
- Inputs read in place, never modified: every subject's `OpenSimData/Mocap/{IK,SO,Model}/*` under
  `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject{2..11}/` (external,
  read-only mount); `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`
  (the already-published anchor this doc cross-checks against, reloaded not recomputed).
- Reused, not modified: `scripts/msk/metabolic_cost.py`, `scripts/msk/cross_subject_validation.py`
  (path-patching technique precedent), `docs/MECHANISM_METABOLIC_CALORIMETRY.md` (the Koelewijn
  anchor value and the Ralston/Browning-Kram citation groundwork this doc builds on),
  `docs/MECHANISM_TRUST_LEDGER.md` (tier vocabulary).

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
