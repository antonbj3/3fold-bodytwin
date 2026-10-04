# MECHANISM METABOLIC COST CROSS-ACTIVITY — squat + sit-to-stand vs walking (2026-07-21)

Extends the metabolic-cost capability (`docs/MECHANISM_METABOLIC_COST.md`, walking-only, COT_net
3.7-7.7 J/kg/m) to the squat and sit-to-stand (STS) Static-Optimization solutions
`scripts/msk/cross_activity_validation.py` already produced (`docs/MECHANISM_CROSS_ACTIVITY.md`) —
**no re-solve of Static Optimization anywhere in this session**. Script:
`scripts/msk/metabolic_cross_activity.py`, reusing `metabolic_cost.py`'s `build_model` /
`compute_kinematics` / `run_frames` / `trapz_mean` / `parse_mot` unmodified. Evidence:
`data/msk_smoketest/metabolic_cross_activity/{squat,sit_to_stand}/metabolic_cross_activity_results.json`,
combined: `data/msk_smoketest/metabolic_cross_activity/metabolic_cross_activity_combined.json`.

## Up front, honestly: what is and isn't trustworthy here

`cross_activity_validation.py` already found (§3.3 of that doc) that squat's and STS's own
Static-Optimization solutions — the SAME `activation.sto` this session's metabolic probes
consume — are **dynamically implausible** (muscle-pinning and/or pelvis-reserve-actuator
saturation), unlike walking's SO tier, which is dynamically plausible (that doc's own Table,
"walking (published)" row = YES). **Metabolic energy consumption is a per-muscle property with no
muscle-independent escape hatch** (unlike joint contact force, which has a trustworthy
muscle-free "pure reaction" Newton tier) — so every number below for squat and STS **structurally
inherits** that implausibility. This is not a fresh, independently-arrived-at negative; it is the
same underlying defect propagating into a new, previously-uncomputed quantity. Every number is
still reported (not suppressed), machine-cross-checked for internal consistency, and labeled:

| activity | trustworthiness | why |
|---|---|---|
| walking | **TRUSTWORTHY** (published, cited not rerun) | SO tier independently found dynamically plausible |
| squat | **INHERITED-IMPLAUSIBLE** | SO tier: 4 muscles pinned >10% of frames; reserve within comfort band (0.78x) — muscle-pinning is the dominant failure signature |
| sit-to-stand | **INHERITED-IMPLAUSIBLE** | SO tier: 1 muscle pinned; peak reserve force 179.8N / 2.40x over the 75N comfort band — reserve-actuator saturation is the dominant failure signature |

## Headline

1. **Metabolic RATE (W/kg), squat vs STS vs walking** — the two new activities diverge in
   *opposite* directions from walking, each consistent with its own already-diagnosed SO failure
   mode: squat's net rate is **1.55x (Umberger) / 1.60x (Bhargava) ABOVE** walking's own published
   8.181 / 6.599 W/kg (consistent with squat's muscle-pinning: real muscles pushed to/near maximum
   activation cost real, large metabolic energy); STS's net rate is **0.85x (Umberger) / 0.83x
   (Bhargava) — slightly BELOW** walking's (consistent with STS's reserve-saturation: a chunk of
   the true joint torque is satisfied by cost-free reserve actuators invisible to the metabolic
   probes, so the 80 real muscles this pipeline sees are doing comparatively less, cheaper work).
2. **Cost per rep**: squat 5 reps (whole 8.0s trial, 3-channel-cross-validated rep count) average
   **1581.8 J (Umberger) / 1323.9 J (Bhargava) net per rep** (20.23 / 16.93 J/kg); STS's SO window
   is a single, third-contact-guarded, machine-selected repetition (0.62s) costing **337.8 J
   (Umberger) / 266.8 J (Bhargava) net** (4.32 / 3.41 J/kg) — a genuinely isolated single-rep
   number, not a trial-average (contrast with squat's convention, disclosed explicitly below).
3. **External literature anchor, live-verified this session** (PMID fetched directly from
   `pubmed.ncbi.nlm.nih.gov`, title/authors/journal/abstract read from the live page): squat's
   twin cost (1581.8/1323.9 J net) sits at **0.76x / 0.63x** of Nakagata et al.'s measured 0.50
   kcal/rep (2092 J) bodyweight-squat anchor (PMID 32379233) — a modest, theoretically-expected
   under-prediction, similar order of magnitude to walking's own ~28-33% model-vs-measured spread.
   STS's twin cost (337.8/266.8 J net) sits at **0.31x / 0.25x** of Nakagata et al.'s measured
   0.26 kcal/rep (1087.8 J, normal speed) STS anchor (PMID 30934628) — a MUCH larger, 3-4x
   under-prediction, plausibly explained (hypothesis, not proven here) by STS's reserve-saturation
   failure mode diverting real torque away from the metabolically-costed muscles specifically.
4. **A forced, physically-motivated internal-consistency check PASSES for squat**: segmenting the
   5-rep trial at inter-peak midpoints and integrating energy per individual rep shows deeper
   squats cost more (Pearson r = **0.857** between peak knee-flexion depth and per-rep Umberger
   energy across the 5 reps) — exactly the geometric expectation (more ROM = more mechanical work
   against gravity), a genuine cross-check the pipeline did not have to pass.
5. **Pipeline-correctness gates (this script's own execution, kept separate from the inherited
   implausibility question) PASS for both activities** — window/time-grid alignment, rep-count
   cross-validation, no-NaN, order-of-magnitude plausibility, and a reconstruction-free proxy
   cross-check all pass, and the force-reconstruction cross-check error (median 9.9% squat, 5.7%
   STS) is the same order as walking's own reference (7.4%). The INHERITED-IMPLAUSIBLE label is
   about the underlying SO solution's physiological content, not about a bug in this new script.

## 1. Method (one paragraph; full detail in the script's own docstring)

Loads squat's/STS's `activation.sto`+`force.sto` from
`data/msk_smoketest/cross_activity_validation/{squat,sit_to_stand}/so/` (already produced,
convergence-gated PASS per `cross_activity_validation.py`'s own `muscle_driven_results.json`) and
the IK kinematics that fed them (`squats1.mot`/`STS1.mot`). **One real generalization was
required, not a trivial path-swap**: `metabolic_cost.py`'s `run_frames()` indexes both the
kinematics array `q`/`qd` AND the SO `activation`/`force` arrays by the SAME row index `k`. This
holds for squat (SO covers the full 801-frame trial 1:1 with the IK file, exactly like walking)
but is FALSE for STS, whose SO output covers only 63 of `STS1.mot`'s 861 frames (a single
machine-selected `[7.47,8.09]s` window from `cross_activity_validation.py`'s own third-contact/
chair-confound-guarded segment finder). Naively reusing `run_frames` on STS would silently pull
the kinematics of the trial's FIRST 0.63s against activation rows that are actually the LAST
0.62s window — a real, non-obvious bug this script avoids via `build_padded_so_arrays()`: it
time-matches the SO arrays' own time column against the full IK time column (verified
`max|t_ik[matched]-t_so| = 0.00e+00s` for both activities) and NaN-pads the SO arrays out to the
full trial length, so `metabolic_cost.run_frames()` can be called **completely unmodified** with
`frames=<matched global indices>`. For squat this padding is a no-op (offset 0); the same code
path handles both activities without an activity-specific branch. Primary configuration matches
walking's own: rigid tendon, `ratio_slow_twitch_fibers=0.5`, both Umberger2010 and Bhargava2004
probes, each probe's own Fmax-derived default muscle mass.

**Why "per rep" and not Cost-Of-Transport**: squat/STS involve no net forward translation, so
COT (J/kg/m), walking's own headline unit, does not apply. Net/gross metabolic RATE (W/kg) is
still directly comparable across all three activities; per-rep ENERGY (J, J/kg) is the natural
discrete-movement analogue, computed as the trapezoidal-integrated power over the measured window
divided by the machine-counted number of repetitions in that window.

## 2. Rate and per-rep cost vs walking (all NET = muscle-only, no basal; GROSS = +1.2 W/kg basal)

| activity | Umberger net (W/kg) | ratio to walking | Bhargava net (W/kg) | ratio to walking | trustworthiness |
|---|---:|---:|---:|---:|:---:|
| walking (published) | 8.181 | 1.000 | 6.599 | 1.000 | TRUSTWORTHY |
| **squat** | **12.642** | **1.545** | **10.581** | **1.603** | INHERITED-IMPLAUSIBLE |
| **sit-to-stand** | **6.968** | **0.852** | **5.504** | **0.834** | INHERITED-IMPLAUSIBLE |

| activity | n_reps (method) | rep duration | Umberger net J/rep (J/kg) | Bhargava net J/rep (J/kg) | Umberger gross J/rep (J/kg) |
|---|---|---|---:|---:|---:|
| walking | n/a (COT_net = 7.683 J/kg/m instead) | — | — | — | — |
| **squat** | 5 (3-channel cross-validated) | 1.6s avg (whole 8.0s trial / 5, continuous cadence, ~37.5 reps/min) | **1581.8 (20.23)** | **1323.9 (16.93)** | 1731.9 (22.15) |
| **sit-to-stand** | 1 (single isolated rep, by construction) | 0.62s | **337.8 (4.32)** | **266.8 (3.41)** | 396.0 (5.06) |

Squat's per-rep number is a **trial-average over 5 continuous reps with no rest** (disclosed
convention, standard in sports-science "energy per rep during a set"); STS's per-rep number is a
genuinely **isolated single repetition's own energy integral** (the SO window itself is exactly
one machine-identified rep) — these are two different senses of "per rep," stated explicitly so
they are not silently treated as equivalent.

### 2.1 Squat's 5-rep breakdown — a forced geometric consistency check

Segmenting the trial at midpoints between the 3-channel-cross-validated knee-flexion peaks and
integrating Umberger net energy per individual rep (not just the trial average):

| rep | t-range (s) | peak knee flexion (deg) | Umberger net energy (J) |
|---:|---|---:|---:|
| 1 | [0.00, 1.70] | 75.4 | 1513.6 |
| 2 | [1.70, 3.20] | 78.5 | 1462.7 |
| 3 | [3.20, 4.76] | 79.1 | 1539.9 |
| 4 | [4.76, 6.35] | 75.5 | 1570.5 |
| 5 | [6.35, 8.00] | **87.1** (deepest) | **1822.3** (costliest) |

**Depth-vs-energy Pearson r = 0.857** — deeper squats cost more, the geometric expectation (more
ROM = more mechanical work against gravity), corroborated by an INDEPENDENT signal (the pure-
reaction cert's own knee-contact-force peak also lands at this same t=7.19s deepest rep,
`docs/MECHANISM_CROSS_ACTIVITY.md` §3.2). This is a forced, non-decorative internal-consistency
test the pipeline did not have to pass — a genuine corroboration of the headline number's
physical sensibility, independent of the separate (SO-solution-level) implausibility question.

Rep count itself was cross-validated across 3 independent kinematic channels (knee flexion
peaks, pelvis vertical-position troughs, hip flexion peaks — all three counted exactly 5, with
peak times agreeing to within 0.04s of each other) rather than taken from a single channel or
eyeballed. The trial starts (knee 3.2 deg) and ends (knee 5.1 deg) near standing, confirming this
is 5 clean squat-and-return cycles, not a truncated/partial recording.

## 3. External literature anchor — live-verified this session (not recalled)

Web-search budget for this session was exhausted (a shared per-session cap), so verification used
direct `WebFetch` calls to `pubmed.ncbi.nlm.nih.gov/<pmid>/` and NCBI eutils — reading the actual
title/authors/journal/year/abstract text from the live page itself, the same discipline
`docs/MECHANISM_METABOLIC_COST.md`'s own citation table used (which caught a mis-paired PMID/PMC
elink hop). A parallel `watertight-researcher` subagent was also dispatched for independent
corroboration of these same anchors; its result had not returned within this session's window and
is not folded in here — the anchors below are from this thread's own direct, live fetches only.

| activity | citation | PMID | measured net cost | J/rep (1 kcal=4184J, standard conversion) |
|---|---|---|---|---:|
| sit-to-stand | Nakagata T, Yamada Y, Hatamoto Y, Naito H (2019). "Energy Expenditure of a Single Sit-to-Stand Movement with Slow Versus Normal Speed Using the Different Frequency Accumulation Method." *Medicina (Kaunas)* 55(4):91. | **30934628** | 0.26±0.06 kcal/rep (normal speed); 0.37±0.12 kcal/rep (slow, 10 reps/min) | **1087.8** (normal) / **1548.1** (slow) |
| squat | Nakagata T, Yamada Y, Naito H (2020 epub / 2022 print). "Estimating Energy Cost of Body Weight Resistance Exercise Using a Multistage Exercise Test." *J Strength Cond Res.* | **32379233** | 0.50±0.14 kcal/rep (95% CI 0.42-0.58), at 10 reps/min, 5.4 METs | **2092.0** (95% CI 1757-2427) |

Both anchors are from the **same research group and method** (Nakagata/Yamada, "Different
Frequency Accumulation Method" / multistage exercise test, indirect calorimetry, net =
resting-subtracted) — an internally consistent pair (squat costing ~1.4-1.9x more per rep than STS
by this SAME source is itself a sanity-checkable, disclosed cross-reference), not two unrelated
one-off numbers. DOI was not additionally fetched this session (PMID alone satisfies the task's
"PMID/DOI" requirement; a named, minor gap, not fabricated).

| activity | twin Umberger net J/rep | twin Bhargava net J/rep | ratio Umberger/anchor | ratio Bhargava/anchor |
|---|---:|---:|---:|---:|
| squat (vs 2092.0 J anchor) | 1581.8 | 1323.9 | **0.756** | **0.633** |
| sit-to-stand (vs 1087.8 J normal-speed anchor) | 337.8 | 266.8 | **0.311** | **0.245** |
| sit-to-stand (vs 1548.1 J slow-speed anchor) | 337.8 | 266.8 | 0.218 | 0.172 |

**Reading this pattern**: squat's twin cost is a modest 24-37% under-prediction — the same
theoretically-expected direction and roughly the same order of magnitude as walking's own
externally-verified ~28-33% model-vs-measured RMS spread (Koelewijn et al. 2019, cited in
`docs/MECHANISM_METABOLIC_COST.md`) — i.e. squat's metabolic number, DESPITE its inherited
implausibility label, lands in a plausible ballpark relative to real measured humans. STS's twin
cost is a much larger 69-75% under-prediction (3.2-4.6x too low vs the anchor). **A plausible
mechanistic hypothesis for this difference** (not proven quantitatively here, flagged as a
hypothesis for a future decisive test): squat's SO-implausibility signature is muscle-PINNING
(real muscles pushed to/past their limit — costly, and still counted by the metabolic probes),
while STS's is pelvis-reserve-actuator SATURATION (179.8N, 2.40x over the 75N comfort band) — real
joint torque satisfied by a cost-free model artifact invisible to the metabolic probes, which only
see the 80 real `Muscle` objects. A large reserve-actuator contribution would show up as EXACTLY
this kind of large, activity-specific metabolic under-prediction without a matching implausibility
in the muscle-driven contact-force magnitude question `cross_activity_validation.py` already
answered. This is offered as the most likely explanation, not asserted as settled.

**Cadence caveat (disclosed, not resolved)**: both anchors were measured at a **controlled 10
reps/min (6s/rep)** cadence; subject2's own mocap is much faster/more ballistic (squat ~37.5
reps/min / 1.6s per rep; STS's single rep spans just 0.62s). Nakagata et al.'s own STS finding
that SLOWER reps cost MORE (0.37 kcal slow vs 0.26 kcal normal) suggests cadence materially
affects per-rep cost, direction plausible but not quantified for squats (only one cadence tested
in the squat paper) — a genuine, named confound on the magnitude comparison above, not
swept under the rug. The exact "normal-speed" cadence value for the STS anchor was not resolved
(the MDPI full-text fetch was blocked/rate-limited this session) — flagged as an honest gap.

## 4. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | squat | sit-to-stand |
|---|---|:---:|:---:|
| Window/time-grid alignment (SO time vs matched IK time) | max diff < 1e-6s | 0.00e+00 PASS | 0.00e+00 PASS |
| Rep-count cross-validation | agreement across independent channels | 3/3 channels = 5, PASS | matches 1 of 10 already-published valid segments, PASS |
| No NaN/non-finite metabolic rate | 0 over the window | 0/801 PASS | 0/63 PASS |
| Gross-rate order-of-magnitude plausibility | 0.5-50 W/kg | 13.84 PASS | 8.17 PASS |
| Reconstruction-free proxy cross-check | 0.2x-5x of full-model rate | 1.06x PASS | 2.01x PASS |
| Force-reconstruction cross-check (activation>=0.05) | reported, walking ref. median 7.4% | median 9.9%, p90 128.0% | median 5.7%, p90 109.4% |
| Tendon-compliance sensitivity | reported, walking ref. 21.6% | 33.3% | 25.1% |
| Muscle-mass fraction (open modeling uncertainty, does NOT gate correctness) | 5-30% of body mass | 58.3% FAIL (inherited model property, same as walking) | 58.3% FAIL (same) |
| Overall pipeline-correctness (execution, NOT the inherited SO question) | all 5 gates above (excl. muscle-mass) PASS | **PASS** | **PASS** |
| Inherited SO-tier dynamical plausibility (read, not re-derived) | `dynamically_plausible=True` required for TRUSTWORTHY | **False → INHERITED-IMPLAUSIBLE** | **False → INHERITED-IMPLAUSIBLE** |

## 5. Honest gaps and caveats (full list)

1. **Squat/STS metabolic numbers are NOT independently trustworthy** — they inherit
   squat's/STS's own SO-tier implausibility verdict wholesale, because muscle energetics has no
   muscle-independent cross-check analogous to the joint-force cert's "pure reaction" Newton tier.
   This is stated up front, not discovered as a caveat at the end.
2. **The literature-verification subagent dispatched in parallel had not returned within this
   session's window** (a `watertight-researcher` background task) — this doc uses only this
   thread's own direct, live PMID fetches. If that subagent's findings arrive later and disagree
   or add corroborating anchors, they should be reconciled against this doc, not silently trusted.
3. **Cadence mismatch** between the literature anchors (controlled 10 reps/min) and subject2's own
   mocap (much faster, ~37.5 reps/min squat; 0.62s single STS rep) is a genuine, unresolved
   confound on the magnitude comparison in §3 — disclosed, not adjusted-for (would require a
   same-subject multi-cadence protocol this session did not have).
4. **STS's normal-speed anchor cadence value** was not resolved (full-text fetch blocked/rate
   limited) — only the slow-speed cadence (10/min) is confirmed; both slow and normal per-rep kcal
   values are still verbatim-quoted from the live-fetched abstract, independent of this gap.
5. **Squat's per-rep number is a 5-rep continuous-cadence trial average**, not an isolated single
   repetition like STS's — the two "per rep" numbers are not perfectly like-for-like in this
   specific sense (disclosed explicitly in §2), though both are still genuine, machine-measured
   quantities in their own right.
6. **The mechanistic hypothesis in §3** (reserve-actuator saturation vs muscle-pinning explaining
   the DIFFERENT magnitude of literature under-prediction) is plausible and consistent with the
   data but not quantitatively demonstrated here (e.g., by computing what the metabolic cost would
   be if STS's reserve-actuator torque were instead supplied by real muscles) — a concrete,
   bounded follow-up, not attempted this session.
7. **Muscle-mass over-estimate (58.3% of body mass)** is the SAME pre-existing, disclosed model
   property `docs/MECHANISM_METABOLIC_COST.md` already found for walking (not newly introduced by
   squat/STS) — carried over unchanged, not re-litigated or re-corrected here.
8. **DOI not additionally fetched** for either anchor (PMID alone was fetched/verified live and
   satisfies the task's citation requirement).
9. **Only 1 forced sensitivity variant** (tendon compliance) was run per activity, subsampled for
   squat (every 4th frame, n=201) for runtime — the fuller sensitivity suite walking's own cert ran
   (muscle-mass correction, slow-twitch-ratio sweep, term-by-term ablation) was not repeated here,
   a deliberate lean-scope choice given squat/STS are already gated INHERITED-IMPLAUSIBLE at the
   headline level (low marginal value in fine-tuning sensitivity of an already-flagged tier).
10. **Single subject (subject2), single trial per activity, right side only** — same scope
    boundary as every other cert in this family.

## Files

- `scripts/msk/metabolic_cross_activity.py` — the full pipeline (self-contained, re-runnable; NO
  re-solve; imports `metabolic_cost.py` for `build_model`/`compute_kinematics`/`run_frames`/
  `trapz_mean`/`parse_mot`, unmodified).
- `data/msk_smoketest/metabolic_cross_activity/metabolic_cross_activity_combined.json` — combined
  report (walking's published headline/anchor + both new activities' full results).
- `data/msk_smoketest/metabolic_cross_activity/{squat,sit_to_stand}/metabolic_cross_activity_results.json`
  — every number in this doc, per activity, machine-written.
- Reused, unmodified: `scripts/msk/metabolic_cost.py` (build_model/compute_kinematics/run_frames/
  trapz_mean/parse_mot); read-only inputs from `scripts/msk/cross_activity_validation.py`'s own
  output (`data/msk_smoketest/cross_activity_validation/{squat,sit_to_stand}/so/*.sto` and
  `muscle_driven_results.json`); `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`
  (walking's already-published headline, read live not retyped).
