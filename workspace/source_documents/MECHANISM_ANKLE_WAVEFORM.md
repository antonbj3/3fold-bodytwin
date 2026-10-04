# MECHANISM ANKLE (TALOCRURAL) CONTACT-FORCE WAVEFORM — does the push-off mechanism resolve at the joint that carries it? (2026-07-21)

Extends `docs/MECHANISM_CONTACT_WAVEFORM.md` (knee/hip contact force over-predicts, CONCENTRATED at the
push-off hump — SO-knee ratio 1.46x/Z=2.67 vs weight-acceptance 1.09x/Z=0.97) and
`docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md` (the ankle-level Achilles-equivalent tendon force at that SAME
instant is only modestly elevated, +37.1%/Z=1.06, and the plantarflexors' own natural peak arrives later,
t=0.59s, than the knee/hip peak, t=0.51s) down to the ankle (talocrural) joint's own full contact-force
WAVEFORM — the most direct joint-level view of the push-off mechanism, since this is the joint the
plantarflexors actually act across.

**Headline: outcome B, not A. The ankle contact force does NOT over-predict at push-off — it already
MATCHES its (weak-tier) cadaveric anchor there (ratio 1.016-1.088x depending on anchor construction, both
comfortably inside the pre-registered 20%-materiality band), and a full-resolution scan finds ZERO %GC
where the model/anchor ratio exceeds that band anywhere across the entire gait cycle.** This further
LOCALIZES the knee/hip push-off excess to the knee-crossing/biarticular routing rather than ankle-level
plantarflexion demand — corroborating, not contradicting, the plantarflexor thread's own finding that the
Achilles force itself is only modestly (not dominantly) elevated. Three convergent, machine-checked findings
support this reading: **(1)** the ankle contact-force curve shows NO genuine early/weight-acceptance hump
in either solver (prominence 0.02-8.9% of the curve's own peak, all far below the pre-registered 15% bar,
robust to smoothing in both directions) — a single dominant late-stance peak, exactly the shape the
cadaveric literature itself describes ("rose... during the **latter** part of stance phase", Stauffer 1977)
and qualitatively different from knee/hip's genuine double hump; **(2)** the ankle's own peak instant
(t=0.58s, 50.1%GC) sits within 10ms/0.75%GC of the plantarflexors' own natural peak (t=0.59s, 50.8%GC,
already established) and 70ms/5.3%GC away from the knee/hip peak (t=0.51s, 44.8%GC) — the pre-registered
phase offset, confirmed; **(3)** the ankle contact-force waveform correlates strongly with the
Achilles-equivalent tendon-force waveform (r=0.953 full-cycle, r=0.930 in the push-off window) — though this
third piece is the weakest discriminator on its own (ankle-vs-knee contact force is nearly as high, r=0.934,
since all gait-driven quantities share a common loading envelope; flagged, not oversold).

**Confidence tier: cadaveric/published-plausibility** (NOT in-vivo-anchored — no instrumented-ankle program
has ever existed; re-verified absent again, consistent with `docs/MECHANISM_ANKLE_FORCE.md` §0).

## Headline numbers

| quantity | value | anchor / comparison | ratio | verdict |
|---|---:|---:|---:|---|
| SO ankle push-off window (30-62%GC) peak | 484.32 %BW @ 50%GC | blended anchor 476.67 %BW (n=3, incl. Giddings talocalcaneal-adjacent) | **1.016x** | not material (<1.20 bar) |
| SO ankle push-off window peak | 484.32 %BW @ 50%GC | talocrural-only anchor 445.00 %BW (n=2, Stauffer+Procter&Paul, stricter) | **1.088x** | not material |
| SO ankle weight-acceptance window (0-30%GC) peak [phase-mismatched anchor, caveated] | 309.90 %BW @ 30%GC | blended anchor | 0.650x | — (no separate WA anchor exists) |
| SO ankle weight-acceptance window peak [phase-mismatched] | 309.90 %BW @ 30%GC | talocrural-only anchor | 0.696x | — |
| Full-waveform material-excess (>1.20x) band, either anchor | **NONE FOUND** (0/101 grid points, both anchors) | — | — | **no material excess anywhere in the cycle** |
| Genuine early-hump prominence, SO (smoothed / raw) | 0.04% / 2.2% of own peak | pre-registered bar: 15% | — | **NO genuine hump** (both far below bar) |
| Genuine early-hump prominence, CMC (smoothed / raw) | 0.5% / 8.9% of own peak | pre-registered bar: 15% | — | **NO genuine hump** (raw is the closer call, still below) |
| Ankle peak instant vs Achilles-equiv peak instant | t=0.58s vs t=0.59s | offset −0.75%GC (−10ms) | — | matches plantarflexor timing |
| Ankle peak instant vs knee/hip peak instant | t=0.58s vs t=0.51s | offset +5.28%GC (+70ms) | — | offset from knee/hip timing |
| Ankle-vs-Achilles-equiv waveform correlation | r=0.953 (full), 0.930 (push-off window) | — | — | strong, but not sharply discriminating (see below) |
| CMC ankle peak vs SO ankle peak | 405.54 vs 484.48 %BW | ratio 0.837 | — | CMC LOWER than SO — opposite of established knee/hip CMC>SO pattern |

All numbers machine-computed by `scripts/msk/ankle_waveform_analysis.py`, written to
`data/msk_smoketest/subject2_walking1/ankle_waveform/ankle_waveform_results.json` — none transcribed from
console prose.

## 1. Pre-registration (stated in the script's module docstring before the final numbers were computed)

**Three separable questions, each with a pre-registered test:**

- **Q1 (shape):** does the ankle contact-force curve show a genuine two-hump structure like knee/hip, or a
  single dominant late-stance peak? Test: a candidate early local maximum (on the smoothed curve) strictly
  before the overall peak, followed by a real dip (also strictly before the overall peak, so the main peak's
  own trailing edge into swing can never be mis-scored as "the dip after an early hump") with prominence
  **≥15% of the curve's own overall peak value** — set deliberately above the already-measured
  smoothing-sensitivity noise ceiling in this repo (`docs/MECHANISM_CONTACT_WAVEFORM.md` §3.2: hump values
  shift only 3-9% under an identical 5-point moving average).
- **Q2 (timing):** does the ankle peak's own instant land closer to the plantarflexors' natural peak
  (t=0.59s) or the knee/hip peak (t=0.51s)? Purely model-internal, no anchor uncertainty.
- **Q3 (magnitude vs anchor):** is the ankle's push-off-window value materially (>20%, this repo's
  established bar) above its cadaveric anchor, while any early-stance value stays non-material — or does
  push-off already match?

**Decision rule (pre-registered):** ratio_pushoff > 1.20 AND weight-acceptance non-material →
**outcome A** (ankle also over-predicts, concentrated at push-off — corroborates gastroc/soleus over-force
at the joint level). ratio_pushoff in [0.80, 1.20] → **outcome B** (ankle matches its anchor — localizes the
knee/hip excess to knee-crossing/biarticular routing instead, consistent with the plantarflexor thread's own
"modest, non-dominant" finding). Both outcomes were named as legitimate by the task before computing;
**measured: outcome B.**

**Forced adversary #1 (caught before trusting the anchor numbers):** an early version of the talocrural-only
anchor filter string-matched `"talocrural" in quantity_text`, which silently DROPPED Stauffer 1977 (whose own
quantity string reads "ankle joint compressive force," not literally "talocrural," even though
`docs/MECHANISM_ANKLE_FORCE.md` §1 explicitly classifies Stauffer as one of "the two studies that measure the
talocrural joint itself"). This collapsed the talocrural-only anchor to n=1/mean=390 instead of the
established n=2/mean=445 — caught by cross-checking the computed anchor against the established doc's own
prose before trusting it, fixed via an explicit PMID allowlist (912978, 7174695) rather than a fragile
substring match. Confirmed fixed: `talocrural_only_n=2`, `talocrural_only_mean_pct_bw=445.00`, matching
`docs/MECHANISM_ANKLE_FORCE.md` exactly.

**Forced adversary #2 (Q1, the "smoothing could hide a real hump" direction):** smoothing before testing for
a genuine hump risks the opposite failure mode from the one it guards against — it could suppress a REAL
early hump present in the raw per-frame solve, not just avoid manufacturing a fake one from frame jitter. The
identical genuine-hump test was therefore also run on the RAW (unsmoothed) curve for both solvers — the
condition most favorable to detecting a hump. Result: SO raw prominence 2.2% (still far below the 15% bar);
CMC raw prominence 8.9% — the closest call in the whole analysis, but still 40% relatively below the bar, not
a knife-edge. The "no genuine early hump" finding survives this adversary in both solvers.

## 2. Method

**Reused, unedited, zero new OpenSim run and zero new code duplication:** `scripts/msk/ankle_waveform_analysis.py`
imports `scripts/msk/contact_waveform_analysis.py` directly (module-level code there is import-safe — only
stdlib/numpy/openpyxl imports and function/constant definitions, no I/O at import time, `main()` guarded by
`__name__` check) rather than re-implementing its parser, %GC phase mapping, or hump-finder a second time, per
the task's own instruction to reuse that machinery. Specifically reused: `parse_mot`, `model_pct_gc_mapping`
(same `T_stride=1.327s`, `R_heel_strike_t=1.2425s` convention as the knee/hip cert — so all %GC values here
are directly comparable to that doc's own numbers), `verify_phase_alignment`, `model_curve_on_grid`,
`find_two_humps`, `smooth5`, and the constants `BW_N=766.88003`, `SO_JR`, `CMC_JR`, `GRID_N=101`.

**Ankle contact-force column:** `ankle_r_on_talus_r_in_talus_r_{fx,fy,fz}` — confirmed present (via direct
header grep) in BOTH `static_optimization/jr/walking1_JointReaction_ReactionLoads.sto` (SO) and
`cmc_second_solve/jr_cmc/walking1_JointReaction_ReactionLoads.sto` (CMC), the SAME already-committed,
already-certified files `contact_waveform_analysis.py` already parses for knee/hip — these columns were
simply never extracted before (the ankle cert, `docs/MECHANISM_ANKLE_FORCE.md`, only ever pulled the SINGLE
peak row, not the full time series). No new OpenSim run. Confirmed bug-free at this column: the ankle-force
cert's sign-bug fix (`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`) affected only the SELF-COMPUTED
(R−Σmuscle-crossing) candidate, never the official `opensim.JointReaction` output used here, which the
corrected self-computed number now matches to 0.0014% anyway.

**Achilles-equivalent tendon-force curve:** `gasmed_r+gaslat_r+soleus_r` SO tension, summed, resampled onto
the identical %GC grid — reads directly from `static_optimization/so/walking1_StaticOptimization_force.sto`,
confirmed to share the IDENTICAL 158-row/0.01s-step time grid as the JointReaction `.sto` files
(`np.allclose(t_force, t_jr) == True`, checked before trusting any cross-file resampling).

**Anchor:** the SAME three cadaveric/model-literature studies already PMID/DOI-verified live in
`docs/MECHANISM_ANKLE_FORCE.md` — pulled directly from that doc's own already-committed evidence JSON
(`ankle_force_validation_results.json`'s `published_ankle_hindfoot_anchor` field), not retyped from prose, to
avoid a transcription-drift error: Stauffer, Chao & Brewster 1977 (PMID 912978, 500 %BW, "ankle joint
compressive force... rose to about 5 times body weight during the latter part of stance phase"); Procter &
Paul 1982 (PMID 7174695, 390 %BW, talocrural joint resultant force, 3-D force-equilibrium model); Giddings et
al. 2000 (PMID 10731005, 540 %BW, talocalcaneal — anatomically adjacent, not the identical joint, disclosed).

**Self-consistency gates (machine-asserted, would raise not silently pass on mismatch):** SO ankle peak
re-derived here = 484.480 %BW @ t=0.580s vs `docs/MECHANISM_ANKLE_FORCE.md`'s published 484.48 @ t=0.58 —
**0.0000% relative difference, exact instant match.** Achilles-equivalent peak re-derived here = 320.265 %BW
@ t=0.590s vs `docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md`'s published 320.2649 @ t=0.59 — **0.0000% relative
difference, exact instant match.** Both PASS before anything new was built on top of them.

**Environment:** identical convention to `contact_waveform_analysis.py` — run with `/usr/bin/python3` (numpy
2.2.6, openpyxl 3.1.5 confirmed present this session), not `.venv-msk`; never imports OpenSim.

## 3. Results in detail

### 3.1 Q1 — no genuine two-hump structure at the ankle

The mechanical 30%GC-split hump-finder (reused verbatim for cross-joint table comparability with knee/hip)
reports an SO "hump1" of 307.5 %BW @ 24%GC — **but this is explicitly NOT a genuine hump.** Sampling the raw
curve at 5%GC steps shows the ankle contact force rising essentially continuously from 61.2 %BW at heel-strike
(0%GC) through 146.4 (10%GC), 266.6 (20%GC), 309.9 (30%GC), 386.1 (40%GC) up to the single peak of 484.3 %BW
at 50%GC, then falling sharply through 201.4 (60%GC) to a swing-phase floor of ~45-51 %BW (65-95%GC) before
rising back into the next cycle. 84% of the frame-to-frame steps between 0-50%GC are non-decreasing (the
rest is small negative jitter, the same per-frame Static-Optimization noise already documented in the knee/hip
cert). The pre-registered genuine-hump test confirms this quantitatively: prominence of the best early-stance
candidate is 0.04% (SO, smoothed) to 8.9% (CMC, raw — the closest call) of the curve's own overall peak, in
every one of the four smoothed/raw × SO/CMC combinations tested — all far below the pre-registered 15% bar.
**The ankle contact force has a single dominant late-stance peak, not a double hump** — qualitatively
different from knee/hip's genuine, large-margin double-hump structure, and directly consistent with the
cadaveric literature's own description (Stauffer 1977's own words: force "rose... during the **latter** part
of stance phase," with no mention of an earlier peak).

### 3.2 Q2 — the ankle peak instant lands with the plantarflexors, not with knee/hip

| peak | t (s) | %GC |
|---|---:|---:|
| Knee/hip (SO/CMC, established) | 0.51 | 44.80 |
| **Ankle contact force (this doc)** | **0.58** | **50.08** |
| Achilles-equivalent tendon force (established) | 0.59 | 50.83 |

The ankle contact-force peak sits **0.75%GC (10ms) from the plantarflexors' own natural peak** and
**5.28%GC (70ms) from the knee/hip peak** — a clean, well-separated confirmation of the phase offset the
plantarflexor thread already found and this task pre-registered as the corroborating signal to check. This is
the single sharpest, most decisive piece of evidence in this analysis, because it requires no anchor at all
(purely a comparison of model-internal, exactly-measured instants) — unlike the magnitude comparisons (§3.3),
it carries no cadaveric-anchor-tier uncertainty.

### 3.3 Q3 — the ankle already matches its anchor at push-off; no material excess anywhere

The push-off-window (30-62%GC) peak value (484.3 %BW) sits at **1.016x** the blended 3-study anchor mean
(476.7 %BW) and **1.088x** the stricter, talocrural-joint-only 2-study anchor mean (445.0 %BW, Stauffer +
Procter & Paul only, excluding Giddings' anatomically-adjacent-but-distinct talocalcaneal number) — both
comfortably inside the pre-registered ±20% materiality band, and this result is **robust to which anchor
construction is used.** A full-resolution scan (every one of the 101 %GC grid points, not just the two
coarse hump windows) finds **zero** %GC locations where the model/anchor ratio exceeds 1.20x, against EITHER
anchor — the strongest, most decisive form of "no over-prediction" available from this data: not merely
"the peak matches," but "nowhere in the entire gait cycle does the ankle contact force materially exceed its
(weak-tier) cadaveric anchor." The weight-acceptance-window value (309.9 %BW @ 30%GC) is well BELOW the anchor
(ratio 0.650-0.696x) — but this comparison is **explicitly phase-mismatched**: no literature source
separately quantifies an early-stance ankle contact-force anchor (all three studies report only a single,
late-stance-understood peak number), so this ratio shows the model's early-stance value is smaller in scale
than the late-stance anchor (expected under any reasonable model), not that the model "matches" or
"under-predicts" a genuine weight-acceptance ground truth that simply does not exist in the published
record.

### 3.4 Waveform correlation — real, but the weakest discriminator here

The ankle contact-force curve correlates strongly with the Achilles-equivalent tendon-force curve (r=0.953
full-cycle, r=0.930 restricted to the 30-65%GC push-off window) — consistent with the mechanistic expectation
that the ankle joint reaction is substantially determined by the plantarflexor tendon tension crossing it.
**Read carefully, not oversold:** the ankle-vs-KNEE contact-force correlation is nearly as high (r=0.934),
because knee, ankle, and plantarflexor-force curves all substantially share the same gait-cycle-driven
loading envelope (rise through stance, fall in swing) — a generic feature any two joint-level quantities in
the same trial would likely share. The correlation DELTA (0.953 vs 0.934, gap=0.019) is small and does not by
itself sharply discriminate "ankle tracks plantarflexors specifically" from "ankle tracks the generic gait
envelope that everything tracks." This piece is reported as consistent, supporting context — the decisive
evidence for the corroboration claim is §3.2's timing offset and §3.1/§3.3's shape/magnitude findings, not
this correlation number in isolation.

### 3.5 A disclosed divergence: CMC ankle is LOWER than SO, opposite the established knee/hip pattern

CMC's ankle peak (405.54 %BW @ t=0.57s) is **16.3% LOWER** than SO's (484.48 %BW), ratio 0.837 — the
**opposite direction** from the already-established knee/hip pattern where CMC reads 8-11% HIGHER than SO
everywhere (`docs/MECHANISM_CONTACT_WAVEFORM.md`). This is new territory (no prior CMC-ankle number exists
anywhere in this repo to self-consistency-check against), reported plainly rather than forced to fit the
knee/hip pattern — a genuine, disclosed cross-solver divergence at this joint, not explained further here (no
per-muscle CMC-ankle decomposition was run). **A separate, independent replication in the SAME direction as
an already-disclosed finding:** CMC's swing-phase floor at the ankle (68.7-105.3 %BW, 65-95%GC) reads
consistently HIGHER than SO's (44.6-51.0 %BW) — reproducing, at a THIRD joint, the already-disclosed
CMC-forward-dynamics-tracking-controller residual-loading-during-swing artifact
(`docs/MECHANISM_CONTACT_WAVEFORM.md` §3.2's "during SWING... CMC-knee reads 79.3 %BW vs OrthoLoad's 26.4%BW
... SO-knee reads only 34.3%BW"; `docs/MECHANISM_RRA_TASK_GAINS.md`; `docs/MECHANISM_ANKLE_RESERVE_FIX.md`) —
strengthening that finding's generality even as the PEAK-level CMC/SO relationship itself flips sign at this
joint.

## 4. Honest gaps

1. **Weak anchor tier (pre-registered, not a new caveat):** no in-vivo instrumented-ankle program has ever
   existed (re-verified absent again this session, consistent with `docs/MECHANISM_ANKLE_FORCE.md` §0 and
   `docs/MECHANISM_MSK_BUILD_PLAN.md`). The anchor here is cadaveric/quasi-static/finite-element literature
   (1977-2000), not OrthoLoad-grade in-vivo instrumented data — a plausibility check, not a validation with
   the knee/hip cert's epistemic weight.
2. **The anchor is a single peak number per study, not a phase-resolved curve.** Unlike OrthoLoad's
   `averDTW`-aligned, hump-specific knee/hip ensemble, none of the three ankle studies separately quantifies a
   weight-acceptance-phase force — so the weight-acceptance ratio (§3.3) is a scale comparison against a
   phase-mismatched anchor, disclosed at every use, not a like-for-like in-vivo check.
3. **Cross-study spread is not a population estimate.** The 3-study SD (77.67 %BW) or 2-study talocrural-only
   SD reflects METHOD heterogeneity (quasi-static cine-film vs 3-D force-equilibrium model vs finite-element
   model) as much as true inter-subject variance — deliberately NOT reported as a Z-score the way the
   OrthoLoad-anchored knee/hip cert's Z-scores are, to avoid implying an epistemically equivalent population
   estimate from n=2-3 heterogeneous studies.
4. **Giddings' 540 %BW figure is the talocalcaneal (subtalar) joint, not talocrural** — disclosed and handled
   by reporting BOTH the blended (with Giddings) and talocrural-only (without) anchor constructions throughout;
   the verdict is robust to either choice (ratio 1.016x vs 1.088x, both non-material).
5. **The 15%-of-own-peak genuine-hump threshold and the 20% materiality bar are both chosen conventions**, not
   physical constants — reported alongside the raw prominence numbers (0.02-8.9%, all well clear of 15%) and
   raw ratios (1.016-1.088x, both well clear of 1.20x) specifically so the reader can see the findings are not
   knife-edge calls sensitive to the exact threshold chosen.
6. **n=1 subject/trial (subject2/walking1), single (right) leg** — same scope caveat inherited from every
   cert in this family; no claim of generality across subjects/speeds.
7. **The CMC-ankle-lower-than-SO divergence (§3.5) is reported, not explained** — no per-muscle CMC-ankle
   decomposition was run this session (out of scope for this waveform-focused task; a natural, cheap follow-up
   flagged, not pursued here per LEAN discipline since it does not bear on this doc's own falsifier).
8. **The %GC mapping (T_stride=1.327s contralateral-symmetry convention) is inherited, re-verified (both
   out-of-sample alignment checks re-run and PASS), but not independently re-derived a third way** — same
   caveat already carried by every doc in this family that uses this mapping.

## 5. Verification (machine-checked, not narrated)

- Self-consistency: SO ankle peak (484.480 %BW @ t=0.580s) vs `docs/MECHANISM_ANKLE_FORCE.md`'s published
  484.48 @ 0.58 — 0.0000% relative difference, exact instant match (asserted in-script). Achilles-equivalent
  peak (320.265 %BW @ t=0.590s) vs `docs/MECHANISM_PUSHOFF_PLANTARFLEXOR.md`'s published 320.2649 @ 0.59 —
  0.0000% relative difference, exact instant match (asserted in-script).
- Phase-alignment: both out-of-sample cross-checks (R toe-off 66.9%GC in [55,68] band; L toe-off 17.0%GC in
  [8,20] band) re-run and PASS (`phase_alignment_verification.all_pass = true`).
- Grid alignment: `force.sto` and `JointReaction.sto` time vectors confirmed bit-identical
  (`np.allclose == True`) before any cross-file resampling was trusted.
- Talocrural-only anchor filter bug (n=1/mean=390 instead of n=2/mean=445) caught by cross-checking against
  `docs/MECHANISM_ANKLE_FORCE.md`'s own prose, fixed via explicit PMID allowlist, re-verified
  (`talocrural_only_n=2`, asserted in-script, would raise otherwise).
- Genuine-hump test run in 4 configurations (SO/CMC × smoothed/raw) — all 4 return `genuine_early_hump=False`,
  all well below the pre-registered 15% bar (max observed: 8.9%, CMC raw).
- Full-resolution material-excess scan run against BOTH anchor constructions (blended, talocrural-only) —
  both return an empty band list.
- All numbers in this document are pulled directly from
  `data/msk_smoketest/subject2_walking1/ankle_waveform/ankle_waveform_results.json`, not transcribed from
  console prose.

## 6. Files

- `scripts/msk/ankle_waveform_analysis.py` (new) — the full, re-runnable pipeline. Run with
  `/usr/bin/python3 scripts/msk/ankle_waveform_analysis.py` (needs numpy+openpyxl transitively via its import
  of `contact_waveform_analysis`; deliberately NOT `.venv-msk`; never imports OpenSim). Imports
  `scripts/msk/contact_waveform_analysis.py` directly rather than re-implementing its parser/%GC-mapping/
  hump-finder — exit 0, all gates PASS.
- `data/msk_smoketest/subject2_walking1/ankle_waveform/ankle_waveform_results.json` (new) — every number in
  this document: the %GC mapping + alignment verification, SO/CMC ankle contact-force curves (full 101-point
  grid) + self-consistency gates, the Achilles-equivalent tendon-force curve + self-consistency gate, the
  genuine-hump test (4 configurations), the peak-timing offsets, the waveform-correlation coefficients, and
  the full anchor comparison (both anchor constructions, both windows, the full-resolution material-excess
  band scan, smoothing-sensitivity check).
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/static_optimization/{so,jr}/`,
  `data/msk_smoketest/subject2_walking1/cmc_second_solve/jr_cmc/` (model force/reaction time series),
  `data/msk_smoketest/subject2_walking1/ankle_force_validation/ankle_force_validation_results.json` (anchor
  literature + established peak numbers, read not retyped),
  `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/ForceData/walking1_forces.mot`
  (GRF, for heel-strike detection, via the imported mapping function).
- **Not touched**: `scripts/msk/contact_waveform_analysis.py` (imported read-only),
  `scripts/msk/validate_ankle_force.py`, `scripts/msk/pushoff_plantarflexor_achilles.py`,
  `scripts/msk/contact_muscle_decomp.py`, and any concurrent, independent instance's in-progress work.

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task). All new files
are untracked, for the coordinator to commit.
