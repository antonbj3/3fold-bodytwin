# MECHANISM METABOLIC COST vs IN-VIVO CALORIMETRY — does the dampened common-mode elasticity mean the metabolic layer is MORE trustworthy than the force certs? (2026-07-21)

**Question** (from `docs/MECHANISM_COMMON_MODE.md`): that probe measured `metabolic_cost`'s Umberger2010
net-power **elasticity to a uniform +20% SO-activation-bias perturbation = 0.55** — dampened, vs e.g.
`muscle_fatigue`'s 15.06× amplification or `tendon_elastic_energy`'s 1.73× amplification. A naive reading
says: "the metabolic layer barely moves when the shared upstream SO solve is wrong, so it deserves MORE
trust than the joint-contact-force certs (which over-predict OrthoLoad in-vivo knee/hip force by
~1.41-1.66×)." **This doc tests that reading directly** by anchoring the twin's metabolic cost of
transport (COT) against measured indirect-calorimetry human-walking energetics, at the matched
definition and (as closely as achievable) matched speed, and comparing the resulting ratio-to-anchor
against the force layer's own already-published ratio-to-anchor, same subject, same trial, same two
solve methods (SO and CMC).

**No re-solve.** Every twin-side number below is reloaded from two already-published, machine-generated
JSONs (`metabolic_cost_results.json`, `cmc_second_solve_results.json`) by
`scripts/msk/metabolic_calorimetry_anchor.py`, which recomputes every ratio from the RAW numbers in
those files (never trusts a pre-stored ratio blindly) and cross-checks the two JSONs agree with each
other before using either. The only new work this session is (a) literature verification and (b) the
cross-layer comparison itself.

## Pre-registration

**C** ("metabolic deserves a HIGHER trust tier than the force certs, consistent with its dampened
elasticity"): the metabolic-layer's ratio-to-external-calorimetry-anchor is materially **tighter**
(lower) than the force-layer's ratio-to-OrthoLoad, on both the mean and the max of a 4-vs-4 comparison
(2 solve methods × 2 sub-metrics per layer), by a non-trivial margin (≥15% relative, fixed before running
the comparison script — see the honesty note below on what "pre-registered" can mean here).

**¬C** ("over-prediction is a broad, systemic feature of the muscle-driven layers generally, not
specific to the free-body joint-contact estimate"): the metabolic ratio is NOT tighter — i.e. C's
margin test fails.

**Falsifier, stated by the task itself before this doc was assembled**: if twin metabolic COT is ~2×
measured (in the same range the force certs already are, ~1.5-2×), that supports ¬C/"broader than
force"; if metabolic matches calorimetry distinctly better than force matches OrthoLoad, that supports
C/"different trust tiers."

**Honesty note on pre-registration**: the component numbers (SO/CMC COT, OrthoLoad force ratios) were
already read into context during this session's orientation pass, before the 15%-margin rule was written
into the script — true blinding of those specific numbers was not possible (they are already-published,
committed facts, not a fresh experiment). What IS true, and reported precisely rather than glossed over:
the discriminating rule was fixed at a generic, round, non-cherry-picked margin, and the actual result
(below) is **not a knife-edge call the exact margin choice could flip** — the sign of the gap between
metabolic and force ratios is the wrong direction for C at ANY positive margin, so this is a robust,
direction-level finding, not a threshold-sensitive one.

## Headline result

| | SO | CMC |
|---|---:|---:|
| Umberger2010 COT_net (J/kg/m) | 7.6834 | 8.3198 |
| Bhargava2004 COT_net (J/kg/m) | 6.1977 | 6.1239 |
| twin trial speed (m/s) | 1.0647 | 1.0605 |

| ratio to external anchor | knee force / OrthoLoad | hip force / OrthoLoad | Umberger COT / Koelewijn | Bhargava COT / Koelewijn |
|---|---:|---:|---:|---:|
| SO | 1.515 | 1.412 | 1.945 | 1.569 |
| CMC | 1.655 | 1.564 | 2.106 | 1.550 |

- **Force layer** (4 values: knee×{SO,CMC}, hip×{SO,CMC}, vs OrthoLoad in-vivo implant telemetry,
  258.22 %BW knee / 273.93 %BW hip): **mean 1.536×, max 1.655×**.
- **Metabolic layer** (4 values: Umberger×{SO,CMC}, Bhargava×{SO,CMC}, vs Koelewijn 2019 in-vivo
  indirect-calorimetry net COT, 3.95 J/kg/m interpolated level-1.3-m/s anchor): **mean 1.793×, max
  2.106×**.
- **Metabolic is 16.7% HIGHER than force on the mean, 27.3% HIGHER on the max.** C is **NOT
  SUPPORTED — falsified**. ¬C is supported: the over-prediction is broader than the free-body
  joint-contact-force estimate specifically; it recurs, at similar-to-larger magnitude, in a
  mechanistically unrelated (muscle-heat-rate-equation-based, not Newton's-law-based) downstream layer.
- **Not driven uniformly across both metabolic models**: sorting all 8 ratios, the two Umberger values
  (1.945, 2.106) sit clearly ABOVE the entire force range [1.412, 1.655]; the two Bhargava values
  (1.569, 1.550) sit WITHIN/at the top of the force range, not distinguishably worse. The
  "metabolic ≥ force" finding is real on the aggregate (mean, max) but is **concentrated in the
  Umberger2010 probe specifically** — consistent with, and reinforcing, `docs/MECHANISM_METABOLIC_COST.md`'s
  own already-disclosed "two-model spread 21.4%" finding (Umberger consistently reads higher than
  Bhargava in this twin). Reported precisely, not smoothed into a single number that would hide this.

**What this means for the elasticity finding**: elasticity (0.55) and this ratio-to-anchor are
**different axes that are shown here to disagree**. Elasticity measures how much a layer's OWN headline
number moves in response to ONE specific, controlled upstream perturbation (a uniform +20% SO-activation
bias) — a propagation-sensitivity property. This doc's ratio-to-anchor measures TOTAL absolute error
against ground truth, accumulated from every source at once (solve-method choice, muscle-model-formulation
choice, muscle-mass estimation, tendon-compliance assumption, the basal constant, and whatever the
+20%-bias probe was itself modeling). A layer can be simultaneously **less sensitive to one named
perturbation** and **no more (here, less) accurate in total** — these are not the same claim, and this
project's own trust-ledger tiering should not conflate them. **The dampened elasticity does NOT imply the
metabolic layer is better-trusted; measured here, it is not.**

## The forced adversaries (both checked, neither explains the gap away)

**1. Definitional mismatch (net vs. gross vs. net-above-standing)** — the task's own named risk, since
these conventions can differ ~2×. Checked explicitly:
- Twin's **"net"** = raw `Umberger2010MuscleMetabolicsProbe`/`Bhargava2004MuscleMetabolicsProbe` output
  (muscle heat-rate + mechanical-work-rate), **zero** basal/resting term added — i.e. cost attributable
  to the 80 modeled muscles' own activity only, nothing else.
- Koelewijn's **"net"** = directly-measured (indirect calorimetry) whole-body metabolic rate during
  walking, MINUS directly-measured quiet-STANDING metabolic rate (verbatim, their Methods: "The resting
  trial was subtracted from each walking trial," where the resting trial was subjects who "stood on the
  treadmill for three minutes"). This is precisely the convention the task brief calls **"net-above-
  standing"** — confirmed by direct quote, not assumed.
- These two "nets" are conceptually matched (both isolate the increment attributable to the exercise of
  walking, over a quiet/basal reference) though mechanistically different (twin: zero baseline by
  omission of non-muscle metabolism entirely; Koelewijn: a real, measured, non-zero standing baseline
  subtracted). **Forced check of the failure mode**: if gross (net + 1.2 W/kg basal add-on) were
  mistakenly compared against Koelewijn's net anchor instead, the ratio would read 2.230 (SO) / 2.393
  (CMC) instead of the correct 1.945 / 2.106 — a **+13.6% to +14.7% spurious inflation**, real but
  **secondary**, not the ~2× a genuine gross/net mismatch can cause in other contexts (e.g. very slow
  walking, where the fixed basal term dominates a near-zero net signal). Net-vs-net (used throughout) is
  confirmed correct and materially different from the mistaken comparison — the observed ~1.6-2.1× gap is
  not a definitional artifact.
- **Second-order, argued-not-corrected point**: the twin's "net" structurally omits 31/111 modeled DOF
  actuators (ideal torque actuators for trunk/arms, pelvis residuals, joint reserves — a disclosed
  under-estimate source from `docs/MECHANISM_METABOLIC_COST.md`) that a real, intact human's whole-body
  calorimetry measurement necessarily includes. If this were quantified and added, it could only push the
  twin's number UP, not down — meaning the true like-for-like gap, if anything, is **understated**, not
  manufactured, by the current comparison.

**2. Speed mismatch (1.06 m/s twin vs. 1.3 m/s Koelewijn anchor)** — checked via Ralston HJ (1958), the
paper that established the existence of an energy-optimal walking speed. Ralston's own primary abstract
is not extractable (see citation table — pre-dates PubMed's abstract-indexing era), but a secondary,
encyclopedia-tier source (Wikipedia "Preferred walking speed," fetched live, attributing the numbers to
Ralston 1958 by name) gives: **gross-COT-minimizing speed ≈ 1.23 m/s; net-COT-minimizing speed ≈ 1.05
m/s.** The twin's own trial speed (1.060-1.065 m/s) sits **~0.01-0.015 m/s from Ralston's own reported
net-COT-minimizing speed** — i.e., essentially on top of the physiologically CHEAPEST point on the
net-COT-vs-speed curve, not in a slower/more-expensive regime. **This adversary does not explain the
over-prediction — if anything it argues the opposite direction**: had the twin walked slower than 1.05
m/s, its TRUE physiological net COT would be expected to rise (U-shaped curve), which would have made the
model-vs-measured gap partly attributable to a genuine speed effect; instead, the twin's trial speed
picks almost exactly the speed at which a real human's net COT is lowest, so the observed over-prediction
cannot be explained away as "the twin was just walking at an inherently more expensive speed than the
anchor." (Argued from the verified curve-SHAPE finding; not independently re-derived from Koelewijn's own
single-speed dataset, which cannot show a speed slope on its own — a disclosed scope limit on this
specific check.)

## Citations — verified LIVE this session (NCBI eutils + Europe PMC, not recalled)

| citation | PMID | verified | role here |
|---|---|---|---|
| Ralston HJ (1958) "Energy-speed relation and optimal speed during level walking." *Int Z Angew Physiol Einschl Arbeitsphysiol* 17(4):277-83 | **13610523** | title/journal/year/pages via NCBI esummary; **no abstract exists in the MEDLINE record** (pre-dates abstract indexing) | Establishes the optimal-speed/U-shaped-curve concept by name; specific numbers (1.23/1.05 m/s) sourced from a secondary, encyclopedia-tier fetch (Wikipedia), flagged at that lower confidence, not primary-verified |
| Browning RC, Kram R (2005) "Energetic cost and preferred speed of walking in obese vs. normal weight women." *Obesity Research* 13(5):891-9 | **15919843** | full abstract verified via Europe PMC REST API | Speed-triangulation only: preferred speed 1.47 m/s (normal-weight), 1.40 m/s (obese). Abstract's own metric is explicitly **GROSS** cost — NOT used as a second net-COT quantitative anchor (would risk the very definitional mismatch this doc guards against) |
| Miller RH (2014) "A comparison of muscle energy models for simulating human walking in three dimensions." *J Biomech* 47(6):1373-81 | **24581797** | full abstract verified via Europe PMC REST API | See correction below — used as an independent inter-model-spread bracket, not a model-vs-calorimetry finding |
| Koelewijn AD, Heinrich D, van den Bogert AJ (2019) "Metabolic cost calculations of gait using musculoskeletal energy models, a comparison study." *PLoS ONE* 14(9):e0222037 | **31532796** (PMC6750598) | PMC full text; this is the 3rd independent fetch overall (2 in the prior session that built `MECHANISM_METABOLIC_COST.md`, 1 more this session targeted specifically at Introduction/Discussion/Fig 5) | **Primary quantitative anchor** — real indirect-calorimetry measurement in living human subjects (n=12), net = standing-subtracted |

**Correction to the task brief's framing (same discipline this repo's own docs already apply — e.g.
`metabolic_cost.py`'s Umberger-PMID catch)**: the task brief characterized Miller 2014 as showing
"Umberger2010/Bhargava2004 muscle-energetics models systematically over-predict whole-body COT vs
calorimetry." **Miller 2014's verified abstract does not say this.** It reports a **model-vs-model**
finding: 5 different Hill-based muscle-energy models, applied to the SAME simulated walking motion,
predict costs "vary[ing] over a three-fold range (2.45-7.15 J/m/kg), with the distinction arising from
whether or not eccentric work was subtracted from the net heat rate" — no calorimetry comparison appears
in the abstract at all. This is still directly useful (see below), just for a different, more precise
reason than the task brief stated.

**Miller 2014 bracket check**: twin Umberger values [7.683 (SO), 8.320 (CMC)] sit **at/above the top** of
Miller's own published 5-model range (2.45-7.15); twin Bhargava values [6.198 (SO), 6.124 (CMC)] sit
**within** Miller's range. Cross-study caveat: Miller's range comes from a different simulated walking
motion/subject/protocol (optimal-control-generated gait, not subject2's mocap-driven SO/CMC), so this is
a magnitude bracket, not a strict replication — but it corroborates, independently, that muscle-energy-
**model-formulation choice alone** is a first-order, ~3× source of variance, separate from (additive to)
the SO-vs-CMC solve-method spread this repo already measured (8.3% Umberger, 1.2% Bhargava,
`docs/MECHANISM_CMC_SECOND_SOLVE.md`), and helps explain why the twin's own Umberger number can land at
the high end of what published models even produce.

**Margaria** — the task brief's "Margaria" citation could not be resolved to a specific PMID. A live
PubMed author search (`Margaria R[Author] AND walking`) returned 5 hits, verified via esummary: none is
the classic 1938 walking-economy-vs-speed curve (that paper, in Italian, pre-dates MEDLINE indexing by a
decade and has no PMID). The 5 hits found are Cavagna & Margaria's **mechanical** (not metabolic) walking-
mechanics papers (1963, 1965, 1966), a subgravity-biomechanics piece (1973), and a clinical foot-lesion
energetics paper (1979) — none gives the classic healthy-adult metabolic-COT-vs-speed number. Separately,
Koelewijn 2019's own reference list includes a "MARG68" model variant (Margaria 1968, a mechanical-
efficiency-based model: "25% efficient shortening, 120% efficient lengthening") — a **different, later**
Margaria reference than the classic 1938 empirical curve the task brief likely intended. Disclosed as an
honest gap, not silently substituted or fabricated.

## Confidence tier

Per **this repo's own established ledger convention** (`docs/MECHANISM_TRUST_LEDGER.md`, tier legend),
**not** the task brief's generic phrasing: the tier legend explicitly reserves **in-vivo-anchored** for
"OrthoLoad implant telemetry... the source for every row in this tier" — a term-of-art for the specific
joint-force-vs-force comparison, not a general "measured in a living human" label. Koelewijn's indirect
calorimetry IS a real in-vivo measurement in living humans, but it is a different modality (whole-body
respiratory gas exchange) compared against a different-methodology model output (muscle heat-rate
equations, not a free-body-diagram joint force) — the trust ledger's own existing row for
`METABOLIC_COST.md` already classifies this identical Koelewijn comparison as **cadaveric-or-published-
plausibility**, and this doc follows that same, already-established convention rather than relitigating
it. So: **cadaveric-or-published-plausibility** (a genuinely external, decorrelated, published human
dataset) — **muscle-energetics-model-limited** (Miller 2014 independently confirms model-formulation
choice alone spans a ~3× range; this twin's own solve-method choice (SO vs CMC) adds a further 1.2-8.3%;
neither of these is a measurement-pipeline defect, both are open, named modeling-choice uncertainties).

## What this does and does not resolve

**Resolved**: the specific hypothesis that dampened common-mode elasticity (0.55) implies the metabolic
layer deserves MORE trust than the force certs is **not supported** — measured here, the metabolic
layer's absolute over-prediction against its own external anchor is comparable to, and by mean/max
slightly larger than, the force layer's. Both forced adversaries named by the task (definitional
mismatch, speed mismatch) were checked explicitly and neither explains the gap away; if anything both
argue the true gap is at least as large as measured, not smaller.

**Not resolved**: this remains a single-subject, single-trial (subject2/`walking1`) comparison, inheriting
every upstream common-mode risk `docs/MECHANISM_COMMON_MODE.md` and `docs/MECHANISM_CMC_SECOND_SOLVE.md`
already disclosed (both SO and CMC share the same scaled model, same subject, same single trial; CMC
carries its own uncleared Hicks-moment-residual defect). This doc does not decompose HOW MUCH of the
metabolic-specific extra gap (beyond what force already shows) traces to muscle-mass over-estimate vs.
tendon-compliance choice vs. Miller's model-formulation axis specifically — `docs/MECHANISM_METABOLIC_COST.md`'s
own sensitivity sweep already shows the combined tendon+mass correction alone can move the ratio from
1.95 down to 0.93, so the "true" gap is highly sensitive to disclosed, still-open modeling choices, not a
single fixed number.

## Honest gaps

1. Ralston 1958's own primary quantitative values are not independently verified this session (no
   abstract in the MEDLINE record); the specific numbers used (1.23/1.05 m/s) come from a secondary,
   encyclopedia-tier source, one citation-hop removed from the primary literature.
2. Margaria's classic 1938 walking-economy value has no PMID and was not obtained from any source this
   session — a disclosed, not fabricated, gap.
3. Browning & Kram 2005 is used only for its speed datapoints, not as a second quantitative net-COT
   anchor — its own actual gross-COT J/kg/m number was not extracted (would require the paywalled full
   text/tables, not attempted this session).
4. Koelewijn's own level-walking (0% grade) net COT is still not extractable as a single verbatim number
   from live-fetchable text (3 independent attempts across 2 sessions, this session's attempt targeted
   specifically at Fig 5) — the 3.95 J/kg/m figure remains OUR OWN linear interpolation of the verified
   downhill/uphill values, not Koelewijn's own reported number. Added this session: the true grade-vs-COT
   relationship is known (general gait-energetics literature) to be non-linear/convex, so the linear
   interpolation carries an unquantified-direction bias — disclosed, not resolved.
5. The 4-vs-4 (metabolic-vs-force) ratio comparison is a comparison of small, internally-correlated
   samples (both "Umberger vs Bhargava" and "knee vs hip" share substantial common-mode through the same
   underlying SO/CMC solve) — not 4 independent replicates in either arm. The qualitative direction
   (metabolic ≥ force) is robust to this (driven mainly by Umberger specifically, per the sorted-ratios
   finding above), but the exact "+16.7%/+27.3%" figures should be read as a modest-sample comparison, not
   a large-N statistical result.
6. Single subject (subject2), single trial (`walking1`), one speed (~1.06 m/s) — no claim of generality
   across subjects, trials, or speeds. Same scope boundary as every other cert in this family.
7. This doc does not re-run or re-verify the force-layer's own OrthoLoad numbers (258.22/273.93 %BW) or
   sign-bug-correction history — those are reloaded read-only from `cmc_second_solve_results.json`,
   itself already cross-checked against `docs/MECHANISM_STATIC_OPT.md`'s correction banner in this
   session's orientation pass (391.10/391.11 %BW match exactly, confirming the corrected, not the stale
   pre-fix 233.20 %BW, pipeline is what feeds this comparison).

## Files

- `scripts/msk/metabolic_calorimetry_anchor.py` — the full, re-runnable comparison script (no OpenSim
  dependency, no re-solve; reloads `metabolic_cost_results.json` + `cmc_second_solve_results.json`,
  recomputes every ratio from raw numbers, runs the pre-registered discriminating test + both forced
  adversaries + the Miller 2014 bracket check).
- `data/msk_smoketest/subject2_walking1/metabolic_calorimetry_anchor/metabolic_calorimetry_anchor_results.json`
  — every number in this doc, machine-written, including the full citation table with PMIDs.
- Inputs read in place, never modified: `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`,
  `data/msk_smoketest/subject2_walking1/cmc_second_solve/cmc_second_solve_results.json`.
- Reused context, not modified: `docs/MECHANISM_COMMON_MODE.md` (the elasticity=0.55 finding this doc
  tests), `docs/MECHANISM_CMC_SECOND_SOLVE.md` (the force-vs-OrthoLoad ratios), `docs/MECHANISM_STATIC_OPT.md`
  (confirms the corrected, non-buggy force numbers), `docs/MECHANISM_TRUST_LEDGER.md` (the tier vocabulary).

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
