# MECHANISM GRIP STRENGTH — the hand/finger-flexor force chain, function ↔ dysfunction (2026-07-22)

Fills the named MSK gap: **no grip-force doc/instrument exists anywhere in this twin** (confirmed
live this session — `docs/MECHANISM_HAND_FIVE_DIGIT.md`L42 states plainly "no live grip-force
instrument exists in this pipeline"; `docs/MECHANISM_ELBOW_WRIST_FORCE.md`L203 confirms "no
fingertip/grip marker exists in this model's `MarkerSet`"). This doc builds the quantitative
force-chain model (PCSA/Fmax → tendon → carpal-tunnel/pulley redirect → joint torque → fingertip
force), states it against the twin's OWN machine-measured model numbers, and forces two adversaries
to their strongest form against real, live-verified external data. **Both forced adversaries FALL,
but in opposite directions**: (1) the classic **inverted-U grip-force-vs-aperture** shape is
confirmed 4/4 across independent large-N studies (peak ≈5–5.5cm, tighter than this task's own
5–7cm prior); (2) a **geometric reconstruction of maximum grip magnitude from the twin's own
Fmax/moment-arm numbers falls well short of real population norms** — under every one of 8 tested
combinations of budget/transmission-ratio assumptions, predicted max grip (45–134N) lands BELOW
even the EWGSOP2 clinical sarcopenia cutoff (157–265N), let alone the healthy population median
(304–500N). This is a genuine, non-redundant finding: the hand model is **under-strong**, the
opposite direction from the knee/hip finding in `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`
(over-strong) — different lineage, different joint, both plausible in a composite twin.

## 0. Headline scorecard

| Claim | Pre-registered threshold | Forced adversary | Result | Anchor |
|---|---|---|---|---|
| **Shape**: grip force vs. grip-span/aperture is an inverted-U, interior peak (not flat, not monotonic) | Peak must be interior to the tested span range in ≥2 independent studies | "flat" (trivial) then "monotonic" (strongest fair form: a pure lever-arm argument predicts monotonic if r(θ) is monotonic in θ) | **Adversary FALLS 4/4** — Firrell&Crain 1996 (n=288, 89% peak at interior setting II), Ruiz-Ruiz 2002 (n=70, optimal span identified both sexes), Lee/Kong/Lowe 2009 (n=46, peak at 50–55mm of 45–65mm tested), Blackwell 1999 (n=18, "middle > small OR large") | Peak ≈5.0–5.5cm (Ruiz-Ruiz men fixed 5.5cm; Lee et al. 50–55mm); tighter than the 5–7cm prior |
| **Peak location numeric** | within task's stated 5–7cm | — | **PASS**, refined to 5.0–5.5cm (lower/tighter than the stated envelope) | Ruiz-Ruiz 2002 PMID 12239682; Lee et al. 2009 PMID 19424925 |
| **Magnitude**: twin's own Fmax budget reaches realistic max-grip norms | predicted grip force (any fair geometric assumption) ≥ EWGSOP2 low-strength cutoff (27kg/16kg) | steelmanned: generous 455N all-closing-muscle budget × most-favorable (DIP-joint, 3.40x) transmission ratio, cross-checked against Schuind 1992's real in-vivo 3.43x | **Adversary FALLS 8/8** (2 budgets × 4 geometric ratio-bases): predicted 44.7–133.9N, all below the female EWGSOP2 cutoff (156.9N), all below the male cutoff (264.8N) | Dodds 2014 PMID 25474696 (norms); EWGSOP2 PMID 30312372 (cutoffs); Schuind 1992 PMID 1564277 (in-vivo tendon force, geometric cross-check) |
| **Dysfunction**: grip strength is a validated mortality/frailty biomarker | HR CI excludes 1.0, independent of confounders, large N | N/A (literature anchor, not re-derived) | **Confirmed**, HR 1.16 (1.13–1.20) per 5kg decrease, all-cause mortality, n=139,691, 21 countries | Leong et al. 2015 Lancet PMID 25982160 |

All PMIDs in this doc were live-verified this session via NCBI eutils (`esearch`/`esummary`/`efetch`
against `eutils.ncbi.nlm.nih.gov`) and, for EWGSOP2's exact cutoff table, a direct fetch of the PMC
full text. Two citation corrections were caught and fixed in the process (not silently propagated):
the correct Leong PURE-mortality paper is PMID 25982160 (a companion paper, PMID 27104109, is a
DIFFERENT PURE norms-only paper, cited separately below); PMID 31081853 for EWGSOP2 is literally an
**erratum notice** for the real paper (PMID 30312372) — the primary citation used throughout is
30312372.

---

## 1. Mechanism: the force chain, in order

Power/precision grip force originates in the **extrinsic finger flexors** — flexor digitorum
profundus (FDP) and superficialis (FDS), forearm-bellied muscles — and is modulated by the
**intrinsics** (lumbricals, interossei, thenar group). The chain, in strict anatomical order:

1. **Muscle belly** (forearm): PCSA × specific tension × activation × force-length ×
   force-velocity multipliers → active fiber force, reduced by `cos(pennation angle)` onto the
   tendon axis.
2. **Tendon** (near-inextensible but with disclosed series compliance — see
   `docs/MECHANISM_TENDON_ELASTIC.md`, `MECHANISM_TENDON_STRAIN.md`, `MECHANISM_TENDON_SLACK_SENSITIVITY.md`,
   not re-litigated here): carries the muscle's tension distally, largely unchanged in magnitude.
3. **Carpal tunnel**: a rigid redirect under the transverse carpal ligament — changes tendon
   *direction*, not magnitude (friction losses disclosed as a simplification, not modeled).
4. **Annular/cruciate pulley system (A1–A5 equivalent)**: sets the **moment arm** at each of
   MCP/PIP/DIP geometrically (twin's own pulleyed model, `docs/MECHANISM_HAND_PULLEYS.md`, confirmed
   this session — see §2). Joint torque τ = F_tendon × r(θ) at each joint the tendon crosses.
5. **Multi-joint chain → fingertip force**: converting 3 joint torques into a single fingertip
   force requires the chain's geometric Jacobian — this is where a **short moment arm vs. long
   phalanx lever-arm** produces a real mechanical disadvantage (same principle the twin's own
   `docs/MECHANISM_ELBOW_WRIST_FORCE.md` already established at the elbow/wrist: "a short
   muscle-tendon moment arm... versus a much longer lever arm... the classic third-class-lever
   mechanical disadvantage" — §4 below re-derives this at the finger joints specifically, using the
   twin's own numbers, and it converges with real in-vivo data).
6. **Summed across digits**: a dynamometer/handle measures the sum of finger (+thumb) contact
   forces. Real per-finger contribution splits (measured, not assumed): index/middle/ring/little =
   30/30/22/18% (Amis 1987, n=17, stable across cylinder diameters) or 20.2/37.5/28.7/13.6%
   (Lee/Kong/Lowe 2009, n=46, at a fixed 45–65mm span range) — the two studies disagree on which
   finger dominates (index≈middle in Amis vs. middle-dominant in Lee et al.), a genuine unresolved
   discrepancy, disclosed not smoothed over.
7. **Length-tension governs the aperture-dependence** (§3): as grip aperture changes, joint angles
   change, MTU length changes, and the fiber's operating point moves along the classic
   ascending-limb/plateau/descending-limb curve — producing the inverted-U.

---

## 2. The twin's own hand model — machine-measured this session, not recalled from prior prose

Two live model forks exist; **neither has both full digit coverage AND the pulley fix** (a genuine,
newly-surfaced gap this session, not previously stated in this form in any prior doc):

| Model file | Digits w/ FDP+FDS | Pulley-corrected moment arms | Muscle count |
|---|---|---|---|
| `data/msk_models/subject2_unified_v2_thumb.osim` | index(2), middle(3) + thumb (FPL) | no | 246 |
| `data/msk_models/subject2_unified_v2_ringlittle.osim` | index,middle,ring(4),little(5) + thumb | **no** (forked from thumb.osim's pre-pulley geometry) | 257 |
| `data/msk_models/subject2_unified_v2_pulleys.osim` | index, middle + thumb ONLY | **yes** (13/13 moment-arm-sign PASS, `docs/MECHANISM_HAND_PULLEYS.md`) | 246 (forked from `thumb.osim`, NOT `ringlittle.osim`) |

Fmax/optimal-fiber-length/tendon-slack-length, extracted live by direct XML parse of the `.osim`
files (not from any doc's prose):

| Muscle | Fmax (N) | Optimal fiber length (mm) | Pennation (deg) |
|---|---:|---:|---:|
| FDP2_r, FDP3_r, FDP4_r, FDP5_r (index/mid/ring/little profundus) | **45** each | 60.0 (identical, all 4) | 8.0 |
| FDS2_r, FDS3_r, FDS4_r, FDS5_r (…superficialis) | **40** each | 60.0 (identical, all 4) | 8.0 |
| FPL_r (thumb long flexor) | 45 | 60.0 | 8.0 |
| APB_r / FPB_r / OP_r / AdP_r (thumb intrinsics) | 11 / 9 / 17 / 30 | — | — |
| LUM1_r / DI1_r (finger intrinsics present in the pulleys fork) | 0.2 / 3.2 | 60.0 / 45.0 | — |

**Two real, specific, citation-backed discrepancies flagged (not fabricated numbers, a genuine
disclosed gap):**
- **Fmax is uniform across all 4 digits** (45N FDP / 40N FDS regardless of index/middle/ring/little)
  — a generic-template simplification. Lieber et al. 1992 (PMID 1401782, n=154 cadaveric forearm
  muscles) found the deep/superficial digital flexors "very similar architecturally, **with the
  exception of the small finger flexor digitorum superficialis, which was much smaller and shorter
  than the rest of the digital flexors**" — i.e., the twin's FDS5_r=40N (same as FDS2-4) is
  specifically contradicted by the primary architecture literature; the uniform-FDP assumption is
  comparatively better supported by the same quote.
- **Optimal fiber length is an identical round 60.0mm for every one of the 5 flexors** — a
  template/generic value, not subject- or muscle-specific architecture (real fiber lengths differ
  measurably by muscle per Lieber 1992's own methodology, laser diffraction on cadaveric
  specimens).

Post-pulley moment arms (index finger FDP2_r, live-measured in `docs/MECHANISM_HAND_PULLEYS.md`,
cross-read not re-derived): MCP 1.29–9.18mm across ROM [−25.7°,86.0°]; PIP 2.14–6.18mm across
[−13.7°,97.2°]; DIP 3.77–5.77mm across [−6.0°,81.6°]. Phalanx segment lengths, extracted live this
session from the joint-frame translations in the same `.osim` (`pip2_flex_r_joint`'s parent frame
on `index_proximal_r`, `dip2_flex_r_joint`'s parent frame on `index_medial_r`): **proximal phalanx
= 51.696mm, middle phalanx = 29.540mm** (both model-measured); distal phalanx length is **not**
model-measured (no fingertip marker exists — disclosed gap, confirmed live in
`MECHANISM_ELBOW_WRIST_FORCE.md`), approximated at 17mm from standard hand-anthropometry tables
(flagged as literature-typical, not model-measured).

---

## 3. Claim 1 — grip force vs. aperture is an inverted-U (length-tension mechanism)

**Mechanism**: FDP/FDS moment-arm sign convention (flexor: flexion shortens the MTU — already
established and pulley-stabilized in `docs/MECHANISM_HAND_PULLEYS.md`) means an **open** hand
(wide aperture) puts the flexor MTU near its longest length; a **fully closed fist** (narrow
aperture) puts it near its shortest. Both extremes push the sarcomere operating point off the
length-tension plateau (descending limb when long/open, ascending limb + reduced mechanical
leverage when short/closed); an intermediate aperture keeps fibers near-optimal — producing the
inverted-U.

**Geometric plausibility check** (using the twin's own numbers, not assumed): total MTU excursion
for FDP2_r across its own full disclosed ROM, computed as Σ(r_avg × ROM_in_radians) per joint:
MCP≈10.20mm + PIP≈8.05mm + DIP≈7.29mm = **25.55mm total**, against OFL=60mm → **42.6% of optimal
fiber length**. This is the right order of magnitude to traverse a meaningful fraction of a
standard skeletal active force-length curve (active-force-bearing width typically ~1.0–1.5×L₀) —
**geometrically consistent with, not proof of**, the mechanism being strong enough to produce a
real inverted-U. Explicitly NOT claimed: an aperture(cm)→joint-angle(deg) kinematic map has not
been built in this twin, so this is a plausibility bound, not a literal reproduction of the
empirical curve (see §6, proposed cell).

**Forced adversary (machine-checked, not narrated)**: is the argmax of grip force vs. span
INTERIOR to the tested range (falsifying both "flat" and the stronger "monotonic" adversary) in
each independent study?

| Study | n | Design | Argmax location | Interior (not at either tested extreme)? |
|---|---:|---|---|---|
| Firrell & Crain 1996 (PMID 8724468) | 288 hands | 5 fixed Jamar settings I–V | 89% (256/288) peak at setting **II** | **YES** — II ≠ {I, V} |
| Ruiz-Ruiz et al. 2002 (PMID 12239682) | 70 (40F/30M) | 5 grip spans × 10 trials | optimal span identified both sexes; men fixed at 5.5cm | **YES** — explicit rejection of "no optimal span exists" null |
| Lee, Kong, Lowe & Song 2009 (PMID 19424925) | 46 males | 5 spans, 45–65mm | peak 50mm (430.8N)/55mm (433.6N); 65mm lowest | **YES** — peak ≠ {45mm, 65mm} |
| Blackwell, Kornatz & Heath 1999 (PMID 10484275) | 18 | 4 grip spans, submax 60–65%MVC | "middle grip sizes > small OR large" | **YES**, by the paper's own explicit finding |

**4/4 PASS** across 4 independent populations/dynamometers (US, Spain, Korea/US-mixed), total
n=422. The monotonic adversary is FORCED (it is the mechanically-motivated stronger alternative to
"flat," since a pure lever-arm argument alone could in principle produce a monotonic curve) and
FALLS in all 4 instances — a diverse instance-space, satisfying the falsifier.

Amis 1987 (PMID 3682795, n=17, cylinders 31–116mm) shows a **monotonic decrease** within its own
tested range — not contradictory: 31mm is already near/past the ~50–55mm empirical peak's
descending side is not directly comparable (cylinder diameter ≠ handle span metric), so this is
read as **corroborating the descending limb**, not as a 5th independent full-U demonstration (kept
separate to avoid overclaiming).

**Peak location, numeric**: converges on **5.0–5.5cm** (Ruiz-Ruiz: men fixed 5.5cm, women
`y=x/5+1.5cm`; Lee et al.: 50–55mm) — this is a **refinement** of the task's stated 5–7cm prior to
a tighter, lower band, backed by 2 independent quantitative sources.

Sancho-Bru, Giurintano, Pérez-González & Vergara 2003 (PMID 14605652) — an independent
**biomechanical-model** (not measurement) prediction of 33mm optimal cylinder diameter — same
order of magnitude, different metric (diameter vs. span), a supporting theoretical cross-check.

---

## 4. Claim 2 — max grip magnitude reconstructed from the twin's own Fmax vs. real norms

**Step 1 — tendon-level budget** (machine-summed from §2's table): extrinsic flexors alone
(4×FDP@45N + 4×FDS@40N + FPL@45N) = **385.0N**. Pennation (8.0°) reduces this by only
`1−cos(8°)=0.97%` — negligible, checked not assumed. A generous "all closing-relevant muscles"
budget (+ thumb intrinsics APB/FPB/OP/AdP + built finger intrinsics LUM1/DI1) = **455.4N**.

**Step 2 — geometric transmission ratio (tendon force needed per unit fingertip force), derived
from the twin's OWN measured r and L, not assumed**: for a single joint, τ=F_tendon·r=F_tip·L (a
lever-balance approximation) ⟹ F_tendon/F_tip = L/r:

| Joint | r (twin-measured, mid-ROM) | L (phalanx length) | F_tendon/F_tip |
|---|---:|---:|---:|
| MCP | ~6.0mm | 51.696mm (model-measured) | **8.62×** |
| PIP | ~5.0mm | 29.540mm (model-measured) | **5.91×** |
| DIP | ~5.0mm | 17mm (literature estimate, not model-measured) | **3.40×** |
| Combined 3-joint (coordinated closure, Σr/ΣL) | 16mm | 98.2mm | **6.14×** |

**External cross-check**: Schuind, Garcia-Elias, Cooney & An 1992 (PMID 1564277) measured REAL
in-vivo tendon force via buckle transducers in 5 patients during tip pinch: tendon force up to
12.0kgf (117.7N) against a mean applied pinch force of 3.5kgf (34.3N) → **ratio 3.43×** — within 1%
of the twin's own **geometrically-derived DIP-joint estimate (3.40×)**, computed independently from
the twin's own moment arm and a literature phalanx length. This convergence (two independent
methods: direct geometric calculation vs. real in-vivo measurement) is a genuine cross-check,
reported with the caveat that some of the closeness may be coincidental given the approximations
involved (a representative mid-ROM r, a non-model-measured L_distal). Amis 1987's finding that
"distal segments of the fingers imposed forces significantly larger than those of the middle and
proximal segments" supports DIP-dominance applying to power grip as well as tip pinch, not just
pinch specifically.

**Step 3 — predicted max grip, all combinations** (2 budgets × 4 ratio-bases, the full sweep, not
a single cherry-picked number):

| Budget \ Ratio-basis | DIP (3.40×) | Combined (6.14×) | PIP (5.91×) | MCP (8.62×) |
|---|---:|---:|---:|---:|
| Extrinsic-only, 385.0N | 113.2N | 62.7N | 65.1N | 44.7N |
| Generous, 455.4N | 133.9N | 74.2N | 77.1N | 52.8N |

**All 8/8 combinations fall in 44.7–133.9N.**

**Step 4 — decorrelated external anchor**: Dodds et al. 2014 (PMID 25474696, PLoS One, n=49,964
participants/60,803 observations, 12 British studies) — peak median grip: **51kg men (500.1N)**,
**31kg women (304.0N)**. EWGSOP2 (Cruz-Jentoft et al. 2019, PMID 30312372, Table 3, verified via
direct PMC full-text fetch) clinical low-strength cutoffs: **<27kg men (264.8N)**, **<16kg women
(156.9N)** — the LOWEST clinically-meaningful real threshold (already-pathological weakness),
the most charitable comparison available for the "twin's budget is sufficient" adversary.

**Verdict**: every one of the 8 tested combinations (44.7–133.9N) is BELOW the EWGSOP2 **female**
cutoff (156.9N), let alone the male cutoff (264.8N) or either sex's healthy median (304–500N). The
adversary "the twin's current hand-flexor Fmax parametrization can plausibly reach realistic max
grip" was steelmanned with its most generous budget and most favorable geometric ratio and still
FALLS, robustly, across the whole tested assumption-space. **This is the opposite-direction finding
from `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`'s knee/hip result** (there: twin over-strong,
mean ratio 1.146 vs. Handsfield-MRI-implied Fmax; here: twin's hand plausibly under-strong by a
factor of ~2–6×) — a genuinely different, non-redundant result for a different body region built
on a different model lineage (uniform round-number Fmax/OFL vs. the knee's Rajagopal/Arnold-lineage
scaling), not a contradiction.

**What this is NOT**: a certified, precise number — no static-optimization simulation was run this
session (see §6). This is a bounded, multiply-cross-checked, geometrically-derived HYPOTHESIS for
independent QC, not a final measurement.

---

## 5. The dysfunction pole — grip strength as a validated function↔dysfunction axis

**Mortality (decorrelated, large-N, prospective)**: Leong, Teo, Rangarajan, Lopez-Jaramillo et al.
2015 (PURE study, Lancet, PMID 25982160), n=139,691 across 17 countries, median follow-up 4.0
years: grip strength (per 5kg **decrease**) — all-cause mortality HR **1.16 (95%CI 1.13–1.20,
p<0.0001)**; cardiovascular mortality **1.17 (1.11–1.24)**; non-CV mortality **1.17 (1.12–1.21)**;
MI **1.07 (1.02–1.11)**; stroke **1.09 (1.05–1.15)**. The paper's own explicit finding: grip
strength was a **stronger predictor of all-cause and CV mortality than systolic blood pressure**.
A companion paper from the same PURE cohort, Leong et al. 2016 (J Cachexia Sarcopenia Muscle,
PMID 27104109, n=125,462), provides cross-country reference ranges (a distinct paper — do not
conflate the two Leong citations, a real distinct-paper trap caught this session).

**Sarcopenia/frailty (diagnostic axis)**: Cruz-Jentoft et al. 2019 (EWGSOP2, Age Ageing, PMID
30312372) — grip strength is the primary case-finding measure for "low muscle strength," the first
diagnostic step; Table 3 cutoffs **<27kg (M) / <16kg (F)** confirmed via direct PMC fetch this
session (PMID 31081853 is the erratum notice for this same paper, not a separate source). Dodds et
al. 2014 (PMID 25474696) is the companion normative-curve anchor: weak grip (≥2.5SD below
sex-specific peak) reaches 23% (M) / 27% (F) prevalence by age 80 — the aging decline this
diagnostic axis is built to catch.

**Carpal tunnel syndrome (neuromuscular dysfunction pole)**: median nerve compression under the
transverse carpal ligament affects the thenar motor branch (± FDS sensory/proprioceptive
involvement depending on compression level) — mechanism reviewed in Katz & Simmons 2002 (NEJM
Clinical Practice, PMID 12050342, a standard clinical review with no numeric abstract body — its
citation is for mechanism, not a quantified number). A real, quantified, but **mechanistically
distinct** CTS-grip finding: Lowe & Freivalds 1999 (Ergonomics, PMID 10204420, n=7 CTS + 7 matched
controls) found CTS patients used **54% higher grip-to-load-force ratio** (p<0.05) and **12% lower**
force modulation (p<0.05); both a **coordination/sensory** deficit, not a raw maximum-strength
deficit — explicitly NOT conflated with motor weakness here. **Honest gap**: a specific, well-powered,
quantitative "CTS reduces max grip force by X kgf/N vs. healthy controls" primary study was
searched for but not found/verified this session — disclosed, not fabricated. Katz et al. 1998
(Maine Carpal Tunnel Study, J Hand Surg Am, PMID 9708386, n=429) shows surgical treatment produces
23–45% improvement in composite symptom/functional-status scores, sustained to 30 months — real
quantitative outcome data, but a composite score, not isolated grip kgf.

---

## 6. Honest gaps (disclosed, not chased further — no p-hacking)

1. **No aperture(cm)→joint-angle(deg) kinematic map exists in this twin.** §3's plausibility check
   uses the twin's own ROM-derived excursion, not an actual span-to-angle function; the exact
   empirical curve (peak at 5.0–5.5cm specifically) was NOT reproduced end-to-end from first
   principles, only shown geometrically plausible in order of magnitude.
2. **No single model file has both full 5-digit coverage AND the pulley-corrected moment arms**
   (newly surfaced this session — `pulleys.osim` forked from the pre-ring/little `thumb.osim`
   lineage, not from `ringlittle.osim`).
3. **Fmax/OFL are uniform generic-template values across all 4 finger-flexor digits** — contradicted
   specifically for FDS5 (little finger) by Lieber et al. 1992's own architecture data (real
   finding, not invented); no hand-specific Fmax/PCSA validation cert exists yet analogous to
   `docs/MECHANISM_FMAX_PCSA_VALIDATION.md` (knee/hip).
4. **Distal phalanx length is not model-measured** (no fingertip marker in this hand model,
   consistent with `MECHANISM_ELBOW_WRIST_FORCE.md`'s own disclosed gap); §4 used a
   literature-typical 17mm estimate, flagged wherever used.
5. **§4's transmission-ratio derivation is a single-joint lever approximation** (and one simple
   coordinated-closure sum), not a true multi-joint Jacobian/static-optimization solve — the
   44.7–134N prediction band is a bounded HYPOTHESIS, not a certified simulation result.
6. **No specific quantitative CTS-vs-healthy max-grip-force deficit citation was verified** this
   session (§5) — a real, disclosed negative, not papered over.
7. **A dedicated CTS/hand-dysfunction grip-force video or sensor dataset was not located** in this
   repo this session (the one grip-relevant video asset found,
   `data/video_index_tables/.../Tian_Tao_Floating_Snatch_Grip_Deficit_Deadlifts...`, is a
   weightlifting grip-failure clip, a function-pole not dysfunction-pole asset — see COUPLES_TO).

## 7. Proposed next cell

`scripts/msk/hand_grip_static_optimization.py`: (a) merge the pulley fix into the `ringlittle`
lineage (or vice versa) to get one model with both full-digit coverage and pulley-corrected moment
arms; (b) add a fingertip contact marker/frame (closing the "no grip marker exists" gap cited 3×
independently this session); (c) build a real aperture(cm)→joint-angle kinematic map (e.g. from
Buchholz & Armstrong-style hand-closure kinematics) and run OpenSim static optimization at 5+
apertures spanning 3.5–9cm to directly reproduce the Ruiz-Ruiz/Lee-Kong-Lowe curve shape AND get a
non-approximated max-grip-magnitude number to replace §4's bounded estimate.

---

## Citations (all PMIDs live-verified via NCBI eutils this session; DOIs as returned by PubMed)

| # | Citation | PMID | DOI |
|---|---|---|---|
| 1 | Dodds RM et al. Grip strength across the life course: normative data from twelve British studies. PLoS One. 2014;9(12):e113637. | 25474696 | 10.1371/journal.pone.0113637 |
| 2 | Leong DP et al. Prognostic value of grip strength: findings from the PURE study. Lancet. 2015;386(9990):266-73. | 25982160 | 10.1016/S0140-6736(14)62000-6 |
| 3 | Leong DP et al. Reference ranges of handgrip strength from 125,462 healthy adults in 21 countries (PURE study). J Cachexia Sarcopenia Muscle. 2016;7(5):535-546. | 27104109 | — |
| 4 | Cruz-Jentoft AJ et al. Sarcopenia: revised European consensus on definition and diagnosis (EWGSOP2). Age Ageing. 2019;48(1):16-31. | 30312372 | 10.1093/ageing/afy169 |
| 5 | (Erratum for #4) Age Ageing. 2019;48(4):601. | 31081853 | 10.1093/ageing/afz046 |
| 6 | Lieber RL, Jacobson MD, Fazeli BM, Abrams RA, Botte MJ. Architecture of selected muscles of the arm and forearm. J Hand Surg Am. 1992;17(5):787-98. | 1401782 | 10.1016/0363-5023(92)90444-t |
| 7 | Amis AA. Variation of finger forces in maximal isometric grasp tests on a range of cylinder diameters. J Biomed Eng. 1987;9(4):313-20. | 3682795 | 10.1016/0141-5425(87)90079-3 |
| 8 | Blackwell JR, Kornatz KW, Heath EM. Effect of grip span on maximal grip force and fatigue of FDS. Appl Ergon. 1999;30(5):401-5. | 10484275 | 10.1016/s0003-6870(98)00055-6 |
| 9 | Ruiz-Ruiz J, Mesa JL, Gutiérrez A, Castillo MJ. Hand size influences optimal grip span in women but not in men. J Hand Surg Am. 2002;27(5):897-901. | 12239682 | 10.1053/jhsu.2002.34315 |
| 10 | Lee SJ, Kong YK, Lowe BD, Song S. Handle grip span for optimising finger-specific force capability as a function of hand size. Ergonomics. 2009;52(5):601-8. | 19424925 | 10.1080/00140130802422481 |
| 11 | Firrell JC, Crain GM. Which setting of the dynamometer provides maximal grip strength? J Hand Surg Am. 1996;21(3):397-401. | 8724468 | 10.1016/S0363-5023(96)80351-0 |
| 12 | Sancho-Bru JL, Giurintano DJ, Pérez-González A, Vergara M. Optimum tool handle diameter for a cylinder grip. J Hand Ther. 2003;16(4):337-42. | 14605652 | 10.1197/s0894-1130(03)00160-1 |
| 13 | Schuind F, Garcia-Elias M, Cooney WP 3rd, An KN. Flexor tendon forces: in vivo measurements. J Hand Surg Am. 1992;17(2):291-8. | 1564277 | 10.1016/0363-5023(92)90408-h |
| 14 | An KN, Ueba Y, Chao EY, Cooney WP, Linscheid RL. Tendon excursion and moment arm of index finger muscles. J Biomech. 1983;16(6):419-425. | 6619158 | 10.1016/0021-9290(83)90074-x |
| 15 | Katz JN, Simmons BP. Clinical practice. Carpal tunnel syndrome. N Engl J Med. 2002;346(23):1807-12. | 12050342 | 10.1056/NEJMcp013018 |
| 16 | Katz JN et al. Maine Carpal Tunnel Study. J Hand Surg Am. 1998;23(4):697-710. | 9708386 | 10.1016/S0363-5023(98)80058-0 |
| 17 | Lowe BD, Freivalds A. Effect of carpal tunnel syndrome on grip force coordination on hand tools. Ergonomics. 1999;42(4):550-64. | 10204420 | 10.1080/001401399185469 |

## COUPLES_TO

- **Wrist/carpal force** (in flight): `docs/MECHANISM_ELBOW_WRIST_FORCE.md` — same finger-flexor
  muscles (FDP2-5_r/FDS2-5_r/FPL_r) already identified there as the ONLY wrist-flexion agonists
  (no dedicated FCR/FCU/palmaris longus in this model); its own 2kg-hold case (248.4N wrist
  contact force, 0.787 max activation, no saturation) is a DIFFERENT task (static support, not
  active power grip) from this doc's grip-force reconstruction — not conflated.
- **Forearm/specific-tension/PCSA**: `docs/MECHANISM_FMAX_PCSA_VALIDATION.md` (knee/hip,
  over-strong finding) — this doc's hand finding is the mirror-image (under-strong), same method
  family (Fmax-vs-real-anchor), opposite direction, different lineage.
- **Hand model lineage**: `docs/MECHANISM_ANATOMICAL_HAND.md`, `MECHANISM_HAND_PULLEYS.md`,
  `MECHANISM_HAND_FIVE_DIGIT.md`, `MECHANISM_HAND_FOREARM_THUMB.md`, `MECHANISM_HAND_INHERIT.md` —
  all read this session; the "no grip-force instrument" gap was independently stated in two of
  them and is the direct motivation for this doc.
- **Operator's lifting grip** (weightlifting): `data/video_index_tables/.../Tian Tao Floating
  Snatch Grip Deficit Deadlifts Almaty 2014 [3af8Ri9t-hw]...` — an existing real-world video asset
  literally titled around grip being the lift-limiting factor; a candidate real-world validation
  clip for a future grip-force-limited-lift cell (function pole, not yet used for that purpose).
- **Sarcopenia/frailty** (dysfunction pole): `data/body_twin/agent_outputs/frailty-multisystem-aging__adba752834918a1ca.json`
  exists as a prior related note (itself an LLM output, not primary data — mentioned as a pointer,
  not cited as evidence).
- **Function↔dysfunction organizing axis**: this doc is structured explicitly around both poles
  (§3-4 function/mechanism, §5 dysfunction/biomarker) per the operator's cross-cutting axis.
