# MECHANISM LANDING / DECELERATION MECHANICS — the eccentric energy-absorption FUNCTION pole, and its ACL-injury DYSFUNCTION pole (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** Literature-forced synthesis + this session's own
machine-computed physics cross-checks (impulse-momentum, energy balance) — **no new OpenSim/JAM
solve executed this pass**. Every number below traces either to a PMID/DOI abstract fetched **live**
this session via NCBI eutils (not recalled — recall drifts) or to
`scripts` output in `data/msk_smoketest/landing_mechanics/landing_mechanics_compute.py`, never hand
arithmetic. Full machine-readable ledger: `docs/MECHANISM_LANDING_MECHANICS_evidence.json`.

This is the **second** of the two canonical non-contact-ACL-injury mechanisms, complementing
`docs/MECHANISM_CUTTING_ACL_MECHANISM.md` (the sidestep-cut mechanism, landed concurrently this
session by the sibling instance). That cert's falsifier is about **plane of loading** (multiplanar
vs. sagittal-alone); this cert's falsifier is about **landing technique/stiffness** (how the same
touchdown momentum gets absorbed) — a different axis of the same FUNCTION↔DYSFUNCTION archetype
(`bt_memory/function-vs-dysfunction-build-the-athlete-baseline-in-parallel-with-disease.md`). The
two share 3 dysfunction-pole anchors (Hewett 2005, Ford 2003, McLean 2007); this document
independently re-fetched all 3 from PubMed rather than copying the sibling cell's numbers, and they
match exactly (a genuine cross-check that neither transcribed a number wrong), plus adds one fresh
decorrelated triangulation (Pappas 2007) the sibling cell did not use.

---

## 0. What this is / is not

**Is**: a geometric/physics-first model of landing as an eccentric energy-absorption task
(ankle→knee→hip), a falsifier-forced test of whether landing **stiffness/technique** (not just
drop height) governs peak force and its joint-by-joint distribution, and the coupling of that
FUNCTION-pole model to the DYSFUNCTION pole (stiff/extended/valgus landing → ACL loading) already
partly anchored in the sibling cutting-ACL cert.

**Is not**: a new subject-specific OpenSim landing simulation on this twin's own geometry (no drop-jump
trial with ground-truth GRF exists in this repo's currently-loaded corpus at the fidelity this task
needs — see §8.1); a re-derivation of the cutting-task's multiplanar-loading mechanism (that is
`MECHANISM_CUTTING_ACL_MECHANISM.md`'s job, not repeated here); a claim that any single-timepoint
landing-stiffness screen reliably predicts *which* individual gets injured (§6 carries the same
contested-classifier caveat the sibling cell already established for KAM, applied consistently here).

---

## 1. Pre-registration (set before tabulating the specific numbers below)

**C1 (impulse-momentum):** for a **fixed** touchdown momentum (fixed drop height), peak vertical GRF
and loading rate are **not** fixed — they scale with landing stiffness/technique (shorter effective
deceleration time → higher peak force), reproducible across ≥2 independent, decorrelated studies at
matched height. **¬C1 (the adversary):** peak GRF is a function of height/momentum **alone**; at
matched height, technique/stiffness makes no measurable difference (within-study effect
non-significant or <15%).

**C2 (joint redistribution):** the eccentric (negative) work absorbing that momentum is shared across
ankle+knee+hip, and the **share** (not just the total) shifts with stiffness — hip+knee dominant when
soft/flexed, ankle relatively more engaged when stiff/erect — replicated across ≥2 decorrelated
labs/years. **¬C2:** the per-joint percentage split is invariant to technique (a fixed proximal-distal
ratio regardless of how the landing is performed).

**C3 (dysfunction pole, extending — not re-deriving — the sibling cert's anchors):** a stiff,
extended-knee, valgus-collapsing landing loads the ACL more than a soft, flexed, neutral one, via
(a) the higher peak force/loading rate C1 predicts, and (b) the knee sitting closer to the
angular range this twin's **own** ligament model (`MECHANISM_KNEE_LIGAMENTS.md`) already found the ACL
is most taut in. **¬C3:** landing kinetics are unrelated to the prospective ACL-injury literature's
own measured group differences.

**Forced adversary:** "landing force is fixed regardless of technique" — forced to its strongest fair
form as a same-height (momentum-matched), multi-technique comparison (§4), not a cross-height
confound. **External anchor (decorrelated, never a tautology gate):** DeVita & Skelly 1992, Zhang/
Bates/Dufek 2000, Self & Paine 2001, Peng 2011, Decker et al. 2003 (FUNCTION pole); Hewett 2005, Ford
2003, McLean 2007, Pappas et al. 2007 (DYSFUNCTION pole) — 9 primary sources, all PMID/DOI-verified
live this session (§10), none of them tuned to or derived from this twin's own model.

---

## 2. Method

No live biomechanics trial in this repo currently has a ground-truth drop-landing GRF trace at
sufficient fidelity (§8.1) — so this cell is a **literature-forced mechanism synthesis**, exactly the
same genre as the sibling cutting-ACL cert, plus this session's own derived-number computations
(never asserted without a script trace):
`data/msk_smoketest/landing_mechanics/landing_mechanics_compute.py` (stdlib-only, re-runnable),
output `data/msk_smoketest/landing_mechanics/landing_mechanics_results.json`. Literature: NCBI eutils
`esearch`→`efetch`(abstract) against PubMed, live this session, `curl --max-time 25` every call — every
PMID below was fetched, not recalled. Two of the three dysfunction-pole citations (Hewett 2005, Ford
2003, McLean 2007) are shared with `MECHANISM_CUTTING_ACL_MECHANISM.md`; this document re-fetched them
independently (fresh esearch/efetch calls, not a copy of the sibling doc's text) as its own
spot-check — see §6 for the cross-check result.

---

## 3. The geometry: touchdown kinetic energy, cross-checked against a real dataset's own numbers

Touchdown vertical speed from free fall: **v = √(2gh)**. Zhang, Bates & Dufek 2000 (PMID 10776901,
DOI 10.1097/00005768-200004000-00014; *Med Sci Sports Exerc*, n=9 male subjects) report **both** the
drop heights **and** the landing velocities they used — an external number this twin did not
generate, so predicting it from `√(2gh)` alone and comparing is a genuine, decorrelated,
machine-computed cross-check, not a tautology:

| height (m) | v = √(2gh), computed | v reported (Zhang 2000) | % diff |
|---:|---:|---:|---:|
| 0.32 | 2.506 m/s | 2.5 m/s | +0.23% |
| 0.62 | 3.488 m/s | 3.5 m/s | −0.35% |
| 1.03 | 4.495 m/s | 4.5 m/s | −0.10% |

**All 3 within 0.4%** — the free-fall model is exactly right at these scales (as it must be — this is
the decisive, load-bearing physics identity underneath everything else in this document, and it is
confirmed against numbers this twin did not choose).

**Energy balance, refined beyond the task's own first-order framing.** KE at touchdown is `½v²` per
kg. But the task's framing ("total negative work = KE at touchdown") implicitly assumes the COM stops
at the *same height* it touched down at. In a real landing the COM continues to descend by some
crouch depth `Δy` while decelerating, so gravity does **additional** positive work `mgΔy` during the
absorption phase that the eccentric muscles must *also* remove:

**|W_eccentric| = ½mv²_touchdown + mgΔy_crouch — always ≥ ½mv²_touchdown, not equal to it.**

This predicts the simple `½v²` figure is a **lower bound** on required eccentric work, and the bound
should be *tighter* (closer to true) for a **stiff** landing (small `Δy`, little knee travel) than for
a **deep, soft** one (large `Δy`) — i.e., a soft landing's muscles do *more* total work than naive KE
alone suggests, not merely the *same* total work spread differently. No trial in the fetched
literature reports `Δy_crouch` directly, so this refinement is not numerically closed here — it is a
geometric correction that sharpens, not overturns, the task's model, and it means the ratio computed
next is an upper bound on the fraction of *true* required work the measured joint quantities capture.

---

## 4. Falsifier — does technique/stiffness modulate peak force at FIXED momentum? Adversary forced and rejected

The **forced adversary** ("peak GRF is fixed by height/momentum alone") is tested at its strongest —
same height, technique varied — in 3 independent, decorrelated studies:

1. **Self & Paine 2001** (PMID 11474336, DOI 10.1097/00005768-200108000-00015; *Med Sci Sports Exerc*)
   — the cleanest test: **one fixed height** (30.48 cm), **4 landing techniques** (bent-knee natural;
   stiff-knee natural; stiff-knee-plantarflexor-absorbed; stiff-knee-heel-absorbed). Peak vertical
   force and peak tibial acceleration were **highest for the stiffest/heel-focused technique**: **2418
   N, 20.7 G** — by the paper's own explicit ranking, not every technique produced the same force at
   the same momentum. **¬C1 falls here directly**: same `mΔv`, different peak F.
2. **Zhang, Bates & Dufek 2000** (PMID 10776901) — a **3-height × 3-technique** factorial (soft/
   normal/stiff at each height): "general increases in peak ground reaction forces, peak joint
   moments, and powers with increases in landing height **and stiffness**" — stiffness is reported as
   its **own** main effect, not absorbed into the height effect.
3. **DeVita & Skelly 1992** (PMID 1548984, DOI 10.1249/00005768-199201000-00018; *Med Sci Sports
   Exerc*, n=8 female athletes, single 59 cm drop height, soft vs. stiff technique only): soft landing
   averaged **117°** knee flexion, stiff landing **77°** — "the stiff landing had larger GRFs" (the
   abstract's own direct statement; exact peak-force values sit in the full-text table, not the
   abstract — disclosed gap, §8.2).
4. **Peng 2011** (PMID 21869631, DOI 10.1519/JSC.0b013e318201bcb3; *J Strength Cond Res*, n=16): peak
   vGRF of the dominant leg **exceeded 3× body weight** at 50 and 60 cm drop heights (significantly
   greater than 20/30/40 cm, p<0.004) — **while** directly-measured leg/knee/ankle **stiffness**
   **decreased** at 60 cm relative to the 3 lower heights (p<0.037): the body's own endogenous
   stiffness is a real, separately-modulated variable that correlates with the force outcome, not an
   epiphenomenon of height alone.

**Verdict: ¬C1 is rejected at 4/4 independent, decorrelated sources; C1 holds.** The impulse-momentum
geometry (`F_peak ∝ 1/T_contact` for fixed `mΔv`, from a half-sine contact-pulse model,
`F_peak = (π/2)·mv/T`) explains *why*: shorter effective contact/deceleration time — the operational
definition of "stiffer" — mechanically requires higher peak force for the same momentum change.
Illustrative-only computation (no measured `T` exists in the fetched abstracts, so absolute
scaling — not the trend — is illustrative), using DeVita & Skelly's 59 cm drop (v=3.402 m/s):

| illustrative contact time T | predicted F_peak / BW |
|---:|---:|
| 0.40 s (soft, long) | 1.36× |
| 0.10 s (stiffer, short) | 5.45× |
| 0.03 s (near-rigid) | 18.16× |

This bounds the task's own framing ("stiff/extended landing → ~3-6×BW") as physically reasonable for
intermediate contact times, and shows a truly **rigid** (near-zero flexion) landing predicts forces
far above any measured soft-landing value — consistent with why case reports of straight-leg landings
associate with fracture, not just soft-tissue strain. A convergent check: back-solving Self & Paine's
own reported 2418 N (given an **assumed**, not reported, 70 kg mass — flagged explicitly, the
abstract does not state subject mass) through the same half-sine model implies **T≈111 ms**, which
sits squarely in the known short-ground-contact-time range (comparable to sprint contact times) — a
physically sensible number recovered independently, not fitted to match.

---

## 5. Where the energy goes — joint redistribution, and how it shifts with stiffness

**Zhang et al. 2000**'s own reported mean eccentric work by joint and height (J/kg), and this
session's own computed percentage split:

| height (m) | ankle (J/kg) | hip (J/kg) | knee (J/kg) | total | %ankle | %hip | %knee | total / (½v²) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.32 | 0.52 | 0.94 | 1.21 | 2.67 | 19.5% | 35.2% | **45.3%** | 0.851 |
| 0.62 | 0.74 | 1.31 | 1.63 | 3.68 | 20.1% | 35.6% | **44.3%** | 0.605 |
| 1.03 | 0.87 | 2.15 | 2.26 | 5.28 | 16.5% | 40.7% | **42.8%** | 0.523 |

Knee is the largest single contributor at every height (42.8–45.3%, stable); hip's share **grows**
with height (35.2%→40.7%) while ankle's *relative* share shrinks slightly (19.5%→16.5%, even though
its *absolute* J/kg still rises) — exactly the paper's own stated conclusion ("a shift from ankle to
hip strategy… as landing height increased"), now reproduced as an explicit percentage rather than
narrated. The last column (measured total ÷ `½v²`) falls from 85%→52% with height — **an honest,
disclosed, only-partially-discriminating finding**: because §3 showed `½v²` is itself an
under-estimate of the true required work (missing `mgΔy_crouch`), a ratio below 100% is close to
guaranteed regardless of whether these sagittal, muscle-only, single-plane joint-power measurements
are capturing everything real — the falling-with-height *trend* is the more informative part
(candidate causes, none confirmed this pass: un-measured frontal/transverse-plane and trunk/arm work,
passive-tissue/SSC energy storage, incomplete "stick" at higher demand — an honest, unresolved gap,
§8).

**Independent replication of the redistribution pattern itself (not the height trend), via a
different manipulation:** **DeVita & Skelly 1992** (soft vs. stiff technique, single height): hip
absorbed more **power** in soft (−0.60 vs. −0.39 W/kg, p<0.01), knee more in soft (−0.89 vs. −0.61
W/kg, p<0.01), **ankle more in stiff** (−0.88 vs. −1.00 W/kg, p<0.05) — the same
proximal-dominant-when-soft / distal-shift-when-stiff pattern Zhang's height manipulation shows,
now confirmed by a **decorrelated** lab, era, and manipulation (technique, not height). Note these
are **power** (W/kg), a different metric type than Zhang's **work** (J/kg) — qualitatively convergent,
not numerically comparable term-for-term.

**A third decorrelated replication, via a sex manipulation — Decker, Torry, Wyland, Sterett & Steadman
2003** (PMID 12880714, DOI 10.1016/s0268-0033(03)00090-1; *Clin Biomech*, n=12 male/9 female, 60 cm
drop): knee is the primary shock absorber for **both** sexes, but the **second-largest** contributor
differs — **ankle** plantarflexors for females, **hip** extensors for males; females land in a more
erect posture. This is the FUNCTION-pole redistribution axis (§C2) converging, a third independent
way, with the DYSFUNCTION pole (§6): the same erect/distally-shifted strategy literature repeatedly
associates with elevated ACL risk in females.

**Verdict: C2 holds** — replicated across 3 decorrelated manipulations (height: Zhang; technique:
DeVita & Skelly; sex: Decker), same qualitative direction each time.

---

## 6. The DYSFUNCTION pole — reused anchors independently re-verified, plus one new triangulation

**Hewett TE et al. 2005** (PMID 15722287, DOI 10.1177/0363546504269591; *Am J Sports Med*, n=205
female athletes, 9 later ACL ruptures) — re-fetched live this session, independently of the sibling
cell: 8° greater knee abduction angle (p<.05), **2.5× greater** knee abduction moment (p<.001), **20%
higher** GRF (p<.05), **16% shorter** stance time in the 9 injured vs. 196 uninjured. **Exact match**
to the sibling `MECHANISM_CUTTING_ACL_MECHANISM.md`'s independently-fetched numbers — two separate
live NCBI pulls landing on identical figures is itself a small but real cross-check that neither
transcription is garbled.

This document does **not** re-litigate the sibling cell's own symmetric-QC finding and defers to it
rather than re-deriving: the **mechanism** (higher force/shorter time/greater valgus in the injured
group) is a real, replicated within-study finding, but **Croncurrent, Creaby & Ageberg 2020** (PMID
32819327, PMC7441716 — a 9-cohort prospective meta-analysis, already independently verified in
`MECHANISM_CUTTING_ACL_MECHANISM.md` §3) found pooled KAM angle/moment does **not** significantly
predict individual future ACL injury. **The same caveat applies here**: this document's C3 rests on
the mechanistic force/kinetics relationship, not on any claim that a single landing-stiffness
screening measurement classifies which individual gets hurt.

**Sex** — **Ford, Myer & Hewett 2003** (PMID 14523314, DOI 10.1249/01.MSS.0000089346.85744.D9, n=81:
47F/34M) re-fetched live, identical to the sibling cell: females land with significantly greater
total valgus motion and maximum valgus angle during a drop vertical jump.

**Fatigue** — **McLean et al. 2007** (PMID 17473777, DOI 10.1249/mss.0b013e3180d47f0, n=20 NCAA
athletes) re-fetched live, identical to the sibling cell: fatigue increases initial-contact and
peak-stance knee abduction/internal-rotation motion **and** moments, more pronounced in females.

**New this session, a fourth decorrelated triangulation not used by the sibling cell — Pappas,
Sheikhzadeh, Hagins & Nordin 2007** (PMID 24149228, PMCID PMC3778703; *J Sports Sci Med*, n=32
recreational athletes, 40 cm bilateral drop-jump landing): females landed with **9° greater** peak
knee valgus (p=0.001) **and 140% greater** maximum vGRF normalized to body weight than males
(p=0.003); fatigue increased peak vGRF by **20%** (p=0.038), foot abduction by 1.7° (p=0.042), and
rectus femoris peak activity by 27% (p=0.018). The valgus-angle sex effect (9°) **replicates** Ford
2003's direction almost exactly (independent lab, cohort, and era); the fatigue-raises-peak-force
finding independently replicates McLean 2007's fatigue direction using a **different** instrument
(peak vGRF + EMG, not 3D moments). The 140% GRF-magnitude sex effect is **flagged, not silently
accepted**: it is a large effect size with no second GRF-magnitude-specific sex-comparison study
cross-checked this pass (unlike the valgus-angle finding, now triangulated 3×) — reported as a
single-source data point pending independent replication, per this repo's own standing discipline
that an anchor is a hypothesis too.

**A geometric convergence with this twin's own already-built ligament model.**
`docs/MECHANISM_KNEE_LIGAMENTS.md` (§4 there) independently found, on this twin's own OpenSim
`Blankevoort1991Ligament` ACL bundles, that **ACL tension is highest near knee extension** (4.318 N
mean over 0–15° vs. 2.143 N over 125–140°, pre-registered and PASS in that document). This document's
own mechanism supplies a **second, independent** reason the same conclusion should hold: a stiff/
erect landing (§4–5) keeps the knee in the same near-extension range that document already measured
as the ACL's most-taut zone, on top of (not instead of) delivering a higher peak force in less time.
Two decorrelated mechanisms — impulse-momentum kinetics here, ligament-strain-vs-angle geometry there
— point the same direction. **Terminology note, not to be conflated**: "landing stiffness" (this
document's central variable — a modulable, technique-dependent mechanical property of the whole limb
during a dynamic task) and "ligament structural stiffness" (`MECHANISM_LIGAMENT_STIFFNESS.md`'s ACL
283 N/mm, a fixed material property) are related but distinct quantities; a stiffer *landing* routes
more load onto a fixed-stiffness *ligament* rather than dissipating it through compliant eccentric
muscle-tendon work — the two "stiffnesses" interact, they are not the same measurement.

**Derived (not literature-reported) loading-rate figure.** Combining Hewett 2005's own two ratios
(peak GRF +20%, stance time −16%, i.e., ×0.84) into a crude average-loading-rate proxy (peak
force/time-to-peak): **1.20 / 0.84 = 1.43×** — the ACL-injured group's landings imply roughly 43%
higher average loading rate than the uninjured group's, **computed by this document from Hewett's
own two published percentages**, not a number Hewett et al. themselves report.

---

## 7. Video-markerless coupling — this task's stiffness signature IS in-plane, unlike the cutting task's KAM

The sibling `MECHANISM_CUTTING_ACL_MECHANISM.md` §6 found an **honest negative** for frontal-plane
(KAM) precision from markerless video (OpenCap: frontal CMC only 0.47–0.78 vs. sagittal CMC>0.94; a
2026 ACL-re-injury-screening validation explicitly concluding markerless video "cannot be recommended
for ACL re-injury risk assessment" at the frontal-plane level) — reused here, not re-verified, per
this repo's own "don't re-litigate" convention.

**This document's own falsifier (§4–5) is built entirely from SAGITTAL-plane quantities**: knee
flexion depth (77–117° range, DeVita & Skelly), stance/contact time (Hewett's 16% figure), and a
GRF-magnitude proxy — precisely the variable class the sibling cell's own §6 already found is
video-tractable **today** (CMC>0.9, few-degree error), in direct contrast to the frontal-plane knee
abduction moment that is not. This is the concrete, actionable difference this document adds to the
graph: a **landing**-based stiffness flag (shallow knee flexion + short stance + a fast-rising
GRF-proxy, all sagittal) is a plausible video-estimable risk marker **today**, whereas a
**cutting**-based KAM flag is not — an asymmetry between the two canonical ACL mechanisms worth
carrying forward explicitly, not assumed away. This is a proposed, unvalidated pathway (§9), not a
built or tested video pipeline output this session.

---

## 8. Honest gaps

1. **No new raw-data solve executed this pass** — no OpenSim drop-landing simulation on this twin's
   own geometry; no ground-truth GRF trial for a controlled drop-landing exists in this repo's
   currently-loaded corpus at this fidelity (the closest asset, `MECHANISM_PLYO_JUMP_CASCADE.md`, uses
   Instagram-sourced ballistic jump clips selected for a *takeoff*-cascade-timing question, not a
   controlled-height drop-landing GRF measurement — a different question, not reused here as if it
   answered this one).
2. **DeVita & Skelly 1992's exact peak-GRF Newton/BW values are not in the fetched abstract** ("larger
   GRFs" is stated qualitatively; the abstract is a pre-2000-era MSSE abstract without a full results
   table) — not invented to fill the gap; the quantitative peak-force claims in this document instead
   rest on Self & Paine 2001, Peng 2011, and the physics-model illustrative table.
3. **The `Δy_crouch` correction (§3) is not numerically closed** — no fetched source reports COM
   crouch depth directly, so the "52–85% of `½v²`" ratio is presented as a bounded, partially
   informative plausibility check, explicitly flagged as having limited discriminating power on its
   own (§5), not as a tight validated energy balance.
4. **Self & Paine's 2418 N/20.7 G-implies-70 kg-subject and the T≈111 ms back-out both use an
   ASSUMED, not reported, subject mass** — flagged inline every time it is used; not to be read as a
   number the paper itself published.
5. **The Pappas 2007 "140% greater" GRF/BW sex effect is single-source this pass** (§6) — large
   effect size, not independently cross-checked against a second GRF-magnitude-specific
   sex-comparison study.
6. **The dysfunction-pole KAM-as-individual-injury-classifier caveat (Croncurrent 2020) is reused, not
   re-derived**, from the sibling cutting-ACL cert — appropriate per this repo's "don't re-litigate"
   convention, but this document did not independently re-fetch that meta-analysis itself.
7. **No frontal-plane (valgus/KAM) measurement exists for a pure drop-landing in this document** —
   the dysfunction-pole anchors (Hewett, Ford, McLean, Pappas) all use a drop-vertical-jump or
   bilateral-drop-jump protocol (landing + an immediate rebound jump), not an isolated "stick the
   landing" trial; this is the standard protocol in this literature, not a mismatch this document
   introduced, but it means the FUNCTION-pole "stick" energy-balance framing (§3) and the
   DYSFUNCTION-pole citations are not from identically-matched tasks.
8. **§7's video-tractability claim is proposed, not tested** — no video-derived knee-flexion/
   stance-time/GRF-proxy was actually computed against a ground-truth landing trial this session.

---

## 9. Couplings

- `docs/MECHANISM_CUTTING_ACL_MECHANISM.md` (+ its `_evidence.json`) — the other canonical non-contact
  ACL mechanism; shares 3 dysfunction-pole anchors (independently re-verified here, §6), and this
  document's stiffness/impulse-momentum falsifier is the complementary axis to that cell's
  plane-of-loading falsifier.
- `docs/MECHANISM_KNEE_LIGAMENTS.md` — this twin's own ACL `Blankevoort1991Ligament` bundles engage
  near extension (§6's geometric convergence).
- `docs/MECHANISM_LIGAMENT_STIFFNESS.md` — ACL structural stiffness 283 N/mm; the terminology
  distinction from "landing stiffness" (§6).
- `docs/MECHANISM_PLYO_JUMP_CASCADE.md` — the jump-cascade cell (the flight/propulsion side of the
  same movement family; this document is the descent/landing side).
- `docs/MECHANISM_JOINT_FORCE_SCORECARD.md`, `docs/MECHANISM_HIP_FORCE.md`,
  `docs/MECHANISM_ANKLE_FORCE.md`, `docs/MECHANISM_KNEE_CARTILAGE_MATERIAL.md` — the walking-anchored
  hip/knee/ankle joint-force family this document's landing-specific joint-redistribution model
  extends into a ballistic-impact regime.
- `docs/MECHANISM_POSE_PIPELINE.md`, `bt_memory/video-to-injury-load-two-barriers-id-double-
  differentiation-and-frontal-plane-knee-dof.md` — the video-markerless pipeline (§7); the same ID
  double-differentiation impact-window degradation that document found for the cutting/drop-jump
  chain applies equally to any future landing-trial video pipeline.
- `bt_memory/function-vs-dysfunction-build-the-athlete-baseline-in-parallel-with-disease.md` — the
  operator's own organizing-axis directive this document, and its sibling, execute.
- Proposed graph node: **`MSK-LANDING-DECELERATION`** (new — not yet in `data/MECHANISM_ANCHOR_GRAPH.json`;
  this session only writes the 4 files listed under **Files** below, per this repo's
  HYPOTHESIS-awaiting-QC convention — the coordinator adds the graph node after QC).

---

## 10. Citations verified live this session (PMID/DOI/n/measure — see evidence.json for exact quoted text)

| # | citation | PMID | DOI | n | measure |
|---|---|---|---|---|---|
| 1 | DeVita & Skelly 1992, *Med Sci Sports Exerc* | 1548984 | 10.1249/00005768-199201000-00018 | 8 (female) | soft (117°) vs. stiff (77°) knee-flexion landing, joint moments/power |
| 2 | Zhang, Bates & Dufek 2000, *Med Sci Sports Exerc* | 10776901 | 10.1097/00005768-200004000-00014 | 9 (male) | eccentric work (J/kg) by joint × 3 heights × 3 techniques |
| 3 | Self & Paine 2001, *Med Sci Sports Exerc* | 11474336 | 10.1097/00005768-200108000-00015 | not stated | peak vGRF/tibial accel across 4 landing techniques, fixed height |
| 4 | Peng 2011, *J Strength Cond Res* | 21869631 | 10.1519/JSC.0b013e318201bcb3 | 16 | peak vGRF & leg/knee/ankle stiffness across 5 drop heights |
| 5 | Decker, Torry, Wyland, Sterett & Steadman 2003, *Clin Biomech* | 12880714 | 10.1016/s0268-0033(03)00090-1 | 12M/9F | sex-based joint energy-absorption strategy |
| 6 | Hewett et al. 2005, *Am J Sports Med* | 15722287 | 10.1177/0363546504269591 | 205 (9 injured) | prospective KAM/valgus → ACL injury |
| 7 | Ford, Myer & Hewett 2003, *Med Sci Sports Exerc* | 14523314 | 10.1249/01.MSS.0000089346.85744.D9 | 47F/34M | sex difference, valgus motion, DVJ |
| 8 | McLean et al. 2007, *Med Sci Sports Exerc* | 17473777 | 10.1249/mss.0b013e3180d47f0 | 10M/10F | fatigue → valgus/IR motion+moment |
| 9 | Pappas, Sheikhzadeh, Hagins & Nordin 2007, *J Sports Sci Med* | 24149228 (PMC3778703) | — | 32 | sex + fatigue, peak valgus & vGRF/BW |

**Provenance note on row 1's DOI**: PubMed's own 1992-era record for DeVita & Skelly does not surface
a DOI in the fetched abstract text (pre-DOI-registration era). That DOI was recovered via a
**separate** live Crossref bibliographic-search API call this session (`api.crossref.org`, exact
title-match, single top result) — a second, decorrelated registry, not invented — and is reported
here distinctly from the other 8 rows, which came directly from the PubMed/NCBI eutils record itself.
Row 9 (Pappas 2007) has no DOI in either registry checked; PMCID is used as the locator instead.

---

## 11. Proposed next cell

Not a video→landing-stiffness regression yet (§7's claim is proposed, untested). Decisive next step:
(1) identify or acquire a controlled drop-landing trial with synchronized GRF + marker/video data
(this repo's own corpus is Instagram-sourced ballistic content, not controlled-height drop landings —
§8.1); (2) if found, reproduce Zhang/DeVita-style joint-power redistribution on this twin's own
geometry as a genuine new measurement (not literature-only, closing this cell's biggest honest gap);
(3) only then attempt the video-tractable sagittal stiffness-proxy screen §7 proposes, gated the same
way the sibling cell gates its own next step — inject realistic video-kinematic noise before trusting
any derived risk flag.

## Files

- `docs/MECHANISM_LANDING_MECHANICS.md` — this document.
- `docs/MECHANISM_LANDING_MECHANICS_evidence.json` — machine-readable mirror: every citation, the
  falsifier table, the computed velocity/energy/joint-redistribution numbers, honest gaps, couplings.
- `data/msk_smoketest/landing_mechanics/landing_mechanics_compute.py` — self-contained stdlib script;
  every "computed"/"derived" number in this document traces to its printed output.
- `data/msk_smoketest/landing_mechanics/landing_mechanics_results.json` — that script's raw output.

ISOLATION (`COORDINATOR.md` §1): bodytwin only; all coupled files read in place, never mutated; no git
commit/push; the 4 files above are the only new files this session.
