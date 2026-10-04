# MECHANISM CUTTING / NON-CONTACT ACL INJURY MECHANISM — the FUNCTION↔DYSFUNCTION archetype (2026-07-22)

Builds a quantitative, falsifier-tested model of the sidestep-cut ACL-injury mechanism: rapid
deceleration + foot plant + redirection, non-contact ACL rupture in early stance, near-extended
knee, under a **multiplanar** load (valgus/abduction + anterior shear + internal tibial rotation).
This is the operator's named FUNCTION↔DYSFUNCTION worked example: a validated **healthy-cut**
baseline (minimized knee-abduction moment via trunk/hip/hamstring control) is exactly what lets a
**dysfunctional** deviation (valgus collapse) be flagged from markerless video. Literature-forced
synthesis — no new OpenSim/JAM solve executed this pass (status **OPEN**, not PROVEN). Full
citation ledger: `docs/MECHANISM_CUTTING_ACL_MECHANISM_evidence.json` (every number below traces to
a PMID/DOI fetched **live** this session via NCBI eutils, not recalled).

**Builds on, does not re-litigate**, 4 existing graph nodes this repo already has: `MSK-CUTTING-ACL-LOAD`
(OPEN/partial — cutting-task KAM magnitude already byte-verified from a raw waveform by a prior
session), `MSK-KNEE-6DOF-JAM-ACL`/`AUTO-DO-NOT-REGRESS-VIDEO-KAM-INJURY-RISK-DIR` (SEED-DESIGN — the
JAM/COMAM ACL-force route + a noise-injection OODA gate, proposed but not executed), `VIDEO-TO-KNEEMOMENT-DROPJUMP-ACL`
(video→sagittal-knee-moment chain, 19-25% peak error), `VIDEO-LINK1-KINEMATICS` (this twin's
video→kinematics front-end). This document supplies the piece those seeds were missing: a
falsifier-forced mechanism synthesis with fresh, independently-verified primary literature.

---

## 0. Headline

**All 4 pre-registered falsifier sub-claims PASS on live-verified primary literature.** ACL
force/strain in a cut is a **multiplanar** (valgus + anterior-shear + internal-rotation) load, not a
sagittal one — confirmed both by direct cadaveric mechanics (Markolf 1995) and a validated
finite-element model (Navacchia/Hewett 2019), and DECISIVELY falsified in the sagittal-only
direction by McLean 2004's own pre-registered 2000N criterion (peak anterior-drawer force **never**
exceeded it across 5000 stochastic perturbations, while valgus loads reached injury-capable
magnitudes). Hamstring co-contraction is confirmed protective at the ACL in-situ-force level (Li
1999: 30-44% force reduction, p<0.05). The DYSFUNCTION-pole prospective anchor (Hewett 2005: 2.5x
greater KAM in later-injured athletes, p<.001) is real and reproduced here from the live abstract —
**but** its population-level replication as an individual-injury CLASSIFIER is genuinely contested:
a 9-cohort meta-analysis (Croncurrent 2020, PMC7441716) found pooled KAM does **not** significantly
predict future ACL injury (mean diff crosses zero). This is reported symmetrically, not buried: the
**mechanism** (combined loading → higher ACL force) is robust; the **single-timepoint-KAM-as-injury-
classifier** claim is not, and the two must not be conflated. Video-markerless coupling is an
**honest negative** for frontal-plane (KAM) precision, confirmed by 3 independently-converging
sources, while sagittal kinematics are video-tractable today.

---

## 1. The geometric mechanism (derived, not just cited)

The ACL is **one** oblique fiber bundle (femoral posterolateral condyle → tibial anteromedial
plateau) — not three independent structures, one per anatomical plane. Two geometric coupling paths
force multiplanar loads onto this single ligament's tension, and a sagittal-only free-body model
sets both to zero **by construction**:

1. **Tibial posterior slope.** The lateral tibial plateau slopes posteriorly-downward. An applied
   knee-abduction (valgus) moment increases lateral-compartment compressive contact force, and the
   slope's geometry converts part of that compression into an **additional anterior-directed shear**
   component at the lateral plateau. Navacchia & Hewett's validated finite-element models measure
   this directly: ACL force correlates at **r=0.99** with the anterior component of lateral-plateau
   contact force, and that component is increased by larger KAM *and* larger internal-rotation
   torque acting **through the posterior slope** (PMID 30977558). This anterior-shear component
   SUMS with whatever sagittal-plane anterior shear GRF/quadriceps pull already produce, on the same
   fiber bundle.
2. **Oblique-fiber wind-up.** Internal tibial rotation twists the tibia under the femoral condyles
   along the ACL's own long axis, tightening (winding) its oblique fibers independent of flexion
   angle — the same geometric fact that makes the clinical **pivot-shift** test (valgus + internal
   rotation + flexion-extension, combined) more sensitive for ACL deficiency than the anterior-drawer
   test (sagittal alone).

**Consequence:** "a sagittal-plane-only model must fail" is not a statistical accident of one
dataset — it is a structural prediction of the joint's own articular geometry. Section 2 confirms
this is exactly what happens when it is measured.

---

## 2. The falsifier (symmetric, pre-registered before reading the numbers)

**C:** ACL force/strain during a cut scales with the **combined** valgus + anterior-shear +
internal-rotation load; a sagittal-plane-only model **fails** to reach ACL-injury force; hamstring
co-contraction **reduces** anterior shear. **¬C:** sagittal loading alone suffices; frontal/transverse
loads are incidental; hamstrings don't matter. **Adversary forced:** valgus is harmless in isolation,
OR baseline single-timepoint KAM reliably predicts which individual gets hurt (forced via the
largest available pooled test, §4).

| # | sub-claim | verdict | decisive evidence (PMID, live-verified) |
|---|---|:-:|---|
| F1 | Peak KAM in a real cut is reproducible & large | **PASS** | Kristianslund 2013 (23287439): 1.58±0.60 Nm/kg sidestep cut, n=120 elite female handball, **6x** the drop-jump value (0.25±0.16). Mai 2022 (36439622): 1.52-1.73 Nm/kg, n=51. |
| F2 | Combined multiplanar loading > sagittal alone, at the ACL-force level | **PASS** | Markolf 1995 (8544031): cadaveric dial-in loading — anterior-tibial-force + **internal torque** produced the **highest** ACL force in the whole study; + valgus raised force in flexion. Navacchia/Hewett 2019 (30977558): validated FE model, ACL force r=0.99 vs the KAM+internal-rotation-driven contact-force component. |
| F3 | Sagittal-plane loading ALONE cannot reach ACL-injury force (the decisive falsifier) | **PASS** | McLean 2004 (15342155): n=20 subject-specific forward-dynamics models, 5000 stochastic perturbations, pre-set criterion (anterior drawer >2000N) **never** exceeded; valgus loads DID reach injury-capable levels, more often in females. |
| F4 | Hamstring co-contraction reduces ACL loading (protective) | **PASS** | Li 1999 (10213029): cadaveric robotic UFS, n=10 — 80N hamstring co-contraction cuts ACL in-situ force by 30-44% (p<0.05) and anterior tibial translation by 18% (p<0.05). |

All four sub-claims pass on primary, live-verified literature — the falsifier the task specified does
**not** kill the multiplanar mechanism. See `..._evidence.json` `falsifier_tests` for full quoted
abstract text, applied-load values, and honest caveats per test.

---

## 3. The DYSFUNCTION pole — two prospective anchors, one contested replication (symmetric QC)

**Anchor 1 — Hewett TE et al. 2005**, *Am J Sports Med* (PMID 15722287, DOI 10.1177/0363546504269591).
n=205 female athletes, 9 later confirmed ACL ruptures. The 9 injured had **8° greater** knee
abduction angle at landing (p<.05), a **2.5x greater** knee abduction moment (p<.001), 20% higher
GRF (p<.05), 16% shorter stance time. KAM classified injury status at 78% sensitivity / 73%
specificity; a dynamic-valgus composite gave predictive r²=0.88.

> **Headline-precision correction** (this repo's established convention — see
> `docs/MECHANISM_KNEE_LIGAMENTS.md` §0, `docs/MECHANISM_LIGAMENT_STIFFNESS.md` §0 for the same
> discipline applied to other citations): the commonly-repeated "4-6 fold" ACL-injury figure is
> Hewett 2005's own **background** citation of the separately-established female:male injury
> *incidence-rate* ratio — **not** this study's own new finding, which is the **2.5x** KAM
> group-difference above. The two are both real but answer different questions; do not conflate
> them. Also flagged: the frequently-quoted "25.25 Nm" cutoff value is a full-text/Table figure, not
> present in the abstract fetched live this session — not independently re-verified here.

**Anchor 2, decorrelated instrument — Zazulak BT, Hewett TE, et al. 2007**, *Am J Sports Med* (PMID
17468378). n=277 collegiate athletes (140F/137M), 25 later knee injuries. A **trunk** (core)
sudden-force-release perturbation test — a completely different measurement instrument and hidden
state from knee-specific 3D kinematics — found **lateral** (not sagittal) trunk displacement was the
**strongest** predictor of ligament injury (p=.009); a logistic model (trunk displacement +
proprioception + LBP history) predicted knee ligament injury at 91% sensitivity/68% specificity, and
ACL injury in females specifically at 91% accuracy. **Common-mode flag:** shares senior author
(Hewett) with Anchor 1 — a real decorrelation on instrument + cohort, not on research program.

**The adversary, forced to its strongest fair form — Croncurrent A, Creaby MW, Ageberg E 2020**, *BMC
Musculoskelet Disord* (PMID 32819327, **PMC7441716**), a PRISMA systematic review/meta-analysis
pooling **9 independent prospective cohorts**. Neither peak knee-abduction angle nor peak knee-
abduction **moment** significantly predicted future ACL injury (moment: mean diff **−10.61** Nm,
95% CI **−26.73 to 5.50**, n_injured=54/n_controls=1330) — the point estimate trends in the
**opposite** direction from Anchor 1's within-study finding, though the wide CI crosses zero.

**Resolution (symmetric QC, not a hidden contradiction):** the **mechanistic** claim (§2: combined
loading generates more ACL force than sagittal alone, and hamstrings protect against it) is a
physics/cadaver-mechanics fact, robust to population sampling noise — it stands. The **prospective-
classifier** claim (a single baseline KAM measurement in a healthy athlete predicts *which specific
individual* tears their ACL years later) is a much harder population-inference claim, and the
current fairly-pooled evidence does **not** support it. These are different claims; conflating them
is exactly the overclaim this document exists to avoid. (This independently reproduces, via a fresh
PMID chase rather than reuse, the same symmetric-QC pattern this repo's sibling `MSK-CUTTING-ACL-LOAD`
cell already found with a single refuted study — here the refutation is stronger: a 9-study
meta-analysis, not one paper.)

---

## 4. Modulators: sex and fatigue (both push toward the dysfunctional pole)

- **Sex** — Ford, Myer, Hewett 2003, *Med Sci Sports Exerc* (PMID 14523314, n=81 HS basketball
  players): females (n=47) land with significantly greater total valgus motion and greater maximum
  valgus angle than males (n=34) during a drop vertical jump, plus significant side-to-side asymmetry.
- **Fatigue** — McLean et al. 2007, *Med Sci Sports Exerc* (PMID 17473777, n=20 NCAA athletes):
  fatigue increases initial-contact and peak-stance knee abduction/internal-rotation motion **and**
  moments; the abduction-moment fatigue effect is more pronounced in females.

---

## 5. Task-specificity and reactive cutting — the "position of no return" is measurable

- **Besier et al. 2001** (PMID 11445765, n=11): unanticipated cutting → varus/valgus and
  internal/external-rotation moments up to **2x** the preplanned magnitude; sagittal moments
  unchanged.
- **Giesche et al. 2021 meta-analysis**, *Br J Sports Med* (PMID 34344709, **25 trials, 485
  participants**): unplanned tasks → significantly higher external knee **abduction** moment (SMD
  0.34, 95% CI 0.16-0.51) **and** tibial **internal-rotation** moment (SMD 0.51, CI 0.23-0.79); **no**
  significant difference for sagittal mechanics. Effect larger in non-professional athletes.
- **Disclosed minority nuance — Mai et al. 2022** (PMID 36439622, n=51): within complexity-matched
  tasks, true unanticipation (1.64 Nm/kg) did **not** exceed a complex-but-preplanned task (1.73
  Nm/kg) — both clearly exceed the simplest preplanned cut (1.52 Nm/kg). A real, smaller-N exception
  reported honestly against the larger pooled direction, not hidden.
- **Scope boundary — Kristianslund & Krosshaug 2013** (PMID 23287439): cross-task rank correlation
  for KAM itself between drop-jump and sidestep-cut is **poor** (ρ=0.135) despite moderate
  correlation for valgus angle (ρ=0.706) — Hewett 2005's own screening task (a drop-vertical-jump)
  does not reliably rank-order an athlete's actual cutting-task KAM.

**Synthesis:** reactive/unanticipated demand robustly (pooled, 25 trials) elevates frontal +
transverse-plane loading specifically, leaving sagittal loading unchanged — mechanistically coherent
with §2's finding that the dangerous load is multiplanar. This is the quantitative version of "the
position of no return": a late-decision cut is measurably, not just anecdotally, loaded differently
and more dangerously than a rehearsed one.

---

## 6. Coupling to the video pipeline — can this twin flag it from markerless video?

**Honest negative for direct frontal-plane (KAM) precision**, confirmed by 3 independently-converging
sources — this is the load-bearing answer for the operator's stated goal, so it is reported plainly,
not softened:

1. **Turner, Chaaban, Padua 2024**, *J Biomech* (PMID 38905926, n=437 trials): OpenCap markerless
   vs marker-based — sagittal knee/hip CMC>0.94 (excellent), but **frontal** plane CMC only
   **0.47-0.78** and transverse 0.51-0.6.
2. **Fisber, Horsak, Paternoster 2026**, *Sci Rep* (PMID 41876778, PMC13018474, n=24/240 drop
   jumps) — a validation built specifically for ACL re-injury screening: frontal-plane knee
   kinematics RMSE **exceeded 6°** despite a strong waveform correlation (r>0.90 — the same
   high-r/poor-absolute-accuracy pattern this repo's own memory already flags as a generic trap).
   Explicit conclusion: **"OpenCap currently cannot be recommended for ACL re-injury risk
   assessment."** (This citation's unusual 2026 date was cross-checked against 2 decorrelated
   registries — NCBI PubMed and CrossRef — before trusting it; both agree it is real, not
   fabricated.)
3. **This twin's own structural gap** (`VIDEO-LINK1-KINEMATICS`, `bt_memory/video-to-injury-load-
   two-barriers...`): the deployed LaiArnold knee is a 1-DOF `CustomJoint` — ab/adduction is a fixed
   polynomial *function* of flexion angle, not an independent coordinate, so frontal knee moment is
   structurally **absent** from this twin's pipeline output regardless of video accuracy. Separately,
   inverse-dynamics double-differentiation degrades video-derived knee-moment accuracy specifically
   in the first ~100ms post-contact (r drops from 0.79 to 0.46-0.54) — **exactly** the ~40-50ms
   early-stance window this task names as the injury window. The video chain is measurably worst
   exactly where the injury happens, for a temporal reason fully independent of the planar one.

**What IS defensibly video-estimable today:** sagittal knee/hip flexion kinematics and moments
(CMC>0.9, few-degree error) — usable for stance-time, loading-rate, and GRF-proxy screening. A gross
valgus **angle** (not moment) flag is plausible at coarse resolution (CMC 0.47-0.78 is degraded, not
zero) but was **not** directly tested by the sources above and must not be assumed validated.

---

## 7. Couplings

- `docs/MECHANISM_KNEE_LIGAMENTS.md` — this twin's own ACL `Blankevoort1991Ligament` bundles (12
  fibers/side) engage near extension in its own flexion sweep, consistent with Markolf's finding
  that the anterior-force+internal-torque combination peaks near extension/hyperextension.
- `docs/MECHANISM_LIGAMENT_STIFFNESS.md` — ACL structural stiffness 283 N/mm, PASS vs young/healthy
  cadaveric literature (218-286 N/mm) — the stiffness this mechanism's forces load against.
- `docs/MECHANISM_QUAD_HAM_FASCICLE.md` — the biarticular hamstring fascicle-timing cert; couples
  F4's hamstring-protective FORCE finding to this twin's own muscle-fiber-kinematics layer.
- `docs/MECHANISM_JOINT_FORCE_VALIDATION.md`, `docs/MECHANISM_WOBBLING_JOINT_FORCE.md` — the
  OrthoLoad-anchored sagittal-plane knee/hip contact-force chain this mechanism extends into the
  frontal/transverse ACL-specific plane.
- `docs/MECHANISM_POSE_PIPELINE.md` — the video→`.trc` pipeline whose frontal-plane and impact-window
  limits §6 quantifies.
- Graph nodes: `MSK-CUTTING-ACL-LOAD`, `MSK-KNEE-6DOF-JAM-ACL`, `AUTO-DO-NOT-REGRESS-VIDEO-KAM-
  INJURY-RISK-DIR`, `VIDEO-TO-KNEEMOMENT-DROPJUMP-ACL`, `VIDEO-LINK1-KINEMATICS`.
- The FUNCTION↔DYSFUNCTION organizing axis (`bt_memory/function-vs-dysfunction-...`): FUNCTION pole
  = minimized KAM via trunk/hip/hamstring control (Zazulak 2007, Li 1999); DYSFUNCTION pole = valgus
  collapse (Hewett 2005, Ford 2003, McLean 2007).

---

## 8. Honest gaps

1. **No new raw-data solve executed this pass** — no OpenSim/JAM run. This is a literature-forced
   MECHANISM synthesis; status **OPEN**, not PROVEN/ASSUMED. The decisive next step (§9) is already
   designed by the sibling `MSK-KNEE-6DOF-JAM-ACL` cell, not repeated here.
2. Hewett 2005's "25.25 Nm" cutoff and the Zazulak/Hewett shared-senior-author common-mode are both
   flagged in §3 — not independently re-derived / not fully decorrelated.
3. The 3 KAM-magnitude cutting-cohort sources (Kristianslund, Mai, and the sibling cell's
   Arefin/PLOS-0297592 number) span **2** labs (1 lineage confirmed twice via shared senior author +
   1 independent source not personally re-verified from primary text this pass, since its abstract
   does not itself surface an explicit abduction-moment figure). A fully independent, freshly
   verified 3rd cutting-KAM cohort would harden this further.
4. No dataset in this repo's archive combines {raw GRF + kinematics + KAM + manipulation +
   prospective ACL outcome} (per the sibling cell's prior search: Dryad ACL-outcome bot-gated,
   Camargo cutting dataset structurally 1-DOF/no-KAM) — this document's falsifier rests on
   published cross-study literature, not this repo's own re-executed data.
5. Markolf 1995's abstract reports **direction and ranking** of combined-load ACL force (highest for
   anterior-force+internal-torque; unloading for anterior-force+external-torque), not exact
   resultant-force Newton values for every combination — those live in the paper's figures/tables,
   not fetched this pass. No number was invented to fill this gap.
6. §6's "what IS video-estimable" (gross valgus angle as a coarse flag) is an explicitly-flagged
   plausibility statement, not directly tested by the cited sources — do not treat as validated.
7. Sex/fatigue modulators (§4) are single studies here, not meta-analyzed, unlike the anticipation
   question (§5) which had a meta-analysis available.

---

## 9. Proposed next cell

Do **not** build a direct video→KAM→injury-risk regression — both links are individually weak/contested
(video's frontal-plane axis per §6; single-timepoint-KAM's injury-classifier power per §3's forced
adversary). Instead, carrying forward the sibling `MSK-KNEE-6DOF-JAM-ACL` cell's own design:

1. Swap this twin's deployed LaiArnold 1-DOF knee for OpenSim-JAM's Lenhart2015 6-DOF+ligament knee
   (reuses this repo's own already-built COMAK toolchain from the knee-flagship close — see
   `bt_memory/knee-comak-flagship-close-and-fork-build.md`) so ACL **force** becomes a native model
   output, anchored on the mechanistically-solid Markolf/Navacchia combined-loading relationship
   rather than the contested KAM-alone epidemiological one.
2. **Before trusting any cut/jump output**, inject synthetic frontal-plane kinematic noise at the
   measured CMC 0.47-0.78 / RMSE>6° level (Turner 2024, Fisber 2026 — both now independently
   verified in this pass, not assumed) into COMAK's secondary-kinematics input on an existing
   bundled trial, and machine-measure ACL-force error propagation vs a noise-free baseline. Small
   propagation → proceed to a video→JAM/COMAK→ACL-force cutting cell. Large propagation → re-scope
   honestly to a sagittal-plane-dominant load estimate with frontal-plane flagged explicitly
   not-computable-from-video.
3. Explore Zazulak 2007's trunk-lateral-displacement instrument as a **second, video-tractable**
   proxy for the FUNCTION/DYSFUNCTION axis — markerless trunk sway is plausibly estimable from a
   single camera, unlike knee frontal moment, and is not yet explored in this graph.

---

## 10. Files

- `docs/MECHANISM_CUTTING_ACL_MECHANISM_evidence.json` — every citation (PMID/DOI/n/exact quoted
  number), the 4 falsifier tests, the geometric-mechanism derivation, the 2 prospective anchors +
  the forced meta-analytic adversary, sex/fatigue modulators, task-specificity data, the video-
  coupling verdict, couplings, honest gaps, and the proposed next cell — the machine-readable mirror
  of this document.
- Read (unmodified, prior art in this same graph): `data/MECHANISM_ANCHOR_GRAPH.json` nodes
  `MSK-CUTTING-ACL-LOAD`, `MSK-KNEE-6DOF-JAM-ACL`, `AUTO-DO-NOT-REGRESS-VIDEO-KAM-INJURY-RISK-DIR`,
  `VIDEO-TO-KNEEMOMENT-DROPJUMP-ACL`, `VIDEO-LINK1-KINEMATICS`;
  `data/body_twin/agent_outputs/6dof-jam-knee-acl__a21389a02e8c77390.json`;
  `bt_memory/video-to-injury-load-two-barriers-id-double-differentiation-and-frontal-plane-knee-dof.md`;
  `bt_memory/function-vs-dysfunction-build-the-athlete-baseline-in-parallel-with-disease.md`;
  `bt_memory/a-high-correlation-over-a-wide-range-can-be-tail-leveraged...md`;
  `bt_memory/knee-comak-flagship-close-and-fork-build.md`; `docs/MECHANISM_KNEE_LIGAMENTS.md`;
  `docs/MECHANISM_LIGAMENT_STIFFNESS.md`.

ISOLATION (`COORDINATOR.md` §1): bodytwin only; all graph/memory files above read in place, never
mutated; no git commit/push; the 2 files listed at the top are the only new files this session.
