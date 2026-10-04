# MECHANISM — Spinal Cord Injury: Ranking the Neural-Bypass Routes to Restored Movement (2026-07-22)

**Status:** HYPOTHESIS awaiting QC — literature-anchored comparative review, no in-repo simulation run
this session. Scope: rank the real, demonstrated routes to restoring movement after spinal cord injury
(SCI) — (a) intracortical BCI + FES reanimation, (b) epidural/transcutaneous spinal cord stimulation
(SCS), (c) the brain-spine "digital bridge," (d) nerve/axon regeneration + regenerative adjuncts — and
force the honest adversary: real daily-life independence, or small-n lab demonstrations with tethered
equipment. Evidence manifest: `docs/MECHANISM_SCI_NEURAL_BYPASS_evidence.json`. **16 PMIDs — 15 in the
core comparison table (§2) + 1 supporting mechanistic citation (§1) — every one live-verified this
session** via direct NCBI eutils `esearch`→`esummary`→`efetch` (not recalled). Of 14 distinct first-pass
literature searches attempted, **5 returned the WRONG paper** (confirmed by journal/year/title mismatch on
`esummary` — Bouton, Ajiboye, Gill's own erratum-vs-original, and the natural-recovery anchor all
mis-hit on the first try, and a 5th, a suspected "take-home cervical-stimulation SCI" citation, resolved
to a real but off-topic post-stroke paper and was abandoned rather than corrected) and **2 returned zero
hits** (Rowald, the INSPIRE scaffold trial), all re-searched with tighter terms until a title/journal/author
match was confirmed or the citation was dropped. Full account in §8 — a concrete demonstration of this
repo's own measured ~62% memory-citation-drift risk, not swept under the rug.

---

## 0. Pre-registered falsifier (stated before the verdict)

**C0** (claim under test): *current neurotechnology / regenerative approaches restore REAL functional
independence in daily life after SCI* (not merely a monitored lab/clinic demonstration).
**Pre-registered threshold for C0** — met only if a method clears **all four**: (i) explicit
documentation of use **outside** a research session (the source text itself says "home" or "community
setting"); (ii) sustained benefit for **≥3 months**; (iii) the function is performed by the **person's
own body effector** (their own limb/muscle), not a surrogate robotic device; (iv) traceable to a
live-verified peer-reviewed primary source, n≥1.
**¬C0:** fails any one of the four, or is demonstrated only within a monitored lab/clinic visit.

**Adversary being forced (the thing tempting a lazy "yes"):** that a single dramatic press-covered case
(digital bridge, BCI+FES) is mistaken for a scalable, reproducible therapy; and conversely, the void-floor
adversary that would make a lazy "no" too easy — that any observed gain is just **spontaneous natural
recovery**, not the intervention. Forced via Fawcett et al. 2007 (PMID 17179973, ICCP panel, live-verified,
direct quote): *"the spontaneous recovery of motor function in patients with motor-complete SCI is fairly
limited and predictable... the vast majority of recovery occurs in the first 3 months, but a small amount
can persist for up to 18 months."* Every trial below enrolled **chronic** patients (2.5–12+ years
post-injury, all past that 18-month window, several explicitly after **failing** conventional locomotor
training first, e.g. Angeli 2018 NEJM), and several build a **within-subject ON/OFF stimulation
contrast** directly into the case report (Harkema 2011: leg control recovered "but only during epidural
stimulation") — a causal design spontaneous recovery cannot mimic, because natural recovery does not
switch off when a stimulator is switched off. The void-floor adversary is forced to its strongest form
(the paper's own natural-history authority, not a strawman) and **it falls** for (b)/(c); it is not even
in contention for (a)/(d), which don't rely on residual spared tissue at all (§1).
**Anchor (decorrelated, not a tautology):** 9 independent research groups/institutions across 3 countries
supply the full evidence base across all 4 method classes (§2) — never a single self-referential source.
Within that, the 4 groups working on method (b) specifically (Louisville, Mayo, Lausanne, Seattle)
converge qualitatively on the SCS finding (§4), while the regeneration finding sharply diverges internally
— a single-arm scaffold signal (Kim 2022) killed by its own field's fair randomized test (Harrop 2024/2025,
§6) — the anchor is the CONTRAST itself, not a uniform consensus manufactured by citing only agreeing
sources.

---

## 1. The pathway geometry — WHERE each method re-inserts a bypass around the same cut edge

Treat the descending motor command as a directed graph (this repo's own convention, already built in
`docs/MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md` §1 for the intact case):

```
M1 (motor cortex)  --CST (corticospinal tract, white matter)-->  LSN (spinal circuitry: interneurons +
motoneuron pool, e.g. lumbosacral)  --ventral root / peripheral nerve-->  NMJ  -->  MUSCLE
```

**A traumatic SCI lesion severs the CST *edge* at the injury level — it does not destroy the M1 node or
the LSN/muscle nodes**, which are each typically anatomically intact but disconnected (exactly the
UMN-lesion mechanism `MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md` §5 already quantifies via Doherty et al.
2002). Each of the four approaches is a different **edge-insertion point** in this same graph, and the
insertion point determines what tissue must still be spared for the method to work at all — this is the
real, falsifiable, non-heuristic reason the four methods have such different eligibility criteria and
different ceilings:

| Method | New edge inserted | What must be spared for it to work | Bypasses the lesion, or boosts the survivor? |
|---|---|---|---|
| **(a) BCI+FES** | `M1 --[decode+wire]--> MUSCLE` (direct) | Nothing below the cut — decode is upstream, drive is downstream of *everything*, including the LSN itself | Full bypass; works even with **zero** spared CST/LSN tissue (both landmark patients were AIS A) |
| **(b) Epidural/transcutaneous SCS** | *No new edge from M1.* Raises the gain/excitability of the existing `LSN` node itself | Some spared fraction of the CST edge (even sub-clinical) **or** an intact afferent-feedback loop into the LSN | Boosts the survivor — does not reconnect to the brain at all; mechanism confirmed by Angeli et al. 2014 (PMID 24713270, live-verified): voluntary movement recovered even in patients "diagnosed with a motor **and sensory** complete lesion," i.e. neuromodulating sub-threshold excitability, not routing a detectable signal |
| **(c) Digital bridge** | `M1 --[decode+wire]--> LSN` (not muscle) | The `LSN --> NMJ --> MUSCLE` machinery downstream must still be intact and drivable by stimulation (same requirement as (b)) | Reconstructs the missing **edge itself**, restoring continuous cortical modulation over the spinal pattern generator — the topologically closest match to the original intact graph among all four methods |
| **(d) Nerve transfer** (peripheral, real axon regrowth) | New edge at the **peripheral-nerve** tier: reroutes a spare/intact donor nerve — itself still driven by an intact (if underused) UMN→CST→LMN path from a **different, uninjured** segment — into the denervated `MUSCLE`'s NMJ | A physically intact donor peripheral nerve *and* an available denervated target muscle at/near the injury level; never touches the CST edge or the lesion at all | Peripheral-only rewire; does not attempt to regrow an axon **through** the lesion/glial scar — that specific edge (true CNS axon regeneration across the injury site) is the one nobody has yet rebuilt in a human at the function level (§6, INSPIRE 2.0 RCT null) |

**A mechanistic side-note this repo can now cross-check against its own already-verified motor-unit
model:** FES in (a) drives peripheral motor axons **electrically**, and extracellular electrical
stimulation activates the **largest-diameter** (lowest-threshold, fastest, most fatigable) axons **first**
— the exact reverse of the natural Henneman size-principle order (smallest/slowest-fatiguing units first)
this repo already live-verified in `docs/MECHANISM_MOTOR_UNIT.md` (Henneman et al. 1965, PMID 14328454;
Fuglevand et al. 1993, PMID 8120594). Classic reference for the reverse-recruitment mechanism: Gorman PH,
Mortimer JT (1983) *IEEE Trans Biomed Eng* 30(7):407-14, **PMID 6604691** — title/journal/DOI confirmed
live this session, but (like Lawrence & Kuypers 1968 in the corticospinal doc) it is a **pre-abstract-era
record with no fetchable abstract text**; the mechanism itself (cable theory: activation threshold ∝
1/diameter under extracellular stimulation) is classically-established textbook electrophysiology, flagged
here as recalled-not-requoted rather than hidden. This is the real, geometric (not narrated) reason every
FES-based method in this review (Bouton, Ajiboye, Inanici) reports session-duration limits from rapid
muscle fatigue — a genuine, falsifiable consequence of *where in the pathway* the bypass injects current,
distinct from and not yet modeled by `docs/MECHANISM_CENTRAL_FATIGUE.md` (which scopes **voluntary
central** fatigue only, not electrically-imposed peripheral recruitment order — an open coupling, §11).

---

## 2. Citation table — 15 PMIDs, every one live-verified this session (NCBI eutils `esearch`→`esummary`→`efetch`)

| # | Citation | PMID / DOI | n / cohort | Key verified content |
|---|---|---|---|---|
| 1 | Bouton CE, Shaikhouni A, Annetta NV, et al. (2016). Restoring cortical control of functional movement in a human with quadriplegia. *Nature* 533(7602):247-50. | **27074513**, `10.1038/nature17435` | n=1, C5-C6 quadriplegia | Chronic intracortical array (M1) → custom NMES → forearm muscles; **6 wrist/hand motions**; "when using the system," motor level improved C5-C6→C7-T1 unilaterally, enabling grasp/manipulate/release. |
| 2 | Ajiboye AB, Willett FR, Young DR, et al. (2017). Restoration of reaching and grasping movements through brain-controlled muscle stimulation in a person with tetraplegia. *Lancet* 389(10081):1821-30. | **28363483**, `10.1016/S0140-6736(17)30601-3` | n=1, 53yo man, C4, **AIS A** | 2 intracortical arrays + 36 percutaneous FES electrodes + motorized mobile arm support (gravity compensation); point-to-point accuracy **80–100%**; drank coffee **11/12** attempts at 463 days post-implant; self-fed at 717 days. BrainGate2 trial NCT00912041. |
| 3 | Harkema S, Gerasimenko Y, Hodes J, et al. (2011). Effect of epidural stimulation of the lumbosacral spinal cord on voluntary movement, standing, and assisted stepping after motor complete paraplegia. *Lancet* 377(9781):1938-47. | **21601270**, `10.1016/S0140-6736(11)60547-3` | n=1 (Rob Summers), C7-T1, motor-complete | 170 locomotor sessions/26 months pre-implant; 16-electrode L1-S1 epidural array; full weight-bearing stand, balance-assist only, **4.25 min**; supraspinal leg control recovered at 7 months **"only during epidural stimulation."** |
| 4 | Angeli CA, Edgerton VR, Gerasimenko YP, Harkema SJ (2014). Altering spinal cord excitability enables voluntary movements after chronic complete paralysis in humans. *Brain* 137(5):1394-409. | **24713270**, `10.1093/brain/awu038` | n=4 (incl. #3's patient + 3 new) | Voluntary movement immediately after implant in 2 patients diagnosed **motor AND sensory complete**; mechanism = neuromodulating sub-threshold spinal excitability, "key to recovery... in four of four." |
| 5 | Angeli CA, Boakye M, Morton RA, et al. (2018). Recovery of Over-Ground Walking after Chronic Motor Complete Spinal Cord Injury. *N Engl J Med* 379(13):1244-50. | **30247091**, `10.1056/NEJMoa1803588` | n=4, chronic 2.5–3.3 yrs, **after failing locomotor training alone** | 2/4 achieved over-ground (not treadmill) walking: **278 sessions/85 weeks** and **81 sessions/15 weeks**; all 4 gained independent standing + trunk stability; 1 hip fracture adverse event. |
| 6 | Gill ML, Grahn PJ, Calvert JS, et al. (2018). Neuromodulation of lumbosacral spinal networks enables independent stepping after complete paraplegia. *Nat Med* 24(11):1677-82. | **30250140**, `10.1038/s41591-018-0175-7` (original; **30353100** is a Publisher Correction, disclosed §8) | n=1 (Jeff Marquis), complete sensorimotor paralysis | 43-week multimodal rehab + EES → independent treadmill stepping (no trainer/BWS); independent **over-ground** stepping with a front-wheeled walker, **trainer assistance at the hips for balance**. |
| 7 | Wagner FB, Mignardot JB, Le Goff-Mignardot CG, et al. (2018). Targeted neurotechnology restores walking in humans with spinal cord injury. *Nature* 563(7729):65-71. | **30382197**, `10.1038/s41586-018-0649-2` | n not stated in abstract (honest gap, §8) | >4 yrs post-injury; spatiotemporal epidural stim, real-time triggered; within 1 week, re-established overground walking; after months, **voluntary control without stimulation** regained + walk/cycle in "ecological settings" **during** stimulation. |
| 8 | Rowald A, Komi S, Demesmaeker R, et al. (2022). Activity-dependent spinal cord neuromodulation rapidly restores trunk and leg motor functions after complete paralysis. *Nat Med* 28(2):260-71. | **35132264**, `10.1038/s41591-021-01663-5` | **n=3**, complete sensorimotor paralysis | New paddle-lead design; within **a single day**: stand, walk, cycle, swim, control trunk; rehab "restore these activities in **community settings**." |
| 9 | Lorach H, Galvez A, Spagnolo V, et al. (2023). Walking naturally after spinal cord injury using a brain-spine interface. *Nature* 618(7963):126-33. | **37225984**, `10.1038/s41586-023-06094-5` | n=1 (Gert-Jan Oskam), chronic tetraplegia | Fully implanted cortex+spine wireless "digital bridge"; calibrated in minutes; stable **>1 year, "including during independent use at home"**; stand/walk/stairs/complex terrain; regained crutch-assisted overground walking **even with BSI switched off**. |
| 10 | Inanici F, Brighton LN, Samejima S, Hofstetter CP, Moritz CT (2021). Transcutaneous Spinal Cord Stimulation Restores Hand and Arm Function After Spinal Cord Injury. *IEEE Trans Neural Syst Rehabil Eng* 29:310-9. | **33400652**, `10.1109/TNSRE.2021.3049133` | n=6, complete paralysis, up to 12 yrs post-injury | **Non-invasive**; rapid + sustained hand/arm recovery, "matched or exceeded" implanted-stimulation results (comparative claim, not head-to-head); gains maintained **3–6 months after stimulation stopped** in all 6; resumed guitar/painting. |
| 11 | van Zyl N, Hill B, Cooper C, Hahn J, Galea MP (2019). Expanding traditional tendon-based techniques with nerve transfers for the restoration of upper limb function in tetraplegia. *Lancet* 394(10198):565-75. | **31280969**, `10.1016/S0140-6736(19)31143-2` | **n=16 (27 limbs)**, early (<18mo) cervical SCI, C5 and below | 59 nerve transfers; 24-month ARAT 16.5→**34.0** (p<0.0001), GRT 35.0→**125.2** (p<0.0001); grasp strength 2.8–3.9 kg; 6 minor AEs, no lasting consequence. |
| 12 | Curtis E, Martin JR, Gabel B, et al. (2018). A First-in-Human, Phase I Study of Neural Stem Cell Transplantation for Chronic Spinal Cord Injury. *Cell Stem Cell* 22(6):941-50. | **29859175**, `10.1016/j.stem.2018.05.014` | n=4, T2-T12 chronic SCI | NSI-566 neural stem cells, 6 stereotactic injections; no serious AEs 18-27mo; **2/4 gained 1-2 ISNCSCI levels**; authors' own words: "lacks statistical power or a control group." |
| 13 | Kim KD, Lee KS, Coric D, et al. (2022). Acute Implantation of a Bioresorbable Polymer Scaffold... 24-Month Follow-up (INSPIRE). *Neurosurgery* 90(6):668-75. | **35442254**, `10.1227/neu.0000000000001932` | n=19 implanted (16 evaluable), single-arm | **7/16 (44%) improved ≥1 AIS grade at 6 months** vs. informal "historical controls." |
| 14 | Harrop JS, Kim KD, Okonkwo DO, et al. (2024/2025). ...A Randomized Controlled Trial (INSPIRE 2.0). *Neurosurgery* 96(4):751-62. | **39471088**, `10.1227/neu.0000000000003180` | **n=20 (10 scaffold / 10 surgery-alone control), randomized, blinded** | At 6 months: scaffold **2/10 (20%)** converted vs. control **3/10 (30%)** — scaffold numerically *worse*; pre-registered success bar (≥20-point advantage) missed; **"did not produce probable clinical benefit."** Study closed. |
| 15 | Fawcett JW, Curt A, Steeves JD, et al. (2007). ICCP guidelines: spontaneous recovery after SCI and statistical power for clinical trials. *Spinal Cord* 45(3):190-205. | **17179973**, `10.1038/sj.sc.3102007` | ICCP panel, Sygen-trial reanalysis | Void-floor anchor (§0): motor-complete spontaneous recovery "fairly limited and predictable"; "vast majority... first 3 months," tail to 18 months. |

---

## 3. Method (a) — intracortical BCI + FES reanimation

**C_a:** cortical decode + FES restores real voluntary control of the person's own paralyzed muscles for
functional tasks. **Threshold:** ≥1 functional ADL-relevant task completed using the participant's own
limb, at above-chance accuracy, in a live-verified primary source. **Cleared** — Bouton 2016 (#1) and
Ajiboye 2017 (#2) both clear it decisively: 6 wrist/hand motions and 80–100% point-to-point accuracy,
self-feeding and drinking coffee using the participant's own reanimated arm, not a robotic surrogate.

**Forced adversary (real independence vs. lab-demo, its strongest form):** both are **n=1**, and both
systems are structurally incapable of home use as built: Ajiboye's rig needs 36 **percutaneous** leads
plus an external **motorized mobile arm support** for gravity compensation (the BCI does not even control
the shoulder's own weight-bearing); Bouton's array fed a benchtop NMES controller, and the reported
functional gain is explicitly qualified **"when using the system"** — not a lasting off-system recovery
the way Wagner 2018/Lorach 2023 report for method (b)/(c). Neither paper's abstract, nor any live search
this session, surfaced a take-home/community-use claim for either system (§8). **Adversary does NOT
fall** — it is the correct description of the current state of the art for this specific method.
**Verdict: C_a TRUE at the task level (real, own-limb, above-chance function), FALSE at the independence
level** — the honest ceiling today is in-lab task restoration under an operator's supervision, full stop.

---

## 4. Method (b) — epidural / transcutaneous spinal cord stimulation

**C_b:** SCS restores standing/stepping/walking in **chronic, motor-complete** SCI beyond the
natural-recovery void floor. **Threshold:** documented in patients >18 months post-injury (past
Fawcett 2007's own recovery window, #15), motor-complete (AIS A/B) status established by standard exam,
with a causal (not merely correlational) link to the stimulation.

| Study | n | Time post-injury | Result | Independence marker |
|---|---:|---|---|---|
| Harkema 2011 (#3) | 1 | chronic | Standing 4.25 min | Leg control present **only stim-ON** (built-in causal contrast) |
| Angeli 2014 (#4) | 4 | chronic | Voluntary movement immediate post-implant, 2/4 sensorimotor-complete | Mechanism = sub-threshold excitability, no signal routing needed |
| Angeli 2018 NEJM (#5) | 4 | 2.5–3.3 yrs, **post-training-failure** | 2/4 over-ground walking after 81–278 sessions | Still needs walker/gait-training support |
| Gill 2018 (#6) | 1 | chronic | Independent treadmill stepping (43 wks); over-ground with walker | **Trainer assistance at hips** for balance even at best |
| Wagner 2018 (#7) | n/s (honest gap) | >4 yrs | Overground walking in 1 week; voluntary stim-OFF control after months | Ecological-setting walk/cycle still **during** stimulation |
| Rowald 2022 (#8) | 3 | chronic | Stand/walk/cycle/swim/trunk control in 1 day | "**Community settings**" — the strongest real-world language in this tier |
| Inanici 2021, transcutaneous (#10) | 6 | up to 12 yrs | Hand/arm function restored, non-invasive | Gains **outlast** stimulation by 3–6 months; resumed hobbies |

**Forced adversary — the void-floor, in its strongest form (§0):** every cohort is chronic, several
explicitly post-training-failure, and 3 of the 7 rows above build an explicit ON/OFF or before/after
contrast that a spontaneous-recovery explanation cannot produce. **The adversary falls** across this
whole, methodologically diverse set (4 independent groups: Louisville/Angeli-Harkema, Mayo/Gill,
Lausanne/Courtine-Wagner-Rowald, Seattle/Moritz-Inanici — implanted **and** non-invasive, 2 countries).
**Verdict: C_b TRUE and the best-replicated finding in this entire review** — but independence is
**partial**: every single row still required an implanted or surface stimulator, and all but the
non-invasive row needed a walker or trainer-assisted balance. "Community settings" (Rowald 2022) is real
progress, not full autonomy from the hardware.

---

## 5. Method (c) — the brain-spine "digital bridge"

**C_c:** a wireless cortex→spinal-stimulation link restores natural walking with real independent,
at-home use. **Threshold:** explicit "home" or "community" language in the primary source, ≥3 months
durability. **Cleared, uniquely strongly** — Lorach 2023 (#9) states directly: reliability "remained
stable over one year, **including during independent use at home**," natural gait over "stairs" and
"complex terrains," **and** partial device-independent carryover (crutch-assisted overground walking with
the BSI **switched off**) — a stronger, more explicit independence claim than any other single source in
this review.

**Forced adversary:** n=1, two fully implanted bespoke systems (cortical recording array + spinal
stimulator + wireless link), a single surgical/engineering team (EPFL/CHUV/Clinatec), and **no
independent replication** was found live this session. **Adversary partially falls** (the specific home-use
claim is textually explicit and decisive for this one participant) but **does not fall** on
generalizability (n=1, unreplicated, extreme cost/complexity). **Verdict: C_c TRUE for n=1, unresolved for
n>1** — the single most complete "real independence" demonstration among all implanted-neurotechnology
methods, honestly bounded to a proof-of-concept scale.

---

## 6. Method (d) — nerve/axon regeneration + regenerative adjuncts

Three sub-routes at very different maturities, forced against each other, not averaged together:

**(d1) Peripheral nerve transfer surgery — van Zyl 2019 (#11), the strongest REAL-independence result in
this entire review.** n=16 (27 limbs), 59 nerve transfers, 24-month durable, statistically **and**
clinically significant (ARAT/GRT/SCIM beyond MDC/MCID), and — critically for C0's threshold (iii) — the
restored grasp is produced by the **person's own reinnervated muscle**, with **zero external hardware,
battery, or wire** after healing. This clears all four C0 gates cleanly (home use is simply "using your
hand," durability confirmed to 24 months, n=16 traceable to a live-verified Lancet paper). **Scope limit,
disclosed plainly, not hidden:** only applies to upper-limb function in **early** (<18 months),
**cervical**, motor-level-C5-or-below injuries with an available donor nerve — it does not touch gait,
does not apply to thoracic/lumbar paraplegia, and (§1) does not reconstruct anything across the actual
CNS lesion.

**(d2) Neural stem cell transplantation — Curtis 2018 (#12), honestly weak.** n=4, Phase I safety trial,
2/4 patients gained 1-2 neurological levels — but the authors themselves disclose it "lacks statistical
power or a control group needed to evaluate functional changes." No forcing needed here: the source
forces its own honesty.

**(d3) Bioresorbable polymer scaffold (INSPIRE) — the cleanest forced-adversary KILL in this whole
review.** The single-arm study (Kim 2022, #13) reported a promising 44% AIS-conversion signal against
informal historical controls. The properly randomized, blinded follow-up (Harrop 2024/2025, #14, **the
fair adversary this signal owed itself**) tested the *same* device against the *correct* control (spine
surgery alone) at the *same* 6-month endpoint: scaffold 20% vs. control 30% — the scaffold arm was
numerically **worse**, the pre-registered ≥20-point benefit bar was missed, and the trial was closed. This
is a textbook symmetric-QC result: a real, published, single-arm positive signal, killed by the exact fair
comparator it was always going to need, at the same measurement the promising result was originally
claimed on. **True axon regeneration across the injury site itself remains unproven at the human function
level as of this review.**

**Verdict on (d), unaveraged:** nerve transfer (peripheral, real, durable, home-independent, scope-limited)
is genuinely excellent; stem cells are an honest early safety signal; scaffold-mediated intra-lesion
regeneration failed its own fair test. Do not report "regenerative medicine" as one number — the three
sub-routes sit at three entirely different maturities.

---

## 7. Ranked comparison — real functional restoration, independence-weighted

| Rank | Method | n | Injury | Durability | Home/lab | What actually returned | Own-body or hardware? |
|---|---|---|---|---|---|---|---|
| **1** | Nerve transfer surgery (#11) | 16 (27 limbs) | Cervical, C5↓, early | 24 mo, durable | **Home** (surgical, permanent) | Grasp/pinch/elbow extension, own hand | **Own body**, zero hardware |
| **2** | Epidural/transcutaneous SCS (#3–8, #10) | 1–6 per study, 4 groups | Chronic motor-complete, thoracic/cervical | Months–years | Mostly lab/clinic; "community settings" in #8 | Standing, stepping, overground walking, hand/arm (transcutaneous) | Hardware-dependent (implant or surface electrodes) |
| **3** | Digital bridge (#9) | 1 | Chronic tetraplegia | **>1 yr**, incl. explicit home use | **Home**, uniquely explicit | Natural gait, stairs, terrain; partial stim-off carryover | Hardware-dependent (2 implants), one own-body carryover result |
| **4** | Intracortical BCI + FES (#1, #2) | 1, 1 | Cervical, incl. AIS A | Session-based only | **Lab only**, no home evidence found | 6-motion grasp, drink, self-feed via own arm | Own limb driven by tethered/benchtop hardware |
| **5** | Regenerative adjuncts (#12–14) | 4 (stem cell), 19/20 (scaffold) | Thoracic, chronic/acute | Safety-tier only; scaffold RCT null | Lab/surgical only | Modest neuro-level gains (stem cell, uncontrolled); **no benefit** (scaffold RCT) | Intended own-body, unproven at function level |

---

## 8. Overall verdict — the true current ceiling, stated plainly

**No approach today gives a person with complete, chronic SCI unassisted, hardware-free, normal walking
or hand function at population scale, on demand.** What is real, forced, and machine-checkable from this
session's live-verified sources:

1. **Nerve-transfer surgery** restores genuinely independent, permanent, hardware-free upper-limb function
   — but only for appropriate cervical candidates, upper limb only, and it never touches gait.
2. **Epidural/transcutaneous spinal cord stimulation** is the **best-replicated** real effect (4 independent
   groups, implanted and non-invasive both converging) — standing and stepping/walking restored well
   beyond the natural-recovery window — but it remains implant/electrode-and-usually-walker-dependent for
   nearly every patient reported.
3. **The digital bridge** delivered the single most complete real-world independence demonstration in the
   literature (explicit >1-year home use) — for **one** person, unreplicated, at enormous
   bespoke cost/complexity.
4. **Intracortical BCI+FES** restores real, richly dexterous, cortically-driven control of a person's own
   paralyzed limb for specific tasks — but has **never**, in either landmark trial, been shown outside a
   monitored research session.
5. **Regenerative medicine** for the spinal cord lesion itself is not yet a functional-restoration method
   — the one rigorously randomized trial in this space found no benefit over standard surgical care.

**Honest gap disclosure, this session's own citation-drift catches (methodologically load-bearing, not
housekeeping):** my first-pass title/author searches for Bouton 2016, Ajiboye 2017, Gill 2018 (original vs.
its own Publisher Correction), and the natural-recovery void-floor anchor all initially returned a
**wrong** PMID (confirmed by `esummary` journal/year/author mismatch) — the first "Bouton 2016" hit was a
2016 IEEE EMBS conference abstract by a co-author, not the Nature paper; the first "Ajiboye 2017" hit was
an unrelated 2020 *Scientific Reports* decoding paper; the first "spontaneous recovery" hit was Dalkilic et
al.'s CSF/MRI biomarker paper, not the intended ICCP guidelines — each was re-searched with tighter terms
until the title matched the primary source exactly, landing on Fawcett et al. 2007 (#15) for the
natural-recovery anchor. Rowald 2022 and the INSPIRE scaffold trial both returned **zero hits** on an
over-specific exact-title-phrase search and needed a looser author+keyword query. A search for a suspected
"Powell et al., take-home cervical epidural stimulation for SCI hand function" paper resolved instead to
Powell et al. 2023 (PMID 36807682) — a real paper, but about **post-stroke**, not SCI, upper-limb paresis;
it is excluded here as off-topic rather than silently substituted, and Inanici 2021 (#10) stands as this
review's non-invasive/take-home-adjacent SCS citation instead. Wagner 2018's participant count is not
stated in its own abstract; the commonly-reported n=3 (STIMO trial) is disclosed as recalled, not
independently reconfirmed via primary text this session. No epidemiological/prevalence figure is included
— a deliberate scoping choice (the task asked for method comparison, not population burden), not a
discovered-then-hidden gap.

---

## 9. Function↔dysfunction summary (closing the requested axis)

The intact pathway is `M1 →(CST)→ LSN →(peripheral nerve)→ MUSCLE`, one continuous, over-determined
circuit. SCI is a **local edge cut**, not a node failure — cortex and muscle are each typically still
alive and excitable on their own. The four methods are four different, mutually non-exclusive edge-repair
strategies (§1's table), which is exactly why they are not really in competition with each other clinically
— a nerve-transfer candidate can also receive epidural stimulation for gait, and a digital-bridge candidate
is, by construction, also an SCS candidate. The one edge nobody has yet rebuilt directly — true axon
regrowth **through** the lesion itself, reconnecting M1 to LSN biologically rather than electronically —
is exactly the one sub-route (scaffold-mediated intraparenchymal regeneration) that failed its fair,
randomized test this session (#14).

---

## 10. Honest gaps (consolidated)

1. Wagner 2018 (#7) participant count not confirmed via its own abstract text this session (commonly
   reported n=3 elsewhere) — disclosed, not asserted as verified.
2. No take-home/community-use BCI+FES (method a) trial was found live this session for either landmark
   system — the absence itself is the finding (§3), but the search is not exhaustive of all grey literature.
3. Curtis 2018 (#12) is explicitly underpowered by its own authors' statement — used here only at the
   confidence tier its authors themselves assign it.
4. INSPIRE 2.0 (#14)'s null result is itself from a small trial (n=20, 10/arm) — a symmetric caveat: the
   null could in principle be under-powered too, though the point estimate ran in the *unfavorable*
   direction for the scaffold, which is the conservative direction for this particular adversary.
5. Gorman & Mortimer 1983 (PMID 6604691, §1) is a pre-abstract-era PubMed record — title/journal/DOI
   confirmed live, no abstract text exists to quote; the reverse-recruitment mechanism it is cited for is
   classically-established electrophysiology, flagged as recalled-not-requoted, matching this repo's own
   precedent for Lawrence & Kuypers 1968 in `MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md`.
6. Inanici 2021 (#10)'s "matched or exceeded" comparison to implanted stimulation is the source paper's own
   comparative textual claim, not an independently-run head-to-head trial against methods (b)'s implanted
   arm — reported here as their claim, not re-verified as a controlled comparison.
7. No in-repo simulation, OpenSim cell, or live-computed number was produced this session — every figure
   above is a literature citation, machine-cross-checked against the source abstract text, not a
   twin-computed result. Tier: literature-anchored HYPOTHESIS, not a MEASURED cell.
8. Global/national SCI prevalence and incidence were not fetched this session (scoping choice, §8) — the
   existing `WAVE_PLAN.md` Wave 20 entry already points to the NSCISC open SCI registry for that axis if a
   future cell needs it.

---

## 11. Proposed cell

`scripts/mechanism/sci_bypass_pathway_channel_probe.py` (not yet built): formalize §1's edge-insertion-point
argument as an explicit information/gain-budget model per method — cortical decode channel capacity at the
BCI read-out (measurable against **DANDI:000147**, the paralyzed-human intracortical BCI finger-decoding
dataset already catalogued in `WAVE_PLAN.md` Wave 20 — a real public dataset, not yet pulled into a cell
this session), spinal-circuit excitability gain at the SCS insertion point, and peripheral axon count
available at the nerve-transfer insertion point — then check which tier is the binding bottleneck for a
GIVEN injury level, the same three-tier channel-bottleneck logic `docs/MECHANISM_BIONIC_HAND_ARCHITECTURE.md`
§4 already used for prosthetic-hand control (decode-channel, not actuator-DOF, as the binding constraint) —
predicting, before any new case is attempted, which of the four bypass points is the best fit for a given
patient's specific spared-tissue profile, checked against real post-intervention outcome scores rather than
assessed post hoc.

---

## 12. Couples to

- `docs/MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md` — this doc's §1 pathway-geometry table is built directly
  on that document's already-verified M1 population code, corticospinal conduction-delay physics, and the
  Doherty et al. 2002 (n=156) UMN-vs-LMN anatomical-gradient table by SCI level — the healthy/dysfunction
  baseline this document's four bypass methods are all trying to reconstruct or route around.
- `docs/MECHANISM_MOTOR_UNIT.md` — the Henneman size-principle recruitment order (Henneman 1965, PMID
  14328454; Fuglevand 1993, PMID 8120594) this document's §1 FES-reverse-recruitment mechanism directly
  contrasts against; not re-derived here.
- `docs/MECHANISM_CENTRAL_FATIGUE.md` — scopes **voluntary central** fatigue only; this document's
  electrically-imposed peripheral/reverse-recruitment-order fatigue mechanism (§1) is a distinct,
  not-yet-modeled mechanism, an open coupling rather than a duplicate.
- `docs/MECHANISM_BIONIC_HAND_ARCHITECTURE.md` — the decode-channel-vs-actuator-DOF bottleneck argument this
  document's §11 proposed cell directly mirrors for the BCI+FES pathway specifically.
- `WAVE_PLAN.md` Wave 20 (`PARALYSIS-STROKE-STAGE-LOCALIZATION` node, status ASSUMED) — a related but
  **distinct** scope (localizing sensor/relay/central failure stage from neural data) that already points
  to **DANDI:000147** (intracortical BCI data) and the **NSCISC** open SCI registry as data sources this
  document's proposed cell (§11) would consume, not duplicate.
