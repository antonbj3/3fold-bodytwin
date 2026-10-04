# MECHANISM CENTRAL FATIGUE — the CNS/voluntary-activation force limit (2026-07-22)

Builds the **CENTRAL** delta this twin is explicitly missing. `docs/MECHANISM_MUSCLE_FATIGUE.md` (this
repo, same wave) already states its own scope limit verbatim: *"No central fatigue, metabolic/pH, or
motor-unit size-principle recruitment order is modeled — this is the peripheral, lumped-compartment
phenomenological layer only."* `MSK-MUSCLE-FATIGUE-EMG` (`docs/MECHANISM_STATE.md`) is a third, still
different capability — a surface-EMG spectral-shift (MDF-down/RMS-up) fatigue *signature* detector, not a
voluntary-drive measurement. This document is the fourth, non-overlapping layer: **how much of the force
lost during fatigue is because the CNS stops driving the muscle fully**, not because the muscle itself
can no longer respond — quantified by the interpolated-twitch technique, and shown to be, in part, an
**actively regulated** protective limit, not passive failure. Literature-forced synthesis, status **OPEN**,
`mechanism_grade: HYPOTHESIS-B`, awaiting independent QC. Every number below traces to a PMID/DOI fetched
**live** this session via NCBI eutils (esearch → esummary → efetch); one citation's own middle initial
("Merton PB") was caught and corrected to the true "Merton PA" only by this live verification — a direct,
disclosed instance of this project's own previously-measured ~62% citation-drift-from-recall rate. Full
falsifier ledger, every citation, and both quantitative models: `docs/MECHANISM_CENTRAL_FATIGUE_evidence.json`.

---

## 0. Headline

**7 falsifier tests: 5 PASS, 1 PASS-as-MIXED (forced, not swept aside), 1 honest OPEN gap.** Voluntary
activation (VA), measured by twitch interpolation, is a real, reproducible, subject-specific quantity
(Allen/Gandevia/McKenzie 1995) that **falls with fatigue** — 90.7% from >99% over a 3-min sustained
maximal contraction (Gandevia et al. 1996) — and this is **not** a peripheral-muscle-property artifact
riding along: a transcranial-magnetic-stimulation (TMS) superimposed-twitch increment grows ~10× with
fatigue (1.0%→9.8%) and, decisively, **does not recover under post-fatigue ischaemia while a decorrelated
cortical-excitability signal recorded in the identical ischaemia does recover** — a forced double
dissociation a shared peripheral confound cannot produce. Separately, and more surprisingly, group III/IV
muscle afferents **actively cap** how much peripheral fatigue is allowed to develop: pharmacologically
blocking them (lumbar intrathecal fentanyl, sparing efferent motor function) **increases** central motor
drive (+9 to +21% EMG/iEMG across three independent trials) **and increases** peripheral quadriceps
fatigue by a wide, reproducible margin (twitch-force decline −44 to −52% vs. −31 to −34% with intact
afferents) — confirmed on a **decorrelated** measurement axis (direct muscle-biopsy metabolite chemistry,
not just EMG/force) that itself correlates with the twitch deficit (r²=0.79–0.92). The naive "it's just
opioid disinhibition, athletes just push harder" adversary is forced onto its own strongest ground and
still gets the **sign wrong**: blockade shortened, not lengthened, time-to-exhaustion (6.8 vs. 8.7 min) —
losing the afferent brake made overall performance **worse**, not better, because the same afferents also
support adequate ventilation/circulation. At the DYSFUNCTION pole, chronic fatigue syndrome (CFS) shows a
clean, axis-orthogonal central-only deficit (voluntary activation down, every peripheral marker normal),
while multiple sclerosis (MS) forces a genuinely **mixed** verdict — a dramatic central (twitch-occlusion)
deficit (normal ≈100% vs. "rarely >60%" in MS) **plus** a real peripheral/excitation-contraction-coupling
abnormality demonstrated with voluntary drive bypassed entirely (electrically-stimulated contractions) —
reported as PARTIAL, not oversold as "purely central." Overtraining syndrome is the honest, OODA-forced
gap: no dedicated twitch-interpolation/VA study in diagnosed overtraining was located after two targeted
search angles; the field's own 2013 consensus statement states no validated marker exists at all.

| # | falsifier | pre-registered threshold | result |
|---|---|---|---|
| F1 | VA is a real, reproducible signal that falls with fatigue (not noise/motivation) | ≥1 independent reliability study (between-subject consistency, no antagonist-coactivation confound) + ≥2 fatigue studies showing significant VA decline from a high baseline | **PASS** (90.7% from >99%; quadriceps EMG-force dissociation in 6/9 subjects) |
| F2 | the deficit has a genuine supraspinal/cortical (not purely peripheral-artifact) component | a TMS-specific superimposed-twitch increment must grow with fatigue AND fail to recover under an ischaemia block that a decorrelated cortical-excitability signal, recorded simultaneously, DOES recover from — a forced double dissociation | **PASS** (1.0%→9.8%, does not recover in ischaemia; MEP excitability changes do recover; independently replicated as a spinal-level time-course dissociation, 30s vs 90s) |
| F3 | group III/IV muscle afferents actively regulate (not merely correlate with) central drive and peripheral-fatigue severity | pharmacological blockade (resting-cardioventilatory-control-verified) must ↑ central drive AND ↑ peripheral fatigue across ≥2 protocols, confirmed on a decorrelated (non-EMG/force) biopsy/metabolite leg, correlated with the twitch deficit | **PASS**, wide margins (iEMG +12–21%; ΔQtw −44 to −52% vs −31 to −34%; metabolite perturbation +27–129% greater, r²=0.79–0.92) |
| F4 | central:peripheral apportionment is task-dependent, not a fixed universal ratio | ≥2 decorrelated modalities/muscle groups must show materially different apportionment | **PASS** (running > cycling/skiing at matched low-damage loads; ~25% quantified for sustained elbow-flexor MVC specifically, not claimed universal) |
| F5 | CFS shows a clean, axis-orthogonal CENTRAL-only deficit | VA significantly lower in patients vs controls WHILE peripheral markers (metabolism, M-wave, twitch tension) are statistically indistinguishable | **PASS** (Kent-Braun et al. 1993, authors' own conclusion: "important central component") |
| F6 | MS central-vs-peripheral verdict, forced against a pure-central adversary | an electrically-stimulated (voluntary-drive-bypassing) protocol must still show excess fatigue/metabolic perturbation in MS vs controls for the "purely central" adversary to fail | **PASS as MIXED** (twitch-occlusion: ~100% normal vs "rarely >60%" MS; but electrically-stimulated 9-min protocol independently shows excess intramuscular fatigue too — both axes real) |
| F7 | overtraining syndrome central-fatigue link | a direct twitch-interpolation/VA study in diagnosed OTS | **OPEN / honest gap** (2 targeted search angles, none found; field's own 2013 consensus states no validated marker exists at all; one decorrelated autonomic/chemoreflex marker is suggestive, not diagnostic) |

---

## 1. Geometric structure — two distinct mechanisms, not one heuristic

**(A) Twitch interpolation is a null-space projection, not an arbitrary index.** The achievable force of a
muscle at any instant is a point on a manifold set by which motor units are recruited and how fast they
fire (the same recruitment-threshold/rate-coding geometry `docs/MECHANISM_MOTOR_UNIT.md` builds explicitly
for soleus). A supramaximal electrical shock delivered *during* a voluntary effort can only add force by
exciting whatever capacity voluntary drive left **unused** — unrecruited motor units, or units firing below
fusion. If voluntary drive already occupies the full achievable manifold (VA = 100%), the shock has
**nothing orthogonal to add** and the superimposed twitch is exactly zero; if drive is submaximal, the
shock's evoked increment is a direct, quantitative readout of the **unused (null-space) dimension** of the
recruitment/rate-coding space. `VA% = [1 − (superimposed twitch / potentiated resting twitch)] × 100` is
the field-standard formula built on exactly this logic (used throughout the twitch-interpolation literature
cited below, e.g. Allen/Gandevia/McKenzie 1995 §2, Gandevia 2001 §2) — an honest gap, not asserted as
independently re-derived from a live-quoted formula this session (§12).

**(B) Group III/IV afferent feedback is a measured negative-feedback loop clamping a state-space
trajectory, not a metaphor.** Treat (central motor drive, accumulated peripheral fatigue) as a
two-dimensional state. Amann & Dempsey (2008, §5) show that **prefatiguing** a subject to two different
degrees before an identical subsequent time trial produces two different *paths* through this space but
the **same convergent endpoint** — end-exercise peripheral quadriceps-twitch decline of −35% to −37%
regardless of the pre-existing fatigue level (P = 0.35, no difference) — while central drive itself is
measurably suppressed in dose-dependent proportion (−23% with the higher pre-fatigue). That is a genuine
control-theoretic signature: an **invariant ceiling reached from different starting trajectories**, which a
passive, unregulated muscle-fatigability limit does not predict (a passive limit would simply track
cumulative peripheral load, not converge). Cutting the feedback loop (lumbar intrathecal fentanyl, §6)
should — and does — remove the clamp: the ceiling rises (to −44% to −52%) and central drive is measurably
disinhibited (+9 to +21%), in three independent replications, one of them using a completely different
observable (intramuscular metabolite chemistry from a muscle biopsy, not EMG or force at all). This is a
real, falsifiable regulator, derived from the geometry of a feedback loop's own defining signature
(convergence-then-overshoot-when-cut), not a rhetorical label glued onto a correlation.

---

## 2. Citations — 23 unique PMIDs, every one verified LIVE this session (NCBI eutils esearch→esummary→efetch)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Merton PA (1954). Voluntary strength and fatigue. *J Physiol* 123(3):553-64. | **13152698**, `10.1113/jphysiol.1954.sp005070`, PMC1366225 | **Founding twitch-interpolation paper.** Bibliographic record confirmed live (pre-abstract era, no abstract text in PubMed; PMC full text publisher-restricted, same restriction seen on other pre-1975/1990s records in this repo's sibling docs). **Self-correction, disclosed**: this session's first 4 esearch attempts under the recalled author "Merton PB" returned **zero** hits; only a volume+page+journal-field search (bypassing the author field entirely) surfaced PMID 13152698, whose esummary shows the true author string is **"MERTON PA."** A direct, live-caught instance of this project's own previously-measured ~62% citation-drift-from-memory rate — corrected here, not silently passed through. |
| 2 | Gandevia SC (2001). Spinal and supraspinal factors in human muscle fatigue. *Physiol Rev* 81(4):1725-89. | **11581501**, `10.1152/physrev.2001.81.4.1725` | **Master review**, full abstract fetched/quoted. States plainly: "voluntary activation of human motoneurons and muscle fibers is suboptimal," twitch interpolation is the revealing technique, VA "usually diminishes during maximal voluntary isometric tasks," TMS reveals cortical and "supraspinal drive" changes, spinal-level changes arise partly from "group III and IV muscle afferents." |
| 3 | Gandevia SC, Allen GM, Butler JE, Taylor JL (1996). Supraspinal factors in human muscle fatigue: evidence for suboptimal output from the motor cortex. *J Physiol* 490(Pt2):529-36. | **8821149**, `10.1113/jphysiol.1996.sp021164`, PMC1158689 | **F1/F2 decisive quantitative + ischaemia double-dissociation.** VA of biceps: >99%→90.7% over a 3-min sustained MVC; TMS-evoked force increment 1.0%→9.8% over a 2-min MVC; post-fatigue ischaemia does NOT restore voluntary force/VA but DOES let TMS-evoked EMG-excitability changes recover — the decisive dissociation (§4). |
| 4 | Belanger AY, McComas AJ (1981). Extent of motor unit activation during effort. *J Appl Physiol Respir Environ Exerc Physiol* 51(5):1131-5. | **7298453**, `10.1152/jappl.1981.51.5.1131` | Classic VA-quantification methodology paper (no abstract text retrievable this session — title/journal/DOI confirmed live). |
| 5 | Allen GM, Gandevia SC, McKenzie DK (1995). Reliability of measurements of muscle strength and voluntary activation using twitch interpolation. *Muscle Nerve* 18(6):593-600. | **7753121**, `10.1002/mus.880180605` | **F1 reliability anchor.** 5 subjects × 5 days × 10 MVCs: subjects failed to drive the elbow flexors maximally in 75% of ALL contractions (even fresh, non-fatigued) — VA<100% is common at baseline, not itself pathological; consistent between-subject differences (P<0.01); good within-subject reproducibility (4/5 subjects); explicitly not explained by antagonist co-contraction. |
| 6 | Bigland-Ritchie B, Jones DA, Hosking GP, Edwards RH (1978). Central and peripheral fatigue in sustained maximum voluntary contractions of human quadriceps muscle. *Clin Sci Mol Med* 54(6):609-14. | **657729**, `10.1042/cs0540609` | **F1 classic quadriceps dissociation.** 9 subjects, tetanic-stimulation-vs-voluntary-force comparison: 3/9 showed force declining in parallel with tetanic (little/no central fatigue, EMG constant while force fell ~60%); 6/9 showed force falling faster than tetanic with parallel EMG decline (central fatigue present), and could transiently boost force + EMG together with extra effort. |
| 7 | McNeil CJ, Martin PG, Gandevia SC, Taylor JL (2009). The response to paired motor cortical stimuli is abolished at a spinal level during human muscle fatigue. *J Physiol* 587(Pt23):5601-12. | **19805743**, `10.1113/jphysiol.2009.180968`, PMC2805373 | **F2 spinal-level double dissociation.** Paired-pulse TMS (long-interval intracortical inhibition) vs. cervicomedullary stimulation (bypasses cortex): both conditioned MEPs AND CMEPs collapse within 30s of a 2-min MVC (spinal, not purely cortical, mechanism); unconditioned responses recover in <30s but conditioned MEPs/CMEPs need 90s — a genuine time-course dissociation between sub-mechanisms. |
| 8 | Taylor JL, Todd G, Gandevia SC (2006). Evidence for a supraspinal contribution to human muscle fatigue. *Clin Exp Pharmacol Physiol* 33(4):400-5. | **16620309**, `10.1111/j.1440-1681.2006.04363.x` | **F2/F4.** Quantifies supraspinal fatigue as accounting for "about one-quarter of the loss of force" in sustained/intermittent elbow-flexor MVCs — a specific, muscle-group-scoped number, not claimed universal. Names fatigue-sensitive muscle afferents as a candidate contributor to supraspinal fatigue itself. |
| 9 | Taylor JL, Gandevia SC (2008). A comparison of central aspects of fatigue in submaximal and maximal voluntary contractions. *J Appl Physiol* 104(2):542-50. | **18032577**, `10.1152/japplphysiol.01053.2007` | **F4.** Central fatigue is demonstrable even in brief maximal efforts nested inside a weak (<15% max) ongoing contraction; group III/IV afferent activity acts supraspinally to limit motor cortical output; best submaximal-task evidence is a disproportionate rise in perceived effort. |
| 10 | Amann M, Dempsey JA (2008). Locomotor muscle fatigue modifies central motor drive in healthy humans and imposes a limitation to exercise performance. *J Physiol* 586(1):161-73. | **17962334**, `10.1113/jphysiol.2007.141838`, PMC2375542 | **F3 regulated-ceiling anchor.** 8 cyclists, two different pre-fatigue protocols (83% vs 67% Wpeak) before an identical 5km TT: central drive (quadriceps EMG) fell −23% with the higher pre-fatigue dose, but end-exercise peripheral fatigue (ΔQtw,pot) converged to the SAME −35% to −37% regardless (P=0.35) — the invariant-ceiling signature. |
| 11 | Amann M, Proctor LT, Sebranek JJ, Pegelow DF, Dempsey JA (2009). Opioid-mediated muscle afferents inhibit central motor drive and limit peripheral muscle fatigue development in humans. *J Physiol* 587(1):271-83. | **19015193**, `10.1113/jphysiol.2008.163303`, PMC2670040 | **F3 THE fentanyl-blockade paper.** Double-blind, placebo-controlled, 8 cyclists, 5km TT: intrathecal fentanyl (L3-L4) vs. placebo. iEMG (central drive) +12±3% in the first half (P<0.05); ΔQtw −46±2% (fentanyl) vs −33±2% (placebo), P<0.001; blood lactate higher (16.3 vs 12.6%), arterial O2 saturation lower (89 vs 94%) with fentanyl. Fentanyl had NO effect on resting quadriceps strength (rules out a direct muscle-level drug effect). |
| 12 | Amann M, Blain GM, Proctor LT, Sebranek JJ, Pegelow DF, Dempsey JA (2011). Implications of group III and IV muscle afferents for high-intensity endurance exercise performance in humans. *J Physiol* 589(Pt21):5299-309. | **21878520**, `10.1113/jphysiol.2011.213769`, PMC3225681 | **F3 exhaustion-protocol + direct-effect control + wrong-sign-performance nuance.** 7 males, constant-load cycling (80% Wpeak) to exhaustion. Explicit control: "unchanged resting cardioventilatory responses" (rules out a generic sedative confound). CMD +9±3% at end-exercise with fentanyl (P<0.05); ΔQtw,pot −44±2% (fentanyl) vs −34±2% (placebo); rate of peripheral fatigue development 67±10% greater (P<0.01); **time to exhaustion SHORTER** with fentanyl (6.8±0.3 vs 8.7±0.3 min) — losing the brake made overall performance worse, the opposite of a naive "opioid disinhibition helps you push harder" story, because O2 transport (ventilation/circulation) was simultaneously compromised (9% lower VE/VCO2, 5 mmHg higher end-tidal PCO2, 3% lower Hb saturation). |
| 13 | Blain GM, Mangum TS, Sidhu SK, Weavil JC, Hureau TJ, Jessop JE, Bledsoe AD, Richardson RS, Amann M (2016). Group III/IV muscle afferents limit the intramuscular metabolic perturbation during whole body exercise in humans. *J Physiol* 594(18):5303-15. | **27241818**, `10.1113/JP272283`, PMC5023715 | **F3 decorrelated leg — muscle biopsy, not EMG/force.** 8 subjects, 5km TT, vastus lateralis biopsies pre/post + surface EMG. Motoneuronal output (EMG) +21±12% with fentanyl (P<0.05). Intramuscular Pi +55%, H+ +62%, ADP +129%, lactate +47% (all P<0.001), PCr depletion +27% (P<0.01) greater with fentanyl. ΔQTsingle −52±2% (fentanyl) vs −31±1% (control), correlated with ΔPi (r²=0.79, P<0.01) and ΔH+ (r²=0.92, P<0.01). Time to complete the (fixed-distance) TT was SIMILAR between conditions (~8.8 min both) — a disclosed, genuine protocol-dependent contrast with #12's fixed-workload-to-exhaustion design, where time-to-exhaustion WAS shortened; task structure (self-paced fixed-distance vs. fixed-workload-to-failure) changes which performance metric the blockade affects. |
| 14 | Millet GY, Lepers R (2004). Alterations of neuromuscular function after prolonged running, cycling and skiing exercises. *Sports Med* 34(2):105-16. | **14965189**, `10.2165/00007256-200434020-00004` | **F4 central-vs-peripheral BY MODALITY.** Names "task dependency of muscle fatigue" as the organizing principle; central fatigue after prolonged running is attenuated in low-muscle-damage modalities (e.g. cycling) at matched duration; after a 30km run, individuals with the greatest knee-extensor strength loss showed the largest activation deficit specifically. |
| 15 | Enoka RM, Duchateau J (2008). Muscle fatigue: what, why and how it influences muscle function. *J Physiol* 586(1):11-23. | **17702815**, `10.1113/jphysiol.2007.139477`, PMC2375565 | **F4 conceptual anchor.** "There is no global mechanism responsible for muscle fatigue... the mechanisms that cause fatigue are specific to the task being performed" — the explicit theoretical grounding for why a single universal central:peripheral ratio should NOT be expected. |
| 16 | Kent-Braun JA, Sharma KR, Weiner MW, Massie B, Miller RG (1993). Central basis of muscle fatigue in chronic fatigue syndrome. *Neurology* 43(1):125-31. | **8423875**, `10.1212/wnl.43.1_part_1.125` | **F5 THE CFS anchor.** Tibialis anterior, intermittent submaximal + sustained maximal protocols. "Extent of fatigue, metabolic response, and changes in both M-wave amplitude and twitch tension... were similar in patients and controls... [and] response to systemic exercise was also normal" — BUT "voluntary activation of the tibialis was significantly lower in the patients during maximal sustained exercise," "well in excess of that found in controls." Authors' own conclusion: "an important central component of muscle fatigue in CFS." |
| 17 | Kent-Braun JA, Sharma KR, Weiner MW, Miller RG (1994). Effects of exercise on muscle activation and metabolism in multiple sclerosis. *Muscle Nerve* 17(10):1162-9. | **7935523**, `10.1002/mus.880171006` | **F6 MS activation-failure-without-proportional-metabolic-demand.** 6 MS + 8 controls, ankle dorsiflexors. MVC fell sooner in MS but similarly overall; end-exercise Pi/PCr rose to 1.86±0.22 in controls vs only 0.66±0.04 in MS (P<0.01); pH fell to 6.75±0.04 in controls vs unchanged 7.06±0.04 in MS (P<0.01) — a smaller metabolic perturbation at the same relative intensity, consistent with a genuine activation failure (not full peripheral engagement) even in MILD MS. |
| 18 | Sharma KR, Kent-Braun J, Mynhier MA, Weiner MW, Miller RG (1995). Evidence of an abnormal intramuscular component of fatigue in multiple sclerosis. *Muscle Nerve* 18(12):1403-11. | **7477063**, `10.1002/mus.880181210`, PMC2733338 | **F6 the forced adversary — MS is NOT purely central.** 9 min of INTERMITTENT ELECTRICAL STIMULATION (explicitly "used to eliminate central sources of muscle fatigue" — voluntary drive fully bypassed): tetanic-force, PCr, and pH decline were all GREATER in MS patients than controls; CMAP amplitude did not fall (rules out neuromuscular-junction failure); tetanic-force recovery was delayed post-exercise. Authors' own conclusion: MS fatigue has "both central (perception, upper motor neuron dysfunction) AND peripheral (impaired metabolism and excitation-contraction coupling) components." |
| 19 | Rice CL, Vollmer TL, Bigland-Ritchie B (1992). Neuromuscular responses of patients with multiple sclerosis. *Muscle Nerve* 15(10):1123-32. | **1406770**, `10.1002/mus.880151011` | **F6 the dramatic MS twitch-occlusion number.** 4 MS patients vs. normal controls, twitch occlusion: normal subjects activated their muscles maximally; MS patients "rarely achieved greater than 60% activation." Motoneuron firing rates in MS rarely exceeded 17 Hz vs ~24 Hz normal — consistent with the reduced activation, not an independent contradiction. |
| 20 | Meeusen R, Duclos M, Foster C, Fry A, Gleeson M, Nieman D, Raglin J, Rietjens G, Steinacker J, Urhausen A; European College of Sport Science; American College of Sports Medicine (2013). Prevention, diagnosis, and treatment of the overtraining syndrome: joint consensus statement. *Med Sci Sports Exerc* 45(1):186-205. | **23247672**, `10.1249/MSS.0b013e318279a10a` | **F7 overtraining consensus — the honest gap, stated by the field itself.** Explicitly: "several markers (hormones, performance tests, psychological tests, and biochemical and immune markers) are used, but none of them meet all the criteria to make their use generally accepted." No twitch-interpolation/VA-specific criterion is named. |
| 21 | Aubry A, Hausswirth C, Louis J, Coutts AJ, Buchheit M, Le Meur Y (2015). The Development of Functional Overreaching Is Associated with a Faster Heart Rate Recovery in Endurance Athletes. *PLoS One* 10(10):e0139754. | **26488766**, `10.1371/journal.pone.0139754`, PMC4619310 | **F7 decorrelated, suggestive (not diagnostic) link.** 31 triathletes; functional-overreaching subgroup (n=10) showed faster heart-rate recovery, explicitly interpreted by the authors as potentially "induced by a decreased central command and by a lower chemoreflex activity" — an autonomic-axis finding directionally consistent with a chronically altered central/afferent-feedback regulatory state, but NOT a voluntary-activation/twitch-interpolation measurement, and not on the same muscle-force construct as F1–F6. Cited as suggestive scope, not folded into the F1-F6 evidence chain. |
| 22 | Noakes TD (2001). Evidence that a central governor regulates exercise performance during acute hypoxia and hyperoxia. *J Exp Biol* 204(Pt18). | **11581338** | **Scoping/disambiguation, not certified here.** Title/journal/year/PMID confirmed live. Noakes' broader "Central Governor Model" extends to brain-mediated pacing/teleoanticipation across the whole exercise-regulation literature — a related but distinct, more theoretically expansive and more contested framework than the mechanistically-tested, pharmacologically-forced Amann group III/IV findings (#10-13) this document actually certifies. Named here explicitly so the two are not conflated. |
| 23 | Noakes TD (2007? — exact year not independently re-verified beyond esummary `pubdate:2007`). The central governor model of exercise regulation applied to the marathon. *Sports Med*. | **17465612** | Same scoping role as #22 — an applied review of the broader, more speculative theory. |

---

## 3. F1 — voluntary activation is a real, reproducible signal, and it falls with fatigue

**Adversary forced**: is the "VA deficit" just noise, submaximal motivation, or antagonist co-contraction
masquerading as a central signal? Allen, Gandevia & McKenzie (1995, #5) ran the decisive reliability check
directly: 5 subjects × 5 days × 10 maximal efforts each. Even in **fresh, non-fatigued** subjects, VA was
not 100% in 75% of all contractions — incomplete activation is a common, real, baseline trait, not a
fatigue artifact by itself — but it was **consistent between subjects** (P<0.01) and **reproducible within
a subject across days** (4/5 subjects), and explicitly **not explained by antagonist co-contraction**
(checked directly in that study). On top of that real, characterized baseline, Gandevia et al. (1996, #3)
measured a genuine **decline with fatigue**: VA of biceps brachii fell from >99% to an average 90.7% over a
3-minute sustained MVC. Bigland-Ritchie et al. (1978, #6) independently found the same phenomenon in a
completely different muscle group (quadriceps) and a different technique (tetanic stimulation vs. voluntary
force comparison, not twitch interpolation): 6 of 9 subjects showed force declining faster than the
electrically-stimulated tetanic reference, with EMG amplitude falling in parallel — and those same subjects
could **transiently boost both force and EMG together** with an extra conscious effort, a within-subject
demonstration that reserve capacity genuinely existed and was not being used. **Gate (VA decline ≥5
percentage points from a ≥95% baseline, ≥2 independent studies/muscle groups, reliability independently
established): PASS**, cleared with margin (8+ points in #3; a qualitative but unambiguous EMG-force
dissociation in 2/3 of the sample in #6).

## 4. F2 — the deficit has a genuine supraspinal (cortical) component, forced against a peripheral-artifact adversary

**Adversary forced**: could the TMS-evoked "extra force" simply be tracking a peripheral muscle-property
change (e.g., potentiation/depression) that has nothing to do with the brain, riding along with fatigue?
Gandevia et al. (1996, #3) built the decisive test directly into the same protocol: a sphygmomanometer cuff
blocked blood flow to the fatigued muscle immediately after the sustained MVC. **Voluntary force and VA did
NOT recover during this ischaemia** — but the TMS-evoked EMG-excitability changes (a decorrelated cortical
signal, recorded in the identical ischaemic block) **did recover rapidly**. A single shared peripheral
confound (e.g., a metabolite washout, contractile-protein state) cannot produce two signals with opposite
recovery behaviour under the identical manipulation — this is a forced double dissociation, not a
correlation. McNeil, Martin, Gandevia & Taylor (2009, #7) independently replicate the same logic with a
different sub-technique (paired-pulse long-interval intracortical inhibition vs. cervicomedullary
stimulation, which bypasses the cortex entirely): both conditioned MEPs and CMEPs collapse within 30s of a
2-minute MVC, but recover on a **different time-course** (90s) than the unconditioned responses (<30s) — a
second, independently-derived dissociation in time rather than in an ischaemic manipulation. Taylor, Todd &
Gandevia (2006, #8) quantify the resulting supraspinal contribution at roughly **one-quarter** of the total
force loss in sustained/intermittent elbow-flexor MVCs — a specific, scoped number, not a universal
constant. **Gate (TMS increment grows with fatigue AND fails to recover under ischaemia while a decorrelated
cortical signal in the identical ischaemia DOES recover): PASS**, independently corroborated by a second
technique/paper showing an analogous but methodologically distinct time-course dissociation.

## 5. F3 — group III/IV muscle afferents actively regulate central drive and peripheral-fatigue severity (the Amann governor)

**Adversary forced to its strongest fair form**: not "does blocking afferents change something" (a weak
test) but the two hardest versions — (a) is the fentanyl effect just a generic sedative/CNS-depressant or
disinhibition/analgesia artifact rather than true afferent-signal removal, and (b) is the resulting
"regulation" claim actually beneficial/adaptive, or could it be dressed-up passive muscle fatigability with
no causal feedback role at all?

Against (b): Amann & Dempsey (2008, #10) show the signature that only a real regulator produces — two
different pre-fatigue *paths* (67% vs 83% Wpeak) converge on the **same** end-exercise peripheral-fatigue
ceiling (−35% to −37%, P=0.35 no difference) while central drive is dose-dependently suppressed in
proportion to how much pre-existing fatigue exists. A passive, unregulated fatigability limit does not
predict convergence from different starting trajectories; a clamped feedback loop does.

Against (a): Amann et al. (2011, #12) built the direct control — "unchanged resting cardioventilatory
responses" to fentanyl, ruling out a generic depressant effect at rest; the measured changes (CMD, ΔQtw,
ventilation, blood pressure) only emerge **during exercise**, exactly when the blocked afferents would
normally be firing. And the sharpest test of all is the **sign** of the performance effect: if fentanyl
simply disinhibited effort/reduced perceived pain generically, overall performance should improve. Instead,
time-to-exhaustion was **shorter**, not longer (6.8 vs 8.7 min) — because losing the afferent brake also
removes the ventilatory/circulatory support those same afferents normally recruit (9% lower VE/VCO2, lower
arterial O2 saturation, lower blood pressure during exercise). This is a forced adversary getting the
**wrong sign**, not a marginal miss.

Numbers, replicated across three independent trials/protocols, one of them (#13) on a completely
**decorrelated measurement axis** (intramuscular metabolite chemistry from muscle biopsy, not EMG or
force):

| study | design | central drive (blocked vs intact) | peripheral fatigue, ΔQtw (blocked vs intact) | other |
|---|---|---|---|---|
| Amann 2009 (#11) | 5km TT | iEMG +12±3% (1st half) | −46±2% vs −33±2% (P<0.001) | lactate 16.3 vs 12.6%; SaO2 89 vs 94% |
| Amann 2011 (#12) | fixed-load-to-exhaustion | CMD +9±3% (end-exercise) | −44±2% vs −34±2%; rate +67±10% (P<0.01) | time-to-exhaustion 6.8 vs 8.7 min (SHORTER) |
| Blain 2016 (#13) | 5km TT + biopsy | EMG +21±12% (P<0.05) | −52±2% vs −31±1% (P<0.001) | Pi/H+/ADP/lactate/PCr +27–129% greater; r²=0.79 (Pi), 0.92 (H+) with ΔQT |

**Gate (blockade must ↑ central drive AND ↑ peripheral fatigue across ≥2 protocols, confirmed on a
decorrelated non-EMG/force leg, correlated with the twitch deficit): PASS**, all three replications land the
same direction with wide margins, and the metabolite leg is measured by an entirely different instrument
(mass spectrometry on biopsy tissue) than the EMG/twitch-force leg.

## 6. F4 — central-vs-peripheral proportion is task-dependent, not a universal ratio

Millet & Lepers (2004, #14) name this explicitly "task dependency of muscle fatigue": central fatigue after
prolonged running is measurable via twitch interpolation, EMG/M-wave ratio, or voluntary-vs-evoked force
comparison, but is **attenuated** in activities inducing less muscle damage (e.g. cycling) at matched
duration; after a 30km run, the subjects with the greatest knee-extensor strength loss specifically showed
the largest activation deficit — a within-modality dose-relationship, not just a between-modality one.
Taylor & Gandevia (2008, #9) show the phenomenon is present even in brief maximal efforts nested inside weak
(<15% max) ongoing contractions, and that its best submaximal-task signature is a **disproportionate rise
in perceived effort** relative to force output — a distinct observable from the maximal-contraction VA
decline. Enoka & Duchateau (2008, #15) state the field's own governing principle directly: "there is no
global mechanism responsible for muscle fatigue... the mechanisms... are specific to the task being
performed." **Gate (≥2 decorrelated modalities/groups show materially different apportionment, defeating a
fixed-universal-ratio adversary): PASS** — running/cycling/skiing differ, and the specific ~25% figure
(#8) for sustained elbow-flexor MVC is explicitly scoped to that muscle group and task, not generalized.

## 7. F5 — chronic fatigue syndrome: a clean, axis-orthogonal CENTRAL-only deficit

Kent-Braun, Sharma, Weiner, Massie & Miller (1993, #16) is the decisive dissociation: in CFS patients vs.
matched controls performing intermittent submaximal and sustained maximal tibialis-anterior protocols, the
**entire peripheral axis is normal** — "extent of fatigue, metabolic response, and changes in both M-wave
amplitude and twitch tension... were similar in patients and controls," and the systemic (whole-body)
exercise response was also normal — while voluntary activation was **specifically and significantly lower**
in patients during the sustained maximal condition, "well in excess of that found in controls." This is a
clean orthogonal-axis result: the peripheral coordinate is unchanged, the central coordinate alone is
shifted, which is exactly the pattern a genuine central-specific deficit predicts and a generalized
deconditioning/peripheral-myopathy account of CFS fatigue does not. **Gate (VA significantly lower WHILE
peripheral markers are statistically indistinguishable from controls): PASS.**

## 8. F6 — multiple sclerosis: forced against a "purely central" adversary, and it does NOT fully fall — a genuine mixed verdict

The central-axis deficit in MS is, if anything, larger and more dramatic than in CFS: Rice, Vollmer &
Bigland-Ritchie (1992, #19) found normal subjects could always activate their muscles maximally under
twitch occlusion, while **MS patients rarely achieved greater than 60% activation** — a huge central
shortfall — and Kent-Braun et al. (1994, #17) independently show MS subjects reach a smaller metabolic
perturbation (Pi/PCr 0.66 vs 1.86, pH unchanged 7.06 vs 6.75 in controls, both P<0.01) at the same relative
exercise intensity, consistent with genuinely incomplete muscle engagement rather than full peripheral
loading. **But the adversary "MS fatigue is purely central/upper-motor-neuron" is then forced onto its
hardest test and fails**: Sharma, Kent-Braun, Mynhier, Weiner & Miller (1995, #18) ran 9 minutes of
intermittent **electrical stimulation** — explicitly designed "to eliminate central sources of muscle
fatigue," i.e. voluntary drive is bypassed by construction — and MS patients STILL showed a greater decline
in tetanic force, phosphocreatine, and pH than controls, with delayed post-exercise force recovery
(consistent with a genuine excitation-contraction-coupling abnormality) while CMAP amplitude stayed stable
(ruling out neuromuscular-junction failure as the cause). The authors' own conclusion, quoted directly:
fatigue in MS has "both central (perception, upper motor neuron dysfunction) and peripheral (impaired
metabolism and excitation-contraction coupling) components." **Gate (pure-central adversary must fail
against a voluntary-drive-bypassing protocol showing excess MS fatigue): PASS-as-adversary-partially-
survives → reported here as a forced, symmetric MIXED verdict**, not oversold as "MS proves central fatigue
alone."

## 9. F7 — overtraining syndrome: the honest, OODA-forced gap

Two targeted search angles were run this session specifically to force this before declaring a gap
(Observe: does a direct twitch-interpolation/VA study in diagnosed overtraining syndrome exist? Orient: the
disease is defined behaviourally/hormonally, not electrophysiologically, in the primary consensus
literature — so the search targeted both "overtraining"+voluntary-activation and "overreaching"+central-
fatigue phrasing directly, Title-restricted to reduce false negatives from generic keyword co-occurrence):
neither surfaced a dedicated VA/twitch-interpolation study in a diagnosed-OTS population. This is not merely
this session's search failing — Meeusen et al. (2013, #20), the field's own joint ECSS/ACSM consensus
statement, states explicitly that of the several markers in current use (hormonal, performance,
psychological, biochemical/immune) "none of them meet all the criteria to make their use generally
accepted" — i.e. the field itself has not converged on a validated marker, let alone specifically a
central-activation one. The single decorrelated finding located, Aubry et al. (2015, #21), shows a faster
heart-rate-recovery signature in functional overreaching that the authors themselves interpret as possibly
reflecting "decreased central command and... lower chemoreflex activity" — directionally consistent with a
chronically altered central/afferent-regulatory state, but measured via autonomic/cardiovascular markers,
not voluntary muscle force or twitch interpolation, and in *functional* overreaching (the mild, reversible
end of the spectrum) rather than diagnosed overtraining syndrome. **Verdict: OPEN, honestly reported as a
gap after being forced once, not a zero-effort "honest negative."**

---

## 10. Scoping note — this document certifies Amann's mechanistic findings, NOT the broader "Central Governor Model"

Noakes' Central Governor Model (#22, #23) is a related, more theoretically expansive claim — that the brain
paces essentially all exercise performance via an anticipatory, subconscious regulatory algorithm across
distances/environments, extending well beyond the specific group III/IV afferent-feedback mechanism tested
here. That broader model remains actively contested in the exercise-physiology literature. This document
certifies specifically what #10-13 measured: a pharmacologically-forced, decorrelated-leg-confirmed,
sign-tested causal role for group III/IV muscle afferent feedback in capping central drive and peripheral
fatigue during a specific class of high-intensity cycling protocols. Naming the broader theory here is a
deliberate scope fence, not a citation of supporting evidence for it.

## 11. Honest gaps (disclosed, not discovered after the fact)

- **The exact VA% formula was not independently re-derived from a live-quoted primary source this
  session** — `VA% = [1 − (SIT/RT)] × 100` is presented as the field-standard convention (consistent with
  every twitch-interpolation study cited), not extracted verbatim from a fetched abstract; unlike the
  cosine-tuning formula in this repo's `MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md`, none of the abstracts
  fetched this session spelled out the algebra explicitly, and full-text access to the classic
  methodology papers (#4, #5) was not attempted this pass.
- **Merton (1954, #1) bibliographic-only.** Pre-abstract-era PubMed record; PMC full-text fetch returned a
  publisher restriction ("does not allow downloading of the full text in XML form" — the identical
  restriction this repo's sibling docs hit on other pre-1990s PMC records). The paper's own original
  quantitative finding (which muscle, which specific numbers) was NOT independently re-extracted this
  session — cited for provenance/priority only, not for a specific number.
- **No new simulation/measurement was run this pass** — this is a literature-forced synthesis (status
  OPEN, `mechanism_grade: HYPOTHESIS-B`), not a from-scratch numeric model; unlike
  `docs/MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md`'s own machine-gated simulation, no code was written here to
  independently reproduce a twitch-interpolation VA measurement or a group-III/IV-afferent feedback-loop
  model on this twin's own data. The proposed_cell (§13) names the decisive, cheap first executable step.
- **F6's "mixed" verdict is not itself quantitatively decomposed** — Sharma 1995 (#18) and Rice 1992 (#19)
  are decorrelated studies (different protocols, likely different specific patients) showing central AND
  peripheral deficits separately in MS; no single study here partitions a given patient's total fatigue into
  an exact central:peripheral percentage split, so "mixed" is a qualitative, forced-adversary-survives
  verdict, not a quantified ratio.
- **F7 is a genuine, currently-open field gap**, not merely a one-session search miss — disclosed at length
  in §9 rather than papered over with a weaker adjacent citation stretched to fit.
- **Group III/IV muscle afferents (Amann's target, metabo-/mechanoreceptive, from working locomotor
  muscle) are anatomically and functionally DISTINCT from the group III/IV joint afferents already cited in
  this repo's `docs/MECHANISM_NOCICEPTION.md`** (Grigg/Schaible/Schmidt 1986, posterior articular nerve,
  cat knee) — same numerical/conduction-velocity classification scheme, different receptive organ and
  reflex target. Not conflated here; flagged explicitly as a couples_to disambiguation (§12).
- **All protocols cited are cycling/isometric-elbow-flexor/isometric-quadriceps/tibialis-anterior
  paradigms in adult humans** (plus rodent-free — no animal data used in this document, unlike the sibling
  eccentric-damage doc). Running/skiing central-fatigue findings are cited via the Millet & Lepers (2004)
  review, not independently re-verified against primary running/skiing studies this session.
- **Sample sizes are small across the decisive Amann trials** (n=7-8 per study, #10-13) — real, published,
  peer-reviewed, replicated across 3 independent papers from the same group with a decorrelated leg, but
  not a large-cohort finding; effect sizes are large and consistent, which is why the falsifier gate is set
  on direction+magnitude-across-replications rather than a single p-value.
- **No twitch-interpolation study cited here was independently re-analyzed on raw trial-level data** — every
  number is taken from the source paper's own reported abstract statistic (means±SD, percentages), not
  recomputed from a downloaded dataset; this document is a literature synthesis, consistent with its
  disclosed OPEN/HYPOTHESIS-B status.

## 12. Couples to

- **`docs/MECHANISM_MUSCLE_FATIGUE.md`** (this repo, same wave) — the PERIPHERAL, motor-unit-compartment
  (Xia & Frey Law) fatigue/recovery ODE this document's central-drive layer sits upstream of; that doc's own
  §10 explicitly names "no central fatigue... modeled" as its scope limit — this document is the designed
  fill for exactly that gap, not a duplicate (a genuinely different mechanism: reduced neural DRIVE vs.
  reduced muscle CAPACITY given a drive).
- **`MSK-MUSCLE-FATIGUE-EMG`** (`docs/MECHANISM_STATE.md`) — a third, non-overlapping capability: a surface-
  EMG spectral-shift (median-frequency-down/RMS-up) fatigue *signature* detector. Detects that fatigue is
  happening; does not partition how much of it is central vs peripheral. This document supplies the missing
  partition.
- **`docs/MECHANISM_MOTOR_UNIT.md`** — supplies the recruitment/rate-coding manifold whose unused (null-space)
  capacity twitch interpolation directly measures (§1A); this document's "central drive" is that document's
  "excitation" input, now shown to be itself regulated and imperfect rather than a free parameter.
- **`docs/MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md`** — the descending-drive layer this document's cortical
  (TMS-revealed, #3/#7) and spinal (#7) fatigue mechanisms directly modulate; that document's population-
  vector/conduction-delay model describes the healthy, unfatigued baseline this document's supraspinal-
  fatigue findings are a *deviation from*, in exactly the function-vs-dysfunction sense that document's own
  stroke section uses for a different pathology.
- **`docs/MECHANISM_ECCENTRIC_MUSCLE_DAMAGE.md`** — the sibling athlete-training-load cert (repeated-bout
  effect, DOMS, rhabdomyolysis dysfunction pole); central fatigue and eccentric damage are DECORRELATED
  history states affecting the same athlete-movement axis (one a same-session neural-drive limit, the other
  a days-to-months structural damage/repair cycle) — both instantiate the FUNCTION-vs-DYSFUNCTION organizing
  axis (`bt_memory/function-vs-dysfunction-build-the-athlete-baseline-in-parallel-with-disease.md`) on the
  same athlete substrate, from opposite ends of the timescale spectrum.
- **`docs/MECHANISM_VO2MAX_AEROBIC_CAPACITY.md`** — Amann's group III/IV findings (#10-13) are measured on
  exactly the cycling-endurance-performance construct that document's Fick-product ceiling models; a
  cardiovascular-delivery limitation and a central-motor-drive limitation are DECORRELATED candidate
  ceilings on the same overall endurance-performance outcome, not yet jointly modeled in one cell.
  `docs/MECHANISM_ATHLETE_HEART_REMODELING.md` shares the same athlete-endurance context.
- **`docs/MECHANISM_NOCICEPTION.md` / `docs/MECHANISM_PAIN_NOCICEPTION.md`** — group III/IV afferent
  nomenclature is shared with this repo's existing joint-nociceptor content, but the receptive organ and
  reflex target are DIFFERENT (muscle metabo-/mechanoreceptors here vs. joint capsule nociceptors there) —
  disambiguated explicitly, not merged (§11).
- **`bt_memory/function-vs-dysfunction-build-the-athlete-baseline-in-parallel-with-disease.md`** — this
  document instantiates the axis directly: §5's Amann governor is the validated FUNCTION-pole mechanism;
  §7-8's CFS/MS findings are the DYSFUNCTION pole (disproportionate/qualitatively-different central-fatigue
  component); §9's overtraining gap is the honestly-unresolved boundary between the two.

## 13. Proposed next cell (not built this pass)

The cheapest, most decisive next executable step: extend `scripts/msk/muscle_fatigue.py` (the existing
peripheral 3-compartment motor-unit model, §12) with a **central-drive multiplier** `CD(t) ∈ (0, 1]` applied
upstream of the existing `TL(t)` target-load input, parametrized by (i) a baseline-VA distribution drawn
from the Allen/Gandevia/McKenzie (1995, #5) reliability data (VA<100% in 75% of fresh contractions) and (ii)
a fatigue-dependent decay term calibrated to the Gandevia et al. (1996, #3) 3-minute-MVC trajectory
(>99%→90.7%). Decisive test: does adding `CD(t)` change which muscle in this twin's own 80-muscle gait
reconstruction is predicted to cross its own peak-demand/MF ceiling first (the mechanism
`MECHANISM_MUSCLE_FATIGUE.md` §7 already derived from its own conservation law) — a cheap, fully
machine-checkable comparison against data already on disk, zero new acquisition required, exactly parallel
to the eccentric-damage document's own proposed next-cell pattern.

## 14. Files

- This document.
- `docs/MECHANISM_CENTRAL_FATIGUE_evidence.json` — machine-readable falsifier ledger (7 tests), the
  geometric derivation, the citation table (23 PMIDs), couples_to, honest gaps, and the proposed next cell.
