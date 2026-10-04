# MECHANISM ALTITUDE ACCLIMATIZATION — the integrated hypoxia response, time-ordered (2026-07-22)

**Status: HYPOTHESIS (literature-synthesis + self-contained arithmetic, `status:OPEN` per this repo's node
schema) — awaiting an independent QC pass.** Models the O2 cascade (PIO2 → PAO2 via the alveolar gas
equation → PaO2 → SaO2 via the O2-Hb curve → CaO2) at extreme altitude, and the TIME-ORDERED acclimatization
response layered on top of it (HVR minutes → 2,3-BPG/EPO-signal hours → renal compensation days → hematocrit
weeks → full RBC turnover ~4 months). Distinct from, and coupled to, this repo's own `MECHANISM_ERYTHROPOIESIS.md`
(the EPO→RBC loop alone, built for a blood-donation/hemorrhage trigger, explicitly scoped **out** altitude).
15 PMIDs live-verified via NCBI eutils this session, 4 Wikipedia pages live-fetched (flagged textbook-tier),
plus verbatim reuse of already-certified constants from 4 sibling docs. Machine-readable companion:
`docs/MECHANISM_ALTITUDE_ACCLIMATIZATION_evidence.json`.

---

## 0. Scope note — the gap this doc closes, checked live before writing a line

`data/MECHANISM_ANCHOR_GRAPH.json` (998+ nodes) was searched live first (whole-word regex
`ALTITUDE|HYPOXIA|HVR|HACKETT|HAPE|HACE|MOUNTAIN SICKNESS`) — **zero** altitude-acclimatization nodes exist
(the only tangential hit, `ONCO-ANGIOGENESIS-HYPOXIA`, is an unrelated cancer-biology node). More tellingly,
**three existing sibling docs each independently disclosed the identical open gap**: `MECHANISM_RESPIRATORY.md`
("NO altitude / inspired-O2 coupling. Sea-level, normoxic constants only"), `MECHANISM_ERYTHROPOIESIS.md`
("no 2,3-DPG/altitude/fetal-Hb coupling"), and `MECHANISM_BLOOD_OXYGEN_TRANSPORT.md` ("no ... altitude
coupling"). This doc is the first to close it — by **reusing**, not duplicating, each sibling's own
already-validated machinery (the Hill/Bohr O2-Hb curve, the Henderson-Hasselbalch acid-base model, the
EPO/marrow delayed-feedback dynamics), rather than re-deriving any of them from scratch.

## 1. Geometric structure — the compositional argument, stated up front

The O2 cascade is a **composition of maps**:

```
PIO2 --(alveolar gas eq: f(Pb, FIO2, PaCO2, R))--> PAO2 --(diffusion)--> PaO2
      --(Hill/Bohr ODC: f(PaO2, pH, PaCO2, T))--> SaO2 --(affine, ×[Hb])--> CaO2 --(×Q, Fick)--> DO2
```

The erythropoiesis/[Hb] variable enters this chain **only at the CaO2 stage** — it is structurally absent
from the *domain* of every upstream map. This is not a timescale finding or a modeling choice; it is a fact
about which variables each stage's function actually takes as arguments. `SaO2 = PaO2^n/(PaO2^n + P50(pH,
PaCO2,T)^n)` simply has no `[Hb]` argument. §5 verifies this by direct computation (not assumption): SaO2 is
**identical to 4 decimal places** across a 10→24 g/dL Hb sweep (severe anemia to CMS-grade polycythemia).
This is the geometric backbone of the whole doc: the ventilatory/acid-base arm (PaCO2, pH) is the **only**
lever that can defend PaO2/SaO2; the erythropoietic arm's lever ([Hb]) is one stage downstream and can only
ever rescale O2 *content*, never saturation or partial pressure — at ANY timescale, however fast (instant
hemoconcentration) or slow (weeks-scale true erythropoiesis) its own mechanism runs.

**Reused, not new**: the Hill-equation form and its P50=26.6 mmHg / n=2.7 parameterization are
`MECHANISM_BLOOD_OXYGEN_TRANSPORT.md`'s own already-validated constants (Collins et al. 2015, PMID 26632351);
the CaO2 formula and its analytic `dCaO2/dHb = 1.34×SaO2` derivative are that same doc's own §10 result,
reproduced here (§5) to 4 significant figures.

## 2. Method, one paragraph

Model the cascade using the alveolar gas equation (`PAO2 = FIO2·(Pb−PH2O) − PaCO2/R`) driven by the real,
independently-measured Everest-summit barometric pressure (West 1999, PMID 10066724). Check it against three
legs of real data (§4): the OEII chamber's own calibration target, West et al. 1983's own derived PAO2, and
an arterial PaO2 convergence between two methodologically decorrelated studies. Force the "EPO/hematocrit-only"
adversary to its **strongest fair form** — not "erythropoiesis is slow" (a weak timescale argument) but the
**fastest possible blood-side confound**, hemoconcentration (§5) — and show it falls via an exact partial-derivative
computation plus a counterfactual/void-floor sweep proving the ventilatory arm is *necessary*, not just faster.
Build the time-ordered cascade (§6) from primary sources, several of which directly reuse this repo's own
already-certified `erythropoiesis.py` dynamical constants. Close with the two symmetric dysfunction poles
(§7 AMS/HAPE/HACE = fast arms incomplete; §8 CMS = slow arm overshooting) and the functional live-high-train-low
pole (§9).

## 3. Citations — 15 PMIDs live-verified via NCBI eutils this session, 4 Wikipedia pages live-fetched (flagged)

| # | Citation | PMID/DOI | Role |
|---|---|---|---|
| 1 | Sutton JR, Reeves JT, Wagner PD, et al. (1988). "Operation Everest II: oxygen transport during exercise at extreme simulated altitude." *J Appl Physiol* 64(4):1309-21. | **3132445**, DOI 10.1152/jappl.1988.64.4.1309 | 40-day decompression-chamber ascent, n=8; VO2max 3.98→1.17 L/min; arterial PO2=28±1, PCO2=11±1 Torr at 60W/PIO2=43 Torr; fourfold ventilation increase; cardiac output rose only at the most extreme tier |
| 2 | West JB, Boyer SJ, Graber DJ, et al. (1983). "Maximal exercise at extreme altitudes on Mount Everest." *J Appl Physiol* 55(3):688-98. | **6415008**, DOI 10.1152/jappl.1983.55.3.688 | AMREE; PIO2=43 Torr "equivalent to that on the summit"; alveolar PCO2 7-8 Torr; VO2max ~1 L/min; diffusion limitation demonstrated |
| 3 | West JB, Hackett PH, Maret KH, et al. (1983). "Pulmonary gas exchange on the summit of Mount Everest." *J Appl Physiol* 55(3):678-87. | **6415007**, DOI 10.1152/jappl.1983.55.3.678 | 34 alveolar samples (4 on the actual summit); PACO2=7.5 Torr; PAO2(calc)=35 Torr; arterial pH>7.7; PaO2(calc)=28 Torr |
| 4 | West JB (1984). "Human physiology at extreme altitudes on Mount Everest." *Science* 223(4638):784-8. | **6364351**, DOI 10.1126/science.6364351 | Summary corroboration of #3's own numbers |
| 5 | West JB (1999). "Barometric pressures on Mt. Everest: new data and physiological significance." *J Appl Physiol* 86(3):1062-6. | **10066724**, DOI 10.1152/jappl.1999.86.3.1062 | Summit Pb=253 Torr, directly measured (1981) + reconfirmed within ~1 Torr (1997) + weather-balloon corroboration |
| 6 | Hackett PH, Roach RC (2001). "High-altitude illness." *N Engl J Med* 345(2):107-14. | **11450659** | Modern AMS/HAPE/HACE framework; no abstract indexed (structured review) — framework citation only |
| 7 | Hackett PH, Rennie D, Levine HD (1976). "The incidence, importance, and prophylaxis of acute mountain sickness." *Lancet* 2(7996):1149-55. | **62991**, DOI 10.1016/s0140-6736(76)91677-9 | PRIMARY, n=278, Pheriche 4243m: 53% AMS, 7 HAPE+5 HACE, ascent-rate-dependent |
| 8 | Roach RC, Hackett PH, Oelz O, et al. (2018). "The 2018 Lake Louise Acute Mountain Sickness Score." *High Alt Med Biol* 19(1):4-6. | **29583031**, DOI 10.1089/ham.2017.0164 | Consensus revision; removed "disturbed sleep" as a scored item |
| 9 | Bistsch P, Swenson ER (2013). "Clinical practice: Acute high-altitude illnesses." *N Engl J Med* 368(24):2294-302. | **23758234**, DOI 10.1056/NEJMcp1214870 | Clinical-vignette format, no extractable numeric abstract — framework citation only |
| 10 | Powell FL, Milsom WK, Mitchell GS (1998). "Time domains of the hypoxic ventilatory response." *Respir Physiol* 112(2):123-34. | **9716296**, DOI 10.1016/s0034-5687(98)00026-7 | Formal seconds-to-years time-domain framework for HVR |
| 11 | Lenfant C, Torrance J, English E, et al. (1968). "Effect of altitude on oxygen binding by hemoglobin and on organic phosphate levels." *J Clin Invest* 47(12):2652-6. | **5725278**, DOI 10.1172/JCI105948 | PRIMARY: 2,3-BPG/Hb-affinity shift within 24h of altitude change |
| 12 | Eckardt KU, Boutellier U, Kurtz A, et al. (1989). "Rate of erythropoietin formation in humans in response to acute hypobaric hypoxia." *J Appl Physiol* 66(4):1785-8. | **2732171**, DOI 10.1152/jappl.1989.66.4.1785 | PRIMARY, n=6: EPO significantly elevated by 84-114 min; 1.8-3.0-fold production-rate rise; half-life 5.2h |
| 13 | León-Velarde F, Maggiorini M, Reeves JT, et al. (2005). "Consensus statement on chronic and subacute high altitude diseases." *High Alt Med Biol* 6(2):147-57. | **16060849**, DOI 10.1089/ham.2005.6.147 | CMS/subacute mountain sickness framework |
| 14 | Beall CM (2007). "Two routes to functional adaptation: Tibetan and Andean high-altitude natives." *PNAS* 104(Suppl 1):8655-60. | **17494744**, DOI 10.1073/pnas.0701985104 | PRIMARY: divergent evolutionary routes to the same O2-delivery outcome |
| 15 | Levine BD, Stray-Gundersen J (1997). "'Living high-training low'..." *J Appl Physiol* 83(1):102-12. | **9216951**, DOI 10.1152/jappl.1997.83.1.102 | PRIMARY, randomized n=39: VO2max +5% ∝ red-cell-mass +9% (r=0.37); 5000m time improved only in high-low group |
| — | Wikipedia "Alveolar gas equation", "Chronic mountain sickness", "Altitude sickness", "Respiratory alkalosis" | — | **TEXTBOOK-GRADE, flagged.** The first page's own quoted "107 mmHg" sea-level worked example was NOT reproduced by this doc's independent recomputation (§4) — caught and rejected, not propagated (§10). |

**Reused verbatim (not re-verified this session)**: `MECHANISM_BLOOD_OXYGEN_TRANSPORT.md`'s Hill/CaO2
constants (Collins 2015, PMID 26632351); `MECHANISM_ACID_BASE_CO2.md`'s Henderson-Hasselbalch form and its
own disclosed-uncertain Bohr-coefficient sweep convention; `MECHANISM_ERYTHROPOIESIS.md`'s marrow-transit
delay (τ=6d), 15d/193d recovery range, and RBC lifespan (Mock 2011, PMID 21062290); `MECHANISM_VO2MAX_AEROBIC_CAPACITY.md`'s
Fick product and its Buick 1980 (PMID 7380690)/Bárány 1993 (PMID 8482049) O2-content-raises-VO2max citations.

## 4. C1 — the O2 cascade reproduces Everest summit and OEII, cross-checked by three legs

**Claim**: the alveolar gas equation, driven by the real measured summit Pb, reproduces (a) the OEII chamber's
own calibration target, (b) West et al. 1983's own derived PAO2, and (c) an arterial PaO2 convergence between
two decorrelated studies.

**Pre-registered threshold**: within ±5% of the independently-published value, on each of 3 legs.

| leg | computed (this doc) | independently published | diff |
|---|---:|---|---:|
| PIO2 at Pb=253 Torr (West 1999) | **43.12 Torr** | OEII's own calibration target: 43 Torr (Sutton 1988/West 1983) | **0.27%** |
| PAO2, using West's own PACO2=7.5 Torr, R=0.85 | **34.29-34.57 Torr** (simplified/full form) | West et al. 1983's own stated 35 Torr | **−1.2% to −2.0%** |
| Arterial PaO2 | measured **28 Torr**, West 1983 (real summit, resting) | measured **28 Torr**, Sutton 1988 (OEII chamber, 60W exercise) | **0.0%** |
| Sea-level self-consistency | PAO2=99.2-101.3 mmHg | universal ~100 mmHg teaching value | **<1.3%** |

**Gate: PASS, all 3 legs.** The PaO2=28-Torr convergence is this doc's **decorrelated anchor**: two
methodologically independent programs — a real Everest summit ascent (elite climbers, alveolar-gas sampling,
rest) vs. a decompression-chamber simulation (healthy volunteers, direct arterial catheter, 60W exercise) —
agree exactly, via the shared physics of the alveolar gas equation, without either number deriving the other.

## 5. C2 — THE DECISIVE FALSIFIER: the ventilatory arm is necessary; erythropoiesis is structurally excluded

**Claim**: immediate-to-day-scale PaO2/SaO2 defense is achieved exclusively by the ventilatory/acid-base arm.
The erythropoiesis/hematocrit arm is not merely *slower* — it is **structurally excluded** from this outcome
variable at any timescale.

**Adversary, forced to its strongest fair form**: not "erythropoiesis is slow" (a weak, timescale-only
argument that a sufficiently patient adversary could dodge) but the **fastest possible blood-side confound**:
hemoconcentration — plasma-volume contraction raising [Hb] within *hours*, requiring no new RBC production
at all. If even this fastest blood-side mechanism cannot move SaO2/PaO2, the exclusion is airtight regardless
of which blood-side mechanism or timescale anyone proposes.

**Test 1 — exact geometric proof (not approximate).** SaO2 and CaO2 computed at the measured summit
PaO2=28 Torr across Hb=10→24 g/dL (severe anemia to CMS-grade polycythemia):

| Hb (g/dL) | SaO2 (%) | CaO2 (mL/dL) |
|---:|---:|---:|
| 10 | 53.4568 | 7.247 |
| 12 | 53.4568 | 8.680 |
| 15 | 53.4568 | 10.829 |
| 18 | 53.4568 | 12.978 |
| 21 | 53.4568 | 15.127 |
| 24 | 53.4568 | 17.276 |

SaO2 is **identical to 4 decimal places** at every Hb value; CaO2 rises monotonically (analytic
`dCaO2/dHb = 1.34×SaO2 = 0.7163` mL/dL per g/dL, matching `MECHANISM_BLOOD_OXYGEN_TRANSPORT.md`'s own formula
exactly). **This is where erythropoiesis DOES matter** (content, delivery, exercise capacity — §9) and
**exactly where it cannot** (saturation, partial pressure).

**Test 2 — counterfactual impossibility, R-swept.** Hold PaCO2 at the sea-level-typical 40 Torr at
Pb=253 Torr (i.e., *no* ventilatory increase): PAO2 = **−3.94 Torr** — negative, physiologically impossible,
regardless of Hb. Swept across the physiological R range:

| R | critical PaCO2 (PAO2=0 boundary) | sea-level-typical PaCO2 (40) exceeds it by |
|---:|---:|---:|
| 0.80 | 34.49 Torr | 5.51 Torr |
| 0.85 | 36.65 Torr | 3.35 Torr |
| 0.90 | 38.80 Torr | 1.20 Torr |

The impossibility holds at every swept R — not an artifact of one chosen value.

**Test 3 — void-floor sweep on big margins (PaCO2 5→45 Torr):**

| PaCO2 | PAO2 @ Pb=760 (sea level) | PAO2 @ Pb=253 (Everest summit) |
|---:|---:|---:|
| 5 | 143.35 | 37.23 |
| 15 | 131.58 | 25.47 |
| 25 | 119.82 | 13.70 |
| 35 | 108.05 | 1.94 |
| **40** | 102.17 | **−3.94 (impossible)** |
| 45 | 96.29 | −9.83 (impossible) |

At sea level the **entire** swept range stays strongly positive — ventilation barely matters for survival
(huge margin). At Everest-summit pressure the identical range crosses zero right at the sea-level-normal
value — a genuine **regime-dependent** finding (ventilation is existentially decisive only in the extreme
regime), not a fixed, tunable threshold.

**Test 4 — Bohr-shift benefit, direction-robust across 6 swept coefficients** (reusing
`MECHANISM_ACID_BASE_CO2.md`'s own disclosed-uncertain-coefficient convention, applied to the alkalosis
direction): at measured summit pH>7.7 (point estimate 7.70), P50 falls from 26.6 to 16.4-23.2 mmHg depending
on the swept coefficient, raising SaO2 at the fixed measured PaO2=28 Torr from 53.46% (unshifted) to
62.5-80.9% — a substantial, monotonic benefit at every one of 6 swept coefficients. This quantifies West's
classic teaching point: the "paradoxical" extreme alkalosis is adaptive, not just a side effect — it buys
back saturation exactly where the ventilatory arm needs it (the lung), at a cost paid elsewhere (tissue
unloading) that matters less at this extremity.

**Gate: PASS on all 4 tests.** The erythropoiesis-only adversary falls in its strongest fair form (fast
hemoconcentration, not slow true erythropoiesis) via an exact, non-approximate structural argument, not a
timescale appeal — and the ventilatory arm is shown *necessary* (not just empirically faster), robust across
the swept R range, with a genuine regime-dependent void floor.

## 6. C3 — the time-ordered cascade, each stage independently anchored

| stage | time constant | source |
|---|---|---|
| HVR (carotid body) onset | **seconds to minutes** | classic teaching / Powell, Milsom & Mitchell 1998 (PMID 9716296) time-domain framework |
| EPO **signal** onset | **84-114 min**, rising continuously over 5.5h, post-hypoxia half-life 5.2h | PRIMARY, Eckardt et al. 1989 (PMID 2732171), n=6 |
| 2,3-BPG shift | **within 24h** | PRIMARY, Lenfant et al. 1968 (PMID 5725278) |
| Renal HCO3⁻ compensation | **3-5 days** | Wikipedia "Respiratory alkalosis" (flagged), live-fetched: "3-5 day delay in kidney compensation" |
| Reticulocyte/hematocrit measurable rise | **days to weeks** (9% red-cell-mass rise over a real 28-day exposure) | PRIMARY, Levine & Stray-Gundersen 1997 (PMID 9216951); REUSES `MECHANISM_ERYTHROPOIESIS.md`'s own τ=6d marrow-transit delay, 15d(feedback)/193d(passive) recovery range |
| Full RBC population turnover | **115-120 days** | REUSED verbatim, `MECHANISM_ERYTHROPOIESIS.md` (Mock et al. 2011, PMID 21062290) |

**Forced (Orient-stage) correction, not an assumption**: "EPO is slow" is imprecise. The EPO **signal** is
fast — comparable to, or faster than, the 2,3-BPG shift. What is slow is EPO's downstream hematological
**payload** (new red-cell mass), gated by the marrow-transit delay and RBC-accumulation kinetics this repo's
own `erythropoiesis.py` already certified for a *different* trigger (hemorrhage). This doc re-triggers the
identical mechanism with a hypoxia stimulus instead — a direct, substantive reuse, not a restatement.
**Gate: PASS** — 6 stages, strictly time-ordered from seconds to ~4 months, each independently sourced.

## 7. C4 — dysfunction pole 1: AMS/HAPE/HACE, ascent RATE dominates over absolute altitude

**Claim**: AMS/HAPE/HACE are the maladaptive signature of the fast arms being incomplete or overwhelmed — and
the primary literature directly forces the naive "altitude-alone" adversary to fail.

**Hackett, Rennie & Levine 1976 (PMID 62991)**, n=278, Pheriche, Nepal (4243m, a **fixed** destination
altitude, ruling out the altitude confound within-study): overall AMS incidence **53%**; incidence *increased*
in those who flew to 2800m, climbed fast, and spent fewer acclimatization nights; severity was "highly
correlated with speed of ascent" and unrelated to sex, prior altitude experience, or load carried. Of the 12
severe cases (7 HAPE + 5 HACE), **11 had flown in, 9 had spent only one night** at Pheriche, none were on
acetazolamide. Acetazolamide reduced AMS incidence/severity specifically in the flew-in group — exactly where
the fast-arm-had-no-time-to-act confound was worst.

**Cross-check (Wikipedia "Altitude sickness", flagged)**: AMS ~20% at 2500m, ~40% at 3000m (rapid ascent),
onset "often within 10 hours." Together with Hackett's 53% at 4243m: a **monotonic incidence gradient across
3 altitude bands from 2 independent sources** — consistent with, not contradicting, the ascent-rate finding.

**Gate: PASS** — the primary 1976 cohort directly falsifies the altitude-alone adversary within its own data
(same destination altitude, ascent-rate-dependent incidence/severity, acetazolamide effective specifically in
the fast-ascent subgroup).

## 8. C5 — dysfunction pole 2: chronic mountain sickness, the slow arm overshooting, population-dependent

**Claim**: CMS (Monge's disease) is the symmetric complement to §7 — not the fast arms failing to keep up,
but the slow (erythropoietic) arm chronically overshooting into maladaptive polycythemia, at a rate that is
population/genotype-dependent.

Qinghai-type diagnostic criteria (Wikipedia, flagged, textbook-tier): hemoglobin ≥21 g/dL (♂) / ≥19 g/dL (♀),
hematocrit >65%, SaO2<85%. Prevalence gradient (same source): Ethiopia 0% (3600-4100m), **Tibet 0.91-1.2%**,
Bolivia 8-10%, Cerro de Pasco Peru (4300m) **14.8-18.2%** — an order-of-magnitude difference between Tibetan
and Andean populations under comparable chronic hypoxia.

**Beall CM 2007 (PMID 17494744, PRIMARY, live-verified)**: "Tibetan and Andean high-altitude natives have
adapted differently... large quantitative differences in numerous physiological traits comprising the oxygen
delivery process... the two followed different routes to the same functional outcome... natural selection is
ongoing in the Tibetan population, where women estimated to have genotypes for high oxygen saturation of
hemoglobin... have higher offspring survival." Two populations under near-identical chronic-hypoxia stress
reach adequate O2 delivery via **different weightings of the same physiological parameters** — Tibetans favor
a strategy that raises O2 saturation without driving Hb as high, Andeans favor a strategy that drives Hb
higher — not a single universal "more hemoglobin" axis. The Tibetan strategy's much lower CMS prevalence is
consistent with it overshooting less often.

**Gate: PASS** on the qualitative/comparative claim (primary-source, live-verified); the specific numeric
prevalence percentages and diagnostic thresholds are Wikipedia-sourced (flagged, not independently re-verified
against the primary Qinghai/2013-prevalence papers this session — §10).

## 9. C6 — function pole: live-high-train-low, the deliberate athletic exploitation of the slow arm

**Levine & Stray-Gundersen 1997 (PMID 9216951, PRIMARY, live-verified)**, randomized, n=39 competitive
runners, 4 weeks: "high-low" group (living 2500m / training 1250m, n=13) — VO2max **+5%**, in direct
proportion to red-cell-mass **+9%** (r=0.37, P<0.05); 5000m time improved **13.4±10s**, in direct proportion
to the VO2max change (r=0.65, P<0.01) — improvement seen **only** in the high-low group, not in matched
"high-high" or "low-low" controls. This is athletes deliberately engaging the **same slow erythropoietic arm**
§6-§8 characterize physiologically, while training at low altitude specifically to dodge the fast ventilatory
arm's acute performance penalty. Directly couples to `MECHANISM_VO2MAX_AEROBIC_CAPACITY.md`'s own C2 falsifier
(Buick 1980 PMID 7380690, Bárány 1993 PMID 8482049 — raising O2 content raises VO2max) — the identical
mechanism, reached here by altitude instead of reinfusion or rHuEPO. **Gate: PASS.**

## 10. Symmetric QC — a caught tool error, and what this does NOT prove

**A caught near-miss, disclosed, not hidden**: a `WebFetch` summary of Wikipedia's "Alveolar gas equation"
page reported a sea-level worked-example PAO2 of "107 mmHg." This doc's own independent recomputation, from
the *same* page's own quoted formula and standard constants, gives 99.2-101.3 mmHg — matching the universal
~100 mmHg teaching consensus. The 107 figure was **not reproduced and is not used anywhere** in this doc's
claims — a small-model summarization artifact, machine-cross-checked and rejected rather than propagated
(exactly the discipline this repo's own erythropoiesis/blood-oxygen-transport docs already modeled for their
own citation-drift and hinge-discontinuity catches).

- **0-D, per-stage-steady-state model** — no continuous multi-week ODE coupling all 6 time-constant stages
  into one simulation; each is anchored independently (consistent with `MECHANISM_RESPIRATORY.md`'s own
  disclosed "0-D, steady-state only" scope). §11's `proposed_cell` names this as the natural next step.
- **R=0.85 is West's own stated value**, not independently re-derived; the PAO2=35 reproduction carries a
  small, disclosed 1.2-2.0% residual, not smoothed to an exact match.
- **The Bohr coefficient is not independently pinned down this session** — the identical disclosed gap
  `MECHANISM_ACID_BASE_CO2.md` and `MECHANISM_BLOOD_OXYGEN_TRANSPORT.md` each already found (a **third**
  independent replication of the same gap-shape). Swept across 6 coefficients, direction-gated, not
  magnitude-gated.
- **Altitude arterial pH is a point estimate (7.70)** from West's own reported floor ("over 7.7") — the true
  Bohr-shift benefit is likely understated, not overstated, by this doc's numbers.
- **No SaO2 was directly measured** on the real 1981 summit — West's team measured/calculated PaO2, pH,
  PaCO2; this doc's SaO2 values are *derived* from those measured inputs via the Hill equation, not an
  independently measured quantity.
- **Hackett 1976 is a single 1970s cohort** (n=278, one site) — modern prophylaxis uptake may shift modern
  incidence; the paper's own acetazolamide finding (effective only in the flew-in subgroup) is reported as a
  real, non-generalized result, not smoothed into "acetazolamide always helps."
- **CMS numeric thresholds/prevalence are Wikipedia-sourced** (flagged) — not independently re-extracted from
  the primary Qinghai-score or 2013 prevalence papers this session; the comparative *mechanism* claim (Beall
  2007) IS primary-source, live-verified.
- **The "fourfold ventilation" cross-check** (this doc's own PaCO2-ratio estimate, 3.64-5.33×) assumes
  constant CO2 production across rest/exercise and chamber/field — a real, disclosed simplification; it
  brackets, not exactly reproduces, Sutton 1988's own measured "fourfold."
- **Eckardt 1989 (n=6) and Lenfant 1968 measured lower simulated altitudes (3000-4000m)** than Everest-summit
  (8848m) — extrapolating their timecourses' *magnitude* (not direction/existence) to extreme altitude is a
  disclosed, unverified extrapolation.
- **No permanent `scripts/msk/*.py` file committed** (literature-synthesis-plus-arithmetic genre, matching
  `MECHANISM_VO2MAX_AEROBIC_CAPACITY.md` / `MECHANISM_CORE_STABILITY_IAP.md` precedent) — no graph-edge write
  performed this session.
- **All population-level, not subject-specific** — no individual human twin subject in this repo has directly
  measured altitude blood gases, consistent with the scope every sibling `MECHANISM_*` doc discloses.

## 11. SCENE-EYES aside — reusing an already-certified geometric fact in a new regime

`MECHANISM_BLOOD_OXYGEN_TRANSPORT.md`'s own §9 established that pulse oximetry (SpO2) sits on the O2-Hb
curve's **flat plateau** at normal sea-level PaO2 (~95 mmHg) — a genuine sensor null-space, ~41× less
sensitive than the curve's steep low-PO2 shoulder. At extreme-altitude PaO2 (28-35 Torr, this doc), SaO2 sits
on or near that **same curve's steep shoulder** — meaning pulse oximetry becomes maximally informative exactly
in the regime this doc studies. Not a new claim — a direct reuse of an already-certified geometric fact,
extended to a new regime.

## 12. couples_to

- **`MECHANISM_ERYTHROPOIESIS.md`** — supplies the marrow-transit delay (τ=6d), the 15d/193d recovery
  timescale, and the RBC lifespan (Mock 2011, PMID 21062290) reused verbatim in §6; this doc re-triggers that
  identical mechanism with a hypoxia stimulus instead of that doc's hemorrhage trigger, and supplies the
  altitude/2,3-DPG coupling that doc's own §14 explicitly disclosed as missing.
- **`MECHANISM_RESPIRATORY.md`** — this doc directly closes that doc's own disclosed gap ("NO altitude /
  inspired-O2 coupling"); the ventilation/VO2 machinery there is the sea-level special case of the alveolar
  gas equation this doc extends to low Pb.
- **`MECHANISM_BLOOD_OXYGEN_TRANSPORT.md`** — supplies the Hill/CaO2 constants (§1, §5) and the SCENE-EYES
  null-space fact (§11) reused here; this doc supplies that doc's own disclosed "no altitude coupling" gap.
- **`MECHANISM_ACID_BASE_CO2.md`** — supplies the Henderson-Hasselbalch form and the disclosed-uncertain
  Bohr-coefficient sweep convention (§5, §10), here extended from that doc's acidosis direction to this doc's
  alkalosis direction; a genuinely new regime, not a re-run.
- **`MECHANISM_VO2MAX_AEROBIC_CAPACITY.md`** — supplies the Fick product and the O2-content-raises-VO2max
  citations (Buick 1980, Bárány 1993) this doc's §9 (LHTL) directly extends to a third manipulation type
  (altitude acclimatization, alongside that doc's reinfusion/hyperoxia/rHuEPO manipulations).
- **Function↔dysfunction organizing axis** — §7 (AMS/HAPE/HACE, fast-arm incomplete) and §8 (CMS, slow-arm
  overshoot) are the two symmetric dysfunction poles; §4-§6 (the cascade defending SaO2/PaO2) plus §9 (LHTL)
  are the function pole — mirroring the identical organizing structure `MECHANISM_VO2MAX_AEROBIC_CAPACITY.md`
  already established for its own Q-side/avO2diff-side split.

(Per `docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/§4, promoting these couplings into canonical graph edges in
`data/MECHANISM_ANCHOR_GRAPH.json` requires the separate `mechanism_fold` pipeline — not performed this session,
consistent with every sibling doc's own stated practice.)

## 13. Pre-registered gates — 6/6 claims PASS (+ disclosed tiering, none laundered)

```
C1_PIO2_matches_OEII_calibration_target:              PASS (0.27% diff)
C1_PAO2_reproduces_west1983_own_derived_value:        PASS (1.2-2.0% diff)
C1_PaO2_convergence_two_decorrelated_studies:         PASS (0.0% diff, exact)
C1_sea_level_self_consistency:                        PASS (<1.3% vs 100mmHg anchor)
C2_geometric_proof_SaO2_Hb_independent:               PASS (identical to 4 decimals, Hb 10-24 g/dL)
C2_geometric_CaO2_Hb_dependent_matches_sibling:        PASS (0.7163 = 1.34*SaO2 exactly)
C2_counterfactual_no_hvr_impossible_R_swept:           PASS (negative at R=0.80/0.85/0.90)
C2_void_floor_regime_dependence:                      PASS (huge margin @ sea level, crosses zero @ summit)
C2_bohr_shift_benefit_direction_robust_6_coefficients: PASS (6/6 coefficients, +9 to +27 pts)
C3_time_ordering_strict_seconds_to_4_months:           PASS (6 stages, each independently sourced)
C3_epo_signal_vs_payload_distinction_forced:           PASS (Orient-stage correction, not assumed)
C4_ascent_rate_falsifies_altitude_alone_adversary:     PASS (Hackett 1976, same-altitude cohort)
C4_monotonic_incidence_gradient_2_independent_sources: PASS (20/40/53% across 2500/3000/4243m)
C5_cms_qualitative_comparative_mechanism:              PASS (Beall 2007, primary, live-verified)
C5_cms_numeric_thresholds_prevalence:                  DISCLOSED WIKIPEDIA-TIER (not primary-re-verified)
C6_lhtl_dose_response_quantified:                      PASS (r=0.37 VO2max~mass, r=0.65 perf~VO2max)

OVERALL: 6/6 claims PASS; 1 sub-item explicitly disclosed at a lower (Wikipedia) evidence tier, not
laundered into the same confidence as the primary-source claims; 1 tool-summary artifact ("107 mmHg") caught
via independent recomputation and explicitly rejected, not propagated.
```

## 14. Repro

All arithmetic (alveolar gas equation, Hill/Bohr O2-Hb curve, geometric Hb-independence proof, counterfactual
and void-floor sweeps) was computed this session in self-contained Python (stdlib `math` only), against the
independently-sourced inputs and live-verified anchors tabulated in §3-§9 and in
`docs/MECHANISM_ALTITUDE_ACCLIMATIZATION_evidence.json`. No OpenSim call, no `.osim`/`.sto` file, no
subject-specific twin data — population-level literature synthesis plus arithmetic cross-checks, consistent
with the `docs/MECHANISM_VO2MAX_AEROBIC_CAPACITY.md` / `docs/MECHANISM_CORE_STABILITY_IAP.md` precedent (no
permanent `scripts/msk/*.py` file created for this subagent-scoped, two-file deliverable).

**Paths**: this doc `docs/MECHANISM_ALTITUDE_ACCLIMATIZATION.md`; evidence
`docs/MECHANISM_ALTITUDE_ACCLIMATIZATION_evidence.json`.
