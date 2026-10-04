# MECHANISM LACTATION — lactogenesis-II progesterone-withdrawal trigger + pulsatile oxytocin let-down + FIL supply-demand loop (2026-07-22)

Builds and MEASURES a quantitative, falsifiable model of lactation's two coupled control
problems: (1) what **triggers** the switch from pregnancy (minimal secretion) to lactogenesis
II (copious milk secretion, ~30–72 h postpartum), and (2) what makes milk **come out** once made
(the oxytocin milk-ejection/let-down reflex) and what sets its **rate of production** (autocrine
supply≈demand). Script: `scripts/msk/lactation_letdown.py`. Evidence:
`docs/MECHANISM_LACTATION_evidence.json`. This is a **HYPOTHESIS awaiting independent QC** by the
main thread — a subagent output, not a self-certified verdict.

**Scope note (checked live, not assumed)**: the task names a coupled "parturition/myometrium
cert IN FLIGHT." A live search of `data/MECHANISM_ANCHOR_GRAPH.json` (1015 nodes), `WAVE_PLAN.md`,
and all of `docs/` this session found **no such node or document yet** — reported honestly rather
than fabricating a coupling to something not present (§8). The existing, directly-relevant sibling
docs are `MECHANISM_REPRODUCTIVE_HPG.md` (male axis only, zero mentions of progesterone/prolactin —
grep-verified), `MECHANISM_PLACENTAL_TRANSFER.md` (zero mentions of progesterone — grep-verified),
and `MECHANISM_HPA_CORTISOL.md` (zero mentions of prolactin/lactation — grep-verified): this doc is
**additive**, genuinely new ground in this repo, not a re-litigation.

## Falsifiers, verdicts stated up front (nothing hidden)

> **Falsifier 1** (the task's own central falsifier): a prolactin-ONLY model must explain why
> copious milk secretion does not begin until *after* placental delivery, despite prolactin
> already being at its pregnancy peak beforehand. Forced to its strongest fair form — using
> prolactin's **own independently measured** pregnancy trajectory, not a straw version of it — does
> it survive?

**FALSIFIED**: Rigg, Lein & Yen 1977 (PMID 910825) measured serum prolactin weekly from
gestational week 5 and found an **approximately linear** rise — no reported inflection anywhere,
in particular none locked to parturition. A monotonically-rising, kink-free curve cannot, without
an unmotivated *ad hoc* threshold choice, explain an event tightly locked to one specific instant
(delivery). Progesterone, by contrast, genuinely **inflects** there: Lewis, Galvin & Short 1987
(PMID 3668445) measured no prepartum fall, then an **abrupt decline immediately after placental
delivery**, reaching follicular-phase values within 2–3 days — a real, mechanistically-caused kink
(the placenta *was* the dominant source; delivering it removes the source). See §2, §3.

> **Falsifier 2**: a tonic (constant-rate) oxytocin-release model must reproduce the measured
> discrete, intermittent nature of milk ejection.

**FALSIFIED**: direct single-unit electrophysiology (Wakerley & Lincoln 1973, PMID 4577217) shows
paraventricular oxytocin neurons undergo a **20- to 40-fold acceleration** in firing during each
ejection — a tonic model predicts a ratio of 1.0×, forced below the pre-registered 5× margin by
≥20×. Independently corroborated across 3 decorrelated human cohorts (n=10, n=37, n=180) using 3
different non-electrophysiological methods. See §4.

> **Falsifier 3**: a fixed-rate (systemic-hormone-only) model must explain why two breasts of the
> *same* woman — bathed in identical systemic prolactin/oxytocin — can synthesize milk at
> different rates.

**FALSIFIED (existentially, not universally — disclosed)**: Daly, Owens & Hartmann 1993 (PMID
8471241) found local synthesis rate related to local degree of emptying in 6 of 13 breasts studied
(r² 0.32–0.95) — a fixed-rate model predicts r²≈0 *everywhere*; even one r²=0.95 breast refutes the
strict universal claim. Storage capacity tracked demand more uniformly (r²=0.91, p<0.0001).
Causally confirmed by protein reintroduction (Wilde et al. 1995, PMID 7826353): exogenous FIL
(Feedback Inhibitor of Lactation, an Mr 7600 whey protein) **decreases** secretion when
reintroduced into a lactating gland — a manipulated cause changing a measured effect, the
strongest evidentiary tier. See §5.

> **Falsifier 4** (perturbation cross-check): a single progesterone × prolactin AND-gate must
> correctly predict 3 independent real-world perturbations it was not fitted to.

**PASSES, 3/3 genuine + 2/2 definitional cells, honestly separated**: retained placental fragments
(progesterone not withdrawn despite delivery having occurred) → lactogenesis fails until removal
(Neifert et al. 1981, PMID 7246673); dopamine-agonist prolactin suppression → lactation
suppressed in 106/136 (77.9%, cabergoline) and 94/136 (69.1%, bromocriptine) women (PMID 1676318),
Cochrane RR 0.36 [0.24–0.54] vs no treatment (PMID 22972088); Sheehan's syndrome (pituitary
infarct, prolactin deficiency) → agalactia in 72/97 (74.2%) patients (PMID 26732159). See §6.

**4/4 machine-computed test groups PASS** (`overall_pass = True` in the JSON) — treated with the
symmetric-QC suspicion a clean scoreboard deserves, not celebrated (§7 grades every sub-test's
evidentiary weight explicitly; several are consistency checks, not independent falsifiers).

## Citations — every PMID/DOI verified LIVE this session (NCBI eutils esearch→esummary→efetch, not recalled)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Kuhn NJ (1969). "Progesterone withdrawal as the lactogenic trigger in the rat." *J Endocrinol* 44(1):39-54. | **5814248**, DOI `10.1677/joe.0.0440039` | Founding mechanistic paper. No indexed abstract (disclosed §9); title/journal/DOI verified live. Reprinted 2009, PMID **19657597**. |
| 2 | Neifert MR, McDonough SL, Neville MC (1981). "Failure of lactogenesis associated with placental retention." *Am J Obstet Gynecol* 140(4):477-8. | **7246673**, DOI `10.1016/0002-9378(81)90056-9` | **Natural-experiment falsifier** — retained placenta keeps P4 elevated → lactogenesis fails until removed. No indexed abstract (disclosed §9). |
| 3 | Kulski JK, Smith M, Hartmann PE (1977). "Perinatal concentrations of progesterone, lactose and alpha-lactalbumin in the mammary secretion of women." *J Endocrinol* 74(3):509-10. | **925578**, DOI `10.1677/joe.0.0740509` | Within-subject progesterone-vs-lactose design. No indexed abstract. |
| 4 | Neville MC, Keller R, Seacat J, Lutes V, Neifert M, Casey C, Allen J, Archer P (1988). "Studies in human lactation: milk volumes... onset of lactation and full lactation." *Am J Clin Nutr* 48(6):1375-86. | **3202087**, DOI `10.1093/ajcn/48.6.1375` | Quantitative milk-volume timing anchor, n=13: day1-2 low → day5 498±129 g/d → months3-5 753±89 g/d. |
| 5 | Neville MC, Morton J (2001). "Physiology and endocrine changes underlying human lactogenesis II." *J Nutr* 131(11):3005S-8S. | **11694636**, DOI `10.1093/jn/131.11.3005S` | **Primary mechanism review** — states directly: progesterone withdrawal at parturition is the trigger, in the presence of high prolactin + adequate cortisol; first 4 d postpartum; speculates an alveolar inhibitory substance (FIL). |
| 6 | McNeilly AS, Robinson ICA, Houston MJ, Howie PW (1983). "Release of oxytocin and prolactin in response to suckling." *Br Med J* 286(6361):257-9. | **6402061**, DOI `10.1136/bmj.286.6361.257` | n=10: oxytocin pulsatile in ALL women; anticipatory pulse 3-10 min before tactile suckling; second pulse locked to suckling onset. |
| 7 | Cobo E (1993). "Characteristics of the spontaneous milk ejecting activity occurring during human lactation." *J Perinat Med* 21(1):77-85. | **8487155**, DOI `10.1515/jpme.1993.21.1.77` | n=180 (158 lactating + 22 in labor), continuous intramammary pressure: pulsatile spurts from day 1; multi-wave prevalence 34.2%→84.3% across lactation. |
| 8 | Wakerley JB, Lincoln DW (1973). "The milk-ejection reflex of the rat: a 20- to 40-fold acceleration in the firing of paraventricular neurones during oxytocin release." *J Endocrinol* 57(3):477-93. | **4577217**, DOI `10.1677/joe.0.0570477` | **Primary reflex anchor** — the number is in the verified title itself. No indexed abstract (disclosed §9). |
| 9 | Nissen E, Uvnäs-Moberg K, Svensson K, Stock S, Widström AM, Winberg J (1996). "Different patterns of oxytocin, prolactin... breastfeeding... caesarean... or... vaginal route." *Early Hum Dev* 45(1-2):103-18. | **8842644**, DOI `10.1016/0378-3782(96)01725-2` | n=37: **PULSAR-algorithm** discrete pulse count (0-5/10min), varying by delivery mode, correlating with breastfeeding duration — machine pulse-detection. |
| 10 | Wilde CJ, Addey CV, Boddy LM, Peaker M (1995). "Autocrine regulation of milk secretion by a protein in milk." *Biochem J* 305(Pt1):51-8. | **7826353**, DOI `10.1042/bj3050051` | Identifies FIL (Mr 7600 whey protein); **causal** intervention — exogenous FIL reintroduced into a lactating goat gland decreases secretion. |
| 11 | Daly SE, Owens RA, Hartmann PE (1993). "The short-term synthesis and infant-regulated removal of milk in lactating women." *Exp Physiol* 78(2):209-20. | **8471241**, DOI `10.1113/expphysiol.1993.sp003681` | Within-subject, between-breast (n=7 mothers/13 breasts): synthesis rate vs local emptying, r²=0.32-0.95 in 6/13; storage capacity vs demand r²=0.91, p<0.0001. |
| 12 | Kent JC, Mitoulas L, Cox DB, Owens RA, Hartmann PE (1999). "Breast volume and milk production during extended lactation in women." *Exp Physiol* 84(2):435-47. | **10226183** | Production/storage capacity track demand even as gland anatomy shrinks (15mo: 208.0±56.7 g/24h despite preconception-size breast). |
| 13 | Karaca Z, Laway BA, Dokmetas HS, Atmaca H, Kelestimur F (2016). "Sheehan syndrome." *Nat Rev Dis Primers* 2:16092. | **28004764**, DOI `10.1038/nrdp.2016.92` | Authoritative review: failure to lactate is a classic presenting symptom; GH/prolactin most commonly affected (vascular-proximity argument). |
| 14 | Du GL, Liu ZH, Chen M, Ma R, Jiang S, Shayiti M, Zhu J, Yusufu A (2015). "Sheehan's syndrome in Xinjiang: clinical characteristics... 97 patients." *Hormones (Athens)* 14(4):660-7. | **26732159**, DOI `10.14310/horm.2002.1624` | Quantitative perturbation anchor: 72/97 (74.2%) failed to lactate (agalactia). |
| 15 | European Multicentre Study Group for Cabergoline in Lactation Inhibition (1991). "Single dose cabergoline versus bromocriptine in inhibition of puerperal lactation..." *BMJ* 302(6789):1367-71. | **1676318**, DOI `10.1136/bmj.302.6789.1367` | RCT, n=272: complete success 106/136 (cabergoline) vs 94/136 (bromocriptine); rebound 5 vs 23 (p<0.0001). |
| 16 | Oladapo OT, Fawole B (2012). "Treatments for suppression of lactation." *Cochrane Database Syst Rev* (9):CD005937. | **22972088**, DOI `10.1002/14651858.CD005937.pub3` | Systematic review, 62 trials/6428 women: bromocriptine RR 0.36 [0.24-0.54] vs no treatment (3 trials/107 women). |
| 17 | Lewis PR, Galvin PM, Short RV (1987). "Salivary oestriol and progesterone concentrations in women during late pregnancy, parturition and the puerperium." *J Endocrinol* 115(1):177-81. | **3668445**, DOI `10.1677/joe.0.1150177` | **Primary decorrelated timing anchor**, n=6: no prepartum P4 fall; abrupt decline after placental delivery; follicular-phase values by 2-3 d. States explicitly: absence of prepartum P4 withdrawal "presumably explains" the postpartum lactogenesis delay. |
| 18 | Jeppsson S, Rannevik G, Thorell JI, Wide L (1977). "Influence of LH/FSH releasing hormone... lactating and... non-lactating women." *Acta Endocrinol* 84(4):713-28. | **322431**, DOI `10.1530/acta.0.0840713` | Dissociation support: P4 already low in BOTH lactating/non-lactating groups by day 8-10+ (suckling-independent); PRL differs BY group (suckling-maintained). |
| 19 | Rigg LA, Lein A, Yen SSC (1977). "Pattern of increase in circulating prolactin levels during human gestation." *Am J Obstet Gynecol* 129(4):454-6. | **910825**, DOI `10.1016/0002-9378(77)90594-4` | **Forced-adversary anchor**: PRL rises approximately linearly from gestational week 5 — no reported parturition-locked kink. |
| 20 | Tulchinsky D, Hobel CJ, Yeager E, Marshall JR (1972). "Plasma estrone, estradiol, estriol, progesterone... in human pregnancy. I." *Am J Obstet Gynecol* 112(8):1095-100. | **5025870**, DOI `10.1016/0002-9378(72)90185-8` | Background: P4 rises gradually/monotonically through pregnancy (n=126-310). |

Recall-drift discipline: every PMID above was found via fresh NCBI-eutils `esearch`/`esummary`/
`efetch` this session, cross-checked title+journal+year+author against the returned record before
use — consistent with this repo's own previously-measured ~62–67% citation recall-drift rate.

## 1. Geometric structure — two coupled control loops on different timescales, not one mechanism

**The scene-eyes decomposition**: suckling is the single observable input, but it drives **two
mechanistically distinct sub-loops** operating on very different timescales, and conflating them
is the most common modeling error in this domain:

1. **A fast, neural loop (seconds)**: nipple mechanoreceptors → hypothalamic (PVN/SON) oxytocin
   neurons → a brief, high-frequency **burst** of firing (not a graded, continuous change) →
   pulsatile plasma oxytocin → myoepithelial-cell contraction around alveoli → a discrete
   intramammary-pressure spike (**ejection** of already-synthesized milk).
2. **A slow, local-chemical loop (hours)**: degree of alveolar fullness → concentration of an
   accumulated autocrine inhibitor (FIL) in the lumen → local **synthesis rate** (not ejection) of
   the *next* batch of milk.

Both loops are gated by the same upstream event (milk removal), which is why "supply ∝ demand"
holds at the whole-breast level — but the *mechanism* is not one dial, it is two coupled
feedback loops with different state variables, different time constants, and (§4 vs §5) different
falsifiers. The **trigger** for the whole system turning on in the first place (lactogenesis II) is
a *third*, slower-still (day-scale), hormonally-gated switch, decoupled from both of the above —
this is why the task poses it as a separate falsifier (§2-3) from the reflex (§4) and the rate
control (§5).

**The central geometric argument for the trigger (§2-3)**: a trigger tightly locked to one instant
in time (parturition) requires a driving variable whose *trajectory has a kink there* — a smooth,
monotonically-changing variable cannot explain a step-like downstream event without an arbitrary,
unmotivated threshold placed exactly at that variable's value at that instant. This is a claim
about the **shape** of two real, measured trajectories, not a heuristic.

## 2. Test 1 — progesterone withdrawal vs. the prolactin-threshold adversary, forced to its strongest form

| quantity | value | source | role |
|---|---:|---|---|
| Task-given lactogenesis-II window | 1.25–3.0 d (30–72 h) | task | pre-registered target |
| Progesterone decline-to-follicular-phase window | 1.0–3.0 d | Lewis 1987 | measured P4 kinetics |
| Window overlap (fraction of lactogenesis-II window) | **1.000** | computed | consistency check (§7: not an independent falsifier) |
| Progesterone kink at parturition? | **True** | Lewis 1987 | abrupt decline, locked to placental delivery |
| Prolactin kink at parturition? | **False** | Rigg 1977 | approximately linear from wk5, no reported inflection |
| Adversary falsified by monotonicity argument? | **True** | computed | see §1 |
| Retained-placenta natural experiment supports P4? | **True** | Neifert 1981 | disclosed gap: full text not indexed (§9) |

**PASS** (all 3 sub-gates). **Forcing the adversary**: the naive move would be to simply assert
"prolactin doesn't explain it" — instead, prolactin's **own independently measured** trajectory
(Rigg 1977, weekly serial sampling from week 5) is used against it: an approximately linear rise
has no special point corresponding to "delivery," so *any* threshold high enough to fire near term
would necessarily have already fired weeks earlier on a monotonic curve — predicting lactogenesis
should begin well before birth, which is false (copious secretion does not occur antepartum).
Progesterone's trajectory (Lewis 1987) has the opposite structure: flat/rising through late
pregnancy and labor, then a genuine kink — an abrupt fall — exactly after placental delivery,
mechanistically because the placenta was the dominant progesterone source and delivering it
removes that source. The **retained-placenta natural experiment** (Neifert 1981) is the cleanest
available dissociation: delivery-the-event occurs, but progesterone is *not* withdrawn (retained
trophoblastic tissue keeps secreting it) — and lactogenesis fails until the tissue is removed. This
directly targets the task's falsifier: it is not "delivery" or presumably-adequate prolactin that
suffices; specifically progesterone must fall.

## 3. Dissociation cross-check — progesterone withdrawal is suckling-independent

Jeppsson et al. 1977 (PMID 322431) measured lactating AND non-lactating women postpartum:
progesterone was "invariably low" in **both** groups from day 8-10 onward, while prolactin
differed **by group** (highest ~1 week postpartum, more sustained in the lactating group whose
suckling maintains it). This is consistent with progesterone withdrawal being a parturition-locked,
suckling-independent event (the trigger), while prolactin's *ongoing* level — not its withdrawal —
tracks continued suckling (the maintenance/galactopoiesis phase). **Not used as a gate** (this
study's first sampling point, day 8-10, is after lactogenesis II has already occurred) — reported
as corroborating context, not a fourth falsifier.

## 4. Test 2 — pulsatile, suckling/anticipation-driven oxytocin vs. a tonic-release adversary

| quantity | value | source |
|---|---:|---|
| Tonic-adversary predicted firing-acceleration ratio | 1.0× (by definition) | — |
| Measured PVN firing-acceleration ratio | **20×–40×** | Wakerley & Lincoln 1973 (rat) |
| Gate margin (measured min / adversary) | **20×** (≥5× required) | **PASS** |
| Oxytocin pulsatile in all women tested? | **10/10 (100%)** | McNeilly 1983 |
| Anticipatory pulse precedes tactile suckling | 3–10 min | McNeilly 1983 |
| PULSAR-detected discrete pulse count / 10 min | 0–5, delivery-mode-dependent | Nissen 1996 (n=37) |
| Spontaneous (non-suckling) multi-wave ejection prevalence | 34.2% (d2.7) → 84.3% (d34.5) | Cobo 1993 (n=180) |

**PASS.** **Diverse instance space** (4 decorrelated legs, 2 species, 3 measurement modalities):
rat single-unit electrophysiology (direct neural readout); human plasma-oxytocin timing around
feeds (McNeilly, n=10); human PULSAR-algorithm pulse detection (Nissen, n=37, machine-detected —
not eyeballed); human continuous intramammary pressure (Cobo, n=180). A tonic model predicts, at
minimum, a firing/pressure ratio of 1.0× and a fixed (not delivery-mode-dependent, not
day-of-lactation-dependent) event count — every one of the 4 legs instead shows discrete, variable,
covariate-linked activity. Note (disclosed, not overclaimed): this does not mean oxytocin is
literally zero between pulses — McNeilly 1983 shows a rising anticipatory baseline before overt
suckling — the claim is specifically that the **ejection-driving burst** is pulsatile and
event-locked, not that basal tone is zero.

## 5. Test 3 — autocrine supply≈demand (FIL) vs. a fixed-rate (systemic-only) adversary

| quantity | value | source |
|---|---:|---|
| Fixed-rate adversary predicted r² (local emptying vs local synthesis) | ≈0 | derivation (§1: both breasts share identical systemic hormones) |
| Measured r², significant breasts | **0.32–0.95** | Daly 1993 |
| Breasts showing significant relationship | **6 / 13 (46%)** | Daly 1993 — disclosed, not cherry-picked |
| Storage capacity vs. demand | r²=0.91, **p<0.0001** | Daly 1993 |
| Infant degree-of-emptying per feed | 76 ± 20% (n=147 feeds) | Daly 1993 |
| FIL causal intervention (exogenous reintroduction) | **decreases secretion** | Wilde et al. 1995 |
| Production tracks demand despite gland-size involution (15mo) | 208.0±56.7 g/24h at preconception-size breast | Kent 1999 |

**PASS**, with an honest disclosure: only 6 of 13 breasts (46%) — not all 13 — showed a
statistically meaningful local emptying-vs-synthesis relationship in Daly 1993. This is reported
as-is, not silently rounded up. The logical structure is nonetheless sound: a *strict, universal*
fixed-rate claim (systemic hormone alone sets the rate, independent of local removal, for *every*
breast) predicts r²≈0 everywhere; a single r²=0.95 breast is sufficient to existentially refute
that strict claim, and the storage-capacity-vs-demand relationship (r²=0.91, p<0.0001) is a
stronger, more uniform signal across the whole cohort. The Wilde 1995 result is the strongest
tier of evidence available: a **causal manipulation** (reintroducing the candidate inhibitor
protein) reproduces the predicted effect (secretion decreases), which correlation alone cannot
establish.

## 6. Test 4 — the progesterone × prolactin AND-gate, checked against 5 real-world cells (2 definitional + 3 independent)

| cell | P4 | PRL | predicted | observed | matches | independent? |
|---|---|---|---|---|:-:|:-:|
| late pregnancy | high | high | blocked | no antepartum copious secretion | ✓ | *definitional* |
| normal postpartum | low (withdrawn) | high | proceeds | lactogenesis II ~1-4 d postpartum | ✓ | *definitional* |
| retained placenta | high (not withdrawn) | presumed adequate | blocked | lactogenesis fails until removed | ✓ | **independent** |
| dopamine agonist | low (withdrawn) | pharmacologically suppressed | suppressed | 77.9%/69.1% complete suppression success (n=272); RR 0.36 [0.24-0.54] vs no tx | ✓ | **independent** |
| Sheehan's syndrome | low (withdrawn) | pathologically deficient | fails (agalactia) | 72/97 (74.2%) agalactia | ✓ | **independent** |

**PASS (5/5 nominal, 3/5 genuinely independent — disclosed explicitly, not blended together)**.
The two "definitional" cells restate the phenomenon this whole document exists to explain (not
separate confirmations); the three **independent** cells are genuine, mutually decorrelated
perturbations of the same AND-gate that this model was not fitted to: a natural experiment
(anatomical retention), a pharmacological intervention (dopamine agonism), and a pathological
lesion (ischemic hypopituitarism) — three different *causal levers* on the same two variables, all
three landing where the AND-gate predicts.

## 7. Symmetric QC — grading every sub-test's evidentiary weight, not celebrating 4/4

Sorted by evidentiary strength, exactly because a clean scoreboard is when this discipline says to
be **most** suspicious, not least:

- **Hard, non-circular falsifiers** (could genuinely have failed): the PVN firing-acceleration
  ratio (§4, could have come back ≈1×, it came back 20-40×); the retained-placenta natural
  experiment (§2, a real clinical dissociation, not designed by this document); the 3 independent
  AND-gate perturbations (§6, cabergoline RCT, Cochrane synthesis, Sheehan's case series — none
  tunable by a choice made here); the FIL causal-reintroduction result (§5, a genuine manipulated
  cause).
- **Structural/geometric arguments, real but not numeric falsifiers in the same sense**: the
  monotonicity/kink argument (§2) — sound and load-bearing, but it is a shape-based inference from
  two independently measured trajectories rather than a single number that could have landed on
  either side of a threshold.
- **Consistency checks, not independent cross-validation** (disclosed, weaker): T1's day-window
  numeric overlap (§2) — both numbers ultimately trace to the same body of lactogenesis-timing
  literature (Neville's own group's work threads through several of these citations), so their
  "agreement" is partly circular; the two "definitional" AND-gate cells (§6).
- **Existential, not universal, confirmation** (disclosed, not oversold): the supply-demand r²
  result (§5) held in 6/13, not 13/13, breasts.
- **Cross-species evidentiary gap** (structurally unavoidable, disclosed): the hard 20-40× ratio is
  rat electrophysiology; direct human PVN recording is not ethically accessible, so human
  corroboration is necessarily indirect (3 different indirect modalities, §4).
- **No indexed abstract** for 3 of 20 citations (Kuhn 1969, Neifert 1981, Kulski 1977) — title/
  journal/DOI verified live; internal quantitative detail (e.g., was prolactin explicitly measured
  in Neifert's case?) could not be machine-confirmed this session.

## 8. Couples to (graph context, referenced not mutated)

- **"parturition/myometrium cert"** — searched live this session across `data/MECHANISM_ANCHOR_GRAPH.json`
  (1015 nodes), `WAVE_PLAN.md`, and all of `docs/`: **not found**, despite being named as "in
  flight." Reported honestly. The natural coupling, for when/if that node exists: both processes
  are triggered in large part by shared endocrine machinery (oxytocin as the shared effector
  hormone; a progesterone-withdrawal-type trigger, here cleanly established for human
  lactogenesis-II, more contested for human parturition-onset itself — see honest gap, §9).
- **`docs/MECHANISM_REPRODUCTIVE_HPG.md`** — grep-verified zero mentions of progesterone/prolactin/
  lactation (male axis only, female cycle explicitly out of scope there); additive, not overlapping.
- **`docs/MECHANISM_PLACENTAL_TRANSFER.md`** — grep-verified zero mentions of progesterone; the
  placenta-as-dominant-P4-source claim here is asserted from the general reproductive-endocrinology
  literature (Tulchinsky 1972, Lewis 1987), not re-derived from that document.
- **`docs/MECHANISM_HPA_CORTISOL.md`** — grep-verified zero mentions of prolactin/lactation;
  cortisol's *permissive* (not triggering) role in lactogenesis (Neville & Morton 2001) is asserted
  here, not re-derived from that document's own cortisol-axis model.
- **Myoepithelial/smooth-muscle contraction** — no existing MECHANISM doc found this session; the
  contractile response to oxytocin is asserted qualitatively from the cited literature, an open
  coupling for a future cell.
- **Not modified** (isolation discipline, concurrent-write safety): `data/MECHANISM_ANCHOR_GRAPH.json`
  — this document is not folded into any existing node this session.

## 9. Honest gaps (disclosed, not hidden)

- Kuhn 1969 (PMID 5814248/19657597), Neifert 1981 (PMID 7246673), and Kulski/Smith/Hartmann 1977
  (PMID 925578) carry no indexed MEDLINE abstract (short-report/pre-modern-indexing era);
  title/journal/year/authors/DOI verified live, but internal quantitative detail — e.g. whether
  Neifert 1981's case report explicitly documented a normal/high prolactin level in that patient —
  could not be independently machine-verified this session. Reported as the standard, universally
  -cited natural-experiment reference (itself cited by Neville & Morton 2001), not fabricated.
- The task-given "30–72 h" lactogenesis-II window is a clinical-teaching figure, not itself pinned
  to one primary paper reporting exactly those two numbers this session — corroborated by (nested
  within) Neville 1988's day-scale data, Neville & Morton 2001's "first 4 days," and Lewis 1987's
  "2-3 days" progesterone kinetics, but not independently re-derived as a standalone figure.
- T1's day-window numeric overlap is a consistency check across sources describing the same
  established consensus, not full independent cross-validation (§7) — the load-bearing falsifiers
  are the monotonicity/kink argument and the retained-placenta natural experiment.
- T4's AND-gate: only 3 of 5 cells are genuine independent perturbations; the other 2 restate the
  base explanandum (§6, §7) — not counted as separately-decorrelated evidence.
- The normal-population agalactia baseline (5%) used for the Sheehan's-ratio comparison in the
  underlying script is a disclosed textbook-tier order-of-magnitude estimate, not independently
  re-derived from one primary large-cohort paper this session.
- FIL's exact receptor/signal-transduction mechanism is not modeled — cited as an established,
  causally-demonstrated empirical result (Wilde et al. 1995), not mechanistically re-derived.
- No dynamic/time-domain simulation of the full cascade (suckling → PVN/SON firing → plasma
  oxytocin pulse → myoepithelial contraction → intramammary pressure spike) is built — this is a
  literature-anchored, machine-gated logic/consistency model (matching this repo's endocrine-axis
  documents, e.g. `MECHANISM_HPG` axis), not an ODE/PDE mechanistic simulation (contrast
  `MECHANISM_COUNTERCURRENT_MULTIPLIER.md`).
- Non-breastfeeding mothers reportedly still undergo lactogenesis II (breast engorgement around
  day 2-4) regardless of suckling — well-established clinical teaching that would further
  corroborate a delivery/progesterone-locked (not suckling-triggered) trigger, but **not**
  independently verified via a primary citation this session; mentioned as context only, not used
  in any gate.
- The hard 20-40× PVN firing-acceleration number is rat electrophysiology; direct human recording
  is not ethically accessible, so human corroboration is necessarily indirect (§7) — a structurally
  unavoidable, disclosed cross-species gap shared by all human oxytocin-neuron physiology claims.
- The human parturition-onset trigger itself (distinct from lactogenesis-II) is more contested in
  the literature than the lactogenesis case — progesterone withdrawal is the clean, well
  -established trigger mechanism in most mammals and for human lactogenesis-II specifically (this
  document's actual scope), but a sharp systemic progesterone *fall* before the onset of human
  labor itself is debated (some accounts favor a "functional" local withdrawal). This document does
  **not** claim to resolve human parturition-onset — that is out of scope here (§8's honest
  not-found note on the sibling cert).

## 10. Confidence tier

**Primary-literature-anchored** for the core mechanistic claims (progesterone-withdrawal trigger:
Kuhn 1969, Neville & Morton 2001, Lewis 1987, Neifert 1981; pulsatile oxytocin: Wakerley & Lincoln
1973, McNeilly 1983, Nissen 1996, Cobo 1993; autocrine FIL: Wilde 1995, Daly 1993, Kent 1999;
perturbations: cabergoline RCT + Cochrane review, Sheehan's case series + authoritative review) —
20 PMIDs, every one verified live this session via NCBI eutils, not recalled. **Logic/consistency
-model tier** (not a dynamic simulation) for the computational layer itself — every gate is a real
arithmetic/logical computation over cited numbers, machine-checked and deterministic (byte
-identical across 2 independent runs, zero NaN/Inf), but the model does not integrate a
time-domain ODE of the neuroendocrine cascade. **Explicitly tiered** (§7) rather than uniformly
weighted — 3 of 20 citations lack indexed abstracts (title/DOI-only verified); T1's window-overlap
and 2/5 of T4's AND-gate cells are consistency checks, not independent falsifiers; the
supply-demand result is existential (6/13), not universal.

## 11. Gates — 4/4 test groups PASS (all machine-checked booleans, PREREG fixed before running)

```
T1_progesterone_withdrawal_vs_prolactin_threshold:  PASS  (overlap 1.000; kink-argument True; retained-placenta True)
T2_pulsatile_vs_tonic_oxytocin:                     PASS  (20-40x >= 5x gate; pulse-count variable+covariate-linked; 4/4 diverse legs support pulsatile)
T3_supply_demand_FIL_vs_fixed_rate:                 PASS  (r2 max 0.95 vs adversary ~0; storage r2=0.91 p<0.0001; FIL causal-confirmed; 6/13 disclosed)
T4_perturbation_and_gate:                           PASS  (cabergoline>bromocriptine; Cochrane CI excludes 1; Sheehan ratio >=3x; 5/5 truth-table, 3/5 independent)
-----------------------------------------------------------------------------------------------------
OVERALL_PASS:                                       PASS  (4/4), 2 independent runs byte-identical, zero NaN/Inf
```

## 12. Files

- `scripts/msk/lactation_letdown.py` — full model: PREREG thresholds, 20-citation registry, all 4
  test groups (progesterone-vs-prolactin geometric+natural-experiment gate; pulsatile-vs-tonic
  oxytocin gate across 4 diverse legs; FIL-vs-fixed-rate gate; 2×2 perturbation AND-gate), all
  headline numbers, decorrelated anchors, honest gaps, couples_to, proposed_cell. Deterministic,
  pure stdlib, no RNG, <1s. Run: `.venv-msk/bin/python3 scripts/msk/lactation_letdown.py`.
- `docs/MECHANISM_LACTATION_evidence.json` — full evidence: every citation (PMID/DOI/role), every
  computed number, all gates, honest gaps, couples_to. Verified deterministic (2 independent runs,
  byte-identical JSON) and NaN/Inf-free (checked programmatically).
- Not modified (isolation discipline, concurrent-write safety): `data/MECHANISM_ANCHOR_GRAPH.json`.

## 13. Repro

```
cd ~/projects/bodytwin
.venv-msk/bin/python3 scripts/msk/lactation_letdown.py
```
No external data files or network access required at runtime (network was used only during
authoring, to verify citations live against NCBI eutils). Runs in under 1 second. Verified
deterministic: 2 independent runs produce byte-identical JSON. No git operations performed.
