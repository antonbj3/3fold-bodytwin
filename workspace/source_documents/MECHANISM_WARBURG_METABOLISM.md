# MECHANISM WARBURG EFFECT / CANCER METABOLIC REPROGRAMMING — a certified model of aerobic glycolysis, the ATP-yield-vs-flux tradeoff, and WHY (2026-07-22)

Builds and MEASURES the Warburg effect (aerobic glycolysis: tumor cells favor glycolysis + lactate
secretion even with O2 present) as a CERTIFIED model: the ATP-yield-vs-flux tradeoff (glycolysis 2
ATP/glucose, fast, no O2 needed; OXPHOS 30-38 ATP/glucose historically, 33.45 modern precise estimate,
~15-17x higher yield, but needs bulky mitochondrial machinery), a geometric solvent-capacity/molecular-
crowding model that DERIVES why a
demand-dependent shift toward glycolysis is optimal, and literature-measured flux/ATP-fraction
numbers — testing explicitly whether the model and the data refute Warburg's own 1956 "damaged
respiration" hypothesis (PMID 13298683), which is measurably FALSE.

**RESOLVES** (informs; does NOT edit — isolation rule "touch only files you create") the pre-existing
SEED-DESIGN node `ONCO-WARBURG-METABOLISM` already in `data/MECHANISM_ANCHOR_GRAPH.json` (998 nodes).
That node's 4 core datapoints were marked `"fetched_live": true` by a prior agent — per this repo's
own discipline (a prior agent's claim is a hypothesis, not truth, until independently re-checked) ALL
FOUR were independently re-verified live this session via fresh NCBI efetch and matched EXACTLY (a
genuine case where the prior claim held up completely — reported here as a real, checked finding, not
assumed).

## Headline (nothing hidden, both halves reported)

**On the geometric mechanism**: YES, forced through two OODA-driven bug-fixes (see §2) — a
2-constraint linear resource-allocation model (glucose-uptake capacity × solvent/molecular-crowding
capacity) reproduces the Vazquez/Shlomi/Molenaar-described demand-dependent threshold switch from
pure OXPHOS to a gradually-more-glycolytic mixed strategy, cross-checked 3 independent computational
ways (hand-derived vertex enumeration / brute-force grid / scipy LP solver) to <1e-6 relative error
everywhere, and the model makes and PASSES its OWN derived, falsifiable, bidirectional prediction: a
genuine mixed strategy is geometrically optimal if and only if κ_OXPHOS/κ_glycolysis (the relative
crowding cost) exceeds Y_OXPHOS/Y_glycolysis ≈ 16.7-16.9 (25/25 swept ratios correct, both directions).

**On "broken mitochondria"**: NO — measurably refuted 3 independent ways: a causal knockdown+rescue
experiment (Fantin 2006), a real human in-vivo 9-patient intraoperative ¹³C-flux study (Hensley
2016), and Koppenol/Bounds/Dang 2011's direct historical correction of Warburg's own claim.

**On elevated-but-not-exclusive glycolytic ATP fraction**: YES, tumor lines show a measured
glycolytic-ATP-fraction range of 1-64% (up to 66% same-O2, genotype-driven) against a normal-cell
baseline of ~30% — genuinely elevated, but NEVER 100% (residual OXPHOS always measured).

**On "nothing proven" (symmetric QC, held OPEN per task instruction)**: the spread by tumor
type/microenvironment is HUGE (1-71% glycolytic-ATP-fraction depending on cell line and O2 tension);
much of the literature's
biggest swings are hypoxia-driven (a classical Pasteur effect), not O2-independent; naive
glutaminase-inhibitor therapy translating the mechanism into a drug is NULL-to-marginal in real
trials; and in 2 independent cancer types the clinically most dangerous subpopulation is the LEAST
glycolytic one, inverting the naive reading. All reported below, not smoothed over.

## Pre-registered falsifiers (stated before the gates were run — verbatim from the task)

1. Does the model reproduce the MEASURED tumor glycolytic flux / lactate production relative to
   OXPHOS (fraction of ATP from glycolysis, tumor vs normal)?
2. Does it correctly predict mitochondria are NOT broken (measured residual OXPHOS fraction)?
3. Decorrelated check: the ATP-per-unit-volume/crowding argument (Vazquez) or the glucose-carbon-to-
   biomass diversion argument.
4. Symmetric QC: hold OPEN that the quantification varies hugely by tumor type/microenvironment; note
   in-vitro Crabtree vs in-vivo differences.

All four are addressed below with machine-computed gates, not narration.

---

## 0. Verified citations (26 sources; every PMID checked LIVE this session via NCBI E-utilities —
none trusted from memory, per this repo's own measured ~62-67% citation-drift-from-recall rate)

| # | Citation | PMID (verify method) | Role |
|---|---|---|---|
| 1 | Warburg O (1956). On the origin of cancer cells. *Science* 123(3191):309-14. | **13298683** (esearch+esummary; disambiguated live from a similarly-titled 1956 *Oncologia* reprint, PMID 13335077, which is a DIFFERENT record) | Original hypothesis (its "damaged respiration" form is measurably refuted below) |
| 2 | Warburg O, Wind F, Negelein E (1927). The metabolism of tumors in the body. *J Gen Physiol* 8(6):519-30. | **19872213**, PMCID PMC2140820 (esearch+esummary) | Original in-vivo quantification. **Primary numbers NOT independently extracted this session** — PDF is a scanned image, no OCR body text via EuropePMC/PMC API (checked, confirmed empty, disclosed gap); relying on citation #3's live-verified modern requotation instead |
| 3 | Koppenol WH, Bounds PL, Dang CV (2011). Otto Warburg's contributions to current concepts of cancer metabolism. *Nat Rev Cancer* 11(5):325-37. | **21508971** (esearch+esummary+efetch abstract) | Decisive refutation of "broken mitochondria"; the ~10-fold glucose-to-lactate flux requotation |
| 4 | Pfeiffer T, Schuster S, Bonhoeffer S (2001). Cooperation and competition in the evolution of ATP-producing pathways. *Science* 292(5516):504-7. | **11283355** (esearch+esummary+efetch abstract) | Foundational rate-vs-yield tradeoff (evolutionary game theory, not cancer-specific) |
| 5 | Vazquez A, Liu J, Zhou Y, Oltvai ZN (2010). Catabolic efficiency of aerobic glycolysis: the Warburg effect revisited. *BMC Syst Biol* 4:58. | **20459610**, PMCID PMC2880972 | The threshold-switch finding this doc's geometric model reproduces |
| 6 | Vazquez A, Oltvai ZN (2011). Molecular crowding defines a common origin for the Warburg effect in proliferating cells and the lactate threshold in muscle physiology. *PLoS ONE* 6(4):e19538. | **21559344**, PMCID PMC3084886 | The ATP-per-volume-density crowding mechanism (**title independently corrected this session** — commonly misremembered as "...and the growth advantage of cancer cells"; live-verified actual title ties to the muscle lactate threshold, a cross-domain corroboration) |
| 7 | Shlomi T, Benyamini T, Gottlieb E, Sharan R, Ruppin E (2011). Genome-scale metabolic modeling elucidates the role of proliferative adaptation in causing the Warburg effect. *PLoS Comput Biol* 7(3):e1002018. | **21423717**, PMCID PMC3053319 | Independent (genome-scale FBA) convergence: biomass-production-rate driver + glutaminolysis preference + "three phase" staged behavior |
| 8 | Vander Heiden MG, Cantley LC, Thompson CB (2009). Understanding the Warburg effect: the metabolic requirements of cell proliferation. *Science* 324(5930):1029-33. | **19460998**, PMCID PMC2849637 | The biosynthetic-precursor-diversion argument |
| 9 | Molenaar D, van Berlo R, de Ridder D, Teusink B (2009). Shifts in growth strategies reflect tradeoffs in cellular economics. *Mol Syst Biol* 5:323. | **19888218**, PMCID PMC2795476 | Cross-SPECIES decorrelated instance: unicellular organisms/yeast/bacteria overflow metabolism obey the same resource-allocation logic |
| 10 | Diaz-Ruiz R, Rigoulet M, Devin A (2011). The Warburg and Crabtree effects. *Biochim Biophys Acta* 1807(6):568-76. | **20804724** | In-vitro Crabtree vs in-vivo Warburg distinction (symmetric-QC hold-open item) |
| 11 | Rich PR (2003). The molecular machinery of Keilin's respiratory chain. *Biochem Soc Trans* 31(6):1095-105. | **14641005** | Background: revised (down from 36-38) modern P/O-ratio ATP yield; no single number in its own abstract, context-only |
| 12 | Zu XL, Guppy M (2004). Cancer metabolism: facts, fantasy, and fiction. *Biochem Biophys Res Commun* 313(3):459-65. | **14697210** | THE forced adversary: hypoxia (Pasteur effect), not cell-intrinsic reprogramming, as an alternative explanation — taken in its strongest form |
| 13 | Moreno-Sánchez R, Rodríguez-Enríquez S, Marín-Hernández A, Saavedra E (2007). Energy metabolism in tumor cells. *FEBS J* 274(6):1393-418. | **17302740** | Confirms cell-line-dependent spread; some lines oxidative-dominant |
| 14 | Fantin VR, St-Pierre J, Leder P (2006). Attenuation of LDH-A expression uncovers a link between glycolysis, mitochondrial physiology, and tumor maintenance. *Cancer Cell* 9(6):425-34. | **16766262** | DECISIVE causal (knockdown + rescue) proof mitochondria remain functionally capable |
| 15 | Gatenby RA, Gillies RJ (2004). Why do cancers have high aerobic glycolysis? *Nat Rev Cancer* 4(11):891-9. | **15516961** | Alternative/complementary hypothesis (acid-mediated invasion) — symmetric QC: multiple non-exclusive drivers plausible |
| 16 | Kubota R, Kubota K, Yamada S, Tada M, Ido T, Tamahashi N (1994). Microautoradiographic study... FDG uptake. *J Nucl Med* 35(1):104-12. | **8271030** | Cellular-level FDG mechanistic validation; discloses macrophages can out-uptake tumor cells (FDG-PET not perfectly tumor-specific) |
| 17 | Fletcher JW et al. (2008). Recommendations on the use of ¹⁸F-FDG PET in oncology. *J Nucl Med* 49(3):480-508. | **18287273** | Clinical consensus guideline; no specific SUV ratio number in its own abstract — disclosed gap |
| 18 | Zheng J (2012). Energy metabolism of cancer: Glycolysis versus oxidative phosphorylation (Review). *Oncol Lett* 4(6):1151-7. | **23226794**, PMCID PMC3506713 | **KEY quantitative synthesis** — see method note below |
| 19 | Mookerjee SA, Gerencser AA, Nicholls DG, Brand MD (2017). Quantifying intracellular rates of glycolytic and oxidative ATP production... *J Biol Chem* 292(17):7189-207. + erratum *J Biol Chem* 293(32):12649-52. | **28270511** + **30097494** (erratum, PMCID PMC6093231, full text independently fetched) | Modern precise ATP-yield figure (33.45/glucose) + the "glycolytic index" concept |
| 20 | Berghmans T et al. (2008). Primary tumor SUVmax... prognostic value for survival in NSCLC: meta-analysis. *J Thorac Oncol* 3(1):6-12. | **18166834** (independently re-verifying the pre-existing seed-node datapoint) | THE clinical-outcome anchor: HR=2.27 (95% CI 1.70-3.02), n=1,474, 13 studies |
| 21 | Hensley CT, Faubert B, Yuan Q, et al. (2016). Metabolic Heterogeneity in Human Lung Tumors. *Cell* 164(4):681-94. | **26853473**, PMCID PMC4752889 (independently re-verified) | DECISIVE human in-vivo refutation of "broken mitochondria": 9/9 real patients |
| 22 | Tannir NM et al. (2022). CANTATA Randomized Clinical Trial (telaglenastat+cabozantinib, RCC). *JAMA Oncol* 8(10):1411-8. | **36048457**, PMCID PMC9437824 (independently re-verified) | Symmetric QC: unselected glutaminase-inhibitor therapy is NULL |
| 23 | Lee CH et al. (2022). ENTRATA Trial (telaglenastat+everolimus, RCC). *Clin Cancer Res* 28(15):3248-55. | **35576438**, PMCID PMC10202043 (independently re-verified) | Symmetric QC: marginal, not significant at conventional two-sided α |
| 24 | Viale A et al. (2014). Oncogene ablation-resistant pancreatic cancer cells depend on mitochondrial function. *Nature* 514(7524):628-32. | **25119024**, PMCID PMC4376130 (erratum *Nature* 2026 Apr;652(8111):E9 exists, not fetched, disclosed) | Polarity-inversion finding #1: relapse-driving cells are HIGH-OXPHOS |
| 25 | Farge T et al. (2017). Chemotherapy-Resistant AML Cells... Require Oxidative Metabolism. *Cancer Discov* 7(7):716-35. | **28416471**, PMCID PMC5501738 | Polarity-inversion finding #2 (independent cancer type/group) |
| 26 | DeBerardinis RJ, Chandel NS (2016). Fundamentals of cancer metabolism. *Sci Adv* 2(5):e1600200. | **27386546**, PMCID PMC4928883 | Broad modern synthesis, topical anchor |

**Method note on #18 (Zheng 2012)** — this is the single most quantitatively load-bearing citation
(the measured glycolytic-ATP-fraction numbers) and was verified BEYOND the usual esearch/esummary/
efetch-abstract pattern: the WebFetch tool's own summarizing model reported specific percentages
(HCT116 40%→66%, HeLa 79%→29%, MCF 91%→36%, "1-64%" range), and — because a summarizing model's
report is itself an unverified claim, not machine truth, per this repo's own standing caution about
a prior "confirmed live" turning out fabricated — the RAW PMC HTML was independently fetched via
`curl` and grepped BY ME directly for every one of those exact strings. All matched byte-for-byte.
Tracing Zheng's own reference list (also independently fetched and read): the "1-64%" figure is
Zheng's citation of Zu & Guppy 2004 (#12 above); the HeLa/MCF hypoxia figure cites Rodríguez-Enríquez
et al. 2010 *Int J Biochem Cell Biol* 42:1744-51 (author/journal/year read directly from Zheng's own
live-fetched bibliography, not independently re-verified via a fresh PMID lookup this session —
disclosed scope boundary); the HCT116 p53 figure cites a reference not relocated in the list this
session (disclosed).

---

## 1. The ATP-yield-vs-flux tradeoff (grounding fact-check, machine cross-checked)

| | glycolysis | OXPHOS (complete oxidation) |
|---|---:|---:|
| ATP yield per glucose | **2** (net, substrate-level phosphorylation) | **33.45** (modern max, P/O=2.79; historical range 30-38) |
| O2 required | no | yes |
| Rate-yield ratio (OXPHOS/glycolysis) | — | **16.7x** (modern) / **15.0-19.0x** (historical spread) |

Source: Mookerjee et al. 2017 (#19), whose 2018 erratum was independently fetched and confirmed to
correct only measured C2C12 flux figures (basal total ATP 41.7→48.2, oxidative ATP 30.8→36.0
pmol/min/µg protein — up to 17% shifts) and NOT the theoretical 33.45 ATP/glucose maximum used here.
**Gate**: rate-yield ratio exceeds 10x under EVERY historical-or-modern estimate — `gate_ratio_exceeds_10x = True`.

---

## 2. The geometric model — WHY, derived from a 2-constraint resource-allocation LP

**Setup** (a real linear program, not rote algebra): two ATP-generating pathways compete for TWO
shared, limited resources — glucose-uptake capacity `G` and solvent/molecular-crowding capacity `V`
(enzyme + mitochondrial volume the cytoplasm can hold; Vazquez et al. 2010/2011, #5/#6). Each pathway
occupies volume `κ` per unit glucose-equivalent flux; mitochondria are bulkier per unit flux than
glycolytic enzymes (`κ_OXPHOS > κ_glycolysis`).

```
maximize   Y_g·J_g + Y_o·J_o                        (total ATP production rate)
subject to J_g + J_o ≤ G                              (glucose-uptake-capacity constraint)
           κ_g·J_g + κ_o·J_o ≤ V                       (solvent/crowding-capacity constraint)
           J_g, J_o ≥ 0
```

By the fundamental theorem of linear programming, the optimum sits at a vertex of the feasible
polytope — a direct geometric fact, not an assumption. Full vertex enumeration gives the general
solution for ANY (κ_g, κ_o), not one cherry-picked regime.

### Two real bugs caught and fixed via forced OODA (not a premature honest-negative)

The first-draft implementation FAILED its own pre-registered cross-check gate on the first run.
Rather than reporting a one-shot "honest negative," both failures were diagnosed to their exact root
cause and fixed:

1. **Brute-force grid bug**: the independent grid-search cross-check swept `J_g` over `[0, G]`
   without checking that `J_g` alone (with `J_o=0`) already satisfies `κ_g·J_g ≤ V` — for
   `G > V/κ_g` this let the grid evaluate physically INFEASIBLE points and report spuriously high
   "ATP" (e.g. 3.0 vs the true 2.0 at one test point — a 50% error). Diagnosed by comparing against
   `scipy.optimize.linprog` (an independent, trusted, off-the-shelf LP solver), which agreed with the
   hand-derived formula to machine epsilon while the old brute force did not. **Fix**: clip `J_g`'s
   own upper bound to `min(G, V/κ_g)` before ever considering `J_o`.
2. **Closed-form regime bug** (the deeper, more scientifically interesting catch): after fixing (1),
   `scipy_lp` and the closed form STILL disagreed at `κ_o/κ_g = 3` (closed-form/scipy_lp both gave
   2.0 vs true optimum 11.15 — a 458% error). Diagnosis: the original hand-derived 3-phase formula
   silently assumed the "Vazquez regime" (`Y_g/κ_g > Y_o/κ_o`, i.e. glycolysis wins per-volume)
   without checking it. Below `κ_o/κ_g ≈ 16.7-16.9` (see below), OXPHOS actually dominates glycolysis
   on BOTH the per-glucose AND per-volume axes — glycolysis should NEVER activate, and the old formula
   was simply wrong outside its own unstated precondition. **Fix**: replaced the regime-restricted
   formula with the fully general vertex-enumeration solver above, correct for any parameter values.

This is exactly the OODA loop this project's own discipline demands (Observe the failure → Orient by
finding the specific root cause → Decide the fix → Act and re-test) — not a shrug, and not achieved by
loosening any tolerance.

### The model's own derived, falsifiable prediction

Fixing bug 2 surfaced a genuine, sharper finding: **a mixed/Warburg-like strategy is geometrically
optimal if and only if κ_OXPHOS/κ_glycolysis exceeds Y_OXPHOS/Y_glycolysis** — mitochondria must be
more than **~16.7x** (modern estimate) to **~15.0-19.0x** (historical spread) bulkier per unit
glucose-equivalent flux than glycolytic enzymes for the crowding argument to overturn OXPHOS's
per-glucose yield advantage. Below that crossover, OXPHOS wins on both axes and glycolysis is
predicted to never activate — a real, testable, bidirectional boundary of the model, not a hand-wave.

### Worked illustrative example (`κ_o=20`, `κ_g=1`, `V=1` — an illustrative sweep, NOT Vazquez's own
internal fitted parameter, which was not available from the abstracts alone; disclosed, not fabricated)

| G (glucose-uptake capacity) | phase | J_glycolysis | J_OXPHOS | glycolytic ATP fraction |
|---:|---|---:|---:|---:|
| 0.01 | pure OXPHOS | 0.000 | 0.010 | 0.0% |
| 0.05 | pure OXPHOS (boundary, `G*_low=V/κ_o=0.05`) | 0.000 | 0.050 | 0.0% |
| 0.20 | gradual mixed | 0.158 | 0.042 | 18.3% |
| 0.50 | gradual mixed | 0.474 | 0.026 | 51.8% |
| 0.75 | gradual mixed | 0.737 | 0.013 | 77.0% |
| 1.00 | pure glycolysis (boundary, `G*_high=V/κ_g=1.0`) | 1.000 | 0.000 | 100.0% |

A smooth, monotonic, demand-driven transition from pure-OXPHOS to pure-glycolysis — exactly the
"gradual activation above a threshold" Vazquez et al. 2010 (#5) describe qualitatively — reproduced
here from first-principles LP geometry, not curve-fit to their data.

### Robustness sweep — pre-registered gate, PASS

Swept `κ_o/κ_g` log-uniformly over **[3, 300]** (25 points), spanning the derived crossover ≈16.7-16.9
on both sides:

| gate | result |
|---|---|
| Universal math-correctness (closed-form = brute-force = scipy-LP, tol 1e-6, at EVERY tested ratio) | **25/25 PASS** (max rel. error 8.2e-7 vs brute-force, 6.3e-10 vs scipy — both after fixing the two bugs above) |
| Regime-prediction correctness (mixed phase present ⟺ κ_ratio > crossover, BOTH directions) | **25/25 PASS** — 9/9 ratios below crossover correctly show NO mixed phase; 16/16 above correctly show one |
| Determinism | 2 independent process re-runs, byte-identical JSON via `diff` |

**`all_gates_pass: True`** in `warburg_metabolism_results.json` — reached honestly, through 2 real
bug-fixes, not by relaxing a threshold.

---

## 3. Falsifier 1 — elevated but NOT exclusive glycolytic ATP fraction (measured, external anchor)

| | glycolytic ATP fraction |
|---|---:|
| Normal cells (baseline, Zheng 2012 citing #18's own synthesis) | ~30% (i.e. "70% OXPHOS") |
| HCT116 colon cancer, p53 wild-type (+/+) | 40% |
| HCT116 colon cancer, p53 null (-/-) — **same O2 condition** | **66%** |
| Cancer cell lines overall (Zu & Guppy 2004 range, via Zheng 2012) | 1-64% ("does not generally exceed 50-60%") |
| HeLa, normoxia | 21% (79% OXPHOS) |
| HeLa, hypoxia | 71% (29% OXPHOS) |
| MCF breast carcinoma, normoxia | 9% (91% OXPHOS) |
| MCF breast carcinoma, hypoxia | 64% (36% OXPHOS) |

**Pre-registered gates**: (a) *elevated* — max same-O2-condition cancer figure exceeds normal baseline
by ≥10 percentage points → **PASS** (66% vs 30%, +36pp). (b) *not exclusive* — max reported cancer
figure stays below 100% → **PASS** (max 66% genotype-driven / 64% Zu-Guppy range ceiling; residual
OXPHOS always present).

**The forced adversary, in its strongest form (not dismissed)**: Zu & Guppy 2004's own conclusion is
"no evidence that cancer cells are inherently glycolytic ... some tumours might indeed be glycolytic
in vivo as a result of their hypoxic environment" — i.e. much of the "aerobic" (O2-independent)
glycolysis literature may just be a classical Pasteur effect (real hypoxia), not a special cancer
reprogramming. The HeLa/MCF numbers above are the adversary's own best ammunition: BOTH lines are
LESS glycolytic under normoxia (9-21%) than the stated normal-cell baseline (30%) and only exceed it
sharply once actual hypoxia is imposed (64-71%) — the adversary genuinely explains much of the
literature's biggest swings. **But it does not fully explain the effect away**: HCT116's 40%→66%
shift is driven by p53 genotype under the SAME O2 condition (per Zheng 2012's own text) — a genuine
O2-independent ("aerobic") component survives the adversary, and Hensley 2016's human in-vivo finding
(§4) independently corroborates that enhanced glycolysis is "common" in real tumors even where
oxidative metabolism is simultaneously present. Reported as a genuine partial-adversary-survival, not
smoothed into either "fully explained by hypoxia" or "hypoxia is irrelevant."

---

## 4. Falsifier 2 — mitochondria are NOT broken (causal + human in-vivo)

**Causal, decisive**: Fantin, St-Pierre, Leder 2006 (#14) knocked down LDH-A (blocking pyruvate→
lactate) and measured **stimulation of mitochondrial respiration** — the reverse of what "permanently
damaged" mitochondria would show — AND showed the phenotype **reverses on re-expressing LDH-A** (a
rescue design that rules out an off-target-toxicity artifact: a genuine off-target effect would not
be cleanly reversed by re-adding the single targeted gene product).

**Human in-vivo, decisive**: Hensley et al. 2016 (#21) infused ¹³C-glucose intraoperatively into 9
real NSCLC patients. "While enhanced glycolysis and glucose oxidation were common among these tumors,
we observed evidence for oxidation of multiple nutrients in **each of them**, including lactate as a
potential carbon source" — **9/9 (100%)** real human tumors showed oxidative use of non-glucose
carbon alongside enhanced glycolysis, in vivo, not a cell-culture artifact.

**Gate**: `LDH-A knockdown increases measured respiration = True` AND `fraction of in-vivo patients
with oxidative multi-substrate use ≥ 0.50 threshold` (measured 1.0/9) → **PASS**.

**Historical correction, directly on point**: Koppenol, Bounds, Dang 2011 (#3), re-examining Warburg's
own data: "this increase in aerobic glycolysis in cancer cells is often erroneously thought to occur
INSTEAD OF mitochondrial respiration and has been misinterpreted as evidence for damage to
respiration... In fact, many cancers exhibit the Warburg effect **while retaining mitochondrial
respiration**." Also requotes Warburg's own original finding in modern terms: tumor tissues metabolize
"approximately tenfold more glucose to lactate in a given time than normal tissues" — the real,
measured flux elevation Warburg himself observed, WITHOUT requiring "damaged" mitochondria to explain
it (this doc's geometric model in §2 explains it instead via a resource-allocation optimum).

**Disclosed limits**: Fantin 2006 is a single paper/model system (n=1 cell-line study, not a broad
survey). Hensley 2016's n=9 is a small, single-center, early-stage-resectable-NSCLC-only cohort —
generalization to advanced/metastatic or non-lung tumors is unmeasured here.

---

## 5. Falsifier 3 (decorrelated clinical anchor) — the metabolic phenotype has real prognostic weight

Berghmans et al. 2008 (#20; independently re-verified this session, not inherited from the prior
seed-node claim on trust) pooled 13 independent NSCLC studies, n=1,474: **combined hazard ratio 2.27
(95% CI 1.70-3.02)** for high-vs-low FDG-PET SUV → death; excluding studies using a post-hoc "best"
cutoff, HR=2.08 (95% CI 1.43-3.04). **Gate**: 95% CI excludes the null (HR=1) → **PASS**. High SUV
(the clinical signature of elevated glucose uptake this doc's model explains mechanistically) is not
just descriptively present, it carries genuine prognostic weight.

A universal clinical tumor:normal SUV RATIO number was **not** pinned to one figure this session
(searched; Fletcher et al. 2008's consensus guideline, #17, gives no such single number in its own
abstract) — disclosed gap. The mechanistic uptake chain (HIF1→GLUT1/GLUT3 overexpression→elevated
FDG trapping) is separately supported by Gatenby & Gillies 2004 (#15: "upregulation of glycolysis...
can be observed with clinical tumour imaging") and Kubota et al. 1994's cellular-level microauto-
radiography (#16), which also discloses FDG-PET is not perfectly tumor-cell-specific (macrophages/
granulation tissue can out-uptake tumor cells in some conditions).

---

## 6. WHY — three decorrelated methodologies converging on "not just ATP rate"

| method | group | finding |
|---|---|---|
| Evolutionary game theory | Pfeiffer, Schuster, Bonhoeffer 2001 (#4) | Rate-yield tradeoffs create a genuine evolutionary dilemma; high-rate/low-yield ATP production is a form of "defection" against cooperative high-yield use of shared resources |
| Genome-scale FBA + solvent-capacity constraint | Vazquez/Oltvai 2010/2011 (#5/#6), Shlomi et al. 2011 (#7) | A resource-allocation optimum favors glycolysis above a threshold flux/growth rate BECAUSE it yields more ATP per unit crowding volume, not per unit glucose; Shlomi's independent, larger genome-scale model additionally finds the shift is "a direct consequence of ... adaptation to increase BIOMASS production rate" |
| Cell-biology / signaling | Vander Heiden, Cantley, Thompson 2009 (#8) | Cancer (and all proliferating) cells divert glycolytic carbon into the biosynthesis of nucleotides, amino acids, and lipids needed to BUILD a new cell — glycolysis is not solely, or even primarily, about ATP rate |

These three genuinely independent routes (game theory / constrained optimization / direct cell
biology) converge on the same qualitative answer: aerobic glycolysis is not a broken-machine artifact,
it is an optimal allocation of a scarce resource (crowding volume, or carbon-for-biomass) under high
proliferative demand — and Vazquez 2011's larger model + Shlomi 2011 both independently flag
**glutaminolysis** activation as accompanying the lactate switch, which is the DIRECT mechanistic
bridge to §7's clinical trials below (glutaminase inhibition is the drugable consequence of this
exact prediction).

Molenaar et al. 2009 (#9) supplies the diverse-instance-space check this project's own discipline
demands: the SAME resource-allocation logic (a tradeoff between enzyme-synthesis investment and
metabolic yield) independently explains growth-rate-dependent "overflow metabolism" in unicellular
organisms/yeast/bacteria — a fully decorrelated biological system, different experimental method,
same geometric conclusion.

---

## 7. Symmetric QC — nothing proven, held OPEN as instructed

1. **Naive drug translation is NOT confirmed.** The mechanistic prediction (glutaminolysis
   dependency, §6) was tested directly in two real trials: CANTATA (#22, n=444, telaglenastat+
   cabozantinib, unselected RCC) HR=0.94 (95% CI 0.74-1.21, P=.65) — **null**. ENTRATA (#23, n=69,
   telaglenastat+everolimus) HR=0.64 (95% CI 0.34-1.20, one-sided P=.079) — **marginal, not
   significant at conventional two-sided α**. The SAME drug with a DIFFERENT combination partner
   flips from null to marginal — no single universal glutaminase-inhibitor effect size exists.
   Honest-negative = PASS per this repo's own discipline; reported, not hidden.
2. **Metabolic-dependency polarity can INVERT.** Two independent cancer types/groups — pancreatic
   cancer oncogene-ablation survivors (Viale et al. 2014, #24) and chemo-resistant AML cells (Farge
   et al. 2017, #25) — both show the clinically most dangerous (relapse-driving / chemo-resistant)
   subpopulation is HIGH-OXPHOS, LOW-glycolysis: the opposite of the bulk-tumor Warburg phenotype.
   This inverts the naive "high FDG uptake = vulnerable to glycolysis blockade" assumption for
   exactly the cells that matter most clinically.
3. **Huge spread by tumor type/microenvironment** (§3): 1-64% (Zu & Guppy range) up to 66-71%
   (genotype- or hypoxia-driven extremes) vs a ~30% normal baseline — reported as a range, never
   forced to one number.
4. **In-vitro Crabtree effect vs in-vivo Warburg effect genuinely differ** (Diaz-Ruiz, Rigoulet,
   Devin 2011, #10): much of the cell-culture literature (standard media glucose ~25 mM, vs
   physiological blood glucose ~5 mM) may reflect the reversible, short-term Crabtree effect rather
   than the sustained in-vivo Warburg phenotype — a real, disclosed, unresolved confound affecting an
   unknown fraction of the cell-line numbers cited in §3 (including Mookerjee's C2C12 measurements,
   which are additionally a non-tumorigenic myoblast line, used here only as a partially-oxidative
   comparator, not a cancer exemplar).
5. **The geometric model's exact crowding-cost parameter (κ_OXPHOS/κ_glycolysis) was not fitted to
   real biological data this session** — Vazquez et al.'s own internal fitted value was not available
   from the abstracts alone. The model is validated QUALITATIVELY (its predicted shape of the
   threshold switch, and its own derived crossover condition) against their independently-published
   conclusion, not quantitatively reproduced.
6. Warburg's own 1927 primary numbers (#2) were not independently re-extracted (scanned PDF, no OCR
   text available this session) — relying on Koppenol 2011's live-verified modern requotation.

---

## 8. Confidence tier (per `docs/MECHANISM_TRUST_LEDGER.md`'s vocabulary — assessed per-claim, not inflated to one blanket tier)

| claim | tier |
|---|---|
| Elevated-but-not-exclusive glycolytic ATP fraction (§3) | **cadaveric-or-published-plausibility** (real, external, decorrelated cell-line literature; not this twin's own measurement) |
| Mitochondria not broken — causal (Fantin 2006) | **cadaveric-or-published-plausibility** (single external causal study) |
| Mitochondria not broken — human in-vivo (Hensley 2016) | **in-vivo-anchored** (real living human tumors, intraoperative measurement) |
| FDG-PET SUV prognostic value (Berghmans 2008) | **in-vivo-anchored** (real patient survival outcomes, n=1,474) |
| Naive drug-translation null result (CANTATA/ENTRATA) | **in-vivo-anchored** (real randomized clinical trials) — reported as a disclosed negative |
| The geometric crowding model itself | **method-only-no-external-anchor** for its exact parameters; qualitatively cross-checked against 3 independently-published models (Vazquez 2010/2011, Shlomi 2011) |

---

## 9. Couples to (existing graph, not edited this session)

`ONCO-WARBURG-METABOLISM` (the SEED-DESIGN node this doc informs), `ONCO-CANCER-HALLMARKS`,
`CANCER-ROGUE-GROWTH-CELL`, `ORG-GLYCEMIC-METABOLIC` + `ORG-PANCREAS-GLUCOSE-INSULIN` (glucose-insulin
uptake side of the glucose/ATP economy this doc's flux numbers feed), `GOAL-CANCER-IMMUNE-WARGAME-TWIN`,
`ORG-METABOLIC-SYNDROME-CLUSTER`, `MOL-MITOPHAGY-QC` (per the seed node's own `couples_to`).

## 10. Files

- `scripts/msk/warburg_metabolism.py` — new, self-contained (numpy + scipy only, no OpenSim; scipy
  used purely as a 3rd independent LP-solver cross-check). Implements the general vertex-enumeration
  geometric model, the 3-way cross-check, the robustness/regime sweep, all literature-anchored
  falsifier gates, and writes the evidence JSON. Re-run with
  `source .venv-msk/bin/activate && python3 scripts/msk/warburg_metabolism.py`.
- `data/msk_smoketest/warburg_metabolism/warburg_metabolism_results.json` — every number in this doc,
  machine-written: the stoichiometry fact-check, the full geometric-model worked example (all 41
  points, all 3 phases), the 25-point robustness/regime sweep (both the universal math cross-check
  and the bidirectional regime-prediction check), all 26 citations with verification-method flags and
  quoted text, and all 5 falsifier gates with pre-registered thresholds and computed booleans.
  Verified deterministic (2 independent process re-runs, byte-identical JSON via `diff`).
- No cellular metabolomics/Seahorse panel exists for this twin's own subject data (same disclosed
  scope as `thyroid_metabolic_axis.py`, `circadian_rhythm.py`) — every number here is population/
  literature-anchored, not re-derived from subject2/walking1.

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task). Only
files created this session were touched.
