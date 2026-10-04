# MECHANISM INVASION-METASTASIS CASCADE + METASTATIC INEFFICIENCY + EMT/MET — a chain-graph model of WHY metastasis is astronomically inefficient and WHICH step is rate-limiting (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** Builds a **NEW, complementary** cell — checked live
this session against `data/MECHANISM_ANCHOR_GRAPH.json` before building (dedup discipline): a
substantial `ONCO-EMT-METASTASIS` node **already exists** (status `OPEN`, `mechanism_grade:
SEED-DESIGN`), built on Fischer 2015/Zheng 2015 lineage-tracing + genetic-KO data (per-animal
survival/incidence computed there) and a sarcoma E-cadherin-prognosis meta-analysis. That node's own
`hidden_state` is **"transient/partial-EMT state conferring invasion+stemness+drug-tolerance"** —
i.e. EMT's necessity for metastatic *incidence* vs its role in *chemoresistance*. **This doc is a
DISTINCT hidden state**: *which stage of the cascade is rate-limiting* (measured via a chain-graph
decomposition of primary cell-tracking data no prior mechanism cell had built) and *why EMT must be
reversible* (a decisive, quantitative, in-vivo forced-constitutive-vs-transient comparison the
existing node does not cover). Not a duplicate — coupled to it, not folded into it (isolation rule:
touch only files created this session). Script: `scripts/msk/metastasis_invasion_cascade.py`.
Evidence (machine-written): `data/msk_smoketest/metastasis_invasion_cascade/metastasis_invasion_cascade_results.json`
(md5 `c1b75da91907826826ba4e35074664ff`, byte-identical across repeated runs — determinism confirmed,
no RNG anywhere in this script).

**Every citation independently verified LIVE this session** via NCBI E-utilities
(esearch/esummary/efetch) and Europe PMC REST full-text (fallback for open-access PMC records). Two
of my own recalled identifiers were **wrong or ambiguous** this session — a fresh, in-session
demonstration of the disclosed ~62% citation recall-drift rate this task was pre-warned about: (1) a
guessed PMID for a Savagner wound-healing/EMT paper (9927757) resolved to an unrelated ribosomal-DNA
paper — corrected via fresh `esearch`; (2) Butler & Gullino 1975's own shedding-rate exponent came back
OCR-garbled from `efetch` plaintext ("10-6", internally inconsistent with the same abstract's own
qualitative claim) — diagnosed and resolved via an independent live secondary source (§6). Six of the
20 citations here are **REUSED** from the existing `ONCO-EMT-METASTASIS` node (Fischer/Zheng/Nieto/
Ye-Weinberg/Kroger/Pastushenko) — every one of the six was **independently re-confirmed this session**
via a fresh `esummary` call (not blindly trusted from the prior agent's own `fetched_live:true` flag),
per this repo's own discipline that a prior agent's claim is a hypothesis, not truth.

**Confidence tier, split precisely, not blended:** the multi-step-bottleneck/chain-graph finding (F1)
and the MET-reversibility falsifier (F4, Tsai 2012 + Ocana 2012) are **in-vivo-anchored** (real mice,
real causal genetic/pharmacological manipulation, machine-recomputed statistics). The human
cross-species check (F2, Coumans 2013) is **in-vivo-anchored at the population-epidemiological level**
(38,715 real patients) but is a fitted MODEL, not a direct cell-tracking measurement. The EMT-TF
mechanism (F3) is **in-vitro/xenograft-anchored** (cell-culture + mouse xenograft causal
gain/loss-of-function). Disclosed per-claim in §8, not smoothed into one blanket tier.

---

## 0. Falsifiers (pre-registered, matching the task) + verdict up front

| # | Falsifier | Verdict |
|---|---|---|
| F1 | The invasion-metastasis cascade, decomposed as a chain graph from Luzzi 1998's own primary cell-tracking data, identifies **POST-extravasation colonization** (not intravasation/extravasation) as the dominant bottleneck | **PASS** — colonization stages carry **97.4%** of total log-loss; naive "extravasation is the bottleneck" adversary directly falls (only 2.6% share); overall computed efficiency **0.02%**, exactly at the task's own stated upper bound |
| F2 | Cross-species/method over-determination: an independent human population-epidemiological model (Coumans 2013) is ALSO extremely inefficient, and Luzzi's point estimate falls inside Coumans' own independently-quoted murine-literature range | **PASS** — human 1.7×10⁻⁶%; Luzzi's 0.02% sits inside the held-out [0.0001%, 0.588%] murine-literature bracket; **honest, disclosed ~4,200-fold human-vs-murine-median gap**, not swept away |
| F3 | EMT-transcription-factor convergence: ≥3 independent genes/labs directly repress E-cadherin | **PASS**, 3/3 (Snail/Cano2000, Twist/Yang2004, ZEB1/Eger2005) |
| F4 | **THE decisive falsifier**: forced-constitutive (irreversible) EMT SUPPRESSES metastatic colonization relative to transient/reversible EMT, causally, in vivo | **PASS** — Tsai 2012: 86% (12/14, reversible) vs 23% (3/13, sustained), Fisher exact OR=20.0, p=0.0018; Ocana 2012 (decorrelated lab/system, same journal issue) independently confirms MET-is-required |
| Decorrelated | CTC shedding (millions/day/g) vs realized colonization rate (<0.02%) | **PASS** — Butler & Gullino 1975 (exponent-disambiguated via Mathias 2020): 1.7–4.6×10⁶ cells/24h/g shed; Meng 2004 (real human patients): CTC half-life 1–2.4h, continuously replenished, mostly dying |

`overall_pass = True` (all 4 falsifier-gates + the decorrelated check). **One real bug was caught and
fixed via OODA this session, disclosed not hidden** (§1): a first-draft exact-boundary float comparison
failed on IEEE-754 representation noise (3.5×10⁻¹⁸ off, not a real discrepancy) — fixed with
`math.isclose`, not by loosening any scientific threshold.

---

## 1. F1 — the cascade as a chain graph: colonization, not extravasation, is the bottleneck (measured)

**Geometric frame (not a metaphor):** the invasion-metastasis cascade is modeled as a **directed,
acyclic chain graph** — source = primary tumor, sink = macrometastasis, each edge weighted by a
measured stage-survival probability. For a serial chain (no parallel paths), overall success
probability is the **product** of edge weights; in log-space this is an **additive decomposition**
(`log(P_total) = Σ log(P_i)`) — the same structure as optical density or impedance in series. The
rate-limiting stage is the edge contributing the **largest magnitude** to that sum — a directly
measured graph bottleneck, not an assumed one.

**Primary data — Luzzi et al 1998 (PMID 9736035),** B16F1 melanoma injected intraportally into mouse
liver, in vivo videomicroscopy + a novel cell-accounting assay + Ki-67(proliferation)/TUNEL(apoptosis)
immunohistochemistry. Quoted verbatim (efetch, this session): **(1)** 80% of injected cells survived
the liver microcirculation and extravasated by day 3; **(2)** of those, only 1-in-40 began to grow,
forming micrometastases (4–16 cells) by day 3; **(3)** of those micrometastases, only 1-in-100
progressed to macroscopic tumors by day 13 (**most micrometastases disappeared**); **(4)** 36% of
injected cells remained at day 13 as solitary, dormant cells (proliferation index 2%, apoptosis index
3% — vs macroscopic-tumor-resident cells at 91% proliferation, 6% apoptosis/necrosis).

| Stage (edge) | Survival fraction | log₁₀ loss | Share of total log-loss |
|---|---:|---:|---:|
| Arrest → survival → extravasation | 0.80 | 0.097 | **2.6%** |
| Growth initiation → micrometastasis (day 3) | 1/40 = 2.5% | 1.602 | **43.3%** |
| Micrometastasis → macroscopic tumor (day 13) | 1/100 = 1.0% | 2.000 | **54.1%** |
| **Overall (product)** | **0.02%** | **3.699 (total)** | **100%** |

**Pre-registered gate**: post-extravasation stages (rows 2+3) carry >50% of total log-loss →
**PASS, 97.4%**. **Forced adversary** (the naive, pre-Luzzi assumption that arrest/extravasation
itself is the dominant bottleneck): tested directly on the *same* primary data — its own log-loss share
is only 2.6%, far under the 50% bar → **adversary falls**. Luzzi's own conclusion, quoted verbatim,
**is** this doc's falsifier claim, independently re-derived here via an explicit chain-graph
computation rather than re-quoted: *"metastatic inefficiency is principally determined by two distinct
aspects of cell growth AFTER extravasation: failure of solitary cells to initiate growth and failure of
early micrometastases to continue growth into macroscopic tumors."*

**A second, independent leg from the SAME paper** (marker-based, not cell-counting-based) converges on
the identical conclusion: the 36% solitary-dormant population shows near-zero proliferation (2%) AND
near-zero apoptosis (3%) — direct IHC evidence that **dormancy, not death**, differentiates the
non-colonizing majority. This is a growth-initiation/microenvironmental-licensing problem, not a
cell-survival problem, at the post-extravasation stage.

**Overall computed efficiency: 0.02%** — landing exactly at the task's own pre-registered "~0.01-0.02%"
upper bound, computed here from Luzzi's raw stage fractions (0.80 × 1/40 × 1/100), not fit to that
target. **Chambers, Groom, MacDonald (2002), Nat Rev Cancer, PMID 12154349** — the same senior authors
(Chambers & Groom) as Luzzi 1998 itself — is the canonical review framing this same conclusion at the
field level; its own abstract (live-fetched this session) is a short summary with no single quotable
number, so it is cited here for identity/framing only, disclosed, not as the source of the 0.02% figure
(which is this doc's own computation from Luzzi's primary data).

**OODA bug-catch, disclosed not hidden**: the first-draft exact `<=0.02` boundary check FAILED on the
first run (Python float arithmetic gives `0.020000000000000004`, not exactly `0.02`). Diagnosed as
IEEE-754 representation noise (diff ≈3.5×10⁻¹⁸, sixteen orders of magnitude below any measurement
precision claimed here) — fixed via `math.isclose(rel_tol=1e-9)`, not by loosening the scientific
threshold. Reported per this repo's own discipline: a one-shot fail is not an honest negative until
Orient is applied.

---

## 2. F2 — cross-species/method over-determination (Coumans 2013), with the gap honestly disclosed

**Coumans, Siesling, Terstappen (2013), BMC Cancer, PMID 23763955** built an independent 5-stage
mathematical model of the metastatic cascade — their own words, quoted verbatim (full text,
live-fetched this session): *"1. local tumor growth; 2. dissemination into circulation; 3. survival in
circulation; 4. extravasation into tissue; and 5. growth into a metastasis"* — fit to **38,715**
T-any-N-any-M0 breast cancer patients (Netherlands Cancer Registry), calibrated on 1,489 T1b patients.
Result, quoted verbatim: **"a metastatic efficiency of 1.7×10⁻⁸ metastases formed per disseminated cell
(range 1.3×10⁻⁸–4.2×10⁻⁸), or approximately 60 million disseminated cells per formed
macrometastasis"** — and, directly comparing to the murine literature: **"substantially less efficient
than the murine model median estimate of 1 metastasis in 14,000 disseminated cells (range 1 in 170 to
1 in 1,000,000)."**

| Source | Method | Species | Metastatic efficiency |
|---|---|---|---:|
| Luzzi 1998 (this doc's F1) | single-cell in vivo videomicroscopy | mouse | 0.02% (1 in 5,000) |
| Coumans 2013's own quoted **murine-literature range** | secondary synthesis | mouse | 0.0001%–0.588% (1 in 1,000,000 to 1 in 170) |
| Coumans 2013's own **human** model | population epidemiological fit | **human** | 0.0000017% (1 in ~60,000,000) |

**Gate**: Luzzi's point estimate falls inside Coumans' own independently-sourced (held-out, not fit to
it) murine-literature range → **PASS**. Both the mouse point-estimate and the human model estimate are
**far below** a pre-registered 1%-inefficiency threshold → **PASS**.

**Symmetric QC, disclosed not swept**: the human estimate is **~4,202-fold** less efficient than
Coumans' own quoted murine-literature *median*, and the mouse-vs-human gap is a genuine, unresolved,
disclosed spread — the two methods (single-cell tracking in one organ/cell-line vs population-scale
epidemiological back-calculation across a real registry) agree only on **order-of-magnitude
extremity** (both ≪1%), which is the actual certified claim here; the exact digit is
species/method-dependent and is NOT collapsed into one number.

---

## 3. F3 — EMT-transcription-factor convergence: Snail, Twist, ZEB1, three independent groups

| Gene | Citation | Sufficiency (ectopic expression → EMT) | Necessity (loss-of-function reverses it) |
|---|---|---|---|
| **Snail** | Cano et al 2000, Nat Cell Biol, PMID 10655586 | YES — *"mouse Snail is a strong repressor of transcription of the E-cadherin gene... epithelial cells that ectopically express Snail adopt a fibroblastoid phenotype and acquire tumorigenic and invasive properties"* | not tested in this paper (disclosed) |
| **Twist** | Yang et al 2004, Cell, PMID 15210113 | YES — *"ectopic expression of Twist results in loss of E-cadherin-mediated cell-cell adhesion, activation of mesenchymal markers, and induction of cell motility"* | YES — *"suppression of Twist expression in highly metastatic mammary carcinoma cells specifically inhibits their ability to metastasize"* |
| **ZEB1 (δEF1)** | Eger et al 2005, Oncogene, PMID 15674322 | YES — *"ectopic expression of deltaEF1...was sufficient to downregulate E-cadherin and to induce EMT"* | YES — *"RNA interference-mediated downregulation of deltaEF1...was sufficient to derepress E-cadherin expression and restore cell to cell adhesion"* |

All three quotes independently live-verified this session (NCBI efetch). **Gate** (pre-registered:
≥3 independent genes/groups converge on direct E-cadherin transcriptional repression as the
EMT-initiating mechanism) → **PASS, 3/3**. This is the mechanistic basis for the cascade's *invasion*
step — upstream of, and mechanistically distinct from, F1's colonization bottleneck.

---

## 4. F4 — THE decisive falsifier: forced-constitutive EMT suppresses colonization (computed, not narrated)

**Tsai, Donaher, Murphy, Chau, Yang (2012), Cancer Cell, PMID 23201165** — a doxycycline-inducible
K5-Twist1 **spontaneous squamous cell carcinoma** mouse model. Two administration routes give the
**same drug, same gene, different reversibility**: **topical** doxycycline (local/transient Twist1
induction — permits reversion, i.e. MET, once cells reach a distant site) vs **oral** doxycycline
(systemic/sustained induction — Twist1 stays ON everywhere, blocking MET). Quoted verbatim (full text,
live-fetched via NCBI PMC this session — Europe PMC's `fullTextXML` 404'd, disclosed fallback route
used): *"12 out of 14 K5-Twist1 mice (86%) receiving topical doxycycline developed distant metastases.
In contrast, only 3 out of 13 K5-Twist1 mice (23%) receiving oral doxycycline developed distant
metastases... all distant metastatic lesions showed no Twist1 expression [regardless of route]...only
reversible, but not irreversible induction of Twist1 significantly promotes distant metastasis."*
Control (no Twist1 induction) baseline: 27–33%.

**This doc's own machine-computed statistics** (Fisher's exact test on the 12/14-vs-3/13 contingency
table, `scipy.stats.fisher_exact` — not narrated from the paper's own reported p-value):

| | Topical (reversible) | Oral (sustained/constitutive) |
|---|---:|---:|
| Metastasis-positive / n | 12/14 | 3/13 |
| % | **85.7%** | **23.1%** |

**Odds ratio = 20.0, Fisher exact two-sided p = 0.00184** (pre-registered α=0.05) → **PASS**. Relative
risk (reversible vs sustained) = **3.71×**. Direction correct (reversible > sustained) → **PASS**.

**Every successful macrometastasis, regardless of route, had already reverted to Twist1-negative** —
i.e. completed MET — by the time it was detected: MET at the colonization site is essentially
**obligate** for outgrowth, not merely correlated with better odds.

**Decorrelated confirmation — Ocaña, Corcoles, Fabra, Moreno-Bueno, Acloque, Vega, Barrallo-Gimeno,
Cano, Nieto (2012), Cancer Cell, PMID 23201163** — a **different lab/system** (Nieto/Cano, Spain, vs
Tsai/Yang, UCSD), published the **same issue** of Cancer Cell (2012 Dec 11) with a joint editorial
cross-reference. Quoted verbatim (abstract, live-fetched; no PMCID found this session, so this citation
is abstract-tier only, disclosed): *"The loss of Prrx1 is required for cancer cells to metastasize in
vivo, which revert to the epithelial phenotype concomitant with the acquisition of stem cell
properties. Thus, unlike the classical EMT transcription factors, Prrx1 uncouples EMT and stemness."*

**Gate**: significant (p<0.05) + correct direction + confirmed in ≥2 independent decorrelated systems →
**PASS, 2/2** (Tsai's SCC model + Ocaña's Prrx1 system). **My own initially-recalled publication date
for Ocaña 2012 was wrong** (recalled as mid-2012; live esummary gives 2012 Dec 11) — corrected, not
asserted from memory, a fresh in-session instance of the disclosed recall-drift risk.

---

## 5. Decorrelated check — CTC shedding (millions/day) vs realized colonization (<0.02%)

**Butler & Gullino (1975), Cancer Res, PMID 1090362** measured cell shedding into efferent blood of rat
MTW9 mammary carcinomas. Live `efetch` returned an OCR-garbled exponent ("3.2×10⁻⁶ / 4.1×10⁻⁶
cells/24h/g") that is **internally inconsistent** with the same abstract's own qualitative claim ("a
2-g MTW9 carcinoma pours enough cells into the host circulation to transplant the tumor every 24 hr" —
impossible at a 10⁻⁶ rate). **Disambiguated** via an independent live secondary source — **Mathias,
Chang, Martin, Vitolo (2020), Cancers, PMID 32245166** — full text quoted verbatim: *"Butler and Gullino
determined that 1.7–4.6×10⁶ tumor cells were shed every 24 hour per gram from growing hormone dependent
rat mammary tumors."* A live-caught-and-resolved digit ambiguity, not silently picked either way.

Per this doc's own F1/F2 computations, producing **one** macrometastasis requires, on average,
~**5,000** disseminated cells (Luzzi mouse model) to ~**59,000,000** (Coumans human model) — against a
shedding rate of **millions of cells per day** from even a small (2 g) tumor. **Meng, Tripathy, Frenkel,
et al (2004), Clin Cancer Res, PMID 15623589** — **real human patients** (13/36 breast-cancer-dormancy
candidates 7–22 years post-mastectomy had detectable CTCs vs 1/26 controls, P=0.0043) — independently
confirms the shedding side is real and sustained, not a rodent-only artifact: *"CTCs...had a half-life
measured in 1 to 2.4 hours...must be replenished every few hours by replicating tumor cells somewhere in
the tissues."* CTCs are a continuously-turned-over, mostly-dying population, not a static reservoir —
millions shed, vanishingly few colonize, across every model checked.

---

## 6. Symmetric QC — the Fischer/Zheng controversy, held OPEN, not resolved

**Fischer et al (2015), Nature, PMID 26560033** and **Zheng et al (2015), Nature, PMID 26560028** —
both **independently re-confirmed live this session** (fresh esummary + efetch, not trusted from the
existing `ONCO-EMT-METASTASIS` node's prior `fetched_live:true` flag) — used EMT lineage-tracing
(breast-to-lung) and genetic Snail/Twist deletion (pancreatic GEMM) respectively to show that
genetically blocking EMT does **NOT** reduce metastatic incidence or systemic dissemination. Zheng's
own words, quoted verbatim: *"EMT suppression in the primary tumour does not alter the emergence of
invasive PDAC, systemic dissemination or metastasis... Snail- or Twist-induced EMT is not rate-limiting
for invasion and metastasis."* The existing graph node's own per-animal computation (not re-derived
here, cited): any-site metastasis incidence unchanged, 54.8% (EMT-intact control) vs 59.1%
(EMT-genetically-blocked), PMID 26560028. **BOTH papers also show EMT drives chemoresistance** — a
completely separate axis from invasion/incidence.

This is a **genuinely, actively contested** finding. **Ye & Weinberg (2015), Trends Cell Biol, PMID
26437589** (identity re-confirmed live this session) explicitly represent the opposing school: EMT/
epithelial-mesenchymal plasticity as *"a central regulator of cancer progression."* **Held OPEN here,
not resolved**, exactly as instructed — this doc does not adjudicate the controversy.

**An interpretive bridge is offered, explicitly flagged as such, NOT as a new certified claim**: this
doc's own F1 finding (colonization carries 97.4% of the total log-loss; the arrest/extravasation stage
carries only 2.6%) is at least *coherent with* Fischer/Zheng's null incidence result — if EMT's
mechanistic role is specifically in the invasion/intravasation step (F3), and that step contributes only
a small share of total log-inefficiency, there may be little headroom for genetically ablating EMT to
move the *overall* incidence number, even if EMT were completely dispensable for that early step. This
is offered for discussion only; it is **not independently statistically tested** by anything in this
doc's own gates, and is not swept into `overall_pass`.

**A separate, important scope caveat**: Luzzi 1998's own experimental design starts **after**
intravasation (B16F1 cells are injected directly into the portal circulation) — it can only speak to
the "back half" of the cascade (circulation survival → extravasation → colonization), not to the
relative efficiency of local invasion/intravasation itself (EMT's own proposed domain, F3). The claim
that intravasation is *also* not rate-limiting rests on Zheng 2015's own separate finding (systemic
dissemination unchanged by EMT-TF deletion), not on Luzzi's data — two different primary sources for
two different halves of the falsifier, not merged into one dataset.

---

## 7. Cross-domain couplings, live-verified this session

- **Cancer/Warburg + angiogenesis** (`docs/MECHANISM_WARBURG_METABOLISM.md`,
  `docs/MECHANISM_ANGIOGENESIS_VEGF.md`): metabolic reprogramming co-occurs with EMT in aggressive
  subclones; hypoxia (angiogenesis doc's own F1 HIF-1α axis) is a well-known EMT inducer upstream of
  this doc's F3 — not independently re-derived here, a prospective coupling pointer.
- **Cytoskeleton** (`docs/MECHANISM_CYTOSKELETON_DYNAMICS.md`): EMT's acquired migratory phenotype is
  mechanistically actin-treadmilling-driven; that doc's own live-verified cofilin-knockdown finding
  (max cell length 52.16→99.33µm, P<0.0001) is cited there, not re-derived here.
- **Wound-healing / MET** (`docs/MECHANISM_WOUND_HEALING.md`): a **genuine, live-verified** cross-domain
  link, not an asserted analogy. **Savagner, Kusewitt, Carver, Magnino, Choi, Gridley, Hudson (2005), J
  Cell Physiol, PMID 15389643** — quoted verbatim: *"During re-epithelialization of cutaneous wounds,
  keratinocytes recapitulate SEVERAL aspects of...EMT, including migratory activity and reduced
  intercellular adhesion... Epithelial cell outgrowth from skin explants of Slug knockout mice was
  severely compromised."* Slug is a **Snail-family member** — the same TF family modeled in F3 drives a
  **partial**-EMT-like program in wound-healing keratinocyte migration. Precisely disclosed, not
  overclaimed: Snail itself (distinct from Slug) was explicitly **not** modulated in this same assay —
  gene-specific, not a blanket EMT program. A corroborating French-language review (Savagner 2009, PMID
  20666012, identity-confirmed live, language barrier disclosed) is cited for context only.
- **`ONCO-EMT-METASTASIS`** (existing graph node): coupled to, not duplicated — see header and §6.
- **`GOAL-CANCER-IMMUNE-WARGAME-TWIN`, `ONCO-THERAPY-RESISTANCE-MECHANISM`**: carried forward from the
  existing node's own coupling list (chemoresistance axis of Fischer/Zheng 2015), not re-derived here.

---

## 8. Confidence tier (per `docs/MECHANISM_TRUST_LEDGER.md` vocabulary, assessed per-claim)

| Claim | Tier |
|---|---|
| F1 cascade chain-graph decomposition / colonization is rate-limiting (Luzzi 1998) | **in-vivo-anchored** (real mouse, in vivo videomicroscopy + IHC, direct primary cell-tracking) |
| F2 human cross-species check (Coumans 2013) | **in-vivo-anchored at population scale** (38,715 real patients) but a **fitted model**, not direct cell-tracking — disclosed distinction |
| F3 EMT-TF convergence (Snail/Twist/ZEB1) | **in-vitro/xenograft-anchored** (cell-culture + mouse xenograft causal gain/loss-of-function) |
| F4 MET-reversibility falsifier (Tsai 2012, Ocaña 2012) | **in-vivo-anchored** (real mice, causal genetic/pharmacological route manipulation, machine-recomputed statistics) |
| Decorrelated CTC-shedding-vs-realized-rate (Butler & Gullino 1975 + Meng 2004) | **in-vivo-anchored**, rat (shedding) + **real human patients** (Meng 2004 CTC dormancy) |
| Fischer/Zheng controversy (§6) | **in-vivo-anchored**, genuinely **contested** — held OPEN, not adjudicated |
| Wound-healing/Slug coupling (§7) | **in-vivo-anchored** (mouse knockout + human keratinocyte in vitro), a genuine but partial (not full-EMT) mechanistic overlap |

---

## 9. Gates — machine-printed, not narrated

```
F1_post_extravasation_dominant_bottleneck                    True   (97.4% share)
F1_naive_extravasation_adversary_falls                        True   (2.6% share, <50% bar)
F1_overall_efficiency_within_task_stated_range                True   (0.02%, float-boundary bug caught+fixed)
F1_pass                                                        True
F2_luzzi_point_estimate_within_independent_murine_lit_range    True
F2_both_species_highly_inefficient                             True
F2_pass                                                        True
F3_emt_tf_convergence_pass                                     True   (3/3: Snail, Twist, ZEB1)
F4_tsai_2012_fisher_exact_significant                          True   (OR=20.0, p=0.00184)
F4_tsai_2012_direction_correct                                 True   (85.7% vs 23.1%)
F4_decorrelated_2_of_2_systems_confirm                         True   (Tsai + Ocana)
F4_pass                                                        True

OVERALL_PASS: True
```
Deterministic: byte-identical md5 `c1b75da91907826826ba4e35074664ff` confirmed across repeated process
runs (no RNG anywhere in this script — closed-form arithmetic + one deterministic `scipy.stats.
fisher_exact` call). All floats independently re-verified finite in a fresh process (0 non-finite
values found, walked recursively over the full parsed JSON).

---

## 10. Honest gaps — symmetric QC: what this does NOT prove

1. **F1's entire chain-graph decomposition rests on ONE study, ONE cell line (B16F1 melanoma), ONE
   target organ (liver via intraportal injection)** — generalization across tumor types/organs/species
   beyond Luzzi 1998 and Coumans 2013's breast-cancer model is not independently verified this session.
2. **Luzzi's experimental design starts after intravasation** (§6) — it cannot itself speak to the
   relative efficiency of local invasion/intravasation; that claim rests on the separate Zheng 2015
   finding, not merged into Luzzi's dataset.
3. **The human-vs-mouse metastatic-efficiency gap is ~4,200-fold and unresolved** (§2) — both numbers
   agree only on order-of-magnitude extremity, not on the exact digit.
4. **F4's decisive falsifier rests on 2 mouse systems** (Tsai 2012 SCC model, Ocaña 2012 — no PMCID
   found this session for Ocaña, abstract-tier only) — independent human confirmation of the specific
   "forced-constitutive EMT suppresses colonization" causal claim was not located this session.
5. **The Fischer/Zheng-vs-Ye/Weinberg controversy is explicitly held OPEN** (§6), not adjudicated; the
   offered reconciliation is an interpretive bridge, not an independently tested claim, and is excluded
   from `overall_pass`.
6. **Butler & Gullino 1975 is a RAT model** — used for order-of-magnitude shedding-rate framing only,
   cross-checked against Meng 2004's independent real-human CTC finding, not treated as a human
   quantitative anchor itself.
7. **Thiery 2002 (PMID 12189386)** has no MEDLINE abstract text retrievable this session — identity/
   citation-only, no number sourced from it.
8. **The wound-healing coupling (§7) is a partial, gene-specific overlap** (Slug, not the full Snail/
   Twist/ZEB program) — disclosed precisely, not overclaimed as "wound healing IS cancer EMT."
9. **Reduced, literature-anchored, closed-form model** — not a spatial/agent-based simulation of the
   physical cascade; F1's chain graph has exactly 3 edges, matching Luzzi's own reported stage
   granularity, not a finer-grained mechanistic subdivision.
10. **No cellular/molecular-scale subject-specific data exists for this twin** — every number here is
    literature/population-anchored, same disclosed scope as the sibling Warburg/angiogenesis/
    cytoskeleton/wound-healing docs.

---

## 11. Files

- `source_repository/scripts/msk/metastasis_invasion_cascade.py` — the model: 20-entry
  `CITATIONS` dict (tiered quant/identity/mechanism/reused, each with PMID+DOI+exact quoted finding),
  the F1 chain-graph decomposition (Luzzi 1998, incl. the disclosed float-boundary OODA fix), F2 the
  cross-species over-determination (Coumans 2013), F3 the EMT-TF convergence table (Snail/Twist/ZEB1),
  F4 the MET-reversibility falsifier (Tsai 2012 Fisher-exact recomputation + Ocaña 2012 decorrelated
  confirmation), the CTC-shedding-vs-realized-rate decorrelated check, the symmetric-QC dict, the gates
  block, and the evidence-JSON writer. Run with `source .venv-msk/bin/activate && python3
  scripts/msk/metastasis_invasion_cascade.py` (<1s wall time, stdlib+scipy only, no OpenSim dependency,
  deterministic — confirmed byte-identical across 2 runs, no RNG anywhere in this script).
- `source_repository/data/msk_smoketest/metastasis_invasion_cascade/metastasis_invasion_cascade_results.json`
  — full machine-written evidence (all 20 citations, F1-F4 computed numbers, the decorrelated CTC
  check, symmetric QC, gates, `overall_pass`). md5 `c1b75da91907826826ba4e35074664ff`; all floats
  independently re-verified finite in a fresh process (0 non-finite values found).
- `source_repository/data/MECHANISM_ANCHOR_GRAPH.json` — the existing `ONCO-EMT-METASTASIS`
  node (`grep -n '"ONCO-EMT-METASTASIS"' data/MECHANISM_ANCHOR_GRAPH.json`), status `OPEN`,
  `mechanism_grade: SEED-DESIGN`, **unchanged this session** (fold not performed, isolation rule —
  touch only files created this session; consistent with every sibling `MECHANISM_*` doc's own
  precedent of stating results in prose/evidence-JSON only, per `docs/MECHANISM_HARDENED_CONVENTIONS.md`
  §2/§4, which requires the separate `mechanism_fold → fold_gate_v2` pipeline to promote anything here
  into a canonical graph node/edge).
- Sibling docs (read for convention/coupling, not re-litigated, not modified):
  `docs/MECHANISM_WARBURG_METABOLISM.md`, `docs/MECHANISM_ANGIOGENESIS_VEGF.md`,
  `docs/MECHANISM_CYTOSKELETON_DYNAMICS.md`, `docs/MECHANISM_WOUND_HEALING.md`,
  `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md`, `docs/MECHANISM_HARDENED_CONVENTIONS.md`,
  `docs/MECHANISM_PREAMBLE.md`.

No git commit, no git push performed (isolation respected). Only files created this session were
touched; `data/MECHANISM_ANCHOR_GRAPH.json` and all sibling docs were read-only.
