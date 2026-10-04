# MECHANISM GUT MICROBIOME — bacterial census, SCFA production, butyrate as colonocyte fuel (2026-07-22)

Builds the GUT MICROBIOME organ-system layer as a **quantitative, mechanistic census + fermentation-
output model** — explicitly NOT the unfalsifiable "dysbiosis" hand-waving that dominates the popular
microbiome literature. Three numbers, each independently live-verified and cross-checked this
session: (1) colonic bacterial density and the corrected bacteria:human-cell ratio (~1:1, refuting
the entrenched 10:1 myth), (2) short-chain-fatty-acid (SCFA) production ratio and total flux, and
(3) butyrate as the dominant colonocyte fuel. Couples to GI (fermentation substrate), immune (SCFA
→ Treg), and metabolic (caloric harvest). Script: `scripts/msk/gut_microbiome.py`. Evidence:
`reports/probes/gut_microbiome.json`.

**Headline: PRIMARY falsifier PASS (4/4 conjunctive gates), secondary corroboration PASS (6/6),
15/16 gates overall** (the 1 designed-to-fail gate is a deliberate regional-specificity probe, not
a defect — see Section B). This is a narrower, more basic quantitative layer than this repo's
existing microbiome cells (`MOL-MICROBIOME-METAGENOMICS-HOST-COUPLING`,
`SYS-GUT-MICROBIOME-HOST-AXIS`, `MOL-MICROBIOME-GUT-BRAIN-HUB`, `MICROBIOME-HOST-METABOLIC` — all
checked read-only in `data/MECHANISM_ANCHOR_GRAPH.json`, none edited) — it does not re-litigate
disease-association or gut-brain causality, it builds the census/flux arithmetic underneath them.

## 0. Scope and disambiguation from existing MICROBIOME-* graph nodes

Checked live (read-only; another instance writes this file concurrently, never edited here):
`MOL-MICROBIOME-METAGENOMICS-HOST-COUPLING` (metagenome→metabolite→IBD/CGM coupling, graded
MEASURED-B, anchor corroboration only MARGINAL — `reports/probes/mt_microbiome.json`),
`SYS-GUT-MICROBIOME-HOST-AXIS` / `MOL-MICROBIOME-GUT-BRAIN-HUB` (gut-brain/immune causal-vs-
correlational triangulation, three asymmetric-confidence edges, psychobiotics flagged
open/contested — `data/body_twin/agent_outputs/*.json`), and `MICROBIOME-HOST-METABOLIC` (function-
vs-outcome cross-cohort triangulation, already flags "fecal SCFA concentration != production flux"
re: Schwiertz 2010). **None of these builds the basic quantitative physiology layer this document
builds**: does the directly-measured bacterial census reproduce the corrected ~1:1 ratio, does SCFA
production/ratio/absorption arithmetic self-consistently close, and is butyrate really the dominant
colonocyte fuel. This document is read-only with respect to all of them — no PMID is re-verified
from those cells' own evidence; where their conclusions are relevant they are cross-referenced by
ID, not re-derived (Section E).

Population-parametrized census/production arithmetic, same class of layer as
`gi_absorption_transit.py` / `renal_filtration.py` / `hepatic_clearance.py` — no per-subject data,
no 16S/metagenomic sequencing re-analysis. Pure Python/numpy, `.venv-msk`.

## 1. Geometric structure (derive, don't assert)

- The bacteria:human-cell ratio is a **ratio of two independent volume-integrals**: bacteria =
  density × colon-content volume; human cells = dominated by hematopoietic-lineage density × blood/
  marrow volume. The historical **10:1 myth is a geometric (volumetric) error, not a density
  error**: Luckey (1972) used the *same* bacterial density this document uses (10^11/g) but
  integrated it over the **wrong volume** — 1 L (the whole alimentary tract) instead of the correct
  0.4 L (colon-content-only). Upstream-GI bacterial density is genuinely 3-8 orders of magnitude
  lower than colonic density (gastric acid + bile + fast transit suppress it) — a real anatomical
  fact, not an arbitrary refit (Section A).
- SCFA absorption is a **mass-balance closure**: production − fecal excretion = absorbed, and the
  resulting fraction must equal the independently-reported ~95% — an over-determined system (3
  separately-measured quantities constrained by 1 conservation law), not a free parameter
  (Section B).
- Colonocyte fuel preference is a **substrate-competition** question, not a single-substrate
  artifact: the forced adversary is whether butyrate's dominance survives when a competing substrate
  (glucose) is added in the same assay (Section C).

## 2. Citations — 15 PMIDs, every one esearch+efetch LIVE this session (NCBI eutils), several via
full PMC HTML fetch (not abstract-only) to pull exact in-text numbers not visible in the abstract

This repo's own measured **~62% citation-drift rate from memory** is why every number below traces
to a live fetch. **Two drifts were caught this session**: Roediger 1982's PMID was initially
misrecalled as 7084603 — the live search result is **7084619**. A prior agent's pass (data/body_twin/
agent_outputs/microbiome-gut-brain-hub__a857dea352c10e12a.json) reported Bäckhed 2004 as "+57% total
body fat, +61% epididymal fat" — this session's own live abstract fetch shows the actual headline is
a single **"60% increase in body fat content"** figure; the 57/61 sub-breakdown is not verified this
session and is not repeated here as independently confirmed.

| # | Citation | PMID | DOI | Role / verbatim number |
|---|---|---|---|---|
| 1 | Sender, Fuchs, Milo 2016, *PLoS Biol* 14(8):e1002533 | **27541692** | 10.1371/journal.pbio.1002533 | PRIMARY census anchor. "total number of bacteria in the 70 kg reference man to be 3.8×10^13 ... revise past estimates to 3.0×10^13 human cells ... updates the widely-cited 10:1 ratio, showing that the number of bacteria in the body is actually of the same order as the number of human cells." Full text (PMC4991899): B/H=1.3, ±25% uncertainty, 53% population CV. Colon density 10^11/mL (Table 1); colon volume 0.4 L (Box 1). Table 2: 12 independent empirical stool-density studies, 1966-2008, range ~0.35-3.2×10^11/g wet stool, CV 24-78%. |
| 2 | Luckey 1972, *Am J Clin Nutr* 25(12):1292-4 | **4639749** | 10.1093/ajcn/25.12.1292 | Sender 2016's own ref [3] — the traceable **origin** of the 10:1 myth: "assumes the volume of the alimentary tract to be 1 liter ... multiplies this volume by the number density of bacteria, known to be about 10^11 bacteria per gram" — same density, wrong (whole-tract, not colon-only) volume. |
| 3 | Savage 1977, *Annu Rev Microbiol* 31:107-33 | **334036** | 10.1146/annurev.mi.31.100177.000543 | Sender 2016's own ref [2] — cited as a/the source that propagated the pre-2016 ~10^14-bacteria consensus figure. |
| 4 | Cummings, Pomare, Branch, Naylor, Macfarlane 1987, *Gut* 28(10):1221-7 | **3678950** | 10.1136/gut.28.10.1221 | PRIMARY **direct luminal** (not fecal) SCFA anchor — autopsy sampling, sudden-death victims, <4h post-mortem. Total SCFA (mmol/kg): 13±6 terminal ileum, **131±9 caecum**, **80±11 descending colon**. Blood (µmol/L): portal 375±70, hepatic 148±42, peripheral 79±22. "Acetate was the principal anion." One of den Besten 2013's own 3 cited sources for the 60:20:20 ratio (confirmed via live full-text reference-list fetch). |
| 5 | den Besten, van Eunen, Groen, Venema, Reijngoud, Bakker 2013, *J Lipid Res* 54(9):2325-40 | **23821742** | 10.1194/jlr.R036012 | PRIMARY SCFA-quantity review, full text (PMC3735932) fetched live. "Acetate, propionate, and butyrate are present in an approximate molar ratio of **60:20:20** in the colon and stool." "total concentration of SCFAs decreases from **70 to 140 mM in the proximal colon** to **20 to 70 mM in the distal colon**." "**95%** of the produced SCFAs are rapidly absorbed by the colonocytes while the remaining **5%** are secreted in the feces." "Fermentation ... yield **400-600 mmol SCFAs/day** ... equivalent to ~10% of the human caloric requirements." Fecal secretion: 10-30 mmol/day (high-fiber) vs 5-15 mmol/day (control diet). |
| 6 | Bergman 1990, *Physiol Rev* 70(2):567-90 | **2181501** | 10.1152/physrev.1990.70.2.567 | PRIMARY independent caloric-contribution + cross-species ratio anchor. "produced in a ratio varying from approximately **75:15:10 to 40:40:20**." "VFA contribute approximately 70% to the caloric requirements of ruminants ... approximately **10% for humans**, and approximately 20-30% for several other omnivorous or herbivorous animals." |
| 7 | Roediger 1980, *Gut* 21(9):793-8 | **7429343** | 10.1136/gut.21.9.793 | PRIMARY **human** colonocyte-fuel anchor (isolated colonocytes, ascending n=7 + descending n=7). Butyrate alone: **73% (ascending) / 75% (descending)** of O2 consumption. WITH 10mM glucose: **59%/72%** — forced adversary, dominance attenuated not abolished. Glucose alone: 85%/30%. |
| 8 | Roediger 1982, *Gastroenterology* 83(2):424-9 | **7084619** | (none on record) | Cross-species (RAT) replication. Butyrate (10mM) alone = **86%** of total O2 consumption, suppressed endogenous-fuel oxidation by 82%; preference order butyrate>acetate>propionate; glucose alone=30%. |
| 9 | Clausen, Mortensen 1994 (rat), *Gastroenterology* 106(2):423-32 | **8299908** | 10.1016/0016-5085(94)90601-7 | DECORRELATED METHOD (Km kinetics, not %-O2). Km order: butyrate (0.184mM) < propionate (0.339) < acetate (0.487) < glucose (0.777) — butyrate has the highest affinity. |
| 10 | Clausen, Mortensen 1995 (human), *Gut* 37(5):684-9 | **8549946** | 10.1136/gut.37.5.684 | Human replication of the Km-ordering result; no significant Vmax/Km difference by UC-disease-status (n=14 UC + n=8 control). |
| 11 | Donohoe et al 2011, *Cell Metab* 13(5):517-26 | **21531334** | 10.1016/j.cmet.2011.02.018 | DECORRELATED METHOD+SPECIES (mouse, transcriptomic/metabolic-state, not O2-tracer or kinetics). Germ-free colonocytes = energy-deprived (↓NADH/NAD+, OXPHOS, ATP → AMPK → autophagy); butyrate rescues, "acting as an energy source rather than as an HDAC inhibitor." |
| 12 | Bäckhed et al 2004, *PNAS* 101(44):15718-23 | **15505215** | 10.1073/pnas.0407076101 | PRIMARY causal germ-free-vs-colonized anchor (binary presence/absence). "conventionalization of adult germ-free (GF) ... mice ... produces a **60% increase in body fat content** and insulin resistance within 14 days **despite reduced food intake**" — rules out an "ate more" confound. |
| 13 | Turnbaugh et al 2006, *Nature* 444(7122):1027-31 | **17183312** | 10.1038/nature05414 | Task's own named citation. Live re-read: **no percentage-of-host-energy figure** in the abstract — a donor-**composition**-dependent result: "colonization of germ-free mice with an obese microbiota results in a significantly greater increase in total body fat than ... a lean microbiota." Qualitative corroboration only (Section D, OODA Fix #1). |
| 14 | Forslund et al 2015, *Nature* 528(7581):262-6 | **26633628** | 10.1038/nature15766 | Freshly-verified (not reused) confound anchor, n=784 metagenomes. "antidiabetic medication confounds" T2D-microbiome association; controlling for metformin, unified T2D signature = "depletion of butyrate-producing taxa." |
| 15 | Human Microbiome Project Consortium 2012, *Nature* 486(7402):207-14 | **22699609** | 10.1038/nature11234 | Inter-individual variation anchor. "healthy individuals differ remarkably in the microbes that occupy habitats ... vary widely even among healthy subjects" BUT "metagenomic carriage of metabolic pathways was stable ... despite variation in community structure." |

## 3. Section A — Bacterial census: colonic density, total count, and the corrected ~1:1 ratio

**Inputs** (Sender/Fuchs/Milo 2016, PMID 27541692): colonic bacterial density **10^11 per mL**
content (their Table 1); colon content volume **0.4 L** (400 mL, their Box 1); their own
independently-stated headline totals: **3.8×10^13 bacteria**, **3.0×10^13 human cells**
(hematopoietic lineage ≈90% of the human-cell count), **B/H = 1.3** (±25% uncertainty, 53%
population CV), total bacterial mass **0.2 kg** (≈0.3% of body weight, updating prior "1-3 kg /
1-3%" claims).

**Arithmetic cross-check (F2, machine-computed, not eyeballed)**: density × volume = 10^11/mL ×
400 mL = **4.0×10^13**, vs. the independently-stated headline 3.8×10^13 — **relative error 5.26%**,
well inside the pre-registered 20% tolerance. This is a genuine internal-consistency check (the
paper's own Table-1 density and Box-1 volume numbers reproduce its own separately-stated headline
total), not a tautology — the two numbers come from different parts of the source paper.

**Ratio gates (F1, the primary falsifier's first conjunct)**: B/H = 1.3 falls inside the
pre-registered "same order of magnitude as 1:1" band **[0.3, 3.3]** — **PASS** — and explicitly
**excludes** the pre-registered old-myth band **[7, 13]** — **PASS**. Both required for the
refutation claim; both hold.

**Forced adversary — is the "corrected" 1:1 figure itself just an equally arbitrary re-guess?**
Falls, on three independent grounds:
1. **Cross-validated against 12 independent empirically-measured stool-density studies** (Sender's
   own Table 2, spanning 1966-2008, multiple methods): range **0.35-3.2×10^11/g wet stool** — the
   10^11 order-of-magnitude figure sits squarely inside this empirically-measured range, not resting
   on the back-of-envelope estimate alone.
2. **Myth-mechanism reconstruction**: replaying Luckey (1972)'s own stated arithmetic — the *same*
   density (10^11/g) × the *wrong* volume (1 L = 1000 mL, whole alimentary tract) — gives
   **1.0×10^14**, matching the ~10^14 figure Sender's own paper identifies (their refs [2],[3]) as
   the traceable historical source of the pre-2016 consensus. The correction is a **volume** fix,
   not a density re-guess.
3. **The volume fix is anatomically grounded, not arbitrary**: upper-GI bacterial density (Sender's
   own Table 1) is **10^3-10^4/mL** (stomach, duodenum, jejunum) to **10^8/mL** (ileum) — **3 to 8
   orders of magnitude below** the colon's 10^11/mL, because gastric acid, bile, and fast transit
   genuinely suppress bacterial growth upstream. Restricting the integration volume to the colon is
   a real anatomical fact, independently measured, not a refit chosen to produce the desired answer.

**Verdict: Section A PASS, 5/5 gates** (`census.*` in the evidence JSON).

## 4. Section B — SCFA molar ratio, total production, luminal concentration, absorption closure

**Molar ratio (F3)**: den Besten 2013's own quoted convention, **60:20:20** (acetate:propionate:
butyrate), citing Cummings 1987 among 3 sources. Sums to 100 exactly (trivial simplex check, PASS).
Checked component-wise against Bergman 1990's **independent**, cross-species/cross-diet reported
band (endpoints 75:15:10 to 40:40:20): acetate 60∈[40,75] ✓, propionate 20∈[15,40] ✓, butyrate
20∈[10,20] ✓ — **all 3 components PASS**, though butyrate sits **exactly at the reported upper
edge** of Bergman's band (disclosed, not hidden — the two source papers agree closely but not with
daylight to spare on butyrate specifically).

**Total production (F5-lite) — a disclosed correction to the task's own band**: den Besten 2013
states **400-600 mmol SCFA/day** (central 500), citing Bergman 1990 among its sources. The task's
own pre-registered band ("~50-100+ mmol/day") **undershoots this by a factor of ~5x** at its own
upper edge. The lenient "+" qualifier is technically not violated (500 ≥ 50 — gate PASSES), but the
headline number reported here is the **measured 400-600 mmol/day**, not the task's own looser band
— stated plainly, not smoothed over by the "+".

**Luminal concentration (F4)**: task's own "~100 mM" point estimate is bracketed by Cummings 1987's
**direct** autopsy measurement (131 mmol/kg caecum, 80 mmol/kg descending colon) and by den Besten's
independent **70-140 mM proximal / 20-70 mM distal** band. Task's ~100mM falls inside the **union**
range [20,140] (PASS) and inside the **proximal-only** range [70,140] (PASS) — but **NOT** inside
the **distal-only** range [20,70] (deliberate FAIL, the one designed-to-fail gate in this document:
`scfa.task_100mM_within_distal_range_ONLY`). This is a genuine, disclosed **regional gradient**, not
a defect: "~100mM" is an accurate single-point summary of the *proximal* colon only; the *distal*
colon runs measurably lower. Reporting a single scalar (as the task does) silently discards this
spatial structure — flagged here, not smoothed over.

**Absorption mass-balance closure (F5, an over-determination check)**: den Besten independently
states (a) production 400-600 mmol/day, (b) fecal excretion 5-15 (control diet) or 10-30 (high-fiber
diet) mmol/day, and (c) "95% absorbed / 5% fecal" — three numbers that were each separately reported
and must be arithmetically mutually consistent under simple mass conservation. Computing
1 − fecal/production from (a)+(b) alone: **control diet → 98% implied absorption; high-fiber diet →
96%** — both within **3 percentage points** of the independently-stated 95% (pre-registered tolerance
5 points) — **PASS**. This is a genuine over-determination check (3 independently-measured/-stated
quantities, 1 conservation law), not a tautology.

**Forced adversary (does this document commit the fecal-proxy error its own repo already flags?)**
No: `MICROBIOME-HOST-METABOLIC` (existing graph cell, read-only) already flags "fecal SCFA
concentration ≠ production flux" (re: Schwiertz 2010). This document's own headline numbers are
**not** fecal-proxy — Cummings 1987 is direct autopsy intestinal-content sampling (not voided
stool), and den Besten/Bergman's 400-600 mmol/day + 95%-absorbed figures are explicit production-
flux estimates. Fecal SCFA (the ~5% non-absorbed remainder) is reported here **separately**
(Section 6) and never substituted for production.

**Verdict: Section B PASS on the primary-falsifier gates (ratio, luminal-union) + absorption-closure
gate; 1 deliberate regional-specificity FAIL disclosed, not a defect.**

## 5. Section C — Butyrate as the dominant colonocyte fuel (~70%)

**Primary human anchor** (Roediger 1980, isolated human colonocytes, ascending n=7 + descending
n=7): butyrate alone accounts for **73% (ascending) / 75% (descending)** of O2 consumption — mean
**74%**, within **4 percentage points** of the task's own "~70%" figure (pre-registered tolerance 10
points) — **PASS**.

**Forced adversary — is this an artifact of a substrate-starved isolated-cell assay?** Tested
directly by Roediger's own experiment: adding 10mM glucose (a competing substrate) drops the
butyrate-attributable fraction to **59% (ascending) / 72% (descending)** — mean 65.5%, an **11.5%**
relative attenuation — but dominance **survives** the pre-registered ≥50% floor at both sites —
**PASS**. Butyrate's primacy is not merely an artifact of substrate starvation; it persists under
competition, just measurably reduced.

**Cross-species replication** (Roediger 1982, rat): butyrate alone = **86%**, ratio-to-human-mean
**1.16x** — same order, rat somewhat higher. Reported as a range (73-86% across species), not forced
to one knife-edge number.

**Decorrelated-method corroborations** (not double-counted as numeric replicates):
- **Kinetics** (Clausen & Mortensen 1994 rat + 1995 human): butyrate has the **lowest Km** (highest
  affinity) among acetate/propionate/butyrate/glucose in both species (rat: butyrate 0.184mM <
  propionate 0.339 < acetate 0.487 < glucose 0.777) — supports preferential use at physiological
  sub-saturating concentrations, via a completely different metric (affinity, not %-of-O2).
- **Genetic/metabolic-state** (Donohoe 2011, mouse): germ-free colonocytes are measurably
  energy-deprived (↓NADH/NAD+, oxidative phosphorylation, ATP → AMPK activation → autophagy);
  butyrate specifically rescues this state, acting as an energy source (not via its separate
  HDAC-inhibitor role) — a qualitative, mechanism-level confirmation via transcriptomic/metabolic
  readout, orthogonal to both the O2-consumption-fraction and the kinetics evidence above.

**Verdict: Section C PASS, 2/2 gates**, plus 2 non-gated decorrelated qualitative corroborations.

## 6. Section D — Decorrelated caloric-harvest check (germ-free vs. colonized)

**OODA Fix #1 (a catch in the task's own framing, diagnosed and fixed, not hidden)**: the task text
bundles "~10% of host energy from SCFA" together with "Turnbaugh 2006 germ-free-transplant
adiposity" as if a single finding. Live re-read of Turnbaugh et al 2006 (PMID 17183312)'s own
abstract: it reports **no percentage-of-host-energy figure at all** — it is a donor-**composition**-
dependent result (obese-donor microbiota → more adiposity than lean-donor microbiota, same diet).
The "~10%" figure is **Bergman 1990's separate**, independently-methodologized (species-comparative
mass-balance) estimate. **Fixed**: reported as two separate citations below, not merged.

**OODA Fix #2 (a near-double-count, caught by tracing citations, not asserted)**: den Besten 2013's
own "~10% of human caloric requirements" sentence was nearly counted as a **third** independent
corroborating source alongside Bergman 1990 and Bäckhed 2004. A live full-text fetch of den Besten's
own reference list shows its citation for that number **is Bergman 1990** (den Besten's own ref
[31]) — i.e., den Besten **repeats** Bergman's estimate, it does not independently re-derive it.
**Fixed**: this document counts only **2** genuinely decorrelated angles, not 3 — avoiding a
shared-source evidence-inflation false positive that would have made the triangulation look
stronger than it is.

**The 2 genuinely decorrelated angles**:
1. **Direct causal manipulation** (Bäckhed et al 2004, PMID 15505215): conventionalizing adult
   germ-free mice with a normal cecal microbiota produces a **60% increase in body fat** within 14
   days **despite reduced food intake** — ruling out an "ate more" confound; mechanism = increased
   monosaccharide absorption + hepatic de novo lipogenesis. Gate: ≥10% (non-trivial-effect floor) —
   **PASS**, and by a wide margin (60% vs. a 10% floor).
2. **Independent mass-balance/stoichiometric estimate** (Bergman 1990, PMID 2181501): SCFA
   contribute **~10%** of human daily caloric requirements (vs. ~70% ruminants, ~20-30% other
   omnivores/herbivores — humans are a genuine outlier, low but non-zero). Gate: within the
   task-adjacent band [5,20]% — **PASS**.

These two angles are methodologically about as decorrelated as this literature gets (an *in vivo*
causal knockout-style manipulation in mice vs. a cross-species stoichiometric energy-balance
estimate in humans) and converge on the same qualitative conclusion — microbiome-derived calories
are real and quantitatively non-trivial — **without being the same measurement counted twice.**

**Verdict: Section D PASS, 2/2 gates**, with the double-count avoided explicitly disclosed.

## 7. Section E — Symmetric QC, held OPEN per task instruction (not resolved away)

- **Individual variation is large, not resolved to a point estimate**: Sender 2016's own B/H ratio
  has **53% population CV**; its own Table 2 shows **24-78% CV** *within each* of 12 independent
  stool-density studies. HMP 2012 (PMID 22699609, freshly verified this session): microbial
  community **structure** "vary[ies] widely even among healthy subjects" — but metagenomic
  **metabolic-pathway carriage was stable** despite that structural variation. This document's
  numbers are population medians/typical orders of magnitude, **not** individual predictions.
- **Diet-driven, not fixed**: den Besten's own numbers show fecal SCFA secretion swings **~2-3x**
  by diet (10-30 mmol/day high-fiber vs. 5-15 mmol/day control), and the luminal-concentration band
  itself is diet-dependent.
- **Most disease associations are correlational, not causal, and often drug/diet-confounded**:
  Forslund et al 2015 (PMID 26633628, n=784 metagenomes, freshly verified this session, not reused
  from a sibling cell) — a widely-reported T2D-microbiome "signature" from two prior unstratified
  studies was substantially a **metformin drug effect**, not a disease effect, until stratified. This
  document does **not** extend its own census/SCFA/fuel numbers into any disease-association or
  causal claim.
- **"Psychobiotic" / gut-brain causal-mood overreach is NOT re-litigated here** — already covered,
  more skeptically and at greater depth, by this repo's own existing (read-only, **not** re-verified
  this session) graph cells: `SYS-GUT-MICROBIOME-HOST-AXIS` flags PMID 29197739's null
  psychobiotic-mood RCT meta-analysis as open/contested; `MOL-MICROBIOME-GUT-BRAIN-HUB`'s own Edge C
  is explicitly "LOW confidence ... explicitly marked hypothesis"; `MOL-MICROBIOME-METAGENOMICS-
  HOST-COUPLING`'s own IBD/CGM anchor is only MARGINAL (AUROC 0.5825, z=1.45). This document is a
  narrower, more basic quantitative layer underneath all of them, and does not contradict or
  re-derive their (more cautious) conclusions.
- **SCFA-fecal-proxy caveat — this document clears the bar the wider field often misses, but the
  general warning still stands**: this repo's own `MICROBIOME-HOST-METABOLIC` cell already flags
  "fecal SCFA concentration ≠ production flux" (re: Schwiertz 2010) as a field-wide hazard. This
  document's own core numbers are **not** fecal-proxy (Section 4), but most published cohort-level
  "SCFA level" claims in the wider literature *are* fecal-concentration stand-ins — the caveat
  remains valid as a general field warning even though this document's own numbers clear it.

## 8. Couplings (prose only — not folded into the graph this session)

Per this session's isolation scope (another instance writes `data/MECHANISM_ANCHOR_GRAPH.json`
concurrently), these are **prose-only**, left for the canonical `mechanism_fold.py` path:

- **GI (fermentation substrate)** / `docs/MECHANISM_GI_ABSORPTION_TRANSIT.md`: that document's
  small-bowel-transit layer (and its own disclosed gap — colonic transit is not modeled there
  either) sets the rate at which fermentable substrate (resistant starch, non-starch polysaccharide)
  reaches the colon. This document takes Western-diet fiber intake (20-25 g/day, den Besten) as a
  disclosed fixed input, not itself modeled — a genuine open extension for a future pass.
- **Immune (SCFA → Treg)**: mechanistically causal in gnotobiotic mice (Furusawa 2013, Smith 2013
  GPR43/FFAR2 pathway) — already live-verified by the existing, read-only `SYS-GUT-MICROBIOME-HOST-
  AXIS` graph cell; **not** re-verified in this session, cross-referenced only.
- **Metabolic (caloric harvest)**: Section 6 above — Bäckhed 2004's causal germ-free manipulation
  (+60% body fat/14d) and Bergman 1990's independent ~10%-of-human-calories mass-balance estimate,
  with the near-triple-counted third angle (den Besten's repeat of Bergman) explicitly caught and
  excluded (OODA Fix #2).

## 9. Honest gaps (disclosed, not hidden)

- No subject-specific data anywhere — population-parametrized census/production arithmetic, matching
  every other MECHANISM_* systemic layer's own disclosed scope.
- The 10:1 myth's correction is verified here on the **numerator** side (bacteria: ~10^14 → 3.8×10^13,
  a volume-misattribution error, directly reconstructed). The myth's **denominator** side (exactly
  what "old" human-cell figure combined with the old ~10^14 bacteria estimate to produce precisely
  "10:1") is not independently re-derived from one single live-verified primary source this
  session — Sender 2016's own explicit statement that it "updates the widely-cited 10:1 ratio" is
  relied on directly, rather than reconstructing every historical step of the denominator's own
  provenance.
- SCFA total-production task band ("~50-100+ mmol/day") undershoots the directly-verified literature
  figure (400-600 mmol/day) by ~5x at the task's own upper edge — disclosed prominently (Section 4),
  not silently absorbed by the lenient "+" qualifier.
- Colonocyte-fuel percentages (Roediger 1980/1982) come from **isolated-cell suspensions ex vivo**,
  not in-vivo flux measurement — the forced glucose-competition adversary (Section 5) partially
  addresses the "substrate-starved-assay-inflates-butyrate" concern but does not fully substitute
  for an in-vivo tracer study (none was found live this session with a directly quotable
  colonocyte-specific in-vivo % figure).
- Both headline numbers (bacteria:host ratio, SCFA figures) are population-level literature
  syntheses — individual variation is large and diet-dependent (Section 7), reported, not resolved.

## 10. Confidence tier

**In-vivo-anchored**: Sender/Fuchs/Milo 2016 census (PMID 27541692, itself cross-validated against
12 independent empirical stool-density studies in their own Table 2); Cummings 1987 direct autopsy
luminal SCFA measurement (PMID 3678950); Roediger 1980 human isolated-colonocyte O2-consumption
(PMID 7429343); Bäckhed 2004 causal germ-free-vs-colonized mouse manipulation (PMID 15505215).
Population-parametrized arithmetic, **not** subject-specific — matching `gi_absorption_transit.py` /
`renal_filtration.py`'s own stated confidence tier. Weaker specifically for the total-production
number (review-synthesized, not a single directly-measured flux) and the colonocyte-fuel-%
numbers (ex-vivo isolated-cell, not in-vivo) — both disclosed in Section 9, not hidden behind the
overall tier label.

## 11. Gates (machine-computed, `reports/probes/gut_microbiome.json`)

```
PRIMARY FALSIFIER (task's own pre-registered conjunction) -- PASS, 4/4:
  census.ratio_same_order_as_1to1:                    PASS (B/H=1.3, band [0.3,3.3])
  census.ratio_excludes_10to1_myth_band:              PASS (1.3 excludes [7,13])
  scfa.ratio_60_20_20_within_bergman_crossspecies_band: PASS (all 3 components; butyrate at edge)
  scfa.task_100mM_within_luminal_union_range:         PASS (100 in [20,140])

SECONDARY CORROBORATION (decorrelated, not gated as harshly, all reported) -- PASS, 6/6:
  census.density_x_volume_reproduces_headline_lt20pct: PASS (5.26% rel. err.)
  scfa.absorption_arithmetic_closes_within_5pts:       PASS (max dev 3.0 pts vs stated 95%)
  fuel.human_mean_within_10pts_of_task_70pct:          PASS (74% vs 70%, diff 4pts)
  fuel.dominance_survives_glucose_competition:         PASS (65.5% mean, floor 50%)
  harvest.backhed_germfree_nontrivial_ge10pct:         PASS (60% vs floor 10%)
  harvest.bergman_estimate_within_5_20pct_band:        PASS (10% in [5,20])

ADDITIONAL DIAGNOSTIC GATES -- 5/6 (1 deliberate, disclosed, non-defect FAIL):
  census.density_cross_validated_by_12_independent_studies: PASS
  census.upper_gi_density_suppressed_ge3_orders:       PASS (3-8 orders)
  scfa.ratio_sums_to_100:                              PASS (trivial simplex check)
  scfa.production_central_ge_task_floor:               PASS (500 >= 50, though 5x above ceiling)
  scfa.task_100mM_within_proximal_range:               PASS (100 in [70,140])
  scfa.task_100mM_within_distal_range_ONLY:            FAIL (100 > 70; regional gradient, disclosed
                                                         Section 4 -- NOT a defect, a designed probe)

overall gates: 15/16 PASS (the 1 FAIL is intentional/informational, see above)
verdict.primary_falsifier_ratio_and_scfa_conjunction: PASS
verdict.secondary_corroboration_fuel_and_harvest:     PASS
verdict.overall:                                      PASS
```

## 12. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/gut_microbiome.py
```
Pure Python/numpy, no OpenSim, no subject data, runs in <1s. Writes
`reports/probes/gut_microbiome.json`. No git operations. Files touched this session:
`scripts/msk/gut_microbiome.py`, `reports/probes/gut_microbiome.json`, this doc.
`data/MECHANISM_ANCHOR_GRAPH.json` (shared, concurrently written by another instance) was read-only
referenced (Section 0), never edited. `data/body_twin/agent_outputs/*.json` (prior agents' passes)
were read-only referenced for context/non-duplication and to catch one prior-session PMID-adjacent
claim (Bäckhed sub-breakdown) that this session's own live fetch could not independently confirm —
never edited.

## 13. Files

- `scripts/msk/gut_microbiome.py` — full model: 15 citations (dict with PMID/DOI/verbatim role),
  census arithmetic (density×volume cross-check, myth reconstruction, upper-GI suppression), SCFA
  arithmetic (ratio-band check, production-vs-task disclosure, luminal bracketing, absorption
  mass-balance closure), colonocyte-fuel forced-adversary (glucose competition), caloric-harvest
  decorrelated triangulation (with the OODA double-count catch), symmetric-QC prose, all gates.
- `reports/probes/gut_microbiome.json` — full evidence: every citation, every computed number, all
  16 gates, the two-tier verdict, honest gaps, confidence tier.
- Read-only, not modified: `data/MECHANISM_ANCHOR_GRAPH.json` (Section 0 disambiguation only),
  `data/body_twin/agent_outputs/microbiome-gut-brain-hub__a857dea352c10e12a.json` and
  `gut-microbiome-brain-axis__a4bc09700fa51fe49.json` (context/non-duplication only),
  `reports/probes/mt_microbiome.json` (context only).
