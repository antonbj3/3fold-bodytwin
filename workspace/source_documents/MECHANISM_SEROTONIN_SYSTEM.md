# MECHANISM SEROTONIN (5-HT) SYSTEM — synthesis (TPH1/TPH2), SERT reuptake, receptor diversity, the SSRI paradox (2026-07-22)

Script: `scripts/msk/serotonin_system.py`. Raw results: `data/serotonin_system/serotonin_system_results.json`.
Evidence (citations): `docs/MECHANISM_SEROTONIN_SYSTEM_evidence.json`. Raw fetch logs (all esearch/esummary/efetch
calls this session, verbatim, plus one PMC full-text pull): `data/raw_fetch/serotonin_*.{json,txt,xml}`.

**Provenance note:** no pre-existing serotonin/5-HT/SERT script or doc was found anywhere in this repo
(repo-wide grep before starting) — this is the **first serotonin-system doc**, written fresh this session,
no collision with concurrent work.

## 0. Why this layer + couples_to

Couples **nerve conduction** (`docs/MECHANISM_NERVE_CONDUCTION.md` — a different neurotransmitter system, same
conduction/transmission-kinetics layer; grep-confirmed zero serotonin content there), **dopamine/monoamine**
(`docs/MECHANISM_DOPAMINE_KINETICS.md` — a genuinely load-bearing coupling, not a topic label: Part 1 below
explicitly contrasts TPH's near-Km substrate-sensitive operation against that doc's TH near-saturation
treatment, **the same governing Michaelis-Menten hyperbola, different physiological operating point**),
**gut/enteric 5-HT** (`docs/MECHANISM_GUT_MICROBIOME.md` — a *disjoint* mechanism, that doc models bacterial
census/SCFA/butyrate fermentation, this doc models host enterochromaffin-cell TPH1 synthesis; grep-confirmed
zero serotonin content there, genuinely complementary not duplicative), and **psychiatric**
(`data/MECHANISM_ANCHOR_GRAPH.json` nodes `PSYCH-DEPRESSION-MECHANISM`, `PSYCH-KYNURENINE-NEUROTOX-BRANCHPOINT`,
`PSYCH-INFLAMMATORY-SUBTYPE-DIMENSION` — all OPEN/autonomous-expansion, operating at the **clinical-
subtype/biomarker** level; **this doc is a distinct, lower layer**: the actual synthesis/reuptake/receptor
mechanism any such clinical model must assume as its physical substrate — the same relationship
`dopamine_kinetics.py` has to `NEU-REWARD-DOPAMINE`). Checked read-only, none edited.

## 1. Citations — 24 PMIDs, every one verified LIVE this session (NCBI eutils esearch/esummary/efetch + 1 PMC full-text pull)

| # | Citation | PMID | Role | Verified this session |
|---|---|---|---|---|
| 1 | Walther DJ et al (2003). *Science* 299(5603):76. | 12511643 | TPH2 discovery | Title/DOI/PMID (no abstract, Brevia format) |
| 2 | Côté F et al (2003). *PNAS* 100(23):13525-30. | 14597720 | **Tph1-KO, peripheral pool** | Full abstract, verbatim |
| 3 | Alenina N et al (2009). *PNAS* 106(25):10332-7. | 19520831 | **Tph2-KO, CNS pool** | Full abstract, verbatim |
| 4 | Savelieva KV et al (2008). *PLoS ONE* 3(10):e3301. | 18923670 | **DKO, decisive dissociation** | Full abstract, verbatim |
| 5 | Fernstrom JD, Wurtman RJ (1971). *Science* 173(3992):149-52. | 5581909 | **Trp-availability sensitivity** | Full abstract, verbatim |
| 6 | Lovenberg W et al (1967). *Science* 155(3759):217-9. | 6015530 | TPH activity assay | Full abstract, verbatim |
| 7 | Blakely RD et al (1991). *Nature* 354(6348):66-70. | 1944572 | SERT cloning | Full abstract, verbatim |
| 8 | Meyer JH et al (2004). *Am J Psychiatry* 161(5):826-35. | 15121647 | **PET occupancy, PRIMARY falsifier** | Full abstract, verbatim |
| 9 | Meyer JH et al (2001). *Am J Psychiatry* 158(11):1843-9. | 11691690 | PET occupancy, corroborating | Full abstract, verbatim |
| 10 | Hoyer D, Hannon JP, Martin GR (2002). *Pharmacol Biochem Behav* 71(4):533-54. | 11888546 | **7 receptor families** | Full abstract, verbatim |
| 11 | Hoyer D et al (1994). *Pharmacol Rev* 46(2):157-203. | 7938165 | Classification history (symmetric QC) | Abstract (truncated 400w) |
| 12 | Blier P, de Montigny C (1994). *Trends Pharmacol Sci* 15(7):220-6. | 7940983 | Context (thin abstract) | Full abstract, verbatim |
| 13 | Invernizzi R et al (1992). *Brain Res* 584(1-2):322-4. | 1515949 | **Raphe/cortex forced adversary** | Full abstract, verbatim |
| 14 | Artigas F et al (1996). *Trends Neurosci* 19(9):378-83. | 8873352 | **Desensitization mechanism, PRIMARY** | Full abstract, verbatim |
| 15 | Pérez V et al (2001). *J Clin Psychopharmacol* 21(1):36-45. | 11199945 | **Pindolol RCT, clinical anchor** | Full abstract, verbatim |
| 16 | Le Poul E et al (1995). *Naunyn Schmiedebergs Arch Pharmacol* 352(2):141-8. | 7477436 | **Desensitization timecourse, PRIMARY** | Full abstract, verbatim |
| 17 | Gershon MD, Tack J (2007). *Gastroenterology* 132(1):397-414. | 17241888 | Gut mechanism (no "90%" in abstract) | Full abstract, verbatim |
| 18 | Yano JM et al (2015). *Cell* 161(2):264-76. | 25860609 | **">90%" gut serotonin, PRIMARY** | Abstract + **PMC full text** |
| 19 | O'Mahony SM et al (2015). *Behav Brain Res* 277:32-48. | 25078296 | Gut-brain-microbiome review | Full abstract, verbatim |
| 20 | Fernstrom JD (2013). *Amino Acids* 45(3):419-30. | 22677921 | **BBB/LNAA mechanism** | Full abstract, verbatim |
| 21 | Ruhé HG et al (2007). *Mol Psychiatry* 12(4):331-59. | 17389902 | **Tryptophan-depletion, PRIMARY** | Full abstract, verbatim |
| 22 | Moncrieff J et al (2022). *Mol Psychiatry* 28(8):3243-3256. | 35854107 | **Symmetric QC, PRIMARY** | Full abstract, verbatim |
| 23 | Bartová L et al (2023). *Mol Psychiatry* 28(8):3153-4. | 37322062 | Rebuttal (contested status) | Title/DOI/PMID (brief reply, no abstract) |
| 24 | Rudnick G (1977). *J Biol Chem* 252(7):2170-4. | 849926 | **SERT Km=0.5µM, PRIMARY** | Full abstract, verbatim |

**Recall-drift caught live this session** (not from memory, per this repo's own ~62%-drift warning):
1. **Meyer 2004** — first esearch returned the *2001* paper (PMID 11691690) by title/topic match; caught by
   esummary date-checking, then separately found the correct 2004 five-SSRI paper (PMID 15121647).
2. **Fernstrom & Wurtman** — commonly mis-cited as 1972; live-verified actual date is **1971** Jul 9.
3. **Le Poul** — initially mis-recalled as "le Poul 2000"; live esearch found the correct paper is **1995**.
4. **Gershon & Tack 2007's own abstract does NOT contain the famous "90%" figure** — checked directly by
   reading the fetched text rather than assuming the citation supports the number widely attributed to it.
   The actual verbatim source (found by pulling **PMC full text**, not stopping at the abstract) is Yano et
   al 2015: *"More than 90% of the body's 5-HT is synthesized in the gut."*

**Verbatim quotes carrying the decisive numbers:**

- **Meyer et al 2004** (PMID 15121647): *"Minimum therapeutic doses of paroxetine and citalopram produce 80%
  occupancy for the serotonin (5-HT) transporter (5-HTT)... Mean occupancy at this dose was 76%-85%...
  Occupancy of 80% across five SSRIs occurs at minimum therapeutic doses... 80% 5-HTT blockade is important
  for therapeutic effect."* n=77, striatal [11C]DASB PET.
- **Meyer et al 2001** (PMID 11691690, independent n=12 cohort): paroxetine 20mg/day → **83%** occupancy;
  citalopram 20mg/day → **77%** occupancy — converges on the same ~80% figure as the larger 2004 study.
- **Le Poul et al 1995** (PMID 7477436): *"the potency of the 5-HT1A receptor agonist, 8-OH-DPAT, to depress
  the firing of serotoninergic neurons... was significantly reduced as early as after a 3-day treatment...
  The proportion of recorded neurons showing desensitization... increased along the treatment from
  approximately 40% on the 3rd day to 60-80% on the 21st day."* Also: specific receptor **binding** was
  unmodified at any timepoint — functional/coupling desensitization, not a receptor-density change.
- **Invernizzi et al 1992** (PMID 1515949): citalopram 1mg/kg i.p. *"significantly increased dialysate
  serotonin in the dorsal raphe, but not in the frontal cortex"* — UNLESS the raphe autoreceptor is blocked
  (methiothepine infusion, which *"by itself did not change cortical serotonin concentrations"*), at which
  point *"the same 1mg/kg citalopram dose significantly increased the extracellular concentration of
  serotonin in the frontal cortex."*
- **Artigas et al 1996** (PMID 8873352): *"SSRIs... increase the extracellular concentration of 5-HT in the
  midbrain raphé nuclei, thereby activating inhibitory somatodendritic 5-HT1A autoreceptors. Consequently,
  the firing activity of 5-HT neurons is reduced and the enhancement of extracellular 5-HT concentration in
  forebrain is dampened."*
- **Pérez et al 2001** (PMID 11199945): *"Median times to sustained response were 19 days for fluoxetine plus
  pindolol (N=55) and 29 days for fluoxetine plus placebo (N=56) (p=0.01)"*; matched-responder subset **18
  vs 10 days, p=0.0002**.
- **Savelieva et al 2008** (PMID 18923670): *"dramatically reduced central 5-HT levels in Tph2 knockout
  (TPH2KO) and Tph1/Tph2 double knockout (DKO) mice; and substantially reduced peripheral 5-HT levels in
  DKO, but NOT TPH2KO mice."*
- **Côté et al 2003** (PMID 14597720): *"the neuronal tph2 is expressed in neurons of the raphe nuclei and of
  the myenteric plexus, whereas the nonneuronal tph1... is in the pineal gland and the enterochromaffin
  cells."* Tph1-/- mice: *"larger heart sizes... abnormal cardiac activity, which ultimately leads to heart
  failure"* — a **peripheral (cardiac)** phenotype, not a CNS one.
- **Alenina et al 2009** (PMID 19520831): Tph2-/- mice *"lack serotonin in the central nervous system...
  growth retardation and 50% lethality in the first 4 weeks... These data confirm that the majority of
  central serotonin is generated by TPH2."*
- **Yano et al 2015** (PMID 25860609, PMC4393509 full text): *"More than 90% of the body's 5-HT is
  synthesized in the gut."*
- **Rudnick 1977** (PMID 849926): *"Kinetic studies reveal a Km of 0.5 microM for the [SERT] transport
  process"* (human platelet plasma-membrane vesicles).
- **Ruhé, Mason, Schene 2007** (PMID 17389902): *"5-HT or NE/DA depletion did not decrease mood in healthy
  controls... In drug-free patients with MDD in remission, a moderate mood decrease was found for ATD...
  ATD induced relapse in patients with MDD in remission who used serotonergic antidepressants... they fail
  to demonstrate a causal relation."*
- **Moncrieff et al 2022** (PMID 35854107): *"The main areas of serotonin research provide no consistent
  evidence of there being an association between serotonin and depression, and no support for the hypothesis
  that depression is caused by lowered serotonin activity or concentrations."* Two of the largest studies
  (genetic association **n=115,257**; collaborative meta-analysis **n=43,165**) *"revealed no evidence of an
  association."* **Nuance, not flattened**: smaller imaging legs (5-HT1A largest n=561; SERT binding largest
  n=1845) showed *"weak and inconsistent evidence of reduced binding... [with] effects of prior
  antidepressant use... not reliably excluded."* Carries **8 published comments** (counted directly), including
  a direct rebuttal (Bartová et al 2023, PMID 37322062) — actively contested, not unanimous.

## 2. Method — geometric, machine-checked (not narrated)

**One governing shape** — the saturating hyperbola `S/(Km+S)` — recurs at three physically distinct scales:
TPH-substrate kinetics (Part 1), SERT-substrate kinetics (Part 2), and SSRI-dose occupancy (Part 2). Its
closed-form **local elasticity**, `Km/(Km+S)`, is exactly 0.5 at `S=Km` and decays toward 0 as `S/Km` grows —
operating near Km is substrate-**sensitive**; operating deep past Km is substrate-**insensitive** (saturated).
This is the *same* hyperbola `dopamine_kinetics.py` uses for DAT/TH — Part 1 below computes the elasticity
contrast directly (5.5× more sensitive at S=Km vs S=10·Km) and anchors the qualitative claim (TPH really
operates near-Km in vivo, unlike TH) via Fernstrom & Wurtman 1971's measured result, not assertion.

**A second, independent saturating shape** — a first-order time-relaxation, `1-exp(-t/tau)` — governs 5-HT1A
autoreceptor desensitization (Part 3), fit directly to le Poul et al 1995's own two reported time points
(day 3 = 40%, day 21 = 60-80%). The **SSRI paradox is a genuinely geometric fact**: it is the clash of two
saturating curves on two different, non-collapsible axes (dose vs time) with very different knee-widths —
occupancy's dose-domain knee is reached within the drug's pharmacokinetic dosing window; desensitization's
time-domain knee sits at days-to-3-weeks. Conflating the two axes (expecting clinical response to track
occupancy 1:1 in time) is exactly the naive model's error.

## 3. Gates — 22/22 pre-registered PASS, but gate STRENGTH is disclosed, not uniform

| Gate | Result |
|---|---|
| G1 Part1 TPH elasticity contrast (S=Km vs S=10Km) ≥5.0× | **5.50×** PASS |
| G2 Part2 SERT synthetic-control Vmax recovery ≤15% | 2.50% PASS |
| G3 Part2 SERT synthetic-control Km recovery ≤15% | 4.63% PASS |
| G4 Part2 MM-vs-exponential functional-form adversary ≥1.5× fold, monotonic | **3.00×**, monotonic PASS† |
| G5 Part2 occupancy hyperbola stays bounded <100% | max 96.97% (at 8×dose) PASS |
| G6 Part2 linear-occupancy adversary exceeds 100% (falsified) | max 640% PASS (adversary correctly falls) |
| G7 Part3 seven receptor families citation-confirmed | Hoyer 2002 verbatim PASS |
| G8 Part3 desensitization τ lands in sane band [1,15]d | 3.56 days PASS‡ |
| G9 Part3 desensitization Dmax lands in sane band [40%,100%] | 70.2% PASS‡ |
| G10 Part3 rodent/human timescale windows overlap | [3,21]d ∩ [10,29]d = **[10,21]d** PASS |
| G11 Part3 Invernizzi autoreceptor-brake required (raw data) | True PASS§ |
| G12 Part3 Invernizzi null (no-brake) model falsified | True PASS§ |
| G13 Part3 pindolol ratio (all responders) ≥1.3× | 1.53× PASS§ |
| G14 Part3 pindolol ratio (matched responders) ≥1.3× | 1.80× PASS§ |
| G15 Part3 pindolol both results significant (p≤0.05) | p=0.01, p=0.0002 PASS§ |
| G16 Part4 genotype×compartment pattern is complementary, not the shared-pool null | True PASS§ |
| G17 Part4 shared-pool null falsified at the decisive cell (Tph2-KO, periphery) | True PASS§ |
| G18 Part4 peripheral-fraction ≥0.85 floor (Yano's >90%) | 0.90 PASS§ |
| G19 Part5 Ruhé group-wise effect monotonic with serotonergic vulnerability | True PASS§ |
| G20 Part5 healthy-control effect < remitted-on-SSRI relapse effect | True PASS§ |
| G21 Part6 both large Moncrieff studies (n=115,257; n=43,165) report null | True PASS§ |
| G22 Part6 occupancy-claim and etiology-claim scopes are logically independent | True PASS¶ |

**† G4 note (bug caught and fixed this session, disclosed, not hidden):** this gate **initially FAILED**
(fold=0.92, non-monotonic) using one fixed absolute time window across a 100× range of C0/Km ratios — the
fixed window silently under-sampled the slow, high-C0 clearance tail. Diagnosed (OODA: observed the actual
per-ratio RMSE values, oriented on the mechanism — window-vs-clearance-timescale mismatch — before touching
code), fixed by making the fit window **per-ratio adaptive** (spanning each ratio's own closed-form
time-to-99%-cleared). Re-ran: fold=3.00, monotonic, PASS. Verified this was a real fix, not a threshold-loosen.

**‡ G8/G9 note:** only 2 literature time points exist for a 2-parameter fit — this exactly reproduces both
points by construction and is **not** an independent validation of the single-first-order functional form; it
is a parameter-**extraction** gated only on landing in a physically sane band (a real check: an inconsistent
input pair could have produced τ<0 or Dmax>100%).

**§ Direct-citation-data gates** (G11-G21): these encode the papers' **own reported** contrasts/numbers as
structured booleans and check they have the claimed shape — they convert narrative into machine-falsifiable
form (a wrong transcription shows up as FAIL) but are **not** a fresh computation the way G1-G10 are. Weaker
evidence than G1-G10, reported honestly, not oversold.

**¶ G22 note:** the single weakest gate — a scope-labeling assertion (occupancy = drug-target engagement;
Moncrieff = baseline etiology, different variables), not a numeric test. Included for explicit auditability
of the non-contradiction argument only.

**Determinism:** 2 independent runs of `scripts/msk/serotonin_system.py` produce byte-identical JSON,
md5sum `d3171308f1f3164e3e00a963a62e97c0`, verified this session.

## 4. Symmetric QC — held OPEN per the task's own explicit instruction, not resolved

- **The "chemical imbalance" theory has no consistent supporting evidence at the largest human sample sizes**
  (Moncrieff et al 2022, two studies n=115,257 and n=43,165 both null) — reported as the headline finding,
  not softened.
- **But Moncrieff's own review does NOT claim uniform-zero-everywhere**: smaller imaging legs show "weak and
  inconsistent" (not flatly null) evidence, confounded by prior antidepressant exposure not reliably excluded
  — steelmanned rather than flattened into a stronger negative than the source supports.
- **Moncrieff 2022 is itself actively disputed** — 8 published comments, a direct rebuttal (Bartová et al
  2023). This doc does not adjudicate that dispute; it reports the dispute's existence.
- **SERT-occupancy (~80%, Meyer 2001/2004) and baseline-serotonin-abnormality-in-depression (Moncrieff) are
  LOGICALLY INDEPENDENT claims** — a drug engaging its target is not evidence the target's baseline state was
  abnormal. Conflating them is exactly the naive-model error this doc must not repeat.
- **The exact "occupancy is fast" (days, not weeks) temporal claim rests on general SSRI pharmacokinetics**
  (steady-state within ~1 week by the standard ~5-half-life rule for most SSRIs), **not** on an independently
  live-fetched serial early-timepoint PET citation this session — a targeted search this session for
  single-dose/early-timepoint SSRI-occupancy PET returned zero PubMed hits. Disclosed as a real,
  structurally-plausible-but-not-independently-re-verified gap, not smoothed over.
- **Pool-separation (Part 4) is mouse-genetic-anchored only** — Cote 2003/Alenina 2009/Savelieva 2008 are all
  from the same overlapping gene-targeting research program (Bader/Vodjdani/Lexicon labs); no human
  causal-manipulation equivalent exists (cannot ethically knock out human TPH1/TPH2). The complementary
  3-genotype structure is strong but not fully dispositive against all conceivable pleiotropy — held OPEN.
- **The 2-point desensitization fit (le Poul 1995) is underdetermined** — 2 parameters, 2 points, exact fit by
  construction. Read as a sane-range parameter extraction, not a validated curve shape.
- **le Poul's own data shows desensitization is functional/coupling-state, not receptor-density** — specific
  radioligand binding was unmodified at every timepoint tested; this doc's "Dmax" state variable should be
  read as "fraction of neurons in a desensitized functional state," not "fraction of receptors internalized."
- **TPH's exact numeric Km for tryptophan was not independently verified live this session** — the geometric
  elasticity argument (Part 1) is citation-independent (depends only on the S/Km ratio); the qualitative
  claim that TPH sits near-Km in vivo is anchored via Fernstrom & Wurtman's functional result, not a Km number.

## 5. Additional convergence found this session

- **A clean 3-genotype double dissociation lives in ONE paper** (Savelieva et al 2008): Tph1-KO alone spares
  the CNS pool; Tph2-KO alone spares the peripheral pool; only the double-KO touches both — a genuinely hard
  pattern for generic off-target pleiotropy to produce by chance, independently corroborating Cote 2003 and
  Alenina 2009's separate single-gene findings from the same broader research program.
- **Two methodologically decorrelated species/methods converge on the SAME 2-4 week window**: le Poul 1995's
  rodent electrophysiological desensitization-maturation window ([3,21] days) and Pérez 2001's human
  clinical-RCT response-latency window ([10,29] days) **overlap** at [10,21] days — different species,
  different instrument (single-unit electrophysiology vs symptom-severity scale), different lab/decade.
- **The Pérez 2001 RCT's placebo-controlled (not untreated-control) design is itself the forced-adversary
  control** for "it's just natural-history/placebo timing" — both arms share that generic confound, so the
  significant between-arm difference isolates the autoreceptor-blockade-specific contribution.
- **Two independent PET cohorts, different years, different N, converge on ~80% occupancy**: Meyer 2001
  (n=12: 83%/77%) and Meyer 2004 (n=77: 76-85%) — not one study's number.

## 6. Honest gaps

- Occupancy-is-fast (days) is inferred from general SSRI pharmacokinetics, not independently re-verified via
  an early-serial-timepoint PET citation this session (Sec. 4).
- SERT's exact Vmax (only Km=0.5µM is citation-verified, Rudnick 1977) is an illustrative, disclosed value —
  same disclosed-gap pattern as `dopamine_kinetics.py`'s own DAT Vmax gap.
- TPH's exact numeric Km for tryptophan is not independently verified live this session (Sec. 4).
- The 2-point le Poul desensitization fit is underdetermined (2 params / 2 points, Sec. 4).
- Pool-separation is mouse-genetic-anchored only, not human-causal-manipulation-anchored (Sec. 4).
- G22 (scope-independence) is a labeling assertion, not a numeric test — the weakest gate in the set.
- This session did not independently re-derive the ODE/curve-fit arithmetic from scratch a second way; one
  real bug (G4's fixed-window artifact) was caught by running the code and inspecting per-ratio output, not
  by a from-scratch parallel implementation.

## 7. Overall result

**22 of 22 pre-registered gates PASS** (11 genuinely computational/forced-adversary gates, one of which
initially FAILED and was fixed after diagnosing the real mechanism — not silently patched; 10 direct-citation-
data gates converting the papers' own reported contrasts into machine-falsifiable form; 1 labeling-assertion
gate), and **all 24 cited PMIDs were verified live this session** (21 with full verbatim abstracts; 3 with
title/DOI/PMID-only confirmation, a disclosed and expected limitation for brief-report/reply formats), finding
**zero fabrications** and **catching 4 live recall-drift instances** before they entered the record (Sec. 1).
The task's two named falsifier anchors were directly confirmed: the SSRI-occupancy paradox (Meyer 2001/2004,
~80-83%, PLUS le Poul 1995's 40%→60-80% desensitization timecourse whose day-21 maturity window overlaps
Pérez 2001's own 10-29 day clinical-response window) and the ~90%-enteric-serotonin fact (Yano 2015 PMC full
text, corroborated by a clean 3-genotype genetic double-dissociation, Cote/Alenina/Savelieva). The
decorrelated tryptophan-depletion check (Ruhé 2007) reproduces the task's own causal-but-conditional framing
exactly (null in healthy controls, relapse-inducing only in remitted patients on serotonergic drugs). The
symmetric-QC leg (Moncrieff 2022) is reported honestly — large-N null, smaller-N nuanced-not-flattened, and
itself contested (8 comments, 1 direct rebuttal) — held OPEN, not resolved in either direction. Confidence
tier: **in-vivo-anchored** (human PET SERT-occupancy + human tryptophan-depletion RCT meta-analysis + a human
pindolol-augmentation RCT) for the SSRI-paradox falsifier; **mouse-genetic-anchored** (three independently-
generated knockout genotypes) plus a PMC-full-text human/mouse quantitative anchor for the pool-separation
falsifier — consistent with, not upgraded beyond, the task's own framing.

## Repro

```
source .venv-msk/bin/activate
python3 scripts/msk/serotonin_system.py
```

Pure closed-form arithmetic + `scipy.optimize.curve_fit`/`fsolve` over literature-anchored parameters — no
external data dependency, <5s wall time. Writes only `data/serotonin_system/serotonin_system_results.json`
(did not exist before this session — first write, not a clobber). This session's independent citation
verification: raw NCBI eutils esearch/esummary/efetch responses (verbatim) in
`data/raw_fetch/serotonin_*.{json,txt}` + one PMC full-text pull (`data/raw_fetch/serotonin_yano2015_pmc.xml`).
No git operations, no edits to any pre-existing file.
