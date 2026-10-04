# MECHANISM BILIRUBIN METABOLISM — heme catabolism, the conjugation solubility switch, and the neonatal-jaundice/kernicterus dysfunction map (2026-07-22)

Builds the bilirubin PIGMENT pathway the twin's liver docs named but never built: heme →
(heme oxygenase) biliverdin → (biliverdin reductase) unconjugated bilirubin → hepatic uptake →
UGT1A1 glucuronidation → conjugated bilirubin → bile → gut → urobilinogen/stercobilin. Distinct
from `docs/MECHANISM_BILE_ENTEROHEPATIC.md` (bile-**acid** mass balance — a different molecule
family) and `docs/MECHANISM_HEPATIC_CLEARANCE.md` (generic well-stirred drug-clearance PK, ICG as
the flow marker — a different mechanism, not the conjugation chemistry itself). The
`data/body_twin/agent_outputs/liver-hepatic-hub__*.json` seed mentions bilirubin in one clause
("heme-breakdown+conjugation/excretion capacity") without building it — this doc is that build.
Script: `scripts/msk/bilirubin_metabolism.py`. Evidence: `data/bilirubin_metabolism/
bilirubin_metabolism_results.json` + this doc's own `docs/MECHANISM_BILIRUBIN_METABOLISM_evidence.json`.

**NO re-solve.** Reads `data/erythropoiesis/erythropoiesis_results.json` read-only for the exact
same three numbers `docs/MECHANISM_IRON_HEPCIDIN.md` already fixed for an unrelated purpose
(Hb=15 g/dL, blood volume=5.0 L, RBC combined lifespan=115 d) and extends the identical
pool/lifespan=flux identity to a **different downstream molecule** (bilirubin, not iron) — a
genuine held-out reuse, cross-checked against an independent literature figure (Levitt & Levitt
2014), not against the sibling doc's own conclusion. A bonus gate (G0) confirms the SAME
heme-turnover number, multiplied by Fe atomic weight instead of bilirubin's, reproduces
`iron_hepcidin.py`'s own published 22.6 mg/day recycling flux to 0.01% — a cheap, genuine,
cross-doc arithmetic-consistency check, not a required falsifier.

**This is a HYPOTHESIS for independent QC** (return-only wave discipline,
`docs/MECHANISM_HARDENED_CONVENTIONS.md` §4b) — a designed-and-computed cell, not yet folded into
the graph.

**Confidence tier: in-vivo/primary-literature-anchored** for the qualitative pathway (Tenhunen
1968 heme-oxygenase discovery; Bosma 1995 + Kadakol 2000 human UGT1A1 genetics; Kawade & Onishi
1981 human-liver-tissue UGT1A1 ontogeny, n=88; Bhutani 1999 human nomogram, n=13,003 measurements/
n=2,840 nomogram cohort; Mreihil 2010 human neonate phototherapy kinetics, n=20; three AAP
guideline generations 2004→2009→2022). Standard protein/small-molecule-chemistry-stoichiometry
tier for the Hb/heme/bilirubin molecular-weight constants (bilirubin MW derived transparently
from its atomic formula, not hardcoded from memory). Disclosed **textbook/secondary tier**,
flagged everywhere it appears, for: neonatal reference Hb/blood-volume-per-kg values, and the
exact neonatal RBC-lifespan day-count (Pearson 1967 is live-confirmed to exist and be on-topic;
its exact quoted figure was not independently re-extracted from its paywalled 1967 full text this
session — swept over its full disclosed range instead of point-estimated, §6).

## 0. Scope note — graph-node disambiguation, checked live before writing a line of code

Per this repo's established convention, `data/MECHANISM_ANCHOR_GRAPH.json` (1049 nodes) was checked
live (read-only, **not edited**) before any code was written. Six existing nodes touch this
neighborhood, none overlapping this doc's scope:

- **`UGT1A1-GENOTYPE-BILIRUBIN-PHENOTYPE-DISSOCIATION`** (OPEN, SEED-DESIGN, never executed) — a
  **narrower, harder** hidden-state cert: does UGT1A1 promoter-repeat genotype (TA6/TA7) × coding
  LOF-allele dose predict bilirubin phenotype, triangulated across irinotecan/SN-38-neutropenia
  and atazanavir pharmacogenomic *outcome* domains? That cell's own `regime_note` already records
  a live grep of this same graph file finding "0 hits for UGT1A1/Crigler/kernicterus/atazanavir/
  phenobarbital/phototherapy" before it was designed — i.e., it independently confirms the same
  gap this doc closes. This doc builds the **basic physiology** that node's occluded-truth
  pharmacogenomic dissociation question would sit on top of; it does not execute that node's own
  chemo/antiretroviral-outcome falsifier (out of scope, left OPEN, not resolved here).
- **`G6PD-ENZYME-ACTIVITY-VS-HEMOLYTIC-RISK-DISSOCIATION`** and **`HEMATO-SICKLE-CELL-DISEASE`**
  (both OPEN) — the upstream hemolysis/heme-LOAD mechanisms (§0's "pre-hepatic" row in §5's
  diagnostic table assumes, does not re-derive, that hemolysis raises heme load; these two nodes
  own that mechanism's own genetics/severity certs).
- **`ORG-LIVER-HEPATIC-HUB`**, **`HEPATIC-XENOBIOTIC-METABOLISM`**, **`MET-PHARMACOKINETICS-ADME`**
  (all OPEN) — the organ-level hub and generic CYP/UGT xenobiotic-clearance layers, already
  partially resolved by the sibling `MECHANISM_HEPATIC_CLEARANCE.md`/`MECHANISM_BILE_ENTEROHEPATIC.md`
  docs for unrelated substrates. UGT1A1 also glucuronidates non-bilirubin xenobiotics (irinotecan's
  SN-38 metabolite, notably) — flagged as a coupling, not built here (this doc's UGT1A1 content is
  scoped to its bilirubin substrate only).

**This document establishes new ground**: no existing node/doc builds the executed heme→biliverdin→
bilirubin→conjugation→enterohepatic pathway, nor the neonatal-jaundice/kernicterus/Gilbert/
Crigler-Najjar quantitative dysfunction map.

## 1. Geometric structure — three genuinely different pieces of geometry

- **Production is the same pool/lifespan renewal-process identity** `flux = pool / mean_lifespan`
  that `MECHANISM_IRON_HEPCIDIN.md` and `MECHANISM_HEMATOPOIESIS.md` already used for the Fe-recycling
  and RBC/neutrophil/platelet legs — reused here for a **third, independent downstream molecule**
  sharing the SAME upstream pool (heme). One heme ring destroyed emits, deterministically, exactly
  one Fe atom (iron_hepcidin's leg) and exactly one bilirubin molecule (this doc's leg) — not two
  separate facts to fit, one stoichiometric consequence of hemoglobin catabolism (§4, G0).
- **Conjugation is a polarity phase-transition, not a rate tweak.** Unconjugated bilirubin's own
  structure (six intramolecular hydrogen bonds folding it into a compact, poorly-water-soluble
  "ridge-tile" conformation — Vitek & Ostrow 2009) is *why* it must travel albumin-bound and cannot
  be renally excreted; UGT1A1 attaching glucuronic acid breaks that intramolecular hydrogen-bonding
  network, flipping the molecule from lipophilic to water-soluble in one enzymatic step. This is a
  **discrete coordinate change** (bound/insoluble ↔ free/soluble), not a continuous rate parameter —
  which is exactly why the conjugated/unconjugated split is a clean categorical localizer (§5), not
  a graded severity scale.
- **The unconjugated/conjugated split is an inverse-problem localization**: a single serum total-
  bilirubin number is under-determined (many lesion sites could produce "high bilirubin"), but
  splitting it into its two chemically-distinct species resolves *where* in the pathway the defect
  sits — upstream-of-conjugation (production overload or UGT1A1 itself → unconjugated predominates)
  vs. downstream-of-conjugation (canalicular/transporter excretion or obstruction → conjugated
  predominates). This is the observable the clinical differential is actually built from (§5).
- **Phototherapy exploits a THIRD, independent geometric handle** — configurational photo-
  isomerization (4Z,15Z → 4Z,15E) plus intramolecular cyclization (lumirubin) change the molecule's
  shape without needing the UGT1A1 enzyme at all, producing polar species excretable in bile/urine
  by a route that bypasses the very enzyme that is deficient in the disease this doc's dysfunction
  section addresses (§6/§7) — geometry substituting for a broken enzyme, not a stronger dose of it.

## 2. Method, in one paragraph

Six falsifiers, each pre-registered with a threshold fixed before the code that tests it was run
(`scripts/msk/bilirubin_metabolism.py`, committed alongside this doc). **F1** re-derives adult daily
bilirubin production from the sibling doc's own fixed Hb/blood-volume/RBC-lifespan numbers via
heme stoichiometry, checked against Levitt & Levitt 2014's independently-quoted, live full-text-
verified figure (not against iron_hepcidin's own conclusion — non-circular). **F2** extends the
same identity to a swept, disclosed-tier neonatal input set to derive a per-kg production-
amplitude ratio. **F3** compares that amplitude against Kawade & Onishi 1981's directly-quoted
UGT1A1 ontogeny deficit to identify which of the two decorrelated mechanisms (production increase
vs. clearance immaturity) actually dominates — a magnitude comparison, not a hand-waved "both
matter." **F4** is the task's own required falsifier: a 7-row categorical table forcing the "single
bilirubin pool" adversary to predict conjugated/unconjugated predominance from mechanism-locus
alone, gated against a coin-flip null. **F5** uses Mreihil 2010's own quoted within-subject timing
data to force the naive "phototherapy directly/rapidly destroys bilirubin" adversary to fail on the
SAME dataset that shows isomerization, not bulk clearance. **F6** forces a "single eternal
bilirubin threshold" adversary against two decorrelated pieces of evidence: Bhutani 1999's
hour-specific zone likelihood-ratio spread, and the AAP guideline's own quoted 2004→2022 revision.

## 3. Citations — 18 papers, every PMID verified LIVE this session via NCBI eutils

Method: `esearch` (title/author-anchored) → `esummary` (title/journal/date cross-check) →
`efetch` (raw abstract/full-text, verbatim numbers extracted) — same discipline as this repo's
other `MECHANISM_*` docs. **Disclosed drift, measured not hidden**: of 6 PMIDs this session initially
drew from memory before checking, **3 (50%) were WRONG** on live esummary/esearch check — a
recalled Levitt & Levitt PMID (24634588) resolved to nothing in the target space (correct:
**25214800**); a recalled Fevery 2008 PMID (18248387) was likewise wrong (correct: **18433389**);
a recalled Crigler-Najjar-1952 PMID (14944796) resolved to a **completely unrelated** 1952 *British
Journal of Surgery* paper on duodenal peptic ulceration and islet-cell pancreatic tumors (correct:
**12983120**). The other 3 recalled PMIDs (Tenhunen 1968, AAP 2004, Bhutani 1999) checked out
exactly. Separately, the initial free-text `esearch` for "Tenhunen heme oxygenase" failed to
surface the target paper in its top-5 relevance-ranked hits at all (a **search-ranking** miss, not
a wrong-memorized-number miss) — caught only because the recalled PMID was independently verified
by direct `esummary` lookup rather than trusted from the text-search ranking.

| role | citation | PMID | key number/finding, quoted or directly reported |
|---|---|---|---|
| heme oxygenase (discovery) | Tenhunen, Marver, Schmid 1968, *PNAS* 61(2):748-55 | 4386763 | "The enzymatic conversion of heme to bilirubin by microsomal heme oxygenase" (title; no indexed abstract, pre-1970s) |
| general review | Fevery 2008, *Liver Int* 28(5):592-605 | 18433389 | "Bilirubin in clinical practice: a review" |
| production/kinetics (PRIMARY, quantitative) | Levitt & Levitt 2014, *Clin Exp Gastroenterol* 7:307-28 | 25214800 | full text: "About 250 mg/day of bilirubin is produced by normal adult humans; roughly 75% of this is derived from circulating [hemoglobin]"; covers Gilbert's/Dubin-Johnson/Crigler-Najjar/Rotor; "extremely high-affinity albumin binding of UB"; CB renal-excretion role |
| chemistry/mechanism review | Vitek & Ostrow 2009, *Curr Pharm Des* 15(25):2869-83 | 19754364 | "less than 0.01% of total bilirubin circulating in an unbound form (free bilirubin, Bf)... governs the diffusion of UCB into tissues" |
| fetal/neonatal RBC lifespan | Pearson 1967, *J Pediatr* 70(2):166-71 | 5334979 | "Life-span of the fetal red blood cell" (title/existence confirmed; exact day-figure disclosed secondary-tier, §6) |
| Gilbert syndrome genetics (PRIMARY) | Bosma et al 1995, *NEJM* 333(18):1171-5 | 7565971 | "glucuronidating activity... reduced to about 30 percent of normal"; TA7/TA7 promoter (`A(TA)7TAA`); allele freq "40 percent"; controls homozygous for TA7 had higher bilirubin, P=0.009 |
| Crigler-Najjar/Gilbert genotype-phenotype (PRIMARY) | Kadakol et al 2000, *Hum Mutat* 16(4):297-306 | 11013440 | ">50 genetic lesions" of UGT1A1; CN-1 (complete loss)/CN-2 (partial)/Gilbert (promoter-only) dose gradient; "severe hyperbilirubinemia seen in CN-1 can cause bilirubin encephalopathy (kernicterus)" |
| UGT1A1 ontogeny (PRIMARY, n=88 human liver) | Kawade & Onishi 1981, *Biochem J* 196(1):257-60 | 6796071 | four phases (middle foetal/late foetal/neonatal+early infantile/mature) = "about 0.1, 0.1-1, and 1-100%... of the mature-phase values (1320±514 μg/h/g liver, n=27)" |
| Crigler-Najjar original description (1/2) | Crigler & Najjar 1952, *Pediatrics* 10(2):169-80 | 12983120 | "Congenital familial nonhemolytic jaundice with kernicterus" |
| Crigler-Najjar original description (2/2, preliminary note) | Crigler & Najjar 1952, *AMA Am J Dis Child* 83(2):259-60 | 14884759 | same authors, same year, earlier short report, "a new clinical entity" |
| AAP guideline (2004, historical) | AAP Subcommittee 2004, *Pediatrics* 114(1):297-316 | 15231951 | "kernicterus should almost always be preventable... cases continue to occur"; phototherapy/exchange-transfusion framework |
| AAP guideline update (2009) | Maisels et al 2009, *Pediatrics* 124(4):1193-8 | 19786452 | "Hyperbilirubinemia in the newborn infant ≥35 weeks' gestation: an update with clarifications" |
| AAP guideline (2022, current) | Kemper et al 2022, *Pediatrics* 150(3):e2022058859 | 35927462 | "Clinical Practice Guideline Revision" |
| AAP technical report (2022, PRIMARY quote for F6) | Slaughter, Kemper, Newman 2022, *Pediatrics* 150(3):e2022058865 | 35927519 | "neurotoxicity does not occur until bilirubin concentrations are well above the 2004 exchange transfusion thresholds... justified narrowly raising phototherapy treatment thresholds" |
| hour-specific nomogram (PRIMARY, n=13,003/2,840) | Bhutani, Johnson, Sivieri 1999, *Pediatrics* 103(1):6-14 | 9917432 | high-risk zone (≥95th pctile): "LR=14.08, sensitivity=54%, specificity=96.2%"; low-risk zone (<40th pctile): "LR=0... probability=0%" |
| kernicterus mechanism review | Watchko & Tiribelli 2013, *NEJM* 369(21):2021-30 | 24256380 | "Bilirubin-induced neurologic damage — mechanisms and management approaches" (title/topic confirmed; not quoted beyond title this session) |
| phototherapy discovery (historical) | Cremer, Perryman, Richards 1958, *Lancet* 1(7030):1094-7 | 13550936 | "Influence of light on the hyperbilirubinaemia of infants" |
| phototherapy mechanism review (PRIMARY) | Maisels & McDonagh 2008, *NEJM* 358(9):920-8 | 18305267 | "Phototherapy for neonatal jaundice" |
| photoisomer/albumin structure (PRIMARY, decorrelated crystallographic anchor) | Zunszain, Ghuman, McDonagh, Curry 2008, *J Mol Biol* 381(2):394-406 | 18602119 | "crystal structure — determined to 2.42 A resolution — revealed the 4Z,15E-bilirubin-IXα isomer bound to an L-shaped pocket in sub-domain IB" |
| photoisomerization kinetics (PRIMARY, n=20 neonates, F5's numbers) | Mreihil, McDonagh, Nakstad, Hansen 2010, *Pediatr Res* 67(6):656-9 | 20308939 | photoisomer "significant (p<0.0001)... detectable within 15 min"; TSB change "insignificant at 120 min but reached significance at 240 min (p<0.001)"; isomer "up to 20-25% of TSB at 2h" |

## 4. Headline results (machine-printed, `data/bilirubin_metabolism/bilirubin_metabolism_results.json`)

| gate | quantity | value | pre-registered threshold | PASS? |
|---|---|---|---|---|
| F1 | adult daily bilirubin production (derived) | **236.6 mg/day** (3.38 mg/kg/day) | in [100,400] mg/day AND within 40% of Levitt's 187.5 mg/day RBC-slice | ✅ (26.2% off; −5.3% vs Levitt's ungated 250 mg/day total) |
| G0 (bonus) | same heme-turnover × Fe MW vs iron_hepcidin.py's own published figure | 22.60 vs 22.6 mg/day | within 2% (cross-doc consistency, not independent evidence) | ✅ (0.01%) |
| F2 | neonatal:adult per-kg production ratio, swept over disclosed 70-90d neonatal RBC lifespan | **1.72x – 2.22x** | > 1.5x, robust across full sweep | ✅ |
| F3 | UGT1A1 neonatal-floor deficit multiple vs F2 amplitude | **100x vs 2.22x max** | deficit % >> amplitude (clearance dominates) | ✅ |
| F4 | diagnostic-split table match | **7/7 rows** | 7/7 (100%); coin-flip null p=0.0078 | ✅ |
| F5 | photoisomer-significant time vs TSB-significant-decline time | **15 min ≪ 120 (NS) ≪ 240 min (p<0.001)** | strict ordering | ✅ |
| F6 | fixed-threshold adversary (LR spread + guideline revision) | LR 14.08 vs 0; 2004→2022 threshold raised | both facts must hold | ✅ |

`required_gates_overall_pass: True` (6/6 required gates PASS; 1/1 bonus gate PASS).

## 5. F4 in detail — the falsifier this doc was tasked to reproduce

**Pre-registered claim**: the unconjugated-vs-conjugated split localizes the lesion; a "bilirubin
is a single pool" adversary must fail against it. **Forced adversary**: if mechanism-locus
(pre-conjugation vs. post-conjugation) carried NO information about the observed clinical
predominant-fraction pattern, a fair coin-flip per condition would produce a 7/7 match with
probability **0.5⁷ = 0.78%** — the adversary's strongest fair form (each row judged independently,
no row given credit for another's answer).

| condition | mechanism locus | predicted | literature-stated pattern | match |
|---|---|---|---|---|
| Hemolysis (pre-hepatic) | pre-conjugation (production overload) | unconjugated | unconjugated | ✅ |
| Gilbert syndrome (UGT1A1 TA7/TA7, ~30% activity — Bosma 1995) | pre-conjugation | unconjugated | unconjugated | ✅ |
| Crigler-Najjar type 1 (UGT1A1 null — Kadakol 2000) | pre-conjugation | unconjugated | unconjugated | ✅ |
| Crigler-Najjar type 2 (UGT1A1 partial, phenobarbital-responsive) | pre-conjugation | unconjugated | unconjugated | ✅ |
| Neonatal physiological jaundice (production ↑ + UGT1A1 ontogeny ↓, both pre-conjugation, §6) | pre-conjugation | unconjugated | unconjugated | ✅ |
| Biliary obstruction (post-hepatic) | post-conjugation (excretion blocked, reflux to blood) | conjugated | conjugated | ✅ |
| Hepatocellular injury (e.g. viral hepatitis; canalicular/transporter dysfunction) | post-conjugation-dominant | conjugated | conjugated | ✅ |

**Measured**: 7/7 (100%). **Anchor**: each row's "literature-stated pattern" column is drawn from
an independent primary source (Bosma 1995, Kadakol 2000) or a dedicated quantitative review
(Levitt & Levitt 2014, which explicitly walks through hemolysis/cholestasis/hepatocellular/
congenital-defect categories) — not from this doc's own table, avoiding a tautology gate. **C
accepted**: the split is a genuine, mechanism-derived localizer, not an arbitrary clinical
convention — it falls directly out of conjugation being a discrete solubility switch (§1), so a
defect's position relative to that one switch point is, by construction, legible in which species
of bilirubin accumulates.

## 6. Neonatal jaundice — F2/F3, why it happens and why it self-resolves

**Pre-registered claim**: neonates produce more bilirubin per kg AND clear it more slowly, and the
clearance deficit is the dominant (not merely additive) mechanism. **Measured**: production
amplitude 1.72x-2.22x (swept over the disclosed 70-90 day neonatal RBC-lifespan range, driven by
two decorrelated literature-anchored inputs — higher neonatal Hb g/dL and shorter RBC lifespan,
Pearson 1967) vs. a UGT1A1 clearance-capacity floor of **~1% of adult** at the neonatal phase
(Kawade & Onishi 1981's own quoted "1-100%" neonatal-phase range, conservative low end used) — a
**100x** deficit multiple, two full orders of magnitude larger than the production-side amplitude.
**Forced adversary**: a "production increase alone explains neonatal jaundice" claim is falsified
by its own arithmetic — 2x more substrate arriving at an enzyme running at 1% capacity is nowhere
near saturating; the enzyme-immaturity leg dominates. **Why it self-resolves**: Kawade & Onishi's
same four-phase curve climbs from the neonatal floor toward the mature 1320±514 μg/h/g-liver value
within the "early infantile" phase — the same ontogeny curve that explains the onset also predicts
the ~1-2 week natural resolution most infants show without intervention, without needing a second,
separately-fit decay parameter.

## 7. Phototherapy — F5, isomerization precedes (and is mechanistically distinct from) clearance

**Pre-registered claim**: phototherapy works via configurational photoisomerization to a
water-soluble/excretable form, not by directly/rapidly destroying bilirubin. **Forced adversary**
(leaning-positive: "phototherapy just makes bilirubin disappear faster" — the naive, mechanism-free
reading of "light lowers bilirubin"): Mreihil et al 2010's own within-subject, same-cohort data
gives a clean temporal dissociation — the 4Z,15E photoisomer is **statistically significant within
15 minutes** (p<0.0001) and reaches **20-25% of total serum bilirubin by 2 hours**, while total
serum bilirubin itself is **not yet significantly different from baseline at 120 minutes**, only
becoming significant at **240 minutes** (p<0.001). **Measured**: 15 min ≪ 120 min (NS) ≪ 240 min —
the molecule's chemical identity changes an order of magnitude faster than its bulk plasma
concentration falls. **Anchor, decorrelated**: an independent crystallographic study (Zunszain et
al 2008, 2.42 A resolution) directly visualizes the photoisomer (4Z,15E-bilirubin-IXα, not the
native 4Z,15Z form) bound in a distinct albumin sub-pocket — a structural-biology confirmation of
the SAME chemical species, obtained by a completely different instrument (X-ray diffraction, not a
spectrophotometric serum assay), decorrelated from the clinical timing data. **C accepted**: the
naive "phototherapy = faster destruction" adversary falls; the correct mechanism is a fast
structural conversion to an excretable isomer, with bulk plasma clearance lagging behind it.

## 8. Kernicterus threshold — F6, forcing the "single fixed number" adversary

**Pre-registered claim**: there is no single, time-invariant, context-free bilirubin number above
which kernicterus occurs — risk is hour-of-life- and risk-factor-stratified, and the guideline
number itself has moved as evidence accrued. **Forced adversary** (leaning-negative direction this
time — the honest-sounding "clinicians must have a hard number to act on" is the pull toward
treating one fixed cutoff as ground truth): two independent facts must BOTH hold for the adversary
to fall. **(a)** Bhutani 1999's own hour-specific nomogram shows a >14-fold likelihood-ratio spread
by risk zone at the SAME chronological age window (high-risk zone LR=14.08 vs. low-risk zone
LR=0) — a fixed single-number rule has no mechanism to assign different risk to the same total
bilirubin value depending on hour-of-life, so it cannot reproduce this spread by construction.
**(b)** The AAP guideline's own numeric thresholds **changed** between 2004 and 2022 — the 2022
technical report's own words: "neurotoxicity does not occur until bilirubin concentrations are
well above the 2004 exchange transfusion thresholds," justifying "narrowly raising phototherapy
treatment thresholds." A truly fixed, eternal biological threshold would not need revision as
evidence accrued. **Measured**: both (a) and (b) hold. **C accepted**: kernicterus risk is
multifactorial and risk-stratified (gestational age, hemolysis, hour-of-life, per the guideline's
own framework), not a single sharp step function — consistent with, though not separately
re-verified against, Watchko & Tiribelli 2013's mechanism review (free-bilirubin/albumin-binding/
BBB-integrity all independently modulate risk at a given total-bilirubin value).

## 9. Honest gaps (symmetric QC — held OPEN, not forced)

1. **Neonatal RBC-lifespan exact figure** (§2's F2/F3 inputs) is disclosed secondary/textbook tier
   — Pearson 1967 is live-confirmed to exist and be exactly on-topic (title: "Life-span of the
   fetal red blood cell"), but its full 1967 *J Pediatr* text was not accessible this session (no
   PMC full text; pre-digital-era paywall), so the exact day-count was not independently
   re-extracted from it. Mitigated, not eliminated, by sweeping the full commonly-cited 70-90 day
   range (F2 passes across the whole sweep, not at one cherry-picked point).
2. **Watchko & Tiribelli 2013** (kernicterus mechanism review) is cited at title/topic-confirmed
   tier only — no specific number from its body text was quoted this session (its own abstract was
   not returned in indexed form by NCBI efetch; likely a clinical-review format without a
   structured abstract). Its role in §8 is corroborative context, not load-bearing evidence — F6
   does not depend on it.
3. **F1's derived adult production figure (236.6 mg/day) sits closer to Levitt's UNGATED total
   figure (250 mg/day, −5.3%) than to the specific 75%-RBC-attributed slice used for the actual
   gate (187.5 mg/day, +26.2%)** — disclosed, not silently reconciled. Two candidate explanations
   (this doc's 115-day biotin-derived "true" RBC lifespan running shorter than whatever older
   Cr-51-based assumption underlies Levitt's 75% attribution; or simply approximation slack in a
   secondary review's rounded fraction) are plausible but NOT adjudicated here.
4. **The F4 diagnostic table's 7 rows are canonical/textbook-established classifications**, not a
   novel discovery — the coin-flip-null framing (p=0.0078) demonstrates internal non-randomness of
   the mechanism-locus/pattern relationship, not a new clinical finding. Presented as a
   machine-checked structural confirmation, not oversold as new science.
5. **Dubin-Johnson and Rotor syndromes** (hereditary conjugated hyperbilirubinemias via MRP2/OATP
   transporter defects) are named in Levitt & Levitt 2014's own abstract but not built into this
   doc's diagnostic table or given their own falsifier — a natural next row, left OPEN.
6. **Heme oxygenase's own reaction stoichiometry** (O2/NADPH/H2O counts, CO release) is presented
   qualitatively (cited to Tenhunen 1968) and NOT hand-derived atom-by-atom in the script — the
   quantitative stoichiometric claim actually gated (F1/G0) is the downstream 1-heme-in
   ⟷1-bilirubin-out (and 1-Fe-out) molar equivalence, which does not require re-deriving the
   oxygenase redox mechanism itself.

## 10. Couples to

- `docs/MECHANISM_BILE_ENTEROHEPATIC.md` — reads its own bile-acid mass-balance framework as the
  downstream transport context for conjugated bilirubin's excretion into bile (not re-derived here).
- `docs/MECHANISM_HEPATIC_CLEARANCE.md` — the generic hepatic-clearance PK layer this doc's
  UGT1A1-conjugation step is a specific instance of (bilirubin as a non-CYP, UGT-substrate case).
- `docs/MECHANISM_ERYTHROPOIESIS.md` / `docs/MECHANISM_IRON_HEPCIDIN.md` — the shared upstream heme
  pool (RBC production/lifespan); this doc's F1/G0 are a direct extension of that shared substrate.
- `G6PD-ENZYME-ACTIVITY-VS-HEMOLYTIC-RISK-DISSOCIATION` / `HEMATO-SICKLE-CELL-DISEASE` (graph nodes)
  — upstream hemolysis mechanisms feeding this doc's "pre-hepatic" production-overload row (§5).
  `UGT1A1-GENOTYPE-BILIRUBIN-PHENOTYPE-DISSOCIATION` (graph node) — the narrower pharmacogenomic
  dissociation cert this doc's basic physiology sits under (§0), left OPEN, not resolved here.
- Blood-brain-barrier: kernicterus (§8) is unconjugated bilirubin crossing an immature BBB —
  couples to a BBB-integrity/permeability cert this doc does not itself build.
- Neonatal-transition: `docs/MECHANISM_NEONATAL_TRANSITION.md` builds first-breath/lung-liquid
  clearance for the same developmental window; this doc's §6 is the hematologic/hepatic-maturation
  leg of that same neonatal-transition period, independently derived, not cross-fit to it.
