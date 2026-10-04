# MECHANISM TASTE TRANSDUCTION — five modalities, two transduction classes, one shared bottleneck (2026-07-22)

Script: `scripts/msk/taste_transduction.py`. Raw results:
`data/taste_transduction/taste_transduction_results.json`. Evidence (citations):
`docs/MECHANISM_TASTE_TRANSDUCTION_evidence.json`. 21 sources (20 PMIDs + 1 population dataset),
every PMID live-verified this session via NCBI eutils (`curl --max-time 25`, esearch→esummary→efetch);
the NHANES dataset independently re-computed from its raw `.XPT` file, not narration.

**Provenance note:** fresh build this session. Repo-wide grep for T1R2/T1R3/T2R/TRPM5/CALHM1/OTOP1/gustducin
found zero pre-existing taste-transduction doc or script. The only pre-existing taste-adjacent artifact is
the anchor-graph node `CHEMOSENSORY-OLFACTION-GUSTATION`, whose own text states *"taste is the cluster's
clearest zero-coverage gap"* and names the exact residual needed: *"acquire taste-bud/receptor molecular
data ... then taste half is a 2-sided cert."* This doc is that molecular-receptor leg. Not a duplicate.

## 0. Why this layer + couples_to

- **`data/MECHANISM_ANCHOR_GRAPH.json` node `CHEMOSENSORY-OLFACTION-GUSTATION`** — this doc supplies the
  taste-receptor molecular leg that node explicitly flags as missing. That node's own `verify` field already
  names `nhanes-taste-csx` as the dataset that "partially unblocks" taste; §4 below independently re-derives
  numbers from that exact dataset (not copied from the node, not narration).
- **`docs/MECHANISM_OLFACTORY_TRANSDUCTION.md`** — sibling chemosensory-transduction mechanism doc, same
  GPCR→Gα→2nd-messenger→channel *shape*, but a genuine point of **divergence**: taste's sweet/umami/bitter
  arm uses gustducin (a Gq-family-like Gα) → PLCβ2 → IP3/Ca²⁺, while olfaction uses Golf → adenylyl cyclase →
  cAMP. Two GPCR-coupled senses, two different second-messenger arms — not a duplicate.
- **`docs/MECHANISM_NERVE_CONDUCTION.md`** — grep-confirmed zero gustatory-afferent (CN VII/IX/X; chorda
  tympani/glossopharyngeal/vagal) content. That doc's generic AP-generation/conduction-velocity machinery is
  the architectural downstream consumer of this doc's output (the CALHM1→ATP→P2X2/P2X3 generator potential
  in the afferent terminal), not yet built for this fiber type — structural coupling, disclosed gap, not built here.
- **GPCR signaling generally** (docs with PLCβ/IP3/Ca²⁺ content, e.g. `MECHANISM_INSULIN_PI3K_AKT.md`,
  `MECHANISM_MAPK_ERK_SIGNALING.md`) — same cascade *topology* reused with a different receptor family and
  physiological role; architectural family resemblance, not content overlap.
- **The sensory/CNS layer + flavor integration** — CT/glossopharyngeal/vagal afferents → nucleus tractus
  solitarius → thalamus (VPM parvocellular) → gustatory cortex (insula/operculum); flavor = this doc's output
  integrated with `MECHANISM_OLFACTORY_TRANSDUCTION.md`'s retronasal-smell output (both certs exist, their
  integration is not attempted here).

## 1. Citations — 20 PMIDs + 1 population dataset (n=3708), all live-verified this session

| # | Citation | PMID | Role |
|---|---|---|---|
| 1 | Adler E et al (2000). *Cell* 100(6):693-702. | 10761934 | T2R gene family discovery |
| 2 | Chandrashekar J et al (2000). *Cell* 100(6):703-11. | 10761935 | T2Rs function as bitter receptors |
| 3 | Nelson G et al (2001). *Cell* 106(3):381-90. | 11509186 | T1R2+T1R3 = sweet receptor |
| 4 | Nelson G et al (2002). *Nature* 416(6877):199-202. | 11894099 | T1R1+T1R3 = umami receptor |
| 5 | Wong GT, Gannon KS, Margolskee RF (1996). *Nature* 381(6585):796-800. | 8657284 | Gustducin-KO: reduced (not abolished) bitter+sweet |
| 6 | Perez CA et al (2002). *Nat Neurosci* 5(11):1169-76. | 12368808 | TRPM5 molecular identification in TRCs |
| 7 | Zhang Y et al (2003). *Cell* 112(3):293-301. | 12581520 | **TRPM5/PLCβ2-KO: DECISIVE double-dissociation** |
| 8 | Damak S et al (2003). *Science* 301(5634):850-3. | 12869700 | **T1R3-KO: diminished not abolished — corrects the premise** |
| 9 | Taruno A et al (2013). *Nature* 495(7440):223-6. | 23467090 | **CALHM1-KO, independent lab, replicates the dissociation** |
| 10 | Ma Z et al (2018). *Neuron* 98(3):547-561. | 29681531 | CALHM3 = obligate heteromer partner (not redundancy) |
| 11 | Meyerhof W et al (2010). *Chem Senses* 35(2):157-70. | 20022913 | 25 hTAS2Rs, combinatorial receptive ranges |
| 12 | Kim UK et al (2003). *Science* 299(5610):1221-5. | 12595690 | **TAS2R38/PTC — decorrelated human-genetics anchor** |
| 13 | Chandrashekar J et al (2010). *Nature* 464(7286):297-301. | 20107438 | **ENaCα-KO: DECISIVE, complete loss of salt taste** |
| 14 | Oka Y et al (2013). *Nature* 494(7438):472-5. | 23407495 | **High salt co-opts bitter+sour cells — corrects the premise** |
| 15 | Tu YH et al (2018). *Science* 359(6379):1047-50. | 29371428 | OTOP1/2/3 = proton-selective channels (discovery) |
| 16 | Teng B et al (2019). *Curr Biol* 29(21):3647-3656. | 31543453 | **OTOP1-KO: DECISIVE, sour-specific attenuation** |
| 17 | Mueller KL et al (2005). *Nature* 434(7030):225-9. | 15759003 | Labeled-line evidence (cell identity sets valence) |
| 18 | Wu A et al (2015). *Nat Commun* 6:8171. | 26373451 | Across-fiber evidence (tuning breadth ∝ concentration) |
| 19 | Mojet J, Christ-Hazelhof E, Heidema J (2001). *Chem Senses* 26(7):845-60. | 11555480 | Psychophysical thresholds, 10 compounds, n=42 |
| 20 | Finger TE et al (2005). *Science* 310(5753):1495-9. | 16322458 | ATP→P2X2/P2X3 afferent step |
| 21 | NHANES 2013-14 CSX_H (Taste & Smell Exam) | n=3708 | **Population behavioral anchor, re-derived from raw data this session** |

## 2. The model

**Sweet + umami (T1R family, class-C GPCR):** T1R2+T1R3 heterodimer = sweet receptor (PMID 11509186);
T1R1+T1R3 heterodimer = umami receptor, broadly tuned to L-amino acids (PMID 11894099). Both couple through
gustducin (a Gα subunit closely related to transducin, PMID 8657284) → PLCβ2 → IP3 → Ca²⁺ release → TRPM5
(a Ca²⁺-activated monovalent cation channel, PMID 12368808) → depolarization → CALHM1 voltage-gated ATP-release
channel (PMID 23467090, obligate CALHM3 heteromer partner per PMID 29681531) → ATP → P2X2/P2X3 on the afferent
terminal (PMID 16322458).

**Bitter (T2R family, ~25 human genes, class-A-like GPCR):** T2Rs (PMID 10761934, 10761935) are exclusively
expressed in gustducin+ cells and converge on the **SAME** downstream cascade (gustducin → PLCβ2 → TRPM5 →
CALHM1) as sweet/umami. A single TRC expresses a large T2R repertoire (PMID 10761934); Meyerhof et al 2010
(PMID 20022913) functionally mapped all 25 against 104 bitter compounds and found combinatorial, unequal
coverage (§3, Gate G4) — the molecular basis for one uniform "bitter" percept from structurally unrelated toxins.

**Salt:** two mechanistically DISTINCT arms, not one graded system. (a) **Appetitive, low-concentration:**
ENaC-expressing TRCs, amiloride-sensitive, genetically dedicated (PMID 20107438). (b) **Aversive,
high-concentration:** NOT a third salt-specific receptor — Oka et al 2013 (PMID 23407495) show high salt
**cross-activates the pre-existing sour- and bitter-sensing cells**; genetically silencing those two
pathways abolishes high-salt aversion while leaving low-salt attraction untouched. This is an **auditable
correction** to a "third amiloride-insensitive salt sensor" framing — the mechanism is cross-wiring, not a
new molecule (see §5).

**Sour:** OTOP1, a proton-selective ion channel (PMID 29371428) in Type III TRCs, conducts H⁺ directly into
the cytosol, acidifying it and driving action-potential firing (PMID 31543453) — no GPCR, no gustducin, no
TRPM5 involvement. Ionotropic, not metabotropic.

**Geometric framing (a convergence bottleneck, not a metaphor):** the sweet/umami/bitter transduction
cascade is a many-to-one converging network. Upstream, MULTIPLE parallel receptor/Gα paths feed the SAME
downstream node: umami has T1R1+T1R3 **and** taste-mGluR4 (Damak 2003's own discussion); the Gα stage has
gustducin **and** transducin (Wong 1996's own discussion). Downstream, TRPM5 and CALHM1 are reported as
**obligate, non-redundant** nodes — no second channel is named that could carry the signal if either is lost.
A node with ≥2 independent parallel inputs loses only a fraction of its signal when ONE input is cut; a node
with exactly one path (no redundancy) loses all of it. This one structural fact — **redundancy upstream,
bottleneck downstream** — is what the data actually show (§3, Gate G2), not an assumption: T1R3-KO and
gustducin-KO (both upstream, both b≥2) give *diminished* phenotypes; TRPM5/PLCβ2-KO and CALHM1-KO (both
downstream, both b=1) give *abolished/severely-impaired* phenotypes. Four-for-four concordant — a real,
falsifiable structural account of *why* some lesions are partial and others total, not decorative language.

## 3. THE falsifier — machine-checked gates (5/5 gated PASS; `taste_transduction.py`)

**Adversary forced:** a "single universal/undifferentiated taste receptor" — no molecular double-dissociation
should be possible if this were true; every lesion should hit every modality roughly equally, and
identification of a stimulus's quality should be near chance.

- **G0 — NHANES CSX_H population behavioral anchor (n=3708, machine-recomputed from the raw `.XPT` this
  session, not narration).** Whole-mouth 1 mM quinine identified as **"Bitter" in 82.66%** of 3114
  non-missing responses, cross-confused as "Salty" in only **1.09%**. Whole-mouth 1 M NaCl identified as
  **"Salty" in 97.23%**, cross-confused as "Bitter" in only **1.06%**. Pre-registered thresholds (≥70% correct,
  ≤10% cross) both cleared by a wide margin. **This is the strongest single anchor in this doc**: human,
  population-scale, self-report psychophysics — decorrelated from every mouse-genetic citation by species,
  method, and institution. Paired within-subject NaCl dose-response (0.32 M → 1 M, a 3.125× concentration
  step): mean intensity (gLMS 0–100) rose 33.08 → 52.06 (1.574×), and **84.44% of n=3214 paired subjects**
  rated the higher concentration as more intense — a real, monotonic, large-N concentration→intensity
  dose-response for one modality, directly measured. **PASS.**
- **G1 — dissociation block-contrast.** A 7-row × 6-column matrix (every cell cited to a verbatim primary
  quote, `evidence.json`) coding impairment on a 0 (normal) → 3 (abolished) scale, honest `None` for
  untested cells. For the 4 rows with BOTH sides measured in their own source paper — TRPM5/PLCβ2-KO,
  CALHM1-KO, gustducin-KO — the margin between the row's own affected block and everything else measured
  is **≥1.0 on the 0–3 scale for all 4** (pre-registered threshold, could have failed). The other 3 rows
  (T1R3, OTOP1, ENaCα) report only their own side in their own abstract — honest gap, not filled in. **PASS
  (4/4 testable).**
- **G2 — funnel-redundancy concordance.** The branching-factor model in §2, checked against measured
  category: 4/4 concordant. **Disclosed as a structural-consistency check, not a blind prediction** — b for
  gustducin/T1R3 was read from the same papers being checked, so a pass is partly expected given a correct
  reading, not fully out-of-sample. **PASS, honestly tiered.**
- **G3 — TAS2R38/PTC decorrelated anchor.** 3 coding SNPs → 5 haplotypes (≤ 2³=8 combinatorial ceiling,
  trivial sanity check) "completely explain the bimodal distribution of PTC taste sensitivity ... 55 to 85%
  of the variance" (Kim et al 2003, verbatim). Independent species (human), independent method (population
  linkage genetics, not transgenics), independent lab (NIH/Drayna) — converges with the mouse-KO
  architecture on the same conclusion: bitter perception is gated by a specific, identifiable gene, not an
  undifferentiated system. **PASS.**
- **G4 — bitter combinatorial overrepresentation.** Meyerhof et al 2010's own numbers: 3 of 25 hTAS2Rs
  jointly cover ≈50% of 104 tested bitter compounds, vs a naive uniform-coverage null of 3/25=12% —
  a **4.17× overrepresentation** (≥2× pre-registered threshold), falsifying "all ~25 receptors contribute
  equally" using the source's own reported numbers. **PASS** (disclosed as an order-of-magnitude comparator,
  not a rigorously derived null distribution).
- **G5 — psychophysics spread (descriptive only, not gated).** Mojet et al 2001's age-related threshold
  fold-changes range from 1.32× (aspartame) to 5.70× (IMP) — a 4.32× spread, sitting alongside (not
  overriding) the source's own more conservative "predominantly generic" conclusion. Only 2 of 10 tested
  compounds' fold-changes are in the abstract; the rest are an honest gap.
- **G6 — labeled-line vs across-fiber-pattern coding: HELD OPEN, not adjudicated (§6).**

## 4. Two auditable corrections to the pre-registered premise

Symmetric QC requires reporting where the primary data does NOT match a pre-registered claim, exactly as
carefully as where it does.

1. **"T1R3-KO must abolish sweet+umami" is FALSIFIED AS STATED.** Damak et al 2003's own data: *"diminished
   but not abolished ... T1r3-independent sweet- and umami-responsive receptors and/or pathways exist."*
   Full abolition is measured only for the artificial-sweetener subset (mice show no preference for
   sucralose/saccharin-type compounds specifically, since those have no non-T1R route). This does **not**
   rescue a "single universal receptor" adversary — the residual response is explained by a SECOND, named,
   modality-specific mechanism (taste-mGluR4 for umami), which is *more* modality-specific molecular
   structure, not less. The corrected, measured claim: T1R3 is *one of at least two* modality-specific
   umami/sweet mechanisms, and removing it produces a strong but partial (not complete) loss.
2. **"An amiloride-insensitive aversive high-salt path" is not a third salt-specific sensor.** Oka et al
   2013's own data: high salt "recruits" and cross-activates the **existing bitter and sour cells**;
   silencing those two (already-characterized) cell populations abolishes high-salt aversion while sparing
   low-salt attraction. There is no separate high-salt-specific transduction molecule — the mechanism is
   cross-wiring of concentration-dependent stimulus strength into pre-existing aversive channels, not a new
   receptor class.

## 5. Concentration→intensity dose-response

Directly measured this session (§3, G0): NaCl whole-mouth intensity rises monotonically with concentration in
84.4% of n=3214 paired human subjects (NHANES CSX_H), 0.32 M→1 M giving a 1.57× intensity increase for a
3.125× concentration increase (compressive, consistent with a saturating/sub-linear psychophysical function,
not tested against a specific power-law exponent this session — an honest scope limit, not a claim). Bitter
coding is combinatorial-and-unequal, not a single dose-response curve: Meyerhof et al 2010's receptive-range
data show individual hTAS2Rs vary from narrowly to extremely broadly tuned (§3, G4). Cross-modality relative
sensitivity (Mojet et al 2001, §3 G5) is reported only at the qualitative/fold-change level this session;
absolute mM thresholds per modality were not extracted from full text — an honest gap, not fabricated from
recall (per this session's own instruction that recall drifts ~62%).

## 6. HELD OPEN — labeled-line vs across-fiber-pattern coding (not adjudicated)

Two real, non-contradicted anchors on different aspects of the question, deliberately not resolved:

- **Labeled-line-supporting:** Mueller et al 2005 (PMID 15759003) — swapping a bitter receptor's ligand
  specificity while keeping expression restricted to gustducin+/bitter-fated cells redirects behavioral
  aversion to the NEW ligand. Cell-TYPE identity, not receptor identity, sets the valence — the strongest
  evidence for a labeled line at the CELL level.
- **Complicating pure labeled-line at the FIBER level:** Wu et al 2015 (PMID 26373451) — individual
  gustatory afferent fibers are narrowly tuned near threshold but become progressively more broadly tuned as
  stimulus concentration rises. Tuning breadth is itself concentration-dependent, not a fixed either/or
  property of a fiber.

These are not measured here to contradict each other: cell-type identity can set the labeled valence at the
periphery while afferent population breadth still varies with concentration — a plausible reconciliation,
**not adjudicated, and not claimed as resolved.**

## 7. Overall result

**5 of 5 gated checks PASS** (G0–G4; G5 descriptive-only, G6 explicitly held open per task instruction).
**21 sources, every PMID live-verified this session** via NCBI eutils; the NHANES CSX_H population dataset
(n=3708) independently re-computed from its raw `.XPT` (pandas.read_sas), cross-checked against — not
blindly trusted from — a WebFetch summary of the CDC codebook, and found to match exactly. **Two auditable
corrections** to the task's own pre-registered premise were found and reported precisely rather than
silently smoothed over (§4) — a rejection of part of the prompt's framing, backed by the primary sources
themselves. The decisive genetic double/triple-dissociation (TRPM5/PLCβ2-KO and, independently, CALHM1-KO in
a different lab, abolish/severely-impair sweet+umami+bitter while sparing sour+salt; OTOP1-KO abolishes sour
specifically; ENaCα-KO abolishes salt-appetitive specifically) falls the "single universal receptor"
adversary across a genuinely diverse instance space (5 distinct genes, 5 distinct labs, 2 species, 2
methodologies — transgenic knockout and human population genetics/behavior). The labeled-line vs
across-fiber-pattern coding debate is reported, not settled, per instruction. Confidence tier:
**in-vivo mouse-genetic-knockout-anchored** (the double-dissociation core) + **molecular/heterologous-expression-anchored**
(receptor identification, combinatorial receptive ranges) + **human-population-genetics-anchored** (TAS2R38)
+ **human-population-behavioral-anchored** (NHANES CSX_H, directly re-derived, not narration).

## Repro

```
python3 scripts/msk/taste_transduction.py
```

Pure stdlib arithmetic over a citation-anchored matrix (no numpy dependency, <1s wall time) for G1–G6.
G0's NHANES numbers are hardcoded from a live fetch + independent pandas re-derivation this session
(URL, fetch method, and exact counts disclosed in `docs/MECHANISM_TASTE_TRANSDUCTION_evidence.json` — the
script does not re-fetch over the network on every run, matching this repo's no-external-runtime-dependency
convention for these cells). Writes only `data/taste_transduction/taste_transduction_results.json` (did not
exist before this session — first write). No git operations, no edits to any pre-existing file, no writes to
`data/MECHANISM_ANCHOR_GRAPH.json` (graph integration is the coordinator's job via the fold path, not this
subagent's).
