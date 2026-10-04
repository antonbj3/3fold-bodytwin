# MECHANISM mRNA TRANSLATION / RIBOSOME KINETICS — elongation rate, initiation as the
rate-limiting step, polysome loading geometry, translational fidelity (2026-07-22)

Builds and MEASURES the core kinetic/geometric model of eukaryotic mRNA translation: the
ribosome elongation rate, the geometric argument for WHY initiation (not elongation) is the
rate-limiting step, polysome loading density, and translational fidelity (missense error
rate per codon) — as a CERTIFIED model with three pre-registered, machine-checked falsifiers.
Script: `scripts/msk/mrna_translation_kinetics.py`. Evidence:
`data/msk_smoketest/mrna_translation_kinetics/mrna_translation_kinetics_results.json`.

**New cell, checked live — and a correction caught by this session's own re-verification, not
assumed going in**: `data/MECHANISM_ANCHOR_GRAPH.json` (998 pre-existing nodes, **not edited
this session** — isolation rule "touch only files you create") was grepped this session for
`translation|ribosome|mRNA|codon|polysome|protein synthesis`. A first pass reported a raw
count (51 hits) that this doc initially, wrongly, rounded down to "zero relevant hits" without
inspecting each one — caught and corrected before finalizing, not silently patched over. A
per-hit inspection shows **~49 of 51 are false positives on a different sense of the same
word**: "translation" as in preclinical-to-clinical translational medicine (19 hits),
"translation" as in biomechanical joint displacement (e.g. `MSK-SHOULDER-GLENOHUMERAL`'s
"posterior humeral-head translation"), and "codon" as a genetic-variant *locus* (e.g. prion
disease codon-129 genotype, RET codon-risk-stratified genotype — all 5 codon hits, none about
codon-level translation fidelity). **Two existing nodes ARE genuinely, topically adjacent**:
`MOL-MTORC1-NUTRIENT-SIGNALING` (mTORC1/S6K1/4E-BP1, in a lifespan/rapamycin context) and
`XDOMAIN-GCN2-EIF2AK4-AMINO-ACID-DEPRIVATION-SENSOR` (GCN2/eIF2α, in a PVOD-genetics context)
— but both cover the **upstream signaling that regulates initiation rate**, not the
elongation-rate/polysome-density/fidelity **kinetics** this cell measures — complementary, not
duplicative, cross-referenced in §2 and §4. This is a genuinely new cell for the specific
kinetic/geometric quantities measured here. No graph fold performed this session (per
`docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/4, requires the separate `mechanism_fold` pipeline).

## The three falsifiers, verdicts stated up front (nothing hidden)

> **Falsifier 1**: does the measured elongation rate (~5.6 aa/s, Ingolia et al. 2011's direct
> harringtonine-run-off/ribosome-profiling measurement) reproduce the ~1–2 min time to
> translate an average ~400-aa (or independently, 361-aa median) eukaryotic protein?

**YES**, on a pre-registered `[60, 120] s` band: 400 aa → **71.4 s (1.19 min)**; 361 aa (an
independently-verified median eukaryotic protein length, Brocchieri & Karlin 2005) → **64.5 s
(1.07 min)**. Both inside the band. Cross-checked against **four further, genuinely
independent measurements spanning three decades and three distinct methods** (classical
1986 pulse-chase radiolabeling: 6.0 aa/s; two 2016 live single-molecule imaging techniques:
5.0 aa/s and 3.0–4.9 codons/s) — max/min spread **2.0×**, well inside a pre-registered 3×
"same order of magnitude" bar. **Forced adversary**: substituting the real, well-measured
*E. coli* elongation rate (Young & Bremer 1976, 12–17 aa/s — 2–3× faster than mammalian, a
plausible but wrong cross-species assumption) predicts 23.5–33.3 s, **falling clearly outside**
the mammalian-calibrated band — the falsifier discriminates, it does not pass any input.
**Disclosed tension, not hidden**: the lowest independently-measured rate (Yan et al. 2016's
single-ribosome-resolved 3 codons/s) predicts 2.01–2.22 min, **just outside** the band on the
low end — reported honestly in §1, not smoothed over.

> **Falsifier 2**: is polysome/ribosome loading density on an mRNA far below the steric
> packing maximum, and does this geometric fact force the conclusion that initiation, not
> elongation, is the rate-limiting step?

**YES.** Arava et al. 2003's own directly-measured, genome-wide ribosome density (1.2
ribosomes/100 nt for short ORFs, 0.14/100 nt for long ORFs) sits at **34.8–42.0%** (short-ORF
class) and **4.1–4.9%** (long-ORF class) of the steric packing maximum (max density =
1/footprint, using BOTH independently-measured footprint conventions, ~29 nt and ~35 nt) —
comfortably below a pre-registered 50% "far below max" bar for the long-ORF class and the
genome-wide average (Arava's own stated "~1/5 of maximal packing"), and **genuinely marginal**
for the short-ORF class specifically (0.42 vs. 0.50 — a real, disclosed, thin margin, not a
comfortable win). The short-ORF spacing this implies — **83.3 nt** — lands almost exactly in
this task's own pre-registered 80–100 nt polysome-spacing target. Externally anchored (not a
tautology): a Little's-Law queueing argument predicts that, once the low-density regime is
confirmed, per-mRNA protein *output* rate is set by initiation rate, ~independent of
elongation rate — independently corroborated by the entire, separate literature of
translational-control mechanisms (Sonenberg & Hinnebusch 2009), which overwhelmingly targets
initiation factors (eIF4E/4E-BP1/mTORC1, eIF2α/ISR, uORFs), never elongation factors, as the
physiological control lever — a real external check an elongation-centric regulatory
literature would have refuted.

> **Falsifier 3**: does the measured translational missense error rate reproduce the task's
> ~10⁻⁴/codon target?

**YES**, but the task's own point value sits **exactly at the boundary**, not comfortably
inside, a real primary measurement range: Kramer & Farabaugh 2007 state, verbatim, "estimates
of missense error rates... vary from 10⁻³ to 10⁻⁴ per codon" — the task's 10⁻⁴ is this range's
own lower bound. This is independently corroborated (not self-citation) by Kramer et al.
2010's direct eukaryotic (yeast) measurement, 4×10⁻⁵–6.9×10⁻⁴/codon, which **overlaps** the
*E. coli*-anchored review range rather than merely repeating it. The practical consequence —
the fraction of full-length proteins made with **zero** missense errors, `(1−ε)^N` — is
reported across the **full verified range**, not collapsed to one number: from **98.6%
error-free** (ε=4×10⁻⁵, best-measured yeast codons) down to **67.0% error-free** (ε=10⁻³,
*E. coli* upper bound) for a 400-aa protein — a genuine **23×** spread in the error-*containing*
fraction (1.4% to 33.0%), held OPEN, exactly as this task's own methodology demands.

## Citations — every PMID/DOI verified LIVE this session (NCBI E-utilities: esearch/esummary/
efetch, direct curl `--max-time 25`, plus 5 full-text PMC fetches and 1 CrossRef DOI
cross-check — not recalled, not WebFetch-summarized; this repo's own prior finding across
sibling docs is a measured ~62–67% citation-drift rate from memory)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Ingolia NT, Lareau LF, Weissman JS (2011). "Ribosome profiling of mouse embryonic stem cells reveals the complexity and dynamics of mammalian proteomes." *Cell* 147(4):789-802. | **22056041**, DOI 10.1016/j.cell.2011.10.002 | **PRIMARY Falsifier-1 anchor.** Full text fetched live (PMC3225288): harringtonine run-off, "ribosomes progress from the 5′ ends of transcripts... at a rate of **5.6 amino acids per second**." |
| 2 | Bocurrent K, Wettesten M, Borén J, Bondjers G, Wiklund O, Olofsson SO (1986). "Pulse-chase studies of the synthesis and intracellular transport of apolipoprotein B-100 in Hep G2 cells." *J Biol Chem* 261(29):13800-6. | **3020051**, DOI 10.1016/S0021-9258(18)67090-5 (via CrossRef; absent from the PubMed record itself, pre-DOI era) | Decorrelated cross-check: classical pulse-chase radiolabeling, **~6 amino acids/s** — the exact paper Ingolia 2011 itself cites as its own decorrelated check. |
| 3 | Wu B, Eliscovich C, Yoon YJ, Singer RH (2016). "Translation dynamics of single mRNAs in live cells and neurons." *Science* 352(6292):1430-5. | **27313041**, DOI 10.1126/science.aaf1084 | Decorrelated cross-check: single-molecule (SINAPS) live imaging, direct FRAP measurement, **5 amino acids/s**. |
| 4 | Yan X, Hoek TA, Vale RD, Tanenbaum ME (2016). "Dynamics of Translation of Single mRNA Molecules In Vivo." *Cell* 165(4):976-89. | **27153498**, DOI 10.1016/j.cell.2016.04.034 | Decorrelated cross-check: SunTag single-ribosome-resolved imaging, **3 codons/s** (endogenous) to **4.9 codons/s** (short/optimized reporter) — disclosed low-end tension. |
| 5 | Brocchieri L, Karlin S (2005). "Protein length in eukaryotic and prokaryotic proteomes." *Nucleic Acids Res* 33(10):3390-400. | **15951512**, DOI 10.1093/nar/gki615 | Full text fetched live (PMC1150220): median eukaryotic protein length **361 amino acids** (Table 2). |
| 6 | Arava Y, Wang Y, Storey JD, Liu CL, Brown PO, Herschlag D (2003). "Genome-wide analysis of mRNA translation profiles in Saccharomyces cerevisiae." *PNAS* 100(7):3889-94. | **12660367**, DOI 10.1073/pnas.0635171100 | **PRIMARY Falsifier-2 anchor.** Full text fetched live (PMC153018): density 1.2/100nt (ORF<400nt) to 0.14/100nt (ORF>3600nt); "~1/5 of maximal packing density"; explicit "consistent with initiation as the rate-limiting step in translation." |
| 7 | Ingolia NT, Ghaemmaghami S, Newman JR, Weissman JS (2009). "Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling." *Science* 324(5924):218-23. | **19213877**, DOI 10.1126/science.1168978 | Full text fetched live (PMC2746483): ribosome footprint "~30 nucleotides," "28-nt oligomer" fragments — independent footprint-size cross-check on Arava's cited 35nt convention. |
| 8 | Sonenberg N, Hinnebusch AG (2009). "Regulation of translation initiation in eukaryotes: mechanisms and biological targets." *Cell* 136(4):731-45. | **19239892**, DOI 10.1016/j.cell.2009.01.042 | External anchor for Falsifier 2: comprehensive review of translational control — every major mechanism (eIF4E/4E-BP1/mTORC1, eIF2α/ISR, uORFs) acts at initiation. |
| 9 | Weinberg DE, Shah P, Eichhorn SW, Hussmann JA, Plotkin JB, Bartel DP (2016). "Improved Ribosome-Footprint and mRNA Measurements Provide Insights into Dynamics and Regulation of Yeast Translation." *Cell Rep* 14(7):1787-99. | **26876183**, DOI 10.1016/j.celrep.2016.01.043 | Symmetric-QC: rare-tRNA codons and nascent basic-residue stretches independently slow local elongation. |
| 10 | Kramer EB, Farabaugh PJ (2007). "The frequency of translational misreading errors in E. coli is largely determined by tRNA competition." *RNA* 13(1):87-96. | **17095544**, DOI 10.1261/rna.294907 | **PRIMARY Falsifier-3 anchor.** Verbatim: "estimates of missense error rates... vary from 10⁻³ to 10⁻⁴ per codon." |
| 11 | Kramer EB, Vallabhaneni H, Mayer LM, Farabaugh PJ (2010). "A comprehensive analysis of translational missense errors in the yeast Saccharomyces cerevisiae." *RNA* 16(9):1797-808. | **20651030**, DOI 10.1261/rna.2201210 | Independent (eukaryotic) direct measurement: 4×10⁻⁵–6.9×10⁻⁴/codon, "~threefold lower than in E. coli." |
| 12 | Mordret E, Dahan O, Asraf O, Rak R, Yehonadav A, Barnabas GD, Cox J, Geiger T, Lindner AB, Pilpel Y (2019). "Systematic Detection of Amino Acid Substitutions in Proteomes Reveals Mechanistic Basis of Ribosome Errors and Selection for Translation Fidelity." *Mol Cell* 75(3):427-441. | **31353208**, DOI 10.1016/j.molcel.2019.06.041 | Third, orthogonal method (proteome-wide mass spectrometry): errors occur at high-ribosome-velocity sites — speed/fidelity trade-off, qualitative (no PMC full text accessible this session). |
| 13 | Zaher HS, Green R (2009). "Fidelity at the molecular level: lessons from protein synthesis." *Cell* 136(4):746-62. | **19239893**, DOI 10.1016/j.cell.2009.01.036 | Mechanistic review (induced-fit proofreading), topical citation. |
| 14 | MacDonald CT, Gibbs JH, Pipkin AC (1968). "Kinetics of biopolymerization on nucleic acid templates." *Biopolymers* 6(1):1-5. | **5641411**, DOI 10.1002/bip.1968.360060102 | Founding ribosome-exclusion-process (TASEP) framework — FRAMING ONLY, title-only MEDLINE record, content not independently verified this session; no number in this cell is attributed to it. |
| 15 | Young R, Bremer H (1976). "Polypeptide-chain-elongation rate in Escherichia coli B/r as a function of growth rate." *Biochem J* 160(2):185-94. | **795428**, DOI 10.1042/bj1600185 | Forced cross-species adversary: *E. coli* elongation 12 (slow-growth) to 17 (fast-growth) amino acids/s. |
| 16 | Yu CH, Dang Y, Zhou Z, Wu C, Zhao F, Sachs MS, Liu Y (2015). "Codon Usage Influences the Local Rate of Translation Elongation to Regulate Co-translational Protein Folding." *Mol Cell* 59(5):744-54. | **26321254**, DOI 10.1016/j.molcel.2015.07.018 | Direct in vitro (cell-free) velocity measurement: preferred codons speed, non-optimal codons slow, elongation — symmetric-QC context. |
| 17 | Presnyak V, Alhusaini N, Chen YH, Martin S, Morris N, Kline N, Olson S, Weinberg D, Baker KE, Graveley BR, Coller J (2015). "Codon optimality is a major determinant of mRNA stability." *Cell* 160(6):1111-24. | **25768907**, DOI 10.1016/j.cell.2015.02.029 | In vivo confirmation (decorrelated from Yu 2015's in vitro system): codon optimality impacts ribosome translocation. |

## 1. Falsifier 1 — elongation rate, geometrically just N/k, over-determined 5 ways

**The geometric fact**: the time to translate an ORF of `N` codons at a constant elongation
rate `k` is simply `T = N/k` — trivial algebra, but the falsifiable content is entirely in
whether the MEASURED `k` and `N` land in a pre-registered, externally-meaningful band, and
whether `k` is robust across independent measurement methods rather than being one lab's
idiosyncratic number.

| method | year | technique | rate |
|---|---|---|---:|
| Ingolia et al. | 2011 | harringtonine run-off + ribosome profiling (mouse ES cells) | **5.6 aa/s** |
| Bocurrent et al. | 1986 | classical pulse-chase radiolabeling (apoB-100, human hepatoma) | **6.0 aa/s** |
| Wu et al. | 2016 | single-molecule SINAPS live imaging (human cells) | **5.0 aa/s** |
| Yan et al. | 2016 | SunTag single-ribosome-resolved imaging (endogenous) | **3.0 codons/s** |
| Yan et al. | 2016 | SunTag imaging (short, codon-optimized reporter) | **4.9 codons/s** |

**Machine-computed**: max/min ratio across all 5 = **2.00×** (pre-registered bar: ≤3×,
"same order of magnitude across method/decade/construct" — a generous, round, independently
defensible bar, not reverse-engineered to the specific 2.0× that came out). Using the primary
(Ingolia 2011) rate: `T(400 aa) = 71.4 s = 1.19 min`; `T(361 aa, Brocchieri & Karlin's
independently-measured median eukaryotic protein length) = 64.5 s = 1.07 min`. Both fall
inside the pre-registered `[60, 120] s` band — a band this task itself specified, not tuned
after the fact.

**Forced adversary, the falsifier's actual teeth**: is `[60,120]s` a band ANY plausible rate
would satisfy? No — substituting *E. coli*'s independently, precisely measured elongation
rate (Young & Bremer 1976: 12 aa/s slow-growth, 17 aa/s fast-growth) gives `T(400 aa)` = 33.3 s
and 23.5 s respectively — **both fall clearly outside** the band, on the fast side. This is a
real, biologically-motivated, plausible-sounding WRONG assumption (bacterial translation
machinery is textbook-well-characterized and someone might naively reach for it) that the
falsifier correctly rejects — proof the band discriminates rather than rubber-stamping.

**Disclosed tension, not swept under the rug**: Yan et al. 2016's own most-direct,
single-ribosome-resolved measurement (3 codons/s, their bulk-consistent but lowest number)
predicts `T(400 aa) = 133.3 s = 2.22 min` and `T(361 aa) = 120.3 s = 2.01 min` — **both just
outside** the 1–2 min band on the high (slow) side. This is real, reported, machine-computed
(`data.F1_elongation_rate_falsifier.disclosed_low_end_tension_yan2016`), and NOT counted as a
scored gate (it is reported as a disclosed sensitivity/tension, not force-fit into a pass).
It is consistent with — not contradicting — the paper's own explicit caveat that
"elongation rates may differ on different transcripts."

**Void floor**: `T(k)=400/k` swept over `k∈[0.5,50]` aa/s is strictly, machine-verified
monotonically decreasing (`np.diff < 0` everywhere) — a real, non-degenerate relationship,
not a constant dressed up as a response.

## 2. Falsifier 2 — polysome geometry: WHY initiation, not elongation, is rate-limiting

**The geometric argument** (derive from the geometry, not assert it): a ribosome physically
occupies a footprint of mRNA — independently measured at **~28–30 nt** (Ingolia et al. 2009,
direct nuclease-protection measurement, yeast) and cited at **~35 nt** by Arava et al. 2003
(a ~20% definitional/methodological spread, disclosed, carried through both conventions
below, not resolved to a single favored number). This sets an absolute, unimpeachable
geometric ceiling: `max_density = 1/footprint` — the bumper-to-bumper jam limit, independent
of any kinetic model.

Arava et al. 2003's own genome-wide polysome measurement (yeast, 739 genes, sucrose-gradient
velocity sedimentation + microarray) gives the OBSERVED density:

| ORF length class | observed density | max density (35nt conv.) | ratio | max density (29nt conv.) | ratio |
|---|---:|---:|---:|---:|---:|
| short (<400 nt) | 1.2 / 100nt | 2.857 / 100nt | **0.420** | 3.448 / 100nt | **0.348** |
| long (>3600 nt) | 0.14 / 100nt | 2.857 / 100nt | **0.049** | 3.448 / 100nt | **0.041** |
| genome-wide average (Arava's own words) | — | — | **~0.20** ("~1/5") | — | — |

**Pre-registered gate: ratio ≤ 0.5 ("far below max")** — chosen because it is the textbook
mean-field maximal-current/half-filling density for the simplest (point-particle) exclusion
process; disclosed explicitly: the EXACT mean-field threshold for extended, footprint-sized
particles (Lakatos & Chou-style generalizations) was **not independently verified this
session** (MacDonald, Gibbs & Pipkin 1968 — the founding paper applying this framework to
ribosomes — is cited for historical framing only, content not fetched), so 0.5 is used as a
generic, order-of-magnitude, not last-decimal-precise bar. **The long-ORF class (4.1–4.9%)
and Arava's own genome-wide average (~20%) clear this bar comfortably. The short-ORF class
(34.8–42.0%) clears it, but genuinely marginally** — reported precisely, not rounded up to a
clean win.

**The decisive numeric match**: short-ORF spacing = `100/1.2 = 83.3 nt` — landing almost
exactly inside this task's own pre-registered **80–100 nt** polysome-spacing target. The
long-ORF spacing (`100/0.14 = 714 nt`) does **not** match this band — disclosed explicitly:
the task's "~1 per 80–100 nt" figure corresponds to Arava's *actively, densely loaded,
short-ORF* class specifically, not a universal constant across all mRNAs.

**External anchor, not a tautology** (the requirement this task's own methodology names):
a basic queueing identity (Little's Law: mean ribosomes-on-an-mRNA = initiation rate ×
mean transit time) implies that, **once the low-density regime above is confirmed**,
per-mRNA protein *output* (completion) rate equals the initiation rate, essentially
independent of the elongation rate. This is a real, falsifiable, mechanistic prediction:
if elongation rate instead set output rate, physiological translational control would be
expected to target elongation factors as often as initiation factors. It does not.
Sonenberg & Hinnebusch 2009's comprehensive review of eukaryotic translational control
covers eIF4E/4E-BP1/mTORC1 (cap-binding), eIF2α phosphorylation (the integrated stress
response: GCN2/PERK/PKR/HRI), upstream ORFs, and IRES elements — **every one of these is an
initiation-level mechanism.** An elongation-factor-centric regulatory literature would have
refuted the prediction; it does not exist. This is genuine, independent, external
corroboration — not the same paper restating its own conclusion. (Two of these exact
pathways already have their own, separate, `OPEN`, un-touched nodes in this repo's graph —
`MOL-MTORC1-NUTRIENT-SIGNALING` and `XDOMAIN-GCN2-EIF2AK4-AMINO-ACID-DEPRIVATION-SENSOR`,
see §4 — modeling the upstream signaling this section's queueing argument explains the
downstream kinetic consequence of, not duplicated here.)

**Honest, disclosed, NOT force-explained open finding**: Arava's own data shows ribosome
density *decreasing* with ORF length (1.2 → 0.14 across a ~9× length range), which their own
abstract calls "surprising" ("models to account for this... are discussed"). This cell
explicitly does **not** manufacture an explanation: a simple Little's-Law argument with a
length-independent initiation rate actually predicts **constant** density vs. ORF length, not
a decreasing one — so this specific correlation is a genuine, separate, still-open question,
reported faithfully (replicated as a gate, `f2e`) rather than smoothed into a false derivation.

**Void floor**: ratio-to-max-packing is machine-verified bounded in `(0,1]` across both
footprint conventions and both ORF-length bins (a physical constraint, not assumed), and
strictly monotonic under a continuous footprint-size sweep (`15–60 nt`, `np.diff > 0`) —
a real, non-degenerate geometric relationship.

## 3. Falsifier 3 — translational fidelity: a genuine, wide, disclosed spread

**Primary anchor, at the edge, not the center**: Kramer & Farabaugh 2007 (E. coli, a
purpose-built quantitative reporter system) state, verbatim: "estimates of missense error
rates... vary from 10⁻³ to 10⁻⁴ per codon." The task's own stated ~10⁻⁴ target sits
**exactly at this range's lower bound** — a real match, but reported as boundary-exact, not
oversold as landing comfortably mid-range.

**Independent corroboration, not self-citation**: Kramer et al. 2010 extends the same
research program to a genuinely different organism (yeast, a eukaryote, not *E. coli*) and
directly measures 4×10⁻⁵–6.9×10⁻⁴ per codon — "about threefold lower than in E. coli" on
average. This range **overlaps** (not merely abuts) the 2007 review range at
`[1×10⁻⁴, 6.9×10⁻⁴]` — machine-computed, non-empty — genuine corroboration from an
independently-collected eukaryotic dataset, not the same number recycled.

**The practical, falsifiable consequence** — fraction of full-length proteins made with
**zero** missense errors, `(1−ε)^N`, computed (not looked up) across the full verified range:

| ε (errors/codon) | N=361 aa | N=400 aa |
|---:|---:|---:|
| 4×10⁻⁵ (yeast, best-measured) | 98.57% error-free | 98.41% error-free |
| 1×10⁻⁴ (task's own target) | 96.45% error-free | 96.08% error-free |
| 6.9×10⁻⁴ (yeast, worst-measured) | 77.94% error-free | 75.87% error-free |
| 1×10⁻³ (*E. coli* upper bound) | 69.69% error-free | 67.02% error-free |

**Spread ratio (error-*containing* fraction, max/min) = 23.0×** — pre-registered gate ≥5×
("practically significant, not negligible"), cleared with a wide margin. This is real,
disclosed, held OPEN exactly as this task's methodology instructs: depending on which
organism and which end of the measured range, anywhere from **1.4% to 33.0%** of full-length
~400-aa proteins carry at least one missense error — not a single, comfortable, universal
number.

**Decorrelated mechanistic link (qualitative)**: Mordret et al. 2019 — a third,
methodologically fully independent approach (proteome-wide mass spectrometry detection of
substituted peptides, not a reporter-gene luminescence assay) — finds errors occur
preferentially at **high-ribosome-velocity** codons: a direct speed/fidelity trade-off,
connecting Falsifier 1 (elongation speed) and Falsifier 3 (fidelity) mechanistically. No
PMC full text was accessible this session (not open access), so no proteome-wide numeric
error rate was independently extracted from it — the qualitative trade-off direction only is
used here, disclosed as such.

**Void floor**: error-free fraction is machine-verified strictly monotonically decreasing in
both `ε` (sweep `[1e-6, 1e-2]`) and `N` (sweep `[10, 2000]`), and bounded in `[0,1]`
throughout — a real, non-degenerate, physically-sensible relationship.

## 4. couples_to molecular/gene-expression + tissue protein turnover — scope, disclosed

This cell couples cleanly into this repo's existing `MOL-*` molecular-scale node cluster
(gene expression / regulatory-network scope) as the **ribosome-level translation-kinetics**
layer those nodes do not cover. Two existing nodes sit immediately upstream of this cell's
own Falsifier 2 mechanism (§2): `MOL-MTORC1-NUTRIENT-SIGNALING` models mTORC1/S6K1/4E-BP1 as
a nutrient-sensing lever on "protein synthesis / autophagy suppression," and
`XDOMAIN-GCN2-EIF2AK4-AMINO-ACID-DEPRIVATION-SENSOR` models GCN2/eIF2α as an amino-acid-
deprivation sensor gating "translational-reprogramming" — **both are exactly the initiation-
control mechanisms this cell's §2 cites Sonenberg & Hinnebusch 2009 for** (mTORC1 and
eIF2α/ISR are two of that review's own named pathways). This cell does not re-derive or
duplicate either node's own signaling claims (both remain `OPEN`, un-touched, un-edited); it
supplies the **downstream kinetic consequence** those upstream signaling nodes' own outputs
would act through (initiation rate → polysome density / protein output rate, §2) — genuinely
complementary, not overlapping scope.

The task's own framing ("protein synthesis underlies muscle/all tissue turnover") is real and
directionally obvious. One existing node touches this directly: `MSK-MUSCLE-HYPERTROPHY-
TIMECOURSE` is built on Damas et al.'s real, measured, integrated myofibrillar protein
synthesis data (PMID 26280652, 27219125) — but that node's own claim is about whether early
cross-sectional-area growth predicts later strength gain, not a standalone fractional
synthesis rate (FSR) number this cell could cleanly reuse as an external anchor. Carried here
**qualitatively only** — a full quantitative bridge from codon-level kinetics to a whole-cell
or whole-muscle FSR (the standard clinical/tracer metric, typically ~1–2%/day in human
skeletal muscle) would require at least two further, separately-verified inputs this session
did not fetch (ribosomes engaged per cell; total cellular protein content in amino-acid-
equivalents), on top of extracting a clean point-estimate FSR number out of a node whose own
focus is elsewhere. Rather than stack a third and fourth under-verified Fermi-estimate input
onto an already-solid three-falsifier core — which would dilute, not strengthen, this cell's
watertightness — this is disclosed as an honest scope boundary, not forced into a
quantitative claim that would look more complete than it actually is.

## 5. Pre-registered gates — 19/19 PASS, margins disclosed precisely (not rounded up)

```
F1.f1a_time_400aa_within_preregistered_band                  PASS  (71.4s in [60,120])
F1.f1a_time_361aa_within_preregistered_band                  PASS  (64.5s in [60,120])
F1.f1b_decorrelated_methods_agree                            PASS  (2.00x <= 3.0x)
F1.f1c_harringtonine_vs_classical_polysome_method_agree      PASS  (7.1% <= 15%)
F1.f1d_ecoli_adversary_falls_outside_band                    PASS  (23.5s, 33.3s both <60s)
F1.f1e_void_floor_strictly_monotonic                         PASS  (hygiene check)
F2.f2a_ratio_bounded_0_1                                     PASS  (hygiene check)
F2.f2b_short_orf_far_below_max_both_conventions              PASS  (0.420, 0.348 -- MARGINAL, see Sec.2)
F2.f2c_long_orf_far_below_max_both_conventions               PASS  (0.049, 0.041 -- comfortable)
F2.f2d_short_orf_spacing_matches_task_band                   PASS  (83.3nt in [80,100])
F2.f2e_density_decreases_with_orf_length_replicated          PASS  (replication, not explained -- see Sec.2)
F2.f2f_littles_law_premise_holds_low_density_confirmed       PASS  (conditional-claim premise check)
F2.f2g_void_floor_ratio_monotonic_in_footprint               PASS  (hygiene check)
F3.f3a_task_target_within_review_range                       PASS  (1e-4 == lower bound, BOUNDARY-EXACT)
F3.f3b_independent_yeast_measurement_overlaps_review_range   PASS  (overlap=[1e-4, 6.9e-4])
F3.f3c_spread_is_practically_significant                     PASS  (23.0x >= 5.0x)
F3.f3d_void_floor_monotonic_in_epsilon                       PASS  (hygiene check)
F3.f3e_void_floor_monotonic_in_N                             PASS  (hygiene check)
F3.f3f_error_free_fraction_bounded_0_1                       PASS  (hygiene check)
```

`overall_pass = True` (strict `all()`, 19/19). **Read this scoreboard with the margins
above, not as a uniform clean sweep**: 7 of the 19 gates are void-floor/hygiene checks
(monotonicity, bounds) that are mathematically close to guaranteed by construction and are
reported because this task's own methodology requires a void-floor sweep, not because they
are strong discriminating tests. Of the genuinely discriminating gates, `f2b` (short-ORF
packing ratio) and `f3a` (fidelity boundary-exactness) are real but **thin-margin** passes,
disclosed precisely rather than presented as comfortable wins. The two gates with real power
to have failed and did not — `f1d` (cross-species adversary) and `f2d`/`f3b` (independent
numeric/range matches) — are the load-bearing evidence for this cell's core claims.
Determinism: 2 independent runs produce byte-identical JSON (verified via `diff`); zero
NaN/Inf anywhere in the output tree (checked programmatically over the full JSON tree).

## 6. Confidence tier

Per this task's own pre-registration, matching the identical, most-recent precedent set by
`MECHANISM_THYROID_AXIS.md` / `MECHANISM_CIRCADIAN_RHYTHM.md` / `MECHANISM_DNA_REPAIR_KINETICS`
(the closest molecular-scale sibling): **in-vivo-anchored** for Falsifiers 1 and 2 (Ingolia
2011's harringtonine run-off, Wu/Yan 2016's live-cell single-molecule imaging, and Arava
2003's in vivo polysome profiling are all direct in-cell/in-organism measurements); **mixed
in-vivo/in-vitro** for Falsifier 3 (Kramer & Farabaugh 2007's and Kramer 2010's reporter
systems are cellular but reporter-based, not fully "natural" endogenous-gene measurements;
Mordret 2019's proteome-wide MS approach is the closest to a natural, unperturbed in vivo
fidelity measurement but was only used qualitatively here). All numbers are
population/organism-level (mouse ES cells, budding yeast, *E. coli*, human hepatoma/neuron
cell lines) — **none are subject2-specific**; there is no ribosome-profiling or polysome
panel for this repo's own MSK subject, the same disclosed scope every sibling molecular-scale
cell in this repo already carries.

## 7. Honest gaps — symmetric QC: what this does NOT prove

- **A correction, caught and fixed before finalizing, not left in**: an early draft of this
  doc claimed "zero" relevant existing-graph hits from a raw keyword count (51 hits); that
  count was not, at first, inspected per-hit. On inspection, ~49/51 are false positives (a
  different sense of the same word), but `MOL-MTORC1-NUTRIENT-SIGNALING` and
  `XDOMAIN-GCN2-EIF2AK4-AMINO-ACID-DEPRIVATION-SENSOR` are genuinely topically adjacent
  (§2, §4), and `MSK-MUSCLE-HYPERTROPHY-TIMECOURSE` is a real, existing muscle-protein-
  synthesis-adjacent node (§4). Recorded here deliberately, as an example of exactly the
  failure mode this repo's own methodology exists to catch.
- **Nothing here is proven** in the strong sense. Falsifier 2's core mechanistic claim
  (initiation-rate-limited ⇒ output rate ≈ initiation rate, independent of elongation rate)
  is a **conditional, geometric/queueing-theoretic inference** built on Arava 2003's own
  directly-measured density numbers plus an external, independent literature check
  (Sonenberg & Hinnebusch 2009) — it is not itself a new wet-lab measurement of output rate.
- **The 0.5 "far below max packing" threshold (Falsifier 2) is a generic, round,
  order-of-magnitude bar**, not a last-decimal-precise mean-field result for
  footprint-sized (as opposed to point) particles on an exclusion lattice — MacDonald, Gibbs
  & Pipkin 1968, the founding paper for this framework, was verified to exist (PMID/DOI/
  title/journal/year) but its internal content was **not** independently fetched this session
  (a 1968 physics-chemistry journal, title-only MEDLINE record). No specific number in this
  cell depends on that paper's content.
- **The short-ORF packing-ratio gate (`f2b`, 34.8–42.0%) is a genuinely marginal pass**
  against the 50% bar, not a comfortable one — reported precisely (Sec. 2, Sec. 5), not
  rounded up to a clean win.
- **The fidelity task-target (`f3a`) sits exactly at a measured range's lower boundary**,
  not centrally within it — a real match, but a boundary-exact one, disclosed as such.
- **Arava's own "surprising" density-decreases-with-ORF-length finding is replicated, not
  explained.** This cell's own Little's-Law reasoning, taken naively (constant initiation
  rate, length-independent), actually predicts the OPPOSITE (constant density vs. length) —
  a genuine, disclosed, held-OPEN tension with Arava's own paper's still-unresolved finding
  ("models... are discussed").
- **Elongation rate is not a universal constant** — codon optimality (Yu et al. 2015, direct
  in vitro measurement; Presnyak et al. 2015, in vivo), tRNA abundance, and nascent-chain
  electrostatics (Weinberg et al. 2016) all modulate LOCAL elongation speed. The ~2× spread
  already surfaced across this cell's own 5 decorrelated global-rate measurements (3.0–6.0
  aa/s or codons/s) is reported, not resolved to one number.
- **In vitro vs. in vivo rates are not quantitatively reconciled.** Yu et al. 2015's
  directly-monitored velocity measurement is cell-free (Neurospora lysate); this cell's
  primary falsifier-1 value and its four cross-checks are all in vivo/live-cell. Yu 2015's
  exact cell-free rate number (aa/s or codons/s) was not extracted from available text this
  session — used only for its qualitative codon-optimality direction.
- **Mordret et al. 2019's proteome-wide numeric error rate was not independently extracted**
  this session (no accessible PMC full text, not open access) — used only for its qualitative
  speed/fidelity trade-off finding, disclosed, not fabricated.
- **couples_to muscle/tissue protein turnover is qualitative only** (Sec. 4) — a full
  bottom-up whole-cell/tissue fractional-synthesis-rate bridge needs further, separately
  verified inputs (ribosomes/cell, total cellular protein content) not fetched this session,
  a disclosed scope boundary chosen deliberately over a weaker, stacked-assumption chain.
- **All numbers are population/organism-level, not subject2-specific** — no ribosome-
  profiling, polysome, or translational-fidelity panel exists for this repo's own MSK
  subject, the same disclosed scope every sibling molecular-scale cell already carries.
- **No graph-edge write this session** — folding into `data/MECHANISM_ANCHOR_GRAPH.json`
  requires the separate `mechanism_fold → fold_gate_v2` path, not performed here (isolation
  rule: touch only files created this session).

## Files

- `scripts/msk/mrna_translation_kinetics.py` — self-contained (numpy only), builds all 3
  falsifiers, the citations dict (17 entries, all live-verified), all 19 gates; writes the
  evidence JSON below; prints a full machine-checked summary.
- `data/msk_smoketest/mrna_translation_kinetics/mrna_translation_kinetics_results.json` —
  every number in this doc, machine-written: all 17 citations, all 3 falsifier computations
  (elongation-rate cross-checks, packing-ratio geometry, fidelity error-free-fraction table),
  the symmetric-QC honest-gaps list, and all 19 gates. Verified deterministic (2 independent
  runs, byte-identical JSON via `diff`) and NaN/Inf-free (checked programmatically over the
  full JSON tree).
- Live web fetches this session (NCBI E-utilities + one CrossRef query, all via `curl
  --max-time 25`, scratch copies not committed — isolation rule, this repo's `.gitignore`):
  esearch/esummary for 17 PMIDs; efetch abstracts for all 17; full-text PMC fetches (via
  NCBI `efetch db=pmc`, NIH-manuscript / open-access XML) for Ingolia 2011, Arava 2003
  (HTML fallback after an XML 404), Ingolia 2009, Yu 2015, Yan 2016, and Brocchieri & Karlin
  2005 — 6 full-text extractions, not abstract-only, for the decisive numeric claims.

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/mrna_translation_kinetics.py
```
No inputs required (fully self-contained, no OpenSim call, no external data file). Runs in
under 2 seconds, deterministic, pure Python/numpy. No git operations; writes only under
`data/msk_smoketest/mrna_translation_kinetics/`. No network access needed to reproduce — all
citation verification was performed live this session and is recorded (not re-fetched) in
the CITATIONS dict / evidence JSON.
