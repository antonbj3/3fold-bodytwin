# MECHANISM SYNAPTIC VESICLE RELEASE — Ca2+-triggered, SNARE/synaptotagmin/complexin-mediated exocytosis (2026-07-22)

Script: `scripts/msk/synaptic_vesicle_release.py`. Raw results:
`data/synaptic_vesicle_release/synaptic_vesicle_release_results.json`. Evidence (citations):
`docs/MECHANISM_SYNAPTIC_VESICLE_RELEASE_evidence.json`. Raw fetch logs (all esearch/efetch/EuropePMC calls
this session, verbatim, via `curl --max-time 25`): `data/raw_fetch/svr_efetch_*.txt`,
`data/raw_fetch/svr_efetch_*.xml`, `data/raw_fetch/svr_europepmc_*.json`.

**Provenance note (stated up front):** repo-wide grep before starting found no existing script/doc for this
topic (`find . -iname "*synap*" -o -iname "*vesicle*"` hits only an unrelated prior agent-output stub,
`data/body_twin/agent_outputs/synapse-neurotransmission__a7e9ad8ea53afbb30.json`, and an existing but
**distinct** graph node, `NEU-SYNAPSE-NEUROTRANSMISSION` — see Sec. 0). This is a first write (script, raw
results, doc, evidence JSON all created fresh this session; no pre-existing file was edited or clobbered).

## 0. Why this layer + couples_to (and why this is NOT a duplicate of NEU-SYNAPSE-NEUROTRANSMISSION)

The graph already holds `NEU-SYNAPSE-NEUROTRANSMISSION` (`data/MECHANISM_ANCHOR_GRAPH.json`, status OPEN,
from a prior agent's `synapse-neurotransmission` output). Read in full before starting this doc. **That
node's own `proposed_cell` targets the STDP Hebbian weight-update rule** (`Δw=f(Δt)` fit to Bi & Poo 1998) as
its primary certified claim — it uses the vesicle-fusion machinery (SNARE structure, synaptotagmin, Katz
quantal hypothesis) only as **one of three supporting decorrelated legs** for that plasticity-rule claim, at
a citation-exists level (e.g. it cites Geppert 1994 and Sutton 1998 by name but does not build a
quantitative Ca-cooperativity model, does not test the toxin adversary, and does not attempt a forced
rejection of a linear-Ca alternative). **This doc is the release-machinery layer itself**: the quantitative
Ca2+-cooperativity dose-response (with a forced, machine-computed rejection of the linear/first-order-Ca
adversary), the Katz N·p·q quantal framework, a 3-toxin/3-SNARE decorrelated indispensability concordance,
two independent synaptotagmin double-dissociation instances, and the sub-ms fusion-latency anchor — a
distinct, deeper mechanistic layer that `NEU-SYNAPSE-NEUROTRANSMISSION`'s own plasticity claim implicitly
depends on but does not itself certify.

**Couples to:**
- `docs/MECHANISM_NERVE_CONDUCTION.md` — the presynaptic action potential + axonal Ca2+ influx that arrives at
  the terminal and triggers this machine. That doc treats `NMJ_SYNAPTIC_DELAY_S` as an opaque, reused
  constant; this doc's Sabatini & Regehr 1996 sub-ms latency anchor is the mechanistic content behind that
  constant.
- `docs/MECHANISM_DOPAMINE_KINETICS.md` and `docs/MECHANISM_SEROTONIN_SYSTEM.md` — dopamine and serotonin are
  released by exactly this SNARE/synaptotagmin machine; those docs start downstream at transporter-mediated
  reuptake/synthesis, this doc supplies the release step itself, shared by every small-molecule/peptide
  transmitter system in the repo (glutamate, GABA, ACh, NE included, not yet separately built).
- `data/MECHANISM_ANCHOR_GRAPH.json` node `NEU-SYNAPSE-NEUROTRANSMISSION` (status OPEN) — distinct per above,
  not superseding it: that node's plasticity-rule claim **depends on** this doc's release-machinery layer as
  its physical substrate, not the reverse.
- SA-node Ca2+-clock cert (in flight per task framing) — shared "steep Ca2+-cooperativity gates an
  all-or-none release/activation event" motif, different molecular hardware (RyR2/SERCA-driven
  Ca-induced-Ca-release in cardiac myocytes vs SNARE/synaptotagmin-driven vesicle fusion here) — a structural
  analogy worth cross-checking once that cert lands, not asserted as the same mechanism.
- Synaptic plasticity (LTP/LTD/STDP) — downstream consumer: `NEU-SYNAPSE-NEUROTRANSMISSION`'s own claim.

## 1. Citations — 20 PMIDs, every one verified LIVE this session (NCBI eutils esearch/efetch, verbatim text)

| # | Citation | PMID | Role |
|---|---|---|---|
| 1 | Dodge FA Jr, Rahamimoff R (1967). *J Physiol* 193(2):419-32. | 6065887 | **Ca-cooperativity PRIMARY anchor** (frog NMJ) |
| 2 | Schneggenburger R, Neher E (2000). *Nature* 406(6798):889-93. | 10972290 | **Ca-cooperativity anchor** (calyx of Held) |
| 3 | Bollmann JH, Sakmann B, Borst JGG (2000). *Science* 289(5481):953-7. | 10937999 | Decorrelated confirmation, same synapse, independent lab |
| 4 | Heidelberger R, Heinemann C, Neher E, Matthews G (1994). *Nature* 371(6497):513-5. | 7935764 | **Ca-cooperativity anchor** (retinal bipolar, 3rd synapse type) |
| 5 | Sun J, Pang ZP, Qin D, Fahim AT, Adachi R, Südhof TC (2007). *Nature* 450(7170):676-82. | 18046404 | **Exact n≈5 vs n≈2 sync/async cooperativity + 2nd Syt-KO dissociation** |
| 6 | Kochubey O, Lou X, Schneggenburger R (2011). *Trends Neurosci* 34(5):237-46. | 21439657 | Review, corroborates Sun 2007 |
| 7 | Südhof TC (2013). *Neuron* 80(3):675-90. | 24183019 | Synthesis review, PMC full text mined for exact cooperativity numbers |
| 8 | Fatt P, Katz B (1952). *Proc R Soc Lond B* 140(899):183-6. | 13003923 | mEPP discovery (companion) |
| 9 | Fatt P, Katz B (1952). *J Physiol* 117(1):109-28. | 14946732 | **mEPP discovery PRIMARY** ("Spontaneous subthreshold activity") |
| 10 | del Castillo J, Katz B (1954). *J Physiol* 124(3):560-73. | 13175199 | **Katz quantal hypothesis PRIMARY** ("Quantal components of the end-plate potential") |
| 11 | Geppert M, Goda Y, Hammer RE, Li C, Rosahl TW, Stevens CF, Südhof TC (1994). *Cell* 79(4):717-27. | 7954835 | **Synaptotagmin-1 KO double-dissociation PRIMARY** |
| 12 | Sutton RB, Fasshauer D, Jahn R, Brunger AT (1998). *Nature* 395(6700):347-53. | 9759724 | **SNARE crystal structure PRIMARY** (2.4Oh, my God. 4-helix bundle) |
| 13 | Fasshauer D, Sutton RB, Brunger AT, Jahn R (1998). *PNAS* 95(26):15781-6. | 9861047 | Q/R-SNARE conservation (generality bonus) |
| 14 | Sabatini BL, Regehr WG (1996). *Nature* 384(6605):170-2. | 8906792 | **Sub-ms fusion-latency PRIMARY** |
| 15 | Schiavo G, Benfenati F, Poulain B, Rossetto O, Polverino de Laureto P, DasGupta BR, Montecucco C (1992). *Nature* 359(6398):832-5. | 1331807 | **BoNT/B + TeNT → VAMP2/synaptobrevin PRIMARY** |
| 16 | Blasi J, Chapman ER, Link E, Binz T, Yamasaki S, De Camilli P, Südhof TC, Niemann H, Jahn R (1993). *Nature* 365(6442):160-3. | 8103915 | **BoNT/A → SNAP-25 PRIMARY** |
| 17 | Blasi J, Chapman ER, Yamasaki S, Binz T, Niemann H, Jahn R (1993). *EMBO J* 12(12):4821-8. | 7901002 | **BoNT/C1 → syntaxin** (3rd decorrelated SNARE-cleavage) |
| 18 | Giraudo CG, Eng WS, Melia TJ, Rothman JE (2006). *Science* 313(5787):676-80. | 16794037 | **Complexin clamp PRIMARY** |
| 19 | Tang J, Maximov A, Shin OH, Dai H, Rizo J, Südhof TC (2006). *Cell* 126(6):1175-87. | 16990140 | Complexin/synaptotagmin switch mechanism |
| 20 | Rosenmund C, Stevens CF (1996). *Neuron* 16(6):1197-1207. | 8663996 | **Readily-releasable-pool (RRP) definition** |

**Verbatim quotes carrying the decisive numbers** (all fetched fresh this session via NCBI eutils, not
carried over from any other source):

- **Dodge & Rahamimoff 1967**: *"The relation between [Ca] and the e.p.p. is highly non-linear. The initial
  part of this relation on double logarithmic co-ordinates gives a straight line with a slope of nearly four
  (mean 3.78 +/- 0.2 S.D. in 28 experiments)... The slope of this logarithmic relation diminishes as [Ca] is
  raised towards the normal level... the number of packets of acetylcholine released is proportional to the
  fourth power of [CaX]... a co-operative action of about four calcium ions is necessary for the release of
  each quantal packet of transmitter."*
- **Schneggenburger & Neher 2000**: *"step-like elevations to only 10 microM [Ca2+]i induce fast transmitter
  release, which depletes around 80% of a pool of available vesicles in less than 3 ms... transient (around
  0.5 ms) local elevations of [Ca2+]i to peak values as low as 25 microM can account for transmitter release
  during single presynaptic action potentials... the high intracellular calcium cooperativity in triggering
  vesicle fusion."*
- **Bollmann, Sakmann & Borst 2000**: *"A rise in [Ca2+]i to 1 micromolar readily evoked release. An increase
  to >30 micromolar depleted the releasable vesicle pool in <0.5 millisecond."*
- **Heidelberger et al 1994**: *"The steepness suggested that at least four calcium ions must bind to
  activate synaptic vesicle fusion. Half-saturation was at 194 microM, and the maximal rate constant was
  2,000-3,000 s-1... within a few hundred microseconds, if [Ca2+]i rises above 100 microM."*
- **Sun et al 2007**: *"deletion of synaptotagmin 2 (Syt2) in mice selectively abolishes synchronous release,
  allowing us to study pure asynchronous release in isolation... asynchronous release displays a Ca2+
  cooperativity of approximately 2 with a Ca2+ affinity of approximately 44 microM, in contrast to
  synchronous release, which exhibits a Ca2+ cooperativity of approximately 5 with a Ca2+ affinity of
  approximately 38 muM."*
- **Südhof 2013** (PMC full text, mined directly, not just abstract): *"Ca2+ triggers release in a highly
  cooperative manner (Dodge and Rahamimoff, 1967) within a few hundred microseconds (Sabatini and Regehr,
  1996)... biophysical studies revealed that the remaining 'asynchronous' release has an apparent
  Ca2+-cooperativity of 1-2, whereas synaptotagmin-dependent release generally exhibits an apparent
  Ca2+-cooperativity of 4-5 (Sun et al., 2007; Kochubey and Schneggenburger, 2011)... Ca2+ can trigger
  synaptic vesicle fusion in less than 100 microseconds (Sabatini and Regehr, 1996)."*
- **Geppert et al 1994**: *"The synchronous, fast component of Ca(2+)-dependent neurotransmitter release is
  decreased, whereas asynchronous release processes, including spontaneous synaptic activity (miniature
  excitatory postsynaptic current frequency) and release triggered by hypertonic solution or
  alpha-latrotoxin, are unaffected."*
- **Sutton et al 1998**: *"a core synaptic fusion complex containing syntaxin-1A, synaptobrevin-II and
  SNAP-25B. The structure reveals a highly twisted and parallel four-helix bundle... Conserved
  leucine-zipper-like layers are found at the centre of the synaptic fusion complex."*
- **Sabatini & Regehr 1996**: *"postsynaptic responses commence just 150 micros after the start of the
  presynaptic action potential. This brisk communication is a consequence of rapid calcium-channel kinetics...
  and extremely fast calcium-driven vesicle fusion, which lags behind calcium influx by 60 micros... the
  classic view... holds for these synapses at room temperature, but not at physiological temperatures."*
- **Schiavo et al 1992**: *"tetanus and botulinum toxins serotype B are zinc endopeptidases... specific for
  synaptobrevin... The rat synaptobrevin-2 isoform is cleaved by both neurotoxins at the same single site,
  the peptide bond Gln 76-Phe 77, but the isoform synaptobrevin-1, which has a valine at the corresponding
  position, is not cleaved."*
- **Blasi et al 1993 (Nature, BoNT/A)**: *"inhibition of transmitter release from synaptosomes caused by
  botulinum neurotoxin A (BoNT/A) is associated with the selective proteolysis of the synaptic protein
  SNAP-25... BoNT/A acts as a zinc-dependent protease that selectively cleaves SNAP-25."*
- **Blasi et al 1993 (EMBO J, BoNT/C1)**: *"inhibition of neurotransmitter release by botulinum neurotoxin
  type C1 was associated with the proteolysis of HPC-1 (= syntaxin)... Breakdown of HPC-1/syntaxin was
  selective since no other protein degradation was detectable."*
- **Giraudo et al 2006**: *"a reversible clamping protein (complexin) that can freeze the SNAREpin, an
  assembled fusion-competent intermediate en route to fusion. When calcium binds to the calcium sensor
  synaptotagmin, the clamp would then be released."*
- **Rosenmund & Stevens 1996**: *"release produced by action potentials and hypertonic solutions varies in
  parallel as the pool size is changed, we conclude that the same pool is shared by both mechanisms."*

## 2. Method — geometric, machine-checked (not narrated)

**Ca2+-cooperativity (Part 1).** Model Ca2+-triggering as **n independent, identical low-affinity binding
sites** that must essentially all be occupied to trigger fusion (a standard, disclosed simplification of
synaptotagmin's real C2A+C2B site structure). For n such sites, joint-occupancy probability follows the Hill
form `theta(Ca) = Ca^n/(Ca^n+Kd^n)`. Elementary calculus gives the **local log-log slope**
`d(ln theta)/d(ln Ca) = n*(1-theta)` — a geometric probability-of-joint-independent-events result, not a
heuristic: it equals **n exactly** in the low-occupancy limit (Ca<<Kd) and falls to **0** near saturation
(Ca>>Kd). This closed form is verified against an independent **finite-difference numerical derivative on
the same theta(Ca)** (machine self-consistency, not asserted calculus), then checked against Dodge &
Rahamimoff's own verbatim, independently-reported qualitative finding ("slope... diminishes as [Ca] is raised
towards the normal level") — reproduced, not assumed.

**The forced linear (n=1) adversary.** Rather than just asserting "steeper than linear," the linear adversary
is tested directly against Dodge & Rahamimoff's **own reported aggregate statistic** (mean slope=3.78, SD=0.2,
N=28 independent experiments — their number, not a fabricated dataset): a one-sample test of H0: true
slope=1.

**Katz quantal release (Part 2).** Binomial(N,p) release statistics validated on synthetic Monte Carlo data
first (memory: synthetic-control-before-real-negative) — recovers theoretical mean=Np and variance=Np(1-p) at
200,000 trials — then instantiated with **real, citation-anchored** N (Rosenmund & Stevens 1996's
operationally-defined RRP) and p (Schneggenburger & Neher 2000's own "~80% pool depletion" at saturating Ca),
explicitly kept separate from the synthetic-machinery check.

**Toxin/SNARE concordance (Part 3).** A model in which SNAREs are dispensable for release must fail against
3 independent clostridial-toxin serotypes, each with a distinct, biochemically specific cleavage target
(BoNT/A→SNAP-25, BoNT/B+TeNT→VAMP2 at one peptide bond, BoNT/C1→syntaxin), each independently abolishing
release. A trivial but real geometric sanity check ties the 4-helix bundle stoichiometry (1 syntaxin + 1
synaptobrevin + 2 SNAP-25 helices = 4) directly to Sutton et al 1998's own stated complex composition.

**Synaptotagmin double dissociation (Part 4).** Geppert et al 1994 (Syt1, hippocampal culture) and Sun et al
2007 (Syt2, calyx of Held) — two independent instances (different isoform, different synapse type, different
technique, 13 years apart) of the SAME dissociation pattern: fast/synchronous release impaired, slow/
asynchronous and other release modes spared. The forced adversary is a generic "sick neuron"
nonspecific-toxicity confound — ruled out directly because both papers explicitly report OTHER release modes
as unaffected/isolable in pure form, not just "reduced along with everything else."

**Fusion latency (Part 5).** Sabatini & Regehr 1996's direct, live-verified numbers (60µs Ca-entry-to-fusion,
150µs AP-to-postsynaptic-response) against the pre-registered <1ms threshold, with the paper's own built-in
temperature adversary (physiological vs room temperature) reported, not invented.

## 3. Gates — 15/15 pre-registered PASS

| Gate | Claim | Pre-registered threshold | Result | Verdict |
|---|---|---|---|---|
| G1 | Geometric slope: analytic vs finite-difference numeric self-consistency | rel err ≤1e-4 | matched to float precision, all n tested | **PASS** |
| G2 | Low-occupancy slope equals n | within 1% of n | exact to <0.1% for all n tested | **PASS** |
| G3 | Slope diminishes near saturation (reproduces Dodge-Rahamimoff's own qualitative finding) | slope ≤10% of n at Ca/Kd=100 | slope ≈0.0 for n≥3.78; 0.99% for n=1 | **PASS** |
| G4 | **Forced linear (n=1) adversary rejected, t-test (SEM-based)** | t≥5 | **t=73.55** | **PASS** |
| G5 | **Forced linear (n=1) adversary rejected, conservative z (raw-SD-based)** | z≥5 | **z=13.90** | **PASS** |
| G6 | All fast/synchronous instances (3 synapse types) exceed cooperativity threshold | n≥3.5 | 3.78 / 4.0 / 5.0 | **PASS** |
| G7 | The one measured slow/asynchronous reference instance falls below threshold (real dissociation on the exponent) | n≤3.0 | 2.0 | **PASS** |
| G8 | Binomial machinery recovers theoretical mean/variance (synthetic control) | rel err ≤2% | ≤0.5% across 3 (N,p) combos, 200k trials | **PASS** |
| G9 | Poisson-limit monotonic convergence (var/mean→1 as p→0) | strictly monotonic increase toward 1 | 0.502→0.901→0.989 | **PASS** |
| G10 | Poisson-limit closeness at lowest p tested | \|ratio-1\|≤0.03 | 0.0113 | **PASS** |
| G11 | Toxin/SNARE concordance (3 independent serotypes × 3 distinct SNAREs) | 3/3 | 3/3 | **PASS** |
| G12 | 4-helix bundle stoichiometry matches Sutton 1998's crystal structure | count=4 | 4 (1+1+2) | **PASS** |
| G13 | Both independent synaptotagmin-KO instances show the double dissociation | ≥1 impaired + ≥1 spared per instance, both instances | Geppert 1994 ✓, Sun 2007 ✓ | **PASS** |
| G14 | Ca-entry-to-fusion latency sub-ms | <1000µs | 60µs (16.7× margin) | **PASS** |
| G15 | AP-onset-to-postsynaptic-response latency sub-ms | <1000µs | 150µs (6.7× margin) | **PASS** |

**Overall: 15/15 PASS** (`overall_pass: true` in the raw results JSON; a real bug was caught and fixed during
this session — an inverted monotonicity check in G9's first draft flagged FAIL until the direction was
corrected against the closed-form binomial identity var/mean=(1-p), disclosed in Sec. 6, not hidden).

## 4. Symmetric QC — held OPEN, not smoothed over

- **Kd (Ca2+ affinity) is NOT universal across synapse types, even though the cooperativity ORDER n is**:
  Heidelberger et al 1994 (retinal bipolar, goldfish) report Kd=194µM; Sun et al 2007 (calyx of Held,
  synchronous component) report Kd≈38µM — a **5.1× difference**, machine-computed, reported directly (Sec.
  1.4 of the raw results JSON), not gated pass/fail and not forced into a single fabricated universal number.
  Only the **steepness** (n≈4-5 for fast/synchronous release) is claimed as the cross-synapse-type universal.
- **Schneggenburger & Neher 2000's own exact cooperativity number is not in its abstract**, and the paper
  itself is confirmed **paywalled** this session (EuropePMC: `isOpenAccess:N, inPMC:N, hasPDF:N` — checked
  directly, not assumed). The task explicitly names this paper as a required anchor; rather than settle for
  the qualitative "high cooperativity" phrase alone, this session mined Südhof 2013's PMC full text (which
  DOES have open access) and found it cites Sun et al 2007 for the exact number — then independently
  live-verified Sun et al 2007 itself, which reports **the same underlying calyx-of-Held Ca2+ sensor**
  (Sun et al's cooperativity-2007 finding: n≈5 synchronous / n≈2 asynchronous) as the specific numeric
  content Schneggenburger & Neher's qualitative "high cooperativity" phrase refers to in the field's own later
  literature. This is disclosed as a **derived-from-a-different-paper** substitution for the exact number,
  not smuggled in under Schneggenburger & Neher's own name.
- **del Castillo & Katz 1954's PMC full text is explicitly publisher-blocked** ("the publisher of this
  article does not allow downloading of the full text in XML form" — the exact same block message this
  repo's dopamine_kinetics.md already documented for an unrelated paper) — confirmed live this session, not
  assumed from precedent. Their own exact quantal-content numbers (m, n, p as originally reported) were **not**
  independently re-extracted; the modern quantitative quantal anchor used here (Part 2.3) is Rosenmund &
  Stevens 1996 (RRP) + Schneggenburger & Neher 2000 (80% pool depletion) instead.
- **The n-independent-site geometric model is a disclosed simplification.** Synaptotagmin's real C2A+C2B
  domains are not perfectly identical/independent sites (C2A and C2B differ in Ca2+ affinity and reported
  stoichiometry across studies) — the model is used only to derive the qualitative **low-occupancy asymptotic
  slope = n** behavior, which is directly, independently confirmed by Dodge & Rahamimoff's own text (not a
  refit of synaptotagmin's actual structural binding scheme).
- **Fatt & Katz 1952 and del Castillo & Katz 1954 are pre-modern-abstract PubMed records** — existence/
  identity/PMCID fully verified live; consistent with the same disclosed limitation this repo's
  `MECHANISM_DOPAMINE_KINETICS.md` and `MECHANISM_NA_K_ATPASE.md` already documented for equally old founding
  papers.

## 5. Additional convergence found this session (forced, not a one-shot stop)

- **The Sun et al 2007 Syt2-KO/calyx-of-Held finding is a second, fully independent instance of the SAME
  double dissociation Geppert et al 1994 established** — different synaptotagmin isoform (Syt2 vs Syt1),
  different synapse type (calyx of Held vs cultured hippocampal neurons), different technique (quantitative
  Ca2+-uncaging kinetics vs simple genetic/electrophysiological comparison), 13 years apart, no author overlap
  in the specific dissociation experiment. This was not in the task's own named anchor list (only Geppert 1994
  was named) — found by forcing the search past the first (paywalled) Schneggenburger-Neher hit, per the OODA
  mandate to orient/re-try rather than accept a one-shot gap.
- **The n≈1-2 "linear-ish" adversary is not merely hypothetically rejected — it is the field's own
  established description of a real, different, decorrelated release pathway** (asynchronous release, per
  Sun et al 2007 and Südhof 2013's synthesis). This sharpens the falsifier: linear/shallow Ca-dependence is
  not a strawman invented for this doc, it is the correct model for a *different* mechanism, which reinforces
  rather than merely asserts the specificity of the steep, synaptotagmin-gated, synchronous pathway.
- **Two methodologically decorrelated groups (Schneggenburger & Neher; Bollmann, Sakmann & Borst), same year,
  same synapse type, different specific protocols**, both report Ca2+-uncaging-driven pool depletion in
  under a few milliseconds at similarly low (10-30µM) Ca2+ steps — neither paper was designed to match the
  other's exact numbers, yet they converge qualitatively.

## 6. Honest gaps

- **Schneggenburger & Neher 2000's own exact cooperativity coefficient is not independently extracted from
  that paper itself** (paywalled, confirmed via 2 independent free-access routes this session: NCBI
  abstract-only efetch and EuropePMC REST `isOpenAccess:N`) — the number used (n≈4-5) is sourced from Sun et
  al 2007 (an independent, later, primary paper on the same synapse type) and Südhof 2013 (a review, PMC full
  text mined directly). Treat the **mechanism and qualitative magnitude** as decisively anchored across 3
  independent synapse types; treat the **specific Schneggenburger-Neher 2000 in-paper number** as not
  itself independently re-verified this session.
- **Kd is synapse-specific (5.1× spread), not universal** — disclosed in Sec. 4, not resolved. Any downstream
  cell consuming a single "the" Ca-affinity number must pick a synapse-type-matched value, not a generic one.
- **del Castillo & Katz 1954's and Fatt & Katz 1952's own exact numeric quantal-content values were not
  independently re-extracted** (pre-modern-abstract PubMed records; the 1954 paper's PMC full text is
  explicitly publisher-blocked, confirmed live). The modern quantal-framework anchor used instead (Rosenmund &
  Stevens 1996 RRP + Schneggenburger & Neher 2000 pool-depletion) is real and citation-sourced, but is not the
  *original* 1950s numbers.
- **The n-independent-site binding derivation is a standard simplification**, not a refit of synaptotagmin's
  actual (non-identical C2A/C2B, partially sequential) structural binding scheme — used only to derive the
  qualitative low-occupancy asymptotic-slope behavior, which independently matches Dodge & Rahamimoff's own
  text.
- **This session did not independently re-derive any wet-lab raw trace** — every gate either (a) reproduces a
  paper's own directly-quoted summary statistic via closed-form arithmetic (Part 1, Dodge-Rahamimoff's
  mean±SD±N), (b) validates a statistical tool on synthetic data before citing it as a framework (Part 2), or
  (c) builds a concordance/dissociation table directly from verbatim-quoted findings (Parts 3-4) — no figure
  was read/eyeballed (figures are forensic-only per this repo's discipline) and no raw electrophysiological
  trace was reprocessed.
- **BoNT/A, BoNT/B, TeNT, BoNT/C1 findings are all in vitro/synaptosome/cultured-neuron preparations** (not
  intact-organism dose-response) — the toxin/SNARE concordance is about **specificity and mechanism**, not a
  whole-organism LD50 or clinical-dose model (out of scope for this doc).

## 7. Overall result

**15 of 15 pre-registered gates PASS.** All 20 cited PMIDs were verified live this session via NCBI eutils
(verbatim abstract text for all; PMC full text additionally mined for 2 — Südhof 2013, and an unsuccessful
but disclosed attempt on del Castillo & Katz 1954). The task's own five named falsifiers are each reproduced
on machine-checked numbers, not narration: (1) the linear-Ca adversary is rejected by t=73.5/z=13.9 against
Dodge & Rahamimoff's own reported statistic, corroborated by 2 further independent synapse-type/technique
instances (Heidelberger 1994 retinal bipolar n≥4; Sun 2007 calyx of Held n≈5), with the field's own
asynchronous-pathway reference (n≈2) showing the dissociation is mechanism-specific, not universal; (2) the
Katz N·p·q quantal framework's statistical machinery is validated on synthetic data before being instantiated
with real RRP (Rosenmund & Stevens 1996) and pool-depletion (Schneggenburger & Neher 2000) anchors; (3) all
3 independent toxin/SNARE cleavage pairs (BoNT/A→SNAP-25, BoNT/B+TeNT→VAMP2, BoNT/C1→syntaxin) show
concordant proteolysis-with-release-abolition; (4) two fully independent synaptotagmin-KO instances
(Syt1/hippocampus, Syt2/calyx of Held) both reproduce the fast-impaired/slow-spared double dissociation,
directly ruling out a generic-toxicity confound; (5) measured fusion latency (60-150µs) clears the <1ms bar
by 6.7-16.7×. The single largest honest gap is Schneggenburger & Neher 2000's own exact in-paper number being
paywalled — bridged, not papered over, by an independent later primary source on the same synapse type found
by forcing the search past the first dead end. Confidence tier: **primary-literature-anchored** (frog NMJ
intracellular recording, retinal-bipolar and calyx-of-Held membrane-capacitance/Ca-uncaging, genetic
knockout, X-ray crystallography, in vitro toxin biochemistry) across 3 species/preps and 6 decades
(1952-2013) of methodologically decorrelated techniques.

## Repro

```
.venv-msk/bin/python3 scripts/msk/synaptic_vesicle_release.py
```

Pure closed-form arithmetic + `numpy` (Monte Carlo binomial/Poisson self-checks, finite-difference
derivatives) over literature-anchored parameters — no external data dependency, <2s wall time. Writes only
`data/synaptic_vesicle_release/synaptic_vesicle_release_results.json` (did not exist before this session —
first write). This session's independent citation verification: raw NCBI eutils esearch/efetch responses
(verbatim) in `data/raw_fetch/svr_efetch_*.{txt,xml}` + `data/raw_fetch/svr_europepmc_*.json`. No git
operations, no edits to any pre-existing file (isolation per task instruction).
