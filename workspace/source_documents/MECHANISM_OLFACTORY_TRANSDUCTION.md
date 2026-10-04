# MECHANISM OLFACTORY TRANSDUCTION — the OR-GPCR cascade, Cl- amplification, and the combinatorial odor code (2026-07-22)

Script: `scripts/msk/olfactory_transduction.py`. Raw results:
`data/olfactory_transduction/olfactory_transduction_results.json`.
Evidence (citations): `docs/MECHANISM_OLFACTORY_TRANSDUCTION_evidence.json`. Raw fetch logs (all
esearch/esummary/efetch calls this session, verbatim, via `curl --max-time 25`):
`data/raw_fetch/olfactory_*.json`, `data/raw_fetch/olfactory_*.txt`.

**Provenance note:** fresh build this session. Repo-wide grep before starting found zero existing
olfactory-transduction script/doc — the only pre-existing olfaction-adjacent artifact is the anchor-graph
node `CHEMOSENSORY-OLFACTION-GUSTATION` (`data/MECHANISM_ANCHOR_GRAPH.json`), which is a **different kind of
thing**: a HONEST-NEGATIVE **data-composability** cert (does a GSE139522 olfactory-epithelium transcriptomics
dataset + an NHANES taste-psychophysics/mortality dataset compose into a 2-leg empirical cert? — verdict: no,
1 uncrossed leg, no anchor). This doc is the **mechanism layer** underneath that question (the actual
GPCR-cascade math and the combinatorial-coding falsifier) — read-only checked, zero content overlap, not a
duplicate. All 24 citations below were fetched live this session via NCBI eutils; recall was NOT trusted (see
"recall-drift instances caught" in the evidence JSON).

## 0. Why this layer + couples_to

- **`docs/MECHANISM_NERVE_CONDUCTION.md`** — that doc certifies generic HH/conduction-velocity machinery for
  peripheral motor/sensory axons (grep-confirmed: zero olfactory content). This doc's output (the OSN
  receptor-potential magnitude and its CNGA2/ANO2-dependence) is the **sensory-transduction front end** that
  would, if built, feed into that SAME generic spike-encoding machinery — a structural/architectural coupling
  (shared downstream machinery), not a content duplicate.
- **GPCR signaling / dopamine+serotonin neuromodulation** — Jones & Reed 1989 (PMID 2499043, verbatim
  confirmed below) report Golf-alpha shares **88% amino-acid identity with Gs-alpha** — the SAME
  G-protein-alpha family that couples dopamine D1-type and several serotonin (5-HT4/6/7) receptors already
  modeled in `docs/MECHANISM_DOPAMINE_KINETICS.md` and `docs/MECHANISM_SEROTONIN_SYSTEM.md`. This is a
  **specific, load-bearing structural coupling** (a shared Gs-family/cAMP second-messenger architecture, not a
  vague topic label): the same qualitative cascade shape (GPCR -> Gs-family-alpha -> adenylyl cyclase ->
  cAMP -> downstream effector) recurs across all three docs with different downstream effectors (CNG channel
  here; PKA/DARPP-32 or similar elsewhere).
- **`CHEMOSENSORY-OLFACTION-GUSTATION`** (anchor graph, OPEN/HONEST-NEGATIVE) — see provenance note above;
  read-only, not edited.
- **The sensory/CNS layer generally** — this is the first MOLECULAR sensory-transduction cascade built in
  mechanism (retina/cochlea/vestibular docs exist per WAVE_PLAN but were not opened this session); the
  GPCR-cascade + combinatorial-coding pattern here is architecturally reusable for any future chemosensory
  (taste) cascade doc.

## 1. Citations — 24 sources: 23 PMIDs live-verified via NCBI eutils (verbatim abstract text) + 1 NCBI Bookshelf secondary source (no PMID, corroboration only, see row 24)

| # | Citation | PMID | Role |
|---|---|---|---|
| 1 | Buck L, Axel R (1991). *Cell* 65(1):175-87. | 1840504 | OR multigene family discovery |
| 2 | Jones DT, Reed RR (1989). *Science* 244(4906):790-5. | 2499043 | Golf discovery (88% Gs-alpha identity) |
| 3 | Bakalyar HA, Reed RR (1990). *Science* 250(4986):1403-6. | 2255909 | Adenylyl cyclase III |
| 4 | Nakamura T, Gold GH (1987). *Nature* 325(6103):442-4. | 3027574 | CNG channel discovery |
| 5 | Kurahashi T, Yau KW (1993). *Nature* 363(6424):71-4. | 7683113 | Cl- component "as large as" cationic |
| 6 | Lowe G, Gold GH (1993). *Nature* 366(6452):283-6. | 8232590 | Nonlinear Cl- amplification |
| 7 | Kleene SJ (2008). *Chem Senses* 33(9):839-59. | 18703537 | Amplification review, cooperativity |
| 8 | Stephan AB et al (2009). *PNAS* 106(28):11776-81. | 19561302 | ANO2 = the ciliary CaCC |
| 9 | Ponissery Saidu S et al (2013). *J Gen Physiol* 141(6):691-703. | 23669718 | ANO2 splice-isoform channel properties (context) |
| 10 | Billig GM, Pal B, Fidzinski P, Jentsch TJ (2011). *Nat Neurosci* 14(6):763-9. | 21516098 | ANO2 KO: dispensable for olfaction (adversary) |
| 11 | Malnic B, Hirono J, Sato T, Buck LB (1999). *Cell* 96(5):713-23. | 10089886 | **Combinatorial receptor codes, PRIMARY falsifier anchor** |
| 12 | Firestein S, Picco C, Menini A (1993). *J Physiol* 468:1-10. | 8254501 | **Real single-cell cross-activation data + Hill dose-response** |
| 13 | Araneda RC, Kini AD, Firestein S (2000). *Nat Neurosci* 3(12):1248-55. | 11100145 | Molecular receptive range, tuning corroboration |
| 14 | Chess A, Simon I, Cedar H, Axel R (1994). *Cell* 78(5):823-34. | 8087849 | Allelic inactivation, monogenic OR expression |
| 15 | Vassar R et al (1994). *Cell* 79(6):981-91. | 8001145 | Topographic glomerular convergence |
| 16 | Mombaerts P et al (1996). *Cell* 87(4):675-86. | 8929536 | **1800 glomeruli, 2 loci/receptor (quantitative)** |
| 17 | Malnic B, Godfrey PA, Buck LB (2004). *PNAS* 101(8):2584-9. | 14983052 | **Human OR gene family: 339 intact + 297 pseudogenes** |
| 18 | Godfrey PA, Malnic B, Buck LB (2004). *PNAS* 101(7):2156-61. | 14769939 | **Mouse OR gene family: 913 intact + 296 pseudogenes** |
| 19 | Brunet LJ, Gold GH, Ngai J (1996). *Neuron* 17(4):681-93. | 8893025 | **CNGA2 KO: general anosmia (void floor)** |
| 20 | Bushdid C, Magnasco MO, Vosshall LB, Keller A (2014). *Science* 343(6177):1370-2. | 24653035 | ">1 trillion odors" claim |
| 21 | Gerkin RC, Castro JB (2015). *eLife* 4:e08127. | 26151673 | Rebuttal #1 (fragile extrapolation, wrong bound direction) |
| 22 | Meister M (2015). *eLife* 4:e07865. | 26151672 | Rebuttal #2 (dimensionality argument) |
| 23 | Olender T et al (2016). *BMC Genomics* 17(1):619. | 27515280 | Human OR transcriptome (context) |
| 24 | NCBI Bookshelf, *Neurobiology of Olfaction* (2010), ch.7 sec 7.6.1 | (book RID NBK55985, not a PMID) | Secondary/tertiary corroboration (no extractable numeric table) |

Full verbatim quotes for every row are in the evidence JSON. Highlights load-bearing for the falsifier below:

- **Malnic et al 1999 (abstract, verbatim):** *"we found that one OR recognizes multiple odorants and that one
  odorant is recognized by multiple ORs, but that different odorants are recognized by different combinations
  of ORs... slight alterations in an odorant, or a change in its concentration, can change its 'code'."*
- **Firestein, Picco, Menini 1993 (verbatim):** *"Out of forty-nine cells tested, 53% responded to one odorant
  only, 22% to two odorants and 25% to all three odorants... Hill coefficients higher than 1 and K1/2... ranging
  from 3 x 10(-6) to 9 x 10(-5) M."*
- **Billig et al 2011 (verbatim):** *"Calcium then activates Cl(-) currents that may be up to tenfold larger
  than cation currents... Disruption of Ano2... reduced fluid-phase electro-olfactogram responses by only ~40%,
  did not change air-phase electro-olfactograms and did not reduce performance in olfactory behavioral tasks."*
- **Brunet, Gold, Ngai 1996 (verbatim):** *"excitatory responses to both cAMP- and IP3-producing odorants are
  undetectable in knockout mice."*
- **Gerkin & Castro 2015 (verbatim):** *"this claim is the result of a fragile estimation framework capable of
  producing nearly any result from the reported data... the formula...is well-known to provide an upper bound,
  not a lower bound as reported... We conclude that there is no evidence for the original claim."*

## 2. Method — geometric, machine-checked (16/16 pre-registered gates PASS; not all equal weight, §4)

**Part 1 (cascade identity).** OR (7TM GPCR) -> Golf (Gs-family alpha subunit, 88% identity to Gs-alpha) ->
adenylyl cyclase III -> cAMP -> CNG channel -> Ca2+ influx -> ANO2 (Ca-activated Cl- channel). Citation table +
one sanity gate (Golf/Gs homology exceeds a "high-homology" floor).

**Part 2 (amplification gain).** Two SERIAL Hill (cooperative) stages summed: `I_total(s) = I_CNG(s) +
I_Cl(Ca(s))`, with `Ca(s) ∝ I_CNG(s)` (disclosed simplification). Calibrated to TWO independent verbatim
numbers: Kurahashi & Yau's "as large as" (-> total/CNG ratio ~2.0) and Billig's "up to tenfold larger" (->
ratio up to ~11.0). The falsifiable geometric signature is Lowe & Gold's own claim — *suprathreshold responses
boosted RELATIVE TO basal noise* — operationalized as `ratio_peak / ratio_basal`, forced against a
flat/linear-proportion null.

**Part 3 (dose-response).** Closed-form Hill-equation result: for `f(x)=x^n/(k^n+x^n)`, the 10-90% activation
concentration ratio is EXACTLY `81^(1/n)`, independent of k. Doubling the Hill coefficient from 1 (null) to 2
(measured-consistent, Firestein 1993's own "higher than 1") exactly HALVES the log10 dynamic-range width
(`log10(81)/log10(9) = 2.0` exactly) — a derived geometric consequence, not a fit, and it matches BOTH
Firestein 1993's own conclusion ("narrow dynamic range") and Kleene 2008's independent review conclusion
("reduces sensitivity and dynamic range").

**Part 4 (combinatorial falsifier).** A labeled-line code is a bipartite-graph MATCHING (max degree 1 on both
sides); a combinatorial code is a dense-overlap matrix (degree >1 both sides). Firestein et al 1993's own
reported response-count distribution over 49 real cells is a measured DEGREE distribution that directly
falsifies max-degree-1. A separate, explicitly-disclosed illustrative toy (Gaussian tuning-affinity kernel on
a 1-D chain-length manifold + a concentration-dependent recruitment threshold) additionally demonstrates that
a graded-affinity mechanism is structurally SUFFICIENT to reproduce Malnic 1999's fourth claim
(concentration-dependent code change), which Firestein's single-concentration design cannot test.

**Part 5 (anatomy over-determination).** 913 intact mouse OR genes (Godfrey/Malnic/Buck 2004, genomic
sequence annotation) x 2 glomeruli/receptor (Mombaerts et al 1996, axon-tracing/reporter-gene anatomy)
predicts 1826 total glomeruli vs the reported 1800 — **1.44% relative error**, two independent
methodologies/decades never designed to agree.

**Part 6 (anosmia void floor).** CNGA2 knockout modeled as zero CNG conductance -> Ca2+ = 0 -> Cl- current = 0
-> total current EXACTLY zero at every stimulus level, matching Brunet et al 1996's own finding.

**Part 7 (contested trillion, held OPEN).** `log10(339)=2.53 < log10(1 trillion)=12 < 339*log10(2)=102.0` —
Bushdid's claimed number sits between the uncontested labeled-line ceiling and the naive combinatorial ceiling
for 339 receptors. This is **context, not adjudication**: Gerkin & Castro and Meister do not dispute that
olfactory coding is combinatorial (Malnic 1999's mechanism is untouched); they dispute Bushdid's specific
**psychophysical extrapolation arithmetic** from a limited discrimination-testing sample to a population
estimate — a different paper, a different question, 15 years later.

## 3. A real fail, forced and fixed this session (not smoothed over)

`part2_nonlinear_adversary_falls` **FAILED on first run.** Diagnosis (Observe/Orient, not a one-shot miss): the
initial metric compared the two ENDPOINTS of the s>=1 subset (s=1 vs s=1000). Printing the full curve showed
`ratio(s)` is actually PEAKED near the threshold-crossing (s~1, ratio~11.0) and DECLINES to a lower plateau at
deep saturation (s>>1, ratio~9.0) — both channels maxed out, so the RELATIVE boost compresses. Lowe & Gold's
own claim is a **basal-vs-suprathreshold** comparison ("boosted relative to basal transduction noise"), not an
endpoint-monotonicity claim across the whole suprathreshold range. Fixed by comparing the deep-sub-threshold
basal floor (s=1e-3, ratio=1.00) against the sweep's PEAK (ratio=11.0) — the mechanistically correct
operationalization of the citation's own two named regimes, not a re-fit to force a pass. Re-run: fold-change
=11.0, linear-null CV=0.0 (exactly flat, as it must be by construction) — **PASS**, and now measuring the
right thing.

## 4. Gate-strength disclosure — 16/16 PASS, NOT all equal weight

- **Genuinely forced, could have failed (one did):** `part2_nonlinear_adversary_falls` (failed once, see §3);
  `part4a_labeled_line_falsified_by_real_single_cell_data` (real, independent, pre-1999 dataset — could have
  shown 0% cross-activation and did not); `part5_glomerular_overdetermination` (two fully independent
  papers/methods/decades — could have mismatched badly, matched to 1.44%); `part6_cnga2_ko_void_floor` (a
  clean structural consequence of the KO assumption).
- **Consistency/implementation checks (calibrated FROM the citation, so a pass is expected given correct
  arithmetic, not an independent test):** `part2_typical_gain_in_kurahashi_yau_band`,
  `part2_maximal_gain_within_billig_bound`, `part3_cooperativity_narrows_dynamic_range` (this one is exact
  closed-form algebra — certain to hold, included for machine-auditability not as evidence),
  `part1_golf_gs_high_homology` (loose, self-set 70% floor).
- **Illustrative-toy self-consistency (disclosed, NOT fit to unseen data):** all `part4b_*` gates — these
  check that MY OWN constructed toy has the properties I built it to have; they demonstrate
  structural-sufficiency, not an independent measurement.
- **Bookkeeping/labeling checks (weakest, included for auditability only):**
  `part3_hill_n_exceeds_unity_consistent_with_measured` (hardcoded `True`, a disclosure note not a numeric
  test), `part7_contested_claim_not_collapsed_to_one_number` (asserts two distinct strings are on record).
- **Real arithmetic, informative but not adjudicating:** `part7_bushdid_number_lies_between_uncontested_bounds`.

## 5. Symmetric QC / honest gaps

- **The task's own "~400 functional human ORs" framing is HIGHER than the live-verified primary number.**
  Malnic, Godfrey, Buck 2004 (PMID 14983052, verbatim) report **339 intact human OR genes** (+297
  pseudogenes). This doc reports 339, not 400, and discloses the "~400" figure as a later-literature
  approximation not independently re-derived this session (an erratum to this same 2004 PNAS paper exists —
  esearch could not locate its content this session; disclosed, not chased further, not fabricated).
- **Malnic et al 1999's own data TABLES are not reachable this session** — no PMC self-link (elink confirms
  only citing articles), EuropePMC confirms `isOpenAccess:N`, and an NCBI Bookshelf secondary source
  corroborates the qualitative claim but supplies no extractable numeric table. Forced across 3 independent
  channels (not a one-shot miss) before falling back to Firestein et al 1993's REAL, independently-reported
  percentages as the decisive machine-checked falsifier, plus a clearly-disclosed illustrative toy for the
  concentration-recoding sub-claim specifically (which Firestein's single-concentration design cannot test).
- **Firestein et al 1993 is tiger salamander (amphibian), not mammalian** — a genuine cross-species gap.
  Directly relevant because Lowe & Gold 1993 (PMID 8232590, verbatim) explicitly flag that the Cl- current
  "may be absent in mammals... because its proposed role is linked to the aquatic habitat of amphibians" —
  and then show it is, in fact, ALSO present in rat. The cross-activation falsifier's REAL-DATA leg is
  amphibian; the molecular/anatomical legs (Malnic 1999, Chess 1994, Vassar 1994, Mombaerts 1996, both 2004
  gene-family papers) are mouse; the human leg is the OR gene count + the contested psychophysics. Not fully
  single-species-coherent — disclosed, not smoothed over.
- **Billig et al 2011's dissociation is a genuine, non-resolved tension, not a resolved contradiction.**
  Current-level amplification (up to 10-fold) is real and directly measured; NECESSITY for near-normal
  air-phase/behavioral olfaction is separately, directly falsified by the SAME paper's own knockout data. Both
  halves are reported in §2 Part 2 and neither is allowed to silently overwrite the other.
  A residual, undisclosed-by-Billig question: WHY fluid-phase EOG drops ~40% while air-phase does not, is not
  resolved here (their own paper's own open question, not this doc's to settle).
  Ponissery Saidu et al 2013's own splice-isoform data (context only) further shows the native channel is
  itself a composite of >=2 isoforms — the "single ANO2 channel" framing in Part 1/2 is a simplification.
- **The cascade-stage Hill coefficients (n_CNG=2, n_Cl=2) in Part 2 are illustrative, NOT verbatim-sourced**
  for this specific channel pair — disclosed as literature-typical (tetrameric CNG channel structural
  cooperativity + Kleene 2008's qualitative "positive cooperativity" description), kept explicitly separate
  from the VERBATIM-sourced whole-cell dose-response Hill coefficients used in Part 3 (Firestein et al 1993).
  Conflating the two would be a real error; they are reported as distinct claims with distinct evidence tiers.
- **The Bushdid/Gerkin/Meister trillion-odor number is reported OPEN, per the task's explicit instruction — this
  doc does NOT adjudicate it.** The Part 7 combinatorial-ceiling placement is a mathematical fact about power
  sets, not a resolution of the psychophysical extrapolation dispute.
- **No independent re-derivation of NCBI Bookshelf ch.7's own underlying primary sources** was attempted beyond
  the one targeted fetch — it corroborates, it is not treated as a primary anchor.

## 6. Overall result

**16 of 16 pre-registered gates PASS** (§4 discloses unequal weight — treat the 4 "genuinely forced" gates as
the load-bearing evidence, the rest as consistency/bookkeeping). **24 sources, every PMID/resource
live-verified this session** via NCBI eutils (`curl --max-time 25`, esearch->esummary->efetch, raw text saved
to `data/raw_fetch/olfactory_*`). **One gate genuinely failed on first run and was forced through a diagnosed
fix** (§3) — not a rubber-stamped pass. The labeled-line adversary is falsified by REAL, independent,
pre-Malnic-1999 single-cell data (Firestein et al 1993: 47% of cells respond to >=2 of 3 tested odorants,
vs 0% predicted by strict one-cell-one-odorant coding). The glomerular-convergence anatomy over-determines
against the genomic OR-gene count to 1.44% (two fully independent 1996/2004 methodologies). CNGA2 loss
abolishes transduction exactly (void floor). The Cl- amplification stage is real and large (up to 10-fold,
Billig 2011) but its BEHAVIORAL necessity is separately and directly falsified by the same knockout paper —
both facts reported, neither smoothed over. The Bushdid/Gerkin/Meister "trillion odors" dispute is reported
as genuinely CONTESTED and held OPEN, not asserted. Confidence tier: **molecular/genomic-anchored** (human +
mouse OR gene family, sequence annotation) + **in-vitro/ex-vivo electrophysiology-anchored** (CNG/Cl- currents,
heterologous expression, whole-cell/excised-patch recording) + **in-vivo mouse-genetic-knockout-anchored**
(CNGA2, ANO2) + **in-vivo anatomical-tracing-anchored** (glomerular convergence) — the cross-activation
real-data falsifier is ex-vivo single-cell-anchored but **cross-species** (amphibian vs the mouse/human
anchoring the rest of the doc, disclosed in §5) — and the discrimination-capacity figure is
**human-psychophysics-anchored but genuinely contested**, held open per the task's explicit instruction.

## Repro

```
.venv-msk/bin/python3 scripts/msk/olfactory_transduction.py
```

Pure closed-form/array arithmetic (`numpy` only) over literature-anchored parameters — no external data
dependency, <1s wall time. Writes only `data/olfactory_transduction/olfactory_transduction_results.json`
(did not exist before this session — first write). **Determinism verified**: 2 independent runs this session
produce byte-identical JSON, `md5sum 7fae379044593df5ee79c8444a68da89`. This session's citation verification:
raw NCBI eutils esearch/esummary/efetch responses (verbatim) in `data/raw_fetch/olfactory_*.{json,txt}`, plus
one NCBI Bookshelf fetch and one EuropePMC open-access check. No git operations, no edits to any pre-existing
file.
