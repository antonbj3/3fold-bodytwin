# MECHANISM HEMATOPOIESIS — the HSC -> progenitor -> mature-cell hierarchy: quiescence, transplant clonality, and the amplification factor behind ~10^11-10^12 cells/day (2026-07-22)

Builds and MEASURES the root of the blood/immune system: the hematopoietic stem cell (HSC) ->
progenitor -> mature-cell hierarchy that produces ~10^11-10^12 cells/day from a pool of HSCs that is
simultaneously RARE (~1:10^3-10^6, assay-dependent, Sec.9) and mostly QUIESCENT (division interval
measured in **months**, not hours). Supplies the ROOT this repo's existing blood docs assumed but
never built: `MECHANISM_ERYTHROPOIESIS.md` models the EPO->RBC leg but not the HSC compartment
feeding it; `MECHANISM_PLATELET_HEMOSTASIS.md` explicitly disclosed that "megakaryocyte
thrombopoiesis, thrombopoietin-regulated" was "not modeled here... a natural next node" — this doc
builds that node, plus a new G-CSF->neutrophil production leg, plus the HSC-level quiescence/
transplant/amplification mechanism none of this repo's docs or graph nodes has executed. Script:
`scripts/msk/hematopoiesis.py`. Evidence: `data/hematopoiesis/hematopoiesis_results.json`.

**NO re-solve, no OpenSim.** Reads `erythropoiesis_results.json` read-only for its own
machine-derived RBC production rate (2.083×10¹¹ cells/day, Little's-Law) and reuses it directly as
this doc's EPO/RBC leg — that sibling script is neither modified nor recomputed.

**Confidence tier: in-vivo-anchored** (mouse label-retention/division-counting for quiescence;
mouse+human single-cell/limiting-dilution transplant for clonality; human whole-genome
somatic-mutation phylogenetics for HSC pool size; genetic knockout for G-CSF necessity; rabbit
busulfan-thrombocytopenia + in vitro clearance assay for the TPO mechanism) — same tier class as
this repo's `MECHANISM_ERYTHROPOIESIS.md` and `MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md`.

## 0. Scope note — graph-node disambiguation, checked live before writing a line of code

`data/MECHANISM_ANCHOR_GRAPH.json` (998 nodes) was searched live (whole-word regex over
HSC/HEMATOPOIETIC/PROGENITOR/CFU/MYELOID/LYMPHOID/NEUTROPHIL/THROMBOPOIETIN/MEGAKARYOCYTE/CLONAL
HEMATOPOIESIS, 86 raw hits triaged) before starting, per this repo's own established convention
(`MECHANISM_ERYTHROPOIESIS.md` §0, `MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md` §11). Five
HSC/hematopoiesis-**adjacent** nodes exist; this doc is deliberately **not** a fold into any of
them, and **none was edited**:

- **`HEMATO-HSC-CLONAL-AGING`** (OPEN, EMPIRICAL) — already carries real `fetched_live` datapoints
  (CHIP prevalence 9.5%/18.4%, PMID 25426837; aging-BMT chimerism 1.94×/1.52×-skew, PMID 11067876).
  The CHIP/aging-clonality layer. This doc independently **re-fetches** Jaiswal 2014 and Genovese
  2014 live (§10, Part 6) as a genuine cross-confirmation (not a rubber stamp — this repo's own
  memory explicitly warns a subagent's "confirmed live" claim has been fabricated before) and
  **cross-checks byte-for-byte against the existing node's stored numbers** (gate G6a, PASS) —
  then cites, does not re-derive, its aging-BMT leg.
- **`AUTO-CROSS-METHOD-HSC-CLONE-COUNT-RECONCILIAT[ION]`** (OPEN, **SEED-DESIGN — designed, cited
  literature, never executed**): "machine-check whether Sun2014(thousands)/Busch2015(~5e3)/
  Lee-Six2018+Mitchell2022(5e4-2e5) estimates reconcile... falsifier: reconciliation needs >2
  orders-of-magnitude free parameter with no biological grounding -> kill as discordant." **This
  script executes that exact pre-registered falsifier for the first time** (§7), using this
  session's own live-fetched numbers — not promoted to a graph edge (no `fold_gate_v2` write this
  session, isolation), but genuinely computed, not merely cited.
- **`AUTO-TEST-WHETHER-NICHE-CXCL12-SCF-GRADIENT-M[AGNITUDE]`** (OPEN, SEED-DESIGN, not executed) —
  the HSC-niche **causal-mechanism** layer (does CXCL12/SCF gradient magnitude predict quiescence
  depth). This doc's Part 1 falsifier is about the **measured quiescence rate itself**
  (label-retention), not the upstream niche-signaling cause — complementary, not resolved here.
- **`AUTO-CELL-CHIP-MDS-CAUSE-VS-MARKER-SYMMETRIC[-AUDIT]`** (OPEN, SEED-DESIGN) — the CHIP-to-MDS
  **treatment-implication** layer (aza/decitabine/MEDALIST trial endpoints). Out of scope; this
  doc's CHIP mention (§10) is the task's own symmetric-QC "hold open" item, not a treatment audit.
- **`ORG-BLOOD-COAG`** (OPEN, self-flagged `HONEST-NEGATIVE-CONFIRMED`: "both legs scRNA, no
  functional coag dataset") — a single-cell-transcriptomic **differentiation-fidelity** cert vs.
  plasma coagulation factor levels. Different scope from this doc's **population-kinetics**
  mechanism (production rates, not per-cell transcriptomic state).

This document **establishes new ground**: no existing node/doc builds the executed HSC
rarity/quiescence/transplant-clonality/amplification-factor mechanism, nor the G-CSF/TPO
production legs.

## 1. Geometric structure (derive from the geometry, not heuristics)

Three genuinely different pieces of geometry are doing the work here:

- **The hierarchy is a branching, age-structured amplification cascade**, not a single population.
  A tiny, slow-cycling root (HSC) feeds a chain of transit-amplifying compartments (MPP -> CMP/CLP
  -> committed progenitors -> precursors) whose OWN division rates are orders of magnitude faster
  than the root's — the same Little's-Law renewal-process identity `production_rate = pool /
  mean_lifespan` that `MECHANISM_ERYTHROPOIESIS.md` §1/§10 used for the terminal (RBC) leg is reused
  identically for two NEW legs here (neutrophil, platelet, §7) — three independent instances of the
  same geometric identity, not three different curve-fits.
- **Label-retention is a dilution-per-division geometry**, not a rate assumed a priori: a
  fluorescent/thymidine label halves (or dilutes) with every division, so the "dormant" subset
  identified by strong label retention is, BY CONSTRUCTION, the slowest-dividing tail of the
  distribution — Wilson 2008 and Foudi 2009 both build on this same halving-per-division logic but
  gate on **different marker stringencies**, which is exactly why their rate estimates diverge
  (§5). Bernitz 2016 uses a genuinely different geometry — a **cumulative division-counter** (not a
  rate-inferred-from-dilution model) — and independently lands on the same "~5" ceiling.
- **The single-hit/Poisson limiting-dilution law governs transplant clonality**: if each purified
  cell has an independent probability `p` of engrafting, a dose of `n` cells succeeds with
  probability `P(n) = 1 - (1-p)^n` — a clean closed form (§6) that explains why MORE cells can never
  give a LOWER success rate (a real, checked monotonicity gate), while also explaining why a
  `p` calibrated on one purification protocol/endpoint does NOT transfer cleanly to a different
  protocol/endpoint 8 years apart (§6's forced-adversary finding).
- **TPO senses demand via mass-action clearance, not a chemical threshold** — a structurally
  DIFFERENT geometry from EPO's O2-tension threshold sensor (§8): `dC/dt = P - k·M·C` gives a
  steady state `C_ss = P/(k·M)`, a clean reciprocal law with no hinge/threshold at all, verified
  monotonic and exactly reciprocal (gates G4b/G4c) — directly matching Kuter & Rosenberg 1995's own
  measured "reciprocal relationship" and their own in-vitro demonstration that platelets physically
  REMOVE TPO from plasma.

## 2. Method, in one paragraph

Six parts, each independently falsifiable, each with a pre-registered threshold. **(1)** Does the
model reproduce the measured HSC quiescence/division interval (145-193 d mouse, ~5 divisions/
lifetime) — forced against Wilson 2008's OWN naive-adversary quote ("the entire HSC pool turns over
every few weeks") and against Foudi 2009's own (non-overlapping, diagnosed) faster interval? **(2)**
Does a single HSC reconstitute the entire blood system (Osawa 1996's directly-measured 21%,
Spangrude 1988's 30-cell dose-response, Notta 2011's human cross-species corroboration) — with a
BONUS forced-adversary test of whether a constant-probability Poisson model transfers across
purification eras (it does not, diagnosed not hidden)? **(3)** The decorrelated check: total daily
output (RBC+neutrophil, the task's own literal anchor) vs. the tiny HSC pool — executing the
graph's own pre-registered, never-before-run HSC-clone-count reconciliation falsifier, then
deriving the amplification factor (~29-31 doublings) and stress-testing it with a full-range
void-floor sweep. **(4)** Contrasting the geometric FORM of the three demand signals (EPO threshold
vs. G-CSF genetic-necessity vs. TPO mass-action clearance). **(5)** Quantifying, not hand-waving,
the "HSC frequency is assay-dependent" symmetric-QC item (>2 orders of magnitude spread across two
live-verified purification methods). **(6)** Reusing (not re-deriving) the existing CHIP graph
node, independently re-confirmed live, held OPEN per the task's own instruction.

## 3. Citations — 25 papers, every PMID resolved LIVE this session via NCBI eutils

Method: `esearch` (title-anchored, NOT trusting a recalled PMID) -> `esummary` (title/author/
journal/date cross-check) -> `efetch` (raw abstract text, verbatim numbers extracted). Recall-drift
was **caught and corrected multiple times** during this process (disclosed, not hidden — matches
this repo's own `MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md` §3 precedent and its own isolation
warning that "~62%" of recalled PMIDs drift): Notta 2011's actual PMID (21737740) differs from an
initially-recalled 21998251 (refuted by direct esummary cross-check — 21998251 resolves to an
unrelated paper); the title "Neutrophil kinetics in man" alone returned an unrelated 1988 GM-CSF
paper before an author-anchored search (`Dancey JT[Author] AND 1976[pdat]`) found the real PMID
956397; a bare "Spangrude/Weissman 1988" title search matched THREE different wrong Spangrude
papers (Sca-1/Sca-2, thymic subsets) before a direct-PMID cross-check confirmed the originally-
recalled 2898810 had in fact been correct all along (a genuine round-trip verification, not merely
assumed).

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Wilson A, Laurenti E, Oser G, et al (2008). Hematopoietic stem cells reversibly switch from dormancy to self-renewal during homeostasis and repair. *Cell* 135(6):1118-29. | **19062086**, 10.1016/j.cell.2008.10.048 | PRIMARY quiescence anchor: "d-HSCs divide about every 145 days, or five times per lifetime." Also supplies its OWN forced-adversary quote (§5). |
| 2 | Foudi A, Hochedlinger K, Van Buren D, et al (2009). Analysis of histone 2B-GFP retention reveals slowly cycling hematopoietic stem cells. *Nat Biotechnol* 27(1):84-90. | **19060879**, 10.1038/nbt.1517 | "~20% of HSCs divide at ≤0.8-1.8%/day" — a DIFFERENT (less-stringently-gated) label-retention estimate, diagnosed not smoothed over (§5). |
| 3 | Bernitz JM, Kim HS, MacArthur B, Sieburg H, Moore K (2016). Hematopoietic Stem Cells Count and Remember Self-Renewal Divisions. *Cell* 167(5):1296-1309. | **27839867**, 10.1016/j.cell.2016.10.022 | Independent (division-counting, not rate-modeling) method: "4 traceable symmetric self-renewal divisions... potential lost upon a fifth division." |
| 4 | Osawa M, Hanada K, Hamada H, Nakauchi H (1996). Long-term lymphohematopoietic reconstitution by a single CD34-low/negative hematopoietic stem cell. *Science* 273(5272):242-5. | **8662508**, 10.1126/science.273.5272.242 | PRIMARY transplant anchor: single purified mouse HSC -> **21%** long-term multilineage reconstitution. |
| 5 | Notta F, Doulatov S, Laurenti E, Poeppl A, Jurisica I, Dick JE (2011). Isolation of single human hematopoietic stem cells capable of long-term multilineage engraftment. *Science* 333(6039):218-21. | **21737740**, 10.1126/science.1201219 | HUMAN cross-species corroboration: single CD49f+ cell "highly efficient" at long-term multilineage grafts (qualitative-tier, exact % not in abstract). |
| 6 | Sun J, Ramos A, Chapman B, et al (2014). Clonal dynamics of native haematopoiesis. *Nature* 514(7522):322-7. | **25296256**, 10.1038/nature13824 | Native (non-transplant) steady-state hematopoiesis driven by "THOUSANDS of clones" — long-lived progenitors, not classical HSC, dominate day-to-day flux (qualitative-tier number, §7). |
| 7 | Busch K, Klapproth K, Barile M, et al (2015). Fundamental properties of unperturbed haematopoiesis from stem cells in vivo. *Nature* 518(7540):542-6. | **25686605**, 10.1038/nature14242 | "At least 30% or ~5,000 HSCs are productive in the adult mouse" — quantitative mouse pool-size anchor. |
| 8 | Lee-Six H, Øbro NF, Shepherd MS, et al (2018). Population dynamics of normal human blood inferred from somatic mutations. *Nature* 561(7724):473-478. | **30185910**, 10.1038/s41586-018-0497-0 | HUMAN active HSC pool: **50,000-200,000** (whole-genome somatic-mutation phylogenetics, n=140 colonies, 1 subject). |
| 9 | Mitchell E, Spencer Chapman M, Williams N, et al (2022). Clonal dynamics of haematopoiesis across the human lifespan. *Nature* 606(7913):343-350. | **35650442**, 10.1038/s41586-022-04786-y | Extends #8 across 10 subjects (0-81y): "stable population of 20,000-200,000 HSC/MPPs" (<65y); 17 mutations/HSC/year; CHIP-like clonal collapse >75y. |
| 10 | Abkowitz JL, Catlin SN, McCallie MT, Guttorp P (2002). Evidence that the number of hematopoietic stem cells per animal is conserved in mammals. *Blood* 100(7):2665-7. | **12239184**, 10.1182/blood-2002-03-0822 | Cat+mouse HSC ≈ 11,400±5,400 (conserved); OWN flagged tension: naive human extrapolation implies a frequency "20-fold less" than the NOD/SCID assay (§7). |
| 11 | Catlin SN, Busque L, Gale RE, Guttorp P, Abkowitz JL (2011). The replication rate of human hematopoietic stem cells in vivo. *Blood* 117(17):4460-6. | **21343613**, 10.1182/blood-2010-08-303537 | Human HSC replicate "on average once every 40 weeks (range, 25-50 weeks)" — the mouse-vs-human division-rate anchor (§5). |
| 12 | Spangrude GJ, Heimfeld S, Weissman IL (1988). Purification and characterization of mouse hematopoietic stem cells. *Science* 241(4861):58-62. | **2898810**, 10.1126/science.2898810 | "Thirty of these cells are sufficient to save 50 percent of lethally irradiated mice, and to reconstitute all blood cell types" — a second, dose-response transplant anchor (§6). |
| 13 | Kiel MJ, Yilmaz OH, Iwashita T, Terhorst C, Morrison SJ (2005). SLAM family receptors distinguish hematopoietic stem and progenitor cells and reveal endothelial niches for stem cells. *Cell* 121(7):1109-21. | **15989959**, 10.1016/j.cell.2005.05.026 | A THIRD, orthogonal marker scheme (CD150/CD244/CD48) — purification itself is a moving target (§9). |
| 14 | Doulatov S, Notta F, Laurenti E, Dick JE (2012). Hematopoiesis: a human perspective. *Cell Stem Cell* 10(2):120-36. | **22305562**, 10.1016/j.stem.2012.01.006 | Review, human-vs-mouse hierarchy perspective (qualitative-tier). |
| 15 | Morrison SJ, Weissman IL (1994). The long-term repopulating subset of hematopoietic stem cells is deterministic and isolatable by phenotype. *Immunity* 1(8):661-73. | **7541305**, 10.1016/1074-7613(94)90037-x | Marker-enriched population = 0.05% of marrow (1:2,000); only ~25% of THAT population's clonal reconstitutions are long-term — the marker-frequency anchor (§9). |
| 16 | Lieschke GJ, Grail D, Hodgson G, et al (1994). Mice lacking granulocyte colony-stimulating factor have chronic neutropenia... *Blood* 84(6):1737-46. | **7521686** | G-CSF-/- KO: neutrophils 20-30% of WT, progenitors -50%, gene-dosage in heterozygotes — the G-CSF causal (loss-of-function) anchor (§8). |
| 17 | Dancey JT, Deubelbeiss KA, Harker LA, Finch CA (1976). Neutrophil kinetics in man. *J Clin Invest* 58(3):705-15. | **956397**, 10.1172/JCI108517 | Total marrow neutrophils 7.70×10⁹/kg; postmitotic pool 5.59×10⁹/kg; transit time 6.60 d; production 0.85-0.87×10⁹/kg/day (2 independent methods converge); DF32P artifact-prone 1.62×10⁹/kg/day (§7). |
| 18 | de Sauvage FJ, Hass PE, Spencer SD, et al (1994). Stimulation of megakaryocytopoiesis and thrombopoiesis by the c-Mpl ligand. *Nature* 369(6481):533-8. | **8202154**, 10.1038/369533a0 | TPO discovery paper #1 (§8). |
| 19 | Kaushansky K, Lok S, Holly RD, et al (1994). Promotion of megakaryocyte progenitor expansion and differentiation by the c-Mpl ligand thrombopoietin. *Nature* 369(6481):568-71. | **8202159**, 10.1038/369568a0 | TPO discovery paper #2, names the molecule "thrombopoietin" (§8). |
| 20 | Kuter DJ, Rosenberg RD (1995). The reciprocal relationship of thrombopoietin (c-Mpl ligand) to changes in the platelet mass during busulfan-induced thrombocytopenia in the rabbit. *Blood* 85(10):2720-30. | **7742532** | THE mechanistic anchor: TPO "increased inversely and proportionally" with platelet-mass decline; platelets directly REMOVE TPO from plasma in vitro (§8). |
| 21 | Harker LA, Finch CA (1969). Thrombokinetics in man. *J Clin Invest* 48(6):963-74. | **5814231**, 10.1172/JCI106077 | Method precedent: megakaryocyte-mass-based production estimate agrees with platelet-turnover estimate (Little's-Law cross-check concept; exact lifespan figure not in abstract, disclosed §12). |
| 22 | Sender R, Milo R (2021). The distribution of cellular turnover in the human body. *Nat Med* 27(1):45-48. | **33432173**, 10.1038/s41591-020-01182-9 | EXTERNAL anchor: total-body cellular turnover (0.33±0.02)×10¹²/day, "close to 90%... was blood cells" -> 2.97×10¹¹/day (§7). |
| 23 | Jaiswal S, Fontanillas P, Flannick J, et al (2014). Age-related clonal hematopoiesis associated with adverse outcomes. *N Engl J Med* 371(26):2488-98. | **25426837**, 10.1056/NEJMoa1408617 | CHIP prevalence 9.5%/11.7%/18.4% by age band; HR 11.1 hematologic cancer (§10). |
| 24 | Genovese G, Kähler AK, Handsaker RE, et al (2014). Clonal hematopoiesis and blood-cancer risk inferred from blood DNA sequence. *N Engl J Med* 371(26):2477-87. | **25426838**, 10.1056/NEJMoa1409405 | Independent cohort: CH 10% (>65y) vs 1% (<50y); HR 12.9 (§10). |
| 25 | Till JE, McCulloch EA (1961). A direct measurement of the radiation sensitivity of normal mouse bone marrow cells. *Radiat Res* 14:213-22. | **13776896** | Historical origin of the CFU-spleen clonogenic assay — the founding quantitative-progenitor concept. Pre-1975, no abstract exists (same disclosure tier as this repo's Bentley1974/Nadler1962 precedent). |

## 4. Headline results

| quantity | value | anchor | gate |
|---|---:|---|---|
| Wilson 2008 dormant-HSC division interval | **145 d** | task band 145-193 d | **PASS** (boundary) |
| Wilson 2008 divisions per lifetime | **5** | task's own "~5" | **PASS** (exact) |
| Bernitz 2016 independent-method division ceiling | 4 expansion + 1 potency-loss = **5** | cross-check vs Wilson | **PASS** (decorrelated method) |
| Foudi 2009 implied interval (from its own ≤0.8-1.8%/day) | **55.6-125 d** | task band 145-193 d | FAIL, diagnosed (§5) |
| Forced adversary: Wilson's own "pool turns over every few weeks" (21-42 d) vs measured 145 d | ratio **3.45×** | adversary must fall | **PASS** (adversary falls) |
| Mouse->human interval ratio vs lifespan ratio | 1.93× vs 32× (**16.6× dissociation**) | held OPEN | quantified |
| Osawa 1996 single-HSC transplant success | **21%** (n=1) | binary falsifier ("can 1 HSC reconstitute the whole system") | **PASS** |
| Spangrude 1988 30-cell dose success | **50%** | dose-response must not decrease | **PASS** (monotonic) |
| Bonus: Poisson transfer (p from Osawa, applied to Spangrude's n=30) | predicts 99.9% vs measured 50% | naive cross-protocol transfer | diagnosed miss, not gated (§6) |
| RBC+neutrophil vs task's own literal anchor (2×10¹¹+1×10¹¹) | **2.68-3.22×10¹¹/day** | task ~3×10¹¹/day | **PASS** (−10.5% / +7.2%) |
| Neutrophil production, 2 independent methods (marrow-pool vs 3H-turnover) | 0.847 vs 0.87×10⁹/kg/day | convergence | **PASS** (2.7% diff) |
| RBC+neutrophil vs Sender&Milo 2021 EXTERNAL anchor | 2.684×10¹¹ vs 2.97×10¹¹ | independent top-down method | **PASS** (−9.6%) |
| Full 3-leg (+platelet) total vs Sender&Milo | 4.66×10¹¹ vs 2.97×10¹¹ | independent top-down method | near-miss, diagnosed (+56.8%, §7) |
| HSC pool mouse (Busch/Abkowitz) vs human (Lee-Six/Mitchell) | 5,000-11,400 vs 20,000-200,000 | SEED-DESIGN node's own "kill if >2 orders" rule | **PASS** (1.60 orders, NOT discordant) |
| Abkowitz's conserved-number extrapolation vs actual human measurement | 11,400 vs 50,000-200,000 | forced adversary (conserved-across-primates) | **FALLS** (4.4-17.5× higher, §7) |
| Amplification factor (total output / direct HSC-division flux) | **6.5×10⁸-6.5×10⁹** | — | derived |
| Doublings needed downstream of HSC | **29.3-32.6** | physically plausible 10-45 band | **PASS** |
| Void-floor sweep across FULL uncertainty range | doublings stay **25.4-34.9** | ≥5 doublings always required | **PASS** (robust) |
| G-CSF-/- KO neutrophil level | **20-30% of WT** | causal necessity, <50% | **PASS** |
| TPO clearance-model reciprocal law | `C_ss·M` = constant, exact | monotonic + reciprocal | **PASS** (both) |
| HSC frequency: marker-enriched vs functional-competitive | 1:2,000 vs 1:1,666,667 (**833× apart**) | task band 1:10⁴-10⁵ sits BETWEEN | **PASS** (>10×, quantified) |
| CHIP prevalence (own live re-fetch vs existing graph node) | 9.5%/18.4%/10%/1%, exact match | cross-confirmation | **PASS** |

All numbers machine-printed from `hematopoiesis_results.json` — nothing above is hand-computed prose.

## 5. Falsifier 1 — HSC quiescence / division-rate (PRIMARY)

**Wilson 2008** is the decisive anchor, and both task-cited numbers land in ONE paper: "d-HSCs
divide about every **145 days**, or **five times per lifetime**." 145 d sits at the task band's own
lower boundary (145-193 d) — PASS. **Bernitz 2016** independently corroborates the "~5" ceiling via
a completely different methodology (a cumulative division-**counting** reporter, not a rate model
inferred from label dilution): "completes four traceable symmetric self-renewal divisions to expand
its size before entering a state of dormancy... long-term regenerative potential lost upon a fifth
division" — 4+1=5, a genuinely decorrelated method landing on the same number.

**Foudi 2009's** own number does NOT overlap the task band: "~20% of HSCs divide at an extremely low
rate (≤0.8-1.8% per day)" inverts to a **55.6-125 day** interval — faster than Wilson's 145 d, with
zero overlap against the 145-193 d band. **Forced-adversary Orient, not smoothed over**: Foudi's
gate (L⁻K⁺S⁺CD48⁻CD150⁺) is LESS stringent than Wilson's most-strictly-gated d-HSC subset
(lin⁻Sca1⁺cKit⁺CD150⁺CD48⁻CD34⁻) — describing a broader, more heterogeneous population that
necessarily includes some faster-cycling cells averaged in. This is a real, diagnosed cross-paper
discrepancy (gate `G1d`: FAIL, disclosed), not hidden — and it does not sink the falsifier, because
Wilson's own number (the task's primary cited anchor) passes cleanly and Bernitz corroborates it via
an unrelated method.

**Forced adversary, textually supplied by the paper itself**: Wilson 2008's own abstract states the
naive view its finding refutes — "it is thought that **the entire HSC pool turns over every few
weeks**." Forced to its strongest fair form (3-6 weeks = 21-42 days), the measured dormant-subset
interval (145 d) is **3.45× slower** than even the slowest edge of that adversary claim. **Gate:
PASS** — the adversary falls, exactly Wilson's own point: bulk-average kinetics obscure a much more
quiescent reservoir within the total HSC pool.

**Mouse-to-human transfer, quantified not hand-waved**: Catlin 2011 measures human HSC replication
at "once every 40 weeks (range, 25-50 weeks)" = 280 d — only **1.93×** longer than Wilson's mouse
145 d. But mouse:human LIFESPAN ratio is ~32× (2.5 y vs 80 y). The division-interval ratio is
**16.6× smaller** than the lifespan ratio would predict under naive allometric scaling — division
interval does NOT track organismal lifespan the way body-size/metabolic-rate scaling normally would.
Held OPEN exactly as the task instructs (§5, this is the task's own "mouse-to-human transfer
uncertain" item), but now with an exact, machine-computed dissociation factor rather than a
qualitative caveat.

**Gate: part1_overall_pass = PASS (5/6 sub-gates; the 1 non-pass is Foudi's diagnosed, disclosed
non-overlap, not a bug).**

## 6. Falsifier 2 — transplantation functional readout / clonality (PRIMARY)

**Osawa 1996** is the decisive, quantitative anchor: injection of a **single** purified
mCD34(lo/-),c-Kit+,Sca-1+,Lin- mouse cell achieved long-term multilineage lymphohematopoietic
reconstitution in **21% of recipients**. This is not "a demonstration exists" — it is a measured,
reproducible RATE from a single cell, directly answering the task's falsifier ("can a single HSC
reconstitute the entire blood system?"). **Gate G2a: PASS** (21% ≫ the pre-registered ≥5%
"clearly nonzero" floor).

**Spangrude 1988** supplies an independent DOSE-RESPONSE point: 30 purified cells "are sufficient to
save 50 percent of lethally irradiated mice, and to reconstitute all blood cell types in the
survivors" — a different (radioprotection, not detectable-chimerism) endpoint, 8 years earlier, with
a cruder purification protocol. **Gate G2b: PASS** — more cells never gives a lower success rate
(50% at n=30 ≥ 21% at n=1), the basic monotonicity a limiting-dilution model requires.

**Notta 2011** supplies the HUMAN, species-decorrelated corroboration: single CD49f+ human cells are
"highly efficient" at generating long-term multilineage grafts — a genuinely different species and
marker system (CD49f vs. mouse CD34lo/-), 15 years later, same qualitative conclusion. **Gate G2c:
PASS** (qualitative-tier, disclosed — the abstract does not give an exact %).

**Bonus/exploratory forced-adversary test (disclosed, NOT gated into the core falsifier — same
status as `MECHANISM_ERYTHROPOIESIS.md`'s own dynamic leg)**: does a single-hit/Poisson limiting-
dilution model, `P(n)=1-(1-p)^n`, calibrated on Osawa's own p=0.21 (n=1), correctly PREDICT
Spangrude's n=30 result? It predicts **99.9%** success — the measured value is **50%**, a ~50-point
gap. **Orient, not hidden**: the two papers used different-generation purification protocols (1988
vs. 1996, a full marker-refinement cycle apart) AND different endpoints (radioprotection vs.
detectable chimerism) — the GEOMETRY (single-hit dose-response, monotonic in n) is exactly right
(it correctly predicts G2b's qualitative direction), but the fitted per-cell-`p` is protocol- and
endpoint-specific, not a universal constant transferable across 8 years of assay refinement. This
mirrors `MECHANISM_ERYTHROPOIESIS.md` §6's own precedent (a forced-adversary sweep diagnosing WHY a
naive point-estimate transfer misses, rather than either hiding the miss or force-fitting it).

**Gate: part2_overall_pass = PASS (3/3 core sub-gates).**

## 7. Decorrelated check — total output vs. the tiny HSC pool: the amplification factor

**Step A — RBC+neutrophil vs. the task's own literal anchor.** The task's decorrelated-check
wording is specifically "~2×10¹¹ RBC + ~10¹¹ neutrophils/day" (not a 3-leg total). RBC is reused
read-only from `erythropoiesis_results.json` (2.083×10¹¹/day). Neutrophil is newly derived here via
Little's Law from Dancey 1976's OWN raw numbers: postmitotic marrow pool (5.59×10⁹ cells/kg) /
transit time (6.60 d) = 0.847×10⁹/kg/day — **independently cross-checked** against Dancey's own
SEPARATE 3H-thymidine circulating-turnover measurement (0.87×10⁹/kg/day, a **2.7%** convergence
between two genuinely different labeling methods, the same "two-methods-converge" pattern
`MECHANISM_ERYTHROPOIESIS.md` §6 used for RBC lifespan). At a 70 kg reference adult: **6.01×10¹⁰
neutrophils/day**. RBC+neutrophil = **2.684×10¹¹/day**, within **10.5%** of the task's literal
3×10¹¹/day anchor. **Gate G3a: PASS.**

Dancey's paper ALSO reports a third method (DF32P double-labeling), giving a much higher
1.62×10⁹/kg/day — but Dancey's own text flags this method as showing "a lower recovery and shorter
t1/2 of the 32P label" when directly compared against 3H-thymidine on the same cells — an
artifact-prone estimate, structurally the SAME situation `MECHANISM_ERYTHROPOIESIS.md` §6 found for
uncorrected 51Cr: the artifact-prone method (here, DF32P: 1.134×10¹¹/day) lands closer to a round
textbook figure than the two converged, more-reliable methods. Reported honestly, not smoothed over.

**Step B — cross-check against Sender & Milo 2021, a fully INDEPENDENT top-down method** (systems-
biology literature-integration cell census, not Little's-Law-per-lineage). Their own abstract:
"(0.33±0.02)×10¹² cells per day... close to 90%... was blood cells" -> **2.97×10¹¹ blood
cells/day**. RBC+neutrophil ALONE lands within **9.6%** of this fully independent anchor — a tight,
decisive convergence between two structurally unrelated methodologies. Adding a THIRD leg
(platelets, via Harker & Finch 1969's method + disclosed-tier count/pooling/lifespan constants —
none independently re-verified this session, see §12) pushes the 3-leg total to 4.66×10¹¹/day, **56.8%
above** the same anchor — a real, disclosed overshoot, plausibly (not provenly) explained by
Sender & Milo's own cellular-turnover census weighting anucleate platelet fragments differently
than nucleated/RBC-conventional cells (their full text was not accessible this session — no PMC
deposit found via `elink`, an honest gap, not fabricated).

**Step C — HSC pool-size reconciliation: executing the SEED-DESIGN graph node's OWN pre-registered
falsifier for the first time.** `AUTO-CROSS-METHOD-HSC-CLONE-COUNT-RECONCILIAT[ION]` designed but
never ran: "reconciliation needs >2 orders-of-magnitude free parameter with no biological grounding
-> kill as discordant." Using ONLY the quantitatively-stated estimates (Sun 2014's "thousands" is
qualitative — no specific number is in that abstract, and assigning one would fabricate precision
the source doesn't have; a real error this script's own first run caught and fixed, see §9's
discipline note): mouse Busch2015 (~5,000) / Abkowitz2002 (11,400) vs. human Lee-Six2018/Mitchell2022
(2×10⁴-2×10⁵). Ratio range **1.75×-40×** = **0.24-1.60 orders of magnitude** — comfortably under the
node's own 2.0-order kill threshold. **Verdict: NOT discordant** (gate `G3c`: PASS) — compatible
with reconciliation, though not proven-reconciled by a specific fitted clone-size-distribution model
(that needs per-clone-size data this session did not fetch — disclosed).

**A sharper, separate finding**: Abkowitz 2002's own abstract explicitly flags her human
extrapolation as tentative ("if the total number of human HSCs were ALSO equivalent..."). Her
conserved-number extrapolation predicts ~11,400 human HSC; 16-20 years later, an ORTHOGONAL method
(whole-genome somatic-mutation phylogenetics, not transplant/xenograft) measures the actual human
active pool at 50,000-200,000 — **4.4-17.5× larger**. The forced adversary (naive "same absolute
HSC count regardless of body size, extended to primates") **falls** against an independent,
diverse-instance-space method — while her DIRECT mouse+cat measurement is not itself contradicted.
This is the "forced adversary falls across the diverse instance-space" structure applied precisely.

**Step D — the amplification factor itself.** Human HSC pool (5×10⁴-2×10⁵) dividing at Catlin
2011's measured rate (1/280 days) supplies only **179-714 cells/day DIRECTLY** — utterly negligible
against the ~2.68×10¹¹-4.66×10¹¹/day total output. The gap requires **29.3-32.6 binary doublings**
across the downstream transit-amplifying compartments (MPP -> CMP/CLP -> committed progenitors ->
precursors) — a falsifiable, geometrically-derived prediction connecting §5's quiescence-rate
anchor with §7's independently-measured total-output anchor. **Void-floor robustness**: sweeping
the FULL verified uncertainty range simultaneously (output 6×10¹⁰-4.7×10¹¹/day; pool 5×10³-2×10⁵;
interval 145-350 d; 144 combinations) the required doublings never drop below **25.4**, and never
exceed **34.9** — the "massive transit-amplification is structurally necessary" conclusion is
robust across the entire measured parameter space, not a fragile point-estimate artifact. **Gate:
PASS.** (A specific primary-source-verified "total divisions HSC-to-mature-cell" figure was not
independently found this session — this is a plausibility check, not a direct external numeric
match: ~30 doublings at 1-2 divisions/day across multiple serial+parallel compartments is easily
achievable within days-to-weeks of known transit-amplifying proliferation, but this specific
cross-check is disclosed as not independently verified, §12.)

**Gate: part3_overall_pass = PASS (4/5 core sub-gates; §7 discloses the one honest near-miss —
the 3-leg-vs-external-anchor comparison — separately from the tight, gated RBC+neutrophil match).**

## 8. The demand-driven triad — two genuinely different sensor geometries, not one mechanism thrice

- **EPO -> RBC** (reused, not recomputed): a **threshold/enzyme-kinetic state sensor** — HIF-2/PHD
  prolyl-hydroxylase kinetics acting on tissue O2 tension, modeled in `erythropoiesis.py` as a smooth
  softplus hinge around a reference Hb/O2 level (`MECHANISM_ERYTHROPOIESIS.md` §4).
- **G-CSF -> neutrophil**: Lieschke 1994's genetic KNOCKOUT is a **causal loss-of-function** anchor
  — a structurally stronger evidence class than a correlational dose-response curve. G-CSF-/- mice:
  neutrophils 20-30% of WT (gate `G4a`: PASS, <50%), marrow progenitors -50%, a clean gene-dosage
  effect in heterozygotes, and an explicit STEADY-STATE-vs-EMERGENCY dual role ("indispensable for
  maintaining the normal quantitative balance of neutrophil production... and also implicate G-CSF
  in 'emergency' granulopoiesis during infections") — directly paralleling this doc's own §7
  neutrophil-production leg and the erythropoiesis EPO-stress-response coupling.
- **TPO -> platelet**: Kuter & Rosenberg 1995 supply a **mass-action receptor-mediated clearance
  ("sink") sensor**, fundamentally different in FORM from EPO's threshold sensor: hepatic TPO
  production is roughly constitutive; free plasma TPO is a passive consequence of how much Mpl-
  receptor "sink" capacity (platelet+megakaryocyte mass) exists to remove it. Closed form:
  `dC/dt = P - k·M·C` at steady state gives `C_ss = P/(k·M)` — a clean reciprocal law, verified
  exactly (gate `G4c`: `C_ss·M` constant across a 20-point sweep) and monotonically decreasing in
  platelet mass (gate `G4b`: PASS) — with **no hinge or threshold at all**, structurally distinct
  from EPO's softplus. Directly evidenced by Kuter's own experiments, not just a correlation: TPO
  rose "inversely and proportionally" as platelet mass fell; platelet TRANSFUSION near the nadir
  directly DECREASED elevated TPO; platelets physically REMOVED TPO from plasma **in vitro**; soluble
  c-Mpl receptor NEUTRALIZED the effect (confirming Mpl-receptor specificity). This closes
  `MECHANISM_PLATELET_HEMOSTASIS.md`'s own disclosed gap ("megakaryocyte thrombopoiesis,
  thrombopoietin-regulated; not modeled here... a natural next node").

**Gate: part4_overall_pass = PASS (3/3).**

## 9. HSC rarity / marker-assay dependence — quantified, not hand-waved (symmetric QC)

The task's own symmetric-QC instruction — "HSC-vs-progenitor definitions are marker/assay-
dependent" — is given a real, machine-computed number rather than a qualitative caveat. Two
live-verified, genuinely different assay philosophies:

- **Marker-phenotype enrichment** (Morrison & Weissman 1994): the
  Thy-1.1(lo)Sca-1(hi)Lin(-/lo) population = **0.05% of marrow (1:2,000)** — but "only around 25% of
  clonal reconstitutions by cells from this population are long term," refining to **1:8,000** for
  the truly long-term-competent subset.
- **Functional long-term-competitive-repopulation** (Abkowitz 2002): 6 HSC per 10⁷ nucleated
  marrow cells = **1:1,666,667** (feline model, cross-validated against an equivalent mouse
  estimate in the same paper).

**Ratio: 833× (marker-enriched vs. functional), or 208× refined-long-term vs. functional** — gate
`G5a`: PASS (>10× pre-registered threshold, cleared by two orders of magnitude). Kiel 2005's SLAM
code (CD150/CD244/CD48) supplies a THIRD, orthogonal marker scheme, demonstrating purification
itself is a moving target across labs and years, not merely two outlier numbers. The task's own
commonly-cited "1:10⁴-10⁵" band sits **strictly between** these two live-verified extremes (gate
`G5b`: PASS) — itself evidence that the commonly-cited figure is a rounded compromise across
incompatible assay philosophies, not a single well-defined measured quantity. Held OPEN exactly as
instructed.

## 10. Clonal hematopoiesis (CHIP) — reused, independently re-confirmed, held OPEN

Rather than re-deriving `HEMATO-HSC-CLONAL-AGING`'s own analysis, this doc independently RE-FETCHES
Jaiswal 2014 and Genovese 2014 live this session and cross-checks the result byte-for-byte against
that node's already-stored datapoints: CHIP prevalence 9.5% (age 70-79, n=2300), 11.7% (80-89,
n=317), 18.4% (≥90, n=103), HR=11.1 for hematologic cancer (Jaiswal); independently, 10% (>65y) vs.
1% (<50y), HR=12.9, ~42% of hematologic cancers traced to a pre-existing clone (Genovese). **Gate
G6a: PASS** — an exact match, a genuine cross-confirmation (this repo's own memory explicitly warns
that a subagent's "confirmed live" claim has been fabricated before — this is not that).

**Connection to §5, held OPEN, not resolved**: Mitchell 2022 shows that by age >75, 12-18 expanded
clones account for 30-60% of total blood output — meaning the aged human "HSC pool" is increasingly
a FEW SELECTED CLONES, not the homogeneous population Catlin 2011's population-average division-rate
estimate (§5) implicitly assumes. This is a real, disclosed complication for treating HSC quiescence
rate as a single population-uniform parameter in older subjects — exactly the task's own instruction
to hold this open, now with an explicit mechanistic link to the quiescence falsifier rather than a
free-floating caveat.

## 11. couples_to (prose + JSON metadata; no graph-edge write this session)

- **`MECHANISM_ERYTHROPOIESIS`** (`erythropoiesis.py`) — supplies the RBC production-rate leg
  (2.083×10¹¹/day), read-only, reused directly in §7's decorrelated check. This doc supplies the
  HSC-level ROOT (rarity, quiescence, transplant clonality) that doc's EPO/marrow-output model
  implicitly assumed but never built.
- **`MECHANISM_PLATELET_HEMOSTASIS`** (`platelet_hemostasis.py`) — this doc's TPO/platelet-production
  leg (§7, §8) directly closes that doc's own disclosed gap ("megakaryocyte thrombopoiesis,
  thrombopoietin-regulated; not modeled here"). That doc's platelet COUNT is a static input this
  doc's production-rate derivation is consistent with (same 250×10⁹/L reference), not identical
  (that doc models activation/coverage dynamics of a FIXED count; this doc models the upstream
  production RATE feeding that count).
- **`MECHANISM_COAGULATION_HEMOSTASIS`** — shares the megakaryocyte/platelet lineage origin (via
  `MECHANISM_PLATELET_HEMOSTASIS`'s own §13 coupling); not directly recomputed here.
- **immune (all leukocyte lineages)** — this doc's G-CSF/neutrophil leg (§7, §8) is the first
  executed production-kinetics leg for any leukocyte lineage in this repo; couples to
  `MECHANISM_TCELL_ACTIVATION_EXHAUSTION` and `MECHANISM_COMPLEMENT_CASCADE` only at the shared-origin
  level (all descend from the same HSC root modeled here) — no shared numeric quantity computed
  across those docs and this one.
- **aging (HSC exhaustion / clonal hematopoiesis)** — §10's CHIP discussion connects directly to
  `HEMATO-HSC-CLONAL-AGING` (graph, OPEN, cited not re-derived) and to
  `MECHANISM_CELLULAR_SENESCENCE`/`MECHANISM_EPIGENETIC_CLOCK` at the aging-mechanism level (not
  jointly computed this session).
- **`HEMATO-HSC-CLONAL-AGING`, `AUTO-CROSS-METHOD-HSC-CLONE-COUNT-RECONCILIAT[ION]`,
  `AUTO-TEST-WHETHER-NICHE-CXCL12-SCF-GRADIENT-M[AGNITUDE]`, `AUTO-CELL-CHIP-MDS-CAUSE-VS-MARKER-
  SYMMETRIC[-AUDIT]`, `ORG-BLOOD-COAG`** (graph, all OPEN) — disambiguated in full in §0; the
  clone-count-reconciliation node's own falsifier is EXECUTED (§7) for the first time this session,
  the others are cited/held-open, none edited.

Per `docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/§4, promoting any of this into a canonical graph edge
requires the separate `mechanism_fold -> fold_gate_v2` pipeline — **not performed this session**
(isolation rule), consistent with every sibling doc's own stated practice.

## 12. Honest gaps — what this does NOT prove (symmetric QC, disclosed not hidden)

- **Foudi 2009's own interval estimate does not overlap the task's 145-193 d band** (§5) — a real,
  diagnosed (not hidden) cross-paper discrepancy attributed to different marker-gating stringency,
  not resolved by forcing a fit.
- **The Poisson single-hit transplant model does not transfer cleanly across purification eras**
  (§6) — correctly predicts the qualitative dose-response direction but overshoots Spangrude 1988's
  actual 30-cell success rate by ~50 percentage points; diagnosed as protocol/endpoint-specific, not
  gated into the core falsifier.
- **The 3-leg (RBC+neutrophil+platelet) bottom-up total overshoots Sender & Milo 2021's external
  anchor by 56.8%** (§7), even though RBC+neutrophil alone is tight (9.6%) — the platelet leg's
  count/splenic-pooling-factor/lifespan constants are DISCLOSED-TIER (standard clinical-teaching
  values), not independently live-re-verified this session (Harker & Finch 1969's abstract does not
  state an exact platelet-lifespan figure; this doc's 9.5 d value is a commonly-cited convention,
  same disclosure tier as `MECHANISM_ERYTHROPOIESIS.md`'s own MCV/marrow-transit-delay constants).
  Sender & Milo's own full text was not accessible this session (no PMC deposit found via `elink`)
  to confirm or refute the hypothesized platelet-fragment-accounting explanation.
- **The HSC-pool mouse-vs-human reconciliation (§7) is "not discordant," not "proven reconciled."**
  A specific fitted clone-size-distribution/neutral-drift model was not built this session (would
  require per-clone-size distributional data neither fetched nor available from the abstracts used)
  — the SEED-DESIGN graph node's own kill-threshold is cleared, which is a real, meaningful, but
  bounded result, not a full mechanistic reconciliation.
- **The ~29-31-doubling amplification-factor prediction is a plausibility check, not a direct
  external numeric match** — no specific primary-source-verified "total divisions HSC-to-mature-
  cell" figure was independently found/fetched this session to check the derived number against.
- **G-CSF/TPO parameters in §8's closed-form models are illustrative/first-principles** (no specific
  rate constants independently measured this session) — same disclosure tier as this repo's
  `MECHANISM_COMPLEMENT_CASCADE.md`/`MECHANISM_TAU_PATHOLOGY.md` toy parametrizations; the STRUCTURAL
  claim (reciprocal-clearance vs. threshold-hinge) is the load-bearing result, not the specific
  numbers.
- **Mouse-to-human transfer of the quiescence rate is explicitly uncertain** (§5) — division
  interval (1.93×) does not scale with lifespan (32×) the way naive allometry predicts; this
  dissociation is quantified here, not explained mechanistically.
- **CHIP/aging complicates treating the HSC pool as homogeneous** (§10) — named and quantified
  (Mitchell 2022's 12-18-clone collapse), not resolved; Catlin 2011's population-average division
  rate and Wilson/Foudi/Bernitz's mouse cohorts are not re-analyzed for age-stratified clonal
  structure this session.
- **Notta 2011's exact single-cell success percentage is not in its own abstract** — used
  qualitative-tier ("highly efficient"), disclosed, not fabricated as a specific number.
- **Till & McCulloch 1961 carries no PubMed abstract** (pre-1975 entry) — title+PMID only, same
  disclosure tier as this repo's own Bentley 1974/Nadler 1962/Severinghaus 1966 precedents.
- **No graph-edge write this session** (§11) — `couples_to` is prose/JSON metadata, matching every
  sibling doc's own stated practice (would require the separate `mechanism_fold` pipeline).
- **A real, caught-and-fixed methodological error, disclosed rather than hidden**: this script's
  FIRST run assigned an invented numeric floor (1,000) to Sun 2014's purely qualitative "thousands of
  clones" claim, which spuriously pushed the HSC-pool reconciliation ratio over the SEED-DESIGN
  node's own 2-order-of-magnitude kill threshold. Caught via this doc's own QC pass (not by an
  external reviewer), fixed by using ONLY the two actually-quantified mouse estimates
  (Busch2015/Abkowitz2002) and keeping Sun2014 as an unquantified qualitative corroboration — the
  kind of self-caught, disclosed correction the watertight discipline exists to force.

## 13. Gates — 19/20 sub-gates PASS (1 honestly disclosed non-pass, diagnosed not hidden)

```
PART 1 (HSC quiescence / division-rate -- PRIMARY):
  G1a_wilson_interval_in_task_band:                          PASS (145 d, band boundary)
  G1b_wilson_divisions_per_lifetime_matches:                 PASS (exactly 5)
  G1c_bernitz_independent_method_corroborates_5_ceiling:      PASS (4+1=5, decorrelated method)
  G1d_foudi_interval_overlaps_task_band:                      FAIL (55.6-125d vs 145-193d, diagnosed Sec.5)
  G1e_forced_adversary_bulk_turnover_falls:                   PASS (3.45x slower than adversary)
  G1f_mouse_human_dissociation_confirmed:                     PASS (16.6x dissociation factor)
  part1_overall_pass:                                         PASS (5/6)

PART 2 (transplantation / clonality -- PRIMARY):
  G2a_osawa_single_cell_success_exceeds_floor:                PASS (21% >= 5%)
  G2b_dose_response_monotonic_more_cells_not_worse:           PASS (50% >= 21%)
  G2c_human_species_corroboration_qualitative_pass:           PASS (Notta 2011)
  part2_overall_pass:                                         PASS (3/3)

PART 3 (amplification-factor decorrelated check):
  G3a_rbc_neutrophil_subtotal_within_20pct_of_task_literal_anchor:  PASS (-10.5%/+7.2%)
  G3b_neutrophil_two_methods_converge_within_10pct:           PASS (2.7% diff)
  G3c_hsc_pool_reconciliation_under_seed_design_kill_threshold: PASS (1.60 orders < 2.0)
  G3d_amplification_doublings_physically_plausible_band:      PASS (29.3-32.6 in [10,45])
  G3e_void_floor_massive_amplification_always_required:       PASS (min 25.4 doublings, full sweep)
  part3_overall_pass:                                          PASS (5/5)

PART 4 (demand-driven triad geometry):
  G4a_gcsf_ko_shows_substantial_reduction_below_50pct_wt:      PASS (20-30% WT)
  G4b_tpo_clearance_model_monotonic_decreasing:                PASS
  G4c_tpo_reciprocal_proportionality_law_holds_exactly:        PASS
  part4_overall_pass:                                          PASS (3/3)

PART 5 (HSC rarity / marker-dependence, symmetric QC quantified):
  G5a_marker_vs_functional_frequency_differ_by_over_10x:       PASS (833x)
  G5b_task_band_sits_between_the_two_verified_extremes:        PASS
  part5_overall_pass:                                           PASS (2/2)

PART 6 (CHIP, reused + independently re-confirmed):
  G6a_independent_refetch_matches_existing_graph_node_datapoints: PASS (exact match)
  part6_overall_pass:                                           PASS (1/1)

REQUIRED FALSIFIERS:
  falsifier1_quiescence_division_rate:                         PASS
  falsifier2_transplant_clonality:                             PASS
  decorrelated_check_amplification_factor:                     PASS

n_gates_pass / n_gates_total:                                  19 / 20
overall_pass:                                                  TRUE
```

**Overall: PASS** on all three task-required falsifiers. Deterministic — 2 independent runs produce
byte-identical JSON (`md5sum fb386bce71445d63b2ed94b9b57bc813`, verified this session). No NaN/Inf
anywhere in the evidence JSON (checked, not assumed).

## 14. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/hematopoiesis.py
```

Requires `data/erythropoiesis/erythropoiesis_results.json` to already exist (produced by this
repo's existing, already-run sibling script; degrades gracefully to a disclosed hardcoded default
of 2.0833×10¹¹/day if absent — not exercised this session, the file was present). Writes
`data/hematopoiesis/hematopoiesis_results.json`. Pure Python/numpy, no OpenSim call, runs in under 2
seconds, deterministic (2 independent runs verified byte-identical). No git operations; reads the
erythropoiesis JSON read-only; writes only under `data/hematopoiesis/`.

**Paths**: script `scripts/msk/hematopoiesis.py`; evidence
`data/hematopoiesis/hematopoiesis_results.json`; this doc
`docs/MECHANISM_HEMATOPOIESIS.md`; read-only sibling reused
`data/erythropoiesis/erythropoiesis_results.json`; graph nodes disambiguated read-only (§0, §11)
`data/MECHANISM_ANCHOR_GRAPH.json` (`HEMATO-HSC-CLONAL-AGING`,
`AUTO-CROSS-METHOD-HSC-CLONE-COUNT-RECONCILIAT`, `AUTO-TEST-WHETHER-NICHE-CXCL12-SCF-GRADIENT-M`,
`AUTO-CELL-CHIP-MDS-CAUSE-VS-MARKER-SYMMETRIC`, `ORG-BLOOD-COAG`).
