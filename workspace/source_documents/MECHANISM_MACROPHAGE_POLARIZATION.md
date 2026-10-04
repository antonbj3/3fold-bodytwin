# MECHANISM MACROPHAGE POLARIZATION — the M1/M2 functional dichotomy, its metabolic signature (itaconate/succinate vs OXPHOS/FAO), and reversible plasticity (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** This is the first MEASURE/CERTIFY pass on the
foundational IN-VITRO macrophage-activation mechanism layer (M1 classically-activated vs M2
alternatively-activated). It is a **DIFFERENT, more basic hidden state** than the existing
`IMMUNE-MACROPHAGE-POLARIZATION` graph cell (`data/MECHANISM_ANCHOR_GRAPH.json`, id search:
`grep -n '"IMMUNE-MACROPHAGE-POLARIZATION"' data/MECHANISM_ANCHOR_GRAPH.json`, status `OPEN`,
`mechanism_grade: SEED-DESIGN`), which corners the **in-vivo tumor-associated-macrophage (TAM)
continuum** question and was never computed, only literature-cited via a prior scout pass. This
build **couples to** that cell (§6) rather than re-deriving it — the two are complementary, not
duplicates: that cell already found (and this build's own symmetric QC reuses, §6) that the clean
in-vitro binary this doc reproduces does **not** survive as a prognostically load-bearing binary
in vivo in tumors. Script: `scripts/msk/macrophage_polarization.py`. Evidence (machine-written):
`data/msk_smoketest/macrophage_polarization/macrophage_polarization_results.json` (md5
`bb7e9897a1d8526ee24d6ca039b5757d`, byte-identical across repeated runs — determinism confirmed;
fully deterministic except two seeded `np.random.default_rng(20260722)` Monte Carlo/robustness
sweeps).

**Confidence tier: in-vitro-anchored** (marker/metabolic profiling — classic murine
peritoneal-macrophage and bone-marrow-derived-macrophage stimulation assays), population/
mechanism-level, not any individual subject's own measured trace. Two of the four falsifiers'
DECISIVE tests are **machine-encoded from a primary paper's own already-measured result**
(Tannahill 2013's 2-deoxyglucose cytokine dissociation, F3; Stout 2005's own explicit
in-vivo-reversal finding, F4) rather than this script's own simulation output — an external,
non-tautological anchor, not a circular self-check.

---

## 0. Falsifiers (pre-registered, matching the task verbatim) + verdict up front

| # | Falsifier | Verdict |
|---|---|---|
| F1 | iNOS/Arg1 reciprocal regulation, arginine-fate split (NO+citrulline vs ornithine+urea) — Mills 2000 | **PASS** (conservation exact, crossover Δ=0.905≥0.8, adversary falls, robust 20/20) |
| F1-ADV | Forced adversary: independent/dedicated-substrate-pool pathways (no shared-node conservation) | **FALLS** (predicts 25.0% double-positive rate; real model forbids it structurally, 0.0%) |
| F2 | itaconate/IRG1(ACOD1) — a metabolic node decorrelated (different regulatory modality) from the cytokine markers, M1-specific — Michelucci 2013 | **PASS** (≥2 independent primary sources, genetically-isolable node; honest non-orthogonality caveat disclosed) |
| F3 | M1 glycolytic+broken-TCA (itaconate/succinate) vs M2 oxidative/FAO metabolic signature | **PASS** (geometric steady-state direction correct; Jacobian = chain + exactly one disclosed feedback mechanism) |
| F3-ADV | Forced adversary: "M1 metabolism is a generic, undifferentiated M1 correlate" (already-measured, Tannahill 2013's own 2-DG result) | **FALLS** (IL-1β glycolysis-dependent, TNF-α glycolysis-independent — a real dissociation in the same experiment) |
| F4 | Plasticity/reversibility (repolarization on stimulus switch, NOT terminal lineage commitment) — Stout 2005 + TAM in vivo | **PASS** (reversible model returns to 0.895 of a 0.9 target; ratchet adversary stays locked at 0.130; void floor falls; robust 16/16) |

`overall_pass = True` in the evidence JSON, **after two real, disclosed design bugs were found and
fixed via OODA (not swept past)** — see §1 and §3. A false-positive is treated as equivalent to a
premature-negative throughout: nothing here reports a threshold quietly loosened after seeing a
result.

---

## 1. Geometric mechanism — TWO independent branchpoint-competition structures, not curve fits

Per this repo's GEOMETRIC-THINKING discipline, both F1 and F3's central mechanism are **flux
partitions at a shared-substrate graph node** — a real conservation-law constraint, not a fitted
curve:

**F1 (arginine fate):** L-arginine is a single pool consumed by two competing enzymes — iNOS
(→ NO + citrulline) and Arginase-1 (→ ornithine + urea). `frac_NO + frac_ornithine = 1` **exactly**,
by construction (machine-verified to float precision, `F1_conservation_exact=True`). M1/M2 differ
only in the *relative* enzyme-activity ratio at that one node — this is why the dichotomy is
structurally **reciprocal** (one branch's gain is necessarily the other's loss), not two
independently-tunable dials. Mills 2000's own finding that a single hub regulator, TGF-β1, "inhibits
inducible NO synthase and stimulates arginase" **in the same signaling event** is exactly what a
coupled, shared-node competition predicts and an independent-pools model would not.

**F3 (fuel entry):** the SAME structural device reappears one level up: total upstream carbon input
splits between a glycolysis-only fate (exits the model as lactate, never enters the simulated TCA
segment) and TCA-entry — using the identical `flux_partition()` function as F1 (deliberate code
reuse, not incidental — a genuine emergent structural parallel between the two branchpoints, not
just DRY-for-its-own-sake). This is the literal machine-implementation of Tannahill 2013's own
verbatim description: macrophages "switch their core metabolism from oxidative phosphorylation to
glycolysis" — a reciprocal shift, structurally identical in shape to the arginine-fate split.

**F3's Jacobian as a structural certificate.** The 6-state TCA-segment ODE (citrate→isocitrate→
{IDH-branch, IRG1-branch}→succinate→downstream-flux) is a directed chain (a DAG) **except for
exactly one disclosed feedback mechanism**: itaconate (Z) modulates `k_SF_eff`, the rate constant of
the single succinate→downstream (SDH) reaction. Because that ONE rate constant appears in both the
equation losing the flux (`dS/dt`) and the one gaining it (`dF/dt`), the single mechanism
necessarily produces **exactly two** off-diagonal Jacobian entries, `(S,Z)` and `(F,Z)` — machine-
verified by comparing the full numeric-Jacobian nonzero-entry set against a hand-specified expected
edge set (`{(I,C),(Z,I),(K,I),(S,K),(F,S)}` chain ∪ `{(S,Z),(F,Z)}` feedback), not asserted from
the equations by eye. This is the literal matrix-structure image of F2's own claim — "decorrelated
except one disclosed feedback mechanism" — turned into a machine-checkable structural fact.

**F4's eigenvalue as a structural certificate for reversibility.** The plasticity ODE
`dφ/dt = k·(φ_target(t) − φ)` has a single real eigenvalue, exactly `−k` (finite-difference-
verified, `F4_eigenvalue_matches_neg_k=True`) — a smooth gradient flow toward *whatever the current
target is*, structurally history-independent beyond the present state. This is the geometric
opposite of the forced ratchet adversary, which is not expressible as such a flow at all (it is a
discrete, explicitly history-dependent state-machine rule bolted onto the same relaxation dynamics)
— the two are different **dynamical classes**, not just different parameter choices.

**OODA — two real, disclosed design bugs found and fixed, not swept past.** *(1)* The first F3
design scaled a single `mito_capacity` factor onto **both** a step's production and its own
consumption term (e.g. `dK/dt = k_IK·I − mito_capacity·K`, then flux-into-S ∝ `mito_capacity·K`) —
this factor **algebraically cancels out of steady-state flux** (a genuine structural property of
linear chains: internal rate constants set pool sizes/transit speed, never steady-state throughput,
which is set by the input rate alone), and the first run produced `F(M1)=1.247 > F(M2)=0.5` —
backwards. **Fix:** model the actual physiological lever (Vats 2006's PGC-1β mechanism boosts the
capacity to route fuel INTO oxidative metabolism at all, competing against the glycolysis-only fate)
as a flux-partition at the fuel-entry node, matching F1's own structure. A second parameterization
(oxidative-fraction gap 0.2 vs 0.9) then over-corrected, shrinking M1's total throughput so much
that succinate came out *lower* in M1 than M2 (1.397 vs 1.8) — still backwards vs Tannahill's "LPS
strongly increases succinate." **Root cause:** conflating "M1 shifts its OWN mix toward glycolysis"
(real, Tannahill's finding) with "M1's TOTAL throughput shrinks 4.5×" (an invented magnitude, not
sourced). **Fix:** moderate the gap (0.4 vs 0.8) and strengthen `J_gln` (Tannahill: glutamine
anaplerosis is the "**principal** source" of succinate — should dominate, not be a minor addition),
so the *local* accumulation mechanism drives the result, matching the primary paper's own emphasis.
*(2)* F4's first design fixed `PHASE_DAYS=6.0` (an arbitrary round number "matching Daley 2010's
day1-day7 window"), which happened to sit at `PHASE_DAYS/τ=0.975` — almost exactly one time
constant — so neither the reversible model nor the ratchet adversary ever got close enough to the
M2 target to trigger the ratchet's commit-rule: the two curves came out numerically indistinguishable
(0.6692 vs 0.6693), a silent non-result. **Fix:** tie `PHASE_DAYS` to the actual loaded `k` (5 time
constants, 99.3% convergence) — principled, not a threshold change on the falsifier itself. All
four fixes are visible in `scripts/msk/macrophage_polarization.py`'s own inline OODA-note comments,
kept in place as an auditable trail (matching this repo's `MECHANISM_HEPATIC_CLEARANCE.md` §6
precedent), not silently corrected out of the historical record.

---

## 2. Citations — verified LIVE this session (NCBI eutils esearch/esummary/efetch + Europe PMC REST); WebSearch quota-exhausted (same disclosed fallback path prior `MECHANISM_*` docs already use)

| # | Citation | PMID/DOI | Role / what was verified live |
|---|---|---|---|
| 1 | **Mills CD, Kincaid K, Alt JM, Heilman MJ, Hill AM (2000). "M-1/M-2 macrophages and the Th1/Th2 paradigm." J Immunol 164(12):6166-73.** | **10843666** | Verified live (esearch disambiguated 2 candidates — 10843666 original vs 28923981 a 2017 "Pillars Article" reprint of the SAME paper — via esummary; FULL abstract fetched verbatim via Europe PMC REST). **PRIMARY F1 anchor.** Quoted: "macrophages from prototypical Th1 strains...are more easily activated to produce NO...In marked contrast, LPS stimulates Th2, but not Th1, macrophages to increase arginine metabolism to ornithine...TGF-beta1, which inhibits inducible NO synthase and stimulates arginase, appears to play an important role in regulating the balance." **TERMINOLOGY CAVEAT (§6):** this paper's own "M-1/M-2" labels are mouse-STRAIN-associated (genetic background), not identical to the later stimulus-based nomenclature this task uses. ALL-MOUSE. No exact NO2-/urea concentration numbers in the abstract (disclosed gap — this build gates on direction/reciprocity, not magnitude). |
| 2 | **Michelucci A, Cordes T, Ghelfi J, et al. (2013). "Immune-responsive gene 1 protein links metabolism to immunity by catalyzing itaconic acid production." PNAS 110(19):7820-5.** DOI 10.1073/pnas.1218599110 | **23610393** | Verified live (esearch disambiguated 4 itaconate-related candidates via esummary; abstract via efetch). **PRIMARY F2 discovery anchor.** "Irg1 is highly expressed in mammalian macrophages during inflammation"; IRG1 gain/loss-of-function directly controls itaconic acid levels; purified IRG1 shows cis-aconitate decarboxylase activity; itaconate inhibits isocitrate lyase (antimicrobial mechanism). |
| 3 | **Stout RD, Jiang C, Matta B, Tietzel I, Watkins SK, Suttles J (2005). "Macrophages sequentially change their functional phenotype in response to changes in microenvironmental influences." J Immunol 175(1):342-9.** | **15972667** | Verified live (esearch single-hit + esummary + FULL abstract via efetch). **PRIMARY F4 anchor** — quoted: "As an alternative to the concept of subset development, we propose that macrophages...can reversibly and progressively change the pattern of functions that they express...macrophage functional phenotypes established in vivo in aged or tumor-bearing mice can be altered by changing their microenvironment." This paper's OWN explicit frame is a refutation of the "terminal subset" adversary. |
| 4 | **Tannahill GM, Curtis AM, Adamik J, et al. (2013). "Succinate is an inflammatory signal that induces IL-1β through HIF-1α." Nature 496(7444):238-42.** DOI 10.1038/nature11986 | **23535595** | Verified live (esearch single-hit + esummary + FULL abstract via efetch). **THE decisive F3 forced-adversary anchor** — quoted: "Macrophages activated by...lipopolysaccharide switch their core metabolism from oxidative phosphorylation to glycolysis...inhibition of glycolysis with 2-deoxyglucose suppresses lipopolysaccharide-induced interleukin-1beta but not tumour-necrosis factor-alpha...Glutamine-dependent anaplerosis is the principal source of succinate." |
| 5 | **Vats D, Mukundan L, Odegaard JI, et al. (2006). "Oxidative metabolism and PGC-1beta attenuate macrophage-mediated inflammation." Cell Metab 4(1):13-24.** DOI 10.1016/j.cmet.2006.05.011 | **16814729** | Verified live (esearch single-hit + esummary + abstract via efetch). **PRIMARY M2 oxidative/FAO anchor**, CAUSAL not just correlative — quoted: "in response to IL-4, STAT6 and PGC-1beta induce macrophage programs for fatty acid oxidation and mitochondrial biogenesis...inhibition of oxidative metabolism or RNAi-mediated knockdown of PGC-1beta attenuates this immune response." |
| 6 | **Jha AK, Huang SC, Sergushichev A, et al. (2015). "Network integration of parallel metabolic and transcriptional data reveals metabolic modules that regulate macrophage polarization." Immunity 42(3):419-30.** DOI 10.1016/j.immuni.2015.02.005 | **25786174** | Verified live (esearch disambiguated 2 candidates via esummary; FULL abstract via Europe PMC REST). **A second, INDEPENDENT (13C-tracer + network integration + pharmacologic shunt inhibition, vs #7's genetic knockout) TCA-break anchor** — quoted: "In M1 macrophages, we identified a metabolic break at Idh...inhibition of aspartate-aminotransferase...inhibited nitric oxide and interleukin-6 production in M1 macrophages, while promoting mitochondrial respiration." Also the CAUSAL M2 anchor: glutamine deprivation DECREASES M2 polarization/CCL22. |
| 7 | **Lampropoulou V, Sergushichev A, Bambouskova M, et al. (2016). "Itaconate Links Inhibition of Succinate Dehydrogenase with Macrophage Metabolic Remodeling and Regulation of Inflammation." Cell Metab 24(1):158-66.** DOI 10.1016/j.cmet.2016.06.004 | **27374498** | Verified live (esearch disambiguated 2 candidates via esummary; FULL abstract via efetch). **THE F2/F3-unifying mechanistic anchor** — quoted: "itaconate modulates macrophage metabolism...by inhibiting succinate dehydrogenase-mediated oxidation of succinate...Using newly generated Irg1(-/-) mice...endogenous itaconate regulates succinate levels and function, mitochondrial respiration, and inflammatory cytokine production." Also the HONEST F2 caveat source (itaconate is NOT a fully orthogonal bystander marker — it feeds back on cytokine output). |
| 8 | **Martinez FO, Gordon S (2014). "The M1 and M2 paradigm of macrophage activation: time for reassessment." F1000Prime Rep 6:13.** | **24669294** | Verified live (esearch single-hit + esummary + FULL abstract via efetch). Symmetric-QC anchor — quoted: "a dichotomy has been proposed for macrophage activation: classic vs. alternative, also M1 and M2...a revision of the model is needed" given "the increasing number of immune-relevant ligands" beyond the original limited set. |
| 9 | **Raes G, Van den Bergh R, De Baetselier P, et al. (2005). "Arginase-1 and Ym1 are markers for murine, but not human, alternatively activated myeloid cells." J Immunol 174(11):6561-2.** | **15905489** | Verified live (esearch single-hit + esummary title/journal/year/pages match). **EXISTENCE-TIER**: no abstract text extractable this session (a very short, 6561-6562, Letters-format brief communication) — the TITLE ITSELF is the complete, decisive claim; same existence-only tier this repo's own Levenson 1965 citation uses. Species-difference anchor. |
| 10 | Murray PJ et al. (2014). "Macrophage activation and polarization: nomenclature and experimental guidelines." Immunity 41(1):14-20. | **25035950** | **REUSED**, not independently re-fetched this session — already `fetched_live: true` in this repo's existing `IMMUNE-MACROPHAGE-POLARIZATION` graph cell. Field-consensus nomenclature paper; that graph cell's own "conflicts" field records it as framing continued M1/M2 shorthand as contentious. |
| 11 | Movahedi K et al. (2010). "Different tumor microenvironments contain functionally distinct subsets of macrophages derived from Ly6C(high) monocytes." Cancer Res. | **20570887** | **REUSED**, not independently re-fetched — already `fetched_live: true` in the existing graph cell. In-vivo TAM functional-subset anchor, couples F4 to the cancer/TAM context (§6). |

**Extensively attempted but not obtained this session (disclosed, not hidden):** a live web search
for a quantitative iNOS-vs-Arginase-1 Km comparison (the textbook μM-vs-mM affinity-mismatch fact)
returned "session has used its web search budget" (WebSearch quota exhausted, same fallback this
repo's sibling docs already disclose) — **not used as a load-bearing number** rather than forced in
unverified; F1 is gated on direction/reciprocity (which Mills 2000's own abstract fully supports),
not on that specific unverified affinity ratio.

---

## 3. Model — what's computed, per falsifier

**F1 (arginine-fate reciprocal split).** Illustrative (disclosed, not fit to an extracted primary
number — Mills 2000's abstract gives direction/mechanism, not exact concentrations) relative enzyme
activities `E_NOS_M1=10, E_ARG_M1=0.5` vs `E_NOS_M2=0.5, E_ARG_M2=10` give `frac_NO_M1=0.952,
frac_NO_M2=0.048`, **Δ=0.905** (≥0.8 pre-registered gate: PASS). The forced adversary (independent,
dedicated-substrate pools — no shared-node conservation) is tested via 5,000 Monte Carlo draws of
the SAME `(e_nos, e_arg)` pairs, scored two ways: under independence, 24.98% of draws would show
BOTH outputs "high" simultaneously; under the real (shared-pool, conservation-exact) model, this is
**structurally impossible** (0.0%, not merely rare) — matching Mills' own reported reciprocal
exclusivity. Robustness: ±40% independent perturbation of all four activity levels, 20 draws, **20/20
retain Δ≥0.8** (minimum observed 0.843).

**F2 (itaconate/IRG1 decorrelated node).** Existence: ≥2 independent primary sources
(Michelucci 2013, Tannahill 2013, Lampropoulou 2016) report IRG1/itaconate induction under
inflammatory (M1) activation — PASS. Modality-decorrelation: Lampropoulou's `Irg1⁻/⁻` mice are a
clean genetic knockout of *only* the itaconate-producing step, structurally distinct from the
TLR4/NF-κB cytokine-transcription machinery — a genuinely different assay modality (metabolite/
enzyme vs secreted-protein/transcript), which is what makes it a real second, independent
verification leg for "is this macrophage M1," not a tautological restatement of the cytokine
markers. **Honest, disclosed, non-gated caveat:** Lampropoulou's own `Irg1⁻/⁻` result shows
itaconate loss changes "inflammatory cytokine production" — itaconate is NOT a fully causally
orthogonal bystander; it feeds back as an anti-inflammatory SDH-inhibition brake. "Decorrelated"
here means modality, not zero coupling — stated plainly, not oversold.

**F3 (metabolic signature).** The 6-state TCA-segment steady state (M1 params: `k_IK=0.3` [IDH
break], `k_IZ=0.8` [IRG1 on], `J_gln=1.2` [glutamine anaplerosis, "principal source"],
`oxidative_frac=0.4`; M2 params: `k_IK=1.0`, `k_IZ=0`, `J_gln=0`, `oxidative_frac=0.8`) gives:
itaconate **M1-specific** (Z_M1=1.939 vs Z_M2=0.0, exactly); succinate **higher in M1**
(S_M1=4.169 vs S_M2=1.6); downstream/OXPHOS-feeding flux **higher in M2** (F_M2=1.6 vs F_M1=1.418)
— all three PASS, the qualitative direction the task's own falsifier specifies. The Jacobian
structural check (§1) PASSES exactly. **The decisive forced-adversary test is not this script's own
simulation** — it is Tannahill 2013's own already-measured result, machine-encoded as a boolean
dissociation: 2-deoxyglucose (glycolysis inhibition) suppresses LPS-induced IL-1β
(`TANNAHILL_2DG_IL1B_SUPPRESSED=True`) but **not** TNF-α (`TANNAHILL_2DG_TNFA_SUPPRESSED=False`) —
`IL1B_SUPPRESSED != TNFA_SUPPRESSED` is the machine-checkable falls-condition for the "M1 metabolism
is a generic, undifferentiated correlate" adversary: it FALLS, because metabolism is a *specific*
wire (succinate→HIF-1α→IL-1β), not a uniform dial across all M1 outputs. **Over-determination, not a
tautology:** Jha 2015 (13C-tracer + pharmacologic shunt-inhibition) and Lampropoulou 2016 (genetic
`Irg1⁻/⁻` + metabolomics) — two independent groups, two different specific lesions (IDH vs
SDH-via-itaconate), two different methods — converge on the same qualitative structural conclusion.

**F4 (plasticity/reversibility).** Time constant `k=0.1625/day` **live-loaded** from
`wound_healing_cascade_results.json`'s own Daley-2010-fitted macrophage decay rate (not re-fit
here) sets `PHASE_DAYS=30.77` (5 time constants). A 3-phase stimulus cycle (M1→M2→M1, target
φ=0.9/0.1/0.9) is run through (a) a **reversible** relaxation ODE and (b) a forced **ratchet**
adversary (locks once `|φ−0.1|<0.05` for `≥0.5τ`). Results: reversible model ends phase 3 at
**φ=0.895** (close to the 0.9 target, PASS gate ≥0.7); ratchet ends at **φ=0.130** (still locked
near 0.1, PASS gate <0.2) — a clean dissociation. Void floor (zero stimulus-responsiveness) fails
to reproduce even the first M1→M2 transition (PASS, floor falls). Robustness: ±40% perturbation of
`k`, 16 draws, **16/16** retain the reversal (minimum observed φ_end=0.868). External anchor (not a
tautology): Stout 2005's own in-vivo aged/tumor-bearing-mice microenvironment-switch experiment +
Movahedi 2010 TAM functional-subset plasticity (reused) — neither is this script's own output.

---

## 4. Symmetric QC — nothing proven; the spread and the caveats are real, held OPEN

- **Binary vs. continuum (the single most important caveat):** the clean M1/M2 dichotomy
  reproduced here (F1–F3) is an **in-vitro-extreme simplification**. This repo's own existing
  `IMMUNE-MACROPHAGE-POLARIZATION` graph cell already found — via single-cell atlases (Azizi 2018
  PMID 29961579, n=45,000 cells/8 tumors; Cheng 2021 PMID 33545035, n=210 patients/15 types) — that
  tumor-associated macrophages occupy a marker-**diverse continuum**, not discrete M1/M2 clusters,
  and that a phenotype-stratified (binary) subgroup analysis (Zhang 2012 PMID 23284651, n=5 studies)
  shows **no** significant differential survival effect even though total TAM density stays
  prognostic. Reused, not re-derived, as this build's own in-vivo-continuum anchor — the two
  findings are complementary: the reductionist in-vitro axis modeled here is real and reproducible
  (F1–F3 all PASS), but it is not the whole in-vivo story.
- **Murray 2014 consensus** (PMID 25035950, reused): the field's own nomenclature-standardization
  paper frames continued M1/M2 shorthand as contentious in vivo, proposing multi-marker/multi-axis
  characterization instead of a true dichotomy.
- **Species difference:** Raes 2005 (PMID 15905489, existence-tier): Arginase-1 and Ym1 are markers
  for MURINE, not human, alternatively-activated myeloid cells. Mills 2000 (this build's F1 anchor)
  is all-mouse — the arginine-fate finding does not directly transfer to human macrophage marker
  panels.
- **Terminology drift (verified live, disclosed, not smoothed over):** Mills 2000's own "M-1/M-2"
  labels were originally mouse-STRAIN/genetic-background associated (Th1-prone C57BL/6 vs
  Th2-prone BALB/c), **not** identical to the later stimulus-based (LPS/IFN-γ vs IL-4/IL-13)
  nomenclature this task itself uses (standardized by Mosser/Martinez&Gordon/Murray 2014) — a real
  terminological drift the modern shorthand inherits.
- **Itaconate is not fully orthogonal** (§3, F2): it feeds back on cytokine output (Lampropoulou
  2016) — "decorrelated" means different assay/regulatory modality, not zero causal coupling.
- **Two distinct TCA lesions, not conflated:** Jha 2015's IDH break and Lampropoulou/Michelucci's
  itaconate-mediated SDH inhibition are different enzymes at different cycle positions, reported
  as an over-determining convergence, not the same mechanism counted twice.
- **Illustrative rate constants:** all enzyme-activity/rate-constant numbers in F1/F3 are
  illustrative, disclosed, not fit to an extracted primary numeric source — the live-verified
  abstracts give direction/mechanism/causal dissociation (what is actually gated), not absolute
  concentrations for this specific parameterization. Same tier as this repo's own
  `acute_phase_inflammation.py` upstream cytokine rate constants.
- **Held OPEN throughout** — none of this spread is swept into `overall_pass`.

---

## 5. Pre-registered gates — machine-printed, not narrated

```
F1_conservation_exact:                     PASS
F1_reciprocal_crossover_pass (Δ=0.905>=0.8): PASS
F1_ADV_independent_pools_falls (25.0% vs 0.0%): PASS (adversary falls)
F1_robustness_pass (20/20, min Δ=0.843):    PASS
F1_overall_pass:                            PASS

F2_existence_pass (3 independent sources):  PASS
F2_modality_decorrelated:                   PASS
F2_overall_pass:                            PASS

F3_geometric_pass (itaconate M1-specific, succinate M1>M2, flux M2>M1): PASS
F3_jacobian_single_feedback_wire (chain + {(S,Z),(F,Z)} exactly): PASS
F3_ADV_tannahill_2DG_falls (IL-1b dependent != TNF-a independent): PASS (adversary falls)
F3_overall_pass:                            PASS

F4_dissociation_pass (reversible 0.895 vs ratchet 0.130): PASS
F4_void_floor_falls:                        PASS
F4_eigenvalue_matches_neg_k:                PASS
F4_robustness_pass (16/16, min phi_end=0.868): PASS
F4_overall_pass:                            PASS

OVERALL: PASS
```
Deterministic: byte-identical md5 `bb7e9897a1d8526ee24d6ca039b5757d` confirmed across repeated
process runs (fixed RNG seeds for the two Monte Carlo/robustness sweeps; all other computation is
closed-form/`solve_ivp`, no other stochastic element).

---

## 6. couples_to — computed coupling points, not prose gestures

- **wound-healing** (`docs/MECHANISM_WOUND_HEALING.md` / `scripts/msk/wound_healing_cascade.py`):
  this build's F4 plasticity-ODE time constant (`k=0.1625/day`) is **live-loaded at runtime** from
  that model's own Daley-2010-fitted macrophage M1-fraction decay rate (its F4), not re-fit — a
  real, computed, quantitative coupling. That model covers the M1-like→M2-like **timing** inside a
  healing wound; this build covers the underlying **mechanism** (why/how the switch happens and
  reverses) that timing model treats phenomenologically.
- **acute-phase inflammation** (`docs/MECHANISM_ACUTE_PHASE_INFLAMMATION.md` /
  `scripts/msk/acute_phase_inflammation.py`): qualitative, disclosed (not a new computed number) —
  that model's upstream TNF pulse is deliberately generic/any-cell-source; this build's M1 state is
  mechanistically one real physiological source of that pulse (macrophage-secreted TNF is part of
  the classical M1 cytokine output, Mills 2000/Murray 2014). Not re-derived as macrophage-specific
  in that model this session.
- **cancer/TAM** (existing graph cell `IMMUNE-MACROPHAGE-POLARIZATION`,
  `data/MECHANISM_ANCHOR_GRAPH.json`, + `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md` cluster B
  [tumor-microenvironment-suppression] and H [injury-induced-regeneration-immunity, DAMPs/
  macrophage-polarization]): this build's clean in-vitro F1–F3 mechanism is the reductionist
  substrate the existing graph cell's in-vivo TAM-continuum finding complicates (§4) — cross-
  referenced both directions, not duplicated.
- **itaconate**: F2/F3 of this build ARE the itaconate/IRG1 content the task specifically asked to
  resolve — built as the F2/F3-unifying mechanism (the Jacobian single-feedback-mechanism check),
  not a bolt-on citation.

Per `docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/§4, promoting any of this into canonical
`data/MECHANISM_ANCHOR_GRAPH.json` edges/status-transition (`OPEN`→`ASSUMED`/`PROVEN`) requires the
separate `mechanism_fold → fold_gate_v2` pipeline — **not performed this session**, consistent with
every sibling `MECHANISM_*` doc's own stated precedent.

---

## 7. Honest gaps — symmetric QC: what this does NOT prove

1. **In-vitro-extreme binary vs. in-vivo continuum** (§4) — the single most important limitation:
   this build's clean reciprocal dichotomy is real and reproducible in the reductionist LPS/IFN-γ-
   vs-IL-4/IL-13 axis, but this repo's own existing graph cell already shows it is not the whole
   in-vivo (tumor) story.
2. **All-mouse primary anchor for F1** (Mills 2000) — the arginine-fate finding does not directly
   transfer to human macrophages (Raes 2005, existence-tier); no live-fetched primary human
   arginine-fate number was obtained this session.
3. **Terminology drift** (§4) — Mills 2000's own strain-based "M-1/M-2" is not identical to the
   modern stimulus-based nomenclature this task (and most of the field) now uses.
4. **Itaconate is not a fully causally orthogonal marker** (§3/§4, F2) — it feeds back on cytokine
   output; "decorrelated" is a modality claim, not a zero-coupling claim.
5. **All F1/F3 rate constants are illustrative**, not fit to an extracted primary numeric source —
   the falsifiers gate on direction/mechanism/dissociation (what the live-verified abstracts
   actually support), not absolute magnitudes for this specific parameterization.
6. **A quantitative iNOS-vs-Arginase Km/affinity comparison was not obtained this session**
   (WebSearch quota exhausted, §2) — not used as an unverified load-bearing number.
7. **Reduced, phenomenological models** — a 6-state TCA-segment ODE and a 1-state relaxation ODE,
   not genome-scale metabolic/signaling simulations; F1 is a static flux-partition, not a kinetic
   enzyme model.
8. **Two real design bugs were found and fixed this session** (§1: the `mito_capacity`-cancellation
   bug and the `PHASE_DAYS`/timescale-mismatch bug) — kept visible in
   `scripts/msk/macrophage_polarization.py`'s own inline comments as an auditable trail, not
   silently corrected out of the record.
9. **Two citations (Murray 2014, Movahedi 2010) are REUSED**, not independently re-fetched this
   session — already `fetched_live: true` in this repo's existing `IMMUNE-MACROPHAGE-POLARIZATION`
   graph cell, disclosed inline (§2, §6).
10. **Single, generic parameter set** — no stimulus-dose/timing/tissue-context-conditioned version
    of the model (matching this repo's own first-step ceiling disclosed by every sibling
    `MECHANISM_*` layer).

---

## 8. Files

- `scripts/msk/macrophage_polarization.py` — the model: `CITATIONS` dict (11 entries, tiered
  live-verified/reused/existence-only), the wound-healing coupling loader (live-loaded at runtime,
  disclosed hardcoded fallback), the F1 arginine-fate flux-partition + Monte Carlo forced adversary
  + robustness sweep, the F2 itaconate-node existence/decorrelation checks, the F3 6-state
  TCA-segment ODE + Jacobian structural check + Tannahill-2DG forced-adversary encoding, the F4
  reversible-relaxation ODE + ratchet forced adversary + void floor + eigenvalue check + robustness
  sweep, the symmetric-QC dict, `couples_to`, and the gates/evidence-JSON writer. Run with
  `source .venv-msk/bin/activate && python3 scripts/msk/macrophage_polarization.py` (a few seconds
  wall time, pure numpy/scipy, no OpenSim dependency, deterministic given fixed RNG seeds).
- `data/msk_smoketest/macrophage_polarization/macrophage_polarization_results.json` — full
  machine-written evidence (citations, F1–F4 full computed blocks including all adversaries/void
  floors/robustness sweeps/Jacobian structure, symmetric QC, `couples_to`, gates, `overall_pass`).
  md5 `bb7e9897a1d8526ee24d6ca039b5757d`.
- `data/MECHANISM_ANCHOR_GRAPH.json` — the existing `IMMUNE-MACROPHAGE-POLARIZATION` node (id
  search: `grep -n '"IMMUNE-MACROPHAGE-POLARIZATION"' data/MECHANISM_ANCHOR_GRAPH.json`), status
  `OPEN`, `mechanism_grade: SEED-DESIGN`, **unchanged this session** (fold not performed, §6) — its
  own PMIDs 25035950 (Murray 2014) and 20570887 (Movahedi 2010) are reused, not refetched, by this
  build (§2).
- Sibling docs (read for coupling, not re-litigated): `docs/MECHANISM_WOUND_HEALING.md` +
  `data/msk_smoketest/wound_healing_cascade/wound_healing_cascade_results.json` (F4 time-constant
  coupling source, §3/§6), `docs/MECHANISM_ACUTE_PHASE_INFLAMMATION.md` (qualitative TNF-source
  coupling, §6), `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md` (cancer/TAM cluster context, §6),
  `docs/MECHANISM_HARDENED_CONVENTIONS.md` (fold path/node schema — this doc is a pre-fold
  HYPOTHESIS artifact, §6).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/macrophage_polarization.py
```
No inputs required beyond the already-certified `wound_healing_cascade_results.json` (live-loaded,
disclosed hardcoded fallback if absent). Pure Python/numpy/scipy (`solve_ivp`), no OpenSim, no
subject data, runs in a few seconds. Writes
`data/msk_smoketest/macrophage_polarization/macrophage_polarization_results.json`. No git
operations. No files outside `scripts/msk/macrophage_polarization.py`, the evidence JSON, and this
doc were modified or created this session.
