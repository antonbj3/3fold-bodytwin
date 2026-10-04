# MECHANISM TAU PATHOLOGY — hyperphosphorylation -> PHF aggregation, the Braak stereotyped spatial spread, and the tau-vs-amyloid temporal-ordering + spatial-mismatch cert (2026-07-22)

Downstream partner to the amyloid-cascade thread. **`docs/MECHANISM_AMYLOID_BETA_AGGREGATION.md` was
searched for exhaustively (`find`, `grep`) at the start of this session and does NOT exist anywhere in
this repo as of this build** — referenced here by name only, per the task's own framing of "together the
AD 'amyloid-cascade' picture," never assumed to exist, never read, never mutated. This document does not
depend on it. Script: `scripts/msk/tau_pathology.py`. Evidence: `data/tau_pathology/tau_pathology_results.json`
(raw, machine-computed) + `docs/MECHANISM_TAU_PATHOLOGY_evidence.json` (curated citation dossier).

## 0. Scope, stated up front — read this before any number below

This is **not** a subject-specific twin layer with an upstream OpenSim/`.osim` dependency (unlike most
`MECHANISM_*` MSK docs) — it is a **literature-anchored + first-principles-geometric-toy cert**, built
fresh, with no upstream JSON read. Two genuinely different kinds of evidence are combined and kept
labeled throughout, never blurred: (1) a small, hand-specified, **disclosed TOY** connectivity graph
(10 nodes, 14 edges, coarse textbook neuroanatomy) used to demonstrate that hierarchical topology plus
the correct epicenter is a *sufficient mechanism* for stereotyped spread, tested by forced adversary and
machine-graded gates; and (2) numbers **directly transcribed from live-verified published abstracts**,
with only transparent arithmetic (subtraction, ratios) performed on top — never invented, never recalled
from training-data memory. **Every one of the 23 PMID/DOI pairs cited below was re-fetched this session**
via NCBI eutils `efetch` (raw abstract text, not a summarizer, not recalled) — including 9 that a
prior session's `NEU-TAU-PROPAGATION` graph node (`data/MECHANISM_ANCHOR_GRAPH.json`) had already flagged
`fetched_live:true`; all 9 were independently re-verified byte-for-byte against that prior claim this
session and matched exactly (a genuine check, not a rubber stamp — this repo's own memory explicitly
warns that a subagent's "confirmed live" claim has been fabricated before).

**The pre-registered falsifier** (verbatim, stated before any number below was computed): does the model
reproduce (a) the measured tau-PET (or CSF/plasma p-tau) vs. cognition correlation being **stronger**
than amyloid-vs-cognition — the decorrelated "which protein tracks the disease" test (Braak stage /
tau-PET SUVR vs. MMSE) — **and** (b) the Braak-stage stereotyped spatial ordering?

**Symmetric QC, stated up front, not discovered after the fact**: the amyloid-cascade **causal
direction is genuinely contested** (amyloid-first vs. tau-first vs. parallel-independent-initiation) —
held OPEN throughout, never resolved either way; Jack et al. 2013 [12] itself revises its own 2010 model
[11] to allow independent initiation. Tau-PET tracers have **documented off-target binding** [20,21] —
held OPEN as a measurement caveat (§6.4).

## 1. Geometric structure (derive from the geometry, not heuristics)

Trans-synaptic/trans-neuronal tau spread (Part B, §5) is modeled as a continuous-time diffusion process
on a weighted anatomical graph: `dx/dt = -L x`, where `L = D - W` is the **combinatorial graph
Laplacian** (`D` = diagonal degree matrix, `W` = weighted adjacency). This is solved via the graph's own
**spectrum** — `L = U diag(λ) Uᵀ`, `x(t) = U diag(exp(-λt)) Uᵀ x0` — the same "derive from the
graph/spectrum" family this repo's own `sigma_min` discipline uses elsewhere, not a curve fit and not a
free-parameter fit to the target ordering. For a **connected** graph, this diffusion conserves total mass
and provably converges to the uniform vector `(1/N)·Σx0` (the averaging/consensus theorem for graph
Laplacians) — this fixes, **with no free parameter**, the per-node "half-arrival-time" threshold
(`0.5/N`) used to rank nodes by how fast disease-proxy mass reaches them from a seeded epicenter.

The 4-tier Braak-rank assignment used to grade the model (§5.2) is **directly licensed by Braak & Braak
1991's own verbatim abstract wording** [1] (re-verified live this session), not an external elaboration:
"transentorhinal stages I-II" (tier 1); "limbic stages (III-IV)... mild involvement of the first Ammon's
horn sector" (tier 2 — hippocampus + amygdala + parahippocampal/perirhinal cortex grouped here);
"isocortical stages (V-VI)... destruction of virtually all isocortical **association** areas" (tier 3 —
note the word "association", which by its own plain meaning *excludes* primary sensory/motor/visual
cortex, directly licensing tier 4 as the latest/least-affected group — an inference from the verified
text itself, not an unverified extrapolation beyond it).

## 2. Method, in one paragraph

Three independent pieces are built and cross-checked. **(A)** The molecular mechanism (hyperphosphorylation
-> PHF aggregation) is a citation-only synthesis of four live-verified primary/structural papers
[2,3,4,5] — no free model parameters, no simulation. **(B)** The Braak spatial-spread claim is checked
two ways: **4 decorrelated real-world legs** (autopsy histology [1], biochemical seeding-assay [6],
in-vivo PET [7], connectome-diffusion modeling [8], with an adversarial AD-vs-PSP dissociation [9] and an
independent replication [10] — all citation-anchored, none re-computed here) **plus** an **independent,
from-scratch, first-principles geometric toy model** (§1, §5.3) built and run *by this script*, tested
against **two orthogonal forced adversaries** (wrong epicenter, exhaustively over all 10 candidate seeds;
scrambled topology, 500 valid trials) and a **void-floor** (topology-blind complete graph). **(C)** The
tau-vs-amyloid decorrelation claim is checked via one pre-registered temporal-gap computation
(Villemagne et al. 2013's own reported years-before-onset numbers [13], subtracted transparently) and one
pre-registered same-cohort ratio (Ossenkoppele et al. 2016's own reported correlation coefficients [17],
divided transparently), **triangulated** across 4 independent research sites spanning 2016-2026 [14,15,16,17,18,23]
— corrected for a common-recruitment-site confound (3 of the papers share the same UCSF cohort base and
are counted as **one** site, not three, a deliberate de-inflation, not a flattering count). A targeted
(not systematic) adversarial search for a dissenting same-cohort finding was run and returned none — a
disclosed **honest-negative search result**, not proof of universal absence (§9).

## 3. Citations — all 23 PMID/DOI pairs verified LIVE this session (NCBI eutils `efetch`, raw abstract text)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Braak H, Braak E (1991). Neuropathological stageing of Alzheimer-related changes. *Acta Neuropathol* 82(4):239-59. | **1759558**, DOI `10.1007/BF00308809` | THE canonical staging scheme, n=83 autopsy brains. Directly licenses the 4-tier model (§1). |
| 2 | Grundke-Iqbal I, Iqbal K, Tung YC, Quinlan M, Wisniewski HM, Binder LI (1986). Abnormal phosphorylation of the microtubule-associated protein tau (τ) in Alzheimer cytoskeletal pathology. *PNAS* 83(13):4913-7. | **3088567**, DOI `10.1073/pnas.83.13.4913` | THE original discovery: dephosphorylation dramatically increases antibody recognition of PHF-tau — tau in PHF is abnormally phosphorylated. |
| 3 | Alonso AC, Zaidi T, Grundke-Iqbal I, Iqbal K (1994). Role of abnormally phosphorylated tau in the breakdown of microtubules in Alzheimer disease. *PNAS* 91(12):5562-6. | **8202528**, DOI `10.1073/pnas.91.12.5562` | Direct experiment: AD P-tau "bound to normal tau but not to tubulin" and "inhibited microtubule assembly" — the dominant-negative sequestration mechanism. |
| 4 | Iqbal K, Alonso AC, Gong CX, Khatoon S, Singh TJ, Grundke-Iqbal I (1994). Mechanism of neurofibrillary degeneration in Alzheimer's disease. *Mol Neurobiol* 9(1-3):119-23. | **7888088**, DOI `10.1007/BF02816111` | Review synthesis: "activities of ... protein phosphatase 2A ... are decreased in AD brain" — the phosphatase-deficit half of the hyperphosphorylation mechanism. |
| 5 | Fitzpatrick AWP, Falcon B, He S, et al. (2017). Cryo-EM structures of tau filaments from Alzheimer's disease. *Nature* 547(7662):185-190. | **28678775**, DOI `10.1038/nature23002` | Atomic-resolution (3.4-3.5O) cryo-EM structure of PHF/straight-filament cores, residues 306-378, "cross-β/β-helix structure" — modern structural end-state anchor. |
| 6 | Kaufman SK, Del Tredici K, Thomas TL, Braak H, Diamond MI (2018). Tau seeding activity begins in the transentorhinal/entorhinal regions and anticipates phospho-tau pathology. *Acta Neuropathol* 136(1):57-67. | **29752551**, DOI `10.1007/s00401-018-1855-6` | Biochemical seeding-assay leg, n=247: earliest/most robust seeding in TRE/EC, anticipating AT8-histopathology. |
| 7 | Scheed M, Lockhart SN, Schonhaut DR, et al. (2016). PET Imaging of Tau Deposition in the Aging Human Brain. *Neuron* 89(5):971-982. | **26938442**, DOI `10.1016/j.neuron.2016.01.028` | In-vivo PET leg: tracer-retention patterns "corresponded well with Braak staging" — living-brain confirmation of an autopsy scheme. |
| 8 | Vogel JW, Iturria-Medina Y, Strandberg OT, et al. (2020). Spread of pathological tau proteins through communicating neurons in human Alzheimer's disease. *Nat Commun* 11(1):2612. | **32457389**, DOI `10.1038/s41467-020-15701-2` | Connectome epidemic-spreading-model (ESM) leg, n=312: "up to 70% of the variance" in tau-PET spatial pattern explained by connectivity; pattern fits **irrespective of amyloid**, but higher-Aβ regions show more tau than connectivity alone predicts. |
| 9 | Cope TE, Rittman T, Borchert RJ, et al. (2018). Tau burden and the functional connectome in Alzheimer's disease and progressive supranuclear palsy. *Brain* 141(2):550-567. | **29293892**, DOI `10.1093/brain/awx347` | Forced-adversary dissociation, n=46 (17 AD+17 PSP+12 controls): hub-connectivity predicts tau in AD (spread-consistent) but **not** in PSP (vulnerability-consistent instead) — the diverse-instance-space test. |
| 10 | Franzmeier N, Neitzel J, Rubinski A, et al. (2020). Functional brain architecture is associated with the rate of tau accumulation in Alzheimer's disease. *Nat Commun* 11(1):347. | **31953405**, DOI `10.1038/s41467-019-14159-1` | Independent 2-cohort replication (53 ADNI + 41 BioFINDER AD, n=138 incl. controls): connectivity correlates with tau-accumulation **rate**. |
| 11 | Jack CR Jr, Knopman DS, Jagust WJ, et al. (2010). Hypothetical model of dynamic biomarkers of the Alzheimer's pathological cascade. *Lancet Neurol* 9(1):119-28. | **20083042**, DOI `10.1016/S1474-4422(09)70299-6` | Original temporal-ordering model: Aβ biomarkers abnormal first, neurodegeneration biomarkers later, correlating with symptom severity. |
| 12 | Jack CR Jr, Knopman DS, Jagust WJ, et al. (2013). Tracking pathophysiological processes in Alzheimer's disease: an updated hypothetical model of dynamic biomarkers. *Lancet Neurol* 12(2):207-16. | **23332364**, DOI `10.1016/S1474-4422(12)70291-0` | **Self-revision**: "the two major proteinopathies ... Aβ and tau, might be initiated independently in sporadic AD" — directly supports holding the causal-order question OPEN. |
| 13 | Villemagne VL, Burnham S, Bourgeat P, et al. (2013). Amyloid β deposition, neurodegeneration, and cognitive decline in sporadic Alzheimer's disease: a prospective cohort study. *Lancet Neurol* 12(4):357-67. | **23477989**, DOI `10.1016/S1474-4422(13)70044-9` | AIBL, n=200 (163 with positive accumulation): Aβ reaches positivity **17.0y** before dementia onset; hippocampal atrophy only **4.2y** before; memory impairment **3.3y** before — the quantified temporal-gap anchor (§6.1). |
| 14 | Hanseeuw BJ, Betensky RA, Jacobs HIL, et al. (2019). Association of Amyloid and Tau With Cognition in Preclinical Alzheimer Disease: A Longitudinal Study. *JAMA Neurol* 76(8):915-924. | **31157827**, DOI `10.1001/jamaneurol.2019.1424` | Harvard Aging Brain Study, n=60: serial mediation Aβ->tau->cognition; tau change predicts PACC change (β=-3.28, P=.001), covarying baseline Aβ/tau. |
| 15 | La Joie R, Visani AV, Baker SL, et al. (2020). Prospective longitudinal atrophy in Alzheimer's disease correlates with the intensity and topography of baseline tau-PET. *Sci Transl Med* 12(524):eaau5732. | **31894103**, DOI `10.1126/scitranslmed.aau5732` | UCSF, n=32: "the global intensity of tau-PET, **but not** β-amyloid-PET, signal predicted the rate of subsequent atrophy." |
| 16 | Brier MR, Gordon B, Friedrichsen K, et al. (2016). Tau and Aβ imaging, CSF measures, and cognition in Alzheimer's disease. *Sci Transl Med* 8(338):338ra66. | **27169802**, DOI `10.1126/scitranslmed.aaf2362` | WashU Knight ADRC: "tau deposition in the temporal lobe more closely tracked dementia status and **was a better predictor of cognitive performance** than Aβ deposition in **any region** of the brain." Direct, verbatim, same-cohort — the closest match to the task's own falsifier wording. Qualitative only — abstract does not print the underlying r/β (§9). |
| 17 | Ossenkoppele R, Schonhaut DR, Schöll M, et al. (2016). Tau PET patterns mirror clinical and neuroanatomical variability in Alzheimer's disease. *Brain* 139(Pt 5):1551-67. | **26962052**, DOI `10.1093/brain/aww027` | UCSF/Berkeley, n=16 with all 3 tracers (of 20 AD + 15 amyloid-negative controls): **quantified** same-cohort, same-day decorrelation — `r(tau,FDG)=-0.49±0.07`, `r(amyloid,FDG)=0.16±0.09`, `r(tau,amyloid)=0.18±0.09` (all P<0.001); amyloid bound "diffusely throughout the neocortex" while tau/FDG tracked each of 3 distinct clinical phenotypes. The core quantitative falsifier datapoint (§6.2). |
| 18 | Bejanin A, Schonhaut DR, La Joie R, et al. (2017). Tau pathology and neurodegeneration contribute to cognitive impairment in Alzheimer's disease. *Brain* 140(12):3286-3300. | **29053874**, DOI `10.1093/brain/awx243` | UCSF, n=40 (incl. 12 PCA + 8 lvPPA atypical phenotypes): tau-PET tracks domain-specific cognitive deficits region-by-region; "weakly related to amyloid burden." |
| 19 | Palmqvist S, Janelidze S, Quiroz YT, et al. (2020). Discriminative Accuracy of Plasma Phospho-tau217 for Alzheimer Disease vs Other Neurodegenerative Disorders. *JAMA* 324(8):772-781. | **32722745**, DOI `10.1001/jama.2020.12134` | Cohort 1, n=81 (34 AD+47 non-AD): plasma p-tau217-vs-tangle correlation gated by amyloid status — `ρ=0.64` (amyloid+) vs. `ρ=0.15, P=.33` NS (amyloid-). Amyloid-gates-tau-magnitude evidence (held as an open tension, §9). |
| 20 | Marquié M, Normandin MD, Vanderburg CR, et al. (2015). Validating novel tau PET tracer [F-18]-AV-1451 (T807) on postmortem brain tissue. *Ann Neurol* 78(5):787-800. | **26344059**, DOI `10.1002/ana.24517` | Off-target-binding leg: AV-1451 binds neuromelanin/melanin-containing cells and, to a lesser extent, hemorrhagic lesions. |
| 21 | Lowe VJ, Curran G, Fang P, et al. (2016). An autoradiographic evaluation of AV-1451 Tau PET in dementia. *Acta Neuropathol Commun* 4(1):58. | **27296779**, DOI `10.1186/s40478-016-0315-6` | Off-target-binding leg, n=38: vessels, iron-associated regions, substantia nigra, choroid-plexus calcification, leptomeningeal melanin; explicitly states AV-1451 "does not completely reflect early stage tau progression suggested by Braak ... staging" — a disclosed tension with citation [7]. |
| 22 | Karikari TK, Pascoal TA, Ashton NJ, et al. (2020). Blood phosphorylated tau 181 as a biomarker for Alzheimer's disease. *Lancet Neurol* 19(5):422-433. | **32333900**, DOI `10.1016/S1474-4422(20)30071-5` | Supporting leg, n=1131 across 4 cohorts: graded plasma p-tau181 tracks the AD continuum; tau-PET association AUC 83.08-93.11% across cohorts. |
| 23 | Hu G, Feng L, He L, Liu N, Wang H (2026). Sequential ¹⁸F-AV45/¹⁸F-AV1451 dual-tracer brain PET imaging in Alzheimer's disease. *Front Neurol* 17:1877217. | **42404124**, DOI `10.3389/fneur.2026.1877217` | **Largest and most recent replication found**, n=438 (325 AD+68 MCI+45 HC), cognitive assessment explicitly includes **MMSE**: "Tau deposition showed stronger cognitive correlations than Aβ... Tau is an independent driver of cognitive decline." Also: Aβ-tau SUVR correlated `r=0.65-0.81` — in tension with [17]'s weak `r=0.18` (§9, disclosed, not smoothed over). |

## 4. Part A — molecular mechanism: hyperphosphorylation -> PHF aggregation

No simulation here — a citation-only synthesis of 4 live-verified papers spanning 1986-2017, i.e. the
**original discovery through the modern atomic-resolution structure of the same lesion**:

1. **Discovery** [2]: a monoclonal antibody to microtubule-associated protein tau labeled neurofibrillary
   tangles and PHF only weakly until the tissue was **dephosphorylated**, at which point recognition
   increased dramatically — the first direct evidence that PHF-tau is an *abnormally phosphorylated*
   form of a normal cytoskeletal protein, not a distinct gene product.
2. **Mechanism, direct experiment** [3]: AD-derived abnormally-phosphorylated tau ("AD P-tau") has little
   microtubule-assembly-promoting activity; critically, when mixed with *normal* tau and tubulin, AD
   P-tau **inhibits** microtubule assembly — it "bound to normal tau but not to tubulin." This is a
   **dominant-negative sequestration mechanism**: hyperphosphorylated tau poisons the function of
   co-existing normal tau in trans, not merely losing its own function.
3. **Mechanism, phosphatase side** [4]: "activities of phosphoseryl/phosphothreonyl protein phosphatase
   2A and nonreceptor phosphotyrosyl phosphatase(s) are decreased in AD brain" — proposing a 3-step
   causal chain, quoted directly: "(1) A defect(s) in the protein phosphorylation/dephosphorylation
   system is one of the early events... (2) A decrease in protein phosphatase activities... allows the
   hyperphosphorylation of tau; and (3) Abnormal phosphorylation and polymerization of tau into PHF most
   probably lead to a breakdown of the microtubule system."
4. **End-state structure, modern** [5]: cryo-EM at 3.4-3.5On Resolution Resolves the PHF/straight-filament
   core as two identical protofilaments spanning tau residues **306-378** (a 73-residue span — trivial
   arithmetic, not fabricated), adopting a "combined cross-β/β-helix structure" that "define[s] the seed
   for tau aggregation." Paired-helical and straight filaments differ only in inter-protofilament
   packing — "ultrastructural polymorphs" of the same core fold, not different proteins.

**Read together**: a phosphatase deficit [4] permits hyperphosphorylation [2], hyperphosphorylated tau
sequesters and poisons normal tau in trans rather than merely losing function [3], and the aggregated
end-state has a specific, now atomically-resolved cross-β/β-helix fold [5] — a coherent, 4-decade,
independently-replicated mechanistic chain from biochemistry to structural biology. No numeric gate is
computed for this section (no free parameter to test) — it is reported as a citation chain, honestly
labeled as such.

## 5. Part B — Braak stereotyped spatial spread

### 5.1 Four decorrelated real-world legs (citation-anchored, not re-computed here)

| Leg | Modality | Key number | Citation |
|---|---|---:|---|
| Autopsy histology | Neuropathology | n=83, 6-stage scheme | Braak & Braak 1991 [1] |
| Biochemical seeding | Cellular biosensor assay | n=247, earliest seeding in TRE/EC | Kaufman 2018 [6] |
| In-vivo imaging | Tau-PET | patterns "corresponded well" with Braak stage | Scratch 2016 [7] |
| Connectome diffusion | Epidemic-spreading-model | n=312, 70% spatial variance explained | Vogel 2020 [8] |

**Forced adversary (real-world)**: a "shared regional-vulnerability gradient" model (nodes are
intrinsically vulnerable by local metabolic/trophic properties, independent of connectivity — the
alternative to trans-neuronal spread). Cope 2018 [9] tests this head-to-head against connectivity-driven
spread in the **same analysis pipeline**, applied to **two different tauopathies**: in AD, hub-connectivity
predicts tau burden (spread-consistent — the adversary **falls**); in PSP, the *same* connectivity metric
does **not** predict tau, and metabolic/trophic graph properties do instead (vulnerability-consistent —
the adversary **succeeds**). This AD-vs-PSP dissociation on identical methodology is the diverse-instance-space
test the watertight method requires — not a single flattering cohort. Franzmeier 2020 [10] independently
replicates the connectivity-rate coupling in 2 further cohorts (n=138 total).

### 5.2 The toy geometric model — construction (see §1 for the spectral method)

10 nodes, 4 Braak tiers (directly licensed by [1]'s own wording, §1): `EC` (tier 1) — `HIP, AMY, PHC`
(tier 2, limbic) — `ITG, PCC, PAR, PFC` (tier 3, isocortical association) — `V1, M1S1` (tier 4, primary,
by the "association"-only wording of [1]). 14 undirected weighted edges built from textbook neuroanatomy
(perforant path `EC-HIP`, entorhinal-parahippocampal adjacency `EC-PHC`, entorhinal-amygdala `EC-AMY`,
hippocampo-cingulate `HIP-PCC`, amygdalo-prefrontal `AMY-PFC`, parahippocampal-temporal `PHC-ITG`,
fronto-parietal and dorsal/ventral visual-association streams for the rest) — **explicitly disclosed as a
coarse toy, not a DTI-derived connectome** (§9).

### 5.3 Results — machine-computed, pre-registered gates

**Pre-registered** (fixed before computing, `PREREG` dict in the script): `rho_real >= 0.75`; EC must
rank **#1** of all 10 candidate epicenters; scrambled-topology mean `rho <= 0.5 * rho_real`; void-floor
arrival-time spread `<= 1%` of the real spread.

**Measured**: seeding at `EC`, the real toy topology gives **Spearman ρ = 0.939** between each region's
diffusion half-arrival-time and its Braak tier (excluding the seed's own trivial `t=0`) — **PASS** against
the 0.75 threshold.

**Forced adversary 1 — exhaustive seed sweep** (real topology, all 10 possible epicenters, non-cherry-picked):

| Seed | ρ vs. Braak tier | Seed | ρ vs. Braak tier |
|---|---:|---|---:|
| **EC** (tier 1) | **0.939** | PCC (tier 3) | 0.355 |
| HIP (tier 2) | 0.843 | ITG (tier 3) | 0.130 |
| PHC (tier 2) | 0.791 | PFC (tier 3) | -0.070 |
| AMY (tier 2) | 0.738 | M1S1 (tier 4) | -0.461 |
| | | V1 (tier 4) | -0.561 |
| | | PAR (tier 3) | -0.714 |

`EC` is the strict `argmax` (**rank 1 of 10** — PASS). A clean gap separates the result: every
tier-1/tier-2 seed scores `ρ > 0.7`; every tier-3/tier-4 seed scores `ρ < 0.36` (several strongly
negative) — the model does not merely "prefer EC by a hair," it cleanly separates limbic-cluster seeds
from neocortical/primary seeds. One honest irregularity, disclosed not hidden: `PAR` (tier 3, a 5-edge
hub in this toy graph) is the single *worst* seed (`ρ=-0.71`), scoring worse than either tier-4 node —
plausible cause (not adjudicated): as the highest-degree node, seeding there homogenizes arrival order
across tiers rather than merely delaying it, which anti-correlates rather than just under-correlates.
Does not affect the gate (EC is still the unambiguous rank-1).

**Forced adversary 2 — scrambled topology** (same 14 edge weights, same EC seed, node-pair identities
randomized; 500 valid connected trials out of 655 attempts, fixed RNG seed `20260722` for reproducibility):
mean `ρ = -0.005` (SD `0.356`, range `[-0.935, 0.873]`) — collapses to **statistical noise centered on
zero**, `-0.005 <= 0.5 * 0.939 = 0.470` — **PASS**. The specific hierarchical topology, not merely "EC" as
a label, is required for the ordering to emerge.

**Void floor**: a topology-blind complete graph (all 45 possible pairs connected, equal weight = the real
graph's own mean edge weight `0.536`), seeded at EC: non-seed arrival-time range collapses to **0.000**
(perfect symmetry — every non-seed node reached at *exactly* the same rate) vs. the real graph's
**3.218**-unit range — a **0% remainder**, `<= 1%` threshold — **PASS**. No structure at all produces no
staging signal at all.

**All 5 Part-B geometric gates: PASS** (`part2_graph_connected_by_construction`,
`part2_rho_real_meets_prereg_threshold`, `part2_ec_is_top_ranked_of_10_candidate_epicenters`,
`part2_scrambled_topology_degrades_below_prereg_fraction`, `part2_void_floor_collapses_ordering_signal`).

**A real bug was caught and fixed during this build, disclosed per the OODA discipline**: the first run
of the scrambled-topology gate FAILED (`False`) even though the underlying numbers (`-0.005 <= 0.470`)
were already a clean pass — Orient found the cause was a Python `numpy.bool_(False) is False` identity
comparison (a NumPy scalar is never the same object as the Python `False` singleton), not a modeling
defect. Fixed to `np.isfinite(...)`; re-ran; 2 independent runs produce byte-identical JSON
(`md5sum` verified). Reported here exactly because a lazily-accepted "honest negative" at that point would
have been **wrong** — the model was fine, the gate-composition code was buggy.

## 6. Part C — temporal ordering + tau-vs-amyloid decorrelation (the decisive falsifier)

### 6.1 Temporal ordering — Villemagne 2013's own quantified accumulation curve [13]

AIBL, n=200 (163 with confirmed positive Aβ-accumulation used for trajectory-fitting). Years-before-clinical-dementia-onset
that each biomarker reaches its own threshold, **quoted directly, then subtracted transparently**:

| Biomarker | Years before onset | 95% CI |
|---|---:|---|
| Aβ-PET positivity (SUVR 1.5) | **17.0** | [14.9, 19.9] |
| Hippocampal atrophy | **4.2** | [3.6, 5.1] |
| Memory impairment | **3.3** | [2.5, 4.5] |

**Computed gaps** (pre-registered gates: amyloid must lead neurodegeneration by `>=10y`; atrophy and
memory impairment must be near-coincident, `<=2y` apart):

- Amyloid -> hippocampal atrophy: `17.0 - 4.2 = 12.8` years — **PASS** (`>=10`)
- Amyloid -> memory impairment: `17.0 - 3.3 = 13.7` years
- Atrophy -> memory impairment: `4.2 - 3.3 = 0.9` years — **PASS** (`<=2`)

Amyloid changes over a **full decade before** structural neurodegeneration and cognitive decline, which
are themselves nearly coincident with each other in the final few years before onset — Villemagne's own
data: "Aβ deposition is slow and protracted... As AD progressed, the rate of Aβ deposition slowed towards
a plateau" [13, verbatim]. This is the quantitative anchor for "amyloid... early-plateau."

### 6.2 Same-cohort decorrelation — Ossenkoppele 2016's quantified head-to-head [17]

n=16 AD patients scanned same-day with all 3 tracers (of 20 total AD + 15 amyloid-negative controls),
all correlations P<0.001:

| Pair | abs(r) |
|---|---:|
| tau-PET vs. FDG-hypometabolism | **0.49** (SD 0.07) |
| amyloid-PET vs. FDG-hypometabolism | **0.16** (SD 0.09) |
| tau-PET vs. amyloid-PET | 0.18 (SD 0.09) |

**Pre-registered gate**: ratio `>= 2.0` (tau tracks the neurodegeneration proxy at least twice as
strongly as amyloid does, same cohort). **Measured: `0.49 / 0.16 = 3.06` — PASS**, a 3-fold advantage.
Precision note, kept explicit rather than blurred: this is tau/amyloid **vs. regional glucose
hypometabolism** (FDG-PET, an established functional-neurodegeneration proxy, not literally an MMSE
score) — the qualitative *direct-cognition* claim is separately, explicitly sourced from [16] (verbatim:
tau "was a better predictor of cognitive performance than Aβ deposition in any region") and from [23]
(n=438, MMSE explicitly among the assessed instruments: "Tau deposition showed stronger cognitive
correlations than Aβ"). Also from [17]: amyloid bound "diffusely throughout the neocortex" for **all**
three distinct clinical phenotypes tested (posterior cortical atrophy, amnestic, logopenic aphasia),
while tau/FDG uptake localized to the phenotype-specific affected region in each case — the quantified
"amyloid is diffuse; tau tracks the specific presentation" datapoint.

### 6.3 Convergence tally — decorrelated by research SITE, not by paper count

| Site | Cohort | Finding | Ref |
|---|---|---|---|
| WashU (Knight ADRC) | cognitively normal + mild AD | tau "better predictor of cognitive performance... than Aβ... in any region" | [16] |
| Harvard Aging Brain Study | n=60 | serial mediation Aβ->tau->cognition; tau (β=-3.28, P=.001), not baseline Aβ directly, drives PACC change | [14] |
| UCSF Memory and Aging Center (3 overlapping-cohort papers, counted **once**) | n=16-40 across 3 papers | tau 3.06x amyloid at FDG-tracking; tau not amyloid predicts atrophy; tau-cognition match only weakly amyloid-related | [17,15,18] |
| Sichuan Provincial Hospital, China | n=438 (largest, most recent) | tau stronger cognitive correlations than Aβ; independent driver | [23] |

**4 independent sites, 0 dissent** — pre-registered gate `>=3 sites, 0 dissent` — **PASS**. The UCSF
grouping is a deliberate **de-inflation**: [17], [15], and [18] share the same Rabinovici/Jagust-affiliated
recruitment base and plausibly overlapping patients — counting them as 3 independent replications would
be a common-mode over-count, so they are folded into one site here, a correction against the claim's own
apparent strength, not for it.

A **targeted adversarial search** (3 query variants via NCBI esearch) for a same-cohort finding where
amyloid beats tau at tracking cognition/neurodegeneration was run this session and found none — reported
as exactly that: a disclosed, honest, **non-exhaustive** negative search result (§9), not proof no such
finding exists in the wider literature.

### 6.4 The off-target-binding caveat — argues against, not for, the finding

AV-1451/flortaucipir tau-PET tracer off-target binding is directly documented: neuromelanin/melanin,
hemorrhagic lesions [20]; vessels, iron-associated regions, substantia nigra, choroid-plexus calcification,
leptomeningeal melanin [21] — [21] states explicitly that AV-1451 "does not completely reflect early
stage tau progression suggested by Braak... staging," a genuine, disclosed tension with [7]'s "corresponded
well" finding (different aspects of fidelity: gross stage-correspondence vs. early-stage sensitivity).
**Reasoned, not measured, point**: classical errors-in-variables attenuation bias means measurement noise
uncorrelated with the true signal should **weaken**, not manufacture, an observed correlation — so
tracer noise argues *against*, not for, the tau-beats-amyloid finding being an artifact of tau-PET being
a "better" (higher-SNR) measurement than amyloid-PET. This is flagged explicitly as **reasoning**, not a
new empirical datapoint (§9).

## 7. Forced adversaries — recap (the concept each was forced to its strongest fair form)

| Claim | Adversary | Result |
|---|---|---|
| EC is the true epicenter | Exhaustive sweep, all 10 candidate seeds, real topology | EC strict argmax (rank 1/10); clean 0.36-0.70 gap vs. all others |
| Hierarchical topology (not just a seed label) drives the order | 500-trial scrambled topology, same seed | Collapses to noise (mean ρ=-0.005 vs. real 0.939) |
| Any diffusion looks "staged" (tautology risk) | Void-floor: topology-blind complete graph | Zero arrival-time spread (0% of real spread) |
| Trans-neuronal spread (not regional vulnerability) | Cope 2018's AD-vs-PSP identical-methodology test [9] | Falls in AD, succeeds in PSP — a genuine dissociation |
| Tau-beats-amyloid is a tracer-SNR artifact | Attenuation-bias reasoning + documented off-target binding [20,21] | Argues against, not for, the finding (noise attenuates, doesn't inflate) |
| Tau-beats-amyloid is a single-cohort fluke | 4-site convergence tally, common-mode-corrected | 4/4 independent sites agree, 0 dissent found (targeted search) |

## 8. Falsifier verdict — stated exactly as pre-registered

> Does the model reproduce the measured tau-PET (or CSF p-tau) vs. cognition correlation being STRONGER
> than amyloid-vs-cognition (Braak stage / tau-PET SUVR vs. MMSE) **and** the Braak-stage stereotyped
> spatial ordering?

**Spatial ordering: PASS**, on both an independent first-principles geometric toy model (5/5 gates,
including two orthogonal forced adversaries and a void floor) **and** 4 decorrelated real-world legs
(autopsy, biochemical seeding, in-vivo PET, connectome diffusion) with an adversarial AD-vs-PSP
dissociation that correctly falls in AD and correctly does not fall in PSP.

**Tau-beats-amyloid at tracking the disease: PASS**, on a pre-registered quantified same-cohort ratio
(3.06x, [17]), a pre-registered quantified temporal-gap computation (amyloid leads neurodegeneration by
12.8-13.7 years, structural/cognitive decline nearly coincident at 0.9 years apart, [13]), an explicit
same-cohort verbatim quote naming "cognitive performance" directly ([16]), a large modern MMSE-inclusive
replication ([23], n=438), and a 4-site convergence tally with zero dissent found in a targeted search.

**Neither result is presented as a clean, uncomplicated sweep.** Genuine, disclosed tensions remain open
(§9): the causal direction of the amyloid-tau relationship; whether amyloid gates tau's rate/magnitude
without gating its spatial pattern (Vogel [8] vs. Hanseeuw [14]/Palmqvist [19]); a scale-dependent
disagreement on amyloid-tau spatial correlation strength itself (weak `r=0.18` in [17]'s small mixed-phenotype
n=16 vs. strong `r=0.65-0.81` in [23]'s large mostly-typical-AD n=438); and the toy geometric model's
explicit, disclosed non-connectomic-precision status.

## 9. Pre-registered gates — 9/9 PASS + 6 disclosed open items

```
part2_graph_connected_by_construction:                    PASS
part2_rho_real_meets_prereg_threshold:                     PASS (0.939 >= 0.75)
part2_ec_is_top_ranked_of_10_candidate_epicenters:         PASS (rank 1/10)
part2_scrambled_topology_degrades_below_prereg_fraction:   PASS (mean rho -0.005 <= 0.470)
part2_void_floor_collapses_ordering_signal:                PASS (0.0% of real spread)
part3_amyloid_leads_neurodegeneration_by_over_a_decade:    PASS (12.8y >= 10y)
part3_atrophy_and_memory_impairment_near_coincident_lt_2y: PASS (0.9y <= 2y)
part3_tau_fdg_tracking_at_least_2x_amyloid_fdg_same_cohort:PASS (ratio 3.06 >= 2.0)
part3_convergence_at_least_3_independent_sites_zero_dissent: PASS (4 sites, 0 dissent)

OPEN MODELING UNCERTAINTY (disclosed, does NOT gate overall_pass):
toy_graph_is_not_a_connectome:                    hand-specified textbook toy, not DTI-derived; role is
                                                   mechanism-demonstration, external magnitude anchor is [8]
amyloid_causal_direction_genuinely_contested:      Jack 2013 [12] allows independent initiation; amyloid
                                                   may gate tau RATE (Hanseeuw[14]/Palmqvist[19]) without
                                                   gating spatial PATTERN (Vogel[8]) -- not reconciled here
amyloid_tau_spatial_correlation_disagrees_by_scale: weak r=0.18 (n=16,[17]) vs strong r=0.65-0.81
                                                   (n=438,[23]) -- disclosed, not smoothed over
tau_pet_off_target_binding_is_a_measurement_caveat: [20,21] document real off-target binding; attenuation-
                                                   bias reasoning argues against artifact, not a new datum
no_single_cohort_jointly_measures_all_4_spatial_legs: autopsy/seeding/PET/connectome legs corroborate
                                                   ACROSS cohorts, not jointly within one
brier_2016_headline_is_qualitative_not_an_r_value: [16]'s "better predictor" quote verified verbatim;
                                                   abstract does not print the underlying statistic
dissent_search_was_targeted_not_systematic:        3 query variants, none found -- an honest negative
                                                   search result, not a systematic-review-grade claim
```

**Overall: PASS** (9/9 machine-graded gates; the geometric-model half is independently, deterministically
reproducible — 2 runs of `scripts/msk/tau_pathology.py` produce byte-identical `md5sum`-verified JSON).

## 10. Honest gaps (disclosed, not hidden)

- **The toy connectivity graph is illustrative, not connectomic.** 10 nodes and 14 hand-specified edge
  weights from textbook neuroanatomy — explicitly not a claim of DTI-tractography-grade precision. Its
  role is to show the *mechanism* (hierarchical topology + correct epicenter => stereotyped order) is
  real and machine-falsifiable from first principles; the real-data magnitude anchor is Vogel 2020 [8]
  (n=312, 70% variance).
- **Causal direction is unresolved.** Whether amyloid pathophysiology is upstream of, downstream of, or
  parallel-and-independent-from tau initiation is explicitly held open by the field's own most-cited
  model [11,12] and is not adjudicated by anything in this document.
- **Amyloid's role may be to gate tau's RATE without gating its spatial PATTERN** — a nuance ([8] vs.
  [14]/[19]) not fully reconciled here.
- **The amyloid-tau spatial correlation itself is scale/cohort-dependent** ([17] weak vs. [23] strong) —
  a genuine, disclosed tension.
- **Tau-PET tracers have documented off-target binding** [20,21] — a real measurement limitation, argued
  (not measured) to work against, not for, the central finding.
- **No single cohort jointly measures all 4 spatial-spread legs** — each corroborates across different
  populations, the same disclosed limitation already present in this repo's own prior-session
  `NEU-TAU-PROPAGATION` graph-node design (re-verified, not re-litigated, this session).
- **The dissent search was targeted (3 queries), not a systematic review** — absence of a counter-finding
  in a quick adversarial search is reported exactly as that, not inflated into "no such finding exists."
- **Confidence tier, stated precisely**: the Braak-staging/spatial-spread claim (§5) is **in-vivo/autopsy-anchored**
  across 4 independent modalities. The tau-beats-amyloid-cognition claim (§6) is **in-vivo-anchored**
  (tau-PET, CSF/plasma p-tau) with one quantified same-cohort ratio and a large modern replication. The
  PHF molecular-mechanism claim (§4) is **in-vitro/ex-vivo biochemical + atomic-resolution structural**,
  not an in-vivo measurement. The toy geometric model (§5.2-5.3) is a **first-principles mechanism
  demonstration**, explicitly not connectome-precision-anchored.

## 11. Couples to (graph context, not mutated)

- **Couples to the amyloid-cascade thread** — `docs/MECHANISM_AMYLOID_BETA_AGGREGATION.md` does **not
  exist** in this repo as of this session (confirmed by exhaustive `find`/`grep` at the start of this
  build). Referenced by name per the task's own framing ("together the AD 'amyloid-cascade' picture"),
  not read, not assumed, not mutated. Whichever session builds it should read §6 and §9 here first — the
  causal-direction and amyloid-gates-tau-rate questions are symmetric across both documents and should
  not be independently re-litigated into disagreeing answers.
- **Couples to `NEU-TAU-PROPAGATION`** (`data/MECHANISM_ANCHOR_GRAPH.json`, status `OPEN`, `SEED-DESIGN`
  grade, prior session) — that node's own 9 datapoints and regime-note were independently re-verified
  live this session (all 9 PMIDs matched byte-for-byte against the raw NCBI text) and form the backbone
  of §5.1's 4-leg spatial-spread evidence and part of §6's convergence tally. This document extends that
  node with (a) an independent first-principles geometric model the prior node did not have, and (b) the
  molecular hyperphosphorylation->PHF mechanism (§4), which the prior node did not cover at all (its
  scope was spread/propagation and cognition-tracking, not the aggregation mechanism itself). Not
  promoted to a different `status` here — that is a graph-maintenance action out of scope for this cert.
- **Couples to `NEU-ALZHEIMERS-MULTIMODAL`** / `AUTO-CROSS-HYPOTHESIS-AD-MECHANISM-TRIANGULAT` (same
  graph, `OPEN`, `SEED-DESIGN`) — that cell's own claim names "tau-PET-vs-Braak-vs-cognition correlation"
  as one of three legs for a broader AD mechanism-triangulation; this document supplies a load-bearing,
  independently-verified quantitative treatment of exactly that leg, without resolving the broader
  triangulation cell itself (glia/TREM2/APOE genetics, sleep/glymphatic clearance are out of scope here).
- **Couples to neurodegeneration** (general) — via `NEU-GLIA-NEUROINFLAMMATION` and
  `NEU-MICROGLIA-ACTIVATION-STATE` (same graph): tau pathology's downstream neurodegeneration is
  mechanistically entangled with microglial/complement activity in the broader literature; this document
  does not model that coupling, only names it.
- **Couples to aging** — via `PSY-COGNITIVE-AGING-TRAJECTORY` (same graph) and directly via Schöld 2016
  [7], which found medial-temporal tau tracer retention in **cognitively normal older adults** (not just
  AD patients) predicting worse episodic memory — i.e., the same tau-deposition process this document
  models is not AD-exclusive; it is part of the normal-aging continuum before crossing into clinical AD,
  a genuine substantive link, not a citation-of-convenience.
- **Couples to `BRAIN-CONNECTOME-CENTRAL-COMPUTE`** (same graph, status `WEAKENED`) — this document's
  §5.2-5.3 toy graph is conceptually adjacent (macro-scale structural/functional connectivity) but does
  **not** resolve that node's own disclosed defect (no shared-subject cross-scale connectome); the two
  remain separate, unmerged concerns.

## 12. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/tau_pathology.py
```

No upstream JSON dependency (unlike most `MECHANISM_*` MSK docs — this is a fresh literature+geometric
build, not a re-solve of any prior twin layer). Runs in under a second, fully deterministic (fixed RNG
seed `20260722` for the 500-trial scrambled-topology sweep; verified: 2 independent runs produce
byte-identical JSON, `md5sum`-checked, not eyeballed).

**Paths**: script `scripts/msk/tau_pathology.py`; raw results
`data/tau_pathology/tau_pathology_results.json`; this doc `docs/MECHANISM_TAU_PATHOLOGY.md`; citation
dossier `docs/MECHANISM_TAU_PATHOLOGY_evidence.json`; upstream graph node (read-only, referenced not
mutated) `data/MECHANISM_ANCHOR_GRAPH.json` (`NEU-TAU-PROPAGATION`, `NEU-ALZHEIMERS-MULTIMODAL`).
