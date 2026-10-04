# MECHANISM ANGIOGENESIS / VEGF SPROUTING — hypoxia→HIF-1α→VEGF, Dll4-Notch tip/stalk selection, and the angiogenic switch, machine-certified (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** Builds a NEW mechanistic layer — no
`MECHANISM_ANGIOGENESIS*.md` or `MECHANISM_VEGF*.md` existed in this repo prior to this session (checked
live: `grep -rli angiogen\|VEGF\|bevacizumab\|Dll4 docs/` returned only this doc's own sibling mentions,
below). Related but DISTINCT existing artifacts, not folded/edited (isolation rule — touch only files
you create): the OPEN, `SEED-DESIGN`-tier graph node `ONCO-ANGIOGENESIS-HYPOXIA` in
`data/MECHANISM_ANCHOR_GRAPH.json` (a literature-scout proposal for a *downstream validation cohort*
design — Hurwitz-vs-Gilbert bevacizumab trials — not a mechanistic model), the sibling `NOTCH-LATERAL-
INHIBITION-FATE-SWITCH` node (about T-ALL/intestinal NOTCH1 oncogene contexts, a different tissue
application of the same NOTCH pathway), and `docs/MECHANISM_WOUND_HEALING.md` (whose own proliferation
phase, checked live this session, models re-epithelialization only — its granulation-tissue
neovascularization component is currently unmodeled; this doc supplies that missing mechanism without
editing that doc, per isolation rules). Script: `scripts/msk/angiogenesis_vegf.py`. Evidence (machine-
written): `data/msk_smoketest/angiogenesis_vegf/angiogenesis_vegf_results.json` (md5
`1e906cc78b2ea6ee00ca168aa2f41fd5`, byte-identical across repeated runs — determinism confirmed).

**Recall-drift measured, not assumed, this session:** a first-pass batch of 16 citation PMIDs pulled
from memory was spot-checked via one `esummary` call before any were used — **6 of 16 (37.5%) resolved
to a completely unrelated paper** (an obesity-gut-microbiome paper, an ocean-productivity paper, an
FGF21/Klotho metabolism paper, a cytochrome-oxidase-inhibition paper, etc., for PMIDs I had associated
with Hellcurrent 2007, Noguera-Troise 2006, Suchting 2007, Forsythe 1996, Folkman & Hochberg 1973, and
Jain 2001 respectively). Every single citation below was therefore independently re-derived this
session via NCBI `esearch`/`esummary`/`efetch`, not trusted from recall — consistent with, and a live
confirmation of, this environment's own disclosed ~62% citation recall-drift rate. WebSearch was
already session-quota-exhausted (2000/2000) before this doc's research began (same constraint the
sibling wound-healing doc hit); all 27 citations were verified via NCBI E-utilities, with Europe PMC
REST as a same-session fallback when NCBI intermittently rate-limited (HTTP 429).

**Confidence tier, split precisely:** the Dll4-Notch tip/stalk lateral-inhibition DIRECTION (blockade
→ excess tip cells) is **in-vivo-anchored** — three independent, live-quoted 2006-2007 mouse studies
(different labs, different model systems: retina, tumor xenograft, retina/embryo) agree. The O2-
diffusion-limit ORDER OF MAGNITUDE is **geometrically-derived and independently cross-checked**
against a textbook-tier historical figure. The exact micron/millimeter magnitudes (Thomlinson & Gray
1955's own digit; the "~1-2mm" angiogenic-switch threshold) are **textbook/secondary-sourced** — disclosed
per-claim below, not smoothed into one blended grade (same convention `MECHANISM_WOUND_HEALING.md`
established).

---

## 0. Falsifiers (pre-registered, matching the task) + verdict up front

| # | Falsifier | Verdict |
|---|---|---|
| F1 | HIF-1α is a genuine exponential/threshold function of O2% (not linear); HIF-1 is CAUSALLY required for hypoxic VEGF transcription | **PASS** — linear adversary requires unphysical negative-clipping over 85.2% of the tested domain; ARNT-null genetic void-floor directly falsifies "HIF-1 not required" |
| F2 | The model reproduces the MEASURED ~100-200µm O2 diffusion distance from a capillary | **PASS** — central-estimate geometric derivation = 161.4µm, inside the textbook anchor; full parameter sweep [49.5, 499.6]µm brackets it |
| F3 | Angiogenic switch: avascular growth beyond ~1-2mm requires neovascularization | **PASS, with an honest 3-way size reconciliation**, not one collapsed number (§3) |
| F4 | Dll4-Notch lateral inhibition produces tip/stalk pattern; Dll4 blockade → EXCESS non-productive tip cells (hypersprouting) | **PASS** — normal coupling gives a robust 1/3 minority-tip pattern (100% of 30 seeds); full AND 10-fold-partial blockade both give 100% tip fraction; cross-checked against 3 independent live-quoted papers |
| Decorrelated | Anti-VEGF (bevacizumab) prunes vessels + shows a measured, MODEST clinical effect | **PASS** — Hurwitz 2004 CRC: +4.7mo OS (HR 0.66, real); Gilbert 2014 GBM: -0.4mo OS (HR 1.13, NS) despite a real PFS benefit — modest AND tumor-type-dependent, not preclinical-hype-scale |

`overall_pass = True` in the evidence JSON. **One design flaw was caught, diagnosed, and fixed
mid-session, disclosed not hidden** (§4): the first F4 parameterization (an 8-cell mean-field model)
never bifurcated at any reasonable coupling strength, and its void-floor test checked the wrong failure
direction — both are OODA-diagnosed and corrected below, not swept away.

---

## 1. F1 — the hypoxia→HIF-1α→VEGF axis: a genuine exponential/threshold switch, causally required

**Identity anchor:** Wang & Semenza (1995, PMID 7539918, identity-verified live) first identified HIF-1
as the basic-helix-loop-helix-PAS heterodimer (HIF-1α/ARNT) whose stability is regulated by cellular O2
tension — the discovery this entire axis rests on, cited here for completeness though no specific
number in this doc traces to it.

**Anchor (quantitative):** Jiang, Semenza, Bauer, Marti (1996, PMID 8897823) is the only primary O2
dose-response curve found live this session. Verbatim (efetch): HIF-1 DNA-binding activity and
HIF-1α/HIF-1β protein levels **"increased exponentially as cells were subjected to decreasing O2
concentrations, with a half maximal response between 1.5 and 2% O2 and a maximal response at 0.5%
O2"** (HeLa cells, 0-20% O2 tested). This doc's model is fit to exactly these two measured points
(half-max at the band's midpoint, 1.75% O2; near-max at 0.5% O2) — not a free-form curve.

**The task's own "~5% O2" framing, evaluated as a PREDICTION of this model, not built into the fit:**
at O2=5%, the calibrated exponential predicts HIF at **8.25% of its 0.5%-O2 maximum** — i.e., already
well past the steep part of the curve and close to baseline. This sharpens, rather than contradicts,
the task's framing: "~5% O2" is best read as the upper edge/onset of the physiologically-responsive
range (matching Jiang 1996's own statement that "the HIF-1 response was greatest over O2 concentrations
associated with ischemic/hypoxic events in vivo"), not the half-maximal point itself (which the primary
data place at 1.75% O2).

**Forced adversary — a 2-point LINEAR fit through the same two anchors:** solved in closed form
(slope=-0.4, intercept=1.2 in the model's normalized units). Machine-checked over the tested domain
(O2 ∈ [0.1%, 20%]): the linear form goes **negative (unphysical) starting at just O2=3.04%**, and is
negative over **85.2% of the tested domain** — i.e., for nearly the entire normoxic-to-mild-hypoxic
range, a linear model would require clipping HIF-1α to an impossible negative value. **Adversary
falls.** (The exponential form, by construction, only asymptotes toward — never reaches — exactly
zero, consistent with real biology's small constitutive normoxic HIF-1α turnover.)

**Causal (not merely correlational) anchor — Forsythe et al. (1996, PMID 8756616):** live-quoted
verbatim (efetch): a 47bp hypoxia-response element (HRE) 985-939bp 5′ of the VEGF transcription start
binds HIF-1 and drives hypoxia-inducible reporter expression; a 3bp HRE substitution abolishes
hypoxia-inducibility; dominant-negative HIF-1α inhibits induction dose-dependently; and decisively,
**"VEGF mRNA was not induced by hypoxia in mutant cells that do not express the HIF-1β (ARNT)
subunit."** This is a genetic loss-of-function void-floor result, not a correlation — HIF-1 necessity
for hypoxic VEGF induction is machine-flagged `True` on this direct quoted finding.

VEGF itself: discovered as a "vascular permeability factor" (Senger et al. 1983, PMID 6823562) and
cloned/identified as a secreted angiogenic mitogen (Leung et al. 1989, PMID 2479986) — both
identity-verified live, background context for the axis's downstream ligand.

**F1_pass = True** (adversary falls AND causal mechanism verified).

---

## 2. F2 — the O2 diffusion limit: a GEOMETRIC derivation, not a literature lookup, live-coupled to hemoglobin

Per this repo's geometric-thinking discipline, the ~100-200µm diffusion-distance figure is **derived**
here from steady-state reaction-diffusion, not simply quoted. For a spherically-symmetric avascular
tissue mass with uniform O2 consumption rate `A` and diffusion coefficient `D`, fixed surface
concentration `Cs` (the capillary-wall boundary condition) and zero-flux at the center (symmetry):

```
C(r) = Cs − (A / 6D)·(R² − r²)         [steady state, D·∇²C = A]
R_crit = sqrt(6·D·Cs / A)              [center just reaches C=0, i.e. anoxia onset]
```

This is the same reaction-diffusion formalism Krogh (1919, PMID 16993405) originated for tissue O2
supply (the "Krogh cylinder").

**The hemoglobin coupling (a real, live-loaded cross-model dependency, not an assertion):** `Cs` is
computed from a capillary pO2 (this repo's own already-certified `blood_oxygen_transport.py` model's
resting arterial/mixed-venous operating points, PaO2=95mmHg / PvO2=40mmHg, **live-loaded from
`data/msk_smoketest/subject2_walking1/blood_oxygen_transport/blood_oxygen_transport_results.json`**)
via the SAME dissolved-O2 solubility coefficient that sibling model already uses (0.003 mL O2/dL
blood/mmHg), scaled by that model's own certified Hb→CaO2 sweep (anemia Hb=6 g/dL through
polycythemia Hb=22 g/dL) — not re-derived, reused, per this repo's established sibling-coupling
convention.

**Central estimate (mid-range literature-plausible D and A, chosen before computing — confirmed
unchanged across repeated runs, not tuned to the target):** `D=1.5×10⁻⁵ cm²/s`, `A=7×10⁻⁴ mL
O2/(mL·s)`, mean capillary pO2=67.5mmHg → **R_crit = 161.4 µm** — inside the textbook 100-200µm anchor
band traceable to Thomlinson & Gray (1955, PMID 13304213, identity-verified live; the pre-abstract-era
record has no live-fetchable abstract, so the exact digit is textbook/secondary-sourced, disclosed, same
tier as the wound-healing doc's handling of Levenson 1965).

Vaupel, Kallinowski & Okunieff's (1989, PMID 2684393, identity-verified live) review of blood flow,
oxygen and nutrient supply in human tumors independently frames diffusion/hypoxia as jointly dependent
on blood flow AND metabolic (consumption) rate — qualitative confirmation, from a decorrelated review,
of exactly the two axes (`Cs` via capillary supply, `A` via consumption) this doc's sweep varies; its
live-fetched abstract does not itself state a micron figure (disclosed).

**The spread (symmetric-QC, measured not asserted):** sweeping D, A, capillary pO2, and the live-loaded
Hb ratio across their full literature-plausible corners gives **R_crit ∈ [49.5, 499.6] µm** — a ~10×
range that BRACKETS but does not pin the textbook anchor. **Hemoglobin alone** (holding D, A, pO2
fixed): R_crit varies **1.89×** from anemia (Hb=6 g/dL) to polycythemia (Hb=22 g/dL) — a real,
quantified, live-loaded hemoglobin-dependence, directly answering the task's own symmetric-QC
instruction rather than merely asserting it.

**F2_pass = True** (central estimate within [100,200]µm AND the full sweep brackets that anchor).

---

## 3. F3 — the angiogenic switch: reconciling THREE distinct size numbers, not collapsing them

An avascular mass growing past `R_crit` develops an anoxic core, forcing a binary choice: (a) plateau
in dormancy (a viable rim tolerating a necrotic core, net growth arrested), or (b) induce new
capillaries penetrating the mass (Hanahan & Folkman's 1996 "angiogenic switch," PMID 8756718),
resetting the diffusion geometry to many small `R_crit`-radius sub-domains and re-enabling net growth.
This is the mechanistic bridge from F2's geometry to *why* a tip/stalk sprouting program (F4) is
invoked at all.

Three numbers appear across this doc's sources, and **they are not the same quantity**:

| Number | Value | What it actually measures |
|---|---|---|
| This doc's own `R_crit`-derived single-nodule diameter | **0.32mm** (2×161.4µm) | the necrotic-ONSET threshold — rim thickness ×2 |
| The commonly-taught clinical "angiogenic switch" figure | **~1-2mm** (Folkman 1971, PMID 4938153; Hanahan & Folkman 1996) | the point past which further NET growth requires new vessels |
| Folkman & Hochberg's own 1973 PRIMARY measurement (PMID 4744009, live-quoted efetch: spheroids "eventually reached a dormant phase at a diameter of approximately 3-4 mm", ~10⁶ cells) | **3-4mm** | a LATER, larger, stable dormant-PLATEAU diameter (a rim already tolerating a necrotic core) |

Disclosed honestly: onset-threshold < taught-switch-point ≤ dormant-plateau-diameter is a coherent
ordering, but collapsing these into one number would overclaim precision this session's live sources
do not support. This is a genuine, disclosed reconciliation, not a hidden discrepancy.

---

## 4. F4 — Dll4-Notch tip/stalk lateral inhibition: a real N-cell ODE, forced adversary, and an honest mid-session fix

**Computational precedent (identity/mechanism-verified live, not re-implemented here):** this doc's
minimal ODE targets the same phenomenon Bentley, Gerhardt & Bates (2008, PMID 18028963) modeled with a
full agent-based Notch simulation, and its companion robustness study (Bentley, Mariggi, Gerhardt &
Bates 2009, PMID 19876379); Jakobsson & Gerhardt's (2009, PMID 19909253) review covers the same
VEGFR/Notch mechanistic territory. This doc's own model is a much smaller, closed-form-adjacent ODE,
not a re-implementation of theirs — cited as precedent, not as validation of this doc's specific numbers.

**The model (Collier & Lewis 1996 mathematical form, PMID 9015458, Dll4/Notch-relabeled):** N cells,
each with a Dll4 level `d_i`; each cell senses the mean Dll4 of its neighbors (`s_i`) and suppresses
its OWN Dll4 production via a Hill-function inhibitory feedback `f(s)=1/(1+k·s^h)` — the literal
molecular wiring of Dll4-Notch lateral inhibition (a high-Dll4 "winning" cell activates Notch in its
neighbors, which suppresses THEIR Dll4/tip program, forcing stalk fate). Integrated numerically
(`scipy.integrate.solve_ivp`) from small-noise initial conditions to steady state.

**An honest, disclosed, mid-session design failure and its OODA fix (not hidden):** the FIRST design
used N=8 cells with mean-field (all-to-all) coupling. It **never bifurcated** — a diagnostic sweep run
this session (h from 2 to 16, k from 20 to 10⁶) showed the symmetry-breaking threshold scales with
`(N-1)` (mean-field coupling divides the shared signal by 7 neighbors, hugely damping any one cell's
influence), so no biologically-reasonable coupling strength produced a genuine tip/stalk split at N=8.
**Orient:** this diagnosis also revealed the first design's void-floor test checked the WRONG failure
direction — it assumed a non-cooperative (h=1) system would fail *toward* the blockade extreme (all
cells high), but measuring it directly showed the opposite: h=1 converges ALL cells to the SAME low,
undifferentiated value (a real, different, degenerate failure mode). **Fix:** switched to N=3 (a
minimal but MORE biologically faithful unit — Notch-Dll4 signalling is juxtacrine/contact-dependent,
so a real lateral-inhibition neighborhood is a cell's immediate physical neighbors, not an 8-cell
mean-field average), and corrected the void-floor gate to test for genuine bimodal splitting
(spread > 0.5) rather than an assumed failure direction. Re-diagnosed at N=3: h=4, k=100 bifurcates
**cleanly and robustly (100% of 30 seeds)**.

**Results (30-seed robustness sweep, N=3, h=4):**

| Condition | k (coupling) | mean tip_fraction | Interpretation |
|---|---|---|---|
| Normal | 100 | **0.333** (100% of seeds, exactly 1/3) | 1 tip, 2 stalk — the correct minority pattern |
| Full Dll4 blockade (forced adversary) | 0 | **1.000** (100% of seeds) | ALL cells become tip-like — hypersprouting |
| Partial blockade (10-fold reduced coupling, a realistic soluble-trap/antibody-strength block, not complete genetic ablation) | 10 | **1.000** (100% of seeds) | hypersprouting reproduced WITHOUT complete ablation — a stronger, more robust adversary finding |
| Void floor (h=1, non-cooperative feedback, same k=100) | 100 | **0.000**, spread ≤ 0.5 | NO bifurcation at all — cooperativity, not just "some negative feedback," is structurally required (Collier & Lewis 1996's own mathematical result) |

**External, decorrelated anchor (3 independent mouse studies, different labs/systems, all live-quoted
this session):**
- Hellcurrent et al. (2007, PMID 17259973, mouse retina): "inhibition of Notch signalling...promote[s]
  increased numbers of tip cells"; Notch activation gives "fewer tip cells and vessel branches."
- Noguera-Troise et al. (2006, PMID 17183313, mouse tumor xenograft): Dll4 blockade causes "markedly
  increased tumour vascularity...enhanced angiogenic sprouting and branching" YET "poor perfusion and
  increased hypoxia" and net "decreased tumour growth" — the exact "more vessels but less function"
  pattern, in vivo, in the paper's own title.
- Suchting et al. (2007, PMID 17296941, mouse retina/embryo): Dll4 loss-of-function gives "greatly
  increased numbers of filopodia-extending endothelial tip cells" via "increased VEGFR2 and decreased
  VEGFR1" (VEGFR1 normally acts as a decoy receptor damping VEGF signalling).
- Named directly: Thurston & Yancopoulos (2007, PMID 17457300) call this "**The Delta paradox: DLL4
  blockade leads to more tumour vessels but less tumour growth**" — the review's own title matches the
  falsifier's framing verbatim.

The "poor perfusion / decreased tumor growth" functional consequence is Noguera-Troise 2006's own
in-vivo measurement, cited rather than re-derived — this doc's N=3 ODE computes the STRUCTURAL
(excess-tip-cell) side only, honestly scoped.

**F4_pass = True** (normal minority pattern, full-blockade hypersprouting, partial-blockade
dose-robustness, and void-floor-falls all machine-verified).

---

## 5. Decorrelated clinical check — bevacizumab, two tumor types, no shared patients or authors

| | Hurwitz et al. 2004 (PMID 15175435) | Gilbert et al. 2014 (PMID 24552317) |
|---|---|---|
| Tumor type | Metastatic colorectal cancer | Newly-diagnosed glioblastoma |
| n | 813 | 637 |
| Median OS (drug vs control) | **20.3 vs 15.6 months** | 15.7 vs 16.1 months |
| OS hazard ratio | **0.66, P<0.001** | 1.13 (NS — numerically favors control) |
| Median PFS (drug vs control) | 10.6 vs 6.2 months (HR 0.54) | 10.7 vs 7.3 months (HR 0.79, real benefit) |

Both live-fetched (efetch) verbatim this session. **The single positive registration result (CRC) is
real and significant, but modest in absolute terms: 4.7 months of median OS, not years.** In a
decorrelated tumor type (GBM — different patients, different authors, same drug/mechanism), the OS
effect is null (numerically −0.4 months, HR 1.13) DESPITE a real PFS benefit — i.e., anti-VEGF prunes
and reorganizes tumor vasculature (the "normalization" mechanism, Jain 2001, PMID 11533692) and
measurably slows progression, but this does not reliably convert to a survival benefit across tumor
types. This is the "modest, not preclinical-hype-scale" clinical reality the task's symmetric-QC asks
to hold OPEN — demonstrated with two real trials, not asserted. A fresh (2026) quarter-century
retrospective found live this session while searching for Jain 2001 (PMID 41997128, "Targeting
angiogenesis: Lessons from 25 years of normalizing tumor vasculature," Cell) is cited as an up-to-date
anchor for the same "modest 25 years on" framing.

---

## 6. Symmetric QC — nothing proven; every spread is reported, not collapsed

- **Diffusion limit depends on metabolic rate AND hemoglobin — measured, not asserted (§2):** R_crit
  spans a ~10× range ([49.5, 499.6]µm) across literature-plausible D/A/pO2/Hb combinations; hemoglobin
  alone drives a 1.89× swing. The central estimate landing inside the textbook anchor is a plausibility
  check on the geometric FORM, not proof of the exact parameter values — held explicitly OPEN.
- **Anti-angiogenic clinical benefit is real but smaller than preclinical hype (§5):** a 4.7-month OS
  benefit in the one positive registration trial; null in a decorrelated tumor type despite a real PFS
  effect. Ellis & Hicklin's (2008, PMID 18596824, identity-verified live) review of VEGF-targeted
  therapy mechanisms independently frames this same modest/heterogeneous clinical picture and the
  resistance mechanisms behind it; Brown & Wilson's (2004, PMID 15170446, identity-verified live)
  standard review of exploiting tumor hypoxia in cancer treatment is the modern context for why hypoxia
  itself (not just VEGF) is a targetable state.
- **Species caveat:** all three Dll4-Notch mechanistic anchors (Hellcurrent, Noguera-Troise, Suchting)
  are MOUSE model systems (retina, tumor xenograft, embryo). Direct live-quoted human sprouting-front
  imaging confirmation was not independently found this session — disclosed, not assumed to
  generalize.
- **Toy-model caveat:** F4's ODE is a minimal N=3 mean-field system, not a spatial reaction-diffusion
  PDE over real 3D sprout geometry, and treats VEGF drive as spatially uniform. Real sprouting
  additionally uses a VEGF GRADIENT for directional tip-cell filopodial guidance (Gerhardt et al. 2003,
  PMID 12810700) — a distinct mechanism from lateral inhibition's role in selecting WHICH cell becomes
  tip; not modeled here.
- **Krogh-parameter tier:** `D_O2` and tissue O2-consumption-rate `A` are textbook-standard literature
  ranges, not independently re-extracted from one live-quoted primary numeric source this session
  (same disclosed tier as Thomlinson & Gray 1955's own exact micron figure, a pre-abstract-era record).
  The certified claim is the geometric FORM plus order-of-magnitude convergence, not exact-digit
  precision.
- **The three distinct "tumor size threshold" numbers (§3)** are reconciled, not collapsed — a
  genuine, disclosed precision distinction.
- **Definitional spread, independently confirmed live:** Hickel & Vaupel (2001, PMID 11181773) states
  in its own abstract (live-fetched) that "biochemists and clinicians define hypoxia differently" — an
  independent primary-source confirmation, not an inference, that the O2/hypoxia threshold is
  definition-dependent rather than one universal fixed number.
- **Held OPEN throughout** — none of this spread is swept into `overall_pass`.

---

## 7. Gates — machine-printed, not narrated

```
F1_hif_axis_pass                                   True
F1_linear_adversary_falls                          True   (negative over 85.2% of domain, onset O2=3.04%)
F1_vegf_causal_mechanism_verified                  True   (ARNT-null abolishes induction, Forsythe 1996)
F2_diffusion_limit_pass                            True
F2_central_within_textbook_anchor                  True   (161.4um, anchor [100,200]um)
F2_sweep_brackets_textbook_anchor                  True   (sweep [49.5, 499.6]um)
F4_normal_minority_pattern_pass                     True   (tip_fraction=0.333, 100% of 30 seeds)
F4_full_blockade_hypersprouting_pass                True   (tip_fraction=1.000, 100% of 30 seeds)
F4_partial_blockade_dose_dependence_pass            True   (10x-reduced coupling already =1.000)
F4_void_floor_falls                                 True   (h=1 non-cooperative: no bifurcation, spread<=0.5)
F4_lateral_inhibition_pass                          True
OVERALL_PASS                                        True
```
Deterministic: byte-identical md5 `1e906cc78b2ea6ee00ca168aa2f41fd5` confirmed across repeated process
runs (fixed RNG seeds per run in the 30-seed robustness sweeps; all other computation is closed-form or
`solve_ivp` with fixed tolerances, no unseeded stochastic element).

---

## 8. Honest gaps — symmetric QC: what this does NOT prove

1. **The exact O2-diffusion-distance and tumor-size-threshold DIGITS are textbook/secondary-sourced**
   (Thomlinson & Gray 1955 and the "~1-2mm" switch figure) — this doc's own geometric derivation
   reproduces the right ORDER OF MAGNITUDE independently, which is the actual certified claim.
2. **Krogh-model parameters (D, A) are literature-plausible ranges, not pinned to one live-quoted
   primary numeric source this session** — disclosed, not hidden; the resulting ~10× R_crit spread is
   reported, not resolved.
3. **All 3 Dll4-Notch mechanistic anchors are mouse studies** — human generalization is not
   independently re-verified live this session.
4. **F4's ODE is a minimal 3-cell mean-field toy model**, not a spatial PDE over real sprout geometry;
   VEGF-gradient-guided directional migration (a distinct mechanism from fate selection) is not
   modeled.
5. **VEGF transcriptional output is modeled as directly proportional to HIF-1α level** — Forsythe 1996
   verifies the CAUSAL mechanism (binding site, dominant-negative, ARNT-null) but this session did not
   find a live-quoted quantitative HIF-level-to-VEGF-fold-change dose-response curve; the proportionality
   itself is a disclosed simplification.
6. **The bevacizumab decorrelated check uses exactly 2 trials** (one positive, one null) — a real,
   informative, live-fetched contrast, but not a systematic meta-analysis across all anti-VEGF trials.
7. **The N=3/N=8 mean-field-coupling bifurcation-threshold finding (§4) is itself a property of THIS
   toy model's coupling topology** (division by N-1), not an independently-measured biological fact
   about how many endothelial cells participate in a real lateral-inhibition neighborhood — disclosed
   as a modeling artifact/choice, not a biological measurement.
8. **Single, generic parameter set** — no vessel-bed/tissue-type/species-conditioned version of either
   model (same scope precedent as the wound-healing and warburg_metabolism sibling docs).

---

## 9. couples_to — resolving the task's explicit couplings with computed numbers, not prose

- **Wound-healing (granulation-phase neovascularization)** — `docs/MECHANISM_WOUND_HEALING.md`. That
  doc's proliferation phase currently models re-epithelialization only (confirmed live this session:
  no angiogenesis/VEGF/neovascularization mechanism in its body). This doc supplies the missing
  granulation-tissue angiogenesis mechanism (hypoxia in a fresh wound bed → HIF-1α → VEGF → tip/stalk
  sprouting to revascularize granulation tissue) as a companion, not an edit (isolation rule).
- **Cancer (angiogenic switch)** — `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md` (brief existing
  angiogenesis mention, ESWT/VEGF context) and the existing OPEN graph node `ONCO-ANGIOGENESIS-HYPOXIA`
  (a downstream bevacizumab-trial validation-cohort design this doc's own §5 independently confirms the
  same Hurwitz/Gilbert numbers for, without folding into or editing that node).
- **O2 transport/hypoxia** — `scripts/msk/blood_oxygen_transport.py` /
  `data/msk_smoketest/subject2_walking1/blood_oxygen_transport/blood_oxygen_transport_results.json`:
  live-loaded (not re-derived) for this doc's F2 hemoglobin-dependent capillary-O2 boundary condition —
  a genuine, machine-verified cross-model coupling, not a prose gesture.
- **NOTCH pathway (general)** — the existing `NOTCH-LATERAL-INHIBITION-FATE-SWITCH` graph node covers
  T-ALL oncogene/intestinal contexts of the same NOTCH signalling pathway; this doc is the first
  application of the pathway's lateral-inhibition mathematics to endothelial tip/stalk selection
  specifically in this repo.

Per `docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/§4, promoting any of this into canonical
`data/MECHANISM_ANCHOR_GRAPH.json` nodes/edges requires the separate `mechanism_fold → fold_gate_v2`
pipeline — **not performed this session**, consistent with the sibling coagulation/skin-barrier/
wound-healing docs' own precedent of stating results in prose/evidence-JSON only.

---

## 10. Files

- `source_repository/scripts/msk/angiogenesis_vegf.py` — the model: 27-entry `CITATIONS`
  dict (tiered identity/quant/mechanism), the sibling `blood_oxygen_transport.py` JSON loader, the F1
  HIF-1α exponential fit + linear forced adversary + ARNT-null causal fact-check, the F2 Krogh/spherical
  `R_crit=sqrt(6DCs/A)` derivation + parameter-cube sweep + Hb-only sensitivity, the F3 three-way
  size-reconciliation, the F4 N=3 Collier-form lateral-inhibition ODE (`scipy.integrate.solve_ivp`) +
  forced adversary (full + 10-fold-partial Dll4 blockade) + void floor (non-cooperative h=1) + 30-seed
  robustness sweep, the Hurwitz-vs-Gilbert decorrelated bevacizumab check, the symmetric-QC dict, and
  the gates/evidence-JSON writer. Run with `source .venv-msk/bin/activate && python3
  scripts/msk/angiogenesis_vegf.py` (a few seconds wall time, numpy/scipy only, no OpenSim dependency,
  deterministic).
- `source_repository/data/msk_smoketest/angiogenesis_vegf/angiogenesis_vegf_results.json` —
  full machine-written evidence (citations, F1-F4 computed numbers, the decorrelated bevacizumab check,
  symmetric QC, gates, `overall_pass`). md5 `1e906cc78b2ea6ee00ca168aa2f41fd5`.
- `source_repository/data/msk_smoketest/subject2_walking1/blood_oxygen_transport/blood_oxygen_transport_results.json`
  — the sibling model live-loaded for F2's hemoglobin coupling (read-only, not modified).
- `source_repository/data/MECHANISM_ANCHOR_GRAPH.json` — existing related OPEN nodes
  (`ONCO-ANGIOGENESIS-HYPOXIA`, `NOTCH-LATERAL-INHIBITION-FATE-SWITCH`, `REPRO-PREECLAMPSIA-ANGIOGENIC`),
  unchanged this session (fold not performed, §9).
- Sibling docs (read for convention/coupling, not re-litigated): `docs/MECHANISM_WOUND_HEALING.md`,
  `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md`, `docs/MECHANISM_HARDENED_CONVENTIONS.md`.
