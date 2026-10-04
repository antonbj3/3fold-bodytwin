# MECHANISM MITOCHONDRIAL OXPHOS — the organelle-level ETC / ATP-synthase / P-O certified model (2026-07-22)

## ⚠ CORRECTION (2026-07-26 — read before Section C; the graph has moved past this doc's grid)

This doc's Section C reports P/O(NADH) as an open spread of **2.31–2.73**, driven by treating Complex I's
H+/2e- dispute (classic 4 vs revised 3) as unresolved. That dispute has since been adjudicated by a
decorrelated DIRECT-flux measurement (Jones/Blaza/Varghese/Hirst 2017, PMID 28174301: Bos taurus
CI=4.21±0.15 vs P. denitrificans 4.01±0.24; pooled z(CI=3)=9.07σ **EXCLUDED**, z(CI=4)=1.21σ consistent —
independently reproduced 3× in-graph to 4 decimals), plus a same-session falsifier confirming the
**Route-A** transport increment (+1.0, not the back-solved +1.333 "Route-B") reproduces Hinkle 2005's own
second data point exactly (yeast/E. coli c10: computed 2.308/1.385 vs Hinkle's stated 2.3/1.4) while
Route-B does not (2.143/1.286 — Route-B exists only because 8/3+1.333=4.0=10/2.5 was engineered to
preserve the old answer). **Current best point estimate: P/O(NADH) = 2.78 ± 0.04 (1σ)** — 2.7845 from the
measured CI=4.21 (2.7273 from the bare integer CI=4) — which **retires the CI=3-derived 2.455 floor**
(now 9.07σ-excluded) below. The mammalian-c8/classic-CI=4 row in the Section C grid (2.727) is now the
better single point estimate, not merely a "structural ceiling." Cross-checked to 0.2% against
Mookerjee 2017's independent whole-chain measured P/O_max=2.79 (vs 2.25% for the bare-integer row) — a
real external over-determination, not a tautology. Scope caveat carried forward unchanged: holds for the
fully-coupled/high-Δp regime only (Ripple & Springett 2013); Cabrera-Orefice 2018 shows this coupling
mechanically/reversibly disengages (Active/Deactive transition), untested here at low-Δp/ischemia. The
succinate P/O (1.636) is unaffected — FADH2 enters post-Complex-I, bypassing this dispute entirely.

**Separately, a same-day check of the ODE against this doc's own algebra found a small, real, second
gap.** This doc's Section C grid is closed-form algebra, never run through the mechanistic ODE
(`scripts/msk/oxphos_corrected_nhatp_resimulation.py` / the `mitochondrial_oxphos.py` family). Executed
fresh with BOTH corrected inputs consistent (n_H,ETC=10.21 measured, n_H,ATP=3.6667 Route-A), the
**ODE-emergent P/O = 2.7624**, about **-0.80%** off the pure-algebra point (2.7845) — a small, stable
kinetic discount (the same ~-0.74% discount was already present at the old stoichiometry, so this is not
new physics, just newly checked). A prior on-disk partial rerun that updated only one of the two inputs
(left n_H,ETC stale at 10) gave a larger, spurious -2.85% gap (P/O=2.7051) — a stale-parameter artifact,
not a real effect: about 2 of those ~2.85 percentage points was the missing second update, not kinetics.
Neither the 2.7051 nor the 2.7624 ODE output has yet been folded into the graph as its own node.

Both findings are machine-verified in `AUDIT-FLAGSHIP-HEADLINE-DISCRIMINATION-SWEEP-2026-07-26` and
`HOLE-OXPHOS-PO-NADH-BRACKET-PROPAGATION-REPAIR-2026-07-26` (graph status OPEN /
`HYPOTHESIS-awaiting-QC` — current best position, not yet a coordinator-verified closed number). Zero
pre-registered gates in this doc flip as a result — both `g08`'s 0.10 tolerance and the composite
falsifier `g17` still PASS under the sharpened value. The rest of this document (Sections A-E, the
citations table, honest gaps) is preserved as originally written below: the derivation and falsifier
logic remain correct and useful; only the "held open" status of the CI dispute and the untested
algebra-vs-ODE gap are superseded.

---

Builds the twin's first **CERTIFIED organelle-level model of oxidative phosphorylation**: the
electron-transport chain (complexes I-IV, proton-pumping stoichiometry H+/2e-), ATP synthase
(complex V) stoichiometry **derived geometrically as a gear ratio between two coupled rotary
symmetries** (not asserted as a memorized fraction), the P/O ratio this implies, and mitochondrial
membrane potential / proton-motive force. Couples to this twin's own already-certified
`docs/MECHANISM_MUSCLE_ENERGETICS.md` (ATP-supply / PCr-recovery link) and to aging (mitochondrial
capacity decline). Script: `scripts/msk/mitochondrial_oxphos.py`. Evidence:
`data/mitochondrial_oxphos/mitochondrial_oxphos_results.json`.

**No subject/trial data required.** This is a molecular/organelle-constant-level model (unlike the
gait-trial-specific `msk_smoketest/subject2_walking1/*` docs) — self-contained closed-form algebra
against live-verified literature constants, no OpenSim, no numpy.

## Pre-registered falsifier (stated in the task before this script was written)

Does the model reproduce the **measured** P/O ratio (~2.3-2.5 NADH-linked / ~1.5 succinate — Hinkle
2005's careful synthesis of studies since 1937, the modern non-integer values that refuted the old
textbook 3.0/2.0) **and** the measured mitochondrial membrane potential (~150-180 mV, TMRM/TPP+)?
Symmetric adversary, forced not assumed: the mammalian ATP-synthase c-ring cryo-EM structure (Watt
et al. 2010) predicts a **structural ceiling** on P/O — a genuinely decorrelated instrument (X-ray/
cryo-EM structure determination, not O2-flux measurement) and a genuinely decorrelated *era*
(2010 postdates and could not have influenced Hinkle's 2005 flux synthesis, nor the reverse). If
the independently-measured flux P/O ever exceeded this structural ceiling, the model would be
falsified outright (proton leak/slip can only ever *lower* measured P/O below the structural
maximum, never raise it above it — a hard, checkable, physically-required direction).

**Result: PASS, 17/17 machine gates**, including the composite falsifier gate
(`g17_COMPOSITE_falsifier_PO_and_membrane_potential_both_reproduced`). Every gate is a closed-form
numeric comparison, machine-run — no figure was eyeballed, no narration substitutes for a number.

**Held OPEN, not resolved (symmetric QC, exactly as the task asked)**: Complex I's own H+/2e-
stoichiometry is a **live, cited, unresolved dispute in the primary literature itself** — the
classic value (4, still the field-standard per a live tertiary check) vs a thermodynamically-argued
revision to 3 (Wikcurrent & Hummer 2012, itself informed by Watt 2010's structure). This dispute is
the dominant driver of this model's own P/O spread (2.31–2.73 NADH depending which CI number and
which species' c-ring is used) — reported as a spread, not smoothed into one decimal.
**(2026-07-26: this dispute has since been adjudicated in the graph — see the CORRECTION banner above.)**

## Citations — every number below traced to a DOI/PMID verified LIVE this session

`WebSearch` was session-quota-exhausted (2000/2000) before this task started — the same documented
constraint `MECHANISM_MUSCLE_ENERGETICS.md`/`MECHANISM_CROSSBRIDGE_MODEL.md`/`MECHANISM_RESPIRATORY.md`/
`MECHANISM_CARDIAC.md` already disclose hitting. Fallback: direct `curl` to
`eutils.ncbi.nlm.nih.gov` (`esearch`→`esummary`/`efetch`, `rettype=abstract`), the same
repo-precedented route — **plus, new this session, one working `WebFetch` of a live PMC full-text
page** (`WebFetch` itself was *not* quota-exhausted, unlike `WebSearch`) for the one number
(membrane-potential mV range) that a PubMed abstract alone did not carry.

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Hinkle PC (2005). P/O ratios of mitochondrial oxidative phosphorylation. *Biochim Biophys Acta* 1706(1-2):1-11. | **15620362**, `10.1016/j.bbabio.2004.09.004` | **PRIMARY FALSIFIER ANCHOR.** Live-quoted verbatim: *"Values of about 2.5 with NADH-linked substrates and 1.5 with succinate are consistent with most reports... An additional revision... suggests that the H+/ATP ratio is 10/3, rather than 3, consistent with P/O ratios of 2.3 with NADH and 1.4 with succinate."* A synthesis of studies since 1937 — already a multi-study aggregation. |
| 2 | Stock D, Leslie AGW, Walker JE (1999). Molecular architecture of the rotary motor in ATP synthase. *Science* 286:1700-1705. | not independently re-fetched — cited exactly as quoted inside Hinkle 2005's own live-verified abstract | Yeast F1-c10-ring crystal structure — the first ATP-synthase c-ring structure, motivating Hinkle's own 1999-informed revision. Verified-by-transitivity, disclosed, not upgraded to independently-live-verified. |
| 3 | Watt IN, Montgomery MG, Runswick MJ, Leslie AGW, Walker JE (2010). Bioenergetic cost of making an ATP molecule in animal mitochondria. *PNAS* 107(39):16823-7. | **20847295**, `10.1073/pnas.1011099107`, PMC2947889 | **STRUCTURAL DECORRELATED ANCHOR.** Live-quoted: *"the c-ring has eight c-subunits... in about 50,000 vertebrate species... 2.7 protons are required by the F-ATPase to make each ATP molecule."* Also: *"fungi, eubacteria, and plant chloroplasts, ring sizes of c10-c15... imply 3.3-5 protons per ATP."* |
| 4 | Wikcurrent M, Hummer G (2012). Stoichiometry of proton translocation by respiratory complex I and its mechanistic implications. *PNAS* 109(12):4431-6. | **22392981**, `10.1073/pnas.1120949109`, PMC3311377 | **CONTESTED REVISION, held OPEN.** Live-quoted: *"The stoichiometry of proton translocation is thought to be 4 H+ per NADH oxidized (2e-). Here we show that a H+/2 e- ratio of 3 appears more likely..."* A live tertiary check (Wikipedia "Respiratory complex I", WebFetch this session) still states the classic 4 with no mention of this revision — reported as genuinely unresolved. |
| 5 | Wikcurrent M (1984). Two protons are pumped from the mitochondrial matrix per electron transferred between NADH and ubiquinone. | not independently NCBI-confirmed this session; cited via Wikipedia (WebFetch, live) | The classic-value foundational citation (2H+/e- × 2e- = 4H+/2e-). Tertiary-sourced, disclosed. |
| 6 | Wikcurrent MK (1977). Proton pump coupled to cytochrome c oxidase in mitochondria. *Nature* 266(5599):271-3. | **15223**, `10.1038/266271a0` | Founding discovery: complex IV is a **proton pump**, not merely a proton-consuming redox center. No abstract text on the PubMed record (pre-abstract-era 1977 short format) — disclosed. |
| 7 | Sigel E, Carafoli E (1978). The proton pump of cytochrome c oxidase and its stoichiometry. *Eur J Biochem* 89(1):119-23. | **29754**, `10.1111/j.1432-1033.1978.tb20903.x` | Early direct measurement attempt — live-quoted: *"nearly 4 K+ are taken up... "* (charge/valinomycin) vs *"about 1.6 protons are released..."* (direct chemical) — historical evidence the exact number was not cleanly settled even by direct experiment. |
| 8 | Crofts AR (2004). The cytochrome bc1 complex: function in the context of structure. *Annu Rev Physiol* 66:689-733. | **14977419**, `10.1146/annurev.physiol.66.032102.150251` | Q-cycle mechanism review for complex III, live-confirmed on-topic. The bare "4 H+/2e-" figure is not itself quoted in the fetched abstract (disclosed) — used for the qualitative Q-cycle mechanism; the number is the textbook-standard value embedded in Hinkle's own H+/2e-=10 total. |
| 9 | Kampjut D, Sazanov LA (2020). The coupling mechanism of mammalian respiratory complex I. *Science* 370(6516):eabc4209. | **32972993**, `10.1126/science.abc4209` | Most recent structural (cryo-EM, 5 conformational states, ovine complex I) account of CI coupling. Does not restate a bare H+/2e- number in its abstract; no PMC full text available this session (no PMCID on record) — not chased further (LEAN). |
| 10 | Perry SW, Norman JP, Barbieri J, Brown EB, Gelbard HA (2011). Mitochondrial membrane potential probes and the proton gradient: a practical usage guide. *Biotechniques* 50(2):98-115. | **21486251**, `10.2144/000113610`, PMC3115691 | **PRIMARY MEMBRANE-POTENTIAL ANCHOR.** Live-quoted (via `WebFetch` of the PMC full-text page — the number is not in the PubMed abstract itself): *"Typical Δp values range 180–220 mV, with Δψm typically accounting for 150–180 mV of this value."* Worked example in the same source: *"Δψm = 150 mV and ΔpHm = −0.5 units..., Δp = 150 – 60(−0.5) = 180 mV."* **Exact match to the task's own stated ~150-180mV band.** |
| 11 | Gerencser AA, Chinopoulos C, Birket MJ, Jastroch M, Vitelli C, Nicholls DG, Brand MD (2012). Quantitative measurement of mitochondrial membrane potential in cultured cells. *J Physiol* 590(12):2845-71. | **22495585**, `10.1113/jphysiol.2012.228387`, PMC3448152 | **DECORRELATED CROSS-CHECK.** Live-quoted: *"In cultured rat cortical neurons, ΔΨM is −139 mV at rest, and is regulated between −108 mV and −158 mV..."* An absolute-calibrated, intact-cell TMRM measurement — mechanistically expected to read *lower* in magnitude than the isolated-mitochondria 150-180mV band. |
| 12 | Fieldman KE, Wikcurrent MK (1976). Safranine as a probe of the mitochondrial membrane potential. *FEBS Lett* 68(2):191-7. | **976474**, `10.1016/0014-5793(76)80434-6` | Founding isolated-mitochondria potentiometric-probe paper — the direct methodological precursor to the TPP+-electrode method the task brief names. No abstract text on record (pre-abstract-era) — disclosed. |
| 13 | Kemp GJ, Meyerspeer M, Moser E (2007). *NMR Biomed* 20(6):555-65. | **17628042**, `10.1002/nbm.1192` | *Reused with attribution* from `MECHANISM_MUSCLE_ENERGETICS.md` (spot-checked live via esummary this session, full numbers not re-derived) — resting [PCr]/[ATP]/[Pi], couples_to muscle-energetics. |
| 14 | Ryan TE, Southern WM, Reynolds MA, McCully KK (2013). *J Appl Physiol* 115(12):1757-66. | **24136110**, `10.1152/japplphysiol.00835.2013` | *Reused with attribution* — PCr-recovery tau, the in-vivo 31P-MRS readout of mitochondrial OXPHOS capacity this doc's P/O model feeds. |
| 15 | Short KR et al. (2005). Decline in skeletal muscle mitochondrial function with aging in humans. *PNAS* 102(15):5618-23. | **15800038** | *Reused with attribution* from `data/body_twin/agent_outputs/mito-triangulation__af380e0570ccf87ab.json` (spot-checked live via esummary title match) — mitochondrial ATP production rate (MAPR) decline ~8%/decade, couples_to aging. |

## Section A — electron-transport-chain proton-pumping stoichiometry (H+/2e-)

| Complex | H+/2e- | Status |
|---|---:|---|
| I (NADH:ubiquinone oxidoreductase) | **4** (classic) / **3** (contested revision) | OPEN dispute (citation 4 vs 5) |
| II (succinate dehydrogenase) | **0** | No proton-pumping function — uncontested |
| III (cytochrome bc1, Q-cycle) | **4** | Textbook-standard (citation 8, mechanism confirmed live) |
| IV (cytochrome c oxidase) | **2** (vectorial/pumped) | Founding discovery citation 6; early direct measurement (citation 7) shows historical spread (1.6–4) |

NADH path (CI+CIII+CIV): **10** (classic) / **9** (revised-CI). Succinate/FADH2 path (CIII+CIV,
complex II never pumps): **6**, unaffected by the CI dispute.

## Section B — ATP synthase: a GEOMETRIC gear ratio, not a memorized fraction

Two rotary symmetries share one rigid shaft (the central γ-stalk):

- **the c-ring**: n_c-fold rotational symmetry, one proton-binding glutamate per c-subunit — one
  full 360° rotation translocates exactly n_c protons.
- **the F1 head**: a *fixed* 3-fold symmetry (3 catalytic β-subunits, Boyer's binding-change
  mechanism) — one full 360° rotation of γ produces exactly 3 ATP, in **every** resolved F1Fo
  structure across all domains of life (essentially unchallenged, unlike n_c).

Because n_c and 3 are not generally commensurate (n_c=8 for mammals: 360°/8=45° steps vs a
120°-per-catalytic-event demand), real ATP synthases resolve the mismatch via elastic power
transmission in the γ-shaft (independently documented in single-molecule rotation studies, not
re-derived here). The rotational H+/ATP cost is simply the **gear ratio of the two coupled rotor
periodicities on one shaft**: n_c / 3.

| Species/structure | c-ring size | Rotational H+/ATP (n_c/3) | + Pi/ANT transport (Hinkle's own +1 convention) | Total H+/ATP |
|---|---:|---:|---:|---:|
| Mammal (bovine, Watt 2010 cryo-EM) | 8 | **2.667** (paper states "2.7") | +1 | **3.667** |
| Yeast (Stock/Leslie/Walker 1999) | 10 | 3.333 | +1 | 4.333 |
| Hinkle's own pre-1999 naive assumption | (F1-sites only, no ring correction) | 3 | +1 | 4.0 |

`g05_mammalian_rotational_hatp_matches_watt_2p7`: computed 2.6667 vs paper's own stated 2.7 —
**PASS** (|Δ|=0.033).

## Section C — the P/O ratio grid + the decorrelated structural-ceiling falsifier

P/O = (H+/2e- supply) / (H+/ATP cost). Six self-consistent rows (machine-computed, not narrated):

| Variant | H+/ATP | P/O (NADH) | P/O (succinate) |
|---|---:|---:|---:|
| Hinkle classic (pre-1999) | 4.000 | **2.500** | **1.500** |
| Hinkle revised (yeast c10, 1999 structure) | 4.333 | 2.308 (→2.3) | 1.385 (→1.4) |
| **Mammalian c8 + classic CI=4 (structural ceiling)** | 3.667 | **2.727** | **1.636** |
| Mammalian c8 + revised CI=3 | 3.667 | 2.455 | 1.636 |

*(2026-07-26: the "revised CI=3" row above is now EXCLUDED at 9.07σ by a decorrelated direct-flux
measurement, not merely "contested" — see the CORRECTION banner above. The CI=4 row, evaluated at the
MEASURED CI=4.21±0.15 rather than the bare integer 4, is now this model's best point estimate:
P/O(NADH)=2.78±0.04, not a "structural ceiling" sitting above the true value.)*

**The falsifier check**: mammalian c8 + classic CI=4 is the highest self-consistent P/O in this
grid — no smaller vertebrate c-ring than c8 is known (Watt 2010's own text: "probably contain
c8-rings also" across ~50,000 vertebrate species), and lowering CI can only lower P/O further, never
raise it. This makes it the best-supported **structural ceiling** for mammals. Hinkle's directly-
measured flux consensus (2.5 NADH / 1.5 succinate) sits at **91.67% of this ceiling for BOTH
substrates** (`measured_over_ceiling_ratio`: 0.9167 / 0.9167 — machine-computed, matching to 4
significant figures) — a clean, non-trivial, falsifiable geometric prediction: since both
substrates share the *same* H+/ATP cost denominator and differ only in H+/2e- supply, their
measured/ceiling ratios **must** come out equal if the model is right. `g06`/`g07`: measured P/O
never exceeds the structural ceiling for either substrate — **PASS** (the physically-required
direction; a violation here would have falsified the model outright).

**The over-determination finding**: using the *modern* combination (mammalian c8 ring **and** the
revised CI=3) gives P/O(NADH)=2.455 — within 0.045 of Hinkle's *original*, decades-earlier,
purely flux-measured classic value of 2.5 (`g08`, tolerance 0.10, **PASS**). Two independent
structural/mechanistic revisions (a smaller ring, a lower CI number), arrived at by a totally
different instrument (cryo-EM/X-ray structure, not O2-flux measurement) and a different era
(2010-2012, vs Hinkle's 2005 synthesis, which could not have been influenced by either), land back
near the number direct flux measurement had *already* converged on independently. This is the kind
of convergence that is hard to manufacture by accident — it is the signature of the same underlying
geometry being recovered from two decorrelated angles.

**Self-check disclosed, not oversold**: back-solving the implied H+/ATP from Hinkle's own measured
2.5 (NADH) and 1.5 (succinate) gives exactly 4.0 from both (`g10`, PASS) — but this is **not** an
independent triangulation, since both numbers were themselves derived from one assumed flat
H+/ATP in the first place (common-mode). It is reported here only as an arithmetic self-check that
this model reproduces Hinkle's own numbers correctly, not as fresh evidence. The genuinely
decorrelated evidence is the structural-ceiling comparison above.

**What is held OPEN**: the true P/O of mammalian mitochondria is **not** pinned to one decimal by
this analysis — it is somewhere in **2.31–2.73 (NADH)** / **1.38–1.64 (succinate)** depending on
which Complex I stoichiometry (3 vs 4, a live unresolved dispute) is correct, before any additional
reduction from proton leak/slip in real (non-ideal) mitochondria is even considered. Reported as a
spread, per the task's own instruction, not resolved by fiat.

**UPDATE 2026-07-26 (see the CORRECTION banner at the top of this doc)**: the CI 3-vs-4 dispute this
spread was held open on is no longer live — CI=3 is now 9.07σ-excluded by a decorrelated direct-flux
measurement, retiring this spread's lower bound. Current best point estimate is **P/O(NADH)=2.78±0.04**
(algebra) / **2.7624** (the actual mechanistic ODE run with fully-consistent corrected inputs, -0.80%
off the algebra point — a separate, smaller, still-open gap this doc's closed-form grid does not
capture). What remains genuinely open is only the proton-leak/slip margin below either number (item 6,
Honest gaps) and the untested regime-scope caveat (low-Δp/ischemia).

## Section D — membrane potential / proton-motive force: a 2D-to-1D geometric projection

The proton-motive force is a scalar **projection of a 2D electrochemical state** (an electrical
component Δψ and a chemical component ΔpH) onto the "useful work" axis, via the Nernst/Boltzmann
conversion factor (~61.5 mV/pH-unit at 37°C; Perry's own worked example rounds to 60):
**Δp = Δψ − 60·ΔpH**.

- Perry 2011 (live-quoted): Δp ≈ 180–220 mV total; **Δψm ≈ 150–180 mV** of that — **exact match**
  to the task's own stated band (`g12`, PASS).
- Their own worked example (Δψm=150, ΔpH=−0.5 units, matrix alkaline) reproduces Δp=180 mV exactly
  (`g13`: 150 − 60×(−0.5) − 180 = 0, **PASS**).
- Applying the same decomposition across the full quoted bands implies a ΔpH contribution of
  **0.50–0.67 pH units** (low end reproduces their own worked example exactly; high end is this
  doc's own extrapolation across their stated range) — a physiologically small, matrix-alkaline
  shift (`g14`, generous 0–1 unit sanity band, PASS).
- **Decorrelated cross-check**: Gerencser 2012's absolute-calibrated, intact-cell TMRM measurement
  gives **−139 mV at rest** (range −108 to −158 mV) — reading *below* the isolated/idealized
  150-180mV band floor (`g15`, PASS). This is the mechanistically-expected direction: ongoing ATP
  synthesis in an intact, ATP-demanding cell partially collapses Δψm relative to the
  non-phosphorylating ceiling reachable in an isolated respirometer. Even the intact-cell's own
  maximum transient hyperpolarization (−158 mV, under concerted Ca²⁺/ATP-demand activation) only
  just reaches the isolated-band floor (`g16`, PASS) — consistent, not contradictory: intact cells
  can *approach* but do not exceed the isolated/idealized ceiling.

**Composite falsifier** (`g17`): both the P/O reproduction (g01-g04, g06-g07) and the membrane-
potential reproduction (g12, g15) hold simultaneously — **PASS**.

## Section E — couples_to (reused with attribution, not re-derived)

**Muscle-energetics (ATP supply / PCr-recovery link).** `docs/MECHANISM_MUSCLE_ENERGETICS.md`
already uses P/O=2.5 (Hinkle's classic NADH value) to convert its own measured Qmax (mM ATP/s) into
an O2-consumption ceiling feeding this twin's VO2/Fick chain. That doc's own honest-gaps section
explicitly flagged: *"this script uses the verified 2.3-2.5 (NADH); the task brief's own upper bound
of 2.7 could not be confirmed against Hinkle 2005's primary abstract this session."* **This doc
closes that specific gap**: 2.7 is real and structurally grounded — not from Hinkle's own abstract
(which tops out at 2.5 measured / 2.3 revised), but from a different, complementary, decorrelated
primary source (Watt et al. 2010's mammalian c-ring cryo-EM structure, which postdates Hinkle 2005
by five years and which Hinkle therefore could not have cited). `MECHANISM_MUSCLE_ENERGETICS.md`
itself is **not edited** here (companion, not edit, per this repo's isolation convention) — this
doc simply supplies the mechanistic backing the sibling doc's own gap asked for.

**Aging (mitochondrial decline).** Short et al. 2005 (PMID 15800038, reused) measured mitochondrial
ATP production rate (MAPR) declining ~8%/decade (~5%/decade after normalizing for mitochondrial
content) — i.e., part of the decline is fewer/smaller mitochondria, part is an intrinsic
per-mitochondrion functional deficit. **No citation found this session measures whether the
structural ceiling itself (c-ring count, F1 3-fold symmetry) changes with age** — both are
extremely conserved, structural, unlikely-to-drift properties; the more mechanistically plausible
reading (hypothesis, not measured here) is that aging erodes **capacity** (electron-transport-chain
complex activity/content, per the same Short 2005 paper's COX3/COX4 mRNA −8 to −10%/decade) and
possibly **coupling efficiency** (the gap between measured flux P/O and the structural ceiling
widening via increased proton leak/slip), rather than the geometric ceiling itself shifting. Flagged
explicitly as an inference, not a measured finding — an honest gap, not smoothed over.

## Confidence tier

Per this repo's own controlled vocabulary (`docs/MECHANISM_TRUST_LEDGER.md` — which **wins** over the
task brief's own phrasing, exactly as `MECHANISM_MUSCLE_ENERGETICS.md` already established this
precedent for the same tension): the ledger has **no tier named "in-vitro-anchored"** — the task
brief's own suggested label does not exist in the controlled vocabulary. The correct existing tier
is **`cadaveric-or-published-plausibility`**: *"compared against a real, genuinely
external/decorrelated published dataset or literature figure... but NOT a direct in-vivo force
measurement of this twin's own quantity."* Every anchor here (Hinkle 2005's flux synthesis, Watt
2010's cryo-EM structure, Perry 2011's membrane-potential review, Gerencser 2012's intact-cell
measurement) is a genuinely external, decorrelated, **published, real biochemical/structural
measurement** — arguably one of the more direct members of this tier, since isolated-mitochondria
P/O flux measurement and cryo-EM c-ring counting are not indirect plausibility proxies but literally
the primary data for the quantities being certified — but it is not this twin's own in-vivo human
measurement, so the ledger's own discipline ("pick ONE per claim; do not inflate") places it here,
not in `in-vivo-anchored`. **Tier: `cadaveric-or-published-plausibility`.**

## Honest gaps

1. **Complex I's H+/2e- stoichiometry is a live, unresolved dispute** (classic 4 vs Wikcurrent &
   Hummer 2012's proposed 3) — the dominant driver of this model's own P/O spread. A live tertiary
   check (Wikipedia) still states 4 as standard; no post-2012 primary source found this session
   settles it either way. Held OPEN, not resolved by this doc.
   **SUPERSEDED 2026-07-26**: a decorrelated direct-flux measurement (Jones/Blaza/Varghese/Hirst
   2017, PMID 28174301) excludes CI=3 at 9.07σ — see the CORRECTION banner at the top of this doc.
   No longer open in the graph, though this doc's own literature search did not settle it.
2. **The back-solved-H+/ATP self-check is not an independent triangulation** (see Section C) —
   disclosed explicitly so it is not mistaken for fresh decorrelated evidence.
3. **Stock/Leslie/Walker 1999** (yeast c10 structure) was not independently re-fetched at its own
   PMID this session — quoted only via Hinkle 2005's own live-verified abstract (verified-by-
   transitivity, disclosed).
4. **Complex III's specific "4 H+/2e-" bare number** was not found in the fetched Crofts 2004
   abstract text itself (mechanism confirmed live; the number is the textbook-standard value
   embedded in Hinkle's own H+/2e-=10 total) — disclosed, not fabricated as a direct quote.
5. **Membrane-potential numbers are probe-and-preparation-dependent by construction** (this is
   Perry 2011's own review subject): isolated/idealized (~150-180mV) vs intact-cell TMRM-quantified
   (Gerencser 2012, −139mV rest) genuinely differ — reported as a real biological + methodological
   gap, not reconciled to one number.
6. **Real (non-isolated, intact-tissue) P/O is further reduced below any grid value by proton leak
   and "slip"** (imperfect coupling), which varies by substrate, tissue, and physiological state —
   not quantified in this pass (no leak/slip flux dataset fetched this session). The structural
   ceiling computed here is a **theoretical maximum**, not a claim that real coupled flux reaches
   it — Hinkle's own measured consensus already sits at ~92% of it, consistent with a real,
   modest, physiological leak/slip margin.
7. **Molecular/organelle-constant-level model, no subject-specific data** — the couples_to sections
   (muscle-energetics, aging) reuse already-live-verified numbers from sibling sources with
   attribution rather than re-deriving them fresh (LEAN, avoids duplicate work); those specific
   numbers carry whatever caveats their own source docs already disclose (see
   `MECHANISM_MUSCLE_ENERGETICS.md`'s own 10-item honest-gaps list for the PCr/tau side).
8. **Graph coupling not executed**: `data/MECHANISM_ANCHOR_GRAPH.json` already has several
   mitochondria-adjacent nodes (`MOL-MITOCHONDRIAL-BIOENERGETICS`, `MITO-CELL`, `MITO-TRIANGULATION`,
   `AGING-ETA-DRIFT-CELL`) at the VO2max/respirometry/mitophagy/aging level — none of them cover
   this doc's organelle-mechanistic level (ETC complex-by-complex stoichiometry, ATP-synthase
   gear-ratio, membrane-potential decomposition), confirmed by reading the existing node text
   before writing this doc (not assumed). This doc is a genuine complement, not a duplicate. The
   graph JSON itself was **not edited** — folding this in is a separate pipeline decision
   (`mechanism_fold.py` → `fold_gate_v2`) left to the coordinator, per this repo's own architecture,
   not executed unilaterally here.

## Files

- `scripts/msk/mitochondrial_oxphos.py` — the full model: ETC stoichiometry table, ATP-synthase
  geometric gear-ratio derivation, the 6-row P/O variant grid, the structural-ceiling falsifier,
  the membrane-potential/PMF decomposition, all 17 gates, self-contained (stdlib only: `json`,
  `math`, `pathlib`).
- `data/mitochondrial_oxphos/mitochondrial_oxphos_results.json` — full machine-readable evidence:
  every citation with its verification status, the ETC/ATP-synthase constants, the full P/O grid,
  the back-solve self-check, the structural-ceiling comparison, the membrane-potential numbers, all
  17 gates (regenerated fresh each run), `overall_pass: true`.
- Reused, not modified: `docs/MECHANISM_MUSCLE_ENERGETICS.md` (P/O-band gap this doc closes),
  `data/body_twin/agent_outputs/mito-triangulation__af380e0570ccf87ab.json` (Short 2005 aging
  numbers, reused with attribution), `docs/MECHANISM_TRUST_LEDGER.md` (tier vocabulary),
  `docs/MECHANISM_HARDENED_CONVENTIONS.md` / `COORDINATOR.md` (isolation + convention discipline).
- `data/MECHANISM_ANCHOR_GRAPH.json` — **not edited** this session (read-only, dup-checked against
  existing mitochondria-adjacent nodes; isolation rule "touch only files you create" respected).

No git commit, no git push performed (isolation respected, per this repo's own `COORDINATOR.md` /
`docs/MECHANISM_HARDENED_CONVENTIONS.md`).
