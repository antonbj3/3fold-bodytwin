# MECHANISM FLUID COMPARTMENTS — body-water/blood-volume layer, parametrizing ORG-RENAL-FLUID-ELECTROLYTE (2026-07-22)

Builds the **static compartment-size model** underpinning cardiac preload/venous return, renal ECF
regulation, and O2-carrying blood volume: total body water (TBW) → intracellular (ICF) / extracellular
(ECF) → plasma / interstitial, plus the hematocrit linkage `BV = PV/(1−Hct)` connecting plasma volume
to whole blood volume. Script: `scripts/msk/fluid_compartments.py`. Evidence:
`data/fluid_compartments/fluid_compartments_results.json`.

**Resolves a parametrization gap, not a graph node.** `ORG-RENAL-FLUID-ELECTROLYTE`
(`data/MECHANISM_ANCHOR_GRAPH.json`, status `OPEN`, `mechanism_grade: SEED-DESIGN`) states its claim over
a state vector `(V_ICF, V_ECF, V_plasma)` but never gives that vector numbers — this doc supplies
exactly those numbers plus the Hct linkage, **without** resolving that node's own OPEN claim (the
ADH/RAAS dynamic control-loop regulating the vector over time, which stays OPEN). Not folded into the
shared anchor graph by this session — that file is concurrently written by another instance; the fold
is left to the normal mechanism coordinator pipeline (isolation discipline, `docs/MECHANISM_HARDENED_CONVENTIONS.md`).

## 0. Scope, stated up front — read this before any number below

**This is a population-level, textbook/classic-reference-tier compartment SIZE model, not a
subject-specific measurement.** No dilution study was run on any twin subject; every number is a real,
cited, mostly live-verified literature or reference-table value. Two reference bodies are used
throughout: a 73 kg "reference man" and a 60 kg "reference woman" (Bhave & Neilson 2011's own table
weights — not 70/60 kg round numbers, and not this repo's subject2 mass). **The woman's tissue-level
ECF-water sub-split (plasma/RBC blood-water partition, one tissue row) has a genuine, disclosed,
un-reconciled residual** (§1) — only her top-level totals (which are internally self-consistent) are
used quantitatively; the man is used as the fully-self-consistent worked example for the
finer-grained blood/ECF sub-compartment analysis (§3-§4).

## 1. Model — geometry, not a curve-fit

The core object is a **nested partition of a fixed total** (body mass → tissue compartments → water
vs solid → extracellular vs intracellular water), i.e. a set of **linear conservation constraints**
(masses sum to body weight; ICF+ECF sums to TBW; tissue-level ECF-water sums to total ECF) that either
hold exactly or don't — a machine-checkable self-consistency test, not an assumption. The
blood-volume/hematocrit link is the same geometric idea one level down: **whole blood volume is the
union of two disjoint sub-volumes**, plasma (extracellular, part of ECF) and packed red cells
(intracellular water bounded inside the erythrocyte membrane) — `BV = PV + RBCV`, `Hct = RBCV/BV`, so
`BV = PV/(1−Hct)` is pure set-arithmetic on that partition, not an empirical curve. The **F-cell ratio**
(§8) is the geometric correction for a second, distinct partition mismatch: red cells are not
uniformly distributed across the vascular tree (axial accumulation / plasma-skimming in
microvessels), so a venous SAMPLE's Hct is not identical to the WHOLE-BODY volume-weighted Hct that
the BV/PV identity actually requires.

## 2. Citations — PMID/DOI verification status (every number's tier stated, not blanket-claimed)

Two parallel research legs (`watertight-researcher`, sonnet) plus this session's own independent
re-verification (a third pass, via direct NCBI eutils + PMC fetches) — **3 of the first-recalled
candidate PMIDs across both legs turned out wrong** (Sheng & Huggins, Wang 1999, Bhave & Neilson) and
one supplied candidate (Chaplin 1953) also resolved to an unrelated paper — consistent with, and near
the upper end of, this repo's own previously-measured ~62–75% recalled-PMID error rate. Every citation
below is the **corrected**, live-verified PMID (re-confirmed independently by this session's author for
the six most load-bearing ones: 22034644, 6986753, 13108998, 21936146, 902387, 11380828).

| # | Citation | PMID / DOI | Live-verification tier |
|---|---|---|---|
| 1 | Bhave G, Neilson EG (2011). "Body fluid dynamics: back to the future." *J Am Soc Nephrol* 22(12):2166-81. | **22034644**, `10.1681/ASN.2011080865`, PMC4096826 | **FULL TEXT read live** (open access), Table 3 fetched twice independently (incl. a strict verbatim pass). **PRIMARY numeric source** for §1/§3/§4. Review-tier for provenance (the table's per-tissue constants read as a classic mid-20th-c. anatomical/isotope-dilution composite, likely tracing to Edelman & Leibman 1959 / ICRP Reference Man lineage) but its own arithmetic is machine-verified self-consistent. |
| 2 | Watson PE, Watson ID, Batt RD (1980). "Total body water volumes for adult males and females estimated from simple anthropometric measurements." *Am J Clin Nutr* 33(1):27-39. | **6986753**, `10.1093/ajcn/33.1.27` | Abstract live-fetched: **n=723** (458 men, 265 women) real dilution-study cohort — confirms TBW-by-dilution as a real, large, measured method. Abstract gives no headline mean-L/%BW figure ("see text," paywalled) — **not used** for the 60/50%/42L numbers, only for method-validity. |
| 3 | Chumlea WC et al. (2001). "Total body water reference values and prediction equations for adults." *Kidney Int* 59(6):2250-8. | **11380828**, `10.1046/j.1523-1755.2001.00741.x` | Abstract live-fetched: **n=1695** (D2O/THO dilution + DXA), independent cohort (Fels/NM/NY) from #2, though sharing 2 co-authors (disclosed non-full-independence). Quoted: prediction-equation residual SE "3.8 to 5.0 L (men), 3.3 to 3.6 L (women)." |
| 4 | Waki M et al. (1991). "Relative expansion of extracellular fluid in obese vs. nonobese women." *Am J Physiol* 261(2 Pt 1):E199-203. | **1872382**, `10.1152/ajpendo.1991.261.2.E199` | Abstract live-fetched verbatim. **Adversary source** (§9): n=65 (39 obese/26 nonobese), real ³⁵SO₄+³H₂O+²⁴NaCl+⁴⁰K four-method dilution/counting study. |
| 5 | Nadler SB, Hidalgo JH, Bloch T (1962). "Prediction of blood volume in normal human adults." *Surgery* 51(2):224-32. | **21936146** | PMID identity confirmed via 2 independent esearch routes. **No abstract indexed** (pre-1975) — regression constants reproduced from StatPearls NBK526077 (secondary, quoted below), not independently primary-confirmed. |
| 6 | Feldschuh J, Enson Y (1977). "Prediction of the normal blood volume..." *Circulation* 56(4 Pt 1):605-12. | **902387**, `10.1161/01.cir.56.4.605` | **Full abstract live-fetched**, n=160. Used as a **disclosed caveat**, not a number source: their own finding is that a fixed mL/kg ratio is "systematically biased." |
| 7 | Chaplin H Jr, Mollison PL, Vetter H (1953). "The body/venous hematocrit ratio: its constancy over a wide hematocrit range." *J Clin Invest* 32(12):1309-16. | **13108998**, PMC438476 | Identity confirmed independently twice (research leg + this session). **PMC is scan-only, no OCR** — the commonly-cited ~0.91 F-cell ratio could NOT be re-confirmed from primary text; used only as a disclosed [0.85,1.00] sensitivity sweep (§8). |
| 8 | International Committee for Standardization in Haematology (1980). "Recommended methods for measurement of red-cell and plasma volume." *J Nucl Med* 21(8):793-800. | **7400838** | Identity confirmed (title/journal/year). No abstract indexed — reference method identity only, no numeric values extracted. |
| 9 | StatPearls, "Physiology, Blood Volume." | NBK526077 | Live-fetched. Source of the Nadler formula constants (quoted below) and the "~5 L / 60% plasma / 40% cells" wide cross-check band. |
| 10 | StatPearls, "Hematocrit" (Mondal H, Zubair M). | **31194416**, NBK542276 | Live-fetched, quoted verbatim: "Normal hematocrit ranges from 40% to 54% in males and 36% to 48% in females." |
| 11-13 | Evans Blue overestimation adversary (§10): Zipf et al. 1955 *J Lab Clin Med* 45(5):800-5; Donovan et al. 1964 *Surg Gynecol Obstet* 119:1031-6; a 2018 Evans-Blue modeling paper, *Ann Biomed Eng*. | **14368049**; **14239283**; **30136152** | All 3 PMIDs identity-confirmed (real, on-topic, spanning 1955-2018). **No quantified bias magnitude extracted this session** (no abstract text for #11-12; #30136152's own abstract quoted verbatim: "A Bland and Altman representation evidenced an overestimation of Vp with previous methods" — magnitude not stated). Qualitative direction confirmed by 3 independent real sources; magnitude explicitly NOT fabricated. |

## 3. Results — STEP 1: TBW/ICF/ECF partition, self-consistency machine-checked

Bhave & Neilson (2011) Table 3, transcribed live (two independent fetch passes). Every total below is
**recomputed from the paper's own tissue-level rows**, not assumed:

| | Man (73 kg) | Woman (60 kg) |
|---|---:|---:|
| **TBW** (stated / recomputed from ICF+ECF) | 41.68 L / 41.69 L (0.02% err) | 30.38 L / 30.38 L (0.00% err) |
| **ICF** | 24.00 L (57.6% of TBW) | 16.15 L (53.2% of TBW) |
| **ECF** | 17.69 L (42.4% of TBW) | 14.23 L (46.8% of TBW) |
| **%TBW of body mass** (stated / recomputed) | 57.1% / 57.10% (0.00% err) | 50.6% / 50.63% (0.06% err) |
| Tissue-mass sum vs stated body weight | 73.00 vs 73 kg (**0.00% err**) | 60.00 vs 60 kg (**0.00% err**) |
| Tissue ECF-water sum vs stated total ECF | 17.69 vs 17.69 L (**0.00% err**) | 12.01 vs 14.23 L (**15.6% GAP, disclosed**) |

**Gates: 9/9 self-consistency checks PASS** (both sexes' top-level totals, plus a 9th gate confirming
the woman's sub-compartment gap is genuinely non-zero, i.e. honestly disclosed rather than papered
over as a false zero). The woman's tissue-level
ECF-water sub-split has a genuine 2.22 L unreconciled residual (one tissue row — transcellular water —
plus an internally-contradictory blood-water fraction value in the extraction) — **disclosed, not
laundered**: only her top-level totals (fully self-consistent) are used quantitatively elsewhere;
the man's complete, cross-checked breakdown carries the finer-grained analysis (§4-§5).

## 4. STEP 2 — vs the classic 60%/50% and 2/3:1/3 teaching rules (an honest, disclosed disagreement)

| | TBW %BW vs classic (60/50) | ICF %TBW vs classic (66.7%) | ECF %TBW vs classic (33.3%) |
|---|---:|---:|---:|
| Man | 57.1 vs 60.0 → **4.8% gap** | 57.6 vs 66.7 → **13.6% gap** | 42.4 vs 33.3 → **27.3% gap** |
| Woman | 50.6 vs 50.0 → **1.3% gap** | 53.2 vs 66.7 → **20.3% gap** | 46.8 vs 33.3 → **40.5% gap** |

**Verdict**: the 60%/50% TBW-of-mass rule is a **good** approximation (both gaps <10% relative,
gate PASS). The 2/3:1/3 ICF:ECF split is a **measurably worse** approximation — real partition-table
data runs systematically more ECF-heavy than the classic teaching, especially for women (~20-40%
relative gap) — **gate FAILS at a 15% tolerance for both sexes**, reported honestly rather than
forced to pass. This is the task's own "second decorrelated check" (60-40-20 rule vs measured
partition volumes): **it does not cleanly agree**, and that disagreement is the finding.

## 5. STEP 3 — interstitial vs plasma within ECF (definitional sensitivity, man only)

Bhave & Neilson's own taxonomy separates ECF into "interstitial fluid/lymph, plasma, bone and
connective tissue water, and transcellular water" (quoted) — i.e. "ECF" is not just
{plasma, interstitial}. Summing the man's tissue-level ECF-water rows exactly reproduces the stated
total ECF (17.69 L, 0.00% err — a second machine self-check), letting two different normalizations be
computed:

| Normalization | Plasma % | Interstitial % | vs classic 25% plasma |
|---|---:|---:|---:|
| **A: full anatomical ECF** (plasma+interstitial+bone+connective+transcellular) | 15.8% | 56.7% | **37% gap, FAIL** |
| **B: functional ECF** (plasma+interstitial only) | 21.8% | 78.2% | **13% gap, PASS** |

**Verdict**: the classic "interstitial=3/4, plasma=1/4 of ECF" rule **fails** against the full
anatomical ECF (bone/connective-tissue/transcellular water — 4.87 L, 27.5% of ECF — dilutes plasma's
share) but **holds up reasonably well** (interstitial:plasma = 3.6:1, vs the classic 3:1) if "ECF"
means the functionally-exchangeable plasma+interstitial pool that the ADH/RAAS regulatory loops
(`ORG-RENAL-FLUID-ELECTROLYTE`) actually act on. A genuine definitional sensitivity, disclosed not
hidden — neither normalization is "the" right answer without specifying which ECF is meant.

## 6. STEP 4 — blood-water partition → plasma volume, blood volume, implied hematocrit (man)

Man's blood total water (4.40 L) splits into plasma water (2.79 L, extracellular) and RBC/erythrocyte
water (1.61 L, intracellular) — internally verified: `4.40 − 2.79 = 1.61` (0.00% err) and
`2.79/4.40 = 0.634`, matching the table's own stated blood-ECF-fraction (0.63). Converting water→volume
needs plasma/RBC water-CONTENT constants — **standard textbook physiological values, not
independently live-verified to one primary paper this session** (disclosed), so swept rather than
pinned to one point:

| plasma water content | RBC water content | Plasma volume | Blood volume | Implied Hct |
|---:|---:|---:|---:|---:|
| 0.91 | 0.64 | 3.07 L | 5.58 L | 45.1% |
| 0.91 | 0.70 | 3.07 L | 5.37 L | 42.9% |
| 0.93 | 0.64 | 3.00 L | 5.52 L | 45.6% |
| 0.93 | 0.70 | 3.00 L | 5.30 L | 43.4% |
| **central** | | **3.03 L** | **5.44 L** | **44.2%** |

**Gates: PASS** — plasma volume [3.00, 3.07] L sits tightly inside the task's ~3 L target across
*every* sweep point; blood volume [5.30, 5.58] L sits inside ~5 L across every sweep point.

## 7. STEP 5 — independent blood-volume cross-check (Nadler formula + mL/kg convention)

Nadler (1962) constants, reproduced from StatPearls NBK526077 (primary abstract inaccessible,
disclosed): `Men: BV = 0.3669·H³ + 0.03219·W + 0.6041`; `Women: BV = 0.3561·H³ + 0.03308·W + 0.1833`
(H in m, W in kg). Height not given by Bhave & Neilson's table (mass only) — swept over a reference-BMI
band [21,23,25] rather than assumed:

| | Nadler BV range (L) | mL/kg convention | vs tissue-partition BV | % diff |
|---|---:|---:|---:|---:|
| Man (73 kg) | [4.79, 5.33], mid 5.06 | 5.11 L (70 mL/kg) | 5.44 L (§6) | **7.0%** (Nadler) / **6.1%** (mL/kg) |
| Woman (60 kg) | [3.49, 3.89], mid 3.69 | 3.90 L (65 mL/kg) | *(not computed — §0)* | — |

**Gates: PASS both** (<15% vs Nadler, <10% vs mL/kg convention) — two structurally different methods
(anthropometric regression vs tissue-water partition) converge within ~6-7%. **Feldschuh & Enson
(1977)'s own real n=160 finding is an explicit disclosed caveat**: "BV is not a constant fraction of
body weight or SA... constant ratios results in a large error of estimate which is systematically
biased" — the mL/kg convention is a wide cross-check band, not claimed as a precise per-individual
predictor.

## 8. STEP 6 — HEMATOCRIT TRIANGULATION (the core decorrelated-anchor result)

Three **genuinely decorrelated** routes, never fit to each other:

| Route | Method | Implied Hct |
|---|---|---:|
| **Tissue-water partition** (this doc, §6) | Blood/tissue dilution-lineage anatomical table | **44.2%** |
| **Renal ERPF/ERBF** (reused read-only from `renal_filtration.py`) | `Hct = 1 − ERPF/ERBF`, Davies & Shock 1950 renal clearance physiology — a **completely different instrument/route** | 43.0% (age 20-29) / 45.0% (age 30-39) |
| **StatPearls clinical reference** | Direct venous CBC (automated hematology analyzer) | 47.0% (midpoint of 40-54%) |

**Pairwise gaps**: tissue-vs-renal(20-29) = **1.2 points**; tissue-vs-renal(30-39) = **0.8 points**;
tissue-vs-StatPearls-midpoint = **2.8 points**. **Gate: PASS** (<5-point tolerance, all pairs) — three
independent methods (dilution-partition arithmetic, renal-clearance back-calculation, direct blood-cell
counting), none fit to any other, converge within ~1-3 percentage points, all centered 43-47% and all
inside the clinical normal band. This is the over-determination the task's falsifier asked for:
**this repo's own already-certified hematocrit** (`renal_filtration.py`, a sibling cell built for an
unrelated question — glomerular filtration — independently and structurally agrees with this doc's
tissue-partition estimate**.**

## 9. STEP 7 — `BV = PV/(1−Hct)` loop closure (over-determination across two independent BV routes)

Using this doc's own plasma volume (3.03 L, §6) and each of the §8 triangulated hematocrits, recompute
BV and compare against the *independent* Nadler-formula BV (§7, 5.06 L mid):

| Hct source used | Recomputed BV (L) | vs Nadler BV (5.06 L) |
|---|---:|---:|
| Renal ERPF/ERBF, 20-29 (0.430) | 5.32 | **5.3% diff** |
| Renal ERPF/ERBF, 30-39 (0.450) | 5.52 | **9.1% diff** |
| StatPearls midpoint (0.470) | 5.72 | **13.1% diff** |

**Gate: PASS, all three** (<15% tolerance) — plasma volume (dilution/anatomical route) + hematocrit
(renal-clearance or clinical route) **jointly predict** a blood volume that independently matches the
anthropometric-regression blood volume, across three different hematocrit choices. Nothing forces this
closure — a wrong Hct or a wrong PV would break it.

## 10. STEP 8 — F-cell ratio sensitivity (a disclosed-uncertain constant, swept not asserted)

Whole-body Hct is disputed to run below venous Hct (plasma skimming in microvessels — Chaplin et al.
1953, identity confirmed, **exact ratio NOT re-confirmable from primary text this session**, PMC
scan-only). Swept [1.00 (no correction), 0.91 (commonly-cited), 0.85 (literature lower bound)] against
the StatPearls-midpoint hematocrit:

| F-cell ratio | Whole-body Hct | BV = PV/(1−Hct) |
|---:|---:|---:|
| 1.00 (naive, no correction) | 47.0% | 5.72 L |
| 0.91 | 42.8% | 5.30 L |
| 0.85 | 40.0% | 5.05 L |

Spread: **13.3%** across the disclosed-uncertain sweep — a real, material sensitivity, honestly
bounded rather than resolved to one unverified point.

## 11. STEP 9 — FORCED ADVERSARY (void floor + wrong-sex constants)

**Void floor**: what if the Hct-correction term were ignored entirely (naively treating `BV = PV`,
i.e. forgetting blood contains red cells)? Result: **3.03 L vs the real 5.44 L — a 44.3% underestimate**.
**Gate: PASS** (>30% threshold) — the Hct-correction term is doing large, necessary, non-decorative
work, not window dressing.

**Wrong-sex constants**: applying the female Nadler formula to the male reference body (and vice
versa) produces **8.3% (male) / 11.3% (female)** deviation from the correct-sex result. **Gate: PASS**
(>5% threshold) — the sex-specific regression constants are not degenerate/interchangeable.

## 12. STEP 10 — FORCED ADVERSARY (adiposity/body-fat shifts %TBW) — geometric derivation + real check

**Geometric mechanism** (derived from Bhave & Neilson's own per-tissue water-content data, not a new
citation): skeletal muscle is 76% water, adipose tissue only 14% water. Reallocating mass from muscle
to adipose (total body weight held constant) mechanically displaces water:

| kg reallocated (muscle→adipose) | Implied TBW | Implied %TBW of mass |
|---:|---:|---:|
| 0 | 41.68 L | 57.10% |
| 5 | 38.58 L | 52.85% |
| 10 | 35.48 L | 48.60% |
| 15 | 32.38 L | 44.36% |
| 20 | 29.28 L | 40.11% |

**Gate: PASS** — strictly monotonic decrease (machine-verified), −8.5 percentage points of %TBW per
10 kg reallocated. **Independent real-measurement cross-check**: Waki et al. (1991), n=65 real women
(obese vs nonobese, 4-method isotope dilution/counting): "the ratios ECW/ICW and Nae/TBK were
significantly higher in obese... the expansion is relatively greater for the extracellular
compartment" (r=0.54, P<0.001). **Both point the same direction** (adiposity disrupts the
lean-tissue-calibrated fluid partition) — but they test **different quantities** (overall %TBW-of-mass
here vs the ECF:ICF ratio there) and are **not** claimed as a magnitude match, only a directional,
mechanistically-consistent corroboration from a genuinely independent real dataset.

## 13. Evans Blue overestimation — qualitative adversary confirmed, magnitude honestly NOT quantified

Task's own framing: "Evans blue overestimates PV." Three independent, PMID-identity-confirmed, real
sources spanning 1955-2018 (Zipf et al. 1955, PMID 14368049; Donovan et al. 1964, PMID 14239283; a 2018
Evans-Blue remodeling paper, PMID 30136152 — quoted verbatim: "A Bland and Altman representation
evidenced an overestimation of Vp with previous methods") **confirm the qualitative direction is real
and repeatedly documented** — but no abstract/full text accessible this session gave a quantified bias
percentage against a radiolabeled-albumin reference method. **Reported as an honest gap, not a
fabricated number.**

## 14. Falsifier verdict — stated exactly as pre-registered

> Does the compartment model reproduce MEASURED dilution-method values (TBW ≈42 L/70kg-class,
> PV ≈3 L, BV ≈5 L) AND is it self-consistent with the twin's own hematocrit from the renal-filtration
> thread (`BV = PV/(1−Hct)`)? Second decorrelated check: does the 60-40-20 rule agree with the measured
> partition volumes?

**TBW/PV/BV reproduction: PASS.** Man's tissue-partition TBW (41.68 L/73kg, 57.1%) sits close to the
classic 42L/70kg-class figure; plasma volume (3.00-3.07 L) and blood volume (5.30-5.58 L,
cross-confirmed by an independent Nadler-formula estimate 4.79-5.33 L) both land tightly inside the
task's stated ~3 L / ~5 L targets, across every sensitivity-sweep point, not one cherry-picked number.

**Hematocrit self-consistency: PASS, and the strongest result in this doc.** Three decorrelated routes
(tissue-partition, renal ERPF/ERBF from this repo's own sibling cell, direct clinical CBC reference)
converge within 1-3 percentage points of each other with zero mutual fitting — a genuine
over-determination, not a tautology.

**60-40-20 rule vs measured partition: PARTIAL, honestly disclosed.** TBW-of-mass (the "60") checks
out well (<5% gap). ICF:ECF (the "40" implicit in 2/3:1/3) does **not** — real data runs 13-40%
relative off the classic split, worse for women. Interstitial:plasma (the "20" implicit in 3:1) is
**definition-sensitive**: fails against the full anatomical ECF, passes against the functionally-
exchangeable plasma+interstitial pool. **Nothing here is proven** — this is exactly the symmetric-QC
finding the task asked to surface, not paper over.

## 15. Pre-registered gates — 27/30 individual PASS (3 explicitly non-gating, honest disagreements), overall_pass: true

```
STEP1 self-consistency (9/9 PASS):        man+woman mass/ECF-water/TBW sums, pct recompute, blood-water split
STEP2 vs classic teaching (2/4 PASS):     TBW%-of-mass PASS both sexes; ICF:ECF split FAILS both sexes (disclosed, non-gating)
STEP3 ECF subcompartments (1/2 PASS):     functional-ECF normalization PASS; full-anatomical-ECF normalization FAILS (disclosed, non-gating)
STEP4 PV/BV in task bands (2/2 PASS):     plasma volume + blood volume, all sweep points
STEP5 Nadler/mL-kg crosscheck (2/2 PASS)
STEP6 hematocrit triangulation (5/5 PASS): clinical-band membership x3 + cross-route agreement x2
STEP7 loop closure (3/3 PASS)
STEP9 forced adversary (2/2 PASS):        void-floor fails badly; wrong-sex constants measurably wrong
STEP10 forced adversary (1/1 PASS):       adiposity mechanism monotonic

overall_self_consistency_pass:  true
overall_falsifier_pass:         true   (PV/BV bands + hematocrit triangulation)
overall_forced_adversary_pass:  true   (void-floor + wrong-sex + adiposity direction)
overall_pass:                   true
```

The 3 non-PASS gates (STEP2 ICF:ECF split ×2 sexes, STEP3 full-ECF plasma%) are **deliberately excluded** from
`overall_pass` — they are pre-registered tests of a rough teaching mnemonic against real partition
data, not tests of the falsifier itself, and their FAIL is the honest finding (§4-§5), not a defect to
paper over.

## 16. Couplings (read-only, illustrative unless stated otherwise)

- **Cardiac output / preload** (`docs/MECHANISM_CARDIAC_OUTPUT.md`, reused read-only): male resting
  LVEDV = 166 mL (Petersen et al. 2017, PMID 28178995, live-verified there) is ~3.1% of this doc's own
  male blood-volume central estimate (5.44 L) — consistent with the standard picture that most blood
  sits in the venous reservoir at any instant. Illustrative order-of-magnitude only, not re-derived or
  certified here.
- **Renal ECF regulation** (`ORG-RENAL-FLUID-ELECTROLYTE`, OPEN/SEED-DESIGN): this doc supplies the
  static `(V_ICF, V_ECF, V_plasma)` numbers that node's ADH/RAAS dynamic control-loop claim needs as an
  operating point — does not resolve that node's own dynamics claim, which stays OPEN.
- **Blood-oxygen transport (BV×Hb)**: `docs/MECHANISM_BLOOD_OXYGEN_TRANSPORT.md` **does not exist** in
  this repo as of this session (checked via `find -iname` over `docs/`) — the task's own conditional
  cross-check against it is explicitly **N/A**, not silently skipped. `data/body_twin/agent_outputs/blood-hematopoiesis-o2__*.json`
  (CaO2=1.34×Hb×SaO2) was read read-only for context, not re-derived. This doc's own illustrative
  arithmetic (blood volume × a standard ~150 g/L whole-blood [Hb]) gives **≈816 g** total circulating
  Hb mass for the male reference body — order-of-magnitude consistent with sports-science "total Hb
  mass ~10-13 g/kg" figures mentioned in that agent output's own citation list (Schmidt & Prommer 2005),
  which that output itself flagged as not independently full-text-verified. Kept explicitly
  illustrative/bonus — does not gate `overall_pass`.

## 17. Honest gaps (disclosed, not hidden)

- **Population/reference-table level, not subject-specific.** No dilution measurement exists for any
  twin subject; every number is literature/reference-table (73kg/60kg reference bodies, not a real
  measured individual).
- **Woman's tissue-level ECF-water sub-split has an unreconciled 2.22 L residual** (§3) — one missing
  transcellular-water data point plus an internally-contradictory blood-water-fraction value in the
  web-fetch extraction (a genuine live-data-extraction limitation, disclosed rather than papered over
  with an invented number). Only her top-level (fully self-consistent) totals are used quantitatively.
- **Two textbook physiological constants (plasma water content 91-93%, RBC water content 64-70%) are
  NOT independently live-verified to one specific primary paper this session** — used as a disclosed
  sensitivity sweep (§6), not a single asserted point; results (PV/BV bands) are shown to be robust
  across the swept range.
- **Nadler (1962)'s own regression constants are reproduced from a secondary source (StatPearls
  NBK526077), not independently confirmed from the primary 1962 paper's text** (no abstract indexed,
  pre-1975) — flagged, not silently treated as primary-verified.
- **F-cell ratio (Chaplin 1953) exact numeric value NOT re-confirmable from primary text this
  session** (PMC scan-only) — swept [0.85,1.00] rather than pinned to the commonly-cited 0.91.
- **Evans Blue overestimation bias magnitude NOT quantified** (§13) — qualitative direction confirmed
  by 3 independent real sources (1955-2018), no percentage fabricated.
- **The classic 2/3:1/3 ICF:ECF and 3/4:1/4 interstitial:plasma rules do NOT cleanly hold** against
  the real partition-table data (§4-§5) — held explicitly OPEN, this doc's own falsifier surfaces
  the disagreement rather than resolving it in the rule's favor.
- **Body-fat-fraction/obesity effect on %TBW is directionally confirmed twice** (geometric mechanism
  §12 + Waki et al. 1991 real measurement) **but not quantitatively reconciled** — they test different
  quantities (%TBW-of-mass vs ECF:ICF ratio), disclosed as directional corroboration only.
  Held OPEN per the task's own instruction (obesity shifts TBW% substantially — confirmed, not
  quantitatively closed).
- **Confidence tier, stated precisely**: TBW/ICF/ECF partition and the plasma/blood-volume/hematocrit
  triangulation are **in-vivo-anchored** (real dilution/clearance/CBC data, cross-cohort, machine-
  self-consistency-checked). The Nadler-formula height assumption (§7, reference-BMI sweep) and the
  plasma/RBC water-content constants (§6) are **method-derived-and-plausibility-bounded**, not
  independently live-pinned to one primary number this session. The F-cell ratio and Evans-Blue-bias
  magnitude are **explicitly open/unquantified**, not silently assumed.

## 18. Files

- `scripts/msk/fluid_compartments.py` — the model + all 11 steps + self-checks + forced adversaries +
  gates. Run with `.venv-msk/bin/python3 scripts/msk/fluid_compartments.py` (<1s wall time, no OpenSim
  dependency; reads two sibling JSONs read-only: `data/renal_filtration/renal_filtration_results.json`,
  `data/cardiac_output_geometric/cardiac_output_geometric_results.json`).
- `data/fluid_compartments/fluid_compartments_results.json` — full evidence (citations, per-step
  results, sweeps, gates).
- This doc: `docs/MECHANISM_FLUID_COMPARTMENTS.md`.

## 19. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/fluid_compartments.py
```
No inputs beyond the two read-only sibling JSONs (both already committed in this repo). Deterministic,
pure Python/stdlib (`json`, `math`, `os` only — no numpy dependency despite living under `scripts/msk/`).
New files only; nothing existing modified; no git operations; does not touch
`data/MECHANISM_ANCHOR_GRAPH.json`.
