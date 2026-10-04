# MECHANISM SPERMATOGENESIS — the Sertoli-supported meiotic production line & blood-testis barrier
# (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** This is the PRODUCTION LINE, not the product and
not the systemic control loop — two sibling docs already cover those and are explicitly NOT
duplicated here:

- `docs/MECHANISM_SPERM_MOTILITY_AXONEME.md` = the flagellum (9+2 axoneme, dynein sliding -> bend
  wave -> swimming). That doc's own honest gap #6 names this exact production-line axis as
  out of scope for it.
- `docs/MECHANISM_REPRODUCTIVE_HPG.md` = the systemic GnRH -> LH/FSH -> testosterone endocrine
  feedback loop (the hormones as free-standing signals, not what they do *inside* the tubule).
- `REPRO-SPERMATOGENESIS` node (`data/MECHANISM_ANCHOR_GRAPH.json`, read-only, status OPEN) —
  despite the name, its actual content is the DNA-fragmentation-index / time-to-pregnancy /
  ART-outcome fertility construct (semen quality as a functional OUTPUT), not the cellular
  production-line mechanism. Confirmed by direct inspection this session (not assumed).
- `REPRO-AZOOSPERMIA-RETRIEVAL-GENOTYPE-GATE` — AZF/Klinefelter genotype -> micro-TESE retrieval
  prognosis. A genotype/retrieval-outcome axis, not the Sertoli-support mechanism.

**This doc's actual delta**: spermatogonial mitosis (self-renewal + differentiation) -> primary
spermatocyte (meiosis I) -> secondary spermatocyte -> spermatid (meiosis II) -> spermiogenesis,
run *inside* Sertoli-cell-nursed, blood-testis-barrier-sequestered seminiferous tubules, timed by
the measured ~16-day epithelial cycle / ~64-day total transit, and gated by intratesticular (not
serum) testosterone. Confirmed FRESH by direct graph inspection (grep'd the full 1049-node
`data/MECHANISM_ANCHOR_GRAPH.json` for `SPERMATOGEN|SERTOLI|SEMINIFEROUS` — 7 hits, none of them
this content; grep'd `bt_memory/` — zero hits) before building, per this repo's own "don't
re-litigate a covered topic" discipline.

**Scope honesty, stated up front**: this is a LITERATURE-SYNTHESIS cell, not a numeric
forward-simulation like its MSK siblings (there is no differential equation to integrate here —
the object is a cellular production line + a barrier, not a continuous dynamical system). The
discipline this doc upholds instead: every quantitative claim below is backed by a
**machine-computed** arithmetic/statistical gate (ratio comparisons, a Fisher exact test computed
from scratch, proportionality checks) run over numbers extracted **verbatim** from live-fetched
primary-source abstracts — never a narrated impression, never a recalled-from-training-data
number presented as verified. Script: `scripts/msk/spermatogenesis_sertoli.py`. Evidence:
`data/msk_smoketest/spermatogenesis_sertoli/spermatogenesis_sertoli_results.json` (byte-identical
across 2 independent process runs, timestamp field excepted — verified via `diff` after stripping
the self-reported generation time; NaN/Inf-free, checked programmatically over the full JSON tree
after one genuine catch, see §7) — also written verbatim to
`docs/MECHANISM_SPERMATOGENESIS_evidence.json` per this task's own requested path.

All 13 PMID/DOI citations below independently live-verified this session via NCBI E-utilities
(esearch then efetch, direct curl, full abstract text fetched and read — not WebFetch-summarized,
not recalled from memory; this repo's own prior finding across sibling docs is a measured
~62–67% citation-drift rate from memory alone).

---

## 0. Pre-registration — falsifiers, thresholds, and the adversary each one forces

| # | Claim | Pre-registered threshold | Adversary forced |
|---|---|---|---|
| F1 | Total spermatogenesis transit ≈ 64 days, organized as ~4 iterations of a 16-day epithelial cycle | 2 methodologically-independent eras/techniques must agree within 20% relative difference | "The 64-day figure is a stale, single-study (1963), small/ethically-fraught-sample artifact, never independently reproduced" |
| F2 | Sertoli cell number sets a hard ceiling on adult germ-cell output — germ cells do NOT autonomously compensate | Experimentally reduced Sertoli population -> proportional (within 20% relative deviation) reduction in germ-cell output, with the general-toxicity confound explicitly ruled out | "Germ cells autonomously divide without Sertoli support" (the task's own named adversary) |
| F3 | Intratesticular testosterone (ITT) is ≥1 order of magnitude above serum in every measured human cohort, and IS the operative local threshold variable (not serum T) for quantitatively normal output | ITT/serum ratio ≥10x in every independently measured cohort; suppressing ITT to serum-comparable levels (even with normal/high serum T) collapses sperm output by ≥50% | "The '100x' figure is a single-study/single-method artifact" AND "serum T, not ITT specifically, is what actually matters" |
| F4 | Blood-testis-barrier tight junctions confer immune privilege on post-meiotic germ cells; barrier-relevant structural lesions associate with measurable antisperm-antibody (ASA) formation, near-absent in matched normals | Control-group ASA rate ≈0%; lesion-group rate significantly >0 (Fisher exact p<0.05, computed fresh, not just quoted); ≥90% of ASA+ cases must localize to an actual local structural/genital lesion (ruling out a generic disease-severity bystander effect) | "ASA formation is a generic autoimmune/inflammatory bystander phenomenon of systemic illness, unrelated to barrier breach specifically" |
| F5 | Sertoli-cell-only syndrome (SCOS/Del Castillo) demonstrates the dependency is ASYMMETRIC: Sertoli cells + organized tubules can exist with zero germ cells; the reverse is not reported | Direct histopathological description: Sertoli cells present, tubules organized, germ cells/spermatogonia absent | The implicit symmetric-mutual-dependency null (germ cells and Sertoli cells are equally load-bearing for each other's presence) |

---

## 1. Citations (13, every PMID/DOI independently live-verified this session)

| # | Citation | PMID / DOI | Tier | Role |
|---|---|---|---|---|
| 1 | Heller CG, Clermont Y (1963). "Spermatogenesis in man: an estimate of its duration." *Science* 140(3563):184-6. | **13953583**, DOI 10.1126/science.140.3563.184 | PRIMARY, founding | F1 anchor: intratesticular ³H-thymidine + radioautography -> 16-day cycle, ~64-day total transit. |
| 2 | Clermont Y (1963). "The cycle of the seminiferous epithelium in man." *Am J Anat* 112:35-51. | **14021715**, DOI 10.1002/aja.1001120103 | EXISTENCE-VERIFIED (no abstract indexed) | Companion founding paper, the staged histological classification underlying the 16-day figure. |
| 3 | Misell LM, Holochwost D, Boban D, Santi N, Shefi S, Hellerstein MK, Turek PJ (2006). "A stable isotope-mass spectrometric method for measuring human spermatogenesis kinetics in vivo." *J Urol* 175(1):242-6. | **16406920**, DOI 10.1016/S0022-5347(05)00053-4 | PRIMARY, **DECORRELATED (method+era)** | F1's independent replication: ²H₂O labeling + GC/MS, n=11, mean 64±8d (range 42-76). |
| 4 | Orth JM, Gunsalus GL, Lamperti AA (1988). "Evidence from Sertoli cell-depleted rats indicates that spermatid number in adults depends on numbers of Sertoli cells produced during perinatal development." *Endocrinology* 122(3):787-94. | **3125042**, DOI 10.1210/endo-122-3-787 | PRIMARY, mechanistic | F2's core anchor: araC-induced Sertoli depletion (-54%) -> spermatid depletion (-55%), confounds controlled. |
| 5 | Valdes-Socin H, Parisel A, Coppens L, Henry L, Gaspard O, Sempels M, Petrossians P (2026). "Sertoli cell-only syndrome (Del Castillo syndrome): Past, present and future." *Ann Endocrinol (Paris)* 87(1):102481. | **41485595**, DOI 10.1016/j.ando.2025.102481 | PRIMARY, **DECORRELATED (human, opposite-direction)** | F2/F5 anchor: Sertoli cells + organized tubules, zero germ cells — traces to the original 1947 Buenos Aires description. |
| 6 | Jarow JP, Chen H, Rosner TW, Trentacoste S, Zirkin BR (2001). "Assessment of the androgen environment within the human testis: minimally invasive method to obtain intratesticular fluid." *J Androl* 22(4):640-5. | **11451361** | PRIMARY quantitative, HUMAN | F3 primary anchor: percutaneous aspiration, n=21, ITT 609±50 ng/mL, ">100-fold" vs serum. |
| 7 | Coviello AD, Bremner WJ, Matsumoto AM, Herbst KL, Amory JK, Anawalt BD, Yan X, Brown TR, Wright WW, Zirkin BR, Jarow JP (2004). "Intratesticular testosterone concentrations comparable with serum levels are not sufficient to maintain normal sperm production in men receiving a hormonal contraceptive regimen." *J Androl* 25(6):931-8. | **15477366**, DOI 10.1002/j.1939-4640.2004.tb03164.x | PRIMARY quantitative, HUMAN, decorrelated cohort | F3's threshold/dose-response anchor: baseline ITT/serum ~36x; suppression to serum-comparable ITT (despite normal serum T) -> sperm count falls 98%. |
| 8 | Jarow JP, Zirkin BR (2005). "The androgen microenvironment of the human testis and hormonal control of spermatogenesis." *Ann N Y Acad Sci* 1061:208-20. | **16467270**, DOI 10.1196/annals.1336.023 | PRIMARY review | States the "100-fold" figure explicitly; also states plainly the field does not yet know the human quantitative-normal-spermatogenesis ITT threshold — quoted, not smoothed over. |
| 9 | Wang RS, Yeh S, Tzeng CR, Chang C (2009). "Androgen receptor roles in spermatogenesis and fertility: lessons from testicular cell-specific androgen receptor knockout mice." *Endocr Rev* 30(2):119-32. | **19176467**, DOI 10.1210/er.2008-0025, PMC2662628 | PRIMARY review of primary knockouts, **DECORRELATED (mouse genetics)** | F3 mechanism anchor: Sertoli-AR-KO arrests at diplotene/pre-MI; Leydig-AR-KO arrests later (round spermatid); germ-cell-AR-KO has NO effect. |
| 10 | Dym M, Fawcett DW (1970). "The blood-testis barrier in the rat and the physiological compartmentation of the seminiferous epithelium." *Biol Reprod* 3(3):308-26. | **4108372**, DOI 10.1093/biolreprod/3.3.308 | EXISTENCE-VERIFIED (no abstract indexed) | F4 founding structural anchor (Sertoli-Sertoli tight junctions, basal/adluminal compartments). |
| 11 | Sinisi AA, D'Apuzzo A, Pasquali D, Venditto T, Esposito D, Pisano G, De Bellis A, Ventre I, Papparella A, Perrone L, Bellastella A (1997). "Antisperm antibodies in prepubertal boys treated with chemotherapy for malignant or non-malignant diseases and in boys with genital tract abnormalities." *Int J Androl* 20(1):23-8. | **9202987**, DOI 10.1046/j.1365-2605.1997.00101.x | PRIMARY quantitative, HUMAN | F4 core anchor: 0/100 controls vs 26/264 (9.8%) diseased ASA+; 26/26 positives localize to a structural testicular/genital lesion. Authors' own attribution: BTB impairment. |
| 12 | Lee R, Goldstein M, Ullery BW, Ehrlich J, Soares M, Razzano RA, Herman MP, Callahan MA, Li PS, Schlegel PN, Witkin SS (2009). "Value of serum antisperm antibodies in diagnosing obstructive azoospermia." *J Urol* 181(1):264-9. | **19013620**, DOI 10.1016/j.juro.2008.09.004 | PRIMARY quantitative, HUMAN, corroborating/adjacent | F4 corroboration (excurrent-duct/vasectomy level, disclosed as anatomically distinct from the seminiferous-tubule BTB proper): IgG sensitivity 85%, specificity 97%, n=484. |
| 13 | Luetjens CM, Weinbauer GF, Wistuba J (2005). "Primate spermatogenesis: new insights into comparative testicular organisation, spermatogenic efficiency and endocrine control." *Biol Rev Camb Philos Soc* 80(3):475-88. | **16094809**, DOI 10.1017/s1464793105006755 | PRIMARY qualitative (review) | **Disclosed caution, not a positive anchor**: explicitly revises the old "human is inefficient" assumption — used to justify NOT asserting a specific human-vs-primate efficiency multiplier. |

---

## 2. F1 — the ~64-day transit is a genuine two-era, two-method over-determination

**The geometric object**: the seminiferous epithelium is not a homogeneous sheet dividing in
lockstep — at any instant, a given tubule cross-section shows one of several fixed germ-cell
*associations* (stages), because successive spermatogonial cohorts enter the pipeline on a fixed
clock and are read out along the tubule length as a **spatial wave**, not a synchronized pulse.
Heller & Clermont's 1963 method measured this directly, not by counting stages at one instant, but
by **injecting a covalent, heritable label (³H-thymidine) into the testis and tracking its most
advanced position over real time** — the label marks a spermatogonial cohort at S-phase, and
serial biopsies at increasing intervals show how far that cohort has advanced along the mitosis ->
meiosis I -> meiosis II -> spermiogenesis sequence. Result, machine-quoted from the fetched
abstract: preleptotene spermatocyte at 1 hour, midpachytene spermatocyte at 16 days, immature
spermatid at 32 days -> **one cycle of the epithelium = 16 days, total spermatogenesis ≈ 64 days**
(4 cycles, consistent with the specific germ-cell transitions being spaced at fixed cycle-length
intervals, not evenly smeared).

**Forced adversary**: a 1963 estimate, from a small sample, using a radioactive tracer that could
no longer ethically be administered to healthy volunteers today, is exactly the kind of number
this discipline flags as a possible stale, unreplicated artifact — real science has plenty of
mid-20th-century point estimates that never got an independent modern check.

**Forced to its strongest form — a genuinely different, non-invasive, non-radioactive, modern
method**: Misell et al. 2006 used oral deuterium oxide (heavy water) as a metabolic label,
detected by gas chromatography/mass spectrometry in sperm DNA, with the "cycle" operationalized as
the **lag time until labeled sperm first appear in ejaculate** — a completely different physical
measurement principle (isotope-ratio mass spec vs. autoradiography), a different era (2006 vs
1963), non-invasive (daily oral dosing + semen collection vs. intratesticular injection +
biopsy), and explicitly framed by its own authors as a replication check: *"this estimate is
derived mainly from a single older, descriptive, kinetic analysis of spermatogenesis... We
confirmed the postulated length of a normal cycle of spermatogenesis."*

**Measured**: Misell 2006, n=11, mean **64±8 days** (range 42-76). Heller & Clermont's 64-day point
estimate sits **exactly at the modern mean**, and comfortably inside the modern study's own
measured range.

**Machine gates** (`spermatogenesis_sertoli.py`, §`G1`-`G2`): relative difference between the
1963 point estimate and the 2006 mean = **0.0%** (threshold <20%) — **PASS**. 1963 estimate falls
within the 2006 study's measured range [42,76] — **PASS**.

**Accept F1**: the "stale single-study artifact" adversary FALLS — a fully decorrelated modern
method, on a different cohort, using a fundamentally different physical labeling/detection
principle, reproduces the same central number. **Honest, disclosed, not smoothed over**: real
individual variation exists (42-76 day range in the modern study) — this is not a hard biological
constant, and the modern paper's own conclusion is that the true duration sits "on the shorter
side of traditional estimates," i.e. genuine science, not a manufactured exact match. The specific
**stage count** (classically 6 in man, per Clermont 1963) is **not independently re-verified from
a live-fetched primary abstract this session** — that companion paper's abstract was not available
via efetch (title/DOI-only tier) — disclosed as an honest gap (§7), not silently upgraded from
recalled knowledge to "verified."

---

## 3. F2 — Sertoli cell number sets the ceiling; the "germ cells autonomously divide" adversary is forced and falls

**The geometric object**: each Sertoli cell physically spans the full radial thickness of the
tubule wall (basal lamina to lumen), and each germ cell at every stage is either directly attached
to or embedded in Sertoli-cell cytoplasm — there is no free-standing germ-cell compartment. If
germ-cell output were autonomously set (Sertoli support merely permissive, not load-bearing), then
perturbing Sertoli cell *number* while leaving germ cells otherwise undamaged should leave adult
output largely unaffected, or at least sub-proportional (germ cells "packing in" around whatever
Sertoli cells remain).

**Forced adversary, and forced to its strongest fair form**: the task's own named adversary,
"germ cells autonomously divide without Sertoli support." The strongest fair test is a *selective*
Sertoli-cell perturbation with the shared-toxicity confound explicitly ruled out — not a generic
whole-testis insult (which would make ANY output drop uninformative about Sertoli-specificity).
Orth, Gunsalus & Lamperti 1988 built exactly this: neonatal rat testes given cytosine arabinoside
(araC) to hit Sertoli cells (still dividing) while germ cells (mitosis not yet begun) are
protected by timing alone, with the drug itself cleared before germ-cell division starts — a
disclosed, machine-checked control, not assumed (`orth_timing_confound_ruled_out`). A second
confound — generic testicular/endocrine toxicity — is also ruled out: Leydig cell volume and
ventral-prostate weight (an androgen-bioassay proxy) were normal in the Sertoli-depleted group
(`orth_leydig_confound_ruled_out`).

**Measured**: Sertoli cell number **-54%**; round spermatid number **-55%** — the ratio of
spermatids per Sertoli cell is, in the paper's own words, "essentially identical" between normal
and depleted animals.

**Machine gates** (§`G3`-`G6`): proportionality relative deviation |54-55|/54 = **1.85%**
(threshold <20%) — **PASS**. Both confound-control booleans — **PASS**. Adversary margin: the
autonomous-adversary's own null prediction is **0%** output change (Sertoli support merely
permissive); the measured change is **55 percentage points** away from that null (threshold >30pts)
— **PASS**, a large, non-knife-edge margin.

**Decorrelated, opposite-direction anchor (human, not rat)**: Sertoli-cell-only syndrome (SCOS /
Del Castillo syndrome) is the anatomical demonstration that the dependency is **asymmetric**.
Valdes-Socin et al. 2026 (tracing the original 1947 Buenos Aires description): testicular biopsy
shows "only Sertoli cells and seminiferous tubules... no spermatozoa or spermatogonia" — i.e.
Sertoli cells CAN exist and organize a normal-looking tubule architecture with **zero** germ
cells. The reverse configuration (organized germ-cell maturation with no Sertoli cells) is not a
recognized clinical or experimental entity in this literature. Machine gate §`G7` — **PASS**
(boolean: Sertoli present AND tubules organized AND germ cells absent).

**Accept F2**: the shared-toxicant-confound adversary is explicitly ruled out at the source
(timing + Leydig-normalcy controls); the "germ cells autonomously compensate" adversary predicts
either no change or a sub-proportional change in output — it is falsified by a near-exact
proportional (54%→55%) loss, over-determined by the qualitative human SCOS asymmetry in the
opposite direction. **Honest limitation**: the quantitative leg (Orth 1988) is a rat perinatal
manipulation model — cross-species, disclosed (§7); the human leg (SCOS) is qualitative/
directional, not a matched quantitative ratio-conservation test.

---

## 4. F3 — intratesticular testosterone: ≥1 order of magnitude above serum in every measured cohort, and the operative LOCAL threshold variable — with a disclosed, forced measurement-method tension

**The geometric object**: Leydig cells (interstitial, outside the tubule) synthesize testosterone;
Sertoli cells are the androgen-receptor-bearing relay INSIDE the tubule (germ cells themselves, by
the receptor-genetics evidence below, do not sense androgen directly). For this relay to work at
the concentration germ cells actually experience, what matters is not circulating (serum) T but
the LOCAL concentration inside the seminiferous compartment — a genuinely different physical
quantity, historically hard to measure without contaminating the sample with blood or interstitial
fluid.

**Forced adversary (the task's own named suspicion) — is "100x" a real constant or an artifact of
one measurement technique?** Forced to its strongest form by comparing **two independent human
cohorts using minimally-invasive percutaneous testicular aspiration**, from an overlapping research
group but genuinely different subjects/timepoints:

- **Jarow et al. 2001** (n=21, fertile men, baseline): mean testicular-fluid T = 609±50 ng/mL,
  stated directly as **"more than 100-fold greater"** than normal serum T.
- **Coviello et al. 2004** (n=7, different cohort, baseline pre-treatment): ITT = 822±136 nmol/L
  vs. serum T = 22.8±1.9 nmol/L -> ratio = **36.05x** (machine-computed, §`G9`), the paper's own
  text rounding this to "~40x."

**Measured disagreement, disclosed not smoothed**: the two studies' own ratios differ by a factor
of **2.77x** (§`G12`, 100/36.05) — a real, method/cohort-dependent spread, not a fixed universal
constant. **This synthesis explicitly does NOT assert a single "always exactly 100x" claim.**
What DOES survive in both, independently: both cohorts clear a conservative, pre-registered
**≥10x** floor (§`G8`, §`G9` — both **PASS**), a demanding bar chosen precisely because it is not
trivially satisfied by any hormone/tissue pairing on its face.

**The threshold-for-function part, not just the ratio**: Coviello 2004's own design is a
within-subject perturbation — 6 months of testosterone-enanthate + levonorgestrel drove serum T to
normal-or-slightly-high (28.7±2.0 nmol/L, unchanged from baseline, P=.12) while **crashing ITT 98%**
to 13.1±4.5 nmol/L — a level at or *below* on-treatment serum T (ratio 0.46x, §`G10`, **PASS**
against the <2x gate — i.e. the local elevation is fully abolished, not merely reduced). Sperm
concentration fell from 65±15 to 1.3±1.3 million/mL — a **98% fall** (§`G11`, threshold ≥50%,
**PASS**) *despite serum T being normal*. This directly dissociates serum T from spermatogenic
output and identifies ITT specifically as the operative local variable.

**Mechanism anchor, fully decorrelated (mouse genetics, not human endocrinology)**: Wang et al.
2009's review of cell-specific androgen-receptor knockout mice gives the receptor-level reason
*why* the signal must route through Sertoli cells. Germ-cell-specific AR knockout: **"does not
affect spermatogenesis and male fertility"** (§`G13` — germ-cell autonomous androgen-sensing is
**FALSE**, machine-encoded). Sertoli-cell-specific AR knockout: arrest at the **diplotene primary
spermatocyte stage, prior to completion of meiosis I** (§`G14` — Sertoli relay is **required**).
Leydig-cell-specific AR knockout arrests at a **different, later** checkpoint (round spermatid
stage) — two decorrelated cell-specific phenotypes, not one smeared effect (§`G15`).

**Accept F3**: the "100x-is-an-artifact" adversary partially WINS on the exact-multiplier question
(disclosed: 36x-100x+, real spread) but FALLS on the qualitative floor (both cohorts clear ≥10x)
and on the operative-threshold question (a controlled within-subject perturbation shows ITT, not
serum T, gates output — and germ-cell AR knockouts show germ cells do not even need to sense
androgen directly, ruling out the "serum T acts directly on germ cells" alternative mechanism).
**Honest field-level limitation, quoted directly, not softened**: Jarow & Zirkin's own 2005 review
states "we do not yet know how much testosterone is required within the human testis to... restore
quantitatively normal spermatogenesis" — the precise human dose-response curve remains an open
question in the primary literature itself, not just in this synthesis.

---

## 5. F4 — the blood-testis barrier confers immune privilege; barrier-relevant lesions associate with measurable, localized antisperm-antibody formation

**The geometric object**: meiosis begins at puberty, roughly a decade after central immune
tolerance is established in infancy — haploid and recombined germ-cell surface antigens are
antigenically **novel** to the adult immune system. Dym & Fawcett's 1970 founding paper (structural
anchor, existence-verified this session, abstract not indexed — §`G20`) established that adjacent
Sertoli cells form tight junctions dividing the tubule into a **basal compartment** (spermatogonia,
early primary spermatocytes — accessible to the interstitium) and an **adluminal compartment**
(later meiotic and post-meiotic stages — sealed off). This is the structural mechanism; the
falsifiable functional CONSEQUENCE is that breaching this seal should expose the systemic immune
system to those novel antigens.

**Forced adversary**: elevated antisperm antibodies (ASA) in sick or anatomically abnormal
patients could simply be a generic autoimmune/inflammatory bystander phenomenon of being
systemically unwell — not evidence of a specific barrier breach.

**Forced to its strongest form**: Sinisi et al. 1997 studied **prepubertal** boys (ages 1.2-13,
i.e. largely BEFORE active spermatogenesis begins — a design that itself tests whether structural
genital lesions matter independent of ongoing meiotic activity), across three groups: malignant
disease (n=42) + nephrotic syndrome (n=10) = Group I; genital tract abnormalities (cryptorchidism,
inguinal hernia, funicular torsion, hypospadias, n=202) + cystic fibrosis (n=10) = Group II; n=100
normal controls = Group III. If ASA were a generic sick-child phenomenon, it should appear across
Group I as readily as Group II.

**Measured**: 0/100 (0%) controls ASA+; 26/264 (9.8%) combined diseased groups ASA+. Of the 26
positive, **24 had genital tract abnormalities** and **2 had leukemia WITH direct testicular
infiltration** — i.e., **26/26 (100%)** of all ASA+ cases had an actual local testicular/genital
structural lesion; **zero** were "systemically diseased with no local lesion." The authors'
own mechanistic attribution: *"potential impairment of the blood testis (Sertoli cell) barrier."*

**Machine gates** (§`G16`-`G19`): control rate exactly 0% — **PASS**. Diseased rate >0% —
**PASS**. Fisher exact test **computed fresh from the raw 2x2 table** (not merely quoted from the
paper's own significance claim) on `[[26,238],[0,100]]`: **p = 3.29×10⁻⁴** (threshold <0.05) —
**PASS**. Localization fraction 26/26 = **100%** (threshold ≥90%) — **PASS**, directly forcing down
the generic-bystander adversary: every single positive case has an identifiable local lesion.

**Corroborating, but honestly disclosed as mechanistically ADJACENT (not the same barrier)**:
Lee et al. 2009 (n=484 adult men) found serum IgG-ASA predicts surgically-confirmed excurrent-duct
obstruction — dominated by **vasectomy** per the authors' own conclusion — with 85% sensitivity,
97% specificity (§`G21`). This is real, decisive evidence for the same *underlying principle*
(sequestered sperm antigens -> immune sensitization if the seal breaks), but vasectomy severs
continuity at the **vas deferens**, anatomically distinct from the Sertoli-Sertoli tight junction
inside the tubule (a closely related but NOT identical "blood-epididymis barrier" mechanism) —
flagged as a distinct locus, not conflated (§7).

**Accept F4**: the generic-autoimmune-bystander adversary is forced via within-cohort structure
(malignancy-alone vs. genital-abnormality) and FALLS — the effect concentrates specifically and
completely (100%) on cases with an actual local lesion, statistically decisive (p<0.001, computed
not quoted), with zero false "background" positives in matched normal controls.

---

## 6. Results summary — 21/21 machine-computed gates PASS

```
G1  cycle_relative_diff_lt_20pct:                       PASS  (0.0%)
G2  cycle_point_within_misell_range:                    PASS  (64 in [42,76])
G3  sertoli_depletion_proportionality_lt_20pct_dev:      PASS  (1.85% deviation)
G4  sertoli_depletion_timing_confound_ruled_out:         PASS
G5  sertoli_depletion_leydig_confound_ruled_out:         PASS
G6  sertoli_adversary_margin_gt_30pts:                   PASS  (55pts vs null's 0)
G7  scos_demonstrates_asymmetric_dependency:             PASS
G8  jarow2001_ratio_gte_floor:                           PASS  (>=100x)
G9  coviello2004_baseline_ratio_gte_floor:                PASS  (36.05x)
G10 coviello2004_ontreatment_ratio_lt_2x:                PASS  (0.456x)
G11 coviello2004_sperm_fall_gte_50pct:                   PASS  (98.0%)
G12 itt_disagreement_disclosed_both_clear_floor:         PASS  (disclosed 2.77x cross-study spread)
G13 scarko_germ_cell_autonomous_sensing_falsified:       PASS  (False, as required)
G14 scarko_sertoli_relay_required:                       PASS
G15 scarko_leydig_decorrelated_later_checkpoint:         PASS
G16 btb_control_group_asa_rate_is_zero:                  PASS  (0/100)
G17 btb_diseased_group_asa_rate_gt_zero:                 PASS  (26/264)
G18 btb_fisher_exact_p_lt_0p05:                          PASS  (p=3.29e-4, computed fresh)
G19 btb_all_positives_localize_to_structural_lesion:     PASS  (26/26 = 100%)
G20 dymfawcett_citation_exists:                          PASS  (existence-only tier)
G21 lee2009_corroborating_adjacent_mechanism:             PASS  (sens 0.85, spec 0.97)
```

`overall_pass_strict_all = True`. Determinism: 2 independent process runs produce identical JSON
content (verified via `diff` after stripping the self-reported `generated_utc` timestamp, which is
expected to vary by construction). NaN/Inf scan: **one genuine catch, fixed at the source, not
hidden** — see §7.

---

## 7. Symmetric QC — what went wrong once, and what this does NOT prove

**A real, disclosed catch during this build**: the first run of the NaN/Inf scanner correctly
flagged a non-finite value — the Fisher-exact odds ratio on `[[26,238],[0,100]]` is mathematically
**+infinity**, because the control arm has zero events (division by zero in the odds-ratio
formula). This is a REAL, EXPECTED statistical fact (a perfectly clean 0% control rate), not a
computational bug — but `+inf` is not valid strict JSON. Fixed at the source: the odds ratio is
now stored as a disclosed string with an explanatory note, and the gate itself keys off the
p-value (always finite here), which is the actually decisive statistic. Reported here as an
example of the discipline in action — the scanner did its job, the fix was to make the honest
statistical fact JSON-safe, not to suppress or round away the infinity.

**Sorted by evidentiary weight**:
- **Hard, non-circular, machine-computed falsifiers** (could genuinely have failed): F1's
  cross-era/cross-method numeric agreement (§2); F2's proportionality-preservation under a
  confound-controlled perturbation (§3); F3's within-subject ITT-suppression/output-collapse
  dissociation from serum T (§4); F4's Fisher-exact test computed fresh from the raw contingency
  table, not merely quoted (§5, §6/`G18`).
- **Disclosed, forced (not hidden) disagreements**: the ITT/serum ratio is genuinely
  method/cohort-dependent (36x vs 100x+, §4/`G12`) — this synthesis explicitly refuses to assert a
  single universal "100x" constant, even though that number is the one most often quoted in
  secondary literature.
- **Existence-verified but not quantitatively re-extracted** (disclosed limitation, matching the
  precedent set by `MECHANISM_SPERM_MOTILITY_AXONEME.md`'s 3 "weaker, existence-only" citations):
  Clermont 1963 (the specific stage count, classically 6, was NOT re-verified from a live-fetched
  abstract this session — no abstract was indexed) and Dym & Fawcett 1970 (the specific
  tracer/EM evidence for the tight junctions was not re-extracted from a live-fetched abstract —
  same reason). Both citations are real and correctly attributed; their DETAILED numeric content
  is a textbook-level gap, not a fabricated one.
- **Cross-species transplants, disclosed**: Orth 1988 (rat, the quantitative Sertoli-ceiling leg)
  and Wang 2009's synthesized knockout lines (mouse, the AR-relay mechanism leg) are not human —
  the human legs (SCOS/Del Castillo, Coviello/Jarow ITT measurements, Sinisi/Lee ASA studies)
  anchor the same claims independently in the actual target species, but the mechanistic
  granularity (exact ratio-conservation; exact knockout phenotype) rests on animal models.
  Testicular architecture and the AR/Sertoli relay are conserved across mammals, but this is not
  itself independently re-verified this session for the human case at that mechanistic grain.
- **A deliberately NOT-made claim**: this synthesis explicitly declines to assert a specific
  "human spermatogenesis is Nx less efficient than other primates" multiplier — Luetjens et al.
  2005's own review argues the older version of that claim is likely wrong (multi-stage tubular
  organisation does not in fact imply low efficiency) — using a citation to justify NOT making a
  folk claim is itself disclosed as a deliberate scope choice, not an oversight.
- **A precision maintained throughout, not conflated**: the blood-testis barrier proper
  (Sertoli-Sertoli tight junctions inside the tubule) is kept distinct from the anatomically
  adjacent but different blood-epididymis barrier (relevant to the vasectomy/Lee 2009
  corroboration) — both are real, both plausibly share the same underlying "sequestered antigen"
  logic, but they are not the same physical structure and this doc does not claim they are.
- **Single confidence tier throughout**: population/literature-anchored (human clinical +
  cross-species experimental/genetic), never subject-specific — no testicular biopsy, aspirate, or
  serology panel exists for subject2, matching every sibling endocrine/reproductive doc's own
  disclosed scope (`MECHANISM_REPRODUCTIVE_HPG.md` §9, `MECHANISM_SPERM_MOTILITY_AXONEME.md` §5).
- **No graph-edge write this session** — folding into `data/MECHANISM_ANCHOR_GRAPH.json` via the
  canonical `mechanism_fold -> fold_gate_v2` path is the coordinator's next step, not performed here
  (isolation rule: touch only files created this session; the anchor graph and `bt_memory/` were
  read, never written).

---

## 8. couples_to

- **`docs/MECHANISM_REPRODUCTIVE_HPG.md`** (`MSK-HPG-RATE-CONSTANT-LAYER` et al.) — that doc's own
  output state vector (`serum_T`, `serum_FSH`, `LH_pulse_frequency_t`) is the SYSTEMIC input this
  doc's Sertoli/tubule machinery consumes; this doc's own finding (§4) that serum T is
  DISSOCIATED from output while ITT is the operative variable is a direct, disclosed refinement of
  what "the hormonal input" actually means for spermatogenesis — serum T is necessary for Leydig
  function and diagnostic convenience, but the tubule itself responds to a different, local
  concentration. No numeric re-derivation performed into that doc's own state vector this session
  (disclosed, matching that doc's own §8 scope note).
- **`docs/MECHANISM_SPERM_MOTILITY_AXONEME.md`** — this doc's production line is the upstream
  process that MAKES the flagellated cell that doc's RFT/dynein model then propels; that doc's own
  honest gap #6 names this exact axis and disclaims coverage of it — now filled.
- **`REPRO-SPERMATOGENESIS`** (`data/MECHANISM_ANCHOR_GRAPH.json`, read-only) — shares the
  organ/process but is a genuinely distinct construct (DNA-fragmentation-index / fertility-outcome
  vs. this doc's cellular-production-line mechanism); not duplicated, not merged, per direct
  content inspection this session.
- **`REPRO-AZOOSPERMIA-RETRIEVAL-GENOTYPE-GATE`** (read-only) — SCOS (§3, §5 here) is one of the
  histological causes of the non-obstructive azoospermia that node's genotype/retrieval-prognosis
  axis addresses from a different angle (which specific genetic lesion, and can sperm still be
  retrieved) — a shared clinical endpoint, different mechanistic question.
- **Mucociliary clearance / PCD cert** (`PCD-CILIARY-GENOTYPE-ULTRASTRUCTURE-PHENOTYPE-GATE`,
  read-only) and `MECHANISM_SPERM_MOTILITY_AXONEME.md`'s own dynein-machinery coupling — indirectly
  related via the shared cytoskeletal-motor biology, not directly coupled by this doc's own
  content (Sertoli/BTB mechanisms are unrelated to axonemal dynein).

---

## 9. Files

- `scripts/msk/spermatogenesis_sertoli.py` — self-contained (numpy/scipy only), the `CITATIONS`
  dict (13 live-verified entries), all raw extracted numbers, all 21 gate computations (including
  a Fisher exact test computed from the raw contingency table, not quoted), the NaN/Inf scanner,
  and the evidence-JSON writer (writes to both the repo-convention path and this task's requested
  path).
- `data/msk_smoketest/spermatogenesis_sertoli/spermatogenesis_sertoli_results.json` /
  `docs/MECHANISM_SPERMATOGENESIS_evidence.json` (identical) — full machine-readable evidence:
  citations, all 5 claims' raw numbers, all 21 gates, the verdict block. Determinism verified
  (content-identical across 2 independent process runs, timestamp field excepted); NaN/Inf-free
  after the disclosed odds-ratio fix (§7).
- Read but NOT modified (isolation: touch only files created this session):
  `data/MECHANISM_ANCHOR_GRAPH.json` (`REPRO-SPERMATOGENESIS`, `REPRO-AZOOSPERMIA-RETRIEVAL-
  GENOTYPE-GATE`, `AUTO-CELL-SPERM-DECLINE-DECORRELATED-CHECK-MA`,
  `AUTO-CELL-PATERNAL-GERMLINE-EPI-BRIDGE-A-MACH` — all confirmed distinct, not duplicated),
  `bt_memory/` (grepped for prior spermatogenesis/Sertoli learnings — zero hits, confirmed fresh),
  `docs/MECHANISM_REPRODUCTIVE_HPG.md`, `docs/MECHANISM_SPERM_MOTILITY_AXONEME.md` (couples_to
  context), `docs/MECHANISM_HARDENED_CONVENTIONS.md` (the authoritative fold-path/schema/isolation
  reference consulted before building).

## 10. Reproduction

```
cd ~/projects/bodytwin
.venv-msk/bin/python3 scripts/msk/spermatogenesis_sertoli.py
```

No args, no network access at runtime (all 13 citations independently live-verified via NCBI
eutils before this script was written, hardcoded with citations in `CITATIONS`). Runtime <1s (pure
arithmetic + one `scipy.stats.fisher_exact` call; no Monte Carlo, no OpenSim). Writes
`data/msk_smoketest/spermatogenesis_sertoli/spermatogenesis_sertoli_results.json` (repo
convention) and `docs/MECHANISM_SPERMATOGENESIS_evidence.json` (this task's requested path) plus a
printed console summary. Exit code 0 iff `overall_pass_strict_all`.

## 11. Roadmap (not attempted here — first cell only)

1. Independently re-extract Clermont 1963's specific stage count and staging criteria from the
   full paper text (only title/DOI verified this session, no abstract indexed) — §7 honest gap.
2. Independently re-extract Dym & Fawcett 1970's specific tracer/EM evidence from the full paper
   text (same reason) — §7 honest gap.
3. Resolve the disclosed ITT/serum ratio spread (36x vs 100x+) with a wider, systematically
   searched panel of independent measurement techniques (microdialysis vs. aspiration vs.
   homogenate) rather than the 2 cohorts available this session.
4. Couple a specific numeric re-derivation into `docs/MECHANISM_REPRODUCTIVE_HPG.md`'s state vector
   (e.g., an explicit serum-T -> ITT conversion/gain term), rather than the qualitative dissociation
   argument made here.
5. Fold this hypothesis into `data/MECHANISM_ANCHOR_GRAPH.json` via the canonical
   `mechanism_fold -> fold_gate_v2` path (not performed this session, isolation rule).
