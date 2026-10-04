# MECHANISM TMJ JOINT — bicondylar jaw lever mechanics: bite force + condylar reaction (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** Confidence tier: **geometrically-forced,
multi-anchor triangulated** for the ORDINAL/ratio claims (molar>incisor bite force; the
molar/incisor ratio and the balancing>working asymmetry are proven algebraically
independent of muscle-force magnitude, not fitted) — **DOMAIN_CONSENSUS/estimated** for
absolute-force magnitudes (several geometric inputs are anatomical-range estimates, not
literature-pinned this session; disclosed in full in §6).

## Headline result

| quantity | model (machine-solved) | measured anchor | verdict |
|---|---:|---:|---|
| Molar/incisor max bite-force ratio | **2.15** (central), sweep range **[1.60, 2.97]** | **2.38** (men), **2.39** (women) — Waltimo & Könönen 1995 | **PASS** (range contains both) |
| ...and that ratio is **provably muscle-force-independent** | `ratio_from_forces` − `d_incisor/d_molar` = **4.4×10⁻¹⁶** (machine-verified identity, n=4000) | geometric necessity, not a fit | **PASS** |
| Condylar reaction force, molar bite | central **1603 N**, sweep **[201, 3885] N** — always **> 0** | Hylander 1975: joint reaction is real, substantial, condylar neck built to withstand it | **PASS** (falsifies the "R=0" adversary) |
| Condylar reaction force, incisor bite | central **2305 N**, sweep **[489, 4800] N** | same qualitative anchor | **PASS** |
| Joint-load-per-unit-bite-force, incisor vs molar | **3.78** vs **1.22** (central) — incisor **3.1×** less efficient, **100%** of sweep | "molar bite → lower joint load per unit bite force" (task's own stated physiology) | **PASS** |
| Balancing-side > working-side condylar reaction (unilateral molar bite) | **100%** of sweep (n=1500), even up to 99.9% working-side muscle-force bias | Hylander 1975 direct EMG-based conclusion; Shi et al. 2012 macaque multibody-dynamics cross-check | **PASS** |
| Interincisal opening (rotation-chord + measured translation) | sweep **[40.5, 76.3] mm**, central **57.9 mm** | **46.0–54.8 mm** (Lewis/Buschang/Throckmorton 2001); **49.9–50.3 mm** (Agrawal 2015) | **PASS** (range overlaps; only 19% of grid points land inside the tight band — disclosed, §6) |
| Muscle-weakness adversary | 19.67% PCSA drop → **exactly** 19.67% bite-force drop (linear identity, 1.1×10⁻¹⁴) | direction matches van Spronsen 1992 (smaller jaw-muscle CSA ↔ smaller max molar bite force) | **PASS** |
| Disc-displacement (delayed-translation) adversary | opening deficit **0→12 mm** as delay_frac 0→0.75, strictly monotonic | Kalaykova et al. 2010: moment-of-disc-reduction shifts later in opening as locking develops (directional match) | **PASS** |

**10/10 machine gates PASS** (`scripts/msk/tmj_lever_model.py`, exit 0). Every equilibrium is
solved as an explicit linear system (`numpy.linalg.solve`, cross-checked against a vectorized
closed-form to **0.0 N** max difference), never hand algebra — see §4 for a case where an
initial **hand-derived** estimate (a "crossover at α≈0.94") was **wrong**, caught only by
letting the code compute the real answer (no crossover at all up to α=0.999).

## 0. What this is / is not

This is a **first-principles static-equilibrium lever model** of the bicondylar jaw,
analogous in spirit to this repo's existing `MECHANISM_MOMENT_ARM_VALIDATION.md` (geometry
decorrelated from force) and `MECHANISM_JOINT_FORCE_VALIDATION.md` (limb joint reaction
force), but built from scratch in plain Python/numpy rather than OpenSim — **no OpenSim
model in this repo includes a mandible/jaw** (verified: `find data/msk_models -iname
"*.osim"` → 20 models, all lower-limb/arm `LaiArnoldModified2017`-lineage, zero
masseter/temporalis/jaw joints). This is the right scope for a first TMJ cert: a lean,
auditable, three-part rigid-body statics model, not a claim to have built a full
3D/dynamic/FE jaw twin.

It is **not** a fit: no parameter in this model was tuned to hit the measured bite-force
ratio or asymmetry. The two headline ordinal claims (molar/incisor ratio; balancing>working)
are **proven algebraically** to be independent of the muscle-force-magnitude parameters
(PCSA, specific tension, `r_muscle`) — verified numerically to machine precision (§3), not
merely argued. Absolute force magnitudes DO depend on those parameters and are reported at a
visibly weaker confidence tier (§6).

## 1. Literature anchors — all live NCBI-eutils-verified this session (DOI/PMID cross-checked)

WebSearch was unavailable (session-wide quota exhausted, consistent with the pre-existing
constraint already logged in `docs/MECHANISM_MUSCLE_AUDIT.md`/`MECHANISM_ANKLE_FORCE.md`).
Verification used `curl` to NCBI eutils (`esearch`/`esummary`/`efetch`) directly against
PubMed, plus Europe PMC's REST API and one PMC full-text fetch — every number below is
machine-fetched abstract/table text with DOI+PMID confirmed, not narrated recall.

| # | Citation | DOI / PMID | What it gives this model |
|---|---|---|---|
| 1 | Waltimo A, Könönen M. 1995. *Acta Odontol Scand* 53(4):254-8 | [10.3109/00016359509005982](https://doi.org/10.3109/00016359509005982) / PMID [7484109](https://pubmed.ncbi.nlm.nih.gov/7484109/) | **PRIMARY bite-force anchor.** n=129. Molar/incisor: men 909±177N / 382±133N; women 777±168N / 325±116N. Verbatim: *"Mean maximal bite force values for men were 909 N... in the molar region and 382 N... in the incisal region."* |
| 2 | Waltimo A, Könönen M. 1993. *Scand J Dent Res* 101(3):171-5 | [10.1111/j.1600-0722.1993.tb01658.x](https://doi.org/10.1111/j.1600-0722.1993.tb01658.x) / PMID [8322012](https://pubmed.ncbi.nlm.nih.gov/8322012/) | 2nd independent device/cohort (n=30): molar 847N(M)/597N(F) — same lab, different recorder, same qualitative regime. |
| 3 | Hylander WL. 1975. *Am J Phys Anthropol* 43(2):227-42 | [10.1002/ajpa.1330430209](https://doi.org/10.1002/ajpa.1330430209) / PMID [1101706](https://pubmed.ncbi.nlm.nih.gov/1101706/) | **PRIMARY asymmetry anchor.** Verbatim: *"during powerful unilateral molar biting the resultant adductor muscle force is passing between the bite point and the balancing (non-biting side) condyle... reaction forces are larger on the balancing side than on the working side."* Also: condylar neck is strong enough to withstand reaction force (falsifies "mandible = pure link, no joint reaction" hypothesis). |
| 4 | Gingerich PD. 1979. *Am J Phys Anthropol* 51(1):135-7 | [10.1002/ajpa.1330510116](https://doi.org/10.1002/ajpa.1330510116) / PMID [453343](https://pubmed.ncbi.nlm.nih.gov/453343/) | Commentary refining Hylander's 1978 incisor-bite data: mandible functions as **both** lever and link simultaneously during incisal biting — nuance disclosed, not flattened. |
| 5 | Osborn JW, Baragar FA. 1992. *J Biomech* 25(9):967-74 | [10.1016/0021-9290(92)90032-v](https://doi.org/10.1016/0021-9290(92)90032-v) / PMID [1517273](https://pubmed.ncbi.nlm.nih.gov/1517273/) | Direct methodological precedent: linear-programming static-equilibrium model, joint-force **direction constrained to be compression-only** (normal to the condylar surface) at bite forces 100–1000N, incisor/premolar/1st/3rd molar. Independently corroborates this model's own emergent tensile/compression subtlety (§4). |
| 6 | Pruim GJ, de Jongh HJ, ten Bosch JJ. 1980. *J Biomech* 13(9):755-63 | [10.1016/0021-9290(80)90237-7](https://doi.org/10.1016/0021-9290(80)90237-7) / PMID [7440590](https://pubmed.ncbi.nlm.nih.gov/7440590/) | Classic static bite-force + joint-reaction model. **Existence/citation-count verified via Europe PMC (162 citations, real DOI)**; abstract text unavailable (pre-abstract-era MEDLINE record) — cited as existence-verified only, not quoted. |
| 7 | van Eijden TM, Korfage JA, Brugman P. 1997. *Anat Rec* 248(3):464-74 | [10.1002/(sici)1097-0185(199707)248:3<464::aid-ar20>3.3.co;2-4](https://doi.org/10.1002/(sici)1097-0185(199707)248:3%3C464::aid-ar20%3E3.3.co;2-4) / PMID [9214565](https://pubmed.ncbi.nlm.nih.gov/9214565/) | n=8 cadavers. Jaw-closers (masseter/temporalis/med.pterygoid) architecturally suited for FORCE (larger PCSA, shorter fibers, higher pennation); jaw-openers suited for velocity. |
| 8 | van Eijden TM, Koolstra JH, Brugman P. 1996. *Anat Rec* 246(4):565-72 | [10.1002/(SICI)1097-0185(199612)246:4<565::AID-AR17>3.0.CO;2-M](https://doi.org/10.1002/(SICI)1097-0185(199612)246:4%3C565::AID-AR17%3E3.0.CO;2-M) / PMID [8955797](https://pubmed.ncbi.nlm.nih.gov/8955797/) | n=8 cadavers, 6 AP temporalis portions, PCSA range **1.82–2.93 cm²/portion**. |
| 9 | van Eijden TM, Koolstra JH, Brugman P. 1995. *J Dent Res* 74(8):1489-95 | [10.1177/00220345950740080901](https://doi.org/10.1177/00220345950740080901) / PMID [7560404](https://pubmed.ncbi.nlm.nih.gov/7560404/) | n=8 cadavers. Medial pterygoid PCSA (ant 2.47±0.57 + post 3.53±0.97) = **6.00 cm² direct sum**, used as-is. |
| 10 | Koolstra JH, van Eijden TM. 2001. *J Biomech* 34(9):1179-88 | [10.1016/s0021-9290(01)00053-7](https://doi.org/10.1016/s0021-9290(01)00053-7) / PMID [11506788](https://pubmed.ncbi.nlm.nih.gov/11506788/) | Masticatory system is kinematically/mechanically **indeterminate** (more muscles than DOF) — methodological precedent for treating this system with explicit linear-system solves. |
| 11 | Koolstra JH, van Eijden TM. 1995. *J Dent Res* 74(9):1564-70 | [10.1177/00220345950740091001](https://doi.org/10.1177/00220345950740091001) / PMID [7560417](https://pubmed.ncbi.nlm.nih.gov/7560417/) | Dynamic 6-DOF model: "swing-slide" condylar movement along the articular eminence generated by masseter + medial pterygoid — mechanistic confirmation of the rotate-then-translate kinematics. |
| 12 | van Spronsen PH et al. 1989. *J Dent Res* 68(12):1765-70 | [10.1177/00220345890680120901](https://doi.org/10.1177/00220345890680120901) / PMID [2600258](https://pubmed.ncbi.nlm.nih.gov/2600258/) | n=12. Masseter + medial-pterygoid CSA correlate significantly with max voluntary bite force; **temporalis CSA did not** — disclosed nuance. |
| 13 | van Spronsen PH et al. 1992. *J Dent Res* 71(6):1279-85 | [10.1177/00220345920710060301](https://doi.org/10.1177/00220345920710060301) / PMID [1613176](https://pubmed.ncbi.nlm.nih.gov/1613176/) | **Muscle-weakness adversary anchor.** Long-face (n=13) vs normal (n=35) adults: masseter/med-pterygoid/ant-temporalis CSA **30%/22%/15% smaller**, with correspondingly smaller max molar bite force (direction only; no inline %). |
| 14 | Weijs WA, Hillen B. 1985. *Acta Morphol Neerl Scand* 23(3):267-74 | PMID [4096273](https://pubmed.ncbi.nlm.nih.gov/4096273/) | Classic jaw-muscle CSA/strength paper. Existence-verified; abstract text unavailable (pre-abstract-era record). |
| 15 | Daboul A et al. 2018. *J Nutr Health Aging* 22(7):829-836 | [10.1007/s12603-018-1029-1](https://doi.org/10.1007/s12603-018-1029-1) / PMID [30080228](https://pubmed.ncbi.nlm.nih.gov/30080228/), PMCID [PMC12880510](https://pmc.ncbi.nlm.nih.gov/articles/PMC12880510/) | **Masseter PCSA anchor.** n=747 (SHIP cohort). Masseter CSA, 30–39y: **5.10±0.93 cm²(M) / 3.90±0.78 cm²(F)** (MRI single-slice anatomical CSA — a likely-conservative proxy for full pennation-corrected PCSA, disclosed §6). |
| 16 | Agrawal J et al. 2015. *Indian J Dent Res* 26(4):361-5 | [10.4103/0970-9290.167638](https://doi.org/10.4103/0970-9290.167638) / PMID [26481881](https://pubmed.ncbi.nlm.nih.gov/26481881/) | n=500. Max mouth opening **50.3±6.26mm(M) / 49.9±6.74mm(F)** — independent ROM cross-check. |
| 17 | Kaboosaya B. 2026. *Arch Oral Biol* 185:106554 | [10.1016/j.archoralbio.2026.106554](https://doi.org/10.1016/j.archoralbio.2026.106554) / PMID [41740348](https://pubmed.ncbi.nlm.nih.gov/41740348/) | n=298, ages 17-85. Very recent independent MIO reference-range cross-check. |
| 18 | Lewis RP, Buschang PH, Throckmorton GS. 2001. *Am J Orthod Dentofacial Orthop* 120(3):294-303 | [10.1067/mod.2001.115612](https://doi.org/10.1067/mod.2001.115612) / PMID [11552129](https://pubmed.ncbi.nlm.nih.gov/11552129/) | **PRIMARY kinematics anchor.** n=56. Incisor opening 52.1mm(M)/46.0mm(F); condylar translation **15.4–17.6mm(M) / 12.4–12.7mm(F)** straight-line (20.5–20.7 / 16.2–17.9mm curvilinear). Explicitly: condylar translation does **not** correlate with incisor opening (disclosed nuance, not oversimplified). |
| 19 | Chen X. 1998. *Am J Phys Anthropol* 106(1):35-46 | [10.1002/(SICI)1096-8644(199805)106:1<35::AID-AJPA3>3.0.CO;2-C](https://doi.org/10.1002/(SICI)1096-8644(199805)106:1%3C35::AID-AJPA3%3E3.0.CO;2-C) / PMID [9590523](https://pubmed.ncbi.nlm.nih.gov/9590523/) | Real jaw opening is **simultaneous** rotation+translation, not a sharp two-phase switch (rotation "somewhat more significant" only in the first 10°) — the textbook two-phase picture is an idealization; folded in as an explicit, disclosed limitation of Part C (§5). |
| 20 | Kalaykova S, Lobbezoo F, Naeije M. 2010. *J Orofac Pain* 24(4):373-8 | PMID [21197509](https://pubmed.ncbi.nlm.nih.gov/21197509/) | **Disc-displacement adversary anchor.** n=55. Moment-of-disc-reduction (MDR) shifts to **later** mouth opening as intermittent locking develops — the directional anchor for the delayed-translation adversary. |
| 21 | Marpaung CM, Kalaykova SI, Lobbezoo F, Naeije M. 2014. *J Oral Rehabil* 41(4):243-9 | [10.1111/joor.12130](https://doi.org/10.1111/joor.12130) / PMID [24533784](https://pubmed.ncbi.nlm.nih.gov/24533784/) | n=53 (referred/clinical sample, NOT general population — disclosed). ADDR: 27.6% clinical exam / 15.2% movement recording / **44.8% MRI**. |
| 22 | Naeije M et al. 2013. *J Oral Rehabil* 40(2):139-58 | [10.1111/joor.12016](https://doi.org/10.1111/joor.12016) / PMID [23199296](https://pubmed.ncbi.nlm.nih.gov/23199296/) | **Systematic review, general population.** Disc displacement prevalence **18–35%**; mostly stable/pain-free/lifelong ("noisy annoyance"); rarely progresses to painful closed lock. |
| 23 | Nickel JC, Iwasaki LR, Walker RD, McLachlan KR, McCall WD Jr. 2003. *J Dent Res* 82(3):212-7 | [10.1177/154405910308200312](https://doi.org/10.1177/154405910308200312) / PMID [12598551](https://pubmed.ncbi.nlm.nih.gov/12598551/) | n=6, EMG-validated. Muscle-force recruitment during static biting is consistent with minimization-of-joint-load OR minimization-of-muscle-effort, depending on individual/bite location — the real system is an OPTIMIZATION over an indeterminate space; this model uses the simpler determinate reduction (single resultant), disclosed as a simplification (§6). |
| 24 | Shi J, Curtis N, Fitton LC, O'Higgins P, Fagan MJ. 2012. *J Theor Biol* 310:21-30 | [10.1016/j.jtbi.2012.06.006](https://doi.org/10.1016/j.jtbi.2012.06.006) / PMID [22721994](https://pubmed.ncbi.nlm.nih.gov/22721994/) | **Cross-species, cross-method corroboration.** Macaque multibody-dynamics + optimization model predicts working:balancing muscle-force ratios and TMJ reaction ratios "comparable to those observed in vivo," peak bite forces matching published data — independent confirmation of the balancing>working finding via a completely different method (optimization-based multibody dynamics vs this model's determinate statics) and species. |
| 25 | Throckmorton GS, Finn RA, Bell WH. 1980. *Am J Orthod* 77(4):410-20 | [10.1016/0002-9416(80)90106-2](https://doi.org/10.1016/0002-9416(80)90106-2) / PMID [6928742](https://pubmed.ncbi.nlm.nih.gov/6928742/) | Direct precedent: 2D mechanical-advantage model of masseter+temporalis, used to compare facial-morphology groups — same model class as Part A. |
| 26 | Finn RA, Throckmorton GS, Bell WH, Legan HL. 1980. *J Oral Surg* 38(4):257-64 | PMID [6928454](https://pubmed.ncbi.nlm.nih.gov/6928454/) | Same 2D mechanical-advantage model, clinical application (surgical correction of mandibular deficiency changes mechanical advantage by 13–21%). |
| 27 | Koc D, Dogan A, Bek B. 2010. *Eur J Dent* 4(2):223-32 | PMID [20396457](https://pubmed.ncbi.nlm.nih.gov/20396457/), PMCID [PMC2853825](https://pmc.ncbi.nlm.nih.gov/articles/PMC2853825/) | Review: bite-force measurement is confounded by pain/TMD/sex/age/craniofacial morphology/occlusal factors/device — the methodological-heterogeneity caveat behind why absolute bite-force numbers vary 2-3× across studies. |

## 2. Model — three sub-models + two perturbation adversaries, all machine-solved

Script: `scripts/msk/tmj_lever_model.py` (self-contained, no OpenSim dependency — pure
numpy). Run: `python3 scripts/msk/tmj_lever_model.py`, exit 0. Every equilibrium is an
explicit 2-equation/2-unknown linear system solved via `numpy.linalg.solve` (never hand
algebra); a vectorized closed-form fast-path is cross-checked against the reference solver
to **0.0 N** max difference (`verify_vectorization()`) before being used in the sweep.

**A. Sagittal single-fulcrum model** (bilateral symmetric clench). Free body: mandible
pivoting at the condyle. Unknowns `[F_bite, R_condyle]`; equations `[moment about condyle,
vertical force]`. Elevator force = `(masseter + temporalis + medial pterygoid) PCSA ×
specific tension`, bilateral. Molar/incisor bite point at `d_molar`/`d_incisor` from the
condyle. **Key algebraic fact, verified not assumed**: `F_bite = F_muscle·r_muscle/d_bite`,
so `F_bite(molar)/F_bite(incisor) = d_incisor/d_molar` **exactly** — `F_muscle` and
`r_muscle` cancel. Verified numerically across n=4000 sweep samples: max
`|ratio_predicted − ratio_from_forces|` = **4.4×10⁻¹⁶** (machine epsilon). `R_condyle =
F_muscle − F_bite`, which is what the "fixed-fulcrum-ignoring-the-joint-reaction" adversary
implicitly sets to zero — this model computes it explicitly and shows it is never zero
(gate F3).

**B. Coronal two-condyle model** (unilateral molar clench). Free body: mandible supported
at two condyles (`x=0` balancing, `x=w` working), carrying the upward muscle resultant and
the downward working-side bite reaction. Unknowns `[R_L, R_R]`; equations `[vertical force,
moment about the balancing condyle]`. The working-side muscle-force fraction `α` is swept
`[0.5, 0.999]` as the forced adversary: can enough working-side EMG bias make the working
condyle's reaction exceed the balancing condyle's? `find_alpha_critical()` does a dense
(20000-point) vectorized root search per geometry sample — **not** a hand-derived estimate
(see §4 for a caught hand-algebra error).

**C. Kinematic two-phase ROM model.** Pure-rotation chord length `2·L·sin(θ/2)` about a
fixed condyle, plus a directly-measured translation distance (Lewis et al. 2001), summed —
an explicitly idealized additive decomposition (Chen 1998's finding that real motion is
simultaneous, not sequential, is disclosed as the reason this idealization runs high, §6).

**Perturbation adversary 1 — muscle weakness.** Applies van Spronsen (1992)'s real measured
CSA reductions (masseter −30%, medial pterygoid −22%, anterior temporalis −15%) to Part A's
PCSA inputs, unchanged code path. Because `F_muscle` is an exactly-linear sum of
`PCSA_i × σ` terms, %bite-force drop must equal the PCSA-weighted-average %drop to machine
precision — verified: **19.67% vs 19.67%**, diff 1.1×10⁻¹⁴.

**Perturbation adversary 2 — disc displacement (delayed translation).** A `delay_frac`
parameter withholds a fraction of the normal translation distance, reusing
`kinematic_opening_mm()` unchanged. Monotonicity (more delay → more opening deficit) is
machine-verified, not assumed. Explicitly **not** a claim to reproduce the discrete click
event itself (§6).

**Perturbation adversary 3 — unilateral asymmetric condylar loading = Part B itself** (the
task's third named adversary is the coronal model's core subject, not a separate add-on).

## 3. Machine-checked gates (PASS/FAIL, not eyeballed)

| gate | pre-registered threshold | measured | verdict |
|---|---|---:|---|
| F1 ratio band | predicted ratio sweep range contains measured 2.38/2.39 | range [1.598, 2.966] | **PASS** |
| F1b identity | `ratio_predicted` == `ratio_from_forces` to <1e-9 | max diff 4.4×10⁻¹⁶ | **PASS** |
| F2 absolute regime | central F_bite within (0.3, 3.0)× measured (men) | molar ratio 1.44×, incisor ratio 1.60× | **PASS** |
| F3 nonzero reaction | R_condyle > 0 N, both bite positions, 100% of sweep | min 201N (molar) / 489N (incisor) | **PASS** |
| F4 efficiency order | load/bite (incisor) > load/bite (molar), 100% of sweep | median 3.77 vs 1.21 | **PASS** |
| F5 balancing>working | \|R_L\|>\|R_R\| at α=0.5, 100% of sweep | 1500/1500 | **PASS** |
| F6 crossover margin | crossover α (if any) > 0.80 for ≥95% of sweep | **0/1500 samples found any crossover** in [0.5, 0.999] | **PASS** |
| F7 ROM overlap | predicted range overlaps measured [46.0, 52.1]mm | range [40.5, 76.3]mm | **PASS** (only 19.4% of grid points land strictly inside the tight band — disclosed) |
| F8 weakness linearity | %bite-force drop == %PCSA drop, <1e-9 | 19.67% vs 19.67%, diff 1.1×10⁻¹⁴ | **PASS** |
| F9 disc monotonicity | opening deficit strictly non-decreasing with delay | 0→4→8→12mm | **PASS** |

**10/10 PASS.**

## 4. Key findings — including one caught hand-algebra error (forced by the discipline itself)

1. **The molar/incisor bite-force ratio is a geometric necessity, not a fitted parameter.**
   `F_bite = F_muscle·r_muscle/d_bite` means the ratio between any two bite positions equals
   the inverse ratio of their distances from the condyle, full stop — independent of how
   strong the muscles are. This is the cleanest possible answer to the task's own
   requirement ("must emerge from the lever geometry... not be assumed"): it is not merely
   consistent with that requirement, it is **mathematically forced** by it.

2. **The condylar reaction is large and directionally informative, falsifying the
   "ignore-the-joint" adversary structurally, not just empirically.** Because `r_muscle <
   d_bite` always (a lever operating at a mechanical disadvantage — the textbook-known
   reason jaw muscles must pull several-fold harder than the resulting bite force), `R =
   F_muscle − F_bite > 0` is **forced by the same inequality**, not an empirical add-on.
   The incisor position (larger `d_bite`) always shows a worse (larger) joint-load-per-bite
   ratio than molar — reproducing "molar bite is more joint-efficient" from the same single
   inequality, at no extra cost.

3. **A genuinely forced adversary surfaced a real physical subtlety: the unconstrained
   coronal solve predicts a *negative* (tensile) working-side reaction in 100% of the
   sweep at α=0.5.** A real TMJ is a compression-only contact joint (it cannot pull) —
   Osborn & Baragar (1992) impose exactly this constraint in their own model. A negative
   `R_R` here signals the true (contact-constrained) solution would clamp `R_R` to 0 and
   push the entire remainder onto `R_L` — i.e. the *unconstrained* linear solve, if
   anything, **understates** the true balancing/working asymmetry. This is disclosed, not
   hidden inside an absolute-value comparison (see `gates.F5_balancing_gt_working.note` in
   the evidence JSON). The exact magnitude is sensitive to un-pinned transverse geometry
   (§6); the sign/ordinal claim is not.

4. **A hand-derived estimate was wrong; the code caught it.** Working the coronal-model
   algebra by hand mid-build suggested a finite crossover near α≈0.94 (working-side force
   fraction at which the asymmetry should flip). The actual dense vectorized sweep found
   **zero crossovers in [0.5, 0.999] across all 1500 geometry samples** — the hand estimate
   had the direction of the α-dependence backwards (increasing working-side dominance pulls
   the muscle resultant *toward* the bite point, which *increases* the far-support reaction
   demand, not decreases it). This is reported explicitly per this task's own instruction
   never to trust hand algebra over a machine check — and it is a real illustration of why.

## 5. Perturbation adversaries (task's own three, addressed explicitly)

- **TMJ disc displacement (clicking/locking).** Modeled as a delayed-translation fraction
  in Part C: opening deficit grows monotonically with delay (0→12mm over delay_frac
  0→0.75), directionally matching Kalaykova et al. (2010)'s clinical finding that the
  moment-of-disc-reduction shifts later in opening as intermittent locking develops.
  Naeije et al. (2013)'s systematic review anchors the population picture: 18–35%
  prevalence, but disc displacement is mostly a benign, pain-free, lifelong "noisy
  annoyance" — only rarely progressing to a painful closed lock. **Explicitly not** a
  reproduction of the discrete click event (a real reciprocal click is a sudden
  jump/snap-through as the condyle passes the posterior band, not a smooth delay) —
  flagged as a simplified mechanistic proxy, §6.
- **Muscle-weakness reduced bite force.** §2/§4 above — exact linear identity, direction
  matches van Spronsen (1992)'s long-face-adult finding (smaller jaw-muscle CSA ↔ smaller
  measured max molar bite force), magnitude not independently pinned (van Spronsen's own
  abstract gives CSA %, not bite-force %).
- **Unilateral-chewing asymmetric condylar loading.** This IS Part B — balancing > working
  in 100% of the swept geometry, with no crossover found even under an extreme (99.9%)
  working-side muscle-force bias, corroborated by both a 50-year-old EMG-based primary
  source (Hylander 1975) and a modern cross-species multibody-dynamics-optimization model
  (Shi et al. 2012, macaque) — two decorrelated methods, two eras, two species, same
  qualitative conclusion.

## 6. Honest gaps / caveats (full list, symmetric QC)

1. **Four geometric nuisance parameters are anatomical-RANGE estimates, not literature-
   pinned this session**: `r_muscle` (elevator resultant moment arm about the condyle,
   swept 15–30mm), intercondylar width `w` (90–120mm), the transverse bite-point/muscle-
   attachment fractions in Part B, and the pure-rotation angle `θ_rot` (15–30°). The
   ordinal claims (§4.1-4.2) don't need these to be exact (they cancel or the sign is
   robust across the whole range); the ABSOLUTE force/ROM numbers do, and are reported at a
   correspondingly weaker (`DOMAIN_CONSENSUS`-like) tier.
2. **Central absolute bite-force estimate runs ~44–60% high** vs measured (gate F2 passes
   on a generous 0.3–3.0× band, the same convention already used in this repo's
   `MECHANISM_ANKLE_FORCE.md`, not a tight match). Driven mainly by the specific-tension
   choice (this model reuses this repo's own already-verified but LIMB-muscle-sourced
   15–100 N/cm² range/convention from `MECHANISM_SPECIFIC_TENSION.md` — a disclosed
   cross-domain proxy, not a jaw-specific verified number this session).
3. **Temporalis total PCSA is an ESTIMATE** (mean-of-reported-range × 6 portions = 14.25
   cm²), not a literature-reported sum — van Eijden (1996)'s abstract gives only the
   per-portion range (1.82–2.93 cm²), not the 6 individual portion values.
4. **Masseter uses an MRI anatomical cross-sectional area** (Daboul 2018, 5.10 cm²), not a
   fiber-architecture-corrected PCSA — likely a modest underestimate for this pennate
   muscle (PCSA ≥ anatomical CSA by construction for any pennation angle > 0).
5. **The kinematic model (Part C) is an idealized additive decomposition that runs
   high**: central/median predicted opening (57.9–58.5mm) exceeds the ~46–52mm measured
   central tendency; only 19.4% of the swept grid lands inside the tight measured band,
   though the full range does overlap it (satisfying the pre-registered but weaker overlap
   gate). Chen (1998)'s finding that real jaw opening is simultaneous rotation+translation,
   not cleanly separable sequential phases, is the mechanistic reason: summing a
   pure-rotation chord and a pure-translation distance double-counts some shared
   displacement. A genuinely simultaneous (not additive) kinematic model is flagged as a
   next step (§7).
6. **The working-side-tensile finding's exact magnitude is model-sensitive** (§4.3) — the
   ordinal claim (balancing>working) is robust to the un-pinned transverse geometry, but
   whether the true working-side reaction is "small and compressive" or "genuinely near
   zero/clamped" depends on parameters not independently verified this session.
7. **Two classic, highly-cited sources (Pruim 1980, 162 citations; Weijs & Hillen 1985)
   are existence-verified only** (DOI/PMID/citation-count confirmed via Europe PMC) —
   their abstract text is unavailable (pre-abstract-era MEDLINE records), so they are cited
   as corroborating context, never quoted or used as a numeric source.
8. **Single maximal-static-clench scenario throughout** — no dynamic chewing-cycle timing,
   no EMG-driven muscle recruitment optimization (Nickel et al. 2003 shows the real system
   optimizes over an indeterminate muscle set; this model uses a single determinate
   resultant, a disclosed simplification, not a claim that real recruitment is this simple).
9. **Marpaung/Kalaykova/Naeije (2014)'s 44.8% MRI disc-displacement prevalence is a
   referred/clinical sample (n=53), not general population** — disclosed distinctly from
   Naeije (2013)'s general-population 18–35% review figure; the two must not be conflated.
10. **No disc/cartilage-contact mechanics modeled here** — this doc is the LEVER/FORCE
    layer; contact pressure/pressure-distribution at the disc-condyle interface is a
    distinct, not-yet-built leg (couples to `MECHANISM_CARTILAGE_CONTACT.md`'s existing
    tibiofemoral/patellofemoral pattern, not yet extended to the TMJ disc).

## 7. Next steps

1. Pin `r_muscle`/`w`/transverse fractions from a real 3D coordinate source (a digitized
   CT/CBCT mandible or an open craniofacial multibody model) rather than swept anatomical
   ranges — would tighten §6.1-2 from `DOMAIN_CONSENSUS` toward `VERIFIED_QUANT`.
2. Replace Part B's unconstrained linear solve with a proper compression-only
   (non-negative reaction / linear-complementarity) contact formulation — resolves §4.3/§6.6
   exactly instead of flagging the sign as informative-but-magnitude-uncertain.
3. Replace Part C's additive rotation+translation idealization with a genuinely
   simultaneous (single time-parameterized) rigid-body kinematic chain, informed by Chen
   (1998)'s instantaneous-center-of-rotation data — should shrink the current ~58mm-vs-~49mm
   central-tendency gap (§6.5).
4. Extend to disc/cartilage contact-pressure mechanics at the condyle-disc-fossa interface
   (couples to `MECHANISM_CARTILAGE_CONTACT.md`'s existing pattern) — the disc-displacement
   adversary here is kinematic-only; a contact-pressure leg would let a future cert test
   whether disc displacement is also a *pressure* (not just a translation-timing) failure
   mode.

## Files

- `scripts/msk/tmj_lever_model.py` — the full pipeline (self-contained, no OpenSim
  dependency; re-runnable, seeded RNG for reproducibility, exit 0).
- `data/msk_smoketest/tmj_lever_model/tmj_lever_model_results.json` — every gate, every
  central estimate, every sweep summary, machine-written.
- `data/msk_smoketest/tmj_lever_model/sagittal_sweep_raw.csv` (n=4000),
  `coronal_sweep_raw.csv` (n=1500), `kinematic_sweep_raw.csv` (n=36) — full per-sample audit
  trail behind every range quoted above.
- `docs/MECHANISM_TMJ_JOINT_evidence.json` — this doc's structured evidence summary.

## Couples to

- **PCSA/muscle-force-capacity certs** (`MECHANISM_FMAX_PCSA_VALIDATION.md`,
  `MECHANISM_SPECIFIC_TENSION.md`) — this doc reuses the repo's own already-verified
  specific-tension range/convention as a disclosed cross-domain (limb→jaw) proxy (§6.2).
- **Moment-arm/lever mechanics** (`MECHANISM_MOMENT_ARM_VALIDATION.md`,
  `MECHANISM_JOINT_FORCE_VALIDATION.md`) — same "geometry decorrelated from force,
  cadaveric/anatomical-range tolerance-band" methodology, extended here to a joint with no
  existing OpenSim representation.
- **Cartilage-contact mechanics** (`MECHANISM_CARTILAGE_CONTACT.md`) — this doc is the
  lever/reaction-FORCE layer; the disc-condyle contact-PRESSURE layer is an open next step
  (§7.4), following the same tibiofemoral/patellofemoral pattern already built for the knee.
- **Cervical/postural chain** — the mandible's muscular sling (masseter/temporalis) and
  hyoid-suspended infrahyoid/suprahyoid muscles mechanically couple jaw posture to cervical
  spine loading; not modeled here, flagged as an open coupling, not fabricated.
