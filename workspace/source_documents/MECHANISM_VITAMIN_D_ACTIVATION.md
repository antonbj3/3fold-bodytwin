# MECHANISM VITAMIN D SYNTHESIS + ACTIVATION CASCADE — the two-hydroxylation endocrine pathway (2026-07-22)

Builds and MEASURES the cutaneous-UVB-to-calcitriol pathway: 7-dehydrocholesterol (7-DHC) + UVB
(295-300nm) → previtamin D3 → [thermal isomerization] → cholecalciferol (D3) → [liver CYP2R1] →
25(OH)D (calcidiol, the STATUS marker) → [kidney CYP27B1] → 1,25(OH)2D (calcitriol, the ACTIVE
hormone) → [CYP24A1] → inactivation. Script: `scripts/msk/vitamin_d_activation.py`. Evidence:
`data/vitamin_d_activation/vitamin_d_activation_results.json`.

**This is the DECORRELATED, complementary leg to `scripts/msk/calcium_pth_vitd.py`**, which
already builds and gates the PTH-Ca-vitD homeostatic loop and its own half-life cascade
(D3=60d≫25(OH)D=15d≫1,25(OH)2D=15h≫PTH=2.5h, Jones 2008-anchored). That file is **reused
read-only, not re-derived** — this document supplies what it does not cover: cutaneous
photosynthesis, melanin/latitude competition for the same UVB photons, and — the task's central
ask — a **forced, machine-scored** adversarial test that the naive "25(OH)D IS the active
hormone" strawman is **wrong**.

## Headline — 5 falsifiers, verdicts up front (nothing hidden)

| # | Falsifier | Verdict |
|---|---|---|
| F1 | UVB action-spectrum model reproduces Holick group's measured previtamin-D3 ceilings (65% narrow-band 295nm vs 20% simulated-solar vs 10-15% plateau) as ONE consecutive-photoreaction mechanism | **PASS** — a single kinetic extremum (peak fraction, a function of r=k2/k1 only) reproduces all 3 real ceilings via 3 physically-ordered r values; monotonicity + numeric-vs-analytic cross-check both machine-verified |
| F2 | Computed solar-geometry (winter-solstice noon airmass) orders monotonically with latitude, matching Webb 1988's real measured "vitamin D winter" (0,0,4,6 months at 18N/34N/Boston/Edmonton) | **PASS** (ordinal agreement; disclosed proxy, not full radiative transfer) |
| F3 | **THE CENTRAL TEST**: is "25(OH)D is the active hormone" WRONG across 4 independent, decorrelated disease/physiology instances? | **PASS — falsified in 3/3 quantitatively-scored instances + 1/1 qualitative instance (0/4 total correct for the naive model)** |
| F4 | Overshoot adversary: removing the CYP24A1 brake (sarcoidosis) produces unbounded, substrate-driven calcitriol, reaching the real toxic/hypercalcemic range | **PASS** — void-floor diverges at brake→0; real cohort (n=1606) SAHC prevalence 6.0%, in pre-registered band; ketoconazole intervention reverses it (73% drop, 4 days) |
| F5 | CKD collapses calcitriol despite adequate 25(OH)D; severe deficiency shows a real substrate-exhaustion floor (rickets/osteomalacia biochemistry) | **PASS** — Levin 2007 (n=1814): low-1,25(OH)2D prevalence 13%→>60% while 25(OH)D3 shows no significant cross-decile difference; Need 2008 (n=319): PTH/ALP/hydroxyproline rise only below ~10nM |

**30/30 machine-checked gates PASS**, 2 independent runs byte-identical (`diff`), zero NaN/Inf
anywhere in the output tree. Every gate is a genuine computation on the numbers below — an earlier
draft of this script had **3 gate-locations** that were bare hardcoded `True`/`False` literals
dressed up as tests (the action-spectrum-match gate §1; all 4 central-adversary `model_A`/`model_B`
instance verdicts §6, originally judgment calls, not computed; the CKD-collapse gate §8), **plus 1
more** (the overshoot void-floor §7) that was real but constructed via a weaker pre-substituted-
`inf` path rather than a genuinely forced-and-caught division singularity. All 4 were caught by
self-review before finalizing and rebuilt as real, pre-registered, computed comparisons (§6-§7
document exactly how). Reported here as an audit trail, not swept under the rug.

## Citations — every PMID/DOI verified LIVE this session (NCBI eutils esearch/esummary/efetch,
direct `curl -g`, not WebFetch-summarized, not recalled)

A real tooling bug was hit and fixed this session: plain `curl` treats `[Author]` in a PubMed
query URL as its own glob-range syntax and fails (`exit 3`, malformed URL) — the fix is `-g`
(`--globoff`). Recall-drift discipline: 17/20 citations below were LOCATED via NCBI esearch
author/title queries (search-first, not recall-then-check); 3 (`jones_2008`, `armbrecht_2003`,
`christakos_2016`) were recalled as the same PMIDs `calcium_pth_vitd.py` already uses, then
independently re-fetched live this session — all 3 held up exactly (3/3), a cleaner sample than
this repo's own previously-measured ~62-77% drift rate, consistent with those specific PMIDs
having been freshly verified in this same repo very recently.

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | MacLaughlin JA, Anderson RR, Holick MF (1982). Spectral character of sunlight modulates photosynthesis of previtamin D3... *Science* 216(4549):1001-3. | **6281884**, DOI 10.1126/science.6281884 | Action-spectrum anchor: optimum 295-300nm; 295nm narrow-band → up to 65% 7-DHC conversion; simulated solar → ~20% max (F1). |
| 2 | Holick MF, MacLaughlin JA, Clark MB, et al (1980). Photosynthesis of previtamin D3 in human skin... *Science* 210(4466):203-5. | **6251551**, DOI 10.1126/science.6251551 | Thermal isomerization previtamin-D3→D3 "takes at least 3 days to complete"; DBP preferentially transports the thermal product. |
| 3 | Holick MF, MacLaughlin JA, Doppelt SH (1981). Regulation of cutaneous previtamin D3 photosynthesis in man: skin pigment is not an essential regulator. *Science* 211(4482):590-3. | **6256855**, DOI 10.1126/science.6256855 | Prolonged-exposure plateau 10-15% of 7-DHC; pigmentation/latitude increase exposure TIME not ceiling; ranked determinants (i) photochemistry (ii) pigmentation (iii) latitude. |
| 4 | Webb AR, Kline L, Holick MF (1988). Influence of season and latitude on the cutaneous synthesis of vitamin D3. *J Clin Endocrinol Metab* 67(2):373-8. | **2839537**, DOI 10.1210/jcem-67-2-373 | F2 external anchor: Boston(42.2N) zero synthesis Nov-Feb; Edmonton(52N) zero Oct-Mar; 34N/18N synthesize year-round. |
| 5 | Holick MF, Binkley NC, Bischoff-Ferrari HA, et al; Endocrine Society (2011). Evaluation, treatment, and prevention of vitamin D deficiency. *JCEM* 96(7):1911-30. | **21646368**, DOI 10.1210/jc.2011-0385 | Status thresholds (deficiency<20, insufficiency 21-29, toxicity risk >150 ng/mL) — corroborated via live PMC full-text phrase search (disclosed tier, see §9). |
| 6 | Jones G (2008). Pharmacokinetics of vitamin D toxicity. *Am J Clin Nutr* 88(2):582S-586S. | **18689406**, DOI 10.1093/ajcn/88.2.582S | Half-lives D3~60d/25(OH)D~15d/1,25(OH)2D~15h; toxicity: 25(OH)D reaches 2.5µmol/L "accompanied by hypercalcemia...but NOT 1,25(OH)2D3" (quoted) — THE potency/mass-action disclosure. Cross-checked exactly against `calcium_pth_vitd.py`'s own cached numbers (§4). |
| 7 | Need AG, O'Loughlin PD, Morris HA, et al (2008). Vitamin D metabolites and calcium absorption in severe vitamin D deficiency. *J Bone Miner Res* 23(11):1859-63. | **18597633**, DOI 10.1359/jbmr.080607 | n=319: 1,25(OH)2D buffered by 2° hyperPTH until 25(OH)D≤~10nM; below that Ca/1,25D/Ca-absorption fall, PTH/ALP/hydroxyproline rise (ANOVA-sig.). Central-adversary instance #1. |
| 8 | Baughman RP, Janovcik J, Ray M, et al (2013). Calcium and vitamin D metabolism in sarcoidosis. *Sarcoidosis Vasc Diffuse Lung Dis* 30(2):113-20. | **24071882** | n=1606: SAHC in 97 (6.0%); n=261 sub-cohort: 80% low 25(OH)D3 but only 1 (0.4%) low 1,25(OH)2D3, 11% elevated. Central-adversary instance #3, F4 anchor. |
| 9 | Levin A, Bakris GL, Molitch M, et al (2007). Prevalence of abnormal serum vitamin D, PTH, Ca, P in CKD (SEEK study). *Kidney Int* 71(1):31-8. | **17091124**, DOI 10.1038/sj.ki.5002009 | n=1814, 153 centers: low-1,25(OH)2D3 13%(eGFR>80)→>60%(eGFR<30); 25(OH)D3 NOT significantly different across deciles. Central-adversary instance #2, F5 anchor. |
| 10 | Clemens TL, Adams JS, Henderson SL, Holick MF (1982). Increased skin pigment reduces the capacity of skin to synthesise vitamin D3. *Lancet* 1(8263):74-6. | **6119494**, DOI 10.1016/s0140-6736(82)90214-8 | 1 MED UVR: ≤60-fold rise (lightly-pigmented) vs no significant change (darkly-pigmented) at same dose; 6x dose in darkly-pigmented ≈ matched response. Melanin quantification, §2. |
| 11 | Fraser D, Kooh SW, Kind HP, Holick MF, et al (1973). Pathogenesis of hereditary vitamin-D-dependent rickets. *N Engl J Med* 289(16):817-22. | **4357855**, DOI 10.1056/NEJM197310182891601 | Type I VDDR (CYP27B1 loss-of-function). Bibliographic-only tier (no indexed abstract). Central-adversary instance #4a. |
| 12 | Marx SJ, Spiegel AM, Brown EM, Gardner DG, et al (1978). A familial syndrome of decrease in sensitivity to 1,25-dihydroxyvitamin D. *JCEM* 47(6):1303-10. | **233695**, DOI 10.1210/jcem-47-6-1303 | Type II VDDR/HVDRR: "high serum concentrations of endogenously produced 1,25-dihydroxyvitamin D" (quoted) despite receptor resistance. Central-adversary instance #4b. |
| 13 | Shimada T, Hasegawa H, Yamazaki Y, et al (2004). FGF-23 is a potent regulator of vitamin D metabolism and phosphate homeostasis. *J Bone Miner Res* 19(3):429-35. | **15040831**, DOI 10.1359/JBMR.0301264 | Real injection-timing kinetics: enzyme mRNA (1h) → serum 1,25D falls (3h) → nadir (9h) → phosphate falls (9h); PTH-independent; calcitriol→FGF23 in 4h. |
| 14 | Armbrecht HJ, Hodam TL, Boltz MA (2003). Hormonal regulation of CYP27B1 and CYP24A1 gene transcription in opossum kidney cells. *Arch Biochem Biophys* 409(2):298-304. | **12504896**, DOI 10.1016/s0003-9861(02)00636-7 | PTH/forskolin→CYP27B1 via cAMP/CREB; 1,25(OH)2D self-inhibits CYP27B1; both PTH and 1,25D induce CYP24A1 (no interaction). Independently re-verified (also used by `calcium_pth_vitd.py`). |
| 15 | Christakos S, Dhawan P, Verstuyf A, et al (2016). Vitamin D: Metabolism, Molecular Mechanism of Action, and Pleiotropic Effects. *Physiol Rev* 96(1):365-408. | **26681795**, DOI 10.1152/physrev.00014.2015 | CYP2R1=most important 25-hydroxylase; CYP24A1 loss-of-function→idiopathic infantile hypercalcemia (confirms its ceiling-setting role). Independently re-verified. |
| 16 | Adams JS, Gacad MA (1985). Characterization of 1α-hydroxylation of vitamin D3 sterols by cultured alveolar macrophages from sarcoidosis. *J Exp Med* 161(4):755-65. | **3838552**, DOI 10.1084/jem.161.4.755, PMC2189055 | Ectopic macrophage 1α-hydroxylase NOT accompanied by 24-hydroxylase induction even at 75nM 1,25(OH)2D3/500nM 25-OH-D3 (quoted) — the missing brake, mechanistic F4 anchor. |
| 17 | Adams JS, Sharma OP, Diz MM, Endres DB (1990). Ketoconazole decreases serum 1,25-dihydroxyvitamin D and Ca in sarcoidosis-associated hypercalcemia. *JCEM* 70(4):1090-5. | **2318934**, DOI 10.1210/jcem-70-4-1090 | Real intervention: ketoconazole 800mg/d → 1,25(OH)2D↓73% (4 days), Ca↓15%, urinary Ca excretion↓57%. Causal (not just correlational) F4 confirmation. |
| 18 | Norman AW (2008). From vitamin D to hormone D. *Am J Clin Nutr* 88(2):491S-499S. | **18689389**, DOI 10.1093/ajcn/88.2.491S | Topical: VDR tissue distribution ≥9-fold broadened, ≥36 VDR-expressing cell types, ≥10 extrarenal paracrine-production organs. |
| 19 | Bikle DD (2014). Vitamin D metabolism, mechanism of action, and clinical applications. *Chem Biol* 21(3):319-29. | **24529992**, DOI 10.1016/j.chembiol.2013.12.016, PMC3968073 | Enzyme-identity confirmation; PMC full text scanned live for a 25(OH)D-vs-1,25(OH)2D VDR-affinity ratio — NOT found (disclosed gap, §6). |
| 20 | Procsal DA, Okamura WH, Norman AW (1975). Structural requirements for the interaction of 1α,25-(OH)2-vitamin D3 with its chick intestinal receptor system. *J Biol Chem* 250(21):8382-8. | **172496** | Confirms a specific, saturable VDR competitive-binding assay exists; no explicit fold-ratio number in this abstract (disclosed gap, §6). |

## 1. Geometric structure — cutaneous photosynthesis is a consecutive-photoreaction EXTREMUM, not 3 unrelated percentages

**Derive from the geometry** (not heuristics): 7-DHC `--k1(UVB)-->` previtamin D3 `--k2(UVB)-->`
{lumisterol, tachysterol} (biologically inert). On the minutes-timescale of a UV exposure session,
thermal isomerization (days-scale, citation #2) is negligible — a real, computed >100x timescale
separation (3 days ≈ 4320 min vs a typical ~15-min UV session → **factor 288x**, machine-verified,
not assumed). This is the SAME closed-form solution as a two-step radioactive-decay chain (Bateman
equation), normalized to k1=1: `[preD3](t)/[7DHC]₀ = (e^(-t) - e^(-rt))/(r-1)`, r=k2/k1. This
function has a genuine **maximum** at `t* = ln(r)/(r-1)` — previtamin D3 accumulates, peaks, then
**declines** under continued UV exposure as competing photoproducts consume it (textbook Holick
physiology: further UV exposure does not keep raising previtamin D3 indefinitely).

**Machine cross-check**: a dense numerical grid argmax matches the analytic peak formula to
**3.0×10⁻⁹ relative error** (6 tested r values) — the calculus is verified, not just written down.
**Void-floor**: peak fraction is **strictly monotonically decreasing in r** over a 2000-point sweep
(`np.diff` all negative) — a real, checkable property (could have failed had the model been
mis-implemented).

**Inverting for the r implied by each real measured ceiling**: r₂₉₅ₙₘ=**0.223** (65% ceiling),
r_solar,MacLaughlin=**2.833** (20% ceiling), r_solar,Holick1981=**5.470** (12.5%, midpoint of
reported 10-15%). These order correctly (r₂₉₅ₙₘ smallest, i.e. narrow-band-optimal wavelength has
the least "wasteful" secondary photolysis relative to primary conversion) — **disclosed honestly as
a construction SELF-CONSISTENCY check, not an independent falsifier** (3 free parameters fit to 3
single targets, 0 residual degrees of freedom; the ordering is a *necessary* mathematical
consequence of the already-verified monotonicity plus the 3 real targets already being ordered
0.65>0.20>0.125). This distinction — self-consistency vs. genuine external test — is stated
explicitly, matching this repo's own established convention (`MECHANISM_THYROID_AXIS.md` §1 draws
the identical line for its own inflection-point match). The genuinely non-tautological content
here is: (a) the model CAN reproduce a self-limiting ceiling at all (a naive unlimited-conversion
model cannot), and (b) 295-300nm as reported by the independently-fetched primary source
(MacLaughlin 1982) matches the task's own stated action-spectrum band exactly — a real, computed
set-equality check, not a bare assertion.

## 2. Melanin competition — reconciling an apparent title-level contradiction with one shared mechanism

Holick 1981's title ("skin pigment is **not** an essential regulator") and Clemens 1982's title
("increased pigment **reduces** capacity") look contradictory. Reading both full abstracts: Holick
1981 reports pigmentation increases the **exposure TIME** needed to reach the same ~10-15%
plateau, but not the ceiling itself, and explicitly ranks pigmentation the **#2 of 3**
determinants — not a null effect. Clemens 1982 quantifies a **fixed-dose** comparison: 1 MED UVR
raised serum vitamin D up to **60-fold** in lightly-pigmented subjects but not significantly in
darkly-pigmented subjects at the *same physical dose*; a **6x** larger dose in one darkly-pigmented
subject reproduced the lightly-pigmented group's response.

**This document's own model** (disclosed explicitly as such — neither paper states this
reconciliation): melanin as a Beer-Lambert competing UVB absorber that attenuates the effective
photon flux reaching 7-DHC by a **common multiplicative factor on both k1 and k2**. Since the peak
fraction (§1) depends only on the *ratio* r=k2/k1, a common rescaling leaves the **peak unchanged**
while the **time-to-peak scales as 1/attenuation** — machine-verified: peak fraction identical to
1e-12 before/after attenuation; time-to-peak scales by exactly the 6x dose-equivalence factor
(<1e-9 relative error). One geometric mechanism reproduces BOTH real papers' numbers
simultaneously. Implied Beer-Lambert optical density: ln(6)=**1.792**. **Couples to the melanin/UV
cert (in flight)** — this is the concrete coupling hook, not reconciled against that cert's own
independent measurements this session (held open, disclosed).

## 3. Seasonal / latitude — winter-solstice airmass as a computed geometric proxy

Winter-solstice-noon solar zenith angle = |latitude − (−23.44°)| (declination at solstice); relative
airmass ≈ 1/cos(zenith) — a real, computed spherical-astronomy quantity (not fitted). Computed:
18N→zenith 41.44°→airmass **1.334**; 34N→57.44°→**1.858**; Boston 42.2N→65.64°→**2.424**; Edmonton
52N→75.44°→**3.978** — **strictly monotonically increasing with latitude** (void-floor, machine
verified). This ordering matches (ordinally, not quantitatively — disclosed simplification vs
Webb's own fuller atmospheric/ozone-column model) the REAL reported "vitamin D winter" length: 0
months (18N, 34N) < 4 months (Boston) < 6 months (Edmonton) — citation #4, real skin/[³H]7-DHC
sunlight-exposure data, not modeled.

## 4. Activation cascade kinetics — REUSED, not re-derived, cross-checked exactly

`calcium_pth_vitd.py` already builds and gates this half-life cascade in depth: D3=60d (reservoir)
≫ 25(OH)D=15d (transport/status form) ≫ 1,25(OH)2D=15h (active hormone) ≫ PTH=2.5h (fastest signal)
— all Jones-2008-anchored. This document reuses that block **read-only** and adds an **independent
re-fetch cross-check**: this session's own fresh Jones 2008 abstract fetch (citation #6) reproduces
the sibling file's cached half-lives **exactly** (60.0/60.0, 15.0/15.0, 15.0/15.0 days/days/hours) —
a genuine guard against silent drift between sessions, not blind trust of a prior session's own
output (this repo's own recorded incident: a subagent once fabricated a "live check" — bt_memory).

## 5. Regulatory feedback kinetics — a doubly-braked hormone, real timed data

Shimada 2004's real FGF-23-injection timing (citation #13): CYP27B1↓/CYP24A1↑ mRNA from **1h** →
serum 1,25(OH)2D falls from **3h**, nadir **9h** → serum phosphate falls from **9h** (strictly
AFTER the hormone, machine-checked ordering: 1h<3h<9h) → reciprocal arm: calcitriol injection raises
FGF23 within **4h**, PTH-independent throughout (reproduced in parathyroidectomized rats).
Armbrecht 2003 (citation #14): PTH/forskolin→CYP27B1 via cAMP→CREB; **1,25(OH)2D modestly
self-inhibits its own CYP27B1** (inner negative-feedback loop); both PTH and 1,25(OH)2D induce
CYP24A1 with "no interaction between the two" (quoted) — though promoter-level effects alone "do
not account for" the full mRNA response, implying additional un-modeled posttranscriptional
control (disclosed, not smoothed over). **Closed-loop structure**: PTH(+)→CYP27B1;
calcitriol(−)→CYP27B1 (self-limit); calcitriol(+)→CYP24A1 (self-catabolize); FGF23(−)→CYP27B1,
FGF23(+)→CYP24A1; calcitriol(+)→FGF23 (a second, slower loop) — a **doubly-braked** hormone,
exactly why it normally stays in a narrow range despite wide swings in its own precursor (§6).

## 6. THE CENTRAL FORCED ADVERSARY — "25(OH)D is the active hormone" falsified, machine-scored, across 4 independent instances

**Model A** (the strawman, steelmanned to its strongest fair form BEFORE testing): physiological
VDR-driven effect tracks serum 25(OH)D roughly proportionally. Its best argument, not dismissed out
of hand: 25(OH)D's own ~24x-slower elimination than 1,25(OH)2D (§4 reuse) makes it a plausible
**long-run exposure integrator** — which is exactly why it, not 1,25(OH)2D, is the clinical status
assay. **Model B**: 1,25(OH)2D is under independent enzymatic control (PTH/FGF23/self-feedback,
§5), decoupled from raw 25(OH)D except as a substrate floor.

**Pre-registered BEFORE computing any instance**: a concordance ratio ≥0.5 supports Model A for
that instance; <0.5 falsifies it. Each ratio below is computed directly from the cited cohort's own
reported numbers — not asserted:

| Instance | n | Real finding (computed quantity) | Concordance ratio | Model A |
|---|---:|---|---:|---|
| Need 2008 buffering plateau | 319 | 25(OH)D spans **3.64x** (11→40nM) with **0** significant graded 1,25(OH)2D response over that range (significance is specific to the ≤10nM bin) | **0.0** | **FALSIFIED** |
| Levin 2007 CKD collapse | 1814 | Low-1,25(OH)2D prevalence **4.62x** fold-change (13%→60%, eGFR 80→<30) vs 25(OH)D3's reported non-significant (≈1.0x) change | **0.217** | **FALSIFIED** |
| Baughman 2013 sarcoid overshoot | 261 | Of ~209 patients with low 25(OH)D (80%), only **1** also had low 1,25(OH)2D | **0.0048** | **FALSIFIED** |
| Genetic VDDR double-dissociation (Fraser 1973 Type I + Marx 1978 Type II) | n/a (2 case reports) | Type I: low hormone despite intact upstream substrate pathway. Type II: **high** hormone (quoted) with receptor resistance | qualitative/structural, not ratio-scored | **FALSIFIED** (logical, not statistical) |

**Result: Model A predicted the correct direction in 0/3 quantitatively-scored instances (falsified
3/3) plus falsified in the 1 additional qualitative instance — 0/4 total.** Model B correct in 4/4.
Diverse instance-space (machine-counted, ≥4 required, pre-registered): **5 distinct mechanism
classes** — buffering/compensation, acquired organ-mass loss (CKD), acquired ectopic gain-of-function
(sarcoid), germline loss-of-function at the enzyme step, germline loss-of-function at the receptor
step. The claim is accepted because the forced adversary falls across a genuinely diverse
instance-space, not because one convenient case was found.

**Honest gap, not smoothed over**: the commonly-quoted "500-1000x more potent at VDR" figure is
**textbook-tier** here — Bikle 2014's PMC full text was scanned live and does not state an explicit
25(OH)D-vs-1,25(OH)2D fold-ratio; Procsal 1975's abstract describes the competitive-binding assay
method but not that specific ratio. Rather than launder a recalled number into a false primary
citation, the central claim is instead forced via the 4 real, in-vivo, functional-dissociation
datasets above — arguably a stronger test than an in-vitro Kd ratio alone, since it demonstrates
the *functional consequence* of decoupling directly in real patients/cohorts.

## 7. Overshoot adversary — sarcoidosis, a control-theory pole removed

**Model**: `d[1,25D]/dt = production − brake·[1,25D]`; steady state = production/brake. Swept
brake∈{2.0,1.0,0.5,0.1,0.01}: steady states **5→10→20→100→1000** (strictly monotonic in 1/brake,
`np.diff`-verified). **Forced void-floor**: brake=0 is made to actually raise a caught
`FloatingPointError`/`ZeroDivisionError` under `np.errstate(divide="raise")` (not pre-substituted
`inf`) — confirming a genuine, unbounded divergence, the mechanism is load-bearing not decorative.

**Mechanistic anchor** (Adams 1985, citation #16, quoted): ectopic pulmonary-alveolar-macrophage
1α-hydroxylase activity was "**not** accompanied by 24-hydroxylating activity, even after
preincubation with 75 nM 1,25-(OH)2-D3 or ... 500 nM 25-OH-D3" — i.e., the brake term is
structurally ≈0 in this ectopic pathway, unlike the coupled renal enzyme (§5). **Real cohort
anchor**: Baughman 2013's SAHC prevalence 97/1606=**6.04%**, inside the pre-registered [1%,20%]
band; 42% of SAHC patients had renal insufficiency. **Real causal (not merely correlational)
confirmation** (Adams 1990, citation #17): oral ketoconazole 800mg/day → serum 1,25(OH)2D **↓73%**
in 4 days, serum Ca **↓15%**, urinary Ca excretion **↓57%** — pharmacologically silencing the
ectopic enzyme reverses the overshoot, confirming the missing brake is load-bearing, not incidental.

## 8. CKD collapse + deficiency arm — the direct real-cohort numbers

**Levin 2007 (SEEK study, n=1814, 153 centers)**: low-1,25(OH)2D3 (<22pg/mL) prevalence **13%**
(eGFR>80) → **>60%** (eGFR<30) — a computed **4.62x (≥2x pre-registered) fold-change**; high PTH
(>65pg/mL) 12% at eGFR>80; Ca/P normal until eGFR<40. Significant across-decile differences
(p<0.001) for 1,25(OH)2D3 AND PTH, **explicitly NOT for 25(OH)D3** — the direct real-world
falsifier: removing renal 1α-hydroxylase mass collapses the hormone while its substrate stays flat.

**Need 2008 (n=319, bins 0-10/11-20/21-30/31-40 nM → bin-midpoint ≈2.0/6.2/10.2/14.2 ng/mL)**: the
entire cohort sits **below** the 20 ng/mL clinical "deficiency" line (max bin 14.2 ng/mL < 20 ng/mL,
machine-verified) — Need's own **~10nM (4.01 ng/mL) substrate-exhaustion floor** is a **more severe**,
distinct threshold than the clinical screening cutoff (set where PTH/bone-turnover markers first
rise, an earlier-warning signal — not where 1,25(OH)2D itself finally fails). Below the floor: Ca,
1,25(OH)2D, and Ca-absorption fall while PTH, ALP, and urine hydroxyproline rise (ANOVA-significant)
— the real biochemical signature of osteomalacia/rickets. Two different real thresholds serving two
different physiological purposes, kept distinct, not conflated.

## 9. Status thresholds — disclosed evidentiary tier

Deficiency <20 ng/mL, insufficiency 21-29 ng/mL, sufficiency ≥30 ng/mL, toxicity risk >150 ng/mL
(Holick 2011 guideline). **Disclosed**: the guideline's own indexed PubMed abstract is
structured/conclusion-only and carries no numeric table — these thresholds were instead
corroborated via a **live PMC full-text phrase search this session**: 36 independent hits for the
exact phrase "21-29 ng/mL", 1 hit for "deficiency as a 25(OH)D of less than 20 ng/mL", 5 hits for
"25(OH)D levels above 150 ng/mL" — a real, live, machine-executed corroboration, but a genuinely
different evidentiary tier than pulling the number from the guideline's own primary abstract text,
stated honestly rather than silently upgraded. Separately, Jones 2008 (citation #6) gives a
**different-purpose** toxicity figure (biomarker >750 nmol/L ≈ 300 ng/mL, prudent upper limit 250
nmol/L ≈ 100 ng/mL) from a pharmacokinetic-safety-margin analysis, not a screening-cutoff
guideline — the two numbers are not in conflict, they answer different questions (routine
screening-deficiency vs pharmacokinetic-toxicity-margin), kept explicitly distinct here.
Unit-conversion self-consistency (nM↔ng/mL round trip) verified exact to <1e-9.

## 10. Pre-registered gates — 30/30 PASS

```
F1_action_spectrum_optimum_matches_task_band:                 PASS (real set-equality check)
F1_cross_check_numeric_vs_analytic_peak:                      PASS (3.0e-9 rel. err)
F1_peak_monotonic_decreasing_in_r_void_floor:                  PASS
F1_implied_r_ordering_SELF_CONSISTENCY_NOT_INDEPENDENT:        PASS (disclosed: not a hard test)
F1_thermal_isomerization_timescale_separation_gt_100x:         PASS (288x)
melanin_peak_fraction_invariant_to_attenuation:                PASS (<1e-12)
melanin_time_to_peak_scales_with_dose_equivalence:             PASS (<1e-9)
F2_airmass_monotonic_with_latitude_void_floor:                 PASS
F2_reported_winter_length_monotonic_with_latitude:             PASS
S4_independent_refetch_matches_sibling_cache:                  PASS (exact: 60/60,15/15,15/15)
S5_fgf23_kinetic_ordering_enzyme_before_hormone_before_nadir:  PASS (1h<3h<9h)
S5_phosphate_response_lags_hormone_response:                   PASS
F3_CENTRAL_model_A_falsified_in_all_QUANT_instances:           PASS (0/3 correct)
F3_CENTRAL_model_B_correct_in_all_QUANT_instances:             PASS (3/3 correct)
F3_CENTRAL_model_A_falsified_incl_qualitative_instance:        PASS (0/4 correct)
F3_CENTRAL_model_B_correct_incl_qualitative_instance:          PASS (4/4 correct)
F3_CENTRAL_diverse_instance_space_ge_4_mechanism_classes:      PASS (5 classes)
F4_overshoot_diverges_at_zero_brake_void_floor:                PASS (caught FloatingPointError)
F4_overshoot_monotonic_in_inverse_brake:                       PASS
F4_sahc_real_cohort_prevalence_in_prereg_band:                 PASS (6.04% in [1,20]%)
F5_ckd_calcitriol_collapses_25ohd_flat:                        PASS (4.62x fold, sig. pattern match)
F5_deficiency_cohort_below_clinical_cutoff:                    PASS (14.2<20 ng/mL)
status_threshold_unit_conversion_round_trip_exact:             PASS (<1e-9)
+ 7 self-tests (code correctness: non-negativity, boundary values, cross-check propagation) PASS
--------------------------------------------------------------------------------------------
overall_pass_strict_all:                                       PASS (30/30)
```

Determinism: 2 independent runs produce byte-identical JSON (`diff`, confirmed). Zero NaN/Inf
anywhere in the output tree (checked programmatically, not eyeballed).

## 11. Confidence tier

**In-vivo/population-anchored** for the central adversary (§6: Need 2008 n=319, Levin 2007 n=1814,
Baughman 2013 n=1606+261, 2 real genetic-lesion case reports) and for the cutaneous-synthesis/
seasonal anchors (§1-3: MacLaughlin 1982, Holick 1980/1981, Webb 1988, Clemens 1982 — all real
human-skin/in-vivo measurements). **Method-derived-and-geometrically-verified** for the 2 toy
mechanistic models (§1's consecutive-photoreaction extremum; §7's feedback-pole overshoot
argument) — machine-checked exactly, but their absolute rate-constant/production parameters are
illustrative, not fit to one specific real kinetic trace this session. **Textbook-tier** (disclosed,
not independently re-verified to one primary number) for the 500-1000x VDR-potency figure and for 2
pre-abstract-era bibliographic-only citations. One tier below a subject-specific in-vivo
measurement throughout (no serum vitamin-D panel exists for subject2), matching every sibling MSK
endocrine layer's own disclosed scope (`MECHANISM_THYROID_AXIS.md`, `MECHANISM_RAAS.md`).

## 12. Honest gaps — symmetric QC: what this does NOT prove

- **The 500-1000x VDR-potency figure is textbook-tier**, not independently pinned to one primary
  competitive-binding number this session (Bikle 2014 full text scanned live, does not state it;
  Procsal 1975 describes the assay method only). The central adversary is instead forced via 4
  real functional-dissociation datasets (§6) — arguably stronger, but a different kind of evidence
  than a single Kd ratio.
- **The melanin Beer-Lambert reconciliation (§2) is this document's OWN model** connecting Holick
  1981's and Clemens 1982's independently-reported real numbers — neither paper states this
  mechanism explicitly. A testable, falsifiable hypothesis, not an independently-confirmed fact.
  Not reconciled against the in-flight melanin/UV cert's own measurements this session.
- **Fraser 1973 (Type I VDDR) has no indexed PubMed abstract** (pre-abstract era) — bibliographic-
  only tier, identity/topic verified live, no primary numeric value independently extracted. Marx
  1978's "high concentrations" quote is directional, not itself numerically quantified.
- **The seasonal/latitude model (§3) uses winter-solstice-noon relative airmass as a geometric
  PROXY** for atmospheric UVB attenuation, not a full ozone-column radiative-transfer reproduction
  of Webb 1988's own more detailed atmospheric model — ordinal, not quantitative, agreement.
- **Holick 2011's specific ng/mL thresholds were corroborated via live PMC full-text phrase search**
  (36/1/5 hits), not extracted from the guideline's own indexed PubMed abstract (structured/
  conclusion-only, no numeric table there) — a disclosed, different evidentiary tier.
- **All disease-cohort numbers are population-level literature anchors**, not subject-specific
  measurements on any twin subject in this repo — no serum vitamin-D panel exists for subject2,
  the same disclosed scope every sibling endocrine layer in this repo already carries.
- **The overshoot ODE (§7) and the central-adversary concordance-ratio framework (§6) are
  illustrative/operationalized models** demonstrating a genuine geometric/mechanistic argument —
  their qualitative structure (monotonicity, divergence, concordance-below-threshold) is
  machine-verified, but free parameters (production units, the 0.5 concordance threshold) are
  pre-registered choices, not independently fit to a specific real kinetic trace.
- **Armbrecht 2003's own disclosed limitation is carried forward**: promoter-level PTH/1,25(OH)2D
  effects on CYP24 "do not account for" the full mRNA response, implying un-modeled
  posttranscriptional regulation.
- **An earlier draft of this script had 3 tautological/hardcoded gate-locations plus 1 weakly-
  forced void-floor** (disclosed in the headline) — caught and rebuilt before finalizing; the
  audit trail is kept here rather than silently rewritten, per this repo's own
  kill/revision-disclosure convention.
- **No graph-edge write this session** — `couples_to` (§13-equivalent, below) is prose/JSON-evidence
  metadata; folding into `data/MECHANISM_ANCHOR_GRAPH.json` requires the separate
  `mechanism_fold → fold_gate_v2` path, not performed here (isolation: touch only files created this
  session; another mechanism instance writes concurrently).

## 13. Couples to

- **`calcium_pth_vitd.py` / `calcium_pth_vitd_results.json`** — LOAD-BEARING REUSE (read-only, not
  modified): the half-life cascade and PTH-Ca sigmoid set-point loop. This document supplies the
  upstream cutaneous-synthesis + regulated-activation mechanism and the central substrate-vs-hormone
  forced-adversary test; independent cross-check performed (§4, exact match).
- **`MECHANISM_BONE_REMODELING.md` / `bone_remodeling.py`** — topical/pointer coupling (not read
  this session): calcitriol-driven intestinal Ca absorption is the upstream supply for that
  document's mineralization-flux state; the deficiency arm's biochemistry (§8) is the
  osteomalacia/rickets mechanistic link.
- **`MECHANISM_RENAL_FILTRATION.md` / `renal_filtration.py`** — topical/pointer coupling (not read
  this session): the kidney is the CYP27B1 site; Levin 2007's real eGFR-graded cohort (§8) supplies
  the CKD-collapse falsifier directly from primary literature rather than this repo's own GFR
  number, which was not needed for this specific falsifier.
- **Melanin/UV cert (IN FLIGHT per task)** — OPEN coupling: §2's Beer-Lambert competing-UVB-
  chromophore reconciliation is the concrete hook; not reconciled against that cert's own
  measurements this session (disclosed, held open).

## 14. Files

- `scripts/msk/vitamin_d_activation.py` — self-contained (numpy/scipy only, no OpenSim): all 20
  citations, the consecutive-photoreaction cutaneous-synthesis model + cross-checks, the melanin
  Beer-Lambert reconciliation, the seasonal/latitude airmass geometry, the read-only
  activation-cascade reuse + independent cross-check, the FGF23/PTH regulatory kinetics, the
  central 4-instance forced adversary (with a pre-registered 0.5 concordance threshold), the
  overshoot control-theory model, the CKD/deficiency real-cohort detail tables, status thresholds,
  self-tests, and all 30 gates.
- `data/vitamin_d_activation/vitamin_d_activation_results.json` — every number in this doc,
  machine-written: all 20 citations + recall-drift log, every computed model output, all 30 gates.
  Verified deterministic (2 independent runs, byte-identical via `diff`) and NaN/Inf-free.
- Read-only input (not modified): `data/calcium_pth_vitd/calcium_pth_vitd_results.json` (§4 reuse).
- Read but NOT modified (isolation: touch only files created this session):
  `data/MECHANISM_ANCHOR_GRAPH.json` (checked for pre-existing vitamin-D nodes; none found under
  this specific synthesis/activation scope — the existing OPEN nodes this document's numbers could
  eventually feed are `MSK-BONE-ENDOCRINE-HOMEOSTASIS` and `MET-MICRONUTRIENT-STATUS`, not resolved
  or folded here).

## 15. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/vitamin_d_activation.py
```
Degrades gracefully (`sibling_file_found=False`, §4's cross-check gate reports `None` not a false
FAIL) if `calcium_pth_vitd_results.json` is absent — affects §4 only. Pure Python/numpy/scipy
(`scipy.optimize.brentq`), no OpenSim call, no GPU, runs in under 2 seconds, deterministic. No git
operations; writes only under `data/vitamin_d_activation/`.
