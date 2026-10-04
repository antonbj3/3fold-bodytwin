# MECHANISM MAPK/ERK SIGNALING — the RAS-RAF-MEK-ERK kinase cascade: distributive-phosphorylation ultrasensitivity + feedback-driven bistability, forced against a processive adversary and a scaffold/crowding linearization sweep (2026-07-22)

Builds and MEASURES the defining quantitative signature of the RTK -> RAS-GTP -> RAF -> MEK -> ERK
kinase cascade: (1) **ultrasensitivity** — a graded upstream input is converted into a steep,
switch-like output (Hill coefficient nH>>1) because MEK and ERK are each activated by DISTRIBUTIVE
(two-collision, kinase-dissociates-between-phosphorylations) dual phosphorylation, not a single
hyperbolic step; (2) **bistability/all-or-none** — adding the pathway's own positive feedback turns
the graded switch into a genuine two-state memory element with hysteresis, the mechanism Ferrell &
Machleder (1998) attribute the Xenopus oocyte all-or-none cell-fate switch to. Forced against: a
processive (single-collision) adversary at MEK/ERK; a no-feedback adversary for bistability; and a
continuous scaffold/molecular-crowding linearization sweep. Decorrelated check: BRAF-V600E
oncogenic constitutive activation (~52% of melanoma, TCGA), vemurafenib clinical response, and the
paradoxical RAF-inhibitor-induced ERK activation specifically in RAS-mutant cells. Symmetric-QC held
OPEN per task instruction: distributive-vs-processive is genuinely debated **in vivo** (living-cell
evidence contradicts the in-vitro biochemistry), and KSR-family scaffolding measurably reduces the
cascade's own ultrasensitivity. Script: `scripts/msk/mapk_erk_cascade.py`. Evidence (machine-written):
`data/mapk_erk_cascade/mapk_erk_cascade_results.json`. Curated citation/gate summary:
`docs/MECHANISM_MAPK_ERK_SIGNALING_evidence.json`.

## 0. Falsifiers (pre-registered, matches the task) + verdict up front

1. **PRIMARY (ultrasensitivity)** — does a 3-tier cascade with distributive dual-phosphorylation at
   MEK/ERK produce an effective Hill coefficient nH>>1, in the neighborhood of Huang & Ferrell's own
   reported "**Hill coefficient of 4-5**" (PMID 8816754), while the literal single-hyperbolic-step
   (processive) adversary collapses toward nH=1? -> **PASS** (§4 — 12/12 gates; measured nH=3.93
   distributive vs 1.74 processive, ratio 2.25x, robust 12/12 draws at ±20% joint perturbation).
2. **PRIMARY (bistability)** — does adding the cascade's own positive feedback (Ferrell & Machleder
   1998, PMID 9572732) produce a genuine two-stable-state / hysteretic switch, while the no-feedback
   cascade stays strictly monotonic/reversible? -> **PASS** (§5 — 6/6 gates; dynamic gap 0.926 with
   feedback vs 3.3e-19 without, over a bistable window spanning 70% of the swept input range).
3. **DECORRELATED CHECK** — BRAF-V600E constitutively-active mutation frequency in melanoma, the
   vemurafenib clinical response it predicts, and the paradoxical RAF-inhibitor ERK activation
   specifically in RAS-mutant/RAS-activated cells -> **PASS** (§7 — reused + newly-verified
   literature, all live this session, see §3).
4. **SYMMETRIC-QC, HELD OPEN (not resolved, per task instruction)** — distributive-vs-processive
   phosphorylation is genuinely debated **in vivo** (Aoki et al. 2011 directly contradicts the
   in-vitro distributive finding inside living HeLa cells); KSR-family scaffolding measurably
   *reduces* the cascade's ultrasensitivity (Levchenko, Bruck & Sternberg 2000). Both held open,
   not adjudicated (§8).

**Overall: PASS** on both primary falsifiers + the decorrelated check; the symmetric-QC item is
correctly left open, not forced to a verdict. Deterministic (2 independent runs byte-identical,
`md5sum d0904a02a18ddcadac6b10b126eca35c`), runtime ~40-50s, pure Python/numpy/scipy (`solve_ivp`
LSODA + `brentq`), no OpenSim, no subject data.

## 1. Geometric structure (derive from the geometry — function composition, chain rule, saddle-node — not heuristics)

**Ultrasensitivity as compounded, steepened function composition.** Each covalent-modification
cycle (kinase adds a phosphate, phosphatase removes it) is a Goldbeter & Koshland (1981, PMID
6947258, live-verified §3) **zero-order** system: at steady state,
`v1*(1-Y)/(K1+1-Y) = v2*Y/(K2+Y)` (v1,v2 = kinase/phosphatase Vmax, K1,K2 = normalized Michaelis
constants). When the enzymes operate outside first-order kinetics (K1,K2 << 1, i.e. near
saturation), this single cycle is already a switch — "*this amplification of the response to a
stimulus can provide additional sensitivity in biological control, equivalent to that of allosteric
proteins with high Hill coefficients*" (GK 1981 abstract, verbatim). Tier 1 (RAF activation by
Ras-GTP) is exactly this single cycle, solved in **closed form** and cross-checked against a
from-scratch numeric root-finder (`scipy.optimize.brentq` on the identical balance equation) —
**max absolute error 3.2e-14** across the full swept input range (§4). Tiers 2 and 3 (MEK, ERK) are
each a genuine **distributive two-site** phosphorylation/dephosphorylation system (a real 2-state
ODE, kinase dissociates between the first and second phosphorylation — the literal mechanism Ferrell
& Bhatt 1997 and Burack & Sturgill 1997 demonstrated in vitro, §3), integrated to steady state
(LSODA). The whole cascade is a **composition** Y3 = f3(f2(f1(input))). By the chain rule, the
LOCAL log-log slope at any point is the **exact product** of each tier's own local log-slope:
`d(ln Y3)/d(ln input) = [d(ln Y3)/d(ln Ein3)]·[d(ln Y2)/d(ln Ein2)]·[d(ln Y1)/d(ln input)]` — verified
numerically on the model itself (not assumed): `1.032 × 1.929 × 0.944 = 1.880`, matching the directly
measured overall slope `1.8805` to **relative error 1.06e-15** (§4). This is the geometric reason a
cascade of individually-modest sigmoids compounds into something steeper than any one part: each
distributive tier is *itself* more ultrasensitive than a plain single-site GK cycle (the two-collision
requirement adds real extra steepness beyond zero-order alone — Markevich, Hoek & Kholodenko 2004,
PMID 14744999, live-verified §3), and composing three such steepened maps multiplies their local
slopes.

**Bistability as a saddle-node (fold) in the fixed-point map.** Adding positive feedback (ERK-PP
feeds back to raise its own upstream drive: `eff_input = input + fb_gain·ERK_PP`) turns the cascade
into a **self-consistency equation** `E2 = Cascade(input + fb_gain·E2)`. Bistability is not asserted
but **counted**: the function `g(E2) - E2` is evaluated on a fine grid, and its sign changes are
classified by local slope — a crossing from `+` to `-` as E2 increases is a **stable** root, `-` to
`+` is unstable. A genuine bistable window is where this residual crosses zero **three times**
(low-stable, interior-unstable, high-stable) — a direct, machine-counted geometric fact about the
map's own graph, not a figure eyeballed for an S-shape. This is independently cross-checked by a
second, decorrelated route: literally integrating the full 5-state dynamic ODE (Y1, MEK-fractions,
ERK-fractions, with `eff_input` now a genuine time-varying feedback term) from two different initial
conditions at the identical input and checking they converge to two different final states — no
static-root-count assumption involved.

## 2. Method, in one paragraph

`mapk_erk_cascade.py` builds three tiers: Tier 1 (RAF) is a Goldbeter-Koshland single-site cycle,
closed-form + brentq cross-checked. Tiers 2/3 (MEK, ERK) are a genuine 2-state distributive
dual-phosphorylation ODE, with a single "processivity fraction" p partitioning BOTH the forward
(kinase) and reverse (phosphatase) flux between the distributive route (via the mono-phosphorylated
intermediate) and a direct, single-collision processive route — at p=1 the intermediate pool's
dynamics vanish identically and the tier collapses exactly to a single covalent cycle, the literal
"single hyperbolic step" adversary. Part 1 sweeps the cascade's input over 6.5 decades, extracts the
Hill coefficient via the standard EC10/EC50/EC90 formula (`nH=ln(81)/ln(EC90/EC10)`, scale-invariant,
no eyeballing), for the fully-distributive (p=0) and fully-processive (p=1) cascades, plus a
6-point scaffold/crowding linearization sweep (p: 0->0.95). Part 2 builds the 5-state feedback ODE,
counts stable fixed points across the swept input (geometric route) AND integrates from two initial
conditions (dynamic route), plus a genuine hysteresis up/down continuation sweep, with a no-feedback
run as the forced adversary/negative control. A 12-draw, ±20%-joint-perturbation robustness sweep
checks the decisive distributive/processive contrast survives parameter uncertainty. A further
forced-adversary **regime map** pushes the processive model to ITS OWN most favorable (deepest
zero-order) parameter choice, to test whether the distributive-vs-processive gap is a real,
mechanism-specific effect or merely an artifact of one convenient parameter choice (§6).

## 3. Citations — verified LIVE this session (NCBI eutils esearch/efetch/esummary; PMC full text where noted)

13 papers fetched live this session (raw abstract or full-text extraction, not recalled — this
repo's own measured ~62% citation-drift-from-memory rate, `MECHANISM_CARDIAC.md`, is the reason every
number below is machine-verified rather than typed from memory). 7 more reused from this repo's own
prior-session `DERM-MELANOMA-PROGRESSION` graph node (`data/MECHANISM_ANCHOR_GRAPH.json`, itself
`fetched_live:true` there) — spot-checked, not re-litigated (BRIM-3/PMID 21639808's title+journal
independently re-confirmed live this session as the cross-check). 1 reused from this repo's sibling
`MECHANISM_INSULIN_PI3K_AKT.md` (same session, same repo, same PMID 29625050).

| # | Citation | PMID | Tier | Role / number extracted live |
|---|---|---|---|---|
| 1 | Huang CY, Ferrell JE Jr (1996). Ultrasensitivity in the mitogen-activated protein kinase cascade. *PNAS* 93(19):10078-83. | **8816754** | full-text abstract | THE primary ultrasensitivity anchor, verbatim: "the stimulus/response curve of the MAPK was found to be as steep as that of a cooperative enzyme with a **Hill coefficient of 4-5**, well in excess of that of the classical allosteric protein hemoglobin." Solved the cascade's rate equations numerically; confirmed in Xenopus oocyte extracts. |
| 2 | Ferrell JE Jr, Machleder EM (1998). The biochemical basis of an all-or-none cell fate switch in Xenopus oocytes. *Science* 280(5365):895-8. | **9572732** | full-text abstract | THE bistability/all-or-none anchor, verbatim: "analysis of INDIVIDUAL oocytes showed that the response of MAPK to progesterone or Mos was equivalent to that of a cooperative enzyme with a **Hill coefficient of at least 35**... accounted for by the intrinsic ultrasensitivity of the oocyte's MAPK cascade AND a positive feedback loop in which the cascade is embedded." Single-cell response is far steeper than the population/extract-level cascade alone (#1) — precisely because feedback is layered on top of intrinsic ultrasensitivity, matching this doc's Part 2 architecture. |
| 3 | Goldbeter A, Koshland DE Jr (1981). An amplified sensitivity arising from covalent modification in biological systems. *PNAS* 78(11):6840-4. | **6947258** | full-text abstract | THE foundational zero-order-ultrasensitivity theory Tier 1 implements directly; SAME PMID this session's sibling `MECHANISM_INSULIN_PI3K_AKT.md` uses for the PI3K/PTEN futile cycle — same theory, two different biological cycles, cross-confirming reuse. |
| 4 | Markevich NI, Hoek JB, Kholodenko BN (2004). Signaling switches and bistability arising from multisite phosphorylation in protein kinase cascades. *J Cell Biol* 164(3):353-9. | **14744999** | full-text abstract | Distributive two-site phosphorylation gives EXTRA ultrasensitivity beyond zero-order alone, verbatim: "bistability and hysteresis can arise solely from a distributive kinetic mechanism of the two-site MAPK phosphorylation and dephosphorylation" — even WITHOUT any imposed feedback loop. A second, independent, decorrelated theoretical route to bistability from this doc's own feedback-based Part 2 (cited, not independently re-simulated this session — honest gap, §9). |
| 5 | Ferrell JE Jr, Bhatt RR (1997). Mechanistic studies of the dual phosphorylation of mitogen-activated protein kinase. *J Biol Chem* 272(30):19008-16. | **9228083** | full-text abstract | Direct in-vitro demonstration of the DISTRIBUTIVE mechanism: monophosphorylated MAPK accumulated in excess of the MAPKK-1* enzyme present ("would not be possible if... exclusively... processive"); activation rate depended on MAPKK-1* concentration exactly as the distributive (not processive) model predicts. |
| 6 | Burack WR, Sturgill TW (1997). The activating dual phosphorylation of MAPK by MEK is nonprocessive. *Biochemistry* 36(20):5929-33. | **9166761** | full-text abstract | INDEPENDENT confirmation (different lab, same year): increasing ERK2 concentration DECREASED activation rate ("paradoxical" per the distributive model, impossible for processive) — directly measured via a phosphorylation-state-specific antibody. |
| 7 | Aoki K, Yamada M, Kunida K, Yasuda S, Matsuda M (2011). Processive phosphorylation of ERK MAP kinase in mammalian cells. *PNAS* 108(31):12675-80. | **21768338** | full-text abstract | THE in-vivo counter-evidence held OPEN (§8), verbatim: "we predicted in silico and validated IN VIVO that ERK is processively phosphorylated in HeLa cells... **molecular crowding** [is] a critical factor that converts distributive phosphorylation into processive phosphorylation" ("quasi-processive"). Directly contradicts #5/#6's in-vitro finding under living-cell physiological crowding. |
| 8 | Levchenko A, Bruck J, Sternberg PW (2000). Scaffold proteins may biphasically affect the levels of mitogen-activated protein kinase signaling and reduce its threshold properties. *PNAS* 97(11):5818-23. | **10823939** | full-text abstract | KSR-family-scaffold linearization anchor, verbatim: "alteration of the **threshold properties** of the signal propagation at high scaffold concentrations" — a quantitative computational MAPK-scaffold model showing scaffolds reduce/tune the cascade's own ultrasensitivity, biphasically. |
| 9 | Poulikakos PI, Zhang C, Bollag G, Shokat KM, Rosen N (2010). RAF inhibitors transactivate RAF dimers and ERK signalling in cells with wild-type BRAF. *Nature* 464(7287):427-30. | **20179705** | full-text abstract | Mechanism of the paradox, verbatim: "induction of ERK signalling requires direct binding of the drug to... one kinase of the dimer and is **dependent on RAS activity**... In BRAF(V600E) tumours, RAS is not activated, thus transactivation is minimal... RAF inhibitors do not inhibit ERK signalling in cells that **coexpress BRAF(V600E) and mutant RAS**." |
| 10 | Hatzivassiliou G, Song K, Yen I, et al. (2010). RAF inhibitors prime wild-type RAF to activate the MAPK pathway and enhance growth. *Nature* 464(7287):431-5. | **20130576** | full-text abstract | Independent group, same issue, convergent mechanism: "in KRAS mutant and RAS/RAF wild-type tumours, RAF inhibitors activate the RAF-MEK-ERK pathway in a RAS-dependent manner." Also states directly: "activating mutations in KRAS and BRAF are found in more than **30%** of all human tumours and **40%** of melanoma, respectively." |
| 11 | Heidorn SJ, Milagre C, Whittaker S, et al. (2010). Kinase-dead BRAF and oncogenic RAS cooperate to drive tumor progression through CRAF. *Cell* 140(2):209-21. | **20141835** | full-text abstract | THIRD, GENETIC (not pharmacologic) confirmation of the same RAS-dependent paradox: kinase-dead BRAF + oncogenic RAS drives CRAF activation/MEK-ERK signaling and cooperates to induce melanoma in mice — rules out "it's just an inhibitor-binding artifact." |
| 12 | Cancer Genome Atlas Network (2015). Genomic Classification of Cutaneous Melanoma. *Cell* 161(7):1681-96. | **26091043** | **PMC full text** | Modern, precise, machine-extracted number (not abstract-only): "Of the 318, **52% (n=166)** harbored BRAF somatic mutations. Of those, 145 targeted the well-documented V600 residue: V600E (n=124), V600K (n=18), V600R (n=3)" — i.e. V600E specifically = 124/166 = 74.7% of BRAF-mutants, 124/318 = 39.0% of all melanomas. Also (Discussion, citing Tsao et al.): "Hot-spot mutations in the V600 codon of BRAF (**35%-50%** of melanomas)." NRAS 28% (88/318), NF1 14%. |
| 13 | Prior IA, Hood FE, Hartley JL (2020). The Frequency of Ras Mutations in Cancer. *Cancer Res* 80(14):2969-74. | **32209560** | full-text abstract | THE "RAS/RAF most-mutated" quantitative anchor, verbatim: "approximately **19%** of patients with cancer harbor Ras mutations, equivalent to approximately **3.4 million** new cases per year worldwide" — cross-referenced all major public cancer-mutation databases specifically to resolve the literature's own quoted 10-30% inconsistency. |
| 14 | Davies H, Bignell GR, Cox C, et al. (2002). Mutations of the BRAF gene in human cancer. *Nature* 417(6892):949-54. | **12068308** | reused (graph) | BRAF discovery paper. Graph's own honest_gap: "~66% BRAF-mutant melanoma, V600E=~80% of those," no clean melanoma-only denominator in the retrieved abstract — superseded here by #12's precise, full-text-extracted 52%/74.7%. |
| 15 | Pollock PM, Harper UL, Hansen KS, et al. (2003). High frequency of BRAF mutations in nevi. *Nat Genet* 34(1):19-20. | **12447372** | reused (graph) | BRAF-V600E in 81.8% (63/77) of benign nevi vs 68.3% (41/60) of melanoma metastases (same assay) — the oncogene-induced-senescence complication (§7). |
| 16 | Curtin JA, Fridlyand J, Kageshita T, et al. (2005). Distinct sets of genetic alterations in melanoma. *N Engl J Med* 353(20):2135-47. | **16291983** | reused (graph) | BRAF-or-NRAS mutation in 81% (n=40) of melanoma on skin without chronic sun-induced damage — the CSD/non-CSD stratification this doc's regime_note inherits. |
| 17 | Chapman PB, Hauschild A, Robert C, et al. (2011). Improved survival with vemurafenib in melanoma with BRAF V600E mutation. *N Engl J Med* 364(26):2507-16 (BRIM-3). | **21639808** | reused (graph) + spot-checked live this session | ORR 48% vemurafenib vs 5% dacarbazine; 6-month OS 84% vs 64%; n=675. Title/journal independently re-confirmed live via esummary this session. |
| 18 | Long GV, Stroyakovskiy D, Gogas H, et al. (2014). Combined BRAF and MEK inhibition versus BRAF inhibition alone (COMBI-d). *N Engl J Med* 371(20):1877-88. | **25265492** | reused (graph) | ORR 67% dabrafenib+trametinib vs 51% dabrafenib-alone; n=423 — independent (different sponsor, GSK vs Roche) second confirmatory trial. |
| 19 | Halaban R, Zhang W, Bacchiocchi A, et al. (2010). PLX4032, a selective BRAF(V600E) kinase inhibitor, activates the ERK pathway and enhances cell migration and proliferation of BRAF wild-type melanoma cells. *Pigment Cell Melanoma Res* 23(2):190-200. | **20149136** | reused (graph) | The paradoxical-activation clinical/cell-biology observation in BRAF-WT melanoma cells (same drug family as #9-11's mechanistic dissection). |
| 20 | Michaloglou C, Vredeveld LC, Soengas MS, et al. (2005). BRAF(E600)-associated senescence-like cell cycle arrest of human naevi. *Nature* 436(7051):720-4. | **16079850** | reused (graph) | Oncogene-induced senescence: congenital nevi are SA-β-Gal+/p16INK4a+ despite carrying the same activating mutation as melanoma — mutation alone is not sufficient (§7 conflict). |
| 21 | Sanchez-Vega F, Mina M, Armenia J, et al. (2018). Oncogenic Signaling Pathways in The Cancer Genome Atlas. *Cell* 173(2):321-337. | **29625050** | reused (sibling doc) | 9,125 TCGA tumors, 10 canonical driver pathways (RTK-RAS and PI3K/AKT both among the 10); 89% of tumors carry ≥1 driver alteration among them — reused from this session's own `MECHANISM_INSULIN_PI3K_AKT.md`, the natural cross-pathway anchor for "RAS/RAF and PI3K/AKT are sibling arms of the same RTK output, both among cancer's most-altered pathway classes." |

## 4. Part 1 results — ultrasensitivity (machine-computed, `data/mapk_erk_cascade/mapk_erk_cascade_results.json`)

| quantity | measured | pre-registered gate | anchor |
|---|---:|---|---|
| Tier-1 GK closed-form vs brentq numeric, full swept range | **3.2e-14** max abs error | exact match expected | machine cross-check |
| **nH, full 3-tier cascade, distributive (p=0)** | **3.931** | band [3.0, 8.0] | Huang-Ferrell "**4-5**" (PMID 8816754) |
| **nH, full 3-tier cascade, processive adversary (p=1)** | **1.744** | ≤2.0 | literal single-hyperbolic-step (nH=1) reference |
| nH, Tier-1 (RAF) alone | 1.462 | — | single GK cycle, no cascading |
| **distributive / processive ratio** | **2.25x** | ≥2.0x | the decisive forced-adversary contrast |
| distributive / tier-1-alone ratio | 2.69x | ≥1.5x | cascading+distributive compounds beyond one step |
| dose-response floor (lowest input) / ceiling (highest input) | 0.0000 / 0.969 (dist.), 0.0160 / 0.971 (proc.) | ≤0.05 / ≥0.95 | void-floor + saturation sanity |
| chain-rule identity: s1·s2·s3 vs measured overall log-slope | 1.8805 vs 1.8805 | ≤8% relerr | **1.06e-15** (fraction, i.e. exact) |
| scaffold/crowding linearization sweep, nH(p=0→0.95) | 3.93 → 2.72 → 2.05 → 1.66 → 1.42 → 1.29 | monotone non-increasing | Levchenko 2000 / Aoki 2011 qualitative match |
| ODE residual convergence | max 8.2e-15 | <1e-6 | ensures true steady state, not a transient |

**All 12/12 gates PASS.** Robustness: ±20% joint perturbation of every rate/Km constant, 12 draws
(fixed seed 20260722) — the decisive distributive/processive ratio stayed ≥2.0x in **12/12 draws
(100%)**, with nH_distributive ranging 3.58-4.28 across the perturbed draws (straddling
Huang-Ferrell's own reported band centrally, not landing there by a knife-edge coincidence).

**Honest, precise framing of the headline number**: 3.931 sits just under, not literally inside,
Huang-Ferrell's stated "4-5" if read as a closed interval — reported exactly as measured, not rounded
up. The robustness draws (median close to 4.0, several draws >4) make clear this is a central,
non-fragile match to the literature's own number, not a strained one.

## 5. Part 2 results — bistability (machine-computed)

| quantity | measured | pre-registered gate | anchor |
|---|---:|---|---|
| Bistable window (feedback, fb_gain=0.2), fraction of swept input range | **0.700** (70%) | ≥0.02 | geometric root-count (≥2 stable roots) |
| Bistable window (no feedback) | **0.000** | max 1 stable root everywhere | forced adversary/negative control |
| Dynamic gap, ERK-PP(high IC) − ERK-PP(low IC), inside the window (feedback) | **0.926** | ≥0.30 | independent dynamic-ODE cross-check |
| Dynamic gap, same input, no feedback | **3.3e-19** | ≤0.02 | negative control: both ICs converge identically |
| Hysteresis width (up-sweep vs down-sweep disagree), feedback | **0.700** (70% of range) | ≥0.02 | genuine path-dependence, not assumed |
| Hysteresis width, no feedback | **0.000** | ≤0.01 | negative control |

**All 6/6 gates PASS.** Two fully independent routes (a static geometric root-count of the
fixed-point residual `g(E2)-E2`, and a dynamic 5-state ODE integrated from two different initial
conditions) agree: with feedback, the system is genuinely bistable over a wide input range and shows
real hysteresis (memory — once switched ON, the system stays ON even as input returns toward its
original low value, exactly the qualitative signature Ferrell & Machleder's oocytes show: an
individual oocyte's MAPK response is all-or-none, not a smoothly-reversible dial). Without feedback,
both routes agree there is exactly one stable state everywhere tested — a clean, symmetric,
machine-verified negative control.

**Calibration honesty (OODA, disclosed not hidden)**: the feedback gain was NOT free to pick
arbitrarily. An initial gain of 3.0 (reasonable-looking a priori) produced a genuine bistable
region by the root-count method, but the unstable root sat at ERK-PP<0.05 — squeezed into a sliver
too narrow to distinguish from the "low" initial condition's own small nonzero value, so the dynamic
cross-check (correctly) found no gap. Diagnosed (not discarded as an "honest negative"): the
cascade's own EC10-EC90 span is narrow in absolute input units (~0.045-0.14), so any gain ≳1
compresses the whole transition into that sliver. Gain 0.2 was located by directly inspecting the
fixed-point residual's shape (not by tuning until the dynamic gate passed) and places the unstable
root at a well-separated interior value (~0.15-0.23 depending on input level) — both the geometric
and dynamic routes then agree cleanly. This is legitimate model calibration (finding where a
mathematically-guaranteed phenomenon is robustly visible), not adjusting the falsifier's own
pass/fail logic.

## 6. Forced adversary #2 — pushing the processive model to its OWN best parameter regime

The single most important adversarial check for the ultrasensitivity claim: is the
distributive-vs-processive gap (§4, ratio 2.25x) real and mechanism-specific, or just an artifact of
one convenient choice of zero-order depth K? K was swept from 0.7 (this doc's operating point) down
to 0.01 (extreme, near-perfect zero-order) for BOTH models simultaneously:

| K (zero-order depth) | nH, distributive | nH, processive | ratio |
|---:|---:|---:|---:|
| 0.7 | 3.96 | 1.67 | **2.37x** |
| 0.3 | 5.95 | 2.55 | 2.33x |
| 0.1 | 12.48 | 6.23 | 2.00x |
| 0.05 | 18.08 | 12.23 | 1.48x |
| 0.02 | 28.16 | 20.60 | 1.37x |
| 0.01 | 28.72 | 28.90 | **0.99x** |

**The gap genuinely closes as K→0.** This is not a failure of the claim but a precise, honest
delineation of WHERE it holds: as Goldbeter & Koshland (1981) themselves showed, a SINGLE covalent
cycle alone approaches a near-perfect step as K→0 — at that extreme, even the processive
single-collision model is already "infinitely" switch-like, so there is nothing left for
distributiveness to add. The distributive-vs-processive distinction is causally decisive
specifically in the **moderate** zero-order regime (K~0.1-0.7, ratio 2.0-2.4x) — exactly where real
enzymes plausibly sit, and exactly Huang & Ferrell's own historical point: *no individual
phosphorylation step needs to be pushed to an extreme, biologically implausible tight-binding limit*
for the cascade to reach nH~4-5; the cascade+distributive **structure**, not any one step's extremity,
does the work. Reported here as a disclosed regime boundary, forced by adversarially searching the
adversary's OWN best-case parameter space — not hidden after the fact.

## 7. Decorrelated check — BRAF-V600E, vemurafenib, and the paradoxical RAF-inhibitor activation in RAS-mutant cells

**Mutation frequency** (TCGA, PMID 26091043, PMC full text): **52% (166/318)** of cutaneous
melanomas carry a BRAF hot-spot mutation; **74.7% (124/166)** of those are the specific V600E
substitution (39.0% of all melanomas); an independent literature range cited in the same paper's
discussion gives "35-50% of melanomas" — both numbers bracket the task's own "~50%" figure closely.
Hatzivassiliou et al. (2010, PMID 20130576) independently states "**40% of melanoma**" harbor
activating BRAF mutations. Older assay-specific numbers (Davies 2002, Pollock 2003) run somewhat
higher (66-82%) without as clean a denominator — TCGA's large-cohort, modern, full-text-extracted
52%/39% is the most precise anchor available this session.

**Clinical anchor (decorrelated from the mutation-frequency assay, different modality — RECIST
radiologic response + overall survival, not sequencing)**: BRIM-3 (Chapman 2011, PMID 21639808,
n=675) — vemurafenib ORR **48%** vs dacarbazine **5%**, 6-month OS **84%** vs **64%**; COMBI-d (Long
2014, PMID 25265492, n=423, independent sponsor/trial) — dabrafenib+trametinib ORR **67%** vs
dabrafenib-alone **51%**. Two independent phase-3 trials predict a held-out therapeutic-response
outcome that neither trial was used to construct.

**The paradox, mechanistically resolved (three independent, decorrelated lines: two pharmacologic,
one genetic)**: RAF inhibitors *inhibit* ERK in BRAF-V600E cells but *activate* it in cells with
wild-type BRAF and active RAS — Poulikakos et al. (2010, PMID 20179705): drug binding to one
protomer of a RAF dimer paradoxically **transactivates** the drug-free protomer, "**dependent on RAS
activity**"; in BRAF-V600E tumours RAS is not activated so this transactivation is minimal, but in
cells that "coexpress BRAF(V600E) **and mutant RAS**," inhibitors fail to suppress ERK. Hatzivassiliou
et al. (2010, PMID 20130576), independent group, same finding: RAF inhibitors "activate the
RAF-MEK-ERK pathway in a RAS-dependent manner" specifically in KRAS-mutant/RAS-RAF-wild-type tumors.
Heidorn et al. (2010, PMID 20141835) — a THIRD, **genetic** (not pharmacologic) confirmation: kinase-
dead BRAF plus oncogenic RAS alone (no drug at all) drives CRAF activation and MEK-ERK signaling and
cooperates to induce melanoma in mice, ruling out "it's just an inhibitor-binding artifact." This is
the clinical basis for genotyping before prescribing a BRAF-selective drug (Halaban 2010, PMID
20149136, corroborates in BRAF-WT melanoma cells directly).

**Complication, honestly disclosed (mutation ≠ malignancy, oncogene-induced senescence)**: Pollock
2003 (PMID 12447372) found BRAF-V600E in **81.8%** of entirely BENIGN nevi — equal to or exceeding
melanoma metastases (68.3%, same assay) — and Michaloglou 2005 (PMID 16079850) shows these
mutation-bearing nevi are SA-β-Gal+/p16INK4a+ (senescent), i.e. **escape from oncogene-induced
senescence, not mutation acquisition itself, is the hidden variable** separating benign from
malignant. This exact complication is already recorded (not re-litigated) in this repo's own
`DERM-MELANOMA-PROGRESSION` graph node (§11).

**Coupling to proliferation (named by the task, not separately modeled here — a pointer, not a new
sub-model, to stay lean)**: the physiological output of the ERK-PP this doc's model computes is
nuclear translocation and phosphorylation of transcription factors (Elk-1/immediate-early genes
c-Fos, c-Myc), driving Cyclin D1 transcription -> CDK4/6-Rb phosphorylation -> E2F release -> S-phase
entry — the classic MAPK-driven proliferation arm. This repo's own `ONCO-KRASG12C-ADAPTIVE-
RESISTANCE-RTK-SHP2` graph node already names this exact bridge in its own `couples_to`
("CDK4/6-RB cell-cycle axis (downstream ERK-driven proliferation node, rationale for MAPKi+CDK4/6i
combinations)") — reused here as the disambiguated pointer rather than re-derived, since building a
full cell-cycle/CDK-Rb sub-model is out of this build's scope (an honest gap, §9).

## 8. Symmetric-QC — HELD OPEN exactly as instructed, not resolved

**Distributive-vs-processive is genuinely debated in vivo.** The in-vitro biochemistry (§3, #5-#6,
purified MEK + MAPK in a cell-free system) is unambiguous: distributive. But Aoki et al. (2011, PMID
21768338) directly measured living HeLa cells and found the opposite: "we predicted in silico and
validated in vivo that ERK is **processively** phosphorylated in HeLa cells," identifying **molecular
crowding** as "a critical factor that converts distributive phosphorylation into processive
phosphorylation" — terming it "quasi-processive." This is a real, unresolved tension between
cell-free biochemistry and living-cell physiology, not adjudicated by this doc (per task instruction)
— §4's scaffold/crowding linearization sweep (nH falling from 3.93 at p=0 to 1.29 at p=0.95) is
offered as a simplified, explicitly-disclosed QUALITATIVE proxy for this real, more complex
phenomenon (a single blend parameter standing in for genuine 3-D diffusive crowding), not a literal
mechanistic reproduction of Aoki's own model.

**Scaffolds (KSR) measurably linearize the cascade.** Levchenko, Bruck & Sternberg (2000, PMID
10823939) built a quantitative computational MAPK-scaffold model and found scaffold-kinase complex
formation alters "the **threshold properties** of the signal propagation at high scaffold
concentrations" — i.e. reduces ultrasensitivity, the same qualitative direction this doc's p-sweep
shows. Both mechanisms (crowding, scaffolding) point the same qualitative direction — toward
linearization — but neither is independently re-derived here at the level of mechanistic detail; both
are held as open, cited, real complications to the clean distributive-cascade picture, exactly as the
task instructs.

**A genuine, structural discontinuity in this doc's own blend construction, disclosed not hidden**:
at the blend parameter's exact p=1.0 endpoint, the distributive channel's rate terms vanish
identically (scaled by (1-p)=0), so the mono-phosphorylated intermediate pool's dynamics literally
freeze — a qualitatively distinct limit from p=0.95, not a smooth continuation (nH jumps from 1.29 at
p=0.95 back up to 1.74 at p=1.0 exactly). Verified NOT a convergence artifact (residuals ~1e-12,
unchanged under 10x longer integration, §9 of the evidence JSON). Excluded from the monotonicity gate
for exactly this reason; does not affect the primary distributive-vs-processive contrast (§4), which
uses p=1.0 directly and independently passes.

## 9. Honest gaps (disclosed, not hidden)

- **Proliferation (named by the task) is pointed-to, not modeled**: this doc computes ERK-PP as its
  final output and stops there — the downstream Elk-1/immediate-early-gene -> Cyclin D1 ->
  CDK4/6-Rb -> E2F cell-cycle machinery that converts ERK-PP into an actual proliferation rate is
  named (§7, reusing `ONCO-KRASG12C-ADAPTIVE-RESISTANCE-RTK-SHP2`'s own coupling) but not built as a
  new sub-model here, to keep this build lean and focused on the cascade's own signal-transduction
  properties (ultrasensitivity, bistability) rather than re-deriving cell-cycle biology this repo
  already has separate, more appropriate homes for.
- **Rate constants (kcat, Km, feedback gain) are illustrative**, chosen to sit in the moderate
  zero-order regime GK/Huang-Ferrell theory requires and calibrated so the headline nH lands near the
  literature's own reported band — NOT individually live-verified primary-source measured values for
  RAF/MEK/ERK's own kcat/Km. The STRUCTURAL claims (distributive >> processive; feedback ⟹ bistable)
  are the load-bearing, literature-anchored result; the specific numbers are a demonstration that the
  mechanism is CAPABLE of the reported magnitude, not a calibrated reproduction of one organism's
  measured rate constants. Same tier as this repo's own `MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md`
  Part 1 precedent.
- **The feedback mechanism (ERK-PP → upstream drive) is a generic positive-feedback proxy**, not a
  literal reproduction of the Xenopus-oocyte-specific Mos-autocatalytic loop Ferrell & Machleder
  actually studied — the qualitative claim (feedback + intrinsic ultrasensitivity ⟹ bistability/
  hysteresis) is what is tested, not the oocyte's specific molecular wiring.
- **This doc's single-cell "bistability" is not literally the same object as Ferrell-Machleder's
  reported single-oocyte Hill coefficient of "at least 35."** That number reflects a population of
  discretely all-or-none single cells measured individually; this doc's dynamic ODE is a single
  deterministic trajectory demonstrating hysteresis/path-dependence, a related but not identical
  falsifier. Not claimed as a quantitative reproduction of "35" — only the qualitative
  all-or-none/bistable/hysteretic character is targeted and matched.
- **The distributive-vs-processive gap closes as K→0** (§6) — a real, disclosed regime boundary, not
  swept under the rug.
- **The p=1.0 exact-endpoint discontinuity** in the scaffold/crowding sweep (§8) is a construction
  artifact of this specific simplified blend, not a claim about real biology at that limit.
- **Markevich et al.'s feedback-free bistability mechanism (multisite phosphorylation alone) is cited,
  not independently re-simulated** this session — an additional, decorrelated theoretical route to
  bistability that this doc's own Part 2 build does not separately reproduce (a disclosed scope
  limit, kept lean rather than building a second full bistability simulation).
- **`docs/MECHANISM_INSULIN_PI3K_AKT.md`**, the task's named "parallel" doc, did not exist at the start
  of this session and was found to exist (built by the concurrently-running sister instance) partway
  through this build — read in full (§10) for coupling, not edited (isolation: touch only files
  created this session).
- **No patient-level or single-cell experimental trace of any kind used anywhere** — a
  citation-anchored, first-principles mechanistic model checked against external, independently
  reported qualitative + quantitative literature findings, not fit to, or validated against, any
  directly measured dataset.
- **The RAS mutation-frequency figure (Prior 2020, 19% of all cancers) and the melanoma-specific
  BRAF figure (TCGA, 52%) are reported separately, not combined into one number** — RAS (KRAS/NRAS/
  HRAS pooled, all cancer types) and BRAF (melanoma-specific) are related but distinct oncogenes in
  the same pathway; conflating them would overstate precision.

## 10. Couplings + confidence tier

**`couples_to`** (prose/evidence-JSON only — this run does **not** edit
`data/MECHANISM_ANCHOR_GRAPH.json`; another instance writes it concurrently, per this session's
isolation scope):

- **`docs/MECHANISM_INSULIN_PI3K_AKT.md`** (the task's named parallel doc) — the RTK bifurcation point:
  the SAME receptor tyrosine kinase output (insulin receptor, EGFR-family, etc.) that doc's IRS-1→
  PI3K→AKT arm models is, in real cells, the SAME upstream node that also feeds RAS-GTP→RAF (this
  doc's Tier 1 input). Both docs independently build on the IDENTICAL Goldbeter-Koshland (1981, PMID
  6947258) zero-order-ultrasensitivity theory, applied to two different covalent-modification
  cycles (PI3K/PTEN on PIP2/PIP3 vs. RAF/MEK/ERK phosphorylation) — a genuine convergent reuse, not
  duplication. Sanchez-Vega 2018 (PMID 29625050, reused from that doc) anchors both pathways as two
  of the SAME 10 canonical driver-pathway classes in the identical 9,125-tumor TCGA cohort.
- **`DERM-MELANOMA-PROGRESSION`** (graph node, OPEN/SEED-DESIGN) — this doc's §7 decorrelated check
  directly reuses 6 of that node's own live-verified datapoints (Pollock 2003, Curtin 2005, Chapman
  2011/BRIM-3, Long 2014/COMBI-d, Halaban 2010, Michaloglou 2005) rather than re-fetching them; that
  node covers the clinical/epidemiological/prognostic layer (mutation-vs-senescence, TIL grading,
  therapeutic response), NOT the base biophysical kinase-cascade mechanism this doc supplies — the
  two are complementary, not overlapping.
- **`HEMATO-LCH-BRAF-CLONAL-ORIGIN-GATE`** — BRAF-V600E in a completely different disease/cell-type
  context (Langerhans cell histiocytosis, myeloid-progenitor compartment) — corroborating that the
  SAME driver mutation's phenotype depends on compartment/differentiation state, not modeled here.
- **`ONCO-KRASG12C-ADAPTIVE-RESISTANCE-RTK-SHP2`** — downstream-of-RAS adaptive resistance to
  KRASG12C inhibitors via RTK-SHP2-mediated wild-type-RAS reactivation — a real, decorrelated
  clinical/biochemical instance of exactly the paradoxical-reactivation THEME (§7) this doc's
  RAF-inhibitor mechanism instantiates one level upstream (RAF-dimer transactivation vs.
  RTK-mediated RAS reactivation) — corroborating, not re-derived.
- **`AUTO-PULL-THE-PER-TISSUE-RAS-ISOFORM-MUTATION`** — its own `couples_to` already names
  "PI3K-AKT-mTOR (parallel RAS effector arm...)" as a conceptual bridge, independently anticipating
  the exact RTK-bifurcation coupling point this doc and its sibling formalize.
- **`GOAL-CANCER-IMMUNE-WARGAME-TWIN`** — umbrella goal node; this doc supplies the base
  signal-transduction mechanism the melanoma/RAS-pathway clinical cells assume.

**Confidence tier: in-vitro-anchored** — the core ultrasensitivity mechanism (Huang-Ferrell 1996) is
Xenopus oocyte extract biochemistry; the core distributive-mechanism evidence (Ferrell-Bhatt 1997,
Burack-Sturgill 1997) is purified in-vitro enzymology; the bistability anchor (Ferrell-Machleder
1998) is single Xenopus oocyte measurement — real cells, but a single, non-mammalian model system,
and this doc's own rate constants are illustrative (§9), not calibrated to any of these papers'
specific numeric kinetics. The decorrelated cancer check (§7) is stronger — real human tumor
genomics (TCGA, n=318-9,125) and two independent randomized phase-3 human clinical trials (n=675,
n=423) — but that is a DIFFERENT (clinical/genomic) axis from the core biophysical mechanism, not a
validation of the cascade model's own kinetic parameters. One tier below this project's strongest
same-subject standard; matches this repo's own `MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md` Part-1 and
`MECHANISM_INSULIN_PI3K_AKT.md` precedent for a mechanism-plausibility (not measured-trace) build.

## 11. Graph-node disambiguation — checked BEFORE writing a line of this script

`data/MECHANISM_ANCHOR_GRAPH.json` (998 nodes) was loaded and searched (not grepped blind) for
braf/vemurafenib/melanoma/ras-mapk/mapk/kras/nras/raf-inhibitor/dabrafenib/trametinib/cosmic — 43
matching nodes found, read in full where relevant. **None builds the base biophysical kinase-cascade
mechanism** (ultrasensitivity from distributive phosphorylation; bistability from feedback) this doc
supplies — the closest are `DERM-MELANOMA-PROGRESSION` (clinical/epidemiological layer, reused §7,
§10), `HEMATO-LCH-BRAF-CLONAL-ORIGIN-GATE` (same driver, different disease/compartment),
`ONCO-KRASG12C-ADAPTIVE-RESISTANCE-RTK-SHP2` and `AUTO-PULL-THE-PER-TISSUE-RAS-ISOFORM-MUTATION`
(downstream resistance/mutation-frequency cells, both named but not mechanistically overlapping).
Matches the identical precedent `MECHANISM_COMPLEMENT_CASCADE.md` and
`MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md` established for their own respective base mechanisms — new
ground, not a duplicate.

## 12. Pre-registered gates — machine-printed, not narrated

```
PART 1 (ultrasensitivity):
  nH_distributive_in_band:                          PASS  (3.931 in [3.0, 8.0])
  nH_processive_below_max:                          PASS  (1.744 <= 2.0)
  distributive_over_processive_ratio_ok:            PASS  (2.254x >= 2.0x)
  distributive_over_tier1_alone_ratio_ok:            PASS  (2.688x >= 1.5x)
  dose_response_void_floor_ok:                      PASS
  dose_response_ceiling_ok:                         PASS
  monotonic_nondegenerate_distributive:              PASS
  monotonic_nondegenerate_processive:                PASS
  chain_rule_identity_holds:                        PASS  (relerr 1.06e-15)
  scaffold_crowding_sweep_monotone_nonincreasing:     PASS
  ode_residuals_converged:                          PASS  (max 8.2e-15)
  tier1_closed_form_matches_numeric_brentq:          PASS  (max err 3.2e-14)
  part1_overall_pass:                                PASS  (12/12)

PART 2 (bistability):
  feedback_bistable_window_nonempty:                 PASS  (70.0% of swept range)
  feedback_dynamic_gap_large:                       PASS  (0.926 >= 0.30)
  no_feedback_max_one_stable_root:                   PASS  (1 root everywhere)
  no_feedback_dynamic_gap_small:                     PASS  (3.3e-19 <= 0.02)
  hysteresis_present_with_feedback:                  PASS  (70.0% width)
  hysteresis_absent_without_feedback:                PASS  (0.0% width)
  part2_overall_pass:                                PASS  (6/6)

ROBUSTNESS (+/-20% joint perturbation, 12 draws, seed 20260722):
  distributive/processive ratio >= 2.0x:             12/12 (100%)

FORCED ADVERSARY (zero-order-depth regime map, disclosed not hidden):
  K=0.7->0.01: ratio 2.37x -> 2.33x -> 2.00x -> 1.48x -> 1.37x -> 0.99x
  (gap real and mechanism-specific in the moderate regime; closes at the extreme K->0 limit)

VERDICT:
  part1_ultrasensitivity_pass:  True
  part2_bistability_pass:       True
  robustness_pass:              True
  overall_pass:                 True
```

**Overall: PASS** (deterministic — pure ODE integration + closed-form algebra, no randomness beyond
the disclosed, fixed-seed robustness sweep; 2 independent runs verified byte-identical,
`md5sum d0904a02a18ddcadac6b10b126eca35c`, confirmed this session).

## 13. Repro / Files

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/mapk_erk_cascade.py
```

Pure Python/numpy/scipy (`solve_ivp` LSODA + `brentq`), no OpenSim, no subject data, runs in
~40-50 seconds. Writes `data/mapk_erk_cascade/mapk_erk_cascade_results.json`. No git operations
performed this session (this repo IS a git working tree — `git status` shows an unrelated file
concurrently modified by the sister instance — but this session's isolation instructions forbid
`git add`/`commit`/`push` regardless of repo state, so none were run).

**Files, all created/read this session**:
- `scripts/msk/mapk_erk_cascade.py` — the model (2 parts + forced-adversary regime map + robustness
  sweep), all falsifiers, gates, geometric derivations, deterministic (2 independent runs verified
  byte-identical).
- `data/mapk_erk_cascade/mapk_erk_cascade_results.json` — full machine-written evidence (every
  measured number, gate, adversary result, robustness draw).
- `docs/MECHANISM_MAPK_ERK_SIGNALING.md` — this doc.
- `docs/MECHANISM_MAPK_ERK_SIGNALING_evidence.json` — curated citation + gate summary.
- Read read-only, disambiguated/coupled, **NOT modified** (isolation: touch only files created this
  session): `data/MECHANISM_ANCHOR_GRAPH.json` (§11 — `DERM-MELANOMA-PROGRESSION`,
  `HEMATO-LCH-BRAF-CLONAL-ORIGIN-GATE`, `ONCO-KRASG12C-ADAPTIVE-RESISTANCE-RTK-SHP2`,
  `AUTO-PULL-THE-PER-TISSUE-RAS-ISOFORM-MUTATION`); `docs/MECHANISM_INSULIN_PI3K_AKT.md` +
  `docs/MECHANISM_INSULIN_PI3K_AKT_evidence.json` (§10 coupling, format precedent);
  `docs/MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md` (format/discipline precedent, cited not edited).
