# MECHANISM PANCREATIC EXOCRINE SECRETION — juice volume, enzyme complement, secretin-driven
# bicarbonate, and the Di Magno nonlinear-reserve threshold (2026-07-22)

Completes the upper-digestion cluster alongside the GI thread (`docs/MECHANISM_GI_ABSORPTION_TRANSIT.md`
— reads its small-bowel transit time read-only) and the bile thread
(`docs/MECHANISM_BILE_ENTEROHEPATIC.md` — fat-micelle solubilization, co-required for lipolysis, prose-
coupled). Builds pancreatic juice volume (~1.5 L/day), the enzyme complement (lipase, amylase,
trypsin(ogen)/proteases) and their normal secretory outputs, secretin-driven bicarbonate secretion (up
to ~120–140 mmol/L, neutralizing gastric acid toward duodenal pH ~6), and the digestive-capacity
**reserve** — the Di Magno 1973 finding that steatorrhea appears only when enzyme output falls below
~5–10% of normal, a genuine nonlinear threshold, not a linear response. Script:
`scripts/msk/pancreatic_exocrine.py`. Evidence: `reports/probes/pancreatic_exocrine.json`.

## 0. The falsifier this doc closes (verbatim from the task, stated before any number below)

*"Does the model reproduce the MEASURED enzyme-output-vs-maldigestion threshold (fat malabsorption/
steatorrhea appears only when lipase output falls below ~5-10% of normal — Di Magno 1973 NEJM, the
classic nonlinear reserve) AND the measured bicarbonate secretory response (secretin-stimulated HCO3
concentration/output)?"*

**Answered: yes, on both legs — 12/12 pre-registered gates PASS (§10)** — with the reserve threshold
itself, the enzyme-output data points' evidentiary tier, and the gastric-acid side of the pH mass
balance all held explicitly OPEN per the task's own symmetric-QC instruction (§11), not forced.

## 1. Scope, stated up front

**Population-parameter forward/threshold model**, same class of layer as `gi_absorption_transit.py` /
`bile_enterohepatic.py` / `glucose_insulin_minimal_model.py` — no OpenSim, no per-subject data, pure
Python/numpy on literature-anchored bands. All thresholds in §10 were fixed **before** the first run
(a genuine pre-registration, matching this repo's own established convention).

## 2. Geometric structure — a capacity/demand bottleneck, not rote algebra

This doc's central object is a **capacity-vs-demand ratio**, the same "queueing/bottleneck" family
`gi_absorption_transit.py` already uses for its own `rho = Vmax_total/(dose·k_liq)` regime split:

- Let `x` = enzyme secretory output as a fraction of normal. Let `f_th` = the threshold fraction (Di
  Magno's 5–10% band) at which secretory **capacity** exactly equals digestive **demand** within the
  available small-intestinal transit window (216 min / 3.6h, reused read-only from
  `gi_absorption_transit.py` — the "buffer time" a backlog would need to clear within). Define
  `rho(x) = x/f_th`.
- For `rho>=1` (capacity exceeds demand), essentially all substrate is digested within the transit
  buffer — undigested fraction stays flat at a small baseline, **insensitive** to reductions in `x`.
  Only once `rho<1` does an ever-larger undigested backlog escape, rising as the capacity shortfall
  `(1-rho)`. This **piecewise-linear-in-rho** ("kinked", not smooth-exponential) shape is the correct
  geometric signature of a bottleneck/rate-limiting-capacity system — and it is exactly what makes a
  ~10× reserve margin (`f_th≈0.10`) look flat, then fall off a cliff, on a linear enzyme-output axis,
  even though nothing in the underlying biochemistry is discontinuous. (A pure single-exponential
  digestion-kinetics model — undigested = `exp(-k·E·T)` — was explicitly tried and rejected first: it
  cannot simultaneously reproduce ~3% undigested at 100% output and a sharp knee near 10%, because its
  curvature is wrong everywhere; the bottleneck/capacity framing is not a decorative choice.)
- **The forced adversary is the identical functional family with one parameter changed**: the naive
  "no reserve, purely proportional" model is `rho(x)=x/f_th` with `f_th` forced to **1.0** (zero
  reserve margin). This removes the objection that the reserve-capacity shape was rigged in its own
  favor — the entire empirical question reduces to the value of **one** structural parameter (the
  reserve ratio `1/f_th`), decided by external data, not by shape-choice (§7).
- The bicarbonate concentration ceiling is **also** geometric/mechanistic, not a free parameter: Sohma,
  Gray, Imai & Argent (2001, PMID 11875259) derive a **two-stage, two-location** secretory mechanism
  along the duct axis — a proximal anion-**exchanger**-dominated segment caps luminal HCO₃⁻ at ~70 mM;
  a more distal CFTR/Cl⁻ **channel**-mediated segment pushes it further to ~150 mM. The final
  concentration is bounded by which transporter population is recruited, not asserted.

## 3. Citations — 15 sources, every one esearch+esummary+efetch (or EuropePMC/StatPearls fetch) LIVE

This repo's own measured ~62% citation-drift rate from memory is why every number below traces to a
live fetch. Citation 3 (Keller & Layer 2005) was fetched as a **full 28-page PDF** via EuropePMC's PMC-
render endpoint and read directly — page images, not a summarized intermediary — the richest single
source this session.

| # | Citation | ID | Role / verbatim number |
|---|---|---|---|
| 1 | DiMagno EP, Go VL, Summerskill WH (1973). *N Engl J Med* 288(16):813-5. | PMID **4693931** | THE original finding. No AbstractText exists in PubMed for this short report (confirmed via raw XML, disclosed). Identity independently confirmed 3 ways: MeSH major-topics (Lipase/metabolism, Trypsin/metabolism, Pancreatitis/complications→Steatorrhea/etiology, Cholecystokinin, Amino Acids, Perfusion, Stimulation-Chemical) topically match an EAA-vs-CCK design; Keller & Layer's own bibliography ref [105] quotes this exact title/journal/pages; that review's Table 4 attributes EAA-stimulus (trypsin 20%, lipase 15%) and CCK-stimulus (trypsin 10%, lipase 10%) rows to "DiMagno et al 1973". |
| 2 | DiMagno EP, Malagelada JR, Go VL, Moertel CG (1977). *N Engl J Med* 296(23):1318-22. | PMID **16213** | Companion paper. A confusable same-title-fragment match (PMID 1777547, Delchier 1991, a different paper) was caught and diffed out before selecting this one. Feeds Table 4's "mixed meal, <1%/<1%" severe-insufficiency row — the cleanest, verbatim clinically-labeled falsifier point. |
| 3 | Keller J, Layer P (2005). *Gut* 54(Suppl 6):vi1-28. | PMID **15951527**, PMCID PMC1867805 | PRIMARY quantitative anchor (full-text). Verbatim: "steatorrhoea and creatorrhoea do not occur until secretion... is decreased below 5-10% of normal" (citing the "large reserve capacity of the exocrine pancreas"); "duodenal pH is around 6 early postprandially, drops towards 5-5.5..."; Tables 1&2 (enzyme outputs/activities); Section 2.6.1 ratios; Section 3.2.6 lipase pH-inactivation. |
| 4 | Hotz J, Goberna R, Clodi PH (1973). *Digestion* 9:212-23. | PMID **4765730** | Cross-SPECIES (rat), cross-METHODOLOGY (95% surgical pancreatectomy, not secretory-output fraction) convergent anchor on the same ~5% figure. No abstract retrievable (disclosed). |
| 5 | Tian C, Ghodeif AO, Arshad S, Gillespie E. *Exocrine Pancreatic Insufficiency.* StatPearls NBK555926, updated 2025-09-14. | PMID **32310386** | Independent modern restatement, 5 decades decorrelated from citations 1/4: verbatim "~1.5 L of pancreatic fluid"; "Fat malabsorption is defined by a decrease in pancreatic lipase and trypsin levels of at least 5% to 10%"; steatorrhea = fecal fat >7g/day on 100g/day fat diet; fecal elastase-1 cutoffs (<200/<100/<50 µg/g). |
| 6 | Sohma Y, Gray MA, Imai Y, Argent BE (2001). *JOP* 2(4 Suppl):198-202. | PMID **11875259** | Mechanistic two-stage bicarbonate derivation: exchanger caps ~70mM; CFTR/Cl⁻ channels push to ~150mM. |
| 7 | Steward MC, Ishiguro H (2009). *Curr Opin Gastroenterol* 25(5):447-53. | PMID **19571747** | Verbatim: duct epithelium secretes HCO₃⁻ "at concentrations as high as 140 mmol/l" — matches the task's own upper band almost exactly. |
| 8 | Denyer ME, Cotton PB (1979). *Gut* 20(2):89-97. | PMID **428831** | THIRD, decorrelated, direct human ductal-cannulation study: peak bicarbonate "in excess of 100 mmol/l" (exact cohort for that sentence not fully disambiguated from the abstract, disclosed). |
| 9 | Chandra R, Liddle RA (2009). *Curr Opin Gastroenterol* 25(5):441-6. | PMID **19535978** | Discloses real regulatory complexity beyond secretin/CCK (vagal, ghrelin, orexin-A, NPY, melatonin, obestatin, leptin) — symmetric QC. |
| 10 | Morisset J (2008). *Pancreas* 37(1):1-12. | PMID **18580437** | Discloses genuine open question: CCK-/secretin-releasing-factor negative feedback, established in animals, is disputed in humans ("supporters and detractors"). |
| 11 | Whitcomb DC, Lowe ME (2007). *Dig Dis Sci* 52(1):1-17. | PMID **17205399** | Canonical review establishing enzyme-complement identity (proteases, lipases, amylase) — brief "meta" abstract, no line-item list extracted (disclosed). |
| 12 | Rune SJ (1972). *Gastroenterology* 62(4):533-9. | PMID **5020865** | Classical primary duodenal-pH anchor, bibliographically confirmed; no abstract retrievable (disclosed); content corroborated via citation 3's independent restatement. |
| 13 | Layer P, Go VL, DiMagno EP (1986). *Am J Physiol* 251(4 Pt 1):G475-80. | PMID **2429560** | Enzyme fate during SI transit: only 1%/74%/22% of lipase/amylase/trypsin activity reaches the terminal ileum after 50g rice starch — same ~3.6h transit window `gi_absorption_transit.py` measures. |
| 14 | Vadukoot Lazar M, C S Menon A, Thomas J (2026). *Cureus* 18(3):e105089. | PMID **41994742** | Verbatim: Type 3c diabetes "is characterized by BOTH endocrine and exocrine pancreatic insufficiency (EPI)" — the shared-organ coupling anchor to `ORG-PANCREAS-GLUCOSE-INSULIN`. |
| 15 | (in-repo) `scripts/msk/gi_absorption_transit.py` / `reports/probes/gi_absorption_transit.json` | n/a | SI mean transit 216 min (3.6h), reused read-only for the transit-buffer framing. |

## 4. Headline results (all machine-computed, `reports/probes/pancreatic_exocrine.json`)

| quantity | computed value | anchor | verdict |
|---|---:|---|---|
| Pancreatic juice volume | **1.5 L/day** | StatPearls (independent) = task central, exact | **PASS** |
| Lipase:amylase output ratio (from Table 1's own central values) | **6.0** | Keller & Layer's own stated 3–6:1 band | **PASS** (exact upper edge — disclosed, not knife-edge-hidden) |
| Lipase:trypsin output ratio | **7.5** | stated 5–10:1 band | **PASS**, comfortably interior |
| Threshold model, undigested fraction at EAA-stimulus data (lipase 15%, trypsin 20% of normal) | **2–5%** (baseline, swept) | steatorrhea cutoff 7% | **PASS** — correctly stays non-steatorrheic |
| Forced adversary (zero-reserve, same family, f_th=1.0) at the SAME EAA data | **80.4–85.75%** undigested | same 7% cutoff | **FALSIFIED** — wrongly predicts severe steatorrhea, **11.49×** over cutoff |
| Threshold model at severe/decompensated data (mixed meal, <1% of normal) | **80.4–90.5%** undigested | 7% cutoff | **PASS** — correctly predicts decisive steatorrhea |
| Bicarbonate ceiling: 3 decorrelated sources (Steward/Ishiguro 140, Sohma 150, Denyer/Cotton >100 mmol/L) | 2/3 within task's [120,140] band or 10% | task band | **PASS** |
| Two-stage bicarbonate mechanism | 70 mM (exchanger) < 150 mM (channel) | monotonic staging | **PASS** |
| Duodenal pH, early postprandial | **6.0** (Keller & Layer, direct measurement) | task target ~6 | **PASS**, exact |
| Cross-species/methodology reserve convergence | rat 95%-resection remnant = **5%** | human secretory band [5%,10%] | **PASS** |
| Void floor (E=0) | undigested = **1.0** exactly | structural | **PASS** |
| Ceiling floor (E=1× or 50× normal) | undigested = baseline exactly (**3%** at central baseline) | structural | **PASS** |

**Gates: 12/12 PASS** (§10). Held explicitly OPEN, not forced: exact reserve-threshold point value, the
EAA/CCK data points' evidentiary tier, and the gastric-acid side of the pH mass balance (§11).

## 5. Volume (§A)

StatPearls' independently-stated "~1.5 L of pancreatic fluid" per day (secretin+CCK stimulated) matches
the task's own central figure **exactly** (`abs_diff_L_day = 0.0`) — not a fitted band, a genuine
cross-source point match.

## 6. Enzyme complement + ratios (§B) — a source-internal-consistency check

Keller & Layer (2005) state TWO independent facts that should agree if their own data is self-
consistent: (1) Table 1's own central postprandial output values (lipase 4500, amylase 750, trypsin
600 U/min — early-postprandial midpoints); (2) a separately-stated ratio range (lipase:amylase ~3–6:1,
lipase:trypsin ~5–10:1, Section 2.6.1). Computing the ratios directly from (1) gives **6.0** and **7.5**
respectively — both fall inside (1) at the exact upper edge, (2) comfortably interior — the same "does
a live-fetched source's own numbers add up" discipline as `bile_enterohepatic.py`'s StatPearls-
arithmetic check.

**Disclosed caveat**: Table 2's concentrations are **duodenal juice** (diluted by co-secreted gastric
acid, bile, Brunner's-gland mucus), not pure pancreatic secretion. Dividing Table 1's output by Table
2's concentration gives an implied duodenal flow of **4.5 mL/min** postprandially — physiologically
unremarkable in isolation, but **not** directly commensurable with the 1.5 L/day pure-pancreatic-juice
figure (naively annualizing 4.5 mL/min over 24h gives 6.48 L/day, ~4× too high) — reported as an
order-of-magnitude plausibility check, not forced into a single-number match.

## 7. THE central falsifier — capacity/demand threshold model, forced adversary (§C)

Sweeping `f_th` over Di Magno's own **[5%, 10%]** band (21 grid points) and a baseline-undigested-
fraction sweep (2%, 3%, 5%; never a single point), evaluated at Di Magno's own **three** stimulation-
mode data points (reused from Keller & Layer's Table 4 / Section 3.2.1, not re-fit):

| data point | x (% of normal) | threshold model: undigested | forced adversary (f_th=1.0): undigested |
|---|---:|---:|---:|
| EAA stimulus, lipase | 15% | 2–5% | **85.3–85.75%** |
| EAA stimulus, trypsin | 20% | 2–5% | **80.4–81.0%** |
| CCK stimulus, lipase & trypsin | 10% | 2–5% | **90.2–90.5%** |
| Mixed meal (1977, severe), lipase & trypsin | <1% (used 1%) | **80.4–90.5%** | **99.02–99.05%** |

**The discriminating test is the EAA/CCK mid-range (10–20% of normal).** The threshold model predicts
near-baseline digestion there (matching the "compensated insufficiency" framing Di Magno's own papers
and Keller & Layer's review consistently use for this range) — **PASS**, robust across the entire
pre-registered sweep. The forced adversary — **the identical functional family with the reserve ratio
set to 1 instead of 10–20** (not a different shape, not a strawman) — wrongly predicts severe
steatorrhea there, **11.49× over the clinical cutoff** at its most favorable (lowest) value — a large,
non-knife-edge falsification margin. Both models correctly predict decisive steatorrhea at the
explicitly-labeled severe/decompensated (<1%) point — a necessary sanity floor, but this point alone
does **not** discriminate the two hypotheses (both agree there); the mid-range is where the reserve-
capacity claim earns its keep.

**CCK point, reported not gated**: at x=10% (Di Magno's own stated *lower* edge of the reserve band),
the threshold model's undigested fraction stays at baseline (2–5%) across the *entire* pre-registered
`f_th` sweep — because x=10% coincides with the sweep's own upper bound, `rho=x/f_th` never drops
below 1.0 for any `f_th` in [5%,10%]. This is a genuine computed result, but it is a **boundary
coincidence** (x sits exactly at the sweep's own edge), not a robust interior margin like the EAA
points enjoy — reported precisely rather than oversold as either "ambiguous" or "robustly confirmed".

## 8. Bicarbonate concentration + two-stage mechanism (§D)

Three **decorrelated** sources — one direct human ductal-cannulation primary study (Denyer & Cotton,
different methodology from the two mechanistic reviews) and two independent mechanistic/computational
reviews (Sohma et al 2001; Steward & Ishiguro 2009) — against the task's own [120,140] mmol/L band:

| source | value (mM) | vs task band |
|---|---:|---|
| Steward & Ishiguro (2009) | 140 | exact match, 0% outside |
| Sohma et al (2001), final channel-mediated stage | 150 | 7.1% over upper edge |
| Denyer & Cotton (1979), direct human ductal juice | >100 (floor) | 16.7% under lower edge |

**2/3 within band or 10%** — gated PASS. The Sohma two-stage mechanism (exchanger-only cap **70 mM** <
channel-mediated final **150 mM**, monotonic, structural not fit) explains *why* 140 mM (Steward &
Ishiguro) sits between the two Sohma stages: 140 mM is a commonly-observed achieved maximum, 150 mM the
modeled asymptotic ceiling of the **same** underlying transporter-recruitment mechanism, not two
competing claims.

## 9. Duodenal pH neutralization (§E)

Keller & Layer's own live-fetched, verbatim statement — "duodenal pH is around 6 early postprandially"
— is a **direct measurement** of the task's target quantity, not a reconstruction: `abs diff = 0.0`,
exact match. Interdigestive pH sits higher (6–7); late-postprandial pH drops to 5–5.5 as gastric acid
delivery outpaces neutralization capacity later in the meal.

**Mechanistic cross-link (why neutralization matters, not just a number)**: pancreatic lipase activity
falls to 50% at pH 7 vs pH 9 (Keller & Layer, Section 3.2.6) and is **irreversibly** inactivated below
pH 4. In severe exocrine insufficiency, intraduodenal pH falls to ~4 late postprandially (Keller & Layer
citing Di Magno's own group's work) — meaning the enzyme-threshold axis (§7) and the bicarbonate/pH
axis (§9) interact **multiplicatively**, not independently, in decompensated disease: the small residual
lipase that IS secreted gets further destroyed by the acid it should never have been exposed to.

**Illustrative, explicitly non-gated**: at the central volume (1.5 L/day) and central concentration
(130 mM, band midpoint), computed HCO₃⁻ output = **195 mmol/day** (≈65 mmol/meal, 3-meal even split).
The matching gastric-acid-output side of this titration (mEq H⁺/hour) was **not** independently live-
verified this session (disclosed gap, §11) — the duodenal-pH gate instead uses a direct measurement of
the target quantity itself, stronger evidence than an indirect mass-balance reconstruction would be.

## 10. Reserve-capacity cross-methodology convergence (§F) + void floor/ceiling (§G)

Two **independent axes**, different species, different measurement methodology, converge on the same
figure: Di Magno's **secretory-output-fraction** threshold in humans (~5–10%) vs Hotz, Goberna & Clodi's
(1973) **anatomical-resection-fraction** threshold in rats (95% pancreatectomy leaves exactly 5% of
tissue). The rat remnant fraction (**5%**) sits inside the human secretory band **[5%,10%]** — gated
PASS, a genuine over-determination (cross-species AND cross-methodology), not a restatement of the same
measurement.

**Structural checks (analytic, not numerically approximate)**: undigested fraction at zero enzyme output
= **1.0 exactly** (total malabsorption, the correct void floor); undigested fraction at 1× or 50×
normal output = the baseline **exactly** (a superphysiological enzyme dose cannot outperform the
residual non-enzymatic digestion floor) — both PASS.

```
A_volume_exact_match_task_vs_statpearls:                   PASS (1.5 == 1.5)
B_ratio_lipase_amylase_internally_consistent:              PASS (6.0, edge of [3,6])
B_ratio_lipase_trypsin_internally_consistent:              PASS (7.5, interior of [5,10])
C_threshold_model_correct_at_EAA_compensated_10_20pct:      PASS (2-5% undigested, cutoff 7%)
C_adversary_falsified_at_EAA_large_margin:                 PASS (adversary 80.4-85.75%, 11.49x cutoff)
C_threshold_model_correct_at_severe_decompensated_lt1pct:   PASS (80.4-90.5% undigested)
D_majority_bicarb_sources_confirm_120_140_band:             PASS (2/3 sources)
D_two_stage_mechanism_monotonic:                            PASS (70 < 150 mM)
E_duodenal_pH_direct_measured_match_to_target_6:            PASS (exact)
F_cross_species_cross_methodology_convergence:              PASS (5% in [5,10]%)
G_void_floor_exact:                                         PASS (1.0 exactly)
G_ceiling_floor_exact:                                      PASS (baseline exactly)

OVERALL: PASS (12/12). Deterministic -- 2 independent runs byte-identical (verified this session,
md5-compared on disk, not just asserted).
```

## 11. Symmetric QC — held OPEN, not forced (per task instruction)

- **Enzyme outputs vary hugely with meal/CCK stimulation and assay**: Keller & Layer's own Table 1
  spans a >3–6× range between interdigestive and postprandial states for the *same* normal subject —
  "normal" is itself state-dependent, never collapsed to one point anywhere in this script (every gate
  uses a swept band).
- **The reserve threshold is itself a range** (5–10%), used throughout as a sweep across 21×3 grid
  points, never one number.
- **Common-mode risk, disclosed, not hidden**: Di Magno's own 3 stimulation-mode data points (EAA/CCK/
  mixed-meal) trace to one research group across 2 papers (1973, 1977). This is only *partially*
  mitigated by Hotz's cross-species/cross-methodology convergence (§10) and by Keller & Layer's
  independent restatement of the *qualitative* rule three decades later (itself further restated by
  StatPearls in 2025) — the qualitative "5–10%" rule has 3 independent restatements; the specific
  quantitative 15–20%/10%/<1% data triplet used in the forced-adversary test does not.
- **A real evidentiary-tier asymmetry between the two falsifier legs**: the EAA/CCK points' "no
  steatorrhea" status is *inferred* from the source's own stated threshold rule and "compensated
  insufficiency" terminology, not from a directly-quoted per-patient outcome sentence for those exact
  rows. The <1% mixed-meal point, by contrast, *does* have an explicit verbatim clinical-outcome
  sentence (steatorrhea present, needing enzyme therapy). Disclosed rather than smoothed over.
- **Denyer & Cotton's ">100 mmol/l" figure's exact cohort** (healthy vs the chronic-pancreatitis
  patients also studied in that paper) is not fully disambiguated from the abstract alone — used only
  as a floor/order-of-magnitude corroboration.
- **The gastric-acid-output side of the duodenal-neutralization mass balance was not independently
  live-verified this session** (multiple targeted PubMed searches returned no directly quotable normal-
  subject mEq/hour figure) — the duodenal-pH gate uses a direct measurement of the target quantity
  instead (§9); the HCO₃⁻ mass-balance number is reported as illustrative and explicitly not gated.
- **Whitcomb & Lowe (2007)'s indexed abstract** is a brief "meta" summary without a line-item enzyme
  list — cited for enzyme-complement *identity* (textbook-grade), not numerically mined.
- **Rune (1972), DiMagno (1973), and Hotz (1973) have no retrievable AbstractText in PubMed** for this
  session (confirmed via raw XML fetch, not assumed) — all 3 cited bibliographically (title/journal/
  year/volume/pages independently verified live), content corroborated via Keller & Layer's independent
  restatement/table data instead.
- **The CCK data point (10%) is a boundary coincidence, not a robust interior result** (§7) — reported
  precisely, not oversold.
- **No direct CCK/secretin-to-insulin-secretion hormonal cross-talk** (distinct from the classical
  GLP-1/GIP incretin axis) was independently verified this session and none is claimed in §12 — only
  the shared-organ/Type-3c-diabetes structural coupling is claimed, with a live citation.
- **No subject-specific data anywhere** — a population-parametrized forward/threshold-model consistency
  check, matching every other `MECHANISM_*` systemic layer's own disclosed scope.

## 12. Couplings (prose only — not folded into the shared graph this session, per isolation scope)

- **`docs/MECHANISM_GI_ABSORPTION_TRANSIT.md` / `gi_absorption_transit.py`**: SI transit time (3.6h,
  reused read-only) is the "buffer window" §2's capacity/demand model implicitly assumes digestion
  must complete within; that script's own SGLT1/GLUT2 model covers glucose absorption only, not yet a
  lipid/protein digestion-completeness pathway this doc could feed.
- **`docs/MECHANISM_BILE_ENTEROHEPATIC.md` / `bile_enterohepatic.py`**: fat digestion requires BOTH this
  doc's pancreatic lipase (triglyceride→FFA+monoglyceride hydrolysis) AND that doc's bile-salt
  micellization (solubilizing the hydrolysis products for mucosal uptake) — neither doc alone completes
  the fat-absorption pathway; that doc's own §8 already flags this exact gap from its side. Neither has
  been quantitatively joined into one forward model this session.
- **`docs/MECHANISM_GLUCOSE_INSULIN.md` / `ORG-PANCREAS-GLUCOSE-INSULIN`**: same organ, different cell
  population (acinar/exocrine here vs islet/endocrine there). Type 3c (pancreatogenic) diabetes is the
  clinical entity at their intersection — Vadukoot Lazar et al (2026, PMID 41994742) verbatim: "is
  characterized by BOTH endocrine and exocrine pancreatic insufficiency (EPI)". Qualitative only — no
  quantitative cross-model built this session; a direct CCK/secretin-to-insulin hormonal cross-talk is
  explicitly NOT claimed (§11).

## 13. Confidence tier

**In-vivo-anchored (pancreatic function tests / Di Magno), tiered — not a single flat confidence
level**: (1) the central reserve-capacity finding traces to 2 primary Di Magno-group human papers
(1973, 1977, both bibliographically confirmed live, neither with a retrievable PubMed abstract this
session) plus 1 independent cross-species/cross-methodology primary paper (Hotz 1973, rat
pancreatectomy) plus 2 independent modern secondary restatements of the qualitative rule (Keller &
Layer 2005's full-text review, StatPearls 2025) — the qualitative 5–10% rule is very strongly
triangulated; the specific EAA/CCK/mixed-meal quantitative data triplet used in the forced-adversary
test all traces to one research lineage (disclosed common-mode risk, §11); (2) the bicarbonate-
concentration ceiling is anchored to 3 independent, decorrelated sources (1 direct human ductal-
cannulation primary study, 2 independent mechanistic/computational reviews spanning 2001–2009); (3)
the duodenal-pH-target claim uses a direct measured quantity from a live full-text-fetched review, not
a reconstruction. Population-parametrized forward/threshold model throughout — no subject-specific
data, matching this repo's other systemic-layer scope (`glucose_insulin_minimal_model.py`,
`renal_filtration.py`, `hepatic_clearance.py`, `gi_absorption_transit.py`, `bile_enterohepatic.py`).

## 14. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/pancreatic_exocrine.py
```
Pure Python/numpy, no OpenSim, no subject data, runs in under a second, deterministic (verified: 2
independent runs byte-identical, md5-compared on disk). Writes `reports/probes/pancreatic_exocrine.json`.
Read-only import of `reports/probes/gi_absorption_transit.json` (SI transit time; falls back to a
hardcoded 3.6h literature value if absent, disclosed via `used_repo_gi_json: false` in the JSON). No
git operations. Files touched this session: `scripts/msk/pancreatic_exocrine.py`,
`reports/probes/pancreatic_exocrine.json`, this doc. `data/MECHANISM_ANCHOR_GRAPH.json` (shared,
concurrently written by another instance) was not read or edited this session — no existing
`ORG-PANCREAS-EXOCRINE`-class node was found in an earlier grep of that file (only the unrelated
endocrine `ORG-PANCREAS-GLUCOSE-INSULIN`), so no disambiguation lookup against it was needed.

## 15. Files

- `scripts/msk/pancreatic_exocrine.py` — full model: 15 citations, the capacity/demand bottleneck
  derivation (§2), the enzyme-complement/ratio internal-consistency check, the central forced-adversary
  falsifier test (Di Magno's own EAA/CCK/mixed-meal data vs the zero-reserve adversary), the 3-source
  decorrelated bicarbonate check + two-stage mechanism, the duodenal-pH direct-anchor + lipase-pH-
  inactivation mechanistic cross-link, the cross-species/cross-methodology reserve-capacity convergence
  (Hotz 1973 vs Di Magno), the structural void-floor/ceiling checks, and all gates.
- `reports/probes/pancreatic_exocrine.json` — full evidence: every citation, every computed number, all
  12 gates, the falsifier verdict, the honest gaps, the confidence tier.
- Read-only, not modified: `reports/probes/gi_absorption_transit.json` (SI transit time).
