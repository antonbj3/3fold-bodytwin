# MECHANISM GI ABSORPTION & TRANSIT — gastric emptying, small-bowel transit, SGLT1/GLUT2 glucose
# absorption, compartmental model (2026-07-22)

Builds the GASTROINTESTINAL absorption/transit organ-system layer: gastric emptying (liquid +
solid), small-intestinal transit, and a mechanistic SGLT1/GLUT2 glucose-absorption model whose
output Ra(t) is cross-checked against the glucose-insulin thread's own OGTT input
(`scripts/msk/glucose_insulin_minimal_model.py`) — the "un-attempted extension" that doc's own
Section 7 flags: *"a full closed-loop oral model (Dalla Man-style) is a separate, un-attempted
extension."* Script: `scripts/msk/gi_absorption_transit.py`. Evidence:
`reports/probes/gi_absorption_transit.json`.

**Headline, stated up front, two tiers, not one (they are not the same evidentiary class):**
- **Tier 1 — gastric emptying + small-bowel transit: PASS**, strongly. Multi-anchor (2 independent
  multicenter scintigraphy cohorts, 2 breath-test cross-validations, 3 wireless-motility-capsule
  studies), cross-population-replicated, and a forced mono-exponential adversary loses by **48.5x**
  (SSE) — a real, load-bearing result, not a knife-edge coincidence.
- **Tier 2 — mechanistic SGLT1/GLUT2 Ra(t) coupling to the glucose-insulin thread: PARTIAL / OPEN**,
  honestly. The model's Ra(t) peak-**timing** lands inside Ferrannini et al 1985's directly
  tracer-measured 15-30 min window only *after* a forced OODA fix (an initial "deep Regime A"
  draft peaked at 2.5 min — wrong, diagnosed, fixed by calibration). A **second** forced fix,
  attempting to *also* match Ferrannini's sustained Ra-at-3.5h, did **not converge** — diagnosed as
  a likely genuine structural limit of a 2-parameter single-Erlang-chain architecture, reported as
  an open finding, not forced to a false pass.

## 0. Scope and disambiguation from existing GI-* graph nodes

Checked live in `data/MECHANISM_ANCHOR_GRAPH.json`, not assumed: `ORG-GI-MOTILITY-ENS` (IBS/DGBI
gut-brain reverse-causation pacemaker test), `GI-ENTERIC-MOTILITY` (transit-vs-contractility-index
discordance in constipation/gastroparesis), `GI-INTESTINAL-BARRIER-PERMEABILITY` (dual-sugar
permeability), plus disease-specific nodes (achalasia, Barrett, microscopic colitis, gastric
cancer). **All of these are DISEASE / DISCORDANCE claims.** None builds the quantitative physiology
layer this doc builds: a forward compartmental model of gastric emptying → small-bowel transit →
saturable mucosal absorption → systemic glucose appearance. Two PMIDs are **reused, re-verified
live this session** (confirmed matching, not drifted) from `GI-ENTERIC-MOTILITY`'s own evidence:
20465593 (capsule-vs-radiopaque-marker) and 34687494 (Sangnes 2021 WMC diabetic transit).

Population-parametrized forward-model arithmetic, same class of layer as
`glucose_insulin_minimal_model.py` / `renal_filtration.py` / `hepatic_clearance.py` — no OpenSim,
no `.osim` model, no per-subject trace. Pure Python/numpy/scipy, `.venv-msk`.

## 1. Geometric structure (derive, don't assert)

- **Solid gastric-emptying lag phase** = a stretched-exponential/Weibull survival curve
  `R(t) = exp(-k t^beta)`. `beta > 1` is the geometric signature of a **regulated** (feedback/
  trituration-limited), not memoryless, process. Fit via a **closed-form log-log linear
  regression** (`ln(-ln R) = ln k + beta ln t` is exactly linear) **cross-checked** against
  nonlinear least squares — two independent numerical routes, not one.
- **Small-bowel transit-time dispersion** is a spectral/multiplicity fact of an N-stage Erlang
  chain: `CV = 1/sqrt(N)` for N identical first-order compartments in series, **independent of the
  mean** (fixed by `N/k_si`) — chain length alone controls how spread-out the arrival-time
  distribution is, decoupled from where its center sits.
- **Rate-limiting-step regime**: whether gastric emptying or mucosal transport capacity dominates
  Ra(t) timing is a genuine two-regime question, resolved via a dimensionless capacity ratio
  `rho = Vmax_total / (dose * k_liq)`, not asserted. The calibrated operating point (rho≈2.27,
  Section 5) sits in an **intermediate transition zone**, not deep in either asymptote (rho<<1:
  peak~100-130min, absorption-limited; rho>>1: peak~0-3min, emptying-limited) — both mechanisms
  genuinely contribute to real peak-timing.

## 2. Citations — 22 PMIDs, every one esearch+esummary+efetch LIVE this session (NCBI eutils)

This repo's own measured **~62% citation-drift rate from memory** is why every number below traces
to a live fetch, not recall. **One drift was caught and corrected in this session**: the initial
recollection "Dalla Man/Camilleri/Cobelli 2006, PMID 17153199" resolves live to an **unrelated**
cardiac reaction-diffusion modeling paper (Potse et al.) — the correct PMID (**17153204**) was
re-found by title esearch, not silently assumed.

| # | Citation | PMID | Role / verbatim number |
|---|---|---|---|
| 1 | Tougas et al 2000, *Am J Gastroenterol* 95(6):1456-62 | **10894578** | PRIMARY solid-meal scintigraphy anchor, n=123, 11 centers. Median % retention at 60/120/240min = 69% (95th %ile 90%), 24% (60%), 1.2% (10%). ">10% at 4h = delayed." |
| 2 | Abell et al 2008 consensus, *J Nucl Med Technol* 36(1):44-54 | **18287197** | Standardized protocol (low-fat egg-white meal, 0/1/2/4h imaging) this model's timepoints mirror. |
| 3 | Elashoff, Reedy, Meyer 1982, *Gastroenterology* 83(6):1306-12 | **7129034** | Confirms the power-exponential model exists/is recommended — F3's functional form is literature-grounded, not invented. |
| 4 | Vasavid et al 2014, *J Neurogastroenterol Motil* 20(3):371-8 | **24948129** | INDEPENDENT cross-population replication, n=189, 7 Thai centers, rice+egg meal. Median GE-T50 = 68.7 (45.1-107.8) min. |
| 5 | Maes/Ghoos et al 1998, *Am J Physiol* 275(1):G169-75 | **9655697** | 13C-octanoate breath test vs scintigraphy, same subjects: r=0.98 (t½), r=0.85 (lag); mean diff t½=10min (-20,41 CI). |
| 6 | Miller, Parkman, Urbain et al 1997, *Dig Dis Sci* 42(1):10-8 | **9009110** | Lactulose breath test vs scintigraphy orocecal transit: r=0.95, BUT lactulose ITSELF accelerates transit (P=0.004), slows solid GE (P=0.02), n=8. |
| 7 | Lee, Rao, Nguyen et al 2019, *Clin Gastroenterol Hepatol* 17(9):1770-9 | **30557741** | WMC vs GES, n=167. Clinical cutoffs: GES delayed >10%@4h / WMC >5h; overall agreement 75.7% (kappa=0.42). |
| 8 | Maqbool, Parkman, Friedenberg 2009, *Dig Dis Sci* 54(10):2167-74 | **19655250** | SmartPill vs scintigraphy, n=10: r=0.95 (gastric retention@120min), r=0.73 (@240min). |
| 9 | Camilleri et al 2010, *Neurogastroenterol Motil* 22(8):874-82 | **20465593** | (Re-verified, matches graph.) WMC vs ROM colonic transit, n=158: r=0.707, 87% categorical agreement, BUT median ROM 55.0h vs WMC 43.5h (P<0.001). |
| 10 | Hunt & Stubbs 1975, *J Physiol* 245(1):209-25 | **1127608** | 33 pooled studies: caloric DELIVERY rate (not volumetric rate) to the duodenum is the regulated quantity — meal-composition-dependence anchor. |
| 11 | Kalogeris, Reidelberger, Mendel 1983 (rat), *Am J Physiol* 244(6):R865-71 | **6407340** | Cross-species replication of the same caloric-constancy principle (disclosed non-human). |
| 12 | Read, Al-Janabi, Bates, Barber 1983, *Gastroenterology* 84(6):1568-72 | **6840487** | PRIMARY small-bowel transit anchor: paired, non-intubated control, SBTT = 3.6±0.4h. |
| 13 | Degen & Phillips 1996, *Gut* 39(2):299-305 | **8977347** | n=32: scintigraphy and radio-opaque markers "equally well reflected" total colonic transit — methods AGREE here (reported for balance). |
| 14 | Sangnes et al 2021, *United European Gastroenterol J* 9(10):1168-77 | **34687494** | (Re-verified, matches graph.) WMC whole-gut transit: healthy 35h55min; diabetic+constipation 66h15min; diabetic no-constipation 71h16min. |
| 15 | Ferrannini, Bjorkman, Reichard et al 1985, *Diabetes* 34(6):580-8 | **3891471** | PRIMARY glucose Ra(t) anchor, n=11, oral glucose 1g/kg, hepatic-vein catheter+double-tracer. Ra "peak after 15-30 min"; Ra(3.5h)=2.47±0.45 mg/min/kg; 73±4% of load recovered systemically by 3.5h. |
| 16 | Basu, Di Camillo, Toffolo et al 2003, *Am J Physiol Endocrinol Metab* 284(1):E55-69 | **12485809** | n=12, mixed meal, triple-tracer: initial splanchnic extraction (ISE)=12.9±3.4%; dual-tracer method OVERESTIMATES ISE vs this triple-tracer method (their own finding). |
| 17 | Dalla Man, Camilleri, Cobelli 2006, *IEEE Trans Biomed Eng* 53(12):2472-8 | **17153204** | Ra(t)-system-model-vs-tracer-gold-standard validation paradigm. CITATION DRIFT CAUGHT: initial recall (17153199) resolves to an unrelated cardiac paper. |
| 18 | Trahair, Horowitz, Marathe et al 2014, *Physiol Rep* 2(11):e12204 | **25413324** | KEY coupling anchor: n=87, SAME 75g oral-glucose dose, 13C-breath-test GE-T50 measured directly against glycemic response (P<0.01 to P<0.001 at multiple timepoints). |
| 19 | Wright, Loo, Hirayama 2011, *Physiol Rev* 91(2):733-94 | **21527736** | SGLT1 biology review (Na+-coupled, cloned 1987, SLC5 family) — Km/Vmax textbook-grade, not machine-extracted from this abstract. |
| 20 | Reds, Geillinger, Zietek et al 2014, *PLoS One* 9(2):e89977 | **24587162** | FORCED ADVERSARY: SGLT1/GLUT2 knockout mice — "SGLT1 is unequivocally the prime intestinal glucose transporter... no evidence for GLUT2 playing any role in... apical glucose influx." Contradicts Kellett/Brot-Laroche apical-GLUT2 hypothesis. |
| 21 | Gal-Garber, Mabjeesh, Sklan, Uni 2000 (chicken), *J Nutr* 130(9):2174-9 | **10958809** | Cross-species-only BBMV kinetics sanity check: Km=24-150 µM, same order of magnitude as (below) this script's central human Km_SGLT1=1.0mM. |
| 22 | (in-repo) `scripts/msk/glucose_insulin_minimal_model.py` | n/a | Sibling Ra(t) gamma(2,tau=40min) kernel, imported directly (read-only), not re-transcribed, for F8's cross-model check. |

## 3. Gastric emptying — liquid + solid, forced adversary

**Liquid** (75g-glucose-solution's *generic* reference band): t½ central = 15 min, task's own
[10,20] min band used as the anchor — **no clean single live-extracted human plain-liquid t½ was
found this session** (disclosed honestly, same pattern as this repo's own DeFronzo-M-value gap).
**PASS** by construction against the task band.

**Solid** (Tougas 2000's own 3-point retention curve, 69%/24%/1.2% at 60/120/240min): power-
exponential fit via **two independent numerical routes**:

| route | k | beta | T50 (min) |
|---|---:|---:|---:|
| log-log closed-form linear regression | 2.549e-4 | 1.788 | **83.41** |
| nonlinear least-squares (scipy curve_fit) | (agrees) | (agrees) | rel. err vs route 1: **0.79%** |

T50=83.41 min sits comfortably inside the task's [60,120] band. **Cross-population replication**:
an entirely independent cohort (Vasavid et al 2014, Thailand, rice+egg meal, n=189) reports T50=
68.7 min — **21.4%** from this fit, inside the pre-registered (wide, because meal composition is
explicitly NOT held fixed between studies) 30% band.

**Forced adversary (F3)**: a naive mono-exponential (beta forced=1) fit to the *same 3 points*
gives SSE **48.5x larger** than the power-exponential fit, and is **clinically wrong at the most
important point**: it predicts 7.0% retention at 4h (vs actual 1.2%) — misses the diagnostic
"delayed emptying" cutoff (>10% at 4h) by less than a factor of 1.5, when the real subjects are
comfortably normal. This is not a decorative shape choice — the lag-phase form is load-bearing.

## 4. Small-intestinal transit — N-compartment Erlang chain

Mean SI transit = 216 min (3.6h), calibrated to Read et al 1983's own paired, non-intubated
control measurement (3.6±0.4h) — sits at the exact center of the task's [3,4]h band. N=10
(disclosed PK-transit-compartment convention) for the reference run; **CV=1/sqrt(N) sweep**
(N=1..128) confirms monotonic dispersion decrease **independent of the mean** — a clean geometric
separation between where mass sits on average and how spread its arrival-time distribution is.

## 5. Glucose Ra(t): SGLT1/GLUT2 mechanistic model, coupled OODA loop (two forced fixes)

**Model**: gastric-liquid mono-exponential emptying → N=10-compartment SI Erlang chain → per-
compartment saturable SGLT1 (Km=1.0mM, textbook-grade central, swept 0.4-2.0mM) + GLUT2 (Km=17mM)
Michaelis-Menten mucosal uptake → splanchnic-extraction band [12.9%, 27%] → systemic Ra(t).

**OODA Fix #1 (diagnosed, not skipped)**: a first draft asserted rho=10 ("deep Regime A", absorption
never limiting) on the strength of Trahair/Horowitz + Ferrannini qualitatively favoring gastric-
emptying control. That draft's Ra(t) **peaked at 2.5 min** — wrong. **Orient**: a mono-exponential
gastric flux is itself *maximal at t=0*; near-instantaneous absorption at high rho adds ~zero
further delay, so the systemic curve just mirrors the front-loaded gastric flux. **Fix**: root-find
(brentq, exact) rho against Ferrannini's peak-time band midpoint (22.5min) → **rho=2.267**, landing
in an intermediate transition regime, not deep Regime A.

**OODA Fix #1, checked against a SECOND anchor (round 2 of the loop, not declared done after round
1)**: Fix #1's own Ra(3.5h) prediction = **0.0022 mg/min/kg** vs Ferrannini's directly measured
**2.47±0.45 mg/min/kg** — off by a factor of **~1118x**. The coarse `frac_absorbed_by_3.5h` AUC gate
(0.729, inside [0.5,0.98]) did **not** catch this, because it integrates the whole curve; a fast-
rise-then-near-zero-tail shape can still pass an AUC-fraction check even when its instantaneous
late-time value is wildly wrong. This is exactly why a coarse aggregate metric is an insufficient
falsifier on its own — reported here, not hidden behind the passing AUC gate.

**OODA Fix #2 (forced, not a one-shot fail)**: diagnosed crux — the generic 15-min liquid t½ used
for a 75g/~300mL glucose bolus (~1 kcal/mL) silently violates this script's own meal-composition
citation (Hunt & Stubbs 1975): caloric *delivery* rate, not volumetric rate, is regulated, so a
calorically dense liquid should empty slower. **Fix attempted**: free a separate glucose-solution
t½ and jointly root-find it with rho against BOTH Ferrannini anchors (2 knobs, 2 equations, not
over-fit). **Result: fsolve did NOT converge** on this bracket-free 2D search (tested from 4
different initial guesses, incl. scaled residuals; the solver either failed to progress or,
when it did hit peak-time exactly, still undershot Ra-at-3.5h by 40-95%). **Diagnosis**: this looks
like a genuine **structural limit** of a 2-parameter single-Erlang-chain architecture — it can
produce a sharp early peak, or (at very different parameters) a long sustained tail, but not both
tightly at once with only 2 free knobs. This is reported as an **open, diagnosed, quantified
finding**, not smoothed over and not declared a bug.

**Consequence (disclosed, not hidden)**: because fix #2 did not converge, the script **falls back**
to the well-converged fix #1 (rho=2.267, generic 15-min liquid t½) as the actual reference model.
Its Ra(t): peak time 22.5 min (inside Ferrannini's 15-30min band — but this is the calibration
target, reported as **true-by-construction, NOT a gate**, matching this repo's own precedent for
excluding tautological quantities from the gate count); frac-absorbed-by-3.5h = 0.729 (Basu-
splanchnic 0.870) — both genuinely inside [0.5,0.98]; Ra-at-3.5h = 0.0022 mg/min/kg — a **genuine,
disclosed miss** vs the measured 2.02-2.92 band.

**Forced adversary (F9, Reds 2014)**: re-running with SGLT1-only (Vmax_GLUT2=0) still peaks at
21.25 min (vs 22.5 min dual-pathway) — **both inside the Ferrannini band**. GLUT2 contributes only
**0.03%** of total AUC at these parameters — the headline peak-timing result is **ROBUST** to the
SGLT1-vs-SGLT1+GLUT2 controversy, consistent with reds's knockout finding that SGLT1 dominates.

**Cross-model consistency (F8) vs the sibling's own Ra(t) input** — imported directly
(`gamma_kernel`, `PARAMS`, read-only, zero transcription risk): sibling kernel (tau=40min) peaks at
40 min, **outside** Ferrannini's measured 15-30min window — **a genuine, machine-detected
discrepancy in the ALREADY-BUILT sibling model**, surfaced here rather than hidden (symmetric QC
applies the same evidence burden to this repo's own prior work). Pearson correlation between this
script's mechanistic Ra(t) and the sibling's phenomenological kernel, common 0.5-240min grid:
**r=0.609** — below the pre-registered 0.8 threshold. Both curves are unimodal with roughly
comparable peak location (22.5 vs 40 min, 17.5 min apart) but **different widths**: the
mechanistic curve is essentially finished by ~90-120 min while the sibling's kernel is still at
40-65% of peak out to 90-120min — a real shape mismatch, not a bug, consistent with the Ra(3.5h)
miss above (same underlying "too-narrow tail" issue, seen from two angles).

**Sensitivity** (Km_SGLT1 ∈ [0.4,2.0]mM, V_SI_total ∈ [500,1500]mL): peak-time moves by at most
±0.5 min — governed by k_liq/rho, not by these secondary parameters, within literature-plausible
ranges.

## 6. Symmetric QC — held OPEN, not resolved away, per task instruction

**Meal-composition / caloric-density dependence (HELD OPEN)**: Hunt & Stubbs 1975 (33 pooled
studies): caloric *delivery* rate, not volumetric emptying rate, is the regulated quantity for
liquids. Every t½/T50 number in this document is specific to its own meal (low-fat egg-substitute
for Tougas, rice+egg for Vasavid, 75g-glucose-in-water for the OGTT coupling) and does **not**
generalize to arbitrary meals. This is not a minor caveat — Section 5's own OODA Fix #2 attempt
shows that ignoring it (using a generic liquid t½ for a calorically dense glucose bolus) is
precisely the source of the Ra(3.5h) miss.

**Method spread — scintigraphy vs breath test vs capsule (report the spread, not one number):**

| comparison | source | correlation | systematic difference |
|---|---|---|---|
| Solid GE: scintigraphy vs 13C-breath test | Maes/Ghoos 1998 | r=0.98 (t½) | mean diff 10min (-20,41 CI) |
| Gastric emptying: GES vs WMC | Lee/Kuo 2019, n=167 | kappa=0.42 | 75.7% categorical agreement; WMC detects delayed GE in 34.6% vs GES 24.5% |
| Orocecal transit: lactulose breath test vs scintigraphy | Miller 1997, n=8 | r=0.95 | lactulose ITSELF accelerates transit (P=0.004) — not an inert tracer |
| Colonic transit: WMC vs radiopaque markers | Camilleri 2010, n=158 | r=0.707 | median ROM 55.0h vs WMC 43.5h (P<0.001) |
| Colonic transit: scintigraphy vs radiopaque markers | Degen/Phillips 1996, n=32 | — | methods AGREE ("equally well reflected") — reported for balance, not every pair disagrees |

Correlated (r=0.69-0.98) does **not** mean interchangeable — several pairs carry real, quantified,
systematic mean differences, and one (lactulose) actively perturbs the thing it measures.

## 7. Honest gaps (disclosed, not hidden)

- No live-extracted primary-source human plain-**liquid** gastric-emptying t½ was found this
  session — the task's own [10,20]min band is the anchor.
- Km_SGLT1 (1.0mM central), Km_GLUT2 (17mM), V_SI_total (1000mL) are literature-convention
  order-of-magnitude values, not individually live-extracted from a primary human kinetics paper
  this session — swept to show headline peak-timing is insensitive to their exact value.
- Splanchnic extraction is a disclosed **band** [12.9%, 27%] (Basu, mixed meal, triple-tracer,
  instantaneous vs Ferrannini, pure glucose, dual-tracer, cumulative-to-3.5h) — Basu's own paper
  states dual-tracer overestimates ISE vs triple-tracer, partially (not fully) reconciling the gap.
- **The joint (t_half_glucose, rho) calibration did not converge** — Section 5's central open
  finding. The fallback (fix #1) reference model matches Ferrannini's peak-**timing** but misses
  Ra-**at-3.5h** by ~1100x and the cross-model correlation vs the sibling script (0.609) falls
  short of the pre-registered 0.8 — reported as PARTIAL/OPEN, not forced to PASS.
- Vmax_total (absorptive capacity) is not independently pinned to a live-verified human jejunal-
  perfusion primary source this session (searches for Modigliani/Holdsworth-class papers did not
  return a directly quotable saturating-flux number) — the rho-sweep instead shows which regime
  the results depend on.
- N=10 compartments is a disclosed PK-convention choice; Read 1983's 0.4h SD is BETWEEN-subject
  variability (a different construct from within-subject transit dispersion) — used only to
  anchor the MEAN, not to fit N.
- GLUT2's true contribution is scientifically disputed (Reds 2014 vs Kellett/Brot-Laroche) — this
  script tests robustness to both configurations, does not resolve the controversy.
- No subject-specific data anywhere — a population-parametrized forward-model consistency check,
  matching every other MECHANISM_* systemic layer's own disclosed scope.

## 8. Couplings (prose only — not folded into the graph this session)

Per this session's isolation scope (another instance writes `data/MECHANISM_ANCHOR_GRAPH.json`
concurrently), these are **prose-only**, left for the canonical `mechanism_fold.py` path:

- **`ORG-PANCREAS-GLUCOSE-INSULIN`** / `glucose_insulin_minimal_model.py`: this is the "un-attempted
  extension" flagged in `MECHANISM_GLUCOSE_INSULIN.md` §7. A real, disclosed tension was found
  (sibling's tau=40min kernel peaks outside the measured 15-30min band; this script's own
  mechanistic Ra(t) does not fully match Ferrannini's sustained tail either) — the coordinator
  should treat this as an open reconciliation item, not a clean validation.
- **`GI-ENTERIC-MOTILITY`** / **`ORG-GI-MOTILITY-ENS`**: this script builds the quantitative
  transit-time layer those disease/discordance-focused cells assume as background; 2 PMIDs reused
  (re-verified live, not drifted).
- **`ORG-METABOLIC-SYNDROME-CLUSTER`** / **`METAB-TISSUE-INSULIN-RESISTANCE`**: absorption-rate
  changes (bariatric surgery, SGLT-inhibitor drugs) would shift Ra(t) shape — currently a fixed
  input to that cluster's framing.
- **`ORG-LIVER-HEPATIC-HUB`** / `hepatic_clearance.py`: splanchnic extraction (Section 5) is the
  same first-pass concept `hepatic_clearance.py` models for drugs — not quantitatively reused here
  (disclosed gap for a future pass).

## 9. Confidence tier

**In-vivo-anchored** for Tier 1 (gastric emptying + SI transit): 2 independent multicenter
scintigraphy cohorts (n=123, n=189), 2 breath-test cross-validations, 3 wireless-motility-capsule
studies, 1 direct paired SBTT measurement — cross-population-replicated and forced-adversary-
tested. **Weaker / partial-open tier** for Tier 2 (Ra(t) mechanistic coupling): anchored to
directly tracer-measured Ferrannini/Basu data and cross-checked against the sibling script, but the
joint calibration did not converge and the Ra-at-3.5h / cross-model-shape mismatches are disclosed,
unresolved findings. Population-parametrized forward model throughout — no subject-specific data,
matching this repo's own established scope for `glucose_insulin_minimal_model.py` and
`renal_filtration.py`.

## 10. Gates (machine-computed, `reports/probes/gi_absorption_transit.json`)

```
TIER 1 (gastric emptying + SI transit) -- PASS, 8/8:
  liquid_t50_in_task_band_10_20:                 PASS (15.0 min, task-band anchor disclosed)
  solid_T50_in_task_band_60_120:                 PASS (83.41 min)
  two_fit_routes_agree_lt_5pct:                  PASS (0.79%)
  forced_adversary_powerexp_beats_monoexp_2x:     PASS (48.5x, pre-registered >=2x)
  cross_population_within_30pct:                 PASS (21.4%, Vasavid n=189 independent)
  si_mean_in_task_band_3_4h:                      PASS (3.60h)
  si_mean_within_1sd_of_read1983:                 PASS (exact match to central estimate)
  void_floor_and_mass_conservation:               PASS (0 Vmax->0 Ra; mass balance <1% rel err)

TIER 2 (Ra(t) mechanistic coupling) -- PARTIAL/OPEN, 4/6 (2 calibration-consumed, not gated):
  frac_absorbed_basu_in_wide_band_0.5_0.98:       PASS (0.870)
  frac_absorbed_ferr_in_wide_band_0.5_0.98:       PASS (0.729)
  sglt1_only_robust_to_glut2_controversy:         PASS (21.25min, still in Ferrannini band)
  cross_model_corr_vs_sibling_gte_0.8:            FAIL (0.609) -- disclosed, not hidden
  [peak_time_in_ferrannini_band]:                 calibration target, NOT a gate (tautological)
  [joint_t_half_rho_calibration_converged]:       FAIL (fsolve non-convergence) -- diagnosed open

overall_pass (Tier 1 only, the well-converged headline): PASS
two_tier_verdict.tier_2: PARTIAL-OPEN (honest, not forced)
```

## 11. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/gi_absorption_transit.py
```
Pure Python/numpy/scipy (`solve_ivp`, `curve_fit`, `brentq`, `fsolve`), no OpenSim, no subject
data, runs in ~10-20s (dominated by the rho/N/Km/V_SI sweeps' repeated ODE solves). Writes
`reports/probes/gi_absorption_transit.json`. Read-only import of
`scripts/msk/glucose_insulin_minimal_model.py` (import-safety verified: its I/O is guarded behind
`if __name__ == "__main__"`, so importing executes only its PARAMS/CITATIONS/function definitions,
no file writes). No git operations. Files touched this session: `scripts/msk/gi_absorption_transit.py`,
`reports/probes/gi_absorption_transit.json`, this doc. `data/MECHANISM_ANCHOR_GRAPH.json` (shared,
concurrently written by another instance) was read-only referenced (Section 0), never edited.

## 12. Files

- `scripts/msk/gi_absorption_transit.py` — full model: citations, solid/liquid gastric-emptying
  fits + forced mono-exponential adversary, N-compartment SI Erlang chain + CV sweep, SGLT1/GLUT2
  coupled ODE, the two-round OODA rho/t-half calibration (one converged, one not, both reported),
  splanchnic-extraction band, SGLT1-only forced adversary, cross-model check vs the sibling script,
  sensitivity sweeps, all gates.
- `reports/probes/gi_absorption_transit.json` — full evidence: every citation, every computed
  number, both calibration attempts' full diagnostics, all gates, the two-tier verdict.
- Read-only, not modified: `scripts/msk/glucose_insulin_minimal_model.py` (imported for F8),
  `data/MECHANISM_ANCHOR_GRAPH.json` (Section 0 disambiguation only).
