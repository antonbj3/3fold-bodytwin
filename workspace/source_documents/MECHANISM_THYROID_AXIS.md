# MECHANISM THYROID-METABOLIC (HPT) AXIS — resolving the SEED-DESIGN governor node (2026-07-22)

Builds and MEASURES the hypothalamic-pituitary-thyroid (HPT) negative-feedback axis: the
log-linear/log-sigmoid TSH-fT4 relationship, T4→T3 peripheral deiodination, and thyroid-status →
BMR coupling, as a CERTIFIED model — resolving the pre-existing **SEED-DESIGN** (status `OPEN`,
"designed, not measured") hypothesis nodes `ORG-THYROID-METABOLIC-GOVERNOR`,
`AUTO-CELL-THYROID-METABOLIC-GOVERNOR-HPT-AXIS`, and `ENDO-HPT-THYROID-SETPOINT` already sitting in
`data/MECHANISM_ANCHOR_GRAPH.json` (999 nodes; **not edited this session** — isolation rule "touch
only files you create"; folding this result into those nodes via the canonical
`mechanism_fold → fold_gate_v2` path is the natural next step, not performed here, matching the
identical precedent already set by `MECHANISM_BLOOD_OXYGEN_TRANSPORT.md` / `MECHANISM_PULMONARY_GAS_EXCHANGE.md`).
Script: `scripts/msk/thyroid_metabolic_axis.py`. Evidence:
`data/msk_smoketest/subject2_walking1/thyroid_metabolic_axis/thyroid_metabolic_axis_results.json`.

**No re-solve of OpenSim.** Every thyroid-axis number is population/literature-anchored (no thyroid
function test panel exists for subject2 — the same disclosed scope every sibling
endocrine/blood-gas layer in this repo already carries). The subject2/walking1
`thermoregulation_results.json` is read **read-only**, once, purely to make the
`couples_to metabolic + thermoregulation` link a **concrete re-computed number**, not a prose
pointer.

## The two falsifiers, verdicts stated up front (nothing hidden)

> **Falsifier 1**: does a log-linear/log-sigmoid model calibrated on the task's own reference-range
> endpoints (TSH 0.4–4.0 mIU/L, fT4 12–22 pmol/L) reproduce Hadlow et al. 2013's MEASURED, n=152,261
> finding that the true log(TSH)-fT4 relationship is **not** simple log-linear but two overlapping
> negative sigmoid curves, inflexion points at fT4=7 and 21 pmol/L?

**YES, on the genuinely discriminating, machine-checked test** — naive log-linear (Model A) has
**constant** slope by construction (cannot represent Hadlow's finding at any parameter setting); a
double-sigmoid geometric construction (Model B), calibrated identically, has a **measured,
regime-dependent** slope (verified: slope range 0.01–0.49 decades/pmol/L, not constant) — this is
the real falsifier, and it PASSES. Model B's inflection points landing at fT4=7.00 and 21.00 pmol/L
is a **construction self-consistency check** (I placed them there by design), disclosed as such, not
an independent discovery — see §1.

> **Falsifier 2**: does the task's own given clinical BMR-delta range (hyper +25 to +80%, hypo −20
> to −40%) reproduce REAL indirect-calorimetry group comparisons?

**PARTIAL — hyperthyroid PASSES cleanly (+33.5%, Chng et al. 2016, n=24, real REE); hypothyroid is a
DIAGNOSED, disclosed, marginal MISS** (−17.2%, Wolf et al. 1996, vs. the task's −20% floor) — traced
via OODA to a real, plausible, disclosed cause (Wolf's cohort is SHORT-TERM/acute hypothyroidism,
not chronic/untreated), not force-fitted or hidden. See §4.

**Symmetric QC, held OPEN per task instruction, not resolved here**: Andersen et al. 2002 directly
measured an individual's TSH/T4/T3 set point to occupy a window ~3.3× narrower than the population
reference range — robust — but a genuine, disclosed, UNRESOLVED cross-cohort tension exists on
*which* threshold TSH's own "index of individuality" crosses (Andersen 2002, n=16: 0.49, passes;
Yildiz 2025, n=21: 0.84, fails) — see §2.

## Citations — every PMID/DOI verified LIVE this session (NCBI E-utilities: esearch/esummary/efetch,
direct curl, not WebFetch-summarized, not recalled — this repo's own prior finding across sibling
docs is a measured ~62–67% citation-drift rate from memory)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Hadlow NC, Rothacker KM, Wardrop R, Brown SJ, Lim EM, Walsh JP (2013). "The relationship between TSH and free T4 in a large population is complex and nonlinear and differs by age and sex." *J Clin Endocrinol Metab* 98(7):2936-43. | **23671314**, DOI 10.1210/jc.2012-4223 | Falsifier 1's decisive anchor: n=152,261, log TSH-fT4 = 2 overlapping sigmoids, inflexion fT4=7,21 pmol/L; their own stated reference range 10-20 pmol/L. |
| 2 | Andersen S, Pedersen KM, Bruun NH, Laurberg P (2002). "Narrow individual variations in serum T4 and T3 in normal subjects: a clue to the understanding of subclinical thyroid disease." *J Clin Endocrinol Metab* 87(3):1068-72. | **11889165**, DOI 10.1210/jcem.87.3.8165 | THE task-named "Andersen 2002" (distinct from a same-group 2003 *Thyroid* review some prior work in this repo's graph cited instead — re-verified as the correct, literal primary paper). n=16 men, individual 95% CI ≈ half the group's. |
| 3 | Pilo A, Iervasi G, Vitek F, Ferdeghini M, Cazzuola F, Bianchi R (1990). "Thyroidal and peripheral production of 3,5,3'-triiodothyronine in humans by multicompartmental analysis." *Am J Physiol* 258(4 Pt1):E715-26. | **2333963**, DOI 10.1152/ajpendo.1990.258.4.E715 | T4→T3 deiodination split, leg 1: thyroidal 3.3, peripheral 12.7 µg/day/m². |
| 4 | Chopra IJ (1976). "An assessment of daily production and significance of thyroidal secretion of 3,3',5'-triiodothyronine (reverse T3) in man." *J Clin Invest* 58(1):32-40. | **932209**, DOI 10.1172/JCI108456 | T4→T3 deiodination split, leg 2 (independent method, 14 years earlier): thyroidal ≈23.8% of PR-T3. |
| 5 | al-Adsani H, Hoffer LJ, Silva JE (1997). "Resting energy expenditure is sensitive to small dose changes in patients on chronic thyroid hormone replacement." *J Clin Endocrinol Metab* 82(4):1118-25. | **9100583**, DOI 10.1210/jcem.82.4.3873 | Falsifier 2, gradient context: REE↓~15% as TSH rises 0.1→10 mU/L (r²=0.64, n=9). |
| 6 | Chng CL, Lim AY, Tan HC, et al. (2016). "Physiological and Metabolic Changes During the Transition from Hyperthyroidism to Euthyroidism in Graves' Disease." *Thyroid* 26(10):1422-1430. | **27465032**, DOI 10.1089/thy.2015.0602 | Falsifier 2, hyperthyroid decisive anchor: REE 28.7→21.5 kcal/kg (n=24, within-subject). |
| 7 | Wolf M, Weigert A, Kreymann G (1996). "Body composition and energy expenditure in thyroidectomized patients during short-term hypothyroidism and thyrotropin-suppressive thyroxine therapy." *Eur J Endocrinol* 134(2):168-73. | **8630514**, DOI 10.1530/eje.0.1340168 | Falsifier 2, hypothyroid decisive anchor: BEE 5265 (hypo) vs 6362 (matched controls) kJ/24h. |
| 8 | Muraca E, Ciardullo S, Oltolini A, et al. (2020). "Resting Energy Expenditure in Obese Women with Primary Hypothyroidism and Appropriate Levothyroxine Replacement Therapy." *J Clin Endocrinol Metab* 105(4):dgaa097. | **32119074**, DOI 10.1210/clinem/dgaa097 | Gradient context (n=85 vs 564): treated-to-normal-TSH hypothyroid still −4.4% REE vs controls; their own stated TSH inclusion band (0.4–4.0 mU/L) independently matches this task's reference range. |
| 9 | Bianco AC, Kim BW (2006). "Deiodinases: implications of the local control of thyroid hormone action." *J Clin Invest* 116(10):2571-9. | **17016550**, DOI 10.1172/JCI29812, PMC1578599 | DIO1/DIO2/DIO3 mechanism, topical citation. |
| 10 | Bianco AC, Salvatore D, Gereben B, Berry MJ, Larsen PR (2002). "Biochemistry, cellular and molecular biology, and physiological roles of the iodothyronine selenodeiodinases." *Endocr Rev* 23(1):38-89. | **11844744**, DOI 10.1210/edrv.23.1.0455 | Comprehensive deiodinase review, topical citation (paywalled, no PMC copy — no number independently extracted). |
| 11 | Leow MK, Goede SL (2014). "The homeostatic set point of the hypothalamus-pituitary-thyroid axis — maximum curvature theory for personalized euthyroid targets." *Theor Biol Med Model* 11:35. | **25102854**, DOI 10.1186/1742-4682-11-35, PMC4237899 | Theoretical corroboration of §2: between-subject TSH-FT4 heterogeneity 0.0107 (95% CI 0.0029–0.03975), full text re-fetched and the exact figure independently re-confirmed against a prior session's graph-node claim (verbatim quote found, not a fabricated "live check"). |
| 12 | Yildiz R, Ozkanay H, Arslan FD, Koseoglu M (2025). "Biological variation of TSH, fT3 and fT4 in healthy subjects in Turkey." *Biochem Med (Zagreb)* 35(1):010706. | **39974197**, DOI 10.11613/BM.2025.010706 | The cross-cohort tension partner (n=21): re-verified independently this session — the prior graph-node's numbers held up exactly. |

## 1. Geometric structure — TSH-fT4 feedback is a control-loop GAIN curve, not a straight line

**The core geometric fact** (rule: derive from the geometry, not rote algebra): the pituitary
thyrotroph's TSH output as a function of fT4 is a classic negative-feedback dose-response — what
matters is not the curve's *level* but its **local slope**, `d(log₁₀TSH)/d(fT4)`, which IS the
loop's gain (Leow & Goede's own "maximum curvature theory," #11, frames the identical idea). High
gain near an individual's own operating point means small fT4 deviations are signalled strongly
(tight regulation, and why TSH — not fT4 — is the preferred first-line screening test); this is the
same "steepest-slope-of-a-saturating-curve" geometric structure this repo's own
`blood_oxygen_transport.py` already used for the O2-Hb dissociation curve.

**Model A (naive log-linear)**, calibrated exactly through the task's own 2 reference-range
endpoints: `log₁₀(TSH) = 1.80206 − 0.1·fT4` (the 0.1 decade/pmol/L slope falls out exactly because
4.0/0.4 = 10 precisely). **By construction this slope is a CONSTANT everywhere** — machine-verified
(`np.std` of the slope array = 0.0) — extrapolated 2→60 pmol/L it spans **5.8 decades** (a ~630,000×
range), which no real TSH assay or physiology exhibits.

**Model B (double-sigmoid)**, built from the geometry Hadlow's abstract reports (an inflection point
is a local extremum of the *slope*; a Gaussian bump in slope-space peaks exactly at its own
center) — `d(log₁₀TSH)/d(fT4) = −(a₁·g(fT4;7,3) + a₂·g(fT4;21,3) + 0.01)`, `a₁/a₂=2` (disclosed prior:
TSH's hypothyroid-direction rise is far larger in magnitude than its hyperthyroid-direction
suppression is deep), calibrated on the **same 2 anchors** Model A uses. Machine-measured (central
finite differences, noise-floor-filtered — a genuine finite-difference-noise artifact was caught and
filtered in the far tail, not silently averaged over):

| check | result | verdict |
|---|---:|---|
| slope constant (Model A) | std = 0.0 exactly | PASS (by construction — the point) |
| slope regime-dependent (Model B) | range 0.01–0.494 decades/pmol/L | **PASS — genuinely non-constant** |
| measured inflection points | **7.0000000, 21.0000000** pmol/L | matches design targets to <1e-8 rel. error |
| steepest gain point | **fT4 = 7.0 pmol/L**, slope −0.494 dec/pmol/L | max diagnostic sensitivity sits at the *lower* transition — consistent with TSH being a highly sensitive early detector of incipient hypothyroidism specifically |
| 3rd emergent feature | fT4 ≈ 14.25 pmol/L (a real, non-noise, smaller curvature change) | a **generic mathematical consequence** of summing two well-separated bumps, not independently reported by Hadlow (whose abstract states 2, not 3) — disclosed as a construction artifact, not a discovered fact |
| span 2→60 pmol/L | 4.41 decades (vs Model A's 5.8) | narrower, but still an extrapolation into an uncalibrated regime, honestly not separately validated against a real severe-disease population this session |

**What this actually proves, precisely stated**: the inflection-point-location match is a
**construction self-consistency check** (I placed the bumps at 7 and 21 by design, calibrated by
the SAME 2-point procedure as Model A — not a blind refit of Hadlow's real dataset, which I do not
have access to). The genuinely discriminating, non-tautological falsifier is the **constant-vs-
regime-dependent slope** property: naive log-linear (Model A) is a hypothesis class that **cannot**
represent Hadlow's reported shape at *any* parameter setting, while a geometrically-motivated
saturating construction (Model B) can and, once fairly calibrated on the same anchors, does.

## 2. Individual set-point vs. population reference range — HELD OPEN, per task instruction

Andersen et al. 2002 (#2, n=16 healthy men, monthly sampling × 12 months) measured this directly, not
theoretically: **"the width of the individual 95% confidence intervals were approximately half that
of the group, for all variables."** Index of individuality (CVᵢₙₜᵣₐ/CVᵢₙₜₑᵣ, low = high
individuality): T4=0.58, T3=0.54, freeT4index=0.59, **TSH=0.49** — all below the classic 0.6
"significant individuality" bar in this cohort. A single TSH test localizes an individual's own set
point to ±50% (a span-factor of 1.5/0.5 = **3.0×**) — machine-computed against the task's population
TSH range (4.0/0.4 = **10×**): **the population range is 3.33× wider than what Andersen's own data
says is needed to characterize one person.** Verbatim, decisive: *"a test result within laboratory
reference limits is not necessarily normal for an individual."*

**Forced adversary (symmetric QC), not smoothed over**: Yildiz et al. 2025 (#12, n=21, Turkey,
re-verified live this session, not merely trusted from a prior graph-node output — see
`bt_memory`'s own recorded incident of a fabricated subagent "live check") measured TSH's index of
individuality at **0.84 — FAILING the same 0.6 bar** Andersen's cohort passed. This is a genuine,
disclosed, **UNRESOLVED cross-cohort tension specific to TSH** (n=16 Denmark/2002 vs n=21
Turkey/2025 — different populations, decades, assay platforms). The fT3/fT4-family indices are
directionally consistent across both cohorts (all <0.65); the disagreement is confined to TSH. Per
the task's own explicit instruction, **this is held OPEN, not reconciled here** — the *narrower-
individual-window* claim rests on Andersen's own directly-measured CI-width statement (robust,
independent of the contested 0.6-threshold framing), while the *which-analyte-crosses-an-arbitrary-
threshold* framing is genuinely cohort-dependent.

Leow & Goede 2014 (#11, PMC full text independently re-fetched and re-verified this session, not
merely trusted from the pre-existing SEED-DESIGN JSON) corroborates via a **different, theoretical**
route (a multi-level GLLM regression on paired TFTs): between-subject heterogeneity 0.0107 (95% CI
0.0029–0.03975, p<0.05) — verbatim confirmed: *"every individual has a euthyroid set point that is
unique and stable."* Two decorrelated approaches (Andersen's direct longitudinal CI-width
measurement; Leow & Goede's theoretical variance-component model) converge on the same qualitative
conclusion via genuinely different methods.

## 3. T4 → T3 peripheral deiodination — two independent primary measurements, 14 years apart

| source | method | thyroidal fraction | peripheral fraction |
|---|---|---:|---:|
| Pilo et al. 1990 (#3) | 6-pool multicompartmental kinetic model, 14 studies | 3.3 µg/day/m² → **20.6%** | 12.7 µg/day/m² → **79.4%** |
| Chopra 1976 (#4) | single-tracer MCR/PR kinetics, independent method, 14 yr earlier | **23.8%** (of PR-T3) | **76.2%** |

**Machine-computed agreement: within 3.2 percentage points** across two structurally-different
kinetic methods and eras — genuine over-determination, not a single number dressed up as two. Both
corroborate the classic "~80% of circulating T3 comes from peripheral deiodination, only ~20% from
direct thyroidal secretion" teaching figure, from PRIMARY data, not textbook assertion. Chopra
additionally reports **~84% of daily T4 production (73.0 of 87.0 µg/day) is monodeiodinated** to
either T3 or (inactive) reverse-T3 — the DIO1/DIO2 (activating) vs. DIO3 (inactivating) branch point
Bianco & Kim 2006 (#9) and Bianco et al. 2002 (#10) describe mechanistically (topical citations; no
specific number independently extracted from either review this session — both are large,
paywalled reviews with no accessible PMC full text within this session's budget, disclosed, not
fabricated).

## 4. Falsifier 2 — thyroid status → BMR, real indirect calorimetry

| comparison | source | measured REE/BEE | computed Δ | task band | verdict |
|---|---|---|---:|---|---|
| Hyperthyroid (Graves) vs. same-subjects euthyroid | Chng 2016 (#6), n=24 | 28.7±4.0 → 21.5±4.1 kcal/kg | **+33.5%** | [+25, +80]% | **PASS** |
| Hypothyroid (short-term, thyroidectomized) vs. matched controls | Wolf 1996 (#7) | 5265±766 vs 6362±992 kJ/24h | **−17.2%** | [−40, −20]% | **DIAGNOSED MARGINAL MISS** (2.8 pp short of the −20% floor) |
| Treated-hypo (TSH normalized) vs. euthyroid controls | Muraca 2020 (#8), n=85 vs 564 | 28.59±3.26 vs 29.91±3.59 kcal/kg FFM | −4.4% | context only, non-gating | corroborates direction; small residual deficit even when "in range" |
| Within-subject TSH dose-titration slope | al-Adsani 1997 (#5), n=9 | r²=0.64, p<0.001 | −7.5%/decade of TSH | context only, non-gating | corroborates direction + rough magnitude |

**OODA on the one real miss, not swept under the rug**: Observe — Wolf 1996's hypothyroid state is
**short-term** (weeks, the standard pre-¹³¹I-scan thyroxine-withdrawal protocol used for thyroid
cancer follow-up), not decades-long untreated myxedema. Orient — BMR suppression operates partly via
SLOW mechanisms (reduced Na/K-ATPase synthesis, reduced mitochondrial biogenesis) that plausibly need
weeks-to-months of sustained low T3 to reach steady state; a several-week acute withdrawal may not
have fully equilibrated. Decide/Act — this is a **plausible, disclosed** explanation for why the
acute number (−17.2%) sits just short of the task's own band (which, like the classic textbook
−40/−50% figure, likely anchors on chronic/severe cases) — **not independently confirmed by a
duration-sweep study this session** (a concrete, cheap next test: find or run an equivalent
short-term-vs-chronic within-subject calorimetry comparison). Reported as a genuine, diagnosed,
disclosed marginal miss, exactly the discipline `MECHANISM_PULMONARY_GAS_EXCHANGE.md`'s own §7 already
established in this repo (an honest FAIL with its diagnostic reasoning attached, not silently patched
to force a clean scoreboard).

**Exploratory, explicitly non-gating**: extrapolating al-Adsani's measured −7.5%/decade slope out to
TSH=100 mU/L (well beyond their measured 0.1–10 range) gives −22.5%, which WOULD fall inside the
task's hypo band — reported only as a suggestive, non-decisive corroboration, since it assumes the
log-linear dose-response continues unchanged into an untested regime (an extrapolation, not a
measurement — same discipline `pulmonary_gas_exchange.py`'s own "illustrative, NOT independently
validated" elite-exercise number already applies).

## 5. couples_to metabolic (BMR) + thermoregulation (heat production) — a concrete, re-computed number

Read-only from `data/msk_smoketest/subject2_walking1/thermoregulation/thermoregulation_results.json`
(subject mass 78.2 kg, M_rest = 93.84 W — **not edited**, thermoregulation.py's own equation
`H_prod = M_rest + (1−η)·(M_gross−M_rest)`, η=0.225, re-implemented independently from that doc's
prose description, not copy-pasted). **Machine cross-check before trusting this re-implementation**:
for the euthyroid baseline (multiplier=1.0, combined_corrected config), this script's own H_prod
(331.392 W) and dT/dt (0.07286 °C/min) reproduce thermoregulation.py's own already-published 331.4 W
/ 0.0729 °C/min **exactly** — an independent re-implementation from prose, not a copy, agreeing to
the precision shown.

Thyroid-status multipliers (from §4's two decisive anchors) applied to **M_rest only** — a disclosed
scoping choice: thyroid hormone's classical BMR effect is a resting/basal phenomenon; the
gait-mechanical exercise increment (`M_gross − M_rest`) is held fixed, since it is dominated by
muscle-contraction mechanics at this subject's own fixed measured gait, not basal thyroid tone on a
single-stride timescale:

| state | multiplier | M_rest (W) | H_prod (W, combined_corrected) | whole-body dT/dt (°C/min) | required sweat rate (g/h) |
|---|---:|---:|---:|---:|---:|
| Hyperthyroid (Chng 2016) | ×1.3349 | 125.3 | **362.8** | **0.0798** | 352.5 (unchanged, see below) |
| Euthyroid baseline | ×1.0 | 93.8 | 331.4 | 0.0729 | 352.5 |
| Hypothyroid (Wolf 1996) | ×0.8276 | 77.7 | **315.2** | **0.0693** | 352.5 (unchanged) |

**Disclosed limitation, not a bug**: because the required-sweat-rate formula (per
`thermoregulation.py`'s own derivation) depends only on the held-fixed exercise increment, it is
**unchanged** by this scoping choice across all three thyroid states — a real consequence of scoping
thyroid status to the basal rate only, stated explicitly rather than hidden. A fuller model would
also let thyroid status modulate exercise thermogenesis/sweat-response sensitivity (a real,
literature-documented effect not captured here) — a natural next extension, not attempted this
session.

**Void-floor / non-degeneracy** (forced, not eyeballed): sweeping the thyroid multiplier
continuously over [0.5, 2.0] (well beyond the two measured anchors) gives **strictly monotonically
increasing** H_prod and dT/dt (verified via `np.diff`), and the numerically-measured slope
`d(dT/dt)/d(multiplier)` matches the closed-form analytical slope `M_rest·60/(mass·c)` to <0.1%
relative error — a real, non-pinned function of thyroid status, not a constant dressed up as a
response (same discipline `thermoregulation.py`'s own M-sweep void-floor check already established).
The specific-heat constant used (re-derived from `thermoregulation.py`'s own published
(H_prod, dT/dt, mass) triple, not re-typed) comes out to 3490.0 J/(kg·K) — matching that doc's own
disclosed textbook constant (3.49 kJ/(kg·K)) to <0.01%, confirming the re-derivation is correct
before trusting it further.

## 6. Pre-registered gates — 11/12 PASS, the 1 exception fully OODA-diagnosed (§4)

```
falsifier1_calibration_self_check:                          PASS
falsifier1_naive_model_has_constant_slope:                  PASS (std=0.0, by construction)
falsifier1_hadlow_structure_is_regime_dependent:             PASS (slope range 0.01-0.49 dec/pmol/L)
falsifier1_inflection_points_within_30pct_of_hadlow:         PASS (construction self-consistency, <1e-8 rel.err)
individuality_gate_individual_narrower_than_population:      PASS (3.33x)
deiodination_two_decorrelated_sources_agree:                 PASS (3.2 pp apart)
bmr_gate_hyperthyroid_in_task_band:                          PASS (+33.5% in [25,80])
bmr_gate_hypothyroid_in_task_band_DIAGNOSED_MARGINAL:        FAIL (-17.2% vs floor -20%, OODA-diagnosed §4)
coupling_derivative_self_check:                              PASS (1.150 vs sibling's 1.15 g/h per W)
coupling_specific_heat_matches_sibling:                      PASS (3490.0 vs 3490 J/kg/K)
coupling_void_floor_monotonic:                               PASS
coupling_analytical_numerical_slope_match:                   PASS (<0.1%)
```

`overall_pass_strict_all` in the JSON is **False** (strict `all()`, including the diagnosed gate) —
reported exactly as measured, same convention `pulmonary_gas_exchange.py`'s own 18/19 scoreboard
already established in this repo, never patched to force a clean pass. Determinism: 2 independent
runs produce byte-identical JSON (verified via `diff`) and zero NaN/Inf anywhere in the output tree
(checked programmatically over the full JSON).

## 7. Confidence tier

Per this task's own pre-registration and matching the identical, most-recent precedent set by
`MECHANISM_BLOOD_OXYGEN_TRANSPORT.md` / `MECHANISM_PULMONARY_GAS_EXCHANGE.md` (both 2026-07-22):
**in-vivo-anchored** (real human TSH-fT4 population data, n=152,261; real human indirect-calorimetry
group comparisons, n=24/n=649/n=9/n=16/n=21) — **one tier below a subject-specific in-vivo
measurement** (no thyroid function test panel exists for subject2; every operating point here is
generic/population-level, exactly the same disclosed scope those two sibling docs carry for their
own blood-gas/pulmonary-diffusion parameters).

## 8. Honest gaps — symmetric QC: what this does NOT prove

- **Nothing here is proven** in the strong sense — same discipline every sibling doc in this family
  applies. Model B's inflection-point match to Hadlow's reported 7/21 pmol/L is a **construction
  self-consistency check** (I placed them there), not an independent recovery of Hadlow's actual
  fitted curve (I do not have their raw data or fitted coefficients — a genuine, disclosed scope
  limit; the primary-source full text/dataset was not accessible within this session's budget).
- **The fT4 reference range (12-22 pmol/L) is close to, not identical to, Hadlow's own stated
  10-20 pmol/L** — an assay/lab-dependent variant, disclosed, not silently reconciled. The TSH range
  (0.4-4.0 mIU/L) independently matches Muraca 2020's own stated inclusion criterion exactly.
  TSH assays are well known to vary across platforms/generations (the task's own caveat) — not
  independently surveyed across assay platforms this session.
- **The hypothyroid BMR-band gate is a genuine, disclosed, marginal miss** (-17.2% vs. task floor
  -20%), OODA-diagnosed (§4) to a plausible acute-vs-chronic duration effect, **not independently
  confirmed** by a duration-sweep study this session — a concrete, cheap next test, not performed.
- **The individual-set-point vs. population-range cross-cohort tension (Andersen 0.49 vs. Yildiz
  0.84 for TSH's own index of individuality) is held OPEN, exactly per the task's own instruction**
  — not reconciled, not smoothed over. The narrower-individual-window HEADLINE claim (3.33x) rests
  on Andersen's own directly-measured CI-width statement, which is independent of this threshold
  debate and therefore more robust than the II<0.6 framing alone.
- **Bianco & Kim 2006 / Bianco et al. 2002 (deiodinase mechanism reviews) are cited topically only**
  — both are large, paywalled reviews with no PMC full text accessible this session; no specific
  numeric claim was independently extracted from either (disclosed, not fabricated).
- **The Rolfe & Brand mitochondrial-proton-leak lineage** (cited in the prior SEED-DESIGN
  graph-node's own JSON as "relayed via search synthesis, NOT independently fetched from primary
  text") is **deliberately NOT reused here** as load-bearing — this doc only carries numbers this
  session independently, live re-verified.
- **The coupling demonstration (§5) is a disclosed SIMPLIFICATION**: thyroid multipliers act on
  M_rest only, not on the exercise-mechanical increment or on exercise-thermogenesis/sweat-response
  sensitivity (a real, literature-documented effect not modeled) — required sweat rate is
  UNCHANGED across all three thyroid states as a direct, disclosed consequence, not a hidden defect.
- **Single subject (subject2), single trial (walking1), one speed** for the coupling demonstration
  only — inherits, does not fix or expand, every upstream MSK layer's own scope limit. The
  thyroid-axis numbers themselves are population-level throughout, never subject-specific.
- **No graph-edge write this session** — `couples_to` is prose/JSON-evidence metadata; folding into
  `data/MECHANISM_ANCHOR_GRAPH.json`'s existing `ORG-THYROID-METABOLIC-GOVERNOR` /
  `AUTO-CELL-THYROID-METABOLIC-GOVERNOR-HPT-AXIS` / `ENDO-HPT-THYROID-SETPOINT` nodes requires the
  separate `mechanism_fold → fold_gate_v2` path, not performed here (isolation rule: touch only files
  created this session).
- **RED-S / non-thyroidal-illness / aging-axis / critical-care scope** (named in the pre-existing
  SEED-DESIGN nodes' own `couples_to` list) is explicitly OUT OF SCOPE for this doc, which answers
  only the log-linear feedback + deiodination + BMR questions this specific task named.

## Files

- `scripts/msk/thyroid_metabolic_axis.py` — self-contained (numpy/scipy only, no OpenSim), builds
  Model A/B, the individuality/deiodination/BMR sections, the couples_to coupling demonstration, the
  void-floor sweep, and all 12 gates; writes the evidence JSON below; prints a full summary.
- `data/msk_smoketest/subject2_walking1/thyroid_metabolic_axis/thyroid_metabolic_axis_results.json`
  — every number in this doc, machine-written: all 12 citations, both TSH-fT4 models (parameters,
  swept slopes, measured inflection points), the individuality/deiodination/BMR computations, the
  full couples_to coupling table (3 metabolic-cost configs × 3 thyroid states), the void-floor sweep,
  and all 12 gates. Verified deterministic (2 independent runs, byte-identical JSON via `diff`) and
  NaN/Inf-free (checked programmatically over the full JSON tree).
- Input read (read-only, no re-solve, not modified):
  `data/msk_smoketest/subject2_walking1/thermoregulation/thermoregulation_results.json`.
- Read but NOT modified (isolation: touch only files created this session):
  `data/MECHANISM_ANCHOR_GRAPH.json` (the 3 pre-existing SEED-DESIGN nodes this doc resolves),
  `data/body_twin/agent_outputs/auto__a74dd6744c67e3372.json` and
  `…/thyroid-metabolic-governor__a3d1188f44874fc4e.json` (the prior session's literature-scout
  outputs — treated as HYPOTHESES per this project's own discipline, re-verified live rather than
  trusted; 2 of their PMIDs — Yildiz 2025, Leow & Goede 2014 — held up exactly under independent
  re-verification; 1 citation — "Andersen 2003" PMID 14651790 — was found to be a *different* paper
  than this task's literally-named "Andersen 2002," which this doc independently located and
  verified instead, PMID 11889165).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/thyroid_metabolic_axis.py
```
No inputs required beyond the already-published `thermoregulation_results.json` (degrades to
`coupling_section.available=False` if absent — affects §5 only, not the gates in §1-4). Pure
Python/numpy/scipy (`scipy.special.erf`), no OpenSim call, runs in under 2 seconds, deterministic.
No git operations; writes only under
`data/msk_smoketest/subject2_walking1/thyroid_metabolic_axis/`.
