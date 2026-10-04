# MECHANISM LYMPHATIC RETURN — closes the microcirculation fluid loop (2026-07-22)

Closes the loop `docs/MECHANISM_CAPILLARY_STARLING.md` opens: that doc does **not exist as a
doc yet** in this repo (checked live via `find -iname` — only its script + evidence JSON exist,
`scripts/msk/capillary_starling.py` / `data/capillary_starling/capillary_starling_results.json`,
`overall_pass: true`), so this document reads that JSON directly, the same way
`docs/MECHANISM_FLUID_COMPARTMENTS.md` already reads `data/renal_filtration/` and
`data/cardiac_output_geometric/` read-only. Builds a certified model of **whole-body lymph
flow**, the **lymphatic pump** (intrinsic contraction + extrinsic muscle/respiratory pump), the
**interstitial-pressure-vs-lymph-flow relationship** (Guyton/Aukland), and **protein return**.
Script: `scripts/msk/lymphatic_return.py`. Evidence: `data/lymphatic_return/lymphatic_return_results.json`.

## 0. The adversary this document is built not to dodge

`capillary_starling.py`'s own `open_modeling_uncertainty` discloses that its `Kf_total` (1.205
L/day/mmHg) was **calibrated from** the Renkin (1986) 8 L/day anchor, not independently measured
(Kf/sigma are tissue-specific and hard to measure). This session independently re-fetched Bhave &
Neilson (2011)'s full text (PMC4096826 — not copy-trusted from either sibling doc's citation
table) and found their 10-g-albumin/hour protein-flux figure (used below in STEP 4) traces to
**the same Renkin (1986) review**, its reference 82. So a mass-balance "closure" test — does
modeled lymph flow match capillary filtration? — that only checks against Renkin's own number
would be a **tautology** (one external source cited on both sides), not the over-determination
the task's falsifier asks for. STEP 1 below machine-flags this shared common-mode explicitly
(never laundered as independent) and builds a **genuinely decorrelated** leg — Mobley et al.
(1989)'s direct dog thoracic-duct-cannulation flow rate, allometrically scaled, sharing no
citation with Renkin — to give the falsifier an actual chance to fail.

**Two real bugs were caught and fixed during this build (OODA, not swept under the rug — see
`scripts/msk/lymphatic_return.py` STEP 3 comments for the full diagnosis):** (1) a first attempt
at the Pi-vs-flow recruitment curve placed its inflection point at +0.5 mmHg, putting the
steepest part of the curve INSIDE the "plateau" sampling window instead of the "rising" one —
caught because the machine-computed steepness ratio came back 0.59 (backwards: the plateau region
was *steeper* than the rise). (2) after fixing that, a self-consistency bisection was written to
solve for `q_max` directly, which has a unique root independent of the search bracket — silently
collapsing an intended 3-point ceiling-multiple sweep into three byte-identical rows (all still
"passing", for a degenerate reason). Both are disclosed in the script and below, not hidden.

## 1. Citations — PMID/DOI verified LIVE this session (17 sources, tier stated per number)

Every PMID/DOI below was verified via raw NCBI eutils JSON (esearch existence → esummary
title/journal/author identity → efetch abstract/full-text content) plus one independent PMC
full-text refetch of Bhave & Neilson (2011). One recall-then-correction this session: Aukland &
Reed (1993)'s PMID was recalled as 8419966 — live search found the real PMID is **8419962** —
consistent with this repo's own previously-measured ~62–75% recalled-PMID drift rate.

| # | Citation | PMID / DOI | Tier / role |
|---|---|---|---|
| 1 | Renkin EM (1986). "Some consequences of capillary permeability to macromolecules: Starling's hypothesis reconsidered." *Am J Physiol* 250(5 Pt 2):H706-10. | **3706547**, `10.1152/ajpheart.1986.250.5.H706` | **REUSED, shared common-mode** with `capillary_starling.py` (§0). Quoted via Bhave & Neilson (2011): 8 L/day capillary filtrate → initial lymphatics, ~4 L/day reabsorbed in lymph nodes, remainder ~4 L/day via thoracic duct. Also the traced source of the 10 g-albumin/hour protein flux (STEP 4). |
| 2 | Bhave G, Neilson EG (2011). "Body fluid dynamics: back to the future." *J Am Soc Nephrol* 22(12):2166-81. | **22034644**, PMC4096826 | **Full text independently re-fetched** this session (not copy-trusted). Source of: Fig.6 compliance-curve quote, the `Jv=Q_L` steady-state identity quote, the 10g/hr albumin quote, and the Pi=-4-to-0mmHg figure. |
| 3 | Aukland K, Reed RK (1993). "Interstitial-lymphatic mechanisms in the control of extracellular fluid volume." *Physiol Rev* 73(1):1-78. | **8419962** (PMID recall-corrected live this session — 8419966 was wrong), `10.1152/physrev.1993.73.1.1` | THE definitive modern review. Abstract fetched live verbatim: Pi consensus "0 to -4 mmHg" in soft tissues; confirms Pi-volume curves have been recorded. **No PMC full text** (elink returns only a citing-articles list) — the primary Pi-vs-lymph-flow curve numbers were not extracted from body text this session, disclosed. |
| 4 | Guyton AC (1963). "A concept of negative interstitial pressure based on pressures in implanted perforated capsules." *Circ Res* 12:399-414. | **13951514**, `10.1161/01.res.12.4.399` | Originating METHOD paper. No abstract indexed (pre-1975); identity confirmed live; present in Bhave & Neilson's own reference list (confirmed live). |
| 5 | Guyton AC, Granger HJ, Taylor AE (1971). "Interstitial fluid pressure." *Physiol Rev* 51(3):527-63. | **4950077**, `10.1152/physrev.1971.51.3.527` | THE classic definitive review. No abstract indexed; identity confirmed live; no PMC full text — curve not independently extracted, disclosed. |
| 6 | Olszewski WL, Engeset A (1980). "Intrinsic contractility of prenodal lymph vessels and lymph flow in human leg." *Am J Physiol* 239(6):H775-83. | **7446752**, `10.1152/ajpheart.1980.239.6.h775` | **DIRECT HUMAN IN-VIVO** measurement — exactly the confidence tier this task asked for. Identity confirmed live (NCBI + EuropePMC, cited-by=137). No abstract accessible either route — cited for existence only, no number asserted from it. |
| 7 | Scallan JP, Zawieja SD, Castorena-Gonzalez JA, Davis MJ (2016). "Lymphatic pumping: mechanics, mechanisms and malfunction." *J Physiol* 594(20):5749-68. | **27219461**, PMC5063934 | Review, abstract fetched live verbatim: pumping regulated by "preload, afterload, spontaneous contraction rate, contractility" — the framework behind STEP 2. |
| 8 | von der Weid PY, Zawieja DC (2004). "Lymphatic smooth muscle: the motor unit of lymph drainage." *Int J Biochem Cell Biol* 36(7):1147-53. | **15109561**, `10.1016/j.biocel.2003.12.008` | Review, verbatim: intrinsic contraction "represents the principal mechanism by which lymph flow is generated." |
| 9 | Davis MJ, Scallan JP, Wolpers JH, Muthuchamy M, Gashev AA, Zawieja DC (2012). "Intrinsic increase in lymphangion muscle contractility in response to elevated afterload." *Am J Physiol Heart Circ Physiol* 303(7):H795-808. | **22886407**, PMC3469705 | **Direct quantitative measurement** (isolated rat mesenteric lymphangion, cannulated ex vivo). Verbatim: afterload ramp → "30% leftward shift in the end-systolic P-V relationship accompanied an 84% increase in dP/dt"; "weaker pumps... pump failure." |
| 10 | Majgaard J, Skov FG, Kim S, Hjortdal VE, Boedtkjer DMB (2022). "Positive chronotropic action of HCN channel antagonism in human collecting lymphatic vessels." *Physiol Rep* 10(16):e15401. | **35980021**, PMC9387113 | **Direct human ex vivo** (thoracic duct/mesenteric segments, surgical patients). Full text fetched live: "baseline CF was 4.4 ± 0.5 min⁻¹." Also discloses HCN channels are NOT the primary chronotropic driver in human vessels (open question vs rodent). |
| 11 | Havas E, Parviainen T, Vuorela J, Toivanen J, Nikula T, Vihko V (1997). "Lymph flow dynamics in exercising human skeletal muscle as detected by scintography." *J Physiol* 504(Pt 1):233-9. | **9350633**, PMC1159951 | **Direct human in-vivo** (99mTc-HSA scintigraphy, n=16). Rest 0.04±0.05 %/min; dynamic/isometric-ext/isometric-flex 0.16/0.20/0.09 %/min; abstract's own words "three- to sixfold." **Corrects** this repo's own prior flagged-uncertain "3-10x" claim (`data/body_twin/agent_outputs/lymphatic-immune-trafficking__a822249e9b42df4ac.json`) to the primary-verified figure. |
| 12 | Mobley WP, Kintner K, Witte CL, Witte MH (1989). "Contribution of the liver to thoracic duct lymph flow in a motionless subject." *Lymphology* 22(2):81-4. | **2770355** | **Direct animal (dog) in-vivo cannulation** — the genuinely independent leg (shares no citation with Renkin). Resting TDL flow 0.60±0.17 mL/min, TDL protein 3.4±0.5 g/dL; liver ≈1/3 of resting TDL flow; TDL/plasma-protein ratio 0.58→0.48 when hepatic inflow clamped. |
| 13 | Moman RN, Gupta N, Singh C, Varacallo MA. "Physiology, Albumin." StatPearls, 2026 Jan. | **29083605** | Verbatim: albumin "3.5-5 g/dL," hepatic synthesis "10 g to 15 g per day" — a DIFFERENT quantity from the transcapillary flux (STEP 4 forces these apart). |
| 14 | Parving HH, Hansen JM, Nielsen SL, Rossing N, Munck O, Lassen NA (1979). "Mechanisms of edema formation in myxedema." *N Engl J Med* 301(9):460-5. | **460364**, `10.1056/NEJM197908303010902` | **Direct human** (¹³¹I-albumin turnover, n=7). Verbatim: hypothyroidism → increased TER-albumin + "inadequate lymphatic drainage may explain" edema — real disease-contrast confirming TER-albumin is measurable and that lymphatic-failure edema is a distinct mechanism. Normal-subject TER%/hr baseline not stated in the abstract — disclosed gap. |
| 15 | Kramer GC, Harms BA, Gunther RA, Renkin EM, Demling RH (1981). "The effects of hypoproteinemia on blood-to-lymph fluid transport in sheep lung." *Circ Res* 49(5):1173-80. | **7296783**, `10.1161/01.res.49.5.1173` | Direct animal (sheep), **different organ/mechanism** (lung, oncotic not pressure-driven): "lymph flows increased to a maximum of 4 times baseline." Corroborating order-of-magnitude only, explicitly flagged as a different mechanism. |
| 16 | Venugopal AM, Stewart RH, Laine GA, Dongaonkar RM, Quick CM (2007). "Lymphangion coordination minimally affects mean flow in lymphatic vessels." *Am J Physiol Heart Circ Physiol* 293(2):H1183-9. | **17468331**, `10.1152/ajpheart.01340.2006` | Computational + bovine-mesenteric-validated. Verbatim: "Coordination of contraction had little impact on mean flow." Independent confirmation of STEP 2's conservation-in-series prediction for equal-capacity segments. |
| 17 | Bertram CD, Macaskill C, Davis MJ, Moore JE Jr (2017). "Valve-related modes of pump failure in collecting lymphatics." *Biomech Model Mechanobiol* 16(6):1987-2003. | **28699120**, PMC5671905 | Ex vivo + numerical. Verbatim: ramping afterload on a single lymphangion "until pumping fails"; "seven different possible modes of pump failure." The opposite-regime complement to Venugopal (2007). |

## 2. STEP 1 — Whole-body lymph flow: mass-balance closure, forced honestly

**Machine-checked**: `CITATIONS["renkin_1986"]["pmid"] == capillary_starling_json["citations"]["renkin_1986"]["pmid"] == "3706547"` → **True**. Both threads share the identical external source for
the 8/4/4 L/day figures — confirmed by direct field comparison, not narration.

**Genuinely independent leg** (Mobley et al. 1989, dog thoracic-duct cannulation, motionless/resting):
0.60 mL/min = **0.864 L/day** (dog). Allometrically scaled to this repo's own 73 kg reference man
(`fluid_compartments.py`), swept over dog mass [15,20,25,30] kg × scaling exponent [0.67
(surface-area), 0.75 (metabolic)] — 8 combinations, no cherry-picking:

| | Value |
|---|---:|
| Scaled human estimate, min | **1.57 L/day** |
| Scaled human estimate, median | **1.99 L/day** |
| Scaled human estimate, max | **2.83 L/day** |
| Task's own pre-registered band | [2.0, 4.0] L/day |

**A resting/motionless condition (no extrinsic muscle pump) lands at or just below the lower edge
of the ambulatory-inclusive Renkin-derived band** — the physiologically-expected direction
(rest ≤ active), from a completely different species, method (surgical cannulation vs
tracer/review synthesis), and condition, sharing zero citations with Renkin (1986). This is the
actual over-determination the falsifier asked for — not a second look at the same number.

**Void-floor forced adversary** — if lymphatic return were zero, filtrate (8 L/day, Renkin's total
capillary-filtrate figure, the correct "inflow" per Bhave's own `Jv=Q_L` identity) accumulates in
the interstitium at its initial rate. Using this repo's own already-certified interstitial volume
(`fluid_compartments.py`, reference man, 10.03 L) and Bhave & Neilson's own quoted edema-onset
threshold (Pi saturates once IFV is 20-50% above euvolemia):

| Threshold | Time (static, first-order) |
|---|---:|
| +20% interstitial volume (edema onset begins) | **6.0 hours** |
| +50% interstitial volume | **15.0 hours** |
| Double interstitial volume | **30.1 hours** |

A real, disclosed simplification (Pi-feedback on filtration is not modeled — this bounds an order
of magnitude, not a time-domain edema trajectory) — but it quantifies, rather than asserts, why
lymphatic return is necessary: without it, clinically-apparent edema is a matter of **hours**, not
weeks.

## 3. STEP 2 — The lymphatic pump: intrinsic + extrinsic, with a bottleneck governor

**Intrinsic**: baseline spontaneous contraction frequency in **human** collecting lymphatics
(thoracic duct, ex vivo) = **4.4 ± 0.5 min⁻¹** (Majgaard et al. 2022). Afterload response
(Davis et al. 2012, isolated rat lymphangion): elevated afterload → **84% increase in dP/dt**,
30% leftward P-V shift — a homeometric-autoregulation analog to the cardiac Anrep effect — until a
vessel-specific ceiling beyond which "weaker pumps" fail (their own finding, not asserted here).

**Extrinsic** (Havas et al. 1997, human scintigraphy, n=16): rest clearance 0.04%/min; exercise
raises it to 0.09-0.20%/min depending on contraction mode. Machine-recomputed multipliers vs rest
(not just repeating the abstract's rounded prose):

| Contraction mode | Clearance (%/min) | Multiplier vs rest |
|---|---:|---:|
| Rest | 0.04 | 1.0x |
| Isometric flexion | 0.09 | **2.25x** |
| Dynamic knee extension | 0.16 | **4.00x** |
| Isometric extension | 0.20 | **5.00x** |

Recomputed range [2.25x, 5.0x] **overlaps** (does not exactly equal) the abstract's own stated
"three- to sixfold" — reported precisely, not silently rounded to match.

**Series-chain bottleneck argument (geometric, derived not curve-fit)**: N lymphangions in a
valved series chain, at periodic steady state, must carry the SAME mean flow at every
cross-section (mass conservation) — chain throughput is therefore governed by the
**weakest/lowest-capacity segment (min), not the average**. Illustrative machine check: segment
capacities [1.0, 1.0, 0.4, 1.0, 1.0] → chain throughput = min = **0.4**, well below the mean
(**0.84**) — confirmed. Two decorrelated, opposite-looking real papers bracket this exactly:
**Venugopal et al. (2007)** (coordination among roughly-equal-capacity segments "had little impact
on mean flow" — no dominant bottleneck when segments are alike) and **Bertram et al. (2017)**
(ramping afterload on a SINGLE lymphangion drives multi-modal pump failure — a bottleneck matters
when one segment is genuinely weaker). The theory correctly predicts both regimes; neither paper
was fit to the other.

## 4. STEP 3 — Interstitial-pressure-vs-lymph-flow curve (Guyton/Aukland)

**Geometric derivation** (not a curve-fit): Bhave & Neilson (2011)'s own Figure 6, quoted verbatim
— *"Pi normally varies with interstitial volume (IFV). At low IFV, compliance is low and pressure
rises significantly. Once IFV increases 20-50% above euvolemia, compliance increases dramatically
and Pi essentially remains near constant allowing for edema formation."* Composing this saturating
Pi(IFV) with a monotonic pump-recruitment law Q_L(Pi) (calibrated so Q_L at the sibling thread's
own central-point Pi = -2.0 mmHg reproduces the STEP 1 baseline, 3.0 L/day) mechanically predicts a
flow-plateaus-above-Pi~0 shape — machine-verified on a numeric grid (-4 to +4 mmHg, 0.05 mmHg
steps), never eyeballed off a textbook figure:

| Ceiling assumption (q_max) | Steep-region slope (Pi -4→-1) | Plateau-region slope (Pi +1→+4) | Steepness ratio | Magnitude gate (≥5x) |
|---|---:|---:|---:|---:|
| 3x baseline (9.0 L/day) | 1.78 | 0.13 | **13.4x** | PASS |
| 6x baseline (18.0 L/day) | 2.39 | 0.62 | **3.8x** | fail |
| 10x baseline (30.0 L/day) | 2.75 | 1.72 | **1.6x** | fail |

**Two separate pre-registered bars, not merged**: the **directional** claim (post-0 region less
steep than pre-0 region AT ALL, ratio > 1) holds at **all 3** sweep points and is what gates
`overall_pass`. The stronger **magnitude** claim (ratio ≥ 5x) clears only the most conservative
ceiling (3x); it does NOT clear at 6x or 10x. This is reported as an honest **sensitivity
boundary**, not p-hacked by lowering the bar after seeing the result: with the recruitment curve's
steepness held fixed, self-consistency at Pi_normal mechanically drags the curve's inflection
toward 0 mmHg as the assumed ceiling grows (pi_mid: -1.43 → -0.70 → -0.23 mmHg across 3x/6x/10x),
weakening the plateau's sharpness within the sampled range. **The qualitative
plateau-above-Pi~0 claim is robust; its degree is sensitive to an unpinned ceiling parameter**,
most dramatic at the conservative end of the range this doc's own real (if
mechanistically-different) corroborations support — Havas 1997 (2-5x), Kramer/Renkin 1981 sheep
lung (4x).

**Steady-state existence/divergence — the task's own "assumes steady state, no edema" caveat,
machine-verified via explicit Euler integration of `dV/dt = Jv_fixed − Q_L(Pi(V))`** (60-day
simulated window, central sweep point):

| Regime | Final IFV/IFV₀ after 60 days |
|---|---:|
| Jv fixed at baseline (3.0 L/day, ≤ pump ceiling) | **0.80** (converges to a bounded fixed point) |
| Jv fixed at 1.5× pump ceiling (filtration exceeds max pump capacity) | **27.24** (grows without bound) |

A fixed point exists when the filtration line intersects the pump-recruitment curve; push
filtration above the curve's supremum and no intersection exists — decompensated edema, exactly
the boundary condition the task's own steady-state caveat names, now a computed result rather than
a prose assertion.

## 5. STEP 4 — Protein return: flux vs turnover, forced apart

Bhave & Neilson (2011), quoting Renkin (1986): **10 g albumin/hour** (= 240 g/day) transits from
plasma to lymph; ~50-60% of total body albumin sits extravascular at any instant. Using this
repo's own certified plasma volume (`fluid_compartments.py`, 3.00-3.07 L) and StatPearls' albumin
concentration band (3.5-5.0 g/dL):

| Quantity | Value |
|---|---:|
| Plasma albumin mass | 105.0 – 153.3 g |
| Derived TER (10 g/hr ÷ mass) | **6.52% – 9.52% per hour** |
| Implied extravascular-pool transit time | 0.44 – 0.96 days |
| Albumin synthesis/degradation turnover (StatPearls) | **10 – 15 g/day** |
| Flux ÷ synthesis ratio | **16× – 24×** |

**Forced adversary**: the one-way transcapillary-escape/lymph-return flux (240 g/day) is **16-24
times larger** than the synthesis/degradation turnover (10-15 g/day) — these are genuinely
different quantities (circulating-and-returning mass vs net metabolic replacement), and conflating
them would be a real, catchable error. The derived TER (6.5-9.5%/hr) sits in the same order of
magnitude as the commonly-cited clinical TER-albumin figure (~5%/hr) without being forced to match
it exactly — an honest, order-of-magnitude cross-check, not a precision claim.

**Liver cross-check** (Mobley et al. 1989): thoracic-duct-lymph/plasma-protein ratio falls from
0.58 to 0.48 when hepatic inflow is surgically clamped — direct confirmation that liver lymph is
protein-enriched relative to the rest of the body's lymph (consistent with hepatic sinusoidal
endothelium's known high permeability), corroborating why a whole-body-average lymph protein
concentration masks large organ-specific heterogeneity.

## 6. Falsifier verdict — stated exactly as pre-registered

> Does the modeled lymph flow MATCH the net capillary filtration rate the Starling thread
> independently computed (mass-balance closure, ~2-4 L/day) AND does the interstitial-
> pressure-vs-lymph-flow curve reproduce the measured Guyton/Aukland relationship (flow plateaus
> as Pi rises above ~0)?

**Mass-balance closure: PARTIAL, honestly qualified.** Against Renkin's own number alone, this
would be a tautology (§0) — flagged, not laundered. Against the **genuinely independent** Mobley
et al. (1989) dog-cannulation leg, allometrically scaled, the estimate (1.57-2.83 L/day) sits at or
just below the task's own 2-4 L/day band, in the physiologically-correct direction (motionless <
ambulatory). **This is real corroboration, not proof** — different species, an approximate
allometric scaling, and a resting (not typical-daily-activity) condition are all genuine, disclosed
limits.

**Pi-vs-lymph-flow curve: PASS on shape (directional), OPEN on magnitude.** The
monotonic-then-plateau shape is machine-verified across every point of the ceiling-multiple sweep
(all steepness ratios > 1). The dramatic 5x+ version of "plateau" holds only for the most
conservative (empirically-best-supported) ceiling assumption — reported as a sensitivity boundary,
not forced to a uniform pass.

## 7. Symmetric QC — nothing proven, spread reported (per the task's own pre-registration)

- **Whole-body lymph flow is genuinely hard to measure and estimates vary widely.** This doc
  deliberately reports TWO numbers rather than collapsing to one: the Renkin/Bhave-anchored
  [2,4] L/day band (shared common-mode with `capillary_starling.py`, disclosed) and the
  Mobley-1989-independent allometrically-scaled [1.57, 2.83] L/day band (different species/
  method/condition). The spread itself — not a false single point estimate — is the honest
  finding.
- **The filtration=lymph-return closure assumes steady state (no edema).** This is not just
  asserted: STEP 3's explicit numerical integration shows the model's own fixed point exists below
  the pump ceiling and provably diverges above it — the exact boundary the task's caveat names.
- **The exact fold-multiplier of the classical Guyton/Aukland curve was not independently
  re-derived from primary full text this session** (Physiol Rev 1971/1993, both paywalled, no PMC
  self-deposit found via elink — only "cited-by" lists). The curve's *shape* is derived from a
  directly-quoted mechanism (Bhave & Neilson's Fig.6) and confirmed by machine computation; its
  *ceiling magnitude* is an illustrative, disclosed sweep, bounded by (not equal to) two real but
  mechanistically-different reserve-capacity numbers (Havas 1997 exercise 2-5x; Kramer/Renkin 1981
  sheep-lung hypoproteinemia 4x). Held explicitly OPEN.
- **Two real bugs were caught mid-build** (§0) — a sign/placement error in the recruitment curve,
  and a degenerate self-consistency solve that silently collapsed a 3-point sweep to one repeated
  value while still "passing" gates. Both are disclosed in the script's own comments, not scrubbed
  from the history.

## 8. Gates — pre-registered, machine-checked, 16/16 PASS (with one explicit non-gating sensitivity finding)

```
shared_common_mode_with_capillary_starling_disclosed:        true  (identity-checked, not narrated)
mobley_independent_crosscheck_same_order_of_magnitude:       true
mobley_consistent_with_rest_lt_active_direction:             true
void_floor_edema_onset_plausible_1_to_48h:                   true  (6.0h)
havas_multiplier_direction_all_gt_1:                         true
havas_recompute_overlaps_reported_3to6x_claim:               true
davis2012_afterload_contractility_response_positive:         true
bottleneck_min_below_mean_confirmed:                         true
qL_vs_Pi_monotonic_nondecreasing_allsweep:                   true
qL_vs_Pi_plateau_direction_pass_allsweep_ratio_gt1:          true  (13.4x, 3.8x, 1.6x -- all >1)
sweep_is_genuinely_varying:                                  true  (post-bugfix machine check)
steady_state_exists_below_ceiling:                           true
steady_state_diverges_above_ceiling:                         true
ter_order_of_magnitude_plausible:                            true
flux_and_synthesis_are_distinct_quantities_confirmed:        true  (16-24x apart)
liver_protein_richness_direction_confirmed:                  true

NON-GATING (disclosed, not forced into overall_pass):
qL_vs_Pi_plateau_magnitude_5x_pass_by_point: {3x: true, 6x: false, 10x: false}

overall_pass: true   (16/16 gating checks; magnitude-sensitivity finding reported separately)
```

## 9. Couplings

- **capillary-Starling** (`scripts/msk/capillary_starling.py`, `data/capillary_starling/`): the
  filtration side this doc's mass-balance closure returns — shared common-mode (Renkin 1986)
  explicitly disclosed, not laundered as independent (§0, STEP 1).
- **fluid-compartments** (`scripts/msk/fluid_compartments.py`, `data/fluid_compartments/`): read
  read-only for interstitial volume (10.03 L) and plasma volume (3.00-3.07 L), used in the
  void-floor and protein-return calculations. Not modified.
- **immune / lymph nodes** (`ORG-LYMPHATIC-IMMUNE-TRAFFICKING`, OPEN): Renkin's own 8/4/4 L/day
  breakdown states ~4 L/day of the 8 L/day total is reabsorbed **in lymph nodes** before reaching
  the thoracic duct — lymph nodes are a real fluid-reabsorbing waypoint, not just an immune
  filter. Not otherwise modeled here.
- **`data/MECHANISM_ANCHOR_GRAPH.json`**: partially addresses `AUTO-CROSS-CHECK-THE-REVISED-
  STARLING-LYMPHAT` (the lymph-return side of the mass-balance question — not that node's own
  multi-tissue-bed tracer-method cross-check) and `AUTO-LYMPHATIC-IMMUNE-TRAFFICKING-ORGAN-CELL`'s
  low-churn "fluid-return Starling range" sub-block (not its node-surveillance-trafficking,
  meningeal-glymphatic, or lymphedema/oncology sub-blocks). **Not folded into the graph this
  session** — isolation discipline, another instance writes concurrently; left to the normal
  mechanism coordinator pipeline.

## 10. Confidence tier, stated precisely

**In-vivo/animal-anchored, with direct human measurements where available** — precisely as the
task specified, not blanket-claimed: Olszewski & Engeset (1980, human leg, in-vivo, intrinsic
contractility + flow — cannulation-adjacent) and Havas et al. (1997, human, in-vivo scintigraphy)
and Parving et al. (1979, human, ¹³¹I-albumin in-vivo turnover) are direct human in-vivo; Majgaard
et al. (2022) is direct human ex vivo (isolated vessel, surgical patients); Mobley et al. (1989,
dog) and Davis et al. (2012, rat) and Kramer et al. (1981, sheep) are direct animal in-vivo/ex-vivo
cannulation/perfusion. The whole-body volumetric anchor (Renkin 1986/Bhave & Neilson 2011) is a
narrative-review synthesis, not a single primary cohort — its shared-common-mode risk with the
sibling capillary-Starling thread is the reason STEP 1 does not rest on it alone.

## 11. Honest gaps (disclosed, not hidden)

- Whole-body lymph flow is NOT resolved to one number — reported as a spread across two
  decorrelated methods (§7), per the task's own instruction.
- Aukland & Reed (1993) and Guyton, Granger & Taylor (1971) — the two most authoritative primary
  sources for the classic Pi-vs-lymph-flow curve — have **no accessible full text** this session
  (both paywalled, no PMC self-deposit; elink returns only citing-article lists). The curve's shape
  is derived from a directly-quoted mechanism in a different, open-access paper (Bhave & Neilson
  2011's Fig.6); its exact magnitude is swept/illustrative, not primary-pinned.
- Olszewski & Engeset (1980)'s own measured numbers (the one truly direct human whole-leg
  intrinsic-contractility-and-flow dataset) were not accessible this session — cited for existence
  only.
- The pump-chain model is validated STRUCTURALLY (bottleneck/conservation argument, machine-
  checked and literature-bracketed) — it does NOT derive an absolute L/day figure, since human
  lymphangion stroke volume was not independently pinned this session (avoiding fabrication of an
  unmeasured constant).
- The Mobley et al. (1989) dog mass is not stated in the accessible abstract — swept over a
  plausible range [15,30] kg rather than assumed at one unverified value.
- Normal (non-myxedema) human TER-albumin %/hour baseline was not extracted from Parving et al.
  (1979)'s abstract (comparison-framed, not stated) — this doc's own 6.5-9.5%/hr figure is a
  computed derivation from a different real number (Bhave/Renkin's 10 g/hr), not independently
  cross-validated against a primary clinical TER cohort this session.
- The void-floor and steady-state dynamics checks are disclosed simplifications (static/first-
  order filtration; a single lumped compartment) — not a full time-domain multi-compartment edema
  model.

## 12. Files

- `scripts/msk/lymphatic_return.py` — the model (4 steps + forced adversaries + gates). Run with
  `.venv-msk/bin/python3 scripts/msk/lymphatic_return.py` (<1s wall time, no OpenSim dependency,
  pure Python stdlib — `json`, `math`, `os`, `statistics`; reads two sibling JSONs read-only:
  `data/capillary_starling/capillary_starling_results.json`,
  `data/fluid_compartments/fluid_compartments_results.json`).
- `data/lymphatic_return/lymphatic_return_results.json` — full evidence (17 citations, per-step
  results, sweeps, gates, non-gating sensitivity findings).
- This doc: `docs/MECHANISM_LYMPHATIC_RETURN.md`.

## 13. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/lymphatic_return.py
```
No inputs beyond the two read-only sibling JSONs (both already present in this repo). Deterministic,
pure Python/stdlib. New files only; nothing existing modified; no git operations; does not touch
`data/MECHANISM_ANCHOR_GRAPH.json`.
