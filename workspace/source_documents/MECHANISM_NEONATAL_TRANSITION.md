# MECHANISM NEONATAL TRANSITION — first-breath opening pressure (Laplace geometry, inverted) + fetal-lung-liquid clearance forced against a "passive-drainage-only" adversary across 3 species/methods (2026-07-22)

Builds and MEASURES the two-part falsifier for birth's cardiopulmonary transition: **(a)** the
first-breath opening pressure (~40-80 cmH2O, far above tidal) is reproduced not by a new curve fit
but by **inverting** the sibling surfactant cert's own verified Laplace relation `P=2*gamma/r` —
asking what effective radius of curvature the MEASURED pressure implies, given the SAME
bare/unspread surface-tension bracket already anchored in `docs/MECHANISM_PULMONARY_SURFACTANT.md`
— and checking that implied radius is geometrically sane (much smaller than the mature alveolus);
**(b)** fetal-lung-liquid clearance is forced against a "passive drainage only" adversary using
DIRECT measured necessity evidence — a genetic knockout (mouse), a pharmacological blockade (fetal
sheep), and a clinical cohort (human) — three decorrelated species/methods, not a single simulated
curve, all independently showing that removing active Na+ transport (not mechanical squeeze, not
Starling forces) makes clearance fail. Script: `scripts/msk/neonatal_transition.py`. Raw evidence:
`data/neonatal_transition/neonatal_transition_results.json`. Curated citation/gate summary:
`docs/MECHANISM_NEONATAL_TRANSITION_evidence.json`.

**This is a HYPOTHESIS for independent QC** (return-only wave discipline, `docs/MECHANISM_HARDENED_CONVENTIONS.md` §4b) — a designed-and-computed cell, not yet folded into the graph.

## 0. Falsifiers (pre-registered, matches the task spec) + verdict up front

1. **Claim A (opening pressure)** — the first-breath transpulmonary opening pressure (~40-80
   cmH2O, task-specified / Karlberg-Koch-anchored) is far above tidal-breathing pressure, and this
   gap is explained by surface tension acting at a much smaller effective radius of curvature
   before the lung is aerated, NOT by chest-wall stiffness or by the mature-alveolus geometry alone
   → **PASS, 5/5 machine gates** (§4 — inverting `P=2*gamma/r` at the measured pressure range,
   using the bare/unspread tension bracket already anchored in the surfactant cert, implies an
   effective radius of **12.7-35.7 micron** — 2.8x-7.9x smaller than the mature 100-micron alveolus
   used throughout that cert; the mature-geometry-alone ceiling, even at ZERO surfactant function,
   tops out at **14.28 cmH2O**, well short of the measured 40 cmH2O floor; infant chest wall is
   **2.86x MORE compliant than the lung** — Papastamelos et al 1995 — ruling out chest-wall
   stiffness as the driver).
2. **Claim B (fluid clearance, PRIMARY)** — a "passive drainage only" adversary (mechanical
   squeeze + Starling forces, present in a normally-born animal, MINUS functional Na+/ENaC
   transport) must FAIL to clear lung liquid at the measured/observed rate → **PASS, 7/7 machine
   gates** (§5 — **genetic**: alpha-ENaC(-/-) mice, Na+ transport abolished, die within **40 h**
   of birth from failure to clear lung liquid, Hummler et al 1996, PMID 8589728; **pharmacological**:
   amiloride at 10⁻⁴ M abolishes the adrenaline-induced reabsorption response in fetal lambs, **25x**
   above its own measured KI, Olver/Ramsden et al 1986, PMID 3795077; **clinical**: bypassing labor
   entirely (elective caesarean) carries an OR **3.9** (95% CI 2.4-6.5) for composite respiratory
   morbidity at 37 weeks, falling monotonically to 1.9 at 39 weeks as maturity/proximity-to-labor
   increases, Hansen et al 2008, PMID 18077440 — **3 decorrelated species/methods** (mouse-genetic,
   sheep-pharmacological, human-clinical-epidemiological), not one).
3. **Mechanism/timing (labor catecholamine switch)** — the fetal lung's Cl⁻-driven SECRETION
   switches to Na+-driven ABSORPTION specifically coincident with the measured labor catecholamine
   surge, not before → **PASS** (§5 — plasma adrenaline rises **78.85x** from early labor (0.087
   ng/mL) to the last 50 min (6.86 ng/mL), Brown/Olver et al 1983, PMID 6655575; the SAME paper's
   adrenaline-infusion dose-response in NON-labouring fetuses reproduces the identical
   secretion→absorption switch, establishing sufficiency, not mere correlation; Olver & Strang 1974,
   PMID 4443921, independently establishes the PRE-labor baseline: Cl⁻ actively secreted, Na⁺ moves
   passively — the switch has a genuine "before" state to switch FROM).
4. **Secondary/illustrative cross-check (viscous resistance)** — NOT independently re-derived this
   session; reused by reference from the surfactant cert's own Hooper et al 2013 anchor (liquid
   resistance "≈100 times greater than air", PMID 24035400). A bulk Newtonian water/air viscosity
   ratio at 37°C (**36.4x**, standard tabulated constants, not a biomedical claim) is directionally
   consistent (>>1) but explicitly NOT claimed to reproduce the measured ~100x figure (§6, §8).

**Overall: PASS, 13/13 machine gates** (§7 prints the full pre-registered gate block; gates are
graded by decisiveness in §7, not presented as uniformly load-bearing — symmetric QC applies to
passes, not just kills). Deterministic (2 independent runs byte-identical,
`md5sum 84b80672de1831a9bd4aec0f9b971dca`), pure Python arithmetic (no simulation, no fitting),
runtime <1 second.

## 1. Geometric structure (derive from the geometry, not heuristics)

**Claim A is an INVERSION of an already-verified relation, not a new curve fit.** The sibling
surfactant cert established `f(r) = 2*gamma/r` (Laplace's law) and evaluated it FORWARD at the
mature alveolar radius (`r=1e-4 m`) to get pressures of 0.41-14.28 cmH2O depending on surfactant
status. This cert runs the SAME relation BACKWARD: given a MEASURED pressure (the first-breath
opening pressure, 40-80 cmH2O) and the SAME bare/unspread tension bracket (50-70 mN/m — physically
appropriate here because a freshly-forming interface has not yet had surfactant spread across it,
the same "not yet functional" state the surfactant cert uses for its own bare-saline comparator),
what radius does the geometry REQUIRE?

```
r_implied = 2*gamma / P
```

This is a genuine falsifier, not a tautology: **r_implied could have come out equal to or larger
than the mature 100-micron alveolar radius**, which would have REFUTED the "tiny, not-yet-opened
interface" mechanism (it would mean the pressure gap has nothing to do with curvature at all).
Instead, across the full 3(gamma)x2(pressure) grid, r_implied lands at **12.7-35.7 micron** —
consistently 2.8x-7.9x SMALLER than the mature radius, i.e. geometrically consistent with "the
interface starts at much higher curvature before the first breath than it settles to once FRC is
established." The direction of this consistency check is exactly the kind of over-determination the
discipline requires: two independent numbers (a literature-measured pressure, a literature-anchored
tension bracket already used for an unrelated purpose in the sibling cert) combine through one
geometric relation to produce a THIRD number (implied radius) that is checked against a FOURTH,
independent anchor (the mature alveolar radius, itself measured/used elsewhere) — no free parameter
was tuned to make this land in the physiologically sane range.

**Why the "mature-geometry-alone" adversary is forced to its strongest form and still fails.** The
most tempting shortcut here is to skip re-deriving anything and just assert "surface tension
matters, so pressure is high." The fair adversary is: *"maybe no radius change is needed at all —
maybe bare (zero-surfactant-function) tension AT THE ALREADY-MATURE alveolar geometry is enough to
explain the measured first-breath pressure."* Computed directly (§4, `mature_bare_max_cmH2O`), the
worst case (70 mN/m, r=100 micron) tops out at **14.28 cmH2O** — this is the IDENTICAL number
already reported in the surfactant cert's own F1 table, reused here as a cross-check, not
recomputed with new inputs — a full **2.8x short of the measured 40 cmH2O floor**. The adversary
is not merely disfavored; it is quantitatively insufficient by nearly 3-fold at the LOW end of the
measured range and by nearly 6-fold at the high end. A radius change (curvature increase) is
therefore NECESSARY, not merely sufficient, to close this gap under the bare-tension ceiling.

**A second, independent adversary (chest-wall stiffness) is forced and also falls.** A different
confound: maybe the high pressure isn't about the LUNG at all, but about a stiff neonatal chest
wall resisting expansion. Papastamelos et al (1995, PMID 7713809, live-verified) directly measured
chest-wall-vs-lung compliance in infants under 1 year: **Cw/Cl = 2.86 ± 1.06** — the chest wall is
nearly 3x MORE compliant (easier to distort) than the lung, the OPPOSITE of what the "stiff chest
wall" adversary requires. This rules out chest-wall stiffness as the source of the high opening
pressure and redirects the explanatory burden back onto the lung itself (surface tension + liquid),
consistent with Claim A's own mechanism.

**Claim B's geometry is a threshold/switch, not a smooth ramp — and the switch is measured, not
assumed.** Brown, Olver, Ramsden, Strang & Walters (1983, PMID 6655575) show the fetal lamb's
response to adrenaline is GESTATION-DEPENDENT and changes QUALITATIVE SIGN: between 120-130 days,
adrenaline only slows secretion; after 130 days, the SAME stimulus produces net absorption — a
genuine bifurcation in the response, structurally the same kind of "control-parameter crosses a
threshold and the qualitative behavior flips" geometry as the surfactant cert's own n=1/2 stability
threshold and this repo's `MECHANISM_PARTURITION_MYOMETRIUM.md` oxytocin-gain threshold — though here
the threshold is crossed by GESTATIONAL MATURATION (the epithelium's own Ai sensitivity falls
**14.83x**, from 0.43 to 0.029 ng/mL, between 132-4 days and >140 days) rather than swept as a
free model parameter.

## 2. Method, in one paragraph

`neonatal_transition.py` (pure Python arithmetic, no simulation) computes three parts. **Part 1**
evaluates the SAME `P=2*gamma/r` relation as the surfactant cert (reused, not re-derived) forward
at the mature radius (reproducing that cert's own 0.41-14.28 cmH2O numbers as an internal
consistency check) and INVERTED at the measured first-breath pressure range to get the implied
opening radius (§1), plus the two forced-adversary checks (mature-geometry-alone ceiling;
chest-wall compliance ratio) and one order-of-magnitude cross-check (a clinical sustained-lung-
inflation protocol pressure, Lista & Castoldi 2010, PMID 21089711). **Part 2** computes fold-changes
and margins-over-threshold directly from literature-extracted numbers (catecholamine surge across
labor; amiloride block margin over its own measured KI; gestational sensitization of the Ai
threshold) and tabulates the qualitative necessity results (genetic knockout, pharmacological
blockade, clinical cohort) as an explicit 3-species/method diversity count. **Part 3** is a
secondary, explicitly-labeled-illustrative bulk-viscosity cross-check that does NOT attempt to
reproduce Hooper et al 2013's measured ~100x liquid-vs-air airway-resistance figure (reused by
reference from the sibling surfactant cert instead). All inputs are either (i) reused, verbatim,
from a citation already live-verified in a sibling doc this repo, or (ii) freshly live-verified via
NCBI eutils this session (§3) — no recalled-and-trusted numbers.

## 3. Citations — verified LIVE this session (NCBI eutils esearch/esummary/efetch)

10 papers freshly fetched live this session (full abstract text extracted for 8 of 10; 2 pre-1975
Karlberg-series records have no MEDLINE-indexed abstract, disclosed not hidden, §8). 2 further
citations are REUSED BY REFERENCE from `docs/MECHANISM_PULMONARY_SURFACTANT.md` (already
live-verified in that sibling session, not re-fetched here).

| # | Citation | PMID / DOI | Tier | Role / number extracted live |
|---|---|---|---|---|
| 1 | Hummler E, Barker P, Gatzy J, Beermann F, Verdumo C, Schmidt A, Boucher R, Rossier BC (1996). Early death due to defective neonatal lung liquid clearance in alpha-ENaC-deficient mice. *Nat Genet* 12(3):325-8. | **8589728**, DOI 10.1038/ng0396-325 | full abstract fetched live | THE decisive genetic-necessity anchor for Claim B: alpha-ENaC(-/-) abolishes amiloride-sensitive Na+ transport in airway epithelia; neonates "died within 40 h of birth from failure to clear their lungs of liquid." |
| 2 | Olver RE, Strang LB (1974). Ion fluxes across the pulmonary epithelium and the secretion of lung liquid in the foetal lamb. *J Physiol* 241(2):327-57. | **4443921**, DOI 10.1113/jphysiol.1974.sp010659 | full abstract fetched live | THE baseline fetal-SECRETION-phase anchor: Cl⁻ actively transported plasma→lung liquid; Na⁺ moves PASSIVELY down its electrochemical gradient during this pre-labor phase — establishes the "before" state Claim B's switch acts on. |
| 3 | Brown MJ, Olver RE, Ramsden CA, Strang LB, Walters DV (1983). Effects of adrenaline and of spontaneous labour on the secretion and absorption of lung liquid in the fetal lamb. *J Physiol* 344:137-52. | **6655575**, DOI 10.1113/jphysiol.1983.sp014929 | full abstract fetched live | THE labor-catecholamine-switch anchor: gestation-dependent sign-flip of the adrenaline response (slows secretion at 120-130d; causes absorption after 130d); measured plasma adrenaline 0.087->6.86->7.17 ng/mL (early labor -> last 50 min -> early postnatal) and noradrenaline 1.71->12.14->9.10 ng/mL; upper-airway one-way-valve mechanism (outflow only). |
| 4 | Olver RE, Ramsden CA, Strang LB, Walters DV (1986). The role of amiloride-blockable sodium transport in adrenaline-induced lung liquid reabsorption in the fetal lamb. *J Physiol* 376:321-40. | **3795077**, DOI 10.1113/jphysiol.1986.sp016156 | full abstract fetched live | THE pharmacological-necessity anchor (Claim B, decorrelated species+method from Hummler): amiloride at 10⁻⁴ M "abolished the changes in p.d. and ion flux induced by adrenaline"; KI (50% inhibition) = 4×10⁻⁶ M. |
| 5 | Hansen AK, Wisborg K, Uldbjerg N, Henriksen TB (2008). Risk of respiratory morbidity in term infants delivered by elective caesarean section: cohort study. *BMJ* 336(7635):85-7. | **18077440**, DOI 10.1136/bmj.39405.539282.BE | full abstract fetched live | THE clinical cross-check anchor (34,458-baby Aarhus cohort): elective CS (labor bypassed) OR for composite respiratory morbidity (TTN + RDS + PPHN named explicitly as the outcome bundle) = 3.9 (37wk), 3.0 (38wk), 1.9 (39wk); serious morbidity OR 5.0 at 37wk. |
| 6 | Hooper SB, Polglase GR, Roehr CC (2015). Cardiopulmonary changes with aeration of the newborn lung. *Paediatr Respir Rev* 16(3):147-50. | **25870083**, DOI 10.1016/j.prrv.2015.03.003 | full abstract fetched live | THE coupling anchor to the fetal-circulation cert: "lung aeration triggers the increase in pulmonary blood flow (PBF) at birth," securing left-ventricular preload — the causal link this doc's aeration/clearance mechanics feeds into the circulatory shunt-closure side. |
| 7 | Steinhorn RH (2016). Advances in Neonatal Pulmonary Hypertension. *Neonatology* 109(4):334-44. | **27251312**, DOI 10.1159/000444895 | full abstract fetched live | Modern clinical anchor for the PPHN dysfunction coupling: "a surprisingly common event in the neonatal intensive care unit," confirms active clinical relevance of failed PVR fall. |
| 8 | Papastamelos C, Panitch HB, England SE, Allen JL (1995). Developmental changes in chest wall compliance in infancy and early childhood. *J Appl Physiol* 78(1):179-84. | **7713809**, DOI 10.1152/jappl.1995.78.1.179 | full abstract fetched live | Adversary-elimination anchor for Claim A: Cw/Cl = 2.86±1.06 in infants <1yr — chest wall nearly 3x MORE compliant than lung, ruling out chest-wall stiffness as the opening-pressure driver. |
| 9 | Lista G, Castoldi F (2010). Alveolar recruitment in the delivery room: sustained lung inflation. *Minerva Pediatr* 62(3 Suppl 1):17-8. | **21089711**, no DOI in record | full abstract fetched live | Order-of-magnitude clinical cross-check: sustained lung inflation protocol uses "a peak pressure of 25-30 cm H2O for 10-20 seconds" — same order of magnitude as the measured spontaneous first-breath peak, disclosed as consistent-not-identical (different protocol: sustained vs brief natural peak). |
| 10 | Karlberg P, Koch G (1962). Respiratory studies in newborn infants. III. Development of mechanics of breathing during the first week of life. *Acta Paediatr Suppl* 135:121-9. + Karlberg P et al (1960). Respiratory studies in newborn infants. I. Apparatus and methods... *Acta Paediatr* 49. | **14453970** (Part III) / **14453970**; **14404498** (Part I) | bibliographic only (no MEDLINE abstract, pre-1975/1960s records) | THE task-named Karlberg newborn-respiratory-mechanics series — bibliographically confirmed live (title/journal/volume/pages/DOI), establishing the series exists and covers exactly this subject (breathing mechanics in the first minutes/days of life). The specific ~40-80 cmH2O opening-pressure figure is the task's own specified, textbook-transmitted magnitude — NOT independently re-extracted from a live full-text fetch of this series this session (honest gap, §8; the SPECIFIC installment reporting the pressure-volume loop numbers, commonly cited as "Part II" of a related Karlberg/Cherry/Escardo series, was searched for directly and not located this session — see §8). |

**Reused by reference (already live-verified in `docs/MECHANISM_PULMONARY_SURFACTANT.md`, NOT
re-fetched this session):**

- **Avery ME, Mead J (1959)**, PMID 13649082 — RDS/surfactant-deficiency dysfunction anchor.
- **Hooper SB, Siew ML, Kitchen MJ, te Pas AB (2013)**, PMID 24035400 — FRC establishment; liquid
  viscosity "≈100 times greater than air"; the viscous-resistance anchor for Part 3 (§6).

**Honest disclosure (searched, not found, not fabricated):** a paper titled with "onset of
respiration" or "first minutes of life" (the specific installment of the Karlberg/Cherry/Escardo
newborn-respiratory series most commonly cited for the actual first-breath pressure-volume loop
data) was searched directly (`Karlberg P[Author] AND "onset of respiration"`, `AND "first breath"`,
`AND Cherry RB[Author]`, `AND Escardo FE[Author]`, plus a broad `Karlberg P[Author]` listing of all
118 of this prolific author's PubMed-indexed papers, scanned by esummary title) — the exact
installment was **not** definitively identified this session. Part I (methods, 1960, PMID 14404498)
and Part III (first-week development, 1962, PMID 14453970) of the closely-related "Respiratory
studies in newborn infants" series WERE found and bibliographically confirmed, establishing the
series and its subject matter beyond reasonable doubt, but the specific numeric P-V loop data are
not independently re-extracted this session — an honest gap, not a fabricated identifier, exactly
the same disclosure pattern as the surfactant cert's von Neergaard/Bachofen gaps.

## 4. Part 1 results — opening-pressure geometry (machine-computed)

| quantity | measured | pre-registered gate | anchor |
|---|---:|---|---|
| mature-radius (r=100 micron) bare-tension pressure, 50/60/70 mN/m | 10.20 / 12.24 / **14.28** cmH2O | — | identical formula+inputs as surfactant cert F1 (internal cross-check) |
| mature-radius surfactant-active pressure, 2/5/10 mN/m | 0.41 / 1.02 / 2.04 cmH2O | — | ditto |
| measured first-breath opening pressure (task-specified/Karlberg) | 40 — 80 cmH2O | — | Karlberg-Koch series (bibliographic, §3#10) |
| **r_implied (inverse Laplace, bare tension bracket)** | **12.7 — 35.7 micron** | < 50 micron (half of r_mature) | **PASS** |
| mature-bare-alone ceiling vs measured floor | 14.28 vs 40 cmH2O | ceiling < floor | **PASS** (2.8x short) |
| first-breath / mature-tidal pressure ratio | **19.6x — 196x** | ≥ 10x | **PASS** |
| chest-wall / lung compliance ratio, infants <1yr | **2.86 ± 1.06** | > 1 (wall more compliant) | **PASS** (Papastamelos 1995) |
| SLI clinical pressure (order-of-magnitude cross-check) | 25 — 30 cmH2O | same order of magnitude as 40-80 | **PASS (weak/supportive only, §7)** |

**All 5 gates PASS.** The inverse-Laplace-implied radius (12.7-35.7 micron) sits well inside the
physiologically sane range — smaller than the mature alveolus by 2.8x-7.9x, consistent with "the
air-liquid interface starts at much higher curvature before aeration than it settles to at FRC" —
and BOTH forced adversaries (mature-geometry-alone; chest-wall stiffness) are quantitatively
insufficient, not merely disfavored.

## 5. Part 2 results — fluid-clearance adversary-forcing across 3 species/methods (PRIMARY claim)

| quantity | measured | pre-registered gate | anchor |
|---|---:|---|---|
| adrenaline fold-change, early labor -> last 50 min | 0.087 -> 6.86 ng/mL = **78.85x** | ≥ 10x | **PASS** (Brown/Olver 1983) |
| noradrenaline fold-change, same window | 1.71 -> 12.14 ng/mL = **7.10x** | — (reported, not gated) | Brown/Olver 1983 |
| Ai-threshold gestational sensitization | 0.43 -> 0.029 ng/mL = **14.83x more sensitive** | > 1x (sensitizes) | **PASS** (Brown/Olver 1983) |
| amiloride block margin over its own KI | 1×10⁻⁴ M applied / 4×10⁻⁶ M KI = **25.0x** | ≥ 5x | **PASS** (Olver/Ramsden 1986) |
| ENaC-KO phenotype | Na+ transport abolished; death within **40 h**, failure to clear liquid | acute (≤72h), categorical necessity | **PASS** (Hummler 1996) |
| elective-CS-without-labor OR, 37wk (95% CI) | **3.9** (2.4 — 6.5) | lower CI > 1 | **PASS** (Hansen 2008) |
| OR trend with gestational maturity | 3.9 (37wk) > 3.0 (38wk) > 1.9 (39wk) | monotonic decrease | **PASS** (Hansen 2008) |
| decorrelated species/methods | mouse-genetic, sheep-pharmacological, sheep-endocrine, human-clinical | ≥ 3 | **PASS** (4 counted) |

**All 7 gates PASS.** The decisive result is qualitative-categorical, not merely a threshold
crossing: in EVERY one of 3 independent species/methods (mouse genetic knockout; fetal-lamb acute
pharmacological blockade; human clinical cohort bypassing labor), removing or bypassing active
Na+-transport machinery — while leaving whatever passive/mechanical mechanisms a normal birth
provides fully intact — causes measurably worse fluid-clearance outcomes (mouse: death from
unresorbed liquid; lamb: abolished reabsorption response; human: 2-6x elevated respiratory
morbidity). The "passive drainage only" adversary is not merely under-favored by a model choice —
it is empirically falsified in three independent experimental systems.

**Honest strength grading (symmetric QC — a pass needs the same scrutiny as a kill, §7):** the
Hummler 1996 abstract does not itself state the delivery mode of the ENaC-KO litters. Standard
gene-targeting mouse phenotyping uses normal (vaginal, dam-delivered) birth unless stated otherwise
— a reasonable, standard inference, but **not independently confirmed from the fetched abstract
text this session** (§8). The core necessity conclusion (active Na+ transport is required; genetic
loss is not compensated) holds regardless of delivery mode; the STRONGEST framing ("mechanical
squeeze present, still fails") is the standard-inference version, disclosed as such, not the
directly-quoted version.

## 6. Part 3 — secondary, explicitly-illustrative viscous-resistance cross-check

| quantity | measured/computed | note |
|---|---:|---|
| mu_water at 37C | 0.6913 mPa·s | standard tabulated physical constant, NOT NCBI-verified (not a biomedical claim) |
| mu_air at 37C | 0.0190 mPa·s | ditto |
| bulk viscosity ratio (water/air) | **36.4x** | simplified Newtonian lower-bound estimate |
| Hooper et al 2013 measured whole-airway resistance ratio | **~100x** | PMID 24035400, reused by reference (surfactant cert) |

**Gate PASS (directional consistency only, explicitly not decisive):** both numbers exceed 10x
(same DIRECTION: liquid resistance >> air resistance, by more than an order of magnitude), but this
cert does **not** claim the bulk-viscosity ratio reproduces Hooper's measured ~100x figure — real
fetal lung liquid is not pure water (it carries protein and other solutes) and real airway
resistance includes geometric and non-Newtonian effects a bulk viscosity ratio cannot capture. This
gap is disclosed, not smoothed over: presenting 36.4x as if it "confirmed" the measured ~100x would
be exactly the kind of self-flattering, uninspected consistency check the discipline warns against.
This part is scoped as illustrative context for Claim A's second (viscous) mechanism, not as
load-bearing evidence.

## 7. Pre-registered gates — machine-printed, not narrated, graded by decisiveness

```
A1_r_implied_lt_half_mature_radius:                     PASS  (35.69 micron < 50)   [DECISIVE]
A2_mature_bare_alone_cannot_reach_measured_min:         PASS  (14.28 < 40 cmH2O)    [DECISIVE -- forced adversary #1]
A3_first_breath_over_mature_tidal_ratio_ge_10x:         PASS  (19.6x - 196.1x)      [DECISIVE]
A4_chest_wall_more_compliant_than_lung:                 PASS  (2.86 > 1)            [DECISIVE -- forced adversary #2]
A5_SLI_same_order_of_magnitude_as_first_breath:         PASS  (25-30 vs 40-80)      [WEAK/SUPPORTIVE ONLY -- disclosed, not decisive]
B1_catecholamine_surge_ge_10x:                          PASS  (78.85x)             [DECISIVE]
B2_ai_threshold_sensitizes_with_gestation:               PASS  (14.83x)             [SUPPORTIVE -- mechanistic nuance, not itself the necessity claim]
B3_amiloride_margin_ge_5x:                              PASS  (25.0x)              [DECISIVE -- forced adversary, pharmacological]
B4_enac_ko_fails_despite_normal_birth_mechanics:        PASS  (40h, categorical)    [DECISIVE -- forced adversary, genetic; delivery-mode inferred not quoted, see honest_gaps]
B5_cs_or_37wk_lower_CI_excludes_1:                      PASS  (CI 2.4-6.5)          [DECISIVE -- clinical cross-check]
B6_cs_or_monotonic_with_gestational_maturity:           PASS  (3.9>3.0>1.9)         [SUPPORTIVE -- dose-response-style consistency]
B7_decorrelated_methods_ge_3:                           PASS  (4 counted)           [STRUCTURAL -- the diverse-instance-space requirement itself]
C1_bulk_viscosity_ratio_gt_10x_consistent_direction:    PASS  (36.4x)               [WEAK/ILLUSTRATIVE ONLY -- explicitly not claimed to match Hooper's ~100x]

VERDICT: 13/13 boolean gates PASS (8 decisive/structural, 3 supportive, 2 explicitly weak/illustrative-only — none hidden, none silently promoted). overall_pass = True
```

Deterministic: 2 independent runs verified byte-identical, `md5sum 84b80672de1831a9bd4aec0f9b971dca`.

## 8. Honest gaps (disclosed, not hidden)

- **The specific Karlberg-series installment reporting the actual first-breath P-V loop numbers**
  (commonly cited, in secondary/textbook sources, as "Part II" of a related Karlberg/Cherry/Escardo
  series) was searched for directly this session and **not located** — Part I (methods, PMID
  14404498) and Part III (first-week development, PMID 14453970) of the closely-related series WERE
  found and bibliographically confirmed. The ~40-80 cmH2O figure remains the task's own specified,
  textbook-transmitted magnitude, cross-checked only indirectly (the inverse-Laplace geometric
  consistency check, §1/§4; the Lista 2010 SLI order-of-magnitude comparison, §4) — not a direct
  live re-extraction of the original tabulated pressure-volume data.
- **The "1/3 mechanical-squeeze : 2/3 active-transport" split** sometimes quoted in older secondary
  /textbook sources for how much of total lung liquid clearance is attributable to each mechanism
  was deliberately **not** used or verified this session. This cert instead relies on the
  Hummler/Olver-Ramsden NECESSITY evidence (removing active transport causes categorical failure,
  regardless of exact partition), which is a stronger and more directly falsifiable form of evidence
  than a quantitative split whose original source and precision this session did not verify — a
  disclosed scope choice, not an oversight.
- **No standalone quantitative "passive-only predicted clearance half-life" was computed.** Unlike
  the surfactant cert's two-alveolus ODE (which swept a free parameter through a predicted
  bifurcation), this cert's Part 2 adversary-forcing rests entirely on REAL measured necessity
  phenotypes (a genetic knockout, a pharmacological block, a clinical cohort) rather than a
  simulated passive-only rate constant — arguably stronger evidence (real measured outcomes, not an
  assumed model), but it means no simulated "how long would passive-only clearance take" number
  exists in this cert.
- **Hummler et al 1996's delivery mode for the ENaC-KO litters is inferred, not directly quoted**
  from the fetched abstract (standard mouse-phenotyping practice = vaginal, dam-delivered, unless
  otherwise stated) — disclosed explicitly in §5 and in the raw JSON's
  `enac_ko_categorical_result` field. The core necessity conclusion is robust to this either way.
- **The inverse-Laplace calculation (§1, §4) assumes the pre-aeration interface's surface tension
  is bare/unspread-film-like** (same 50-70 mN/m bracket as the surfactant cert's own bare-saline
  comparator) rather than some intermediate, partially-spread value — a physiologically reasonable
  but not directly-measured-this-session assumption, explicitly disclosed (mirrors the surfactant
  cert's own disclosed assumptions about its gamma(A) calibration).
- **Part 3's bulk-viscosity-ratio cross-check (36.4x) is explicitly NOT claimed to reproduce**
  Hooper et al 2013's measured ~100x whole-airway liquid-vs-air resistance figure — the gap
  (attributable to real airway geometry, non-Newtonian effects, and fetal lung liquid's protein
  content vs. pure water) is disclosed, not smoothed over (§6).
- **Papastamelos et al 1995's chest-wall-compliance measurement** was made in sedated infants 2
  weeks to 3.5 years old, not specifically at the moment of birth/first breath — the closest
  live-verified developmental data available, not itself a birth-moment measurement (a disclosed
  extrapolation, though the direction of the ratio — chest wall more compliant than lung — is not
  expected to reverse in the first minutes of life).
- **This cert does not quantitatively partition** how much of the measured 40-80 cmH2O first-breath
  pressure is attributable to surface tension (Part 1) vs. viscous liquid-displacement resistance
  (Part 3) — both are real, independently documented, decorrelated contributors; no claim is made
  about their relative shares.
- **PPHN and meconium aspiration** (task-named dysfunctions) are covered here only at the level of
  a single modern clinical-review citation each (Steinhorn 2016 for PPHN; meconium aspiration is
  named in the task's dysfunction list but no dedicated citation was fetched this session for it,
  given this cert's scope is the aeration/clearance mechanism rather than meconium-specific
  pathophysiology) — flagged as a lighter-touch coupling, not a gap in the PRIMARY falsifier.

## 9. Couplings + confidence tier

**`couples_to`:**

- **Fetal-circulation cert (IN FLIGHT per task spec)** — checked live this session
  (`docs/`, `data/MECHANISM_ANCHOR_GRAPH.json`): **no doc or graph node yet exists** for
  fetal/ductus-arteriosus/foramen-ovale/PPHN circulatory-transition content as of this session (not
  a duplicate; a genuine open coupling point). This doc supplies the RESPIRATORY half of the
  transition — Hooper et al 2015 (§3#6) directly states lung aeration TRIGGERS the pulmonary-blood-
  flow increase that the circulatory cert's own PVR-fall/shunt-closure mechanism depends on; PPHN
  (§3#7, Steinhorn 2016) is the shared dysfunction-mode where this coupling fails.
- **`docs/MECHANISM_PULMONARY_SURFACTANT.md`** — direct numeric reuse, not just topical adjacency:
  this doc's entire Part 1 (§1, §4) is built by INVERTING that cert's own `P=2*gamma/r` relation and
  bare/surfactant tension brackets; the mature-radius bare-tension ceiling (14.28 cmH2O) is the
  IDENTICAL number from that cert's own F1 table, reused as a forced-adversary cross-check here, not
  independently recomputed with new inputs. Also reuses Hooper 2013 (FRC/viscosity) and Avery-Mead
  1959 (RDS) by reference.
- **`docs/MECHANISM_PARTURITION_MYOMETRIUM.md`** (checked read-only, not modified) — that cert models
  the OXYTOCIN/Ferguson-reflex mechanism driving myometrial contraction; this cert's labor-switch
  mechanism (§1, §5) is driven by a DIFFERENT hormone axis entirely (fetal adrenal-medulla
  catecholamines, Brown/Olver 1983) acting on a different target organ (fetal lung epithelium, not
  myometrium). The shared coupling point is the EVENT of labor itself (both mechanisms are
  labor-triggered, on different molecular pathways) — not a shared molecular pathway, disclosed
  precisely rather than overclaimed as a direct mechanistic link.
- **`docs/MECHANISM_MUCOCILIARY_CLEARANCE.md`** (checked read-only) — a DIFFERENT clearance
  modality (ciliary-driven mucus transport in the mature/pediatric airway) vs. this cert's
  transepithelial Na+-driven LIQUID absorption at birth — adjacent topic (airway clearance,
  broadly), decorrelated mechanism, no direct numeric dependency.
- **`docs/MECHANISM_RESPIRATORY.md`, `docs/MECHANISM_CARDIAC_OUTPUT.md`** (checked read-only) — both
  are adult/pediatric steady-state organ-system layers (VO2/ventilation; Frank-Starling stroke
  volume), operating one developmental stage above this doc's birth-transition mechanics; no direct
  numeric dependency, same organ systems.
- **Dysfunction coupling (function<->dysfunction, task-named):** RDS (surfactant deficiency, reused
  from the surfactant cert), TTN (delayed clearance, this cert's §5 Hansen-2008 anchor is the
  composite outcome TTN is typically the largest component of, disclosed as composite not
  TTN-specific, §8), meconium aspiration (named, not separately anchored this session, §8), PPHN
  (failed PVR fall, Steinhorn 2016 + Hooper 2015 coupling to the in-flight fetal-circulation cert).

**Confidence tier: mechanism-plausibility, geometrically self-consistent (Claim A) plus
multi-species-necessity-anchored (Claim B), not a novel wet-lab or simulated measurement.** Claim
A's inverse-Laplace calculation is arithmetic on two independently-anchored inputs (a measured
pressure range, a tension bracket already verified for an unrelated purpose in the sibling cert) —
as solid as that arithmetic and as those two anchors, with one disclosed assumption (bare/unspread
tension applies pre-aeration). Claim B rests on THREE independently live-verified primary-literature
necessity results across three species/methods (strongest form: a genetic knockout, a pharmacological
blockade, and a clinical cohort all converging on the same conclusion) — this is a comparatively
STRONG evidentiary base (real measured phenotypes, not simulated parameters), one tier above a
single-species or single-method literature anchor. Weakest links, disclosed: the specific Karlberg
first-breath P-V numbers were not independently re-extracted (bibliographic-only, pre-1975/1960s);
the viscous-resistance cross-check (Part 3) is explicitly illustrative, not decisive. Comparable to
this repo's own `MECHANISM_PULMONARY_SURFACTANT.md` precedent (geometrically-derived relation +
directly-measured calibration anchors, some pre-1975 bibliographic-only gaps disclosed).

## 10. Repro / Files

```
cd ~/projects/bodytwin
python3 scripts/msk/neonatal_transition.py
```

Pure Python arithmetic (no numpy/scipy dependency, no simulation, no subject data), runs in
<1 second. Writes `data/neonatal_transition/neonatal_transition_results.json`. No git operations
performed (isolation: never commit/push/add in this session regardless of repo state).

**Files, all created this session:**
- `scripts/msk/neonatal_transition.py` — the model (3 parts: opening-pressure inverse-Laplace
  geometry + 2 forced adversaries; fluid-clearance 3-species/method necessity gates; secondary
  viscous cross-check), deterministic (2 independent runs verified byte-identical).
- `data/neonatal_transition/neonatal_transition_results.json` — full machine-written evidence.
- `docs/MECHANISM_NEONATAL_TRANSITION.md` — this doc.
- `docs/MECHANISM_NEONATAL_TRANSITION_evidence.json` — curated citation + gate summary.
- Read read-only, NOT modified (isolation: touch only files created this session):
  `docs/MECHANISM_PULMONARY_SURFACTANT.md` + its `_evidence.json` (direct numeric reuse source),
  `docs/MECHANISM_PARTURITION_MYOMETRIUM.md`, `docs/MECHANISM_RESPIRATORY.md`,
  `docs/MECHANISM_CARDIAC_OUTPUT.md`, `docs/MECHANISM_MUCOCILIARY_CLEARANCE.md`,
  `docs/MECHANISM_HARDENED_CONVENTIONS.md`, `COORDINATOR.md`, `data/MECHANISM_ANCHOR_GRAPH.json`
  (checked for an existing fetal-circulation node — none found, §9).
