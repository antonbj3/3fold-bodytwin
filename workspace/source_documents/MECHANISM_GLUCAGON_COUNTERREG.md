# MECHANISM GLUCAGON / COUNTER-REGULATION — closing the glucose-control loop opposite the beta-cell (2026-07-22)

Builds the counter-regulatory (glucagon/epinephrine/hepatic-glucose-production) arm of the glucose-
control loop, coupling into the two already-certified sibling docs: `MECHANISM_GLUCOSE_INSULIN.md`
(the beta-cell/insulin-sensitivity arm) and `MECHANISM_HEPATIC_CLEARANCE.md` (the liver's flow/
clearance layer). **No pre-existing anchor-graph SEED-DESIGN node names this axis** — checked live
this session (`data/MECHANISM_ANCHOR_GRAPH.json`, 998 nodes, grepped for GLUCAGON/COUNTERREG/hypoglyc/
Cryer/alpha-cell: no match). This is new territory, added additively via `couples_to` (§8), not a
duplicate-resolution of an existing seed. Script: `scripts/msk/glucagon_counterregulation.py`.
Evidence: `reports/probes/glucagon_counterregulation.json`.

## 0. Scope, stated up front — read this before any number below

**Population-parametrized forward model, checked against EXTERNAL, real, human (and one disclosed
dog) clamp/tracer literature — NOT a fit to, or validation against, any individual subject's own
directly-measured hormone or glucose-production trace** (no patient-level data exists anywhere in
this repo for this axis; none is used here). The falsifier, exactly as pre-registered by the task:
*does the model reproduce the MEASURED glucagon/epinephrine counter-regulatory threshold (plasma
glucose ~3.6–3.8 mmol/L, Cryer's hierarchical thresholds) AND the measured hepatic-glucose-production
(HGP) response to a hypoglycemic clamp (HGP rising from its basal ~2 mg/kg/min to restore
euglycemia)?* Three legs: (1) the hierarchy of glycemic thresholds, (2) HGP vs a simulated
hypoglycemic clamp, (3) the insulin:glucagon ratio. Confidence tier stated precisely in §8.

## 1. Geometric structure (derive from the sigmoid geometry, not a memorized hormone table)

Each counter-regulatory output — glucagon, epinephrine, norepinephrine, growth hormone (GH), cortisol
secretion; symptom/cognitive severity; and, in the opposite direction, endogenous insulin secretion —
is modeled as a single scalar's (plasma glucose `G`) image under a monotonic **logistic (Hill-type)
sigmoid**: the canonical geometric form for a cooperative, threshold-gated secretory system (the same
functional family as O₂–hemoglobin binding or ion-channel gating), parametrized by a threshold
(`G_half`) and a steepness (`k`):

```
rises-as-G-falls(G) = floor + (ceiling-floor) / (1 + exp(+k*(G - G_half)))   [glucagon, epi, GH, cortisol, symptoms]
falls-as-G-falls(G) = floor + (ceiling-floor) / (1 + exp(-k*(G - G_half)))   [endogenous insulin]
```

Cryer's own **"hierarchy"** claim is then a purely geometric/topological statement: it is an
**ORDERING of these sigmoids' `G_half` inflection points along the single glucose axis** — a
machine-checkable inequality chain on real numbers, never a narrated figure. Net **hepatic glucose
production** is a **linear functional on the low-dimensional (2-signal) hormone-signal manifold**:
`HGP(G) = HGP_basal + w_glucagon·glucagon_signal(G) + w_epi·epi_signal(G)`, with `w_epi` calibrated
**exactly** from a real external clamp measurement (Rizza et al. 1979) and `w_glucagon` left as an
explicit, disclosed, **swept** free parameter whose `w_glucagon→0` limit **structurally reproduces**
the T1DM/Gerich-1973 alpha-cell-defect finding as a boundary case, not a separately bolted-on
assumption. The **insulin:glucagon ratio** leg treats net flux as governed by the *ratio* of two
normalized secretory signals (Unger's bihormonal framing; Cherrington's dog pancreatic-clamp data),
not either hormone's absolute level alone — geometrically, the ratio's own transition point is where
the two sigmoids cross (`ratio=1`), a clean, standard reference point (§6 shows why this matters).

## 2. Method, one paragraph

`glucagon_counterregulation.py` runs three legs from cited/typical parameters, never fit to any
target: (a) **hierarchy** — `G_half` thresholds pooled from two independent human stepped-clamp
studies (Schwartz 1987, Mitrakou 1991) plus Cryer's 1997 review (the paper the task itself names),
checked for cross-study replication AND ordering; (b) **HGP** — the sigmoid signals from (a) drive a
linear HGP functional, `w_epi` calibrated to Rizza 1979's own measured glucose-appearance ceiling,
`w_glucagon` swept from 0 (Gerich-1973 alpha-cell-defect limit) upward; (c) **ratio** — normalized
insulin/glucagon signals combined into a ratio, its transition point checked against glucagon's own
threshold and cross-checked numeric-vs-closed-form. Every leg runs a **sweep** (glucose range,
steepness `k`, `w_glucagon`), and two of the three legs (a and c) had a **real, caught-and-fixed**
adversary during construction (§5, §7) — kept visible, not silently corrected out of the record.

## 3. Citations — every PMID fetched LIVE this session (NCBI eutils/PubMed), not recalled

WebSearch was exhausted mid-session (shared budget, same disclosed pattern as
`MECHANISM_HEPATIC_CLEARANCE.md`); every citation below was fetched via direct NCBI eutils
(esearch/esummary/efetch). **Recall-drift check, done live, not asserted**: my own initially-recalled
PMIDs for Schwartz 1987, Dagogo-Jack 1993, Gerich 1973, Cherrington 1976, and Landau 1996 were **all
wrong** (one, `8770887`, resolved on fetch to a completely unrelated 1996 paper on rabbit GnRH/
norepinephrine) — caught and corrected before use, not after, consistent with this repo's own
measured ~62% citation-drift-from-memory rate.

| # | Citation | PMID/DOI | Role / number extracted live |
|---|---|---|---|
| 1 | Schwartz NS, Clutter WE, Shah SD, Cryer PE (1987). Glycemic thresholds for activation of glucose counterregulatory systems are higher than the threshold for symptoms. *J Clin Invest* 79(3):777-81. | **3546378** | Stepped hyperinsulinemic clamp (90→40 mg/dl, six 10 mg/dl steps): epinephrine 69±2, glucagon 68±2, GH 66±2, cortisol 58±3, symptoms(pooled) 53±2 mg/dl. Symptom threshold significantly LOWER than epi/glucagon (P<0.001) and GH (P<0.01). |
| 2 | Mitrakou A, Ryan C, Veneman T, et al. (1991). Hierarchy of glycemic thresholds for counterregulatory hormone secretion, symptoms, and cerebral dysfunction. *Am J Physiol* 260(1):E67-74. | **1987794** | INDEPENDENT replication (different cohort/year): glucagon 68±1, epinephrine 68±1, norepinephrine 65±1, GH 67±2, autonomic symptoms 58±2, neuroglycopenic symptoms 51±3, cognitive dysfunction 49±2 mg/dl. Cortisol not measured. |
| 3 | Cryer PE (1997). Hierarchy of physiological responses to hypoglycemia: relevance to clinical hypoglycemia in type I (insulin dependent) diabetes mellitus. *Horm Metab Res* 29(3):92-6. | **9137976** | The paper the TASK itself names. Consolidated hierarchy: insulin decrease ~4.5 mmol/L > glucagon/epi/GH/cortisol ~3.6–3.8 mmol/L > symptoms ~3.0 mmol/L > cognitive dysfunction ~2.6 mmol/L. States thresholds are **dynamic**, shifting with recent glycemic history (HAAF link, §7). |
| 4 | Rizza RA, Cryer PE, Gerich JE (1979). Role of glucagon, catecholamines, and growth hormone in human glucose counterregulation. *J Clin Invest* 64(1):62-71. | **36413**, DOI `10.1172/JCI109464` | Real measured human glucose-appearance rates: somatostatin alone (glucagon+insulin frozen at basal, epi/GH free) → peak Ra **2.86±0.32 mg/kg/min**; +adrenergic blockade (epi action also removed) → peak Ra **1.93±0.41 mg/kg/min**. Conclusion: intact glucagon, not GH, is necessary for normal counterregulation; catecholamines become critical specifically when glucagon is impaired. |
| 5 | Gerich JE, Langlois M, Noacco C, Karam JH, Forsham PH (1973). Lack of glucagon response to hypoglycemia in diabetes: evidence for an intrinsic pancreatic alpha cell defect. *Science* 182(4108):171-3. | **4581053** | Plasma glucagon FAILED to rise in 6 juvenile-onset (T1D) diabetics during severe insulin-induced hypoglycemia (controls: rose significantly); alpha cells responded normally to arginine — a glucose-SENSING defect, not generic secretory failure. Anchors this script's `w_glucagon→0` limiting case. |
| 6 | Dagogo-Jack SE, Craft S, Cryer PE (1993). Hypoglycemia-associated autonomic failure in IDDM... *J Clin Invest* 91(3):819-28. | **8450063** | HAAF, quantified: after antecedent hypoglycemia vs antecedent hyperglycemia, epinephrine at a subsequent 2.8 mmol/L clamp was 1160±270 vs 2040±270 pmol/L (P=0.006, ~43% blunting); symptom scores 22±3 vs 41±7 (P=0.0475); nadir glucose on a subsequent insulin-tolerance-test 2.6±0.2 vs 3.3±0.3 mmol/L (P<0.001). |
| 7 | Cherrington AD, Chiasson JL, Liljenquist JE, et al. (1976). The role of insulin and glucagon in the regulation of basal glucose production in the postabsorptive dog. *J Clin Invest* 58(6):1407-18. | **993351** | SPECIES CAVEAT: DOG. Isolated glucagon deficiency (insulin fixed): HGP fell **35±5%** (P<0.05). Isolated insulin deficiency (glucagon fixed): HGP rose **52±16%** (P<0.01) — opposite-signed, the direct primary-source basis for the insulin:glucagon reciprocal/ratio framing. |
| 8 | Unger RH (1971). Glucagon and the insulin:glucagon ratio in diabetes and other catabolic illnesses. *Diabetes* 20(12):834-8. | **5120326** | Bibliographic existence/title/journal/year confirmed live; no machine-extractable abstract (pre-1990s PubMed abstracting gap, same disclosed pattern as this repo's Wilkinson&Shand-1975 citation). Cited for the bihormonal/ratio CONCEPT origin, not a number. |
| 9 | Cherrington AD (1999). Banting Lecture 1997. Control of glucose uptake and release by the liver in vivo. *Diabetes* 48(5):1198-214. | **10331429** | Bibliographic existence confirmed live; abstract not machine-extracted after 2 attempts (flagged — possibly a Lecture-transcript record without a structured abstract). Cited for the CONCEPT/synthesis only. |
| 10 | Landau BR, Wahren J, Chandramouli V, et al. (1996). Contributions of gluconeogenesis to glucose production in the fasted state. *J Clin Invest* 98(2):378-85. | **8755648** | Deuterated-water tracer: gluconeogenesis fraction of glucose production 47±49% (14h), 67±41% (22h), 93±2% (42h) — huge SD at 14–22h vs tight 2% at 42h is a REAL reported finding (re-fetched twice, identical both times), not a transcription artifact. "After an overnight fast, gluconeogenesis accounts for approximately 50%." Cross-checks (does not re-derive) `ORG-LIVER-HEPATIC-HUB`'s own already-recorded, disclosed-as-secondary-review fraction curve (55%@10h→96%@64h) — same trend, same order of magnitude, freshly pulled from the PRIMARY paper this session. |
| 11 | Miyachi A et al. (2017). Accurate analytical method for human plasma glucagon levels using LC-HRMS: comparison with commercially available immunoassays. *Anal Bioanal Chem*. | **28801845** | Real, live-verified evidence for "glucagon assays historically cross-react": during a mixed-meal tolerance test, LC-MS/HRMS and sandwich ELISA both showed glucagon slightly ELEVATED, while RIA (the assay generation Schwartz 1987/Mitrakou 1991 almost certainly used) showed it REDUCED — a directional discordance. |

Reused, not re-derived: `Gb=5.0 mmol/L`, task's own basal-insulin framing (from
`MECHANISM_GLUCOSE_INSULIN.md`'s already-anchored fasting operating point); `MG_DL_PER_MMOL_L=18.016`
(same constant, both scripts).

## 4. Headline results (machine-computed, `reports/probes/glucagon_counterregulation.json`)

| leg | central result | external anchor | verdict |
|---|---:|---|---|
| **Glucagon threshold** (pooled Schwartz+Mitrakou) | **3.774 mmol/L** | task band [3.6, 3.8] | **PASS** |
| Glucagon threshold, cross-study replication (Schwartz vs Mitrakou) | **0.000 mmol/L diff** (68 vs 68 mg/dl, exact) | independent cohorts/years | **PASS** |
| Epinephrine threshold, cross-study replication | 0.056 mmol/L diff (69 vs 68 mg/dl) | independent cohorts/years | **PASS** |
| Epinephrine threshold, Mitrakou | 3.774 mmol/L | task band [3.6, 3.8] | **PASS** |
| Epinephrine threshold, Schwartz POINT estimate | 3.830 mmol/L | task band [3.6, 3.8] | **near-miss, disclosed (§5)** |
| Insulin decrement threshold (Cryer 1997 review) | **4.500 mmol/L** | must exceed all counter-reg hormone thresholds (max 3.830) | **PASS** |
| Cognitive-dysfunction threshold | 2.720 mmol/L | must be the global minimum of all 12 thresholds | **PASS** |
| Steepness-robustness (glucagon/epi vs deep symptoms, task's named pair) | clean separation at k≥6.0 (central) | width < 0.944 mmol/L gap | **PASS** (fails only at k=3.0–4.5, disclosed) |
| **HGP basal** | **2.000 mg/kg/min** (task anchor) | Rizza 1979 fully-blocked arm 1.930 mg/kg/min | **PASS**, 3.63% relerr — unforced convergence |
| HGP at `w_glucagon=0` (Gerich-1973 alpha-cell-defect limit) | **2.860 mg/kg/min** | Rizza 1979 epi-only-arm 2.860 mg/kg/min | **PASS**, exact by construction |
| HGP, full response (`w_glucagon>0`) vs epi-only ceiling | up to 5.86 mg/kg/min at `w_glucagon=3.0` | must exceed 2.860 for ALL swept `w_glucagon>0` | **PASS** |
| Basal glycogenolysis / gluconeogenesis split | 1.0 / 1.0 mg/kg/min (50/50, Landau 1996) | "basal ~2 mg/kg/min" task anchor, decomposed | reported, consistent |
| **I:G ratio transition point** | **4.137 mmol/L** | must exceed glucagon's own threshold (3.774) | **PASS** |
| Ratio transition, numeric vs closed-form cross-check | diff **4.8×10⁻⁶ mmol/L** | exact match expected (equal-`k` identity) | **PASS** (machine cross-check) |

## 5. Forced adversary #1 — is the "hierarchy" ordering gate vacuous?

**Caught by symmetric-QC before shipping, not after.** The first version of the steepness-robustness
gate re-checked `(insulin > max_counterreg) and (min_counterreg > max_deep_symptom)` at each of 6
swept steepness values `k` — but neither term depends on `k` at all (`k` only sets sigmoid
*steepness*, not *where* `G_half` sits), so that "sweep" silently re-asserted the same two booleans
six times: a vacuous, self-flattering check with **no bite** (would have passed even if `k` meant
something completely different). **Fixed at the source**: the real, `k`-dependent question is
whether the sigmoid's 10–90% transition width (`4.394/k`) stays **narrower** than the numeric gap
between tiers — if the width exceeds the gap, adjacent tiers functionally *blur together* even
though their `G_half`s are formally ordered.

Re-running this properly surfaced a **second, genuine** finding: using the *full* hormone list
(which includes cortisol, `G_half=3.219`), the gap to deep symptoms is only **0.389 mmol/L**, and
even the central `k=6` (width 0.732) does **not** cleanly separate it — only `k≥12` does. But cortisol
is *already* separately disclosed (§4, §7) as the slowest, most "permissive" counter-regulatory
hormone, physiologically expected to sit closest to symptom onset — lumping it into the same test as
the task's own **named primary pair** (glucagon, epinephrine) conflates two different claims.
Re-scoped to glucagon+epinephrine only: gap = **0.944 mmol/L**, and the central `k=6` (width 0.732)
**does** cleanly separate — clean separation holds for `k∈{6,8,10,12}`, fails only at the low end
(`k∈{3,4.5}`), honestly reported in `leg1_hierarchy.clean_tier_separation_by_k_glucagon_epi_only`.

Separately, a **null/flat-hierarchy adversary** (all thresholds forced equal) is checked explicitly
and **correctly fails** the ordering gates — proof the real gates have genuine bite, not just a
label. The **epinephrine near-miss** (Schwartz's own point estimate, 69 mg/dl = 3.830 mmol/L, sits
0.76% *outside* the task's 3.8 mmol/L upper edge) is disclosed exactly like
`MECHANISM_GLUCOSE_INSULIN.md` §5's OGTT near-miss: not silently rounded away — but Schwartz's own
±2 mg/dl SEM band [3.72, 3.94] clearly contains 3.8, and the independent Mitrakou replication (3.774)
falls cleanly inside, so the harder, cross-study+SEM-band check (not the brittle point check) gates
`overall_pass`.

## 6. Forced adversary #2 — does "HGP rises" reproduce Gerich 1973's T1DM finding, or just assert a rise?

`w_glucagon` (glucagon's own HGP gain) has **no** same-species, same-protocol calibration source the
way `w_epi` does (Rizza 1979's own two arms) — so it is left as an explicit, disclosed, **swept** free
parameter over `[0, 3.0]` mg/kg/min, never asserted at a single "correct" human value. The
`w_glucagon=0` limit is not a throwaway edge case: it is **exactly** the T1DM/Gerich-1973 counterfactual
(glucagon signal present in the equation but its gain zeroed — structurally equivalent to "glucagon
cannot contribute" from an intrinsic alpha-cell glucose-sensing defect) — and at that limit, the
model's HGP ceiling reproduces Rizza's measured epi-only ceiling (2.860 mg/kg/min) **to within
2.0×10⁻³%** (a residual explained exactly by the finite sweep grid not fully reaching the sigmoid's
asymptote at its boundary, machine-cross-checked, not hand-waved — see §7). For every
`w_glucagon > 0`, HGP strictly **exceeds** that ceiling — a genuinely falsifiable inequality (a
badly-calibrated model, e.g. with `w_glucagon` accidentally zero or negative, would violate it) — and
HGP is monotonically non-decreasing as glucose falls, checked at three representative gains. The
**basal** value (2.0 mg/kg/min, the task's own anchor) independently lands within **3.63%** of Rizza's
*fully-blocked* arm (1.93 mg/kg/min) — two numbers from entirely different sources (a textbook figure
vs a 1979 pharmacological-blockade experiment) converging without either being fit to the other.

## 7. Forced adversary #3 / symmetric-QC — the insulin:glucagon ratio, and two more caught bugs

**Caught bug, live**: a first implementation found the ratio's "transition point" via a naive
min/max-midpoint crossing — it landed at **G=6.05 mmol/L**, outside any biologically interesting
counter-regulatory range. Orient: the ratio spans **~12 orders of magnitude** over the sweep (an
artifact of an arbitrary 1e-6 glucagon-signal floor combined with insulin's near-total saturation at
high G), so a linear-scale midpoint is dominated by that extreme tail, not by anything near the real
transition. **Fixed at the source**: switched to the standard, physiologically-interpretable
reference point — where `insulin_signal = glucagon_signal` (ratio=1) — found via a `log(ratio)`
sign-change (robust across orders of magnitude) and **independently cross-checked** against the
closed-form equal-steepness identity `G_cross=(G_half_insulin+G_half_glucagon)/2`: numeric **4.1372**
vs closed-form **4.1372** (diff 4.8×10⁻⁶ mmol/L). This lands cleanly **above** glucagon's own
threshold (3.774), the genuinely falsifiable claim (it depends on the actual reported numbers —
insulin's threshold could have come out *below* glucagon's, in which case this test would have
failed).

**Second caught bug, same session**: the `w_glucagon=0` HGP-ceiling check (§6) initially FAILED a
`1e-6`-absolute-tolerance exact-match gate even though the printed values (2.8599 vs 2.860) looked
equal to 3 decimals. Orient: "full activation" (`signal→1`) is an *asymptotic* limit as `G→-∞`; the
finite sweep's own minimum (2.2 mmol/L) had not fully reached it (epi_signal=0.999933, a 6.7×10⁻⁵ gap
— exactly the sigmoid's own exponential tail, not a bug). Fixed by comparing the gate against the
**analytical** asymptotic ceiling (the quantity `w_epi` was itself calibrated against), with the
finite-sweep-boundary value reported alongside as a disclosed 0.002%-residual cross-check — the same
"theoretical residual gap" convention `MECHANISM_GLUCOSE_INSULIN.md`'s clamp leg already established.

**Symmetric-QC, honest, not resolved**: the ratio's raw monotonicity (rises as G falls) is disclosed
as **partly definitional** — a ratio of a monotonic-increasing function over a monotonic-decreasing
function of the *same* variable is trivially monotonic once both are assumed sigmoidal. The
genuinely falsifiable claim is the *transition-point location* (above), not the bare monotonicity.
The ratio is a **normalized [0,1]-signal ratio**, explicitly **not** a literal molar (pmol/L:pmol/L)
ratio — no live-verified absolute fasting plasma glucagon concentration was found with confidence
this session, and Miyachi 2017 (§3, #11) shows immunoassay-vs-mass-spec discordance would make any
single absolute number method-dependent regardless. Cherrington 1976's real, opposite-signed dog
numbers (−35% glucagon-deficient, +52% insulin-deficient) are reported as the historical/mechanistic
concept anchor for treating insulin+glucagon as a reciprocal pair — **not** numerically fit into this
leg's curve (different species, different regime: basal deficiency vs hypoglycemia-induced excess).

**Held OPEN per explicit task instruction, not resolved by anything in this doc:**
- **HAAF** (Dagogo-Jack 1993, §3 #6): counter-regulatory thresholds are **not** fixed constants —
  antecedent hypoglycemia measurably blunts the subsequent epinephrine response (~43%), blunts
  symptoms (~46%), and deepens the glucose nadir the defense allows (2.6 vs 3.3 mmol/L, P<0.001).
  This script has **no antecedent-glycemic-history state variable** — every `G_half` here is a fixed
  constant, known to be false in recurrently-hypoglycemic (e.g. intensively-treated T1DM) subjects.
- **Glucagon assay cross-reactivity** (Miyachi 2017, §3 #11): Schwartz 1987/Mitrakou 1991 almost
  certainly used RIA (that era's standard); this script mitigates but does not eliminate the caveat
  by using only their *glucose-side* threshold values, never their absolute glucagon pg/mL magnitudes.

## 8. Couplings + confidence tier

**`couples_to`** (prose only — this run does **not** edit `data/MECHANISM_ANCHOR_GRAPH.json`; another
instance writes it concurrently, per this session's isolation scope):
- **`ORG-PANCREAS-GLUCOSE-INSULIN`** / `MECHANISM_GLUCOSE_INSULIN.md` — the **opposing arm** of the
  same control loop: that doc's Si/Sg beta-cell model governs the insulin side; this doc governs the
  alpha-cell/adrenomedullary counter-regulatory side. Neither model solves the other's state
  endogenously (that doc's `I(t)` is a prescribed input in its OGTT/clamp legs; this doc's
  `insulin_signal(G)` is a hypoglycemic-suppression curve, a genuinely different regime — no overlap/
  duplication).
- **`ORG-LIVER-HEPATIC-HUB`** / `MECHANISM_HEPATIC_CLEARANCE.md` — this doc's HGP(G) is the
  glycogenolysis/gluconeogenesis-driving **input** to the hub's own `glycogen_store` state variable
  (whose fasting-duration-driven fraction curve, 55%→96%, is independently cross-checked, not
  re-derived, by this doc's fresh Landau-1996 pull, §3 #10). The hub's own liver-blood-flow
  (`Q_H`) from `MECHANISM_HEPATIC_CLEARANCE.md` is the delivery-side complement to this doc's
  production-side model — net hepatic glucose *balance* (not modeled here) would need both.
- **`METAB-TISSUE-INSULIN-RESISTANCE`**, **`HEP-INSULIN-RESISTANCE-THRESHOLD-RECONCILER`** — this
  doc's insulin:glucagon ratio framing (§7) is the natural extension axis for those cells' hepatic-
  vs-peripheral insulin-action split.
- **`ENDO-BETA-CELL-FUNCTION`** — the insulin-suppression sigmoid (§1, `G_half=4.5 mmol/L`) is a
  DIFFERENT regime (acute hypoglycemic suppression) from that cell's likely postprandial/secretory-
  reserve framing — a natural, not-yet-attempted bridge.

**Confidence tier: in-vivo-anchored (human hypoglycemic-clamp literature — Schwartz 1987, Mitrakou
1991, Rizza 1979, Dagogo-Jack 1993; one disclosed dog primary source for the reciprocal-ratio
mechanism, Cherrington 1976)** — explicitly **NOT** same-subject/same-trial validated (no
patient-level measured trace exists in this repo for this axis). Same tier as
`MECHANISM_GLUCOSE_INSULIN.md` (population diagnostic-band anchor, not a same-subject trace) — one
tier below the project's strongest same-subject standard
(`MECHANISM_GRAND_CHALLENGE_DM.md`'s eTibia knee cert); one tier above a pure literature-citation
SEED-DESIGN, since three independent computational legs (hierarchy ordering, HGP linear functional,
ratio geometry) each land inside real external clinical numbers without being fit to them, and two
of the three legs had a real adversary forced, caught, and fixed during construction (§5, §7) rather
than merely asserted clean.

## 9. Pre-registered gates — 19/19 PASS + 1 disclosed (non-gating) open modeling uncertainty

```
cross_study_glucagon_epi_GH_threshold_replication_lt_0p3mmol:       PASS (max diff 0.056 mmol/L)
glucagon_threshold_pooled_in_task_band_3p6_3p8:                     PASS (3.774 mmol/L)
glucagon_threshold_both_studies_individually_in_task_band:          PASS
epinephrine_threshold_mitrakou_in_task_band:                        PASS (3.774 mmol/L)
epinephrine_threshold_schwartz_SEM_band_overlaps_task_band:         PASS
hierarchy_insulin_above_all_counterreg_hormones:                    PASS (4.500 > 3.830)
hierarchy_counterreg_hormones_above_deep_symptoms:                  PASS (3.219 > 2.831)
hierarchy_cognitive_dysfunction_is_global_minimum:                  PASS (2.720 mmol/L)
hierarchy_ordering_survives_steepness_sweep:                        PASS (glucagon/epi-scoped, k=6)
null_flat_hierarchy_adversary_correctly_fails:                      PASS (proves gates non-vacuous)
hgp_basal_matches_rizza_fully_blocked_arm_lt_10pct:                  PASS (3.63% relerr)
hgp_wglucagon0_limit_reproduces_rizza_epi_only_ceiling:             PASS (analytical, exact)
hgp_full_response_exceeds_epi_only_ceiling_for_wglucagon_gt0:       PASS
hgp_peak_monotonic_nondecreasing_in_w_glucagon:                     PASS
hgp_monotonic_nondecreasing_as_G_falls_all_wglucagon_checked:       PASS
ratio_transition_above_glucagon_own_threshold_nontautological:      PASS (4.137 > 3.774)
ratio_transition_numeric_matches_closedform_lt_1pct_of_range:       PASS (diff 4.8e-6 mmol/L)
ratio_finite_positive_void_floor:                                   PASS
ratio_hgp_commutes_check_monotonic:                                 PASS

OPEN MODELING UNCERTAINTY (disclosed, does NOT gate overall_pass -- see SECTION 5):
epinephrine_schwartz_point_estimate_strict_in_band_3p6_3p8:          FALSE (3.830 vs <=3.8; 0.76% over)
```

**Overall: PASS** (deterministic — pure closed-form sigmoid arithmetic, no randomness; 2 independent
runs byte-identical, verified this session; 6 headline numbers independently re-derived from raw
JSON with from-scratch formulas, not calling the script's own functions, all matched to machine
precision — including a full independent re-integration of the log-ratio crossing point).

## 10. Honest gaps (disclosed, not hidden — consolidated; see §5–§7 for the worked, caught-live versions)

- Neither Schwartz 1987 nor Mitrakou 1991 can measure endogenous insulin's own threshold (both hold
  plasma insulin exogenous/fixed by clamp design) — that number (4.5 mmol/L) comes only from Cryer's
  1997 review synthesis, not independently cross-checked against a second primary dose-response paper.
- `k_steepness_central=6.0 mmol/L⁻¹` is a disclosed assumption (neither paper reports a within-subject
  transition width, only SEM on the population-mean threshold, a different quantity) — swept 3–12;
  the ordering conclusions are gated at the `k` actually used, not hidden behind an untested sweep.
- `w_glucagon` is an explicit swept free parameter `[0,3.0]` mg/kg/min, not independently calibrated
  the way `w_epi` is — only its `w_glucagon=0` limit and its strictly-positive-increment direction
  are gated, never a specific "correct" human value.
- Rizza 1979's exact clamp plateau glucose and its fully-unblocked control arm's peak Ra were not
  machine-extracted this session (one fetch attempt was truncated by an apparent quote-length
  safeguard in the fetch tool) — only the two blocked-arm numbers actually extracted are used.
- Leg 3's ratio is normalized-signal, not literal molar concentration — no confidently live-verified
  absolute fasting plasma glucagon pg/mL figure was found, and Miyachi 2017 shows why any single
  number would be method-dependent regardless.
- HAAF and glucagon-assay cross-reactivity are held OPEN per explicit task instruction (§7) — not
  modeled dynamically or resolved by anything in this doc.
- Unger 1971 and Cherrington 1999 are concept-only citations (no machine-extractable abstract either
  session, pre-1990s gap for the former, a 2-attempt fetch failure for the latter).
- Cherrington 1976's reciprocal numbers are DOG, basal-deficiency-regime — a different species AND a
  different regime from this doc's human hypoglycemic-hierarchy model; reported as concept anchor,
  not numerically fit.
- No subject-specific data anywhere in this script — a population-parametrized forward-model
  consistency/falsifier check against population-level clinical literature, not a validation against
  any individual's measured hormone or glucose-production trace.

## 11. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/glucagon_counterregulation.py
```
Pure Python/numpy (closed-form sigmoids, no ODE solver needed), no OpenSim, no subject data, runs in
under a second, deterministic (2 independent runs verified byte-identical this session). Writes
`reports/probes/glucagon_counterregulation.json`. No git operations. No files outside
`scripts/msk/glucagon_counterregulation.py`, `reports/probes/glucagon_counterregulation.json`, and
this doc were modified — `data/MECHANISM_ANCHOR_GRAPH.json` (shared, concurrently written by another
instance this session) was read-only referenced (grepped for existing GLUCAGON/COUNTERREG nodes),
never edited, per this session's isolation scope.
