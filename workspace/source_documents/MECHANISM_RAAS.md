# MECHANISM RAAS — the renin-angiotensin-aldosterone system, the slow hormonal arm coupling renal Na handling, arterial pressure, and fluid volume (2026-07-22)

Resolves the RAAS-specific slice of this twin's cardiorenal-pressure axis: renin (JG cells) →
angiotensinogen → Ang I → [ACE] → Ang II → dual effector (AT1-vasoconstriction raising TPR +
zona-glomerulosa aldosterone raising Na/water retention → ECFV/plasma volume). Script:
`scripts/msk/raas.py`. Evidence: `data/raas/raas_results.json`.

**Headline result:** the pre-registered falsifier survives on BOTH legs, plus a forced,
non-tautological decorrelation test. **F1** (pressure-natriuresis / Guyton "infinite gain"): a real
population salt-intake range of **1210x** (Intersalt, n=10,079, 52 centres) is, in that same
primary source's own words, related to the *slope* of BP-with-age but explicitly **NOT** to
*median* BP across 48/52 centres — the qualitative and quantitative signature of a high-gain
regulator, not an assumption. **F2** (ACE-inhibitor/ARB blockade): standard-dose BP drop
**9.1/5.5 mmHg** (Law et al., 354 trials, 56,000 subjects) plus a directly-quoted **paradoxical
renin RISE + aldosterone FALL** on ACE inhibition (Atlas et al. 1984); decorrelated by a genuinely
different, modern, real-patient RCT running the intervention in **reverse** — exogenous Ang II
raises BP in vasodilatory shock (ATHOS-3, n=321, odds ratio **7.95**, 95% CI 4.76–13.3). The
**decorrelation** (primary hyperaldosteronism): a real n=251 cohort (PAPY study) shows the
aldosterone-renin ratio discriminates autonomous aldosteronomas from everyone else at **AUC
0.973** — and a machine-verified geometric argument shows *why*: an "autonomous" (renin-independent)
aldosterone source has an ARR-vs-renin log-log slope of **exactly −1**, a full order of magnitude
(**10x**) more divergent than the simplest renin-driven default (slope **0**), as renin is
suppressed. **17/17 gates PASS**, 2 independent runs byte-identical. **In-vivo-anchored** overall
(RAAS-blockade trials + real ARR diagnostic data); several mechanistic pieces are honestly
disclosed as textbook/order-of-magnitude tier, never laundered into a false precision.

## 0. Scope, stated up front — which graph nodes this is, and is NOT

This document builds a **standalone, quantitative, falsifiable model of the classical circulating
RAAS cascade itself** (mechanism + 2 falsifiers + 1 decorrelation test), and supplies it as a
**literature-anchored input** to five already-OPEN, much larger graph nodes in
`data/MECHANISM_ANCHOR_GRAPH.json` — it resolves none of them in full, and is not folded into the
graph this session (isolation discipline, `docs/MECHANISM_HARDENED_CONVENTIONS.md`; another mechanism
instance writes concurrently):

- **`INT-CARDIORENAL-PRESSURE-AXIS`** (OPEN): the full Guyton cardiorenal integrator over state
  `{MAP, ECFV, GFR, CO, TPR, PRA/AngII, aldosterone, ANP, BNP, renal_sympathetic_tone}`, including a
  **fitted neurogenic-reset weight term** (per its own Osborn-critique note). This doc supplies the
  `PRA/AngII` + `aldosterone` state variables' own mechanism and falsifier, but does **not** fit or
  calibrate the neurogenic-reset weight, ANP/BNP antagonism, or the full closed-loop ODE.
- **`ORG-ADRENAL-STRESS-HORMONES`** (OPEN): a 4-node adrenal state machine (zG-aldosterone,
  zF-cortisol, zR-DHEA, medulla-catecholamine). This doc supplies the zG (aldosterone) node's
  **upstream driver** (renin/Ang II) and downstream renal effector in depth, but builds none of the
  other 3 nodes.
- **`ORG-RENAL-FLUID-ELECTROLYTE`** (OPEN, SEED-DESIGN): its own claim names "(b) volume→RAAS/
  aldosterone loop gated by an eGFR-decline modifier" as one of two coupled control loops. This doc
  supplies that loop's mechanism/falsifier; the eGFR-decline gating and the (a) ADH/thirst loop are
  not built here.
- **`ORG-KIDNEY-NEPHRON`** (OPEN): a segmented nephron transport ODE that takes RAAS as an
  **input**. This doc supplies exactly that input's real magnitude/mechanism, not the segmental
  transport ODE itself.
- **`ENDO-ALDOSTERONE-HIDDEN-HTN`** (OPEN, SEED-DESIGN): a much larger claim about **hard clinical
  outcomes** (MI/stroke/HF/AF/CKD/mortality) in primary aldosteronism vs essential HTN, via
  claims/registry Cox/meta-analytic ascertainment. This doc's own ARR decorrelation (§7) is a
  **mechanism/diagnosis**-level result (real AUC data), explicitly **not** a hard-outcomes claim —
  does not resolve that node.
- **`XDOMAIN-ACE2-RAAS-COVID-RCT-NULL`** (OPEN): the ACE2/Ang-(1-7)/Mas **counter-regulatory** arm
  and RAAS-inhibitor COVID RCT null findings. **Explicitly out of scope** — this document models the
  classical ACE/Ang II/AT1 axis only, the same scope boundary the [8],[9] tissue-local-RAS
  citations already draw a line around (§9).

## 1. Geometric structure (derive from the geometry, not heuristics)

1. **Pressure-natriuresis "infinite gain" is a steady-state mass-balance identity, not a mystical
   property.** At steady state, Na-excretion-rate = Na-intake-rate. If the renal function curve
   (excretion vs MAP) has local slope `S = dExcretion/dMAP`, a step change in intake `dI` forces a
   new steady state at `dMAP_ss = dI/S` — exact, not fitted (verified analytically AND via numpy
   finite-difference, §5). As `S → 0`, `dMAP_ss` diverges for *any* nonzero `dI` — the forced
   adversary/void-floor showing the mechanism is load-bearing (§5).
2. **The aldosterone-renin ratio's diagnostic power is a causal-structure fact, not an arbitrary
   cutoff.** Model `aldo(renin) = k·renin^p`. Wherever aldosterone is renin/Ang-II-**driven**
   (`p>0`), `ARR = aldo/renin = k·renin^(p-1)` has log-log slope `(p-1)` — bounded, flat near `p=1`.
   The `p=0` case **is**, by construction, the definition of "autonomous" secretion — and for
   `p=0` the log-log slope of ARR-vs-renin is **exactly −1**, always, independent of `k`: as renin
   is suppressed by the resulting volume expansion (mechanism 1, acting through the same
   pressure-natriuresis loop), ARR diverges. A genuine, exponent-based structural discriminator
   (§7), not a rate/heuristic threshold.
3. **Two arms, one existing pressure law.** Ang II's vasoconstriction is a TPR-arm input;
   aldosterone's Na/water retention is a CO-arm (preload) input. Both multiply into this twin's own
   `MAP = CO × TPR` identity (`arterial_pressure.py`) — RAAS is a source term on both factors of an
   already-certified product, not a new law (§8).
4. **Timescale separation** (the 0-D pole/eigenvalue argument this repo's other MSK layers already
   use for the baroreflex/Windkessel): the baroreflex's own implied sympathetic loop delay
   (2.5–5.0 s, `baroreflex.py`) is the FAST pole; the renal-fluid/RAAS volume-accumulation arm is
   SLOW (hours-to-days, Guyton's [1] own words) — a real, >1000x, computed timescale separation
   (§8), not merely asserted.

## 2. Method, in one paragraph

Reads FOUR already-computed sibling JSONs read-only (`arterial_pressure_results.json`'s MAP/implied
TPR; `fluid_compartments_results.json`'s plasma/blood volume; `renal_filtration_results.json`'s
filtration fraction; `baroreflex_results.json`'s implied sympathetic delay) — no re-solve, no
OpenSim, no new subject data. Builds the RAAS cascade's geometric structure, then runs two
pre-registered falsifiers and one decorrelation test, each with an analytic identity machine-checked
against a numeric cross-check, a forced adversary/void-floor, and an external, non-tautological
real-data anchor — never a figure eyeballed, never a tautology gate.

## 3. Citations — every PMID verified LIVE this session (NCBI eutils + EuropePMC), not recalled

**Recall-drift measured explicitly, not swept under the rug**: of 9 PMIDs this session first tried
to recall from memory *before* searching, only **2 were correct** on the first guess (Guyton 1972,
Funder 2016) — a **7/9 = 77.8%** wrong-on-first-recall rate, at the high end of this repo's own
previously-measured ~62–67% citation-drift finding, and worse than `arterial_pressure.py`'s own
2/2=100% sample. Full log in `data/raas/raas_results.json["recall_drift_log"]`. Every number below
uses the live-corrected PMID.

| # | Citation | PMID/DOI | Tier / role |
|---|---|---|---|
| 1 | Guyton AC (1991). "Blood pressure control—special role of the kidneys and body fluids." *Science* 252(5014):1813-6. | **2063193** (live, full abstract) | Quoted: kidney system acts "within hours or days," is "the dominant method of establishing long-term pressure control"; BP "rarely deviates from normal by more than 10 to 15 percent." THE %-constancy anchor for F1. |
| 2 | Guyton AC, Coleman TG, Cowley AV Jr, et al (1972). "Arterial pressure regulation. Overriding dominance of the kidneys..." *Am J Med* 52(5):584-94. | **4337474** (identity live-verified; recalled guess CORRECT) | Bibliographic-only (no abstract indexed). THE original "infinite gain" dog-experiment paper. |
| 3 | Cowley AW Jr (1992). "Long-term control of arterial blood pressure." *Physiol Rev* 72(1):231-300. | **1731371** (live, full abstract; recalled guess 1731369 WRONG by 2 digits) | Quoted: baroreceptors "reset... cannot provide a sustained negative feedback signal" long-term, unlike the kidney. |
| 4 | Montani JP, Van Vliet BN (2009). "Understanding the contribution of Guyton's large circulatory model..." *Exp Physiol* 94(4):382-8. | **19286638** (live, full abstract; recalled guess WRONG) | Defends the Guyton-Coleman renal-dominant model; quoted on the pressure-natriuresis relationship. |
| 5 | Osborn JW, Averina VA, Fink GD (2009). "Current computational models do not reveal the importance of the nervous system..." *Exp Physiol* 94(4):389-96. | **19286640** (live, full abstract; companion point-counterpoint to [4]; recalled guess WRONG) | Quoted: the model "overestimates the importance of renal control." THE symmetric-QC decorrelated adversary, held OPEN. |
| 6 | Intersalt Cooperative Research Group (1988). *BMJ* 297(6644):319-28. | **3416162** (live, full abstract) | REAL n=10,079, 52 centres. Quoted: Na excretion 0.2–242 mmol/24h; 48/52 centres show NO significant relation to median BP. THE external F1 anchor. |
| 7 | Castrop H, Höcherl K, Kurtz A, et al (2010). "Physiology of kidney renin." *Physiol Rev* 90(2):607-73. | **20393195** (live, full abstract) | Renin-release cAMP(+)/Ca2+(−) signaling mechanism. |
| 8 | Paul M, Poyan Mehr A, Kreutz R (2006). "Physiology of local renin-angiotensin systems." *Physiol Rev* 86(3):747-803. | **16816138** (live, full abstract) | Tissue-local RAS, real. Held-OPEN symmetric-QC anchor (§9). |
| 9 | Yang T, Xu C (2017). "...Intrarenal Renin-Angiotensin System: An Update." *JASN* 28(4):1040-9. | **28255001** (live, full abstract; recalled guess 28904124 WRONG) | Quoted: intrarenal RAS is "a unique entity SEPARATE from systemic angiotensin II generation." |
| 10 | Law MR, Wald NJ, Morris JK, Jordan RE (2003). "Value of low dose combination treatment..." *BMJ* 326(7404):1427. | **12829555** (live, full abstract) | REAL 354-trial meta-analysis, 40,000 treated+16,000 placebo. Quoted: "9.1 mm Hg systolic and 5.5 mm Hg diastolic at standard dose." THE quantified BP-drop anchor. |
| 11 | Law MR, Morris JK, Wald NJ (2009). "...meta-analysis of 147 randomised trials." *BMJ* 338:b1665. | **19454737** (live, full abstract) | 464,000 participants; 5 drug classes "similarly effective" — corroborates [10]. |
| 12 | Atlas SA, Case DB, Sealey JE, Laragh JH, McKinstry DN (1979). *Hypertension* 1(3):274-80. | **399239** (identity verified 2 independent routes; NO abstract either place; recalled guess 396467 WRONG) | Title: captopril → sustained reduction in aldosterone, K+ retention, natriuresis. |
| 13 | Case DB, Atlas SA, Laragh JH, et al (1978). *Prog Cardiovasc Dis* 21(3):195-206. | **214819** (identity live-verified, no abstract) | Early captopril clinical-experience paper. |
| 14 | Atlas SA, Case DB, Yu ZY, Laragh JH (1984). "Hormonal and metabolic effects of ACE inhibitors..." *Am J Med* 77(2A):13-7. | **6089554** (live, full abstract) | Quoted verbatim: ACE inhibition "causes an INCREASE in plasma renin levels and a FALL in plasma and urine aldosterone." THE direct paradoxical-renin-rise confirmation. |
| 15 | Khanna A, English SW, Wang XS, et al; ATHOS-3 (2017). "Angiotensin II for the Treatment of Vasodilatory Shock." *NEJM* 377(5):419-30. | **28528561** (live, full abstract; recalled guess 28528550 WRONG by 11) | REAL RCT, n=321. BP-response 114/163 (69.9%) Ang II vs 37/158 (23.4%) placebo, OR 7.95 (95% CI 4.76–13.3). THE decorrelated converse confirmation. |
| 16 | Funder JW, Carey RM, Mantero F, et al (2016). Endocrine Society PA Guideline. *JCEM* 101(5):1889-1916. | **26934393** (live, full abstract; recalled guess CORRECT) | Establishes ARR as case-detection standard; abstract itself carries no numeric cutoff (disclosed — [17] supplies the numbers). |
| 17 | Rossi GP, Barisa M, Belfiore A, et al; PAPY (2010). *J Hypertens* 28(9):1892-9. | **20683340** (live, full abstract) | REAL n=251. PA prevalence 13.2% (6.4% APA+6.8% IHA); AUC 0.973 (PRA-ARR)/0.870 (DRA-ARR); DRA cutoff 27.3 ng/mIU. THE Conn's/ARR anchor. |
| 18 | Jacob G, Ertl AC, Shannon JR, et al (1998). "Effect of standing on neurohumoral responses..." *J Appl Physiol* 84(3):914-21. | **9480952** (live, full abstract) | REAL n=10. Standing: PV falls 13% (375±35mL); PRA rise correlated with aldosterone rise; kidney Na↓/K↑. Activation-arm confirmation. |
| 19 | Elijovich F, Weinberger MH, Anderson CA, et al (2016). AHA Scientific Statement. *Hypertension* 68(3):e7-e46. | **27443572** (identity live-verified; statement-format, no short abstract) | Salt sensitivity is a real, individually-variable phenotype spectrum. Held-OPEN anchor (§9). |

Reused, not re-verified (already live-verified in the sibling docs that produced them): MAP=93.33
mmHg / implied TPR=1278–1300 dyn·s·cm⁻⁵ (`arterial_pressure.py`, Razminia 2004, PMID 15558774);
plasma volume=3.033 L / blood volume=5.441 L (`fluid_compartments.py`); filtration fraction=20.1%
(`renal_filtration.py`, Davies & Shock 1950, PMID 15415454); implied sympathetic loop delay
2.5–5.0 s (`baroreflex.py`, deBoer 1987, PMID 3631301).

## 4. Cascade geometry

`renin (JG cells, rate-limiting) → angiotensinogen (hepatic, near-constant substrate excess) →
Ang I (inactive decapeptide) → [ACE, pulmonary/vascular endothelium] → Ang II (AT1 receptor) →
{vasoconstriction → raises TPR} + {zona glomerulosa → aldosterone → principal-cell ENaC/Na-K-ATPase
→ Na/water retention → raises ECFV/plasma volume}`.

**Renin-release triggers** (3, converging on the same cAMP(+)/Ca2+(−) JG-cell signaling machinery,
Castrop 2010 [7]): macula densa distal-NaCl sensing (low delivery → renin ↑), renal perfusion
pressure via the afferent-arteriole baroreceptor (low pressure → renin ↑), renal sympathetic tone
via β1 (SNS ↑ → renin ↑). **Short-loop negative feedback**: Ang II, via AT1 receptors *on the JG
cells themselves*, directly suppresses further renin release — the mechanism whose *interruption*
(ACE inhibition) produces the paradoxical renin rise (§6). Angiotensinogen's vast hepatic-secreted
excess relative to renin's own catalytic capacity (renin, not substrate, is rate-limiting) is
disclosed as textbook-tier, not individually live-pinned to a Km/substrate-concentration paper this
session.

## 5. Falsifier 1 — pressure-natriuresis / Guyton "infinite gain"

**Self-consistency** (analytic vs numpy finite-difference on the exact identity `dMAP_ss/dIntake =
1/S`): match to **<0.01%** diff, strictly monotonic — **PASS**.

**Forced adversary / void floor** (`S → 0`): raises a genuine division singularity under
`np.errstate(raise)` — confirms **no finite steady state exists for any nonzero intake change**
without a nonzero pressure-natriuresis slope. The mechanism is necessary, not decorative — **PASS**.

**Backward-implied minimum slope** (disclosed, NOT itself the gate — that would be circular):
combining Guyton's own 15% BP-constancy band with Intersalt's real intake range gives an implied
minimum `S ≈ 17.3 mmol/mmHg` (extreme 1210x range) or `≈14.3 mmol/mmHg` (typical 50–250 mmol/day
dietary sub-range) — reported as a derived quantity, not a fabricated absolute physiological
constant.

**External, non-tautological anchor** (the actual falsifier): Intersalt's [6] **own** quoted
finding — n=10,079, 52 centres, Na excretion range 0.2–242 mmol/24h (**1210x**), and "across the
other 48 centres sodium was significantly related to the SLOPE of blood pressure with age but NOT
to MEDIAN blood pressure" — **PASS**, a real, live-quoted, external primary finding, not asserted.

**Disclosed exception, not laundered**: 4 of the 52 centres (extreme-low-salt societies) *did* show
low BP and no age-related BP rise — a genuine, real tension (symmetric QC, held open per Elijovich
2016 [19], §9).

**Suppression-direction real data**: Jacob et al. 1998 [18] directly measured the **activation**
arm (n=10 healthy subjects, standing → 13% plasma-volume fall, PRA rise correlated with aldosterone
rise, kidney Na↓/K↑ response) — real, quantified, in-vivo. The task's own stated **suppression**
direction (salt loading suppresses renin/aldosterone) is the logically-entailed reverse of this same
mechanism; live-extracted numeric suppression values were not found this session (**disclosed gap**)
— but the entire Funder 2016 [16]/Rossi 2010 [17] confirmatory-testing paradigm (saline-suppression
test, captopril-challenge test) is *built* on the premise that normal aldosterone/renin *does*
suppress under volume/ACEi challenge, and primary aldosteronism is *defined* by failure to suppress
— real, live-verified, clinical-standard corroboration of the direction, if not this session's own
numeric extraction of it.

## 6. Falsifier 2 — ACE-inhibitor/ARB blockade + paradoxical renin rise + Ang II converse

**BP drop**: Law et al. 2003 [10], 354 trials, 40,000 treated + 16,000 placebo — standard-dose drop
**9.1 mmHg systolic / 5.5 mmHg diastolic**, inside the pre-registered [4,15]/[2,10] mmHg bands —
**PASS** (both). Corroborated by Law 2009 [11]'s independent 147-trial/464,000-participant finding
that the 5 drug classes are "similarly effective." Disclosed: 5-class average, ACEi not individually
isolated in this abstract.

**Void-floor on the band itself** (does it have real teeth?): the null/no-effect value (0 mmHg —
literally the trials' own placebo-arm reference, not a hypothetical) and an implausibly large value
(30 mmHg) both fall **outside** the pre-registered band — **PASS**, the band is not vacuously wide.

**Paradoxical renin rise**: Atlas et al. 1984 [14], quoted verbatim — ACE inhibition "causes an
increase in plasma renin levels and a fall in plasma and urine aldosterone, which can be sustained
for many years" — **direction PASS** (a direct primary-literature statement, not inferred).
Magnitude not quantified in this abstract (**disclosed gap**).

**Decorrelated converse** (a genuinely different, real, modern clinical scenario — *adding* Ang II
rather than removing it): Khanna et al. 2017, ATHOS-3 [15], n=321 vasodilatory/septic-shock
patients — BP-response endpoint reached by **114/163 (69.9%)** Ang II patients vs **37/158
(23.4%)** placebo, rate ratio **2.99x** (pre-registered ≥1.5x — **PASS**), odds ratio **7.95** (95%
CI 4.76–13.3, excludes 1 — **PASS**).

**Diverse instance-space tally** (machine-counted from the source list, not hardcoded): 3 distinct
study designs (pooled meta-analysis / comparative pharmacology study / single dedicated RCT), 2
distinct clinical contexts (chronic outpatient hypertension / acute ICU shock), 2 distinct
intervention directions (remove / add Ang II) — **3/3 axes** show ≥2 distinct values, meeting the
pre-registered bar — **PASS**. (An earlier draft of this table used free-text labels that
accidentally differed between legs meant to share a category, inflating the counts to 3/3/3; caught
by the script's own `len(set(...))` machine check before finalizing — fixed, not hidden.)

## 7. Decorrelation — Conn's/primary hyperaldosteronism, the ARR diagnostic dissociation

**The geometric argument** (§1.2): model `aldo(renin) = k·renin^p`. The **autonomous** case (`p=0`
— aldosterone secretion structurally independent of renin, i.e. primary hyperaldosteronism by
construction) gives `ARR = k/renin`, log-log slope **exactly −1**, machine-verified (analytic vs
`np.polyfit` regression, diff `<1e-6` for every tested `p` — **PASS**). The **simplest, non-cherry-
picked default** renin-driven case (`p=1`, aldosterone directly proportional to renin) gives `ARR =
k`, log-log slope **exactly 0** (flat). Evaluating both at a physiological renin range (suppressed
0.2 vs normal 2.0 ng/mL/h): the autonomous case's ARR ratio is **10.0x**; the `p=1` case's is
**1.0x** — a full order-of-magnitude gap, **PASS** against the pre-registered ≥3x bar.

**Self-caught sensitivity, disclosed not gated**: an earlier draft of this script gated on the
*maximum* divergence gap across a swept range of `p` values — a gameable design, since the ratio is
a *continuous* function of `p` with `p=0` as its supremum (never reached for any `p>0`, but
approached arbitrarily closely as `p→0+`). Fixed: the primary gate now uses only `p=1` (the obvious
default), and the full sweep (`p∈{0.5,0.8,1.0,1.3}`) is reported as a disclosed sensitivity — ratios
range 0.50x–3.16x, always strictly less divergent than the `p=0` case, but the *margin* is a
real, disclosed, un-resolved-by-the-toy-model question. The toy model's job is the **qualitative**
structural claim (p=0 is a genuine exponent threshold); the **quantitative** "how separated are real
patients" question is answered empirically by Rossi's real data, next.

**External, real, in-vivo anchor**: Rossi et al. 2010 [17], PAPY study, n=251 — primary
aldosteronism prevalence **13.2%** (within the pre-registered [5%,20%] band — **PASS**), ARR AUC for
identifying an aldosterone-producing adenoma **0.973** (PRA-based) / **0.870** (direct-renin-based),
both against the pre-registered >0.80 "excellent discrimination" bar — **PASS**, **PASS**. Real
patients separate on this ratio exactly as the causal-structure argument predicts.

## 8. Coupling — reused, not re-derived (4 sibling MSK JSONs, read-only)

| quantity | value | source |
|---|---:|---|
| MAP (resting) | 93.33 mmHg | `arterial_pressure.py` / `baroreflex.py` |
| Implied TPR (geo / subj) | 1300.2 / 1278.3 dyn·s·cm⁻⁵ | `arterial_pressure.py` |
| Plasma volume | 3.033 L | `fluid_compartments.py` |
| Blood volume | 5.441 L | `fluid_compartments.py` |
| Filtration fraction | 20.1% | `renal_filtration.py` |
| RBF | 1.077 L/min | `renal_filtration.py` |
| Implied sympathetic delay | 2.5–5.0 s | `baroreflex.py` |

**Two-arm re-entry into `MAP = CO × TPR`** (no new pressure law needed): Ang II's AT1-mediated
vasoconstriction is a source term on the **TPR** factor — the implied TPR reused here (1278–1300
dyn·s·cm⁻⁵) already sits *above* `arterial_pressure.py`'s own literal 900–1200 task-band, consistent
with (not proof of) nonzero resting Ang-II-driven tone. Aldosterone-driven Na/water retention is a
source term on plasma volume → venous return/preload → stroke volume (Frank-Starling,
`cardiac_output_geometric.py`, not re-derived) → cardiac output — the **CO** factor. **Efferent-
arteriole note**: Ang II preferentially constricts the efferent arteriole, raising filtration
fraction (this repo's own measured FF=20.1%) for a given RBF — the same mechanism by which ACEi/ARB
can drop GFR specifically in bilateral renal-artery stenosis, a well-known clinical consequence of
the identical mechanism modeled here.

**Timescale separation** (computed, not asserted): baroreflex fast pole 2.5–5.0 s vs a disclosed
order-of-magnitude RAAS/renal-fluid slow arm of 24–96 hours (Guyton's [1] own "within hours or
days") → ratio range **17,280x–138,240x** — **PASS** (>100x pre-registered bar), consistent with
this repo's own `ORG-ADRENAL-STRESS-HORMONES` graph node framing ("baroreflex fast loop handing off
to RAAS slow loop").

## 9. Symmetric QC — held OPEN, not resolved

- **Tissue-local vs circulating RAAS**: Paul et al. 2006 [8] and Yang & Xu 2017 [9] (quoted:
  intrarenal RAS is "a unique entity SEPARATE from systemic angiotensin II generation") both
  establish real, tissue-resident paracrine/autocrine RAS components. This document models the
  **circulating/systemic axis only** — tissue-local RAS's independent regulation is real,
  established, and **not modeled here**.
- **Pressure-natriuresis gain, genuinely debated in the primary literature itself**: Montani & Van
  Vliet 2009 [4] defend the Guyton-Coleman renal-dominant model; Osborn, Averina & Fink 2009 [5] — a
  real point-counterpoint in the **same journal issue** — argue it "overestimates the importance of
  renal control of body fluids and total blood volume" and that the sympathetic nervous system plays
  an important, independent long-term role. **Not resolved by this document.**
- **Salt-sensitivity is a spectrum, not a universal constant**: Elijovich, Weinberger et al. 2016
  [19] (AHA Scientific Statement) frame salt sensitivity of blood pressure as an actively-researched,
  individually-variable phenotype — the `S` parameter in §5's own identity is not one fixed number
  across people. Intersalt's [6] own 4-centre exception (§5) is a second, independent piece of the
  same open question.

## 10. Pre-registered gates — 17/17 PASS

```
f1_selfconsistency_analytic_numeric_deriv_match:    PASS (diff <0.01%)
f1_void_floor_s_zero_necessary:                     PASS (division singularity confirmed)
f1_external_intersalt_anchor:                       PASS (48/52 centres, quoted verbatim)
f2_sbp_drop_in_prereg_band:                          PASS (9.1 in [4,15] mmHg)
f2_dbp_drop_in_prereg_band:                          PASS (5.5 in [2,10] mmHg)
f2_bp_drop_band_void_floor_has_teeth:                PASS (0 and 30 mmHg both excluded)
f2_renin_rise_aldo_fall_direction:                   PASS (Atlas 1984, quoted verbatim)
f2_khanna_rate_ratio_prereg:                         PASS (2.99x >= 1.5x)
f2_khanna_or_ci_excludes_1:                          PASS (OR 7.95, CI 4.76-13.3)
f2_diverse_instance_space_ge_3_axes:                 PASS (3/3 axes >=2 distinct values)
dec_powerlaw_slope_analytic_numeric_match_all_p:     PASS (diff <1e-6, all 5 tested p)
dec_structural_discriminator_autonomous_more_divergent: PASS (10.0x vs 1.0x, p=1 default)
dec_rossi_auc_pra_pass:                              PASS (0.973 > 0.80)
dec_rossi_auc_dra_pass:                              PASS (0.870 > 0.80)
dec_rossi_prevalence_in_band:                        PASS (13.2% in [5,20]%)
coupling_no_missing_sibling_json:                    PASS (4/4 sibling JSONs found)
coupling_timescale_separation_fast_slow:             PASS (17,280x-138,240x > 100x)
--------------------------------------------------------------------------------
overall_pass:                                        PASS (17/17)
```

**Overall: PASS.** 2 independent runs produce byte-identical JSON (deterministic, pure
Python/numpy, no RNG).

## 11. Confidence tier

**In-vivo-anchored** for both falsifiers' external anchors (Intersalt's real 10,079-subject/
52-centre population data; Law's real 354-trial/56,000-subject meta-analysis; ATHOS-3's real
321-patient RCT) and for the decorrelation's external anchor (Rossi's real 251-patient PAPY-study
AUC/prevalence data). **Method-derived-and-geometrically-verified** for the two structural
arguments (the pressure-natriuresis steady-state identity; the ARR power-law exponent argument) —
machine-checked exactly, but their real-world parameter values (the true renal-function-curve slope
`S`; the true physiological aldosterone-vs-renin exponent `p`) are not independently live-pinned to
one primary number this session. **Bibliographic-only tier** for 3 pre-1985 citations with no
indexed abstract ([2],[12],[13]) — identity confirmed live, content not independently re-verified
beyond title/journal/author/year.

## 12. Honest gaps (disclosed, not hidden)

- **Salt-loading-suppresses-renin/aldosterone was not independently measured in that exact direction
  with live-extracted numeric values this session** (§5) — the activation-arm reverse (Jacob 1998)
  and the entire PA confirmatory-testing paradigm (Funder 2016/Rossi 2010) corroborate the direction,
  but no numeric supine/high-salt-diet PRA-and-aldosterone suppression values were extracted.
- **ACE-inhibitor-induced renin-rise magnitude is not quantified this session** (§6) — Atlas 1984's
  direction ("increase," "sustained for many years") is a direct quote; the fold-magnitude is not.
- **The renal-function-curve slope `S` and the aldosterone-renin exponent `p` are both, honestly,
  toy-model parameters** — the pressure-natriuresis and ARR sections demonstrate the *geometric
  mechanism* (an exact identity; an exponent threshold) rather than independently live-pinning the
  real physiological values of `S` or `p` to one primary measurement paper this session.
- **3 citations are bibliographic-only** ([2] Guyton 1972, [12] Atlas 1979, [13] Case 1978) — no
  abstract indexed/found live in either NCBI eutils or EuropePMC; identity confirmed by
  title/journal/author/year match only.
- **Funder 2016's [16] own abstract carries no numeric ARR cutoff** — the actual numbers used (§7)
  come from Rossi 2010 [17] instead, a stronger (in-vivo, not consensus-committee) evidentiary tier,
  disclosed rather than silently substituted.
- **Brunner & Laragh 1972 (renin-sodium profiling, PMID 4257928, identity-confirmed live) was
  checked this session but not ultimately used as a citation** — Guyton's own papers [1],[2] were
  found to cover the same "renin tracks sodium status" ground more directly relevant to this
  document's specific falsifiers; retained only in the recall-drift log (§3) for honest disclosure
  of the research process, not cited in the numbered list.
- **Tissue-local/intrarenal RAS is real and explicitly not modeled** (§9) — this document is
  circulating-axis-only throughout.
- **The pressure-natriuresis "gain" debate (Montani/Van Vliet vs Osborn) and the salt-sensitivity
  spectrum are both genuine, live, unresolved scientific questions** (§9) — reported, not resolved.
- **Single, population-level, resting-state model** — no per-subject trace, no dynamic/time-domain
  ODE simulation of the full loop (unlike, e.g., `arterial_pressure.py`'s Windkessel), no diurnal/
  postural/exercise regime beyond the single Jacob 1998 standing snapshot.
- **ACE2/Ang-(1-7)/Mas counter-regulatory axis not modeled** — explicit scope boundary vs the
  existing `XDOMAIN-ACE2-RAAS-COVID-RCT-NULL` graph node (§0).

## 13. Couples to (graph context, referenced not mutated)

- **Supplies an input to** `INT-CARDIORENAL-PRESSURE-AXIS`, `ORG-ADRENAL-STRESS-HORMONES`,
  `ORG-RENAL-FLUID-ELECTROLYTE`, `ORG-KIDNEY-NEPHRON` — see §0 for exactly which sub-claim of each
  is (and is not) resolved.
- **Couples to arterial-pressure** (`docs/MECHANISM_ARTERIAL_PRESSURE.md`): the TPR arm re-enters
  `MAP=CO*TPR` directly (§8), load-bearing (implied TPR reused, not re-derived).
- **Couples to fluid-compartments** (`docs/MECHANISM_FLUID_COMPARTMENTS.md`): the aldosterone arm's
  Na/water retention target is that document's own plasma-volume/ECFV state (§8), load-bearing.
- **Couples to renal-filtration** (`docs/MECHANISM_RENAL_FILTRATION.md`): the efferent-arteriole/
  filtration-fraction mechanism note (§8) directly reuses that document's own measured FF=20.1%.
- **Couples to baroreflex** (`docs/MECHANISM_BAROREFLEX.md`): the fast/slow timescale-separation
  argument (§8) directly reuses that document's own implied sympathetic delay bracket, load-bearing.
- **Not modified** (isolation discipline, concurrent-write safety): `data/MECHANISM_ANCHOR_GRAPH.json`
  — §0 explains this document's relationship to the five existing RAAS-adjacent nodes there.

## 14. Files

- `scripts/msk/raas.py` — full computation: cascade geometry, both falsifiers, the decorrelation's
  geometric power-law argument, coupling (4 sibling JSONs read read-only), all gates. Run with
  `.venv-msk/bin/python3 scripts/msk/raas.py` (<1s, no OpenSim, no GPU, pure Python/numpy).
- `data/raas/raas_results.json` — full evidence: every citation/PMID, every computed number
  (self-consistency checks, void-floors, implied slopes, power-law sweep, diverse-instance-space
  tally, coupling values), the recall-drift log, all 17 gates.
- Reused read-only (not modified): `data/arterial_pressure/arterial_pressure_results.json`,
  `data/fluid_compartments/fluid_compartments_results.json`,
  `data/renal_filtration/renal_filtration_results.json`, `data/baroreflex/baroreflex_results.json`.
- Not modified (isolation discipline): `data/MECHANISM_ANCHOR_GRAPH.json` — §0 explains why.

## 15. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/raas.py
```
Requires the 4 sibling JSONs listed in §14 to already exist (all read-only, never modified; run
their own scripts first if missing). No OpenSim call, no new subject-trial data, runs in under a
second, deterministic (verified: 2 independent runs produce byte-identical JSON). No git operations.

**Paths**: script `scripts/msk/raas.py`; evidence `data/raas/raas_results.json`; this doc
`docs/MECHANISM_RAAS.md`; upstream (read-only) inputs `data/arterial_pressure/arterial_pressure_results.json`,
`data/fluid_compartments/fluid_compartments_results.json`,
`data/renal_filtration/renal_filtration_results.json`, `data/baroreflex/baroreflex_results.json`.
