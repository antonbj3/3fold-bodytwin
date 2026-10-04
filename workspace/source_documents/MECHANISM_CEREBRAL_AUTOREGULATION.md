# MECHANISM CEREBRAL BLOOD FLOW AUTOREGULATION — brain-perfusion layer, coupling arterial pressure (CPP) and acid-base (CO2 reactivity) (2026-07-22)

Extends the cardiovascular chain (`docs/MECHANISM_CARDIAC_OUTPUT.md` → `docs/MECHANISM_ARTERIAL_PRESSURE.md`,
`MAP=CO*TPR`) and the acid-base chain (`docs/MECHANISM_ACID_BASE_CO2.md`, PaCO2) into **cerebral blood
flow (CBF) regulation**: resting CBF (~50 mL/100g/min whole-brain, ~80 grey, ~20 white), the
pressure-autoregulation plateau (CBF held roughly constant over CPP~60-150mmHg via myogenic+metabolic+
neurogenic mechanisms), `CPP = MAP - ICP`, and CO2 reactivity (CBF rises with PaCO2). Script:
`scripts/msk/cerebral_autoregulation.py`. Evidence:
`data/cerebral_autoregulation/cerebral_autoregulation_results.json`.

**Does NOT re-solve any upstream layer.** No OpenSim, no `.osim`, no `.sto` file. Reads TWO already-
computed JSON outputs read-only — `arterial_pressure_results.json` (this twin's own MAP, both the
classic-formula and Windkessel-true-mean estimates) and `acid_base_co2_results.json` (this twin's own
acute-respiratory PaCO2 perturbation scenario, 40→50mmHg) — and does pure-Python/numpy closed-form
arithmetic on top.

## 0. Scope, stated up front — read this before any number below

This is a **lumped, 0-D, STATIC-autoregulation-only** brain-perfusion layer: one scalar CBF state, a
closed-form pressure-flow curve, a closed-form CO2-reactivity curve, no vascular tree, no regional
blood-flow map beyond a grey/white split, **no dynamic autoregulation** (the transient response to a
sudden BP change — thigh-cuff release, Tiecks ARI 0-9 — is a *different*, faster-timescale mechanism,
explicitly not modeled here; Rangel-Castilla et al 2008 [15] state directly that "assessment of cerebral
autoregulation is best achieved with dynamic autoregulation methods," held OPEN per the task's own
framing, §13). **The classic Lassen (1959) [1] "flat 60-150mmHg" plateau is the textbook idealization
this document explicitly tests against modern literature, not an assumed truth** — Willie et al 2014 [5]
state, as one of their own review's four central theses, that "cerebral autoregulation does **not**
maintain constant perfusion through a mean arterial pressure range of 60-150mmHg"; Numan et al 2014 [4]
quantify a real, nonzero "leak" in the plateau. Both tensions are held OPEN below (§8), not resolved
either way. No subject-specific CBF/ICP measurement exists for subject2 — every external anchor is
population-level literature, the same first-step scope every other `MECHANISM_*` layer discloses for
itself.

## 1. Geometric structure (derive from the geometry, not heuristics)

Pressure autoregulation is modeled as an actively-**regulated vascular resistance** `R(P)` that tracks
`R_ideal(P) = P/F0` (the resistance value that would hold flow EXACTLY at `F0` at every pressure `P`)
but is physically bounded, `R ∈ [R_min, R_max]` (a vessel cannot dilate/constrict without limit).
Setting `R_min = LLA/F0`, `R_max = UBP/F0` (continuity at the knees, by construction) and substituting
into `F(P) = P/R(P)` gives, with **no free curve-fitting parameter**:

```
F(P) = F0                for LLA <= P <= UBP     (the idealized flat "plateau")
F(P) = F0 * (P/LLA)      for P < LLA               (pressure-passive line through the origin)
F(P) = F0 * (P/UBP)      for P > UBP               (pressure-passive line through the origin)
```

a clipped-line / saturating-gain-scheduled-controller construction — the same geometric family
`arterial_pressure.py` already uses for its own R/C void-floor clamps. The **modern correction**
(§8) replaces the flat middle segment with a power law obtained by directly integrating Numan et al
2014's **own measured regression elasticities** `e = dlnF/dlnP` (their %ΔCBF/%ΔMAP regression
coefficients) outward from a reference pressure: `F(P) = F0*(P/P_ref)^e` — the textbook definition of
"integrate a measured elasticity," not a heuristic curve shape. **CO2 reactivity (§10) uses the exact
same clipped-line family a third time**: linear above a threshold near the eucapnic reference (matching
Battisti-Charbonney et al 2011's [7] own reported linear supra-threshold regime), clipped to a floor
below it (matching their own reported sigmoidal/saturating sub-threshold shape) — never an unbounded
line, because an unbounded line goes **negative** (unphysical) at realistic hypocapnic PaCO2 (§10's
forced adversary). Three different physiological curves in this document, built from the **same**
underlying geometric primitive (a regulation line clipped to a physical floor/ceiling) — not three
unrelated heuristics.

## 2. Method, in one paragraph

Two upstream twin layers are read, no re-solve: `arterial_pressure.py`'s own MAP (classic-formula
93.33mmHg; Windkessel-true-mean 91.89mmHg) and `acid_base_co2.py`'s own acute-respiratory scenario
(PaCO2 40→50mmHg, its `delta_paco2=10` entry, reused verbatim, not invented). The resting CBF baseline
(whole/grey/white) is tested against two independent modern measurement methods — PET [9] and ASL-MRI
[10] — that the task's own classic 50/80/20 figures were never fit to. The autoregulation curve is built
**two ways**: the idealized classic (Lassen [1]) flat plateau, and a modern "leaky-plateau" curve that
literally integrates Numan et al 2014's [4] own measured elasticities — the falsifiable, quantified
comparison the task asks be held open. Both curve variants are tested against a **forced adversary**
("no autoregulation," a fixed resistance) and against the task's own explicit falsifier (a MAP=50
perturbation, below the classic knee, must drop CBF). CO2 reactivity is calibrated to Battisti-
Charbonney et al 2011's [7] live-quoted TCD slope (5.5%/mmHg) and cross-checked against an independent
BOLD-context source (Mutch et al 2012 [8], 5-11%/mmHg) and the task's own stated 3-4%/mmHg band — and
tested against a forced unbounded-linear adversary that goes physically negative at extremes. `CPP=
MAP-ICP` (StatPearls [12][13]) is tested by feeding this twin's OWN already-computed MAP through a
normal-ICP sweep (cross-layer over-determination with the arterial-pressure layer) and by a forced
"ICP-blind" adversary at a pathologically raised ICP. The acid-base coupling is tested by feeding
`acid_base_co2.py`'s OWN PaCO2 perturbation through this model.

## 3. Citations — every PMID/DOI verified LIVE this session (NCBI eutils + EuropePMC), not recalled

**Two numbers are flagged at a lower confidence tier** (marked ⚠WebFetch below): they were obtained via
a WebFetch summarization pass over PMC full text, not independently raw-machine-verified, because a
direct `curl` fetch of PMC full-text pages returned a Google reCAPTCHA challenge page this session (a
real, disclosed tool limitation — confirmed by inspecting the raw returned HTML, not assumed). Every
other number below was verified via NCBI's own `efetch` abstract text or EuropePMC's own machine-readable
`abstractText` field — genuinely live, not recalled.

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Lassen NA (1959). "Cerebral blood flow and oxygen consumption in man." *Physiol Rev* 39(2):183-238. | **13645234**, DOI `10.1152/physrev.1959.39.2.183` (verified live) | THE classic paper compiling pooled data into the textbook "flat 60-150mmHg" curve. Pre-abstracting era (no PubMed abstract) — **provenance only**; the 60-150mmHg numbers are used here as the task's own stated textbook target, not re-extracted from this paper's own tables. |
| 2 | Kety SS, Schmidt CF (1948). "The nitrous oxide method for the quantitative determination of cerebral blood flow in man: theory, procedure and normal values." *J Clin Invest* 27(4):476-83. | **16695568**, DOI `10.1172/JCI101994`, PMC439518 (verified live) | THE original quantitative resting-CBF method paper. Pre-abstracting era — provenance only. |
| 3 | Kety SS, Schmidt CF (1948). "The effects of altered arterial tensions of carbon dioxide and oxygen on cerebral blood flow and cerebral oxygen consumption of normal young men." *J Clin Invest* 27(4):484-92. | **16695569**, DOI `10.1172/JCI101995`, PMC439519 (verified live) | THE original CO2-reactivity quantification paper. Pre-abstracting era — provenance only. |
| 4 | Numan T, Bain AR, Hoiland RL, Smirl JD, Lewis NC, Ainslie PN (2014). "Static autoregulation in humans: a review and reanalysis." *Med Eng Phys* 36(11):1487-95. | **25205587**, DOI `10.1016/j.medengphy.2014.08.001` (verified live, full abstract fetched) | **THE key modern reanalysis.** Quoted verbatim: "based on 40 studies (49 individual experimental protocols)... linear regression coefficient between MAP and CBF... of 0.82±0.77%ΔCBF/%ΔMAP during **decreases** in MAP (n=23)... 0.21±0.47%ΔCBF/%ΔMAP during **increases** (n=26; p<0.001). After correction for... PaCO2, the slopes were not significantly different: 0.64±1.16 (n=16) and 0.39±0.30 (n=12)... (p=0.60)." A real, quantified, nonzero "leak" — much of the apparent hysteresis is itself a PaCO2 confound, disclosed not hidden. |
| 5 | Willie CK, Tzeng YC, Fisher JA, Ainslie PN (2014). "Integrative regulation of human brain blood flow." *J Physiol* 592(5):841-59. | **24396059**, PMC3948549, DOI `10.1113/jphysiol.2013.268953` (verified live via **EuropePMC's own machine-readable `abstractText` field**, not a summarizer) | Abstract states, verbatim, as **one of four key theses** of the whole review: cerebral autoregulation "does **not** maintain constant perfusion through a mean arterial pressure range of 60-150mmHg" — named among "evidence against several canonized paradigms of CBF control." ⚠WebFetch (PMC full text, unverified raw): CO2 reactivity "approximate 3-6% increase and/or 1-3% decrease in flow per mmHg... above and below eupnoeic PaCO2"; "CO2 reactivity of... grey matter is greater than... white matter." |
| 6 | Willie CK, Colino FL, Bailey DM, et al (2011). "Utility of transcranial Doppler ultrasound for the integrative assessment of cerebrovascular function." *J Neurosci Methods* 196(2):221-37. | **21276818**, DOI `10.1016/j.jneumeth.2011.01.011` (verified live) | TCD methodology review — topical support for TCD velocity as a CBF proxy throughout this script. |
| 7 | Battisti-Charbonney A, Fisher J, Duffin J (2011). "The cerebrovascular response to carbon dioxide in humans." *J Physiol* 589(Pt 12):3039-48. | **21521758**, PMC3139085, DOI `10.1113/jphysiol.2011.206052` (verified live, full abstract fetched) | **PRIMARY, quantitative, TCD.** Quoted verbatim: "the MCAv response to CO2 was **sigmoidal** below a discernible threshold CO2 tension... centred at a CO2 tension close to normal resting values (overall mean 36mmHg)... **Above** this threshold both MCAv and MAP increased **linearly** with CO2 tension... overall mean slopes **5.5% mmHg⁻¹** and 2.1mmHg mmHg⁻¹, respectively." The primary anchor for this script's hypercapnic slope and two-regime shape. |
| 8 | Mutch WA, Mandell DM, Fisher JA, Mikulis DJ, Crawley AP, Pucci O, Duffin J (2012). "Approaches to brain stress testing: BOLD magnetic resonance imaging with computer-controlled delivery of carbon dioxide." *PLoS One* 7(11):e47443. | **23139743**, PMC3489910, DOI `10.1371/journal.pone.0047443` (verified live, full abstract fetched) | Independent (BOLD-context) secondary source, quoted verbatim: CO2 "changing cerebral blood flow (CBF) by up to **5-11 percent/mmHg** in normal adults" — corroborates that real CO2 reactivity runs **higher** than the task's own stated ~3-4%/mmHg. |
| 9 | Ito H, Kanno I, Kato C, Sasaki T, et al (2004). "Database of normal human cerebral blood flow... a multicentre study in Japan." *Eur J Nucl Med Mol Imaging* 31(5):635-43. | **14730405**, DOI `10.1007/s00259-003-1430-8` (verified live, full abstract fetched) | **PRIMARY, quantitative, PET.** Quoted verbatim: "Subjects comprised 70 healthy volunteers... Overall mean±SD values for cerebral **cortical** regions were: CBF=**44.4±6.5** ml·100ml⁻¹·min⁻¹." A modern cross-check against the task's stated "grey~80." |
| 10 | Alisch JSR, Khattar N, Kim RW, et al (2021). "Sex and age-related differences in cerebral blood flow... pseudo-continuous arterial spin labeling MRI." *Aging (Albany NY)* 13(4):4911-4925. | **33596183**, PMC7950235, DOI `10.18632/aging.202673` (verified live, full abstract fetched — N=80, "reference CBF values for the... ISMRM Perfusion Study Group") | ⚠WebFetch (PMC Table 3, unverified raw): overall cohort mean GM CBF=**33.0±6.17**, WM CBF=**18.0±4.13** ml/100g/min. The abstract's own qualitative finding (GM>WM, GM falls / WM rises with age) IS independently, primary-verified at full confidence; only the exact decimals carry the flag. |
| 11 | Czosnyka M, Pickard JD (2004). "Monitoring and interpretation of intracranial pressure." *J Neurol Neurosurg Psychiatry* 75(6):813-21. | **15145991**, PMC1739058, DOI `10.1136/jnnp.2003.033126` (verified live) | Topical/provenance anchor: "Information which can be derived from ICP and its waveforms includes cerebral perfusion pressure (CPP), regulation of cerebral blood flow and volume..." — CPP-guided therapy as a real, measured clinical quantity. |
| 12 | Mount CA, Das JM. "Cerebral Perfusion Pressure." *StatPearls* [Internet]. | **30725956** (verified live, full text fetched) | Quoted verbatim: CPP "is the difference between the mean arterial pressure (MAP) and the intracranial pressure (ICP)... Normal CPP lies between **60 and 80mmHg**." |
| 13 | Pinto VL, Adeyinka A. "Increased Intracranial Pressure." *StatPearls* [Internet]. | **29489250** (verified live, full text fetched) | Quoted verbatim: "Normal intracranial pressure (ICP) in adults typically ranges from **7 to 15mmHg** in the supine position. Values above **20 to 25mmHg** are generally considered pathological." States the Monro-Kellie doctrine and that raised ICP "can reduce cerebral perfusion pressure (CPP)... ultimately cause ischemia." |
| 14 | Silverman A, Petersen NH. "Physiology, Cerebral Autoregulation." *StatPearls* [Internet]. | **31985976** (verified live, full text fetched) | Quoted verbatim: distinguishes autoregulation from "carbon dioxide reactivity" and "flow-metabolism coupling" as **three separate** mechanisms; resistance follows "the Hagen-Poiseuille equation"; flags "individualizing cerebral perfusion pressure targets" as an active research frontier — supports holding LLA/UBP individual variability open. |
| 15 | Rangel-Castilla L, Gasco J, Nauta HJ, Okonkwo DO, Robertson CS (2008). "Cerebral pressure autoregulation in traumatic brain injury." *Neurosurg Focus* 25(4):E7. | **18828705**, DOI `10.3171/FOC.2008.25.10.E7` (verified live, full abstract fetched) | Quoted verbatim: "cerebral autoregulation is more a **concept** than a physically measurable entity... assessment of cerebral autoregulation is best achieved with **dynamic** autoregulation methods" — the static-vs-dynamic distinction held open (§0, §13). |

## 4. Inputs — the twin's own MAP and PaCO2 perturbation, two decorrelated upstream legs, no re-solve

| upstream layer | source | value used |
|---|---|---:|
| Arterial pressure (classic formula) | `arterial_pressure_results.json` `map_approximation.classic_map_mmhg` | 93.33 mmHg |
| Arterial pressure (Windkessel true mean) | `arterial_pressure_results.json` `pulse_pressure.windkessel_sim.mean_p_analytical` | 91.89 mmHg |
| Acid-base (acute respiratory scenario) | `acid_base_co2_results.json` `step3_acute_respiratory[1]` | PaCO2 40→50 mmHg |

Both reused verbatim, not re-derived or invented — the `delta_paco2=10` entry is `acid_base_co2.py`'s
own first acute-perturbation scenario, asserted equal at runtime (`assert abs(delta_paco2_from_sibling
- PACO2_COUPLED_DELTA) < 1e-9`).

## 5. Resting CBF baseline — task's stated bands vs two independent modern methods (PET, ASL)

| quantity | task's stated target | pre-registered band (±15%) | modern measurement | verdict |
|---|---:|---|---|---|
| Whole-brain CBF | 50 mL/100g/min | [42.5, 57.5] | Ito 2004 PET [9] cortical=44.4±6.5 lands **inside** | interesting coincidence, not a grey-matter match |
| Grey matter CBF | 80 mL/100g/min | [68.0, 92.0] | Ito 2004 PET [9]=44.4±6.5 | **MISS** |
| Grey matter CBF | 80 mL/100g/min | [68.0, 92.0] | Alisch 2021 ASL [10, ⚠WebFetch]=33.0±6.17 | **MISS** |
| White matter CBF | 20 mL/100g/min | [17.0, 23.0] | Alisch 2021 ASL [10, ⚠WebFetch]=18.0±4.13 | **PASS** |

**Grey-mass-fraction self-consistency** (does the task's own triad of numbers reconcile?): solving
`f·80 + (1-f)·20 = 50` gives **f_grey = 0.500** — comfortably inside a wide plausibility sanity band
[0.25, 0.75]: **PASS**. **Cross-METHOD reconciliation** (does any grey-mass-fraction blend Alisch's OWN
GM/WM ASL numbers into Ito's OWN PET cortical number?): solving `f·33.0 + (1-f)·18.0 = 44.4` gives
**f = 1.76** — **outside** [0,1], i.e. **no valid blend exists**: a genuine, disclosed inter-method
(PET-vs-ASL) absolute-quantification gap, machine-computed, not hidden. Plausible (not adjudicated)
explanations: partial-volume/ROI-averaging dilution in modern cortical ROIs vs older, more localized
measurement traditions; a documented tendency for ASL to under-quantify absolute CBF vs PET/Kety-Schmidt.

**Honest reading**: the task's white-matter figure (~20) replicates cleanly against modern ASL. The
grey-matter figure (~80) does **not** replicate against either modern method tested (PET: 44.4; ASL:
33.0) — both land 45-60% below the task's stated value. This is reported exactly as measured, not
smoothed into a false consistency (§13).

## 6. Classic (idealized Lassen) autoregulation curve — flat by construction, non-falsifying

Built via §1's clipped-resistance geometry (`LLA=60`, `UBP=150`, `F0=50`): **max deviation from
perfectly flat within [60,150] = 0.000000%** — PASS, but explicitly flagged **non-falsifying** (true by
construction, the same discipline `acid_base_co2.py` applies to its own reference-point sanity check).
Continuous at both knees (PASS/PASS); strictly monotonic falling below LLA and rising above UBP
(PASS/PASS). This is the textbook shape the task states as the baseline model — reproduced exactly, then
immediately tested against modern data below.

## 7. The decisive falsifier — modern "leaky-plateau" curve (Numan 2014's own elasticities integrated)

**Pre-registered prediction**, stated before computing: the modern/Numan-anchored curve **will** deviate
from perfectly flat by more than 5% somewhere strictly inside the nominal 60-150mmHg plateau.

**Measured: max deviation-from-flat inside [60,150] = 30.4%** — confirmed. Worked examples (closer to
the reference point, more representative than the extreme knee-edge value): at CPP=65mmHg (5mmHg inside
the lower knee), classic predicts F=50.0 (flat); the Numan-anchored model predicts **F=37.17 (−25.7%)**.
At CPP=140mmHg (10mmHg inside the upper knee): **F=54.44 (+8.9%)**. Continuous at both knees (power-law
meets passive line, PASS/PASS).

**Symmetric-QC caveat on the 30.4% figure itself**: this is a **power-law extrapolation** of Numan et
al's measured LOCAL regression slope over the full 33mmHg half-width of the nominal plateau (93.33→60).
The primary paper's own exclusion criterion was only "ΔMAP<5%" (a *lower* bound on the ΔMAP studies had
to clear to be included) — it does not establish that the elasticity stays *constant* all the way out to
a 36%-of-MAP swing. The 30.4% figure at the exact knee-edge is therefore likely an **upper bound** on
what a constant-elasticity model implies, not a precision measurement; the more moderate ±26%/+9% values
at CPP=65/140 (5-10mmHg inside the knees) are the better-supported headline numbers. Either way, the
qualitative conclusion is unambiguous and externally corroborated twice over (Numan's own regression AND
Willie's own review thesis): **the classic flat plateau is a useful order-of-magnitude idealization, not
a precise description of the measured relationship.**

## 8. Forced adversary (void-floor): "no autoregulation" vs the regulated model

Constant resistance fixed at the resting operating point (`R=P_ref/F0`), extended to every pressure —
the sharpest, cleanest test that autoregulation is doing real, necessary work, not decorative. Total CBF
variation across [60,150]mmHg: **WITH** autoregulation = 0.0000 mL/100g/min (0, by construction) vs
**WITHOUT** (void-floor) = **48.22** mL/100g/min: **PASS** (void-floor swings far more). At CPP=60, the
void-floor predicts F=32.14 vs regulated F=50.0; at CPP=150, void-floor predicts F=80.36 vs regulated
F=50.0 — **a 150% swing with no mechanism, vs 0% with one.**

## 9. The task's explicit falsifier — MAP=50 (below the classic knee) must drop CBF

**Pre-registered gate**: CBF must drop by ≥10% under BOTH curve variants, under BOTH plausible
interpretations of "MAP=50" (literal, treated as if it were CPP directly; and MAP=50 run through
`CPP=MAP-ICP` at a typical mid-normal ICP=11mmHg → CPP=39mmHg).

| interpretation | classic-model drop | leaky-model drop |
|---|---:|---:|
| Literal MAP=50 as CPP | **16.7%** | **42.0%** |
| MAP=50 → CPP=39 (via ICP=11) | **35.0%** | **54.8%** |

**All four ≥10%: PASS**, and robust to curve-variant choice — the qualitative conclusion (CBF drops)
holds whether the idealized or the modern-leaky curve is used, and whether ICP is subtracted first or
not. This is the single most decisive, unambiguous result in this document.

## 10. CO2 reactivity — calibrated to TCD (Battisti-Charbonney), cross-checked vs BOLD (Mutch) and the task's own band

Reference point `g(PaCO2=40)=1.000000` exactly (non-falsifying sanity check, PASS). The model's own
hypercapnic slope (a calibration check, **not** independent evidence, since it is the anchor itself) =
**5.50%/mmHg**, matching Battisti-Charbonney [7] exactly by construction. The genuine, independent
cross-checks:

- vs Mutch et al [8]'s stated 5-11%/mmHg range: **PASS** (5.5 falls inside) — real external corroboration.
- vs the **task's own** pre-registered 3-4%/mmHg band: **MISS** — a disclosed tension, not force-fit.
  Both live-verified sources measuring CO2 reactivity (TCD [7], BOLD-context [8]) independently find the
  real value running **higher** than the task's stated textbook figure.

**Forced adversary — the naive symmetric-unbounded-linear model** (the same hypercapnic slope, applied
in both directions, no floor, no sigmoid — the "simplest possible" model without the correction this
document actually uses): over PaCO2∈[15,70]mmHg, the naive adversary's minimum multiplier = **−0.375**
(**negative — unphysical**, PASS: the adversary correctly falls). The saturating model's minimum =
**0.500** (stays physically bounded, PASS). **Void-floor (zero CO2 reactivity)**: 0.0% swing (by
definition) vs the real model's **215.0%** swing across the same sweep — PASS by a wide margin.

## 11. CPP = MAP − ICP — cross-layer over-determination with `arterial_pressure.py`'s own MAP

Using StatPearls [12]'s formula directly on this twin's own already-computed MAP, swept over the normal
ICP range [7,15]mmHg [13]:

| ICP source | CPP from classic MAP (93.33) | CPP from Windkessel-true-mean (91.89) |
|---|---:|---:|
| ICP=7 (low-normal) | 86.3 | 84.9 |
| ICP=11 (mid-normal) | 82.3 | 80.9 |
| ICP=15 (high-normal) | 78.3 | 76.9 |

**All six values land inside the autoregulation plateau [60,150]: PASS** — a genuine cross-layer
over-determination spanning **three** already-built twin layers (`cardiac_output.py` →
`arterial_pressure.py` → this cerebral-autoregulation model): they agree the resting operating point is
safely inside the regulated range, as a healthy resting subject's must be. One disclosed boundary note:
at mid-normal ICP, CPP=82.3mmHg sits just above StatPearls' own stated normal-CPP upper bound of 80mmHg
— a small, real mismatch, not forced to fit (§13).

**Forced adversary — "ICP-blind" model** (uses MAP directly as if it were CPP, ignoring ICP entirely),
tested at a pathologically raised ICP (StatPearls [13]'s own 25mmHg threshold), with a full OODA trail:

- **Observe** (first pass, reported not discarded): at the twin's own high-headroom resting MAP
  (93.33mmHg) + ICP=25: true CPP=68.33, *still* inside the plateau (>60) — both models trivially agree
  (F=50, 0% error). A real, informative finding: a high-MAP patient has enough "plateau headroom" that
  ICP-blindness alone does not yet matter.
- **Orient**: the adversary needs a MAP that *looks* unremarkable on its own (still ≥ LLA, so a clinician
  reading MAP alone would not be alarmed) but where the TRUE CPP, once ICP is correctly subtracted,
  falls below LLA — isolating ICP-blindness as the error source, not confounding it with an
  already-abnormal MAP.
- **Decide/Act**: re-test at a borderline-but-still-≥LLA MAP=70mmHg + the same pathological ICP=25:
  true CPP=45.00mmHg (now *below* LLA=60) → correct model predicts F=37.50; the ICP-blind model (reads
  MAP=70 directly, sees it as ≥LLA, predicts flat/normal) predicts F=50.00 — **error=33.3%: PASS**
  (ignoring ICP produces a clinically meaningful, *falsely reassuring* CBF prediction).
- **Convergence**: error at high-headroom MAP=93.3 = 0.0% vs error at borderline MAP=70.0 = 33.3% —
  ICP-blindness is not a fixed error; its danger **depends on how much plateau headroom the baseline MAP
  has**, a genuine geometric consequence of the clipped-resistance structure (§1), not an arbitrary
  scenario choice.

## 12. Cross-layer acid-base coupling — `acid_base_co2.py`'s own PaCO2 perturbation propagated through

`acid_base_co2.py`'s own acute-respiratory scenario (PaCO2 40→50mmHg, its own `delta_paco2=10`, reused
verbatim) propagated through this model predicts a CBF multiplier **g=1.5500 (+55.0%)**: **PASS**
(positive, bounded, order-of-magnitude-plausible — a 10mmHg PaCO2 rise at 5.5%/mmHg is exactly 55% by
construction; the value of this check is that it is a real, load-bearing data read from a sibling twin
layer, not a fresh invented number).

## 13. Combined joint model — F(CPP, PaCO2), sanity sweep

`F(P,C) = F_autoregulation(P) · g(C)` — both mechanisms independently adjust a shared resistance state
(a disclosed structural simplification, see Honest Gaps).

| scenario | CPP (mmHg) | PaCO2 (mmHg) | predicted CBF | vs rest |
|---|---:|---:|---:|---:|
| rest | 93.3 | 40.0 | 50.00 | +0.0% |
| hypotensive + hypercapnic stress | 50.0 | 55.0 | 52.93 | +5.9% |
| hypertensive + hypocapnic stress | 160.0 | 28.0 | 44.78 | −10.4% |
| severe hypotension, below LLA | 40.0 | 40.0 | 23.20 | −53.6% |

All four scenarios stay positive/bounded: **PASS**. The rest scenario reproduces `F0` exactly by
construction. Note the "hypotensive+hypercapnic" scenario shows CO2-driven vasodilation *partially
offsetting* a pressure-driven drop (net +5.9%, despite CPP=50 being below LLA on its own) — a real,
qualitatively-expected interaction of the two mechanisms even under this document's disclosed
multiplicative-independence simplification.

## 14. Falsifier verdict — stated exactly as pre-registered

> Does the model reproduce the measured autoregulation curve (CBF ~flat 60-150mmHg, falling steeply
> outside — Lassen's classic curve, re-examined by modern PET/TCD which show a NARROWER true plateau
> than the textbook) AND the measured CO2-reactivity slope (~3-4%/mmHg by TCD/BOLD)? A perturbation at
> MAP=50 (below the knee) must drop CBF.

**Autoregulation curve: PASS, with the modern-vs-classic tension explicitly quantified, not smoothed
over.** The classic flat-plateau shape is reproduced exactly (§6, non-falsifying by construction). The
task's own instruction to hold open "a narrower true plateau than the textbook" is directly confirmed,
twice, by independent modern sources: Willie et al 2014 [5] state this as one of their review's own four
central theses; Numan et al 2014 [4] quantify it as a real, nonzero regression slope inside the nominal
plateau, reproduced here as a 26-30% deviation-from-flat near the knees (§7, with the appropriate
extrapolation-range caveat). **The MAP=50 falsifier, the task's single most explicit requirement, PASSES
decisively and robustly**: CBF drops 16.7-54.8% depending on curve variant and MAP/CPP interpretation,
all four combinations clearing a pre-registered ≥10% threshold by a wide margin (§9).

**CO2 reactivity: PASS on external cross-check (Mutch), MISS on the task's own literal 3-4%/mmHg band —
disclosed, not laundered.** The primary TCD literature (Battisti-Charbonney [7], 5.5%/mmHg) and an
independent BOLD-context source (Mutch [8], 5-11%/mmHg) **agree with each other** and both sit **above**
the task's stated band. This is reported honestly as a real tension between the task's textbook framing
and live-verified modern TCD/BOLD measurements, not forced into a false match (§10).

**Symmetric QC, stated plainly**: nothing here is proven. The classic Lassen plateau is real,
useful, and order-of-magnitude correct, but modern reanalysis shows it is neither as flat nor as
universal as commonly taught — held open, quantified two independent ways (§7, §8). The task's own
grey-matter (~80) and CO2-reactivity (~3-4%/mmHg) figures do **not** replicate cleanly against live
modern measurements, while white-matter CBF (~20) and the MAP=50 falsifier both replicate cleanly. Static
and dynamic autoregulation are explicitly **not** the same thing, and only static autoregulation is
modeled here (§0, §13). What earns more than "lands in a wide band" here is the **forced adversary set**
(no-autoregulation void-floor, unbounded-CO2-adversary, ICP-blind adversary — all three fall cleanly,
§8/§10/§11) and the **cross-layer over-determination** (three already-built twin layers agreeing the
resting operating point sits safely inside the regulated plateau, §11).

## 15. Pre-registered gates — 19/19 PASS + 7 disclosed open items

```
grey_white_fraction_selfconsistent_plausible:              PASS (f_grey=0.500 in [0.25,0.75])
classic_curve_flat_in_plateau_nonfalsifying:                PASS (0.000000% deviation, by construction)
classic_curve_continuous_at_knees:                          PASS
classic_curve_falls_outside_plateau:                        PASS
leaky_plateau_confirms_modern_vs_classic_tension_gt5pct:    PASS (30.4% max deviation, pre-registered >5%)
leaky_plateau_continuous_at_knees:                          PASS
void_floor_no_autoregulation_adversary_falls:               PASS (48.22 vs 0.00 mL/100g/min variation)
map50_perturbation_drops_cbf_classic:                       PASS (16.7-35.0% drop, gate >=10%)
map50_perturbation_drops_cbf_leaky:                         PASS (42.0-54.8% drop, gate >=10%)
co2_reference_point_nonfalsifying:                          PASS (g(40)=1.000000)
co2_slope_matches_calibration_source:                       PASS (5.50 vs 5.5%/mmHg, by construction)
co2_slope_within_independent_mutch_range:                   PASS (5.5 in [5,11])
co2_naive_unbounded_adversary_falls_negative:               PASS (min=-0.375)
co2_saturating_model_stays_physical:                        PASS (min=0.500)
co2_void_floor_zero_reactivity_beaten:                      PASS (215.0% vs 0.0% swing)
cpp_crosslayer_overdetermination_in_plateau:                PASS (all 6 CPP values in [60,150])
icp_blind_adversary_falls_at_pathological_icp:              PASS (33.3% error at borderline MAP)
acid_base_crosslayer_coupling_positive_bounded:             PASS (+55.0%)
combined_model_bounded:                                     PASS

OPEN MODELING UNCERTAINTY (disclosed, does NOT gate overall_pass):
classic_plateau_is_textbook_idealization_not_measured_reality:  Willie[5]+Numan[4], held OPEN
co2_reactivity_measured_higher_than_tasks_stated_band:          5.5-11 vs stated 3-4%/mmHg
grey_matter_cbf_measured_lower_than_tasks_stated_80:            PET 44.4 / ASL 33.0 vs stated 80
static_vs_dynamic_autoregulation_not_modeled_as_separate:       only static modeled, per task's own framing
hypocapnic_slope_and_alisch_values_are_webfetch_summary_tier:   PMC direct-fetch reCAPTCHA-blocked
multiplicative_pressure_co2_independence_disclosed_simplification: no interaction term modeled
cpp_normal_band_boundary_case:                                  82.3 vs stated upper bound 80
```

**Overall: PASS** (19/19 pipeline-correctness/forced-adversary gates; 2 independent runs produce
byte-identical JSON, machine-verified, not eyeballed).

## 16. Honest gaps (disclosed, not hidden)

- **Static autoregulation only.** No dynamic/transient-response model (Tiecks ARI, thigh-cuff-release
  transient) exists in this script. Rangel-Castilla et al 2008 [15]: "cerebral autoregulation is more a
  concept than a physically measurable entity... best achieved with dynamic autoregulation methods."
  Held explicitly OPEN, per the task's own instruction.
- **The classic 60-150mmHg plateau is a textbook idealization**, not a precise description of the
  measured relationship — confirmed two independent ways (§7, §8), held open rather than resolved.
- **Grey-matter CBF (~80) does not replicate** against either PET [9] (44.4) or ASL [10] (33.0); the two
  modern methods do not even fully reconcile with each other (§5's cross-method check, required
  fraction f=1.76, outside [0,1]).
- **CO2 reactivity (~3-4%/mmHg) does not replicate** against live TCD [7] or BOLD-context [8] literature
  (5.5-11%/mmHg) — both independently run higher than the task's stated band.
- **Two numbers are WebFetch-summary-tier, not independently raw-verified**: the 1-3%/mmHg hypocapnic
  slope [5] and the exact GM=33.0/WM=18.0 decimal values [10] — direct `curl` access to PMC full text
  was reCAPTCHA-blocked this session (confirmed from the raw returned HTML, a genuine tool limitation).
  Both papers' own abstracts ARE independently, live, primary-verified at full confidence; only these
  specific decimals carry the flag.
- **Multiplicative pressure×CO2 independence is a disclosed simplification.** Real cerebrovascular
  physiology has documented pressure-CO2 interaction effects (e.g. hypercapnic vasodilation consuming
  some of the vessel's dilatory reserve, narrowing the effective autoregulation range) that this
  model's structural form does not capture — no specific interaction citation was live-verified this
  session for a quantified magnitude, so this is an open structural gap, not an uncited assertion.
- **`LLA=60`/`UBP=150` are the classic population-level constants**, not subject-specific — StatPearls
  [14] itself flags "individualizing cerebral perfusion pressure targets" as an active research
  frontier; individual LLA is known to vary substantially, especially in TBI/critical illness [15].
- **No subject-specific CBF, ICP, or TCD measurement exists for subject2** — every external anchor here
  is population-level literature, the same first-step scope every other `MECHANISM_*` layer discloses.
- **Confidence tier, stated precisely**: the `CPP=MAP-ICP` arithmetic and the MAP=50 falsifier (§9) are
  **in-vivo-anchored** (StatPearls clinical reference ranges + this twin's own already-certified
  upstream MAP). The autoregulation-curve shape (§6-§8) is **in-vivo-anchored with a disclosed,
  quantified, modern-vs-classic tension** (Numan's real patient-study reanalysis). The CO2-reactivity
  leg (§10) is **in-vivo-anchored via TCD** (Battisti-Charbonney, primary) **with a disclosed band
  mismatch** against the task's own stated figure. The grey/white regional baseline (§5) is the
  **weakest** leg — population PET/ASL literature that does not itself agree at the ~80/33/44 mL/100g/min
  level, an honest, disclosed, unresolved tension.

## 17. Couples to (graph context, not mutated)

- **Couples to arterial-pressure** (`scripts/msk/arterial_pressure.py`,
  `docs/MECHANISM_ARTERIAL_PRESSURE.md`): `CPP=MAP-ICP` directly consumes this twin's own already-computed
  MAP (both the classic-formula and Windkessel-true-mean estimates, §4, §11) — a real, load-bearing data
  read, not a prose reference. The cross-layer over-determination (§11) is load-bearing for this
  document's own confidence tier.
- **Couples to acid-base** (`scripts/msk/acid_base_co2.py`, `docs/MECHANISM_ACID_BASE_CO2.md`): the CO2-
  reactivity leg (§10, §12) directly consumes that layer's own `normal_paco2_mmhg=40` reference point and
  its own acute-respiratory PaCO2 perturbation scenario (40→50mmHg) — again a real data read, not prose.
- **`data/MECHANISM_ANCHOR_GRAPH.json` node referenced, not edited** (out of scope for this cert — the
  graph is a large, separately-maintained structure; this doc reports against it in prose only):
  `NEU-NEUROVASCULAR-COUPLING` (status `OPEN`, SEED-DESIGN tier, claim: "NVC-reserve-gate... require an
  independent CBF or CVR channel... before treating a BOLD-fMRI amplitude as a neural-activity readout")
  — that node's own `couples_to` list explicitly names "respiratory/CO2-chemoreflex physiology: PaCO2 is
  the dominant CBF driver, shared substrate with breath-hold BOLD calibration" — this document supplies a
  first quantitative mechanism for exactly that shared substrate (the CO2-reactivity curve, §10), but
  does **not** resolve `NEU-NEUROVASCULAR-COUPLING` itself (a different claim, about BOLD-fMRI
  interpretation validity, not brain perfusion per se) — remains **OPEN**, untouched by this document.
- **Couples to** (downstream, not executed here): any future oxygen-delivery/metabolism layer would need
  this document's CBF output as its flow term, the same pattern `blood_oxygen_transport.py` establishes
  for the arterial-pressure layer's MAP output.

## 18. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/cerebral_autoregulation.py
```
Requires `data/arterial_pressure/arterial_pressure_results.json` and
`data/msk_smoketest/acid_base_co2/acid_base_co2_results.json` to already exist (both read-only, never
modified). No OpenSim call, no new subject-trial data, runs in under a second, deterministic (verified:
2 independent runs produce byte-identical JSON).

**Paths**: script `scripts/msk/cerebral_autoregulation.py`; evidence
`data/cerebral_autoregulation/cerebral_autoregulation_results.json`; this doc
`docs/MECHANISM_CEREBRAL_AUTOREGULATION.md`; upstream (read-only) inputs
`data/arterial_pressure/arterial_pressure_results.json` and
`data/msk_smoketest/acid_base_co2/acid_base_co2_results.json`.
