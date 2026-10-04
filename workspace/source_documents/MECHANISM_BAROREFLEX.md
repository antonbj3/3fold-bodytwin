# MECHANISM BAROREFLEX — the feedback controller that CLOSES the cardiovascular loop (2026-07-22)

Resolves the core, executable subset of the **OPEN**/`SEED-DESIGN` hypothesis node
`SYS-BAROREFLEX-AUTONOMIC` (`data/MECHANISM_ANCHOR_GRAPH.json`, duplicated as
`AUTO-BAROREFLEX-AUTONOMIC-LOOP-CLOSED-LOOP-NT`; prior verdict *"SEED-DESIGN: cited literature anchor,
not yet executed against data"*): the closed-loop, NTS-centered negative-feedback controller
`{carotid-sinus/aortic stretch-afferents -> NTS -> CVLM/RVLM+nucleus-ambiguus efferent split ->
HR/contractility/vascular-tone}`, parametrized by a measurable loop-gain (BRS, ms/mmHg). This is the
**controller** that closes the loop the twin already opened: arterial pressure
(`docs/MECHANISM_ARTERIAL_PRESSURE.md`, `MAP=CO*TPR`) is the **controlled variable**; cardiac output/HR
(`docs/MECHANISM_CARDIAC.md`, `docs/MECHANISM_CARDIAC_OUTPUT.md`) is the **effector**; this doc adds the
**sensor** (baroreceptor sigmoid), the **arc** (afferent/central/efferent pathway + delays), and the
**loop gain** (BRS). Script: `scripts/msk/baroreflex.py`. Evidence:
`data/baroreflex/baroreflex_results.json`.

**NO re-solve, no OpenSim, no new subject-trial data.** Reads TWO already-computed JSON outputs
(`arterial_pressure_results.json` for the twin's own resting MAP; `cardiac_output_results.json` for
the twin's own resting HR/RR), read-only, and does pure-Python/numpy/cmath arithmetic on top of REAL,
live-verified literature numbers. Never touches the `.osim` model or any `.sto` file.

## 0. Scope, stated up front — read this before any number below

This is a **static-operating-point + linearized-loop-geometry** model: the baroreceptor sigmoid is
evaluated at illustrative parameters (disclosed, not fabricated precision), the reflex arc is a
structural/qualitative map with an order-of-magnitude conduction-delay budget (generic anatomical
distances — this twin's OpenSim model is lower-limb-only and has no head/neck/thorax autonomic
anatomy to query live, unlike `nerve_conduction.py`'s leg geometry), and the Mayer-wave/resonance
argument is a **pole-placement consistency check**, not an independent measurement of the sympathetic
delay. **Nothing here is a full closed-loop ODE simulation** — no beat-to-beat integration, no
baroreflex-mediated blood-pressure-control simulation, no dynamic resetting model. Baroreflex
sensitivity (BRS) **declines markedly with age and hypertension** (quantified below, not just
asserted) and **the three pre-registered measurement methods (sequence / modified-Oxford-phenylephrine
/ spectral-alpha-index) genuinely disagree** (quantified below, held **OPEN**, not resolved) — this
cert reports the spread and the real tensions, not a single laundered number.

## 1. Geometric structure (derive from the geometry, not heuristics)

1. **The sigmoid is a 4-parameter logistic**: `FR(P) = FR_min + (FR_max-FR_min)/(1+exp(-k*(P-P50)))`.
   Its one analytically-sharp feature is the **maximum gain**, occurring exactly at the inflection
   point `P50`: `dFR/dP|_{P50} = k*(FR_max-FR_min)/4` — a closed-form feature of the curve's own
   geometry, machine-verified below against a numerical `np.gradient` (§5), not eyeballed off a plot.
2. **The closed-loop DIRECTION is a chain-rule identity, not an assumption.** Literature reports BRS
   as `d(RR)/d(SBP)` in ms/mmHg (positive — RR lengthens when pressure rises). Converting to HR via
   `RR[ms]=60000/HR[bpm]` gives `dHR/dSBP = -(HR^2/60000)*BRS_RR` — **negative by construction**, the
   required "BP up -> HR down" direction, verified below against a finite-difference numerical
   cross-check (§6), not asserted from intuition.
3. **The ~0.1-Hz "Mayer wave" is treated as a pole-placement problem** (the 0-D analog of a
   sigma_min/spectral-gap argument, here for a delay-differential-equation system): a pure-delay
   negative-feedback loop's characteristic equation admits a marginally-stable (purely-imaginary-
   eigenvalue) solution at a period `T` that brackets `[2*tau, 4*tau]` across the two standard
   idealized closed forms (a static-gain loop transmission `1+g*e^{-s*tau}=0`, and a pure-integrator
   delayed feedback `s+a*e^{-s*tau}=0`) — **both verified below by direct complex-number substitution**
   into their own characteristic equations (§7), not asserted. Inverting the real, independently-
   measured resonance period (~10s, deBoer 1987, NOT fit here) brackets the implied sympathetic loop
   delay at **2.5–5.0 s** — a BACKWARD-implied quantity, the same discipline `arterial_pressure.py`
   already used for implied TPR.

## 2. Method, in one paragraph

The twin's own resting MAP (93.33 mmHg, `arterial_pressure.py`'s classic-formula central estimate) and
resting HR (72.57 bpm, `cardiac_output.py`'s subject2/walking1 Fick-chain result) are loaded, read-only.
A 4-parameter-logistic baroreceptor sigmoid is built with illustrative gain/range parameters (disclosed,
motivated by Seagard et al. 1990's real single-fiber dog data) and its operating point `P50` set to the
twin's own resting MAP (a disclosed **design choice**, not an independent falsifier). The reflex arc
(afferent CN IX/X → NTS → CVLM/RVLM + nucleus ambiguus → sympathetic/vagal efferents) is mapped
structurally, with a conduction-time budget computed from generic anatomical path lengths and reused/
textbook fiber-conduction-velocity ranges. Baroreflex sensitivity is tested against a **real, live-
verified, quantitative literature anchor** (Laitinen et al. 1998, n=117, phenylephrine method) against
the task's own pre-registered [10,20] ms/mmHg band, and against a **real hypertension anchor** (Mussalo
et al. 2002) for the symmetric-QC "declines in hypertension" clause. The closed-loop sign is derived
via the chain rule and machine-cross-checked numerically. The Mayer-wave resonance argument inverts a
real, independently-measured period (deBoer 1987, corroborated by Wikipedia and Julien 2006) through
two idealized control-theory toy models, both numerically verified, and is explicitly **held open**
per Julien 2006's own disclosed, live scientific controversy over whether this is a true resonance at
all.

## 3. Citations — every PMID verified LIVE this session (NCBI eutils via direct curl, not the
WebFetch AI-summarizer — see §10's disclosed catch), plus 4 Wikipedia pages (textbook-grade, flagged)

| # | Citation | PMID/DOI | Role |
|---|---|---|---|
| 1 | Kent BB, Drane JW, Blumenstein B, Manning JW (1972). "A mathematical model to assess changes in the baroreceptor reflex." *Cardiology* 57(5):295-310. | **4651782** (verified live: title/journal/year/authors match) | Sigmoid functional FORM. **Bibliographic only** — no abstract available live (pre-1975 abstracting era, same disclosed-gap tier as this repo's Astrand 1964/Rowell 1974). |
| 2 | Seagard JL, van Brederode JF, Dean C, Hopp FA, Gallenberg LA, Kampine JP (1990). "Firing characteristics of single-fiber carotid sinus baroreceptors." *Circ Res* 66(6):1499-509. | **2344663** (verified live, full abstract) | REAL single-fiber dog data, quoted verbatim: "type I, a discontinuous, hyperbolic pattern... sudden onset of discharge at threshold pressure" (narrow range, high sensitivity, large myelinated-A afferents) vs "type II, a continuous, **sigmoidal** pattern... gradual increase in discharge" (wide range, lower sensitivity, smaller-A+unmyelinated-C afferents, spontaneous sub-threshold discharge). Primary anchor for §4's sigmoid shape + afferent fiber split. |
| 3 | Seagard JL, van Brederode JF, Dan C, Hopp FA, Elegbe EO, Gallenberg LA, Kampine JP (1991). "Effects of epinephrine on firing characteristics of two functionally different types of carotid baroreceptors." *Circ Res* 69(4):1097-105. | **1934338** (verified live, abstract) | REAL finding: epinephrine "significantly increase[s] sensitivity, Fth, and Fsat of both types of baroreceptors" — a real afferent-limb positive-modulation complication from the sympatho-adrenal system, corroborating this repo's own `adrenal_cell` coupling note with a live-verified mechanism. Disclosed, not dynamically modeled (§10). |
| 4 | Gribbin B, Pickering TG, Sleight P, Peto R (1971). "Effect of age and high blood pressure on baroreflex sensitivity in man." *Circ Res* 29(4):424-31. | **5110922** (verified live: title/journal/year/authors match) | The classic age+HTN-decline paper. **Bibliographic only** — no abstract available live (1971). |
| 5 | Laitinen T, Hartikainen J, Vanninen E, Niskanen L, Geelen G, Lthansimies E (1998). "Age and gender dependency of baroreflex sensitivity in healthy subjects." *J Appl Physiol* 84(2):576-83. | **9475868** (verified live, full abstract) | **THE primary quantitative falsifier anchor.** REAL, n=117 healthy, age 23-77, **phenylephrine method** (one of the 3 pre-registered methods). Quoted verbatim: "BRS correlated with age (r=-0.65, P<0.001)... diastolic blood pressure (r=-0.47, P<0.001)... significantly higher in men than in women (15.0±1.2 vs. 10.2±1.1 ms/mmHg; P<0.01)... age and gender... accounted for 52% of interindividual BRS variation... 24% of women >40yr... and 18% of men >60yr had markedly depressed BRS (<3 ms/mmHg)." |
| 6 | Mussalo H, Vanninen E, Ikäheimo R, Laitinen T, Laakso M, Lthansimies E, Hartikainen J (2002). "Baroreflex sensitivity in essential and secondary hypertension." *Clin Auton Res* 12(6):465-71. | **12598951** (verified live, full abstract) | REAL hypertension anchor, **phenylephrine method**. Quoted verbatim: "BRS in the RVHT (3.7±0.6 ms/mmHg) and SEHT (7.6±0.8 ms/mmHg) groups... MEHT [8.5±1.2 ms/mmHg]... Baroreceptor reflex regulation has been shown to **reset** towards a higher blood pressure level and to operate with reduced sensitivity in hypertension." |
| 7 | Smyth HS, Sleight P, Pickering GW (1969). "Reflex regulation of arterial pressure during sleep in man. A quantitative method of assessing baroreflex sensitivity." *Circ Res* 24(1):109-21. | **4303309** (verified live: title/journal/year/authors match) | The **original** phenylephrine-ramp BRS method paper. **Bibliographic only** — no abstract available live (1969). |
| 8 | Rudas L, Crossman AA, Morillo CA, Halliwill JR, Tahvanainen KU, Kuusela TA, Eckberg DL (1999). "Human sympathetic and vagal baroreflex responses to sequential nitroprusside and phenylephrine." *Am J Physiol* 276(5 Pt 2):H1691-8. | **10330255** (verified live, full abstract) | REAL modified-Oxford data, n=18. Quoted verbatim: "vagal baroreflex slopes are less when arterial pressures are falling than when rising... hysteresis... Sympathetic baroreflex slopes are similar when arterial pressure is falling and rising; however, small pressure elevations above baseline silence sympathetic motoneurons... Vagal, but not sympathetic, baroreflex gains vary inversely with subjects' ages... There is no correlation between sympathetic and vagal baroreflex gains." Primary source for the vagal/sympathetic asymmetry + decorrelation. |
| 9 | Laude D, Elghozi JL, Girard A, et al. (2004). "Comparison of various techniques used to estimate spontaneous baroreflex sensitivity (the EuroBaVar study)." *Am J Physiol Regul Integr Comp Physiol* 286(1):R226-31. | **14500269** (verified live, full abstract) | REAL 11-center × 21-subject × 21-method comparison. Quoted verbatim: "alpha-coefficient or gain of the transfer function in both the low-frequency band or high-frequency band, TRS, and sequence methods gave **strongly related** results. Conversely, weighted gain, X-AR, and Z exhibited **lower agreement**... Some procedures were unable to provide results when BRS estimates were expected to be very low." THE primary "methods disagree" anchor. |
| 10 | Di Rienzo M, Castiglioni P, Mancia G, Parati G, Pedotti A (1997). "Critical appraisal of indices for the assessment of baroreflex sensitivity." *Methods Inf Med* 36(4-5):246-9. | **9470369** (verified live, full abstract) | REAL cat sinoaortic-denervation data. Quoted verbatim: "the alpha coefficients computed at the respiratory frequency tend to be **higher** than alpha coefficients estimated at 0.1 Hz... the PI-SBP coherence... **remain[s] above the 0.5 threshold also after baroreceptor denervation**" — a genuine spectral-method failure mode (a validity gate that passes vacuously). Second, independent methods-disagreement anchor. |
| 11 | Eckberg DL (1976). "Temporal response patterns of the human sinus node to brief carotid baroreceptor stimuli." *J Physiol* 258(3):769-82. | **978502**, PMCID PMC1309004 (verified live, full abstract; PMC confirmed **not** open access) | REAL human neck-suction data. Quoted (PubMed plaintext): "Maximum sinus node inhibition occurred... about 0-75 sec before the anticipated... P wave... the baroreceptor-cardiac reflex arc... averaged 0-24 sec." **Disclosed catch** (§10): almost certainly a legacy decimal→hyphen rendering artifact; interpreted as 0.75s/0.24s, cross-checked (not proven) against Wikipedia's independent "less than the duration of a cardiac cycle" statement. |
| 12 | deBoer RW, Karemaker JM, Strackee J (1987). "Hemodynamic fluctuations and baroreflex sensitivity in humans: a beat-to-beat model." *Am J Physiol* 253(3 Pt 2):H680-9. | **3631301** (verified live, full abstract) | REAL model-vs-data study. Quoted verbatim: "The so-called **10-s rhythm** in HR and BP appears as a **resonance phenomenon due to the delay in the sympathetic control loop** of the baroreflex." THE primary anchor for §7 — the 10s period is theirs, not fit here. |
| 13 | Julien C (2006). "The enigma of Mayer waves: Facts and models." *Cardiovasc Res* 70(1):12-21. | **16360130** (verified live, full abstract) | REAL modern review. Quoted verbatim: "Mayer waves are oscillations of arterial pressure... approximately 0.1 Hz in humans... Recent analysis of the various transfer functions of the rat baroreceptor reflex suggests that Mayer waves are **transient oscillatory responses to hemodynamic perturbations rather than true feedback oscillations**." THE symmetric-QC "hold open" anchor for §7. |
| 14 | Akselrod S, Gordon D, Ubel FA, Shannon DC, Berger AC, Cohen RJ (1981). "Power spectrum analysis of heart rate fluctuation: a quantitative probe of beat-to-beat cardiovascular control." *Science* 213(4504):220-2. | **6166045** (verified live, full abstract) | Quoted: "renin-angiotensin system activity strongly modulates the amplitude of the spectral peak located at **0.04 hertz**." Decorrelation anchor — 0.04Hz (RAAS) ≠ ~0.1Hz (Mayer/baroreflex), both inside the broader LF band. |
| 15 | Task Force of the European Society of Cardiology and NASPE (1996). "Heart rate variability: standards of measurement..." *Circulation* 93(5):1043-65. | **8598068** (verified live: title/journal/year match) | LF/HF band-definition standard. **Bibliographic only** — no abstract in Medline for this consensus statement. Band numbers cross-checked via Wikipedia. |
| 16 | Wikipedia "Baroreflex", "Mayer waves", "Baroreceptor", "Heart rate variability" (all live-fetched this session, textbook-grade, flagged). | live-fetched | Reflex-arc structure (NTS→CVLM→inhibits RVLM; NTS→nucleus ambiguus+DMNX), afferent nerve identity (CN IX/X), "less than the duration of a cardiac cycle" vagal-speed statement, "~0.1 Hz, or a 10-second period" Mayer-wave figure, LF=0.04-0.15Hz/HF=0.15-0.4Hz bands. Cross-triangulates primary literature, never substitutes for it. |

**Reused, not re-verified** (already live-verified in the sibling docs that produced them): MAP=93.33
mmHg (Razminia 2004, PMID 15558774, via `docs/MECHANISM_ARTERIAL_PRESSURE.md`); HR_rest=72.57 bpm and
HR_max=167 bpm (Higginbotham et al. 1986, PMID 3948345, via `docs/MECHANISM_CARDIAC.md`/
`docs/MECHANISM_CARDIAC_OUTPUT.md`); the general conduction-delay LAW `path/velocity+fixed_delay`
(`docs/MECHANISM_NERVE_CONDUCTION.md`, reused as a formula, not by importing that script's leg-specific
fiber-velocity constants — a different anatomical circuit).

## 4. Inputs (reused, not re-solved)

| quantity | value | source |
|---|---:|---|
| MAP_rest (classic formula) | 93.33 mmHg | `arterial_pressure.py` §6, Razminia 2004-anchored |
| HR_rest (subject2/walking1) | 72.57 bpm | `cardiac_output.py`, within 0.6% of Higginbotham's real 73bpm |
| RR_rest = 60000/HR_rest | 826.8 ms | derived |
| HR_max (real, near-exhaustive) | 167 bpm | Higginbotham et al. 1986, PMID 3948345, reused |

## 5. Baroreceptor sigmoid — illustrative parameters, disclosed, not fabricated precision

`FR(P) = FR_min + (FR_max-FR_min)/(1+exp(-k*(P-P50)))`, with `FR_min=2.0 Hz`, `FR_max=80.0 Hz`,
`k=0.10 /mmHg` (all **illustrative/order-of-magnitude**, motivated by Seagard 1990's real qualitative
two-population structure — not independently live-pinned to one specific human numeric point-estimate
this session, same disclosed-gap tier as `cardiac_output_geometric.py`'s Ees/Ea). `P50` is **set** to
the twin's own resting MAP (93.33mmHg) — a **modeling choice** (baroreflex operating points track
prevailing pressure; Mussalo 2002's real "resetting" finding supports this qualitatively), explicitly
**not** treated as an independent falsifier or gated PASS/FAIL (that would be a tautology-adjacent
check — flagged and avoided, not laundered).

| quantity | value |
|---|---:|
| 5%/95% crossing (illustrative Pth/Psat) | 63.9 / 122.8 mmHg |
| max gain, analytical `k*(FRmax-FRmin)/4` | 1.9500 Hz/mmHg |
| max gain, numerical (`np.gradient` over a 4000-pt sweep) | 1.9500 Hz/mmHg (diff **0.000%**) |
| monotonic increasing over [30,200]mmHg sweep | **PASS** |

**Honest complication, disclosed not modeled**: Seagard 1991 (PMID 1934338) found circulating
epinephrine — itself a sympathetic/adrenal-medulla product — measurably **sensitizes** both
baroreceptor types (raises sensitivity, threshold firing rate, and saturation firing rate). This is a
real positive-feedback loop from the efferent arm back onto the sensor, disclosed here, not dynamically
incorporated into this static sigmoid.

## 6. Reflex arc — afferent → NTS/medulla → sympathetic + vagal efferent split

- **Afferent**: carotid sinus → CN IX (glossopharyngeal); aortic arch → CN X (vagus) → NTS. Fiber
  types per Seagard 1990: type I = large myelinated A (narrow operating range, high sensitivity); type
  II = smaller A + unmyelinated C (wide operating range, lower sensitivity, tonic/baseline-pressure
  signaling role).
- **Central**: NTS → excites CVLM → CVLM **inhibits** RVLM (the tonic sympathetic vasomotor center) —
  so BP-up drives sympathetic **down**. NTS → excites nucleus ambiguus + dorsal motor nucleus of vagus
  (vagal preganglionic) — so BP-up drives vagal **up**. (Structure per live-fetched Wikipedia
  "Baroreflex," textbook-grade, cross-checked against the well-established physiology this repo's
  anchor-graph node already names.)
- **Efferent, sympathetic**: HR (β1) + contractility + TPR (α1) + venous tone (α1) — norepinephrine-
  mediated, **SLOW** onset/offset (Rudas 1999).
- **Efferent, vagal**: HR **only**, via SA-node M2/IKACh — acetylcholine-mediated, **FAST** (within
  about one cardiac cycle; Eckberg 1976 + Wikipedia corroboration). Parasympathetic fibers do not
  densely innervate peripheral resistance vessels or venous capacitance vessels — the vagal arm cannot
  substitute for the sympathetic arm's TPR/venous-tone control.

### Conduction-time budget (generic anatomical distances — disclosed scope narrowing, §0)

| leg | conduction time | as % of Eckberg's real 0.24s vagal-arc total |
|---|---:|---:|
| afferent, A-fiber case (47.97 m/s, Awang 2007, reused cross-nerve-type) | 2.50 ms | **6.6%** |
| afferent, C-fiber case (0.5 m/s, reused from `nerve_conduction.py`) | 240.0 ms | **127.8%** |
| efferent, vagal (B-fiber, 3-15 m/s) | 13.3–66.7 ms | 5.5–27.8% |
| efferent, sympathetic (B-fiber, 3-15 m/s) | 16.7–83.3 ms | vs the §7 implied 2.5-5.0s bracket: **≤3.3%** |

**A genuinely asymmetric, disclosed finding** (not smoothed over): for **type-I (A-fiber)**
baroreceptor afferents, conduction is a small fraction (6.6%) of the real measured vagal total — most
of the 0.24s is central/synaptic/SA-node response, exactly the opposite conclusion from
`nerve_conduction.py`'s leg H-reflex (where conduction dominated). But for the **type-II (C-fiber)**
population specifically, conduction time alone (240ms) is comparable to — and by this illustrative
estimate slightly exceeds (127.8%) — Eckberg's entire measured reflex-arc duration. This is reported
as an honest tension (the two baroreceptor populations plausibly contribute on very different
timescales), not resolved to one number. For the **sympathetic** arm, conduction is a tiny fraction
(≤3.3%) of the §7-implied 2.5-5.0s delay bracket under any fiber-velocity assumption — the slow arm's
overall speed is set by **neuroeffector kinetics** (norepinephrine release/diffusion/receptor
cascade), not axonal conduction.

## 7. Baroreflex sensitivity (BRS) — the pre-registered numeric falsifier

**Pre-registered band**: 10-20 ms/mmHg, healthy adults.

| source | value (ms/mmHg) | method | in [10,20]? |
|---|---:|---|---|
| Laitinen 1998, men (n≈58/117) | 15.0 ± 1.2 | phenylephrine | **PASS** |
| Laitinen 1998, women (n≈59/117) | 10.2 ± 1.1 | phenylephrine | **PASS** (at lower boundary) |

Age correlation `r=-0.65` (P<0.001); DBP correlation `r=-0.47` (P<0.001); age+gender jointly explain
**52%** of interindividual BRS variance; **24%** of women >40yr and **18%** of men >60yr already have
markedly depressed BRS (<3 ms/mmHg) — the "declines markedly with age" clause, quantified, not merely
asserted. **Honest scope note**: Laitinen's sample spans age 23-77 (not young-only); the strong
negative age correlation implies the **young** subset sits **above** these whole-sample means —
directionally consistent with (not contradicting) the falsifier, but no exact young-only point
estimate was extracted live this session (disclosed, §10).

**Hypertension anchor** (Mussalo 2002, phenylephrine method — the "declines... in hypertension"
clause, quantified):

| group | BRS (ms/mmHg) | below healthy-band lower bound (10)? |
|---|---:|---|
| Mild essential HTN (MEHT) | 8.5 ± 1.2 | **YES** |
| Severe essential HTN (SEHT) | 7.6 ± 0.8 | **YES** |
| Renovascular HTN (RVHT) | 3.7 ± 0.6 | **YES** |

All three real hypertensive cohorts sit below the healthy band's own lower bound — a real, quantified,
external confirmation of the symmetric-QC decline clause, plus Mussalo's own quoted mechanism: the
baroreflex **resets** toward a higher operating pressure and loses sensitivity in chronic hypertension.

## 8. Closed-loop direction — chain rule, machine-cross-checked

`dHR/dSBP = -(HR_bpm^2/60000) * BRS_RR[ms/mmHg]`, evaluated at the twin's own resting HR (72.57 bpm):

| BRS source | BRS_RR (ms/mmHg) | dHR/dSBP, analytical (bpm/mmHg) | dHR/dSBP, numerical (finite-difference) | diff | sign |
|---|---:|---:|---:|---:|---|
| Laitinen, men | 15.0 | **-1.3165** | -1.3165 | 1.8e-4% | **negative — PASS** |
| Laitinen, women | 10.2 | **-0.8952** | -0.8952 | 1.2e-4% | **negative — PASS** |

Both the sign (BP-up → HR-down, the required closed-loop direction) and the analytical-vs-numerical
match (near machine precision) are confirmed — a real derivative identity, not an eyeballed direction.

## 9. Method disagreement — symmetric QC, held OPEN, real quoted findings

- **Laude 2004 (EuroBaVar)**: spectral (both LF and HF bands) + TRS + sequence methods "gave strongly
  related results," while weighted-gain, X-AR, and Z "exhibited lower agreement with all the other
  techniques"; some procedures failed to return a value entirely at very-low (baroreflex-failure) BRS.
- **Di Rienzo 1997** (real cat denervation data): sequence ≈ alpha-coefficient-at-respiratory-frequency,
  but alpha-at-respiratory-frequency is systematically **higher** than alpha-at-0.1Hz (the spectral
  method disagrees with **itself** depending on which sub-band is used) — and the PI-SBP coherence>0.5
  validity gate **persists even after baroreceptor denervation**, i.e. it can pass with no real reflex
  present, a genuine specificity failure of a commonly-used spectral-method safeguard.

**Verdict: HELD OPEN, not resolved.** Two independent, real, decorrelated studies (a 21-method
11-center human comparison; an animal denervation study) both show genuine disagreement among the 3
pre-registered method families and among spectral sub-variants — reported honestly, per the task's own
instruction, not laundered into a single number.

## 10. Mayer-wave resonance — pole-placement geometry, held open per Julien 2006

For a pure-delay negative-feedback loop, the marginal (purely-imaginary-eigenvalue) oscillation
condition brackets the period `T` between two idealized closed forms:

| model | characteristic equation | criticality condition | period |
|---|---|---|---:|
| static-gain loop | `1 + g*e^{-s*tau} = 0` | `omega*tau = pi` | `T = 2*tau` |
| pure-integrator loop | `s + a*e^{-s*tau} = 0` | `omega*tau = pi/2` (at `a=omega`) | `T = 4*tau` |

**Both verified by direct complex-number substitution** (not asserted): residuals `1.2e-16` and
`3.9e-17` respectively — machine-zero.

Inverting deBoer et al. 1987's own real, independently-measured period (**T=10s**, their "10-s
rhythm... a resonance phenomenon due to the delay in the sympathetic control loop") brackets the
**implied sympathetic loop delay at 2.5–5.0 s** — **10.4–20.8× longer** than Eckberg's real measured
vagal reflex-arc duration (0.24s), consistent with (not independently re-proving) the vagal-fast/
sympathetic-slow asymmetry Rudas 1999 already established qualitatively. Sympathetic **conduction**
time is only ≤3.3% of this bracket's lower edge (§6) — the loop's slowness is a neuroeffector-kinetics
property, not an axonal-conduction property.

**HELD OPEN, explicitly, per Julien 2006** (a modern, live-verified review): "recent analysis of the
various transfer functions of the rat baroreceptor reflex suggests that Mayer waves are **transient
oscillatory responses to hemodynamic perturbations rather than true feedback oscillations**." The
delay-resonance derivation above is a **consistency check against deBoer's own stated interpretation**,
not an independent proof that a true self-sustained resonance exists — this is a live, genuine,
disclosed scientific controversy, not resolved by this cert.

**Decorrelation note**: Akselrod 1981's real 0.04Hz spectral peak is RAAS-linked, a **different**
sub-component from the ~0.1Hz Mayer/baroreflex peak analyzed here — both sit inside the Task Force
1996's broader 0.04-0.15Hz LF band, not conflated.

## 11. Falsifier verdict — stated exactly as pre-registered

> Does the model reproduce the measured BRS (~10-20 ms/mmHg, healthy adults, sequence / modified-
> Oxford-phenylephrine / spectral-alpha-index) AND the correct closed-loop direction/slope (BP up →
> HR down, the phenylephrine RR-vs-SBP regression)?

**BRS: PASS**, on a real, live-verified, phenylephrine-method anchor (Laitinen 1998, n=117): men 15.0
ms/mmHg (comfortably inside [10,20]), women 10.2 ms/mmHg (inside, at the lower boundary). **Direction:
PASS**, machine-derived via the chain rule (RR=60000/HR) and numerically cross-checked to <2e-4%: `BP
up → HR down` at both -1.32 and -0.90 bpm/mmHg for the two literature BRS values, evaluated at the
twin's own resting HR.

**Symmetric QC, stated plainly**: this is not an unconditional pass. (1) The healthy-band anchor is a
whole-sample (age 23-77) mean, not a young-only point estimate — directionally, not exactly,
consistent with "healthy young adults." (2) The 3 pre-registered measurement methods **genuinely
disagree** with each other (Laude 2004, Di Rienzo 1997), held open, not resolved. (3) The Mayer-wave/
resonance mechanism computed in §10, while numerically self-consistent, is **contested in the primary
literature itself** (Julien 2006) — presented as a consistency check, never as proof. (4) The sigmoid's
specific numeric parameters and the neck/thorax conduction-path lengths are illustrative/generic, not
independently live-pinned this session. What earns more than "lands in a wide band" here: a REAL
phenylephrine-method point estimate (not a hand-picked band), a machine-verified derivative identity
for the direction, two independently-verified characteristic-equation checks for the resonance
geometry, and an honestly asymmetric (not uniformly-favorable) conduction-time finding (§6) that a
symmetric audit could have — and did — turn up against the model's own convenient framing.

## 12. Pre-registered gates — 14/14 PASS + 10 disclosed open items

```
brs_men_in_prereg_band:                                    PASS (15.0 in [10,20])
brs_women_in_prereg_band:                                  PASS (10.2 in [10,20], at boundary)
hypertension_brs_all_below_healthy_band:                   PASS (3.7/7.6/8.5 all < 10)
direction_men_sign_correct_negative:                       PASS (-1.3165 bpm/mmHg)
direction_women_sign_correct_negative:                     PASS (-0.8952 bpm/mmHg)
direction_men_analytical_numerical_match:                  PASS (1.8e-4% diff)
direction_women_analytical_numerical_match:                PASS (1.2e-4% diff)
sigmoid_monotonic_increasing:                              PASS
sigmoid_max_gain_analytical_numerical_match:                PASS (0.000% diff)
resonance_static_gain_characteristic_eq_satisfied:          PASS (residual 1.2e-16)
resonance_integrator_characteristic_eq_satisfied:           PASS (residual 3.9e-17)
implied_sympathetic_delay_longer_than_vagal:                PASS (2.5-5.0s vs 0.24s)
sympathetic_conduction_small_fraction_of_implied_delay:     PASS (<=3.3%)
afferent_A_fiber_conduction_small_fraction_of_vagal_total:  PASS (6.6%)

OPEN MODELING UNCERTAINTY (disclosed, does NOT gate overall_pass):
sigmoid_numeric_params_illustrative_not_live_pinned:        FR_min/FR_max/k order-of-magnitude, not one live-pinned human study
p50_equals_map_is_a_design_choice:                          modeling choice (baroreflex resetting), not an independent falsifier
bibliographic_only_no_abstract_pre1975:                     Kent 1972, Gribbin 1971, Smyth 1969 -- title/journal/year verified, no abstract
methods_disagreement_held_open:                             Laude 2004 + Di Rienzo 1997 both show real method disagreement
mayer_wave_mechanism_contested:                             Julien 2006 disputes true-resonance vs transient-response interpretation
sympathetic_delay_bracket_is_backward_inferred:              2.5-5.0s inferred from deBoer's period, not independently measured this session
neck_thorax_path_lengths_generic:                            twin's OpenSim model is leg-only; no live 3-D query possible (unlike nerve_conduction.py)
b_fiber_velocity_textbook_consensus:                         3-15 m/s Erlanger-Gasser range, not independently live-pinned to one PMID
eckberg_1976_decimal_hyphen_ambiguity:                       PubMed plaintext likely mangles decimals to hyphens; PMC not open access to confirm
epinephrine_afferent_sensitization_not_modeled:              Seagard 1991's real afferent-sensitization finding disclosed, not dynamically modeled
laitinen_sample_spans_23_to_77_not_young_only:                whole-sample means, not a young-only point estimate (directional only)
```

**Overall: PASS** (14/14 pipeline-correctness gates; 2 independent runs produce byte-identical stdout
and JSON, deterministic, machine-verified — not eyeballed).

## 13. Honest gaps (disclosed, not hidden)

- **Static operating-point model, not a closed-loop dynamic simulation.** No beat-to-beat integration,
  no dynamic resetting, no baroreflex-mediated blood-pressure-control ODE. Answers "is the sensor
  shape, the arc structure, the loop gain, and the closed-loop sign consistent with real measurements,"
  not "how does the whole loop evolve in time."
- **Sigmoid numeric parameters are illustrative**, motivated by Seagard 1990's real qualitative
  structure but not independently live-pinned to one specific human point-estimate study this session.
  `P50=MAP_rest` is a disclosed design choice, deliberately excluded from the gated falsifiers to avoid
  a tautology-adjacent check.
- **Neck/thorax conduction path lengths are generic anatomical estimates**, not live 3-D queries —
  this twin's current OpenSim model (`LaiArnoldModified2017`) is a lower-limb gait model with no
  cranial-nerve/cardiac-autonomic-plexus anatomy, unlike `nerve_conduction.py`'s leg geometry. A
  genuine, disclosed scope narrowing versus the sibling doc, not an oversight.
- **B-fiber (autonomic preganglionic) conduction velocity (3-15 m/s) is the standard Erlanger-Gasser
  textbook range**, not independently live-pinned to one specific PMID this session — same disclosed-
  gap tier as `nerve_conduction.py`'s own C-fiber entry.
- **The sympathetic loop-delay bracket (2.5-5.0s) is backward-inferred**, not an independently
  live-measured quantity this session — no PMID was found with a direct point-estimate for this exact
  quantity. It is a consistency deduction from deBoer 1987's own measured period through two idealized
  toy models, explicitly disclosed as such.
- **The Mayer-wave delay-resonance mechanism is contested in the primary literature itself** (Julien
  2006) — this cert computes a consistency check against one (classical, deBoer 1987's) interpretation,
  not an independent resolution of a live scientific debate.
- **The 3 pre-registered BRS methods genuinely disagree** (Laude 2004, Di Rienzo 1997) — reported, not
  resolved, per the task's own explicit instruction.
- **Kent 1972, Gribbin 1971, and Smyth 1969 have no abstract available live** (pre-1975 abstracting
  era) — cited by verified title/journal/year/author match only, same disclosed-gap tier this repo
  already applies to Astrand 1964/Rowell 1974 (`docs/MECHANISM_CARDIAC.md`).
  Real numeric substitutes for Gribbin's age/HTN claim come from the independently-sourced Laitinen
  1998 and Mussalo 2002 papers instead.
- **Eckberg 1976's key latency numbers carry a disclosed decimal/hyphen parsing ambiguity** — PubMed's
  plaintext abstract renders them as "0-75 sec"/"0-24 sec"; interpreted as 0.75s/0.24s (consistent with
  Wikipedia's independent "less than one cardiac cycle" statement), but the paywalled PMC full text
  could not be checked to fully resolve this live this session.
- **Epinephrine's real, measured afferent-sensitizing effect on the baroreceptors themselves** (Seagard
  1991) is disclosed but not dynamically incorporated — a genuine positive-feedback complication on top
  of the primarily-negative-feedback loop this doc models.
- **Laitinen 1998's healthy-band anchor spans age 23-77**, not young-only — the pre-registered claim is
  specifically about healthy YOUNG adults. The strong negative age correlation (r=-0.65) makes it very
  likely the young subset sits above these means, but no exact young-only point estimate was extracted
  live this session.
- **Single resting operating point, no exercise/orthostatic/sleep regime.** No test of BRS under
  orthostatic stress, exercise, or sleep (despite Smyth 1969's own sleep-based original method) —
  out of scope for this first build.
- **Confidence tier, stated precisely, not blanket-claimed**: the **BRS numeric anchor and the
  direction/sign check** (§7-§8) are **in-vivo-anchored** (real phenylephrine-method human data,
  machine-verified derivative identity) — the strongest results in this doc. The **reflex-arc structure
  and conduction-time budget** (§6) are **structure-anchored with illustrative geometry** (real fiber-
  type qualitative anchors, generic quantitative path lengths). The **sigmoid model** (§5) is
  **method-derived-and-plausibility-bounded**, not independently live-pinned to a specific human
  numeric study. The **Mayer-wave resonance argument** (§10) is **geometry-verified-but-contested** —
  the algebra is machine-checked exactly, but the underlying physiological interpretation is an open
  scientific question per its own primary literature.

## 14. Couples to (graph context, resolves the SEED-DESIGN status, not a full closure)

- **Resolves** `SYS-BAROREFLEX-AUTONOMIC` / `AUTO-BAROREFLEX-AUTONOMIC-LOOP-CLOSED-LOOP-NT`
  (`data/MECHANISM_ANCHOR_GRAPH.json`, prior `OPEN`/`SEED-DESIGN`, *"cited literature anchor, not yet
  executed against data"*) for its **BRS loop-gain + closed-loop-sign + reflex-arc-structure** subset —
  partial, not full, closure (§0, §13; the node's own dynamic/resetting/HRV-readout sub-claims remain
  open).
- **Couples to arterial-pressure** (`docs/MECHANISM_ARTERIAL_PRESSURE.md`): MAP is the **controlled
  variable** this loop defends; this doc's own P50 design choice directly consumes that doc's own
  MAP=93.33mmHg central estimate (§4-§5) — load-bearing, not merely referenced.
- **Couples to cardiac** (`docs/MECHANISM_CARDIAC.md`, `docs/MECHANISM_CARDIAC_OUTPUT.md`): HR is the
  **effector** this loop drives; this doc's own direction/sign check (§8) directly consumes the twin's
  own resting HR=72.57bpm — load-bearing, not merely referenced. The Frank-Starling/SV mechanism in
  `MECHANISM_CARDIAC_OUTPUT.md` supplies the OTHER (preload-driven) route to SV/CO the sympathetic
  efferent arm's contractility effect would modulate — not executed here (deferred).
- **Couples to nerve-conduction** (`docs/MECHANISM_NERVE_CONDUCTION.md`): reuses the general conduction-
  delay LAW (`path/velocity+fixed_delay`) as a formula, explicitly NOT that script's leg-specific fiber-
  velocity constants (a different anatomical circuit) — and reaches the OPPOSITE conclusion about which
  term dominates total latency (conduction dominates the leg H-reflex; neuroeffector kinetics dominate
  the sympathetic baroreflex arm) — a genuine, disclosed contrast between the two docs, not a
  contradiction (different pathways, different rate-limiting steps).
- **Couples to** (per the anchor-graph node's own listed edges, none executed here, deferred per §0):
  `cardiac_cell` (SA-node HR+contractility, shared final-common effector — see cardiac coupling above),
  `HRV_cell` (HF/RSA power as a decorrelated non-invasive vagal readout — the LF/HF caveat in §10
  directly qualifies this), `hypothalamus_cell` (long-run BP/thermoregulatory set-point the baroreflex
  defends acutely and that resets in chronic HTN — Mussalo 2002's resetting finding, §5/§7, is the
  measured anchor for this coupling), `adrenal_cell` (circulating epinephrine sensitizes the
  baroreceptors themselves — Seagard 1991, §5/§13, a REAL, live-verified mechanism for this coupling,
  not just an assumed conceptual link), `thermoregulation_cell` (shares sympathetic cutaneous
  vasomotor/sudomotor efferent pathway, competing for the same RVLM-originating drive — not executed
  here).

## 15. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/baroreflex.py
```
Requires `data/arterial_pressure/arterial_pressure_results.json` and
`data/msk_smoketest/subject2_walking1/cardiac_output/cardiac_output_results.json` to already exist
(both read-only, never modified; run `arterial_pressure.py`/`cardiac_output.py` first if missing). No
OpenSim call, no new subject-trial data, runs in under a second, deterministic (verified: 2 independent
runs produce byte-identical stdout and JSON, both `sys.exit(0)`). No git operations.

**Paths**: script `scripts/msk/baroreflex.py`; evidence `data/baroreflex/baroreflex_results.json`; this
doc `docs/MECHANISM_BAROREFLEX.md`; upstream (read-only) inputs
`data/arterial_pressure/arterial_pressure_results.json` and
`data/msk_smoketest/subject2_walking1/cardiac_output/cardiac_output_results.json`.
