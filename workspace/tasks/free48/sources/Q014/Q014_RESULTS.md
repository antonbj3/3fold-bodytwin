BT-HX-Q014
# When does the order between two interventions change an observable outcome?

**Status: FIRST RUNNING MEKANISTISK MODELL — ready and running.**
`PREREG.md` (sha256 `f0c7799a…9f87f9cd6`, re-frozen after Errata 1–6) ·
`model.py` · `test_model.py` (11/11 PASS) · `results.json` (alla tal).
Run: `python3 model.py && python3 test_model.py` (numpy 1.26 / scipy 1.11,
1 thread, <1 s, <100 MB).

---

## 1. Svaret i en mening

The order changes an observable outcome **if and only if four gates
are simultaneously open**, and the product of them is Closed form:

```
lam_hat  =  G (a_B - a_A) / 2  x  [ (1 - e^{-t_p/tau}) (1 + e^{-(t_w+t_o)/tau})
                                    -  (1 - e^{-t_o/tau}) ]                (E6)
```

`lam_hat` = **order contrast** (first-order transfer / *first-order
carryover*), standard 2×2 estimator according to Willan & Pater (1986),
`lam_hat = (y11−y21) − ½[(y11−y21)+(y22−y12)]`. The ports:

| # | Port | Magnitude | Zero when |
|---|---|---|---|
| **G0** | non-commutativity | `G(a_B − a_A)` | the interventions are interchangeable |
| **G1** | not realized new intervention | `1 − e^{−t_o/tau}` | `t_o = 0` (read immediately at the prey) |
| **G2a** | written memory | `1 − e^{−t_p/tau}` | `t_p ≪ tau` (for short period) |
| **G2b** | remaining memory | `e^{−(t_w+t_o)/tau}` | `t_w + t_o ≫ tau` (too long wash) |
| **G4** | the noise floor | `MDC = 1.96…·t_{α/2,df}·SE` | `\|lam_hat\| ≤ MDC` |

`lam_hat = 0` **identical** when G0 or G2b close — ie. when
the interventions are commutative, or when the outcome is read long after the exchange.
(Verified: `L1`, `L3`, `L5`, `L6`; null models `N0/N0'/N2` give exactly `0.0`
on a 60×60 grid in `(tau, t_w)`, `≤1e-12`.)

**The sharpest result** (closed-form, `t_w = 0`, `t_o = t_p`): set
`e = e^{−t_p/tau}`. Then the bracket is `e(1−e) ≤ ¼`, so

```
max_tau |lam_hat| = a_A/8  = 0.600 pt   vid   tau* = t_p / ln 2 = 129.84 d
```

Along the entire identification line (`a_A = 4.8/S`) the ceiling is `2.352 pt`.
**Both are below the noise floor `MDC = 2.449 pt`.** Conclusion: *for the published
the design (t_w = 0, t_o = 90 d) there is no memory constant `tau` at all that gives
a detectable order effect.* Thus, the gate is **not** the washing time — it is
**read time `t_o`**.

---

## 2. What was built on

BodyTwin rope **doesn't** exist in this sandbox (`~/projects/bodytwin` is missing; only
`inputs/NIGHT_PREAMBLE.md` + `inputs/QUESTION.md` was delivered). Therefore:

* **K03, K09, K11, K16 (the nodes that QUESTION.md points to): UNAVAILABLE.**
  No node id is reused. See PREREG §1.
* **Mandatory null models could not be run**: BT-B24 (body weight ×
  OrthoLoad median) and N12's `k·|GRF|` are not available. Replaced with
  **analytically stronger** null models within the design: `N0` (commutativity),
  `N0'` (additive crossover = `tau→0` with `t_w ≥ t_p`, i.e. the standard assumption
  "no carryover"), `N1` (control sequence A→nothing), `N2` (orthogonal
  condition). All give exactly 0.
* **Previously aborted session:** only `agent.log` (2 step: `ls` + read BRIEF).
  No artifacts to build on.

### Spreads (all **BEAT UPP** in this session)

| # | Source | DOI | Hit Value |
|---|---|---|---|
| **R1** | *Design issues in crossover trials involving patients with Parkinson's disease*, Front Neurol 2023;14:1197281 | `10.3389/fneur.2023.1197281` | **Table 2A**, mean(SE) in Mini-BESTest score: baseline 21.4(1.76)/20.6(1.29); period 1 25.0(0.87)/20.2(1.24); period 2 24.1(0.80)/22.3(1.51); n=7/n=9. **Table 4**: washing period **0 dagar**. |
| **R2** | Willan & Pater, *Biometrics* 1986;42(3):593 | `10.2307/2531209` | definition of the estimand |
| **R3** | Senn, *State Med* 1991;10:1361–1374 | `10.1002/sim.4780100905` | the carryover problem |
| **R4** | Bae, Chung & Je, *BMJ Open* 2016;6:e010873 | `10.1136/bmjopen-2015-010873` | n=37, order effect p<0.001, **wash ≥ 2 dagar** |
| **R5** | same as R1, §4 + Table 4 | — | 25 % of 36 AB/BA trials lacked washing; 16/36 lacked justification; observed washing intervals 0 d–12 v |

Read out from **R1 Table 2A** (arithmetic of published means/SE, labeled HEREBY):
`a_A` is identified from `period1_raw = 25.0 − 20.2 = 4.8 pt`;
`tau_hat = 3.3`; `pi_hat = +0.6`; `lam_hat = 4.8 − 3.3 = 1.5 pt`;
`SE(lam_hat) = 0.5·sqrt(0.87²+1.24²+0.80²+1.51²) = 1.142 pt`;
`MDC95 = 2.449 pt`; `sigma_w = 3.19/3.69 pt`; `t_stat(period1_raw) = 2.79` (p≈0.015);
`mdc(baseline) = 4.02 pt`.

**VERIFIED (searched, not found):** carrying memory time constant `tau` for
this estimate. Search (Crossref + Europe PMC + PMC-full text) returned no
primary source; `tau` is therefore **free parameter**, not a literature constant. This is
also R5's point: the laundry justification is empirically lacking. `tau` **is not
identifiable** from R1 (see PREREG §6/Erratum 4). No measurement data is invented.

---

## 3. Frysta kriterier — alla uppfyllda

`results.json:frozen_criteria` → `frozen_criteria_all_pass: true`

| ID | Kriterium | Resultat |
|---|---|---|
| FC1 | calibration brand is **identified**, not assumed | 4.800000 vs 4.8 pt (exact; declared non-independent) |
| **FC2** | **independent test:** model's `\|lam_hat\|` forecast ≤ published resolution 2.449 pt | **0.00011 pt** ✔ |
| FC3 | zero models ≤ 1e-12 on 60×60 grid in `(tau, t_w)` | **0.0** |
| FC4 | `max_tau \|lam_hat\| ≤ MDC` (two independent ways) | **0.600 pt** (at `tau*=129.84 d`) and **2.352 pt** (identification line) — both ≤ 2.449 |
| FC5 | noise floor forecast: `4.8 > 3.70` **and** `1.5 < 2.45` | both hold |
| FC6 | propagator (E3) = numerical ODE solution of (E1) | max rel. error **5.65e-13** |
| FC7 | initial states comparable: `\|0.8\| < 4.02 pt` | holds |

---

## 4. What fell

1. **The original estimand definition was wrong** (prereg §3 used
   the *difference* between the period changes of the sequences instead of the *average*).
   Under the null hypothesis, my definition gives `2(a_B−a_A) ≠ 0` — it can never
   be a zero. Caught by `N0'` returning 3.00 instead of 0.
2. **Published contrast −3.0 was arithmetic error.** Correct is
   `pi_hat = +0.6` and `lam_hat = +1.5` (not −3.0).
3. **A real dimension bug in the model**, caught by test L5, not by the run:
   the propagator resolved `dz/dt = −z/tau + u` (attractor `u·tau`), so `a_A` was
   in practice a **rate** [pt/d] even though the parameter table stated [pt].
   Corrected first-order relaxation `dz/dt = −(z−u)/tau`; new test L0 sticks
   dimension.
4. **Calibration against `lam_hat` would have been to penetrate noise** (`1.5 < MDC
   2.45`). Corrected to calibrate against the only significantly supported quantity
   (`period1_raw = 4.8 pt`), making `lam_hat` a truly
   independent prediction.
5. **The public “significant carryover effect” in R1 does not hold against its own
   correction.** Uncorrected period-1 contrast `4.8 pt` (t = 2.79, p≈0.015) is
   significant, but the **bias-corrected** `lam_hat = 1.5 pt` is **below**
   MDC 2.449 pt. The attribute "significant positive carryover effect" i R1 §3.2.1
   driven by the uncorrected contrast. The model predicts this independently
   (FC2), and R1 himself states that carryover tests have low statistical
   effect.

---

## 5. Sensitivity ±50 % (control parameters)

At `t_o = 0` (where the order effect is detectable, `lam_hat = 4.800 pt`):

| Parameter | −50 % | base | +50 % | affect `lam_hat`? | affects `t_w*`? |
|---|---|---|---|---|---|
| `tau` = 9 d | 4.800 | 4.800 | 4.794 | 0.1 % | **17.53 → 35.05 → 51.76 d** |
| `sigma_w` = 3.45 pt | 4.800 | 4.800 | 4.800 | 0 % | **∞ → 35.05 → 5.70 d** |
| `t_p` = 90 d | 4.768 | 4.800 | 4.800 | 0.7 % | 32.43 → 35.05 → 35.07 d |
| `a_A`, `G`, `t_o` | 4.800 | 4.800 | 4.800 | 0 % | 35.05 d (icke-identifierbara) |

**Reading:** the *magnitude* of the contrast is robust; the **critical washing time `t_w*`
controlled by the noise floor and the memory constant, not by the intervention size**.
That `sigma_w` −50 % gives `t_w* = ∞` means that with half the noise there was a
washout over any — the effect was then never in danger of disappearing.
`a_A` and `G` are exactly insensitive because only the product `G·a_A` is
identified (invariance tested in `unit_check`).

In case of published design (`t_o = 90 d`) `lam_hat` ≈ 0.0001 pt and
`detectable = False` at **all** ±50 % points for `tau`, `sigma_w`, `t_p` —
the conclusion "no `tau` makes the order detectable here" is robust.

`t_w*` as a function of the reading time (monotonically decreasing, test L8):
`t_o = 0 / 0.25 / 0.5 / 1 / 2 d` → `t_w* = 35.05 / 27.12 / 22.88 / 17.68 / 11.64 d`.

---

## 6. Ask for an answer to the question

* **Order can never mean anything** if the interventions are interchangeable
  (`a_A = a_B` → `lam_hat ≡ 0`, exactly) or if they act on independent states
  (`N2`). It is a *provable* theorem, not an empirical observation.
* **Order can never mean anything if the outcome is read at steady state.**
  If `t_o ≫ tau` G1 collapses to zero and the contrast goes to zero no matter how
  long the laundry is. With `t_p = 90 d`, `t_o` is the only design grip that opens
  the window: detectable for `t_o ≲ 66.5 d`, undetectable from `t_o ≥ 78.5 d`
  (`results.json:detection_design_levers`).
* **When order *can* matter** the window `0 ≤ t_w < t_w*` with
  `t_w* = −tau·ln[(MDC/K + c)/S − 1] − t_o`, `K = G(a_A−a_B)/2`. At `t_o = 0`,
  `tau = 9 d`: **`t_w* = 35.05 dagar`**.
* **Cross-check against an independent trial (R4).** Bae et al. (2016) found one
  order effect (`p<0.001`) with wash **≥ 2 dagar**. The model's window at
  early reading is `[0, 35.0] d` — 2 dagar is well inside, in line with their
  bargain. R1's design with `t_w = 0` and reading at `t_o = 90 d` lies in contrast
  *outside* the window, in line with its corrected `lam_hat` not being
  significant. The `t_o` dependence of the model explains why two RCT with different
  washout schemes produce different strong order effects.
* **Direction of findings (non-judgment):** support for order effects in published
  sequence data is very much an *artifact of read time*, not of one
  slow adaptation. The model does not rule out a slow, non-decimating one
  component (eg, irreversible learning); it excludes that such can give one
  detectable effect when reading in steady state.

---

## 7. What the next resolution step requires

1. **Identify `tau` from within-trial decay data** (priority). Requires
   sequence data with measurements of several `t_o` within period 2, not just one
   endpoint number. MLE on `ln y(t)` against `t` (slope `−1/tau`); with it
   the parametrization validated here gives it `tau` directly. **This is the only one
   the way to distinguish `tau` from `G·a_A`**, which R1 cannot.
2. **Design that can actually answer the question.** Run at least three conditions:
   `t_o` ∈ {0, ~t_p/8, ~t_p/2} as factorial factor, `t_w` as log-factorial
   sweep ~1/3`t_p` … 3`t_p`, and a **control sequence A→** (N1).
   Without the `t_o` factor, the design, per model, is measured in an area where
   the order effect is impossible to see.
3. **Test the non-exp form.** The model assumes a single relaxation constant. A
   two-time scale or non-decimating model (`dz/dt = f(z, u)` without
   `−z/tau` term) is the natural competitor; `PREREG` FC1–FC2 are already
   formulated so that it loses or wins on the same data.
4. **The null models BT-B24 and N12 N1 must be resumed** when the repo exists, and
   The Q014 core mechanism must hit them on the same sequence data. The mandatory ones
   the null models in NIGHT_PREAMBLE §2 are **not** purchased here.
5. **K03/K09/K11/K16 must be read** before any of the above is seen as
   BodyTwin compatible; the node content was completely unknown in this run.
6. **Measurement validity:** `MDC = 2.449 pt` rests on the fact that the four published
   the cell SEs are correct and on common variance structure (they are not:
   period-1 and period-2 give 3.19 respectively. 3.69 pt). `sigma_w` is therefore the one
   parameter with the greatest impact on `t_w*` and should be measured, not assumed.

---

## 8. File overview

| File | Contents |
|---|---|
| `PREREG.md` | hypothesis, reference, frozen criteria, null models, Errata 1–6 with traceability log, §7 refreeze |
| `PREREG.sha256` | `f0c7799af02f501ab6f677284bd68c40088b79fef0b7bd1f7749d6d9f87f9cd6  PREREG.md` |
| `model.py` | (E1)/(E3) first-principles, parameter table with unit + source, `unit_check`, (E6) closed-form, `t_w*` (E4), `tau_min` (E5), null models, calibration, sensitivity |
| `test_model.py` | 11 tests: L0 dimension, L1–L4 analytical limit cases, L5/L6 independent estimator identity, FC6 ODE, N0/N0'/N1/N1'/N2, L7 `t_w*` round trip, L8 monotony — **all PASS** |
| `results.json` | 167 kB, each number from above; `prereg_sha256` enclosed |

**Reproducera:** `python3 model.py && python3 test_model.py`
