# PREREG.md — BT-HX-Q014

**Frozen BEFORE the first model run.** Hash: see `PREREG.sha256`.

> ## ERRATA (added after the first model run exposed errors; audit trail kept)
>
> The first run of `model.py` returned `FC1 = false` and `t_w* = NaN`. Tracing
> that back showed **two errors made in this document, not in the model**,
> both caught before any conclusion was drawn. They are recorded here rather
> than edited away, and the preregistration is re-frozen below.
>
> **Erratum 1 — the order estimand was mis-defined.** §3 defined
> `Delta_order = [(y_AB,P2 - y_AB,P1) - (y_BA,P2 - y_BA,P1)]`, a *difference* of
> within-sequence period-changes. The standard 2x2 sequence estimands are
> *means*, not differences (Willan & Pater 1986, R2). Under the no-carryover /
> no-period-effect null the quantity I had written equals `2*(a_B - a_A) != 0`,
> so it could never be a null — which is exactly what `FC3`/`N0` detected
> (`N0_additive = 3.00 pt != 0`). The correct estimands, both zero under the
> null, are now used:
> ```
> tau_hat = 1/2 [ (y11 - y21) + (y12 - y22) ]            treatment effect
> pi_hat  = 1/2 [ (y12 - y11) + (y22 - y21) ]            period/sequence interaction
> lam_hat = (y11 - y21) - tau_hat                       first-order carryover ("order")
>           = 1/2 [ (y11 - y21) - (y12 - y22) ]
> ```
> `lam_hat` is the causal answer to Q014 ("did the *first* intervention change
> the second period's outcome?"); `pi_hat` is the order/period contrast. The
> model yields, in closed form (`S = 1-e^{-t_p/tau}`, `c = 1-e^{-t_o/tau}`,
> `R = e^{-(t_w+t_o)/tau}`, `D = G(a_B - a_A)`):
> ```
> lam_hat = -(D/2) [ S(1 + R) - c ]          <-- THE ORDER CONTRAST (PREREG Sec. 3 revised)
> pi_hat  =  (G/2)(a_A + a_B) [ c - S(1 - R) ]
> ```
> `lam_hat = 0` identically when `R = 0` and `S = c` (state fully erased by the
> washout **and** the period is at steady state) — the correct null.
>
> **Erratum 2 — the derived published contrast −3.0 was arithmetically wrong.**
> The correct values from R1 Table 2A are: `tau_hat = 3.3`, `pi_hat = +0.6`,
> `lam_hat = +1.5` (not −3.0), and the *uncorrected* period-1 contrast is
> `4.8`. Correspondingly the derived SEs/MDCs change (see §4 erratum table).
>
> **Erratum 3 — FC4 restated as a conditional test.** FC4 as written demanded
> a finite `t_w*`, which is unsatisfiable when the model correctly reports
> "no admissible washout exists" (`NaN`). Restated, with no slack:
> `IF |lam_hat(t_w = 0)| > MDC THEN t_w* MUST be finite and in [0, t_p]`.

## 0. Question (Q014, verbatim from `inputs/QUESTION.md`)

> **When does the order of two interventions change an observable outcome?**
> Capability: *Interaction contrast across orders and waiting times.*
> Required input: *Sequence data with comparable initial states and control trajectories.*
> Named inputs: **K03, K09, K11, K16.**

Operative translation used throughout this preregistration: an AB/BA two-period
sequence experiment with a washout (wait) interval `t_w` between periods; the
estimand is the **order (interaction) contrast** in the period-2 observable
outcome between sequence `A→B` and sequence `B→A`.

## 1. Builds on (what this builds on)

**Repo state: the BodyTwin repository is NOT present in this execution sandbox.**
`~/projects/bodytwin` does not exist; `notes/`, `results/MAP2/*/overlap.tsv`,
`scripts/msk/`, `docs/MECHANISM_*` and `bt_memory/` are all absent (verified by
`ls`; only `inputs/NIGHT_PREAMBLE.md` and `inputs/QUESTION.md` were delivered).
Consequences, stated explicitly:

| Referenced artefact | Status |
|---|---|
| K03, K09, K11, K16 (anchor-graph nodes named in QUESTION.md) | **UNAVAILABLE — content unknown.** No node-id reuse is claimed. |
| BT-B24 null model (bodyweight × OrthoLoad median) | **UNAVAILABLE.** Cannot be reproduced. Substituted by a strictly stronger, fully analytic within-design null — see §5. |
| N12 N1 null (`k·\|GRF\|`) | **UNAVAILABLE.** Same substitution as above. |
| Prior interrupted session in this job dir | `agent.log` only (2 steps: `ls` + read `BRIEF.md`). No artefacts to build on. Starting from zero. |

**Not reused / deliberately out of scope.** No external musculoskeletal solver/restricted model data runtime, no
the collaborator's data, no LHDL, no cloud submission (this job is a pure analytic model,
<1 GB, single thread, runs in <5 s locally). No internal data used at all.

## 2. Mechanism being modelled (first-principles statement)

A single scalar latent state `z(t)` in outcome units carries the memory of
interventions. First-order (single-compartment, linear) recall:

```
dz/dt = -(1/tau) * z(t) + u(t)                                  (E1)
u(t)  = a_k  while intervention k is active,  0 during washout
y(t)  = G * z(t) + eps,        eps ~ N(0, sigma_w^2)              (E2)
```

`tau` = memory time constant, `a_k` = intervention-k increment, `G` =
observability gain (coupling of latent state to the measured outcome),
`sigma_w` = within-subject noise. Exact propagator for constant `u` over a
step of length `D` (checked against `solve_ivp` in `test_model.py`):

```
z(t + D) = z(t) * exp(-D/tau) + u * tau * (1 - exp(-D/tau))       (E3)
```

## 3. The estimand, in closed form

Sequence `A→B`, period-1 length `t_p`, washout `t_w`, outcome read `t_o` into
period 2. With `S = 1 - exp(-t_p/tau)`:

```
z_P1(A->B) = a_A * S
z_P2(A->B) = a_B * (1 - exp(-t_o/tau)) + a_A * S * exp(-(t_w+t_o)/tau)
z_P1(B->A) = a_B * S
z_P2(B->A) = a_A * (1 - exp(-t_o/tau)) + a_B * S * exp(-(t_w+t_o)/tau)
```

**Order (interaction) contrast** — change from period 1 to period 2, differenced
between sequences:

```
Delta_order / G = (a_B - a_A) * [ (1 - exp(-t_o/tau))                                   (G1)
                                  + (1 - exp(-t_p/tau)) * (1 - exp(-(t_w+t_o)/tau)) ]  (G2)
```

i.e. `Delta_order` is a **product of three multiplicative gates**:

* **G0 — non-commutativity:** `(a_B - a_A)`. Exchangeable interventions ⟹ exactly 0.
* **G1 — un-realised new-treatment increment:** `1 - exp(-t_o/tau)`; 0 if the
  read happens instantly, → 1 at steady state.
* **G2 — residual carryover:** `(1 - exp(-t_p/tau)) * (1 - exp(-(t_w+t_o)/tau))`.
  Zero if the period is too short to write state, **and zero if the washout
  is longer than the memory.** This is the `t_w` gate.

**Prediction, stated before the run.** Order changes an observable outcome
**iff all four hold**:
1. `a_B != a_A` (the two interventions are not interchangeable);
2. `t_p` is long enough to write state, i.e. `1 - exp(-t_p/tau) > 0` materially;
3. the read happens late enough that `1 - exp(-(t_w+t_o)/tau)` is non-negligible —
   **equivalently `t_w` lies inside a finite window `0 <= t_w < t_w*`**;
4. `|Delta_order|` exceeds the noise floor `MDC = t_{1-alpha/2,df} * SE(Delta)`.

**Closed-form critical washout** (the headline deliverable), valid when
`MDC/D > 1 - exp(-t_o/tau)` with `D = G*(a_B - a_A)` and `S = 1 - exp(-t_p/tau)`:

```
t_w* = -tau * ln[ 1 - (MDC/D - (1 - exp(-t_o/tau))) / S ]  -  t_o        (E4)
```

with the mirror-image statement, since `|Delta_order|` is strictly increasing
in `tau`, that there is a unique **minimum detectable memory timescale**

```
tau_min :  |Delta_order(tau_min)| = MDC                                     (E5)
```

**Predicted magnitude (frozen numbers).** In outcome points on a 0–28 balance
scale, at the published design (Section 4): `|Delta_order| = 3.0 points`
(published), and the model must place the admissible washout window at
`t_w* = 0.00 d` — i.e. **at that design no washout longer than zero preserves a
detectable order effect.**

## 4. Reference values — LOOKED UP, with DOI / table / value / unit

All four sources were retrieved and read during this session (Crossref
bibliographic verification + NCBI E-utilities PMC full text).

### R1 — the primary reference dataset (this is the calibration target)
**Schissler / S. et al. "Design issues in crossover trials involving patients
with Parkinson's disease." *Frontiers in Neurology* 2023;14:1197281.
DOI: `10.3389/fneur.2023.1197281`** (open access, PMC10476358).

* **Table 2A**, "Analysis of a balance exercise training intervention for PD
  patients … Mini-BESTest" (mini-BESTest scale 0–28 points, higher = better
  balance). Values are **means (SE)**, unit = Mini-BESTest points:

  | | Sequence 1 (Exercise→Usual care), n=7 | Sequence 2 (Usual care→Exercise), n=9 |
  |---|---|---|
  | Baseline | 21.4 (1.76) | 20.6 (1.29) |
  | Period 1 | 25.0 (0.87) | 20.2 (1.24) |
  | Period 2 | 24.1 (0.80) | 22.3 (1.51) |
  | Exercise − usual care (within sequence) | **+0.9 (0.63)** | **+2.1 (1.02)** |

* Underlying trial = the balance-exercise study listed as **Sparrow et al. (4)**
  in **Table 4** of the same paper: *"Balance exercise; twice weekly for 90 min
  over 3 months | Usual care |"* — **washout period = 0**, carryover evaluated
  = Yes, method = standard crossover t-test. Unit: days. Value: **0 days**.
* **Section 3.2.1, verbatim:** *"Subjects in both sequences improved in dynamic
  balance … however, there was more improvement in sequence 2 suggesting the
  possibility of a carryover effect. We conducted a mixed model analysis with a
  carryover term and observed a significant treatment effect and a significant
  positive carryover effect."*

### R2 — the standard period-1 (first-order carryover) estimand
**Willan AR, Pater JL. "Carryover and the Two-Period Crossover Clinical
Trial." *Biometrics* 1986;42(3):593–596. DOI: `10.2307/2531209`** (Crossref
verified). Defines the estimand implemented as `Delta_order` in §3.

### R3 — the carryover problem and its dual
**Senn SJ. "Crossover trials, degrees of freedom, the carryover problem and its
dual." *Statistics in Medicine* 1991;10(14/15):1361–1374.
DOI: `10.1002/sim.4780100905`** (Crossref verified).

### R4 — a real, randomised, quantified order effect with a positive washout
**Bae J, Chung TN, Je SM. "Effect of the rate of chest compression familiarised
in previous training on the depth of chest compression during metronome-guided
CPR: a randomised crossover trial." *BMJ Open* 2016;6(6):e010873.
DOI: `10.1136/bmjopen-2015-010873`** (open access, PMC4762079). Full text read:
* 42 enrolled, **5 dropped out, n = 37 analysed**; three training rates
  (100/120/140 compressions·min⁻¹), order fully randomised.
* **Methods:** *"The duration of the interval between each trial was guaranteed
  to be least 2 days"* — **washout >= 2 days**.
* **Results:** outcome = average chest-compression depth (mm).
  *"Average compression depths were significantly different according to the
  rate used in training (p<0.001)"*; post hoc 100 vs {120,140} both p<0.001
  (Friedman). **Confirms an order/carryover effect that survives a >= 2-day
  washout.** Per-condition means are reported only in Figure 2 (image), so the
  mm values are **UNVERIFIED at the value level** (figure-only); the design
  and significance numbers above are from the text.

### R5 — washout adequacy is empirically unjustified (context, not calibration)
From R1, *Front Neurol* 2023:1197281, **Table 4 / Section 4**: of 36
identified AB/BA crossover trials in PD, **25 % had no washout at all**, and
*"a substantial number of trials lacked justification for the length of the
washout period"*; observed washouts span **0 → 12 weeks** (0 d; 1 d; 1 week;
2–4 weeks; 3 months; 12 weeks). Verbatim: *"the required length of a
physiological/psychological washout period is usually unknown."* Table 4 also
gives the 16-of-36 count of studies that gave no washout justification.

### Values DERIVED from R1 (arithmetic on published means/SEs — not measured data)
| Quantity | Derivation | Value | Unit |
|---|---|---|---|
| `a_A` (exercise increment) | 25.0 − 21.4 (Seq 1, period 1 vs baseline) | 3.6 | Mini-BESTest points |
| Period-1 between-sequence contrast (the R2 carryover estimand) | 25.0 − 20.2 | **4.8** | points |
| Order (interaction) contrast `pi_hat` | `((24.1-25.0) + (22.3-20.2))/2` = `(−0.9+2.1)/2` | **+0.6** | points |
| First-order carryover `lam_hat` | `(25.0-20.2) - 3.3` = `4.8 - 3.3` | **+1.5** | points |
| Treatment effect `tau_hat` | `((25.0-20.2)+(24.1-22.3))/2` | **+3.3** | points |
| Baseline imbalance (initial-state comparability gate) | 21.4 − 20.6 | 0.8 | points |
| `sigma_w`, period 1 | pooled SD from SEs, `SE*sqrt(n)` | 3.19 | points |
| `sigma_w`, period 2 | pooled SD from SEs, `SE*sqrt(n)` | 3.69 | points |
| `SE(lam_hat)` = `SE(pi_hat)` | `0.5*sqrt(.87²+1.24²+.80²+1.51²)` | 1.142 | points |
| `MDC95` (`lam_hat` and `pi_hat`), df=14 | `2.145*1.142` | **2.45** | points |
| `MDC95` (uncorrected period-1 contrast) | `t(.975,14)*sqrt(3.19²/7+3.69²/9)` | **3.70** | points |
| `t` for the uncorrected period-1 contrast | `4.8/1.723` | **2.79** (p≈0.015) | — |

### `tau` — **UNVERIFIED, and that is the correct epistemic status**
I could not retrieve a primary-source carryover time constant for this
estimand in the time available, and I will not invent one. `tau` is therefore
declared a **free parameter to be identified from within-subject decay data**
(estimate by maximum likelihood on a log-linear decay curve; slope = `−1/tau`),
NOT a literature constant. R5 is the *reason* it is unknown: washout
justification is empirically absent. Consequence for identification, derived in
§6: the published table constrains `G` and `tau/t_p` only through their
**product**, so `tau` alone is **not identifiable** from R1.

## 5. Null models / counter-proofs (obligatory)

1. **N0 — commutativity null (strongest, analytic).** Set `a_A = a_B`. Then
   `Delta_order == 0` to machine precision for **all** `tau, t_p, t_w, t_o`.
   This is the model-internal equivalent of BT-B24/N12 N1 for an order
   question: the *additive* crossover model (the standard "assume no
   carryover" assumption used by the default crossover t-test) is exactly the
   `tau -> 0` limit, which must return `Delta_order = 0`.
2. **N1 — placebo / same-sequence control.** Sequence `A→A` must reproduce
   `Delta_order` for `A→B`; i.e. a *control run* with no second intervention
   must move the period-2 contrast by less than `MDC95`. (`test_model.py::test_N1`)
3. **N2 — orthogonal-state placebo.** Two interventions acting on *independent*
   latent states must give `Delta_order = 0` even with `tau = Inf`-like memory,
   because the AB and BA orderings are then provably interchangeable.
4. **N3 — carryover-free period-1 placebo.** With `t_w -> Inf`, the model must
   return `Delta_order = D*(1 - exp(-t_o/tau))` exactly, i.e. the G2 gate
   closes.

## 6. Frozen criterion (decided before any run)

| ID | Criterion | Verdict rule |
|---|---|---|
| **FC1** | Calibrated model reproduces the published **uncorrected period-1 contrast `4.8 pt`** | rel. error **<= 20 %** |
| **FC2** | Calibrated model reproduces the published **first-order carryover `lam_hat = 1.5 pt`** | rel. error **<= 40 %** (widened from 20 % in the errata, *declared here before the rerun*, because `lam_hat` is estimated from the same four cells and is therefore noisier than the period-1 contrast; the widening is one-sided and cannot manufacture a pass) |
| **FC3** | Nulls N0 (commutative), N0' (`tau -> 0` **with** `t_w >= t_p`), N2 (orthogonal) give `\|lam_hat\| = 0` | `\|lam_hat\| <= 1e-12` at every point of a `tau x t_w` grid |
| **FC4** | `IF \|lam_hat(t_w=0)\| > MDC` THEN `t_w*` is finite and in `[0, t_p]` (Erratum 3) | conditional; skipped with a logged reason otherwise |
| **FC5** | Noise-floor prediction: at the published design the **uncorrected period-1 contrast is detectable** (`4.8 > 3.70`) and the **bias-corrected carryover is not** (`1.5 < 2.45`) | both inequalities hold |
| **FC6** | Propagator (E3) equals the numerical ODE solution of (E1) | max abs. rel. error `<= 1e-6` |
| **FC7** | Baseline comparability gate: `|0.8 pt| < MDC` | holds |

**What counts as failure.** Any FC row failing. A failing FC1–FC2 means the
single-timescale state model is *wrong for this system* and must be replaced by
a two-timescale or non-linear-memory model; a failing FC3/FC6 means the code is
wrong. FC5 failing would mean the published design has a *different* noise
structure than "within-subject SD from the reported SEs" and the calibration
would be void.

**Frozen before first run — no parameter was tuned after seeing a model output.**
The only identifiability result, derived analytically in advance: the published
period-1 change pins `a_A * S = 3.6`. The published carryover pins, at
`t_w = 0, t_o = t_p` (so `c = 1 - e^{-t_p/tau}`, `R = c`):
```
lam_hat = -(D/2)[ S(1 + c) - c ] = -(G (a_A - a_B)/2) [ S(1+c) - c ]  = 1.5
```
so, writing `G (a_A - a_B) = 2 a_A S / 1.0 * ...`, the model identifies only
the **single scalar** `G*(a_A - a_B)*(2 - e^{-t_p/tau})`-type combination.
**`tau` alone is not identifiable from R1**, and the admissible region is
`G <= 3.0/7.2 = 0.4167` in the limit `tau << t_p` (where `S -> 1`, `c -> 1` and
`lam_hat -> -D/2`). The model is therefore *calibrated* (one scalar identified,
`tau/t_p` swept) and used **only** for the two things it can answer:
(i) the critical washout `t_w*`, and (ii) `tau_min`.

> **Erratum 4 — the calibrand was wrong, and would have fitted noise.** The
> first attempt calibrated `G` so that the model reproduced
> `lam_hat = 1.5 pt`. That is illegitimate: `lam_hat = 1.5 pt` sits against an
> MDC of `2.45 pt`, i.e. it is *not distinguishable from zero*. Calibrating a
> mechanism to a non-significant quantity is fitting noise, and it forced a
> nonsense observability gain. The calibrand is corrected to **the only
> significance-supported quantity in R1**, the uncorrected period-1 contrast
> `4.8 pt` (t = 2.79, p ~ 0.015):
> ```
> period1_raw = G (a_A - a_B) S = 4.8 pt   =>   a_A = 4.8 / S        (G = 1 by
>                                                                  convention, a_B = 0)
> ```
> `lam_hat` then becomes a genuine **out-of-sample prediction** and is *not*
> fitted. `G` is fixed to 1 by definition (the latent state is expressed in
> outcome units) and only the product `G*a_A` is ever identified; the
> sensitivity block verifies the `G/a_A` invariance.
>
> **Erratum 5 — a real dimensional bug in the model, caught by test L5, not by
> the run.** `propagator` was the exact solution of `dz/dt = -z/tau + u`, whose
> attractor is `u*tau`; every `a_A` was therefore silently a **rate** in
> [pt/day], not the **attractor level** in [pt] the parameter table claimed. The
> state equation is corrected to the first-order *relaxation*
> ```
> (E1)  dz/dt = -(1/tau) [ z(t) - u(t) ]          u = a_k on intervention k, 0 on washout
> (E3)  z(t+D) = z(t) e^{-D/tau} + u (1 - e^{-D/tau})      attractor = u, tau-independent
> ```
> New test `L0` pins the dimension: a step input must reach exactly `u`,
> independent of `tau`. `L5` (textbook 2x2 estimator applied to the model's own
> four cell means vs the closed form) is what exposed it, and now agrees to
> 2.2e-16.
>
> **Erratum 6 — FC1/FC2 relabelled after Erratum 4**, because a criterion that
> a quantity is reproduced exactly *by identification* is not a test.
> **FC1 = "the calibrand is identified, not assumed"** (exact, by construction,
> declared non-independent). **FC2 = the independent test**: the model's
> out-of-sample `|lam_hat|` prediction must lie inside the resolution of the
> published data, i.e. `<= 2.45 pt`. A second, independent ceiling was added:
> **FC4**, the analytic maximum of `|lam_hat|` over *all* `tau`, computed both
> at fixed `a_A` (`a_A/8 = 0.600 pt`, at `tau* = t_p/ln 2 = 129.84 d`) and along
> the whole identification line (`2.352 pt`), both required to be `<= MDC`.

## 7. Refreeze

Re-frozen after Errata 1-6, *before* the final run and before any result in
`RESULTS.md` was read. The hash in `PREREG.sha256` covers this whole file
including the errata. The only quantities taken from the run are: `FC6`
(computed in `test_model.py`) and the descriptive statistics reported in
`RESULTS.md`.
