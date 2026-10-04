# MECHANISM CROSSBRIDGE MODEL — a sub-cellular mechanistic contraction model, one muscle, first step (2026-07-21)

**Scope, stated up front (do not read past this as more than it is):** this is a **research PROTOTYPE on
ONE muscle** (soleus_r, right leg), a **first step** toward a mechanistic contraction model for the mechanism
biology substrate — explicitly **NOT** a full-body swap of the twin's existing Hill-type muscle model.
Parameter identification is **generic** (canonical/illustrative rate constants + geometry, not subject-specific),
the model is **decoupled from the full-body simulation** (a standalone sarcomere/half-sarcomere, not wired into
the `.osim` musculoskeletal chain), and several numeric choices are explicitly flagged as illustrative rather
than independently re-verified against primary sources this session (see §7). **ISOLATION** (bodytwin
`COORDINATOR.md` §1): `data/msk_models/*.osim` was read-only (loaded via the OpenSim API, never mutated or saved);
no git commit, no git push, performed by this work.

Code: `scripts/msk/crossbridge_contraction_model.py` (self-contained, numpy-only for the mechanistic core;
the "vs Millard" comparisons additionally require the repo's real OpenSim environment,
`.venv-msk/bin/python3`). Evidence: `data/msk_smoketest/crossbridge_contraction_model/crossbridge_contraction_model_evidence.json`
(regenerated every run — a HYPOTHESIS artifact, not a certified cell, per this repo's node-schema conventions).

---

## 1. Why (the gap this fills)

The twin's muscle model everywhere in this repo is `Millard2012EquilibriumMuscle` (verified live: 168 instances
in `data/msk_models/LaiArnoldModified2017_erector_spinae_subject2_scaled.osim`, including `soleus_r`/`soleus_l`).
This is a **Hill-type phenomenological model**: three curve-fits (active-force-length, force-velocity,
tendon/passive-fiber elasticity) reproduce whole-muscle force output but do not represent the underlying
sub-cellular contractile *mechanism* — sarcomere half-length tension generation via filament overlap,
actin-myosin cross-bridge attachment/detachment cycling, or calcium excitation-contraction coupling. This
script builds that missing mechanistic layer for one representative muscle and asks two **symmetric**
questions, exactly as posed: (1) does the mechanistic model reproduce the canonical, falsifiable cross-bridge
signatures (isometric force-length plateau + descending limb; a force-velocity Hill hyperbola *emerging* from
the kinetics, not imposed; an isometric twitch); (2) does it reproduce Hill/Millard's own behavior on the
**same real muscle, same inputs**, and does it add anything Hill **cannot structurally have** — history
dependence. Per the task's own framing, a cross-bridge model that merely re-fits Hill with more parameters and
adds nothing would be a **weak result** — the history-dependence test is the real discriminator, and is
reported honestly below, not as a "we found the expected positive" narrative.

---

## 2. Citations (all live-verified this session — DOI/PMID/PMCID machine-cross-checked, not recalled)

WebSearch was unavailable this session (session-wide quota exhausted — a pre-existing constraint already
flagged in `docs/MECHANISM_MSK_BUILD_PLAN.md`/`docs/MECHANISM_MUSCLE_AUDIT.md`); verification instead used the
Crossref (`api.crossref.org`) and NCBI (`eutils`/`idconv`) **public structured-metadata APIs directly** via
`curl` — the same fallback pattern this repo's `scripts/msk/validate_ankle_force.py` already uses. Two numbers
were pulled as **direct quotes from the primary source's own abstract** (Gordon-Huxley-Julian 1966, via
`WebFetch` on its open PMC page) — genuine primary-source data, not a secondary paraphrase.

| # | Reference | DOI | PMID | Role in this script |
|---|---|---|---|---|
| 1 | Huxley AF (1957). *Muscle Structure and Theories of Contraction.* Prog Biophys Biophys Chem 7:255-318. | `10.1016/S0096-4174(18)30128-8` | 13485191 | **Primary governing physics** — the two-state cross-bridge attachment/detachment PDE and its piecewise-linear rate functions f(x), g(x), implemented directly. |
| 2 | Gordon AM, Huxley AF, Julian FJ (1966). *The variation in isometric tension with sarcomere length in vertebrate muscle fibres.* J Physiol 184(1):170-192. | `10.1113/jphysiol.1966.sp007909` | 5921536 (PMCID PMC1357553) | **External anchor** for the force-length plateau+descending-limb shape. Plateau **2.05–2.2 µm** and ascending-limb "shoulder" at **~1.67 µm** quoted directly from the abstract text this session. |
| 3 | Zahalak GI (1981). *A distribution-moment approximation for kinetic theories of muscular contraction.* Math Biosci 55(1-2):89-114. | `10.1016/0025-5564(81)90014-6` | not found in PubMed's index (Crossref-verified real; honest note in §7) | Literature precedent for a moment-ODE (Q0,Q1,Q2) reduction of the same PDE — **cited, not implemented** (§7). |
| 4 | Zahalak GI, Ma SP (1990). *Muscle Activation and Contraction: Constitutive Relations Based Directly on Cross-Bridge Kinetics.* J Biomech Eng 112(1):52-62. | `10.1115/1.2891126` | 2308304 | Structural precedent for **gating the Huxley attachment rate by a calcium-activation variable** — the Hill-Huxley-hybrid structure implemented here (§3.3), via direct PDE gating rather than their specific moment-ODE mechanics. |
| 5 | Razumova MV, Bukatina AE, Campbell KB (1999). *Stiffness-distortion sarcomere model for muscle simulation.* J Appl Physiol 87(5):1861-1876. | `10.1152/jappl.1999.87.5.1861` | 10562631 | The task's **alternative candidate approach** (multi-state, adds a weakly-bound pre-power-stroke state). **Not implemented** — flagged as the road not taken (§7), directly motivated by the eccentric-force finding in §5.2. |
| 6 | Millard M, Uchida T, Seth A, Delp SL (2013, presented ASME 2012). *Flexing Computational Muscle.* J Biomech Eng 135(2):021005. | `10.1115/1.4023390` | 23445050 (PMCID PMC3705831) | **The exact source paper** for `Millard2012EquilibriumMuscle` — the phenomenological model compared against, on the real `soleus_r` instance, via live OpenSim API queries (not re-implemented). |
| 7 | Arnold EM, Ward SR, Lieber RL, Delp SL (2010). *A Model of the Lower Limb for Analysis of Human Movement.* Ann Biomed Eng 38(2):269-279. | `10.1007/s10439-009-9852-5` | 19957039 | Model-lineage paper for this exact `soleus_r` instance. **Honest gap**: full text paywalled/CAPTCHA-blocked this session (PMC + SpringerLink both blocked `WebFetch`) — the paper's own exact "assumed optimal sarcomere length" number could not be independently pulled; see §3.1 for how this gap is handled. |
| 8 | Felder A, Ward SR, Lieber RL (2005). *Sarcomere length measurement permits high resolution normalization…* J Exp Biol 208(17):3275-3279. | `10.1242/jeb.01763` | 16109889 | Methodology precedent for sarcomere-length-based fiber normalization in this model lineage (referenced by #7's own reference list, confirmed via WebFetch). Context only. |
| 9 | Burkholder TJ, Lieber RL (2001). *Sarcomere Length Operating Range of Vertebrate Muscles During Movement.* J Exp Biol 204(9):1529-1536. | `10.1242/jeb.204.9.1529` | not found in PubMed's index (Crossref-verified real) | Qualitative context: mammalian/human sarcomere operating ranges sit at longer absolute lengths than frog (#2) — consistent with, not a numeric derivation of, the human 2.6–2.8 µm target. |

---

## 3. The model

### 3.1 Sarcomere filament-overlap geometry — derived from geometry, not fit (mandate: geometric thinking)

The fraction of myosin heads with reachable actin ("available" cross-bridges) is the 1-D interval intersection
of the myosin-head zone `[c, b]` (M-line to thick-filament tip, minus the bare zone) and the actin span
`[x-a, x]` (Z-line at half-sarcomere length `x=SL/2`, thin filament length `a`):

```
overlap(x) = max(0, min(x, b) - max(x - a, c))
```

This is **exact** for the plateau + descending limb (`x ≥ a`, no thin-thin double overlap) — precisely the
regime the task asks to validate. Below `x=a` (deep ascending limb) a documented thin-thin interference effect
exists (GHJ1966's own "shoulder" at 1.67 µm, quoted above) that this simple formula does not capture — an
honest, explicit scope limit (not silently smoothed over), and not needed since the task specifies only
"plateau + descending limb."

**Two calibrations of the same formula, both machine-checked (`experiment_force_length`):**
- **Frog cross-check**: `(a,b,c)` solved from the plateau bounds **quoted directly** from GHJ1966 (2.05, 2.2 µm)
  → reproduces those bounds to floating-point precision (a pipeline-correctness check, not a novel claim).
- **Human target** (this task's spec): `(a,b,c)` solved so plateau = **[2.6, 2.8] µm**, bare-zone half-length
  held fixed at an illustrative 0.085 µm (bare-zone length is the most conserved of the three geometric
  parameters — set by thick-filament backbone packing — so it is the natural one to hold fixed rather than fit).
  **Honest note**: the exact human filament-length numbers are canonical/illustrative, chosen to hit the
  operator-specified plateau target, *consistent with* (not independently re-derived from) Burkholder & Lieber
  (2001)'s documented longer mammalian operating range — the Arnold et al. (2010)/Ward et al. (2009) lineage
  paper's own exact assumed value could not be pulled from paywalled full text this session (§7).

### 3.2 Huxley (1957) two-state kinetics

Dimensionless cross-bridge distortion `ξ = x/h` (`h` = power-stroke distance): `ξ∈[0,1]` = attached,
force-generating; `ξ<0` = attached but past the power stroke (must detach fast); `ξ>1` = beyond the attachment
range (cannot newly attach). Piecewise-linear rate functions (Huxley's original form):

```
f(ξ) = f1·ξ   for 0≤ξ≤1, else 0        (attachment)
g(ξ) = g1·ξ   for ξ≥0                   (detachment, still rises past ξ=1 — no separate cutoff)
g(ξ) = g2     for ξ<0                   (fast detachment once past the power stroke)
```

`f1:g1:g2 = 43.3 : 10.0 : 209.0` — the canonical illustrative ratio commonly reproduced in muscle-mechanics
teaching material representing Huxley's 1957 fit; **not independently re-extracted from the original paper's
typeset tables this session** (pre-PMC-era Elsevier archive, paywalled). The absolute magnitude is then
calibrated (§3.5), so the exact historical digits do not matter — what is load-bearing (and undisputed
regardless of precision) is the ratio *structure*: `g2 ≫ f1 > g1`.

**Numerical engine — a deliberate, disclosed implementation choice**: the literature's standard efficient
reduction of this exact PDE is Zahalak's (1981) distribution-moment method (3 ODEs instead of an N-point PDE).
This script instead solves the PDE **directly** via upwind finite-difference with an exact-exponential reaction
split (Strang splitting: the reaction sub-step is solved in closed form, `n_eq + (n−n_eq)·exp(−(gain+loss)·dt)`,
unconditionally stable regardless of rate-constant stiffness) — more implementation steps, but zero
closure-approximation risk, and directly verifiable against Huxley's own closed-form steady-state solution
(§4.1) rather than trusting a hand-re-derived moment closure. This is a first-step engineering choice, not a
claim of having implemented Zahalak's method.

### 3.3 Calcium-activation dynamics (the Hill-Huxley-hybrid layer)

Per Zahalak & Ma (1990)'s structural precedent: a calcium transient gates the attachment rate only (detachment
kinetics are not calcium-gated — the standard assumption in this model class). `Ca(t)` is solved **analytically**
(piecewise-exponential release/uptake ODE, exact — no integration error in the forcing signal): during a
stimulus pulse, `dCa/dt = k_release·(1−Ca) − k_uptake·Ca`; between pulses, `dCa/dt = −k_uptake·Ca`. Thin-filament
activation is a Hill-cooperative function of Ca, `a(Ca) = Ca^n / (Ca^n + K^n)`. Illustrative rate constants
(`k_release=400/s, k_uptake=60/s, K=0.4, n=3`) — not fit to a specific published soleus Ca-transient dataset
this session (§7).

### 3.4 Comparison target — the REAL Millard2012 muscle, not a re-implementation

All "vs Millard" numbers come from **live OpenSim API calls** on the actual `soleus_r`
`Millard2012EquilibriumMuscle` instance (`.venv-msk`, OpenSim `4.6-2026-06-22-85aaf64`):
`ActiveForceLengthCurve.calcValue()`, `ForceVelocityCurve.calcValue()`, etc., evaluated at matched normalized
lengths/velocities — genuinely the same real muscle used elsewhere in this twin, not a re-derivation.
Real, verified parameters: `max_isometric_force=6194.84 N`, `optimal_fiber_length=55.164 mm`,
`tendon_slack_length=352.73 mm`, `max_contraction_velocity=10.0 L0/s` (→ `Vmax=0.5516 m/s`).

### 3.5 Calibration — two deliberately separated free knobs

- **`rate_scale`** (absolute cycling-time magnitude): set so the attachment/detachment time constant
  `τ_xb = 1/((f1+g1)·rate_scale) ≈ 3 ms`, fast relative to calcium kinetics (tens of ms) and the whole twitch —
  the physiologically correct ordering (XB cycling ≪ Ca dynamics ≪ twitch shape), illustrative (§7).
  Result: `rate_scale = 6.254/s`.
- **`velocity_gain`** (the length↔time↔velocity axis scale only): calibrated so the shortening branch reaches
  a **practical 5%-of-isometric threshold** at `norm_v=−1`, matching Millard's convention (its curve is *defined*
  to hit exactly 0 there — live-confirmed: `fvc.calcValue(-1.0)==0.0`). **Not an exact zero-crossing target**,
  because (honest, derived, not asserted): the pure 2-state Huxley model's theoretical F-V curve is
  *asymptotic* — it approaches but never exactly reaches zero force at any finite shortening velocity (shown
  directly by the closed-form solution, §4.1) — a genuine structural difference from Hill's own equation, which
  *does* define an exact finite Vmax. Result: `velocity_gain = 39.83 /µm`.

---

## 4. Validation — pre-registered thresholds, machine-checked, not eyeballed

### 4.1 Self-check: numerics vs. Huxley's own closed-form analytic steady state (PASS, 4/4)

A genuine independent mathematical anchor — the steady-state ODE at constant velocity is solved in closed form
(both `v>0`/lengthening and `v<0`/shortening branches, **across the full domain including the outside-[0,1]
tails**, not just the naive `[0,1]`-only approximation). **OODA note, kept in the record rather than erased**:
the first version of this analytic reference assumed zero population outside `[0,1]`; the self-check initially
**FAILED by 3–40×** (`rel_err` up to 40) and even produced a *negative* total force. Diagnosis (Observe→Orient,
not a one-shot "numerics are broken"): instrumenting the raw `n(ξ)` profile showed the PDE-FD numerics were
actually correct and matched a properly-derived extended analytic tail almost exactly (predicted vs. measured
`n(ξ=3)`: 0.139 vs. 0.1394) — the bug was in the incomplete analytic *reference*, not the PDE. After deriving
and implementing the complete three-region closed form, and separately fixing a real bug (the PDE self-check
inadvertently let sarcomere length drift away from the optimal length during the constant-velocity ramp,
conflating a force-length effect into what should have been a pure-kinetics check):

| norm_v | analytic | PDE-FD | rel. err | |
|---|---|---|---|---|
| −0.60 | 0.06901 | 0.06937 | 0.0052 | PASS |
| −0.30 | 0.14501 | 0.14647 | 0.0101 | PASS |
| +0.30 | 1.42897 | 1.43771 | 0.0061 | PASS |
| +0.60 | 1.38471 | 1.38567 | 0.0007 | PASS |

### 4.2 Isometric force-length: plateau + descending limb (PASS)

- Frog geometry pipeline check (reproduces the GHJ1966-quoted plateau 2.05–2.2 µm exactly): **PASS**.
- Human target plateau **[2.6, 2.8] µm flat** (machine-checked, max−min < 0.02 normalized force across the
  band): **PASS**. Descending limb monotonically decreasing beyond 2.8 µm: **PASS**.
- PDE-vs-analytic isometric cross-check (3 test points): **PASS**.
- **Honest, measured (not assumed) nuance**: the isometric F-L shape is *not* purely a rescaled copy of the
  geometric overlap curve. Steady-state attached fraction at fixed velocity=0 is
  `n_flat(overlap) = overlap·f1/(overlap·f1+g1)` — a *saturating* (Michaelis-Menten-shaped) function of overlap,
  not linear — so the descending limb, expressed as force, is slightly *convex* relative to the exactly-linear
  geometric overlap. Measured: linear-fit R² of the descending limb **in force space = 0.9143** (not 1.0 —
  measured, not asserted). Plateau location/width and the zero-force point are still exactly geometric; only
  the descending limb's fine shape is "softened" by the attachment-kinetics nonlinearity.

### 4.3 Force-velocity: does a Hill hyperbola emerge? (PASS on the pre-registered gate — with a major honest asterisk)

Hill's 1938 equation was never coded into this model — a good fit is a genuine emergent-behavior claim.
**Pre-registered gate: R²>0.95 on the concentric (shortening) branch.**

- Concentric-branch Hill-hyperbola fit: **R² = 0.9987** (`a/F0=0.215, b/Vmax=0.281`) → **PASS**.
- vs. the REAL Millard `ForceVelocityCurve`, split honestly (a single pooled number would hide this):
  - **Concentric (shortening) RMS diff = 0.0537** — good agreement.
  - **Eccentric (lengthening) RMS diff = 1.8804** — poor agreement; peak eccentric force is **2.57× higher**
    than Millard's (which caps at 1.4× isometric, matching real physiology).
- **This over-prediction is a genuine, robust structural finding, not a fixable artifact of one parameter
  choice** — checked directly (OODA, forced before accepting it): an alternative variant with `g(ξ)` *bounded*
  (constant) for `ξ>1` instead of continuing to rise linearly was tested; it made the over-prediction **worse**
  (up to 5.4× isometric, vs. 3.6× for the unbounded version already in the model), confirming the effect is not
  an easy artifact of this one modeling choice. It is a known characteristic of the plain two-state formulation
  (motivating why more detailed multi-state models — Razumova/Campbell 1999, §2 — add extra states) — reported
  honestly, not silently patched.
- PDE-vs-analytic cross-check (3 velocities): **PASS**.

**Reading the two together**: the cross-bridge model reproduces Hill's hyperbola **on the physiologically
dominant shortening side** (where Huxley's 1957 model was originally validated) but substantially
over-predicts on the lengthening side — an honest, scoped, symmetric result.

### 4.4 Isometric twitch (PASS)

Single calcium stimulus at optimal length, full Ca-transient + cross-bridge PDE (not tetanic):
**time-to-peak = 13.2 ms, half-relaxation time = 25.0 ms, twitch/tetanus ratio = 0.684.** Unimodal (clean
rise-then-fall, machine-checked monotonicity on both limbs): **PASS**. Time-to-peak inside a deliberately wide
10–300 ms plausibility window (no independently live-verified soleus-specific twitch-time citation this
session, so no tighter numeric gate is pre-registered — honest, not overclaimed): **PASS**. Ratio < 1: **PASS**.

---

## 5. Same-inputs comparison vs. the real Millard2012 (Experiment 4)

Force-length, matched normalized-length grid, real Millard `ActiveForceLengthCurve` vs. this model's
tetanic-isometric normalized force: **RMS diff = 0.265, Pearson r = 0.954** (good shape correlation, moderate
absolute offset). **Honest, structural, not-a-bug difference**: Millard's AFL curve is a smoothed single-peak
fit (value 1.0 exactly at norm_length=1.0, no flat region, by construction — `OpenSim`/Millard(2013)'s own
curve design). The cross-bridge model instead produces a genuine **flat plateau** from filament-overlap
geometry (§3.1/§4.2) — consistent with the *actual experimentally measured shape* (GHJ1966's real data shows a
true plateau, not a single peak). Both shapes are reported; agreement was not forced.

A segfault was hit and fixed during this comparison's development: OpenSim's SWIG-wrapped curve objects
(`ActiveForceLengthCurve`, etc.) do not by themselves keep their parent `opensim.Model` alive in Python's
reference counting — the model was a local variable in the loader function, and once garbage-collected
(non-deterministic timing; it happened to survive long enough for Experiment 2's calls but not Experiment 4's,
minutes later) the underlying C++ curve objects became dangling pointers and `calcValue()` segfaulted with no
Python exception at all. Diagnosed precisely via `faulthandler` (C-level traceback pointed at the exact
`calcValue` call site) rather than guessed at; fixed by stashing a keep-alive reference tuple
`(model, muscle, curves...)` on the returned curves dict.

---

## 6. History dependence — the real discriminator (Experiment 5)

### 6.1 Residual Force Enhancement / Force Depression: honest, *proven* absence — not a lazy negative

Protocol: tetanic isometric at `L1=1.1` (descending limb) → ramp stretch to `L2=1.5` at 15% Vmax → hold, force
compared to a from-rest tetanic isometric reference at `L2`. Measured: **RFE = −0.085%** (essentially zero,
below the pre-registered 2% detection threshold) → `detected=False`. Symmetric shortening (FDE) protocol:
**+0.005%** → `detected=False`.

**Per the mandate that an honest-negative is not a free pass**: this was *not* accepted at face value. Forced
via OODA — Orient (why would this be zero?): the model's `v=0` steady state is the *unique, globally-stable
fixed point* of a per-`ξ` **linear** ODE (`dn/dt = a·f(ξ)(1−n) − g(ξ)·n`) depending only on the instantaneous
`(activation, SL)` — there is no other state variable that could carry history once `n(ξ)` has re-equilibrated.
This means RFE/FDE are **provably** exactly zero at *true* steady state for this model class — not merely
expected to be small. Verified directly, not just asserted: the residual was measured across increasing hold
durations and **shrinks monotonically toward exactly zero**:

| hold duration | RFE residual |
|---|---|
| 0.25 s | −0.065% |
| 0.50 s | −0.013% |
| 1.00 s | −0.0024% |
| 2.00 s | −0.00038% |
| 4.00 s | −0.00001% |

This confirms the −0.085%/+0.005% headline numbers are the incomplete-relaxation floor of a finite hold
duration, not a genuine (even small) enhancement — a mathematically necessary, not merely measured, negative.
Consistent with the literature's own account of why real RFE requires sarcomere non-uniformity ("popping
sarcomeres") or titin engagement — both explicitly out of scope for this one-muscle, uniform-half-sarcomere
first-step prototype (§7).

### 6.2 Finite relaxation time constant — the real positive discriminator (detected)

Complementary test: a 2% length step at fixed full activation, then hold. The force **does not jump
instantaneously** to its new isometric value — it relaxes with a **fitted exponential time constant of
5.4 ms** (`detected=True`), a direct, unavoidable consequence of having *any* finite attachment/detachment rate
constants.

### 6.3 Millard/Hill idempotency — machine-verified, not argued

The same protocol run through the REAL Millard curves: `F(L,v) = activation·AFL(L)·FV(v)·Fmax`, a pure
function of *instantaneous* `(L,v)` with no internal state. Measured directly (not merely argued from the
equation's form): **RFE_Hill = 0.0000000000%** — exactly zero to floating-point precision, and (mathematically)
the relaxation time constant is exactly zero too, since there is no ODE state to relax. Hill/Millard **cannot**
show RFE, FDE, or a finite relaxation time, **regardless of parameter values** — a structural property of the
model *class*, not a fitting failure.

### 6.4 The honest verdict on the real discriminator

**Mixed, symmetric, and reported as such — not spun positive:**
- RFE/FDE: **absent**, and now *provably* absent for this specific model structure (uniform half-sarcomere,
  plain 2-state kinetics) — a genuine limitation matching known literature explanations for why real RFE needs
  additional structure this first-step prototype does not have.
- Finite relaxation time constant: **present** (5.4 ms), and Hill/Millard structurally cannot have one at all.

The cross-bridge model **does** add one class of history-dependence (rate-limited, non-instantaneous relaxation
after a perturbation) that Hill/Millard cannot represent by construction, and does **not** (in its current
one-muscle, uniform-sarcomere form) add the other (RFE/FDE) — both conclusions are machine-verified, not
asserted, and neither is hidden to make the result look better.

---

## 7. Honest limitations — first-step scope (read before citing this elsewhere)

1. **Generic parameter identification, not subject-specific.** Rate constants, calcium kinetics, and the
   bare-zone geometric parameter are canonical/illustrative (§3.2, §3.3, §3.1), calibrated only against two
   real, verified numbers: the soleus_r `max_contraction_velocity` (10 L0/s, live OpenSim query) and the
   task-specified 2.6–2.8 µm plateau target.
2. **One muscle, decoupled from the full-body sim.** This is a standalone half-sarcomere model, not wired into
   the `.osim` musculoskeletal chain or any joint-force cert in this repo.
3. **Zahalak's (1981) distribution-moment method and Razumova/Campbell (1999)'s multi-state model are cited,
   not implemented** — both are real, correctly-attributed alternatives (§2), deliberately not built this
   session (§3.2, §6.4).
4. **Eccentric force is substantially over-predicted** (§4.3) — a genuine, checked (not hand-waved) limitation
   of the plain two-state formulation at these calibrated rates.
5. **RFE/FDE genuinely absent** in this model's current (uniform-sarcomere, no titin) structure (§6.1) — by
   design/proof, not a bug, but a real scope boundary of what a "first step" can show.
6. **Arnold et al. (2010)'s own exact assumed optimal-sarcomere-length number could not be independently pulled**
   from paywalled/CAPTCHA-blocked full text this session (§2, ref. 7) — the human plateau target instead comes
   directly from the operator's task specification (2.6–2.8 µm), which the geometry was calibrated to hit.
7. **Twitch and calcium-transient timing are illustrative**, not fit to a specific live-verified soleus dataset
   this session (§3.3, §4.4) — the validated claim is the qualitative shape (unimodal, ratio<1, plausible
   order-of-magnitude timing), not a specific millisecond number.
8. **WebSearch was unavailable this session** (quota exhausted) — all citation verification used the Crossref +
   NCBI eutils/idconv public APIs directly via `curl`/`WebFetch`, a documented, repo-precedented fallback
   (§2), not a weaker form of verification for the specific claims it was used for (DOI/PMID/title/author
   existence, and — for GHJ1966 — direct abstract quotation).

---

## 8. Files

- `scripts/msk/crossbridge_contraction_model.py` — the model + all 5 experiments + self-check. Run with
  `.venv-msk/bin/python3 scripts/msk/crossbridge_contraction_model.py` for the full live-Millard comparison
  (falls back to a cached real-parameter snapshot, clearly labeled, for the mechanistic-core-only sections if
  `opensim` is not importable).
- `data/msk_smoketest/crossbridge_contraction_model/crossbridge_contraction_model_evidence.json` — full
  machine-readable evidence (regenerated each run).
- This doc.
