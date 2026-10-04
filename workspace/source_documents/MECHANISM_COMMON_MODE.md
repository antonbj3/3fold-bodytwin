# MECHANISM COMMON-MODE EXPOSURE PROBE — quantifying the shared-SO-solve risk the trust ledger named (2026-07-21)

**Question.** `docs/MECHANISM_TRUST_LEDGER.md`'s cross-cutting caveat (read once above §1, restated
in §10) names a real structural risk: `metabolic_cost`, `muscle_fatigue`, `tendon_elastic_energy`,
`validate_emg_timing`, `muscle_perfusion`, `thermoregulation`, `respiratory`, `cardiac_output`, and
`motor_unit_recruitment` all sound like independent physiology layers, but most of them read the
SAME single subject2/`walking1` Static-Optimization (SO) solve
(`data/msk_smoketest/subject2_walking1/static_optimization/so/`). Their mutual agreement is
internal consistency, not external validation. **This doc does not fix that** (only a second
independent solve or real EMG-drive does — `docs/MECHANISM_EMG_DRIVEN.md` is this repo's one
genuine decorrelated input). It **quantifies the exposure**: inject two known, controlled
perturbations into the SO output, propagate each through every downstream layer's own real code
(imported, never reimplemented), and measure Δ-downstream-headline per Δ-SO-input.

**Bottom line, stated before the detail so it cannot be missed: "how wrong would the stack be if
the one SO solve is off by X%" has NO single answer.** Measured sensitivity spans more than three
orders of magnitude across the 9 layers — from **exactly 0** (two layers structurally immune) to a
**~15× amplifier** (one layer). Whether a downstream number is safe or dangerous under a shared-SO
bias is governed entirely by that layer's own internal nonlinear structure (a conservation ceiling,
a convex material curve, a fixed additive floor, a scale-invariant normalization) — not by "how
much SO error" there is in the abstract.

## Method

Two perturbations, both built on the SO solve's own 158-frame grid, both consumed by every
downstream layer through that layer's **own, unmodified functions** (imported, called with the
perturbed arrays substituted in-memory for the real SO output — never a hand-rolled reimplementation
of any model):

- **(A) Uniform +20% activation bias** — every one of the 80 real Muscle activation columns in
  `walking1_StaticOptimization_activation.sto` scaled ×1.20 (clipped to [0,1] — realized exactly
  20.00% mean change, **0 of 12,640 entries** hit the ceiling, so this bias is clean and undamped by
  clipping); the matching `_force.sto` muscle columns scaled by the same *realized*, per-entry,
  post-clip factor, so activation and force stay mutually consistent. Models "what if the whole SO
  solve is uniformly too strong/weak."
- **(B) Realistic, EMG-calibrated co-contraction offset** — not an invented number.
  `docs/MECHANISM_EMG_DRIVEN.md` already built two method-controlled forward-integration
  constructions on this exact trial: `naive_fwd_only` (SO's own smoothed activation as excitation,
  zero real EMG) and `emg_hybrid` (7/80 muscles driven by real measured EMG instead). Their
  difference isolates the true co-contraction *content* SO misses (the primary, decorrelated
  comparison that doc's own method section identifies), and was independently found there to raise
  the knee joint contact force **+48.4%** (material) and change the hip **−12.2%** (modest, opposite
  direction) at SO's own peak instant. This probe adds that **same real delta**
  (`data/msk_smoketest/subject2_walking1/emg_driven/_diag_cache.npz`, `a_hybrid − a_naive_only` /
  `f_hybrid − f_naive_only`) onto the SO baseline for the 7 EMG-covered muscles only (zero elsewhere
  — disclosed: no real-EMG content exists to inject for the other 73/80 muscles).

**Six independent baseline cross-checks confirm the reuse is faithful** (not a rewrite that happens
to look plausible) before any sensitivity number is trusted:

| layer | this probe's baseline | already-published number | match |
|---|---:|---:|---|
| metabolic_cost | COT_net = 7.6834 J/kg/m | 7.683 J/kg/m (`MECHANISM_METABOLIC_COST.md`) | exact |
| tendon_elastic_energy | Σ(soleus_r+gasmed_r+gaslat_r own-peak) = 11.1803 J | 11.180 J (`MECHANISM_TENDON_ELASTIC.md` §5) | exact |
| muscle_fatigue | worst-muscle 60-min decline = 7.964% (`glmed1_r`) | 7.96% at `glmed1_r` (`MECHANISM_MUSCLE_FATIGUE.md` §7) | exact |
| validate_emg_timing | 0 PASS / 3 PARTIAL / 5 FAIL | 0/3/5 (`MECHANISM_EMG_TIMING.md`) | exact |
| muscle_perfusion | peak settled flow = 131.72 ml/min/100g (`glmed1_r`) | 131.72, `glmed1_r` (`MECHANISM_VASCULATURE.md`) | exact |
| knee joint force (calibration check, §below) | self-computed peak = 391.10 %BW @ t=0.51s | 391.10/391.11 %BW (`MECHANISM_STATIC_OPT.md`) | exact |

**Calibration check on perturbation (B)** (an over-determination, not a tautology: computed via a
**different code path** — `static_opt_knee.py`'s geometric muscle-crossing projection — than
`docs/MECHANISM_EMG_DRIVEN.md` used, `opensim.JointReaction`): at SO's own peak instant (t=0.51s,
the fair, non-edge-flagged reading), baseline self-computed knee contact force = 391.10 %BW,
perturbation-(B) = 568.11 %BW, a **+45.3%** change — matching the referenced **+48.4%** to within 3
points via an independent method and a different baseline (raw SO here vs. `naive_fwd_only` there).
Same direction, same order of magnitude: perturbation (B) is a faithful, calibrated proxy for the
real co-contraction content, not an arbitrary injected number.

## Per-layer sensitivity — the core result

Elasticity = (% change in downstream headline) / (% change in the raw SO input, 20.00% realized for
bias A). Ranked by |elasticity| under bias (A):

| rank | layer | headline | baseline | bias A (+20% activation): value / %Δ / **elasticity** | bias B (EMG co-contraction): value / %Δ | classification |
|---|---|---|---:|---|---|---|
| 1 | **muscle_fatigue** | 60-min sustained peak force-decline %, worst muscle | 7.964% (`glmed1_r`) | 31.96% / +301.3% / **15.06** | 50.06% / +528.6% | **MASSIVE AMPLIFIER** |
| 2 | **tendon_elastic_energy** | Σ Achilles-group own-peak stored energy | 11.180 J | 15.051 J / +34.6% / **1.73** | 26.325 J / +135.5% | **AMPLIFIER** |
| 3 | muscle_perfusion | peak settled flow, max/80 muscles | 131.72 ml/min/100g | 156.87 / +19.1% / **0.95** | 181.15 / +37.5% | ~neutral |
| 4 | metabolic_cost | Umberger2010 net power | 8.1807 W/kg | 9.0763 / +10.9% / **0.55** | 11.113 / +35.8% | dampened |
| 5 | cardiac_output (Q) | cardiac output | 20.087 L/min | 22.005 / +9.5% / **0.48** | / +31.3% | dampened |
| 5 | respiratory (VO2) | O2 uptake | 26.930 mL/kg/min | 29.502 / +9.5% / **0.48** | / +31.3% | dampened |
| 6 | thermoregulation (dT/dt) | whole-body heat-rate | 0.1296 °C/min | 0.1415 / +9.2% / **0.46** | / +30.1% | dampened (double) |
| 6 | thermoregulation (sweat) | required sweat rate | 735.7 g/h | 816.3 / +10.9% / **0.55** | / +35.8% | = net-power (identity) |
| 7 | cardiac_output (HR) | heart rate | 189.70 bpm | 202.7 / +6.9% / **0.34** | / +21.2% | most-dampened chain quantity |
| 8 | **validate_emg_timing** | 8 muscle-group verdicts + Jaccard | 0P/3PA/5F, J̄=0.5138 | **0 flips, J̄ unchanged to every digit** / **elasticity = 0 EXACTLY** | **2/8 flip** (`glmed` FAIL→PARTIAL, `soleus` PARTIAL→FAIL), J̄=0.5210 | ROBUST-to-A, exposed-to-B |
| 9 | **motor_unit_recruitment** | Fuglevand pool recruitment/force curve | — | **not a consumer at all** (audited) | **not a consumer at all** | **ZERO COUPLING** |

## Why each layer moves the way it does (derived from structure, then verified — not asserted)

**muscle_fatigue is the outlier, and it is a real, well-characterized nonlinearity, not a fluke of
one data point.** A dose-response sweep (4 magnitudes: +5/+10/+20/+30% realized activation bias, all
independently re-run) gives elasticity **15.84 / 17.15 / 15.06 / 13.64** — a stable ~14-17× band
that mildly *softens* at the largest bias, not a divergence. The worst-affected muscle (`glmed1_r`)
is identical across all four magnitudes. The mechanism is the one the original script's own STEP 5
already names: the fatigue model's state obeys `MR+MA+MF=1` (a simplex conservation law); *visible*
force decline requires a muscle's own peak task-level demand plus its accumulated fatigue fraction
to collide with that `1.0` ceiling. A muscle sitting close to, but under, that ceiling shows a
**disproportionate** jump in visible decline for a modest activation increase — the same
"geometric insight" mechanism `muscle_fatigue.py` itself uses to explain why `glmed1_r` shows large
visible decline while a higher-hidden-fatigue muscle does not. **Symmetric caveat**: the layer's OWN
external anchor (Frey Law & Avin 2010 match, 35/35 points pass) is measured from a *separate*
synthetic constant-task-level sweep that never touches this trial's SO activation at all — that
specific validated claim is untouched by the common-mode risk; only the real-gait 60-min-decline
headline (the number actually exposed here) is.

**tendon_elastic_energy amplifies because the tendon force-strain curve is convex.** Predicted
*before* measuring, from the curve's own geometry: near the low-strain "toe" region, `F ~ x²` gives
`E = ∫F dx ~ F^1.5`; in the stiffer, near-linear region, `F ~ x` gives `E ~ F²`. Either way, `dE/E =
n · dF/F` with `n ∈ [1.5, 2]`. Measured elasticity **1.73** sits exactly inside that predicted band
— a geometric prediction made from the curve shape, then confirmed on raw data, not the reverse.

**metabolic_cost, respiratory, cardiac_output(Q), thermoregulation all dampen, and by a common,
exactly-verifiable mechanism: a FIXED additive floor dilutes relative sensitivity.**
`gross = net + basal` (basal = 1.2 W/kg, a constant, unaffected by the SO bias) →
`elasticity_gross = net/(net+basal)`. Measured exactly: `8.1807/9.3807 = 0.8721`; `0.8721 × 10.947%
(net's own elasticity·input) = 9.547%` — matching the measured gross-COT change (9.5470%) to 4
significant figures. VO2 and Q are pure proportional conversions of gross power (elasticity
*exactly* 1.0 relative to gross, confirmed: VO2's %change = gross-COT's %change to 6 sig figs) — so
they inherit metabolic_cost's own dampening, not more.

**thermoregulation's own two headline numbers behave DIFFERENTLY — one layer, two exposures.**
Required sweat rate is driven by `(1−η)·(gross−rest) = (1−η)·net_power` — an algebraic identity
(`gross−rest` cancels exactly to `net`, by construction, since `gross` is *defined* as `net+rest`),
so sweat rate's elasticity is *identical* to net power's own elasticity (confirmed: 10.9474% vs
10.9474%, matching to 6 sig figs — not amplified, not dampened relative to net power). Whole-body
dT/dt instead uses the η-*blended* `H_prod = η·rest + (1−η)·gross = rest + (1−η)·net_power` — a
fixed term plus a throttled variable term — giving a **doubly**-dampened elasticity:
`(1−η)·(net/H_prod)·elasticity_net = 0.775 × (639.7/589.6) × 0.547 = 0.460`, matching the measured
0.4603 exactly. (An earlier pass through this analysis predicted amplification here by conflating
the bare exercise-increment with the η-blended `H_prod` that actually drives dT/dt — caught when the
numbers disagreed with the prediction, corrected, not smoothed over; the corrected mechanism is what
is reported above and matches to 4 significant figures.) **The lesson: "one elasticity per layer" is
itself an oversimplification a lazy pass would have missed** — two headline numbers from the same
script can have different exposure.

**cardiac_output's HR is the most-dampened chain quantity**, because stroke volume (`SV`) is itself
an increasing, nonzero-intercept affine function of the same VO2 that drives `Q` — `HR = Q/SV(Q)`
with `SV` rising alongside `Q` mathematically *guarantees* sub-unit elasticity of `HR` with respect
to `Q` (confirmed ratio 0.7176, so HR's elasticity relative to the raw SO bias, 0.343, is lower than
Q's own 0.477).

**validate_emg_timing is EXACTLY invariant to the uniform bias, and this is a provable geometric
fact, not an empirical coincidence.** Its ON/OFF classification uses `(activation − floor)/(peak −
floor)`, a min-max normalization; scaling every value by the same constant scales floor and peak by
that constant too, leaving the ratio bit-for-bit unchanged — **unless clipping at the 1.0 activation
ceiling breaks the symmetry** (confirmed: 0/12,640 entries clipped under bias A, hence the measured
mean Jaccard is unchanged to *every reported digit*, 0.5137863018973969 before and after). This same
normalization offers **no** protection against a *non-uniform*, phase-localized perturbation: bias
(B) — which raises specific muscles at specific times, not everything by the same factor — flips
2/8 verdicts (`glmed` FAIL→PARTIAL, `soleus` PARTIAL→FAIL). **This layer is robust to one entire
CLASS of common-mode bias (uniform miscalibration) and exposed to another (missing co-contraction
timing) — "robust" must specify which kind of SO error.**

**motor_unit_recruitment has ZERO coupling to the SO solve, full stop.** Its own docstring frames
the motor-unit pool as sitting conceptually "underneath" SO's activation output, but that is
motivating prose, not a data dependency — verified directly by source inspection (not assumed
either way): zero references to `ACT_STO`, `SO_DIR`, `static_optimization`, or any `.sto` file
anywhere in the module's *executable* code (the one docstring mention is prose, confirmed by
checking the string occurs only inside `__doc__`). Its own numbers (Fuglevand pool recruitment
order, the force-vs-excitation curve shape) are driven entirely by a synthetic excitation ramp. This
is not "robust" — it is **not a consumer of this SO solve at all**, despite reading like one of the
9 layers the trust ledger's caveat names. A real, disclosed, symmetric finding: the audit could just
as easily have found real coupling and did not.

## Systemic exposure quantification — the honest answer to "how wrong would the stack be"

There is no single multiplier. Two concrete, evidence-based scenarios:

- **If the SO solve is uniformly biased +20%** (a global miscalibration, e.g. a systematic scaling
  error in the objective or solver tolerance): net metabolic power moves +10.9%, tendon elastic
  energy +34.6%, muscle fatigue's 60-min decline moves +301% (a ~15× amplification), perfusion
  +19.1%, VO2/Q +9.5%, HR +6.9% — **and** the EMG-timing PASS/PARTIAL/FAIL verdicts and the
  motor-unit layer do not move AT ALL. A reader who only checked "SO passed its convergence gate"
  would have no way to guess this 30×-plus spread in downstream consequences from that gate alone.
- **If the SO solve is missing real co-contraction at the magnitude the one genuinely decorrelated
  measurement in this repo (real EMG) shows** (not a hypothetical — this is what
  `docs/MECHANISM_EMG_DRIVEN.md` already measured): muscle fatigue's headline would be understated by
  **>6×** (7.96%→50.06%), tendon elastic energy by **~2.4×** (11.18J→26.32J), metabolic power by
  **~36%**, and 2 of 8 EMG-timing verdicts would flip in the direction of appearing *worse*
  (`glmed`) or *better* (`soleus`) than the true SO-only number suggests.

**The exposure ranges over three-plus orders of magnitude (0 to >500%) depending entirely on which
downstream layer and which kind of SO bias** — governed by each layer's own nonlinear structure
(a conservation ceiling, a convex material law, a fixed additive floor, a scale-invariant
normalization), not by a single "propagation factor." The trust ledger's caveat was correct that
agreement among these layers is internal consistency, not external validation — this probe adds the
missing number: for the two most amplifying layers (`muscle_fatigue`, `tendon_elastic_energy`), that
internal consistency is worth very little, because a modest, plausible SO bias moves their headline
numbers by more than the number itself. For the two immune layers (`validate_emg_timing` under
uniform bias, `motor_unit_recruitment` always), it is worth correspondingly more.

## Honest gaps (pre-registered items first)

1. **This does not fix the common-mode risk.** Only a second, genuinely independent SO-equivalent
   solve, or driving more of the model from real EMG (`docs/MECHANISM_EMG_DRIVEN.md`, currently 7/80
   muscles), actually decorrelates these numbers. This probe measures exposure, nothing more.
2. **Two perturbation shapes tested, not the space of all possible SO biases.** A uniform
   multiplicative bias and one EMG-calibrated co-contraction offset are informative, plausible
   probes — not proof about, e.g., a bias correlated with gait phase, a systematic sign error in one
   muscle group, or a solver-tolerance-driven bias. The muscle_fatigue dose-response sweep (4
   magnitudes) shows the elasticity is a stable local property over 5-30% activation bias, not
   necessarily beyond that range.
3. **Single trial, single subject** (subject2/`walking1`) — the same scope limit as every other cert
   in this family. No claim of generality across subjects, trials, or gait speeds.
4. **Perturbation (B) reuses only the 7/80 EMG-covered muscles' delta**, zero elsewhere — a
   disclosed floor on how much of the true co-contraction gap this specific perturbation can
   represent (the other 73 muscles may have their own unmeasured co-contraction content this probe
   cannot inject, since no real EMG covers them).
5. **The "elasticity" framing assumes local, roughly stable sensitivity.** Confirmed reasonable for
   muscle_fatigue via the 4-point sweep; NOT independently re-verified for the other 8 layers at a
   second magnitude (a disclosed, not fatal, gap — the chain layers' elasticities are exact algebraic
   consequences of closed-form formulas, verified to 4+ significant figures, so they are not at risk
   of being an artifact of the one magnitude tested; `muscle_perfusion` and `tendon_elastic_energy`
   were not swept beyond the two variants reported here).
6. **`validate_emg_timing`'s exact scale-invariance is a property of the min-max-normalization
   metric, not evidence the underlying activation trace is unbiased** — a badly biased SO solve
   could pass this layer's own gates unchanged while being wrong in absolute terms; robustness here
   is about verdict STABILITY under one class of common-mode error, not correctness.
7. **The knee-force calibration check (Sec. "Method" above) confirms direction + order of magnitude
   (+45.3% vs. +48.4%), not exact numerical agreement** — a real, disclosed ~3-point gap, expected
   given the different baseline (raw SO here, `naive_fwd_only` in `docs/MECHANISM_EMG_DRIVEN.md`) and
   different code path (geometric projection here, `opensim.JointReaction` there).

## Files

- `scripts/msk/common_mode_probe.py` — the full, self-contained, re-runnable probe (imports
  `metabolic_cost.py`, `muscle_fatigue.py`, `tendon_elastic_energy.py`, `validate_emg_timing.py`,
  `muscle_perfusion.py`, `thermoregulation.py`, `respiratory.py`, `cardiac_output.py`,
  `motor_unit_recruitment.py`, `static_opt_knee.py`, `validate_joint_force.py` — every physics
  computation reused from its own validated module, nothing reimplemented). ~150s wall time (9
  layers × 3 variants + the calibration check + the 4-point dose-response sweep).
- `data/msk_smoketest/subject2_walking1/common_mode_probe/common_mode_probe_results.json` — every
  number in this doc, machine-written (raw per-variant, per-layer results + the sensitivity table +
  the dose-response sweep).
- Real bug found and fixed while building this probe (disclosed, not hidden): the first version of
  `probe_tendon_energy` cached `TendonForceLengthCurve`/`Muscle` SWIG objects in a dict returned from
  a function whose local `model` variable then went out of scope — a use-after-free that
  intermittently raised a caught OpenSim C++ exception (`Property<T>::getValue()...`) and, in the
  full pipeline, an outright segfault (classic undefined-behavior symptoms of a dangling reference,
  not a math bug). Fixed by making `probe_tendon_energy` build its own model fresh every call
  (cheap: ~0.05s), never letting curve objects outlive their parent model's Python scope.
- Reused, not re-solved: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_{activation,force}.sto`; reused, not reprocessed:
  `data/msk_smoketest/subject2_walking1/emg_driven/_diag_cache.npz` (`a_hybrid`, `a_naive_only`,
  `f_hybrid`, `f_naive_only`, `muscle_names` — order-verified live to match the SO columns exactly
  before use, not assumed).
- Isolation (`COORDINATOR.md` §1): bodytwin only; every input read in place, never modified; all new
  outputs under `data/msk_smoketest/subject2_walking1/common_mode_probe/`; no git commit, no git
  push.
