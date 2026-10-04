# MECHANISM SUBJECT-SPECIFIC PATH — verifying and prototyping fidelity Ceiling #1 (2026-07-21)

**Purpose.** The fidelity architecture audit (`docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`) names
"ONE real subject-specific 3D imaging dataset" as the highest-ceiling-clearing migration (§5 item 2),
and cites a striking, load-bearing claim: muscle strength (`max_isometric_force`, Fmax) is
**bit-identical** between the generic template and subject2's scaled model — i.e. OpenSim's scaling
pipeline *structurally* never touches strength. This doc (1) independently re-verifies that claim,
(2) maps exactly what OpenSim's `Scale` tool does and does not touch, grounded in the actual C++
source (not recalled), and (3) prototypes ONE concrete improvement — a segment-size-derived PCSA
proxy for Fmax — and measures its effect on the knee/hip contact-force cert already in this repo.
Executed via `scripts/msk/subject_specific_scaling.py` (real run, exit 0; console + full JSON at
`data/msk_smoketest/subject2_walking1/pcsa_fmax_proxy/subject_specific_scaling_results.json`).

## 0. Executive verdict

**Claim 1 (Fmax bit-identical): CONFIRMED**, independently, via two decorrelated extraction methods
(0/480 mismatches) AND grounded in the actual OpenSim C++ scaling mechanism — it is a **structural**
property of the tool, not a coincidence of this one export. **Claim 2 (map what Scale touches):**
geometry (segment length, muscle/ligament path-derived lengths, joint frames) is rescaled; gross
body mass is rescaled only globally, and — a new, unplanned finding — subject2's own pipeline used
`preserve_mass_distribution=true`, which measurably makes body MASS carry **zero** per-segment shape
signal for this subject (uniform 1.038003 ratio across all 22 bodies); the model's separate,
per-body mesh geometry *does* vary genuinely per segment (1.20–1.99× volume). **Claim 3 (prototype):**
built two Fmax proxies and reran the existing validated Static-Optimization + dual-method
contact-force pipeline on both. Result, after catching and fixing a **stale-cache bug in my own
first pass** (below): the literal "mass-informed" proxy is measurably near-vacuous for this subject
(an honest, disclosed finding, not swept under the rug) and, coincidentally, moves the knee estimate
closer to the OrthoLoad anchor; the more geometrically-defensible proxy (real per-segment volume
scale factors) moves it further away. **Neither proxy is a validated fidelity improvement** — the
exercise's real value is a working, mechanistically-understood pipeline for testing Fmax proxies at
all, and a clear diagnosis of *why* the knee is far more sensitive to this choice than the hip.

---

## 1. VERIFY: is Fmax really bit-identical, generic vs scaled?

**Models compared** (same pair the muscle audit used): GENERIC =
`.../subject2/OpenSimData/Mocap/Model/LaiArnoldModified2017_poly_withArms_weldHand_generic.osim`;
SCALED = `..._scaled.osim` (both on the read-only `mechanism_data` NTFS drive). 80 muscles each
(`Millard2012EquilibriumMuscle`), same name set.

**Method (two independent extraction paths, cross-checked against each other):**
1. Raw XML parse (`xml.etree.ElementTree`, no OpenSim API involved).
2. OpenSim 4.6 Python API (`opensim.Model(...).getMuscles()`, `getMaxIsometricForce()` etc.).

| check | result |
|---|---|
| XML-parse vs API-parse agreement (Fmax/Lopt/Lts, 80 muscles × 2 models = 480 comparisons) | **0 mismatches** — methods agree |
| `max_isometric_force`, generic vs scaled, exact (tol=0) | **0/80 differ → CONFIRMED bit-identical** |
| `pennation_angle_at_optimal`, generic vs scaled | 0/80 differ |
| `max_contraction_velocity`, generic vs scaled | 0/80 differ (not reliably serialized in the raw XML — 0/80 muscles have an explicit tag, confirmed by grep; compared via the API's default-resolved value instead, disclosed methodological wrinkle) |
| `optimal_fiber_length`, generic vs scaled | 80/80 muscles change, ratio ∈ **[1.0635, 1.2539]**× (pre-registered sanity band [0.5, 2.0]×: 0/80 outside) |
| `tendon_slack_length` ratio vs `optimal_fiber_length` ratio, same muscle | **identical for 80/80 muscles** (to float precision) |

The last row is the empirical fingerprint of a single shared scale factor driving both quantities —
explained, not just observed, in §2.

**Verdict: CONFIRMED.** This independently reproduces `docs/MECHANISM_MUSCLE_AUDIT.md` §3's claim with
a different extraction method (raw XML vs that audit's OpenSim-API-only approach) and adds a
mechanism-level anchor (§2) the original audit did not have.

---

## 2. MAP: what does `OpenSim::Model::scale()` actually touch?

Traced directly in the C++ source (`opensim_jam_build/opensim-core-jam`, an independently-built
OpenSim-JAM fork of OpenSim 4.x — a different artifact from the compiled pip package `opensim==4.6`
that ran the checks above, i.e. two decorrelated copies of the "same" library, not one).
`Model::scale()` (`Simulation/Model/Model.cpp:1590`) runs three passes over every `ModelComponent`:
`preScale` → `scale` (`extendScale`) → `postScale` (`extendPostScale`), then a separate body-mass step.

| component / property | touched? | mechanism (file:line) |
|---|---|---|
| `Body.mass_center` | **YES, always** | `Body::extendScale` (`SimbodyEngine/Body.cpp:223`) — repositioned by the body's own 3-axis scale factors |
| `Body.mass` (per-body) | **CONDITIONAL** | `Body::scaleInertialProperties(scaleFactors, scaleMass=!preserve_mass_distribution)` — subject2's `setup_scale.xml` has `<preserve_mass_distribution>true</preserve_mass_distribution>` (grepped directly), so this per-body step is **skipped** |
| `Body.mass` (ALL bodies, global) | **YES** | `Model::scale()`'s `finalMass` step re-multiplies every body's mass by ONE constant = `target_mass / current_total_mass`, regardless of the row above. **Measured live: exactly 1.038003 for all 22 bodies** (pelvis through hand_l) — mass carries **zero** per-segment shape signal for this subject |
| `Body` mesh `scale_factors` (visual/geometric, separate XML field) | **YES, genuinely per-segment** | Not part of the physics scaling pass above — a per-body recorded 3-axis geometric ratio. Measured range **1.20–1.99× volume** across segments (femur_r 1.72×, tibia_r 1.99×, torso 1.20×, pelvis 1.39× — anisotropic) — sane for subject2 (**1.96 m tall**, per `sessionMetadata`, vs a shorter nominal template) |
| `Body.inertia` | **YES, always** | `scaleInertialProperties` — length² term always applied, ×mass-factor too if mass is also being scaled |
| `Muscle.optimal_fiber_length` | **YES** | `Millard2012EquilibriumMuscle::extendPostScale` (`Actuators/Millard2012EquilibriumMuscle.cpp:810-825`): `scaleFactor = path.getLength(s) / path.getPreScaleLength(s)`; `optimal_fiber_length *= scaleFactor` |
| `Muscle.tendon_slack_length` | **YES** | same function, same call: `tendon_slack_length *= scaleFactor` — **the same variable**, which is exactly why Lopt-ratio == Lts-ratio for 80/80 muscles in §1 |
| `Muscle.max_isometric_force` | **NO — never** | not referenced anywhere in `extendPostScale`; the base `Muscle` class (`Simulation/Model/Muscle.cpp`/`.h`) defines **no** `extendScale`/`extendPreScale`/`extendPostScale` of its own — structurally absent from the class hierarchy, not merely unset by this one pipeline run |
| `Muscle.pennation_angle_at_optimal`, `.max_contraction_velocity` | **NO — never** | same absence |
| `Ligament.resting_length` / `Blankevoort1991Ligament.slack_length` | **YES** | `Ligament::extendPostScale` / `Blankevoort1991Ligament::extendPostScale` (`Simulation/Model/Blankevoort1991Ligament.cpp:130-143`) — **identical pattern**: `*= scaleFactor` |
| `Blankevoort1991Ligament.linear_stiffness` (the ligament's Fmax-analogue) | **NO — never** | not referenced — an independent confirmation of the same structural-absence pattern, on a *different* tissue class |
| `CustomJoint`/`EllipsoidJoint`/`ScapulothoracicJoint`/`ConstantCurvatureJoint` frames | **YES** | `extendScale` repositions parent/child joint frames to track their (now-scaled) bodies |

**Reading:** OpenSim's scale mechanism systematically rescales **geometry** (segment dimensions,
path-derived lengths, joint frames) and, separately, **gross body mass** (globally, not necessarily
per-segment — gated by `preserve_mass_distribution`) — but **never** rescales tissue **strength**
(muscle Fmax, ligament stiffness) on any tissue class anywhere in the source tree. This generalizes
the audit's claim from "true for subject2" to "true for any OpenSim `ScaleTool` run on this model
family" — a structural fact, not a one-off.

---

## 3. PROTOTYPE: a PCSA proxy for Fmax

**Geometric derivation.** PCSA = Volume / OptimalFiberLength (standard identity; Fmax = specific
tension × PCSA). `OptimalFiberLength`'s own scaled/generic ratio is already known exactly (§2's
`scaleFactor`) — only a **volume** proxy is needed, since no muscle-tissue imaging exists.

- **Proxy-A ("mass-informed", built exactly as asked):** volume proxy = the muscle's origin-body
  MASS ratio (scaled/generic). **Measured to be the constant 1.038003 for every body** (§2) — this
  proxy is honestly near-**vacuous** for this subject/pipeline: `Fmax_new ≈ Fmax_gen × 1.038 /
  Lopt_ratio`, i.e. it just anti-tracks each muscle's own (already-known) Lopt ratio through a fixed
  constant, contributing no real per-segment differential information.
- **Proxy-B ("geometric segment-volume", the forced adversary on Proxy-A's degeneracy):** volume
  proxy = the origin body's own recorded mesh-geometry volume ratio (§2's `Sx·Sy·Sz`), which **does**
  vary genuinely per segment.
- `Fmax_new = Fmax_generic × (Volume_ratio_origin_body / Lopt_ratio_muscle)`, applied to a **clone**
  of the scaled model (only Fmax overwritten — everything else, including the already-scaled Lopt/Lts/
  paths, stays exactly as the audited pipeline produced it). New files: `data/msk_models/
  subject2_scaled_pcsa_fmax_proxyA_mass.osim`, `..._proxyB_geom.osim`.

| muscle | origin | Proxy-A ratio | Fmax gen→A | Proxy-B ratio | Fmax gen→B |
|---|---|---:|---|---:|---|
| soleus_r | tibia_r | 0.828 | 6194.8 → 5129.0 N | 1.590 | 6194.8 → 9849.1 N |
| gasmed_r | femur_r | 0.828 | 3115.5 → 2579.2 N | 1.374 | 3115.5 → 4281.3 N |
| glmax1_r | pelvis | 0.940 | 983.8 → 924.6 N | 1.259 | 983.8 → 1238.1 N |
| recfem_r | pelvis | 0.870 | 2191.7 → 1906.6 N | 1.165 | 2191.7 → 2553.2 N |
| psoas_r | pelvis | 0.952 | 1426.8 → 1358.5 N | 1.275 | 1426.8 → 1819.2 N |
| tibant_r | tibia_r | 0.836 | 1227.5 → 1026.1 N | 1.605 | 1227.5 → 1970.4 N |

Full-model distribution (n=80, pre-registered plausibility band [0.3, 4.0]×, 0/80 outside for
either proxy): **Proxy-A ratio ∈ [0.828, 0.976], mean 0.879** (Fmax **decreases** for every single
muscle — mechanically inevitable once mass-ratio is a constant < every Lopt-ratio); **Proxy-B ratio
∈ [1.109, 1.648], mean 1.320** (Fmax increases 11–65%, physiologically the expected *direction* for
a 1.96 m subject vs a shorter template).

**Sensitivity check** (origin-body vs mean-of-all-crossed-bodies convention, the other reasonable
way to assign a segment to a multi-body-spanning muscle): Proxy-A max relative difference = **0.00%**
(mechanically trivial — the mass ratio is a global constant, so the assignment convention cannot
matter); Proxy-B max relative difference = **22.45%** (mean 12.09%) across the 80/80 muscles that
cross >1 body — a genuine, disclosed sensitivity of the geometric proxy to this modeling choice.

---

## 4. Effect on knee/hip contact force

### 4.1 A stale-cache bug caught before trusting the comparison

My first full run compared the two proxies against a **cached** baseline number
(`data/msk_smoketest/subject2_walking1/static_optimization/static_opt_knee_results.json`,
written 2026-07-21 13:06): knee self-computed = 233.20 %BW. But `scripts/msk/static_opt_knee.py`'s
*current* source already contains a sign-bug fix to `knee_crossing_muscles_and_forces` (its own
inline comment documents it — see `docs/MECHANISM_SIGN_BUG_REMEDIATION.md`), and the fidelity audit's
own text independently states the **corrected** value is 391.10 %BW, converging with
`opensim.JointReaction` to <0.1%. The cached JSON was never regenerated after the fix landed —
233.20 is the *pre-fix* number. **Forced check:** recomputed the self-cross-check fresh, from the
same underlying (unaffected — Fmax-independent) SO output, using the current code:

| | cached JSON (13:06) | fresh recompute (current code) | official JointReaction |
|---|---:|---:|---:|
| knee_r self-computed | 233.20 %BW | **391.10 %BW** | 391.11 %BW |

Fresh matches JointReaction to <0.01% — exactly the convergence the audit's text describes; the
cache was stale. (Hip was **not** stale: 386.77 cached == 386.77 fresh, because
`validate_hip_force.py`'s cross-check was written after the fix was already in place.) Fixed at the
source: `subject_specific_scaling.py`'s `compute_baseline_fresh()` now **always** recomputes through
the identical, current, live-executed code path shared with both proxies — it never reads a cached
`self_computed` number again. All numbers below use the corrected baseline.

### 4.2 Corrected comparison (Static Optimization rerun fresh for both proxies; SO convergence PASS for both; self-computed and official-JointReaction agree to <0.1% in every row, same discipline as the existing pipeline)

| model | knee_r %BW | knee vs OrthoLoad (258.22) | hip_r %BW | hip vs OrthoLoad (273.93) |
|---|---:|---:|---:|---:|
| baseline (generic Fmax) | 391.10 | 1.515× | 386.77 | 1.412× |
| Proxy-A (mass-informed) | 348.08 | 1.348× | 385.93 | 1.409× |
| Proxy-B (geometric volume) | 430.79 | 1.668× | 392.51 | 1.433× |
| OrthoLoad anchor (median) | 258.22 | 1.000 | 273.93 | 1.000 |

Baseline (zero strength-scaling — the audited status quo) **already overshoots** both anchors by
40–52%. Proxy-A moves the knee *closer* (1.515×→1.348×); Proxy-B moves it *further away*
(1.515×→1.668×). Both proxies leave the hip almost unchanged (1.41×→1.41–1.43×).

### 4.3 Forcing the mechanism (Orient, not just report the aggregate)

Naive expectation was ambiguous going in; the measured direction (Proxy-B *increasing* Fmax makes
the knee estimate *worse*, not better) is not self-evidently a bug, so it was forced rather than
reported flat. Per-muscle SO force+activation at the knee's peak instant (t=0.51s, identical across
all three runs) shows why:

- **The knee-crossing force is dominated by ONE muscle**: `gasmed_r` (medial gastrocnemius) supplies
  **45.7–47.1% of the total 13-muscle knee-crossing force sum**, stably across 5 nearby frames
  (t=0.47–0.55s) — not a one-frame artifact. Proxy-B raises `gasmed_r`'s own Fmax by 37% (origin body
  femur_r, volume ratio 1.72×).
- Static Optimization minimizes **Σactivation²**, not force — it does not conserve total force when
  Fmax changes. At the peak instant, `gasmed_r`'s activation **fell** under Proxy-B (0.441→0.384, as
  naively expected for a now-stronger muscle) but its **force still rose** (1140 N→1364 N, +20%),
  because `force = activation × Fmax × f(l,v)`: a 37% Fmax increase outweighs a 13% activation drop.
  Since routing force through a now-individually-stronger actuator is *activation-cheaper* per unit
  force, the quadratic-cost optimum does not "give back" the full capacity increase as reduced force.
- **The hip is comparatively robust to exactly the same perturbation** because its crossing-force sum
  is spread over 25 muscles with **no single muscle exceeding 24%** of the total (glmed1_r 23.2–23.8%,
  iliacus_r/recfem_r/psoas_r 13–20% each) — an individual muscle's Fmax change gets diluted across a
  much more even distribution, rather than dominating the aggregate the way `gasmed_r` does at the knee.

**Geometric reading:** a joint's contact-force estimate's sensitivity to any single muscle's Fmax
scales with how *concentrated* the load-bearing is across that joint's crossing-muscle set at the
loaded instant — an effective-number-of-contributors argument, not a per-joint coincidence.

### 4.4 Honest verdict on the prototype

**Neither proxy is a validated fidelity improvement.** Proxy-A's apparent knee improvement is
coincidental: it comes from a proxy already shown (§3) to carry no real per-segment signal for this
subject, which happens to point the ratio-correcting direction for the knee's specific overshoot —
not because reduced-strength is subject-specifically correct. Proxy-B, the more geometrically
defensible construction (real, subject-measured segment size variation, physiologically the right
*direction* for a very tall subject), makes the already-overshooting knee number worse. The
honest conclusion is not "which proxy wins" — it is that **the joint-force cert's response to a
strength perturbation is now measured, mechanistically explained, and shown to be joint-specific**
(knee: high-sensitivity/single-muscle-dominated; hip: low-sensitivity/distributed) — the necessary
groundwork before any subject-specific strength signal (real or proxied) can be trusted at either
joint.

---

## 5. Honest ceiling — the ceiling this does NOT clear

Both proxies are **bone-segment** geometry/mass signals; **neither ever measures muscle tissue**.
Concretely, not closed by this work:
- **No imaging of any kind.** Individual hypertrophy/atrophy, sarcopenia, or asymmetric training
  effects (e.g. a dominant-leg athlete) are invisible to a bone/segment-geometry proxy — bone
  geometry does not change with muscle-specific gain or loss.
- **Single origin-body assignment is a simplification**, crudest for biarticular muscles (measured
  up to 22.45% sensitivity to the origin-vs-mean-crossed-body convention, §3) — a real, disclosed
  source of the proxy's own uncertainty.
- **Specific tension is assumed constant** (population-literature value) — never subject-measured;
  it cancels out of the *ratio* used here but is still an unverified absolute assumption.
- **Pennation angle** is untouched by either proxy (as it is by the base scaling pipeline, §1/§2).

**Concrete path to real subject-specificity** (ranked by cost, cheapest first):
1. **Published anthropometry-based muscle-volume regressions** (e.g. Handsfield et al. 2014-style
   equations fit to actual cadaver/MRI individual-muscle volumes as a function of height × mass) —
   still not imaging, but replaces this doc's crude single-segment-scale proxy with a proxy fit to
   real per-muscle (not per-bone) data. No new acquisition needed, a software-only upgrade.
2. **Ultrasound fascicle-length + pennation measurement**, per muscle, this subject — directly
   replaces the assumed-generic Lopt/pennation with a subject-measured value (cheap, common in sports
   science, no MRI needed).
3. **MRI/CT volumetric segmentation of individual muscles** → true per-muscle volume → PCSA =
   volume / (measured) fiber length → Fmax = specific_tension × PCSA. This is the audit's own
   Ceiling #1 migration — the only step that replaces a *proxy* with a *measurement*.
4. Dynamometry-based isokinetic torque-velocity curves (this subject's own) as the discriminating
   validation anchor for whichever Fmax source is used — never measured for subject2, noted as a
   falsifier already in `MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`'s Muscles row and still true here.

This work is scoped to the **muscle-strength** slice of Ceiling #1 only. Lopt/Lts are already
touched by the existing pipeline (via geometry, not imaging — a real but different signal); the
other three ceilings (Hill→cross-bridge, rigid-body→deformable, monocular→medical-grade acquisition)
are untouched and remain exactly as ranked in the architecture audit.

---

## Falsifiers (pre-registered, machine-checked, all in the script's own PASS/FAIL output)

- Fmax bit-identical claim: falsified if any of 80 muscles differ by >0.0 N between generic/scaled,
  or if the XML-parse and API-parse extraction methods disagree on any of 480 comparisons. **0/80,
  0/480 — claim survives.**
- Mechanism map: falsified if `Millard2012EquilibriumMuscle::extendPostScale` or the base `Muscle`
  class referenced `max_isometric_force` anywhere. **Grepped directly — 0 references.**
- PCSA-proxy plausibility: falsified if any of 80 muscles' new Fmax ratio fell outside the
  pre-registered [0.3, 4.0]× band. **0/80 for either proxy.**
- SO convergence: falsified if `convergence_pass=False` for either proxy run, or if self-computed
  and official JointReaction disagreed by >1% at either joint. **Both PASS; max disagreement 0.07%
  (Proxy-A hip: 385.93 vs 386.21).**
- Stale-cache check: falsified (i.e. "no bug") if fresh recompute matched the cached JSON. **It did
  not (233.20 vs 391.10) — the bug is real, disclosed, and fixed at the source in this script.**

## Evidence index

Script: `scripts/msk/subject_specific_scaling.py` (real run, exit 0). Full JSON evidence:
`data/msk_smoketest/subject2_walking1/pcsa_fmax_proxy/subject_specific_scaling_results.json`. New
models: `data/msk_models/subject2_scaled_pcsa_fmax_proxyA_mass.osim`,
`..._proxyB_geom.osim` (Fmax overwritten in-XML, spot-verified against console output post-write).
New SO/JR runs: `data/msk_smoketest/subject2_walking1/pcsa_fmax_proxy/{proxyA_mass,proxyB_geom}/{so,jr}/`.
Baseline (generic-Fmax) SO/JR cache: `data/msk_smoketest/subject2_walking1/static_optimization/`
(read, never overwritten — mtimes verified unchanged, still 2026-07-21 13:06). Prior docs relied on:
`docs/MECHANISM_MUSCLE_AUDIT.md`, `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`,
`docs/MECHANISM_JOINT_FORCE_VALIDATION.md`, `docs/MECHANISM_STATIC_OPT.md`,
`docs/MECHANISM_SIGN_BUG_REMEDIATION.md`. External C++ source read (read-only, isolation-respected):
`/media/anton/8838D60F38D5FBDE/mechanism_data/opensim_jam_build/opensim-core-jam/OpenSim/{Simulation/Model/Model.cpp,Simulation/SimbodyEngine/Body.cpp,Simulation/Model/GeometryPath.cpp,Actuators/Millard2012EquilibriumMuscle.cpp,Simulation/Model/Ligament.cpp,Simulation/Model/Blankevoort1991Ligament.cpp}`.

Isolation respected throughout: bodytwin only; the two external NTFS drives (`mechanism_data`,
OrthoLoad) and the opensim-core-jam source tree were read-only; all new files are under
`scripts/msk/`, `docs/`, `data/msk_models/`, `data/msk_smoketest/` inside this repo. No git commit,
no git push.
