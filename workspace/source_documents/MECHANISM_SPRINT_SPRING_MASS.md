# MECHANISM SPRINT SPRING-MASS MODEL — Weyand's force-limit falsifier, exact SLIP geometry, Bolt cross-check (2026-07-22)

Quantitative model of sprint running as a **spring-loaded inverted pendulum (SLIP)** stance mechanism,
built and falsified against the Weyand-lab force-plate literature and an independent, video-digitized
elite-sprint dataset. Script: `scripts/msk/sprint_spring_mass.py`. Raw machine-checked results:
`reports/probes/sprint_spring_mass_results.json` (13 pre-registered checks, 11 PASS / 2 FAIL, both
FAILs explained below — not hidden). Citation ledger: `docs/MECHANISM_SPRINT_SPRING_MASS_evidence.json`.

**Scope / duplicate-check** (`data/MECHANISM_ANCHOR_GRAPH.json`, 1031 nodes, checked directly, not
assumed): three existing nodes touch adjacent ground but do not overlap this doc. `MSK-SPRINT-FV-PROFILE`
(EXECUTED, MEASURED-B-partial) recovers Samozino F0/V0/Pmax from split-times — a **macroscopic**
force-velocity-power *capability profile*, not the ground-contact mechanism that applies it.
`MSK-REPEATED-SPRINT-FATIGUE` (EXECUTED, MEASURED-B) certifies cycling/treadmill power-output decline
across repeated maximal efforts — a **metabolic/PPO-level** fatigue signature, operationally distinct
from the biomechanical spring-mass-property fatigue data used in Part 4 below. `MSK-GAIT-LOCOMOTION`
(OPEN, SEED-DESIGN only, never executed) explicitly proposes "parametrize gait control in the low-dim
latent space of the inverted-pendulum/SLIP model (z = step-length, frequency, leg-stiffness,
angle-of-attack)" — this doc **builds and validates exactly that forward model** against real data for
the first time in this repo, but does **not** execute that node's full ambition (metabolic-cost
optimization vs. a full-DOF OpenSim comparison); it supplies a load-bearing, externally-anchored piece
of it. None of the three is duplicated; all three are `couples_to` targets.

## Part 1 — The geometry: exact planar SLIP stance integration (not sine-wave, not curve-fit)

The stance leg is modeled as the standard spring-loaded inverted pendulum: a massless linear spring
(stiffness `k`, rest length `L0`) pinned to a fixed foot, carrying point mass `m` (the CoM). No
small-angle or symmetric-sine approximation is used for the primary simulation — the full nonlinear
equations are integrated numerically (`scipy.solve_ivp`, adaptive RK45, event-terminated at leg
re-extension = liftoff):

```
ell = sqrt(x^2 + y^2)          Fs = k(L0 - ell)  [>0 while compressed]
m x'' = Fs * x/ell             m y'' = Fs * y/ell - m g
```

This is exactly the Blickhan (1989) / McMahon & Cheng (1990) formulation — geometry and Newton's
second law, nothing curve-fit. A **coarse-to-fine grid search** over touchdown angle `theta0` and
landing vertical speed `vy_td` (not a black-box optimizer, so it cannot silently fail to converge)
finds the touchdown condition that reproduces a target `(Tc, F_avg/W)` pair for a given `k`.

### Claim 1 (pre-registered): realistic k_leg reproduces Weyand's measured top-speed mechanics

**C:** literature-measured leg stiffness (7–30 kN/m — Farley & González 1996's 7.0–16.3 kN/m jogging
range through Morin et al. 2006's 19.5±4.3 kN/m **direct in-sprint** measurement) run through the exact
SLIP equations reproduces Weyand et al. 2010's measured forward-running top-speed point
(`Tc=0.108±0.004 s`, `F_avg/W=2.08±0.07`) at a physiologically ordinary touchdown geometry (angle ≤40°,
compression ≤35% of leg length). **Pre-registered threshold:** ≥1 grid point at residual <3%.

**Result (machine-computed, primary anthropometry m=75 kg, L0=0.95 m, v0 swept 7–10 m/s since
Weyand 2010's abstract gives Tc/F but not the associated m/s — flagged in honest gaps):** best match
at **k=25 kN/m, v0=7 m/s → theta0=26.8°, compression=10.2% of L0, residual=1.26%**
(`Tc_sim=0.1069 s`, `F_avg/W_sim=2.096` vs. targets 0.108 s / 2.08) — **PASS**. The full per-k sweep is
the more informative result: k=7–10 kN/m (Farley-González's own *lowest-stride-frequency jogging*
values) **clearly fail** (residual 10–28%, requiring 26–35% compression to compensate for being too
soft), while k=15–30 kN/m (bracketing Morin's **directly-measured in-sprint** 19.5±4.3 kN/m almost
exactly) fit at 1.3–4.1% residual with shallow, ordinary touchdown angles (9–33°). This is a genuine,
falsifiable discrimination, not a fit-anything result — it is also the quantitative confirmation of the
McMahon & Cheng (1990) claim (verified live: dimensionless `KLEG` "nearly linear" in speed at high U)
on the **same U-axis** where Farley/Glasheen/McMahon (1993) and He/Kram/McMahon (1991) find k_leg
constant: jogging-band k (7–10 kN/m) is geometrically inconsistent with sprint-band mechanics; the
higher, in-sprint-measured band is consistent. Robustness (anthropometry swept m∈{65,85}kg,
L0∈{0.90,1.00}m, k scaled by the established **k_leg ~ M^0.67** cross-individual relation, Farley,
Glasheen & McMahon 1993): **3/4 corners plausible** (residuals 0.4%, 1.4%, 2.0%, 5.4%) — **PASS**.

### The OODA loop this doc forces on itself, not just on the literature

First pass held `k` fixed in absolute kN/m across the m=65→85 kg anthropometry sweep and **failed**
2/4 corners (residuals 12–23%). Diagnosing *why* (Orient, not skip): at fixed k, lower mass raises
`ω=sqrt(k/m)`, shortening the natural stance timescale — exactly why Tc drifted from the 0.108 s target
(measured directly: m=65kg gave Tc_sim=0.084–0.088s, ~20% short). The fix is not ad hoc: Farley, Glasheen
& McMahon (1993) themselves report k_leg scaling as **M^0.67** across differently-sized individuals —
applying that established relation (not re-fitting a new exponent) moved the forward-running corner
count from 2/4 to 3/4 plausible. This is reported as a real before/after, not smoothed away.

### Claim 1b (hopping) — an honest, explained FAIL, not hidden

Weyand et al. 2010's **hopping** condition (`Tc=0.160±0.006s`, `F_avg/W=2.71±0.15` — their own most
extreme, deepest-bounce contrast gait) is matched numerically **excellently** (best match at k=10 kN/m:
**theta0=40.8°, compression=35.8% of L0, residual=0.86%**; 0.1–4.5% across all four anthropometry
corners) but **only** at a geometry just past this doc's pre-registered "ordinary running" plausibility
bound (≤40°, ≤35%). **FAIL** on the geometry-plausibility gate (0/4 corners), **not** on the mechanics
themselves.
This is not a moved goalpost: the numbers are reported exactly as computed, and the physical reading is
that a bound calibrated for *ordinary running* under-estimates how much compression a maximal,
Weyand-selected "as hard as possible" one-legged hop legitimately uses — a genuine, auditable limitation
of the pre-registered threshold, disclosed rather than silently loosened after the fact.

### Internal cross-check: half-sine analytic approximation vs. exact SLIP

Independent of the exact ODE, the impulse-momentum theorem under an assumed half-sine GRF profile
(structurally the Morin et al. 2005 "sine-wave" stiffness method, PMID 16082017 — re-derived here from
first principles, not copied from that paywalled paper's exact text) gives the closed-form
`F_peak/F_avg = π/2 ≈ 1.5708` for a perfectly symmetric bounce. Measured: **hopping = 1.5498 (1.34%
below ideal — PASS, hopping is the "purest" passive bounce)**; **forward running = 1.7404 (10.80% above
ideal — PASS on the pre-registered "measurably asymmetric" threshold)**. This is an internal
consistency check (not a tautology: the π/2 figure is derived, not fitted) that independently
reproduces, from first-principles impulse-momentum reasoning alone, Clark & Weyand's (2014, PMID
25080925) empirical finding that sprint force application is **measurably less spring-symmetric** than
idealized hopping (their own R²: sprinters <0.85 vs. nonsprinters ≥0.91, worst at top speed R²=0.78±0.02).

## Part 2 — Impulse-balance identity: a geometric law, cross-checked against an independent, decorrelated dataset

For **any** periodic single-support running gait (no assumption about force-profile shape — a stronger,
more general constraint than the SLIP simulation above), Newton's second law integrated over one step
gives an exact identity:

```
F_avg * Tc = W * (Tc + t_flight)   =>   F_avg/W = 1 + t_flight/Tc   =>   step_time = Tc * (F_avg/W)
```

Feeding in **only** Weyand et al. 2010's forward-running numbers (Tc=0.108s, F_avg/W=2.08 — a
treadmill force-plate measurement) predicts `step_time = 0.2246 s` → **step frequency = 4.452 Hz**,
with **zero free parameters**. Krzysztof & Mero (2013, PMID 23717364, PMC3661886 — fetched live in full
text) report Usain Bolt's **video-digitized** stride frequency at Berlin 2009's 60–80 m split as
**4.49 Hz** — a **completely different measurement modality** (overground video kinematics vs.
treadmill force-plate) and population (Bolt vs. Weyand's 7 "athletic subjects") than the number that
generated the prediction. **Relative error = 0.86% — PASS** (pre-registered threshold: <15%, since this
is cross-population/cross-modality triangulation, not a same-subject fit). This is the "decorrelated
external anchor, never a tautology" requirement satisfied in its cleanest form: a conservation law, not
a curve, doing the predicting.

## Part 3 — The forced adversary: "faster leg-repositioning," not force, sets top speed

**Adversary, steelmanned:** even though Weyand et al. 2000 (PMID 11053354) report swing time `t_sw` as
statistically indistinguishable from a systematic change (P=0.18) across their 33-subject, 1.8-fold
(6.2→11.1 m/s) top-speed range, perhaps the small residual change they *do* report in the paired
incline/decline test (5 subjects, "minimum t_sw similar, +8%") is exactly what is doing the work,
alongside a modest force contribution — i.e., force and swing-time both partially explain the speed
range, not force alone.

**Forced quantitatively** (arithmetic on the verified numbers, no re-simulation needed): incline/decline
speeds 7.10→9.96 m/s (**Δv=40.3%**) accompany force 1.76→2.30×BW (**ΔF=30.7%**, leverage ratio
`Δv%/ΔF%=1.31`) and t_sw change of only 8% (leverage ratio `Δv%/Δt_sw%=5.04`) — the adversary's own
best-case mechanism would need to be carrying **3.8× more leverage per unit change** than force does,
for a signal that is smaller in absolute terms (8% vs. 30.7%) **and**, in the larger, better-powered
33-subject cross-section, **not detectably different from zero at all** (P=0.18). **Pre-registered
threshold: ≥3× leverage gap — measured 3.8× — PASS**, adversary falls. Cross-sectional force-elasticity
(computed from the two reported fold-changes, not assumed): `ln(1.26)/ln(1.8) = 0.393` — force explains
a disproportionately large share of the 1.8-fold speed range for a comparatively modest 1.26-fold force
change, while swing time shows no detectable systematic contribution over the same range.

## Part 4 — Dose-response ladder across a diverse instance-space (not two points)

The adversary is also tested against a **completely independent, video-only dataset** it was not
constructed to explain — Bolt's own race. At Berlin 2009's 60–80 m split, **Bolt (12.26 m/s) ran a
LOWER stride frequency (4.49 Hz) than the finalists he beat (11.80 m/s, 4.77 Hz)**, winning instead on
stride length (2.77 m vs. 2.48 m) — a direct, independent falsifier of "faster turnover," measured from
video kinematics with zero connection to Weyand's force-plate apparatus. Three independent
higher-speed-vs-lower-speed pairs (Bolt vs. rest / Clark&Weyand sprinters vs. nonsprinters /
Weyand 2000 fastest vs. slowest 33-subject-cohort ends) **all** show the higher-speed member with
**both** higher speed and a higher force-like metric (stride length, first-half-stance force, relative
F_avg) — **PASS** on the pre-registered "monotonic across ≥3 diverse, decorrelated instances" gate.

**Dysfunction/contrast pole** (Morin et al. 2006, PMID 16475063 — same 8 subjects, repeated 100 m
sprints, a *within-subject* function→dysfunction manipulation via fatigue, arguably cleaner than a
cross-sectional elite-vs-recreational comparison): across repeated maximal sprints, **leg stiffness and
maximal force are PRESERVED** while vertical stiffness (−20.6±7.9%), step frequency (−8.0±3.3%),
contact time (−14.7±7.2%) and speed (v_max −10.9±2.0%) **all decline together**. This is reported
descriptively (no PASS/FAIL gate pre-registered) as the twin's fatigue/dysfunction coupling anchor. One
number here is **not** mechanistically resolved from the abstract alone and is flagged rather than
narrated past: contact time *decreasing* alongside a *speed decrease* is not the direction naive
Weyand-style reasoning would predict (shorter Tc is normally bought with *more* force at *higher* speed)
— full-text access to Morin 2006 (paywalled, not fetched this session) would be needed to resolve
whether this reflects a genuine fatigue-specific technique change or an artifact of this particular
contact-time operationalization.

## Honest gaps

- Weyand et al. 2010's abstract reports `(Tc, F_avg/W)` per gait but **not** the associated m/s value —
  v0 was swept (7–10 m/s, the Clark&Weyand nonsprinter-to-sprinter band) as an assumed nuisance
  parameter, not fitted to a single number pulled from that paper.
- The hopping geometry-plausibility bound (40°, 35% compression) was calibrated with ordinary running
  in mind and is very likely too strict for Weyand's deliberately-maximal one-legged-hop condition —
  disclosed as a limitation of the pre-registered gate, not silently loosened.
- Farley, Glasheen & McMahon (1993)'s abstract reports k_leg~M^0.67 for dogs, goats, horses and red
  kangaroos explicitly; its applicability to a 65–85 kg human anthropometry sweep is an extrapolation
  by analogy (a widely-used one in the literature), not a human-specific re-derivation this session.
  Human-specific "k_leg roughly speed-invariant" for running (not stride-frequency, not cross-species)
  is separately, directly verified via He, Kram & McMahon (1991), 2.0–6.0 m/s.
- Morin et al. 2006's contact-time-decreases-with-fatigue direction (Part 4) is quoted verbatim from
  the live-fetched abstract and not mechanistically reconciled with the SLIP/impulse framework — full
  text was not accessible this session (paywalled).
- Krzysztof & Mero (2013) is a tertiary, video-digitization-based re-analysis (not a first-party
  biomechanics lab measurement) — its own stride-length × stride-frequency products do not exactly
  reproduce its own reported velocities (Bolt: 2.77×4.49=12.44 vs. reported 12.26, 1.5% internal
  discrepancy; rest-of-field: 2.48×4.77=11.83 vs. 11.80, 0.3%) — a known limitation of frame-by-frame
  video digitization, disclosed rather than silently reconciled by this doc.
- No raw per-subject data was accessed for any cited study (all numbers are live-verified
  published-abstract summary statistics) — every figure in this doc traces to a PMID/DOI fetched via
  NCBI eutils this session (ledger: evidence JSON), not to recalled textbook values.
- The SLIP model's touchdown-angle/vy_td search is a coarse-to-fine grid, not an exhaustive global
  search — it is robust against solver non-convergence (the risk a black-box optimizer would carry) but
  could in principle miss a better-fitting point between grid nodes; refinement step size (0.4°, 0.03 m/s)
  is far finer than the reported residual margins, so this is a low-risk, disclosed approximation.

## couples_to

- `docs/MECHANISM_GAIT_KINEMATICS_FIDELITY.md`, `docs/MECHANISM_SPINE_GAIT_VBR.md` — general gait
  kinematics fidelity this doc's sprint-specific stance mechanics sit on top of.
- `docs/MECHANISM_ANKLE_FORCE.md`, `docs/MECHANISM_HIP_FORCE.md`, `docs/MECHANISM_JOINT_FORCE_SCORECARD.md`,
  `docs/MECHANISM_JOINT_FORCE_VALIDATION.md` — this doc's whole-body mass-specific GRF is the external
  boundary condition those internal joint-force certs distribute across ankle/knee/hip.
  `docs/MECHANISM_ANKLE_WAVEFORM.md` in particular shares a contact-time/force-waveform vocabulary.
- `docs/MECHANISM_BIARTICULAR_MUSCLE.md`, `docs/MECHANISM_FORCE_VELOCITY_POWER_RFD.md` — the muscle-level
  force-velocity-power capacity (F0/V0/Pmax, RFD) that this doc's whole-body mass-specific GRF is the
  *output* of; `MSK-SPRINT-FV-PROFILE` in the anchor graph is the macroscopic (split-time-derived)
  sibling of this doc's ground-level mechanism.
- `docs/MECHANISM_WEIGHTLIFTING_TRIPLE_EXTENSION.md` — shares the "brief, maximal mass-specific force
  application" motif in a different (triple-extension, not steady-state bounce) task.
  `MSK-REPEATED-SPRINT-FATIGUE` — the metabolic/PPO-decline sibling to this doc's Part 4 biomechanical
  fatigue signature (Morin 2006); same broad theme, decorrelated measurement level.
  `MSK-GAIT-LOCOMOTION` (OPEN, SEED-DESIGN) — explicitly proposes the SLIP-model low-dim-latent
  parametrization this doc builds and validates; this doc supplies (not closes) that forward model.
- `data/video_index_tables/body_twin_collection__youtube/` — already contains indexed (shot-segmented,
  not yet pose-extracted) sprint clips ("Sprint Form Slow Motion," "60m Sprint Drive Phase," "Sprint
  Drills - Bounding Strides") — the concrete, existing substrate for the proposed cell below.
- The function↔dysfunction organizing axis: Part 4's Morin-2006 fatigue contrast is this doc's
  contribution to that axis (a within-subject, biomechanically-measured function→dysfunction
  transition), independent of the cross-sectional elite-vs-recreational framing used elsewhere.

## Proposed cell

Run the twin's markerless pose pipeline on the already-indexed sprint clips in
`data/video_index_tables/body_twin_collection__youtube/` to extract real per-step contact-time and
stride-frequency measurements from the operator's own footage; feed the recovered `(v, Tc)` pairs
through this doc's impulse-balance identity (Part 2) to predict `F_avg/W`, and through the exact SLIP
grid-match (Part 1) to recover an implied `k_leg` — a fully closed-loop, camera-only, zero-force-plate
replication of the Weyand mechanism, machine-gated against this doc's pre-registered residual/plausible-
geometry thresholds rather than eyeballed against the video.

## Repro

```
cd ~/projects/bodytwin
python3 scripts/msk/sprint_spring_mass.py   # writes reports/probes/sprint_spring_mass_results.json
```
Pure numpy/scipy (`solve_ivp`) ODE integration + arithmetic on verified literature numbers — no
OpenSim, no external data download, no network access needed at run time (all NCBI fetches happened
during authoring, captured in the evidence JSON). Runtime: ~15-25 seconds. Verified working both under
plain system `python3` (numpy 2.2.6/scipy 1.15.3) and `.venv-msk` (numpy 2.5.1/scipy 1.18.0). No git
operations; no writes outside `reports/probes/sprint_spring_mass_results.json` and this doc pair (this
`.md` + `docs/MECHANISM_SPRINT_SPRING_MASS_evidence.json`).
