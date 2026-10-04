# MECHANISM METABOLIC COST — muscle-energetics capability (2026-07-21)

Couples the musculoskeletal twin to whole-body energetics/calorimetry (COORDINATOR.md Sec.5: "optics
= the eye biology's forward model, **battery/thermal verticals = metabolism/calorimetry
cert-twins**, compute = the mind"). Computes the twin's metabolic cost of transport (COT) and
gross metabolic power for the subject2/walking1 trial, from the **already-validated** Static
Optimization output — **no re-solve**. Script: `scripts/msk/metabolic_cost.py`. Evidence:
`data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`.

## Method, in one paragraph

Loads `data/msk_smoketest/subject2_walking1/static_optimization/so/` (activation.sto +
force.sto, convergence-gated PASS per `static_opt_knee_results.json`) and the IK kinematics that
fed it (`walking1.mot`, 158 rows — time grid matches the SO output **exactly**, max abs diff
0.0s, verified not assumed). At each of the 158 frames: sets coordinate values/speeds (deg→rad,
central-difference) and each of the model's 80 real `Muscle` actuators' activation directly from
SO's own activation.sto (never re-derived); calls `model.equilibrateMuscles(state)` to solve
fiber length for quasi-static force equilibrium (the same per-instant, no-cross-frame-continuity
assumption Static Optimization itself makes); realizes Velocity+Dynamics; queries OpenSim's own
`Umberger2010MuscleMetabolicsProbe` and `Bhargava2004MuscleMetabolicsProbe` (vendored in
`.venv-msk` opensim 4.6 — used directly, no custom energetics equations written). Only the 80
genuine `Muscle` objects are registered with the probes; the model's other 31 DOF actuators
(13 ideal torque actuators for lumbar+arms, 6 pelvis residuals, 12 joint reserves — enumerated
from force.sto's own 130-column header) carry no metabolic cost here — a disclosed under-
estimate source (real trunk/arm muscles cost real energy; this model represents them as ideal
actuators).

## Citations — verified LIVE, not recalled (a real catch)

| model | citation as verified | note |
|---|---|---|
| Umberger 2010 | Umberger BR (2010) "Stance and swing phase costs in human walking," *J R Soc Interface* 7(50):1329–40. **PMID 20356877**, PMC2894890 | **The task brief's PMID (20356878) is WRONG** — resolves to an unrelated paper ("Intercellular mechanotransduction during multicellular morphodynamics," Kim/Dooling/Asthagiri 2010, same journal). Off-by-one-digit error, verified via NCBI eutils esummary, corrected here, not propagated. |
| Bhargava 2004 | Bhargava LJ, Pandy MG, Anderson FC (2004) "A phenomenological model for estimating metabolic energy consumption in muscle contraction," *J Biomech* 37(1):81–8. **PMID 14672571** | Task brief's PMID verified CORRECT. |
| External anchor | Koelewijn AD, Heinrich D, van den Bogert AJ (2019) "Metabolic cost calculations of gait using musculoskeletal energy models, a comparison study," *PLoS ONE* 14(9):e0222037. **PMID 31532796**, PMC6750598 | Open-access, fetched full text (not just abstract). An automated PMC-elink pass initially **mis-paired** this PMID with PMC2894890 (Umberger's own paper) — caught by fetching both PMC landing pages directly and reading their actual titles, not by trusting the first machine answer. |

Recall discipline: this project's own prior finding is a **measured ~62% citation-drift rate**
from memory. Every number/PMID used below was fetched live (NCBI eutils, PMC full text, or direct
Python object introspection of the installed OpenSim build) this session, not recalled.

## Headline result (primary, pre-registered configuration)

Primary configuration: rigid tendon (`ignore_tendon_compliance=True` on all muscles — chosen
*before* seeing the final COT number, because it reproduces SO's own force.sto more closely than
the model-default elastic tendon, see Force cross-check below) + each probe's own Fmax-derived
default muscle mass + uniform `ratio_slow_twitch_fibers=0.5`.

| quantity | Umberger2010 | Bhargava2004 |
|---|---:|---:|
| net rate (muscle-only, no basal) | 8.181 W/kg | 6.599 W/kg |
| gross rate (+1.2 W/kg basal) | 9.381 W/kg | 7.799 W/kg |
| **COT net** | **7.683 J/kg/m** | **6.198 J/kg/m** |
| COT gross | 8.810 J/kg/m | 7.325 J/kg/m |

Trial: subject2/walking1, mass 78.2 kg, whole-body-COM-derived forward speed **1.065 m/s** over
1.671 m / 1.57 s (cross-checked: gravity is (0,−g,0) in this model, pelvis_tz net drift only
0.025 m vs pelvis_tx's 1.66 m over the same window — confirms the forward-axis assumption).
Two-model (Umberger vs Bhargava, same reconstructed states) spread: **21.4%** — comparable in
order of magnitude to Koelewijn's own externally-reported ~28–33% inter-model RMS error (below).

## External anchor comparison — and an honest SURPRISE

Koelewijn et al. 2019 (healthy adults, n=12, 6F/6M, age 24±5y, mass 70±12kg), treadmill walking
at 1.3 m/s, **NET** (indirect-calorimetry-measured, resting-subtracted) cost of transport —
quoted verbatim from their abstract, consistently across two independent fetch passes:

- downhill (−8% grade): **2.0 ± 0.40 J/kg/m**
- uphill (+8% grade): **5.9 ± 0.34 J/kg/m**
- level (0% grade): not cleanly extractable as a single verbatim number from two automated fetch
  passes (got 1.1 and ~3.5 J/kg/m inconsistently) — rather than pass off an uncertain extraction
  as a quote, level COT here is **our own linear interpolation** of the two verbatim, twice-
  consistent grade values: (2.0+5.9)/2 = **3.95 J/kg/m**, flagged explicitly as an interpolation,
  not Koelewijn's own reported number. Sits inside the task brief's own prior range (2–4 J/kg/m).
- Their own model-vs-measured RMS error at level/1.3 m/s: UMBE03 (Umberger's 2003 predecessor)
  1.11 J/kg/m, BHAR04 (Bhargava 2004) 1.31 J/kg/m → **~28–33% relative error**, independently
  verifying the task brief's own "±30% spread" prior from the primary comparison-study source.

**PRE-REGISTERED gate (set before computing the trial-averaged number): twin COT_net should be
BELOW the measured anchor** (SO's effort-minimization discards antagonist co-contraction, which
costs real energy — Koelewijn et al. state this explicitly: *"Most models... generally
underestimate the metabolic cost. A possible reason is the approach used to find the muscle
activations by minimizing muscular effort. Then, muscular co-contraction was disregarded."*).

**Result: SURPRISE, not suppressed.** Twin COT_net (7.68 J/kg/m Umberger, 6.20 Bhargava) is
**ABOVE** the 3.95 J/kg/m anchor — ratio 1.95 / 1.57 — the opposite of the pre-registered
expected direction. Symmetric disclosure: Umberger's OWN 2010 paper (PMC2894890, direct fetch)
independently reports the SAME opposite-direction pattern for its own (non-SO, dynamic-
optimization-generated) muscle activations: *"Total muscle energy consumption in the model was
higher than the net energy expenditure in the subjects."* I.e. the same cost equations can over-
or under-shoot depending on how the driving muscle activations were generated — the under-
prediction direction is a property of SO/effort-minimization specifically, not guaranteed by the
Umberger/Bhargava heat-rate equations themselves, and this twin's SURPRISE result had to be
force-diagnosed, not hand-waved away either as confirming or refuting the prior.

## OODA-forced diagnosis of the surprise (the actual investigative work)

Ruled OUT (checked directly, not assumed):
- **Activation-state corruption** through `equilibrateMuscles`: verified activation is preserved
  bit-identical before/after (0.3 in, 0.3 out) — not a bug.
- **`muscle_effort_scaling_factor` misconfiguration**: default is 1.0 (verified by direct object
  introspection) — not the issue.
- **`aerobic_factor` misconfiguration**: default 1.5 is documented (OpenSim header, fetched) as
  *"S=1.0 for primarily anaerobic conditions and S=1.5 for primarily aerobic conditions"* — 1.5
  is the physiologically CORRECT choice for walking (aerobic), not an inflation bug.
- **Per-frame reconstruction noise in fiber velocity**: curvature-ratio smoothness check across 5
  muscles shows fiber_length(t) tracks the purely-geometric MTU_length(t)'s curvature almost
  exactly (e.g. vaslat_r: 0.0101 vs 0.0101; tibpost_l: 0.0804 vs 0.0805) — the equilibrium solve
  is well-conditioned and smooth, not jagged/noisy. PASS.
- **Force-reconstruction fragility (Step 3) driving the total**: a reconstruction-FREE proxy
  (activation × muscle_mass × Bhargava's own act/maintenance constants, no `equilibrateMuscles`
  call at all) gives 5.65 W/kg — same order of magnitude as the full pipeline's 8.18/6.60 W/kg —
  confirms the headline isn't an artifact of the force-mismatch named below.

Identified and QUANTIFIED (the real drivers):
1. **Muscle-mass over-estimate.** OpenSim's standard Fmax/specific_tension×density×optimal_fiber
   _length mass formula (the field-standard technique, and what these probes use internally by
   default) gives 45.6 kg summed over the 80 registered muscles — **58.3% of this subject's body
   mass**, exceeding even whole-body skeletal muscle mass estimates (~40% of body mass, Janssen
   et al. 2000) for just a lower-limb+hip subset. Most likely mechanism: this (pre-existing,
   unchanged) model's per-muscle Fmax values were inflated by the subject-specific OpenSim
   scaling step. Correcting summed mass down to a literature-reasoned ~17.2 kg (whole-body-%
   derived, not directly cited — flagged): **−23.2%** on the headline number.
2. **Tendon-compliance treatment.** Rigid tendon (primary choice, justified by the force cross-
   check below) vs the model's default elastic tendon: **−21.6%**.
3. **Combined** (both corrections applied together, run jointly not just multiplied): net rate
   drops to **3.92 W/kg → COT 3.68 J/kg/m → ratio to the Koelewijn anchor = 0.93** — landing
   almost exactly on the external anchor, 7% below it, which IS the theoretically-expected
   direction. This does not mean "3.68 J/kg/m is the corrected true answer" (the 17.2 kg mass
   target is itself a rough, disclosed estimate) — it means the primary configuration's surprise
   is **explained**, not mysterious: two identified, individually well-motivated modeling choices
   account for essentially all of the gap. Honest range to carry forward: **COT_net spans ~3.7–
   7.7 J/kg/m depending on these two disclosed choices, bracketing the 3.95 J/kg/m external
   anchor from both sides.**
4. **Term-by-term ablation** (which heat-rate term actually drives the total): activation+
   maintenance heat rate = 48.6% of the total, basal = 29.3%, shortening/lengthening = 22.0%,
   mechanical work rate = 0.2% (negligible). This is *why* the mass correction (which only scales
   the mass-proportional maintenance term) has a real but partial effect, not a 1:1 effect.
5. **`ratio_slow_twitch_fibers` sensitivity measured at ~0** (4.2e-10, machine precision) across
   a 0.3/0.5/0.7 sweep, in this probe configuration (`use_Bhargava_recruitment_model=True`
   default). Verified exactly, not just "small" — reported honestly as an unresolved-but-
   precisely-measured curiosity; the mechanism (why it is *exactly* insensitive, not just
   weakly sensitive) was not fully traced into the C++ recruitment-model implementation within
   this task's budget.

## Force-reconstruction cross-check (named limitation, not hidden)

Does the per-frame quasi-static state built here reproduce SO's own reported force.sto? Median
relative error (activation ≥ 0.05, the metabolically-meaningful regime): **7.4%**, p90 65.7%.
Error concentrates in specific muscles: `tibpost_l` (467% median rel. err.), `addbrev_l` (219%),
`perlong_r` (178%), `addbrev_r` (170%), `grac_l` (155%) — likely a wrapping-surface path-length
sensitivity between SO's own fast/simplified internal force model and the full nonlinear
`Millard2012EquilibriumMuscle` equilibrium used here; not resolved further within this task's
scope, but shown (Step 4's reconstruction-free proxy) not to dominate the aggregate total.

## Scope, honest gaps

- Single subject (subject2), single trial (walking1), one speed (1.065 m/s) — not a generality
  claim.
- 31 non-muscle DOF actuators (trunk/arm ideal torque actuators, pelvis residuals, joint
  reserves) carry zero metabolic cost here — a real, disclosed under-estimate source.
- "Gross" here = probe net output + OpenSim's own default basal term (1.2 W/kg, cross-checked
  against 1 MET = 3.5 mL O2/kg/min = 1.162 W/kg, Wikipedia "Metabolic equivalent of task" —
  textbook-grade, not primary-literature, source, flagged as such). This is NOT the same
  "resting baseline" Koelewijn measured directly per-subject; the gross-vs-gross comparison
  carries this additional caveat on top of the speed mismatch (twin 1.065 m/s vs anchor's
  1.3 m/s — COT is the speed-normalized quantity and is compared directly; W/kg gross rate is
  speed-dependent and not directly comparable across speeds without this caveat).
- Whole-body muscle-energetics models generally carry a documented ±30% spread vs measured
  indirect calorimetry (now independently verified via Koelewijn's own RMS-error table, not just
  assumed from the task prior) — this twin's numbers should be read at that resolution, not as
  point predictions.

## Repro

```
cd ~/projects/bodytwin
.venv-msk/bin/python scripts/msk/metabolic_cost.py
```
Writes `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json`.
Runtime ~2–3 min (primary 158-frame pass ×2 probes + 8 additional sensitivity/ablation passes on
a 53-frame subsample). No git operations; reads the shared `LabValidation_withVideos` /
OrthoLoad-adjacent data tree read-only; writes only under `data/msk_smoketest/subject2_walking1/`.
