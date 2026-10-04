# MECHANISM POSTURAL BALANCE CONTROL — ankle/hip/step strategy hierarchy + the function-dysfunction fall-risk axis (2026-07-22)

Companion to (**not** a rebuild of) `docs/MECHANISM_VESTIBULAR_BALANCE.md` Part 2 and the executed
`ORG-POSTURAL-BALANCE` graph node. This doc's job is the two pieces the task brief asks for that genuinely
do **not** already exist in this repo: (1) a from-first-principles geometric proof that passive ankle
stiffness alone cannot hold quiet stance, forced against a second, independent stiffness measurement; (2)
the ankle-hip-step **strategy hierarchy** as a function of perturbation size; (3) the function-vs-dysfunction
**fall-risk** axis with real prospective cohorts. All PMIDs below were checked live this session via NCBI
eutils (esearch -> esummary -> efetch abstract) — none taken from recall. No new script was written; the one
new numeric result (Part 1) is a six-line, dependency-free arithmetic check on numbers already sitting in
this repo's own JSON, reproduced verbatim in the Repro section.

## Scope note — what already exists, re-verified not rebuilt

`docs/MECHANISM_VESTIBULAR_BALANCE.md` Part 2 (`scripts/msk/postural_sway_pendulum.py`,
`reports/probes/postural_sway_pendulum_results.json`) already built and machine-verified: the single-link
inverted-pendulum-plus-delay model (`I*theta''=mgh*theta-Kp*theta_hat(t-delay)-Kd*thetadot(t-delay)+xi`),
Loram & Lakie's insufficiency number, Peterka's delay-not-damping mechanism, a genuine stability-margin-
governs-amplitude structure (27.6x variance growth from 50ms to 200ms delay, machine-crosschecked
exact-frequency-domain-vs-Monte-Carlo after a real instability bug was caught and fixed), the internal
`ORG-POSTURAL-BALANCE` real-HBEDB-data anchor, and an honest negative on quantitatively reproducing the
Romberg velocity/RMS metric-specific dissociation. **This is not redone here** (isolation discipline: re-verify
the load-bearing citations, don't re-run the simulation). The two most decisive citations that model depends
on (Loram & Lakie 2002, Peterka 2002) were independently spot-rechecked this session via NCBI esummary — both
bibliographically match the existing doc's citation exactly, zero drift found.

## Part 1 — Why active feedback is load-bearing: the passive-only limit is provably, quickly unstable

**Claim (C) / threshold**: a purely passive ankle spring (no neural delay, no active gain — i.e.
`T_A(t)=Kp_passive*theta(t)` instantaneously, nothing else) at either literature-measured intrinsic-stiffness
fraction of critical (`mgh`) produces a real, strictly-positive characteristic root (unconditional exponential
instability), pre-registered threshold: `Kp_passive/mgh < 1.0`. **Adversary (leaning-positive claim, so the
adversary is the confound that could rescue passive-sufficiency)**: co-contraction or a different measurement
paradigm might push the *true* intrinsic stiffness at or above criticality, making Loram & Lakie's 91% an
under-estimate. Forced via an independent replication using a **different lab, different disturbance-size
regime, different device**:

| source | method | intrinsic stiffness / critical | verified |
|---|---|---|---|
| Loram ID, Lakie M (2002). *J Physiol* 545(3):1041-53. **PMID 12482906**, DOI 10.1113/jphysiol.2002.025049, PMC2290720 | very small (imperceptible) perturbations, sway-referenced measurement, Birmingham UK | **91 +/- 23%** | Bibliographic + numbers re-confirmed live this session (matches existing in-repo citation exactly, no drift) |
| Casadio M, Morasso PG, Sanguineti V (2005). *Gait Posture* 21(4):410-24. **PMID 15886131**, DOI 10.1016/j.gaitpost.2004.05.005 | motorised footplate, 1-degree step disturbances, EMG-checked to exclude reflex contamination, Genova IT | **64 +/- 8%** | Abstract fetched live this session — NEW citation, not in the existing doc |

The adversary is **forced to its strongest fair form and falls harder, not softer**: an independent group,
independent decade, independent device, and an explicitly *larger* perturbation (1 degrees vs Loram & Lakie's
imperceptible probe) — the direction a genuine adversary would need to push the number *up* toward or past
100% — instead measures a **lower** fraction (64% vs 91%). Casadio et al. themselves read the two numbers as
the lower/upper bound of one real range (smaller test-disturbance -> higher apparent stiffness), not as
contradictory estimates — and **both bounds of that range are still below the critical value**. Two
decorrelated real measurements, same conclusion, one even more decisive than the other.

### Geometric derivation (eigenvalue sign, not narrative) + machine-computed consequence

For the homogeneous passive-only system, `I*theta'' = (mgh - Kp_passive)*theta` is a linear 2nd-order ODE
whose state-matrix `[[0,1],[omega^2,0]]` has eigenvalues `+/-omega`, `omega = sqrt((mgh-Kp_passive)/I)`. Since
both measured fractions give `Kp_passive < mgh`, `omega` is **real** (not imaginary) — an unconditionally
growing mode exists, not a bounded oscillation. Using this repo's own already-verified real subject parameters
(`mass=78.2kg` from `metabolic_cost_results.json`, `h=0.969m` = Winter's 0.56x1.73m stature fraction,
`mgh=743.2 Nm` — all three numbers reused verbatim from `postural_sway_pendulum_results.json`, not
re-derived):

| stiffness estimate | omega (1/s) | doubling time (s) | time to grow 0.005rad (~0.3deg, realistic quiet-stance amplitude) to 0.2rad (~11.5deg, an unrecoverable lean) |
|---|---:|---:|---:|
| Loram & Lakie 91% | 0.954 | 0.73 | **3.87 s** |
| Casadio et al. 64% | 1.909 | 0.36 | **1.93 s** |

**PASS** (claim confirmed, both independent measurements): a passive-only ankle would collapse within 2-4
seconds with zero active correction — nowhere near the minutes-long duration of real quiet standing. This is a
genuinely new, machine-computed number (not printed by the existing pendulum script, which reports variances
under the *active* model only) — six lines of dependency-free Python, shown in Repro, using only numbers
already sitting in this repo's own JSON plus the two independent literature fractions.

## Part 2 — Sway magnitude/frequency reproduction: what's PASS, what's OPEN (reused, not rerun)

The existing model's **frequency-domain** plausibility check already PASSES, unfitted: predicted PSD peaks at
0.20 Hz with 99.9% of power below 1 Hz, matching the well-documented real-world fact that quiet-stance COP
sway concentrates below ~1-2 Hz. The **absolute magnitude** (COP RMS in real mm) is explicitly an **OPEN,
disclosed gap in the existing model** — its own honest-gaps list states "this script does not calibrate
absolute COP RMS in cm against a literature value" because the noise input is in arbitrary units (only
EO/EC *ratios* are asserted as falsifier-relevant). This is re-stated here, not re-litigated: I did **not**
attempt a new absolute-unit calibration of that simulation. Instead, two real-unit anchors close to (but not
identical to) the RMS-in-mm target were live-verified this session:

- **Internal, already byte-verified** (`ORG-POSTURAL-BALANCE`, PhysioNet HBEDB, n=163/644 trials): mean COP
  **velocity** young=0.905 cm/s, old=1.274 cm/s (ratio 1.409, cross-checked by direct division =1.4077,
  matches the reported rounded value).
- **Prieto TE, Myklebust JB, Hoffmann RG, Lovett EG, Myklebust BM (1996).** *IEEE Trans Biomed Eng*
  43(9):956-66. **PMID 9214811**, DOI 10.1109/10.532130 — abstract fetched live this session (a NEW citation).
  n=20 healthy young (21-35y) vs n=20 healthy elderly (66-70y), EO/EC quiet stance. Headline finding, quoted:
  "Mean velocity of the COP was the **only** measure that identified age-related changes in both eye
  conditions, and differences between eye conditions in both groups."

**Convergence, not tautology**: this is an independent dataset, decade, cohort, and research group finding the
*exact same metric-specificity pattern* — velocity-type measures carry the age/sensory signal, other COP
measures do not — that the internal `ORG-POSTURAL-BALANCE` cell already found on real HBEDB data (RMS/area
Romberg effect null, p=0.58/0.22, while velocity/path-length are significant at p<1e-4). Two decorrelated real
datasets, 20+ years apart, agreeing on *which metric is the reliable one* is a stronger anchor than either
alone. **Honest gap, explicitly not closed**: Prieto 1996's exact mm/cm-per-second numbers are in its tables,
not its (paywalled) abstract — not extracted live; the pendulum model's own absolute-unit calibration remains
open exactly as the existing doc already disclosed it.

## Part 3 — Ankle -> hip -> step: a torque continuum shaped by biomechanical constraints, not three discrete boxes

**Claim / pre-registered threshold**: postural strategy shifts from ankle-dominant to hip-recruiting to
stepping as perturbation magnitude/task difficulty increases, verified via (a) a real EMG/kinetic dissociation
across an independent perturbation-like variable, and (b) a distinct, age-linked third (stepping) tier.

| source | independent variable | finding |
|---|---|---|
| Horak FB, Nashner LM (1986). *J Neurophysiol* 55(6):1369-81. **PMID 3734861**, DOI 10.1152/jn.1986.55.6.1369 | support-surface **length/geometry** (normal vs short-relative-to-foot-length) | Ankle strategy (normal surface): 73-110ms EMG latency, activation begins ankle and radiates distal-to-proximal, produces ankle torque, restores balance by moving COM about the ankle. Hip strategy (short surface): same latencies, but antagonist muscles activate proximal-to-distal (opposite sequence), little/no ankle torque, produces horizontal shear force, motion focused at the hip. Intermediate surface lengths -> complex combined patterns. |
| Runge CF, Shupert CL, Horak FB, Zajac FE (1999). *Gait Posture* 10(2):161-70. **PMID 10502650**, DOI 10.1016/s0966-6362(99)00032-6 | backward translation **velocity** on a FLAT surface, 5-55 cm/s | Joint-torque analysis: a hip flexor torque **adds onto** the ankle plantarflexor torque as velocity increases — a continuum. **"Hip torque without accompanying ankle torque (pure hip strategy) was not observed."** |
| McIlroy WE, Maki BE (1996). *J Gerontol A Biol Sci Med Sci* 51(6):M289-96. **PMID 8914501**, DOI 10.1093/gerona/51a.6.m289 | perturbation exceeding what ankle/hip torque alone can recover (600ms platform translation) | 98% of all trials required a step at all. Elderly (n=9, 65-81y) needed **multiple** steps in **63%** of trials vs young (n=5, 22-28y) described as "twice as likely" (implies young ~31%, my arithmetic from the stated ratio, not a number printed verbatim in the abstract — flagged as such) — third tier engaged more, and with age-linked deficit, exactly where ankle+hip torque is exhausted. |

**Forced adversary and honest refinement**: the naive textbook reading — three cleanly discrete strategy
"boxes" selected by perturbation size — is the version I was tempted to just accept. Runge et al. 1999 is
itself an explicit test of this adversary (does hip strategy appear as a genuinely separate, discrete mode, or
as a continuum?) on a *fixed, flat* surface, varying only translation velocity. Result: **the adversary
partially falls** — hip torque is real and does scale with perturbation size, but it never appears alone
(no "pure hip" case with zero ankle torque was found); this is a **continuum**, with the surface-geometry
manipulation in Horak & Nashner 1986 acting as a *separate, additional* constraint that can suppress the
ankle-torque option almost entirely (short/compliant surface). **PASS on the qualitative structure** (a
real, EMG/kinetic-confirmed transition exists, across two independent variables — surface geometry and
translation velocity — from two different decades/groups, plus a genuinely distinct, age-linked stepping
tier); the brief's own "hierarchy" framing is **refined, not simply confirmed**: better described as a
torque continuum whose *discreteness is constraint-dependent* (surface geometry can force a sharp
ankle-torque cutoff; velocity alone does not).

## Part 4 — The function<->dysfunction fall-risk axis: sway prospectively predicts falls, but not as one universal clean law

**Claim / pre-registered threshold**: elevated sway (esp. velocity-type measures) prospectively predicts
future falls in older adults above chance, in a genuinely prospective (not merely concurrent/retrospective)
design. **Adversary (leaning-positive — "more sway = more falls" is tempting to accept uncritically as a
universal law)**: force with a **more frail, different-instance-space cohort** — does the relationship hold
with the same strength there too?

| source | cohort | design | finding |
|---|---|---|---|
| Maki BE, Holliday PJ, Topper AK (1994). *J Gerontol* 49(2):M72-84. **PMID 8126355**, DOI 10.1093/geronj/49.2.m72 | n=100, ambulatory/independent elderly, age 62-96 | 1-year **prospective** fall monitoring after force-plate balance testing | Lateral spontaneous-sway amplitude (eyes-closed) = single best predictor of future falling, "moderate accuracy," predictive **even in subjects with no recent fall history** (a genuinely prospective, not reverse-causation, signal) |
| Fernie GR, Gryfe CI, Holliday PJ, Llewellyn A (1982). *Age Ageing* 11(1):11-6. **PMID 7072557**, DOI 10.1093/ageing/11.1.11 | n=205, **institutionalized** geriatric, mean age 81.8y (more frail instance-space) | 1-year prospective fall monitoring | Mean sway speed significantly greater in fallers (p<0.05) — direction confirmed — **but** "the difference was less than might have been expected" and **no dose-response trend** between sway magnitude and fall frequency was found |

**The adversary partially falls, honestly reported, not force-fit to a clean PASS**: direction replicates in
both a clean (ambulatory-independent) and a harder (institutionalized-frail) cohort, but effect cleanliness
degrades in the frailer population — exactly the kind of forced, diverse-instance-space test this method
requires, and exactly the kind of result that should **not** be silently smoothed into one uniform "sway
predicts falls" headline. **Convergent internal + cross-decade anchor**: the internal `ORG-POSTURAL-BALANCE`
cell's own real, byte-verified HBEDB measurement (age-sway COP-velocity ratio old/young **1.409**,
Mann-Whitney p=1.9e-6) agrees on direction and on Prieto 1996's independent finding (Part 2) that velocity
specifically, not every COP measure, is what tracks the effect — three decorrelated real datasets across four
decades (1982/1994/1996/HBEDB-era), same direction, same metric-specificity texture.

**Honest gap, not closed**: none of these studies (nor this doc) establish that sway is a *mechanistic*
readout of the delayed-feedback control loop itself rather than an epiphenomenal marker of general frailty
that independently also causes falls — Maki 1994's "predictive even with no fall history" result argues
against pure reverse-causation but does not by itself prove the mechanistic-control-loop link. A plausible,
disclosed-not-verified actuator-level mechanism this repo already has real data for: `MSK-SARCOPENIA-COMPOSITE`
(MEASURED-B, EXECUTED 2026-07-18) found a real dynapenia dissociation (grip strength declines faster than
lean mass with age, NHANES n=2810) — reduced torque-generating capacity would directly cap the active-feedback
gain terms (`Kp`,`Kd`) in the Part-1/existing pendulum model, a concrete, testable but **not yet tested**
coupling between two already-measured cells.

## Part 5 — Sensory removal / Romberg generalizes across channels (reused, not re-verified this session)

Not re-verified here (already real, byte-verified, in-repo): `ORG-POSTURAL-BALANCE`'s own Romberg ratios
(eyes-closed/open) 1.155 (firm surface, p=1.2e-4) and 1.304 (foam surface, p=4e-21), with the foam/firm ratio
(3.109, p=5e-28) showing the surface (proprioceptive) manipulation has ~10x the effect of the visual one — a
real sensory hierarchy, not just "removing any one sense increases sway" flattened to a single number. Peterka
2002 (**PMID 12205132**, spot-rechecked live this session, matches existing citation) supplies the general
mechanism: sensory reweighting — reliance shifts toward the remaining channels as one is removed/degraded, and
stiffness/gain measurably changes with the reweighting, not just noise amplitude (already used by the existing
pendulum model to motivate its Kp-change mechanism test).

## Part 6 — Video-markerless coupling: real anchor confirmed, not yet built here for quiet stance

This repo's own `scripts/msk/pose_extract.py` (`docs/MECHANISM_POSE_PIPELINE.md`) produces real per-frame 3D
landmarks (33 MediaPipe points incl. ankle/hip/shoulder, 228/264 clips passing its own gate suite) but does
**not** currently compute COP/COM or run on quiet-stance clips — it targets gait/plyometric clips for OpenSim
IK. The graph's own `SNS-VESTIBULAR-BALANCE` node (SEED-DESIGN, **not yet executed**) already proposes
exactly this extension — camera-only recovery of Peterka sensory weights via a markerless COP pipeline citing
**PMC12189434**. That PMC ID was independently confirmed live this session (not taken on the existing node's
word): Feng R, Ugbolue UC, Yang C, Liu H (2025). "Estimation of Three-Dimensional Ground Reaction Force and
Center of Pressure During Walking Using a Machine-Learning-Based Markerless Motion Capture System."
*Bioengineering (Basel)* 12(6) — a real, live-verified paper confirming markerless-video COP/GRF estimation is
methodologically real (for walking), not vaporware. **Honest gap**: neither that paper nor this repo's
pipeline has been run on quiet standing; `SNS-VESTIBULAR-BALANCE` remains SEED-DESIGN, correctly so.

## Symmetric QC — what is and is not claimed

**Claimed (PASS, machine-checked, externally anchored)**:
1. Passive-only ankle stiffness at either of two independent, live-verified measurements (91%, 64% of
   critical) yields an unconditionally unstable linear system (real positive eigenvalue) that collapses within
   2-4 seconds — a new, geometric, machine-computed number using this repo's own already-verified parameters.
2. A real ankle-to-hip torque continuum exists, confirmed across two independent perturbation-defining
   variables (surface geometry, translation velocity) from two decades/groups; a third, distinct, age-linked
   stepping tier is engaged when the first two are exhausted.
3. Elevated sway prospectively predicts future falls in a clean ambulatory-elderly cohort, converging in
   direction (not in dose-response cleanliness) with a frailer institutionalized cohort and with this repo's
   own internal, byte-verified age-sway measurement and an independent cross-decade metric-specificity finding.

**NOT claimed (explicitly OPEN or partial)**:
4. Absolute COP RMS-in-mm reproduction from the mechanistic model — the existing pendulum model's frequency
   shape passes, its noise-scale is admittedly arbitrary; not newly calibrated here.
5. The ankle/hip/step "hierarchy" as three cleanly discrete boxes — refined to a torque continuum whose
   discreteness is constraint- (surface-) dependent, per Runge et al.'s own forced test.
6. Sway as a *mechanistic*, not merely statistical/epiphenomenal, readout of the control loop — plausible
   (sarcopenia/dynapenia coupling identified) but not tested.
7. Camera-only quiet-stance COP/COM recovery — the literature anchor for the general method is real and
   live-confirmed; this repo has not run it on quiet stance; `SNS-VESTIBULAR-BALANCE` correctly stays
   SEED-DESIGN.

`couples_to`: `docs/MECHANISM_VESTIBULAR_BALANCE.md` (Part 2 — the model this doc extends, not duplicates),
`ORG-POSTURAL-BALANCE` (real HBEDB anchor, unedited), `SNS-VESTIBULAR-BALANCE` (SEED-DESIGN camera-COP sibling,
same video-markerless framing), `SENS-VESTIBULAR-BALANCE` (clinical UVL/BVL laterality fall-risk — a different
specific mechanism, same fall-risk theme), `MSK-SARCOPENIA-COMPOSITE` (plausible actuator-gain coupling,
untested), `docs/MECHANISM_POSE_PIPELINE.md`/`scripts/msk/pose_extract.py` (the raw video landmark stream any
future COP-from-video extension would consume), `docs/MECHANISM_PROPRIOCEPTION.md` (the afferent/Ia channel
conceptually feeding the pendulum model's noisy position estimate).

## Honest gaps (consolidated)

- Absolute-unit (mm) COP-RMS calibration of the mechanistic model remains open (inherited gap, not newly
  closed).
- The ~150-200ms delay parameter in the existing pendulum model is still an unpinned round estimate; this
  session did surface a more precise, live-verified, directly relevant number — Horak & Nashner 1986's 73-110ms
  automatic-postural-response EMG latency — but that measures EMG *onset*, not the same quantity as a
  modeled total loop delay (which plausibly includes electromechanical delay on top); offered as a candidate
  refinement, not substituted in.
- McIlroy & Maki 1996's exact young-group step-repeat percentage is my own arithmetic from the stated "twice
  as likely" ratio (~31%), not a number printed verbatim in the abstract.
- Prieto 1996's exact mm/cm-per-second values are in its tables, not its abstract — not extracted live
  (paywalled).
- The sway -> prospective-fall link is statistical/prospective, not demonstrated to be mechanistically caused
  by the specific delayed-feedback control loop modeled in Part 1/2 (vs. general frailty/sarcopenia acting as
  a shared upstream cause) — a concrete, testable, but unexecuted coupling is named (Part 4).
- Fernie 1982's institutionalized cohort shows a real but weak, non-dose-responsive sway-fall relationship —
  reported as a partial confirmation, not smoothed into the cleaner Maki 1994 result.

## Repro

The one new numeric result (Part 1) is reproduced exactly by:

```python
import math
mass_kg, h_m, g = 78.2, 0.969, 9.80665          # reused verbatim from this repo's postural_sway_pendulum_results.json
mgh = mass_kg * g * h_m
I = mass_kg * h_m**2
for frac in (0.91, 0.64):                        # Loram & Lakie 2002 (PMID 12482906); Casadio et al 2005 (PMID 15886131)
    omega = math.sqrt((mgh - frac*mgh) / I)
    print(frac, omega, math.log(2)/omega, math.log(0.2/0.005)/omega)
```

No git operations; no new files outside this doc pair. Every PMID above was checked live this session via
`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{esearch,esummary,efetch}.fcgi` (esearch to find candidates by
author/year/topic where the number was not already pinned, esummary/efetch to confirm bibliographic identity
and pull real abstract numbers) — none taken from memory.
