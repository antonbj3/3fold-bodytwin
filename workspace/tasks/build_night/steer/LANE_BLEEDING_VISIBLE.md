# Steering LANE_BLEEDING_VISIBLE

# Round 2 (the coordinator, 2/10 22:50) — fix the optics FIRST, then calculate the reregistrations
The falsifier did not trigger: the field is obscured, 146,15–317,31 mm² in the first minute at 2 mm depth and 30°
tilt, against a critical film thickness of 50,87 µm at 550 nm. The chain is thus justified.
- **Obstacle 1, what makes the number unusable:** `held_roggan_optics = FAIL`, 6 of 10 anchor consistencies,
  maximum relative error 0,304, `compiled_independent = False`. The hydrostatics are consistent to
  3,3e-13, so the arithmetic is not the error — the optical coupling is. Show WHY it failed: is the
  wavelength dependence, hematocrit dependence, or scattering versus absorption modeled incorrectly?
  Report which of the ten anchor consistencies fails and by what number.
- **Obstacle 2, the decision-relevant number was missing:** `clinical_extra_reregistrations = None`. A
  declared policy range 0–4 is not a calculated value. This is the whole reason the lane exists: a
  blood pool that looks real changes no decision, an obscured field fraction that forces a 150-second
  reregistration changes the robot branch's calculation.
- **Changed operation:** calculate reregistrations per procedure as a function of occluded area per
  minute and the landmark's area. A landmark is hidden when occluded area covers it, and
  registration must be repeated. At 150 s per registration (124 s solution + 21 s resampling) and
  7,25 s per update, 1192 updates fit in a procedure of published mean duration — say how
  many of them are lost per extra reregistration.
- **Falsifier:** if the number of extra reregistrations is zero even at the high area 317 mm²/min,
  bleeding does not bind the robot and the chain should not be built further.
**Forbidden:** rendering; reporting the area as a result while the optics gate is FAIL.

# Round 3 (the coordinator, 2/10 23:20) — the new finding about timing comes BEFORE the optics
You calculated the decision chain conditionally: 1 to 3 extra reregistrations, 21 to 62 lost updates
of 1 192, and marker pose error 1,07e-14 mm so the geometry is clean. But you also found something that weighs
more than your own lane: **`published_timing_excludes_spatial_registration = True`**.
- **Take that first.** The whole robot branch rests on registration taking 150 s, of which 124 s solution and
  21 s resampling, and on 1 192 updates fitting in a procedure. If the published timing
  excludes spatial registration, the 150 s are scoped incorrectly, and then the binding term is incorrectly
  set in four lanes. Say what the timing ACTUALLY includes, with a quote and locator, and what it
  excludes. Recalculate 150 s and 1 192 if the scope changes.
- **Then the optics.** `optical_gate = FAIL` remains, and your two repairs went the wrong way: removing
  packaging raised maximum relative error from 0,304 to 0,718. Both pass 6 of 10. Show WHICH
  of the ten anchor consistencies fails and by what number, instead of switching variants. Is it
  wavelength dependence, hematocrit dependence, or scattering versus absorption?
- **Falsifier:** if no optical variant within published coefficients passes more than 6 of 10, the
  occluded area is not certifiable from the optical data we have, and then the decision number is conditional
  forever — say it plainly and leave it as an acquisition item.
**Forbidden:** rendering; reporting the area as certified while the gate is FAIL.

# Round 4 (the coordinator, 3/10 00:00) — YOUR OPTICS GATE HAS A SOLUTION THAT ALREADY EXISTS
You have stayed at `optical_gate = FAIL` with 6 of 10 anchor consistencies for two rounds, and your own
repair attempts went the wrong way (removing packaging raised the maximum error 0,304 → 0,718).
**The laser lane measured tonight that a HISTORICAL optical set passes 9 of 10.** It therefore
already exists in the project.
- **Do this first:** read `results/LANE_LASER_SURGERY/` — r1.json, RESULTS.md and
  `SKIN_OPTICAL_ANCHORS_R1_V1.json` with 27 held anchors across three tissue compartments. Identify
  the historical set that passes 9, run your gate with it, and report the outcome against your
  current 6.
- **Two warnings from there that apply directly to you:** anisotropy g = 0,9 is an ASSUMPTION, not a
  measurement, and `BALLISTIC_IDENTIFIABILITY_R1_V1.json` shows **380× difference in unscattered signal across
  g sensitivity** — so your occluded area may depend on an assumed number. And **bulk dermis must not
  be copied to papillary and reticular layers**, which is exactly what `SURG_HEMOSTASIS_model.py`
  does with its STANDARD values per layer. Say which of your layers have their own support and which
  inherit bulk.
- **The absorption length 46,5 µm** is a number you can use directly against your critical film thickness
  50,87 µm — they are almost equal, which means a film that just obscures is also about one
  absorption length thick. That is no coincidence and it can simplify the whole chain. Calculate it.
- **Falsifier:** if the historical set passes 9 of 10 in the laser lane but fewer than 6 with
  you, the gates differ and then the gates should be compared, not the coefficients.

# New round (the coordinator, 3/10 00:35)
You went from 6 to 9 of 10 and replaced the threshold with photon transport — it is now a chain. Changed operation: take the marker return all the way to the decision. At what film thickness does the return fall below what the tracking system requires to maintain registration, and how many extra reregistrations does that give per procedure? You have 0,484-0,666 at 50 um and 0,068-0,210 at 200 um, so the threshold lies between them. And solve the ONLY remaining optical error: 650 nm oxygenated with 25,3 percent error. Falsifier: if the marker return never falls below the system requirement at the occluded areas you calculated, bleeding does not bind the robot and the chain should close with that.

# New round (the coordinator, 3/10 01:20)
DROP ISOTROPY. The gate's own name says where the error sits: FAIL_ISOTROPIC_ALL10_PRESERVED. Fixing
650 nm oxygenated to 19,3 % breaks 450 nm deoxygenated, and THREE independent routes give −23,0 %,
−20,2 % and −21,2 %. Three different constructions landing within two percentage points of one another is
the signature of missing physics, not a fit deficiency — a fit deficiency would have spread out.

Changed operation: a red blood cell is a biconcave disc, not a sphere, and it orients itself in
shear flow. Shape consistency already passes 10 of 10 at both 0 and 30 degrees for the external
fixture, so the geometry is clean — it is the ISOTROPIC averaging over orientation that must
go. Calculate scattering with an orientation distribution instead of an isotropic mean, and say whether
the 450 deoxy error shrinks. The laser lane measured 380× difference in unscattered signal across anisotropy sensitivity
and that g = 0,9 is an assumption — that number is your strongest reason to suspect that particular axis.

Strongest control: the historical set that passes 9 of 10 isotropically. The gain is if an
orientation-resolved variant passes 10 of 10 without any other anchor falling.

Falsifier: if an orientation distribution does not reduce the 450 deoxy error below 10 %, the missing
physics is not anisotropy, and then the next candidate must be named — dependent scattering at high hematocrit
is the natural one, and you already have a dependent-dense error of −30 % pointing there.

# New round (the coordinator, 3/10 02:05)
You refuted my hypothesis with an argument stronger than a failed test: a distribution-independent
molecular absorption bound excludes orientation-based repair, and 0 of 48 evaluations pass
the ten-percent requirement. The ceiling sits in the molecule, not the geometry. I was wrong about anisotropy.

The next candidate is the one I named in the falsifier and it now has support in your own numbers: dependent
scattering at high hematocrit. Your dense phase gave the anisotropy parameter 0,0196 against held 0,028, thus
−30 percent, and that is the direction in which independent scattering would be wrong. At physiological hematocrit
the cells are not optically isolated from one another, and the independent-scattering assumption is then invalid.

Changed operation: calculate scattering with a structure factor for dense packing instead of as the sum of
independent scatterers. The gate is now called FAIL_ORIENTATION_ONLY_STRAIGHT_RAY_10PERCENT, so the
straight-ray assumption is also part of the error — say whether the structure factor and multiple scattering are separable
or whether they must be treated together.

Strongest control: the historical set that passes 9 of 10 under independent scattering.
Falsifier: if a structure factor does not reduce the 450-deoxy error below 10 percent either, the
optical coefficients are not certifiable from available data and the occluded area remains
conditional. Then say it plainly as an acquisition item and move on to the marker return, which is already a
usable number.

# New round (the coordinator, 3/10 04:20)
You have drifted into library forensics. Three hypothesis gates fail on eigenvalue ratios,
centroid deviation at 0,05 pixel and an external vision library's contour numerics, and
paper_exact_source_version is unknown. That is not your question.

Your question is whether bleeding obscures the field so a robot must reregister, and you already have what
is needed: independent optical validation against MCML with 8 million photons and 38 of 38 angles,
marker return 0,48 to 0,67 at 50 micrometers of film and 0,07 to 0,21 at 200, occluded area 146
to 317 square millimeters in the first minute, and a reregistration that costs 21 of 1 192 updates.

Changed operation, and it is the last question the lane needs to answer: at what film thickness does
the marker return fall below what tracking requires, and how many extra reregistrations per procedure does that give?
The threshold lies between 50 and 200 micrometers because the return goes from 0,6 to 0,14 between them.
Deliver a number with its interval and close the chain.

Leave the contour operator. If it is needed to know what tracking requires, instead state the requirement as
a return threshold with justification and declare it as an assumption — a declared assumption with a
number is worth more than an unresolved library comparison.

Falsifier: if the marker return never falls below a reasonable tracking threshold at the areas you calculated,
bleeding does not bind the robot, and then the lane closes with that answer. That is a valid outcome.

# Round 15 (the coordinator, 3/10 07:20) — you retrieve sources instead of calculating the decision number
This round retrieved a primary document with fifteen blood-related and thirteen
blood-occlusion rows, and then calculated nothing: new_transport_solves = 0, new_photons = 0,
optical_model_reestimated = False. And the three decision numbers remain None for the fifth round in
a row — clinical_extra_reregistrations, actual_sensor_hcrit_um, new_certified_occult_area.

You have had everything you need for three rounds: external optical validation against MCML with 8 million photons
and 38 of 38 angles, an external gate that passes within 7 percent with a control that is 31
percent wrong, marker return 0,48 to 0,67 at 50 micrometers of film and 0,07 to 0,21 at 200,
occluded area 146 to 317 square millimeters in the first minute, and a reregistration that costs 21 of
1 192 updates.

Changed operation, and it is the ONLY task: calculate the sensor's threshold thickness. At what
film thickness does the marker return fall below what tracking requires? The threshold lies between 50 and 200
micrometers because the return goes from about 0,6 to 0,14 between them. If you do not know what tracking
requires, DECLARE a requirement with justification and calculate with it — a declared assumption with a number
is worth more than a fifth None. Then state the number of extra reregistrations per procedure.

Forbidden: retrieving another document before the three numbers are calculated. The sources are not
the bottleneck.

Falsifier: if the marker return never falls below a reasonable tracking threshold at the areas you calculated,
bleeding does not bind the robot and the lane CLOSES with that answer. That is a valid outcome and it is
better than a sixth round without a decision number.
