# Steering LANE_NEEDLE_FRICTION_ANCHOR — round 1 (the coordinator, 1/10 13:15)

Anton 1/10: Sol runs sparingly, and the swarm (The swarm and reserve_worker) takes the width. This is the only Sol lane now. Therefore, put the budget where a weaker model will not do: derivation from first principles and source criticism, not cataloguing.

- **Hinder:** the effective toughness of the needle is a tool-specific lump, and friction and breakage are not separated by any independent measure.
- **Modified operation:** friction is derived from tissue stiffness and hole geometry (cavity recoil) and judged against measured retraction force or force in prepunched hole. The crime term will be the rest.
- **Strongest control:** INCISION R4's per-tool toughness (25 %) and Field breakdown with fitted friction (26 %).
- **Falsifier:** if derived friction deviates from measured by more than a factor 2 in the same tissue, or if the tool dependence of the residual term does not shrink, the partition is wrong. Say which term carries the error and change the operation in the same round.

Add follow-ups that the swarm can run (narrow calculations, data checks) in FOLLOWUPS.json with external_referent filled in: locator, compared_quantity, our_value and command/excerpt. The bonus now follows delivery.

## Coordination with Fields (1/10 13:35) — effective now
Field's SOL_FALT_NALFRIKTION_20261001 owns the data card: retraction in the same hole and insertion in prepunched hole, digitized and source coated, in ../3fold-motion-engine/_private/romi_collab/build/SOL_FALT_NALFRIKTION_20261001/raw/anchor_curves/ with INDEX.json. Build NOT your own data hunt for the same curves. Read their INDEX.json; catalog is not yet available, do the derivation and pre-register the comparison against their curves. The line's part is the derivation: friction per length out of skin stiffness + hole geometry (cavity return) and the fracture as a residue. Deliver PORT.json with friction_per_length(stiffness, r_needle, r_hole) and residual_fracture_energy to Field.

## Field anchor curves now exist (1/10 13:45) — and their value
INDEX.json supplied: khadem2016 fig7, langevin2001 fig3 (BI/NO/UNI), zoudani2025 fig3b + fig4c (pigskin). Reading errors in raw/digitization.json, sources in raw/sources.json. Field verdict: assumption of LIKA friction on reinsertion is UNSAFE — no series separates forward friction in first vs. repeat pass in same hole, and pullout force requires handling of direction and contact position change.

Impact on the lane: compare NOT mot en "same-hole"-series that does not exist. Do this instead:
1. Derive friction per length from skin stiffness + Hole geometry (cavity return) and judge against the curves as EXISTS, using Field's reading error as uncertainty. Specify which curve measures what (forward, extended, suction cup fixed) before the comparison.
2. If the derivation can only be tested against tensile force: model the change of direction explicitly (contact position, possible adhesion) instead of assuming symmetry. It is then part of the derivation, not a source of error.
3. Report what would be required to make the assumption decisive (Chu 2019 full protocol according to Field). Say it as a measurement specification, not as an obstacle.

# Round 2 (the coordinator, 1/10 13:45) — replaces the points above as missions

## Review of round 1: carrying, and a failed assumption
Approved as real design change. Crucial: Field speed fit (trained on total max forces 16G/21G) missed the held-out shaft force after penetration (Barnett 2016 Fig11, 25G, pigskin) by factor 3,6–9,9, at least 2,96 after 0,05 N read margin; 4/4 speeds failed factor-2-the gate. I recalculated the quotas myself from r1/CONTROL_COMPONENT_V1.json: 9,88 / 4,96 / 3,94 / 3,57.

Stronger than the error size is RIKTNINGEN: measured shaft force rises 2,95× from 1 to 80 mm/s, while the cavity image with fixed contact geometry, constant shear and pressure out of decreasing modulus yields almost flat (1,07×). Disruption, not calibration error. Fixed contact + velocity-independent shear is therefore pigskin in the range — the lane's own assumption failed, which is the right kind of result.

Disclaimer: the held-out anchor is FIELD's digitization (SOL_FALT_NALKRAFT_20261001/raw/rate_pilot.json, sha256 posted), not our own reading, and it's PENDING_INDEPENDENT_REVIEW on theirs. Treat as semi-external until someone reads Fig11 independently. Type it in RESULTS.md.

## Changed operation: change the interface law, not the parameters
1. **Velocity-dependent interface attraction from first principles.** The shear stress must come from a mechanism that CAN provides increasing force with velocity: viscous film/grinterstitial fluid lubrication, poroelastic pore pressure build-up in the contact zone (fluid does not escape 80 mm/s), or velocity-dependent adhesion/stick-slip. Derive which gives 2,95× over 1→80 mm/s with the skin's own transport parameters (permeability, viscosity, contact length) and which cannot. Discriminating: several mechanisms provide rising force but different EXPONENT and different dependence on needle diameter.
2. **Finite strain and slot geometry.** The hole is a slot, not a circular cavity; the printed circular cavity (fixed b/R, c, xi) gave 53,2 % error against Montanari bulk slopes. Calculate the contact pressure for a slot at finite strain and state how the contact area differs. Montanari's prehole calibration is INDATA, never reference.
3. **First then the fracture residue.** With an interface layer that can handle velocity dependence: remove friction and deformation and see if the rest becomes tool common across INCISION R4's held-out tool (where the error was 22–55 %).

## Control, reference, falsifier
- Control: Field frozen fits (now refuted at component level), INCISION R4's per-tool toughness, Montanari's published cavity calibration.
- Reference: same held-out Barnett shaft force plus at least one new independent from Field's anchor_curves/INDEX.json. State per curve what it measures before comparison.
- Falsifier: if the derived mechanism does not give 2,95× ±factor 2 with published skin parameters, it is the wrong mechanism. Name the term that carries the error and try the next mechanism in SAME round. Don't end on a negative outcome.
- Forbidden: fitting a speed exponent against the held-out points. The exponent should fall out of the mechanism.

## Deliverable
PORTS_R1.json has the correct form but all null. Fill native_skin_friction_prediction_N, and tool_shared_Gamma_J_m2 if step 3 is reached. Narrow follow-ups to the swarm in FOLLOWUPS.json with external_referent complete (locator, compared_quantity, our_value, command/excerpt).

# Round 3 (the coordinator, 1/10 14:45) — replaces the mission above

## Review of round 2: bearing and discriminating. Two mechanisms excluded.
The round did exactly what the round steer 2 asked for: the exponent, not the amplitude, decided. Verified by me in r2/mechanismsv1/ANALYSIS.json and ROUND_SUMMARY.json (gate 1,475–5,9 around the held-out 2,95):
- **Viscous film / Brinkman / Couette: UTESLUTEN.** All film thicknesses and both permeability scales give velocity ratio exactly 80,0, i.e. linear in velocity, against measured 2,95. Furthermore, the force level 54–162× is too high. Linear viscous shear is the wrong mechanism, regardless of film thickness — that's the strength of trying the exponent.
- **Poroelastic pore pressure: UTESLUTEN in the central case.** Ratio 1,03–1,38 (optimistic casing up to 2,24) to gate bottom 1,475. Too weak, not too strong.
- **Binding kinetics (Schallamach/Bell): KVAR, but conditional.** 58 of 105 grid instances are consistent with both the strength increase and the quotient. No forward prediction is possible: skin-metal's N0, lambda, stiffness and binding-/lossningsdata are missing.
- **Free slot contact at finite strain: numerically clear.** Convergence went from 0,161 (old fixed area) to 0,0028 at mesh refinement. So the old fixed incompressible contact was not converged, and that is fixed.
- **New independent anchor, and it drops another one:** Langevin 2001 (doi 10.1152/jappl.2001.91.6.2471) gives UNI/NO-kvot 2,67, while a history-free contact law predicts 1. Thus, a history-free law cannot explain the direction dependence. Note that the anchor is total skin plus subcutaneous pullout, not pure shaft friction; write that disclaimer.

## Changed operation: make bond law measured, not conditional
The obstacle is named and right: the binding mechanism survives only as a grid of possibilities. Round 3 should make it a prediction.
1. **Get an independent pressure dependence binding-/lossningslag for skin against metal.** Search published measurements of adhesion, shear strength per area or friction against steel/stainless steel for skin or skin imitation, with pressure dependence. Accept related systems (gelatin, PDMS, agar, imitation epidermis) but declare interspecies and state why the transfer is reasonable. Freeze N0, lambda and stiffness with locator. If no such law exists for skin: say it exactly, and use the closest one with declared crossing — build NOT further on a grid.
2. **Measure the reference geometry under load.** The true loaded slot or anisotropic law is missing. Derive what the contact fraction 0,715 and P_GR 7,38 will be at the loaded geometry instead of the unloaded one.
3. **The directional dependence as its own test.** With a pressure-dependent binding law with sign: predict Langevin's UNI/NO 2,67 without farting against it. It is now the lane's sharpest reference, as it is a quota and thus insensitive to calibration level.

## Control, counterfeiters, bans
- **Control:** Field frozen total force fit, INCISION R4's per-tool toughness, Montanaris cavity calibration, and known Schallamach/Bell-frenew with same sub-information (no algorithm gain claimed).
- **Falsifier:** if the frozen bond law gives velocity ratio outside 1,475–5,9 or UNI/NO outside 1,8–3,8, the bond image is refuted. Then name what carries the error and try the next mechanism in the same round. Don't end on a negative outcome.
- **Forbidden:** to pussy lambda or N0 against 2,95 or 2,67. Both must fall out.
- **Delivery:** PORTS is still null. Fill native_skin_friction_prediction_N this round, also as intervals with declared uncertainty.
# Round 4 (the coordinator, 2/10 16:10) — replaces the mission in round 3

## CORRECTION FIRST (the coordinator, 2/10 16:30): my own exclusions were too strongly written
I wrote below that three mechanisms are excluded. I checked against `r2/mechanisms_v1/ANALYSIS.json` and it is true for ONE out of three. Correct position, with the gate [1,475; 5,9] around the held-out 2,95:
- **Viscous film: UTESLUTEN, by a large margin.** Ratio exact 80,0 at all film thicknesses and power level 54–162× too high. It's an order mistake, not a calibration mistake.
- **Poroelastic pore pressure: UTESLUTEN ONLY I CENTRALFALLET.** `transport_optimistic_pore_ratio_minmax` = [1,033; 2,241], and 2,241 are INSIDE the gate. Thus, the optimistic envelope is not excluded, and an excluded central case is not an excluded mechanism. What would determine that is a permeability-FREE upper limit, namely the undrained/drained modulus ratio — it would have determined the exclusion without requiring the permeability, and it is not counted.
- **Binding Kinetics: ALDRIG UTESLUTEN.** 58 of 105 grid instances are consistent with both the strength gain and the quotient. The specific FRYSTA PASSED law gives 1,455 against the gate 1,475, thus a rejection by 1,4 % margin — it is within every reasonable uncertainty and must not be booked as an exclusion.

That I wrote "uteslutna" about all three is the same error class that the lane itself should avoid: a central case and a cover are not the same statement. Treat the below with that fix.

## Review of round 3 and the only question left
R1–R3 excluded viscous film, constrained the poroelasticity to its central case and left the binding kinetics conditional, and the directional difference failed on elastic memory. What remains is not another mechanism to test against the same number. It's that the **number itself is an average that hides a cover**, and that's my own finding on lane's data, not lane's:

The held-out rate dependence 2,95× over 1→80 mm/s corresponds to a dimensionless exponent m = ln(2,95)/ln(80) = 0,2469. But our own segment exponents over the same interval are **0,235, 0,356 and 0,188** — a dispersion of factor 1,9. One thus, single exponent cannot describe the interval, and the one fitted to 2,95 will lie between two regimes and belong to neither of them. Compare: viscous film gives m = 1, poroelastic m = 0,023–0,074, the frozen bonding law m = 0,086. Neither is close to 0,2469 — but now that's a weaker objection than 0,2469 not being a constant.

## CORRECTION of my own wording, written before I counted (2/10 16:15)
I recalculated the segments myself from `r1/CONTROL_COMPONENT_V1.json` (speed_mm_s [1, 20, 40, 80], held_phase3_N [0,0816; 0,1652; 0,2114; 0,2408]) and got 0,2353 / 0,3558 / 0,1877, 1,896, total m = 0,24687 and quotient dispersion exactly 2,95. The numbers are correct. **But the pattern is not what I first wrote.** m STIGER and FALLER then. A cover between two regimes gives a monotone m, so the two-regime picture cannot produce this. Therefore, do not seek a cover rate as the first hypothesis. What can produce a maximum in the middle is (a) read error — three segments out of four points, and an error at the midpoint moves two segments in the opposite direction at the same time, which is exactly the signature we see; or (b) a mechanism with its own internal maximum in speed. Try (a) first, because it's cheaper and more likely.

## Changed operation: first determine if the dispersion is read error
1. **Determine if the dispersion survives the read error.** Propagate Field read error from `raw/digitization.json` point by point to each segment exponent and specify confidence intervals for 0,2353, 0,3558 and 0,1877. Test in particular the hypothesis that the pattern points to: that ONE read error at the midpoint (20 or 40 mm/s) explains the entire deviation, since an incorrect center point moves the two surrounding segments in opposite directions. Calculate how large that error must be to make m constant, and compare with the declared reading error. If it is less: the dispersion is an artifact, m is a constant within the error, and IT is the round's result — report it and move on to point 3 in the same round.
2. **Only if the dispersion survives: a mechanism with inner maximum, not two regimes.** A monotonic cover is already excluded by the shape of the pattern. Name a mechanism whose m has a maximum in speed, and say which quantity sets the position of the maximum. Each regime should get its own mechanism with its own scaling, and the wrap rate should drop out of the comparison between their time scales — not set. Name for each regime the quantity that sets the time scale, and calculate at what rate the two are equal. That number is the prediction.
3. **Discriminating, and that's the point of dimensionless:** two mechanisms that happen to give the same m in a range are distinguished by how m depends on NEEDLE DIAMETER and on tissue permeability. Enter the predicted derivative of m with respect to the diameter for each candidate, as it is a measurement found in the literature for several diameters.

## Control, counterfeiters, bans
- **Control:** a single fitted power law over the entire range, so exactly m = 0,2469. It has three free numbers against our three segments and will look good; the gain must be that the wrap rate is PREDICTED and that the diameter dependence points correctly, not that the residuals become smaller.
- **Falsifier:** if the segment spread is within the reading error, there is no wrap — report it as an honest negative and move on to the diameter dependency in the same round. If the wrap rate falls outside the tried range 1–80 mm/s, the dual-mode image is not identified by this anchor; tell me what speed would identify it.
- **Forbidden:** to fit an exponent against 2,95; to report a single m as the lane result; to choose the mechanism according to which produces the least residual.

## Deliverable
PORTS is still null. Supply the wrapping speed with intervals, the m of the two regimes with their mechanisms, and dm/d(diameter) per candidate. If step 1 falls: deliver the segment exponents with propagated read error and the statement that the interval does not separate them.

# ADDENDUM to round 4 (the coordinator, 2/10 16:30) — a published in vivo measure of MOTSATT TECKEN

This preempts the segment exponents if you only manage one. A crossing of today's agent waves found the only hard external contradiction to a BodyTwin result found in the entire material.

**The measure:** needle-tissue friction shear stress FAILED from 0,58 ± 0,27 to 0,16 ± 0,08 kPa as the insertion speed increased from 0,2 to 10 mm/s, in vivo in rat brain. Regional: cortex 0,227, external capsule 0,222, caudate-putamen 0,383 kPa. Locator according to the agent: PMID 25151066, doi 10.1016/j.jneumeth.2014.08.012. **Verify the locator yourself before using the numbers** — it comes from an agent's reading, not our own, and in the same material only 14 of 43 source records were page verified. If it doesn't match: say so and report what you searched for.

Our number rises 2,95× over 1→80 mm/s. Two margin reports disagreeing on THE SIGN is exactly what `margin_net` and `claim_federation` exist to decide.

**OThe operation:** determines if this is a CONTRADICTION or a REGIME LIMIT, and count it, don't describe it. The validity boxes look disjoint on at least three axes and each one must be declared explicitly:
1. **Species and tissue:** rat brain versus pig skin. Two completely different nets and water levels.
2. **Magnitude:** shear stress per area against total axis force. It is not the same quantity, and a force can rise while the voltage falls if the contact area grows faster than the voltage falls. **Calculate if that very thing can explain both signs at the same time** — it's the sharpest hypothesis in the entire supplement and it can be tested with our own slot geometry.
3. **Hastighetsintervall:** 0,2–10 mot 1–80 mm/s. They just overlap between 1 and 10 mm/s, and our first segment, 1→20 mm/s med m = 0,2353, is right there.

**The outcomes, and both are worth delivering.** If the boxes are disjoint, it is a regime boundary and BOTH outcomes stand — and it will be the first time our 2,95× is placed in a declared regime instead of just asserted. If they overlap, it is a contradiction, and then `hidden_variable` should name the axis that carries it, with the common sign gate (pinned commit `352c6d3`), not with the axis-parallel sweep.

**Forbidden:** to dismiss the measure because it is rat; species difference is a declared axis, not a reason to throw away a measure. To merge kPa and N to a comparison without the area. To report "border of regime" without the declared boxes.

# Round 5 (the coordinator, 2/10 16:45)

## Review of round 4: you refuted my hypothesis, and it was right
Steg 1 were made exactly as ordered and the outcome went towards the client. `constant_slope_feasible = True` med `minimum_common_error_N = 0,005933`, which is 7,3 / 3,6 / 2,8 / 2,5 % of the four forces held, while Field declared replicate scattering is 15–27 %So the required error lies well in the uncertainty, `crossover_identified = False`; and **m = 0,2469 is thus the correct description**. Min "1,9× dissemination conceals a cover" was an artifact of treating digitized points as accurate. Don't look further for a cover in this anchor, and don't build anything on the segment differences.

Two things you supplied are more important than they look: the conditional diameter derivatives dm/dd and the amplitude-free double ratio 0,933 / 0,789. Both are independent of the calibration level, and that's exactly what is needed now.

## The mission: the addition about it MOTSATTA THE SIGN, which r4 never saw
The addendum was written 16:30, five minutes after your round ended, so it's untouched. **Read it in this file and make it the main task of the round.** Brief: a published in vivo measurement in rat brain has the frictional stress FALLANDE 0,58 ± 0,27 → 0,16 ± 0,08 kPa over 0,2 → 10 mm/s, against ours 2,95× rising above 1 → 80 mm/s. Verify the locator yourself (stated PMID 25151066, doi 10.1016/j.jneumeth.2014.08.012).

In addition, now there is the right instrument to determine it, which did not exist when the addendum was written:
1. **Area is the sharpest hypothesis and you can now test it.** A KRAFT can rise while a VOLTAGE falls if the contact area grows faster than the voltage falls. You have a converged free slot contact at finite strain. Calculate the contact area as a function of speed in our case, and determine if the two signs can be the same physics expressed in two quantities. If yes, the contradiction is apparent and that is the strongest possible outcome.
2. **Use the double ratio, not the amplitude, for the comparison.** It is amplitude-free and thus survives that the two measurements have different calibration, species and geometry.
3. **Declare the boxes and count.** `margin_net` and `claim_federation` decide REGIME LIMIT against CONTRADICT, and the shared character gate is on pinned commit `352c6d3`, branch `research/typed-throws-20261001`. Disjoint boxes ⇒ regime boundary and both stand, and our 2,95× is placed for the first time in a declared regime. Overlapping ⇒ contradiction, and `hidden_variable` should name the axis with the common the gate, never with the axis-parallel sweep.

## Control, falsifier, prohibitions
- **Control:** to treat the two measurements as incompatible without counting area or boxes, i.e. choosing a side. That's what we did before.
- **Falsifierare:** if the area cannot explain the character difference AND the boxes overlap, are at least one of the two measurements wrong or measure something other than it says. Name which and what would determine it.
- **Forbidden:** to dismiss the rat measurement because it is a rat — species difference is a declared axis, not a reason to throw away a measurement. To compare kPa against N without the area. To report regime limit without the declared boxes.

# Round 6 (the coordinator, 2/10 17:35) — turn the instrument towards our OWN number

## Review of round 5: PASS, and you refuted my hypothesis by several orders of magnitude
The regime limit is determined and calculated with the pinned engine, not described. Our 2,95× is now in a declared regime. My area hypothesis was wrong by a margin: the area changes 0,0313 % over 80× rate and the ceiling is 1,4142×, so it can't flip a sign. Your explanation is stronger than mine because it is measured in the same preparation: the rat protocol has **120 s dwell** at 0,2 mm/s, and deposit against post-relaxation separates **8,143×** in the same brain.

## The quest: is our own 2,95× partly a relaxation artifact?
It follows directly from your own finding and it attacks our speech, not theirs. If a state after relaxation separates 8,14× from insertion in the same tissue, **the time under load** is a bearing axis — and then we need to know where our own anchor lies on it.

1. **Determine our own protocol's time under load.** At 1 mm/s a given insertion distance takes 80× longer than at 80 mm/s. Calculate the actual contact time per speed for Barnett's geometry. If the tissue relaxes on a time scale that lies INOM that interval, the four points measure different degrees of relaxation and not just different rates.
2. **Calculate how much of 2,95× can be relaxation.** With a measured skin relaxation time scale and the calculated contact time per point: how much of the force difference between 1 and 80 mm/s is due to the fast point simply not having time to relax? Enter it as a factor with intervals. If it becomes comparable to 2,95, our velocity dependence is largely a time artifact, and that is a more important result than anything we found in five rounds.
3. **Discriminatory test:** a pure speed effect and a relaxation effect are distinguished by the fact that the latter depends on the INSERTION DEPTH, since the time under load grows with the distance at fixed speed. Give predicted depth dependence for both and say which Barnett's curves can differ.
4. **Find the relaxation time scale independently**, not from our own curves. If it is not available for skin: say so and enter the measurement.

## Control, counterfeiters, prohibition
- **Control:** our own interpretation of 2,95× as a pure speed dependence, so m = 0,2469. The gain should be a breakdown of speed and time-under-load, not a better exponent.
- **Falsifier:** if the skin relaxation time scale is outside the contact time range of all four points, relaxation cannot contribute, and 2,95× is purely velocity dependent. Report it — it STRENGTHENS our speech and is a perfectly good outcome.
- **Forbidden:** to estimate the relaxation time scale from the curves it is supposed to explain; mixing dwell and deposit time in the same term without saying which is which.

# Round 7 (the coordinator, 2/10 18:20) — try if the exclusions survive the corrected target

## Round Review 6: you attacked our own speech, and it didn't hold
The contact times 4 / 0,2 / 0,1 / 0,05 s against porcent relaxation span 0,515–7,625 s: the relaxation is inside the range, so the points measure different degree of relaxation. The proportion is 28–76 % of the logarithm, thus factor 1,351–2,273, and the remainder becomes 1,30–2,18× instead of 2,95×. Proper reticence not to attribute when contact age is unknown.

## The mission, and it's a consequence you didn't count on
I recalculated the dimensionless requirement with your relaxation rate: **m falls from 0,2469 to 0,0595–0,1782.** Towards that interval two of the mechanisms lane EXCLUDED in round 2–3 lie within:
- poroelastic pore pressure in its upper envelope, m = 0,074 — inside
- the frozen transferred bond law, m = 0,086 — inside
- viscous film, m = 1,0 — still 5,6× too high, remains excluded

The exclusions were thus tested against a target inflated by a time artifact. **Retry them against the corrected range**, one at a time, and tell each one if it survives, using the same gate logic as before but with the new threshold. It is not an overtaking: the threshold is a different quantity now.

1. **The poroelastic first**, because it has a permeabilityFREE upper limit that was never counted: the ratio undrained/drained modulus. It determines the exclusion without requiring the permeability, and it stood in play 3:s notes as "andra starka" without being executed. Make it now.
2. **The bond law**: its m = 0,086 lies inside, but it was conditional on 58 of 105 grid instances. With the new threshold: how many instances are compatible, and does the grid shrink to a prediction?
3. **Keep the relaxation part as an interval all the way through.** You have 1,351–2,273, not a number. Each exclusion should be tested against BOTH ends, and a mechanism that survives only at one end should be reported as conditional on the relaxation fraction — which then becomes the deciding measurement.

## Control, counterfeiters, prohibition
- **Control:** round 2–3:s own exclusions against m = 0,2469. The payoff is to say which of them were artifacts of an inflated target.
- **Falsifier:** if both mechanisms also fall toward the corrected range, the exclusions remain and the relaxation finding only changes the number, not the conclusions. Report it straight.
- **Forbidden:** to treat the relaxation fraction as a point value; to revive viscous film; to fit the relaxation part so that a favorite mechanism fits.

# Round 8 (the coordinator, 2/10 18:50) — separate the two candidates, stop widening the gate

## Review of round 7: PASS, and your own obstacle formulate the next task correctly
The reopening is confirmed: the poroelastic overlaps the target in 1,2978–1,3839, the frozen bond layer lies inside at 1,4552, and bond plus contact age spans 2,0144–3,5369, thus the observed 2,95. The permeability free ceiling is calculated and I verified the half-space relation for all three ν_d (1,25 × (1−0,2²)/(1−0,5²) = 1,6000 exactly).

But you write yourself that wide gate authorization does not identify the physics, and that is correct. **Two mechanisms that both fit the interval is a discrimination task, not an outcome.** Don't widen the gate any more.

## The mission: use the amplitude-free instruments you already built
You have two quantities that are independent of calibration level and can therefore distinguish the candidates without us knowing the absolute level:
1. **The diameter derivative dm/dd from round 4** (micro −0,0609 / −0,0202 / −0,0141; nano −0,2354 / −0,0630 / …). Poroelastic pore pressure and bond kinetics have different diameter dependences: one is set by drainage length versus contact width, the other by bond density per area. **Predict dm/dd for each candidate separately and compare with the measured ones.** It is a ratio to a ratio, thus independent of the calibration level.
2. **The double ratio 0,9328 / 0,7886 from round 5.** Enter what each candidate predicts for it.

Report per candidate: predicted dm/dd, predicted double ratio, and whether the measured is within the respective prediction. A candidate that passes the range but misses the diameter dependence is refuted — and it's the kind of refutation that's worth something, because it doesn't go via the amplitude.

## If both can handle this too
Say it plainly and then state **the only measurement** that would distinguish them, with quantity, unit and preparation. After eight rounds, a named tiebreaker is a better result than a ninth candidate.

## Control, counterfeiters, bans
- **Control:** the gate authorization from round 7, so that both fit the interval. The gain should be to distinguish them or to name the measurement that does so.
- **Falsifier:** if the diameter dependence cannot be predicted for any of the candidates without additional unknowns, the instrument is not discriminating and should be stated rather than forced.
- **Forbidden:** to widen the residual interval to fit a candidate; to compare absolute levels when an amplitude-free ratio exists; to introduce a third candidate before the two are separated.

# Round 9 (the coordinator, 2/10 19:15) — AVSLUTANDE: write down the measurement and gate, nothing more

## Round Review 8: You guys were right and I was wrong about the instrument
`measured_diameter_derivative`, `measured_Q`, and `measured_Cread` are all None, and that's correct. I named the diameter derivative and double quotient "instrument du redan byggt" — but the r4 field is named `conditional_25G_segment_dm_dd_per_mm`, with *conditional* in the name, and I read predictions as measurements. The fault was mine, not the round's.

**But you delivered something better than a discrimination: a pure separation.** `Cread` 1,0 for pure pressure against 1,4552 for the frozen bond layer is 45,5 %, against Field's replicate spread 15–27 %. This means that ONE measurement determines eight rounds of investigation.

This is the lane's **last round**. Don't budget for a ninth candidate or a wider gate.

## Deliver three items and close
1. **The measurement specification for Cread**, written so that it can be purchased or searched: exactly what quantity, in what units, on what preparation, at what geometry and speed, with what replication is required to distinguish 1,0 from 1,455 given 15–27 % dispersion. Calculate the number of replicates explicitly — with 45,5 % separation and 27 % dispersion, it's a small number, and specifying it makes the difference between a wish and an order.
2. **PORT.json to Field**, which consumes our Γ per layer. It shall bear: that native friction and Γ remain UNKNOWN; that the velocity dependence 2,95× to 28–76 % may be contact age and that the pure velocity effect is thus 1,30–2,18×; that two mechanisms remain and which; the permeability-free roof E*u/E*d 1,2–1,6; and the Cread measurement as the deciding factor. Each number with its distribution type named.
3. **EOne line about what is EXCLUDED and survives the correction:** viscous film, 5,6× too high even against the corrected target. It's the only clean cut after eight rounds, and it should stand as such.

## Forbidden
To open a new mechanism. To report UNKNOWN without the measurement that would make it KNOWN. To hide that two candidates remain.


## PIN UPDATED (the coordinator 2/10 21:05)
Pinned commit moved from `cf9c1e9` to `352c6d3`, three commits later, seven files changed. The reason is that the three contain precisely the corrections we have pursued tonight: *Abstain when the downstream cert says nothing, instead of reading silence as a pass*, *Stop a validity filter from passing a point it never checked*, and *Say whether the declaration could separate two discordant reports, and count unconsumed verdicts*. The older pin would have let silence pass as approved.
The engine is now also on both cloud hosts under `/opt/agents/graph_engine`, with the cited path symlinked there, so a job no longer needs to reimplement it. Verified: import through the cited path works on both.
