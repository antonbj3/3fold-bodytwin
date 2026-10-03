# Control LANE_FIBRIL_NETWORK_FRACTURE — round 1 (the coordinator, 1/10 15:30)

This is the only Sol lane (Anton: run Sol sparingly, the swarm takes the width). Put the budget where a weaker model is not enough: the net solver, the breaking criterion and the source criticism.

- **Obstacle (inherited, named):** three lanes have stayed that a summary toughness does not move between cases, and the needle friction eliminated the interface mechanisms one by one. What remains is the material: fiber recruitment and bond breakage.
- **Changed operation:** a bond-resolved mesh that delivers solidification, tearing energy and velocity dependence out of the SAME object.
- **Strongest control:** continuum hyperelastic model with the same information. The profit must be of a different kind, not just lower errors.
- **Forgers:** see the brief's four points. Change the supporting condition in the same round in the event of an error.

**Order is important.** Don't build skin first. Step 1 is Licup's collagen network: if the solver cannot handle the stress hardening, everything after that is pointless. Report step 1 separately even if the budget ends there — a verified net is a delivery.

**What's already there, read before building:** results/LANE_SURGICAL_BINDINGS (fibril inventories, Γ budgets per interface, pull-out energies), LANE_SKIN_TOUGHNESS_GAP (bridging candidate and poroelastic limit), LANE_NEEDLE_FRICTION_ANCHOR r1–r3 (the excluded mechanisms and the held speed points). DO NOT rebuild the fibril inventory; it is frozen with springs.

# Round 2 (the coordinator, 1/10 16:30)

## Review of round 1: step 1 passed in the central case
Verified by me in METRICS_R1.json: the mesh is bending dominated at insertion (bending part 0,917) and stiffens with slope **1,039** over a stress span of 5,6 with 12 points, and K/G0 reaches 12,5 while bending dominates. Licup's bending-dominated regime gives slope 1, so 1,039 is the right power law, not just the right direction. Both factor 2 gates PASS. It is a verified net solver, and it was the delivery for steps 1.

What doesn't hold: the small box falls on all caps (small_box_all_kappa_onset_rows FAIL), and the extended bend-dominated slope is UNKNOWN with only 2 points after insertion into the W20 geometry. The control (homogenized hyperelastic with the same information) gives TIE at 161 split states — expected, since step 1 is not that query. PORT.json has the correct form but all Γ are null.

## Changed operation: step 2, and separate the two interfaces
Lanen's own diagnosis is correct and should be followed: **breaking in a fibril is not the same interface as sliding between fibrils.** Do not mix them in one criterion.
1. **Separate the two.** Intrafibrillar: crosslinks inside the fibril break (parameters from Gautieri 2011 doi 10.1021/nl103943u, Depalle 2015 JMBBM 52:1–13). Interfibrillar: anchoring is lost and fibrils slide against each other (separate work per surface). Report which carries the fracture in which stress range.
2. **Validate the breach against a mesh reference BEFORE dermis.** Lane named Burla 2020 (mesh fracture in collagen) — check the locator and cite the figure. Same order as step 1: if the breaking criterion does not pass a published mesh, a dermis map is meaningless.
3. **Fix the small box or explain it.** A FAIL on all coat in small boxes is either finite size or an error in the criterion. Decide which with a size range; report the convergence as you did for the slot in the friction lane.
4. **Extend the bend-dominated span in W20** so the slope gets more than 2 points after the insertion, otherwise the anisotropy in step 3 cannot be trusted.

## Control, falsifiers, bans
- **Control:** same homogenized hyperelastic energy. The gain must be of a different kind: a failure criterion per interface that the continuum model does not have, not lower solidification errors.
- **Falsifier:** if the breaking criterion misses Burla's mesh reference by more than factor 2, or if the same criterion must be given different parameters for the two interfaces to hit, the split is wrong. Name the term and switch in the same round.
- **Forbidden:** to calibrate bond parameters against skin tear resistance or against 2,95. Step 3 (Annaidh anisotropy and velocity ratio) is not touched until the fracture passes a mesh reference.
- **Delivery:** fill Gamma_by_layer at least for dermis as range, or write exactly what measurement is missing to be able to do it.

# Round 3 (the coordinator, 1/10 18:00)

## Review of batch 2: a real structure find, and a correction of batch 1
**Finding, verified in METRICS_R2_V2.json:** the two interfaces are separated and the intersection is found. With the quotient between junction strength and fibril strength (eta), the fracture distribution becomes interfibrillar versus intrafibrillar: eta 0,25 → **56 versus 0**, eta 1 → **1 against 22**, eta 4 → **0 against 23**. The cover is thus close to eta ≈ 1, and which interface fails is completely controlled by that ratio. It is a solid result and it was the right division to make.

The span in the bending slope was extended as the control required: 18 points over stress span 15,6.

**Correction of my own accounting of round 1:** I wrote slope 1,039 as if step 1 was done. Run 2 shows that the slope varies between realizations (0,796, 0,880, 0,919 in different runs) and that `size_convergence` is **FAIL**: finite size dependence remains and is not resolved by larger W. Also, small boxes give ghost peaks (0,325 and 0,31 at L16/L24 vs. 0,513 in the matched one). So the power law is correct in magnitude, but the ensemble is not converged and 1,039 was a realization, not an ensemble value. Report the slope as a range of realizations starting now.

## Changed operation for round 3: convergence first, not Γ
The obstacle you named is correct but is shared now: **the data laws are laid out on the swarm** (BT-FIBX-FIBRIL-RUPTURE-WORK and BT-FIBX-FIBRIL-JUNCTION-WORK retrieve full fibril curve to break with hysteresis and separation work per node and node density in dermis). Do not download the same data yourself. Read their results.json if they exist; otherwise work with eta as parameter.

Sol's part is what the swarm can't handle:
1. **Solve the ensemble convergence.** Report peak stress and bending slope averaged with spread over at least 8 realizations per size, for at least three sizes. Determine if the ghost peaks at L16/L24 are finite size (scales off) or a criterion error (remains). That is the question that determines whether everything downstream matters.
2. **Close terminal separation.** `terminal_separation` is UNKNOWN and that is what prevents Γ. With eta as parameter: derive Γ(eta) as a curve instead of a number, with the active crack area measured and not assumed. Then the swarm's node work can be deployed as it comes, without a new round.
3. **Burla reference: do it unconditionally or tell me why it doesn't work.** `held_peak_central_factor2` is PASS_CONDITIONAL_ANALOGUE_TRANSFER and `held_peak_uniform_band` is UNKNOWN. Name exactly which transfer is conditional and what would make it unconditional.

## Control, falsifiers, bans
- **Control:** same homogenized hyperelastic energy. Γ(eta) as a curve is a guarantee the control cannot give; it is the profit, not lower error.
- **Falsifier:** if the ghost peaks do not scale with size, the breaking criterion is wrong, not the mesh. If Γ(eta) is not monotonic in eta, the coupling between the interfaces is incorrect. Name the term and switch in the same round.
- **Forbidden:** reporting a single realization as a net value (it was my mistake in round 1); to calibrate against the tear resistance of the skin or against 2,95.
- **Delivery:** PORT.json with Γ_dermis as FUNKTION of eta, plus the spread over realizations for each number.

## Addition during round 3 (the coordinator, 1/10 18:30) — what energy Γ should be compared with

Field's Web Search (~/research/HUD_NAL_HYPOTESER_20261001/HYPOTESER.md and KALLOR.md) Redefining the target. Read finds 1–4 where before you compare Γ The core, which I read in their file:

- **20–30 kJ/m² is a FAR FIELD J**, measured on tensile tested skin strips with 25 mm precut crack, mod III 20,4 and mod I 30,4 kJ/m² at 0,3 mm/s, juvenile pig. That number includes the entire mm-scale dissipation zone around the tip, not a local cutting energy.
- **Local and blade driven energies are 1–3 orders of magnitude lower:** intrinsic shear energy 276 ± 17 J/m² (bovine collagen membrane, EJ skin), J_IC 76–186 J/m² from needle penetration in pig liver in vivo, scissor cutting 1,77 kJ/m² (second hand quote, not opened).
- Consequence for our network: **Γ(eta) must be reported as two separate quantities**, not a number. (a) The local separation energy provided by the network when the breaking occurs at a bond or junction. (b) The far-field J that a sample with the same mesh would show, thus (a) plus the dissipation produced by the mesh itself in the zone around the tip (fiber straightening, rotation, sliding). The difference between (a) and (b) is a prediction, and that's the interesting one: if the mesh reproduces the ratio of local to far-field energy, about 10²–10³, without tuning there, it is a stronger result than hitting 20–30 kJ/m².
- **So never compare the local network Γ direkt mot 20–30 kJ/m².** There was an error in my previous briefing, and it's the same kind of error that Field found in the interlacing area.
- The honest limitation of the field, carries it further: even with 276 J/m² and DEJ 1,5–3 J/m² the quotient becomes DEJ versus dermis 0,005–0,01, well below the threshold of about 0,25 mentioned for the He–Hutchinson criterion. That threshold is VERIFIED in primary source. Lowering the toughness of the dermis is therefore not sufficient as an explanation.

## Parked 1/10 20:45, with a sharper target for round 4 (Field BLADBANA)
The lane is parked because the obstacle is DATA which the swarm picks up (BT-FIBX-FIBRIL-RUPTURE-WORK: full fibril curve to break with hysteresis; BT-FIBX-FIBRIL-JUNCTION-WORK: separation work per node and node density in the dermis). Read their results.json first on reboot.

**The target is now bracketed by an independent measurement.** Field SOL_FALT_BLADBANA_20261001 (PORT.json in romi_collab/build/) measured that the local capsule energy lies **70–117× under the far-field J of the skin**. I previously wrote 10²–10³ as expected ratio; Field range is at the lower end of that and is measured, not estimated. Samples and tissues differ (capsule vs. skin), so carry it as a direction of declared transfer, not as a reference. When the network provides local separation and far-field J, the ratio must be compared against 70–117 first, and against 10²–10³ as a further frame.

Two more things from the same delivery:
- **Blade steering does not remove the branches.** A cohesive model gives branches of about 1,6–1,75 µm per page despite fixed head tip, and the sharp 1 µm case is subcritical in the scenarios tested. So don't assume that a controlled trajectory will produce a single clean crack.
- **Peel is force per width and must be work balanced** before it becomes breaking energy. Same error as I did with the DEJ values.

# Round 4 (coordinator, 2/10 19:50) — resumed, and your step 3 target changed by factor 2,3

The lane was parked. It now takes over the slot from LANE_HIDDEN_AXIS_CERT, which closed after seven rounds.

## Why you're restarting: the target you were supposed to hit was wrong
Step 3 asked you to predict the rate ratio of needle force **2,95×** over 1→80 mm/s from fiber recruitment kinetics, with gate 1,475–5,9. The pin track has since shown that that number is largely a TIDSARTEFAKT, not a speed effect: the contact time per speed is 4 / 0,2 / 0,1 / 0,05 s, porcine skin relaxation is at 0,515–7,625 s, thus in the middle of the span, and the relaxation carries 28–76 % of the logarithm.

**The corrected target is thus a pure speed effect of 1,30–2,18×, corresponding to m = 0,0595–0,1782.** Do not fart against 2,95. And keep the interval as an interval — fiber recruitment that passes one end but not the other should be reported as conditional on the relaxation fraction, which then becomes the measurement that counts.

The same patch reopened two mechanisms in the pin track previously ruled out as "too weak" against the inflated target. Your recruitment mechanism was tested against the same incorrect threshold.

## What is resolved and what is still blocking
- Step 1 converged: the bending slope goes 1,186 → 1,033 → 0,980 with the size, and the two largest enclose 1,0. Don't run over it.
- Γ is still null, and the blocker is the two data laws. Both are now queued to the swarm as `BT-FIBX-FIBRIL-RUPTURE-WORK` and `BT-FIBX-FIBRIL-JUNCTION-WORK`. If their results are available: use them. If they are not: say it as `missing_prerequisite` and do step 3 anyway, since the speed ratio does not require an absolute Γ.
- Frozen microstructure numbers in `results/LANE_SURGICAL_BINDINGS/SOURCES_R2.json`: fibril diameter 82 ± 14 nm, fiber diameter 2,23 ± 0,96 µm, radius of curvature 6,56 µm, opening angle 32°. Do not remeasure them.
- **Elastic straightening and rotation is recoverable and must never be recorded as dissipation.** It was my mistake in a previous control and it was dropped.

## Control, falsifiers, bans
- **Control:** a continuum hyperelastic model with the same information. The gain should be a prediction of the velocity ratio from recruitment kinetics, not a calibrated parameter.
- **Falsifier:** if the recruitment falls outside 1,30–2,18× at both ends of the relaxation fraction, the mechanism is folded towards the corrected target. Name the term that carries the error and change a carrying condition in the same round.
- **Forbidden:** to fart against 2,95 or against 1,30–2,18; to treat the relaxation rate as a point value; to record elastic recovery as dissipation.

# Next round (the coordinator, 2/10 22:35) — DIN EXPONENT READING IS NOT ONE CERTIFIKAT
I have posted that the mesh bending slope 1,186 → 1,033 → 0,980 converges to the published 1,0, and
that the segment exponents 0,2353 / 0,3558 / 0,1877 s have something about regime. A graded find in old
the project traps that read, and I've verified the file:
`~/projects/bodytwin/docs/inherited_memory/exponent-cert-by-self-similar-collapse-and-causal-coupling-not-3pt-fit.md`
- A three-point log-log slope is **almost tautological**, and a prefactor that itself depends on
  the control parameter gives the correct slope with the WRONG mechanism.
- That a fitted exponent converges over several runs is **no certificate**.
- A stable exponent can be a stable FEL exponent: a quantity where was xmin-stable with drift 0,12
  and still had KS = 0,37, so bad absolute fit.

**Changed operation — two certificates instead of pussy:**
1. **Self-collapse.** If the quantity reads y(t;R) ~ Rⁿ·f(t/R^a) y/Rⁿ against t/R^a must fall on EN
   R-independent master curve. Measure peak normalized maximum deviation, and measure it ONLY in the adult range
   (the master curve above ~10 % of its peak) — small-t is often non-self and a relative measure
   explodes there. Use the high R end; self-similarity is asymptotic and the finite-R correction
   seen as operating in y_max/Rⁿ. Benchmark to beat: peak deviation 0,0068 over a factor 4 in R.
2. **Causal Coupling Button.** Scale the coupling of the mechanism and show that the exponent shifts by one
   PREDICTED amount. In our case: the angular distribution of fibril recruitment or the crosslink density
   is the button. A button that moves the exponent as predicted proves the mechanism; one that doesn't
   drops it.
3. Report **KS goodness** alongside each exponent, not just exponent stability.

**Check:** your current fit. **Falsifier:** if the collapse breaks while the slope persists, is
the exponent a prefactor effect and the mechanism interpretation falls — that is a perfectly good and important outcome.
**Prohibited:** reporting a converging fit as a certificate.

# Round 7 (the coordinator, 2/10 23:25) — the collapse fell per realization, so the KAUSALA KNAPPEN decide now
You ran the collapse and got the discrimination test exists for: mean deviates 4,34 % and
passes, the realizations deviate 18,84 % and fall. My posted "converges towards 1,0" was thus a
aggregate artifact, and the fix is in the log as mine. And you correctly rejected my imported
reference level — the tube transient's 0,0068 is not a universal fiber threshold, and I transferred a number
from another physical system.
- **Hinder:** with the collapse dropped per realization, only the second certificate remains, and that's it
  stronger: **the causal coupling button.** Scale the coupling of the mechanism and show that the exponent
  moves by a PREDICTED amount. Candidate buttons with us: fibril recruitment
  angular distribution, the crosslink density, and the network connectivity (you have 6 complete
  connectivity paths and 3 certified curves — use them).
- **Changed operation:** select EN button, predict exponent shift before run, run, and report
  predicted versus measured. A button that moves the exponent as predicted proves the mechanism even when
  the collapse is weak. One that does not fail the mechanism interpretation, and that is a perfectly good outcome.
- **And explain the dispersion instead of smoothing it:** why do the realizations diverge 18,84 %
  when the mean deviates 4,34 %? If the spread depends on the realization of the network rather than on
  the load, it is a structural characteristic of the network and should be reported as such, not as noise.
  Measure the correlation between deviation and connectivity per realization.
- **Report KS goodness** alongside each exponent, not just its stability. A stable exponent
  can be a stable error exponent.
- **Falsifier:** if no key within published range moves the exponent measurably, the exponent is
  not mechanism bearing and the three segment values 0,2353 / 0,3558 / 0,1877 should not be interpreted as regime.
**Prohibited:** reporting the average collapse as a collapse; to import a reference level from
another physical system without showing that it applies to fibers.

# New round (the coordinator, 3/10 00:25)
The collapse is saved: 95,25 % of the variation is within magnitude, and matched stress and magnitude yields 2,69 % versus 18,84 %. Rewrite the statement to be true — the collapse holds if conditions are matched, and the previous 18,84 % was an unsorted comparison. But the causal key yielded predicted_delta_m = 0,0, so no predicted shift, and the connectivity correlates weakly (−0,28), so the mechanism is not proven. Changed operation: select a button that ACTUALLY predicts a nonzero displacement — crosslink density is the candidate because it sits in the slip term. Predict the displacement before driving. Falsifiers: if no key within published ranges yields a predicted nonzero displacement, the exponent is not mechanism bearing and the three segment values should not be interpreted as regime.

## ADDENDUM — 62 literature records on bending versus stretching, verified
`source_repository/data/tissue_lit_refs/measurements.bending_vs_stretch_mode.jsonl`
(readable, read-only) carries **62 records** of exactly your lane's magnitude, with `unit` = *ratio
(bending/stretch at matched peak tensile strain)*. The verification distributes as **42
VERIFIED_PRIMARY_ABSTRACT, 10 COMPUTED_FROM_VERIFIED_INPUTS, 3 DERIVED_FROM_VERIFIED_PRIMARY and 3
VERIFIED_ABSENCE_OF_EVIDENCE** — the last class is rare and valuable: it says that someone has
searched and not found.

A floor is already calculated: `R_voidfloor_pointwise_linear_compression_inert` = **0,25** at e0 = 0,1 for
pure bending moment with the neutral axis in the middle of the wall. The tissues are named (collagen membrane, constructed
heart valve tissue, bovine late fibrocartilage, thoracic aorta), so the scope can be matched against
ours.

1. **Set the net's bending slope against these entries instead of against our own fixture.** The lane has one
   convergence to published 1,0 which holds only at matched stress and magnitude (2,69 %) and not per
   realization (18,84 %) — the 62 entries are the independent comparison it requires.
2. **Count the three VERIFIED_ABSENCE_OF_EVIDENCE entries separately.** If our open question coincides with
   one of them is the question unanswered in the literature and should be recorded as an acquisition item, not as our deficiency.
3. **And specify the scope match explicitly** per consumed item: tissue, strain level,
   geometry. A ratio measured on a heart valve must not be entered into a fibril model without the difference being noted.

**Prohibited:** Entries naming a device, device, or private target tissue analog — tissue generic
mechanics only.

## OBLIGATORISKT FILTER — read tasks/build_night/LIT_REFS_FILTER.md before the first entry
Two files in `tissue_lit_refs` are not used at all, and the other eight are read with row filters. After the filter
remaining 549 of 671 records. Your lane loses almost nothing: the thermal file has an affected entry of 25 and
bend/stretch file zero of 62. Report how many records the filter removed.
