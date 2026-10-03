# LANE_ZERO_FLUX_AUDIT — a flow that doesn't die when the drive dies

Resultatmapp `results/LANE_ZERO_FLUX_AUDIT/`.

## The error that motivates the lane, measured last night
Cellane's round 3 found why a cell exceeded the molecular oxygen ceiling at **78 of 78**
physical points: the native ATP synthase produces **56 882 ATP per second at ZERO proton-motive
force**. A flux that does not go to zero when its driving force goes to zero makes work out of nothing, and
that's why the ceiling could be broken with median 2,18× and max 12,8×.

Same error in the adenine nucleotide translocator: **12 of 15 equilibrium points wrong natively, 0 after repair**,
where the repaired shape satisfies a **published zero flow equation to 9,77e-17**. And
synthase equilibrium failed in **9 of 9** points.

**It is a FELKLASS and not a bug.** Any constitutive relationship in the project that writes a flow
as a function of a state can have the same defect, and it can be tested mechanically without new physics:
set the thrust to zero and see if the flow does.

## Why it deserves its own lane
The 43 cells carry **1 336 parameter families and 155 state families** over 1 491 magnitude families
(unit census, same night). A fraction of them are flows, and of these only those with one
identifiable driving force to which the test applies. But each hit is a number that is wrong in a direction that
not seen in any unit check and not in any dimensional analysis — only against thermodynamics.
The unit invariance error in the oxygen case was exactly 0,0, so dimensional checks don't catch it.

## Do like this
1. **Find the flows.** Go through `tasks/free48/sources/*/` and list each expression that returns a
   flow, a speed or a current. For each: what is the magnitude of its DRIVKRAFT? Candidates are
   differences — concentration difference, pressure difference, potential difference, temperature difference,
   chemical affinity. A flow whose driving force cannot be named is already a flag and should
   reported as such.
2. **Run the zero flow test.** Set the thrust to zero and evaluate. The flow must be exactly zero
   or zero within machine precision. Report per flow: operating force name, flow value at zero
   driving force, and if it is zero within 1e-12 relative to its typical value.
3. **And run the REVERSE test, which is stronger.** Reverse the sign of the driving force. The flow should reverse sign.
   A flux that is positive in both directions breaks microscopic reversibility though
   happens to be zero at zero. Report the number that pass the null test but fail the sign test —
   that number is the actual find of the lane, since the null test is the easy one.
4. **Measure the distribution before any system statement is made.** How many flows are there, how many have one
   named driving force, how many pass the null test, how many pass the sign test. Four numbers. A
   individual violation is not a statement about the system.
5. **For each violation: is it a miswritten relation or a deliberate approximation?** A linearized
   relation around a working point SHOULD not be valid at zero driving force, and then the error is not i
   the relationship but in the fact that its scope is not declared. Separate the two. The first group
   are bugs, the other are undeclared scopes, and they require different action.
6. **The repair pattern already exists.** Cellane's repaired form meets the published one
   the zero flow equation to 9,77e-17. Read how it is written before suggesting your own form, and
   use the same pattern where appropriate. Their external form is a published equation 14 for one
   zero flow surface.

## Strongest control and falsifier
- **Check:** unit and dimension checks, which all 43 cells already pass. The point is that
  they don't capture this — the unit invariance error in the oxygen case was exactly 0,0. The profit is the number of crimes
  that the unit checks let through.
- **Falsifier:** if no flow outside the already known cell breaks the zero flow test, the error is
  isolated and not an error class, and the lane closes with it. It is a perfectly good and soothing outcome
  and should be reported outright, not downplayed.
- **Prohibited:** fixing a relationship without specifying whether it was a bug or an undeclared one
  linearization; reporting a violation without the driving force's name and unit; changing a cell's
  default values so that a previous result moves quietly — declare the repair next to and
  measure what it moves.

## Delivery
`PORT.json`: the four distribution numbers, one row per flow with driving force, zero flow value and
character test result, as well as the breakdown of bug versus undeclared scope. Plus, for each offense,
what past results rest on it.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.
