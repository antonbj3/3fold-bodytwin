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

# New round (the coordinator, 3/10 02:05)
The error class is confirmed and the distribution is measured: 32 flows; 29 klarar nolltestet, 26 teckentestet,
and 2 pass zero but fall on signs — Just what I predicted would be the actual find.
Starkast av allt: 647 av 2 000 instanser med negativ dissipation nativt mot 0 after repair, i.e.
ett andra-lagsbrott i 32,4 % of cases, and a closed system that produced 2,95 mikrojoule ur
ingenting.

Changed operation, in that order:
1. De 11 unsolved families of 43 Each unsolved family can carry the same error and the distribution is
   Not complete until they're there. — is not driven, the flow is not
   Isolable, or isn't the cell running?
2. De 6 which fall within the character test, each shall be classified as: BUGG eller ODEKLARERAD LINARISATION.
   linjarisering kring en arbetspunkt SKA does not apply in the case of reverse thrust, and then the error is that:
   the validity range is not declared, not that the relationship is wrong. The two groups require different
   action and shall not be mixed.
3. Track the consequences. The repair moves net glucose in the brain with 19 percent. For each
   ‘repaired flow' means the past results of which rest on it; ореpaired the speech? That's the most expensive part
   and the person who makes the audit a correction instead of an observation.

Strongest control: the unit and dimensional controls as all 43 cells are already passing through.
that they let through all these errors.

Falsifierare: om de 11 unsolved families turn out to lack identifiable flows completely, is 32 Total flows
population and distribution is complete — Then say it, it makes the statement stronger and not weaker.

# New Round (koordinatorn, 3/10 04:05)
You measured the third independent instance of the night by the same principle: lumped model 22,90 procent fel mot
axially dissolved 0,126, that is 181 together with the heat dose 254 times and
the reversal of the deflection verdict is now overdetermined. I started a dedicated lane on this,
LANE_COARSEST_SUFFICIENT, so leave the general question there and stick to your own.

Your question is the wrong class, and it has four numbers: 32 flows; 32 med namngiven drivkraft, 29 klarar
nolltestet, 27 teckentestet, med 5 real character error after fixing another inactive one.
remains is the most expensive and most important: track the consequences. A repair moved a net uptake
70 percent and another net glucose in the brain 19 per cent. For each repaired flow: which earlier
results rest on the unrepaired number, and do they move?

And the classification bug against undeclared linearization is still not delivered for the 5. En
linjarisering kring en arbetspunkt SKA does not apply in the case of reverse thrust; then the error is that:
the field of validity is not declared; the two groups require different measures and shall not be mixed.

Falsifier: if no previous result moves by any repair, the error class is real but without
Consistency, and then it's a cleaning issue and not a bargain.

# New round (coordinator, 3/10 04:35)
The classification is done and the answer is that all five are real: two bugs, one topology, two unknown physics,
and zero Taylor linearizations, nothing can be explained as undeclared validity.
impact detection indicates that the cost is concentrated: renal outcomes do not move at all;
gastriskt 3,6 percent, blood brain barrier 9 till 48 % local and up to 70 i netto.

Changed operation, in that order:
1The two unknown physics cases are the most expensive and should be described so that someone else can work on them.
   Each: what quantity, what momentum, and what would a proper relationship have to fulfill?
   cases marked unknown physics without that description are not transferable.
2Repair the blood brain barrier first, because the cost is sitting there. And report which previous ones
   results based on the unrepaired number — a move on 70 procent i netto kan ha burit
   slutsatser.
3. The gastric historical anchor is not verified and its screening falls. An anchor that does not
   is verified must not carry a 3,6-percentage adjustment. Verify or mark the outcome as
   villkorligt.

Falsifier: if the two unknown-physics cases on closer examination prove to be the same error in two
shape, there are four faults and not five — and it is a correction worth making of your own account.

# New Round (koordinatorn, 3/10 05:20)
You verified the gastric anchor and first dropped a miscitation. — the pair used was ten-percent
dextrose and not water, and the verified values are: 15 and 70 minuter. Bra. Men sedan blir resultatet
Inconvenient and that is now the most important question of the lan.

Den NATIVA the model, which violates thermodynamics, lies 19 % low on liquid and 13 percent high
It's stuck. REPARERADE the mixing model lies 170 percentage of fluid failure and at least 157 on a solid. Thus:
The repair makes the compliance dramatically worse.

Changed operation, and determines NOT This is by choosing the appropriate model: two interpretations are open and
they require different measures.
One, compensation: the thermodynamic error of the native model compensated another error in the opposite direction; and
then there is a third error in the model that is now visible only when the first is gone. The test is to look for the second
The fault, not to reverse the repair.
Two, parametric: repair is correct in principle but its parameters are wrong, and then it must:
the repaired form is calibrated to the verified anchor and see if it KAN reach 15 and 70 inom rimliga
parametrar.
Decide which, and tell me how you decided it. A repair that impairs compliance is not
automatic error, but it must not remain unexplained.

And the attachment of solid content falls from 0,984 till 0,684 In mixing. There is a big change in
a quantity that should be closely preserved. Tell us if it is a consequence of the repair or a fault of its own.

Falsifier: if the repaired form cannot reach 15 and 70 for ANY parameter set within published
interval, repair is wrong in shape and not only in parameters — And that's the heaviest outcome.

# OPERATOR DIRECTIV 3/10 09:25: EMOVE THE SOLUTION OF HELA BANAN
Operator's directive: Raise the resolution all over the track. You have the hardest measured case — klumpad modell
22,90 % error against axially dissolved 0,126 percent, that is 181 times, over 205 instanser.

Changed operation: run the repaired flows at the UPLOST setting and report what is
move. You have already shown that the repair of a flow moves net removal 70 percent and
net blood glucose in the brain 19; the question now is how much of what is repair and how much is
The resolution. The two have been mixed in each number you reported, and none of them can be attributed until
they are separated. Report both effects separately for at least three flows.

Falsifier: if the dissolution carries most of the 70 percent, flow repairs are less
significant than they look and the error class is more expensive in accounting than in physics.
