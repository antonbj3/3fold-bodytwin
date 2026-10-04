# Control LANE_CELLS_FROM_DOMAIN_DATA — after r48

## Calculate in MOLEKYLER, not in µM — then the obstacle is immediately visible
`target_gap_uM = 0,20756740959193906` is bit-identical in r41, r43, r45, r46, r47 and r48, and in r45
also says 0,4151348191838781 = exactly double. I divided: the number is **number of molecules ×
0,1037837048 µM**, i.e. 2,0000 molecules in r48 and 4,0000 in r45, in a declared volume of
1,6000e-17 L = 0,0160 µm³. The gap is not a measured concentration but a unit conversion of one
integers.

And then it appears that six rounds have hidden: `budget_uM = 0,0915193951833` is **0,8818 molecules**.
The budget thus requires resolution finer than EN molecule in the declared volume. It is impossible of
construction, not because the instrument is bad. Rewrite the budget and gap in molecules in each
artifact hereafter; in µM they are read as measured quantities and they are not.

## The obstacle, reformulated
There is nothing to retrieve that makes 0,88 molecules measurable. Either the budget is increased to at least one
molecule (and then the question is whether a decision that tolerates the uncertainty of a molecule is useful), or so
a larger volume is declared so that a molecule becomes a smaller concentration. Both are choices, aren't they
measurements, and they must be reported as choices.

## The second in r48 to be said outright
`modeled_ATP_ligands = 0`, while `ADP_ligands = 17` and `ANP_ligands = 4` (non-hydrolyzable
analogue). The chain derives `net_native_ATP_consumption = [2, 3, 4]` from structures that do not contain
a single modeled ATP. ADP bound and analog bound states frame the catalysis cycle but do not show
ATP, and it should be written in plain text next to the consumption number — otherwise read 21 nucleotide ligands as
support for a ATP occupation that no structure carries.

## Strongest control
The same 2 266 heteroresidy controls, but with the ATP occupation set to zero everywhere. If
the consumption figure does not change, it is not structurally borne, and then the structural acquisition is not what it is
carry the result.

## The falsifier
If the resolution of one molecule is sufficient for the decision, write which decision can still be made
at ±1 molecule. If no decision can be made there, the entire chain is below its own resolution limit and that
is RESULTATET, not a deficit to fix.


## Addition after r49 — you already have a decision, expressed in molecules
I recalculated r49's three radii in molecular units (0,1037837048 µM per molecule):

| quantity | µM | molecules |
|---|---|---|
| budget | 0,0915193952 | **0,88183** |
| singleton radius | 0,0206063822 | 0,19855 |
| unknown-bound radius | 0,1243900870 | **1,19855** |
| target gap | 0,2075674096 | 2,00000 |

The unknown-bound radius is the singleton radius plus **exactly 1,00000 molecule**, and the target gap is exactly
2 molecules. Everything in the lane is quantized at one molecule; The µM suit hides it.

And then the decision is already in the numbers: the budget 0,882 molecules lies BETWEEN the two radii.
The chain passes its budget when the species is singleton (0,199 < 0,882) and misses it with **1,359×** so
soon the bond is unknown (1,199 > 0,882). Formulate it as the result of the lane: *determines if the species is
single bound, and the chain keeps its budget; leave the bond unknown, and it cannot hold it.*
It's an ability with a falsifier, not a deficit.

## The acquisition exchange must be recorded
`thesis_pages = 117` gave `thesis_nonmissing_bound_data = 4` and `thesis_calibration_amounts = 4`,
thus 0,034 usable numbers per page. Compare with the structure acquisition in r48 (4 structures → 21 ligands
with 583 heavy atoms). Print the yield per source type in the next round, so that the next acquisition goes where the numbers are
actually exists.
