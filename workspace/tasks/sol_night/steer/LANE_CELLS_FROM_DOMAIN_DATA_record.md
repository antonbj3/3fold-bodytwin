# Styrning LANE_CELLS_FROM_DOMAIN_DATA — efter r48

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
