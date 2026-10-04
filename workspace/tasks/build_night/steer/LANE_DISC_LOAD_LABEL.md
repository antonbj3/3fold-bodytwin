# Steering LANE_DISC_LOAD_LABEL — after r9

## You closed the alias, and that was what I asked for
The previous steering required the electrical alias to be closed first, because everything computed on top of
a summary that does not separate states inherits the gap. That is done:
`primary_electrical_mechanical_mean_gap_MPa = 0,0006816163471549052` and max
**0,0017345916495278926**. Against the previous round's alias gap of 0,08108314647592979 MPa, it is **46,7×**
better on max and **119×** on mean, and the max gap is now **1,445 %** of the reference band. Zero fitted
parameters.

## But a new and larger gap has taken its place
`full_salt_curve_ion_summary_identity_error_exact = 0` while
`matrix_pressure_gap_MPa = **0,11155306633210749**` across **99** matrix cases, of which **92** are above
2 kPa — that is, 92,9 %. The gap is **93,0 %** of the entire reference band, and it sits behind an
ion summary that is **exactly identical**. It is the eighth insufficiency case in the net and
the second largest in absolute terms.

The changed-volume case is milder: 0,032555646878362245 MPa across 861 cases.

## The obstacle
`changed_state_pressure_enclosure_MPa = [0,16663797048690213; 0,24982681048629915]`, width 0,083189.
The reference band 0,53–0,65 lies **entirely outside** that enclosure. So the new state is not what Wilke
measured, and `exact_primary_Ringers_composition = NOT_ACQUIRED` says why: the bath composition has
not been acquired, and it sets the ionic strength on which the entire matrix pressure depends.

And `primary_supports_are_statistical_or_physical = False` — the supports are declared assumptions, not
measurements. That must appear in every row that cites the enclosure.

## Changed operation
Acquire the bath's exact composition. It is a single entry from the source, and it determines whether the matrix gap of
0,1116 MPa is physics or a free parameter. Everything else you compute before then inherits the same uncertainty
as the alias did.

## Strongest control
The electrical-mechanical agreement you just established, 0,0017 MPa on max. Every new
matrix result must be measured against it, not against the reference — because the reference lies outside the enclosure and therefore cannot
separate two matrix readings.

## Falsifier
If the matrix gap of 0,1116 MPa does not shrink when the bath composition is acquired, then
ionic strength does not carry it, and the edge must say that matrix pressure is undetermined by what we can measure.
