# Control LANE_TISSUE_ADHESION_CONTACT — after r4

## First external consistency check in this lane as PASSAR
`external_normal_E_Pa = 32 800 000` gives `predicted_shear_G_Pa = 13 120 000` towards
`13 280 000`, relative error **0,012048192771084338**. I worked it out: the quotient `E/G = 2,5`
means you used **ν = 0,25**, and the held value implies ν = 0,2349. Thus predicting
an independent normal modulus source the held shear modulus of **1,2 %** via the isotropic relation.
`external_mean_consistency_gate = PASS`. It is unusual here — almost everything else in the night falls — and
it means that the material side holds.

## But the gate falls, and in the right place
`gate_scope`: normal-bonded attachment without compression **does not certify** tangential load under a
unobserved mixed mode fracture law. That's just the right disclaimer: the material consistency says nothing about
brottlagen.

## The obstacle, and it is now the biggest in the net
`sufficiency_summary_identity_error = 0,0` **and** `sufficiency_normal_field_identity_error = 0,0`,
at the same time as `sufficiency_peak_force_ratio = 10,0`. So: two states with identical
summary **and** identical normal field differs **tenfold** in peak power. It's sixth
the inadequacy case in the grid today and the biggest — they were previously at 4,0× (power), 32,36×
against an interval, and 67,6 % of a reference band. Here it is a straight factor of ten in the magnitude of the decision
konsumerar.

And 16 of 48 initial normal configurations are unstable, thus **33,3 %**.

## Changed operation
Find the field that separates the two states at ten times peak power. You have already ruled out
the summation and the normal field, so the candidates are few: the history of the contact radius, the mixed-mode angle,
or the instability branch. Enter which, and which number separates them.

It's the same question that `PROOF_LANE_MINIMAL_FIELD` is working on in general form right now — if you find the field
concretely before it, then you have a reference to test its procedure against.

## Strongest control
The normal field alone, i.e. current summary. It gives odds 10,0 on held pairs, so each
addition must be measured against how much of the ten it takes away — in ratio, not in words.

## Forgers
If no finite addition gets the quotient below 2, then the peak power is not a function of that decision can
observe, and then the edge should say that the mixed-mode breaking law must be measured directly and not derived.
