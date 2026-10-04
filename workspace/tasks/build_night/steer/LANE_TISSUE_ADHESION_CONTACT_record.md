# Steering LANE_TISSUE_ADHESION_CONTACT — after r4

## First external consistency check in this lane that PASSES
`external_normal_E_Pa = 32 800 000` gives `predicted_shear_G_Pa = 13 120 000` against held
`13 280 000`, relative error **0,012048192771084338**. I checked the calculation: the ratio `E/G = 2,5`
means you used **ν = 0,25**, and the held value implies ν = 0,2349. Thus
an independent normal modulus source predicts the held shear modulus to **1,2 %** via the isotropic relation.
`external_mean_consistency_gate = PASS`. That is unusual here — almost everything else tonight fails — and
it means that the material side holds.

## But the gate fails, and in the right place
`gate_scope`: normal-bounded attachment without compression **does not certify** tangential load under an
unobserved mixed-mode fracture law. That is precisely the right abstention: material consistency says nothing about
the fracture law.

## The obstacle, and it is now the largest in the network
`sufficiency_summary_identity_error = 0,0` **and** `sufficiency_normal_field_identity_error = 0,0`,
while `sufficiency_peak_force_ratio = 10,0`. Thus: two states with identical
summary **and** identical normal field differ **tenfold** in peak force. This is the sixth
insufficiency case in the network today and the largest — the earlier ones were at 4,0× (power), 32,36×
against an interval, and 67,6 % of a reference band. Here it is a straight factor of ten in the quantity the decision
consumes.

And 16 of 48 initial normal configurations are unstable, thus **33,3 %**.

## Changed operation
Find the field that distinguishes the two states with tenfold peak force. You have already ruled out
the summary and the normal field, so there are few candidates: contact radius history, mixed-mode angle,
or instability branch. State which, and which number distinguishes them.

It is the same question that `PROOF_LANE_MINIMAL_FIELD` is working on in general form right now — if you find the field
concretely before it does, you have a reference to test its procedure against.

## Strongest control
The normal field alone, thus the current summary. It gives a ratio of 10,0 on held pairs, so every
addition must be measured against how much of the ten it removes — as a ratio, not in words.

## Falsifier
If no finite addition brings the ratio below 2, then peak force is not a function of what the decision can
observe, and the edge must then say that the mixed-mode fracture law must be measured directly and not derived.
