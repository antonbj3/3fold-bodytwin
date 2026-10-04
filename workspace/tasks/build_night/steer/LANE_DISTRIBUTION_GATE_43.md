# Styrning LANE_DISTRIBUTION_GATE_43 — efter r9

## The proof holds, but its terms disagree by a margin
`exact_class_proof = PASS_CONDITIONAL_26_INNER_KERNEL_TWO_MODE_BETA_LT1_ETA2_NE0`, noll motexempel
over 108 profile instances and 221 root attempts. But the condition is `β < 1`, and I calculated the margin:
`common_beta = 0,99999995` lies **5,0e−08** from 1, while `common_ratio_halfwidth_unknown = 1e−07`.
The case `[0,99999985; 1,00000005]` ** thus contains 1**. The proof applies in an area which is
narrower than the uncertainty in the parameter the condition is about.

And third mode is not a technicality: **12 by 12 third mode witnesses are above the requested
the error**. The two-mode assumption thus carries the proof, and every witness you tried outside of it breaks it.

## Hindret
Two things that look like precision but aren't. `independent_unhmixed_held_gap` is
`[0,28411910877330215; 0,28411910877330215]` — both endpoints equal, so a **point** printed
in casing form. And `initial_total_field_identity_exact = 0` at the same time as
`initial_total_field_held_gap = 0,148502787956807`: another case where identical summaries
hides different states. In addition, 113 retained failed starts and
`failed_verifier_cpu_s = UNKNOWN`.

## Changed operation
Enter `β` with the uncertainty it actually has, and determine if the condition `β < 1` can **be determined**
or just assumed. If the envelope crosses 1, the proof is not conditional but indefinite, and that is a difference.
Print what measurement would shrink the `β` envelope below 5e−08.

And make the third mode case the main track, not a witness. Twelve out of twelve above tolerance is one
result: state how big the error is in the quantity the decision consumes.

## Starkaste kontrollen
The two-mode solution as it stands, with β assumed strictly under 1. Any expansion must beat it on the
held-out cases, not on the synthetic ones.

## Falsifierare
If the `β` case cannot be shrunk below its own distance to 1,, then `PASS_CONDITIONAL` is effectively
`UNKNOWN` and shall be recorded as such. Say it with both numbers.
