# PROOF_LANE_MINIMAL_FIELD — which SINGLE field makes an insufficient summary sufficient?

## The pattern, five cases measured today by me
Five times today, a summary with **identity error exactly 0** has hidden states giving different
decisions. I have calculated every pair myself:

| var | identiska sammanfattningar | dolt gap | gapet mot beslutets egen skala |
|---|---|---|---|
| biofilm plateau | ordinates 0, ¼, ½, ¾ at 0/10/20/30 s | 3/16 = 0,1875 in normalized signal | the decision’s entire interval |
| immune window | both identity errors 0,0 | 17,60797843701942 nM MAC | **32,36×** the interval [0; 0,5441770553588867] |
| distribution gate | `initial_total_field_identity_exact = 0` | 0,148502787956807 | — |
| disc’s electrical alias | `electrical_summary_identity_error_exact = 0` | **0,08108314647592979 MPa** | 67,6 % of the ground truth band 0,53–0,65, **and the diagnosis changes** |
| diffusivity | identical mean V and mean I | **20,0 W** | 4,0× the 5,0 W implied by the terminals |

You have already given me the two-copy contract: `d` factors through `σ` exactly when a formula is
unsatisfiable. It **decides** whether a summary suffices. It does not say what is missing.

## The question
Given a summary `σ`, a decision function `d` and a witness pair that `σ` does not distinguish — is there
a **constructive** procedure returning a minimal additional field `f` such that `(σ, f)`
suffices? And is the minimal field **unique** or only minimal in size?

Theorems, not reasoning, and a measurable quantity per condition. Especially:

- when `f` can be chosen as a **function of what `σ` already carries** (i.e. no new measurement required, only
  another projection) and when it requires a new observable;
- a lower bound: how many fields are required **at minimum** for a given family of witness pairs, so "add
  a field" does not become an infinite ladder;
- whether order matters — whether `f₁` then `f₂` may give a different final result than `f₂` then `f₁`;
- and the negative case: when **no** finite addition suffices, with a criterion I can run.

## The strongest control
Adding the entire raw state, which always suffices and is never useful. Your minimal field
should be measured against it in number of numbers that must be carried, not elegance.

## Falsifier
Construct a case yourself where your procedure returns a field that does not suffice, or where a
smaller field suffices. If you succeed, the criterion is wrong and should be rejected, not patched.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no biological validation statement, no excluded_category material, nothing like
names an individual. The endpoint rule applies to any precision.
