# PROOF_LANE_SUFFICIENT_BY_CONSTRUCTION — a summary that CANNOT be insufficient

## The pattern, six cases measured by me today, all with identity error exactly zero
Each row is a pair of states that our summary **does not distinguish** and that the decision distinguishes:

| kedja | vad som var identiskt | dolt gap | mot beslutets egen skala |
|---|---|---|---|
| adhesive contact | summary **and** normal field, both 0,0 | **factor 10,0** in peak force | a straight tenfold in the quantity the decision consumes |
| immune window | both identity errors 0,0 | 17,60797843701942 nM MAC | **32,36×** the interval [0; 0,5441770553588867] |
| the disc's electrical alias | `electrical_summary_identity_error_exact = 0` | 0,08108314647592979 MPa | **67,6 %** of the ground-truth band, **and the diagnosis changes** |
| diffusivity | identical mean-V and mean-I | 20,0 W | **4,0×** the 5,0 W the terminals give |
| the biofilm plateau | every observed point | 3/16 in normalised signal | the entire decision interval |
| the distribution gate | `initial_total_field_identity_exact = 0` | 0,148502787956807 | — |

Six independent chains, six times the same error. It is not six accidents.

## The question I want answered, and it is not "what is missing"
`PROOF_LANE_MINIMAL_FIELD` runs in parallel to find the smallest addition **per case**. It repairs one
row of the table at a time. This question is the other one: is there a **construction** that makes the repair
unnecessary?

My suspicion, and I want to know whether it holds: every row above is a summary that is **invariant
under a transformation the decision is not invariant under**. Diffusivity is the cleanest case —
the mean of a product is not the product of the means, so the summary quotients out
the time correlation while the decision reads it. If this holds generally, the elegant
summary is the **complete invariant for the decision's own equivalence relation**: quotient out
exactly what the decision cannot see, neither more nor less.

So: construct it, or show that it does not exist. Theorems, not reasoning. In particular:

- existence: for which class of decision functions does a complete invariant exist that is **finite**,
  and what breaks finiteness?
- **canonicity**: is it unique up to bijection, or are there several non-equivalent complete
  invariants for the same decision?
- **cost**: how many numbers must it carry compared with the raw state, and is there a case where it
  is no smaller — thus where no compression is possible?
- **composition**: if two decisions read the same summary, is the complete invariant for
  the pair the union of their invariants or something strictly larger?
- and the **negative** case: a decision function whose complete invariant is provably
  the raw state itself, so that every summary is insufficient. I want to recognise that
  situation before I build on it.

## The test
Your construction must be run against all six rows above and for each row say either which field it
would have required from the start, or that the row belongs to the negative case. **A construction that does not
reproduce at least four of six is not an answer.**

## The strongest control
The raw state as summary. It is always sufficient, never compressing, and it is exactly
what we do today when we give up. Your construction is measured in how many numbers it saves against it, per row.

## Falsifier
Construct a decision function and a state pair yourself where your complete invariant is equal but
the decisions differ. If you succeed, the construction is not complete and must be rejected, not patched.

## Warning
A width-based rounding rule fails — your own `PROOF_LANE_PRECISION_PROP` made one fail on
all 1 200 constructed boundary cases. The endpoint rule applies to all precision you state.

## Regler
`PENDING_INDEPENDENT_REVIEW`, ingen utsaga om biologisk validering, inget excluded_category material, inget som
namnger en enskild person.
