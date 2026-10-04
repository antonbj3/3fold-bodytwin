# Steering LANE_PRECISION_LABEL — new lane 2026-10-04

## The finding came from the swarm and is verified
A creative job on edge H-E29 took the number `2.220261437908497` from our own network and reconstructed it
as exactly `3397/1530`. I verified: `Fraction(x).limit_denominator(10**4)` gives 3397/1530 and
`float(3397/1530) == x` is true. The number therefore carries FOUR significant digits and the remaining eleven are
floating-point representation. Quoting it with fifteen digits is false precision.

It is a class, not a case. A chain that divides measured quantities produces long decimals from
short measurements, and a long decimal is read as precise. The same error occurs in my own work today twice:
a residual reported as −0,000000 when the inputs carried 3e-5, and an interval quoted with
ten digits that was 1,81 float32-ULP wide.

## My started instrument and why it is inadequate
`tasks/assembly/precision_audit.py` exists and RUNS, but the filter is wrong in two steps that are both instructive:
- The first version caught 24 026 numbers, almost all 1/6, 4/3 and 8/3 printed as long floating-point numbers.
  A simple fraction in long format is not false precision, it is the opposite problem.
- The second version required denominators ≥ 50 and caught 1/86400 and 1/60000 instead, meaning seconds per
  day and minutes — unit conversions, not measurements. It also reads its own output file, which is
  a feedback loop I have not closed.

## The task
Build the discriminator that separates three things, because that is the whole question:
1. **Two short measurements divided** — dangerous, must be labelled. Signature: both numerator and denominator are
   credible readings, three to five digits, neither a round unit factor.
2. **Simple fraction or unit conversion** — harmless. 1/6, 4/3, 1/86400, 1/60.
3. **Actual long calculation or irrational number** — harmless, leave it.

Exclude output files under `results/ASSEMBLY_PRECISION_AUDIT` and `ASSEMBLY_QUANTITY_INDEX` from the input.

## Strongest control
The swarm's own case: 2.220261437908497 should end up in class 1. And at least ten numbers you select from the internet
that you first classify by hand — if the discriminator does not match your manual classification on them, it is not
ready.

## Falsifier
Print before the run how many of the network's numbers you expect in class 1. If the discriminator catches more
than a fifth of all numbers it is too broad and catches unit conversions; if it catches zero it is too
narrow and misses the swarm's own case.

## Why it is worth a round
Every number we publish with more digits than it carries invites the next reader to compare on digits that
do not exist. It is the same error as the convention collisions, one level down: wrong frame replaced by wrong
resolution. A precision label per number makes it impossible to repeat.
