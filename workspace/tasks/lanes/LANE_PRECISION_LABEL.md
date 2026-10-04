# Control LANE_PRECISION_LABEL — new lane 2026-10-04

## The find came out of the swarm and is verified
A creative job on edge H-E29 took the number `2.220261437908497` from our own network and reconstructed it
as exactly `3397/1530`. I verified: `Fraction(x).limit_denominator(10**4)` gives 3397/1530 and
`float(3397/1530) == x` is true. The number thus carries FYRA value digits and the remaining eleven are
floating point representation. Quoting it in fifteen figures is false precision.

It is a class, not a case. A chain that divides measured quantities produces long decimals out of them
short measurements, and a long decimal is read as exact. The same error is found in my own work on day two
times: a residual reported as −0,000000 when the inputs carried 3e-5, and an interval quoted with
ten digits that were 1,81 float32-ULP wide.

## My started instrument and why it doesn't work
`tasks/assembly/precision_audit.py` exists and RUNS, but the filter is wrong in two steps, both instructive:
- First version caught 24 026 numbers, almost all 1/6, 4/3 and 8/3 printed as long floating-point numbers.
  A simple fraction in long format is not false precision, it's the opposite problem.
- The second version required denominators ≥ 50 and instead captured 1/86400 and 1/60000, thus seconds per
  days and minutes — unit conversions, not measurements. It also reads its own output file, which is
  a feedback loop I didn't clog.

## The task
Build the discriminator that separates three things, because that's the whole question:
1. **Two short measurements divided** — dangerous, must be labeled. Signature: both numerator and denominator are
   credible readings, three to five digits, none of them a round factor of unity.
2. **Simple fraction or unit conversion** — harmless. 1/6, 4/3, 1/86400, 1/60.
3. **Real long calculation or irrational number** — harmless, leave.

Exclude output files under `results/ASSEMBLY_PRECISION_AUDIT` and `ASSEMBLY_QUANTITY_INDEX` from the input.

## Strongest control
The swarm's own case: 2.220261437908497 should end up in class 1. And at least ten numbers you choose from the internet
that you first grade by hand — if the discriminator doesn't match your hand grade on them, it isn't
clear.

## Forgers
Before the run, print out how many of the grid's numbers you expect in class 1. Catches the discriminator more
than a fifth of all numbers it is too wide and captures unit conversions; catches it zero is it for
tight and misses the swarm's own fall.

## Why it's worth a spin
Every number we publish with more numbers than it carries invites the next reader to compare on numbers like
does not exist. It's the same error as the convention collisions, one level down: not wrong frame, but wrong
resolution. A precision marking per number makes it impossible to redo.
