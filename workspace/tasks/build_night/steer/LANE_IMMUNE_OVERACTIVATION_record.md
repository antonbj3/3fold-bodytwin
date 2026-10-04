# Steering LANE_IMMUNE_OVERACTIVATION — after r1

## You refuted the axis, and I have lowered the edges
`rho_status = UNMEASURED_MODEL_CLOSURE` . The recycling rate is a model closure, not a measured
quantity, and the lysis thresholds are synthetic. The falsifier I wrote in the brief thus triggered, and
you were right to leave the shared graph untouched and leave the decision with me. `T-E14` and
`T-E15` are now **UNKNOWN instead of TIGHT**, and `T-E30` is qualified as a numeric result.

And the inadequacy is grosser than I thought: two states with **identical** summations (both
identity errors exactly 0,0) differ by **17,60797843701942 nM** in MAC. I set it against your
primary range `[0; 0,5441770553588867] nM` — the gap is **32,36 times the full range** the
decision is working in. The stressed range is 0,4324131012 nM, 79,5% of the primary.

## Obstacle
The entire chain hangs on a axis that no one has measured. Anything
you calculate on top of it is correct and useless for a decision.

## Changed operation
Change axis. Find a quantity in the same chain that **is** measured and for which ρ is a function of
— surface-bound C3b density, regulator concentration, convertase half-life (we have 129,84 s surface-bound
and 111,09 s in fluid phase from our own complement decision). Express the window at that size instead.
If none exists, it's the result: print what measurement would make the axis measurable.

## Strongest control
The one-sided reading as it stands now, with ρ as a free parameter. Every
new axis should beat it by being measurable, not by being slimmer.

## Falsifier
If the new axis also does not have a published measurement with a locator, then the
entire complement regulation in the network is a model statement, and then all six immune
edges should be UNKNOWN with that reason — not just the two I lowered today.
