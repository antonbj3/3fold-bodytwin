# LANE_DISTRIBUTIONAL_GATES — build the acceptance gate our cells lack

Resultatmapp `results/LANE_DISTRIBUTIONAL_GATES/`.

## The measurement that makes the lane necessary
A swarm job replaced the hemostasis cell’s point-in-band gate with an **omitted distribution test** and kept
everything else equal — same model, parameters, integrator, 15 000-step grid and the same published
mean and spread, only the decision statistic differed. **The control PASSES, the treatment FAILS.**
I opened the raw data and checked:

| quantity | value |
|---|---|
| mass that never closes | **0,07305** against threshold 0,05 |
| same mass at timestep 0,020 / 0,010 / 0,004 | **0,069 / 0,069 / 0,069** |
| median | 2,6793 min against the referent’s 2,79 ± 0,78 |
| q95 and tail | 4,8392 min and up to 12,7420 min respectively |
| overlap | OVL 0,7726, BC 0,9185 |

Non-closure being **identical at three timesteps** makes it structural and not numerical. And
the job itself declares that the FAIL verdict is prior-dependent: with K_OFF_P fixed, the test passes at
±20 % and ±40 %. **The conclusion is therefore the gate’s FORM, not the model’s failure:** a scalar acceptance test
cannot see that 7 % of cases never finish.

And the distribution across our cells, measured before the statement: **21 of 43 cells carry a scalar point or
band gate, only 7 carry any distribution statistic.** Eight cells have a scalar gate without any
distribution statistic at all: `IMMUNITY`, `MITOSTRESS`, `Q009`, `Q031`, `SOLBENCH`, `SURG_HEALING`,
`SURG_HEMOSTASIS`, `SURG_INCISION`.

## Do this
1. **Build ONE reusable gate module** in the workspace that takes a cell’s predictive outcomes and returns:
   quantiles (q05, median, q95), tail mass above a declared bound, **the mass that never reaches the target**,
   and overlap against a reference distribution when one exists. One module, not eight copies.
2. **Run it on the eight cells** and report per cell: did it pass the scalar gate, does it pass the
   distribution gate, and **does the verdict reverse**. That is the lane’s main number: how many of eight reverse.
3. **Separate prior dependence from structural failure.** For every rejection: rerun with the most uncertain
   parameter fixed at its delivered value. If the cell still fails, it is structural; if it passes,
   it is prior-dependent and should be reported as such. That is exactly what the swarm job did right and
   is the difference between a finding and an alarm.
4. **And check timestep stability** for every non-closure mass, with at least three steps as in
   the measurement above. A tail mass that moves with timestep is numerical and not a tissue property.
5. **For cells without a reference distribution: say so** and report the quantiles anyway. A cell that
   cannot be compared against anything is an acquisition entry, with quantity and unit, in the format of
   `notes/ACQUISITION_TARGETS.json`.

## Strongest control and falsifier
- **Control:** the cell’s CURRENT scalar gate, run on the same outcomes. It is an equally informed
  control by construction — same model, same data, only the decision statistic differs — and that is
  why the result can be claimed at all.
- **Falsifier:** if none of the eight verdicts reverse, the scalar gate suffices in practice
  and the hemostasis case is an exception. That would be a reassuring result and should be reported plainly, not
  downplayed.
- **Prohibited:** changing any cell’s physics or parameters — only the gate is added; reporting
  a rejection without the prior check in point 3; calling a tail mass structural without
  the timestep check.

## Delivery
`PORT.json`: the gate module’s path, one row per cell with both verdicts and whether it reverses, for each rejection
its prior versus structural classification, and the timestep table for each non-closure mass.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.

## ADDENDUM — extend to all 43 cells, and take the asymmetry seriously
Eight cells was a sample, not an answer. **3 of 8 verdicts flipped**, and 21 of 43 cells accept on a
scalar point-or-band test, so there are 13 cells we know nothing about. Extend to all 43 and report the
three numbers per cell: passes the scalar gate, passes the distributional gate, verdict flips.

And one thing the sample did not separate, which the field lane's inspection makes concrete. Their
deflection decision reports **74 of 100 boxes fully decided against 13 falsely safe answers under corner
testing**. Falsely safe is the dangerous direction and a symmetric error count hides it. So for every
flip, say WHICH WAY it flips: did the scalar gate pass something the distribution rejects (falsely safe),
or reject something the distribution accepts (falsely cautious)? Those are not the same finding and only
the first one is a hazard.

Report falsely-safe and falsely-cautious counts separately for all 43. A gate that is wrong in the
cautious direction costs compute; one that is wrong in the safe direction costs a claim.
