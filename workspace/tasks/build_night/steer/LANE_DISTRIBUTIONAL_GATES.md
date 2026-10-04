# LANE_DISTRIBUTIONAL_GATES — build the acceptance gate our cells lack

Resultatmapp `results/LANE_DISTRIBUTIONAL_GATES/`.

## The measurement that makes the lane necessary
A swarm job swapped the hemostasis cell's point-in-band gate for a **distribution test omitted** and held
everything else the same — same model, parameters, integrator, 15 000 step mesh and same published
mean and spread, only the decision statistics are different. **The control PASSERAR, the treatment FALLER.**
I opened the raw data and checked:

| magnitude | value |
|---|---|
| mass that never ends | **0,07305** towards the threshold 0,05 |
| same mass at time step 0,020 / 0,010 / 0,004 | **0,069 / 0,069 / 0,069** |
| median | 2,6793 min against the referent's 2,79 ± 0,78 |
| q95 and tail | 4,8392 min respectively up to 12,7420 min |
| overlap | OVL 0,7726, BC 0,9185 |

That the non-inference is **identical at three time steps** makes it structural and not numerical. And
the job declares itself that the FAIL verdict is prior-dependent: with K_OFF_P fixed, the test passes at
±20 % and ±40 %. **The conclusion is therefore the gate's FORM, not the model's case:** a scalar acceptance test can
don't see that 7 % of the cases never end.

And the distribution over our cells, measured before the statement: **21 of 43 cells carry a scalar point or
bandgate, only 7 carries any distribution statistics.** Eight cells have scalar gate without any
distribution statistics total: `IMMUNITY`, `MITOSTRESS`, `Q009`, `Q031`, `SOLBENCH`, `SURG_HEALING`,
`SURG_HEMOSTASIS`, `SURG_INCISION`.

## Do like this
1. **Build EN reusable gate module** in the workspace that takes a cell's predictive outcome and returns:
   quantiles (q05, median, q95), tail mass above a declared limit, **the mass that never reaches the target**,
   and overlap against a reference distribution when one exists. One module, not eight copies.
2. **Run it on the eight cells** and report per cell: passed the scalar gate, it passes
   the distribution gate, and **reverses the value**. It is the main number of the lane: how many out of eight flips.
3. **Distinguish prior dependence from structural case.** For each trap: rerun with the most uncertain one
   the parameter fixed at its supplied value. If the cell still falls, it is structural; passes it is
   the prior dependent and should be reported as such. That was exactly what the swarm job did right and that
   is the difference between a find and an alarm.
4. **And check the time step stability** for each non-closing mass, with at least three steps as in
   the measurement above. A tail mass that moves with the time step is numerical and not a tissue property.
5. **For the cells without a reference distribution: say so** and report the quantiles anyway. A cell that does not
   can be compared against something is an acquisition item, with greatness and unity, in the form i
   `notes/ACQUISITION_TARGETS.json`.

## Strongest control and falsifier
- **Control:** cell's NUVARANDE scalar gate, run on the same outcome. It is an equally informed
  control by design — same model, same data, only the decision statistics differ — and that's it
  therefore the result can be claimed at all.
- **Falsifier:** if none of the eight verdicts reverse, the scalar gate is sufficient in practice
  and the hemostasis case is an exception. That would be a reassuring result and should be reported outright, not
  tonas ner.
- **Forbidden:** to change any cell's physics or parameters — only the gate is added; to report
  a precipitation without the prior tuning in point 3; to call a tail mass structural without
  tidsstegskontrollen.

## Delivery
`PORT.json`: gate module path, one row per cell with both vertices and if it reverses, per fold
its prior and structure classification, respectively, and the time step table for each non-closure mass.

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
