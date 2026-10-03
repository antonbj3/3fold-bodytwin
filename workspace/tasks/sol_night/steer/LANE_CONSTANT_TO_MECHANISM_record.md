# Direction LANE_CONSTANT_TO_MECHANISM — round 2 (the coordinator, 2/10 16:15)

The quota is open for Sol again: Anton's 75 % limit for lane_runner no longer applies, and last night's STOP file was set to "save Sol until tomorrow", which is now.

- **Obstacle:** the measured enzyme activity 0,9625 is below the published 2,9–9,1 nmol/min/mg, and the lane has treated it as a capacity gap that must be carried by a tissue port.
- **Changed operation:** decide FIRST whether it is a gap or three normalizations. It is not a theoretical objection, it is a measurement I made today on an adjacent case and it changed the conclusion by a factor 160.

## Why normalization comes before the port
Today I retrieved Zhang & Wright 2022 (10.3390/ijms23031472) to decide a claimed gap of 243× in OCT2's turnover number. The paper turned out to state the number in **three** places that disagree: the abstract 12,6 ± 2,4 s⁻¹ at ~510 fmol/cm², Results 2.4 gives 70 ± 12, and Table 2 gives 62,2 ± 20,5 at 154 ± 70,7 fmol/cm². The abstract and table thus use **different surface densities** for the same normalization. After temperature correction from ~23 °C the claimed gap fell from 243× to 1,5–0,75×. The gap was a literature artifact, not a capacity deficit.

The same questions apply to 0,9625 against 2,9–9,1 nmol/min/mg, and they must be answered before a port is built to carry the difference:
1. **Per what?** nmol/min/mg microsomal protein, mg total protein, mg membrane protein or nmol enzyme are four different numbers. State for EVERY source, verbatim, what the denominator is. If they differ the numbers are not comparable and the gap is partly or entirely apparent.
2. **Active fraction.** An activity per mg protein measures the preparation, not the enzyme. If the active fraction differs between preparations the ratio of activities is not the ratio of catalytic capacities.
3. **Temperature and pH.** State each source's values. Report the correction for Q10 = 2, 2,5 and 3 **separately** and say that the choice of Q10 is an assumption, not a measurement.
4. **Species, preparation, substrate concentration and saturation.** An activity below K_m scales with substrate, not with capacity.
5. **Errata.** Search for corrections and corrigenda for every source, and report which phrases you searched with even when you find nothing.

## Strongest control and falsifier
- **Control:** conventional IVIVE with the same input, plus the lane's own earlier port. The gain shall be a statement about whether the gap EXISTS, not a better fit across it.
- **Falsifier:** if all sources share denominator, active fraction, temperature and degree of saturation, the gap is real and shall be carried by the mechanism. Then say so plainly and build the port — but with individual ratios and an uncertainty budget, never a fitted organ-Vmax.
- **Forbidden:** choosing the normalization that makes the gap smallest; treating a single Q10 as measured; fitting organ-Vmax against the endpoint. The anchor is held-out data and must never be model input.

## Delivery
A table source × denominator × active fraction × temperature × saturation, and a verdict: the gap is apparent (with the factor by which it shrinks), partly real, or real. Narrow follow-ups to the swarm in FOLLOWUPS.json with external_referent complete: kind, locator, compared_quantity, our_value, command or excerpt, refutes_us.

Allt PENDING_INDEPENDENT_REVIEW. Inga interna data.

# Round 3 (the coordinator, 2/10 16:45)

## Review of round 2: the question was worth asking, but was not the answer
`normalization_verdict = PARTIALLY_REAL_UNRESOLVED_CAUSE`. The assay method carried 1,17951 and raised the prediction 0,9625 → 1,13528; I recalculated it and the arithmetic holds, and `remaining_held_over_prediction` × the prediction gives exactly the held-out band 2,9–9,1. The temperature factor became 1, thus the opposite of the OCT2 case where temperature was the entire explanation. **Normalization explained one fifth of an order of magnitude, not the gap.** Remaining: 2,55–8,02× and `intrinsic_capacity_deficit = UNKNOWN`.

Therefore stop searching in normalization. Your own obstacle names the next term correctly: different cohorts, counting metrology and **competence**, meaning what proportion of the inventoried enzyme is actually catalytically active.

## Changed operation: active fraction as a measured quantity, not a residual
1. **Active fraction must not be the term that absorbs the gap.** If it is calculated as the ratio between measured activity and inventoried amount it becomes by construction exactly as small as the gap requires, and then nothing has been explained — it is the same circularity that another of our lanes was just rejected for, where a drag coefficient was back-solved from the speed it was supposed to explain. Seek an **independent** measurement of active fraction or specific activity per mol enzyme: titration with a covalent or tight-binding inhibitor, active-site titration, or specific activity of purified recombinant enzyme. Freeze it with a locator before you calculate.
2. **If no such measurement exists: state it as a measurement specification, not an obstacle.** Name the titration, which preparation it requires and which quantity it would freeze. It is a fully adequate outcome for the round.
3. **Calculate an upper bound that needs no active fraction at all.** If every inventoried enzyme molecule is assumed 100 % active you get an upper bound on capacity. If the held-out activity 2,9–9,1 is ABOVE that bound the inventory or amount measurement is wrong, not competence, and this is decided without knowing the active fraction. **Do that calculation first** — it is cheap and can move the question alone.
4. **The cohort question separately:** if the held-out activities come from different individuals than the amount measurements, their ratio is not a ratio within the same individual. State per source whether amount and activity are PAIRED in the same preparation, and calculate the gap only on paired rows if such exist.

## Control, falsifier, prohibitions
- **Control:** conventional IVIVE with the same input, plus your own linear programming comparator. The gain shall be a statement about WHERE the gap is, not a better fit across it.
- **Falsifier:** if the 100 %-active bound is above the held-out activities, competence is a possible explanation and shall be measured; if it is below, the inventory is rejected. Report which.
- **Forbidden:** calculating active fraction as a residual; multiplying several corrections so none can be attributed; treating a single Q10 as measured; fitting organ-Vmax against the endpoint.

# Round 4 (the coordinator, 2/10 17:15)

## Review of round 3: the cheapest calculation decided most, and it excluded an entire class
The 100 %-active ceiling became 1,0164 pooled and 1,1984 radiometric against a held-out band 2,9–9,1, thus a violation of 2,85–8,95× and 2,42–7,59× respectively (recalculated by me). Since a lower active fraction can only LOWER the prediction, competence is thus excluded without any active fraction being measured. It is the right kind of result: an exclusion that did not require the missing measurement.

And you yourself found what is sharper than the main result: **pooled versus paired**. On paired rows the violation is 1,1533× (k_cat 91) and 1,0946× (k_cat 96) against held-out 3,3016–3,4709. Pooled, the gap is multiples; paired, it is around ten percent. The large part of the gap thus looks like a cohort artifact from amount and activity coming from different individuals.

## The assignment: decide whether ANY gap survives on paired data
This is the round's only question, and it is decidable.
1. **Recalculate everything only on paired rows**, with uncertainty propagated per row and the type of dispersion named per source. Report the paired gap as an interval, not a point. If one at 1,0 lies within the interval there is no gap to explain on paired data, and then the entire original discrepancy is a cohort and metrology question — say so plainly.
2. **Quantify the cohort contribution explicitly.** The difference between the pooled and paired gap is a measure of how much of the discrepancy is created by the lack of pairing. State it as a factor with uncertainty. It is a statement about our own measurement setup and it is reusable for every other family where we pair amount from one cohort with activity from another — and we do that in more places.
3. **Attack the amount measurement, not competence.** Your own conclusion points there. What sets the inventoried count: which standard, which peptide, which normalization? Is there an INDEPENDENT amount measurement with another standard or another method? A ten-percent deviation on paired data is precisely the order of magnitude of a standard or normalization difference, so it may well be the entire remainder.
4. **The 1,1× violation must not be rounded away.** If it survives propagated uncertainty it is real and small, and a small real violation of a 100 %-active ceiling is more interesting than a large one, because the ceiling is a hard ceiling: it says the inventory is too low by a measurable factor.

## Control, falsifier, prohibitions
- **Control:** round 3's pooled figures, thus the comparison that gave multiples. The gain shall be to say how much of them was lack of pairing.
- **Falsifier:** if the paired interval contains 1,0 there is no capacity gap on paired data, and the lane's original obstacle was wrongly posed. Report it as a fully adequate outcome — it is a resolution gain, not a failure.
- **Forbidden:** calculating active fraction as a residual; mixing paired and unpaired rows in the same ratio; treating a single Q10 as measured; rounding away 1,09× without propagated uncertainty.

# Round 5 (the coordinator, 2/10 18:00) — CLOSE the track with a port, and explain the regime factor

## Review of round 4: the question dissolved, which is the best kind of outcome
The cohort contrast 2,1690–7,7631 reproduces the held-out band exactly (2,169 × 1,3370 = 2,9; 7,763 × 1,1722 = 9,1). Pooling alone thus carried most of a gap we have treated as physics for five rounds. The paired residual gap 1,0437–1,4742 at 1,96 SE requires only 2,1 % extra symmetric error to reach 1,0, and at 3 SE or Student df 3 the interval already contains 1,0. And the independent amount method's cross-basis gives [0,9400; 1,0599], which contains 1,0.

**The track is thus practically exhausted, and in the right way: the question was largely wrongly posed.** This is the closing round. Spend no budget chasing the residual further.

## Two things to deliver, nothing more
1. **PORT.json with the verdict, quantified.** The gap 2,9–9,1 resolves into: cohort artifact 2,17–7,76×, assay normalization 1,18× (round 2), and a residual 1,04–1,47 whose existence depends on the error model. State for each term what it rests on and which type of dispersion was used. Write explicitly that the residual is not robust to the choice of degrees of freedom, and that it is decided by a number the sources do not report.
2. **Explain `kinetics_regime_factor = 2,5069`.** It stands unexplained and it worries me, because exactly that error class was rejected in another of our chains today: a requirement calculated in a saturated regime was compared against a measured value in a linear regime, and they differ by Km/C. Decide which regime each part of your comparison is in: is substrate concentration far below, near or above K_m? If any term is calculated as V_max while being compared against something that is in practice V_max/K_m, the factor 2,5 is a regime artifact and not biology. Show the arithmetic both ways.

## The measurement specification is a REPORTING SPECIFICATION this time
What would decide the residual is not a new experiment but a number that was produced and not printed: the raw degrees of freedom and covariance between amount and activity determination in the same preparation. Formulate it as a short, concrete request that can be put to an author or retrieved from supplementary material, and say which of our conclusions reverses depending on the answer.

## Control, falsifier, prohibitions
- **Control:** round 1's original claim of a 3–9× capacity gap. The gain shall be its resolution, not a smaller residual.
- **Falsifier:** if regime analysis shows that 2,5069 IS biology and not an artifact, a real unexplained factor remains and the port shall say so plainly.
- **Forbidden:** chasing the residual with more error models; choosing the error model that gives the desired answer; closing the lane without deciding the regime question.
