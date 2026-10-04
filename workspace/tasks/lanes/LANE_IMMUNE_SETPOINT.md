# LANE_IMMUNE_SETPOINT

What sets the boundary between underactivation and overactivation in the immune response, and is it a point or an interval? Results directory `results/LANE_IMMUNE_SETPOINT/`.

## Why this lane

Anton explicitly seeded the immune system’s **underactivity and overactivity** as a track where recent scientific advances should be sought as data anchors. The track has received swarm jobs but no Sol lane, and it is the only one of his seeded directions lacking a deriving lane.

Internally there is already much to stand on, so this is not an empty field:

- `BT-RESEARCH-IMMUNITY` is one of our largest nodes with about 390 completed results.
- A swarm job today used CRP’s **half-life invariance** as a conservation constraint and rejected four of seven arms in an acute-phase model: all four were clearance-channel arms, and the three production-channel arms did not fail (but were not thereby validated). The anchor was radioiodinated human CRP with half-life 18,8 ± 3,9 h, constant across inflammatory states. It is a strong pattern: an invariance rejects an entire class of mechanisms without measuring the mechanism.
- Nine `BT-IMMX` jobs (IL6 buffer, IL2-Treg dose, PD1 occupancy, LPS dose tolerance, IgE threshold, IgG pneumonia, HLA-DR state, inflammaging, two-axis map) died at a provider limit 00:04–00:15 without producing anything. Their questions are thus formulated but unanswered.

## The question, and why it is a derivation and not cataloging

"Underactivation" and "overactivation" are used as if they were two sides of a scale with a point between them. If instead there is an **interval** where both errors are possible simultaneously, or where neither is defined, the entire clinical decision logic is wrongly formulated. This is a structural question that can be decided from a model, not from more measurements.

Test this concretely: an immune response that both fails to clear a pathogen AND damages the host is not a point on a scale between too little and too much. If that configuration is possible in a mechanistic model with published parameters, "setpoint" is an incorrect concept and should be replaced with a region.

## Do this

1. **Choose ONE mechanism where both errors can be expressed in the same quantities.** Candidates to test, not assume: pathogen clearance versus tissue damage in an acute response; Treg-mediated suppression versus effector expansion; endotoxin tolerance where a second stimulus gives a weaker response. Say why the chosen one can carry both errors simultaneously.
2. **Derive whether the two errors are separated by a point or a region.** Calculate, with published parameters and their spread, whether there are parameter combinations where both error criteria are satisfied simultaneously. Report the measure of that set, not only whether it is empty.
3. **Use an invariance as a constraint, as the CRP job did.** A quantity constant across states rejects an entire mechanism class without measuring it. Seek at least one such invariance with locator, and say which class it rejects. It is the lane’s strongest possible move and is already tested in our own material.
4. **Keep the anchor outside the model.** Every published measurement is held data and must never be input. State the type of spread per source — SD, SE, 95 % CI, IQR and range are not interchangeable.

## Strongest control and falsifier

- **Control:** the usual description with one threshold, i.e. a point between underactivation and overactivation, calibrated on the same data. The gain should be a statement about the set’s shape, not a better threshold value.
- **Falsifier:** if the set where both errors are satisfied is empty for every published parameter combination, the point picture is right and the lane has shown it. Report that just as clearly — a confirmed point picture with a computed margin is a result.
- **Prohibited:** calibrating against the outcome metric; treating a clinical cutoff as a measurement (it is often a committee convention); comparing a concentration against a flux without the denominator; counting an active fraction as a remainder.

## Delivery

`PORT.json` with the set’s shape (point or region) and its measure, the invariance used and which mechanism class it rejects, and the missing quantities with the measurement that would freeze each one. Narrow follow-ups to the swarm in FOLLOWUPS.json with external_referent complete — and **no template repeated per node**, because another lane today produced 24 identical follow-up proposals and only one of 31 carried external ground truth.

No internal data, no patient data. Everything PENDING_INDEPENDENT_REVIEW.
