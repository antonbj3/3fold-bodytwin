# the proof lane: when is improving a component counterproductive?

## The measured pattern, seven chains in two sessions

Seven independent chains carry a favourable number that requires CANCELLATION to be retained. Each time,
an isolated improvement of a component makes the system WORSE.

**My four, with numbers:**

1. **Toric system margin.** Two stages, per-stage standard deviation 0,45562 and 0,59698 D, measured
   correlation ρ = −0,7961 (95 % interval −0,8565 to −0,7254, n = 69). Combined sd 0,36178 against
   independent 0,75098, so cancellation removes 51,8 % of the error. Repairing stage 1 completely gives
   combined sd 0,59698, a factor **1,6501** worse; stage 2 completely gives 0,45562, factor **1,2594**.
   The interior optimum for the two-stage case is σ₂* = −ρ·σ₁ = 0,36273, and stage 1 is already at 0,9586 of its
   optimum, thus marginally BETTER than benefits the system.
   The gate that certified this is called "decorrelation" and is ONE-SIDED: the rule is that the correlation's
   upper confidence bound lies below 0,15, which lets strong anticorrelation through under a name
   that means absence of correlation.

2. **Laser.** Sharp endpoint bounds exist only "after spatial cancellation is retained".
   The aggregate strictly underestimates in **14 412 of 15 413** positive instances, thus 93,51 %.
   The witness is exact: with the weight pair q = (1/2, −1/2), which sums to zero, and the summary
   S_N = Y = 0 exactly, the source release is 0 against 4/3 Pa·s with threshold 2/3 exactly halfway between.

3. **The swarm's error budget.** Reported spreads underestimated the truth 1,90× to 2,74× through
   superposition of opposite signs. Verified on one case: 0,01172341/0,00586116 = 2,00019.

4. **The field lane's three**, their own words: the needle force's effective J is "first minus repeated" and
   the friction correction cancels exactly only if contact is the same on reinsertion; the incision's
   guaranteed energy difference is narrow only when shared errors between endpoints are kept correlated;
   and the guarantee for force per tooth holds only when a shared geometry error is carried through the entire
   coupling.

## The question

There are two readings and I cannot decide between them:

**A. Property of the systems.** Biological and mechanical chains work at compensated operating points,
where terms with opposite signs balance. Then cancellation is physical, the margin is real, and
the rule is never to optimise a part without the pair.

**B. Property of our error composition.** We estimate correlations from the same finite sample from which we estimate
the spreads, and a negative estimate is favourable by construction. Then the margin is partly an
artifact and would shrink with more data.

Which? And if it is A, what is the exact theorem?

## What I want built

1. **Generalise the interior optimum to K stages.** For the two-stage case, σⱼ* = −ρ·σ_andra is exact. For K stages
   with covariance matrix Σ, what is the condition for a reduction of σⱼ to INCREASE 1ᵀΣ1? State it exactly,
   not as a rule of thumb.

2. **State when a component improvement is counterproductive**, as a condition on measurable quantities.
   Which numbers must be known to know it in advance?

3. **Distinguish A from B with a test.** Construct a test that on a FINITE sample distinguishes physical
   anticorrelation from estimation artifact. It may use the fact that we have n = 69 for the toric case
   with a CI on ρ, and n = 15 413 instances for the laser case.

4. **Judge the one-sided gate.** Is "upper CI bound below tol" a defensible certification rule for
   a system where negative correlation is favourable, or must it be two-sided? If two-sided, what is
   the right tolerance and why?

## Requirements

- Every claim must have a condition, not a feeling. A theorem must carry its assumption.
- Control: the reading that treats the stages as independent. It is what a reader uses today.
- Falsifier written out before every derivation.
- Say what you CANNOT determine. A boundary on the claim is a result.
- No claim about real tissue. This is a question about error composition and certification.
- PENDING_INDEPENDENT_REVIEW.
