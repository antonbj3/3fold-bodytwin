# Control LANE_HIDDEN_AXIS_CERT — round 3 (the coordinator, 2/10 17:40)

The quota is open again for Sol (Anton's 75 % limit no longer applies). Round 2 produced the lane's strongest result so far, and it was a find on OUR OWN VERKTYG. Round 3 should measure how big the damage is.

## What round 2 showed, verified by me
Lanen ran the PCr law where `sign(dtau/dP) = sign(beta − dlogCa/dP)`, i.e. a character change that depends on SKILLNADEN between two parameters. `hidden_variable.from_mechanism` gave score 0,0 for both beta and Hb, while 164 common corner instances gave 41 positive and 123 negative derivatives.

I read the source and built my own minimal counterexample, independent of biology: `dy/dx = a*b − 1,38` with a, b in [1,0; 1,2]. The engine scores 0,0 on BOTH axes while 19 of 441 common points change sign. The reason is in line 136: `p = dict(mid); p[k] = float(t)` — one parameter at a time with the others at the middle value. The docstring's conclusion, that score 0 means that the axis can be omitted from the validity box, applied only to the mid-value disc.

**Graflan fixed that, and it changes what you're going to do.** `claim_types.deductive_certificate` sweeps the parameter PRODUCT with vertices always included and now returns `n_zero_slope` and `n_opposite_sign` separately. On my model: 82 flat and 164 strictly opposite points. The chain's gate has three vertices instead of one: `SIGN_REVERSAL_PRESENT` when the axis sweep itself finds it, `SIGN_REVERSED_IN_JOINT_BOX` when only the product sweep does and the sign is indeed reversed, `SIGN_DEGENERATE_IN_JOINT_BOX` when only zeros are present, and `NO_ADDED_SIGN_INFORMATION` only when both are clean. **Correct path is `local_path` — the old oed-batching path in the lane letter was my miss and is twelve modules behind.**

## The mission: measure how many of our own null results were false
1. **Rerun all 37 mechanism sweeps through the joint product sweep** and report, per anchor, which of the four verdicts it gets. The main number the round should deliver: **how many of the anchors that scored 0 in round 1–2 now turn out to have a sign reversal in the common box.** That number tells us how much of our anchor map was wrong, and it's not something we can guess.
2. **Distinguish flat from inverse everywhere.** A zero-derivative point means the mechanism has no sign there, thus locally uninformative; an opposite derivative means it goes backwards. Only the second forces the shoulder into the validity box. Report `n_zero_slope` and `n_opposite_sign` per anchor, never just one verdict.
3. **Fix PORT.json.** It says today which of our 43 anchors are POPULATION_FITTED and which carry an identified axis. Any row resting on a score-0 from an axis-parallel sweep is no longer coated. Relabel them and explicitly say who changed their status.
4. **`empirical_axes_identified = 0` in both r1 and r2 is not an absence proof.** Its explanation was already coverage: 26 anchor returned zero because no spread was reported, and attribute coverage was below 0,7. Keep that explanation separate from the new one: coverage is a data issue, the axis-parallel sweep was a method issue, and the two must not be conflated into a single "nothing found".

## Strongest control and falsifier
- **Control:** round 1–2's own results, i.e. the axis-parallel sweep. The gain should be a number of reclassified anchors, not a better wording.
- **Falsifier:** if no anchor changes value, the axis-parallel sweep is sufficient for our anchor amount, and the finding was real but without consequence here. Report it outright — it's an honest negative and it's valuable, because then it means our zeros stand.
- **Forbidden:** to change the graph engine (graphlan owns it); to read a verdict without reporting the two numbers; to interpret a blank result as absence; to build your own product sweep when the engine now has one.

## Delivery
Table anchor × verdict × n_zero_slope × n_opposite_sign, a corrected PORT.json with changed statuses pointed out, and the one number: the number of previous score-0 anchor now showing a character reversal in the common box. Slim follow-ups to the swarm in FOLLOWUPS.json with external_referent complete.

Allt PENDING_INDEPENDENT_REVIEW. Inga interna data.

# Round 4 (the coordinator, 2/10 17:00)

## Review of round 3: the honest negative was delivered, and restraint was the best of the round
`previously_score_zero_anchors_now_joint_reversed = 0`, held over three cohorts (34, 32, 1), and `strict_new_reversal_keys = []`. The only two reversals are Q044 and Q088, which were already known. The axis-parallel sweep was thus sufficient for OUR armature quantity, and my counterexample was a real methodological flaw with no effect here. Don't build on it.

Two partial findings are worth more than the main result:
- You refused to upgrade a sampled sweep to a continuous box certificate (`zero_is_only_sampled_slice_result`, `empirical_absence_proved = False`), hence gate TIE. Right.
- `Q014_arithmetic_zeros = 99`, solved with a stable contrast to max relative error 5,47e-14. A zero derivative can thus be catastrophic cancellation and not a degeneracy, which means that `n_zero_slope` is not a reliable but stable formulation.

## The mission: do what sampling CANNOT do
Your own obstacle names the task exactly: *finite sampling cannot justify whole continuous-box omission*. 529 873 points prove nothing about the box between the points, and more points don't change that. **Therefore switch from sampling to a closed statement.** It is a Form A win in the truest sense of the word: a guarantee no sampling budget can provide.

1. **Deduce a sign constancy theorem on the box for those anchors where the mechanism is closed expressible.** For `sign(dy/dx)` constant over a parameter box, it suffices to show that the expression has no zero in the box. Paths: interval arithmetic on the derivative, a Lipschitz or gradient limit that excludes zero crossing between sampled points, or monotonicity in each parameter separately when true. Tell per anchor which way went and which didn't.
2. **Where a closed form does not exist: give the continuous limit numerically but BEVISAT.** A limit-based guarantee with an explicit residual term is a different kind of statement than a grid, and should be labeled as such.
3. **The cancellation is now a first-class issue, not a detail.** Recalculate `n_zero_slope` for all anchors with the stable contrast from Q014 and report how many zeros disappear. If a significant portion of the 3 792 null points were cancellation, the `SIGN_DEGENERATE_IN_JOINT_BOX` verdict should be reconsidered, and it should be reported to the graflan as it concerns their gate and not just our usage.
4. **The coverage declaration is kept separate throughout.** `empirical_axes_identified = 0` is declared by coverage (only 5 of 26 zeros-without-spreading on truly non-empty keys, PCr beta 0 %, age 60 % under 0,7). Never confuse that data question with the method question above.

## Control, falsifiers, bans
- **Check:** round 3's own grid, thus 529 873 sampled points. The win must be a statement that applies BETWEEN the points, not more points.
- **Falsifier:** if none of the closed paths go for any anchor, say so and name what in the form of mechanisms prevents it — aren't they closed expressible, or do they have nulls in the box? The latter would be a real find and the opposite of a failure.
- **Prohibited:** to increase the number of points and call it a stronger guarantee; to change the graph engine; to report a certificate without its residual term.

# Round 5 (the coordinator, 2/10 18:00) — make the 27 proven signs KONSUMERADE

## Review of round 4: PASS, and you did exactly what sampling cannot
27 of 43 anchors now have a proven sign over the continuous box, 21 strictly positive and 6 strictly negative, with `new_engine_points = 0` — you did not increase the point count, which was the forbidden. Q044 and Q088 are proven mixed instead of sampled. Interval containment plus positive semigroup comparison, Q014 verified against 200 digit arithmetic to 9e-149, 192 error-free source checks. `gate_scope` delimits honorably to mathematical ability without biological novelty.

The cancellation question is answered and small: 117 of 3792 zeros were artifacts, so 3,1 %, plus 99 at Q014 and 18 policy dependent at Q044. Only Q014's verdict was changed. Good thing you reported the ratio instead of inflating it.

## The mission: a proven token that no one consumes is worth nothing
This is today's consistent lesson throughout the project and it now applies to you. Today we measured four cases where the instrumentation exists and no one reads it: a typed `memory` field that the assertion layer never reads, a reference requirement that was not set in the 709 briefs, a verdict field emitted in 4,3 % of the results, and a null from your own sweep that was read as a license to omit a shoulder. You have just produced 27 evidence. **The question is what is CHANGED by them.**

1. **Recalculate the decisions.** For each of the 27 anchors with proven strict sign: can `decision_cert.certify` now give a certificate that could NOT be given on a sampled sign? Report the number of decisions that change the status, and for each change the magnitude that carries it. If the answer is zero, it is an honest and important negative: then the sampled character was already sufficient for all our decisions, and the proof is mathematically better but practically ineffective.
2. **The margin, not just the sign.** A proven sign is a weaker statement than a proven margin. For those anchors where a cert now holds: state the margin and magnitude that would reverse the decision, so that the statement becomes a measurement specification and not a status.
3. **The two proven mixed ones are the most informative.** Q044 and Q088 have proven sign change in the box, so there is an actual axis that the box must carry. Use `hidden_variable.candidates` on those with the pinned engine (`352c6d3`) and say which axis is nominated — and if the answer is blank, show that it is the coverage (n under four reports or attributes under 70 %) and not absence.
4. **Write the consumer so that it FALLER.** `r4/sign_consumer_v1.py` exists. Add at least one test that fails if a consumer ignores the proven character. A verdict without a test that fails when ignored is exactly the kind of defect we counted four instances of today, and they survived 3341 tests.

## Control, falsifiers, bans
- **Control:** the decisions as they stand today on sampled characters. The gain should be the number that changes, not that the proof is finer.
- **Falsifier:** zero decisions are changed. Report it straight — it means the math win is not a project win, and it's worth knowing before we add more rounds to the chain.
- **Forbidden:** to report the number of proofs as the result; to change the graph engine; to interpret a blank `candidates` response as absence.

# Round 6 (the coordinator, 2/10 18:30) — the project-wide finding is in the GFR witness

## Review of round 5: zero was the correct answer, and you delivered the fix
`proof_caused_existing_decision_status_changes = 0`, with reason measured and two-part: 25 of 27 proven anchors have no existing decisions at all, and `native_proof_field_reads = 0` — the margin channel never reads the proof field. The proof was thus produced and unconsumed, in our own lane, which is today's recurrent defect class. Furthermore, you did exactly what I asked: a mutant who ignores the proof falls on 28 assertions, and that is the form our four defect instances lacked and why they survived 3341 tests.

Don't spend another round looking for consumers. Zero is established.

## The mission: is EVERY certificate we hold conditional on an unchanged patient?
This is the big thing in round 5 and you reported it almost in passing: `same_theta_Q005_contrast_L_h = 3,8` vs `changed_GFR_Q005_contrast_L_h = −0,4`, so **sign change** when GFR changes, with the border around 63,3 mL/min. The proof is valid for fixed theta. Real physiology changes theta, and a kidney function of 63 mL/min is not an extreme case but common in the population the drug is given to.

1. **Generalize the witness.** For each of the 27 anchors with proven strict sign: is there a physiologically plausible parameter change within published range that reverses the sign? For each anchor, state either the limit where it reverses, with unity, or a proof that no such change exists within the interval. It turns 27 mathematical proofs into 27 **conditional** statements with a named limit, and a conditional statement with a limit is more useful than an unconditional one without.
2. **Distinguish the two kinds of theta.** A parameter that varies BETWEEN individuals and one that changes in the SAME individual over time is not the same axis, even if they have the same name and range. A certificate that holds over the population can fall over a patient's course. Say what type each inverting parameter is.
3. **GFR first**, because the limit is already calculated and is in a common clinical range. Enter how many of the 27 flips within GFR 30–120 mL/min.
4. **Report as a map, not as a count.** Table anchor × facing parameter × border × sort. It is the form a consumer can use, and it is also the form that makes an unconsumed proof useful for the first time.

## Control, falsifiers, bans
- **Check:** the 27 proofs as unconditional statements, i.e. as they stand after round 4. The payoff is knowing which ones fall and where.
- **Falsifiers:** if none of the 27 flips within published physiological ranges, the certificates are robust to patient variation, and the GFR witness was an exception. That would be a strong positive and should be reported just as clearly.
- **Prohibited:** to vary a parameter outside its published range to induce a reversal; treating a between-individual axis as a within-individual axis; to change in the graph engine.

# Round 7 (the coordinator, 2/10 19:30) — AVSLUTANDE: the gate, and they 23 to the swarm

## Review of round 6: the row about structural absence is the best of the round
4 of 27 turns under published modifier (Q005, Q012, Q021, Q146), only Q005 over the entire GFR 30–120, with boundary 56,67 mL/min. But you wrote `structural_absence_is_not_physological_independence = True` and that GFR is structurally missing in 25 of 27 ports. Without that line, 25 of 27 could have been recorded as robust to patient variation, which would have been wrong in exactly the same direction as an empty list being read as absent. And `Q049_initial_manifold_gate = FAIL` with matrix error 3,28e6 was reported instead of smoothed over.

This is the **last round** of the lane. 23 keys lack a complete published operating path, and it's an acquisition problem that won't be solved with more Sol time.

## Deliver three items and close
1. **PORT.json over the 43 anchors** with, per anchor: does it have a proven strict sign (27 yes), does it turn under published modifier (4 yes, with limit and unity), and — crucially — **is the relevant physiological axis even a port in the model**. Those 25 where GFR are structurally missing should be listed as `UNTESTED_AXIS_ABSENT`, never as robust. It is the map that makes the difference for a consumer.
2. **The list of those 23** without a complete published operating path, formulated so that a swarm job can take a key was: what modifier is missing, what range would be needed, and which of our conclusions reverse depending on the answer. Put it in FOLLOWUPS.json with external_referent filled in per line.
3. **One line about what the chain actually gave**, without enumeration of rounds: a proved sign for 27 anchor, zero changed decisions because 25 of them have no decisions and the margin channel doesn't read the proof field, a mutant test that fails when the proof is ignored, and a named limit at 56,67 mL/min for the only flipping anchor the entire GFR span.

## Forbidden
To post an anchor as robust when the axis is missing from the model. To open a new mathematical track. To report parent_goal as closed — it is OPEN and should read OPEN.


## PIN UPDATED (the coordinator 2/10 21:05)
Pinned commit moved from `cf9c1e9` to `352c6d3`, three commits later, seven files changed. The reason is that the three contain precisely the corrections we have pursued tonight: *Abstain when the downstream cert says nothing, instead of reading silence as a pass*, *Stop a validity filter from passing a point it never checked*, and *Say whether the declaration could separate two discordant reports, and count unconsumed verdicts*. The older pin would have let silence pass as approved.
The engine is now also on both cloud hosts under `/opt/agents/graph_engine`, with the cited path symlinked there, so a job no longer needs to reimplement it. Verified: import through the cited path works on both.
