# Direction — round 1 (coordinator, 2/10 21:30)

New lane in a Sol wave. Anton 2/10: "it is up to you to fix all this and make it work, you know the goal is comparative improvement". The wave builds connections; it revises no numbers.

Read the lane brief in `tasks/lanes/` for the task. These are the additions that apply to all lanes in the wave.

## The graph engine is now in the cloud
It is at `/opt/agents/graph_engine` on both cloud hosts, with the quoted path symlinked there, pinned commit `352c6d3`. Import via the path is verified. **Do not reimplement it** — four jobs today did so silently because it was missing, and one gave a result I had to qualify because its "engine" was its own rewrite.

## Five error classes that cost us today, each with its price
1. **Wrong regime.** A requirement was calculated as flux = k_cat × N, i.e. saturated, and compared against a value in the linear regime where clearance = V_max/K_m. The factor was over a hundred and it rejected today’s largest claim. State which regime and denominator every number belongs to before comparing.
2. **Pooled is not paired.** A gap of 2,9–9,1× turned out to be a 2,17–7,76× cohort artifact because amount and activity came from different individuals.
3. **Structural absence is not independence.** 25 of 27 anchors did not change sign under a physiological modifier — because the modifier was not even a port in the model. Label such things untested, never robust.
4. **Prediction is not measurement.** A field with *conditional* in its name was read as measured, by me. Label every number MEASURED or DERIVED.
5. **A verdict nobody reads is worthless.** Four instances today. Every verdict you produce should have a test that FAILS if a consumer ignores it.

## Missing file
`missing_prerequisite` in the FIRST paragraph of RESULTS.md, not a footnote. Never silently reconstruct a missing input — four jobs did that today and one produced an incorrect headline that I amplified further.

Everything PENDING_INDEPENDENT_REVIEW. No internal data leaves the machine.

# Round 2 (coordinator, 2/10 22:15) — bisect the error, do not touch the chain
The chain RUNS. 18,0009 µM against held median 3,7 = 4,87× above, with 7 added assumptions. It is the project’s first working connection and must not be rebuilt.
- **Obstacle:** seven assumptions were added and we do not know which carries 4,87×.
- **Changed operation:** bisect. Set one assumption at a time to its most conservative value and measure the ratio. Report error contribution per assumption, sorted.
- **Control:** the two cells separately against their own anchors. The gain is knowing WHICH term carries the error.
- **Falsifier:** if no single assumption carries more than a factor of 1,3, the error is distributed and the chain is structurally wrong, not miscalibrated. Say so.
Prohibited: adjusting several assumptions simultaneously; calibrating against 3,7.

# Round 4 (coordinator, 2/10 22:50) — the missing clearance is not a parameter, so seek the MECHANISM
You have dissected the error and the answer is structural, not numerical. Active transport alone diverges
(`active_only_finite_steady = False`, minimum growth 66,14 µM/h, 22 079 µM at 240 h). With a passive
fraction, equilibrium becomes finite but 5,81–7,81× above held median 3,7, and the span shrinks only 1,34×
when the passive fraction grows fifty times. **Thus the missing clearance is not the passive
fraction.**
- **Changed operation:** name the candidate mechanisms for the missing sink and compute what
  order of magnitude each must have to close 5,81×. Candidates to test, not assume:
  microbial degradation before absorption, hepatic first pass, renal secretion beyond
  filtration, and tissue binding. For each: what quantity would it need to have, and is it
  measured? Search ALL FOUR pools (see COMMON.md) before saying something is missing.
- **The 21,48× segment gap** (29,04 µmol/h held against 1,352 our synthesis) is correctly labeled
  `CONDITIONAL_SEGMENT_SUBSET_ONLY`. Keep that label. But calculate what the whole-body value would
  be if the segment’s flux per length applied to the entire intestine, and say whether it is physiologically possible —
  if not, the segment measurement is not scalable and that is a separate finding.
- **Strongest control:** held median 3,7 µM with its own spread. Report 5,81× against it.
- **Falsifier:** if no named mechanism can carry 5,81× within published orders of magnitude,
  the chain’s STRUCTURE is wrong and not its parameters — and that is a stronger result.

# Round 5 (coordinator, 2/10 23:20) — now it must CLOSE, and be falsified hard
You found that the 5,81× error is two errors in opposite directions: our source 4,13× too low (5,378 against published
22,208 µmol/h) and our renal clearance 26,16× too low (4,167 against 109 mL/min). With both published values,
the chain lands at 3,3958 µM against held median 3,7, IQR 2,4–6,2. I have independently recalculated it from
your numbers and reproduce your `renal109_forward_TMAO_uM` to the last digit, so the scaling is valid.
**This is not a prediction, and the next round should treat it as suspicious.**
- **Strongest falsifier, and it must run first:** two errors in opposite directions that almost cancel
  are exactly the form of a post hoc fit. Test it: are the published 109 mL/min and
  22,208 µmol/h taken from the SAME population and the same state, or from two different studies? If
  they come from different cohorts, the composition is not licensed and the hit may be chance.
  Report source, population and state for both.
- **Second falsifier:** how large an interval of source and clearance lands within IQR 2,4–6,2? If
  almost any combination within published intervals lands there, the IQR is too broad to
  test anything, and the hit says nothing. Calculate the fraction of the published parameter rectangle that
  falls within the IQR. That number decides whether the result carries weight.
- **Third:** why was our own source term 4,13× too low? This is a model question, not a
  parameter question. Which pathway in the source term is missing?
- **Retain and strengthen your own negatives:** the microbial assay at 7 672× physiological concentration does
  not suffice as support, and uniform extrapolation to 506 cm gives 489,8 µmol/h = 16,9× the input, so
  the segment measurement is not scalable. Both should remain in the port.
- **Control:** held median with its IQR, never against zero.
**Prohibited:** reporting 0,918× as a prediction or model validation; choosing
the combination of published values that hits best.

# Round 6 (coordinator, 3/10 00:15) — the base-rate test rejected my headline number, build on what HELD
You ran the falsifier I set and it hit me: **69,9 % of the published rectangle hits
the IQR** (67,8 % loguniform), so the hit at 3,3958 µM against median 3,7 carries almost no information. I
have corrected my own log row. But **with the source’s own marginal, the fraction is 17,2 % and 15,0 %**, and there
it carries weight.
- **Changed operation:** make the stricter distribution the main statement. What is the source’s own marginal
  and why is it the right distribution rather than a rectangle? If the source’s marginal is the distribution
  a real population has, 17,2 % is the number to report and the rectangle is a straw man.
  Report both, with the stricter first.
- **Next real track:** `source_gap_over_microbe_ceiling = 10,21 µmol/h`. The source gap exceeds
  the microbial capacity ceiling, so microbial production cannot carry it. What can? Calculate what
  order of magnitude each remaining pathway must have, as the previous round did for clearance.
- **And explain the recovery:** at clearance 109, the direct pathway gives 0,99999 in 24-hour
  recovery against held 0,96. A value of 0,99999 where the measurement says 0,96 means that something
  in the model loses nothing at all. Where do the 4 percent go in reality?
- **Falsifier:** if the source’s own marginal cannot be justified as a population distribution,
  the 17,2 % statement also fails and the entire composition is then diagnostic and nothing more.

# New round (coordinator, 3/10 00:35)
The dietary pathway is the answer and must now be tested as a hypothesis, not assumed. A portion of cod gives 528,9 mg against the requirement 18,4 mg per day, i.e. 28,7 times margin. Changed operation: can the chain reproduce the held median with a REALISTIC dietary distribution instead of our endogenous source term? It requires an intake distribution in the population, not one portion. If yes, the model’s error is identified as a missing PATHWAY and not a missing parameter, and that is a structural finding. And the base-rate test must be corrected: the marginal witnesses giving 17,2 percent are not the source’s patients, so even the stricter distribution rests on an unmatched population. Report it. Falsifier: if a realistic dietary distribution gives a wider spread than the held IQR, the dietary pathway is too coarse an explanation.

# Round 17 (coordinator, 3/10 06:50) — WRITE THE PORT, the chain of rejections IS the result
Sixteen rounds have built a chain of rejections, each parameter-free or nearly so, and
it is time to write it as a result instead of adding a seventeenth on top.

The chain, as it stands: starting point 18,0 µM against held median 3,7, i.e. 4,87× above. The error turned out
to be TWO errors in opposite directions — source 4,13× too low and renal clearance 26,16× too low — and with
both published values the chain lands at 3,3958 µM, inside the held interquartile range. But the base-rate test
showed that 70 % of the published rectangle hits the same range, so the hit carries information only
under the source’s own marginal, and those marginal witnesses are not the source’s patients. The dietary pathway explains
the source gap with a large margin — one portion of cod covers the daily requirement 29 times over. Urine data require
46,5 % of the dose to lie in a non-plasma pool at six hours while storage can carry at most 8 %. And now:
a source-only explanation is rejected by a factor of 1,71 and an additive apical one by 1,36, without a single new
fit.

Changed operation: write PORT.json as that chain, in order, with each rejection’s factor and its
status as parameter-free or conditional. It is the project’s most complete example of how a gap
is dissected, and it has value as a method example regardless of whether the final answer is reached.

And name the only measurement remaining: `actual_matched_TMAO_pairs = 0` through sixteen rounds.
Paired measurements in the same individuals are what the chain has lacked all along. Formulate it as an entry in
notes/ACQUISITION_TARGETS.json format — node ID, quantity, unit, what changes — so it can be acquired by
the data acquisition channel.

Falsifier: if the chain cannot be written without some link turning out to rest on a fit I
believe is parameter-free, say which. A chain with a hidden fit in the middle is worse than no chain.
