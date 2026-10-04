# Steering — round 1 (the coordinator, 2/10 21:30)

New lane in a Sol wave. Anton 2/10: "it is up to you to fix all this and make it work, you know the goal is comparative improvement". The wave builds couplings; it revises no numbers.

Read the lane brief in `tasks/lanes/` for the task. These are the additions that apply to all lanes in the wave.

## The graph engine is now in the cloud
It is at `/opt/agents/graph_engine` on both cloud hosts, with the cited path symlinked there, pinned commit `352c6d3`. Import via the path is verified. **Do not reimplement it** — four jobs today did so silently because it was missing, and one of them gave a result I had to qualify because its "engine" was its own rewrite.

## Five error classes that cost us today, each with its price
1. **Wrong regime.** A requirement was calculated as flux = k_cat × N, thus saturated, and compared against a value in the linear regime where clearance = V_max/K_m. The factor was over a hundred and it refuted today's largest claim. Say which regime and which denominator every number is in before comparing.
2. **Pooled is not paired.** A gap of 2,9–9,1× turned out to be 2,17–7,76× cohort artifact from amount and activity coming from different individuals.
3. **Structural absence is not independence.** 25 of 27 anchors did not reverse sign under a physiological modifier — because the modifier was not even a port in the model. Mark that untested, never robust.
4. **Prediction is not measurement.** A field with *conditional* in its name was read as measured, by me. Mark every number MEASURED or DERIVED.
5. **A verdict nobody reads is worthless.** Four instances today. Every verdict you produce must have a test that FAILS if a consumer ignores it.

## Missing file
`missing_prerequisite` in the FIRST paragraph of RESULTS.md, not in a footnote. Never silently reconstruct missing input — four jobs did so today and one produced an incorrect headline that I amplified further.

Everything PENDING_INDEPENDENT_REVIEW. No internal data leaves the machine.

# Round 2 (the coordinator, 2/10 21:45) — you refuted my premise, build the cheaper answer

## Review of round 1: TIE, and that was the right answer
I built the lane on the box form being unable to carry our 11,27× case. Your obstacle says: *"Observed synthetic history compresses into scalar internal state; unconsumed/missing state is the defect."* History thus compresses into a scalar internal state — exactly the condition my own criterion gave for when the box suffices. **It suffices.** No new path identity should be built.

Verified by me: `order_work_bound_J_m2` [605,9237662351542; 7494,591721487457] is exactly the numbers from `PORTS_R4_FINAL.json` — no invented values. And with the state declared, CONTRADICTION becomes OK/VIOLATED, with the scalar control at REGIME-BOUNDARY. Five mutants rejected, IGNORE_IDENTITY 5, IGNORE_LEDGER 4, IGNORE_MEMORY 3.

## The task: make the state declared and consumed
The defect is not representation but that the state is neither declared nor read. It is the fifth instance today of the same thing, and this time we know exactly which scalar it is.

1. **Name the scalar and its update rule**, in a form a cell can declare as a port: its name, its unit, how it is updated per event, and what it is initialized to. "Number of previous passages" and "accumulated dose" are candidates — say which our case actually requires.
2. **Show what changes in a consumer that reads it** against one that does not. You already have CONTRADICTION → OK/VIOLATED; make that transition an executable example, not an observation.
3. **Test whether that scalar suffices for the SURGICAL case**, not just for heat opening: "third passage through the same tissue" should differ from the first even when all other declared quantities are equal. If a scalar suffices there too, the question is closed and that is a strong result.
4. **Hand the proposal to the graph lane, do not build in the engine.** They own it and have four verdicts in the shared sign box since tonight.

**Falsifier:** if the surgical case requires more than one scalar, my original premise remains for that particular case — say so, with which part of history does not compress.
**Forbidden:** building a path identity when a scalar suffices; changing the graph engine; reporting a state without a test that fails when it is ignored.

# Round 3 (the coordinator, 2/10 22:45) — FINAL: deliver the port, the question is closed
`ordinal_single_scalar_suffices = True`, `physical_count_only_suffices = False`, and two local modes suffice under the tested Prony law without path identity. You have refuted my premise twice and the answer is complete.
- **Deliver:** PORT.json with the ordinal's name, unit, update rule and initial value, plus the definition of the two modes. Formulate it as a port contract a cell can declare.
- **To the graph lane as a PROPOSAL, not as a change** — they own the engine and have four verdicts in the shared sign box since tonight. Say that an ordinal plus two modes suffice for our path dependence, so they can decide whether `Box` needs anything at all.
- **Keep the 192 counterexamples** as regression cases: identical present state, different future, gap 0,0038–0,0705. It is the evidence that the state is needed and that the difference is small enough to miss.
Forbidden: reopening path identity; reporting without the test that fails when the state is ignored.
