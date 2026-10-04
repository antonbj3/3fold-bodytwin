# Steering LANE_CAPACITY_FROM_TMAO — round 1 (the coordinator, 2/10 20:50)

New lane. It takes the slot from LANE_ROBOT_PORT_CONNECT, which closed after five rounds with the conclusion that the robot branch's computation term is not the limit.

- **Obstacle:** the capacity chain has asked about metformin clearance for months, which a certificate today shows is the only one of three substances whose extraction ratio (E = 0,698) makes capacity unbounded above.
- **Changed operation:** change substance. Bound capacity with TMAO (E = 0,368, separation succeeds) from a PAIRED measurement, and check against a thermodynamic ceiling that does not go through flow at all.
- **Strongest control:** the metformin chain that gave 243×. The gain should be a finite interval where it gave an unbounded one.
- **Falsifier:** if the TMAO bound also becomes unbounded after propagated spread, separation does not hold in practice for our anchor — say so, that would be a negative about a result we just booked.

## Start with the regime, because that was today's most costly error
Our largest claim today failed because a requirement was calculated as flux = k_cat × N, thus saturated regime, and compared against a measured value in the linear regime where clearance = V_max/K_m. The factor between them is K_m/C, and for metformin it is over a hundred. **Write out, for every link, whether C lies far below, near or above K_m before comparing anything.** Do that first; everything else rests on it.

## The paired anchor is unusual and should be treated as such
219 ± 78 mL/min TMAO against 119 ± 21 creatinine and 55 ± 14 urea in the SAME individuals. Pooled versus paired refuted two conclusions today — a gap of 2,9–9,1× turned out to be 2,17–7,76× cohort artifact. Here the pairing exists, so use it and say what it is worth compared with pooling.

## If a path is missing where you run
`local_path` does not exist on cloud hosts. Four jobs today silently reimplemented missing modules or reconstructed missing constants, and one of them produced an incorrect headline that I amplified further. If something is missing: write `missing_prerequisite` in the FIRST paragraph of RESULTS.md, not in a footnote.

Allt PENDING_INDEPENDENT_REVIEW. Inga interna data, inga patientdata.

# Round 2 (the coordinator, 2/10 23:00) — the finite bound exists, qualify it

## Review of round 1: you delivered what the lane existed for
`effective_H_point = 346,68 mL/min` and `conditional_H_1SD = [166,09; 1225,52]` with `upper_unbounded: False` — thus a **FINITE** capacity bound, where the metformin chain gave unbounded. At 2SD it becomes unbounded above, however. Gate FAIL is right because the native turnover interval is null, but the finite 1SD interval is the result.

## The task
1. **Say exactly which condition makes the bound finite.** It holds at 1SD and fails at 2SD. What in the spread reverses it, and is 1SD defensible for this anchor or is it a convenience? Answer with propagated uncertainty, not with a choice.
2. **Search the index first.** `notes/OLD_CELL_INDEX.json` has 1309 completed cells with decisive numbers, 696 MEASURED. Run `python3 tasks/build_night/index_old_cells.py --quantity clearance` and `--quantity flow` before deriving anything new. If the quantity is measured there, consume it and cite the cell ID.
3. **Your obstacle names the next step:** "scalar conductance + thermodynamic ceiling do not identify serial rates". Two serial steps cannot be separated by a scalar conductance. Which observable separates them? That is the question, and the thermodynamic bound was meant to be the independent control — say whether it can become that.

**Forbidden:** choosing spread level according to the answer it gives; reporting 346,68 as a number without its condition.

# Round 3 (the coordinator, 2/10 22:15) — FINAL: name the observable or close
You refuted your own result: `one_SD_empirical_support_defensible = False`. The finite bound does not stand. That was the right answer to my question.
- **Obstacle:** no clearance-based route gives a defensible finite capacity bound, because D_2SD = −16,75 and worst_new_person95 = −48,20.
- **Changed operation:** last round. Name the observable that WOULD give a defensible finite bound, with quantity, unit and preparation — concentration ratio at steady state is the candidate because it is set by thermodynamics and not by delivery. Deliver PORT.json and the measurement specification, not a new estimate.
- **Falsifier:** if no observable can give a finite bound, capacity is structurally unidentifiable from what we can measure, and that must stand as the lane's conclusion.
Forbidden: choosing spread level according to the answer; reporting 346,68 without it being refuted.
