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

# Round 2 (the coordinator, 2/10 21:45) — scale the content and root the reported

## Review of batch 1: PASS, and the measurement was correct
Ten entries with geometric condition, 6 practitioner reported and 4 measured, and `baseline_questions_with_this_content = 0` — a planner without the entries returns zero. It is the information link check done correctly. Three source cases waived instead of guessed, 21 of 21 adversa cases passed, and the caveat "structural N_eff not empirical independence" is exactly the difference that makes root tracing meaningful.

Anton has said that the surgeon part is important and that it should not be left out. It scales now.

## The task
1. **From ten to as many as the budget will carry, without losing the field requirement.** Each entry should still have tissue, action, geometric condition, expected behavior, locator, and `evidence_class`. An entry without a geometry condition is not counted — the robot only knows its own geometry.
2. **Broaden the tissue families.** You have two. Surgery involves all forms of tissue: skin and its layers, fascia, muscle, tendon, bone, vessel wall, nerve, intestinal wall, liver, brain. Say who you covered and who you didn't reach.
3. **Rottraca the six reported.** `claim_federation` on pinned commit `352c6d3` traces each source to its root sources and effectively counts independent counts. You flagged yourself that your N_eff is structural — make it empirical where the sources allow. A trick-of-the-trade that fifty sources repeat but traces to one observation should carry its N_eff, not its citation count. The scale is what makes this necessary: at thousands of sources, the nominal number is meaningless.
4. **Search the index first.** `notes/OLD_CELL_INDEX.json` has 1309 ready cells, 696 MEASURED and 431 REFUTED. Run `python3 tasks/build_night/index_old_cells.py --quantity <storhet>` before deriving a number. An already measured quantity should be consumed with cell-ID quoted, and an already refuted question should not be rerun.
5. **Contradictions are the most valuable items.** Where a statement of practice conflicts with a measured quantity we hold, one of them is wrong in a context where both are used. Lift them separately.

**Falsifier:** if every new entry that can be formulated with a geometric condition is already a measured quantity in the index, the practice layer adds nothing new and it should be said.
**Forbidden:** entries without a geometric condition; citations counted as independent support; to mix MEASURED and PRACTITIONER_REPORTED in an unlabeled statement.

# Round 3 (the coordinator, 2/10 22:45) — rottraca the seven, and the source is now questionable
33 records, 10 tissue families, 25 MEASURE against 7 practitioner reported. The distribution turned in the right direction.
**Your records are now in a queryable source:** `notes/SURGICAL_EXPERT_SOURCE.json`, queried with `python3 tasks/build_night/surgical_expert.py --quantity <klass> --trust high`. 155 entries total over four blocks, keyed per STORHET because edges are formed between magnitudes and not anatomy. Ask it before you harvest more — duplicates cost money.
- **Hinders:** seven practitioner entries are untraced, and the RSTL genealogy is incomplete, so some entries rest on a tension line map whose lineage is not traced.
- **Changed operation:** root the seven with `claim_federation` (pinned `352c6d3`) and report N_eff against nominal count. Trace RSTL to its primary source — a map that everyone cites but traces to an observation should bear it.
- **Calibration you should know:** another block had **14 of 69 hard-error entries** on resolution pass — ten misattributed quotes and four numbers that were not in the source. Expect ~20 % in memory generated material and verify accordingly.
- **Falsifier:** if RSTL traces to a single observation, every record resting on it should be labeled with that N_eff.

# Round 4 (the coordinator, 2/10 22:40) — same test on the others 26
The result is sharp: **all seven practitioner records have structurally N_eff = 1,0**, and bibliographically
fanout overestimated independence by **factor 5**. It is exactly the calibration the source needed.
- **Hinder:** 26 of 33 entries have not received the same treatment, and `empirical_N_eff_resolved = 0` —
  structurally N_eff is counted, empirically not.
- **Changed operation:** run the same root tracing on the remaining 26 and report the distribution of
  N_eff, not an average value. If the fanout factor 5 his over the entire set, it is a property of
  the literature we harvest from, and then it must be entered into the metadata of the source so that every future question
  see it. And solve the unresolved genealogy plus the RSTL root.
- **Write back in the questionable source:** `notes/SURGICAL_EXPERT_SOURCE.json` carries today by mail
  `trust` but nothing N_eff. An entry with N_eff = 1 and `trust: high` is not the same as five
  independent high records, and whoever asks the source must see the difference. Put N_eff by post.
- **Falsifier:** if the fanout factor varies greatly between tissue families, it is not one
  characteristic of the literature but of our harvest, and then it is the harvesting method that should be reported.
**Check:** nominal number of entries, i.e. 33 against the sum of N_eff.

# Round 5 (the coordinator, 2/10 23:05) — my quintuple fanout fell, so the question is the HARVEST METHOD
The result is uniform and strong: **all 33 entries have N_eff = 1,0**, sum 33,0, and you refuse
read the sum as jointly independent — right. 35 distinct publications, 17 source study families.
And `fivefold_universal_fanout_gate = FAIL` fell my own book five factor as universal; the
was specific to the first seven entries. That fix is yours, not mine, and it's in the log.
- **Obstacles:** `empirical_N_eff_resolved = 0` and an unresolved genealogy. Structurally N_eff is counted
  over the entire set, empirically not for a single item.
- **Changed operation:** 33 of 33 at N_eff = 1 is an outcome that can also be a property of
  THE HARVEST and not with literature. Try it: if every record was picked up by a specialist lens
  asked to name a threshold, each entry per construct comes with a source. Measure it — how many of
  the 17 source study families have MER than one independent study behind the same magnitude, and how many
  record could have gotten N_eff > 1 if harvest looked for replicates instead of thresholds?
  That number determines whether 33-of-33 sowns anything about the world or about us.
- **The next harvest round should search REPLIKAT:** for the magnitude classes that have the most records
  (clearance_distance 44, pressure 24, time 21, flow 16, force 11 in `SURGICAL_EXPERT_SOURCE.json`),
  seek a second independent study for the same magnitude and same tissue family. An entry with N_eff = 2 is
  host more than five new with N_eff = 1.
- **Falsifier:** if none of the 17 families has a second independent study for the same magnitude, is
  the whole area single-sourced and it is a heavy finding about the tissue literature — then report it
  as such, not as a crop failure.
**Check:** nominal 33 entries against the sum of N_eff, so 33 against 33 today.

# Round 6 (the coordinator, 2/10 23:35) — enter the N_eff mode in the source and stop harvesting
15 of 17 families lack a verified independent comparison in your search, 0 exact protocol replicates,
and you put `negative_is_world_absence = False` — right, and better discipline than myself kept four
times tonight.
- **The conclusion is now strong enough to stop harvesting.** The source can RIKTA a measurement, LIMIT
  a range and CONTRADICT a cell, but not certify a decision. That is the limit, and it must stand
  in the source's metadata so any future questioner sees it without reading a log.
- **Changed operation:** finish writing `notes/SURGICAL_EXPERT_SOURCE.json` as a consumable resource.
  Four defects are known and must be fixed: (1) 50 of 188 entries have `unit = None` because the unit
  is embedded in `limiting_quantity`; (2) metrics lie as JSON strings inside `threshold` in
  instead of structured fields with value, SD and n; (3) magnitude class is derived from free text and
  fail — a cufflink entry was classified as `time` because the threshold text contained "time to
  onset of pain", which is an error in MIN `derive_class` in `tasks/build_night/surgical_expert.py` where
  time is tried before pressure; (4) locator notes have been copied between entries. Fix (1), (2) and (4) in
  JSON. For (3): suggest the fix but don't change the script — it's mine, I'll take it.
- **Add the fields a cell needs to read a record without a human in the middle:**
  `quantity_role` (parameter, boundary condition, validation target or validity condition), `maps_to` with
  cell and field, `comparator`, `applicable_range`, and species/temperature/ex-vivo.
- **Forger:** if no record can be provided with a `maps_to` pointing to a field that actually
  is in an executable cell, the source is not an information source but a reading list, and it should be said.
**Forbidden:** more entries before the schedule is consumable; to raise a post's trust because it received
a N_eff.

# New round (coordinator, 3/10 00:35)
The schedule is now consumable and the call is fired correctly. The obstacle is the 155 records without value: a record without value can neither be consumed as a parameter nor falsify a cell, so they are dead weight in an information call. Second operation: for the 155, determine who HAS a food item in their locator but when it has not been extracted, against who never had one. The first is extraction work, the second is a property of the cold material. Report the two numbers. And connect at least three entries to a field in a crossable cell with maps_to, said call proves it can be consumed. Forgers: if none of the 155 have an extractable food item, the older half of the summon is a load list and not an information summon, and it should be said outright.

# New round (the coordinator, 3/10 02:05)
The falsifier fell in favor of the source: 82 of 155 older records are primary numerical observations,
thus 52,9 percent, plus 90 with positive source measurement and 97 verified numerical lines. The older one
half is not a reading list and the work is extraction, not pruning.

But after nine rounds, the big obstacle remains unchanged: native_predictive_bindings = 0. The source is
well-ordered, queryable, N_eff-labeled, and correctly abstaining from 940 questions — and no cell consumes a single
mail. A source of information that no one asks is not a source of information.

Changed operation, and it's the only one that counts now: get EN record consumed by EN executable cell, whole
the road. Select from the 82 primary numeric the entry whose magnitude matches a field actually present in
tasks/free48/sources/SURG_HEMOSTASIS/ or SURG_COLLAGEN/, type maps_to, run the cell with the entry as
input, and report what moved. If nothing moves, the post is not load-bearing and the next one should
be chosen.

The 58 with unresolved measurement status and the 6 which are protocols and not outcomes should remain as they are. To
sorting them further is cheaper than connecting a record, and therefore the wrong priority.

Forger: if none of the 82 has a magnitude that matches a field in any executable cell, the source and
the cells keyed to different magnitudes, and then that's the key problem to solve — not more records.
It would be a heavy and useful find after nine rounds.
