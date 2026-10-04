# Analysis-2: contradictions, gaps and evidence integrity

You're the lane runner `lane-model` with full mandate. You are one of three OBEROENDE analytiker av samma material (the other two you don't see)Write your own assessment.

## Goal
Read the four outputs and actively look for **contradictions and gaps**: where do boundary/diseases/maps contradict each other, where are claims unfounded, where are there silent assumptions, and where is evidence integrity broken (definition objects that look like results)?

## Input (read all)
- `external_research_path` + `..._BOUNDARY_MUTATIONS.json`
- `external_research_path`
- `external_research_path` + `.json`
- `external_research_path`
- `external_research_path`
- `tasks/lanes/LANE_RUNNER_*_RESULT.md`

## To produce
1. **Contradiction list:** exactly where A/B/C clash (mechanisms, ports, spaces, priorities), with quote/reference.
2. **Unfounded claims:** where more is claimed than the material supports (especially "solves disease", "map exists", coverage).
3. **Evidence integrity:** verify that everything is `NOT_ADMITTED`/definitions and that no measured data have been fabricated.
4. **Gaps against what you check in the graph yourself:** run `./graph working rank/packet` for spot checks; do claims match the graph's own notes?
5. **Ranked risks** for further development.

## Utdata
`external_research_path` + JSON. Avsluta med `tasks/lanes/LANE_RUNNER_ANALYSIS2_RESULT.md`.

## Constraints
No fabricated data; the graph's own notes = ground truth. No emails/pushes/credentials/SEED touching.
