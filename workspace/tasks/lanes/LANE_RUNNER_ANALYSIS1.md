# Analysis-1: overlap, common maps, density — and 300-vs-deep

You're the lane runner `lane-model` with full mandate. You are one of three OBEROENDE analytiker av samma material (the other two you don't see)Write your own assessment.

## Goal
Read the four outputs below and analyze **overlap and common maps**: which maps/operators are REUSED between the boundary mutations, the diseases and the maps library. Calculate an honest density (reuse proposal, not observed run). Then answer the question: **does the whole 300+ need new objects, or should it deepen on the existing ones?**

## Input (read all)
- `~/research/BODYTWIN_BOUNDARY_REPORT.md` + `~/research/BODYTWIN_BOUNDARY_MUTATIONS.json`
- `~/research/BODYTWIN_DISEASE_LIBRARY.md` (+ ev. JSON-bibliotek)
- `~/research/BODYTWIN_MAPS_OPERATORS.md` + `~/research/BODYTWIN_MAPS_OPERATORS.json`
- `~/research/BODYTWIN_MUTATIONS_IMPROVED_v6.md` (de 500) som referens
- `~/research/ORGANIZATION_AUDIT.md`
- `tasks/lanes/LANE_RUNNER_BOUNDARY_RESULT.md`, `LANE_RUNNER_DISEASE_RESULT.md`, `LANE_RUNNER_MAPS_RESULT.md`

## To produce
1. **Reuse map:** which maps are reused over ≥2 diseases/mutations. Count honestly; distinguish between "proposed reuse" and "observed".
2. **Density measure:** number of unique maps vs number of uses; which families carry the most.
3. **Value assessment:** which of the 12 boundary + 27 maps are the most load-bearing (unlocks the most downstream).
4. **300-vs-deepening:** a clear recommended answer with justification (number ≠ value).
5. **Gap remaining** to proceed.

## Utdata
`~/research/ANALYSIS_1_OVERLAP.md` + machine-readable reuse map (JSON). Avsluta med `tasks/lanes/LANE_RUNNER_ANALYSIS1_RESULT.md`.

## Constraints
No invented data; reuse is suggestion; graph's notes = reference. No email/pushes/credentials/touching SEED.
