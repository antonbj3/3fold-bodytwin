# CX-SURGERYFE — Parametric osteotomy to tissue risk with honest uncertainty

Based on `results/CX-WHATIF2` (195/216 KKT approved points, shape/active mass gap), `results/BT-N28` (pressure model), `results/BT-N33` (inertial), `results/BT-C1b` (UNKNOWN refusal) and `results/CX-PARAM-DENTAL` (transfer of parameter→design, proposal only). This is a mechanical what-if, not a clinical treatment recommendation.

## Task
1. Select three approved CCD/AV cases and construct anatomically consistent before/after bone geometries, with an actual HJC/lever arm/wrap update. Reject 21 KKT gaps; do not fill them with interpolation as a reference.
2. Run the load through regional field/FE or documented proxy with density, material and contact tape. Calculate separate force, surface pressure and maximum local stress. Describe if external tissue geometry/material measurement is missing.
3. Compute a Pareto front with form, load and material insecurity, as well as certificates that stop undefined density, wrong frame or invalid mask. Compare against dental chain traceability scheme as method transfer.

## PREREG and gate
≥3 forms × at least 5 accepted parameter points; all force calculations 141/141 KKT ≤1e-10; FE/proxy solution converges within 5% at halved cell size; at least one robust improvement ≥5% in both force and stress outside the 95% uncertainty band for an actual claim. Counter test: zero change, person-changed geometry, constant-radius pressure and B24/weight × OrthoLoad. If the last gate is not passed: UNKNOWN for planning benefit.

## Gemensamma regler
- Skriv endast under `results/CX-SURGERYFE/`, samt stora mellanresultat i `external_media`Read it. `tasks/NIGHT_PREAMBLE.md`, and relevant A-rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first computation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, countertests and error definition. Document `Builds on` with graph node and source files plus `Not redone`. Preserve negative results and `UNKNOWN`.
- Provide `results/CX-SURGERYFE/RESULTS.md` with first line `# CX-SURGERYFE`, `results.json`, executable code, provenance/hashes and meaningful checks. Report both number of valid and lapsed units.
