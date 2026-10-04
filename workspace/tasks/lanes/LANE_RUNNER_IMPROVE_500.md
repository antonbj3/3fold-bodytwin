# Assignment: improve all 500 mutations (M1–M500) — more value, more innovation potential

You're a lane runner. `lane-model`Agent with full mandate from the coordinator. (Anton). Free hand to improve, but do not touch the evidence integrity. Perform without asking; account for your standard choices.

## Goal
Go through **all 500 mutations** and improve them so they become **more value-creating** and create **greater innovation potential** — without destroying traceability, numbering or evidence.

## What "better" means (concretely)
For each mutation, sharpen, do not blur:
1. **Decision & observability:** which decision changes if the mechanism holds; measurable quantity, units, time window, regime, and who consumes the result.
2. **Discriminating test:** must distinguish the mutation from its outperformed baseline at matched (setup+update) cost. Make the test sharper and cheaper where possible.
3. **Evaluation target:** raise the bar so it beats the strongest known baseline at matched accuracy, not a weak one.
4. **Innovation potential:** connect to the graph (contradiction engine T43, spectral unity T41, method transfer T44, operator composition T45, emergent combination T48) and to new decision boundaries. Prioritize mutations that unlock downstream capability or remove a decisive uncertainty.
5. **First principles:** every mutation should have a decomposition according to `FIRST_PRINCIPLES.md` (constituents, interactions, governing relations, leaf status: `DERIVED_UNDER_ASSUMPTIONS` / `EXTERNALLY_MEASURED` / `CONSTITUTIVE_CLOSURE` / `UNKNOWN`, stopping argument, discriminating test).
6. **Remove redundancy:** merge near duplicates into distinct mechanisms; one well-formed mutation > several almost identical ones.

## Inputs (read everything first)
- Mutations (5 batches, M1–M500): `external_research_path`, `_v2.md`, `_v3.md`, `_v4.md`, `_v5.md`.
- Canonical program: `tasks/free48/SEED_PROGRAM_500_20260929.json`.
- Run context: `.../tasks/free48/INITIAL_JOBS.json`, `STATE.json`, `CATALOG.json`; job directories `BT-FW48-SEED-001..500` in `external_mount` (each job has `BRIEF.md`, `inputs/`, `JOB.json`).
- Rules: `START.md` and `local_path` plus `external_research_path`.
- Graph: use the project's own `./graph working rank` and `./graph working packet --id <ID>` to ground the improvements before choosing/reshaping mutations (see `notes/GRAPH_WORKFLOW.md`).

## Freedom and boundaries
- **Keep the numbering M1–M500** so results keep mapping to the same job ID (`BT-FW48-SEED-###`). Improve the content; do not rename.
- **Do not touch running/completed jobs.** Update `BRIEF.md` only for jobs that **do not** have `.ovh_claim`, **are not** running and **do not** have `RESULTS.md`. In-flight and completed jobs are left untouched.
- Put new/updated definitions in **new files** (see Output). Do not change `GRAPH.json`/generated files by hand.
- No automatic evidence admission; retain negative results and their validity limits. No invented data or citations.
- No mail, pushes, publication, credential changes.

## Output
1. `external_research_path` — all 500 improved, M1–M500, with explicit "what was improved and why" per mutation.
2. `tasks/free48/SEED_PROGRAM_500_IMPROVED_20260929.json` — machine-readable, same 500 IDs.
3. `external_research_path` — the 40 with most value creation/innovation potential, ranked with justification.
4. For unstarted jobs: updated `BRIEF.md` (only within the boundary above), and note which were updated.
5. `.../tasks/lanes/LANE_RUNNER_IMPROVE_500_RESULT.md` — short summary: method, count changed, redundancy merged, estimated value increase, what was NOT changed and why.

## Working method
- Start by reading all 5 batches + JSON and build a map (theme T1–T50, domain, baseline, test).
- Work in batches (100 at a time), save intermediate results continuously.
- Stay within reasonable compute time; this is design/review, not new physical measurements.

## Definition of done
- All 500 exist in v6 + JSON, with sharpened test/comparative improvement/innovation connection and first-principles decomposition.
- TOP40 ranking complete.
- Redundancy consolidated and reported.
- No running/completed SEED jobs affected; no evidence integrity broken.
