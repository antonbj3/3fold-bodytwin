# Dental — BodyTwin + tillverkning

Arbetsmapp: `local_path`.

**Research work 2026-09-23:** use `./graph working rank` and [the graph workflow](notes/GRAPH_WORKFLOW.md) for dental nodes, new results, selection of the next lane and feedback. The older `./graph rank` below concerns the imported source structure.

The shared read view combines BodyTwin with manufacturing's dependency graph, three constraint nets, stress map, source registry and measurement logs.

```bash
./graph status
./graph rank --limit 10
./graph unbound --limit 5
```

Start with `MERGED_GRAPH.json` for the full structure and `VALIDATION.json` for separate layer counts. `./graph show --id <exakt-ID>` shows a record and its connections. [The shared README file](../README.md) describes all detail layers and [FORMAT.md](../FORMAT.md) defines the contract.

- `notes/`, `tasks/`, `results/`: the project's own continuing work.
- `references/current_bodytwin`: current published BodyTwin code and geometry implementation.
- `references/field_engine`: the field engine's public code.

The sources' own statuses are preserved. Statistical fusion of raw margins requires additional uncertainty, normalisation and provenance. Unbound records, conflicts and older evidence gaps remain explicit in the review view.

`data/corpus` now points to the checksum-verified data copy on `external_mount`. `references/dataset_mapping` contains the collection mapping. The move covered 52 005 files; this does not mean all originally planned dataset downloads are complete.

The project context is possible experimental validation with Osayd Alawawda around dental geometry through mesh → field → simulation → optimization → manufacturing. The exact experiment and prioritization are selected in the next step.

For a new session: [complete startup message with goals and research mandate](STARTUP_MESSAGE.md).

The assignment includes expanding the entire chain as graph nodes and relations: anatomy/tissues, materials, physics, generative design, manufacturing and experiment. The factors shall get an explicit mathematical role in the design problem. This domain expansion is planned continuing research; it is not yet included as new nodes in the imported read view.

Material knowledge is its own innovation track: composition, structure, process, properties and material models. Simulation shall also cover the operations themselves and their course, with linkage between tool, tissue and construction before, during and after relevant steps. Both tracks shall feed knowledge back to the graph and generative design.

[The data collection addendum from anton-02](notes/ADDENDUM_DATA_COLLECTION_2026-09-22.md) preserves the handoff and its checking needs.

[Existing computational cells and local experiments](notes/ADDENDUM_COMPUTE_CELLS_2026-09-22.md) points to published BodyTwin code, private trials, cad-to-simulation's cell registry and the engines' codebases.
