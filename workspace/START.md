# BodyTwin

Arbetsmapp: `.`.

**Context check 2026-09-23:** Read [the check of existing work](notes/CONTEXT_COVERAGE_AUDIT_20260923.md) before new experiments. The working view's nine selected packets do not cover the private implementation. `MERGED_GRAPH.json` is an index; code references and previous corrections are in the detailed graph `GRAPH.json`. New coordinator lanes wait until the plan is ready according to Anton.

**Research work 2026-09-23:** use `./graph working rank` and [the graph workflow](notes/GRAPH_WORKFLOW.md) for separate the collaborator/BodyTwin targets, new results, selection of the next lane and feedback. The older `./graph rank` below concerns the imported source structure.

The private BodyTwin graph, its source records and inherited measurement logs are loaded.

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

`references/john_context.md` contains verified contexts about the collaborative Rasmussen and his questions. It is a part of the general BodyTwin orientation.

For a new session: [complete startup message — BodyTwin and the collaborator, first 24 hours](STARTUP_MESSAGE.md).

The session has two parallel tracks: the collaborator's actual needs/software and a relevant geometry demonstration, and broad improvement and innovation hunting in BodyTwin. The first day must provide executable results where possible, tested hypotheses and a basis for further prioritisation. Materials, tissues, computational cells and geometric relationships are part of the broader direction.

The broader direction also includes [researcher roles, multimodality and combined interventions](notes/RESEARCH_USERS_AND_MODALITIES.md), and [Mechanism, movement, video, splats and camera twin](notes/MECHANISM_MODALITY_SOURCES.md). Relevant public data and tools may be acquired according to the startup message. Trial registration and Wine experiments are postponed until the discussion with Anton tomorrow.

## Added research clusters — 2026-09-26

[First-principles mitochondrial stress and immunity clusters](notes/PHYSIOLOGY_CLUSTERS_20260926/README.md): 16 mitochondrial/stress questions, 24 immune questions and 4 coupling/audit tasks. Dedicated working targets: `BT-RESEARCH-MITOSTRESS` and `BT-RESEARCH-IMMUNITY`. Actual source code, prior negative results and primary references are supplied through `tasks/free48` catalog keys `MITOSTRESS` and `IMMUNITY`; free-model planners can expand them. Graph proposal and definition packets are not scientific admission. `ACTIVATION.json` records queued IDs; process state must be checked afresh.
