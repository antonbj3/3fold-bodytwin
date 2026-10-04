# BodyTwin

Arbetsmapp: ``.

**Context check 2026-09-23:** Read [the check of existing work](notes/CONTEXT_COVERAGE_AUDIT_20260923.md) before new experiments. Arbetsvyn's nine selected packages do not cover the private implementation. `MERGED_GRAPH.json` is an index; code references and previous corrections can be found in the detail graph `GRAPH.json`. New coordinator lanes wait until the plan is ready according to Anton.

**Research work 2026-09-23:** use `./graph working rank` and [graph workflow](notes/GRAPH_WORKFLOW.md) for separate the collaborator/BodyTwin targets, new results, next lane selection and feedback. The older `./graph rank` below applies to the imported source structure.

The private BodyTwin graph, its source records and inherited measurement logs are loaded.

```bash
./graph status
./graph rank --limit 10
./graph unbound --limit 5
```

Start with `MERGED_GRAPH.json` for the entire structure and `VALIDATION.json` for separate layer counts. `./graph show --id <exakt-ID>` shows a record and its connections. [The common README-file](../README.md) describes all detail layers and [FORMAT.md](../FORMAT.md) defines the contract.

- `notes/`, `tasks/`, `results/`: projektets eget fortsatta work.
- `references/current_bodytwin`: current published BodyTwin code and geometry implementation.
- `references/field_engine`: field engine public code.

The sources' own statuses are preserved. Statistical fusion of raw margins requires complementary uncertainty, normalization and provenance. Unbound records, conflicts and older evidence gaps remain explicitly in the review view.

`references/collaborator_context.md` contains verified context about the external collaborator and his questions. It is a partial task within the general BodyTwin orientation.

For a new session: [Complete start message — BodyTwin and the collaborator, first 24 hours](STARTUP_MESSAGE.md).

The session has two parallel tracks: the collaborator's actual behov/software and a relevant geometry demonstration, as well as broad improvement and innovation shunt in BodyTwin. The first day should provide executable results where possible, tested hypotheses and basis for continued prioritization. Materials, tissues, calculation cells and geometric relationships are part of the broader orientation.

The broader orientation also includes [researcher roles, multimodality and combined interventions](notes/RESEARCH_USERS_AND_MODALITIES.md), as well as [Mechanism, movement, video, space and camera twin](notes/MECHANISM_MODALITY_SOURCES.md). Relevant public data and tools may be downloaded according to the start notice. Trial registration and Wine trials are postponed until the discussion with Anton tomorrow.

## Added research clusters — 2026-09-26

[First-principles mitochondrial stress and immunity clusters](notes/PHYSIOLOGY_CLUSTERS_20260926/README.md): 16 mitochondrial/stress questions, 24 immune questions and 4 coupling/audit tasks. Dedicated working targets: `BT-RESEARCH-MITOSTRESS` and `BT-RESEARCH-IMMUNITY`. Actual source code, prior negative results and primary references are supplied through `tasks/free48` catalog keys `MITOSTRESS` and `IMMUNITY`; free-model planners can expand them. Graph proposal and definition packets are not scientific admission. `ACTIVATION.json` records queued IDs; process state must be checked afresh.
