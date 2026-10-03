# BodyTwin and dental — shared private graph view

Two workspaces use the same version-locked Graph Engine and the same reading contract:

| Workspace | Anchor graph | Other layers |
|---|---|---|
| [BodyTwin](bodytwin/START.md) | 3 950 nodes, 543 dependencies | source registries, explicit quotations and preserved inherited logs |
| [Dental: BodyTwin + manufacturing](dental/START.md) | 4 047 nodes, 624 dependencies | three constraint nets, stress map, all configured measurement logs and sources |

**The anchor count describes goals and claims.** Variables, constraint relations, fold records, file artifacts and citation records are counted separately. A larger total count does not mean more verified results.

## Start working

```bash
cd local_path
./graph status
./graph show --id manufacturing:VEH-CONSTRAINT-NET --limit 10
./graph show --id manufacturing:constraint:vehicle:edge:V-E16-weight-vs-overdimensioning
./graph unbound --limit 5
./graph rank --limit 10
```

The same commands exist in `bodytwin`. `rank` is the existing engine's structural work order based on inherited dependencies, goal weights and costs or cost proxies; it is not calibrated information value or approval of the sources' claims.

## Filerna som visar grafen

| File | Contents |
|---|---|
| `MERGED_GRAPH.json` | shared index with typed nodes and links across all imported layers |
| `GRAPH.json` | lossless working view of the anchor graphs, with namespaces and recalculated load |
| `RELATION_GRAPH.json` | variables, 59 constraint relations in dental, raw margins and stress comparison |
| `MEASUREMENTS.json` | fold records, all original occurrences, exact bindings, corrections and ID conflicts |
| `SOURCES_VIEW.json` | source registries and explicit citation records with provenance |
| `EVIDENCE_BINDINGS.json` | source files, hash checks, resolved JSON keys and comparison of cited values |
| `ENGINE_MARGIN_VIEW.json` | the actual engine's reading of constraints; statistical NO-DATA where supporting data are missing |
| `OPEN_INTEGRATION_ITEMS.json` | unbound references, ambiguous IDs, differences and information needed for numerical fusion |
| `SOURCE_AUDIT.json`, `VALIDATION.json` | older evidence gaps and the import's checks and separate layer counts, respectively |
| `PROFILE.json`, `PROFILE_CONFLICTS.json` | engine profile merged with existing `EngineProfile.merge`; default values are uncalibrated |

## Verified contents in dental

- Three nets: vehicle, projector and factory. 48 variable names (28 explicitly declared, 20 only named endpoints), 59 constraints, of which 10 cover more than two variables.
- 13 measurement log files: 3 799 occurrences become 2 772 unique JSON records after merging 1 027 identical copies. All original rows, files and line numbers are preserved.
- 1 116 source registry records. BodyTwin's 18 397 explicit source references are represented by 16 166 unique citation records. These are neither 16 166 independent measurements nor a counted number of articles.
- 170 explicit references from constraints have been examined. 34 distinct file contents are frozen; 56 cited values match, five differ, 53 file references lack a resolved file. Other references are preserved with their more specific resolution statuses.
- One fold-ID conflict and three differences against the older stress map are exposed. 2 345 fold records lack exact binding to any imported anchor or constraint node; they are searchable as unbound records.

## Readiness and limits

**The layers are imported and connected in a shared reading view.** Original statuses are preserved. No source, dependency edge or measurement has become more certain through the import.

The 26 numerical source margins have different definitions and lack explicit σ. Automatic statistical fusion is therefore not complete: raw values are shown, but no invented σ or independence is used. NO-DATA shows `null` for probabilities and estimates. Shared variables and sources across domain boundaries need reviewed identities and validity domains; the same word is not enough.

The manufacturing anchor graph's older evidence check passes. BodyTwin has 4 793 evidence references that the engine's older validator does not recognize and three decided nodes without evidence. Many references are free text. These are reported format/traceability gaps, not a count of false claims. Existing coupling notes and statuses remain for review.

## Update from the sources

```bash
./graph refresh
./graph validate
```

A new generation is created and activated atomically. Source files are write-protected for the reader; nothing is written back. Own notes, tasks and results are in `notes/`, `tasks/`, `results/` and are preserved. Changes in generated files stop refresh. Changed reader code must be reviewed and then followed by refresh; validate detects code drift.

The source map is `~/projects/graph_workspace/GRAPH_MAP.md`. `SOURCES.json` and `LEDGERS.json` specify the concrete reading. CS-I is the source for the manufacturing graph and nets; other worktrees contribute measurement logs. Seed→commit index, machine telemetry and superseded CS-D graphs are not included.

The engine is frozen at research commit `73e76dd83a601ddf2ccb5bab041512dd426ddb2c`; the canonical public head during preparation was `3928720730b3348cd2c410288f023363a23ada73`. Exact code hashes are in `ENGINE.json` and the generation manifest. The graph session's checkout is untouched.

See [FORMAT.md](FORMAT.md) for the contract. The workspace is private, without publishing or mailing functions.
