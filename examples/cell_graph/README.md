# Cell graph
A claims graph over 35 of the physiology cells in `src/bodytwin/cells/` (the first public set; the cells added later are not catalogued here): one node per cell module (id = module name, cluster = its subsystem package), plus one aggregate node per subsystem that depends on the cells beneath it.
Node status comes from each cell's own selftest run: PROVEN where the run prints an overall-pass token, ASSUMED where it prints a mixed summary, REFUTED where the claimed number is not reproduced, OPEN for the shared library that has no gates.
Edges follow one mechanical rule, stated in full at the top of `build_cell_graph.py`: module imports first, then shared-quantity edges (ATP, O2, membrane potential, Ca) from each consuming cell to the cell that measures that quantity, ordered by supply direction so the graph stays acyclic.
`designs/cell_designs.json` is the catalog the framework builder reads; `ANCHOR_GRAPH.json` is what it writes through the guarded writer (its companion inventory markdown targets a different cluster vocabulary and is written to a scratch path instead).
`FOLD_LEDGER.jsonl` holds one entry per cell run: the decisive number, the command that reproduces it relative to the repository root, and the nodes that consume it.
Counts: 45 nodes (35 cells, 10 subsystems), 44 edges (9 cell-to-cell, 35 aggregate).
Rebuild from the repository root: `.venv/bin/python examples/cell_graph/build_cell_graph.py`.
Check: `.venv/bin/python -m pytest -q tests/test_cell_graph.py`.
Schema check with the graph engine: `anchor_graph_tools.py validate --graph examples/cell_graph/ANCHOR_GRAPH.json --ledger examples/cell_graph/FOLD_LEDGER.jsonl`.
