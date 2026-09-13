"""Cell-graph tests: rebuild examples/cell_graph from the cell modules and validate it.

The rebuild is deterministic (it reads the 35 cell modules plus the recorded run outcomes), so
it is run in place and the regenerated files must be byte-identical to the committed ones.
Schema validation uses the graph engine's anchor_graph_tools when that repository is reachable
(GRAPH_ENGINE_TOOLS, or a sibling checkout); otherwise the structural checks below still run.
"""
import json
import os
import pathlib
import subprocess
import sys

import pytest

from conftest import CELLS, EXAMPLES

GRAPH_DIR = EXAMPLES / "cell_graph"
BUILDER = GRAPH_DIR / "build_cell_graph.py"
SUBSYSTEMS = sorted(p.name for p in CELLS.iterdir() if p.is_dir() and (p / "__init__.py").exists())


def _tools_dir():
    env = os.environ.get("GRAPH_ENGINE_TOOLS")
    candidates = [pathlib.Path(env)] if env else []
    repo = pathlib.Path(__file__).resolve().parents[2]
    candidates.append(repo / "3fold-graph-engine" / "src" / "graph_engine" / "tools")
    for c in candidates:
        if (c / "anchor_graph_tools.py").exists():
            return c
    return None



def catalogued_cells():
    """The cells the example catalogues (a fixed subset of src/bodytwin/cells); each must exist on disk."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("build_cell_graph", BUILDER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cells = sorted(r["mod"] for r in mod.CELL_ROWS)
    on_disk = {p.stem for p in CELLS.rglob("*.py") if p.name != "__init__.py"}
    assert set(cells) <= on_disk
    return cells

@pytest.fixture(scope="module")
def rebuilt():
    # Include the persistent lock file in both complete byte snapshots.
    (GRAPH_DIR / ".graph.lock").touch(exist_ok=True)
    before = {p.name: p.read_bytes() for p in GRAPH_DIR.iterdir() if p.is_file()}
    proc = subprocess.run([sys.executable, str(BUILDER)], capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr
    after = {p.name: p.read_bytes() for p in GRAPH_DIR.iterdir() if p.is_file()}
    assert after == before, "rebuild is not reproducible: " + str(
        sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k)))
    return json.loads((GRAPH_DIR / "ANCHOR_GRAPH.json").read_text())


def test_one_node_per_cell_and_subsystem(rebuilt):
    cells = catalogued_cells()
    subsystems = sorted({p.parent.name for p in CELLS.rglob("*.py") if p.stem in cells})
    ids = {n["id"] for n in rebuilt["nodes"]}
    assert set(cells) <= ids
    assert {f"{s}_subsystem" for s in subsystems} <= ids
    assert len(rebuilt["nodes"]) == len(cells) + len(subsystems)
    for n in rebuilt["nodes"]:
        assert n["cluster"] in SUBSYSTEMS


def test_ledger_has_one_entry_per_cell(rebuilt):
    rows = [json.loads(line) for line in
            (GRAPH_DIR / "FOLD_LEDGER.jsonl").read_text().splitlines() if line.strip()]
    cells = catalogued_cells()
    assert sorted(r["node"] for r in rows) == cells
    for r in rows:
        assert r["decisive_number"] and r["gate_cmd"].endswith(r["source"])
        assert (pathlib.Path(__file__).resolve().parents[1] / r["source"]).exists()


def test_graph_is_acyclic_and_loads_match(rebuilt):
    nodes = {n["id"]: n for n in rebuilt["nodes"]}
    load = {i: 0 for i in nodes}
    for n in rebuilt["nodes"]:
        for d in n["depends_on"]:
            assert d in nodes, d
            load[d] += 1
    for n in rebuilt["nodes"]:
        assert n["load"] == load[n["id"]]
    seen, done = set(), set()

    def walk(i):
        assert i not in seen, f"cycle at {i}"
        if i in done:
            return
        seen.add(i)
        for d in nodes[i]["depends_on"]:
            walk(d)
        seen.discard(i)
        done.add(i)

    for i in nodes:
        walk(i)


def test_validate_reports_no_errors(rebuilt):
    tools = _tools_dir()
    if tools is None:
        pytest.skip("graph engine tools not reachable")
    sys.path.insert(0, str(tools))
    import anchor_graph_tools as agt
    report = agt.validate(rebuilt, ledger_path=GRAPH_DIR / "FOLD_LEDGER.jsonl", repo_root=GRAPH_DIR)
    assert report["error_count"] == 0, report["errors"]
    assert agt.next_actions(rebuilt)
    assert agt.priority(rebuilt)
