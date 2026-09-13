"""Fold-layer tests: graph build, guarded write, hardening check, skip rules, extraction."""
import json
import os
import shutil
import subprocess
import sys

import pytest

from conftest import EXAMPLES, FRAMEWORK, framework_env

import bodytwin_graph_io as gio
import bodytwin_fold_skip_rules as skips
import mt_extract
import mt_id_alias

DESIGNS = EXAMPLES / "mini_graph" / "designs" / "cell_designs.json"


def build_graph(tmp_path):
    env = framework_env(tmp_path)
    env["BODYTWIN_DESIGNS"] = str(DESIGNS)
    proc = subprocess.run([sys.executable, str(FRAMEWORK / "bodytwin_build_anchor_graph.py")],
                          env=env, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    return env, json.loads((tmp_path / "ANCHOR_GRAPH.json").read_text())


def test_build_anchor_graph(tmp_path):
    _, graph = build_graph(tmp_path)
    assert set(graph) >= {"_meta", "nodes"}
    ids = [n["id"] for n in graph["nodes"]]
    assert len(ids) == len(set(ids)) == 2
    assert all(n["status"] == "OPEN" for n in graph["nodes"])
    assert all(n["cert_design"]["legs"] for n in graph["nodes"])


def test_hardening_check_passes_on_the_example_graph(tmp_path):
    env, _ = build_graph(tmp_path)
    proc = subprocess.run([sys.executable, str(FRAMEWORK / "bodytwin_hardening_check.py")],
                          env=env, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PASS" in proc.stdout


def test_hardening_check_fails_on_a_duplicate_id(tmp_path):
    env, graph = build_graph(tmp_path)
    graph["nodes"].append(dict(graph["nodes"][0]))
    (tmp_path / "ANCHOR_GRAPH.json").write_text(json.dumps(graph))
    proc = subprocess.run([sys.executable, str(FRAMEWORK / "bodytwin_hardening_check.py")],
                          env=env, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 1
    assert "duplicate" in proc.stdout.lower()


def test_scorecard_runs(tmp_path):
    env, _ = build_graph(tmp_path)
    proc = subprocess.run([sys.executable, str(FRAMEWORK / "bodytwin_scorecard.py")],
                          env=env, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0
    assert "total cells: 2" in proc.stdout


def test_guarded_writer_is_atomic_and_verifies(tmp_path):
    path = tmp_path / "ANCHOR_GRAPH.json"
    graph = {"_meta": {"name": "test"}, "nodes": [{"id": "A"}, {"id": "B"}]}
    gio.write_graph(graph, str(path))
    reread = json.loads(path.read_text())
    assert [n["id"] for n in reread["nodes"]] == ["A", "B"]
    assert not list(tmp_path.glob("*.tmp*"))


def test_extract_cell_json_roundtrip(tmp_path):
    payload = {"proposed_cell": {"cell_id": "MODEL-TEST-CELL", "claim": "x"}}
    raw = tmp_path / "raw.txt"
    raw.write_text("prose\nMT-JSON-BEGIN " + json.dumps(payload) + " MT-JSON-END\nprose\n")
    got = mt_extract.extract_cell(str(raw))
    assert got["proposed_cell"]["id"] == "MODEL-TEST-CELL"


def test_id_alias_normalisation():
    for key in mt_id_alias.ID_KEYS:
        payload = {"proposed_cell": {key: "MODEL-X"}}
        assert mt_id_alias.declared_id(mt_id_alias.normalize_id_alias(payload)) == "MODEL-X"


def test_skip_rules_flag_a_design_stub_and_pass_a_real_cell():
    stub = {"id": "MODEL-STUB", "claim": "Build a cartilage thickness map", "cert_design": {},
            "data_sources": []}
    assert skips.skip_reason(stub) == "design-stub"
    real = {"proposed_cell": {"id": "MODEL-REAL",
                               "claim": "GHK with measured concentrations reproduces a resting potential of -61 mV",
                               "cert_design": {"legs": [{"name": "GHK"}, {"name": "pump"}],
                                                "anchor": "held-out measurement"}},
            "decorrelated_anchor": "an independently measured potential",
            "honest_gaps": "intracellular concentrations are preparation-dependent",
            "data_sources": ["published electrophysiology"]}
    assert skips.skip_reason(real) is None


def test_fold_module_imports_without_the_private_gate():
    import bodytwin_fold
    assert hasattr(bodytwin_fold, "extract_cell_json")
    assert hasattr(bodytwin_fold, "fold_one")
