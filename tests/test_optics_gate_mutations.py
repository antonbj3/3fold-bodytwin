"""Mutation checks for the optics probe gates that used to be literal or vacuous
(from the 2026-09 review).

Each check evaluates the gate expression exactly as committed in the probe (extracted
with ``ast``, so the probe's ``main`` and its report writes are not run) on a correct
input and on mutated inputs, and asserts that the gate can turn false. No GPU or
Modal call is made; photon_modal_suite is never imported.
"""
import ast
import importlib.util

import numpy as np
import pytest

from conftest import REPO

OPTICS = REPO / "probes" / "optics"


def _source(name):
    return (OPTICS / name).read_text(encoding="utf-8")


def _gate_keyword(tree, key):
    """The value expression of ``key=...`` inside the ``gates=dict(...)`` call."""
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "gates"
                and isinstance(node.value, ast.Call)
                and getattr(node.value.func, "id", None) == "dict"):
            for kw in node.value.keywords:
                if kw.arg == key:
                    return kw.value
    raise AssertionError(f"gate {key!r} not found")


def _eval(expr, namespace):
    code = compile(ast.Expression(body=expr), "<gate>", "eval")
    return eval(code, dict(namespace))  # noqa: S307 - expression taken from the repo source


def _load(name):
    path = OPTICS / name
    spec = importlib.util.spec_from_file_location("gate_mutation_" + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------ tetra_ray_grid_probe ----
@pytest.fixture(scope="module")
def grid():
    probe = _load("tetra_ray_grid_probe.py")
    v, t, labels = probe.fixture(2)
    _, arrays = probe.observe(v, t, labels)
    return probe, v, t, labels, arrays


def test_partition_gate_uses_the_predicate_not_a_literal():
    expr = _gate_keyword(ast.parse(_source("tetra_ray_grid_probe.py")), "partition_accepted")
    assert not isinstance(expr, ast.Constant)
    assert "conforming_partition" in ast.unparse(expr)


def test_partition_predicate_accepts_the_conforming_fixture(grid):
    probe, v, t, labels, arrays = grid
    assert probe.conforming_partition(v, t, labels, arrays) is True


def test_partition_predicate_rejects_a_missing_cell(grid):
    probe, v, t, labels, arrays = grid
    assert probe.conforming_partition(v, t[1:], labels[1:], arrays) is False


def test_partition_predicate_rejects_a_label_count_mismatch(grid):
    probe, v, t, labels, arrays = grid
    assert probe.conforming_partition(v, t, labels[:-1], arrays) is False


def test_partition_predicate_rejects_a_gap_in_a_ray_sequence(grid):
    probe, v, t, labels, arrays = grid
    mutated = {k: a.copy() for k, a in arrays.items()}
    offsets = mutated["offsets"]
    ray = next(i for i in range(len(offsets) - 1) if offsets[i + 1] - offsets[i] > 1)
    mutated["intervals"][offsets[ray] + 1, 1] += 1e-3
    assert probe.conforming_partition(v, t, labels, mutated) is False


def test_partition_predicate_rejects_a_repeated_cell(grid):
    probe, v, t, labels, arrays = grid
    mutated = {k: a.copy() for k, a in arrays.items()}
    offsets = mutated["offsets"]
    ray = next(i for i in range(len(offsets) - 1) if offsets[i + 1] - offsets[i] > 1)
    mutated["intervals"][offsets[ray] + 1, 0] = mutated["intervals"][offsets[ray], 0]
    assert probe.conforming_partition(v, t, labels, mutated) is False


# ------------------------------------------------------ tetra_interface_mechanism ----
def test_positive_orientation_gate_flips_on_a_flipped_cell():
    probe = _load("tetra_interface_mechanism.py")
    expr = _gate_keyword(ast.parse(_source("tetra_interface_mechanism.py")),
                         "positive_oriented_volumes")
    rows = probe.measure()
    assert _eval(expr, {"a": rows}) is True
    # the declared inverted control is negative and stays excluded
    assert any(x < 0 for r in rows if r["case"] == "inverted" for x in r["cell_six_volumes"])
    flipped = [dict(r) for r in rows]
    target = next(r for r in flipped if r["case"] == "layers")
    target["cell_six_volumes"] = [-x for x in target["cell_six_volumes"]]
    assert _eval(expr, {"a": flipped}) is False
    zero = [dict(r) for r in rows]
    next(r for r in zero if r["case"] == "single")["cell_six_volumes"] = [0.0]
    assert _eval(expr, {"a": zero}) is False


# ----------------------------------------------------------------- nirs_spec_probe ----
def _nirs_rows(box=(True, True), curved=(True, True), narrow=(False, False), full=True):
    return [dict(name=n, frozen_observer_bytes_exact=f, frozen_spec_match=s,
                 full_records_reference_exact=full)
            for n, (f, s) in (("box", box), ("curved", curved), ("narrow", narrow))]


def test_frozen_two_exact_gate_mutations():
    expr = _gate_keyword(ast.parse(_source("nirs_spec_probe.py")), "frozen_two_exact")
    assert _eval(expr, {"rows": _nirs_rows()}) is True
    # a non-matching spec that still reproduces the frozen record (bands ignored)
    assert _eval(expr, {"rows": _nirs_rows(narrow=(True, False))}) is False
    # a matching spec that no longer reproduces the frozen record
    assert _eval(expr, {"rows": _nirs_rows(box=(False, True))}) is False
    # no row matches the frozen parameters at all
    none = _nirs_rows(box=(False, False), curved=(False, False))
    assert _eval(expr, {"rows": none}) is False
    # records differ from the reference propagation
    assert _eval(expr, {"rows": _nirs_rows(full=False)}) is False
    # a missing field is an error, not a silent pass
    broken = _nirs_rows()
    del broken[2]["frozen_observer_bytes_exact"]
    with pytest.raises(KeyError):
        _eval(expr, {"rows": broken})


# -------------------------------------------------------------- photon_modal_suite ----
def _photon_namespace(folder):
    tree = ast.parse(_source("photon_modal_suite.py"))
    families = None
    artifact_names = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id == "FAMILIES":
                families = ast.literal_eval(node.value)
            if node.targets[0].id == "artifact_names":
                artifact_names = node.value
    assert families is not None and artifact_names is not None
    ns = {"FAMILIES": families, "CAPTURES": {"L4": None, "A10G": None, "H100": None},
          "folder": folder}
    ns["artifact_names"] = _eval(artifact_names, ns)
    return tree, ns


def _write_capture_tree(folder, names):
    for backend in ("L4", "A10G", "H100"):
        (folder / backend).mkdir(parents=True, exist_ok=True)
        for name in names:
            (folder / backend / name).write_bytes(b"x")


def test_artifact_integrity_gate_mutations(tmp_path):
    tree, ns = _photon_namespace(tmp_path)
    expr = _gate_keyword(tree, "artifact_integrity")
    assert not isinstance(expr, ast.Constant)
    _write_capture_tree(tmp_path, ns["artifact_names"])
    assert _eval(expr, ns) is True
    missing = sorted(ns["artifact_names"])[0]
    (tmp_path / "A10G" / missing).unlink()
    assert _eval(expr, ns) is False
    (tmp_path / "A10G" / missing).write_bytes(b"x")
    (tmp_path / "H100" / "stray.txt").write_bytes(b"x")
    assert _eval(expr, ns) is False
    (tmp_path / "H100" / "stray.txt").unlink()
    (tmp_path / "L4" / "capture.json").write_bytes(b"")
    assert _eval(expr, ns) is False
    (tmp_path / "L4" / "capture.json").write_bytes(b"x")
    assert _eval(expr, ns) is True
    import shutil
    shutil.rmtree(tmp_path / "A10G")
    assert _eval(expr, ns) is False  # a missing backend directory is a failed gate, not a crash


def test_committed_tetra_report_is_a_run_of_the_current_predicate(grid):
    """reports/tetra_ray_grid.json and its arrays belong together and satisfy the computed
    partition predicate (not the former literal)."""
    import hashlib
    import json

    probe = grid[0]
    report = json.loads((REPO / "reports" / "tetra_ray_grid.json").read_text())
    arrays = np.load(REPO / "reports" / "tetra_ray_grid_arrays.npz")
    assert report["gates"] == {"full_sequences": True, "analytic_exit": True,
                               "full_repeat": True, "partition_accepted": True}
    v, t, labels = probe.fixture()
    for i, row in enumerate(report["rows"]):
        stored = {k: arrays[f"{k}_{i}"] for k in row["hashes"]}
        for k, h in row["hashes"].items():
            assert hashlib.sha256(stored[k].tobytes()).hexdigest() == h, (i, k)
        assert row["cells"] == len(t) and row["refusals"] == [] and row["sequence_errors"] == 0
        assert probe.conforming_partition(v, t, labels, stored) is True
