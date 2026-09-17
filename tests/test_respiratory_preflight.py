"""Acceptance tests for the shared preflight's third consumer: respiratory.

The respiratory cell reads the SAME ``<BODYTWIN_OUT>/metabolic_cost/metabolic_cost_results.json``
as thermoregulation and cardiac_output (a SUBSET of the nine shared quantities: mass, the four
``headline.*`` w/kg values and the combined tendon+mass correction), so the metabolic-cost
fixtures are reused byte-for-byte (no duplicated inputs). The tests mirror the cardiac tests:

  (a) ``tests/fixtures/i1/synthetic_demo`` (default mode, no registry) -> exit 0, tagged a non-scientific
      synthetic demo (this is what the legacy ``tests/test_cells.py`` environment requires);
  (b) ``tests/fixtures/i1/variants/wrong_unit`` -> exit 2, ``respiratory_refused.json`` with a
      ``unit mismatch`` and NO ``respiratory_results.json``;
  (c) ``tests/fixtures/i1/variants/stale_input`` with strict mode + registry -> exit 2;
  (d) ``tests/fixtures/i1/public_control`` (the real adapted producer result, whose own gates
      failed) -> refused by the preflight, exit 2, no results file;
  (d') ``tests/fixtures/i1/contract_gates_pass`` (contract-test fixture, not evidence) -> the
      preflight itself accepts (``ACCEPT_SCIENTIFIC``). The
      respiratory cell's OWN ten gates are independent of the preflight and one of them
      (``primary_configs_exceed_published_range_expected_direction``) FAILS on the adapted public
      payload, so the process exits 2 while still writing its report. This test therefore asserts
      the preflight state from the written report and documents the exit-2 outcome explicitly (it
      does NOT pretend the cell's own gates pass);
  (e) a ``test_cells``-style legacy env run (cwd = cell dir, no PYTHONPATH, no registry) on the
      synthetic fixture -> exit 0, banner + ``synthetic_demo`` report;
  (f) the generic ``preflight_metabolic_cost`` and the ``preflight_thermoregulation`` wrapper agree
      on the same input (the wrapper is a thin, behaviour-preserving forwarder);
  (g) a REFUSE removes a prior accepted ``respiratory_results.json`` (same hardening as thermo).

Stdlib + pytest only. Subprocesses inherit this interpreter and run with
``BLAS_NUM_THREADS=1 OMP_NUM_THREADS=1``.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
FIXTURES = REPO / "tests" / "fixtures" / "i1"
REGISTRY_PATH = FIXTURES / "registry" / "thermoregulation_registry.json"
RESPIRATORY = SRC / "bodytwin" / "cells" / "respiratory" / "respiratory.py"
SYNTHETIC_INPUT = FIXTURES / "synthetic_demo" / "metabolic_cost" / "metabolic_cost_results.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import bodytwin.framework.consumer_preflight_v1 as cp  # noqa: F401


# --------------------------------------------------------------------------- helpers --------
def _copy_fixture(src_dir, tmp_path):
    if not src_dir.exists():
        pytest.fail(f"fixture directory not built yet: {src_dir}")
    dst = tmp_path / src_dir.name
    shutil.copytree(src_dir, dst)
    return dst


def _run_consumer(out_root, evidence_mode=None, use_registry=True):
    """Run the patched respiratory cell standalone (PYTHONPATH=src) against ``out_root``.

    ``use_registry=False`` is the no-registry path the synthetic-demo test requires; the sealed
    scientific/stale paths need the registry pin and use ``use_registry=True``.
    """
    env = dict(os.environ)
    env["BODYTWIN_OUT"] = str(out_root)
    env["BLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("BODYTWIN_I1_REGISTRY", None)
    if evidence_mode is not None:
        env["BODYTWIN_EVIDENCE_MODE"] = evidence_mode
    if use_registry and REGISTRY_PATH.exists():
        env["BODYTWIN_I1_REGISTRY"] = str(REGISTRY_PATH)
    return subprocess.run(
        [sys.executable, str(RESPIRATORY)],
        env=env,
        cwd=str(REPO),
        capture_output=True,
        text=True,
        timeout=900,
    )


def _report_path(out_root):
    return out_root / "respiratory" / "respiratory_results.json"


def _refused_path(out_root):
    return out_root / "respiratory" / "respiratory_refused.json"


def _read_report(out_root):
    path = _report_path(out_root)
    assert path.exists(), f"consumer did not write report: {path}"
    return json.loads(path.read_text(encoding="utf-8"))


def _registry():
    if not REGISTRY_PATH.exists():
        return None
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


# --------------------------------------------------------------- end-to-end subprocess ------
def test_consumer_synthetic_demo_default_mode(tmp_path):
    """Legacy default ``auto`` on the committed synthetic fixture: exit 0, honest demo tag."""
    out = _copy_fixture(FIXTURES / "synthetic_demo", tmp_path)
    proc = _run_consumer(out, use_registry=False)  # the task's "no registry" path
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "EVIDENCE: SYNTHETIC DEMO" in proc.stdout, proc.stdout

    report = _read_report(out)
    assert report["evidence_class"] == "synthetic_demo"
    assert report["scientific_claim"] is False
    assert report["preflight"]["state"] == "ACCEPT_SYNTHETIC_DEMO"
    assert report["preflight"]["mode"] == "auto"


def test_consumer_wrong_unit_writes_refusal_and_exits_2(tmp_path):
    """A sealed envelope with a wrong unit is refused by the shared preflight before any number."""
    out = _copy_fixture(FIXTURES / "variants" / "wrong_unit", tmp_path)
    proc = _run_consumer(out)
    assert proc.returncode == 2, proc.stdout + proc.stderr

    refused = _refused_path(out)
    assert refused.exists(), f"expected refusal record at {refused}"
    record = json.loads(refused.read_text(encoding="utf-8"))
    assert any("unit mismatch" in f for f in record["failures"]), record["failures"]
    assert not _report_path(out).exists(), "REFUSE must not leave respiratory_results.json"


def test_consumer_stale_input_strict_exits_2(tmp_path):
    """Strict + registry: an envelope whose declared input_sha256 does not match the pin refuses."""
    out = _copy_fixture(FIXTURES / "variants" / "stale_input", tmp_path)
    proc = _run_consumer(out, evidence_mode="strict")
    assert proc.returncode == 2, proc.stdout + proc.stderr

    refused = _refused_path(out)
    assert refused.exists(), f"expected refusal record at {refused}"
    record = json.loads(refused.read_text(encoding="utf-8"))
    assert any("stale input" in f for f in record["failures"]), record["failures"]


def test_consumer_contract_fixture_preflight_state_from_report(tmp_path):
    """The shared preflight ACCEPTS the contract-test fixture (same payload numbers, passing gate
    outcome) as scientific, but the respiratory cell's
    OWN gates do NOT all pass on the adapted public payload.

    Honest, explicit outcome: the cell writes ``respiratory_results.json`` with
    ``preflight.state == "ACCEPT_SCIENTIFIC"`` / ``evidence_class == "calibrated_public_model"``
    AND ``overall_pass is False`` (its ``primary_configs_exceed_published_range_expected_direction``
    gate fails because the adapted payload's primary configs land INSIDE the published walking VO2
    range, i.e. the expected direction is not observed on this payload), so the process exits 2.
    This test asserts the preflight state from the report (as the task directs) and pins the
    independent-gate outcome so exit-2 is not mistaken for a preflight refusal.
    """
    out = _copy_fixture(FIXTURES / "contract_gates_pass", tmp_path)
    proc = _run_consumer(out)  # auto -> strict for a sealed envelope, registry pins hashes
    # A report IS written even when the cell's own gates fail (base behaviour), so we can read
    # the preflight verdict from it.
    report = _read_report(out)
    assert report["preflight"]["state"] == "ACCEPT_SCIENTIFIC"
    assert report["evidence_class"] == "calibrated_public_model"
    assert report["scientific_claim"] is True
    assert "EVIDENCE: SCIENTIFIC" in proc.stdout, proc.stdout
    # Independent respiratory gates fail on this payload -- documented, not hidden.
    assert report["overall_pass"] is False
    assert report["gates"]["primary_configs_exceed_published_range_expected_direction"] is False
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert not _refused_path(out).exists(), "the preflight accepted; there must be no refusal record"


def test_consumer_public_control_is_refused(tmp_path):
    """The real adapted producer result (its own gates failed) is refused before any number is
    used: exit 2, a refusal record, no results file, no scientific banner."""
    out = _copy_fixture(FIXTURES / "public_control", tmp_path)
    proc = _run_consumer(out)
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "REFUSED: producer gates failed" in proc.stdout
    assert "EVIDENCE: SCIENTIFIC" not in proc.stdout
    assert not _report_path(out).exists()
    record = json.loads(_refused_path(out).read_text(encoding="utf-8"))
    assert "producer gates failed: overall_pass is not true" in record["failures"]


def test_refuse_removes_stale_accepted_output(tmp_path):
    """A later REFUSE must not leave a prior accepted respiratory_results.json for readers."""
    out = _copy_fixture(FIXTURES / "synthetic_demo", tmp_path)
    first = _run_consumer(out, use_registry=False)
    assert first.returncode == 0, first.stdout + first.stderr
    assert _report_path(out).exists(), "accepted run did not write a report"

    bad = FIXTURES / "variants" / "wrong_unit" / "metabolic_cost" / "metabolic_cost_results.json"
    shutil.copyfile(bad, out / "metabolic_cost" / "metabolic_cost_results.json")
    second = _run_consumer(out)
    assert second.returncode == 2, second.stdout + second.stderr
    assert not _report_path(out).exists(), "REFUSE left a stale accepted respiratory_results.json"
    assert _refused_path(out).exists()


def test_legacy_test_cells_env_runs_synthetic_as_demo(tmp_path):
    """Old passing examples must keep working. ``tests/test_cells.py`` runs every cell as
    ``subprocess.run([python, script.name], cwd=script.parent, env=BODYTWIN_OUT=...)`` with NO
    PYTHONPATH and NO registry. In that exact environment the committed synthetic fixture must
    still exit 0 and be tagged a non-scientific synthetic demo."""
    out = _copy_fixture(FIXTURES / "synthetic_demo", tmp_path)
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("BODYTWIN_I1_REGISTRY", None)
    env.pop("BODYTWIN_EVIDENCE_MODE", None)
    env["BODYTWIN_OUT"] = str(out)
    env["BLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.run(
        [sys.executable, RESPIRATORY.name],
        cwd=str(RESPIRATORY.parent),
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert proc.stdout.strip(), "a cell must print its own numbers"
    assert "EVIDENCE: SYNTHETIC DEMO" in proc.stdout, proc.stdout
    report = _read_report(out)
    assert report["evidence_class"] == "synthetic_demo"
    assert report["scientific_claim"] is False
    assert report["preflight"]["state"] == "ACCEPT_SYNTHETIC_DEMO"


# ----------------------------------------------------------- generic API ---------------------
def test_generic_preflight_and_wrapper_agree_and_spec_alias_is_shared(tmp_path):
    """The extracted generic preflight and the thermoregulation wrapper are behaviour-identical,
    and the exported metabolic-cost spec alias names the same object (the respiratory cell routes
    through the generic entry point)."""
    if not SYNTHETIC_INPUT.exists():
        pytest.fail(f"fixture not built yet: {SYNTHETIC_INPUT}")
    assert cp.METABOLIC_COST_QUANTITY_SPECS is cp.THERMO_QUANTITY_SPECS

    raw = json.loads(SYNTHETIC_INPUT.read_text(encoding="utf-8"))
    generic_payload, generic_verdict = cp.preflight_metabolic_cost(raw, mode="auto")
    wrapper_payload, wrapper_verdict = cp.preflight_thermoregulation(raw, mode="auto")

    assert generic_verdict.state == wrapper_verdict.state == "ACCEPT_SYNTHETIC_DEMO"
    assert generic_verdict.failures == wrapper_verdict.failures
    assert generic_verdict.flags == wrapper_verdict.flags
    assert generic_verdict.evidence_class == wrapper_verdict.evidence_class
    assert generic_payload == wrapper_payload == raw


def test_generic_preflight_refuses_wrong_unit_fixture():
    """The generic entry point applies the same shared contract to the same bad fixture the
    respiratory cell consumes."""
    path = FIXTURES / "variants" / "wrong_unit" / "metabolic_cost" / "metabolic_cost_results.json"
    if not path.exists():
        pytest.fail(f"fixture not built yet: {path}")
    raw = json.loads(path.read_text(encoding="utf-8"))
    registry = _registry()
    if registry is None:
        pytest.fail("registry unavailable")
    rentry = cp.registry_entry(registry, "metabolic_cost_results.json")
    current = hashlib.sha256(json.dumps(raw, indent=2).encode("utf-8")).hexdigest()
    payload, verdict = cp.preflight_metabolic_cost(
        raw, mode="auto", registry_entry=rentry, current_input_sha256=current,
    )
    assert payload is None
    assert verdict.state == "REFUSE"
    assert any("unit mismatch" in f for f in verdict.failures), verdict.failures
