"""Acceptance tests for the I1 consumer preflight + patched thermoregulation consumer.

Two layers:

1. **In-process** drive of every ``tests/fixtures/i1/manifest.json`` entry through
   ``consumer_preflight_v1.preflight_thermoregulation``, asserting the declared
   ``expected_state`` and ``expected_failure_substring``.

2. **Subprocess** end-to-end drive of the patched
   ``src/bodytwin/cells/organ_systems/thermoregulation.py`` with ``BODYTWIN_OUT`` pointed at a
   throw-away copy of each fixture, asserting exit codes and the reports it writes:
     (a) ``tests/fixtures/i1/contract_gates_pass``       -> exit 0, scientific evidence, ACCEPT_SCIENTIFIC
     (a') ``tests/fixtures/i1/public_control``           -> exit 2, producer gates failed
     (b) ``tests/fixtures/i1/synthetic_demo`` (default)  -> exit 0, evidence_class synthetic_demo
     (c) ``tests/fixtures/i1/variants/wrong_unit``       -> exit 2, thermoregulation_refused.json written
     (d) ``tests/fixtures/i1/variants/stale_input`` strict -> exit 2
     (e) ``tests/fixtures/i1/variants/wrong_regime``     -> exit 2, refusal contains ``regime mismatch``

Missing framework modules, fixtures or the test registry are errors, not skips.

Stdlib + pytest only.  The subprocess inherits this interpreter (the venv python) and runs with
``BLAS_NUM_THREADS=1`` / ``OMP_NUM_THREADS=1`` and ``PYTHONPATH=<repo>/src``.
"""
import copy
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
CONSUMER = SRC / "bodytwin" / "cells" / "organ_systems" / "thermoregulation.py"
MANIFEST_PATH = FIXTURES / "manifest.json"

# The registered adapted producer the sealed fixtures were built against. Re-sealing a mutated
# envelope with any other file must be refused at receipt level (hardening H2).
# Not part of this repository: located by BODYTWIN_REF_ADAPTED_PRODUCER (tests/research_data.py).
SYNTHETIC_INPUT = FIXTURES / "synthetic_demo" / "metabolic_cost" / "metabolic_cost_results.json"
# Contract-test fixture (gate outcome replaced so the accepting path can be exercised; not
# evidence) and the real adapted producer result, whose own gates failed.
VALID_INPUT = FIXTURES / "contract_gates_pass" / "metabolic_cost" / "metabolic_cost_results.json"
CONTROL_INPUT = FIXTURES / "public_control" / "metabolic_cost" / "metabolic_cost_results.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import bodytwin.framework.consumer_preflight_v1 as cp  # noqa: F401

import bodytwin.framework.result_envelope_v1 as rev  # noqa: F401

_MANIFEST = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
_MANIFEST_ITEMS = sorted(_MANIFEST.items())

SCIENTIFIC_TERMS = {"measured", "calibrated", "calibrated_public_model"}


# --------------------------------------------------------------------------- helpers --------
def _fixture_path(entry, name):
    raw = entry.get("path") or name
    path = Path(raw)
    if not path.is_absolute():
        path = REPO / path
    if path.is_dir():
        path = path / "metabolic_cost" / "metabolic_cost_results.json"
    return path


def _registry():
    if not REGISTRY_PATH.exists():
        return None
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


def _registry_entry_or_fail():
    registry = _registry()
    if registry is None:
        pytest.fail("registry unavailable")
    rentry = cp.registry_entry(registry, "metabolic_cost_results.json")
    if rentry is None:
        pytest.fail("registry entry unavailable")
    return rentry


def _seal_with_registered_producer(payload):
    from research_data import optional_path

    return rev.seal(payload, optional_path("BODYTWIN_REF_ADAPTED_PRODUCER"))


def _preflight_strict(sealed, registry_entry):
    current = hashlib.sha256(json.dumps(sealed, indent=2).encode("utf-8")).hexdigest()
    return cp.preflight_thermoregulation(
        sealed,
        mode="strict",
        registry_entry=registry_entry,
        current_input_sha256=current,
        producer_path=None,
    )


def _copy_fixture(src_dir, tmp_path):
    if not src_dir.exists():
        pytest.fail(f"fixture directory not built yet: {src_dir}")
    dst = tmp_path / src_dir.name
    shutil.copytree(src_dir, dst)
    return dst


def _run_consumer(out_root, evidence_mode=None):
    env = dict(os.environ)
    env["BODYTWIN_OUT"] = str(out_root)
    env["BLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    if evidence_mode is not None:
        env["BODYTWIN_EVIDENCE_MODE"] = evidence_mode
    if REGISTRY_PATH.exists():
        env["BODYTWIN_I1_REGISTRY"] = str(REGISTRY_PATH)
    return subprocess.run(
        [sys.executable, str(CONSUMER)],
        env=env,
        cwd=str(REPO),
        capture_output=True,
        text=True,
        timeout=900,
    )


def _read_report(out_root):
    path = out_root / "thermoregulation" / "thermoregulation_results.json"
    assert path.exists(), f"consumer did not write report: {path}"
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------- in-process manifest ---
@pytest.mark.parametrize(
    "name,entry", _MANIFEST_ITEMS, ids=[n for n, _ in _MANIFEST_ITEMS]
)
def test_manifest_entry_preflight(name, entry):
    path = _fixture_path(entry, name)
    if not path.exists():
        pytest.fail(f"fixture not built yet: {path}")

    raw = json.loads(path.read_text(encoding="utf-8"))
    kind = cp.classify_input(raw)
    mode = entry.get("mode") or "auto"

    registry = _registry()
    rentry = None
    if registry is not None:
        rentry = cp.registry_entry(registry, path.name)
    if kind == cp.SEALED and rentry is None:
        pytest.fail("registry entry unavailable; cannot pin producer/input hashes")

    current_input_sha256 = (rentry or {}).get("input_sha256")
    payload, verdict = cp.preflight_thermoregulation(
        raw,
        mode=mode,
        registry_entry=rentry,
        current_input_sha256=current_input_sha256,
        producer_path=None,
    )

    assert verdict.state == entry["expected_state"], (name, verdict)
    expected_sub = entry.get("expected_failure_substring")
    if expected_sub:
        assert any(expected_sub in f for f in verdict.failures), (
            name, expected_sub, verdict.failures,
        )
    if verdict.state in ("ACCEPT_SCIENTIFIC", "ACCEPT_SYNTHETIC_DEMO"):
        assert payload is not None, (name, "accepted but no payload")


def test_manifest_covers_the_required_defect_classes():
    expected = {
        "unit mismatch", "region mismatch", "unknown region_id", "frame mismatch",
        "time mismatch", "regime mismatch", "missing regime", "stale input",
        "producer changed", "disallowed provenance", "stale result", "producer gates failed",
    }
    covered = {
        entry.get("expected_failure_substring")
        for entry in _MANIFEST.values()
        if entry.get("expected_failure_substring")
    }
    missing = expected - covered
    assert not missing, f"manifest does not cover defect substrings: {sorted(missing)}"


def test_public_control_with_failed_producer_gates_is_refused():
    """The real adapted public gait2392 result reports ``overall_pass=false``
    (``speed_sanity_ok=false``, ``expected_direction_pass=false``). Scientific acceptance refuses
    it by default (``producer gates failed``); strict mode cannot waive the gates; the explicit
    demo mode accepts it only as a non-scientific run and still shows the failed gates."""
    raw = json.loads(CONTROL_INPUT.read_text(encoding="utf-8"))
    assert raw["payload"]["overall_pass"] is False
    rentry = _registry_entry_or_fail()
    current = hashlib.sha256(json.dumps(raw, indent=2).encode("utf-8")).hexdigest()

    for mode in ("auto", "strict"):
        payload, verdict = cp.preflight_thermoregulation(
            raw, mode=mode, registry_entry=rentry, current_input_sha256=current,
            producer_path=None,
        )
        assert payload is None
        assert verdict.state == "REFUSE"
        assert "producer gates failed: overall_pass is not true" in verdict.failures

    with pytest.raises(ValueError, match="strict mode cannot waive producer gates"):
        cp.preflight_thermoregulation(
            raw, mode="strict", registry_entry=rentry, current_input_sha256=current,
            producer_path=None, require_producer_gates_pass=False,
        )

    payload, demo = cp.preflight_thermoregulation(
        raw, mode="synthetic_demo", registry_entry=rentry, current_input_sha256=current,
        producer_path=None,
    )
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert "producer gates failed: overall_pass is not true" in demo.flags


def test_contract_fixture_is_accepted_only_with_the_registry_pin():
    raw = json.loads(VALID_INPUT.read_text(encoding="utf-8"))
    rentry = _registry_entry_or_fail()
    assert rentry["payload_sha256"] == rev.payload_sha256(raw["payload"])
    payload, verdict = cp.preflight_thermoregulation(
        raw, mode="strict", registry_entry=rentry, producer_path=None)
    assert verdict.state == "ACCEPT_SCIENTIFIC", verdict
    unpinned = {k: v for k, v in rentry.items() if k != "payload_sha256"}
    payload, verdict = cp.preflight_thermoregulation(
        raw, mode="strict", registry_entry=unpinned, producer_path=None)
    assert verdict.state == "REFUSE"
    assert ("no trusted artifact pin: expected_payload_sha256 is required for scientific "
            "acceptance") in verdict.failures


# --------------------------------------------------------------- end-to-end subprocess ------
def test_consumer_public_control_is_refused(tmp_path):
    out = _copy_fixture(FIXTURES / "public_control", tmp_path)
    proc = _run_consumer(out)  # default mode = auto -> strict for a sealed envelope
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "REFUSED: producer gates failed" in proc.stdout
    assert not (out / "thermoregulation" / "thermoregulation_results.json").exists()
    refused = json.loads(
        (out / "thermoregulation" / "thermoregulation_refused.json").read_text(encoding="utf-8"))
    assert "producer gates failed: overall_pass is not true" in refused["failures"]


def test_consumer_contract_fixture_is_accepted(tmp_path):
    out = _copy_fixture(FIXTURES / "contract_gates_pass", tmp_path)
    proc = _run_consumer(out)  # default mode = auto -> strict for a sealed envelope
    assert proc.returncode == 0, proc.stdout + proc.stderr

    report = _read_report(out)
    assert report["evidence_class"] in SCIENTIFIC_TERMS, report["evidence_class"]
    assert report["preflight"]["state"] == "ACCEPT_SCIENTIFIC"
    assert report.get("scientific_claim") is True
    assert not any("producer gates failed" in f for f in report["preflight"]["flags"])


def test_consumer_synthetic_demo_default_mode(tmp_path):
    out = _copy_fixture(FIXTURES / "synthetic_demo", tmp_path)
    proc = _run_consumer(out)  # default/auto: raw synthetic -> synthetic_demo
    assert proc.returncode == 0, proc.stdout + proc.stderr

    report = _read_report(out)
    assert report["evidence_class"] == "synthetic_demo"
    assert report["scientific_claim"] is False
    assert report["preflight"]["state"] == "ACCEPT_SYNTHETIC_DEMO"


def test_consumer_wrong_unit_writes_refusal_and_exits_2(tmp_path):
    out = _copy_fixture(FIXTURES / "variants" / "wrong_unit", tmp_path)
    proc = _run_consumer(out)
    assert proc.returncode == 2, proc.stdout + proc.stderr

    refused = out / "thermoregulation" / "thermoregulation_refused.json"
    assert refused.exists(), f"expected refusal record at {refused}"
    record = json.loads(refused.read_text(encoding="utf-8"))
    assert any("unit mismatch" in f for f in record["failures"]), record["failures"]


def test_consumer_stale_input_strict_exits_2(tmp_path):
    out = _copy_fixture(FIXTURES / "variants" / "stale_input", tmp_path)
    proc = _run_consumer(out, evidence_mode="strict")
    assert proc.returncode == 2, proc.stdout + proc.stderr


def test_consumer_wrong_regime_writes_refusal_and_exits_2(tmp_path):
    """End-to-end regime gate: an out-of-domain regime (anaerobic sprint) is not scientific."""
    out = _copy_fixture(FIXTURES / "variants" / "wrong_regime", tmp_path)
    proc = _run_consumer(out)  # default mode = auto -> strict for a sealed envelope
    assert proc.returncode == 2, proc.stdout + proc.stderr

    refused = out / "thermoregulation" / "thermoregulation_refused.json"
    assert refused.exists(), f"expected refusal record at {refused}"
    record = json.loads(refused.read_text(encoding="utf-8"))
    assert any("regime mismatch" in f for f in record["failures"]), record["failures"]


# ---------------------------------------------- hardening regressions (H1-H5) ----------------
def test_sealed_bare_raw_payload_is_refused():
    """H1/H5: a sealed payload that is a raw result (no ``envelope_version``) must be refused;
    it must not be re-wrapped with registry metadata and promoted to scientific evidence."""
    if not SYNTHETIC_INPUT.exists():
        pytest.fail(f"fixture not built yet: {SYNTHETIC_INPUT}")
    raw = json.loads(SYNTHETIC_INPUT.read_text(encoding="utf-8"))
    sealed = _seal_with_registered_producer(raw)
    payload, verdict = _preflight_strict(sealed, _registry_entry_or_fail())
    assert payload is None
    assert verdict.state == "REFUSE"
    assert any(
        "sealed payload is not a ResultEnvelope" in f for f in verdict.failures
    ), verdict.failures


def test_sibling_envelope_does_not_masquerade_for_the_payload():
    """H5: a valid sibling top-level ``envelope`` must not be validated while a different bare
    payload is consumed. The payload here is the attacker's bare result; the sibling is a clean,
    valid envelope. The whole input must be refused."""
    if not VALID_INPUT.exists():
        pytest.fail(f"fixture not built yet: {VALID_INPUT}")
    valid_payload = json.loads(VALID_INPUT.read_text(encoding="utf-8"))["payload"]
    attacker_b = {
        "synthetic": True,
        "note": "attacker-controlled payload that is NOT a ResultEnvelope",
        "muscle_mass": {"total_body_mass_kg": 999.0, "summed_muscle_mass_kg": 888.0},
        "headline": {
            "umberger2010_gross_w_per_kg": 7.2,
            "umberger2010_net_w_per_kg": 1.162,
            "bhargava2004_gross_w_per_kg": 7.0,
            "bhargava2004_net_w_per_kg": 1.2,
        },
        "sensitivity": {
            "combined_tendon_and_mass_correction_w_per_kg": 5.0,
            "muscle_mass_correction_scale_applied": 1.0,
        },
        "distance_speed": {"speed_mps": 1.4},
    }
    sealed = _seal_with_registered_producer(attacker_b)
    sealed["envelope"] = copy.deepcopy(valid_payload)  # the sibling A
    payload, verdict = _preflight_strict(sealed, _registry_entry_or_fail())
    assert payload is None
    assert verdict.state == "REFUSE"
    assert any(
        "sealed payload is not a ResultEnvelope" in f for f in verdict.failures
    ), verdict.failures


def test_foreign_producer_receipt_is_refused(tmp_path):
    """H2: a valid envelope sealed against a foreign producer file must be refused because the
    receipt producer hash is now bound to the registry pin."""
    if not VALID_INPUT.exists():
        pytest.fail(f"fixture not built yet: {VALID_INPUT}")
    valid_payload = json.loads(VALID_INPUT.read_text(encoding="utf-8"))["payload"]
    foreign = tmp_path / "foreign_producer.py"
    foreign.write_text("# a foreign producer, never the registered one\n", encoding="utf-8")
    sealed = rev.seal(copy.deepcopy(valid_payload), foreign)
    payload, verdict = _preflight_strict(sealed, _registry_entry_or_fail())
    assert payload is None
    assert verdict.state == "REFUSE"
    assert any("producer changed" in f for f in verdict.failures), verdict.failures


def test_physically_impossible_muscle_mass_is_refused():
    """H4: the physical-impossibility check is wired through the consumer preflight. A local
    working-muscle mass above the total body mass must be refused."""
    if not VALID_INPUT.exists():
        pytest.fail(f"fixture not built yet: {VALID_INPUT}")
    envelope = copy.deepcopy(json.loads(VALID_INPUT.read_text(encoding="utf-8"))["payload"])
    envelope["quantities"]["muscle_mass.summed_muscle_mass_kg"]["value"] = 100.0
    envelope["quantities"]["muscle_mass.total_body_mass_kg"]["value"] = 64.0
    # The real fixture's correction scale (0.574) would make the local mass 57.4 kg < 64 kg;
    # pin it to 1.0 so the mutated local mass is genuinely impossible (local > total).
    envelope["quantities"]["sensitivity.muscle_mass_correction_scale_applied"]["value"] = 1.0
    sealed = _seal_with_registered_producer(envelope)
    payload, verdict = _preflight_strict(sealed, _registry_entry_or_fail())
    assert payload is None
    assert verdict.state == "REFUSE"
    assert any("physically impossible mass" in f for f in verdict.failures), verdict.failures


def test_refuse_removes_stale_accepted_output(tmp_path):
    """H3: a REFUSE must not leave a prior ACCEPT_SCIENTIFIC result for a downstream reader."""
    out = _copy_fixture(FIXTURES / "contract_gates_pass", tmp_path)
    first = _run_consumer(out)
    assert first.returncode == 0, first.stdout + first.stderr
    results = out / "thermoregulation" / "thermoregulation_results.json"
    assert results.exists(), "accepted run did not write a results file"

    bad = FIXTURES / "variants" / "wrong_unit" / "metabolic_cost" / "metabolic_cost_results.json"
    shutil.copyfile(bad, out / "metabolic_cost" / "metabolic_cost_results.json")
    second = _run_consumer(out)
    assert second.returncode == 2, second.stdout + second.stderr
    assert not results.exists(), "REFUSE left a stale accepted thermoregulation_results.json"
    assert (out / "thermoregulation" / "thermoregulation_refused.json").exists()


def test_invalid_evidence_mode_refuses_without_traceback(tmp_path):
    """An unknown BODYTWIN_EVIDENCE_MODE is a refusal (exit 2), not an uncaught traceback."""
    out = _copy_fixture(FIXTURES / "synthetic_demo", tmp_path)
    proc = _run_consumer(out, evidence_mode="banana")
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "Traceback" not in (proc.stderr or ""), proc.stderr
    refused = out / "thermoregulation" / "thermoregulation_refused.json"
    assert refused.exists(), f"expected refusal record at {refused}"
    record = json.loads(refused.read_text(encoding="utf-8"))
    assert any("invalid evidence mode" in f for f in record["failures"]), record["failures"]


def test_legacy_test_cells_env_runs_synthetic_as_demo(tmp_path):
    """Old passing examples must keep working. tests/test_cells.py runs every cell as
    ``subprocess.run([python, script.name], cwd=script.parent, env=BODYTWIN_OUT=...)`` with NO
    PYTHONPATH and NO registry. In that exact environment the committed synthetic fixture must
    still exit 0, print numbers, and be tagged as a non-scientific synthetic demo."""
    out = _copy_fixture(FIXTURES / "synthetic_demo", tmp_path)
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("BODYTWIN_I1_REGISTRY", None)
    env.pop("BODYTWIN_EVIDENCE_MODE", None)
    env["BODYTWIN_OUT"] = str(out)
    env["BLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.run(
        [sys.executable, CONSUMER.name],
        cwd=str(CONSUMER.parent),
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
