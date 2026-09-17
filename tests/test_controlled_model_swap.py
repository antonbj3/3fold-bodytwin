"""Controlled model swap: a producer change is accepted ONLY via an explicit pin.

This is a **contract-mechanics** test. It copies or claims **no** physiology: it takes the
*envelope payload* of the existing healthy control and re-seals it against two deterministic
producer stubs (``tests/fixtures/i1/model_swap/producer_v1.py`` / ``producer_v2.py``), then shows that:

  * a payload sealed by v2 (and declaring v2's ``source_sha256``) is REFUSED as ``producer
    changed`` while the registry pins v1;
  * the *same* sealed payload is ACCEPTED as ``ACCEPT_SCIENTIFIC`` once the registry pin is
    explicitly updated to v2 -- the pin is the change-control point;
  * a cross-swap (receipt ``producer_sha256`` = v2 but the envelope forges
    ``producer.source_sha256`` back to v1) is refused, because the receipt pin and the envelope
    declaration disagree.

The fixtures are built deterministically by ``tests/fixtures/i1/build_model_swap_fixtures.py`` (stdlib only,
no network, no repo dependency). The two producer stub hashes are recorded below and printed by
``test_producer_stub_hashes_are_pinned_and_distinct``.
"""
import hashlib
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
SWAP_DIR = REPO / "tests" / "fixtures" / "i1" / "model_swap"
MANIFEST_PATH = SWAP_DIR / "manifest.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import bodytwin.framework.consumer_preflight_v1 as cp

import bodytwin.framework.result_envelope_v1 as rev

_MANIFEST = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

# Recorded deterministic stub hashes (the pin identities). They are reproducible from the fixed
# literals in tests/fixtures/i1/build_model_swap_fixtures.py; a change here means the fixture must be re-pinned.
PRODUCER_V1_SHA256 = "23332ca5b3a85ced2a23cfe6f007a1a1bba9cc4faf4a7165ad3f64b566b23eae"
PRODUCER_V2_SHA256 = "ac4477151fed282aa4912e53d7d9bb2d957ef99dc01ce402df3fbde00658b23a"


# --------------------------------------------------------------------------- helpers -------
def _file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load_fixture(name):
    rel = _MANIFEST["fixtures"][name]["path"]
    path = REPO / rel
    assert path.exists(), f"model-swap fixture missing: {path}"
    return json.loads(path.read_text(encoding="utf-8"))


def _preflight(name, pin):
    """Drive the shared thermo preflight with the registry pin ``pin``."""
    raw = _load_fixture(name)
    return cp.preflight_thermoregulation(
        raw,
        mode="strict",
        registry_entry=_MANIFEST["registry_pins"][pin],
        current_input_sha256=_MANIFEST["input_sha256"],
        producer_path=None,
    )


def _has(failures, needle):
    return any(needle in f for f in failures), tuple(failures)


# ------------------------------------------------------------------- fixture integrity --------
def test_producer_stub_hashes_are_pinned_and_distinct():
    v1 = _file_sha256(SWAP_DIR / "producer_v1.py")
    v2 = _file_sha256(SWAP_DIR / "producer_v2.py")
    assert v1 != v2
    assert v1 == PRODUCER_V1_SHA256, f"producer_v1.py hash moved: {v1}"
    assert v2 == PRODUCER_V2_SHA256, f"producer_v2.py hash moved: {v2}"
    assert _MANIFEST["producers"]["producer_v1.py"]["sha256"] == v1
    assert _MANIFEST["producers"]["producer_v2.py"]["sha256"] == v2
    # Record the two stub hashes in the test output (visible with ``pytest -s``) and in the docs.
    print(f"producer_v1.py sha256={v1}")
    print(f"producer_v2.py sha256={v2}")


def test_fixture_receipts_and_declarations_match_their_registry_pins():
    for name in ("sealed_by_v1", "sealed_by_v2", "cross_swap"):
        sealed = _load_fixture(name)
        entry = _MANIFEST["fixtures"][name]
        assert sealed["receipt"]["schema"] == rev.RECEIPT_SCHEMA
        assert sealed["receipt"]["producer_sha256"] == entry["receipt_producer_sha256"]
        assert sealed["payload"]["producer"]["source_sha256"] == entry["declared_source_sha256"]
        # The receipt binds the payload it actually sealed.
        assert rev.verify_receipt(
            sealed, expected_producer_sha256=entry["receipt_producer_sha256"]
        ) == sealed["payload"]

    # The two stubs are genuinely different identities.
    assert (
        _MANIFEST["fixtures"]["sealed_by_v1"]["declared_source_sha256"]
        != _MANIFEST["fixtures"]["sealed_by_v2"]["declared_source_sha256"]
    )


# --------------------------------------------------------------- controlled model swap --------
def test_pin_v1_accepts_the_v1_sealed_payload():
    payload, verdict = _preflight("sealed_by_v1", "pin_v1")
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.accepted is True
    assert verdict.failures == ()
    assert payload is not None


def test_v2_sealed_payload_is_refused_while_the_registry_still_pins_v1():
    """The swap is not silently accepted: the v1 pin refuses the v2 producer at the receipt."""
    payload, verdict = _preflight("sealed_by_v2", "pin_v1")
    assert verdict.state == "REFUSE"
    assert verdict.accepted is False
    assert payload is None
    ok, failures = _has(verdict.failures, "producer changed")
    assert ok, failures


def test_v2_sealed_payload_is_accepted_once_the_pin_is_explicitly_updated_to_v2():
    """The SAME sealed payload becomes scientific evidence after the pin update -- and only then."""
    payload, verdict = _preflight("sealed_by_v2", "pin_v2")
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.accepted is True
    assert verdict.failures == ()
    assert payload is not None
    assert verdict.evidence_class == "calibrated_public_model"


def test_pin_v2_refuses_the_v1_sealed_payload():
    """The control is symmetric: a rollback is also an explicit pin change, never silent."""
    payload, verdict = _preflight("sealed_by_v1", "pin_v2")
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "producer changed")
    assert ok, failures


def test_cross_swap_forged_envelope_source_is_refused():
    """Receipt pins v2 (the file that sealed it) but the envelope forges ``source_sha256`` to v1.

    The receipt pin is satisfied, so the failure must come from the envelope's own
    ``producer.source_sha256`` check -- the receipt and the declaration disagree.
    """
    payload, verdict = _preflight("cross_swap", "pin_v2")
    assert verdict.state == "REFUSE"
    assert payload is None
    ok, failures = _has(verdict.failures, "producer changed: source_sha256")
    assert ok, failures

    # And a v1 pin refuses it one step earlier, at the receipt (receipt producer_sha256 = v2).
    _payload_v1_pin, verdict_v1 = _preflight("cross_swap", "pin_v1")
    assert verdict_v1.state == "REFUSE"
    ok, failures = _has(verdict_v1.failures, "producer changed")
    assert ok, failures
