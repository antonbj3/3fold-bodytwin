"""Chained (cross-consumer) receipts.

A downstream consumer (e.g. blood_oxygen_transport <- cardiac_output) seals its output with
``upstream=<sealed upstream artifact>`` and a verifier proves which upstream payload it used.
"""
import copy
import hashlib
import json

import pytest

import bodytwin.framework.result_envelope_v1 as rev

BASE_KEYS = {"payload_sha256", "producer_sha256", "producer_path", "sealed_utc", "schema"}


def _hash(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


@pytest.fixture
def producers(tmp_path):
    up = tmp_path / "upstream_producer.py"
    up.write_text("# upstream consumer (cardiac_output)\n", encoding="utf-8")
    down = tmp_path / "downstream_producer.py"
    down.write_text("# downstream consumer (blood_oxygen_transport)\n", encoding="utf-8")
    return up, down


def _chain(producers):
    up, down = producers
    a = rev.seal({"cardiac_output_l_min": 5.1}, up)
    b = rev.seal({"do2_ml_min": 1000.0}, down, upstream=a)
    return a, b


def test_two_level_chain_verifies(producers):
    up, down = producers
    a, b = _chain(producers)
    link = b["receipt"]["upstream"]
    assert set(b["receipt"]) == BASE_KEYS | {"upstream"}
    assert link == {
        "name": str(up),
        "payload_sha256": _hash(a["payload"]),
        "producer_sha256": a["receipt"]["producer_sha256"],
        "schema": "result_receipt_v1",
    }
    got = rev.verify_receipt(
        b,
        producer_path=down,
        require_upstream=True,
        expected_upstream_payload_sha256=a["receipt"]["payload_sha256"],
    )
    assert got == b["payload"]


def test_wrong_expected_upstream_is_stale(producers):
    _, b = _chain(producers)
    with pytest.raises(ValueError, match="stale upstream"):
        rev.verify_receipt(b, expected_upstream_payload_sha256="0" * 64)


def test_require_upstream_on_plain_receipt_refused(producers):
    a, _ = _chain(producers)
    with pytest.raises(ValueError, match="missing upstream receipt"):
        rev.verify_receipt(a, require_upstream=True)
    with pytest.raises(ValueError, match="missing upstream receipt"):
        rev.verify_receipt(a, expected_upstream_payload_sha256="0" * 64)
    # no chain arguments -> plain receipt still verifies
    assert rev.verify_receipt(a) == a["payload"]


def test_tampered_upstream_recorded_and_refused(producers):
    up, down = producers
    a, _ = _chain(producers)
    original = a["receipt"]["payload_sha256"]
    tampered = copy.deepcopy(a)
    tampered["payload"]["cardiac_output_l_min"] = 9.9  # receipt left stale
    b = rev.seal({"do2_ml_min": 1800.0}, down, upstream=tampered)
    recorded = b["receipt"]["upstream"]["payload_sha256"]
    assert recorded == _hash(tampered["payload"]) != original  # recomputed, not copied
    with pytest.raises(ValueError, match="stale upstream"):
        rev.verify_receipt(b, require_upstream=True, expected_upstream_payload_sha256=original)
    # the tampered upstream itself fails its own receipt
    with pytest.raises(ValueError, match="payload hash mismatch"):
        rev.verify_receipt(tampered, producer_path=up)


def test_receipt_only_upstream_and_malformed(producers):
    _, down = producers
    a, _ = _chain(producers)
    b = rev.seal({"x": 1}, down, upstream=a["receipt"])
    assert b["receipt"]["upstream"]["payload_sha256"] == a["receipt"]["payload_sha256"]
    bad = dict(a["receipt"])
    del bad["payload_sha256"]
    with pytest.raises(ValueError, match="invalid upstream receipt"):
        rev.seal({"x": 1}, down, upstream=bad)
    with pytest.raises(ValueError, match="invalid upstream receipt"):
        rev.seal({"x": 1}, down, upstream={"payload": {}, "receipt": "nope"})
    with pytest.raises(ValueError, match="invalid upstream receipt"):
        rev.seal({"x": 1}, down, upstream="not-a-dict")
    with pytest.raises(ValueError, match="invalid upstream receipt"):
        rev.seal_chain({"x": 1}, down, None)


def test_seal_chain_read_verified_chain_roundtrip(producers, tmp_path):
    up, down = producers
    a = rev.seal({"ve_l_min": 8.0}, up)
    b = rev.seal_chain({"muscle_w": 120.0}, down, a)
    path = tmp_path / "chained.json"
    path.write_text(json.dumps(b, indent=2) + "\n", encoding="utf-8")
    exp = a["receipt"]["payload_sha256"]
    assert rev.read_verified_chain(
        path, down, expected_upstream_payload_sha256=exp
    ) == b["payload"]
    assert rev.read_verified_chain(
        path, expected_upstream_payload_sha256=exp,
        expected_producer_sha256=b["receipt"]["producer_sha256"],
    ) == b["payload"]
    # a bare upstream link (no independently expected upstream identity) is refused
    with pytest.raises(ValueError, match="expected upstream identity required"):
        rev.read_verified_chain(path, expected_producer_sha256=b["receipt"]["producer_sha256"])
    with pytest.raises(ValueError, match="stale upstream"):
        rev.read_verified_chain(path, expected_upstream_payload_sha256="f" * 64)
    plain = tmp_path / "plain.json"
    rev.write_sealed(plain, {"ve_l_min": 8.0}, up)
    with pytest.raises(ValueError, match="missing upstream receipt"):
        rev.read_verified_chain(plain, up, expected_upstream_payload_sha256=exp)


def test_chaining_on_an_upstream_envelope_with_failed_gates_is_refused(producers, tmp_path):
    up, down = producers
    failed = {"envelope_version": 1, "quantities": {}, "gates": {"g": False},
              "overall_pass": False}
    a = rev.seal(failed, up)
    with pytest.raises(ValueError, match="producer gates failed: upstream envelope"):
        rev.seal_chain({"x": 1}, down, a)
    # an explicit non-scientific chain records the outcome, and the consumer reader refuses it
    b = rev.seal_chain({"x": 1}, down, a, allow_failed_upstream=True)
    assert b["receipt"]["upstream"]["upstream_overall_pass"] is False
    path = tmp_path / "failed_chain.json"
    path.write_text(json.dumps(b), encoding="utf-8")
    with pytest.raises(ValueError, match="producer gates failed: the chained upstream"):
        rev.read_verified_chain(path, down,
                                expected_upstream_payload_sha256=rev.payload_sha256(failed))
    passed = dict(failed, gates={"g": True}, overall_pass=True)
    c = rev.seal_chain({"x": 1}, down, rev.seal(passed, up))
    assert c["receipt"]["upstream"]["upstream_overall_pass"] is True


def test_upstream_less_receipt_keys_unchanged(producers):
    up, _ = producers
    sealed = rev.seal({"a": 1}, up)
    assert set(sealed) == {"payload", "receipt"}
    assert set(sealed["receipt"]) == BASE_KEYS
    assert "upstream" not in sealed["receipt"]
    assert rev.verify_receipt(sealed, producer_path=up) == {"a": 1}


def test_chain_refuses_named_failure_despite_passing_aggregate(producers, tmp_path):
    up, down = producers
    upstream = rev.seal({"envelope_version": "1", "overall_pass": True,
                         "gates": {"mass_balance": False}}, up)
    with pytest.raises(ValueError, match="producer gates failed"):
        rev.seal_chain({"value": 1}, down, upstream)
    # Explicit demo sealing must not hide the failure from a scientific chain reader.
    demo = rev.seal_chain({"value": 1}, down, upstream, allow_failed_upstream=True)
    path = tmp_path / "chain.json"
    path.write_text(json.dumps(demo))
    with pytest.raises(ValueError, match="producer gates failed"):
        rev.read_verified_chain(path, down,
            expected_upstream_payload_sha256=upstream["receipt"]["payload_sha256"])
