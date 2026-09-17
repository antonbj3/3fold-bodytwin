#!/usr/bin/env python
"""Build the I1 controlled-model-swap fixtures (deterministic, stdlib only).

This demonstrates that a producer/model change is accepted ONLY through an explicit registry
pin update, never silently. This tool is pure contract mechanics -- it makes **no** physiological
claim. It takes the *envelope payload* of the contract-test fixture
(``tests/fixtures/i1/contract_gates_pass/metabolic_cost/metabolic_cost_results.json``, a contract-test fixture whose gate outcome passes) and re-seals it against
two distinct deterministic producer stubs so a caller can exercise the "producer changed" pin:

Run from the repository root (after ``build_fixtures.py``)::

    PYTHONPATH=src python tests/fixtures/i1/build_model_swap_fixtures.py

Output (below ``tests/fixtures/i1/model_swap/``)::

    producer_v1.py        deterministic stub producer, version 1
    producer_v2.py        deterministic stub producer, version 2 (different bytes/hash)
    sealed_by_v1.json     envelope declares v1, receipt producer_sha256 == v1
    sealed_by_v2.json     envelope declares v2, receipt producer_sha256 == v2
    cross_swap.json       receipt producer_sha256 == v2, envelope source_sha256 forged back to v1
    manifest.json         stub hashes + per-fixture declared/receipt hashes + registry pins

Determinism: the producer scripts are fixed literals and each receipt's ``sealed_utc`` is pinned;
the receipt hash covers only the payload, so overriding the timestamp does not weaken the
payload/producer binding. The two stub hashes are recorded in the manifest and in
``tests/test_controlled_model_swap.py``.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]  # repository root
sys.path.insert(0, str(ROOT / "src"))

from bodytwin.framework.result_envelope_v1 import file_sha256, seal  # noqa: E402

FIXTURES = HERE
SWAP_DIR = FIXTURES / "model_swap"
VALID_INPUT = FIXTURES / "contract_gates_pass" / "metabolic_cost" / "metabolic_cost_results.json"

FIXED_SEALED_UTC = "2026-09-16T00:00:00+00:00"
PRODUCER_NAME = "model_swap_stub"

# Deterministic stub producers. They are intentionally trivial: their ONLY role is to be two
# distinct, reproducible byte-streams whose sha256 the registry can pin. Deliberately different
# content (a different version constant and body) so the hashes differ.
PRODUCER_V1 = '''#!/usr/bin/env python
"""Deterministic I1 model-swap producer stub, version 1 (contract mechanics only)."""
PRODUCER_VERSION = 1
MODEL_ID = "bodytwin-model-swap-stub-v1"


def produce():
    """Return the stub's (non-physiological) result marker."""
    return {"model_id": MODEL_ID, "version": PRODUCER_VERSION}


if __name__ == "__main__":
    print(produce())
'''

PRODUCER_V2 = '''#!/usr/bin/env python
"""Deterministic I1 model-swap producer stub, version 2 (contract mechanics only)."""
PRODUCER_VERSION = 2
MODEL_ID = "bodytwin-model-swap-stub-v2"
CHANGED = "v2 changes the model identity, so its sha256 differs from v1"


def produce():
    """Return the stub's (non-physiological) result marker."""
    return {"model_id": MODEL_ID, "version": PRODUCER_VERSION, "note": CHANGED}


if __name__ == "__main__":
    print(produce())
'''


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _seal_fixed(payload, producer_path, when=FIXED_SEALED_UTC):
    """Seal and pin ``sealed_utc`` so the fixture bytes are reproducible."""
    sealed = seal(payload, producer_path)
    sealed["receipt"]["sealed_utc"] = when
    return sealed


def _write_json(path: Path, obj) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    return file_sha256(path)


def _payload_sha(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _declaring_payload(base_payload, *, source_sha256, producer_path):
    """Deep-copy the contract-test envelope and declare ``producer_path`` as its producer.

    The deep copy carries the base payload's ``gates`` / ``overall_pass`` along verbatim; only the
    producer identity is rewritten here.
    """
    payload = copy.deepcopy(base_payload)
    payload["producer"]["name"] = PRODUCER_NAME
    payload["producer"]["source_sha256"] = source_sha256
    payload["producer"]["command"] = str(producer_path)
    return payload


def main():
    if not VALID_INPUT.exists():
        raise SystemExit(f"healthy control fixture missing: {VALID_INPUT}")
    base_payload = json.loads(VALID_INPUT.read_text(encoding="utf-8"))["payload"]
    input_sha256 = base_payload["producer"]["input_sha256"]
    # Record the producer gate outcome the re-sealed stubs preserve (the contract-test base's
    # passing gate). A missing field is normalised to the make_envelope default ({} / True)
    # exactly as the envelope writer does.
    gate_outcome = {
        "gates": base_payload.get("gates") or {},
        "overall_pass": base_payload.get("overall_pass", True),
    }

    SWAP_DIR.mkdir(parents=True, exist_ok=True)
    v1_path = SWAP_DIR / "producer_v1.py"
    v2_path = SWAP_DIR / "producer_v2.py"
    v1_path.write_text(PRODUCER_V1, encoding="utf-8", newline="\n")
    v2_path.write_text(PRODUCER_V2, encoding="utf-8", newline="\n")

    v1_sha = file_sha256(v1_path)
    v2_sha = file_sha256(v2_path)
    # Record producer paths relative to the repository root (never an absolute host path):
    # ``seal`` stores the path it is given, so seal from the repository root.
    os.chdir(ROOT)
    v1_path = v1_path.relative_to(ROOT)
    v2_path = v2_path.relative_to(ROOT)
    assert v1_sha != v2_sha, "the two producer stubs must have distinct sha256"
    assert _sha256_text(PRODUCER_V1) == v1_sha and _sha256_text(PRODUCER_V2) == v2_sha

    # sealed_by_v1: the baseline pin. Envelope declares v1, sealed against v1.
    sealed_v1 = _seal_fixed(
        _declaring_payload(base_payload, source_sha256=v1_sha, producer_path=v1_path), v1_path
    )
    # sealed_by_v2: the swapped-in model. Envelope declares v2, sealed against v2.
    sealed_v2 = _seal_fixed(
        _declaring_payload(base_payload, source_sha256=v2_sha, producer_path=v2_path), v2_path
    )
    # cross_swap: the envelope forges the model identity back to v1 while the receipt still pins
    # v2 (the file that actually sealed it). The two identities disagree.
    sealed_cross = _seal_fixed(
        _declaring_payload(base_payload, source_sha256=v1_sha, producer_path=v1_path), v2_path
    )

    # The re-sealed stubs must carry the base payload's real gate outcome unchanged; this
    # is the guarantee the tests rely on. Check it before writing so a regression fails the build.
    for sealed in (sealed_v1, sealed_v2, sealed_cross):
        assert sealed["payload"].get("gates", {}) == gate_outcome["gates"]
        assert sealed["payload"].get("overall_pass", True) == gate_outcome["overall_pass"]

    sealed_v1_sha = _write_json(SWAP_DIR / "sealed_by_v1.json", sealed_v1)
    sealed_v2_sha = _write_json(SWAP_DIR / "sealed_by_v2.json", sealed_v2)
    sealed_cross_sha = _write_json(SWAP_DIR / "cross_swap.json", sealed_cross)

    manifest = {
        "note": (
            "Controlled-model-swap fixtures: the registry pin is the change-control "
            "point. No physiology is copied or claimed; this is contract mechanics only."
        ),
        "producer_name": PRODUCER_NAME,
        "input_sha256": input_sha256,
        # The producer gate outcome the re-sealed envelopes preserve (the contract-test base's
        # passing gate), recorded so the fixture's meaning is inspectable.
        "producer_gate_outcome": gate_outcome,
        "producers": {
            "producer_v1.py": {"path": "tests/fixtures/i1/model_swap/producer_v1.py", "sha256": v1_sha},
            "producer_v2.py": {"path": "tests/fixtures/i1/model_swap/producer_v2.py", "sha256": v2_sha},
        },
        "fixtures": {
            "sealed_by_v1": {
                "path": "tests/fixtures/i1/model_swap/sealed_by_v1.json",
                "sha256": sealed_v1_sha,
                "declared_source_sha256": v1_sha,
                "receipt_producer_sha256": v1_sha,
            },
            "sealed_by_v2": {
                "path": "tests/fixtures/i1/model_swap/sealed_by_v2.json",
                "sha256": sealed_v2_sha,
                "declared_source_sha256": v2_sha,
                "receipt_producer_sha256": v2_sha,
            },
            "cross_swap": {
                "path": "tests/fixtures/i1/model_swap/cross_swap.json",
                "sha256": sealed_cross_sha,
                "declared_source_sha256": v1_sha,
                "receipt_producer_sha256": v2_sha,
            },
        },
        "registry_pins": {
            "pin_v1": {
                "producer_name": PRODUCER_NAME,
                "producer_sha256": v1_sha,
                "input_sha256": input_sha256,
                "payload_sha256": _payload_sha(sealed_v1["payload"]),
            },
            "pin_v2": {
                "producer_name": PRODUCER_NAME,
                "producer_sha256": v2_sha,
                "input_sha256": input_sha256,
                "payload_sha256": _payload_sha(sealed_v2["payload"]),
            },
        },
    }
    manifest_path = SWAP_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"producer_v1.py sha256={v1_sha}")
    print(f"producer_v2.py sha256={v2_sha}")
    print(f"sealed_by_v1.json sha256={sealed_v1_sha}")
    print(f"sealed_by_v2.json sha256={sealed_v2_sha}")
    print(f"cross_swap.json sha256={sealed_cross_sha}")
    print(f"wrote {manifest_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
