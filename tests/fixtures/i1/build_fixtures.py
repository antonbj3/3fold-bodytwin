#!/usr/bin/env python
"""Build the I1 preflight fixtures from the registered read-only bases (deterministic).

Run from the repository root::

    BODYTWIN_REF_ADAPTED_PRODUCER=<path to metabolic_cost_adapted.py> PYTHONPATH=src \
        python tests/fixtures/i1/build_fixtures.py

Inputs (READ-ONLY):
  * the payload of the adapted public gait2392 metabolic-cost producer (not part of this
    repository): ``out/metabolic_cost_results.json`` next to the producer script
    ``metabolic_cost_adapted.py`` (sha256 6bf7599c...), which ``BODYTWIN_REF_ADAPTED_PRODUCER``
    points to. The receipt records only the script's file name, never a host path;
  * committed synthetic fixture ``examples/synthetic_inputs/metabolic_cost/metabolic_cost_results.json``
    (sha256 4753cf5c...), kept as an *unsealed* raw synthetic demo;
  * ``tests/fixtures/i1/registry/thermoregulation_registry.json`` (producer/input test pins).

Output layout (below ``tests/fixtures/i1/``)::

    public_control/metabolic_cost/metabolic_cost_results.json       # SEALED, real producer gates
    contract_gates_pass/metabolic_cost/metabolic_cost_results.json  # SEALED, contract test only
    synthetic_demo/metabolic_cost/metabolic_cost_results.json       # raw synthetic, unsealed
    variants/<name>/metabolic_cost/metabolic_cost_results.json      # one mutation each
    manifest.json                                                   # name -> expectations

``public_control`` carries the adapted producer's real gate outcome (``overall_pass=false``) and is
refused by a scientific consumer. ``contract_gates_pass`` is the same payload with the producer
gate outcome REPLACED by a single passing ``contract_test_fixture`` gate. It exists only so the
tests can exercise the accepting path of the contract; it is not the producer's result and not
evidence. The variants are mutations of ``contract_gates_pass`` so each carries one defect.

Registry admission: this builder is the only writer of the test registry's ``payload_sha256``
(the approved artifact, here ``contract_gates_pass``); the change becomes effective only through a
reviewed commit of ``registry/thermoregulation_registry.json``.

Determinism: every sealed file uses a fixed ``sealed_utc`` (the receipt hash covers only the
payload, so overriding the timestamp does not weaken the payload/sealed-producer binding).
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

from bodytwin.framework.consumer_preflight_v1 import THERMO_QUANTITY_SPECS  # noqa: E402
from bodytwin.framework.result_envelope_v1 import make_envelope, seal  # noqa: E402

# Registered read-only bases.
_PRODUCER = os.environ.get("BODYTWIN_REF_ADAPTED_PRODUCER")
if not _PRODUCER or not Path(_PRODUCER).is_file():
    raise SystemExit("set BODYTWIN_REF_ADAPTED_PRODUCER to the adapted producer script")
PRODUCER_DIR = Path(_PRODUCER).resolve().parent
ADAPTED_OUT = PRODUCER_DIR / "out" / "metabolic_cost_results.json"
# Relative to PRODUCER_DIR; ``seal`` records it verbatim, so the build seals with that cwd.
ADAPTED_PRODUCER = Path(Path(_PRODUCER).name)
SYNTHETIC_FIXTURE = ROOT / "examples" / "synthetic_inputs" / "metabolic_cost" / "metabolic_cost_results.json"
REGISTRY_PATH = HERE / "registry" / "thermoregulation_registry.json"
FIXTURES = HERE

FIXED_SEALED_UTC = "2026-09-16T00:00:00+00:00"

REGISTRY = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
ENTRY = REGISTRY["entries"]["metabolic_cost_results.json"]

PRODUCER_NAME = ENTRY["producer_name"]
PRODUCER_SHA256 = ENTRY["producer_sha256"]
INPUT_SHA256 = ENTRY["input_sha256"]
PROVENANCE = ENTRY["provenance"]


def _sha256_payload(payload) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _seal_fixed(payload, producer_path, when=FIXED_SEALED_UTC):
    """Seal and pin ``sealed_utc`` so the fixture bytes are reproducible."""
    cwd = os.getcwd()
    os.chdir(PRODUCER_DIR)
    try:
        sealed = seal(payload, producer_path)
    finally:
        os.chdir(cwd)
    sealed["receipt"]["sealed_utc"] = when
    return sealed


def _write_json(path: Path, obj) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    return _sha256_file(path)


def _build_base_envelope():
    raw = json.loads(ADAPTED_OUT.read_text(encoding="utf-8"))
    # Propagate the adapted producer's OWN gate outcome. The healthy public control
    # reports ``overall_pass=false`` (``speed_sanity_ok=false``, ``expected_direction_pass=false``);
    # ``make_envelope`` used to hardcode ``gates={}`` / ``overall_pass=true``, silently reporting
    # success. Read the payload's real ``gates``/``overall_pass`` so ``public_control`` carries the
    # true outcome. ``make_envelope`` keeps its historical
    # defaults for payloads that lack the fields (``None`` -> ``{}`` / ``True``).
    raw_gates = raw.get("gates")
    raw_overall_pass = raw.get("overall_pass")
    return make_envelope(
        raw,
        producer_name=PRODUCER_NAME,
        producer_sha256=PRODUCER_SHA256,
        provenance=PROVENANCE,
        input_sha256=INPUT_SHA256,
        quantity_specs=THERMO_QUANTITY_SPECS,
        command=str(ADAPTED_PRODUCER),
        seed=None,
        regime="level_walking_steady_state",
        gates=raw_gates if isinstance(raw_gates, dict) else None,
        overall_pass=raw_overall_pass if isinstance(raw_overall_pass, bool) else None,
    )


# --------------------------------------------------------------------------- mutations -----
def _mutate_wrong_unit(env):
    env["quantities"]["headline.umberger2010_gross_w_per_kg"]["units"] = "W"


def _mutate_wrong_region(env):
    env["quantities"]["muscle_mass.summed_muscle_mass_kg"]["region_id"] = "WHOLE-BODY"


def _mutate_unknown_region(env):
    env["quantities"]["muscle_mass.summed_muscle_mass_kg"]["region_id"] = "MSK-NOT-A-REGION"


def _mutate_wrong_frame(env):
    env["quantities"]["headline.umberger2010_gross_w_per_kg"]["frame"] = "instantaneous"


def _mutate_wrong_time(env):
    env["quantities"]["headline.umberger2010_gross_w_per_kg"]["time"] = None


def _mutate_stale_input(env):
    env["producer"]["input_sha256"] = "0" * 64


def _mutate_changed_producer(env):
    env["producer"]["source_sha256"] = "f" * 64


def _mutate_missing_provenance(env):
    env["producer"]["provenance"] = None
    for quantity in env["quantities"].values():
        quantity["provenance"] = None


def _mutate_synthetic_in_strict(env):
    env["producer"]["provenance"] = "synthetic"
    for quantity in env["quantities"].values():
        quantity["provenance"] = "synthetic"


def _mutate_stale_check_disabled(env):
    env["validity"]["stale_if_input_sha256_mismatch"] = False


def _mutate_wrong_regime(env):
    # A recognised out-of-domain exercise regime: the contract knows it, but it is not the
    # steady-state walking regime the thermoregulation consumer accepts.
    env["validity"]["regime"] = "sprint_anaerobic"


def _mutate_missing_regime(env):
    env["validity"]["regime"] = None


# name -> (mutator, expected state, expected failure substring, needs_reseal)
VARIANTS = (
    ("wrong_unit", _mutate_wrong_unit, "REFUSE", "unit mismatch", True),
    ("wrong_region", _mutate_wrong_region, "REFUSE", "region mismatch", True),
    ("unknown_region", _mutate_unknown_region, "REFUSE", "unknown region_id", True),
    ("wrong_frame", _mutate_wrong_frame, "REFUSE", "frame mismatch", True),
    ("wrong_time", _mutate_wrong_time, "REFUSE", "time mismatch", True),
    ("wrong_regime", _mutate_wrong_regime, "REFUSE", "regime mismatch", True),
    ("missing_regime", _mutate_missing_regime, "REFUSE", "missing regime", True),
    ("stale_input", _mutate_stale_input, "REFUSE", "stale input", True),
    ("changed_producer", _mutate_changed_producer, "REFUSE", "producer changed", True),
    ("missing_provenance", _mutate_missing_provenance, "REFUSE", "disallowed provenance", True),
    ("changed_payload", None, "REFUSE", "stale result", False),
    ("synthetic_in_strict", _mutate_synthetic_in_strict, "REFUSE", "disallowed provenance", True),
    ("stale_check_disabled", _mutate_stale_check_disabled, "REFUSE", "stale check disabled", True),
)


def main():
    FIXTURES.mkdir(parents=True, exist_ok=True)
    manifest = {}

    # ---- public_control: the real adapted producer outcome (its own gates failed) ----
    real_env = _build_base_envelope()
    real_sealed = _seal_fixed(real_env, ADAPTED_PRODUCER)
    control_path = FIXTURES / "public_control" / "metabolic_cost" / "metabolic_cost_results.json"
    manifest["public_control"] = {
        "path": str(control_path.relative_to(ROOT)),
        "sha256": _write_json(control_path, real_sealed),
        "mode": "auto",
        "expected_state": "REFUSE",
        "expected_failure_substring": "producer gates failed",
        "note": "real producer gate outcome (overall_pass=false); refused for scientific use",
    }

    # ---- contract_gates_pass: contract-mechanics fixture, gate outcome replaced ----
    base_env = copy.deepcopy(real_env)
    base_env["gates"] = {"contract_test_fixture": True}
    base_env["overall_pass"] = True
    base_sealed = _seal_fixed(base_env, ADAPTED_PRODUCER)
    valid_path = FIXTURES / "contract_gates_pass" / "metabolic_cost" / "metabolic_cost_results.json"
    manifest["contract_gates_pass"] = {
        "path": str(valid_path.relative_to(ROOT)),
        "sha256": _write_json(valid_path, base_sealed),
        "mode": "auto",
        "expected_state": "ACCEPT_SCIENTIFIC",
        "expected_failure_substring": None,
        "note": ("contract test only: the producer gate outcome is replaced by a passing "
                 "contract_test_fixture gate; not the producer's result, not evidence"),
    }

    # ---- registry admission: pin the approved artifact ----
    REGISTRY["entries"]["metabolic_cost_results.json"]["payload_sha256"] = \
        _sha256_payload(base_sealed["payload"])
    REGISTRY_PATH.write_text(json.dumps(REGISTRY, indent=2) + "\n", encoding="utf-8")

    # ---- variants ----
    for name, mutator, state, substring, reseal in VARIANTS:
        if name == "changed_payload":
            # mutate one quantity value inside the sealed payload and DO NOT re-seal:
            # the receipt payload hash must catch it.
            sealed = copy.deepcopy(base_sealed)
            sealed["payload"]["quantities"]["muscle_mass.total_body_mass_kg"]["value"] += 1.0
        else:
            envelope = copy.deepcopy(base_env)
            mutator(envelope)
            sealed = _seal_fixed(envelope, ADAPTED_PRODUCER) if reseal else sealed
        path = FIXTURES / "variants" / name / "metabolic_cost" / "metabolic_cost_results.json"
        sha = _write_json(path, sealed)
        manifest[f"variants/{name}"] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": sha,
            "mode": "auto",
            "expected_state": state,
            "expected_failure_substring": substring,
        }

    # ---- synthetic_demo: byte copy of the committed synthetic fixture, unsealed ----
    if _sha256_file(SYNTHETIC_FIXTURE) != "4753cf5cef0738515bd6ffa32a54d5e041eb6e24103e5d15680d5b93a7c89357":
        raise SystemExit("committed synthetic fixture hash does not match the registered base")
    synth_path = FIXTURES / "synthetic_demo" / "metabolic_cost" / "metabolic_cost_results.json"
    synth_path.parent.mkdir(parents=True, exist_ok=True)
    synth_path.write_bytes(SYNTHETIC_FIXTURE.read_bytes())
    manifest["synthetic_demo"] = {
        "path": str(synth_path.relative_to(ROOT)),
        "sha256": _sha256_file(synth_path),
        "mode": "auto",
        "expected_state": "ACCEPT_SYNTHETIC_DEMO",
        "expected_failure_substring": None,
    }

    manifest_path = FIXTURES / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    for name, entry in manifest.items():
        print(f"{name:28s} {entry['expected_state']:22s} sha256={entry['sha256'][:16]}...")
    print(f"wrote {manifest_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
