"""Cross-module coupling checks for the 2026-09 BodyTwin build packages.

The integrated modules are exercised together, without mocks:

* I2 (``bodytwin.geometry``): region/material/mass payload built from the repo's layered box,
  structurally validated and registered against the material roster;
* I1 (``bodytwin.framework``): the payload is wrapped into a ResultEnvelope, sealed with a receipt
  bound to the producing module's source, verified, and preflighted by a consumer contract;
* the three migrated consumers (thermoregulation, cardiac_output, respiratory) run as real cells
  on the committed synthetic input and on sealed I1 fixtures;
* chained receipts: blood_oxygen_transport's output is sealed with the cardiac_output
  output it actually read as upstream, and verified against that upstream payload hash;
* I3 (``bodytwin.cells.musculoskeletal``): a paired relaxation fit is sealed in the same envelope
  and consumed under a unit/region/material contract.

Each negative case asserts the specific refusal text. Scientific acceptance needs a trusted
artifact pin, passing producer gates and scientific provenance; the helpers below pass the pin of
exactly the payload under test, and the pin itself is tested separately.

These are contract/integrity checks on synthetic or public-model inputs. They say nothing about
biological validity.
"""
import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from conftest import CELLS, EXAMPLES, REPO

SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bodytwin.cells.musculoskeletal import tendon_paired_analysis_v1 as paired  # noqa: E402
from bodytwin.framework import consumer_preflight_v1 as cp  # noqa: E402
from bodytwin.framework import region_envelope_v1 as renv  # noqa: E402
from bodytwin.framework import region_vocabulary_v1 as rvocab  # noqa: E402
from bodytwin.framework import result_envelope_v1 as rev  # noqa: E402
from bodytwin.geometry import material_roster_v1 as roster  # noqa: E402
from bodytwin.geometry.i2_payload_validate_v1 import validate_payload  # noqa: E402
from bodytwin.geometry.i2_region_registration_v1 import registration_issues  # noqa: E402

I1_FIXTURES = REPO / "tests" / "fixtures" / "i1"
REGISTRY_PATH = I1_FIXTURES / "registry" / "thermoregulation_registry.json"
EXAMPLE_BUILDER = REPO / "examples" / "geometry" / "i1_region_example.py"
I2_PRODUCER = SRC / "bodytwin" / "geometry" / "region_mass_v1.py"
I3_PRODUCER = SRC / "bodytwin" / "cells" / "musculoskeletal" / "tendon_paired_analysis_v1.py"
REGION_IDS = ("SYNTH-BOX-1", "SYNTH-BOX-2")
MATERIAL_MAP = {"SYNTH-BOX-1": "MSK-MUSCLE", "SYNTH-BOX-2": "MSK-MUSCLE"}


def _payload_sha(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------------ I2 payload -----
@pytest.fixture(scope="module")
def i2_payload():
    """Region payload built by the committed example builder (I2 modules, repo layered box)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("i1_region_example_for_xmod", EXAMPLE_BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    payload = json.loads(json.dumps(module.build()))  # plain JSON types, as consumers see it
    committed = json.loads((EXAMPLES / "geometry" / "region_example.json").read_text())
    assert payload == committed, "examples/geometry/region_example.json is stale vs its builder"
    return payload


def _calibrated(payload):
    """Same geometry, provenance relabelled calibrated, so the strict scientific path applies."""
    p = copy.deepcopy(payload)
    p["synthetic"] = False
    for region in p["regions"]:
        region["mass"]["provenance"] = "calibrated"
        region["heat_capacity"]["provenance"] = "calibrated"
    p["totals"]["mass"]["provenance"] = "calibrated"
    p["totals"]["heat_capacity"]["provenance"] = "calibrated"
    return p


def _region_preflight(payload, **kw):
    args = dict(region_ids=REGION_IDS, material_region_map=MATERIAL_MAP,
                structural_validator=validate_payload, mode="strict",
                producer_name="region_mass_v1", producer_sha256=_file_sha(I2_PRODUCER),
                input_sha256="1" * 64, current_input_sha256="1" * 64,
                expected_producer_sha256=_file_sha(I2_PRODUCER),
                expected_payload_sha256=rev.payload_sha256(payload))
    args.update(kw)
    return renv.preflight_region_payload(payload, **args)


def test_i2_payload_is_structurally_valid_and_registered(i2_payload):
    assert validate_payload(i2_payload) == []
    assert registration_issues(i2_payload) == []
    # the framework region vocabulary is the geometry roster in this repository, not a copy
    assert tuple(roster.REGION_IDS) == rvocab.I2_REGION_IDS
    assert rvocab.assert_i2_parity()[0] is True
    assert cp.THERMO_REGION_VOCABULARY == frozenset(roster.REGION_IDS)


def test_i2_payload_sealed_in_i1_envelope_verifies_and_is_accepted(i2_payload):
    envelope, verdict = _region_preflight(_calibrated(i2_payload))
    assert verdict.state == "ACCEPT_SCIENTIFIC", verdict
    sealed = rev.seal(envelope, I2_PRODUCER)
    assert rev.verify_receipt(sealed, producer_path=I2_PRODUCER,
                              expected_producer_sha256=_file_sha(I2_PRODUCER)) == envelope
    q = envelope["quantities"]["regions.SYNTH-BOX-2.mass"]
    region2 = next(r for r in i2_payload["regions"] if r["region_id"] == "SYNTH-BOX-2")
    assert q["value"] == region2["mass"]["value"] and q["units"] == "kg"
    assert q["material_region_id"] == "MSK-MUSCLE"


def test_synthetic_i2_payload_is_only_a_demo(i2_payload):
    _, strict = _region_preflight(i2_payload)
    assert strict.state == "REFUSE"
    assert any("disallowed provenance" in f for f in strict.failures), strict.failures
    _, demo = _region_preflight(i2_payload, mode="synthetic_demo")
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO", demo


# --------------------------------------------------------------- I2 -> I1 negatives -----
def test_wrong_unit_in_region_payload_is_refused(i2_payload):
    bad = _calibrated(i2_payload)
    bad["regions"][0]["mass"]["units"] = "g"
    envelope, verdict = _region_preflight(bad, structural_validator=None)
    assert verdict.state == "REFUSE"
    assert any(f.startswith("unit mismatch: regions.SYNTH-BOX-1.mass") for f in verdict.failures), \
        verdict.failures


def test_gram_mass_is_refused_by_the_geometry_validator(i2_payload):
    # closed gap 1: the validator enforces the declared kg mass contract (no silent conversion)
    bad = _calibrated(i2_payload)
    bad["regions"][0]["mass"]["units"] = "g"
    errors = validate_payload(bad)
    assert errors == ["regions[0].mass.units: expected 'kg', got 'g'"]
    total = _calibrated(i2_payload)
    total["totals"]["mass"]["units"] = "g"
    assert validate_payload(total) == ["totals.mass.units: expected 'kg', got 'g'"]
    envelope, verdict = _region_preflight(bad)
    assert envelope is None and verdict.state == "REFUSE"
    assert list(verdict.failures) == errors


def test_wrong_coordinate_frame_is_refused(i2_payload):
    framed = _calibrated(i2_payload)
    frame = {"name": "lab", "units": "mm", "origin_mm": [0.0, 0.0, 0.0]}
    for region in framed["regions"]:
        region["frame"] = frame
    # a declared frame is refused under the default scalar contract ...
    _, verdict = _region_preflight(framed, structural_validator=None)
    assert verdict.state == "REFUSE"
    assert any(f.startswith("frame mismatch: regions.SYNTH-BOX-1.mass") for f in verdict.failures)
    # ... accepted when the consumer requires exactly that frame ...
    frames = {rid: frame for rid in REGION_IDS}
    _, ok = _region_preflight(framed, structural_validator=None, frames=frames)
    assert ok.state == "ACCEPT_SCIENTIFIC", ok
    # ... and refused when the consumer requires a different frame.
    other = {rid: dict(frame, name="body") for rid in REGION_IDS}
    _, wrong = _region_preflight(framed, structural_validator=None, frames=other)
    assert any("frame mismatch" in f for f in wrong.failures), wrong.failures


def test_unknown_region_id_is_refused(i2_payload):
    payload = _calibrated(i2_payload)
    _, verdict = _region_preflight(payload, region_ids=("SYNTH-BOX-1",),
                                   material_region_map={"SYNTH-BOX-1": "MSK-MUSCLE"})
    assert verdict.state == "REFUSE"
    assert any(f.startswith("unknown region_id: regions.SYNTH-BOX-2") for f in verdict.failures), \
        verdict.failures


def test_region_material_mismatch_is_refused(i2_payload):
    payload = _calibrated(i2_payload)
    swapped = dict(MATERIAL_MAP, **{"SYNTH-BOX-1": "MSK-BONE"})
    _, verdict = _region_preflight(payload, material_region_map=swapped)
    assert verdict.state == "REFUSE"
    assert ("material region mismatch: regions.SYNTH-BOX-1.mass expected 'MSK-BONE' "
            "got 'MSK-MUSCLE'") in verdict.failures


def test_unregistered_material_region_is_refused_at_the_seam(i2_payload):
    payload = _calibrated(i2_payload)
    payload["regions"][0]["material_region_id"] = "MSK-NOPE"
    envelope, verdict = _region_preflight(payload, structural_validator=None)
    assert envelope is None and verdict.state == "REFUSE"
    assert any("material region" in f and "MSK-NOPE" in f for f in verdict.failures)


def _literal_cited(payload):
    """Same geometry, non-synthetic, every value labelled a repo literal (``literal_cited``)."""
    p = copy.deepcopy(payload)
    p["synthetic"] = False
    for region in p["regions"]:
        region["mass"]["provenance"] = "literal_cited"
        region["heat_capacity"]["provenance"] = "literal_cited"
    p["totals"]["mass"]["provenance"] = "literal_cited"
    p["totals"]["heat_capacity"]["provenance"] = "literal_cited"
    return p


def test_literal_cited_region_payload_is_not_accepted_as_scientific(i2_payload):
    # closed gap 7: a repository literal is its own provenance category, never calibrated
    assert rev.normalize_provenance("literal_cited") == "literal_cited"
    assert "literal_cited" not in rev.DEFAULT_SCIENTIFIC_ALLOWED
    assert validate_payload(_literal_cited(i2_payload)) == []
    _, verdict = _region_preflight(_literal_cited(i2_payload))
    assert verdict.state == "REFUSE"
    assert any(f.startswith("disallowed provenance: producer") and "literal_cited" in f
               for f in verdict.failures), verdict.failures
    _, demo = _region_preflight(_literal_cited(i2_payload), mode="synthetic_demo")
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"


def test_abstained_i2_mass_is_not_accepted(i2_payload):
    # The payload reports an unknown mass (e.g. missing density) as value None + missing; the envelope must refuse
    payload = _calibrated(i2_payload)
    payload["regions"][0]["mass"]["value"] = None
    payload["regions"][0]["mass"]["missing"] = ["density"]
    _, verdict = _region_preflight(payload, structural_validator=None)
    assert verdict.state == "REFUSE"
    assert "abstain: required value missing for regions.SYNTH-BOX-1.mass" in verdict.failures


def test_tampered_sealed_region_payload_is_refused(i2_payload):
    envelope, verdict = _region_preflight(_calibrated(i2_payload))
    assert verdict.accepted
    sealed = rev.seal(envelope, I2_PRODUCER)
    sealed["payload"]["quantities"]["regions.SYNTH-BOX-1.mass"]["value"] *= 2.0
    with pytest.raises(ValueError, match="stale result: payload hash mismatch"):
        rev.verify_receipt(sealed, producer_path=I2_PRODUCER)


def test_changed_producer_source_after_sealing_is_refused(i2_payload, tmp_path):
    envelope, _ = _region_preflight(_calibrated(i2_payload))
    producer = tmp_path / "region_mass_v1.py"
    shutil.copyfile(I2_PRODUCER, producer)
    sealed = rev.seal(envelope, producer)
    producer.write_text(producer.read_text() + "\n# edited after sealing\n")
    with pytest.raises(ValueError, match="stale result: producer code changed since seal"):
        rev.verify_receipt(sealed, producer_path=producer)


def test_stale_upstream_input_is_refused(i2_payload):
    # the envelope declares the input it was computed from; the consumer now sees another input
    _, verdict = _region_preflight(_calibrated(i2_payload), current_input_sha256="2" * 64)
    assert verdict.state == "REFUSE"
    assert any(f.startswith("stale input: producer input_sha256") for f in verdict.failures)


def test_failed_upstream_gates_are_refused_by_the_region_preflight(i2_payload):
    # closed gap 2: scientific acceptance refuses an upstream whose own gates failed
    for kw in ({}, {"require_producer_gates_pass": True}):
        _, verdict = _region_preflight(_calibrated(i2_payload), overall_pass=False,
                                       gates={"geometry_certificate": False}, **kw)
        assert verdict.state == "REFUSE"
        assert "producer gates failed: overall_pass is not true" in verdict.failures
    with pytest.raises(ValueError, match="strict mode cannot waive producer gates"):
        _region_preflight(_calibrated(i2_payload), overall_pass=False,
                          gates={"geometry_certificate": False},
                          require_producer_gates_pass=False)
    _, demo = _region_preflight(_calibrated(i2_payload), overall_pass=False,
                                gates={"geometry_certificate": False}, mode="synthetic_demo")
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert "producer gates failed: overall_pass is not true" in demo.flags


def test_region_scientific_acceptance_requires_the_trusted_payload_pin(i2_payload):
    payload = _calibrated(i2_payload)
    _, unpinned = _region_preflight(payload, expected_payload_sha256=None)
    assert unpinned.state == "REFUSE"
    assert ("no trusted artifact pin: expected_payload_sha256 is required for scientific "
            "acceptance") in unpinned.failures
    edited = copy.deepcopy(payload)
    edited["regions"][0]["mass"]["value"] *= 10.0
    _, verdict = _region_preflight(edited, expected_payload_sha256=rev.payload_sha256(payload))
    assert verdict.state == "REFUSE"
    assert any(f.startswith("artifact not approved: region payload sha256")
               for f in verdict.failures), verdict.failures


# ------------------------------------------------ I1 fixtures through the consumers -----
def _registry_entry():
    return cp.registry_entry(cp.load_registry(REGISTRY_PATH), "metabolic_cost_results.json")


def _fixture(name):
    return json.loads((I1_FIXTURES / name / "metabolic_cost" / "metabolic_cost_results.json")
                      .read_text())


def test_public_control_with_failed_producer_gates_is_refused():
    # closed gap 3: the real adapted producer result (overall_pass=false) is refused by default
    raw = _fixture("public_control")
    assert raw["payload"]["overall_pass"] is False
    for mode in ("auto", "strict"):
        _, verdict = cp.preflight_metabolic_cost(raw, mode=mode, registry_entry=_registry_entry())
        assert verdict.state == "REFUSE"
        assert "producer gates failed: overall_pass is not true" in verdict.failures
    _, profile = cp.preflight_metabolic_cost(raw, mode="strict", registry_entry=_registry_entry(),
                                             **rev.SCIENTIFIC_REQUIREMENTS)
    assert "producer gates failed: overall_pass is not true" in profile.failures


def test_contract_fixture_is_the_only_approved_artifact():
    raw = _fixture("contract_gates_pass")
    entry = _registry_entry()
    assert entry["payload_sha256"] == rev.payload_sha256(raw["payload"])
    _, verdict = cp.preflight_metabolic_cost(raw, mode="auto", registry_entry=entry)
    assert verdict.state == "ACCEPT_SCIENTIFIC", verdict


@pytest.mark.parametrize("variant,substring", [
    ("wrong_unit", "unit mismatch"),
    ("wrong_frame", "frame mismatch"),
    ("wrong_region", "region mismatch"),
    ("unknown_region", "unknown region_id"),
    ("stale_input", "stale input"),
    ("changed_payload", "stale result: payload hash mismatch"),
    ("changed_producer", "producer changed"),
])
def test_consumer_preflight_refuses_each_variant(variant, substring):
    _, verdict = cp.preflight_metabolic_cost(_fixture(f"variants/{variant}"), mode="strict",
                                             registry_entry=_registry_entry())
    assert verdict.state == "REFUSE"
    assert any(substring in f for f in verdict.failures), verdict.failures


def test_registry_input_pin_change_is_refused_as_stale():
    """The registry pins a different upstream input than the sealed envelope declares (the
    envelope-side mutation is the ``stale_input`` variant above)."""
    entry = dict(_registry_entry(), input_sha256="3" * 64)
    _, verdict = cp.preflight_metabolic_cost(_fixture("contract_gates_pass"), mode="strict",
                                             registry_entry=entry)
    assert verdict.state == "REFUSE"
    assert any(f.startswith("stale input") for f in verdict.failures), verdict.failures


def _edited_with_recomputed_hash():
    sealed = _fixture("contract_gates_pass")
    sealed["payload"]["quantities"]["muscle_mass.total_body_mass_kg"]["value"] = 726.0
    sealed["receipt"]["payload_sha256"] = _payload_sha(sealed["payload"])
    return sealed


def test_receipt_is_integrity_only_an_edit_with_recomputed_hash_verifies():
    """The receipt is an unkeyed sha256: it detects accidental change or a stale payload, not a
    deliberate edit whose hash is recomputed (no authenticity guarantee)."""
    sealed = _edited_with_recomputed_hash()
    assert rev.verify_receipt(sealed)["quantities"]["muscle_mass.total_body_mass_kg"]["value"] == 726.0
    untouched = _fixture("contract_gates_pass")
    untouched["payload"]["quantities"]["muscle_mass.total_body_mass_kg"]["value"] = 726.0
    with pytest.raises(ValueError, match="stale result: payload hash mismatch"):
        rev.verify_receipt(untouched)


def test_edited_payload_with_recomputed_hash_is_refused_by_the_registry_pin():
    # closed gap 8: the receipt still verifies (integrity only), but the edited payload is not
    # the artifact the trusted registry approved, so scientific acceptance refuses it
    edited = _edited_with_recomputed_hash()
    assert rev.verify_receipt(edited) is edited["payload"]
    _, verdict = cp.preflight_metabolic_cost(edited, mode="strict",
                                             registry_entry=_registry_entry())
    assert verdict.state == "REFUSE", verdict
    assert any(f.startswith("artifact not approved: payload_sha256") for f in verdict.failures), \
        verdict.failures
    assert not any("producer gates failed" in f for f in verdict.failures)


# --------------------------------------------------------------- real cell runs ------
CONSUMERS = {
    "thermoregulation": CELLS / "organ_systems" / "thermoregulation.py",
    "cardiac_output": CELLS / "cardiovascular" / "cardiac_output.py",
    "respiratory": CELLS / "respiratory" / "respiratory.py",
}


def _cell_env(out, registry=False, mode=None):
    env = dict(os.environ, BODYTWIN_OUT=str(out))
    env.pop("BODYTWIN_I1_REGISTRY", None)
    env.pop("BODYTWIN_EVIDENCE_MODE", None)
    if registry:
        env["BODYTWIN_I1_REGISTRY"] = str(REGISTRY_PATH)
    if mode:
        env["BODYTWIN_EVIDENCE_MODE"] = mode
    return env


def _run_cell(script, env, timeout=300):
    return subprocess.run([sys.executable, script.name], cwd=script.parent, env=env,
                          capture_output=True, text=True, timeout=timeout)


def _stage(out, fixture_dir):
    shutil.copytree(fixture_dir, out, dirs_exist_ok=True)


@pytest.mark.parametrize("variant,substring", [
    ("wrong_unit", "unit mismatch"),
    ("wrong_frame", "frame mismatch"),
    ("unknown_region", "unknown region_id"),
])
@pytest.mark.parametrize("name", sorted(CONSUMERS))
def test_consumer_cell_refuses_bad_contract_before_using_numbers(name, variant, substring,
                                                                 tmp_path):
    _stage(tmp_path, I1_FIXTURES / "variants" / variant)
    proc = _run_cell(CONSUMERS[name], _cell_env(tmp_path, registry=True))
    assert proc.returncode == 2, (proc.stdout + proc.stderr)[-2000:]
    assert f"REFUSED: {substring}" in proc.stdout
    assert not (tmp_path / name / f"{name}_results.json").exists()
    refused = json.loads((tmp_path / name / f"{name}_refused.json").read_text())
    assert any(substring in f for f in refused["failures"])


@pytest.mark.parametrize("name", sorted(CONSUMERS))
def test_consumer_cell_refuses_stale_input_strict(name, tmp_path):
    _stage(tmp_path, I1_FIXTURES / "variants" / "stale_input")
    proc = _run_cell(CONSUMERS[name], _cell_env(tmp_path, registry=True, mode="strict"))
    assert proc.returncode == 2, (proc.stdout + proc.stderr)[-2000:]
    assert "REFUSED: stale input" in proc.stdout


@pytest.mark.parametrize("name", sorted(CONSUMERS))
def test_consumer_cell_labels_synthetic_demo(name, tmp_path):
    _stage(tmp_path, EXAMPLES / "synthetic_inputs")
    proc = _run_cell(CONSUMERS[name], _cell_env(tmp_path))
    assert proc.returncode == 0, (proc.stdout + proc.stderr)[-2000:]
    assert "EVIDENCE: SYNTHETIC DEMO" in proc.stdout
    report = json.loads((tmp_path / name / f"{name}_results.json").read_text())
    assert report["evidence_class"] == "synthetic_demo"
    assert report["scientific_claim"] is False


@pytest.mark.parametrize("name", sorted(CONSUMERS))
def test_consumer_cell_refuses_failed_producer_gates(name, tmp_path):
    # closed gap 4: no consumer cell makes a scientific claim on a gate-failed upstream
    _stage(tmp_path, I1_FIXTURES / "public_control")
    proc = _run_cell(CONSUMERS[name], _cell_env(tmp_path, registry=True))
    assert proc.returncode == 2, (proc.stdout + proc.stderr)[-2000:]
    assert "REFUSED: producer gates failed: overall_pass is not true" in proc.stdout
    assert "EVIDENCE: SCIENTIFIC" not in proc.stdout
    assert not (tmp_path / name / f"{name}_results.json").exists()
    refused = json.loads((tmp_path / name / f"{name}_refused.json").read_text())
    assert "producer gates failed: overall_pass is not true" in refused["failures"]


# ------------------------------------------------------- chained receipts ---------------
@pytest.fixture(scope="module")
def cell_chain(tmp_path_factory):
    """cardiac_output -> blood_oxygen_transport on the committed synthetic inputs."""
    out = tmp_path_factory.mktemp("xmod_chain")
    _stage(out, EXAMPLES / "synthetic_inputs")
    env = _cell_env(out)
    runs = {}
    for rel in ("cardiovascular/cardiac_output.py", "respiratory/respiratory.py",
                "respiratory/pulmonary_gas_exchange.py", "haematology/blood_oxygen_transport.py"):
        proc = _run_cell(CELLS / rel, env)
        runs[rel] = proc
        assert proc.returncode == 0, (rel, (proc.stdout + proc.stderr)[-2000:])
    return out


def _seal_chain_files(out, sealed_dir):
    co_path = out / "cardiac_output" / "cardiac_output_results.json"
    bo_path = out / "blood_oxygen_transport" / "blood_oxygen_transport_results.json"
    co_script = CELLS / "cardiovascular" / "cardiac_output.py"
    bo_script = CELLS / "haematology" / "blood_oxygen_transport.py"
    upstream = rev.seal(json.loads(co_path.read_text()), co_script)
    downstream = rev.seal_chain(json.loads(bo_path.read_text()), bo_script, upstream)
    sealed_dir.mkdir(parents=True, exist_ok=True)
    (sealed_dir / "cardiac_output.sealed.json").write_text(json.dumps(upstream))
    down_file = sealed_dir / "blood_oxygen_transport.sealed.json"
    down_file.write_text(json.dumps(downstream))
    return co_path, bo_script, upstream, downstream, down_file


def test_chain_verifies_against_the_upstream_actually_read(cell_chain, tmp_path):
    co_path, bo_script, upstream, downstream, down_file = _seal_chain_files(cell_chain, tmp_path)
    read_hash = _payload_sha(json.loads(co_path.read_text()))
    assert downstream["receipt"]["upstream"]["payload_sha256"] == read_hash
    payload = rev.read_verified_chain(down_file, bo_script,
                                      expected_upstream_payload_sha256=read_hash)
    assert payload == downstream["payload"]


def test_chain_refuses_an_upstream_changed_after_sealing(cell_chain, tmp_path):
    co_path, bo_script, _, _, down_file = _seal_chain_files(cell_chain, tmp_path)
    changed = json.loads(co_path.read_text())
    changed["i1_changed_after_seal"] = True
    with pytest.raises(ValueError, match="stale upstream"):
        rev.read_verified_chain(down_file, bo_script,
                                expected_upstream_payload_sha256=_payload_sha(changed))


def test_chain_refuses_a_tampered_downstream_payload(cell_chain, tmp_path):
    co_path, bo_script, _, downstream, down_file = _seal_chain_files(cell_chain, tmp_path)
    tampered = copy.deepcopy(downstream)
    tampered["payload"]["i1_tampered"] = True
    down_file.write_text(json.dumps(tampered))
    with pytest.raises(ValueError, match="stale result: payload hash mismatch"):
        rev.read_verified_chain(down_file, bo_script,
                                expected_upstream_payload_sha256=_payload_sha(
                                    json.loads(co_path.read_text())))


def test_chain_refuses_a_receipt_without_upstream(cell_chain, tmp_path):
    bo_path = cell_chain / "blood_oxygen_transport" / "blood_oxygen_transport_results.json"
    bo_script = CELLS / "haematology" / "blood_oxygen_transport.py"
    unchained = tmp_path / "unchained.json"
    unchained.write_text(json.dumps(rev.seal(json.loads(bo_path.read_text()), bo_script)))
    with pytest.raises(ValueError, match="missing upstream receipt"):
        rev.read_verified_chain(unchained, bo_script, expected_upstream_payload_sha256="0" * 64)


def test_chain_refuses_an_upstream_with_failed_gates(tmp_path):
    # closed gap 5: a chained result cannot build on an upstream whose own gates failed
    raw = _fixture("public_control")
    upstream_producer = tmp_path / "upstream.py"
    upstream_producer.write_text("# upstream stand-in\n")
    upstream = rev.seal(raw["payload"], upstream_producer)
    assert upstream["payload"]["overall_pass"] is False
    down_producer = tmp_path / "downstream.py"
    down_producer.write_text("# downstream stand-in\n")
    with pytest.raises(ValueError, match="producer gates failed: upstream envelope"):
        rev.seal_chain({"value": 1.0}, down_producer, upstream)
    sealed = rev.seal_chain({"value": 1.0}, down_producer, upstream, allow_failed_upstream=True)
    path = tmp_path / "down.json"
    path.write_text(json.dumps(sealed))
    with pytest.raises(ValueError, match="producer gates failed: the chained upstream"):
        rev.read_verified_chain(path, down_producer,
                                expected_upstream_payload_sha256=_payload_sha(raw["payload"]))


def test_chain_without_expected_upstream_identity_is_refused(cell_chain, tmp_path):
    # closed gap 9: a bare upstream link is not evidence; a foreign upstream is caught by the
    # independently expected upstream hash
    co_path, bo_script, upstream, _, down_file = _seal_chain_files(cell_chain, tmp_path)
    with pytest.raises(ValueError, match="expected upstream identity required"):
        rev.read_verified_chain(down_file, bo_script)
    foreign = copy.deepcopy(upstream)
    foreign["payload"]["i1_foreign_upstream"] = True
    bo_path = cell_chain / "blood_oxygen_transport" / "blood_oxygen_transport_results.json"
    foreign_file = tmp_path / "foreign_chain.json"
    foreign_file.write_text(json.dumps(rev.seal_chain(json.loads(bo_path.read_text()), bo_script,
                                                      foreign)))
    with pytest.raises(ValueError, match="stale upstream"):
        rev.read_verified_chain(foreign_file, bo_script,
                                expected_upstream_payload_sha256=_payload_sha(
                                    json.loads(co_path.read_text())))


# ------------------------------------------------------------ I3 paired analysis ------
AREA_MM2 = 57.5
REF_MM = 200.0
TAU_TRUE = 30.0
TAU_SPEC = {"tau_s": {"units": "s", "region_id": "MSK-TENDON", "frame": "scalar",
                      "time": None, "required": True, "material_region_id": "MSK-TENDON"}}


def _paired_inputs(t, eps):
    stress = paired._sls_stress_from_strain(t, eps, 0.8e9, 0.4e9, TAU_TRUE)
    return paired.PairedInputs(
        time_s=t, force_gf=stress * AREA_MM2 * 1e-6 / paired.m.GF_TO_N,
        displacement_mm=eps * REF_MM, reference_length_mm=REF_MM,
        area_mm2=AREA_MM2, specimen="syn", donor="syn")


def _clean_step():
    t = np.linspace(0.0, 300.0, 3001)
    return t, np.where(t <= 1.0, 0.02 * t, 0.02)


def _unloading_step():
    t = np.linspace(0.0, 60.0, 6001)
    eps = np.where(t <= 1.0, 0.02 * t, 0.02)
    return t, np.where(t > 40.0, 0.02 * np.clip((60.0 - t) / 20.0, 0, 1), eps)


def _failure_ramp_step():
    t = np.linspace(0.0, 25.0, 2501)
    return t, 0.04 * t / 25.0


def _tau_envelope(fit, provenance="calibrated"):
    result = {"tau_s": fit.get("tau_s")}
    spec = {k: dict(v) for k, v in TAU_SPEC.items()}
    for v in spec.values():
        v.pop("required")
    return rev.make_envelope(
        result, producer_name="tendon_paired_analysis_v1",
        producer_sha256=_file_sha(I3_PRODUCER), provenance=provenance,
        input_sha256="4" * 64, quantity_specs=spec, regime="resting",
        gates={"step_valid": fit["state"] == "PASS"}, overall_pass=fit["state"] == "PASS")


def _tau_verdict(envelope, **kw):
    args = dict(expected_quantities=TAU_SPEC,
                allowed_provenance=rev.allowed_provenance_for(), mode="strict",
                current_input_sha256="4" * 64,
                expected_producer_sha256=_file_sha(I3_PRODUCER),
                expected_regions=frozenset(roster.REGION_IDS),
                expected_material_regions={"tau_s": "MSK-TENDON"},
                material_region_vocabulary=frozenset(rvocab.I2_MATERIAL_REGION_IDS),
                expected_payload_sha256=rev.payload_sha256(envelope))
    args.update(kw)
    return rev.validate_for_consumer(envelope, **args)


def test_i3_paired_fit_sealed_in_envelope_is_accepted():
    fit = paired.fit_tau_from_pairs(_paired_inputs(*_clean_step()))
    assert fit["state"] == "PASS" and abs(fit["tau_s"] - TAU_TRUE) / TAU_TRUE < 1e-6
    envelope = _tau_envelope(fit)
    sealed = rev.seal(envelope, I3_PRODUCER)
    payload = rev.verify_receipt(sealed, producer_path=I3_PRODUCER,
                                 expected_producer_sha256=_file_sha(I3_PRODUCER))
    verdict = _tau_verdict(payload, **rev.SCIENTIFIC_REQUIREMENTS | {
        "require_uncertainty": False, "require_provenance_detail": False})
    assert verdict.state == "ACCEPT_SCIENTIFIC", verdict


def test_i3_tau_with_wrong_unit_or_region_is_refused():
    fit = paired.fit_tau_from_pairs(_paired_inputs(*_clean_step()))
    envelope = _tau_envelope(fit)
    envelope["quantities"]["tau_s"]["units"] = "ms"
    verdict = _tau_verdict(envelope)
    assert "unit mismatch: tau_s expected 's' got 'ms'" in verdict.failures
    envelope = _tau_envelope(fit)
    envelope["quantities"]["tau_s"]["region_id"] = "MSK-MUSCLE"
    envelope["quantities"]["tau_s"]["material_region_id"] = "MSK-MUSCLE"
    verdict = _tau_verdict(envelope)
    assert ("region mismatch: tau_s expected 'MSK-TENDON' got 'MSK-MUSCLE'" in verdict.failures)
    assert ("material region mismatch: tau_s expected 'MSK-TENDON' got 'MSK-MUSCLE'"
            in verdict.failures)


@pytest.mark.parametrize("step,reason", [
    (_unloading_step, "unloading_after_peak"),
    (_failure_ramp_step, "no_hold"),
])
def test_i3_bad_step_is_rejected_and_never_reaches_a_consumer_as_a_value(step, reason):
    fit = paired.fit_tau_from_pairs(_paired_inputs(*step()))
    assert fit["state"] == "REJECT", fit
    assert reason in fit["reasons"]
    assert "tau_s" not in fit
    envelope = _tau_envelope(fit)
    assert envelope["overall_pass"] is False
    verdict = _tau_verdict(envelope)
    assert verdict.state == "REFUSE"
    assert "abstain: required value missing for tau_s" in verdict.failures


def test_i3_unguarded_bad_step_is_refused_by_the_producer_gate():
    # guards=False reproduces the pre-guard behaviour (a tau on an unloading step); the
    # envelope carries the step verdict, so scientific acceptance refuses it by default.
    inputs = _paired_inputs(*_unloading_step())
    unguarded = paired.fit_tau_from_pairs(inputs, guards=False)
    assert unguarded["state"] == "PASS" and "tau_s" in unguarded
    guarded = paired.fit_tau_from_pairs(inputs)
    envelope = _tau_envelope(dict(unguarded, state=guarded["state"]))
    strict = _tau_verdict(envelope)
    assert strict.state == "REFUSE"
    assert "producer gates failed: overall_pass is not true" in strict.failures


def test_i3_tampered_sealed_tau_is_refused():
    fit = paired.fit_tau_from_pairs(_paired_inputs(*_clean_step()))
    sealed = rev.seal(_tau_envelope(fit), I3_PRODUCER)
    sealed["payload"]["quantities"]["tau_s"]["value"] = 52.44
    with pytest.raises(ValueError, match="stale result: payload hash mismatch"):
        rev.verify_receipt(sealed, producer_path=I3_PRODUCER)


# ------------------------------------------ material citations vs the cited cells ------
def _citation_rows():
    from bodytwin.geometry import material_roster_v1 as mr
    from bodytwin.geometry import region_mass_v1 as rm
    from bodytwin.geometry import thermal_capacity_v1 as tc

    rows = [(r.source, r.source_sha256, r.value)
            for r in list(rm.DENSITY_REGISTRY.values()) + list(rm.SPECIFIC_HEAT_REGISTRY.values())]
    for source in (tc.LATENT_HEAT_SOURCE, tc.SPECIFIC_HEAT_SOURCE):
        head = source.split(" ", 1)[0]
        value = (tc.LATENT_HEAT_VAPORIZATION_SWEAT_J_PER_G if "LATENT" in head
                 else rm.SPECIFIC_HEAT_REGISTRY["WHOLE-BODY"].value)
        rows.append((head, None, value))
    absent = [mr.lookup("WHOLE-BODY", p) for p in ("thermal_conductivity",
                                                   "perfusion_heat_transport")]
    return rows, absent


def _resolve_symbol(tree, symbol):
    """Value of ``NAME``, ``NAME["k"]["f"]`` (dict literal / dict(...)) or ``func(param)``."""
    import ast
    import re

    m = re.fullmatch(r"(\w+)\((\w+)\)", symbol)
    if m:
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == m.group(1):
                args = node.args.args[len(node.args.args) - len(node.args.defaults):]
                for arg, default in zip(args, node.args.defaults):
                    if arg.arg == m.group(2):
                        return ast.literal_eval(default)
        raise AssertionError(f"no default {symbol}")
    name = re.match(r"\w+", symbol).group(0)
    keys = re.findall(r'\["(\w+)"\]', symbol)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets):
            value = node.value
            for key in keys:
                if isinstance(value, ast.Call):
                    value = {kw.arg: kw.value for kw in value.keywords}[key]
                else:
                    value = dict(zip([ast.literal_eval(k) for k in value.keys],
                                     value.values))[key]
            return ast.literal_eval(value)
    raise AssertionError(f"symbol {symbol} not found")


def test_material_citations_resolve_in_the_current_tree():
    """Each ``<path>::<symbol>`` citation names a symbol in the current file whose literal equals
    the registered value, and the pinned sha256 matches the current file. Works without any
    git history (the snapshot repository has no base commit)."""
    import ast

    rows, absent = _citation_rows()
    assert len(rows) == 8
    for source, sha, value in rows:
        path, symbol = source.split("::", 1)
        data = (REPO / path).read_bytes()
        if sha is not None:
            assert hashlib.sha256(data).hexdigest() == sha, f"re-pin {source}"
        got = _resolve_symbol(ast.parse(data.decode()), symbol)
        assert got == value, (source, got, value)
    for record in absent:
        path = record["source"].split(" ", 1)[0]
        assert hashlib.sha256((REPO / path).read_bytes()).hexdigest() == record["source_sha256"]
