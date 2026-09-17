"""Acceptance tests for the region-envelope adapter (``region_envelope_v1``).

These tests cover the region-payload -> envelope seam: ``bodytwin.framework`` owns the
``ResultEnvelope``; ``bodytwin.geometry`` owns the
``i2_region_example_v1`` payload. They drive I2's example payload (``examples/geometry/region_example.json``)
and I2's structural validator (``bodytwin.geometry.i2_payload_validate_v1``) from this
repository. Both are part of the same commit series, so the recorded hashes are asserted.
"""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
EXAMPLE_PATH = REPO / "examples" / "geometry" / "region_example.json"
VALIDATOR_PATH = REPO / "src" / "bodytwin" / "geometry" / "i2_payload_validate_v1.py"

# Recorded I2 source hashes. In this repository the I2 files are committed together with the
# seam, so a hash change is a deliberate re-pin and the pin tests fail instead of skipping.
EXAMPLE_SHA256 = "9bdc4be96021fe8c7776f12b8d24a4f23d6263d6f5a131e61ca3a93b6ed05e9c"
VALIDATOR_SHA256 = "268ca3854ac28671fd7cacd1d82a04f1d08f49af427ce2d9d6df8a330ec87042"

import bodytwin.framework.region_envelope_v1 as rev2

import bodytwin.framework.result_envelope_v1 as base  # noqa: F401

import bodytwin.framework.consumer_preflight_v1 as cp  # noqa: F401

import bodytwin.framework.region_vocabulary_v1 as rvf

REGION_IDS = ("SYNTH-BOX-1", "SYNTH-BOX-2")
BARE_NUMBER_SUBSTRING = "bare numeric mass/heat value is not allowed"


# --------------------------------------------------------------------------- helpers -------
def _load_example():
    if not EXAMPLE_PATH.exists():
        pytest.fail(f"region example missing: {EXAMPLE_PATH}")
    return json.loads(EXAMPLE_PATH.read_text(encoding="utf-8"))


def _load_i2_validator():
    if not VALIDATOR_PATH.exists():
        pytest.fail(f"structural validator missing: {VALIDATOR_PATH}")
    from bodytwin.geometry.i2_payload_validate_v1 import validate_payload

    return validate_payload


def _frozen_hash_matches(path, expected):
    """True iff ``path`` still hashes to the recorded freeze (and exists)."""
    if not path.exists():
        return False
    return hashlib.sha256(path.read_bytes()).hexdigest() == expected


def _fail_if_pin_moved(path, expected, label):
    """Fail when a pinned region-payload file no longer matches its recorded sha256.

    The example and the validator are part of this repository, so a changed hash is a real change
    that needs a deliberate re-pin, never a skip."""
    if not _frozen_hash_matches(path, expected):
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
        pytest.fail(
            f"I2 {label} changed since the recorded hash "
            f"({actual[:12] if actual else 'absent'} != {expected[:12]}); re-pin deliberately"
        )


def _all_calibrated_copy():
    payload = _load_example()
    payload["synthetic"] = False
    for region in payload["regions"]:
        region["mass"]["provenance"] = "calibrated"
        region["heat_capacity"]["provenance"] = "calibrated"
    payload["totals"]["mass"]["provenance"] = "calibrated"
    payload["totals"]["heat_capacity"]["provenance"] = "calibrated"
    return payload


def _preflight(payload, **kwargs):
    """``preflight_region_payload`` with the trusted pin of exactly this payload (scientific
    acceptance requires a pin; the pin itself is tested separately)."""
    kwargs.setdefault("expected_payload_sha256", base.payload_sha256(payload))
    return rev2.preflight_region_payload(payload, **kwargs)


def _make(payload, **overrides):
    args = dict(
        producer_name="i2_region_example",
        producer_sha256=None,
        input_sha256=None,
    )
    args.update(overrides)
    return rev2.make_region_envelope(payload, **args)


# ------------------------------------------------------------------- frozen source evidence -
def test_i2_example_hash_matches_the_recorded_frozen_source():
    _fail_if_pin_moved(EXAMPLE_PATH, EXAMPLE_SHA256, "example")


def test_i2_validator_hash_matches_the_recorded_frozen_source():
    _fail_if_pin_moved(VALIDATOR_PATH, VALIDATOR_SHA256, "structural validator")


# --------------------------------------------------------------------- make_region_envelope -
def test_make_region_envelope_provenance_and_quantity_contract():
    envelope, errors = _make(_load_example())
    assert errors == []
    assert envelope is not None
    assert envelope["producer"]["provenance"] == "synthetic"
    assert envelope["validity"]["regime"] is None

    quantities = envelope["quantities"]
    assert set(quantities) == {
        "regions.SYNTH-BOX-1.mass",
        "regions.SYNTH-BOX-1.heat_capacity",
        "regions.SYNTH-BOX-2.mass",
        "regions.SYNTH-BOX-2.heat_capacity",
        "totals.mass",
        "totals.heat_capacity",
        # The payload's area block is carried additively under ``area.<kind>``.
        "area.material_interface_area",
    }

    mass = quantities["regions.SYNTH-BOX-1.mass"]
    assert mass["units"] == "kg"
    assert mass["region_id"] == "SYNTH-BOX-1"
    assert mass["frame"] == "scalar"
    assert mass["time"] is None
    assert mass["value"] == pytest.approx(0.03816)
    assert quantities["regions.SYNTH-BOX-1.heat_capacity"]["units"] == "J/K"
    assert quantities["regions.SYNTH-BOX-2.heat_capacity"]["region_id"] == "SYNTH-BOX-2"

    assert quantities["totals.mass"]["units"] == "kg"
    assert quantities["totals.mass"]["region_id"] is None
    assert quantities["totals.heat_capacity"]["units"] == "J/K"
    assert quantities["totals.heat_capacity"]["region_id"] is None

    # The mirrored area block is a region-less, scalar quantity tagged by its kind in the
    # dotted key; the frozen I2 example declares no ``area.provenance``, so the adapter stores the
    # canonical "unknown" (which a consumer that explicitly expects the area then refuses) without
    # touching the top-level provenance (still "synthetic" here).
    area = quantities["area.material_interface_area"]
    assert area["value"] == pytest.approx(3600.0)
    assert area["units"] == "mm2"
    assert area["region_id"] is None
    assert area["frame"] == "scalar"
    assert area["time"] is None
    assert area["provenance"] == "unknown"
    assert envelope["producer"]["provenance"] == "synthetic"


def test_wrong_schema_is_refused_with_a_schema_message():
    envelope, errors = _make({"schema": "not_i2_region_example_v1"})
    assert envelope is None
    assert errors and any("schema" in e for e in errors), errors

    envelope2, errors2 = _make(["not", "an", "object"])
    assert envelope2 is None
    assert any("schema" in e for e in errors2), errors2

    envelope3, verdict = _preflight(
        {"schema": "wrong"}, region_ids=REGION_IDS, mode="synthetic_demo"
    )
    assert envelope3 is None
    assert verdict.state == "REFUSE"
    assert any("schema" in f for f in verdict.failures), verdict.failures


# ---------------------------------------------------------------- strict vs synthetic demo --
def test_synthetic_example_refuses_strict_and_accepts_synthetic_demo():
    payload = _load_example()

    strict_env, strict = _preflight(payload, region_ids=REGION_IDS, mode="strict")
    assert strict.state == "REFUSE"
    assert any("disallowed provenance" in f for f in strict.failures), strict.failures
    assert any("synthetic" in f for f in strict.failures), strict.failures
    # the envelope is returned even on a post-wrap REFUSE, so a caller can inspect it.
    assert strict_env is not None

    demo_env, demo = _preflight(
        payload, region_ids=REGION_IDS, mode="synthetic_demo"
    )
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert demo.accepted is True
    assert demo.failures == ()
    assert demo_env is not None


def test_all_calibrated_copy_accepts_scientific():
    envelope, verdict = _preflight(
        _all_calibrated_copy(), region_ids=REGION_IDS, mode="strict"
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.accepted is True
    assert verdict.evidence_class == "calibrated"
    assert envelope is not None
    assert envelope["producer"]["provenance"] == "calibrated"


# --------------------------------------------------------------------------- region gates ---
def test_unregistered_region_id_is_refused_as_unknown_region_id():
    payload = _load_example()
    payload["regions"][0]["region_id"] = "MSK-LOWERLIMB-HIP"
    _env, verdict = _preflight(
        payload, region_ids=REGION_IDS, mode="synthetic_demo"
    )
    assert verdict.state == "REFUSE"
    assert any(
        "unknown region_id" in f and "MSK-LOWERLIMB-HIP" in f for f in verdict.failures
    ), verdict.failures


def test_requested_subset_refuses_the_extra_region():
    payload = _load_example()
    _env, verdict = _preflight(
        payload, region_ids=("SYNTH-BOX-1",), mode="synthetic_demo"
    )
    assert verdict.state == "REFUSE"
    assert any(
        "unknown region_id" in f and "SYNTH-BOX-2" in f for f in verdict.failures
    ), verdict.failures


def test_msk_lowerlimb_hip_refused_as_unknown_region_id_by_the_validator():
    """Validator unit test: the new I2-sourced vocabulary refuses the old placeholder."""
    assert "MSK-LOWERLIMB-HIP" not in cp.THERMO_REGION_VOCABULARY
    region_id = "MSK-LOWERLIMB-HIP"
    key = f"regions.{region_id}.mass"
    envelope = base.make_envelope(
        {"regions": {region_id: {"mass": {"value": 1.0}}}},
        producer_name="x",
        producer_sha256=None,
        provenance="synthetic",
        input_sha256=None,
        quantity_specs={
            key: {
                "units": "kg",
                "region_id": region_id,
                "frame": "scalar",
                "time": None,
            }
        },
    )
    verdict = base.validate_for_consumer(
        envelope,
        expected_quantities={
            key: {
                "units": "kg",
                "region_id": region_id,
                "frame": "scalar",
                "time": None,
                "required": True,
            }
        },
        allowed_provenance=base.allowed_provenance_for(opt_in=("synthetic",)),
        mode="synthetic_demo",
        expected_regions=cp.THERMO_REGION_VOCABULARY,
    )
    assert verdict.state == "REFUSE"
    assert any(
        "unknown region_id" in f and region_id in f for f in verdict.failures
    ), verdict.failures


# ------------------------------------------------------------------- real I2 validator ------
def test_i2_structural_validator_accepts_example_and_rejects_a_bare_number():
    validate_payload = _load_i2_validator()
    payload = _load_example()
    # The frozen validator hash is asserted by the dedicated frozen-source test above; here we
    # require only that the real structural validator (whatever revision) accepts the example and
    # rejects a bare number, so the functional seam stays covered while the validator evolves.
    assert validate_payload(payload) == []

    bad = copy.deepcopy(payload)
    bad["regions"][0]["mass"] = 0.03816  # bare number where a quantity object is required
    errors = validate_payload(bad)
    assert errors
    assert any(BARE_NUMBER_SUBSTRING in e for e in errors), errors

    envelope, wrapper_errors = _make(bad, structural_validator=validate_payload)
    assert envelope is None
    assert wrapper_errors == errors

    envelope2, verdict = _preflight(
        bad, region_ids=REGION_IDS, mode="synthetic_demo", structural_validator=validate_payload
    )
    assert envelope2 is None
    assert verdict.state == "REFUSE"
    assert any(BARE_NUMBER_SUBSTRING in f for f in verdict.failures), verdict.failures


# --------------------------------------------------------------------- pure unit helpers ----
def test_derive_provenance_returns_the_least_trustworthy_term():
    assert rev2.derive_provenance([]) == "unknown"
    assert rev2.derive_provenance([None]) == "unknown"
    assert rev2.derive_provenance(["not-a-term"]) == "unknown"
    assert rev2.derive_provenance(["measured"]) == "measured"
    assert rev2.derive_provenance(["calibrated", "literal_cited"]) == "literal_cited"
    assert rev2.derive_provenance(["literal_cited", "synthetic"]) == "synthetic"
    assert rev2.derive_provenance(["measured", "synthetic"]) == "synthetic"
    assert rev2.derive_provenance(["calibrated_public_model", "tuned"]) == "tuned"
    assert rev2.derive_provenance(["measured", "assumption"]) == "unknown"


def test_region_expected_quantities_shape():
    expected = rev2.region_expected_quantities(("A", "B"))
    assert set(expected) == {
        "regions.A.mass",
        "regions.A.heat_capacity",
        "regions.B.mass",
        "regions.B.heat_capacity",
        "totals.mass",
        "totals.heat_capacity",
    }
    assert expected["regions.A.mass"] == {
        "units": "kg",
        "region_id": "A",
        "frame": "scalar",
        "time": None,
        "required": True,
    }
    assert expected["regions.B.heat_capacity"] == {
        "units": "J/K",
        "region_id": "B",
        "frame": "scalar",
        "time": None,
        "required": True,
    }
    assert expected["totals.mass"]["region_id"] is None
    assert expected["totals.heat_capacity"]["region_id"] is None

    without_totals = rev2.region_expected_quantities(("A",), require_totals=False)
    assert "totals.mass" not in without_totals
    assert "regions.A.mass" in without_totals


def test_computed_provenance_maps_to_unknown_and_is_not_promoted():
    assert base.normalize_provenance("computed") == "unknown"
    assert rev2.derive_provenance(["computed", "measured"]) == "unknown"

    payload = _load_example()
    payload["synthetic"] = False
    for region in payload["regions"]:
        region["mass"]["provenance"] = "computed"
        region["heat_capacity"]["provenance"] = "computed"
    payload["totals"]["mass"]["provenance"] = "computed"
    payload["totals"]["heat_capacity"]["provenance"] = "computed"

    _env, strict = _preflight(payload, region_ids=REGION_IDS, mode="strict")
    assert strict.state == "REFUSE"
    assert any("disallowed provenance" in f for f in strict.failures), strict.failures

    # It remains admissible only via the explicit synthetic/unknown opt-in demo path.
    demo_env, demo = _preflight(
        payload, region_ids=REGION_IDS, mode="synthetic_demo"
    )
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert demo_env["producer"]["provenance"] == "unknown"


# --------------------------------------------------------- abstention ------------
def test_region_abstention_is_refused_by_the_i1_envelope_though_i2_accepts_it():
    """The payload schema permits an explicit abstention (``value=null`` + non-empty ``missing``); the result
    envelope must refuse a required quantity whose value is None. The payload's own structural validator
    accepts the payload, so this pins the gap the result envelope closes (an unknown mass is not a
    mass).

    The copy is all-``calibrated``/``synthetic=False`` so the abstention is the ONLY failure in
    strict mode (provenance would otherwise dominate)."""
    payload = _all_calibrated_copy()
    payload["regions"][0]["mass"]["value"] = None
    payload["regions"][0]["mass"]["missing"] = ["mass"]

    validate_payload = _load_i2_validator()
    assert validate_payload(payload) == []  # the structural validator itself accepts the explicit abstention

    envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        mode="strict",
        structural_validator=validate_payload,
    )
    assert verdict.state == "REFUSE"
    assert any("abstain:" in f for f in verdict.failures), verdict.failures
    # the envelope is still returned on a post-wrap REFUSE so a caller can inspect it.
    assert envelope is not None
    assert envelope["quantities"]["regions.SYNTH-BOX-1.mass"]["value"] is None


def test_region_payload_valid_provenance_detail_is_accepted_when_required():
    """A region payload with a real region provenance as ``data`` passes the required detail."""
    payload = _all_calibrated_copy()
    detail = {
        # the frozen example's mass provenance is ``literal_cited`` -> normalises to calibrated.
        "data": "literal_cited",
        "calibration": "literal_cited",
        "model": "region_mass_v1",
        "execution": {"code_sha256": "a" * 64, "command": "python make_region.py"},
    }
    envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        mode="strict",
        provenance_detail=detail,
        require_provenance_detail=True,
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert envelope is not None
    assert envelope["producer"]["provenance_detail"] == detail


def test_region_payload_without_provenance_detail_is_refused_when_required():
    _envelope, verdict = _preflight(
        _all_calibrated_copy(),
        region_ids=REGION_IDS,
        mode="strict",
        require_provenance_detail=True,
    )
    assert verdict.state == "REFUSE"
    assert any("missing provenance_detail" in f for f in verdict.failures), verdict.failures


# ------------------------------------------------ uncertainty at the region-payload seam ----------
def _all_calibrated_with_uncertainty():
    """``_all_calibrated_copy()`` with an ``uncertainty`` field on every I2 quantity object."""
    payload = _all_calibrated_copy()
    for region in payload["regions"]:
        region["mass"]["uncertainty"] = 0.01
        region["heat_capacity"]["uncertainty"] = 0.02
    payload["totals"]["mass"]["uncertainty"] = 0.0  # 0.0 is a valid declared uncertainty
    payload["totals"]["heat_capacity"]["uncertainty"] = 0.04
    return payload


def test_region_uncertainty_is_threaded_verbatim_into_the_envelope():
    envelope, errors = _make(_all_calibrated_with_uncertainty())
    assert errors == []
    quantities = envelope["quantities"]
    assert quantities["regions.SYNTH-BOX-1.mass"]["uncertainty"] == 0.01
    assert quantities["regions.SYNTH-BOX-1.heat_capacity"]["uncertainty"] == 0.02
    assert quantities["totals.mass"]["uncertainty"] == 0.0
    # Absent on the original example -> None, never an invented value.
    plain, _ = _make(_load_example())
    assert plain["quantities"]["regions.SYNTH-BOX-1.mass"]["uncertainty"] is None


def test_region_payload_uncertainty_accepted_under_the_strict_profile():
    payload = _all_calibrated_with_uncertainty()
    _envelope, verdict = _preflight(
        payload, region_ids=REGION_IDS, mode="strict", require_uncertainty=True
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_region_payload_without_uncertainty_refuses_under_the_strict_profile():
    payload = _all_calibrated_copy()
    _envelope, verdict = _preflight(
        payload, region_ids=REGION_IDS, mode="strict", require_uncertainty=True
    )
    assert verdict.state == "REFUSE"
    assert any("missing uncertainty" in f for f in verdict.failures), verdict.failures


# ------------------------------------------ producer gate outcome at the seam ---------
def test_region_envelope_gate_defaults_keep_the_historical_envelope():
    """An I2 payload carries no gates by contract; the wrapped envelope keeps ``{}`` / ``True``."""
    envelope, errors = _make(_load_example())
    assert errors == []
    assert envelope["gates"] == {}
    assert envelope["overall_pass"] is True


def test_region_envelope_with_overall_pass_false_flags_or_refuses_per_the_parameter():
    """At the seam: a region-derived envelope can carry an upstream producer's real
    ``overall_pass=false``. Strict (scientific) mode refuses it; the explicit demo mode flags it
    unless ``require_producer_gates_pass=True``."""
    gates = {"region_mass_sanity": False, "regions_present": True}
    envelope, errors = _make(_all_calibrated_copy(), gates=gates, overall_pass=False)
    assert errors == []
    assert envelope["gates"] == gates
    assert envelope["overall_pass"] is False

    common = dict(
        expected_quantities=rev2.region_expected_quantities(REGION_IDS),
        allowed_provenance=base.allowed_provenance_for(),
        expected_regions=frozenset(REGION_IDS),
        expected_payload_sha256=base.payload_sha256(envelope),
    )
    strict = base.validate_for_consumer(envelope, mode="strict", **common)
    assert strict.state == "REFUSE"
    assert "producer gates failed: overall_pass is not true" in strict.failures

    demo = base.validate_for_consumer(envelope, mode="synthetic_demo", **common)
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert "producer gates failed: overall_pass is not true" in demo.flags

    required = base.validate_for_consumer(
        envelope, mode="synthetic_demo", require_producer_gates_pass=True, **common
    )
    assert required.state == "REFUSE"
    assert "producer gates failed: overall_pass is not true" in required.failures


def test_preflight_region_payload_threads_the_producer_gate_outcome():
    """``preflight_region_payload`` forwards ``gates``/``overall_pass``; strict mode refuses."""
    payload = _all_calibrated_copy()

    _env, strict = _preflight(
        payload, region_ids=REGION_IDS, mode="strict",
        gates={"region_mass_sanity": False}, overall_pass=False,
    )
    assert strict.state == "REFUSE"
    assert "producer gates failed: overall_pass is not true" in strict.failures
    with pytest.raises(ValueError, match="strict mode cannot waive producer gates"):
        _preflight(payload, region_ids=REGION_IDS, mode="strict", gates={"g": False},
                   overall_pass=False, require_producer_gates_pass=False)


def test_region_scientific_acceptance_requires_the_payload_pin():
    payload = _all_calibrated_copy()
    _env, unpinned = rev2.preflight_region_payload(payload, region_ids=REGION_IDS, mode="strict")
    assert unpinned.state == "REFUSE"
    assert ("no trusted artifact pin: expected_payload_sha256 is required for scientific "
            "acceptance") in unpinned.failures
    _env, other = rev2.preflight_region_payload(
        payload, region_ids=REGION_IDS, mode="strict", expected_payload_sha256="0" * 64)
    assert other.state == "REFUSE"
    assert any(f.startswith("artifact not approved: region payload sha256") for f in other.failures)
    _env, ok = _preflight(payload, region_ids=REGION_IDS, mode="strict")
    assert ok.state == "ACCEPT_SCIENTIFIC"


# ------------------------------------- area-kind tagging at the region-payload seam -----------------
# An area is not interchangeable with a different area kind. I2's ``area.kind`` is
# carried into the envelope key (``area.<kind>``); a consumer must require the exact key, so a
# wall (``exterior_area``) or a shared face (``material_interface_area``) can never silently stand
# in for a ``functional_exchange_area``.
def _calibrated_with_area_provenance(provenance="calibrated"):
    """``_all_calibrated_copy()`` with a declared ``area.provenance`` (I2's example has none)."""
    payload = _all_calibrated_copy()
    payload["area"]["provenance"] = provenance
    return payload


def test_area_kind_is_carried_as_a_kind_tagged_quantity():
    payload = _calibrated_with_area_provenance()
    envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        area_kinds=["material_interface_area"],
        mode="strict",
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert envelope is not None
    area = envelope["quantities"]["area.material_interface_area"]
    assert area["units"] == "mm2"
    assert area["value"] == pytest.approx(3600.0)
    assert area["region_id"] is None
    assert area["frame"] == "scalar"
    assert area["time"] is None
    assert area["provenance"] == "calibrated"


def test_the_wall_area_cannot_stand_in_for_a_functional_exchange_area():
    """The example's area is ``material_interface_area``; a consumer that needs an exchange area
    must get ``missing quantity: area.functional_exchange_area`` -- the key encodes the kind, so an
    area/volume (or wall/exchange) swap cannot be silently absorbed."""
    payload = _calibrated_with_area_provenance()
    _envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        area_kinds=["functional_exchange_area"],
        mode="strict",
    )
    assert verdict.state == "REFUSE"
    assert "missing quantity: area.functional_exchange_area" in verdict.failures, verdict.failures
    # The material-interface area is present but is a DIFFERENT quantity and must not satisfy it.
    assert not any(
        "area.material_interface_area" in f and "missing quantity" in f for f in verdict.failures
    )


def test_area_with_wrong_units_is_refused():
    """A length (``mm``) or a volume (``mm3``) must not pass as an area (``mm2``).

    The unit is pinned, so an area/volume swap cannot be silently absorbed by the envelope.
    """
    for wrong_units in ("mm", "mm3"):
        payload = _calibrated_with_area_provenance()
        payload["area"]["units"] = wrong_units
        _envelope, verdict = _preflight(
            payload,
            region_ids=REGION_IDS,
            area_kinds=["material_interface_area"],
            mode="strict",
        )
        assert verdict.state == "REFUSE", (wrong_units, verdict.failures)
        assert any(
            "unit mismatch" in f and "area.material_interface_area" in f
            for f in verdict.failures
        ), (wrong_units, verdict.failures)


def test_unprovenanced_area_is_refused_when_expected_but_ignored_otherwise():
    """I2's frozen example declares no ``area.provenance``.

    A consumer that explicitly expects the area quantity must not accept an unprovenanced area
    (``disallowed provenance`` on the per-quantity ``None`` -> ``unknown``); a consumer that does
    NOT expect area is unchanged and the undeclared area provenance does not drag the whole
    envelope's top-level provenance to ``unknown``.
    """
    payload = _all_calibrated_copy()
    assert "provenance" not in payload["area"]  # the frozen I2 example indeed declares none

    # Explicitly expected -> refused as an unprovenanced scientific quantity.
    _envelope, expected = _preflight(
        payload,
        region_ids=REGION_IDS,
        area_kinds=["material_interface_area"],
        mode="strict",
    )
    assert expected.state == "REFUSE"
    assert any("disallowed provenance" in f for f in expected.failures), expected.failures
    assert all("area.material_interface_area" in f for f in expected.failures if "provenance" in f)

    # Not expected -> the existing calibrated-copy acceptance is unchanged (no top-level regression).
    envelope, verdict = _preflight(
        payload, region_ids=REGION_IDS, mode="strict"
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert envelope is not None
    assert envelope["producer"]["provenance"] == "calibrated"
    assert envelope["quantities"]["area.material_interface_area"]["provenance"] == "unknown"


def test_area_expected_quantities_shape():
    assert rev2.area_expected_quantities(("material_interface_area",)) == {
        "area.material_interface_area": {
            "units": "mm2",
            "region_id": None,
            "frame": "scalar",
            "time": None,
            "required": True,
        }
    }
    optional = rev2.area_expected_quantities(("exterior_area",), required=False)
    assert optional == {
        "area.exterior_area": {
            "units": "mm2",
            "region_id": None,
            "frame": "scalar",
            "time": None,
            "required": False,
        }
    }
    assert rev2.area_expected_quantities([]) == {}


# ----------------------------------- material-region binding at the seam ------------
# The region payload carries a per-region ``material_region_id`` (an explicit region -> material
# binding, e.g. SYNTH-BOX-1 -> MSK-MUSCLE). The adapter copies it into BOTH region quantity specs
# and refuses an unregistered id outright; a consumer enforces an exact binding via
# ``preflight_region_payload(material_region_map=...)``.
MATERIAL_REGION_MAP = {"SYNTH-BOX-1": "MSK-MUSCLE", "SYNTH-BOX-2": "MSK-MUSCLE"}


def test_region_material_region_ids_are_copied_to_both_quantities_and_are_registered():
    payload = _load_example()
    # the frozen example's bindings are all members of the I2 material roster.
    for region in payload["regions"]:
        assert region["material_region_id"] in rvf.I2_MATERIAL_REGION_IDS

    envelope, errors = _make(payload)
    assert errors == []
    quantities = envelope["quantities"]
    for region_id, material_region_id in (
        ("SYNTH-BOX-1", "MSK-MUSCLE"),
        ("SYNTH-BOX-2", "MSK-MUSCLE"),
    ):
        assert quantities[f"regions.{region_id}.mass"]["material_region_id"] == material_region_id
        assert (
            quantities[f"regions.{region_id}.heat_capacity"]["material_region_id"]
            == material_region_id
        )
    # no binding is invented on the aggregates.
    assert "material_region_id" not in quantities["totals.mass"]
    assert "material_region_id" not in quantities["totals.heat_capacity"]


def test_preflight_region_payload_accepts_the_correct_material_region_map():
    payload = _all_calibrated_copy()
    envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        mode="strict",
        material_region_map=MATERIAL_REGION_MAP,
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert envelope is not None
    assert (
        envelope["quantities"]["regions.SYNTH-BOX-1.mass"]["material_region_id"] == "MSK-MUSCLE"
    )


def test_preflight_region_payload_refuses_a_swapped_material_region_map():
    payload = _all_calibrated_copy()
    swapped = {"SYNTH-BOX-1": "MSK-BONE", "SYNTH-BOX-2": "MSK-MUSCLE"}
    _envelope, verdict = _preflight(
        payload, region_ids=REGION_IDS, mode="strict", material_region_map=swapped
    )
    assert verdict.state == "REFUSE"
    assert any(
        "material region mismatch: regions.SYNTH-BOX-1.mass" in f
        and "'MSK-BONE'" in f
        and "'MSK-MUSCLE'" in f
        for f in verdict.failures
    ), verdict.failures
    assert any(
        "material region mismatch: regions.SYNTH-BOX-1.heat_capacity" in f
        for f in verdict.failures
    ), verdict.failures
    # the correctly-bound region is untouched.
    assert not any("SYNTH-BOX-2" in f and "material region" in f for f in verdict.failures)


def test_unregistered_material_region_id_is_refused_by_the_adapter():
    """The adapter (``make_region_envelope``), not the consumer, catches an unregistered id: it
    refuses to wrap the payload at all (envelope is ``None``) with a stable ``material region``
    substring, independent of any consumer contract."""
    payload = _all_calibrated_copy()
    payload["regions"][0]["material_region_id"] = "MSK-NOPE"

    envelope, errors = _make(payload)
    assert envelope is None
    assert any("material region" in e for e in errors), errors

    env2, verdict = _preflight(payload, region_ids=REGION_IDS, mode="strict")
    assert env2 is None
    assert verdict.state == "REFUSE"
    assert any("material region" in f for f in verdict.failures), verdict.failures


def test_unknown_material_region_vocabulary_is_enforced_by_the_shared_validator():
    """The validator-level ``unknown material region`` rule (the second layer): an envelope built
    directly (bypassing the adapter) with an unregistered binding is refused when the consumer
    passes the material vocabulary."""
    key = "regions.SYNTH-BOX-1.mass"
    envelope = base.make_envelope(
        {"regions": {"SYNTH-BOX-1": {"mass": {"value": 1.0}}}},
        producer_name="x",
        producer_sha256=None,
        provenance="synthetic",
        input_sha256=None,
        quantity_specs={
            key: {
                "units": "kg",
                "region_id": "SYNTH-BOX-1",
                "frame": "scalar",
                "time": None,
                "material_region_id": "MSK-NOPE",
            }
        },
    )
    verdict = base.validate_for_consumer(
        envelope,
        expected_quantities={
            key: {
                "units": "kg",
                "region_id": "SYNTH-BOX-1",
                "frame": "scalar",
                "time": None,
                "required": True,
            }
        },
        allowed_provenance=base.allowed_provenance_for(opt_in=("synthetic",)),
        mode="synthetic_demo",
        expected_material_regions=frozenset(rvf.I2_MATERIAL_REGION_IDS),
        material_region_vocabulary=frozenset(rvf.I2_MATERIAL_REGION_IDS),
    )
    assert verdict.state == "REFUSE"
    assert any(
        "unknown material region: regions.SYNTH-BOX-1.mass 'MSK-NOPE'" in f
        for f in verdict.failures
    ), verdict.failures


def test_no_material_region_map_does_not_enforce_the_binding_but_the_roster_applies():
    payload = _all_calibrated_copy()

    # no map supplied -> the registered binding is accepted (not enforced).
    _env, accepted = _preflight(payload, region_ids=REGION_IDS, mode="strict")
    assert accepted.state == "ACCEPT_SCIENTIFIC"
    assert accepted.failures == ()

    # a *different but registered* binding is also accepted when no map is supplied.
    swapped = copy.deepcopy(payload)
    swapped["regions"][0]["material_region_id"] = "MSK-BONE"
    _env2, swapped_verdict = _preflight(
        swapped, region_ids=REGION_IDS, mode="strict"
    )
    assert swapped_verdict.state == "ACCEPT_SCIENTIFIC"
    assert swapped_verdict.failures == ()

    # the roster vocabulary still applies: an unregistered id is refused by the adapter.
    bad = copy.deepcopy(payload)
    bad["regions"][0]["material_region_id"] = "MSK-NOPE"
    bad_env, errors = _make(bad)
    assert bad_env is None
    assert any("material region" in e for e in errors), errors


def test_region_expected_quantities_material_region_map_shape():
    expected = rev2.region_expected_quantities(
        ("A", "B"), material_region_map={"A": "MSK-BONE"}
    )
    assert expected["regions.A.mass"]["material_region_id"] == "MSK-BONE"
    assert expected["regions.A.heat_capacity"]["material_region_id"] == "MSK-BONE"
    assert "material_region_id" not in expected["regions.B.mass"]
    assert "material_region_id" not in expected["regions.B.heat_capacity"]
    assert "material_region_id" not in expected["totals.mass"]

    # no map -> the historical expected-quantity shape is unchanged (no extra key).
    without_map = rev2.region_expected_quantities(("A",))
    assert "material_region_id" not in without_map["regions.A.mass"]

    # a map entry for an id the caller did not request is a clear caller error.
    with pytest.raises(ValueError, match="material region map"):
        rev2.region_expected_quantities(("A",), material_region_map={"TYPO": "MSK-BONE"})


def test_preflight_rejects_a_material_region_map_with_an_unknown_local_id():
    _envelope, verdict = _preflight(
        _all_calibrated_copy(),
        region_ids=REGION_IDS,
        mode="strict",
        material_region_map={"TYPO": "MSK-BONE"},
    )
    assert verdict.state == "REFUSE"
    assert any("material region map" in f for f in verdict.failures), verdict.failures


# ----------------------------------- per-region / area frame at the seam ----------
# The payload may declare an optional per-region ``frame`` (and an ``area.frame``): a frame descriptor
# dict (here ``{origin_mm, axes, pitch_mm, units, declared}``) or a string. The adapter consumes it
# verbatim; absent/``None`` keeps ``scalar``. The frozen example declares no frame at all, so the
# default contract is unchanged and an undeclared frame is never invented.
FRAME_DICT = {
    "origin_mm": [10.0, 20.0, 30.0],
    "axes": [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]],
    "pitch_mm": 2.0,
    "units": "mm",
    "declared": True,
}


def test_region_frame_is_consumed_when_declared():
    payload = _all_calibrated_copy()
    payload["regions"][0]["frame"] = FRAME_DICT
    # The framed copy is accepted by the REAL structural validator (exercised when present).
    envelope, errors = _make(payload, structural_validator=_load_i2_validator())
    assert errors == []
    assert envelope is not None
    for name in ("mass", "heat_capacity"):
        assert envelope["quantities"][f"regions.SYNTH-BOX-1.{name}"]["frame"] == FRAME_DICT
    # a region with no declared frame and the aggregates keep the scalar convention.
    assert envelope["quantities"]["regions.SYNTH-BOX-2.mass"]["frame"] == "scalar"
    assert envelope["quantities"]["regions.SYNTH-BOX-2.heat_capacity"]["frame"] == "scalar"
    assert envelope["quantities"]["totals.mass"]["frame"] == "scalar"


def test_consumer_expecting_scalar_frame_refuses_a_declared_region_frame():
    payload = _all_calibrated_copy()
    payload["regions"][0]["frame"] = FRAME_DICT
    _envelope, verdict = _preflight(
        payload, region_ids=REGION_IDS, mode="strict"
    )
    assert verdict.state == "REFUSE"
    assert any(
        "frame mismatch: regions.SYNTH-BOX-1.mass" in f for f in verdict.failures
    ), verdict.failures
    assert any(
        "frame mismatch: regions.SYNTH-BOX-1.heat_capacity" in f for f in verdict.failures
    ), verdict.failures
    # the region without a declared frame is untouched.
    assert not any("SYNTH-BOX-2" in f and "frame mismatch" in f for f in verdict.failures)


def test_consumer_expecting_the_declared_region_frame_accepts_canonically():
    payload = _all_calibrated_copy()
    payload["regions"][0]["frame"] = FRAME_DICT
    # A different dict key ORDER must still match: the shared validator compares sorted JSON.
    reordered = {key: FRAME_DICT[key] for key in reversed(list(FRAME_DICT))}
    assert list(reordered) != list(FRAME_DICT)
    envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        mode="strict",
        frames={"SYNTH-BOX-1": reordered},
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert envelope["quantities"]["regions.SYNTH-BOX-1.mass"]["frame"] == FRAME_DICT


def test_region_frame_null_or_absent_keeps_scalar():
    # the frozen example declares no per-region frame at all ...
    payload = _load_example()
    for region in payload["regions"]:
        assert "frame" not in region
    envelope, errors = _make(payload)
    assert errors == []
    assert envelope["quantities"]["regions.SYNTH-BOX-1.mass"]["frame"] == "scalar"

    # ... and an explicit ``frame: null`` keeps the scalar convention too (never invented).
    payload_null = _load_example()
    for region in payload_null["regions"]:
        region["frame"] = None
    envelope_null, errors_null = _make(payload_null)
    assert errors_null == []
    for name in ("mass", "heat_capacity"):
        assert (
            envelope_null["quantities"][f"regions.SYNTH-BOX-2.{name}"]["frame"] == "scalar"
        )


def test_area_frame_is_consumed_when_declared():
    payload = _calibrated_with_area_provenance()
    payload["area"]["frame"] = FRAME_DICT

    # a scalar expectation refuses the declared area frame ...
    _envelope, refused = _preflight(
        payload,
        region_ids=REGION_IDS,
        area_kinds=["material_interface_area"],
        mode="strict",
    )
    assert refused.state == "REFUSE"
    assert any(
        "frame mismatch: area.material_interface_area" in f for f in refused.failures
    ), refused.failures

    # ... requiring the declared area frame accepts ...
    envelope, verdict = _preflight(
        payload,
        region_ids=REGION_IDS,
        area_kinds=["material_interface_area"],
        mode="strict",
        frames={"material_interface_area": FRAME_DICT},
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert envelope["quantities"]["area.material_interface_area"]["frame"] == FRAME_DICT

    # ... and an explicit null area frame keeps the scalar convention.
    payload_null = _calibrated_with_area_provenance()
    payload_null["area"]["frame"] = None
    envelope_null, errors = _make(payload_null)
    assert errors == []
    assert envelope_null["quantities"]["area.material_interface_area"]["frame"] == "scalar"


def test_region_and_area_expected_frames_shape():
    frame = {"origin_mm": [0.0, 0.0, 0.0]}
    expected = rev2.region_expected_quantities(("A", "B"), frames={"A": frame})
    assert expected["regions.A.mass"]["frame"] == frame
    assert expected["regions.A.heat_capacity"]["frame"] == frame
    assert expected["regions.B.mass"]["frame"] == "scalar"
    assert expected["totals.mass"]["frame"] == "scalar"

    # frames keyed by a region id do not leak into the area contract and vice versa.
    area_expected = rev2.area_expected_quantities(
        ("material_interface_area",), frames={"A": frame}
    )
    assert area_expected["area.material_interface_area"]["frame"] == "scalar"
    framed_area = rev2.area_expected_quantities(
        ("material_interface_area",), frames={"material_interface_area": frame}
    )
    assert framed_area["area.material_interface_area"]["frame"] == frame

    # default shapes are unchanged (no ``frames`` -> scalar).
    assert rev2.region_expected_quantities(("A",))["regions.A.mass"]["frame"] == "scalar"
    assert rev2.region_expected_quantities(("A",), frames={"A": None})["regions.A.mass"][
        "frame"
    ] is None
    assert rev2.area_expected_quantities(()) == {}

    # the region call can build the whole region+area contract in one go (area_kinds convenience).
    combined = rev2.region_expected_quantities(
        ("A",),
        area_kinds=["material_interface_area"],
        frames={"material_interface_area": frame},
    )
    assert combined["area.material_interface_area"]["frame"] == frame
    assert combined["regions.A.mass"]["frame"] == "scalar"
