"""Acceptance tests for I1 ``bodytwin.framework.result_envelope_v1``.

Written to the exact public API + stable failure substrings of the result envelope.
Stdlib + pytest only.  The module is imported from ``src`` (``PYTHONPATH=src``); an import
error fails the test run.

Coverage:
  * ``make_envelope`` shape / defaults / only-present-quantities
  * every ``validate_for_consumer`` failure class (unit, region, unknown region, frame,
    time, absent provenance, changed producer, name mismatch, stale input)
  * strict vs opt-in provenance semantics
  * receipt seal/verify round-trip + tamper + producer hash mismatch
  * no universal mass-fraction bound (0.585 accepted/flagged), physical impossibility refused
  * ``allowed_provenance_for`` opt-in semantics
"""
import pytest

import bodytwin.framework.result_envelope_v1 as rev  # noqa: F401

PRODUCER_SHA = "6bf7599c63a2203b175c43749a20d4a52ccf99fabefb8303c15009243142fe38"
INPUT_SHA = "a" * 64
SCIENTIFIC = {"measured", "calibrated", "calibrated_public_model"}


# --------------------------------------------------------------------------- helpers -------
def _result_obj():
    return {
        "muscle_mass": {"total_body_mass_kg": 72.6, "summed_muscle_mass_kg": 24.0},
        "headline": {"gross_w_per_kg": 5.0, "net_w_per_kg": 4.5},
    }


def _specs():
    return {
        "muscle_mass.total_body_mass_kg": {
            "units": "kg", "region_id": "WHOLE-BODY", "frame": "scalar",
            "time": None, "uncertainty": 0.4,
            "provenance": "calibrated_public_model",
        },
        "muscle_mass.summed_muscle_mass_kg": {
            "units": "kg", "region_id": "MSK-MUSCLE", "frame": "scalar",
            "time": None, "uncertainty": 0.4,
            "provenance": "calibrated_public_model",
        },
        "headline.gross_w_per_kg": {
            "units": "W/kg", "region_id": "WHOLE-BODY", "frame": "scalar",
            "time": "steady_state",
        },
        "headline.net_w_per_kg": {
            "units": "W/kg", "region_id": "WHOLE-BODY", "frame": "scalar",
            "time": "steady_state",
        },
    }


def _expected():
    return {
        "muscle_mass.total_body_mass_kg": {
            "units": "kg", "region_id": "WHOLE-BODY", "frame": "scalar",
            "time": None, "required": True,
        },
        "muscle_mass.summed_muscle_mass_kg": {
            "units": "kg", "region_id": "MSK-MUSCLE", "frame": "scalar",
            "time": None, "required": True,
        },
        "headline.gross_w_per_kg": {
            "units": "W/kg", "region_id": "WHOLE-BODY", "frame": "scalar",
            "time": "steady_state", "required": True,
        },
        "headline.net_w_per_kg": {
            "units": "W/kg", "region_id": "WHOLE-BODY", "frame": "scalar",
            "time": "steady_state", "required": True,
        },
    }


def _envelope(**overrides):
    args = dict(
        result_obj=_result_obj(),
        producer_name="metabolic_cost",
        producer_sha256=PRODUCER_SHA,
        provenance="calibrated_public_model",
        input_sha256=INPUT_SHA,
        quantity_specs=_specs(),
        command="python metabolic_cost_adapted.py",
        seed=7,
        regime="walking",
    )
    args.update(overrides)
    return rev.make_envelope(**args)


def _validate(env=None, **overrides):
    """Validate a valid envelope, with the stale-input check satisfied by default."""
    if env is None:
        env = _envelope()
    args = dict(
        expected_quantities=_expected(),
        allowed_provenance=rev.allowed_provenance_for(),
        current_input_sha256=INPUT_SHA,
        # the trusted registry pin of exactly this envelope (scientific acceptance requires it)
        expected_payload_sha256=rev.payload_sha256(env),
    )
    args.update(overrides)
    return rev.validate_for_consumer(env, **args)


def _has(failures, needle):
    return any(needle in f for f in failures), tuple(failures)


# --------------------------------------------------------------------------- constants -----
def test_constants_match_spec():
    assert rev.ENVELOPE_VERSION == 1
    assert rev.SCALAR_FRAME == "scalar"
    assert set(rev.PROVENANCE_VOCABULARY) == {
        "measured", "calibrated", "calibrated_public_model", "literal_cited", "tuned",
        "synthetic", "unknown",
    }
    assert rev.DEFAULT_SCIENTIFIC_ALLOWED == frozenset(
        {"measured", "calibrated", "calibrated_public_model"}
    )
    assert rev.OPT_IN_REQUIRED == frozenset({"literal_cited", "tuned", "synthetic", "unknown"})
    assert set(rev.REGIME_VOCABULARY) == {
        "resting", "level_walking_steady_state", "load_carrying_steady_state",
        "synthetic_demo", "unknown",
    }
    for key in (
        "resting", "rest", "basal", "level_walking_steady_state",
        "level-walking-steady-state", "level_walking", "steady_state_walking",
        "load_carrying_steady_state", "load_carrying", "synthetic_demo", "synthetic", "unknown",
    ):
        assert key in rev.REGIME_ALIASES, key


def test_normalize_provenance_aliases():
    assert rev.normalize_provenance("calibrated-public-model") == "calibrated_public_model"
    assert rev.normalize_provenance("calibrated_from_public_model") == "calibrated_public_model"
    assert rev.normalize_provenance("public_model") == "calibrated_public_model"
    # a repository literal is its own category, never promoted to calibrated
    assert rev.normalize_provenance("literal_cited") == "literal_cited"
    assert rev.normalize_provenance("assumption") == "unknown"
    assert rev.normalize_provenance("UNKNOWN") == "unknown"
    assert rev.normalize_provenance(" calibrated ") == "calibrated"
    assert rev.normalize_provenance(None) == "unknown"
    assert rev.normalize_provenance(5) == "unknown"
    assert rev.normalize_provenance("not-a-term") == "unknown"
    for term in rev.PROVENANCE_VOCABULARY:
        assert rev.normalize_provenance(term) == term


def test_normalize_regime_aliases():
    assert rev.normalize_regime("rest") == "resting"
    assert rev.normalize_regime(" basal ") == "resting"
    assert rev.normalize_regime("RESTING") == "resting"
    assert rev.normalize_regime("level-walking-steady-state") == "level_walking_steady_state"
    assert rev.normalize_regime("level_walking") == "level_walking_steady_state"
    assert rev.normalize_regime("steady_state_walking") == "level_walking_steady_state"
    assert rev.normalize_regime("load_carrying") == "load_carrying_steady_state"
    assert rev.normalize_regime("synthetic") == "synthetic_demo"
    assert rev.normalize_regime(" unknown ") == "unknown"

    # absent / non-string / unrecognised nonsense -> "unknown" sentinel.
    assert rev.normalize_regime(None) == "unknown"
    assert rev.normalize_regime(5) == "unknown"
    assert rev.normalize_regime("") == "unknown"
    assert rev.normalize_regime("banana_regime") == "unknown"

    # a recognised but out-of-domain exercise regime canonicalises to itself (known, not unknown).
    assert rev.normalize_regime("sprint_anaerobic") == "sprint_anaerobic"

    for term in rev.REGIME_VOCABULARY:
        assert rev.normalize_regime(term) == term


def test_allowed_provenance_for_opt_in_semantics():
    assert rev.allowed_provenance_for() == rev.DEFAULT_SCIENTIFIC_ALLOWED
    assert "synthetic" not in rev.allowed_provenance_for()
    assert "tuned" not in rev.allowed_provenance_for()
    assert "unknown" not in rev.allowed_provenance_for()

    opted = rev.allowed_provenance_for(opt_in=("synthetic", "unknown"))
    assert {"synthetic", "unknown"} <= opted
    assert rev.DEFAULT_SCIENTIFIC_ALLOWED <= opted

    # opt-in accepts aliases and canonicalises them.
    assert "tuned" in rev.allowed_provenance_for(opt_in=("tuned",))

    # scientific=False is the tooling vocabulary: every declared non-scientific term.
    tooling = rev.allowed_provenance_for(scientific=False)
    assert rev.OPT_IN_REQUIRED <= tooling
    assert rev.DEFAULT_SCIENTIFIC_ALLOWED <= tooling


def test_get_dotted_and_missing_raises():
    assert rev.get_dotted(_result_obj(), "muscle_mass.total_body_mass_kg") == 72.6
    with pytest.raises(KeyError):
        rev.get_dotted(_result_obj(), "muscle_mass.not_a_key")
    with pytest.raises(KeyError):
        rev.get_dotted(_result_obj(), "no.such.path")


# --------------------------------------------------------------------------- envelope ------
def test_make_envelope_top_level_shape():
    env = _envelope()
    assert set(env) == {
        "envelope_version", "producer", "quantities", "validity", "gates", "overall_pass",
    }
    assert env["envelope_version"] == rev.ENVELOPE_VERSION
    assert set(env["producer"]) == {
        "name", "source_sha256", "command", "seed", "provenance", "input_sha256",
    }
    assert env["producer"]["name"] == "metabolic_cost"
    assert env["producer"]["source_sha256"] == PRODUCER_SHA
    assert env["producer"]["input_sha256"] == INPUT_SHA
    assert set(env["validity"]) == {
        "regime", "valid_from_utc", "valid_until_utc", "stale_if_input_sha256_mismatch",
    }
    assert env["validity"]["stale_if_input_sha256_mismatch"] is True


def test_make_envelope_quantity_shape_and_defaults():
    env = _envelope()
    assert set(env["quantities"]) == set(_expected())
    q = env["quantities"]["headline.gross_w_per_kg"]
    assert set(q) == {
        "value", "units", "region_id", "frame", "time", "uncertainty", "provenance",
    }
    assert q["value"] == 5.0
    assert q["units"] == "W/kg"
    assert q["frame"] == "scalar"
    assert q["time"] == "steady_state"
    # per-quantity provenance written where declared
    assert env["quantities"]["muscle_mass.total_body_mass_kg"]["provenance"] == (
        "calibrated_public_model"
    )
    # missing provenance -> producer provenance
    assert env["quantities"]["headline.gross_w_per_kg"]["provenance"] == (
        "calibrated_public_model"
    )


def test_make_envelope_missing_frame_time_region_defaults_and_present_only():
    specs = {
        "muscle_mass.total_body_mass_kg": {"units": "kg"},  # no frame/time/region
        "headline.gross_w_per_kg": {"units": "W/kg", "time": "steady_state"},
        "not.in_result": {"units": "x", "region_id": "WHOLE-BODY"},
    }
    env = _envelope(quantity_specs=specs)
    assert set(env["quantities"]) == {
        "muscle_mass.total_body_mass_kg", "headline.gross_w_per_kg",
    }
    total = env["quantities"]["muscle_mass.total_body_mass_kg"]
    assert total["frame"] == rev.SCALAR_FRAME
    assert total["time"] is None
    assert total["region_id"] is None


# --------------------------------------------------------------------------- validation ----
def test_valid_envelope_accepts_as_scientific():
    verdict = _validate()
    assert isinstance(verdict, rev.Verdict)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.accepted is True
    assert verdict.failures == ()
    assert verdict.evidence_class == "calibrated_public_model"


def test_missing_required_quantity_is_refused():
    obj = _result_obj()
    del obj["headline"]["net_w_per_kg"]
    env = _envelope(result_obj=obj)
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing quantity: headline.net_w_per_kg")
    assert ok, failures


def test_non_required_missing_quantity_is_not_a_failure():
    obj = _result_obj()
    del obj["headline"]["net_w_per_kg"]
    env = _envelope(result_obj=obj)
    expected = _expected()
    expected["headline.net_w_per_kg"]["required"] = False
    verdict = _validate(env, expected_quantities=expected)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert not any("missing quantity" in f for f in verdict.failures)


def test_unit_mismatch_is_refused():
    env = _envelope()
    env["quantities"]["muscle_mass.total_body_mass_kg"]["units"] = "g"
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "unit mismatch: muscle_mass.total_body_mass_kg")
    assert ok, failures


def test_region_mismatch_is_refused():
    env = _envelope()
    env["quantities"]["muscle_mass.summed_muscle_mass_kg"]["region_id"] = "MSK-LOWERLIMB-HIP"
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "region mismatch: muscle_mass.summed_muscle_mass_kg")
    assert ok, failures


def test_missing_region_id_is_refused():
    env = _envelope()
    env["quantities"]["muscle_mass.total_body_mass_kg"]["region_id"] = None
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing region_id: muscle_mass.total_body_mass_kg")
    assert ok, failures


def test_unknown_region_id_is_refused():
    # expected region_id None (so the mismatch branch is skipped) + a vocabulary gate.
    env = _envelope()
    env["quantities"]["muscle_mass.total_body_mass_kg"]["region_id"] = "BOGUS-REGION"
    expected = _expected()
    expected["muscle_mass.total_body_mass_kg"]["region_id"] = None
    verdict = _validate(
        env, expected_quantities=expected,
        expected_regions=frozenset({"WHOLE-BODY", "MSK-MUSCLE"}),
    )
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "unknown region_id: muscle_mass.total_body_mass_kg")
    assert ok, failures


def test_frame_mismatch_is_refused():
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["frame"] = "sagittal"
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "frame mismatch: headline.gross_w_per_kg")
    assert ok, failures


def test_frame_dict_is_canonicalised_by_sorted_json():
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["frame"] = {"b": 2, "a": 1}
    expected = _expected()
    expected["headline.gross_w_per_kg"]["frame"] = {"a": 1, "b": 2}
    verdict = _validate(env, expected_quantities=expected)
    assert not any("frame mismatch" in f for f in verdict.failures), verdict.failures


def test_time_mismatch_including_required_rate_with_none_is_refused():
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["time"] = None
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "time mismatch: headline.gross_w_per_kg")
    assert ok, failures


# ------------------------------- absent vs explicit expected frame/time ------------
# Previously ``spec.get("frame")``/``spec.get("time")`` made an OMITTED key indistinguishable
# from an explicit ``None``: a consumer that simply never declared a frame was demanding a null
# frame. Omitting the key now means "do not check"; ``"frame": None`` (or ``"time": None``) still
# means "expect null"; any other explicit value is still compared canonically.
def test_absent_expected_frame_and_time_keys_are_not_checked():
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["frame"] = "sagittal"
    env["quantities"]["headline.gross_w_per_kg"]["time"] = "instantaneous"
    expected = _expected()
    del expected["headline.gross_w_per_kg"]["frame"]
    del expected["headline.gross_w_per_kg"]["time"]
    verdict = _validate(env, expected_quantities=expected)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert not any("frame mismatch" in f for f in verdict.failures), verdict.failures
    assert not any("time mismatch" in f for f in verdict.failures), verdict.failures


def test_explicit_none_frame_and_time_still_expect_null():
    # An EXPLICIT null expectation is satisfied by a null quantity ...
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["frame"] = None
    env["quantities"]["headline.gross_w_per_kg"]["time"] = None
    expected = _expected()
    expected["headline.gross_w_per_kg"]["frame"] = None
    expected["headline.gross_w_per_kg"]["time"] = None
    verdict = _validate(env, expected_quantities=expected)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()

    # ... and a non-null quantity against an explicit-null expectation still mismatches.
    concrete = _envelope()  # frame "scalar", time "steady_state"
    verdict2 = _validate(concrete, expected_quantities=expected)
    assert verdict2.state == "REFUSE"
    ok, failures = _has(verdict2.failures, "frame mismatch: headline.gross_w_per_kg")
    assert ok, failures
    ok, failures = _has(verdict2.failures, "time mismatch: headline.gross_w_per_kg")
    assert ok, failures


def test_explicit_expected_frame_and_time_values_still_mismatch():
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["frame"] = "sagittal"
    env["quantities"]["headline.gross_w_per_kg"]["time"] = "instantaneous"
    verdict = _validate(env)  # _expected() declares scalar / steady_state explicitly
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "frame mismatch: headline.gross_w_per_kg")
    assert ok, failures
    ok, failures = _has(verdict.failures, "time mismatch: headline.gross_w_per_kg")
    assert ok, failures


# --------------------------------------------------------------------------- regime --------
SCI_REGIMES = frozenset({"level_walking_steady_state"})


def test_expected_regimes_none_is_a_noop():
    """Backwards compatibility: with no ``expected_regimes`` the regime is not inspected."""
    env = _envelope(regime="banana_regime")
    verdict = _validate(env)  # expected_regimes defaults to None
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert not any("regime" in f for f in verdict.failures)


def test_allowed_regime_is_accepted():
    env = _envelope(regime="level_walking_steady_state")
    verdict = _validate(env, expected_regimes=SCI_REGIMES)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert not any("regime" in f for f in verdict.failures)

    # an alias normalises onto the allowed term and is accepted.
    env2 = _envelope(regime="level-walking-steady-state")
    verdict2 = _validate(env2, expected_regimes=SCI_REGIMES)
    assert verdict2.state == "ACCEPT_SCIENTIFIC"
    assert verdict2.failures == ()


def test_wrong_regime_is_refused():
    # a known-but-not-allowed regime (and a recognised out-of-domain one) -> regime mismatch.
    for raw in ("resting", "load_carrying_steady_state", "sprint_anaerobic"):
        env = _envelope(regime=raw)
        verdict = _validate(env, expected_regimes=SCI_REGIMES)
        assert verdict.state == "REFUSE", raw
        ok, failures = _has(verdict.failures, f"regime mismatch: {raw!r}")
        assert ok, (raw, failures)
        assert any("not in" in f for f in verdict.failures), (raw, failures)


def test_missing_regime_is_refused():
    env = _envelope(regime=None)
    verdict = _validate(env, expected_regimes=SCI_REGIMES)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing regime: envelope declares no regime")
    assert ok, failures

    # an absent ``validity.regime`` key is also a missing regime.
    env2 = _envelope()
    del env2["validity"]["regime"]
    verdict2 = _validate(env2, expected_regimes=SCI_REGIMES)
    assert verdict2.state == "REFUSE"
    assert any("missing regime" in f for f in verdict2.failures), verdict2.failures


def test_unknown_regime_is_refused():
    for raw in ("banana_regime", "unknown", ""):
        env = _envelope(regime=raw)
        verdict = _validate(env, expected_regimes=SCI_REGIMES)
        assert verdict.state == "REFUSE", raw
        ok, failures = _has(verdict.failures, f"unknown regime: {raw!r}")
        assert ok, (raw, failures)

    # a non-string regime normalises to unknown and is refused.
    env_nonstr = _envelope(regime=5)
    verdict_nonstr = _validate(env_nonstr, expected_regimes=SCI_REGIMES)
    assert verdict_nonstr.state == "REFUSE"
    assert any("unknown regime: 5" in f for f in verdict_nonstr.failures), verdict_nonstr.failures


def test_regime_check_precedes_provenance_in_failure_order():
    """The regime check sits immediately after the per-quantity time check, before provenance."""
    env = _envelope(provenance=None, regime=None)
    verdict = _validate(env, expected_regimes=SCI_REGIMES)
    assert verdict.state == "REFUSE"
    idx_regime = next(i for i, f in enumerate(verdict.failures) if "missing regime" in f)
    idx_prov = next(i for i, f in enumerate(verdict.failures) if "disallowed provenance" in f)
    assert idx_regime < idx_prov, verdict.failures


def test_absent_producer_provenance_is_refused():
    env = _envelope(provenance=None)
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "disallowed provenance: producer")
    assert ok, failures

    env2 = _envelope()
    del env2["producer"]["provenance"]
    verdict2 = _validate(env2)
    assert verdict2.state == "REFUSE"
    assert any("disallowed provenance: producer" in f for f in verdict2.failures)


def test_per_quantity_disallowed_provenance_is_refused():
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["provenance"] = "synthetic"
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(
        verdict.failures, "disallowed provenance: quantity headline.gross_w_per_kg"
    )
    assert ok, failures


def test_producer_changed_is_refused():
    verdict = _validate(expected_producer_sha256="b" * 64)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "producer changed: source_sha256")
    assert ok, failures


def test_producer_name_mismatch_is_refused():
    verdict = _validate(expected_producer_name="someone_else")
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "producer name mismatch:")
    assert ok, failures


def test_stale_input_is_refused():
    verdict = _validate(current_input_sha256="c" * 64)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "stale input:")
    assert ok, failures


def test_opting_out_of_stale_check_is_refused():
    # A producer must not be able to disable the staleness contract: declaring the flag false
    # is itself a refusal, even when the input hash happens to match.
    env = _envelope()
    env["validity"]["stale_if_input_sha256_mismatch"] = False
    verdict = _validate(env, current_input_sha256=INPUT_SHA)
    assert verdict.state == "REFUSE"
    assert any("stale check disabled" in f for f in verdict.failures), verdict.failures

    # Absent flag is also refused (must be explicitly true).
    env2 = _envelope()
    del env2["validity"]["stale_if_input_sha256_mismatch"]
    verdict2 = _validate(env2)
    assert verdict2.state == "REFUSE"
    assert any("stale check disabled" in f for f in verdict2.failures), verdict2.failures

    # A changed input with the flag still true is refused as stale.
    env3 = _envelope()
    verdict3 = _validate(env3, current_input_sha256="c" * 64)
    assert verdict3.state == "REFUSE"
    assert any("stale input" in f for f in verdict3.failures), verdict3.failures


def test_disallowed_provenance_strict_vs_opt_in():
    env = _envelope(provenance="synthetic")

    strict = _validate(env)
    assert strict.state == "REFUSE"
    assert any("disallowed provenance: producer" in f for f in strict.failures)

    opted = _validate(env, allowed_provenance=rev.allowed_provenance_for(opt_in=("synthetic",)))
    assert opted.state == "ACCEPT_SYNTHETIC_DEMO"
    assert opted.accepted is True
    assert opted.failures == ()

    # "unknown" (absent provenance) can only be admitted by explicit opt-in.
    unknown_env = _envelope(provenance=None)
    assert _validate(unknown_env).state == "REFUSE"
    opted_unknown = _validate(
        unknown_env, allowed_provenance=rev.allowed_provenance_for(opt_in=("unknown",))
    )
    assert opted_unknown.state == "ACCEPT_SYNTHETIC_DEMO"


def test_synthetic_demo_mode_labels_evidence_honestly():
    env = _envelope(provenance="synthetic")
    allowed = rev.allowed_provenance_for(opt_in=("synthetic", "unknown"))
    demo = _validate(env, mode="synthetic_demo", allowed_provenance=allowed)
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert demo.evidence_class == "synthetic_demo"

    # a scientific provenance consumed in demo mode keeps its real label.
    sci = _validate(mode="synthetic_demo")
    assert sci.state == "ACCEPT_SYNTHETIC_DEMO"
    assert sci.evidence_class == "calibrated_public_model"


def test_invalid_mode_raises_value_error():
    with pytest.raises(ValueError):
        _validate(mode="bogus")


def test_malformed_envelope_is_refused():
    for bad in ("not a dict", {}, {"quantities": {}}, [1, 2, 3]):
        verdict = _validate(bad)
        assert verdict.state == "REFUSE", bad
        assert any("malformed envelope" in f for f in verdict.failures), bad


# --------------------------------------------------------------------------- plausibility --
def test_mass_fraction_0585_without_source_is_not_refused():
    """No universal mass-fraction bound: an unusual fraction with no declared source is flagged,
    never refused.  0.585 (local mass / total mass) is physiologically possible."""
    env = _envelope()
    plausibility = [{
        "name": "muscle_mass_fraction",
        "kind": "region_mass_fraction",
        "value": 0.585,
        "fraction": 0.585,
        "local_mass_kg": 0.585 * 72.6,
        "total_mass_kg": 72.6,
        "max_without_source_ok": True,
        # deliberately NO "source"
    }]
    verdict = _validate(env, plausibility=plausibility)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.accepted is True
    assert verdict.failures == ()
    assert not any("physically impossible" in f for f in verdict.failures)


def test_physically_impossible_mass_is_refused():
    env = _envelope()
    plausibility = [{
        "name": "local_muscle_mass",
        "kind": "region_mass_fraction",
        "value": 1.15,
        "fraction": 1.15,
        "local_mass_kg": 84.0,
        "total_mass_kg": 72.6,
        "max_without_source_ok": True,
    }]
    verdict = _validate(env, plausibility=plausibility)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "physically impossible mass")
    assert ok, failures


# ----------------------------------------------- plausibility rule 13 (hardening H4) --------
def test_plausibility_bare_local_total_check_is_quiet_when_possible():
    """A check that declares ONLY local/total (no ``max_fraction``, no ``source``) stays quiet
    when ``local <= total``: it must not append the unsourced-bound flag."""
    env = _envelope()
    plausibility = [{
        "name": "local_working_muscle_mass",
        "kind": "region_mass_fraction",
        "local_mass_kg": 17.2,
        "total_mass_kg": 72.6,
    }]
    verdict = _validate(env, plausibility=plausibility)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert verdict.flags == ()


def test_plausibility_unsourced_flag_requires_declared_max_fraction():
    env = _envelope()
    plausibility = [{
        "name": "muscle_mass_fraction",
        "kind": "region_mass_fraction",
        "local_mass_kg": 0.585 * 72.6,
        "total_mass_kg": 72.6,
        "max_fraction": 0.6,  # declared (sourceless) bound -> non-blocking flag only
    }]
    verdict = _validate(env, plausibility=plausibility)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert any("unsourced plausibility bound not enforced" in f for f in verdict.flags), verdict.flags


def test_plausibility_sourced_bound_flag_is_kept():
    env = _envelope()
    plausibility = [{
        "name": "muscle_mass_fraction",
        "kind": "region_mass_fraction",
        "local_mass_kg": 0.585 * 72.6,
        "total_mass_kg": 72.6,
        "max_fraction": 0.6,
        "source": "some-paper-2001",
    }]
    verdict = _validate(env, plausibility=plausibility)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert any("sourced bound" in f for f in verdict.flags), verdict.flags
    assert not any("unsourced" in f for f in verdict.flags)


# --------------------------------------------------------------------------- receipt --------
def test_receipt_seal_verify_roundtrip(tmp_path):
    producer = tmp_path / "producer.py"
    producer.write_text("# adapted public producer\n", encoding="utf-8")
    payload = {"a": 1, "nested": {"b": [1, 2, 3]}}

    sealed = rev.seal(payload, producer)
    assert set(sealed) == {"payload", "receipt"}
    receipt = sealed["receipt"]
    assert receipt["schema"] == "result_receipt_v1"
    assert set(receipt) == {
        "payload_sha256", "producer_sha256", "producer_path", "sealed_utc", "schema",
    }
    assert rev.verify_receipt(sealed, producer_path=producer) == payload

    sealed_path = tmp_path / "sealed.json"
    rev.write_sealed(sealed_path, payload, producer)
    assert sealed_path.exists()
    assert rev.read_verified(sealed_path, producer_path=producer) == payload


def test_receipt_tampered_payload_is_refused(tmp_path):
    producer = tmp_path / "producer.py"
    producer.write_text("# producer\n", encoding="utf-8")
    sealed = rev.seal({"x": 1}, producer)
    sealed["payload"]["x"] = 2  # tamper, keep the old receipt
    with pytest.raises(ValueError) as excinfo:
        rev.verify_receipt(sealed, producer_path=producer)
    assert "stale result" in str(excinfo.value)
    assert "payload hash mismatch" in str(excinfo.value)


def test_receipt_producer_file_change_is_refused(tmp_path):
    producer = tmp_path / "producer.py"
    producer.write_text("# producer v1\n", encoding="utf-8")
    sealed = rev.seal({"x": 1}, producer)
    producer.write_text("# producer v2 -- changed after seal\n", encoding="utf-8")
    with pytest.raises(ValueError):
        rev.verify_receipt(sealed, producer_path=producer)


def test_receipt_expected_producer_sha_mismatch_is_refused(tmp_path):
    producer = tmp_path / "producer.py"
    producer.write_text("# producer\n", encoding="utf-8")
    sealed = rev.seal({"x": 1}, producer)
    with pytest.raises(ValueError):
        rev.verify_receipt(sealed, expected_producer_sha256="0" * 64)
    # correct sha is accepted
    assert rev.verify_receipt(
        sealed, expected_producer_sha256=sealed["receipt"]["producer_sha256"]
    ) == {"x": 1}


def test_receipt_missing_is_refused():
    with pytest.raises(ValueError):
        rev.verify_receipt({"payload": {"x": 1}})
    with pytest.raises(ValueError):
        rev.verify_receipt(None)
    with pytest.raises(ValueError):
        rev.verify_receipt({"payload": {"x": 1}, "receipt": "not-a-dict"})


# ------------------------------------------------------- abstention is not acceptance ------
# A present quantity whose value is None is an explicit unknown number, never a value.
def test_required_quantity_with_none_value_is_refused_as_abstention():
    env = _envelope()
    env["quantities"]["muscle_mass.total_body_mass_kg"]["value"] = None
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    ok, failures = _has(
        verdict.failures,
        "abstain: required value missing for muscle_mass.total_body_mass_kg",
    )
    assert ok, failures


def test_optional_quantity_with_none_value_is_flagged_not_refused():
    env = _envelope()
    env["quantities"]["headline.net_w_per_kg"]["value"] = None
    expected = _expected()
    expected["headline.net_w_per_kg"]["required"] = False
    verdict = _validate(env, expected_quantities=expected)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert not any("abstain:" in f for f in verdict.failures), verdict.failures
    assert any(
        "abstain: optional value missing for headline.net_w_per_kg" in f
        for f in verdict.flags
    ), verdict.flags


def test_abstention_is_refused_in_synthetic_demo_mode_too():
    """An unknown number is not a number: demo mode is not exempt from the abstain rule."""
    env = _envelope(provenance="synthetic")
    env["quantities"]["muscle_mass.summed_muscle_mass_kg"]["value"] = None
    allowed = rev.allowed_provenance_for(opt_in=("synthetic", "unknown"))
    verdict = _validate(env, mode="synthetic_demo", allowed_provenance=allowed)
    assert verdict.state == "REFUSE"
    assert any("abstain:" in f for f in verdict.failures), verdict.failures


def test_abstention_is_independent_of_units_region_frame_time():
    """A None value is refused even when units/region/frame/time are otherwise all correct."""
    env = _envelope()
    env["quantities"]["headline.gross_w_per_kg"]["value"] = None
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    assert any("abstain:" in f for f in verdict.failures), verdict.failures
    assert not any(
        any(marker in f for marker in ("unit mismatch", "region mismatch", "frame mismatch",
                                       "time mismatch", "missing region_id"))
        for f in verdict.failures
    ), verdict.failures


# ---------------------------------------------------- structured provenance detail --------
# Optional producer.provenance_detail {data, calibration, model, execution}.
def _provenance_detail(**overrides):
    detail = {
        "data": "calibrated",
        "calibration": "literal_cited",  # alias -> calibrated
        "model": "gait2392",
        "execution": {
            "code_sha256": PRODUCER_SHA,
            "command": "python metabolic_cost_adapted.py",
            "seed": 7,
            "host": "cpu",
            "utc": "2026-09-16T00:00:00Z",
        },
    }
    detail.update(overrides)
    return detail


def test_make_envelope_stores_provenance_detail_and_omits_it_when_none():
    env = _envelope(provenance_detail=_provenance_detail())
    assert env["producer"]["provenance_detail"] == _provenance_detail()
    # absent by default: the key is omitted, not None (historical producer shape preserved).
    assert "provenance_detail" not in _envelope()["producer"]


def test_valid_provenance_detail_is_accepted():
    verdict = _validate(_envelope(provenance_detail=_provenance_detail()))
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_absent_provenance_detail_is_accepted_by_default():
    verdict = _validate(_envelope(), require_provenance_detail=False)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_absent_provenance_detail_is_refused_when_required():
    verdict = _validate(_envelope(), require_provenance_detail=True)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing provenance_detail")
    assert ok, failures


def test_present_provenance_detail_satisfies_the_requirement():
    verdict = _validate(
        _envelope(provenance_detail=_provenance_detail()), require_provenance_detail=True
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_unknown_provenance_detail_field_is_refused():
    verdict = _validate(_envelope(provenance_detail=_provenance_detail(bogus=1)))
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "provenance detail: unknown field 'bogus'")
    assert ok, failures


def test_non_object_provenance_detail_block_is_refused():
    verdict = _validate(_envelope(provenance_detail="calibrated"))
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "provenance detail: block must be an object")
    assert ok, failures


@pytest.mark.parametrize("bad", ["assumption", "unknown", None, 5, "not-a-term"])
def test_provenance_detail_bad_data_is_refused(bad):
    verdict = _validate(_envelope(provenance_detail=_provenance_detail(data=bad)))
    assert verdict.state == "REFUSE"
    assert any(
        "provenance detail: data" in f and "is not a declared provenance" in f
        for f in verdict.failures
    ), (bad, verdict.failures)


def test_provenance_detail_bad_calibration_is_refused():
    verdict = _validate(
        _envelope(provenance_detail=_provenance_detail(calibration="assumption"))
    )
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "provenance detail: calibration")
    assert ok, failures


def test_provenance_detail_execution_code_sha256_mismatch_is_refused():
    detail = _provenance_detail(execution={"code_sha256": "b" * 64, "command": "run"})
    verdict = _validate(_envelope(provenance_detail=detail))
    assert verdict.state == "REFUSE"
    assert any(
        "provenance detail: execution.code_sha256" in f
        and "!= producer.source_sha256" in f
        for f in verdict.failures
    ), verdict.failures


def test_provenance_detail_matching_code_sha256_is_accepted():
    detail = _provenance_detail(execution={"code_sha256": PRODUCER_SHA, "command": "run"})
    verdict = _validate(_envelope(provenance_detail=detail))
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_provenance_detail_execution_missing_required_fields_is_refused():
    for bad_exec in ({}, {"code_sha256": PRODUCER_SHA}, {"command": "run"}):
        verdict = _validate(
            _envelope(provenance_detail=_provenance_detail(execution=bad_exec))
        )
        assert verdict.state == "REFUSE", bad_exec
        assert any("provenance detail: execution" in f for f in verdict.failures), (
            bad_exec, verdict.failures,
        )


def test_provenance_detail_model_string_is_accepted_and_empty_dict_refused():
    accepted = _validate(_envelope(provenance_detail=_provenance_detail(model="gait2392")))
    assert accepted.state == "ACCEPT_SCIENTIFIC"
    assert accepted.failures == ()

    refused = _validate(_envelope(provenance_detail=_provenance_detail(model={})))
    assert refused.state == "REFUSE"
    ok, failures = _has(refused.failures, "provenance detail: model")
    assert ok, failures


def test_provenance_detail_model_object_sha256_validation():
    ok = _validate(
        _envelope(
            provenance_detail=_provenance_detail(
                model={"name": "gait2392", "sha256": "a" * 64}
            )
        )
    )
    assert ok.state == "ACCEPT_SCIENTIFIC"

    bad = _validate(
        _envelope(
            provenance_detail=_provenance_detail(model={"name": "gait2392", "sha256": "xyz"})
        )
    )
    assert bad.state == "REFUSE"
    ok_fail, failures = _has(bad.failures, "provenance detail: model.sha256")
    assert ok_fail, failures


# ------------------------------------------------- uncertainty requirement ------------
def _specs_with_uncertainty(value=0.4):
    """``_specs()`` with a declared ``uncertainty`` on every quantity spec."""
    specs = {key: dict(spec) for key, spec in _specs().items()}
    for spec in specs.values():
        spec["uncertainty"] = value
    return specs


def test_require_uncertainty_accepts_when_every_required_quantity_declares_one():
    env = _envelope(quantity_specs=_specs_with_uncertainty())
    verdict = _validate(env, require_uncertainty=True)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert not any("missing uncertainty" in f for f in verdict.failures)


def test_require_uncertainty_refuses_absent_uncertainty():
    # ``_envelope()`` declares no uncertainty on the two headline rates -> the envelope value is
    # None, which under the strict profile is missing.
    verdict = _validate(_envelope(), require_uncertainty=True)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing uncertainty: headline.gross_w_per_kg")
    assert ok, failures
    assert any("missing uncertainty: headline.net_w_per_kg" in f for f in verdict.failures), failures


def test_require_uncertainty_refuses_explicit_none_uncertainty():
    specs = _specs_with_uncertainty()
    specs["muscle_mass.total_body_mass_kg"]["uncertainty"] = None
    verdict = _validate(_envelope(quantity_specs=specs), require_uncertainty=True)
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing uncertainty: muscle_mass.total_body_mass_kg")
    assert ok, failures


def test_require_uncertainty_zero_is_a_valid_declared_uncertainty():
    verdict = _validate(
        _envelope(quantity_specs=_specs_with_uncertainty(0.0)), require_uncertainty=True
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_uncertainty_is_not_required_by_default():
    # The permissive default still accepts an envelope whose quantities carry no uncertainty.
    verdict = _validate(_envelope(), require_uncertainty=False)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert not any("missing uncertainty" in f for f in verdict.failures)


def test_require_uncertainty_ignores_optional_quantities():
    # Only *required* quantities are gated: an optional present quantity with no uncertainty is
    # left alone (no failure).
    specs = _specs_with_uncertainty()
    specs["headline.net_w_per_kg"]["uncertainty"] = None
    expected = _expected()
    expected["headline.net_w_per_kg"]["required"] = False
    verdict = _validate(
        _envelope(quantity_specs=specs),
        expected_quantities=expected,
        require_uncertainty=True,
    )
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()


def test_scientific_requirements_is_a_strict_non_default_profile():
    assert rev.SCIENTIFIC_REQUIREMENTS == {
        "require_uncertainty": True,
        "require_provenance_detail": True,
        "require_producer_gates_pass": True,
    }
    # NOT the default: the historical envelope (no uncertainty, no detail, overall_pass true) is
    # accepted...
    assert _validate(_envelope()).state == "ACCEPT_SCIENTIFIC"
    # ...but refused once the strict profile is spread into the validator.
    verdict = _validate(_envelope(), **rev.SCIENTIFIC_REQUIREMENTS)
    assert verdict.state == "REFUSE"
    assert any("missing uncertainty" in f for f in verdict.failures), verdict.failures
    assert any("missing provenance_detail" in f for f in verdict.failures), verdict.failures
    # This envelope declares overall_pass=True, so the producer-gate branch is satisfied and does
    # NOT fire; the profile still bundles it for envelopes that report a real failure.
    assert not any("producer gates failed" in f for f in verdict.failures), verdict.failures

    # A payload that satisfies the strict profile is accepted.
    strict_env = _envelope(
        quantity_specs=_specs_with_uncertainty(),
        provenance_detail=_provenance_detail(),
    )
    strict_ok = _validate(strict_env, **rev.SCIENTIFIC_REQUIREMENTS)
    assert strict_ok.state == "ACCEPT_SCIENTIFIC"
    assert strict_ok.failures == ()


# ------------------------------------------------ producer gate outcome ------------
def _adapted_style_gates():
    """The healthy adapted public producer's own five gates (two fail -> overall_pass false)."""
    return {
        "peak_rot_speed_plausible": True,
        "speed_sanity_ok": False,
        "no_nan_metabolic_rate": True,
        "proxy_same_order_of_magnitude": True,
        "expected_direction_pass": False,
    }


def test_make_envelope_gate_defaults_are_unchanged():
    """A caller that supplies neither ``gates`` nor ``overall_pass`` gets the historical envelope."""
    env = _envelope()
    assert env["gates"] == {}
    assert env["overall_pass"] is True


def test_make_envelope_stores_the_supplied_gate_outcome():
    gates = _adapted_style_gates()
    env = _envelope(gates=gates, overall_pass=False)
    assert env["gates"] == gates
    assert env["overall_pass"] is False


def test_overall_pass_false_is_refused_by_default_in_strict_mode():
    env = _envelope(gates=_adapted_style_gates(), overall_pass=False)
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    assert "producer gates failed: overall_pass is not true" in verdict.failures
    with pytest.raises(ValueError, match="strict mode cannot waive producer gates"):
        _validate(env, require_producer_gates_pass=False)
    # the individual failing gates are listed by name (sorted), as a second non-blocking flag.
    assert any(
        "producer gates failed: expected_direction_pass, speed_sanity_ok" in f
        for f in verdict.flags
    ), verdict.flags


def test_overall_pass_false_is_refused_when_producer_gates_are_required():
    env = _envelope(gates=_adapted_style_gates(), overall_pass=False)
    verdict = _validate(env, require_producer_gates_pass=True)
    assert verdict.state == "REFUSE"
    assert verdict.accepted is False
    assert any(
        "producer gates failed: overall_pass is not true" in f for f in verdict.failures
    ), verdict.failures
    # the per-gate listing remains a non-blocking flag (it is diagnostic, not a second refusal).
    assert any("producer gates failed:" in f for f in verdict.flags), verdict.flags


def test_overall_pass_true_with_a_failing_gate_is_refused_in_strict_mode():
    """An inconsistent producer (overall_pass=true, a named falsy gate) is refused."""
    env = _envelope(gates={"expected_direction_pass": False, "ok": True}, overall_pass=True)
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    assert "producer gates failed: expected_direction_pass" in verdict.failures
    assert not any("overall_pass is not true" in f for f in verdict.failures + verdict.flags)
    demo = _validate(env, mode="synthetic_demo",
                     allowed_provenance=rev.allowed_provenance_for(opt_in=("synthetic",)))
    assert demo.state == "ACCEPT_SYNTHETIC_DEMO"
    assert "producer gates failed: expected_direction_pass" in demo.flags


def test_missing_overall_pass_is_not_assumed_to_pass():
    """The validator never invents a passing gate set for an envelope that declares none."""
    env = _envelope()
    del env["overall_pass"]
    verdict = _validate(env)
    assert verdict.state == "REFUSE"  # strict mode requires a declared passing gate outcome
    assert "producer gates failed: overall_pass is not true" in verdict.failures

    refused = _validate(env, require_producer_gates_pass=True)
    assert refused.state == "REFUSE"
    assert any("producer gates failed" in f for f in refused.failures), refused.failures


def test_overall_pass_is_independent_of_units_region_frame_time():
    """A gate failure is not a substitute for a contract failure: both are reported."""
    env = _envelope(gates=_adapted_style_gates(), overall_pass=False)
    env["quantities"]["headline.gross_w_per_kg"]["units"] = "W"
    verdict = _validate(env, require_producer_gates_pass=True)
    assert verdict.state == "REFUSE"
    assert any("unit mismatch" in f for f in verdict.failures), verdict.failures
    assert any("producer gates failed" in f for f in verdict.failures), verdict.failures


def test_strict_profile_refuses_the_adapted_style_envelope_for_all_three_reasons():
    """The REAL adapted public control: calibrated_public_model provenance, no declared
    uncertainty, no provenance_detail, and its own producer reports overall_pass=false. The
    strict scientific profile refuses it for all three (evidence) reasons."""
    env = _envelope(gates=_adapted_style_gates(), overall_pass=False)
    verdict = _validate(env, **rev.SCIENTIFIC_REQUIREMENTS)
    assert verdict.state == "REFUSE"
    assert any("missing uncertainty" in f for f in verdict.failures), verdict.failures
    assert any("missing provenance_detail" in f for f in verdict.failures), verdict.failures
    assert any("producer gates failed" in f for f in verdict.failures), verdict.failures


def test_strict_profile_refuses_an_otherwise_complete_envelope_only_for_producer_gates():
    """Satisfy uncertainty + provenance_detail: the producer gate outcome alone still refuses."""
    env = _envelope(
        quantity_specs=_specs_with_uncertainty(),
        provenance_detail=_provenance_detail(),
        gates=_adapted_style_gates(),
        overall_pass=False,
    )
    verdict = _validate(env, **rev.SCIENTIFIC_REQUIREMENTS)
    assert verdict.state == "REFUSE"
    assert not any("missing uncertainty" in f for f in verdict.failures), verdict.failures
    assert not any("missing provenance_detail" in f for f in verdict.failures), verdict.failures
    assert any("producer gates failed" in f for f in verdict.failures), verdict.failures


def test_gate_flags_are_surfaced_even_when_other_failures_refuse():
    """Flags accompany a REFUSE too: a refusing envelope still reports the producer's gates."""
    env = _envelope(gates=_adapted_style_gates(), overall_pass=False)
    env["producer"]["provenance"] = None  # a separate, dominant refusal
    verdict = _validate(env)
    assert verdict.state == "REFUSE"
    assert any("disallowed provenance" in f for f in verdict.failures), verdict.failures
    assert any("producer gates failed" in f for f in verdict.flags), verdict.flags


# ------------------------------------------- material-region binding ---------------
# A quantity may carry an optional ``material_region_id`` (an explicit region -> material binding).
# ``make_envelope`` copies it only when the spec provides the key; ``validate_for_consumer`` accepts
# either a MAPPING (exact per-key binding) or a SET/frozenset (a vocabulary over all quantities).
def _material_specs(material_region_id="MSK-MUSCLE"):
    specs = {key: dict(spec) for key, spec in _specs().items()}
    for spec in specs.values():
        spec["material_region_id"] = material_region_id
    return specs


def test_make_envelope_copies_material_region_id_only_when_the_spec_provides_it():
    # A spec that provides the key gets it copied; one that does not keeps the historical shape.
    specs = _specs()
    specs["muscle_mass.total_body_mass_kg"]["material_region_id"] = "MSK-BONE"
    env = _envelope(quantity_specs=specs)
    assert env["quantities"]["muscle_mass.total_body_mass_kg"]["material_region_id"] == "MSK-BONE"
    assert "material_region_id" not in env["quantities"]["headline.gross_w_per_kg"]

    # no spec entry -> the key is OMITTED (historical seven-key quantity shape preserved).
    plain = _envelope()
    assert "material_region_id" not in plain["quantities"]["muscle_mass.total_body_mass_kg"]

    # an explicit spec entry of None is still a *provided* key and is copied as None.
    specs = _specs()
    specs["headline.gross_w_per_kg"]["material_region_id"] = None
    explicit = _envelope(quantity_specs=specs)
    assert "material_region_id" in explicit["quantities"]["headline.gross_w_per_kg"]
    assert explicit["quantities"]["headline.gross_w_per_kg"]["material_region_id"] is None


def test_expected_material_regions_mapping_match_is_accepted():
    env = _envelope(quantity_specs=_material_specs("MSK-MUSCLE"))
    expected_map = {key: "MSK-MUSCLE" for key in _expected()}
    verdict = _validate(env, expected_material_regions=expected_map)
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert not any("material region" in f for f in verdict.failures)


def test_expected_material_regions_mapping_mismatch_is_refused():
    env = _envelope(quantity_specs=_material_specs("MSK-MUSCLE"))
    expected_map = {key: "MSK-MUSCLE" for key in _expected()}
    expected_map["headline.gross_w_per_kg"] = "MSK-BONE"  # swapped binding
    verdict = _validate(env, expected_material_regions=expected_map)
    assert verdict.state == "REFUSE"
    ok, failures = _has(
        verdict.failures,
        "material region mismatch: headline.gross_w_per_kg expected 'MSK-BONE' got 'MSK-MUSCLE'",
    )
    assert ok, failures


def test_expected_material_regions_mapping_absent_is_missing():
    # The envelope declares no material binding at all.
    verdict = _validate(
        _envelope(),
        expected_material_regions={"headline.gross_w_per_kg": "MSK-MUSCLE"},
    )
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "missing material region: headline.gross_w_per_kg")
    assert ok, failures


def test_expected_material_regions_set_vocabulary_accepts_members_and_refuses_unknown():
    vocabulary = frozenset({"MSK-MUSCLE", "MSK-BONE"})
    accepted = _validate(
        _envelope(quantity_specs=_material_specs("MSK-MUSCLE")),
        expected_material_regions=vocabulary,
    )
    assert accepted.state == "ACCEPT_SCIENTIFIC"
    assert accepted.failures == ()

    refused = _validate(
        _envelope(quantity_specs=_material_specs("MSK-NOPE")),
        expected_material_regions=vocabulary,
    )
    assert refused.state == "REFUSE"
    ok, failures = _has(
        refused.failures,
        "unknown material region: muscle_mass.total_body_mass_kg 'MSK-NOPE'",
    )
    assert ok, failures


def test_expected_material_regions_set_applies_to_quantities_outside_the_contract():
    """Mirroring the region-vocabulary rule: the set form also covers quantities the consumer did
    not declare in ``expected_quantities``."""
    env = _envelope(quantity_specs=_material_specs("MSK-MUSCLE"))
    env["quantities"]["extra.unlisted"] = {
        "value": 1.0,
        "units": "x",
        "region_id": None,
        "frame": "scalar",
        "time": None,
        "uncertainty": None,
        "provenance": "calibrated_public_model",
        "material_region_id": "MSK-NOPE",
    }
    verdict = _validate(env, expected_material_regions=frozenset({"MSK-MUSCLE"}))
    assert verdict.state == "REFUSE"
    ok, failures = _has(verdict.failures, "unknown material region: extra.unlisted 'MSK-NOPE'")
    assert ok, failures


def test_expected_material_regions_none_is_a_noop():
    """Backwards compatibility: with no ``expected_material_regions`` the binding is not inspected,
    even when a quantity declares an otherwise unregistered id."""
    verdict = _validate(_envelope(quantity_specs=_material_specs("MSK-NOPE")))
    assert verdict.state == "ACCEPT_SCIENTIFIC"
    assert verdict.failures == ()
    assert not any("material region" in f for f in verdict.failures)


def test_material_region_vocabulary_applies_independently_of_the_mapping():
    env = _envelope(quantity_specs=_material_specs("MSK-NOPE"))
    verdict = _validate(env, material_region_vocabulary=frozenset({"MSK-MUSCLE"}))
    assert verdict.state == "REFUSE"
    ok, failures = _has(
        verdict.failures, "unknown material region: muscle_mass.total_body_mass_kg 'MSK-NOPE'"
    )
    assert ok, failures


# ------------------------------------------------ trusted pin / literal_cited / demo path ------
def test_strict_acceptance_requires_the_trusted_artifact_pin():
    env = _envelope()
    missing = _validate(env, expected_payload_sha256=None)
    assert missing.state == "REFUSE"
    assert ("no trusted artifact pin: expected_payload_sha256 is required for scientific "
            "acceptance") in missing.failures
    other = _validate(env, expected_payload_sha256="0" * 64)
    assert other.state == "REFUSE"
    assert any(f.startswith("artifact not approved: payload_sha256") for f in other.failures)
    assert _validate(env).state == "ACCEPT_SCIENTIFIC"


def test_edited_payload_with_recomputed_receipt_is_refused_by_the_pin():
    env = _envelope()
    pin = rev.payload_sha256(env)
    sealed = {"payload": env, "receipt": {"payload_sha256": pin}}
    sealed["payload"]["quantities"]["muscle_mass.total_body_mass_kg"]["value"] = 726.0
    sealed["receipt"]["payload_sha256"] = rev.payload_sha256(sealed["payload"])
    assert rev.verify_receipt(sealed) is sealed["payload"]  # integrity check only: passes
    verdict = _validate(sealed["payload"], expected_payload_sha256=pin)
    assert verdict.state == "REFUSE"
    assert any(f.startswith("artifact not approved") for f in verdict.failures)


def test_literal_cited_is_not_scientific_provenance():
    env = _envelope(provenance="literal_cited")
    for spec in env["quantities"].values():
        spec["provenance"] = "literal_cited"
    refused = _validate(env)
    assert refused.state == "REFUSE"
    assert any("disallowed provenance: producer" in f and "literal_cited" in f
               for f in refused.failures), refused.failures
    opted = _validate(env, allowed_provenance=rev.allowed_provenance_for(opt_in=("literal_cited",)))
    assert opted.state == "ACCEPT_SYNTHETIC_DEMO"


def test_demo_mode_never_yields_scientific_acceptance():
    env = _envelope()
    verdict = _validate(env, mode="synthetic_demo")
    assert verdict.state == "ACCEPT_SYNTHETIC_DEMO"
    assert verdict.state != "ACCEPT_SCIENTIFIC"
