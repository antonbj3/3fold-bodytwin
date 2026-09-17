"""Acceptance tests for the framework-side region-vocabulary mirror (``region_vocabulary_v1``).

The declared vocabulary is only useful if it matches the I2 roster in this repository
(``bodytwin.geometry.material_roster_v1``). The parity check must hold; an import failure of
the roster is a test failure.
"""
import pytest

import bodytwin.framework.region_vocabulary_v1 as rvf

import bodytwin.framework.consumer_preflight_v1 as cp  # noqa: F401

# The registered region ids (src/bodytwin/geometry/material_roster_v1.py).
I2_ROSTER_SHA256 = "1998eec14cd3195f35705fb3d08324a765b02e23e9fde31a0fabbfb69e08ae2a"
FROZEN_IDS = (
    "CARD-VESSEL",
    "GENERIC",
    "MSK-BONE",
    "MSK-MUSCLE",
    "MSK-TENDON",
    "NEU-TISSUE",
    "WHOLE-BODY",
)


def test_declared_ids_are_the_frozen_i2_tuple():
    assert rvf.I2_REGION_IDS == FROZEN_IDS
    assert rvf.I2_ROSTER_SHA256 == I2_ROSTER_SHA256
    assert rvf.I2_ROSTER_PATH.name == "material_roster_v1.py"
    assert rvf.I2_ROSTER_PATH.as_posix().endswith("src/bodytwin/geometry/material_roster_v1.py")
    assert rvf.I2_ROSTER_MODULE == "bodytwin.geometry.material_roster_v1"


def test_in_repo_roster_file_matches_the_recorded_sha256():
    import hashlib
    assert rvf.I2_ROSTER_PATH.is_file()
    assert hashlib.sha256(rvf.I2_ROSTER_PATH.read_bytes()).hexdigest() == I2_ROSTER_SHA256


def test_declared_ids_match_i2_roster():
    status, message = rvf.assert_i2_parity()
    assert status is True, message
    assert rvf.load_i2_region_ids() == FROZEN_IDS


def test_thermo_region_vocabulary_equals_the_i2_declared_set():
    assert cp.THERMO_REGION_VOCABULARY == frozenset(rvf.I2_REGION_IDS)
    assert cp.THERMO_REGION_VOCABULARY == frozenset(FROZEN_IDS)


def test_msk_lowerlimb_hip_is_not_a_registered_region():
    assert "MSK-LOWERLIMB-HIP" not in cp.THERMO_REGION_VOCABULARY
    assert "MSK-LOWERLIMB-HIP" not in rvf.I2_REGION_IDS
    # kept only as a named placeholder for the placeholder-era tests.
    assert cp.THERMO_DECLARED_SYNTHETIC_REGION == "MSK-LOWERLIMB-HIP"


def test_load_returns_none_for_a_missing_or_broken_path(tmp_path):
    assert rvf.load_i2_region_ids(tmp_path / "does_not_exist.py") is None
    assert rvf.assert_i2_parity(tmp_path / "does_not_exist.py") == (None, "I2 roster not importable")
