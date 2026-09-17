#!/usr/bin/env python3
"""Tests for the I3 source-receipt module.

Plain python3:  OMP_NUM_THREADS=1 PYTHONPATH=src python tests/test_tendon_source_receipt_v1.py
Also pytest-compatible (functions named test_*).

Covers the acceptance items that can be tested in isolation:
  * receipt field completeness, frozen-ness and internal consistency
  * verify_receipt all_ok on the real UniPD pack (hashes + byte count; needs
    BODYTWIN_REF_UNIPD_PACK, skipped otherwise)
  * a corrupted expected hash is detected and the strict wrapper raises
  * load_normalized_holds reproduces 18 specimens / 91 holds / 4 donors
  * gf -> N is exact (1.0 gf == 9.80665e-3 N) through the receipt path
  * absolute_force_newton refuses a "N"-declared receipt
  * absolute material parameters ABSTAIN with the four missing prerequisites
  * the module self_test is fully green
"""
import dataclasses
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULE = os.path.join(ROOT, "src", "bodytwin", "cells", "musculoskeletal",
                      "tendon_source_receipt_v1.py")
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
import numpy as np  # noqa: E402

from bodytwin.cells.musculoskeletal import tendon_source_receipt_v1 as m  # noqa: E402


def _require_pack():
    """The parsed UniPD pack (BODYTWIN_REF_UNIPD_PACK; see tests/research_data.py)."""
    from research_data import optional_path

    return str(optional_path("BODYTWIN_REF_UNIPD_PACK"))


# --------------------------------------------------------------- receipt pin ----
def test_receipt_fields_complete_and_consistent():
    r = m.UNIPD_RECEIPT
    for name in m.RECEIPT_FIELDS:
        assert getattr(r, name) is not None, name
    assert r.dataset_id == "unipd_foot_tendon_1897"
    assert r.doi == "10.25430/researchdata.cab.unipd.it.00001897"
    assert r.url == "https://researchdata.cab.unipd.it/1897/"
    assert r.licence == "CC BY 4.0"
    assert r.file_name == "experimental_data.xlsx"
    assert r.file_bytes == 21717298
    assert r.force_unit == "gf" and r.time_unit == "s" and r.position_unit == "mm"
    assert r.n_specimens == 18 and r.n_holds == 91 and r.n_donors == 4
    assert abs(r.force_quantization_gf_median_min_nonzero - 1.26312) < 1e-12
    assert r.geometry_present is False
    assert r.displacement_present is False
    # the dataclass is frozen: mutation must be refused
    try:
        r.n_holds = 0  # type: ignore[misc]
        raise AssertionError("expected FrozenInstanceError")
    except dataclasses.FrozenInstanceError:
        pass


def test_receipt_hashes_are_hex64():
    r = m.UNIPD_RECEIPT
    for value in (r.file_sha256, r.extracted_npz_sha256, r.extracted_meta_sha256):
        assert len(value) == 64 and all(c in "0123456789abcdef" for c in value)


# ------------------------------------------------------------- verification ----
def test_verify_receipt_all_ok_on_real_pack():
    PACK_ROOT = _require_pack()
    verdict = m.verify_receipt(m.UNIPD_RECEIPT, PACK_ROOT)
    assert verdict["all_ok"] is True, verdict
    assert verdict["missing_files"] == []
    assert verdict["checks"]["xlsx_bytes_match"] is True
    assert verdict["checks"]["xlsx_sha256_match"] is True
    assert verdict["checks"]["npz_sha256_match"] is True
    assert verdict["checks"]["meta_sha256_match"] is True


def test_corrupted_hash_is_detected_and_strict_wrapper_raises():
    PACK_ROOT = _require_pack()
    corrupted = dataclasses.replace(m.UNIPD_RECEIPT, file_sha256="0" * 64)
    verdict = m.verify_receipt(corrupted, PACK_ROOT)
    assert verdict["all_ok"] is False
    assert verdict["checks"]["xlsx_sha256_match"] is False
    # verify must NOT raise; only the strict wrapper does
    try:
        m.raise_if_invalid(verdict)
        raise AssertionError("expected ValueError from raise_if_invalid")
    except ValueError:
        pass
    # and the valid verdict passes the strict wrapper
    m.raise_if_invalid(m.verify_receipt(m.UNIPD_RECEIPT, PACK_ROOT))


# ------------------------------------------------------------------- loading ----
def test_load_normalized_holds_counts():
    PACK_ROOT = _require_pack()
    holds, info, receipt = m.load_normalized_holds(PACK_ROOT)
    assert receipt is m.UNIPD_RECEIPT
    assert len(holds) == 91
    assert info["n_specimens"] == 18
    assert info["n_holds"] == 91
    assert info["n_donors"] == 4
    assert len(set(info["donors"])) == 4  # donor labels are not repeated here
    # the adapter's decimation appends the final sample if the stride misses it,
    # so a hold can carry at most nmax+1 points.
    assert all(h["x"].size <= 301 for h in holds)


# ------------------------------------------------------------ absolute force ----
def test_absolute_force_newton_one_gf_exact():
    got = m.absolute_force_newton(1.0, m.UNIPD_RECEIPT)
    assert float(got) == 9.80665e-3
    assert float(got) == m.m.GF_TO_N
    # a real trace takes the frame-guarded path and still scales exactly
    trace = m.absolute_force_newton(np.array([1000.0, 500.0]), m.UNIPD_RECEIPT)
    assert float(trace[0]) == 1000.0 * m.m.GF_TO_N


def test_absolute_force_newton_refuses_n_receipt():
    bad = dataclasses.replace(m.UNIPD_RECEIPT, force_unit="N")
    try:
        m.absolute_force_newton(np.array([1000.0, 500.0]), bad)
        raise AssertionError("expected refusal for a 'N'-declared receipt")
    except ValueError:
        pass


# --------------------------------------------------------- abstain material ----
def test_absolute_parameter_prerequisites_abstains():
    v = m.absolute_parameter_prerequisites(m.UNIPD_RECEIPT)
    assert v["state"] == "ABSTAIN"
    assert set(v["missing"]) == {"area_mm2", "reference_length_mm",
                                 "strain_protocol_verified", "displacement_channel"}
    assert v["present"] == []
    assert v["receipt_geometry_present"] is False
    assert v["receipt_displacement_present"] is False


# ------------------------------------------------------------------ self-test ----
def test_module_self_test_all_true():
    out = m.self_test()
    assert out and all(out.values()), out


def _run_all():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    results, failed = [], []
    for t in tests:
        try:
            t()
            results.append((t.__name__, "PASS"))
        except Exception as exc:  # noqa: BLE001
            results.append((t.__name__, f"FAIL: {exc}"))
            failed.append(t.__name__)
    for name, res in results:
        if res == "PASS":
            print(f"{res:60s} {name}")
        else:
            print(f"{name}: {res}")
    print(f"\n{len(tests) - len(failed)}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_run_all())
