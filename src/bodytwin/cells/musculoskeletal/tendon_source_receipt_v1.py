"""Provenance + unit binding for the parsed UniPD tendon data.

Why this module exists
----------------------
The relaxation analysis (``tendon_relaxation_v1``) consumes a parsed copy of the UniPD
human foot-tendon series. Those numbers are only meaningful if the *source* they came from, the
*exact bytes* that were parsed, and the *declared units* are pinned to the
analysis. This module freezes that binding as a single immutable
:class:`SourceReceipt` and provides the only sanctioned path for converting the
recorded gram-force values to absolute newtons.

Honesty boundary
----------------
A verified receipt licenses *absolute force* in newtons. It does NOT license
absolute *material* parameters (elastic modulus, viscosity, physical relaxation
time). Those additionally require specimen geometry (cross-sectional area,
reference/gauge length) and a verified displacement/strain protocol. The UniPD
release carries neither (``geometry_present=False``, ``displacement_present=False``),
so :func:`absolute_parameter_prerequisites` returns ABSTAIN and must keep doing
so. No default geometry is ever substituted.

Scope
-----
Pure Python/NumPy. Streams file hashes in chunks. Performs no network access and
writes nothing; the caller supplies the pack root explicitly (a directory holding
``data/`` with the workbook and the parsed npz, and ``evidence/`` with the parse metadata).
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
import sys
import tempfile
from dataclasses import dataclass, fields, replace
from typing import Any

import numpy as np

# ------------------------------------------------------------- sibling module ----
# Reuse the relaxation analyzer (read-only). As a package module the sibling is
# imported normally; when this file is run as a standalone script (the repo's
# cell runner does that) the src root is put on sys.path and the package import
# is used as well, so both paths share one module object.
_HERE = os.path.dirname(os.path.abspath(__file__))
if __package__:
    from . import tendon_relaxation_v1 as m
else:  # standalone script run
    _SRC_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
    if _SRC_ROOT not in sys.path:
        sys.path.insert(0, _SRC_ROOT)
    from bodytwin.cells.musculoskeletal import tendon_relaxation_v1 as m

__all__ = [
    "SourceReceipt", "UNIPD_RECEIPT", "RECEIPT_FIELDS",
    "sha256", "verify_receipt", "raise_if_invalid", "load_normalized_holds",
    "absolute_force_newton", "absolute_parameter_prerequisites", "self_test",
]


# ------------------------------------------------------------------- receipt ----
@dataclass(frozen=True)
class SourceReceipt:
    """Immutable provenance+units binding for one parsed data release.

    Every field is either an identity/attribution string, a byte/hash pin, a
    declared unit, or an observed count/quality statistic. ``geometry_present``
    and ``displacement_present`` are explicit capability flags: False means the
    release cannot support absolute material parameters, not that the value is
    merely unknown.
    """

    dataset_id: str
    title: str
    doi: str
    url: str
    licence: str
    retrieved: str
    file_name: str
    file_bytes: int
    file_sha256: str
    extracted_npz_name: str
    extracted_npz_sha256: str
    extracted_meta_name: str
    extracted_meta_sha256: str
    force_unit: str
    time_unit: str
    position_unit: str
    n_specimens: int
    n_holds: int
    n_donors: int
    force_quantization_gf_median_min_nonzero: float
    geometry_present: bool
    displacement_present: bool
    notes: str


RECEIPT_FIELDS = tuple(f.name for f in fields(SourceReceipt))


# Frozen literal values for the UniPD release. These are transcribed, never
# derived at import time, so a mismatch is detectable rather than self-consistent.
UNIPD_RECEIPT = SourceReceipt(
    dataset_id="unipd_foot_tendon_1897",
    title="Mechanical tests on human fresh-frozen foot tendons",
    doi="10.25430/researchdata.cab.unipd.it.00001897",
    url="https://researchdata.cab.unipd.it/1897/",
    licence="CC BY 4.0",
    retrieved="2026-09-16",
    file_name="experimental_data.xlsx",
    file_bytes=21717298,
    file_sha256="55e367b6fc3f204b03a6f749ffcfe3495a1c0a765e71005fc0ddfa42f1e680c3",
    extracted_npz_name="unipd_series.npz",
    extracted_npz_sha256="b670d43dd6e3d15e737d4ce4e7fbff68f858f6cf22bc6abc0c857b7301d0a0af",
    extracted_meta_name="unipd_series_meta.json",
    extracted_meta_sha256="8fa730806a0f30954ea3128f431ff8ad63d5e768567fc9ce8fcf10f8ab6708cd",
    force_unit="gf",
    time_unit="s",
    position_unit="mm",
    n_specimens=18,
    n_holds=91,
    n_donors=4,
    force_quantization_gf_median_min_nonzero=1.26312,
    geometry_present=False,
    displacement_present=False,
    notes=("UniPD human foot-tendon mechanical tests (Fontanella et al.). Force is "
           "gram-force, quantised with a median minimum non-zero step of ~1.263 gf. "
           "Failure sheets carry actuator position but no time; stress-relaxation "
           "sheets carry time but no displacement. No cross-sectional area, gauge "
           "length or strain protocol -> absolute material parameters ABSTAIN; only "
           "the normalized relaxation shape is identifiable."),
)


# ---------------------------------------------------------------------- hash ----
def sha256(path: str) -> str:
    """Streamed SHA-256 of a file (1 MiB chunks), returned as lowercase hex."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------- verification ----
def _pack_paths(receipt: SourceReceipt, pack_root: str) -> dict[str, str]:
    """Pack layout: ``<pack_root>/data/<xlsx, npz>`` and ``<pack_root>/evidence/<meta json>``."""
    data_dir = os.path.join(pack_root, "data")
    evidence_dir = os.path.join(pack_root, "evidence")
    return {
        "xlsx": os.path.join(data_dir, receipt.file_name),
        "npz": os.path.join(data_dir, receipt.extracted_npz_name),
        "meta": os.path.join(evidence_dir, receipt.extracted_meta_name),
    }


def verify_receipt(receipt: SourceReceipt, pack_root: str) -> dict:
    """Check a receipt against the on-disk pack. Never raises on mismatch.

    Returns ``{"checks": {name: bool}, "all_ok": bool, "missing_files": [...]}``.
    A mismatch is reported as ``False``; only :func:`raise_if_invalid` turns that
    into an exception.
    """
    paths = _pack_paths(receipt, pack_root)
    checks: dict[str, bool] = {
        "xlsx_present": os.path.isfile(paths["xlsx"]),
        "npz_present": os.path.isfile(paths["npz"]),
        "meta_present": os.path.isfile(paths["meta"]),
        "xlsx_bytes_match": False,
        "xlsx_sha256_match": False,
        "npz_sha256_match": False,
        "meta_sha256_match": False,
    }
    missing = [p for name in ("xlsx", "npz", "meta")
               for p in (paths[name],) if not os.path.isfile(p)]
    try:
        if checks["xlsx_present"]:
            checks["xlsx_bytes_match"] = bool(
                os.path.getsize(paths["xlsx"]) == receipt.file_bytes)
            checks["xlsx_sha256_match"] = bool(
                sha256(paths["xlsx"]) == receipt.file_sha256)
        if checks["npz_present"]:
            checks["npz_sha256_match"] = bool(
                sha256(paths["npz"]) == receipt.extracted_npz_sha256)
        if checks["meta_present"]:
            checks["meta_sha256_match"] = bool(
                sha256(paths["meta"]) == receipt.extracted_meta_sha256)
    except OSError as exc:  # pragma: no cover - race/IO; still report, not raise
        checks["io_error"] = False
        missing.append(str(exc))
    return {
        "checks": checks,
        "all_ok": bool(all(checks.values())),
        "missing_files": missing,
        "dataset_id": receipt.dataset_id,
        "pack_root": pack_root,
    }


def raise_if_invalid(verdict: dict) -> dict:
    """Strict wrapper: raise ValueError unless ``verdict['all_ok']`` is True."""
    if not verdict.get("all_ok"):
        failed = [k for k, v in verdict.get("checks", {}).items() if not v]
        raise ValueError(
            "source receipt verification failed: "
            f"failed_checks={failed} missing_files={verdict.get('missing_files')}")
    return verdict


# ------------------------------------------------------------------- loading ----
def load_normalized_holds(pack_root: str, receipt: SourceReceipt = UNIPD_RECEIPT,
                          nmax: int = 300) -> tuple[list, dict, SourceReceipt]:
    """Verify the receipt strictly, then load the normalized holds via the adapter.

    Returns ``(holds, info, receipt)`` where ``holds``/``info`` are exactly what
    ``tendon_relaxation_v1.load_unipd_holds`` returns. The declared force unit must
    be ``"gf"``; anything else is refused rather than silently interpreted.
    """
    raise_if_invalid(verify_receipt(receipt, pack_root))
    if receipt.force_unit != "gf":
        raise ValueError(
            f"load_normalized_holds requires a 'gf'-declared receipt, got "
            f"{receipt.force_unit!r}")
    paths = _pack_paths(receipt, pack_root)
    holds, info = m.load_unipd_holds(paths["npz"], paths["meta"], nmax=nmax)
    return holds, info, receipt


# ------------------------------------------------------------- absolute force ----
def absolute_force_newton(force_gf: Any,
                          receipt: SourceReceipt = UNIPD_RECEIPT) -> np.ndarray:
    """The ONLY sanctioned absolute-force path: verified receipt -> newtons.

    An absolute force in newtons requires the verified source receipt (exact
    bytes, declared ``gf`` frame). Absolute *material* parameters additionally
    require geometry (cross-sectional area), a reference length and a verified
    strain/displacement protocol; for UniPD those are absent, so material
    parameters remain ABSTAIN regardless of this conversion.

    The declared frame must be ``"gf"`` (ValueError otherwise). For a real force
    *trace* (more than one sample) the module's magnitude band guard
    ``assert_declared_force_frame`` is applied to catch a gf/N frame swap. A lone
    reference scalar is not a trace and is converted directly with the exact
    factor.
    """
    if receipt.force_unit != "gf":
        raise ValueError(
            f"absolute_force_newton requires a 'gf'-declared receipt, got "
            f"{receipt.force_unit!r}")
    arr = np.asarray(force_gf, dtype=float)
    if arr.ndim == 0 or arr.size == 1:
        return m.gf_to_newton(arr)
    m.assert_declared_force_frame(arr, receipt.force_unit)
    return m.gf_to_newton(arr)


# -------------------------------------------------------- material prereqs gate ----
def absolute_parameter_prerequisites(receipt: SourceReceipt = UNIPD_RECEIPT) -> dict:
    """ABSTAIN gate for absolute material parameters, driven by the receipt flags.

    The receipt stores capability flags, not geometry values. Because
    ``geometry_present=False`` (no area, no reference length) and
    ``displacement_present=False`` (no strain protocol/channel), every input to
    :func:`tendon_relaxation_v1.assess_geometry_prerequisites` is None and the
    verdict is ABSTAIN. It must stay ABSTAIN for UniPD.
    """
    geometry = bool(receipt.geometry_present)
    displacement = bool(receipt.displacement_present)
    # The receipt deliberately carries no geometry numbers. If geometry were
    # present we would still have to refuse, because a flag is not a measurement;
    # the caller would need a receipt variant holding the actual values.
    area_mm2 = None
    reference_length_mm = None
    strain_protocol_verified = None
    displacement_channel = None
    if geometry or displacement:  # pragma: no cover - never true for UNIPD_RECEIPT
        raise ValueError(
            "receipt claims geometry/displacement present but stores no values; "
            "extend the receipt before assessing absolute parameters")
    verdict = dict(m.assess_geometry_prerequisites(
        area_mm2=area_mm2,
        reference_length_mm=reference_length_mm,
        strain_protocol_verified=strain_protocol_verified,
        displacement_channel=displacement_channel,
    ))
    verdict["receipt_dataset_id"] = receipt.dataset_id
    verdict["receipt_geometry_present"] = geometry
    verdict["receipt_displacement_present"] = displacement
    return verdict


# ------------------------------------------------------------------ self-test ----
def self_test() -> dict:
    """Hermetic self-test. No network, no writes outside a temp dir."""
    out: dict[str, bool] = {}
    r = UNIPD_RECEIPT

    out["receipt_fields_present"] = bool(
        all(getattr(r, name) is not None for name in RECEIPT_FIELDS))
    out["receipt_counts_consistent"] = bool(
        r.n_specimens == 18 and r.n_holds == 91 and r.n_donors == 4
        and r.force_unit == "gf" and r.time_unit == "s"
        and r.position_unit == "mm"
        and abs(r.force_quantization_gf_median_min_nonzero - 1.26312) < 1e-12
        and r.geometry_present is False and r.displacement_present is False)
    out["receipt_fields_expected"] = bool(
        set(RECEIPT_FIELDS) == {
            "dataset_id", "title", "doi", "url", "licence", "retrieved",
            "file_name", "file_bytes", "file_sha256", "extracted_npz_name",
            "extracted_npz_sha256", "extracted_meta_name",
            "extracted_meta_sha256", "force_unit", "time_unit", "position_unit",
            "n_specimens", "n_holds", "n_donors",
            "force_quantization_gf_median_min_nonzero", "geometry_present",
            "displacement_present", "notes"})

    # exact gf -> N for 1.0 through both the primitive and the receipt path
    out["gf_to_newton_one_is_9p80665e_minus3"] = bool(
        float(m.gf_to_newton(1.0)) == m.GF_TO_N == 9.80665e-3)
    out["absolute_force_newton_one_exact"] = bool(
        float(absolute_force_newton(1.0, r)) == 9.80665e-3)

    # a "N"-declared receipt is refused
    bad_unit = replace(r, force_unit="N")
    try:
        absolute_force_newton([1.0, 100.0], bad_unit)
        out["refuses_n_receipt"] = False
    except ValueError:
        out["refuses_n_receipt"] = True

    # corrupted expected hash is detected; strict wrapper raises
    with tempfile.TemporaryDirectory() as td:
        data_dir = os.path.join(td, "data")
        ev_dir = os.path.join(td, "evidence")
        os.makedirs(data_dir)
        os.makedirs(ev_dir)
        xlsx = os.path.join(data_dir, r.file_name)
        npz = os.path.join(data_dir, r.extracted_npz_name)
        meta = os.path.join(ev_dir, r.extracted_meta_name)
        for p, blob in ((xlsx, b"xlsx-bytes" * 10), (npz, b"npz-bytes"),
                        (meta, b"{}")):
            with open(p, "wb") as fh:
                fh.write(blob)
        good = replace(r, file_bytes=os.path.getsize(xlsx),
                       file_sha256=sha256(xlsx),
                       extracted_npz_sha256=sha256(npz),
                       extracted_meta_sha256=sha256(meta))
        v_good = verify_receipt(good, td)
        corrupted = replace(good, file_sha256="0" * 64)
        v_bad = verify_receipt(corrupted, td)
        out["verify_receipt_all_ok_on_valid"] = bool(v_good["all_ok"])
        out["corrupted_hash_detected"] = bool(
            v_bad["all_ok"] is False
            and v_bad["checks"]["xlsx_sha256_match"] is False
            and v_bad["checks"]["xlsx_bytes_match"] is True)
        try:
            raise_if_invalid(v_bad)
            out["raise_on_invalid"] = False
        except ValueError:
            out["raise_on_invalid"] = True
    return out


if __name__ == "__main__":  # pragma: no cover - manual probe
    import json
    print(json.dumps(self_test(), indent=2))
