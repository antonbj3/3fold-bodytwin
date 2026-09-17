"""ResultEnvelope v1 -- smallest reusable cached-result hand-off wrapper + consumer validator.

It combines three earlier designs from the 2026-09 BodyTwin review (none of them part of
this repository):

* an envelope prototype (envelope shape, provenance vocabulary, validate-before-use rule);
* a result-receipt prototype (payload + producer hash receipt);
* a geometry contract (a frame is either ``scalar`` or a JSON-canonicalisable descriptor;
  units/region are explicit and never inferred).

Acceptance vs integrity: :func:`validate_for_consumer` in ``strict`` mode is the only path that
can return ``ACCEPT_SCIENTIFIC``. It requires (among the checks below) that the producer's own
gates passed and that the envelope is the artifact an independently trusted registry approved
(``expected_payload_sha256``, compared with :func:`payload_sha256`). ``synthetic_demo`` mode never
returns ``ACCEPT_SCIENTIFIC``. The receipt functions (:func:`seal`, :func:`verify_receipt`,
:func:`read_verified`) check integrity only: an unkeyed sha256 detects accidental change or a stale
payload, but anyone can recompute it, so a receipt alone never establishes that an artifact was
approved. :func:`read_verified_chain` is the consumer-side chain reader and requires the
independently expected upstream payload hash.

Scope, deliberately small and stdlib-only (``json, hashlib, dataclasses, datetime, pathlib``):

* :func:`make_envelope` wraps a raw (possibly dotted/nested) result dict into the envelope.
* :func:`validate_for_consumer` is the consumer preflight: pure, no I/O, no framework.
* :func:`seal` / :func:`verify_receipt` / :func:`write_sealed` / :func:`read_verified` add a
  tamper-evident receipt around any payload; :func:`seal_chain` / :func:`read_verified_chain`
  bind a downstream receipt to a verified upstream artifact (``receipt["upstream"]``).

Schema::

    ResultEnvelope:
      envelope_version: int
      producer:  {name, source_sha256, command, seed, provenance, input_sha256,
                  provenance_detail?}            (optional; omitted when None)
      quantities:{key: {value, units, region_id|null, frame|null, time|null,
                        uncertainty|null, provenance, material_region_id?}}
      validity:  {regime, valid_from_utc, valid_until_utc|null,
                  stale_if_input_sha256_mismatch}
      gates:     {...}
      overall_pass: bool

The envelope is intentionally *strict*: a synthetic/unprovenanced cached result is refused by
a default scientific consumer and only accepted when the consumer explicitly opts in.
There is NO universal mass-fraction bound here: a plausibility check may only add a
``flag`` (a physically impossible ``local_mass_kg > total_mass_kg`` is the sole failure).

Five further rules:

* **Abstention is not acceptance.** A quantity that is *present* but whose ``value is None``
  is an explicit abstention (an unknown number), not a value. A consumer must not accept it:
  a required quantity fails with a stable ``abstain:`` substring; an optional one is flagged
  only. This holds in ``strict`` *and* ``synthetic_demo`` mode.
* **Structured provenance detail.** An optional ``producer.provenance_detail`` block separates
  ``data`` / ``calibration`` / ``model`` / ``execution`` instead of a single provenance label.
  When present it is always validated (any mode); a consumer may additionally require it via
  ``validate_for_consumer(..., require_provenance_detail=True)``.
* **Uncertainty requirement.** A scientific claim must state its uncertainty: the optional
  ``require_uncertainty=True`` refuses any *required* quantity that is present but declares no
  uncertainty (``uncertainty is None``). A value of ``0.0`` is a valid declared uncertainty.
* **Producer gate outcome.** A producer publishes its own ``gates`` / ``overall_pass``
  through :func:`make_envelope` (defaults ``{}`` / ``True`` keep historical callers byte-identical).
  When ``overall_pass is not True`` a ``strict`` (scientific) consumer ALWAYS refuses the envelope
  (stable substring ``producer gates failed``); strict mode cannot waive this. In the explicit
  ``synthetic_demo`` mode it is a non-blocking flag unless the caller passes
  ``require_producer_gates_pass=True``. Any individual falsy gate is additionally listed by name.
  :data:`SCIENTIFIC_REQUIREMENTS` additionally requires uncertainty and a provenance detail block.
* **Material-region binding.** A quantity may carry an optional ``material_region_id``
  (an explicit region -> material binding). ``make_envelope`` copies it only when the quantity spec
  provides it, so historical callers that do not supply one keep the exact quantity shape (the key
  is omitted, never added as ``None``). ``validate_for_consumer(..., expected_material_regions=...)``
  accepts a MAPPING ``{quantity_key: material_id}`` (exact binding: a wrong value is
  ``material region mismatch``, an absent one is ``missing material region``) or a SET/frozenset
  vocabulary (any present id outside it is ``unknown material region``, applied to ALL quantities
  including ones outside ``expected_quantities``). ``None`` is a no-op. An optional
  ``material_region_vocabulary`` applies the vocabulary independently of the mapping, so a consumer
  can enforce the exact binding *and* the registered vocabulary in one pass.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

ENVELOPE_VERSION = 1

# --------------------------------------------------------------------------- vocabulary ----
PROVENANCE_VOCABULARY = (
    "measured",
    "calibrated",
    "calibrated_public_model",
    "literal_cited",
    "tuned",
    "synthetic",
    "unknown",
)

# Ad-hoc / alias producer labels -> canonical vocabulary term. Anything unrecognised (and any
# non-string / absent label) normalises to ``"unknown"``, which a default consumer refuses.
PROVENANCE_ALIASES = {
    "measured": "measured",
    "calibrated": "calibrated",
    "calibrated_public_model": "calibrated_public_model",
    "calibrated-public-model": "calibrated_public_model",
    "calibrated_from_public_model": "calibrated_public_model",
    "public_model": "calibrated_public_model",
    # A citation of a repository literal is its own category: it records where a number came
    # from, not that the number was calibrated or measured. It never grants scientific
    # acceptance by itself (it is opt-in only, like ``tuned``).
    "literal_cited": "literal_cited",
    "tuned": "tuned",
    "synthetic": "synthetic",
    "assumption": "unknown",
    # Provisional mapping of the region-payload vocabulary: the payload validator lists
    # ``computed`` among its provenance terms (e.g. a derived mass fraction). A derived value
    # whose inputs are not declared must not be silently promoted to a scientific term, so
    # ``computed`` maps to ``unknown`` (opt-in only). This mapping is provisional until both
    # vocabularies are frozen together.
    "computed": "unknown",
    "UNKNOWN": "unknown",
    "unknown": "unknown",
}

# A scientific consumer trusts these by default. literal_cited / tuned / synthetic / unknown are
# never implicit, and an opted-in non-scientific term never yields ``ACCEPT_SCIENTIFIC``.
DEFAULT_SCIENTIFIC_ALLOWED = frozenset({"measured", "calibrated", "calibrated_public_model"})
OPT_IN_REQUIRED = frozenset({"literal_cited", "tuned", "synthetic", "unknown"})

# --------------------------------------------------------------------------- regime --------
# Regimes the contract knows about. A consumer accepts only the regimes it declares via
# ``validate_for_consumer(expected_regimes=...)``; ``unknown`` is the sentinel for an
# absent / non-string / nonsense label (never implicit).
REGIME_VOCABULARY = (
    "resting",
    "level_walking_steady_state",
    "load_carrying_steady_state",
    "synthetic_demo",
    "unknown",
)

# Ad-hoc / alias regime labels -> canonical vocabulary term. Anything unrecognised (and any
# non-string / absent label) normalises to ``"unknown"``, which a consumer that declares
# ``expected_regimes`` refuses. A *named* out-of-domain exercise regime (e.g. an anaerobic
# sprint) canonicalises to itself: the contract knows it exists but it is never in a scientific
# consumer's allowed set, so the refusal is the more informative ``regime mismatch`` rather than
# ``unknown regime``.
REGIME_ALIASES = {
    "resting": "resting",
    "rest": "resting",
    "basal": "resting",
    "level_walking_steady_state": "level_walking_steady_state",
    "level-walking-steady-state": "level_walking_steady_state",
    "level_walking": "level_walking_steady_state",
    "steady_state_walking": "level_walking_steady_state",
    "load_carrying_steady_state": "load_carrying_steady_state",
    "load_carrying": "load_carrying_steady_state",
    "synthetic_demo": "synthetic_demo",
    "synthetic": "synthetic_demo",
    "unknown": "unknown",
    # Recognised out-of-domain / non-steady-state regimes (known to exist, never accepted here).
    "sprint_anaerobic": "sprint_anaerobic",
}

# A scalar has no geometry frame; a non-scalar frame is canonicalised via sorted JSON.
SCALAR_FRAME = "scalar"

# Sentinel distinguishing an OMITTED expected ``frame``/``time`` key from an explicit
# ``None``. A caller that omits the key asks for NO check; a caller that writes ``"frame": None``
# asks for a NULL frame. A module-private ``object()`` is unambiguous (it can never equal a JSON
# value a spec could carry) and keeps the default lookup a single ``spec.get(..., _ABSENT)``.
_ABSENT = object()

RECEIPT_SCHEMA = "result_receipt_v1"

# The only fields a structured ``producer.provenance_detail`` block may carry; ``data`` is
# required inside the block, the rest are optional (see ``_validate_provenance_detail``).
PROVENANCE_DETAIL_FIELDS = ("data", "calibration", "model", "execution")

# Strict scientific profile: the full evidential bar a *scientific* claim must clear. A caller
# opts in explicitly by spreading it into ``validate_for_consumer``::
#
#     validate_for_consumer(envelope, ..., **SCIENTIFIC_REQUIREMENTS)
#
# Each key is a keyword argument of :func:`validate_for_consumer`. Uncertainty and the provenance
# detail block are not required by default (the historical consumer controls carry neither); a
# declared uncertainty of ``0.0`` is valid — only an absent / ``None`` uncertainty is "missing".
# The producer-gate requirement is already the strict-mode default; it is listed here for
# explicitness.
SCIENTIFIC_REQUIREMENTS = {
    "require_uncertainty": True,
    "require_provenance_detail": True,
    "require_producer_gates_pass": True,
}
"""Strict scientific profile (uncertainty + provenance detail + producer gates). Uncertainty and
provenance detail are opt-in via ``**``; producer gates are enforced in strict mode regardless."""


# --------------------------------------------------------------------------- provenance -----
def normalize_provenance(label) -> str:
    """Map an ad-hoc provenance label onto the canonical ``PROVENANCE_VOCABULARY``.

    ``None``, non-strings and unrecognised labels normalise to ``"unknown"`` (which a default
    scientific consumer refuses rather than silently trusting).
    """
    if not isinstance(label, str):
        return "unknown"
    return PROVENANCE_ALIASES.get(label.strip(), PROVENANCE_ALIASES.get(label.strip().lower(), "unknown"))


def normalize_regime(label) -> str:
    """Map an ad-hoc regime label onto the canonical ``REGIME_VOCABULARY``.

    ``None``, non-strings and unrecognised nonsense labels normalise to ``"unknown"``; whitespace
    and case are tolerated. A recognised but out-of-domain regime (e.g. ``sprint_anaerobic``)
    canonicalises to itself and is refused as a ``regime mismatch`` by a consumer that declares
    ``expected_regimes``.
    """
    if not isinstance(label, str):
        return "unknown"
    stripped = label.strip()
    return REGIME_ALIASES.get(stripped, REGIME_ALIASES.get(stripped.lower(), "unknown"))


def allowed_provenance_for(*, opt_in=(), scientific=True) -> frozenset:
    """Allowed provenance set for a consumer.

    * ``scientific=True`` (default): ``{measured, calibrated, calibrated_public_model}``.
    * ``opt_in`` adds canonical terms the caller explicitly accepts (e.g. ``{"synthetic"}``);
      ``tuned`` / ``synthetic`` / ``unknown`` are never implicit.
    * ``scientific=False``: the whole declared vocabulary (use only for tooling).
    """
    base = DEFAULT_SCIENTIFIC_ALLOWED if scientific else frozenset(PROVENANCE_VOCABULARY)
    return frozenset(base) | frozenset(normalize_provenance(p) for p in opt_in)


# --------------------------------------------------------------------------- helpers -------
def get_dotted(mapping, dotted):
    """Return ``mapping[a][b][c]`` for ``"a.b.c"``, raising ``KeyError(dotted)`` when absent."""
    current = mapping
    for part in dotted.split("."):
        if not isinstance(current, dict) or part not in current:
            raise KeyError(dotted)
        current = current[part]
    return current


def _canonical_frame(value):
    """Canonicalise a frame for comparison: dict -> sorted JSON, None -> None, str -> str."""
    if value is None:
        return None
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True)
    return json.dumps(value, sort_keys=True)


def _prov_label(raw, normalized):
    """Human-readable provenance for failure messages (shows the normalisation)."""
    return repr(raw) if raw == normalized else f"{raw!r} (normalized {normalized!r})"


def _is_sha256_hex(value) -> bool:
    """True iff ``value`` is a 64-character hex string (case-insensitive)."""
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(c in "0123456789abcdefABCDEF" for c in value)
    )


def _validate_provenance_detail(detail, producer, failures):
    """Validate a structured ``producer.provenance_detail`` block into ``failures``.

    Pure and additive: only called when the block is present. Every failure string starts with
    the stable substring ``"provenance detail:"``.
    """
    if not isinstance(detail, dict):
        failures.append(
            f"provenance detail: block must be an object (got {type(detail).__name__})"
        )
        return

    for key in detail:
        if key not in PROVENANCE_DETAIL_FIELDS:
            failures.append(f"provenance detail: unknown field {key!r}")

    # data: required inside the block, must be a declared (non-unknown) provenance term.
    if "data" not in detail:
        failures.append("provenance detail: data is required")
    else:
        raw_data = detail.get("data")
        norm_data = normalize_provenance(raw_data)
        if norm_data not in PROVENANCE_VOCABULARY or norm_data == "unknown":
            failures.append(
                f"provenance detail: data {raw_data!r} is not a declared provenance"
            )

    # calibration: optional; None or a declared (non-unknown) provenance term.
    if "calibration" in detail and detail.get("calibration") is not None:
        raw_cal = detail.get("calibration")
        norm_cal = normalize_provenance(raw_cal)
        if norm_cal not in PROVENANCE_VOCABULARY or norm_cal == "unknown":
            failures.append(
                f"provenance detail: calibration {raw_cal!r} is not a declared provenance"
            )

    # model: optional; a nonempty string, or an object with a nonempty ``name`` string and an
    # optional 64-hex ``sha256``.
    if "model" in detail:
        raw_model = detail.get("model")
        if isinstance(raw_model, str):
            if not raw_model.strip():
                failures.append("provenance detail: model name is empty")
        elif isinstance(raw_model, dict):
            name = raw_model.get("name")
            if not isinstance(name, str) or not name.strip():
                failures.append("provenance detail: model.name must be a nonempty string")
            sha = raw_model.get("sha256")
            if sha is not None and not _is_sha256_hex(sha):
                failures.append(
                    f"provenance detail: model.sha256 {sha!r} is not 64 hex chars"
                )
        else:
            failures.append(
                f"provenance detail: model {raw_model!r} must be a nonempty string "
                "or an object with a name"
            )

    # execution: optional; when present it must carry ``code_sha256`` and ``command`` (nonempty
    # strings) and, when the producer pins a ``source_sha256``, they must agree.
    if "execution" in detail:
        raw_exec = detail.get("execution")
        if not isinstance(raw_exec, dict):
            failures.append(f"provenance detail: execution {raw_exec!r} must be an object")
        else:
            code_sha = raw_exec.get("code_sha256")
            command = raw_exec.get("command")
            if not isinstance(code_sha, str) or not code_sha.strip():
                failures.append(
                    "provenance detail: execution.code_sha256 must be a nonempty string"
                )
            elif not _is_sha256_hex(code_sha):
                failures.append(
                    f"provenance detail: execution.code_sha256 {code_sha!r} is not 64 hex chars"
                )
            else:
                source_sha = producer.get("source_sha256")
                if source_sha is not None and code_sha != source_sha:
                    failures.append(
                        f"provenance detail: execution.code_sha256 {code_sha!r} != "
                        f"producer.source_sha256 {source_sha!r}"
                    )
            if not isinstance(command, str) or not command.strip():
                failures.append(
                    "provenance detail: execution.command must be a nonempty string"
                )


# --------------------------------------------------------------------------- envelope ------
def make_envelope(
    result_obj,
    *,
    producer_name,
    producer_sha256,
    provenance,
    input_sha256,
    quantity_specs,
    command=None,
    seed=None,
    regime=None,
    valid_from_utc=None,
    valid_until_utc=None,
    provenance_detail=None,
    gates=None,
    overall_pass=None,
) -> dict:
    """Wrap a raw result dict into a ResultEnvelope.

    ``quantity_specs`` maps a dotted key to ``{"units", "region_id", "frame", "time",
    "uncertainty", "provenance"}`` and an optional ``"material_region_id"``; only keys present in
    ``result_obj`` become quantities. A missing ``frame`` defaults to :data:`SCALAR_FRAME`; missing
    ``time``/``region_id`` become ``None``; a missing per-quantity ``provenance`` falls back to the
    producer provenance. The optional ``material_region_id`` is copied verbatim when the spec
    provides the key and the key is **omitted** otherwise, so callers that supply no
    binding stay byte-identical to the historical quantity shape.

    ``provenance_detail`` is an optional structured block (``data`` / ``calibration`` /
    ``model`` / ``execution``); when supplied it is stored as
    ``producer["provenance_detail"]``. When ``None`` the key is omitted entirely, so the
    historical producer shape is unchanged.

    ``gates`` / ``overall_pass`` carry the **producer's own gate outcome** into the envelope
   , so a producer that reports ``overall_pass=false`` is no longer silently reported
    as a success. ``gates`` defaults to ``{}`` and ``overall_pass`` defaults to ``True`` when
    ``None``, so a caller that supplies neither gets the historical envelope **byte-identically**
    (``"gates": {}``, ``"overall_pass": true``). The values are stored verbatim; the validator
    interprets a non-``True`` ``overall_pass`` as a flag (default) or a failure when the consumer
    opts into ``require_producer_gates_pass=True``.
    """
    quantities = {}
    for key, spec in quantity_specs.items():
        try:
            raw = get_dotted(result_obj, key)
        except KeyError:
            continue
        if isinstance(raw, dict) and "value" in raw:
            value = raw.get("value")
        else:
            value = raw
        frame = spec.get("frame")
        quantity = {
            "value": value,
            "units": spec.get("units"),
            "region_id": spec.get("region_id"),
            "frame": SCALAR_FRAME if frame is None else frame,
            "time": spec.get("time"),
            "uncertainty": spec.get("uncertainty"),
            "provenance": spec.get("provenance", provenance),
        }
        # Copy the optional per-quantity ``material_region_id`` binding only when the
        # spec actually provides the key. The key is OMITTED otherwise, so a caller that supplies
        # no binding keeps the historical seven-key quantity shape byte-identically (adding
        # ``material_region_id: None`` unconditionally would break the exact-shape contract).
        if "material_region_id" in spec:
            quantity["material_region_id"] = spec.get("material_region_id")
        quantities[key] = quantity

    producer = {
        "name": producer_name,
        "source_sha256": producer_sha256,
        "command": command,
        "seed": seed,
        "provenance": provenance,
        "input_sha256": input_sha256,
    }
    if provenance_detail is not None:
        producer["provenance_detail"] = provenance_detail

    # Carry the producer's own gate outcome. ``gates`` stays the historical empty dict
    # and ``overall_pass`` the historical ``True`` when the caller supplies neither, so existing
    # callers are byte-identical; a producer that passes its real ``overall_pass=false`` (or its
    # real gate dict) through here no longer has that silently overwritten with success.
    if gates is None:
        gates = {}
    if overall_pass is None:
        overall_pass = True

    return {
        "envelope_version": ENVELOPE_VERSION,
        "producer": producer,
        "quantities": quantities,
        "validity": {
            "regime": regime,
            "valid_from_utc": valid_from_utc,
            "valid_until_utc": valid_until_utc,
            "stale_if_input_sha256_mismatch": True,
        },
        "gates": gates,
        "overall_pass": overall_pass,
    }


# --------------------------------------------------------------------------- verdict -------
@dataclass(frozen=True)
class Verdict:
    """Consumer verdict over an envelope.

    ``state`` is one of ``ACCEPT_SCIENTIFIC`` / ``ACCEPT_SYNTHETIC_DEMO`` / ``REFUSE``.

    Scope of ``ACCEPT_SCIENTIFIC``: the envelope passed this consumer's contract checks (units,
    region, frame, time, regime label, provenance label, staleness, producer pin), its producer's
    own gates passed, it is the artifact a trusted registry pinned, and it declares a
    measured/calibrated provenance. It is not a biological or physiological validation.
    ``ACCEPT_SYNTHETIC_DEMO`` is any accepted non-scientific run (explicit demo mode, or an opted-in
    non-scientific provenance such as ``literal_cited``).
    """

    state: str
    failures: tuple = ()
    flags: tuple = ()
    evidence_class: str = "refused"

    @property
    def accepted(self) -> bool:
        return self.state in ("ACCEPT_SCIENTIFIC", "ACCEPT_SYNTHETIC_DEMO")


# --------------------------------------------------------------------------- validate ------
def validate_for_consumer(
    envelope,
    *,
    expected_quantities,
    allowed_provenance,
    mode="strict",
    current_input_sha256=None,
    expected_producer_sha256=None,
    expected_producer_name=None,
    expected_regions=None,
    expected_material_regions=None,
    material_region_vocabulary=None,
    expected_regimes=None,
    plausibility=None,
    require_provenance_detail=False,
    require_uncertainty=False,
    require_producer_gates_pass=None,
    expected_payload_sha256=None,
) -> Verdict:
    """Validate an envelope for one consumer. Pure: no I/O, no framework.

    ``strict`` mode is scientific acceptance. It refuses an envelope whose producer gates did not
    pass (``require_producer_gates_pass`` may not be ``False`` there; ``ValueError``) and one that
    is not the artifact an independently trusted registry approved: ``expected_payload_sha256``
    must be given and equal :func:`payload_sha256` of the envelope (``no trusted artifact pin`` /
    ``artifact not approved``). ``synthetic_demo`` mode never returns ``ACCEPT_SCIENTIFIC``; there
    failed producer gates are a flag unless ``require_producer_gates_pass=True``.

    Failure strings are stable substrings (asserted in tests/test_result_envelope_v1.py). Plausibility checks are
    NON-BLOCKING: they only append to ``flags``; the single exception is a physical
    impossibility (``local_mass_kg > total_mass_kg``), which is a failure.

    ``require_provenance_detail=True`` additionally refuses an envelope whose producer carries no
    structured ``provenance_detail``; a present block is always validated (any mode).

    ``require_uncertainty=True`` additionally refuses (substring ``missing uncertainty``) every
    *required* quantity that is present but declares ``uncertainty is None`` (an absent key and an
    explicit ``None`` are both missing; ``0.0`` is a valid declared uncertainty). Off by default.

    ``require_producer_gates_pass`` governs the producer's own ``overall_pass``. When
    the envelope's ``overall_pass`` is not ``True`` -- including when the key is absent, because
    the validator never invents a passing gate set -- the validator ALWAYS records
    ``"producer gates failed: overall_pass is not true"``: as a FAILURE in ``strict`` mode (always)
    or in ``synthetic_demo`` mode when this is ``True``, otherwise as a non-blocking flag.
    Independently, every falsy entry of the envelope's ``gates`` dict is listed by name in
    ``"producer gates failed: <sorted names>"``: a failure where gates are required and the
    producer claimed ``overall_pass=true`` (inconsistent producer), otherwise a flag next to the
    ``overall_pass`` failure.

    **Absent vs explicit expected ``frame``/``time``.** For every expected quantity a
    spec that OMITS the ``frame`` (or ``time``) key is not checked at all; a spec that explicitly
    carries ``"frame": None`` (or ``"time": None``) demands a null value, and any other explicit
    value is compared canonically. This mirrors the region seam's absent-frame semantics (an
    undeclared frame is never invented) and means omitting a field cannot silently mean "expect
    null". Existing callers always pass ``frame``/``time`` explicitly, so their behaviour is
    unchanged.

    ``expected_material_regions`` validates the optional per-quantity
    ``material_region_id`` binding. A MAPPING ``{quantity_key: material_id}`` demands the exact id
    per key: a different id is ``f"material region mismatch: {key} expected {want!r} got {got!r}"``
    and an absent one is ``f"missing material region: {key}"``. A SET/frozenset is a vocabulary:
    any *present* id outside it fails with ``f"unknown material region: {key} {got!r}"`` and the
    rule is applied to ALL quantities, including ones outside ``expected_quantities`` (mirroring the
    region-vocabulary rule; a quantity that declares no binding is exempt). ``None`` is a no-op.
    ``material_region_vocabulary`` applies a vocabulary independently of the mapping form, so a
    caller can enforce an exact binding and the registered vocabulary in one pass.
    """
    if mode not in ("strict", "synthetic_demo"):
        raise ValueError(f"unknown mode {mode!r}; expected 'strict' or 'synthetic_demo'")
    if mode == "strict":
        if require_producer_gates_pass is False:
            raise ValueError(
                "strict mode cannot waive producer gates: use mode='synthetic_demo' for a "
                "non-scientific run"
            )
        require_producer_gates_pass = True
    else:
        require_producer_gates_pass = bool(require_producer_gates_pass)

    if not isinstance(envelope, dict) or "envelope_version" not in envelope:
        return Verdict("REFUSE", ("malformed envelope",), (), "refused")

    quantities = envelope.get("quantities") or {}
    producer = envelope.get("producer") or {}
    validity = envelope.get("validity") or {}
    allowed = frozenset(normalize_provenance(p) for p in allowed_provenance)

    failures = []
    flags = []

    # 3. required quantity present.
    for key, spec in expected_quantities.items():
        if spec.get("required", True) and key not in quantities:
            failures.append(f"missing quantity: {key}")

    # 3a. abstention is not acceptance. A quantity that is PRESENT but carries ``value is None``
    # is an explicit "unknown number", not a value; it must not be accepted as scientific or as
    # demo evidence. A required quantity is refused with a stable ``abstain:`` substring; an
    # optional one is only flagged. This applies in BOTH modes and is independent of whether the
    # units/region/frame/time are otherwise correct.
    for key, spec in expected_quantities.items():
        if key not in quantities:
            continue
        quantity = quantities.get(key) or {}
        if quantity.get("value") is None:
            if spec.get("required", True):
                failures.append(f"abstain: required value missing for {key}")
            else:
                flags.append(f"abstain: optional value missing for {key}")

    # 3b. uncertainty requirement. A strict scientific consumer can require that every
    # *required* quantity that is present declares an uncertainty. An absent key and an explicit
    # ``None`` are both "missing"; a value of ``0.0`` is a valid declared uncertainty. This is off
    # by default: the healthy thermo/region controls carry no uncertainty and must keep passing.
    if require_uncertainty:
        for key, spec in expected_quantities.items():
            if not spec.get("required", True) or key not in quantities:
                continue
            quantity = quantities.get(key) or {}
            if quantity.get("uncertainty") is None:
                failures.append(f"missing uncertainty: {key}")

    # 4-7. units / region / frame / time for every declared quantity that is present.
    for key, spec in expected_quantities.items():
        if key not in quantities:
            continue
        quantity = quantities.get(key) or {}

        want_units = spec.get("units")
        got_units = quantity.get("units")
        if want_units is not None and got_units != want_units:
            failures.append(f"unit mismatch: {key} expected {want_units!r} got {got_units!r}")

        want_region = spec.get("region_id")
        got_region = quantity.get("region_id")
        if want_region is not None and got_region is None:
            failures.append(f"missing region_id: {key}")
        elif got_region is not None and got_region != want_region:
            failures.append(f"region mismatch: {key} expected {want_region!r} got {got_region!r}")
        # ``expected_regions`` is a region vocabulary: it applies only to quantities that actually
        # declare a region. A quantity whose contract declares ``region_id=None`` (e.g. a
        # ``totals.*`` aggregate) legitimately has no region and must not be reported as an
        # unknown region; a quantity that *requires* a region is already failed above as
        # ``missing region_id`` when it carries ``None``, so this is not an evasion path.
        if (
            expected_regions is not None
            and got_region is not None
            and got_region not in expected_regions
        ):
            failures.append(f"unknown region_id: {key} {got_region!r}")

        # An OMITTED ``frame`` key means "do not check this quantity's frame", while
        # an explicit ``"frame": None`` means "expect a null frame" (the same absent-vs-explicit
        # distinction the region seam uses for an undeclared frame). Before this, a spec without a
        # ``frame`` key was silently treated as ``frame is None`` and demanded a null frame from a
        # quantity the caller never asked about. Existing callers always pass frame/time explicitly,
        # so their behaviour is unchanged.
        want_frame = spec.get("frame", _ABSENT)
        if want_frame is not _ABSENT:
            want_frame = _canonical_frame(want_frame)
            got_frame = _canonical_frame(quantity.get("frame"))
            if want_frame != got_frame:
                failures.append(
                    f"frame mismatch: {key} expected {want_frame!r} got {got_frame!r}"
                )

        # Same absent-vs-explicit rule for ``time``; an omitted key is not checked.
        want_time = spec.get("time", _ABSENT)
        if want_time is not _ABSENT:
            got_time = quantity.get("time")
            if want_time != got_time:
                failures.append(f"time mismatch: {key} expected {want_time!r} got {got_time!r}")

    # 7a. region vocabulary across ALL declared quantities, not only the required contract: a
    # consumer that declares ``expected_regions`` must not accept an envelope that also carries a
    # quantity bound to a region outside that vocabulary (e.g. a region the caller never
    # requested). ``region_id=None`` is a legitimately global quantity (an aggregate) and is
    # exempt, as is a quantity whose key the contract already covers (checked in 4-7 above).
    if expected_regions is not None:
        for key, quantity in quantities.items():
            if key in expected_quantities:
                continue
            got_region = quantity.get("region_id") if isinstance(quantity, dict) else None
            if got_region is not None and got_region not in expected_regions:
                failures.append(f"unknown region_id: {key} {got_region!r}")

    # 7a-material. Material-region binding / vocabulary. A quantity may declare an
    # explicit region -> material binding (``material_region_id``). ``expected_material_regions``
    # is either a MAPPING ``{quantity_key: material_id}`` (exact binding: a present-but-different
    # id is a mismatch, an absent id is missing) or a SET/frozenset vocabulary (any *present* id
    # outside it is unknown -- applied to ALL quantities, including ones outside
    # ``expected_quantities``, mirroring the region-vocabulary rule above). ``None`` is a no-op, so
    # historical callers are unaffected. ``material_region_vocabulary`` applies the vocabulary
    # independently of the mapping form, so a caller can enforce the exact binding AND the
    # registered vocabulary in one pass (``preflight_region_payload`` does exactly that).
    if expected_material_regions is not None:
        if isinstance(expected_material_regions, Mapping):
            for key, want in expected_material_regions.items():
                quantity = quantities.get(key)
                got = quantity.get("material_region_id") if isinstance(quantity, dict) else None
                if got is None:
                    failures.append(f"missing material region: {key}")
                elif got != want:
                    failures.append(
                        f"material region mismatch: {key} expected {want!r} got {got!r}"
                    )
        else:
            vocabulary = frozenset(expected_material_regions)
            for key, quantity in quantities.items():
                got = quantity.get("material_region_id") if isinstance(quantity, dict) else None
                if got is not None and got not in vocabulary:
                    failures.append(f"unknown material region: {key} {got!r}")
    if material_region_vocabulary is not None:
        vocabulary = frozenset(material_region_vocabulary)
        for key, quantity in quantities.items():
            got = quantity.get("material_region_id") if isinstance(quantity, dict) else None
            if got is not None and got not in vocabulary:
                failures.append(f"unknown material region: {key} {got!r}")

    # 7b. regime, only when the consumer declares which regimes it accepts. An absent regime
    # (None) is missing, an unrecognised label is unknown, and a recognised but disallowed regime
    # is a mismatch; none of the three may be accepted as scientific evidence. When
    # ``expected_regimes`` is None the check is a no-op (backwards compatible).
    if expected_regimes is not None:
        raw_regime = validity.get("regime")
        if raw_regime is None:
            failures.append("missing regime: envelope declares no regime")
        else:
            norm_regime = normalize_regime(raw_regime)
            if norm_regime == "unknown":
                failures.append(f"unknown regime: {raw_regime!r}")
            elif norm_regime not in expected_regimes:
                failures.append(
                    f"regime mismatch: {raw_regime!r} (normalized {norm_regime!r}) "
                    f"not in {sorted(expected_regimes)}"
                )

    # 8. producer provenance allowed (absent -> "unknown" -> refuses for a scientific consumer).
    producer_prov_raw = producer.get("provenance")
    producer_prov = normalize_provenance(producer_prov_raw)
    if producer_prov not in allowed:
        failures.append(
            "disallowed provenance: producer "
            f"{_prov_label(producer_prov_raw, producer_prov)} not in {sorted(allowed)}"
        )

    # 9. per-quantity provenance, when present, allowed.
    for key in expected_quantities:
        if key not in quantities:
            continue
        qprov_raw = (quantities.get(key) or {}).get("provenance")
        if qprov_raw is None:
            continue
        qprov = normalize_provenance(qprov_raw)
        if qprov not in allowed:
            failures.append(
                "disallowed provenance: quantity "
                f"{key} {_prov_label(qprov_raw, qprov)} not in {sorted(allowed)}"
            )

    # 9b. structured provenance detail (data / calibration / model / execution). Absent is fine
    # by default; a present block is ALWAYS validated (any mode). A consumer may require it.
    detail = producer.get("provenance_detail")
    if detail is None:
        if require_provenance_detail:
            failures.append("missing provenance_detail")
    else:
        _validate_provenance_detail(detail, producer, failures)

    # 10. producer code hash pinned by the registry.
    if expected_producer_sha256 is not None:
        got_sha = producer.get("source_sha256")
        if got_sha != expected_producer_sha256:
            failures.append(
                f"producer changed: source_sha256 {got_sha!r} != expected {expected_producer_sha256!r}"
            )

    # 11. producer name pinned by the registry.
    if expected_producer_name is not None:
        got_name = producer.get("name")
        if got_name != expected_producer_name:
            failures.append(f"producer name mismatch: {got_name!r} != {expected_producer_name!r}")

    # 12. stale input. A scientific/consumer path must not accept an envelope that opts OUT of
    # staleness checking: a producer could otherwise flip the flag and silently ship a changed
    # input (acceptance rule 1). Declaring the flag false is itself a contract violation.
    if validity.get("stale_if_input_sha256_mismatch") is not True:
        failures.append(
            "stale check disabled: envelope declares stale_if_input_sha256_mismatch "
            f"{validity.get('stale_if_input_sha256_mismatch')!r} (must be true)"
        )
    elif producer.get("input_sha256") != current_input_sha256:
        failures.append(
            "stale input: producer input_sha256 "
            f"{producer.get('input_sha256')!r} != current {current_input_sha256!r}"
        )

    # 12a. approved artifact. Scientific acceptance requires that this exact envelope is the
    # artifact an independently trusted registry pinned; a matching receipt is not enough,
    # because anyone can recompute an unkeyed digest after editing the payload.
    if mode == "strict":
        if expected_payload_sha256 is None:
            failures.append(
                "no trusted artifact pin: expected_payload_sha256 is required for scientific "
                "acceptance"
            )
        else:
            got_payload_sha = payload_sha256(envelope)
            if got_payload_sha != expected_payload_sha256:
                failures.append(
                    f"artifact not approved: payload_sha256 {got_payload_sha!r} != trusted pin "
                    f"{expected_payload_sha256!r}"
                )

    # 12b. producer's own gate outcome. A producer may publish its ``gates`` /
    # ``overall_pass`` via ``make_envelope``; the validator must not silently report success when
    # the producer itself reports failure. ``envelope.get("overall_pass") is not True`` covers
    # both an explicit ``False`` and a missing key (an envelope that never declared a gate outcome
    # is not proven to pass). The message is the same in both roles; only its severity depends on
    # the consumer's explicit opt-in. The per-gate listing is always a non-blocking flag and is
    # independent of ``overall_pass`` (a producer could self-report ``overall_pass=true`` while a
    # named gate is falsy -- that inconsistency must still be visible).
    overall_pass = envelope.get("overall_pass")
    if overall_pass is not True:
        _producer_gates_msg = "producer gates failed: overall_pass is not true"
        if require_producer_gates_pass:
            failures.append(_producer_gates_msg)
        else:
            flags.append(_producer_gates_msg)
    envelope_gates = envelope.get("gates")
    if isinstance(envelope_gates, dict):
        failing_gate_names = sorted(name for name, passed in envelope_gates.items() if not passed)
        if failing_gate_names:
            # A named failing gate is a refusal wherever the gate outcome is required, even when
            # the producer inconsistently reports overall_pass=true.
            (failures if require_producer_gates_pass and overall_pass is True else flags).append(
                "producer gates failed: " + ", ".join(failing_gate_names))

    # 13. plausibility: NON-BLOCKING flags; only physical impossibility is a failure.
    if plausibility:
        for check in plausibility:
            if not isinstance(check, dict):
                continue
            name = check.get("name") or check.get("kind") or "plausibility"
            local = check.get("local_mass_kg")
            total = check.get("total_mass_kg")
            if local is not None and total is not None and float(local) > float(total):
                failures.append(
                    f"physically impossible mass: {name} "
                    f"local_mass_kg={local!r} > total_mass_kg={total!r}"
                )
                continue
            source = check.get("source")
            if source is not None:
                flags.append(f"{name}: sourced bound {source!r} checked (non-blocking)")
            elif check.get("max_fraction") is not None:
                # Only a check that declares a (sourceless) fraction bound is flagged; a bare
                # local/total check stays quiet when local <= total.
                flags.append(
                    f"{name}: unsourced plausibility bound not enforced "
                    "(max_without_source_ok=True; never refuses)"
                )

    if failures:
        return Verdict("REFUSE", tuple(failures), tuple(flags), "refused")

    if mode == "synthetic_demo":
        if producer_prov in DEFAULT_SCIENTIFIC_ALLOWED:
            return Verdict("ACCEPT_SYNTHETIC_DEMO", (), tuple(flags), producer_prov)
        return Verdict("ACCEPT_SYNTHETIC_DEMO", (), tuple(flags), "synthetic_demo")

    # strict: scientific terms accept as scientific; an explicitly opted-in non-scientific term
    # (literal_cited / tuned / synthetic / unknown) is accepted only as a non-scientific run.
    if producer_prov in DEFAULT_SCIENTIFIC_ALLOWED:
        return Verdict("ACCEPT_SCIENTIFIC", (), tuple(flags), producer_prov)
    if producer_prov in allowed:
        return Verdict("ACCEPT_SYNTHETIC_DEMO", (), tuple(flags), "synthetic_demo")
    # Unreachable: rule 8 already refused a provenance outside ``allowed``.
    return Verdict(
        "REFUSE",
        (f"disallowed provenance: producer {_prov_label(producer_prov_raw, producer_prov)}",),
        tuple(flags),
        "refused",
    )


# --------------------------------------------------------------------------- receipt -------
def _sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha256(path) -> str:
    """sha256 of a file's bytes (the producer / payload artifact identity)."""
    return _sha256_hex(Path(path).read_bytes())


def _payload_bytes(payload) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")


def payload_sha256(payload) -> str:
    """Canonical sha256 of a JSON payload (sorted keys, compact); the value receipts and registry
    pins use. Integrity identity only -- it says which bytes, not who approved them."""
    return _sha256_hex(_payload_bytes(payload))


def _upstream_link(upstream, allow_failed_upstream=False) -> dict:
    """Build the ``receipt["upstream"]`` link from a sealed artifact or a receipt.

    When the sealed upstream payload is a ResultEnvelope, its ``overall_pass`` is recorded as
    ``upstream_overall_pass``; an upstream whose own gates did not pass is refused
    (``producer gates failed``) unless ``allow_failed_upstream=True`` (explicit non-scientific use,
    which :func:`read_verified_chain` still refuses).

    * Full sealed artifact ``{"payload", "receipt"}``: the upstream ``payload_sha256`` is
      RECOMPUTED from ``upstream["payload"]`` (a caller-supplied hash is not trusted). The
      upstream receipt itself is not verified here; a mutated upstream payload therefore records
      its new hash, which a verifier holding the original expected hash refuses as stale.
    * Receipt dict: its ``payload_sha256`` is copied as-is; the caller is responsible for having
      verified that receipt (``verify_receipt``) before chaining on it.
    """
    is_envelope = False
    upstream_pass = None
    if isinstance(upstream, Mapping) and "payload" in upstream and "receipt" in upstream:
        receipt = upstream["receipt"]
        if not isinstance(receipt, Mapping):
            raise ValueError("invalid upstream receipt: 'receipt' is not a dict")
        payload_sha = _sha256_hex(_payload_bytes(upstream["payload"]))
        upstream_payload = upstream["payload"]
        if isinstance(upstream_payload, Mapping) and "envelope_version" in upstream_payload:
            is_envelope = True
            upstream_pass = upstream_payload.get("overall_pass")
            upstream_gates = upstream_payload.get("gates")
            if isinstance(upstream_gates, Mapping) and any(
                not passed for passed in upstream_gates.values()
            ):
                # A contradictory aggregate flag cannot override a failed individual gate.
                # Preserve the effective failure in the link even for an explicit demo seal.
                upstream_pass = False
            if upstream_pass is not True and not allow_failed_upstream:
                raise ValueError(
                    "producer gates failed: upstream envelope overall_pass is not true or a named gate failed; a chained "
                    "result cannot build on it"
                )
    elif isinstance(upstream, Mapping):
        receipt = upstream
        payload_sha = receipt.get("payload_sha256")
        if payload_sha is None:
            raise ValueError("invalid upstream receipt: missing payload_sha256")
        if not _is_sha256_hex(payload_sha):
            raise ValueError(
                f"invalid upstream receipt: payload_sha256 {payload_sha!r} is not 64 hex chars"
            )
    else:
        raise ValueError(
            f"invalid upstream receipt: expected a sealed artifact or receipt dict, "
            f"got {type(upstream).__name__}"
        )
    link = {
        "name": str(receipt.get("producer_path")),
        "payload_sha256": payload_sha,
        "producer_sha256": receipt.get("producer_sha256"),
        "schema": RECEIPT_SCHEMA,
    }
    if is_envelope:
        link["upstream_overall_pass"] = upstream_pass
    return link


def seal(payload, producer_path, *, upstream=None, allow_failed_upstream=False) -> dict:
    """Return ``payload`` wrapped with a receipt binding payload + producer code version.

    ``upstream`` (optional): a sealed artifact ``{"payload", "receipt"}`` or a receipt
    dict of the upstream consumer this payload was computed from; recorded as
    ``receipt["upstream"]`` (see :func:`_upstream_link`). ``upstream=None`` yields exactly the
    receipt without chaining (no ``upstream`` key).
    """
    receipt = {
        "payload_sha256": _sha256_hex(_payload_bytes(payload)),
        "producer_sha256": file_sha256(producer_path),
        "producer_path": str(producer_path),
        "sealed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "schema": RECEIPT_SCHEMA,
    }
    if upstream is not None:
        receipt["upstream"] = _upstream_link(upstream, allow_failed_upstream=allow_failed_upstream)
    return {"payload": payload, "receipt": receipt}


def seal_chain(payload, producer_path, upstream_sealed, *, allow_failed_upstream=False) -> dict:
    """Seal a downstream payload bound to a verified upstream artifact (``upstream`` required).

    An upstream envelope whose own gates did not pass is refused (``producer gates failed``)
    unless ``allow_failed_upstream=True``."""
    if upstream_sealed is None:
        raise ValueError("invalid upstream receipt: seal_chain requires an upstream artifact")
    return seal(payload, producer_path, upstream=upstream_sealed,
                allow_failed_upstream=allow_failed_upstream)


def verify_receipt(
    sealed,
    *,
    producer_path=None,
    expected_producer_sha256=None,
    expected_upstream_payload_sha256=None,
    require_upstream=False,
):
    """Return the payload, raising ``ValueError`` on a missing/tampered/stale receipt.

    INTEGRITY ONLY. An unkeyed sha256 detects accidental change and a stale payload; it does not
    show who wrote the payload (an edit with a recomputed hash verifies) and it is not evidence
    acceptance. Scientific acceptance is :func:`validate_for_consumer` in strict mode with a
    trusted ``expected_payload_sha256`` pin.

    Chain checks run after the existing payload/producer checks:
    ``require_upstream`` refuses a receipt without ``upstream``; ``expected_upstream_payload_sha256``
    refuses a receipt whose recorded upstream payload hash differs (stale upstream).
    """
    if not isinstance(sealed, dict) or "receipt" not in sealed or "payload" not in sealed:
        raise ValueError("unsealed result: no receipt (treat as stale)")
    receipt = sealed["receipt"]
    if not isinstance(receipt, dict):
        raise ValueError("unsealed result: no receipt (treat as stale)")
    if _sha256_hex(_payload_bytes(sealed["payload"])) != receipt.get("payload_sha256"):
        raise ValueError("stale result: payload hash mismatch")
    if producer_path is not None and file_sha256(producer_path) != receipt.get("producer_sha256"):
        raise ValueError("stale result: producer code changed since seal")
    if expected_producer_sha256 is not None and receipt.get("producer_sha256") != expected_producer_sha256:
        raise ValueError(
            "producer changed: producer_sha256 "
            f"{receipt.get('producer_sha256')!r} != expected {expected_producer_sha256!r}"
        )
    link = receipt.get("upstream")
    if (require_upstream or expected_upstream_payload_sha256 is not None) and link is None:
        raise ValueError("missing upstream receipt (chain not sealed)")
    if link is not None and not isinstance(link, dict):
        raise ValueError("invalid upstream receipt: 'upstream' is not a dict")
    if expected_upstream_payload_sha256 is not None:
        got = link.get("payload_sha256")
        if got != expected_upstream_payload_sha256:
            raise ValueError(
                f"stale upstream: upstream payload_sha256 {got!r} != "
                f"expected {expected_upstream_payload_sha256!r}"
            )
    return sealed["payload"]


def write_sealed(path, payload, producer_path) -> None:
    """Seal ``payload`` and write ``{"payload", "receipt"}`` to ``path`` as indented JSON."""
    Path(path).write_text(
        json.dumps(seal(payload, producer_path), indent=2) + "\n", encoding="utf-8"
    )


def read_verified(path, producer_path=None, expected_producer_sha256=None):
    """Read a sealed file and return its payload after receipt verification."""
    sealed = json.loads(Path(path).read_text(encoding="utf-8"))
    return verify_receipt(
        sealed, producer_path=producer_path, expected_producer_sha256=expected_producer_sha256
    )


def read_verified_chain(
    path,
    producer_path=None,
    *,
    expected_upstream_payload_sha256=None,
    require_upstream=True,
    **kw,
):
    """Consumer-side chain reader: verify the receipt and the upstream link.

    ``expected_upstream_payload_sha256`` is REQUIRED: the caller must supply the upstream identity
    it independently expects (e.g. the payload hash of the upstream it actually read or a registry
    pin). A bare upstream link is refused, and so is a link that records an upstream whose own
    gates did not pass.
    """
    if expected_upstream_payload_sha256 is None:
        raise ValueError(
            "expected upstream identity required: pass expected_upstream_payload_sha256 "
            "(a bare upstream link is not evidence)"
        )
    sealed = json.loads(Path(path).read_text(encoding="utf-8"))
    link = (sealed.get("receipt") or {}).get("upstream") if isinstance(sealed, dict) else None
    if isinstance(link, dict) and "upstream_overall_pass" in link \
            and link["upstream_overall_pass"] is not True:
        raise ValueError("producer gates failed: the chained upstream did not pass its own gates")
    return verify_receipt(
        sealed,
        producer_path=producer_path,
        expected_upstream_payload_sha256=expected_upstream_payload_sha256,
        require_upstream=require_upstream,
        **kw,
    )
