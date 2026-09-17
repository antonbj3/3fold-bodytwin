"""Consumer preflight for cached metabolic-cost results (thermoregulation and siblings).

This is the ONE wired consumer of the :mod:`bodytwin.framework.result_envelope_v1` contract.
It decides, before any numbers are used, whether a cached ``metabolic_cost_results.json`` is:

* a sealed, registry-pinned, calibrated-public-model result -> scientific evidence;
* an explicitly declared ``synthetic: true`` stand-in -> a labelled synthetic demo only;
* anything else (unsealed envelope, raw unknown, wrong units/region/frame/time, wrong/missing/
  unknown regime, changed producer, stale input, absent provenance) -> REFUSE.

The regime gate accepts only :data:`THERMO_ACCEPTED_REGIMES` for scientific evidence (plus
``synthetic_demo`` on the explicit demo path); it is enforced by the shared
:func:`validate_for_consumer(expected_regimes=...)` check.

Region vocabulary is the **roster-registered** set, mirrored by
:mod:`bodytwin.framework.region_vocabulary_v1` from the frozen ``bodytwin.geometry.material_roster_v1``
(sha256 ``1998eec14cd3195f35705fb3d08324a765b02e23e9fde31a0fabbfb69e08ae2a``):
``CARD-VESSEL``, ``GENERIC``, ``MSK-BONE``, ``MSK-MUSCLE``, ``MSK-TENDON``, ``NEU-TISSUE``,
``WHOLE-BODY``. This module does NOT modify the region payload schema. ``THERMO_DECLARED_SYNTHETIC_REGION``
(``MSK-LOWERLIMB-HIP``) is kept only as a named placeholder the roster never registered, so a
payload declaring it is refused as ``unknown region_id``.

An ``ACCEPT_SCIENTIFIC`` verdict is a contract result, not a biological validation. It requires
that the producer's own gates passed and that the sealed payload is the registry-pinned artifact.
The adapted public gait2392 control fixture, whose producer gates report
``overall_pass=false``, is refused (``producer gates failed``). A synthetic stand-in is only ever
``ACCEPT_SYNTHETIC_DEMO``. See :class:`bodytwin.framework.result_envelope_v1.Verdict`.

Stdlib only; no network; the only I/O is reading a registry JSON supplied by the caller.
"""
from __future__ import annotations

import json
from pathlib import Path

from bodytwin.framework.region_vocabulary_v1 import I2_REGION_IDS
from bodytwin.framework.result_envelope_v1 import (
    SCALAR_FRAME,
    Verdict,
    allowed_provenance_for,
    get_dotted,
    make_envelope,
    validate_for_consumer,
    verify_receipt,
)

# --------------------------------------------------------------------------- vocabulary ----
# The vocabulary is exactly the region ids the roster registers (single source mirror, see
# region_vocabulary_v1); it is deliberately NOT {WHOLE-BODY, MSK-MUSCLE, MSK-LOWERLIMB-HIP}.
THERMO_REGION_VOCABULARY = frozenset(I2_REGION_IDS)

# Named placeholder, deliberately NOT in THERMO_REGION_VOCABULARY: the roster never registers
# ``MSK-LOWERLIMB-HIP``, so a payload declaring it is refused as ``unknown region_id``.
THERMO_DECLARED_SYNTHETIC_REGION = "MSK-LOWERLIMB-HIP"

# Derived through the shared vocabulary helper (never a raw string comparison).
THERMO_ALLOWED_PROVENANCE = allowed_provenance_for() | frozenset({"calibrated_public_model"})

# Scientific thermoregulation evidence is only meaningful for a steady-state walking regime.
# The explicit synthetic-demo path additionally declares ``synthetic_demo`` (see preflight).
THERMO_ACCEPTED_REGIMES = frozenset({"level_walking_steady_state"})

# Exact per-quantity contract for the nine quantities thermoregulation.py reads.
THERMO_QUANTITY_SPECS = {
    "muscle_mass.total_body_mass_kg": {
        "units": "kg", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": None, "required": True,
    },
    "muscle_mass.summed_muscle_mass_kg": {
        "units": "kg", "region_id": "MSK-MUSCLE", "frame": SCALAR_FRAME,
        "time": None, "required": True,
    },
    "headline.umberger2010_gross_w_per_kg": {
        "units": "W/kg", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": "steady_state", "required": True,
    },
    "headline.umberger2010_net_w_per_kg": {
        "units": "W/kg", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": "steady_state", "required": True,
    },
    "headline.bhargava2004_gross_w_per_kg": {
        "units": "W/kg", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": "steady_state", "required": True,
    },
    "headline.bhargava2004_net_w_per_kg": {
        "units": "W/kg", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": "steady_state", "required": True,
    },
    "sensitivity.combined_tendon_and_mass_correction_w_per_kg": {
        "units": "W/kg", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": "steady_state", "required": True,
    },
    "sensitivity.muscle_mass_correction_scale_applied": {
        "units": "dimensionless", "region_id": "MSK-MUSCLE", "frame": SCALAR_FRAME,
        "time": None, "required": True,
    },
    "distance_speed.speed_mps": {
        "units": "m/s", "region_id": "WHOLE-BODY", "frame": SCALAR_FRAME,
        "time": "steady_state", "required": True,
    },
}

# The nine metabolic-cost quantities are shared verbatim between the thermoregulation
# and cardiac_output consumers (same producer result file, same units/regions/frames/times). The
# spec set is therefore shared; this documented alias names it from the second consumer's side so
# neither consumer has to reach for a thermoregulation-branded constant.
METABOLIC_COST_QUANTITY_SPECS = THERMO_QUANTITY_SPECS

# Input classification tags (also the ``kind`` argument of :func:`resolve_mode`).
SEALED = "SEALED"
ENVELOPE_UNSEALED = "ENVELOPE_UNSEALED"
RAW_SYNTHETIC = "RAW_SYNTHETIC"
RAW_UNKNOWN = "RAW_UNKNOWN"

_UNSEALED_FAILURE = "unsealed result: no receipt (treat as stale)"
_RAW_UNKNOWN_FAILURE = "unknown input: no provenance"
_SEALED_NOT_ENVELOPE = "malformed envelope: sealed payload is not a ResultEnvelope"


# --------------------------------------------------------------------------- registry ------
def load_registry(path) -> dict:
    """Load a registry JSON file (trusted pins: producer, input and approved payload hashes).

    Admission and update rules: the registry is read-only at run time; no code in this package
    writes it. A pin is admitted or changed only by an explicit, reviewed change of the registry
    file (for the test registry: ``tests/fixtures/i1/build_fixtures.py`` writes the approved
    ``payload_sha256`` and the change is committed). A consumer's scientific acceptance requires
    the sealed payload to equal the registry's ``payload_sha256``; a recomputed receipt for an
    edited payload does not match that pin and is refused.
    """
    return json.loads(Path(path).read_text(encoding="utf-8"))


def registry_entry(registry, result_name):
    """Return the registry entry for ``result_name`` (or ``None``)."""
    if not isinstance(registry, dict):
        return None
    return (registry.get("entries") or {}).get(result_name)


# --------------------------------------------------------------------------- classify ------
def classify_input(obj) -> str:
    """Classify a loaded JSON object into one of the four input kinds."""
    if isinstance(obj, dict) and "payload" in obj and "receipt" in obj:
        return SEALED
    if isinstance(obj, dict) and "envelope_version" in obj:
        return ENVELOPE_UNSEALED
    if isinstance(obj, dict) and obj.get("synthetic") is True:
        return RAW_SYNTHETIC
    return RAW_UNKNOWN


def resolve_mode(requested, kind) -> str:
    """Resolve the effective mode: an explicit mode is honoured, ``auto`` is classified.

    ``auto``/``None`` -> ``synthetic_demo`` for an explicit synthetic payload, else ``strict``.
    """
    if requested in ("strict", "synthetic_demo"):
        return requested
    if requested is None or (isinstance(requested, str) and requested.strip().lower() in ("auto", "")):
        return "synthetic_demo" if kind == RAW_SYNTHETIC else "strict"
    raise ValueError(f"unknown mode {requested!r}; expected 'strict', 'synthetic_demo' or 'auto'")


# --------------------------------------------------------------------------- helpers -------
def _unnest_quantities(quantities):
    """Rebuild the raw nested result dict a consumer reads from flat dotted quantities."""
    out = {}
    for dotted, quantity in (quantities or {}).items():
        value = quantity.get("value") if isinstance(quantity, dict) else quantity
        cursor = out
        parts = dotted.split(".")
        for part in parts[:-1]:
            nxt = cursor.get(part)
            if not isinstance(nxt, dict):
                nxt = {}
                cursor[part] = nxt
            cursor = nxt
        cursor[parts[-1]] = value
    return out


def _envelope_and_result(payload):
    """Return ``(envelope, raw_result)`` for a verified sealed payload.

    The sealed ``payload`` MUST itself be a ResultEnvelope: a dict with ``envelope_version`` and
    a ``quantities`` dict. The consumed raw result is reconstructed ONLY from those validated
    quantities (``_unnest_quantities``), so validation and consumption are literally the same
    object. There is no fallback: a sibling top-level ``envelope`` key is ignored, and an
    arbitrary raw payload is never re-wrapped with registry metadata (both were soundness holes
    that let a validated envelope coexist with, and legitimise, a different consumed payload).

    Returns ``(None, None)`` when the payload is not a ResultEnvelope; the caller turns that
    into a REFUSE (``_SEALED_NOT_ENVELOPE``).
    """
    if (
        isinstance(payload, dict)
        and "envelope_version" in payload
        and isinstance(payload.get("quantities"), dict)
    ):
        return payload, _unnest_quantities(payload["quantities"])
    return None, None


def thermo_plausibility(result):
    """Physical-impossibility check for the local working-muscle mass.

    Returns a one-element plausibility list when the three mass keys are present, else ``[]``.
    The check declares only local/total masses (no ``max_fraction`` and no ``source``), so the
    validator stays quiet for a possible-but-unusual fraction and fails only when the local
    working-muscle mass exceeds the total body mass.
    """
    try:
        summed = get_dotted(result, "muscle_mass.summed_muscle_mass_kg")
        scale = get_dotted(result, "sensitivity.muscle_mass_correction_scale_applied")
        total = get_dotted(result, "muscle_mass.total_body_mass_kg")
    except (KeyError, TypeError):
        return []
    try:
        local = float(summed) * float(scale)
        total_kg = float(total)
    except (TypeError, ValueError):
        return []
    return [
        {
            "name": "local_working_muscle_mass",
            "kind": "region_mass_fraction",
            "local_mass_kg": local,
            "total_mass_kg": total_kg,
        }
    ]


# --------------------------------------------------------------------------- preflight -----
def preflight_metabolic_cost(
    raw_obj,
    *,
    mode="auto",
    registry_entry=None,
    current_input_sha256=None,
    producer_path=None,
    provenance_detail=None,
    require_provenance_detail=False,
    require_uncertainty=False,
    require_producer_gates_pass=None,
):
    """Consumer-agnostic preflight of a ``metabolic_cost_results.json`` input.

    This is the shared implementation extracted from the original
    thermoregulation-specific preflight. It decides, before any number is used, whether a cached
    ``metabolic_cost_results.json`` is a sealed/registry-pinned/scientific result, an explicitly
    declared ``synthetic: true`` stand-in (a labelled demo only), or a REFUSE. Every
    metabolic-cost consumer reads the same nine quantities, so the quantity contract
    (:data:`METABOLIC_COST_QUANTITY_SPECS`), region vocabulary
    (:data:`THERMO_REGION_VOCABULARY`) and accepted regimes (:data:`THERMO_ACCEPTED_REGIMES`)
    are shared; the constant names retain their historical ``THERMO_`` prefix so the
    thermoregulation consumer's behaviour is unchanged.

    Returns ``(payload | None, Verdict)``; ``payload`` is the raw result dict the consumer then
    uses (or ``None`` on REFUSE).

    ``provenance_detail`` / ``require_provenance_detail`` thread the structured provenance block
    into the SEALED validation: a caller may annotate an envelope for a
    validation exercise and/or require that a block is present. The production consumers pass
    neither, so the sealed trust boundary is unchanged. The RAW_SYNTHETIC stand-in envelope is
    never given a detail block (a demo must not fabricate provenance); if a caller nevertheless
    requires one on a demo path, that refusal is their explicit choice.

    ``require_uncertainty`` threads the strict uncertainty requirement into both the
    SEALED and the RAW_SYNTHETIC validation. It is off by default because the healthy
    metabolic-cost control carries no per-quantity uncertainty; a caller treating the result as a
    *scientific claim* can require it (or spread ``result_envelope_v1.SCIENTIFIC_REQUIREMENTS``).

    Scientific acceptance (the strict path, which ``auto`` resolves to for a sealed envelope)
    refuses a producer whose own gates failed (``producer gates failed``; the adapted public
    control with ``speed_sanity_ok`` / ``expected_direction_pass`` false is refused) and requires
    the sealed payload to equal the registry's ``payload_sha256`` pin (``no trusted artifact
    pin`` / ``artifact not approved``). ``require_producer_gates_pass`` only matters on the
    explicit synthetic demo path. The RAW_SYNTHETIC stand-in builds its envelope with no producer
    gates (``{}`` / ``True``) and can only be accepted as ``ACCEPT_SYNTHETIC_DEMO``.
    """
    kind = classify_input(raw_obj)
    resolved = resolve_mode(mode, kind)

    entry = registry_entry if isinstance(registry_entry, dict) else {}
    expected_producer_sha256 = entry.get("producer_sha256")
    expected_producer_name = entry.get("producer_name")
    registry_input = entry.get("input_sha256")
    registry_payload_pin = entry.get("payload_sha256")
    effective_input = registry_input if registry_input is not None else current_input_sha256

    allowed = allowed_provenance_for(
        opt_in=("synthetic", "unknown") if resolved == "synthetic_demo" else ()
    )

    # The explicit demo path may consume a ``synthetic_demo`` envelope; the scientific strict
    # path accepts only steady-state walking. Both paths validate the regime.
    accepted_regimes = THERMO_ACCEPTED_REGIMES | (
        frozenset({"synthetic_demo"}) if resolved == "synthetic_demo" else frozenset()
    )

    if kind == SEALED:
        try:
            # Bind the receipt's producer hash to the registry pin as well as the payload hash.
            payload = verify_receipt(
                raw_obj,
                producer_path=producer_path,
                expected_producer_sha256=expected_producer_sha256,
            )
        except ValueError as exc:
            return None, Verdict("REFUSE", (str(exc),), (), "refused")
        envelope, result = _envelope_and_result(payload)
        if envelope is None:
            return None, Verdict("REFUSE", (_SEALED_NOT_ENVELOPE,), (), "refused")
        if provenance_detail is not None:
            # Caller-supplied annotation (validation exercise): shallow-copy the producer so the
            # sealed payload object itself is not mutated.
            envelope = dict(envelope)
            envelope["producer"] = dict(envelope.get("producer") or {})
            envelope["producer"]["provenance_detail"] = provenance_detail
        verdict = validate_for_consumer(
            envelope,
            expected_quantities=THERMO_QUANTITY_SPECS,
            allowed_provenance=allowed,
            mode=resolved,
            current_input_sha256=effective_input,
            expected_producer_sha256=expected_producer_sha256,
            expected_producer_name=expected_producer_name,
            expected_regions=THERMO_REGION_VOCABULARY,
            expected_regimes=accepted_regimes,
            plausibility=thermo_plausibility(result),
            require_provenance_detail=require_provenance_detail,
            require_uncertainty=require_uncertainty,
            require_producer_gates_pass=require_producer_gates_pass,
            expected_payload_sha256=registry_payload_pin,
        )
        if verdict.state == "REFUSE":
            return None, verdict
        return result, verdict

    if kind in (ENVELOPE_UNSEALED, RAW_UNKNOWN):
        failure = _UNSEALED_FAILURE if kind == ENVELOPE_UNSEALED else _RAW_UNKNOWN_FAILURE
        return None, Verdict("REFUSE", (failure,), (), "refused")

    # RAW_SYNTHETIC: build the declared envelope (values from the raw object) and classify.
    envelope = make_envelope(
        raw_obj,
        producer_name="metabolic_cost_synthetic_standin",
        producer_sha256=None,
        provenance="synthetic",
        input_sha256=None,
        quantity_specs=THERMO_QUANTITY_SPECS,
        regime="synthetic_demo",
    )
    verdict = validate_for_consumer(
        envelope,
        expected_quantities=THERMO_QUANTITY_SPECS,
        allowed_provenance=allowed,
        mode=resolved,
        current_input_sha256=None,
        expected_regions=THERMO_REGION_VOCABULARY,
        expected_regimes=accepted_regimes,
        plausibility=thermo_plausibility(raw_obj),
        require_provenance_detail=require_provenance_detail,
        require_uncertainty=require_uncertainty,
        require_producer_gates_pass=require_producer_gates_pass,
    )
    if verdict.state == "REFUSE":
        return None, verdict
    return raw_obj, verdict


def preflight_thermoregulation(
    raw_obj,
    *,
    mode="auto",
    registry_entry=None,
    current_input_sha256=None,
    producer_path=None,
    provenance_detail=None,
    require_provenance_detail=False,
    require_uncertainty=False,
    require_producer_gates_pass=None,
):
    """Thin wrapper over :func:`preflight_metabolic_cost` for the thermoregulation consumer.

    Kept as the stable, named entry point used by ``thermoregulation.py`` and its tests. It
    forwards every argument verbatim, so the thermoregulation behaviour is byte-identical to the
    historical implementation (only where the shared body lives changed, not what it does).
    """
    return preflight_metabolic_cost(
        raw_obj,
        mode=mode,
        registry_entry=registry_entry,
        current_input_sha256=current_input_sha256,
        producer_path=producer_path,
        provenance_detail=provenance_detail,
        require_provenance_detail=require_provenance_detail,
        require_uncertainty=require_uncertainty,
        require_producer_gates_pass=require_producer_gates_pass,
    )
