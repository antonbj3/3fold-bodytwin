"""Adapter: wrap a region payload (``bodytwin.geometry``) into a :class:`ResultEnvelope`.

Division of ownership: ``bodytwin.framework`` owns the result envelope;
``bodytwin.geometry`` owns the region payload schema (``i2_region_example_v1``). This module is the seam between them. It is stdlib +
:mod:`bodytwin.framework.result_envelope_v1` only, performs no I/O, and never modifies the geometry package.

Region-payload sources in this repository (read-only; hashes recorded when the seam was
built)::

    examples/geometry/region_example.json          sha256 9bdc4be96021fe8c7776f12b8d24a4f23d6263d6f5a131e61ca3a93b6ed05e9c
    src/bodytwin/geometry/i2_payload_validate_v1.py sha256 268ca3854ac28671fd7cacd1d82a04f1d08f49af427ce2d9d6df8a330ec87042
    src/bodytwin/geometry/material_roster_v1.py     sha256 1998eec14cd3195f35705fb3d08324a765b02e23e9fde31a0fabbfb69e08ae2a

``examples/geometry/i1_region_example.py`` is the builder of the example payload; it is recorded
for provenance only and is never imported by this module.

I2's structural validator (``validate_payload(payload) -> list[str]``) is *injected* by the
caller (``bodytwin.geometry.i2_payload_validate_v1.validate_payload``); it is pure and
imports only ``math``. Passing it here refuses a malformed I2 payload *before* it is sealed into an envelope.

Mapping (I2 payload -> I1 quantities)::

    regions[<rid>].mass           -> regions.<rid>.mass          (kg,   region_id <rid>)
    regions[<rid>].heat_capacity  -> regions.<rid>.heat_capacity (J/K,  region_id <rid>)
    totals.mass                   -> totals.mass                 (kg,   region_id None)
    totals.heat_capacity          -> totals.heat_capacity        (J/K,  region_id None)
    area.<kind>                   -> area.<kind>                 (mm2,  region_id None)

``regions[<rid>].frame``, when present and not ``None``, is copied verbatim onto BOTH
of that region's quantities; ``area.frame`` likewise onto ``area.<kind>``. ``None``/absent keeps
the scalar convention (:data:`...result_envelope_v1.SCALAR_FRAME`) -- the adapter never invents a
frame. A declared frame may be a descriptor dict (``{origin_mm, axes, pitch_mm, units, declared}``
or a string; the shared validator compares frames canonically (a dict via sorted
JSON), so a consumer that requires a different frame is refused ``frame mismatch``.

``area.<kind>`` encodes I2's ``area.kind`` in the dotted quantity key itself, one of
``exterior_area`` / ``material_interface_area`` / ``functional_exchange_area`` (the identical names
used by the geometry contract). A consumer must require the exact ``area.<kind>`` key
(:func:`area_expected_quantities`), so an ``exterior_area`` (a wall) or a
``material_interface_area`` (a shared face) can never silently stand in for a
``functional_exchange_area``. The geometry package itself refuses to return a functional exchange area without a
caller declaration; this adapter mirrors that refusal at the envelope seam.

Quantities default to scalar static geometry (``frame='scalar'``, ``time=None``) unless the payload declares
a frame; the envelope declares ``validity.regime = None`` because no exercise regime is
claimed for a geometry record.

Each region may declare an explicit ``regions[<rid>].material_region_id`` (a
region -> material binding, e.g. ``SYNTH-BOX-1 -> MSK-MUSCLE``). The adapter copies it verbatim
into BOTH region quantity specs (``regions.<rid>.mass`` and ``regions.<rid>.heat_capacity``) so a
region/material swap is detectable, and ALWAYS refuses a present id that is not a member of
``I2_MATERIAL_REGION_IDS`` (the frozen material roster) with a ``"material region ..."`` error --
independent of any consumer. A consumer enforces the exact binding via
``preflight_region_payload(..., material_region_map={local_rid: material_id})`` (or, at the
envelope level, ``validate_for_consumer(..., expected_material_regions=...)``). When no map is
supplied the binding is not enforced, but the roster vocabulary still applies.

Interop gaps declared here (not hidden):
  * I2's ``provenance`` vocabulary (``measured``/``calibrated``/``literal_cited``/``assumption``/
    ``computed``/``UNKNOWN``) is normalised onto I1's canonical vocabulary;
    ``computed`` -> ``unknown`` and ``assumption`` -> ``unknown`` (never silently promoted).
  * ``totals.*`` are aggregate quantities with no single region; their ``region_id`` is ``None``
    by contract. The shared validator's ``expected_regions`` vocabulary therefore only applies to
    quantities that actually declare a region.
  * The payload may carry an optional per-region / per-area ``frame`` (a descriptor dict or a string); it is
    now consumed verbatim and an absent/``None`` frame is recorded as ``scalar`` rather
    than invented. The payload may
    carry an optional per-quantity ``uncertainty`` (the payload schema allows it); when present it is
    copied verbatim into the envelope quantity spec, and when absent it stays ``None`` -- the
    adapter never invents an uncertainty value.
  * The payload carries no ``area.provenance`` and its structural validator does not require one. The adapter
    therefore stores the canonical ``"unknown"`` for an area whose block declares no provenance
    (rather than leaving it ``None``, which the shared validator treats as "not checked"), so a
    consumer that *explicitly expects* the area quantity refuses an unprovenanced area as
    ``disallowed provenance``. An undeclared area provenance is deliberately NOT fed into the
    envelope's top-level :func:`derive_provenance` input: a structural area must not drag an
    otherwise calibrated payload down to ``unknown`` when no one asked for the area.
"""
from __future__ import annotations

from bodytwin.framework.region_vocabulary_v1 import I2_MATERIAL_REGION_IDS
from bodytwin.framework.result_envelope_v1 import (
    SCALAR_FRAME,
    Verdict,
    allowed_provenance_for,
    make_envelope,
    normalize_provenance,
    payload_sha256,
    validate_for_consumer,
)

# The frozen I2 payload schema id this adapter accepts.
I2_SCHEMA_ID = "i2_region_example_v1"

# Least-trustworthy -> most-trustworthy. ``derive_provenance`` returns the worst of its inputs so
# a derived value can never inherit more trust than its weakest input.
_PROVENANCE_TRUST_ORDER = (
    "unknown",
    "tuned",
    "synthetic",
    "literal_cited",
    "calibrated_public_model",
    "calibrated",
    "measured",
)

# Canonical unit expected per mapped quantity (also used by ``region_expected_quantities``).
_REGION_QUANTITY_UNITS = (("mass", "kg"), ("heat_capacity", "J/K"))

# Canonical unit expected for an area-kind-tagged quantity (``area.<kind>``). The payload uses mm2; the
# kind (not the unit) is what distinguishes a wall from a material interface from an exchange area,
# so the unit is still pinned to make a length/area (mm vs mm2) swap a hard ``unit mismatch``.
_AREA_QUANTITY_UNITS = "mm2"


def derive_provenance(provenances) -> str:
    """Return the least-trustworthy canonical provenance over ``provenances``.

    Each input is mapped through :func:`normalize_provenance`; the result is the lowest term in
    the trust order ``unknown < tuned < synthetic < literal_cited < calibrated_public_model <
    calibrated < measured`` (assume the worst). An empty input yields ``"unknown"``.
    """
    canonical = [normalize_provenance(p) for p in (provenances or ())]
    if not canonical:
        return "unknown"
    return min(canonical, key=lambda term: _PROVENANCE_TRUST_ORDER.index(term))


def region_expected_quantities(
    region_ids,
    *,
    require_totals=True,
    material_region_map=None,
    frames=None,
    area_kinds=None,
) -> dict:
    """Build the expected-quantity contract for a region payload.

    For every region id: ``regions.<rid>.mass`` (kg) and ``regions.<rid>.heat_capacity`` (J/K),
    both ``region_id=<rid>``/``time=None``/required. When ``require_totals`` (default):
    ``totals.mass`` (kg) and ``totals.heat_capacity`` (J/K) with ``region_id=None``. When
    ``area_kinds`` (a convenience) is given, the matching
    :func:`area_expected_quantities` contract is merged in as well, so one call can build the whole
    region + area contract; ``preflight_region_payload`` uses exactly that.

    The expected ``frame`` is :data:`...result_envelope_v1.SCALAR_FRAME` by default. ``frames``
    is an optional ``{local_id: expected_frame}`` mapping (the same shape threaded
    through :func:`preflight_region_payload`); a requested region present in it gets that frame on
    BOTH of its expected quantity specs, so a consumer can require a declared non-scalar region
    frame. A frame value is used verbatim (a dict is compared canonically by the shared validator's
    sorted-JSON rule); a region absent from the mapping keeps ``SCALAR_FRAME``, and an explicit
    ``{rid: None}`` entry is honoured as "expect null". Unknown keys in ``frames`` are ignored (the
    same mapping is also passed to :func:`area_expected_quantities`, whose keys are area kinds), so
    there is deliberately no "unknown frame key" ``ValueError``.

    ``material_region_map`` is an optional ``{local_region_id: material_id}`` mapping:
    each requested region that has an entry gets ``"material_region_id": <material_id>`` added to
    BOTH of its expected quantity specs, so the contract itself records the expected region ->
    material binding. A map entry whose local id is NOT one of the requested ``region_ids`` is a
    caller error (a silent typo would silently drop the binding), so it raises a clear
    :class:`ValueError`; a requested region absent from the map simply gets no binding. The key is
    omitted when there is no binding, keeping the historical expected-quantity shape.
    """
    unknown = sorted(set(material_region_map or ()) - set(region_ids))
    if unknown:
        raise ValueError(
            f"material region map: unknown local region id(s) {unknown} not in the requested "
            f"region_ids {list(region_ids)}"
        )
    frames = frames or {}
    expected = {}
    for region_id in region_ids:
        material_region_id = (material_region_map or {}).get(region_id)
        frame = frames[region_id] if region_id in frames else SCALAR_FRAME
        for name, units in _REGION_QUANTITY_UNITS:
            spec = {
                "units": units,
                "region_id": region_id,
                "frame": frame,
                "time": None,
                "required": True,
            }
            if material_region_id is not None:
                spec["material_region_id"] = material_region_id
            expected[f"regions.{region_id}.{name}"] = spec
    if require_totals:
        for name, units in _REGION_QUANTITY_UNITS:
            expected[f"totals.{name}"] = {
                "units": units,
                "region_id": None,
                "frame": SCALAR_FRAME,
                "time": None,
                "required": True,
            }
    if area_kinds:
        # Convenience: one call can build the region+area contract. The same ``frames``
        # mapping is forwarded so an area kind present in it gets its expected frame.
        expected.update(area_expected_quantities(area_kinds, frames=frames))
    return expected


def area_expected_quantities(area_kinds, *, required=True, frames=None) -> dict:
    """Build the expected-quantity contract for I2's area block.

    For every kind in ``area_kinds`` emit ``area.<kind>`` with units ``mm2``, ``region_id=None``,
    ``time=None`` and ``required``. The kind lives in the key, so a consumer's contract names the
    exact area semantics it wants: requiring ``area.functional_exchange_area`` is NOT satisfied by
    an ``area.material_interface_area`` (the shared validator then reports the stable ``missing
    quantity: area.functional_exchange_area`` for the wall-area envelope). ``required=False`` makes
    the area contract optional (an abstention is then flagged rather than failed).

    ``frames`` is the same optional ``{local_id: expected_frame}`` mapping as
    :func:`region_expected_quantities`; here a kind present in it gets that frame on the
    ``area.<kind>`` spec, so a consumer can require a declared non-scalar area frame. A kind absent
    from the mapping keeps ``SCALAR_FRAME``; frames keyed by a region id are simply ignored.
    """
    frames = frames or {}
    return {
        f"area.{kind}": {
            "units": _AREA_QUANTITY_UNITS,
            "region_id": None,
            "frame": frames[kind] if kind in frames else SCALAR_FRAME,
            "time": None,
            "required": required,
        }
        for kind in area_kinds
    }


def _quantity_units_and_provenance(quantity):
    """``(units, provenance)`` of an I2 quantity object (or ``(None, None)`` for a non-object)."""
    if isinstance(quantity, dict):
        return quantity.get("units"), quantity.get("provenance")
    return None, None


def _quantity_uncertainty(quantity):
    """Optional ``uncertainty`` of an I2 quantity object (or ``None`` for a non-object/absent).

    The payload schema may add an ``uncertainty`` field to a quantity object; when it is present it is
    copied verbatim (including ``0.0``), and when it is absent it stays ``None``. The adapter
    never invents an uncertainty value.
    """
    if isinstance(quantity, dict):
        return quantity.get("uncertainty")
    return None


def make_region_envelope(
    i2_payload,
    *,
    producer_name,
    producer_sha256,
    input_sha256,
    structural_validator=None,
    command=None,
    seed=None,
    provenance_detail=None,
    gates=None,
    overall_pass=None,
):
    """Wrap an I2 region payload into a ResultEnvelope.

    Returns ``(envelope, [])`` on success or ``(None, errors)`` where ``errors`` is a non-empty
    list of stable messages. A non-object / wrong-``schema`` payload is refused with a message
    containing ``"schema"``; when ``structural_validator`` is supplied (I2's real
    ``validate_payload``) its errors are returned verbatim.

    The payload's ``area`` block (when it carries a nonempty ``kind`` string) is carried as ONE
    quantity under the dotted key ``area.<kind>`` (units ``mm2``, ``region_id=None``), so a
    consumer can require the exact area semantics it wants. An undeclared area
    provenance becomes the per-quantity ``"unknown"`` but is not folded into the top-level
    provenance; see the module docstring.

    A region's optional ``regions[<rid>].frame`` (a descriptor dict or string) is copied
    verbatim onto BOTH of the region's quantities, and ``area.frame`` onto ``area.<kind>``; a
    ``None``/absent frame keeps :data:`...result_envelope_v1.SCALAR_FRAME`, so a frame is never
    invented and the frozen I2 example (no frame) is byte-identical to before.

    ``provenance_detail`` is passed through to :func:`make_envelope` (optional structured
    ``data`` / ``calibration`` / ``model`` / ``execution`` block); ``None`` omits the key.

    ``gates`` / ``overall_pass`` carry an upstream producer's own gate outcome into the wrapped
    envelope. I2 payloads carry no such block by default, so both are ``None`` and
    :func:`make_envelope` stores the historical ``{}`` / ``True``.

    A region's optional ``material_region_id`` is copied verbatim into BOTH of its
    quantity specs (``regions.<rid>.mass`` and ``regions.<rid>.heat_capacity``); a present id that
    is not a member of :data:`bodytwin.framework.region_vocabulary_v1.I2_MATERIAL_REGION_IDS`
    is refused with ``(None, ["material region ... not registered ..."])`` (stable substring
    ``"material region"``). An absent id is never invented.
    """
    if not isinstance(i2_payload, dict):
        return None, [
            f"region payload: schema must be {I2_SCHEMA_ID!r} (got {type(i2_payload).__name__})"
        ]
    if i2_payload.get("schema") != I2_SCHEMA_ID:
        return None, [f"region payload: schema must be {I2_SCHEMA_ID!r}"]

    if structural_validator is not None:
        errors = list(structural_validator(i2_payload))
        if errors:
            return None, errors

    regions = i2_payload.get("regions")
    totals = i2_payload.get("totals")
    if not isinstance(totals, dict):
        totals = {}

    region_map = {}
    region_material_ids = {}
    region_frames = {}
    provenance_inputs = []
    if isinstance(regions, list):
        for region in regions:
            if not isinstance(region, dict):
                continue
            region_id = region.get("region_id")
            if not isinstance(region_id, str) or not region_id:
                continue
            # The region -> material binding is enforced HERE, independent of any
            # consumer. A present ``material_region_id`` must be a member of the frozen material
            # roster; an unregistered (or non-string) id is refused with a stable "material region"
            # substring rather than being copied into the envelope. An absent id stays absent (the
            # adapter never invents a binding).
            material_region_id = region.get("material_region_id")
            if material_region_id is not None:
                if (
                    not isinstance(material_region_id, str)
                    or material_region_id not in I2_MATERIAL_REGION_IDS
                ):
                    return None, [
                        f"material region {material_region_id!r} is not registered in the I2 "
                        f"material roster (registered material regions: "
                        f"{list(I2_MATERIAL_REGION_IDS)})"
                    ]
                region_material_ids[region_id] = material_region_id
            region_map[region_id] = {
                "mass": region.get("mass"),
                "heat_capacity": region.get("heat_capacity"),
            }
            # The payload may declare an optional per-region ``frame`` (a frame descriptor dict
            # or a string). It is captured verbatim; ``None``/absent means the scalar convention,
            # which is what the envelope records -- a frame is never invented.
            region_frames[region_id] = region.get("frame")
            for name, _ in _REGION_QUANTITY_UNITS:
                provenance_inputs.append(region.get(name, {}).get("provenance")
                                         if isinstance(region.get(name), dict) else None)
    for name, _ in _REGION_QUANTITY_UNITS:
        provenance_inputs.append(totals.get(name, {}).get("provenance")
                                 if isinstance(totals.get(name), dict) else None)

    # Carry the payload's single area block as ONE area-kind-tagged quantity. The kind is encoded
    # in the dotted key (``area.<kind>``) so a consumer must require the exact semantics it wants;
    # an ``exterior_area``/``material_interface_area`` can never satisfy an
    # ``area.functional_exchange_area`` contract. The area is additive: when no consumer supplies an
    # area contract the quantity is ignored (it declares ``region_id=None``, so the shared
    # ``expected_regions`` check exempts it).
    area = i2_payload.get("area")
    area_kind = None
    if isinstance(area, dict):
        kind = area.get("kind")
        if isinstance(kind, str) and kind:
            area_kind = kind
            # The provenance is fed to ``derive_provenance`` ONLY when the block actually declares
            # one. The region example declares no ``area.provenance``; including an absent/None one would
            # normalise to ``unknown`` and silently drag an otherwise calibrated envelope down.
            if area.get("provenance") is not None:
                provenance_inputs.append(area.get("provenance"))

    # ``result_obj`` is keyed by region id so the dotted quantity keys (``regions.<rid>.mass``)
    # resolve through the shared ``get_dotted`` traversal (the payload stores ``regions`` as a list). The
    # area block is nested under its kind so ``area.<kind>`` resolves to the block (which carries
    # ``value``/``units``/...).
    result_obj = {"regions": region_map, "totals": totals}
    if area_kind is not None:
        result_obj["area"] = {area_kind: area}

    quantity_specs = {}
    for region_id, quantities in region_map.items():
        material_region_id = region_material_ids.get(region_id)
        # Consume the region's declared frame; an absent/None frame keeps the scalar
        # convention. ``make_envelope`` stores the value verbatim (a dict is compared canonically
        # by the shared validator's sorted-JSON rule).
        region_frame = region_frames.get(region_id)
        if region_frame is None:
            region_frame = SCALAR_FRAME
        for name, _ in _REGION_QUANTITY_UNITS:
            units, provenance = _quantity_units_and_provenance(quantities.get(name))
            spec = {
                "units": units,
                "region_id": region_id,
                "frame": region_frame,
                "time": None,
                # Copied verbatim from the I2 quantity object; ``None`` when absent (never invented).
                "uncertainty": _quantity_uncertainty(quantities.get(name)),
                "provenance": provenance,
            }
            # The SAME region -> material binding is copied into BOTH region quantities
            # (mass and heat_capacity), so a consumer can require either one and a swapped
            # region/material binding is detectable. Omitted when the region declares none.
            if material_region_id is not None:
                spec["material_region_id"] = material_region_id
            quantity_specs[f"regions.{region_id}.{name}"] = spec
    for name, _ in _REGION_QUANTITY_UNITS:
        units, provenance = _quantity_units_and_provenance(totals.get(name))
        quantity_specs[f"totals.{name}"] = {
            "units": units,
            "region_id": None,
            "frame": SCALAR_FRAME,
            "time": None,
            "uncertainty": _quantity_uncertainty(totals.get(name)),
            "provenance": provenance,
        }

    if area_kind is not None:
        # The payload's ``area`` block has no declared ``provenance``. An undeclared provenance is stored as
        # the canonical ``"unknown"`` (rather than left ``None``, which the shared validator treats
        # as "not present, nothing to check"), so a consumer that explicitly expects the area
        # quantity refuses an unprovenanced area as ``disallowed provenance`` -- an unprovenanced
        # area is not scientific evidence. A *declared* provenance is copied verbatim (e.g.
        # ``calibrated``) and normalised by the validator like every other mapping. This per-quantity
        # ``"unknown"`` does NOT affect ``provenance_inputs`` above (only a declared area provenance
        # is folded into the top-level derive), so it cannot drag the envelope's top-level
        # provenance.
        area_provenance = area.get("provenance")
        if area_provenance is None:
            area_provenance = normalize_provenance(None)
        # The payload's area block may also declare an optional frame (a descriptor dict or a
        # string); it is consumed verbatim, with absent/None keeping the scalar convention.
        area_frame = area.get("frame")
        if area_frame is None:
            area_frame = SCALAR_FRAME
        quantity_specs[f"area.{area_kind}"] = {
            "units": area.get("units"),
            "region_id": None,
            "frame": area_frame,
            "time": None,
            # The payload's area block may grow an optional ``uncertainty``; copied verbatim when present and
            # left ``None`` when absent, exactly like the region/totals quantity objects.
            "uncertainty": area.get("uncertainty"),
            "provenance": area_provenance,
        }

    top_provenance = (
        "synthetic"
        if i2_payload.get("synthetic") is True
        else derive_provenance(provenance_inputs)
    )

    envelope = make_envelope(
        result_obj,
        producer_name=producer_name,
        producer_sha256=producer_sha256,
        provenance=top_provenance,
        input_sha256=input_sha256,
        quantity_specs=quantity_specs,
        command=command,
        seed=seed,
        regime=None,
        provenance_detail=provenance_detail,
        gates=gates,
        overall_pass=overall_pass,
    )
    return envelope, []


def preflight_region_payload(
    i2_payload,
    *,
    region_ids,
    area_kinds=None,
    material_region_map=None,
    frames=None,
    mode="strict",
    structural_validator=None,
    allowed_provenance=None,
    current_input_sha256=None,
    expected_producer_sha256=None,
    expected_producer_name=None,
    producer_name="i2_region_example",
    producer_sha256=None,
    input_sha256=None,
    provenance_detail=None,
    require_provenance_detail=False,
    require_uncertainty=False,
    require_producer_gates_pass=None,
    gates=None,
    overall_pass=None,
    expected_payload_sha256=None,
):
    """Validate an I2 region payload for a caller that expects ``region_ids``.

    Scientific acceptance (``mode="strict"``) requires ``expected_payload_sha256``: the trusted
    pin of the region payload itself (:func:`...result_envelope_v1.payload_sha256` of
    ``i2_payload``). A missing pin is refused (``no trusted artifact pin``) and a different payload
    is refused (``artifact not approved``). Failed producer gates are always refused in strict
    mode. ``literal_cited`` values are not scientific provenance.

    Returns ``(envelope | None, Verdict)``. The envelope is ``None`` only when the payload could
    not be wrapped (bad schema or structural errors) — a REFUSE *after* wrapping still returns the
    envelope, so a caller can inspect exactly what was refused. ``region_ids`` is both the
    expected-quantity contract and the allowed region vocabulary.

    ``area_kinds``, when given, adds :func:`area_expected_quantities` to the contract:
    the caller must name the exact area semantics it wants (e.g. ``["functional_exchange_area"]``),
    so an ``area.exterior_area`` or ``area.material_interface_area`` envelope is refused with the
    stable ``missing quantity: area.<kind>``. ``None`` (default) expects no area and leaves the
    additive ``area.<kind>`` quantity ignored, exactly as before.

    ``frames`` is an optional ``{local_id: expected_frame}`` mapping threaded into
    both :func:`region_expected_quantities` and :func:`area_expected_quantities`: a region id (or
    area kind) present in it makes the consumer require that exact frame on the corresponding
    expected quantity, compared canonically by the shared validator (a dict is sorted-JSON
    compared). A local id absent from the mapping keeps the scalar convention, so the default
    contract is unchanged and an undeclared frame is never invented. Keys for entities the caller
    did not request are ignored (see :func:`region_expected_quantities`), so the single mapping can
    carry both region-id and area-kind entries.

    ``material_region_map`` is an optional ``{local_region_id: material_id}`` binding
    the caller requires. When given, each requested region's expected quantities carry the
    expected ``material_region_id`` (via :func:`region_expected_quantities`) and the exact binding
    is enforced by the shared validator (``expected_material_regions``), so a swapped region ->
    material binding is refused with ``material region mismatch`` and an absent one with
    ``missing material region``. The frozen material vocabulary
    (:data:`...region_vocabulary_v1.I2_MATERIAL_REGION_IDS`) is ALWAYS passed as well, so an
    unregistered ``material_region_id`` anywhere is refused with ``unknown material region`` (an
    unregistered id in the payload is already refused earlier by :func:`make_region_envelope`).
    When ``material_region_map`` is ``None`` (default) the binding is NOT enforced -- but the
    roster vocabulary still applies. A map entry for a local id outside ``region_ids`` is a caller
    error and is surfaced as a REFUSE (``material region map: unknown local region id(s) ...``).

    ``provenance_detail`` is written into the wrapped envelope's producer; a caller may also
    ``require_provenance_detail=True`` to refuse an envelope that carries no block. A caller may
    ``require_uncertainty=True`` to refuse any required quantity whose I2 quantity
    object carried no ``uncertainty`` (the adapter threads a present field through verbatim; an
    absent one stays ``None`` and is refused under this strict profile).

    ``gates`` / ``overall_pass`` are forwarded to :func:`make_region_envelope` so a
    caller can carry an upstream producer's own gate outcome into the region envelope; both
    default to ``None`` (the historical ``{}`` / ``True``). An ``overall_pass`` other than ``True``
    is refused in strict mode; in ``synthetic_demo`` mode it is a flag unless
    ``require_producer_gates_pass=True``.
    """
    envelope, errors = make_region_envelope(
        i2_payload,
        producer_name=producer_name,
        producer_sha256=producer_sha256,
        input_sha256=input_sha256,
        structural_validator=structural_validator,
        provenance_detail=provenance_detail,
        gates=gates,
        overall_pass=overall_pass,
    )
    if envelope is None:
        return None, Verdict("REFUSE", tuple(errors), (), "refused")

    allowed = allowed_provenance or allowed_provenance_for(
        opt_in=("synthetic", "unknown", "literal_cited") if mode == "synthetic_demo" else ()
    )
    # Trusted pin of the region payload itself. The envelope is derived deterministically from the
    # payload and the caller's arguments, so a matching payload pin is passed on as the pin of the
    # derived envelope; a missing or different pin leaves the envelope unpinned (refused).
    pin_failures = []
    envelope_pin = None
    if expected_payload_sha256 is not None:
        got_payload_sha = payload_sha256(i2_payload)
        if got_payload_sha == expected_payload_sha256:
            envelope_pin = payload_sha256(envelope)
        elif mode == "strict":
            pin_failures.append(
                f"artifact not approved: region payload sha256 {got_payload_sha!r} != trusted pin "
                f"{expected_payload_sha256!r}"
            )
    try:
        expected_quantities = region_expected_quantities(
            region_ids,
            material_region_map=material_region_map,
            frames=frames,
            # The area contract is additive and carries the same ``frames`` mapping.
            area_kinds=area_kinds,
        )
    except ValueError as exc:
        # A material_region_map entry keyed by a local id the caller did not request is a
        # contract-construction error; surface it as a REFUSE (keeping the wrapped envelope) rather
        # than letting a validation entry point raise.
        return envelope, Verdict("REFUSE", (str(exc),), (), "refused")
    # Derive the exact expected binding map from the annotated expected quantities (only
    # requested regions that had a map entry), and ALWAYS pass the frozen material vocabulary so an
    # unregistered id is refused even when no explicit binding map is supplied.
    expected_material_regions = (
        {
            key: spec["material_region_id"]
            for key, spec in expected_quantities.items()
            if "material_region_id" in spec
        }
        if material_region_map is not None
        else None
    )
    verdict = validate_for_consumer(
        envelope,
        expected_quantities=expected_quantities,
        allowed_provenance=allowed,
        mode=mode,
        current_input_sha256=current_input_sha256,
        expected_producer_sha256=expected_producer_sha256,
        expected_producer_name=expected_producer_name,
        expected_regions=frozenset(region_ids),
        expected_material_regions=expected_material_regions,
        material_region_vocabulary=frozenset(I2_MATERIAL_REGION_IDS),
        require_provenance_detail=require_provenance_detail,
        require_uncertainty=require_uncertainty,
        require_producer_gates_pass=require_producer_gates_pass,
        expected_payload_sha256=envelope_pin,
    )
    if pin_failures:
        verdict = Verdict("REFUSE", tuple(pin_failures) + tuple(verdict.failures), verdict.flags,
                          "refused")
    return envelope, verdict
