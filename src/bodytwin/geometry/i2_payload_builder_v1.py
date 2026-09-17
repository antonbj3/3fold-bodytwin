"""Builder for the ``i2_region_example_v1`` region payload.

Additive; changes no existing symbol; pure, no I/O (it does not read
examples/geometry/region_example.json).

A reusable, deterministic builder for the frozen ``i2_region_example_v1`` payload. The example
script (``examples/geometry/i1_region_example.py``) hand-rolls the payload; this module factors
that shape into one pure function so other fixtures can produce the same record without copying
the hand-rolled assembly.

It reuses, unmodified:
  * ``region_mass_v1``       -- Region, MaterialProperty, integrate_region, area_mm2,
                                interface_area_by_pair_mm2, region_boundary_volume_mm3,
                                mass_moments_surface, total_quantity, MASS_BASIS;
  * ``i2_payload_validate_v1`` -- validate_payload (called on the finished payload);
  * ``i2_region_registration_v1`` -- registration_issues (only from ``build_and_check``);
  * ``material_roster_v1``   -- lookup / is_known (resolve a ``*_key`` to a KNOWN value or
                                record the declared absence; never a default).

Design:
  * ``region_specs`` are matched to the interface region labels in sorted-label order. Every
    label selectable by the data is consumed once; a spec may pin ``region_label`` explicitly
    (for a subset), but then all specs must pin it. A count/label mismatch is refused.
  * A material input must be declared for every region: either an explicit ``MaterialProperty``
    (``density``/``specific_heat``) or a roster ``*_key``. A registered-but-unknown key (e.g.
    ``MSK-TENDON`` density) is an explicit declaration of absence and yields null plus a
    matching ``missing`` entry; a missing declaration is refused.
  * ``missing != default``: an unresolved value is null with a nonempty ``missing`` list.
  * The finished payload is passed through ``validate_payload``; a nonempty error list is
    raised as ``PayloadBuildError``. Duplicate region ids therefore surface through the
    validator. Registration is opt-in via ``build_and_check``.

No physiological claim is made: the built values are the caller's declared synthetic inputs.
"""
from __future__ import annotations

import math

import numpy as np

from . import (
    frame_v1,
    i2_payload_validate_v1,
    i2_region_registration_v1,
    material_roster_v1,
    region_mass_v1,
    thickness_v1,
    uncertainty_v1,
)

__all__ = ['PayloadBuildError', 'build_and_check', 'build_payload']

SCHEMA_ID = 'i2_region_example_v1'
AREA_UNITS = 'mm2'
DEFAULT_BASIS = 'region_volume'
DEFAULT_VOLUME_METHOD = 'outward boundary integral of a synthetic closed surface'

_AREA_METHODS = {
    'material_interface_area': 'shared-face count once, back->front',
    'exterior_area': 'exterior label-0 face area, summed once',
}


class PayloadBuildError(ValueError):
    """The frozen ``i2_region_example_v1`` payload could not be built or did not validate."""


def _quantity(q) -> dict:
    return {'name': q.name, 'value': q.value, 'units': q.units,
            'provenance': q.provenance, 'missing': list(q.missing)}


def _material(m) -> dict | None:
    """Material object shaped exactly like the frozen example (or an explicit null)."""
    if m is None:
        return None
    return {'name': m.name, 'value': m.value, 'units': m.units, 'provenance': m.provenance,
            'source': m.source, 'source_sha256': m.source_sha256}


def _frame_block(raw, region_id) -> dict:
    """Normalise a ``frame_v1.Frame`` or a plain frame dict into the payload frame block.

    Accepts a ``frame_v1.Frame`` (origin/axes/pitch/declared are read from it) or a dict
    already in the emitted shape ``{origin_mm, axes, pitch_mm, units, declared}``. Everything
    else, and any frame that violates the mm/orthonormal right-handed contract, raises
    ``PayloadBuildError``.
    """
    if isinstance(raw, frame_v1.Frame):
        frame = raw
    elif isinstance(raw, dict):
        for key in ('origin_mm', 'axes', 'pitch_mm'):
            if key not in raw:
                raise PayloadBuildError(
                    f'region {region_id!r}: frame dict requires {key!r}')
        units = raw.get('units', 'mm')
        if units != 'mm':
            raise PayloadBuildError(
                f"region {region_id!r}: frame units must be 'mm', got {units!r}")
        declared = raw.get('declared', False)
        if not isinstance(declared, bool):
            raise PayloadBuildError(
                f'region {region_id!r}: frame declared must be bool, got '
                f'{type(declared).__name__}')
        try:
            frame = frame_v1.Frame(origin_mm=raw['origin_mm'], axes=raw['axes'],
                                   pitch_mm=raw['pitch_mm'], pitch_declared=declared,
                                   units='mm')
        except frame_v1.FrameError as exc:
            raise PayloadBuildError(f'region {region_id!r}: frame: {exc}') from exc
    else:
        raise PayloadBuildError(
            f'region {region_id!r}: frame must be a frame_v1.Frame or dict, got '
            f'{type(raw).__name__}')
    axes = np.asarray(frame.axes, dtype=np.float64)
    return {
        'origin_mm': [float(value) for value in np.asarray(frame.origin_mm, dtype=np.float64)],
        'axes': [[float(value) for value in row] for row in axes],
        'pitch_mm': float(frame.pitch_mm),
        'units': 'mm',
        'declared': bool(frame.pitch_declared),
    }


def _coerce_uncertainty(name, raw, region_id) -> dict:
    """Validate one uncertainty entry and return its payload object."""
    if not isinstance(name, str) or not name:
        raise PayloadBuildError(
            f'region {region_id!r}: uncertainty name must be a nonempty string, got {name!r}')
    if isinstance(raw, uncertainty_v1.Uncertainty):
        entry = raw
    elif isinstance(raw, dict):
        fields = {key: value for key, value in raw.items()
                  if key not in ('name', 'uncertainty')}
        try:
            entry = uncertainty_v1.Uncertainty(**fields)
        except TypeError as exc:
            raise PayloadBuildError(
                f'region {region_id!r}: uncertainty {name!r} malformed: {exc}') from exc
        except uncertainty_v1.UncertaintyError as exc:
            raise PayloadBuildError(
                f'region {region_id!r}: uncertainty {name!r}: {exc}') from exc
    else:
        raise PayloadBuildError(
            f'region {region_id!r}: uncertainty {name!r} must be an '
            f'uncertainty_v1.Uncertainty or dict, got {type(raw).__name__}')
    return {'name': name, 'value': float(entry.value), 'units': entry.units,
            'kind': entry.kind, 'provenance': entry.provenance}


def _uncertainty_block(raw, region_id) -> list:
    """Normalise a spec ``uncertainty`` into a list of payload uncertainty objects.

    Accepts a ``{name: entry}`` mapping, a single uncertainty-shaped dict (with ``name``), or
    a list of entries where an entry is an uncertainty dict (with ``name``), a
    ``{name, uncertainty}`` wrapper, a ``(name, entry)`` pair, or an ``Uncertainty``. A bare
    ``Uncertainty`` has no name of its own, so it must arrive through one of the named forms.
    """
    entries = []
    if isinstance(raw, dict):
        if 'value' in raw or 'kind' in raw:
            entries.append((raw.get('name'), raw))
        else:
            for name in sorted(raw):
                entries.append((name, raw[name]))
    elif isinstance(raw, (list, tuple)):
        for index, item in enumerate(raw):
            if isinstance(item, (list, tuple)) and len(item) == 2 and isinstance(item[0], str):
                entries.append((item[0], item[1]))
            elif isinstance(item, dict) and 'name' in item:
                if 'uncertainty' in item:
                    entries.append((item['name'], item['uncertainty']))
                else:
                    entries.append((item['name'], item))
            elif isinstance(item, uncertainty_v1.Uncertainty):
                name = getattr(item, 'name', None)
                if not isinstance(name, str) or not name:
                    raise PayloadBuildError(
                        f'region {region_id!r}: uncertainty[{index}] is an Uncertainty with no '
                        'name; pass a {name: Uncertainty} mapping or a (name, Uncertainty) pair')
                entries.append((name, item))
            else:
                raise PayloadBuildError(
                    f'region {region_id!r}: uncertainty[{index}] must be an Uncertainty, a named '
                    f'dict or a (name, entry) pair, got {type(item).__name__}')
    else:
        raise PayloadBuildError(
            f'region {region_id!r}: uncertainty must be a dict or list, got '
            f'{type(raw).__name__}')
    return [_coerce_uncertainty(name, value, region_id) for name, value in entries]


def _thickness_block(raw, region_id) -> dict:
    """Normalise a ``thickness_v1.ThicknessDistribution`` or dict into the payload block."""
    if isinstance(raw, thickness_v1.ThicknessDistribution):
        dist = raw
    elif isinstance(raw, dict):
        units = raw.get('units', thickness_v1.UNITS)
        if units != thickness_v1.UNITS:
            raise PayloadBuildError(
                f'region {region_id!r}: thickness units must be {thickness_v1.UNITS!r}, '
                f'got {units!r}')
        fields = {key: value for key, value in raw.items() if key != 'units'}
        try:
            dist = thickness_v1.ThicknessDistribution(**fields)
        except TypeError as exc:
            raise PayloadBuildError(
                f'region {region_id!r}: thickness malformed: {exc}') from exc
        except thickness_v1.ThicknessError as exc:
            raise PayloadBuildError(f'region {region_id!r}: thickness: {exc}') from exc
    else:
        raise PayloadBuildError(
            f'region {region_id!r}: thickness must be a thickness_v1.ThicknessDistribution or '
            f'dict, got {type(raw).__name__}')
    return {
        'samples_mm': [float(sample) for sample in dist.samples_mm],
        'units': dist.units,
        'method': dist.method,
        'provenance': dist.provenance,
        'source': dist.source,
    }


def _require_interfaces(interfaces) -> None:
    if not isinstance(interfaces, dict):
        raise PayloadBuildError(f'interfaces: dict required, got {type(interfaces).__name__}')
    for key in ('vertices', 'faces', 'front', 'back'):
        if key not in interfaces:
            raise PayloadBuildError(f'interfaces: missing required key {key!r}')


def _interface_labels(interfaces) -> list:
    back = np.asarray(interfaces['back'])
    front = np.asarray(interfaces['front'])
    labels = {int(value) for value in back.tolist()} | {int(value) for value in front.tolist()}
    labels.discard(0)
    return sorted(labels)


def _outward_faces(interfaces, label) -> np.ndarray:
    faces = np.asarray(interfaces['faces'])
    back = np.asarray(interfaces['back'])
    front = np.asarray(interfaces['front'])
    return np.concatenate([faces[back == label], faces[front == label][:, ::-1]])


def _resolve_labels(interfaces, specs) -> list:
    labels = _interface_labels(interfaces)
    if not labels:
        raise PayloadBuildError('interfaces: no nonzero region labels found')
    explicit = [spec.get('region_label') for spec in specs]
    has_explicit = [value is not None for value in explicit]
    if any(has_explicit) and not all(has_explicit):
        raise PayloadBuildError('region_label: give it for every region spec or for none')
    if all(has_explicit):
        try:
            chosen = [int(value) for value in explicit]
        except (TypeError, ValueError) as exc:
            raise PayloadBuildError('region_label: integer labels required') from exc
        if len(set(chosen)) != len(chosen):
            raise PayloadBuildError(f'region_label: duplicate labels {chosen}')
        unknown = sorted(set(chosen) - set(labels))
        if unknown:
            raise PayloadBuildError(
                f'region_label: {unknown} not present in interface labels {labels}')
        return chosen
    if len(labels) != len(specs):
        raise PayloadBuildError(
            f'{len(specs)} region specs but the interfaces carry {len(labels)} region labels '
            f'{labels}; pin region_label on every spec to select a subset')
    return labels


def _resolve_material(spec, property_name):
    """Return a MaterialProperty, or None for an explicitly declared roster absence.

    ``spec[property_name]`` (a MaterialProperty) and ``spec[property_name + '_key']`` (a roster
    key) are mutually exclusive. An unregistered key is refused; a registered key with no audited
    value keeps the declared absence (null + matching ``missing``), never a generic default.
    """
    region_id = spec.get('region_id')
    explicit = spec.get(property_name)
    key = spec.get(f'{property_name}_key')
    if explicit is not None and key is not None:
        raise PayloadBuildError(
            f'region {region_id!r}: give {property_name!r} or {property_name}_key, not both')
    if explicit is not None:
        if not isinstance(explicit, region_mass_v1.MaterialProperty):
            raise PayloadBuildError(
                f'region {region_id!r}: {property_name!r} must be a MaterialProperty, '
                f'got {type(explicit).__name__}')
        return explicit
    if key is None:
        raise PayloadBuildError(
            f'region {region_id!r}: missing {property_name!r} or {property_name}_key; a '
            f'material input must be declared (missing != default)')
    try:
        row = material_roster_v1.lookup(key, property_name)
    except KeyError as exc:
        raise PayloadBuildError(str(exc)) from exc
    return row if material_roster_v1.is_known(row) else None


def build_payload(interfaces, region_specs, *, area_kind='material_interface_area',
                  synthetic=True, source_sha256=None, area_to_volume=False) -> dict:
    """Build the frozen ``i2_region_example_v1`` payload; raise ``PayloadBuildError`` otherwise.

    ``interfaces`` is a label/tetra interface dict (vertices/faces/front/back). ``region_specs``
    is an iterable of dicts with ``region_id``, ``material_region_id`` and either an explicit
    ``density``/``specific_heat`` MaterialProperty or a roster ``density_key``/
    ``specific_heat_key``; optional ``basis``, ``volume_method``, ``region_label`` and
    ``with_moments``. Deterministic: specs are matched to sorted interface labels in order and
    the finished payload is passed through ``i2_payload_validate_v1.validate_payload``.

    Optional blocks are emitted only when requested or present, never by default:

    * ``area_to_volume=True`` adds ``area.area_to_volume`` = ``{value, units: '1/mm'}`` from the
      total material-interface area over the total built volume;
    * a spec ``frame`` adds ``regions[].frame`` (``frame_v1.Frame`` or a plain frame dict);
    * a spec ``uncertainty`` adds ``regions[].uncertainty`` (``Uncertainty``/dict entries);
    * a spec ``thickness`` adds ``regions[].thickness_distribution`` (``ThicknessDistribution``
      or a dict).

    A malformed optional spec raises ``PayloadBuildError``; the required-field behaviour and the
    existing defaults are unchanged.
    """
    _require_interfaces(interfaces)
    if not isinstance(synthetic, bool):
        raise PayloadBuildError(f'synthetic: bool required, got {type(synthetic).__name__}')
    if source_sha256 is not None and not isinstance(source_sha256, dict):
        raise PayloadBuildError('source_sha256: dict of {filename: sha256} or None required')
    try:
        specs = list(region_specs)
    except TypeError as exc:
        raise PayloadBuildError('region_specs: iterable of dicts required') from exc
    if not specs:
        raise PayloadBuildError('region_specs: at least one region spec is required')
    for index, spec in enumerate(specs):
        if not isinstance(spec, dict):
            raise PayloadBuildError(
                f'region_specs[{index}]: dict required, got {type(spec).__name__}')

    labels = _resolve_labels(interfaces, specs)

    try:
        area_value = region_mass_v1.area_mm2(interfaces, area_kind)
    except (region_mass_v1.UnitError, region_mass_v1.UndeclaredError) as exc:
        raise PayloadBuildError(f'area_kind {area_kind!r}: {exc}') from exc
    method = _AREA_METHODS.get(area_kind, f'{area_kind} area, summed once')
    by_pair = {f'{a}-{b}': value for (a, b), value in
               sorted(region_mass_v1.interface_area_by_pair_mm2(interfaces).items())}

    states = []
    regions = []
    moments_blocks = {}
    for spec, label in zip(specs, labels):
        region_id = spec.get('region_id')
        material_region_id = spec.get('material_region_id')
        if not isinstance(region_id, str) or not region_id:
            raise PayloadBuildError('region_specs: region_id must be a nonempty string')
        if not isinstance(material_region_id, str) or not material_region_id:
            raise PayloadBuildError(
                f'region {region_id!r}: material_region_id must be a nonempty roster id')

        basis = spec.get('basis', DEFAULT_BASIS)
        if basis not in region_mass_v1.MASS_BASIS:
            raise PayloadBuildError(
                f'region {region_id!r}: basis {basis!r} not in {list(region_mass_v1.MASS_BASIS)}')
        volume_method = spec.get('volume_method', DEFAULT_VOLUME_METHOD)
        if not isinstance(volume_method, str) or not volume_method:
            raise PayloadBuildError(
                f'region {region_id!r}: volume_method must be a nonempty string')

        density = _resolve_material(spec, 'density')
        specific_heat = _resolve_material(spec, 'specific_heat')

        try:
            volume = region_mass_v1.region_boundary_volume_mm3(interfaces, label)
            geometry = region_mass_v1.Region(
                region_id, volume, 'mm3', volume_method,
                area_kind=area_kind, area_mm2=area_value)
            state = region_mass_v1.integrate_region(
                geometry, density=density, specific_heat=specific_heat, basis=basis)
        except (region_mass_v1.UnitError, ValueError) as exc:
            raise PayloadBuildError(f'region {region_id!r}: {exc}') from exc

        states.append(state)
        region = {
            'region_id': region_id,
            'material_region_id': material_region_id,
            'volume_mm3': state.geometry.volume_mm3,
            'volume_method': state.geometry.volume_method,
            'mass_basis': state.basis,
            'density': _material(state.density),
            'mass': _quantity(state.mass_kg),
            'specific_heat': _material(state.specific_heat),
            'heat_capacity': _quantity(state.heat_capacity_j_per_k),
            'missing': list(state.missing),
        }

        if 'frame' in spec:
            region['frame'] = _frame_block(spec['frame'], region_id)
        if 'uncertainty' in spec:
            region['uncertainty'] = _uncertainty_block(spec['uncertainty'], region_id)
        if 'thickness' in spec:
            region['thickness_distribution'] = _thickness_block(spec['thickness'], region_id)
        regions.append(region)

        if spec.get('with_moments', False):
            block = region_mass_v1.mass_moments_surface(
                interfaces['vertices'], _outward_faces(interfaces, label),
                density=state.density, region_id=region_id, basis=basis)
            moments_blocks[f'mass_moments_{region_id}'] = {
                'region_id': block.region_id,
                'basis': block.basis,
                'mass_kg': block.mass_kg,
                'centroid_mm': None if block.centroid_mm is None else block.centroid_mm.tolist(),
                'inertia_kg_mm2': (None if block.inertia_kg_mm2 is None
                                   else block.inertia_kg_mm2.tolist()),
                'units': {'mass': 'kg', 'centroid': 'mm', 'inertia': 'kg*mm2'},
                'inertia_SI_units': block.inertia_SI_units,
                'missing': list(block.missing),
            }

    try:
        totals = {
            'mass': _quantity(region_mass_v1.total_quantity(states, 'mass_kg')),
            'heat_capacity': _quantity(
                region_mass_v1.total_quantity(states, 'heat_capacity_j_per_k')),
        }
    except (region_mass_v1.UnitError, ValueError) as exc:
        raise PayloadBuildError(f'totals: {exc}') from exc

    area_block = {'kind': area_kind, 'value': area_value, 'units': AREA_UNITS,
                  'method': method, 'by_pair': by_pair}
    if area_to_volume:
        total_volume = float(math.fsum(state.geometry.volume_mm3 for state in states))
        if not math.isfinite(total_volume) or total_volume <= 0.0:
            raise PayloadBuildError(
                f'area_to_volume: total region volume must be positive, got {total_volume!r}')
        try:
            interface_area = region_mass_v1.area_mm2(interfaces, 'material_interface_area')
        except (region_mass_v1.UnitError, region_mass_v1.UndeclaredError) as exc:
            raise PayloadBuildError(f'area_to_volume: {exc}') from exc
        area_block['area_to_volume'] = {
            'value': float(interface_area) / total_volume, 'units': '1/mm'}

    payload = {
        'schema': SCHEMA_ID,
        'synthetic': synthetic,
        'area': area_block,
        'regions': regions,
        'totals': totals,
    }
    if source_sha256 is not None:
        payload['source_sha256'] = dict(source_sha256)
    payload.update(moments_blocks)

    errors = i2_payload_validate_v1.validate_payload(payload)
    if errors:
        raise PayloadBuildError('i2 payload invalid:\n- ' + '\n- '.join(errors))
    return payload


def build_and_check(interfaces, region_specs, **kwargs) -> dict:
    """``build_payload`` then require no registration issue; raise ``PayloadBuildError`` otherwise.

    ``require_known_density`` (default True) is consumed here; every other keyword is forwarded
    to ``build_payload`` (``area_kind``, ``synthetic``, ``source_sha256``, ``area_to_volume``).
    """
    require_known_density = kwargs.pop('require_known_density', True)
    payload = build_payload(interfaces, region_specs, **kwargs)
    issues = i2_region_registration_v1.registration_issues(
        payload, require_known_density=require_known_density)
    if issues:
        raise PayloadBuildError('i2 region registration invalid:\n- ' + '\n- '.join(issues))
    return payload
