"""Region-bound geometry -> material -> mass and thermal capacity, with explicit units,
declared area kind/measurement method, bound provenance and abstention on missing data.

Additive only: it imports the installed seam ``bodytwin.geometry.mesh_ingest_v1``. Its
area-kind semantics were checked against an earlier area-descriptor prototype (optional
parity test in tests/test_region_mass_v1.py).

Location: src/bodytwin/geometry/region_mass_v1.py
No existing symbol is modified and nothing is written.

Rules enforced here (geometry and constitutive contract):
  * lengths are mm, volumes mm3, areas mm2, density kg/mm3, heat capacity J/K;
  * a missing density or specific heat yields UNKNOWN, never a silent default;
  * an undeclared functional exchange area is refused, not equated with a wall area;
  * a unit tag that is not in the contract raises instead of being guessed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from bodytwin.geometry.mesh_ingest_v1 import surface_mesh

PROVENANCE = ('measured', 'calibrated', 'literal_cited', 'assumption', 'UNKNOWN')
AREA_KINDS = ('exterior_area', 'material_interface_area', 'functional_exchange_area')

_LENGTH_TO_MM = {'mm': 1.0, 'cm': 10.0, 'm': 1000.0}
_VOLUME_TO_MM3 = {'mm3': 1.0, 'cm3': 1.0e3, 'm3': 1.0e9, 'L': 1.0e6, 'mL': 1.0e3}
_DENSITY_TO_KG_PER_MM3 = {'kg/L': 1.0e-6, 'g/cm3': 1.0e-6, 'g/mL': 1.0e-6,
                          'kg/m3': 1.0e-9, 'g/L': 1.0e-9}
_CP_TO_J_PER_KG_K = {'J/(kg*K)': 1.0, 'J/(g*K)': 1.0e3, 'kJ/(kg*K)': 1.0e3}


class UnitError(ValueError):
    """A declared unit is absent or outside the contract."""


class UndeclaredError(ValueError):
    """A descriptor kind needs an explicit caller declaration."""


class NonConformingError(ValueError):
    """Coincident boundaries reveal a non-conforming assembled mesh, not exterior area."""


# Active muscle mass, total muscle mass and total body mass must never be conflated. A region
# volume integral is its own basis. Scalar values such as the 17.2 kg active muscle mass used
# by the muscle_energetics and thermoregulation cells are separate: that value is
# model-relative (it comes from one reference producer run), not a literature anchor and not
# a fixed constant.
MASS_BASIS = ('region_volume', 'active_muscle', 'total_muscle', 'total_body')


def _convert(value, units, table, what):
    if units not in table:
        raise UnitError(f'{what}: unsupported unit {units!r}; contract units {sorted(table)}')
    v = float(value)
    if not np.isfinite(v) or v <= 0:
        raise UnitError(f'{what}: finite positive value required, got {value!r}')
    return v * table[units]


def to_mm(value, units):
    return _convert(value, units, _LENGTH_TO_MM, 'length')


def to_mm3(value, units):
    return _convert(value, units, _VOLUME_TO_MM3, 'volume')


def to_kg_per_mm3(value, units):
    return _convert(value, units, _DENSITY_TO_KG_PER_MM3, 'density')


def to_J_per_kgK(value, units):
    return _convert(value, units, _CP_TO_J_PER_KG_K, 'specific heat')


@dataclass(frozen=True)
class Quantity:
    name: str
    value: float | None
    units: str
    provenance: str
    missing: tuple[str, ...] = ()

    @property
    def known(self):
        return self.value is not None


@dataclass(frozen=True)
class MaterialProperty:
    name: str
    value: float
    units: str
    provenance: str
    source: str | None = None
    source_sha256: str | None = None

    def __post_init__(self):
        if type(self.name) is not str or not self.name:
            raise ValueError('material property needs a nonempty name')
        if self.provenance not in PROVENANCE:
            raise ValueError(f'provenance {self.provenance!r} not in {PROVENANCE}')
        if not np.isfinite(self.value) or self.value <= 0:
            raise ValueError('material property needs a finite positive value')

    def as_density_kg_per_mm3(self):
        return Quantity(self.name, to_kg_per_mm3(self.value, self.units), 'kg/mm3',
                        self.provenance)

    def as_cp_J_per_kgK(self):
        return Quantity(self.name, to_J_per_kgK(self.value, self.units), 'J/(kg*K)',
                        self.provenance)


# Seeded from an audited property registry of the repository (values are
# literals cited in the staged repo, not individual measurements). ``literal_cited`` means
# "repo literal with a source pointer"; the primary literature source (DOI/PMID) is NOT
# verified here, and species, age and hydration regime are undeclared. ``source`` is
# ``<repo path>::<symbol>`` and ``source_sha256`` is the sha256 of that file; both are checked
# by tests/test_cross_module_integration.py, so an edit of a cited file needs a deliberate re-pin. No tendon density
# exists in the repository, so MSK-TENDON is deliberately absent -> UNKNOWN.
DENSITY_REGISTRY = {
    'MSK-MUSCLE': MaterialProperty('density_muscle', 1.06, 'kg/L', 'literal_cited',
                                   'src/bodytwin/cells/energy/muscle_energetics.py::MUSCLE_DENSITY_KG_PER_L',
                                   'd08a611c6549d3de86304e6b4fa74525471d7f44650af5ec8a82c0f574e0ed0c'),
    'NEU-TISSUE': MaterialProperty('density_brain', 1.04, 'kg/L', 'literal_cited',
                                   'src/bodytwin/cells/energy/krogh_bmr_voidfloor_gate.py::tissues["brain"]["density"]',
                                   '688b25fc6d61a070188b4e57d509966d760d2ddc25448103d8c90bf8e7651b40'),
    'CARD-VESSEL': MaterialProperty('density_heart', 1.05, 'kg/L', 'literal_cited',
                                    'src/bodytwin/cells/energy/krogh_bmr_voidfloor_gate.py::tissues["heart"]["density"]',
                                    '688b25fc6d61a070188b4e57d509966d760d2ddc25448103d8c90bf8e7651b40'),
    'MSK-BONE': MaterialProperty('apparent_density_bone', 1.9, 'g/cm3', 'literal_cited',
                                 'src/bodytwin/cells/musculoskeletal/bone_calcium_flux_reconciliation.py::APPARENT_DENSITY_G_CM3',
                                 '0701e8eb8bf783ec75ed0a6b0a4081ccc8461c9a9e0d1004d7ab346f8c6422c7'),
    'GENERIC': MaterialProperty('density_default_generic', 1.05, 'kg/L', 'assumption',
                                'src/bodytwin/cells/energy/krogh_bmr_voidfloor_gate.py::M_from_VO2(density_kg_L)',
                                '688b25fc6d61a070188b4e57d509966d760d2ddc25448103d8c90bf8e7651b40'),
}
GENERIC_KEY = 'GENERIC'

# Only a whole-body lumped value is audited (audited = present in the repo with a pointer; the
# primary source is not verified); there is no per-tissue specific heat in either repo, so a
# tissue region must receive an explicit declared assumption or stay UNKNOWN.
SPECIFIC_HEAT_REGISTRY = {
    'WHOLE-BODY': MaterialProperty('specific_heat_body', 3490.0, 'J/(kg*K)', 'literal_cited',
                                   'src/bodytwin/cells/organ_systems/thermoregulation.py::SPECIFIC_HEAT_BODY_J_PER_KG_K',
                                   '4474bfbf7b2bec7209d2283b658ad573a88ef00d94843e7c47870f7d1387dbc8'),
}


def lookup_density(region_id):
    """Return the registered density or None; never substitute a generic default."""
    return DENSITY_REGISTRY.get(region_id)


def lookup_specific_heat(region_id):
    return SPECIFIC_HEAT_REGISTRY.get(region_id)


@dataclass(frozen=True)
class Region:
    region_id: str
    volume: float
    volume_units: str
    volume_method: str
    area_kind: str | None = None
    area_mm2: float | None = None
    source_sha256: str | None = None
    volume_mm3: float = field(init=False)

    def __post_init__(self):
        if type(self.region_id) is not str or not self.region_id:
            raise ValueError('region needs a nonempty region_id')
        if type(self.volume_method) is not str or not self.volume_method:
            raise ValueError('region needs a declared volume_method')
        object.__setattr__(self, 'volume_mm3', to_mm3(self.volume, self.volume_units))
        if self.area_kind is not None and self.area_kind not in AREA_KINDS:
            raise ValueError(f'area_kind {self.area_kind!r} not in {AREA_KINDS}')
        if self.area_mm2 is not None and (not np.isfinite(self.area_mm2) or self.area_mm2 < 0):
            raise ValueError('area_mm2 must be finite and nonnegative')


@dataclass(frozen=True)
class RegionState:
    region_id: str
    geometry: Region
    density: MaterialProperty | None
    specific_heat: MaterialProperty | None
    mass_kg: Quantity
    heat_capacity_j_per_k: Quantity
    basis: str = 'region_volume'

    def __post_init__(self):
        if self.basis not in MASS_BASIS:
            raise UnitError(f'mass basis {self.basis!r} not in {MASS_BASIS}')

    @property
    def missing(self):
        return tuple(sorted(set(self.mass_kg.missing) | set(self.heat_capacity_j_per_k.missing)))


def _mass_quantity(volume_mm3, density):
    if density is None:
        return Quantity('mass', None, 'kg', 'UNKNOWN', ('density',))
    v = to_mm3(volume_mm3, 'mm3')
    q = density.as_density_kg_per_mm3()
    return Quantity('mass', v * q.value, 'kg', q.provenance)


def _cp_quantity(mass_kg, specific_heat):
    if specific_heat is None:
        return Quantity('heat_capacity', None, 'J/K', 'UNKNOWN', ('specific_heat',))
    if mass_kg is None:
        return Quantity('heat_capacity', None, 'J/K', 'UNKNOWN', ('mass',))
    cp = specific_heat.as_cp_J_per_kgK()
    return Quantity('heat_capacity', mass_kg * cp.value, 'J/K', cp.provenance)


def integrate_region(region, *, density=None, specific_heat=None,
                     density_key=None, specific_heat_key=None, basis='region_volume'):
    """Bind one region to a material and integrate mass plus thermal capacity.

    Density/specific heat come from an explicit MaterialProperty or a registry key. Only an
    explicit ``density_key='GENERIC'`` may use the generic assumption, and it keeps the
    ``assumption`` provenance. ``basis`` names which mass the region integral represents
    (active muscle vs total muscle vs total body vs region volume); it is required metadata
    and is never inferred. Nothing is defaulted silently.
    """
    if basis not in MASS_BASIS:
        raise UnitError(f'mass basis {basis!r} not in {MASS_BASIS}')
    if density is None and density_key is not None:
        density = lookup_density(density_key)
    if specific_heat is None and specific_heat_key is not None:
        specific_heat = lookup_specific_heat(specific_heat_key)
    if density is not None and not isinstance(density, MaterialProperty):
        raise UnitError('density must be a MaterialProperty with explicit units')
    if specific_heat is not None and not isinstance(specific_heat, MaterialProperty):
        raise UnitError('specific heat must be a MaterialProperty with explicit units')
    mass = _mass_quantity(region.volume_mm3, density)
    heat = _cp_quantity(mass.value, specific_heat)
    return RegionState(region.region_id, region, density, specific_heat, mass, heat, basis)


def mass_fraction(numerator_kg, denominator_kg, *, numerator_basis, denominator_basis,
                  tol=1e-9):
    """Declared-basis mass fraction; refuses identical or inconsistent bases.

    No universal physiological fraction band is applied (an ad hoc active-muscle threshold
    is not a universal criterion). Only arithmetic consistency is enforced: the two bases must differ and the
    fraction must not exceed one.
    """
    for basis in (numerator_basis, denominator_basis):
        if basis not in MASS_BASIS:
            raise UnitError(f'mass basis {basis!r} not in {MASS_BASIS}')
    if numerator_basis == denominator_basis:
        raise UnitError('a mass fraction needs two distinct declared bases')
    n, d = float(numerator_kg), float(denominator_kg)
    if not np.isfinite(n) or not np.isfinite(d) or d <= 0 or n < 0:
        raise UnitError('finite nonnegative numerator and positive denominator required')
    value = n / d
    if value > 1.0 + tol:
        raise UnitError(f'fraction {value:.6g} exceeds 1: {numerator_basis!r} and '
                        f'{denominator_basis!r} are inconsistent')
    return Quantity('mass_fraction', value, '1', 'computed')


def total_quantity(states, attribute):
    """Strict partition total: known only when every region contributes the same basis."""
    if states and len({s.basis for s in states}) != 1:
        raise UnitError(f'mixed mass bases in one total: {sorted({s.basis for s in states})}')
    quantities = [getattr(s, attribute) for s in states]
    missing = tuple(s.region_id for s, q in zip(states, quantities) if not q.known)
    units = quantities[0].units if quantities else None
    if missing:
        return Quantity(attribute, None, units or '', 'UNKNOWN', missing)
    provenance = 'literal_cited' if all(q.provenance == 'literal_cited' for q in quantities) else 'assumption'
    return Quantity(attribute, float(sum(q.value for q in quantities)), units, provenance)


def _triangle_areas_mm2(vertices, faces):
    v = np.asarray(vertices, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    if v.ndim != 2 or v.shape[1] != 3 or f.ndim != 2 or f.shape[1] != 3:
        raise ValueError('expected (N,3) vertices and (M,3) faces')
    tri = v[f]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    return 0.5 * np.linalg.norm(n, axis=1)


def area_mm2(interfaces, kind='material_interface_area'):
    """Area of a label/tetra interface dict, counted once per shared face.

    Mirrors ``interface_area_mm2`` of the earlier area-descriptor prototype. ``functional_exchange_area`` is deliberately
    not reachable here; use ``functional_exchange_area_mm2`` with a declaration.
    """
    front = np.asarray(interfaces['front'])
    back = np.asarray(interfaces['back'])
    a = _triangle_areas_mm2(interfaces['vertices'], interfaces['faces'])
    if kind == 'exterior_area':
        mask = (front == 0) | (back == 0)
    elif kind == 'material_interface_area':
        mask = (front != 0) & (back != 0)
    elif kind == 'functional_exchange_area':
        raise UndeclaredError('declare the exchange pairs; a wall area is not exchange area')
    else:
        raise UnitError(f'unknown area kind {kind!r}; contract kinds {AREA_KINDS}')
    return float(a[mask].sum())


def functional_exchange_area_mm2(interfaces, *, declared_pairs=None):
    """Area of the caller-declared material pairs, or refuse.

    ``declared_pairs`` is a nonempty iterable of unordered {front, back} label pairs with a
    stated mechanism. Requiring it prevents an exchange term silently reusing a wall area.
    """
    if declared_pairs is None:
        raise UndeclaredError('functional exchange area must be declared, not assumed')
    pairs = {frozenset(p) for p in declared_pairs}
    if not pairs or any(len(p) != 2 or 0 in p for p in pairs):
        raise UndeclaredError('declared pairs must be nonempty {front, back} with nonzero labels')
    front = np.asarray(interfaces['front'])
    back = np.asarray(interfaces['back'])
    a = _triangle_areas_mm2(interfaces['vertices'], interfaces['faces'])
    mask = np.array([(f != 0 and b != 0 and frozenset((int(f), int(b))) in pairs)
                     for f, b in zip(front, back)], dtype=bool)
    return float(a[mask].sum())


def coincident_vertex_groups(vertices, *, tolerance_mm=1e-6):
    """Count coordinate groups holding more than one distinct vertex index.

    A conforming boundary has unique coordinates (the repo's seams weld or share indices).
    A non-conforming assembly of coincident vertices with distinct indices leaves the shared
    split plane reported as exterior by index-based extractors, which is the non-conforming failure.
    """
    if not np.isfinite(tolerance_mm) or tolerance_mm <= 0:
        raise UnitError('finite positive tolerance_mm required')
    v = np.asarray(vertices, dtype=np.float64)
    if v.ndim != 2 or v.shape[1] != 3 or not np.isfinite(v).all():
        raise ValueError('expected finite (N,3) vertices')
    rounded = np.round(v / tolerance_mm).astype(np.int64)
    _, counts = np.unique(rounded, axis=0, return_counts=True)
    duplicated = counts[counts >= 2]
    return int(len(duplicated)), int(duplicated.sum())


def assert_conforming_interfaces(interfaces, *, tolerance_mm=1e-6):
    """Refuse a plausible exterior area that hides a coincident internal boundary."""
    groups, vertices = coincident_vertex_groups(interfaces['vertices'], tolerance_mm=tolerance_mm)
    if groups:
        raise NonConformingError(
            f'{groups} coincident vertex groups ({vertices} distinct indices) at '
            f'tolerance_mm={tolerance_mm}: the assembled mesh is non-conforming, so the shared '
            'boundary is reported as exterior. Supply a conforming mesh or declare a coordinate '
            'weld before trusting the interface area.')
    return True


def material_interface_area_checked(interfaces, *, tolerance_mm=1e-6):
    """Material-material area after the conforming certificate; never a silent wall area."""
    assert_conforming_interfaces(interfaces, tolerance_mm=tolerance_mm)
    return area_mm2(interfaces, 'material_interface_area')


def closed_surface_volume_mm3(vertices, faces):
    """Region volume from a closed, consistently oriented mm surface (reuses mesh_ingest)."""
    return surface_mesh(vertices, faces, units='mm').volume_mm3


def region_boundary_volume_mm3(interfaces, region_label):
    """Volume of one label region from its outward boundary; shared faces counted once."""
    if region_label == 0:
        raise ValueError('label 0 is exterior, not a material region')
    faces = np.asarray(interfaces['faces'])
    front = np.asarray(interfaces['front'])
    back = np.asarray(interfaces['back'])
    outward = np.concatenate([faces[back == region_label], faces[front == region_label][:, ::-1]])
    return closed_surface_volume_mm3(interfaces['vertices'], outward)


def voxel_region_volume_mm3(labels, region_label, *, pitch_mm=1.0, pitch_declared=False):
    """Voxel-count volume; the repo convention is corners with origin 0 and pitch 1."""
    arr = np.asarray(labels)
    if arr.ndim != 3:
        raise ValueError('expected a 3-D label array')
    if pitch_mm != 1.0 and not pitch_declared:
        raise UnitError('non-unit pitch needs an explicit transform declaration')
    if not np.isfinite(pitch_mm) or pitch_mm <= 0:
        raise UnitError('finite positive pitch required')
    count = int(np.count_nonzero(arr == region_label))
    if count == 0:
        raise ValueError(f'region label {region_label} absent from the grid')
    return count * (pitch_mm ** 3)


def tetra_volume_mm3(vertices, cells):
    """Analytic tetrahedral-complex volume; no surface extraction, no trimesh."""
    v = np.asarray(vertices, dtype=np.float64)
    t = np.asarray(cells)
    if v.ndim != 2 or v.shape[1] != 3 or t.ndim != 2 or t.shape[1] != 4 or len(t) == 0:
        raise ValueError('expected (N,3) vertices and (M,4) cells')
    if t.min() < 0 or t.max() >= len(v):
        raise ValueError('cell vertex index outside vertices')
    p = v[t]
    six = np.einsum('ij,ij->i', p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
    if not np.isfinite(six).all():
        raise ValueError('nonfinite cell volume')
    return float(np.abs(six).sum() / 6.0)


def interface_area_by_pair_mm2(interfaces):
    """Material-material area per unordered {label, label} pair, each face counted once.

    Region mapping needs to know which two materials meet, not only the total; the pair sums
    equal ``area_mm2(..., 'material_interface_area')``.
    """
    front = np.asarray(interfaces['front'])
    back = np.asarray(interfaces['back'])
    a = _triangle_areas_mm2(interfaces['vertices'], interfaces['faces'])
    mask = (front != 0) & (back != 0)
    pairs = {}
    for f, b, area in zip(front[mask].tolist(), back[mask].tolist(), a[mask].tolist()):
        key = (min(f, b), max(f, b))
        pairs[key] = pairs.get(key, 0.0) + float(area)
    return pairs


def assert_unique_regions(states):
    """Refuse two states that claim the same region id: a region maps to one material."""
    ids = [s.region_id for s in states]
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate region_id in one binding: a region maps to one material')
    return True


@dataclass(frozen=True, eq=False)
class MassMoments:
    region_id: str
    mass_kg: float | None
    centroid_mm: np.ndarray | None
    inertia_kg_mm2: np.ndarray | None
    density: MaterialProperty | None
    basis: str = 'region_volume'
    missing: tuple[str, ...] = ()

    @property
    def known(self):
        return self.mass_kg is not None

    @property
    def inertia_SI_units(self):
        return 'kg*m2 (multiply kg*mm2 by 1e-6)'


def _assemble_moments(region_id, V, M1, M2, density, basis):
    if basis not in MASS_BASIS:
        raise UnitError(f'mass basis {basis!r} not in {MASS_BASIS}')
    if density is None:
        return MassMoments(region_id, None, None, None, None, basis, ('density',))
    if not np.isfinite(V) or V <= 0:
        raise ValueError('positive finite region volume required')
    rho = density.as_density_kg_per_mm3().value
    centroid = M1 / V
    covariance = M2 - V * np.outer(centroid, centroid)
    inertia = rho * (np.trace(covariance) * np.eye(3) - covariance)
    return MassMoments(region_id, rho * V, centroid, inertia, density, basis, ())


def mass_moments_surface(vertices, faces, *, density=None, region_id='', basis='region_volume'):
    """Mass, centroid and inertia tensor (kg*mm2) of a closed homogeneous mm region.

    Reuses the strict ``mesh_ingest`` seam (watertight, consistently oriented, positive
    volume). Missing density yields unknown moments, never a default.
    """
    mesh = surface_mesh(vertices, faces, units='mm')
    v = mesh.vertices_mm
    f = mesh.faces
    tri = v[f]
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    six = np.einsum('ij,ij->i', a, np.cross(b, c))
    V = float(six.sum() / 6.0)
    if abs(V - mesh.volume_mm3) > 1e-9 * max(1.0, abs(V)):
        raise ValueError('surface moment volume disagrees with the seam volume')
    u = a + b + c
    M1 = np.einsum('i,ij->j', six / 24.0, u)
    S = np.einsum('ij,ik->ijk', u, u) + np.einsum('ij,ik->ijk', a, a) \
        + np.einsum('ij,ik->ijk', b, b) + np.einsum('ij,ik->ijk', c, c)
    M2 = np.einsum('i,ijk->jk', six / 120.0, S)
    return _assemble_moments(region_id, V, M1, M2, density, basis)


def mass_moments_tetra(vertices, cells, *, density=None, region_id='', basis='region_volume'):
    """Mass, centroid and inertia tensor of a positive tetrahedral complex (analytic path)."""
    v = np.asarray(vertices, dtype=np.float64)
    t = np.asarray(cells)
    if v.ndim != 2 or v.shape[1] != 3 or t.ndim != 2 or t.shape[1] != 4 or len(t) == 0:
        raise ValueError('expected (N,3) vertices and (M,4) cells')
    if t.min() < 0 or t.max() >= len(v):
        raise ValueError('cell vertex index outside vertices')
    p = v[t]
    six = np.einsum('ij,ij->i', p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
    if not np.isfinite(six).all() or np.any(six == 0):
        raise ValueError('zero or nonfinite tetrahedron volume')
    V = float(six.sum() / 6.0)
    weights = np.repeat(six / 6.0, 4).reshape(-1, 4)
    M1 = (weights[:, :, None] * p).sum(axis=(0, 1)) / 4.0
    origin = p[:, 0]
    a, b, c = p[:, 1] - origin, p[:, 2] - origin, p[:, 3] - origin
    local = _local_tetra_second_moment(a, b, c, six / 6.0)
    shifted = (np.einsum('i,ij,ik->ijk', six / 6.0, origin, origin) + local
               + np.einsum('i,ij,ik->ijk', six / 6.0, origin, (a + b + c) / 4.0)
               + np.einsum('i,ij,ik->ijk', six / 6.0, (a + b + c) / 4.0, origin))
    M2 = shifted.sum(axis=0)
    return _assemble_moments(region_id, V, M1, M2, density, basis)


def _local_tetra_second_moment(a, b, c, volumes):
    u = a + b + c
    S = (np.einsum('ij,ik->ijk', u, u) + np.einsum('ij,ik->ijk', a, a)
         + np.einsum('ij,ik->ijk', b, b) + np.einsum('ij,ik->ijk', c, c))
    return (np.einsum('i,ijk->ijk', volumes / 20.0, S))


def anatomical_cross_section_mm2(volume_mm3, length_mm, *, length_units='mm'):
    length = to_mm(length_mm, length_units)
    return to_mm3(volume_mm3, 'mm3') / length


def physiological_cross_section_mm2(volume_mm3, fibre_length_mm, pennation_angle,
                                    *, angle_units='rad', length_units='mm'):
    """PCSA = V cos(theta) / L_fibre. Angle must be declared rad or deg."""
    length = to_mm(fibre_length_mm, length_units)
    angle = float(pennation_angle)
    if angle_units == 'deg':
        angle = np.deg2rad(angle)
    elif angle_units != 'rad':
        raise UnitError("pennation angle units must be 'rad' or 'deg'")
    return to_mm3(volume_mm3, 'mm3') * np.cos(angle) / length
