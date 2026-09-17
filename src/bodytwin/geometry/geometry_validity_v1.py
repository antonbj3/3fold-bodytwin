"""Geometry validity certificate for region geometry: a surface is only accepted when it
is watertight, winding-consistent, positive-volume, conforming (no coincident vertices with
distinct indices) AND has an *evaluated* negative self-intersection result.

Key rule: an UNEVALUATED self-intersection check (``self_intersections=None``) is an explicit
abstention, NOT an acceptance. The geometry contract requires to
reject self-intersecting/open input rather than return a plausible area from it.

Location: src/bodytwin/geometry/geometry_validity_v1.py
It is additive, changes no existing symbol, and writes nothing.

Reuses the installed seams unchanged:
  * ``bodytwin.geometry.mesh_ingest_v1.surface_mesh`` for the strict mm/watertight/oriented
    positive-volume gate;
  * ``region_mass_v1.coincident_vertex_groups`` for the non-conforming assembly detection.
Pure NumPy/trimesh; no I/O and no writes.
"""
from __future__ import annotations

import numpy as np
import trimesh

from bodytwin.geometry.mesh_ingest_v1 import surface_mesh

from .region_mass_v1 import coincident_vertex_groups

__all__ = [
    'GeometryContractError',
    'evaluate_self_intersections',
    'geometry_certificate',
    'require_geometry_certificate',
]


class GeometryContractError(ValueError):
    """The geometry does not satisfy the validity contract (or the check was skipped)."""


def _reject_reason(cert):
    """Human-readable reason for a non-accepted certificate."""
    parts = []
    if not cert['watertight']:
        parts.append('surface is not watertight')
    if not cert['winding_consistent']:
        parts.append('winding is inconsistent')
    if not cert['positive_volume']:
        parts.append('volume is not finite/positive')
    if not cert['conforming']:
        parts.append('non-conforming assembled mesh: coincident vertices with distinct '
                     'indices')
    if cert['self_intersections'] is True:
        parts.append('self-intersections detected')
    elif cert['self_intersections'] is None:
        parts.append('self-intersection check unevaluated (None): abstaining, not accepting')
    return 'rejected: ' + '; '.join(parts) if parts else 'rejected'


# --- self-intersection evaluation -----------------------------------------------------------
# trimesh 5.1.0 ships no public self-intersection predicate (``Trimesh.intersection`` is a
# boolean-volume op requiring Manifold3D/Blender; ``trimesh.intersections`` only slices). The
# evaluator therefore prefers a public API when a future trimesh version provides one, and
# otherwise falls back to an exact triangle/triangle intersection test built only on public
# trimesh data (``Trimesh.triangles`` + ``trimesh.triangles.bounds_tree``) and NumPy. It never
# guesses: any failure or missing broad-phase backend yields ``None`` (abstention).

_TOL = 1e-9  # geometric tolerance for the exact tri/tri test (mm-scale coordinates)


def _public_self_intersection(mesh):
    """Return a bool from a public trimesh self-intersection API, else ``None``."""
    for name in ('is_self_intersecting', 'self_intersection'):
        attr = getattr(mesh, name, None)
        if attr is None:
            continue
        value = attr() if callable(attr) else attr
        return _as_intersection_bool(value)
    for name in ('mesh_self_intersection', 'self_intersection'):
        fn = getattr(trimesh.intersections, name, None)
        if fn is not None:
            return _as_intersection_bool(fn(mesh))
    return None


def _as_intersection_bool(value):
    """Normalize a public-API result to a bool (face-pair arrays count as intersections)."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    arr = np.asarray(value)
    if arr.dtype == bool:
        return bool(arr.any())
    return bool(arr.size > 0)


def _segment_triangle_hit(p0, p1, v0, v1, v2, tol=_TOL):
    """True when segment ``p0-p1`` meets triangle ``v0-v1-v2`` (Moller-Trumbore)."""
    edge1 = v1 - v0
    edge2 = v2 - v0
    pvec = np.cross(p1 - p0, edge2)
    det = edge1 @ pvec
    if abs(det) <= tol:
        return False  # parallel/coplanar; coplanar overlap is handled separately
    inv_det = 1.0 / det
    tvec = p0 - v0
    u = (tvec @ pvec) * inv_det
    if u < -tol or u > 1.0 + tol:
        return False
    qvec = np.cross(tvec, edge1)
    v = ((p1 - p0) @ qvec) * inv_det
    if v < -tol or u + v > 1.0 + tol:
        return False
    t = (edge2 @ qvec) * inv_det
    return -tol <= t <= 1.0 + tol


def _segment_segment_hit_2d(p0, p1, q0, q1, tol=1e-12):
    """True when 2-D segments ``p0-p1`` and ``q0-q1`` intersect."""
    r = p1 - p0
    s = q1 - q0
    den = r[0] * s[1] - r[1] * s[0]
    if abs(den) <= tol:
        return False
    qp = q0 - p0
    t = (qp[0] * s[1] - qp[1] * s[0]) / den
    u = (qp[0] * r[1] - qp[1] * r[0]) / den
    return (-tol <= t <= 1.0 + tol) and (-tol <= u <= 1.0 + tol)


def _point_in_triangle_2d(point, tri, tol=1e-12):
    """True when 2-D ``point`` lies inside or on triangle ``tri``."""

    def _side(a, b, c):
        return (a[0] - c[0]) * (b[1] - c[1]) - (b[0] - c[0]) * (a[1] - c[1])

    d1 = _side(point, tri[0], tri[1])
    d2 = _side(point, tri[1], tri[2])
    d3 = _side(point, tri[2], tri[0])
    has_neg = (d1 < -tol) or (d2 < -tol) or (d3 < -tol)
    has_pos = (d1 > tol) or (d2 > tol) or (d3 > tol)
    return not (has_neg and has_pos)


def _coplanar_triangles_overlap(tri1, tri2, normal):
    """True when two coplanar triangles overlap (edge crossing or containment)."""
    drop = int(np.argmax(np.abs(normal)))
    keep = [axis for axis in (0, 1, 2) if axis != drop]
    a = tri1[:, keep]
    b = tri2[:, keep]
    for i in range(3):
        p0, p1 = a[i], a[(i + 1) % 3]
        for j in range(3):
            q0, q1 = b[j], b[(j + 1) % 3]
            if _segment_segment_hit_2d(p0, p1, q0, q1):
                return True
    return _point_in_triangle_2d(a[0], b) or _point_in_triangle_2d(b[0], a)


def _triangles_intersect(tri1, tri2, tol=_TOL):
    """Exact triangle/triangle intersection test (handles coplanar overlap)."""
    a, b, c = tri1
    d, e, f = tri2
    n1 = np.cross(b - a, c - a)
    n2 = np.cross(e - d, f - d)
    n1_len = np.linalg.norm(n1)
    n2_len = np.linalg.norm(n2)
    if n1_len <= tol or n2_len <= tol:
        return False  # degenerate triangle is not a crossing
    dist2 = (np.array([d, e, f]) - a) @ n1 / n1_len
    if np.all(dist2 > tol) or np.all(dist2 < -tol):
        return False
    dist1 = (np.array([a, b, c]) - d) @ n2 / n2_len
    if np.all(dist1 > tol) or np.all(dist1 < -tol):
        return False
    if np.all(np.abs(dist2) <= tol) and np.all(np.abs(dist1) <= tol):
        return _coplanar_triangles_overlap(tri1, tri2, n1)
    for p0, p1 in ((a, b), (b, c), (c, a)):
        if _segment_triangle_hit(p0, p1, d, e, f):
            return True
    for p0, p1 in ((d, e), (e, f), (f, d)):
        if _segment_triangle_hit(p0, p1, a, b, c):
            return True
    return False


def evaluate_self_intersections(vertices, faces):
    """Evaluate whether a triangle surface intersects itself.

    Returns ``True`` when intersecting face pairs are found, ``False`` when the surface is
    evaluated clean, and ``None`` when no supported evaluation is available (missing public
    API *and* unusable broad-phase/data) or when any check raises. ``None`` is an abstention,
    never an acceptance.
    """
    try:
        mesh = trimesh.Trimesh(vertices, faces, process=False)
        public = _public_self_intersection(mesh)
        if public is not None:
            return public
        triangles = np.asarray(mesh.triangles, dtype=np.float64)
        n_faces = len(triangles)
        if n_faces < 2:
            return False
        tree = trimesh.triangles.bounds_tree(triangles)
        face_vertices = np.asarray(mesh.faces)
        candidates = set()
        for i in range(n_faces):
            bounds = np.concatenate([triangles[i].min(axis=0), triangles[i].max(axis=0)])
            for j in tree.intersection(tuple(bounds)):
                j = int(j)
                if j <= i:
                    continue
                # Faces sharing a vertex (index) are legitimately adjacent, not crossing.
                if set(face_vertices[i].tolist()) & set(face_vertices[j].tolist()):
                    continue
                candidates.add((i, j))
        for i, j in sorted(candidates):
            if _triangles_intersect(triangles[i], triangles[j]):
                return True
        return False
    except Exception:
        return None


def geometry_certificate(vertices, faces, *, tolerance_mm=1e-6, self_intersections=None,
                         declared_units='mm', evaluate=False):
    """Return an explicit validity certificate for a candidate mm surface.

    ``self_intersections`` is tri-state: True (found), False (evaluated, none found), or None
    (not evaluated). Only ``False`` can satisfy the acceptance contract; ``None`` abstains.

    When ``evaluate=True`` and no ``self_intersections`` declaration is supplied, the tri-state
    is filled by :func:`evaluate_self_intersections`; if that evaluation is unavailable it
    returns ``None`` and the certificate abstains rather than accepting. With ``evaluate=False``
    (the default) behaviour is unchanged and ``self_intersections`` passes through verbatim.

    On a rejecting ``surface_mesh`` ValueError the certificate is returned (not raised) with
    ``accepted=False``, ``reason=str(exc)`` and explicit False validity fields.
    """
    try:
        tolerance_ok = bool(np.isfinite(tolerance_mm)) and float(tolerance_mm) > 0.0
    except (TypeError, ValueError):
        tolerance_ok = False
    if not tolerance_ok:
        raise GeometryContractError(
            f'finite positive tolerance_mm required, got {tolerance_mm!r}')

    if evaluate and self_intersections is None:
        self_intersections = evaluate_self_intersections(vertices, faces)

    cert = {
        'units': 'mm',
        'watertight': False,
        'winding_consistent': False,
        'positive_volume': False,
        'conforming': False,
        'self_intersections': self_intersections,
        'volume_mm3': None,
        'tolerance_mm': float(tolerance_mm),
        'accepted': False,
        'reason': '',
    }

    # Conformity is determinable from the raw arrays even when surface_mesh refuses them.
    try:
        groups, _indices = coincident_vertex_groups(vertices, tolerance_mm=tolerance_mm)
        cert['conforming'] = (groups == 0)
    except ValueError:
        cert['conforming'] = False

    try:
        mesh = surface_mesh(vertices, faces, units=declared_units)
    except ValueError as exc:
        cert['reason'] = str(exc)
        return cert

    # surface_mesh enforces explicit mm coordinates, watertightness, consistent winding and
    # finite positive volume, so reaching here establishes each of those guarantees.
    cert['units'] = 'mm'
    cert['watertight'] = True
    cert['winding_consistent'] = True
    cert['positive_volume'] = True
    cert['volume_mm3'] = float(mesh.volume_mm3)

    cert['accepted'] = (
        cert['watertight']
        and cert['winding_consistent']
        and cert['positive_volume']
        and cert['conforming']
        and cert['self_intersections'] is False
    )
    cert['reason'] = ('accepted: closed, winding-consistent, positive-volume, conforming mm '
                      'surface with self_intersections=False'
                      if cert['accepted'] else _reject_reason(cert))
    return cert


def require_geometry_certificate(vertices, faces, **kwargs):
    """Return the certificate or raise GeometryContractError when it is not accepted."""
    cert = geometry_certificate(vertices, faces, **kwargs)
    if cert['accepted'] is not True:
        raise GeometryContractError(cert['reason'])
    return cert
