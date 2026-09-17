"""Frame binding for region mass moments.

Binds a region's mass moments (``region_mass_v1.MassMoments``) to a frame descriptor
(``frame_v1.Frame``), expressing the centroid and the inertia tensor in that frame. Without it a
mass-moments block cannot state the frame its numbers live in.
The returned ``'frame'`` field is that provenance tag.

Location: src/bodytwin/geometry/frame_region_bridge_v1.py
It is additive, changes no existing symbol, and is pure (no I/O, no writes, no
third-party imports beyond numpy, no imports beyond the sibling geometry modules).

Semantics (matching ``frame_v1``): a frame is an origin plus an orthonormal right-handed
basis. ``from_global`` maps shared-mm coordinates into the frame, and a tensor rotates as
``axes.T @ I @ axes``. The inertia returned by ``region_mass_v1`` is about the centroid, so
only the basis affects it; a pure origin translation leaves ``inertia_kg_mm2`` unchanged.
Missing inputs (``moments.known`` is False) yield ``None`` numbers and forward the missing
reasons; a missing datum is never a hard error and never silently substituted. No
physiological claim is made here.
"""
from __future__ import annotations

import numpy as np

from . import frame_v1, region_mass_v1

__all__ = ['FrameBindingError', 'moments_in_frame', 'bind_region_moments']

UNITS = {'mass': 'kg', 'centroid': 'mm', 'inertia': 'kg*mm2'}


class FrameBindingError(ValueError):
    """A frame cannot be bound to region moments (wrong type or malformed input)."""


def _require_frame(frame):
    if not isinstance(frame, frame_v1.Frame):
        raise FrameBindingError(
            f'frame must be a frame_v1.Frame, got {type(frame).__name__}')
    return frame


def moments_in_frame(moments, frame):
    """Express ``moments`` in ``frame`` as an additive mass-moments block.

    Returns a dict with keys ``region_id``, ``mass_kg``, ``centroid_mm`` (3-vector in
    ``frame``), ``inertia_kg_mm2`` (3x3 in ``frame``), ``frame`` (the descriptor), ``units``
    (exact contract units) and ``missing`` (forwarded reasons). When ``moments.known`` is
    False the three numeric fields are ``None`` and ``missing`` is nonempty; this never
    raises for missing data.
    """
    frame = _require_frame(frame)

    missing = list(moments.missing)
    if not moments.known:
        return {
            'region_id': moments.region_id,
            'mass_kg': None,
            'centroid_mm': None,
            'inertia_kg_mm2': None,
            'frame': frame,
            'units': dict(UNITS),
            'missing': missing,
        }

    centroid_global = np.asarray(moments.centroid_mm, dtype=np.float64).reshape(3)
    inertia_global = np.asarray(moments.inertia_kg_mm2, dtype=np.float64)

    centroid_local = frame_v1.from_global(
        frame, centroid_global.reshape(1, 3))[0]
    axes = np.asarray(frame.axes, dtype=np.float64)
    # Equivalent to frame_v1.transfer_tensors(Frame(), frame, inertia) (verified in tests).
    inertia_local = axes.T @ inertia_global @ axes

    return {
        'region_id': moments.region_id,
        'mass_kg': float(moments.mass_kg),
        'centroid_mm': centroid_local,
        'inertia_kg_mm2': inertia_local,
        'frame': frame,
        'units': dict(UNITS),
        'missing': missing,
    }


def bind_region_moments(region_id, vertices, faces, density, frame, *,
                        basis='region_volume'):
    """Compute surface mass moments for a region, then express them in ``frame``.

    ``density`` is a ``region_mass_v1.MaterialProperty`` or ``None``; ``None`` (or an
    unregistered material) abstains, producing ``None`` moments and ``missing`` reasons.
    """
    frame = _require_frame(frame)
    moments = region_mass_v1.mass_moments_surface(
        vertices, faces, density=density, region_id=region_id, basis=basis)
    return moments_in_frame(moments, frame)
