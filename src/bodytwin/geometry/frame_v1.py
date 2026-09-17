"""Minimal frame / coordinate descriptor for region geometry.

It is additive
plumbing intended for later force / optics coupling: it carries an origin, an orthonormal
right-handed basis and a voxel pitch, and provides pure coordinate transfers. It does NOT
implement a fibre / anisotropy field (out of scope here), and makes no
physiological claim.

Location: src/bodytwin/geometry/frame_v1.py
It is additive, changes no existing symbol, and is pure (no I/O, no writes, no imports
from the repo).

Conventions (matching the repo convention origin=0, pitch=1):
  * lengths are millimetres; ``units`` must be ``'mm'`` (the mm/m guard lives here);
  * ``axes`` is a 3x3 orthonormal, right-handed basis (det = +1 within 1e-9);
  * ``origin_mm`` is the local origin expressed in the shared mm frame;
  * ``local @ axes.T`` maps local coordinates to the shared mm frame and
    ``(global - origin) @ axes`` is its inverse;
  * any pitch other than 1.0 must be explicitly declared via ``pitch_declared=True``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    'FrameError',
    'Frame',
    'to_global',
    'from_global',
    'transfer_points',
    'transfer_vectors',
    'transfer_tensors',
]

_ORTHO_TOL = 1e-9


class FrameError(ValueError):
    """A frame descriptor is malformed or a transfer violates the declared basis."""


@dataclass
class Frame:
    """An orthonormal, right-handed mm frame with an origin and a declared pitch."""

    origin_mm: object = (0.0, 0.0, 0.0)
    axes: object = None  # default = 3x3 identity, applied in __post_init__
    pitch_mm: float = 1.0
    pitch_declared: bool = False
    units: str = 'mm'

    def __post_init__(self):
        if self.units != 'mm':
            raise FrameError(f"frame units must be 'mm', got {self.units!r} (mm/m guard)")

        try:
            origin = np.asarray(self.origin_mm, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise FrameError(
                f'origin_mm must be a finite length-3 vector, got {self.origin_mm!r}') from exc
        if origin.shape != (3,) or not np.all(np.isfinite(origin)):
            raise FrameError(
                f'origin_mm must be a finite length-3 vector, got {self.origin_mm!r}')

        raw_axes = np.eye(3) if self.axes is None else self.axes
        try:
            axes = np.array(raw_axes, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise FrameError(f'axes must be a 3x3 numeric matrix: {exc}') from exc
        if axes.shape != (3, 3):
            raise FrameError(f'axes must have shape (3, 3), got {axes.shape}')
        if not np.all(np.isfinite(axes)):
            raise FrameError('axes must be finite')

        ortho_err = float(np.max(np.abs(axes @ axes.T - np.eye(3))))
        if not np.isfinite(ortho_err) or ortho_err > _ORTHO_TOL:
            raise FrameError(
                'axes must be orthonormal '
                f'(max|axes@axes.T - I| = {ortho_err:.3e} > {_ORTHO_TOL:.0e})')

        det = float(np.linalg.det(axes))
        if not np.isfinite(det) or abs(det - 1.0) > _ORTHO_TOL:
            raise FrameError(
                'axes must be right-handed '
                f'(det = {det:.6g}, expected +1 within {_ORTHO_TOL:.0e})')

        try:
            pitch = float(self.pitch_mm)
        except (TypeError, ValueError) as exc:
            raise FrameError(
                f'pitch_mm must be a finite positive number, got {self.pitch_mm!r}') from exc
        if not np.isfinite(pitch) or pitch <= 0.0:
            raise FrameError(f'pitch_mm must be finite and positive, got {self.pitch_mm!r}')
        if pitch != 1.0 and not self.pitch_declared:
            raise FrameError(
                f'pitch_mm={pitch!r} != 1.0 requires pitch_declared=True '
                '(repo convention is origin=0, pitch=1)')

        axes.flags.writeable = False
        self.axes = axes
        self.pitch_mm = pitch


def _require_frame(frame, name):
    if not isinstance(frame, Frame):
        raise FrameError(f'{name} must be a Frame, got {type(frame).__name__}')
    return frame


def _as_vectors(values, name):
    try:
        arr = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise FrameError(f'{name} must be an (N, 3) numeric array: {exc}') from exc
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise FrameError(f'{name} must have shape (N, 3), got {arr.shape}')
    if not np.all(np.isfinite(arr)):
        raise FrameError(f'{name} must be finite')
    return arr


def _as_tensors(values, name):
    try:
        arr = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise FrameError(f'{name} must be an (N, 3, 3) numeric array: {exc}') from exc
    if arr.ndim != 3 or arr.shape[1:] != (3, 3):
        raise FrameError(f'{name} must have shape (N, 3, 3), got {arr.shape}')
    if not np.all(np.isfinite(arr)):
        raise FrameError(f'{name} must be finite')
    return arr


def to_global(frame, local_points):
    """Map ``local_points`` from ``frame``'s local coordinates to the shared mm frame."""
    frame = _require_frame(frame, 'frame')
    local = _as_vectors(local_points, 'local_points')
    origin = np.asarray(frame.origin_mm, dtype=np.float64)
    return origin + local @ frame.axes.T


def from_global(frame, global_points):
    """Inverse of :func:`to_global`: map shared-mm points into ``frame`` local coordinates."""
    frame = _require_frame(frame, 'frame')
    glob = _as_vectors(global_points, 'global_points')
    origin = np.asarray(frame.origin_mm, dtype=np.float64)
    return (glob - origin) @ frame.axes


def transfer_points(frame_a, frame_b, points):
    """Express points given in ``frame_a`` coordinates in ``frame_b`` coordinates."""
    frame_a = _require_frame(frame_a, 'frame_a')
    frame_b = _require_frame(frame_b, 'frame_b')
    pts = _as_vectors(points, 'points')
    return from_global(frame_b, to_global(frame_a, pts))


def transfer_vectors(frame_a, frame_b, vectors):
    """Rotate vectors (no translation) from ``frame_a`` coordinates into ``frame_b``."""
    frame_a = _require_frame(frame_a, 'frame_a')
    frame_b = _require_frame(frame_b, 'frame_b')
    vecs = _as_vectors(vectors, 'vectors')
    return vecs @ (frame_b.axes.T @ frame_a.axes).T


def transfer_tensors(frame_a, frame_b, tensors):
    """Rotate (N, 3, 3) tensors from ``frame_a`` into ``frame_b`` as ``R @ T @ R.T``.

    ``R = frame_b.axes.T @ frame_a.axes`` is the same rotation used for vectors, so a
    tensor and a vector transform consistently under the declared bases.
    """
    frame_a = _require_frame(frame_a, 'frame_a')
    frame_b = _require_frame(frame_b, 'frame_b')
    tens = _as_tensors(tensors, 'tensors')
    R = frame_b.axes.T @ frame_a.axes
    return R @ tens @ R.T
