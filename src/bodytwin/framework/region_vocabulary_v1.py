"""Framework-side copy of the region vocabulary (provenance + parity check).

``bodytwin.framework`` owns the result *envelope*; ``bodytwin.geometry`` owns the region
*payload schema*. This module only mirrors the region ids the material roster registers so
that the consumer preflight can refuse a payload that declares a region id that was never
registered (``unknown region_id``), instead of trusting it.

Provenance (I2 source in this repository, read-only)::

    src/bodytwin/geometry/material_roster_v1.py
    sha256 1998eec14cd3195f35705fb3d08324a765b02e23e9fde31a0fabbfb69e08ae2a
    REGION_IDS == ('CARD-VESSEL', 'GENERIC', 'MSK-BONE', 'MSK-MUSCLE',
                   'MSK-TENDON', 'NEU-TISSUE', 'WHOLE-BODY')
    published alongside ``examples/geometry/region_example.json``.

The roster is a single frozen vocabulary used for **two** distinct concepts, so it is mirrored under
two names rather than keeping a second list that could drift:

* :data:`I2_REGION_IDS` -- the payload's ``regions[].region_id`` (a geometry label, which may be a
  synthetic id such as ``SYNTH-BOX-1``); and
* :data:`I2_MATERIAL_REGION_IDS` -- the payload's per-region ``material_region_id`` (the explicit
  region -> material binding, whose values are the roster's material regions, e.g. ``MSK-MUSCLE``).

This module does **not** own or modify the payload schema: this file is a read-only mirror plus a parity check.
The roster is imported lazily from :mod:`bodytwin.geometry.material_roster_v1` (its import chain
pulls in ``numpy``/``trimesh``), so importing this module stays cheap for the consumer cells.

``MSK-LOWERLIMB-HIP`` is deliberately **not** part of this vocabulary: it is a placeholder that
the roster never registered (see :mod:`bodytwin.framework.consumer_preflight_v1`).
"""
from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path

# ------------------------------------------------------------------ I2 source in this repo -
I2_ROSTER_MODULE = "bodytwin.geometry.material_roster_v1"
I2_ROSTER_PATH = Path(__file__).resolve().parents[1] / "geometry" / "material_roster_v1.py"
I2_ROSTER_SHA256 = "1998eec14cd3195f35705fb3d08324a765b02e23e9fde31a0fabbfb69e08ae2a"

# Framework-side declared copy of ``material_roster_v1.REGION_IDS`` (kept in lock-step by
# :func:`assert_i2_parity`; a mismatch is a failing test, never a silent divergence).
I2_REGION_IDS = (
    "CARD-VESSEL",
    "GENERIC",
    "MSK-BONE",
    "MSK-MUSCLE",
    "MSK-TENDON",
    "NEU-TISSUE",
    "WHOLE-BODY",
)

# The same frozen roster is also the *material-region* vocabulary. An I2 region payload
# carries a per-region ``material_region_id`` (an explicit region -> material binding, e.g.
# ``SYNTH-BOX-1 -> MSK-MUSCLE``); that id must be a member of the roster. It is deliberately the
# SAME declared tuple object (not a second parallel list) so the two can never drift: the roster
# *is* the material-region vocabulary.
I2_MATERIAL_REGION_IDS = I2_REGION_IDS

# Unique module name for an explicit by-path import (never installed, never cached across calls).
_I2_ROSTER_BY_PATH_NAME = "_i1_by_path_material_roster_v1"


def load_i2_region_ids(roster_path=None):
    """Return I2's ``material_roster_v1.REGION_IDS`` tuple, or ``None`` when not importable.

    With ``roster_path=None`` the in-repo module :data:`I2_ROSTER_MODULE` is imported. An explicit
    ``roster_path`` is loaded by path (used by tests for the missing/broken-path case); a missing
    file or an ``ImportError`` inside the module returns ``None`` so a caller can degrade
    gracefully. Nothing is written.
    """
    if roster_path is None:
        try:
            module = importlib.import_module(I2_ROSTER_MODULE)
        except ImportError:
            return None
    else:
        path = Path(roster_path)
        if not path.exists():
            return None
        spec = importlib.util.spec_from_file_location(_I2_ROSTER_BY_PATH_NAME, path)
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
        except ImportError:
            return None

    region_ids = getattr(module, "REGION_IDS", None)
    if region_ids is None:
        return None
    return tuple(region_ids)


def assert_i2_parity(roster_path=None):
    """Return ``(loaded == declared, message)``; ``(None, ...)`` when I2 is not importable.

    A ``None`` status means the roster could not be imported; a ``False`` status is a real
    divergence between the declared ids and the material roster.
    """
    loaded = load_i2_region_ids(roster_path)
    if loaded is None:
        return (None, "I2 roster not importable")
    if loaded == I2_REGION_IDS:
        return (
            True,
            f"I1 declared ids match I2 material_roster_v1.REGION_IDS (sha256 {I2_ROSTER_SHA256})",
        )
    return (
        False,
        f"I2 region ids diverged: loaded {loaded!r} != declared {I2_REGION_IDS!r}",
    )
