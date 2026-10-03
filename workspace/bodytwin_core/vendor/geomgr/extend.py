"Extend a geometry registry with tissue regions, footprints and source-supplied landmarks."
from __future__ import annotations

import numpy as np

from . import core as C
from . import identity as I
from . import measures as MS
from . import wrap as W

FOOT_R = 8.0   # N7b e_regions FOOT_R (assumed footprint radius)


def tissue_regions(TV, reg):
    """results/N7b/e_regions.py regions() (copied logic): head articular cap, neck, diaphysis (25-75 % L_mech),
    distal epiphysis (< 12 %), muscle footprints (attachment +- 8 mm). On the template (vertex ids are shared)."""
    vals = I.evaluate(reg, TV)
    L = MS.landmarks_of(vals, reg.part)
    m = C.measures(L, TV)
    head_c, knee_c, mech = m['_head_c'], m['_knee_c'], m['_mech']
    r_head = m['D_head'] / 2
    neck_c = np.mean([L[k] for k in MS.NECK4], 0)
    nax = C.unit(head_c - neck_c)
    h = (TV - knee_c) @ mech / m['L_mech']
    dist_head = np.abs(np.linalg.norm(TV - head_c, axis=1) - r_head)
    beyond = (TV - neck_c) @ nax > 0.35 * np.linalg.norm(head_c - neck_c)
    R = {}
    R['head_articular'] = np.flatnonzero((dist_head < 2.0) & beyond & ((TV - head_c) @ nax > -0.2 * r_head))
    t = np.clip((TV - neck_c) @ nax, -15, np.linalg.norm(head_c - neck_c))
    d_ax = np.linalg.norm(TV - (neck_c + np.outer(t, nax)), axis=1)
    R['neck'] = np.setdiff1d(np.flatnonzero((d_ax < 20) & ((TV - neck_c) @ nax > -15) & ~beyond), R['head_articular'])
    R['diaphysis'] = np.flatnonzero((h > 0.25) & (h < 0.75))
    R['distal_epiphysis'] = np.flatnonzero(h < 0.12)
    foot = {}
    att = [e for e in reg.by_kind('attachment')]
    muscles = sorted({e.name.split('|')[0].rsplit('_', 1)[0] for e in att})
    for mus in muscles:
        pts = np.array([vals[e.base] for e in att if e.name.split('|')[0].rsplit('_', 1)[0] == mus])
        d = np.min(np.linalg.norm(TV[:, None, :] - pts[None], axis=2), axis=1)
        foot[mus] = np.flatnonzero(d < FOOT_R)
    return R, foot


def extend_registry(reg, TV, tlem_native_V, tlem_points_native):
    """Return (new registry, diagnostics). tlem_points_native: {construct name: xyz in the TLEM member's native
    frame of the population}."""
    import copy
    r = copy.deepcopy(reg)
    r.version_note = (reg.version_note + ' | N7c: + N7b tissue regions/footprints, '
                      'TLEM femur constructs')
    R, foot = tissue_regions(TV, reg)
    for k, idx in R.items():
        r.add(I.Entity(r.part, 'region', f'tissue.{k}', indices=[int(i) for i in idx], source='derived',
                       note='results/N7b/e_regions.py regions() on the template'))
    for k, idx in foot.items():
        if len(idx):
            r.add(I.Entity(r.part, 'region', f'footprint.{k}', indices=[int(i) for i in idx], source='atlas_internal',
                           license='internal', note=f'TLEM attachments of {k} +- {FOOT_R} mm (radius assumed; N7b e)'))
    diag = {}
    names = sorted(tlem_points_native)
    ents, cp, Wt = W.construct_entities(r.part, tlem_native_V, np.array([tlem_points_native[n] for n in names]), names,
                                        'TLEM 2.0 table point carried by TPS from the TLEM femur in the population')
    for e in ents:
        r.add(e)
    diag['tlem_constructs'] = dict(n=len(names), affine_reproduction_err_mm=float(
        np.abs(Wt @ np.asarray(tlem_native_V)[cp] - np.array([tlem_points_native[n] for n in names])).max()))
    diag['n_entities'] = len(r.entities)
    diag['manifest_sha256'] = r.manifest_hash()
    return r, diag
