"""geomgr.certify: validity certificate per operation (interpretation f).

Numeric part (numpy/scipy, runs on the cloud box):
  * input contract: known feature names, units, finite values, plausibility (population mean +- 6 SD);
  * statistical consistency: stage-A innovation chi2 (features vs population), stage-B landmark residual chi2
    after pose+shape fit, per-landmark outlier test (Bonferroni);
  * validity domain: posterior coefficient Mahalanobis vs chi2_r(0.999) -> extrapolation FLAG (not a rejection);
  * geometry: X1b rotation-invariant fold check against the population mean, positive signed volume.
Local part (needs trimesh/rtree and the workspace): BodyTwin geometry_certificate(evaluate=True) (closed,
winding-consistent, no self-intersection), BT-LV1 outward-rounded interval volume, C1 contract quantities.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from . import core as C

P_REJECT = 1e-3
SD_PLAUSIBLE = 6.0


class InputError(ValueError):
    pass


def check_inputs(model, feat_obs, feat_units, lm_obs, lm_units, known_points):
    """Contract checks before any numerics. Raises InputError (certificate: rejected, reason)."""
    for k, v in (feat_obs or {}).items():
        if k not in model.feat_names:
            raise InputError(f'unknown feature {k!r}')
        u = (feat_units or {}).get(k)
        if u != C.feature_units(k):
            raise InputError(f'feature {k}: units {u!r}, expected {C.feature_units(k)!r}')
        if not np.isfinite(v):
            raise InputError(f'feature {k}: non-finite value')
    for k, v in (lm_obs or {}).items():
        if k not in known_points:
            raise InputError(f'unknown landmark {k!r}')
        if lm_units != 'mm':
            raise InputError(f'landmark units {lm_units!r}, expected mm')
        if not np.all(np.isfinite(v)) or np.shape(v) != (3,):
            raise InputError(f'landmark {k}: non-finite or not 3D')


def plausibility(model, feat_obs):
    out = {}
    for k, v in (feat_obs or {}).items():
        j = model.feat_names.index(k)
        mu, sd = model.M[:, j].mean(), model.M[:, j].std(ddof=1)
        out[k] = float((v - mu) / sd)
    return out


def numeric_certificate(model, post, X, F, feat_obs=None):
    """Certificate for 'instantiate'. X: shape-frame vertices of the instance (model frame)."""
    checks = {}
    reasons = []
    pz = plausibility(model, feat_obs)
    checks['plausibility_z'] = pz
    if any(abs(z) > SD_PLAUSIBLE for z in pz.values()):
        reasons.append('feature outside population mean +- 6 SD')
    sa = post.stage_a
    checks['stage_a'] = dict(d2=sa['d2'], dof=sa['dof'], p=sa['p'])
    if sa['dof'] and sa['p'] < P_REJECT:
        reasons.append(f"features inconsistent with population (chi2 p={sa['p']:.2e})")
    if post.stage_b is not None:
        sb = post.stage_b
        L = len(sb['per_landmark_d2'])
        worst = max(sb['per_landmark_d2'].items(), key=lambda kv: kv[1])
        p_worst = float(stats.chi2.sf(worst[1], 3))
        checks['stage_b'] = dict(ssr=sb['ssr'], dof=sb['dof'], p=sb['p'], worst_landmark=worst[0],
                                 worst_p=p_worst, residual_rms_mm=sb['residual_rms_mm'],
                                 pose_identified=bool(post.pose_identified))
        if sb['p'] < P_REJECT:
            reasons.append(f"landmarks inconsistent with any rigidly placed population shape (chi2 p={sb['p']:.2e})")
        if p_worst < P_REJECT / L:
            reasons.append(f'landmark outlier {worst[0]} (p={p_worst:.2e})')
    d2b = float(np.sum(post.b ** 2 / model.lam))
    lim = float(stats.chi2.ppf(0.999, model.r))
    checks['extrapolation'] = dict(mahalanobis2=d2b, limit_chi2_0999=lim, inside=bool(d2b <= lim))
    fold = C.fold_check(model.mu, X, F)
    vol = C.signed_volume(X, F)
    checks['fold'] = fold
    checks['signed_volume_mm3'] = vol
    if fold['n_folded_faces'] > 0:
        reasons.append(f"{fold['n_folded_faces']} folded faces")
    if not vol > 0:
        reasons.append('non-positive volume')
    return dict(op='instantiate', accepted=not reasons, reasons=reasons, checks=checks,
                flags=[] if checks['extrapolation']['inside'] else ['outside population (extrapolation)'])


def geometry_certificate_local(V, F):
    """BodyTwin geometry_certificate(evaluate=True): closed, winding-consistent, positive volume, conforming,
    no self-intersection. Local only."""
    from ._local import bodytwin_geometry_certificate
    return bodytwin_geometry_certificate(V, F)


def interval_volume(V, F):
    """BT-LV1 outward-rounded binary64 interval of the enclosed volume (mm^3) and centroid (mm)."""
    from .vendor import interval_moments as im
    Vl = np.asarray(V, float).tolist()
    Fl = np.asarray(F, int).tolist()
    mom, ref = im.shell(Vl, Fl)
    vol = mom[0]
    fl = C.signed_volume(np.asarray(V, float), np.asarray(F))
    com = [(mom[1 + i] / vol) + im.IV(int(ref[i])) for i in range(3)] if vol.lo > 0 else None
    return dict(volume_mm3=[vol.lo, vol.hi], width_mm3=vol.width(), float_volume_mm3=fl,
                float_inside=bool(vol.lo <= fl <= vol.hi),
                centroid_mm=None if com is None else [[c.lo, c.hi] for c in com])
