"""geomgr.noise: MEASURED landmark noise (BT-B30) for averaged landmark observations (N7c).

N7a's per-landmark covariance C_lm (anatomical frame) is empirical for ONE rating: second moment of the VSD rating
deviations around each subject's 20-rating mean + landmark-definition noise (LOO consensus). When a user supplies
the MEAN of k ratings, the usual assumption is C/k. BT-B30 measured the variance components (VSD 19 x 5 raters
x 4 trials): intra-rater (trial) and inter-rater variance per landmark. The mean of n_raters x n_trials ratings has
rating variance  v_k = s_inter^2 / n_raters + s_intra^2 / (n_raters n_trials)  instead of (s^2)/k, so

    C_k = f_k * C_rating + C_definition,     f_k = v_k / (s_inter^2 + s_intra^2)

(definition noise does not average out). API:
    B = B30Noise.load(extra_dir / 'b30_femur_noise.json')
    f = B.factor('SGT', n_raters=1, n_trials=4)            # -> scalar
    C_lm_k = B.scale_C(gm, n_raters, n_trials)             # dict name -> 3x3, for gm.instantiate(..., C_lm override)
numpy only.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import core as C


class B30Noise:
    def __init__(self, d):
        self.d = d
        self.lm = d['landmarks']

    @classmethod
    def load(cls, path):
        return cls(json.loads(Path(path).read_text()))

    def components(self, name):
        names = C.HEAD6 if name == 'HC6' else (name,)
        intra = np.mean([self.lm[n]['intra_sd_3d_mm'] ** 2 for n in names])
        inter = np.mean([self.lm[n]['inter_sd_3d_mm'] ** 2 for n in names])
        return float(intra), float(inter)

    def factor(self, name, n_raters=1, n_trials=1):
        intra, inter = self.components(name)
        return (inter / n_raters + intra / (n_raters * n_trials)) / (inter + intra)


def rating_and_definition_cov(gm_obs, keep, name_index):
    """Split N7a's C_lm into the rating part and the definition part (same data as api.GeometryManager)."""
    O = gm_obs
    valid = np.isfinite(O['obs_lm'][:, :, 0, 0])
    D = np.array([O['dev_anat'][s, r, name_index] for s in keep for r in np.flatnonzero(valid[s])])
    Nd = O['ndef_all'][keep, name_index]
    return np.cov(D.T), Nd.T @ Nd / len(Nd)


def scale_C(gm, O, keep, n_raters, n_trials, model='b30', b30=None):
    """C_lm for the mean of n_raters x n_trials ratings. model: 'b30' (measured components), 'naive' (C/k),
    'single' (k = 1, N7a)."""
    out = {}
    k = n_raters * n_trials
    for j, n in enumerate(C.PT_NAMES):
        Cr, Cd = rating_and_definition_cov(O, keep, j)
        f = {'b30': (b30.factor(n, n_raters, n_trials) if b30 is not None else None), 'naive': 1.0 / k, 'single': 1.0}[model]
        out[n] = f * Cr + Cd
    return out
