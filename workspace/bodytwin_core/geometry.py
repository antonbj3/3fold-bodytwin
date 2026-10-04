"""N7c geometry manager and generic geometry exports (vendor/geomgr/export.py)."""
from .vendor.geomgr.api import GeometryManager
from .vendor.geomgr import core, export, identity, ops, certify2
from .config import paths
import numpy as np
from .geometry_priors import (prior_from_anthropometry, contralateral,
                              prior_from_tibia, scalar_population_prior)

def instantiate(subject='z001', scenario='S6', config=None):
    p = paths(config)
    gm = GeometryManager(p['geometry_data'], exclude_vsd=subject,
                         extra_dir=p['geometry_extra'], cache_dir='external_media')
    obs = np.load(p['geometry_data'] / 'vsd_obs.npz')
    si = list(map(str, obs['subj'])).index(subject)
    ri = int(np.flatnonzero(np.isfinite(obs['obs_lm'][si, :, 0, 0]))[0])
    names = {'S6': (core.X5, ['HC6', 'SGT', 'LT', 'MEC', 'LEC']),
             'S2': (core.X5, []), 'S4': ([], ['HC6', 'SGT', 'LT', 'MEC', 'LEC'])}
    fn, ln = names[scenario]
    feat = {k: float(obs['obs_feat'][si, ri, core.FEMUR_FEATURES.index(k)]) for k in fn}
    lm = {k: obs['obs_lm'][si, ri, core.PT_NAMES.index(k)] for k in ln}
    return gm.instantiate(features=feat, feature_units={k: core.feature_units(k) for k in feat},
                          landmarks=lm, subject=subject, scenario=scenario)
