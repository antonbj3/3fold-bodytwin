"""Trivial LP certificates; general cases still execute the inherited HiGHS LP.

For sum(w_i |z_i-prior_i|), w_i>=0, z_i>=lower_i, the independent
minimum is max(prior,lower). If it also satisfies all coupling constraints,
it is a global optimum. This changes no objective or feasible set.
"""
import numpy as np
from legacy.generators import optimize as lp_optimize, submit
from legacy.checks import feasibility as lp_feasibility
from legacy.geometry import section_weights

def special_ok(t, z, inner):
    r = t['requirements']
    if t['family'] == 'implant_crown':
        a = np.asarray(t['implant_axis'])
        angle = np.degrees(np.arccos(np.clip(a[2] / np.linalg.norm(a), -1, 1)))
        if angle > r['channel_max_deg']:
            return False
    if t['family'] == 'lattice_onlay' and t['minimum_strut_mm'] < r['strut_mm']:
        return False
    if t['family'] == 'bridge3':
        for x in t['connector_x']:
            if section_weights(t['xy'], t['faces'], x) @ (z - inner) < r['connector_mm2']:
                return False
    return not t['A'].shape[0] or bool(np.all(t['A'] @ z <= t['obstacle_b']))

def optimize(t):
    inner = np.asarray(t['preparation_z']) + t['requirements']['film_min_mm']
    lower = inner + t['requirements']['wall_mm']
    z = np.maximum(t['prior'], lower)
    if np.all(t['weights'] >= 0) and special_ok(t, z, inner):
        return submit(t, z, inner, dict(operation='trivial global L1 LP certificate'))
    return lp_optimize(t)

def feasibility(t):
    inner = np.asarray(t['preparation_z']) + t['requirements']['film_min_mm']
    lower = inner + t['requirements']['wall_mm']
    if special_ok(t, lower, inner):
        return dict(status='FEASIBLE', witness='explicit lower-bound vector satisfies all locked inequalities', z=lower.tolist())
    return lp_feasibility(t)
