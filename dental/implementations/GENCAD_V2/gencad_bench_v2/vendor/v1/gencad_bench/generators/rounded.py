"""Round-2 changed CAD operation; not a new material or physical model."""
from copy import deepcopy
import numpy as np
from ..geometry import THETA, SECTORS, rounded
from .baselines import parametric

def round_intaglio(profile):
    R = float(np.max(profile['radii_mm'][0]))
    H = max(float(profile['apex_mm']), R)
    straight = H - R
    zs = [0.0]
    radii = [R]
    if straight > 0:
        zs.append(straight)
        radii.append(R)
    for theta in THETA[1:]:
        zs.append(straight + R * np.sin(theta))
        radii.append(R * np.cos(theta))
    return dict(z_mm=rounded(zs), radii_mm=rounded(np.repeat(np.array(radii)[:, None], SECTORS, axis=1)), apex_mm=round(H, 6))

def generate(task):
    t = deepcopy(task)
    t['public']['intaglio_template'] = round_intaglio(t['public']['intaglio_template'])
    d = parametric(t)
    d['generator'] = 'rounded_access_IFU'
    return d

def scale_profile(profile, factor):
    return dict(z_mm=rounded(np.asarray(profile['z_mm']) * factor), radii_mm=rounded(np.asarray(profile['radii_mm']) * factor), apex_mm=round(profile['apex_mm'] * factor, 6))
