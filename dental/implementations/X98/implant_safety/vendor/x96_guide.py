import numpy as np
from scipy.optimize import brentq

def tail(t, mu, sd):
    t = np.asarray(t, dtype=float)
    return np.where(t > mu, np.minimum(1, np.minimum(mu / np.maximum(t, 1e-300), sd * sd / (sd * sd + np.maximum(t - mu, 0) ** 2))), 1.0)

def bound(clearance, p, h, r=2.0):
    d = np.atleast_1d(clearance).astype(float)
    rad = r + h
    z = np.linspace(0, 2 * rad, 2001)
    theta = 2 * np.arcsin(np.clip(z / (2 * rad), 0, 1))
    q = tail(theta, *np.radians(p['angle_mean_sd_deg']))
    c = d[:, None] - 2 - z
    return np.minimum(1, (tail(c, *p['entry_mean_sd_mm']) + tail(c, *p['apex_mean_sd_mm']) + q).min(1))

def displacement95(p, h):
    return float(brentq(lambda d: float(bound(d, p, h)[0]) - 0.05, 2, 100, xtol=1e-10) - 2)

def scalar_control(d, p, h):

    def t(x, m, s):
        return 1.0 if x <= m else min(1.0, m / x, s * s / (s * s + (x - m) ** 2))
    best = 1.0
    for z in np.linspace(0, 2 * (2 + h), 2001):
        th = 2 * np.arcsin(z / (2 * (2 + h)))
        c = d - 2 - z
        best = min(best, t(c, *p['entry_mean_sd_mm']) + t(c, *p['apex_mean_sd_mm']) + t(th, *np.radians(p['angle_mean_sd_deg'])))
    return best
