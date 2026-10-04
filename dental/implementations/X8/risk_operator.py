"""Geometric research operators. None returns clinical nerve-injury probability."""
import numpy as np
from scipy.special import gammaincc
from scipy.optimize import brentq
from numpy.polynomial.hermite import hermgauss

def moment_tail(t, mean, sd):
    """Conservative upper P(X>t), X>=0, given exact population moments.

    Sample moments are NOT exact population moments. Caller must retain this condition.
    """
    t = np.asarray(t, float)
    markov = np.divide(mean, t, out=np.ones_like(t), where=t > 0)
    cantelli = sd * sd / (sd * sd + np.maximum(t - mean, 0) ** 2)
    return np.where(t > mean, np.minimum(1, np.minimum(markov, cantelli)), 1.0)

def gamma_signed_survival(s, mean, sd):
    """P(R*U>s), R gamma(mean,sd), U uniform[-1,1] independent.

    For s>=0: 1/2*[Q(k,s/theta)-s/(theta*(k-1))*Q(k-1,s/theta)].
    The gamma/uniform law is an explicit closure, not recovered individual data.
    """
    k = (mean / sd) ** 2
    scale = sd * sd / mean
    if k <= 1:
        raise ValueError('This closed form requires gamma shape>1')
    s = np.asarray(s, float)
    v = np.abs(s) / scale
    p = 0.5 * (gammaincc(k, v) - v / (k - 1) * gammaincc(k - 1, v))
    p = np.clip(p, 0, 0.5)
    return np.where(s >= 0, p, 1 - p)

def scenario_probability(distance, profile, bias=0.059, sd_canal=np.hypot(0.06, 0.22), order=96):
    (nodes, weights) = hermgauss(order)
    d = np.asarray(distance, float)
    shifts = bias + np.sqrt(2) * sd_canal * nodes
    p = gamma_signed_survival(d[..., None] - 1 - shifts, *profile['apex_mean_sd_mm'])
    return p @ weights / np.sqrt(np.pi)

def scenario_margin(profile, x=0.01):
    return brentq(lambda d: float(scenario_probability(d, profile)) - x, 0, 40, xtol=1e-10)

def body_bound(distance, profile, boundary_budget=0.3, radius=2, extra_budget=0, n_allocations=1001):
    """Whole finite cylinder Hausdorff/union bound, no joint independence required.

    d_true>=d_nom-max(Rentry,Rapex)-2r*sin(angle/2)-B-extra.
    Supplied moments and B MUST be valid for the actual pose/anatomy. Here they
    are sensitivity conditions. `extra`=1.5 is a drill-path budget, not measured.
    """
    d = np.asarray(distance, float)
    z = np.linspace(0, 2, n_allocations)
    th = 2 * np.arcsin(np.clip(z / (2 * radius), 0, 1))
    (mu, sd) = np.radians(profile['angle_mean_sd_deg'])
    angular = moment_tail(th, mu, sd)
    slack = d[..., None] - 1 - boundary_budget - extra_budget - z
    p = moment_tail(slack, *profile['entry_mean_sd_mm']) + moment_tail(slack, *profile['apex_mean_sd_mm']) + angular
    return np.minimum(1, np.min(p, axis=-1))

def bound_margin(profile, x=0.01, boundary_budget=0.3, extra_budget=0):
    return brentq(lambda d: float(body_bound(d, profile, boundary_budget, extra_budget=extra_budget)) - x, 0, 100, xtol=1e-09)

def scalar_control(distance, profile, B=0.3, radius=2):
    """Independently enumerated scalar conventional union/Cantelli comparator."""

    def tail(t, mu, sd):
        if t <= mu:
            return 1.0
        return min(1.0, mu / t, sd ** 2 / (sd ** 2 + (t - mu) ** 2))
    best = 1.0
    for z in np.linspace(0, 2, 1001):
        theta = 2 * np.arcsin(z / (2 * radius))
        a = tail(theta, *np.radians(profile['angle_mean_sd_deg']))
        s = distance - 1 - B - z
        best = min(best, a + tail(s, *profile['entry_mean_sd_mm']) + tail(s, *profile['apex_mean_sd_mm']))
    return best
