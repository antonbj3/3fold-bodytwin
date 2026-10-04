"""Paired observation likelihood, joint (physical ratio, scanner bias) state.
Exact finite Gaussian arithmetic for declared additive observation closure.
No physical additive-error enclosure is claimed.
"""
import math
import numpy as np
from scipy.stats import norm

class PairedPort:
    A = np.array([[1.0, 1.0], [1.0, 0.0]])

    def __init__(self, mean=(0.8, 0.0), sd=(0.02, 0.02)):
        self.prec = np.diag(1 / np.array(sd) ** 2)
        self.h = self.prec @ np.array(mean)
        self.n = 0

    def update(self, scan, gauge, scan_sd, gauge_sd, corr):
        if not -1 < corr < 1 or min(scan_sd, gauge_sd) <= 0:
            raise ValueError('Paired metrology covariance not positive definite')
        cov = np.array([[scan_sd ** 2, corr * scan_sd * gauge_sd], [corr * scan_sd * gauge_sd, gauge_sd ** 2]])
        V = np.linalg.inv(cov)
        self.prec += self.A.T @ V @ self.A
        self.h += self.A.T @ V @ np.array([scan, gauge])
        self.n += 1

    def moments(self):
        C = np.linalg.inv(self.prec)
        return (C @ self.h, C)

    def report(self):
        (mu, C) = self.moments()
        q = norm.ppf(0.975)
        out = dict(observations=self.n, mean=mu.tolist(), covariance=C.tolist(), resolution='POPULATION', status='UPDATED_PAIRED_CONDITIONAL' if self.n else 'UNKNOWN_NO_PAIRED_MEASUREMENT')
        for (i, k) in enumerate(('sinter_factor', 'scanner_bias')):
            sd = math.sqrt(C[i, i])
            out[k] = dict(mean=float(mu[i]), sd=sd, credible_95=[float(mu[i] - q * sd), float(mu[i] + q * sd)], unit='post/pre ratio', observations=self.n, resolution='POPULATION', status=out['status'])
        return out
