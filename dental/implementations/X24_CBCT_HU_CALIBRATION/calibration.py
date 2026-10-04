"""X24 region-level reference calibration. H_ref is empirical reference CT#, not universally true HU."""
import numpy as np

def fit_affine(gray, href):
    g = np.asarray(gray, float)
    h = np.asarray(href, float)
    if len(g) != 2 or not np.isfinite(g).all() or (not np.isfinite(h).all()) or (g[1] <= g[0]) or (h[1] <= h[0]):
        raise ValueError('ordered finite distinct anchors required')
    b = (h[1] - h[0]) / (g[1] - g[0])
    a = h[0] - b * g[0]
    return np.array([a, b])

def apply(c, g):
    return np.polynomial.polynomial.polyval(np.asarray(g, float), c)

def bone_branch_b(h):
    h = np.asarray(h, float)
    app = 2 * (0.001 * h + 0.19)
    ash = 0.598 * np.maximum(app, 0)
    return {'rho_app': app, 'rho_ash': ash, 'E_MPa': np.maximum(1, 10500 * ash ** 2.57), 'density_valid': app > 0, 'published_NewTom_range': (h >= -744.4) & (h <= 230), 'Keller_ash_range': (ash >= 0.03) & (ash <= 1.22)}
