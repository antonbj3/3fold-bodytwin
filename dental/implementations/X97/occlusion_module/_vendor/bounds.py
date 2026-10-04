import numpy as np

def support(t):
    ff = t['grid_faces']
    p = t['original_gap']
    r = t['reference_gap']
    return ff[np.isfinite(p[ff]).all(1) & np.isfinite(r[ff]).all(1)]
