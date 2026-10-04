"""R9 add linear axial membrane strips, without external supports or target fit."""
import numpy as np
from scipy.linalg import eigh
from whole_arch import build as bending_build
from whole_arch_refined import solve_floating
from contact_model import thickness, skew

def build(arch, records):
    s = bending_build(arch, records)
    K = s['K'].copy()
    nd = len(K)
    count = 0

    def displacement(j, k):
        p = s['patch'][5 * j + k]
        n = p['normal']
        r = p['lever']
        P = np.eye(3) - np.outer(n, n)
        D = np.zeros((3, nd))
        D[:, 11 * j:11 * j + 3] = P
        D[:, 11 * j + 3:11 * j + 6] = -P @ skew(r)
        D[:, 11 * j + 6 + k] = n
        return D
    regions = ['Buccogingival', 'Buccal', 'Incisal/occlusal', 'Palatal', 'Palatogingival']
    from contact_model import REGIONS
    for j in range(len(s['ids']) - 1):
        H0 = arch['teeth'][str(s['ids'][j])]['crown_height_mm']
        H1 = arch['teeth'][str(s['ids'][j + 1])]['crown_height_mm']
        for (k, region) in enumerate(REGIONS):
            p0 = s['patch'][5 * j + k]
            p1 = s['patch'][5 * (j + 1) + k]
            delta = p1['point'] - p0['point']
            L = np.linalg.norm(delta)
            if L < 1e-06:
                raise ValueError('Coincident membrane endpoints')
            t = delta / L
            v = t @ (displacement(j + 1, k) - displacement(j, k))
            h = thickness(region, records, 'regional', 1)
            w = (H0 + H1) / 10
            K += 2746 * h * w / L * np.outer(v, v)
            count += 1
    (vals, V) = eigh(K)
    positive = vals > vals[-1] * 1e-10
    N = V[:, ~positive]
    P = V[:, positive] / vals[positive] @ V[:, positive].T
    s.update(K=K, N=N, P=P, W=s['A'] @ P @ s['A'].T, C=N.T @ s['A'].T, null_modes=N.shape[1], minimum_eigenvalue=float(vals[0]), membrane_edges=count)
    return s
