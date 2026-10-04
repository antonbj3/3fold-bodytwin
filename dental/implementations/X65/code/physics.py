"""Legacy physical closure, vector roots and an independent scalar check."""
import importlib.util, json, math, functools
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('legacy', ROOT / 'inputs/run_k1b_eval_legacy.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid
PROFILE = max((json.loads(l) for l in (ROOT / 'inputs/thr_harm.jsonl').read_text().splitlines() if abs(json.loads(l)['geo']['rho'] - 0.08) < 1e-09), key=lambda r: r['level'])
(FAC, _) = legacy.path_factors(json.loads((ROOT / 'inputs/g_med_post.json').read_text()))

def vector_root(g, s0, HV, sqa, machine_R=0.1, strength=1.0):
    (g, s0, HV, sqa) = np.broadcast_arrays(g, s0, HV, sqa)
    if np.any(sqa <= 0):
        raise ValueError('invalid defect input')
    lo = np.ones(g.shape)
    hi = np.full(g.shape, 5000.0)
    base = np.minimum(1.43 * (HV + 120) / sqa ** (1 / 6), 1.6 * HV) * strength
    alpha = 0.226 + HV * 0.0001
    for _ in range(80):
        F = (lo + hi) * 0.5
        mx = s0 + g * F
        mn = s0 + g * machine_R * F
        R = np.maximum(np.divide(mn, mx, out=np.zeros_like(mx), where=mx != 0), -1.0)
        sw = base * np.maximum((1 - R) * 0.5, 0) ** alpha
        high = (mx > 0) & (0.5 * (mx - mn) > sw)
        hi = np.where(high, F, hi)
        lo = np.where(high, lo, F)
    return np.where(g > 0, (lo + hi) * 0.5, np.inf)

def scalar_root(g, s0, HV, sqa, machine_R=0.1, strength=1.0):

    def fn(F):
        mx = s0 + g * F
        mn = s0 + g * machine_R * F
        if mx <= 0:
            return -1.0
        local_R = max(mn / mx, -1.0)
        base = min(1.43 * (HV + 120) / sqa ** (1 / 6), 1.6 * HV) * strength
        return 0.5 * (mx - mn) - base * ((1 - local_R) * 0.5) ** (0.226 + HV * 0.0001)
    return brentq(fn, 1.0, 5000.0, xtol=1e-10)

def draws(rows, n=2000):
    rng = np.random.default_rng(20260923)
    out = {}
    for row in rows:
        D = float(row['D_mm'])
        band = (0.55, 0.8) if D <= 3 else (0.7, 1.0) if D <= 3.6 else (0.8, 1.1) if D <= 4.2 else (0.9, 1.2)
        mc = legacy.mat_class(row['material'], row['manufacturing'])
        hvband = legacy.HV_BAND[mc]
        u = rng.random((n, 6))
        h = 0.25 + 0.2 * u[:, 0]
        rb = band[0] + (band[1] - band[0]) * u[:, 1]
        rho = 0.04 + 0.08 * u[:, 2]
        Sz = 10 + 20 * u[:, 3]
        HV = hvband[0] + (hvband[1] - hvband[0]) * u[:, 4]
        mu = 0.12 + 0.38 * u[:, 5]
        out[row['id']] = {'h': h, 'rb': rb, 'rho': rho, 'sqa': 1.95 * Sz, 'HV': HV, 'Fp': legacy.torque_Nmm(row['abutment_torque']) / (0.16 * 0.4 + 0.58 * mu * 1.74 + mu * 2.05 / 2)}
    return out

@functools.lru_cache(maxsize=256)
def profile_means(sizes):
    return np.array([legacy.line_avg(np.array(PROFILE['depth_profile']['d_mm']), np.array(PROFILE['depth_profile']['sigma_over_peak']), x * 0.001) for x in sizes])

def forward(row, d, kt_scale=1.0, strength_scale=1.0, preload_scale=1.0, lever_scale=1.0, defect_scale=1.0, residual_MPa=0.0):
    D = float(row['D_mm'])
    l = float(row['lever_l_mm']) if row['lever_l_mm'][0].isdigit() else 11.0
    r = D / 2 - d['h']
    rb = d['rb']
    I = math.pi / 4 * (r ** 4 - rb ** 4)
    A = math.pi * (r ** 2 - rb ** 2)
    if np.any(I <= 0):
        raise ValueError('bore exceeds critical root')
    sc = (d['rho'] / 0.08) ** (-0.4)
    kb = PROFILE['Kt_bend'] * sc * kt_scale
    ka = PROFILE['Kt_axial'] * sc * kt_scale
    sqa = d['sqa'] * defect_scale
    avg = profile_means(tuple(sqa))
    forces = []
    ports = []
    for (z, f) in FAC.items():
        M = f['eta_M'] * math.sin(math.pi / 6) * (l * lever_scale - f['h_above_emb'])
        N = f['eta_A'] * math.cos(math.pi / 6)
        g = (kb * M * r / I - ka * N / A) * avg
        s0 = -ka * f['kappa'] * d['Fp'] * preload_scale / A * avg + residual_MPa
        forces.append(vector_root(g, s0, d['HV'], sqa, strength=strength_scale))
        ports.append((g, s0))
    return (np.min(forces, axis=0), ports)
