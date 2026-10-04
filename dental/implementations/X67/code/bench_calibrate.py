"""Integral energy calibration for an explicitly declared two-volume observer.

This is a lab prototype, not a validated drilling temperature field. Capacity bounds
must come from independent volume/density/calorimetry, not this same temperature fit.
"""
import csv, argparse, numpy as np
from scipy.optimize import lsq_linear
COLS = ['time_s', 'T1_C', 'T2_C', 'coolant_C', 'heater_power_W', 'mechanical_power_W']

def read_csv(path):
    with open(path) as f:
        r = csv.DictReader(f)
        if r.fieldnames != COLS:
            raise ValueError('CSV units/schema mismatch; expected ' + ','.join(COLS))
        a = np.array([[float(row[k]) for k in COLS] for row in r])
    if a.ndim != 2 or len(a) < 3 or (not np.isfinite(a).all()) or np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError('Finite strictly increasing time required')
    if np.any(a[:, 4:] < 0):
        raise ValueError('Negative power forbidden')
    return a

def check_metadata(cal_meta, drill_meta):
    for k in ['specimen_id', 'material_id', 'geometry_id', 'region_ids', 'observer_model']:
        if cal_meta.get(k) != drill_meta.get(k):
            raise ValueError('Calibration/drilling support mismatch: ' + k)
    if cal_meta.get('observer_model') != 'two_volume_temperatures':
        raise ValueError('Point sensors require a separately calibrated observation operator')
    if not cal_meta.get('known_electrical_power'):
        raise ValueError('Known heater power required to break scale ambiguity')
    if cal_meta.get('capacity_bounds_provenance') not in ['independent_measurement', 'our_own_fixture']:
        raise ValueError('Independent capacity-bounds provenance missing')
    if drill_meta.get('mechanical_work_provenance') not in ['independent_measurement', 'our_own_fixture']:
        raise ValueError('Measured mechanical work required; motor torque limit is not measured drilling torque')

def intervals(a):
    t = a[:, 0]
    dt = np.diff(t)
    x1 = a[:, 1] - a[:, 3]
    x2 = a[:, 2] - a[:, 3]
    if not np.all(a[:, 3] == a[0, 3]):
        raise ValueError('Changing coolant baseline needs expanded energy equation')
    return (dt, np.diff(x1), np.diff(x2), 0.5 * (x1[:-1] + x1[1:]) * dt, 0.5 * (x2[:-1] + x2[1:]) * dt)

def calibrate(a, capacity_bounds):
    (dt, d1, d2, i1, i2) = intervals(a)
    diff = i1 - i2
    n = len(dt)
    X = np.zeros((2 * n, 5))
    Y = np.zeros(2 * n)
    X[0::2, 0] = d1
    X[1::2, 1] = d2
    X[0::2, 2] = diff
    X[1::2, 2] = -diff
    X[0::2, 3] = i1
    X[1::2, 4] = i2
    Y[0::2] = a[:-1, 4] * dt
    lo = np.array([capacity_bounds[0][0], capacity_bounds[1][0], 0.0, 0.0, 0.0])
    hi = np.array([capacity_bounds[0][1], capacity_bounds[1][1], np.inf, np.inf, np.inf])
    fit = lsq_linear(X, Y, bounds=(lo, hi), tol=1e-12, lsmr_tol=1e-12, max_iter=200)
    residual = float(np.max(np.abs(X @ fit.x - Y)))
    unconstrained = np.linalg.lstsq(X, Y, rcond=None)[0]
    if not fit.success or residual > 0.001:
        raise ValueError('Calibration energy residual exceeded frozen 1e-3 J: ' + str(residual))
    if not np.all((unconstrained[:2] >= lo[:2]) & (unconstrained[:2] <= hi[:2])):
        raise ValueError('Independent capacity interval refutes fitted heater scale')
    return dict(parameters={k: float(v) for (k, v) in zip(['C1_J_K', 'C2_J_K', 'G_W_K', 'H1_W_K', 'H2_W_K'], fit.x)}, max_interval_residual_J=residual, rank=int(np.linalg.matrix_rank(X)), condition_number=float(np.linalg.cond(X)), guarantee='Linear-interpolant integral residual only; rigorous measurement-noise, intersample and spatial-closure enclosure absent')

def source_inverse(a, fit):
    p = fit['parameters']
    (dt, d1, d2, i1, i2) = intervals(a)
    energy = float(np.sum(a[:-1, 5] * dt))
    if energy <= 0:
        raise ValueError('Positive independently measured drilling work required')
    q1 = p['C1_J_K'] * d1.sum() + p['G_W_K'] * (i1 - i2).sum() + p['H1_W_K'] * i1.sum()
    q2 = p['C2_J_K'] * d2.sum() - p['G_W_K'] * (i1 - i2).sum() + p['H2_W_K'] * i2.sum()
    eta = float((q1 + q2) / energy)
    chi = float(q2 / (q1 + q2))
    if not (0 <= eta <= 1 and 0 <= chi <= 1):
        raise ValueError('Energy/source fraction inadmissible: eta=' + str(eta) + ' chi=' + str(chi))
    return dict(mechanical_work_J=energy, source_heat_J=float(q1 + q2), eta=eta, chi_region2=chi, source_region_heat_J=[float(q1), float(q2)], protocol_by_bone_class='UNKNOWN', CEM43_safety='UNKNOWN: point-to-wall and empirical bone outcome closure unvalidated')

def run(cal_csv, drill_csv, cal_meta, drill_meta):
    check_metadata(cal_meta, drill_meta)
    fit = calibrate(read_csv(cal_csv), cal_meta['capacity_bounds_J_K'])
    inv = source_inverse(read_csv(drill_csv), fit)
    return dict(calibration=fit, source=inv, external_measurement=cal_meta.get('kind') == 'independent_measurement', physical_validation='UNKNOWN until independent held-out same-material and observer/outcome measurements')
if __name__ == '__main__':
    import json
    ap = argparse.ArgumentParser()
    ap.add_argument('heater_csv')
    ap.add_argument('drill_csv')
    ap.add_argument('heater_meta')
    ap.add_argument('drill_meta')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = run(args.heater_csv, args.drill_csv, json.load(open(args.heater_meta)), json.load(open(args.drill_meta)))
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=2)
