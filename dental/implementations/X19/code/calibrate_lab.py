"""Usage: python3 code/calibrate_lab.py measurements.json --out fit.json

A real six-axis, same-protocol input is required; bundled published scalar data
are not labelled as local bench measurements. No full 6x6 tooth stiffness claim.
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear

def calibrate(x):
    required = ['arch_geometry_sha256', 'tooth_fdi', 'movement_unit', 'moment_origin_xyz_mm', 'frame', 'seating_protocol', 'temperature_C', 'wet_condition', 'elapsed_seconds', 'material_product', 'forming_protocol', 'replicates']
    missing = [k for k in required if k not in x or x[k] is None]
    if missing or len(x.get('replicates', [])) < 2:
        return {'status': 'NEEDS_MATCHED_MEASUREMENT', 'missing': missing, 'reason': 'Two measured formed-thickness states and signed six-axis zero/intermediate/target-displacement measurements required'}
    axis = np.array(x['movement_unit'], float)
    axis /= np.linalg.norm(axis)
    rs = x['replicates']
    fits = []
    for r in rs:
        d = np.array([a['displacement_mm'] for a in r['samples']])
        W = np.array([a['wrench_N_Nmm'] for a in r['samples']])
        h = r['measured_effective_thickness_mm']
        E = r['measured_E_MPa']
        if len(d) < 3 or W.shape != (len(d), 6) or (not np.any(d == 0)) or (h <= 0) or (E <= 0):
            raise ValueError('Missing zero/intermediate/target signed six-axis physical data')
        A = np.column_stack([d, np.ones(len(d))])
        (slope, intercept) = np.linalg.lstsq(A, W, rcond=None)[0]
        k = float(-slope[:3] @ axis)
        if k <= 0:
            raise ValueError('Non-restoring primary measured force slope')
        fit = A @ np.stack([slope, intercept])
        res = float(np.max(np.abs(fit - W)))
        fits.append({'specimen': r['specimen'], 'measured_h_mm': h, 'E_MPa': E, 'K_N_per_mm': k, 'slope_N_Nmm_per_mm': slope.tolist(), 'intercept_N_Nmm': intercept.tolist(), 'max_6axis_linearity_residual': res})
    X = np.column_stack([[1 / (a['E_MPa'] * a['measured_h_mm'] ** 3) for a in fits], np.ones(len(fits))])
    y = 1 / np.array([a['K_N_per_mm'] for a in fits])
    ab = np.linalg.lstsq(X, y, rcond=None)[0]
    if min(ab) < 0:
        raise ValueError('Measured states refute passive shell-series closure')
    ctrl = lsq_linear(X, y, bounds=(0, np.inf), tol=1e-14).x
    return {'status': 'CALIBRATED_CONDITIONAL_LOCAL_PORT', 'claim_type': 'capability', 'a_geometry_mm2': float(ab[0]), 'b_effective_mm_per_N': float(ab[1]), 'fit': fits, 'same_information_control_relative_difference': float(np.linalg.norm(ab - ctrl) / max(np.linalg.norm(ab), 1e-12)), 'context': {k: x[k] for k in required if k != 'replicates'}, 'uncertainty': 'Replicate and independent holdout measurements required; no material-level or clinical transfer inferred', 'query_equation': 'K(E,h)=1/(a_geometry/(E*h^3)+b_effective); wrench slopes/intercepts remain specimen-conditioned until separate component/holdout checks pass'}

def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('measurements', type=Path)
    a.add_argument('--out', type=Path, required=True)
    o = a.parse_args()
    x = json.load(open(o.measurements))
    r = calibrate(x)
    o.out.write_text(json.dumps(r, indent=2, allow_nan=False) + '\n')
    print(r['status'])
if __name__ == '__main__':
    main()
