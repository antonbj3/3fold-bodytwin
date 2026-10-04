"""Research laboratory query over measured local calibration and a bore menu.

Usage: python3 code/query.py --input raw/LAB_QUERY_INPUT.json
Two points fit a linear transverse-area proxy; three points fit a quadratic.
Every physical guarantee remains null without a physical response error enclosure.
"""
import argparse, json, math, pathlib
import numpy as np
from calibration import monotone_inverse
from enclosures import certified_derivative

def query(data):
    if data.get('same_regime_confirmed') is not True:
        return {'status': 'UNKNOWN_REGIME_COMPATIBILITY', 'guaranteed_fixed_bore_mm': None, 'radial_pressure_MPa': None}
    if data.get('measurement_level') not in ['POPULATION', 'PER_TOOTH']:
        raise ValueError('Explicit measurement_level required; no fabricated patient-level data')
    if data.get('units') != 'Ncm':
        raise ValueError('Measured torque units must be Ncm')
    points = sorted(data['calibration'], key=lambda p: p['Df_mm'])
    if len(points) not in [2, 3]:
        raise ValueError('Exactly two or three measured calibration conditions required')
    ds = np.array([p['Df_mm'] for p in points], dtype=float)
    ts = np.array([p['torque_Ncm'] for p in points], dtype=float)
    D = float(data['implant_outer_D_mm'])
    target = float(data['target_Ncm'])
    if not (np.all(np.isfinite(ds)) and np.all(np.isfinite(ts)) and np.all(ds > 0) and np.all(ts > 0) and math.isfinite(D) and (D > 0) and math.isfinite(target) and (target > 0)):
        raise ValueError('Finite positive diameters, torque and target required')
    inv = monotone_inverse(list(zip(ds, ts)), target)
    if inv['status'] in ['INCONSISTENT', 'OUTSIDE_CALIBRATION_RANGE']:
        return dict(status=inv['status'], inverse=inv, guaranteed_fixed_bore_mm=None, radial_pressure_MPa=None)
    u = D * D - ds * ds
    if len(points) == 2:
        b = (ts[0] - ts[1]) / (u[0] - u[1])
        poly = np.array([0.0, b, ts[0] - b * u[0]])
    else:
        poly = np.polyfit(u, ts, 2)
    derivative = certified_derivative(poly, D, float(ds[0]), float(ds[-1]))
    if not derivative['certified_nonincreasing_fitted_torque']:
        return dict(status='UNKNOWN_FITTED_MONOTONICITY', fitted_derivative_enclosure=derivative, inverse=inv, guaranteed_fixed_bore_mm=None, radial_pressure_MPa=None)
    menu = list(data['available_final_bore_mm'])
    valid = [float(d) for d in menu if math.isfinite(float(d)) and ds[0] <= float(d) <= ds[-1]]
    dropped = len(menu) - len(valid)
    preds = [{'Df_mm': d, 'predicted_Ncm': float(np.polyval(poly, D * D - d * d))} for d in sorted(set(valid))]
    if not preds:
        return dict(status='UNKNOWN_NO_AVAILABLE_OPTION_IN_CALIBRATION_RANGE', inverse=inv, guaranteed_fixed_bore_mm=None, radial_pressure_MPa=None, rejected_menu_options=dropped)
    if any((p['predicted_Ncm'] <= 0 for p in preds)):
        return dict(status='INCONSISTENT_NONPOSITIVE_MODEL', inverse=inv, guaranteed_fixed_bore_mm=None, radial_pressure_MPa=None)
    chosen = min(preds, key=lambda p: (abs(math.log(p['predicted_Ncm'] / target)), -p['Df_mm']))
    return dict(status='CONDITIONAL_RESEARCH_CANDIDATE', measurement_level=data['measurement_level'], model_level='PHENOMENOLOGICAL', source_regime=data.get('regime'), target_Ncm=target, selected=chosen, menu_predictions=preds, model_coefficients_c_b_a=poly.tolist(), fitted_derivative_enclosure=derivative, inverse=inv, rejected_menu_options=dropped, nominal_menu_not_actual_profile='Actual bore/contact metrology still required', guaranteed_fixed_bore_mm=None, physical_response_error_enclosure='MISSING', radial_pressure_MPa=None, note='Calibration means/individual specimens do not guarantee a new specimen. This query does not predict patient outcomes.')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=pathlib.Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    items = data if isinstance(data, list) else [data]
    print(json.dumps([query(v) for v in items], indent=2, allow_nan=False))
if __name__ == '__main__':
    main()
