"""Future same-specimen pressure measurement port. Missing data remain NOT_RUN."""
from common import *

def compare(path):
    measurement = json.loads(Path(path).read_text())
    case = int(measurement['case'])
    fdi = int(measurement['fdi'])
    rows = json.loads((H / 'rounds/R2.json').read_text())['stress_results']
    r = next((r for r in rows if r['case'] == case and r['fdi'] == fdi))
    z = np.load(r['basis']['path'])
    n = len(z['area_shares'])
    required = ['independent_measurement_locator', 'measurement_sha256', 'same_specimen_registration_locator', 'total_force_N', 'contact_forces_N', 'force_tolerance_N', 'resolution', 'recording_time', 'force_calibration_locator', 'force_sd_N']
    for k in required:
        if k not in measurement:
            raise ValueError('Missing measurement field: ' + k)
    if measurement['resolution'] != 'PER_SURFACE_REGION':
        raise ValueError('Contact forces must be resolved per original pressure region')
    if not measurement['independent_measurement_locator'] or not measurement['same_specimen_registration_locator']:
        raise ValueError('Independent source and same-specimen registration required')
    forces = np.asarray(measurement['contact_forces_N'], float)
    if forces.shape != (n,) or not np.isfinite(forces).all() or forces.min() < 0:
        raise ValueError('Patch force shape/nonnegativity failed')
    total = float(measurement['total_force_N'])
    tol = float(measurement['force_tolerance_N'])
    if total <= 0 or tol < 0:
        raise ValueError('Invalid force resultant/tolerance')
    if abs(forces.sum() - total) > tol:
        raise ValueError('Measured patch forces do not integrate to independent resultant')
    if len(measurement['force_sd_N']) != n:
        raise ValueError('Uncertainty must be per contact')
    from pressure import peak
    tensor = np.einsum('j,jea->ea', forces / 100, z['stress_tensors_MPa'])
    result = dict(case=case, fdi=fdi, resolution='PER_TOOTH', conditional_predicted_tensile_peak_MPa=peak(tensor), within_original_simplex_bound=peak(tensor) <= r['old_sharp_upper_MPa'] * total / 100 * (1 + 1e-09), force_integral_error_N=float(abs(forces.sum() - total)), measurement=measurement, prediction_freeze_sha256=sha(H / 'FROZEN_PREDICTIONS_R2.json'), scope='Force observation updates the fixed X18 roof only; pressure geometry/strain/support model needs independent validation. Statistical stress uncertainty not inferred from unreported force covariance.')
    return result
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(json.dumps(dict(status='NOT_RUN', reason='No independent registered same-specimen patch pressure measurement exists locally', specification=str(H / 'MEASUREMENT_SPEC.md')), indent=2))
    else:
        result = compare(sys.argv[1])
        dump(H / 'raw/INDEPENDENT_MEASUREMENT_COMPARISON.json', result)
        print(json.dumps(result, indent=2))
