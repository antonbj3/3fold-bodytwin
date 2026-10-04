"""One-command regional query from a supplied local-reference measurement JSON.

Input is a lab mechanical calibration contract, not patient assessment.
All numerical interval propagation is rational before outward float conversion.
"""
import json, sys
from pathlib import Path
from fractions import Fraction as Q
from round3 import response_bounds, directed
from motion_port import local_threshold

def main():
    d = json.loads(Path(sys.argv[1]).read_text())
    if d.get('units') != {'position': 'mm', 'force': 'N', 'moment': 'Nmm', 'displacement': 'um'}:
        raise ValueError('Explicit units must be mm, N, Nmm, um')
    if d.get('bone_marker_readings_um') is None:
        print(json.dumps({'status': 'UNKNOWN_MISSING_LOCAL_BONE_REFERENCE'}))
        return
    for key in ['calibration_wrenches_columns_N_Nmm', 'implant_marker_readings_um', 'bone_marker_readings_um']:
        if len(d[key]) != 2 or any((len(row) != 2 for row in d[key])):
            raise ValueError(key + ' must be a 2-by-2 matrix')
    if len(d['marker_positions_mm']) != 2 or len(d['query_wrench_N_Nmm']) != 2:
        raise ValueError('This in-plane port needs exactly two marker positions and a force/moment pair')
    if not d['regional_positions_mm']:
        raise ValueError('At least one requested region is required')
    scale = lambda m: [[Q(str(x)) / 1000 for x in row] for row in m]
    zi = [Q(str(x)) for x in d['marker_positions_mm']]
    W = [[Q(str(x)) for x in row] for row in d['calibration_wrenches_columns_N_Nmm']]
    w = [Q(str(x)) for x in d['query_wrench_N_Nmm']]
    zq = [Q(str(x)) for x in d['regional_positions_mm']]
    yi = scale(d['implant_marker_readings_um'])
    yb = scale(d['bone_marker_readings_um'])
    error = Q(str(d['absolute_error_um_each_reading'])) / 1000
    if error < 0:
        raise ValueError('error bound must be nonnegative')
    r = response_bounds(zi, W, w, zq, yi, yb, error)
    if not r['status'].startswith('BOUNDED'):
        print(json.dumps(r))
        return
    validation = {'status': 'UNKNOWN_NO_HELD_WRENCH_OBSERVATIONS'}
    held = d.get('held_wrench_validation')
    if held is not None:
        if held.get('bone_readings_um') is None:
            print(json.dumps({'status': 'UNKNOWN_MISSING_HELD_LOCAL_BONE_REFERENCE'}))
            return
        hz = [Q(str(x)) for x in held['positions_mm']]
        hw = [Q(str(x)) for x in held['wrench_N_Nmm']]
        hy = [Q(str(x)) / 1000 for x in held['implant_readings_um']]
        hb = [Q(str(x)) / 1000 for x in held['bone_readings_um']]
        he = Q(str(held['absolute_error_um_each_reading'])) / 1000
        if not hz or len(hw) != 2 or len(hy) != len(hz) or (len(hb) != len(hz)) or (he < 0):
            raise ValueError('Held validation needs positions, matching implant/bone vectors, a force/moment pair and nonnegative error')
        bounds = response_bounds(zi, W, hw, hz, yi, yb, error)['rational_intervals']
        failures = [i for (i, ((lo, hi), a, b)) in enumerate(zip(bounds, hy, hb)) if a - b + 2 * he < lo or a - b - 2 * he > hi]
        validation = {'status': 'REJECT_HELD_WRENCH_OUTSIDE_READING_BOXES' if failures else 'PASS_CONSISTENT_WITH_DECLARED_BOXES', 'rejected_regions': failures, 'scope': 'necessary held-load consistency; does not certify all future loads or nonlinear/spatial residuals'}
        if failures:
            print(json.dumps({'status': validation['status'], 'held_wrench_validation': validation, 'mechanical_query': 'UNKNOWN', 'biological_outcome': 'UNKNOWN'}, indent=2))
            return
    intervals = [[directed(a * 1000, 'lo'), directed(b * 1000, 'hi')] for (a, b) in r['rational_intervals']]
    lo = max((Q(0) if a <= 0 <= b else min(abs(a), abs(b)) for (a, b) in r['rational_intervals']))
    hi = max((max(abs(a), abs(b)) for (a, b) in r['rational_intervals']))
    peak = [directed(lo * 1000, 'lo'), directed(hi * 1000, 'hi')]
    print(json.dumps({'status': 'CONDITIONAL_ON_CALIBRATED_LINEAR_REVERSIBLE_BRANCH', 'regional_intervals_um': intervals, 'sampled_regional_peak_interval_um': peak, 'mechanical_research_band': local_threshold(peak, True), 'biological_outcome': 'UNKNOWN', 'whole_interface_peak': 'UNKNOWN_WITHOUT_SPATIAL_RESIDUAL_BOUND', 'resolution': 'PER_SURFACE_REGION', 'time_scale': 'SIMULTANEOUS', 'measurement_provenance': d.get('measurement_provenance', 'UNKNOWN'), 'held_wrench_validation': validation, 'required_assumptions': ['marker positions and force/moment known exactly in this enclosure', 'same reversible linear branch as calibration', 'rigid two-mode kinematics at sampled regions', 'declared absolute reading error is a bound, not sensor resolution']}, indent=2))
if __name__ == '__main__':
    main()
