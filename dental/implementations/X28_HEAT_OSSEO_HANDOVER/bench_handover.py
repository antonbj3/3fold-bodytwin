"""Unit-explicit lab input. Exact integral of a measured linear interpolant.

This computes thermal screening dose, not necrosis, healing delay or ISQ.
"""
import argparse
import csv
import math
from pathlib import Path
import numpy as np
from common import ROOT, dump, cem_rate
THERMAL_COLUMNS = ['specimen_id', 'region_id', 'point_id', 'tissue', 'population', 'r_from_wall_mm', 'axial_depth_mm', 'time_s', 'temp_C', 'calibration_bound_C', 'retained_at_handover', 'source_locator']

def dose_segment(a, b, dt_s):
    if dt_s < 0:
        raise ValueError('Time must increase')
    if a < 43 < b or b < 43 < a:
        f = (43 - a) / (b - a)
        return dose_segment(a, 43, dt_s * f) + dose_segment(43, b, dt_s * (1 - f))
    r = 0.25 if (a + b) / 2 < 43 else 0.5
    c = math.log(1 / r) * (b - a)
    ratio = math.expm1(c) / c if abs(c) > 1e-12 else 1.0 + c / 2.0 + c * c / 6.0
    return dt_s / 60.0 * r ** (43 - a) * ratio

def integrate(t, temp):
    if len(t) < 2 or len(t) != len(temp) or (not np.isfinite(t).all()) or (not np.isfinite(temp).all()):
        raise ValueError('A peak is not a full finite temperature history')
    if any((b <= a for (a, b) in zip(t, t[1:]))):
        raise ValueError('Strictly increasing seconds required')
    return float(sum((dose_segment(a, b, dt) for (a, b, dt) in zip(temp, temp[1:], np.diff(t)))))

def time_above(t, temp, threshold=47.0):
    out = 0.0
    for (a, b, dt) in zip(temp, temp[1:], np.diff(t)):
        if a >= threshold and b >= threshold:
            out += dt
        elif max(a, b) > threshold:
            out += dt * (max(a, b) - threshold) / abs(b - a)
    return float(out)

def longest_time_above(t, temp, threshold=47.0):
    windows = []
    for (t0, t1, a, b) in zip(t, t[1:], temp, temp[1:]):
        if a >= threshold and b >= threshold:
            (lo, hi) = (t0, t1)
        elif a < threshold < b:
            (lo, hi) = (t0 + (t1 - t0) * (threshold - a) / (b - a), t1)
        elif b < threshold < a:
            (lo, hi) = (t0, t0 + (t1 - t0) * (a - threshold) / (a - b))
        else:
            continue
        if windows and abs(windows[-1][1] - lo) <= 1e-12:
            windows[-1][1] = hi
        else:
            windows.append([lo, hi])
    return float(max((hi - lo for (lo, hi) in windows), default=0.0))

def import_thermal(path):
    grouped = {}
    with Path(path).open() as f:
        reader = csv.DictReader(f)
        if set(THERMAL_COLUMNS) - set(reader.fieldnames or []):
            raise ValueError('Missing units/calibration/location columns')
        for row in reader:
            if any((not row[k].strip() for k in THERMAL_COLUMNS)):
                raise ValueError('Missing thermal metadata')
            if row['retained_at_handover'] not in ['true', 'false']:
                raise ValueError('Retained mask must be explicit')
            key = tuple((row[k] for k in ['specimen_id', 'region_id', 'point_id']))
            grouped.setdefault(key, []).append(row)
    if not grouped:
        raise ValueError('No lab measurement rows supplied')
    ports = []
    for ((spec, region, point), rows) in grouped.items():
        static = ['tissue', 'population', 'r_from_wall_mm', 'axial_depth_mm', 'retained_at_handover', 'source_locator']
        if any((len({r[k] for r in rows}) != 1 for k in static)):
            raise ValueError('Point metadata changed across time')
        rows.sort(key=lambda r: float(r['time_s']))
        t = np.array([float(r['time_s']) for r in rows])
        y = np.array([float(r['temp_C']) for r in rows])
        e = np.array([float(r['calibration_bound_C']) for r in rows])
        if not np.isfinite(e).all() or (e < 0).any():
            raise ValueError('Invalid calibration bound')
        if float(rows[0]['r_from_wall_mm']) < 0:
            raise ValueError('A retained-wall probe needs nonnegative distance from wall')
        if rows[0]['retained_at_handover'] == 'false':
            raise ValueError('Removed tissue cannot be healing initial state')
        doses = [integrate(t, y + offset * e) for offset in [-1, 0, 1]]
        if not doses[0] <= doses[1] <= doses[2]:
            raise ValueError('Dose interval not ordered')
        duration = time_above(t, y)
        longest = longest_time_above(t, y)
        ports.append(dict(specimen_id=spec, region_id=region, point_id=point, tissue=rows[0]['tissue'], population=rows[0]['population'], source_locator=rows[0]['source_locator'], point_r_from_wall_mm=float(rows[0]['r_from_wall_mm']), point_axial_depth_mm=float(rows[0]['axial_depth_mm']), resolution='PER_POINT', timescale='HANDOVER', retained_at_handover=True, observed_window_s=[float(t[0]), float(t[-1])], full_surgery_coverage='UNKNOWN', CEM43_min=dict(lower=doses[0], nominal=doses[1], upper=doses[2]), interpolant='piecewise-linear T between measured samples', sub_sample_peak_uncertainty='UNKNOWN: calibration interval does not bound unobserved peaks', within_window_seconds_at_or_above47=duration, within_window_longest_continuous_seconds_at_or_above47=longest, within_window_binary_flag=longest >= 60, threshold_flag_note='continuous47C exposure in observed window; cumulative/intermittent dose kept separately', full_surgery_binary_flag=True if longest >= 60 else 'UNKNOWN', thermal_to_viability_closure='PHENOMENOLOGICAL: uncalibrated', necrosis_width_mm='UNKNOWN', healing_delay_weeks='UNKNOWN', ISQ_prediction='UNKNOWN'))
    return ports

def join_injury(ports, record):
    required = ['specimen_id', 'region_id', 'tissue', 'population', 'observation_day', 'injury_width_mm', 'source_locator']
    if any((k not in record or record[k] in [None, ''] for k in required)):
        raise ValueError('Incomplete injury record')
    if float(record['observation_day']) < 0 or float(record['injury_width_mm']) < 0:
        raise ValueError('Invalid injury/time units')
    matches = [p for p in ports if all((p[k] == record[k] for k in ['specimen_id', 'region_id', 'tissue', 'population']))]
    if not matches:
        raise ValueError('No matched specimen + region + tissue + population')
    if record.get('clinical_safe') is True or record.get('ISQ_prediction') is not None:
        raise ValueError('Observed injury does not supply clinical safety or an ISQ prediction')
    return dict(edge_kind='HANDOVER', resolution='PER_SURFACE_REGION', specimen_id=record['specimen_id'], region_id=record['region_id'], thermal_points=matches, injury_observation=record, point_to_region_operator='retain all point histories; no average or spatial interpolation', causal_model_status='UNKNOWN', healing_prediction='UNKNOWN')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--thermal', required=True)
    p.add_argument('--out', default='BENCH_HANDOVER.json')
    a = p.parse_args()
    ports = import_thermal(a.thermal)
    Path(a.out).write_text(__import__('json').dumps(ports, indent=2) + '\n')
    print(f'{len(ports)} observed-point ports; no inferred healing model')
if __name__ == '__main__':
    main()
