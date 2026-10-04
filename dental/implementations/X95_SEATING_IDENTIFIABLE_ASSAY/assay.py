"""Research bench port. No synthetic observation can receive physical PASS.

freeze --acquisition acquisition.json --training lab.csv --output frozen.json
score --predictions frozen.json --validation held.csv --output score.json
Incomplete physical logs return UNKNOWN with named missing observations.
"""
import argparse, csv, json, math
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from mpmath import iv
from common import dump, sha, utc
iv.dps = 50
KEYS = ['specimen_id', 'case_id', 'region_id', 'frame_id', 'pose_id', 'material_batch']
COLUMNS = KEYS + ['source_kind', 'time_s', 'age_s', 'measured_at_utc', 'force_N', 'h0_m', 'h1_m', 'h2_m', 'outflow_m3']
UNITS = dict(time='s', height='m', force='N', outflow='m^3', fluidity='1/(Pa*s)')

class Unknown(ValueError):
    pass

def timestamp(s):
    d = datetime.fromisoformat(s.replace('Z', '+00:00'))
    if d.tzinfo is None:
        raise Unknown('UTC time lacks timezone')
    return d.astimezone(timezone.utc)

def hull(x):
    return [float(np.nextafter(float(x.a), -np.inf)), float(np.nextafter(float(x.b), np.inf))]

def interval(lo, hi=None):
    return iv.mpf([lo, lo if hi is None else hi])

def centered(value, bound):
    return interval(value) + interval(-bound, bound)

def read_rows(path):
    with Path(path).open(newline='') as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != COLUMNS:
            raise Unknown('Missing or incorrect SI CSV columns')
        rows = list(reader)
    if not rows:
        raise Unknown('Missing synchronized hydraulic rows')
    return rows

def load_acquisition(path):
    a = json.loads(Path(path).read_text())
    if a.get('units') != UNITS:
        raise Unknown('Incorrect SI unit contract')
    for k in KEYS:
        if not a.get(k):
            raise Unknown('Missing ' + k)
    for k in ['dry_geometry_sha256', 'coordinate_frame_sha256', 'mixing_anchor_sha256', 'calibration_sha256']:
        if not a.get(k) or len(a[k]) != 64:
            raise Unknown('Missing calibrated source hash: ' + k)
    if a.get('mixing_start_utc') != a.get('independent_mixing_start_utc'):
        raise Unknown('Mixing clock disagrees with independent anchor')
    if a.get('clock_model') != 'exact_recorded_time_conditional' or a.get('unmodelled_time_uncertainty_s') != 0:
        raise Unknown('Timestamp uncertainty lacks a rigorous enclosure')
    if a.get('geometry') != 'parallel_disk' or a.get('parallelism_bound_rad') != 0:
        raise Unknown('Tilt/crown resistance has no enclosure in R2 disk port')
    if a.get('constitutive_contract') != 'uniform_newtonian_exponential_age' or not a.get('full_fill_independently_verified'):
        raise Unknown('Missing full-fill or rheology branch observations')
    if a.get('force_attribution') != 'hydraulic_only_no_contact_no_gel':
        raise Unknown('Force attribution unknown: dry contact/gel must be separated')
    if a.get('source_kind') not in ['our_own_fixture', 'independent_measurement']:
        raise Unknown('Unknown observation provenance')
    if a['source_kind'] == 'independent_measurement':
        for stem in ['dry_geometry', 'coordinate_frame', 'mixing_anchor', 'calibration']:
            path = a.get(stem + '_path')
            if not path or not Path(path).is_file() or sha(path) != a[stem + '_sha256']:
                raise Unknown('Missing physical source object or hash mismatch: ' + stem)
        cal = json.loads(Path(a['calibration_path']).read_text())
        for k in ['height_bound_m', 'force_bound_N', 'outflow_bound_m3', 'radius_bound_m', 'unmodelled_time_uncertainty_s']:
            if cal.get(k) != a.get(k):
                raise Unknown('Instrument calibration does not bind ' + k)
        if not a.get('data_object_license'):
            raise Unknown('Missing rights for the actual measurement data object')
    for k in ['height_bound_m', 'force_bound_N', 'outflow_bound_m3', 'radius_bound_m']:
        v = a.get(k)
        if v is None or not math.isfinite(v) or v < 0:
            raise Unknown('Missing calibrated bound ' + k)
    if not a.get('radius_m') or a['radius_m'] <= a['radius_bound_m']:
        raise Unknown('Invalid radius')
    xy = np.asarray(a.get('height_probe_xy_m', []), float)
    if xy.shape != (3, 2) or np.linalg.matrix_rank(np.column_stack([np.ones(3), xy])) != 3:
        raise Unknown('Three noncollinear regional height probes required')
    if not a.get('held_pulse_segments'):
        raise Unknown('Missing held force-pulse plan')
    return a

def validate_rows(rows, a, held=False, frozen=None):
    out = []
    last = -math.inf
    mix = timestamp(a['mixing_start_utc'])
    for r in rows:
        if any((r[k] != a[k] for k in KEYS)):
            raise Unknown('Specimen/case/region/frame/pose/batch mismatch')
        if r['source_kind'] != a['source_kind']:
            raise Unknown('Source-kind mismatch')
        try:
            v = {k: float(r[k]) for k in ['time_s', 'age_s', 'force_N', 'h0_m', 'h1_m', 'h2_m', 'outflow_m3']}
        except ValueError:
            raise Unknown('Missing numerical hydraulic observation')
        if not all((math.isfinite(x) for x in v.values())):
            raise Unknown('Nonfinite observation')
        if v['time_s'] <= last or v['age_s'] < 0:
            raise Unknown('Nonmonotone synchronized clock')
        last = v['time_s']
        expected = v['time_s'] + a['age_at_time_zero_s']
        if abs(v['age_s'] - expected) > 1e-09:
            raise Unknown('Mixing age mismatch')
        measured = timestamp(r['measured_at_utc'])
        if abs((measured - mix).total_seconds() - v['age_s']) > 2e-06:
            raise Unknown('Observation UTC and mixing clock disagree')
        if held and measured <= timestamp(frozen):
            raise Unknown('Held measurement occurred before prediction freeze')
        if a['source_kind'] == 'independent_measurement' and measured > timestamp(utc()):
            raise Unknown('Physical measurement timestamp is in the future')
        if v['force_N'] < 0 or min((v[k] for k in ['h0_m', 'h1_m', 'h2_m'])) <= a['height_bound_m'] or v['outflow_m3'] < 0:
            raise Unknown('Invalid positive-height/force/inventory observation')
        hh = np.array([v['h0_m'], v['h1_m'], v['h2_m']])
        if float(hh.max() - hh.min()) > 2 * a['height_bound_m']:
            raise Unknown('Local heights refute locked parallel pose')
        v['h_m'] = float(hh.mean())
        v['measured_at_utc'] = r['measured_at_utc']
        out.append(v)
    return out

def inventory(rows, a, initial_h):
    radius = centered(a['radius_m'], a['radius_bound_m'])
    area = iv.pi * radius ** 2
    h0 = centered(initial_h, a['height_bound_m'])
    residual = []
    for r in rows:
        expected = hull(area * (h0 - centered(r['h_m'], a['height_bound_m'])))
        observed = hull(centered(r['outflow_m3'], a['outflow_bound_m3']))
        ok = not (observed[1] < expected[0] or observed[0] > expected[1])
        residual.append(dict(time_s=r['time_s'], expected_m3=expected, observed_m3=observed, consistent=ok))
    return residual

def ik(k, dt):
    if k == 0:
        return interval(dt)
    return (1 - iv.exp(-interval(k) * dt)) / interval(k)

def fit_intervals(rows, a):
    if len(rows) != 3:
        raise Unknown('R2 requires three training endpoints; final gap alone cannot identify age law')
    ts = np.array([r['time_s'] for r in rows])
    dt = ts[1] - ts[0]
    if dt <= 0 or abs(ts[2] - ts[1] - dt) > 1e-12:
        raise Unknown('Training intervals must have equal known duration')
    (F1, F2) = (rows[0]['force_N'], rows[1]['force_N'])
    if min(F1, F2) <= a['force_bound_N']:
        raise Unknown('Positive training force required')
    if a.get('training_force_segments') != [[float(ts[0]), float(ts[1]), F1], [float(ts[1]), float(ts[2]), F2]]:
        raise Unknown('Independent constant-force history missing')
    R = centered(a['radius_m'], a['radius_bound_m'])
    C = 4 / (3 * iv.pi * R ** 4)
    ys = [centered(r['h_m'], a['height_bound_m']) ** (-2) for r in rows]
    es = [(ys[i + 1] - ys[i]) / C for i in range(2)]
    if min((hull(e)[0] for e in es)) <= 0:
        raise Unknown('Height resolution cannot establish positive exposure')
    ff = [centered(f, a['force_bound_N']) for f in [F1, F2]]
    kval = hull(-iv.ln(es[1] / ff[1] / (es[0] / ff[0])) / dt)
    if kval[1] < 0:
        raise Unknown('Training refutes k>=0 age contract')
    kval[0] = max(0.0, kval[0])
    ilower = hull(ik(kval[1], dt))[0]
    iupper = hull(ik(kval[0], dt))[1]
    age0 = rows[0]['age_s']
    integral = interval(ilower, iupper)
    aval = hull(es[0] * iv.exp(interval(*kval) * age0) / (ff[0] * integral))
    cn = 4 / (3 * np.pi * a['radius_m'] ** 4)
    en = [(rows[i + 1]['h_m'] ** (-2) - rows[i]['h_m'] ** (-2)) / cn for i in range(2)]
    k = -math.log(en[1] / F2 / (en[0] / F1)) / dt
    I = dt if abs(k) < 1e-12 else -math.expm1(-k * dt) / k
    a0 = en[0] * math.exp(k * age0) / (F1 * I)
    return dict(a0_nominal=a0, k_nominal=k, a0_interval=aval, k_interval=kval, C_interval=hull(C), exposure_intervals=[hull(e) for e in es], training_end_time_s=float(ts[2]), age0=age0, dt=dt, enclosure='Rigorous interval enclosure conditional on disk, exact timing and exponential law; fixture bounds are not instrument calibrations')

def integrate_nominal(segments, start, end, k, age_shift):
    val = 0.0
    for (l, u, F) in segments:
        lo = max(start, l)
        hi = min(end, u)
        if hi <= lo:
            continue
        val += F * (hi - lo) if abs(k) < 1e-12 else F * math.exp(-k * (lo + age_shift)) * -math.expm1(-k * (hi - lo)) / k
    return val

def segment_coverage(segments, start, end):
    cursor = start
    for (l, u, F) in segments:
        if l != cursor or u <= l or F < 0:
            raise Unknown('Incomplete or invalid held pulse segments')
        cursor = u
    if cursor != end:
        raise Unknown('Held pulse plan missing interval')

def predict(t, model, rows, a):
    T = model['training_end_time_s']
    hT = rows[-1]['h_m']
    shift = a['age_at_time_zero_s']
    pulse = a['held_pulse_segments']
    C = 4 / (3 * np.pi * a['radius_m'] ** 4)
    dose = integrate_nominal(pulse, T, t, model['k_nominal'], shift)
    nominal = (hT ** (-2) + C * model['a0_nominal'] * dose) ** (-0.5)
    e = interval(0.0)
    for (l, u, F) in pulse:
        lo = max(T, l)
        hi = min(t, u)
        if hi <= lo:
            continue
        kval = interval(*model['k_interval'])

        def integral_at(k):
            return iv.exp(-interval(k) * (lo + shift)) * ik(k, hi - lo)
        i_lo = hull(integral_at(model['k_interval'][1]))[0]
        i_hi = hull(integral_at(model['k_interval'][0]))[1]
        (fl, fu) = hull(centered(F, a['force_bound_N']))
        fl = max(0.0, fl)
        e += interval(fl, fu) * interval(i_lo, i_hi)
    y = centered(hT, a['height_bound_m']) ** (-2) + interval(*model['C_interval']) * interval(*model['a0_interval']) * e
    bounds = hull(1 / iv.sqrt(y))
    return dict(time_s=t, height_nominal_m=nominal, height_interval_m=bounds, resolution='PER_SURFACE_REGION')

def freeze(acquisition, training, output):
    out = Path(output)
    if out.exists():
        raise Unknown('Preserve existing frozen predictions')
    a = load_acquisition(acquisition)
    raw = read_rows(training)
    rows = validate_rows(raw, a)
    inv = inventory(rows, a, rows[0]['h_m'])
    if not all((x['consistent'] for x in inv)):
        raise Unknown('Outflow refutes full-filled film inventory')
    model = fit_intervals(rows, a)
    times = a.get('held_observation_times_s', [])
    if not times or min(times) <= model['training_end_time_s']:
        raise Unknown('Missing later held observations')
    segment_coverage(a['held_pulse_segments'], model['training_end_time_s'], max(times))
    predictions = [predict(t, model, rows, a) for t in times]
    f = dict(claim_type='capability', created_utc=utc(), acquisition_path=str(Path(acquisition).resolve()), implementation_sha256=sha(Path(__file__)), acquisition_sha256=sha(acquisition), training_path=str(Path(training).resolve()), training_sha256=sha(training), acquisition=a, model=model, predictions=predictions, training_inventory=inv, initial_h_m=rows[0]['h_m'], last_h_m=rows[-1]['h_m'], physical_status='UNKNOWN_SYNTHETIC_LOGS' if a['source_kind'] == 'our_own_fixture' else 'CONDITIONAL_AWAITING_HELD_DATA')
    dump(out, f)
    out.with_suffix(out.suffix + '.sha256').write_text(sha(out) + '\n')
    return f

def score(predictions, validation, output):
    path = Path(predictions)
    if sha(path) != path.with_suffix(path.suffix + '.sha256').read_text().strip():
        raise Unknown('Prediction hash drift')
    f = json.loads(path.read_text())
    a = f['acquisition']
    load_acquisition(f['acquisition_path'])
    if f.get('implementation_sha256') != sha(Path(__file__)):
        raise Unknown('Frozen implementation missing or changed')
    if sha(f['acquisition_path']) != f['acquisition_sha256'] or sha(f['training_path']) != f['training_sha256']:
        raise Unknown('Frozen acquisition or training changed')
    rows = validate_rows(read_rows(validation), a, held=True, frozen=f['created_utc'])
    if [r['time_s'] for r in rows] != [p['time_s'] for p in f['predictions']]:
        raise Unknown('Incomplete held pulse observations')
    inv = inventory(rows, a, f['initial_h_m'])
    scores = []
    for (r, p) in zip(rows, f['predictions']):
        force = next((F for (l, u, F) in a['held_pulse_segments'] if l < r['time_s'] <= u))
        if abs(r['force_N'] - force) > a['force_bound_N']:
            raise Unknown('Held force differs from frozen pulse')
        (lo, hi) = p['height_interval_m']
        obs = hull(centered(r['h_m'], a['height_bound_m']))
        intersects = obs[1] >= lo and obs[0] <= hi
        scores.append(dict(**p, observed_height_m=r['h_m'], observed_height_interval_m=obs, interval_consistent=intersects, absolute_nominal_error_m=abs(r['h_m'] - p['height_nominal_m'])))
    numerical = all((s['interval_consistent'] for s in scores)) and all((s['consistent'] for s in inv))
    tolerance = a.get('physical_acceptance_tolerance_m')
    physical = 'UNKNOWN_SYNTHETIC_LOGS' if a['source_kind'] == 'our_own_fixture' else 'UNKNOWN_MISSING_CALIBRATED_PHYSICAL_TOLERANCE'
    if a['source_kind'] == 'independent_measurement' and tolerance is not None:
        physical = 'PENDING_INDEPENDENT_REVIEW' if numerical else 'PHYSICAL_CONTRACT_REFUTED'
    result = dict(numerical_contract='CONSISTENT' if numerical else 'REFUTED', physical_status=physical, predictions_sha256=sha(path), validation_sha256=sha(validation), scores=scores, inventory=inv, decision_scope='Interval consistency can refute the model but does not prove the law or empirical coverage', review_state='PENDING_INDEPENDENT_REVIEW')
    if Path(output).exists():
        raise Unknown('Preserve previous score')
    dump(output, result)
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    s = parser.add_subparsers(dest='cmd', required=True)
    a = s.add_parser('freeze')
    a.add_argument('--acquisition', required=True)
    a.add_argument('--training', required=True)
    a.add_argument('--output', required=True)
    a = s.add_parser('score')
    a.add_argument('--predictions', required=True)
    a.add_argument('--validation', required=True)
    a.add_argument('--output', required=True)
    args = vars(parser.parse_args())
    cmd = args.pop('cmd')
    try:
        result = freeze(**args) if cmd == 'freeze' else score(**args)
    except (Unknown, KeyError, TypeError, FileNotFoundError) as err:
        result = dict(status='UNKNOWN', missing_observation_or_refutation=str(err), physical_status='UNKNOWN')
        if not Path(args['output']).exists():
            dump(args['output'], result)
    print(json.dumps(dict(status=result.get('numerical_contract', result.get('status', 'FROZEN')), physical_status=result.get('physical_status'), output=args['output'])))
if __name__ == '__main__':
    main()
