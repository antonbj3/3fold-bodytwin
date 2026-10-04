import argparse, csv, math, hashlib, json, datetime
from pathlib import Path
from statistics import mean, stdev, NormalDist
from analyze_r1 import R, read, dump, sha
PROTOCOLS = {'C0_PROTT': dict(crown='inCoris TZI', angle=0, E=8600, diameter=6.36, speed=1.5, cement='RelyX Unicem 2 Automix', die='Tetric Evo Ceram A2 Dentin', surface='MDP_no_abrasion'), 'C30_CHEN': dict(crown='Katana HT', angle=30, E=2100, diameter=3.5, speed=0.5, cement='Rely X Luting Cement', die='NextDent C&B MFH', surface='Al2O3_50um_0.2MPa_15s')}
W = math.log(1 / 0.8) / math.log(1.5 / 0.8)
H = 0.2
Z = NormalDist().inv_cdf(0.9875)
REQUIRED = ['specimen_id', 'protocol_id', 'design_id', 'force_unit', 'fracture_force_N', 'die_E_MPa', 'indenter_diameter_mm', 'crosshead_mm_min', 'load_angle_deg', 'material_batch', 'crown_product', 'cement_product', 'cement_batch', 'die_product', 'die_batch', 'surface_protocol', 'storage_days', 'interlayer_spec', 'manufacturing_spec_sha256', 'failure_mode', 'fracture_origin', 'force_trace_sha256', 'measured_utc', 'status']

def canonical(rows):
    return hashlib.sha256(json.dumps(sorted(rows, key=lambda r: r['specimen_id']), sort_keys=True, separators=(',', ':')).encode()).hexdigest()

def load(path):
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        missing = set(REQUIRED) - set(r.fieldnames or [])
        if missing:
            raise ValueError('missing columns: ' + str(sorted(missing)))
        return list(r)

def integrity():
    f = read('FROZEN_PREDICTIONS.json')
    if sha('FROZEN_PREDICTIONS.json') != (R / 'FROZEN_PREDICTIONS.sha256').read_text().strip():
        raise ValueError('frozen prediction hash drift')
    for (p, h) in f['locked_files'].items():
        if sha(p) != h:
            raise ValueError('locked artifact drift: ' + p)
    return f

def validate(rows, designs):
    if len({r['specimen_id'] for r in rows}) != len(rows):
        raise ValueError('duplicate specimen IDs')
    groups = {}
    problems = []
    for r in rows:
        p = PROTOCOLS.get(r['protocol_id'])
        d = r['design_id']
        if p is None or d not in designs:
            raise ValueError('unknown protocol or wrong calibration/validation design')
        if r['force_unit'] != 'N':
            raise ValueError('force unit must be N')
        for (field, target, tol) in [('die_E_MPa', p['E'], p['E'] * 0.1), ('indenter_diameter_mm', p['diameter'], p['diameter'] * 0.01), ('crosshead_mm_min', p['speed'], p['speed'] * 0.05), ('load_angle_deg', p['angle'], 1.0)]:
            if not math.isfinite(float(r[field])) or abs(float(r[field]) - target) > tol:
                raise ValueError('protocol mismatch ' + field)
        for (field, key) in [('crown_product', 'crown'), ('cement_product', 'cement'), ('die_product', 'die'), ('surface_protocol', 'surface')]:
            if r[field] != p[key]:
                raise ValueError('protocol mismatch ' + field)
        if not 7 <= float(r['storage_days']) <= 8:
            raise ValueError('storage period not matched')
        for field in ['material_batch', 'cement_batch', 'die_batch', 'interlayer_spec']:
            if not r[field]:
                raise ValueError('missing metadata ' + field)
        for field in ['force_trace_sha256', 'manufacturing_spec_sha256']:
            if len(r[field]) != 64 or any((c not in '0123456789abcdef' for c in r[field])):
                raise ValueError('missing/invalid hash ' + field)
        date = datetime.datetime.fromisoformat(r['measured_utc'].replace('Z', '+00:00'))
        if date.tzinfo is None:
            raise ValueError('measurement timestamp needs timezone')
        if r['status'] != 'measured':
            problems.append('censored/unmeasured ' + r['specimen_id'])
        else:
            f = float(r['fracture_force_N'])
            if not math.isfinite(f) or f <= 0:
                raise ValueError('force must be positive finite N')
        groups.setdefault((r['protocol_id'], d), []).append(r)
    for prot in PROTOCOLS:
        rr = [r for r in rows if r['protocol_id'] == prot]
        for field in ['material_batch', 'cement_batch', 'die_batch', 'manufacturing_spec_sha256', 'interlayer_spec']:
            if len({r[field] for r in rr}) > 1:
                raise ValueError('batch/spec change within protocol ' + field)
        for d in designs:
            if len(groups.get((prot, d), [])) < 12:
                problems.append('insufficient n ' + prot + ' ' + d)
    return (groups, problems)

def summary(rows):
    fs = [float(r['fracture_force_N']) for r in rows]
    m = mean(fs)
    s = stdev(fs)
    return dict(mean_N=m, sd_N=s, n=len(fs), log_mean_variance=(s / m) ** 2 / len(fs), CV=s / m)

def calibration(rows):
    (groups, problems) = validate(rows, {'D1', 'M2'})
    if problems:
        raise ValueError('UNKNOWN calibration: ' + '; '.join(problems))
    p = {}
    for prot in PROTOCOLS:
        a = summary(groups[prot, 'D1'])
        b = summary(groups[prot, 'M2'])
        force = a['mean_N'] ** (1 - W) * b['mean_N'] ** W
        p[prot] = dict(D1=a, M2=b, M1_point_N=force, M1_operational_window_N=[force * math.exp(-H), force * math.exp(H)], local_exponent=math.log(b['mean_N'] / a['mean_N']) / math.log(1.5 / 0.8))
    return dict(schema='x1c-physical-calibration-v1', claim_type='information_link', frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), anchor_rows_sha256=canonical(rows), formula_frozen_sha256=sha('FROZEN_PREDICTIONS.json'), predictions=p, measurement_status='USER_SUPPLIED_PHYSICAL_MEASUREMENTS; NOT INDEPENDENTLY_VERIFIED', resolution='PER_TOOTH')

def compare(rows, cal, check_time=True):
    (groups, problems) = validate(rows, {'D1', 'M1', 'M2'})
    anchors = [r for r in rows if r['design_id'] != 'M1']
    if canonical(anchors) != cal['anchor_rows_sha256']:
        raise ValueError('calibration rows changed after freeze')
    if cal['formula_frozen_sha256'] != sha('FROZEN_PREDICTIONS.json'):
        raise ValueError('wrong formula freeze')
    if check_time:
        stamp = datetime.datetime.fromisoformat(cal['frozen_utc'])
        if any((datetime.datetime.fromisoformat(r['measured_utc'].replace('Z', '+00:00')) <= stamp for r in rows if r['design_id'] == 'M1')):
            raise ValueError('M1 measured before numeric prediction freeze')
    if problems:
        return dict(outcome='UNKNOWN', reasons=problems, submitted_rows=len(rows), dropout_fraction=sum((r['status'] != 'measured' for r in rows)) / len(rows))
    out = []
    for prot in PROTOCOLS:
        a = summary(groups[prot, 'D1'])
        b = summary(groups[prot, 'M2'])
        v = summary(groups[prot, 'M1'])
        pred = a['mean_N'] ** (1 - W) * b['mean_N'] ** W
        e = math.log(v['mean_N'] / pred)
        se = math.sqrt(v['log_mean_variance'] + (1 - W) ** 2 * a['log_mean_variance'] + W ** 2 * b['log_mean_variance'])
        ci = [e - Z * se, e + Z * se]
        verdict = 'REJECTED' if ci[1] < -H or ci[0] > H else 'SUPPORTED_OPERATIONAL_TEST' if -H <= ci[0] and ci[1] <= H else 'UNKNOWN'
        mech = []
        for d in ['D1', 'M1', 'M2']:
            rr = groups[prot, d]
            bad = sum((r['failure_mode'] in ['die_fracture', 'cement_debonding'] for r in rr))
            unknown = sum((r['failure_mode'] not in ['crown_tensile_fracture', 'contact_damage', 'die_fracture', 'cement_debonding'] or r['fracture_origin'] == 'unknown' for r in rr))
            m = 'FAIL_MECHANISM' if bad / len(rr) > 0.5 else 'UNKNOWN' if bad or unknown else 'CROWN_FAILURE_OBSERVABLE_ONLY'
            mech.append(dict(design=d, status=m, incompatible_fraction=bad / len(rr), unknown_fraction=unknown / len(rr)))
        if any((z['CV'] > 0.5 for z in [a, b, v])):
            verdict = 'UNKNOWN_HIGH_SCATTER'
        out.append(dict(protocol=prot, group_summaries=dict(D1=a, M1=v, M2=b), predicted_M1_N=pred, residual_log=e, measurement_CI_log=ci, CI_level=0.975, operational_halfwidth_log=H, force_verdict=verdict, mechanism_gates=mech, resolution='PER_TOOTH'))
    return dict(claim_type='information_link', outcome='MEASURED_COMPARISON', protocols=out, submitted_rows=len(rows), rejected_rows=0, dropout_fraction=0, scientific_status='PENDING_INDEPENDENT_REVIEW', note='delta-method measurement CI; operational shape window is not a population prediction interval')

def main():
    pa = argparse.ArgumentParser()
    pa.add_argument('action', choices=['calibrate', 'compare'])
    pa.add_argument('csv')
    pa.add_argument('--output', required=True)
    pa.add_argument('--calibration')
    a = pa.parse_args()
    integrity()
    try:
        rows = load(a.csv)
        if a.action == 'calibrate':
            if Path(a.output).exists():
                raise ValueError('refuse overwrite of frozen calibration')
            result = calibration(rows)
        else:
            if not a.calibration:
                raise ValueError('calibration file required')
            cf = Path(a.calibration)
            ch = cf.with_suffix('.sha256')
            if not ch.exists() or ch.read_text().strip() != hashlib.sha256(cf.read_bytes()).hexdigest():
                raise ValueError('calibration hash drift')
            result = compare(rows, json.loads(cf.read_text()))
    except (ValueError, KeyError, TypeError) as e:
        result = dict(outcome='FAIL_INPUT', reason=str(e))
        print(json.dumps(result))
        raise SystemExit(2)
    f = Path(a.output)
    f.write_text(json.dumps(result, indent=2) + '\n')
    if a.action == 'calibrate':
        f.with_suffix('.sha256').write_text(hashlib.sha256(f.read_bytes()).hexdigest() + '\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
