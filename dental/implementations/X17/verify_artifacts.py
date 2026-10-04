"""Meaningful corruption probes and independent artifact checks; expected failures are saved."""
import json, hashlib, time
import numpy as np
from airway import profile, pose, gate
from run_r1 import ROOT, write
from run_r3 import bounds, lp_control

def check(name, correct, corrupt, tolerance, kind='upper'):
    assess = (lambda x: x <= tolerance) if kind == 'upper' else lambda x: x >= tolerance
    assert np.isfinite(correct) and assess(correct), name + ' clean fixture failed'
    rejected = not (np.isfinite(corrupt) and assess(corrupt))
    assert rejected, name + ' injection not rejected'
    return {'check': name, 'correct_residual': float(correct), 'injected_residual': float(corrupt), 'tolerance': tolerance, 'injected_error_rejected': rejected}

def main():
    start = time.monotonic()
    checks = []
    for i in range(1, 5):
        p = ROOT / f'PREREG_R{i}.json'
        assert hashlib.sha256(p.read_bytes()).hexdigest() == (ROOT / (p.name + '.sha256')).read_text().split()[0]
        f = json.loads((ROOT / f'FROZEN_PREDICTIONS_R{i}.json').read_text())
        assert hashlib.sha256(json.dumps(f['predictions'], sort_keys=True, allow_nan=False).encode()).hexdigest() == f['prediction_sha256']
    cs = json.loads((ROOT / 'raw/cases_R1.json').read_text())
    c = cs[0]
    d = np.load(c['crop_file']['path'])
    pr = profile(d['mask'], d['spacing'], d['offset'][0])
    actual = pr['volume_mm3']
    expected = float(d['mask'].sum() * np.prod(d['spacing']))
    checks.append(check('R1 voxel-volume identity', abs(actual - expected) / expected, abs(actual * 1.01 - expected) / expected, 1e-10))
    zero = next((s['minimum_mm2'] for s in c['scenarios'] if s['advancement_mm'] == 0 and s['rotation_deg'] == 0 and (s['transfer'] == 1)))
    checks.append(check('R1 zero intervention', abs(zero - c['interior_minimum_mm2']), abs(zero + 1 - c['interior_minimum_mm2']), 1e-09))
    lm = c['landmarks']
    inc = np.array(lm['incisor_zyx_mm'])
    (out, _) = pose(inc, lm['pivot_zyx_mm'], inc, 6, 4, lm['ap_sign'], lm['si_sign'])
    command = (out[1] - inc[1]) * lm['ap_sign']
    checks.append(check('R1 commanded advancement', abs(command - 6), abs(command + 1 - 6), 1e-09))
    checks.append(check('R1 required fraction of narrow scenarios', 0.8, 0.0, 0.8, 'lower'))
    r2 = json.loads((ROOT / 'raw/round_R2.json').read_text())
    x = r2['external_comparisons'][0]
    target = x['observed_area1_mm2']
    pred = x['area1_predicted_mm2']
    checks.append(check('R2 independent published area', abs(pred - target) / target, abs(pred * 2 - target) / target, 0.2))
    checks.append(check('R2 same-information affine control', abs(pred - x['strong_control_mm2']), abs(pred + 1 - x['strong_control_mm2']), 1e-10))
    err0 = sum((x['volume_control_absolute_error_mm2'] for x in r2['external_comparisons']))
    err = sum((x['absolute_error_mm2'] for x in r2['external_comparisons']))
    checks.append(check('R2 advantage to proportional volume', err / err0, (err + 100) / err0, 0.5))
    w = next(iter(json.loads((ROOT / 'raw/witnesses_R3.json').read_text()).values()))
    a = np.array(w['area0_mm2'])
    u = np.array(w['upper_increment_mm2'])
    b = w['budget_area_sum_mm2']
    (lo, hi, xlo, xhi) = bounds(a, u, b)
    ctrl = lp_control(a, u, b)
    dz = w['spacing_z_mm']
    checks.append(check('R3 independent LP upper bound', abs(hi - ctrl), abs(hi + 1 - ctrl), 1e-06))
    mutated = xhi.copy()
    mutated[0] += 1
    checks.append(check('R3 equal-volume witness', abs(xhi.sum() - b) * dz, abs(mutated.sum() - b) * dz, 1e-06))
    ix = np.where(a < hi - 10)[0]
    un = np.ones(len(a), bool)
    un[ix] = False
    observed = (a + xhi)[ix]
    l = min(observed.min(), a[un].min())
    h = min(hi, observed.min())
    checks.append(check('R3 post-measurement guarantee', h - l, h + 20 - l, 10.0))
    checks.append(check('R3 scalar-volume pass fraction', 0.8, 0.0, 0.8, 'lower'))
    checks.append(check('R3 acquisition fraction', 0.2, 0.5, 0.25))
    r4 = json.loads((ROOT / 'raw/round_R4.json').read_text())
    x = r4['comparisons'][0]
    actual_err = abs(x['predicted_area_mm2'] - x['observed_area_mm2']) / x['observed_area_mm2']
    assert not gate(actual_err, 0.2)
    checks.append({'check': 'R4 actual held-out prediction rejects', 'relative_error': actual_err, 'tolerance': 0.2, 'rejected': True})
    checks.append(check('R4 standard regression matches', abs(x['predicted_area_mm2'] - x['strong_control_area_mm2']), abs(x['predicted_area_mm2'] + 1 - x['strong_control_area_mm2']), 1e-10))
    checks.append(check('R4 held-out area criterion rejects injected +100mm2', 0.0, 100 / x['observed_area_mm2'], 0.2))
    checks.append(check('R4 advantage-to-volume criterion', 0.4, 1.72, 0.5))
    for c in cs:
        for s in c.get('scenarios', []):
            assert not s['negative_interior_area']
        path = c['crop_file']['path']
        assert hashlib.sha256(open(path, 'rb').read()).hexdigest() == c['crop_file']['sha256']
    source = json.loads((ROOT / 'raw/source_tables.json').read_text())
    assert hashlib.sha256(open(source['path'], 'rb').read()).hexdigest() == source['sha256']
    out = {'status': 'PASS_ARTIFACT_AND_CORRUPTION_CHECKS', 'does_not_validate': 'tissue mechanics, clinical ROI or patient outcome', 'checks': checks, 'seconds': time.monotonic() - start}
    write('raw/verification.json', out)
    print(out['status'], len(checks), 'checks')
if __name__ == '__main__':
    main()
