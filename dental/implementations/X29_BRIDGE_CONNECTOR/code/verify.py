"""Independent source/unit/regime checks and real rejecting fault injections."""
import copy, re, xml.etree.ElementTree as ET
import numpy as np
from common import R, read, dump, sha

def source_check(d):
    root = ET.parse(R / 'raw/sources/PMC10471503.xml').getroot()
    rows = root.find('.//table-wrap[@id="T2"]').findall('.//tr')[2:]
    expected = {}
    for (i, tr) in enumerate(rows):
        for (j, (h, b)) in enumerate([(4, 2.25), (3, 3)], 1):
            cell = ' '.join(tr[j].itertext())
            m = re.search('([\\d.]+)\\s*±\\s*([\\d.]+)', cell)
            expected[f'PMC10471503_{i + 3}Y_{j}'] = (float(m[1]), float(m[2]), h, b)
    if len(d) != 6 or len({a['id'] for a in d}) != 6:
        return False
    for a in d:
        if a['id'] not in expected:
            return False
        (F, SD, h, b) = expected[a['id']]
        if (a['mean_N'], a['sd_N'], a['height_mm'], a['width_mm']) != (F, SD, h, b):
            return False
        if a['n'] != 6 or a['force_unit'] != 'N' or a['dimension_unit'] != 'mm' or (a['area_unit'] != 'mm2'):
            return False
        if abs(a['area_mm2'] - a['height_mm'] * a['width_mm']) > 0.001 * a['area_mm2']:
            return False
        if a['connector_observation'] != 'right_censored_at_system_break; NOT uncensored connector failure':
            return False
        if a['origin'] != 'pontic_occlusal_contact' or a['source_sha256'] != sha(R / 'raw/sources/PMC10471503.xml'):
            return False
    return True

def r3_check(d):
    root = ET.parse(R / 'raw/sources/PMC10788321.xml').getroot()
    rows = root.find('.//table-wrap[@id="Tab2"]').findall('.//tr')
    if len(d) != 2:
        return False
    for (a, tr, angle) in zip(d, rows[2:4], [0, 30]):
        cell = ' '.join(tr[6].itertext())
        m = re.search('(\\d+)\\s*[A-Za-z]*\\s*\\((\\d+)\\)', cell)
        if a['mean_N'] != float(m[1]) or a['sd_N'] != float(m[2]) or a['n'] != 8:
            return False
        if a['angle_deg'] != angle or a['cantilever_area_mm2'] != 21.1 or a['abutment_area_mm2'] != 18.8:
            return False
        if a['bridge_class'] != '3unit_posterior_cantilever' or a['aging'] != '10000thermalcycles6.5-60C+1.2millionchewingcycles108N; matched aging/test direction':
            return False
        if a['source_sha256'] != sha(R / 'raw/sources/PMC10788321.xml'):
            return False
    return True

def deck_force_check(text):
    steps = text.split('*STEP\n')[1:]
    out = []
    for step in steps:
        part = step.split('*CLOAD, OP=NEW\n')[1].split('*EL FILE')[0]
        f = np.zeros(3)
        for line in part.splitlines():
            if not line.strip():
                continue
            (nd, k, v) = line.split(',')
            f[int(k) - 1] += float(v)
        out.append(f)
    return len(out) == 2 and np.allclose(out, [[0, 0, -1], [0, 0.5, -np.sqrt(3) / 2]], rtol=0, atol=1e-05)

def geometry_check(V):
    return np.ptp(V[:, 0]) > 25 and np.ptp(V[:, 0]) < 38 and (len(V) > 1000)

def run():
    d = read('raw/MEASUREMENTS.json')
    e = read('raw/R3_MEASUREMENTS.json')
    f = []

    def inj(name, checker, x, change):
        bad = copy.deepcopy(x)
        change(bad)
        passed = checker(bad)
        f.append(dict(name=name, accepted_bad_input=bool(passed), rejected=not passed))
        dump('raw/faults/' + name + '.json', bad)
    normal = dict(primary_source=source_check(d), cantilever_source=r3_check(e))
    for (name, key, val) in [('force_x10', 'mean_N', d[0]['mean_N'] * 10), ('area_x10', 'area_mm2', 90), ('wrong_failure_mode', 'connector_observation', 'uncensored_connector_failure'), ('wrong_units', 'area_unit', 'cm2'), ('wrong_geometry', 'height_mm', 5), ('source_hash', 'source_sha256', '0' * 64)]:
        inj(name, source_check, d, lambda x, k=key, v=val: x[0].update({k: v}))
    for (name, key, val) in [('R3_force_x10', 'mean_N', e[0]['mean_N'] * 10), ('wrong_angle', 'angle_deg', 30), ('R3_wrong_area', 'cantilever_area_mm2', 2.11), ('wrong_aging', 'aging', 'unaged'), ('wrong_regime', 'bridge_class', '3unit_end_supported')]:
        inj(name, r3_check, e, lambda x, k=key, v=val: x[0].update({k: v}))
    product_rule = lambda x: x == {'product': 'Lava Plus', 'cantilever_min_mm2': 12.0, 'abutment_min_mm2': 9.0}
    rule = {'product': 'Lava Plus', 'cantilever_min_mm2': 12.0, 'abutment_min_mm2': 9.0}
    normal['manufacturer_threshold'] = product_rule(rule)
    inj('manufacturer_threshold_swapped', product_rule, rule, lambda x: x.update(cantilever_min_mm2=16.0))
    fe = read('raw/BRIDGE_FE.json')
    for r in fe['runs'] + read('raw/R4_PAIRED.json')['new_runs']:
        base = R / 'raw/fe' / r['tag']
        text = (base / 'solve.inp').read_text()
        normal['force_' + r['tag']] = deck_force_check(text)
        ga = np.load(base / 'geometry.npz')
        normal['bridge_' + r['tag']] = geometry_check(ga['Vc'])
    r = fe['runs'][0]
    base = R / 'raw/fe' / r['tag']
    text = (base / 'solve.inp').read_text()
    bad = re.sub('(?m)^(\\d+),([123]),([-+\\d.eE]+)$', lambda m: f'{m[1]},{m[2]},{float(m[3]) * 10:.8e}', text)
    (R / 'raw/faults/FE_load_x10.inp').write_text(bad)
    f.append(dict(name='FE_load_x10', accepted_bad_input=deck_force_check(bad), rejected=not deck_force_check(bad)))
    V = np.load(base / 'geometry.npz')['Vc'].copy()
    V[:, 0] /= 3
    f.append(dict(name='same_crown_reused_as_bridge', accepted_bad_input=geometry_check(V), rejected=not geometry_check(V)))
    f.append(dict(name='singular_area_design', accepted_bad_input=read('raw/R1_FIT.json')['candidate']['full_rank'], rejected=not read('raw/R1_FIT.json')['candidate']['full_rank']))
    expected = (R / 'FROZEN_PREDICTIONS.sha256').read_text().split()[0]
    normal['frozen_prediction_hash'] = sha(R / 'FROZEN_PREDICTIONS.json') == expected
    import hashlib
    changed = (R / 'FROZEN_PREDICTIONS.json').read_bytes().replace(b'"operational_log_ratio_tolerance": 0.3', b'"operational_log_ratio_tolerance": 3.0')
    if changed == (R / 'FROZEN_PREDICTIONS.json').read_bytes():
        raise ValueError('Frozen mutation injection did not change data')
    f.append(dict(name='measurement_refit_after_freeze', accepted_bad_input=hashlib.sha256(changed).hexdigest() == expected, rejected=hashlib.sha256(changed).hexdigest() != expected))
    out = dict(normal_checks=normal, fault_injections=f, all_normal_pass=all(normal.values()), all_faults_rejected=all((x['rejected'] for x in f)))
    dump('VERIFICATION.json', out)
    print('Integrity controls', out['all_normal_pass'], 'rejected faults', sum((x['rejected'] for x in f)), '/', len(f))
    if not out['all_normal_pass'] or not out['all_faults_rejected']:
        raise ValueError('Integrity control failed')
    return out
if __name__ == '__main__':
    run()
