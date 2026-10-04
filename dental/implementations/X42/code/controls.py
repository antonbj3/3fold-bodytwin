"""Reject real-development geometry faults. No original benchmark file mutated."""
from common import *
import copy, tempfile
from integrity import Bundle, IntegrityError

def run():
    root = ROOT / 'raw/control_sandbox'
    root.mkdir(exist_ok=True)
    cases = split()['dev_round1']
    chosen = None
    bridge = None
    for key in cases:
        for orig in tasks(key):
            if orig['status'] != 'READY':
                continue
            t = scene(orig)
            d = standard_lp(t)
            if all_checks(t, d)['validity'] == 'PASS':
                if chosen is None and t['A'].shape[0]:
                    chosen = (t, d)
                if t['family'] == 'bridge3' and bridge is None:
                    bridge = (t, d)
        if chosen and bridge:
            break
    if not chosen:
        raise RuntimeError('no real dev positive control')
    (t, d) = chosen
    out = []

    def probe(name, mutate, expected):
        dd = copy.deepcopy(d)
        mutate(dd)
        r = all_checks(t, dd)
        ok = r['validity'] != 'PASS' and (r['validity'] == 'INVALID' if expected == 'schema' else r['checks'].get(expected, {}).get('status') != 'PASS')
        out.append({'name': name, 'positive_source_task': t['task_id'], 'rejected': bool(ok), 'verdict': r['validity'], 'expected_check': expected, 'checks': r['checks']})
    out.append({'name': 'positive_restored', 'passed': all_checks(t, d)['validity'] == 'PASS'})
    probe('wall_low', lambda a: a.update(outer_vertices=np.c_[t['xy'], np.asarray(a['inner_vertices'])[:, 2] + t['requirements']['wall_mm'] - 0.1].tolist()), 'wall')
    probe('film_too_small', lambda a: a.update(inner_vertices=np.c_[t['xy'], t['preparation_z'] + t['requirements']['film_min_mm'] - 0.02].tolist()), 'film_min')
    probe('film_too_large', lambda a: a.update(inner_vertices=np.c_[t['xy'], t['preparation_z'] + t['requirements']['film_max_mm'] + 0.02].tolist()), 'film_max')
    probe('antagonist_penetration', lambda a: a.update(outer_vertices=np.c_[t['xy'], np.full(len(t['xy']), 50.0)].tolist()), 'antagonist')
    probe('nesting_inverted', lambda a: a.update(outer_vertices=np.c_[t['xy'], np.asarray(a['inner_vertices'])[:, 2] - 0.1].tolist()), 'nesting')
    probe('coordinate_shape', lambda a: a.update(outer_vertices=a['outer_vertices'][:-1]), 'schema')
    probe('nonfinite_coordinate', lambda a: a['outer_vertices'][0].__setitem__(2, float('nan')), 'schema')
    probe('face_removed', lambda a: a.update(faces=a['faces'][:-1]), 'schema')
    probe('wrong_units', lambda a: a.update(units='micrometres'), 'schema')
    if bridge:
        (bt, bd) = bridge
        low = submit(bt, np.asarray(bd['inner_vertices'])[:, 2] + bt['requirements']['wall_mm'], np.asarray(bd['inner_vertices'])[:, 2])
        r = all_checks(bt, low)
        out.append({'name': 'connector_area_low', 'rejected': r['checks']['connector_area']['status'] == 'FAIL', 'positive_source_task': bt['task_id'], 'value': r['checks']['connector_area']['value']})
    key = cases[0]
    name = 'task.json'
    dest = root / name
    src = BENCH / 'payload/public/tasks' / (key + '.json')
    dest.write_bytes(src.read_bytes())
    probe_b = object.__new__(Bundle)
    probe_b.root = root
    probe_b.payload = root
    probe_b.files = {name: {'sha256': sha(dest), 'bytes': dest.stat().st_size}}
    positive = probe_b.bytes(name) == src.read_bytes()
    a = json.loads(dest.read_text())
    a[0]['requirements']['clearance_mm'] = -1
    dump(dest, a)
    try:
        probe_b.bytes(name)
        reject = False
    except IntegrityError:
        reject = True
    dest.write_bytes(src.read_bytes())
    restore = probe_b.bytes(name) == src.read_bytes()
    out.append({'name': 'C02_clearance_input_change', 'positive': positive, 'rejected': reject, 'restored': restore})
    key = t['case_key']
    ref = dev_reference(key)[t['family']]
    ok = np.isfinite(ref)
    w = t['weights']
    z = ref.copy()
    z[~ok] = t['prior'][~ok]
    m0 = metric(t, z, ref)
    m1 = metric(t, z + 1, ref)
    out.append({'name': 'native_reference_height_fault', 'positive_rmse_mm': m0['rmse_mm'], 'injected_rmse_mm': m1['rmse_mm'], 'rejecting_fixed_check': 'known +1mm shift must give RMSE1 within1e-10mm', 'rejected': abs(m0['rmse_mm']) < 1e-10 and abs(m1['rmse_mm'] - 1) < 1e-10, 'external_referent': 'published native dev surface; this tests metric arithmetic, not a new measured accuracy claim'})
    good = all((x.get('rejected', x.get('passed', False)) for x in out))
    dump(ROOT / 'raw/CONTROL_REPORT.json', {'all_rejected_and_positive_restored': good, 'controls': out, 'original_files_mutated': False, 'scope': 'V3 roof; control performance never establishes physical milling/fit/strength'})
    if not good:
        raise RuntimeError('fault check failed')
    print('faults', len(out), 'PASS')
if __name__ == '__main__':
    run()
