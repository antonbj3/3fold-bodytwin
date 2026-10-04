"""R2 independent registration, source/output corruption and relocation tests."""
from dental_release.paths import expand as _release_expand
import copy
import datetime
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import numpy as np
from implant_safety import SafetyModule, ContractError
from implant_safety.module import load, sha
from verify import verify
ROOT = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X98-implant-safety-module'))
DATA.mkdir(exist_ok=True, parents=True)

def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def run():
    start = time.perf_counter()
    m = SafetyModule(ROOT)
    out = ROOT / 'raw'
    checks = []
    baseline = verify(out)
    q = load(ROOT / 'examples/query.json')
    s = m.sites[q['site_id']]
    base = m.query(**q)

    def reg(name, **changes):
        inp = dict(site_id=name, points_path=ROOT / s['points_local'], points_sha256=s['points_sha256'], frame=s['frame'], nominal_key='tf2_voxels_zyx', source_locator=s['external_source_bindings'][0]['archive'])
        inp.update(changes)
        return m.register_site(**inp)

    def query(name, pose):
        v = dict(q)
        v.update(site_id=name, cbct_pose=pose)
        return m.query(**v)
    reg('registered_same_source')
    r = query('registered_same_source', q['cbct_pose'])
    err = max((abs(x - y) for (x, y) in zip(r['digital_geometry']['tool_envelope']['interval_mm'], base['digital_geometry']['tool_envelope']['interval_mm'])))
    assert err <= 1e-06
    checks.append(dict(name='registered_source_matches_bundled', result='PASS', error_mm=err))
    frame = copy.deepcopy(s['frame'])
    delta = np.array([16.0, 32.0, 64.0])
    frame['origin_mm'] = delta.tolist()
    pose = copy.deepcopy(q['cbct_pose'])
    pose['frame'] = frame
    pose['entry_zyx_mm'] = (np.array(pose['entry_zyx_mm']) + delta).tolist()
    reg('translated_origin', frame=frame)
    r = query('translated_origin', pose)
    err = max((abs(x - y) for (x, y) in zip(r['digital_geometry']['tool_envelope']['interval_mm'], base['digital_geometry']['tool_envelope']['interval_mm'])))
    assert err <= 1e-06
    checks.append(dict(name='origin_and_pose_translate_together', result='PASS', error_mm=err))
    for (name, change) in [('wrong_hash', dict(points_sha256='0' * 64)), ('missing_locator', dict(source_locator='')), ('wrong_nominal_key', dict(nominal_key='absent')), ('unsupported_rotation', dict(frame=s['frame'] | dict(direction=[[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])))]:
        try:
            reg(name, **change)
        except ContractError as e:
            checks.append(dict(name=name, result='PASS', rejection=str(e)))
        else:
            raise AssertionError('registered source injection accepted: ' + name)
    single = DATA / 'single_source.npz'
    with np.load(ROOT / s['points_local']) as z:
        np.savez_compressed(single, nominal_voxels_zyx=z['tf2_voxels_zyx'])
    reg('single_mask', points_path=single, points_sha256=sha(single), nominal_key='nominal_voxels_zyx')
    rr = query('single_mask', q['cbct_pose'])
    assert rr['terms']['annotation_revision_loss']['interval_mm'] is None
    assert rr['terms']['annotation_revision_loss']['status'] == 'UNKNOWN_NO_REVISION_MEASUREMENT'
    checks.append(dict(name='single_mask_revision_is_unknown', result='PASS'))
    bad = DATA / 'tampered_source.npz'
    shutil.copyfile(single, bad)
    reg('warm_then_mutated', points_path=bad, points_sha256=sha(bad), nominal_key='nominal_voxels_zyx')
    query('warm_then_mutated', q['cbct_pose'])
    with np.load(bad) as z:
        a = z['nominal_voxels_zyx'].copy()
    a[0, 0] += 1
    np.savez_compressed(bad, nominal_voxels_zyx=a)
    try:
        query('warm_then_mutated', q['cbct_pose'])
    except ContractError as e:
        checks.append(dict(name='warm_cache_source_mutation', result='PASS', rejection=str(e)))
    else:
        raise AssertionError('warm cache ignored source mutation')
    files = ['TABLE.csv', 'QUERIES.jsonl', 'CONTROLS.json', 'CONTRACT_TESTS.json', 'SUFFICIENCY.json', 'SOURCE_REPLAY.json']
    for name in ['control_distance_plus1', 'control_gate_false', 'raw_distance_plus1', 'raw_guide_entry_plus1', 'physical_unknown_erased', 'table_distance_plus1']:
        dst = DATA / 'mutation' / name
        dst.mkdir(exist_ok=True, parents=True)
        for f in files:
            shutil.copyfile(out / f, dst / f)
        if name.startswith('control'):
            v = load(dst / 'CONTROLS.json')
            if name == 'control_distance_plus1':
                v['geometry'][0]['control_mm'] += 1
            else:
                v['control_gate'] = 'FAIL'
            dump(dst / 'CONTROLS.json', v)
        elif name.startswith('raw') or name == 'physical_unknown_erased':
            v = [json.loads(l) for l in (dst / 'QUERIES.jsonl').read_text().splitlines()]
            if name == 'raw_distance_plus1':
                v[0]['digital_geometry']['tool_envelope']['interval_mm'][0] += 1
            elif name == 'raw_guide_entry_plus1':
                v[0]['terms']['guide_combination']['entry_mm'] += 1
            else:
                v[0]['physical_safety']['status'] = 'CERTIFIED'
            (dst / 'QUERIES.jsonl').write_text(''.join((json.dumps(x) + '\n' for x in v)))
        else:
            p = dst / 'TABLE.csv'
            p.write_text(p.read_text().replace(str(base['reference_rule']['nominal_interval_mm'][0]), '999', 1))
        try:
            verify(dst)
        except AssertionError as e:
            checks.append(dict(name=name, result='PASS', rejection=str(e), artifact_directory=str(dst)))
        else:
            raise AssertionError('FINAL verifier accepted injected error: ' + name)
    portable = DATA / 'portable_bundle'
    portable.mkdir(exist_ok=True)
    paths = [ROOT / 'SOURCE_MANIFEST.json', ROOT / 'PREREG_R1.json', ROOT / 'PREREG_R1.json.sha256', ROOT / 'PREREG_R2.json', ROOT / 'PREREG_R2.json.sha256', ROOT / 'FROZEN_PREDICTIONS.json', ROOT / 'FROZEN_PREDICTIONS.json.sha256', ROOT / 'FROZEN_PREDICTIONS_R2.json', ROOT / 'FROZEN_PREDICTIONS_R2.json.sha256', ROOT / 'PREREG_R3.json', ROOT / 'PREREG_R3.json.sha256', ROOT / 'CODE_MANIFEST.json', ROOT / 'make_demo.py', ROOT / 'validation.py', ROOT / 'verify.py', ROOT / 'plot_demo.py', ROOT / 'run_all.sh', ROOT / 'release_checks.py']
    paths += list((ROOT / 'implant_safety').rglob('*.py')) + list((ROOT / 'data').iterdir())
    paths += [ROOT / r['local'] for r in load(ROOT / 'SOURCE_MANIFEST.json') if r.get('local')]
    for p in paths:
        d = portable / p.relative_to(ROOT)
        d.parent.mkdir(exist_ok=True, parents=True)
        shutil.copyfile(p, d)
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[k] = '1'
    before = sha(portable / 'FROZEN_PREDICTIONS_R2.json')
    with (DATA / 'PORTABLE_RUN.log').open('w') as log:
        run = subprocess.run(['bash', 'run_all.sh'], cwd=portable, env=env, stdout=log, stderr=subprocess.STDOUT)
    assert run.returncode == 0, 'relocated one-command demo failed: ' + str(DATA / 'PORTABLE_RUN.log')
    assert sha(portable / 'FROZEN_PREDICTIONS_R2.json') == before
    result = subprocess.run([sys.executable, 'verify.py'], cwd=portable, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    checks.append(dict(name='relocated_one_command_and_final_gate', result='PASS', frozen_prediction_identity_error=0, log_sha256=sha(DATA / 'PORTABLE_RUN.log')))
    total = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    assert total < 3000000000
    res = dict(round='R3', claim_type='capability', outcome='REGISTERED_LOCAL_SOURCE_AND_PORTABLE_API_AVAILABLE; PHYSICAL_UNKNOWN', final_gate=baseline, checks=checks, all_passed=all((x['result'] == 'PASS' for x in checks)), external_referent=dict(kind='published_dataset', locator='SOURCE_MANIFEST.json; frozen source-bound R1 and R2 predictions', compared_quantity='Minimum gap after explicit local registration and joint origin/pose translation', refutes_us=False), cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, threads=1, gpu=False, scratch_bytes=total), mutations=dict(final_verifier_rejections=6, registered_source_rejections=5), scratch_artifacts=[dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for p in [single, bad, DATA / 'PORTABLE_RUN.log']], semantic_correction='R1 broad guide scenario class UNRESOLVED_NUMERIC_BRACKET now INTERVAL_CROSSES_REFERENCE. Same endpoints, tolerance and gate; original frozen prediction remains immutable.')
    dump(ROOT / 'rounds/R3.json', res)
    print(json.dumps({k: v for (k, v) in res.items() if k not in ['checks', 'scratch_artifacts']}))
    return res
if __name__ == '__main__':
    run()
