from core import *
import argparse, shutil

def verify():
    st = time.perf_counter()
    checks = []

    def check(name, ok):
        checks.append(dict(check=name, passed=bool(ok)))
        if not ok:
            raise ValueError('FAILED ' + name)
    for p in sorted(R.glob('PREREG*.json')):
        check(p.name, sha(p) == p.with_suffix('.sha256').read_text().split()[0])
    inputs = read(R / 'INPUTS.json')['records']
    for r in inputs:
        for field in ['public', 'private']:
            check(r['key'] + '_' + field, sha(r[field + '_path']) == r[field + '_sha256'])
        check(r['key'] + '_mirrored_outer', sha(r['outer']['mesh_path']) == r['outer']['mesh_sha256'])
    rows = []
    for tag in ['R1', 'R2', 'R3']:
        p = R / f'FROZEN_PREDICTIONS_{tag}.json'
        check(p.name, sha(p) == p.with_suffix('.sha256').read_text().split()[0])
        fr = read(p)
        check(tag + '_same36', len(fr['rows']) == 36)
        if tag != 'R1':
            rows.extend(read(R / f'RESULTS_{tag}.json')['rows'])
    for row in rows:
        if row['status'] != 'GENERATED':
            continue
        check(row['key'] + '_' + row['round'] + '_mesh', sha(row['mesh_path']) == row['mesh_sha256'])
        m = load(row['mesh_path'])
        P = Distance(m['prep_vertices'], m['prep_faces'])
        inn = m['vertices'][m['faces'][m['roles'] == 1]]
        gap = P.query(sample(inn, 128))[0]
        check(row['key'] + '_actual_gap', bool(abs(gap - 0.05).max() <= 0.01))
        certlist = dict(row['certificates'])
        if row.get('actual_mesh_wall_proof'):
            certlist['actual_emitted_outer'] = row['actual_mesh_wall_proof']['exact_outer_separation']
        if row.get('native_containment_obstruction'):
            certlist['outside_witness_in_Q'] = row['native_containment_obstruction']['Q_membership']
        for (name, cc) in certlist.items():
            check(name + '_hash', sha(cc['path']) == cc['sha256'])
            c = read(cc['path'])
            check(name + '_surface_hash', sha(c['support_path']) == c['support_sha256'])
            check(name + '_queries_hash', sha(c['query_path']) == c['query_sha256'])
            exe = 'exact_surface' if 'inside_checked' in cc else 'exact_clearance'
            cmd = [str(D / exe), c['support_path'], c['query_path']]
            if cc.get('inside_checked') is False:
                cmd.append('surface')
            z = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
            check(name + '_backend_exit', z.returncode == 0)
            computed = json.loads(z.stdout)['rows']
            check(name + '_all_exact_replay', computed == c['rows'])
        for field in ['stl', 'three_mf']:
            if row.get('exports'):
                check(row['key'] + '_' + field, sha(row['exports'][field]) == row['exports'][field + '_sha256'])
    print('Validated', len(rows), 'case/material rows; exact digital certificates replayed.', flush=True)
    from controls import run
    controls = run()
    check('fault_controls', controls['all_pass'])
    out = dict(status='PASS', checks=len(checks), seconds=time.perf_counter() - st, all_pass=all((x['passed'] for x in checks)), complete_crowns=read(R / 'results.json')['complete_crowns_per_tooth_type'], control_passes=controls['passed'], control_count=controls['count'], detail=checks, peak_process_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(R / 'raw/REPLAY_VALIDATION.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'detail'}, indent=2))

def recompute():
    size = sum((p.stat().st_size for p in D.rglob('*') if p.is_file()))
    if size > 1.5 * 1024 ** 3:
        raise ValueError('Recompute would risk3GB lane cap; preserved data not deleted automatically')
    runroot = D / ('replay_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S'))
    runroot.mkdir()
    (runroot / 'raw').mkdir()
    (runroot / 'figures').mkdir()
    art = runroot / 'artifacts'
    art.mkdir()
    for p in list(R.glob('PREREG*')) + list(R.glob('DECOMPOSITION*')) + [R / 'INPUTS.json', R / 'INPUTS.sha256']:
        shutil.copy2(p, runroot / p.name)
    for n in ['PROOF.md', 'LITERATURE.md', 'EXISTING_WORK.md']:
        shutil.copy2(R / n, runroot / n)
    for n in ['INITIAL_LOCAL_MEASUREMENT.json', 'CERVICAL_FRAME_DIAGNOSTIC.json']:
        shutil.copy2(R / 'raw' / n, runroot / 'raw' / n)
    for n in ['exact_clearance', 'exact_surface']:
        (art / n).symlink_to(D / n)
    env = dict(os.environ, DENT_FIX_RUN_ROOT=str(runroot), DENT_FIX_DATA_ROOT=str(art))
    code = Path(__file__).resolve().parent
    for (script, args) in [('generate.py', ['R1']), ('regularize.py', []), ('sleeve.py', []), ('evaluate.py', ['R1']), ('evaluate.py', ['R2']), ('evaluate.py', ['R3']), ('parity.py', []), ('controls.py', []), ('report.py', []), ('figure.py', [])]:
        print('RECOMPUTE', script, args, flush=True)
        subprocess.run([sys.executable, str(code / script), *args], env=env, check=True)
    print('Fresh run preserved in', runroot)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--recompute', action='store_true')
    a = p.parse_args()
    recompute() if a.recompute else verify()
