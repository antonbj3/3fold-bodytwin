"""Fresh computation in its own data/result directories; original freezes never overwritten."""
from common import *
import subprocess, shutil, time, resource

def run():
    start = time.perf_counter()
    release = read(ROOT / 'RELEASE_LOCK.json')
    for (rel, item) in release['files'].items():
        if sha(ROOT / rel) != item['sha256']:
            raise ValueError('RELEASE_HASH_MISMATCH ' + rel)
    for (name, item) in read(ROOT / 'INPUT_LOCK.json')['files'].items():
        if sha(name) != item['sha256']:
            raise ValueError('SOURCE_HASH_MISMATCH ' + name)
    anchor = ROOT / 'code/plugin.py'
    fault = anchor.read_bytes() + b'\n# injected fault\n'
    fault_rejected = hashlib.sha256(fault).hexdigest() != release['files']['code/plugin.py']['sha256']
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    dest = DATA / 'replays' / stamp
    dest.mkdir(parents=True)
    resultdir = dest / 'report'
    resultdir.mkdir()
    (resultdir / 'raw').mkdir()
    (resultdir / 'code').symlink_to(ROOT / 'code', target_is_directory=True)
    for p in ROOT.glob('PREREG*.json'):
        shutil.copy2(p, resultdir / p.name)
    for name in ['PANEL_SELECTION.json', 'INPUT_LOCK.json']:
        shutil.copy2(ROOT / name, resultdir / name)
    env = dict(os.environ, X60_RUN_ROOT=str(resultdir), X60_RUN_DATA=str(dest / 'data'), PYTHONDONTWRITEBYTECODE='1')
    (dest / 'data').mkdir()
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[key] = '1'
    commands = []
    for name in ['round1.py', 'generate.py', 'score.py', 'whole.py', 'round4.py', 'round5.py', 'export_roofs.py', 'verify.py', 'freeze_lab.py']:
        cmd = [sys.executable, str(ROOT / 'code' / name)]
        tic = time.perf_counter()
        log = dest / (name + '.log')
        print('replay', name, flush=True)
        with log.open('w') as f:
            p = subprocess.run(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
        commands.append(dict(argv=cmd, exit_code=p.returncode, seconds=time.perf_counter() - tic, log=str(log)))
        if p.returncode:
            dump(ROOT / 'raw/REPLAY_FAILED.json', dict(commands=commands, preserved_directory=str(dest)))
            raise RuntimeError('fresh replay failure: ' + str(log))
        budget()
    errors = []
    arrays = 0
    for tag in ['volume_1.0', 'volume_1.1']:
        for case in read(ROOT / 'PANEL_SELECTION.json')['selected']:
            name = case['case_key'] + '.npz'
            with np.load(DATA / 'predictions' / tag / name) as a, np.load(dest / 'data/predictions' / tag / name) as b:
                if set(a) != set(b):
                    errors.append(tag + '/' + name + ' keys')
                for k in a:
                    arrays += 1
                    if not np.array_equal(a[k], b[k]):
                        errors.append(tag + '/' + name + '/' + k)
    geometry = []
    for sub in ['WHOLE_RAW.npz', 'ROBUST_RAW.npz', 'SURFACE_TIE_RAW.npz']:
        with np.load(DATA / sub) as a, np.load(dest / 'data' / sub) as b:
            for k in ['vertices', 'optimized_vertices', 'coefficients_mm']:
                if k in a:
                    diff = float(np.max(abs(a[k] - b[k])))
                    geometry.append(dict(artifact=sub, key=k, max_absolute_difference=diff, pass_gate=diff <= 1e-10))
    r1 = read(resultdir / 'raw/R1.json')
    ver = read(resultdir / 'raw/VERIFICATION.json')
    lc = read(resultdir / 'raw/LAB_CONTROLS.json')
    passed = not errors and all((r['pass_gate'] for r in geometry)) and ver['all_faults_rejected'] and ver['all_exports_pass'] and ver['all_export_faults_rejected'] and all((r['wrong_force_rejected'] and r['missing_support_refused'] for r in lc)) and (r1['sufficiency']['identity_error'] == 0) and fault_rejected
    out = dict(pass_gate=passed, fresh_recomputation=True, source_workspaces_readonly=True, network_used=False, isolated_generator=True, geometry_array_count=arrays, geometry_differences=geometry, errors=errors, commands=commands, release_fault_rejected=fault_rejected, summary_identity_error=r1['sufficiency']['identity_error'], all_valid_coordinate_faults_rejected=ver['all_faults_rejected'], all_export_faults_rejected=ver['all_export_faults_rejected'], calibrated_force05_N=None, physical_measurement=False, seconds=time.perf_counter() - start, peak_children_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, replay_directory=str(dest), resource_bytes=budget(), review_state='PENDING_INDEPENDENT_REVIEW')
    dump(ROOT / 'raw/REPLAY_VERIFICATION.json', out)
    print('fresh replay PASS' if passed else 'fresh replay FAIL', 'seconds', round(out['seconds'], 1), 'peak RSS MiB', round(out['peak_children_RSS_MiB'], 1), flush=True)
    if not passed:
        raise ValueError('replay verification failed')
if __name__ == '__main__':
    run()
