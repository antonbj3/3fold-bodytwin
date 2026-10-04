from common import *
import argparse, subprocess, platform, importlib.metadata, resource, shutil

def run():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    start = time.perf_counter()
    commands = []
    initial = ROOT / 'history/initial_measurements'
    initial.mkdir(exist_ok=True, parents=True)
    for name in ['LOCAL_PAIRS.json', 'ANNOTATED_PAIRS.json', 'RESCORE.json']:
        dest = initial / name
        if not dest.exists():
            shutil.copyfile(ROOT / 'raw' / name, dest)
    assert sha(initial / 'ANNOTATED_PAIRS.json') == read(ROOT / 'FROZEN_RESCORE_PREDICTIONS.json')['threshold_source_sha256']
    for name in ['PREREG_LOCAL', 'PREREG_ANNOTATED', 'PREREG_RESCORE', 'FROZEN_RESCORE_PREDICTIONS']:
        assert sha(ROOT / (name + '.json')) == (ROOT / (name + '.sha256')).read_text().split()[0]
    for (p, h) in read(ROOT / 'PREREG_LOCAL.json')['input_hashes'].items():
        assert sha(p) == h, p
    for (p, h) in read(ROOT / 'PREREG_RESCORE.json')['input_hashes'].items():
        assert sha(p) == h, p
    for src in read(ROOT / 'raw/ANNOTATED_PAIRS.json')['sources']:
        assert sha(src['path']) == src['obj_sha256']
        assert sha(src['labels_path']) == src['labels_sha256']
    packages = {p: importlib.metadata.version(p) for p in ['numpy', 'scipy', 'trimesh', 'rtree', 'matplotlib', 'beautifulsoup4']}
    dump(ROOT / 'raw/RUNTIME.json', dict(python=sys.version, executable=sys.executable, platform=platform.platform(), packages=packages, threads={k: os.environ.get(k) for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']}, clean_install_tested=False))
    steps = [] if args.verify else ['local_pairs', 'annotated_pairs', 'rescore']
    steps += ['component_audit', 'evidence', 'validation', 'report']
    for step in steps:
        cmd = [PYTHON, str(ROOT / 'code' / (step + '.py'))]
        t = time.perf_counter()
        with (ROOT / 'raw' / ('demo_' + step + '.log')).open('w') as log:
            proc = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(argv=cmd, returncode=proc.returncode, seconds=time.perf_counter() - t))
        dump(ROOT / 'raw/DEMO.json', dict(mode='verify' if args.verify else 'full', commands=commands, seconds=time.perf_counter() - start))
        if proc.returncode:
            raise RuntimeError('Failed ' + step + '; inspect raw/demo_' + step + '.log')
        print(step + ' PASS', flush=True)
    parity = {}
    for name in ['LOCAL_PAIRS.json', 'ANNOTATED_PAIRS.json']:
        before = read(initial / name)['rows']
        after = read(ROOT / 'raw' / name)['rows']
        assert len(before) == len(after)
        errors = [abs(x['metrics']['8192']['p95_mm'] - y['metrics']['8192']['p95_mm']) for (x, y) in zip(before, after)]
        parity[name] = max(errors)
        assert max(errors) <= 1e-09
    b = read(initial / 'RESCORE.json')['rows']
    a = read(ROOT / 'raw/RESCORE.json')['rows']
    assert len(a) == len(b) == 90
    parity['rescore_max_shape_p95_difference_mm'] = max((abs(x['oracle_aligned']['8192']['p95_mm'] - y['oracle_aligned']['8192']['p95_mm']) for (x, y) in zip(a, b) if x['status'] == 'RESCORED'))
    assert parity['rescore_max_shape_p95_difference_mm'] <= 1e-09
    assert all((x.get('shape_envelope_pass') == y.get('shape_envelope_pass') and x.get('contact_gates') == y.get('contact_gates') for (x, y) in zip(a, b)))
    dump(ROOT / 'raw/REPLAY_PARITY.json', dict(mode='verify_only_no_remeasurement' if args.verify else 'full_remeasurement', maximum_errors=parity, all_shape_and_contact_gates_identical=True))
    print('All stages PASS; see RESULTS.md and figures/diagnosis.png', flush=True)
if __name__ == '__main__':
    run()
