from dental_release.paths import expand as _release_expand
import os, sys, subprocess, time, shlex, resource, platform, importlib.metadata
from shared import *
start = time.perf_counter()
stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
run = R / 'runs' / stamp
run.mkdir(parents=True)
DATA.mkdir(parents=True, exist_ok=True)
budget()
roots = ['LANE_X11_TOOTH_SEG', _release_expand('X21'), _release_expand('X7'), _release_expand('X18'), 'LANE_X18B_PREP_SURFACE', _release_expand('X34'), _release_expand('X38'), _release_expand('X12'), _release_expand('X8'), _release_expand('X15'), _release_expand('X20'), 'LANE_X26_DECIDABILITY', _release_expand('X33')]
files = []
for root in roots:
    p = RESULTS / root
    files.extend(p.glob('*.py'))
    files.extend((p / 'code').glob('*.py'))
files.extend((RESULTS / 'LANE_X34_STL_DESIGN_GATE/code/designgate').rglob('*.py'))
for p in [RESULTS / 'LANE_X7_BITE2TEXT/raw/models.pkl', RESULTS / 'LANE_X7_BITE2TEXT/raw/THRESHOLDS_R2.json', Path(_release_expand('@DENTAL_WORK_ROOT@/X18_crown_antagonist/donor_template.npz'))]:
    files.append(p)
for pattern in ['*global_nn.pt', '*matched_context_tree.pkl', '*SELECTION.json']:
    files.extend(Path(_release_expand('@DENTAL_WORK_ROOT@/X11')).glob(pattern))
files.extend((R / 'code').glob('*.py'))
files.extend([R / 'PREREG_R1.json', R / 'PREREG_R2.json', R / 'PREREG_R3.json', R / 'PREREG_R4.json'])
dump(run / 'ENVIRONMENT.json', dict(python=sys.executable, python_sha256=sha(Path(sys.executable).resolve()), python_version=platform.python_version(), packages={p: importlib.metadata.version(p) for p in ['numpy', 'scipy', 'trimesh', 'numba', 'nibabel', 'scikit-learn']}, threads={k: os.environ.get(k) for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']}, torch='X11 isolated torch_cpu runtime, reused read-only', plotting='system Python -s Matplotlib'))
manifest = [artifact(p) for p in sorted(set(files)) if p.is_file()]
dump(run / 'SOURCE_LOCK.json', dict(frozen_utc=now(), files=manifest))
ledger = []

def launch(script, *args, system=False):
    cmd = ['/usr/bin/python3'] + (['-s'] if system else []) + [str(R / 'code' / script), str(run), *map(str, args)]
    st = time.perf_counter()
    env = os.environ.copy()
    log = run / (Path(script).stem + '_' + ('_'.join(map(str, args)) or 'run') + '.log')
    if system:
        env['PYTHONNOUSERSITE'] = '1'
    state(script, 'Previous completed stage saved', 'Execute ' + script + ' ' + ' '.join(map(str, args)))
    with log.open('w') as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    row = dict(command=cmd, shell_display=shlex.join(cmd), wall_s=time.perf_counter() - st, exit_code=p.returncode, log=str(log))
    ledger.append(row)
    dump(run / 'COMMAND_LEDGER.json', ledger)
    print(script, *args, 'exit', p.returncode, 'seconds', round(row['wall_s'], 2), flush=True)
    if p.returncode:
        state('EXECUTION_FAILURE', row, 'Inspect exact log; do not rewrite frozen scientific gates')
        raise SystemExit(p.returncode)
    budget()
launch('ios.py')
launch('crown.py')
for c in [_release_expand('@DENTAL_CASE_A@'), _release_expand('@DENTAL_CASE_B@')]:
    launch('cbct.py', c)
launch('round2.py')
launch('round3.py')
launch('round4.py')
launch('verify.py')
launch('package.py')
launch('figure.py', system=True)
drift = [f['path'] for f in manifest if sha(f['path']) != f['sha256']]
if drift:
    dump(run / 'SOURCE_DRIFT.json', drift)
    raise RuntimeError('Read-only source drift detected')
dump(run / 'RUNTIME.json', dict(wall_s=time.perf_counter() - start, stages=ledger, peak_child_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, data_bytes=budget(), source_lock=artifact(run / 'SOURCE_LOCK.json'), source_drift=[]))
subprocess.run([sys.executable, str(R / 'code/package.py'), str(run)], check=True)
dump(R / 'LATEST_RUN.json', dict(path=str(run), runtime=artifact(run / 'RUNTIME.json'), results=artifact(R / 'results.json')))
state('FOUR_ROUNDS_COMPLETE', 'PARTIAL_CAPABILITY; computational controls PASS; clinical/manufacturing release ABSTAIN', 'Independent review, paired IOS+CBCT and actual preparation/contact/thermal measurements')
print('Complete:', R / 'README_DEMO.md', flush=True)
