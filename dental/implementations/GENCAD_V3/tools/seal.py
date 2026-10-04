"""Author release construction. Replay never invokes this or auto-refreezes."""
import sys, os, shutil, argparse, subprocess, importlib.metadata as metadata
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *

def write_anchor(benchmark, predictions):
    (ROOT / 'code/release_anchor.py').write_text('"""Trusted version configuration; retain its digest independently. Not participant input."""\nBENCHMARK_SHA256 = ' + repr(benchmark) + '\nPREDICTIONS_SHA256 = ' + repr(predictions) + '\n')

def runtime():
    for p in sorted((PAYLOAD / 'runtime').rglob('__pycache__'), reverse=True):
        shutil.rmtree(p)
    packages = {p: metadata.version(p) for p in ['numpy', 'scipy', 'scikit-image', 'trimesh', 'rtree', 'matplotlib', 'contourpy', 'cycler', 'fonttools', 'kiwisolver', 'pyparsing', 'python-dateutil', 'six', 'pillow', 'networkx', 'packaging', 'lazy-loader', 'imageio', 'tifffile']}
    files = {str(p.relative_to(PAYLOAD / 'runtime')): sha(p) for p in sorted((PAYLOAD / 'runtime').rglob('*')) if p.is_file()}
    dump(ROOT / 'RUNTIME_LOCK.json', dict(python_version='3.10.12', interpreter='payload/runtime/bin/python3', interpreter_sha256=sha(PAYLOAD / 'runtime/bin/python3'), platform='Linux x86_64, glibc>=2.35, unprivileged user/mount/PID/network namespaces; OS ELF loader/libraries remain declared dependencies', packages=packages, files=files))
    (ROOT / 'requirements.lock').write_text('CPython==3.10.12\n' + '\n'.join((k + '==' + v for (k, v) in packages.items())) + '\n')
    (ROOT / 'RUNTIME_EXEC.sha256').write_text(sha(PAYLOAD / 'runtime/bin/python3') + '  payload/runtime/bin/python3\n' + sha(PAYLOAD / 'runtime/bin/bwrap') + '  payload/runtime/bin/bwrap\n')

def benchmark():
    if (ROOT / 'BENCHMARK_LOCK.json').exists():
        raise RuntimeError('benchmark already sealed; no automatic refreeze')
    from sandbox import public_code
    for base in ['code', 'vendor']:
        for p in sorted((ROOT / base).rglob('__pycache__'), reverse=True):
            shutil.rmtree(p)
    public_code()
    runtime()
    cohort = read(PAYLOAD / 'private/COHORT.json')['payload']['cases']
    expected = {r['case_key'] for r in cohort if r['split'] != 'train'}
    actual = {p.stem for p in (PAYLOAD / 'public/tasks').glob('*.json')}
    if expected != actual:
        raise RuntimeError('cohort build incomplete')
    closed = ['code', 'vendor', 'review_evidence'] + ['payload/' + d for d in ['public', 'private', 'preparations', 'preparations_r3', 'crowns', 'preparations_r5', 'crowns_r5', 'preparations_r6', 'crowns_r6', 'preparations_r7', 'crowns_r7', 'runtime', 'participant_code']]
    files = {}
    for prefix in closed:
        d = ROOT / prefix
        d.mkdir(parents=True, exist_ok=True)
        for p in sorted(d.rglob('*')):
            if p.is_symlink():
                raise RuntimeError('internal symlink not allowed ' + str(p))
            if p.is_file() and p.name != 'release_anchor.py':
                files[str(p.relative_to(ROOT))] = dict(sha256=sha(p), bytes=p.stat().st_size)
    fixed = ['SCORER_PARAMETERS.json', 'RUNTIME_LOCK.json', 'RUNTIME_EXEC.sha256', 'SOURCE_REUSE.json', 'DECOMPOSITION.json', 'FROZEN_PREPARATIONS.json', 'FROZEN_PREPARATIONS_R3.json', 'FROZEN_CROWNS.json', 'FROZEN_CROWNS_R5.json', 'FROZEN_CROWNS_R6.json', 'FROZEN_CROWNS_R7.json', 'DECOMPOSITION_MAP.json', 'raw/COHORT_SUMMARY.json', 'raw/BUILD_COST.json'] + [p.name for p in sorted(ROOT.glob('PREREG_*.json'))]
    for name in fixed:
        p = ROOT / name
        files[name] = dict(sha256=sha(p), bytes=p.stat().st_size)
    lock = dict(schema='DentalGenCAD-Bench-v3', claim_type='capability', frozen_utc=now(), closed_directories=closed, files=files, trust_boundary='Evaluator plus installed release_anchor.py; retain independent digest from RELEASE_ID.txt. Every task, scene, reference, policy and executable dependency under these closed trees is bound.', task_cases=len(actual), task_count=sum((len(read(p)) for p in (PAYLOAD / 'public/tasks').glob('*.json'))))
    dump(ROOT / 'BENCHMARK_LOCK.json', lock)
    write_anchor(sha(ROOT / 'BENCHMARK_LOCK.json'), 'UNSEALED_NOT_SCORABLE')
    state('BENCHMARK_SEALED', 'Every operative input has relative path/hash', 'Generate isolated baselines, freeze, then score', locked_files=len(files))

def generate():
    from integrity import Bundle
    from sandbox import generate as run
    b = Bundle()
    dest = PAYLOAD / 'predictions'
    if (dest / 'FROZEN_PREDICTIONS.json').exists():
        raise RuntimeError('predictions already frozen')
    run(dest)
    participants = ['population', 'parametric', 'constraint_optimizer']
    cases = {Path(p).stem for (p, _) in b.tasks()}
    expected = {name + '/' + case + '.npz' for name in participants for case in cases}
    actual = {str(p.relative_to(dest)) for p in dest.rglob('*.npz')}
    if expected != actual:
        raise RuntimeError('prediction set differs from locked tasks')
    files = {name: dict(sha256=sha(dest / name), bytes=(dest / name).stat().st_size) for name in sorted(expected)}
    pred = dict(schema='DentalGenCAD-Bench-v3-predictions', frozen_utc=now(), benchmark_sha256=sha(ROOT / 'BENCHMARK_LOCK.json'), participants=participants, files=files, reference_access='public-only namespace, no scorer/source archives/private targets', controller_cost_sha256=sha(ROOT / 'raw/GENERATION_PROCESS.json'))
    dump(dest / 'FROZEN_PREDICTIONS.json', pred)
    shutil.copy2(dest / 'FROZEN_PREDICTIONS.json', ROOT / 'FROZEN_PREDICTIONS.json')
    write_anchor(sha(ROOT / 'BENCHMARK_LOCK.json'), sha(dest / 'FROZEN_PREDICTIONS.json'))
    (ROOT / 'RELEASE_ID.txt').write_text('DentalGenCAD-Bench v3\nbenchmark ' + sha(ROOT / 'BENCHMARK_LOCK.json') + '\npredictions ' + sha(dest / 'FROZEN_PREDICTIONS.json') + '\ntrusted_anchor ' + sha(ROOT / 'code/release_anchor.py') + '\n')
    state('PREDICTIONS_FROZEN', 'Public-only generation complete; predictions hashed before score', 'Run injection suite and separate evaluator')
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('phase', choices=['benchmark', 'generate'])
    a = p.parse_args()
    {'benchmark': benchmark, 'generate': generate}[a.phase]()
