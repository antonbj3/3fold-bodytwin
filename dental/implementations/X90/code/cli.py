"""Portable consumer entry point; data/runtime locations are external configuration."""
import argparse
import datetime
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
import hashlib
from common import ROOT, read, dump, sha, freeze
from render import render, verify_report

def execute(args):
    config = read(args.site_config)
    config_root = Path(args.site_config).resolve().parent

    def configured(path):
        p = Path(path)
        return p if p.is_absolute() else config_root / p
    package = (Path(args.package) if args.package else configured(config['package'])).resolve()
    nerve = configured(config['nerve']).resolve()
    patient = args.patient
    if patient.isdigit():
        patient = 'Bits2Bites:' + patient
    name = patient.replace(':', '_')
    if args.manifest:
        manifest = Path(args.manifest).resolve()
    else:
        manifest = ROOT / 'inputs' / (name + '.json')
        if not manifest.exists():
            raise SystemExit('No case manifest. Use --manifest; examples are in inputs/.')
    if read(manifest)['patient'] != patient:
        raise ValueError('CLI_SUBJECT_MANIFEST_MISMATCH')
    epoch = ROOT / 'inputs/BUILD_EPOCH.json'
    if not epoch.exists():
        freeze(epoch, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), purpose='Immutable report build epoch; never a claim of blind physical prediction'))
    expected = read(ROOT / 'inputs/TOOL_LOCK.json')
    for row in expected['files']:
        root = package if row['root'] == 'package' else nerve
        if sha(root / row['path']) != row['sha256']:
            raise ValueError('TOOL_SOURCE_DRIFT: ' + row['path'])
    prereg = ROOT / 'PREREG_R1.json'
    if sha(prereg) != (ROOT / 'PREREG_R1.json.sha256').read_text().split()[0]:
        raise ValueError('PREREG_DRIFT')
    out = Path(args.out or ROOT / 'reports' / name).resolve()
    out.mkdir(parents=True, exist_ok=True)
    code_hash = hashlib.sha256(''.join((sha(ROOT / 'code' / n) for n in ['common.py', 'operators.py', 'case_worker.py'])).encode()).hexdigest()
    marker_key = hashlib.sha256((patient + code_hash + sha(manifest) + sha(epoch)).encode()).hexdigest()
    marker = ROOT / 'inputs/frozen_case_records' / (marker_key + '.json')
    marker.parent.mkdir(parents=True, exist_ok=True)
    prior_prediction = out / 'FROZEN_PREDICTIONS.json'
    if not marker.exists() and prior_prediction.exists():
        previous = read(prior_prediction)
        if previous.get('adapter_code_sha256') != code_hash or previous.get('input_manifest_sha256') != sha(manifest):
            raise ValueError('EXISTING_FREEZE_VERSION_MISMATCH: use a new output directory')
        if sha(prior_prediction) != prior_prediction.with_suffix('.json.sha256').read_text().strip():
            raise ValueError('EXISTING_FREEZE_HASH_MISMATCH')
        freeze(marker, previous)
    context = read(epoch)
    if marker.exists():
        context['per_case_frozen_utc'] = read(marker)['frozen_utc']
    prediction_epoch = out / 'PREDICTION_EPOCH.json'
    dump(prediction_epoch, context)
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg', NUMBA_CACHE_DIR=str(out / 'numba_cache'))
    env['MPLCONFIGDIR'] = str(out / 'plot_cache')
    if config.get('runtime_deps'):
        env['PYTHONPATH'] = str(configured(config['runtime_deps']))
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS']:
        env[key] = '4'
    cmd = [str(configured(config.get('python', sys.executable))), str(ROOT / 'code/case_worker.py'), '--manifest', str(manifest), '--out', str(out), '--package', str(package), '--nerve', str(nerve), '--dataset-root', str(configured(config['dataset_root'])), '--labels-root', str(configured(config['labels_root'])), '--epoch', str(prediction_epoch)]
    cmd = ['bwrap', '--ro-bind', '/', '/', '--bind', str(out), str(out), '--dev', '/dev', '--proc', '/proc', '--tmpfs', '/tmp', '--unshare-net', '--die-with-parent', '--'] + cmd
    t0 = time.perf_counter()
    with (out / 'execution.log').open('w') as f:
        proc = subprocess.run(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
    dump(out / 'EXECUTION.json', dict(command=cmd, shell_display=shlex.join(cmd), exit_code=proc.returncode, wall_seconds=time.perf_counter() - t0, network='UNSHARED', upstream_writable=False))
    if proc.returncode:
        raise SystemExit('Case calculation failed; see ' + str(out / 'execution.log'))
    freeze(marker, read(out / 'FROZEN_PREDICTIONS.json'))
    render(out)
    verify_report(out)
    print(out / 'report.html')
    return out

def main():
    ap = argparse.ArgumentParser(description='Offline patient-bound dental laboratory research report')
    ap.add_argument('patient', help='Bits2Bites:13 (or 13); explicit case manifest defines available modalities')
    ap.add_argument('--out')
    ap.add_argument('--manifest')
    ap.add_argument('--package', help='Reviewed DEMO48_PACKAGE root; tools must match TOOL_LOCK')
    ap.add_argument('--site-config', default=str(ROOT / 'SITE_CONFIG.json'))
    execute(ap.parse_args())
if __name__ == '__main__':
    main()
