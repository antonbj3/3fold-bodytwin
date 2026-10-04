"""One-command replay always writes to a new private audit copy."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, sys, shutil, datetime, json, hashlib, subprocess
ROOT = Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    lock = json.loads((ROOT / 'RELEASE_LOCK.json').read_text())
    for (p, h) in lock['code_and_contracts'].items():
        if sha(ROOT / p) != h:
            raise RuntimeError('Executed source/copy differs from release: ' + p)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    base = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/X95-crown-prep-envelope')) / stamp
    base.mkdir(parents=True)
    copy = base / 'package'
    shutil.copytree(ROOT, copy, ignore=shutil.ignore_patterns('__pycache__', '.git', 'prior_runs', 'AUDIT_REPLAY_RECEIPT.json'))
    env = os.environ.copy()
    env.update(X95_DATA_DIR=str(base / 'data'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    cmd = [sys.executable, str(copy / 'code/replay_worker.py')]
    if '--full-source' in sys.argv:
        cmd.append('--full-source')
    run = subprocess.run(cmd, env=env, cwd=copy)
    receipt = {'audit_copy': str(copy), 'private_data': str(base / 'data'), 'exit_code': run.returncode, 'copied_release_lock_sha256': sha(copy / 'RELEASE_LOCK.json'), 'result_sha256': sha(copy / 'results.json') if (copy / 'results.json').exists() else None, 'replay_comparison': json.loads((copy / 'raw/REPLAY_COMPARISON.json').read_text()) if (copy / 'raw/REPLAY_COMPARISON.json').exists() else None}
    (ROOT / 'AUDIT_REPLAY_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))
    sys.exit(run.returncode)
if __name__ == '__main__':
    main()
