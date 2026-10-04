"""Run an external Python generator in the same public-only sandbox and score its frozen output."""
import argparse, subprocess, sys, re
from pathlib import Path
from .common import *
from .sandbox import run

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plugin', required=True)
    ap.add_argument('--name', required=True)
    a = ap.parse_args()
    if not re.fullmatch('[A-Za-z][A-Za-z0-9_]{0,39}', a.name):
        raise ValueError('safe unique participant name required')
    (path, func) = a.plugin.rsplit(':', 1)
    path = Path(path).resolve()
    tag = 'external_' + a.name
    if (ROOT / f'FROZEN_PREDICTIONS_{tag}.json').exists():
        raise ValueError('participant already frozen; choose a new version name')
    dest = DATA / 'external_predictions' / a.name
    run(a.name, output=dest, plugin=(str(path), func))
    files = {str(p.relative_to(DATA)): sha(p) for p in sorted(dest.glob('*.json')) if p.name != '_cost.json'}
    freeze(ROOT / f'FROZEN_PREDICTIONS_{tag}.json', dict(participants=[a.name], prediction_dirs={a.name: str(dest.relative_to(DATA))}, files=files, plugin_sha256=sha(path), benchmark_sha256=sha(ROOT / 'data/BENCHMARK_LOCK.json')))
    p = subprocess.run([sys.executable, '-m', 'gencad_bench_v2.score_cli', '--round', tag], cwd=ROOT)
    raise SystemExit(p.returncode)
if __name__ == '__main__':
    main()
