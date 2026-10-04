import json, time, datetime, resource, sys, subprocess, platform
from pathlib import Path
from lab_port import ROOT, sha
from verify import source_check, run as verify
from model import run_full
from round2 import run as run_reduced
from lab_port import run as run_port
from report import run as report

def main():
    start = time.monotonic()
    lock = json.loads((ROOT / 'INPUT_LOCK.json').read_text())
    for r in lock['files']:
        p = Path(r['path'])
        if sha(p) != r['sha256']:
            raise ValueError('INPUT_HASH_CHANGED:' + str(p))
    rows = json.loads((ROOT / 'raw/RETENTION_GROUPS.json').read_text())
    source_check(rows)
    run_full(rows)
    run_reduced(rows)
    run_port()
    verify(rows)
    out = report()
    receipt = {'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'wall_s': time.monotonic() - start, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'interpreter': sys.executable, 'python_version': platform.python_version(), 'threads_limit': 4, 'GPU': False, 'physical_measurements': 0, 'results_sha256': sha(ROOT / 'results.json'), 'input_lock_sha256': sha(ROOT / 'INPUT_LOCK.json')}
    (ROOT / 'RUN_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('DONE: source-bound retention endpoint; absolute/aging calibration UNKNOWN; three rounds replayed; ' + str(round(receipt['wall_s'], 2)) + ' s')
if __name__ == '__main__':
    main()
