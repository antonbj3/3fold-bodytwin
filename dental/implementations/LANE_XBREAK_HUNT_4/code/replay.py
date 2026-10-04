"""Replay own small experiment inputs into a fresh directory, preserving first results."""
import json, hashlib, subprocess, datetime, os, sys, shutil, time
from pathlib import Path

def scientific_projection(d):
    if isinstance(d, dict):
        return {k: scientific_projection(v) for (k, v) in d.items() if k not in ['cost', 'analytic_query_seconds']}
    if isinstance(d, list):
        return [scientific_projection(v) for v in d]
    return d

def compare(a, b, path=''):
    if isinstance(a, dict):
        assert set(a) == set(b), (path, 'keys')
        for k in a:
            compare(a[k], b[k], path + '/' + k)
    elif isinstance(a, list):
        assert len(a) == len(b), (path, 'length')
        for (i, (x, y)) in enumerate(zip(a, b)):
            compare(x, y, path + '/' + str(i))
    elif isinstance(a, (int, float)) and (not isinstance(a, bool)):
        assert abs(a - b) <= 1e-07, (path, a, b)
    else:
        assert a == b, (path, a, b)

def main():
    root = Path(__file__).resolve().parent.parent
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    target = root / 'replay' / stamp
    for p in root.glob('PREREG*.json'):
        assert hashlib.sha256(p.read_bytes()).hexdigest() == p.with_suffix('.sha256').read_text().strip(), p
    x = json.loads((root / 'raw/R1_SOURCE_COUNTS.json').read_text())
    s = json.loads((root / 'raw/R5_TRIAL_INPUTS.json').read_text())
    assert hashlib.sha256((root / 'sources/pbrn.pdf').read_bytes()).hexdigest() == x['source']['source_sha256']
    assert hashlib.sha256((root / 'sources/PMC11629050.html').read_bytes()).hexdigest() == s['source_sha256']
    target.mkdir(parents=True)
    (target / 'raw').mkdir()
    shutil.copytree(root / 'code', target / 'code', ignore=shutil.ignore_patterns('__pycache__'))
    for p in root.glob('PREREG*'):
        shutil.copy2(p, target / p.name)
    for name in ['R1_SOURCE_COUNTS.json', 'R3_PAIRED_REPORTS.json', 'R5_TRIAL_INPUTS.json']:
        shutil.copy2(root / 'raw' / name, target / 'raw' / name)
    env = os.environ.copy()
    for var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[var] = '4'
    results = []
    for i in range(1, 8):
        t = time.perf_counter()
        run = subprocess.run([sys.executable, 'code/r' + str(i) + '.py'], cwd=target, env=env, capture_output=True, text=True)
        (target / f'R{i}_stdout.txt').write_text(run.stdout)
        (target / f'R{i}_stderr.txt').write_text(run.stderr)
        if run.returncode:
            raise RuntimeError(f'R{i} failed: ' + run.stderr)
        a = json.loads((root / f'raw/R{i}_RESULTS.json').read_text())
        b = json.loads((target / f'raw/R{i}_RESULTS.json').read_text())
        compare(scientific_projection(a), scientific_projection(b))
        results.append({'round': f'R{i}', 'scientific_compare': 'PASS', 'full_process_seconds': time.perf_counter() - t, 'raw_sha256': hashlib.sha256((target / f'raw/R{i}_RESULTS.json').read_bytes()).hexdigest()})
    summary = {'directory': str(target), 'results': results, 'source_and_protocol_hashes': 'PASS', 'all_scientific_outcomes': 'PASS', 'total_full_process_seconds': sum((r['full_process_seconds'] for r in results))}
    (target / 'REPLAY_COMPARE.json').write_text(json.dumps(summary, indent=2) + '\n')
    (root / 'REPLAY_LATEST.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    main()
