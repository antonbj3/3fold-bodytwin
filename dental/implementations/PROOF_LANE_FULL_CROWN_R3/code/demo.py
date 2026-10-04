from common_r3 import *
import subprocess, time, resource

def stripped(x):
    if isinstance(x, dict):
        return {k: stripped(v) for (k, v) in x.items() if k not in ['score_seconds', 'score_peak_rss_MiB']}
    if isinstance(x, list):
        return [stripped(v) for v in x]
    return x

def main():
    start = time.perf_counter()
    rows = []
    for name in ['FROZEN_GENERATOR_A.json', 'FROZEN_GENERATOR_B.json', 'FROZEN_GENERATOR_C.json', 'FROZEN_MEMORY_REPAIR.json']:
        for (rel, h) in read(ROOT / name)['files'].items():
            assert sha(ROOT / 'code' / rel) == h, rel
    for (tag, folder) in [('A', 'inputs'), ('C', 'inputs_C')]:
        for (rel, m) in read(ROOT / f'FROZEN_ACQUISITION_{tag}.json')['files'].items():
            assert sha(DATA / folder / rel) == m['sha256'], rel
    for r in read(ROOT / 'results.json')['source_manifest'].values():
        assert sha(r['path']) == r['sha256'], r['path']
    from regenerate_selected import run as regenerate
    regenerate()
    for p in ROOT.glob('*.sha256'):
        expected = p.read_text().split()[0]
        target = ROOT / p.read_text().split()[1]
        assert sha(target) == expected, str(p)
    for tag in ['A', 'B', 'C']:
        frozen = read(ROOT / f'FROZEN_PREDICTIONS_{tag}.json')
        for (rel, m) in frozen['files'].items():
            assert sha(DATA / tag / rel) == m['sha256'], rel
        expected = read(ROOT / 'rounds' / f'{tag}.json')['rows']
        for r in expected:
            ident = tag + '_' + r['key'] + '__' + r['participant']
            dest = ROOT / 'raw/recheck' / (ident + '.json')
            dest.parent.mkdir(exist_ok=True)
            cmd = ['/usr/bin/time', '-v', '-o', str(dest.with_suffix('.memory.txt')), PYTHON, str(ROOT / 'code/score_worker.py'), tag, r['key'], r['participant'], str(dest)]
            proc = subprocess.run(cmd)
            if proc.returncode:
                raise RuntimeError('Recheck process failed:' + ident)
            observed = read(dest)
            same = stripped(observed) == stripped(r)
            if not same:
                dump(ROOT / 'raw/RECHECK_MISMATCH.json', dict(expected=r, observed=observed))
                raise AssertionError('scientific score changed:' + ident)
            rows.append(dict(task=ident, identical_scientific_values=True, peak_rss_MiB=observed.get('score_peak_rss_MiB')))
    from validation import run as vc
    from validate_memory import rows as memory_rows
    vc()
    subprocess.run(['/usr/bin/time', '-v', '-o', str(ROOT / 'raw/SCORER_CONTROLS_MEMORY.txt'), PYTHON, str(ROOT / 'code/scorer_controls.py')], check=True)
    from lab_export_r3 import run as lab
    lab()
    from report_r3 import run as report
    report()
    dump(ROOT / 'raw/DEMO_RECHECK.json', dict(status='PASS', recomputed=len(rows), rows=rows, seconds=time.perf_counter() - start, parent_peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Recomputed all72 construction records and actual exported crowns; three selected generators actually re-executed separately; no physical test claimed'))
    print('DEMO RECHECK PASS', len(rows), flush=True)
if __name__ == '__main__':
    main()
