from dental_release.paths import expand as _release_expand
from common import *
import subprocess, resource
PY = _release_expand('@DENTAL_PYTHON@')

def main():
    start = time.perf_counter()
    checks = 0
    for name in ['INPUT_LOCK_R1.json', 'INPUT_LOCK_R2.json', 'INPUT_LOCK_R3.json']:
        for (p, r) in read(ROOT / name)['files'].items():
            if sha(p) != r['sha256']:
                raise RuntimeError('SOURCE_DRIFT ' + p)
            checks += 1
    for p in ROOT.glob('PREREG_*.json'):
        if sha(p) != p.with_suffix('.sha256').read_text().split()[0]:
            raise RuntimeError('PREREG_DRIFT ' + str(p))
    if sha(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext/PMC10958416.xml')) != read(ROOT / 'PREREG_R4.json')['source_sha256']:
        raise RuntimeError('SOURCE_TABLE_DRIFT')
    for r in read(ROOT / 'FROZEN_MEASUREMENT_PREDICTIONS.json')['pair']:
        for field in ['mesh', 'stl']:
            if sha(r[field]) != r[field + '_sha256']:
                raise RuntimeError('EXPORT_DRIFT ' + r[field])
    log = ROOT / 'raw/RUN_LEDGER.jsonl'
    steps = [(PY, 'verify_contact.py'), (PY, 'r1.py'), (PY, 'r2.py'), ('/usr/bin/python3', 'r3.py'), ('/usr/bin/python3', 'verify_load.py'), (PY, 'r4.py'), (PY, 'r5.py'), (PY, 'verify_outputs.py'), (PY, 'report.py')]
    for (python, script) in steps:
        s = time.perf_counter()
        cmd = [python, '-B', str(ROOT / 'code' / script)]
        print('Executing', script, flush=True)
        with (ROOT / 'raw' / ('demo_' + script + '.log')).open('w') as out:
            r = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, cwd=ROOT)
        rec = dict(utc=now(), argv=cmd, returncode=r.returncode, seconds=time.perf_counter() - s)
        with log.open('a') as f:
            f.write(json.dumps(rec) + '\n')
        if r.returncode:
            raise RuntimeError('STEP_FAILED ' + script + '; see raw/demo_' + script + '.log')
    cost = dict(seconds=time.perf_counter() - start, peak_child_rss_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, source_files_checked=checks, threads_max=4, GPU=False, preparation='Included current hash/read; inherited preparation UNKNOWN', fit='none', discovery='Implementation/source search time UNKNOWN; no method speedup claim', validation='included controls', queries='all rounds freshly executed', fallback='Explicit unknowns/rejected rows; no physical substitute', upstream_training_and_data_acquisition='UNKNOWN; not included in replay time')
    dump(ROOT / 'raw/RUN_COST.json', cost)
    subprocess.run([PY, '-B', str(ROOT / 'code/report.py')], check=True)
    state('DEMO_COMPLETE', 'All frozen rounds and actual-output fault checks completed', 'Independent review and matched loaded contact/force acquisition; see HANDOFF.md')
    print('DEMO_COMPLETE', cost['seconds'], flush=True)
if __name__ == '__main__':
    main()
