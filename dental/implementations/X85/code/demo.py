from common import *
import time, subprocess, resource

def main():
    from review_integrity import verify_bundle
    verify_bundle(ROOT)
    start = time.perf_counter()
    events = []
    for (p, d) in read(ROOT / 'DEPENDENCY_LOCK.json')['files'].items():
        if sha(p) != d['sha256']:
            raise ValueError('Read-only helper drift: ' + p)
    if (ROOT / 'SOURCE_STATE_LOCK.json').exists():
        for (p, d) in read(ROOT / 'SOURCE_STATE_LOCK.json')['files'].items():
            if sha(p) != d['sha256']:
                raise ValueError('X82 source state drift: ' + p)
    for script in ['round01.py', 'round2.py', 'round3.py', 'round4.py', 'round5.py', 'round6.py', 'round7.py', 'external.py']:
        tic = time.perf_counter()
        log = ROOT / 'raw' / ('replay_' + script.replace('.py', '.log'))
        with log.open('w') as f:
            p = subprocess.run([sys.executable, '-B', str(ROOT / 'code' / script)], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        events.append(dict(script=script, seconds=time.perf_counter() - tic, exit_code=p.returncode, log=str(log)))
        print(script, p.returncode, round(events[-1]['seconds'], 2), flush=True)
        if p.returncode:
            raise RuntimeError('Preserved failure ' + str(log))
    tic = time.perf_counter()
    with (ROOT / 'raw/TESTS.log').open('w') as f:
        p = subprocess.run([sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests', '-v'], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
    dump(ROOT / 'raw/TEST_RECEIPT.json', dict(exit_code=p.returncode, seconds=time.perf_counter() - tic, log_sha256=sha(ROOT / 'raw/TESTS.log'), count=17))
    if p.returncode:
        raise RuntimeError('Contract tests failed; raw/TESTS.log preserved')
    dump(ROOT / 'raw/DEMO_RECEIPT.json', dict(status='REPLAY_AND_CONTRACTS_PASS_SCIENTIFIC_UNKNOWNS_RETAINED', events=events, elapsed_seconds=time.perf_counter() - start, peak_children_rss_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, finished_utc=now(), python=sys.version, executable=sys.executable, physical_measurement='NOT_RUN'))
    from report import run
    run()
    print('Demo completed; results.json and figures/demo.png', flush=True)
if __name__ == '__main__':
    main()
