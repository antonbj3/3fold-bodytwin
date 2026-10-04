"""One-command replay, no downloads and no mutations outside this lane/data."""
from common import *
import shutil, subprocess, sys, time, shlex

def main():
    diskcheck()
    runid = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out = DATA / 'runs' / runid
    out.mkdir(parents=True)
    started = time.perf_counter()
    processes = []
    previous = out / 'previous_outputs'
    previous.mkdir()
    for p in (ROOT / 'raw').glob('*RESULTS.json'):
        shutil.copy2(p, previous / p.name)
    if (ROOT / 'CODE_LOCK.json').exists():
        for x in read(ROOT / 'CODE_LOCK.json')['files']:
            assert sha(x['path']) == x['sha256'], x['path']
    for file in ('run_r1.py', 'run_r2.py', 'run_r3.py', 'run_r4.py', 'report.py', 'verify.py'):
        argv = [sys.executable] + (['-s'] if file == 'report.py' else []) + [str(ROOT / 'code' / file)]
        if file in ('run_r2.py', 'run_r4.py'):
            argv = [str(DENT / 'tasks/heavy_run.sh'), '3', '0'] + argv
        t = time.perf_counter()
        log = out / (file + '.log')
        with log.open('w') as f:
            ret = subprocess.run(argv, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        row = dict(argv=argv, exit_code=ret.returncode, wall_seconds=time.perf_counter() - t, log=str(log))
        processes.append(row)
        with (ROOT / 'COMMANDS.md').open('a') as f:
            f.write('\n- `' + shlex.join(argv) + '` → ' + str(ret.returncode) + ', log ' + str(log) + '\n')
        print(file, 'exit', ret.returncode, flush=True)
        if ret.returncode:
            dump(out / 'FAILURE.json', row)
            state('REPLAY_FAILED', file + ' exit ' + str(ret.returncode), 'Inspect preserved log; do not modify frozen gates')
            raise SystemExit(ret.returncode)
    for p in sorted((ROOT / 'raw').glob('*RESULTS.json')):
        shutil.copy2(p, out / p.name)
    for name in ('results.json', 'GRAPH_FEEDBACK.json', 'VERIFICATION.json'):
        shutil.copy2(ROOT / name, out / name)
    receipt = dict(run_id=runid, wall_seconds=time.perf_counter() - started, processes=processes, data_bytes=diskcheck(), source_sha256=sha(ROOT / 'SOURCE_LOCK.json'), result_sha256=sha(ROOT / 'results.json'), artifact_hashes=[dict(path=str(p), sha256=sha(p)) for p in sorted(out.glob('*.json'))], threads_max=4, gpu=False)
    dump(out / 'RUN_RECEIPT.json', receipt)
    dump(ROOT / 'LATEST_RUN.json', receipt)
    state('READY_FOR_INDEPENDENT_REVIEW', 'One-command R1-R4 replay and producer verification PASS', 'Acquire same-Demo1 tissue/pose observations; no further same-information model variants can identify them')
    print('One-command replay PASS', receipt['wall_seconds'], 's', flush=True)
if __name__ == '__main__':
    main()
