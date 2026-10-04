from common import *
import subprocess, sys, shlex

def main():
    stages = [['source_gold.py'], ['sufficiency.py'], ['experiment.py', 'evaluate'], ['material_port.py'], ['absence.py', 'evaluate'], ['spacing.py'], ['evidence_summary.py'], ['replacement_correction.py'], ['measurement_port.py'], ['verify.py'], ['demo_case.py'], ['figure.py'], ['package.py']]
    runid = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S')
    for args in stages:
        cmd = [sys.executable, str(ROOT / 'code' / args[0]), *args[1:]]
        st = time.perf_counter()
        started = now()
        log = ROOT / 'raw' / ('REPLAY_' + runid + '_' + args[0].replace('.py', '') + '.log')
        with log.open('w') as f:
            r = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env=os.environ.copy())
        record = dict(started_utc=started, command=shlex.join(cmd), cwd=str(ROOT), exit_code=r.returncode, wall_s=time.perf_counter() - st, log=str(log), log_sha256=sha(log))
        with (ROOT / 'raw/COMMAND_LEDGER.jsonl').open('a') as f:
            f.write(json.dumps(record) + '\n')
        print(args[0], 'PASS' if r.returncode == 0 else 'FAILED', flush=True)
        if r.returncode:
            raise SystemExit('See ' + str(log))
    state('DEMO_VERIFIED', 'Executable/source checks PASS; scientific gates preserved', 'Independent review then fill frozen complete-negative native/site/material observation panel; no same-geometry material refit', active_background_jobs=0, review_state='PENDING_INDEPENDENT_REVIEW')
if __name__ == '__main__':
    main()
