"""Record exact argv, complete stderr, cost and failure without shell interpretation."""
from common import *
import sys, subprocess, time, resource

def main():
    argv = sys.argv[1:]
    st = time.perf_counter()
    started = now()
    log = P / 'raw/COMMAND_LOG.jsonl'
    log.parent.mkdir(exist_ok=True)
    cp = subprocess.run(argv, cwd=P, text=True, capture_output=True)
    entry = dict(started_utc=started, argv=argv, exit_code=cp.returncode, wall_s=time.perf_counter() - st, child_peak_rss_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, stdout=cp.stdout, stderr=cp.stderr)
    with log.open('a') as f:
        f.write(json.dumps(entry) + '\n')
    print(cp.stdout, end='')
    print(cp.stderr, end='', file=sys.stderr)
    sys.exit(cp.returncode)
if __name__ == '__main__':
    main()
