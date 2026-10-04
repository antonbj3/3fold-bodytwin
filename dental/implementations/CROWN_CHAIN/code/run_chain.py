from common_chain import *
import subprocess, time, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--limit', type=int)
    ap.add_argument('--phase', choices=['r1', 'r2', 'r3', 'referents', 'report', 'all'], default='all')
    args = ap.parse_args()
    if args.phase in ('referents', 'all'):
        from referents import run
        run()
        if args.phase == 'referents':
            return
    if args.phase == 'r2':
        from exact_repair import run
        run()
        return
    if args.phase == 'r3':
        from round3 import run
        run()
        return
    if args.phase == 'report':
        from report_chain import run
        run()
        return
    for name in ['PREREG_R1.json', 'DECOMPOSITION_R1.json', 'FROZEN_INPUTS.json', 'FROZEN_PREDICTIONS_R1.json']:
        verify_freeze(ROOT / name)
    start = time.perf_counter()
    inp = read(ROOT / 'FROZEN_INPUTS.json')
    rows = []
    for (i, r) in enumerate(inp['cohort'][:args.limit]):
        dest = ROOT / 'raw/R1' / f"{r['key']}.json"
        if not (args.resume and dest.exists()):
            log = ROOT / 'raw' / f'worker_{i:02}.log'
            with log.open('w') as f:
                code = subprocess.run([sys.executable, str(ROOT / 'code/worker_r1.py'), str(i)], stdout=f, stderr=subprocess.STDOUT).returncode
            if code:
                state('R1_IMPLEMENTATION_FAILURE', dict(index=i, returncode=code, log=str(log)), 'Preserve failure and repair without changing prereg')
                raise RuntimeError(str(log))
        rows.append(read(dest))
        state('R1_RUNNING', dict(completed=len(rows), requested=len(inp['cohort'])), 'Complete unchanged denominator')
        print(rows[-1]['key'], rows[-1]['design_gate']['verdict'], flush=True)
    dump(ROOT / 'raw/R1_ALL.json', dict(rows=rows, requested=len(inp['cohort']), seconds=time.perf_counter() - start, peak_child_rss_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024))
    if args.phase == 'all' and args.limit is None:
        from exact_repair import run as run_r2
        from round3 import run as run_r3
        from report_chain import run as report
        run_r2()
        run_r3()
        report()
if __name__ == '__main__':
    main()
