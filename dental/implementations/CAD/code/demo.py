import thread_guard
from cadlib import *
import argparse, subprocess, resource

def verify_dependencies():
    if not (ROOT / 'DEPENDENCIES.json').exists():
        return
    for r in read(ROOT / 'DEPENDENCIES.json')['files']:
        if sha(r['path']) != r['sha256']:
            raise ValueError('Read-only dependency changed: ' + r['path'])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--recompute', action='store_true')
    ap.add_argument('--case')
    args = ap.parse_args()
    start = time.perf_counter()
    verify_dependencies()
    if args.case:
        dest = DATA / 'case_export'
        subprocess.run([sys.executable, str(ROOT / 'code/cad_tool.py'), '--manifest', args.case, '--output', str(dest)], check=True)
        return
    commands = [['code/controls.py'], ['code/reproduce_margin.py'] + (['--all'] if args.recompute else []), ['code/reproduce_geometry.py'] + (['--all'] if args.recompute else []), ['code/round5_score.py', 'score'], ['code/report.py']]
    for (i, cmd) in enumerate(commands):
        print('DEMO_STAGE', i + 1, '/', len(commands), cmd[0], flush=True)
        log = ROOT / 'raw' / ('DEMO_' + Path(cmd[0]).stem + '.log')
        with log.open('w') as f:
            subprocess.run([sys.executable, *cmd], cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, check=True)
        print('DEMO_STAGE_COMPLETE', cmd[0], flush=True)
    receipt = dict(command='./run_all.sh' + (' --recompute' if args.recompute else ''), finished_utc=now(), seconds=time.perf_counter() - start, full_regeneration=args.recompute, peak_child_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, results_sha256=sha(ROOT / 'results.json'), controls_sha256=sha(ROOT / 'raw/CONTROLS.json'), status='REPRODUCED_PARTIAL_CAPABILITY', physical_complete=False)
    dump(ROOT / 'raw/DEMO_RECEIPT_FULL.json' if args.recompute else ROOT / 'raw/DEMO_RECEIPT.json', receipt)
    state('DELIVERY_REPRODUCED', read(ROOT / 'results.json')['summary'], 'Full mirrored axial exterior with an independently annotated preparation margin; LAB_PROTOCOL.md supplies missing physical measurements')
    print(json.dumps(receipt), flush=True)
if __name__ == '__main__':
    main()
