"""One-command offline benchmark: controls, local build, isolated generation, freeze, separate score."""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'THREADS']:
    os.environ[key] = '4'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys, time, argparse, subprocess
from pathlib import Path
from gencad_bench_v2.common import *

def main():
    if sys.version_info[:3] != (3, 10, 12):
        raise RuntimeError('Locked runtime CPython 3.10.12 required')
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=['all', 'controls', 'build', 'generate', 'score', 'report'], default='all')
    ap.add_argument('--recompute', action='store_true', help='Rerun all isolated participants to separate files and require identical frozen designs')
    a = ap.parse_args()
    start = time.perf_counter()
    phase = a.phase
    if phase in ['all', 'controls']:
        from gencad_bench_v2.audit import run
        state('TESTING_CONTROLS', 'No numerical result admitted yet', 'Positive fixtures and injected faults')
        run()
    if phase in ['all', 'build']:
        from gencad_bench_v2.build import run
        run()
    if phase in ['all', 'generate']:
        from gencad_bench_v2.measurements import predict
        from gencad_bench_v2.sandbox import run
        from gencad_bench_v2.generators import PARTICIPANTS
        lock = ROOT / 'data/IMPLEMENTATION_LOCK_R1.json'
        if lock.exists():
            for (rel, h) in read(lock)['payload'].items():
                if sha(ROOT / rel) != h:
                    raise ValueError('R1 implementation drift: ' + rel)
        else:
            code = {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / 'gencad_bench_v2').rglob('*.py')) if p.name not in ['report.py', 'plot.py']}
            freeze(lock, code)
        participants = list(PARTICIPANTS)
        dirs = {}
        predfiles = {}
        for name in participants:
            dest = DATA / 'predictions' / name
            dirs[name] = str(dest.relative_to(DATA))
            if not (ROOT / 'FROZEN_PREDICTIONS_R1.json').exists():
                state('GENERATING', 'Public-only sandbox active', 'Participant ' + name)
                run(name, output=dest)
            for p in sorted(dest.glob('*.json')):
                if p.name != '_cost.json':
                    predfiles[str(p.relative_to(DATA))] = sha(p)
        freeze(ROOT / 'FROZEN_PREDICTIONS_R1.json', dict(participants=participants, prediction_dirs=dirs, files=predfiles, implementation_sha256=sha(ROOT / 'data/IMPLEMENTATION_LOCK_R1.json'), benchmark_sha256=sha(ROOT / 'data/BENCHMARK_LOCK.json'), reference_read_by_generators=False))
        if a.recompute:
            recomputed = []
            for name in participants:
                dest = DATA / 'recomputed_predictions' / name
                state('RECOMPUTING_PARTICIPANTS', 'Frozen originals retained', 'Same public sandbox, participant ' + name)
                run(name, output=dest)
                bad = []
                expected = DATA / 'predictions' / name
                for p in sorted(expected.glob('*.json')):
                    if p.name == '_cost.json':
                        continue
                    q = dest / p.name
                    if not q.exists() or sha(p) != sha(q):
                        bad.append(p.name)
                recomputed.append(dict(participant=name, identical=not bad, differences=bad))
            dump(ROOT / 'raw/RECOMPUTATION_RECEIPT.json', dict(participants=recomputed, all_identical=all((r['identical'] for r in recomputed))))
            if not all((r['identical'] for r in recomputed)):
                raise ValueError('Recomputed designs differ; originals and differences preserved')
        freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(round='R1', manifest='FROZEN_PREDICTIONS_R1.json', sha256=sha(ROOT / 'FROZEN_PREDICTIONS_R1.json')))
        predict()
        state('PREDICTIONS_FROZEN', 'All task predictions hashed before scorer', 'Separate scorer opens withheld references')
    if phase in ['all', 'score']:
        cmd = [sys.executable, '-m', 'gencad_bench_v2.score_cli']
        p = subprocess.run(cmd, cwd=ROOT)
        if p.returncode:
            raise RuntimeError('separate scorer failed')
        for (module, args) in [('axial_round', ['all']), ('axial_round3', ['all']), ('witness_replay', []), ('adversarial_integration', [])]:
            p = subprocess.run([sys.executable, '-m', 'gencad_bench_v2.' + module] + args, cwd=ROOT)
            if p.returncode:
                raise RuntimeError('verification/extension failed: ' + module)
    if phase in ['all', 'report']:
        from gencad_bench_v2.report import run
        run()
    dump(ROOT / 'raw/LAST_COMMAND.json', dict(argv=sys.argv, seconds=time.perf_counter() - start, data_bytes=budget(), python=sys.version, threads=4))
if __name__ == '__main__':
    main()
