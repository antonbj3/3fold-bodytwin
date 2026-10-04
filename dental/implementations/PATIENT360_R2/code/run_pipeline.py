import argparse, resource, subprocess, sys, time, traceback
from common import *
import pose, field_round, robust_design, contracts, certificate_check

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--package-existing', help='Regenerate reports from a completed own run, without recomputing')
    args = parser.parse_args()
    if args.package_existing:
        from package import publish_local
        publish_local(Path(args.package_existing))
        return
    started = time.perf_counter()
    cpu = time.process_time()
    budget()
    DATA.mkdir(parents=True, exist_ok=True)
    for p in ROOT.glob('PREREG_*.json'):
        if sha(p) != p.with_suffix('.sha256').read_text().strip():
            raise ValueError('PREREG_HASH_CHANGED:' + str(p))
    lock = ROOT / 'SOURCE_LOCK.json'
    if lock.exists():
        for row in load(lock)['files']:
            if sha(row['path']) != row['sha256']:
                raise ValueError('SOURCE_DRIFT:' + row['path'])
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    run = ROOT / 'runs' / stamp
    run.mkdir(parents=True)
    dump(run / 'CODE_LOCK.json', {str(p.relative_to(ROOT)): sha(p) for p in (ROOT / 'code').glob('*.py')})
    try:
        state('C1_RUNNING', 'Source/prereg hashes verified', 'Pose events, FE counterexample and preserved STL precision failure', run=str(run))
        c1 = pose.run(run / 'C1')
        print('C1:', c1['contract_status'], flush=True)
        state('C2_RUNNING', 'C1 support loss decided; derived STL precision failure preserved', 'Original source-addressed IOS/CBCT fields and registration connector', run=str(run))
        (c2, fields) = field_round.run(run / 'C2', DATA / stamp / 'fields')
        print('C2:', len(c2['patients']), 'patients', flush=True)
        state('C4_RUNNING', 'Source and crop-coordinate gates executed', 'Exact source-facet witness for robust rigid crown offset', run=str(run))
        c4 = robust_design.run(run / 'C4', fields)
        print('C4:', c4['status'], flush=True)
        state('C3_RUNNING', 'Exact robust-offset inequalities frozen', 'Per-decision X26 classes, shared-theta execution and export contracts', run=str(run))
        c3 = contracts.run(run / 'C3', c2, fields, c1, c4)
        print('C3:', c3['decision_count'], 'typed decisions', flush=True)
        independent = certificate_check.run(run / 'CHECK', c4, fields)
        expected_failed = {'original_triangle_ceiling'}
        c1_fail = {c['name'] for c in c1['checks'] if not c['pass_']}
        unexpected = [c for res in (c2, c3, c4, independent) for c in res['checks'] if not c['pass_']]
        verification = dict(delivery_pass=c1_fail == expected_failed and (not unexpected), expected_preserved_failure=sorted(c1_fail), unexpected_failures=unexpected, checks=[dict(round=name, **c) for (name, res) in [('C1', c1), ('C2', c2), ('C3', c3), ('C4', c4), ('CHECK', independent)] for c in res['checks']], review_state='PENDING_INDEPENDENT_REVIEW', physical_validation='NOT_RUN')
        dump(run / 'VERIFICATION.json', verification)
        external = [v for (k, v) in SOURCE_FILES.items() if not Path(k).is_relative_to(ROOT) and (not Path(k).is_relative_to(DATA))]
        for row in external:
            if sha(row['path']) != row['sha256']:
                raise ValueError('MID_RUN_SOURCE_DRIFT:' + row['path'])
        if not lock.exists():
            freeze(lock, dict(files=external))
        dump(run / 'SOURCE_FILES.json', list(SOURCE_FILES.values()))
        runtime = dict(wall_seconds=time.perf_counter() - started, cpu_seconds=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, GPU_seconds=0, threads_max=4, own_data_bytes=budget(), preparation='Included: raw source reads, fields and imports', fitting='No new fit; inherited X11/prior fit cost UNKNOWN', inherited_parent_execution_cost='UNKNOWN, excluded from current replay time; no speedup claim', validation='Included in round timings; no physical measurement', future_acquisition_cost='UNKNOWN', independent_review_cost='UNKNOWN')
        dump(run / 'RUNTIME.json', runtime)
        if not verification['delivery_pass']:
            raise RuntimeError('UNEXPECTED_VERIFICATION_FAILURE; artifacts preserved')
        from package import publish_local
        publish_local(run)
        print('DEMO COMPLETE:', run, flush=True)
    except Exception as exc:
        dump(run / 'FAILURE.json', dict(error=str(exc), traceback=traceback.format_exc(), time=now()))
        state('FAILED_CHECKPOINT', str(exc), 'Inspect frozen failure, repair operation without changing gates', run=str(run))
        raise
if __name__ == '__main__':
    main()
