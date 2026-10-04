from common import *
import time, resource

def run():
    start = time.perf_counter()
    check_lock()
    from acquire import run as acquire
    from construct import run as construct
    from validate import run as validate
    from reachability import run as reachability
    from local_boundary import run as local_boundary
    from local_uniform import run as local_uniform
    from report import run as report
    acquire()
    construct()
    validate()
    reachability()
    local_boundary()
    local_uniform()
    result = report()
    receipt = dict(started_utc=datetime.datetime.fromtimestamp(time.time() - (time.perf_counter() - start), datetime.timezone.utc).isoformat(), finished_utc=now(), wall_seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads=1, physical_measurement='NOT_RUN', R1_gate=result['decision']['inherited_50pct_3of10_gate'], R2_possible_cases=result['decision']['R2_all_family_possible_cases'], data_root=str(DATA), results_sha256=sha(ROOT / 'results.json'))
    dump(ROOT / 'raw/DEMO_RECEIPT.json', receipt)
    state('DEMO_COMPLETE_PENDING_REVIEW', result['decision'], 'Independent copied replay and local source/solid/physical calibration prerequisites in HANDOFF.md')
    print(json.dumps(receipt), flush=True)
if __name__ == '__main__':
    run()
