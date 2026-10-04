"""One-command bounded replay. Original preregs and prediction freezes persist."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '4'
import time, resource, json, datetime
from pathlib import Path
P = Path(__file__).resolve().parent

def run():
    start = time.perf_counter()
    steps = []
    from extract_measurements import extract
    from predict_r1 import run as p1
    from score_r1 import run as s1
    from predict_r2 import run as p2
    from score_r2 import run as s2
    from predict_r3 import run as p3
    from score_r3 import run as s3
    from cement_port import run as cp
    from inverse_query import demo as iq
    from build_demo import reports
    for (name, f) in [('extract', extract), ('R1 predict', p1), ('R1 score', s1), ('R2 predict', p2), ('R2 score', s2), ('R3 predict', p3), ('R3 score', s3), ('physical port', cp), ('inverse queries', iq), ('report + figure', reports)]:
        t = time.perf_counter()
        f()
        steps.append(dict(operation=name, seconds=time.perf_counter() - t))
    account = dict(replay_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), wall_seconds=time.perf_counter() - start, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, cpu_user_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime, cpu_system_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_stime, threads=4, gpu=False, intermediate_limit_bytes=3000000000, own_files_bytes=sum((f.stat().st_size for f in P.rglob('*') if f.is_file())), steps=steps, historical_preparation_discovery_reasoning_cost='UNKNOWN', physical_acquisition_cost='UNKNOWN')
    from predict_r1 import dump
    dump(P / 'RUN_ACCOUNTING.json', account)
    from verify_demo import verify
    verify()
    dump(P / 'CURRENT_WORK_STATE.json', dict(lane='X13-cement-gap', phase='DEMO_READY', latest_gate='R1 FAIL; R2 FAIL; R3 LOCAL PILOT PASS/CLASSICAL TIE; R4 mean-field obstruction', original_full_regional_inverse='UNKNOWN', next_operation='Two-setting same-system calibration + held dose/batch and same-object spatial gap field; freeze future predictions before measurement', review_state='PENDING_INDEPENDENT_REVIEW', updated_utc=account['replay_utc'], result_file='results.json', graph_binding='BLOCKED_OLD_RESULT_HASH_MISMATCH'))
    print('DEMO_READY', round(account['wall_seconds'], 2), 's;', round(account['peak_rss_mib'], 1), 'MiB; full inverse UNKNOWN')
if __name__ == '__main__':
    run()
