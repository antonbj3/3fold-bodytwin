"""One-command rebuild, recording negative scientific gates without hiding them."""
import datetime
import os
import resource
import subprocess
import sys
import time
from common import *

def main():
    started = time.perf_counter()
    records = []
    env = os.environ.copy()
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[key] = '1'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['MPLBACKEND'] = 'Agg'
    (HERE / 'raw').mkdir(exist_ok=True)
    commands = [['python3', 'prepare_inputs.py'], ['python3', 'build_net.py'], ['python3', 'validate.py', '--mutations'], ['python3', 'stress_and_queries.py'], ['python3', 'verify_queries.py'], ['python3', 'residual_query.py'], ['/usr/bin/python3', 'make_figure.py']]
    for (i, cmd) in enumerate(commands):
        step_env = dict(env)
        if 'make_figure.py' in cmd:
            step_env['PYTHONNOUSERSITE'] = '1'
        t = time.perf_counter()
        p = subprocess.run(cmd, cwd=HERE, env=step_env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rec = {'argv': cmd, 'wall_s': time.perf_counter() - t, 'exit_code': p.returncode, 'log': f'raw/run_{i:02}.log', 'environment_threads': 1}
        (HERE / rec['log']).write_text(p.stdout)
        records.append(rec)
        print(' '.join(cmd) + f' -> {p.returncode}', flush=True)
        if p.returncode:
            if 'validate.py' not in cmd:
                dump('FAILED_EXECUTION.json', records)
                print(p.stdout)
                return p.returncode
            v = read(HERE / 'VALIDATION.json')
            if not v['integrity_pass'] or not v['fault_controls']['all_rejected']:
                dump('FAILED_EXECUTION.json', records)
                print(p.stdout)
                return 2
        state = read(HERE / 'CURRENT_WORK_STATE.json')
        state.update(status='REPLAY_RUNNING', latest_replay_step=rec, next_operation=commands[i + 1] if i + 1 < len(commands) else 'report')
        dump('CURRENT_WORK_STATE.json', state)
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    dump('COST_RUN.json', {'measured_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'preparation': 'Source reading and research elapsed not instrumented; not zero', 'fit_s': 0, 'discovery': 'Two frozen exact arithmetic/representation constructions; reasoning cost UNKNOWN', 'validation_and_queries_wall_s': time.perf_counter() - started, 'child_cpu_user_s': usage.ru_utime, 'child_cpu_system_s': usage.ru_stime, 'peak_child_RSS_MiB': usage.ru_maxrss / 1024, 'gpu_hours': 0, 'physical_measurements': 0, 'threads_per_process': 1, 'fallback': 'Strict completeness FAIL retained; no substitute physical values', 'steps': records, 'packaging_time_excluded': True, 'historical_source_preparation_and_fit_cost': 'UNKNOWN, not claimed free'})
    p = subprocess.run(['python3', 'report.py'], cwd=HERE, env=env)
    if p.returncode:
        return p.returncode
    dump('ARTIFACT_MANIFEST.json', {p.name: {'sha256': sha(p), 'bytes': p.stat().st_size} for p in HERE.iterdir() if p.is_file() and p.name not in ['ARTIFACT_MANIFEST.json', 'GRAPH_FEEDBACK_RECEIPT.json', 'CURRENT_WORK_STATE.json']})
    v = read(HERE / 'VALIDATION.json')
    print(f"Saved all artifacts. Strict numeric completeness: {v['numeric_open_unknown_edges']}/{v['numeric_open_unknown_edges'] + v['missing_numeric_open_unknown_edges']} — FAIL.", flush=True)
    return 0 if v['pass'] else 1
if __name__ == '__main__':
    sys.exit(main())
