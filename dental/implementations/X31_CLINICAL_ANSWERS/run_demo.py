"""One-command runner with actual stage timings and retained failed logs."""
import datetime, hashlib, json, resource, subprocess, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    logdir = ROOT / 'logs' / stamp
    logdir.mkdir(parents=True)
    before = {str(p.relative_to(ROOT)): sha(p) for p in sorted(ROOT.glob('*.py'))}
    before['run_all.sh'] = sha(ROOT / 'run_all.sh')
    (logdir / 'CODE_HASHES.json').write_text(json.dumps(before, indent=2) + '\n')
    for name in before:
        (logdir / name).write_bytes((ROOT / name).read_bytes())
    started = time.perf_counter()
    records = []
    for command in [[sys.executable, 'analysis.py'], [sys.executable, 'measurements.py', '--self-check'], [sys.executable, 'minimum_observation.py']]:
        t = time.perf_counter()
        p = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        duration = time.perf_counter() - t
        (logdir / (Path(command[1]).stem + '.stdout')).write_text(p.stdout)
        (logdir / (Path(command[1]).stem + '.stderr')).write_text(p.stderr)
        print(p.stdout, end='')
        print(p.stderr, end='', file=sys.stderr)
        records.append({'argv': command, 'wall_seconds': duration, 'exit_code': p.returncode})
        if p.returncode:
            (logdir / 'FAILED_RUN.json').write_text(json.dumps({'stage_records': records, 'frozen_inputs_preserved': True}, indent=2) + '\n')
            return p.returncode
    r = json.loads((ROOT / 'results.json').read_text())
    r['rounds']['R4'] = json.loads((ROOT / 'raw/R4_EVENT_RESULTS.json').read_text())
    used = resource.getrusage(resource.RUSAGE_CHILDREN)
    r['cost']['executed_stage_records'] = records
    r['cost']['analysis_report_guard_wall_seconds'] = time.perf_counter() - started
    r['cost']['children_cpu_seconds'] = used.ru_utime + used.ru_stime
    r['cost']['children_maxrss_MiB'] = used.ru_maxrss / 1024
    r['cost']['final_manifest_serialization_seconds'] = 'Stored in run receipt after finalization; excluded from prior-stage sum'
    (ROOT / 'results.json').write_text(json.dumps(r, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    t = time.perf_counter()
    p = subprocess.run([sys.executable, 'finalize.py'], cwd=ROOT, capture_output=True, text=True)
    print(p.stdout, end='')
    print(p.stderr, end='', file=sys.stderr)
    records.append({'argv': [sys.executable, 'finalize.py'], 'wall_seconds': time.perf_counter() - t, 'exit_code': p.returncode})
    (logdir / 'finalize.stdout').write_text(p.stdout)
    (logdir / 'finalize.stderr').write_text(p.stderr)
    receipt = {'run_utc': stamp, 'stage_records': records, 'full_command_wall_seconds': time.perf_counter() - started, 'code_unchanged': all((sha(ROOT / name) == h for (name, h) in before.items())), 'physical_measurements_n': 0, 'results_sha256': sha(ROOT / 'results.json')}
    (logdir / 'RUN_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (ROOT / 'LAST_RUN_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return p.returncode
if __name__ == '__main__':
    sys.exit(main())
