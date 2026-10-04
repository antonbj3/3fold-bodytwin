"""Sequential, low-memory replay with full cost and hash receipt."""
import subprocess, sys, time, resource, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
R = Path(__file__).resolve().parent
start = time.perf_counter()
commands = []
for args in [['r1_identify.py'], ['science.py', '2'], ['science.py', '3'], ['science.py', '4'], ['source_controls.py'], ['validation.py'], ['freeze_predictions.py'], ['report.py']]:
    t = time.perf_counter()
    cmd = [sys.executable, '-s'] + args
    p = subprocess.run(cmd, cwd=R, check=False)
    commands.append({'argv': cmd, 'exit_code': p.returncode, 'wall_seconds': time.perf_counter() - t})
    if p.returncode:
        (R / 'raw/PIPELINE_FAILURE.json').write_text(json.dumps({'commands': commands, 'failed_utc': datetime.now(timezone.utc).isoformat()}, indent=2) + '\n')
        sys.exit(p.returncode)
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
receipt = {'completed_utc': datetime.now(timezone.utc).isoformat(), 'wall_seconds': time.perf_counter() - start, 'max_child_rss_kib': usage.ru_maxrss, 'child_user_cpu_seconds': usage.ru_utime, 'child_system_cpu_seconds': usage.ru_stime, 'commands': commands, 'threads': 1, 'gpu': False, 'questions': 0, 'network_calls_during_replay': 0, 'preparation_fit_discovery_included': 'Local input parsing, frozen logic and calibration sort, validation and figure generation included; prior manual source search and upstream data creation reported separately', 'manual_source_discovery_seconds': 'UNKNOWN'}
(R / 'raw/RUN_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
results = json.loads((R / 'results.json').read_text())
results['cost']['full_replay'] = receipt
(R / 'results.json').write_text(json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
feedback = json.loads((R / 'GRAPH_FEEDBACK.json').read_text())
feedback['sha256'] = hashlib.sha256((R / 'results.json').read_bytes()).hexdigest()
(R / 'GRAPH_FEEDBACK.json').write_text(json.dumps(feedback, ensure_ascii=False, indent=2) + '\n')
manifest = []
for p in sorted(R.glob('*.py')) + sorted(R.glob('*.sh')):
    manifest.append({'path': p.name, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'bytes': p.stat().st_size})
(R / 'CODE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
state = json.loads((R / 'CURRENT_WORK_STATE.json').read_text())
state['last_replay'] = receipt
state['graph_feedback'] = 'Prepared with current result hash; dispatch blocked by older Result hash mismatch; feedback blocked by missing dispatch'
(R / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, ensure_ascii=False, indent=2) + '\n')
print(json.dumps({'full_replay_seconds': receipt['wall_seconds'], 'max_rss_MiB': usage.ru_maxrss / 1024, 'result_sha256': feedback['sha256']}))
