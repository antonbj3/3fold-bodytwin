from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, sys, json, subprocess, datetime, shutil, time, hashlib, resource
SRC = Path(__file__).resolve().parents[1]
base = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/X95-native-support-shape'))
if len(sys.argv) > 1:
    target = Path(sys.argv[1]).resolve()
else:
    target = base / ('REPLAY_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%f'))
if target.exists():
    raise ValueError('Replay destination must be new: ' + str(target))
target.mkdir(parents=True)
shutil.copytree(SRC / 'code', target / 'code', ignore=shutil.ignore_patterns('__pycache__'))
for pattern in ['INPUT_LOCK.json*', 'DEPENDENCY_LOCK.json*', 'SOURCE_NATIVE_REFERENCE_LOCK.json*', 'PREREG*.json*', 'FROZEN_PREDICTIONS.json*']:
    for p in SRC.glob(pattern):
        shutil.copy2(p, target / p.name)
env = dict(os.environ, X95_OUTPUT_ROOT=str(target), X95_DATA_ROOT=str(target / 'data'), PYTHONDONTWRITEBYTECODE='1')
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    env[k] = '4'
interpreter = _release_expand('@DENTAL_PYTHON@')
steps = ['dependency_lock', 'native_port', 'support_coverage', 'source_solver', 'compiled_support', 'port_certificate', 'certificate_decimal_check', 'source_facit', 'verify', 'report']
t0 = time.monotonic()
receipts = []
for step in steps:
    cmd = [interpreter, str(target / 'code' / f'{step}.py')]
    if step == 'source_solver':
        gate_cmd = ['python3', _release_expand('@DENTAL_EXTERNAL_ROOT@/projects/bodytwin/scripts/resource_gate.py'), 'check', '--ram-gb', '2', '--vram-gb', '0', '--max-load-per-core', '0']
        gate = subprocess.run(gate_cmd, capture_output=True, text=True)
        (target / 'raw/RESOURCE_PREFLIGHT.txt').write_text(gate.stdout + gate.stderr)
        if gate.returncode:
            cmd += ['--cached-reference']
        else:
            cmd = [_release_expand('@DENTAL_INPUT_ROOT@/workspace/tasks/heavy_run.sh'), '2', '0'] + cmd
    log = target / 'raw' / f'{step}_REPLAY.log'
    log.parent.mkdir(parents=True, exist_ok=True)
    ts = time.monotonic()
    with log.open('w') as f:
        r = subprocess.run(cmd, env=env, cwd=target, stdout=f, stderr=subprocess.STDOUT)
    receipts.append({'step': step, 'command': cmd, 'returncode': r.returncode, 'seconds_including_wait': time.monotonic() - ts})
    if r.returncode:
        (target / 'REPLAY_FAILURE.json').write_text(json.dumps(receipts, indent=2) + '\n')
        print(log.read_text()[-2000:])
        raise SystemExit(r.returncode)
expected = json.loads((SRC / 'results.json').read_text())
got = json.loads((target / 'results.json').read_text())
comparisons = {}

def add(name, a, b, tol=0):
    if isinstance(a, (int, float)):
        err = abs(a - b)
        ok = err <= tol
    else:
        err = 0 if a == b else None
        ok = a == b
    comparisons[name] = {'original': a, 'replayed': b, 'absolute_error': err, 'tolerance': tol, 'pass': ok}
add('outcome', expected['outcome'], got['outcome'])
add('native_max_projection', expected['rounds']['R1']['interface_projection_mm']['max'], got['rounds']['R1']['interface_projection_mm']['max'])
add('qualified_count', expected['rounds']['R2']['native_coverage']['qualified_nodes'], got['rounds']['R2']['native_coverage']['qualified_nodes'])
add('omitted_count', expected['rounds']['R3']['counts']['compiled_omitted_nodes'], got['rounds']['R3']['counts']['compiled_omitted_nodes'])
add('sufficiency_identity', expected['sufficiency']['identity_error'], got['sufficiency']['identity_error'])
add('sufficiency_energy_difference', expected['sufficiency']['relative_downstream_difference'], got['sufficiency']['relative_downstream_difference'], 1e-10)
for name in ['native_all', 'all_master_all', 'coverage_qualified']:
    for i in range(2):
        add(f'energy_{name}_{i}', expected['rounds']['R2']['contrasts_P1_conditional'][name][i]['strain_energy_N_mm_per_unit_load'], got['rounds']['R2']['contrasts_P1_conditional'][name][i]['strain_energy_N_mm_per_unit_load'], 1e-13)
for i in range(2):
    add(f'certified_bound_{i}', expected['rigorous_enclosure']['geometry'][i], got['rigorous_enclosure']['geometry'][i])
add('source_ratio', expected['source_protocol_contrast']['ratio_3Y_5Y'], got['source_protocol_contrast']['ratio_3Y_5Y'])
add('validation_passed', expected['validation']['passed'], got['validation']['passed'])
add('shape_optimizer_runs', 0, got['shape_optimizer_runs'])
receipt = {'status': 'REPLAY_PASS' if all((v['pass'] for v in comparisons.values())) else 'REPLAY_FAIL', 'target': str(target), 'steps': receipts, 'comparisons': comparisons, 'total_wall_seconds_including_gate_wait': time.monotonic() - t0, 'scope': 'Fresh local geometry/FE/KKT/enclosure/source-table/fault replay. Native step is fresh only if gate allows; otherwise hash-verified actual original native reference is reused, explicitly labeled. No clean OS installation, physical measurement or independent scientific review', 'native_replay_status': json.loads((target / 'raw/NATIVE_SOLVER.json').read_text()).get('replay_native_status', 'FRESH_NATIVE_SOLVE'), 'code_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((target / 'code').glob('*.py'))}}
(target / 'REPLAY_VERIFICATION.json').write_text(json.dumps(receipt, indent=2) + '\n')
(SRC / 'raw/REPLAY_VERIFICATION.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'status': receipt['status'], 'target': str(target), 'wall_seconds': receipt['total_wall_seconds_including_gate_wait']}))
if receipt['status'] != 'REPLAY_PASS':
    raise SystemExit(1)
