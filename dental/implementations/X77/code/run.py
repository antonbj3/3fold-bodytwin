"""Read-only package composition; each run is preserved, no mutable scientific cache."""
from dental_release.paths import expand as _release_expand
import argparse, datetime, hashlib, json, os, resource, shlex, subprocess, time
from pathlib import Path
R = Path(__file__).resolve().parents[1]
D = R.parent.parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X77_same_patient_flow'))
LOGICAL = DATA / 'work'
PY = _release_expand('@DENTAL_PYTHON@')

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def write(p, x):
    Path(p).write_text(json.dumps(x, indent=2, sort_keys=True, ensure_ascii=False) + '\n')

def state(phase, gate, nextop):
    write(R / 'CURRENT_WORK_STATE.json', dict(lane='X77-same-patient-flow', claim_type='capability', phase=phase, latest_gate=gate, next_operation=nextop, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), review_state='PENDING_INDEPENDENT_REVIEW'))

def check():
    for name in ['PATIENT_SELECTION.json', 'PACKAGE_MAP.json', 'PREREG_R1_SAME_PATIENT.json', 'PREREG_R2_CORRELATED_DECISIONS.json', 'PREREG_R1B_COLLISION_INSTRUMENT.json', 'PREREG_R3_FULL_PRESCAN.json', 'PREREG_R4_LOCAL_COMPLETION.json', 'PREREG_R5_CANONICAL_RECEIPT.json', 'DECOMPOSITION.json', 'DECOMPOSITION_R3_R4.json', 'SOURCE_LOCK.json']:
        if sha(R / name) != (R / (name + '.sha256')).read_text().strip():
            raise ValueError('Frozen input drift ' + name)
    for q in json.loads((R / 'SOURCE_LOCK.json').read_text())['files']:
        if sha(q['path']) != q['sha256']:
            raise ValueError('Source drift ' + q['path'])
    for q in json.loads((R / 'PREREG_R3_FULL_PRESCAN.json').read_text())['source_lock']:
        if sha(q['path']) != q['sha256']:
            raise ValueError('R3 source drift ' + q['path'])

def budget():
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if n > 3000000000:
        raise RuntimeError('Own intermediates exceed 3GB')
    return n

def command(stage, out):
    cmd = ['/usr/bin/bwrap', '--die-with-parent', '--unshare-net', '--ro-bind', '/', '/', '--bind', str(out), str(LOGICAL), '--tmpfs', '/tmp', '--proc', '/proc', '--dev', '/dev', '--chdir', str(R)]
    for (logical, pkg) in json.loads((R / 'PACKAGE_MAP.json').read_text())['mounts'].items():
        cmd += ['--ro-bind', pkg, logical]
    py = PY if not stage.startswith('full_crown') else _release_expand('@DENTAL_PYTHON@')
    worker = R / 'code/worker.py' if not stage.startswith('full_crown') else R / 'code/full_crown.py'
    cmd += ['--setenv', 'X60_RUN_ROOT', str(LOGICAL), '--setenv', 'X60_RUN_DATA', str(LOGICAL), '--setenv', 'PYTHONDONTWRITEBYTECODE', '1', py, '-B', str(worker), stage, str(LOGICAL)]
    return cmd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--single-run', action='store_true')
    a = ap.parse_args()
    full_start = time.perf_counter()
    check()
    preflight_s = time.perf_counter() - full_start
    DATA.mkdir(parents=True, exist_ok=True)
    LOGICAL.mkdir(exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    root = DATA / 'runs' / stamp
    root.mkdir(parents=True)
    st = time.perf_counter()
    runs = []
    try:
        for tag in ['A'] if a.single_run else ['A', 'B']:
            out = root / tag
            out.mkdir()
            (out / 'raw').mkdir()
            rows = []
            for stage in ['source', 'sdf', 'generate', 'gate', 'milling', 'physical', 'uncertainty', 'full_crown', 'full_crown_R4']:
                budget()
                state('RUN_' + tag + '_' + stage.upper(), 'RUNNING', 'Execute ' + stage)
                cmd = command(stage, out)
                tic = time.perf_counter()
                with (out / 'raw' / (stage + '.log')).open('w') as log:
                    p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
                row = dict(stage=stage, argv=cmd, exit_code=p.returncode, wall_seconds=time.perf_counter() - tic)
                rows.append(row)
                write(out / 'PROCESS_RECEIPT.json', rows)
                with (R / 'COMMANDS.md').open('a') as f:
                    f.write('\n`' + shlex.join(cmd) + '`\n\nExit ' + str(p.returncode) + '; log ' + str(out / 'raw' / (stage + '.log')) + '\n')
                print(tag, stage, 'exit', p.returncode, 's', round(row['wall_seconds'], 2), flush=True)
                if p.returncode:
                    raise RuntimeError('Preserved stage failure ' + str(out / 'raw' / (stage + '.log')))
            runs.append(dict(tag=tag, path=str(out), processes=rows))
            write(root / ('COMPLETED_' + tag + '.json'), dict(path=str(out), registration=json.loads((out / 'REGISTRATION.json').read_text()), anatomical_status='UNKNOWN missing same-patient tissues and independent bounds'))
        check()
        names = ['TASK.json', 'TASK_ARRAYS.npz', 'SCAN.json', 'scene.npz', 'segmented_tooth.stl', 'antagonist.stl', 'measured_capped_roof.stl', 'REGISTRATION.json', 'REGISTRATION_PER_REGION.csv', 'prep.stl', 'SDF.npz', 'SDF.json', 'crown.stl', 'GENERATION.json', 'DESIGN_FIELDS.npz', 'CONTRACT.json', 'REGIONS.json', 'GENERATOR_FAULTS.json', 'GATE.json', 'GATE_FAULTS.json', 'EXPORT.json', 'EXPORT_FRAME.json', 'MILLING.json', 'PHYSICAL.json', 'UNCERTAINTY.json', 'ANATOMICAL_DECISIONS.json', 'export/model.3mf', 'export/model.3mf.json', 'export/DESIGN_GATE.json']
        names += ['GENERATOR_FAULT_REPAIR.json', 'FULL_CROWN.json', 'FULL_CROWN_R4.json', 'DECISION_ERROR_MAP.json']
        names += [str(p.relative_to(root / 'A')) for p in sorted((root / 'A' / 'full_crown').rglob('*')) if p.is_file() and p.suffix in ('.stl', '.npz', '.3mf', '.json')]
        names += [str(p.relative_to(root / 'A')) for p in sorted((root / 'A' / 'full_crown_R4').rglob('*')) if p.is_file() and p.suffix in ('.stl', '.npz', '.3mf', '.json')]
        files = {n: dict(sha256=sha(root / 'A' / n), bytes=(root / 'A' / n).stat().st_size, second_sha256=sha(root / 'B' / n) if not a.single_run else None) for n in names}
        different = [n for (n, q) in files.items() if q['second_sha256'] is not None and q['second_sha256'] != q['sha256']]
        receipt = dict(run_id=stamp, runs=runs, scientific_artifacts=files, byte_identity_gate='NOT_RUN_SINGLE' if a.single_run else 'PASS' if not different else 'FAIL', different_files=different, wall_seconds=time.perf_counter() - st, preflight_hash_seconds=preflight_s, complete_workflow_wall_seconds=time.perf_counter() - full_start, max_child_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, data_bytes=budget(), runtime_telemetry_excluded=['PROCESS_RECEIPT.json', 'raw/*.log', 'GATE_TELEMETRY.json'], code_sha256={str(p.relative_to(R)): sha(p) for p in sorted((R / 'code').glob('*.py'))})
        write(root / 'RECEIPT.json', receipt)
        write(R / 'LATEST_RUN.json', receipt)
        if different:
            raise RuntimeError('Scientific byte drift: ' + str(different))
        subprocess.run([PY, '-B', str(R / 'code/report.py')], check=True)
        state('DIGITAL_CHAIN_AND_CONDITIONAL_DECISIONS_EXECUTED', 'Partial same-patient capability; physical surgery UNKNOWN', 'Same-patient CBCT semantic tissue and independent directed bounds; full axial specimen construction')
    except Exception as e:
        write(root / 'FAILURE.json', dict(error=str(e), runs=runs, wall_seconds=time.perf_counter() - st))
        state('FAILED', str(e), 'Preserve failure and change the diagnosed operation')
        raise
if __name__ == '__main__':
    main()
