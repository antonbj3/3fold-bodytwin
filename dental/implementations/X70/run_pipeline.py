"""Two fresh runs, read-only package calls and exact scientific artifact comparison."""
from dental_release.paths import expand as _release_expand
import argparse, datetime, hashlib, json, os, resource, shutil, subprocess, time
from pathlib import Path
R = Path(__file__).resolve().parent
D = R.parent.parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X70_end_to_end'))
LOGICAL = DATA / 'work'
PYTHON = _release_expand('@DENTAL_PYTHON@')
STAGES = ['scan', 'sdf', 'generate', 'gate', 'milling', 'physical', 'metrology', 'lab']

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def write(p, x):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(x, indent=2, sort_keys=True, ensure_ascii=False) + '\n')

def state(phase, gate, next_operation):
    write(R / 'CURRENT_WORK_STATE.json', dict(lane='X70-end-to-end', claim_type='capability', phase=phase, latest_gate=gate, next_operation=next_operation, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), review_state='PENDING_INDEPENDENT_REVIEW'))

def check_sources():
    for name in ['SOURCE_LOCK.json', 'PATIENT_SELECTION.json', 'PREREG_X70_R1.json', 'PREREG_X70_R2.json', 'PREREG_RUNTIME_SEPARATION.json', 'PREREG_FAULT_REPAIR.json', 'DECOMPOSITION.json', 'PACKAGE_MAP.json']:
        p = R / name
        if sha(p) != (R / (name + '.sha256')).read_text().strip():
            raise ValueError('Frozen contract changed: ' + name)
    for row in json.loads((R / 'SOURCE_LOCK.json').read_text())['files']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Read-only package/input changed: ' + row['path'])

def command(stage, work):
    mounts = json.loads((R / 'PACKAGE_MAP.json').read_text())['mounts']
    cmd = ['/usr/bin/bwrap', '--die-with-parent', '--unshare-net', '--ro-bind', '/', '/', '--bind', str(R), str(R), '--bind', str(work), str(LOGICAL), '--tmpfs', '/tmp', '--proc', '/proc', '--dev', '/dev', '--chdir', str(R)]
    for (logical, template) in mounts.items():
        cmd += ['--ro-bind', template, logical]
    cmd += ['--setenv', 'X60_RUN_ROOT', str(LOGICAL), '--setenv', 'X60_RUN_DATA', str(LOGICAL), '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
    cmd += [PYTHON, '-B', str(R / 'worker.py'), stage, str(LOGICAL)]
    return cmd

def budget():
    n = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    if n > 3000000000:
        raise RuntimeError('Own data exceeds frozen 3 GB')
    return n

def frozen_predictions(work):
    gen = json.loads((work / 'GENERATION.json').read_text())
    phy = json.loads((work / 'PHYSICAL.json').read_text())
    science_files = ['crown.stl', 'prep.stl', 'antagonist.stl', 'DESIGN_FIELDS.npz', 'REGIONS.json', 'CONTRACT.json', 'export/model.3mf', 'export/model.3mf.json', 'export/DESIGN_GATE.json']
    payload = dict(claim_type='capability', patient_case='27', case_key='14ceebb036f59e99', frame=gen['frame'], files={s: dict(sha256=sha(work / s), bytes=(work / s).stat().st_size) for s in science_files}, conditional_model_force_ratio=gen['diagnostics']['model_force_quantile_ratio'], model_ratio_status='Unvalidated rigid-support X60 closure, not measured absolute force', calibrated_patient_force05_N=None, calibrated_patient_fracture_mean_N=None, source_support_status=phy['calibrated_patient_force']['status'], nominal_spatial_film_port='DESIGN_FIELDS.npz inner_z_mm - virtual_prep_z_mm', seated_spatial_cement_film_um=None, physical_measurement_status='NOT_RUN', diagnostic_digital_metrology_status='Synthetic controls already executed; no prospective physical data observed', full_crown_status='UNKNOWN_AXIAL_ANATOMY_ABSENT', missing_affine_enclosure='RIGOROUS_SENSITIVITY_AND_MESH_CONVERGENCE_ENCLOSURES_MISSING', freeze_rule='Before any physical lab outcome, metrology scan or batch calibration. Historical source group comparisons are retrospective; no source-trained patient likelihood transferred.')
    dest = R / 'FROZEN_PREDICTIONS.json'
    if dest.exists():
        old = json.loads(dest.read_text())
        if old['payload'] != payload or sha(dest) != (R / 'FROZEN_PREDICTIONS.json.sha256').read_text().strip():
            raise ValueError('Frozen patient prediction drift')
    else:
        write(dest, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), payload=payload))
        (R / 'FROZEN_PREDICTIONS.json.sha256').write_text(sha(dest) + '\n')
    shutil.copyfile(dest, work / 'FROZEN_PATIENT_PREDICTIONS.json')
    shutil.copyfile(dest, work / 'export/X70_FROZEN_PATIENT_PREDICTIONS.json')
    write(work / 'export/EXPORT_PROVENANCE.json', dict(patient_predictions=dict(path='X70_FROZEN_PATIENT_PREDICTIONS.json', sha256=sha(dest)), parent_gate_verification_suite=dict(path='FROZEN_PREDICTIONS.json', sha256=sha(work / 'export/FROZEN_PREDICTIONS.json'), role='Immutable X49 rule/transport verification suite; not patient27 outcome predictions'), release_status='FAIL/UNKNOWN upstream design retained; transport only'))

def run(work, tag):
    work.mkdir(parents=True)
    (work / 'raw').mkdir()
    src = D / 'results/DEMO48_PACKAGE/batch15/demos/X60/PREREG_R1.json'
    shutil.copyfile(src, work / 'PREREG_R1.json')
    records = []
    for stage in STAGES:
        budget()
        state('RUN_' + tag + '_' + stage.upper(), 'RUNNING', stage)
        cmd = command(stage, work)
        start = time.perf_counter()
        with (work / 'raw' / ('PROCESS_' + stage + '.log')).open('w') as log:
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        row = dict(stage=stage, argv=cmd, exit_code=proc.returncode, wall_seconds=time.perf_counter() - start)
        records.append(row)
        write(work / 'PROCESS_RECEIPT.json', records)
        print(tag, stage, 'exit', proc.returncode, 'seconds', round(row['wall_seconds'], 2), flush=True)
        with (R / 'COMMANDS.md').open('a') as f:
            f.write('\n`' + __import__('shlex').join(cmd) + '`\n\nExit ' + str(proc.returncode) + ', log: ' + str(work / 'raw' / ('PROCESS_' + stage + '.log')) + '\n')
        if proc.returncode:
            state('FAILED_' + tag + '_' + stage.upper(), row, 'Preserve failed attempt; inspect exact log and change construction')
            raise RuntimeError('Stage failed; preserved ' + str(work / 'raw' / ('PROCESS_' + stage + '.log')))
        if stage == 'physical':
            frozen_predictions(work)
    return records
SCIENCE = ['upper.stl', 'lower.stl', 'segmented_tooth.stl', 'antagonist.stl', 'ANTAGONIST_FACET_MAP.npz', 'ANTAGONIST_DOMAIN.json', 'measured_capped_roof.stl', 'scene.npz', 'SCAN.json', 'SDF.npz', 'SDF.json', 'prep.stl', 'crown.stl', 'DESIGN_FIELDS.npz', 'GENERATION.json', 'GENERATOR_FAULTS.json', 'CONTRACT.json', 'REGIONS.json', 'GATE_PREEXPORT.json', 'GATE.json', 'EXPORT.json', 'GATE_FAULTS.json', 'MILLING.json', 'PHYSICAL.json', 'METROLOGY.json', 'metrology/point_fields.npz', 'METROLOGY_FAULTS.json', 'LAB_ALARM.json', 'LAB_MEASUREMENTS_TEMPLATE.csv', 'export/model.3mf', 'export/model.3mf.json', 'export/DESIGN_GATE.json', 'FROZEN_PATIENT_PREDICTIONS.json', 'export/X70_FROZEN_PATIENT_PREDICTIONS.json', 'export/EXPORT_PROVENANCE.json']

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--single-run', action='store_true', help='Development attempt; not full reproducibility gate')
    a = ap.parse_args()
    check_sources()
    DATA.mkdir(parents=True, exist_ok=True)
    LOGICAL.mkdir(exist_ok=True)
    if (R / 'CODE_LOCK.json').exists():
        for q in json.loads((R / 'CODE_LOCK.json').read_text())['files']:
            if sha(R / q['path']) != q['sha256']:
                raise ValueError('Sealed X70 implementation changed: ' + q['path'])
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    attempt = DATA / 'attempts' / stamp
    attempt.mkdir(parents=True)
    start = time.perf_counter()
    runs = []
    try:
        for tag in ['A'] if a.single_run else ['A', 'B']:
            runs.append(dict(directory=str(attempt / tag), processes=run(attempt / tag, tag)))
        check_sources()
        science = {}
        for name in SCIENCE:
            p = attempt / 'A' / name
            if not p.exists():
                raise RuntimeError('Missing scientific output ' + name)
            science[name] = dict(sha256=sha(p), bytes=p.stat().st_size, second_sha256=sha(attempt / 'B' / name) if not a.single_run else None)
        diffs = [name for (name, x) in science.items() if x['second_sha256'] is not None and x['sha256'] != x['second_sha256']]
        receipt = dict(schema='X70-run-v1', attempt=str(attempt), runs=runs, scientific_artifacts=science, byte_identity_gate='NOT_RUN_SINGLE' if a.single_run else 'PASS' if not diffs else 'FAIL', different_files=diffs, wall_seconds=time.perf_counter() - start, max_child_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, data_bytes=budget(), telemetry_excluded_before_run=['PROCESS_RECEIPT.json', 'raw/PROCESS_*.log', 'GATE_TELEMETRY.json', 'raw/R1.json runtime fields', 'metrology/report.json runtime fields', 'nested gate query_s and query_wall_s'], canonical_packaging='ZIP member dates fixed to 1980, sorted member order; floating geometry unchanged. Same stable logical work path inside read-only mounts.')
        write(attempt / 'REPRODUCIBILITY.json', receipt)
        write(R / 'LATEST_RUN.json', receipt)
        state('DIGITAL_CHAIN_EXECUTED', receipt['byte_identity_gate'], 'Report matched-source support, unknown physical ports and freeze lab predictions before measurements')
        if diffs:
            raise RuntimeError('Scientific byte drift in ' + str(diffs))
        if not a.single_run:
            subprocess.run([PYTHON, '-B', str(R / 'report.py')], check=True)
            result = json.loads((R / 'results.json').read_text())
            if not result['all_faults_rejected']:
                raise RuntimeError('An injected false value survived a control')
            state('ROUND_COMPLETE_PARTIAL_CAPABILITY', 'Scientific byte identity PASS; actual design FAIL; physical ports UNKNOWN', 'Matched X70 lab measurement + separately calibrated likelihood; reviewed complete axial generator')
            print('Scientific rerun PASS; patient design FAIL; physical calibration UNKNOWN. Read README_DEMO.md.', flush=True)
    except Exception as e:
        write(attempt / 'FAILURE.json', dict(error=str(e), runs=runs, wall_seconds=time.perf_counter() - start))
        raise
if __name__ == '__main__':
    main()
