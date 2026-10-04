"""One serial command, all frozen identities, original validators, retained failures."""
from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path
from prepare_replay import HERE, PACKAGE, DATA, ALIAS, sha, write, execution_manifest
from receipt_contract import validator_pass
from release_bindings import install as install_release_binding
ARTIFACT_ROOT = DATA / 'artifacts'

def size(root):
    total = 0
    seen = set()
    for (base, dirs, files) in os.walk(root, followlinks=False):
        for name in files:
            p = Path(base) / name
            if not p.is_symlink():
                s = p.stat()
                identity = (s.st_dev, s.st_ino)
                if identity not in seen:
                    total += s.st_size
                    seen.add(identity)
    return total

def state(**changes):
    p = HERE / 'CURRENT_WORK_STATE.json'
    obj = json.loads(p.read_text())
    obj.update(changes)
    obj['updated_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write(p, obj)

def original_index():
    cache = HERE / 'raw/IMMUTABLE_CONTENT_INDEX.json'
    if cache.exists():
        index = json.loads(cache.read_text())
        for manifest in list((HERE / 'raw/prior_runs').glob('*/output_manifests/*.json')) + list((HERE / 'raw/output_manifests').glob('*.json')):
            for r in json.loads(manifest.read_text()):
                if r.get('sha256') and r.get('locator') and Path(r['locator']).is_file():
                    index.setdefault(r['sha256'], {'path': r['locator'], 'bytes': r['bytes'], 'own_output': r['locator'].startswith(str(DATA))})
        return index
    index = {}
    for (base, dirs, files) in os.walk(PACKAGE, followlinks=False):
        dirs[:] = [d for d in dirs if d not in ['runs', 'work', '__pycache__', 'logs']]
        for name in files:
            p = Path(base) / name
            if p.is_symlink():
                continue
            digest = sha(p)
            index.setdefault(digest, {'path': str(p), 'bytes': p.stat().st_size})
    write(cache, index)
    return index

def archive(ident, clone, originals):
    records = []
    retained = ARTIFACT_ROOT / ident
    retained.mkdir(parents=True, exist_ok=True)
    roots = [clone / 'runs'] + list(clone.glob('batch*/runs')) + list(clone.glob('batch*/work'))
    roots += list((DATA / 'mapped_data').iterdir())
    if ident in ['X1C', 'X41', 'X44', 'X43_STATUS']:
        roots += [clone / 'batch8/demos' / ident, clone / 'batch8/logs']
    for root in roots:
        for (base, dirs, files) in os.walk(root, followlinks=False):
            for name in files:
                p = Path(base) / name
                rel = str(p.relative_to(DATA))
                if p.is_symlink():
                    records.append({'run_path': rel, 'symlink_target': os.readlink(p), 'kind': 'shared_input_or_closed_output'})
                    continue
                digest = sha(p)
                n = p.stat().st_size
                if digest in originals:
                    kind = 'reused_byte_exact_own_output' if originals[digest].get('own_output') else 'unchanged_immutable_input'
                    records.append({'run_path': rel, 'sha256': digest, 'bytes': n, 'locator': originals[digest]['path'], 'kind': kind})
                else:
                    target = retained / rel
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(p), target)
                    records.append({'run_path': rel, 'sha256': digest, 'bytes': n, 'locator': str(target), 'kind': 'new_run_output'})
                    originals[digest] = {'path': str(target), 'bytes': n, 'own_output': True}
        if root.exists():
            for q in list(root.iterdir()):
                if q.is_dir() and (not q.is_symlink()):
                    shutil.rmtree(q)
                else:
                    q.unlink()
            if ident in ['X1C', 'X41', 'X44', 'X43_STATUS'] and root == clone / 'batch8/demos' / ident:
                root.rmdir()
    write(HERE / 'raw/output_manifests' / f'{ident}.json', records)
    return {'files': len(records), 'new_bytes': sum((r.get('bytes', 0) for r in records if r['kind'] == 'new_run_output')), 'manifest': str(HERE / 'raw/output_manifests' / f'{ident}.json')}

def fresh_receipt(row, clone, start_ns):
    batch = row['batch']
    candidates = [clone / batch / 'LATEST_VERIFICATION.json'] if batch != 'base' else [clone / 'LATEST_SUBSET_VERIFICATION.json']
    if batch == 'batch8':
        candidates = [clone / 'batch8/REPLAY.json']
    if batch == 'batch9':
        candidates += [clone / 'batch9/LATEST_PATIENT360.json']
    for q in candidates:
        if q.exists() and q.stat().st_mtime_ns >= start_ns:
            return (json.loads(q.read_text()), str(q), sha(q))
    return (None, None, None)

def mismatch_values(obj):
    findings = []

    def visit(x):
        if isinstance(x, dict):
            for (k, v) in x.items():
                if k in ['mismatches', 'headline_mismatches', 'scientific_mismatches', 'differences'] and isinstance(v, list):
                    findings.extend(v)
                else:
                    visit(v)
        elif isinstance(x, list):
            for v in x:
                visit(v)
    visit(obj)
    return findings

def main():
    global ARTIFACT_ROOT
    manifest = execution_manifest()
    frozen = json.loads(manifest.read_text())
    assert sha(manifest) == Path(str(manifest) + '.sha256').read_text().split()[0]
    rows = frozen['executions']
    requested = sys.argv[1:]
    if requested:
        assert len(requested) == len(set(requested)) and set(requested) <= set((r['id'] for r in rows))
        rows = [r for r in rows if r['id'] in requested]
    clone = DATA / 'clean_package'
    (HERE / 'replay_logs').mkdir(exist_ok=True)
    ledger_file = HERE / 'raw/ALL_REPLAY.json'
    tag = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    if ledger_file.exists():
        previous = json.loads(ledger_file.read_text())
        if previous.get('status') == 'RUNNING':
            raise RuntimeError('A run is active; do not start a concurrent replay')
        prior = HERE / 'raw/prior_runs' / tag
        prior.mkdir(parents=True)
        shutil.copy2(ledger_file, prior / 'ALL_REPLAY.json')
        for name in ['receipts', 'output_manifests']:
            if (HERE / 'raw' / name).exists():
                shutil.move(HERE / 'raw' / name, prior / name)
        for name in ['ONE_COMMAND_TIME.json', 'LAST_LEDGER_PRE_TIMER_ENRICHMENT.json']:
            if (HERE / 'raw' / name).exists():
                shutil.move(HERE / 'raw' / name, prior / name)
        if (HERE / 'replay_logs').exists():
            shutil.copytree(HERE / 'replay_logs', prior / 'replay_logs')
        write(prior / 'RELOCATION.json', {'reason': 'Retain each completed attempt before a new one-command replay', 'old_raw_paths': 'raw/prior_runs/' + tag + '/original directory', 'artifact_paths': 'Original data-artifact locators remain unchanged'})
        ledger_file.unlink()
    ARTIFACT_ROOT = DATA / 'artifacts_runs' / tag
    for r in json.loads((HERE / 'CLEAN_COPY_MANIFEST.json').read_text())['copied_files']:
        src = DATA / 'metadata_snapshot' / r['relative_path']
        if src.exists():
            assert sha(src) == r['sha256'], 'Metadata snapshot drift'
            dst = clone / r['relative_path']
            assert not dst.is_symlink(), 'Copied metadata unexpectedly became a source link'
            if dst.exists():
                dst.unlink()
            shutil.copy2(src, dst)
    for q in clone.rglob('LATEST*VERIFICATION*.json'):
        if not q.is_symlink():
            q.unlink()
    originals = original_index()
    records = []
    total_start = time.monotonic()
    state(phase='ALL_DEMO_REPLAY', next_operation='Execute frozen identities serially with protected source paths')
    for row in rows:
        ident = row['id']
        print('START', ident, flush=True)
        state(current_demo=ident, next_operation='Execute ' + ident + ' with its existing validator')
        binding = install_release_binding(clone)
        if row['batch'] == 'batch8':
            parent = clone / 'batch8/demos'
            if parent.is_symlink():
                parent.unlink()
                parent.mkdir()
            assert not (parent / ident).exists(), 'Unarchived batch8 work'
            shutil.copytree(PACKAGE / 'batch8/demos' / ident, parent / ident, symlinks=True)
        immutable = Path(_release_expand('@DENTAL_WORK_ROOT@/PKG_demo48/PKG2_immutable_X1B')).resolve()
        command = ['bwrap', '--die-with-parent', '--ro-bind', '/', '/', '--proc', '/proc', '--dev', '/dev', '--ro-bind', str(PACKAGE), str(ALIAS), '--bind', str(clone), str(PACKAGE), '--bind', str(DATA / 'mapped_data/PKG_demo48'), str(Path(_release_expand('@DENTAL_WORK_ROOT@/PKG_demo48')).resolve()), '--bind', str(DATA / 'mapped_data/DEMO48_PACKAGE_BATCH4'), str(Path(_release_expand('@DENTAL_WORK_ROOT@/DEMO48_PACKAGE_BATCH4')).resolve()), '--ro-bind', str(immutable), str(immutable), '--ro-bind', str((PACKAGE / 'batch4/demos/GENCAD_V3').resolve()), str((PACKAGE / 'batch4/demos/GENCAD_V3').resolve()), '--chdir', str(PACKAGE)]
        route = list(row['route'])
        if route[0] == 'run_demos.sh':
            route[0] = './run_demos.sh'
        command += route
        ram = '4' if ident in ['PROOF_LANE_FULL_CROWN_R3', 'PATIENT360', 'PATIENT360_R2', 'GENCAD_V4', 'X44'] else '3'
        command = [str(PACKAGE.parent.parent / 'tasks/heavy_run.sh'), ram, '0'] + command
        env = dict(os.environ, OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(HERE / 'py_adapter'))
        env['DENTAL_X86_REGULAR_COPY_BYTES'] = str(json.loads((HERE / 'raw/BATCH4_TEMPLATE_COPY_COST.json').read_text())['source_regular_bytes'])
        env['DENTAL_X86_REMAINING_INTERMEDIATE_BYTES'] = str(3000000000 - size(DATA) - size(HERE))
        log = HERE / 'replay_logs' / f'{ident}.log'
        start = time.monotonic()
        start_ns = time.time_ns()
        with log.open('w') as handle:
            result = subprocess.run(command, env=env, stdout=handle, stderr=subprocess.STDOUT)
        (receipt, receipt_path, receipt_sha) = fresh_receipt(row, clone, start_ns)
        mismatches = mismatch_values(receipt) if receipt is not None else None
        receipt_ok = validator_pass(receipt, ident)
        rec = {'demo_id': ident, 'batch': row['batch'], 'command': command, 'exit_code': result.returncode, 'wall_seconds': time.monotonic() - start, 'log': str(log), 'execution_scope': row['scope'], 'validator_receipt': receipt, 'validator_path': receipt_path, 'validator_sha256': receipt_sha, 'mismatches': mismatches, 'mismatch_count': len(mismatches) if mismatches is not None else None, 'status': 'PASS' if result.returncode == 0 and receipt_ok and (not mismatches) else 'MISMATCH' if mismatches else 'EXECUTION_FAILURE', 'captured_at': datetime.datetime.now(datetime.timezone.utc).isoformat()}
        rec['operational_release_binding'] = binding
        if 'RESOURCE_LIMIT:' in log.read_text():
            rec['failure_kind'] = 'RESOURCE_LIMIT'
            rec['resource_limit'] = {'required_regular_input_copy_bytes': int(env['DENTAL_X86_REGULAR_COPY_BYTES']), 'available_retained_plus_intermediate_bytes': int(env['DENTAL_X86_REMAINING_INTERMEDIATE_BYTES']), 'budget_bytes': 3000000000, 'scientific_comparison_reached': False}
        rec['output_archive'] = archive(ident, clone, originals)
        rec['retained_bytes_after'] = size(DATA) + size(HERE)
        records.append(rec)
        write(HERE / 'raw/receipts' / f'{ident}.json', rec)
        write(ledger_file, {'requested': [r['id'] for r in rows], 'records': records, 'completed': len(records), 'wall_seconds_so_far': time.monotonic() - total_start, 'status': 'RUNNING'})
        state(latest_gate=f"{ident}: {rec['status']}", completed=len(records), next_operation='Next frozen execution; retain every failure')
        print('END', ident, rec['status'], round(rec['wall_seconds'], 3), 'seconds', flush=True)
        if rec['retained_bytes_after'] >= 3000000000:
            raise RuntimeError('Own intermediate budget reached; no further writes authorized')
    result = {'requested': [r['id'] for r in rows], 'records': records, 'completed': len(records), 'wall_seconds': time.monotonic() - total_start, 'pass_count': sum((r['status'] == 'PASS' for r in records)), 'mismatch_count': sum((r['status'] == 'MISMATCH' for r in records)), 'execution_failure_count': sum((r['status'] == 'EXECUTION_FAILURE' for r in records)), 'status': 'COMPLETE_WITH_FAILURES' if any((r['status'] != 'PASS' for r in records)) else 'PASS', 'clean_copy': 'Fresh writable package state; shared local inputs read-only; installed environments retained', 'operational_release_binding': 'R7 frozen batch25; four operational file bindings versioned, original69 code/scientific entries and all numerical comparisons retained', 'scope_cutoff': '97 identities from initial batch<=25 registry plus earlier indexed executables; live batch26 excluded', 'physical_measurements': 0, 'cold_OS_installation_tested': False}
    write(ledger_file, result)
    state(phase='REPLAY_COMPLETE', latest_gate=result['status'], next_operation='Write front and per-demo evidence table')
    print('ALL', len(records), 'PASS', result['pass_count'], 'failures', result['execution_failure_count'] + result['mismatch_count'], flush=True)
    subprocess.run([sys.executable, str(HERE / 'write_reports.py')] + (['--partial'] if requested else []), check=True)
    subprocess.run(['/usr/bin/python3', '-s', str(HERE / 'plot_package.py')], check=True)
    raise SystemExit(0 if result['status'] == 'PASS' else 1)
if __name__ == '__main__':
    main()
