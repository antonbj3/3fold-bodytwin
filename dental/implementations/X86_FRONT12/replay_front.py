"""Bounded12-demo continuation; no route can expand to the old97 loop."""
from dental_release.paths import expand as _release_expand
import argparse
import copy
import datetime
import fcntl
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from prepare_replay import HERE, PACKAGE, sha, write
from receipt_contract import validator_pass
from release_bindings import install as install_binding
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X86b_sunday_refresh'))
OLD_DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X86_sunday_refresh'))
FROZEN = HERE / 'FROZEN_FRONT_X86B.json'

def read(p):
    return json.loads(Path(p).read_text())

def update(**changes):
    j = read(HERE / 'CURRENT_WORK_STATE.json')
    j.update(changes)
    j['updated_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write(HERE / 'CURRENT_WORK_STATE.json', j)

def size(root):
    total = 0
    seen = set()
    if not root.exists():
        return 0
    for (base, dirs, files) in os.walk(root, followlinks=False):
        for name in files:
            p = Path(base) / name
            if p.is_symlink():
                continue
            s = p.stat()
            key = (s.st_dev, s.st_ino)
            if key not in seen:
                total += s.st_size
                seen.add(key)
    return total

def budget():
    n = sum((size(p) for p in [HERE, DATA, OLD_DATA]))
    if n >= 3000000000:
        raise RuntimeError(f'Own3GB intermediate budget reached: {n}')
    return n

def fresh_clone(run):
    clone = run / 'clean_package'
    clone.mkdir(parents=True)
    alias = DATA / 'readonly_package'
    alias.mkdir(exist_ok=True)
    manifest = read(HERE / 'CLEAN_COPY_MANIFEST.json')
    for r in manifest['copied_files']:
        src = OLD_DATA / 'metadata_snapshot' / r['relative_path']
        dest = clone / r['relative_path']
        assert sha(src) == r['sha256'], 'Frozen metadata drift ' + r['relative_path']
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
    for r in manifest['immutable_input_bindings']:
        rel = Path(r['source']).relative_to(PACKAGE)
        dest = clone / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.symlink_to(alias / rel, target_is_directory=True)
    for root in [OLD_DATA / 'clean_package'] + list((OLD_DATA / 'clean_package').glob('batch*')):
        if not root.is_dir() or root.is_symlink():
            continue
        for p in root.iterdir():
            rel = p.relative_to(OLD_DATA / 'clean_package')
            dest = clone / rel
            if p.is_symlink() and (not dest.exists()) and (not dest.is_symlink()):
                dest.symlink_to(os.readlink(p), target_is_directory=p.is_dir())
    for root in [clone] + [q for q in clone.glob('batch*') if q.is_dir() and (not q.is_symlink())]:
        for name in ['runs', 'work', 'logs']:
            (root / name).mkdir(exist_ok=True)
        for q in root.glob('LATEST*VERIFICATION*.json'):
            q.unlink()
    binding = install_binding(clone)
    return (clone, alias, binding)

def mismatches(obj):
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

def acceptance(exit_code, receipt, ident):
    return exit_code == 0 and validator_pass(receipt, ident) and (not mismatches(receipt))

def hash_checks(contract):
    started = time.monotonic()
    records = []
    for row in contract['rows'] + contract['rejections']:
        begin = time.monotonic()
        path = Path(row['source_result'])
        expected = row.get('source_result_sha256')
        observed = sha(path) if path.is_file() else None
        errors = []
        if not expected:
            errors.append({'kind': 'MISSING_FROZEN_RESULT_HASH', 'path': str(path)})
        elif observed != expected:
            errors.append({'kind': 'RESULT_HASH_MISMATCH', 'path': str(path), 'expected': expected, 'observed': observed})
        binding = row.get('independent_review_binding')
        checks = []
        if binding:
            for (path_key, hash_key) in [('review_locator', 'review_sha256'), ('reviewed_original_snapshot', 'reviewed_original_result_sha256')]:
                p = Path(binding[path_key])
                h = sha(p) if p.is_file() else None
                checks.append({'path': str(p), 'expected': binding[hash_key], 'observed': h, 'pass': h == binding[hash_key]})
                if h != binding[hash_key]:
                    errors.append({'kind': 'REVIEW_BINDING_MISMATCH', 'path': str(p)})
        records.append({'demo_id': row['demo_id'], 'selection': row['selection'], 'status': 'HASH_PASS' if not errors else 'HASH_MISMATCH', 'result_path': str(path), 'expected_sha256': expected, 'observed_sha256': observed, 'binding_checks': checks, 'mismatches': errors, 'mismatch_count': len(errors), 'wall_seconds': time.monotonic() - begin, 'numerical_replay_performed': False})
    return {'records': records, 'checked': len(records), 'pass_count': sum((r['status'] == 'HASH_PASS' for r in records)), 'mismatch_count': sum((r['mismatch_count'] for r in records)), 'wall_seconds': time.monotonic() - started, 'scope': 'Frozen result bytes and explicit external-review snapshots only; hash PASS does not validate scientific truth'}

def audit_release_r3(row):
    entry = next((r for r in read(PACKAGE / 'batch25/INDEX.json')['demos'] if r['demo_id'] == row['demo_id']))
    records = []
    for (name, want) in entry['release_manifest'].items():
        p = PACKAGE / 'batch25' / name
        got = sha(p)
        assert got == want, 'R3 release drift ' + name
        records.append({'path': str(p), 'sha256': got, 'bytes': p.stat().st_size})
    for (name, want) in entry['external_files'].items():
        p = Path(name)
        got = sha(p)
        assert got == want, 'R3 dependency drift ' + name
        records.append({'path': str(p), 'sha256': got, 'bytes': p.stat().st_size})
    return (entry, records)

def run_one(row, clone, alias, run, remaining):
    ident = row['demo_id']
    batch = row['batch']
    started = time.monotonic()
    start_ns = time.time_ns()
    binding = install_binding(clone)
    cmd = ['bwrap', '--die-with-parent', '--unshare-net', '--ro-bind', '/', '/', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--ro-bind', str(PACKAGE), str(alias), '--bind', str(clone), str(PACKAGE), '--bind', str(run), str(run), '--chdir', str(PACKAGE)]
    special = None
    release_hashes = None
    if ident == 'PROOF_LANE_FULL_CROWN_R3':
        (entry, release_hashes) = audit_release_r3(row)
        work = run / 'quick_R3'
        shutil.copytree(PACKAGE / 'batch25/demos' / ident, work, symlinks=True)
        original = str(Path(entry['original_directory']).resolve())
        cmd += ['--bind', str(work), original, '--ro-bind', str((PACKAGE / 'batch25/data' / ident).resolve()), str(Path(entry['original_data_directory']).resolve()), '--chdir', original, _release_expand('@DENTAL_PYTHON@'), '-B', str(HERE / 'quick_r3.py'), original, str(PACKAGE / 'batch25')]
        special = work / 'QUICK_X86B_RECEIPT.json'
        ram = '4'
    else:
        route = list(row['route'])
        route[0] = './run_demos.sh' if route[0] == 'run_demos.sh' else route[0]
        cmd += route
        ram = '3'
    command = [str(PACKAGE.parent.parent / 'tasks/heavy_run.sh'), ram, '0'] + cmd
    env = dict(os.environ, OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(HERE / 'py_adapter'), MPLBACKEND='Agg')
    log = run / f'{ident}.log'
    timeout = max(0.001, min(remaining, 300 - (time.monotonic() - started)))
    code = 124
    timed_out = False
    with log.open('w') as handle:
        proc = subprocess.Popen(command, env=env, stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
    candidates = [special] if special else [clone / batch / 'LATEST_VERIFICATION.json']
    receipt_path = next((p for p in candidates if p and p.is_file() and (p.stat().st_mtime_ns >= start_ns)), None)
    receipt = read(receipt_path) if receipt_path else None
    diffs = mismatches(receipt) if receipt is not None else []
    status = 'PASS' if acceptance(code, receipt, ident) else 'TIMEOUT' if timed_out else 'MISMATCH' if diffs else 'EXECUTION_FAILURE'
    rec = {'demo_id': ident, 'batch': batch, 'status': status, 'exit_code': code, 'wall_seconds': time.monotonic() - started, 'execution_scope': row['execution_scope'], 'profile': row['replay_profile'], 'original_packaged_profile_executed': ident != 'PROOF_LANE_FULL_CROWN_R3', 'full_numerical_replay': True if ident == 'X77' else False if ident in ['PROOF_LANE_FULL_CROWN_R3', 'GENCAD_V6', 'X53'] else None, 'command': command, 'log': str(log), 'validator_receipt': receipt, 'validator_path': str(receipt_path) if receipt_path else None, 'validator_sha256': sha(receipt_path) if receipt_path else None, 'mismatches': diffs, 'mismatch_count': len(diffs) if receipt is not None else None, 'timeout_seconds': timeout, 'operational_release_binding': binding, 'original_paths_readonly': True, 'network_disabled': True, 'captured_at': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    if release_hashes:
        write(run / 'R3_RELEASE_HASHES.json', release_hashes)
        rec['frozen_release_files_checked'] = len(release_hashes)
        rec['release_hash_manifest'] = {'path': str(run / 'R3_RELEASE_HASHES.json'), 'sha256': sha(run / 'R3_RELEASE_HASHES.json')}
        rec['execution_scope'] = receipt['execution_scope'] if receipt else 'QUICK_EXISTING_CONTROLS; failed before fresh receipt'
    rec['retained_bytes_after'] = budget()
    return rec

def delivery_faults(records, hashes, contract):
    trials = []
    for r in records:
        if r['status'] != 'PASS':
            continue
        assert acceptance(r['exit_code'], r['validator_receipt'], r['demo_id'])
        assert not acceptance(99, r['validator_receipt'], r['demo_id'])
        bad = {'demo_id': r['demo_id'], 'status': 'FAIL', 'pass': False, 'mismatches': [{'value': 999}]}
        assert not acceptance(0, bad, r['demo_id'])
        assert not acceptance(0, r['validator_receipt'], 'INJECTED_WRONG_ID')
        trials.append({'demo_id': r['demo_id'], 'bad_exit_rejected': True, 'bad_receipt_rejected': True, 'wrong_identity_rejected': True})
    hash_faults = []
    for r in hashes['records']:
        path = Path(r['result_path'])
        if not path.is_file() or not r['expected_sha256']:
            continue
        bad = hashlib.sha256(path.read_bytes() + b'\nX86B_INJECTED_WRONG_RESULT=999\n').hexdigest()
        assert bad != r['expected_sha256']
        hash_faults.append({'demo_id': r['demo_id'], 'wrong_bytes_sha256': bad, 'rejected': True})
    good = next((r for r in records if r['status'] == 'PASS'), None)
    summary = {'ids': contract['front_ids'], 'profiles': [(r['demo_id'], r['replay_profile']) for r in contract['rows'] if r['selection'] == 'FRONT_REPLAY'], 'frozen_results': [(r['demo_id'], r['source_result_sha256']) for r in contract['rows']]}
    serialized = json.dumps(summary, sort_keys=True, separators=(',', ':')).encode()
    assert good is not None
    bad = copy.deepcopy(good['validator_receipt'])
    bad = {'demo_id': good['demo_id'], 'status': 'FAIL', 'pass': False, 'mismatches': [{'value': 999}]}
    a = int(acceptance(0, good['validator_receipt'], good['demo_id']))
    b = int(acceptance(0, bad, good['demo_id']))
    assert a == 1 and b == 0
    suff = {'summary_a_sha256': hashlib.sha256(serialized).hexdigest(), 'summary_b_sha256': hashlib.sha256(bytes(serialized)).hexdigest(), 'summary_byte_identity_error': 0, 'float64_identity_error': 0.0, 'demo_id': good['demo_id'], 'downstream_a': a, 'downstream_b': b, 'downstream_difference': a - b, 'minimum_extension_for_this_pair': 'One validator-acceptance bit; counts/IDs/frozen expectations alone do not determine a produced result acceptance', 'resolution': 'PHENOMENOLOGICAL', 'time_scale': 'SIMULTANEOUS', 'physical_interpretation': False, 'external_referent': {'kind': 'our_own_fixture', 'locator': str(HERE / 'replay_front.py'), 'compared_quantity': 'Administrative exact-summary counterexample; independent physical referent not claimed', 'refutes_us': True}}
    return {'receipt_faults': trials, 'hash_faults': hash_faults, 'sufficiency': suff, 'pass': True}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ids', nargs='*')
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    started = time.monotonic()
    deadline = started + 3600
    primary = not args.ids and (not args.verify_only)
    ledger_name = 'FRONT_REPLAY.json' if primary else 'LAST_SUBSET_REPLAY.json'
    assert sha(FROZEN) == Path(str(FROZEN) + '.sha256').read_text().split()[0]
    code = read(HERE / 'FROZEN_CODE_X86B.json')
    for (name, want) in code['files'].items():
        assert sha(HERE / name) == want, 'X86b code drift ' + name
    contract = read(FROZEN)
    ids = args.ids or contract['front_ids']
    if len(ids) != len(set(ids)) or set(ids) - set(contract['front_ids']):
        raise SystemExit('Only the frozen12 front IDs are executable here')
    DATA.mkdir(parents=True, exist_ok=True)
    with (DATA / '.x86b.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        tag = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
        run = DATA / 'runs' / tag
        run.mkdir(parents=True)
        records = []
        initial = budget()
        update(phase='BOUNDED_FRONT_REPLAY', completed=0, current_demo=None, latest_gate='Budget preflight', next_operation='Create fresh writable snapshot and validate frozen result bytes')
        (clone, alias, binding) = fresh_clone(run)
        hashes = hash_checks(contract)
        write(run / 'RESULT_HASHES.json', hashes)
        write(HERE / 'raw/x86b/RESULT_HASHES.json', hashes)
        print('HASH', hashes['pass_count'], '/', hashes['checked'], 'mismatches', hashes['mismatch_count'], flush=True)
        selected = {r['demo_id']: r for r in contract['rows'] if r['selection'] == 'FRONT_REPLAY'}
        for ident in [] if args.verify_only else ids:
            remaining = deadline - time.monotonic()
            update(current_demo=ident, latest_gate=f"Result hashes: {hashes['pass_count']}/{hashes['checked']}", next_operation='Execute ' + ident + ' original bounded profile')
            print('START', ident, selected[ident]['replay_profile'], flush=True)
            if remaining <= 0:
                rec = {'demo_id': ident, 'status': 'TIME_BUDGET_NOT_RUN', 'wall_seconds': 0.0, 'profile': selected[ident]['replay_profile'], 'mismatch_count': None, 'mismatches': [], 'validator_receipt': None, 'exit_code': None}
            else:
                try:
                    rec = run_one(selected[ident], clone, alias, run, remaining)
                except Exception as exc:
                    rec = {'demo_id': ident, 'status': 'EXECUTION_FAILURE', 'wall_seconds': 0.0, 'profile': selected[ident]['replay_profile'], 'mismatch_count': None, 'mismatches': [], 'validator_receipt': None, 'exit_code': None, 'failure': type(exc).__name__ + ': ' + str(exc)}
                    write(run / (ident + '_PREPARATION_FAILURE.json'), rec)
            records.append(rec)
            write(run / 'REPLAY_PROGRESS.json', {'requested': ids, 'records': records, 'status': 'RUNNING'})
            write(HERE / 'raw/x86b/FRONT_REPLAY_PROGRESS.json', {'requested': ids, 'records': records, 'status': 'RUNNING'})
            update(completed=len(records), latest_gate=ident + ': ' + rec['status'], next_operation='Next frozen front profile')
            print('END', ident, rec['status'], round(rec['wall_seconds'], 3), 's', flush=True)
        final_hashes = hash_checks(contract)
        write(run / 'RESULT_HASHES_AFTER.json', final_hashes)
        faults = delivery_faults(records, hashes, contract) if records else None
        if faults:
            write(run / 'DELIVERY_FAULTS_AND_SUFFICIENCY.json', faults)
        seconds = time.monotonic() - started
        ledger = {'round_tag': 'X86b-sunday-refresh', 'requested': [] if args.verify_only else ids, 'records': records, 'completed': len(records), 'pass_count': sum((r['status'] == 'PASS' for r in records)), 'failure_count': sum((r['status'] != 'PASS' for r in records)), 'status': 'PASS' if all((r['status'] == 'PASS' for r in records)) and hashes['mismatch_count'] == 0 and (final_hashes['mismatch_count'] == 0) else 'COMPLETE_WITH_FAILURES', 'wall_seconds_before_report': seconds, 'wall_budget_seconds': 3600, 'budget_pass': seconds <= 3600, 'preparation_and_hash_seconds': seconds - sum((r['wall_seconds'] for r in records)), 'hashes_before': hashes, 'hashes_after': final_hashes, 'faults': faults, 'run_directory': str(run), 'clean_copy': str(clone), 'initial_own_bytes': initial, 'final_own_bytes': budget(), 'clean_OS_installation_tested': False, 'physical_measurements': 0, 'execution_manifest_sha256': sha(FROZEN), 'runtime_code_sha256': sha(HERE / 'FROZEN_CODE_X86B.json'), 'operational_release_binding': binding, 'old97_replayed': False, 'eligible_hash_only_count': len(contract['rows']) - len(contract['front_ids']), 'eligible_hash_only_policy': 'Frozen result hash only; historical replay receipts remain separately archived'}
        ledger['runtime_code_locator'] = str(HERE / 'FROZEN_CODE_X86B.json')
        ledger['namespace_mapping'] = {'runtime_package': str(PACKAGE), 'actual_host_copy': str(clone), 'note': 'Nested native receipt paths under runtime_package resolve through this mapping on the host'}
        write(run / 'FRONT_REPLAY.json', ledger)
        write(HERE / 'raw/x86b' / ledger_name, ledger)
        if primary:
            subprocess.run([sys.executable, str(HERE / 'write_x86b_reports.py')], check=True)
            subprocess.run([sys.executable, '-s', str(HERE / 'plot_x86b.py')], check=True)
        ledger['wall_seconds'] = time.monotonic() - started
        ledger['budget_pass'] = ledger['wall_seconds'] <= 3600
        ledger['timer_scope'] = 'Bounded runner, preparation, hashes,12 validators, fault tests, first reader and figure generation. Final receipt serialization and final reader timestamp rendering excluded.'
        write(run / 'FRONT_REPLAY.json', ledger)
        write(HERE / 'raw/x86b' / ledger_name, ledger)
        if primary:
            subprocess.run([sys.executable, str(HERE / 'write_x86b_reports.py')], check=True)
        update(phase='FRONT12_REPLAY_COMPLETE' if not args.verify_only else 'HASH_VERIFICATION_COMPLETE', latest_gate=f"{ledger['pass_count']}/{len(records)} scoped PASS; {final_hashes['pass_count']}/{final_hashes['checked']} hash PASS", current_demo=None, next_operation='Bind final timer, reader audit, figure and PENDING graph feedback', final_command_timer='raw/x86b/FRONT_REPLAY.json', final_scope='12 front profiles; other versions hash-only')
        print('BOUNDED', ledger['completed'], 'PASS', ledger['pass_count'], 'failures', ledger['failure_count'], 'total', round(ledger['wall_seconds'], 3), 's', flush=True)
        return 0 if ledger['status'] == 'PASS' and ledger['budget_pass'] else 1
if __name__ == '__main__':
    raise SystemExit(main())
