"""A real wrong-hash call through the original global code-lock loop."""
from dental_release.paths import expand as _release_expand
import json, os, subprocess, time
from pathlib import Path
from prepare_replay import HERE, DATA, PACKAGE, ALIAS, sha, write
from release_bindings import install
clone = DATA / 'clean_package'
assert not (HERE / 'raw/ALL_REPLAY.json').exists() or json.loads((HERE / 'raw/ALL_REPLAY.json').read_text())['status'] != 'RUNNING'
installed = install(clone)
original = (clone / 'PACKAGE_CODE_LOCK.json').read_bytes()
bad = json.loads(original)
for r in bad['files']:
    if r['path'] == 'demos.json':
        r['sha256'] = '0' * 64
(clone / 'PACKAGE_CODE_LOCK.json').unlink()
write(clone / 'PACKAGE_CODE_LOCK.json', bad)
immutable = Path(_release_expand('@DENTAL_WORK_ROOT@/PKG_demo48/PKG2_immutable_X1B')).resolve()
command = [str(PACKAGE.parent.parent / 'tasks/heavy_run.sh'), '3', '0', 'bwrap', '--die-with-parent', '--ro-bind', '/', '/', '--proc', '/proc', '--dev', '/dev', '--ro-bind', str(PACKAGE), str(ALIAS), '--bind', str(clone), str(PACKAGE), '--bind', str(DATA / 'mapped_data/PKG_demo48'), str(Path(_release_expand('@DENTAL_WORK_ROOT@/PKG_demo48')).resolve()), '--bind', str(DATA / 'mapped_data/DEMO48_PACKAGE_BATCH4'), str(Path(_release_expand('@DENTAL_WORK_ROOT@/DEMO48_PACKAGE_BATCH4')).resolve()), '--ro-bind', str(immutable), str(immutable), '--chdir', str(PACKAGE), '/usr/bin/python3', 'package_runner.py', 'X14']
start = time.monotonic()
env = dict(os.environ, OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMEXPR_NUM_THREADS='4', PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(HERE / 'py_adapter'))
try:
    r = subprocess.run(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (HERE / 'raw/RELEASE_LOCK_INJECTED_HASH.log').write_text(r.stdout)
    result = {'command': command, 'exit_code': r.returncode, 'wall_seconds': time.monotonic() - start, 'injected_value': '0' * 64, 'injected_field': 'demos.json expected SHA256 in own copied release-lock candidate', 'same_original_loop': 'package_runner.py main PACKAGE_CODE_LOCK check', 'rejected': r.returncode != 0 and 'Package code or lock drift: demos.json' in r.stdout, 'original_scientific_comparisons_reached': False, 'installed_contract': installed, 'scope': 'Administrative byte-identity control, own fault fixture; no physical facit'}
    write(HERE / 'raw/RELEASE_LOCK_INJECTED_HASH.json', result)
    print(result['exit_code'], result['rejected'])
    assert result['rejected']
finally:
    (clone / 'PACKAGE_CODE_LOCK.json').unlink()
    (clone / 'PACKAGE_CODE_LOCK.json').write_bytes(original)
    assert sha(clone / 'PACKAGE_CODE_LOCK.json') == installed['new_lock_sha256']
