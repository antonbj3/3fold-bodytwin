"""Make a standalone own-lane delivery without duplicating the large payload."""
from dental_release.paths import expand as _release_expand
import sys, shutil, subprocess, time, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
source = PAYLOAD.resolve()
base = source.parent
delivery = base / 'DELIVERY'
if delivery.exists():
    raise RuntimeError('delivery exists: inspect prior attempt; never silently replace')
original_digest = read(PAYLOAD / 'evaluation/SCORE_DIGEST.json')['sha256']
delivery.mkdir()
for p in ROOT.iterdir():
    if p.name == 'payload':
        continue
    if p.is_dir():
        shutil.copytree(p, delivery / p.name, ignore=shutil.ignore_patterns('__pycache__'))
    else:
        shutil.copy2(p, delivery / p.name)
source.rename(delivery / 'payload')
PAYLOAD.unlink()
PAYLOAD.symlink_to(delivery / 'payload', target_is_directory=True)
frozen = [str(p.relative_to(delivery)) for p in delivery.rglob('*') if p.is_symlink()]
if frozen:
    raise RuntimeError('standalone package contains symlink ' + str(frozen))
checks = 'from pathlib import Path; import json; paths=' + repr([str(ROOT), str(base), _release_expand('@DENTAL_IMPLEMENTATIONS@/GENCAD_V2'), _release_expand('@DENTAL_WORK_ROOT@/X11')]) + "; assert all(not Path(p).exists() for p in paths); assert len(Path('/proc/net/route').read_text().splitlines())==1; print(json.dumps({'unavailable_original_paths':paths,'network_route_rows':0}))"
prefix = [str(delivery / 'payload/runtime/bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--bind', str(delivery), '/demo', '--chdir', '/demo', '--setenv', 'PATH', '/usr/bin:/bin', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
probe = subprocess.run(prefix + ['/demo/payload/runtime/bin/python3', '-s', '-E', '-B', '-c', checks], capture_output=True, text=True)
if probe.returncode:
    raise RuntimeError('relocation isolation probe failed ' + probe.stderr)
start = time.perf_counter()
with (ROOT / 'raw/relocation.log').open('w') as f:
    p = subprocess.run(prefix + ['/usr/bin/bash', '/demo/run_all.sh', '--recompute'], stdout=f, stderr=subprocess.STDOUT, timeout=14400)
receipt = dict(status='PASS' if p.returncode == 0 else 'FAIL', utc=now(), exit_code=p.returncode, seconds=time.perf_counter() - start, standalone_delivery=str(delivery), creation='Own lane code/docs copied; own payload moved atomically to avoid a second multi-GB data copy. Real directories, no symlinks; unrelated /demo visible execution path.', isolation=json.loads(probe.stdout), argv=prefix + ['/usr/bin/bash', '/demo/run_all.sh', '--recompute'], original_score_digest=original_digest)
if p.returncode == 0:
    receipt['relocated_score_digest'] = read(delivery / 'payload/evaluation/SCORE_DIGEST.json')['sha256']
    receipt['same_score_digest'] = receipt['relocated_score_digest'] == original_digest
    receipt['regeneration'] = read(delivery / 'raw/RECOMPUTATION_RECEIPT.json')
    receipt['status'] = 'PASS' if receipt['same_score_digest'] and receipt['regeneration']['identical'] else 'FAIL'
    scratch = delivery / 'payload/recomputed'
    receipt['removed_recomputation_bytes'] = sum((p.stat().st_size for p in scratch.rglob('*') if p.is_file()))
    shutil.rmtree(scratch)
dump(ROOT / 'raw/RELOCATION_RECEIPT.json', receipt)
dump(delivery / 'raw/RELOCATION_RECEIPT.json', receipt)
for name in ['results.json', 'LEADERBOARD.md', 'RESULTS.md']:
    if (delivery / name).exists():
        shutil.copy2(delivery / name, ROOT / name)
for name in ['CONTROL_REPORT.json', 'SOURCE_REFERENCE_AUDIT.json', 'RECOMPUTATION_RECEIPT.json', 'RECOMPUTATION_PROCESS.json', 'LAST_RUN.json']:
    if (delivery / 'raw' / name).exists():
        shutil.copy2(delivery / 'raw' / name, ROOT / 'raw' / name)
for root in [ROOT, delivery]:
    if (root / 'results.json').exists():
        r = read(root / 'results.json')
        r['relocation'] = receipt
        dump(root / 'results.json', r)
    report = root / 'RESULTS.md'
    if report.exists():
        report.write_text(report.read_text().replace('C01 relocation status: PENDING_UNTIL_RELOCATION_RECEIPT.', 'C01 relocation status: ' + receipt['status'] + '.'))
print(json.dumps(receipt, indent=2))
if receipt['status'] != 'PASS':
    raise RuntimeError('relocation did not pass; receipt retained')
