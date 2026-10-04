"""Fresh complete R8 copied run; disjoint outputs allow parallel reference replay."""
from dental_release.paths import expand as _release_expand
import sys, subprocess, shutil, time, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from util import *
D = PAYLOAD.resolve().parent
workspace = D.parent / 'R8_relocated_outputs'
assert read(ROOT / 'raw/CONTROL_REPORT.json')['integrity']['count'] == 55
assert not workspace.exists()
workspace.mkdir()
prefix = [str(PAYLOAD / 'runtime/bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--bind', str(D), '/demo', '--chdir', '/demo', '--setenv', 'PATH', '/usr/bin:/bin', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
for name in ['evaluation', 'round3/evaluation', 'round5/evaluation', 'round6/evaluation', 'round7/evaluation']:
    p = workspace / name
    p.mkdir(parents=True)
    prefix += ['--bind', str(p), '/demo/payload/' + name]
script = 'from pathlib import Path; import json; p=' + repr([str(ROOT), str(D.parent), _release_expand('@DENTAL_IMPLEMENTATIONS@/GENCAD_V2'), _release_expand('@DENTAL_WORK_ROOT@/X11')]) + "; assert all(not Path(x).exists() for x in p); assert len(Path('/proc/net/route').read_text().splitlines())==1; print(json.dumps({'unavailable_original_paths':p,'network_route_rows':0}))"
probe = subprocess.run(prefix + ['/demo/payload/runtime/bin/python3', '-s', '-E', '-B', '-c', script], capture_output=True, text=True)
assert probe.returncode == 0, probe.stderr
start = time.perf_counter()
with (ROOT / 'raw/relocation_R8.log').open('w') as f:
    p = subprocess.run(prefix + ['/usr/bin/bash', '/demo/run_all.sh', '--recompute'], stdout=f, stderr=subprocess.STDOUT, timeout=14400)
receipt = dict(status='FAIL', revision='R8', utc=now(), exit_code=p.returncode, seconds=time.perf_counter() - start, standalone_delivery=str(D), creation='Release1 own-lane standalone tree explicitly revised under PREREG_R8. Output mounts are disjoint during concurrent reference replay; code/runtime/inputs are unchanged package bytes.', isolation=json.loads(probe.stdout), argv=prefix + ['/usr/bin/bash', '/demo/run_all.sh', '--recompute'], prior_failed_receipt='revisions/RELEASE_1/raw/RELOCATION_RECEIPT.json')
if p.returncode == 0:
    assert 'DONE all' in (ROOT / 'raw/release_replay_R8.log').read_text(), 'reference replay must finish first'
    original = read(PAYLOAD / 'evaluation/SCORES.json')
    assert original['benchmark_sha256'] == sha(ROOT / 'BENCHMARK_LOCK.json')
    receipt.update(original_score_digest=read(PAYLOAD / 'evaluation/SCORE_DIGEST.json')['sha256'], relocated_score_digest=read(workspace / 'evaluation/SCORE_DIGEST.json')['sha256'], regeneration=read(D / 'raw/RECOMPUTATION_RECEIPT.json'))
    receipt['same_score_digest'] = receipt['original_score_digest'] == receipt['relocated_score_digest']
    receipt['status'] = 'PASS' if receipt['same_score_digest'] and receipt['regeneration']['identical'] else 'FAIL'
    dump(ROOT / 'raw/R8_ORIGINAL_SCORE_DIGEST.json', read(PAYLOAD / 'evaluation/SCORE_DIGEST.json'))
    for folder in ['evaluation', 'round3/evaluation', 'round5/evaluation', 'round6/evaluation', 'round7/evaluation']:
        shutil.copytree(workspace / folder, PAYLOAD / folder, dirs_exist_ok=True)
    receipt['removed_disposable_recomputation_bytes'] = sum((q.stat().st_size for q in (PAYLOAD / 'recomputed').rglob('*') if q.is_file()))
    shutil.rmtree(PAYLOAD / 'recomputed')
    shutil.rmtree(workspace)
dump(ROOT / 'raw/RELOCATION_RECEIPT.json', receipt)
dump(D / 'raw/RELOCATION_RECEIPT.json', receipt)
for name in ['results.json', 'LEADERBOARD.md', 'RESULTS.md']:
    if (D / name).exists():
        shutil.copy2(D / name, ROOT / name)
for name in ['CONTROL_REPORT.json', 'SOURCE_REFERENCE_AUDIT.json', 'RECOMPUTATION_RECEIPT.json', 'RECOMPUTATION_PROCESS.json', 'LAST_RUN.json']:
    shutil.copy2(D / 'raw' / name, ROOT / 'raw' / name)
for root in [ROOT, D]:
    if (root / 'results.json').exists():
        out = read(root / 'results.json')
        out['relocation'] = receipt
        out['release_revision'] = read(ROOT / 'raw/R8_REBINDING.json')
        dump(root / 'results.json', out)
    report = root / 'RESULTS.md'
    if report.exists():
        report.write_text(report.read_text().replace('C01 relocation status: FAIL.', 'C01 relocation status: ' + receipt['status'] + '.'))
print(json.dumps(receipt, indent=2))
if receipt['status'] != 'PASS':
    raise RuntimeError('R8 copied replay failed; receipt retained')
