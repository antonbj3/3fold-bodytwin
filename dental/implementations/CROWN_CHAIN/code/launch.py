"""Read-only reviewed modules, own outputs only; no source-tree copies."""
from dental_release.paths import expand as _release_expand
import os, sys, json, hashlib, subprocess
from pathlib import Path
R = Path(__file__).resolve().parents[1]
B = R.parent
P = B / 'DEMO48_PACKAGE'
os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:4])
lock = R / 'DEPENDENCIES.json'
if lock.exists():
    if hashlib.sha256(lock.read_bytes()).hexdigest() != lock.with_suffix('.sha256').read_text().split()[0]:
        raise SystemExit('Dependency seal drift')
    for (path, record) in json.loads(lock.read_text())['files'].items():
        h = hashlib.sha256()
        with Path(path).open('rb') as f:
            for block in iter(lambda : f.read(1048576), b''):
                h.update(block)
        if h.hexdigest() != record['sha256']:
            raise SystemExit('Pinned dependency drift: ' + path)
mapping = {'PROOF_LANE_FULL_CROWN_R4': 'batch27/demos/PROOF_LANE_FULL_CROWN_R4', 'PROOF_LANE_FULL_CROWN_R3': 'batch27/demos/PROOF_LANE_FULL_CROWN_R3', _release_expand('GENCAD_V4'): 'batch27/demos/PROOF_LANE_GENCAD_V4', _release_expand('GENCAD_V6'): 'batch27/demos/PROOF_LANE_GENCAD_V6', _release_expand('X53'): 'batch13/demos/X53', _release_expand('X49'): 'batch12/demos/X49', _release_expand('X82'): 'batch23/demos/X82', _release_expand('X85'): 'batch25/demos/X85', _release_expand('X13'): 'demos/X13', 'LANE_X23_CROWN_WEIBULL_FLOOR': 'batch3/demos/X23', _release_expand('X1B'): 'demos/X1B', _release_expand('X61'): 'batch15/demos/X61'}
data = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_COMBINE_CROWN_CHAIN'))
data.mkdir(exist_ok=True)
if sum((p.stat().st_size for p in data.rglob('*') if p.is_file())) > 3000000000:
    raise SystemExit('Own data cap exceeded')
cmd = ['bwrap', '--die-with-parent', '--unshare-net', '--ro-bind', '/', '/', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--bind', str(R), str(R), '--bind', str(data), str(data)]
for (original, rel) in mapping.items():
    if not (P / rel).is_dir():
        raise SystemExit('Missing reviewed part ' + rel)
    cmd += ['--ro-bind', str(P / rel), str(B / original)]
cmd += ['--ro-bind', str(B / 'PROOF_LANE_GENCAD_V4/payload'), '/tmp/chain_v4_payload']
cmd += ['--chdir', str(R), _release_expand('@DENTAL_PYTHON@'), str(R / 'code/run_chain.py'), *sys.argv[1:]]
(R / 'raw/READONLY_BINDINGS.json').write_text(json.dumps(mapping, indent=2) + '\n')
os.execvp(cmd[0], cmd)
