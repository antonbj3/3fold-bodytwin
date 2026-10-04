from dental_release.paths import expand as _release_expand
from fc_common import *
import subprocess
out = DATA / 'R2_predictions'
out.mkdir(exist_ok=True)
runtime = V4 / 'payload/runtime'
cmd = [str(runtime / 'bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--ro-bind', str(runtime), '/runtime', '--ro-bind', _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5/deps'), '/deps', '--ro-bind', str(DATA / 'public'), '/inputs', '--ro-bind', str(ROOT / 'code'), '/runner', '--bind', str(out), '/output', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    cmd += ['--setenv', k, '4']
cmd += ['/runtime/bin/python3', '-s', '-E', '-B', '/runner/homolog_plugin.py']
freeze(ROOT / 'FROZEN_GENERATOR_R2.json', dict(code_sha256=sha(ROOT / 'code/homolog_plugin.py'), public_sha256=sha(ROOT / 'FROZEN_PUBLIC_INPUTS.json'), prereg_sha256=sha(ROOT / 'PREREG_R2.json'), argv=cmd, withheld_reference_mounted=False))
st = time.perf_counter()
with (ROOT / 'raw/generate_R2.log').open('w') as f:
    r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
dump(ROOT / 'raw/PROCESS_R2.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - st))
if r.returncode:
    raise RuntimeError('R2 failed, preserve log')
freeze(ROOT / 'FROZEN_PREDICTIONS_R2.json', dict(files={str(p.relative_to(out)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / 'FROZEN_GENERATOR_R2.json'), reference_queries_before_freeze=0))
state('R2_PREDICTIONS_FROZEN', 'Public-only complete surfaces produced', 'Score unchanged v5 per-tooth vector; preserve all abstentions')
