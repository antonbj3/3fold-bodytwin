from common import *
import subprocess
base = DATA / 'PREDICTIONS/control_full_kernel'
base.mkdir(exist_ok=True)
b = BENCH / 'payload'
outcost = ROOT / 'raw/full_kernel_cost'
outcost.mkdir(exist_ok=True)
cmd = [str(b / 'runtime/bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--dir', '/work', '--chdir', '/work', '--ro-bind', str(b / 'runtime'), '/runtime', '--ro-bind', str(b / 'participant_code'), '/baseline', '--dir', '/candidate', '--ro-bind', str(ROOT / 'code/full_kernel_participant.py'), '/candidate/participant.py', '--ro-bind', str(ROOT / 'code/public_byte_guard.py'), '/candidate/guard.py', '--ro-bind', str(ROOT / 'PUBLIC_INPUT_LOCK.json'), '/public_lock.json', '--ro-bind', str(ROOT / 'code/generator.py'), '/candidate/generator.py', '--ro-bind', str(b / 'public'), '/inputs', '--ro-bind', str(DATA / 'FULL_KERNEL_CONTROL'), '/models', '--ro-bind', str(DATA / 'FINAL_MODELS/kernel_contact'), '/contacts', '--bind', str(base), '/output', '--bind', str(outcost), '/cost', '--setenv', 'HOME', '/tmp', '--setenv', 'PYTHONWARNINGS', 'ignore']
for var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    cmd += ['--setenv', var, '1']
cmd += ['/runtime/bin/python3', '-s', '-E', '-B', '/candidate/guard.py', '--inputs', '/inputs', '--models', '/models', '--contacts', '/contacts', '--output', '/output', '--cost-output', '/cost/COST.json']
tic = time.perf_counter()
with (ROOT / 'raw/full_kernel_generation.log').open('w') as f:
    p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
dump(ROOT / 'raw/FULL_KERNEL_ISOLATION_RECEIPT.json', dict(argv=cmd, exit_code=p.returncode, seconds=time.perf_counter() - tic, reference_mount=False, scorer_mount=False, network='unshare-all', candidate_modification=False))
if p.returncode:
    raise RuntimeError('additional standard control generation failed')
