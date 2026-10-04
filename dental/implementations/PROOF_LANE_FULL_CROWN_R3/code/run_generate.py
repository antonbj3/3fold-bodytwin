from dental_release.paths import expand as _release_expand
from common_r3 import *
import subprocess, time, argparse
ap = argparse.ArgumentParser()
ap.add_argument('tag')
ap.add_argument('script')
args = ap.parse_args()
out = DATA / args.tag
out.mkdir(exist_ok=True)
runtime = V4 / 'payload/runtime'
cmd = [str(runtime / 'bin/bwrap'), '--die-with-parent', '--unshare-all', '--new-session', '--clearenv', '--ro-bind', '/usr', '/usr', '--ro-bind', '/lib', '/lib', '--ro-bind', '/lib64', '/lib64', '--proc', '/proc', '--dev', '/dev', '--tmpfs', '/tmp', '--ro-bind', str(runtime), '/runtime', '--ro-bind', str(DATA.parent / 'PROOF_LANE_GENCAD_V5/deps'), '/deps', '--ro-bind', str(DATA / 'inputs'), '/inputs', '--ro-bind', str(ROOT / 'code'), '/runner', '--ro-bind', str(OLD / 'code'), '/oldcode', '--bind', str(out), '/output', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1']
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    cmd += ['--setenv', k, '4']
cmd += ['--ro-bind', _release_expand('@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation/.venv-newton/lib/python3.10/site-packages/shapely'), '/extra/shapely', '--ro-bind', _release_expand('@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation/.venv-newton/lib/python3.10/site-packages/shapely.libs'), '/extra/shapely.libs', '--setenv', 'PYTHONPATH', '/extra', '/runtime/bin/python3', '-s', '-B', '/runner/' + args.script]
freeze(ROOT / f'FROZEN_GENERATOR_{args.tag}.json', dict(argv=cmd, files={p.name: sha(p) for p in (ROOT / 'code').glob('*.py')}, prereg_sha256=sha(ROOT / f'PREREG_{args.tag}.json'), whole_target_mounted=False, added_contact_data_mounted=True))
st = time.perf_counter()
state(args.tag + '_GENERATING', 'frozen construction; target meshes not mounted', 'Hash predictions before scorer')
with (ROOT / 'raw' / f'generate_{args.tag}.log').open('w') as f:
    r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
dump(ROOT / 'raw' / f'PROCESS_{args.tag}.json', dict(argv=cmd, exit_code=r.returncode, seconds=time.perf_counter() - st))
if r.returncode:
    raise RuntimeError('generation failed; log retained')
freeze(ROOT / f'FROZEN_PREDICTIONS_{args.tag}.json', dict(files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / f'FROZEN_GENERATOR_{args.tag}.json'), target_mesh_queries_by_generator=0, contact_prescription='retrospectively observed from native scan; explicit added input', physical_measurement='NOT_RUN'))
state(args.tag + '_FROZEN', 'generation complete and predictions frozen', 'Score full anatomy and function')
