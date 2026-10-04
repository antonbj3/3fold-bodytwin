"""Select, vendor and hash the original CPU solver and its MKL run-time files."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import shutil
from common import R, read, dump, sha
W = R.parents[1]
src = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/shared-storage/scratch/dental_F2/opt'))
manifest = []
items = [(src / 'src/CalculiX/ccx_2.23/src/ccx_2.23_pardiso', 'vendor/solver/ccx')]
for name in ['libmkl_intel_lp64.so.2', 'libmkl_gnu_thread.so.2', 'libmkl_core.so.2', 'libmkl_avx2.so.2', 'libmkl_avx512.so.2', 'libmkl_def.so.2', 'libmkl_mc3.so.2', 'libmkl_vml_avx2.so.2', 'libmkl_vml_avx512.so.2', 'libmkl_vml_def.so.2', 'libmkl_vml_mc3.so.2', 'libgomp.so.1', 'libgcc_s.so.1', 'libquadmath.so.0']:
    items.append((src / 'mklenv/lib' / name, 'vendor/solver/lib/' + name))
for name in ['crown_design_fe.py']:
    items.append((W / 'cells/design' / name, 'vendor/' + name))
items.append((W / 'cells/solvers/ccx_highres.py', 'vendor/ccx_highres.py'))
for (p, target) in items:
    q = R / target
    q.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(p, q)
    manifest.append(dict(source=str(p), local=target, sha256=sha(q), source_sha256=sha(p), bytes=q.stat().st_size))
dump('SOLVER_MANIFEST.json', manifest)
lock = read('RUNTIME_LOCK.json')
lock['packages']['skimage'] = '0.26.0'
dump('RUNTIME_LOCK.json', lock)
(R / 'requirements.lock').write_text('\n'.join((('scikit-image' if k == 'skimage' else k) + '==' + v for (k, v) in lock['packages'].items())) + '\n')
print('Solver vendor MB', sum((x['bytes'] for x in manifest)) / 1000000.0)
