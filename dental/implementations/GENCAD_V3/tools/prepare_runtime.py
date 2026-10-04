"""Copy only the installed locked interpreter, standard library and named local dependencies."""
from dental_release.paths import expand as _release_expand
import shutil, pathlib, json, hashlib, subprocess, sys
R = pathlib.Path(__file__).resolve().parents[1]
D = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V3'))
D.mkdir(exist_ok=True)
P = D / 'payload'
P.mkdir(exist_ok=True)
link = R / 'payload'
if not link.exists():
    link.symlink_to(P, target_is_directory=True)
rt = P / 'runtime'
(rt / 'bin').mkdir(parents=True, exist_ok=True)
shutil.copy2('/usr/bin/python3.10', rt / 'bin/python3')
ignore = shutil.ignore_patterns('__pycache__', 'tests', 'test', 'idlelib', 'tkinter', 'ensurepip')
if not (rt / 'lib/python3.10').exists():
    shutil.copytree('/usr/lib/python3.10', rt / 'lib/python3.10', ignore=ignore)
site = rt / 'lib/python3/dist-packages'
site.mkdir(parents=True, exist_ok=True)
old = rt / 'lib/python3.10/site-packages'
if old.exists():
    for p in old.iterdir():
        shutil.move(str(p), site / p.name)
    old.rmdir()
src = pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/.local/lib/python3.10/site-packages'))
names = ['numpy', 'scipy', 'skimage', 'trimesh', 'networkx', 'lazy_loader', 'packaging', 'imageio', 'tifffile', 'PIL']
packages = {}
for name in names:
    paths = [src / name] + list(src.glob(name + '.libs')) + list(src.glob(('scikit_image' if name == 'skimage' else 'pillow' if name == 'PIL' else name) + '-*.dist-info'))
    for p in paths:
        if p.is_dir():
            shutil.copytree(p, site / p.name, dirs_exist_ok=True, ignore=shutil.ignore_patterns('__pycache__'))
    packages[name] = str(src / name)
subprocess.run([str(rt / 'bin/python3'), '-s', '-E', '-c', 'import sys,numpy,scipy,skimage,trimesh; print(sys.version); print(numpy.__version__,scipy.__version__,skimage.__version__,trimesh.__version__);print(sys.prefix)'], check=True)
files = {str(p.relative_to(rt)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(rt.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
(R / 'RUNTIME_LOCK.json').write_text(json.dumps(dict(python_version='3.10.12', platform='Linux x86_64, glibc>=2.35; installed OS dynamic loader and glibc required', packages=packages, files=files), indent=2) + '\n')
print('runtime_bytes', sum((p.stat().st_size for p in rt.rglob('*') if p.is_file())))
