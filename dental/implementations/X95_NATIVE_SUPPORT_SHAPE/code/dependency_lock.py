from common import *
import sys, importlib
X1C = D / 'results/LANE_X1C_CROWN_LITERATURE'

def capture():
    paths = [X1 / 'inputs/literature/PMC10817558.xml', X1C / 'PREREG_R2_BRACKET_CALIBRATION.json', X1C / 'raw/R2_RESULTS.json', X1 / 'RUNTIME_LOCK.json', X1 / 'SOLVER_MANIFEST.json', X1 / 'vendor/ccx_highres.py', X1 / 'vendor/crown_fit_geometry.py', D / 'cells/design/crown_case_sts3d.py', D / 'tasks/heavy_run.sh', SRC / 'code/mechanics_read_reference.txt']
    for row in read(X1 / 'SOLVER_MANIFEST.json'):
        paths.append(X1 / row['local'])
    paths = list(dict.fromkeys(paths))
    paths = [p for p in paths if p.exists()]
    obj = {'frozen_utc': now(), 'scope': 'Supplemental inherited runtime/source-object hashes; no numerical, metric or input-selection change', 'files': [{'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in paths], 'runtime': {'python': sys.version.split()[0], 'executable': sys.executable, 'packages': {m: importlib.import_module(m).__version__ for m in ['numpy', 'scipy', 'trimesh', 'matplotlib', 'skimage']}}}
    freeze(SRC / 'DEPENDENCY_LOCK.json', obj)
    print(len(paths))

def verify():
    obj = read(SRC / 'DEPENDENCY_LOCK.json')
    for r in obj['files']:
        if sha(r['path']) != r['sha256']:
            raise ValueError('Supplemental dependency drift ' + r['path'])
    if sys.version.split()[0] != obj['runtime']['python']:
        raise ValueError('Python drift')
    for (m, v) in obj['runtime']['packages'].items():
        if importlib.import_module(m).__version__ != v:
            raise ValueError('Package drift ' + m)
    print('All supplemental dependency hashes and runtime versions match')
if __name__ == '__main__':
    capture() if '--capture' in sys.argv else verify()
