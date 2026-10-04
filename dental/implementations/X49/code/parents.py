"""Read-only, content-pinned predecessor imports. Never run their main()."""
from dental_release.paths import expand as _release_expand
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parents[1]
DENTAL = HERE.parents[1]
RESULTS = DENTAL / 'results'

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda : f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

def verify_dependencies():
    lock = json.loads((HERE / 'DEPENDENCIES.json').read_text())
    for row in lock['files']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Pinned predecessor changed: ' + row['path'])
sys.path.insert(0, str(RESULTS / 'LANE_X34_STL_DESIGN_GATE' / 'code'))
from designgate import gate as original_gate
from designgate import geometry, distance
sys.path.insert(0, str(RESULTS / 'LANE_X38_EXPORT_GATE' / 'code'))
import mesh_transport as transport
import exportgate
sys.path.insert(0, str(RESULTS / _release_expand('PROOF_LANE')))
from gencad_bench.checks import milling, inherited

def load_cement():
    path = RESULTS / _release_expand('X13') / 'cement_port.py'
    spec = importlib.util.spec_from_file_location('x49_cement', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
