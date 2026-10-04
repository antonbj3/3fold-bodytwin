"""Scoped imports retain predecessor files read-only and avoid common.py collisions."""
from dental_release.paths import expand as _release_expand
from cadlib import *
from contextlib import contextmanager

@contextmanager
def scope(path, names):
    oldpath = list(sys.path)
    saved = {n: sys.modules.get(n) for n in names}
    try:
        for n in names:
            sys.modules.pop(n, None)
        sys.path.insert(0, str(path))
        yield
    finally:
        sys.path[:] = oldpath
        for (n, x) in saved.items():
            if x is None:
                sys.modules.pop(n, None)
            else:
                sys.modules[n] = x

def parents():
    sys.path.insert(0, str(BASE / _release_expand('GENCAD_V2')))
    from gencad_bench_v2 import geometry, axial_round3, tasks, generators, build
    with scope(BASE / 'PROOF_LANE_FULL_CROWN_R4/code', ['common', 'geometry', 'generate', 'generate_c', 'generate_d']):
        r4g = module('geometry', BASE / 'PROOF_LANE_FULL_CROWN_R4/code/geometry.py')
    with scope(BASE / 'PROOF_LANE_GENCAD_V6/code', ['common', 'contact']):
        v6 = module('contact', BASE / 'PROOF_LANE_GENCAD_V6/code/contact.py')
    return dict(v2g=geometry, cap=axial_round3, tasks=tasks, gen=generators, build=build, r4=r4g, v6=v6)

def cone_module():
    return module('insertion_parent', BASE / 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')

def force_refusal(task_id):
    with scope(BASE / 'LANE_X85_OCCLUSAL_ADJUST/code', ['common', 'regional', 'geometry', 'plugin']):
        p = module('plugin', BASE / 'LANE_X85_OCCLUSAL_ADJUST/code/plugin.py')
        return p.generate(dict(task_id=task_id))

def milling(crown, contract):
    with scope(BASE / 'LANE_X89_MILLING_CODESIGN/code', ['codesign', 'milling_gate', 'calibration', 'sliver_scene', 'collision']):
        m = module('milling_gate', BASE / 'LANE_X89_MILLING_CODESIGN/code/milling_gate.py')
        return m.evaluate(str(crown), str(contract))
