"""Read-only adapters. No predecessor main() or output-writing code is run."""
import importlib.util
import sys
from ..io import ROOT, DENTAL, sha, result
sys.dont_write_bytecode = True

def load(name, relative):
    path = ROOT / '_vendor/cone.py' if relative == 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py' else DENTAL / 'results' / relative
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def insertion_cone(normals):
    module = load('gencad_parent_cone', 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
    return module.classify(normals, rounding_envelope=True)

def verify_cone_negative(normals, record):
    """Strengthen parent replay by checking nonempty six-chart/input binding."""
    from .exact import q
    import numpy as np
    try:
        if record.get('status') != 'NO_NONZERO_DIRECTION':
            return False
        certs = record['certificates']
        if len(certs) != 6 or {(c['axis'], c['sign']) for c in certs} != {(j, s) for j in range(3) for s in [-1, 1]}:
            return False
        m = load('gencad_bound_cone', 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
        unit = np.array(normals, float)
        unit = unit / np.linalg.norm(unit, axis=1)[:, None]
        N = np.unique(np.rint(unit * m.SCALE).astype(np.int64), axis=0)
        for c in certs:
            (A, b) = m.chart(N, c['axis'], c['sign'], q(2) / m.SCALE)
            if not c['rows'] or len(c['A']) != len(c['rows']) or len(c['b']) != len(c['rows']):
                return False
            if any((i < 0 or i >= len(A) for i in c['rows'])):
                return False
            if [[q(v) for v in a] for a in c['A']] != [A[i] for i in c['rows']]:
                return False
            if [q(v) for v in c['b']] != [b[i] for i in c['rows']]:
                return False
            if m.check_dual('inequality', c['A'], c['b'], c['y']).status != 'NEJ':
                return False
        return True
    except (KeyError, ValueError, TypeError, IndexError, ZeroDivisionError):
        return False

def insertion(planes):
    from .exact import q
    if not planes:
        return result('INVALID', 'convex cavity translational insertion')
    normals = [[float(q(x)) for x in p[:3]] for p in planes]
    parent = insertion_cone(normals)
    if all((q(p[2]) >= 0 for p in planes)):
        return result('PASS', '+z extraction / reversed -z insertion relative to contained preparation; no remote obstacles', proof='exact recession inequality a_i dot(-z) <= 0 for every cavity wall', reused_finite_cone=parent, lean_replayed_for_this_design=False)
    return result('UNKNOWN', 'convex cavity insertion', reused_finite_cone=parent, reason='Finite normal cone is not a global seating certificate')

def projected_extremum(U, L):
    import numpy as np
    module = load('gencad_parent_gap', 'LANE_NEXT_P_OCCLUSION_VALIDATION/code/continuous_gap.py')
    pairs = np.array([(i, j) for i in range(len(U)) for j in range(len(L))], dtype=int)
    a = module.extremum(np.asarray(U, float), np.asarray(L, float), pairs, True)
    b = module.extremum(np.asarray(U, float), np.asarray(L, float), pairs, False)
    return dict(candidate=a, full_control=b, parity_mm=abs(a['minimum_gap_mm'] - b['minimum_gap_mm']), scope='floating continuous projected height gap; not general 3D nonpenetration')
